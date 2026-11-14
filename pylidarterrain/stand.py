"""Stand state and dieback from the vertical distribution of LiDAR returns.

Two products, at two different units of analysis, both built from the same
weighted height distribution:

  `dieback_index`  a continuous index on ANY unit (a grid cell or a crown
                   object).  High only where a tall mode and a dense low layer
                   coexist - the died-back-and-reshot form.
  `stand_state`    a four-state classification of a grid cell, the landscape
                   product: open, low shrubland, intact tall forest, or tall
                   with a dense low layer.

Everything takes plain arrays of heights above ground `r` and pulse weights
`w`, and returns plain arrays.  No shared state, no file IO, no projection
handling, and every threshold is a keyword argument - none of the Avonside
values are baked in.

WEIGHTS.  `w` is the per-return pulse weight, conventionally 1/number_of_returns
so that one emitted pulse contributes 1.0 however many times it is recorded.
Pass `w=None` for unweighted counts.

NORMALISATION, which is the one thing to get right.  `dieback_index` divides
the low-layer weight by TOTAL pulse weight, including returns from bare
ground - not by canopy weight.  A canopy-weight share is scale-free and so
cannot tell a sparse tall stand from a dense one; at Avonside that is exactly
why the earlier share-based score failed at the field-confirmed dieback site.

VALIDATION STATUS.  See `docs/methods.md`.  In short: the index separates
died-back from shrubland strongly, and from intact tall forest only weakly and
only once stand height is controlled for.  It is not a validated classifier.
"""
import numpy as np

__all__ = ["top_height", "band_weight", "dieback_index", "stand_state",
           "profile_null", "MIN_RETURNS"]

#: Default floor on vegetation returns per unit for a profile to be reported.
#: Below this the vertical distribution is not estimable and the unit must be
#: reported as NOT ASSESSED, never as a zero score - the units that fall below
#: it are systematically the low, dense ones, i.e. the died-back form itself.
MIN_RETURNS = 30


def _as_arrays(r, w):
    r = np.asarray(r, dtype=float)
    w = np.ones_like(r) if w is None else np.asarray(w, dtype=float)
    if r.shape != w.shape:
        raise ValueError("r and w must have the same shape")
    return r, w


def top_height(r, w=None, n_top=3, min_height=0.5):
    """Top height as the mean of the `n_top` highest vegetation returns.

    Averaging the highest few rather than taking the single maximum guards
    against one spurious high return.  At Avonside the spread between the
    highest and third-highest return was under 1 m in 44 of the 50 tallest
    cells, so the mean of three is close to the true crown top there.

    Returns nan when fewer than `n_top` returns clear `min_height`.
    """
    r, w = _as_arrays(r, w)
    veg = r[r > min_height]
    if veg.size < n_top:
        return np.nan
    return float(np.sort(veg)[-n_top:].mean())


def band_weight(r, w=None, band=(0.5, 4.0), normalise="total"):
    """Fraction of pulse weight falling in a height band.

    normalise="total"   divide by the weight of ALL returns, ground included.
                        A cover/density: "how much of this patch has material
                        in the band".  Use this for `dieback_index`.
    normalise="canopy"  divide by the weight of returns above the band floor.
                        A share: "where in the canopy does the weight sit".
                        Scale-free, and therefore blind to a thinning canopy.
    """
    r, w = _as_arrays(r, w)
    lo, hi = band
    num = float(w[(r >= lo) & (r < hi)].sum())
    if normalise == "total":
        den = float(w.sum())
    elif normalise == "canopy":
        den = float(w[r >= lo].sum())
    else:
        raise ValueError("normalise must be 'total' or 'canopy'")
    return num / den if den > 0 else np.nan


def dieback_index(r, w=None, low_band=(0.5, 4.0), height_ref=15.0, power=2.0,
                  n_top=3, min_returns=MIN_RETURNS):
    """Continuous dieback index for one unit.  Range 0 to 1, higher = more
    died-back in form.

        index = min(top_height / height_ref, 1) ** power  *  low-layer density

    The first factor gates on the ABSOLUTE presence of a tall mode, so plain
    shrubland cannot score highly however dense its low layer.  The second is
    the low-layer density, normalised by total pulse weight (see module
    docstring).  High only where both hold: relict tall stems with a thicket
    beneath them.

    Returns nan when fewer than `min_returns` vegetation returns are present.
    That is a REFUSAL TO MEASURE and must be carried through to any summary as
    an unassessed fraction, not replaced by zero.

    PARAMETERS THAT WILL NOT TRANSFER UNCHANGED.  `height_ref` is the height at
    which a stand counts as fully "tall"; 15 m suits the montane eucalypt forest
    it was built on and should be set from the mature stand height of the target
    forest.  `low_band` is the regrowth layer; its floor must clear the local
    grass and ground-noise height.  Both were chosen on a 0.8 pulses/m2 capture
    against a single reference site - treat them as starting values to be
    re-swept, not as constants.
    """
    r, w = _as_arrays(r, w)
    if int((r > low_band[0]).sum()) < min_returns:
        return np.nan
    h = top_height(r, w, n_top=n_top, min_height=low_band[0])
    if not np.isfinite(h):
        return np.nan
    gate = min(h / float(height_ref), 1.0) ** float(power)
    low = band_weight(r, w, band=low_band, normalise="total")
    if not np.isfinite(low):
        return np.nan
    return float(gate * low)


def stand_state(top_h, low_cover, is_tree=True, tall_cut=10.0, low_cut=0.06):
    """Four-state classification of a grid cell from two scalars.

    Two stages, in this order - the order matters, because the low-layer term
    is only meaningful once height has separated shrubland out:

      not a tree cell              -> "open"
      top height  < tall_cut       -> "low shrubland"
      tall, low_cover <  low_cut   -> "intact tall forest"
      tall, low_cover >= low_cut   -> "tall with dense low layer"

    `top_h` may be an empirical top height or a model-fitted canopy top (see
    `canopy.canopy_top`); `low_cover` is `band_weight(..., normalise="total")`.
    `is_tree` is a per-cell flag that a canopy exists at all - e.g. the dBIC
    test in `canopy.fit_cell`, or simply a minimum canopy cover.

    `tall_cut` and `low_cut` are JUDGEMENTS, not natural breaks: the underlying
    distributions are unimodal.  Report their sensitivity on any new dataset
    rather than quoting a single split.  Accepts scalars or arrays.
    """
    top_h = np.asarray(top_h, dtype=float)
    low_cover = np.asarray(low_cover, dtype=float)
    is_tree = np.broadcast_to(np.asarray(is_tree, dtype=bool), top_h.shape)
    out = np.full(top_h.shape, "open", dtype=object)
    scored = is_tree & np.isfinite(top_h)
    tall = scored & (top_h >= tall_cut)
    out[scored & (top_h < tall_cut)] = "low shrubland"
    out[tall & (low_cover < low_cut)] = "intact tall forest"
    out[tall & (low_cover >= low_cut)] = "tall with dense low layer"
    return out if out.ndim else out.item()


def profile_null(r_by_unit, w_by_unit, n_iter=200, rng=None, **index_kw):
    """Sampling-noise null for `dieback_index`.  NOT OPTIONAL - run it before
    reporting any index distribution or any clustering of one.

    For each unit, its vegetation returns are replaced by draws from the POOLED
    relative-height distribution of all units, rescaled to that unit's own
    observed top height.  Unit size, return count, pulse weights and top height
    are all preserved, so the null differs from the data only in having no real
    structural difference between units.

    At Avonside this control changed two conclusions in opposite directions: it
    confirmed 85% of the index's variance as real structure, and it showed that
    the apparent clusters in the index were reproduced by pure noise.  It also
    showed that a difference between two SITES can be produced by their height
    composition alone, so compare observed contrasts against this null rather
    than against 0.5.

    WHAT THE NULL PRESERVES, and therefore cannot test: unit area, return
    count, pulse weights and top height.  Any contrast arising from those is
    reproduced by the null - a height-only difference between two sites comes
    back at roughly half its observed size, and a difference in return density
    at about a third.  Only the excess of the observed contrast OVER the null
    is attributable to vertical arrangement.  This is not a pure height
    control, and it is not a significance test on its own.

    Returns an array of shape (n_iter, n_units).
    """
    rng = np.random.default_rng() if rng is None else rng
    n = len(r_by_unit)
    low0 = index_kw.get("low_band", (0.5, 4.0))[0]
    n_top = index_kw.get("n_top", 3)
    tops, pool = [], []
    for r in r_by_unit:
        r = np.asarray(r, dtype=float)
        h = top_height(r, None, n_top=n_top, min_height=low0)
        tops.append(h)
        if np.isfinite(h) and h > 0:
            veg = r[r > low0]
            if veg.size:
                pool.append(np.clip(veg / h, 0.0, 1.2))
    if not pool:
        raise ValueError("no unit has enough vegetation returns to build a null")
    pool = np.concatenate(pool)
    out = np.full((n_iter, n), np.nan)
    for it in range(n_iter):
        for i, (r, w) in enumerate(zip(r_by_unit, w_by_unit)):
            r = np.asarray(r, dtype=float)
            h = tops[i]
            if not np.isfinite(h):
                continue
            veg = r > low0
            rn = np.zeros_like(r)
            rn[veg] = rng.choice(pool, int(veg.sum())) * h
            out[it, i] = dieback_index(rn, w, **index_kw)
    return out
