"""Canopy structure from classified airborne LiDAR, one vertical distribution per cell.

Instead of segmenting individual crowns, every grid cell gets its own ground
PLANE and a pulse-weighted mixture of a Gaussian ground component and a
skew-normal canopy component fitted to heights above that plane.  The plane
matters on steep ground: at 845 Avonside the tall-tree cells sit at a median
26.6 deg, so the bare earth falls a median 2.8 m (worst case 7.5 m) across a
10 m cell.  A plane absorbs that fall, leaving a fitted ground roughness of
~0.12-0.19 m instead of metres, so canopy heights are measured from the local
ground surface rather than from a cell-mean elevation.

Returns are weighted by 1 / NumberOfReturns, so each emitted PULSE contributes
total weight 1.0 and is split between the strata it hit.  Treating returns as
independent points would count a multi-return pulse several times and inflate
the canopy against the ground.  Drop overlap (class 12) and noise (class 7)
before calling anything here: on an ELVIS tile class 12 can be 22% of returns,
and a pulse-normalised metric computed against that denominator is wrong.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, skewnorm

__all__ = ["ground_plane", "heights_above_plane", "CellFit", "fit_cell",
           "canopy_top", "stratum_cover", "canopy_skewness", "dieback_score",
           "fit_cloud"]


def ground_plane(x, y, z, cx=None, cy=None, iters=4, tukey=4.685):
    """Robust ground plane through classified ground returns, by IRLS with a
    Tukey bisquare loss.  Returns (b0, bx, by) where the surface elevation is
    ``b0 + bx*(x-cx) + by*(y-cy)`` and b0 is therefore the ground elevation at
    the cell centre.  Needs >= 3 points; >= 8 is advisable."""
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    cx = float(np.mean(x)) if cx is None else float(cx)
    cy = float(np.mean(y)) if cy is None else float(cy)
    A = np.c_[np.ones(len(x)), x - cx, y - cy]
    beta = np.linalg.lstsq(A, z, rcond=None)[0]
    for _ in range(iters):
        res = z - A @ beta
        s = 1.4826 * np.median(np.abs(res - np.median(res))) + 1e-3
        u = res / (tukey * s)
        wt = np.where(np.abs(u) < 1, (1 - u ** 2) ** 2, 0.0)
        if wt.sum() < 3:
            break
        beta = np.linalg.lstsq(A * wt[:, None], z * wt, rcond=None)[0]
    return tuple(float(v) for v in beta)


def heights_above_plane(x, y, z, plane, cx, cy):
    """Height of every return above the fitted ground plane."""
    b0, bx, by = plane
    return np.asarray(z, float) - (b0 + bx * (np.asarray(x, float) - cx)
                                   + by * (np.asarray(y, float) - cy))


class CellFit:
    """Fitted two-component mixture for one cell.

    pi_c      canopy weight fraction (an interception/gap-fraction cover)
    mu_g, sigma_g   ground component; sigma_g is roughness ABOUT the plane
    xi, omega, alpha  skew-normal canopy location, scale and shape
    dbic      BIC improvement over a ground-only fit; > 0 means trees present
    """
    __slots__ = ("pi_c", "mu_g", "sigma_g", "xi", "omega", "alpha", "dbic", "n_eff")

    def __init__(self, pi_c, mu_g, sigma_g, xi, omega, alpha, dbic, n_eff):
        self.pi_c, self.mu_g, self.sigma_g = pi_c, mu_g, sigma_g
        self.xi, self.omega, self.alpha = xi, omega, alpha
        self.dbic, self.n_eff = dbic, n_eff

    @property
    def tree(self):
        return self.dbic > 0

    def __repr__(self):
        return (f"CellFit(pi_c={self.pi_c:.3f}, sigma_g={self.sigma_g:.3f}, "
                f"xi={self.xi:.2f}, omega={self.omega:.2f}, alpha={self.alpha:.2f}, "
                f"dbic={self.dbic:.1f})")


def _nll_mix(th, r, w):
    lp, mg, lsg, xi, lom, al = th
    pc = 1.0 / (1.0 + np.exp(-lp))
    f = (1 - pc) * norm.pdf(r, mg, np.exp(lsg)) + pc * skewnorm.pdf(r, al, xi, np.exp(lom))
    return -np.sum(w * np.log(f + 1e-300))


def _nll_ground(th, r, w):
    return -np.sum(w * np.log(norm.pdf(r, th[0], np.exp(th[1])) + 1e-300))


def fit_cell(r, w, canopy_cut=1.5, xi_min=0.5, sigma_g_max=0.60, alpha_max=20.0):
    """Fit the ground + skew-normal canopy mixture to one cell's heights.

    `xi_min` floors the canopy location.  Keep it well BELOW the height of the
    lowest stratum you care about: a 2 m regrowth mode fitted against a 1.5 m
    floor sits on its constraint and is biased upward.  `sigma_g_max` caps
    ground roughness; raising it from 0.40 to 0.60 m changes little (1% of
    tree cells press the bound at Avonside) but costs nothing.

    Returns a CellFit, or None if the cell holds too little to fit.
    """
    r, w = np.asarray(r, float), np.asarray(w, float)
    keep = np.isfinite(r) & np.isfinite(w)
    r, w = r[keep], w[keep]
    if len(r) < 10:
        return None
    n_eff = float(w.sum())
    bounds = [(-8, 8), (-1, 1), (np.log(0.03), np.log(sigma_g_max)),
              (xi_min, 40), (np.log(0.3), np.log(15)), (-alpha_max, alpha_max)]
    g = minimize(_nll_ground, [0.0, np.log(0.15)], args=(r, w), method="L-BFGS-B",
                 bounds=[(-1, 1), (np.log(0.03), np.log(sigma_g_max))])
    hi = r[r > canopy_cut]
    w_can = float(w[r > canopy_cut].sum())
    if w_can <= 0:
        return CellFit(0.0, float(g.x[0]), float(np.exp(g.x[1])),
                       np.nan, np.nan, np.nan, -np.inf, n_eff)
    xi0 = max(xi_min, float(np.percentile(hi, 25))) if len(hi) > 3 else 2.0
    best = None
    for a0 in (-5.0, 0.0, 5.0):
        th0 = [np.log(max(w_can, 1e-3) / max(n_eff - w_can, 1e-3)), 0.0,
               np.log(0.15), xi0, np.log(2.0), a0]
        th0 = [min(max(v, b[0]), b[1]) for v, b in zip(th0, bounds)]
        o = minimize(_nll_mix, th0, args=(r, w), method="L-BFGS-B", bounds=bounds)
        if best is None or o.fun < best.fun:
            best = o
    lp, mg, lsg, xi, lom, al = best.x
    dbic = (2 * g.fun + 2 * np.log(n_eff)) - (2 * best.fun + 6 * np.log(n_eff))
    return CellFit(1.0 / (1.0 + np.exp(-lp)), float(mg), float(np.exp(lsg)),
                   float(xi), float(np.exp(lom)), float(al), float(dbic), n_eff)


def canopy_top(fit, q=0.95):
    """Canopy top as a high quantile of the fitted canopy component.

    Do NOT report the canopy component's MEAN as a tree height: it is a crown
    CENTROID and sits a median 4.8 m below the highest return at Avonside, a
    gap that widens with stand height.  The 0.95 quantile tracks the empirical
    95th percentile closely (farm medians 12.48 m vs 12.69 m).
    """
    if fit is None or not np.isfinite(fit.alpha):
        return np.nan
    return float(skewnorm.ppf(q, fit.alpha, fit.xi, fit.omega))


def stratum_cover(r, w, bands=((1.5, 4.0), (4.0, 8.0), (8.0, np.inf))):
    """Pulse-weight cover in each height band, as a fraction of TOTAL cell
    weight.  The default bands are the low regrowth layer (1.5-4 m), the
    midstorey (4-8 m) and the RELICT COVER above 8 m - the surviving tall
    stems and dead branches of a stand that has died back.  Band covers are interception (gap-fraction) measures: within-crown
    gaps are excluded, so they run lower than a crown-projection cover."""
    r, w = np.asarray(r, float), np.asarray(w, float)
    tot = max(float(w.sum()), 1e-9)
    return [float(w[(r >= lo) & (r < hi)].sum()) / tot for lo, hi in bands]


def canopy_skewness(r, w, canopy_cut=1.5, min_sd=0.2):
    """Pulse-weighted skewness of canopy returns.

    An intact closed crown surface is strongly NEGATIVELY skewed - weight packed
    just under the top with a thin tail downward.  A died-back stand is closer
    to symmetric, with weight spread down the profile (Avonside medians: -0.59
    intact against -0.26 at the dieback reference site).
    """
    r, w = np.asarray(r, float), np.asarray(w, float)
    m = r > canopy_cut
    if w[m].sum() < 3:
        return np.nan
    r2, w2 = r[m], w[m]
    mu = np.average(r2, weights=w2)
    sd = np.sqrt(np.average((r2 - mu) ** 2, weights=w2))
    if sd <= min_sd:
        return np.nan
    return float(np.average(((r2 - mu) / sd) ** 3, weights=w2))


def dieback_score(r, w, method="lowlayer", shrub_band=(1.5, 4.0), relict_cut=8.0,
                  canopy_cut=1.5):
    """EXPERIMENTAL. A per-cell score for forest that has died back and reshot.

    Neither method is a validated classifier; both are published so they can be
    tested on captures dense enough to resolve relict stems.  Measured at
    Avonside against a field-identified tall-dieback site, on TALL cells only
    (canopy top >= 14 m, so height is not doing the work):

      method="lowlayer"     AUC 0.692 vs intact forest - low-stratum cover,
                            the strongest single indicator found.
      method="conjunction"  AUC 0.582.  2*sqrt(shrub_dominance * relict_presence)
                            as originally specified: high only where a dense low
                            layer and tall elements coexist.  It FAILS because
                            relict presence is essentially a height measure, and
                            within tall stands it is slightly LOWER at the
                            dieback site (0.248) than in intact forest (0.266),
                            so multiplying by it destroys the signal.

    Separate low shrubland from tall stands with `canopy_top` FIRST (AUC 0.734);
    this score is only meaningful within the tall class.
    """
    r, w = np.asarray(r, float), np.asarray(w, float)
    tot = max(float(w.sum()), 1e-9)
    w_can = max(float(w[r > canopy_cut].sum()), 1e-9)
    lo, hi = shrub_band
    w_low = float(w[(r >= lo) & (r < hi)].sum())
    w_up = float(w[r >= relict_cut].sum())
    if method == "lowlayer":
        return w_low / tot
    if method == "conjunction":
        return 2.0 * np.sqrt((w_low / w_can) * (w_up / w_can))
    raise ValueError("method must be 'lowlayer' or 'conjunction'")


def fit_cloud(x, y, z, number_of_returns, classification, res=10.0,
              ground_class=2, canopy_cut=1.5, min_ground=8, **fit_kw):
    """Fit every populated cell of a point cloud on a `res`-metre grid.

    Pass returns with overlap (12) and noise (7) ALREADY removed.  Returns a
    dict of numpy arrays, one entry per populated cell, carrying cell centres,
    the fitted parameters, canopy top, stratum covers and canopy skewness.
    """
    x, y, z = (np.asarray(v, float) for v in (x, y, z))
    w = 1.0 / np.asarray(number_of_returns, float)
    cl = np.asarray(classification)
    ix, iy = np.floor(x / res).astype(np.int64), np.floor(y / res).astype(np.int64)
    key = ix * np.int64(10 ** 9) + iy
    uniq, inv = np.unique(key, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    x, y, z, w, cl, inv = (a[order] for a in (x, y, z, w, cl, inv))
    n = len(uniq)
    starts = np.searchsorted(inv, np.arange(n))
    ends = np.searchsorted(inv, np.arange(n), side="right")
    cols = {k: np.full(n, np.nan) for k in
            ("easting", "northing", "ground_z", "pi_c", "sigma_g", "xi", "omega",
             "alpha", "dbic", "n_eff", "canopy_top", "cover", "low_cover",
             "mid_cover", "relict_cover", "skew")}
    cols["tree"] = np.zeros(n, bool)
    for i in range(n):
        s, e = starts[i], ends[i]
        cx = (uniq[i] // np.int64(10 ** 9)) * res + res / 2
        cy = (uniq[i] % np.int64(10 ** 9)) * res + res / 2
        cols["easting"][i], cols["northing"][i] = cx, cy
        gm = cl[s:e] == ground_class
        if gm.sum() < min_ground:
            continue
        plane = ground_plane(x[s:e][gm], y[s:e][gm], z[s:e][gm], cx, cy)
        r = heights_above_plane(x[s:e], y[s:e], z[s:e], plane, cx, cy)
        ww = w[s:e]
        cols["ground_z"][i] = plane[0]
        low, mid, relict = stratum_cover(r, ww)
        cols["low_cover"][i], cols["mid_cover"][i], cols["relict_cover"][i] = low, mid, relict
        cols["cover"][i] = float(ww[r > canopy_cut].sum()) / max(float(ww.sum()), 1e-9)
        cols["skew"][i] = canopy_skewness(r, ww, canopy_cut)
        f = fit_cell(r, ww, canopy_cut=canopy_cut, **fit_kw)
        if f is None:
            continue
        for k, v in (("pi_c", f.pi_c), ("sigma_g", f.sigma_g), ("xi", f.xi),
                     ("omega", f.omega), ("alpha", f.alpha), ("dbic", f.dbic),
                     ("n_eff", f.n_eff)):
            cols[k][i] = v
        cols["tree"][i] = bool(f.tree)
        cols["canopy_top"][i] = canopy_top(f)
    return cols
