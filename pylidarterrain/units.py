"""Tree-sized scoring units from a canopy height model.

The unit of scoring, not a claim about individual trees.  The pipeline is

    chm            -> canopy height model on a raster
    cover_regime   -> per-cell {open, scattered, continuous}
    crowns         -> watershed segments inside the scattered regime
    cap_units      -> split any segment wider than a cap onto the analysis grid
    unit_returns   -> the returns belonging to each unit, for `stand.dieback_index`

Arrays in, arrays out.  Rasters are indexed [row, col] with row 0 at the
SOUTHERN edge (origin="lower"), matching `extent` in matplotlib; this is the
opposite of the north-up convention used by `dem.py`, so do not mix them
without flipping.

THE CAP IS A RULE ABOUT THE UNIT OF SCORING, not a delineation-accuracy
guideline.  A 40x40 m crown object is really many large trees, and each part
should carry its own score.  Seed separation cannot enforce it: seeds control
how many segments exist, not how far each one grows.  It must be a separate
step, and it is not free - see `cap_units`.
"""
import numpy as np
from scipy import ndimage as ndi

__all__ = ["chm", "cover_regime", "crowns", "cap_units", "unit_returns",
           "DISC"]

#: 4-connected structuring element.  Morphology is done on the BINARY canopy
#: mask, never on the height surface - smoothing heights merges adjacent crowns.
DISC = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)


def _grid_index(x, y, x0, y0, res, ny, nx):
    r = np.clip(((np.asarray(y) - y0) // res).astype(np.int64), 0, ny - 1)
    c = np.clip(((np.asarray(x) - x0) // res).astype(np.int64), 0, nx - 1)
    return r, c


def chm(x, y, height, x0, y0, res=1.5, shape=None, median_size=3):
    """Canopy height model: maximum height per cell, nearest-valid gap fill,
    then a median filter.

    CHOOSE `res` SO THAT MOST CELLS ALREADY CARRY A RETURN.  At 0.8 pulses/m2
    a 1.0 m grid leaves 40% of cells empty, and gap-filling that many cells
    smears canopy across the gaps - at Avonside a 1 m subcircle CHM inflated
    canopy cover from 27% to 61% against 38% measured at point level, and
    merged separate crowns.  At 1.5 m, 90% of cells were already populated.
    Check the populated fraction on any new capture before trusting the result.
    """
    x, y, height = (np.asarray(v, float) for v in (x, y, height))
    if shape is None:
        ny = int(np.ceil((y.max() - y0) / res)) + 1
        nx = int(np.ceil((x.max() - x0) / res)) + 1
    else:
        ny, nx = shape
    r, c = _grid_index(x, y, x0, y0, res, ny, nx)
    out = np.full(ny * nx, np.nan)
    order = np.argsort(height)                 # last write wins = maximum
    np.put(out, (r * nx + c)[order], height[order])
    out = out.reshape(ny, nx)
    filled = np.isfinite(out)
    if not filled.all():
        idx = ndi.distance_transform_edt(~filled, return_distances=False,
                                         return_indices=True)
        out = np.where(filled, out, out[tuple(idx)])
    if median_size:
        out = ndi.median_filter(out, size=median_size)
    return out, float(filled.mean())


def cover_regime(x, y, height, weight, x0, y0, cell=10.0, shape=None,
                 canopy_cut=1.5, cover_split=0.6, min_top=2.0, n_top=3):
    """Per-cell structural regime from canopy cover.

    Returns (regime, cover, top), regime coded 0 open / 1 scattered /
    2 continuous.  Cover is pulse weight above `canopy_cut` over total pulse
    weight.  Crown delineation is only meaningful in the scattered regime.

    DO NOT USE CONNECTED-PATCH SIZE FOR THIS.  Above roughly 40% cover the
    canopy percolates: at Avonside a single connected component of 26.6 ha held
    96% of all canopy area, so "a patch far larger than a crown" classified
    nearly the whole farm as continuous.  Per-cell cover does not percolate.

    `cover_split` IS A JUDGEMENT, NOT A NATURAL BREAK - cover is unimodal.  At
    Avonside 0.6 gave 12% continuous, 0.5 gave 25% and 0.7 gave 4%.  Report the
    sensitivity on any new dataset; it will not transfer unchanged.
    """
    x, y, height, weight = (np.asarray(v, float) for v in (x, y, height, weight))
    if shape is None:
        ny = int(np.ceil((y.max() - y0) / cell)) + 1
        nx = int(np.ceil((x.max() - x0) / cell)) + 1
    else:
        ny, nx = shape
    r, c = _grid_index(x, y, x0, y0, cell, ny, nx)
    flat = r * nx + c
    tot = np.bincount(flat, weights=weight, minlength=ny * nx)
    can = np.bincount(flat, weights=weight * (height > canopy_cut),
                      minlength=ny * nx)
    cover = np.where(tot > 0, can / np.maximum(tot, 1e-9), np.nan)
    order = np.lexsort((height, flat))
    fs, hs = flat[order], height[order]
    lo = np.searchsorted(fs, np.arange(ny * nx), "left")
    hi = np.searchsorted(fs, np.arange(ny * nx), "right")
    top = np.full(ny * nx, np.nan)
    ok = np.where(hi - lo >= n_top)[0]
    if ok.size:
        top[ok] = np.mean([hs[hi[ok] - k] for k in range(1, n_top + 1)], axis=0)
    cover, top = cover.reshape(ny, nx), top.reshape(ny, nx)
    regime = np.zeros((ny, nx), np.int8)
    has = (top > min_top) & np.isfinite(cover)
    regime[has & (cover < cover_split)] = 1
    regime[has & (cover >= cover_split)] = 2
    return regime, cover, top


def crowns(chm_raster, mask, res, min_separation=5.0, min_area=15.0,
           canopy_cut=2.0, closing=True, closing_mask=None):
    """Watershed crown segmentation inside `mask`.

    `mask` should be the scattered-regime area.  Canopy is thresholded at
    `canopy_cut`, closed as a BINARY mask, seeded with local maxima at least
    `min_separation` apart, and grown by watershed on the inverted CHM.
    Segments under `min_area` are absorbed into the neighbour with the longest
    shared border.  Returns an int label raster, 0 = background.

    CLOSE OVER THE WHOLE STUDY AREA, THEN RESTRICT.  Pass `closing_mask` as the
    full area of interest (e.g. the property boundary) when `mask` is only the
    scattered regime.  Closing inside the regime mask instead lets the regime
    boundary cut the structuring element, which moves crown edges and changes
    the segment count - at Avonside by about 0.2%, small but not nothing, and
    the difference is concentrated exactly at regime boundaries where crowns
    are most ambiguous.
    """
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed
    region = mask if closing_mask is None else closing_mask
    canopy = (chm_raster > canopy_cut) & region
    if closing:
        canopy = ndi.binary_closing(canopy, structure=DISC) & region
    canopy = canopy & mask
    surface = np.where(canopy, chm_raster, 0.0)
    sep = max(1, int(round(min_separation / res)))
    peaks = peak_local_max(surface, min_distance=sep, labels=canopy,
                           exclude_border=False)
    markers = np.zeros(surface.shape, dtype=np.int32)
    markers[tuple(peaks.T)] = np.arange(1, len(peaks) + 1)
    seg = watershed(-surface, markers, mask=canopy)
    return _absorb_small(seg, res, min_area)


def _absorb_small(seg, res, min_area):
    cell = res * res
    counts = np.bincount(seg.ravel())
    small = np.where((counts * cell < min_area) &
                     (np.arange(len(counts)) > 0))[0]
    for s in small:
        m = seg == s
        if not m.any():
            continue
        nb = seg[ndi.binary_dilation(m, structure=DISC) & ~m]
        nb = nb[nb > 0]
        if nb.size:
            seg[m] = np.bincount(nb).argmax()
    return _renumber(seg)


def _renumber(seg):
    vals = np.unique(seg[seg > 0])
    remap = np.zeros(int(seg.max()) + 1, np.int32)
    remap[vals] = np.arange(1, len(vals) + 1)
    return remap[seg]


def cap_units(seg, res, cap_diameter=10.0, grid=10.0, min_area=15.0):
    """Split every segment wider than `cap_diameter` onto a `grid` lattice.

    Any segment whose area exceeds a circle of `cap_diameter` is intersected
    with the grid; the resulting fragments under `min_area` are absorbed back
    into the largest neighbour.  Returns a new label raster.

    THE COST, WHICH MUST BE STATED WHEREVER THE CAP IS USED.  Splitting divides
    each object's returns among its parts.  At Avonside the median vegetation
    returns per object fell from 45 to 24, so the share of objects clearing a
    30-return profile floor dropped from about 70% to 40%.  The objects pushed
    below the floor are systematically the low, dense ones - which is the
    died-back form - so the cap trades scoring resolution for a larger
    unassessed fraction, biased against the signal of interest.
    """
    cap_area = np.pi * (cap_diameter / 2.0) ** 2
    cell = res * res
    counts = np.bincount(seg.ravel())
    ids = np.where(counts > 0)[0]
    ids = ids[ids > 0]
    big = ids[counts[ids] * cell > cap_area]
    rows, cols = np.mgrid[0:seg.shape[0], 0:seg.shape[1]]
    gy = ((rows * res) // grid).astype(np.int64)
    gx = ((cols * res) // grid).astype(np.int64)
    key = np.where(np.isin(seg, big),
                   seg.astype(np.int64) * 1_000_000 + gy * 1000 + gx,
                   seg.astype(np.int64))
    key[seg == 0] = 0
    _, unit = np.unique(key, return_inverse=True)
    unit = unit.reshape(seg.shape).astype(np.int32)
    unit[seg == 0] = 0
    return _absorb_small(_renumber(unit), res, min_area)


def unit_returns(unit, x, y, height, weight, x0, y0, res):
    """Group returns by the unit whose cell they fall in.

    Returns (unit_ids, r_by_unit, w_by_unit) ready for `stand.dieback_index`
    and `stand.profile_null`.  Returns outside every unit are dropped.
    """
    ny, nx = unit.shape
    r, c = _grid_index(x, y, x0, y0, res, ny, nx)
    uid = unit[r, c]
    keep = uid > 0
    uid, h, w = uid[keep], np.asarray(height, float)[keep], np.asarray(weight, float)[keep]
    order = np.argsort(uid, kind="stable")
    uid, h, w = uid[order], h[order], w[order]
    ids = np.unique(uid)
    bounds = np.searchsorted(uid, np.concatenate([ids, ids[-1:] + 1]))
    return (ids,
            [h[bounds[i]:bounds[i + 1]] for i in range(len(ids))],
            [w[bounds[i]:bounds[i + 1]] for i in range(len(ids))])
