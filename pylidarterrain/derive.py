"""Slope, aspect, flow routing and the topographic wetness index."""
from __future__ import annotations
import numpy as np
from collections import deque

# 8 neighbour offsets, (drow, dcol); row increases SOUTHWARD on a north-up grid.
OFFSETS = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
_OPPOSITE = {k: 7 - k for k in range(8)}


def fill_depressions(dem, res=None, epsilon=True, eps=1e-3):
    """Hydrologically condition the surface by PRIORITY FLOOD (Barnes, Lehman &
    Mulla 2014), so every cell has a descending path to the grid edge.

    `epsilon=True` raises each filled cell a hair above its spill point, which
    imposes a minimal gradient across flats so flow routing has somewhere to go.
    Input must be void-free; run `fill_voids` first.
    """
    import heapq
    d = np.asarray(dem, float).copy()
    if np.isnan(d).any():
        raise ValueError("fill_depressions needs a void-free DEM; call fill_voids first")
    ny, nx = d.shape
    closed = np.zeros((ny, nx), bool)
    heap = []
    for r in range(ny):
        for c in (0, nx - 1):
            heapq.heappush(heap, (d[r, c], r, c)); closed[r, c] = True
    for c in range(1, nx - 1):
        for r in (0, ny - 1):
            heapq.heappush(heap, (d[r, c], r, c)); closed[r, c] = True
    step = eps if epsilon else 0.0
    while heap:
        z, r, c = heapq.heappop(heap)
        for dr, dc in OFFSETS:
            rr, cc = r + dr, c + dc
            if 0 <= rr < ny and 0 <= cc < nx and not closed[rr, cc]:
                closed[rr, cc] = True
                if d[rr, cc] <= z:
                    d[rr, cc] = z + step
                heapq.heappush(heap, (d[rr, cc], rr, cc))
    return d


# Tarboton (1997) D-infinity facets: (cardinal neighbour, diagonal neighbour).
# Row index increases SOUTHWARD, so (-1, 0) is north.
_DINF_FACETS = [((0, 1), (-1, 1)), ((-1, 0), (-1, 1)), ((-1, 0), (-1, -1)),
                ((0, -1), (-1, -1)), ((0, -1), (1, -1)), ((1, 0), (1, -1)),
                ((1, 0), (1, 1)), ((0, 1), (1, 1))]


def dinf_proportions(dem, res):
    """D-infinity flow partition (Tarboton 1997).

    Flow leaves a cell at a continuous angle and is split between the two
    neighbours bounding that angle, which avoids the grid-aligned striping of
    D8 without the over-dispersion of multiple-flow-direction.  Returns
    (j1, p1, j2, p2) as flat receiver indices and proportions; j = -1 where a
    cell has no downslope neighbour.
    """
    z = np.asarray(dem, float)
    ny, nx = z.shape
    big = np.inf
    pad = np.pad(z, 1, constant_values=big)
    core = pad[1:-1, 1:-1]
    limit = np.arctan2(res, res)          # pi/4 on a square grid
    best_s = np.full((ny, nx), -np.inf)
    best_r = np.zeros((ny, nx))
    best_f = np.full((ny, nx), -1, dtype=np.int8)
    for f, ((dr1, dc1), (dr2, dc2)) in enumerate(_DINF_FACETS):
        e1 = pad[1 + dr1:1 + dr1 + ny, 1 + dc1:1 + dc1 + nx]
        e2 = pad[1 + dr2:1 + dr2 + ny, 1 + dc2:1 + dc2 + nx]
        s1 = (core - e1) / res
        s2 = (e1 - e2) / res
        r = np.arctan2(s2, s1)
        s = np.hypot(s1, s2)
        lo = r < 0
        r = np.where(lo, 0.0, r)
        s = np.where(lo, s1, s)
        hi = r > limit
        r = np.where(hi, limit, r)
        s = np.where(hi, (core - e2) / (res * np.sqrt(2.0)), s)
        upd = np.isfinite(s) & (s > best_s)
        best_s[upd], best_r[upd], best_f[upd] = s[upd], r[upd], f
    rows = np.repeat(np.arange(ny), nx).reshape(ny, nx)
    cols = np.tile(np.arange(nx), ny).reshape(ny, nx)
    j1 = np.full((ny, nx), -1, dtype=np.int64)
    j2 = np.full((ny, nx), -1, dtype=np.int64)
    p2 = best_r / limit
    p1 = 1.0 - p2
    for f, ((dr1, dc1), (dr2, dc2)) in enumerate(_DINF_FACETS):
        m = best_f == f
        if not m.any():
            continue
        for (dr, dc), jj in (((dr1, dc1), j1), ((dr2, dc2), j2)):
            rr, cc = rows[m] + dr, cols[m] + dc
            inside = (rr >= 0) & (rr < ny) & (cc >= 0) & (cc < nx)
            idx = np.full(rr.shape, -1, dtype=np.int64)
            idx[inside] = rr[inside] * nx + cc[inside]
            jj[m] = idx
    dead = (best_s <= 0) | (best_f < 0)
    j1[dead] = -1; j2[dead] = -1
    p1 = np.where(dead, 0.0, p1); p2 = np.where(dead, 0.0, p2)
    return j1.ravel(), p1.ravel(), j2.ravel(), p2.ravel()


def flow_accumulation(dem, res, method="dinf", p=1.1):
    """Upslope contributing area in CELLS.

    `method="dinf"` (default) is Tarboton D-infinity, the usual choice for a
    wetness index.
    `method="mfd"` is multiple-flow-direction (Freeman 1991): a cell divides its
    flow between ALL lower neighbours in proportion to slope**p.  It is the right
    default for a wetness index on a hillslope, where water spreads rather than
    following a single thread.  `method="d8"` sends everything to the steepest
    neighbour, which is better for extracting a channel network but draws
    stripey hillslopes.

    Expects a depression-filled surface.  Replaces an earlier richdem-backed
    D-infinity implementation; the package carries no GPL dependency.
    """
    d = np.asarray(dem, float)
    ny, nx = d.shape
    n = ny * nx
    flat = d.ravel()
    acc = np.ones(n)
    nb_off = np.array([dr * nx + dc for dr, dc in OFFSETS])
    length = np.array([np.hypot(dr, dc) * res for dr, dc in OFFSETS])
    rows = np.repeat(np.arange(ny), nx)
    cols = np.tile(np.arange(nx), ny)
    order = np.argsort(flat)[::-1]
    if method == "dinf":
        j1, p1, j2, p2 = dinf_proportions(d, res)
        for i in order:
            a = acc[i]
            if j1[i] >= 0 and p1[i] > 0:
                acc[j1[i]] += a * p1[i]
            if j2[i] >= 0 and p2[i] > 0:
                acc[j2[i]] += a * p2[i]
        return acc.reshape(ny, nx)
    d8 = method == "d8"
    for i in order:
        r, c = rows[i], cols[i]
        zi = flat[i]
        best, best_g, tot = -1, 0.0, 0.0
        ws = []
        for k, (dr, dc) in enumerate(OFFSETS):
            rr, cc = r + dr, c + dc
            if rr < 0 or rr >= ny or cc < 0 or cc >= nx:
                continue
            j = i + nb_off[k]
            g = (zi - flat[j]) / length[k]
            if g <= 0:
                continue
            if d8:
                if g > best_g:
                    best_g, best = g, j
            else:
                wk = g ** p
                ws.append((j, wk)); tot += wk
        if d8:
            if best >= 0:
                acc[best] += acc[i]
        elif tot > 0:
            a = acc[i]
            for j, wk in ws:
                acc[j] += a * wk / tot
    return acc.reshape(ny, nx)


def slope_aspect(dem, res):
    """Slope in DEGREES and aspect in degrees clockwise from north.

    Horn's (1981) 3x3 finite difference, the method behind ArcGIS and GDAL.
    Aspect is the direction a slope FACES, i.e. the downhill direction: 0 = north,
    90 = east.  On a north-up grid row index increases SOUTHWARD, which is the
    sign that most often gets flipped - `assert_north_up` and the synthetic-plane
    tests exist to catch exactly that.  Flat cells return aspect 0.
    """
    z = np.asarray(dem, float)
    p = np.pad(z, 1, mode="edge")
    a, b, c = p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:]      # row to the NORTH
    d, f = p[1:-1, :-2], p[1:-1, 2:]
    g, h, i = p[2:, :-2], p[2:, 1:-1], p[2:, 2:]         # row to the SOUTH
    dz_east = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * res)
    dz_south = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * res)
    dz_north = -dz_south
    slope = np.degrees(np.arctan(np.hypot(dz_east, dz_north)))
    aspect = np.degrees(np.arctan2(-dz_east, -dz_north)) % 360.0
    aspect[slope == 0] = 0.0
    return slope, aspect


def slope_percent(slope_deg):
    """Grade in percent. 100% is 45 degrees -- worth stating on every figure."""
    return 100.0 * np.tan(np.radians(slope_deg))


def flow_dirs(dem_filled, res):
    """Steepest-descent (D8) receiver index per cell, -1 where none."""
    d = np.asarray(dem_filled, float)
    ny, nx = d.shape
    pad = np.pad(d, 1, constant_values=np.inf)
    best = np.full((ny, nx), -1, dtype=np.int8)
    drop = np.zeros((ny, nx))
    for k, (dy, dx) in enumerate(OFFSETS):
        nb = pad[1+dy:1+dy+ny, 1+dx:1+dx+nx]
        length = np.hypot(dy, dx) * res
        g = (d - nb) / length
        upd = g > drop
        best[upd] = k; drop[upd] = g[upd]
    return best


def twi(flow_acc_cells, slope_deg, res, min_slope_deg=0.1):
    """ln(a / tan beta) with a = upslope area per unit contour width (m)."""
    a = np.maximum(flow_acc_cells, 1.0) * res
    tanb = np.tan(np.radians(np.maximum(slope_deg, min_slope_deg)))
    return np.log(a / tanb)


def upslope_mask(flow_dirs_arr, row, col):
    """All cells draining through (row, col), by breadth-first search upstream."""
    best = flow_dirs_arr
    ny, nx = best.shape
    m = np.zeros((ny, nx), bool); m[row, col] = True
    q = deque([(row, col)])
    while q:
        cy, cx = q.popleft()
        for k, (dy, dx) in enumerate(OFFSETS):
            uy, ux = cy + dy, cx + dx
            if 0 <= uy < ny and 0 <= ux < nx and not m[uy, ux] and best[uy, ux] == _OPPOSITE[k]:
                m[uy, ux] = True; q.append((uy, ux))
    return m
