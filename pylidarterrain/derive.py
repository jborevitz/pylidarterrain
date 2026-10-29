"""Slope, aspect, flow routing and the topographic wetness index."""
from __future__ import annotations
import numpy as np
from collections import deque

# 8 neighbour offsets, (drow, dcol); row increases SOUTHWARD on a north-up grid.
OFFSETS = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
_OPPOSITE = {k: 7 - k for k in range(8)}


def _rd(dem, res):
    import richdem as rd
    a = rd.rdarray(np.ascontiguousarray(dem, dtype=float), no_data=-9999)
    a.geotransform = (0.0, res, 0.0, 0.0, 0.0, -res)
    a.projection = ""
    return a


def fill_depressions(dem, res, epsilon=True):
    """Hydrologically condition the surface. `epsilon=True` imposes a minimal
    gradient across flats so flow routing has somewhere to go."""
    import richdem as rd
    return np.array(rd.FillDepressions(_rd(dem, res), epsilon=epsilon, in_place=False))


def slope_aspect(dem, res):
    """Slope in DEGREES and aspect in degrees clockwise from north."""
    import richdem as rd
    a = _rd(dem, res)
    return (np.array(rd.TerrainAttribute(a, attrib="slope_degrees")),
            np.array(rd.TerrainAttribute(a, attrib="aspect")))


def slope_percent(slope_deg):
    """Grade in percent. 100% is 45 degrees -- worth stating on every figure."""
    return 100.0 * np.tan(np.radians(slope_deg))


def flow_accumulation(dem, res, method="Dinf"):
    """Upslope contributing area in CELLS."""
    import richdem as rd
    return np.array(rd.FlowAccumulation(_rd(dem, res), method=method))


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
