"""Terrain horizon angles and the sky view factor."""
from __future__ import annotations
import numpy as np


def _shift(a, dy, dx):
    ny, nx = a.shape
    out = np.full_like(a, np.nan)
    ys = slice(max(dy, 0), ny + min(dy, 0)); yd = slice(max(-dy, 0), ny + min(-dy, 0))
    xs = slice(max(dx, 0), nx + min(dx, 0)); xd = slice(max(-dx, 0), nx + min(-dx, 0))
    if (ys.stop - ys.start) <= 0 or (xs.stop - xs.start) <= 0:
        return out
    out[yd, xd] = a[ys, xs]
    return out


def horizon_angles(dem, res, n_azimuth=32, max_distance=5000.0, n_steps=52):
    """Horizon elevation angle (degrees) per azimuth, shape (n_azimuth, ny, nx).

    Azimuth 0 is north, increasing clockwise.  The search is clipped to the
    grid, so a small tile truncates distant horizons -- `reach` reports how far
    each azimuth actually looked.
    """
    d = np.asarray(dem, float); ny, nx = d.shape
    max_cells = int(max_distance / res)
    radii = np.unique(np.round(np.geomspace(1, max(max_cells, 2), n_steps)).astype(int))
    hor = np.zeros((n_azimuth, ny, nx), dtype=np.float32)
    reach = np.zeros(n_azimuth)
    for k in range(n_azimuth):
        ar = np.radians(k * 360.0 / n_azimuth)
        ux, uy = np.sin(ar), np.cos(ar)
        best = np.zeros((ny, nx), dtype=np.float32); far = 0.0
        for r in radii:
            dx = int(round(ux * r)); dy = int(round(-uy * r))   # row grows southward
            if (dx == 0 and dy == 0) or abs(dy) >= ny or abs(dx) >= nx:
                continue
            dist = np.hypot(dx, dy) * res
            ang = (_shift(d, dy, dx) - d) / dist
            np.maximum(best, np.nan_to_num(ang, nan=-9.0), out=best)
            far = max(far, dist)
        hor[k] = np.degrees(np.arctan(best)); reach[k] = far
    return hor, reach


def sky_view_factor(horizon):
    """Isotropic sky view factor from horizon angles: mean of cos^2(horizon)."""
    return np.mean(np.cos(np.radians(np.clip(horizon, 0, None))) ** 2, axis=0)
