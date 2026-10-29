"""Grid classified LiDAR ground returns into a bare-earth DEM."""
from __future__ import annotations
from dataclasses import dataclass
import glob
import numpy as np
from scipy import ndimage

GROUND_CLASS = 2


@dataclass(frozen=True)
class Grid:
    """A north-up raster definition. Row 0 is the NORTHERN edge."""
    e0: float; n0: float; e1: float; n1: float; res: float

    @property
    def shape(self):
        return int(round((self.n1 - self.n0) / self.res)), int(round((self.e1 - self.e0) / self.res))

    @property
    def transform(self):
        """GDAL-order geotransform, negative pixel height (north-up)."""
        return (self.e0, self.res, 0.0, self.n1, 0.0, -self.res)

    def rowcol(self, easting, northing):
        return ((self.n1 - np.asarray(northing)) / self.res).astype(int), \
               ((np.asarray(easting) - self.e0) / self.res).astype(int)

    def centres(self, rows, cols):
        return self.e0 + (np.asarray(cols) + 0.5) * self.res, \
               self.n1 - (np.asarray(rows) + 0.5) * self.res


def build_dem(laz_paths, grid: Grid, ground_class: int = GROUND_CLASS, reducer: str = "mean"):
    """Bin ground returns onto `grid`. Returns (dem, count); dem is NaN where empty.

    `reducer` is "mean" (default, best for bare-earth LiDAR) or "min" (a
    bare-ground proxy for photogrammetric surface clouds; still sits on canopy
    under trees).
    """
    import laspy
    if isinstance(laz_paths, str):
        laz_paths = sorted(glob.glob(laz_paths))
    ny, nx = grid.shape
    acc = np.zeros((ny, nx)); cnt = np.zeros((ny, nx), dtype=np.int32)
    mn = np.full((ny, nx), np.inf)
    for p in laz_paths:
        las = laspy.read(p)
        m = np.asarray(las.classification) == ground_class if ground_class is not None \
            else np.ones(len(las.x), bool)
        x = np.asarray(las.x)[m]; y = np.asarray(las.y)[m]; z = np.asarray(las.z)[m]
        iy, ix = grid.rowcol(x, y)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        np.add.at(acc, (iy[ok], ix[ok]), z[ok])
        np.add.at(cnt, (iy[ok], ix[ok]), 1)
        np.minimum.at(mn, (iy[ok], ix[ok]), z[ok])
        del las, x, y, z
    if reducer == "min":
        dem = np.where(cnt > 0, mn, np.nan)
    else:
        dem = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    return dem, cnt


def fill_voids(dem, max_iter: int = 25):
    """Fill NaN cells by iterative 3x3 neighbourhood mean."""
    d = dem.copy()
    k = np.ones((3, 3))
    for _ in range(max_iter):
        nan = np.isnan(d)
        if not nan.any():
            break
        f = np.where(nan, 0.0, d); w = (~nan).astype(float)
        sm = ndimage.convolve(f, k, mode="nearest")
        wt = ndimage.convolve(w, k, mode="nearest")
        d = np.where(nan & (wt > 0), sm / np.maximum(wt, 1e-9), d)
    return d
