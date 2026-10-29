"""Guards against the orientation error that silently mirrors aspect.

Binning point data the natural way -- (y - y_min)/res -- puts row 0 at the
SOUTH edge.  GIS rasters put row 0 at the NORTH edge.  Mixing the two mirrors
the grid north-south, which inverts aspect.  In the Northern Hemisphere the
result still "looks right" (sun-facing reads as south-facing), so the error
survives review; in the Southern Hemisphere it inverts every insolation,
temperature and evapotranspiration conclusion drawn from aspect.
"""
import numpy as np


def assert_north_up(dem, aspect, res=1.0, tol=0.02):
    """Raise if `aspect` disagrees with the north-south gradient of `dem`.

    With row 0 at the north edge, the row index increases SOUTHWARD, so a
    positive d(elev)/d(row) means the ground rises to the south -- a
    north-facing slope.  No such cell may be labelled south-facing.
    """
    gy, _ = np.gradient(np.asarray(dem, float), res, res)
    north_facing = gy > 0
    if not north_facing.any():
        return True
    a = np.asarray(aspect)
    south_labelled = (a > 135) & (a < 225)
    bad = south_labelled[north_facing].mean()
    if bad > tol:
        raise AssertionError(
            f"{bad:.1%} of north-facing cells are labelled south-facing: "
            "the DEM is probably south-up while the geotransform says north-up."
        )
    return True
