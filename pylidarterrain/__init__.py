"""Bare-earth DEMs and terrain derivatives from classified airborne LiDAR.

Built for ICSM Level-3 classified point clouds (e.g. Geoscience Australia ELVIS
tiles), where ground returns carry classification 2.  Everything downstream of
`build_dem` works on plain numpy arrays in **north-up** order: row 0 is the
NORTHERN edge, matching the GIS/GDAL convention and a negative pixel height in
the geotransform.
"""
from .dem import build_dem, fill_voids, Grid
from .derive import slope_aspect, flow_dirs, flow_accumulation, twi, fill_depressions
from .horizon import horizon_angles, sky_view_factor
from .checks import assert_north_up

__all__ = ["build_dem","fill_voids","Grid","slope_aspect","flow_dirs",
           "flow_accumulation","twi","fill_depressions","horizon_angles",
           "sky_view_factor","assert_north_up"]
__version__ = "0.1.0"
