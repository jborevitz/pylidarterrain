# pylidarterrain

Bare-earth DEMs and terrain derivatives from classified airborne LiDAR.

Written for ICSM Level-3 classified clouds (Geoscience Australia ELVIS tiles and
equivalents), where ground returns carry classification 2.

```python
from pylidarterrain import Grid, build_dem, fill_voids, slope_aspect, horizon_angles

grid = Grid(e0=652000, n0=5966000, e1=656000, n1=5968000, res=3.0)   # EPSG:28355
dem, counts = build_dem("tiles/*.laz", grid)
dem = fill_voids(dem)
slope, aspect = slope_aspect(dem, grid.res)
hor, reach = horizon_angles(dem, grid.res, n_azimuth=32, max_distance=5000)
```

## The orientation trap

`assert_north_up(dem, aspect)` exists because of a real error. Binning points the
natural way — `(y - y_min)/res` — puts row 0 at the **south** edge; GIS rasters put
it at the **north** edge. Mixing the two mirrors the grid and inverts aspect.

In the Northern Hemisphere the mistake hides, because sun-facing ground still reads
as south-facing. In the Southern Hemisphere it inverts every insolation, temperature
and evapotranspiration conclusion you draw. Call the guard after computing aspect.

## Resolution and what the data supports

Ground-return density is the binding constraint, not the gridding. At the ~0.44
ground points/m² of a standard wide-area survey a 3 m cell holds about four returns
and roughly 4% voids; 1 m holds less than one and is mostly interpolation. Check
`counts` before trusting a fine grid.

## Worked example — 47 ha in the Snowy Monaro, NSW

Built from two 2018 Geoscience Australia ELVIS tiles (ICSM Level 3, 0.44 ground
points/m², 0.3 m vertical accuracy at 95% CI). 3.32 M ground returns, a 3 m grid with
a median of 4 returns per cell and 3.7% voids.

![Elevation and topographic wetness](docs/farm_terrain_wetness.png)

*Bare-earth elevation (888–1196 m AHD) and the topographic wetness index. The wettest
decile picks out 4.6 ha of drainage hollows.*

![Grade in percent](docs/farm_grade_percent.png)

*Grade in percent with degrees on the colourbar. This block reaches 100% — a true 45°.
Reporting slope in one unit without the other invites exactly the confusion the
colourbar is labelled to prevent. Note also that the 10 m grid capped at 36.8°: the
steep ground only appears once the cells are small enough to resolve the gully walls.*

![Grade at 3 m](docs/grade_3m.png)

*Grade across the block at 3 m: median 29%, reaching 100% on the gully walls. The 10 m
grid capped at 36.8 degrees -- the steepest ground only appears once cells are small
enough to resolve the walls themselves.*
