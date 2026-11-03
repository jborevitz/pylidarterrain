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

## What it found at 845 Avonside

Two 2018 Geoscience Australia ELVIS tiles (ICSM Level 3, 0.44 ground points/m²,
0.3 m vertical accuracy at 95% CI). 3.32 M ground returns gridded at **3 m** over the
full 4 × 2 km tile union: 887,778 cells, a median of 4 returns per cell, 3.7% voids.
The block itself falls 309 m, from 887.9 to 1197.0 m AHD.

![Elevation with 10 m contours](docs/elevation_contours_3m.png)

*The land itself: hillshaded elevation with contours every 10 m, labelled every 50 m.
Contours are drawn on the 3 m surface after a 4.5 m Gaussian smooth, which suppresses
return-density noise without moving the lines anywhere a walker would notice. The
contour interval is the honest one for this terrain — at a median grade of 29% the
10 m lines sit about 35 m apart on the ground, close enough to read the shape of a
spur and far enough apart to stay legible on the steep faces, where they crowd to
25 m or less.*

![Grade and wetness as a bivariate map](docs/slope_twi_bivariate.png)

*Grade and the topographic wetness index on one map: red is steep, blue is wet, dark
is both, pale is flat and dry. Drawn across the whole tile with the property boundary
and the 88.7 ha creek bowl catchment overlaid, rather than clipped to the title — a
boundary is a legal object, not a hydrological one, and three quarters of the bowl
lies outside it.*

*The two variables are worth combining because they are largely independent here:
grade and TWI correlate at only −0.35, since TWI is driven by upslope area (r = 0.92)
rather than by the slope term in its own denominator (r = 0.36). The dark threads
where steep meets wet are the gully walls above active channels.*

![Grade at 3 m](docs/grade_3m.png)

*Grade alone: median 29%, reaching 100% — a true 45° — on those gully walls. The 10 m
grid capped at 36.8°, so the steepest ground is not new terrain but terrain the
coarser cells averaged away. The same mechanism drops the shaded annual shortwave
minimum from 4841 MJ/m² at 10 m to 3477 at 3 m.*

Grade is labelled in percent with degrees on the colourbar. 100% is 45°, and
reporting one unit without the other invites exactly the confusion the colourbar is
there to prevent.

### The block divides into four parts along its own terrain

![The four parts](docs/farm_four_parts.png)

*The property is a narrow rectangle, 1914 m on a bearing of 108°, so the useful
coordinate is distance along that axis rather than easting. Three terrain features
cut it into four parts, each boundary a line in that axis frame.*

| part | area | elevation | median grade | % north-facing | potential ET |
|---|---|---|---|---|---|
| Top farm | 11.6 ha | 1114–1197 m | 14% | 22% | 899 mm/yr |
| Upper forest | 10.2 ha | 962–1169 m | 48% | 1.5% | 653 mm/yr |
| Northerly woodland | 11.8 ha | 958–1051 m | 30% | 85% | 1018 mm/yr |
| Lower forest flats | 12.7 ha | 888–1051 m | 29% | 5.5% | 809 mm/yr |

The parts came out of the terrain nearly even in area, 10.2 to 12.7 ha, which nobody
chose. What they are not even in is demand: **365 mm/yr separates the northerly
woodland from the upper forest, 68% of annual rainfall, and the two share a
boundary.** The slope break between the top farm and the upper forest is the sharpest
line on the block, 29% grade to 53% across it. The southern creek enters from the
neighbour's 28 ha catchment and runs diagonally northeast to the pond junction. The
eastern crest is the real divide — east of it water leaves toward the flats, west of
it toward the creek — and it is the least certain of the three, since the ridge is
crisp at both ends and broad in the middle.

## Methods

The conventions this package enforces, the traps it exists to prevent, and what it
deliberately does not do are in **[`docs/methods.md`](docs/methods.md)**.

Methods and results are kept apart on purpose. The methods travel to any site and
will be re-checked against new ground and simulated data; the findings above are one
47 ha property in the Snowy Monaro, some of them specific to it and some likely to
generalise. Either can then be argued with without disturbing the other.
