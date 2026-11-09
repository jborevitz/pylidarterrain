# pylidarterrain — methods

Conventions, traps and limitations. Findings from the block this was developed on are in the [README](../README.md).


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

## Dividing a block into parts

On a long narrow property, northings and eastings are the wrong frame: every
boundary that matters runs across the long axis or along it. Fit the axis first
(principal component of the boundary vertices), then work in distance along it `s`
and perpendicular offset `p`. A boundary that is "roughly east-west across the
block" is then a line of constant `s - m·p`, and the single number `m` is its
obliquity — easy to fit, easy to adjust by hand, and easy to state.

Three ways to place such a line, in increasing order of trust:

- **A break in a terrain statistic.** Bin grade against `s` and look for the step.
  Reliable where the break is sharp; it gave 600 m here against a hand-placed 580.
- **A drainage line.** Take it from the flow grid, not from the DEM by eye.
- **A crest.** Fit the ridge as a line through the per-bin maxima of elevation and
  report the rms residual: 87 m here, against a block 1914 m long, which is honest
  about a ridge that is crisp at both ends and broad in the middle.

Report the residual and let the owner move the line. A boundary that the person
walking the ground disagrees with is wrong regardless of what the algorithm says,
and the parts exist to be managed, not to be optimal.

## Grade, not degrees — but label both

Report slope as a percentage and put degrees on the colourbar. 100% is 45°, and
quoting one unit without the other invites exactly the confusion the second axis is
there to prevent: a "50% slope" sounds like half of vertical and is 26.6°.

## Contours on a 3 m surface need smoothing; the lines do not need moving

Return-density noise puts kinks in a contour that no walker would recognise. A 4.5 m
Gaussian smooth before contouring removes them without displacing a line anywhere
visible on the ground. Choose the interval from the terrain rather than by habit: at
a median grade of 29% a 10 m interval puts lines about 35 m apart, legible on the
steep faces where they crowd to 25 m and still informative on the benches.
