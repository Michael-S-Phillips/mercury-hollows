# g7_terrain_fix: closing the unclassified gaps in the terrain classes (6 Nov 2026)

Michael noticed that Figure S2 panels b and c leave a ring of ground between the peak ring and the wall in no class.

## Cause
`peak_wall_test/peak_wall.py` (and every copy, including the Figure S2 script) uses one floor level for the whole crater
(the 20th percentile of elevation at 0.15-0.6 R) and fixed radial limits: peak within 0.45 R (0.6 R for peak rings), wall
at 0.55 (0.7)-1.05 R, both more than 250 m above the floor level, and floor within 0.8 R and within 250 m of it. Ground
more than 250 m above the floor between the two limits falls in no class, and so does low ground beyond 0.8 R. In
peak-ring craters the annular floor between the ring and the wall stands 200-600 m above the inner floor
(`diag_profiles.py`, `diag_maps.png`), so it was either unclassified (0.6-0.7 R) or counted as wall (beyond 0.7 R).
In v3, 290 of the 3,760 hollows in 30-250 km hosts were unclassified, 180 of them in peak-ring craters.

## Fix (`terrain_classes.py`, used by `peak_wall_fix.py` and `figS2_terrain_examples_fix.py`)
Every pixel inside the rim gets one class. In 24 azimuth sectors, the radial profile of height sets a boundary radius
between the peak side and the wall side. In peak-ring craters it is the low point between the ring crest and 0.85 R, if
the ring stands at least 50 m above it and the wall rises at least 100 m beyond it; beyond the crest, heights are then
measured from that annular floor level. Where fewer than a quarter of the sectors show such a low point (for example
Scarlatti, 1 of 24), the boundary is 0.65 R and heights are measured from the inner floor. In other craters the boundary
is the low point of the profile at 0.35-0.65 R, clipped to 0.45-0.55 R. Boundaries are interpolated in azimuth. Peak and
wall are ground more than 250 m above the local floor level on their side of the boundary; everything else within 1.0 R is
floor. Unclassified hollows fall from 290 to 56 (those left lie 1.0-1.05 R from the center).

## Results (`changes_g7.csv`, `terrain_fix_comparison.csv`)
Configurations: `out_v3_gridv3` (primary: v3 catalog, v3 refit completeness grid), `out_v3_gridv2` and `out_v2_gridv2`
(main definition only, to separate the fix from the catalog and grid changes).
- Peaks to walls: 6.62 (4.25-10.13) with the old rules on v3 becomes 5.72 (3.74-8.57); on v2 the fix alone gives 5.72.
- Floors to walls: 3.00 to 2.64 (1.76-3.87).
- Equator-facing to pole-facing walls: 1.56 (0.94-2.48) becomes 1.74 (1.12-2.64), so the interval no longer includes 1.
- Peaks to walls in nested hosts: 2.08 (1.01-4.17) becomes 1.67 (0.80-3.32), so the interval now includes 1.
- Peak density in basin-hosted craters falls from 12.2 to 7.9 per 1,000 km2, level with their floors (7.8).
- Ring-crater peak/wall: 14.2 to 12.1. Hollows move from the ring's lower outer flank to the annular floor.

## Figures
`figS2_terrain_examples_v3cat_fix.png` (Figure S2 with the new classes) and `fig1_exposure_v9_v3cat_fix.png` (Figure 2,
panel d from the fixed terrain table).

## Text to change at v6.4
S4 class definitions, the Figure S2 caption (no unclassified ground inside the rim except 1.0-1.05 R), Section 3.2
(6.7 and 3.1 become 5.7 and 2.6), Table S2, Table S3, and the S8 wall-orientation sentence (1.52 becomes 1.74, 1.12-2.64).
