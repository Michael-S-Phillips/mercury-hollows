# g9_draft_reruns: four draft analyses rerun on v3r + terrain fix (6 Nov 2026)

This folder reruns the four workshop analyses that `model_rewrite_draft_v1.tex` quotes, so the draft can cite current numbers.
The current configuration ("v3r + terrain fix") is the v3 catalog (`out_v3/hollows_canonical_v3.csv`, non-trace hollows, and
`host_assignments_v3.csv`), the v3 refit completeness grid (`g4_completeness/completeness_10ppd_v3.npz`, Dref = 1 km) and the g7
terrain classes (`g7_terrain_fix/terrain_classes.classify`). No existing file was changed, and nothing was written outside this folder.
All numbers are in `results.json`; it is written by `collect_results.py` from the per-configuration outputs.

## How the scripts were copied
`make_copies.py` makes every copy with exact string replacements, each asserted to occur once, so the changes can be read there.
The first line of each copy names its source and what changed. Environment settings select the configuration:
`G9_CAT` (v2 | v3), `G9_RULES` (fix | old), `G9_GRID` (v3 | v2) and, for the size-transition test, `G9_ST_CAT` (v2 | v3 | v3r).
Outputs go to `<analysis>/out_<CAT>_<RULES>_grid<GRID>/` (positions, slope), `knobs/out_<CAT>_grid<GRID>/` and `size_transition/out_<CAT>/`.

| Analysis | Copied from | Changes |
|---|---|---|
| Peak and floor positions | g3 `melt_source/positions/positions_build.py`, `positions_analyze.py` | class lines replaced by `classify` (az = degrees(arctan2(dy, dx)) % 360); grid; paths; radial bins gain 0.8-1.0 R; new measures (below) |
| Slope-controlled positions | workshop `positions/positions_slope_build.py`, `positions_slope.py` (g3 made no copy) | peak line replaced by `classify`; catalog, grid, paths; optional START STOP for the build (not needed: it runs in 5 s) |
| Knobs around Caloris | g3 `melt_source/knobs/common.py`, `t2_knob_hollows.py`, `t2_peak_hollows.py`, `summarize.py` | grid in `load_C`; output folder. Neither script uses terrain classes (`t2_peak_hollows` uses a geometric peak disk, rho <= 0.15 R), so `classify` was not applied |
| Size transition | workshop `size_transition/st_common.py`, `st1`-`st6` | crater table (v3r: `g4_completeness/g1/summary_model/v3r/crater_features.csv`, v3 hosts with C1 from the v3 grid), hollow catalog, output folder |

New scripts: `positions/merge_chunks.py` (joins build chunks; not needed in the end, since full builds take 20-25 s),
`size_transition/isolated_small_check.py` (counts isolated hosts 5-12 km and those resting only on crater-sized detections, with the
rule of `st8_isolated_small_hosts.py`), `collect_results.py`, and `_scratch/ring_outer_flank.py` (one caveat number).

**New measures in `positions_analyze.py` (`REC["g9_additions"]`).** Michael asked why the floor-margin test used an absolute distance
when host sizes differ. The workshop scripts and README give no reason for 3 km. The copy adds floor within 0.05 R and within 0.1 R
of wall- or peak-class cells (distance divided by the host radius) against farther floor, for hollows and fields, each also with the
crater x radius-bin null. It also adds floor within 0.2 R of the center against the rest of the floor as a split ratio. The draft's
"1.50 (1.10-1.93) times the density of the rest of the floor" is in fact the observed/expected ratio of the 0-0.2 R bin against the
whole floor; both are reported.

**d_nonfloor.** The copy keeps its definition (nearest unmasked non-floor cell within 1.05 R). Under the new classes the floor reaches
1.0 R, so the radial cut inside this measure moves from 0.8 R to the rim; non-floor cells are then peak, wall and the 1.0-1.05 R ring.
No requested number uses it; the margin tests use the distances to wall- and peak-class cells only.

## Checks of the copied chains (old configuration against the existing record)
- Positions: the build of hosts 0-40 (v3, old classes, v2 grid) matches the g3 cache array for array. The analysis run on the g3 v3
  cache and on the workshop v2 cache reproduces both `positions_record.json` files with 0 differences.
- Slope: the v2 build matches the workshop `slope_H250.npz`, and `positions_slope_record.json` is reproduced with 0 differences.
- Knobs: v3 with the v2 grid reproduces g3's `knob_hollows.csv`, `region_hollows.json` and peak chunk byte for byte, and `knob_tests_record.json` exactly.
- Size transition: v2 reproduces `size_transition_record.json` exactly (except the created/folder/inputs fields). The isolated-host check
  reproduces 29 isolated hosts 5-12 km, 16 of them crater-sized only, and 85 of 128 nested hosts.

## Results (95% intervals as the scripts compute them)
Columns: v2 record | v3, old classes, v2 grid (g3) | v3, old classes, v3 grid | **v3r + fix**. The third column separates the grid
change from the class fix; the grid alone moves no number by more than 0.01.

### 1. Peak and floor positions (crater-matched ratios, 2,000 crater bootstrap resamples)
| Measure | v2 | v3 old (g3) | v3 old, v3 grid | v3r + fix |
|---|---|---|---|---|
| Peak hollows / fields / craters | 933 / 127 / 106 | 920 / 118 / 99 | 920 / 118 / 99 | 810 / 116 / 99 |
| Lowest third vs upper two thirds, hollows | 0.61 (0.44-0.83) | 0.59 (0.42-0.80) | 0.59 (0.42-0.80) | 0.62 (0.39-0.90), 63 in band |
| Highest third vs lower two thirds, hollows | 1.29 (0.96-1.76) | 1.29 (0.98-1.81) | 1.29 (0.97-1.81) | 1.38 (1.11-1.75) |
| Lowest third, fields | 0.61 (0.18-1.13) | 0.46 (0.11-0.89) | 0.46 (0.11-0.89) | 0.59 (0.22-1.06) |
| Highest third, fields | 1.94 (1.43-2.67) | 2.15 (1.57-3.02) | 2.15 (1.57-3.02) | 1.61 (1.15-2.29) |
| Peak within 2 km of floor vs rest of peak | 0.79 (0.58-1.04) | 0.79 (0.56-1.05) | 0.79 (0.56-1.05) | 0.73 (0.55-0.91) |
| Floor hollows / fields / craters | 1,796 / 287 / 223 | 1,727 / 243 / 197 | 1,727 / 243 / 197 | 2,063 / 276 / 213 |
| Floor within 0.2 R: O/E vs whole floor (draft's number) | 1.50 (1.10-1.93) | 1.51 (1.10-2.00) | 1.52 (1.10-2.01) | 1.54 (1.15-2.03) |
| Floor within 0.2 R vs rest of floor, hollows | 1.59 (1.12-2.15), 272 | 1.60 (1.12-2.25), 265 | 1.61 (1.12-2.26) | 1.62 (1.17-2.25), 265 |
| Same, fields | 2.43 (1.78-3.23) | 2.62 (1.89-3.61) | 2.62 (1.89-3.62) | 2.50 (1.84-3.38) |
| Floor within 3 km of wall or peak, hollows | 0.90 (0.73-1.09) | 0.91 (0.74-1.10) | 0.91 (0.74-1.10) | 0.92 (0.76-1.11) |
| Same, fields | 0.77 (0.58-0.99) | 0.82 (0.61-1.08) | 0.82 (0.61-1.08) | 0.88 (0.66-1.14) |
| Floor within 0.05 R of wall or peak, hollows | 1.13 (0.84-1.45) | 1.16 (0.85-1.47) | 1.16 (0.85-1.47) | 1.11 (0.88-1.34) |
| Same, fields | 0.89 (0.63-1.18) | 1.05 (0.74-1.41) | 1.05 (0.74-1.41) | 1.12 (0.81-1.48) |
| Floor within 0.1 R of wall or peak, hollows | 1.20 (0.89-1.61) | 1.22 (0.88-1.65) | 1.22 (0.88-1.66) | 1.20 (0.89-1.54) |
| Same, fields | 0.88 (0.67-1.13) | 0.96 (0.72-1.24) | 0.96 (0.71-1.24) | 1.03 (0.78-1.32) |
| Radius-stratified null, hollows: 3 km / 0.05 R / 0.1 R | 0.95 / 1.08 / 1.15 | 0.96 / 1.11 / 1.17 | 0.96 / 1.11 / 1.18 | 1.00 (0.87-1.16) / 1.14 (0.94-1.32) / 1.26 (0.98-1.58) |

With the new classes, floor within 3 km of the wall alone is 0.63 (0.47-0.82) and within 3 km of the peak alone 1.59 (1.20-2.11).
Floor density by radius (O/E): 1.54, 1.02, 1.15, 0.83 and 0.37 (0.21-0.62) for 0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8 and 0.8-1.0 R.
The 0.8-1.0 R bin is new: the old floor class stopped at 0.8 R. Mean height quantile of peak hollows: 0.57 (0.53-0.60).
**Reading.** With the class fix, peak hollows fall by 110 (peak-ring hosts 711 to 593) and floor hollows rise by 336, which come from the outer
flanks of peak rings, from ground the old classes counted as wall, and from formerly unclassified ground.
The low-flank deficit holds, and the high-flank excess for hollows now excludes 1. No margin measure, absolute or radius-scaled, shows
an excess at the floor margins once radius is stratified; the radius-scaled bands lean above 1 because in large hosts they reach farther
onto the central floor beside the peak, where density is highest. The wall foot itself is depleted.

### 2. Slope-controlled positions (Mantel-Haenszel, crater x slope bin, 1-cell slope)
| Measure | v2 | v3 old, v2 grid | v3 old, v3 grid | v3r + fix |
|---|---|---|---|---|
| Hosts passing: loose (15, 500 m) | 158 (74 with peak hollows) | 145 (70) | 145 (70) | 150 (73) |
| Main (30, 750 m): hosts; peak hollows; fields | 89 (50); 675; 65 | 84 (46); 666; 59 | 84 (46); 666; 59 | 88 (48); 591; 61 |
| Strict (60, 1,000 m) | 43 (23) | 41 (22) | 41 (22) | 45 (24) |
| eta-slope r, within crater (main) | -0.014 | -0.011 | -0.011 | -0.014 |
| Lowest third vs upper two thirds, hollows: main | 0.57 (0.40-0.80) | 0.54 (0.37-0.76) | 0.55 (0.37-0.76) | 0.54 (0.34-0.80) |
| loose / strict | 0.57 (0.41-0.79) / 0.46 (0.32-0.63) | 0.54 / 0.45 | 0.55 / 0.45 | 0.54 (0.35-0.81) / 0.47 (0.30-0.71) |
| Highest third vs lower two thirds, fields: main | 2.15 (1.31-3.73) | 2.68 (1.56-4.82) | 2.68 (1.56-4.83) | 2.08 (1.25-3.70) |
| loose / strict | 2.26 (1.51-3.54) / 2.16 (1.01-4.72) | 2.62 / 2.80 | 2.63 / 2.80 | 2.08 (1.34-3.41) / 1.92 (0.81-4.44) |
| Highest third, hollows: main | 1.32 (0.91-2.12) | 1.32 (0.91-2.19) | 1.32 (0.90-2.18) | 1.39 (1.16-1.79) |

g3 did not rerun the slope test, so the v3 old-class columns were computed here.

**Why 750 m.** Peak cells are ground more than H = 250 m above the floor level, so the normalized height eta = h/(z_peak95 - z_floor)
cannot fall below 250/(z_peak95 - z_floor). The lowest third of eta (eta < 1/3) therefore exists only on peaks with more than
3 x 250 = 750 m of relief (positions README, caveats). README_slope adds that small peaks are the most smoothed by the 665 m DEM;
the 30-cell part of the screen addresses that. `positions_slope.py` gives no other reason. The g7 classes keep h > 250 m for peak
cells, so the bound still holds.

### 3. Knobs around Caloris (completeness-weighted densities; knob bootstrap, B = 4,000)
| Measure | v2 | v3, v2 grid (g3) | v3r |
|---|---|---|---|
| Knobs with hollows; hollows; fields | 20; 25; 16 | 19; 24; 15 | 19; 24; 15 |
| Knob vs terrain 1.5-3 knob radii away | 4.04 (1.99-7.72) | 3.98 (1.91-7.77) | 3.98 (1.91-7.76) |
| block bootstrap / exact interval | 1.83-12.4 / 2.33-6.86 | 1.82-12.8 / 2.28-6.83 | 1.82-12.8 / 2.28-6.83 |
| Knob vs within 0.15 R of Kinczyk Class 1-2 centers | 2.35 (1.04-6.51) | 2.37 (1.01-7.29) | 2.35 (1.005-7.24) |
The new grid raises completeness by about 10% on both terrains, so densities fall (knobs 5.04 to 4.59 per 1,000 km2 from the g3 run to v3r) and the ratios do not change.

### 4. Size transition (logistic, spatial-block cluster-robust intervals)
| Measure | v2 | v3, v2 grid | v3r |
|---|---|---|---|
| Isolated craters 5-60 km; hosts | 15,292; 117 | 15,292; 84 | 15,292; 84 |
| Step at 12 km, odds ratio | 0.51 (0.22-1.22) | 0.73 (0.22-2.43) | 0.73 (0.22-2.43); p = 0.60 |
| Calibrated p for any rise over 8-20 km | 0.87 | 0.82 | 0.79 |
| Craters 12-60 km, central peak vs none | 0.79 (0.43-1.47); 2,951 craters, 162 hosts | 0.73 (0.35-1.50); 135 hosts | 0.73 (0.36-1.50); 135 hosts |
| Isolated hosts 5-12 km (crater-sized detections only) | 29 (16) | 13 (2) | 13 (2) |
| Nested hosts 5-12 km (crater-sized only) | 128 (85) | 47 (11) | 47 (11) |
v3r variants of the 12 km step: spline 0.66 (0.19-2.34), logit completeness 0.74 (0.22-2.48), 5-100 km 0.93 (0.28-3.08),
10-60 km 0.63 (0.19-2.14), C1 >= 0.1 0.53 (0.15-1.96). The best step is a drop at 12.75 km (Delta AIC -0.26). The two isolated hosts
that still rest on crater-sized detections (603224, 6.6 km; 103519, 10.6 km) were screened by eye as hollows (class H).

## Caveats
- The 12 km step interval now reaches 2.43, against 1.22 in v2. The test still finds no rise, but it can exclude only a rise of more than
  about 2.4 times, not the 20% the v2 README stated. The "pending rerun" note can go: v3 flags the crater-sized detections.
- In the 35 peak-ring hosts whose annular floor is resolved, eta is still measured from the inner floor level. There 91 of 411 peak
  hollows and 26% of peak cells lie beyond the ring crest, where the local base is the annular floor, so the base of the outer flank sits
  at mid eta. Peak-ring hosts carry 593 of the 810 peak hollows.
- At 30 km, 0.05 R is 0.75 km, about one DEM cell, so in the smallest hosts the 0.05 R band holds only cells next to the boundary.
- The knob-to-Class 1-2 lower bound is 1.005; it rounds to 1.0 as in the draft but sits on the edge.
- The radius-stratified 0.1 R floor-margin ratio, 1.26 (0.98-1.58), is near significance. It reflects floor beside the peak, not the wall foot.
- Configurations were run in the foreground; nothing ran in the background. `_scratch/test_chunks/` holds two 40-host test builds and a
  5-simulation timing run that are not used.

## Rerun
```
export PYTHONDONTWRITEBYTECODE=1 G9_CAT=v3 G9_RULES=fix G9_GRID=v3
cd positions; for H in 250 150 400; do python3 positions_build.py 0 354 $H; done; python3 positions_analyze.py
python3 positions_slope_build.py; python3 positions_slope.py
cd ../knobs; python3 t2_knob_hollows.py; python3 t2_peak_hollows.py 0 4824; python3 summarize.py
cd ../size_transition; export G9_ST_CAT=v3r; for s in st1_bins st2_breaks st2c_blurred st3_nested st4_morph st5_caloris; do python3 $s.py; done
python3 st2b_null.py isolated 250 11; python3 st2b_null.py isolated 250 12; python3 st2b_null.py nested 250 21; python3 st2b_null.py nested 250 22
python3 st6_record.py; python3 isolated_small_check.py
cd ..; python3 collect_results.py
```
For the checks, run the analyses with G9_RULES=old G9_GRID=v2 (and G9_CAT=v2, G9_ST_CAT=v2); positions_analyze.py reads another cache through G9_CACHE_DIR.
