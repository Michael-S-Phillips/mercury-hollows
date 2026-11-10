# g3_catalog: catalog counts, size distribution, Caloris, terrain classes, and workshop tests rerun on v3/v3b (6 Nov 2026)

This folder reruns the catalog-dependent numbers of Sections 2.1, 3.1, 3.2 (Caloris and terrain classes) and 4.4 of
`Hollow_Formation_Model_Paper_v6_3.tex`, and of Sections S1, S2, S4 and S7 of `supplementary-materials_v13.tex`, on the
screened catalogs of `code/rerun/out_v3` (v3: 368 spurious detections removed; v3b: also 232 uncertain ones). No existing file
was changed. Every number with its v2 record, v3 and v3b value is in `changes_g3.csv` (85 rows).

## How the scripts were copied
Each copied script carries a first-line comment naming its source. Only input paths, output paths and environment settings were
changed: `G3_CAT` = v2 | v3 | v3b selects the catalog, and `G3_THETA` = fixed | refit selects the detection model for the size
distribution. The `make_copies_*.py` scripts make the copies with exact string replacements (each asserted to occur once), so
the changes can be read there. Every copied chain was also run on v2. Each v2 run reproduces its existing record exactly,
except where noted below.

| Folder | Copied scripts (source) | New scripts |
|---|---|---|
| `inputs/` | | `prep_inputs.py`: pairs restricted to non-trace hollows of each catalog (the pairs column `h` is the row index of the v2 non-trace hollows, not `hollow_id`); containment-graph nodes with `n_hollows` and `hosts` recomputed from `host_assignments_{v}.csv` |
| `counts/` | | `catalog_counts.py`: hollows, fields (existing ids and rebuilt 20 km chains), pairs, size counts |
| `sfd/` | `refit_detection.py`, `fit_sfd.py`, `sfd_windows.py`, `sfd_numbers.py` (analysis_workshop/completeness_sfd); `sfd_fig_v1.py` (figscripts); `size_law.py` (code/rerun/closure) | `run_sfd.sh`; `sfd_derived.py` (Section 4.4 predictions, exponent range, pair shares, largest-detection excess); `refit_detection_resumable.py` (not needed in the end; the direct refit ran in 30 s) |
| `caloris/` | `caloris_v9.py`, `caloris_v9_controls.py` (analysis_workshop/caloris_v9), `caloris_test.py` (peak_wall_test) | |
| `terrain/` | `peak_wall.py`, `analyze.py` (peak_wall_test); `terrain_table.py` (analysis_v3) | `merge_rows.py` (joins peak_wall chunks); `terrain_obs_exp.py` (Table S2 observed/expected column and 6.1 density, for which no script was found); `run_terrain.sh` |
| `melt_source/positions/` | `positions_build.py`, `positions_analyze.py` | |
| `melt_source/knobs/` | `common.py`, `t2_knob_hollows.py`, `t2_peak_hollows.py`, `summarize.py` | |

The paper's Caloris numbers (179 craters, 0 of 144, 14 of 35, and so on) come from `analysis_workshop/caloris_v9/caloris_v9.py`.
The script named in the request, `peak_wall_test/caloris_test.py`, gives the older 0/128 and 14/32 counts, so both were rerun.
Outputs are in `*/out_{v2,v3,v3b}[_fixed|_refit]/`. `build_changes.py` writes `changes_g3.csv` from those outputs.

## Headline results
- **Catalog.** Hollows: 6,731 to 6,363 (v3) and 6,131 (v3b). Non-trace hollows: 6,667 to 6,299 and 6,067. Fields: 1,228 to 973 and 797
  (978 and 803 when the 20 km chains are rebuilt, because removals split a few fields). Pairs: 92,382 to 85,602 and 81,484.
- **Size distribution.** Removing spurious detections takes out large hollows: 185 to 96 (v3) and 73 (v3b) cataloged hollows
  larger than 5 km. The 0.3-3 km exponents rise from 1.4-1.6 to 1.5-1.6 (v3) and 1.6-1.7 (v3b). The exponent above 3 km rises
  from 2.68 to 3.29 (2.94-3.65) and 3.39. The corrected count larger than 5 km falls from 214 to 107 (86-131) and 82.
  The power law predicts 380-450 (v3), so the deficit of large hollows is stronger, not weaker. The steepening now starts nearer 2 km
  (2-6 km exponent 2.01 against 1.68). Section 4.4 becomes 5.2-5.8 times and 160,000-200,000 hollows larger than 0.1 km (v3), and
  "about 110 against 96" hollows larger than 5 km.
- **Detection model refit (sensitivity).** Refit on the v3 pairs, p_max rises from 0.243 to 0.272, and the upper turnover moves
  from 196 to 187 px with a wider decline. The exponents are then 1.43-1.55 (still 1.4-1.6), 3.26 above 3 km, and the corrected
  count larger than 5 km is 105. Completeness grids built from a refit model would rise about 12-14% everywhere. Other groups
  use the unchanged out_v2 grids.
- **Caloris.** The main test is unchanged: 0 of 144, 14 of 35 (p = 1.0e-11), 2 of 762, and 3 of 5 against 1 of 52. Seven shallow
  craters with one or two hollows become 6 (v3) and 4 (v3b). In S7, 14 of the 20 small basin-parent hosts lose their hollows,
  which the screen classified as bowls or walls. The nested comparison holds: 3 of 12 against 6 of 452 (p = 0.001). The outer-fifth
  comparison does not: 4 of 6 hosts (p = 0.41; controlled odds ratio 2.5, p = 0.30). v3b gives 3 of 4 (p = 0.32).
- **Terrain classes.** The host count falls from 352 to 321 (v3) and 262 (v3b). Peak/wall is 6.62 (4.25-10.13) and 6.92 (4.44-10.44),
  and floor/wall is 3.00 and 3.09. Every row of Table S3 moves by 0.1 or less (v3) or 0.4 or less (v3b).
- **Workshop tests.** The peak-position and floor-margin results are unchanged. The floor-margin result for fields (0.77, 0.58-0.99) becomes
  0.82 (0.61-1.08). For knobs of 5 km or more under v3b, the 3-degree block-bootstrap interval reaches 0.75.

## Conclusion changes (changes_g3.csv, conclusion_change = YES)
- S7 and the last sentence of the Caloris paragraph in 3.2: the claim that small Caloris craters with hollows lie mostly in the
  outer fifth of the basin is no longer significant under v3 or v3b. The nested-parent half of that sentence holds.
- Workshop only: the floor-margin interval for fields now includes 1. Under v3b, one of the three intervals for knobs of 5 km or more includes 1.
- No conclusion of 3.1 or 4.4 changes. "Large hollows are rarer than a power law predicts" holds and strengthens. The text should
  also say that screened crater bowls and walls are excluded, alongside the 64 detections larger than 10 km.

## Caveats
- Table S2 intervals for the isolated, basin-hosted and nested columns differ slightly in the supplement from `terrain_table_record.json`
  (for example 14.7-27.7 against 15.0-27.6). The v2 rerun reproduces the record, not the supplement. The observed/expected column
  and the 6.1 density were re-derived (v2: 2.99, 1.30, 0.45 and 6.2 against 2.98, 1.32, 0.45 and 6.1). Their v2-to-v3 change is
  what to use.
- "14-15%" is reproduced only when one detection per hollow is dropped at random (seed 1). The value depends slightly on the seed.
- `sfd_fig_v1.py` hard-codes the v2 exponent lines (1.4-1.6 and 2.68). The v3 figure shows new points with the old lines.
- Not rerun: the S9 color mask of the 6,731 hollows (color_mdr), the E-based hosting fits of S2 (exponent 1.7 instead of 1.6;
  S2 reports less than 5% change for exponents 1.3-1.9), and the knob color tests (t1_*; summarize.py reads their v2 outputs).
  `composition_v6/depth_lrm.py` reads no catalog file. The LRM share in the joint model belongs to the joint-model rerun.
- Background jobs do not survive a device shell call, so every step was run in the foreground in pieces of 3 minutes or less.
  `terrain/_scratch/` holds two superseded test chunks, which could not be deleted.

## Caloris tests on the v3r grid (9 Nov 2026)
`G3_GRID=v3r python3 caloris/caloris_v9.py` and `G3_GRID=v3r python3 caloris/caloris_v9_controls.py` read C1 from
g4_completeness/g1/v3r/crater_dag_nodes.csv (identical to inputs/crater_dag_nodes_v3.csv except C1) and write caloris/out_v3r
(logs out_v3r.log, out_v3r_controls.log). Without G3_GRID the scripts behave as before. Changes: 2 of 762 -> 2 of 840; 1 of 52 -> 1 of 58
(p = 0.001); 53 -> 56 craters, 0.14 -> 0.13; medians 0.32/0.40 -> 0.36/0.40; 0.86/0.35 -> 0.89/0.34; 5 of 44 -> 5 of 46; OR 69 -> 68.
Grid-independent counts (0 of 144, 14 of 35, 464, 3 of 12, 6 of 452) are unchanged.
