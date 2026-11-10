# out_v3: catalog with spurious detections removed (6 Nov 2026)

## Why
The v2 catalog flags crater traces only when D_eq > 10 km. Visual screening showed that many large detections are
crater bowls, sunlit walls, or other non-hollows. Many small "host" craters (5-12 km) rested on such a detection.

## Screen (bowl_candidates.csv -> screen_set.csv -> screen_classes.csv)
EDR crops (qc/tiles, qc/edr_<set>_<k>.png) were classified H (hollow), B (crater bowl), W (wall, rim or slope),
O (other) or U (uncertain).
- **A:** D_eq >= 0.35 D, centroid within 0.5 R of a Herrick crater. 61 detections: 56 B, 2 W, 3 H.
- **B:** other detections with D_eq >= 1.5 km in fields of <= 3 hollows. 629 detections: 190 B, 94 W, 16 O, 219 U, 110 H.
  Classified by three subagents.
- **C:** D_eq 1.0-1.5 km singletons whose host has <= 2 hollows. 30 detections: 5 B, 2 W, 8 U, 15 H.
- **D (control):** 48 random detections with D_eq >= 1.5 km in fields of >= 10 hollows: 2 B, 1 W, 5 U, 40 H.
  This gives a spurious rate of ~6% in large fields; these were not screened beyond the sample.

## Inputs for reruns
- `hollows_canonical_v3.csv`: v2 columns. is_crater_trace = v2 trace OR class in {B, W, O}. 6,299 non-trace hollows (v2: 6,667).
- `hollows_canonical_v3b.csv`: also excludes U. 6,067 non-trace hollows.
- Both files keep is_crater_trace_v2, screen_class, is_spurious_v3 and is_uncertain_v3.
- `host_assignments_v3.csv` and `_v3b.csv`: GRL_insolation_hollows/data/host_assignments.csv with the excluded hollows removed.
  Host assignment is per hollow, so filtering equals reassigning.
  - Any-hollow hosts: v2 707, v3 534, v3b 441.
  - h09 hosts: v2 608, v3 448, v3b 371.

## Not changed
- Completeness grids (code/rerun/out_v2/completeness_10ppd.npz).
- Field ids (field_id). Fields left with no hollows simply disappear.
- Blewett depth measurements.
- Undetected spurious hollows remain in large fields (~6%) and among small detections, which were not screened.

## Consolidated changes (2026-11-05)

`rerun/changes_all_v3.csv` stacks the three group tables (g1 graph and joint model, g2 insolation, g3 catalog and size-frequency).
It has 397 rows, one per quantity quoted in v6.3 or supplement v13 or in the workshop drafts.
71 rows change a conclusion: 41 in both catalogs, 19 in v3b only, 3 in v3 only, and the rest are qualified in the `conclusion_change` column.
Rows are sorted with conclusion changes first.

## Catalog decision and follow-up reruns (2026-11-06)

Michael chose v3 as the primary catalog after reviewing `qc/decide_U1.png`, `decide_U2.png` and `decide_REF.png` (made by
`qc/make_u_sheets.py`; tile ids in `qc/decide_ids.csv`). v3b is kept as a sensitivity check.
Three more groups ran in `rerun/`:
- `g4_completeness/`: completeness grids rebuilt from the v3 and v3b refit detection models (`completeness_10ppd_v3.npz`,
  `_v3b.npz`), full g1 and g2 suites rerun on them (catalog tags v3r and v3br), the S2 reference-size checks, and the primary
  size-distribution numbers under the refit model.
- `g5_figures/`: every catalog-dependent figure redrawn on v3 (`*_v3cat.png`), including a new panel e in Figure S5.
- `g6_remaining/`: Stage 2/3 against incidence, the S9/S16 color mask, knob colors, the DEM aspect check, nested_test2, and
  a scan for uncovered numbers.
`rerun/changes_all_v3r.csv` stacks all six groups' change tables. g4 rows carry `v3r` (new grid) beside `v3` (old grid).
The step-versus-linear decision figure is `rerun/g2_insolation/out/compare/fig_step_vs_linear.png` (`fig_step_vs_linear.py`).

## g13_traced (6 Nov 2026)
Terrain test, positions and slope, knob comparison, Figures S2 and 2d rerun with Michael's traced outlines for 62 of the 64 main hosts
(curvature v3 for the other 259). See rerun/g13_traced/README.md, results.json and changes_traced.csv. Main changes from g11: individual
hollows no longer favor the top of peaks (fields still do); floor beside the peak (0.05 R, radius-stratified) now includes 1; the nested
class order becomes peak > floor > wall. Peaks/walls 7.19 (4.00-12.80), floors/walls 2.05 (1.29-3.37), equator/pole walls 1.57 (0.83-2.74).
