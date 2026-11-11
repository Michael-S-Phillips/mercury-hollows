# g1_graph_joint: containment graph, vents, degradation state and joint model on the v3/v3b catalogs (6 Nov 2026)

Reruns of the analyses behind main-text Sections 3.3 (nesting), 3.4 (vents), 3.6 (hosting by degradation state; depth
results use the Blewett depths and were not rerun), 3.8 (joint model) and supplement S11, S14, S15, S16, with the
spurious detections removed (`../../hollows_canonical_v3.csv`, `_v3b.csv`; see `../../README.md`).
**Every number checked, with the v2 record, v3, v3b and whether a conclusion changes: `changes_g1.csv` (218 rows).**

## Inputs
- Host assignments: `code/rerun/out_v3/host_assignments_v3.csv` / `_v3b.csv` (`hostfiles.py` selects them with the
  environment variable `CATALOG` = v3 or v3b; `v2` points to the original files and was used only to check the copies).
- The crater table `insolation_thermal/data/test2_craters_h18.csv` carries the v2 hosting flag `hosts` (any assigned
  hollow; it equals the host set of host_assignments.csv exactly). `inputs/make_inputs.py` writes
  `inputs/test2_craters_h18_v3.csv` / `_v3b.csv` with `hosts` recomputed (534 and 441 hosts) and every other column unchanged.
- Hosts within 0.9 radii (h09) in the tests: 606 (v2), 446 (v3), 370 (v3b); any-hollow hosts in the tests: 705, 532, 440.

## Scripts (copies; only input paths, output paths and environment settings changed)
| Folder / file | Copied from | Outputs |
|---|---|---|
| `crater_dag.py`, `crater_dag_tests.py`, `crater_dag_sens.py`, `run_dag.sh` | analysis_v9 (run_v9.sh variants) | `<CAT>/` records and nodes; `<CAT>/uplift_lunar`, `kring_0p8`, `kring_1p3`, `caloris_orgel` |
| `compare_v8_v9.py` | analysis_v9 (run list changed to v2 vs this catalog) | `<CAT>/compare_v2_<CAT>_record.json` |
| `summary_model/` build_features, add_xrs, common, step3, step4, step7, step7_sizestep_check, step8_figure_improved_v2, step9_color | analysis_workshop/summary_model | `summary_model/<CAT>/crater_features*.csv`, `summary_model/<CAT>/out/` (improved.json, dominance.json, joint_model.csv, color_record.json, trees.json, sizestep_check.json, xrs_rate_ratios.csv, figS9_joint_model_<CAT>.png) |
| `degradation/` host_age_tests, match_kinczyk, kinczyk_tests, stage_by_class | analysis_v3, analysis_workshop/kinczyk_degradation | `degradation/<CAT>/` |
| `vents/` test8_vents, vent_assignment_check | analysis_workshop, vent_check | `vents/<CAT>/` |
| `depth_nested/depth_nested.py` | analysis_v8/depth_nested | `depth_nested/<CAT>/` (unchanged results) |

Two copies needed more than a path change: in `host_age_tests.py` the depth block reads a cloud-only file and now runs only
with `RUN_DEPTH=1` (depths do not depend on the catalog); `compare_v8_v9.py` compares v2 with v3/v3b instead of v8 with v9.

New scripts: `parents_and_examples.py` (S14 parents without hollows and the Figure 3 caption counts; no script was kept,
its definition reproduces 149/75/0.08/0.21/22/4 exactly), `degradation/kinczyk_match_sensitivity.py` (reproduces
match_sensitivity.json), `summary_model/cached_pool.py` + `run_step.sh` (runner, below), `summary_model/count_model_check.py`,
`summary_model/boot_failure_check.py`, `make_changes.py` (writes changes_g1.csv).

**Check of the copies.** With `CATALOG=v2` the copies reproduce the v2 records exactly (crater_dag, tests, sens, kinczyk_tests,
stage_by_class, match files, test8_vents, vent_assignment_check, depth_nested; host_age_record to 1e-15),
crater_features.csv and crater_features_xrs.csv are byte-identical, the tree SHAP shares are identical, and the count-model
check gives alpha 9.907 and the recorded shares.

## Running on this machine
Shell calls end after ~1-2.5 min and kill background jobs. `summary_model/cached_pool.py` runs a copied script unchanged
but replaces `multiprocessing.Pool` with a pool that caches each bootstrap job (`summary_model/<CAT>/cache/`) and stops
when its time budget is spent; calling it again resumes. Each job has its own seed, so the result does not depend on the
split. Bootstrap sizes are the originals (step3 B = 200, step7 B = 100, step9 B = 60).

## Numerical problems found (they affect intervals, not point estimates)
1. **Count model.** statsmodels `NegativeBinomial.fit` in step7 does not converge for v3b (llf nan; alpha stays at its start
   value 1.05), and for v3 only with single-threaded BLAS. v3 was rerun with default threads (alpha 12.4). For v3b, and for
   all count rows of changes_g1.csv, `count_model_check.py` estimates alpha from Poisson start values (v2 9.91, v3 12.41,
   v3b 16.33) and reruns step7's own dominance code. The count bars in `figS9_joint_model_v3b.png` come from the failed fit;
   use `summary_model/v3b/out/count_model_check.json` instead.
2. **Thinned-occurrence bootstrap.** The cloglog IRLS of step7/step9 fails in some bootstrap draws (a group then lowers R^2).
   The proportion of failed draws is 2 of 60 (v2), 12 of 100 (v3) and 31 of 100 (v3b). These draws widen the v3/v3b share
   intervals and push some below zero, e.g. size 13 (8-38)% in v3. `summary_model/<CAT>/out/boot_failure_check.json` gives
   intervals without the failed draws; for v3 these are size 9-21, insolation 11-22, degradation 31-48 and nesting 12-20%.
   The logistic (step3) bootstrap has no failures.
3. The S16 check on craters >= 10 km alone fails to converge in v3b (negative shares). The 46% for degradation state is not
   usable there.

## Headline changes (v2 -> v3 / v3b)
- Nesting OR 1.85 (1.54-2.22) -> 1.68 (1.35-2.09) / 1.78 (1.40-2.26). The OR for craters inside hollow-bearing parents
  rises from 2.23 to 3.02 / 3.08. With degradation state controlled it is 2.31 (p = 0.006) / 1.96 (p = 0.043), against
  1.68 (p = 0.055) in v2.
- Became non-significant: fill-only craters (2.09 -> 1.55 (0.82-2.94) / 1.71 (0.87-3.37)) and craters inside basins
  (1.39 -> 1.22 / 1.33). In v3b only, craters inside parents >= 300 km (1.40, 0.93-2.10) and the peak position (2.21, 0.98-4.98).
- AIC: nesting lowers the AIC by 38 rather than 61. Every depth-of-origin model is still >= 25 AIC above nesting alone (v2: 56).
  The position model now has a lower AIC than nested alone (by 2.2 / 0.8), and adding ln z_max to ln z_own gains 22 / 17
  rather than 11.
- Vents OR 4.88 -> 3.48 (1.97-6.17) / 3.22 (1.71-6.07). Isolated craters with vents: 2.53 (1.20-5.32) / 2.22 (0.96-5.13; v3b n.s.).
  In the logistic joint model the block-bootstrap vent interval now includes 1.
- Herrick states: fresh 4.87 / 5.30 and standard 1.94 / 2.04 (v2 4.86, 1.83). Kinczyk per class: 2.89 / 2.84 (v2 2.77).
  In v3b, Class 5 hosts less often than Class 4 (22.7% vs 24.7%), and Kinczyk no longer beats Herrick by AIC with
  diameter and latitude alone (-1.9).
- Joint model: g 0.88 -> 0.955 / 1.033 and R^2 0.25 -> 0.22 / 0.23. Shares for degradation state, insolation, nesting and
  size: 39, 18, 16 and 13% in v3 (37, 18, 17 and 13% in v3b), against 31, 25, 16 and 17% in v2. Nesting and size therefore
  swap ranks, and the well-imaged order (size before nesting) no longer matches. In the logistic and tree versions,
  insolation is no longer among the top two.

## Not rerun
- `peak_wall_test/nested_test2.py`: none of its numbers appear in v6.3 (the 300 km matched comparison in the paper comes
  from crater_dag_sens.py and was rerun). Its per-crater hollow counts (morph_all.csv, n_all) come from the DEM terrain
  pipeline measure_all.py, which uses another host-assignment file.
- Depth results (Blewett et al. 2016 depths): unchanged by construction. depth_nested was rerun as a check and is identical.
- Figure 3 (fig_nested_examples.png) must be redrawn: in both v3 and v3b, 3 of the 5 craters inside Lermontov host hollows
  (v2: 4 of 5). Figure S9 was regenerated in `summary_model/<CAT>/out/`. The v3b count bars in that figure come from the
  failed fit (see above).
