# Crater containment graph and hollow hosting (analysis_v6, 29 Oct 2026)

Question (Michael): does a crater that forms inside an earlier crater, and so excavates or uplifts material the earlier
crater already brought up, reach the hollow source more often? Method adapted from mars-crater-graph (containment edges only).

Scripts (run in order, argument = hollows-emissivity root):
- crater_dag.py: containment forest of 31,571 Herrick et al. (2018) craters plus 53 Orgel et al. (2020) certain/probable
  basins (duplicates of cataloged craters removed); immediate parent = smallest feature >= 2x larger whose rim contains the
  child's center; position in parent (peak rho < 0.2, floor 0.2-0.6, wall >= 0.6); first-order depth of origin
  z_max = z_base + max(0.1 D_tc, 0.086 D^1.03); fill from the depth deficit. Writes crater_dag_nodes.csv, crater_dag_record.json.
- crater_dag_tests.py: stricter hosting (hollows inside 0.9 radii), regional hosting rate within 500 km as a covariate,
  parent type, inheritance (parent hosts), conditional logit within parents. crater_dag_tests_record.json.
- crater_dag_sens.py: subsets (completeness >= 0.5, 10-50 km, preservation class, |lat| <= 60, parents < 300 km) and a
  matched comparison with isolated craters within 300 km. crater_dag_sens_record.json.

Controls in every logistic fit: log10 D, |lat|, completeness for 1 km hollows at the center, northern-plains flag, dT.

Results (hosting = >= 1 hollow inside 0.9 radii unless stated):
- 14,916 of 31,571 craters (47%) lie inside a feature >= 2x larger; chains reach 3 ancestors.
- Nested craters host hollows 1.85 (1.54-2.22) times as often (odds), with the regional rate controlled. Same with
  >= 3 hollows (1.59, p = 0.008), completeness >= 0.5 (1.75), 10-50 km (1.93), preservation class controlled (1.52),
  |lat| <= 60 (1.84). Matched isolated craters within 300 km of similar size and completeness: rate ratio 2.1.
- Parent type: inside craters < 300 km, 2.88x (3.1x when only such parents are compared); inside craters >= 300 km, 1.40x
  (p = 0.06); inside basins, 1.47x.
- Inheritance: children of hosting parent craters host 2.39 (1.73-3.29) times as often; 2.1 with preservation class controlled.
- Depth of origin NOT supported: within the same parent, peak/floor children vs wall children 1.23 (p = 0.20); children
  whose excavation should reach only fill host as often as other nested children; ln z_max predicts hosting worse than the
  child's own ln z_own (AIC 5,573 vs 5,111).
Reading: nesting matters, but through the material the earlier impact processed (exposed, brecciated, and, in the model,
charged with sulfur), not through an additive depth of origin. Every child (excavation >= 0.42 km) cuts through a charged
layer 15-46 m thick wherever it sits in the parent.
Caveats: MESSENGER targeted known hollow sites, so children of hosting parents may be better imaged than completeness at
the center captures; inheritance can also reflect local source richness shared by parent and child.

---
# v8 (30 Oct 2026): one Caloris node; lunar-uplift sensitivity

The text above is the v6 README. v8 changes crater_dag.py only in how Caloris enters the graph (plus options):
- v6 kept both the Herrick et al. (2018) crater 200153 (D 1595 km, 31.93 N 164.10 E) and the Orgel et al. (2020) basin
  Caloris (D 1556 km, 29.91 N 161.43 E). Their centers lie 130 km apart, beyond the 0.15-radius (117 km) duplicate tolerance.
- v8 keeps one node with the Orgel center and diameter. The kept record is the catalog crater 200153 (kind crater, its
  hollow-host flags, which host_assignments.csv keys to the catalog id); the Orgel row is dropped, as the duplicate step does
  for every other cataloged basin. Children are all features inside either footprint (Orgel or catalog circle), so the
  children of both v6 nodes (437 + 50 = 487) attach to it; rho is measured from the Orgel center, so the 50 children beyond
  the Orgel rim have rho >= 1 (wall, as in v6). The Caloris fill stays 3.0 km.
- Options (environment): CRATER_DAG_CALORIS=orgel keeps the Orgel record instead (kind basin, no host flags);
  CRATER_DAG_SU=lunar uses SU = 0.022 D^1.45 (Cintala and Grieve 1998, eq. 12; checked by refitting their Table 4:
  0.0219 D^1.454, r = 0.982), double uses 2 x 0.086 D^1.03; CRATER_DAG_OUT sets the output folder (all three scripts).
- crater_dag_tests.py records coefficients and SEs unrounded (v6 rounded them to 3 decimals); fits are unchanged.

Runs (argument ~/mnt/hollows-emissivity):
- analysis_v8/              default: merged Caloris (catalog record), terrestrial uplift; all three scripts.
- alt_caloris_orgel/        CRATER_DAG_CALORIS=orgel; all three scripts.
- uplift_lunar/             CRATER_DAG_SU=lunar, merged Caloris (catalog record); crater_dag.py only.
- v6_rerun_fullprecision/   v6 scripts rerun unchanged except unrounded test coefficients; identical to the analysis_v6
                            records to their rounding (the "before" values).
- compare_v6_v8.py -> compare_v6_v8_record.json: before/after table (ORs exp(b) with Wald 95% intervals).
- detection_function_v2_check.py -> detection_function_v2_check_record.json: refit of the v2 detection model (92,382 pairs)
  and check of the 1 km completeness grid (used by figscripts/make_figS_completeness_v8.py).

Effect of the merge: counts and all hosts-outcome fits change by <= 0.4 AIC. With the catalog record, the 437 children that
v6 placed in the Orgel basin become children of a >= 300 km crater that hosts hollows, so they enter the inheritance tests:
inheritance 2.39 (1.73-3.29) -> 2.23 (1.62-3.07); preservation controlled 2.10 (1.25-3.55) -> 1.68 (0.99-2.86), p = 0.055;
C1 >= 0.5 5.98 -> 4.64; parent kind >= 300 km crater 1.40 (0.99-1.99) -> 1.56 (1.16-2.10), basin 1.47 -> 1.39. Keeping the
Orgel record instead leaves these near v6 (2.36, 2.05, 6.20; 1.31, 1.50).
Lunar uplift: ln z_own improves (AIC 5,111 -> 5,088) and z_base weakens (p 0.02 -> 0.51), but every depth-of-origin model
still fits worse than the nesting indicator (dAIC vs nested +37 to +469).
