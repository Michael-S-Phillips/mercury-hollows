# g4_completeness: completeness grids from the refit detection model, and the reruns that use them (6 Nov 2026)

The v3 catalog is now the primary catalog, and the detection model refit on the v3 hollow-image pairs (g3_catalog, "Detection
model refit") replaces the v2 model, because the v2 model was fit on pairs that included the spurious detections. This folder
rebuilds the completeness grids and every completeness-derived input from the refit model and reruns the analyses that use them.
Catalog tags: `v3r` = v3 catalog with the v3 refit model (primary); `v3br` = v3b catalog with the v3b refit model (sensitivity).
`v3` and `v3b` below always mean the earlier reruns, which used the old (v2-model) grid. No file outside this folder was changed.

Every quantity, with v2, old-grid v3/v3b and new-grid v3r/v3br, is in `changes_g4.csv` (368 rows). The column
`conclusion_change` compares old-grid v3 with new-grid v3r (and v3b with v3br). The column `conclusion_change_v2_to_v3r` repeats the
earlier tables' comparison with v2, using v3r/v3br in place of v3/v3b.

## How the grids were built
- `build_completeness.py` is a copy of the coverage, completeness-grid and attach blocks of `code/rerun/build_hollow_products.py`.
  Only the detection-model parameters (read from `g3_catalog/sfd/out_<v3|v3b>_refit/theta_out_v2.json` instead of the in-script
  fit), the input catalog and the paths were changed. The grid formula is unchanged: C = 1 - prod over scale bands of
  (1 - p)^n, with p at 60 deg incidence and the band's representative pixel scale, at 11 reference diameters on the 10 ppd grid.
- Outputs: `completeness_10ppd_v3.npz` (refit p_max 0.272, x50 35.3 px, sigma_lo 0.585, x_hi 187.2 px, sigma_hi 0.346,
  b -0.104) and `completeness_10ppd_v3b.npz` (0.278, 34.3, 0.578, 186.9, 0.351, -0.107).
- The attach step writes the per-hollow and per-field columns (cov_*, C_D0.5km, C_D1km, C_D2km, C_D5km, C_own).
- `g2/prep_inputs.py` builds `g2/inputs/v3r` and `v3br` from the g2 inputs of v3 and v3b. Only the completeness columns change:
  the hollow and field columns come from the attach step, and C1 of `crater_features.csv` is resampled at the crater centers
  exactly as `crater_dag.py` samples it. The file `crater_dag_nodes_C1.csv` gives C1 at the analysis_v9 graph nodes for
  `dT_hosting_completeness.py`. In g1, C1 enters through the copied `crater_dag.py`, which reads the new grid when CATALOG is
  v3r or v3br; the joint-model tables are then rebuilt by the copied `build_features.py` and `add_xrs.py`.

**Checks of the chain.** The rebuilt coverage grid is identical to `out_v2/coverage_10ppd.npz`. The v2 configuration
reproduces `out_v2/completeness_10ppd.npz` to 2.4e-6 at every reference size; the difference is not exactly zero because the
in-script fit of the v2 model was never saved, and `theta_out_v2.json` is a refit within the Nelder-Mead tolerance. The attach
code applied to the stored v2 grid returns the C columns of `hollows_canonical_v3.csv` and the g2 fields exactly. The C1 sampling
applied to the stored grid returns the stored C1 of `crater_features.csv` and of `analysis_v9/crater_dag_nodes.csv` exactly. The
C1 that `crater_completeness.py` computes independently from the refit model matches the new grid to 1e-7. In g1, only the C1
column differs between the v3 and v3r feature tables.

## How much the grid changes (`compare_grids.py`, `grid_comparison.json`, `grid_comparison.png`)
New/old over covered cells, 5th / 50th / 95th percentiles:

| Reference | v3 refit grid | v3b refit grid |
|---|---|---|
| 0.5 km | 1.034 / 1.066 / 1.122 | 1.055 / 1.116 / 1.200 |
| 1 km | 1.043 / 1.094 / 1.124 (area-weighted 1.079 / 1.098 / 1.122) | 1.062 / 1.159 / 1.196 |
| 2 km | 1.000 / 1.105 / 1.125 | 1.013 / 1.167 / 1.206 |
| 5 km | 1.000 / 1.057 / 2.739 | 1.005 / 1.079 / 2.948 |

The rise is not uniform. Completeness rises about 9-11% in cells that are poorly imaged or imaged at 10-100 m per pixel, but only
about 1% where the finest images are below 10 m per pixel or where the old C exceeds 0.95. In those cells a 1 km hollow spans
about 160 pixels, which is near the refit model's lower upper turnover (187 px instead of 196 px). At 5 km the wider decline of
the refit model raises C by up to about 2.7 times (95th percentile) in cells covered only by the finest images. At crater centers, C1 rises by a factor of
1.086 / 1.096 / 1.121 (5th / 50th / 95th percentiles over all craters), and by 1.030 / 1.101 / 1.123 over the h09 hosts.

**Rerun depth.** At 1 km, the 5th and 95th percentiles lie -4.6% and +2.7% from the median, which is outside the +-3% criterion,
so the full suites were rerun instead of spot checks. (The area-weighted range, -1.7% / +2.2%, would have passed, but the
departure is systematic: it is largest in the best-imaged cells, which hold many hosts.)

## What was rerun, and the check of each copy
| Folder | Copies (source) | Runs | Check of the copy |
|---|---|---|---|
| `g2/` | all 22 g2 scripts and `run_g2.sh` (`make_copies_g4.py`: folder depth; grid path and C1 source for v3r/v3br); `compile_changes.py` | full suite for v2, v3, v3b, v3r, v3br | v2, v3, v3b records identical to `g2_insolation/out` (23 files each, `out/repro_*.json`); pass A of `compile_changes.py` rebuilds `changes_g2.csv` byte for byte |
| `g1/` | DAG, tests, sensitivity and variant runs, parents, degradation, vents, summary model (`make_copies_g1.py`: CATALOG v3r/v3br and the grid path); `make_changes_g4.py` | v2, v3, v3r, v3br; summary model (bootstraps) for v3r only | v2 and v3 records, nodes, parents, degradation and vents files identical to `g1_graph_joint`; v2/v3 feature tables identical; pass A rebuilds `changes_g1.csv` byte for byte |
| `s2/` | `ref_size_insolation.py`, `ref_size_check.py`, `crater_completeness.py`, `hosting_fits.py`, `joint_thinned.py` (`make_copies_s2.py`) | v2, v3, v3b, v3r, v3br | v2 outputs identical to the stored workshop records |
| `sfd/` | `size_law.py` (g3 copy) | v2, v3, v3b, v3r, v3br | v2, v3, v3b records identical to g3 |
| this folder | `figS_completeness_v3cat.py` (from `figscripts/make_figS_completeness_v8.py`) | v3 refit | the v2 configuration reproduces `figS_completeness_v8.png` pixel for pixel |

The size-distribution numbers under the refit model were already computed by g3 (`g3_catalog/sfd/out_v3_refit`, `out_v3b_refit`,
`sfd_derived_refit.json`); they are collected in `changes_g4.csv`, and every size-distribution quantity of `changes_g3.csv` has a
refit value. The only size-distribution script that reads the 10 ppd grid, `size_law.py` (closure, not in v6.3), was rerun here.

## Results
**Effect on v3.** The new grid moves numbers slightly and changes no conclusion for v3r. In g2, 34 of 94 rows change value and
none changes a conclusion. In g1, 128 of 218 rows change value (mostly in the second digit) and none changes a conclusion. Examples (old-grid v3
-> v3r):
- Field hot-to-warm ratio at 0-15N: 1.71 (1.15-2.55; p = 0.0095) -> 1.75 (1.22-2.60; p = 0.0045). The 15-30N crater odds ratio:
  2.64 (1.55-4.48) -> 2.60 (1.53-4.40). The 30-45N crater odds ratio: 1.49 (0.93-2.40) -> 1.49 (0.93-2.39; still includes 1).
- Delta T_w x log D interaction p with completeness: 0.18 -> 0.21 (AIC 3,326 -> 3,315); best single step: AIC -1.6 -> -2.1
  against the linear term.
- Nesting odds ratio 1.68 (1.35-2.09), unchanged. Parent-host inheritance 3.02 (2.00-4.55) -> 2.98 (1.98-4.50); with degradation
  state controlled 2.31 (p = 0.006) -> 2.28 (p = 0.007). Vents 3.48 (1.97-6.17) -> 3.46 (1.95-6.13).
- Subsets defined by completeness grow: craters inside parents with C >= 0.5, 365 -> 414 (inheritance 4.53, p = 0.0064 ->
  0.0023); nesting with C >= 0.5, 1,119 -> 1,252 craters; well-imaged craters (C >= 0.2) in the joint model, 3,835 -> 4,158.
- Joint model: g 0.955 -> 0.947; shares for degradation state, insolation, nesting and size 39.1, 18.4, 16.1, 13.0 -> 39.4,
  18.1, 16.2, 13.1%. All ranks unchanged. Every containment-test AIC falls by about 15, and the AIC differences between models change by 0.6 or
  less.

**Effect on v3b (sensitivity).** Two results change between v3b and v3br. Parent-host inheritance with degradation state controlled goes
from 1.96 (1.02-3.76; p = 0.043) to 1.90 (0.99-3.65; p = 0.054), and the three rows that report this test are flagged. The
interval of the 60-75N field ratio now excludes 1 (0.67, 0.43-0.99; p = 0.056), as it already did in v3.

**S2, reference size.** The grid changes nothing here; the changes from the paper come from the catalog. With completeness for 0.3, 0.5 or 1 km hollows, the hot
vs warm odds ratios are 2.19-2.63 within 15 deg, 2.60-3.13 at 15-30N and 1.49-1.69 at 30-45N in v3r (paper: 2.2-2.6, 3.2-3.7 and
1.76-1.98; old-grid v3: 2.19-2.63, 2.62-3.11, 1.49-1.69). The odds per 10 K of Delta T_w are 1.36-1.41 (paper 1.43-1.47). The 30-45N excess is
significant only with 0.3 or 0.5 km completeness (p = 0.035, 0.028) and not with 1 km (p = 0.10), 2 km or 5 km, so the S2
sentence that it "loses significance only with completeness for 2 or 5 km hollows" no longer holds; this was already true for
old-grid v3. With 2 km completeness, the excesses within 15 deg (p = 0.057) and at 15-30N (p = 0.063) also lose significance. The
fit comparison also changes: in the containment fit, 1 km now fits better than 2 km (log-likelihood -1,493.8 against -1,508.1;
paper: -2,100 against -2,093). In the joint model, 0.7 km fits slightly better than 1 km (-1,357.0 against -1,362.2; paper:
1 km best).

**S2, size-weighted completeness E.** The exponent for v3r is 1.6-1.7. Nesting 1.75 (1.42-2.17), vents 4.3 (2.5-7.3), fresh vs
degraded 3.3 (2.1-5.2); hot vs warm odds ratios 2.07 (1.41-3.06), 2.30 (1.38-3.81) and 1.47 (0.91-2.36; p = 0.12); 1.38 per 10 K
(exponent 1.7; exponent 1.6 differs by at most 0.02). The 30-45N odds ratio with E no longer excludes 1 (paper: 1.70, 1.08-2.68),
as in old-grid v3. Exponents 1.3 and 1.9 change these odds ratios by at most 4.1% from the 1.6 fit (5.1% from the 1.7 fit). In
the joint model, E moves group shares by up to 1.6 percentage points (paper: 1.2) and still fits worse than C1.

**Size distribution with the refit model (now primary).** Corrected N(>5 km) 105 (84-128) against 96 cataloged; N(>0.3 km) 21,281
(20,426-22,078); N(>1 km) 4,674. Exponents: 1.55 (0.3-0.7 km), 1.43 (0.7-3 km), 1.52 (1-3 km), 3.26 (2.91-3.62) above 3 km, 1.96
over 2-6 km. The range 0.3-3 km is 1.4-1.6 again (old-model v3: 1.5-1.6). Median detection probability at cataloged sites: 54.2% at
1 km, 4.3% at 0.3 km, and 84.4% and 86.3% at 3 and 5 km. The power law predicts 360-490 hollows larger than 5 km, so the deficit of
large hollows holds. The catalog holds one in 4.4 to one in 5.6 hollows larger than 0.3 km (paper: one in five to one in six).
Section 4.4: 4.7-5.8 times as many hollows larger than 0.1 km as larger than 0.3 km, and 117,000-186,000 on the searched surface.
The best-imaged sites (C(0.3 km) >= 0.3) number 275, with 149 hollows larger than 0.3 km, and the exponent over 0.3-1 km is
2.34 (1.80-2.90). Exponents above 0.2, 0.3 and 0.5 km are 1.57, 1.67 and 1.68.

**Figure.** `figS_completeness_v3cat.png` is the supplement figure redrawn with the v3 refit grid and model.

## Caveats
- **Count model.** The negative binomial fit of `step7_improved.py` does not converge for v3r (alpha 1.05 with single-threaded
  BLAS, infinite with default threads). The count-model rows use `count_model_check.py` (alpha 12.52), as g1 did. The count bars
  of `g1/summary_model/v3r/out/figS9_joint_model_v3r.png` come from the failed fit.
- **S16 check for craters >= 10 km alone.** This fit fails to converge in v3r (negative cover share), so its 41.5% is not usable.
- **Bootstrap failures.** In the thinned-occurrence bootstrap, 13 of 100 draws fail in v3r (v3: 12). Intervals without the failed
  draws are in `boot_failure_check.json`.
- **v3br coverage.** The joint model (summary model) was not rerun for v3br; its rows say "not run".
- **Figure scale.** In the figure, the panel b color scale stops at 0.30, but the refit model reaches 0.325 below about 20 deg
  incidence, so those values are clipped.
- **Not rerun.** Other completeness-dependent analyses were not rerun with the new grid and should be when convenient: the g3
  terrain table (area weighted by C), the Caloris tests, catalog_counts, the melt-source positions and knob tests, and the
  workshop scripts outside g1-g3.
- **Text updates for S2.** S2 should quote the refit model: 85,602 pairs, 10,214 of them detections, from 6,299 hollows.
- **Scratch files.** `_scratch/` holds the rebuilt coverage grid, the v2-check grid (25 MB) and the attach outputs. They could not
  be deleted here.

## Files
`build_completeness.py`, `compare_grids.py`, `make_changes_g4.py`, `figS_completeness_v3cat.py`; `completeness_10ppd_v3.npz`,
`completeness_10ppd_v3b.npz`, `grid_comparison.{json,png}`, `figS_completeness_v3cat.{png,pdf}`, `changes_g4.csv`;
`g2/` (inputs/v3r, inputs/v3br, out/<cat>/, changes_g2_pass{A,B,C,D}.csv); `g1/` (<cat>/, degradation/, vents/,
summary_model/v3r/out/, changes_g1_pass{A,B,C,D}.csv); `s2/out/`; `sfd/size_law_record_*.json`.

- **Caloris tests rerun (9 Nov 2026).** g3_catalog/caloris/caloris_v9.py and caloris_v9_controls.py with G3_GRID=v3r (C1 from
  g1/v3r/crater_dag_nodes.csv) write g3_catalog/caloris/out_v3r; g5_figures/caloris_test_record_v9_v3cat.py with G5_GRID=v3r writes
  _scratch/caloris_v3r, from which Figure 2 (fig1_exposure_v9_v3cat_traced_r3_v3r) was redrawn. Used in v6.14 final and supplement v25 final.
