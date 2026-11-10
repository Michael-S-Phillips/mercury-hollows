# g5_figures: paper figures redrawn on the v3 catalog (6 Nov 2026)

Every figure of `Hollow_Formation_Model_Paper_v6_3.tex` and `supplementary-materials_v13.tex` that plots catalog hollows or host counts, redrawn on the screened v3 catalog (`code/rerun/out_v3`). The size distribution uses the refit detection model (`g3_catalog/sfd/out_v3_refit`). Each copied script starts with a comment that names its source and lists what changed. No existing file was modified, and nothing was written into `manuscripts/`.

## Figures

| Figure (file in v6.3 / v13) | New file | What changed |
|---|---|---|
| Fig. 1 `fig_sfd_v1.png` | `fig_sfd_v1_v3cat.png/.pdf` (`sfd_fig_v1_v3cat.py`) | The figure uses the v3 refit outputs. The exponent lines are now read from `sfd_numbers.json` instead of being typed in. The band is still 1.4–1.6 (fitted range 1.43–1.55). The line above 3 km is now 3.26 (2.91–3.62), labeled 3.3 (v2: 2.7). The y-axis top of panel b was raised from 3.5 to 4.0, because the 3.62 bound was clipped. |
| Fig. 2 `fig1_exposure_v9.png` | `fig1_exposure_v9_v3cat.png/.pdf` (`make_figs_pr_v9_v3cat.py`; `caloris_test_record_v9_v3cat.py`) | Panel a shows the 935 v3 fields (g2 rule), sized by their v3 hollow counts. Panel b uses the g2 v3 densities: inside/outside is 0.35 at 40–60°N, and the 0–30°N outside density is 52 (v2: 61). Panel c changes only in two isolated bins: 9/78 and 38/137 (v2: 10/78 and 41/137); the Caloris counts are unchanged. Panel d uses g3 `terrain_table_record.json`, because `peak_wall_by_setting.json` has no script; the point estimates match it exactly for v2. Hosts: 321 (190/68/63), where the caption says 352. |
| Fig. 3 `fig_nested_examples.png` | `fig_nested_examples_v3cat.png/.pdf` (`nested_examples_v3cat.py`, new) | No script had been saved, so I wrote a new one with the same layout, mosaic and symbols. Lermontov now has 3 of 5 hosts. The 5.7 km child at 14.5°N, 311.9°E keeps its place but is drawn as a non-host (blue ring), because its only hollow was screened as a crater bowl. Four screened hollows leave the Lermontov view. Praxiteles is unchanged (4 of 6). The base rate is 1.4% (446 of 31,554), down from 1.9%. The arrows now end at the child rims. |
| Fig. 4 `fig3_insolation_v9.png` | `fig3_insolation_v9_v3cat.png/.pdf` (copy of g2 `out/v3`) | The figure renders correctly; see the note on panel b below. Caption values: 217 fields; fit 0.99 + 0.49 cos 2λ (0.33–0.65); 485 hosts. |
| Fig. S5 `figS5_insolation_v9.png` | `figS5_insolation_v9_v3cat.png/.pdf` (`make_figs_pr_insolation_v9_v3cat.py`) | Panels a–d are g2's v3 panels. Panel e is new (see below). The record of panel e is `figS5_panel_e_v3cat_record.json`. |
| Fig. S1 `figS1_positions_v3.png` | not redrawn | It plots detections, not the screened catalog. No plotted detection, whether in the frame or from other images, belongs to a screened-out hollow. |
| Fig. S2 `figS2_terrain_examples.png` | `figS2_terrain_examples_v3cat.png/.pdf` | Redrawn with v3 host assignments. The image is byte-identical to the original: 83, 343, 136 and 80 hollows, all unchanged. |
| Fig. S3 `figS3_caloris_map_v9.jpg` | `figS3_caloris_map_v9_v3cat.jpg/.png/.pdf` | Redrawn on v3. The image is identical: the classes stay 144/21/14, and two craters change their counts without changing class. |
| Fig. S4 `figS4_caloris_examples.png` | `figS4_caloris_examples_v3cat.png` (`nested_figs_v3cat.py`) | Redrawn on v3. The image is identical: 0, 4, 71 and 0 hollows. |
| Fig. S7 `figS7_color_v1.png` | `figS7_color_v1_v3cat.png` | Panel c rings the 442 v3 hosts (v2: 563); panels a and b are unchanged. |
| Fig. S8 `figS8_joint_model_v2.png` | `figS8_joint_model_v2_v3cat.png` (copy of g1 `summary_model/v3/out/figS9_joint_model_v3.png`) | v3 shares, with the count model converged (alpha 12.4). |

Not catalog-dependent (not redrawn): main-text Fig. 5 `fig4_depth.png` (Blewett depths), plus `figS_tbase_phi_v1` and `figS_completeness_v8` (agent g4), which were skipped as instructed.

## Figure S5, new panel e
Caption fragment: "(e) Hosting probability for a crater of median size, latitude, and completeness (10.6 km, 29°, C = 0.026) against $\Delta T_\mathrm{w}$, for the 27,257 craters $\geq$ 5~km within 60$^\circ$ of the equator (414 hosts), from logistic fits with log diameter, absolute latitude, and completeness as controls: $\Delta T_\mathrm{w}$ bins with 95\% intervals, the steady (linear) rise, and the best single step, at 34~K, whose AIC is 1.6 lower than that of the linear fit."
The source panel also has a dotted 44 K line (the Stage 1 step), which I left out.

## Figure 4b under v3
Panel b still draws the fitted interaction line, which falls from 1.45 at 5 km to about 1.05 at 800 km. The panel gives no p value, so the line implies a size trend even though the trend is not significant (interaction p = 0.061). The band can contain a flat line near 1.25. The six bin points (1.35, 1.21, 1.60, 1.30, 1.13, 1.03) do not decline steadily. The same panel appears as S5b.

## Gaps and caveats
- Fig. S7 panels a and b still use the v2 hollow mask (`color_mdr`), which no group has rerun.
- Fig. S8 intervals include the 12 failed bootstrap draws of 100: size 8–38% and volcanic cover −5.5 to 10.8%. Without the failed draws, size is 9–21% (`boot_failure_check.json`).
- Fig. 2d intervals come from the bootstrap in `terrain_table.py`. For v2 they differ from the original file by up to 0.5.
- Stray writes, moved here because files cannot be deleted:
  - `_scratch/figure_numbers_pr_stray.json` was briefly written to `rerun/`.
  - `_scratch/stray_figS_positions_v3.cpython-310.pyc` came from `manuscripts/PlanetaryResearch_hollows/analysis_v3/__pycache__/`. An empty `__pycache__/` folder remains there.

## Messages that change under v3
- Fig. 1: the size distribution is steeper above 3 km (3.3, against 2.7), so the deficit of large hollows is stronger.
- Fig. 3: Lermontov now has 3 of 5 hosts, and the base rate falls to 1.4%. The figure still shows the point.
- Fig. 4b / S5b: the size dependence is no longer significant.
- S5e: a step fits slightly better than a steady rise.
- S8: insolation falls to 18%, and nesting and size swap ranks.

## 9 Nov 2026 additions
- `caloris_test_record_v9_v3cat.py`: env G5_GRID=v3r reads C1 from the v3r nodes table and the g3 out_v3r Caloris craters and writes
  _scratch/caloris_v3r/caloris_test_record_v9.json (Figure 2c of v6.14 final).
- `color_fig_v1_v3cat.py`: env COLOR_SRC=screened draws all panels from g6_remaining/color/out_v3 (screened-mask PCA, 31,412 craters;
  fit from depth_color_fits.json). `COLOR_SRC=screened python3 color_fig_v1_v3cat.py ~/mnt/hollows-emissivity/manuscripts figS7_color_v2_v3cat.png`
  wrote figS7_color_v2_v3cat.png (Figure S9 of supplement v25 final); quintiles and fit are identical to v1.
