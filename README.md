# Mercury hollows: code and derived data for a hollow formation model

This repository holds the analysis code, derived data, analysis records, and figure scripts behind the paper on how hollows form on Mercury (Phillips, in preparation; working title of manuscript v6: "Hollows on Mercury Form by Daily Cycling of Impact-Charged Sulfur"). It stays private until the paper is submitted; a restricted Zenodo record gives reviewers access during review.

The crater containment graph used in the paper (craters that formed inside older craters, and the interactive map) is a separate public repository with its own DOI: [mercury-crater-graph](https://github.com/Michael-S-Phillips/mercury-crater-graph).

## Layout

The folders keep the layout in which the analyses were run, so the scripts find each other's outputs by the same relative paths. Every script locates the repository root from its own location, or from `$MERCURY_HOLLOWS_ROOT` if that is set.

| Folder | Contents |
|---|---|
| `georeference/` | Recomputation of hollow-detection positions through the MDIS camera model (ISIS), and a check against the archived positions |
| `GIS/vector/` | `bickel2024_hollows_corrected.gpkg`: all 19,110 detections of the 0.7-confidence release at camera-model positions, with its README |
| `data/` | Bickel catalog position-error record; Herrick et al. (2018) crater catalog (`data/ME6M300TGT/herrick2018_all.csv`) |
| `code/rerun/` | Hollow products (`build_hollow_products.py`, `out_v2/`): merged hollows, hollow fields, detection QC, imaging coverage, detection model, completeness; host assignment, hot/warm-longitude, wall-orientation, incidence, and province analyses (`latitude/`) |
| `hollow-networks/` | Fitted per-frame detection function and declustering settings written by `build_hollow_products.py` |
| `manuscripts/GRL_insolation_hollows/` | Thermal model, crater hosting test, plains densities, depth analyses (`scripts/`, `depthctl/`), their inputs (`data/`) and records (`records/`) |
| `manuscripts/PlanetaryResearch_hollows/` | Analyses added for the Planetary Research paper (`analysis_v3/`), figure scripts (`figscripts/`), figures, and the data behind each figure (`figures/data/`) |
| `peak_wall_test/` | Peak, floor, and wall terrain densities; Caloris penetration test; Caloris example figure |
| `slope_test/` | Hollow occurrence against regional slope; depth against slope with shadow geometry corrected |
| `run_all.sh` | Reruns every analysis and figure that needs only the data in this repository |

## Where each result comes from

Section numbers refer to manuscript v6 and its supplement (S1–S16); the formation-model draft uses the same supplement sections. Paths are relative to the repository root; `G/` is `manuscripts/GRL_insolation_hollows/` and `P/` is `manuscripts/PlanetaryResearch_hollows/`.

| Result | Script(s) | Record or output |
|---|---|---|
| Recomputed positions; 1.90 → 0.41 km cross-image separation (S1) | `georeference/georeference_bickel.py` | `GIS/vector/bickel2024_hollows_corrected.gpkg` and README; `data/bickel-2024/bickel_position_error_record.json` |
| Relative offsets between overlapping images, Table S1; Figure S1 | `P/analysis_v3/pair_offsets.py`, `figS_positions.py` | `P/analysis_v3/pair_offsets_record.json` |
| 6,731 hollows, 1,228 fields, 92,382 hollow–image pairs, detection model, completeness (S2) | `code/rerun/build_hollow_products.py corrected code/rerun/out_v2` | `code/rerun/out_v2/`, `hollow-networks/results/per_frame_detection_function.json`, `code/rerun/dedup_before_after.json` |
| Host craters (707 among 31,571) | `code/rerun/latitude/aspect/assign_hosts.py` | `code/rerun/latitude/aspect/host_assignments.csv` |
| Thermal model and ΔT (S5) | `G/scripts/thermal.py` | `G/data/thermal_grid.npz` |
| Plains densities 0.41 and 0.13; Figure 1a–b (S6) | `G/scripts/nsp_mask2.py` (mask), `latcorr.py` → `latcorr3.py`, `latcorr5.py`; `P/analysis_v3/dT_plains.py` | `G/data/nsp_mask_10ppd_r40.npy`, `G/records/latcorr_density_record.json`, `latcorr_inc60_75_record.json`, `P/analysis_v3/dT_plains_record.json` |
| Caloris: 0 of 128, 14 of 32, 3 of 5 vs 1 of 57; Figure 1c (S7); Figures S3, S4 | `peak_wall_test/caloris_test.py`; `P/analysis_v3/figS_caloris_map.py`; `peak_wall_test/nested_figs.py` (Figure S4 is its panel B) | `peak_wall_test/caloris_test_record.json` |
| Peak 6.7 and floor 3.1 times wall density; Figure 1d; Table S2; Figure S2 (S4) | `peak_wall_test/peak_wall.py` (DEM, writes `rows_*.csv`), `analyze.py`; `P/analysis_v3/terrain_table.py`, `figS_terrain_examples.py` | `peak_wall_test/peak_wall_record.json`, `peak_wall_by_setting.json`, `P/analysis_v3/terrain_table_record.json` |
| Hot/warm longitudes 1.9 and 2.1; cos 2λ fit; Figure 3a (S8) | `code/rerun/latitude/hotwarm/hotwarm.py` | `code/rerun/latitude/hotwarm/hotwarm_record.json` (copy in `G/records/`) |
| ΔT × size: 1.60 and 1.18 per 10 K; Figure 3b (S3) | `G/scripts/test2_h18.py`, `P/analysis_v3/dT_hosting.py` | `G/records/test2_h18_record.json`, `P/analysis_v3/dT_record.json` |
| Equator-facing walls (S8) | `code/rerun/latitude/aspect/aspect_stats.py`, `aspect_robustness/aspect_robustness.py` | `aspect_record.json`, `aspect_robustness_record.json` (copies in `G/records/`) |
| Slope occurrence and slope figure (S9) | `slope_test/slope_occurrence.py` (DEM), `occ_analysis.py`; `P/figscripts/slope_fig.py` | `slope_test/occ_analysis_record.json` |
| Layer crossing 368 K: 22 and 10 cm (S15) | `P/analysis_v3/transition_layer.py`, `transition_layer_hilat.py` | `transition_layer_record.json`, `transition_layer_hilat_record.json`, `transition_layer_grid.npz` |
| Fresh craters 4.9 (3.4–7.0); depth ratio 1.02 (S13) | `P/analysis_v3/host_age_tests.py` | `P/analysis_v3/host_age_record.json` |
| Depth against temperature: 1.5-fold vs 794-fold; 0.77 poleward; Figure 4a (S10) | `G/scripts/test1.py`, `depth_reanalysis.py`; `slope_test/subsolar_for_depths.py`, `depth_points_slope.py`, `depth_slope_geometry.py`, `depth_slope_eqface.py` | `G/data/test1_rows.csv`, `test1_rows_inc.csv`, `G/records/depth_reanalysis_record.json`, `slope_test/depth_slope_record.json` |
| Depth against host size (Figure 4b); within-field correlation r = 0.69; r = 0.01 between craters (S11) | `G/depthctl/depth_structure.py`, `depth_structure_extras.py` | `G/depthctl/depth_structure_record.json`, `crater_median_depths.csv` |
| Depth against X-ray ratios, at most 11% (S12) | `G/depthctl/depth_vs_xrs.py`, `depth_vs_xrs_tiles.py` | `G/depthctl/depth_vs_xrs_record.json` |
| Sulfur release, vapor pressure, 15–46 m, iron budget (S16) | `P/analysis_v3/sulfur_charge.py`, `sulfide_s2.py` | `sulfur_charge_record.json`, `sulfide_s2_record.json` |
| Overburden and strength | `G/depthctl/strength_calc.py` | `G/depthctl/strength_calc_record.json` |
| Craters inside older craters (formation-model draft) | [mercury-crater-graph](https://github.com/Michael-S-Phillips/mercury-crater-graph) | that repository's `outputs/` |

| Figure | Script | Data table |
|---|---|---|
| Exposure (v6 Figure 2; draft Figure 1) | `P/figscripts/make_figs_pr.py` | `P/figures/data/fig1_exposure_data.csv` |
| Insolation, full (Figure S5) | `make_figs_pr.py` (saved as `fig2_insolation`, archived as `fig2_insolation_v3`) | `fig2_insolation_v3_data.csv` |
| Insolation, panels a–b (v6 and draft Figure 3) | `P/figscripts/fig3_insolation_v6.py` (crop of the above) | `fig2_insolation_v3_data.csv`, panels a–b |
| Depth (v6 and draft Figure 4) | `make_figs_pr.py` | `fig4_depth_data.csv` |
| Life cycle (v6 Figure 1; draft Figure 5) | `P/figscripts/lifecycle_fig_v5.py G/scripts P/figures/fig5_lifecycle_v5.png` | `fig5_lifecycle_v5_data.csv` (panel 7 is data; panels 1–6 are schematic) |
| Slope (supplement) | `P/figscripts/slope_fig.py` | `fig3_slope_data.csv` |
| Craters inside two hollow-bearing craters (draft Figure 2) | `scripts/make_figures.py` in mercury-crater-graph | that repository's `outputs/` |
| Figures S1–S4 | `P/analysis_v3/figS_*.py`, `peak_wall_test/nested_figs.py` | maps drawn from the MDIS mosaics and the USGS DEM; plotted hollows and craters are in `code/rerun/out_v2/` and `peak_wall_test/` |

`P/figure_numbers_pr.json` lists the numbers printed in the figures. The data tables come from `P/figscripts/figure_data.py`, which reads the plotted lines, points, bars, and interval segments of each panel just before the figure is saved.

## Running

Python 3.10 or later with the packages in `requirements.txt`.

```bash
bash run_all.sh          # 29 analyses and figures, about 3 minutes
git status; git diff --stat
```

`run_all.sh` lists the order and the arguments each script takes. Scripts that read files relative to their own folder (`peak_wall_test/`, `slope_test/`) are run from that folder. The scripts that need external data are run the same way once the data are in place (next section).

**Reproduction check (29 October 2026).** On a fresh copy, `run_all.sh` rewrites 21 analysis records. Nineteen match the archived records to within 10⁻¹¹, the CSV tables match to within floating-point printing, and the figures match to within font rendering (1–2.5% of pixels differ; `fig3_insolation_v6.png` is identical). Two records differ slightly:
- `code/rerun/latitude/aspect/aspect_record.json`: bootstrap intervals differ in the third significant figure under the current NumPy, and the archived record has three text fields (`interpretation`, `limitation`, `verdict`) that were added by hand.
- `G/depthctl/depth_structure_record.json`: `R2_crater_image` is 0.616 against 0.621 archived. That fit has collinear image and crater dummies, so its R² depends on the linear-algebra library; it is not cited.

## External data (not included)

| Data | Used by | Source |
|---|---|---|
| Bickel et al. hollow catalog, 0.7-confidence release (`Bickel_et_al_2024_MDIS_hollows_07_v1.0.csv`) | `build_hollow_products.py`, `georeference_bickel.py` | data release of Bickel et al. (2025); place in `data/bickel-2024/` |
| MDIS CDR cumulative index | `build_hollow_products.py` (downloads it), `depth_reanalysis.py` and `subsolar_for_depths.py` (`external/mdis_index/cdr_index.pkl`, `CDR_CUMINDEX.TAB`) | PDS Imaging Node, `msgrmds_2001/INDEX/` |
| MDIS EDR images | `georeference_bickel.py` (downloads them) | PDS Imaging Node, `msgrmds_1001/DATA/` |
| USGS MESSENGER global DEM, 665 m/pixel (`data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif`) | `peak_wall.py`, `measure_all.py`, `nested_build.py`, `slope_occurrence.py`, `aspect_robustness.py`, `figS_terrain_examples.py` | USGS Astrogeology (Becker et al., 2016) |
| MDIS global mosaics: 166 m/pixel monochrome (`data/`), 665 m/pixel enhanced color (`GIS/raster/`) | `figS_caloris_map.py`, `nested_figs.py` | USGS Astrogeology (Denevi et al., 2018) |
| Map-projected MDIS frames for Figure S1 (`data/mdis_frames/`) | `figS_positions.py` | made with ISIS from PDS EDRs |
| Northern smooth plains lineament shapefile (`data/data/NSP_map_data/`) | `nsp_mask.py`, `nsp_mask2.py` | source to be cited (open item in the manuscript); the resulting mask is included |
| ISIS 8 with MESSENGER kernels | `georeference_bickel.py` | USGS Astrogeology |

## Changes made for the archive

The scripts are the ones that produced the records, with these changes only:
1. Absolute paths from the machines where the analyses ran (`/Users/...`, `/home/claude/...`, `/mnt/user-data/...`, `/sessions/...`) are replaced by paths relative to the repository root; a two-line header in each script defines that root.
2. Inputs and outputs that were read or written in the working directory now point to their folders in the repository: `test1.py`, `test2_h18.py`, `latcorr.py`, `latcorr3.py`, `latcorr5.py`, `depth_vs_xrs.py`, `transition_layer_hilat.py`, `hotwarm.py` (reads `code/rerun/out_v2/`, identical to the former `hollow-networks/data_v2/`), `make_figs_pr.py`.
3. `dT_hosting.py` writes `craters_dT.csv`, the input of `transition_layer.py`, which was previously passed by hand. `make_figs_pr.py`, `slope_fig.py`, and `lifecycle_fig_v5.py` also write the data behind each figure (`figure_data.py`); the figures are unchanged.
4. New scripts recreate results first computed inline, and each reproduces the archived values exactly: `depth_structure_extras.py`, `depth_vs_xrs_tiles.py`, `depth_slope_eqface.py`, `fig3_insolation_v6.py`.
5. `georeference/georeference_bickel.py` re-implements the camera-model georeferencing documented in `GIS/vector/bickel2024_hollows_corrected_README.md`; the original per-image driver script was not saved. Its reader and the comparison script are tested, but the ISIS steps have not been run from this script. Before the code is cited, run it on a few images with ISIS and check the result with `georeference/compare_to_archive.py` (median difference below 0.05 km).

Scripts kept for completeness but not needed for the paper's numbers: `G/scripts/make_figs_v2.py` and `lifecycle_fig.py` (GRL draft figures), `latcorr2.py` and `latcorr4.py` (superseded by `latcorr3.py` and `latcorr5.py`), `peak_wall_test/nested_test*.py` (earlier nested-crater tests), and `code/rerun/latitude/province/` (its analysis script read temporary files from an earlier session and does not run as archived).

## Licenses and credits

Code is under the MIT License (`LICENSE`). Derived data are under CC BY 4.0 (`LICENSE-DATA.md`), which also lists the source data sets. The derived files redistribute fields from the hollow catalog of Bickel et al. (2025), the crater catalog of Herrick et al. (2018), the basin catalog of Orgel et al. (2020), the hollow depths of Blewett et al. (2016), and the MESSENGER X-Ray Spectrometer maps (Nittler et al., 2020); please cite them with this record.
