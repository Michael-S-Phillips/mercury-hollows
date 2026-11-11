# Mercury hollows: code and derived data for a hollow formation model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24121932.svg)](https://doi.org/10.5281/zenodo.24121932)

This repository holds the analysis code, derived data, analysis records, hand-measurement data, and figure scripts behind the paper "Toward a Geologically Plausible Model of Hollow Formation on Mercury" (Phillips, submitted to *Planetary Research*). 

The crater containment graph (craters that formed inside older craters, and an interactive map) is also a separate public repository with its own DOI: [mercury-crater-graph](https://github.com/Michael-S-Phillips/mercury-crater-graph), version 0.3.0, which reproduces the graph results of this paper.

## Running the scripts

The analyses ran in a working folder mounted at `~/mnt/hollows-emissivity`, and the scripts in this repository are byte-identical to the scripts that wrote the archived records, except for one line of `R/g4_completeness/g2/threshold.py` changed in v0.3.1 and the folder and variable names changed in v0.4.2 (see Versions). Most scripts therefore locate their inputs under `~/mnt/hollows-emissivity`. Link the repository there before running anything:

```bash
bash setup_root.sh        # creates ~/mnt/hollows-emissivity -> this repository (stops if that path already exists)
pip install -r requirements.txt
bash run_all.sh           # reruns the analyses that need only the data in this repository (a few minutes)
```

Scripts kept from version 0.1.0 instead find the repository root from their own location or from `$MERCURY_HOLLOWS_ROOT`. Twelve scripts also contain absolute paths of the machines where earlier versions ran (`/Users/...`, `/home/claude/...`, `/mnt/user-data/...`, `/sessions/...`); they are listed in `provenance/hardcoded_paths.txt`. In these scripts the absolute paths belong to helpers that copied scripts between rerun folders, to optional figure panels, to defaults that environment variables override, or to blocks that are off by default; none of the steps in the reproduction check below depends on them.

## What the paper uses

**Catalog.** The paper uses the screened hollow catalog v3: 368 hollow entries that visual screening classed as crater bowls, crater walls, or other features are flagged as crater traces, together with the 64 detections larger than 10 km that outline whole craters. `code/rerun/out_v3/README.md` describes the screen; `code/rerun/out_v3/build_v3_catalog.py` rebuilds the catalog files from the v2 products and the visual classes.

**Detection model and completeness.** Completeness comes from the detection model refit on the screened catalog (tag **v3r**: `code/rerun/out_v3/rerun/g3_catalog/sfd/out_v3_refit/`, grid `g4_completeness/completeness_10ppd_v3.npz`). Results tagged v3b or v3br use the stricter catalog that also drops uncertain entries; they are sensitivity checks.

**Terrain classes.** Peak, floor, and wall classes come only from hand tracings of the host craters (`g12_tracer` and `g17_tracer_round2` hold the tracing pages' inputs; `g13_traced`, `g18_traced_r3`, and `g19_knobs_traced` hold the exported tracings and the analyses). The traced-class scripts import two modules of automated curvature classes (`g10_threshold/curv_classes.py`, `curv_classes_v3.py`) and one gap-filling module (`g7_terrain_fix/terrain_classes.py`). Those modules assign classes to hosts that were not traced, in intermediate tables only; every reported terrain number uses the traced hosts alone (`rows_tracedonly_r3_v3.csv`). The S8 share of wall area with regional slope below 5° (61%) also uses the traced hosts (`R/g18_traced_r3/terrain/out_v3/peak_wall_record_rows_tracedonly_r3_v3.json`, read by `R/g27_text_numbers/text_numbers.py`); supplement v25 had taken it from the automated relief-threshold classes, whose record has left the archive.

## Layout

Folders keep the layout in which the analyses ran, so scripts find one another's outputs by the same relative paths. `R/` below is `code/rerun/out_v3/rerun/`, `I/` is `manuscripts/insolation_thermal/`, and `P/` is `manuscripts/PlanetaryResearch_hollows/`. `I/` and `P/` hold only the submitted manuscript and the scripts, inputs, and records that the reproduction check and the provenance tables use; drafts and superseded analyses are in v0.2.1.

| Folder | Contents |
|---|---|
| `georeference/`, `GIS/vector/`, `data/bickel-2024/` | Camera-model positions of the 19,110 hollow detections of Bickel et al. (2025) and the position-error record (Section S1) |
| `code/rerun/build_hollow_products.py`, `code/rerun/out_v2/`, `hollow-networks/` | Merged hollows, hollow fields, imaging coverage, the first detection model (v2 products) |
| `code/rerun/out_v3/` | Screened catalog v3/v3b, host assignments, screening sets and classes (`qc/`), `build_v3_catalog.py` |
| `R/g1_graph_joint/`, `R/g4_completeness/g1/` | Crater containment graph, degradation state (Herrick and Kinczyk), vents, joint model of hosting (`summary_model/`); `g4_completeness` holds the v3r runs |
| `R/g2_insolation/`, `R/g4_completeness/g2/` | Hot longitudes, ΔT_w, latitude, wall orientation, Stage 1 tests, plains hosts, Figure 4; `g4_completeness/g2/out/v3r/` holds the records |
| `R/g3_catalog/` | Size distribution and detection-model refit (`sfd/`), catalog counts, Caloris penetration test (`caloris/`) |
| `R/g4_completeness/` | Completeness grids from the refit model, Figure S2, and the reference-size checks (`s2/`) |
| `R/g5_figures/`, `R/g8_figS5/` | Figures redrawn on the screened catalog, including Figures 1, 3, S4, S5, and S9 (`g5_figures`) and Figure S6 (`g8_figS5`) |
| `R/g6_remaining/color/` | MDIS eight-filter color of crater interiors (Section S9) |
| `R/g9_draft_reruns/knobs/`, `R/g19_knobs_traced/` | Hollows on the knobs around Caloris (Section S17) |
| `R/g12_tracer/`, `R/g13_traced/`, `R/g17_tracer_round2/`, `R/g18_traced_r3/` | Hand tracings of host terrain and the peak, floor, and wall tests (Section S4, Figure 2d, Figure S3) |
| `R/g20_regional_clusters/`, `R/g25_composition/`, `R/g22_province_screen/` | Regional excess of host craters over the joint model, the composition group, and the visual screen of the excess regions (Section 3.8, Figure 6, Section S16) |
| `R/g23_smallest_hollows/`, `R/g24_blewett_compare/`, `R/g26_measurement_error/` | Automatic shadow depths, hand depths, the blind check of Blewett et al. (2016) depths, and measurement uncertainty (Section 3.7, Figure 5, Section S9, Figures S7 and S8) |
| `I/scripts/`, `I/depthctl/`, `I/data/`, `I/records/` | Thermal model, plains mask, Blewett depth tables, XRS maps |
| `P/analysis_v3/`, `P/analysis_v7/`, `P/analysis_v8/`, `P/analysis_v9/` | Sulfur release and vapor pressure, transition layer, collapse depth, model-logic calculations, first graph versions |
| `P/analysis_workshop/` | Volatile fits, sulfur redox, Kinczyk classes, vent checks, depth pilot, province basins, Caloris |
| `P/figscripts/` | Life-cycle figure (Figure 7, `lifecycle_fig_v12.py`) and base temperature against φ (Figure 8) |
| `P/submission/` | The submitted manuscript and supplement (`Hollow_Formation_Model_Paper_final.tex` = v6.17, `supplementary-materials_final.tex` = v29) and their figures (`figures/`) |
| `peak_wall_test/` | DEM morphometry of the cataloged craters (`morph_all.csv`) and the Caloris children |
| `R/g27_text_numbers/` | Numbers that the text computes from stated inputs (S2, S8, S13, Section 3.6), recorded in `text_numbers.json` |
| `R/g28_fig2a/` | Figure 2 with panel a on the MDIS mosaic (`layers.py` needs the external 166 m mosaic; `make_figs_pr_v10.py` draws the figure) |
| `provenance/` | Where each number and figure comes from (next section) |

## Where each result comes from

`provenance/main_text.csv` (151 rows) and `provenance/supplement.csv` (118 rows) map every quantitative result of the manuscript and the supplement to the script that computed it and the record that holds it, with the JSON key or table column where one was found. The column `value_found` says whether the printed value appears in the record exactly (`yes`), after rounding or a simple derivation (`approx`), or not at all (`no`, with a note). `provenance/figures.csv` gives the script, the data, and any external data behind each of the 17 figures; every submitted figure is byte-identical to its script's output. These tables were compiled on 9 November 2026 by matching the printed values against the records; the row for Equation S1 and Table S21 was added on 11 November 2026, and the table numbers follow supplement v29.

## Reproduction check (9 November 2026)

The following analyses were rerun from this repository, with `~/mnt/hollows-emissivity` pointing at it, and every record they write is byte-identical to the archived record:
- the screened catalog (`build_v3_catalog.py`);
- the containment graph, its tests, and its sensitivity runs (`R/g4_completeness/g1/run_dag.sh v3r main`);
- degradation state and vents (`run_degr_vents.sh v3r`);
- all 22 insolation analyses (`R/g4_completeness/g2/run_g2.sh` with `G2_CAT=v3r`), including Figure 4, which is byte-identical to the submitted figure;
- the size distribution and the detection-model refit (`R/g3_catalog/sfd/run_sfd.sh v3 refit`, `sfd_derived.py`);
- the depth analyses `depth_regress.py`, `within_field.py`, `bench_analyze.py`, `fig_depth_v3.py` (Figure 5, byte-identical), `cmp_analyze.py`, `me.py`, and `no_shadow_check.py`;
- the traced terrain tests (`R/g18_traced_r3/terrain/analyze.py out_v3/rows_tracedonly_r3_v3.csv`) and the knob ratio (`knob_ratio_design.py`);
- the model calculations in `P/analysis_v3` (sulfur charge, FeS vapor), `P/analysis_v7`, `P/analysis_v8`, and `P/analysis_workshop/volatile_fit`;
- the Caloris tests on the refit grid (`G3_CAT=v3 G3_GRID=v3r` with `R/g3_catalog/caloris/caloris_v9.py` and `caloris_v9_controls.py`, records in `caloris/out_v3r/`), which differ only in the recorded paths of their inputs;
- the coefficients of the joint model, Equation S1 and Table S21 (`R/g25_composition/coef_table.py` with `SM_DIR` set to `R/g4_completeness/g1/summary_model/v3r`; rerun on 11 November 2026);
- the depth regressions on the automatic depths (`R/g23_smallest_hollows/depth_paper.py`, which reads the external map of low-reflectance material of Klima et al., 2018), the depth fits on the screened-mask colors (`R/g6_remaining/color/depth_color_fits.py`), and the text numbers (`R/g27_text_numbers/text_numbers.py`).

A rerun on 11 November 2026 (v0.3.1) with Python 3.11.13, NumPy 2.3.3, pandas 2.3.2, SciPy 1.16.2, statsmodels 0.14.6, and Matplotlib 3.10.6 completes every step above except `depth_paper.py`, whose external map was not at hand. Its records agree with the archived records to floating-point precision (the last digit) rather than byte for byte, with one exception: in `R/g4_completeness/g2/out/v3r/threshold.json`, the pointwise 95% band around the spline (`spline_shape` `lo` and `hi`) differs by up to 2.4%, while every estimate in the record agrees to 10⁻¹¹. The paper does not use the band. Figure 4 differs from the submitted figure only in how its text is rendered.

Partial checks:
- **Joint model with the composition group** (`R/g25_composition/comp_model.py`). A run with 2 bootstrap draws reproduces every point estimate (R² 0.226, the shares, and the coefficients) to 10⁻¹⁵. The archived intervals come from 100 draws, which were not rerun with `comp_model.py`. `coef_table.py` fits the same model to the same 100 block resamples and reproduces the composition factors and intervals of `comp.json` to 10⁻¹⁴.
- **Records that need Cantera or external data.** `P/analysis_workshop/sulfur_redox/redox_check.py` (Cantera; `redox_check_record.json`) and `I/scripts/nsp_mask2.py` (the northern smooth plains lineament map; `nsp_mask2_record.json`) were not rerun here.

Not rerun, because the steps need external data: georeferencing (ISIS and MDIS EDRs), the catalog screen and the automatic shadow depths (MDIS EDRs), the completeness grid (MDIS CDR index), DEM morphometry, terrain measurements and hollow positions on peaks (USGS 665 m DEM), MDIS color (MDR tiles), and the map figures (MDIS mosaics). Their outputs are archived.

## Hand-measurement and visual-classification data

The author measured and classified features by hand in browser pages; the exported results are data in this repository.
- **Terrain tracings** of 111 host craters: `R/g13_traced/traces/`, `R/g18_traced_r3/traces/`, and, for the Caloris knobs, `R/g19_knobs_traced/traces/` (raw exports `traces_raw_artifact_*.json` and cleaned `traces_clean*.json`).
- **Hand shadow depths** of 135 hollows: `R/g23_smallest_hollows/bench_shadows_2026-11-09.json`, `bench_depths.csv`.
- **Blind measurement** of 178 points of Blewett et al. (2016): `R/g24_blewett_compare/blewett_cmp_2026-11-09.json`, `cmp_points.csv`.
- **Visual screen of the regions with excess hosts**: `R/g22_province_screen/calls_2026-11-09.json`.
- **Catalog screen**: `code/rerun/out_v3/screen_set.csv`, `screen_classes.csv`, `qc/class_*.csv`, `qc/decide_ids.csv`.

## External data (not included)

| Data | Used by | Source |
|---|---|---|
| Bickel et al. hollow catalog, releases 0.7 and 0.8 confidence (`data/bickel-2024/Bickel_et_al_2024_MDIS_hollows_07_v1.0.csv`, `_08_v1.0.csv`) | `build_hollow_products.py`, `georeference_bickel.py`, the Stage tests in `g2` | Bickel et al. (2025) data release |
| MDIS CDR cumulative index | `build_hollow_products.py`, `g4_completeness/build_completeness.py`, depth scripts | PDS Imaging Node, `msgrmds_2001/INDEX/` |
| MDIS EDR images | georeferencing, catalog screen, automatic depths, province screen | PDS Imaging Node, `msgrmds_1001/DATA/` |
| MDIS MDR eight-filter mosaics (`data/mdis_mdr/`, fetched by `data/mdis_mdr/_download.sh`) | `g6_remaining/color/` | PDS Imaging Node |
| USGS MESSENGER global DEM, 665 m/pixel (`data/DEMs/`) | `peak_wall_test/measure_all.py`, `nested_build.py`, traced terrain and position builders, knob scripts | USGS Astrogeology (Becker et al., 2016) |
| MDIS global mosaics, 166 m/pixel monochrome and 665 m/pixel enhanced color | map figures, tracing pages | USGS Astrogeology (Denevi et al., 2018) |
| Low-reflectance material map (`data/Concentrated_LRM_560_params_BW.tif`) | `depth_paper.py`, composition steps | Klima et al. (2018) |
| Volcanic vent catalog (`data/data/MER_Vent_Locations_Classifications.csv`) | vent tests | Pegg et al. (2021) |
| Northern smooth plains lineament map (`data/data/NSP_map_data/`) | `I/scripts/nsp_mask2.py` | Klimczak et al. (2025); the resulting mask is included |
| Kinczyk et al. (2020) degradation table | `P/analysis_workshop/kinczyk_degradation/` | Kinczyk et al. (2020) supplement; the matched classes are included |
| ISIS 8 with MESSENGER kernels; Cantera (sulfur-redox checks) | georeferencing; `P/analysis_workshop/sulfur_redox/` | USGS Astrogeology; cantera.org |

## Versions

**v0.4.2 (November 2026).** Renames the folder that holds the thermal model, plains mask, Blewett depth tables, XRS maps, and their records to `manuscripts/insolation_thermal/` (alias `I/` in this README). Every path to the folder is changed to match, in the scripts, the provenance tables, and the three records that store input paths (`depth_nested_record.json` in `R/g1_graph_joint/depth_nested/v3/` and `P/analysis_v8/depth_nested/`, and `R/g3_catalog/caloris/out_v2/caloris_v9_record.json`). The variables and environment variables that point to the folder are renamed to match (`INSOL`, `INSOL_DIR`, `INSOL_SCRIPTS`, `INSOL_WORK`, `INSOL_FIGS`, and `G2_INSOLDATA`). No number, figure, or manuscript file changes.

**v0.4.1 (November 2026).** Matches the submitted manuscript v6.17 and supplement v29. v6.17 changes wording only, and both files no longer carry the version comments at the top. v29 adds to Section S16 the fitted thinned occurrence model with every term and coefficient (Equation S1) and a table of the factor by which each term changes the density of hollows, with 95% intervals (Table S21). `R/g25_composition/coef_table.py` (new) computes them and writes `coef_table.csv` and `coef_table.json`; `run_all.sh` reruns it, and `provenance/supplement.csv` has a row for it. In `provenance/supplement.csv`, the rows for Tables S10 to S20 gave table numbers one too high; they now match the supplement. The figures are unchanged. The Zenodo record labeled v0.4.0 (10.5281/zenodo.24146453) was archived by mistake from the v0.3.1 commit and holds the same files as v0.3.1; v0.4.0 has no other release.

**v0.3.1 (November 2026).** Two steps of `run_g2.sh` failed in v0.2.0 to v0.3.0. `latcorr5.py` (Section S6) reads `code/rerun/latitude/incidence_null/out_inc60_75/coverage_10ppd.npz`, which v0.2.0 removed; the file is restored from v0.1.0, and the record it writes is byte-identical to the archived one. `threshold.py` (Section 3.5 and S8) read the design of a formula model as `model_spec`, which statsmodels 0.14 calls `design_info`; it now reads either. The vent catalog is credited to Pegg et al. (2021).

**v0.3.0 (November 2026).** Removes from `manuscripts/` the manuscript drafts (v6.13 to v6.15, supplement v24 to v27), the document builders of an earlier manuscript, and the analyses and figures that neither the reproduction check nor the provenance tables use (195 files). The submitted manuscript and supplement, and every script, input, and record that `run_all.sh` or `provenance/` names, are unchanged. The removed files remain in v0.2.1.

**v0.2.0 (November 2026).** Matches the submitted manuscript v6.16 and supplement v28, held as `manuscripts/PlanetaryResearch_hollows/submission/Hollow_Formation_Model_Paper_final.tex` and `supplementary-materials_final.tex` (the submission folder holds only the current version; v6.13, v6.14, v6.15, v24, v25, v26 and v27 are one level up). v6.15 redraws Figure 2a on the MDIS mosaic (`g28_fig2a`) and letters the panels of Figure 7 (`lifecycle_fig_v12.py`); v26 takes the S8 wall-area share from the traced classes, and v27 fixes the panel titles of Figure S1 (`analysis_v3/figS_positions_v4.py`). v6.16 and v28 rename the solar-day average "mean surface temperature" (Figure 5 v4, `g23_smallest_hollows/fig_depth_v4.py`) and add the pressure-driven flow of S2 to S13 (`g29_pressure_flow/pressure_flow.py`). The corrections came from the provenance cross-check: crater counts and the Caloris tests on the refit completeness grid (`g3_catalog/caloris/out_v3r`, Figure 2c), the S2 pair ranges and rounding, the S9 color fits on the screened mask (`g6_remaining/color/depth_color_fits.py`) with Figure S9 redrawn from it, and new records for numbers that had none (`g27_text_numbers`, `redox_check_record.json`, `nsp_mask2_record.json`); `depth_paper.py` now writes every key of its record. The provenance tables describe v6.15 and v27. Adds the screened catalog and every rerun on it (`code/rerun/out_v3/`), the refit detection model and completeness grid, the hand-traced terrain classes, the automatic and hand shadow depths and their uncertainty, MDIS color, the Kinczyk degradation classes, the joint model of hosting with the composition group, the regional excess, the knobs, the provenance tables, and the submitted figures. The automated relief-threshold terrain classes and the analyses that the paper no longer uses (slope tests, the version-6 figures and their data tables, superseded latitude records) are removed. `build_v3_catalog.py` is new: it rebuilds the screened catalog files, which were first made by commands run outside a script, and it reproduces them exactly.

**v0.1.0 (29 October 2026).** First archive, for manuscript draft v6.

## Licenses and credits

Code is under the MIT License (`LICENSE`). Derived data are under CC BY 4.0 (`LICENSE-DATA.md`), which also lists the source data sets. The derived files redistribute fields from the hollow catalog of Bickel et al. (2025), the crater catalog of Herrick et al. (2018), the basin catalog of Orgel et al. (2020), the degradation classes of Kinczyk et al. (2020), and the hollow depths of Blewett et al. (2016), so the record credits them. Much of the analysis code was written with Claude (Anthropic) and checked by the author.
