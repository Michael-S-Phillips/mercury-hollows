# g6_remaining: the remaining catalog-dependent items rerun on the v3 and v3b catalogs (6 Nov 2026)

This folder covers the items that groups g1, g2 and g3 left out: the Stage 2/3 classification against solar incidence (S8), the color
analysis with the screened hollow mask (S9, S16, main-text Sections 3.7 and 3.8), the knob color tests (workshop), the DEM aspect check,
the check of `peak_wall_test/nested_test2.py`, and a scan of v6.3 and supplement v13 for catalog-dependent numbers that no rerun covered.
The scan found the reference-size checks of S2, which are rerun here as well. The catalogs are those of `code/rerun/out_v3`
(v3: 368 spurious detections removed; v3b: also 232 uncertain ones). No existing file was changed, nothing was written into
`manuscripts/`, and every Python call ran with `PYTHONDONTWRITEBYTECODE=1`, so no bytecode was written into the original folders.

**Every number, with its v2 record, v3 and v3b value and whether a conclusion changes: `changes_g6.csv` (53 rows), written by
`build_changes.py` from the outputs below.**

## Files
| Folder / file | What it holds |
|---|---|
| `prep_inputs.py`, `inputs/` | Hollow masks for the color and knob reruns (`hollows_mask_v3.csv`: 6,363 rows; `_v3b`: 6,131; the 64 v2 crater traces stay in the mask, as in the original 6,731), and the aspect host assignments filtered to each catalog (`aspect_host_assignments_v3.csv`: 5,104 of 5,409 rows; `_v3b`: 4,918) |
| `make_copies_color.py`, `make_copies_knobs.py`, `make_copies_aspect.py`, `make_copies_refsize.py` (+ `.log`) | Make every copy with exact string replacements, each asserted to occur once; each copy's first line names its source and what changed |
| `stage/stage_incidence.py`, `stage/out_<cat>/` | Item 1 (new script; see below) |
| `color/` | Copies of `analysis_workshop/color_mdr/{mdr_georef, extract_tiles, merge_spectra, pca_color, diagnostics, depth_color_lrm}.py`; `check_extract_v2.py` (v2 tile check); `depth_sensitivity.py` (new, three S9 sensitivity fits); outputs in `color/out_<cat>/` (`tiles_qc` there links to `_scratch/tiles_<cat>`) |
| `summary_model/` | Copy of `summary_model/step9_color.py` (color tables from `color/out_<cat>`), with g1's copies of `common.py`, `step7_improved.py`, `cached_pool.py`, `run_step.sh` and `boot_failure_check.py`; `summary_model/<cat>/` links to g1's crater tables; outputs in `summary_model/<cat>/out/` |
| `knobs/` | Copies of `melt_source_tests/knobs/{common, t1_knob_color, t1_peak_color, t1_montes_color, summarize}.py`; outputs in `knobs/out_<cat>/` |
| `aspect/` | `prep_aspect_inputs.py` (builds the 8 ppd DEM, below), copy of `code/rerun/latitude/aspect_robustness/aspect_robustness.py`, outputs in `aspect/out_<cat>/` |
| `nested/search_nested_test2.py`, `nested_test2_search.json` | Item 5 |
| `refsize/` | Copies of `completeness_ref/{ref_size_check, ref_size_insolation}.py` and `completeness_sfd/{hosting_fits, joint_thinned}.py`; outputs in `refsize/out_<cat>/` |
| `audit/scan_numbers.py`, `audit/numbers_scan.csv` | Item 6: every number in v6.3 and v13 with its section and whether a changes table contains it |
| `_scratch/` | MDR tile extractions for v3 and v3b (1.2 GB each), per-file links to the original v2 tiles, the emptied v2 check files, and the 8 ppd DEM (`dem8_gis.npz`). Files cannot be deleted from here; `_scratch/` can be removed once the results are accepted |

## Check of the copies on v2
Each copied chain was run on the v2 inputs before it was trusted.
- **Color.** `extract_tiles.py` reproduces all 54 original tile files array for array (`color/out_v2/extract_check.json`; the compared copies were then emptied to save space). On the original tiles, `merged_raw.csv`, `crater_color.csv`, `hollow_color.csv`, `pca_record.json`, `host_vs_nonhost.csv`, `diagnostics_record.json` and `depth_color_lrm_record.json` are byte-identical. `step9_color.py` reproduces `summary_model/out/color_record.json` to a relative difference of 1e-6 or less (differences in the convergence of the iterative fits; identical at the printed precision). `depth_sensitivity.py` gives the three S9 values (0.90, 0.89, 0.91).
- **Knobs.** The three `t1_*` outputs and `knob_tests_record.json` are byte-identical.
- **Aspect.** With the rebuilt DEM, `aspect_robustness_record.json` is reproduced to a relative difference of 1e-9 or less, including the DEM part (agreement 0.789 / 0.638, median 28.3 deg).
- **Reference size.** `ref_size_record.json`, `ref_size_insolation.json` and `hosting_fits.json` are byte-identical, and `joint_thinned.json` agrees to 5e-13.
- **Stage.** The script behind `stage_vs_incidence.json` is not in the connected folder (see item 1). The rebuilt release fit gives 2.03 per 10 deg against the recorded "x2.0".

## Results
**1. Stage 2/3 classification against incidence (S8).** Neither `stage_vs_incidence.json` nor its script is in `hollows-emissivity`.
`PROJECT_STATE.md` (Section 9.7) places them in the proposal folder and gives the model: the 3,268 labeled detections of the Bickel
et al. 0.8 release, with |lat|, pixel scale and emission as controls, odds x2.0 per 10 deg. `stage/stage_incidence.py` rebuilds this fit
(2.03, 1.88-2.20) and fits the same model to the labeled hollows of each catalog, matched as in `stage_test.py`: 1.93 (1.50-2.48) on
2,468 hollows (v2), 1.93 (1.49-2.51) on 2,405 (v3) and 1.93 (1.49-2.50) on 2,372 (v3b). The incidence term of the S8 model itself gives
1.97 in all three. "Double per 10 deg" holds. The Stage 1 fractions and the 44 K step were already rerun by g2; its outputs give
56-71% / 96-100%, a step at 44 K of 3.05 (1.50-6.20) and an AIC lower by 27.5 in v3 (55.5-71.1% / 95.9-100%, 3.07, 27.6 in v3b).

**2. Color with the screened mask (S9, S16, 3.7, 3.8).** The screened mask frees pixels in 291 (v3) and 408 (v3b) crater interiors, and
39 craters reach 20 valid pixels, so the PCA uses 31,412 craters (99.55%) instead of 31,373 (99.4%). No measured hollow's surroundings
change in v3, and only 3 change in v3b, by at most 0.003 SD. Every S9 number and every color number of Section 3.7 is unchanged at the
printed precision: 88% / 5% of the variance, r = 0.93, rho 0.23 / -0.22 / 0.29, depth factors 0.90 (0.86-0.95), 1.01 (0.97-1.06), 0.92,
0.91, 0.90, 0.89, 0.91, 0.78-0.82, rho 0.65 and -0.46, 0.93 (0.89-0.97) and 1.08 (0.91-1.28), 17% and 76%. In S16 the color subset has
442 (v3) and 366 (v3b) hosts, with 4 small hosts lost to the mask instead of 43. The S16 color numbers differ from g1's v3/v3b values,
which used the v2 mask, only in the second decimal: the low-reflectance share with color is 5.5% (2.4-10.7%) in v3 and 4.7% (0.6-11.2%)
in v3b (main text "6% (3-10%)"), the gain is 3.2 (1.4-7.3) and 2.6 (1.2-6.5) points, and the reduction of the regional increment is 0.0013
and 0.0007. In v3b the share interval no longer reaches below zero (g1: -1.0-10.1%), so g1's v3b "yes" for that row comes from the old
mask. As g1 found, nesting moves ahead of size among the other groups.

**3. Knob color tests (workshop).** Unchanged. Knob interiors are brighter and bluer than their matrix in all three catalogs: dln R750
+0.061 [0.056, 0.067], dln(R1000/R480) -0.033 [-0.037, -0.029], 82% brighter and bluer, 8% darker and bluer. The peak, wall and Montes
contrasts change by 0.001 or less. Test 2 (hollow density) comes from g3's outputs; in v3b its bootstrap intervals differ slightly from
g3's because the summary script draws both tests from one random stream (knobs >= 5 km, block interval 0.77-28.6 against g3's 0.75-27.6;
both include 1).

**4. DEM aspect check.** `dempath.txt` names `vfb855e76_mercury_grids_8ppd.npz` on another machine (`.claude-science` artifacts). No file of
that name, and no `.npz` with a DEM array, is in `hollows-emissivity`. The search did find the USGS 665 m DEM (`data/DEMs/`) and its 8 ppd
block average (`GIS/raster/elevation_8ppd.tif`). The GIS grid starts at 180 deg E, whereas the script indexes columns from 0 deg E (the 665 m
DEM has central meridian 180 deg), so `prep_aspect_inputs.py` rolls it by 1,440 columns and checks it against the 665 m DEM (r = 0.99999998).
With this grid the v2 record is reproduced exactly, so it is the grid the script used. On v3 / v3b, the number quoted in the paper, the
equator-facing fraction at 15-30N with Caloris excluded, is 0.44 (0.41-0.48) per hollow and 0.53 (0.40-0.65) per crater, 731 hollows in 47
craters (v3), and 0.44 (0.40-0.47) and 0.45 (0.32-0.58), 717 in 41 (v3b), the same as g2's `wall_fractions.py`. The DEM agreement (not
quoted) stays 0.79.

**5. `nested_test2.py`.** None of its numbers appear in v6.3 or v13. `nested/search_nested_test2.py` searched both files for all 298
values of `nested_test2_record.json` (counts, k/N strings, rates as fractions and percentages, coefficients, p values, exp(coefficients));
the 583 string hits are all coincidental (generic values such as 0.1 or 1.2 in other contexts). Not rerun.

**6. Numbers no group had rerun.** `audit/scan_numbers.py` lists the 3,709 numbers of v6.3 and v13 with their sections and marks those
contained in a changes table; the remaining numbers were read section by section. Most are not catalog-dependent (thermal and sulfur
calculations, Blewett et al. depths, crater-graph counts, image geometry). The catalog-dependent ones not covered before were:
- **S2, Reference size (rerun here, `refsize/`).** In v3 and v3b, completeness for 0.7 km hollows fits the joint model better than
  completeness for 1 km hollows (log-likelihood -1,357.8 against -1,361.8 in v3; -1,147.2 against -1,153.2 in v3b), so the stated reason
  for using 1 km ("describes where host craters are found best") no longer holds strictly; 1 km is second best. In the containment tests
  1 km now fits better than 2 km (by 10-12) rather than "about equally". The ranges become: nesting 1.68-1.76 / 1.78-1.87, vents 3.6-5.2
  / 3.3-5.4, fresh 2.9-4.3 / 2.9-4.5, hot longitudes 2.2-2.6 / 2.1-2.6 within 15 deg, 2.6-3.1 / 2.4-2.9 at 15-30N, 1.49-1.69 / 1.42-1.62 at
  30-45N, and 1.36-1.41 / 1.33-1.38 per 10 K. With the size-weighted measure E: 1.75 (1.42-2.17), 4.3 (2.5-7.2), 3.3 (2.1-5.2) in v3, and
  for the hot longitudes 2.0 (1.4-3.0), 2.3 (1.4-3.8) and 1.45 (0.90-2.33). The 30-45N excess is therefore not significant with 1 km,
  2 km, 5 km or E in v3 (and also with 0.3 km in v3b), which contradicts "loses significance only with completeness for 2 or 5 km".
  Exponents 1.3 and 1.9 change the E odds ratios by up to 4.2% (v3) and 5.0% (v3b), and the joint-model shares by up to 1.7 and 2.1
  points instead of 1.2.
- **Vent depth (3.4, S15).** 1.10 (0.96-1.28), 13 hosts, 405 hollows: unchanged in g1's `test8_vents` outputs, but not tabulated by g1.
- **S8 wall slopes ("65% of wall area has regional slopes < 5 deg").** The script was not found. The wall area below 5 deg in the
  terrain test is 68.8% (352 hosts), 69.1% (321) and 69.0% (262), so the screen barely moves it.
- **Not rerun: S1 georeferencing checks** (97.7% of centers inside their image; repeat detections 0.41 km apart against 1.90 km; 52%
  within 0.5 km; 1,725 image pairs; Table S1). They describe the positions of all detections, including those the screen removes, and no
  conclusion depends on them.

## Conclusion changes (changes_g6.csv, conclusion_change = yes)
- S2 "Reference size": completeness for 0.7 km hollows now fits the joint model better than 1 km (v3, v3b); the 30-45N hot-longitude
  excess is not significant at 1 km or with E (v3, v3b; the 1 km case was already flagged by g2).
- S16 / 3.8: the order of the other groups with color added (nesting ahead of size), as g1 found.
- Not a change: g1's v3b "yes" for the low-reflectance share with color does not hold with the screened mask (interval 0.6-11.2%).

## Caveats
- **Bootstrap failures.** As g1 found, the cloglog fits of `step9_color.py` fail in some bootstrap draws: 5-6 of 60 in v3 and 17-20 of 60
  in v3b (1 in v2). These widen the share intervals; `summary_model/<cat>/out/boot_failure_check.json` gives intervals without them (v3
  low-reflectance share 2.5-10.8%, v3b 2.2-11.1%).
- **Stage fit.** The release fit is a reconstruction from the description in `PROJECT_STATE.md`; it matches the recorded "x2.0" but could
  not be checked against the original JSON.
- **DEM.** The 8 ppd grid is rebuilt from `GIS/raster/elevation_8ppd.tif`; its agreement with the record is exact, but it is not the
  file named in `dempath.txt`.
- **Knob PCs.** The v3/v3b knob runs use the v3/v3b PCA loadings; the quoted contrasts (R750, R1000/R480) do not depend on them.
- **Completeness inputs.** The reference-size checks use the unchanged v2 completeness grids and `crater_completeness.csv`, as the other
  groups did (g3: a detection model refit on v3 pairs would raise completeness 12-14% everywhere).
- **Overlap with g4.** While this folder was being written, `../g4_completeness` appeared beside it. It also reruns the S2 reference-size
  checks, and also with completeness grids rebuilt from the screened catalogs. Its v3 values with the unchanged grids equal the ones here
  (for example E nesting 1.75 (1.42-2.17) and the same joint-model log-likelihoods), so for S2 its rebuilt-grid results take precedence.
- **Disk.** `_scratch/` holds 2.4 GB of tile extractions that could not be deleted from here.

## depth_color_fits.py (9 Nov 2026)
color/depth_color_fits.py runs the depth part of summary_model/step9_color.py on the screened-mask hollow colors
(color/out_v3/hollow_color.csv) and writes color/out_v3/depth_color_fits.json. The factors equal the version-2-mask values to three
decimals (|lat| control 0.903, 0.862-0.945). Supplement v25 final quotes S9 from this record and depth_sensitivity.json.
