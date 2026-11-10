# Depth-expansion pilot: automatic hollow depths from shadows in raw MDIS NAC frames

Pilot run on 2 November 2026. The pilot tests whether hollow depths can be measured automatically from shadow lengths in raw
MESSENGER MDIS narrow-angle (NAC) EDR frames at the Bickel et al. (2024) detection boxes. It validates the measurements
against the manual shadow-length depths of Blewett et al. (2016), and then measures new depths for hollows in small host
craters and poleward of 60 degrees.

**Bottom line.** The validation is a marginal pass. Automatic and manual depths correlate well (r of ln depth = 0.70,
n = 353 boxes in 110 images), but the automatic depths run about 24% deeper than the manual depths (median ratio 1.24,
95% interval 1.18-1.27). Poleward of 60 degrees the overestimate is larger (1.62, from only 23 boxes in 12 images). The
method misses many shallow hollows whose shadows are 3 pixels long or shorter. The new depths are therefore useful for
relative comparisons made with the same method. They are not a drop-in extension of the Blewett catalog.

## Data and sources
- MDIS EDRs (PDS MSGRMDS_1001) came from the PDS Imaging Node mirror at
  `https://planetarydata.jpl.nasa.gov/img/data/messenger/MSGRMDS_1001/DATA/YYYY_DOY/EN<id>M.IMG`. The day-of-year folder
  comes from the Bickel 'Earth datetime'. The pds-geosciences.wustl.edu paths for MDIS returned 404. The EDRs are cached on
  the VM (`~/edr_cache`), not in the project. In total, 622 frames were read and 2,285 detections were measured.
- The manual depths are the 2,518 Blewett et al. (2016) depths in `GRL_insolation_hollows/data/test1_rows_inc.csv`.
- The detections come from the Bickel release CSV (0.7 confidence). Camera-model positions and box corners come from
  `GIS/vector/bickel2024_hollows_corrected.gpkg`.
- Host craters come from Herrick et al. (2018) (`test2_craters_h18.csv`), assigned with the depth_nested.py rule.

## Method (scripts in this folder)
1. **Reading the frames (`dep_common.py`).** The script reads the PDS3 label (`^IMAGE`, `LINES`, `LINE_SAMPLES`,
   `SAMPLE_BITS = 8`, `SUB_SOLAR_AZIMUTH`, `SUB_SPACECRAFT_AZIMUTH`, incidence, emission, pixel scale and sub-solar point)
   and loads the 8-bit companded DN array in raw orientation. Binned 512 x 512 frames are handled through the pointer. The
   release box coordinates (normalized 0-1, full frame) sit on the hollow clusters in the raw frames
   (`qc/explore_sun_az.png`).
2. **Sun direction.** The label `SUB_SOLAR_AZIMUTH` matches the release 'Sub sol azimuth' to within 0.01 deg. It is
   measured clockwise from +sample, with lines increasing downward, and it points toward the sun. Three checks established
   this convention:
   - In 24 boxes in 6 frames, inspected visually (`qc/explore_sun_az.png`), the cast shadows lie on the sunward (left)
     inner wall and the sun-facing wall is bright.
   - In 888 boxes (`sun_convention.json`), the sunward half of the box is darker than the anti-sunward half in 80% of
     cases. The mean profile across the box dips to 0.96 of the surroundings at u = 0.1 and peaks at 1.10 at u = 0.8.
   - The NORTH_AZIMUTH convention verified earlier in `feature_survey_hires/nsp_gallery/make_panels.py` agrees.
   In the frames used here, the azimuth is 175.4-179.6 deg (1st to 99th percentile), so the sun is toward -sample and
   shadows run along image rows. The general-azimuth code path was therefore barely exercised. Across the whole release,
   about 5% of frames have azimuths below 159 deg.
3. **Shadow measurement (`shadow_measure.py`, v2).**
   - The script crops the box plus a margin and resamples it on a sun-aligned grid with 0.25 px steps.
   - It sets two brightness levels. B is the surroundings median and S is the 2nd-percentile DN inside the box. The
     threshold is the half-contrast level, T = (B + S)/2.
   - On every line it finds the runs darker than T. A run counts only if it meets all four of these conditions:
     - it starts inside the box;
     - it ends inside the box;
     - it is followed within 3 px by a surface at least as bright as B (the lit far wall of a depression);
     - it does not lie on the outer 10% of lines across the box.
   - The shadow length is the maximum of the per-line lengths, after a 1.25 px running median across lines.
   - The ground length is L = L_px x px x f. The factor f corrects foreshortening along the sun direction. The spacecraft
     azimuth is near 0/180 deg, so the emission tilt lies along the shadow. The camera-model affine confirms this: the ground
     scale along samples correlates with f at r = 0.96.
   - The depth is h = L / tan(i), with i the incidence at the detection on a sphere.
   - The widths are the box extents along and across the sun direction, plus the release polygon's equivalent diameter.
   - The flags are `multi`, `no_shadow`, `edge_sun`, `edge_far`, `not_depression`, `short` and `frame_edge`.
4. **Validation (`validate.py`).** Blewett coordinates sit a median 2.2 km from the camera-model positions. Without
   correction, only 3% of Blewett points fall inside a detection box, and 30% fall outside the frame. Each shared image was
   therefore registered with a translation:
   - An affine lon/lat-to-pixel map is fitted to the box corners (median rms 0.3 px in the accepted frames).
   - An exhaustive shift search (within +-8 km) finds the shift that puts the most Blewett points inside boxes. The depths
     are not used in this search.
   - A frame is accepted when it has at least 3 points, at least 50% of its points are captured, and the best shift is
     unique and within 0.5 km.
   Of the 416 shared images, 130 were accepted (median shift 2.6 km). In those frames, 817 Blewett points matched 469
   boxes, and 434 of those boxes have emission <= 30 deg. Each box's automatic depth (the longest shadow in the box) is
   compared with the deepest Blewett depth in that box.
5. **Targets (`targets.py`).** The base set is detections with confidence >= 0.75 in frames with px <= 20 m, incidence
   60-80 deg and emission <= 30 deg: 1,866 detections in 583 frames. Two target groups were drawn from it:
   - (a) host craters smaller than 20 km: 75 detections, 32 hosts, 42 frames;
   - (b) |lat| >= 60 deg: 203 detections, 134 frames.
   For comparisons, depths are corrected to the Blewett median px (14.24 m) and incidence (67.6 deg) with the paper's OLS
   fit, ln h ~ ln px + ln tan(i), fitted on the 2,518 Blewett depths (b = 0.78, -0.73). The same fit on the automatic base
   set gives 0.91 +- 0.06 and -0.79 +- 0.09, so the correction is appropriate. Intervals come from 2,000 bootstrap
   resamples of images.

## Validation result (`validation.csv`, `validation_record.json`, `validation_sensitivity.json`, `qc/validation_scatter.png`)
| Set (emission <= 30, clean) | n (images) | r(ln) | median auto/manual [95%] | MAD ln ratio |
|---|---|---|---|---|
| v1 (first pass) | 353 (110) | 0.68 | 1.30 [1.25-1.38] | 0.25 |
| **v2 (adopted)** | **353 (110)** | **0.70** | **1.24 [1.18-1.27]** | **0.21** |
| v2, one Blewett point per box | 189 (70) | 0.61 | 1.24 [1.17-1.30] | 0.20 |
| v2, \|lat\| >= 60 | 23 (12) | 0.76 | 1.62 [1.23-1.74] | |

The v1 run failed the criterion (median ratio 1.30). The v1 contact sheets (`qc/val_*_v1.png`) showed two failure modes:
long runs on the box's top and bottom edge lines over dark terrain, and dark runs with no lit wall after them. The two v2
rules address these failure modes. They were chosen from the images and were not fitted to the depths.

The edge convention controls the absolute bias. The half-contrast threshold (50%) gives 1.24. A threshold at 40% gives 1.10
(r = 0.78), and one at 60% gives 1.40. The Blewett measurements behave like a darker edge, about 40%. I kept 50% because it
was set before the comparison, and choosing 40% would calibrate the method on the validation set.

The yield is 81% (353 of 434 boxes). By Blewett shadow length, the yield is 57% for shadows of 3 px or less, 75% for 3-5 px,
92% for 5-8 px and 100% for longer shadows. The shallow end is therefore truncated (`qc/val_noshadow_v2.png`). Contact
sheets: `qc/validation_sheet_random_v2.png` (24 random pairs) and `qc/validation_sheet_extremes_v2.png` (worst ratios: runs on
crater-wall shading or across dark floors).

**Judgment: a marginal pass.** The point estimate (1.24) and r (0.70) meet the criterion, but the interval reaches 1.27.

## New depths (`new_depths.csv`, `targets_record.json`)
Of the 254 target detections (24 are in both groups), 186 gave clean measurements.

**Small hosts (< 20 km).** There are 51 clean depths in 25 hosts. Fifteen of these hosts have no Blewett depth at all, but
each of the 15 contributes only one hollow.
- The raw median is 23.0 m (IQR 16.6-33.1).
- The corrected median is 34.6 m [26.0-44.6], against 35.7 [33.8-37.9] for automatic depths in hosts of 20 km or larger.
  The ratio is 0.97 [0.73-1.26].
- In the Blewett set, the small-host median is 22.1 [19.3-26.0] (112 depths, 22 hosts), against 25.3 for larger hosts
  (ratio 0.87).
- The calibrated automatic median (divided by 1.236) is 28.0 [21.1-36.1], which is 1.27 [0.89-1.66] times the Blewett
  small-host median.
- In short, the new depths show no detectable small-host effect, and the interval is wide.

**Poleward (>= 60 deg).** There are 148 clean depths in 103 frames, and 63 of those frames are new relative to the Blewett
set.
- The raw median is 21.8 m (IQR 14.3-30.5).
- The corrected median is 34.2 [31.3-37.0], against 35.9 [34.2-38.0] for automatic depths below 60 deg (ratio 0.95
  [0.85-1.04]) and 38.3 for 50-60 deg (ratio 0.89 [0.80-0.99]).
- In the Blewett set, the poleward median is 20.1 [17.8-22.3] against 25.7 below 60 deg (ratio 0.78).
- The calibrated automatic median (divided by 1.236) is 27.7 [25.3-29.9], which is 1.38 [1.21-1.58] times the Blewett
  poleward median. With the poleward validation factor (1.62), the median is 21.1 [19.3-22.8], and the ratio is 1.05
  [0.92-1.20].
- The method over-reads more at high latitude, and that excess would hide a poleward shallowing of the size Blewett found.
  The new depths therefore do not independently test the latitude trend. They are consistent with the Blewett trend once
  that latitude-dependent bias is allowed for.

**Visual QC.** Of 30 random clean target measurements (`qc/targets_sheet_random.png`, `qc/targets_visual_qc.csv`), 12 are
good, 8 are plausible, 7 are doubtful (2 px runs or runs off the main shadow) and 3 are wrong. The wrong cases are
crater-wall shadow or shading on hollows that sit on crater rims (58 m and 105 m) and one run on the anti-sunward side.

## Caveats
- **Systematic overestimate.** The automatic depths run about 24% deep overall and about 60% deep poleward (a small sample).
  Absolute depths depend on the edge convention. Compare automatic depths with automatic depths, or calibrate them.
- **Truncation.** Shadows of 3-5 px or less are often missed, so the medians are biased high. Comparing the overall
  automatic median with the overall Blewett median gives 1.42, larger than the per-hollow 1.24, and this truncation accounts
  for the difference. Depth scales with px^0.91, so the measurements are resolution-limited, as the manual set is.
- **Boxes, not hollows.** 65-78% of the boxes hold several shadows (`multi`). The depth is that of the longest shadow in
  the box, and the widths are the dimensions of the box (cluster scale), not of single hollows.
- **Registration.** Blewett positions are offset 1-5 km from the camera-model positions. Only 130 of the 416 shared frames
  could be registered, and the matching is at box level.
- **Image processing.** The DN values were not decompanded. The incidence is computed on a sphere, with no local slopes.
  The camera-model positions are uncontrolled (about 0.3-1 km), which matters for host assignment near crater rims.
- **Recommended next steps.** Segment individual hollows inside the boxes, add PSF and decompanding, calibrate the threshold
  on held-out images, and screen out hollows on crater walls and rims.

## Files
- Scripts:
  - `dep_common.py`: loaders, EDR access, frame maps and host assignment.
  - `shadow_measure.py`: the shadow measurement.
  - `validate.py`: register, fetch, measure and stats stages.
  - `sensitivity.py`: the setting variants.
  - `sun_convention.py`: the sun-direction check.
  - `targets.py`: select, run, compare and qc stages.
  - `plot_validation.py`, `qc_render.py`, `qc_validation_sheet.py`: figures and contact sheets.
- Tables:
  - `registration_frames.csv`, `registration_points.csv`: the registration of Blewett points.
  - `validation_measurements.csv`, `validation.csv`: validation measurements and matched pairs.
  - `targets_selection.csv`, `base_measurements.csv`, `new_depths.csv`: target selection and new depths.
- Records: `validation_record.json`, `validation_sensitivity.json`, `sun_convention.json`, `targets_record.json`.
- The first-pass (v1) outputs are kept as `*_v1.*`. `qc/explore_*.png` is early exploration made with v1 code.

## Hollows on circum-Caloris knobs (5 Nov 2026): knob_measure.py
Targets: Bickel detections on knobs of Wright et al. (2020) (supplementary shapefile, 2,011 knobs; matching in
../province_basins/knob_depths.py), in frames with incidence >= 55 deg. There are 26 detections on 14 knobs in 23 frames,
at 12-64 m/pixel.
- Of the 26, 6 are clean, on 4 knobs: raw depths 25-52 m and corrected depths 29-44 m. The corrected median is 36 m,
  against 32 m for the clean automatic base set.
- Visual QC (qc/knob_sheet.png) shows that most of these do not measure hollows on knobs:
  - The clean depths on knobs 51 (3 frames) and 400 are round depressions with bright far rims. They look more like small
    craters than hollows.
  - The clearest hollows on knobs are on knob 235 (Wright et al.'s image EN1045703448M), 239 and 1436. Those boxes are
    dominated by the knob's own shadow (155-180 m "depths") at 29-61 m/pixel.
  - Plausible hollow depths come only from knob 273 (25 m raw at 15 m/pixel) and knob 905 (27 m raw at 24 m/pixel).
- No MDIS NAC frame of 15 m/pixel or better with incidence 55-82 deg covers any of the 32 knobs with detections (CDR index
  search; ../province_basins/wright2020_knobs/knob_site_frames.csv).
- Conclusion: MESSENGER images do not support an automatic depth test on knobs. At best a few manual measurements are
  possible.
