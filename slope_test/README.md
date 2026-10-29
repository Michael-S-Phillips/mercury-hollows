# Hollows and regional slope

Question: if overburden sets the depth of collapse, slope should matter. On a slope the normal stress falls as cos²θ and the shear stress rises. Above the friction angle of the devolatilized breccia (about 30–40°), the uncemented column should slide instead of collapsing in place.

Topography: the USGS 665 m MESSENGER global DEM (`../data/DEMs/..._665m_v2_cog.tif`). The east-west spacing is scaled by cos(latitude).
- s1: central differences on the native grid, a 1.33 km baseline.
- s3: after a 3 × 3 mean filter, about 2.7 km.

## Main limitation

The DEM is too smooth to reach the friction angle. Inside the 352 host craters 30–250 km across:
- 80% of the area is below 5°, and 0.5% is at 20° or more.
- 65% of the wall area is below 5°.

A cutoff near 30–35° cannot be tested with this DEM.
- At 25° or more there are 557 km² of crater area and no hollows.
- The crater-average density predicts 3.9 hollows there.
- The density trend at 15–25° predicts about 1.2, so the absence is not informative.

## Occurrence (slope_occurrence.py → occ_by_crater_bin.csv; occ_analysis.py → occ_analysis_record.json)

Method:
- Same hosts, terrain classes, superposed-crater mask, and completeness weighting as `../peak_wall_test/peak_wall.py`.
- Observed and expected hollow counts are compared within each crater, so differences between craters cancel.
- 95% intervals come from bootstrapping over craters.

Hollow density per area falls with DEM slope. Observed/expected (s1):

| Slope | Observed / expected | 95% interval |
|---|---|---|
| < 5° | 1.15 | 1.11–1.19 |
| 5–10° | 0.69 | 0.54–0.85 |
| 10–15° | 0.40 | 0.27–0.53 |
| 15–20° | 0.29 | 0.09–0.57 |

- Hollow fields show the same trend, weaker: 1.09, 0.93, 0.58, 0.40.
- The s3 baseline gives the same result.
- A Poisson GLM with crater fixed effects agrees: rate ratios against < 5° are 0.57, 0.31, 0.22.

The decline depends on which way the slope faces. Observed/expected at 5–10°, 10–15°, and 15–20°:

| Facing | 5–10° | 10–15° | 15–20° |
|---|---|---|---|
| Equator | 0.96 (0.64–1.26) | 0.71 (0.34–1.14) | 0.67 (0.08–1.66) |
| East or west | 0.70 | 0.31 | 0.19 |
| Pole | 0.37 (0.24–0.53) | 0.27 | 0.08 |

A mechanical limit acts the same on every aspect. The deficit on resolved slopes therefore follows sunlight, or illumination-dependent detection, not slope itself. Equator-facing slopes keep near-average density up to 15–20°; the only hint of a drop is at ≥ 20° (2 hollows against 4.8 expected).

## Depth (depth_points_slope.py on the Mac; subsolar_for_depths.py and depth_slope_geometry.py in the cloud → depth_slope_record.json)

Shadow-length depths assume level ground. If the ground under the shadow tilts by β along the shadow (positive when facing the Sun), d_true = d_meas × (1 + tan β / tan e).
- Solar geometry at each hollow comes from the image's sub-solar point in the CDR index. Its incidence matches the index center incidence to within 0.11° (SD).
- DEM slopes at the 2,518 depth sites are small: median 1.6°, 90th percentile 6°.
- The correction factor has a median of 1.00, and 90% of values fall between 0.83 and 1.27.

Results:
- **The DEM slope is not the slope under the shadow.** Measured ln depth against ln(correction) gives −0.10 ± 0.15, where −1 is expected if the DEM slope were the slope under the shadow. Measurements were probably taken on locally level floors.
- **The headline depth results are unchanged by the correction.**
  - Temperature slope: −288 ± 100 K becomes −284 ± 103 K.
  - Fold change: 1.53 becomes 1.52.
  - Poleward of 60°: 0.76 becomes 0.78 (0.69–0.87).
- **Depth against slope cannot test the normal-stress prediction here.** The slopes are too gentle.
  - Coefficient on ln(1/cos²θ): 2.4 ± 2.0, against +1 predicted by normal stress alone and 0 for no slope effect.
  - Depth relative to 0–2° slopes: 1.03 at 2–4°, 1.04 at 4–6°, 1.11 (0.96–1.30) at 6–10°, and 1.08 (0.89–1.31) at ≥ 10°.
- **Facing direction (tentative).** On slopes ≥ 2°, hollows facing the equator are deeper than those facing the pole: coefficient on cos(facing relative to equator) 0.096 ± 0.042, p = 0.02 (1,053 hollows, 295 images), a ratio of about 1.2 from fully pole-facing to fully equator-facing. On slopes ≥ 4° the coefficient is 0.11 ± 0.06 (p = 0.07). The residual shadow bias would make sun-facing hollows look shallower, so this direction is conservative.

## What would test the friction cutoff

Slopes at the scale of hollow fields (≤ 0.5 km baseline) are needed:
- MLA profiles in the north.
- MDIS stereo or photoclinometry DEMs at hollow sites.
- Future BepiColombo topography.
