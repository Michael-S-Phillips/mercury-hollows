# g23: shadow depths of the smallest cataloged hollows (9 Nov 2026)

Michael asked for a look at the smallest hollows in the catalog and a shadow depth for each, to test the idea that a column of
charged regolith must be wider than some minimum to collapse to full depth.

## Set
All 21 hollows with equivalent diameter < 100 m in the screened v3 catalog used in the paper (g2_insolation/inputs/v3/hollows_canonical.csv,
crater traces removed). They are 39-95 m across, each detected once, in 15 MDIS NAC frames at 1.3-6.9 m/pixel and incidence 58-80 deg,
at 37-63 N and 221-318 E (none in the province). Most frames come from MESSENGER's final low-altitude images (CN1066-1067...).

## Method
- small_shadow.py: the automatic shadow length of the depth-expansion pilot (analysis_workshop/depth_expansion_pilot/shadow_measure.py, v2)
  on the detection box; ground length L = L_px x px x f (foreshortening along the sun); depth h = L / tan(i), with i the local incidence.
  h_fill = (box extent along the sun) / tan(i) is the depth at which the shadow would fill the box. Writes small_shadow.csv, tiles/<id>_z.png
  and _c.png, sheet_1-4.png.
- small_shadow_v2.py: the same with the box widened toward the sun by 0.75 box lengths, so a shadow the detection box left out is counted
  (2496: the box covers only the lit far half of the depression). Writes small_shadow_v2.csv, profiles.json, tiles/<id>_v2.png,
  sheet_v2_1-3.png (yellow = detection box, magenta = widened box, red = measured shadow).
- Visual classes below are my reading of the tiles (sun from the left in all frames). Michael has not screened these.
- Caveats: in the pilot's validation the automatic depths ran 1.24 times the manual depths of Blewett et al. (2016), 1.62 times
  poleward of 60 deg. For a bowl-shaped pit, h = L / tan(i) underestimates depth by up to about 25%. Raw EDRs show column striping.

## Result
| Class | n | Hollows | Depth |
|---|---|---|---|
| Depression with a shadow and a lit far wall | 11 | 4184, 2948, 2771, 2949, 5608, 4142, 3084, 2496, 865, 6230, 3970 | 4.1-11 m, median 6.7 m |
| Ambiguous (crater-like bowl, or small relief inside larger landforms) | 4 | 5606 (0-3.1), 597 (5.2-6.4), 601 (3.1-3.3), 1598 (0-5.1) | 0-6.4 m |
| No measurable shadow (bright patch, low mound, bright crater wall) | 6 | 2275, 2276, 2772, 3752, 4896, 5832 | < 1.0-1.9 m if depressions |

Depths per hollow (m): 4184 6.6; 2948 9.2-9.4; 2771 4.1; 2949 6.9; 5608 6.3-7.1; 4142 4.7-5.2 (frame edge); 3084 10.5-11.2 (three pits
in a bright patch); 2496 10.5 (widened box); 865 5.4-5.7; 6230 10.0; 3970 6.4. Ranges span the two box versions.

- None of the 21 is fully shadowed (shadow at most half the box along the sun). At these incidences a depression 24-28 m deep, the median
  hollow depth, would cast a shadow 39-163 m long and fill a box of this size; h_fill is 11-38 m.
- The 11 depressions are 4-11 m deep, a quarter to a half of the median hollow depth, with depth/diameter 0.05-0.2. Depth does not rise
  with diameter across 42-95 m (Spearman 0.15, n = 11).
- Depth/diameter of 0.05-0.2 overlaps that of small fresh and degraded impact craters, so morphology (irregular outline, flat floor,
  no raised rim, halo) has to decide which are hollows.
- Detection bias: the detector finds bright features. A small hollow at full depth would be mostly in shadow at these incidences and
  could be missed, so these frames cannot show whether full-depth hollows smaller than about 100 m exist.

## Step 3-5: depth against size at fixed pixel scale (9 Nov)
- depth_vs_size.py: the same automatic measurement on the finest detection of every non-trace hollow in frames with incidence 58-81 deg and
  emission <= 50 deg, in two tiers: <= 7 m/pixel (135 hollows, 40 m-1 km) and 7-15 m/pixel (641 hollows, 0.1-3 km). Detections from
  out_v2/detections_qc.csv (same hollow ids as v3). Output depth_vs_size.csv; zoom tiles of the <= 7 m tier in tiles_dvs/.
- dvs_analyze.py: classes (clean; upper limit where no shadow component exists, h < 1.5 px; lower limit where >= 85% of the box is dark;
  excluded for edge flags or shadows that no lit far wall bounds). Output dvs_summary.json, depth_vs_size_classed.csv.
- dvs_checks.py: subsets and a fit of the wall-support (Janssen) curve. Output dvs_checks.json.

Median automatic depth of clean measurements (m):
| D | <= 7 m/px | 7-15 m/px |
|---|---|---|
| 40-100 m | 5.3 (n 16) | - |
| 100-150 m | 10.0 (22) | 11.7 (11) |
| 150-200 m | 11.6 (29) | 15.7 (29) |
| 200-300 m | 15.7 (28) | 20.7 (107) |
| 300-500 m | 16.8 (17) | 25.1 (207) |
| 0.5-1 km | 34 (2) | 37.5 (108) |
| 1-3 km | - | 37.6 (26) |

- Censored fit (Tobit, log space, with ln px and ln tan i): ln depth rises 0.72 +- 0.07 per ln D. Clean-only OLS within each tier, and within
  single-depression boxes, gives 0.64-0.92. No step between 100 and 300 m. 13-35% of boxes show no shadow at all, falling with size.
- Wall-support curve with load-proportional compaction: transition width D_t = 1.9 km (1.2-3.3 km), better than a power law by 4 AIC.
  That would require a charged layer roughly 1-3.6 km thick, far thicker than the 60-280 m of S13, so it is not a good explanation as fitted.
- QC (qc_dvs_le7_sample.png, 12 random clean measurements at 100-500 m, <= 7 m/px): only 4 measure the main shadow of the depression
  (2497 9 m, 5441 21 m, 1933 36 m, 922 42 m; the last two look crater-like). In 7-8 the method measured a small shadow inside a complex
  box (knobby clusters, troughs) or missed the main shadow (2997). The automatic depths therefore underestimate many hollows, most in
  complex boxes, and the size trend is not reliable without manual measurement.
- Side observation: at fixed pixel scale, hollows poleward of 60N are smaller (median 0.39 vs 0.47 km at <= 20 m/px; 0.53 vs 0.99 km at
  20-40 m/px). If depth rises with size, part of the poleward decrease in depth could be a size effect. Not tested.

## Step 6: manual shadow measurements (Hollow Shadow Bench, 9 Nov)
- Artifact "Hollow Shadow Bench": https://claude.ai/artifact/FkT7M8UFS1ZeMka8h2nM8t (db collection "shadows", one doc per hollow id).
- Images and geometry from measure_tool_assets.py: tiles_meas/<id>.png (native EDR crop + 300 px context) and manifest_meas.json
  (crop origin X0, Y0 in frame pixels; box and automatic line in crop pixels; px, f, local incidence, sun azimuth). Bundle: tiles_meas_v1.tgz.
- Michael clicks the start (sunward rim) and tip of the main shadow. The page stores the two points in crop pixels (ax, ay, bx, by),
  L_px = |(B - A) . u| with u the light-travel direction (-cos az, -sin az), off_px (deviation across the sun), and
  depth_m = L_px x px x f / tan(i). Other fields: cls (hollow, crater, notdep, unsure), shadow (measured, full = lower limit,
  none = shallower than 1.5 px x px x f / tan i), note.

## Step 7: result of the hand measurements (9 Nov; bench_analyze.py -> bench_summary.json, bench_depths.csv)
- Michael's calls: 128 hollows, 5 unsure, 1 crater-like, 1 unclassed; 132 measured shadows, 1 lower limit (872), 2 with no shadow (3083, 251).
  18 measurements flagged uncertain from his notes. Two notes say depth varies within the hollow (1151, 5438).
- Median hand depth of measured hollows: 40-100 m 5.2 m (n 21); 100-150 m 11.1 (25); 150-200 m 11.8 (31); 200-300 m 15.9 (32);
  300-500 m 18.9 (19); 0.5-1 km 23.2 (3). Depth/diameter falls from 0.08-0.09 to 0.04-0.05.
- With pixel scale and incidence controlled, ln depth rises 0.57 +- 0.07 per ln D (Tobit with the limits 0.58 +- 0.08; certain-only 0.57 +- 0.08).
  A broken power law and a power law with a plateau both fit worse (AIC) than one power law: no step and no plateau between 40 m and 1 km.
- Hand vs automatic depths (114 with both): median ratio 1.01, 80% within 30%, Spearman 0.88. The automatic method was better than the
  12-tile QC suggested (Michael may have had the automatic line on for some).
- Automatic sample (602 clean, both tiers): without ln D, ln px takes 0.79 +- 0.06, the same as the Blewett et al. depth regression (0.78);
  with ln D it falls to 0.30 +- 0.07. Most of the pixel-scale dependence of measured depth is hollow size.
- Poleward of 60N in the automatic sample (127 of 602; median size 344 vs 346 m): depth ratio 0.90 (0.79-1.03) without ln D, 0.98 (0.87-1.10)
  with it. The pilot found automatic depths 1.62x manual poleward vs 1.24x overall, so this does not test the Blewett poleward decrease.
  The bench set has only a few hollows poleward of 60N.

## Step 8: common cause or size? (within_field.py -> within_field.json, 9 Nov)
- Michael: smaller hollows may form in less volatile-rich material (lower phi), which would make them shallower without size mattering.
- Test: slope of ln depth on ln D within groups that share a field (or host) and a frame (so pixel scale and incidence are shared),
  against the slope between group means.
| Data | Grouping | n (groups) | Within slope | Between slope | Pooled |
|---|---|---|---|---|---|
| Hand | field x frame | 94 (24) | 0.37 +- 0.10 | 0.74 +- 0.15 | 0.57 |
| Hand | host x frame | 80 (21) | 0.43 +- 0.09 | 0.78 +- 0.20 | |
| Automatic | field x frame | 422 (121) | 0.49 +- 0.06 | 0.84 +- 0.10 | 0.60 |
| Automatic | host x frame | 259 (67) | 0.49 +- 0.08 | 0.85 +- 0.12 | |
- Both effects are present. Fields of larger hollows are deeper than size alone predicts (between slope about twice the within slope),
  as a field-level factor such as phi would make them. Within a field, depth still rises with size at 0.4-0.5, about two-thirds to
  four-fifths of the pooled slope. That part is size itself or variation in phi below the scale of a field; these data cannot separate the two.
- For depth regressions on latitude or color: if phi sets both size and depth, controlling hollow size also removes part of the phi signal,
  so such checks should be run with and without hollow size.

## Step 9: does sunlight set hollow size? (size_sun.py -> size_sun.json, 9 Nov)
- Michael: hollows form with the crater and keep growing by sun-driven activity such as scarp retreat, so the model can explain smaller,
  shallower hollows poleward. Test on catalog hollows (finest detection <= 40 m/px), ln D_eq with ln px controlled.
- Hot longitudes: hollows are 1.07 (1.04-1.10) times larger per unit cos(2 lon), about 1.15 times from warm to hot longitude; 1.04 (1.00-1.07)
  within 45 deg of the equator without the province. Sunlight modulates size, as sun-driven growth predicts.
- Poleward of 60N: 0.85 (0.80-0.91) times the size, but 71% of poleward hollows are in the province; without it, 0.94 (0.83-1.07).
  The poleward size deficit is mostly the province, not a general effect of latitude.

## Steps 10-11: latitude and color depth regressions on our own depths (depth_all.py, depth_regress.py -> depth_regress.json, 9 Nov)
- Automatic depths for every suitable hollow: finest detection with px <= 20 m, incidence 58-81, emission <= 50 (2,012 hollows, 713 frames;
  1,558 clean). Same models as the paper (OLS of ln depth, SEs clustered by image), beside the 2,518 Blewett et al. (2016) depths.
  Color here is the host crater's interior color (crater_color.csv PC1/PC2), not the paper's 10 km surroundings; the Blewett depths are
  fitted with host color too as a bridge.
| Result | Blewett depths | Automatic depths | Automatic + ln D |
|---|---|---|---|
| Coefficient on 1/T | -291 +- 99 K | +105 +- 91 K | +89 +- 65 K |
| Poleward depth ratio | 0.77 (0.69-0.86) | 0.93 (0.83-1.03) | 0.98 (0.90-1.07) |
| Poleward, province removed | 0.72 (0.59-0.87) | 1.03 (0.83-1.28) | 1.05 (0.87-1.26) |
| Host PC1 (slope), per SD | 0.93 (0.88-0.99) | 1.02 (0.95-1.10) | 1.06 (1.01-1.12) |
- ln D slope 0.61 +- 0.03 in every automatic model. Poleward hollows in this sample: median 363 m vs 485 m; 72% are in the province
  (69% of Blewett's poleward depths are too).
- Checks: Blewett's poleward ratio stays 0.76 (0.68-0.85) in the same geometry (incidence 58-81, px <= 20), so the difference is not
  incidence or pixel scale. Michael's hand depths agree with the automatic depths poleward (median ratio 0.98, n = 16) as well as
  equatorward (1.02, n = 99).
- So our depths do not reproduce the temperature, poleward, or color dependence of the Blewett depths, and hollow size accounts for most
  of the small poleward decrease they do show. The difference lies in which hollows were measured or how. A direct test: hand-measure on
  the bench the hollows Blewett measured in the pilot's registered frames (validation.csv), poleward and equatorward.
## Step 12: the paper's depth results on our depths (depth_paper.py -> depth_paper.json, depth_paper_rows.csv; fig_depth_v2.py, 9 Nov)
- depth_paper.py reruns every depth regression of the paper (temperature, poleward, hot longitudes, host diameter, nesting, Herrick and Kinczyk states, vents, LRM, XRS, host color, between-host ANOVA) on the 1,558 clean automatic depths, without and with ln hollow diameter. Kinczyk merge deduplicated (202 unmatched rows without an id); merges validated 1:1 / m:1.
- Main results (ratios, 95%; without / with ln D): 324 K to 219 K 0.86 (0.66-1.12) / 0.88 (0.73-1.06); poleward 0.93 / 0.98 (0.90-1.07); per unit cos 2 lon 0.94 / 0.96; host diameter per tenfold 1.16 / 1.13 (1.02-1.25); nested 0.94 / 1.01; fresh vs degraded 1.18 / 1.01; Kinczyk 3 and 4-5 vs 1-2: 1.21, 1.30 / 1.03, 0.99 (Class 1-2 hollows smaller: median 340 vs 529 and 625 m); vents 1.14 / 1.16; LRM >= half 0.99 / 0.95; S/Si 0.93 / 1.02; host PC1 per SD 1.02 / 1.06. ln D slope 0.59-0.62 in every model.
- fig_depth_v2.py writes main-text Figure 5 of v6.11 (copied to manuscripts/.../figures/fig5_depth_v2.png). 9 Nov evening: panel (a) legend moved off the sulfur-front line, which is labeled directly.
## Step 13: supplement figure of automatic shadows (figS_shadow_examples.py -> figS_shadow_examples_v1.png/.pdf/.json, 9 Nov)
- Michael asked for examples of automatically detected shadows in S9 (now Figure S7 of supplement v23).
- (a) Six hand-measured hollows (<= 7 m/px), one per log-size bin 40 m-1 km, each the hollow whose |ln(hand/auto)| is closest to the set median (0.107); uncertain measurements and edge-flagged automatic runs excluded. (b) Six accepted automatic depths at 10-20 m/px drawn at random (seed 7) within 150-300, 300-500 (2), 500-800 (2), 800-1,500 m.
- Recomputed depths match depth_all.csv exactly. In three of the six panels in (a) the automatic and hand lines fall on different shadows of a multi-hollow box, though the depths agree within 9-11%.
## Step 14: depth against diameter for S9 (figS_depth_size.py -> figS_depth_size_v1.png/.pdf, 9 Nov)
- Michael's v22 comment: scatter plot in "Depth and hollow size" (Figure S8 of supplement v24). (a) hand depths corrected to 5 m/px and 70 deg with the hand fit, limits marked; (b) the 1,558 automatic depths corrected to 10 m/px and 70 deg, the 320 no-shadow boxes as upper limits, bin medians, least squares (0.61) and the diameter-scatter-corrected exponents from g26 (0.73, 0.90).

## depth_paper.py keys (9 Nov 2026)
fresh_vs_degraded(+lnD), kinczyk_n, kinczyk_hosts, kinczyk_median_D_m and kinczyk_median_depth_m were first added to depth_paper.json
outside the script. depth_paper.py now computes them (fresh + abslat on fresh and degraded hosts; counts and medians by Kinczyk class,
numpy rounding); a rerun reproduces every value of the stored record (depth_paper_rows.csv unchanged, md5 e70de884).
