# g13_traced: terrain analyses with Michael's traced outlines (6 Nov 2026)

Michael traced the 64 host craters that hold 80% of terrain-test hollows in the Hollow Host Terrain Tracer
(artifact Hnvf7kvMZBqPpkaKvRX6So; page source and assets in g12_tracer/). He marked 62 done and 2 "can't tell" (110709, Namatjira).
He drew his own peak outlines in most craters, kept the curvature v3 wall foot except in Sholem Aleichem and Eminescu (moved inward deliberately, confirmed by Michael,
over wall slump, which he counts as wall), and kept v3 outlines in a few craters (notes: "auto did good enough", Durer, Aksakov, Wang Meng).

## Inputs and rules
- traces/traces_raw_artifact_2026-11-06.json: the 64 database documents as read. traces/traces_clean.json: per crater, peak geometry
  (peak outlines minus floor-inside-ring outlines) and wall foot, as lon/lat rings. Scarlatti's ring was drawn as one self-intersecting
  loop (outer then inner); it is read as an annulus (the inner loop is the hole; a sliver hole < 5% of the outline area at the joining
  stroke is filled), as Michael's note asked.
- traced_classes.py: peak = inside peak outlines; wall = outside the wall foot out to 1.05 R, not peak; floor = inside the wall foot,
  not peak. Hosts not traced (259, including the 2 "can't tell") use curvature v3 (g10 curv_classes_v3.py).
- Everything else as g11: v3 catalog, v3 completeness grid (Dref 1 km), 321 hosts 30-250 km, 3,761 hollows, superposed-crater mask.

## Scripts (copies carry their source and changes on line 1)
- terrain/peak_wall_traced.py (copy of g11 peak_wall_curv.py; env CLS v3|v2, NOTRACE=1 for curvature alone) -> terrain/out_v3/rows_<cfg>.csv
  via merge_rows.py; also per-hollow classes under the primary rule, v3 and v2 (_chunks/hollow_classes_*) and traced/curvature
  area overlaps (_chunks/area_overlap_*). analyze.py, terrain_table.py, terrain_obs_exp.py are g11 copies with paths changed.
  Check: NOTRACE=1 CLS=v2 reproduces g11 rows_curv_n24_h100_s0.03.csv exactly (33 area and count columns, max difference 0).
- terrain/make_figs_pr_v9_v3cat_traced.py -> fig1_exposure_v9_v3cat_traced.png/.pdf (Figure 2; panel d from the traced table).
- terrain/figS2_terrain_examples_traced.py -> figS2_terrain_examples_traced.png/.pdf (Figure S2 on the MDIS 166 m mosaic, traced classes).
- positions/ (g11 copies, G11_RULES=traced) -> positions/out_traced/, positions/results.json (collect_results.py).
- knobs/class12_peaks_traced.py, class12_knob_ratio_traced.py -> knobs/class12_knob_ratio_traced.json (10 of 2,160 Class 1-2 craters traced).
- build_summary.py -> results.json and changes_traced.csv (129 rows, g11 curvature v2 vs g13; conclusion change = interval crosses 1,
  presence test crosses p = 0.05, or class order changes).

## Main results (g11 curvature v2 -> g13 traced + curvature v3)
- Hollows per class peak/floor/wall/other/masked: 1,588/1,607/511/54/1 -> 1,493/1,613/600/54/1. Weighted area (1,000 km2) 121/295/245 -> 86/326/248.
- Peaks/walls 6.26 (4.05-9.50) -> 7.19 (4.00-12.80); fields 2.63 -> 3.44 (2.67-4.50); floors/walls 2.61 -> 2.05 (1.29-3.37);
  inner floor 3.59 -> 2.44; outer floor 2.23 -> 1.88; equator-/pole-facing walls 1.22 (0.81-1.82) -> 1.57 (0.83-2.74), still includes 1.
- Groups, peaks/walls: ring 14.9 -> 28.2; nested 3.01 -> 2.06 (1.07-4.52); unclassified 2.07 -> 1.28 and filled 1.58 -> 0.76 (both include 1, as before).
- Table S2 densities (peak/floor/wall per 1,000 km2): all 13.1/5.5/2.1 -> 17.4/4.9/2.4; nested 6.7/12.7/4.9 -> 10.2/8.5/7.8, so the nested
  order becomes peak > floor > wall (was floor > peak > wall). Observed/expected peak 2.10 -> 3.45, floor 1.02 -> 0.88 (0.71-1.05), wall 0.37 -> 0.42.
- Peak height (positions): individual hollows no longer favor the top of peaks. Lowest third 0.69 (0.47-0.91) -> 1.09 (0.78-1.37); highest third
  1.50 (1.08-2.23) -> 1.06 (0.73-1.59); same with eta_b and with the slope control. Fields still do: highest third 2.17 -> 1.89 (1.45-2.47),
  lowest third 0.55 -> 0.64 (0.45-0.88). The deficit at the peak base strengthens: peak within 2 km of floor 0.75 -> 0.49 (0.41-0.60).
- Floor positions: within 0.2 R of the center 1.91 -> 1.58 (1.08-2.22); within 0.05 R of the peak (radius-stratified) 1.42 (1.12-1.74) ->
  1.16 (0.94-1.36), now includes 1; within 0.1 R of the wall (radius-stratified) 0.78 (0.65-0.93) -> 0.82 (0.71-0.96).
- Knob density / Class 1-2 peak density 2.88 (1.45-5.56) -> 3.24 (1.58-6.61).
- Decomposition: curvature v3 alone gives peaks/walls 6.44, floors/walls 2.79, equator/pole 1.23, so the changes come from the tracing.
  Traced + curvature v2 elsewhere: 6.76 / 2.06 / 1.56. Without the 13 traced hosts whose notes flag poor images, uncertain peaks or
  (Abu Nuwas) poor hollow detection: 8.24 (4.31-16.09) / 2.24 (1.34-4.11) / 1.63 (0.78-3.20).

## Curvature classes against the tracings (2,945 hollows in the 62 traced hosts)
- Hollow classes agree for 75% with curvature v3 (73% with v2). 345 hollows on traced peaks are v3 floor; 241 hollows on traced floor are
  v3 peak; 113 hollows on traced wall are v3 floor (the slump Michael added to the wall in Sholem Aleichem and Eminescu).
- Peak area: traced 39,700 km2, v3 55,200 km2, overlap 30,900 km2 (IoU 0.48); 78% of traced peak lies in v3 peak, 56% of v3 peak in traced peak.

## Michael's notes worth carrying forward
- Abu Nuwas: most detections look like small craters and the visible hollows are missed (catalog check).
- Dario: hollows on the infill and in small craters, possibly venting through volcanic infill.
- Sholem Aleichem: no peak complex; all hollows on a slumped section of wall.

## Revision r2 (7 Nov 2026): Michael's refined wall feet and slump
traces/traces_clean_r2.json (clean_traces.py from the 7 Nov database snapshot, traces_raw_artifact_2026-11-07.json): all 64 hosts have
the wall-foot check done; feet drawn or edited in 56; 13 slump outlines in 11 hosts; peaks as before (62 done, 2 can't tell).
Rule in r2: slump counts as wall (slump inside the wall foot becomes wall; peak still wins). Env G13_REV=r2 is the default in every g13
script; G13_REV=r1 reproduces the 6 Nov outputs (checked: maximum difference 9e-13 in the rows table). r2 outputs: terrain/out_v3/*traced_r2_v3*,
terrain/fig1_exposure_v9_v3cat_traced_r2.png, terrain/figS2_terrain_examples_traced_r2.png, positions/out_traced_r2, positions/results_r2.json,
knobs/*_r2*.
- Hollows peak/floor/wall/other/masked: 1,493/1,613/600/54/1 -> 1,499/1,484/723/54/1 (123 floor hollows become wall).
- Peaks/walls 7.19 -> 6.44 (3.72-10.89); fields 3.44 -> 3.37; floors/walls 2.05 -> 1.77 (1.14-2.81); inner 2.44 -> 2.12; outer 1.88 -> 1.62 (1.02-2.60).
- Equator-/pole-facing walls 1.57 (0.83-2.74) -> 1.72 (1.03-2.84): the interval now excludes 1 (equator-/pole-facing wall hollows 465/258).
- Ring craters peaks/walls 28.2 -> 18.5; nested 2.06 -> 1.94 (1.04-3.85); central-peak 8.02 -> 6.59.
- Table S2 all 17.5/4.8/2.7 (peak/floor/wall); nested 10.0/7.8/8.6 (order peak > wall > floor); observed/expected peak 3.47, floor 0.91 (0.71-1.12), wall 0.44.
- Positions essentially unchanged (peak heights: hollows 1.07 lowest / 1.08 highest third; fields 0.66 / 1.91). Floor within 0.1 R of the wall
  (radius-stratified) 0.82 -> 0.86 (0.71-1.03), now includes 1; floor within 0.05 R of the peak 1.15 (0.94-1.35) unchanged.
- Knobs / Class 1-2 peaks 3.24 (1.58-6.61), unchanged.
- slump_descriptive.py -> slump_descriptive.json: in the 11 hosts with slump, hollows per 1,000 km2 (completeness weighted): slump 12.9 (158
  hollows), rest of wall 1.1 (18), floor 12.0 (453), peak 53.6 (644). 111 of the 158 slump hollows are in Sholem Aleichem and 31 in Sander;
  7 of the 11 hosts have 0-3 hollows on slump. Slump was not mapped blind to hollows. Descriptive only.
