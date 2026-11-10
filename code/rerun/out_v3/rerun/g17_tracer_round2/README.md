# g17_tracer_round2: second tracing round (7 Nov 2026)

Michael (7 Nov): trace the next 50 hosts and the peaks of degraded craters, so the terrain results and the knob comparison rest on hand labels.
The Hollow Host Terrain Tracer (artifact Hnvf7kvMZBqPpkaKvRX6So, versions 5-6) now lists four groups; round-1 traces in the "hosts"
collection are untouched (64 documents checked after the publish). New documents carry group and task fields.

## Groups (select_craters.py -> selection.json)
- A, hosts round 1: the 64 craters of g12 (62 done, 2 can't tell).
- B, hosts round 2: the 50 untraced terrain-test hosts with the most hollows (383 hollows, at least 5 each; the 12 hosts with 5 hollows all
  included). With A, the traced hosts would hold 3,328 of 3,760 terrain-test hollows (88.5%). Full task; outlines start from curvature v3, as in round 1.
- C, degraded craters with hollows: the 84 Kinczyk et al. (2020) Class 1-2 craters with assigned hollows that are not in A or B
  (143 hollows; 14 more Class 1-2 hosts are in B and 10 in A). Peaks only; outlines start empty (curvature v3 shown dashed, copy on request).
- D, degraded sample: 60 Class 1-2 craters without hollows, PPS on the completeness-weighted area within 0.15 R (5 taken with certainty;
  inclusion probabilities in selection.json and data.json "pi"). Peaks only; outlines start empty.
Knob comparison once traced: hollows on traced peaks of A, B, C Class 1-2 craters / (traced peak area of those craters + Horvitz-Thompson
estimate of the peak area of all 2,052 Class 1-2 craters without hollows from D). With curvature areas as a stand-in, the relative standard
error of the D total is 0.23.

## Files
- select_craters.py, selection.json; build_tracer_assets.py and build_color_assets.py (g12 copies; craters from selection.json) -> assets/;
  make_packs.py -> pk/<uid>.json (img, dem, col, ecol as data URLs; one published file per crater keeps the artifact within its file limit);
  round2_bundle.tgz (selection, meta, packs; staged to the cloud); build_data_r2.py (cloud) -> data.json; tracer_page_r2.html (published page).
