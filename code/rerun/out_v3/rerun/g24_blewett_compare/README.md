# g24: blind check of Blewett et al. (2016) depths against Michael's hand depths (9 Nov 2026)

Why: our automatic depths (validated against Michael's hand depths, g23) show no temperature or color dependence and a poleward depth
ratio of 0.93 (0.98 with hollow size), while the Blewett depths give -291 K, 0.90 per SD of color, and 0.77 poleward, also in the same
incidence and pixel-scale range. This compares the two measurements on the same hollows in the same images.

- cmp_assets.py: from the depth pilot's registration (accepted frames, Blewett points inside Bickel boxes), every poleward box (28, 14 frames)
  and 72 random equatorward boxes (at most 3 per frame, seed 24). tiles_cmp/<fid>.png, manifest_cmp.json (page data: no Blewett depths,
  no latitude), key_cmp.csv (Blewett max and median depth per box, latitude, the pilot's automatic depth). Bundle cmp_v1.tgz.
- Artifact "Blewett Depth Check": https://claude.ai/artifact/ErgMkKroiHNaYNjWYuK3dk, db collection "blewett_cmp" (doc id = fid; fields as
  in the g23 bench: ax, ay, bx, by in crop pixels, L_px, depth_m, cls, shadow, note). Michael measures the longest shadow in each box,
  blind to Blewett's values; magenta crosses mark the registered Blewett points.
- Planned comparison: Michael's depth against Blewett's deepest point in each box, poleward and equatorward (ratio and its 95% interval
  by bootstrap over frames); also the implied shadow lengths, and the pilot's automatic depth.
- Key medians before measuring: Blewett max 24.4 m equatorward, 17.5 m poleward; pilot automatic 26.0 and 21.4 m.
- Version 2 of the page (9 Nov): one line per Blewett point. Boxes hold 178 points (42 boxes have more than one, up to 7). Doc fields
  now: lines = {k: {ax, ay, bx, by, L_px, depth_m, off_px} or {skip: true}}, with k the 1-based point number in manifest_cmp.json
  "marks"; depth_m = deepest line; n_lines, n_points. marks_key.csv maps (fid, k) to the Blewett point id and depth. Two docs saved under
  version 1 (1118, 9226) carry a single line in ax..by; the page assigns it to the nearest point.

## Result (cmp_analyze.py -> cmp_summary.json, cmp_points.csv; blewett_registered_size.json; 9 Nov)
- Michael measured all 178 points in the 100 boxes, blind. Same features: Spearman 0.86-0.89 between his depths and Blewett's.
- Blewett depths are shallower than Michael's by a similar factor everywhere (Blewett/Michael, median, 95% by frame bootstrap):
  points 0.83 (0.77-0.88); equatorward 0.85 (0.78-0.91), poleward 0.81 (0.70-0.92). Shadow-length ratios are the same, and the
  incidence angles agree (median difference 0.03 deg), so the offset is where the shadow edges were placed. Poleward factor on the
  ratio: 0.95 (0.82-1.10), p = 0.48: no poleward-specific difference in measurement.
- On these same points, the poleward depth ratio is 0.84 (0.70-1.02) from Blewett's depths and 0.88 (0.74-1.05) from Michael's.
- All 817 registered Blewett points (59 poleward): poleward ratio 0.89 (0.76-1.04); with ln box size added 0.95 (0.82-1.10). Depth rises
  with box size in Blewett's own depths (0.34 +- 0.04 per ln D), and the ln px coefficient falls from 0.85 to 0.55 with size.
- Pilot automatic / Michael on these boxes: 0.96 equatorward, 1.09 poleward.
- Reading: the measurement method does not make the poleward decrease. The decrease is 0.77 in the full Blewett set, 0.84-0.89 in the
  registered subsets and 0.93 in our automatic sample, and 0.95-0.98 wherever hollow size is controlled. Blewett-style depths run about
  15% shallower than Michael's (and the automatic method's), uniformly.
