# g22_province_screen: screening the hollows of the province excess (9 Nov 2026)

Michael wants to know whether the excess of host craters near 65N 125E (region 1 of g20) comes from false detections. 17 of the region's
23 hosts qualify only by hollows < 1 km, and the v3 screen covered detections >= 1 km.

- select_province.py -> province_screen_set.csv (81 hollows within 0.9 R of the 23 hosts; representative detection and bbox from the
  0.8 release) and province_hosts.csv.
- tiles_province.py -> tiles/<hollow_id>_z.png (box + 1.2 box widths) and _c.png (~8 km), raw EDR of the representative detection;
  EDRs fetched from the PDS Imaging Node and deleted after tiling. province_tiles_v1.tgz bundles the tiles and CSVs.
- Screening page: artifact "Province Hollow Screen" (https://claude.ai/artifact/LqEcq5GKsjPBWqT6wgaucM), db collection "calls",
  one document per hollow_id: call (hollow | not | unsure), kind (crater, shadow, slope, bright, other), note, updated.
- Next: read the calls back, recompute hosting for the 23 craters, and rerun g20 clusters.py with the false detections removed.

## Result (9 Nov 2026)
- Michael's screen: all 81 are bona fide hollows except 1802 and 1810, marked unsure (calls_2026-11-09.json). Each unsure hollow is in a
  host with two other confirmed hollows (1810 in 200761, 102 km; 1802 in Kofi, 200532, 128 km), so all 23 craters remain hosts and the
  g20 excess map is unchanged (23 hosts within 500 km of 65N 125E, 4.9 predicted).
- Imaging check (g20_regional_clusters/province_fine_imaging_check.json): among craters 5-30 km at 45-75N whose best image is <= 20 m/pixel,
  13 of 143 in the province host hollows (9.1%), against 9 of 831 elsewhere (1.1%). Fine imaging does not explain the excess.
