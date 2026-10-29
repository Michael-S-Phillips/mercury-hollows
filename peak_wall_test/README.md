# Peak vs wall vs floor hollow density inside hollow-hosting craters

Question: are hollows richer on central peaks than on crater walls per exposed area? If not, exposure alone could explain hollows on peaks.

Method (peak_wall.py, analyze.py): 352 hollow-hosting Herrick (2018) craters 30-250 km across. Terrain classes from the USGS 665 m DEM:
peak = within 0.45 crater radii (0.6 for peak-ring classes) and >250 m above the floor (20th percentile height at 0.15-0.6 radii);
wall = 0.55-1.05 radii (0.7-1.05 for rings) and >250 m above the floor; floor = within 0.8 radii and <=250 m above it.
Areas are weighted by catalog completeness for 1 km hollows (10 ppd grid). Superposed cataloged craters >=5 km are masked,
because hollows in them are assigned to them. Hollows counted only if assigned to the host (smallest containing crater, 1.1 R).
"Nested" = center inside a cataloged crater at least twice as large (68 craters), or inside an Orgel et al. certain/probable basin (73 more).

Files: rows_H*_r*.csv (per crater; H = relief threshold, r = peak radius), peak_wall_record.json (all statistics and splits),
peak_wall_by_setting.json and .png (densities by terrain for isolated vs nested craters), analyze_log.txt.

## Nested-crater penetration test (nested_build.py, nested_test2.py, nested_test3.py, caloris_test.py, nested_figs.py)
Question: inside a larger crater, do only children that cut through the parent's floor fill host hollows?
- morph_all.csv: DEM depth and terrain areas for all cataloged craters >= 10 km (measure_all.py).
- nested_children*.csv: children >= 10 km whose center lies inside a crater >= 2.5x larger; parent fill from depth deficit against the
  90th-percentile depth of isolated craters (fresh_envelope.json); child excavation depth 0.1 x transient diameter (Croft 1985, Dt = 12 km).
- Parents 60-300 km: inconclusive (fill estimates rough, few hollow fields): fill-only 0/253 children host a field of >=3 hollows, through-fill 6/561.
- Caloris (plains 2.5-3.5 km thick, Ernst et al. 2015): children excavating < 2.5 km: 0/128 host a field; >= 2.5 km: 14/32 (Fisher p = 1e-11);
  with completeness >= 0.4: 0/13 vs 9/12. caloris_test_record.json; figures fig_nested_A/B/C.
