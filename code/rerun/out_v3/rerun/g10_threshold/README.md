# g10_threshold: is the 250 m peak/wall relief threshold defensible? (6 Nov 2026)

No README, script, or note records why 250 m was chosen (peak_wall_test/README.md states the rule only).

## relief_survey.py -> relief_survey.csv
2,055 Herrick craters 30-250 km (1,220 with cataloged central peaks, cp/mp; 694 filled craters with no cataloged central
structure; 141 peak-ring craters). Same floor level as the terrain classes. Medians by diameter bin:
- Central relief (max height above the floor within 0.3 R): craters with cataloged peaks 246, 303, 429, 547 m at 30-45,
  45-60, 60-90, 90-130 km; filled craters without cataloged peaks 195, 239, 288, 332 m. A fit to the peak craters gives
  central relief ~ D^0.80 (204 m at 30 km, 742 m at 150 km) and depth ~ D^0.74.
- So a fixed 250 m is about the whole median peak at 30-45 km and less than half of it at 90-130 km, and at 30-60 km the
  665 m topography separates cataloged peaks from peakless floors by only 50-65 m.

## peak_wall_scaled.py (copy of g7 peak_wall_fix.py) -> out_v3_<mode>/, threshold_comparison.csv
v3 catalog, v3 grid, gap-free classes. Thresholds: fixed 150/250/400 m (from g7), 0.10/0.15/0.20 x crater depth (>= 100 m),
and hybrid max(250 m, 0.15 x depth).
- Peaks to walls: 4.98 (150 m), 5.72 (250 m), 7.91 (400 m); 5.29, 6.29, 6.43 (0.10, 0.15, 0.20 x depth); 6.90 (4.37-10.38) hybrid.
- Floors to walls: 2.53-3.25 for every threshold. Equator- to pole-facing walls: 1.50-1.74, lower bounds 1.01-1.12.
- By size: 30-60 km 2.2-3.5 with wide intervals (lower bound 0.44-2.07); 60-120 km 7.0-11.9; 120-250 km 5.5-9.1.

## class12_peaks.py -> class12_peaks.csv, class12_knob_ratio.json
Knob density against the peak terrain of Kinczyk Class 1-2 craters 30-250 km (2,160 craters, hosts and non-hosts), with the
terrain classes in place of the 0.15 R proxy: 3.24 (1.59-6.81) fixed 250 m; 3.50 (1.66-7.71) 0.15 x depth; 3.07 (1.44-6.77)
hybrid; the old proxy gives 2.37 (1.01-7.35) (g9: 2.35, 1.005-7.24).
