# g19_knobs_traced: knobs against the peaks of Class 1-2 craters traced by hand (7 Nov 2026)
- traces/traces_raw_artifact_2026-11-07c.json: tracer database after Michael finished groups C and D. clean_traces.py (g18 copy, all groups) ->
  traces/traces_clean_r3k.json. traced_classes.py (g18 copy) reads it.
- knobs/class12_peaks_traced.py (g13 copy; 2,160 Class 1-2 craters 30-250 km; chunks 0-1200, 1200-2200) -> knobs/_chunks_r3k.
- knobs/knob_ratio_design.py (new) -> class12_peaks_traced_r3k.csv, knob_ratio_design.json: hosts exact (108; 7 can't tell count as no peak),
  non-hosts by Horvitz-Thompson from the g17 PPS sample (60; 1 can't tell dropped, non-certainty part rescaled); 4,000 bootstrap replicates.
- Result: knobs / degraded peaks 0.76 (0.32-1.63); can't-tell craters with curvature peaks 0.81 (0.35-1.81); curvature everywhere 2.54 (1.24-4.90).
  Design check: HT of curvature area from the sample 38,800 km2 against the full curvature total 46,700 km2.

- 201427 (group C, Class 1, 3 hollows): the tracer note reads "no peak obvious. using automated", but the crater is status done with no shapes, so every estimate except curvature_all counts it as a crater with no peak (0 peak hollows, 0 km2). Michael confirmed on 7 Nov 2026 that it has no real peak and that the central topography is a different feature, so this treatment stands and no number changes. It is not in the 111-host terrain test (fewer than 5 hollows).
- knobs/knob_checks_nocurv.py (8 Nov 2026) -> knob_checks_nocurv.json: checks with no curvature classes. The 7 can't-tell hosts hold 13 hollows (1 within 0.15 R); counting hollows and area within 0.15 R of the 8 can't-tell craters as peak gives 0.79 (0.34-1.76). One-component bootstrap intervals: knobs 0.42-1.18, hosts 0.44-1.46, non-host sample 0.60-0.96. Supplement v18 reports the 13 hollows in place of the curvature sensitivity.
