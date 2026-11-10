# g18_traced_r3: terrain and positions from the hosts traced in rounds 1 and 2 (7 Nov 2026)

Michael finished round 2 of the tracer for hosts (group B: 49 done, 1 can't tell). The 114 hosts traced in rounds 1 and 2 are all host craters
30-250 km with at least 5 hollows; 111 are done and 3 can't tell. They hold 3,323 of the 3,761 terrain-test hollows (88%).

## Steps
- traces/traces_raw_artifact_2026-11-07b.json: all 203 tracer documents (ArtifactData read); traces/tracer_data_r2.json: the tracer's data.json.
- clean_traces.py (g13 copy; groups A and B only; a round-2 host counts once its wall-foot check is done) -> traces/traces_clean_r3.json.
- traced_classes.py, terrain/peak_wall_traced.py, merge_rows.py, analyze.py, terrain_table.py, terrain_obs_exp.py, positions/positions_build.py,
  positions_slope_build.py, merge_chunks.py: g13 copies reading traces_clean_r3.json (G13_REV r3; slump counts as wall, as in r2). Traced hosts use
  the tracings, others curvature v3 (rows_traced_r3_v3.csv); rows_tracedonly_r3_v3.csv keeps the 111 traced hosts, and every reported number uses it.
- positions/positions_analyze.py, positions_slope.py, collect_results.py: g16 copies on the g18 caches, traced hosts only -> positions/out_tracedonly_r3,
  positions/results_tracedonly_r3.json (side by side with the g16 62-host values).
- terrain/make_figs_pr_v9_v3cat_traced_r3.py -> fig1_exposure_v9_v3cat_traced_r3.png (Figure 2, panel d from the 111-host table).

## Results (62 hosts, g16 -> 111 hosts)
- Peaks/walls 6.36 (3.50-12.16) -> 7.17 (4.33-12.17); fields 3.05 -> 3.46 (2.53-4.93); floors/walls 1.52 (0.85-2.84) -> 1.63 (0.98-2.73);
  equator-/pole-facing walls 2.12 (1.20-3.92) -> 1.80 (1.08-3.11); without top 3 5.6 (3.6-9.5); Mantel-Haenszel by slope 6.5; isolated 14.1 (7.4-28.4).
- Table S2 (peak/floor/wall per 1,000 km2, all): 29.8/6.8/4.2; observed/expected 3.66 / 0.81 / 0.46; all-terrain density 8.6.
- Peak within 2 km of floor 0.51 -> 0.48 (0.39-0.59). Highest third: hollows 0.99 (0.69-1.48); fields 1.26 -> 1.51 (1.08-2.07) (eta_b 1.83; with slope
  controlled 1.37, 0.92-2.01).
- Floor within 0.2 R of center 0.68 -> 0.73 (0.42-1.09). Within distance bins: margin 1.04 (0.85-1.27), wall 0.87 (0.69-1.08), peak 1.18 (1.00-1.39),
  fields beside the peak 1.50 (1.16-1.96).
- Curvature v3 classes agree with the tracings for 70% of 3,288 hollows in the traced hosts.
- Not yet: knobs vs Class 1-2 peaks (group D sample 5 of 60 traced).

## Figure 2 with the v3r Caloris record (9 Nov 2026)
`CALORIS_RECORD=.../g5_figures/_scratch/caloris_v3r/caloris_test_record_v9.json FIG1_NAME=fig1_exposure_v9_v3cat_traced_r3_v3r python3
terrain/make_figs_pr_v9_v3cat_traced_r3.py` wrote terrain/fig1_exposure_v9_v3cat_traced_r3_v3r.png/.pdf; only panel c differs from the
previous figure (figure_numbers_*_v3r.json). The four Figure S3 craters are identical in the round-2 and round-3 traces.
