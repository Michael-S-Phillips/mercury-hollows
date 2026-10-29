#!/bin/bash
# Rerun every analysis and figure that needs only the data in this repository, in dependency order.
# Usage: bash run_all.sh            (then: git status / git diff --stat to see what changed)
# Steps that need external data (ISIS, MDIS EDRs and mosaics, the USGS DEM, the PDS image index, the Bickel et al. release)
# are listed in README.md and are not run here.
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"; export MERCURY_HOLLOWS_ROOT="$ROOT" MPLBACKEND=Agg
G="$ROOT/manuscripts/GRL_insolation_hollows"; P="$ROOT/manuscripts/PlanetaryResearch_hollows"
fail=0
run() { local dir=$1; shift; printf '%-75s' "$(basename "$dir")/$1"; if (cd "$dir" && python3 "$@" > /dev/null 2>&1); then echo ok; else echo FAILED; fail=1; fi; }
# hosting, plains, insolation
run "$G/scripts" test2_h18.py
run "$G/scripts" latcorr.py
run "$G/scripts" latcorr3.py
run "$G/scripts" latcorr5.py
run "$P/analysis_v3" dT_hosting.py
run "$P/analysis_v3" dT_plains.py
run "$P/analysis_v3" transition_layer.py "$G/scripts" "$P/analysis_v3/craters_dT.csv" "$P/analysis_v3"
run "$P/analysis_v3" transition_layer_hilat.py "$G/scripts" "$G/data/test2_craters_h18.csv" "$P/analysis_v3/transition_layer_grid.npz"
run "$ROOT/code/rerun/latitude/hotwarm" hotwarm.py
run "$ROOT/code/rerun/latitude/aspect" aspect_stats.py
# terrain, Caloris, host age, positions
run "$ROOT/peak_wall_test" analyze.py rows_H150_r0.45.csv rows_H250_r0.45.csv rows_H400_r0.45.csv rows_H250_r0.35.csv
run "$ROOT/peak_wall_test" caloris_test.py
run "$P/analysis_v3" terrain_table.py
run "$P/analysis_v3" host_age_tests.py
run "$P/analysis_v3" pair_offsets.py
# depth
run "$G/scripts" test1.py
run "$G/depthctl" depth_structure.py
run "$G/depthctl" depth_structure_extras.py
run "$G/depthctl" depth_vs_xrs.py
run "$G/depthctl" depth_vs_xrs_tiles.py
run "$G/depthctl" strength_calc.py
run "$ROOT/slope_test" occ_analysis.py
run "$ROOT/slope_test" depth_slope_geometry.py
run "$ROOT/slope_test" depth_slope_eqface.py
# sulfur charging
run "$P/analysis_v3" sulfur_charge.py
run "$P/analysis_v3" sulfide_s2.py
# figures (writes figures/*.png|pdf and figures/data/*_data.csv; make_figs_pr.py names panel set 2 'fig2_insolation')
run "$P/figscripts" make_figs_pr.py
run "$P/figscripts" slope_fig.py
run "$P/figscripts" lifecycle_fig_v5.py "$G/scripts" "$P/figures/fig5_lifecycle_v5.png"
exit $fail
