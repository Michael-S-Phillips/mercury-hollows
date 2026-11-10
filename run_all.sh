#!/bin/bash
# Reruns the analyses of the hollow formation paper that need only the data in this repository (README, "Reproduction check").
# Run setup_root.sh first. Records are rewritten in place; compare with `git status` and `git diff --stat`, and restore with `git checkout -- .`.
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"; export MERCURY_HOLLOWS_ROOT="$ROOT" MPLBACKEND=Agg
if [ "$(cd "$HOME/mnt/hollows-emissivity" 2>/dev/null && pwd -P)" != "$(cd "$ROOT" && pwd -P)" ]; then echo "Run setup_root.sh first."; exit 1; fi
R="$ROOT/code/rerun/out_v3/rerun"; P="$ROOT/manuscripts/PlanetaryResearch_hollows"
fail=0
run() { local dir=$1; shift; printf '%-80s' "${dir#$ROOT/}: $*"; if (cd "$dir" && bash -c "$*" > /dev/null 2>&1); then echo ok; else echo FAILED; fail=1; fi; }
run "$ROOT/code/rerun/out_v3" "python3 build_v3_catalog.py"
run "$R/g4_completeness/g1" "bash run_dag.sh v3r main"
run "$R/g4_completeness/g1" "bash run_degr_vents.sh v3r"
G2="dT_hosting dT_hosting_completeness dT_plains transition_layer_hilat latcorr latcorr3 latcorr5 aspect_stats fractions_g2 hollow_rich_hosts band_crater_hotwarm tm_vs_dT_check cos2lon_lat15 size_quirk hotwarm threshold wall_fractions plains_hosts make_figs_pr_insolation_v9"
[ -f "$ROOT/data/bickel-2024/Bickel_et_al_2024_MDIS_hollows_08_v1.0.csv" ] && G2="$G2 stage_test stage_test_lowinc stage_fieldlevel"
run "$R/g4_completeness/g2" "G2_CAT=v3r bash run_g2.sh $G2; ! grep -v 'exit=0\|done' out/v3r/run_status.txt"
run "$R/g3_catalog/sfd" "bash run_sfd.sh v3 refit"
run "$R/g3_catalog/sfd" "python3 sfd_derived.py"
run "$R/g3_catalog/caloris" "G3_CAT=v3 G3_GRID=v3r python3 caloris_v9.py"
run "$R/g3_catalog/caloris" "G3_CAT=v3 G3_GRID=v3r python3 caloris_v9_controls.py"
run "$R/g6_remaining/color" "python3 depth_color_fits.py v3"
run "$R/g27_text_numbers" "python3 text_numbers.py"
for s in depth_regress within_field bench_analyze fig_depth_v3; do run "$R/g23_smallest_hollows" "python3 $s.py"; done
[ -f "$ROOT/data/Concentrated_LRM_560_params_BW.tif" ] && run "$R/g23_smallest_hollows" "python3 depth_paper.py"
run "$R/g24_blewett_compare" "python3 cmp_analyze.py"
run "$R/g26_measurement_error" "python3 me.py"
run "$R/g26_measurement_error" "python3 no_shadow_check.py"
run "$R/g18_traced_r3/terrain" "python3 analyze.py out_v3/rows_tracedonly_r3_v3.csv"
run "$R/g19_knobs_traced/knobs" "python3 knob_ratio_design.py"
run "$P/analysis_v3" "python3 sulfur_charge.py"
run "$P/analysis_v3" "python3 sulfide_s2.py"
run "$P/analysis_v7" "python3 collapse_depth.py"
run "$P/analysis_v7" "python3 supplement_v7_calcs.py"
run "$P/analysis_v8" "python3 model_logic_v8_calcs.py"
run "$P/analysis_workshop/volatile_fit" "python3 depth_T_fit.py"
run "$P/analysis_workshop/volatile_fit" "python3 volatile_candidates.py"
exit $fail
