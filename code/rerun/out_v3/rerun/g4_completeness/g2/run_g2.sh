# g4 copy of code/rerun/out_v3/rerun/g2_insolation/run_g2.sh (6 Nov 2026). Changed: ROOT depth 5 -> 6.
#!/bin/bash
# g2 rerun driver (6 Nov 2026). Usage: G2_CAT=<v2|v3|v3b> bash run_g2.sh [script ...]
# Runs the copied scripts in dependency order; logs go to out/<cat>/logs/.
cd "$(dirname "$0")"
CAT=${G2_CAT:-v3}; export G2_CAT=$CAT
L=out/$CAT/logs; mkdir -p $L
ROOT=$(cd ../../../../../.. && pwd)
ALL="dT_hosting dT_hosting_completeness dT_plains transition_layer_hilat latcorr latcorr3 latcorr5 aspect_stats fractions_g2 hollow_rich_hosts band_crater_hotwarm tm_vs_dT_check stage_test stage_test_lowinc stage_fieldlevel cos2lon_lat15 size_quirk hotwarm threshold wall_fractions plains_hosts make_figs_pr_insolation_v9"
LIST=${@:-$ALL}
for s in $LIST; do
  t0=$(date +%s)
  case $s in
    dT_hosting_completeness) python3 $s.py "$ROOT" - > $L/$s.log 2>&1 ;;
    transition_layer_hilat) python3 $s.py "$ROOT/manuscripts/insolation_thermal/scripts" "inputs/$CAT/test2_craters_h18.csv" "$ROOT/manuscripts/PlanetaryResearch_hollows/analysis_v3/transition_layer_grid.npz" > $L/$s.log 2>&1 ;;
    *) python3 $s.py > $L/$s.log 2>&1 ;;
  esac
  echo "$CAT $s exit=$? $(( $(date +%s) - t0 ))s" >> out/$CAT/run_status.txt
done
echo "$CAT done" >> out/$CAT/run_status.txt
