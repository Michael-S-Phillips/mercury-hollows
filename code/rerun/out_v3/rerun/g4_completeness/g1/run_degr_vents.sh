#!/bin/bash
# g4 (6 Nov 2026), new driver: degradation-state and vent tests of g1 for one catalog: bash run_degr_vents.sh <v2|v3|v3b|v3r|v3br>
# (order as in g1_graph_joint/README.md; the summary_model/<CAT>/crater_features.csv must exist).
H=$(cd "$(dirname "$0")" && pwd); export CATALOG=$1; mkdir -p $H/degradation/$1 $H/vents/$1
for s in degradation/match_kinczyk degradation/kinczyk_tests degradation/stage_by_class degradation/kinczyk_match_sensitivity degradation/host_age_tests vents/test8_vents vents/vent_assignment_check; do
  t0=$(date +%s); python3 $H/$s.py > $H/${s%/*}/$1/$(basename $s).log 2>&1; echo "$1 $s rc=$? $(( $(date +%s) - t0 ))s"
done
