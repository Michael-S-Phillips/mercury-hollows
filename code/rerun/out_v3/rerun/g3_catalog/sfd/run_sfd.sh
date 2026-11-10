#!/bin/bash
# Runs the copied size-distribution chain for one catalog and detection-model setting. Usage: run_sfd.sh CAT THETA
cd "$(dirname "$0")"
export G3_CAT=$1 G3_THETA=$2
mkdir -p out_${1}_${2}
if [ "$2" = "refit" ]; then python3 refit_detection.py > out_${1}_${2}/refit_detection.log 2>&1; fi
python3 fit_sfd.py > out_${1}_${2}/fit_sfd.log 2>&1 && python3 sfd_windows.py > out_${1}_${2}/sfd_windows.log 2>&1 && python3 sfd_numbers.py > out_${1}_${2}/sfd_numbers.log 2>&1
echo done > out_${1}_${2}/DONE
