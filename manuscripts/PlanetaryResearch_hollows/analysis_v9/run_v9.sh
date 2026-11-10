#!/bin/bash
# v9 runs (argument: hollows-emissivity root). Each writes into its own folder.
R=${1:-$HOME/mnt/hollows-emissivity}; H=$(cd "$(dirname "$0")" && pwd)
CRATER_DAG_BASIN=v8 CRATER_DAG_OUT=$H/check_v8_rule python3 $H/crater_dag.py $R > $H/check_v8_rule/crater_dag.log 2>&1
python3 $H/crater_dag.py $R > $H/crater_dag_v9.log 2>&1
CRATER_DAG_SU=lunar CRATER_DAG_OUT=$H/uplift_lunar python3 $H/crater_dag.py $R > $H/uplift_lunar/crater_dag.log 2>&1
CRATER_DAG_KRING=0.8 CRATER_DAG_OUT=$H/kring_0p8 python3 $H/crater_dag.py $R > $H/kring_0p8/crater_dag.log 2>&1
CRATER_DAG_KRING=1.3 CRATER_DAG_OUT=$H/kring_1p3 python3 $H/crater_dag.py $R > $H/kring_1p3/crater_dag.log 2>&1
python3 $H/crater_dag_tests.py $R > $H/crater_dag_tests_v9.log 2>&1
python3 $H/crater_dag_sens.py $R > $H/crater_dag_sens_v9.log 2>&1
echo ALL_DONE > $H/run_v9.done
