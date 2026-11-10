# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/run_dag.sh (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/bin/bash
# g1 rerun of the containment graph (analysis_v9/run_v9.sh variants) for one catalog: bash run_dag.sh v3|v3b|v2 [main]
# Outputs: <CAT>/ (terrestrial uplift), <CAT>/uplift_lunar, <CAT>/kring_0p8, <CAT>/kring_1p3, <CAT>/caloris_orgel.
R=$HOME/mnt/hollows-emissivity; H=$(cd "$(dirname "$0")" && pwd); C=$1; export CATALOG=$C
mkdir -p $H/$C
CRATER_DAG_OUT=$H/$C python3 $H/crater_dag.py $R > $H/$C/crater_dag.log 2>&1
CRATER_DAG_OUT=$H/$C python3 $H/crater_dag_tests.py $R > $H/$C/crater_dag_tests.log 2>&1
CRATER_DAG_OUT=$H/$C python3 $H/crater_dag_sens.py $R > $H/$C/crater_dag_sens.log 2>&1
[ "$2" == "main" ] && exit 0
for V in uplift_lunar kring_0p8 kring_1p3 caloris_orgel; do mkdir -p $H/$C/$V; done
CRATER_DAG_SU=lunar CRATER_DAG_OUT=$H/$C/uplift_lunar python3 $H/crater_dag.py $R > $H/$C/uplift_lunar/crater_dag.log 2>&1
CRATER_DAG_KRING=0.8 CRATER_DAG_OUT=$H/$C/kring_0p8 python3 $H/crater_dag.py $R > $H/$C/kring_0p8/crater_dag.log 2>&1
CRATER_DAG_KRING=1.3 CRATER_DAG_OUT=$H/$C/kring_1p3 python3 $H/crater_dag.py $R > $H/$C/kring_1p3/crater_dag.log 2>&1
CRATER_DAG_CALORIS=orgel CRATER_DAG_OUT=$H/$C/caloris_orgel python3 $H/crater_dag.py $R > $H/$C/caloris_orgel/crater_dag.log 2>&1
CRATER_DAG_OUT=$H/$C/caloris_orgel python3 $H/crater_dag_tests.py $R > $H/$C/caloris_orgel/crater_dag_tests.log 2>&1
CRATER_DAG_OUT=$H/$C/caloris_orgel python3 $H/crater_dag_sens.py $R > $H/$C/caloris_orgel/crater_dag_sens.log 2>&1
echo DONE > $H/$C/run_dag.done
