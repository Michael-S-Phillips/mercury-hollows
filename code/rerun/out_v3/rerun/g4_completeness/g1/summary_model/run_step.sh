# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/run_step.sh (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/bin/bash
# bash run_step.sh <CATALOG> <script.py> [args]: one resumable call of cached_pool.py (exit 3 = call again; 0 = done).
H=$(cd "$(dirname "$0")" && pwd); C=$1; shift; export SM_DIR=$H/$C CATALOG=$C; cd $H; mkdir -p $SM_DIR/out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BUDGET=${BUDGET:-90} NPROC=${NPROC:-3}
S=${1%.py}; timeout 125 python3 cached_pool.py "$@" > $SM_DIR/out/$S.log.part 2> $SM_DIR/out/$S.err.part; rc=$?
cat $SM_DIR/out/$S.err.part >> $SM_DIR/out/$S.err; [ $rc -eq 0 ] && mv $SM_DIR/out/$S.log.part $SM_DIR/out/$S.log
echo "$(date +%T) $C $S rc=$rc" >> $SM_DIR/out/progress.log; echo "rc=$rc"; tail -n 2 $SM_DIR/out/$S.err.part
