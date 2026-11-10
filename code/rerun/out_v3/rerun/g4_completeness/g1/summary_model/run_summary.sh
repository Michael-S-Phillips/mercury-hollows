# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/run_summary.sh (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/bin/bash
# g1 rerun of the summary (joint) model for one catalog: bash run_summary.sh v3|v3b. Run after build_features.py and add_xrs.py.
# Same order and bootstrap sizes as the originals: step7 (B=100), size-step check, step8 figure, step3 (B=200), step9 (B=60), step4.
H=$(cd "$(dirname "$0")" && pwd); C=$1; export SM_DIR=$H/$C CATALOG=$C; cd $H; mkdir -p $SM_DIR/out
for s in "step7_improved.py 100" "step7_sizestep_check.py" "step8_figure_improved_v2.py" "step3_dominance.py 200" "step9_color.py 60" "step4_trees.py"; do
  set -- $s; echo "$(date +%T) start $1" >> $SM_DIR/out/progress.log
  python3 $s > $SM_DIR/out/${1%.py}.log 2>&1; echo "$(date +%T) end $1 rc=$?" >> $SM_DIR/out/progress.log
done
echo DONE >> $SM_DIR/out/progress.log
