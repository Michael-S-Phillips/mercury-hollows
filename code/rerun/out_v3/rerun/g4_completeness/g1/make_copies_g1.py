#!/usr/bin/env python3
"""g4 (6 Nov 2026): copy the g1 rerun scripts (code/rerun/out_v3/rerun/g1_graph_joint, themselves path-only copies of the
original analysis scripts) into this folder, keeping the subfolders. Changes, each asserted to match exactly once:
  hostfiles.py : CATALOG may also be v3r (v3 host assignments and crater table) or v3br (v3b ones); the crater tables are read
                 from g1_graph_joint/inputs/ (input path only).
  crater_dag.py: the completeness grid is ../completeness_10ppd_v3.npz for CATALOG = v3r and _v3b.npz for v3br (the refit
                 detection model); v2, v3, v3b keep code/rerun/out_v2/completeness_10ppd.npz.
Everything else is verbatim (outputs go to this folder because the scripts write next to themselves). Not copied:
make_changes.py (see make_changes_g4.py), compare_v8_v9.py and show_dag.py (no quoted numbers), inputs/make_inputs.py
(inputs are reused), depth_nested/ (Blewett depths; no completeness)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__)); G1 = os.path.abspath(os.path.join(HERE, "..", "..", "g1_graph_joint"))
FILES = ["crater_dag.py", "crater_dag_tests.py", "crater_dag_sens.py", "run_dag.sh", "hostfiles.py", "parents_and_examples.py",
         "degradation/host_age_tests.py", "degradation/kinczyk_match_sensitivity.py", "degradation/kinczyk_tests.py",
         "degradation/match_kinczyk.py", "degradation/stage_by_class.py", "vents/test8_vents.py", "vents/vent_assignment_check.py"] + \
        [f"summary_model/{f}" for f in ("add_xrs.py", "boot_failure_check.py", "build_features.py", "cached_pool.py", "common.py",
         "count_model_check.py", "run_step.sh", "run_summary.sh", "step3_dominance.py", "step4_trees.py", "step7_improved.py",
         "step7_sizestep_check.py", "step8_figure_improved_v2.py", "step9_color.py")]
SPEC = {
 "hostfiles.py": [
  ('CAT = os.environ.get("CATALOG", "v3"); assert CAT in ("v2", "v3", "v3b"), CAT',
   'CAT = os.environ.get("CATALOG", "v3"); assert CAT in ("v2", "v3", "v3b", "v3r", "v3br"), CAT\n'
   'BASE = {"v3r": "v3", "v3br": "v3b"}.get(CAT, CAT)   # g4: v3r / v3br use the v3 / v3b catalog inputs'),
  ('"v3": f"{ROOT}/code/rerun/out_v3/host_assignments_v3.csv", "v3b": f"{ROOT}/code/rerun/out_v3/host_assignments_v3b.csv"}[CAT]',
   '"v3": f"{ROOT}/code/rerun/out_v3/host_assignments_v3.csv", "v3b": f"{ROOT}/code/rerun/out_v3/host_assignments_v3b.csv"}[BASE]'),
  ('CRATERS = f"{HERE}/inputs/test2_craters_h18_{CAT}.csv"', 'CRATERS = f"{ROOT}/code/rerun/out_v3/rerun/g1_graph_joint/inputs/test2_craters_h18_{BASE}.csv"')],
 "crater_dag.py": [
  ('Cz = np.load(f"{ROOT}/code/rerun/out_v2/completeness_10ppd.npz")',
   'Cz = np.load({"v3r": f"{ROOT}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz", "v3br": f"{ROOT}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3b.npz"}.get(os.environ.get("CATALOG", "v3"), f"{ROOT}/code/rerun/out_v2/completeness_10ppd.npz"))')],
}
log = []
for f in FILES:
    txt = open(f"{G1}/{f}").read(); reps = SPEC.get(f, [])
    for o, n in reps:
        assert txt.count(o) == 1, (f, o, txt.count(o)); txt = txt.replace(o, n)
    what = {"hostfiles.py": "v3r/v3br catalogs map to the v3/v3b inputs; crater tables read from g1_graph_joint/inputs",
            "crater_dag.py": "completeness grid of the refit model for CATALOG = v3r/v3br"}.get(f, "none (verbatim)")
    txt = f"# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/{f} (6 Nov 2026). Changed: {what} (see make_copies_g1.py).\n" + txt
    os.makedirs(os.path.dirname(f"{HERE}/{f}"), exist_ok=True); open(f"{HERE}/{f}", "w").write(txt); log.append(f"{f}: {what}")
open(f"{HERE}/make_copies_g1.log", "w").write("\n".join(log) + "\n"); print("\n".join(log))
