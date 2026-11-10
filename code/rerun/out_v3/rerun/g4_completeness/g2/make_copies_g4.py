#!/usr/bin/env python3
"""g4 (6 Nov 2026): copy the g2 rerun scripts (code/rerun/out_v3/rerun/g2_insolation, themselves path-only copies of the
original analysis scripts) into this folder. Changes, each asserted to match exactly the stated number of times:
  (1) folder depth: G2_ROOT is six levels up from here instead of five (path only);
  (2) completeness grid: common_null.py and hotwarm.py load ../completeness_10ppd_v3.npz (G2_CAT = v3r) or _v3b.npz (v3br)
      instead of code/rerun/out_v2/completeness_10ppd.npz; other catalogs keep the out_v2 grid;
  (3) C1 source: dT_hosting_completeness.py reads C1 from inputs/<cat>/crater_dag_nodes_C1.csv for v3r / v3br instead of
      analysis_v9/crater_dag_nodes.csv (same node ids; C1 resampled on the new grid by prep_inputs.py).
Everything else is verbatim. The first line of each copy names its source and the changes."""
import os
HERE = os.path.dirname(os.path.abspath(__file__)); G2 = os.path.abspath(os.path.join(HERE, "..", "..", "g2_insolation"))
DEPTH = ('*[".."] * 5', '*[".."] * 6')
GRID_OLD = 'np.load(DD + "completeness_10ppd.npz")'
GRID_NEW = ('np.load({"v3r": f"{G2_ROOT}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz", '
            '"v3br": f"{G2_ROOT}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3b.npz"}.get(G2_CAT, DD + "completeness_10ppd.npz"))')
NODES_OLD = 'pd.read_csv(f"{G2_ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_v9/crater_dag_nodes.csv", dtype={"id": str}, low_memory=False)'
NODES_NEW = 'pd.read_csv(f"{IN}/crater_dag_nodes_C1.csv" if G2_CAT in ("v3r", "v3br") else f"{G2_ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_v9/crater_dag_nodes.csv", dtype={"id": str}, low_memory=False)'
BLOCKED = ["dT_hosting", "dT_hosting_completeness", "dT_plains", "transition_layer_hilat", "latcorr", "latcorr3", "latcorr5",
           "aspect_stats", "hotwarm", "common_null", "make_figs_pr_insolation_v9"]
PLAIN = ["fractions_g2", "hollow_rich_hosts", "band_crater_hotwarm", "tm_vs_dT_check", "stage_test", "stage_test_lowinc",
         "stage_fieldlevel", "cos2lon_lat15", "size_quirk", "threshold", "wall_fractions", "plains_hosts"]
SPEC = {n: [DEPTH + (1,)] for n in BLOCKED}
SPEC.update({n: [] for n in PLAIN})
SPEC["common_null"].append((GRID_OLD, GRID_NEW, 1)); SPEC["hotwarm"].append((GRID_OLD, GRID_NEW, 1))
SPEC["dT_hosting_completeness"].append((NODES_OLD, NODES_NEW, 1))
log = []
for n, reps in SPEC.items():
    txt = open(f"{G2}/{n}.py").read()
    for old, new, k in reps:
        assert txt.count(old) == k, (n, old, txt.count(old)); txt = txt.replace(old, new)
    what = "; ".join({DEPTH[0]: "G2_ROOT depth 5 -> 6", GRID_OLD: "completeness grid for v3r/v3br", NODES_OLD: "C1 source for v3r/v3br"}[o] for o, _, _ in reps) or "none (verbatim)"
    txt = f"# g4 copy of code/rerun/out_v3/rerun/g2_insolation/{n}.py (6 Nov 2026). Changed: {what} (see make_copies_g4.py).\n" + txt
    open(f"{HERE}/{n}.py", "w").write(txt); log.append(f"{n}.py: {what}")
sh = open(f"{G2}/run_g2.sh").read(); o = "ROOT=$(cd ../../../../.. && pwd)"; assert sh.count(o) == 1
sh = sh.replace(o, "ROOT=$(cd ../../../../../.. && pwd)")
open(f"{HERE}/run_g2.sh", "w").write("# g4 copy of code/rerun/out_v3/rerun/g2_insolation/run_g2.sh (6 Nov 2026). Changed: ROOT depth 5 -> 6.\n" + sh)
log.append("run_g2.sh: ROOT depth 5 -> 6")
open(f"{HERE}/make_copies_g4.log", "w").write("\n".join(log) + "\n"); print("\n".join(log))
