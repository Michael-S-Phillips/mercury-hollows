#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script. Copies the color chain and applies only input/output path changes (each asserted once).
"""Copies analysis_workshop/color_mdr/{mdr_georef, extract_tiles, merge_spectra, pca_color, diagnostics, depth_color_lrm}.py into
color/, analysis_workshop/summary_model/step9_color.py into summary_model/, and g1's adapted copies of common.py,
step7_improved.py, cached_pool.py and run_step.sh (g1_graph_joint/summary_model) into summary_model/ unchanged.
The scripts are run with the working directory color/out_<cat>/, where tiles_qc is a link to _scratch/tiles_<cat>/.
Env: G6_MASK (hollow mask, relative to the project root), G6_TILES (tile output folder), G6_CF (crater table with h09),
G6_COLOR (folder with crater_color.csv and hollow_color.csv for step9_color)."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); AW = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop"
G1 = f"{HE}/code/rerun/out_v3/rerun/g1_graph_joint/summary_model"; HERE = os.path.dirname(os.path.abspath(__file__))
LOG = []
def copy(src, dst, reps, tag):
    t = open(src).read()
    for a, b in reps:
        assert t.count(a) == 1, (src, a, t.count(a)); t = t.replace(a, b)
    head = ""
    if t.startswith("#!"): head, t = t.split("\n", 1); head += "\n"
    cm = "#"
    open(dst, "w").write(head + f"{cm} g6_remaining copy of {src.replace(HE + '/', '')} (6 Nov 2026). Changed: {tag}\n" + t)
    LOG.append(f"{dst.replace(HERE + '/', '')} <- {src.replace(HE + '/', '')} ({len(reps)} replacements)")
C = f"{AW}/color_mdr"; D = f"{HERE}/color"
copy(f"{C}/mdr_georef.py", f"{D}/mdr_georef.py", [], "nothing (imported module).")
copy(f"{C}/extract_tiles.py", f"{D}/extract_tiles.py",
     [("OUT = 'tiles_qc'", "OUT = os.environ.get('G6_TILES', 'tiles_qc')"),
      ("CAN = 'code/rerun/out_v2/hollows_canonical.csv'", "CAN = os.environ.get('G6_MASK', 'code/rerun/out_v2/hollows_canonical.csv')")],
     "hollow-mask input (G6_MASK) and tile output folder (G6_TILES) only.")
copy(f"{C}/merge_spectra.py", f"{D}/merge_spectra.py", [], "nothing (reads and writes tiles_qc/ in the working directory, a link to _scratch/tiles_<cat>).")
copy(f"{C}/pca_color.py", f"{D}/pca_color.py",
     [("CF = os.path.join(BASE, 'manuscripts/PlanetaryResearch_hollows/analysis_workshop/summary_model/crater_features.csv')",
       "CF = os.environ.get('G6_CF', os.path.join(BASE, 'manuscripts/PlanetaryResearch_hollows/analysis_workshop/summary_model/crater_features.csv'))")],
     "crater table input (G6_CF: the catalog's crater_features.csv, for h09); outputs go to the working directory as before.")
copy(f"{C}/diagnostics.py", f"{D}/diagnostics.py",
     [("import json, numpy as np, pandas as pd", "import os, json, numpy as np, pandas as pd"),
      ("cf = pd.read_csv('../summary_model/crater_features.csv')",
       "cf = pd.read_csv(os.environ.get('G6_CF', os.path.expanduser('~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/summary_model/crater_features.csv')))"),
      ("pd.read_csv('../../../insolation_thermal/data/test1_rows_inc.csv')",
       "pd.read_csv(os.path.expanduser('~/mnt/hollows-emissivity/manuscripts/insolation_thermal/data/test1_rows_inc.csv'))")],
     "input paths (relative paths made absolute; G6_CF for h09); outputs go to the working directory as before.")
copy(f"{C}/depth_color_lrm.py", f"{D}/depth_color_lrm.py",
     [("HERE = os.path.dirname(os.path.abspath(__file__))", "HERE = os.getcwd()")],
     "hollow_color.csv input and record output in the working directory (color/out_<cat>) instead of the script folder.")
S = f"{HERE}/summary_model"
copy(f"{AW}/summary_model/step9_color.py", f"{S}/step9_color.py",
     [('U = "/mnt/user-data/uploads/hollows-emissivity/manuscripts"', 'U = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts")\nCOLOR = os.environ.get("G6_COLOR", f"{U}/PlanetaryResearch_hollows/analysis_workshop/color_mdr")'),
      ('pd.read_csv(f"{U}/PlanetaryResearch_hollows/analysis_workshop/color_mdr/crater_color.csv"', 'pd.read_csv(f"{COLOR}/crater_color.csv"'),
      ('pd.read_csv(f"{U}/PlanetaryResearch_hollows/analysis_workshop/color_mdr/hollow_color.csv")', 'pd.read_csv(f"{COLOR}/hollow_color.csv")')],
     "input paths only (project root as in g1's copy; G6_COLOR selects the color tables); output to SM_DIR/out as in g1's copy.")
for f in ["common.py", "step7_improved.py", "cached_pool.py"]:
    copy(f"{G1}/{f}", f"{S}/{f}", [], "nothing (g1's adapted copy, used unchanged).")
copy(f"{G1}/run_step.sh", f"{S}/run_step.sh", [], "nothing (g1's runner, used unchanged).")
open(f"{HERE}/make_copies_color.log", "w").write("\n".join(LOG) + "\n"); print("\n".join(LOG))
