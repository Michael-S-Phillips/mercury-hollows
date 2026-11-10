#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script. Copies the S2 reference-size checks and changes only input/output paths (each asserted once).
"""Copies analysis_workshop/completeness_ref/{ref_size_check, ref_size_insolation}.py and analysis_workshop/completeness_sfd/
{hosting_fits, joint_thinned}.py into refsize/. Env G6_CAT = v2 | v3 | v3b: the crater table with the hosting flags is g1's
g1_graph_joint/summary_model/<cat>/crater_features.csv (v2: byte-identical to the original), outputs go to refsize/out_<cat>/.
Completeness grids (code/rerun/out_v2/completeness_10ppd.npz) and completeness_sfd/crater_completeness.csv are read unchanged
(they do not depend on the screen; the E_a measures are given for exponents 1.3, 1.6 and 1.9, which bracket the v3 exponent 1.7).
joint_thinned.py imports this folder's summary_model/step7_improved.py (g1's copy), which reads SM_DIR = summary_model/<cat>."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); AW = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop"
HERE = os.path.dirname(os.path.abspath(__file__)); DST = f"{HERE}/refsize"; LOG = []
CF = 'os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g1_graph_joint/summary_model/" + os.environ.get("G6_CAT", "v3") + "/crater_features.csv")'
OUT = 'os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_" + os.environ.get("G6_CAT", "v3"))'
def copy(src, name, reps, tag):
    t = open(src).read()
    for a, b in reps:
        assert t.count(a) == 1, (name, a, t.count(a)); t = t.replace(a, b)
    head = ""
    if t.startswith("#!"): head, t = t.split("\n", 1); head += "\n"
    open(f"{DST}/{name}", "w").write(head + f"# g6_remaining copy of {src.replace(HE + '/', '')} (6 Nov 2026). Changed: {tag}\n" + t)
    LOG.append(f"refsize/{name} <- {src.replace(HE + '/', '')} ({len(reps)} replacements)")
copy(f"{AW}/completeness_ref/ref_size_check.py", "ref_size_check.py",
     [('d = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", dtype={"id": str})', f'd = pd.read_csv({CF}, dtype={{"id": str}})'),
      ('json.dump(R, open(f"{HERE}/ref_size_record.json", "w"), indent=1)', f'os.makedirs({OUT}, exist_ok=True); json.dump(R, open(os.path.join({OUT}, "ref_size_record.json"), "w"), indent=1)')],
     "crater table input (G6_CAT) and output folder only.")
copy(f"{AW}/completeness_ref/ref_size_insolation.py", "ref_size_insolation.py",
     [('HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, f"{HERE}/../insolation_v6")', f'HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, "{AW}/insolation_v6")'),
      ('c = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})', f'c = pd.read_csv({CF}, low_memory=False, dtype={{"id": str}})'),
      ('json.dump(R, open(f"{HERE}/ref_size_insolation.json", "w"), indent=1)', f'os.makedirs({OUT}, exist_ok=True); json.dump(R, open(os.path.join({OUT}, "ref_size_insolation.json"), "w"), indent=1)')],
     "import path of common_null (original folder), crater table input (G6_CAT) and output folder only.")
copy(f"{AW}/completeness_sfd/hosting_fits.py", "hosting_fits.py",
     [('HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, f"{HERE}/../insolation_v6")', f'HERE = "{AW}/completeness_sfd"; sys.path.insert(0, "{AW}/insolation_v6")'),
      ('d = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})', f'd = pd.read_csv({CF}, low_memory=False, dtype={{"id": str}})'),
      ('json.dump(R, open(f"{HERE}/hosting_fits.json", "w"), indent=1)', f'os.makedirs({OUT}, exist_ok=True); json.dump(R, open(os.path.join({OUT}, "hosting_fits.json"), "w"), indent=1)')],
     "input folder (crater_completeness.csv read from the original folder), crater table input (G6_CAT) and output folder only.")
copy(f"{AW}/completeness_sfd/joint_thinned.py", "joint_thinned.py",
     [('HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, f"{HERE}/../summary_model")',
       f'HERE = "{AW}/completeness_sfd"; sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "summary_model"))'),
      ('out_p = f"{HERE}/joint_thinned.json"', f'os.makedirs({OUT}, exist_ok=True); out_p = os.path.join({OUT}, "joint_thinned.json")')],
     "step7_improved import (this folder's summary_model/, g1's copy; SM_DIR selects the catalog), input folder of crater_completeness.csv and output path only.")
open(f"{HERE}/make_copies_refsize.log", "w").write("\n".join(LOG) + "\n"); print("\n".join(LOG))
