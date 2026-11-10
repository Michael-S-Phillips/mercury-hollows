#!/usr/bin/env python3
"""Copies the Caloris scripts into caloris/ and changes only input paths, output paths and env settings.
Env G3_CAT = v2 | v3 | v3b. Host assignments: GRL copy (v2) or code/rerun/out_v3/host_assignments_{v}.csv; containment-graph nodes:
inputs/crater_dag_nodes_{v}.csv (prep_inputs.py). Outputs: caloris/out_{G3_CAT}/."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); DST = os.path.dirname(os.path.abspath(__file__)) + "/caloris"
HDR = ('import os as _os\nG3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")\n'
       'CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/caloris/out_{CAT}"; _os.makedirs(ODIR, exist_ok=True)\n'
       'HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/GRL_insolation_hollows/data/host_assignments.csv") if CAT == "v2"\n'
       '         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))\n')
def copy(src, name, reps, after):
    t = open(src).read()
    for a, b in reps:
        assert t.count(a) == 1, (name, a, t.count(a)); t = t.replace(a, b)
    assert t.count(after) == 1, after; t = t.replace(after, after + HDR, 1)
    open(f"{DST}/{name}", "w").write(f"# g3_catalog copy of {src.replace(HE + '/', '')} (6 Nov 2026): only input/output paths and env settings changed.\n" + t)
C9 = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/caloris_v9"
copy(f"{C9}/caloris_v9.py", "caloris_v9.py",
     [('OUT = os.path.dirname(os.path.abspath(__file__))', 'OUT = ODIR'),
      ('n = pd.read_csv(f"{P}/analysis_v9/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str})', 'n = pd.read_csv(f"{G3}/inputs/crater_dag_nodes_{CAT}.csv", dtype={"id": str, "parent_id": str})'),
      ('h = pd.read_csv(f"{G}/data/host_assignments.csv", usecols=["host_uid", "rho"])', 'h = pd.read_csv(HOSTS, usecols=["host_uid", "rho"])'),
      ('host_assignments=f"{G}/data/host_assignments.csv"', 'host_assignments=HOSTS')],
     "from scipy.stats import fisher_exact, mannwhitneyu\n")
copy(f"{C9}/caloris_v9_controls.py", "caloris_v9_controls.py",
     [('H = os.path.dirname(os.path.abspath(__file__))', 'H = ODIR')],
     "from scipy.stats import fisher_exact\n")
PW = f"{HE}/peak_wall_test"
copy(f"{PW}/caloris_test.py", "caloris_test.py",
     [('ch=pd.read_csv("nested_children_classified.csv"); m=pd.read_csv("morph_all_nested.csv")',
       'ch=pd.read_csv(_os.path.expanduser("~/mnt/hollows-emissivity/peak_wall_test/nested_children_classified.csv")); m=pd.read_csv(_os.path.expanduser("~/mnt/hollows-emissivity/peak_wall_test/morph_all_nested.csv"))'),
      ('ha=pd.read_csv("../code/rerun/latitude/aspect/host_assignments.csv")', 'ha=pd.read_csv(HOSTS)'),
      ('cal.to_csv("caloris_children.csv",index=False)', 'cal.to_csv(f"{ODIR}/caloris_children.csv",index=False)'),
      ('open("caloris_test_record.json","w")', 'open(f"{ODIR}/caloris_test_record.json","w")')],
     "import statsmodels.formula.api as smf\n")
print("copied")
