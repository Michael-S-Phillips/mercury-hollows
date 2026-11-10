#!/usr/bin/env python3
"""Copies melt_source_tests/positions/positions_build.py and positions_analyze.py into melt_source/positions/ and changes only input
paths, output paths and env settings (env G3_CAT = v2 | v3 | v3b; outputs in melt_source/positions/out_{G3_CAT}/)."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); SRC = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/melt_source_tests/positions"
DST = os.path.dirname(os.path.abspath(__file__)) + "/melt_source/positions"; os.makedirs(DST, exist_ok=True)
HDR = ('import os as _os\nG3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")\n'
       'CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/melt_source/positions/out_{CAT}"; _os.makedirs(ODIR + "/cache", exist_ok=True)\n'
       'HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"\n'
       '         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))\n'
       'CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"\n'
       '         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))\n')
def copy(name, reps, after):
    t = open(f"{SRC}/{name}").read()
    for a, b in reps:
        assert t.count(a) == 1, (name, a, t.count(a)); t = t.replace(a, b)
    assert t.count(after) == 1, after; t = t.replace(after, after + HDR, 1)
    open(f"{DST}/{name}", "w").write(f"# g3_catalog copy of analysis_workshop/melt_source_tests/positions/{name} (6 Nov 2026): only input/output paths and env settings changed.\n" + t)
copy("positions_build.py", [('OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")', 'OUT = ODIR + "/cache"'),
                            ('ha = pd.read_csv(f"{HE}/code/rerun/latitude/aspect/host_assignments.csv")', 'ha = pd.read_csv(HOSTS)'),
                            ('hc = pd.read_csv(f"{HE}/code/rerun/out_v2/hollows_canonical.csv")[["hollow_id", "field_id", "is_crater_trace"]]',
                             'hc = pd.read_csv(CANON)[["hollow_id", "field_id", "is_crater_trace"]]')],
     "from scipy.spatial import cKDTree\n")
copy("positions_analyze.py", [('HERE = os.path.dirname(os.path.abspath(__file__)); HE', 'HERE = ODIR; HE')], "import os, json, numpy as np, pandas as pd\n")
print("copied")
