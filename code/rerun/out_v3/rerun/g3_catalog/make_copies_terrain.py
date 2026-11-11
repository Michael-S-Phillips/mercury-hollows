#!/usr/bin/env python3
"""Copies peak_wall.py, analyze.py (peak_wall_test) and terrain_table.py (analysis_v3) into terrain/ and changes only input paths,
output paths and env settings. Env G3_CAT = v2 | v3 | v3b. Host assignments: code/rerun/latitude/aspect/host_assignments.csv (v2;
same bytes as the insolation_thermal copy) or code/rerun/out_v3/host_assignments_{v}.csv; field ids from the matching canonical file.
Outputs: terrain/out_{G3_CAT}/."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); DST = os.path.dirname(os.path.abspath(__file__)) + "/terrain"
HDR = ('import os as _os\nG3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")\n'
       'CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/terrain/out_{CAT}"; _os.makedirs(ODIR, exist_ok=True)\n'
       'HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"\n'
       '         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))\n'
       'CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"\n'
       '         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))\n')
def copy(src, name, reps, after):
    t = open(src).read()
    for a, b in reps:
        assert t.count(a) == 1, (name, a, t.count(a)); t = t.replace(a, b)
    assert t.count(after) == 1, after; t = t.replace(after, after + HDR, 1)
    open(f"{DST}/{name}", "w").write(f"# g3_catalog copy of {src.replace(HE + '/', '')} (6 Nov 2026): only input/output paths and env settings changed.\n" + t)
copy(f"{HE}/peak_wall_test/peak_wall.py", "peak_wall.py",
     [('ha=pd.read_csv(f"{HE}/code/rerun/latitude/aspect/host_assignments.csv")', 'ha=pd.read_csv(HOSTS)'),
      ('hc=pd.read_csv(f"{HE}/code/rerun/out_v2/hollows_canonical.csv")[["hollow_id","field_id"]]', 'hc=pd.read_csv(CANON)[["hollow_id","field_id"]]'),
      ('out=f"{HE}/peak_wall_test/rows_H{int(HT)}_r{RPK}.csv"', '_os.makedirs(f"{ODIR}/_chunks", exist_ok=True); out=f"{ODIR}/_chunks/rows_H{int(HT)}_r{RPK}_{a0:04d}_{a1:04d}.csv"')],
     "from rasterio.windows import from_bounds\n")
copy(f"{HE}/peak_wall_test/analyze.py", "analyze.py",
     [('open("peak_wall_record.json","w")', 'open(f"{ODIR}/peak_wall_record.json","w")')], "from scipy import stats\n")
copy(f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_v3/terrain_table.py", "terrain_table.py",
     [('pd.read_csv(f"{HE}/peak_wall_test/rows_H250_r0.45.csv")', 'pd.read_csv(f"{ODIR}/rows_H250_r0.45.csv")'),
      ('open(f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_v3/terrain_table_record.json", "w")', 'open(f"{ODIR}/terrain_table_record.json", "w")')],
     "import os, json, numpy as np, pandas as pd\n")
print("copied")
