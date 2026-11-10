# g13 copy of g11_curvature/terrain/merge_rows.py (6 Nov 2026). Changed: folder g13_traced/terrain; rows file rows_traced_<CLS>.csv (env CURV_CFG, default traced_v3).
# g11 copy of g10_threshold/merge_rows.py (6 Nov 2026). Changed: folder g13_traced/terrain/out_<CAT> and chunk names rows_curv_<cfg>_START_STOP.csv; G10_MODE removed.
# g10 copy of g7_terrain_fix/merge_rows.py: output folder g10_threshold/out_<CAT>_<G10_MODE> only.
# g7 copy of g3_catalog/terrain/merge_rows.py (6 Nov 2026): output folder g7_terrain_fix/out_<CAT>_grid<G7_GRID> only.
#!/usr/bin/env python3
"""g3_catalog: concatenates peak_wall.py chunks (out_{cat}/_chunks/rows_H*_r*_START_STOP.csv) into out_{cat}/rows_H*_r*.csv, the
file name analyze.py and terrain_table.py expect. Env G3_CAT. New script."""
import os, glob, re, pandas as pd
G3 = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced/terrain"); CAT = os.environ.get("G3_CAT", "v3"); O = f"{G3}/out_{CAT}"
groups = {}
for f in sorted(glob.glob(f"{O}/_chunks/rows_*_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].csv")):
    groups.setdefault(re.sub(r"_\d{4}_\d{4}\.csv$", ".csv", os.path.basename(f)), []).append(f)
for k, fs in groups.items():
    d = pd.concat([pd.read_csv(f) for f in fs]).drop_duplicates("uid").sort_values("uid"); d.to_csv(f"{O}/{k}", index=False); print(k, len(fs), len(d))
