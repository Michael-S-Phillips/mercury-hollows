# g18 copy of g13_traced/terrain/terrain_table.py (7 Nov 2026). Changed: folder g18_traced_r3; traces from traces/traces_clean_r3.json (round-2 snapshot, hosts of tracer groups A and B); G13_REV default r3 (slump counts as wall, as in r2); outputs named r3.
# g13 copy of g11_curvature/terrain/terrain_table.py (6 Nov 2026). Changed: folder g13_traced/terrain; rows file rows_traced_<CLS>.csv (env CURV_CFG, default traced_v3).
# g11 copy of g10_threshold/terrain_table.py (6 Nov 2026). Changed: folder g13_traced/terrain/out_<CAT>; reads rows_curv_<cfg>.csv (cfg from env CURV_CFG, default n24_h100_s0.03); G10_MODE removed.
# g10 copy of g7_terrain_fix/terrain_table.py: output folder g10_threshold/out_<CAT>_<G10_MODE> only.
# g7 copy of g3_catalog/terrain/terrain_table.py (6 Nov 2026): output folder g7_terrain_fix/out_<CAT>_grid<G7_GRID> only.
# g3_catalog copy of manuscripts/PlanetaryResearch_hollows/analysis_v3/terrain_table.py (6 Nov 2026): only input/output paths and env settings changed.
"""Hollow density by terrain class and crater setting for the 352 host craters 30-250 km across (Table 2).
Reads peak_wall_test/rows_H250_r0.45.csv (per-crater completeness-weighted areas and hollow counts from peak_wall.py).
Settings: nested = center inside a cataloged crater >= 2x larger; basin-hosted = otherwise inside an Orgel et al. (2020)
certain/probable basin; isolated = neither. 95% intervals from 2,000 bootstrap resamples of craters (seed 1)."""
import os, json, numpy as np, pandas as pd
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g18_traced_r3/terrain")
CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/out_{CAT}"; CFG = _os.environ.get("CURV_CFG", "traced_v3"); _os.makedirs(ODIR, exist_ok=True)
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))
HE = os.path.expanduser("~/mnt/hollows-emissivity")
r = pd.read_csv(f"{ODIR}/rows_{CFG}.csv"); r = r[(r.D >= 30) & (r.D <= 250)]
def dens(s, t): return s[f"n_{t}"].sum() / s[f"A_{t}"].sum() * 1000
nested = r.parent_D > 0; basin = (~nested) & r.in_orgel_basin; iso = (~nested) & (~r.in_orgel_basin)
rng = np.random.default_rng(1); out = {}
for nm, s in [("all", r), ("isolated", r[iso]), ("basin", r[basin]), ("nested", r[nested])]:
    o = {"N": int(len(s))}
    for t in ("peak", "floor", "wall", "wall_eq", "wall_pole"):
        bs = [dens(s.iloc[rng.integers(0, len(s), len(s))], t) for _ in range(2000)]
        o[t] = dict(density=dens(s, t), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)),
                    n=int(s[f"n_{t}"].sum()), area_km2=float(s[f"A_{t}"].sum()))
    out[nm] = o
json.dump(out, open(f"{ODIR}/terrain_table_record.json" if CFG == "traced_v3" else f"{ODIR}/terrain_table_record_{CFG}.json", "w"), indent=1)
print(json.dumps({k: {t: round(v[t]["density"], 2) for t in ("peak", "floor", "wall")} for k, v in out.items()}))
