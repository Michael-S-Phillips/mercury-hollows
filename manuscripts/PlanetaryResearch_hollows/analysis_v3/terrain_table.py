import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Hollow density by terrain class and crater setting for the 352 host craters 30-250 km across (Table 2).
Reads peak_wall_test/rows_H250_r0.45.csv (per-crater completeness-weighted areas and hollow counts from peak_wall.py).
Settings: nested = center inside a cataloged crater >= 2x larger; basin-hosted = otherwise inside an Orgel et al. (2020)
certain/probable basin; isolated = neither. 95% intervals from 2,000 bootstrap resamples of craters (seed 1)."""
import os, json, numpy as np, pandas as pd
HE = os.path.expanduser(_REPO + "")
r = pd.read_csv(f"{HE}/peak_wall_test/rows_H250_r0.45.csv"); r = r[(r.D >= 30) & (r.D <= 250)]
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
json.dump(out, open(f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_v3/terrain_table_record.json", "w"), indent=1)
print(json.dumps({k: {t: round(v[t]["density"], 2) for t in ("peak", "floor", "wall")} for k, v in out.items()}))
