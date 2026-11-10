# g13 copy of g11_curvature/terrain/terrain_obs_exp.py (6 Nov 2026). Changed: folder g13_traced/terrain; rows file rows_traced_<CLS>.csv (env CURV_CFG, default traced_v3).
# g11 copy of g10_threshold/terrain_obs_exp.py (6 Nov 2026). Changed: folder g13_traced/terrain/out_<CAT>; reads rows_curv_<cfg>.csv (cfg from env CURV_CFG, default n24_h100_s0.03); G10_MODE removed.
# g10 copy of g7_terrain_fix/terrain_obs_exp.py: output folder g10_threshold/out_<CAT>_<G10_MODE> only.
# g7 copy of g3_catalog/terrain/terrain_obs_exp.py (6 Nov 2026): output folder g7_terrain_fix/out_<CAT>_grid<G7_GRID> only.
#!/usr/bin/env python3
"""g3_catalog (6 Nov 2026): the last column of Table S2 (observed/expected hollows by terrain, each crater's peak+floor+wall hollows
spread over its completeness-weighted peak+floor+wall area in proportion to area; 95% intervals from 2,000 bootstrap resamples of
craters, seed 1) and the all-terrain density used for the shading (v2: 6.1 per 1,000 km2), plus the setting ratios quoted in S4
(isolated peak/wall; nested wall / isolated wall). No script for these was found in the project, so this one was written and
checked against the v2 table. New script. Env G3_CAT. Reads terrain/out_{cat}/rows_H250_r0.45.csv; writes terrain/out_{cat}/terrain_obs_exp.json."""
import os, json, numpy as np, pandas as pd
G3 = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced/terrain"); CAT = os.environ.get("G3_CAT", "v3"); O = f"{G3}/out_{CAT}"; CFG = os.environ.get("CURV_CFG", "traced_v3")
r = pd.read_csv(f"{O}/rows_{CFG}.csv"); r = r[(r.D >= 30) & (r.D <= 250)].reset_index(drop=True)
T = ("peak", "floor", "wall"); n = r[[f"n_{t}" for t in T]].values.astype(float); A = r[[f"A_{t}" for t in T]].values
ok = A.sum(1) > 0; E = np.zeros_like(n); E[ok] = n[ok].sum(1, keepdims=True) * A[ok] / A[ok].sum(1, keepdims=True)
def oe(i): return n[i].sum(0) / E[i].sum(0)
rng = np.random.default_rng(1); bs = np.array([oe(rng.integers(0, len(r), len(r))) for _ in range(2000)])
nested = r.parent_D > 0; iso = (~nested) & (~r.in_orgel_basin)
dens = lambda s, t: 1000 * s[f"n_{t}"].sum() / s[f"A_{t}"].sum()
out = dict(N=int(len(r)), all_terrain_density_per_1000km2=float(1000 * n.sum() / A.sum()),
           obs_over_exp={t: dict(ratio=float(oe(np.arange(len(r)))[k]), lo=float(np.percentile(bs[:, k], 2.5)), hi=float(np.percentile(bs[:, k], 97.5))) for k, t in enumerate(T)},
           isolated_peak_over_wall=float(dens(r[iso], "peak") / dens(r[iso], "wall")),
           nested_wall_over_isolated_wall=float(dens(r[nested], "wall") / dens(r[iso], "wall")),
           nested_densest=max(T, key=lambda t: dens(r[nested], t)))
json.dump(out, open(f"{O}/terrain_obs_exp.json" if CFG == "traced_v3" else f"{O}/terrain_obs_exp_{CFG}.json", "w"), indent=1); print(CAT, json.dumps(out))
