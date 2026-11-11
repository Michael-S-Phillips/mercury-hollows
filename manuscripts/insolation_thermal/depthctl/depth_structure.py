#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Where does hollow depth vary? Crater-level vs regional ('provincial') structure, and position within the host crater.
Inputs: data/test1_rows_inc.csv (Blewett 2016 depths + incidence), Herrick 2018 craters. Writes depthctl/depth_structure_record.json."""
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.spatial import cKDTree
G = _REPO + "/manuscripts/insolation_thermal"
d = pd.read_csv(_REPO + "/manuscripts/insolation_thermal/data/test1_rows_inc.csv")
d["ltan"] = np.log10(np.tan(np.radians(d.inc))); d["lp"] = np.log10(d.px); d["ld"] = np.log10(d.depth)
c = smf.ols("ld ~ lp + ltan", data=d).fit(); d["res"] = c.resid           # corrected log10 depth (anomaly)
cr = pd.read_csv(f"{G}/data/test2_craters_h18.csv"); cr = cr[cr.diameter >= 5].reset_index(drop=True)
R_M = 2439.7
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
tree = cKDTree(xyz(cr.latitude.values, cr.lon360.values)); P = xyz(d.latitude.values, d.longitude.values % 360)
Dkm = cr.diameter.values; host = np.full(len(d), -1); rho = np.full(len(d), np.nan)
for i, cl in enumerate(tree.query_ball_point(P, r=0.36)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(P[i] @ tree.data[cl].T, -1, 1)) * R_M
    ok = ang <= 1.1 * Dkm[cl] / 2
    if ok.any():
        j = np.argmin(np.where(ok, Dkm[cl], np.inf)); host[i] = cl[j]; rho[i] = ang[j] / (Dkm[cl[j]] / 2)
d["host"] = host; d["rho"] = rho; d["hostD"] = np.where(host >= 0, Dkm[np.maximum(host, 0)], np.nan)
R = {"n_hollows": int(len(d)), "n_assigned": int((d.host >= 0).sum())}
# ---- crater medians
h = d[d.host >= 0]; g = h.groupby("host")
pc = g.agg(n=("res", "size"), med=("res", "median"), nimg=("filename", "nunique"), D=("hostD", "first"))
pc = pc[pc.n >= 3].copy()
pc["lat"] = cr.latitude.values[pc.index]; pc["lon"] = cr.lon360.values[pc.index]
R["n_craters"] = int(len(pc)); R["crater_median_sd_log10"] = float(pc.med.std())
R["crater_median_range_m"] = [float(10 ** (pc.med.quantile(.1))), float(10 ** (pc.med.quantile(.9)))]   # relative factor
# ---- variance shares (hollow level): crater, image-within-crater, residual
def r2(f, data): return float(smf.ols(f, data=data).fit().rsquared)
hh = h[h.host.isin(pc.index)].copy()
R["R2_crater"] = r2("res ~ C(host)", hh); R["R2_crater_image"] = r2("res ~ C(host) + C(filename)", hh)
R["R2_image_only"] = r2("res ~ C(filename)", hh)
# ---- spatial structure of crater medians: semivariogram and permutation null
X = xyz(pc.lat.values, pc.lon.values); ang = np.arccos(np.clip(X @ X.T, -1, 1)) * R_M
iu = np.triu_indices(len(pc), 1); dist = ang[iu]; m = pc.med.values
sv = 0.5 * (m[iu[0]] - m[iu[1]]) ** 2
bins = [0, 300, 600, 1000, 2000, 4000, 8000]
rng = np.random.default_rng(3)
def vario(mm):
    s = 0.5 * (mm[iu[0]] - mm[iu[1]]) ** 2
    return np.array([s[(dist >= a) & (dist < b)].mean() for a, b in zip(bins[:-1], bins[1:])])
obs = vario(m); null = np.array([vario(rng.permutation(m)) for _ in range(2000)])
R["variogram"] = [dict(bin=f"{a}-{b} km", n_pairs=int(((dist >= a) & (dist < b)).sum()), semivar=float(o),
                       null_mean=float(nl.mean()), p_lower=float((nl <= o).mean())) for (a, b), o, nl in zip(zip(bins[:-1], bins[1:]), obs, null.T)]
R["total_var_crater_medians"] = float(m.var())
# nearest-neighbor crater correlation within 600 km
near = (dist < 600)
R["pairs_lt600_corr"] = float(np.corrcoef(m[iu[0]][near], m[iu[1]][near])[0, 1]) if near.sum() > 5 else None
# ---- regional tiles (30 x 30 deg): share of crater-median variance explained
pc["tile"] = (np.floor(pc.lat / 30)).astype(int).astype(str) + "_" + (np.floor(pc.lon / 30)).astype(int).astype(str)
ft = smf.ols("med ~ C(tile)", data=pc).fit()
R["tiles"] = dict(n_tiles=int(pc.tile.nunique()), R2=float(ft.rsquared), R2_adj=float(ft.rsquared_adj), F_p=float(ft.f_pvalue))
# ---- position within host crater (craters >= 30 km), crater fixed effects
big = h[(h.hostD >= 30) & h.rho.notna()].copy()
big["pos"] = pd.cut(big.rho, [-0.01, 0.25, 0.55, 0.85, 1.1], labels=["center", "inner_floor", "outer_floor_lower_wall", "wall_rim"]).astype(str)
fe = smf.ols("res ~ C(pos, Treatment('inner_floor')) + C(host)", data=big).fit()
R["position_effects_vs_inner_floor"] = {k.split("[T.")[1].rstrip("]"): [float(10 ** fe.params[k]), float(10 ** (fe.params[k] - 1.96 * fe.bse[k])), float(10 ** (fe.params[k] + 1.96 * fe.bse[k]))]
                                         for k in fe.params.index if k.startswith("C(pos")}
R["position_n"] = big.pos.value_counts().to_dict()
# ---- does depth depend on host age class or setting? (Herrick preservation if present)
json.dump(R, open(_REPO + "/manuscripts/insolation_thermal/depthctl/depth_structure_record.json", "w"), indent=1)
print(json.dumps(R, indent=1))
pc.to_csv(_REPO + "/manuscripts/insolation_thermal/depthctl/crater_median_depths.csv")
