#!/usr/bin/env python3
"""Do craters superposed on ghost craters host hollows when they dig through the ghost crater's volcanic fill? (5 Nov 2026)
Ghost craters: Kinczyk et al. (2020) code 6 (140 craters >= 40 km, entirely buried by lava; 70% on the northern smooth plains).
Superposed craters: craters >= 5 km of crater_features.csv whose centre lies inside a ghost crater at least twice their size
(smallest such ghost crater), rho = distance / ghost radius.
Fill thickness over the ghost floor, lower bound: the fresh rim-to-floor depth d = 0.241 D^0.486 (90th percentile of isolated
craters, Section S7 of the supplement; plains surface level with the buried rim crest), capped at 3 km (Caloris interior plains)
for ghost basins; + 0.5 and + 1 km if the plains surface lies above the rim crest. Interior only (rho < 0.6) in the main test;
on the wall (0.6-1) the fill thins toward the rim, taken here as linear to zero at rho = 1.
Reach: excavation depth 0.1 D_tc (D_tc = 12^0.15 D^0.85 for D > 12 km, D/1.2 below), as in the Caloris test, and z_own =
max(excavation, uplift 0.086 D^1.03) as in the containment graph. Hosting = h09. Writes ghost_superposed.json, ghost_superposed.csv."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
from scipy.spatial import cKDTree
from scipy import stats
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); R_M = 2439.4
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
k = pd.read_csv(f"{HERE}/data/crater_degradation_database_MJK_2019.txt", sep="\t")
g = k[k.Class == 6].rename(columns={"x_coord": "glon", "y_coord": "glat", "Diameter_km": "gD"}).reset_index(drop=True)
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
D = f.D.values
f["d_ex"] = 0.1 * np.where(D > 12, 12 ** 0.15 * D ** 0.85, D / 1.2)
f["z_own"] = np.maximum(f.d_ex, np.where(D > 12, 0.086 * D ** 1.03, 0.0))
U = xyz(f.lat.values, f.lon.values % 360); G = xyz(g.glat.values, g.glon.values % 360)
gi = np.full(len(f), -1); rho = np.full(len(f), np.nan)
for j in np.argsort(-g.gD.values):                       # larger first; smaller containing ghost overwrites
    R = g.gD.values[j] / 2
    hits = np.asarray(cKDTree(U).query_ball_point(G[j], 2 * np.sin(R / R_M / 2)) if False else [])
for j in np.argsort(-g.gD.values):
    R = g.gD.values[j] / 2
    a = np.arccos(np.clip(U @ G[j], -1, 1)) * R_M
    ok = (a < R) & (D <= g.gD.values[j] / 2)
    gi[ok] = j; rho[ok] = a[ok] / R
s = f[gi >= 0].copy(); s["g"] = gi[gi >= 0]; s["rho"] = rho[gi >= 0]
s["gD"] = g.gD.values[s.g]; s["glat"] = g.glat.values[s.g]
s["t0"] = np.minimum(0.241 * s.gD ** 0.486, 3.0)
taper = np.where(s.rho < 0.6, 1.0, (1 - s.rho) / 0.4)
R = {"ghost_n": int(len(g)), "superposed_n": int(len(s)), "superposed_hosts": int(s.h09.sum()),
     "ghosts_with_superposed": int(s.g.nunique()), "interior_n": int((s.rho < 0.6).sum()), "interior_hosts": int(s[s.rho < 0.6].h09.sum())}
rows = []
for dz in [0.0, 0.5, 1.0]:
    for reachvar in ["d_ex", "z_own"]:
        for zone, sub in [("interior", s.rho < 0.6), ("all_taper", s.rho < 1.0)]:
            x = s[sub].copy(); t = (x.t0 + dz) * (taper[sub.values] if zone == "all_taper" else 1.0)
            x["pierce"] = (x[reachvar] > t).astype(int)
            a, b = x[x.pierce == 1], x[x.pierce == 0]
            tab = [[int(a.h09.sum()), int(len(a) - a.h09.sum())], [int(b.h09.sum()), int(len(b) - b.h09.sum())]]
            p = float(stats.fisher_exact(tab)[1]) if len(a) and len(b) else None
            rows.append(dict(extra_fill_km=dz, reach=reachvar, zone=zone, pierce_n=int(len(a)), pierce_hosts=int(a.h09.sum()),
                             no_pierce_n=int(len(b)), no_pierce_hosts=int(b.h09.sum()), fisher_p=p,
                             pierce_D_range=[round(float(a.D.min()), 1), round(float(a.D.max()), 1)] if len(a) else None,
                             no_pierce_D_max=round(float(b.D.max()), 1) if len(b) else None))
R["tests"] = rows
# size-matched comparison: hosting by size among (i) ghost-superposed, (ii) northern smooth plains outside ghosts, (iii) elsewhere
f["grp"] = np.where(gi >= 0, "ghost", np.where(f.nsp == 1, "nsp_other", "elsewhere"))
bins = [5, 10, 20, 40, 80, 2000]; f["Dbin"] = pd.cut(f.D, bins, right=False)
R["by_size"] = {str(b): {gname: [int(len(x)), int(x.h09.sum())] for gname, x in gg.groupby("grp")} for b, gg in f.groupby("Dbin", observed=True)}
# logistic within ghost-superposed craters: piercing (lower-bound fill, excavation) at fixed size
x = s[s.rho < 1.0].copy(); x["pierce"] = (x.d_ex > x.t0 * np.where(x.rho < 0.6, 1, (1 - x.rho) / 0.4)).astype(int)
try:
    m = smf.logit("h09 ~ pierce + logD", x).fit(disp=0); b_, ci = m.params["pierce"], m.conf_int().loc["pierce"]
    R["logit_pierce_given_logD"] = [round(float(np.exp(b_)), 2), round(float(np.exp(ci[0])), 2), round(float(np.exp(ci[1])), 2), float(m.pvalues["pierce"])]
except Exception as e:
    R["logit_pierce_given_logD"] = str(e)
# hosting vs the same-size, non-ghost NSP craters (logit with group dummies)
y = f[f.grp.isin(["ghost", "nsp_other"])].copy(); y["ghost"] = (y.grp == "ghost").astype(int)
m = smf.logit("h09 ~ ghost + logD + abslat + C1", y).fit(disp=0); b_, ci = m.params["ghost"], m.conf_int().loc["ghost"]
R["ghost_vs_nsp_other"] = [round(float(np.exp(b_)), 2), round(float(np.exp(ci[0])), 2), round(float(np.exp(ci[1])), 2), float(m.pvalues["ghost"])]
host_list = s[s.h09 == 1][["id", "lat", "lon", "D", "d_ex", "z_own", "gD", "t0", "rho", "nsp", "n09"]].round(3)
R["hosts_detail"] = host_list.to_dict(orient="records")
s.to_csv(f"{HERE}/ghost_superposed.csv", index=False)
json.dump(R, open(f"{HERE}/ghost_superposed.json", "w"), indent=1, default=str)
print(json.dumps({k_: v for k_, v in R.items() if k_ != "tests"}, indent=1, default=str))
print(pd.DataFrame(rows).to_string())
