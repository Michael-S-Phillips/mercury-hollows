#!/usr/bin/env python3
"""Hosting and hollows against distance from the center of b59 (Orgel et al. 2020; 701 km, 50.1N 119.3E), in basin radii,
for craters at 35-70N (5 Nov 2026). Logistic fits with the paper's controls and distance bins; top host craters of the
province with their distance. Writes b59_profile.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
B = pd.read_csv(f"{ROOT}/peak_wall_test/orgel_basins.csv"); b = B[B.Basin_Name == "b59"].iloc[0]; c = xyz([b.y_coord], [b.x_coord % 360])[0]
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
h = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)]
f["r"] = np.arccos(np.clip(xyz(f.lat.values, f.lon.values % 360) @ c, -1, 1)) * R_M / (b.Diam_km / 2)
h["r"] = np.arccos(np.clip(xyz(h.lat.values, h.lon.values % 360) @ c, -1, 1)) * R_M / (b.Diam_km / 2)
edges = [0, 0.5, 1.0, 1.25, 1.5, 2.0, 3.0, 99]; labs = ["0-0.5", "0.5-1", "1-1.25", "1.25-1.5", "1.5-2", "2-3", ">3"]
d = f[f.lat.between(35, 70)].copy(); d["rb"] = pd.Categorical(pd.cut(d.r, edges, labels=labs, right=False), categories=labs[::-1])
m = smf.logit("h09 ~ C(rb) + logD + abslat + C1 + nsp + dT", d).fit(disp=0)
R = {"bins": {}}
for l in labs:
    g = d[d.rb == l]; t = f"C(rb)[T.{l}]"
    o = [round(float(np.exp(m.params[t])), 2), round(float(np.exp(m.conf_int().loc[t][0])), 2), round(float(np.exp(m.conf_int().loc[t][1])), 2)] if t in m.params else "reference"
    R["bins"][l] = {"craters": int(len(g)), "hosts": int(g.h09.sum()), "hollows": int(h[(h.lat.between(35, 70)) & (pd.cut(h.r, edges, labels=labs, right=False) == l)].shape[0]), "odds_ratio_vs_>3": o}
# top hosts in the province box with their distance from b59
hp = h[(h.lat.between(40, 70)) & ((h.lon % 360).between(90, 150)) & (h.rho <= 0.9)].copy()
top = hp.groupby("host_uid").agg(n=("lat", "size"), lat=("host_lat", "first"), lon=("host_lon", "first"), D=("host_D_km", "first")).sort_values("n", ascending=False).head(10)
top["r_b59"] = np.arccos(np.clip(xyz(top.lat.values, top.lon.values % 360) @ c, -1, 1)) * R_M / (b.Diam_km / 2)
R["top_hosts_province"] = top.round(2).reset_index().to_dict("records")
R["province_hollows_within_1.5R"] = int(((h.lat.between(40, 70)) & ((h.lon % 360).between(90, 150)) & (h.r < 1.5)).sum())
json.dump(R, open(f"{HERE}/b59_profile.json", "w"), indent=1, default=str); print(json.dumps(R, indent=1, default=str))
