#!/usr/bin/env python3
"""How unusual is b59? For each Orgel et al. (2020) basin >= 300 km other than Caloris, the odds ratio of hosting for
craters within 1.5 basin radii of its center against craters 1.5-4 radii away, with the paper's controls (5 Nov 2026).
Writes basin_rank.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
B = pd.read_csv(f"{ROOT}/peak_wall_test/orgel_basins.csv"); B = B[(B.Basin_Name != "Caloris") & (B.Diam_km >= 300)].reset_index(drop=True)
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str}); U = xyz(f.lat.values, f.lon.values % 360)
rows = []
for _, b in B.iterrows():
    r = np.arccos(np.clip(U @ xyz([b.y_coord], [b.x_coord % 360])[0], -1, 1)) * R_M / (b.Diam_km / 2)
    d = f[r < 4].copy(); d["near"] = (r[r < 4] < 1.5).astype(int)
    if d.h09.sum() < 5 or d.near.sum() < 20: continue
    try:
        m = smf.logit("h09 ~ near + logD + abslat + C1 + nsp + dT", d).fit(disp=0)
        rows.append(dict(name=b.Basin_Name, D=round(b.Diam_km), lat=round(b.y_coord, 1), lon=round(b.x_coord % 360, 1), conf=b.Confidence,
                         n_near=int(d.near.sum()), hosts_near=int(d[d.near == 1].h09.sum()), OR=round(float(np.exp(m.params["near"])), 2),
                         lo=round(float(np.exp(m.conf_int().loc["near"][0])), 2), hi=round(float(np.exp(m.conf_int().loc["near"][1])), 2), p=float(m.pvalues["near"])))
    except Exception as e:
        pass
T = pd.DataFrame(rows).sort_values("OR", ascending=False).reset_index(drop=True)
R = {"n_basins_tested": int(len(T)), "b59_rank": int(T.index[T.name == "b59"][0]) + 1, "top10": T.head(10).to_dict("records"),
     "median_OR": round(float(T.OR.median()), 2), "n_OR_gt_2_p_lt_0.01": int(((T.OR > 2) & (T.p < 0.01)).sum())}
json.dump(R, open(f"{HERE}/basin_rank.json", "w"), indent=1); print(json.dumps(R, indent=1))
