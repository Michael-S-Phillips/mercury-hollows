#!/usr/bin/env python3
"""Match the Kinczyk et al. (2020) degradation catalog (craters >= 40 km, classes 1-5 plus codes 0 and 6) to the
paper's crater table (summary_model/crater_features.csv, Herrick et al. 2018 craters >= 5 km). 5 Nov 2026.
Match: nearest crater_features crater within 0.3 of the Kinczyk radius... specifically centre distance <= 0.3 * D_K / 2
... widened below; diameter ratio 0.67-1.5; one-to-one (closest pairs first). Writes kinczyk_matched.csv and match_summary.json."""
import os, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
R_M = 2439.4
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
k = pd.read_csv(f"{HERE}/data/crater_degradation_database_MJK_2019.txt", sep="\t")
k = k.rename(columns={"x_coord": "klon", "y_coord": "klat", "Diameter_km": "kD", "Class": "kclass"})
fc = f[f.D >= 25].reset_index(drop=True)
t = cKDTree(xyz(fc.lat.values, fc.lon.values % 360))
P = xyz(k.klat.values, k.klon.values % 360)
cand = []
for i, (p, D) in enumerate(zip(P, k.kD.values)):
    r = 0.5 * D / R_M  # chord ~ angle for small angles; search one radius
    for j in t.query_ball_point(p, r):
        ang = np.arccos(np.clip(np.dot(p, t.data[j]), -1, 1)) * R_M
        ratio = fc.D.values[j] / D
        if 0.67 <= ratio <= 1.5:
            cand.append((ang / (D / 2) + abs(np.log(ratio)), i, j, ang, ratio))
cand.sort(); used_k, used_f, rows = set(), set(), []
for s, i, j, ang, ratio in cand:
    if i in used_k or j in used_f: continue
    used_k.add(i); used_f.add(j); rows.append((i, j, ang, ratio))
m = pd.DataFrame(rows, columns=["ki", "fj", "sep_km", "D_ratio"])
out = k.copy(); out["id"] = None; out["sep_km"] = np.nan; out["D_ratio"] = np.nan
out.loc[m.ki, "id"] = fc.id.values[m.fj]; out.loc[m.ki, "sep_km"] = m.sep_km.values; out.loc[m.ki, "D_ratio"] = m.D_ratio.values
out = out.merge(f, on="id", how="left")
out.to_csv(f"{HERE}/kinczyk_matched.csv", index=False)
mm = out[out.id.notna()]
S = {"kinczyk_n": int(len(k)), "kinczyk_class_counts": {int(a): int(b) for a, b in k.kclass.value_counts().sort_index().items()},
     "matched": int(len(mm)), "sep_over_R_median": round(float((mm.sep_km / (mm.kD / 2)).median()), 3),
     "sep_over_R_p95": round(float((mm.sep_km / (mm.kD / 2)).quantile(0.95)), 3),
     "D_ratio_median": round(float(mm.D_ratio.median()), 3), "D_ratio_p5_p95": [round(float(x), 3) for x in mm.D_ratio.quantile([0.05, 0.95])],
     "unmatched_by_class": {int(a): int(b) for a, b in out[out.id.isna()].kclass.value_counts().sort_index().items()},
     "features_D40_n": int((f.D >= 40).sum()), "features_D40_matched": int(mm[mm.D >= 40].shape[0]),
     "crosstab_kclass_vs_herrick_age": pd.crosstab(mm.kclass, mm.age).to_dict(),
     "hosts_h09_by_kclass": {int(a): [int(len(g)), int(g.h09.sum())] for a, g in mm.groupby("kclass")}}
json.dump(S, open(f"{HERE}/match_summary.json", "w"), indent=1); print(json.dumps(S, indent=1))
