#!/usr/bin/env python3
"""(v2, 5 Nov 2026: all cataloged hollows from hollows_canonical.csv; v1 used host_assignments.csv, hollows in host craters only.)
Orientation of elongated hollow fields NW/N of Caloris (5 Nov 2026). For each field with >= 5 hollows (host_assignments.csv,
crater traces dropped) at 35-75N, 90-170E, the long axis from the principal components of its hollow positions in a local
east-north frame (km). Compared with (a) the direction radial to Caloris (Orgel center 29.9N 161.4E), (b) the tangent to
candidate basin A (Michael's possiblebasins.gpkg; 143.19E 56.96N), (c) the tangent to b59 (Orgel; 119.26E 50.05N).
Statistic: mean |angle difference| (0-90 deg; 45 if unrelated) and the fraction within 30 deg (1/3 if unrelated), with a
permutation null that shuffles orientations among fields. Writes field_alignment.json."""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)]
h["lon360"] = h.lon % 360
def bearing(la1, lo1, la2, lo2):   # initial bearing from point 1 to point 2, degrees from north
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    return np.degrees(np.arctan2(np.sin(lo2 - lo1) * np.cos(la2), np.cos(la1) * np.sin(la2) - np.sin(la1) * np.cos(la2) * np.cos(lo2 - lo1))) % 180
rows = []
for fid, g in h[(h.lat.between(35, 75)) & (h.lon360.between(90, 170))].groupby("field_id"):
    if len(g) < 5: continue
    la0, lo0 = g.lat.mean(), g.lon360.mean()
    x = (g.lon360 - lo0).values * np.cos(np.radians(la0)) * np.pi / 180 * R_M; y = (g.lat - la0).values * np.pi / 180 * R_M
    C = np.cov(np.c_[x, y].T); w, v = np.linalg.eigh(C)
    if w[0] <= 0 or np.sqrt(w[1] / w[0]) < 1.5: continue           # keep elongated fields (axis ratio >= 1.5)
    ax = np.degrees(np.arctan2(v[0, 1], v[1, 1])) % 180              # long axis, degrees from north
    rows.append(dict(field=int(fid), n=len(g), lat=la0, lon=lo0, ratio=float(np.sqrt(w[1] / w[0])), axis=ax,
                     calrad=bearing(la0, lo0, 29.91, 161.43), Atan=(bearing(la0, lo0, 56.96, 143.19) + 90) % 180,
                     b59tan=(bearing(la0, lo0, 50.05, 119.26) + 90) % 180))
F = pd.DataFrame(rows)
def dang(a, b): d = np.abs(a - b) % 180; return np.minimum(d, 180 - d)
rng = np.random.default_rng(0); R = {"n_fields": int(len(F)), "median_axis_ratio": round(float(F.ratio.median()), 2)}
for ref in ["calrad", "Atan", "b59tan"]:
    dd = dang(F.axis.values, F[ref].values); null = np.array([dang(rng.permutation(F.axis.values), F[ref].values).mean() for _ in range(5000)])
    R[ref] = {"mean_abs_diff_deg": round(float(dd.mean()), 1), "frac_within_30": round(float((dd < 30).mean()), 2),
              "perm_p_aligned": round(float((null <= dd.mean()).mean()), 4), "null_mean": round(float(null.mean()), 1)}
F.round(2).to_csv(f"{HERE}/field_alignment_fields.csv", index=False)
json.dump(R, open(f"{HERE}/field_alignment.json", "w"), indent=1); print(json.dumps(R, indent=1))
