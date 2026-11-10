#!/usr/bin/env python3
"""Hollows (host_assignments.csv) that lie inside ghost craters (Kinczyk code 6), with their host crater (h09 rule, rho <= 0.9)
and whether that host's excavation reaches below the ghost-crater fill (lower bound 0.241 D^0.486, capped at 3 km). 5 Nov 2026.
Writes ghost_hollows.json."""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
k = pd.read_csv(f"{HERE}/data/crater_degradation_database_MJK_2019.txt", sep="\t"); g = k[k.Class == 6].reset_index(drop=True)
h = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", low_memory=False)
h = h[~h.is_crater_trace.astype(bool)]
H = xyz(h.lat.values, h.lon.values % 360); G = xyz(g.y_coord.values, g.x_coord.values % 360)
gi = np.full(len(h), -1); grho = np.full(len(h), np.nan)
for j in np.argsort(-g.Diameter_km.values):
    a = np.arccos(np.clip(H @ G[j], -1, 1)) * R_M; ok = a < g.Diameter_km.values[j] / 2; gi[ok] = j; grho[ok] = a[ok] / (g.Diameter_km.values[j] / 2)
x = h[gi >= 0].copy(); x["g"] = gi[gi >= 0]; x["grho"] = grho[gi >= 0]; x["gD"] = g.Diameter_km.values[x.g]
x["t0"] = np.minimum(0.241 * x.gD ** 0.486, 3.0)
hd = x.host_D_km.where(x.rho <= 0.9)
x["host_in_ghost"] = (hd.notna() & (hd <= x.gD / 2))
x["host_d_ex"] = 0.1 * np.where(hd > 12, 12 ** 0.15 * hd ** 0.85, hd / 1.2)
R = {"hollows_total": int(len(h)), "in_ghosts": int(len(x)), "ghosts_with_hollows": int(x.g.nunique()),
     "fields": int(x.field_id.nunique()),
     "by_ghost": [dict(ghost_D=round(float(gg.gD.iat[0]), 1), lat=round(float(g.y_coord[j]), 2), lon=round(float(g.x_coord[j]), 2),
                       fill_lb_km=round(float(gg.t0.iat[0]), 2), hollows=int(len(gg)),
                       hosts=sorted({(round(float(a), 1), round(float(b), 2)) for a, b in zip(gg.host_D_km.where(gg.rho <= 0.9).fillna(-1), gg.host_d_ex.fillna(-1))}),
                       median_grho=round(float(gg.grho.median()), 2))
                  for j, gg in x.groupby("g")]}
json.dump(R, open(f"{HERE}/ghost_hollows.json", "w"), indent=1); print(json.dumps(R, indent=1))
