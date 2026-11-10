#!/usr/bin/env python3
"""g20 (9 Nov 2026): two quick checks on the excess regions. (1) Region 1: how many hosts owe their status only to hollows < 1 km, and
the detection pixel scale of those hollows. (2) Observed / predicted hosts (joint model) in rings of distance from the Caloris center,
outside the basin, and by azimuth in the 1,500-2,200 km ring. Writes region_checks.json. Not for the paper."""
import os, json, numpy as np, pandas as pd
HE = os.path.expanduser("~/mnt/hollows-emissivity"); RR = f"{HE}/code/rerun/out_v3/rerun"; HERE = os.path.dirname(os.path.abspath(__file__)); RM = 2439.7
C = pd.read_csv(f"{HERE}/craters_P.csv", dtype={"id": str})
H = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool)]
ha = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/host_assignments.csv")[["hollow_id", "host_uid", "rho"]]; H = H.merge(ha, on="hollow_id")
H["host_uid"] = H.host_uid.astype("Int64").astype(str); H = H[H.rho <= 0.9]
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
def dist(la, lo, la0, lo0): return RM * np.arccos(np.clip(xyz(la, lo) @ xyz(np.array([la0]), np.array([lo0]))[0], -1, 1))
out = {}
for nm, la0, lo0 in (("1_province", 65, 125), ("2_abu_nuwas_sinan", 13, 329), ("3_tyagaraja", 7, 203), ("all", None, None)):
    s = C if la0 is None else C[dist(C.lat.values, C.lon.values % 360, la0, lo0) <= 500]
    hosts = s[s.h09 == 1].id; g = H[H.host_uid.isin(hosts)].groupby("host_uid").agg(maxD=("D_eq_km", "max"), n=("D_eq_km", "size"), px=("px_best", "median"))
    out[nm] = dict(hosts=int(len(hosts)), hosts_only_lt1km=int((g.maxD < 1).sum()), hosts_only_lt0p5km=int((g.maxD < 0.5).sum()),
                   median_px_of_hosts_only_lt1km=(round(float(g[g.maxD < 1].px.median()), 1) if (g.maxD < 1).any() else None))
d = dist(C.lat.values, C.lon.values % 360, 29.91, 161.43); C["dcal"] = d
az = np.degrees(np.arctan2(np.sin(np.radians(C.lon % 360 - 161.43)) * np.cos(np.radians(C.lat)),
                           np.cos(np.radians(29.91)) * np.sin(np.radians(C.lat)) - np.sin(np.radians(29.91)) * np.cos(np.radians(C.lat)) * np.cos(np.radians(C.lon % 360 - 161.43)))) % 360
C["az"] = az; rings = []
for a, b in ((0, 778), (778, 1200), (1200, 1600), (1600, 2000), (2000, 2400), (2400, 2800), (2800, 3500), (3500, 8000)):
    t = C[(C.dcal >= a) & (C.dcal < b)]; rings.append(dict(km=[a, b], n=int(len(t)), hosts=int(t.h09.sum()), predicted=round(float(t.P.sum()), 1), ratio=round(float(t.h09.sum() / t.P.sum()), 2)))
out["caloris_rings"] = rings; sec = []
t = C[(C.dcal >= 1500) & (C.dcal < 2200)]
for a in range(0, 360, 45):
    u = t[(t.az >= a) & (t.az < a + 45)]; sec.append(dict(az=[a, a + 45], n=int(len(u)), hosts=int(u.h09.sum()), predicted=round(float(u.P.sum()), 1)))
out["ring_1500_2200_by_azimuth"] = sec
json.dump(out, open(f"{HERE}/region_checks.json", "w"), indent=1); print(json.dumps(out, indent=0)[:3000])
