#!/usr/bin/env python3
"""Hollows and the mapped tectonic landforms NW of Caloris (5 Nov 2026). Lines: the global tectonic map of Klimczak et al.
(2025 data; GIS/vector/klimczak_tectonic_map_lines.gpkg: Contractional Landform, Extensional Landform, Trough). Region:
38-62N, 98-145E (the basemap frame: b59, the NW-SE band, candidate A). (1) Orientation of the mapped segments relative to
the direction to the Caloris center (Orgel et al. 2020, 29.91N 161.43E), length-weighted, by type. (2) Distance from
each hollow field (non-trace hollows in hollows_canonical.csv, grouped by field_id) to the nearest mapped landform of each
type, against the fraction of the region's area within the same distance (area-weighted grid), binomial test.
Writes tectonic_lines.json and fig_tectonic_lines.png."""
import os, sys, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import binomtest
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
sys.path.insert(0, HERE); from gpkg_wkb import load
LA0, LO0, LO_MIN, LO_MAX, LA_MIN, LA_MAX = 29.91, 161.43, 98, 145, 38, 62
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
def bearing(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    return np.degrees(np.arctan2(np.sin(lo2-lo1)*np.cos(la2), np.cos(la1)*np.sin(la2)-np.sin(la1)*np.cos(la2)*np.cos(lo2-lo1))) % 360
F = load(f"{ROOT}/GIS/vector/klimczak_tectonic_map_lines.gpkg")
TYPES = ["Contractional Landform", "Extensional Landform", "Trough"]
verts = {t: [] for t in TYPES}; segs = []
for fid, ft, L, parts in F:
    for a in parts:
        lo, la = a[:, 0] % 360, a[:, 1]
        if not ((lo > LO_MIN - 2) & (lo < LO_MAX + 2) & (la > LA_MIN - 2) & (la < LA_MAX + 2)).any(): continue
        P = xyz(la, lo); d = np.arccos(np.clip((P[:-1] * P[1:]).sum(1), -1, 1)) * R_M
        for i in range(len(P) - 1):                                 # densify to <= 0.5 km
            k = max(1, int(np.ceil(d[i] / 0.5))); w = np.linspace(0, 1, k, endpoint=False)[:, None]
            verts[ft].append(P[i] * (1 - w) + P[i + 1] * w)
            mla, mlo = (la[i] + la[i + 1]) / 2, (lo[i] + lo[i + 1]) / 2
            if LA_MIN <= mla <= LA_MAX and LO_MIN <= mlo <= LO_MAX and d[i] > 0:
                az = bearing(la[i], lo[i], la[i + 1], lo[i + 1]) % 180; rad = bearing(mla, mlo, LA0, LO0) % 180
                segs.append((ft, fid, d[i], min(abs(az - rad), 180 - abs(az - rad)), mla, mlo))
        verts[ft].append(P[-1:])
S = pd.DataFrame(segs, columns=["type", "fid", "len_km", "ang_radial", "lat", "lon"])
R = {"region": "38-62N, 98-145E", "caloris_center": [LA0, LO0], "landforms": {}}
for t in TYPES:
    s = S[S.type == t]; w = s.len_km.sum()
    R["landforms"][t] = {"features": int(s.fid.nunique()), "length_km": round(float(w)),
        "frac_len_within30_of_radial": round(float(s.len_km[s.ang_radial <= 30].sum() / w), 3) if w else None,
        "frac_len_within30_of_concentric": round(float(s.len_km[s.ang_radial >= 60].sum() / w), 3) if w else None,
        "uniform_expectation": 0.333}
trees = {t: cKDTree(np.vstack(verts[t])) for t in TYPES if len(verts[t])}
trees["any"] = cKDTree(np.vstack([np.vstack(verts[t]) for t in TYPES if len(verts[t])]))
def dist_km(tree, P): return 2 * np.arcsin(np.clip(tree.query(P)[0] / 2, 0, 1)) * R_M
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)].copy(); h["lon"] = h.lon % 360
h = h[h.lat.between(LA_MIN, LA_MAX) & h.lon.between(LO_MIN, LO_MAX)]
ha = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", usecols=["hollow_id"])
h["in_host"] = h.hollow_id.isin(ha.hollow_id)
fl = h.groupby("field_id").agg(lat=("lat", "mean"), lon=("lon", "mean"), n=("hollow_id", "size"), in_host=("in_host", "mean")).reset_index()
FP = xyz(fl.lat.values, fl.lon.values); HP = xyz(h.lat.values, h.lon.values)
g_la, g_lo = np.meshgrid(np.arange(LA_MIN + 0.025, LA_MAX, 0.05), np.arange(LO_MIN + 0.025, LO_MAX, 0.05), indexing="ij")
GP = xyz(g_la.ravel(), g_lo.ravel()); gw = np.cos(np.radians(g_la.ravel())); gw /= gw.sum()
R["hollows"] = int(len(h)); R["fields"] = int(len(fl)); R["fields_mostly_outside_hosts"] = int((fl.in_host < 0.5).sum()); R["proximity"] = {}
for t in [x for x in TYPES if x in trees] + ["any"]:
    df, dh, dg = dist_km(trees[t], FP), dist_km(trees[t], HP), dist_km(trees[t], GP)
    out = {"median_field_distance_km": round(float(np.median(df)), 1), "median_area_distance_km": round(float(np.median(dg)), 1)}
    for X in (10, 25):
        a = float(gw[dg < X].sum()); k = int((df < X).sum()); ko = int(((df < X) & (fl.in_host.values < 0.5)).sum()); no = int((fl.in_host.values < 0.5).sum())
        out[f"within_{X}km"] = {"area_fraction": round(a, 3), "fields": f"{k}/{len(fl)}", "field_fraction": round(k / len(fl), 3),
            "p_more_than_area": round(binomtest(k, len(fl), a, alternative="greater").pvalue, 4),
            "hollow_fraction": round(float((dh < X).mean()), 3),
            "fields_outside_hosts": f"{ko}/{no}", "p_outside_hosts": round(binomtest(ko, no, a, alternative="greater").pvalue, 4) if no else None}
    R["proximity"][t] = out
json.dump(R, open(f"{HERE}/tectonic_lines.json", "w"), indent=1); print(json.dumps(R, indent=1))
S.to_csv(f"{HERE}/tectonic_segments_region.csv", index=False); fl.to_csv(f"{HERE}/tectonic_fields.csv", index=False)
