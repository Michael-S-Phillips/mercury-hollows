#!/usr/bin/env python3
"""Assets for the ridge-tracing page (5 Nov 2026). Frame 38-62N, 98-145E in plate carree with x scaled by cos(50 deg),
so image x = (lon - 98) / 47 * W and y = (62 - lat) / 24 * H exactly. Writes basemap_color.jpg (the pan-sharpened MDR
color of fig_hrp_mask_v2.py), basemap_gray.jpg (BDR only), and overlays.json (hollows, Klimczak et al. 2025 tectonic
lines, b59 rim, candidate basins A and B, radial lines from the Caloris center)."""
import os, sys, json, sqlite3, struct, numpy as np, pandas as pd
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); PB = os.path.dirname(HERE); ROOT = os.path.abspath(f"{PB}/../../../.."); R_M = 2439.4
sys.path.insert(0, PB); from gpkg_wkb import load
LO0, LO1, LA0, LA1 = 98.0, 145.0, 38.0, 62.0; W = 2400; H = int(round(W / ((LO1 - LO0) * np.cos(np.radians(50)) / (LA1 - LA0))))
B = np.load(f"{PB}/basemap_cache.npz"); pan, col = B["pan"], B["col"]
col = np.repeat(np.repeat(col, 2, axis=1), 2, axis=2)[:, :pan.shape[0], :pan.shape[1]]
def stretch(a, p=(1.5, 98.5)):
    lo_, hi_ = np.nanpercentile(a, p); return np.clip((a - lo_) / (hi_ - lo_), 0, 1)
rgb = np.stack([stretch(col[2]), stretch(col[1]), stretch(col[0])]); I = rgb.mean(0); valid = np.isfinite(I) & (pan > 0)
P = (pan - pan[valid].mean()) / pan[valid].std() * np.nanstd(I[valid]) + np.nanmean(I[valid])
out = np.clip(rgb * (P / np.where(I > 0.02, I, 0.02))[None], 0, 1); out = np.nan_to_num(out) ** 0.9; out = np.moveaxis(out, 0, -1); out[~valid] = 0.15
Image.fromarray((out * 255).astype(np.uint8)).resize((W, H), Image.LANCZOS).save(f"{HERE}/basemap_color.jpg", quality=86, optimize=True)
g = stretch(np.where(pan > 0, pan, np.nan), (0.5, 99.5)); g = np.nan_to_num(g, nan=0.15)
Image.fromarray((g * 255).astype(np.uint8)).resize((W, H), Image.LANCZOS).save(f"{HERE}/basemap_gray.jpg", quality=86, optimize=True)
inside = lambda la, lo: (la > LA0 - 1) & (la < LA1 + 1) & (lo > LO0 - 1.5) & (lo < LO1 + 1.5)
O = {"frame": {"lon": [LO0, LO1], "lat": [LA0, LA1], "image_px": [W, H], "x_scale_lat_deg": 50},
     "caloris_center": [29.91, 161.43], "b59": {"center": [50.05, 119.26], "radius_km": 350.3}}
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)].copy(); h["lon"] = h.lon % 360
h = h[inside(h.lat, h.lon)]; O["hollows"] = [[round(a, 3), round(b, 3)] for a, b in zip(h.lon, h.lat)]
lines = []
for fid, ft, L, parts in load(f"{ROOT}/GIS/vector/klimczak_tectonic_map_lines.gpkg"):
    for a in parts:
        lo, la = a[:, 0] % 360, a[:, 1]
        if not inside(la, lo).any(): continue
        keep = [0]
        for i in range(1, len(lo)):
            j = keep[-1]
            if np.hypot((lo[i] - lo[j]) * np.cos(np.radians(la[i])), la[i] - la[j]) * np.pi / 180 * R_M > 1.5 or i == len(lo) - 1: keep.append(i)
        lines.append({"type": ft, "pts": [[round(float(lo[i]), 3), round(float(la[i]), 3)] for i in keep]})
O["tectonic"] = lines
def circle(clat, clon, Rkm, n=360):
    d = Rkm / R_M; b = np.radians(np.linspace(0, 360, n)); p1, l1 = np.radians(clat), np.radians(clon)
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(b))
    l2 = l1 + np.arctan2(np.sin(b) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2)); return np.degrees(p2), np.degrees(l2) % 360
y, x = circle(50.05, 119.26, 350.3); O["b59"]["pts"] = [[round(float(a), 3), round(float(b), 3)] for a, b in zip(x, y)]
c = sqlite3.connect(f"{ROOT}/GIS/vector/possiblebasins.gpkg"); cands = []
for fid, gb in c.execute("select fid, geom from possiblebasins"):
    env = (gb[3] >> 1) & 7; w = gb[8 + {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]:]
    bo = "<" if w[0] == 1 else ">"; n = struct.unpack(bo + "I", w[9:13])[0]
    xy = np.frombuffer(w[13:13 + 16 * n], dtype=bo + "f8").reshape(n, 2); rho = np.hypot(xy[:, 0], xy[:, 1]); Rs = 2440530.0
    lat = 90 - np.degrees(2 * np.arctan(rho / (2 * Rs))); lon = np.degrees(np.arctan2(xy[:, 0], -xy[:, 1])) % 360
    clat, clon = float(np.mean(lat)), float(np.mean(lon))
    for name, (tl, tn) in {"A": (56.96, 143.19), "B": (52.45, 125.68)}.items():
        if abs(clat - tl) < 2 and abs(clon - tn) < 3 and name not in [x["name"] for x in cands]: cands.append({"name": name, "pts": [[round(float(a), 3), round(float(b), 3)] for a, b in zip(lon, lat)]})
O["candidates"] = cands
rad = []
for az in np.arange(280, 352.5, 2.5):
    pts = []
    for dkm in np.arange(500, 3200, 15):
        d = dkm / R_M; p1, l1, b = np.radians(29.91), np.radians(161.43), np.radians(az)
        p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(b)); l2 = l1 + np.arctan2(np.sin(b) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2))
        la_, lo_ = np.degrees(p2), np.degrees(l2) % 360
        if inside(la_, lo_): pts.append([round(float(lo_), 3), round(float(la_), 3)])
    if len(pts) > 1: rad.append({"az": float(az), "pts": pts})
O["radials"] = rad
json.dump(O, open(f"{HERE}/overlays.json", "w"), separators=(",", ":"))
print(W, H, len(O["hollows"]), len(lines), [x["name"] for x in cands], len(rad), {f: os.path.getsize(f"{HERE}/{f}") for f in ("basemap_color.jpg", "basemap_gray.jpg", "overlays.json")})
