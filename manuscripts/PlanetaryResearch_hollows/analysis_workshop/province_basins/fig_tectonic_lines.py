#!/usr/bin/env python3
"""Figure (5 Nov 2026): mapped tectonic landforms of Klimczak et al. (2025 data) and hollows in 38-62N, 98-145E, on the
color basemap of fig_hrp_mask_v2.py, with b59, the NW-SE band (band_fit.json) and radial directions to the Caloris center.
Writes fig_tectonic_lines.png."""
import os, sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
sys.path.insert(0, HERE); from gpkg_wkb import load
B = np.load(f"{HERE}/basemap_cache.npz"); pan, col, ext = B["pan"], B["col"], B["extent"]
col = np.repeat(np.repeat(col, 2, axis=1), 2, axis=2)[:, :pan.shape[0], :pan.shape[1]]
def stretch(a, p=(1.5, 98.5)):
    lo_, hi_ = np.nanpercentile(a, p); return np.clip((a - lo_) / (hi_ - lo_), 0, 1)
rgb = np.stack([stretch(col[2]), stretch(col[1]), stretch(col[0])]); I = rgb.mean(0); valid = np.isfinite(I) & (pan > 0)
P = (pan - pan[valid].mean()) / pan[valid].std() * np.nanstd(I[valid]) + np.nanmean(I[valid])
out = np.clip(rgb * (P / np.where(I > 0.02, I, 0.02))[None], 0, 1); out = np.nan_to_num(out) ** 0.9; out = np.moveaxis(out, 0, -1); out[~valid] = 0.15
def circle(clat, clon, Rkm, n=720):
    d = Rkm / R_M; b = np.radians(np.linspace(0, 360, n)); p1, l1 = np.radians(clat), np.radians(clon)
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(b))
    l2 = l1 + np.arctan2(np.sin(b) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2)); return np.degrees(p2), np.degrees(l2) % 360
def gc(la1, lo1, la2, lo2, n=200):
    a = np.radians([la1, lo1, la2, lo2]); p = np.array([np.cos(a[0])*np.cos(a[1]), np.cos(a[0])*np.sin(a[1]), np.sin(a[0])]); q = np.array([np.cos(a[2])*np.cos(a[3]), np.cos(a[2])*np.sin(a[3]), np.sin(a[2])])
    w = np.arccos(p @ q); t = np.linspace(0, 1, n)[:, None]; v = (np.sin((1 - t) * w) * p + np.sin(t * w) * q) / np.sin(w)
    return np.degrees(np.arcsin(v[:, 2])), np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360
fig, ax = plt.subplots(figsize=(12, 6.6), dpi=200)
ax.imshow(out, extent=list(ext), aspect=1 / np.cos(np.radians(50)), interpolation="lanczos")
for k in range(-6, 7):                                          # radial directions to Caloris, every 2 deg of azimuth around 315
    az = np.radians(315 + 2.5 * k); d = 2600 / R_M; p1, l1 = np.radians(29.91), np.radians(161.43)
    p2 = np.arcsin(np.sin(p1)*np.cos(d) + np.cos(p1)*np.sin(d)*np.cos(az)); l2 = l1 + np.arctan2(np.sin(az)*np.sin(d)*np.cos(p1), np.cos(d) - np.sin(p1)*np.sin(p2))
    y, x = gc(29.91, 161.43, np.degrees(p2), np.degrees(l2) % 360); ax.plot(x, y, color="white", lw=0.4, alpha=0.45)
cols = {"Contractional Landform": "#ff9f1c", "Extensional Landform": "#2ec4b6", "Trough": "#ff4fd8"}
for fid, ft, L, parts in load(f"{ROOT}/GIS/vector/klimczak_tectonic_map_lines.gpkg"):
    for a in parts:
        lo = a[:, 0] % 360
        if ((lo > 96) & (lo < 147) & (a[:, 1] > 36) & (a[:, 1] < 64)).any(): ax.plot(lo, a[:, 1], color=cols[ft], lw=1.0)
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)]; h = h.assign(lon=h.lon % 360)
ax.scatter(h.lon, h.lat, s=2.5, c="#fff200", edgecolors="none", zorder=5)
y, x = circle(50.05, 119.26, 350.3); ax.plot(x, y, color="#ff5a5a", lw=1.4, ls="--")
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_xlabel("East longitude (°)"); ax.set_ylabel("Latitude (°N)")
hd = [Line2D([], [], color=cols["Contractional Landform"], lw=1.2, label="contractional landform (Klimczak et al. 2025 map)"),
      Line2D([], [], color=cols["Extensional Landform"], lw=1.2, label="extensional landform"),
      Line2D([], [], color="white", lw=0.6, alpha=0.6, label="radial to the Caloris center"),
      Line2D([], [], color="#ff5a5a", lw=1.4, ls="--", label="b59 rim (Orgel et al. 2020)"),
      Line2D([], [], marker="o", ls="", mfc="#fff200", mec="none", ms=3, label="cataloged hollow")]
ax.legend(handles=hd, loc="lower right", fontsize=7.5, framealpha=0.85, facecolor="#1a1a19", edgecolor="#6b6b67", labelcolor="white")
ax.set_title("Mapped tectonic landforms (no troughs mapped in this frame) and hollows NW of Caloris", fontsize=9)
fig.savefig(f"{HERE}/fig_tectonic_lines.png", bbox_inches="tight"); print("ok")
