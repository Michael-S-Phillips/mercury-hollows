#!/usr/bin/env python3
"""Figure (5 Nov 2026): hollow fields and the higher-albedo smooth plains inside b59, on a color basemap. Basemap: MDIS MDR
color (R 1000, G 750, B 480 nm) sharpened with the USGS MDIS BDR 166 m mosaic (Brovey-style ratio), as Fassett et al.
(2009, Fig. 9) merged color with a NAC mosaic. Overlays: the automated bright-plains mask of hrp_edge_test.py (white
outline), b59 (Orgel et al. 2020) and the 1.2-radius test domain, and hollow fields (hollows_canonical.csv) colored by
their position relative to the mask and sized by number of hollows. Fassett arrows dropped: read off their conic projection by eye, they were too uncertain to place. Writes fig_hrp_mask_v2.png."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__)); R_M = 2439.4
B = np.load(f"{HERE}/basemap_cache.npz"); pan, col, ext = B["pan"], B["col"], B["extent"]
col = np.repeat(np.repeat(col, 2, axis=1), 2, axis=2)[:, :pan.shape[0], :pan.shape[1]]
def stretch(a, p=(1.5, 98.5)):
    lo_, hi_ = np.nanpercentile(a, p); return np.clip((a - lo_) / (hi_ - lo_), 0, 1)
rgb = np.stack([stretch(col[2]), stretch(col[1]), stretch(col[0])])          # R 1000, G 750, B 480
I = rgb.mean(0); valid = np.isfinite(I) & (pan > 0)
P = (pan - pan[valid].mean()) / pan[valid].std() * np.nanstd(I[valid]) + np.nanmean(I[valid])
out = np.clip(rgb * (P / np.where(I > 0.02, I, 0.02))[None], 0, 1); out = np.nan_to_num(out) ** 0.9
out = np.moveaxis(out, 0, -1); out[~valid] = 0.15
M = np.load(f"{HERE}/hrp_mask.npz"); F = pd.read_csv(f"{HERE}/hrp_fields.csv")
def circle(clat, clon, Rkm, n=720):
    d = Rkm / R_M; b = np.radians(np.linspace(0, 360, n)); p1, l1 = np.radians(clat), np.radians(clon)
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(b))
    l2 = l1 + np.arctan2(np.sin(b) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2)); return np.degrees(p2), np.degrees(l2) % 360
fig, ax = plt.subplots(figsize=(12, 6.6), dpi=200)
ax.imshow(out, extent=list(ext), aspect=1 / np.cos(np.radians(50)), interpolation="lanczos")
ax.contour(M["lo"], M["la"], M["mask"].astype(float), levels=[0.5], colors="white", linewidths=1.1)
y, x = circle(50.05, 119.26, 350.3); ax.plot(x, y, color="#ff5a5a", lw=1.5, ls="--")
y, x = circle(50.05, 119.26, 1.2 * 350.3); ax.plot(x, y, color="white", lw=0.8, ls=":")
colors = {"on": "#ff4fd8", "near": "#ffd23f", "far": "#4cc9f0"}
for cat, g_ in F.groupby("cat"):
    ax.scatter(g_.lon, g_.lat, s=12 + 9 * np.sqrt(g_.n), c=colors[cat], edgecolors="black", linewidths=0.6, zorder=5)
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_xlabel("East longitude (°)"); ax.set_ylabel("Latitude (°N)")
h = [Line2D([], [], color="white", lw=1.1, label="higher-albedo plains (automated mask)"),
     Line2D([], [], color="#ff5a5a", lw=1.5, ls="--", label="b59 rim (Orgel et al. 2020)"),
     Line2D([], [], color="white", lw=0.8, ls=":", label="test area (1.2 b59 radii)"),
     Line2D([], [], marker="o", ls="", mfc=colors["on"], mec="k", ms=7, label=f"field on the plains ({(F.cat=='on').sum()})"),
     Line2D([], [], marker="o", ls="", mfc=colors["near"], mec="k", ms=7, label=f"field within 25 km of their edge ({(F.cat=='near').sum()})"),
     Line2D([], [], marker="o", ls="", mfc=colors["far"], mec="k", ms=7, label=f"field farther away ({(F.cat=='far').sum()})")]
leg = ax.legend(handles=h, loc="lower right", fontsize=7.5, framealpha=0.85, facecolor="#1a1a19", edgecolor="#6b6b67", labelcolor="white")
ax.set_title("Hollow fields and the higher-albedo smooth plains inside b59 (symbol area scales with the number of hollows)", fontsize=9)
fig.savefig(f"{HERE}/fig_hrp_mask_v2.png", bbox_inches="tight"); print("ok")
