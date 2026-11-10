#!/usr/bin/env python3
"""g20 paper figure (9 Nov 2026): main-text Figure 6 of v6.9. Same data as fig_clusters.py (grid.npz, clusters.json, craters_P.csv;
500 km windows), drawn at the paper's 15 cm width in the style of make_figs_pr_v9 (DejaVu Sans, INK/MUTED). Hillshade of the USGS 665 m
DEM (Becker et al. 2016) at overview level 8. Caloris rim as in Figure 2a (Orgel et al. 2020: 1,556 km, 29.9 N, 161.4 E).
Writes fig6_regional_excess_v1.png (400 dpi) and .pdf."""
import os, json, numpy as np, pandas as pd, rasterio
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LightSource
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5})
INK, MUTED = "#1a1a19", "#6b6b67"; R = 2439.7
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
G = np.load(f"{HERE}/grid.npz"); S = json.load(open(f"{HERE}/clusters.json")); C = pd.read_csv(f"{HERE}/craters_P.csv")
with rasterio.open(f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif") as s:
    z = s.read(1, out_shape=(1440, 2880)).astype(float)
z = np.roll(z, 1440, axis=1)
hs = LightSource(azdeg=315, altdeg=35).hillshade(z, vert_exag=0.02, dx=1, dy=1)
glat, glon = G["glat"], G["glon"]; O, E, Z, ok = G["O"], G["E"], G["Z"], G["ok"]
ratio = np.where(ok, np.log2(np.clip((O + 0.5) / (E + 0.5), 0.25, 4)), np.nan)
cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
fig = plt.figure(figsize=(5.9, 3.25)); ax = fig.add_axes([0.075, 0.25, 0.9, 0.69])
ax.imshow(hs, extent=[0, 360, -90, 90], cmap="gray", vmin=-0.2, vmax=1.2, interpolation="bilinear", rasterized=True)
pm = ax.pcolormesh(np.arange(0, 361, 2.0), np.arange(-90, 91, 2.0), ratio, cmap=cmap, vmin=-2, vmax=2, alpha=0.62, shading="flat", edgecolors="none", rasterized=True)
ax.contour(glon, glat, np.where(ok, Z, -10), levels=[S["z_threshold_excess"]], colors=INK, linewidths=0.9)
h = C[C.h09 == 1]; ax.scatter(h.lon % 360, h.lat, s=1.6, c=INK, lw=0, zorder=4)
la0, lo0, ang = np.radians(29.91), np.radians(161.43), (1555.9 / 2) / R; t = np.linspace(0, 2 * np.pi, 721)
lat = np.arcsin(np.sin(la0) * np.cos(ang) + np.cos(la0) * np.sin(ang) * np.cos(t))
lon = lo0 + np.arctan2(np.sin(t) * np.sin(ang) * np.cos(la0), np.cos(ang) - np.sin(la0) * np.sin(lat))
ax.plot(np.degrees(lon), np.degrees(lat), ls=(0, (4, 1.5)), c=INK, lw=0.8, zorder=5)
ax.text(161.4, 29.9, "Caloris", ha="center", va="center", fontsize=6.2, color=INK, zorder=5)
for s_, x, y_ in (("1", 133, 81), ("2", 316, 25), ("3", 216, 17)):
    ax.text(x, y_, s_, fontsize=7, fontweight="bold", ha="center", va="center", color=INK,
            bbox=dict(boxstyle="circle,pad=0.2", fc="#fcfcfb", ec=INK, lw=0.6), zorder=6)
ax.set_xlim(0, 360); ax.set_ylim(-60, 90); ax.set_xticks(range(0, 361, 30)); ax.set_yticks(range(-60, 91, 30))
ax.set_xticklabels([f"{v}°E" for v in range(0, 361, 30)], fontsize=6.2); ax.set_yticklabels([f"{v}°" for v in range(-60, 91, 30)], fontsize=6.2)
for v, nm in ((0, "hot"), (90, "warm"), (180, "hot"), (270, "warm"), (360, "hot")):
    ax.text(v, 91.5, nm, ha="center", va="bottom", fontsize=6.2, color=MUTED)
cax = fig.add_axes([0.075, 0.085, 0.3, 0.03]); cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
cb.set_ticks([-2, -1, 0, 1, 2]); cb.set_ticklabels(["¼", "½", "1", "2", "4"]); cb.ax.tick_params(labelsize=6.2, length=2); cb.outline.set_linewidth(0.5)
fig.text(0.39, 0.10, "Observed / predicted host craters within 500 km", fontsize=6.6, va="center", color=INK)
fig.savefig(f"{HERE}/fig6_regional_excess_v1.png", dpi=400); fig.savefig(f"{HERE}/fig6_regional_excess_v1.pdf", dpi=400); print("ok")
