#!/usr/bin/env python3
"""g20 figure (9 Nov 2026): observed / expected host craters within 500 km of each 2-degree node, expected from the joint model
(clusters.py), on a hillshade of the USGS 665 m DEM (overview level 16). Outlines: family-wise significant excess (z >= threshold
in clusters.json). Writes fig_regional_clusters.png."""
import os, json, numpy as np, pandas as pd, rasterio
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LightSource
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
G = np.load(f"{HERE}/grid.npz"); S = json.load(open(f"{HERE}/clusters.json")); C = pd.read_csv(f"{HERE}/craters_P.csv")
with rasterio.open(f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif") as s:
    z = s.read(1, out_shape=(720, 1440)).astype(float)
z = np.roll(z, 720, axis=1)                                            # -180..180 -> 0..360 E
hs = LightSource(azdeg=315, altdeg=35).hillshade(z, vert_exag=0.02, dx=1, dy=1)
glat, glon = G["glat"], G["glon"]; O, E, Z, ok = G["O"], G["E"], G["Z"], G["ok"]
ratio = np.where(ok, np.log2(np.clip((O + 0.5) / (E + 0.5), 0.25, 4)), np.nan)
cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
fig = plt.figure(figsize=(12, 7.2), dpi=200); ax = fig.add_axes([0.06, 0.25, 0.9, 0.64])
ax.imshow(hs, extent=[0, 360, -90, 90], cmap="gray", vmin=-0.2, vmax=1.2, interpolation="bilinear")
lon_e = np.arange(0, 361, 2.0); lat_e = np.arange(-90, 91, 2.0)
pm = ax.pcolormesh(lon_e, lat_e, ratio, cmap=cmap, vmin=-2, vmax=2, alpha=0.62, shading="flat", edgecolors="none")
thr = S["z_threshold_excess"]; Zc = np.where(ok, Z, -10)
ax.contour(glon, glat, Zc, levels=[thr], colors="#1a1a19", linewidths=1.3)
h = C[C.h09 == 1]; ax.scatter(h.lon % 360, h.lat, s=4, c="#1a1a19", lw=0, zorder=4, label="host crater")
cal = (29.23, 161.40, 695.6); t = np.linspace(0, 2 * np.pi, 361); ang = cal[2] / 2439.7
la0, lo0 = np.radians(cal[0]), np.radians(cal[1])
lat = np.arcsin(np.sin(la0) * np.cos(ang) + np.cos(la0) * np.sin(ang) * np.cos(t))
lon = lo0 + np.arctan2(np.sin(t) * np.sin(ang) * np.cos(la0), np.cos(ang) - np.sin(la0) * np.sin(lat))
ax.plot(np.degrees(lon), np.degrees(lat), ls=(0, (4, 3)), c="#1a1a19", lw=0.9); ax.text(161.4, 29.2, "Caloris", ha="center", va="center", fontsize=8, color="#1a1a19")
labels = [("1", 132, 79), ("2", 318, 24), ("3", 214, 17)]
for s_, x, y_ in labels:
    ax.text(x, y_, s_, fontsize=11, fontweight="bold", ha="center", va="center", color="#1a1a19",
            bbox=dict(boxstyle="circle,pad=0.25", fc="#fcfcfb", ec="#1a1a19", lw=0.8), zorder=6)
ax.set_xlim(0, 360); ax.set_ylim(-60, 90); ax.set_xticks(range(0, 361, 30)); ax.set_yticks(range(-60, 91, 30))
ax.set_xticklabels([f"{v}°E" for v in range(0, 361, 30)], fontsize=8); ax.set_yticklabels([f"{v}°" for v in range(-60, 91, 30)], fontsize=8)
ax.grid(color="#fcfcfb", lw=0.4, alpha=0.5)
for v, nm in ((0, "hot"), (90, "warm"), (180, "hot"), (270, "warm"), (360, "hot")):
    ax.text(v, 91.5, nm, ha="center", va="bottom", fontsize=8, color="#555")
fig.text(0.06, 0.965, "Hosting is clustered beyond the joint model in three regions", fontsize=13, fontweight="bold", ha="left")
fig.text(0.06, 0.937, f"Observed / expected host craters within 500 km, expected from the joint model of Section 3.8 ({S['n_hosts']} hosts among {S['n_craters']:,} craters). "
         f"Black outlines: excess at family-wise p < 0.05 (z ≥ {thr:.1f}).", fontsize=8.5, ha="left", color="#333")
cax = fig.add_axes([0.06, 0.165, 0.3, 0.022]); cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
cb.set_ticks([-2, -1, 0, 1, 2]); cb.set_ticklabels(["¼", "½", "1", "2", "4"]); cb.ax.tick_params(labelsize=8); cb.outline.set_linewidth(0.5)
cb.set_label("Observed / expected hosts within 500 km (uncolored: < 2 expected)", fontsize=8)
ex = S["excess"]
def agg(sel):
    a = [c for c in ex if sel(c)]; return a
r1 = agg(lambda c: c["peak_lat"] > 40); r2 = agg(lambda c: 300 < c["peak_lon"] < 360 and c["peak_lat"] < 40); r3 = agg(lambda c: 180 < c["peak_lon"] < 240)
txt = ["1  Northern province; 45–75°N, 95–145°E: 23 hosts within 500 km of 65°N 125°E, 4.9 expected (z = 9.4).",
       "2  Abu Nuwas, Sinan, Ts'ai Wen-chi; 5–25°N, 325–340°E: 23 hosts near 13°N 329°E, 8.0 expected (z = 6.1).",
       "3  Around Tyagaraja; 0–10°N, 195–210°E: 19 hosts near 7°N 203°E, 6.2 expected (z = 6.0).",
       f"No region has significantly fewer hosts than expected (lowest z = −3.19; threshold −{abs(S['z_threshold_deficit']):.2f})."]
for i, s_ in enumerate(txt): fig.text(0.40, 0.185 - i * 0.032, s_, fontsize=8, ha="left", va="center", color="#1a1a19")
fig.savefig(f"{HERE}/fig_regional_clusters_v2.png", dpi=200); print("ok")
