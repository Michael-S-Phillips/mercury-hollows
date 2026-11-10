#!/usr/bin/env python3
"""(v2, 5 Nov 2026: all 6,667 cataloged hollows from hollows_canonical.csv; v1 used host_assignments.csv, which holds only the 5,409 hollows inside host craters.)
Map of the region of Fassett et al. (2009, EPSL 285:297) Figure 9, northwest of Caloris (5 Nov 2026): MDIS MDR color
(1000/750/480 nm as R/G/B, as in their figure) with the cataloged hollows, Orgel et al. (2020) basins b59 and Caloris,
Michael's candidate basins (GIS/vector/possiblebasins.gpkg, north polar stereographic), and the approximate targets of
Fassett's white arrows read from their figure. Writes fig_fassett_region.png and fig_fassett_region_record.json."""
import os, sys, json, sqlite3, struct, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
sys.path.insert(0, f"{HERE}/../color_mdr"); import mdr_georef as g
LAT0, LAT1, LON0, LON1, STEP = 36.0, 62.0, 100.0, 170.0, 0.05
la = np.arange(LAT1, LAT0, -STEP); lo = np.arange(LON0, LON1, STEP); LO, LA = np.meshgrid(lo, la)
G = g.Globe()
bands = {b: G.sample(LA.ravel(), LO.ravel(), band=b).reshape(LA.shape) for b in (1, 4, 7)}   # 480, 750, 1000 nm
def stretch(a, p=(2, 98)):
    lo_, hi_ = np.nanpercentile(a, p); return np.clip((a - lo_) / (hi_ - lo_), 0, 1)
rgb = np.dstack([stretch(bands[7]), stretch(bands[4]), stretch(bands[1])]); rgb[np.isnan(rgb)] = 0
# candidate basins from the GeoPackage (north polar stereographic, sphere R = 2,440,530 m)
def wkb_polygon(w):
    bo = '<' if w[0] == 1 else '>'; nr = struct.unpack(bo + 'I', w[5:9])[0]; off = 9; rings = []
    for _ in range(nr):
        n = struct.unpack(bo + 'I', w[off:off + 4])[0]; off += 4
        pts = np.frombuffer(w[off:off + 16 * n], dtype=bo + 'f8').reshape(n, 2); off += 16 * n; rings.append(pts)
    return rings
c = sqlite3.connect(f"{ROOT}/GIS/vector/possiblebasins.gpkg"); cands = []
for fid, gb in c.execute("select fid, geom from possiblebasins"):
    env = (gb[3] >> 1) & 7; w = gb[8 + {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]:]
    xy = wkb_polygon(w)[0]; rho = np.hypot(xy[:, 0], xy[:, 1]); Rs = 2440530.0
    lat = 90 - np.degrees(2 * np.arctan(rho / (2 * Rs))); lon = np.degrees(np.arctan2(xy[:, 0], -xy[:, 1])) % 360
    cands.append(dict(fid=fid, lat=lat, lon=lon, clat=float(np.mean(lat)), clon=float(np.degrees(np.arctan2(np.mean(np.sin(np.radians(lon))), np.mean(np.cos(np.radians(lon))))) % 360)))
def circle(clat, clon, Rkm, n=360):
    d = Rkm / R_M; b = np.radians(np.linspace(0, 360, n)); p1, l1 = np.radians(clat), np.radians(clon)
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(b))
    l2 = l1 + np.arctan2(np.sin(b) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2))
    return np.degrees(p2), np.degrees(l2) % 360
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)]
h["lon360"] = h.lon % 360; hh = h[h.lat.between(LAT0, LAT1) & h.lon360.between(LON0, LON1)]
arrows = [(45.6, 119.5), (48.3, 122.5), (49.8, 128.5), (48.4, 130.5), (52.3, 129.5)]   # approximate, read from Fassett et al. Fig. 9
fig, ax = plt.subplots(figsize=(13, 5.6), dpi=170)
ax.imshow(rgb, extent=[LON0, LON1, LAT0, LAT1], aspect=1 / np.cos(np.radians(49)), origin="upper")
ax.scatter(hh.lon360, hh.lat, s=4, c="#ffd23f", edgecolors="none", alpha=0.9, label=f"cataloged hollows ({len(hh):,})", zorder=4)
for name, (cl, cn, D, col) in {"b59 (Orgel et al. 2020)": (50.05, 119.26, 700.6, "#ff5a5a"), "Caloris rim (Orgel et al. 2020)": (29.91, 161.43, 1555.9, "#ffffff")}.items():
    y, x = circle(cl, cn, D / 2); ax.plot(x, y, color=col, lw=1.4, ls="--", label=name, zorder=5)
for k, cnd in enumerate(cands):
    if LAT0 - 10 < cnd["clat"] < LAT1 + 10 and LON0 - 10 < cnd["clon"] < LON1 + 10:
        ax.plot(cnd["lon"], cnd["lat"], color="#5ad1ff", lw=1.2, label="candidate basins (possiblebasins.gpkg)" if k == 0 or "cand" not in [l.get_label()[:4] for l in ax.lines] else None, zorder=5)
for (a, b_) in arrows: ax.annotate("", xy=(b_, a), xytext=(b_ - 2.5, a + 1.6), arrowprops=dict(arrowstyle="->", color="white", lw=1.6), zorder=6)
ax.set_xlim(LON0, LON1); ax.set_ylim(LAT0, LAT1); ax.set_xlabel("East longitude (°)"); ax.set_ylabel("Latitude (°N)")
hl, lb = ax.get_legend_handles_labels(); seen = {}; [seen.setdefault(l, hdl) for hdl, l in zip(hl, lb) if l]
ax.legend(seen.values(), seen.keys(), loc="lower left", fontsize=7, framealpha=0.85)
ax.set_title("MDIS MDR color (R 1000, G 750, B 480 nm), the region of Fassett et al. (2009) Fig. 9; white arrows approximate theirs", fontsize=9)
fig.savefig(f"{HERE}/fig_fassett_region.png", bbox_inches="tight")
rec = {"candidates": [dict(fid=x["fid"], center_lat=round(x["clat"], 2), center_lon=round(x["clon"], 2)) for x in cands],
       "arrows_approx": arrows, "n_hollows_in_frame": int(len(hh))}
json.dump(rec, open(f"{HERE}/fig_fassett_region_record.json", "w"), indent=1); print(json.dumps(rec, indent=1))
