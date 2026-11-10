# g13 revision (7 Nov 2026): G13_REV=r2 (default) uses the r2 traces (through traced_classes.py) and writes *_r2 figures.
#!/usr/bin/env python3
"""g13 (6 Nov 2026), from g10_threshold/figS2_bdr_v3.py: Figure S2 on the MDIS BDR 166 m mosaic with a light hillshade of the 665 m
topography, terrain classes from Michael's traced outlines (g13_traced/traced_classes.py; all four craters are traced). Same four
craters and drawing code (g10 compare_curv_v2v3.py load/draw). Writes figS2_terrain_examples_traced.png/.pdf here."""
import sys, os, numpy as np, matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt, rasterio
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
OUTD = os.path.dirname(os.path.abspath(__file__)); G10 = os.path.join(os.path.dirname(os.path.dirname(OUTD)), "g10_threshold"); sys.path.insert(0, G10)
sys.path.insert(0, os.path.dirname(OUTD)); from traced_classes import classes_for, TR
sys.argv = [sys.argv[0], "_none", "Tyagaraja"]
src = open(f"{G10}/compare_curv_v2v3.py").read(); src = src[:src.index("out, names = sys.argv[1]")]; exec(src)
picks = ["Tyagaraja", "Praxiteles", "Scarlatti", "Sander"]
mpl.rcParams.update({"font.size": 7.5})
fig, axs = plt.subplots(2, 2, figsize=(5.9, 6.3))
with rasterio.open(DEM) as ds, rasterio.open(BDR) as bs:
    for a, nm in zip(axs.ravel(), picks):
        c = cr[cr.name == nm].sort_values("diam").iloc[-1]; assert int(c.unique_id) in TR, nm
        z, LON, LAT, dx, dy, T, b, bext, R, cl = load(c, ds, bs); rho = np.hypot(dx, dy) / R; az = np.degrees(np.arctan2(dy, dx)) % 360; cent = str(c.cent_str).strip()
        kind = "peak-ring crater" if cent in ("pr", "rpc", "pb") else "central-peak crater"
        draw(a, c, z, LON, LAT, dx, dy, T, b, bext, R, cl, classes_for(c.unique_id, z, rho, az, R, cent, LON, LAT), f"{nm}, {c.diam:.0f} km, {kind}")
        a.set_xlabel("east (km)", fontsize=7); a.set_ylabel("north (km)", fontsize=7); a.title.set_ha("left"); a.title.set_position((0, 1.0))
h = [Patch(fc=COL[k], ec=COL[k] if k in ("peak", "wall") else "none", alpha=.5, label=l) for k, l in [("peak", "peak or ring"), ("floor", "floor"), ("wall", "wall"), ("masked", "superposed crater (masked)")]]
h += [Line2D([], [], marker="o", ms=3.5, mfc="#ffd84a", mec="#1a1a19", mew=.4, ls="none", label="hollow assigned to host"), Line2D([], [], color="#bbbbbb", lw=.8, ls=(0, (3, 2)), label="cataloged rim")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=6.6)
fig.subplots_adjust(left=.09, right=.99, top=.93, bottom=.14, wspace=.30, hspace=.42)
for k, a in zip("abcd", axs.ravel()): a.text(-0.2, 1.12, k, transform=a.transAxes, fontsize=10, fontweight="bold")
SUF = "" if os.environ.get("G13_REV", "r2") == "r1" else "_r2"; fig.savefig(f"{OUTD}/figS2_terrain_examples_traced{SUF}.png", dpi=400); fig.savefig(f"{OUTD}/figS2_terrain_examples_traced{SUF}.pdf"); print("wrote")
