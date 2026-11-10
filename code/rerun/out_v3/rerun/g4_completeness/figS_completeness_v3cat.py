# g4 copy of manuscripts/PlanetaryResearch_hollows/figscripts/make_figS_completeness_v8.py (6 Nov 2026). Changed: inputs (completeness grid completeness_10ppd_v3.npz built from the v3 refit detection model; theta from g3_catalog/sfd/out_v3_refit/theta_out_v2.json) and outputs (figS_completeness_v3cat.png/.pdf in this folder). Run: python3 figS_completeness_v3cat.py ~/mnt/hollows-emissivity
#!/usr/bin/env python3
"""Supplementary figure (v8, 30 Oct 2026): detection probability and catalog completeness.
(a) Completeness C for hollows 1 km across on the 10 pixel-per-degree grid (code/rerun/out_v2/completeness_10ppd.npz, key C at
    Dref = 1.0 km): global equirectangular map and north polar Lambert azimuthal equal-area map from 30 N (0 E at the bottom,
    as in Figure 1a). C = 1 - prod_k (1 - p_k) over searched MDIS NAC images covering the cell, built by
    code/rerun/build_hollow_products.py with each image's pixel scale replaced by its band's representative scale and incidence
    fixed at 60 deg.
(b) Per-image detection probability p(x, i) of the v2 detection model against hollow diameter in pixels x and solar incidence
    angle i, with the parameters refitted (and checked against the recorded 3-decimal values and against the stored grid) by
    analysis_v8/detection_function_v2_check.py, read from analysis_v8/detection_function_v2_check_record.json.
House style: DejaVu Sans, ink #1a1a19, muted #6b6b67; 15 cm wide. Writes figures/figS_completeness_v8.png (300 dpi) and .pdf.
Usage: python3 make_figS_completeness_v8.py [~/mnt/hollows-emissivity]"""
import os, sys, json
import numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pyproj
from scipy.special import erf

HERE = os.path.dirname(os.path.abspath(__file__)); PR = os.path.abspath(os.path.join(HERE, ".."))
ROOT = os.path.expanduser(sys.argv[1]) if len(sys.argv) > 1 else os.path.abspath(os.path.join(PR, "..", ".."))
OUT = HERE
names = [os.path.join(OUT, f"figS_completeness_v3cat.{e}") for e in ("png", "pdf")]
assert os.environ.get("REPLACE_V8") == "1" or not any(os.path.exists(n) for n in names), "outputs exist; set REPLACE_V8=1 to regenerate"

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "mathtext.fontset": "dejavusans",
    "axes.linewidth": .6, "axes.edgecolor": "#5a5a57", "xtick.color": "#5a5a57", "ytick.color": "#5a5a57",
    "xtick.major.width": .6, "ytick.major.width": .6, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.labelcolor": "#1a1a19", "text.color": "#1a1a19", "axes.labelsize": 7.5,
    "figure.facecolor": "white", "axes.facecolor": "white", "pdf.fonttype": 42})
INK, MUTED = "#1a1a19", "#6b6b67"

# ------------------------------------------------------------------ data
Z = np.load(f"{ROOT}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz")
C = Z["C"][list(Z["Dref"]).index(1.0)].astype(float)                     # (1800, 3600), row 0 = 89.95 N, col 0 = 0.05 E
TH = json.load(open(f"{ROOT}/code/rerun/out_v3/rerun/g3_catalog/sfd/out_v3_refit/theta_out_v2.json"))
def Phi(z): return 0.5 * (1 + erf(z / np.sqrt(2)))
def pdet(x, inc):
    p = TH["p_max"] * Phi((np.log(x) - np.log(TH["x50_px"])) / TH["sigma_lo"]) * (1 - Phi((np.log(x) - np.log(TH["x_hi_px"])) / TH["sigma_hi"]))
    return np.clip(p * np.exp(TH["b_inc"] * (inc - 60) / 25), 1e-6, 0.999)

# ------------------------------------------------------------------ layout (inches)
W, H = 15 / 2.54, 4.3
fig = plt.figure(figsize=(W, H))
def ax_in(x0, y0, w, h): return fig.add_axes([x0 / W, y0 / H, w / W, h / H])
ae = ax_in(0.48, 2.55, 2.90, 1.45)                     # equirectangular
ap = ax_in(3.44, 2.37, 1.66, 1.81)                     # north polar
ca = ax_in(5.17, 2.55, 0.08, 1.45)                     # colorbar (a)
ab = ax_in(0.48, 0.40, 4.55, 1.45)                     # heat map
cb = ax_in(5.17, 0.40, 0.08, 1.45)                     # colorbar (b)
fig.text(0.02 / W, 4.06 / H, "a", fontsize=10, fontweight="bold", va="bottom", ha="left")
fig.text(0.02 / W, 1.91 / H, "b", fontsize=10, fontweight="bold", va="bottom", ha="left")
cmapA = mpl.colormaps["cividis"]; normA = mpl.colors.Normalize(0, 1)

# (a) global
ae.imshow(C, extent=[0, 360, -90, 90], origin="upper", cmap=cmapA, norm=normA, interpolation="antialiased", aspect="auto")
ae.set_xlim(0, 360); ae.set_ylim(-90, 90)
ae.set_xticks([0, 90, 180, 270, 360]); ae.set_yticks([-90, -45, 0, 45, 90])
ae.set_xlabel("longitude (°E)", labelpad=2); ae.set_ylabel("latitude (°N)", labelpad=1)
ae.tick_params(length=2.5, pad=1.5)

# (a) north polar, Lambert azimuthal equal area from 30 N
Rm = 2439400.0
ps = pyproj.CRS.from_proj4(f"+proj=laea +lat_0=90 +lon_0=0 +R={Rm} +no_defs"); geo = pyproj.CRS.from_proj4(f"+proj=longlat +R={Rm} +no_defs")
half = 2 * Rm * np.sin(np.radians(60) / 2) * 1.0005; npx = 1400; pix = 2 * half / npx
xx = (np.arange(npx) + .5) * pix - half; Xg, Yg = np.meshgrid(xx, -xx)
lo_, la_ = pyproj.Transformer.from_crs(ps, geo, always_xy=True).transform(Xg, Yg); lo_ = np.mod(lo_, 360)
rr = ((90 - la_) * 10).astype(int).clip(0, 1799); cc = (lo_ * 10).astype(int).clip(0, 3599)
img = np.where(la_ >= 30, C[rr, cc], np.nan)
hk = half / 1e3
ap.imshow(img, extent=[-hk, hk, -hk, hk], origin="upper", cmap=cmapA, norm=normA, interpolation="antialiased", zorder=0)
th = np.linspace(0, 2 * np.pi, 400)
for lc in (30, 45, 60, 75):
    rad = 2 * Rm * np.sin(np.radians(90 - lc) / 2) / 1e3
    ap.plot(rad * np.cos(th), rad * np.sin(th), color=MUTED if lc == 30 else "white", lw=.6 if lc == 30 else .35,
            ls="-" if lc == 30 else (0, (2, 2)), alpha=1 if lc == 30 else .55, zorder=2)
    if lc > 30:
        ap.text(rad * np.cos(np.radians(-58)) + 40, rad * np.sin(np.radians(-58)) - 40, f"{lc}°", fontsize=5.6, color="white", ha="left", va="top", zorder=3)
rad30 = 2 * Rm * np.sin(np.radians(60) / 2) / 1e3
for lg_ in (0, 90, 180, 270):
    ang = np.radians(lg_ - 90)
    ap.plot([0, rad30 * np.cos(ang)], [0, rad30 * np.sin(ang)], color="white", lw=.35, ls=(0, (2, 2)), alpha=.55, zorder=2)
    f = 1.09 if lg_ in (0, 180) else 1.155
    ap.text(rad30 * f * np.cos(ang), rad30 * f * np.sin(ang), f"{lg_}°E", fontsize=5.6, rotation=90 if lg_ in (90, 270) else 0, color=MUTED, ha="center", va="center")
ap.text(rad30 * np.cos(np.radians(-50)) + 40, rad30 * np.sin(np.radians(-50)) - 40, "30°N", fontsize=5.6, color=MUTED, ha="left", va="top")
lim = rad30 * 1.27
ap.set_xlim(-lim, lim); ap.set_ylim(-rad30 * 1.385, rad30 * 1.385); ap.set_aspect("equal"); ap.axis("off")
cbA = fig.colorbar(mpl.cm.ScalarMappable(norm=normA, cmap=cmapA), cax=ca)
cbA.set_label("completeness, $C$\n(1 km hollows)", labelpad=3, linespacing=1.1); cbA.set_ticks([0, .2, .4, .6, .8, 1]); cbA.outline.set_linewidth(.6)
ca.tick_params(length=2.5, pad=1.5)

# (b) detection probability
x = np.logspace(np.log10(2), np.log10(1000), 500); inc = np.linspace(0, 85, 341)
Xm, Im = np.meshgrid(x, inc); Pm = pdet(Xm, Im)
cmapB = mpl.colormaps["viridis"]; normB = mpl.colors.Normalize(0, 0.3)
ab.pcolormesh(x, inc, Pm, cmap=cmapB, norm=normB, shading="auto", rasterized=True)
LEV = [0.05, 0.10, 0.15, 0.20, 0.25]
cs = ab.contour(x, inc, Pm, levels=LEV, colors="white", linewidths=.45, alpha=.75)
ab.set_xscale("log"); ab.set_xlim(2, 1000); ab.set_ylim(0, 85)
ab.set_xticks([2, 5, 10, 20, 50, 100, 200, 500, 1000]); ab.set_xticklabels(["2", "5", "10", "20", "50", "100", "200", "500", "1000"]); ab.minorticks_off()
ab.set_yticks([0, 15, 30, 45, 60, 75])
ab.set_xlabel("hollow diameter (pixels)", labelpad=2); ab.set_ylabel("solar incidence angle (°)", labelpad=2)
ab.tick_params(length=2.5, pad=1.5)
cbB = fig.colorbar(mpl.cm.ScalarMappable(norm=normB, cmap=cmapB), cax=cb)
cbB.set_label("detection probability\nper image, $p$", labelpad=3, linespacing=1.1); cbB.set_ticks([0, .05, .1, .15, .2, .25, .3]); cbB.outline.set_linewidth(.6); cbB.add_lines(cs)
cb.tick_params(length=2.5, pad=1.5)

fig.savefig(names[0], dpi=300); fig.savefig(names[1], dpi=300)
print("wrote", *names, "p range", float(Pm.min()), float(Pm.max()))
