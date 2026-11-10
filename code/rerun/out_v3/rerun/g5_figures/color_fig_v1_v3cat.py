# g5 copy of manuscripts/PlanetaryResearch_hollows/figscripts/color_fig_v1.py (6 Nov 2026). Changed: host flag h09 from g1_graph_joint/summary_model/v3/crater_features_xrs.csv and the depth fit from g1_graph_joint/summary_model/v3/out/color_record.json (v3 catalog); run with output g5_figures/figS7_color_v1_v3cat.png. The crater and hollow colors (analysis_workshop/color_mdr) still use the v2 hollow mask (not rerun by any group).
#!/usr/bin/env python3
"""Figure S10 (v6 / SM v11, 3 Nov 2026): MDIS eight-filter color of crater interiors and hollow settings.
(a) Loadings of the first three principal components of brightness-normalized spectra (each spectrum divided by its
    geometric mean) of 31,373 crater interiors (hollows masked), from color_mdr/pca_record.json.
(b) PC1 (spectral slope) of each crater interior; host craters (a hollow within 0.9 R) ringed.
(c) Hollow depth against PC1 of the surroundings within 10 km (hollows masked), in quintiles: depth corrected for
    ln pixel scale, ln tan(incidence), and |latitude| (OLS without color), median per quintile relative to the overall
    median, 95% intervals from 1,000 bootstrap resamples of images; line = fitted 0.903 per SD (color_record.json).
Inputs: analysis_workshop/color_mdr/{pca_record.json, crater_color.csv, hollow_color.csv},
summary_model/crater_features_xrs.csv (h09), GRL_insolation_hollows/data/test1_rows_inc.csv.
Usage: python3 color_fig_v1.py <manuscripts dir> <output png>"""
import sys, os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf  # g5: os imported here (used above its original import)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

M = sys.argv[1] if len(sys.argv) > 1 else "/mnt/user-data/uploads/hollows-emissivity/manuscripts"
OUTP = sys.argv[2] if len(sys.argv) > 2 else "figS_color_v1.png"
W = "%s/PlanetaryResearch_hollows/analysis_workshop" % M
CDIR = "%s/color_mdr" % W
# COLOR_SRC=screened (added 9 Nov 2026): PCA, crater colors and hollow colors from the screened-mask run
# (g6_remaining/color/out_v3, 31,412 craters) and the depth fit from its depth_color_fits.json, in place of the version-2 mask.
SCREENED = os.environ.get("COLOR_SRC") == "screened"
if SCREENED: CDIR = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g6_remaining/color/out_v3")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#555555", "axes.labelcolor": "#222222",
                     "xtick.color": "#555555", "ytick.color": "#555555", "axes.spines.top": False, "axes.spines.right": False})
INK, MUTED = "#222222", "#6b6b66"
SER = ["#2a78d6", "#eb6834", "#1baf7a"]          # categorical slots 1-3 (validated, light)
DIV = LinearSegmentedColormap.from_list("bluered", ["#0d366b", "#2a78d6", "#86b6ef", "#f0efec", "#f2a2a1", "#e34948", "#8f2423"])

pr = json.load(open("%s/pca_record.json" % CDIR))
wl = [430, 480, 560, 630, 750, 830, 900, 1000]
c = pd.read_csv("%s/crater_color.csv" % CDIR, dtype={"id": str}, usecols=["id", "lat", "lon", "D", "ok", "pc1"])
G1S = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g1_graph_joint/summary_model/v3")  # g5: g1 joint-model rerun, v3 catalog
f = pd.read_csv("%s/crater_features_xrs.csv" % G1S, dtype={"id": str}, usecols=["id", "h09"])  # g5: v3 h09 (was analysis_workshop/summary_model/crater_features_xrs.csv)
c = c[c.ok.fillna(False).astype(bool)].merge(f, on="id", how="left")
t = pd.read_csv("%s/GRL_insolation_hollows/data/test1_rows_inc.csv" % M); t["row"] = np.arange(len(t))
h = pd.read_csv("%s/hollow_color.csv" % CDIR, usecols=["row", "set_pc1", "set_pc2"])
t = t.merge(h, on="row", how="left").dropna(subset=["set_pc1"]).reset_index(drop=True)
t["ltan"] = np.log(np.tan(np.radians(t.inc)))

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.42, wspace=0.32)
axa, axc, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

# (a) loadings
vr = pr["variance_ratio"]
for k, (pc, col, ls) in enumerate(zip(["PC1", "PC2", "PC3"], SER, ["-", "-", "--"])):
    y = [pr["loadings"][pc][str(w)] for w in wl]
    axa.plot(wl, y, ls, color=col, lw=2, marker="o", ms=4, mec="white", mew=0.8, zorder=3)
    lab = {"PC1": "PC1, %.0f%% (slope)", "PC2": "PC2, %.0f%% (curvature)", "PC3": "PC3, %.0f%%"}[pc] % (100 * vr[k])
    axa.lines[-1].set_label(lab)
axa.legend(loc="lower right", frameon=False, fontsize=7.5, handlelength=2.2, borderaxespad=0.1)
axa.axhline(0, color="#cccccc", lw=0.8, zorder=1)
axa.axvspan(600, 650, color="#eeeeee", zorder=0); axa.text(625, 0.92, "600 nm\nband", ha="center", va="top", fontsize=7, color=MUTED)
axa.set_xlim(400, 1010); axa.set_ylim(-1.02, 0.93); axa.set_xticks([400, 600, 800, 1000])
axa.set_xlabel("wavelength (nm)"); axa.set_ylabel("loading")
axa.set_title("a", loc="left", fontweight="bold")

# (c) depth vs setting PC1
g = pd.factorize(t.filename)[0]
base = smf.ols("lnh ~ lnpx + ltan + abslat", t).fit()
t["res"] = base.resid
sd = t.set_pc1.std(); t["z"] = (t.set_pc1 - t.set_pc1.mean()) / sd
t["q"] = pd.qcut(t.z, 5, labels=False)
cent = t.groupby("q").z.median().values
med0 = np.median(t.res)
pt = np.exp(t.groupby("q").res.median().values - med0)
rng = np.random.default_rng(7); imgs = np.unique(g); B = []
for _ in range(1000):
    pick = rng.choice(imgs, len(imgs)); cnt = np.bincount(pick, minlength=g.max() + 1); w = cnt[g]
    rows = np.repeat(np.arange(len(t)), w); tt = t.iloc[rows]
    B.append(np.exp(tt.groupby("q").res.median().values - np.median(tt.res)))
B = np.array(B); lo, hi = np.percentile(B, 2.5, axis=0), np.percentile(B, 97.5, axis=0)
import os
_r = "%s/out/color_record.json" % G1S  # g5: v3 record (same depth fit, 0.903 per SD)
rec = json.load(open(_r if os.path.exists(_r) else "/home/claude/summary_model/out/color_record.json"))
b = rec["depth"]["fits"]["abslat"]["set_pc1_z"][0]
if SCREENED: b = json.load(open("%s/depth_color_fits.json" % CDIR))["fits"]["abslat"]["set_pc1_z"][0]
zz = np.linspace(cent[0] - 0.3, cent[-1] + 0.3, 50)
axc.plot(zz, b ** zz, color=MUTED, lw=1.2, ls="--", zorder=2)
axc.text(zz[-1], b ** zz[-1] * 0.985, "fit: %.2f per SD" % b, ha="right", va="top", fontsize=8, color=MUTED)
axc.vlines(cent, lo, hi, color=SER[0], lw=2, zorder=3)
axc.plot(cent, pt, "o", color=SER[0], ms=8, mec="white", mew=1.5, zorder=4)
axc.axhline(1, color="#cccccc", lw=0.8, zorder=1)
axc.set_xlabel("PC1 of surroundings within 10 km (SD)\nbluer  ←  →  redder")
axc.set_ylabel("corrected depth / median")
axc.set_title("b", loc="left", fontweight="bold")
axc.set_ylim(0.7, 1.3)

# (b) map
cc = c.sort_values("D")
norm = TwoSlopeNorm(vcenter=0, vmin=-0.15, vmax=0.15)
s = 1.5 + 6 * (np.log10(cc.D) - np.log10(5))
sc = axb.scatter(cc.lon % 360, cc.lat, c=cc.pc1, s=s, cmap=DIV, norm=norm, lw=0, rasterized=True, zorder=2)
hh = cc[cc.h09 == 1]
axb.scatter(hh.lon % 360, hh.lat, s=14, facecolors="none", edgecolors=INK, lw=0.6, zorder=3)
axb.set_xlim(0, 360); axb.set_ylim(-90, 90); axb.set_xticks(range(0, 361, 60)); axb.set_yticks(range(-90, 91, 30))
axb.set_xlabel("east longitude (°)"); axb.set_ylabel("latitude (°)")
axb.set_aspect("equal")
cb = fig.colorbar(sc, ax=axb, fraction=0.025, pad=0.015, extend="both")
cb.set_label("crater-interior PC1\n(bluer ← → redder)", fontsize=8); cb.outline.set_edgecolor("#555555")
axb.scatter([], [], s=14, facecolors="none", edgecolors=INK, lw=0.6, label="host crater")
axb.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=8, handletextpad=0.3, borderaxespad=0.2)
axb.set_title("c", loc="left", fontweight="bold")

fig.savefig(OUTP, dpi=300, bbox_inches="tight", facecolor="white")
json.dump({"quintile_center_sd": cent.tolist(), "ratio": pt.tolist(), "lo": lo.tolist(), "hi": hi.tolist(), "fit_per_sd": b,
           "n_craters": int(len(c)), "n_hosts": int((c.h09 == 1).sum()), "n_hollows": int(len(t))},
          open(OUTP.replace(".png", "_record.json"), "w"), indent=1)
print("wrote", OUTP)
