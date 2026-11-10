#!/usr/bin/env python3
"""g4 (6 Nov 2026): compare the completeness grid built from the refit detection model (completeness_10ppd_v3.npz, v3b alike)
with the v2 grid (code/rerun/out_v2/completeness_10ppd.npz) over covered cells (>= 1 searched image and old C > 0).
Per reference diameter: 5th/50th/95th percentiles of new/old over cells (unweighted and area-weighted), and the same by
finest-pixel-scale class of the cell and by old-C class. Writes grid_comparison.json and grid_comparison.png."""
import os, json, numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); RT = os.path.expanduser("~/mnt/hollows-emissivity")
O = np.load(f"{RT}/code/rerun/out_v2/completeness_10ppd.npz"); Dref = [float(x) for x in O["Dref"]]
cv = np.load(f"{H}/_scratch/coverage_10ppd_rebuilt.npz"); nimg = cv["nimg"]; bestpx = cv["bestpx"]
lat = 90 - (np.arange(1800) + .5) / 10; W = np.repeat(np.cos(np.radians(lat))[:, None], 3600, 1)
def wq(x, w, q):
    o = np.argsort(x); cw = np.cumsum(w[o]); cw /= cw[-1]; return [float(np.interp(qq / 100, cw, x[o])) for qq in q]
PXC = [(0, 10), (10, 20), (20, 35), (35, 60), (60, 100), (100, 150.01)]
CC = [(0, .05), (.05, .2), (.2, .5), (.5, .8), (.8, .95), (.95, 1.01)]
R = {"definition": "ratio = C_new / C_old over cells with >= 1 searched image and C_old > 0; percentiles 5/50/95; rel = (p - median) / median",
     "criterion_1km": "uniform if p5 and p95 lie within +-3% of the median at Dref = 1 km"}
for tag, fn in (("v3", "completeness_10ppd_v3.npz"), ("v3b", "completeness_10ppd_v3b.npz")):
    N = np.load(f"{H}/{fn}"); R[tag] = {}
    for k, D in enumerate(Dref):
        o = O["C"][k]; n = N["C"][k]; m = (nimg > 0) & (o > 0); r = (n[m] / o[m]).astype(float); w = W[m]
        q = [float(v) for v in np.percentile(r, [5, 50, 95])]; qw = wq(r, w, [5, 50, 95])
        e = dict(cells=int(m.sum()), p5_p50_p95=q, rel_p5_p95=[(q[0] - q[1]) / q[1], (q[2] - q[1]) / q[1]],
                 area_weighted_p5_p50_p95=qw, rel_area_weighted=[(qw[0] - qw[1]) / qw[1], (qw[2] - qw[1]) / qw[1]],
                 mean_C_old=float(np.average(o[m], weights=w)), mean_C_new=float(np.average(n[m], weights=w)),
                 max_abs_dC=float(np.abs(n[m] - o[m]).max()), median_dC=float(np.median(n[m] - o[m])))
        if D in (0.5, 1.0, 2.0, 5.0):
            bp = bestpx[m]; e["by_finest_px_class"] = {}
            for a, b in PXC:
                s = (bp >= a) & (bp < b)
                if s.sum(): e["by_finest_px_class"][f"{a:g}-{b:g} m"] = dict(cells=int(s.sum()), p5_p50_p95=[float(v) for v in np.percentile(r[s], [5, 50, 95])])
            e["by_old_C_class"] = {}
            for a, b in CC:
                s = (o[m] >= a) & (o[m] < b)
                if s.sum(): e["by_old_C_class"][f"{a:g}-{min(b, 1):g}"] = dict(cells=int(s.sum()), p5_p50_p95=[float(v) for v in np.percentile(r[s], [5, 50, 95])])
        R[tag][str(D)] = e
    e = R[tag]["1.0"]; R[tag]["uniform_at_1km"] = bool(max(abs(e["rel_p5_p95"][0]), abs(e["rel_p5_p95"][1])) <= 0.03)
json.dump(R, open(f"{H}/grid_comparison.json", "w"), indent=1)
for tag in ("v3", "v3b"):
    print(tag, "uniform at 1 km:", R[tag]["uniform_at_1km"])
    for D in Dref:
        e = R[tag][str(D)]; print(f"  {D:5.2f} km  p5/50/95 {np.round(e['p5_p50_p95'], 4)}  rel {np.round(e['rel_p5_p95'], 4)}  area-w {np.round(e['area_weighted_p5_p50_p95'], 4)}  meanC {e['mean_C_old']:.4f}->{e['mean_C_new']:.4f}")
    for D in ("1.0",):
        print("   by px class", {k: np.round(v["p5_p50_p95"], 3).tolist() for k, v in R[tag][D]["by_finest_px_class"].items()})
        print("   by old C", {k: np.round(v["p5_p50_p95"], 3).tolist() for k, v in R[tag][D]["by_old_C_class"].items()})

# ---------------- figure: (a) map of new/old at 1 km (v3), (b) median and 5-95% band of new/old by reference diameter,
# (c) median new/old at 1 km by finest pixel-scale class of the cell
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": .6, "axes.edgecolor": "#5a5a57",
                     "xtick.color": "#5a5a57", "ytick.color": "#5a5a57", "text.color": "#1a1a19", "axes.labelcolor": "#1a1a19"})
INK, MUTED = "#1a1a19", "#6b6b67"
N = np.load(f"{H}/completeness_10ppd_v3.npz"); k1 = Dref.index(1.0)
o = O["C"][k1]; n = N["C"][k1]; rmap = np.where((nimg > 0) & (o > 0), n / np.where(o > 0, o, 1), np.nan)
fig = plt.figure(figsize=(7.1, 5.0))
a = fig.add_axes([0.07, 0.53, 0.60, 0.41]); ca = fig.add_axes([0.685, 0.53, 0.012, 0.41])
cmap = mpl.colormaps["RdBu_r"].copy(); cmap.set_bad("#e8e8e6"); norm = mpl.colors.TwoSlopeNorm(vmin=0.85, vcenter=1.0, vmax=1.15)
im = a.imshow(rmap, extent=[0, 360, -90, 90], origin="upper", cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
a.set_xticks([0, 90, 180, 270, 360]); a.set_yticks([-90, -45, 0, 45, 90]); a.set_xlabel("longitude (°E)"); a.set_ylabel("latitude (°N)")
a.set_title("a  new / old completeness, 1 km hollows (v3 refit model / v2 model); gray: not covered", loc="left", fontsize=7.5)
cb = fig.colorbar(im, cax=ca); cb.set_label("C$_{new}$ / C$_{old}$"); cb.outline.set_linewidth(.6)
b = fig.add_axes([0.07, 0.08, 0.40, 0.31])
for tag, col, off in (("v3", "#2a6fb0", -0.012), ("v3b", "#c7792b", 0.012)):
    x = np.array(Dref); q = np.array([R[tag][str(D)]["p5_p50_p95"] for D in Dref])
    b.vlines(x * (1 + off), q[:, 0], q[:, 2], color=col, lw=1.2, alpha=.8)
    b.plot(x * (1 + off), q[:, 1], "o", ms=3.5, color=col, mec="white", mew=.5, label=f"{tag} refit grid")
b.axhline(1, color=MUTED, lw=.6, ls=":"); b.set_xscale("log"); b.set_ylim(0.95, 1.25)
b.set_xticks([0.1, 0.2, 0.5, 1, 2, 5]); b.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"]); b.minorticks_off()
b.set_xlabel("reference diameter (km)"); b.set_ylabel("C$_{new}$ / C$_{old}$, cells")
b.set_title("b  median and 5-95% range over covered cells", loc="left", fontsize=7.5); b.legend(frameon=False, fontsize=6.5, loc="upper left")
b.text(5.0, 1.245, "5 km p95 = %.2f (v3)" % R["v3"]["5.0"]["p5_p50_p95"][2], fontsize=6, color=MUTED, ha="right", va="top")
c = fig.add_axes([0.57, 0.08, 0.40, 0.31])
cls = list(R["v3"]["1.0"]["by_finest_px_class"].keys())
for j, (tag, col, off) in enumerate((("v3", "#2a6fb0", -0.15), ("v3b", "#c7792b", 0.15))):
    q = np.array([R[tag]["1.0"]["by_finest_px_class"][k]["p5_p50_p95"] for k in cls])
    xx = np.arange(len(cls)) + off
    c.vlines(xx, q[:, 0], q[:, 2], color=col, lw=1.2, alpha=.8); c.plot(xx, q[:, 1], "o", ms=3.5, color=col, mec="white", mew=.5)
c.axhline(1, color=MUTED, lw=.6, ls=":"); c.set_xticks(range(len(cls))); c.set_xticklabels([k.replace(" m", "").replace("150.01", "150") for k in cls], fontsize=6.5)
c.set_xlabel("finest pixel scale in the cell (m)"); c.set_ylabel("C$_{new}$ / C$_{old}$, 1 km"); c.set_ylim(0.95, 1.25)
c.set_title("c  1 km ratio by finest pixel-scale class", loc="left", fontsize=7.5)
for ax in (a, b, c): ax.tick_params(length=2.5, pad=1.5)
fig.savefig(f"{H}/grid_comparison.png", dpi=200); print("wrote grid_comparison.png")
