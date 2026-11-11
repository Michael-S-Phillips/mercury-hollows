# v3 (9 Nov 2026, Michael's v6.12 comment): panel (a) legend moved to the lower right, where it is not mistaken for data, and the
# sulfur-front label moved left of the line. Otherwise as v2. Output fig_depth_v3.*
#!/usr/bin/env python3
"""g23 (9 Nov 2026): main-text depth figure for v6.11 (replaces fig4_depth of v6.10). Blue: shadow-length depths of Blewett et al. (2016);
orange: our automatic depths (depth_paper_rows.csv, 1,558 hollows); dark gray: Michael's hand depths (bench_depths.csv).
(a) Median depth in bands of absolute latitude against band-median mean annual temperature: Blewett corrected to their median pixel
    scale and incidence (depth_reanalysis_record.json 'hc'); ours corrected to px = 10 m, i = 70 deg and D = 0.5 km with the coefficients
    of ln h ~ ln px + ln tan i + ln D fitted to our depths; 95% intervals from 1,000 bootstrap resamples of images. Dashed: sulfur front
    limited by vapor diffusion (-4,514 K), anchored at Blewett's 0-30 deg band.
(b) Depth against hollow diameter: our depths in diameter bins (median, corrected to px = 10 m and i = 70 deg; image bootstrap), the
    fitted power law, and the hand depths corrected to the same reference with the same coefficients.
(c) Median depth per host crater (>= 3 hollows) against host diameter: Blewett as in v6.10 (corrected to median px, incidence and |lat|),
    ours corrected also to D = 0.5 km; lines are the regressions on log host diameter.
Writes fig_depth_v3.png (400 dpi) and .pdf, 15 cm wide, and fig_depth_v3_numbers.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
HE = os.path.expanduser("~/mnt/hollows-emissivity"); INSOL = f"{HE}/manuscripts/insolation_thermal"
INK, MUTED, GRID, BLUE, ORANGE = "#1a1a19", "#6b6b67", "#e4e3df", "#2a78d6", "#eb6834"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.8, "axes.linewidth": .6, "axes.edgecolor": "#5a5a57", "xtick.labelsize": 6.2, "ytick.labelsize": 6.2})
a = pd.read_csv("depth_paper_rows.csv"); a["img"] = a.img.astype(str)
m = smf.ols("lnh ~ lnpx + ltan + lnD", a).fit(); bp, bt, bD = m.params["lnpx"], m.params["ltan"], m.params["lnD"]
PX0, LT0, LD0 = np.log(10.0), np.log(np.tan(np.radians(70.0))), np.log(0.5)
a["lhc"] = a.lnh - bp * (a.lnpx - PX0) - bt * (a.ltan - LT0) - bD * (a.lnD - LD0)          # px, i and size corrected
a["lhs"] = a.lnh - bp * (a.lnpx - PX0) - bt * (a.ltan - LT0)                             # px and i corrected (for panel b)
rng = np.random.default_rng(5); imgs = a.img.unique(); G = {k: v for k, v in a.groupby("img")}
def boot_median(sel, col, nb=1000):
    s = a[sel]; ims = s.img.unique(); g = {k: v[col].values for k, v in s.groupby("img")}; meds = []
    for _ in range(nb):
        pick = rng.choice(ims, len(ims)); meds.append(np.median(np.concatenate([g[i] for i in pick])))
    return float(np.exp(np.median(s[col]))), float(np.exp(np.percentile(meds, 2.5))), float(np.exp(np.percentile(meds, 97.5))), int(sel.sum())
NUM = {"coef": dict(lnpx=round(bp, 3), ltan=round(bt, 3), lnD=round(bD, 3))}
fig, (e1, e2, e3) = plt.subplots(1, 3, figsize=(5.9, 2.35))
# (a)
DR = json.load(open(f"{INSOL}/records/depth_reanalysis_record.json")); bk = ["0-30", "30-45", "45-60", "60-75", "75-90"]
T_ = np.array([DR["bands"][k]["T_med"] for k in bk]); y = np.array([DR["bands"][k]["hc"] for k in bk])
e1.errorbar(T_ - 0.8, y[:, 0], yerr=[y[:, 0] - y[:, 1], y[:, 2] - y[:, 0]], fmt="o", ms=3.8, color=BLUE, mfc="white", mew=1.0, lw=.8, capsize=1.8, zorder=5, label="Blewett et al. (2016)")
ours = []
for lo, hi in [(0, 30), (30, 45), (45, 60), (60, 75), (75, 90)]:
    sel = ((a.abslat >= lo) & (a.abslat < hi)).values
    if sel.sum() >= 15: ours.append((float(np.median(a.loc[sel, "T"])),) + boot_median(sel, "lhc") + (f"{lo}-{hi}",))
ours = np.array([o[:5] for o in ours], float); NUM["a_ours"] = ours.round(2).tolist()
e1.errorbar(ours[:, 0] + 0.8, ours[:, 1], yerr=[ours[:, 1] - ours[:, 2], ours[:, 3] - ours[:, 1]], fmt="o", ms=3.8, color=ORANGE, mfc=ORANGE, mew=1.0, lw=.8, capsize=1.8, zorder=5, label="this study")
Tg = np.linspace(236, 316, 100); T0, h0 = DR["bands"]["0-30"]["T_med"], DR["bands"]["0-30"]["hc"][0]
e1.plot(Tg, h0 * np.exp(-9028.6679 / 2 * (1 / Tg - 1 / T0)), ls=(0, (5, 2)), color=INK, lw=1.0, zorder=3)
e1.text(237, 3.0, "sulfur front,\ndiffusion-limited", fontsize=5.9, color=INK, ha="left", va="center")  # v3: label left of the line, clear of the legend
e1.set_yscale("log"); e1.set_ylim(.3, 120); e1.set_xlim(234, 318)
e1.set_yticks([.5, 1, 2, 5, 10, 20, 50, 100]); e1.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50", "100"]); e1.minorticks_off()
e1.set_xlabel("Mean annual surface temperature (K)"); e1.set_ylabel("Median hollow depth (m)")
e1.legend(frameon=False, fontsize=5.9, loc="lower right", handletextpad=.3, borderpad=.2)
# (b)
bins = [0.1, 0.2, 0.3, 0.5, 1.0]; rows = []   # the automatic method is validated to 1 km (hand depths); larger boxes mostly hold several hollows
mb_ = smf.ols("lnh ~ lnpx + ltan + lnD", a[a.D_can_m <= 1000]).fit(); NUM["b_fit_le1km"] = dict(lnD=round(mb_.params["lnD"], 3), lnpx=round(mb_.params["lnpx"], 3), ltan=round(mb_.params["ltan"], 3))
bp, bt, bD = mb_.params["lnpx"], mb_.params["ltan"], mb_.params["lnD"]; a["lhs"] = a.lnh - bp * (a.lnpx - PX0) - bt * (a.ltan - LT0); m = mb_
for lo, hi in zip(bins[:-1], bins[1:]):
    sel = ((a.D_can_m / 1000 >= lo) & (a.D_can_m / 1000 < hi)).values
    if sel.sum() >= 10: rows.append((float(np.exp(np.median(np.log(a.loc[sel, "D_can_m"] / 1000)))),) + boot_median(sel, "lhs"))
rows = np.array(rows); NUM["b_bins"] = rows.round(3).tolist()
hb = pd.read_csv("bench_depths.csv"); hb = hb[(hb.shadow == "measured") & (hb.depth_m > 0) & ~hb.cls.isin(["crater", "notdep"])].copy()
hb["lhs"] = np.log(hb.depth_m) - bp * (np.log(hb.px) - PX0) - bt * (np.log(np.tan(np.radians(hb.inc))) - LT0)
e2.scatter(hb.D_m / 1000, np.exp(hb.lhs), s=5, c=INK, alpha=.55, lw=0, zorder=3, label="hand-measured")
e2.errorbar(rows[:, 0], rows[:, 1], yerr=[rows[:, 1] - rows[:, 2], rows[:, 3] - rows[:, 1]], fmt="s", ms=3.8, color=ORANGE, mfc=ORANGE, lw=.8, capsize=1.8, zorder=5, label="automatic, binned")
Dg = np.logspace(np.log10(0.04), np.log10(1.0), 50); e2.plot(Dg, np.exp(m.params["Intercept"] + bp * PX0 + bt * LT0 + bD * np.log(Dg)), color=ORANGE, lw=1.0, zorder=4)
e2.text(.04, .96, f"depth ∝ D$^{{{bD:.2f}}}$", transform=e2.transAxes, fontsize=6.3, va="top")
e2.set_xscale("log"); e2.set_yscale("log"); e2.set_xlim(.03, 1.3); e2.set_ylim(1, 100)
e2.set_xticks([.05, .1, .2, .5, 1]); e2.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1"]); e2.set_yticks([1, 2, 5, 10, 20, 50, 100]); e2.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"]); e2.minorticks_off()
e2.set_xlabel("Hollow diameter (km)"); e2.set_ylabel("Depth, corrected (m)"); e2.legend(frameon=False, fontsize=5.9, loc="lower right", handletextpad=.3, borderpad=.2)
# (c) Blewett per-host medians (as v6.10) and ours
d = pd.read_csv(f"{INSOL}/data/test1_rows_inc.csv"); d["ltan"] = np.log(np.tan(np.radians(d.inc)))
cr = pd.read_csv(f"{INSOL}/data/test2_craters_h18.csv"); cr = cr[cr.diameter >= 5]
def xyz(la, lo): la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
tree = cKDTree(xyz(cr.latitude.values, cr.lon360.values)); P = xyz(d.latitude.values, d.longitude.values % 360); Dkm = cr.diameter.values; host = np.full(len(d), -1)
for i, cl in enumerate(tree.query_ball_point(P, r=0.36)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(P[i] @ tree.data[cl].T, -1, 1)) * 2439.7; ok = cl[ang <= 1.1 * Dkm[cl] / 2]
    if len(ok): host[i] = ok[np.argmin(Dkm[ok])]
d["host"] = host; d["hostD"] = np.where(host >= 0, Dkm[np.maximum(host, 0)], np.nan); dd = d[d.host >= 0].copy()
dd["ld"] = np.log10(dd.depth); dd["lp"] = np.log10(dd.px); dd["li"] = dd.ltan / np.log(10)
def per_host(df, key, extra):
    agg = dict(n=("ld", "size"), D=("hostD", "first"), ld=("ld", "median"), lp=("lp", "median"), li=("li", "median"), al=("abslat", "median"), **{k: (k, "median") for k in extra})
    pc = df.groupby(key).agg(**agg); pc = pc[pc.n >= 3].copy(); pc["logD"] = np.log10(pc.D)
    mm = smf.ols("ld ~ logD + lp + li + al" + "".join(" + " + k for k in extra), pc).fit(); pc["ldc"] = pc.ld
    for k in ["lp", "li", "al"] + extra: pc["ldc"] -= mm.params[k] * (pc[k] - pc[k].median())
    return pc, mm
pcb, mb = per_host(dd, "host", [])
a2 = a.dropna(subset=["hostD"]).copy(); a2["ld"] = a2.lnh / np.log(10); a2["lp"] = a2.lnpx / np.log(10); a2["li"] = a2.ltan / np.log(10); a2["lD"] = a2.lnD / np.log(10)
pco, mo = per_host(a2, "host_uid", ["lD"])
for pc, mm, col, fc, lab in [(pcb, mb, BLUE, "white", "Blewett et al. (2016)"), (pco, mo, ORANGE, ORANGE, "this study")]:
    e3.scatter(pc.D, 10 ** pc.ldc, s=np.clip(pc.n, 3, 80) * .5 + 3, facecolor=fc, edgecolor=col, lw=.7, alpha=.85, zorder=4, label=lab)
    Dg = np.logspace(np.log10(pc.D.min() * .9), np.log10(pc.D.max() * 1.1), 50)
    e3.plot(Dg, 10 ** (pc.ldc.median() + mm.params["logD"] * (np.log10(Dg) - pc.logD.median())), color=col, lw=1.0, zorder=5)
NUM["c"] = {"blewett": dict(hosts=int(len(pcb)), slope=round(mb.params["logD"], 3), se=round(mb.bse["logD"], 3)), "ours": dict(hosts=int(len(pco)), slope=round(mo.params["logD"], 3), se=round(mo.bse["logD"], 3))}
e3.set_xscale("log"); e3.set_yscale("log"); e3.set_ylim(3, 150)
e3.set_xticks([5, 10, 30, 100, 300, 1000]); e3.set_xticklabels(["5", "10", "30", "100", "300", "1000"]); e3.set_yticks([5, 10, 20, 50, 100]); e3.set_yticklabels(["5", "10", "20", "50", "100"]); e3.minorticks_off()
e3.set_xlabel("Host crater diameter (km)"); e3.set_ylabel("Median depth per crater (m)"); e3.legend(frameon=False, fontsize=5.9, loc="lower right", handletextpad=.3, borderpad=.2)
for k, lab in zip((e1, e2, e3), "abc"):
    k.grid(True, color=GRID, lw=.4); k.set_axisbelow(True); [k.spines[s].set_visible(False) for s in ("top", "right")]
    k.text(-0.2, 1.03, lab, transform=k.transAxes, fontsize=9, fontweight="bold", va="bottom")
fig.subplots_adjust(left=.075, right=.99, top=.9, bottom=.18, wspace=.42)
fig.savefig("fig_depth_v3.png", dpi=400); fig.savefig("fig_depth_v3.pdf"); json.dump(NUM, open("fig_depth_v3_numbers.json", "w"), indent=1); print(json.dumps(NUM))
