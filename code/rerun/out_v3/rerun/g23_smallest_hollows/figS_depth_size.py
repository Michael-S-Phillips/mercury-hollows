#!/usr/bin/env python3
"""g23 (9 Nov 2026), Michael's v22 comment: scatter plot of depth against hollow size for S9 (supplement v24).
(a) Hand measurements (bench_depths.csv), corrected to 5 m/pixel and 70 deg incidence with the coefficients of the hand fit
    (ln h ~ ln D + ln px + ln tan i; bench_summary.json): measured shadows (filled; uncertain open), the shadow that fills the hollow
    (lower limit, up triangle) and no measurable shadow (upper limit, down triangle). Line: the fitted power law.
(b) Automatic depths (depth_paper_rows.csv, 1,558), corrected to 10 m/pixel and 70 deg with the coefficients of the fit to them
    (fig_depth_v2_numbers.json 'coef'), and the 320 boxes without a shadow component as upper limits (depth_all.csv, n_runs == 0, at
    h_upper). Medians in diameter bins with 95% intervals from image bootstraps; solid line, least squares (exponent 0.61); dashed lines,
    exponents corrected for the scatter of measured diameters with half and all of the repeat-detection variance (g26 me.json).
Writes figS_depth_size_v1.png (400 dpi) and .pdf, 15 cm wide."""
import os, json, numpy as np, pandas as pd, matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); G26 = os.path.join(os.path.dirname(HERE), "g26_measurement_error")
INK, MUTED, GRID, ORANGE, GRAY = "#1a1a19", "#6b6b67", "#e4e3df", "#eb6834", "#4a4a47"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.8, "axes.linewidth": .6, "axes.edgecolor": "#5a5a57", "xtick.labelsize": 6.2, "ytick.labelsize": 6.2})
S = json.load(open(f"{HERE}/bench_summary.json"))["ols_all_measured_hollows"]; bD, bp, bt = S["lnD"][0], S["lnpx"][0], S["lntan_i"][0]
b = pd.read_csv(f"{HERE}/bench_depths.csv"); b = b[b.cls == "hollow"].copy()
corr = lambda h, px, inc, p0, i0, cp, ct: h * (p0 / px) ** cp * (np.tan(np.radians(i0)) / np.tan(np.radians(inc))) ** ct
b["hc"] = corr(b.depth_m, b.px, b.inc, 5.0, 70.0, bp, bt); b["hu"] = corr(b.h_up, b.px, b.inc, 5.0, 70.0, bp, bt)
fig, (e1, e2) = plt.subplots(1, 2, figsize=(5.9, 2.75))
for k in (e1, e2):
    k.set_xscale("log"); k.set_yscale("log"); k.grid(True, color=GRID, lw=.5); k.set_axisbelow(True)
    for s_ in ("top", "right"): k.spines[s_].set_visible(False)
m = b[(b.shadow == "measured") & (b.depth_m > 0)]; cu = m.uncertain.astype(bool)
e1.scatter(m.D_m[~cu] / 1000, m.hc[~cu], s=9, c=GRAY, lw=0, label="measured")
e1.scatter(m.D_m[cu] / 1000, m.hc[cu], s=9, facecolor="none", edgecolor=GRAY, lw=.6, label="measured, uncertain")
lo_ = b[b.shadow == "full"]; up_ = b[b.shadow == "none"]
e1.scatter(lo_.D_m / 1000, lo_.hc, marker="^", s=14, c=INK, lw=0, label="lower limit")
e1.scatter(up_.D_m / 1000, up_.hu, marker="v", s=14, c=INK, lw=0, label="upper limit")
xx = np.exp(np.linspace(np.log(0.035), np.log(1.1), 50)); x0 = np.exp(np.median(np.log(m.D_m / 1000))); y0 = np.exp(np.median(np.log(m.hc)))
e1.plot(xx, y0 * (xx / x0) ** bD, color=INK, lw=1.0, label=f"depth ∝ D$^{{{bD:.2f}}}$")
e1.set_xlim(0.03, 1.2); e1.set_ylim(1, 100); e1.set_xlabel("Hollow diameter (km)"); e1.set_ylabel("Depth, corrected to 5 m/pixel and 70° (m)")
e1.legend(frameon=False, fontsize=5.8, loc="upper left", handletextpad=.3, borderpad=.2)
cf = json.load(open(f"{HERE}/fig_depth_v2_numbers.json"))["coef"]; ap, at, aD = cf["lnpx"], cf["ltan"], cf["lnD"]
a = pd.read_csv(f"{HERE}/depth_paper_rows.csv"); a["hc"] = corr(a.depth_m, a.px, a.inc, 10.0, 70.0, ap, at); a["D"] = a.D_can_m / 1000; a["img"] = a.img.astype(str)
al = pd.read_csv(f"{HERE}/depth_all.csv"); fl = al["flags"].fillna("")
u = al[fl.str.contains("no_shadow") & (al.n_runs == 0) & (al.h_upper_m > 0) & ~fl.str.contains("edge_sun|frame_edge|no_edr")].copy()
u["hu"] = corr(u.h_upper_m, u.px, u.inc, 10.0, 70.0, ap, at)
e2.scatter(a.D, a.hc, s=3, c=ORANGE, alpha=.35, lw=0, label=f"automatic ({len(a):,})", rasterized=True)
e2.scatter(u.D_can_m / 1000, u.hu, marker="v", s=5, c=MUTED, alpha=.45, lw=0, label=f"upper limits, no shadow ({len(u)})", rasterized=True)
rng = np.random.default_rng(5); edges = np.exp(np.linspace(np.log(0.1), np.log(4.0), 9)); rows = []
for lo, hi in zip(edges[:-1], edges[1:]):
    s = a[(a.D >= lo) & (a.D < hi)]
    if len(s) < 15: continue
    ims = s.img.unique(); g = {k_: v.hc.values for k_, v in s.groupby("img")}
    bs = [np.median(np.concatenate([g[i] for i in rng.choice(ims, len(ims))])) for _ in range(500)]
    rows.append((np.exp(np.median(np.log(s.D))), np.median(s.hc), *np.percentile(bs, [2.5, 97.5])))
r = np.array(rows); e2.errorbar(r[:, 0], r[:, 1], yerr=[r[:, 1] - r[:, 2], r[:, 3] - r[:, 1]], fmt="s", ms=3.6, color=INK, mfc=INK, lw=.8, capsize=1.6, zorder=5, label="bin medians")
x0 = np.exp(np.median(np.log(a.D))); y0 = np.exp(np.median(np.log(a.hc))); xx = np.exp(np.linspace(np.log(0.08), np.log(5), 50))
e2.plot(xx, y0 * (xx / x0) ** aD, color=INK, lw=1.0, label=f"least squares, {aD:.2f}")
me = json.load(open(f"{G26}/me.json"))["C_diameter_error"]["paper_quantities"]
for tag, ls in [("half", (0, (4, 2))), ("full", (0, (1.5, 1.5)))]:
    ex = me[tag]["size_exponent"][1]; e2.plot(xx, y0 * (xx / x0) ** ex, color=INK, lw=0.9, ls=ls, label=f"corrected, {ex:.2f}")
e2.set_xlim(0.08, 5); e2.set_ylim(1, 1500); e2.set_xlabel("Hollow diameter (km)"); e2.set_ylabel("Depth, corrected to 10 m/pixel and 70° (m)")
e2.legend(frameon=False, fontsize=5.4, loc="upper left", ncol=2, handletextpad=.3, borderpad=.2, columnspacing=0.8, handlelength=1.8)
for k, tk in [(e1, [0.05, 0.1, 0.2, 0.5, 1]), (e2, [0.1, 0.2, 0.5, 1, 2, 5])]:
    k.set_xticks(tk); k.set_xticklabels([f"{v:g}" for v in tk]); k.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    k.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
for k, lab in zip((e1, e2), "ab"): k.text(-0.16, 1.02, lab, transform=k.transAxes, fontsize=9, fontweight="bold", va="bottom")
fig.subplots_adjust(left=0.09, right=0.99, top=0.93, bottom=0.15, wspace=0.3)
fig.savefig(f"{HERE}/figS_depth_size_v1.png", dpi=400); fig.savefig(f"{HERE}/figS_depth_size_v1.pdf")
print(dict(hand_exp=round(bD, 3), auto_exp=aD, n_auto=len(a), n_upper=len(u), bins=r.round(2).tolist()))
