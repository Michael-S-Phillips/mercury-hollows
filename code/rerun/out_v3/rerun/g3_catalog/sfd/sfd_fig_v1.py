# g3_catalog copy of figscripts/sfd_fig_v1.py (6 Nov 2026): only input/output paths and env settings changed (env G3_CAT, G3_THETA).
# NOTE: the power-law lines are hard-coded in the original at the v2 exponents (1.4-1.6 and 2.68); the points and bands are the rerun's.
"""Main-text figure: size distribution of hollows (4 Nov 2026). Inputs: analysis_workshop/completeness_sfd/sfd_grid.npz (detection
probability C_i(D) at each cataloged hollow's site, fit_sfd.py) and sfd_numbers.json (maximum-likelihood exponents, sfd_numbers.py).
Writes figures/fig_sfd_v1.png (400 dpi) and .pdf, 15 cm wide."""
import os, json, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
PR = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows"); W = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog/sfd/out_" + os.environ.get("G3_CAT", "v3") + "_" + os.environ.get("G3_THETA", "fixed"))
Z = np.load(f"{W}/sfd_grid.npz"); Dg, Ci, D = Z["Dg"], Z["Ci"].astype(float), Z["D"]; N = json.load(open(f"{W}/sfd_numbers.json")); lDg = np.log(Dg)
Cown = np.array([np.interp(np.log(x), lDg, r) for x, r in zip(D, Ci)]); w = 1 / np.clip(Cown, 1e-3, 1)
x = np.exp(np.linspace(np.log(0.15), np.log(12), 50))
obs = np.array([(D > v).sum() for v in x]); corr = np.array([w[D > v].sum() for v in x])
rng = np.random.default_rng(1); B = []
for _ in range(300):
    i = rng.integers(0, len(D), len(D)); B.append([w[i][D[i] > v].sum() for v in x])
lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
Cx = np.array([[np.interp(np.log(v), lDg, r) for v in x] for r in Ci[::2]]); med, q25, q75 = np.median(Cx, 0), np.percentile(Cx, 25, 0), np.percentile(Cx, 75, 0)
N1 = w[D > 1].sum(); N3 = w[D > 3].sum(); xx = np.exp(np.linspace(np.log(0.15), np.log(12), 80)); x3 = xx[xx >= 3]
BLUE, ORANGE, INK, INK2, GRID = "#2a78d6", "#eb6834", "#1a1a19", "#5a5a57", "#e4e3df"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": .6, "axes.edgecolor": INK2, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "axes.labelcolor": INK})
fig = plt.figure(figsize=(5.9, 4.9)); gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.42, wspace=0.3)
A = fig.add_subplot(gs[0, :]); S = fig.add_subplot(gs[1, 0]); C = fig.add_subplot(gs[1, 1])
for k in (A, S, C):
    k.set_xscale("log"); k.set_xlim(0.15, 12); k.grid(True, color=GRID, lw=0.5); k.set_axisbelow(True); k.set_xlabel("Hollow diameter D (km)")
    k.set_xticks([0.2, 0.5, 1, 2, 5, 10]); k.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10"]); k.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    for s in ("top", "right"): k.spines[s].set_visible(False)
A.set_yscale("log")
A.fill_between(xx, N1 * xx ** -1.4, N1 * xx ** -1.6, color=BLUE, alpha=0.12, lw=0)
A.plot(xx, N1 * xx ** -1.5, color=BLUE, lw=1.4, ls="--", label="Power law, exponent 1.4–1.6 (fits over 0.3–3 km)")
A.plot(x3, N3 * (x3 / 3) ** -2.68, color=BLUE, lw=1.4, ls=":", label="Exponent 2.7 (fit above 3 km)")
A.fill_between(x, lo, hi, color=BLUE, alpha=0.22, lw=0)
A.plot(x, corr, "o", ms=3.4, mfc=BLUE, mec="white", mew=0.6, label="Corrected for missed hollows (lower limit at small D)")
A.plot(x, obs, "o", ms=3.4, mfc=ORANGE, mec="white", mew=0.6, label="Cataloged hollows")
A.set_ylim(80, 2e5); A.set_ylabel("Number of hollows larger than D"); A.legend(frameon=False, fontsize=6.6, loc="upper right")
for (a_, b_), (s, l, h, n) in [((0.3, 0.7), N["slope_0.3_0.7"]), ((0.7, 3.0), N["slope_0.7_3"]), ((3.0, 12.0), N["slope_3_12"])]:
    S.fill_between([a_, b_], [l, l], [h, h], color=BLUE, alpha=0.22, lw=0); S.plot([a_, b_], [s, s], color=BLUE, lw=1.8)
s, l, h, n = N["sites_C03_ge_0.3"]["slope_03_1"]
S.fill_between([0.3, 1.0], [l, l], [h, h], color=ORANGE, alpha=0.2, lw=0); S.plot([0.3, 1.0], [s, s], color=ORANGE, lw=1.8)
S.text(0.32, 3.15, "best-imaged sites", color=ORANGE, fontsize=6.5); S.text(1.1, 1.05, "all sites", color=BLUE, fontsize=6.5)
S.axhline(0, color=INK2, lw=0.6); S.set_ylim(-0.2, 3.5); S.set_ylabel("Exponent of N(>D)")
C.fill_between(x, q25, q75, color=ORANGE, alpha=0.2, lw=0); C.plot(x, med, color=ORANGE, lw=1.6)
C.set_ylim(0, 1); C.set_ylabel("Detection probability")
for k, lab in zip((A, S, C), "abc"): k.text(-0.08 if k is A else -0.2, 1.04, lab, transform=k.transAxes, fontsize=9.5, fontweight="bold", va="bottom")
fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.09)
fig.savefig(f"{W}/fig_sfd_v1.png", dpi=200); print("wrote")
