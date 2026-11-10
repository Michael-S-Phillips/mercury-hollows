# v3 (9 Nov 2026, Michael): two panels. The panel of fitted exponents (old b) is dropped; the exponents are given in the text and
# in Section S2. Old c (detection probability) becomes b, beside a. Output fig_sfd_v3_v3cat.
# v2 (9 Nov 2026, Michael's v6.10 comment 5): shaded regions between the cataloged and corrected counts (hollows the detector
# missed) and between the corrected counts and the power law (missing if the power law continues to small sizes; shortfall at
# large sizes, a possible upper limit on hollow size); panel b y label restored (it sat inside a comment). Output fig_sfd_v2_v3cat.
# g5 copy of manuscripts/PlanetaryResearch_hollows/figscripts/sfd_fig_v1.py (6 Nov 2026). Changed: inputs read from g3_catalog/sfd/out_v3_refit (v3 catalog, refit detection model); the power-law lines use the fitted exponents from sfd_numbers.json instead of the hard-coded 1.4-1.6 and 2.68; panel b y-limit raised from 3.5 to 4.0 so the 3-12 km interval (2.91-3.62) is not clipped; output to g5_figures/fig_sfd_v1_v3cat.png/.pdf.
"""Main-text figure: size distribution of hollows (4 Nov 2026). Inputs: analysis_workshop/completeness_sfd/sfd_grid.npz (detection
probability C_i(D) at each cataloged hollow's site, fit_sfd.py) and sfd_numbers.json (maximum-likelihood exponents, sfd_numbers.py).
Writes figures/fig_sfd_v1.png (400 dpi) and .pdf, 15 cm wide."""
import os, json, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
PR = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows"); W = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog/sfd/out_v3_refit")  # g5: v3 catalog, refit detection model
G5 = os.path.dirname(os.path.abspath(__file__))  # g5: output folder
Z = np.load(f"{W}/sfd_grid.npz"); Dg, Ci, D = Z["Dg"], Z["Ci"].astype(float), Z["D"]; N = json.load(open(f"{W}/sfd_numbers.json")); lDg = np.log(Dg)
Cown = np.array([np.interp(np.log(x), lDg, r) for x, r in zip(D, Ci)]); w = 1 / np.clip(Cown, 1e-3, 1)
x = np.exp(np.linspace(np.log(0.15), np.log(12), 50))
obs = np.array([(D > v).sum() for v in x]); corr = np.array([w[D > v].sum() for v in x])
rng = np.random.default_rng(1); B = []
for _ in range(300):
    i = rng.integers(0, len(D), len(D)); B.append([w[i][D[i] > v].sum() for v in x])
lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
Cx = np.array([[np.interp(np.log(v), lDg, r) for v in x] for r in Ci[::2]]); med, q25, q75 = np.median(Cx, 0), np.percentile(Cx, 25, 0), np.percentile(Cx, 75, 0)
N1 = w[D > 1].sum(); N3 = w[D > 3].sum()
# g5: exponent lines from the fitted values in sfd_numbers.json (were typed in as 1.4, 1.5, 1.6 and 2.68). The band uses the
# min and max of the 0.3-0.7, 0.7-3 and 1-3 km exponents rounded to one decimal, as the paper states the range (sfd_derived.py).
_a3 = [N["slope_0.3_0.7"][0], N["slope_0.7_3"][0], N["slope_1_3"][0]]; A_LO, A_HI = round(min(_a3), 1), round(max(_a3), 1); A_MID = (A_LO + A_HI) / 2; A_STEEP = N["slope_3_12"][0]
xx = np.exp(np.linspace(np.log(0.15), np.log(12), 80)); x3 = xx[xx >= 3]
BLUE, ORANGE, INK, INK2, GRID = "#2a78d6", "#eb6834", "#1a1a19", "#5a5a57", "#e4e3df"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": .6, "axes.edgecolor": INK2, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "axes.labelcolor": INK})
fig = plt.figure(figsize=(5.9, 3.95)); gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=0.28)
A = fig.add_subplot(gs[0, 0]); C = fig.add_subplot(gs[0, 1])
for k in (A, C):
    k.set_xscale("log"); k.set_xlim(0.15, 12); k.grid(True, color=GRID, lw=0.5); k.set_axisbelow(True); k.set_xlabel("Hollow diameter D (km)")
    k.set_xticks([0.2, 0.5, 1, 2, 5, 10]); k.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10"]); k.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    for s in ("top", "right"): k.spines[s].set_visible(False)
A.set_yscale("log")
A.fill_between(xx, N1 * xx ** -A_LO, N1 * xx ** -A_HI, color=BLUE, alpha=0.12, lw=0)
A.plot(xx, N1 * xx ** -A_MID, color=BLUE, lw=1.4, ls="--", label=f"Power law, exponent {A_LO:.1f}–{A_HI:.1f} (fits over 0.3–3 km)")
A.plot(x3, N3 * (x3 / 3) ** -A_STEEP, color=BLUE, lw=1.4, ls=":", label=f"Exponent {A_STEEP:.1f} (fit above 3 km)")
A.fill_between(x, lo, hi, color=BLUE, alpha=0.22, lw=0)
pl = N1 * x ** -A_MID; mpl.rcParams["hatch.linewidth"] = 0.5
A.fill_between(x, obs, corr, where=corr > obs, color=ORANGE, alpha=0.14, lw=0, interpolate=True, label="Missed by the detector (estimated)")
ms_, ml_ = x <= 1.0, x >= 3.0
A.fill_between(x[ms_], corr[ms_], pl[ms_], where=pl[ms_] > corr[ms_], facecolor="none", edgecolor=BLUE, hatch="////", lw=0, interpolate=True, label="Missing if the power law continues to small sizes")
A.fill_between(x[ml_], corr[ml_], pl[ml_], where=pl[ml_] > corr[ml_], facecolor="none", edgecolor=INK2, hatch="\\\\\\\\", lw=0, interpolate=True, label="Shortfall from the power law: a possible upper size limit")
A.plot(x, corr, "o", ms=3.4, mfc=BLUE, mec="white", mew=0.6, label="Corrected for missed hollows (lower limit at small D)")
A.plot(x, obs, "o", ms=3.4, mfc=ORANGE, mec="white", mew=0.6, label="Cataloged hollows")
A.set_ylim(80, 2e5); A.set_ylabel("Number of hollows larger than D"); hh, ll = A.get_legend_handles_labels(); order = [6, 5, 0, 1, 2, 3, 4]
fig.legend([hh[i] for i in order], [ll[i] for i in order], frameon=False, fontsize=6.4, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0), handlelength=2.4, columnspacing=1.6)
C.fill_between(x, q25, q75, color=ORANGE, alpha=0.2, lw=0); C.plot(x, med, color=ORANGE, lw=1.6)
C.set_ylim(0, 1); C.set_ylabel("Detection probability")
for k, lab in zip((A, C), "ab"): k.text(-0.13 if k is A else -0.22, 1.02, lab, transform=k.transAxes, fontsize=9.5, fontweight="bold", va="bottom")
fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.255)
fig.savefig(f"{G5}/fig_sfd_v3_v3cat.png", dpi=400); fig.savefig(f"{G5}/fig_sfd_v3_v3cat.pdf")  # g5: output
print("wrote", dict(exp_band=[A_LO, A_HI], exp_above_3km=A_STEEP, N1=round(N1), N3=round(N3), n_hollows=len(D), obs_gt_12km=int((D > 12).sum()), obs_min_plotted=int(obs.min())))
