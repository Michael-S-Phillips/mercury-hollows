#!/usr/bin/env python3
"""Decision figure (6 Nov 2026): does hosting rise with Delta T_w as a step or a steady rise?

Same sample and controls as threshold.py (lat60_main): craters >= 5 km within 60 deg of the equator, h09 = a hollow
within 0.9 radii, controls logD + abslat + C1. For v2, v3 and v3b:
  top: hosting probability for a crater at the median diameter, latitude and completeness, from a model with
       Delta T_w in bins (points, 95% intervals), and from the linear, best-step and step-plus-linear models (lines);
  bottom: AIC of a step at threshold t (+2 for choosing t) minus the AIC of the linear model, from threshold.json.
The Stage 1 step at 44 K (supplement S8) is marked for reference. Writes out/compare/fig_step_vs_linear.png and
step_vs_linear_record.json.
"""
import json, os, numpy as np, pandas as pd, patsy, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out/compare"; os.makedirs(OUT, exist_ok=True)
CTRL = "logD + abslat + C1"; BINS = [0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 63]
C_LIN, C_STEP, C_SL, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#8a8984"
fit = lambda f, d: smf.logit(f, data=d).fit(disp=0, maxiter=200)

def prob_curve(m, ref):
    X = np.asarray(patsy.build_design_matrices([m.model.data.model_spec], ref)[0])
    eta = X @ m.params.values; se = np.sqrt(np.einsum("ij,jk,ik->i", X, m.cov_params().values, X))
    ex = lambda z: 100 / (1 + np.exp(-z)); return ex(eta), ex(eta - 1.96 * se), ex(eta + 1.96 * se)

rec = {"description": __doc__}
fig, ax = plt.subplots(2, 3, figsize=(10.5, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[1.25, 1]))
for j, cat in enumerate(["v2", "v3", "v3b"]):
    c = pd.read_csv(f"{HERE}/inputs/{cat}/crater_features.csv", low_memory=False, usecols=["id", "dT", "abslat", "logD", "C1", "h09", "lat", "lon"])
    d = c[c.abslat <= 60].copy(); R = json.load(open(f"{HERE}/out/{cat}/threshold.json"))["lat60_main"]; t = R["step_dT_t"]
    d["dT10"] = d.dT / 10; d["st"] = (d.dT >= t).astype(float); d["bin"] = pd.cut(d.dT, BINS, include_lowest=True).astype(str)
    med = {k: float(d[k].median()) for k in ("logD", "abslat", "C1")}
    mb = fit(f"h09 ~ {CTRL} + C(bin)", d); ml = fit(f"h09 ~ {CTRL} + dT10", d)
    ms = fit(f"h09 ~ {CTRL} + st", d); msl = fit(f"h09 ~ {CTRL} + st + dT10", d)
    cats = sorted(d["bin"].unique(), key=lambda s: float(s.strip("([").split(",")[0]))
    ref = pd.DataFrame({**{k: [v] * len(cats) for k, v in med.items()}, "bin": cats})
    p, lo, hi = prob_curve(mb, ref); mid = np.array([(float(s.strip("([").split(",")[0]) + float(s.strip("])").split(",")[1])) / 2 for s in cats])
    g = np.linspace(0, 63, 631); rg = pd.DataFrame({**{k: [v] * len(g) for k, v in med.items()}, "dT10": g / 10, "st": (g >= t).astype(float)})
    a = ax[0, j]
    a.axvline(44, color=MUTED, ls=":", lw=1)
    a.errorbar(mid, p, yerr=[p - lo, hi - p], fmt="o", ms=4.5, color=INK, ecolor=INK, elinewidth=1, capsize=0, zorder=5, label="ΔT_w bins (95%)")
    a.plot(g, prob_curve(ml, rg)[0], color=C_LIN, lw=2, label="steady rise (linear)")
    a.plot(g, prob_curve(ms, rg)[0], color=C_STEP, lw=2, ls="--", label=f"step at {t:.0f} K")
    a.plot(g, prob_curve(msl, rg)[0], color=C_SL, lw=2, ls="-.", label="step + linear")
    aic = R["AIC"]; dl = lambda k: aic[k] - aic["linear_dT"]
    a.text(0.03, 0.97, f"{cat}: {R['hosts']} hosts of {R['n']:,}\nΔAIC vs linear:\n  step {dl('step_dT'):+.1f}\n  step + linear {dl('step_plus_linear_dT'):+.1f}",
           transform=a.transAxes, va="top", fontsize=7.5, color=INK)
    a.set_ylim(0, 1.08 * max(hi.max(), 0.01)); a.set_title(f"{cat} catalog", fontsize=10, color=INK)
    if j == 0: a.set_ylabel(f"hosting probability (%)\nmedian crater, D = {10**med['logD']:.1f} km", fontsize=8.5)
    b = ax[1, j]; sc = R["scan"]; ts = np.array(sc["t"]); dA = np.array(sc["step"]) + 2 - aic["linear_dT"]
    b.axhspan(-10, 0, color=C_STEP, alpha=.08, lw=0); b.axhline(0, color=C_LIN, lw=2); b.axvline(44, color=MUTED, ls=":", lw=1)
    b.plot(ts, dA, color=C_STEP, lw=2, marker="o", ms=3)
    b.annotate(f"best {t:.0f} K ({dA.min():+.1f})", (t, dA.min()), (t + 4, 22), fontsize=7.5, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=.8))
    b.text(1, 2, "linear model", color=INK, fontsize=7.5, ha="left", va="bottom"); b.text(1, -2, "step fits better", color=INK, fontsize=7.5, ha="left", va="top")
    b.set_ylim(-12, 100); b.set_xlabel("ΔT_w (K)", fontsize=8.5)
    if j == 0: b.set_ylabel("AIC(step at t) − AIC(linear)", fontsize=8.5)
    # where the step location comes from: craters and hosts just below / above t
    w = {f"[{t-2:.0f},{t:.0f})": d[(d.dT >= t - 2) & (d.dT < t)], f"[{t:.0f},{t+2:.0f})": d[(d.dT >= t) & (d.dT < t + 2)]}
    rec[cat] = dict(t=t, hosts=R["hosts"], n=R["n"], dAIC_step=dl("step_dT"), dAIC_step_plus_linear=dl("step_plus_linear_dT"),
                    step_plus_linear=R["step_plus_linear"], bins=[dict(bin=s, p=float(x), lo=float(y), hi=float(z)) for s, x, y, z in zip(cats, p, lo, hi)],
                    near_t={k: dict(n=int(len(v)), hosts=int(v.h09.sum()), lat_med=float(v.lat.median()), lon_q=[float(q) for q in v.lon.quantile([.1, .5, .9])]) for k, v in w.items()},
                    median_controls=med)
    print(cat, "t", t, "near t:", rec[cat]["near_t"])
    print("   bins", [(r["bin"], round(r["p"], 2)) for r in rec[cat]["bins"]])
for a in ax.flat:
    a.tick_params(labelsize=8, colors=INK); [s.set_color(MUTED) for s in a.spines.values()]; a.spines[["top", "right"]].set_visible(False)
    a.grid(axis="y", color="#e6e5e0", lw=.6); a.set_axisbelow(True)
h, l = ax[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, frameon=False)
ax[0, 0].text(44.5, ax[0, 0].get_ylim()[1] * .99, "Stage 1\nstep (44 K)", fontsize=7, color=MUTED, va="top")
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(f"{OUT}/fig_step_vs_linear.png", dpi=200); fig.savefig(f"{OUT}/fig_step_vs_linear.pdf")
json.dump(rec, open(f"{OUT}/step_vs_linear_record.json", "w"), indent=1)
