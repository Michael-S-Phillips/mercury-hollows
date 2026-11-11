# g4 copy of code/rerun/out_v3/rerun/g2_insolation/threshold.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g4.py).
#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/insolation_v6/threshold.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Task (c), 2 Nov 2026: is the dependence of hosting on insolation better described as a threshold or a continuum?

Logistic models of h09 (>= 1 hollow within 0.9 radii) for craters >= 5 km within 60 deg of the equator (as in the paper)
and for all latitudes. Base controls: logD, abslat, C1 (completeness for 1 km hollows at the crater center).
Insolation terms compared by AIC:
  linear Delta T_w; step I(Delta T_w >= t) with t scanned 5-60 K (+2 AIC charged for choosing t);
  ramp min(Delta T_w, t) (linear rise then plateau, t scanned, +2 AIC); natural (restricted) cubic splines of
  Delta T_w with 3-5 df; mean annual temperature Tm (linear, spline) in place of Delta T_w; cos 2 lon.
Each also with the x logD interaction of the paper. Sensitivity: completeness as logit(C1) and a logD^2 term.
Tm comes from analysis_v9/crater_dag_nodes.csv. Writes threshold.json and threshold_shape.png.
"""
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common_null import *

c = pd.read_csv(f"{IN}/crater_features.csv", low_memory=False, dtype={"id": str})
n = pd.read_csv(f"{PR}/analysis_v9/crater_dag_nodes.csv", low_memory=False, dtype={"id": str})[["id", "Tm"]]
c = c.merge(n, on="id", how="left"); assert c.Tm.notna().all()
c["lc"] = np.log((c.C1 + 0.01) / (1.01 - c.C1)); c["logD2"] = c.logD ** 2
c["c2l"] = np.cos(np.radians(2 * c.lon)); c["Tm10"] = c.Tm / 10; c["dT10"] = c.dT / 10
OUT = {"description": __doc__}

def fit(f, d):
    return smf.logit(f, data=d).fit(disp=0, maxiter=200)

def run(d, ctrl, tag):
    R = {"n": int(len(d)), "hosts": int(d.h09.sum()), "controls": ctrl}
    A = {}
    A["base"] = fit(f"h09 ~ {ctrl}", d).aic
    A["linear_dT"] = fit(f"h09 ~ {ctrl} + dT10", d).aic
    A["linear_dT_x_logD"] = fit(f"h09 ~ {ctrl} + dT10 * logD", d).aic
    A["cos2lon"] = fit(f"h09 ~ {ctrl} + c2l", d).aic
    A["linear_Tm"] = fit(f"h09 ~ {ctrl} + Tm10", d).aic
    A["linear_Tm_x_logD"] = fit(f"h09 ~ {ctrl} + Tm10 * logD", d).aic
    for k in (3, 4, 5):
        A[f"spline_dT_df{k}"] = fit(f"h09 ~ {ctrl} + cr(dT, df={k})", d).aic
        A[f"spline_Tm_df{k}"] = fit(f"h09 ~ {ctrl} + cr(Tm, df={k})", d).aic
    ts = np.arange(5, 61, 1.0); step, ramp, stepx, rampx = [], [], [], []
    for t in ts:
        d["st"] = (d.dT >= t).astype(float); d["rp"] = np.minimum(d.dT, t) / 10
        step.append(fit(f"h09 ~ {ctrl} + st", d).aic); ramp.append(fit(f"h09 ~ {ctrl} + rp", d).aic)
        stepx.append(fit(f"h09 ~ {ctrl} + st * logD", d).aic); rampx.append(fit(f"h09 ~ {ctrl} + rp * logD", d).aic)
    step, ramp, stepx, rampx = map(np.array, (step, ramp, stepx, rampx))
    def best(a, nm):
        i = int(np.argmin(a)); ok = ts[a - a.min() <= 3.84]          # profile-likelihood 95% set for t (Delta deviance 3.84)
        A[nm] = float(a[i] + 2); R[nm + "_t"] = float(ts[i]); R[nm + "_t_95set"] = [float(ok.min()), float(ok.max())]
        return ts[i]
    tS = best(step, "step_dT"); tR = best(ramp, "ramp_dT"); best(stepx, "step_dT_x_logD"); best(rampx, "ramp_dT_x_logD")
    # Tm step scan (on the residual scale it is a latitude-longitude mix; reported for completeness)
    tms = np.arange(np.ceil(d.Tm.quantile(.05)), np.floor(d.Tm.quantile(.95)), 2.0); stm = []
    for t in tms:
        d["sm"] = (d.Tm >= t).astype(float); stm.append(fit(f"h09 ~ {ctrl} + sm", d).aic)
    stm = np.array(stm); A["step_Tm"] = float(stm.min() + 2); R["step_Tm_t"] = float(tms[np.argmin(stm)])
    # combined: step plus linear above (does dT matter above the threshold?)
    d["st"] = (d.dT >= tS).astype(float)
    m_sl = fit(f"h09 ~ {ctrl} + st + dT10", d); A["step_plus_linear_dT"] = float(m_sl.aic + 2)
    R["step_plus_linear"] = dict(step_OR=float(np.exp(m_sl.params["st"])), dT10_OR=float(np.exp(m_sl.params["dT10"])),
                                 dT10_p=float(m_sl.pvalues["dT10"]), step_p=float(m_sl.pvalues["st"]))
    ms = fit(f"h09 ~ {ctrl} + st", d); b, se = ms.params["st"], ms.bse["st"]
    R["step_OR"] = [float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))]
    R["frac_craters_above_step"] = float((d.dT >= tS).mean()); R["frac_hosts_above_step"] = float(d[d.h09 == 1].dT.ge(tS).mean())
    d["rp"] = np.minimum(d.dT, tR) / 10; mr = fit(f"h09 ~ {ctrl} + rp", d)
    R["ramp_OR_per10K_below_t"] = [float(np.exp(mr.params["rp"])), float(np.exp(mr.params["rp"] - 1.96 * mr.bse["rp"])), float(np.exp(mr.params["rp"] + 1.96 * mr.bse["rp"]))]
    ml = fit(f"h09 ~ {ctrl} + dT10", d)
    R["linear_OR_per10K"] = [float(np.exp(ml.params["dT10"])), float(np.exp(ml.params["dT10"] - 1.96 * ml.bse["dT10"])), float(np.exp(ml.params["dT10"] + 1.96 * ml.bse["dT10"]))]
    R["AIC"] = {k: round(float(v), 1) for k, v in A.items()}
    amin = min(A.values()); R["dAIC"] = {k: round(float(v - amin), 1) for k, v in sorted(A.items(), key=lambda kv: kv[1])}
    R["scan"] = dict(t=ts.tolist(), step=step.tolist(), ramp=ramp.tolist())
    # spline shape: odds ratio vs dT relative to dT = 0 (df=4), with pointwise 95% band
    m4 = fit(f"h09 ~ {ctrl} + cr(dT, df=4)", d)
    g = np.linspace(0, d.dT.max(), 63); ref = d.iloc[[0] * len(g)].copy(); ref["dT"] = g
    for k in ("logD", "abslat", "C1", "lc", "logD2"): ref[k] = d[k].median()
    spec = getattr(m4.model.data, "model_spec", None) or m4.model.data.design_info  # v0.3.1: statsmodels 0.14 names it design_info
    import patsy; X = np.asarray(patsy.build_design_matrices([spec], ref)[0])
    X0 = X - X[0]; eta = X0 @ m4.params.values; se_ = np.sqrt(np.einsum("ij,jk,ik->i", X0, m4.cov_params().values, X0))
    R["spline_shape"] = dict(dT=g.tolist(), OR=np.exp(eta).tolist(), lo=np.exp(eta - 1.96 * se_).tolist(), hi=np.exp(eta + 1.96 * se_).tolist())
    # empirical: host rate / C1-adjusted host rate in dT bins
    e = pd.cut(d.dT, [0, 10, 20, 30, 40, 50, 55, 60, 63])
    R["binned"] = d.groupby(e, observed=True).agg(n=("h09", "size"), hosts=("h09", "sum"), C1sum=("C1", "sum")).assign(
        hosts_per_C1sum=lambda x: x.hosts / x.C1sum).reset_index().astype({"dT": str}).to_dict("records")
    return R

samples = {"lat60": c[c.abslat <= 60].copy(), "global": c.copy()}
for nm, d in samples.items():
    for ctag, ctrl in (("main", "logD + abslat + C1"), ("sens_logitC_logD2", "logD + logD2 + abslat + lc")):
        R = run(d.copy(), ctrl, f"{nm}_{ctag}"); OUT[f"{nm}_{ctag}"] = R
        print(f"=== {nm} {ctag}: n={R['n']} hosts={R['hosts']}")
        print("   dAIC:", R["dAIC"])
        print(f"   step t={R['step_dT_t']} (95% set {R['step_dT_t_95set']}), OR {np.round(R['step_OR'], 2)}, frac craters above {R['frac_craters_above_step']:.2f}, hosts above {R['frac_hosts_above_step']:.2f}")
        print(f"   ramp t={R['ramp_dT_t']} (95% set {R['ramp_dT_t_95set']}), OR/10K below t {np.round(R['ramp_OR_per10K_below_t'], 2)}; linear OR/10K {np.round(R['linear_OR_per10K'], 2)}")
        print(f"   step+linear: {R['step_plus_linear']}   step x logD t={R['step_dT_x_logD_t']} ramp x logD t={R['ramp_dT_x_logD_t']} step Tm t={R['step_Tm_t']}")
        s = R["spline_shape"]; idx = [0, 10, 20, 30, 40, 50, 60]
        print("   spline OR vs dT:", [(round(s['dT'][i]), round(s['OR'][i], 2)) for i in idx])
        print("   binned:", [(r['dT'], r['n'], r['hosts'], round(r['hosts_per_C1sum'], 3)) for r in R["binned"]])
json.dump(OUT, open(f"{G2_OUT}/threshold.json", "w"), indent=1)

# figure: spline shape and step/ramp profiles (lat60 main)
R = OUT["lat60_main"]; s = R["spline_shape"]
fig, ax = plt.subplots(1, 2, figsize=(7.5, 3))
ax[0].fill_between(s["dT"], s["lo"], s["hi"], color="#1F6FB2", alpha=.2, lw=0); ax[0].plot(s["dT"], s["OR"], color="#1F6FB2")
ax[0].axvline(R["step_dT_t"], color="#D4541E", ls="--", lw=1); ax[0].set_yscale("log")
ax[0].set_xlabel("ΔT_w (K)"); ax[0].set_ylabel("odds ratio of hosting vs ΔT_w = 0"); ax[0].set_title("|lat| ≤ 60°, natural spline df 4", fontsize=8)
sc = R["scan"]; mn = min(min(sc["step"]), min(sc["ramp"]))
ax[1].plot(sc["t"], np.array(sc["step"]) - mn, label="step at t"); ax[1].plot(sc["t"], np.array(sc["ramp"]) - mn, label="ramp to t, then flat")
ax[1].axhline(R["AIC"]["linear_dT"] - mn, color="k", ls=":", label="linear"); ax[1].set_ylim(-1, 40)
ax[1].set_xlabel("threshold t (K)"); ax[1].set_ylabel("AIC − min (before +2 for t)"); ax[1].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{G2_OUT}/threshold_shape.png", dpi=200)
