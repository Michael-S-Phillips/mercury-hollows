#!/usr/bin/env python3
"""g26 (9 Nov 2026), for Michael: is the measurement uncertainty of the automatic shadow depths propagated into the paper's depth results?
(A) Differential error against Michael's hand depths: ln(hand/auto) regressed on hollow size, pixel scale, incidence, latitude, and the
    multi flag (g23 bench, <= 7 m/px), SEs clustered by image.
(B) Error of each method from the three measurements of the same 100 boxes in g24 (Blewett max, Michael max, pilot automatic):
    three-cornered hat on ln depth, with plain and robust (MAD) variances; ln(auto/Michael) on pole, ln px, ln tan i, ln D by frame.
(C) Error in hollow diameter, a regressor: within-hollow scatter of ln D among repeat detections at <= 20 m/px (out_v2/detections_qc.csv),
    the reliability lambda = 1 - var_err / var(ln D | ln px, ln tan i) in the automatic sample, and method-of-moments errors-in-variables
    corrections of the paper's size-controlled coefficients (depth_paper_rows.csv; 200 image bootstraps; full and half error variance).
    (A SIMEX run with quadratic extrapolation was unstable at this reliability and was replaced.)
(D) Censoring: boxes with no shadow component (n_runs == 0) are upper limits h < h_upper; Tobit refit of the size slope and the poleward
    ratio with them included, against OLS on accepted depths; intervals from 100 image bootstraps.
Writes me.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats, optimize
HERE = os.path.dirname(os.path.abspath(__file__)); RR = os.path.dirname(HERE); G23 = f"{RR}/g23_smallest_hollows"; G24 = f"{RR}/g24_blewett_compare"
out = {}
def ci(m, k): return [round(float(m.params[k]), 3), round(float(m.bse[k]), 3), round(float(m.pvalues[k]), 3)]
# ---- (A) bench
b = pd.read_csv(f"{G23}/bench_depths.csv"); man = json.load(open(f"{G23}/manifest_meas.json"))
b = b[(b.cls == "hollow") & (b.shadow == "measured") & b.auto_h.notna() & (b.auto_h > 0) & (b.depth_m > 0)].copy()
b["lat"] = [man[str(int(h))]["lat"] for h in b.hollow_id]; b["img"] = [man[str(int(h))]["img"] for h in b.hollow_id]
b["y"] = np.log(b.depth_m / b.auto_h); b["lnD"] = np.log(b.D_m / 1000); b["lnpx"] = np.log(b.px); b["ltan"] = np.log(np.tan(np.radians(b.inc)))
b["pole"] = (b.lat >= 60).astype(int); fl = b.auto_flags.fillna(""); b["multi"] = fl.str.contains("multi").astype(int); b["edge"] = fl.str.contains("edge").astype(int)
b["unc"] = b.uncertain.astype(bool).astype(int)
mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
A = dict(n=len(b), images=int(b.img.nunique()), median_ratio=round(float(np.exp(b.y.median())), 3), sd_ln=round(float(b.y.std()), 3), robust_sd_ln=round(float(mad(b.y)), 3), n_pole=int(b.pole.sum()))
m = smf.ols("y ~ lnD + lnpx + ltan + pole + multi + edge + unc", b).fit(cov_type="cluster", cov_kwds=dict(groups=b.img))
A["fit_all"] = {k: ci(m, k) for k in ["lnD", "lnpx", "ltan", "pole", "multi", "edge", "unc"]}
c = b[(b.edge == 0) & (b.unc == 0)]
m2 = smf.ols("y ~ lnD + lnpx + ltan + pole + multi", c).fit(cov_type="cluster", cov_kwds=dict(groups=c.img))
A["fit_clean"] = dict(n=len(c), **{k: ci(m2, k) for k in ["lnD", "lnpx", "ltan", "pole", "multi"]})
# effect on the hand-data size slope if hand depths are replaced by automatic ones on the same hollows
for nm, col in [("hand", "depth_m"), ("auto", "auto_h")]:
    mm = smf.ols(f"np.log({col}) ~ lnD + lnpx + ltan", c).fit(cov_type="cluster", cov_kwds=dict(groups=c.img)); A[f"size_slope_same_hollows_{nm}"] = ci(mm, "lnD")
out["A_bench_hand_vs_auto"] = A
# ---- (B) three measurements of the same boxes
pts = pd.read_csv(f"{G24}/cmp_points.csv"); key = pd.read_csv(f"{G24}/key_cmp.csv"); mc = json.load(open(f"{G24}/manifest_cmp.json"))
ok = pts[(pts.skip != True) & (pts.m_depth > 0)]
k = key.merge(ok.groupby("fid").m_depth.max().rename("M"), on="fid").rename(columns={"b_max": "B", "auto_depth_pilot": "Au"})
k = k[(k.Au > 0) & (k.B > 0) & (k.M > 0)].copy(); k["D"] = [mc[str(f)]["D"] for f in k.fid]; k["px"] = [mc[str(f)]["px"] for f in k.fid]; k["inc"] = [mc[str(f)]["inc"] for f in k.fid]
lB, lM, lA = np.log(k.B), np.log(k.M), np.log(k.Au)
def hat(var):
    vMB, vMA, vBA = var(lM - lB), var(lM - lA), var(lB - lA)
    sM = (vMB + vMA - vBA) / 2; sB = (vMB + vBA - vMA) / 2; sA = (vMA + vBA - vMB) / 2
    return {kk: round(float(np.sqrt(v)) if v > 0 else float("nan"), 3) for kk, v in [("hand", sM), ("blewett", sB), ("auto", sA)]}
rng = np.random.default_rng(3); boots = []
fr = k.num.unique()
for _ in range(1000):
    pick = rng.choice(fr, len(fr)); s = pd.concat([k[k.num == f] for f in pick])
    l1, l2, l3 = np.log(s.B), np.log(s.M), np.log(s.Au); vMB, vMA, vBA = np.var(l2 - l1, ddof=1), np.var(l2 - l3, ddof=1), np.var(l1 - l3, ddof=1)
    boots.append([(vMB + vMA - vBA) / 2, (vMB + vBA - vMA) / 2, (vMA + vBA - vMB) / 2])
boots = np.array(boots); q = lambda j: [round(float(np.sqrt(max(x, 0))), 3) for x in np.percentile(boots[:, j], [2.5, 97.5])]
B = dict(n_boxes=len(k), frames=int(k.num.nunique()), sd_plain=hat(lambda v: np.var(v, ddof=1)), sd_robust=hat(lambda v: mad(v) ** 2),
         ci_plain=dict(hand=q(0), blewett=q(1), auto=q(2)))
k["y"] = lA - lM; k["lnD"] = np.log(k.D / 1000); k["lnpx"] = np.log(k.px); k["ltan"] = np.log(np.tan(np.radians(k.inc))); k["pole"] = k.pole.astype(int)
m3 = smf.ols("y ~ pole + lnD + lnpx + ltan", k).fit(cov_type="cluster", cov_kwds=dict(groups=k.num))
B["auto_over_hand"] = dict(median=round(float(np.exp(k.y.median())), 3), **{kk: ci(m3, kk) for kk in ["pole", "lnD", "lnpx", "ltan"]})
out["B_three_methods_same_boxes"] = B
# ---- (C) error in hollow diameter, and SIMEX
d = pd.read_csv(f"{RR}/../../out_v2/detections_qc.csv"); d = d[d.px_m <= 20]; d["l"] = np.log(d.D_eq_km)
g = d.groupby("hollow_id").l; nn = g.size(); w = g.var(ddof=1)[nn >= 2]; nw = nn[nn >= 2]
var_err = float((w * (nw - 1)).sum() / (nw - 1).sum())
a = pd.read_csv(f"{G23}/depth_paper_rows.csv"); a["img"] = a.img.astype(str)
r = smf.ols("lnD ~ lnpx + ltan", a).fit(); var_cond = float(r.resid.var())
C = dict(hollows_with_repeats=int(len(w)), detections=int(nw.sum()), sd_err_lnD=round(np.sqrt(var_err), 3), sd_lnD_sample=round(float(a.lnD.std()), 3),
         sd_lnD_given_px_inc=round(np.sqrt(var_cond), 3), reliability=round(1 - var_err / var_cond, 3))
# errors-in-variables by the method of moments (classical error in ln D only): beta = (X'X - n S_uu)^-1 X'y; intervals from 300
# bootstrap resamples of images; also with half the error variance, in case repeat detections overstate the error of one diameter.
models = {"size": ["lnD", "lnpx", "ltan"], "pole": ["pole", "lnD", "lnpx", "ltan"], "T": ["x", "lnD", "lnpx", "ltan"],
          "nested": ["nested", "lnD", "lnpx", "ltan"], "kinczyk": ["k3", "k45", "lnD", "lnpx", "ltan"], "fresh": ["fresh", "degraded", "lnD", "lnpx", "ltan"]}
def mom_xy(X, y, j, ve):
    S = np.zeros((X.shape[1], X.shape[1])); S[j, j] = ve * len(y)
    return np.linalg.solve(X.T @ X - S, X.T @ y)[1:], np.linalg.lstsq(X, y, rcond=None)[0][1:]
rng = np.random.default_rng(11); eiv = {}
for nm, cols in models.items():
    s = a.dropna(subset=cols + ["lnh"]); X = np.column_stack([np.ones(len(s))] + [s[c_].values.astype(float) for c_ in cols]); y = s.lnh.values
    j = 1 + cols.index("lnD"); ims = s.img.values; uim = np.unique(ims); gi = {g_: np.where(ims == g_)[0] for g_ in uim}
    picks = [np.concatenate([gi[g_] for g_ in rng.choice(uim, len(uim))]) for _ in range(200)]
    res = {}
    for tag, ve in [("full", var_err), ("half", var_err / 2)]:
        b_eiv, b_ols = mom_xy(X, y, j, ve); bs = np.array([mom_xy(X[ii], y[ii], j, ve)[0] for ii in picks])
        lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        res[tag] = {c_: dict(ols=round(float(b_ols[k_]), 4), eiv=round(float(b_eiv[k_]), 4), ci=[round(float(lo[k_]), 4), round(float(hi[k_]), 4)]) for k_, c_ in enumerate(cols) if c_ not in ("lnpx", "ltan")}
    eiv[nm] = res
C["eiv"] = eiv
ex = lambda v: round(float(np.exp(v)), 3)
C["paper_quantities"] = {tag: dict(
    size_exponent=[eiv["size"][tag]["lnD"]["ols"], eiv["size"][tag]["lnD"]["eiv"], eiv["size"][tag]["lnD"]["ci"]],
    poleward_ratio=[ex(eiv["pole"][tag]["pole"]["ols"]), ex(eiv["pole"][tag]["pole"]["eiv"]), [ex(v) for v in eiv["pole"][tag]["pole"]["ci"]]],
    kinczyk3=[ex(eiv["kinczyk"][tag]["k3"]["ols"]), ex(eiv["kinczyk"][tag]["k3"]["eiv"]), [ex(v) for v in eiv["kinczyk"][tag]["k3"]["ci"]]],
    kinczyk45=[ex(eiv["kinczyk"][tag]["k45"]["ols"]), ex(eiv["kinczyk"][tag]["k45"]["eiv"]), [ex(v) for v in eiv["kinczyk"][tag]["k45"]["ci"]]],
    nested=[ex(eiv["nested"][tag]["nested"]["ols"]), ex(eiv["nested"][tag]["nested"]["eiv"]), [ex(v) for v in eiv["nested"][tag]["nested"]["ci"]]],
    fresh_vs_degraded=[ex(eiv["fresh"][tag]["fresh"]["ols"] - eiv["fresh"][tag]["degraded"]["ols"]), ex(eiv["fresh"][tag]["fresh"]["eiv"] - eiv["fresh"][tag]["degraded"]["eiv"])],
    T_ratio_324_219=[ex(eiv["T"][tag]["x"]["ols"] * (1 / 324.3 - 1 / 219.0)), ex(eiv["T"][tag]["x"]["eiv"] * (1 / 324.3 - 1 / 219.0)), [ex(v * (1 / 324.3 - 1 / 219.0)) for v in eiv["T"][tag]["x"]["ci"]][::-1]])
    for tag in ["full", "half"]}
out["C_diameter_error"] = C
# ---- (D) censoring: boxes without a shadow component
al = pd.read_csv(f"{G23}/depth_all.csv"); flg = al["flags"].fillna("")
acc = al[(al.depth_m > 0) & ~flg.str.contains("no_shadow|edge_sun|frame_edge|no_edr")].copy(); acc["k"] = "clean"
up = al[flg.str.contains("no_shadow") & (al.n_runs == 0) & (al.h_upper_m > 0) & ~flg.str.contains("edge_sun|frame_edge|no_edr")].copy(); up["k"] = "upper"
u = pd.concat([acc, up]); u["lnD"] = np.log(u.D_can_m / 1000); u["lnpx"] = np.log(u.px); u["ltan"] = np.log(np.tan(np.radians(u.inc))); u["pole"] = (u.lat >= 60).astype(int)
u["yy"] = np.where(u.k == "clean", np.log(u.depth_m.clip(lower=1e-3)), np.log(u.h_upper_m.clip(lower=1e-3)))
X = np.column_stack([np.ones(len(u)), u.lnD, u.lnpx, u.ltan, u.pole]); yy = u.yy.values; kk = u.k.values
def nll(th):
    mu = X @ th[:-1]; s_ = np.exp(th[-1]); z = (yy - mu) / s_
    return -(np.sum(stats.norm.logpdf(z[kk == "clean"]) - th[-1]) + np.sum(stats.norm.logcdf(z[kk == "upper"])))
o = smf.ols("np.log(depth_m) ~ lnD + lnpx + ltan + pole", u[u.k == "clean"]).fit()
f_ = optimize.minimize(nll, np.r_[o.params.values, np.log(np.sqrt(o.scale))], method="BFGS")
Hi = f_.hess_inv; se = np.sqrt(np.diag(Hi))
u["img"] = u.num.astype(str); ims = u.img.unique(); gi = {g_: np.where(u.img.values == g_)[0] for g_ in ims}; tb = []
for _ in range(100):
    ii = np.concatenate([gi[g_] for g_ in rng.choice(ims, len(ims))]); Xb, yb, kb = X[ii], yy[ii], kk[ii]
    def nllb(th):
        mu = Xb @ th[:-1]; s_ = np.exp(th[-1]); z = (yb - mu) / s_
        return -(np.sum(stats.norm.logpdf(z[kb == "clean"]) - th[-1]) + np.sum(stats.norm.logcdf(z[kb == "upper"])))
    tb.append(optimize.minimize(nllb, f_.x, method="BFGS").x)
tb = np.array(tb)
out["D_censoring"] = dict(tobit_boot=dict(lnD=[round(float(v), 3) for v in np.percentile(tb[:, 1], [2.5, 97.5])], pole_ratio=[round(float(np.exp(v)), 3) for v in np.percentile(tb[:, 4], [2.5, 97.5])]), n_clean=int((kk == "clean").sum()), n_upper=int((kk == "upper").sum()), frac_upper=round(float((kk == "upper").mean()), 3),
                          ols=dict(lnD=round(float(o.params.lnD), 3), pole_ratio=round(float(np.exp(o.params.pole)), 3)),
                          tobit=dict(lnD=[round(float(f_.x[1]), 3), round(float(se[1]), 3)], pole_ratio=round(float(np.exp(f_.x[4])), 3), pole_se=round(float(se[4]), 3)))
json.dump(out, open(f"{HERE}/me.json", "w"), indent=1); print(json.dumps(out, indent=1))
