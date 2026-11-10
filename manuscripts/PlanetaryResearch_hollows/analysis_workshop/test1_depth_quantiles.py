#!/usr/bin/env python3
"""Workshop test 1 (2 Nov 2026): does the poleward decrease in hollow depth shift the whole depth distribution
(imaging bias, or a depth that scales with temperature) or only its lower tail (hollows still growing toward a
stopping depth, version E)?

Input: GRL_insolation_hollows/data/test1_rows_inc.csv (2,518 shadow-length depths of Blewett et al. 2016, with px, inc,
T, image), as in depth_reanalysis.py and Section S9. Poleward = |lat| >= 60.
Analyses
 A. OLS replication: ln h ~ pole + ln px + ln tan(inc), SE clustered by image (Section S9 ratio 0.77).
 B. Quantile regression of ln h on pole + ln px + ln tan(inc) at tau = 0.1, 0.25, 0.5, 0.75, 0.9; 95% intervals from
    400 bootstrap resamples of images.
 C. Detection-threshold check: the shortest measurable shadow sets a minimum depth h_min = k px / tan(inc) that differs
    between groups. The shortest shadows in the data are ~2 pixels (5th percentile 2.2 px), so k = 2 is the empirical
    floor. Keep only hollows deeper than a common threshold D0 at which >= 95% of measurements in BOTH groups could have
    been measured with a shadow of k pixels (k = 2, 3); compare quantiles of h | h >= D0 (with and without the px/inc
    correction), bootstrap by image.
 D. Equatorial-imaging resampling: each poleward depth is kept only if it would cast a shadow of >= 2 pixels in an
    equatorward image (px, inc) drawn at random; 2,000 draws; quantiles of the kept poleward depths against the
    equatorward depths.
Writes test1_depth_quantiles_record.json.
"""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
ROOT = os.path.expanduser("~/mnt/hollows-emissivity")
F = f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv"
d = pd.read_csv(F)
d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["img"] = d.filename.astype(str)
d["pole"] = (d.abslat >= 60).astype(int)
d["shadow_px"] = d.depth * np.tan(np.radians(d.inc)) / d.px
R = {"input": F, "n": int(len(d)), "n_pole": int(d.pole.sum()), "n_img": int(d.img.nunique())}
R["groups"] = {g: dict(n=int(len(s)), n_img=int(s.img.nunique()), px_med=float(s.px.median()), inc_med=float(s.inc.median()),
                      depth_q=[float(x) for x in s.depth.quantile([.1, .25, .5, .75, .9])],
                      shadow_px_med=float(s.shadow_px.median()))
               for g, s in [("equatorward", d[d.pole == 0]), ("poleward", d[d.pole == 1])]}
# A
m = smf.ols("lnh ~ pole + lnpx + ltan", d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.img)[0]})
b, se = m.params["pole"], m.bse["pole"]
R["A_ols"] = dict(ratio=float(np.exp(b)), lo=float(np.exp(b - 1.96 * se)), hi=float(np.exp(b + 1.96 * se)))
# B
taus = [0.1, 0.25, 0.5, 0.75, 0.9]
def qfit(df):
    out = []
    for t in taus:
        r = smf.quantreg("lnh ~ pole + lnpx + ltan", df).fit(q=t, max_iter=5000)
        out.append(r.params["pole"])
    return np.array(out)
est = qfit(d)
rng = np.random.default_rng(42); imgs = d.img.unique(); groups = {k: v.index.values for k, v in d.groupby("img")}
B = 400; boots = []
for _ in range(B):
    pick = rng.choice(imgs, len(imgs), replace=True)
    idx = np.concatenate([groups[i] for i in pick]); s = d.loc[idx].reset_index(drop=True)
    if s.pole.sum() < 20: continue
    try: boots.append(qfit(s))
    except Exception: pass
boots = np.array(boots)
R["B_quantreg"] = {f"tau_{t}": dict(ratio=float(np.exp(e)), lo=float(np.exp(np.percentile(boots[:, i], 2.5))),
                                    hi=float(np.exp(np.percentile(boots[:, i], 97.5))))
                   for i, (t, e) in enumerate(zip(taus, est))}
R["B_quantreg"]["n_boot"] = int(len(boots))
# slope difference across quantiles: tau 0.9 minus tau 0.1
diff = boots[:, -1] - boots[:, 0]
R["B_quantreg"]["ratio_tau0.9_over_tau0.1"] = dict(est=float(np.exp(est[-1] - est[0])),
    lo=float(np.exp(np.percentile(diff, 2.5))), hi=float(np.exp(np.percentile(diff, 97.5))))
# C
qs = [0.1, 0.25, 0.5, 0.75, 0.9]
mc = smf.ols("lnh ~ lnpx + ltan", d).fit()
d["h_corr"] = np.exp(d.lnh - mc.params["lnpx"] * (d.lnpx - d.lnpx.median()) - mc.params["ltan"] * (d.ltan - d.ltan.median()))
def qratio(s, col):
    a, p = s[s.pole == 0][col], s[s.pole == 1][col]
    return np.array([p.quantile(q) / a.quantile(q) for q in qs])
R["C_threshold"] = {}
for k in [2, 3]:
    hk = k * d.px / np.tan(np.radians(d.inc))
    D0 = float(max(hk[d.pole == 0].quantile(0.95), hk[d.pole == 1].quantile(0.95)))
    t = d[d.depth >= D0].copy(); tg = {k_: v.index.values for k_, v in t.groupby("img")}; timgs = t.img.unique()
    rec = dict(D0=D0, n_equatorward=int((t.pole == 0).sum()), n_poleward=int(t.pole.sum()))
    for col in ["depth", "h_corr"]:
        est_c = qratio(t, col); bb = []
        for _ in range(B):
            pick = rng.choice(timgs, len(timgs), replace=True)
            s = t.loc[np.concatenate([tg[i] for i in pick])]
            if s.pole.sum() < 8 or (s.pole == 0).sum() < 8: continue
            bb.append(qratio(s, col))
        bb = np.array(bb)
        rec[col] = {str(q): dict(ratio=float(e), lo=float(np.percentile(bb[:, i], 2.5)), hi=float(np.percentile(bb[:, i], 97.5)))
                    for i, (q, e) in enumerate(zip(qs, est_c))}
    R["C_threshold"][f"k{k}"] = rec
# D
eq = d[d.pole == 0]; po = d[d.pole == 1]
thr_eq = (2 * eq.px / np.tan(np.radians(eq.inc))).values
kept = []
for _ in range(2000):
    th = rng.choice(thr_eq, len(po), replace=True)
    kept.append(po.depth.values[po.depth.values >= th])
allk = np.concatenate(kept)
R["D_resample"] = dict(frac_kept=float(len(allk) / (2000 * len(po))),
    quantile_ratio_pole_over_eq={str(q): float(np.quantile(allk, q) / eq.depth.quantile(q)) for q in qs},
    pole_raw_over_eq={str(q): float(po.depth.quantile(q) / eq.depth.quantile(q)) for q in qs})
json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test1_depth_quantiles_record.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
