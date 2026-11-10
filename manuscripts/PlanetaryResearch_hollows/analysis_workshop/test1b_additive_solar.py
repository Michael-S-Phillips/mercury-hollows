#!/usr/bin/env python3
"""Workshop test 1b (2 Nov 2026): is hollow depth an impact-set floor plus a sunlight-driven addition?
If so, the shallowest hollows (low quantiles) should not change with temperature, and the deepest (high quantiles)
should deepen with mean annual temperature T and toward the hot longitudes (cos 2*lon), at fixed pixel scale and incidence.
Quantile regression of ln h on T (per 10 K) [or cos 2*lon] + ln px + ln tan(inc); 300 bootstrap resamples of images.
Samples: all 2,518 depths; equatorward of 60 deg only. Writes test1b_additive_solar_record.json.
"""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
ROOT = os.path.expanduser("~/mnt/hollows-emissivity")
d = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")
d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["img"] = d.filename.astype(str)
d["T10"] = d["T"] / 10.0; d["c2"] = np.cos(np.radians(2 * d.longitude))
taus = [0.1, 0.25, 0.5, 0.75, 0.9]; rng = np.random.default_rng(7); R = {}
def fit(df, var):
    return np.array([smf.quantreg(f"lnh ~ {var} + lnpx + ltan", df).fit(q=t, max_iter=5000).params[var] for t in taus])
for name, df in [("all", d), ("equatorward", d[d.abslat < 60].reset_index(drop=True))]:
    gi = {k: v.index.values for k, v in df.groupby("img")}; imgs = df.img.unique()
    for var in ["T10", "c2"]:
        est = fit(df, var); bb = []
        for _ in range(300):
            pick = rng.choice(imgs, len(imgs), replace=True)
            try: bb.append(fit(df.loc[np.concatenate([gi[i] for i in pick])].reset_index(drop=True), var))
            except Exception: pass
        bb = np.array(bb)
        R[f"{name}_{var}"] = {str(t): dict(factor=float(np.exp(e)), lo=float(np.exp(np.percentile(bb[:, i], 2.5))),
                                          hi=float(np.exp(np.percentile(bb[:, i], 97.5)))) for i, (t, e) in enumerate(zip(taus, est))}
        dd = bb[:, -1] - bb[:, 0]
        R[f"{name}_{var}"]["tau0.9_minus_tau0.1"] = dict(factor=float(np.exp(est[-1] - est[0])),
            lo=float(np.exp(np.percentile(dd, 2.5))), hi=float(np.exp(np.percentile(dd, 97.5))))
json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test1b_additive_solar_record.json"), "w"), indent=1)
for k, v in R.items():
    print(k, {t: (round(x["factor"], 2), round(x["lo"], 2), round(x["hi"], 2)) for t, x in v.items()})
