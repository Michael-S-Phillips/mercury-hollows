#!/usr/bin/env python3
"""g26 (9 Nov 2026), Michael: the boxes without a shadow plot as 'upper limits' at the bottom of the depth cloud; are they shallow hollows?
For every box (accepted automatic depths and boxes with no shadow component, n_runs == 0), the probability that a hollow following the
fitted relation ln h ~ ln D + ln px + ln tan i (least squares on accepted depths) falls below its image's threshold h_upper, with the
residual scatter of the fit and with 30% more scatter (the fit to measured depths is censored, so its scatter is low). Sum by size bin
against the observed number of no-shadow boxes. Writes no_shadow_check.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); G23 = os.path.join(os.path.dirname(HERE), "g23_smallest_hollows")
a = pd.read_csv(f"{G23}/depth_all.csv"); fl = a["flags"].fillna(""); base = ~fl.str.contains("edge_sun|frame_edge|no_edr")
acc = a[base & (a.depth_m > 0) & ~fl.str.contains("no_shadow")].copy(); up = a[base & fl.str.contains("no_shadow") & (a.n_runs == 0) & (a.h_upper_m > 0)].copy()
for d in (acc, up): d["lnD"] = np.log(d.D_can_m / 1000); d["lnpx"] = np.log(d.px); d["ltan"] = np.log(np.tan(np.radians(d.inc)))
m = smf.ols("np.log(depth_m) ~ lnD + lnpx + ltan", acc).fit(); sd = float(np.sqrt(m.scale))
b = pd.concat([acc.assign(k="acc"), up.assign(k="up")]); mu = m.predict(b); thr = np.log(b.h_upper_m)
b["p"] = stats.norm.cdf((thr - mu) / sd); b["p2"] = stats.norm.cdf((thr - mu) / (1.3 * sd)); b["bin"] = pd.cut(b.D_can_m, [0, 150, 300, 500, 800, 1500, 5000])
rows = []
for k, x in b.groupby("bin", observed=True):
    rows.append(dict(bin=str(k), n=len(x), observed_noshadow=int((x.k == "up").sum()), expected=round(float(x.p.sum()), 1), expected_wider=round(float(x.p2.sum()), 1),
                     median_threshold_m=round(float(x.h_upper_m.median()), 1), median_predicted_m=round(float(np.exp(m.predict(x).median())), 1)))
out = dict(resid_sd=round(sd, 3), bins=rows, total=dict(observed=int((b.k == "up").sum()), expected=round(float(b.p.sum()), 1), expected_wider=round(float(b.p2.sum()), 1)))
json.dump(out, open(f"{HERE}/no_shadow_check.json", "w"), indent=1); print(json.dumps(out["total"]), out["resid_sd"])
