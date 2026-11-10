#!/usr/bin/env python3
"""g23 step 4: depth against size from depth_vs_size.csv. Classes: clean measurement (depth > 0, no pilot BAD flag: no_shadow, edge_sun,
frame_edge); upper limit (no shadow: h < 1.5 px); lower limit (box >= 85% dark: h >= full_frac x h_fill). Bins of canonical diameter
within each pixel-scale tier; OLS of ln h on ln D, ln px and ln tan i for clean measurements; and a censored (Tobit) fit in log space
that also uses the limits. Writes dvs_summary.json."""
import json, numpy as np, pandas as pd
from scipy import optimize, stats
r = pd.read_csv("depth_vs_size.csv"); fl = r["flags"].fillna("")
r = r[~fl.str.contains("no_edr")].copy(); fl = r["flags"].fillna("")
bad = fl.str.contains("edge_sun|frame_edge")
r["cls"] = "excluded"
r.loc[(r.full_frac >= 0.85) & ~bad, "cls"] = "lower"
# "no_shadow" also marks boxes whose dark runs were all rejected (not bounded by a lit far wall); those have shadow components
# (n_runs >= 1) and are excluded rather than counted as upper limits
r.loc[fl.str.contains("no_shadow") & ~bad & (r.cls != "lower") & (r.n_runs == 0), "cls"] = "upper"
r.loc[(r.depth_m > 0) & ~bad & ~fl.str.contains("no_shadow") & (r.cls != "lower"), "cls"] = "clean"
r["h_lo"] = (r.full_frac * r.h_fill_m).where(r.cls == "lower")
bins = [30, 100, 150, 200, 300, 500, 1000, 3000]; r["Dbin"] = pd.cut(r.D_can_m, bins)
out = {}
for t in ["le7", "7to15"]:
    q = r[r.tier == t]; rows = []
    for b, g in q.groupby("Dbin", observed=True):
        c = g[g.cls == "clean"]
        rows.append(dict(D=str(b), n=len(g), clean=len(c), upper=int((g.cls == "upper").sum()), lower=int((g.cls == "lower").sum()),
                         excl=int((g.cls == "excluded").sum()), px_med=round(g.px.median(), 1), inc_med=round(g.inc.median(), 1),
                         h_med=round(c.depth_m.median(), 1) if len(c) else None,
                         h_q25=round(c.depth_m.quantile(.25), 1) if len(c) else None, h_q75=round(c.depth_m.quantile(.75), 1) if len(c) else None,
                         h_upper_med=round(g[g.cls == "upper"].h_upper_m.median(), 1) if (g.cls == "upper").any() else None,
                         h_lo_med=round(g.h_lo.median(), 1) if (g.cls == "lower").any() else None))
    out[t] = rows; print(t); print(pd.DataFrame(rows).to_string(index=False))
# OLS on clean
c = r[r.cls == "clean"].copy(); c["lh"] = np.log(c.depth_m); c["lD"] = np.log(c.D_can_m); c["lpx"] = np.log(c.px); c["lti"] = np.log(np.tan(np.radians(c.inc)))
X = np.c_[np.ones(len(c)), c.lD, c.lpx, c.lti]; beta, *_ = np.linalg.lstsq(X, c.lh, rcond=None)
res = c.lh - X @ beta; s2 = res.var(ddof=4); cov = s2 * np.linalg.inv(X.T @ X); se = np.sqrt(np.diag(cov))
out["ols_clean"] = dict(n=len(c), coef=dict(zip(["const", "lnD", "lnpx", "lntan_i"], beta.round(3).tolist())), se=dict(zip(["const", "lnD", "lnpx", "lntan_i"], se.round(3).tolist())))
print("OLS clean", out["ols_clean"])
for t in ["le7", "7to15"]:
    q = c[c.tier == t]; rho, p = stats.spearmanr(q.D_can_m, q.depth_m); out[f"spearman_{t}"] = dict(n=len(q), rho=round(rho, 3), p=float(p)); print(t, "spearman D-h", round(rho, 3), p, len(q))
# Tobit in log space: clean = exact, upper = left-censored at h_upper, lower = right-censored at h_lo
u = r[r.cls.isin(["clean", "upper", "lower"])].copy()
y = np.where(u.cls == "clean", np.log(u.depth_m), np.where(u.cls == "upper", np.log(u.h_upper_m), np.log(u.h_lo.clip(lower=0.1))))
Xu = np.c_[np.ones(len(u)), np.log(u.D_can_m), np.log(u.px), np.log(np.tan(np.radians(u.inc)))]; k = u.cls.values
def nll(th):
    b, ls = th[:-1], th[-1]; s = np.exp(ls); mu = Xu @ b; z = (y - mu) / s
    return -(np.sum(stats.norm.logpdf(z[k == "clean"]) - ls) + np.sum(stats.norm.logcdf(z[k == "upper"])) + np.sum(stats.norm.logsf(z[k == "lower"])))
th0 = np.r_[beta, np.log(np.sqrt(s2))]; fit = optimize.minimize(nll, th0, method="BFGS")
H = fit.hess_inv; se_t = np.sqrt(np.diag(H))
out["tobit"] = dict(n=len(u), n_clean=int((k == "clean").sum()), n_upper=int((k == "upper").sum()), n_lower=int((k == "lower").sum()),
                    coef=dict(zip(["const", "lnD", "lnpx", "lntan_i", "ln_sigma"], fit.x.round(3).tolist())), se=dict(zip(["const", "lnD", "lnpx", "lntan_i", "ln_sigma"], se_t.round(3).tolist())))
print("Tobit", out["tobit"])
# Tobit with a break: slope of ln D below and above D_b (profile over D_b)
def tobit_break(Db):
    X2 = np.c_[np.ones(len(u)), np.minimum(np.log(u.D_can_m), np.log(Db)), np.maximum(np.log(u.D_can_m) - np.log(Db), 0), np.log(u.px), np.log(np.tan(np.radians(u.inc)))]
    def nll2(th):
        b, ls = th[:-1], th[-1]; s = np.exp(ls); mu = X2 @ b; z = (y - mu) / s
        return -(np.sum(stats.norm.logpdf(z[k == "clean"]) - ls) + np.sum(stats.norm.logcdf(z[k == "upper"])) + np.sum(stats.norm.logsf(z[k == "lower"])))
    f = optimize.minimize(nll2, np.r_[fit.x[0], fit.x[1], fit.x[1], fit.x[2], fit.x[3], fit.x[4]], method="BFGS")
    return f.fun, f.x, np.sqrt(np.diag(f.hess_inv))
prof = []
for Db in [80, 100, 150, 200, 300, 500, 700]:
    f, x, se2 = tobit_break(Db); prof.append(dict(Db=Db, nll=round(f, 2), slope_below=round(x[1], 3), se_below=round(se2[1], 3), slope_above=round(x[2], 3), se_above=round(se2[2], 3)))
out["tobit_break_profile"] = prof; out["tobit_nll_nobreak"] = round(fit.fun, 2); print(pd.DataFrame(prof).to_string(index=False), "no-break nll", round(fit.fun, 2))
json.dump(out, open("dvs_summary.json", "w"), indent=1, default=str)
r.to_csv("depth_vs_size_classed.csv", index=False)
