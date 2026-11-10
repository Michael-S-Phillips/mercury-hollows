#!/usr/bin/env python3
"""g23 step 5: checks on the depth-size relation. (1) Subsets: single shadow runs only (n_runs <= 1, no 'multi'); emission <= 30;
confidence >= 0.8. (2) Tobit fit of the arching curve: h = h_inf * g(D_t / D) with g(x) = 2 (x - 1 + e^-x) / x^2, the mean load on a
column D wide relative to an unconfined column (Janssen; load grows as D for D << D_t and is unchanged for D >> D_t), plus ln px and
ln tan i terms, against the power-law Tobit of dvs_analyze.py. Writes dvs_checks.json."""
import json, numpy as np, pandas as pd
from scipy import optimize, stats
r = pd.read_csv("depth_vs_size_classed.csv"); u = r[r.cls.isin(["clean", "upper", "lower"])].copy()
def ycls(u):
    k = u.cls.values; y = np.where(k == "clean", np.log(u.depth_m.where(k == "clean", 1.0)), np.where(k == "upper", np.log(u.h_upper_m), np.log(u.h_lo.fillna(1).clip(lower=0.1))))
    return y, k
def tobit(mu_fn, th0, u):
    y, k = ycls(u)
    def nll(th):
        mu = mu_fn(th[:-1], u); s = np.exp(th[-1]); z = (y - mu) / s
        return -(np.sum(stats.norm.logpdf(z[k == "clean"]) - th[-1]) + np.sum(stats.norm.logcdf(z[k == "upper"])) + np.sum(stats.norm.logsf(z[k == "lower"])))
    f = optimize.minimize(nll, th0, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-6, fatol=1e-6))
    f = optimize.minimize(nll, f.x, method="BFGS"); return f
lin = lambda b, u: b[0] + b[1] * np.log(u.D_can_m) + b[2] * np.log(u.px) + b[3] * np.log(np.tan(np.radians(u.inc)))
def g(x): x = np.asarray(x, float); return np.where(x < 1e-4, 1 - x / 3, 2 * (x - 1 + np.exp(-x)) / np.maximum(x, 1e-12) ** 2)
arch = lambda b, u: b[0] + np.log(g(np.exp(b[1]) / u.D_can_m)) + b[2] * np.log(u.px) + b[3] * np.log(np.tan(np.radians(u.inc)))
out = {}
for name, sub in [("all", u), ("single_run", u[(u.n_runs <= 1) & ~u["flags"].fillna("").str.contains("multi")]), ("emi_le30", u[u.emi <= 30]), ("conf_ge0.8", u[u.conf >= 0.8])]:
    f = tobit(lin, np.r_[-1, 0.6, 0.1, -0.3, -0.1], sub); se = np.sqrt(np.diag(f.hess_inv))
    c = sub[sub.cls == "clean"]; rho = stats.spearmanr(c.D_can_m, c.depth_m)[0]
    out[name] = dict(n=len(sub), n_clean=len(c), slope_lnD=round(f.x[1], 3), se=round(se[1], 3), spearman_clean=round(rho, 3), nll=round(f.fun, 2))
    print(name, out[name])
fa = tobit(arch, np.r_[np.log(40), np.log(400), 0.1, -0.3, -0.1], u); sea = np.sqrt(np.diag(fa.hess_inv))
fl = tobit(lin, np.r_[-1, 0.6, 0.1, -0.3, -0.1], u)
# profile likelihood interval for D_t
grid = np.exp(np.linspace(np.log(100), np.log(5000), 40)); prof = []
y, k = ycls(u)
for Dt in grid:
    def nll_fix(th):
        mu = th[0] + np.log(g(Dt / u.D_can_m)) + th[1] * np.log(u.px) + th[2] * np.log(np.tan(np.radians(u.inc))); s = np.exp(th[3]); z = (y - mu) / s
        return -(np.sum(stats.norm.logpdf(z[k == "clean"]) - th[3]) + np.sum(stats.norm.logcdf(z[k == "upper"])) + np.sum(stats.norm.logsf(z[k == "lower"])))
    ff = optimize.minimize(nll_fix, np.r_[fa.x[0], fa.x[2], fa.x[3], fa.x[4]], method="BFGS"); prof.append(ff.fun)
prof = np.array(prof); ok = prof <= prof.min() + 1.92; Dt_ci = (float(grid[ok].min()), float(grid[ok].max()))
# implied h_inf at the median px and incidence of the 7-15 m tier, and the width scale
pxm, im = u[u.tier == "7to15"].px.median(), u[u.tier == "7to15"].inc.median()
h_inf = float(np.exp(fa.x[0] + fa.x[2] * np.log(pxm) + fa.x[3] * np.log(np.tan(np.radians(im)))))
out["arching_fit"] = dict(h_inf_m_at_tier2_median=round(h_inf, 1), Dt_m=round(float(np.exp(fa.x[1])), 0), Dt_ci95_profile=[round(Dt_ci[0]), round(Dt_ci[1])],
                          lnpx=round(fa.x[2], 3), lntan_i=round(fa.x[3], 3), sigma=round(float(np.exp(fa.x[4])), 3), nll=round(fa.fun, 2), nll_powerlaw=round(fl.fun, 2),
                          aic_arch=round(2 * fa.fun + 10, 1), aic_powerlaw=round(2 * fl.fun + 10, 1), px_med=round(pxm, 1), inc_med=round(im, 1))
print(out["arching_fit"])
# predicted h/h_inf at sizes
Dt = np.exp(fa.x[1]); out["arching_curve"] = {str(D): round(float(g(Dt / D)), 3) for D in [50, 100, 200, 300, 500, 1000, 2000, 5000]}; print(out["arching_curve"])
json.dump(out, open("dvs_checks.json", "w"), indent=1)
