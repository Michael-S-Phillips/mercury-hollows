#!/usr/bin/env python3
"""g23 step 7: Michael's hand measurements from the Hollow Shadow Bench (bench_shadows_2026-11-09.json, 135 hollows at <= 7 m/pixel).
Depth = L_px x px x f / tan(i) along the sun. Classes: measured (exact), full (lower limit), none (upper limit, 1.5 px).
Uncertain = a note that calls the measurement unsure, unreliable, uncertain, diffuse, blurry, tough, or poorly resolved.
Fits: median by size bin; OLS of ln h on ln D, ln px, ln tan i; Tobit with the limits; a broken power law (profiled break);
a saturating curve h = h_max D^b / (D^b + D_s^b) ... kept simple: power law below D_s, flat above (profiled D_s).
Comparison with the automatic depths (manifest_meas.json auto_h). Writes bench_summary.json and bench_depths.csv."""
import json, re, numpy as np, pandas as pd
from scipy import optimize, stats
b = pd.DataFrame(json.load(open("bench_shadows_2026-11-09.json"))); man = json.load(open("manifest_meas.json"))
b["auto_h"] = [man[str(int(i))]["auto_h"] for i in b.hollow_id]; b["auto_flags"] = [man[str(int(i))]["auto_flags"] for i in b.hollow_id]
b["note"] = b["note"].fillna("")
b["uncertain"] = b.note.str.contains(r"unsure|unreliable|uncertain|diffuse|blurry|tough|poor|bad res|can't|not quite|might fill", case=False, regex=True) | (b.cls == "unsure")
b["h_up"] = 1.5 * b.px * b.f / np.tan(np.radians(b.inc))
b.to_csv("bench_depths.csv", index=False)
out = {"n": len(b), "cls": b.cls.fillna("none").value_counts().to_dict(), "shadow": b.shadow.value_counts().to_dict(), "uncertain": int(b.uncertain.sum())}
hol = b[b.cls.isin(["hollow", "unsure"]) | b.cls.isna()]
m = hol[hol.shadow == "measured"]
bins = [0, 100, 150, 200, 300, 500, 1e4]; lab = ["40-100", "100-150", "150-200", "200-300", "300-500", "500+"]
m = m.assign(Dbin=pd.cut(m.D_m, bins, labels=lab))
t = m.groupby("Dbin", observed=True).agg(n=("depth_m", "size"), h_med=("depth_m", "median"), q25=("depth_m", lambda x: x.quantile(.25)), q75=("depth_m", lambda x: x.quantile(.75)),
    hD=("depth_m", lambda x: np.nan), px=("px", "median"), inc=("inc", "median")).round(2)
t["hD_med"] = m.groupby("Dbin", observed=True).apply(lambda g: (g.depth_m / g.D_m).median()).round(3); t = t.drop(columns="hD")
c = m[~m.uncertain]; t["n_certain"] = c.groupby(pd.cut(c.D_m, bins, labels=lab), observed=True).size(); t["h_med_certain"] = c.groupby(pd.cut(c.D_m, bins, labels=lab), observed=True).depth_m.median().round(2)
print(t.to_string()); out["bins"] = t.reset_index().astype(str).to_dict("records")
def ols(d):
    X = np.c_[np.ones(len(d)), np.log(d.D_m), np.log(d.px), np.log(np.tan(np.radians(d.inc)))]; y = np.log(d.depth_m)
    bb, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ bb; cov = r.var(ddof=4) * np.linalg.inv(X.T @ X); return bb, np.sqrt(np.diag(cov)), r.std(ddof=4)
for nm, d in [("all_measured_hollows", m), ("certain_only", c)]:
    bb, se, sd = ols(d); rho, p = stats.spearmanr(d.D_m, d.depth_m)
    out[f"ols_{nm}"] = dict(n=len(d), lnD=[round(bb[1], 3), round(se[1], 3)], lnpx=[round(bb[2], 3), round(se[2], 3)], lntan_i=[round(bb[3], 3), round(se[3], 3)], resid_sd=round(sd, 3), spearman=[round(rho, 3), float(p)])
    print(nm, out[f"ols_{nm}"])
# D-only fit (no px/inc terms) for reference
bb = np.polyfit(np.log(m.D_m), np.log(m.depth_m), 1); out["powerlaw_D_only"] = dict(slope=round(bb[0], 3), h_at_100=round(float(np.exp(np.polyval(bb, np.log(100)))), 2), h_at_500=round(float(np.exp(np.polyval(bb, np.log(500)))), 2))
print("D only", out["powerlaw_D_only"])
# Tobit with limits (all classed hollow/unsure/none; crater excluded)
u = hol.copy(); k = u.shadow.values
y = np.where(k == "measured", np.log(u.depth_m.where(k == "measured", 1.0)), np.where(k == "none", np.log(u.h_up), np.log(u.depth_m.where(k == "full", 1.0).clip(lower=0.1))))
X = np.c_[np.ones(len(u)), np.log(u.D_m), np.log(u.px), np.log(np.tan(np.radians(u.inc)))]
def nll(th, X=X):
    mu = X @ th[:-1]; s = np.exp(th[-1]); z = (y - mu) / s
    return -(np.sum(stats.norm.logpdf(z[k == "measured"]) - th[-1]) + np.sum(stats.norm.logcdf(z[k == "none"])) + np.sum(stats.norm.logsf(z[k == "full"])))
b0, *_ = np.linalg.lstsq(X[k == "measured"], y[k == "measured"], rcond=None)
f = optimize.minimize(nll, np.r_[b0, -0.5], method="Nelder-Mead", options=dict(maxiter=40000, xatol=1e-9, fatol=1e-9))
def hess(fn, th, h=1e-4):
    n = len(th); H = np.zeros((n, n)); E = np.eye(n) * h
    for i in range(n):
        for j in range(n): H[i, j] = (fn(th + E[i] + E[j]) - fn(th + E[i] - E[j]) - fn(th - E[i] + E[j]) + fn(th - E[i] - E[j])) / (4 * h * h)
    return H
se = np.sqrt(np.diag(np.linalg.inv(hess(nll, f.x)))); out["tobit"] = dict(n=len(u), lnD=[round(f.x[1], 3), round(se[1], 3)], lnpx=[round(f.x[2], 3), round(se[2], 3)], lntan_i=[round(f.x[3], 3), round(se[3], 3)], nll=round(f.fun, 2))
print("tobit", out["tobit"])
# broken power law and plateau, on measured hollows with px and inc terms; profile the break
lD = np.log(m.D_m.values); Xo = np.c_[np.log(m.px), np.log(np.tan(np.radians(m.inc)))]; yy = np.log(m.depth_m.values)
def fit_piece(Db, flat=False):
    x1 = np.minimum(lD, np.log(Db)); x2 = np.maximum(lD - np.log(Db), 0)
    XX = np.c_[np.ones(len(yy)), x1, Xo] if flat else np.c_[np.ones(len(yy)), x1, x2, Xo]
    bb, *_ = np.linalg.lstsq(XX, yy, rcond=None); r = yy - XX @ bb; return (r ** 2).sum(), bb
rss0 = ((yy - np.c_[np.ones(len(yy)), lD, Xo] @ np.linalg.lstsq(np.c_[np.ones(len(yy)), lD, Xo], yy, rcond=None)[0]) ** 2).sum()
grid = np.exp(np.linspace(np.log(60), np.log(700), 60)); n = len(yy)
prof_b = [fit_piece(D)[0] for D in grid]; prof_f = [fit_piece(D, True)[0] for D in grid]
ib, iff = int(np.argmin(prof_b)), int(np.argmin(prof_f))
aic = lambda rss, kpar: n * np.log(rss / n) + 2 * kpar
out["piecewise"] = dict(n=n, aic_powerlaw=round(aic(rss0, 5), 2), broken=dict(Db=round(grid[ib]), aic=round(aic(prof_b[ib], 7), 2), coefs=np.round(fit_piece(grid[ib])[1], 3).tolist()),
                        plateau=dict(Ds=round(grid[iff]), aic=round(aic(prof_f[iff], 6), 2), coefs=np.round(fit_piece(grid[iff], True)[1], 3).tolist()))
# profile interval for plateau D_s (delta AIC-like: n ln RSS within 3.84)
ok = [grid[i] for i in range(len(grid)) if n * np.log(prof_f[i] / prof_f[iff]) <= 3.84]; out["piecewise"]["plateau"]["Ds_range95"] = [round(min(ok)), round(max(ok))]
print("piecewise", out["piecewise"])
# manual vs automatic
both = m[m.auto_h.notna() & (m.auto_h > 0)]; r = both.depth_m / both.auto_h
out["vs_auto"] = dict(n=len(both), median_manual_over_auto=round(float(r.median()), 3), frac_auto_low_by_2x=round(float((r > 2).mean()), 3), frac_within_30pct=round(float(((r > 0.77) & (r < 1.3)).mean()), 3),
                      spearman=round(float(stats.spearmanr(both.depth_m, both.auto_h)[0]), 3))
print("vs auto", out["vs_auto"])
out["notes_multi_point"] = b[b.note.str.contains("multi", case=False)].hollow_id.tolist()
out["limits"] = b[b.shadow != "measured"][["hollow_id", "D_m", "cls", "shadow", "depth_m", "h_up", "note"]].round(2).to_dict("records")
out["off_axis_gt_25pct"] = int((m.off_px > 0.25 * m.L_px).sum())
json.dump(out, open("bench_summary.json", "w"), indent=1, default=str)
print("limits", out["limits"]); print("off-axis", out["off_axis_gt_25pct"], "multi", out["notes_multi_point"])
