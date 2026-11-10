#!/usr/bin/env python3
"""Observed vs expected hosts for craters superposed on ghost craters, split by whether excavation reaches below the
lower-bound fill. Expected counts come from a hosting model (h09 ~ logD, logD^2, |lat|, |lat|^2, C1) fitted to craters
off the smooth plains (nsp = 0, outside Caloris, not on a ghost crater). 5 Nov 2026. Writes ghost_expected.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
s = pd.read_csv(f"{HERE}/ghost_superposed.csv", dtype={"id": str})
ref = f[(f.nsp == 0) & (f.in_cal == 0) & (~f.id.isin(s.id))].copy()
m = smf.logit("h09 ~ logD + I(logD**2) + abslat + I(abslat**2) + C1", ref).fit(disp=0)
s["p"] = m.predict(s)
taper = np.where(s.rho < 0.6, 1.0, (1 - s.rho) / 0.4)
R = {"ref_n": int(len(ref)), "ref_hosts": int(ref.h09.sum())}
def summ(x):
    o, e = int(x.h09.sum()), float(x.p.sum())
    lo, hi = (stats.chi2.ppf(0.025, 2 * o) / 2 if o > 0 else 0.0), stats.chi2.ppf(0.975, 2 * o + 2) / 2
    return dict(n=int(len(x)), observed=o, expected=round(e, 2), obs_over_exp=round(o / e, 2) if e > 0 else None,
                ratio_95=[round(lo / e, 2), round(hi / e, 2)], poisson_p_le_obs=round(float(stats.poisson.cdf(o, e)), 4))
for zone, sub in [("interior", s.rho < 0.6), ("all_taper", s.rho < 1.0)]:
    for dz in [0.0, 0.5]:
        x = s[sub].copy(); t = (x.t0 + dz) * (taper[sub.values] if zone == "all_taper" else 1.0)
        pr = x.d_ex > t
        R[f"{zone}_fill+{dz}"] = {"pierce": summ(x[pr]), "no_pierce": summ(x[~pr])}
nsp = f[(f.nsp == 1) & (~f.id.isin(s.id))].copy(); nsp["p"] = m.predict(nsp); R["nsp_not_ghost"] = summ(nsp)
R["all_ghost_superposed"] = summ(s)
json.dump(R, open(f"{HERE}/ghost_expected.json", "w"), indent=1); print(json.dumps(R, indent=1))
