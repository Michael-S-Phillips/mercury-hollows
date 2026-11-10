# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/count_model_check.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""g1 rerun (6 Nov 2026), new diagnostic for the count model of step7_improved.py. step7 estimates the negative binomial
alpha with statsmodels NegativeBinomial(...).fit(disp=0, maxiter=300) and then runs the dominance analysis with alpha fixed.
For v3b (and for v3 with single-threaded BLAS) that fit does not converge (llf = nan) and leaves alpha at its start value
1.051, so the count-model shares in improved.json are not usable. Here alpha is estimated from Poisson GLM start values
(statsmodels NegativeBinomial with start_params, method 'nm' then 'bfgs'), and the dominance analysis is rerun with
step7_improved.make(..., 'nb', alpha) exactly as in step7. Writes SM_DIR/out/count_model_check.json."""
import os, json, numpy as np, statsmodels.api as sm
import step7_improved as s
d = s.d; groups = {**s.PG, "composition": ["lrm"]}; cols = [c for g in groups for c in groups[g]]
X = sm.add_constant(d[["lnC"] + cols].astype(float)); y = d.n09.values
po = sm.GLM(y, X, family=sm.families.Poisson(), offset=d.lnA.values).fit()
nb = sm.NegativeBinomial(y, X, offset=d.lnA.values)
r0 = nb.fit(start_params=np.r_[po.params.values, 1.0], method="nm", maxiter=20000, disp=0)
r1 = nb.fit(start_params=r0.params.values, method="bfgs", maxiter=2000, disp=0)
alpha = float(r1.params["alpha"])
yy, off, fit, r2_all, dominance = s.make(d, groups, "nb", alpha); p = dominance(r2_all(np.ones(len(d))))
R = dict(alpha=alpha, converged=bool(r1.mle_retvals.get("converged")), llf=float(r1.llf), point=p)
json.dump(R, open(f"{s.OUT}/count_model_check.json", "w"), indent=1)
print(os.environ["CATALOG"], "alpha", round(alpha, 3), "converged", R["converged"], "R2", round(p["full"], 4), {k: round(100 * v, 1) for k, v in p["share"].items()})
