#!/usr/bin/env python3
"""g25 (9 Nov 2026): the joint model of hosting with a composition group up front (Michael's v6.10 comment 32): the fraction of
concentrated low-reflectance material plus the first two principal components of MDIS eight-filter color (crater interior within
0.9 R, hollows masked; g6_remaining/color/out_v3/crater_color.csv). The 181 craters without valid color keep PC values of 0 (the
mean) and a no-color indicator, so all 31,554 craters and 446 hosts stay in the model. Everything else as step7_improved.py
('thinned_presence_all'): thinned occurrence model (cloglog, offset ln A, ln C), dominance shares, regional increment, 100 block
bootstraps (~640 km); plus point estimates for the count model (negative binomial) and the logistic model of well-imaged craters,
and the regional increment with and without XRS on the XRS subset. Run with the summary model's runner:
  SM_DIR=.../summary_model/v3r CATALOG=v3r BUDGET=110 NPROC=3 python3 .../summary_model/cached_pool.py comp_model.py 100
Writes comp.json and comp_coefs.csv here."""
import os, sys, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
SMD = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g4_completeness/g1/summary_model"); sys.path.insert(0, SMD)
import step7_improved as s7
from multiprocessing import Pool
B = int(sys.argv[1]) if len(sys.argv) > 1 else 100
col = pd.read_csv(os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g6_remaining/color/out_v3/crater_color.csv"), dtype={"id": str}, usecols=["id", "ok", "pc1", "pc2"])
d = s7.d.copy(); d["id"] = d.id.astype(str); d = d.merge(col, on="id", how="left")
ok = d.ok.fillna(False).astype(bool).values; d["nocolor"] = (~ok).astype(int)
for k in ["pc1", "pc2"]:
    m, s = d.loc[ok, k].mean(), d.loc[ok, k].std(); d[k + "_z"] = np.where(ok, (d[k] - m) / s, 0.0)
COMP = ["lrm", "pc1_z", "pc2_z", "nocolor"]; GR = {**s7.PG, "composition": COMP}
XRS = s7.XRS
def boot(seed):
    y, off, fit, r2_all, dom = s7.make(d, GR, "cloglog")
    rng = np.random.default_rng(seed); blocks = d.block.unique()
    cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); w = d.block.map(cnt).fillna(0).values.astype(float)
    return dom(r2_all(w))
if __name__ == "__main__":
    y, off, fit, r2_all, dom = s7.make(d, GR, "cloglog")
    jobs = list(range(1, B + 1))
    with Pool(3) as p: res = p.map(boot, jobs, chunksize=1)
    R = {"B": B, "n": int(len(d)), "events": int(y.sum()), "n_nocolor": int(d.nocolor.sum()), "point": dom(r2_all(np.ones(len(d))))}
    cols = [c for g in GR for c in GR[g]]; b, _ = fit(cols, np.ones(len(d))); R["coefs"] = dict(zip(["Intercept", "lnC"] + cols, [float(x) for x in b]))
    q = lambda f: [float(np.percentile([f(x) for x in res], 2.5)), float(np.percentile([f(x) for x in res], 97.5))]
    R["ci"] = {g: q(lambda x, g=g: x["share"][g]) for g in res[0]["share"]}; R["ci"]["full"] = q(lambda x: x["full"]); R["ci"]["region_add"] = q(lambda x: x["region_add"])
    # rate ratios per SD (per unit for lrm and nocolor) with block-bootstrap intervals of the full model coefficients
    rng = np.random.default_rng(21); blocks = d.block.unique(); BB = []
    for _ in range(100):
        cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); w = d.block.map(cnt).fillna(0).values.astype(float); BB.append(fit(cols, w)[0])
    BB = np.array(BB); names = ["Intercept", "lnC"] + cols
    R["rate_ratios"] = {c: [float(np.exp(b[names.index(c)])), float(np.exp(np.percentile(BB[:, names.index(c)], 2.5))), float(np.exp(np.percentile(BB[:, names.index(c)], 97.5)))] for c in COMP}
    # XRS on the XRS subset: regional increment with and without XRS (point)
    sub = d[d[XRS].notna().all(axis=1)].reset_index(drop=True)
    for nm, gg in [("xrs_subset_comp", GR), ("xrs_subset_comp_xrs", {**s7.PG, "composition": COMP + XRS})]:
        _, _, _, r2a, dm = s7.make(sub, gg, "cloglog"); R[nm] = dm(r2a(np.ones(len(sub)))); R[nm]["n"] = int(len(sub))
    R["region_add_reduction_by_xrs"] = R["xrs_subset_comp"]["region_add"] - R["xrs_subset_comp_xrs"]["region_add"]
    # count model and logistic of well-imaged craters, point estimates
    import statsmodels.api as sm
    X = sm.add_constant(d[["lnC"] + cols].astype(float)); nbm = sm.NegativeBinomial(d.n09.values, X, offset=d.lnA.values).fit(disp=0, maxiter=300)
    alpha = float(nbm.params["alpha"]); _, _, _, r2a, dm = s7.make(d, GR, "nb", alpha); R["counts"] = dm(r2a(np.ones(len(d)))); R["counts"]["alpha"] = alpha
    import step3_dominance as s3
    s3.GROUPS = {**s3.GROUPS, "composition": COMP} if hasattr(s3, "GROUPS") else None
    keep = (d.C1 >= 0.2).values; s3.d = d[keep].reset_index(drop=True); s3.y = s3.d.h09.values.astype(float)
    try:
        pt = s3.dominance(s3.r2_all(None)); R["logistic_C_ge_0.2"] = {"n": int(keep.sum()), "events": int(s3.y.sum()), "share": {g: pt["general"][g] / pt["full"] for g in pt["general"]}, "full": pt["full"]}
    except Exception as e: R["logistic_C_ge_0.2"] = {"error": str(e)}
    json.dump(R, open(f"{HERE}/comp.json", "w"), indent=1)
    pd.DataFrame([dict(term=k, coef=v, rate_ratio=float(np.exp(v))) for k, v in R["coefs"].items()]).to_csv(f"{HERE}/comp_coefs.csv", index=False)
    d[["id", "pc1_z", "pc2_z", "nocolor"]].to_csv(f"{HERE}/crater_color_z.csv", index=False)
    print(json.dumps({k: R[k] for k in ["n", "events", "point", "ci", "rate_ratios", "region_add_reduction_by_xrs"]}, indent=1))
