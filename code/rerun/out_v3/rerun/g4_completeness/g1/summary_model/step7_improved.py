# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/step7_improved.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Summary model, step 7 (3 Nov 2026): three improvements on step 3.
(1) Detection as thinning, area as exposure. Hollows are treated as a Poisson process over the crater floor whose
    intensity depends on the crater's properties and which the catalog sees with probability tied to completeness C.
    The probability that a crater is observed to host is then 1 - exp(-A r(x) C^g), a complementary log-log model with
    offset ln A (A ~ D^2) and ln(C + 0.01) always in the model (g = 1 would be pure thinning; 18 hosts have C = 0 at the
    center, hence the 0.01). Dominance analysis then splits only the physical groups, and "size" is size per unit area.
(2) The same model for hollow counts (negative binomial, alpha fixed at its full-model value), as a check.
(3) Composition: LRM plus MESSENGER XRS S/Si, Mg/Si, Ca/Si, Al/Si (per standard deviation), on the 27,810 craters with
    all four ratios; the regional hosting rate is added last with and without XRS to see whether composition explains
    the leftover regional clustering.
Also: the step 3 logistic model restricted to well-imaged craters (C >= 0.2).
Intervals: B bootstraps over ~640 km blocks (argv[1], default 100). Writes out/improved.json, out/improved_coefs.csv."""
import sys, json, itertools, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.special import gammaln
from common import load, GROUPS, OUT, fit_logit, design, ll_null
import common
B = int(sys.argv[1]) if len(sys.argv) > 1 else 100
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/step7_improved.py; only input/output paths changed. SM_DIR (env) is the
# folder holding this catalog's crater_features*.csv and out/ (g1_graph_joint/summary_model/<CATALOG>).
import os
_loc = os.path.join(os.environ["SM_DIR"], "crater_features_xrs.csv")
common.DATA = _loc
d = load()
for c in ["s_si", "mg_si", "ca_si", "al_si"]: d[c + "_z"] = (d[c] - d[c].mean()) / d[c].std()
d["lnA"] = 2 * np.log(d.D); d["lnC"] = np.log(d.C1 + 0.01)
PG = {g: GROUPS[g] for g in GROUPS if g != "detection"}            # physical groups; size = logD, logD2 (per unit area now)
XRS = ["s_si_z", "mg_si_z", "ca_si_z", "al_si_z"]

def fit_cloglog(X, y, off, w, iters=60):
    n, p = X.shape; b = np.zeros(p); b[0] = np.log(-np.log(1 - y.mean())) - off.mean()
    def ll(bb):
        e = np.clip(X @ bb + off, -30, 5); mu = -np.expm1(-np.exp(e)); mu = np.clip(mu, 1e-12, 1 - 1e-12)
        return np.sum(w * (y * np.log(mu) + (1 - y) * np.log1p(-mu)))
    cur = ll(b)
    for _ in range(iters):
        e = np.clip(X @ b + off, -30, 5); ee = np.exp(e); mu = np.clip(-np.expm1(-ee), 1e-12, 1 - 1e-12)
        dm = ee * np.exp(-ee); W = w * dm ** 2 / (mu * (1 - mu)); z = (y - mu) / dm
        H = X.T @ (X * W[:, None]) + 1e-6 * np.eye(p); g = X.T @ (W * z); step = np.linalg.lstsq(H, g, rcond=None)[0]; t = 1.0
        while t > 1e-6:
            new = ll(b + t * step)
            if new >= cur - 1e-10: break
            t /= 2
        b = b + t * step; done = abs(new - cur) < 1e-9; cur = new
        if done: break
    return b, cur

def fit_nb(X, y, off, w, alpha, iters=60):
    n, p = X.shape; b = np.zeros(p); b[0] = np.log(np.average(y, weights=w) + 1e-9) - off.mean()
    def ll(bb):
        mu = np.exp(np.clip(X @ bb + off, -30, 20)); r = 1 / alpha
        return np.sum(w * (gammaln(y + r) - gammaln(r) - gammaln(y + 1) + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu))))
    cur = ll(b)
    for _ in range(iters):
        mu = np.exp(np.clip(X @ b + off, -30, 20)); W = w * mu / (1 + alpha * mu); z = (y - mu) / mu
        H = X.T @ (X * W[:, None]) + 1e-6 * np.eye(p); g = X.T @ (W * z); step = np.linalg.lstsq(H, g, rcond=None)[0]; t = 1.0
        while t > 1e-6:
            new = ll(b + t * step)
            if new >= cur - 1e-10: break
            t /= 2
        b = b + t * step; done = abs(new - cur) < 1e-9; cur = new
        if done: break
    return b, cur

def make(dd, groups, family, alpha=None):
    y = (dd.h09 if family == "cloglog" else dd.n09).values.astype(float); off = dd.lnA.values
    G = list(groups); K = len(G)
    def fit(cols, w):
        X = np.column_stack([np.ones(len(dd)), dd.lnC.values] + [dd[c].values.astype(float) for c in cols])
        return (fit_cloglog(X, y, off, w) if family == "cloglog" else fit_nb(X, y, off, w, alpha))
    def r2_all(w):
        l0 = fit([], w)[1]; out = {}
        for k in range(K + 1):
            for S in itertools.combinations(range(K), k):
                cols = [c for g in S for c in groups[G[g]]]
                out[frozenset(S)] = (1 - fit(cols, w)[1] / l0) if cols else 0.0
        full = [c for g in G for c in groups[g]]
        out["region"] = 1 - fit(full + ["ln_local"], w)[1] / l0
        return out
    def dominance(r2):
        full = r2[frozenset(range(K))]; gd = {}
        for g in range(K):
            o = [x for x in range(K) if x != g]
            gd[G[g]] = float(np.mean([np.mean([r2[frozenset(S) | {g}] - r2[frozenset(S)] for S in itertools.combinations(o, k)]) for k in range(K)]))
        return dict(full=float(full), general=gd, share={g: gd[g] / full for g in G}, region_add=float(r2["region"] - full))
    return y, off, fit, r2_all, dominance

ANALYSES = {
    "thinned_presence_all": (lambda: d, {**PG, "composition": ["lrm"]}, "cloglog"),
    "thinned_presence_xrs_subset_lrm_only": (lambda: d[d[XRS].notna().all(axis=1)].reset_index(drop=True), {**PG, "composition": ["lrm"]}, "cloglog"),
    "thinned_presence_xrs_subset_lrm_xrs": (lambda: d[d[XRS].notna().all(axis=1)].reset_index(drop=True), {**PG, "composition": ["lrm"] + XRS}, "cloglog"),
}

def run_boot(args):
    name, seed = args
    get, groups, fam = ANALYSES[name]; dd = get(); y, off, fit, r2_all, dominance = make(dd, groups, fam)
    rng = np.random.default_rng(seed); blocks = dd.block.unique()
    cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); w = dd.block.map(cnt).fillna(0).values.astype(float)
    return name, dominance(r2_all(w))

if __name__ == "__main__":
    R = {"B": B}
    for name, (get, groups, fam) in ANALYSES.items():
        dd = get(); y, off, fit, r2_all, dominance = make(dd, groups, fam)
        R[name] = {"n": int(len(dd)), "events": int(y.sum()), "point": dominance(r2_all(np.ones(len(dd))))}
        cols = [c for g in groups for c in groups[g]]; b, _ = fit(cols, np.ones(len(dd)))
        R[name]["coefs"] = dict(zip(["Intercept", "lnC"] + cols, [float(x) for x in b]))
        print(name, R[name]["n"], R[name]["events"], json.dumps(R[name]["point"]["share"]), "g =", round(b[1], 3), "region add", round(R[name]["point"]["region_add"], 4), flush=True)
    jobs = [(name, s) for name in ANALYSES for s in range(1, B + 1)]
    with Pool(2) as p: res = p.map(run_boot, jobs, chunksize=4)
    for name in ANALYSES:
        bs = [r for n_, r in res if n_ == name]
        q = lambda f: [float(np.percentile([f(b) for b in bs], 2.5)), float(np.percentile([f(b) for b in bs], 97.5))]
        R[name]["ci"] = {g: q(lambda b, g=g: b["share"][g]) for g in bs[0]["share"]}
        R[name]["ci"]["full"] = q(lambda b: b["full"]); R[name]["ci"]["region_add"] = q(lambda b: b["region_add"])
    # regional increment with and without XRS, same bootstrap draws (paired by seed)
    a = [r for n_, r in res if n_ == "thinned_presence_xrs_subset_lrm_only"]; b2 = [r for n_, r in res if n_ == "thinned_presence_xrs_subset_lrm_xrs"]
    diff = [x["region_add"] - y_["region_add"] for x, y_ in zip(a, b2)]
    R["region_add_reduction_by_xrs"] = dict(point=R["thinned_presence_xrs_subset_lrm_only"]["point"]["region_add"] - R["thinned_presence_xrs_subset_lrm_xrs"]["point"]["region_add"],
                                            ci=[float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))])
    # (2) counts, negative binomial with alpha from the full model (statsmodels), point estimates only
    import statsmodels.api as sm
    groups = {**PG, "composition": ["lrm"]}; cols = [c for g in groups for c in groups[g]]
    X = sm.add_constant(d[["lnC"] + cols].astype(float)); nbm = sm.NegativeBinomial(d.n09.values, X, offset=d.lnA.values).fit(disp=0, maxiter=300)
    alpha = float(nbm.params["alpha"]); y, off, fit, r2_all, dominance = make(d, groups, "nb", alpha)
    R["thinned_counts_all"] = {"n": int(len(d)), "alpha": alpha, "point": dominance(r2_all(np.ones(len(d))))}
    print("counts", json.dumps(R["thinned_counts_all"]["point"]["share"]), flush=True)
    # (4) logistic, well-imaged craters only
    import step3_dominance as s3
    keep = (d.C1 >= 0.2).values; s3.d = d[keep].reset_index(drop=True); s3.y = s3.d.h09.values.astype(float)
    pt = s3.dominance(s3.r2_all(None)); R["logistic_C_ge_0.2"] = {"n": int(keep.sum()), "events": int(s3.y.sum()), "point": pt,
                                                                    "share": {g: pt["general"][g] / pt["full"] for g in pt["general"]}}
    print("C>=0.2", json.dumps(R["logistic_C_ge_0.2"]["share"]))
    json.dump(R, open(f"{OUT}/improved.json", "w"), indent=1)
    rows = []
    for name in ANALYSES:
        for k, v in R[name]["coefs"].items(): rows.append(dict(analysis=name, term=k, coef=v, rate_ratio=float(np.exp(v))))
    pd.DataFrame(rows).to_csv(f"{OUT}/improved_coefs.csv", index=False)
    print(json.dumps({k: R[k] for k in R if k not in ANALYSES and k != "thinned_counts_all"}, indent=1))
