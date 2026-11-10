# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/step3_dominance.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Summary model, step 3: one logistic model of hosting (h09) with every group of variables (common.py), and general
dominance analysis of McFadden R^2 by group: for each group, its added R^2 averaged over all 2^7 subsets of the other
groups (first within each subset size, then across sizes), so variance shared between groups is split evenly.
The pre-specified combinations and the regional hosting rate are added after the eight groups, and their increments
are reported separately. Intervals: bootstrap over ~640 km spatial blocks (B set by argv[1], default 200).
Also: the joint model's odds ratios, and a check that the paper's nesting specification is reproduced.
Writes out/dominance.json, out/joint_model.csv."""
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/step3_dominance.py; unchanged (paths come from common.py: SM_DIR).
import sys, json, itertools, numpy as np, pandas as pd
from multiprocessing import Pool
from common import load, OUT, GROUPS, INTER, REGION, LABEL, fit_logit, design, ll_null
d = load(); y = d.h09.values.astype(float); G = list(GROUPS); K = len(G)
B = int(sys.argv[1]) if len(sys.argv) > 1 else 200

def r2_all(w):
    """McFadden R^2 for every subset of groups (dict keyed by frozenset), plus +combinations and +regional, for weights w."""
    l0 = ll_null(y, w); out = {}
    for k in range(K + 1):
        for S in itertools.combinations(range(K), k):
            cols = [c for g in S for c in GROUPS[G[g]]]
            out[frozenset(S)] = 1 - fit_logit(design(d, cols), y, w)[1] / l0 if cols else 0.0
    full = [c for g in G for c in GROUPS[g]]
    out["inter"] = 1 - fit_logit(design(d, full + INTER), y, w)[1] / l0
    out["region"] = 1 - fit_logit(design(d, full + INTER + REGION), y, w)[1] / l0
    out["region_only_add"] = 1 - fit_logit(design(d, full + REGION), y, w)[1] / l0
    return out

def dominance(r2):
    gd = {}
    for g in range(K):
        others = [o for o in range(K) if o != g]; per_k = []
        for k in range(K):
            inc = [r2[frozenset(S) | {g}] - r2[frozenset(S)] for S in itertools.combinations(others, k)]
            per_k.append(np.mean(inc))
        gd[G[g]] = float(np.mean(per_k))
    full = r2[frozenset(range(K))]
    alone = {G[g]: float(r2[frozenset([g])]) for g in range(K)}
    last = {G[g]: float(full - r2[frozenset(range(K)) - {g}]) for g in range(K)}
    return dict(general=gd, alone=alone, last=last, full=float(full), inter_add=float(r2["inter"] - full),
                region_add=float(r2["region"] - r2["inter"]), region_add_no_inter=float(r2["region_only_add"] - full))

def boot(seed):
    rng = np.random.default_rng(seed); blocks = d.block.unique()
    pick = rng.choice(blocks, size=len(blocks), replace=True)
    cnt = pd.Series(pick).value_counts(); w = d.block.map(cnt).fillna(0).values.astype(float)
    return dominance(r2_all(w))

if __name__ == "__main__":
    point = dominance(r2_all(None))
    with Pool(2) as p: bs = p.map(boot, range(1, B + 1))
    def ci(path):
        v = np.array([path(b) for b in bs]); return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    R = {"n": int(len(d)), "events": int(y.sum()), "n_blocks": int(d.block.nunique()), "B": B, "point": point, "ci": {}}
    for g in G:
        R["ci"][g] = dict(general=ci(lambda b: b["general"][g]), share=ci(lambda b: b["general"][g] / b["full"]),
                          alone=ci(lambda b: b["alone"][g]), last=ci(lambda b: b["last"][g]))
    R["ci"]["full"] = ci(lambda b: b["full"]); R["ci"]["inter_add"] = ci(lambda b: b["inter_add"]); R["ci"]["region_add"] = ci(lambda b: b["region_add"])
    R["share"] = {g: point["general"][g] / point["full"] for g in G}
    R["rank_stability"] = {g: float(np.mean([sorted(b["general"], key=lambda x: -b["general"][x]).index(g) + 1 for b in bs])) for g in G}
    # joint model odds ratios (all groups + combinations), Wald and block-bootstrap intervals
    cols = [c for g in G for c in GROUPS[g]] + INTER
    X = design(d, cols); b0, ll, cov = fit_logit(X, y)
    rng = np.random.default_rng(99); blocks = d.block.unique(); BB = []
    for _ in range(B):
        cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); w = d.block.map(cnt).fillna(0).values.astype(float)
        BB.append(fit_logit(X, y, w)[0])
    BB = np.array(BB); se = np.sqrt(np.diag(cov))
    jm = pd.DataFrame({"term": ["Intercept"] + cols, "label": ["Intercept"] + [LABEL[c] for c in cols], "coef": b0, "se_wald": se,
                       "OR": np.exp(b0), "OR_lo_wald": np.exp(b0 - 1.96 * se), "OR_hi_wald": np.exp(b0 + 1.96 * se),
                       "OR_lo_block": np.exp(np.percentile(BB, 2.5, axis=0)), "OR_hi_block": np.exp(np.percentile(BB, 97.5, axis=0))})
    jm["group"] = ["-"] + [g for g in G for c in GROUPS[g]] + ["combinations"] * len(INTER)
    jm.to_csv(f"{OUT}/joint_model.csv", index=False)
    # check: paper nesting specification (analysis_v9/crater_dag_tests.py nested_local_h09: logD + abslat + C1 + nsp + dT + ln_local + nested)
    pc = ["logD", "abslat", "C1", "nsp", "dT", "ln_local", "nested"]; bp, _, cp = fit_logit(design(d, pc), y, ridge=0)
    R["check_paper_nesting_OR"] = [float(np.exp(bp[-1])), float(np.exp(bp[-1] - 1.96 * np.sqrt(cp[-1, -1]))), float(np.exp(bp[-1] + 1.96 * np.sqrt(cp[-1, -1])))]
    json.dump(R, open(f"{OUT}/dominance.json", "w"), indent=1)
    print(json.dumps({"full": point["full"], "share": R["share"], "inter_add": point["inter_add"], "region_add": point["region_add"],
                      "check": R["check_paper_nesting_OR"], "rank": R["rank_stability"]}, indent=1))
    print(jm[["label", "OR", "OR_lo_block", "OR_hi_block"]].round(2).to_string(index=False))
