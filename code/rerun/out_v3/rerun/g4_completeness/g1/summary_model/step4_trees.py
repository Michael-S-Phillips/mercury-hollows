# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/step4_trees.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Summary model, step 4: gradient-boosted trees (LightGBM) as a check on the logistic model. The trees find nonlinear
effects and combinations without being told. (1) Spatial 5-fold cross-validation over the ~640 km blocks: AUC and log loss
of the logistic model (eight groups; + combinations; + regional rate) and of the trees (without and with the regional rate).
(2) SHAP values of the trees fitted to all craters, summed within each group (|sum| averaged over craters), as a second
ranking of the groups. (3) SHAP interaction values on all hosts plus 6,000 random non-hosts, ranked by mean |interaction|,
for single variables and for pairs of groups. Writes out/trees.json."""
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/step4_trees.py; unchanged (paths come from common.py: SM_DIR).
import json, numpy as np, pandas as pd, lightgbm as lgb, shap, warnings; warnings.filterwarnings("ignore")
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, log_loss
from common import load, OUT, GROUPS, INTER, REGION, fit_logit, design
d = load(); y = d.h09.values
d["age_code"] = d.age.map({"unrated": 0, "filled": 1, "degraded": 2, "standard": 3, "fresh": 4}).astype(int)
d["zone_code"] = d.zone.map({"none": 0, "wall": 1, "floor": 2, "peak": 3}).astype(int)
TG = {"detection": ["C1"], "size": ["logD"], "insolation": ["abslat", "dT"], "age": ["age_code", "rayed"],
      "nesting": ["nested", "zone_code", "parent_hosts", "basin_parent", "ln_zmax"], "cover": ["nsp", "in_cal", "fill_only"],
      "volcanism": ["vent", "parent_vent"], "composition": ["lrm"]}
F = [c for g in TG for c in TG[g]]; G_of = {c: g for g in TG for c in TG[g]}
P = dict(objective="binary", learning_rate=0.03, num_leaves=15, min_child_samples=40, subsample=0.8, subsample_freq=1,
         colsample_bytree=0.8, reg_lambda=1.0, verbose=-1, n_jobs=2)
LOG = {"logistic: 8 groups": [c for g in GROUPS for c in GROUPS[g]]}
LOG["logistic: + combinations"] = LOG["logistic: 8 groups"] + INTER
LOG["logistic: + combinations + regional"] = LOG["logistic: + combinations"] + REGION
gkf = GroupKFold(n_splits=5); R = {"cv": {}}
folds = list(gkf.split(d, y, groups=d.block))
def cv_scores(pred): return dict(auc=float(roc_auc_score(y, pred)), logloss=float(log_loss(y, np.clip(pred, 1e-6, 1 - 1e-6))))
for name, cols in LOG.items():
    pr = np.zeros(len(d))
    for tr, te in folds:
        b = fit_logit(design(d.iloc[tr], cols), y[tr].astype(float))[0]; pr[te] = 1 / (1 + np.exp(-design(d.iloc[te], cols) @ b))
    R["cv"][name] = cv_scores(pr)
best = None
for feats, tag in [(F, "trees"), (F + ["ln_local"], "trees + regional")]:
    for ne in [150, 300, 600]:
        pr = np.zeros(len(d))
        for tr, te in folds:
            m = lgb.LGBMClassifier(n_estimators=ne, **P).fit(d.iloc[tr][feats], y[tr]); pr[te] = m.predict_proba(d.iloc[te][feats])[:, 1]
        R["cv"][f"{tag}: {ne} trees"] = cv_scores(pr)
        if tag == "trees" and (best is None or R["cv"][f"{tag}: {ne} trees"]["logloss"] < R["cv"][f"trees: {best} trees"]["logloss"]): best = ne
R["n_trees"] = best
m = lgb.LGBMClassifier(n_estimators=best, **P).fit(d[F], y)
ex = shap.TreeExplainer(m); sv = ex.shap_values(d[F]); sv = sv[1] if isinstance(sv, list) else sv
grp = {g: np.abs(sv[:, [F.index(c) for c in TG[g]]].sum(axis=1)).mean() for g in TG}; tot = sum(grp.values())
R["shap_group"] = {g: float(v) for g, v in grp.items()}; R["shap_share"] = {g: float(v / tot) for g, v in grp.items()}
R["shap_feature"] = {c: float(np.abs(sv[:, i]).mean()) for i, c in enumerate(F)}
# direction of the main effects in the trees: mean SHAP among craters with the variable "on" (binary) or in the top decile
dirn = {}
for c in F:
    x = d[c].values
    on = x > 0 if set(np.unique(x)) <= {0, 1} else x >= np.quantile(x, 0.9)
    dirn[c] = float(sv[on, F.index(c)].mean()) if on.any() else None
R["shap_mean_when_high_or_on"] = dirn
rng = np.random.default_rng(7); idx = np.r_[np.where(y == 1)[0], rng.choice(np.where(y == 0)[0], 6000, replace=False)]
iv = ex.shap_interaction_values(d.iloc[idx][F]); iv = iv[1] if isinstance(iv, list) else iv
M = np.abs(iv).mean(axis=0); pairs = []
for i in range(len(F)):
    for j in range(i + 1, len(F)): pairs.append((F[i], F[j], G_of[F[i]], G_of[F[j]], float(2 * M[i, j])))
pairs.sort(key=lambda p: -p[4]); R["top_pairs"] = pairs[:15]
R["main_vs_interaction"] = dict(main=float(np.trace(M)), interaction=float(M.sum() - np.trace(M)))
gp = {}
for a, b, ga, gb, v in pairs:
    k = tuple(sorted((ga, gb))); gp[k] = gp.get(k, 0) + v
R["group_pairs"] = sorted([(a, b, v) for (a, b), v in gp.items()], key=lambda p: -p[2])[:12]
json.dump(R, open(f"{OUT}/trees.json", "w"), indent=1)
print(json.dumps({k: R[k] for k in ["cv", "n_trees", "shap_share", "top_pairs", "group_pairs", "main_vs_interaction"]}, indent=1))
