#!/usr/bin/env python3
"""Completeness controls for task (b) of caloris_v9.py (run after it). Writes caloris_v9_controls_record.json."""
import json, os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import fisher_exact
H = os.path.dirname(os.path.abspath(__file__))
c = pd.read_csv(f"{H}/caloris_v9_small_craters.csv", dtype={"id": str, "parent_id": str})
c["nested"] = (c.grp == "nested").astype(int); c["logD"] = np.log10(c.D)
d = c[c.grp == "direct"].copy(); d["outer"] = (d.rho_cal >= 0.8).astype(int)
R = {}
m = smf.logit("anyall ~ nested + C1 + logD", data=c).fit(disp=0)
R["nested_logit"] = {k: [float(m.params[k]), float(m.bse[k]), float(m.pvalues[k])] for k in m.params.index}
R["nested_OR"] = float(np.exp(m.params["nested"]))
m = smf.logit("anyall ~ outer + C1 + logD", data=d).fit(disp=0)
R["outer_logit_direct"] = {k: [float(m.params[k]), float(m.bse[k]), float(m.pvalues[k])] for k in m.params.index}
R["outer_OR"] = float(np.exp(m.params["outer"]))
R["C1_median_outer_inner"] = [float(d[d.outer == 1].C1.median()), float(d[d.outer == 0].C1.median())]
dd = d[d.C1 >= 0.2]; a, b = dd[dd.outer == 1], dd[dd.outer == 0]
R["C1_ge_0p2_outer_vs_inner"] = dict(outer=f"{int(a.anyall.sum())}/{len(a)}", inner=f"{int(b.anyall.sum())}/{len(b)}",
    fisher_p=float(fisher_exact([[a.anyall.sum(), len(a) - a.anyall.sum()], [b.anyall.sum(), len(b) - b.anyall.sum()]])[1]))
ne = c[c.nested == 1]; ne5 = ne[ne.parent_id.isin(["200520", "201095", "201324", "201995", "109885", "111583"])]
R["nested_parents_dex_gt_2p5"] = dict(nested=f"{int(ne5.anyall.sum())}/{len(ne5)}", direct=f"{int(d.anyall.sum())}/{len(d)}",
    fisher_p=float(fisher_exact([[ne5.anyall.sum(), len(ne5) - ne5.anyall.sum()], [d.anyall.sum(), len(d) - d.anyall.sum()]])[1]))
json.dump(R, open(f"{H}/caloris_v9_controls_record.json", "w"), indent=1); print(json.dumps(R, indent=1))
