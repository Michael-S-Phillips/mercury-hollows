# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/degradation/kinczyk_match_sensitivity.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""g1 rerun (6 Nov 2026), new script: reproduces analysis_workshop/kinczyk_degradation/match_sensitivity.json (no script was
kept for it): Kinczyk per-class odds ratio of hosting (h09) for all matches and for matches with center offset <= 0.3 radius,
with log D + |lat| and with completeness added. Reads degradation/<CATALOG>/kinczyk_matched.csv; writes match_sensitivity.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
G1 = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, G1); from hostfiles import CAT
O = f"{G1}/degradation/{CAT}"; m = pd.read_csv(f"{O}/kinczyk_matched.csv", dtype={"id": str}, low_memory=False)
k = m[m.id.notna() & m.kclass.between(1, 5)].copy(); k["kc"] = k.kclass.astype(int); k["r"] = k.sep_km / (k.kD / 2)
R = {}
for lab, s in [("all", k), ("r_le_0.3", k[k.r <= 0.3])]:
    for ctl in ["logD + abslat", "logD + abslat + C1"]:
        mm = smf.logit(f"h09 ~ kc + {ctl}", s).fit(disp=0); ci = mm.conf_int().loc["kc"]
        R[f"{lab} | {ctl}"] = [int(len(s)), int(s.h09.sum()), round(float(np.exp(mm.params["kc"])), 2), round(float(np.exp(ci[0])), 2), round(float(np.exp(ci[1])), 2)]
json.dump(R, open(f"{O}/match_sensitivity.json", "w"), indent=1); print(json.dumps(R))
