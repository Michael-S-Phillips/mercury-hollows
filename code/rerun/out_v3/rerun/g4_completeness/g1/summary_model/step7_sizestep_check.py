# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/step7_sizestep_check.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Check (5 Nov 2026, reviewer's flag): does the degradation-state ("age") share in the thinned occurrence model absorb a
size step at 10 km? Herrick et al. (2018) rate only craters > 10 km, so the reference category (unrated) is mostly 5-10 km.
Point estimates of the dominance shares for (a) the paper's model, (b) size group + indicator D >= 10 km, (c) craters >= 10 km
only. Also reports g. Writes out/sizestep_check.json."""
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/step7_sizestep_check.py; unchanged (paths come from common.py: SM_DIR).
import json, numpy as np
import step7_improved as s
d = s.d.copy(); d["ge10"] = (d.D >= 10).astype(float)
out = {}
for name, dd, groups in [("paper", d, {**s.PG, "composition": ["lrm"]}),
                         ("size_plus_ge10", d, {**{k: v for k, v in s.PG.items() if k != "size"}, "size": ["logD", "logD2", "ge10"], "composition": ["lrm"]}),
                         ("D_ge10_only", d[d.D >= 10].reset_index(drop=True), {**s.PG, "composition": ["lrm"]})]:
    y, off, fit, r2_all, dominance = s.make(dd, groups, "cloglog")
    dm = dominance(r2_all(np.ones(len(dd))))
    full = [c for g in groups for c in groups[g]]; b = fit(full, np.ones(len(dd)))[0]
    out[name] = {"n": int(len(dd)), "hosts": int(y.sum()), "R2": round(dm["full"], 4), "share_pct": {g: round(100 * v, 1) for g, v in dm["share"].items()}, "g": round(float(b[1]), 3)}
    print(name, out[name], flush=True)
json.dump(out, open(f"{s.OUT}/sizestep_check.json", "w"), indent=1)
