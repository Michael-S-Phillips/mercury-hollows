#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/insolation_v6/tm_vs_dT_check.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Add-on to threshold.py, 2 Nov 2026: does mean annual temperature Tm beat Delta T_w only because, with |lat| linear in
the controls, Tm also absorbs curvature in latitude? Refit with abslat^2 (and a latitude spline) in the controls.
Writes tm_vs_dT_check.json."""
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from common_null import *
c = pd.read_csv(f"{IN}/crater_features.csv", low_memory=False, dtype={"id": str})
c = c.merge(pd.read_csv(f"{PR}/analysis_v9/crater_dag_nodes.csv", low_memory=False, dtype={"id": str})[["id", "Tm"]], on="id")
c["abslat2"] = c.abslat ** 2; c["Tm10"] = c.Tm / 10; c["dT10"] = c.dT / 10
out = {}
for nm, d in (("lat60", c[c.abslat <= 60]), ("global", c)):
    for cn, ctrl in (("abslat", "logD + abslat + C1"), ("abslat_quad", "logD + abslat + abslat2 + C1"), ("abslat_spline", "logD + cr(abslat, df=5) + C1")):
        r = {}
        for v in ("dT10", "Tm10", "dT10 * logD", "Tm10 * logD", "dT10 + Tm10"):
            m = smf.logit(f"h09 ~ {ctrl} + {v}", data=d).fit(disp=0, maxiter=200); r[v] = round(float(m.aic), 1)
            if v == "dT10 + Tm10": r["joint_p"] = {k: float(m.pvalues[k]) for k in ("dT10", "Tm10")}
        out[f"{nm}_{cn}"] = r; print(nm, cn, r)
json.dump(out, open(f"{G2_OUT}/tm_vs_dT_check.json", "w"), indent=1)
