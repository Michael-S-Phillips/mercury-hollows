#!/usr/bin/env python3
"""Sensitivity of the nesting and inheritance results (run after crater_dag.py and crater_dag_tests.py).
v8: reads/writes CRATER_DAG_OUT if set."""
import sys, os, json
# g1 rerun (6 Nov 2026): copy of analysis_v9/crater_dag_sens.py; only input/output paths changed. CATALOG (env, v3 or v3b) selects
# code/rerun/out_v3/host_assignments_<CATALOG>.csv and inputs/test2_craters_h18_<CATALOG>.csv (hosts recomputed); see hostfiles.py.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from hostfiles import HOSTS
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf
ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/mnt/hollows-emissivity")
OUT = os.environ.get("CRATER_DAG_OUT", os.path.dirname(os.path.abspath(__file__))); G = f"{ROOT}/manuscripts/GRL_insolation_hollows"
n = pd.read_csv(f"{OUT}/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str})
h = pd.read_csv(HOSTS, usecols=["host_uid", "rho"]).dropna(); h["host_uid"] = h.host_uid.astype(int).astype(str)
n["h09"] = n.id.map(h[h.rho <= 0.9].host_uid.value_counts()).fillna(0).gt(0).astype(int)
s = n[(n.kind == "crater") & (~n.possible_ring)].copy(); s["nested"] = (s.parent >= 0).astype(int)
la, lo = np.radians(s.lat.values), np.radians(s.lon.values); u = np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
nb = cKDTree(u).query_ball_point(u, 2*np.sin(500/2439.4/2)); hv = s.hosts.values
s["ln_local"] = np.log(np.array([(hv[j].sum()-hv[i]+0.5)/(len(j)-1+1.0) for i, j in enumerate(nb)]))
pid = dict(zip(n.id, n.h09)); s["parent_h09"] = s.parent_id.map(pid).fillna(0).astype(int)
s["pclass"] = s.preserved.fillna("").replace({"": "u"})
R = {}
def fit(f, d, name):
    try:
        m = smf.logit(f, data=d).fit(disp=0)
    except Exception as e:
        R[name] = str(e); return
    keys = [k for k in m.params.index if k.startswith(("nested", "parent_h09"))]
    R[name] = dict(n=int(m.nobs), events=int(d.h09.sum()), **{k: dict(OR=round(float(np.exp(m.params[k])), 2),
               ci=[round(float(np.exp(m.params[k]-1.96*m.bse[k])), 2), round(float(np.exp(m.params[k]+1.96*m.bse[k])), 2)],
               p=float(m.pvalues[k])) for k in keys})
ctrl = "logD + abslat + C1 + nsp + dT + ln_local"
fit(f"h09 ~ {ctrl} + nested", s, "all")
fit(f"h09 ~ {ctrl} + nested", s[s.C1 >= 0.5], "C1_ge_0.5")
fit(f"h09 ~ {ctrl} + nested", s[(s.D >= 10) & (s.D < 50)], "D_10_50")
fit(f"h09 ~ {ctrl} + nested + C(pclass)", s[(s.D >= 10) & s.pclass.isin(["f", "s", "d"])], "D_ge10_preservation_controlled")
fit(f"h09 ~ {ctrl} + nested", s[s.nsp == 0], "outside_northern_plains")
fit(f"h09 ~ {ctrl} + nested", s[(s.parent < 0) | (s.parent_D < 300)], "parents_lt300_only")
fit(f"h09 ~ {ctrl} + nested", s[s.abslat <= 60], "abslat_le60")
c = s[(s.nested == 1) & (s.parent_kind == "crater")]
fit(f"h09 ~ {ctrl} + parent_h09", c, "inherit_all")
fit(f"h09 ~ {ctrl} + parent_h09", c[c.C1 >= 0.5], "inherit_C1_ge_0.5")
fit(f"h09 ~ {ctrl} + parent_h09 + C(pclass)", c[(c.D >= 10) & c.pclass.isin(["f", "s", "d"])], "inherit_preservation_controlled")
# matched check: nested vs isolated craters within 300 km, same size within 25%, completeness within 0.15
iso = s[s.nested == 0]; ui = u[(s.nested == 0).values]; ti = cKDTree(ui)
ne = s[s.nested == 1]; une = u[(s.nested == 1).values]
pairs = []
for (idx, row), uu in zip(ne.iterrows(), une):
    j = ti.query_ball_point(uu, 2*np.sin(300/2439.4/2))
    if not j: continue
    cand = iso.iloc[j]; cand = cand[(abs(np.log(cand.D/row.D)) < np.log(1.25)) & (abs(cand.C1-row.C1) < 0.15)]
    if len(cand): pairs.append((row.h09, cand.h09.mean(), len(cand)))
P = np.array(pairs)
R["matched_300km"] = dict(nested_with_matches=int(len(P)), nested_rate=float(P[:, 0].mean()), matched_isolated_rate=float(P[:, 1].mean()),
                          ratio=float(P[:, 0].mean()/P[:, 1].mean()))
json.dump(R, open(f"{OUT}/crater_dag_sens_record.json", "w"), indent=1)
print(json.dumps(R, indent=1))
