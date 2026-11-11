#!/usr/bin/env python3
"""Robustness tests for the nesting result of crater_dag.py.
1. Stricter hosting: hollows inside 0.9 radii of the crater center (rho from host_assignments.csv), so hollows on the
   parent's surface just outside a child's rim are not counted for the child.
2. Regional context: the hosting rate of other craters within 500 km (leave-one-out) as a covariate.
3. Within-parent comparison: conditional logit with parent groups (peak/floor vs wall children of the same parent).
4. Inheritance: do children of hosting parents host more often?
Usage: python3 crater_dag_tests.py ~/mnt/hollows-emissivity   (after crater_dag.py)
v8: reads/writes CRATER_DAG_OUT if set; coefficients and SEs are recorded unrounded (v6 rounded them to 3 decimals)."""
import sys, os, json
# g1 rerun (6 Nov 2026): copy of analysis_v9/crater_dag_tests.py; only input/output paths changed. CATALOG (env, v3 or v3b) selects
# code/rerun/out_v3/host_assignments_<CATALOG>.csv and inputs/test2_craters_h18_<CATALOG>.csv (hosts recomputed); see hostfiles.py.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from hostfiles import HOSTS
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf
from statsmodels.discrete.conditional_models import ConditionalLogit
ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/mnt/hollows-emissivity")
OUT = os.environ.get("CRATER_DAG_OUT", os.path.dirname(os.path.abspath(__file__))); G = f"{ROOT}/manuscripts/insolation_thermal"
n = pd.read_csv(f"{OUT}/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str})
h = pd.read_csv(HOSTS, usecols=["host_uid", "rho"]).dropna()
h["host_uid"] = h.host_uid.astype(int).astype(str)
k09 = h[h.rho <= 0.9].host_uid.value_counts()
n["n09"] = n.id.map(k09).fillna(0).astype(int); n["h09"] = (n.n09 > 0).astype(int); n["h09_3"] = (n.n09 >= 3).astype(int)
s = n[(n.kind == "crater") & (~n.possible_ring)].copy()
s["nested"] = (s.parent >= 0).astype(int)
s["zone_c"] = pd.Categorical(np.where(s.fill_only, "fill_only", s.zone), categories=["none", "peak", "floor", "wall", "fill_only"])
# regional hosting rate within 500 km, leave-one-out
la, lo = np.radians(s.lat.values), np.radians(s.lon.values)
u = np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
t = cKDTree(u); nb = t.query_ball_point(u, 2 * np.sin(500 / 2439.4 / 2))
hv = s.hosts.values
loc = np.array([(hv[j].sum() - hv[i] + 0.5) / (len(j) - 1 + 1.0) for i, j in enumerate(nb)])
s["ln_local"] = np.log(loc)
R = {}
def fit(formula, d, name, y="hosts"):
    mm = smf.logit(formula.replace("Y ~", f"{y} ~"), data=d).fit(disp=0)
    R[name] = dict(n=int(mm.nobs), events=int(d[y].sum()), aic=float(mm.aic),
                   params={k: [float(v), float(mm.bse[k]), float(mm.pvalues[k])] for k, v in mm.params.items()})
    return mm
ctrl = "logD + abslat + C1 + nsp + dT"
for y in ["hosts", "h09", "h09_3"]:
    fit(f"Y ~ {ctrl} + nested", s, f"nested_{y}", y)
    fit(f"Y ~ {ctrl} + ln_local + nested", s, f"nested_local_{y}", y)
    fit(f"Y ~ {ctrl} + ln_local + C(zone_c)", s, f"zone_local_{y}", y)
    fit(f"Y ~ {ctrl} + ln_local + C(chain.clip(upper=2))", s, f"chain_local_{y}", y)
# parent kind and size: nesting in a basin vs in a crater
s["pk"] = pd.Categorical(np.where(s.nested == 0, "none", np.where(s.parent_kind == "basin", "basin", np.where(s.parent_D >= 300, "crater_ge300", "crater_lt300"))),
                         categories=["none", "crater_lt300", "crater_ge300", "basin"])
fit(f"Y ~ {ctrl} + ln_local + C(pk)", s, "parentkind_local_h09", "h09")
# inheritance: parent hosts (parents that are cataloged craters)
pid = dict(zip(n.id, n.hosts))
c = s[(s.nested == 1) & (s.parent_kind == "crater")].copy()
c["parent_hosts"] = c.parent_id.map(pid).fillna(0).astype(int)
fit(f"Y ~ {ctrl} + ln_local + parent_hosts", c, "inherit_h09", "h09")
R["inherit_rates_h09"] = {str(k): f"{int(v.h09.sum())}/{len(v)}" for k, v in c.groupby("parent_hosts")}
# within-parent: conditional logit on children, groups = parent; zone effect relative to wall
w = s[(s.nested == 1) & (~s.fill_only) & (s.zone.isin(["peak", "floor", "wall"]))].copy()
grp = w.groupby("parent_id").h09.agg(["sum", "count"])
keep = grp[(grp["sum"] > 0) & (grp["sum"] < grp["count"])].index
w = w[w.parent_id.isin(keep)]
X = pd.DataFrame({"inner": w.zone.isin(["peak", "floor"]).astype(float), "logD": w.logD, "C1": w.C1})
cl = ConditionalLogit(w.h09.values, X.values, groups=w.parent_id.values).fit(disp=0)
R["within_parent_clogit_h09"] = dict(parents=int(len(keep)), children=int(len(w)), events=int(w.h09.sum()),
    params={nm: [float(b), float(se), float(p)] for nm, b, se, p in zip(X.columns, cl.params, cl.bse, cl.pvalues)},
    rates={"inner": f"{int(w[X.inner == 1].h09.sum())}/{int((X.inner == 1).sum())}", "wall": f"{int(w[X.inner == 0].h09.sum())}/{int((X.inner == 0).sum())}"})
json.dump(R, open(f"{OUT}/crater_dag_tests_record.json", "w"), indent=1)
for k, v in R.items():
    if isinstance(v, dict) and "params" in v:
        ps = {kk: vv for kk, vv in v["params"].items() if kk not in ("Intercept", "abslat", "C1", "nsp", "dT")}
        print(k, "n", v.get("n", v.get("children")), "events", v["events"], "AIC", round(v.get("aic", float("nan")), 1), ps)
    else:
        print(k, v)
