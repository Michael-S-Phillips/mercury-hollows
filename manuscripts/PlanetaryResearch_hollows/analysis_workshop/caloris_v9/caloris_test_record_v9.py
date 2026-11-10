#!/usr/bin/env python3
"""v9-membership version of peak_wall_test/caloris_test.py, written in the same JSON structure (caloris_test_record_v9.json)
so figscripts/make_figs_pr_v9.py can draw Figure 1c from it. Differences from caloris_test.py:
  - Caloris set = the 179 crater descendants (D >= 10 km) of node 200153 in analysis_v9/crater_dag_nodes.csv (was: direct
    children named Caloris in nested_children_classified.csv, 160 craters);
  - isolated set = craters >= 10 km with no parent in the v9 graph, D <= 250 km (was: morph_all_nested parent < 0, not in an
    Orgel basin, D <= 250);
  - completeness = C1 (1-km hollow, at the crater center; nodes table) in place of C_mean, for both sets. The C_cal/C_iso
    columns are medians, as in caloris_test.py;
  - basin-relative position rho_p = distance from the Orgel et al. (2020) Caloris center in Orgel radii.
Hosting: GRL_insolation_hollows/data/host_assignments.csv (same bytes as the file caloris_test.py read); n_in = hollows within
0.9 radii; field = n_in >= 3; any = n_in >= 1. Run after caloris_v9.py (reads caloris_v9_craters_ge10km.csv)."""
import os, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, statsmodels.formula.api as smf
H = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{H}/../../../..")
n = pd.read_csv(f"{ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_v9/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str})
n["parent_id"] = n.parent_id.fillna("")
h = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", usecols=["host_uid", "rho"]).dropna()
h["host_uid"] = h.host_uid.astype(int).astype(str)
n["n_in"] = n.id.map(h[h.rho <= 0.9].host_uid.value_counts()).fillna(0).astype(int)
cal = pd.read_csv(f"{H}/caloris_v9_craters_ge10km.csv", dtype={"id": str, "parent_id": str})
iso = n[(n.kind == "crater") & (n.D >= 10) & (n.parent_id == "") & (n.D <= 250)].copy()
cal["rho_p"] = cal.rho_cal; cal["n_in"] = cal.n09
for d in (cal, iso):
    d["field"] = (d.n_in >= 3).astype(int); d["any"] = (d.n_in >= 1).astype(int); d["logD"] = np.log10(d.D); d["abslat"] = d.lat.abs()
assert len(cal) == 179 and (cal.n_in == cal.id.map(n.set_index("id").n_in)).all()
d_ex = lambda D: 0.1 * 12 ** 0.15 * D ** 0.85
Dthr = {t: (t / (0.1 * 12 ** 0.15)) ** (1 / 0.85) for t in (2.5, 3.0, 3.5)}
cal["set"] = "caloris"; iso["set"] = "isolated"
R = {"Ernst2015_plains_km": [2.5, 3.5], "D_for_excavation_depth": {str(k): round(v, 1) for k, v in Dthr.items()},
     "membership": "v9 crater containment graph: descendants of node 200153 (179 craters >= 10 km); isolated = no parent, D <= 250 km",
     "completeness": "C1 (1-km hollow, crater center) for both sets; C_cal and C_iso are medians", "n_caloris": len(cal), "n_isolated": len(iso)}
bins = [10, 15, 20, 25, 30, 35, 40, 50, 70, 250]
for cmin in (0.0, 0.2):
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        a = cal[(cal.D >= lo) & (cal.D < hi) & (cal.C1 >= cmin)]; b = iso[(iso.D >= lo) & (iso.D < hi) & (iso.C1 >= cmin)]
        rows.append(dict(D=f"{lo}-{hi}", caloris_field=f"{int(a.field.sum())}/{len(a)}", caloris_any=f"{int(a['any'].sum())}/{len(a)}",
                         iso_field=f"{int(b.field.sum())}/{len(b)}", iso_any=f"{int(b['any'].sum())}/{len(b)}",
                         C_cal=round(a.C1.median(), 2) if len(a) else None, C_iso=round(b.C1.median(), 2) if len(b) else None))
    R[f"rates_Cmin{cmin}"] = rows
A = pd.concat([cal, iso])[["field", "any", "logD", "C1", "abslat", "set", "D"]]
A = A[(A.D >= 10) & (A.D <= 100)]
for out in ("field", "any"):
    mo = smf.logit(f"{out} ~ logD*C(set, Treatment('isolated')) + C1 + abslat", data=A).fit(disp=0, maxiter=300)
    R[f"logit_{out}"] = {k: [round(float(mo.params[k]), 3), round(float(mo.bse[k]), 3), round(float(mo.pvalues[k]), 5)] for k in mo.params.index}
    pr = {}
    for s in ("caloris", "isolated"):
        X = pd.DataFrame({"logD": np.log10([15, 25, 35, 45]), "C1": 0.3, "abslat": 32, "set": s}); pr[s] = [round(float(p), 4) for p in mo.predict(X)]
    R[f"pred_{out}_C0.3"] = pr
cal["through"] = (d_ex(cal.D) >= 2.5).astype(int)
mo = smf.logit("field ~ through + C1", data=cal).fit(disp=0)
R["caloris_step_field"] = {k: [round(float(mo.params[k]), 3), round(float(mo.bse[k]), 3), round(float(mo.pvalues[k]), 6)] for k in mo.params.index}
R["caloris_counts"] = {"d_ex<2.5": f"{int(cal[cal.through == 0].field.sum())}/{int((cal.through == 0).sum())}",
                       "d_ex>=2.5": f"{int(cal[cal.through == 1].field.sum())}/{int((cal.through == 1).sum())}"}
s = cal[(cal.D >= 20) & (cal.D < 40)]
R["caloris_20_40_edge_vs_center"] = {"rho_p>=0.7": f"{int(s[s.rho_p >= 0.7].field.sum())}/{int((s.rho_p >= 0.7).sum())}",
                                     "rho_p<0.7": f"{int(s[s.rho_p < 0.7].field.sum())}/{int((s.rho_p < 0.7).sum())}"}
cal.to_csv(f"{H}/caloris_children_v9.csv", index=False)
json.dump(R, open(f"{H}/caloris_test_record_v9.json", "w"), indent=1)
print(json.dumps({k: R[k] for k in ("rates_Cmin0.2", "caloris_counts", "caloris_20_40_edge_vs_center", "D_for_excavation_depth")}))
