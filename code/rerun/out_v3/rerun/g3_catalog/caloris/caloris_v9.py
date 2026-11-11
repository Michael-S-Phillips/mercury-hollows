# g3_catalog copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/caloris_v9/caloris_v9.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Caloris penetration test on the v9 crater containment graph (2 Nov 2026).
Reruns peak_wall_test/caloris_test.py with membership from analysis_v9/crater_dag_nodes.csv: every crater descendant of
the merged Caloris node 200153 (Orgel center and diameter; footprint = union of Orgel and catalog circles). Isolated =
craters with no parent in the v9 graph, D <= 250 km. Hosting from insolation_thermal/data/host_assignments.csv
(identical bytes to the code/rerun/latitude/aspect copy used by v5). Completeness: C1 (1-km hollow, at crater center,
nodes table) is primary because C_mean (morph_all.csv) is missing for 47 isolated craters 10-40 km (mostly at lon ~ 0/360);
C_mean results are given as a sensitivity on the craters that have it.
Also: (b) small craters inside Caloris, direct children vs nested in an intermediate crater; (c) Stage 3 candidates.
Usage: python3 caloris_v9.py ~/mnt/hollows-emissivity"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")
CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/caloris/out_{CAT}"
# G3_GRID=v3r (added 9 Nov 2026): C1 from the refit completeness grid, read from the v3r crater-graph nodes table
# (g4_completeness/g1/v3r/crater_dag_nodes.csv; identical to inputs/crater_dag_nodes_v3.csv except C1); output in out_v3r.
GRID = _os.environ.get("G3_GRID", "")
if GRID: ODIR = f"{G3}/caloris/out_{CAT}{GRID[2:] if GRID.startswith(CAT) else GRID}"
_os.makedirs(ODIR, exist_ok=True)
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/insolation_thermal/data/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
sys.setrecursionlimit(20000)
ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/mnt/hollows-emissivity")
P = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"; G = f"{ROOT}/manuscripts/insolation_thermal"
OUT = ODIR
NODES = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g4_completeness/g1/v3r/crater_dag_nodes.csv")
         if GRID == "v3r" else f"{G3}/inputs/crater_dag_nodes_{CAT}.csv")
n = pd.read_csv(NODES, dtype={"id": str, "parent_id": str})
n["parent_id"] = n.parent_id.fillna("")
m = pd.read_csv(f"{ROOT}/peak_wall_test/morph_all.csv"); m["id"] = m.uid.astype(str)
n = n.merge(m[["id", "C_mean"]], on="id", how="left")
names = pd.read_csv(f"{G}/data/test2_craters_h18.csv", usecols=["unique_id", "name"])
names["id"] = names.unique_id.astype(str); n = n.merge(names[["id", "name"]], on="id", how="left")
n["name"] = n.name.fillna("").astype(str).str.strip()
h = pd.read_csv(HOSTS, usecols=["host_uid", "rho"]).dropna(subset=["host_uid"])
h["host_uid"] = h.host_uid.astype(int).astype(str)
n["n_all"] = n.id.map(h.host_uid.value_counts()).fillna(0).astype(int)
n["n09"] = n.id.map(h[h.rho <= 0.9].host_uid.value_counts()).fillna(0).astype(int)
n["field"] = (n.n09 >= 3).astype(int); n["any09"] = (n.n09 >= 1).astype(int); n["anyall"] = (n.n_all >= 1).astype(int)
assert (n.anyall == (n.n_hollows > 0)).all()
d_ex = lambda D: np.where(D > 12, 0.1 * 12 ** 0.15 * np.asarray(D) ** 0.85, 0.1 * np.asarray(D) / 1.2)
n["d_ex"] = d_ex(n.D.values); n["through"] = (n.d_ex >= 2.5).astype(int)
# Caloris membership: follow parent links up the chain
idx = {k: i for i, k in enumerate(n.id)}; par = n.parent_id.values; memo = {}
def inside(i):
    if i in memo: return memo[i]
    p = par[i]; r = False if p == "" else (True if p == "200153" else inside(idx[p])); memo[i] = r; return r
n["in_cal"] = [inside(i) for i in range(len(n))]
cal0 = n.loc[n.id == "200153"].iloc[0]
def gc(lat1, lon1, lat2, lon2):
    a, b, c_, d = map(np.radians, (lat1, lon1, lat2, lon2))
    return 2439.4 * np.arccos(np.clip(np.sin(a) * np.sin(c_) + np.cos(a) * np.cos(c_) * np.cos(b - d), -1, 1))
n["rho_cal"] = gc(n.lat, n.lon, cal0.lat, cal0.lon) / (cal0.D / 2)
f = lambda d, col="field": f"{int(d[col].sum())}/{len(d)}"
def fisher(a, b, col="field"):
    return float(fisher_exact([[a[col].sum(), len(a) - a[col].sum()], [b[col].sum(), len(b) - b[col].sum()]])[1])
R = {"settings": dict(caloris_center_lat_lon_D=[float(cal0.lat), float(cal0.lon), float(cal0.D)],
                      host_assignments=HOSTS, completeness_primary="C1", nodes=NODES)}

# ------------------------------------------------------------------ (a)
cr = n[(n.kind == "crater") & (n.D >= 10)]
cal = cr[cr.in_cal].copy(); iso = cr[(cr.parent_id == "") & (cr.D <= 250)].copy()
R["a_membership"] = dict(caloris_ge10=len(cal), direct_children=int((cal.parent_id == "200153").sum()),
                         nested_deeper=int((cal.parent_id != "200153").sum()), possible_ring=int(cal.possible_ring.sum()),
                         beyond_orgel_rim_rho_ge1=int((cal.rho_cal >= 1).sum()), isolated=len(iso),
                         old_v5_ids_not_in_v9=None)
old = pd.read_csv(f"{ROOT}/peak_wall_test/caloris_children.csv"); old["id"] = old.uid.astype(str)
R["a_membership"]["old_v5_ids_not_in_v9"] = {i: n.loc[n.id == i, "parent_id"].iat[0] for i in old.id[~old.id.isin(cal.id)]}
R["a_membership"]["v9_not_in_v5"] = int((~cal.id.isin(old.id)).sum())
sh, dp = cal[cal.through == 0], cal[cal.through == 1]
R["a_step"] = dict(shallow=f(sh), deep=f(dp), fisher_p=fisher(sh, dp))
R["a_step_rho_lt1"] = dict(shallow=f(sh[sh.rho_cal < 1]), deep=f(dp[dp.rho_cal < 1]), fisher_p=fisher(sh[sh.rho_cal < 1], dp[dp.rho_cal < 1]))
s12 = sh[sh.any09 == 1]
R["a_shallow_with_1_2"] = s12[["id", "name", "D", "d_ex", "n09", "n_all", "rho_cal", "parent_id"]].round(3).to_dict("records")
R["a_shallow_any_all_rho"] = f(sh, "anyall")
for C in ("C1", "C_mean"):
    a_, b_ = cal[cal[C].notna()], iso[iso[C].notna()]
    r = {"iso_missing": int(iso[C].isna().sum())}
    A = a_[(a_.D >= 10) & (a_.D < 30) & (a_[C] >= 0.2)]; B = b_[(b_.D >= 10) & (b_.D < 30) & (b_[C] >= 0.2)]
    rate = B.field.sum() / len(B)
    r["iso_10_30_C02"] = f(B); r["cal_10_30_C02"] = f(A); r["expected_cal"] = round(float(rate * len(A)), 2)
    r["cal_10_30_C02_field_ids"] = A[A.field == 1][["id", "D", "d_ex"]].round(2).to_dict("records")
    for lo, hi in [(25, 30), (30, 35), (35, 40), (40, 50), (28.5, 40)]:
        A = a_[(a_.D >= lo) & (a_.D < hi) & (a_[C] >= 0.2)]; B = b_[(b_.D >= lo) & (b_.D < hi) & (b_[C] >= 0.2)]
        r[f"{lo}-{hi}"] = dict(cal=f(A), iso=f(B), C_cal_median=round(float(A[C].median()), 3), C_iso_median=round(float(B[C].median()), 3),
                               C_cal_mean=round(float(A[C].mean()), 3), C_iso_mean=round(float(B[C].mean()), 3), fisher_p=fisher(A, B))
    # binned rates as in v5 record
    bins = [10, 15, 20, 25, 30, 35, 40, 50, 70, 250]
    r["rates_C02"] = [dict(D=f"{lo}-{hi}", cal=f(a_[(a_.D >= lo) & (a_.D < hi) & (a_[C] >= 0.2)]),
                           iso=f(b_[(b_.D >= lo) & (b_.D < hi) & (b_[C] >= 0.2)])) for lo, hi in zip(bins[:-1], bins[1:])]
    R[f"a_{C}"] = r
s = cal[(cal.D >= 20) & (cal.D < 40)]
R["a_edge_20_40"] = dict(outer_rho_ge0p7=f(s[s.rho_cal >= 0.7]), inner=f(s[s.rho_cal < 0.7]),
                         fisher_p=fisher(s[s.rho_cal >= 0.7], s[s.rho_cal < 0.7]))
cal.to_csv(f"{OUT}/caloris_v9_craters_ge10km.csv", index=False)

# ------------------------------------------------------------------ (b)
c5 = n[(n.kind == "crater") & (n.D >= 5) & n.in_cal & (n.d_ex < 2.5)].copy()
c5["grp"] = np.where(c5.parent_id == "200153", "direct", "nested")
dr, ne = c5[c5.grp == "direct"], c5[c5.grp == "nested"]
pr = n.set_index("id").loc[sorted(set(ne.parent_id))]
R["b"] = dict(n_direct=len(dr), n_nested=len(ne),
              any_all=dict(direct=f(dr, "anyall"), nested=f(ne, "anyall"), fisher_p=fisher(ne, dr, "anyall")),
              any09=dict(direct=f(dr, "any09"), nested=f(ne, "any09"), fisher_p=fisher(ne, dr, "any09")),
              parents=pr.reset_index()[["id", "name", "D", "d_ex", "n_all", "n09", "parent_id", "C1"]].round(3).to_dict("records"),
              n_parents=len(pr), parents_all_dex_gt_2p5=bool((pr.d_ex > 2.5).all()),
              all_nested_parents_are_craters=bool((pr.kind == "crater").all()))
R["b"]["nested_children"] = ne[["id", "D", "parent_id", "n_all", "n09", "C1", "rho"]].round(3).to_dict("records")
hs, nh = dr[dr.anyall == 1], dr[dr.anyall == 0]
R["b"]["direct_hosts_outer_fifth"] = dict(hosts=f"{int((hs.rho_cal >= 0.8).sum())}/{len(hs)}",
                                          all_direct=f"{int((dr.rho_cal >= 0.8).sum())}/{len(dr)}",
                                          rate_outer=f(dr[dr.rho_cal >= 0.8], "anyall"), rate_inner=f(dr[dr.rho_cal < 0.8], "anyall"),
                                          fisher_p=fisher(dr[dr.rho_cal >= 0.8], dr[dr.rho_cal < 0.8], "anyall"),
                                          hosts_beyond_orgel_rim=int((hs.rho_cal >= 1).sum()))
R["b"]["direct_completeness"] = {C: dict(hosts_median=round(float(hs[C].median()), 3), nonhosts_median=round(float(nh[C].median()), 3),
                                         hosts_mean=round(float(hs[C].mean()), 3), nonhosts_mean=round(float(nh[C].mean()), 3),
                                         n_hosts=int(hs[C].notna().sum()), n_nonhosts=int(nh[C].notna().sum()),
                                         mwu_p=float(mannwhitneyu(hs[C].dropna(), nh[C].dropna()).pvalue))
                                 for C in ("C1",)}
R["b"]["nested_completeness_C1"] = dict(median=round(float(ne.C1.median()), 3), mean=round(float(ne.C1.mean()), 3))
R["b"]["direct_C1_ge_nested_min"] = None
# completeness-restricted comparison: direct children with C1 >= min C1 of nested hosts
cmin = float(ne[ne.anyall == 1].C1.min()); d2 = dr[dr.C1 >= cmin]
R["b"]["direct_C1_ge_nested_min"] = dict(C1_min=round(cmin, 3), direct=f(d2, "anyall"), fisher_p=fisher(ne, d2, "anyall"))
R["b"]["direct_hosts"] = hs[["id", "name", "D", "d_ex", "n_all", "n09", "rho_cal", "C1"]].sort_values("D").round(3).to_dict("records")
c5.to_csv(f"{OUT}/caloris_v9_small_craters.csv", index=False)

# ------------------------------------------------------------------ (c)
ch = n[(n.parent_id != "") & (n.kind == "crater")]
chh = ch[(ch.anyall == 1) & (~ch.possible_ring)]
k = chh.groupby("parent_id").agg(n_host_children=("id", "size"), n_hollows_children=("n_all", "sum"),
                                 max_child_D=("D", "max"))
ring_drop = ch[(ch.anyall == 1) & ch.possible_ring].parent_id.unique()
P_ = n.set_index("id").join(k, how="inner")
P_ = P_[P_.kind == "crater"]
cand = P_[P_.hosts == 0].copy()
both = P_[P_.hosts == 1]
cand["n_children"] = cand.index.map(ch.parent_id.value_counts()).fillna(0).astype(int)
cand = cand.reset_index().rename(columns={"index": "id"})
cols = ["id", "name", "lat", "lon", "D", "C1", "C_mean", "n_host_children", "n_hollows_children", "max_child_D", "n_children",
        "parent_id", "chain", "nsp", "dT", "preserved", "filled", "cent"]
cand = cand[cols].sort_values("C1", ascending=False)
cand.to_csv(f"{OUT}/stage3_candidate_parents.csv", index=False)
R["c"] = dict(parents_with_hosting_child=len(P_), candidates_hosts0=len(cand), parents_also_hosting=len(both),
              ring_children_excluded_parents=[str(x) for x in ring_drop],
              D_quantiles=cand.D.quantile([0, .25, .5, .75, 1]).round(1).tolist(),
              D_bins=pd.cut(cand.D, [10, 20, 40, 80, 160, 2000]).value_counts().sort_index().astype(int).rename(str).to_dict(),
              C1_quantiles=cand.C1.quantile([0, .25, .5, .75, 1]).round(3).tolist(),
              C1_bins=pd.cut(cand.C1, [-0.01, 0.2, 0.4, 0.6, 0.8, 1.0]).value_counts().sort_index().astype(int).rename(str).to_dict(),
              n_C1_ge_0p5=int((cand.C1 >= 0.5).sum()), n_C1_ge_0p2=int((cand.C1 >= 0.2).sum()),
              both_C1_median=round(float(both.C1.median()), 3), cand_C1_median=round(float(cand.C1.median()), 3),
              both_D_median=round(float(both.D.median()), 1),
              top10=cand.head(10)[["id", "name", "lat", "lon", "D", "C1", "n_host_children"]].round(3).to_dict("records"))
json.dump(R, open(f"{OUT}/caloris_v9_record.json", "w"), indent=1, default=str)
print(json.dumps(R, indent=1, default=str))
