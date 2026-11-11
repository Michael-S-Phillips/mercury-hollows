#!/usr/bin/env python3
"""Depth of hollows in nested host craters (host formed inside an older crater or basin) vs isolated host craters.

Question (Michael, 30 Oct 2026): do hollows in host craters that formed inside an older crater or basin have different
depths from hollows in isolated host craters, with the controls the paper uses for depth (Section S10)?

Inputs (argv[1] = hollows-emissivity root; default ~/mnt/hollows-emissivity):
- manuscripts/insolation_thermal/data/test1_rows_inc.csv: the 2,518 shadow-length depths of Blewett et al. (2016) with pixel scale (px, m), depth (m),
  image (filename), mean annual surface temperature (T, K) and solar incidence (inc, deg) from the MDIS CDR index
  (written by depth_reanalysis.py; identical to the copy in the mercury-hollows repository).
- manuscripts/insolation_thermal/data/test2_craters_h18.csv: Herrick et al. (2018) catalog (the catalog crater_dag.py reads).
- P/analysis_v8/crater_dag_nodes.csv: containment graph, merged Caloris (v8 default run).
Host assignment, as in manuscripts/insolation_thermal/depthctl/depth_structure.py (Section S11) and Section S3: the smallest catalog crater >= 5 km whose
1.1 radii contain the hollow (R_Mercury = 2439.7 km, catalog geometry). Hosts are joined to the graph by catalog id.
Regression, as in depth_reanalysis.py (Section S10, model T_px_inc): ln depth ~ group + 1/T + ln px + ln tan(inc), OLS with
standard errors clustered by image. Depth ratio = exp(b). Sensitivity: SEs clustered by host crater; possible ring
artifacts (crater_dag possible_ring) excluded.
Median corrected depth: depth corrected to the median pixel scale and incidence with lnh ~ lnpx + ltan fitted to all 2,518
(as the band table of Section S10), 95% interval from 2,000 bootstrap resamples of images.
Writes depth_nested_record.json next to this script.
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/mnt/hollows-emissivity")
G = f"{ROOT}/manuscripts/insolation_thermal"
HERE = os.path.dirname(os.path.abspath(__file__))
NODES = os.path.join(os.path.dirname(HERE), "crater_dag_nodes.csv")
F_DEPTH, F_CAT = f"{G}/data/test1_rows_inc.csv", f"{G}/data/test2_craters_h18.csv"

# ---------------------------------------------------------------- depths and controls (depth_reanalysis.py)
d = pd.read_csv(F_DEPTH)
d["x"] = 1.0 / d["T"]; d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["img"] = d.filename.astype(str)
assert np.allclose(d.lnh, np.log(d.depth)) and np.allclose(d.lnpx, np.log(d.px))

# ---------------------------------------------------------------- host assignment (depth_structure.py)
cr = pd.read_csv(F_CAT); cr = cr[cr.diameter >= 5].reset_index(drop=True)
R_M = 2439.7
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
tree = cKDTree(xyz(cr.latitude.values, cr.lon360.values)); Pt = xyz(d.latitude.values, d.longitude.values % 360)
Dkm = cr.diameter.values; host = np.full(len(d), -1); rho = np.full(len(d), np.nan)
for i, cl in enumerate(tree.query_ball_point(Pt, r=0.36)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(Pt[i] @ tree.data[cl].T, -1, 1)) * R_M
    ok = ang <= 1.1 * Dkm[cl] / 2
    if ok.any():
        j = np.argmin(np.where(ok, Dkm[cl], np.inf)); host[i] = cl[j]; rho[i] = ang[j] / (Dkm[cl[j]] / 2)
d["host_id"] = np.where(host >= 0, cr.unique_id.values[np.maximum(host, 0)].astype(str), "")
d["hostD"] = np.where(host >= 0, Dkm[np.maximum(host, 0)], np.nan); d["host_rho"] = rho

# ---------------------------------------------------------------- graph attributes of the host
n = pd.read_csv(NODES, dtype={"id": str, "parent_id": str}, keep_default_na=False, na_values=[""])
n["parent_id"] = n.parent_id.fillna("")
hosts_flag = dict(zip(n.id, n.hosts))
cols = ["id", "kind", "D", "hosts", "parent", "parent_id", "parent_kind", "parent_D", "chain", "zone", "possible_ring", "n_hollows"]
a = d.merge(n[cols].rename(columns={"id": "host_id", "D": "node_D", "hosts": "host_hosts_flag"}), on="host_id", how="left")
assert len(a) == len(d)
R = dict(inputs=dict(depths=F_DEPTH, catalog=F_CAT, graph_nodes=NODES,
                     graph_columns=dict(id="id", nesting="parent (>= 0 nested) / parent_id", n_ancestors="chain",
                                        position="zone (peak rho<0.2, floor 0.2-0.6, wall >=0.6)", parent_kind="parent_kind",
                                        parent_hosts="hosts of the parent row (catalog host flag; NaN for Orgel basins)",
                                        ring_artifact="possible_ring"),
                     depth_columns=dict(depth="depth (m), lnh", pixel_scale="px (m), lnpx", incidence="inc (deg)",
                                        temperature="T (K, mean annual)", image="filename")))
R["counts"] = dict(hollows=int(len(a)), images=int(a.img.nunique()), assigned=int((a.host_id != "").sum()),
                   unassigned=int((a.host_id == "").sum()), unassigned_images=int(a[a.host_id == ""].img.nunique()),
                   assigned_not_in_graph=int(((a.host_id != "") & a.kind.isna()).sum()),
                   host_craters=int(a[a.host_id != ""].host_id.nunique()),
                   host_flag_0_in_catalog=int(a[a.host_id != ""].drop_duplicates("host_id").host_hosts_flag.eq(0).sum()),
                   caloris_as_host=int((a.host_id == "200153").sum()))
s = a[(a.host_id != "") & a.kind.notna()].copy()
s["nested"] = (s.parent >= 0).astype(int)
s["chain3"] = pd.Categorical(np.clip(s.chain, 0, 2).astype(int).astype(str), categories=["0", "1", "2"])
ph = s.parent_id.map(hosts_flag)
s["ptype"] = pd.Categorical(np.where(s.nested == 0, "isolated", np.where(s.parent_kind == "basin", "in_basin",
                            np.where(ph == 1, "parent_hosts", "parent_crater_no_hollows"))),
                            categories=["isolated", "parent_hosts", "parent_crater_no_hollows", "in_basin"])
# median corrected depth (Section S10 band table): to median pixel scale and incidence
c = smf.ols("lnh ~ lnpx + ltan", data=d).fit()
mp, mt = d.lnpx.median(), d.ltan.median()
for z in (d, s):
    z["hc"] = np.exp(z.lnh - c.params["lnpx"] * (z.lnpx - mp) - c.params["ltan"] * (z.ltan - mt))
R["correction"] = dict(b_lnpx=float(c.params["lnpx"]), b_ltan=float(c.params["ltan"]), px_median_m=float(np.exp(mp)),
                       inc_median_deg=float(np.degrees(np.arctan(np.exp(mt)))))
rng = np.random.default_rng(20261030)
def med_ci(sub, nb=2000):
    imgs = sub.img.unique(); g = {i: sub.hc.values[sub.img.values == i] for i in imgs}
    bt = [np.median(np.concatenate([g[i] for i in rng.choice(imgs, len(imgs))])) for _ in range(nb)]
    return [float(sub.hc.median()), float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))]
def groups(sub, col):
    out = {}
    for k, v in sub.groupby(col, observed=True):
        out[str(k)] = dict(hollows=int(len(v)), images=int(v.img.nunique()), host_craters=int(v.host_id.nunique()),
                           median_corrected_depth_m=med_ci(v), median_raw_depth_m=float(v.depth.median()),
                           median_host_D_km=float(v.drop_duplicates("host_id").node_D.median()))
    return out
def fit(sub, term, name, ref=None):
    rhs = f"C({term}, Treatment('{ref}'))" if ref else term
    f = f"lnh ~ {rhs} + x + lnpx + ltan"
    out = dict(formula=f, n=int(len(sub)), n_img=int(sub.img.nunique()), n_host=int(sub.host_id.nunique()), terms={})
    for cl_name, grp in (("image", sub.img), ("host", sub.host_id)):
        m = smf.ols(f, data=sub).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(grp)[0]})
        for k in m.params.index:
            if not k.startswith(term) and not k.startswith("C("): continue
            b, e, p = float(m.params[k]), float(m.bse[k]), float(m.pvalues[k])
            lab = k.split("[T.")[1].rstrip("]") if "[T." in k else k
            out["terms"].setdefault(lab, {})[f"cluster_{cl_name}"] = dict(b=b, se=e, p=p, ratio=float(np.exp(b)),
                                                                          ci95=[float(np.exp(b - 1.96 * e)), float(np.exp(b + 1.96 * e))])
        if cl_name == "image":
            out["x_coef_K"] = [float(m.params["x"]), float(m.bse["x"])]
            out["lnpx"] = float(m.params["lnpx"]); out["ltan"] = float(m.params["ltan"])
    R[name] = out
    return out
s["logD"] = np.log10(s.node_D)
sets = dict(all=s, ge10=s[s.node_D >= 10])
for tag, sub in sets.items():
    R[f"groups_nested_{tag}"] = groups(sub.assign(grp=np.where(sub.nested == 1, "nested", "isolated")), "grp")
    fit(sub, "nested", f"nested_{tag}")
    R[f"groups_parent_{tag}"] = groups(sub, "ptype"); fit(sub, "ptype", f"parent_{tag}", "isolated")
    R[f"groups_chain_{tag}"] = groups(sub, "chain3"); fit(sub, "chain3", f"chain_{tag}", "0")
R["zone_counts_nested_hollows"] = s[s.nested == 1].zone.value_counts().to_dict()
R["parent_kind_counts_nested_hollows"] = s[s.nested == 1].parent_kind.value_counts().to_dict()
R["possible_ring_hollows"] = int(s.possible_ring.astype(bool).sum())   # 0: no ring-artifact sensitivity needed
# sensitivity: host size as an added control (nested hosts are smaller; Section S10 finds no size dependence)
m = smf.ols("lnh ~ nested + logD + x + lnpx + ltan", data=s).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(s.img)[0]})
b, e = m.params["nested"], m.bse["nested"]
R["nested_all_hostlogD"] = dict(ratio=float(np.exp(b)), ci95=[float(np.exp(b - 1.96 * e)), float(np.exp(b + 1.96 * e))],
                                p=float(m.pvalues["nested"]), b_logD=[float(m.params["logD"]), float(m.bse["logD"])])
# parents of nested hosts: hollows, host craters, parent D and hosting flag
pp = s[s.nested == 1].groupby("parent_id").agg(hollows=("lnh", "size"), hosts=("host_id", "nunique"), images=("img", "nunique"),
                                              parent_kind=("parent_kind", "first"), parent_D=("parent_D", "first"), ptype=("ptype", "first"))
R["parents_of_nested_hosts"] = {k: dict(hollows=int(v.hollows), host_craters=int(v.hosts), images=int(v.images), kind=v.parent_kind,
                                        D_km=float(v.parent_D), group=str(v.ptype)) for k, v in pp.sort_values("hollows", ascending=False).iterrows()}
R["note"] = ("Groups with fewer than 3 host craters (>= 2 ancestors: 2 hosts, 14 hollows, 3 images; 1 hollow with hosts >= 10 km) "
             "are not interpretable; their coefficients are recorded but not reported.")
# check against the paper's crater-median table (Section S11: 87 craters with >= 3 measured hollows, 1,894 hollows)
k3 = a[a.host_id != ""].host_id.value_counts(); k3 = k3[k3 >= 3]
R["check_S11_craters_ge3"] = dict(craters=int(len(k3)), hollows=int(k3.sum()))
json.dump(R, open(f"{HERE}/depth_nested_record.json", "w"), indent=1)
def show(name):
    for lab, v in R[name]["terms"].items():
        ci, cj = v["cluster_image"], v["cluster_host"]
        print(f"  {name:22s} {lab:28s} {ci['ratio']:.3f} ({ci['ci95'][0]:.3f}-{ci['ci95'][1]:.3f}) p={ci['p']:.3g} | host-cl p={cj['p']:.3g} ({cj['ci95'][0]:.3f}-{cj['ci95'][1]:.3f})")
print(json.dumps(R["counts"]), R["check_S11_craters_ge3"], R["correction"])
for tag in sets:
    print(tag, R[f"nested_{tag}"]["n"], R[f"nested_{tag}"]["n_img"], R[f"nested_{tag}"]["n_host"], "x", [round(v) for v in R[f"nested_{tag}"]["x_coef_K"]])
    for nm in ("nested", "parent", "chain"):
        show(f"{nm}_{tag}")
        for g_, v in R[f"groups_{nm}_{tag}"].items():
            print(f"     {g_:26s} n={v['hollows']:5d} img={v['images']:4d} hosts={v['host_craters']:3d} med={[round(z,1) for z in v['median_corrected_depth_m']]} medD={v['median_host_D_km']:.0f}")
print("logD ctrl", R["nested_all_hostlogD"])
for k, v in list(R["parents_of_nested_hosts"].items())[:12]: print("  parent", k, v)
print(R["zone_counts_nested_hollows"], R["parent_kind_counts_nested_hollows"], "ring", R["possible_ring_hollows"])
