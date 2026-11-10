#!/usr/bin/env python3
"""Workshop test 8 (2 Nov 2026): do volcanic vents mark a second heat source for hollow hosts?
If intrusions charge some hosts, craters with vents host hollows more often at fixed size and nesting. Under version C,
hollows in vent craters should also differ in depth (gentler gradient above a deep intrusion); under version E, depth
should not differ.
Inputs: data/data/MER_Vent_Locations_Classifications.csv (309 vents; Lat, Lon, Compound, Image Quality, Facula, In Crater);
analysis_v9/crater_dag_nodes.csv (craters and basins with hosting flag, nesting, and the controls of the paper's base model:
logD, abslat, C1 completeness, nsp smooth plains, dT = Delta T_w); GRL_insolation_hollows/data/test1_rows_inc.csv (depths).
Vent assignment: a crater 'has a vent' if a vent lies within 1.0 of its radius (great circle, R = 2439.7 km).
Hosting: logit hosts ~ vent + nested + logD + abslat + C1 + nsp + dT, on craters (kind == crater) >= 5 km with a hosting flag;
also isolated craters only, and craters >= 20 km only.
Depth: each depth assigned to the smallest crater >= 5 km whose 1.1 radii contain it (as depth_nested.py); OLS
ln h ~ host_vent + ln px + ln tan(inc) + pole, SE clustered by image (and by host).
Writes test8_vents_record.json.
"""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.spatial import cKDTree
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); P = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"
v = pd.read_csv(f"{ROOT}/data/data/MER_Vent_Locations_Classifications.csv")
v.columns = [c.strip() for c in v.columns]
# The file's 'Lat' and 'Lon' headers are swapped: 'Lat' spans -179..180 and 'Lon' spans -76..62, and the swapped reading puts
# vents at Praxiteles (27N, 300E), NE Rachmaninoff (36N, 64E) and Tyagaraja (4N, 211E). Use lat = Lon, lon = Lat.
v["vlat"], v["vlon"] = v["Lon"], v["Lat"]
n = pd.read_csv(f"{P}/analysis_v9/crater_dag_nodes.csv", dtype={"id": str}, low_memory=False)
R_M = 2439.7
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
R = {"n_vents": int(len(v)), "header_swap_corrected": True, "vent_lat_range": [float(v.vlat.min()), float(v.vlat.max())],
     "node_lon_range": [float(n.lon.min()), float(n.lon.max())]}
c = n[(n.kind == "crater") & (n.D >= 5) & n.hosts.notna()].reset_index(drop=True).copy()
tree = cKDTree(xyz(c.lat.values, c.lon.values % 360))
V = xyz(v.vlat.values, v.vlon.values % 360)
c["vent"] = 0; vent_host = np.full(len(v), -1)
for i, cl in enumerate(tree.query_ball_point(V, r=0.4)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(V[i] @ tree.data[cl].T, -1, 1)) * R_M
    ok = ang <= 1.0 * c.D.values[cl] / 2
    if ok.any():
        j = cl[np.argmin(np.where(ok, c.D.values[cl], np.inf))]; vent_host[i] = j
c.loc[np.unique(vent_host[vent_host >= 0]), "vent"] = 1
c["nested"] = (c.parent >= 0).astype(int); c["hosts"] = c.hosts.astype(int)
R["vents_in_catalog_craters"] = int((vent_host >= 0).sum()); R["craters_with_vent"] = int(c.vent.sum())
R["hosting_rate"] = {"vent": float(c[c.vent == 1].hosts.mean()), "no_vent": float(c[c.vent == 0].hosts.mean()),
                     "n_vent_craters": int(c.vent.sum()), "n_vent_hosts": int(c[c.vent == 1].hosts.sum())}
R["vent_crater_D_quantiles"] = [float(x) for x in c[c.vent == 1].D.quantile([.1, .5, .9])]
def logit(df, f, name):
    m = smf.logit(f, data=df).fit(disp=0); b, se = m.params["vent"], m.bse["vent"]
    R[name] = dict(n=int(len(df)), n_vent=int(df.vent.sum()), odds_ratio=float(np.exp(b)),
                   lo=float(np.exp(b - 1.96 * se)), hi=float(np.exp(b + 1.96 * se)), p=float(m.pvalues["vent"]))
base = "hosts ~ vent + nested + logD + abslat + C1 + nsp + dT"
logit(c, base, "hosting_all")
logit(c[c.nested == 0], base.replace(" + nested", ""), "hosting_isolated")
logit(c[c.D >= 20], base, "hosting_D_ge20")
# Paper specification (analysis_v9/crater_dag_tests.py): host = >= 1 hollow within 0.9 radii (h09), possible rings dropped,
# controls + log regional hosting rate within 500 km (leave-one-out).
hh = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", usecols=["host_uid", "rho"]).dropna()
hh["host_uid"] = hh.host_uid.astype(int).astype(str); k09 = hh[hh.rho <= 0.9].host_uid.value_counts()
c["h09"] = (c.id.map(k09).fillna(0) > 0).astype(int)
cs = c[~c.possible_ring.astype(bool)].copy()
la_, lo_ = np.radians(cs.lat.values), np.radians(cs.lon.values)
u = np.c_[np.cos(la_) * np.cos(lo_), np.cos(la_) * np.sin(lo_), np.sin(la_)]
nb = cKDTree(u).query_ball_point(u, 2 * np.sin(500 / 2439.4 / 2)); hv = cs.hosts.values
cs["ln_local"] = np.log(np.array([(hv[j].sum() - hv[i] + 0.5) / (len(j) - 1 + 1.0) for i, j in enumerate(nb)]))
R["paper_spec_counts"] = dict(n=int(len(cs)), h09_hosts=int(cs.h09.sum()), vent_craters=int(cs.vent.sum()), vent_h09_hosts=int(cs[cs.vent == 1].h09.sum()))
pspec = "h09 ~ vent + nested + logD + abslat + C1 + nsp + dT + ln_local"
logit(cs, pspec, "paperspec_all")
logit(cs[cs.nested == 0], pspec.replace(" + nested", ""), "paperspec_isolated")
logit(cs[cs.D >= 20], pspec, "paperspec_D_ge20")
logit(cs, pspec.replace(" + ln_local", ""), "paperspec_no_local")
# depth
d = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")
d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["img"] = d.filename.astype(str); d["pole"] = (d.abslat >= 60).astype(int)
cr = n[(n.kind == "crater") & (n.D >= 5)].reset_index(drop=True)
t2 = cKDTree(xyz(cr.lat.values, cr.lon.values % 360)); Pt = xyz(d.latitude.values, d.longitude.values % 360); host = np.full(len(d), -1)
for i, cl in enumerate(t2.query_ball_point(Pt, r=0.36)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(Pt[i] @ t2.data[cl].T, -1, 1)) * R_M; ok = ang <= 1.1 * cr.D.values[cl] / 2
    if ok.any(): host[i] = cl[np.argmin(np.where(ok, cr.D.values[cl], np.inf))]
vent_ids = set(c[c.vent == 1].id)
d["host_id"] = np.where(host >= 0, cr.id.values[np.maximum(host, 0)], "")
s = d[d.host_id != ""].copy(); s["host_vent"] = s.host_id.isin(vent_ids).astype(int)
R["depth_counts"] = dict(assigned=int(len(s)), in_vent_hosts=int(s.host_vent.sum()), vent_hosts=int(s[s.host_vent == 1].host_id.nunique()))
for grp, gcol in [("image", "img"), ("host", "host_id")]:
    m = smf.ols("lnh ~ host_vent + lnpx + ltan + pole", s).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(s[gcol])[0]})
    b, se = m.params["host_vent"], m.bse["host_vent"]
    R[f"depth_ratio_cluster_{grp}"] = dict(ratio=float(np.exp(b)), lo=float(np.exp(b - 1.96 * se)), hi=float(np.exp(b + 1.96 * se)))
R["vent_hosts_with_depths"] = s[s.host_vent == 1].groupby("host_id").agg(n=("depth", "size"), D=("host_id", lambda x: float(cr.set_index("id").loc[x.iloc[0], "D"])), lat=("latitude", "median"), lon=("longitude", "median"), depth_med=("depth", "median")).reset_index().to_dict("records")
json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test8_vents_record.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
