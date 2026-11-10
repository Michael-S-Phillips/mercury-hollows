# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/build_features.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Summary model, step 1 (3 Nov 2026): one table of every crater-level variable tested in paper v5, for a joint model of
hosting. Unit: craters >= 5 km of Herrick et al. (2018) in the v9 containment graph (kind == crater, possible peak rings
dropped): 31,554 craters. Writes crater_features.csv.
Sources:
  analysis_v9/crater_dag_nodes.csv  size, completeness C1, latitude, dT (Delta T_w), northern smooth plains flag (nsp),
                                    filled/preserved, nesting (parent, zone, chain), depth of origin z_max, fill_only, hosts
  GRL_insolation_hollows/data/host_assignments.csv   hollows within 0.9 radii (h09, the paper's host definition)
  GRL_insolation_hollows/data/test2_craters_h18.csv  ejecta code (rayed = contains 'r')
  data/data/MER_Vent_Locations_Classifications.csv   vents (headers swapped: lat = 'Lon', lon = 'Lat'); smallest crater >= 5 km
                                                     whose rim contains the vent (as in test 8)
  data/Concentrated_LRM_560_params_BW.tif            LRM raster (no georeference; read as equirectangular, -180..180 E, 90N at
                                                     top, as in data/von_mises_fit.ipynb; Tolstoj and Rachmaninoff fall on LRM);
                                                     lrm = fraction of pixels > 0 within one crater radius
Caloris: craters descended from node 200153; cal_shallow / cal_deep split at excavation depth 2.5 km
(0.1 D_sc^0.15 D^0.85, D_sc = 12 km; D/1.2 for simple craters)."""
import os, numpy as np, pandas as pd, rasterio, warnings; warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/build_features.py; only input/output paths changed. SM_DIR (env) is the
# folder holding this catalog's crater_features*.csv and out/ (g1_graph_joint/summary_model/<CATALOG>).
# Graph nodes from g1_graph_joint/<CATALOG>/crater_dag_nodes.csv; host assignments from hostfiles.HOSTS.
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); PR = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"
import sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from hostfiles import HOSTS, CAT
HERE = os.environ["SM_DIR"]; os.makedirs(HERE, exist_ok=True); R_M = 2439.4
n = pd.read_csv(f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/{CAT}/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str}, low_memory=False)
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
U = xyz(n.lat.values, n.lon.values % 360); par = n.parent.values.astype(int)
# hollows within 0.9 radii
h = pd.read_csv(HOSTS, usecols=["host_uid", "rho"]).dropna()
h["host_uid"] = h.host_uid.astype(int).astype(str)
n["n09"] = n.id.map(h[h.rho <= 0.9].host_uid.value_counts()).fillna(0).astype(int); n["h09"] = (n.n09 > 0).astype(int)
# rayed
t2 = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test2_craters_h18.csv", usecols=["unique_id", "ejecta"], dtype={"unique_id": str})
n["rayed"] = n.id.map(dict(zip(t2.unique_id, t2.ejecta.fillna("").astype(str).str.contains("r")))).fillna(False).astype(int)
# vents
v = pd.read_csv(f"{ROOT}/data/data/MER_Vent_Locations_Classifications.csv"); v.columns = [c.strip() for c in v.columns]
V = xyz(v.Lon.values, v.Lat.values % 360); tr = cKDTree(U); vent = np.zeros(len(n), int)
for i, cl in enumerate(tr.query_ball_point(V, 0.4)):
    if not cl: continue
    cl = np.array(cl); cl = cl[n.kind.values[cl] == "crater"]
    if not len(cl): continue
    a = np.arccos(np.clip(V[i] @ U[cl].T, -1, 1)) * 2439.7; ok = a <= n.D.values[cl] / 2
    if ok.any(): vent[cl[np.argmin(np.where(ok, n.D.values[cl], np.inf))]] = 1
n["vent"] = vent
n["parent_hosts"] = np.where(par >= 0, n.hosts.fillna(0).values[np.maximum(par, 0)], 0).astype(int)
n["parent_vent"] = np.where(par >= 0, vent[np.maximum(par, 0)], 0)
# Caloris
cal = int(np.where(n.id.values == "200153")[0][0])
def in_cal(i):
    while par[i] >= 0:
        i = par[i]
        if i == cal: return True
    return False
n["in_cal"] = np.array([in_cal(i) for i in range(len(n))]).astype(int)
n["d_ex"] = 0.1 * np.where(n.D > 12, 12 ** 0.15 * n.D ** 0.85, n.D / 1.2)
n["cal_shallow"] = ((n.in_cal == 1) & (n.d_ex < 2.5)).astype(int); n["cal_deep"] = ((n.in_cal == 1) & (n.d_ex >= 2.5)).astype(int)
# LRM fraction within one radius
L = rasterio.open(f"{ROOT}/data/Concentrated_LRM_560_params_BW.tif").read(1) > 0
ny, nx = L.shape; ppd = nx / 360.0
lrm = np.zeros(len(n))
for i, (la, lo, D) in enumerate(zip(n.lat.values, ((n.lon.values + 180) % 360) - 180, n.D.values)):
    rdeg = D / 2 / (R_M * np.pi / 180); cy = (90 - la) * ppd; cx = (lo + 180) * ppd
    ry = max(rdeg * ppd, 0.5); rx = min(max(rdeg * ppd / max(np.cos(np.radians(la)), 0.05), 0.5), nx / 2)
    ys = np.arange(max(int(np.floor(cy - ry)), 0), min(int(np.ceil(cy + ry)) + 1, ny))
    xs = np.arange(int(np.floor(cx - rx)), int(np.ceil(cx + rx)) + 1)
    yy, xx = np.meshgrid(ys, xs, indexing="ij")
    d2 = ((yy + 0.5 - cy) / ry) ** 2 + ((xx + 0.5 - cx) / rx) ** 2; m = d2 <= 1
    if not m.any(): m = d2 == d2.min()
    lrm[i] = L[yy[m], xx[m] % nx].mean()
n["lrm"] = lrm
# age class
n["age"] = np.where(n.filled == "y", "filled", n.preserved.map({"f": "fresh", "s": "standard", "d": "degraded"}).fillna("unrated"))
n["nested"] = (par >= 0).astype(int); n["basin_parent"] = (n.parent_kind == "basin").astype(int)
n["ln_zmax"] = np.log(n.z_max.clip(lower=0.1))
s = n[(n.kind == "crater") & (~n.possible_ring.astype(bool))].copy()
s["hosts"] = s.hosts.astype(int)
# regional hosting rate within 500 km, leave-one-out (as in analysis_v9/crater_dag_tests.py)
Us = xyz(s.lat.values, s.lon.values % 360); nb = cKDTree(Us).query_ball_point(Us, 2 * np.sin(500 / R_M / 2)); hv = s.hosts.values
s["ln_local"] = np.log(np.array([(hv[j].sum() - hv[i] + 0.5) / (len(j) - 1 + 1.0) for i, j in enumerate(nb)]))
cols = ["id", "lat", "lon", "D", "logD", "abslat", "C1", "dT", "nsp", "age", "rayed", "nested", "zone", "chain", "parent_hosts", "basin_parent",
        "ln_zmax", "fill_only", "in_cal", "cal_shallow", "cal_deep", "vent", "parent_vent", "lrm", "ln_local", "hosts", "h09", "n09"]
s[cols].to_csv(f"{HERE}/crater_features.csv", index=False)
print(len(s), "craters;", int(s.h09.sum()), "h09 hosts;", int(s.hosts.sum()), "hosts;", int(s.vent.sum()), "vent craters;", int(s.rayed.sum()), "rayed")
print(s.age.value_counts().to_dict()); print("lrm mean all", round(s.lrm.mean(), 3), "h09 hosts", round(s[s.h09 == 1].lrm.mean(), 3))
