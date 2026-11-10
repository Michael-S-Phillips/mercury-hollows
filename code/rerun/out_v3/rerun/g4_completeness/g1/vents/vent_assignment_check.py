# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/vents/vent_assignment_check.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Vent check, part 3 (2 Nov 2026).
(1) Agreement between the catalog's own 'In Crater' flag and whether the vent position (swapped reading: lat = file 'Lon',
    lon = file 'Lat') falls within 1.0 radius of a catalog crater (crater_dag_nodes, kind == crater, D >= 5 km). Compared with
    the as-labeled reading (valid rows only) and positions shifted 100 km in a random direction.
(2) Test 8 paper specification (h09 ~ vent + nested + logD + abslat + C1 + nsp + dT + ln_local, possible rings dropped)
    re-run on subsets of vents: In Crater = yes; image quality good or reasonable; Facula = yes; good quality only;
    and a stricter assignment radius of 0.8 R.
(3) Basemap cutouts (300 km) around the four named faculae with no catalog vent within 50 km: Suge, Madu, Coatl, Bibilava.
Writes vent_assignment_check.json and vent_outlier_faculae.png."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.spatial import cKDTree
# g1 rerun (6 Nov 2026): copy of analysis_workshop/vent_check/vent_assignment_check.py; only input/output paths changed (graph nodes from g1_graph_joint/<CATALOG>/,
# host assignments from hostfiles.HOSTS, outputs in vents/<CATALOG>/).
import sys as _sys; _G1 = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); _sys.path.insert(0, _G1)
from hostfiles import HOSTS, CAT; OUTDIR = f"{_G1}/vents/{CAT}"; os.makedirs(OUTDIR, exist_ok=True)
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); P = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"
HERE = OUTDIR; R_M = 2439.7
v = pd.read_csv(f"{ROOT}/data/data/MER_Vent_Locations_Classifications.csv"); v.columns = [c.strip() for c in v.columns]
v["vlat"], v["vlon"] = v["Lon"], v["Lat"] % 360
n = pd.read_csv(f"{_G1}/{CAT}/crater_dag_nodes.csv", dtype={"id": str}, low_memory=False)
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
OUT = {}

# (1) In Crater agreement
allc = n[n.kind == "crater"].reset_index(drop=True)
T = cKDTree(xyz(allc.lat.values, allc.lon.values % 360)); Dall = allc.D.values
def inside(lat, lon, frac=1.0):
    V = xyz(lat, lon); res = np.zeros(len(V), bool)
    for i, cl in enumerate(T.query_ball_point(V, r=0.4)):
        if cl:
            cl = np.array(cl); ang = np.arccos(np.clip(V[i] @ T.data[cl].T, -1, 1)) * R_M
            res[i] = bool((ang <= frac * Dall[cl] / 2).any())
    return res
flag = (v["In Crater"].str.strip() == "yes").values
ins_s = inside(v.vlat.values, v.vlon.values)
rng = np.random.default_rng(3); th = rng.uniform(0, 2*np.pi, len(v))
dl = 100 / (R_M * np.pi / 180)
ins_n = inside(np.clip(v.vlat.values + dl*np.sin(th), -89, 89), (v.vlon.values + dl*np.cos(th)/np.maximum(np.cos(np.radians(v.vlat.values)), 0.2)) % 360)
ok = v.Lat.abs() <= 90
ins_a = inside(v.Lat[ok].values, v.Lon[ok].values % 360)
def tab(ins, fl):
    return dict(incrater_yes_inside=float(ins[fl].mean()), incrater_no_inside=float(ins[~fl].mean()),
                n_yes=int(fl.sum()), n_no=int((~fl).sum()), agreement=float((ins == fl).mean()))
OUT["in_crater_agreement"] = {"swapped": tab(ins_s, flag), "shifted_100km_null": tab(ins_n, flag),
                              "as_labeled_valid_rows": tab(ins_a, flag[ok.values])}

# (2) Test 8 paper spec on vent subsets
c = n[(n.kind == "crater") & (n.D >= 5) & n.hosts.notna()].reset_index(drop=True).copy()
c["nested"] = (c.parent >= 0).astype(int); c["hosts"] = c.hosts.astype(int)
hh = pd.read_csv(HOSTS, usecols=["host_uid", "rho"]).dropna()
hh["host_uid"] = hh.host_uid.astype(int).astype(str); k09 = hh[hh.rho <= 0.9].host_uid.value_counts()
c["h09"] = (c.id.map(k09).fillna(0) > 0).astype(int)
tree = cKDTree(xyz(c.lat.values, c.lon.values % 360))
def assign(sub, frac=1.0):
    V = xyz(sub.vlat.values, sub.vlon.values); host = np.full(len(V), -1)
    for i, cl in enumerate(tree.query_ball_point(V, r=0.4)):
        if not cl: continue
        cl = np.array(cl); ang = np.arccos(np.clip(V[i] @ tree.data[cl].T, -1, 1)) * R_M
        okk = ang <= frac * c.D.values[cl] / 2
        if okk.any(): host[i] = cl[np.argmin(np.where(okk, c.D.values[cl], np.inf))]
    x = np.zeros(len(c), int); x[np.unique(host[host >= 0])] = 1; return x
cs_mask = ~c.possible_ring.astype(bool).values
cs = c[cs_mask].copy()
la_, lo_ = np.radians(cs.lat.values), np.radians(cs.lon.values)
u = np.c_[np.cos(la_)*np.cos(lo_), np.cos(la_)*np.sin(lo_), np.sin(la_)]
nb = cKDTree(u).query_ball_point(u, 2*np.sin(500/2439.4/2)); hv = cs.hosts.values
cs["ln_local"] = np.log(np.array([(hv[j].sum() - hv[i] + 0.5) / (len(j) - 1 + 1.0) for i, j in enumerate(nb)]))
pspec = "h09 ~ vent + nested + logD + abslat + C1 + nsp + dT + ln_local"
q = v["Image Quality"].str.strip(); fa = v.Facula.str.strip() == "yes"
subsets = {"all_vents": (v, 1.0), "in_crater_yes": (v[flag], 1.0), "quality_good_or_reasonable": (v[q != "poor"], 1.0),
           "facula_yes": (v[fa], 1.0), "quality_good": (v[q == "good"], 1.0),
           "in_crater_yes_and_not_poor": (v[flag & (q != "poor")], 1.0), "all_vents_0.8R": (v, 0.8)}
OUT["test8_paperspec_subsets"] = {}
for name, (sub, frac) in subsets.items():
    cs["vent"] = assign(sub, frac)[cs_mask]
    m = smf.logit(pspec, data=cs).fit(disp=0); b, se = m.params["vent"], m.bse["vent"]
    OUT["test8_paperspec_subsets"][name] = dict(n_vents=int(len(sub)), vent_craters=int(cs.vent.sum()), vent_hosts=int(cs[cs.vent == 1].h09.sum()),
        odds_ratio=round(float(np.exp(b)), 2), lo=round(float(np.exp(b - 1.96*se)), 2), hi=round(float(np.exp(b + 1.96*se)), 2), p=float(m.pvalues["vent"]))
json.dump(OUT, open(f"{HERE}/vent_assignment_check.json", "w"), indent=1); print(json.dumps(OUT, indent=1))

# (3) cutouts of the outlier faculae
import rasterio, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rasterio.windows import Window
r = rasterio.open(f"{ROOT}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif"); px = r.transform.a; RR = 2439400.0
FAC = [("Suge Facula (IAU 26.10N, 300.40W; in Rachmaninoff)", 26.10, 59.60, 32), ("Madu Facula (IAU 21.24S, 74.87W)", -21.24, 285.13, 70),
       ("Coatl Facula (IAU 29.75S, 216.55W)", -29.75, 143.45, 48), ("Bibilava Faculae (IAU 16.40N, 202.80W)", 16.40, 157.20, 100)]
fig, ax = plt.subplots(1, 4, figsize=(16, 4.4))
for a, (lab, la, lo, dkm) in zip(ax, FAC):
    half = 150; lon = ((lo + 180) % 360) - 180
    col = (RR*np.radians(lon) - r.transform.c) / px; row = (r.transform.f - RR*np.radians(la)) / px
    hw = int(half*1000/px/np.cos(np.radians(la))); hh_ = int(half*1000/px)
    img = r.read(1, window=Window(int(col-hw), int(row-hh_), 2*hw, 2*hh_), out_shape=(600, 600)).astype(float)
    lo_, hi_ = np.percentile(img[img > 0], [1, 99.5]) if (img > 0).any() else (0, 1)
    a.imshow(img, cmap="gray", vmin=lo_, vmax=hi_, extent=[-half, half, -half, half])
    a.add_patch(plt.Circle((0, 0), dkm/2, fill=False, ec="cyan", lw=1, ls="--"))
    dx = (((v.vlon - lo + 180) % 360) - 180) * np.cos(np.radians(la)) * RR/1000*np.pi/180; dy = (v.vlat - la) * RR/1000*np.pi/180
    s = (dx.abs() < half) & (dy.abs() < half); a.plot(dx[s], dy[s], "r+", ms=10, mew=1.5)
    for i in np.where(s)[0]: a.text(dx.iloc[i]+4, dy.iloc[i]+4, str(v["Vent ID"].iloc[i]), color="r", fontsize=8)
    a.set_title(lab, fontsize=9); a.set_xlabel("km east"); a.set_ylabel("km north")
fig.suptitle("Named faculae with no catalog vent within 50 km (cyan: IAU center and diameter; red: catalog vents, swapped reading)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{HERE}/vent_outlier_faculae.png", dpi=110)
