#!/usr/bin/env python3
"""g23 step 11 (9 Nov): the paper's latitude and color depth regressions on automatic depths (depth_all.csv), with and without hollow size,
with and without the province, beside the same models on the 2,518 Blewett et al. (2016) depths (test1_rows_inc.csv).
- Clean automatic depths: depth > 0 and no pilot flag no_shadow, edge_sun, frame_edge (as dvs_analyze.py).
- T: mean annual surface temperature from the paper's thermal grid (GRL_insolation_hollows thermal.lookup). Poleward = lat >= 60N.
- Color: the host crater's interior color (crater_color.csv pc1 = spectral slope, pc2 = curvature; within 0.9 R, hollows masked), z-scored
  over the rows of each data set. The paper's depth-color result uses the surroundings within 10 km of each depth point; that map
  extraction is not rerun here, so the Blewett depths are also fitted with host color as a bridge.
- OLS of ln depth, standard errors clustered by image, as in the paper. Writes depth_regress.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings; warnings.filterwarnings("ignore")
HE = os.path.expanduser("~/mnt/hollows-emissivity"); GRL = f"{HE}/manuscripts/GRL_insolation_hollows"
sys.path.insert(0, f"{GRL}/scripts"); from thermal import lookup
sys.path.insert(0, f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot"); import dep_common as C
g = np.load(f"{GRL}/data/thermal_grid.npz")
col = pd.read_csv("../g6_remaining/color/out_v3/crater_color.csv")[["id", "pc1", "pc2"]].rename(columns={"id": "host_uid"}); col["host_uid"] = col.host_uid.astype(str)
def xyz(la, lo): la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
prov = lambda la, lo: (xyz(la, lo) @ xyz(np.array([65.0]), np.array([125.0]))[0]) >= np.cos(500 / 2439.7)
# automatic depths
a = pd.read_csv("depth_all.csv"); fl = a["flags"].fillna("")
a = a[(a.depth_m > 0) & ~fl.str.contains("no_shadow|edge_sun|frame_edge|no_edr")].copy()
ha = pd.read_csv("../g2_insolation/inputs/v3/host_assignments.csv")[["hollow_id", "host_uid"]]; ha["host_uid"] = ha.host_uid.astype("Int64").astype(str)
a = a.merge(ha, on="hollow_id", how="left").merge(col, on="host_uid", how="left")
a["lnh"] = np.log(a.depth_m); a["lnD"] = np.log(a.D_can_m / 1000); a["lnpx"] = np.log(a.px); a["ltan"] = np.log(np.tan(np.radians(a.inc)))
a["T"] = lookup(g["la"], g["lo"], g["G"], a.lat.values, a.lon.values); a["x"] = 1 / a["T"]; a["pole"] = (a.lat >= 60).astype(int); a["abslat"] = a.lat.abs()
a["prov"] = prov(a.lat.values, a.lon.values); a["img"] = a.num.astype(str)
# Blewett depths
b = pd.read_csv(f"{GRL}/data/test1_rows_inc.csv"); b["x"] = 1 / b["T"]; b["ltan"] = np.log(np.tan(np.radians(b.inc))); b["pole"] = (b.abslat >= 60).astype(int)
b["img"] = b.filename.astype(str); b["prov"] = prov(b.latitude.values, b.longitude.values % 360)
hb = C.assign_hosts(b.latitude.values, b.longitude.values % 360); b["host_uid"] = hb[0] if isinstance(hb, tuple) else hb
b["host_uid"] = b.host_uid.astype(str); b = b.merge(col, on="host_uid", how="left")
for d in (a, b):
    for k in ("pc1", "pc2"): d[k + "z"] = (d[k] - d[k].mean()) / d[k].std()
def fit(d, f, keys):
    d = d.dropna(subset=[c for c in ["pc1z", "pc2z"] if c in f] + ["lnh"]); m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.img)[0]})
    o = {"n": int(m.nobs), "n_img": int(d.img.nunique())}
    for k in keys:
        if k == "x": o["coef_invT_K"] = [round(m.params[k], 1), round(m.bse[k], 1), round(float(m.pvalues[k]), 4)]
        else: o[k] = [round(float(np.exp(m.params[k])), 3), round(float(np.exp(m.params[k] - 1.96 * m.bse[k])), 3), round(float(np.exp(m.params[k] + 1.96 * m.bse[k])), 3)]
    if "lnD" in m.params: o["lnD_slope"] = [round(m.params["lnD"], 3), round(m.bse["lnD"], 3)]
    return o
R = {}
base = "lnh ~ {} + lnpx + ltan"
specs = [("T", "x", ["x"]), ("pole", "pole", ["pole"]), ("color", "pc1z + pc2z + abslat", ["pc1z", "pc2z"])]
for nm, d in [("blewett", b), ("auto", a), ("auto_no_province", a[~a.prov]), ("blewett_no_province", b[~b.prov])]:
    for sn, rhs, keys in specs:
        R[f"{nm}|{sn}"] = fit(d, base.format(rhs), keys)
        if nm.startswith("auto"): R[f"{nm}|{sn}+lnD"] = fit(d, base.format(rhs) + " + lnD", keys)
for k, v in R.items(): print(f"{k:32s}", v)
# context: sizes and color by pole in the automatic sample
R["auto_context"] = dict(n=len(a), n_pole=int(a.pole.sum()), prov_share_of_pole=round(float(a[a.pole == 1].prov.mean()), 3),
    medD_pole_m=round(float(a[a.pole == 1].D_can_m.median()), 0), medD_eq_m=round(float(a[a.pole == 0].D_can_m.median()), 0),
    rho_pc1_lnD=round(float(a[["pc1", "lnD"]].corr("spearman").iloc[0, 1]), 3), rho_abslat_pc1=round(float(a[["abslat", "pc1"]].corr("spearman").iloc[0, 1]), 3),
    blewett_prov_share_of_pole=round(float(b[b.pole == 1].prov.mean()), 3), blewett_n_pole=int(b.pole.sum()))
print(R["auto_context"])
json.dump(R, open("depth_regress.json", "w"), indent=1)
