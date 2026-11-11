#!/usr/bin/env python3
"""g23 step 12 (9 Nov 2026): the depth results of the paper (Sections 3.6-3.7, S9-S11, S15) recomputed on our automatic depths
(depth_all.csv: 1,558 clean hollows, finest detection at <= 20 m/px, incidence 58-81 deg, emission <= 50 deg; validated against
Michael's hand depths), each without and with ln hollow diameter. OLS of ln depth with ln pixel scale and ln tan(incidence),
standard errors clustered by image, as in the paper. Host = the v3 host assignment (smallest crater >= 5 km whose 1.1 R contain the
hollow); host attributes from crater_features_xrs.csv (v3r) and kinczyk_matched.csv; LRM within 10 km from the Klima et al. map as in
depth_color_lrm.py; T from the paper's thermal grid. Writes depth_paper.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, rasterio, warnings; warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
HE = os.path.expanduser("~/mnt/hollows-emissivity"); INSOL = f"{HE}/manuscripts/insolation_thermal"; RR = f"{HE}/code/rerun/out_v3/rerun"
sys.path.insert(0, f"{INSOL}/scripts"); from thermal import lookup
g = np.load(f"{INSOL}/data/thermal_grid.npz")
a = pd.read_csv("depth_all.csv"); fl = a["flags"].fillna("")
a = a[(a.depth_m > 0) & ~fl.str.contains("no_shadow|edge_sun|frame_edge|no_edr")].copy()
ha = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/host_assignments.csv")[["hollow_id", "host_uid", "rho"]]; ha["host_uid"] = ha.host_uid.astype("Int64").astype(str)
cf = pd.read_csv(f"{RR}/g4_completeness/g1/summary_model/v3r/crater_features_xrs.csv", dtype={"id": str})
cf = cf[["id", "D", "nested", "age", "vent", "parent_vent", "lrm", "s_si", "mg_si", "ca_si", "al_si"]].rename(columns={"id": "host_uid", "D": "hostD"})
kz = pd.read_csv(f"{RR}/g4_completeness/g1/degradation/v3r/kinczyk_matched.csv", dtype={"id": str})[["id", "kclass"]].rename(columns={"id": "host_uid"})
kz = kz.dropna(subset=["host_uid"]).drop_duplicates("host_uid")   # 202 unmatched Kinczyk craters have no id
col = pd.read_csv(f"{RR}/g6_remaining/color/out_v3/crater_color.csv", dtype={"id": str}, usecols=["id", "ok", "pc1", "pc2"]).rename(columns={"id": "host_uid"})
a = a.merge(ha, on="hollow_id", how="left", validate="1:1").merge(cf, on="host_uid", how="left", validate="m:1").merge(kz, on="host_uid", how="left", validate="m:1").merge(col, on="host_uid", how="left", validate="m:1")
a["lnh"] = np.log(a.depth_m); a["lnD"] = np.log(a.D_can_m / 1000); a["lnpx"] = np.log(a.px); a["ltan"] = np.log(np.tan(np.radians(a.inc))); a["img"] = a.num.astype(str)
a["T"] = lookup(g["la"], g["lo"], g["G"], a.lat.values, a.lon.values); a["x"] = 1 / a["T"]; a["pole"] = (a.lat >= 60).astype(int); a["abslat"] = a.lat.abs()
a["c2l"] = np.cos(np.radians(2 * a.lon)); a["loghostD"] = np.log10(a.hostD)
L = rasterio.open(f"{HE}/data/Concentrated_LRM_560_params_BW.tif").read(1) > 0; ny, nx = L.shape; ppd = nx / 360; R = 2439.4; lr = []
for la, lo in zip(a.lat.values, ((a.lon.values + 180) % 360) - 180):
    r = 10 / (R * np.pi / 180); cy, cx = (90 - la) * ppd, (lo + 180) * ppd; ry = r * ppd; rx = ry / max(np.cos(np.radians(la)), 0.05)
    ys = np.arange(max(int(cy - ry) - 1, 0), min(int(cy + ry) + 2, ny)); xs = np.arange(int(cx - rx) - 1, int(cx + rx) + 2)
    yy, xx = np.meshgrid(ys, xs, indexing="ij"); m = ((yy + .5 - cy) / ry) ** 2 + ((xx + .5 - cx) / rx) ** 2 <= 1
    lr.append(L[yy[m], xx[m] % nx].mean() if m.any() else L[int(cy), int(cx) % nx])
a["lrm10"] = lr; a["lrm_half"] = (a.lrm10 >= 0.5).astype(int)
okc = a.ok.fillna(False).astype(bool)
for k in ["pc1", "pc2"]: a[k + "z"] = np.where(okc, (a[k] - a.loc[okc, k].mean()) / a.loc[okc, k].std(), np.nan)
for k in ["s_si", "mg_si", "ca_si", "al_si"]: a[k + "_z"] = a[k]
a["fresh"] = (a.age == "fresh").astype(int); a["degraded"] = (a.age == "degraded").astype(int)
a["k3"] = (a.kclass == 3).astype(int); a["k45"] = a.kclass.isin([4, 5]).astype(int)
a.to_csv("depth_paper_rows.csv", index=False)
def fit(d, rhs, keys, size=False):
    f = "lnh ~ " + rhs + " + lnpx + ltan" + (" + lnD" if size else ""); d = d.dropna(subset=[c for c in keys if c in d.columns])
    m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.img)[0]}); o = {"n": int(m.nobs), "n_img": int(d.img.nunique())}
    for k in keys:
        o[k] = [round(float(m.params[k]), 4), round(float(m.bse[k]), 4), round(float(m.pvalues[k]), 4)]
        o[k + "_ratio"] = [round(float(np.exp(m.params[k])), 3), round(float(np.exp(m.params[k] - 1.96 * m.bse[k])), 3), round(float(np.exp(m.params[k] + 1.96 * m.bse[k])), 3)]
    if size: o["lnD"] = [round(float(m.params["lnD"]), 3), round(float(m.bse["lnD"]), 3)]
    return o
R = {"n": int(len(a)), "n_img": int(a.img.nunique()), "n_hosts": int(a.host_uid.nunique()), "T_range": [round(a["T"].min(), 1), round(a["T"].max(), 1)],
     "median_depth_m": round(float(a.depth_m.median()), 1), "median_D_m": round(float(a.D_can_m.median()), 0)}
eq = a[a.abslat < 60]
tests = {"T": (a, "x", ["x"]), "pole": (a, "pole", ["pole"]), "hot_lon_lt60": (eq, "c2l + abslat", ["c2l"]), "T_lt60": (eq, "x", ["x"]),
         "hostD": (a, "loghostD + abslat", ["loghostD"]), "nested": (a.dropna(subset=["nested"]), "nested + abslat", ["nested"]),
         "lrm_half": (a, "lrm_half + abslat", ["lrm_half"]), "color": (a, "pc1z + pc2z + abslat", ["pc1z", "pc2z"]),
         "herrick": (a[a.age.isin(["fresh", "standard", "degraded"])], "fresh + degraded + abslat", ["fresh", "degraded"]),
         "kinczyk": (a[a.kclass.notna()], "k3 + k45 + abslat", ["k3", "k45"]), "vent": (a.dropna(subset=["vent"]), "vent + abslat", ["vent"])}
for nm, (d, rhs, keys) in tests.items():
    R[nm] = fit(d, rhs, keys); R[nm + "+lnD"] = fit(d, rhs, keys, size=True)
# temperature: change across the T range and apparent enthalpy
for k in ["T", "T+lnD"]:
    b, se = R[k]["x"][0], R[k]["x"][1]; t0, t1 = 219.0, 324.0; dx = 1 / t0 - 1 / t1
    R[k]["depth_ratio_324_over_219"] = [round(float(np.exp(-b * dx)), 3), round(float(np.exp(-(b + 1.96 * se) * dx)), 3), round(float(np.exp(-(b - 1.96 * se) * dx)), 3)]
    R[k]["apparent_dH_kJ"] = [round(-2 * 8.314 * b / 1000, 2), round(-2 * 8.314 * (b + 1.96 * se) / 1000, 2), round(-2 * 8.314 * (b - 1.96 * se) / 1000, 2)]
# XRS: change across the 10th-90th percentile of each ratio (host level), one ratio at a time
R["xrs"] = {}
for k in ["s_si", "mg_si", "ca_si", "al_si"]:
    d = a.dropna(subset=[k]).copy(); lo, hi = d[k].quantile([.1, .9]); d["z"] = (d[k] - lo) / (hi - lo)
    for sz in (False, True):
        o = fit(d, "z + abslat", ["z"], size=sz); R["xrs"][k + ("+lnD" if sz else "")] = dict(n=o["n"], ratio_10_90=o["z_ratio"])
# depth differs between hosts; spatial coherence of host means (>= 3 hollows), depth corrected for px, tan i and ln D
m = smf.ols("lnh ~ lnpx + ltan + lnD", a).fit(); a["res"] = m.resid
hm = a.groupby("host_uid").agg(n=("res", "size"), res=("res", "mean"), lat=("lat", "first"), lon=("lon", "first")); hm = hm[hm.n >= 3]
from scipy import stats
grp = [a[a.host_uid == h].res.values for h in hm.index]; F, p = stats.f_oneway(*grp)
def xyz(la, lo): la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
P = xyz(hm.lat.values, hm.lon.values); D = 2439.7 * np.arccos(np.clip(P @ P.T, -1, 1)); iu = np.triu_indices(len(hm), 1)
near = D[iu] < 600; x1, x2 = hm.res.values[iu[0]], hm.res.values[iu[1]]
R["between_hosts"] = dict(hosts=int(len(hm)), anova_F=round(float(F), 2), anova_p=float(p), r_pairs_lt600=round(float(np.corrcoef(x1[near], x2[near])[0, 1]), 3), n_pairs_lt600=int(near.sum()))
# Fresh against degraded hosts only, and hollows, hosts and medians by Kinczyk class (Section S11). Added 9 Nov 2026: these keys were
# first computed outside this script; the code below reproduces the stored values exactly.
fd = a[a.age.isin(["fresh", "degraded"])]
for sz in (False, True):
    o = fit(fd, "fresh + abslat", ["fresh"], size=sz)
    R["fresh_vs_degraded" + ("+lnD" if sz else "")] = dict(n=o["n"], hosts=int(fd.host_uid.nunique()), ratio=o["fresh_ratio"])
kc = a[a.kclass.notna()].copy(); kc["kc"] = np.where(kc.kclass <= 2, "1-2", np.where(kc.kclass == 3, "3", "4-5")); kg = kc.groupby("kc")
R["kinczyk_n"] = {k: int(v) for k, v in kg.size().items()}
R["kinczyk_hosts"] = {k: int(v) for k, v in kg.host_uid.nunique().items()}
R["kinczyk_median_D_m"] = {k: float(np.round(v, 0)) for k, v in kg.D_can_m.median().items()}   # numpy rounding, as stored
R["kinczyk_median_depth_m"] = {k: float(np.round(v, 1)) for k, v in kg.depth_m.median().items()}
json.dump(R, open("depth_paper.json", "w"), indent=1, default=str)
for k, v in R.items(): print(k, v)
