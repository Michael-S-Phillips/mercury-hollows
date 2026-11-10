#!/usr/bin/env python3
"""Workshop test 2 (2 Nov 2026): does hollow depth change with hollow size (version E: small, young hollows shallower until
they reach a stopping depth; version C: no trend), and are shadow-length depths capped by hollow width (shadow saturation)?

Inputs: GRL_insolation_hollows/data/test1_rows_inc.csv (2,518 depths, Blewett et al. 2016);
GIS/vector/bickel2024_hollows_corrected.gpkg, layer hollows_points_corrected (19,110 detections, ISIS camera-model positions);
runs/bickel_hollow_areas.csv (same row order: release polygon area, equivalent diameter, passes_QC).
Match: each depth to corrected detections (confidence >= 0.75, incidence < 85) within R km (great circle, R_M = 2439.4 km);
size W = median equivalent diameter of the matched detections; R = 0.5 km (sensitivity 0.3, 1.0).
Regression: ln h ~ ln W + ln px + ln tan(inc) [+ pole], SE clustered by image.
Saturation: shadow length L = h tan(inc); a shadow cannot run longer than the hollow is wide, so L/W near 1 flags capped depths.
Writes test2_depth_size_record.json.
"""
import os, json, sqlite3, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.spatial import cKDTree
ROOT = os.path.expanduser("~/mnt/hollows-emissivity")
d = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")
d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["img"] = d.filename.astype(str); d["pole"] = (d.abslat >= 60).astype(int)
c = sqlite3.connect(f"{ROOT}/GIS/vector/bickel2024_hollows_corrected.gpkg")
b = pd.read_sql("select fid, mdis_id, confidence, incidence, px_scale_m, lat, lon180, lat_orig, lon180_orig from hollows_points_corrected order by fid", c)
a = pd.read_csv(f"{ROOT}/runs/bickel_hollow_areas.csv")  # one row fewer than the release; join on release position + pixel scale
key = lambda la, lo, px: (np.round(la, 4).astype(str) + "_" + np.round(lo % 360, 4).astype(str) + "_" + np.round(px, 2).astype(str))
b["k"] = key(b.lat_orig, b.lon180_orig, b.px_scale_m); a["k"] = key(a.lat, a.lon_E, a.pixel_scale_m)
a = a.drop_duplicates("k")
b = b.merge(a[["k", "equiv_diam_m", "area_km2", "passes_QC"]], on="k", how="left")
R0 = {"n_release": int(len(b)), "n_without_area": int(b.equiv_diam_m.isna().sum())}
b = b.dropna(subset=["equiv_diam_m"]); b["W"] = b.equiv_diam_m
b = b[(b.confidence >= 0.75) & (b.incidence < 85)].reset_index(drop=True)
RM = 2439.4
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
tree = cKDTree(xyz(b.lat.values, b.lon180.values)); P = xyz(d.latitude.values, ((d.longitude.values + 180) % 360) - 180)
R = {"n_depths": int(len(d)), "n_detections_used": int(len(b)), **R0}
out = {}
for rk in [0.3, 0.5, 1.0]:
    W = np.full(len(d), np.nan); nm = np.zeros(len(d), int); Wpx = np.full(len(d), np.nan)
    for i, cl in enumerate(tree.query_ball_point(P, r=rk / RM)):
        if cl:
            W[i] = np.median(b.W.values[cl]); nm[i] = len(cl); Wpx[i] = np.median(b.px_scale_m.values[cl])
    s = d.assign(W=W, nmatch=nm, Wpx=Wpx).dropna(subset=["W"]).copy()
    s["lnW"] = np.log(s.W); s["L"] = s.depth * np.tan(np.radians(s.inc)); s["LW"] = s.L / s.W
    g = pd.factorize(s.img)[0]
    rec = dict(n_matched=int(len(s)), frac_matched=float(len(s) / len(d)), n_pole=int(s.pole.sum()),
               W_quantiles_m=[float(x) for x in s.W.quantile([.1, .25, .5, .75, .9])])
    for f in ["lnh ~ lnW + lnpx + ltan", "lnh ~ lnW + lnpx + ltan + pole", "lnh ~ pole + lnpx + ltan"]:
        m = smf.ols(f, s).fit(cov_type="cluster", cov_kwds={"groups": g})
        rec[f] = {k: dict(coef=float(m.params[k]), lo=float(m.params[k] - 1.96 * m.bse[k]), hi=float(m.params[k] + 1.96 * m.bse[k]))
                  for k in m.params.index if k in ("lnW", "pole")}
    rec["LW_quantiles"] = {grp: [float(x) for x in s[s.pole == p].LW.quantile([.5, .75, .9, .95])] for grp, p in [("eq", 0), ("pole", 1)]}
    rec["frac_LW_gt_0.5"] = {grp: float((s[s.pole == p].LW > 0.5).mean()) for grp, p in [("eq", 0), ("pole", 1)]}
    rec["W_median_m"] = {grp: float(s[s.pole == p].W.median()) for grp, p in [("eq", 0), ("pole", 1)]}
    # depth by size tercile, corrected to median px and inc
    mc = smf.ols("lnh ~ lnpx + ltan", s).fit()
    s["hc"] = np.exp(s.lnh - mc.params["lnpx"] * (s.lnpx - s.lnpx.median()) - mc.params["ltan"] * (s.ltan - s.ltan.median()))
    s["Wbin"] = pd.qcut(s.W, 4)
    rec["corrected_depth_by_W_quartile"] = {str(k): dict(n=int(len(v)), W_med=float(v.W.median()), h_med=float(v.hc.median()))
                                            for k, v in s.groupby("Wbin", observed=True)}
    out[f"R{rk}"] = rec
    if rk == 0.5: s.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test2_matched_R0.5.csv"), index=False)
R["by_radius"] = out
json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test2_depth_size_record.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
