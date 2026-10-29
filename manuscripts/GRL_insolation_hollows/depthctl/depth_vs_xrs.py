#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Does hollow depth follow regional composition? Samples MESSENGER XRS ratio maps (MESS-H-XRS-3-RDR-MAPS-V1.0, messxrs_3001,
2015-04-24 release; 4 pix/deg simple cylindrical, rows 90N->90S, columns 180W->180E, 8-bit with label scaling) at each depth
measurement and regresses corrected log depth on each ratio. Writes depthctl/depth_vs_xrs_record.json."""
import json, re, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
H = _REPO + "/manuscripts/GRL_insolation_hollows/depthctl"
def load(kind, el):
    lbl = open(f"{H}/../data/xrs_maps/xrs_{kind}_{el}_20150424.lbl").read()
    sf = float(re.search(r"SCALING_FACTOR\s*=\s*([0-9.eE+-]+)", lbl).group(1))
    a = np.fromfile(f"{H}/../data/xrs_maps/xrs_{kind}_{el}_20150424.img", dtype=np.uint8).reshape(720, 1440).astype(float)
    a[a == 0] = np.nan; return a * sf
els = ["s__si", "mg_si", "ca_si", "al_si", "fe_si"]
M = {e: load("map", e) for e in els}; U = {e: load("unc", e) for e in els}; F = {e: load("res", e) for e in els}
d = pd.read_csv(_REPO + "/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")
d["ltan"] = np.log10(np.tan(np.radians(d.inc))); d["lp"] = np.log10(d.px); d["ld"] = np.log10(d.depth)
d["img"] = d.filename.astype("category").cat.codes
lon180 = ((d.longitude + 180) % 360) - 180
r = np.clip(((90 - d.latitude) * 4).astype(int), 0, 719).values; c = np.clip(((lon180 + 180) * 4).astype(int), 0, 1439).values
for e in els:
    d[e] = M[e][r, c]; d[e + "_unc"] = U[e][r, c]; d[e + "_fp"] = F[e][r, c]
R = {"n": int(len(d))}
for e in els:
    s = d[d[e].notna()].copy(); s["z"] = (s[e] - s[e].mean()) / s[e].std()
    m = smf.ols("ld ~ z + lp + ltan", data=s).fit(cov_type="cluster", cov_kwds={"groups": s.img})
    m2 = smf.ols("ld ~ z + lp + ltan + abslat", data=s).fit(cov_type="cluster", cov_kwds={"groups": s.img})
    lo, hi = s[e].quantile(.1), s[e].quantile(.9)
    slope_per_unit = m.params["z"] / s[e].std()
    R[e] = dict(n=int(len(s)), n_img=int(s.img.nunique()), mean=float(s[e].mean()), p10=float(lo), p90=float(hi),
                median_footprint_km=float(s[e + "_fp"].median()), median_unc=float(s[e + "_unc"].median()),
                coef_per_sd_log10=[float(m.params["z"]), float(m.bse["z"]), float(m.pvalues["z"])],
                depth_ratio_p90_vs_p10=[float(10 ** (slope_per_unit * (hi - lo))),
                                        float(10 ** ((m.params["z"] - 1.96 * m.bse["z"]) / s[e].std() * (hi - lo))),
                                        float(10 ** ((m.params["z"] + 1.96 * m.bse["z"]) / s[e].std() * (hi - lo)))],
                with_latitude_coef_per_sd=[float(m2.params["z"]), float(m2.bse["z"]), float(m2.pvalues["z"])])
# collapse prediction for sulfur: depth proportional to S content -> ratio p90/p10 of S/Si
e = "s__si"; R["collapse_prediction_if_phi_proportional_to_S"] = float(R[e]["p90"] / R[e]["p10"])
json.dump(R, open(f"{H}/depth_vs_xrs_record.json", "w"), indent=1)
for e in els:
    x = R[e]; print(f"{e:6s} n={x['n']} fp={x['median_footprint_km']:.0f} km range {x['p10']:.3f}-{x['p90']:.3f}  coef/sd {x['coef_per_sd_log10'][0]:+.3f}±{x['coef_per_sd_log10'][1]:.3f} p={x['coef_per_sd_log10'][2]:.3f}  depth ratio p90/p10 {x['depth_ratio_p90_vs_p10'][0]:.2f} ({x['depth_ratio_p90_vs_p10'][1]:.2f}-{x['depth_ratio_p90_vs_p10'][2]:.2f})  w/lat {x['with_latitude_coef_per_sd'][0]:+.3f}±{x['with_latitude_coef_per_sd'][1]:.3f}")
print("collapse prediction (S_p90/S_p10):", R["collapse_prediction_if_phi_proportional_to_S"])
