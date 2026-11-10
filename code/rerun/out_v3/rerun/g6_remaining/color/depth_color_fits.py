#!/usr/bin/env python3
"""Depth of the Blewett et al. (2016) hollows against the color of their surroundings, on the screened-mask colors (9 Nov 2026).
Same specification as part (2) of g4_completeness/g1/summary_model/step9_color.py (OLS of ln depth on standardized PC1 and PC2
of the surroundings within 10 km, ln pixel scale and ln tan(incidence), with |latitude|, the poleward indicator, or no latitude
control; standard errors clustered by image; 10th-90th percentile range in SD units), but with the hollow colors of
out_<cat>/hollow_color.csv (hollows of the screened v3 catalog masked) in place of analysis_workshop/color_mdr/hollow_color.csv
(version-2 mask). Usage: depth_color_fits.py [cat=v3]; writes out_<cat>/depth_color_fits.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
cat = sys.argv[1] if len(sys.argv) > 1 else "v3"; HERE = os.path.dirname(os.path.abspath(__file__))
t = pd.read_csv(os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")); t["row"] = np.arange(len(t))
h = pd.read_csv(f"{HERE}/out_{cat}/hollow_color.csv")
t = t.merge(h[["row", "set_pc1", "set_pc2"]], on="row", how="left").dropna(subset=["set_pc1", "set_pc2"])
t["ltan"] = np.log(np.tan(np.radians(t.inc))); t["pole"] = (t.abslat >= 60).astype(int); g = pd.factorize(t.filename)[0]
for k in ["set_pc1", "set_pc2"]: t[k + "_z"] = (t[k] - t[k].mean()) / t[k].std()
D = {}
for name, f in [("pole", "lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan + pole"), ("abslat", "lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan + abslat"),
                ("none", "lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan")]:
    m = smf.ols(f, t).fit(cov_type="cluster", cov_kwds={"groups": g})
    D[name] = {k: [float(np.exp(m.params[k])), float(np.exp(m.params[k] - 1.96 * m.bse[k])), float(np.exp(m.params[k] + 1.96 * m.bse[k])), float(m.pvalues[k])] for k in ["set_pc1_z", "set_pc2_z"]}
    D[name]["p10_p90_sd"] = {k: float((t[k].quantile(.9) - t[k].quantile(.1)) / t[k].std()) for k in ["set_pc1", "set_pc2"]}
    s = D[name]["p10_p90_sd"]["set_pc1"]; D[name]["pc1_factor_over_p10_p90"] = [float(x ** s) for x in D[name]["set_pc1_z"][:3]]
R = {"hollow_color": f"out_{cat}/hollow_color.csv", "n": int(len(t)), "n_img": int(len(np.unique(g))), "fits": D}
json.dump(R, open(f"{HERE}/out_{cat}/depth_color_fits.json", "w"), indent=1)
print(json.dumps(R, indent=1))
