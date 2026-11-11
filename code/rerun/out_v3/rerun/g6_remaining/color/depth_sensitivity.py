#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script (no script was kept for these three S9 numbers); same rows, variables and fit as step9_color.py part (2).
"""S9 sensitivity fits of ln depth on the setting PC1 (z) and PC2 (z) with ln pixel scale and ln tan(incidence), SEs clustered by
image: (a) |lat| control, measurements within 60 deg of the equator; (b) |lat| and ln geometric-mean reflectance of the setting;
(c) 10 deg latitude bands (categorical) instead of |lat|. Usage: depth_sensitivity.py <cat>; writes out_<cat>/depth_sensitivity.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
cat = sys.argv[1]; HERE = os.path.dirname(os.path.abspath(__file__))
t = pd.read_csv(os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/insolation_thermal/data/test1_rows_inc.csv")); t["row"] = np.arange(len(t))
h = pd.read_csv(f"{HERE}/out_{cat}/hollow_color.csv")
t = t.merge(h[["row", "set_pc1", "set_pc2", "set_albedo_gm"]], on="row", how="left").dropna(subset=["set_pc1", "set_pc2"]).reset_index(drop=True)
t["ltan"] = np.log(np.tan(np.radians(t.inc))); t["lnalb"] = np.log(t.set_albedo_gm); t["band"] = (t.latitude // 10).astype(int).astype(str)
for k in ["set_pc1", "set_pc2"]: t[k + "_z"] = (t[k] - t[k].mean()) / t[k].std()
def fit(f, d):
    g = pd.factorize(d.filename)[0]; m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": g}); b, s = m.params["set_pc1_z"], m.bse["set_pc1_z"]
    return dict(n=int(len(d)), factor=float(np.exp(b)), lo=float(np.exp(b - 1.96 * s)), hi=float(np.exp(b + 1.96 * s)), p=float(m.pvalues["set_pc1_z"]))
R = {"within_60": fit("lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan + abslat", t[t.abslat <= 60]),
     "with_ln_albedo": fit("lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan + abslat + lnalb", t),
     "lat_bands_10deg": fit("lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan + C(band)", t)}
json.dump(R, open(f"{HERE}/out_{cat}/depth_sensitivity.json", "w"), indent=1)
print(cat, {k: [round(v[x], 3) for x in ("factor", "lo", "hi")] + [v["n"]] for k, v in R.items()})
