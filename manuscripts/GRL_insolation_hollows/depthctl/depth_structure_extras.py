#!/usr/bin/env python3
"""Additions to depthctl/depth_structure_record.json: pair correlations of corrected depth and the red-province depth test.
Recreated for the archive (29 Oct 2026). These two blocks were first computed inline and the code was not saved; this script
reproduces the archived values exactly. Run after depth_structure.py, which rewrites the record without them.
Corrected depth is the residual of log10 depth on log10 pixel scale and log10 tan(incidence), as in depth_structure.py."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
# Reuse the depth table, corrected depth (d.res), and host-crater assignment (d.host) of depth_structure.py.
exec(open(os.path.join(HERE, "depth_structure.py")).read().replace("__file__", repr(os.path.join(HERE, "depth_structure.py"))).split('R = {"n_hollows"')[0])
REC = os.path.join(HERE, "depth_structure_record.json")
R = json.load(open(REC))
# ---- pairs of measurements within 60 km
X = xyz(d.latitude.values, d.longitude.values % 360)
pairs = np.array(sorted(cKDTree(X).query_pairs(2 * np.sin(60 / R_M / 2))))
i, j = pairs[:, 0], pairs[:, 1]
sep = np.arccos(np.clip((X[i] * X[j]).sum(1), -1, 1)) * R_M
same_img = d.filename.values[i] == d.filename.values[j]
same_crater = (d.host.values[i] == d.host.values[j]) & (d.host.values[i] >= 0)
r = d.res.values
def corr(m): return round(float(np.corrcoef(r[i][m], r[j][m])[0, 1]), 2)
P = {"note": "Pearson r of pixel- and incidence-corrected log10 depth between pairs of measurements, by separation and by image/crater membership",
     "diff_image_lt0.2km": corr(~same_img & (sep < 0.2)), "same_image_lt0.2km": corr(same_img & (sep < 0.2))}
for a, b in [(0, 0.5), (0.5, 2), (2, 5), (5, 10), (10, 20), (20, 60)]:
    m = (sep >= a) & (sep < b); row = {}
    for k, mm in (("diff_img_same_crater", m & ~same_img & same_crater), ("diff_img_diff_crater", m & ~same_img & ~same_crater), ("same_img", m & same_img)):
        if mm.sum() >= 20: row[k] = corr(mm)          # bins with fewer than 20 pairs are left out
    P[f"{a:g}-{b:g}km"] = row
dd = r[i][~same_img & (sep < 0.2)] - r[j][~same_img & (sep < 0.2)]
P["per_measurement_sd_dex_from_repeat_pairs"] = round(float(dd.std(ddof=1) / np.sqrt(2)), 3)
R["pair_correlations_corrected_depth"] = P
# ---- red province (90-150E, 40-70N) against the rest of 40-70N
lon360 = d.longitude % 360
s = d[(d.latitude >= 40) & (d.latitude <= 70)].copy(); s["inprov"] = ((lon360 >= 90) & (lon360 <= 150)).astype(int)
m = smf.ols("ld ~ inprov + lp + ltan", data=s).fit(cov_type="cluster", cov_kwds={"groups": s.filename.astype("category").cat.codes})
b, e = m.params["inprov"], m.bse["inprov"]
R["red_province_90_150E_40_70N_vs_rest_40_70N"] = dict(depth_ratio=round(float(10 ** b), 3), ci95=[round(float(10 ** (b - 1.96 * e)), 3), round(float(10 ** (b + 1.96 * e)), 3)],
    p=round(float(m.pvalues["inprov"]), 2), n_in=int(s.inprov.sum()), n_out=int((1 - s.inprov).sum()),
    note="province is S/Si-, Mg/Si-, Ca/Si-depleted and Al/Si-enriched (earlier XRS analysis); pixel scale and incidence controlled, image-clustered SE")
json.dump(R, open(REC, "w"), indent=1)
print(json.dumps({k: R[k] for k in ("pair_correlations_corrected_depth", "red_province_90_150E_40_70N_vs_rest_40_70N")}, indent=1))
