#!/usr/bin/env python3
"""Depth vs setting color and LRM together (3 Nov 2026), for v6 Section 3.6. Same rows and LRM fraction as
composition_v6/depth_lrm.py; setting PC1/PC2 from hollow_color.csv (within 10 km, hollows masked). ln depth on
z(set_pc1), z(set_pc2), lrm_local, ln pixel scale, ln tan(incidence), and a poleward flag or |lat|; SEs clustered by image.
Also the correlation of setting PC1 with lrm_local. Writes depth_color_lrm_record.json."""
import os, json, numpy as np, pandas as pd, rasterio, statsmodels.formula.api as smf, warnings; warnings.filterwarnings("ignore")
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv"); d["row"] = np.arange(len(d))
d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["pole"] = (d.abslat >= 60).astype(int); d["img"] = d.filename.astype(str)
L = rasterio.open(f"{ROOT}/data/Concentrated_LRM_560_params_BW.tif").read(1) > 0; ny, nx = L.shape; ppd = nx / 360
R = 2439.4; out = []
for la, lo in zip(d.latitude.values, ((d.longitude.values + 180) % 360) - 180):
    r = 10 / (R * np.pi / 180); cy, cx = (90 - la) * ppd, (lo + 180) * ppd; ry = r * ppd; rx = ry / max(np.cos(np.radians(la)), 0.05)
    ys = np.arange(max(int(cy - ry) - 1, 0), min(int(cy + ry) + 2, ny)); xs = np.arange(int(cx - rx) - 1, int(cx + rx) + 2)
    yy, xx = np.meshgrid(ys, xs, indexing="ij"); m = ((yy + .5 - cy) / ry) ** 2 + ((xx + .5 - cx) / rx) ** 2 <= 1
    out.append(L[yy[m], xx[m] % nx].mean() if m.any() else L[int(cy), int(cx) % nx])
d["lrm_local"] = out
h = pd.read_csv(f"{HERE}/hollow_color.csv")
d = d.merge(h[["row", "set_pc1", "set_pc2"]], on="row", how="left").dropna(subset=["set_pc1", "set_pc2"]).reset_index(drop=True)
for k in ["set_pc1", "set_pc2"]: d[k + "_z"] = (d[k] - d[k].mean()) / d[k].std()
g = pd.factorize(d.img)[0]; Rr = {"n": int(len(d)), "rho_pc1_lrm": float(d[["set_pc1", "lrm_local"]].corr("spearman").iloc[0, 1])}
for name, f in [("pole", "lnh ~ set_pc1_z + set_pc2_z + lrm_local + lnpx + ltan + pole"), ("abslat", "lnh ~ set_pc1_z + set_pc2_z + lrm_local + lnpx + ltan + abslat")]:
    m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": g})
    Rr[name] = {k: [float(np.exp(m.params[k])), float(np.exp(m.params[k] - 1.96 * m.bse[k])), float(np.exp(m.params[k] + 1.96 * m.bse[k])), float(m.pvalues[k])] for k in ["set_pc1_z", "set_pc2_z", "lrm_local"]}
json.dump(Rr, open(f"{HERE}/depth_color_lrm_record.json", "w"), indent=1); print(json.dumps(Rr, indent=1))
