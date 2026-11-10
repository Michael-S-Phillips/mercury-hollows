#!/usr/bin/env python3
"""Depth vs LRM (3 Nov 2026), for v6 Section 3.6 / S9. LRM: data/Concentrated_LRM_560_params_BW.tif (no georeference; read as
-180..180 E, 90N at top; alignment checked visually against the MDIS basemap, summary_model/_scratch/lrm_alignment_check.png).
lrm_local = fraction of LRM pixels (value > 0) within 10 km of the hollow. ln depth regressed on lrm_local (and on an LRM
indicator, lrm_local >= 0.5), ln pixel scale, ln tan(incidence), and poleward-of-60 flag; SEs clustered by image; also
host-level LRM fraction (crater_features.csv) for hollows assigned to hosts. Writes depth_lrm_record.json."""
import os, json, numpy as np, pandas as pd, rasterio, statsmodels.formula.api as smf, warnings; warnings.filterwarnings("ignore")
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")
d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["pole"] = (d.abslat >= 60).astype(int); d["img"] = d.filename.astype(str)
L = rasterio.open(f"{ROOT}/data/Concentrated_LRM_560_params_BW.tif").read(1) > 0; ny, nx = L.shape; ppd = nx / 360
R = 2439.4; out = []
for la, lo in zip(d.latitude.values, ((d.longitude.values + 180) % 360) - 180):
    r = 10 / (R * np.pi / 180); cy, cx = (90 - la) * ppd, (lo + 180) * ppd; ry = r * ppd; rx = ry / max(np.cos(np.radians(la)), 0.05)
    ys = np.arange(max(int(cy - ry) - 1, 0), min(int(cy + ry) + 2, ny)); xs = np.arange(int(cx - rx) - 1, int(cx + rx) + 2)
    yy, xx = np.meshgrid(ys, xs, indexing="ij"); m = ((yy + .5 - cy) / ry) ** 2 + ((xx + .5 - cx) / rx) ** 2 <= 1
    out.append(L[yy[m], xx[m] % nx].mean() if m.any() else L[int(cy), int(cx) % nx])
d["lrm_local"] = out; d["lrm_bin"] = (d.lrm_local >= 0.5).astype(int)
Rr = {"n": int(len(d)), "n_img": int(d.img.nunique()), "lrm_local_quantiles": [float(x) for x in d.lrm_local.quantile([.1, .25, .5, .75, .9])],
      "share_lrm_bin": float(d.lrm_bin.mean())}
g = pd.factorize(d.img)[0]
for name, f in [("continuous", "lnh ~ lrm_local + lnpx + ltan + pole"), ("binary", "lnh ~ lrm_bin + lnpx + ltan + pole")]:
    m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": g}); k = f.split("~")[1].split("+")[0].strip()
    b, se = m.params[k], m.bse[k]
    Rr[name] = dict(ratio_full_range=float(np.exp(b)), lo=float(np.exp(b - 1.96 * se)), hi=float(np.exp(b + 1.96 * se)), p=float(m.pvalues[k]))
# bootstrap by image for the binary effect
rng = np.random.default_rng(3); imgs = d.img.unique(); bs = []
for _ in range(1000):
    pick = rng.choice(imgs, len(imgs)); dd = pd.concat([d[d.img == i] for i in pick[:0]]) if False else None
    cnt = pd.Series(pick).value_counts(); w = d.img.map(cnt).fillna(0).values
    m = smf.wls("lnh ~ lrm_bin + lnpx + ltan + pole", d, weights=w).fit(); bs.append(m.params["lrm_bin"])
Rr["binary"]["boot_ci"] = [float(np.exp(np.percentile(bs, 2.5))), float(np.exp(np.percentile(bs, 97.5)))]
json.dump(Rr, open(f"{HERE}/depth_lrm_record.json", "w"), indent=1); print(json.dumps(Rr, indent=1))
