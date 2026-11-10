#!/usr/bin/env python3
"""(v2, 5 Nov 2026: all cataloged hollows from hollows_canonical.csv; v1 used host_assignments.csv, which holds only hollows inside host craters.)
Hollow fields and the higher-albedo smooth plains of Fassett et al. (2009, Fig. 9) inside b59 (5 Nov 2026). Mask of
bright, red plains in the MDIS MDR color within 1.2 b59 radii: 750-nm reflectance and 1000/480-nm ratio both above their
60th percentiles in the area, median-filtered (9 px) and blobs < 2,000 km^2 removed. For hollow fields (field_id; one point
per field at its mean position) and for random points in the same area, the fraction inside the mask, within 25 km of its
edge outside it, and farther away. Writes hrp_edge_test.json and fig_hrp_mask.png."""
import os, sys, json, numpy as np, pandas as pd
from scipy import ndimage as ndi
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
sys.path.insert(0, f"{HERE}/../color_mdr"); import mdr_georef as g
c_la, c_lo, Rb = 50.05, 119.26, 700.6 / 2
LAT0, LAT1, LON0, LON1, STEP = 38.0, 62.0, 98.0, 141.0, 0.04
la = np.arange(LAT1, LAT0, -STEP); lo = np.arange(LON0, LON1, STEP); LO, LA = np.meshgrid(lo, la)
G = g.Globe(); r750 = G.sample(LA.ravel(), LO.ravel(), band=4).reshape(LA.shape); r480 = G.sample(LA.ravel(), LO.ravel(), band=1).reshape(LA.shape); r1000 = G.sample(LA.ravel(), LO.ravel(), band=7).reshape(LA.shape)
def dist_km(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2); dl = np.radians(lo2 - lo1)
    return 2 * R_M * np.arcsin(np.sqrt(np.clip(np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2, 0, 1)))
rb = dist_km(LA, LO, c_la, c_lo) / Rb; area = rb < 1.2
ratio = r1000 / r480; ok = area & np.isfinite(r750) & np.isfinite(ratio)
m = ok & (r750 > np.nanpercentile(r750[ok], 60)) & (ratio > np.nanpercentile(ratio[ok], 60))
m = ndi.median_filter(m.astype(np.uint8), size=9).astype(bool) & ok
pix_km2 = (STEP * np.pi / 180 * R_M) ** 2 * np.cos(np.radians(LA))
lab, n = ndi.label(m); sizes = ndi.sum(pix_km2, lab, range(1, n + 1)); keep = np.isin(lab, 1 + np.where(sizes >= 2000)[0]); m = keep
# distance (km) of every pixel to the mask edge, outside the mask
dpix = ndi.distance_transform_edt(~m, sampling=(STEP * np.pi / 180 * R_M, STEP * np.pi / 180 * R_M * np.cos(np.radians(50))))
def classify(lat, lon):
    i = np.clip(((LAT1 - lat) / STEP).astype(int), 0, LA.shape[0] - 1); j = np.clip(((lon - LON0) / STEP).astype(int), 0, LA.shape[1] - 1)
    inside = m[i, j]; near = (~inside) & (dpix[i, j] <= 25); return inside, near, ok[i, j]
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)]
h["lon360"] = h.lon % 360; fld = h.groupby("field_id").agg(lat=("lat", "mean"), lon=("lon360", "mean"), n=("lat", "size")).reset_index()
fld = fld[dist_km(fld.lat.values, fld.lon.values, c_la, c_lo) / Rb < 1.2]
fi, fn, fok = classify(fld.lat.values, fld.lon.values); fld = fld[fok]; fi, fn = fi[fok], fn[fok]
w = pix_km2 * ok; A = w.sum()
R = {"mask_area_frac": round(float((w * m).sum() / A), 3), "near_area_frac": round(float((w * ((~m) & (dpix <= 25))).sum() / A), 3),
     "fields": int(len(fld)), "hollows": int(fld.n.sum()),
     "fields_inside": int(fi.sum()), "fields_near_edge_25km": int(fn.sum()), "fields_farther": int((~fi & ~fn).sum()),
     "hollows_inside": int(fld.n[fi].sum()), "hollows_near": int(fld.n[fn].sum()), "hollows_farther": int(fld.n[~fi & ~fn].sum())}
from scipy.stats import binomtest
R["binom_p_fields_near_vs_area"] = float(binomtest(int(fn.sum()), int(len(fld)), R["near_area_frac"]).pvalue)
R["binom_p_fields_inside_vs_area"] = float(binomtest(int(fi.sum()), int(len(fld)), R["mask_area_frac"]).pvalue)
fig, ax = plt.subplots(figsize=(9, 5.2), dpi=150)
ax.imshow(np.where(ok, r750, np.nan), extent=[LON0, LON1, LAT0, LAT1], cmap="gray", aspect=1 / np.cos(np.radians(50)))
ax.contour(lo, la, m.astype(float), levels=[0.5], colors="#ff9f1c", linewidths=1.0)
ax.scatter(fld.lon, fld.lat, s=10 + 2 * fld.n, c=np.where(fi, "#e63946", np.where(fn, "#ffd23f", "#4cc9f0")), edgecolors="k", linewidths=0.4, zorder=5)
ax.set_xlabel("East longitude (°)"); ax.set_ylabel("Latitude (°N)")
ax.set_title("750-nm reflectance within 1.2 b59 radii; orange: bright, red plains mask; fields red = inside, yellow = within 25 km, blue = farther", fontsize=7.5)
fig.savefig(f"{HERE}/fig_hrp_mask.png", bbox_inches="tight")
np.savez_compressed(f"{HERE}/hrp_mask.npz", mask=m, ok=ok, lo=lo, la=la)
fld.assign(cat=np.where(fi, "on", np.where(fn, "near", "far"))).to_csv(f"{HERE}/hrp_fields.csv", index=False)
json.dump(R, open(f"{HERE}/hrp_edge_test.json", "w"), indent=1); print(json.dumps(R, indent=1))
