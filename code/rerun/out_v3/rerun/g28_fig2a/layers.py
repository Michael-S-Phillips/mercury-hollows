#!/usr/bin/env python3
"""g28 (10 Nov 2026): layers for a redesigned Figure 2a (north polar map), on a Lambert azimuthal equal-area grid centered on the
north pole (2 km pixels, +-2,600 km, as nsp_laea.npz but finer).
(1) basemap: the MDIS 166 m BDR mosaic (data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif) north of 27 N, read at every
    4th pixel, averaged 2 x 2 (1.33 km), and sampled at each grid pixel (nearest).
(2) footprints: best (finest) pixel scale of the images that cover each grid pixel, from the coverage grid of the analysis
    (code/rerun/out_v2/coverage_10ppd.npz, bestpx: MDIS NAC CDR, orbital, <= 150 m/pixel, incidence < 85 deg; each image a latitude-
    longitude box; rebuilt identically by g4_completeness/build_completeness.py), sampled at each grid pixel (nearest 0.1 deg cell). The
    MDIS index on disk has no corner coordinates for NAC images (placeholders -99.999), so true footprint outlines are not available.
Writes _cache/basemap_laea2km.npy, _cache/bestpx_laea2km.npy, layers_record.json."""
import os, json, time, numpy as np, pandas as pd, pyproj, rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window
from PIL import Image, ImageDraw
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__)); C = f"{HERE}/_cache"; os.makedirs(C, exist_ok=True)
R = 2439400.0; PIX = 2000.0; HALF = 2.6e6; N = int(2 * HALF / PIX)
ps = pyproj.CRS.from_proj4(f"+proj=laea +lat_0=90 +lon_0=0 +R={R} +no_defs"); geo = pyproj.CRS.from_proj4(f"+proj=longlat +R={R} +no_defs")
xx = (np.arange(N) + .5) * PIX - HALF; X, Y = np.meshgrid(xx, -xx)
lon, lat = pyproj.Transformer.from_crs(ps, geo, always_xy=True).transform(X, Y)
rec = dict(grid=dict(pix_m=PIX, half_m=HALF, n=N, proj="laea lat_0=90 lon_0=0 R=2439400"))
# (1) basemap
t0 = time.time()
ds = rasterio.open(f"{HE}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif"); ppd = ds.width / 360
rows = int((90 - 27.0) * ppd); f = 4
a = ds.read(1, window=Window(0, 0, ds.width, rows), out_shape=(rows // f, ds.width // f), resampling=Resampling.nearest).astype(np.float32)
a[a == 0] = np.nan
a = np.nanmean(a[: a.shape[0] // 2 * 2, : a.shape[1] // 2 * 2].reshape(a.shape[0] // 2, 2, a.shape[1] // 2, 2), axis=(1, 3))
pp = ppd / f / 2   # pixels per degree of the reduced mosaic (32)
r = np.clip(((90 - lat) * pp).astype(int), 0, a.shape[0] - 1); c = np.clip(((np.mod(lon + 180, 360)) * pp).astype(int), 0, a.shape[1] - 1)
bm = np.where(lat >= 27.0, a[r, c], np.nan).astype(np.float32)
np.save(f"{C}/basemap_laea2km.npy", bm); rec["basemap"] = dict(source="MDIS BDR 166 m mosaic", reduced_ppd=pp, seconds=round(time.time() - t0, 1))
# (2) footprints
cov = np.load(f"{HE}/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"].astype(np.float32)
rr = np.clip(((90 - lat) * 10).astype(int), 0, 1799); cc = np.clip((np.mod(lon, 360) * 10).astype(int), 0, 3599)
bp = np.where(lat >= 27.0, cov[rr, cc], np.nan).astype(np.float32)
np.save(f"{C}/bestpx_laea2km.npy", bp)
g30 = lat >= 30
rec["footprints"] = dict(source="code/rerun/out_v2/coverage_10ppd.npz bestpx", frac_ge30N_imaged=float(np.isfinite(bp[g30]).mean()),
                         frac_ge30N_le40=float((bp[g30] <= 40).mean()), median_px_ge30N=float(np.nanmedian(bp[g30])))
json.dump(rec, open(f"{HERE}/layers_record.json", "w"), indent=1); print(json.dumps(rec, indent=1))
