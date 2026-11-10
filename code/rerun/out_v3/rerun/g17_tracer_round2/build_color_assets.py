# g17 copy of g12_tracer/build_color_assets.py (7 Nov 2026). Changed: craters from selection.json; writes g17_tracer_round2/assets/col and ecol.
#!/usr/bin/env python3
"""g12 (7 Nov 2026): MDIS color basemaps for the tracer, at Michael's request (fewer shadows than the 166 m monochrome mosaic).
For each tracer crater (assets/order.json), samples the MD3 color and enhanced color global mosaics (GIS/raster, 665 m, 3 bands,
lon_0 = 0) bilinearly onto the tracer's local azimuthal-equidistant grid (half-width 1.3 R) at 333 m pixels, stretches the three bands
together between their joint 1st and 99.5th percentiles (keeping each product's color balance), and writes assets/col/<uid>.jpg (MD3)
and assets/ecol/<uid>.jpg (enhanced color). Grid geometry matches build_tracer_assets.py. Usage: build_color_assets.py"""
import os, json, numpy as np, rasterio
from rasterio.windows import Window
from scipy.ndimage import map_coordinates
from PIL import Image
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/assets"
SRC = {"col": f"{HE}/GIS/raster/Mercury_MESSENGER_MDIS_Basemap_MD3Color_Mosaic_Global_665m_cog.tif",
       "ecol": f"{HE}/GIS/raster/Mercury_MESSENGER_MDIS_Basemap_EnhancedColor_Mosaic_Global_665m_cog.tif"}
RM = 2439.4; RMm = 2439400.0; PX = 0.333
def aeqd_inv(x, y, la0, lo0):
    rho = np.hypot(x, y); c = rho / RM; p0 = np.radians(la0)
    with np.errstate(invalid="ignore", divide="ignore"):
        lat = np.arcsin(np.cos(c) * np.sin(p0) + np.where(rho > 0, y * np.sin(c) * np.cos(p0) / rho, 0))
        lon = lo0 + np.degrees(np.arctan2(x * np.sin(c), rho * np.cos(p0) * np.cos(c) - y * np.sin(p0) * np.sin(c)))
    return np.degrees(lat), lon
def read_wrapped(ds, band, r0, r1, c0, c1):
    W = ds.width; r0c, r1c = max(r0, 0), min(r1, ds.height); parts = []; c = c0
    while c < c1:
        cm = c % W; n = min(c1 - c, W - cm); parts.append(ds.read(band, window=Window(cm, r0c, n, r1c - r0c))); c += n
    a = np.concatenate(parts, axis=1)
    if r0c > r0 or r1c < r1: a = np.pad(a, ((r0c - r0, r1 - r1c), (0, 0)))
    return a
SEL = json.load(open(f"{HERE}/selection.json")); order = [g["uid"] for k in ("group_B", "group_C", "group_D") for g in SEL[k]]; done = 0
srcs = {k: rasterio.open(v) for k, v in SRC.items()}
for uid in order:
    m = json.load(open(f"{OUT}/meta/{uid}.json")); la0, lo0, H = m["lat"], m["lon"], m["H_km"]
    N = int(np.ceil(2 * H / PX)); px = 2 * H / N; g = -H + (np.arange(N) + .5) * px; X, Y = np.meshgrid(g, -g)
    lat, lon = aeqd_inv(X, Y, la0, lo0); lp = ((lo0 + 180) % 360) - 180; lonc = lp + (lon - lo0)
    for key, ds in srcs.items():
        T = ds.transform; col = (np.radians(lonc) * RMm - T.c) / T.a - 0.5; row = (np.radians(lat) * RMm - T.f) / T.e - 0.5
        c0, c1 = int(np.floor(col.min())) - 2, int(np.ceil(col.max())) + 3; r0, r1 = int(np.floor(row.min())) - 2, int(np.ceil(row.max())) + 3
        bands = [map_coordinates(read_wrapped(ds, b, r0, r1, c0, c1).astype(np.float32), [row - r0, col - c0], order=1, mode="nearest") for b in (1, 2, 3)]
        A = np.stack(bands, -1); ok = A.max(-1) > 0.5
        lo_, hi_ = np.percentile(A[ok], [1, 99.5]) if ok.any() else (0, 255)
        im = np.clip((A - lo_) / max(hi_ - lo_, 1) * 255, 0, 255).astype(np.uint8); im[~ok] = 0
        Image.fromarray(im, "RGB").save(f"{OUT}/{key}/{uid}.jpg", quality=88, optimize=True)
    done += 1
for ds in srcs.values(): ds.close()
print("craters", done)
