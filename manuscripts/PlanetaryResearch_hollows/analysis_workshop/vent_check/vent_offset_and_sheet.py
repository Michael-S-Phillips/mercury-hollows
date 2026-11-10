#!/usr/bin/env python3
"""Positional check of the swapped vent reading (2 Nov 2026): for each vent, the offset (within +-20 km) of the brightest
8-km-radius disk relative to its 25-45 km annulus on the MDIS 166 m basemap; and a contact sheet of 42 random vents
(+-40 km, crosshair at the catalog position). Writes vent_offsets.csv/.json and vent_sheet.png."""
import os, json, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
from PIL import Image, ImageDraw
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
r = rasterio.open(f"{ROOT}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif"); R = 2439400.0; px = r.transform.a
v = pd.read_csv(f"{ROOT}/data/data/MER_Vent_Locations_Classifications.csv"); v.columns = [c.strip() for c in v.columns]
v["lat"], v["lon"] = v.Lon, v.Lat
def read(lat, lon, half_km, dec):
    lon = ((lon + 180) % 360) - 180; cl = np.cos(np.radians(lat))
    col = (R * np.radians(lon) - r.transform.c) / px; row = (r.transform.f - R * np.radians(lat)) / px
    hw = int(half_km * 1000 / px / cl); hh = int(half_km * 1000 / px)
    c0, r0 = int(col - hw), int(row - hh)
    if c0 < 0 or r0 < 0 or c0 + 2 * hw >= r.width or r0 + 2 * hh >= r.height: return None
    return r.read(1, window=Window(c0, r0, 2 * hw, 2 * hh), out_shape=(2 * hh // dec, 2 * hw // dec)).astype(float)
out = []
for _, q in v.iterrows():
    if abs(q.lat) > 80: out.append(dict(vent_id=q["Vent ID"], dx_km=np.nan, dy_km=np.nan)); continue
    a = read(q.lat, q.lon, 70, 6)
    if a is None: out.append(dict(vent_id=q["Vent ID"], dx_km=np.nan, dy_km=np.nan)); continue
    ny, nx = a.shape; kx, ky = 140 / nx, 140 / ny
    yy, xx = np.mgrid[0:ny, 0:nx]; best = (-1, 0, 0)
    for oy in np.arange(-20, 21, 2.0):
        for ox in np.arange(-20, 21, 2.0):
            rr = np.hypot((xx - nx / 2) * kx - ox, (yy - ny / 2) * ky - oy)
            c = a[(rr < 8) & (a > 0)]; s = a[(rr > 25) & (rr < 45) & (a > 0)]
            if len(c) > 20 and len(s) > 50:
                val = c.mean() / s.mean()
                if val > best[0]: best = (val, ox, oy)
    out.append(dict(vent_id=q["Vent ID"], best_contrast=best[0], dx_km=best[1], dy_km=-best[2]))
o = pd.DataFrame(out); o["dist_km"] = np.hypot(o.dx_km, o.dy_km); o.to_csv(f"{HERE}/vent_offsets.csv", index=False)
x = o.dropna(subset=["dist_km"])
S = dict(n=int(len(x)), median_offset_km=float(x.dist_km.median()), frac_within_4km=float((x.dist_km <= 4).mean()),
         frac_within_8km=float((x.dist_km <= 8).mean()), frac_at_search_edge=float((x.dist_km >= 19).mean()),
         mean_dx_km=float(x.dx_km.mean()), mean_dy_km=float(x.dy_km.mean()))
json.dump(S, open(f"{HERE}/vent_offsets.json", "w"), indent=1); print(json.dumps(S, indent=1))
# contact sheet
smp = v[(v.lat.abs() < 75)].sample(42, random_state=3); tiles = []
for _, q in smp.iterrows():
    a = read(q.lat, q.lon, 40, 3)
    if a is None: continue
    lo, hi = np.percentile(a[a > 0], [1, 99.5]) if (a > 0).any() else (0, 1)
    im = Image.fromarray((np.clip((a - lo) / (hi - lo + 1e-6), 0, 1) * 255).astype(np.uint8)).resize((240, 240)).convert("RGB")
    d = ImageDraw.Draw(im); d.line([110, 120, 100, 120], fill=(255, 0, 0)); d.line([130, 120, 140, 120], fill=(255, 0, 0))
    d.line([120, 110, 120, 100], fill=(255, 0, 0)); d.line([120, 130, 120, 140], fill=(255, 0, 0))
    d.text((3, 3), f"V{q['Vent ID']} {q.lat:.1f},{q.lon % 360:.1f}E {q.Facula[:1]}/{q['Image Quality'][:4]}", fill=(255, 255, 0))
    tiles.append(im)
S2 = Image.new("RGB", (7 * 240, 6 * 240))
for i, t in enumerate(tiles): S2.paste(t, ((i % 7) * 240, (i // 7) * 240))
S2.save(f"{HERE}/vent_sheet.png"); print("sheet", len(tiles))
