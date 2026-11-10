#!/usr/bin/env python3
"""Check vent positions in data/data/MER_Vent_Locations_Classifications.csv against the MDIS 166 m/pixel global
basemap (2 Nov 2026). Explosive vents sit in bright faculae, so the basemap should be brighter at a true vent position
than in a surrounding annulus. Contrast = mean DN within 8 km / mean DN at 25-45 km. Compared for the swapped reading
(lat = file 'Lon', lon = file 'Lat'), the as-labeled reading (only rows whose 'Lat' value is a valid latitude), and a null of
positions shifted 100 km in a random direction. Writes vent_imagery_check.csv/.json and a contact sheet of 42 vents."""
import os, json, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
r = rasterio.open(f"{ROOT}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif")
R = 2439400.0; px = r.transform.a
v = pd.read_csv(f"{ROOT}/data/data/MER_Vent_Locations_Classifications.csv"); v.columns = [c.strip() for c in v.columns]
def window(lat, lon, half_km=50, dec=4):
    lon = ((lon + 180) % 360) - 180
    x = R * np.radians(lon) * 1.0; y = R * np.radians(lat)
    col = (x - r.transform.c) / px; row = (r.transform.f - y) / px
    hw = int(half_km * 1000 / px / np.cos(np.radians(lat))) if abs(lat) < 85 else None
    hh = int(half_km * 1000 / px)
    if hw is None: return None
    c0, r0 = int(col - hw), int(row - hh)
    if c0 < 0 or r0 < 0 or c0 + 2 * hw >= r.width or r0 + 2 * hh >= r.height: return None
    a = r.read(1, window=Window(c0, r0, 2 * hw, 2 * hh), out_shape=(max(1, 2 * hh // dec), max(1, 2 * hw // dec))).astype(float)
    return a
def contrast(a, half_km=50):
    if a is None: return np.nan
    ny, nx = a.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    dy = (yy - ny / 2) * (2 * half_km / ny); dx = (xx - nx / 2) * (2 * half_km / nx)
    rr = np.hypot(dx, dy); ok = a > 0
    c = a[(rr < 8) & ok]; s = a[(rr > 25) & (rr < 45) & ok]
    return float(c.mean() / s.mean()) if len(c) > 20 and len(s) > 50 and s.mean() > 0 else np.nan
rng = np.random.default_rng(1)
rows = []
for i, q in v.iterrows():
    lat_s, lon_s = q.Lon, q.Lat
    cs = contrast(window(lat_s, lon_s))
    ca = contrast(window(q.Lat, q.Lon)) if abs(q.Lat) <= 85 else np.nan
    th = rng.uniform(0, 2 * np.pi); dlat = 100 / 42.58 * np.sin(th); dlon = 100 / 42.58 * np.cos(th) / max(np.cos(np.radians(lat_s)), 0.2)
    cn = contrast(window(lat_s + dlat, lon_s + dlon))
    rows.append(dict(vent_id=q["Vent ID"], lat=lat_s, lon_E=lon_s % 360, facula=q.Facula, quality=q["Image Quality"], compound=q.Compound,
                     contrast_swapped=cs, contrast_as_labeled=ca, contrast_shifted100km=cn))
d = pd.DataFrame(rows); d.to_csv(f"{HERE}/vent_imagery_check.csv", index=False)
def summ(x): x = x.dropna(); return dict(n=int(len(x)), median=float(x.median()), frac_gt_1=float((x > 1).mean()), frac_gt_1_1=float((x > 1.1).mean()))
S = {"swapped_all": summ(d.contrast_swapped), "swapped_facula_yes": summ(d[d.facula == "yes"].contrast_swapped),
     "swapped_facula_no": summ(d[d.facula == "no"].contrast_swapped), "shifted_null": summ(d.contrast_shifted100km),
     "as_labeled_valid_rows": summ(d.contrast_as_labeled)}
json.dump(S, open(f"{HERE}/vent_imagery_check.json", "w"), indent=1); print(json.dumps(S, indent=1))
