#!/usr/bin/env python3
"""Contact sheets for the crater-bowl screen (6 Nov 2026). Each tile: MDIS BDR 166 m mosaic around a canonical hollow,
yellow circle = detection equivalent diameter at its centroid, orange = Herrick craters >= 5 km. Usage: qc_sheet.py SET N SEED"""
import sys, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
from PIL import Image, ImageDraw
ROOT = "/sessions/" if False else None
import os; HE = os.path.expanduser("~/mnt/hollows-emissivity"); OUT = f"{HE}/code/rerun/out_v3"
c = pd.read_csv(f"{OUT}/bowl_candidates.csv"); cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon"] = cr.lon % 360
which, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
sel = {"flag": c[c.bowl_cr >= 0], "border": c[c.borderline], "big": c[c.big_unmatched], "small": c[(c.D_eq_km < 1) & (c.bowl_cr < 0)]}[which]
if len(sel) > N: sel = sel.sample(N, random_state=seed)
KD = 2439.4 * np.pi / 180; T = 240; tiles = []
with rasterio.open(f"{HE}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif") as ds:
    for _, r in sel.iterrows():
        half = max(4 * r.D_eq_km, 8.0, 0.8 * (r.bowl_D if r.bowl_D == r.bowl_D else 0)); lat, lon = r.lat, ((r.lon + 180) % 360) - 180
        dla = half / KD; dlo = half / (KD * max(np.cos(np.radians(lat)), 0.05))
        b = ds.bounds
        x0 = b.left + (lon - dlo + 180) / 360 * (b.right - b.left); x1 = b.left + (lon + dlo + 180) / 360 * (b.right - b.left)
        y0 = b.bottom + (lat - dla + 90) / 180 * (b.top - b.bottom); y1 = b.bottom + (lat + dla + 90) / 180 * (b.top - b.bottom)
        try: a = ds.read(1, window=from_bounds(x0, y0, x1, y1, ds.transform), boundless=True, fill_value=0).astype(float)
        except Exception: continue
        if a.size == 0: continue
        lo_, hi_ = np.percentile(a[a > 0], [1, 99.5]) if (a > 0).any() else (0, 1)
        im = Image.fromarray((np.clip((a - lo_) / (hi_ - lo_ + 1e-6), 0, 1) * 255).astype(np.uint8)).resize((T, T)).convert("RGB"); d = ImageDraw.Draw(im)
        s = T / (2 * half)
        def circ(la, lo_c, Dk, col, w):
            dx = ((lo_c - r.lon + 180) % 360 - 180) * KD * np.cos(np.radians(lat)); dy = (la - lat) * KD
            cx, cy = T / 2 + dx * s, T / 2 - dy * s; rr = Dk / 2 * s; d.ellipse([cx - rr, cy - rr, cx + rr, cy + rr], outline=col, width=w)
        near = cr[(abs(cr.lat - lat) < 2 * dla) & (((cr.lon - r.lon + 180) % 360 - 180).abs() < 2 * dlo)]
        for _, q in near.iterrows(): circ(q.lat, q.lon, q.diam, (255, 120, 0), 1)
        circ(r.lat, r.lon, r.D_eq_km, (255, 240, 0), 2)
        T2 = Image.new("RGB", (T, T + 30)); T2.paste(im, (0, 0)); dd = ImageDraw.Draw(T2)
        dd.text((2, T + 2), f"{int(r.hollow_id)} Deq {r.D_eq_km:.1f} km px {r.px_best:.0f} m", fill=(0, 255, 0))
        dd.text((2, T + 15), f"crD {r.bowl_D:.1f} rel {r.rel:.2f} n_det {int(r.n_det)}" if r.bowl_cr >= 0 else f"n_det {int(r.n_det)} field {int(r.field_id)}", fill=(0, 255, 0))
        tiles.append(T2)
nc = 8; nr = int(np.ceil(len(tiles) / nc)); S = Image.new("RGB", (nc * T, nr * (T + 30)))
for i, t in enumerate(tiles): S.paste(t, ((i % nc) * T, (i // nc) * (T + 30)))
S.save(f"{OUT}/qc/sheet_{which}.png"); sel.to_csv(f"{OUT}/qc/sheet_{which}_ids.csv", index=False); print(which, len(tiles), S.size)
