#!/usr/bin/env python3
"""Cache a basemap for 38-62N, 98-145E (5 Nov 2026): the USGS MDIS BDR 166 m mosaic averaged to 1/64 deg, and MDIS MDR
color (480, 750, 1000 nm) sampled at 1/32 deg. Writes basemap_cache.npz."""
import os, sys, numpy as np, rasterio
from rasterio.windows import Window
from rasterio.enums import Resampling
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../..")
sys.path.insert(0, f"{HERE}/../color_mdr"); import mdr_georef as g
LAT0, LAT1, LON0, LON1 = 38.0, 62.0, 98.0, 145.0
with rasterio.open(f"{ROOT}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif") as s:
    r0, r1 = int((90 - LAT1) * 256), int((90 - LAT0) * 256); c0, c1 = int((LON0 + 180) * 256), int((LON1 + 180) * 256)
    pan = s.read(1, window=Window(c0, r0, c1 - c0, r1 - r0), out_shape=((r1 - r0) // 4, (c1 - c0) // 4), resampling=Resampling.average).astype(np.float32)
step = 1 / 32; la = LAT1 - step * (np.arange(int((LAT1 - LAT0) / step)) + 0.5); lo = LON0 + step * (np.arange(int((LON1 - LON0) / step)) + 0.5)
LO, LA = np.meshgrid(lo, la); G = g.Globe()
col = np.stack([G.sample(LA.ravel(), LO.ravel(), band=b).reshape(LA.shape) for b in (1, 4, 7)])
np.savez_compressed(f"{HERE}/basemap_cache.npz", pan=pan, col=col, extent=np.array([LON0, LON1, LAT0, LAT1]))
print(pan.shape, col.shape)
