#!/usr/bin/env python3
"""Compare positions from georeference_bickel.py with the archived camera-model positions in
GIS/vector/bickel2024_hollows_corrected.gpkg (layer hollows_points_corrected), matched by image and bounding-box center.
Usage: python3 compare_to_archive.py <work_dir>/recomputed_positions.csv"""
import os, sys, sqlite3, numpy as np, pandas as pd
GPKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "GIS", "vector", "bickel2024_hollows_corrected.gpkg")
new = pd.read_csv(sys.argv[1])
new["img_x"] = (new.xmin + new.xmax) / 2; new["img_y"] = (new.ymin + new.ymax) / 2
old = pd.read_sql("select mdis_id, img_x, img_y, lat, lon180 from hollows_points_corrected", sqlite3.connect(GPKG))
old = old[old.mdis_id.isin(new.mdis_id)]
for c in ("img_x", "img_y"): new[c + "_k"] = new[c].round(5); old[c + "_k"] = old[c].round(5)
m = new.merge(old, on=["mdis_id", "img_x_k", "img_y_k"], suffixes=("_new", "_archive"))
R = 2439.4
la1, la2 = np.radians(m.lat_new), np.radians(m.lat_archive); dlo = np.radians(m.lon180_new - m.lon180_archive)
km = R * np.arccos(np.clip(np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(dlo), -1, 1))
print(f"matched {len(m)} of {len(new)} detections in {m.mdis_id.nunique()} images")
print(f"difference from archive: median {np.median(km):.3f} km, 95th percentile {np.percentile(km, 95):.3f} km, max {km.max():.3f} km")
print("PASS" if len(m) and np.median(km) < 0.05 else "CHECK: differences exceed 0.05 km; see the VERIFY notes in georeference_bickel.py")
