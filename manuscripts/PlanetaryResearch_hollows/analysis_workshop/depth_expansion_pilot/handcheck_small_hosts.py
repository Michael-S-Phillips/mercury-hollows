#!/usr/bin/env python3
"""Render the clean small-host (< 20 km) target measurements for a hand check (2 Nov 2026).
Writes qc/smallhost_sheet_{1,2,3}.png (18 tiles each, numbered k = 0..n-1 in new_depths order) and smallhost_index.csv.
Tile: Bickel box (orange), sun direction (yellow, toward the sun), measured shadow run (cyan)."""
import os, sys, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dep_common as C, qc_render as Q, targets as T
d = pd.read_csv(f"{HERE}/new_depths.csv", dtype={"image": str})
s = d[(d.qc_clean) & (d.group.str.contains("small_host"))].reset_index(drop=True)
s[["fid", "image", "group", "lat", "lon", "host_id", "host_D_km", "px_m", "inc_deg", "shadow_len_px", "depth_m", "depth_corr_cal_m", "flags"]].to_csv(f"{HERE}/smallhost_index.csv", index_label="k")
tiles = []
for k, r in s.iterrows():
    num = r.image[2:12]; p = C.edr_path(num)
    if not os.path.exists(p): C.fetch_edr(num, T.s_dt(num))
    img, info = C.read_edr(p)
    cap = f"k{k} f{r.fid} host {r.host_D_km:.0f} km\nh={r.depth_m:.0f} m ({r.shadow_len_px:.1f} px) px{r.px_m:.0f} i{r.inc_deg:.0f}\n{r['flags'] if isinstance(r['flags'], str) else ''}"
    tiles.append(Q.tile(img, info, (r.x0, r.x1, r.y0, r.y1), dict(run_u0=r.run_u0, run_u1=r.run_u1, run_w=r.run_w), cap, size=300))
for i in range(0, len(tiles), 18):
    print(Q.sheet(tiles[i:i + 18], ncol=6, size=300, path=f"{HERE}/qc/smallhost_sheet_{i // 18 + 1}.png"))
print(len(s))
