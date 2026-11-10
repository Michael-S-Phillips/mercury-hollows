#!/usr/bin/env python3
"""g23 step 10 (9 Nov): automatic shadow depths for every suitable hollow, for the latitude and color depth regressions with hollow size.
Selection: non-trace v3 hollows; detections (out_v2/detections_qc.csv) with pixel scale <= 20 m, incidence 58-81 deg, emission <= 50 deg;
the finest such detection of each hollow (2,012 hollows, 713 frames). Same measurement as depth_vs_size.py (pilot shadow_measure v2);
detections already measured there are copied, not re-measured. Output depth_all.csv. Resumable: run until it reports 0 left."""
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import depth_vs_size as DV
OUT = f"{HERE}/depth_all.csv"
def select():
    H = pd.read_csv(f"{DV.RR}/g2_insolation/inputs/v3/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool)]
    d = pd.read_csv(f"{DV.RUN}/out_v2/detections_qc.csv"); d = d[d.hollow_id.isin(H.hollow_id)]
    d = d.merge(H[["hollow_id", "D_eq_km"]].rename(columns={"D_eq_km": "D_can_km"}), on="hollow_id")
    d = d[(d.px_m <= 20) & d.inc.between(58, 81) & (d.emi <= 50)].sort_values("px_m").drop_duplicates("hollow_id").copy(); d["tier"] = "all"
    b = DV.C.load_bickel(with_corners=False)[["thumb", "num", "datetime", "x0", "x1", "y0", "y1", "lat", "lon180"]].drop_duplicates("thumb")
    s = d.drop(columns=["lat"]).merge(b.rename(columns={"lat": "lat_cm"}), on="thumb", how="left"); s["num"] = s.num.astype(str).str.zfill(10)
    return s.reset_index(drop=True)
if not os.path.exists(OUT):
    s = select(); old = pd.read_csv(f"{HERE}/depth_vs_size.csv"); old = old[old.thumb.isin(s.thumb)].drop_duplicates("thumb").copy(); old["tier"] = "all"
    old.to_csv(OUT, index=False); print("copied", len(old))
DV.select = select; DV.OUT = OUT
DV.run()
