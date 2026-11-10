#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script. Writes the hollow masks used by the color and knob reruns.
"""inputs/hollows_mask_{v3,v3b}.csv: the cataloged hollows that remain in each screened catalog, with the columns that
color_mdr/extract_tiles.py reads (lat, lon, D_eq_km). The v2 mask is the 6,731 rows of code/rerun/out_v2/hollows_canonical.csv
(all rows, including the 64 v2 crater traces), so the screened masks keep those 64 traces too and drop only the screened
detections: v3 drops the 368 spurious detections (B, W, O) and v3b also the 232 uncertain ones (U)."""
import os, json, numpy as np, pandas as pd
R = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun")
HERE = os.path.dirname(os.path.abspath(__file__))
v2 = pd.read_csv(f"{R}/out_v2/hollows_canonical.csv")
rec = {"v2_rows": int(len(v2))}
for cat in ["v3", "v3b"]:
    c = pd.read_csv(f"{R}/out_v3/hollows_canonical_{cat}.csv")
    assert (c.hollow_id.values == v2.hollow_id.values).all() and np.allclose(c.lat, v2.lat) and np.allclose(c.lon, v2.lon)
    assert np.allclose(c.D_eq_km.fillna(-1), v2.D_eq_km.fillna(-1))
    drop = c.is_spurious_v3.astype(bool) | (c.is_uncertain_v3.astype(bool) if cat == "v3b" else False)
    m = c.loc[~drop, ["hollow_id", "lat", "lon", "D_eq_km", "is_crater_trace"]]
    m.to_csv(f"{HERE}/inputs/hollows_mask_{cat}.csv", index=False)
    rec[cat] = dict(mask_rows=int(len(m)), dropped=int(drop.sum()), non_trace=int((~m.is_crater_trace.astype(bool)).sum()))
json.dump(rec, open(f"{HERE}/inputs/prep_summary.json", "w"), indent=1); print(rec)
