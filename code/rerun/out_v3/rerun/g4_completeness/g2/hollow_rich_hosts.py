# g4 copy of code/rerun/out_v3/rerun/g2_insolation/hollow_rich_hosts.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g4.py).
#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/insolation_v6/hollow_rich_hosts.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Supporting check for comment 1, 2 Nov 2026: at 15-75 N, do a few hollow-rich host craters at the warm longitudes make
individual hollow counts favor the warm longitudes while host craters favor the hot longitudes? Uses n09 (hollows within
0.9 radii) in crater_features.csv. Writes hollow_rich_hosts.json."""
import json, numpy as np, pandas as pd
from common_null import *
c = pd.read_csv(f"{IN}/crater_features.csv", low_memory=False)
c["hot"] = hot_mask(c.lon.values).astype(int); c["prov"] = ((c.lon >= 90) & (c.lon <= 150) & (c.lat >= 40) & (c.lat <= 70)).astype(int)
h = c[(c.h09 == 1) & (c.lat >= 15) & (c.lat < 75)].sort_values("n09", ascending=False)
out = dict(n_hosts=int(len(h)), hosts_hot=int(h.hot.sum()), hosts_warm=int((1 - h.hot).sum()),
           hollows_hot=int(h[h.hot == 1].n09.sum()), hollows_warm=int(h[h.hot == 0].n09.sum()),
           top10=h.head(10)[["id", "lat", "lon", "D", "n09", "hot", "prov"]].round(2).to_dict("records"),
           top10_share_of_hollows=float(h.head(10).n09.sum() / h.n09.sum()),
           without_top10=dict(hollows_hot=int(h.iloc[10:][h.iloc[10:].hot == 1].n09.sum()), hollows_warm=int(h.iloc[10:][h.iloc[10:].hot == 0].n09.sum())),
           median_n09_hot=float(h[h.hot == 1].n09.median()), median_n09_warm=float(h[h.hot == 0].n09.median()))
print(json.dumps(out, indent=0, default=float))
json.dump(out, open(f"{G2_OUT}/hollow_rich_hosts.json", "w"), indent=1, default=float)
