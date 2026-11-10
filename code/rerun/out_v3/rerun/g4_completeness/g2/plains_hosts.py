# g4 copy of code/rerun/out_v3/rerun/g2_insolation/plains_hosts.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g4.py).
#!/usr/bin/env python3
"""g2 (6 Nov 2026), new script: Section S6 counts of host craters on the northern smooth plains, from the crater table of the
containment tests (inputs/<cat>/crater_features.csv: 31,554 craters >= 5 km; nsp = center on the plains mask; h09 = >= 1
hollow within 0.9 radii). Writes out/<cat>/plains_hosts.json."""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CAT = os.environ.get("G2_CAT", "v3")
c = pd.read_csv(f"{HERE}/inputs/{CAT}/crater_features.csv", low_memory=False)
p = c[c.nsp == 1]; h = p[p.h09 == 1]; o = c[c.nsp == 0]
R = dict(n_plains=int(len(p)), plains_h09=int(len(h)), plains_h09_ge10km=int((h.D >= 10).sum()),
         median_D_plains_h09=float(h.D.median()), median_D_plains_no_h09=float(p[p.h09 == 0].D.median()),
         plains_hosts_any=int(p.hosts.sum()), plains_hosts_any_ge10km=int((p[p.hosts == 1].D >= 10).sum()),
         by_size={f"{lo}-{hi}": dict(plains=f"{int(p[(p.D >= lo) & (p.D < hi)].h09.sum())}/{int(((p.D >= lo) & (p.D < hi)).sum())}",
                                     off=f"{int(o[(o.D >= lo) & (o.D < hi)].h09.sum())}/{int(((o.D >= lo) & (o.D < hi)).sum())}")
                  for lo, hi in ((5, 10), (10, 20), (20, 40), (40, 1e5))},
         plains_h09_list=h.sort_values("D")[["id", "lat", "lon", "D", "n09"]].round(2).to_dict("records"))
os.makedirs(f"{HERE}/out/{CAT}", exist_ok=True); json.dump(R, open(f"{HERE}/out/{CAT}/plains_hosts.json", "w"), indent=1)
print({k: v for k, v in R.items() if k != "plains_h09_list"})
