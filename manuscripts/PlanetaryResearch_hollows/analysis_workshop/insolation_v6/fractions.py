#!/usr/bin/env python3
"""Task (e), 2 Nov 2026: share of hollow fields, host craters, and hollows within 15 deg of the equator (where the
field-level hot-longitude excess holds) and at 15-75 N. Fields: code/rerun/out_v2/fields.csv; hollows: non-trace canonical
hollows; hosts: crater_features.csv (h09, and hosts within 1.1 radii). Also with Caloris < 1,200 km excluded.
Writes fractions.json."""
import json, numpy as np, pandas as pd
from common_null import *
cov, comp, fields, hollows = load()
H = hollows[~hollows.is_crater_trace.astype(bool)]
c = pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv", low_memory=False)
bands = {"|lat|<15": lambda la: np.abs(la) < 15, "0-15N": lambda la: (la >= 0) & (la < 15), "0-15S": lambda la: (la < 0) & (la > -15),
         "15-30N": lambda la: (la >= 15) & (la < 30), "15-75N": lambda la: (la >= 15) & (la < 75), "30-75N": lambda la: (la >= 30) & (la < 75),
         "15-75S": lambda la: (la <= -15) & (la > -75), ">=75 both": lambda la: np.abs(la) >= 75}
out = {}
for nm, df, la, lo in (("fields", fields, fields.lat.values, fields.lon.values), ("hollows", H, H.lat.values, H.lon.values),
                       ("h09_hosts", c[c.h09 == 1], c[c.h09 == 1].lat.values, c[c.h09 == 1].lon.values),
                       ("hosts_1.1R", c[c.hosts == 1], c[c.hosts == 1].lat.values, c[c.hosts == 1].lon.values),
                       ("craters_ge5km", c, c.lat.values, c.lon.values)):
    nc = gc_dist_km(CALORIS[0], CALORIS[1], la, lo) > CALORIS_KM
    out[nm] = {"total": int(len(la)), "total_noCaloris": int(nc.sum())}
    for b, f in bands.items():
        m = f(la); out[nm][b] = dict(n=int(m.sum()), frac=round(float(m.mean()), 3), n_noCal=int((m & nc).sum()), frac_noCal=round(float((m & nc).sum() / nc.sum()), 3))
    print(nm, out[nm]["total"], {b: (v["n"], v["frac"], v["frac_noCal"]) for b, v in out[nm].items() if isinstance(v, dict)})
json.dump(out, open(f"{HERE}/fractions.json", "w"), indent=1)
