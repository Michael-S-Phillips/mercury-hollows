#!/usr/bin/env python3
"""g2 add-on (8 Nov 2026), for Michael's question about the latitude dependence of the hot-longitude preference: the Stage 1 vs
Stage 2 hot-sector test of stage_test.py by 15-degree latitude band (north), with labels from images at incidence < 70 deg
(the labels the paper uses) and, for comparison, all labels. Runs stage_test.py's setup by exec, as stage_test_lowinc.py does.
Writes out/<cat>/stage_by_band.json. Not for the paper."""
import json, numpy as np, pandas as pd
src = open(__file__.replace("stage_by_band.py", "stage_test.py")).read().split("OUT[\"results\"] = R")[0]
exec(src)
out = {}
for lab, cut in (("inc_lt70", 70), ("all_labels", 999)):
    out[lab] = {}
    for lo, hi in ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (30, 60), (30, 75), (0, 30)):
        d = S[(S.lat >= lo) & (S.lat < hi) & (S.inc8_m < cut)].copy()
        if lo == 0: d = S[(S.abslat < hi) & (S.inc8_m < cut)].copy() if hi <= 15 else S[(S.lat > -15) & (S.lat < hi) & (S.inc8_m < cut)].copy()
        r = dict(n=int(len(d)), n_fields=int(d.field_id.nunique()),
                 s1_hot=int(d[(d.hot == 1)].s1.sum()), n_hot=int((d.hot == 1).sum()), s1_warm=int(d[(d.hot == 0)].s1.sum()), n_warm=int((d.hot == 0).sum()),
                 mean_dT_hot=float(d[d.hot == 1].dT.mean()) if (d.hot == 1).any() else None, mean_dT_warm=float(d[d.hot == 0].dT.mean()) if (d.hot == 0).any() else None,
                 mean_inc_hot=float(d[d.hot == 1].inc8_m.mean()) if (d.hot == 1).any() else None, mean_inc_warm=float(d[d.hot == 0].inc8_m.mean()) if (d.hot == 0).any() else None)
        try: r["hot"] = glm(f"s1 ~ hot + {CTRL}", d, "hot")
        except Exception as e: r["hot"] = str(e)[:120]
        out[lab][f"{lo}-{hi}" if lo else (f"|lat|<{hi}" if hi <= 15 else f"-15-{hi}")] = r
json.dump(out, open(f"{G2_OUT}/stage_by_band.json", "w"), indent=1, default=float)
for lab in out:
    print("==", lab)
    for k, r in out[lab].items():
        h = r["hot"]; s = f"OR {h['OR']:.2f} ({h['lo']:.2f}-{h['hi']:.2f}) p={h['p']:.3f}" if isinstance(h, dict) else h
        f = lambda a, b: f"{a}/{b}={a / b:.2f}" if b else "-"
        print(f"  {k:9s} n={r['n']:5d} fields={r['n_fields']:4d} S1 hot {f(r['s1_hot'], r['n_hot'])}  warm {f(r['s1_warm'], r['n_warm'])}  dT hot/warm {r['mean_dT_hot'] or 0:.0f}/{r['mean_dT_warm'] or 0:.0f}  inc {r['mean_inc_hot'] or 0:.0f}/{r['mean_inc_warm'] or 0:.0f}  {s}")
