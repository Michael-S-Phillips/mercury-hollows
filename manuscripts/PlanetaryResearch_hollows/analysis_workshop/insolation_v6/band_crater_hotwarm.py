#!/usr/bin/env python3
"""Supporting check for comment 1, 2 Nov 2026: crater-level hot/warm odds ratio of hosting (h09) in 15-deg latitude bands,
with completeness controlled, beside the field-level ratios of Table S6. Odds ratios for hot sector (within 45 deg of
0/180 E) from h09 ~ hot + logD + C1 + abslat in each band; variants: completeness as logit, well-imaged craters only
(C1 >= 0.2), the 90-150E / 40-70N province and Caloris < 1,200 km excluded, and per 10 K of Delta T_w. Also the
Delta T_w x logD interaction in each band. Writes band_crater_hotwarm.json."""
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from common_null import *
c = pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv", low_memory=False)
c["hot"] = hot_mask(c.lon.values).astype(int); c["dcal"] = gc_dist_km(CALORIS[0], CALORIS[1], c.lat.values, c.lon.values)
c["prov"] = ((c.lon >= 90) & (c.lon <= 150) & (c.lat >= 40) & (c.lat <= 70)).astype(int)
c["lc"] = np.log((c.C1 + 0.01) / (1.01 - c.C1)); c["dT10"] = c.dT / 10
def orr(f, d, t="hot"):
    try:
        m = smf.logit(f, data=d).fit(disp=0); b, se = m.params[t], m.bse[t]
        return [round(float(np.exp(b)), 2), round(float(np.exp(b - 1.96 * se)), 2), round(float(np.exp(b + 1.96 * se)), 2), float(m.pvalues[t])]
    except Exception as ex: return str(ex)
out = {}
bands = [(-30, -15), (-15, 0), (0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (15, 75), (-15, 15), (-60, 60)]
for lo, hi in bands:
    d = c[(c.lat >= lo) & (c.lat < hi)]
    dx = d[(d.prov == 0) & (d.dcal > CALORIS_KM)]
    r = dict(n=int(len(d)), h09=int(d.h09.sum()), h09_hot=int(d[d.hot == 1].h09.sum()), h09_warm=int(d[d.hot == 0].h09.sum()),
             meanC1_hot=round(float(d[d.hot == 1].C1.mean()), 3), meanC1_warm=round(float(d[d.hot == 0].C1.mean()), 3),
             raw=orr("h09 ~ hot + logD + abslat", d),
             C1=orr("h09 ~ hot + logD + C1 + abslat", d),
             logitC=orr("h09 ~ hot + logD + lc + abslat", d),
             C1ge02=orr("h09 ~ hot + logD + C1 + abslat", d[d.C1 >= 0.2]),
             n_C1ge02=int((d.C1 >= 0.2).sum()), h09_C1ge02=int(d[d.C1 >= 0.2].h09.sum()),
             C1_noProvCal=orr("h09 ~ hot + logD + C1 + abslat", dx), h09_noProvCal=int(dx.h09.sum()),
             hosts11_C1=orr("hosts ~ hot + logD + C1 + abslat", d),
             per10K_C1=orr("h09 ~ dT10 + logD + C1 + abslat", d, "dT10"),
             dTxlogD=orr("h09 ~ dT10 * logD + C1 + abslat", d, "dT10:logD"))
    out[f"{lo}..{hi}"] = r
    print(f"{lo:4d}..{hi:3d} n={r['n']:5d} h09={r['h09']:3d} (hot {r['h09_hot']}, warm {r['h09_warm']}) C1 hot/warm {r['meanC1_hot']}/{r['meanC1_warm']}\n"
          f"      raw {r['raw'][:3]}  C1 {r['C1'][:3]} p={r['C1'][3]:.2g}  logitC {r['logitC'][:3]}  C1>=0.2 {r['C1ge02'][:3]} ({r['h09_C1ge02']} hosts)  noProv/Cal {r['C1_noProvCal'][:3]} ({r['h09_noProvCal']})\n"
          f"      hosts1.1 {r['hosts11_C1'][:3]}  per10K {r['per10K_C1'][:3]}  dTxlogD coef-exp {r['dTxlogD'][:3]} p={r['dTxlogD'][3]:.2g}")
json.dump(out, open(f"{HERE}/band_crater_hotwarm.json", "w"), indent=1)
