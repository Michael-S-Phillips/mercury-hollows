# g4 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_ref/ref_size_insolation.py (6 Nov 2026). Changed: input/output paths by S2_CAT, added exponents/measures 1.7 and 1.8 (see make_copies_s2.py).
#!/usr/bin/env python3
"""Insolation results of v6.1 Section 3.4 with completeness at reference sizes 0.3-5 km in place of C1 (3 Nov 2026).
(a) crater-level hot/warm odds ratio of h09 in latitude bands (insolation_v6/band_crater_hotwarm.py spec:
h09 ~ hot + logD + C + abslat); (b) per 10 K of Delta T_w within 60 deg (h09 ~ dT10 + logD + C + abslat); (c) the same
without completeness; (d) mean C at hot and warm longitudes. Writes ref_size_insolation.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings; warnings.filterwarnings("ignore")
import os as _o, sys as _s; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__))); from s2cfg import *   # g4: inputs by S2_CAT
HERE = os.path.dirname(os.path.abspath(__file__)); pass
HE = os.path.expanduser("~/mnt/hollows-emissivity")
c = pd.read_csv(FEAT, low_memory=False, dtype={"id": str})
c["hot"] = hot_mask(c.lon.values).astype(int); c["dT10"] = c.dT / 10
Cz = np.load(GRID); Dref = list(Cz["Dref"])
rr = ((90 - c.lat.values) * 10).astype(int).clip(0, 1799); cc = ((c.lon.values % 360) * 10).astype(int).clip(0, 3599)
def orr(f, d, t):
    m = smf.logit(f, data=d).fit(disp=0); b, se = m.params[t], m.bse[t]
    return [round(float(np.exp(b)), 2), round(float(np.exp(b - 1.96 * se)), 2), round(float(np.exp(b + 1.96 * se)), 2), float(m.pvalues[t]), round(float(m.llf), 1)]
bands = [(-15, 15), (15, 30), (30, 45), (45, 60), (-60, 60)]
R = {}
for D0 in ["none", 0.3, 0.5, 1.0, 2.0, 5.0]:
    if D0 != "none": c["C"] = np.nan_to_num(Cz["C"][Dref.index(D0)][rr, cc])
    ctl = "logD + abslat" + ("" if D0 == "none" else " + C")
    out = {}
    for lo, hi in bands:
        d = c[(c.lat >= lo) & (c.lat < hi)]
        out[f"{lo}..{hi}"] = {"hot_OR": orr(f"h09 ~ hot + {ctl}", d, "hot")}
        if D0 != "none": out[f"{lo}..{hi}"]["meanC_hot_warm"] = [round(float(d[d.hot == 1].C.mean()), 3), round(float(d[d.hot == 0].C.mean()), 3)]
    d = c[c.abslat < 60]
    out["dT_per10K_lt60"] = orr(f"h09 ~ dT10 + {ctl}", d, "dT10")
    R[str(D0)] = out
    print(D0, " | ".join(f"{k}: {v['hot_OR'][:3]} p={v['hot_OR'][3]:.2g}" + (f" C h/w {v['meanC_hot_warm']}" if 'meanC_hot_warm' in v else "") for k, v in out.items() if k != "dT_per10K_lt60"), "| dT/10K<60:", out["dT_per10K_lt60"][:4], flush=True)
json.dump(R, open(f"{HERE}/out/ref_size_insolation_{CAT}.json", "w"), indent=1)
