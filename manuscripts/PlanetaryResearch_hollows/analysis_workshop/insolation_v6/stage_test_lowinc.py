#!/usr/bin/env python3
"""Add-on to stage_test.py, 2 Nov 2026: the same Stage 1 vs Stage 2 tests restricted to labels from images with incidence
< 60 deg (halo contrast visible; stage labels less affected by shadowing). Runs stage_test.py's setup by exec, then
fits linear, hot-sector, and step forms of Delta T_w, and the wall test. Writes stage_test_lowinc.json."""
import json, numpy as np, pandas as pd
src = open(__file__.replace("stage_test_lowinc.py", "stage_test.py")).read().split("OUT[\"results\"] = R")[0]
exec(src)
out = {}
for cut in (50, 60, 70):
    d = S[(S.abslat <= 60) & (S.inc8_m < cut)].copy()
    r = dict(n=int(len(d)), s1_frac=float(d.s1.mean()), dT10=glm(f"s1 ~ dT10 + {CTRL}", d, "dT10"), hot=glm(f"s1 ~ hot + {CTRL}", d, "hot"))
    a = {"base": smf.glm(f"s1 ~ {CTRL}", data=d, family=sm.families.Binomial()).fit().aic,
         "linear": smf.glm(f"s1 ~ dT10 + {CTRL}", data=d, family=sm.families.Binomial()).fit().aic}
    st = {}
    for t in range(20, 61, 2):
        d["st"] = (d.dT >= t).astype(float)
        if d.st.nunique() > 1: st[t] = smf.glm(f"s1 ~ st + {CTRL}", data=d, family=sm.families.Binomial()).fit().aic
    tb = min(st, key=st.get); a["step_best"] = st[tb] + 2; d["st"] = (d.dT >= tb).astype(float)
    r.update(AIC={k: round(v, 1) for k, v in a.items()}, best_t=tb, step=glm(f"s1 ~ st + {CTRL}", d, "st"),
             step_noCaloris=glm(f"s1 ~ st + {CTRL}", d[d.dcal > CALORIS_KM], "st"))
    e = pd.cut(d.dT, [0, 20, 30, 40, 50, 55, 60, 63])
    r["s1_share_by_dT"] = d.groupby(e, observed=True).agg(n=("s1", "size"), s1=("s1", "mean")).reset_index().astype({"dT": str}).to_dict("records")
    ww = w[w.inc8_m < cut]; r["walls_north_eqf"] = glm(f"s1 ~ eqf + {CTRL}", ww, "eqf")
    out[f"inc_lt{cut}"] = r
    print(f"== inc<{cut}: n={r['n']} s1 {r['s1_frac']:.2f}; dT10 OR {r['dT10']['OR']:.2f} ({r['dT10']['lo']:.2f}-{r['dT10']['hi']:.2f}) p={r['dT10']['p']:.2g};"
          f" hot {r['hot']['OR']:.2f} ({r['hot']['lo']:.2f}-{r['hot']['hi']:.2f}) p={r['hot']['p']:.2g}; AIC {r['AIC']} t={tb};"
          f" step {r['step']['OR']:.2f} ({r['step']['lo']:.2f}-{r['step']['hi']:.2f}) p={r['step']['p']:.2g}; noCal {r['step_noCaloris']['OR']:.2f} ({r['step_noCaloris']['lo']:.2f}-{r['step_noCaloris']['hi']:.2f})")
    print("    share", [(x['dT'], x['n'], round(x['s1'], 2)) for x in r["s1_share_by_dT"]], " walls eqf", {k: round(v, 3) if isinstance(v, float) else v for k, v in r["walls_north_eqf"].items()})
json.dump(out, open(f"{HERE}/stage_test_lowinc.json", "w"), indent=1, default=float)
