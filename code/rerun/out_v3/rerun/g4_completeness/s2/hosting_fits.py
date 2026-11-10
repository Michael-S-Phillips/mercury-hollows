# g4 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/hosting_fits.py (6 Nov 2026). Changed: input/output paths by S2_CAT, added exponents/measures 1.7 and 1.8 (see make_copies_s2.py).
#!/usr/bin/env python3
"""The paper's hosting fits with each completeness measure in place of C1 (4 Nov 2026). Specifications as in
completeness_ref/ref_size_check.py (nesting, vents, Delta T_w in the nesting fit, fresh vs degraded >= 10 km) and
ref_size_insolation.py (hot vs warm longitudes by latitude band; Delta T_w per 10 K within 60 deg). Measures:
C1 (paper), lnC1 = ln(C1 + 0.01) (same measure, log form), C1_act (own incidence and pixel scale), lnE_a = ln(E_a + 0.01)
for size-distribution slopes a = 1.3, 1.6, 1.9 (crater_completeness.py), F_0.3. Writes hosting_fits.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings; warnings.filterwarnings("ignore")
import os as _o, sys as _s; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__))); from s2cfg import *   # g4: inputs by S2_CAT
HERE = os.path.dirname(os.path.abspath(__file__)); pass
d = pd.read_csv(FEAT, low_memory=False, dtype={"id": str})
m = pd.read_csv(CCOMP, dtype={"id": str}); d = d.merge(m.drop(columns=["C1"]), on="id", how="left")
assert d.C1_old.notna().all()
d["lnC1"] = np.log(d.C1 + 0.01)
for a in ("1.3", "1.6", "1.7", "1.8", "1.9"): d[f"lnE_{a}"] = np.log(d[f"E_{a}"] + 0.01)
d["hot"] = hot_mask(d.lon.values).astype(int); d["dT10"] = d.dT / 10
def OR(model, t, scale=1.0):
    b, ci = model.params[t], model.conf_int().loc[t]; return [round(float(np.exp(scale * b)), 3), round(float(np.exp(scale * ci[0])), 3), round(float(np.exp(scale * ci[1])), 3)]
MEAS = ["C1", "lnC1", "C1_act", "lnE_1.6", "lnE_1.3", "lnE_1.9", "F_0.3", "lnE_1.7", "lnE_1.8"]
R = {}
for M in MEAS:
    d["C"] = d[M]; o = {}
    ctrl = "logD + abslat + C + nsp + dT + ln_local"
    m1 = smf.logit(f"h09 ~ {ctrl} + nested", d).fit(disp=0); o["nested_OR"] = OR(m1, "nested"); o["dT_OR_per10K_nestfit"] = OR(m1, "dT", 10); o["llf_nestfit"] = round(float(m1.llf), 1)
    m2 = smf.logit(f"h09 ~ {ctrl} + vent", d).fit(disp=0); o["vent_OR"] = OR(m2, "vent")
    a = d[(d.D >= 10) & d.age.isin(["fresh", "standard", "degraded"])].copy(); a["fresh"] = (a.age == "fresh").astype(int); a["standard"] = (a.age == "standard").astype(int)
    m3 = smf.logit("h09 ~ fresh + standard + logD + abslat + C", a).fit(disp=0); o["fresh_vs_degraded_OR"] = OR(m3, "fresh")
    for lo, hi in [(-15, 15), (15, 30), (30, 45), (45, 60)]:
        b = d[(d.lat >= lo) & (d.lat < hi)]; mm = smf.logit("h09 ~ hot + logD + C + abslat", b).fit(disp=0)
        o[f"hot_OR_{lo}..{hi}"] = OR(mm, "hot") + [float(mm.pvalues["hot"])]
    b = d[d.abslat < 60]; mm = smf.logit("h09 ~ dT10 + logD + C + abslat", b).fit(disp=0); o["dT_OR_per10K_lt60"] = OR(mm, "dT10"); o["llf_dT_lt60"] = round(float(mm.llf), 1)
    R[M] = o
    print(M, o, flush=True)
json.dump(R, open(f"{HERE}/out/hosting_fits_{CAT}.json", "w"), indent=1)
