# g4 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_ref/ref_size_check.py (6 Nov 2026). Changed: input/output paths by S2_CAT, added exponents/measures 1.7 and 1.8 (see make_copies_s2.py).
#!/usr/bin/env python3
"""Does the choice of 1 km as the reference hollow size for completeness change the hosting results? (3 Nov 2026,
Michael's comment on v6.1 Section 2.1.) Completeness at the crater center is re-sampled from completeness_10ppd.npz at
reference sizes 0.3, 0.5, 1, 2, and 5 km, exactly as analysis_v9/crater_dag.py samples C1, and the paper's hosting fits are
refitted with each in place of C1: nesting (h09 ~ logD + abslat + C + nsp + dT + ln_local + nested), vents (same controls +
vent), Delta T_w (coefficient per 10 K in the nesting fit), and fresh vs degraded hosts (D >= 10 km with a state;
h09 ~ age + logD + abslat + C). Also the Spearman correlation of each C with C1 and its median. Writes ref_size_record.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings; warnings.filterwarnings("ignore")
import os as _o, sys as _s; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__))); from s2cfg import *   # g4: inputs by S2_CAT
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(FEAT, dtype={"id": str})
Cz = np.load(GRID); Dref = list(Cz["Dref"])
rr = ((90 - d.lat.values) * 10).astype(int).clip(0, 1799); cc = ((d.lon.values % 360) * 10).astype(int).clip(0, 3599)
R = {"n": int(len(d)), "hosts_h09": int(d.h09.sum()), "check_C1_max_abs_diff": None, "by_ref": {}}
for D0 in [0.3, 0.5, 1.0, 2.0, 5.0]:
    d["C"] = np.nan_to_num(Cz["C"][Dref.index(D0)][rr, cc])
    if D0 == 1.0: R["check_C1_max_abs_diff"] = float(np.abs(d.C - d.C1).max())
    out = {"median_C": float(d.C.median()), "median_C_hosts": float(d.C[d.h09 == 1].median()),
           "spearman_with_C1": float(d[["C", "C1"]].corr("spearman").iloc[0, 1])}
    ctrl = "logD + abslat + C + nsp + dT + ln_local"
    m = smf.logit(f"h09 ~ {ctrl} + nested", d).fit(disp=0)
    ci = m.conf_int().loc["nested"]; out["nested_OR"] = [float(np.exp(m.params["nested"])), float(np.exp(ci[0])), float(np.exp(ci[1]))]
    ci = m.conf_int().loc["dT"]; out["dT_OR_per10K"] = [float(np.exp(10 * m.params["dT"])), float(np.exp(10 * ci[0])), float(np.exp(10 * ci[1]))]
    out["C_coef"] = float(m.params["C"]); out["llf"] = float(m.llf)
    m2 = smf.logit(f"h09 ~ {ctrl} + vent", d).fit(disp=0)
    ci = m2.conf_int().loc["vent"]; out["vent_OR"] = [float(np.exp(m2.params["vent"])), float(np.exp(ci[0])), float(np.exp(ci[1]))]
    a = d[(d.D >= 10) & d.age.isin(["fresh", "standard", "degraded"])].copy()
    a["fresh"] = (a.age == "fresh").astype(int); a["standard"] = (a.age == "standard").astype(int)
    m3 = smf.logit("h09 ~ fresh + standard + logD + abslat + C", a).fit(disp=0)
    ci = m3.conf_int().loc["fresh"]; out["fresh_vs_degraded_OR"] = [float(np.exp(m3.params["fresh"])), float(np.exp(ci[0])), float(np.exp(ci[1]))]
    R["by_ref"][str(D0)] = out
    print(D0, {k: (round(v, 3) if isinstance(v, float) else [round(x, 2) for x in v]) for k, v in out.items()}, flush=True)
json.dump(R, open(f"{HERE}/out/ref_size_record_{CAT}.json", "w"), indent=1)
print("C1 check", R["check_C1_max_abs_diff"])
