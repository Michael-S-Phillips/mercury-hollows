#!/usr/bin/env python3
"""Hosting and depth by host-crater degradation class of Kinczyk et al. (2020), compared with Herrick et al. (2018)
states on the same craters (5 Nov 2026). Hosting = h09 (crater_features.csv). Kinczyk classes 1 (most degraded) to 5
(fresh, rayed); code 6 (140 craters, probably the paper's 'ghost' craters) and 0 (3) are left out of the class fits.
Writes kinczyk_tests.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
from scipy.spatial import cKDTree
from scipy import stats
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../..")
R = {}
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
m = pd.read_csv(f"{HERE}/kinczyk_matched.csv", low_memory=False, dtype={"id": str})
k = m[m.id.notna() & m.kclass.between(1, 5)].copy()
k["kc"] = k.kclass.astype(int); k["kcat"] = pd.Categorical(k.kc, categories=[1, 2, 3, 4, 5])
def OR(model, t, scale=1.0):
    b, ci = model.params[t], model.conf_int().loc[t]
    return [round(float(np.exp(scale * b)), 2), round(float(np.exp(scale * ci[0])), 2), round(float(np.exp(scale * ci[1])), 2), float(model.pvalues[t])]
def fit(formula, data):
    return smf.logit(formula, data).fit(disp=0, maxiter=200)

# 0. counts
R["kinczyk_counts_hosts"] = {int(c): [int(len(g)), int(g.h09.sum()), round(float(g.h09.mean()), 3)] for c, g in k.groupby("kc")}
g6 = m[m.id.notna() & (m.kclass == 6)]; R["code6_counts_hosts"] = [int(len(g6)), int(g6.h09.sum())]

# 1. Herrick baseline in the same data (paper spec: state + logD + |lat|, craters >= 10 km), and with completeness
a = f[(f.D >= 10) & f.age.isin(["fresh", "standard", "degraded"])].copy()
a["pres"] = pd.Categorical(a.age, categories=["degraded", "standard", "fresh"])
for lab, form in [("paper_spec", "h09 ~ C(pres) + logD + abslat"), ("with_C1", "h09 ~ C(pres) + logD + abslat + C1")]:
    mm = fit(form, a); R[f"herrick_ge10_{lab}"] = {"n": int(len(a)), "hosts": int(a.h09.sum()),
        "standard_vs_degraded": OR(mm, "C(pres)[T.standard]"), "fresh_vs_degraded": OR(mm, "C(pres)[T.fresh]")}
# Herrick restricted to >= 40 km
a40 = a[a.D >= 40]
mm = fit("h09 ~ C(pres) + logD + abslat", a40)
R["herrick_ge40_paper_spec"] = {"n": int(len(a40)), "hosts": int(a40.h09.sum()), "standard_vs_degraded": OR(mm, "C(pres)[T.standard]"), "fresh_vs_degraded": OR(mm, "C(pres)[T.fresh]")}

# 2. Kinczyk classes (>= 40 km), paper spec and with completeness; ordinal trend
for lab, ctl in [("paper_spec", "logD + abslat"), ("with_C1", "logD + abslat + C1"),
                 ("full_controls", "logD + abslat + C1 + nsp + dT + ln_local")]:
    mm = fit(f"h09 ~ C(kcat) + {ctl}", k)
    R[f"kinczyk_{lab}"] = {"n": int(len(k)), "hosts": int(k.h09.sum()), **{f"class{c}_vs_1": OR(mm, f"C(kcat)[T.{c}]") for c in [2, 3, 4, 5]}, "llf": round(float(mm.llf), 1)}
    mo = fit(f"h09 ~ kc + {ctl}", k)
    R[f"kinczyk_{lab}"]["per_class_step"] = OR(mo, "kc"); R[f"kinczyk_{lab}"]["llf_ordinal"] = round(float(mo.llf), 1)
    R[f"kinczyk_{lab}"]["ordinal_vs_categorical_LR_p"] = float(stats.chi2.sf(2 * (mm.llf - mo.llf), 3))

# 3. Head-to-head on the same craters (>= 40 km, matched, Kinczyk 1-5, Herrick state incl. filled)
b = k[k.age.isin(["fresh", "standard", "degraded", "filled"])].copy()
b["hcat"] = pd.Categorical(b.age, categories=["filled", "degraded", "standard", "fresh"])
H2H = {"n": int(len(b)), "hosts": int(b.h09.sum())}
for lab, ctl in [("simple", "logD + abslat"), ("with_C1", "logD + abslat + C1"), ("full", "logD + abslat + C1 + nsp + dT + ln_local")]:
    m0 = fit(f"h09 ~ {ctl}", b); mh = fit(f"h09 ~ C(hcat) + {ctl}", b); mk = fit(f"h09 ~ C(kcat) + {ctl}", b)
    mb = fit(f"h09 ~ C(hcat) + C(kcat) + {ctl}", b)
    H2H[lab] = {"llf_none": round(float(m0.llf), 1), "llf_herrick4": round(float(mh.llf), 1), "llf_kinczyk5": round(float(mk.llf), 1), "llf_both": round(float(mb.llf), 1),
                "AIC_herrick4": round(float(mh.aic), 1), "AIC_kinczyk5": round(float(mk.aic), 1),
                "pseudoR2_gain_herrick": round(float(1 - mh.llf / m0.llf - (1 - m0.llf / m0.llf)), 4),
                "LR_kinczyk_given_herrick_p": float(stats.chi2.sf(2 * (mb.llf - mh.llf), 4)),
                "LR_herrick_given_kinczyk_p": float(stats.chi2.sf(2 * (mb.llf - mk.llf), 3)),
                "dev_explained_share_herrick": round(float((mh.llf - m0.llf) / (mb.llf - m0.llf)), 3),
                "dev_explained_share_kinczyk": round(float((mk.llf - m0.llf) / (mb.llf - m0.llf)), 3)}
# within Herrick 'filled' craters: does Kinczyk class separate hosting?
fl = b[b.age == "filled"].copy(); fl["old"] = (fl.kc == 1).astype(int)
H2H["filled_by_kclass"] = {int(c): [int(len(g)), int(g.h09.sum())] for c, g in fl.groupby("kc")}
mf = fit("h09 ~ kc + logD + abslat", fl[fl.kc <= 3]); H2H["filled_per_class_step_k1to3"] = OR(mf, "kc")
R["head_to_head_ge40"] = H2H

# 4. Delta T_w per 10 K within Kinczyk groups (|lat| <= 60) and interaction
s2 = k[k.abslat <= 60].copy(); s2["grp"] = pd.Categorical(np.where(s2.kc <= 2, "C1-2", np.where(s2.kc == 3, "C3", "C4-5")), categories=["C1-2", "C3", "C4-5"])
s2["dT10"] = s2.dT / 10
for gname, g in s2.groupby("grp", observed=True):
    mm = fit("h09 ~ dT10 + logD + abslat", g); R[f"dT_per10K_{gname}"] = OR(mm, "dT10") + [int(g.h09.sum())]
mi = fit("h09 ~ dT10 * C(grp) + logD + abslat", s2); R["dT_interaction_p"] = {kk: float(v) for kk, v in mi.pvalues.items() if "dT10:" in kk}
mp = fit("h09 ~ dT10 + C(grp) + logD + abslat", s2); R["dT_per10K_pooled_ge40"] = OR(mp, "dT10")

# 5. Depth by host class (hosts as in step5_depth.py: smallest crater >= 5 km whose 1.1 radii contain the hollow)
R_M = 2439.4
d = pd.read_csv(f"{ROOT}/manuscripts/insolation_thermal/data/test1_rows_inc.csv")
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
t = cKDTree(xyz(f.lat.values, f.lon.values % 360)); P = xyz(d.latitude.values, d.longitude.values % 360); host = np.full(len(d), -1)
for i, cl in enumerate(t.query_ball_point(P, r=0.36)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(P[i] @ t.data[cl].T, -1, 1)) * R_M; ok = ang <= 1.1 * f.D.values[cl] / 2
    if ok.any(): host[i] = cl[np.argmin(np.where(ok, f.D.values[cl], np.inf))]
d = d[host >= 0].copy(); d["host"] = f.id.values[host[host >= 0]]
d["ltan"] = np.log(np.tan(np.radians(d.inc)))
kk = k.set_index("id")
d = d[d.host.isin(kk.index)].copy(); d["kc"] = kk.loc[d.host, "kc"].values; d["age"] = kk.loc[d.host, "age"].values
d["kgrp"] = pd.Categorical(np.where(d.kc <= 2, "C1-2", np.where(d.kc == 3, "C3", "C4-5")), categories=["C1-2", "C3", "C4-5"])
md = smf.ols("lnh ~ C(kgrp) + lnpx + ltan + abslat", d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.filename)[0]})
DEP = {"n_hollows": int(len(d)), "n_hosts": int(d.host.nunique()),
       "by_class": {int(c): [int(len(g)), int(g.host.nunique())] for c, g in d.groupby("kc")}}
for c in ["C3", "C4-5"]:
    b_, ci = md.params[f"C(kgrp)[T.{c}]"], md.conf_int().loc[f"C(kgrp)[T.{c}]"]
    DEP[f"{c}_vs_C1-2"] = [round(float(np.exp(b_)), 2), round(float(np.exp(ci[0])), 2), round(float(np.exp(ci[1])), 2)]
mo = smf.ols("lnh ~ kc + lnpx + ltan + abslat", d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.filename)[0]})
b_, ci = mo.params["kc"], mo.conf_int().loc["kc"]; DEP["per_class_step"] = [round(float(np.exp(b_)), 3), round(float(np.exp(ci[0])), 3), round(float(np.exp(ci[1])), 3)]
R["depth"] = DEP
json.dump(R, open(f"{HERE}/kinczyk_tests.json", "w"), indent=1); print(json.dumps(R, indent=1))
