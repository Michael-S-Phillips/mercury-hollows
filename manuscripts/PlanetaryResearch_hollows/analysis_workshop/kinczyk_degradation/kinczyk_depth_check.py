#!/usr/bin/env python3
"""Depth by host class: Kinczyk vs Herrick on the same hollows (hosts >= 40 km with a Kinczyk class 1-5), clustered by
image and by host; with and without Herrick-filled hosts (5 Nov 2026). Writes kinczyk_depth_check.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
from scipy.spatial import cKDTree
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../..")
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
m = pd.read_csv(f"{HERE}/kinczyk_matched.csv", low_memory=False, dtype={"id": str})
k = m[m.id.notna() & m.kclass.between(1, 5)].copy(); k["kc"] = k.kclass.astype(int); kk = k.set_index("id")
R_M = 2439.4
d = pd.read_csv(f"{ROOT}/manuscripts/insolation_thermal/data/test1_rows_inc.csv")
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
t = cKDTree(xyz(f.lat.values, f.lon.values % 360)); P = xyz(d.latitude.values, d.longitude.values % 360); host = np.full(len(d), -1)
for i, cl in enumerate(t.query_ball_point(P, r=0.36)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(P[i] @ t.data[cl].T, -1, 1)) * R_M; ok = ang <= 1.1 * f.D.values[cl] / 2
    if ok.any(): host[i] = cl[np.argmin(np.where(ok, f.D.values[cl], np.inf))]
d = d[host >= 0].copy(); d["host"] = f.id.values[host[host >= 0]]; d["ltan"] = np.log(np.tan(np.radians(d.inc)))
d = d[d.host.isin(kk.index)].copy(); d["kc"] = kk.loc[d.host, "kc"].values; d["age"] = kk.loc[d.host, "age"].values
d["logD"] = kk.loc[d.host, "logD"].values
d["kgrp"] = pd.Categorical(np.where(d.kc <= 2, "C1-2", np.where(d.kc == 3, "C3", "C4-5")), categories=["C1-2", "C3", "C4-5"])
d["hcat"] = pd.Categorical(d.age, categories=["filled", "degraded", "standard", "fresh"])
R = {"crosstab": pd.crosstab(d.kc, d.age).to_dict(), "hosts_crosstab": pd.crosstab(d.drop_duplicates("host").kc, d.drop_duplicates("host").age).to_dict()}
def ratios(model, terms):
    o = {}
    for tt in terms:
        b, ci = model.params[tt], model.conf_int().loc[tt]; o[tt] = [round(float(np.exp(b)), 2), round(float(np.exp(ci[0])), 2), round(float(np.exp(ci[1])), 2)]
    return o
for cl in ["filename", "host"]:
    g = pd.factorize(d[cl])[0]
    for lab, sub in [("all", d), ("no_filled", d[d.age != "filled"])]:
        gg = pd.factorize(sub[cl])[0]
        mk = smf.ols("lnh ~ C(kgrp) + lnpx + ltan + abslat", sub).fit(cov_type="cluster", cov_kwds={"groups": gg})
        R[f"kinczyk_{lab}_cl_{cl}"] = ratios(mk, ["C(kgrp)[T.C3]", "C(kgrp)[T.C4-5]"])
        mk2 = smf.ols("lnh ~ C(kgrp) + lnpx + ltan + abslat + logD", sub).fit(cov_type="cluster", cov_kwds={"groups": gg})
        R[f"kinczyk_{lab}_cl_{cl}_plus_logD"] = ratios(mk2, ["C(kgrp)[T.C3]", "C(kgrp)[T.C4-5]"])
    mh = smf.ols("lnh ~ C(hcat) + lnpx + ltan + abslat", d).fit(cov_type="cluster", cov_kwds={"groups": g})
    R[f"herrick_all_cl_{cl}"] = ratios(mh, ["C(hcat)[T.degraded]", "C(hcat)[T.standard]", "C(hcat)[T.fresh]"])
R["n"] = {"hollows": int(len(d)), "hosts": int(d.host.nunique()), "hollows_no_filled": int((d.age != "filled").sum()), "hosts_no_filled": int(d[d.age != "filled"].host.nunique())}
json.dump(R, open(f"{HERE}/kinczyk_depth_check.json", "w"), indent=1); print(json.dumps(R, indent=1))
