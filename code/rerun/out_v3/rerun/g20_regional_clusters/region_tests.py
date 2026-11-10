#!/usr/bin/env python3
"""g20 (9 Nov 2026): tests of two explanations for the excess regions. (A) Detection of sub-km hollows: refit the joint model
(thinned presence, cloglog, offset 2 ln D, seven groups) with completeness for 0.3 km or 0.5 km hollows, or both 0.3 and 1 km, in
place of 1 km, and recompute observed / predicted hosts within 500 km of the three peaks. (B) A threshold response to insolation:
observed / predicted hosts by Delta T_w bin within 30 deg of the equator. Writes region_tests.json. Not for the paper."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.api as sm
HE = os.path.expanduser("~/mnt/hollows-emissivity"); RR = f"{HE}/code/rerun/out_v3/rerun"; SM = f"{RR}/g4_completeness/g1/summary_model"; HERE = os.path.dirname(os.path.abspath(__file__)); RM = 2439.7
os.environ["SM_DIR"] = f"{SM}/v3r"; sys.path.insert(0, SM); import common
common.DATA = f"{SM}/v3r/crater_features_xrs.csv"; d = common.load(); d["lnA"] = 2 * np.log(d.D)
G = np.load(f"{RR}/g4_completeness/completeness_10ppd_v3.npz"); Cv = G["C"]; Dref = list(G["Dref"])
iy = np.clip(((90 - d.lat.values) * 10).astype(int), 0, 1799); ix = ((d.lon.values % 360) * 10).astype(int) % 3600
for s in (0.3, 0.5, 1.0): d[f"lnC{s}"] = np.log(Cv[Dref.index(s)][iy, ix] + 0.01)
d["lnC_feat"] = np.log(d.C1 + 0.01)
cols = [c for g in common.GROUPS if g != "detection" for c in common.GROUPS[g]]
fam = sm.families.Binomial(link=sm.families.links.CLogLog()); y = d.h09.values.astype(float)
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
X3 = xyz(d.lat.values, d.lon.values % 360)
peaks = {"1_province": (65, 125), "2_abu_nuwas_sinan": (13, 329), "3_tyagaraja": (7, 203)}
win = {k: (X3 @ xyz(np.array([a]), np.array([b]))[0]) >= np.cos(500 / RM) for k, (a, b) in peaks.items()}
out = {"A_detection_size": {}}
for nm, det in (("C1_features (paper)", ["lnC_feat"]), ("C1_grid", ["lnC1.0"]), ("C0.5_grid", ["lnC0.5"]), ("C0.3_grid", ["lnC0.3"]), ("C0.3_and_C1", ["lnC0.3", "lnC1.0"])):
    X = sm.add_constant(d[det + cols].astype(float)); m = sm.GLM(y, X, family=fam, offset=d.lnA.values).fit(maxiter=200)
    P = m.predict(X, offset=d.lnA.values); r = dict(loglik=round(float(m.llf), 1), sumP=round(float(P.sum()), 1), g={k: round(float(m.params[k]), 3) for k in det})
    for k, w in win.items():
        O, E, V = y[w].sum(), P[w].sum(), (P[w] * (1 - P[w])).sum(); r[k] = dict(O=int(O), E=round(float(E), 1), z=round(float((O - E) / np.sqrt(V)), 1))
    out["A_detection_size"][nm] = r; print(nm, json.dumps(r))
X = sm.add_constant(d[["lnC_feat"] + cols].astype(float)); P = sm.GLM(y, X, family=fam, offset=d.lnA.values).fit(maxiter=200).predict(X, offset=d.lnA.values)
d["P"] = P; t = d[d.abslat <= 30]; B = []
for a, b in ((0, 20), (20, 35), (35, 45), (45, 50), (50, 55), (55, 60), (60, 63)):
    u = t[(t.dT >= a) & (t.dT < b)]; B.append(dict(dT=[a, b], n=int(len(u)), hosts=int(u.h09.sum()), predicted=round(float(u.P.sum()), 1)))
out["B_dT_bins_lat_le30"] = B; print(json.dumps(B))
for k in peaks: print(k, "median dT in window", round(float(d[win[k]].dT.median()), 1))
json.dump(out, open(f"{HERE}/region_tests.json", "w"), indent=1)
