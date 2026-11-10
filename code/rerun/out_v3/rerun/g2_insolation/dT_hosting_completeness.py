#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_v7/dT_hosting_completeness.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Hot-longitude (Delta T) hosting fit of analysis_v3/dT_hosting.py with imaging completeness added as a control
(29 Oct 2026, in response to review). Completeness C1 (probability that the imaging detects a 1 km hollow at the crater
center) is taken from the crater-graph node table (mercury-crater-graph, outputs/crater_dag_nodes.csv.gz), joined by
Herrick unique_id. Usage: python3 dT_hosting_completeness.py <mercury-hollows repo root> <mercury-crater-graph repo root>"""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_GRLDATA = f"{G2_ROOT}/manuscripts/GRL_insolation_hollows/data"
# --- end g2 block ---
import sys, os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
REPO, MCG = sys.argv[1], sys.argv[2]
G = f"{REPO}/manuscripts/GRL_insolation_hollows"; sys.path.insert(0, f"{G}/scripts"); import thermal
z = np.load(f"{G}/data/thermal_grid.npz"); la, lo, T = z["la"], z["lo"], z["G"]
s = pd.read_csv(f"{IN}/test2_craters_h18.csv"); s = s[(s.diameter >= 5) & (s.abslat <= 60)].copy()
s["dT"] = thermal.lookup(la, lo, T, s.latitude.values, s.longitude.values) - thermal.lookup(la, lo, T, s.latitude.values, np.full(len(s), 90.0))
n = pd.read_csv(f"{G2_ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_v9/crater_dag_nodes.csv", dtype={"id": str}, low_memory=False)
s["id"] = s.unique_id.astype(int).astype(str); s = s.merge(n[["id", "C1"]], on="id", how="left")
R = {"n": int(len(s)), "hosts": int(s.hosts.sum()), "C1_missing": int(s.C1.isna().sum())}
s = s.dropna(subset=["C1"])
def summarize(m, name):
    b, bi = m.params["dT"], m.params["dT:logD"]; V = m.cov_params().loc[["dT", "dT:logD"], ["dT", "dT:logD"]].values
    out = dict(aic=float(m.aic), p_interaction=float(m.pvalues["dT:logD"]), zero_km=float(10 ** (-b / bi)))
    for D in (10, 100):
        x = np.array([1, np.log10(D)]); c = b + bi * np.log10(D); se = float(np.sqrt(x @ V @ x))
        out[f"OR10K_{D}km"] = [float(np.exp(10 * c)), float(np.exp(10 * (c - 1.96 * se))), float(np.exp(10 * (c + 1.96 * se)))]
    if "C1" in m.params: out["C1"] = [float(m.params["C1"]), float(m.bse["C1"]), float(m.pvalues["C1"])]
    R[name] = out
summarize(smf.logit("hosts ~ dT * logD + abslat", data=s).fit(disp=0), "without_C1")
summarize(smf.logit("hosts ~ dT * logD + abslat + C1", data=s).fit(disp=0), "with_C1")
json.dump(R, open(f"{G2_OUT}/dT_hosting_completeness_record.json", "w"), indent=1)
print(json.dumps(R, indent=1))
