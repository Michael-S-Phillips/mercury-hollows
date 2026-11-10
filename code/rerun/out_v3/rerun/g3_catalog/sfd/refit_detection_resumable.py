#!/usr/bin/env python3
"""g3_catalog (6 Nov 2026): the detection-model refit of refit_detection.py (same likelihood, bounds and Nelder-Mead options), run
in pieces because a device shell call is limited to ~180 s and background jobs do not survive it. Each call advances one start by
at most MAXITER iterations and saves the simplex in sfd/out_{cat}_refit/refit_state.json; rerun until "converged". Starts: the
v2 optimum (theta_out_v2.json) and the original grid start that reached it (p0 0.2, x0 25). Writes sfd/out_{cat}_refit/theta_out_v2.json
(same keys as the original) when every start has converged. New script. Env G3_CAT; arg MAXITER (default 400)."""
import os, sys, json, numpy as np, pandas as pd
from scipy.special import erf
from scipy.optimize import minimize
G3 = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog"); CAT = os.environ.get("G3_CAT", "v3")
O = f"{G3}/sfd/out_{CAT}_refit"; os.makedirs(O, exist_ok=True); ST = f"{O}/refit_state.json"
MAXITER = int(sys.argv[1]) if len(sys.argv) > 1 else 400
CV = pd.read_csv(f"{G3}/inputs/completeness_pairs_{CAT}.csv").sort_values("h", kind="stable").reset_index(drop=True)
Phi = lambda z: 0.5 * (1 + erf(z / np.sqrt(2)))
def pdet(x, inc_v, th):
    pmax, lx50, ls1, lxhi, ls2, binc = th
    p = pmax * Phi((np.log(x) - lx50) / np.exp(ls1)) * (1 - Phi((np.log(x) - lxhi) / np.exp(ls2)))
    return np.clip(p * np.exp(binc * (inc_v - 60) / 25), 1e-6, 0.999)
hid = CV.h.values; X = CV.D_px.values; I = CV.inc.values; Y = CV.detected.values.astype(float)
uh, starts = np.unique(hid, return_index=True)
def nll(th):
    if not (0.005 < th[0] < 0.95) or th[1] > th[3]: return 1e18
    p = pdet(X, I, th); S1 = np.add.reduceat(Y * np.log(p) + (1 - Y) * np.log(1 - p), starts); S0 = np.add.reduceat(np.log(1 - p), starts)
    return -np.sum(S1 - np.log(1 - np.exp(S0) + 1e-300))
v2 = json.load(open(os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json")))
S = json.load(open(ST)) if os.path.exists(ST) else {"starts": {"v2_optimum": {"x0": v2["theta_raw"], "simplex": None, "nit": 0, "converged": False},
                                                                 "grid_p0.2_x25": {"x0": [0.2, float(np.log(25)), float(np.log(0.5)), float(np.log(170)), float(np.log(0.35)), 0.0], "simplex": None, "nit": 0, "converged": False}}}
for k, s in S["starts"].items():
    if s["converged"]: continue
    opt = dict(maxiter=MAXITER, xatol=1e-5, fatol=1e-3)
    if s["simplex"] is not None: opt["initial_simplex"] = np.array(s["simplex"])
    r = minimize(nll, s["x0"] if s["simplex"] is None else np.array(s["simplex"])[0], method="Nelder-Mead", options=opt)
    s.update(simplex=r.final_simplex[0].tolist(), x=r.x.tolist(), fun=float(r.fun), nit=s["nit"] + int(r.nit), converged=bool(r.success))
    print(k, "nit", s["nit"], "fun", round(r.fun, 3), "converged", r.success, np.round(r.x, 4)); json.dump(S, open(ST, "w"), indent=1)
    break
if all(s["converged"] for s in S["starts"].values()):
    b = min(S["starts"].values(), key=lambda s: s["fun"]); th = b["x"]
    json.dump({"p_max": float(th[0]), "x50_px": float(np.exp(th[1])), "sigma_lo": float(np.exp(th[2])), "x_hi_px": float(np.exp(th[3])), "sigma_hi": float(np.exp(th[4])), "b_inc": float(th[5]), "n_pairs": int(len(CV)), "theta_raw": [float(v) for v in th]}, open(f"{O}/theta_out_v2.json", "w"), indent=1)
    print("all converged; theta written", np.round(th, 3), "nll v2-theta at this catalog:", round(nll(v2["theta_raw"]), 3), "best:", round(b["fun"], 3))
