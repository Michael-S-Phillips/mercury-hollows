# g3_catalog copy of analysis_workshop/completeness_sfd/refit_detection.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Refit the per-image detection model on code/rerun/out_v2/completeness_pairs.csv exactly as build_hollow_products.py does, because
hollow-networks/results/per_frame_detection_function.json was overwritten by a later run (147,704 pairs, p_max 0.055). The out_v2
log gives theta [0.243, 3.57, -0.527, 5.278, -1.322, -0.102]. Writes theta_out_v2.json."""
import os, json, numpy as np, pandas as pd
from scipy.special import erf
from scipy.optimize import minimize
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")
CAT = _os.environ.get("G3_CAT", "v3"); THETA = _os.environ.get("G3_THETA", "fixed")
ODIR = f"{G3}/sfd/out_{CAT}_{THETA}"; _os.makedirs(ODIR, exist_ok=True)
THETA_JSON = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json")
              if THETA == "fixed" else f"{G3}/sfd/out_{CAT}_refit/theta_out_v2.json")
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
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
best = None
for p0 in [0.05, 0.1, 0.2]:
    for x0 in [12, 18, 25]:
        r = minimize(nll, [p0, np.log(x0), np.log(0.5), np.log(170), np.log(0.35), 0.0], method='Nelder-Mead', options=dict(maxiter=10000, xatol=1e-5, fatol=1e-3))
        if best is None or r.fun < best.fun: best = r
th = best.x; print("theta", np.round(th, 3), "nll", best.fun)
json.dump({"p_max": float(th[0]), "x50_px": float(np.exp(th[1])), "sigma_lo": float(np.exp(th[2])), "x_hi_px": float(np.exp(th[3])), "sigma_hi": float(np.exp(th[4])), "b_inc": float(th[5]), "n_pairs": int(len(CV)), "theta_raw": [float(v) for v in th]}, open(f"{ODIR}/theta_out_v2.json", "w"), indent=1)
