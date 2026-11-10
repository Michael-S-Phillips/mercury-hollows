#!/usr/bin/env python3
"""Check of the v2 per-frame detection function and completeness grid (for Figure S-completeness, v8).
The fitted parameters of the v2 detection model (92,382 hollow-image pairs, code/rerun/out_v2/completeness_pairs.csv) were
recorded only to 3 decimals (code/rerun/dedup_before_after.json, out_v2/log.txt); hollow-networks/results/
per_frame_detection_function.json is the stale v1 fit (147,704 pairs) and is not used. This script refits the model with the
likelihood of code/rerun/build_hollow_products.py (binomial over hollow-image pairs, conditional on >= 1 detection per hollow;
Nelder-Mead), starting from the recorded v2 values and from the first start of the original grid, and checks that the stored
completeness grid for 1 km hollows (out_v2/completeness_10ppd.npz) is reproduced from out_v2/coverage_10ppd.npz with the
refitted parameters. Writes detection_function_v2_check_record.json next to this script.
Usage: python3 detection_function_v2_check.py ~/mnt/hollows-emissivity"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.special import erf
from scipy.optimize import minimize
ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/mnt/hollows-emissivity")
OUT = os.path.dirname(os.path.abspath(__file__))
V2 = f"{ROOT}/code/rerun/out_v2"
def Phi(z): return 0.5 * (1 + erf(z / np.sqrt(2)))
def pdet(x, inc_v, th):
    pmax, lx50, ls1, lxhi, ls2, binc = th
    p = pmax * Phi((np.log(x) - lx50) / np.exp(ls1)) * (1 - Phi((np.log(x) - lxhi) / np.exp(ls2)))
    return np.clip(p * np.exp(binc * (inc_v - 60) / 25), 1e-6, 0.999)
CV = pd.read_csv(f"{V2}/completeness_pairs.csv")
hid = CV.h.values; X = CV.D_px.values; I = CV.inc.values; Y = CV.detected.values.astype(float)
uh, starts = np.unique(hid, return_index=True)
def nll(th):
    if not (0.005 < th[0] < 0.95) or th[1] > th[3]: return 1e18
    p = pdet(X, I, th)
    S1 = np.add.reduceat(Y * np.log(p) + (1 - Y) * np.log(1 - p), starts); S0 = np.add.reduceat(np.log(1 - p), starts)
    return -np.sum(S1 - np.log(1 - np.exp(S0) + 1e-300))
rec = json.load(open(f"{ROOT}/code/rerun/dedup_before_after.json"))["detection_function_theta"]["v2"]
fits = []
for name, x0 in [("from_recorded_v2", rec), ("from_original_start_1", [0.05, np.log(12), np.log(0.5), np.log(170), np.log(0.35), 0.0])]:
    r = minimize(nll, x0, method="Nelder-Mead", options=dict(maxiter=10000, xatol=1e-5, fatol=1e-3))
    fits.append(dict(start=name, nll=float(r.fun), nit=int(r.nit), theta_log=[float(v) for v in r.x]))
best = min(fits, key=lambda f: f["nll"]); th = np.array(best["theta_log"])
# completeness grid for 1 km hollows, as build_hollow_products.py builds it (band pixel scales, incidence 60 deg)
cov = np.load(f"{V2}/coverage_10ppd.npz"); BANDS = cov["bands"]; nband = cov["nband"]
bpx = np.array([np.sqrt(max(a, 4) * b) for a, b in BANDS])
Cz = np.load(f"{V2}/completeness_10ppd.npz"); Dref = list(Cz["Dref"]); C1 = Cz["C"][Dref.index(1.0)]
Cr = np.ones_like(C1)
for bi in range(len(BANDS)):
    Cr *= (1 - pdet(1.0 * 1000 / bpx[bi], 60.0, th)) ** nband[bi].astype(np.float32)
Cr = 1 - Cr
lat = 90 - (np.arange(C1.shape[0]) + 0.5) / 10; w = np.cos(np.radians(lat))[:, None] * np.ones_like(C1)
R = dict(n_pairs=int(len(CV)), n_hollows=int(len(uh)), n_detected_pairs=int(Y.sum()),
         form="p(x,i) = p_max Phi((ln x - ln x50)/sigma_lo) [1 - Phi((ln x - ln x_hi)/sigma_hi)] exp(b_inc (i - 60)/25), x = D/pixel scale, clipped to [1e-6, 0.999]",
         recorded_v2_theta_log=rec, nll_at_recorded=float(nll(np.array(rec))), fits=fits,
         theta=dict(p_max=float(th[0]), x50_px=float(np.exp(th[1])), sigma_lo=float(np.exp(th[2])), x_hi_px=float(np.exp(th[3])),
                    sigma_hi=float(np.exp(th[4])), b_inc=float(th[5])),
         data_ranges=dict(D_px_p01_p99=[float(np.percentile(X, 1)), float(np.percentile(X, 99))], D_px_min_max=[float(X.min()), float(X.max())],
                          inc_p01_p99=[float(np.percentile(I, 1)), float(np.percentile(I, 99))], inc_min_max=[float(I.min()), float(I.max())],
                          detected_D_px_p01_p99=[float(np.percentile(X[Y == 1], 1)), float(np.percentile(X[Y == 1], 99))]),
         grid_check_C1km=dict(max_abs_diff=float(np.nanmax(np.abs(Cr - C1))), bands_m=BANDS.tolist(), band_px_m=bpx.tolist(), incidence_deg=60.0),
         C1km_summary=dict(area_frac_C_ge_0_5=float((w * (C1 >= 0.5)).sum() / w.sum()), area_frac_C_ge_0_2=float((w * (C1 >= 0.2)).sum() / w.sum()),
                           area_frac_C_gt_0=float((w * (C1 > 0)).sum() / w.sum()), max=float(np.nanmax(C1)),
                           area_frac_N_ge_0_5=float((w * (C1 >= 0.5))[lat > 0].sum() / w[lat > 0].sum()),
                           area_frac_S_ge_0_5=float((w * (C1 >= 0.5))[lat < 0].sum() / w[lat < 0].sum())))
json.dump(R, open(f"{OUT}/detection_function_v2_check_record.json", "w"), indent=1)
print(json.dumps(R, indent=1))
