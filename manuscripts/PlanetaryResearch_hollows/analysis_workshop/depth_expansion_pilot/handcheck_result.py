#!/usr/bin/env python3
"""Hand-check result for the 51 clean small-host (< 20 km) pilot measurements (Claude, 2 Nov 2026, visual review of
qc/smallhost_sheet_{1,2,3}.png). Classes: good (shadow cast by the sunward wall onto a hollow floor), plausible, doubtful,
wrong (crater or scarp shadow, shadow on the anti-sunward side, or a ring-shaped crater). Accepted = good + plausible.
Compares accepted small-host depths (corrected to the Blewett median pixel scale and incidence, divided by the 1.24
calibration) with Blewett depths in hosts >= 20 km (same correction), both restricted to |lat| < 60 deg (the poleward
calibration differs). Bootstrap by image. Writes smallhost_handcheck.csv and smallhost_handcheck_record.json."""
import os, sys, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dep_common as C
cls = {}
for k in [2, 14, 15, 21, 22, 31, 32, 34, 41, 44, 50]: cls[k] = "good"
for k in [0, 8, 16, 18, 28, 29, 33, 35, 42]: cls[k] = "plausible"
for k in [1, 5, 6, 9, 10, 12, 17, 20, 24, 25, 27, 30, 37, 38, 39, 40, 45, 46, 47]: cls[k] = "doubtful"
for k in [3, 4, 7, 11, 13, 19, 23, 26, 36, 43, 48, 49]: cls[k] = "wrong"
idx = pd.read_csv(f"{HERE}/smallhost_index.csv"); assert len(idx) == 51 and set(cls) == set(range(51))
idx["handcheck"] = idx.k.map(cls); idx["accepted"] = idx.handcheck.isin(["good", "plausible"])
idx.to_csv(f"{HERE}/smallhost_handcheck.csv", index=False)
R = {"counts": idx.handcheck.value_counts().to_dict(), "median_depth_by_class_m": idx.groupby("handcheck").depth_m.median().round(1).to_dict()}
# Blewett reference, corrected as in targets.py
tr = json.load(open(f"{HERE}/targets_record.json"))["correction"]
b = pd.read_csv(os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv"))
b["hc"] = np.exp(np.log(b.depth) - tr["b_lnpx"] * (np.log(b.px) - np.log(tr["px_median_m"])) - tr["b_ltan"] * (np.log(np.tan(np.radians(b.inc))) - np.log(np.tan(np.radians(tr["inc_median_deg"])))))
hid, hD, _ = C.assign_hosts(b.latitude.values, b.longitude.values % 360)
b["hostD"] = hD
ref = b[(b.hostD >= 20) & (b.abslat < 60)]
a = idx[idx.accepted & (idx.lat.abs() < 60)]
allsm = idx[idx.lat.abs() < 60]
rng = np.random.default_rng(3)
def med_ratio(x, xi, y, yi, nb=2000):
    xg = {k: v.values for k, v in x.groupby(xi)}; yg = {k: v.values for k, v in y.groupby(yi)}
    xk, yk = list(xg), list(yg); out = []
    for _ in range(nb):
        xs = np.concatenate([xg[k] for k in rng.choice(xk, len(xk))]); ys = np.concatenate([yg[k] for k in rng.choice(yk, len(yk))])
        out.append(np.median(xs) / np.median(ys))
    return dict(ratio=float(np.median(x) / np.median(y)), lo=float(np.percentile(out, 2.5)), hi=float(np.percentile(out, 97.5)))
R["blewett_ref_hosts_ge20_lat_lt60"] = dict(n=int(len(ref)), median_m=float(ref.hc.median()))
R["accepted_small_lat_lt60"] = dict(n=int(len(a)), n_hosts=int(a.host_id.nunique()), n_images=int(a.image.nunique()), median_m=float(a.depth_corr_cal_m.median()),
                                    ratio_to_blewett_ge20=med_ratio(a.depth_corr_cal_m, a.image, ref.hc, ref.filename))
R["all51_small_lat_lt60"] = dict(n=int(len(allsm)), median_m=float(allsm.depth_corr_cal_m.median()),
                                 ratio_to_blewett_ge20=med_ratio(allsm.depth_corr_cal_m, allsm.image, ref.hc, ref.filename))
json.dump(R, open(f"{HERE}/smallhost_handcheck_record.json", "w"), indent=1); print(json.dumps(R, indent=1))
