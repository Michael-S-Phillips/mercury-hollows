#!/usr/bin/env python3
"""g3_catalog (6 Nov 2026): numbers of Sections 3.1, 4.4 and S2 that are derived from the size-distribution outputs or from the
detection pairs but have no script in analysis_workshop/completeness_sfd. New script; each definition is checked against the
v2 text. Env G3_THETA = fixed | refit selects the sfd/out_{cat}_{theta} folders. Writes sfd/sfd_derived_{theta}.json.
- Exponent range: the paper rounds the min and max of the 0.3-0.7, 0.7-3 and 1-3 km exponents to one decimal (v2: 1.4-1.6) and
  applies it with N1, the corrected count > 1 km: N(>0.3 km) power law = N1 0.3^-a, N(>5 km) = N1 5^-a, N(>0.1 km) = N1 0.1^-a,
  and N(>0.1)/N(>0.3) = 3^a (Section 4.4). Values with the unrounded range are also given.
- Pairs (inputs/completeness_pairs_{cat}.csv): share of the other images that detected a hollow (all pairs minus one detection
  per hollow; v2 4.7%); the same with one detection per hollow dropped at random (seed 1) by D_px bin 60-100 and 100-200
  (v2 15.1% and 14.0%, the "14-15%"); share of pairs in which a 0.3 km hollow spans < 10 pixels (v2 80%) and median pixel scale.
- Largest detection: for non-trace hollows with >= 2 detections in code/rerun/out_v2/detections_qc.csv, median of
  max(D_eq)/D_eq(canonical) - 1 by canonical D bin (v2: 5.1 and 6.1% at 3-5 and 5-10 km; 7.3, 10.9 and 8.7% at 0.5-1, 1-2, 2-3 km)."""
import os, json, numpy as np, pandas as pd
G3 = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog"); HE = os.path.expanduser("~/mnt/hollows-emissivity")
TH = os.environ.get("G3_THETA", "fixed")
canon = {"v2": f"{HE}/code/rerun/out_v2/hollows_canonical.csv", "v3": f"{HE}/code/rerun/out_v3/hollows_canonical_v3.csv", "v3b": f"{HE}/code/rerun/out_v3/hollows_canonical_v3b.csv"}
dq = pd.read_csv(f"{HE}/code/rerun/out_v2/detections_qc.csv"); gq = dq.groupby("hollow_id").agg(Dmax=("D_eq_km", "max"), n=("D_eq_km", "size"))
R = {}
for v in ("v2", "v3", "v3b"):
    th = "fixed" if v == "v2" else TH
    f = f"{G3}/sfd/out_{v}_{th}/sfd_numbers.json"
    if not os.path.exists(f): continue
    N = json.load(open(f)); a3 = [N["slope_0.3_0.7"][0], N["slope_0.7_3"][0], N["slope_1_3"][0]]
    N1 = N["counts"]["1.0"]["corrected"]; cat03 = N["counts"]["0.3"]["catalog"]
    r = {"exponents_0.3-0.7_0.7-3_1-3": a3, "exponent_3_12": N["slope_3_12"], "N1_corrected": N1}
    for nm, (lo, hi) in {"rounded": (round(min(a3), 1), round(max(a3), 1)), "exact": (min(a3), max(a3))}.items():
        r[nm] = dict(range=[lo, hi], ratio_0p1_over_0p3=[round(3 ** lo, 2), round(3 ** hi, 2)],
                     N_gt_0p1km=[round(N1 * 10 ** lo, -3), round(N1 * 10 ** hi, -3)],
                     N_gt_0p3km_powerlaw=[round(N1 * 0.3 ** -lo, -3), round(N1 * 0.3 ** -hi, -3)],
                     catalog_fraction_0p3=[round(cat03 / (N1 * 0.3 ** -lo), 3), round(cat03 / (N1 * 0.3 ** -hi), 3)],
                     N_gt_5km_powerlaw=[round(N1 * 5 ** -hi, -1), round(N1 * 5 ** -lo, -1)])
    r["N_gt_5km"] = dict(catalog=N["counts"]["5.0"]["catalog"], corrected=N["counts"]["5.0"]["corrected"], ci=N["counts"]["5.0"]["corrected_ci"])
    r["N_gt_0p3km"] = dict(catalog=cat03, corrected=N["counts"]["0.3"]["corrected"], ci=N["counts"]["0.3"]["corrected_ci"])
    r["median_C_at_size"] = N["median_C_at_size"]
    r["best_imaged_C03_ge_0.3"] = dict(n=N["sites_C03_ge_0.3"]["n"], n_ge03=N["sites_C03_ge_0.3"]["n_ge03"], slope_03_1=N["sites_C03_ge_0.3"]["slope_03_1"])
    F = json.load(open(f"{G3}/sfd/out_{v}_{th}/sfd_fit.json")); r["power_above_Dmin"] = {k: round(p["slope_cumulative"], 2) for k, p in F["power"].items()}
    r["theta"] = {k: F["theta"][k] for k in ("p_max", "x50_px", "sigma_lo", "x_hi_px", "sigma_hi", "b_inc", "n_pairs")}
    cv = pd.read_csv(f"{G3}/inputs/completeness_pairs_{v}.csv"); nh = cv.h.nunique()
    r["pairs_other_images_detected"] = round(float((cv.detected.sum() - nh) / (len(cv) - nh)), 4)
    o = cv.drop(cv[cv.detected].groupby("h").sample(1, random_state=1).index)
    r["pairs_detected_Dpx_60_100_100_200"] = [round(float(o[(o.D_px >= lo) & (o.D_px < hi)].detected.mean()), 4) for lo, hi in ((60, 100), (100, 200))]
    r["pairs_0p3km_lt_10px"] = round(float((300 / cv.px_m < 10).mean()), 3); r["pairs_median_px_m"] = round(float(cv.px_m.median()), 1)
    c = pd.read_csv(canon[v]); m = c[~c.is_crater_trace.astype(bool)].merge(gq, on="hollow_id"); m = m[m.n >= 2]; m["x"] = m.Dmax / m.D_eq_km - 1
    r["largest_detection_excess_median"] = {f"{lo}-{hi}": [round(float(m[(m.D_eq_km >= lo) & (m.D_eq_km < hi)].x.median()), 3), int(((m.D_eq_km >= lo) & (m.D_eq_km < hi)).sum())]
                                            for lo, hi in ((0.5, 1), (1, 2), (2, 3), (0.5, 3), (3, 5), (5, 10), (3, 10))}
    R[v] = r
json.dump(R, open(f"{G3}/sfd/sfd_derived_{TH}.json", "w"), indent=1); print(json.dumps(R, indent=1))
