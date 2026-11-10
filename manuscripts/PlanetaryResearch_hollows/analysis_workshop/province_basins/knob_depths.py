#!/usr/bin/env python3
"""Hollows on circum-Caloris knobs (Wright et al. 2020 supplementary shapefile, 2,011 knobs) (5 Nov 2026).
A hollow is on a knob if it lies within the knob's half long axis (or 3.5 km for unmeasured certain/uncertain knobs,
2.0 km for knobs < 5 km) plus 0.5 km of the knob point. (1) Cataloged hollows (hollows_canonical.csv, non-trace) on
knobs, and whether they lie in host craters (host_assignments.csv). (2) Depths (test1_rows_inc.csv, 2,518 measured
hollows): ln depth on knob indicator, ln pixel scale, ln tan incidence, and |latitude|, standard errors clustered by
image, against all other measured hollows and against measured hollows in host craters only. Writes knob_depths.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
sys.path.insert(0, f"{HERE}/wright2020_knobs"); from read_dbf import read_dbf
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
k = read_dbf(f"{HERE}/wright2020_knobs/CircumCalorisKnobs_SupplementaryShapefile.dbf"); k["lon"] = k.Longitude % 360
k["r_km"] = np.where(k.LongAxis_k > 0, k.LongAxis_k / 2, np.where(k.KnobType == "Knob (<5 km)", 2.0, 3.5)) + 0.5
KP = xyz(k.Latitude.values, k.lon.values); kt = cKDTree(KP)
def on_knob(la, lo):
    P = xyz(la, lo); out = np.full(len(P), -1)
    for i, idx in enumerate(kt.query_ball_point(P, r=2 * np.sin(11.0 / R_M / 2))):
        if not idx: continue
        d = 2 * np.arcsin(np.clip(np.linalg.norm(KP[idx] - P[i], axis=1) / 2, 0, 1)) * R_M
        ok = [j for j, dd in zip(idx, d) if dd <= k.r_km.values[j]]
        if ok: out[i] = ok[int(np.argmin([d[idx.index(j)] for j in ok]))]
    return out
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)].copy(); h["lon360"] = h.lon % 360
ha = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", usecols=["hollow_id"])
h["knob"] = on_knob(h.lat.values, h.lon360.values); h["in_host"] = h.hollow_id.isin(ha.hollow_id)
hk = h[h.knob >= 0]
R = {"knobs": int(len(k)), "knob_types": k.KnobType.value_counts().to_dict(),
     "cataloged_hollows_on_knobs": int(len(hk)), "knobs_with_hollows": int(hk.knob.nunique()),
     "knobs_with_hollows_by_type": k.loc[hk.knob.unique(), "KnobType"].value_counts().to_dict(),
     "knob_hollows_in_host_craters": int(hk.in_host.sum())}
d = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv"); d["lon360"] = d.longitude % 360
d["knob_i"] = on_knob(d.latitude.values, d.lon360.values); d["knob"] = (d.knob_i >= 0).astype(int)
# host status of measured hollows: nearest cataloged hollow within 1 km
ct = cKDTree(xyz(h.lat.values, h.lon360.values)); dist, j = ct.query(xyz(d.latitude.values, d.lon360.values))
d["in_host"] = np.where(2 * np.arcsin(np.clip(dist / 2, 0, 1)) * R_M < 1.0, h.in_host.values[j], np.nan)
d["lntan"] = np.log(np.tan(np.radians(d.inc)))
R["measured_on_knobs"] = int(d.knob.sum()); R["measured_knobs"] = int(d.loc[d.knob == 1, "knob_i"].nunique())
R["measured_on_knobs_in_host"] = int((d.knob.eq(1) & d.in_host.eq(1)).sum())
R["measured_on_knobs_images"] = int(d.loc[d.knob == 1, "filename"].nunique())
def fit(sub, label):
    m = smf.ols("lnh ~ knob + lnpx + lntan + abslat", sub).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(sub.filename)[0]})
    c = m.params["knob"]; lo_, hi_ = m.conf_int().loc["knob"]
    R[label] = {"n": int(len(sub)), "n_knob": int(sub.knob.sum()), "depth_ratio_knob": [round(float(np.exp(c)), 2), round(float(np.exp(lo_)), 2), round(float(np.exp(hi_)), 2)], "p": round(float(m.pvalues["knob"]), 3)}
fit(d, "vs_all_other_measured")
fit(d[(d.knob == 1) | (d.in_host == 1)], "vs_measured_in_host_craters")
dk = d[d.knob == 1]
R["knob_depth_m_median_raw"] = round(float(dk.depth.median()), 1); R["other_depth_m_median_raw"] = round(float(d[d.knob == 0].depth.median()), 1)
R["knob_depth_corr_median_m"] = round(float(np.exp(dk.lnh_corr.median())), 1); R["other_depth_corr_median_m"] = round(float(np.exp(d[d.knob == 0].lnh_corr.median())), 1)
R["knob_measured_lat_range"] = [round(float(dk.latitude.min()), 1), round(float(dk.latitude.max()), 1)]
json.dump(R, open(f"{HERE}/knob_depths.json", "w"), indent=1); print(json.dumps(R, indent=1))
d[d.knob == 1].to_csv(f"{HERE}/knob_depths_measured.csv", index=False)
