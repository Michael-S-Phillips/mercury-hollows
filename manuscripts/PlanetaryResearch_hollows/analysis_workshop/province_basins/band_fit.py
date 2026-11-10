#!/usr/bin/env python3
"""The NW-SE band of hollows between b59 and Caloris (5 Nov 2026). Hollows (hollows_canonical.csv, non-trace) at 45-60N,
121-142E, east of b59's interior cluster. Fit a great circle (smallest-eigenvalue normal of the hollow unit vectors,
robust: refit on hollows within 60 km of the first fit), then report the band's length and width, its azimuth at the
midpoint against the azimuth to the Caloris center, the distance of the great circle from the Caloris center, and how
many band hollows lie inside host craters (host_assignments.csv). Writes band_fit.json."""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
def ll(v): v = v / np.linalg.norm(v); return float(np.degrees(np.arcsin(v[2]))), float(np.degrees(np.arctan2(v[1], v[0])) % 360)
def bearing(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    return float(np.degrees(np.arctan2(np.sin(lo2-lo1)*np.cos(la2), np.cos(la1)*np.sin(la2)-np.sin(la1)*np.cos(la2)*np.cos(lo2-lo1))) % 360)
h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", low_memory=False); h = h[~h.is_crater_trace.astype(bool)].copy(); h["lon360"] = h.lon % 360
ha = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", usecols=["hollow_id"])
sel = h[h.lat.between(45, 60) & h.lon360.between(121, 142)].copy()
P = xyz(sel.lat.values, sel.lon360.values)
def fit(P):
    w, v = np.linalg.eigh(P.T @ P); return v[:, 0]                    # normal of the best-fit great circle
n = fit(P); d = np.abs(np.degrees(np.arcsin(P @ n))) * np.pi / 180 * R_M
for _ in range(3):
    keep = d < 60; n = fit(P[keep]); d = np.abs(P @ n) * R_M
band = sel[d < 40].copy(); B = xyz(band.lat.values, band.lon360.values)
mid = B.mean(0); mla, mlo = ll(mid)
# direction of the great circle at the midpoint: tangent = n x mid
t = np.cross(n, mid / np.linalg.norm(mid)); tla, tlo = ll(mid / np.linalg.norm(mid) + 0.01 * t)
az_band = bearing(mla, mlo, tla, tlo) % 180; az_cal = bearing(mla, mlo, 29.91, 161.43) % 180
cal = xyz([29.91], [161.43])[0]; dist_cal_gc = abs(float(np.degrees(np.arcsin(cal @ n)))) * np.pi / 180 * R_M
along = np.degrees(np.arctan2(B @ np.cross(n, mid / np.linalg.norm(mid)), B @ (mid / np.linalg.norm(mid)))) * np.pi / 180 * R_M
R = {"hollows_in_box": int(len(sel)), "band_hollows_within_40km": int(len(band)), "band_fields": int(band.field_id.nunique()),
     "band_length_km_5_95": round(float(np.percentile(along, 95) - np.percentile(along, 5))), "band_halfwidth_km_median": round(float(np.median(np.abs(B @ n) * R_M)), 1),
     "midpoint": [round(mla, 2), round(mlo, 2)], "band_azimuth_deg": round(az_band, 1), "azimuth_to_Caloris_center_deg": round(az_cal, 1),
     "angle_band_vs_radial_deg": round(min(abs(az_band - az_cal), 180 - abs(az_band - az_cal)), 1),
     "great_circle_distance_from_Caloris_center_km": round(dist_cal_gc), "Caloris_radius_km": 778,
     "band_hollows_in_host_craters": int(band.hollow_id.isin(ha.hollow_id).sum())}
json.dump(R, open(f"{HERE}/band_fit.json", "w"), indent=1); print(json.dumps(R, indent=1))
