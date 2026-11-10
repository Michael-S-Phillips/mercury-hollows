#!/usr/bin/env python3
"""g20 (9 Nov 2026), for Michael's question about causes of the three excess regions. Profiles the craters within 500 km of each
region's peak node (clusters.json) against all other craters: observed and predicted hosts by crater size; imaging (best pixel scale,
number of images, completeness for 0.3 and 1 km hollows); degradation states; nesting; vents; LRM, XRS ratios and MDIS color;
distance from the Caloris center; and the hollows of the region's hosts (size, detection pixel scale, Stage 1 share at incidence < 70).
Writes region_profiles.json. Not for the paper."""
import os, sys, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
HE = os.path.expanduser("~/mnt/hollows-emissivity"); RR = f"{HE}/code/rerun/out_v3/rerun"; HERE = os.path.dirname(os.path.abspath(__file__)); RM = 2439.7
S = json.load(open(f"{HERE}/clusters.json")); C = pd.read_csv(f"{HERE}/craters_P.csv", dtype={"id": str})
F = pd.read_csv(f"{RR}/g4_completeness/g1/summary_model/v3r/crater_features_xrs.csv", dtype={"id": str})
d = C.merge(F[["id", "age", "rayed", "nested", "parent_hosts", "vent", "lrm", "s_si", "mg_si", "ca_si", "al_si", "nsp", "in_cal"]], on="id", how="left")
cc = pd.read_csv(f"{RR}/g6_remaining/color/out_v3/crater_color.csv", dtype={"id": str})[["id", "pc1", "pc2"]]
d = d.merge(cc, on="id", how="left")
for k in ("s_si", "mg_si", "ca_si", "al_si", "pc1", "pc2", "lrm"): d[k + "_z"] = (d[k] - d[k].mean()) / d[k].std()
G = np.load(f"{RR}/g4_completeness/completeness_10ppd_v3.npz"); Cv = G["C"]; Dref = list(G["Dref"]); cov = np.load(f"{RR}/g4_completeness/_scratch/coverage_10ppd_rebuilt.npz")
iy = np.clip(((90 - d.lat.values) * 10).astype(int), 0, 1799); ix = ((d.lon.values % 360) * 10).astype(int) % 3600
d["C03"] = Cv[Dref.index(0.3)][iy, ix]; d["C1g"] = Cv[Dref.index(1.0)][iy, ix]; d["bestpx"] = cov["bestpx"][iy, ix]; d["nimg"] = cov["nimg"][iy, ix]
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
def gc(la1, lo1, la2, lo2):
    a = xyz(np.atleast_1d(la1), np.atleast_1d(lo1)); b = xyz(np.atleast_1d(la2), np.atleast_1d(lo2)); return RM * np.arccos(np.clip((a * b).sum(1), -1, 1))
v = pd.read_csv(f"{HE}/data/data/MER_Vent_Locations_Classifications.csv"); v["vlat"], v["vlon"] = v["Lon"], v["Lat"] % 360
H = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool)]
ha = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/host_assignments.csv")[["hollow_id", "host_uid", "rho"]]; H = H.merge(ha, on="hollow_id", how="left")
H["host_uid"] = H.host_uid.astype("Int64").astype(str)
tree = cKDTree(xyz(d.lat.values, d.lon.values % 360)); chord = 2 * np.sin(500 / RM / 2)
def profile(sel, label, lat0=None, lon0=None):
    s = d[sel]; h = s[s.h09 == 1]
    out = dict(label=label, n_craters=int(len(s)), hosts=int(s.h09.sum()), predicted=round(float(s.P.sum()), 1))
    for nm, a, b in (("D5_10", 5, 10), ("D10_30", 10, 30), ("D30plus", 30, 1e9)):
        t = s[(s.D >= a) & (s.D < b)]; out[f"hosts_{nm}"] = [int(t.h09.sum()), round(float(t.P.sum()), 1)]
    out.update(median_C1=round(float(s.C1.median()), 3), median_C03=round(float(np.median(s.C03)), 3), median_bestpx_m=round(float(np.median(s.bestpx[np.isfinite(s.bestpx)])), 1),
               median_nimg=float(np.median(s.nimg)), host_median_D=round(float(h.D.median()), 1) if len(h) else None, host_median_P=round(float(h.P.median()), 3) if len(h) else None)
    r10 = s[s.D >= 10]; out["age_frac_D10"] = {k: round(float((r10.age == k).mean()), 3) for k in ("fresh", "standard", "degraded", "filled")}
    out["rayed"] = int(s.rayed.sum()); out["nested_frac"] = round(float(s.nested.mean()), 3); out["vent_craters"] = int(s.vent.sum())
    out["mean_z"] = {k: (round(float(s[k + "_z"].mean()), 2) if s[k + "_z"].notna().any() else None) for k in ("lrm", "s_si", "mg_si", "ca_si", "al_si", "pc1", "pc2")}
    out["mean_z_hosts"] = {k: (round(float(h[k + "_z"].mean()), 2) if len(h) and h[k + "_z"].notna().any() else None) for k in ("lrm", "pc1", "pc2", "mg_si", "s_si")}
    hh = H[H.host_uid.isin(h.id)]
    if len(hh):
        lab = hh[(hh.inc < 70) & hh.stage.isin([1, 2])]
        out["hollows"] = dict(n=int(len(hh)), median_D_km=round(float(hh.D_eq_km.median()), 2), frac_lt_0p5km=round(float((hh.D_eq_km < 0.5).mean()), 2),
                              median_px_best=round(float(hh.px_best.median()), 1), stage1_share_inc_lt70=(round(float((lab.stage == 1).mean()), 2) if len(lab) else None), n_staged=int(len(lab)),
                              frac_in_host_wall_rim=round(float((hh.rho >= 0.55).mean()), 2))
    if lat0 is not None:
        out["dist_caloris_center_km"] = round(float(gc(lat0, lon0, 29.91, 161.43)[0])); nv = int((gc(v.vlat.values, v.vlon.values, np.full(len(v), lat0), np.full(len(v), lon0)) <= 500).sum())
        out["vents_within_500km"] = nv; out["vents_expected_by_area"] = round(len(v) * (1 - np.cos(500 / RM)) / 2, 1)
    return out
R = {"global": profile(np.ones(len(d), bool), "all craters")}
regions = [("1_province", 65, 125), ("2_abu_nuwas_sinan", 13, 329), ("3_tyagaraja", 7, 203)]
allin = np.zeros(len(d), bool)
for nm, la, lo in regions:
    idx = tree.query_ball_point(xyz(np.array([la]), np.array([lo]))[0], chord); m = np.zeros(len(d), bool); m[idx] = True; allin |= m
    R[nm] = profile(m, nm, la, lo)
R["outside_regions"] = profile(~allin, "outside the three 500 km windows")
json.dump(R, open(f"{HERE}/region_profiles.json", "w"), indent=1, default=float)
for k, r in R.items():
    print("==", k); print("  ", {x: r[x] for x in r if x not in ("label", "mean_z", "mean_z_hosts", "hollows", "age_frac_D10")})
    print("   age", r["age_frac_D10"], "| z", r["mean_z"], "| z hosts", r["mean_z_hosts"]); print("   hollows", r.get("hollows"))
