#!/usr/bin/env python3
"""g22 (9 Nov 2026), step 1: the hollows that make craters within 500 km of 65N 125E (region 1 of g20) hosts, for Michael to screen
as hollow or not. All hollows within 0.9 radii of the 23 hosts (h09 definition). Each row carries its representative detection
(finest image) and every other detection of the same hollow. Writes province_screen_set.csv and province_hosts.csv."""
import os, sys, numpy as np, pandas as pd
HE = os.path.expanduser("~/mnt/hollows-emissivity"); RR = f"{HE}/code/rerun/out_v3/rerun"; HERE = os.path.dirname(os.path.abspath(__file__)); RM = 2439.7
sys.path.insert(0, f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot"); import dep_common as C
cr = pd.read_csv(f"{RR}/g20_regional_clusters/craters_P.csv", dtype={"id": str})
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
w = (xyz(cr.lat.values, cr.lon.values % 360) @ xyz(np.array([65.0]), np.array([125.0]))[0]) >= np.cos(500 / RM)
hosts = cr[w & (cr.h09 == 1)].copy()
H = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool)]
ha = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/host_assignments.csv")[["hollow_id", "host_uid", "rho"]]
H = H.merge(ha, on="hollow_id"); H["host_uid"] = H.host_uid.astype("Int64").astype(str); H = H[H.host_uid.isin(hosts.id) & (H.rho <= 0.9)]
g = H.groupby("host_uid").agg(n_scr=("hollow_id", "size"), maxD=("D_eq_km", "max"))
hosts = hosts.merge(g, left_on="id", right_index=True); hosts["only_subkm"] = hosts.maxD < 1
b = C.load_bickel(with_corners=False)[["thumb", "num", "datetime", "x0", "x1", "y0", "y1", "px_rel", "incidence_rel"]].drop_duplicates("thumb")
s = H.merge(b, on="thumb", how="left")
s = s.merge(hosts[["id", "name", "D", "lat", "lon", "n_scr", "only_subkm", "P"]].rename(columns={"id": "host_uid", "name": "host_name", "D": "host_D", "lat": "host_lat", "lon": "host_lon", "P": "host_P"}), on="host_uid")
s = s.sort_values(["only_subkm", "host_D", "D_eq_km"], ascending=[False, True, True]).reset_index(drop=True)
s.to_csv(f"{HERE}/province_screen_set.csv", index=False); hosts.to_csv(f"{HERE}/province_hosts.csv", index=False)
print("hosts", len(hosts), "only sub-km", int(hosts.only_subkm.sum()), "| hollows", len(s), "in sub-km-only hosts", int(s.only_subkm.sum()), "| images", s.num.nunique(), "| missing bbox", int(s.x0.isna().sum()))
print(hosts[["id", "name", "D", "lat", "lon", "n09", "n_scr", "maxD", "only_subkm", "P"]].sort_values("D").to_string(index=False))
print(s.D_eq_km.describe().round(2).to_dict(), s.px_best.describe().round(1).to_dict(), "n_det", s.n_det.value_counts().to_dict())
