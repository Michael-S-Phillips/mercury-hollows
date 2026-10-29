#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
"""Assign corrected canonical (non-trace) hollows to Herrick host craters.

Host crater = the SMALLEST Herrick crater whose 1.1 * (crater radius) contains
the hollow (great-circle distance <= 1.1 * D/2).  Computes:
  rho     = great-circle distance / crater radius (km/km)
  az      = azimuth crater-center -> hollow, 0=north, clockwise (deg)
The wall a hollow sits on FACES azimuth az+180 (the facing aspect).

Writes host_assignments.csv in this folder.
"""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

R_M = 2439.7  # km, Mercury radius

BASE = _REPO + ""
CANON = f"{BASE}/hollow-networks/data_v2/hollows_canonical.csv"
HERR = [f"{BASE}/data/herrick-craters/2017je005516r_5-10_kilometer_diameter_craters.xlsx",
        f"{BASE}/data/herrick-craters/2017je005516r_craters_over_10_km_in_diameter.xlsx"]

def unit_vec(lat_deg, lon_deg):
    la = np.radians(lat_deg); lo = np.radians(lon_deg)
    return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])

def main():
    h = pd.read_csv(CANON)
    h = h[~h["is_crater_trace"]].reset_index(drop=True)
    print("non-trace hollows:", len(h))

    cr = pd.concat([pd.read_excel(f)[["unique_id","lat","lon","diam"]] for f in HERR],
                   ignore_index=True)
    cr["lon"] = cr["lon"] % 360.0
    cr = cr.reset_index(drop=True)
    print("craters:", len(cr), "diam range", cr.diam.min(), cr.diam.max())

    hv = unit_vec(h["lat"].values, h["lon"].values)
    cv = unit_vec(cr["lat"].values, cr["lon"].values)
    htree = cKDTree(hv)

    crad_km = cr["diam"].values / 2.0
    ang_rad = 1.1 * crad_km / R_M            # search half-angle per crater (rad)
    chord = 2.0 * np.sin(np.clip(ang_rad, 0, np.pi)/2.0)  # chord threshold

    # best (smallest-diam) host per hollow
    best_diam = np.full(len(h), np.inf)
    best_ci = np.full(len(h), -1, dtype=int)
    order = np.argsort(cr["diam"].values)    # ascending: first hit that contains wins
    for ci in order:
        idx = htree.query_ball_point(cv[ci], chord[ci])
        if not idx: continue
        idx = np.asarray(idx)
        take = idx[best_ci[idx] == -1]       # only hollows without a (smaller) host yet
        if take.size:
            best_ci[take] = ci
            best_diam[take] = cr["diam"].values[ci]
    print("assigned:", int((best_ci>=0).sum()), "/", len(h))

    assigned = best_ci >= 0
    hs = h[assigned].copy().reset_index(drop=True)
    ci = best_ci[assigned]
    clat = cr["lat"].values[ci]; clon = cr["lon"].values[ci]; cd = cr["diam"].values[ci]
    cuid = cr["unique_id"].values[ci]

    # spherical distance + bearing (crater center -> hollow)
    la1 = np.radians(clat); lo1 = np.radians(clon)
    la2 = np.radians(hs["lat"].values); lo2 = np.radians(hs["lon"].values)
    dlo = lo2 - lo1
    central = np.arccos(np.clip(np.sin(la1)*np.sin(la2) +
                                np.cos(la1)*np.cos(la2)*np.cos(dlo), -1, 1))
    dist_km = R_M * central
    bear = np.degrees(np.arctan2(np.sin(dlo)*np.cos(la2),
                                 np.cos(la1)*np.sin(la2) - np.sin(la1)*np.cos(la2)*np.cos(dlo)))
    az = bear % 360.0
    crad = cd/2.0
    rho = dist_km / crad

    hs["host_uid"] = cuid; hs["host_lat"] = clat; hs["host_lon"] = clon
    hs["host_D_km"] = cd; hs["dist_km"] = dist_km; hs["rho"] = rho; hs["az_deg"] = az
    hs["facing_aspect_deg"] = (az + 180.0) % 360.0
    hs.to_csv("host_assignments.csv", index=False)
    print("wrote host_assignments.csv", hs.shape)
    print("rho distribution:", np.percentile(rho,[5,25,50,75,95]).round(3))

if __name__ == "__main__":
    main()
