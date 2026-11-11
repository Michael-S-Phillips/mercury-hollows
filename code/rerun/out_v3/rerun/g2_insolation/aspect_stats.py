#!/usr/bin/env python3
# g2 copy of code/rerun/latitude/aspect/aspect_stats.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Scarp-aspect statistics on corrected catalogue (position-as-aspect proxy).

LIMITATION (stated prominently): the 665 m DEM cannot resolve hollow-scale
slopes, so we use each hollow's POSITION on the host-crater wall as the aspect
proxy.  A hollow at azimuth A from the crater center sits on the wall that FACES
azimuth A+180.  Equator-facing wall => hollow on the POLEWARD side of the crater.

Reads host_assignments.csv (from assign_hosts.py); writes aspect_record.json,
aspect_wall_hollows.csv.
"""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_INSOLDATA = f"{G2_ROOT}/manuscripts/insolation_thermal/data"
# --- end g2 block ---
import json, numpy as np, pandas as pd
from scipy import stats

RNG = np.random.default_rng(20260911)
NB = 10000

# ---------- circular helpers ----------
def wrap180(a):  # -> (-180,180]
    return (a + 180.0) % 360.0 - 180.0

def rayleigh(ang_deg):
    a = np.radians(ang_deg); n = len(a)
    if n == 0: return dict(n=0, Rbar=np.nan, p=np.nan)
    C, S = np.cos(a).sum(), np.sin(a).sum()
    Rbar = np.hypot(C, S) / n
    Z = n * Rbar**2
    p = np.exp(-Z) * (1 + (2*Z - Z**2)/(4*n) - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4)/(288*n**2))
    return dict(n=int(n), Rbar=float(Rbar), mean_deg=float(np.degrees(np.arctan2(S, C)) % 360), p=float(min(max(p,0),1)))

def vtest(ang_deg, mu0_deg=0.0):
    """V-test: concentration toward mu0. Returns one-sided p."""
    a = np.radians(ang_deg); n = len(a)
    if n == 0: return dict(n=0, V=np.nan, p=np.nan)
    mu0 = np.radians(mu0_deg)
    Vbar = np.mean(np.cos(a - mu0))
    u = Vbar * np.sqrt(2*n)
    p = 1 - stats.norm.cdf(u)
    return dict(n=int(n), Vbar=float(Vbar), u=float(u), p=float(p))

def circ_median(ang_deg):
    """Circular median (Fisher): angle minimizing mean circular distance."""
    a = np.radians(ang_deg); n = len(a)
    if n == 0: return np.nan
    # candidate medians = data points; pick minimizing sum of circular abs dist
    diff = np.abs(wrap180(np.degrees(a)[:, None] - np.degrees(a)[None, :]))
    med = ang_deg[np.argmin(diff.sum(axis=1))]
    return float(med % 360)

def circ_median_boot(ang_deg, nb=2000, cap=250):
    ang_deg = np.asarray(ang_deg)
    if len(ang_deg) < 3: return (np.nan, np.nan)
    meds = np.empty(nb)
    for i in range(nb):
        s = RNG.choice(ang_deg, min(len(ang_deg), cap), replace=True)
        meds[i] = circ_median(s)
    return (float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5)))

# ---------- poleward / equator-facing scoring ----------
def poleward_flag(lat, az_deg):
    """True if hollow sits on the poleward side of the crater (equator-facing wall)."""
    return np.sign(lat) * np.cos(np.radians(az_deg)) > 0

def equator_dev(lat, facing_deg):
    """Deviation of facing aspect from the equator-facing direction (deg, wrapped).
    North: equator dir = 180; South: equator dir = 0/360."""
    eq_dir = np.where(lat >= 0, 180.0, 0.0)
    return wrap180(facing_deg - eq_dir)

def binom_ci(k, n):
    if n == 0: return (np.nan, np.nan, np.nan)
    p = k/n
    bt = stats.binomtest(k, n, 0.5)
    lo, hi = bt.proportion_ci(confidence_level=0.95)
    return p, (lo, hi), bt.pvalue

def per_crater_boot(df, nb=NB):
    """Resample craters (each weighted equally); return mean-of-crater-mean poleward frac + CI."""
    g = df.groupby("host_uid")["poleward"].mean().values
    if len(g) == 0: return dict(mean=np.nan, ci=(np.nan, np.nan), n_craters=0)
    boots = np.array([RNG.choice(g, len(g), replace=True).mean() for _ in range(nb)])
    return dict(mean=float(g.mean()), ci=(float(np.percentile(boots,2.5)), float(np.percentile(boots,97.5))),
                n_craters=int(len(g)))

def summarize(df, label):
    n = len(df)
    if n == 0:
        return dict(label=label, n=0)
    pole = df["poleward"].values
    k = int(pole.sum())
    p, ci, pval = binom_ci(k, n)
    pcb = per_crater_boot(df)
    dev = df["equator_dev"].values
    ray = rayleigh(df["facing_aspect_deg"].values)
    vt = vtest(dev, 0.0)  # toward equator-facing
    # median facing aspect, north only meaningful vs Thomas; report per subset anyway
    med_face = circ_median(df["facing_aspect_deg"].values)
    med_ci = circ_median_boot(df["facing_aspect_deg"].values)
    med_dev = float(np.median(np.abs(dev)))
    return dict(label=label, n=int(n), n_craters=int(df["host_uid"].nunique()),
                k_poleward=k, frac_poleward=float(p),
                frac_poleward_ci95=[float(ci[0]), float(ci[1])], binom_p=float(pval),
                per_crater_frac_poleward=pcb["mean"], per_crater_ci95=list(pcb["ci"]),
                per_crater_n=pcb["n_craters"],
                rayleigh=ray, vtest_equatorfacing=vt,
                median_facing_aspect_deg=med_face, median_facing_aspect_ci95=list(med_ci),
                median_abs_equator_dev_deg=med_dev)

def main():
    df = pd.read_csv(f"{IN}/host_assignments.csv")
    df["hemi"] = np.where(df["lat"] >= 0, "N", "S")
    df["abslat"] = df["lat"].abs()
    df["poleward"] = poleward_flag(df["lat"].values, df["az_deg"].values)
    df["equator_dev"] = equator_dev(df["lat"].values, df["facing_aspect_deg"].values)
    # hot / warm pole by hollow longitude
    lon = df["lon"].values % 360
    d_hot = np.minimum(np.abs(wrap180(lon - 0)), np.abs(wrap180(lon - 180)))
    df["polepref"] = np.where(d_hot < 45, "hot", "warm")

    rec = {"note": "position-as-aspect proxy; 665 m DEM cannot resolve hollow-scale slopes",
           "thomas2014": {"n": 31, "median_aspect_deg": 187, "desc": "equator-facing in the north"},
           "n_assigned_hollows": int(len(df)),
           "rng_seed": 20260911, "n_bootstrap": NB}

    # ---- WALL hollows ----
    wall = df[(df["rho"] > 0.55) & (df["rho"] < 1.05) & (df["host_D_km"] >= 15)].copy()
    wall.to_csv(f"{G2_OUT}/aspect_wall_hollows.csv", index=False)
    rec["wall_definition"] = "0.55 < rho < 1.05 and host_D_km >= 15"
    rec["n_wall_hollows"] = int(len(wall))
    rec["wall_all"] = summarize(wall, "wall_all")
    rec["wall_by_hemisphere"] = {h: summarize(wall[wall.hemi == h], f"wall_{h}") for h in ["N", "S"]}

    # 15-deg |lat| bands, per hemisphere
    bands = {}
    edges = list(range(0, 90, 15))
    for h in ["N", "S"]:
        sub = wall[wall.hemi == h]
        bb = {}
        for e in edges:
            m = sub[(sub.abslat >= e) & (sub.abslat < e+15)]
            bb[f"{e}-{e+15}"] = summarize(m, f"wall_{h}_{e}-{e+15}")
        bands[h] = bb
    rec["wall_by_band"] = bands

    # hot vs warm pole (per hemisphere)
    hw = {}
    for h in ["N", "S"]:
        sub = wall[wall.hemi == h]
        hw[h] = {pp: summarize(sub[sub.polepref == pp], f"wall_{h}_{pp}") for pp in ["hot", "warm"]}
    rec["wall_hotwarm"] = hw

    # ---- FLOOR control (rho < 0.55) ----
    floor = df[(df["rho"] < 0.55) & (df["host_D_km"] >= 15)].copy()
    rec["floor_definition"] = "rho < 0.55 and host_D_km >= 15 (control: no expected preference)"
    rec["n_floor_hollows"] = int(len(floor))
    rec["floor_all"] = summarize(floor, "floor_all")
    rec["floor_by_hemisphere"] = {h: summarize(floor[floor.hemi == h], f"floor_{h}") for h in ["N", "S"]}

    # ---- Peak-ring / central-peak (rho < 0.3 in D > 100 km) ----
    peak = df[(df["rho"] < 0.3) & (df["host_D_km"] > 100)].copy()
    rec["peak_definition"] = "rho < 0.3 and host_D_km > 100 (central-peak/peak-ring)"
    rec["n_peak_hollows"] = int(len(peak))
    rec["peak_all"] = summarize(peak, "peak_all") if len(peak) else dict(n=0)

    # ---- Thomas 2014 comparison (north wall hollows) ----
    wn = wall[wall.hemi == "N"]
    med_n = circ_median(wn["facing_aspect_deg"].values)
    lo, hi = circ_median_boot(wn["facing_aspect_deg"].values)
    thomas_consistent = bool((lo <= 187 <= hi)) if not np.isnan(lo) else None
    rec["thomas_comparison"] = {
        "north_wall_median_facing_aspect_deg": med_n,
        "north_wall_median_ci95": [lo, hi],
        "thomas_median_deg": 187,
        "consistent_187_in_ci": thomas_consistent}

    with open(f"{G2_OUT}/aspect_record.json", "w") as f:
        json.dump(rec, f, indent=1)
    # console summary
    print("n_wall =", rec["n_wall_hollows"])
    for h in ["N", "S"]:
        w = rec["wall_by_hemisphere"][h]
        print(f"  {h}: frac_pole={w['frac_poleward']:.3f} CI{np.round(w['frac_poleward_ci95'],3)} "
              f"binom_p={w['binom_p']:.2e} percrater={w['per_crater_frac_poleward']:.3f}"
              f"{np.round(w['per_crater_ci95'],3)} med_face={w['median_facing_aspect_deg']:.0f} "
              f"Vp={w['vtest_equatorfacing']['p']:.2e}")
    print("floor N/S frac_pole:",
          rec["floor_by_hemisphere"]["N"].get("frac_poleward"),
          rec["floor_by_hemisphere"]["S"].get("frac_poleward"))
    print("thomas:", rec["thomas_comparison"])
    print("wrote aspect_record.json, aspect_wall_hollows.csv")

if __name__ == "__main__":
    main()
