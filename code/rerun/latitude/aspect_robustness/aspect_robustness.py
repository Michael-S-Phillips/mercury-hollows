#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
"""Aspect robustness: decompose the 15-30N frac_poleward reversal and re-check
the 30-75N trend on the corrected (v2) catalogue.

Wall/host definitions taken from ../aspect/aspect_record.json:
  wall  = 0.55 < rho < 1.05 and host_D_km >= 15
  host  = smallest Herrick crater whose 1.1 R contains the hollow
  poleward flag = sign(lat)*cos(radians(az_deg)) > 0  (hollow on poleward
                  half of crater = sits on the equator-facing wall)
  frac_poleward = fraction of wall hollows that are poleward (equator-facing).

Produces aspect_robustness_record.json + aspect_robustness.png.
"""
import json, numpy as np, pandas as pd
from scipy import stats

RNG = np.random.default_rng(20260911)
NB = 10000
R_M = 2439.7  # km
CAL_LAT, CAL_LON = 31.9, 161.9  # Caloris center

BASE = _REPO + ""
HA = f"{BASE}/code/rerun/latitude/aspect/host_assignments.csv"
CANON = f"{BASE}/hollow-networks/data_v2/hollows_canonical.csv"
HERR = [f"{BASE}/data/herrick-craters/2017je005516r_5-10_kilometer_diameter_craters.xlsx",
        f"{BASE}/data/herrick-craters/2017je005516r_craters_over_10_km_in_diameter.xlsx"]
DEMPATH = open(f"{BASE}/code/rerun/latitude/aspect_robustness/dempath.txt").read().strip()

# ---------- helpers ----------
def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0

def gc_dist_km(lat1, lon1, lat2, lon2):
    la1, la2 = np.radians(lat1), np.radians(lat2)
    dlo = np.radians(lon2 - lon1)
    c = np.arccos(np.clip(np.sin(la1)*np.sin(la2) +
                          np.cos(la1)*np.cos(la2)*np.cos(dlo), -1, 1))
    return R_M * c

def poleward_flag(lat, az_deg):
    return np.sign(lat) * np.cos(np.radians(az_deg)) > 0

def per_crater_boot(df, nb=NB):
    g = df.groupby("host_uid")["poleward"].mean().values
    if len(g) == 0:
        return dict(mean=np.nan, ci=[np.nan, np.nan], n_craters=0)
    boots = np.array([RNG.choice(g, len(g), replace=True).mean() for _ in range(nb)])
    return dict(mean=float(g.mean()),
                ci=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                n_craters=int(len(g)))

def summarize(df, label):
    n = len(df)
    if n == 0:
        return dict(label=label, n=0)
    pole = df["poleward"].values.astype(bool)
    k = int(pole.sum())
    p = k / n
    bt = stats.binomtest(k, n, 0.5)
    lo, hi = bt.proportion_ci(confidence_level=0.95)
    pcb = per_crater_boot(df)
    return dict(label=label, n=int(n), n_craters=int(df["host_uid"].nunique()),
                k_poleward=k, frac_poleward=float(p),
                frac_poleward_ci95=[float(lo), float(hi)], binom_p=float(bt.pvalue),
                per_crater_frac_poleward=pcb["mean"], per_crater_ci95=pcb["ci"],
                per_crater_n=pcb["n_craters"])

# ---------- load & derive (matches aspect_stats.py) ----------
df = pd.read_csv(HA)
df["hemi"] = np.where(df["lat"] >= 0, "N", "S")
df["abslat"] = df["lat"].abs()
df["poleward"] = poleward_flag(df["lat"].values, df["az_deg"].values)
lon = df["lon"].values % 360
d_hot = np.minimum(np.abs(wrap180(lon - 0)), np.abs(wrap180(lon - 180)))
df["polepref"] = np.where(d_hot < 45, "hot", "warm")
df["cal_dist_km"] = gc_dist_km(df["lat"].values, df["lon"].values % 360, CAL_LAT, CAL_LON)

WALL = (df["rho"] > 0.55) & (df["rho"] < 1.05) & (df["host_D_km"] >= 15)
wall = df[WALL].copy()

rec = {"note": "aspect robustness: decompose 15-30N reversal, recheck 30-75N trend",
       "wall_definition": "0.55 < rho < 1.05 and host_D_km >= 15",
       "poleward_definition": "sign(lat)*cos(az_deg)>0 = equator-facing wall",
       "caloris_center": [CAL_LAT, CAL_LON], "R_Mercury_km": R_M,
       "rng_seed": 20260911, "n_bootstrap": NB,
       "n_wall_all": int(len(wall))}

# ---- reproduce recorded N bands (validation) ----
repro = {}
for e in [0, 15, 30, 45, 60]:
    m = wall[(wall.hemi == "N") & (wall.abslat >= e) & (wall.abslat < e+15)]
    repro[f"{e}-{e+15}"] = dict(n=int(len(m)), frac_poleward=float(m.poleward.mean()),
                                n_craters=int(m.host_uid.nunique()))
rec["reproduce_recorded_N_bands"] = repro
RECORDED = {"0-15": 0.4794007490636704, "15-30": 0.3901601830663616,
            "30-45": 0.6151368760064412, "45-60": 0.706060606060606, "60-75": 0.75}
rec["reproduce_match"] = {b: bool(abs(repro[b]["frac_poleward"] - RECORDED[b]) < 1e-6) for b in RECORDED}

print("REPRODUCE:", {b: round(repro[b]["frac_poleward"], 4) for b in repro})
print("MATCH:", rec["reproduce_match"])

# =====================================================================
# (1) DECOMPOSE the 15-30N reversal
# =====================================================================
b1530 = wall[(wall.hemi == "N") & (wall.abslat >= 15) & (wall.abslat < 30)].copy()
dec = {"band": "15-30N", "n": int(len(b1530)), "overall": summarize(b1530, "15-30N_all")}

# --- Caloris proximity ---
near = b1530[b1530.cal_dist_km < 1200]
far = b1530[b1530.cal_dist_km >= 1200]
dec["by_caloris"] = {"lt1200km": summarize(near, "15-30N_caloris_lt1200"),
                     "ge1200km": summarize(far, "15-30N_caloris_ge1200")}

# --- host D bins ---
dec["by_hostD"] = {}
for lo_d, hi_d, lab in [(15, 30, "D15-30"), (30, 100, "D30-100"), (100, 1e9, "D>100")]:
    m = b1530[(b1530.host_D_km >= lo_d) & (b1530.host_D_km < hi_d)]
    dec["by_hostD"][lab] = summarize(m, f"15-30N_{lab}")

# --- peak-ring basin hosts (D>100 km): list them ---
pr = b1530[b1530.host_D_km > 100]
prlist = (pr.groupby("host_uid")
          .agg(host_lat=("host_lat", "first"), host_lon=("host_lon", "first"),
               host_D_km=("host_D_km", "first"), n_wall_hollows=("poleward", "size"),
               k_poleward=("poleward", "sum")).reset_index())
prlist["frac_poleward"] = prlist.k_poleward / prlist.n_wall_hollows
dec["peakring_basin_hosts"] = prlist.round(3).to_dict(orient="records")
dec["by_peakring"] = {"is_peakring_D>100": summarize(pr, "15-30N_peakring"),
                      "not_peakring_D<=100": summarize(b1530[b1530.host_D_km <= 100], "15-30N_notpeakring")}

# --- hot/warm longitude ---
dec["by_hotwarm"] = {pp: summarize(b1530[b1530.polepref == pp], f"15-30N_{pp}")
                     for pp in ["hot", "warm"]}

rec["decompose_15_30N"] = dec
print("\n15-30N overall frac_pole=%.3f n=%d ncr=%d" %
      (dec["overall"]["frac_poleward"], dec["overall"]["n"], dec["overall"]["n_craters"]))
print("  Caloris <1200:", round(dec["by_caloris"]["lt1200km"]["frac_poleward"], 3),
      "n=", dec["by_caloris"]["lt1200km"]["n"], "| >=1200:",
      round(dec["by_caloris"]["ge1200km"].get("frac_poleward", float("nan")), 3),
      "n=", dec["by_caloris"]["ge1200km"].get("n"))
for lab in dec["by_hostD"]:
    s = dec["by_hostD"][lab]
    print(f"  {lab}: frac={s.get('frac_poleward')} n={s.get('n')} ncr={s.get('n_craters')} percrater={s.get('per_crater_frac_poleward')}")
print("  peakring hosts (D>100):", len(prlist), "craters, n_hollows=", int(pr.shape[0]),
      "frac=", round(dec["by_peakring"]["is_peakring_D>100"].get("frac_poleward", float('nan')), 3))
print("  hot:", round(dec["by_hotwarm"]["hot"].get("frac_poleward", float('nan')), 3),
      "n=", dec["by_hotwarm"]["hot"].get("n"),
      "| warm:", round(dec["by_hotwarm"]["warm"].get("frac_poleward", float('nan')), 3),
      "n=", dec["by_hotwarm"]["warm"].get("n"))

# =====================================================================
# (2) HOST-ASSIGNMENT check + reassignment with hosts restricted to D 20-100 km
# =====================================================================
from scipy.spatial import cKDTree

def unit_vec(lat_deg, lon_deg):
    la = np.radians(lat_deg); lo = np.radians(lon_deg)
    return np.column_stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])

# load craters (same source as assign_hosts.py)
cr = pd.concat([pd.read_excel(f)[["unique_id", "lat", "lon", "diam"]] for f in HERR],
               ignore_index=True)
cr["lon"] = cr["lon"] % 360.0
cr = cr.reset_index(drop=True)

# --- 2a: for the 15-30N reversal hollows, fraction assigned to D<30 host that is
#         superposed on a LARGER crater whose 1.1R also contains the hollow ---
cv_all = unit_vec(cr["lat"].values, cr["lon"].values)
crtree = cKDTree(cv_all)
crad_all = cr["diam"].values / 2.0
def contains_larger(row):
    """does a Herrick crater with diam > host_D contain this hollow within 1.1R?"""
    hv = unit_vec(np.array([row["lat"]]), np.array([row["lon"]]))[0]
    # candidate craters within generous radius
    idx = crtree.query_ball_point(hv, 2.0 * np.sin(1.1 * (cr["diam"].max()/2.0)/R_M/2.0))
    if not idx: return False
    idx = np.asarray(idx)
    d = gc_dist_km(cr["lat"].values[idx], cr["lon"].values[idx], row["lat"], row["lon"])
    larger = (cr["diam"].values[idx] > row["host_D_km"]) & (d <= 1.1 * crad_all[idx])
    return bool(larger.any())

small_host = b1530[b1530.host_D_km < 30].copy()
if len(small_host):
    small_host["superposed_on_larger"] = small_host.apply(contains_larger, axis=1)
    n_small = int(len(small_host))
    n_super = int(small_host["superposed_on_larger"].sum())
else:
    n_small = n_super = 0
rec["hostcheck_15_30N"] = {
    "n_band": int(len(b1530)),
    "n_host_D_lt30": n_small,
    "frac_of_band_host_D_lt30": n_small / len(b1530),
    "n_host_D_lt30_superposed_on_larger": n_super,
    "frac_D_lt30_superposed": (n_super / n_small) if n_small else None,
    "frac_of_band_smallhost_superposed": n_super / len(b1530)}
print("\n2a small-host check: D<30 hosts =%d/%d (%.3f); of those superposed on larger =%d (%.3f)"
      % (n_small, len(b1530), n_small/len(b1530), n_super, (n_super/n_small) if n_small else float('nan')))

# --- 2b: reassign ALL non-trace hollows using hosts restricted to 20<=D<=100 km ---
h = pd.read_csv(CANON)
h = h[~h["is_crater_trace"]].reset_index(drop=True)
crf = cr[(cr["diam"] >= 20) & (cr["diam"] <= 100)].reset_index(drop=True)
hv = unit_vec(h["lat"].values, h["lon"].values)
cvf = unit_vec(crf["lat"].values, crf["lon"].values)
htree = cKDTree(hv)
crad_km = crf["diam"].values / 2.0
ang = 1.1 * crad_km / R_M
chord = 2.0 * np.sin(np.clip(ang, 0, np.pi)/2.0)
best_ci = np.full(len(h), -1, dtype=int)
order = np.argsort(crf["diam"].values)      # smallest containing wins
for ci in order:
    idx = htree.query_ball_point(cvf[ci], chord[ci])
    if not idx: continue
    idx = np.asarray(idx)
    take = idx[best_ci[idx] == -1]
    if take.size:
        best_ci[take] = ci
assigned = best_ci >= 0
hs = h[assigned].copy().reset_index(drop=True)
ci = best_ci[assigned]
clat = crf["lat"].values[ci]; clon = crf["lon"].values[ci]; cd = crf["diam"].values[ci]
cuid = crf["unique_id"].values[ci]
la1 = np.radians(clat); lo1 = np.radians(clon)
la2 = np.radians(hs["lat"].values); lo2 = np.radians(hs["lon"].values)
dlo = lo2 - lo1
central = np.arccos(np.clip(np.sin(la1)*np.sin(la2) +
                            np.cos(la1)*np.cos(la2)*np.cos(dlo), -1, 1))
dist_km = R_M * central
bear = np.degrees(np.arctan2(np.sin(dlo)*np.cos(la2),
                             np.cos(la1)*np.sin(la2) - np.sin(la1)*np.cos(la2)*np.cos(dlo)))
az = bear % 360.0
hs["host_uid"] = cuid; hs["host_D_km"] = cd; hs["dist_km"] = dist_km
hs["rho"] = dist_km / (cd/2.0); hs["az_deg"] = az
hs["hemi"] = np.where(hs["lat"] >= 0, "N", "S")
hs["abslat"] = hs["lat"].abs()
hs["poleward"] = poleward_flag(hs["lat"].values, hs["az_deg"].values)
wall2 = hs[(hs["rho"] > 0.55) & (hs["rho"] < 1.05) & (hs["host_D_km"] >= 15)].copy()
print("2b reassigned (D 20-100 km hosts): wall hollows =", len(wall2),
      "assigned =", int(assigned.sum()), "/", len(h))

reassign = {"host_diam_restriction_km": [20, 100],
            "n_assigned": int(assigned.sum()), "n_wall": int(len(wall2)),
            "N_bands": {}, "S_bands": {}}
for e in [0, 15, 30, 45, 60, 75]:
    for hemi, key in [("N", "N_bands"), ("S", "S_bands")]:
        m = wall2[(wall2.hemi == hemi) & (wall2.abslat >= e) & (wall2.abslat < e+15)]
        reassign[key][f"{e}-{e+15}"] = summarize(m, f"reassign_{hemi}_{e}-{e+15}")
rec["reassign_D20_100"] = reassign
print("  N bands (D20-100 hosts):")
for e in [0, 15, 30, 45, 60]:
    s = reassign["N_bands"][f"{e}-{e+15}"]
    if s.get("n", 0):
        print(f"    {e}-{e+15}: frac={s['frac_poleward']:.3f} n={s['n']} ncr={s['n_craters']} "
              f"percrater={s['per_crater_frac_poleward']:.3f} CI{np.round(s['per_crater_ci95'],3)}")

# =====================================================================
# (3) ALTERNATIVE ASPECT PROXY: local DEM slope azimuth vs position facing
# =====================================================================
zz = np.load(DEMPATH, allow_pickle=True)
DEM = zz["DEM"].astype(np.float64)
ppd = int(zz["ppd"])
nrow, ncol = DEM.shape
# gradients in index space
grow, gcol = np.gradient(DEM)  # d/drow (south+), d/dcol (east+)
deg_per_px = 1.0 / ppd
def dem_facing_az(lat, lon):
    lon = lon % 360.0
    r = np.clip(((90.0 - lat) * ppd).astype(int), 0, nrow-1)
    c = np.clip(((lon) * ppd).astype(int), 0, ncol-1)
    # dz/dnorth: north = -row direction
    dz_dlat = -grow[r, c] / deg_per_px            # per degree lat
    dz_dnorth = dz_dlat / (R_M * np.radians(1.0)) # per km
    dz_dlon = gcol[r, c] / deg_per_px             # per degree lon (east+)
    cosl = np.cos(np.radians(lat))
    dz_deast = dz_dlon / (R_M * cosl * np.radians(1.0))
    # slope faces downhill: facing vector = -gradient (east, north)
    az = np.degrees(np.arctan2(-dz_deast, -dz_dnorth)) % 360.0
    return az, np.hypot(dz_deast, dz_dnorth)

waz, wslope = dem_facing_az(wall["lat"].values, wall["lon"].values)
wall["dem_facing_az"] = waz
wall["dem_slope"] = wslope
diff = np.abs(wrap180(wall["dem_facing_az"].values - wall["facing_aspect_deg"].values))
wall["dem_pos_diff_deg"] = diff
agree90 = diff < 90.0
agree45 = diff < 45.0
proxy = {"n": int(len(wall)),
         "agreement_rate_lt90deg": float(np.mean(agree90)),
         "agreement_rate_lt45deg": float(np.mean(agree45)),
         "median_abs_diff_deg": float(np.median(diff)),
         "note": ("DEM is 8 ppd (~5.3 km/px at equator); much coarser than hollow-scale "
                  "walls, so this bounds proxy reliability, not validates it at hollow scale.")}
# by band (N)
proxy["by_N_band"] = {}
for e in [0, 15, 30, 45, 60]:
    m = wall[(wall.hemi == "N") & (wall.abslat >= e) & (wall.abslat < e+15)]
    if len(m):
        proxy["by_N_band"][f"{e}-{e+15}"] = dict(
            n=int(len(m)), agreement_lt90=float((m.dem_pos_diff_deg < 90).mean()),
            median_diff=float(m.dem_pos_diff_deg.median()))
rec["dem_proxy_check"] = proxy
print("\n3 DEM slope-azimuth vs position facing: agree<90deg=%.3f  agree<45deg=%.3f  median|diff|=%.1f (n=%d)"
      % (proxy["agreement_rate_lt90deg"], proxy["agreement_rate_lt45deg"],
         proxy["median_abs_diff_deg"], proxy["n"]))

# =====================================================================
# (4) 30-75N trend with hosts D 20-100 km only (already in reassign N_bands)
# =====================================================================
trend = {"hosts": "D 20-100 km, smallest containing", "bands": {}}
for e in [30, 45, 60]:
    trend["bands"][f"{e}-{e+15}"] = reassign["N_bands"][f"{e}-{e+15}"]
# monotonic check per-hollow and per-crater
fh = [trend["bands"][b].get("frac_poleward") for b in ["30-45", "45-60", "60-75"] if trend["bands"][b].get("n")]
pc = [trend["bands"][b].get("per_crater_frac_poleward") for b in ["30-45", "45-60", "60-75"] if trend["bands"][b].get("n")]
trend["per_hollow_fracs_30_45_60_75"] = fh
trend["per_crater_fracs_30_45_60_75"] = pc
trend["per_hollow_increasing"] = bool(all(x < y for x, y in zip(fh, fh[1:]))) if len(fh) > 1 else None
trend["per_crater_increasing"] = bool(all(x < y for x, y in zip(pc, pc[1:]))) if len(pc) > 1 else None
rec["trend_30_75N_D20_100"] = trend
print("\n4 30-75N trend (D20-100 hosts) per-hollow:", [round(x,3) for x in fh],
      "per-crater:", [round(x,3) for x in pc])

with open(f"{BASE}/code/rerun/latitude/aspect_robustness/aspect_robustness_record.json", "w") as f:
    json.dump(rec, f, indent=1)
wall.to_csv(f"{BASE}/code/rerun/latitude/aspect_robustness/wall_with_dem_aspect.csv", index=False)
print("\nwrote aspect_robustness_record.json")
