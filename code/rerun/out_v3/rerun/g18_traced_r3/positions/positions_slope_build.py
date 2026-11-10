# g18 copy of g13_traced/positions/positions_slope_build.py (7 Nov 2026). Changed: folder g18_traced_r3; traces from traces/traces_clean_r3.json (round-2 snapshot, hosts of tracer groups A and B); G13_REV default r3 (slump counts as wall, as in r2); outputs named r3.
# g13 revision (7 Nov 2026): with G11_RULES=traced and G13_REV=r2 (default) outputs go to out_traced_r2; G13_REV=r1 keeps out_traced.
# g13 copy of g11_curvature/positions/positions_slope_build.py (6 Nov 2026). Changed: folder g13_traced/positions; G11_RULES default traced (traced outlines where available, curvature v3 elsewhere via g13_traced/traced_classes.py; env CLS=v2 for curvature v2 elsewhere).
# g11 copy of g9_draft_reruns/positions/positions_slope_build.py (6 Nov 2026). Changed: env G11_RULES (curv | fix; v3 catalog and grid fixed), output folder g11_curvature/positions/out_<RULES>/cache, the peak class (RULES=curv: g10 curv_classes.classify_curv), and new arrays peak_etab, hollows_etab, fields_etab (eta_b = (z - z_base)/(z_peak95 - z_base), z_base = 5th percentile of valid peak-cell elevations) and zbase_rel (z_base - z_floor per crater).
# g9 copy of analysis_workshop/melt_source_tests/positions/positions_slope_build.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID, host and catalog paths, completeness grid, output folder positions/out_<CAT>_<RULES>_grid<GRID>/cache, the peak-class line (RULES=fix: g7 terrain_classes.classify), and optional START STOP arguments (chunks in cache/_chunks, joined by merge_chunks.py) for the 180 s limit.
#!/usr/bin/env python3
"""Slope follow-up, step 1 (5 Nov 2026): peak-cell slope for the melt-source position test.
Same hosts, floor level, terrain classes (H = 250 m), superposed-crater masks and completeness weights as
positions_build.py / peak_wall_test/peak_wall.py. For every valid peak cell, peak hollow and peak hollow field (centroid of the
field's peak hollows, snapped to the nearest member hollow if the centroid cell is not peak) it records the weighted area a,
eta = h/(z_peak95 - z_floor), and DEM slope in degrees from central differences at a 1-cell baseline (+-1 cell, 1.33 km) and a
2-cell baseline (+-2 cells, 2.66 km); east-west spacing = 665.24 m x cos(latitude of the row). Per crater: number of valid peak
cells and z_peak95 - z_floor for the resolution screen. Writes cache/slope_H250.npz."""
import os, numpy as np, pandas as pd, rasterio
import os as _os, sys
import sys as _sys; _sys.dont_write_bytecode = True
G11 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g18_traced_r3")
CAT = "v3"; RULES = _os.environ.get("G11_RULES", "traced"); GRID = "v3"; assert RULES in ("curv", "fix", "traced"), RULES
ODIR = f"{G11}/positions/out_{RULES}" + ("_r3" if RULES == "traced" and _os.environ.get("G13_REV", "r3") == "r3" else ""); _os.makedirs(ODIR + "/cache/_chunks", exist_ok=True)
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))
sys.path.insert(0, _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g7_terrain_fix")); from terrain_classes import classify
_sys.path.insert(0, _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g10_threshold")); from curv_classes import classify_curv
_sys.path.insert(0, _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g18_traced_r3")); from traced_classes import classes_for
from rasterio.windows import from_bounds
HE = os.path.expanduser("~/mnt/hollows-emissivity"); OUT = ODIR + "/cache"
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM = 2439.4; KD = np.pi*RM/180; PX = 665.24
cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"] = cr.lon % 360
ha = pd.read_csv(HOSTS)
hc = pd.read_csv(CANON)[["hollow_id", "field_id", "is_crater_trace"]]
ha = ha.merge(hc, on="hollow_id", how="left", suffixes=("", "_c")); ha = ha[~ha.is_crater_trace.astype(bool)]
C = np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz" if GRID == "v2" else f"{HE}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz"); Cg = C["C"][list(C["Dref"]).index(1.0)]
hosts = cr[cr.unique_id.isin(set(ha.host_uid.dropna().astype(int))) & (cr.diam >= 30)].sort_values("unique_id").reset_index(drop=True)
hosts = hosts[hosts.diam <= 250].reset_index(drop=True)
HT = 250.0; RPK = 0.45
def gcd(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    return RM*np.arccos(np.clip(np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2-lo1)), -1, 1))
def slope(z, LAT, k):
    zp = np.pad(z, k, constant_values=np.nan); nr, nc = z.shape
    gx = (zp[k:k+nr, 2*k:2*k+nc]-zp[k:k+nr, 0:nc])/(2*k*PX*np.cos(np.radians(LAT)))
    gy = (zp[2*k:2*k+nr, k:k+nc]-zp[0:nr, k:k+nc])/(2*k*PX)
    return np.degrees(np.arctan(np.hypot(gx, gy)))
CR, PK, HO, FD = [], [], [], []; PKB, HOB, FDB, ZB = [], [], [], []
with rasterio.open(DEM) as ds:
  A0 = int(sys.argv[1]) if len(sys.argv) > 2 else 0; A1 = int(sys.argv[2]) if len(sys.argv) > 2 else len(hosts)
  for _, c in hosts.iloc[A0:A1].iterrows():
    R = c.diam/2; la0, lo0 = c.lat, c.lon360; cl = np.cos(np.radians(la0))
    dla = 1.25*R/KD; dlo = min(1.25*R/(KD*max(cl, 0.05)), 179)
    xs = lambda lo: np.radians(lo-180)*2439400.0; ys = lambda la: np.radians(la)*2439400.0
    lo_a, lo_b = lo0-dlo, lo0+dlo
    if lo_a < 0 or lo_b > 360 or la0+dla > 89 or la0-dla < -89: continue
    win = from_bounds(xs(lo_a), ys(la0-dla), xs(lo_b), ys(la0+dla), ds.transform)
    z = ds.read(1, window=win, boundless=True, fill_value=-32768).astype(float); z[z == -32768] = np.nan; z *= 0.5
    T = ds.window_transform(win); nr, nc = z.shape
    xc = T.c+T.a*(np.arange(nc)+.5); yc = T.f+T.e*(np.arange(nr)+.5)
    LON = 180+np.degrees(xc/2439400.0)[None, :]*np.ones((nr, 1)); LAT = np.degrees(yc/2439400.0)[:, None]*np.ones((1, nc))
    dx = (LON-lo0)*KD*cl; dy = (LAT-la0)*KD; rho = np.hypot(dx, dy)/R
    ring = c.cent_str in ("pr", "rpc", "pb")
    fz = z[(rho > 0.15) & (rho < 0.6)]; fz = fz[np.isfinite(fz)]
    if len(fz) < 20: continue
    floor = np.percentile(fz, 20); h = z-floor; H = HT
    if RULES == "fix":
        az = np.degrees(np.arctan2(dy, dx)) % 360
        KC = classify(z, rho, az, ring, HT, RPK); assert KC["F0"] == floor; peak = KC["peak"]
    elif RULES == "curv":
        az = np.degrees(np.arctan2(dy, dx)) % 360
        peak = classify_curv(z, rho, az, R, str(c.cent_str).strip())["peak"]
    elif RULES == "traced":   # g13
        az = np.degrees(np.arctan2(dy, dx)) % 360
        peak = classes_for(c.unique_id, z, rho, az, R, str(c.cent_str).strip(), LON, LAT, _os.environ.get("CLS", "v3"))["peak"]
    else:
        peak = (rho <= (0.6 if ring else RPK)) & (h > H)
    sub = cr[(cr.diam >= 5) & (cr.diam < 0.5*c.diam) & (np.abs(cr.lat-la0) < dla)]
    sub = sub[gcd(la0, lo0, sub.lat.values, sub.lon360.values) < 1.05*R]
    m = np.zeros_like(rho, bool)
    for s in sub.itertuples():
        m |= np.hypot((LON-s.lon360)*KD*np.cos(np.radians(s.lat)), (LAT-s.lat)*KD) < 1.1*s.diam/2
    fin = np.isfinite(z); valid = fin & ~m & (rho <= 1.05)
    ci = np.clip(((90-LAT)*10).astype(int), 0, 1799); cj = np.clip((LON*10).astype(int), 0, 3599)
    w = Cg[ci, cj]; w = np.where(np.isfinite(w), w, 0); a = (0.66524**2)*np.cos(np.radians(LAT))*w
    Pv = peak & valid; uid = int(c.unique_id)
    z95 = float(np.percentile(z[Pv], 95)) if Pv.any() else np.nan
    eta = (h/(z95-floor)) if Pv.any() else np.full(z.shape, np.nan)
    z5 = float(np.percentile(z[Pv], 5)) if Pv.any() else np.nan
    etab = ((z-z5)/(z95-z5)) if (Pv.any() and z95 - z5 > 0) else np.full(z.shape, np.nan)
    s1, s2 = slope(z, LAT, 1), slope(z, LAT, 2)
    CR.append((uid, c.diam, la0, lo0, str(c.cent_str).strip(), int(Pv.sum()), z95-floor if Pv.any() else np.nan, float(a[Pv].sum())))
    ZB.append(z5-floor if Pv.any() else np.nan); ip = np.where(Pv); PKB.append(etab[ip].astype(np.float32)); PK.append(np.column_stack([np.full(len(ip[0]), uid), a[ip], eta[ip], s1[ip], s2[ip], h[ip]]))
    hs = ha[ha.host_uid == c.unique_id]
    if len(hs) == 0: continue
    def cells(lat, lon):
        i = np.clip(np.round((T.f-ys(lat))/(-T.e)-.5).astype(int), 0, nr-1)
        j = np.clip(np.round((xs(lon % 360)-T.c)/T.a-.5).astype(int), 0, nc-1); return i, j
    hlat, hlon = hs.lat.values, hs.lon.values % 360; i0, j0 = cells(hlat, hlon)
    isp = peak[i0, j0] & ~m[i0, j0] & fin[i0, j0] & (rho[i0, j0] <= 1.05)
    if isp.any():
        HOB.append(etab[i0[isp], j0[isp]]); HO.append(np.column_stack([np.full(isp.sum(), uid), hs.hollow_id.values[isp], hs.field_id.values[isp], eta[i0[isp], j0[isp]], s1[i0[isp], j0[isp]], s2[i0[isp], j0[isp]]]))
        hdx = (hlon-lo0)*KD*cl; hdy = (hlat-la0)*KD
        for fid in np.unique(hs.field_id.values[isp]):
            ss = isp & (hs.field_id.values == fid); cx, cy = hdx[ss].mean(), hdy[ss].mean()
            fi, fj = cells(np.array([la0+cy/KD]), np.array([lo0+cx/(KD*cl)])); fb = 0
            if not (peak[fi[0], fj[0]] and not m[fi[0], fj[0]] and fin[fi[0], fj[0]]):
                k = np.argmin(np.hypot(hdx[ss]-cx, hdy[ss]-cy)); fi, fj = np.array([i0[ss][k]]), np.array([j0[ss][k]]); fb = 1
            FD.append([uid, fid, int(ss.sum()), fb, eta[fi[0], fj[0]], s1[fi[0], fj[0]], s2[fi[0], fj[0]]]); FDB.append(float(etab[fi[0], fj[0]]))
craters = np.array(CR, dtype=[("uid", "i8"), ("D", "f8"), ("lat", "f8"), ("lon", "f8"), ("cent_str", "U8"), ("n_peak_cells", "i8"), ("z95_rel", "f8"), ("A_peak", "f8")])
np.savez_compressed(f"{OUT}/slope_H250.npz" if (A0, A1) == (0, len(hosts)) else f"{OUT}/_chunks/slope_H250_{A0:03d}_{A1:03d}.npz", craters=craters, peak=np.vstack(PK).astype(np.float32) if PK else np.zeros((0, 6), np.float32), hollows=np.vstack(HO) if HO else np.zeros((0, 6)), fields=np.array(FD, float) if FD else np.zeros((0, 7)),
         peak_etab=np.concatenate(PKB) if PKB else np.zeros(0, np.float32), hollows_etab=np.concatenate(HOB) if HOB else np.zeros(0), fields_etab=np.array(FDB, float), zbase_rel=np.array(ZB, float))
print(len(craters), "craters;", sum(len(x) for x in PK), "peak cells;", sum(len(x) for x in HO), "peak hollows;", len(FD), "peak fields")
