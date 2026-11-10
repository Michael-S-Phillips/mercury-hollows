# g13 revision (7 Nov 2026): with G11_RULES=traced and G13_REV=r2 (default) outputs go to out_traced_r2; G13_REV=r1 keeps out_traced.
# g13 copy of g11_curvature/positions/positions_build.py (6 Nov 2026). Changed: folder g13_traced/positions; G11_RULES default traced (traced outlines where available, curvature v3 elsewhere via g13_traced/traced_classes.py; env CLS=v2 for curvature v2 elsewhere).
# g11 copy of g9_draft_reruns/positions/positions_build.py (6 Nov 2026). Changed: env G11_RULES (curv | fix; v3 catalog and v3 grid fixed), output folder g11_curvature/positions/out_<RULES>/cache, terrain classes (RULES=curv: g10 curv_classes.classify_curv(z, rho, az, R km, cent_str stripped), az as in g9), and new arrays peak_etab, hollows_etab, jitter_etab, fields_etab: eta_b = (z - z_base)/(z_peak95 - z_base), z_base = 5th percentile of the crater's valid peak-cell elevations (NaN if z_peak95 <= z_base); craters gain z5_rel = z_base - z_floor. Existing arrays keep their layout.
# g9 copy of g3_catalog/melt_source/positions/positions_build.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID, output folder positions/out_<CAT>_<RULES>_grid<GRID>/cache, completeness grid (v3 = g4 completeness_10ppd_v3.npz), and the terrain-class lines (RULES=fix: g7 terrain_classes.classify, az = degrees(arctan2(dy, dx)) % 360). d_nonfloor keeps its definition (nearest non-floor cell within 1.05 R); under the new classes the floor reaches 1.0 R, so its radial cut moves from 0.8 R to 1.0 R.
# g3_catalog copy of analysis_workshop/melt_source_tests/positions/positions_build.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Melt-source position test, step 1 (5 Nov 2026): per-cell and per-hollow positions within peak and floor terrain.
Hosts, floor level, terrain classes, superposed-crater masks and completeness weights are copied from
peak_wall_test/peak_wall.py (unchanged logic): 352 hollow-hosting Herrick (2018) craters 30-250 km across, USGS 665 m DEM.
For each host this script records
  peak cells (valid, unmasked): weighted area a = pixel area x completeness(1 km hollows), h = z - z_floor,
      eta = h / (z_peak95 - z_floor) with z_peak95 the 95th percentile of the crater's peak-cell elevations,
      q = completeness-weighted height quantile of the cell within the crater's peak terrain (0 = lowest, 1 = highest),
      d_base = distance (km) to the nearest unmasked floor-class cell (the peak base);
  floor cells: a, rho, d_wall / d_peak = distance (km) to the nearest unmasked wall- / peak-class cell,
      d_elev = distance to the nearest unmasked cell within 1.05 R that is > H above the floor (any non-floor high ground),
      d_nonfloor = distance to the nearest unmasked non-floor cell of any kind within 1.05 R (includes the 0.8 R radial cut);
  hollows assigned to the host: class and the same quantities at the hollow's DEM cell, plus K = 25 jittered copies
      (Gaussian 0.5 km per axis, seed fixed per crater) to propagate position uncertainty;
  hollow fields: for peak and floor separately, the centroid of the field's in-class hollows (if the centroid cell is
      not in the class, the member hollow nearest the centroid is used).
Distances are between cell centers in the local tangent plane used by peak_wall.py (dx scaled by cos(lat0)).
Usage: python3 positions_build.py START STOP [H]   (indices into the 352-host list; writes cache/H{H}_{START}_{STOP}.npz)"""
import sys, os, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
from scipy.spatial import cKDTree
import os as _os
import sys as _sys; _sys.dont_write_bytecode = True
G11 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced")
CAT = "v3"; RULES = _os.environ.get("G11_RULES", "traced"); GRID = "v3"; assert RULES in ("curv", "fix", "traced"), RULES
ODIR = f"{G11}/positions/out_{RULES}" + ("_r2" if RULES == "traced" and _os.environ.get("G13_REV", "r2") == "r2" else ""); _os.makedirs(ODIR + "/cache", exist_ok=True)
sys.path.insert(0, _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g7_terrain_fix")); from terrain_classes import classify
_sys.path.insert(0, _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g10_threshold")); from curv_classes import classify_curv
_sys.path.insert(0, _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced")); from traced_classes import classes_for
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))
HE = os.path.expanduser("~/mnt/hollows-emissivity")
OUT = ODIR + "/cache"
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM = 2439.4; KD = np.pi*RM/180
cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"] = cr.lon % 360
ha = pd.read_csv(HOSTS)
hc = pd.read_csv(CANON)[["hollow_id", "field_id", "is_crater_trace"]]
ha = ha.merge(hc, on="hollow_id", how="left", suffixes=("", "_c"))
ha = ha[~ha.is_crater_trace.astype(bool)]
C = np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz" if GRID == "v2" else f"{HE}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz"); Cg = C["C"][list(C["Dref"]).index(1.0)]
hosts = cr[cr.unique_id.isin(set(ha.host_uid.dropna().astype(int))) & (cr.diam >= 30)].sort_values("unique_id").reset_index(drop=True)
hosts = hosts[hosts.diam <= 250].reset_index(drop=True)
a0, a1 = int(sys.argv[1]), int(sys.argv[2]); HT = float(sys.argv[3]) if len(sys.argv) > 3 else 250.0; RPK = 0.45
K = 25; SIG = 0.5
def gcd(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    return RM*np.arccos(np.clip(np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2-lo1)), -1, 1))
def dist_to(src_mask, tgt_mask, X, Y):
    """distance (km) from every src cell to the nearest tgt cell (inf if no tgt cells)"""
    out = np.full(src_mask.shape, np.nan, np.float32)
    if not src_mask.any(): return out
    if not tgt_mask.any(): out[src_mask] = np.inf; return out
    t = cKDTree(np.column_stack([X[tgt_mask], Y[tgt_mask]]))
    d, _ = t.query(np.column_stack([X[src_mask], Y[src_mask]])); out[src_mask] = d; return out
CR, PK, FL, HO, FD = [], [], [], [], []; PKB, HOB, FDB = [], [], []
with rasterio.open(DEM) as ds:
  for _, c in hosts.iloc[a0:a1].iterrows():
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
        KC = classify(z, rho, az, ring, HT, RPK); assert KC["F0"] == floor
        peak, wall, flo = KC["peak"], KC["wall"], KC["flo"]
    elif RULES == "curv":
        az = np.degrees(np.arctan2(dy, dx)) % 360
        KC = classify_curv(z, rho, az, R, str(c.cent_str).strip())
        peak, wall, flo = KC["peak"], KC["wall"], KC["flo"]
    elif RULES == "traced":   # g13: traced outlines where available, curvature v3 elsewhere
        az = np.degrees(np.arctan2(dy, dx)) % 360
        KC = classes_for(c.unique_id, z, rho, az, R, str(c.cent_str).strip(), LON, LAT, _os.environ.get("CLS", "v3"))
        peak, wall, flo = KC["peak"], KC["wall"], KC["flo"]
    else:
        peak = (rho <= (0.6 if ring else RPK)) & (h > H)
        wall = (rho >= (0.7 if ring else 0.55)) & (rho <= 1.05) & (h > H)
        flo = (rho <= 0.8) & (h <= H)
    sub = cr[(cr.diam >= 5) & (cr.diam < 0.5*c.diam) & (np.abs(cr.lat-la0) < dla)]
    sub = sub[gcd(la0, lo0, sub.lat.values, sub.lon360.values) < 1.05*R]
    m = np.zeros_like(rho, bool)
    for s in sub.itertuples():
        m |= np.hypot((LON-s.lon360)*KD*np.cos(np.radians(s.lat)), (LAT-s.lat)*KD) < 1.1*s.diam/2
    fin = np.isfinite(z); valid = fin & ~m & (rho <= 1.05)
    ci = np.clip(((90-LAT)*10).astype(int), 0, 1799); cj = np.clip((LON*10).astype(int), 0, 3599)
    w = Cg[ci, cj]; w = np.where(np.isfinite(w), w, 0)
    pxa = (0.66524**2)*np.cos(np.radians(LAT)); a = pxa*w
    # targets (unmasked, finite)
    tP = peak & fin & ~m; tW = wall & fin & ~m; tF = flo & fin & ~m
    tE = (h > H) & fin & ~m & (rho <= 1.05); tN = (~flo) & fin & ~m & (rho <= 1.05)
    Pv = peak & valid; Fv = flo & valid
    d_base = dist_to(Pv, tF, dx, dy)
    d_wall = dist_to(Fv, tW, dx, dy); d_peak = dist_to(Fv, tP, dx, dy)
    d_elev = dist_to(Fv, tE, dx, dy); d_nonf = dist_to(Fv, tN, dx, dy)
    zp = z[Pv]; z95 = float(np.percentile(zp, 95)) if Pv.sum() >= 1 else np.nan
    z5 = float(np.percentile(zp, 5)) if Pv.sum() >= 1 else np.nan; etab = np.full(z.shape, np.nan, np.float32)
    if Pv.sum() >= 1 and z95 - z5 > 0: etab[Pv] = ((z[Pv]-z5)/(z95-z5)).astype(np.float32)
    eta = np.full(z.shape, np.nan, np.float32); qq = np.full(z.shape, np.nan, np.float32)
    if Pv.sum() >= 1:
        eta[Pv] = (h[Pv]/(z95-floor)).astype(np.float32)
        hv = h[Pv]; av = a[Pv]; o = np.argsort(hv, kind="stable"); cs = np.cumsum(av[o]); At = cs[-1]
        if At > 0:
            qv = np.empty(len(hv)); qv[o] = (cs-0.5*av[o])/At; qq[Pv] = qv
    uid = int(c.unique_id)
    CR.append(dict(uid=uid, D=c.diam, lat=la0, lon=lo0, cent_str=str(c.cent_str).strip(), floor_m=floor, z95_rel=z95-floor,
                   A_peak=float(a[Pv].sum()), A_floor=float(a[Fv].sum()), A_wall=float(a[wall & valid].sum()),
                   cl=cl, H=HT, z5_rel=z5-floor))
    ip = np.where(Pv); PKB.append(etab[ip]); PK.append(np.column_stack([np.full(len(ip[0]), uid), a[ip], h[ip], eta[ip], qq[ip], d_base[ip], rho[ip]]))
    ifl = np.where(Fv); FL.append(np.column_stack([np.full(len(ifl[0]), uid), a[ifl], h[ifl], rho[ifl], d_wall[ifl], d_peak[ifl], d_elev[ifl], d_nonf[ifl]]))
    # hollows
    hs = ha[ha.host_uid == c.unique_id]
    if len(hs) == 0: continue
    def cells(lat, lon):
        i = np.clip(np.round((T.f-ys(lat))/(-T.e)-.5).astype(int), 0, nr-1)
        j = np.clip(np.round((xs(lon % 360)-T.c)/T.a-.5).astype(int), 0, nc-1); return i, j
    def clsof(i, j):
        k = np.where(peak[i, j], 1, np.where(wall[i, j], 2, np.where(flo[i, j], 3, 0)))   # 1 peak 2 wall 3 floor 0 other
        k = np.where(m[i, j] | ~fin[i, j] | (rho[i, j] > 1.05), -1, k); return k
    def vals(i, j):
        return [clsof(i, j), eta[i, j], qq[i, j], d_base[i, j], d_wall[i, j], d_peak[i, j], d_elev[i, j], d_nonf[i, j], rho[i, j]]
    hlat, hlon = hs.lat.values, hs.lon.values % 360
    i0, j0 = cells(hlat, hlon); v0 = vals(i0, j0)
    hdx = (hlon-lo0)*KD*cl; hdy = (hlat-la0)*KD
    rng = np.random.default_rng(uid)
    ex = rng.normal(0, SIG, (K, len(hs))); ey = rng.normal(0, SIG, (K, len(hs)))
    jl = la0+(hdy[None, :]+ey)/KD; jo = lo0+(hdx[None, :]+ex)/(KD*cl)
    ij, jj = cells(jl.ravel(), jo.ravel()); vj = vals(ij, jj)
    rec = np.column_stack([np.full(len(hs), uid), hs.hollow_id.values, hs.field_id.values] + [np.asarray(x, float) for x in v0])
    jrec = np.stack([np.asarray(x, float).reshape(K, len(hs)) for x in vj], axis=0)   # (9, K, n)
    HO.append((rec, jrec)); HOB.append((etab[i0, j0].astype(float), etab[ij, jj].reshape(K, len(hs)).astype(np.float32)))
    # fields: centroid of in-class hollows
    for kc in (1, 3):
        sel = v0[0] == kc
        if not sel.any(): continue
        for fid in np.unique(hs.field_id.values[sel]):
            ss = sel & (hs.field_id.values == fid)
            cx, cy = hdx[ss].mean(), hdy[ss].mean()
            fi, fj = cells(np.array([la0+cy/KD]), np.array([lo0+cx/(KD*cl)])); fb = 0
            if clsof(fi, fj)[0] != kc:
                k = np.argmin(np.hypot(hdx[ss]-cx, hdy[ss]-cy)); fi, fj = np.array([i0[ss][k]]), np.array([j0[ss][k]]); fb = 1
            vv = vals(fi, fj)
            FD.append([uid, fid, kc, int(ss.sum()), fb] + [float(np.asarray(x)[0]) for x in vv[1:]]); FDB.append(float(etab[fi, fj][0]))
cr_df = pd.DataFrame(CR)
pk = np.vstack(PK).astype(np.float32) if PK else np.zeros((0, 7), np.float32)
fl = np.vstack(FL).astype(np.float32) if FL else np.zeros((0, 8), np.float32)
ho = np.vstack([r for r, _ in HO]) if HO else np.zeros((0, 12))
jh = np.concatenate([j for _, j in HO], axis=2).astype(np.float32) if HO else np.zeros((9, K, 0), np.float32)
fd = np.array(FD, float) if FD else np.zeros((0, 13))
os.makedirs(OUT, exist_ok=True)
np.savez_compressed(f"{OUT}/H{int(HT)}_{a0:03d}_{a1:03d}.npz", craters=cr_df.to_records(index=False), peak=pk, floor=fl, hollows=ho, jitter=jh, fields=fd,
         peak_etab=np.concatenate(PKB).astype(np.float32) if PKB else np.zeros(0, np.float32), hollows_etab=np.concatenate([a for a, _ in HOB]) if HOB else np.zeros(0),
         jitter_etab=np.concatenate([b for _, b in HOB], axis=1) if HOB else np.zeros((K, 0), np.float32), fields_etab=np.array(FDB, float))
print(len(hosts), "hosts;", len(cr_df), "processed;", pk.shape, fl.shape, ho.shape, fd.shape)
