import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Hollow occurrence against regional slope and slope aspect inside hollow-hosting craters 30-250 km across.
Same host set, terrain classes, superposed-crater mask, and completeness weighting as ../peak_wall_test/peak_wall.py.
Slope and aspect from the USGS 665 m MESSENGER DEM, with the east-west pixel spacing scaled by cos(latitude):
  s1 = central differences on the native grid (1.33 km baseline);
  s3 = after a 3 x 3 mean filter (about 2.7 km baseline).
Aspect classes use the downslope (facing) azimuth: eq = within 45 deg of the equator direction, pole = within 45 deg of the pole
direction, ew = east- or west-facing.
Output (long format): one row per crater x baseline x terrain x aspect x slope bin, with completeness-weighted area (km2),
hollow count, and number of distinct hollow fields. Usage: python3 slope_occurrence.py START STOP OUTDIR"""
import sys, os, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import uniform_filter
HE = os.path.expanduser(_REPO + "")
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM = 2439.4; KD = np.pi * RM / 180; PX = 665.24
SB = [0, 5, 10, 15, 20, 25, 30, 90]
cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"] = cr.lon % 360
ha = pd.read_csv(f"{HE}/code/rerun/latitude/aspect/host_assignments.csv")
C = np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz"); Cg = C["C"][list(C["Dref"]).index(1.0)]
hosts = cr[cr.unique_id.isin(set(ha.host_uid.dropna().astype(int))) & (cr.diam >= 30) & (cr.diam <= 250)].sort_values("unique_id").reset_index(drop=True)
a0, a1, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
HT, RPK = 250.0, 0.45
def gcd(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    return RM * np.arccos(np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(np.radians(lo2 - lo1)), -1, 1))
def slope_aspect(z, lat_rows):
    gr, gc = np.gradient(z)                                      # per pixel; rows run north to south
    dzdn = -gr / PX; dzde = gc / (PX * np.cos(np.radians(lat_rows))[:, None])
    s = np.degrees(np.arctan(np.hypot(dzde, dzdn)))
    az = np.degrees(np.arctan2(-dzde, -dzdn)) % 360               # azimuth the slope faces (steepest descent)
    return s, az
rows = []
with rasterio.open(DEM) as ds:
    for _, c in hosts.iloc[a0:a1].iterrows():
        R = c.diam / 2; la0, lo0 = c.lat, c.lon360; cl = np.cos(np.radians(la0))
        dla = 1.25 * R / KD; dlo = min(1.25 * R / (KD * max(cl, 0.05)), 179)
        xs = lambda lo: np.radians(lo - 180) * 2439400.0; ys = lambda la: np.radians(la) * 2439400.0
        lo_a, lo_b = lo0 - dlo, lo0 + dlo
        if lo_a < 0 or lo_b > 360 or la0 + dla > 89 or la0 - dla < -89: continue
        win = from_bounds(xs(lo_a), ys(la0 - dla), xs(lo_b), ys(la0 + dla), ds.transform)
        z = ds.read(1, window=win, boundless=True, fill_value=-32768).astype(float); z[z == -32768] = np.nan; z *= 0.5
        T = ds.window_transform(win); nr, nc = z.shape
        xc = T.c + T.a * (np.arange(nc) + .5); yc = T.f + T.e * (np.arange(nr) + .5)
        latr = np.degrees(yc / 2439400.0)
        LON = 180 + np.degrees(xc / 2439400.0)[None, :] * np.ones((nr, 1)); LAT = latr[:, None] * np.ones((1, nc))
        dx = (LON - lo0) * KD * cl; dy = (LAT - la0) * KD; rho = np.hypot(dx, dy) / R
        ring = c.cent_str in ("pr", "rpc", "pb")
        fz = z[(rho > 0.15) & (rho < 0.6)]; fz = fz[np.isfinite(fz)]
        if len(fz) < 20: continue
        floor = np.percentile(fz, 20); h = z - floor
        peak = (rho <= (0.6 if ring else RPK)) & (h > HT)
        wall = (rho >= (0.7 if ring else 0.55)) & (rho <= 1.05) & (h > HT)
        flo = (rho <= 0.8) & (h <= HT)
        terr = np.where(peak, 1, np.where(wall, 2, np.where(flo, 3, 0)))      # 0 other, 1 peak, 2 wall, 3 floor
        sub = cr[(cr.diam >= 5) & (cr.diam < 0.5 * c.diam) & (np.abs(cr.lat - la0) < dla)]
        sub = sub[gcd(la0, lo0, sub.lat.values, sub.lon360.values) < 1.05 * R]
        m = np.zeros_like(rho, bool)
        for s_ in sub.itertuples():
            m |= np.hypot((LON - s_.lon360) * KD * np.cos(np.radians(s_.lat)), (LAT - s_.lat) * KD) < 1.1 * s_.diam / 2
        valid = np.isfinite(z) & ~m & (rho <= 1.05)
        ci = np.clip(((90 - LAT) * 10).astype(int), 0, 1799); cj = np.clip((LON * 10).astype(int), 0, 3599)
        w = Cg[ci, cj]; w = np.where(np.isfinite(w), w, 0)
        pxa = (PX / 1000) ** 2 * np.cos(np.radians(LAT)) * w
        hs = ha[ha.host_uid == c.unique_id]
        hi = np.clip(np.round((T.f - ys(hs.lat.values)) / (-T.e) - .5).astype(int), 0, nr - 1)
        hj = np.clip(np.round((xs(hs.lon.values % 360) - T.c) / T.a - .5).astype(int), 0, nc - 1)
        hvalid = valid[hi, hj]
        eq_az = 180.0 if la0 >= 0 else 0.0
        zf = np.where(np.isfinite(z), z, np.nanmean(z))
        zs = uniform_filter(zf, 3, mode="nearest"); zs[~np.isfinite(z)] = np.nan
        for bl, zz in (("s1", z), ("s3", zs)):
            s, az = slope_aspect(zz, latr)
            rel = (az - eq_az + 180) % 360 - 180
            asp = np.where(np.abs(rel) <= 45, 1, np.where(np.abs(rel) >= 135, 2, 3))   # 1 eq, 2 pole, 3 ew
            sb = np.clip(np.digitize(s, SB) - 1, 0, len(SB) - 2); sb[~np.isfinite(s)] = -1
            key = (terr * 4 + asp) * 10 + sb
            ok = valid & (sb >= 0)
            A = np.bincount(key[ok].ravel(), weights=pxa[ok].ravel(), minlength=160)
            hk = key[hi, hj]; hok = hvalid & (sb[hi, hj] >= 0)
            N = np.bincount(hk[hok], minlength=160)
            F = {}
            for k_, f_ in zip(hk[hok], hs.field_id.values[hok]): F.setdefault(int(k_), set()).add(f_)
            for k_ in np.where((A > 0) | (N > 0))[0]:
                rows.append(dict(uid=int(c.unique_id), D=c.diam, lat=la0, lon=lo0, baseline=bl, terrain=["other", "peak", "wall", "floor"][k_ // 40],
                                 aspect=["", "eq", "pole", "ew"][(k_ // 10) % 4], sbin=SB[k_ % 10], A_km2=float(A[k_]), n=int(N[k_]),
                                 nf=len(F.get(int(k_), ())), n_host_total=int(len(hs)), n_host_valid=int(hvalid.sum())))
pd.DataFrame(rows).to_csv(f"{OUT}/occ_{a0:04d}_{a1:04d}.csv", index=False)
print(len(hosts), "hosts;", len(set(r["uid"] for r in rows)), "processed in", a0, a1)
