#!/usr/bin/env python3
"""g13 (7 Nov 2026), new script: in the hosts where Michael drew slump (traces_clean_r2.json), completeness-weighted area and hollows on
slump, on the rest of the wall (outside the wall foot to 1.05 R) and on floor, with the terrain test's mask and weights. Descriptive:
11 hosts, no interval. Writes slump_descriptive.json."""
import os, sys, json, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.environ["G13_REV"] = "r2"
from traced_classes import TR, aeqd, _inside
HE = os.path.expanduser("~/mnt/hollows-emissivity"); RM = 2439.4; KD = np.pi * RM / 180
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"] = cr.lon % 360
ha = pd.read_csv(f"{HE}/code/rerun/out_v3/host_assignments_v3.csv", low_memory=False)
G = np.load(f"{HE}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz"); Cg = G["C"][list(G["Dref"]).index(1.0)]
gcd = lambda la1, lo1, la2, lo2: RM * np.arccos(np.clip(np.sin(np.radians(la1))*np.sin(np.radians(la2)) + np.cos(np.radians(la1))*np.cos(np.radians(la2))*np.cos(np.radians(lo2 - lo1)), -1, 1))
xs = lambda lo: np.radians(lo - 180) * 2439400.0; ys = lambda la: np.radians(la) * 2439400.0; rows = []
with rasterio.open(DEM) as ds:
    for uid, t in TR.items():
        if not t.get("slump"): continue
        c = cr[cr.unique_id == uid].iloc[0]; R = c.diam / 2; la0, lo0 = c.lat, c.lon360; cl = np.cos(np.radians(la0)); dla = 1.25 * R / KD; dlo = 1.25 * R / (KD * cl)
        win = from_bounds(xs(lo0 - dlo), ys(la0 - dla), xs(lo0 + dlo), ys(la0 + dla), ds.transform); z = ds.read(1, window=win, boundless=True, fill_value=-32768).astype(float); z[z == -32768] = np.nan
        T = ds.window_transform(win); nr, nc = z.shape; xc = T.c + T.a * (np.arange(nc) + .5); yc = T.f + T.e * (np.arange(nr) + .5)
        LON = 180 + np.degrees(xc / 2439400.0)[None, :] * np.ones((nr, 1)); LAT = np.degrees(yc / 2439400.0)[:, None] * np.ones((1, nc))
        rho = np.hypot((LON - lo0) * KD * cl, (LAT - la0) * KD) / R; X, Y = aeqd(LON, LAT, t["lat0"], t["lon0"])
        pk = _inside(t["peak"], X, Y, t); ft = _inside(t["foot"], X, Y, t); sl = _inside(t["slump"], X, Y, t) & ~pk & (rho <= 1.05)
        wall = (~ft) & ~sl & ~pk & (rho <= 1.05); flo = ft & ~sl & ~pk
        sub = cr[(cr.diam >= 5) & (cr.diam < 0.5 * c.diam) & (np.abs(cr.lat - la0) < dla)]; sub = sub[gcd(la0, lo0, sub.lat.values, sub.lon360.values) < 1.05 * R]
        m = np.zeros_like(rho, bool)
        for q in sub.itertuples(): m |= np.hypot((LON - q.lon360) * KD * np.cos(np.radians(q.lat)), (LAT - q.lat) * KD) < 1.1 * q.diam / 2
        valid = np.isfinite(z) & ~m; w = Cg[np.clip(((90 - LAT) * 10).astype(int), 0, 1799), np.clip((LON * 10).astype(int), 0, 3599)]; w = np.where(np.isfinite(w), w, 0)
        pxa = (0.66524 ** 2) * np.cos(np.radians(LAT)) * w; A = lambda k: float((pxa * (k & valid)).sum())
        hs = ha[ha.host_uid == uid]; hi = np.clip(np.round((T.f - ys(hs.lat.values)) / (-T.e) - .5).astype(int), 0, nr - 1); hj = np.clip(np.round((xs(hs.lon.values % 360) - T.c) / T.a - .5).astype(int), 0, nc - 1)
        n = lambda k: int((k[hi, hj] & valid[hi, hj]).sum()); sl_in = sl & ft
        rows.append(dict(uid=uid, name=t["name"], A_slump=A(sl), n_slump=n(sl), A_slump_inside_foot=A(sl_in), n_slump_inside_foot=n(sl_in), A_wall_other=A(wall), n_wall_other=n(wall), A_floor=A(flo), n_floor=n(flo), A_peak=A(pk), n_peak=n(pk)))
d = pd.DataFrame(rows); tot = d.drop(columns=["uid", "name"]).sum()
dens = {k: round(1000 * tot[f"n_{k}"] / tot[f"A_{k}"], 2) for k in ("slump", "wall_other", "floor", "peak")}
json.dump(dict(hosts=d.round(1).to_dict("records"), totals=tot.round(1).to_dict(), density_per_1000km2=dens), open(f"{HERE}/slump_descriptive.json", "w"), indent=1)
print(d[["name", "n_slump", "A_slump", "n_wall_other", "A_wall_other", "n_floor", "n_peak"]].round(0).to_string(index=False)); print("density per 1,000 km2", dens, "| hollows", {k: int(tot[f"n_{k}"]) for k in ("slump", "wall_other", "floor", "peak")})
