#!/usr/bin/env python3
"""g12 (6 Nov 2026): assets for the host-crater tracing page. For the top hosts by terrain-test hollow count
(rows_curv_n24_h100_s0.03.csv, n_peak + n_wall + n_floor + n_other), writes per crater:
  img/<uid>.jpg  MDIS BDR 166 m mosaic resampled (bilinear) to a local azimuthal-equidistant grid, half-width 1.3 R, x east, y north,
                 pixel max(166.31 m, 2.6 R / 1600); stretched to the 1-99.5 percentiles of nonzero values.
  dem/<uid>.png  USGS 665 m DEM on the same projection at 665.24 m pixels; value = round(z - zmin) in R*256+G (m), B = 255 where valid.
  meta/<uid>.json  crater record, grid geometry, curvature-v3 peak/ring mask polygons (DEM pixel edges, lon/lat, unsimplified),
                 v3 wall-foot line (radius per azimuth, 2 deg steps), hollows assigned to the host (v3), superposed craters masked.
Curvature v3 is run exactly as compare_curv_v2v3.py load() does (1.25 R window on the native DEM grid).
Usage: build_tracer_assets.py [N_TOP=64] [time budget s=150]; skips craters whose three files exist."""
import sys, os, json, time, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window, from_bounds
from rasterio.features import shapes
from scipy.ndimage import map_coordinates
from PIL import Image
T0 = time.time(); NTOP = int(sys.argv[1]) if len(sys.argv) > 1 else 64; BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 150
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__)); RR = os.path.dirname(HERE)
sys.path.insert(0, f"{RR}/g10_threshold"); import curv_classes_v3 as V3
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"; BDR = f"{HE}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif"
RM = 2439.4; RMm = 2439400.0; KD = np.pi * RM / 180; OUT = f"{HERE}/assets"
cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"] = cr.lon % 360
ha = pd.read_csv(f"{HE}/code/rerun/out_v3/host_assignments_v3.csv", low_memory=False)
rows = pd.read_csv(f"{RR}/g11_curvature/terrain/out_v3/rows_curv_n24_h100_s0.03.csv")
rows["n_tot"] = rows.n_peak + rows.n_wall + rows.n_floor + rows.n_other
rows = rows.sort_values("n_tot", ascending=False).reset_index(drop=True); NALL = int(rows.n_tot.sum())
top = rows.head(NTOP)
def gcd(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2); return RM * np.arccos(np.clip(np.sin(p1)*np.sin(p2) + np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2 - lo1)), -1, 1))
def aeqd_inv(x, y, la0, lo0):
    """local azimuthal equidistant (km, x east, y north) -> lat, lon (lon continuous around lo0)"""
    rho = np.hypot(x, y); c = rho / RM; p0 = np.radians(la0)
    with np.errstate(invalid="ignore", divide="ignore"):
        lat = np.arcsin(np.cos(c) * np.sin(p0) + np.where(rho > 0, y * np.sin(c) * np.cos(p0) / rho, 0))
        lon = lo0 + np.degrees(np.arctan2(x * np.sin(c), rho * np.cos(p0) * np.cos(c) - y * np.sin(p0) * np.sin(c)))
    return np.degrees(lat), lon
def read_wrapped(ds, r0, r1, c0, c1, fill):
    """rows r0..r1-1, columns c0..c1-1 (may run past either edge; wraps in longitude)"""
    W = ds.width; r0c, r1c = max(r0, 0), min(r1, ds.height); parts = []; c = c0
    while c < c1:
        cm = c % W; n = min(c1 - c, W - cm)
        parts.append(ds.read(1, window=Window(cm, r0c, n, r1c - r0c))); c += n
    a = np.concatenate(parts, axis=1)
    if r0c > r0 or r1c < r1: a = np.pad(a, ((r0c - r0, r1 - r1c), (0, 0)), constant_values=fill)
    return a
def sample(ds, X, Y, fill, order=1):
    """X, Y in the dataset's projected metres (X continuous); bilinear sample"""
    T = ds.transform; col = (X - T.c) / T.a - 0.5; row = (Y - T.f) / T.e - 0.5
    c0, c1 = int(np.floor(col.min())) - 2, int(np.ceil(col.max())) + 3; r0, r1 = int(np.floor(row.min())) - 2, int(np.ceil(row.max())) + 3
    a = read_wrapped(ds, r0, r1, c0, c1, fill)
    return a, map_coordinates(a.astype(np.float32), [row - r0, col - c0], order=order, mode="nearest"), (row - r0, col - c0)
def one(c, ds, bs):
    uid = int(c.unique_id); R = c.diam / 2; la0, lo0 = float(c.lat), float(c.lon360); cl = np.cos(np.radians(la0)); cent = str(c.cent_str).strip() if isinstance(c.cent_str, str) else ""
    # ---- curvature v3 on the native DEM grid, as compare_curv_v2v3.load()
    xs = lambda lo: np.radians(lo - 180) * RMm; ys = lambda la: np.radians(la) * RMm
    dla = 1.25 * R / KD; dlo = 1.25 * R / (KD * cl)
    win = from_bounds(xs(lo0 - dlo), ys(la0 - dla), xs(lo0 + dlo), ys(la0 + dla), ds.transform)
    z = ds.read(1, window=win, boundless=True, fill_value=-32768).astype(float); z[z == -32768] = np.nan; z *= 0.5
    T = ds.window_transform(win); nr, nc = z.shape; xc = T.c + T.a * (np.arange(nc) + .5); yc = T.f + T.e * (np.arange(nr) + .5)
    LON = 180 + np.degrees(xc / RMm)[None, :] * np.ones((nr, 1)); LAT = np.degrees(yc / RMm)[:, None] * np.ones((1, nc))
    dx = (LON - lo0) * KD * cl; dy = (LAT - la0) * KD; rho = np.hypot(dx, dy) / R; az = np.degrees(np.arctan2(dy, dx)) % 360
    K = V3.classify_curv(z, rho, az, R, cent)
    polys = []
    for geom, val in shapes(K["peak"].astype(np.uint8), mask=K["peak"], transform=T):
        rings = [[[round(180 + float(np.degrees(x / RMm)), 5), round(float(np.degrees(y / RMm)), 5)] for x, y in ring] for ring in geom["coordinates"]]
        polys.append(rings)
    rwall = V3._med3(np.array([r["wall"] for r in K["info"]])); azg = np.arange(0, 360, 2.0) + 1.0; rw = V3._az(rwall, azg, 24) * R
    foot = [[round(lo0 + float(r * np.cos(np.radians(a)) / (KD * cl)), 5), round(la0 + float(r * np.sin(np.radians(a)) / KD), 5)] for a, r in zip(azg, rw)]
    # ---- display grids (local azimuthal equidistant)
    H = 1.3 * R; NB = int(min(1600, np.ceil(2 * H / 0.16631))); ND = int(np.ceil(2 * H / 0.66524))
    out = {}
    for name, N, src in (("img", NB, bs), ("dem", ND, ds)):
        px = 2 * H / N; g = -H + (np.arange(N) + .5) * px; X, Y = np.meshgrid(g, -g)   # row 0 at the north edge
        lat, lon = aeqd_inv(X, Y, la0, lo0)
        if name == "img":
            lp = ((lo0 + 180) % 360) - 180; lonc = lp + (lon - lo0)          # continuous pm180
            a, v, _ = sample(bs, np.radians(lonc) * RMm, np.radians(lat) * RMm, 0)
            nz = v[v > 0.5]; lo_, hi_ = np.percentile(nz, [1, 99.5]) if nz.size else (0, 255)
            im = np.clip((v - lo_) / max(hi_ - lo_, 1) * 255, 0, 255).astype(np.uint8)
            Image.fromarray(im, "L").save(f"{OUT}/img/{uid}.jpg", quality=86, optimize=True, progressive=True)
            out["img"] = dict(n=N, px_km=px, stretch=[float(lo_), float(hi_)])
        else:
            a, v, (rr, cc) = sample(ds, np.radians(lon - 180) * RMm, np.radians(lat) * RMm, -32768)
            ok = map_coordinates((a != -32768).astype(np.float32), [rr, cc], order=0, mode="nearest") > 0.5
            af = a.astype(np.float32); af[a == -32768] = np.nan; fillv = np.nanmedian(af); af[np.isnan(af)] = fillv
            zz = map_coordinates(af, [rr, cc], order=1, mode="nearest") * 0.5
            zmin = float(np.floor(zz[ok].min())) if ok.any() else 0.0; q = np.clip(np.round(zz - zmin), 0, 65535).astype(np.uint32)
            rgb = np.stack([(q >> 8).astype(np.uint8), (q & 255).astype(np.uint8), np.where(ok, 255, 0).astype(np.uint8)], -1)
            Image.fromarray(rgb, "RGB").save(f"{OUT}/dem/{uid}.png", optimize=True)
            out["dem"] = dict(n=N, px_km=px, zmin=zmin, zmax=float(zz[ok].max()) if ok.any() else 0.0)
    hs = ha[ha.host_uid == uid]; hol = [[round(float(lo % 360), 5), round(float(la), 5)] for lo, la in zip(hs.lon.values, hs.lat.values)]
    sub = cr[(cr.diam >= 5) & (cr.diam < 0.5 * c.diam) & (np.abs(cr.lat - la0) < 1.25 * R / KD)]
    sub = sub[gcd(la0, lo0, sub.lat.values, sub.lon360.values) < 1.05 * R]
    rr_ = rows[rows.uid == uid].iloc[0]
    meta = dict(uid=uid, name=(c["name"] if isinstance(c["name"], str) else ""), D=round(float(c.diam), 2), lat=round(la0, 5), lon=round(lo0, 5), cent=cent,
                filled=str(rr_.filled), floor_m=float(rr_.floor_m), depth_m=float(rr_.depth_m),
                n=dict(peak=int(rr_.n_peak), wall=int(rr_.n_wall), floor=int(rr_.n_floor), other=int(rr_.n_other), masked=int(rr_.n_masked), tot=int(rr_.n_tot)),
                H_km=H, grids=out, v3=dict(peak=polys, foot=foot, n_ring_sectors=int(K["n_ring_sectors"])),
                hollows=hol, masked=[[round(float(s.lon360), 5), round(float(s.lat), 5), round(float(s.diam), 2)] for s in sub.itertuples()])
    json.dump(meta, open(f"{OUT}/meta/{uid}.json", "w"), separators=(",", ":"))
    return uid
done = []
with rasterio.open(DEM) as ds, rasterio.open(BDR) as bs:
    for k, r in top.iterrows():
        uid = int(r.uid)
        if all(os.path.exists(f"{OUT}/{d}/{uid}.{e}") for d, e in (("img", "jpg"), ("dem", "png"), ("meta", "json"))): continue
        if time.time() - T0 > BUDGET: break
        c = cr[cr.unique_id == uid].iloc[0]; one(c, ds, bs); done.append(uid); print(k + 1, uid, r["name"], f"{c.diam:.0f} km", f"{time.time() - T0:.0f}s", flush=True)
left = [int(u) for u in top.uid if not os.path.exists(f"{OUT}/meta/{int(u)}.json")]
json.dump(dict(order=[int(u) for u in top.uid], n_all=NALL, n_hosts=len(rows), cum=[int(x) for x in top.n_tot.cumsum()]), open(f"{OUT}/order.json", "w"))
print("this call", len(done), "remaining", len(left), "n_all", NALL)
