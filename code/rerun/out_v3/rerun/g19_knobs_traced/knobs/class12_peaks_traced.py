# g19 copy of g13_traced/knobs/class12_peaks_traced.py (7 Nov 2026). Changed: traced peaks from g19 traced_classes (traces_clean_r3k.json: rounds 1 and 2, all groups); chunks to _chunks_r3k.
# g13 revision (7 Nov 2026): chunks for G13_REV=r2 (default) go to _chunks_r2.
# g13 copy of g11_curvature/terrain/class12_peaks_curv.py (6 Nov 2026). Changed: only two peak columns are computed, curvature v3 (A_peak_v3, n_peak_v3; g10 curv_classes_v3) and the g13 classes (A_peak_tr, n_peak_tr; Michael's traced outlines where traced, curvature v3 elsewhere, g13_traced/traced_classes.py); chunks go to g13_traced/knobs/_chunks/c12_START_STOP.csv.
# g11 copy of g10_threshold/class12_peaks.py (6 Nov 2026). Changed: added the curvature column (A_peak_curv, n_peak_curv) from curv_classes.classify_curv (g10_threshold, defaults nsec 24, hmin 100 m, sig_frac 0.03) next to the existing f250, d015, h015 and p015 columns, which are computed as before; chunks go to g11_curvature/terrain/_chunks/c12_START_STOP.csv.
#!/usr/bin/env python3
"""g10 (6 Nov 2026): peak-terrain hollow density of Kinczyk et al. (2020) Class 1-2 craters 30-250 km with the DEM terrain classes
of the terrain test (g7 terrain_classes.classify), so the knob comparison uses the same peak definition as Section 3.2 in place
of the 0.15 R radial proxy of t2_peak_hollows.py. Hosts and non-hosts alike; hollows count only if assigned to the crater
(host_assignments_v3); superposed Herrick craters >= 5 km masked at 1.1 radii; area weighted by the v3 completeness grid (1 km).
Thresholds: fixed 250 m; 0.15 x depth (>= 100 m); max(250 m, 0.15 x depth). Usage: class12_peaks.py START STOP"""
import sys, os, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{HE}/code/rerun/out_v3/rerun/g7_terrain_fix"); from terrain_classes import classify
sys.path.insert(0, f"{HE}/code/rerun/out_v3/rerun/g19_knobs_traced"); from traced_classes import classes_for, TR   # g13
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"; RM = 2439.4; KD = np.pi * RM / 180
cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"] = cr.lon % 360
kz = pd.read_csv(f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/kinczyk_degradation/kinczyk_matched.csv")
K = kz[(kz.D >= 30) & (kz.D <= 250) & kz.kclass.isin([1, 2])].merge(cr[["unique_id", "cent_str"]].rename(columns={"unique_id": "id"}), on="id", how="left").sort_values("id").reset_index(drop=True)
ha = pd.read_csv(f"{HE}/code/rerun/out_v3/host_assignments_v3.csv", low_memory=False)
G = np.load(f"{HE}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz"); Cg = G["C"][list(G["Dref"]).index(1.0)]
def gcd(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2); return RM * np.arccos(np.clip(np.sin(p1)*np.sin(p2) + np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2 - lo1)), -1, 1))
xs = lambda lo: np.radians(lo - 180) * 2439400.0; ys = lambda la: np.radians(la) * 2439400.0
a0, a1 = int(sys.argv[1]), int(sys.argv[2]); rows = []
with rasterio.open(DEM) as ds:
    for _, c in K.iloc[a0:a1].iterrows():
        R = c.D / 2; la0, lo0 = c.lat, c.lon % 360; cl = np.cos(np.radians(la0)); dla = 1.25 * R / KD; dlo = min(1.25 * R / (KD * max(cl, 0.05)), 179)
        if lo0 - dlo < 0 or lo0 + dlo > 360 or la0 + dla > 89 or la0 - dla < -89: continue
        win = from_bounds(xs(lo0 - dlo), ys(la0 - dla), xs(lo0 + dlo), ys(la0 + dla), ds.transform)
        z = ds.read(1, window=win, boundless=True, fill_value=-32768).astype(float); z[z == -32768] = np.nan; z *= 0.5
        T = ds.window_transform(win); nr, nc = z.shape; xc = T.c + T.a * (np.arange(nc) + .5); yc = T.f + T.e * (np.arange(nr) + .5)
        LON = 180 + np.degrees(xc / 2439400.0)[None, :] * np.ones((nr, 1)); LAT = np.degrees(yc / 2439400.0)[:, None] * np.ones((1, nc))
        dx = (LON - lo0) * KD * cl; dy = (LAT - la0) * KD; rho = np.hypot(dx, dy) / R; az = np.degrees(np.arctan2(dy, dx)) % 360
        fz = z[(rho > 0.15) & (rho < 0.6)]; fz = fz[np.isfinite(fz)]
        if len(fz) < 20: continue
        F0 = np.percentile(fz, 20); rim = [np.nanmax(np.where((rho >= 0.85) & (rho <= 1.2) & (az >= k * 10) & (az < k * 10 + 10), z, np.nan)) for k in range(36)]
        depth = float(np.nanmedian(rim) - F0); ring = str(c.cent_str).strip() in ("pr", "rpc", "pb")
        sub = cr[(cr.diam >= 5) & (cr.diam < 0.5 * c.D) & (np.abs(cr.lat - la0) < dla)]; sub = sub[gcd(la0, lo0, sub.lat.values, sub.lon360.values) < 1.05 * R]
        m = np.zeros_like(rho, bool)
        for s in sub.itertuples(): m |= np.hypot((LON - s.lon360) * KD * np.cos(np.radians(s.lat)), (LAT - s.lat) * KD) < 1.1 * s.diam / 2
        valid = np.isfinite(z) & ~m & (rho <= 1.05)
        w = Cg[np.clip(((90 - LAT) * 10).astype(int), 0, 1799), np.clip((LON * 10).astype(int), 0, 3599)]; w = np.where(np.isfinite(w), w, 0)
        pxa = (0.66524 ** 2) * np.cos(np.radians(LAT)); A = lambda msk: float(np.sum(pxa * w * (msk & valid)))
        hs = ha[ha.host_uid == c.id]
        hi = np.clip(np.round((T.f - ys(hs.lat.values)) / (-T.e) - .5).astype(int), 0, nr - 1); hj = np.clip(np.round((xs(hs.lon.values % 360) - T.c) / T.a - .5).astype(int), 0, nc - 1)
        row = dict(uid=int(c.id), D=c.D, kclass=int(c.kclass), cent_str=str(c.cent_str).strip(), depth_m=depth, n_assigned=len(hs))
        cent = str(c.cent_str).strip() if isinstance(c.cent_str, str) else ""
        pk3 = classes_for(-1, z, rho, az, R, cent, LON, LAT, "v3")["peak"]; pkt = classes_for(c.id, z, rho, az, R, cent, LON, LAT, "v3")["peak"]   # g13
        for tag, pk in (("v3", pk3), ("tr", pkt)):
            row[f"A_peak_{tag}"] = A(pk); row[f"n_peak_{tag}"] = int((pk[hi, hj] & ~m[hi, hj]).sum()) if len(hs) else 0
        row["traced"] = int(c.id) in TR
        row["A_p015"] = A(rho <= 0.15); row["n_p015"] = int(((rho[hi, hj] <= 0.15) & ~m[hi, hj]).sum()) if len(hs) else 0
        rows.append(row)
CH = f"{HERE}/_chunks_r3k"; os.makedirs(CH, exist_ok=True); pd.DataFrame(rows).to_csv(f"{CH}/c12_{a0:05d}_{a1:05d}.csv", index=False); print(len(K), "craters;", len(rows), "rows")   # g11: chunk folder
