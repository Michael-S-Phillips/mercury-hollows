#!/usr/bin/env python3
"""g23 step 3 (9 Nov 2026): shadow depth against hollow size at fixed pixel scale, to see whether depth rises to its full value
somewhere above 100 m. Detections (out_v2/detections_qc.csv; the v3 hollow ids) of non-trace hollows in frames with incidence 58-81 deg
and emission <= 50 deg, in two pixel-scale tiers, <= 7 m and 7-15 m; in each tier the finest detection of each hollow. Same measurement
as small_shadow.py (pilot shadow_measure v2 on the detection box), plus:
  - full_frac: the longest dark run inside the box on any inner line, as a fraction of the box length along the sun. A box that is
    >= 85% dark gives a lower limit, h >= full_frac x along_m / tan(i).
  - no shadow: an upper limit, h < 1.5 px x px x f / tan(i).
Resumable: `python3 depth_vs_size.py run` until it reports 0 left (EDRs fetched in chunks and deleted). Zoom tiles for the <= 7 m tier
in tiles_dvs/. Output depth_vs_size.csv."""
import sys, os, time, numpy as np, pandas as pd
from PIL import Image, ImageDraw
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import small_shadow as G
C, V, M = G.C, G.V, G.M
RR = G.RR; RUN = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun"); TD = f"{HERE}/tiles_dvs"; os.makedirs(TD, exist_ok=True)
OUT = f"{HERE}/depth_vs_size.csv"
def select():
    H = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool)]
    d = pd.read_csv(f"{RUN}/out_v2/detections_qc.csv"); d = d[d.hollow_id.isin(H.hollow_id)]
    d = d.merge(H[["hollow_id", "D_eq_km"]].rename(columns={"D_eq_km": "D_can_km"}), on="hollow_id")
    d = d[d.inc.between(58, 81) & (d.emi <= 50) & (d.px_m <= 15)].copy(); d["tier"] = np.where(d.px_m <= 7, "le7", "7to15")
    d = d.sort_values("px_m").drop_duplicates(["hollow_id", "tier"])
    b = C.load_bickel(with_corners=False)[["thumb", "num", "datetime", "x0", "x1", "y0", "y1", "lat", "lon180"]].drop_duplicates("thumb")
    s = d.drop(columns=["lat"]).merge(b.rename(columns={"lat": "lat_cm"}), on="thumb", how="left"); s["num"] = s.num.astype(str).str.zfill(10)
    return s.reset_index(drop=True)
def full_frac(o, dbg):
    if dbg is None: return np.nan
    sh, inbox, Wv = dbg["sh"], dbg["inbox"], dbg["W"]; wb = Wv[inbox.any(1)]
    if len(wb) == 0: return np.nan
    w0, w1 = wb.min() + 0.1 * (wb.max() - wb.min()), wb.max() - 0.1 * (wb.max() - wb.min()); best = 0.0
    for j in np.where((Wv >= w0) & (Wv <= w1))[0]:
        m = (sh[j] & inbox[j]).astype(np.int8)
        if m.any():
            dd = np.diff(np.r_[0, m, 0]); st, en = np.where(dd == 1)[0], np.where(dd == -1)[0]; best = max(best, (en - st).max() * M.STEP)
    return best / o["along_px"] if o["along_px"] > 0 else np.nan
def run():
    s = select(); done = pd.read_csv(OUT) if os.path.exists(OUT) else pd.DataFrame(columns=["thumb"])
    todo = s[~s.thumb.isin(done.thumb)]; nums = list(todo.num.unique()); t0 = time.time(); rows = []
    dts = todo.drop_duplicates("num").set_index("num").datetime
    while nums and time.time() - t0 < 90:
        chunk, nums = nums[:10], nums[10:]; V.fetch_many([(n, dts[n]) for n in chunk])
        for n in chunk:
            p = C.edr_path(n)
            if not os.path.exists(p):
                for _, r in todo[todo.num == n].iterrows(): rows.append(dict(thumb=r.thumb, hollow_id=r.hollow_id, flags="no_edr"))
                continue
            try: img, info = G.read_any(p)
            except Exception as e:
                print("read fail", n, e); os.remove(p); continue
            Hh, Ww = img.shape; fs = M.foreshortening(info)
            for _, r in todo[todo.num == n].iterrows():
                box = (r.x0 * Ww, r.x1 * Ww, r.y0 * Hh, r.y1 * Hh)
                o = M.measure(img, info, box, return_debug=True); dbg = o.pop("_dbg", None)
                inc = float(C.local_incidence(r.lat_cm, r.lon180, info["subsol_lat"], info["subsol_lon"])); ti = np.tan(np.radians(inc))
                ff = full_frac(o, dbg); along_m = o["along_px"] * info["px"] * fs
                Lm = o["L_px"] * info["px"] * fs if np.isfinite(o["L_px"]) else np.nan
                rows.append(dict(thumb=r.thumb, hollow_id=int(r.hollow_id), tier=r.tier, D_can_m=round(r.D_can_km * 1000, 1), D_det_m=round(r.D_eq_km * 1000, 1),
                                 num=n, px=info["px"], inc=round(inc, 2), emi=info["emi"], f=round(fs, 3), L_px=o["L_px"], along_px=round(o["along_px"], 1),
                                 L_m=round(Lm, 2) if np.isfinite(Lm) else np.nan, depth_m=round(Lm / ti, 2) if np.isfinite(Lm) else np.nan,
                                 along_m=round(along_m, 1), full_frac=round(ff, 3) if np.isfinite(ff) else np.nan,
                                 h_upper_m=round(1.5 * info["px"] * fs / ti, 2), h_fill_m=round(along_m / ti, 2),
                                 contrast=round((o["B"] - o["S"]) / max(o["B"], 1), 3) if np.isfinite(o["B"]) else np.nan,
                                 flags=o["flags"], n_runs=o["n_runs"], lat=r.lat_cm, lon=r.lon180 % 360, conf=r.conf))
                if r.tier == "le7":
                    runxy = None
                    if dbg is not None and np.isfinite(o.get("run_u0", np.nan)):
                        du, dw = dbg["du"], dbg["dw"]; c0 = np.array([dbg["cx"], dbg["cy"]])
                        runxy = (c0 + o["run_u0"] * du + o["run_w"] * dw, c0 + o["run_u1"] * du + o["run_w"] * dw)
                    bw = max(box[1] - box[0], box[3] - box[2]); h_ = rows[-1]["depth_m"]
                    G.tile(img, info, box, max(1.0 * bw, 25), f"{TD}/{int(r.hollow_id)}.png",
                           f"{int(r.hollow_id)} D={r.D_can_km*1000:.0f} m px={info['px']:.1f} i={inc:.0f} h={h_} ff={ff:.2f}", runxy, sbar=20 if info["px"] < 3 else 50)
            os.remove(p)
        if rows: done = pd.concat([done, pd.DataFrame(rows)]); done.to_csv(OUT, index=False); rows = []
    print("left", len(s) - (len(pd.read_csv(OUT)) if os.path.exists(OUT) else 0), "of", len(s), "| frames left", len(nums))
if __name__ == "__main__":
    run()
