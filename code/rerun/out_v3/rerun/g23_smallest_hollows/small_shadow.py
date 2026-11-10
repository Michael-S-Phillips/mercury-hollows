#!/usr/bin/env python3
"""g23 (9 Nov 2026): shadow depths of the smallest cataloged hollows (D_eq < 100 m, the 21 hollows of the screened v3 catalog
used in the paper; all single detections in MDIS NAC frames at 1.4-6.9 m/pixel). For each: the automatic shadow length of the
depth-expansion pilot (shadow_measure.measure, v2), ground length L = L_px * px * f (f = foreshortening along the sun),
depth h = L / tan(i) with i the local incidence at the hollow. Also the depth that would just fill the box with shadow,
h_fill = (box extent along the sun, ground) / tan(i): a hollow fully in shadow is at least this deep, and a hollow whose far
floor is lit is shallower. Tiles: tiles/<id>_z.png (zoom, measured run in red) and _c.png (context); yellow = detection box,
cyan = toward the Sun. Resumable: `python3 small_shadow.py run` until it reports 0 left; then `python3 small_shadow.py sheets`."""
import sys, os, re, time, json, numpy as np, pandas as pd
from PIL import Image, ImageDraw
HE = os.path.expanduser("~/mnt/hollows-emissivity"); RR = f"{HE}/code/rerun/out_v3/rerun"; HERE = os.path.dirname(os.path.abspath(__file__))
TD = f"{HERE}/tiles"; os.makedirs(TD, exist_ok=True)
sys.path.insert(0, f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot")
import dep_common as C, validate as V, shadow_measure as M
DMAX = 0.1; T = 360
def read_any(path):
    try: return C.read_edr(path)
    except AssertionError:
        raw = open(path, "rb").read(); end = raw.find(b"\r\nEND\r\n"); end = raw.find(b"\nEND\n") if end < 0 else end
        lab = raw[:end].decode("latin1"); kw = {}
        for m in re.finditer(r"^\s*([A-Z0-9_:\^]+)\s*=\s*(.+?)\s*$", lab, re.M): kw.setdefault(m.group(1), m.group(2))
        num = lambda k: float(re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", kw[k])[0]) if k in kw else np.nan
        rb = int(num("RECORD_BYTES")); ptr = int(num("^IMAGE")); obj = lab[lab.find("OBJECT = IMAGE"):]
        L = int(re.search(r"LINES\s*=\s*(\d+)", obj).group(1)); S = int(re.search(r"LINE_SAMPLES\s*=\s*(\d+)", obj).group(1))
        st = re.search(r"SAMPLE_TYPE\s*=\s*(\S+)", obj).group(1)
        img = np.frombuffer(raw, ">u2" if "MSB" in st else "<u2", count=L * S, offset=(ptr - 1) * rb).reshape(L, S).astype(np.float32)
        info = dict(lines=L, samples=S, subsol_az=num("SUB_SOLAR_AZIMUTH"), subsc_az=num("SUB_SPACECRAFT_AZIMUTH"), inc=num("INCIDENCE_ANGLE"),
                    emi=num("EMISSION_ANGLE"), px=num("HORIZONTAL_PIXEL_SCALE"), subsol_lat=num("SUB_SOLAR_LATITUDE"),
                    subsol_lon=num("SUB_SOLAR_LONGITUDE"), dark=np.nan, start=kw.get("START_TIME", ""))
        return img, info
def select():
    H = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool) & (H.D_eq_km < DMAX)]
    ha = pd.read_csv(f"{RR}/g2_insolation/inputs/v3/host_assignments.csv")[["hollow_id", "host_uid", "rho"]]
    H = H.merge(ha, on="hollow_id", how="left")
    b = C.load_bickel(with_corners=False)[["thumb", "num", "datetime", "x0", "x1", "y0", "y1", "lat", "lon180"]].drop_duplicates("thumb")
    s = H.merge(b.rename(columns={"lat": "lat_cm"}), on="thumb", how="left"); s["num"] = s.num.astype(str).str.zfill(10)
    return s.sort_values("D_eq_km").reset_index(drop=True)
def tile(img, info, box, margin_px, path, label, run=None, sbar=None):
    H, W = img.shape; x0, x1, y0, y1 = box; cx, cy = (x0 + x1) / 2, (y0 + y1) / 2; half = max(x1 - x0, y1 - y0) / 2 + margin_px
    X0, X1, Y0, Y1 = int(max(0, cx - half)), int(min(W, cx + half)), int(max(0, cy - half)), int(min(H, cy + half))
    sub = img[Y0:Y1, X0:X1]; lo_, hi_ = np.percentile(sub, [0.5, 99.5]); z = T / max(sub.shape)
    im = Image.fromarray((np.clip((sub - lo_) / (hi_ - lo_ + 1e-6), 0, 1) * 255).astype(np.uint8)).resize(
        (max(1, int(sub.shape[1] * z)), max(1, int(sub.shape[0] * z))), Image.LANCZOS).convert("RGB")
    d = ImageDraw.Draw(im); d.rectangle([(x0 - X0) * z, (y0 - Y0) * z, (x1 - X0) * z, (y1 - Y0) * z], outline=(255, 200, 0), width=1)
    if run is not None: (ax, ay), (bx, by) = run; d.line([(ax - X0) * z, (ay - Y0) * z, (bx - X0) * z, (by - Y0) * z], fill=(255, 40, 40), width=2)
    t = Image.new("RGB", (T, T + 16), (20, 20, 20)); t.paste(im, ((T - im.width) // 2, (T - im.height) // 2)); dd = ImageDraw.Draw(t)
    a = np.radians(info["subsol_az"]); ax_, ay_ = T - 22, 18; dd.line([ax_, ay_, ax_ + 14 * np.cos(a), ay_ + 14 * np.sin(a)], fill=(0, 255, 255), width=3)
    dd.ellipse([ax_ + 14 * np.cos(a) - 3, ay_ + 14 * np.sin(a) - 3, ax_ + 14 * np.cos(a) + 3, ay_ + 14 * np.sin(a) + 3], fill=(0, 255, 255))
    if sbar:  # scale bar of sbar metres
        L = sbar / info["px"] * z; dd.line([8, T - 8, 8 + L, T - 8], fill=(0, 255, 0), width=3); dd.text((12 + L, T - 14), f"{sbar} m", fill=(0, 255, 0))
    dd.text((4, T + 2), label, fill=(230, 230, 230)); t.save(path)
def run():
    s = select(); out_csv = f"{HERE}/small_shadow.csv"; done = pd.read_csv(out_csv) if os.path.exists(out_csv) else pd.DataFrame(columns=["hollow_id"])
    todo = s[~s.hollow_id.isin(done.hollow_id)]; t0 = time.time(); rows = []
    nums = list(todo.num.unique()); added = []
    while nums and time.time() - t0 < 140:
        chunk, nums = nums[:8], nums[8:]; dts = todo.drop_duplicates("num").set_index("num").datetime
        V.fetch_many([(n, dts[n]) for n in chunk])
        for n in chunk:
            p = C.edr_path(n)
            if not os.path.exists(p): print("no EDR", n); continue
            img, info = read_any(p); Hh, Ww = img.shape; fs = M.foreshortening(info)
            for _, r in todo[todo.num == n].iterrows():
                box = (r.x0 * Ww, r.x1 * Ww, r.y0 * Hh, r.y1 * Hh)
                o = M.measure(img, info, box, return_debug=True); dbg = o.pop("_dbg", None)
                inc = float(C.local_incidence(r.lat_cm, r.lon180, info["subsol_lat"], info["subsol_lon"])); ti = np.tan(np.radians(inc))
                Lm = o["L_px"] * info["px"] * fs if np.isfinite(o["L_px"]) else np.nan
                along_m = o["along_px"] * info["px"] * fs; across_m = o["across_px"] * info["px"]
                runxy = None
                if dbg is not None and np.isfinite(o.get("run_u0", np.nan)):
                    du, dw = dbg["du"], dbg["dw"]; c0 = np.array([dbg["cx"], dbg["cy"]])
                    runxy = (c0 + o["run_u0"] * du + o["run_w"] * dw, c0 + o["run_u1"] * du + o["run_w"] * dw)
                bw = max(box[1] - box[0], box[3] - box[2])
                tile(img, info, box, max(1.5 * bw, 25), f"{TD}/{int(r.hollow_id)}_z.png",
                     f"{int(r.hollow_id)} zoom  D={r.D_eq_km*1000:.0f} m  px={info['px']:.2f} m  i={inc:.1f}", runxy, sbar=20 if info["px"] < 3 else 50)
                tile(img, info, box, max(300 / info["px"], 4 * bw), f"{TD}/{int(r.hollow_id)}_c.png",
                     f"{int(r.hollow_id)} context  {r.lat:.2f}N {r.lon:.2f}E", None, sbar=100 if info["px"] < 3 else 250)
                rows.append(dict(hollow_id=int(r.hollow_id), D_eq_m=round(r.D_eq_km * 1000, 1), maxext_m=round(r.max_ext_km * 1000, 1), num=n,
                                 px=info["px"], inc_label=info["inc"], inc=round(inc, 2), emi=info["emi"], f=round(fs, 3), subsol_az=info["subsol_az"],
                                 B=o["B"], S=o["S"], contrast=round((o["B"] - o["S"]) / max(o["B"], 1), 3) if np.isfinite(o["B"]) else np.nan,
                                 L_px=o["L_px"], along_px=round(o["along_px"], 1), L_m=round(Lm, 1) if np.isfinite(Lm) else np.nan,
                                 along_m=round(along_m, 1), across_m=round(across_m, 1), shadow_frac=round(o["L_px"] / o["along_px"], 2) if o["along_px"] > 0 else np.nan,
                                 depth_m=round(Lm / ti, 1) if np.isfinite(Lm) else np.nan, h_fill_m=round(along_m / ti, 1), flags=o["flags"],
                                 n_runs=o["n_runs"], host_uid=r.host_uid, rho=r.rho, lat=r.lat, lon=r.lon, conf=r.conf))
            added.append(p)
    if rows: pd.concat([done, pd.DataFrame(rows)]).to_csv(out_csv, index=False)
    for p in added:
        if os.path.exists(p): os.remove(p)
    print("left", len(s) - len(pd.read_csv(out_csv)), "of", len(s))
def sheets():
    r = pd.read_csv(f"{HERE}/small_shadow.csv").sort_values("D_eq_m"); ids = list(r.hollow_id)
    for k in range(0, len(ids), 6):
        sh = Image.new("RGB", (4 * T, 3 * (T + 16)), (0, 0, 0))
        for j, h in enumerate(ids[k:k + 6]):
            for q, sfx in enumerate("zc"):
                p = f"{TD}/{h}_{sfx}.png"
                if os.path.exists(p): sh.paste(Image.open(p), ((j % 2) * 2 * T + q * T, (j // 2) * (T + 16)))
        sh.save(f"{HERE}/sheet_{k // 6 + 1}.png")
    print("sheets", (len(ids) + 5) // 6)
if __name__ == "__main__":
    {"run": run, "sheets": sheets}[sys.argv[1]]()
