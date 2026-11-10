#!/usr/bin/env python3
"""g23 step 6: images and manifest for the "Hollow Shadow Bench" artifact, where Michael measures by hand the main shadow of each of the
135 hollows in the <= 7 m/pixel tier of depth_vs_size.csv. One PNG per hollow (tiles_meas/<id>.png): left, the raw EDR crop at native
resolution (box plus max(1 box size, 25 px) on each side, 0.5-99.5% stretch, nothing drawn on it); right, a 300 px context view
(at least 8 box sizes and 0.8 km across) with the detection box in yellow. manifest_meas.json carries, per hollow, the crop geometry,
the box and the automatic shadow line in crop pixels, and px, f (foreshortening along the sun), local incidence and sun azimuth, so
the page converts a clicked shadow into depth = L_px x px x f / tan(i). Resumable: run until it reports 0 left."""
import sys, os, json, time, numpy as np, pandas as pd
from PIL import Image, ImageDraw
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import small_shadow as G, depth_vs_size as DV
C, V, M = G.C, G.V, G.M
TD = f"{HERE}/tiles_meas"; MF = f"{HERE}/manifest_meas.json"; CW = 300
s = DV.select(); s = s[s.tier == "le7"].reset_index(drop=True)
man = json.load(open(MF)) if os.path.exists(MF) else {}
todo = s[[str(int(h)) not in man for h in s.hollow_id]]; nums = list(todo.num.unique()); t0 = time.time()
dts = todo.drop_duplicates("num").set_index("num").datetime
def stretch(a):
    lo, hi = np.percentile(a, [0.5, 99.5]); return (np.clip((a - lo) / (hi - lo + 1e-6), 0, 1) * 255).astype(np.uint8)
while nums and time.time() - t0 < 100:
    chunk, nums = nums[:10], nums[10:]; V.fetch_many([(n, dts[n]) for n in chunk])
    for n in chunk:
        p = C.edr_path(n)
        if not os.path.exists(p): print("no EDR", n); continue
        img, info = G.read_any(p); H, W = img.shape; fs = M.foreshortening(info)
        for _, r in todo[todo.num == n].iterrows():
            hid = int(r.hollow_id); x0, x1, y0, y1 = r.x0 * W, r.x1 * W, r.y0 * H, r.y1 * H; bw = max(x1 - x0, y1 - y0); cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            half = bw / 2 + max(bw, 25); X0, X1, Y0, Y1 = int(max(0, cx - half)), int(min(W, cx + half)), int(max(0, cy - half)), int(min(H, cy + half))
            z = stretch(img[Y0:Y1, X0:X1]); zh, zw = z.shape
            hc = max(4 * bw, 400 / info["px"]); A0, A1, B0, B1 = int(max(0, cx - hc)), int(min(W, cx + hc)), int(max(0, cy - hc)), int(min(H, cy + hc))
            csub = stretch(img[B0:B1, A0:A1]); sc = CW / max(csub.shape)
            cim = Image.fromarray(csub).resize((max(1, int(csub.shape[1] * sc)), max(1, int(csub.shape[0] * sc))), Image.LANCZOS).convert("RGB")
            ImageDraw.Draw(cim).rectangle([(x0 - A0) * sc, (y0 - B0) * sc, (x1 - A0) * sc, (y1 - B0) * sc], outline=(242, 194, 48), width=1)
            Hc = max(zh, cim.height); comp = Image.new("RGB", (zw + 8 + cim.width, Hc), (17, 20, 23))
            comp.paste(Image.fromarray(z).convert("RGB"), (0, 0)); comp.paste(cim, (zw + 8, 0)); comp.save(f"{TD}/{hid}.png", optimize=True)
            o = M.measure(img, info, (x0, x1, y0, y1), return_debug=True); dbg = o.pop("_dbg", None); auto = None
            if dbg is not None and np.isfinite(o.get("run_u0", np.nan)) and o["L_px"] > 0:
                du, dw = dbg["du"], dbg["dw"]; c0 = np.array([dbg["cx"], dbg["cy"]])
                a = c0 + o["run_u0"] * du + o["run_w"] * dw; b = c0 + o["run_u1"] * du + o["run_w"] * dw
                auto = [round(a[0] - X0, 2), round(a[1] - Y0, 2), round(b[0] - X0, 2), round(b[1] - Y0, 2)]
            inc = float(C.local_incidence(r.lat_cm, r.lon180, info["subsol_lat"], info["subsol_lon"]))
            Lm = o["L_px"] * info["px"] * fs if np.isfinite(o["L_px"]) else np.nan
            man[str(hid)] = dict(id=hid, D=round(r.D_can_km * 1000, 1), px=round(info["px"], 4), f=round(fs, 4), inc=round(inc, 3), az=round(info["subsol_az"], 3),
                                 emi=round(info["emi"], 2), img=f"EN{n}M", lat=round(r.lat_cm, 3), lon=round(r.lon180 % 360, 3), zw=zw, zh=zh, X0=X0, Y0=Y0,
                                 box=[round(x0 - X0, 2), round(x1 - X0, 2), round(y0 - Y0, 2), round(y1 - Y0, 2)], cw=cim.width, ch=cim.height,
                                 ctx_km=round((A1 - A0) * info["px"] / 1000, 2), auto=auto,
                                 auto_h=round(Lm / np.tan(np.radians(inc)), 2) if np.isfinite(Lm) and Lm > 0 else None, auto_flags=o["flags"] or "")
        os.remove(p)
    json.dump(man, open(MF, "w"))
print("left", len(s) - len(man), "of", len(s))
