#!/usr/bin/env python3
"""g23 step 2: the same measurement with the detection box widened toward the Sun by 0.75 box lengths along the sun direction,
so that a shadow cast inside a depression whose sunward part the detection box left out (e.g. 2496) is counted. Also writes,
per hollow, the mean brightness profile along the sun direction through the 5 lines around the longest run (profiles.json),
and zoom tiles with both boxes (tiles/<id>_v2.png: yellow = detection box, magenta = widened box, red = measured run).
Run: python3 small_shadow_v2.py (re-fetches the EDRs and deletes them afterwards)."""
import sys, os, json, numpy as np, pandas as pd
from PIL import Image, ImageDraw
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import small_shadow as G
C, V, M = G.C, G.V, G.M
s = G.select(); r1 = pd.read_csv(f"{HERE}/small_shadow.csv").set_index("hollow_id"); rows = []; prof = {}
nums = list(s.num.unique()); V.fetch_many([(n, s.drop_duplicates("num").set_index("num").datetime[n]) for n in nums])
for n in nums:
    p = C.edr_path(n); img, info = G.read_any(p); Hh, Ww = img.shape; fs = M.foreshortening(info)
    for _, r in s[s.num == n].iterrows():
        x0, x1, y0, y1 = r.x0 * Ww, r.x1 * Ww, r.y0 * Hh, r.y1 * Hh; az = np.radians(info["subsol_az"]); su = np.array([np.cos(az), np.sin(az)])
        ext = 0.75 * max(x1 - x0, y1 - y0)   # widen toward the sun
        bx = (min(x0, x0 + ext * su[0]), max(x1, x1 + ext * su[0]), min(y0, y0 + ext * su[1]), max(y1, y1 + ext * su[1]))
        o = M.measure(img, info, bx, return_debug=True); dbg = o.pop("_dbg", None)
        inc = r1.loc[r.hollow_id, "inc"]; ti = np.tan(np.radians(inc)); Lm = o["L_px"] * info["px"] * fs if np.isfinite(o["L_px"]) else np.nan
        runxy = None
        if dbg is not None and np.isfinite(o.get("run_u0", np.nan)):
            du, dw = dbg["du"], dbg["dw"]; c0 = np.array([dbg["cx"], dbg["cy"]])
            runxy = (c0 + o["run_u0"] * du + o["run_w"] * dw, c0 + o["run_u1"] * du + o["run_w"] * dw)
            j = int(np.argmin(np.abs(dbg["W"] - o["run_w"]))); A = dbg["A"][max(0, j - 2):j + 3].mean(0)
            prof[int(r.hollow_id)] = dict(u_px=dbg["U"][::4].round(2).tolist(), dn=A[::4].round(1).tolist(), B=o["B"], S=o["S"], T=o["T"], u0=o["run_u0"], u1=o["run_u1"])
        # tile with both boxes
        bw = max(bx[1] - bx[0], bx[3] - bx[2]); cx, cy = (bx[0] + bx[1]) / 2, (bx[2] + bx[3]) / 2; half = bw / 2 + max(bw, 20)
        X0, X1, Y0, Y1 = int(max(0, cx - half)), int(min(Ww, cx + half)), int(max(0, cy - half)), int(min(Hh, cy + half))
        sub = img[Y0:Y1, X0:X1]; lo_, hi_ = np.percentile(sub, [0.5, 99.5]); T = 420; z = T / max(sub.shape)
        im = Image.fromarray((np.clip((sub - lo_) / (hi_ - lo_ + 1e-6), 0, 1) * 255).astype(np.uint8)).resize((int(sub.shape[1] * z), int(sub.shape[0] * z)), Image.LANCZOS).convert("RGB")
        d = ImageDraw.Draw(im)
        d.rectangle([(bx[0] - X0) * z, (bx[2] - Y0) * z, (bx[1] - X0) * z, (bx[3] - Y0) * z], outline=(255, 0, 255), width=1)
        d.rectangle([(x0 - X0) * z, (y0 - Y0) * z, (x1 - X0) * z, (y1 - Y0) * z], outline=(255, 200, 0), width=1)
        if runxy is not None: (ax, ay), (bx_, by_) = runxy; d.line([(ax - X0) * z, (ay - Y0) * z, (bx_ - X0) * z, (by_ - Y0) * z], fill=(255, 40, 40), width=2)
        t = Image.new("RGB", (T, T + 18), (20, 20, 20)); t.paste(im, ((T - im.width) // 2, (T - im.height) // 2)); dd = ImageDraw.Draw(t)
        dd.line([T - 26, 18, T - 26 + 16 * su[0], 18 + 16 * su[1]], fill=(0, 255, 255), width=3)
        sb = 20 if info["px"] < 3 else 50; L = sb / info["px"] * z; dd.line([8, T - 8, 8 + L, T - 8], fill=(0, 255, 0), width=3); dd.text((12 + L, T - 14), f"{sb} m", fill=(0, 255, 0))
        dd.text((4, T + 3), f"{int(r.hollow_id)}  D={r.D_eq_km*1000:.0f} m  i={inc:.1f}  v1 h={r1.loc[r.hollow_id,'depth_m']} m  v2 h={Lm/ti:.1f} m", fill=(230, 230, 230))
        t.save(f"{HERE}/tiles/{int(r.hollow_id)}_v2.png")
        rows.append(dict(hollow_id=int(r.hollow_id), L_px_v2=o["L_px"], L_m_v2=round(Lm, 1) if np.isfinite(Lm) else np.nan,
                         depth_m_v2=round(Lm / ti, 1) if np.isfinite(Lm) else np.nan, flags_v2=o["flags"], contrast_v2=round((o["B"] - o["S"]) / max(o["B"], 1), 3) if np.isfinite(o["B"]) else np.nan,
                         run_starts_outside_box=bool(np.isfinite(o.get("run_u0", np.nan)) and runxy is not None and
                                                     ((runxy[0][0] - x0) * -su[0] + 0) > 0.5 and not (x0 - 0.5 <= runxy[0][0] <= x1 + 0.5))))
    os.remove(p)
v2 = pd.DataFrame(rows); out = r1.reset_index().merge(v2, on="hollow_id"); out.to_csv(f"{HERE}/small_shadow_v2.csv", index=False)
json.dump(prof, open(f"{HERE}/profiles.json", "w"))
pd.set_option("display.width", 250)
print(out[["hollow_id", "D_eq_m", "px", "inc", "L_px", "depth_m", "flags", "L_px_v2", "depth_m_v2", "flags_v2", "contrast_v2", "run_starts_outside_box", "h_fill_m"]].sort_values("D_eq_m").to_string(index=False))
# sheets of v2 tiles, 3 x 3
ids = list(out.sort_values("D_eq_m").hollow_id)
for k in range(0, len(ids), 9):
    sh = Image.new("RGB", (3 * 420, 3 * 438), (0, 0, 0))
    for j, h in enumerate(ids[k:k + 9]): sh.paste(Image.open(f"{HERE}/tiles/{h}_v2.png"), ((j % 3) * 420, (j // 3) * 438))
    sh.save(f"{HERE}/sheet_v2_{k // 9 + 1}.png")
