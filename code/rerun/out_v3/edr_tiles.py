#!/usr/bin/env python3
"""Crater-bowl screen, step 3: per-detection EDR tiles, made in batches (fetch, tile, delete the EDR) because the VM disk
is small. Tiles go to qc/tiles/<hollow_id>.png. Resumable; run repeatedly until it reports 0 left. Then: edr_tiles.py sheets"""
import sys, os, time, numpy as np, pandas as pd
from PIL import Image, ImageDraw
HE = os.path.expanduser("~/mnt/hollows-emissivity"); OUT = f"{HE}/code/rerun/out_v3"; TD = f"{OUT}/qc/tiles"; os.makedirs(TD, exist_ok=True)
sys.path.insert(0, f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot")
import dep_common as C, validate as V

import re
def read_any(path):
    try: return C.read_edr(path)
    except AssertionError:
        raw = open(path, "rb").read(); end = raw.find(b"\r\nEND\r\n"); end = raw.find(b"\nEND\n") if end < 0 else end
        lab = raw[:end].decode("latin1"); kw = {}
        for m in re.finditer(r"^\s*([A-Z0-9_:\^]+)\s*=\s*(.+?)\s*$", lab, re.M): kw.setdefault(m.group(1), m.group(2))
        num = lambda k: float(re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", kw[k])[0]) if k in kw else np.nan
        rb = int(num("RECORD_BYTES")); ptr = int(num("^IMAGE")); obj = lab[lab.find("OBJECT = IMAGE"):]
        L = int(re.search(r"LINES\s*=\s*(\d+)", obj).group(1)); S = int(re.search(r"LINE_SAMPLES\s*=\s*(\d+)", obj).group(1))
        bits = int(re.search(r"SAMPLE_BITS\s*=\s*(\d+)", obj).group(1)); st = re.search(r"SAMPLE_TYPE\s*=\s*(\S+)", obj).group(1)
        dt = (">u2" if "MSB" in st else "<u2") if bits == 16 else np.uint8
        img = np.frombuffer(raw, dt, count=L * S, offset=(ptr - 1) * rb).reshape(L, S).astype(np.float32)
        return img, dict(lines=L, samples=S, subsol_az=num("SUB_SOLAR_AZIMUTH"), inc=num("INCIDENCE_ANGLE"), px=num("HORIZONTAL_PIXEL_SCALE"))

s = pd.read_csv(f"{OUT}/screen_set.csv", dtype={"num": str}); s["num"] = s.num.str.zfill(10); T = 240
def tile(r, img, info):
    H, W = img.shape; x0, x1, y0, y1 = r.x0 * W, r.x1 * W, r.y0 * H, r.y1 * H; bw = max(x1 - x0, y1 - y0); m = max(1.2 * bw, 40)
    X0, X1, Y0, Y1 = int(max(0, x0 - m)), int(min(W, x1 + m)), int(max(0, y0 - m)), int(min(H, y1 + m))
    sub = img[Y0:Y1, X0:X1].astype(float); lo_, hi_ = np.percentile(sub, [1, 99.5]); z = (T - 4) / max(sub.shape)
    im = Image.fromarray((np.clip((sub - lo_) / (hi_ - lo_ + 1e-6), 0, 1) * 255).astype(np.uint8)).resize((max(1, int(sub.shape[1] * z)), max(1, int(sub.shape[0] * z)))).convert("RGB")
    d = ImageDraw.Draw(im); d.rectangle([(x0 - X0) * z, (y0 - Y0) * z, (x1 - X0) * z, (y1 - Y0) * z], outline=(255, 200, 0), width=2)
    t = Image.new("RGB", (T, T + 26), (30, 30, 30)); t.paste(im, (0, 0)); dd = ImageDraw.Draw(t)
    a = np.radians(info["subsol_az"]); cx, cy = T - 16, 14; dd.line([cx, cy, cx + 11 * np.cos(a), cy + 11 * np.sin(a)], fill=(0, 255, 255), width=3)
    dd.text((2, T + 1), f"id {int(r.hollow_id)}  Deq {r.D_eq_km:.1f} km  set {r.set}", fill=(0, 255, 0))
    dd.text((2, T + 13), f"px {info['px']:.0f} m  i {info['inc']:.0f}  crop {sub.shape[1]*info['px']/1000:.1f} km", fill=(0, 255, 0))
    t.save(f"{TD}/{int(r.hollow_id)}.png")
if sys.argv[1] == "run":
    t0 = time.time()
    todo = s[[not os.path.exists(f"{TD}/{int(h)}.png") for h in s.hollow_id]]
    nums = list(todo.num.unique())
    while nums and time.time() - t0 < 140:
        chunk, nums = nums[:24], nums[24:]
        dts = todo.drop_duplicates("num").set_index("num").datetime
        V.fetch_many([(n, dts[n]) for n in chunk])
        for n in chunk:
            p = C.edr_path(n)
            if not os.path.exists(p): continue
            try: img, info = read_any(p)
            except Exception as e: print('read fail', n, e); os.remove(p); continue
            for _, r in todo[todo.num == n].iterrows():
                try: tile(r, img, info)
                except Exception as e: print("tile fail", r.hollow_id, e)
            os.remove(p)
    left = sum(not os.path.exists(f"{TD}/{int(h)}.png") for h in s.hollow_id); print("left", left, "of", len(s))
elif sys.argv[1] == "sheets":
    per = 48
    for which in ("A", "B", "C", "D"):
        q = s[s.set == which].reset_index(drop=True)
        for k in range(int(np.ceil(len(q) / per))):
            part = q.iloc[k * per:(k + 1) * per]; nc = 8; nr = int(np.ceil(len(part) / nc))
            S = Image.new("RGB", (nc * T, nr * (T + 26)), (60, 60, 60)); d = ImageDraw.Draw(S)
            for i, (_, r) in enumerate(part.iterrows()):
                f = f"{TD}/{int(r.hollow_id)}.png"; x, y = (i % nc) * T, (i // nc) * (T + 26)
                if os.path.exists(f): S.paste(Image.open(f), (x, y))
                else: d.text((x + 5, y + 100), f"no image id {int(r.hollow_id)}", fill=(255, 80, 80))
                d.text((x + 3, y + 3), f"{i + 1}", fill=(255, 60, 255))
            S.save(f"{OUT}/qc/edr_{which}_{k}.png"); part.assign(tile=np.arange(1, len(part) + 1))[["tile", "hollow_id", "lat", "lon", "D_eq_km", "set"]].to_csv(f"{OUT}/qc/edr_{which}_{k}_ids.csv", index=False)
            print(which, k, len(part))
