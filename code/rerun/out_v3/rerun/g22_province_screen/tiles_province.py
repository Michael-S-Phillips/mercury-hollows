#!/usr/bin/env python3
"""g22 step 2: for each hollow in province_screen_set.csv, two tiles from the EDR of its representative (finest) detection:
a zoom (detection box plus 1.2 box widths, at least 40 pixels, as in the v3 screen's edr_tiles.py) and a context view about 8 km
across (at least 6 box widths). Yellow = detection box; cyan = toward the Sun. EDRs are fetched, tiled, and deleted in batches
(small VM disk). Resumable: run until it reports 0 left. Tiles: tiles/<hollow_id>_z.png and _c.png (raw EDR DN, 1-99.5% stretch)."""
import sys, os, re, time, numpy as np, pandas as pd
from PIL import Image, ImageDraw
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__)); TD = f"{HERE}/tiles"; os.makedirs(TD, exist_ok=True)
sys.path.insert(0, f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot")
import dep_common as C, validate as V
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
s = pd.read_csv(f"{HERE}/province_screen_set.csv", dtype={"num": str}); s["num"] = s.num.str.zfill(10); T = 300
def crop(img, info, r, margin_px, path, label):
    H, W = img.shape; x0, x1, y0, y1 = r.x0 * W, r.x1 * W, r.y0 * H, r.y1 * H; cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    half = max(x1 - x0, y1 - y0) / 2 + margin_px
    X0, X1, Y0, Y1 = int(max(0, cx - half)), int(min(W, cx + half)), int(max(0, cy - half)), int(min(H, cy + half))
    sub = img[Y0:Y1, X0:X1].astype(float); lo_, hi_ = np.percentile(sub, [1, 99.5]); z = T / max(sub.shape)
    im = Image.fromarray((np.clip((sub - lo_) / (hi_ - lo_ + 1e-6), 0, 1) * 255).astype(np.uint8)).resize((max(1, int(sub.shape[1] * z)), max(1, int(sub.shape[0] * z))), Image.LANCZOS).convert("RGB")
    d = ImageDraw.Draw(im); d.rectangle([(x0 - X0) * z, (y0 - Y0) * z, (x1 - X0) * z, (y1 - Y0) * z], outline=(255, 200, 0), width=2)
    t = Image.new("RGB", (T, T), (20, 20, 20)); t.paste(im, ((T - im.width) // 2, (T - im.height) // 2)); dd = ImageDraw.Draw(t)
    a = np.radians(info["subsol_az"]); ax, ay = T - 18, 16; dd.line([ax, ay, ax + 12 * np.cos(a), ay + 12 * np.sin(a)], fill=(0, 255, 255), width=3)
    dd.text((4, T - 14), f"{label}  {sub.shape[1] * info['px'] / 1000:.1f} km across", fill=(0, 255, 0)); t.save(path)
    return dict(px=info["px"], inc=info["inc"])
meta = []
if sys.argv[1] == "run":
    t0 = time.time(); todo = s[[not os.path.exists(f"{TD}/{int(h)}_c.png") for h in s.hollow_id]]; nums = list(todo.num.unique())
    while nums and time.time() - t0 < 140:
        chunk, nums = nums[:12], nums[12:]
        dts = todo.drop_duplicates("num").set_index("num").datetime; V.fetch_many([(n, dts[n]) for n in chunk])
        for n in chunk:
            p = C.edr_path(n)
            if not os.path.exists(p): print("no EDR", n); continue
            try: img, info = read_any(p)
            except Exception as e: print("read fail", n, e); os.remove(p); continue
            for _, r in todo[todo.num == n].iterrows():
                try:
                    bw = max(r.x1 - r.x0, r.y1 - r.y0) * img.shape[1]
                    crop(img, info, r, max(1.2 * bw, 40), f"{TD}/{int(r.hollow_id)}_z.png", "zoom")
                    crop(img, info, r, max(4000 / info["px"], 3 * bw), f"{TD}/{int(r.hollow_id)}_c.png", "context")
                except Exception as e: print("tile fail", r.hollow_id, e)
            os.remove(p)
    print("left", sum(not os.path.exists(f"{TD}/{int(h)}_c.png") for h in s.hollow_id), "of", len(s))
