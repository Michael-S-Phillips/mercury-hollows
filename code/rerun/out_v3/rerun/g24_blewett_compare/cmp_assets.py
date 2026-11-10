#!/usr/bin/env python3
"""g24 (9 Nov 2026): blind comparison of Michael's hand shadow depths with Blewett et al. (2016) on the same hollows.
Source: the depth pilot's registration (analysis_workshop/depth_expansion_pilot/registration_points.csv): Blewett points shifted into the
raw EDR frames of the images they were measured in, in the 130 accepted frames, and matched to Bickel detection boxes (fid >= 0).
Set: every box poleward of 60N (28 boxes, 14 frames) and a random 72 of the 441 equatorward boxes, at most 3 per frame (seed 24).
Per box: tiles_cmp/<fid>.png (native crop, box plus max(1 box size, 25 px) each side, and a 300 px context view), manifest_cmp.json for the
page (crop geometry, box, Blewett point positions, px, f, local incidence, sun azimuth; no Blewett depths and no latitude), and
key_cmp.csv (Blewett depths, latitude and the automatic depth, kept off the page). Resumable."""
import sys, os, json, time, numpy as np, pandas as pd
from PIL import Image, ImageDraw
HE = os.path.expanduser("~/mnt/hollows-emissivity"); P = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot"
sys.path.insert(0, P); import dep_common as C, validate as V, shadow_measure as M
sys.path.insert(0, f"{HE}/code/rerun/out_v3/rerun/g23_smallest_hollows"); import small_shadow as G
HERE = os.path.dirname(os.path.abspath(__file__)); TD = f"{HERE}/tiles_cmp"; MF = f"{HERE}/manifest_cmp.json"; KF = f"{HERE}/key_cmp.csv"; CW = 300
fr = pd.read_csv(f"{P}/registration_frames.csv", dtype={"num": str}); pt = pd.read_csv(f"{P}/registration_points.csv", dtype={"num": str})
bl = C.load_blewett(); k = C.load_bickel(with_corners=False)
p = pt[pt.num.isin(fr[fr.accepted].num) & (pt.fid >= 0)].merge(bl[["bid", "latitude", "longitude", "depth", "px", "inc"]].rename(columns={"px": "b_px", "inc": "b_inc", "depth": "b_depth"}), on="bid")
bx = p.groupby("fid").agg(num=("num", "first"), b_lat=("latitude", "median"), n_b=("bid", "size"), b_max=("b_depth", "max"), b_med=("b_depth", "median"),
                          b_px=("b_px", "median"), b_inc=("b_inc", "median")).reset_index()
bx["pole"] = bx.b_lat >= 60
eq = bx[~bx.pole].sample(frac=1, random_state=24).groupby("num").head(3).head(72)
sel = pd.concat([bx[bx.pole], eq]).merge(k[["fid", "datetime", "x0", "x1", "y0", "y1", "lat", "lon180", "poly_eqd_km"]], on="fid")
val = pd.read_csv(f"{P}/validation.csv")[["fid", "depth_m", "flags"]].rename(columns={"depth_m": "auto_depth_pilot", "flags": "auto_flags_pilot"})
sel = sel.merge(val, on="fid", how="left"); sel.to_csv(KF, index=False)
man = json.load(open(MF)) if os.path.exists(MF) else {}
todo = sel[[str(int(f)) not in man for f in sel.fid]]; nums = list(todo.num.unique()); t0 = time.time(); dts = todo.drop_duplicates("num").set_index("num").datetime
def stretch(a): lo, hi = np.percentile(a, [0.5, 99.5]); return (np.clip((a - lo) / (hi - lo + 1e-6), 0, 1) * 255).astype(np.uint8)
while nums and time.time() - t0 < 100:
    chunk, nums = nums[:10], nums[10:]; V.fetch_many([(n, dts[n]) for n in chunk])
    for n in chunk:
        pth = C.edr_path(n)
        if not os.path.exists(pth): print("no EDR", n); continue
        img, info = G.read_any(pth); H, W = img.shape; fs = M.foreshortening(info)
        for _, r in todo[todo.num == n].iterrows():
            fid = int(r.fid); x0, x1, y0, y1 = r.x0 * W, r.x1 * W, r.y0 * H, r.y1 * H; bw = max(x1 - x0, y1 - y0); cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            half = bw / 2 + max(bw, 25); X0, X1, Y0, Y1 = int(max(0, cx - half)), int(min(W, cx + half)), int(max(0, cy - half)), int(min(H, cy + half))
            z = stretch(img[Y0:Y1, X0:X1]); zh, zw = z.shape
            hc = max(4 * bw, 400 / info["px"]); A0, A1, B0, B1 = int(max(0, cx - hc)), int(min(W, cx + hc)), int(max(0, cy - hc)), int(min(H, cy + hc))
            csub = stretch(img[B0:B1, A0:A1]); sc = CW / max(csub.shape)
            cim = Image.fromarray(csub).resize((max(1, int(csub.shape[1] * sc)), max(1, int(csub.shape[0] * sc))), Image.LANCZOS).convert("RGB")
            ImageDraw.Draw(cim).rectangle([(x0 - A0) * sc, (y0 - B0) * sc, (x1 - A0) * sc, (y1 - B0) * sc], outline=(242, 194, 48), width=1)
            comp = Image.new("RGB", (zw + 8 + cim.width, max(zh, cim.height)), (17, 20, 23)); comp.paste(Image.fromarray(z).convert("RGB"), (0, 0)); comp.paste(cim, (zw + 8, 0))
            comp.save(f"{TD}/{fid}.png", optimize=True)
            mk = p[p.fid == fid][["u", "v"]].values; inc = float(C.local_incidence(r.lat, r.lon180, info["subsol_lat"], info["subsol_lon"]))
            man[str(fid)] = dict(id=fid, px=round(info["px"], 4), f=round(fs, 4), inc=round(inc, 3), az=round(info["subsol_az"], 3), emi=round(info["emi"], 2),
                                 img=f"EN{n}M", zw=zw, zh=zh, X0=X0, Y0=Y0, box=[round(x0 - X0, 2), round(x1 - X0, 2), round(y0 - Y0, 2), round(y1 - Y0, 2)],
                                 marks=[[round(u - X0, 1), round(v - Y0, 1)] for u, v in mk], cw=cim.width, ch=cim.height, ctx_km=round((A1 - A0) * info["px"] / 1000, 2),
                                 D=round(r.poly_eqd_km * 1000, 0) if np.isfinite(r.poly_eqd_km) else None)
        os.remove(pth)
    json.dump(man, open(MF, "w"))
print("left", len(sel) - len(man), "of", len(sel), "| poleward", int(sel.pole.sum()))
