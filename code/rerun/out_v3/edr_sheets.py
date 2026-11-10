#!/usr/bin/env python3
"""Crater-bowl screen, step 2 (6 Nov 2026): EDR crops of candidate spurious detections for visual classification.
Sets: A = rule A (D_eq >= 0.35 D, centroid within 0.5 R of a Herrick crater); B = other detections with D_eq >= 1.5 km in
fields of <= 3 hollows; C = D_eq 1.0-1.5 km singletons whose host has <= 2 hollows within 0.9 radii; D = control, 48 random
detections with D_eq >= 1.5 km in fields of >= 10 hollows. Usage: edr_sheets.py select | fetch | render"""
import sys, os, time, numpy as np, pandas as pd
HE = os.path.expanduser("~/mnt/hollows-emissivity"); OUT = f"{HE}/code/rerun/out_v3"
sys.path.insert(0, f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot")
import dep_common as C
step = sys.argv[1]
if step == "select":
    c = pd.read_csv(f"{OUT}/bowl_candidates.csv"); c["field_n"] = c.field_id.map(c.groupby("field_id").size())
    ha = pd.read_csv(f"{HE}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", usecols=["hollow_id", "host_uid", "rho"])
    n09 = ha[ha.rho <= 0.9].groupby("host_uid").size(); ha["host_n09"] = ha.host_uid.map(n09).fillna(0)
    c = c.merge(ha[["hollow_id", "host_uid", "rho", "host_n09"]].rename(columns={"rho": "host_rho"}), on="hollow_id", how="left")
    A = c.bowl_cr >= 0
    B = (~A) & (c.D_eq_km >= 1.5) & (c.field_n <= 3)
    Cc = (~A) & (~B) & c.D_eq_km.between(1.0, 1.5) & (c.field_n == 1) & (c.host_n09 <= 2)
    D = (~A) & (c.D_eq_km >= 1.5) & (c.field_n >= 10)
    c["set"] = np.where(A, "A", np.where(B, "B", np.where(Cc, "C", "")))
    d = c[D].sample(48, random_state=7).index; c.loc[d, "set"] = "D"
    s = c[c.set != ""].copy()
    b = C.load_bickel(with_corners=False)[["thumb", "num", "datetime", "x0", "x1", "y0", "y1", "px_rel", "incidence_rel"]]
    s = s.merge(b, on="thumb", how="left"); s = s.sort_values(["set", "hollow_id"]).reset_index(drop=True)
    s.to_csv(f"{OUT}/screen_set.csv", index=False); print(s.set.value_counts().to_dict(), "images", s.num.nunique())
elif step == "fetch":
    s = pd.read_csv(f"{OUT}/screen_set.csv", dtype={"num": str}); s["num"] = s.num.str.zfill(10)
    sys.path.insert(0, f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/depth_expansion_pilot"); import validate as V
    t0 = time.time(); todo = [(n, d) for n, d in s.drop_duplicates("num")[["num", "datetime"]].itertuples(index=False) if not os.path.exists(C.edr_path(n))]
    while todo and time.time() - t0 < 150:
        chunk, todo = todo[:16], todo[16:]; V.fetch_many(chunk)
    print("remaining", len([1 for n in s.num.unique() if not os.path.exists(C.edr_path(n))]), "of", s.num.nunique())
elif step == "render":
    from PIL import Image, ImageDraw
    s = pd.read_csv(f"{OUT}/screen_set.csv", dtype={"num": str}); s["num"] = s.num.str.zfill(10)
    which = sys.argv[2]; k0 = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    q = s[s.set == which].reset_index(drop=True); T = 240; per = 48
    os.makedirs(f"{OUT}/qc", exist_ok=True)
    for k in range(k0, int(np.ceil(len(q) / per))):
        part = q.iloc[k * per:(k + 1) * per]; tiles = []
        for _, r in part.iterrows():
            p = C.edr_path(r.num)
            if not os.path.exists(p): tiles.append(None); continue
            img, info = C.read_edr(p); H, W = img.shape
            x0, x1, y0, y1 = r.x0 * W, r.x1 * W, r.y0 * H, r.y1 * H; bw = max(x1 - x0, y1 - y0); m = max(1.2 * bw, 40)
            X0, X1, Y0, Y1 = int(max(0, x0 - m)), int(min(W, x1 + m)), int(max(0, y0 - m)), int(min(H, y1 + m))
            sub = img[Y0:Y1, X0:X1].astype(float); lo_, hi_ = np.percentile(sub, [1, 99.5])
            z = (T - 4) / max(sub.shape)
            im = Image.fromarray((np.clip((sub - lo_) / (hi_ - lo_ + 1e-6), 0, 1) * 255).astype(np.uint8)).resize((max(1, int(sub.shape[1] * z)), max(1, int(sub.shape[0] * z)))).convert("RGB")
            d = ImageDraw.Draw(im); d.rectangle([(x0 - X0) * z, (y0 - Y0) * z, (x1 - X0) * z, (y1 - Y0) * z], outline=(255, 200, 0), width=2)
            t = Image.new("RGB", (T, T + 26)); t.paste(im, (0, 0)); dd = ImageDraw.Draw(t)
            a = np.radians(info["subsol_az"]); cx, cy = T - 16, 14; dd.line([cx, cy, cx + 11 * np.cos(a), cy + 11 * np.sin(a)], fill=(0, 255, 255), width=3)
            dd.text((2, T + 1), f"#{_ + k * 0} id {int(r.hollow_id)} Deq {r.D_eq_km:.1f} km", fill=(0, 255, 0))
            dd.text((2, T + 13), f"px {info['px']:.0f} m  i {info['inc']:.0f}  crop {sub.shape[1]*info['px']/1000:.1f} km", fill=(0, 255, 0))
            tiles.append(t)
        nc = 8; nr = int(np.ceil(len(tiles) / nc)); S = Image.new("RGB", (nc * T, nr * (T + 26)), (40, 40, 40))
        for i, t in enumerate(tiles):
            if t is not None: S.paste(t, ((i % nc) * T, (i // nc) * (T + 26)))
        S.save(f"{OUT}/qc/edr_{which}_{k}.png"); part[["hollow_id", "lat", "lon", "D_eq_km", "num", "set"]].to_csv(f"{OUT}/qc/edr_{which}_{k}_ids.csv", index=True)
        print(which, k, len(part))
