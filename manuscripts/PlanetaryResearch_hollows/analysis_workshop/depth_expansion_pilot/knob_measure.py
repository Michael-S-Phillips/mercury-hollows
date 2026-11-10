#!/usr/bin/env python3
"""Shadow depths of hollows on circum-Caloris knobs (5 Nov 2026), with the pilot's automatic v2 method.
Targets: Bickel et al. (2024) detections that lie on a knob of Wright et al. (2020) (within the knob's half long axis, or
3.5 km for unmeasured certain/uncertain knobs and 2.0 km for knobs < 5 km, plus 0.5 km), in frames with incidence >= 55 deg
(list: ../province_basins/wright2020_knobs/knob_bickel_detections.csv). Frames come from the PDS Imaging Node (dep_common).
Writes knob_measurements.csv, knob_measure_record.json and qc/knob_sheet.png (close and context tiles per detection)."""
import sys, os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from PIL import Image, ImageDraw
import dep_common as C, shadow_measure as M, validate as V, qc_render as Q
HERE = C.HERE; PB = f"{HERE}/../province_basins/wright2020_knobs"
t = pd.read_csv(f"{PB}/knob_bickel_detections.csv", dtype={"num": str})
t = t[t.incidence_rel >= 55].copy(); t["num"] = t.num.str.zfill(10)
print("targets", len(t), "frames", t.num.nunique(), "knobs", t.knob.nunique())
print("have", V.fetch_many(list(t.drop_duplicates("num")[["num", "datetime"]].itertuples(index=False, name=None))), "of", t.num.nunique())
out = f"{HERE}/knob_measurements.csv"
# measure_frames resumes: frames already in the output are skipped
k = C.load_bickel(with_corners=False); k = k[k.fid.isin(t.fid)]
V.measure_frames(list(t.num.unique()), k, out, tmax=1e9)
m = pd.read_csv(out, dtype={"num": str}).merge(t[["fid", "knob", "ktype", "kaxis", "shadow_px_25m"]], on="fid")
# correction to the Blewett median pixel scale and incidence, with the paper's OLS fit on the 2,518 Blewett depths
bl = pd.read_csv(C.F_BLEWETT); bl["lntan"] = np.log(np.tan(np.radians(bl.inc)))
f = smf.ols("lnh ~ lnpx + lntan", bl).fit(); bpx, bti = f.params["lnpx"], f.params["lntan"]
px0, i0 = np.exp(np.median(bl.lnpx)), np.median(bl.inc)
corr = lambda d: np.exp(np.log(d.depth_m) - bpx * (np.log(d.px_m) - np.log(px0)) - bti * (np.log(np.tan(np.radians(d.inc_local))) - np.log(np.tan(np.radians(i0)))))
m["depth_corr_m"] = corr(m)
base = pd.read_csv(f"{HERE}/base_measurements.csv", dtype={"num": str})
base = base[np.isfinite(base.depth_m) & (base["flags"].fillna("") == "")]; base["depth_corr_m"] = corr(base)
clean = m[np.isfinite(m.depth_m) & (m["flags"].fillna("") == "")]
R = {"targets": int(len(t)), "frames": int(t.num.nunique()), "knobs": int(t.knob.nunique()),
     "measured_clean": int(len(clean)), "knobs_clean": int(clean.knob.nunique()),
     "flags": m["flags"].fillna("clean").value_counts().to_dict(),
     "correction": {"b_lnpx": round(bpx, 3), "b_lntan": round(bti, 3), "px0": round(px0, 2), "inc0": round(i0, 1)},
     "knob_clean": clean[["knob", "ktype", "num", "px_m", "inc_local", "emission", "shadow_px", "depth_m", "depth_corr_m"]].round(2).to_dict("records"),
     "knob_median_corr_m": round(float(clean.depth_corr_m.median()), 1) if len(clean) else None,
     "base_auto_median_corr_m": round(float(base.depth_corr_m.median()), 1), "base_n": int(len(base)),
     "base_auto_corr_m_p10_p90": [round(float(base.depth_corr_m.quantile(q)), 1) for q in (0.1, 0.9)],
     "auto_over_manual_validation": 1.236}
json.dump(R, open(f"{HERE}/knob_measure_record.json", "w"), indent=1); print(json.dumps({k_: v for k_, v in R.items() if k_ != "knob_clean"}, indent=1))
print(clean[["knob", "ktype", "num", "px_m", "inc_local", "emission", "shadow_px", "depth_m", "depth_corr_m"]].round(1).to_string(index=False))
print(m[~m.index.isin(clean.index)][["knob", "num", "px_m", "inc_local", "shadow_px", "depth_m", "flags"]].round(1).to_string(index=False))
# contact sheet: close tile (method QC) and context tile (is the box on the knob?)
tiles = []
for _, r in m.sort_values(["knob", "px_m"]).iterrows():
    img, info = C.read_edr(C.edr_path(r.num)); box = (r.x0, r.x1, r.y0, r.y1)
    fl = r["flags"] if isinstance(r["flags"], str) else "clean"
    cap = f"knob {r.knob} EN{r.num}M\n{r.px_m:.0f} m/px i={r.inc_local:.0f} d={r.depth_m:.0f} m\n{fl}"
    tiles.append(Q.tile(img, info, box, r.to_dict(), cap, size=300))
    H, W = img.shape; half = max(60, (r.kaxis if r.kaxis > 0 else 5.0) * 1000 / r.px_m * 0.9); cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
    X0, X1, Y0, Y1 = int(max(0, cx - half)), int(min(W, cx + half)), int(max(0, cy - half)), int(min(H, cy + half))
    sub = img[Y0:Y1, X0:X1]; lo, hi = np.percentile(sub, [1, 99.5]); s = np.clip((sub - lo) / (hi - lo + 1e-6), 0, 1)
    z = 260 / max(sub.shape); im = Image.fromarray((s * 255).astype(np.uint8)).resize((max(1, int(sub.shape[1] * z)), max(1, int(sub.shape[0] * z))), Image.BILINEAR).convert("RGB")
    d = ImageDraw.Draw(im); d.rectangle([(r.x0 - X0) * z, (r.y0 - Y0) * z, (r.x1 - X0) * z, (r.y1 - Y0) * z], outline=(255, 140, 0), width=2)
    T = Image.new("RGB", (300, 300)); T.paste(im, (0, 0)); dd = ImageDraw.Draw(T); dd.text((2, 262), f"context {2*half*r.px_m/1000:.1f} km\n{r.ktype}", fill=(0, 255, 0))
    tiles.append(T)
os.makedirs(f"{HERE}/qc", exist_ok=True); print(Q.sheet(tiles, ncol=6, size=300, path=f"{HERE}/qc/knob_sheet.png"))
