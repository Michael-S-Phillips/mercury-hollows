"""Contact sheet of validation pairs: box (orange), measured shadow (cyan), registered Blewett points (magenta dots),
caption: image fid / auto depth vs Blewett depths (m) / flags. argv: out.png seed n [sort: random|ratio_hi|ratio_lo]"""
import sys, numpy as np, pandas as pd, dep_common as C, qc_render as Q
from PIL import ImageDraw
out, seed, n, mode = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
v = pd.read_csv("validation.csv", dtype={"num": str})
v = v[v.emi_le30 & ~v["flags"].fillna("").str.contains("no_shadow")] if mode != "no_shadow" else v[v.emi_le30 & v["flags"].fillna("").str.contains("no_shadow")]
if mode != "no_shadow": v = v[v.clean]
pts = pd.read_csv("registration_points.csv", dtype={"num": str}); b = C.load_blewett().set_index("bid")
if mode == "random": s = v.sample(n, random_state=seed)
elif mode == "ratio_hi": s = v.sort_values("ratio_max", ascending=False).head(n)
elif mode == "ratio_lo": s = v.sort_values("ratio_max").head(n)
elif mode == "poleward": s = v[v.lat >= 60].sort_values("num").head(n)
else: s = v.sample(min(n, len(v)), random_state=seed)
tiles = []
for _, r in s.iterrows():
    img, info = C.read_edr(C.edr_path(r.num)); box = (r.x0, r.x1, r.y0, r.y1)
    meas = dict(run_u0=r.run_u0, run_u1=r.run_u1, run_w=r.run_w)
    bd = "/".join(f"{x:.0f}" for x in sorted(b.loc[list(map(int, str(r.blewett_ids).split(";"))), "depth"]))
    t = Q.tile(img, info, box, meas, f"{r.num} f{r.fid} px{r.px_m:.0f} i{r.inc_local:.0f}\nauto {r.depth_m:.0f} ({r.shadow_px:.1f}px) r{r.ratio_max:.2f} | B16 {bd}\n{r['flags']}")
    # Blewett points
    H, W = img.shape; bw = max(box[1] - box[0], box[3] - box[2]); m = max(bw * 0.6, 15)
    X0 = int(max(0, box[0] - m)); Y0 = int(max(0, box[2] - m)); X1 = int(min(W, box[1] + m)); Y1 = int(min(H, box[3] + m))
    z = (260 - 40) / max(Y1 - Y0, X1 - X0); d = ImageDraw.Draw(t)
    for _, p in pts[pts.fid == r.fid].iterrows():
        d.ellipse([(p.u - X0) * z - 3, (p.v - Y0) * z - 3, (p.u - X0) * z + 3, (p.v - Y0) * z + 3], outline=(255, 0, 255), width=2)
    tiles.append(t)
print(Q.sheet(tiles, ncol=6, path=out))
