"""Sensitivity of the validation statistics to the measurement settings (run on the matched boxes of validation.csv's
registration). Variants: v1 (no edge exclusion, tail flag only), v2 (edge exclusion + lit-tail requirement; the adopted
method), v2_p90 (90th percentile instead of the maximum of the smoothed per-line lengths), v2_frac0.4 / v2_frac0.6
(threshold 40% / 60% of the way from the box's darkest level to the surroundings), v2_nofs (no foreshortening correction).
Writes validation_sensitivity.json."""
import json, numpy as np, pandas as pd, dep_common as C, shadow_measure as M, validate as V
pts = pd.read_csv("registration_points.csv", dtype={"num": str}); b = C.load_blewett().set_index("bid")
mp = pts[(pts.fid >= 0) & pts.frame_accepted].copy(); mp["bd"] = b.loc[mp.bid, "depth"].values
g = mp.groupby("fid").agg(nb=("bid", "size"), bmax=("bd", "max")).reset_index()
k = C.load_bickel(with_corners=False).set_index("fid").loc[g.fid]
variants = {"v1": dict(EDGE_FRAC=0.0, REQUIRE_TAIL=False, STAT="max", frac=0.5, fs=True),
            "v2": dict(EDGE_FRAC=0.10, REQUIRE_TAIL=True, STAT="max", frac=0.5, fs=True),
            "v2_p90": dict(EDGE_FRAC=0.10, REQUIRE_TAIL=True, STAT="p90", frac=0.5, fs=True),
            "v2_frac0.4": dict(EDGE_FRAC=0.10, REQUIRE_TAIL=True, STAT="max", frac=0.4, fs=True),
            "v2_frac0.6": dict(EDGE_FRAC=0.10, REQUIRE_TAIL=True, STAT="max", frac=0.6, fs=True),
            "v2_nofs": dict(EDGE_FRAC=0.10, REQUIRE_TAIL=True, STAT="max", frac=0.5, fs=False)}
cache = {}
R = {}
for name, p in variants.items():
    M.EDGE_FRAC, M.REQUIRE_TAIL, M.STAT = p["EDGE_FRAC"], p["REQUIRE_TAIL"], p["STAT"]
    rows = []
    for fid, r in k.iterrows():
        if r.num not in cache: cache[r.num] = C.read_edr(C.edr_path(r.num))
        img, info = cache[r.num]; N, L = info["samples"], info["lines"]
        o = M.measure(img, info, (r.x0 * N, r.x1 * N, r.y0 * L, r.y1 * L), frac=p["frac"])
        inc = C.local_incidence(r.lat, r.lon180, info["subsol_lat"], info["subsol_lon"])
        fs = M.foreshortening(info) if p["fs"] else 1.0
        rows.append(dict(fid=fid, image=r.num, emission=info["emi"], flags=o["flags"],
                         depth_m=o["L_px"] * info["px"] * fs / np.tan(np.radians(inc))))
    d = pd.DataFrame(rows).merge(g, on="fid")
    clean = ~d["flags"].fillna("").str.contains("|".join(V.BAD)) & (d.depth_m > 0) & (d.emission <= 30)
    d = d[clean]
    R[name] = dict(params=p, all=V.summarize(d, "depth_m", "bmax"), single=V.summarize(d[d.nb == 1], "depth_m", "bmax"),
                   n_boxes_emi_le30=int((pd.DataFrame(rows).emission <= 30).sum()))
    s = R[name]["all"]; print(f"{name:12s} n={s['n']} r={s['pearson_r_log']:.3f} ratio={s['median_ratio']:.3f} {np.round(s['median_ratio_ci95'],3)} mad={s['mad_log_ratio']:.3f} | single n={R[name]['single']['n']} ratio={R[name]['single']['median_ratio']:.3f} r={R[name]['single']['pearson_r_log']:.3f}")
json.dump(R, open("validation_sensitivity.json", "w"), indent=1)
