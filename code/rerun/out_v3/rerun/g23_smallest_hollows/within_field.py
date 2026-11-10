#!/usr/bin/env python3
"""g23 step 8 (9 Nov): is depth-size a common-cause effect? Michael: smaller hollows may form in less volatile-rich material (lower phi),
which would make them shallower without size itself mattering. Hollows in the same field, measured in the same frame, share local
material as far as fields resolve it, and share pixel scale and incidence. The within-group slope of ln depth on ln D (group fixed
effects) is free of differences in phi between fields; the between-group slope (group means) carries them.
Groups: field_id x frame, and host_uid x frame. Data: hand depths (bench_depths.csv) and automatic depths (depth_vs_size_classed.csv).
Writes within_field.json."""
import json, numpy as np, pandas as pd
RR = "../g2_insolation/inputs/v3"
H = pd.read_csv(f"{RR}/hollows_canonical.csv")[["hollow_id", "field_id"]]; ha = pd.read_csv(f"{RR}/host_assignments.csv")[["hollow_id", "host_uid"]]
man = json.load(open("manifest_meas.json"))
def prep(d, frame):
    d = d.merge(H, on="hollow_id", how="left").merge(ha, on="hollow_id", how="left"); d["frame"] = frame(d)
    d["lh"] = np.log(d.depth_m); d["lD"] = np.log(d.D); d["lpx"] = np.log(d.px); d["lti"] = np.log(np.tan(np.radians(d.inc))); return d
b = pd.read_csv("bench_depths.csv"); b = b[b.cls.isin(["hollow", "unsure"]) | b.cls.isna()]; b = b[(b.shadow == "measured") & (b.depth_m > 0)].copy(); b["D"] = b.D_m
b = prep(b, lambda d: [man[str(int(i))]["img"] for i in d.hollow_id])
a = pd.read_csv("depth_vs_size_classed.csv"); a = a[a.cls == "clean"].copy(); a["D"] = a.D_can_m; a = prep(a, lambda d: d.num.astype(str))
def within(d, gcol):
    g = d.dropna(subset=[gcol]).copy(); g["grp"] = g[gcol].astype(str) + "|" + g.frame.astype(str)
    cnt = g.groupby("grp").lD.transform("size"); rng = g.groupby("grp").lD.transform(lambda x: x.max() - x.min()); g = g[(cnt >= 2) & (rng > 0)]
    if len(g) < 6: return dict(n=len(g))
    ng = g.grp.nunique(); dm = lambda c: g[c] - g.groupby("grp")[c].transform("mean")
    x, y = dm("lD"), dm("lh"); b1 = float((x * y).sum() / (x * x).sum()); r = y - b1 * x
    # cluster-robust SE by group
    s = g.assign(xr=x * r).groupby("grp").xr.sum(); se = float(np.sqrt((s ** 2).sum()) / (x * x).sum())
    # between: group means with px and tan i
    m = g.groupby("grp")[["lh", "lD", "lpx", "lti"]].mean(); X = np.c_[np.ones(len(m)), m.lD, m.lpx, m.lti]
    bb, *_ = np.linalg.lstsq(X, m.lh, rcond=None); rr = m.lh - X @ bb; cov = rr.var(ddof=4) * np.linalg.inv(X.T @ X) if len(m) > 5 else np.full((4, 4), np.nan)
    return dict(n=len(g), groups=ng, D_range_within_median=round(float(np.exp(g.groupby("grp").lD.agg(lambda x: x.max() - x.min()).median())), 2),
                within_slope=round(b1, 3), within_se=round(se, 3), between_slope=round(float(bb[1]), 3), between_se=round(float(np.sqrt(cov[1, 1])), 3))
out = {}
for nm, d in [("hand", b), ("auto", a)]:
    for gcol in ["field_id", "host_uid"]:
        out[f"{nm}_{gcol}"] = within(d, gcol); print(nm, gcol, out[f"{nm}_{gcol}"])
    X = np.c_[np.ones(len(d)), d.lD, d.lpx, d.lti]; bb, *_ = np.linalg.lstsq(X, d.lh, rcond=None); out[f"{nm}_pooled_slope"] = round(float(bb[1]), 3); print(nm, "pooled", out[f"{nm}_pooled_slope"], "n", len(d))
json.dump(out, open("within_field.json", "w"), indent=1)
