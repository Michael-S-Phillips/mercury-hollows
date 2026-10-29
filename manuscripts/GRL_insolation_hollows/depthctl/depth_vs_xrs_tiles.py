#!/usr/bin/env python3
"""Adds 'tile15_clustered' to depthctl/depth_vs_xrs_record.json: the depth-composition regressions of depth_vs_xrs.py with
standard errors clustered by 15 x 15 degree tile (instead of by image). Recreated for the archive (29 Oct 2026); the block
was first computed inline and the code was not saved. This script reproduces the archived values exactly.
Run after depth_vs_xrs.py, which rewrites the record without this block."""
import os, json, numpy as np, statsmodels.formula.api as smf
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, "depth_vs_xrs.py")
exec(open(SRC).read().replace("__file__", repr(SRC)).split('R = {"n": int(len(d))}')[0])   # loads d (with sampled ratios) and els
REC = os.path.join(HERE, "depth_vs_xrs_record.json"); R = json.load(open(REC))
d["tile"] = (np.floor(d.latitude / 15)).astype(int).astype(str) + "_" + (np.floor((d.longitude % 360) / 15)).astype(int).astype(str)
T = {}
for e in els:
    s = d[d[e].notna()].copy(); s["z"] = (s[e] - s[e].mean()) / s[e].std()
    m = smf.ols("ld ~ z + lp + ltan", data=s).fit(cov_type="cluster", cov_kwds={"groups": s.tile.astype("category").cat.codes})
    lo, hi, sd = s[e].quantile(.1), s[e].quantile(.9), s[e].std(); b, se = m.params["z"], m.bse["z"]
    T[e] = dict(n_tiles=int(s.tile.nunique()), coef=float(b), se=float(se), p=float(m.pvalues["z"]),
                ratio=[float(10 ** (b / sd * (hi - lo))), float(10 ** ((b - 1.96 * se) / sd * (hi - lo))), float(10 ** ((b + 1.96 * se) / sd * (hi - lo)))])
R["tile15_clustered"] = T
json.dump(R, open(REC, "w"), indent=1)
print(json.dumps(T, indent=1))
