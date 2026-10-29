#!/usr/bin/env python3
"""Adds 'eq_face_slope_ge2' and 'eq_face_slope_ge4' to slope_test/depth_slope_record.json: the equator-facing term of
ln depth ~ slope + eq_face + 1/T + ln px + ln tan(i), restricted to measurements on slopes >= 2 and >= 4 degrees, with
image-clustered standard errors. Values are [coef, se, p, n measurements, n images]. Recreated for the archive (29 Oct 2026);
the block was first computed inline and the code was not saved. This script reproduces the archived values exactly.
Run from slope_test/ after depth_slope_geometry.py, which rewrites the record without this block."""
import os, json, statsmodels.formula.api as smf
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, "depth_slope_geometry.py"); os.chdir(HERE)
ns = {"__name__": "depth_slope_geometry_head"}
exec(open(SRC).read().replace("__file__", repr(SRC)).split('json.dump(R, open("depth_slope_record.json"')[0], ns)
d = ns["d"]; REC = os.path.join(HERE, "depth_slope_record.json"); R = json.load(open(REC))
for th in (2, 4):
    sub = d[d.slope_s1 >= th]
    m = smf.ols("lnh_s1 ~ slope_s1 + eq_face + x + lnpx + ltan", data=sub).fit(cov_type="cluster", cov_kwds={"groups": sub.filename.astype("category").cat.codes})
    R[f"eq_face_slope_ge{th}"] = [float(m.params["eq_face"]), float(m.bse["eq_face"]), float(m.pvalues["eq_face"]), int(len(sub)), int(sub.filename.nunique())]
    print(th, R[f"eq_face_slope_ge{th}"])
json.dump(R, open(REC, "w"), indent=1)
