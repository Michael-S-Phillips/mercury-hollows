# g19 copy of g18_traced_r3/traced_classes.py (7 Nov 2026). Changed: reads traces/traces_clean_r3k.json (all tracer groups; status done).
# g18 copy of g13_traced/traced_classes.py (7 Nov 2026). Changed: folder g18_traced_r3; traces from traces/traces_clean_r3.json (round-2 snapshot, hosts of tracer groups A and B); G13_REV default r3 (slump counts as wall, as in r2); outputs named r3.
# g13 revision (7 Nov 2026): env G13_REV (r1 | r2, default r2): r2 reads traces/traces_clean_r2.json (Michael's refined wall feet and slump outlines, 7 Nov) and counts slump as wall (slump inside the wall foot becomes wall, peak still wins); r1 reproduces the 6 Nov results.
"""g13 (6 Nov 2026), new module: terrain classes from Michael's traced outlines (traces/traces_clean.json, status done) where they
exist, curvature classes elsewhere (cls = v3: g10 curv_classes_v3, default; v2: g10 curv_classes as in g11). Traced rules:
peak = inside peak outlines minus floor-inside-ring outlines; wall = outside the wall foot out to 1.05 R, not peak;
floor = inside the wall foot, not peak. Point-in-polygon in a local azimuthal-equidistant frame about the catalog center,
the frame the tracer used. Same rules as peak_wall_traced.py."""
import os, sys, json, numpy as np
from matplotlib.path import Path
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(HERE), "g10_threshold"))
from curv_classes_v3 import classify_curv as _v3
from curv_classes import classify_curv as _v2
RM = 2439.4
REV = os.environ.get("G13_REV", "r3")
TR = {int(k): v for k, v in json.load(open(os.path.join(HERE, "traces", "traces_clean_r3k.json")))["craters"].items() if v["status"] == "done"}
def aeqd(lon, lat, la0, lo0):
    p, p0, dl = np.radians(lat), np.radians(la0), np.radians(lon - lo0)
    cc = np.clip(np.sin(p0)*np.sin(p) + np.cos(p0)*np.cos(p)*np.cos(dl), -1, 1); c = np.arccos(cc)
    k = np.where(c > 1e-12, c/np.sin(np.maximum(c, 1e-12)), 1.0)
    return RM*k*np.cos(p)*np.sin(dl), RM*k*(np.cos(p0)*np.sin(p) - np.sin(p0)*np.cos(p)*np.cos(dl))
def _inside(polys, X, Y, t):
    pts = np.column_stack([X.ravel(), Y.ravel()]); m = np.zeros(len(pts), bool)
    for rings in polys:
        rings = [r for r in rings if len(r) >= 3]
        if not rings: continue
        P = [np.column_stack(aeqd(np.array(r)[:, 0], np.array(r)[:, 1], t["lat0"], t["lon0"])) for r in rings]
        q = Path(P[0]).contains_points(pts)
        for h in P[1:]: q &= ~Path(h).contains_points(pts)
        m |= q
    return m.reshape(X.shape)
def classes_for(uid, z, rho, az, R, cent, LON, LAT, cls="v3"):
    uid = int(uid)
    if uid in TR:
        t = TR[uid]; X, Y = aeqd(LON, LAT, t["lat0"], t["lon0"]); fin = np.isfinite(z)
        pk = _inside(t["peak"], X, Y, t) & fin; ft = _inside(t["foot"], X, Y, t)
        sl = _inside(t.get("slump", []), X, Y, t) if REV != "r1" else np.zeros_like(ft)
        return dict(peak=pk, wall=((~ft) | sl) & (rho <= 1.05) & ~pk & fin, flo=ft & ~sl & ~pk & fin, source="traced")
    cent = "" if cent in ("nan", "None") else cent
    K = (_v3 if cls == "v3" else _v2)(z, rho, az, R, cent); K["source"] = f"curv_{cls}"; return K
