"""Assemble data.json for the host tracer from the g12 assets: simplify curvature-v3 peak/ring polygons
(local azimuthal-equidistant km, shapely simplify 0.5 km, drop parts and holes < 2 km2) and copy images."""
import json, os, shutil, numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid
S = "/tmp/claude-0/-home-claude/40c2499e-3a3d-5b4e-8387-2d699ac58079/scratchpad"; A = f"{S}/htracer_src/assets"; O = f"{S}/htracer"
RM = 2439.4
def fwd(lon, lat, la0, lo0):
    p, p0, dl = np.radians(lat), np.radians(la0), np.radians(np.asarray(lon) - lo0)
    cc = np.clip(np.sin(p0)*np.sin(p) + np.cos(p0)*np.cos(p)*np.cos(dl), -1, 1); c = np.arccos(cc); k = np.where(c > 1e-12, c / np.sin(np.maximum(c, 1e-12)), 1.0)
    return RM*k*np.cos(p)*np.sin(dl), RM*k*(np.cos(p0)*np.sin(p) - np.sin(p0)*np.cos(p)*np.cos(dl))
def inv(x, y, la0, lo0):
    x, y = np.asarray(x, float), np.asarray(y, float); rho = np.hypot(x, y); c = rho / RM; p0 = np.radians(la0)
    with np.errstate(invalid="ignore", divide="ignore"):
        lat = np.arcsin(np.cos(c)*np.sin(p0) + np.where(rho > 0, y*np.sin(c)*np.cos(p0)/rho, 0))
        lon = lo0 + np.degrees(np.arctan2(x*np.sin(c), rho*np.cos(p0)*np.cos(c) - y*np.sin(p0)*np.sin(c)))
    return lon, np.degrees(lat)
def ll(ring, la0, lo0):
    xy = np.asarray(ring)[:-1] if np.allclose(ring[0], ring[-1]) else np.asarray(ring)
    lo, la = inv(xy[:, 0], xy[:, 1], la0, lo0); return [[round(float(a), 5), round(float(b), 5)] for a, b in zip(lo, la)]
order = json.load(open(f"{A}/order.json")); out = []
for d in ("img", "dem"): os.makedirs(f"{O}/{d}", exist_ok=True)
for rank, uid in enumerate(order["order"]):
    m = json.load(open(f"{A}/meta/{uid}.json")); la0, lo0 = m["lat"], m["lon"]; shapes = []
    for rings in m["v3"]["peak"]:
        xy = [np.column_stack(fwd(np.array(r)[:, 0], np.array(r)[:, 1], la0, lo0)) for r in rings]
        g = make_valid(Polygon(xy[0], xy[1:]))
        parts = [g] if isinstance(g, Polygon) else [q for q in getattr(g, "geoms", []) if isinstance(q, Polygon)]
        for q in parts:
            q = q.simplify(0.5, preserve_topology=True)
            if q.is_empty or q.area < 2.0: continue
            for qq in ([q] if isinstance(q, Polygon) else list(q.geoms)):
                shapes.append(dict(kind="peak", ll=ll(list(qq.exterior.coords), la0, lo0)))
                for h in qq.interiors:
                    if Polygon(h).area >= 2.0: shapes.append(dict(kind="hole", ll=ll(list(h.coords), la0, lo0)))
    fxy = np.column_stack(fwd(np.array(m["v3"]["foot"])[:, 0], np.array(m["v3"]["foot"])[:, 1], la0, lo0))
    fpoly = Polygon(fxy).simplify(0.15, preserve_topology=True)
    shapes.append(dict(kind="foot", ll=ll(list(fpoly.exterior.coords), la0, lo0)))
    out.append(dict(uid=uid, rank=rank + 1, name=(m["name"] or "").strip(), D=m["D"], lat=la0, lon=lo0, cent=m["cent"], filled=m["filled"], floor_m=m["floor_m"], depth_m=m["depth_m"],
                    n=m["n"], H=round(m["H_km"], 4), img=m["grids"]["img"]["n"], dem=dict(n=m["grids"]["dem"]["n"], zmin=m["grids"]["dem"]["zmin"]),
                    ring_sectors=m["v3"]["n_ring_sectors"], v3=shapes, hollows=m["hollows"], masked=m["masked"]))
    shutil.copy(f"{A}/img/{uid}.jpg", f"{O}/img/{uid}.jpg"); shutil.copy(f"{A}/dem/{uid}.png", f"{O}/dem/{uid}.png")
json.dump(dict(version="g12-v1 (curvature v3, 6 Nov 2026)", n_all=order["n_all"], n_hosts=order["n_hosts"], craters=out), open(f"{O}/data.json", "w"), separators=(",", ":"))
nv = [sum(len(s["ll"]) for s in c["v3"]) for c in out]
print("craters", len(out), "data.json KB", os.path.getsize(f"{O}/data.json") // 1024, "vertices/crater median", int(np.median(nv)), "max", max(nv))
print([(c["name"] or c["uid"], [s["kind"] for s in c["v3"]].count("peak"), [s["kind"] for s in c["v3"]].count("hole")) for c in out[:12]])
