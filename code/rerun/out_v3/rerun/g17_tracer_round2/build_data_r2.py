"""g17 (7 Nov 2026): data.json for round 2 of the tracer. Keeps the 64 round-1 craters as they are (group A) and appends groups B
(50 more hosts), C (Class 1-2 craters with hollows, peaks only) and D (PPS sample of Class 1-2 craters without hollows, peaks only)
from g17_tracer_round2/assets/meta, with curvature-v3 outlines simplified as in g12 build_tracer_data.py. New craters load their
images from pk/<uid>.json (pack = true)."""
import json, os, shutil, numpy as np
from shapely.geometry import Polygon
from shapely.validation import make_valid
T = "/tmp/claude-0/-home-claude/40c2499e-3a3d-5b4e-8387-2d699ac58079/scratchpad/tracer_r2"; A = f"{T}/bundle/assets"; O = f"{T}/site"
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
def entry(m):
    la0, lo0 = m["lat"], m["lon"]; shapes = []
    for rings in m["v3"]["peak"]:
        xy = [np.column_stack(fwd(np.array(r)[:, 0], np.array(r)[:, 1], la0, lo0)) for r in rings]
        g = make_valid(Polygon(xy[0], xy[1:])); parts = [g] if isinstance(g, Polygon) else [q for q in getattr(g, "geoms", []) if isinstance(q, Polygon)]
        for q in parts:
            q = q.simplify(0.5, preserve_topology=True)
            if q.is_empty or q.area < 2.0: continue
            for qq in ([q] if isinstance(q, Polygon) else list(q.geoms)):
                shapes.append(dict(kind="peak", ll=ll(list(qq.exterior.coords), la0, lo0)))
                for h in qq.interiors:
                    if Polygon(h).area >= 2.0: shapes.append(dict(kind="hole", ll=ll(list(h.coords), la0, lo0)))
    fxy = np.column_stack(fwd(np.array(m["v3"]["foot"])[:, 0], np.array(m["v3"]["foot"])[:, 1], la0, lo0))
    shapes.append(dict(kind="foot", ll=ll(list(Polygon(fxy).simplify(0.15, preserve_topology=True).exterior.coords), la0, lo0)))
    return dict(uid=m["uid"], name=(m["name"] or "").strip(), D=m["D"], lat=la0, lon=lo0, cent=m["cent"], filled=m["filled"], floor_m=round(m["floor_m"], 1), depth_m=round(m["depth_m"], 1),
                n=m["n"], H=round(m["H_km"], 4), img=m["grids"]["img"]["n"], dem=dict(n=m["grids"]["dem"]["n"], zmin=m["grids"]["dem"]["zmin"]),
                ring_sectors=m["v3"]["n_ring_sectors"], v3=shapes, hollows=m["hollows"], masked=m["masked"])
live = json.load(open(f"{T}/live/data.json")); sel = json.load(open(f"{T}/bundle/selection.json"))
for c in live["craters"]: c["group"] = "A"
seen = {c["uid"] for c in live["craters"]}; out = list(live["craters"])
for g, key in (("B", "group_B"), ("C", "group_C"), ("D", "group_D")):
    for k, s in enumerate(sel[key]):
        assert s["uid"] not in seen, s["uid"]; seen.add(s["uid"])
        e = entry(json.load(open(f"{A}/meta/{s['uid']}.json"))); e.update(group=g, rank=k + 1, pack=True)
        if g in ("C", "D"): e["kclass"] = s["kclass"]
        if g == "D": e["pi"] = s["pi"]
        out.append(e); shutil.copy(f"{T}/bundle/pk/{s['uid']}.json", f"{O}/pk/{s['uid']}.json")
GROUPS = {"A": dict(title="Hosts, round 1", task="full", start="v3"), "B": dict(title="Hosts, round 2", task="full", start="v3"),
          "C": dict(title="Degraded craters with hollows", task="peaks", start="empty"), "D": dict(title="Degraded craters, random sample", task="peaks", start="empty")}
json.dump(dict(version="g17-r2 (7 Nov 2026): round 1 hosts plus 50 hosts, 84 Class 1-2 craters with hollows and 60 sampled Class 1-2 craters",
               n_all=live["n_all"], n_hosts=live["n_hosts"], groups=GROUPS, craters=out), open(f"{O}/data.json", "w"), separators=(",", ":"))
from collections import Counter
print(Counter(c["group"] for c in out), "data.json KB", os.path.getsize(f"{O}/data.json") // 1024, "packs", len(os.listdir(f"{O}/pk")))
