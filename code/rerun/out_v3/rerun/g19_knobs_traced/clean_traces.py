# g19 copy of g18_traced_r3/clean_traces.py (7 Nov 2026). Changed: keeps all four tracer groups (peaks of groups C and D are needed for the knob comparison) and records each crater's group.
# g18 copy of g13_traced/clean_traces.py (7 Nov 2026). Changed: reads one snapshot JSON (docs keyed by uid) in place of a folder; keeps tracer groups A and B (hosts) only; a round-2 host counts as done only when its wall-foot check is done. Usage: clean_traces.py SNAPSHOT_JSON DATA_JSON OUT
#!/usr/bin/env python3
"""g13 (7 Nov 2026), run where the tracer database snapshot is: cleans the traced outlines into per-crater geometry, in the
tracer's local azimuthal-equidistant frame and back to lon/lat rings. peak = peak outlines minus floor-inside-ring outlines;
foot = wall-foot outline(s); slump = slump outlines (union). Self-intersecting outlines are repaired with shapely make_valid, and holes
smaller than 5% of the outline's area are filled (Scarlatti's ring drawn as one loop becomes an annulus). Writes JSON with status,
pass2 (wall-foot check), note, and geometry as [exterior, holes...] rings of [lon, lat]. Usage: clean_traces.py SNAPSHOT_DIR DATA_JSON OUT"""
import json, glob, sys, numpy as np
from shapely.geometry import Polygon
from shapely.validation import make_valid
from shapely.ops import unary_union
snap, data, out = sys.argv[1:4]; RM = 2439.4
CR = {c["uid"]: c for c in json.load(open(data))["craters"]}
RAW = json.load(open(snap))["docs"]   # g18: one JSON with every database document (ArtifactData read), not a folder
def fwd(C, lon, lat):
    p, p0, dl = np.radians(lat), np.radians(C["lat"]), np.radians(np.asarray(lon) - C["lon"])
    cc = np.clip(np.sin(p0)*np.sin(p) + np.cos(p0)*np.cos(p)*np.cos(dl), -1, 1); c = np.arccos(cc); k = np.where(c > 1e-12, c/np.sin(np.maximum(c, 1e-12)), 1)
    return RM*k*np.cos(p)*np.sin(dl), RM*k*(np.cos(p0)*np.sin(p) - np.sin(p0)*np.cos(p)*np.cos(dl))
def inv(C, x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); rho = np.hypot(x, y); c = rho/RM; p0 = np.radians(C["lat"])
    lat = np.arcsin(np.cos(c)*np.sin(p0) + np.where(rho > 0, y*np.sin(c)*np.cos(p0)/np.maximum(rho, 1e-12), 0))
    lon = C["lon"] + np.degrees(np.arctan2(x*np.sin(c), rho*np.cos(p0)*np.cos(c) - y*np.sin(p0)*np.sin(c)))
    return lon, np.degrees(lat)
polys = lambda g: [g] if g.geom_type == "Polygon" else [q for q in getattr(g, "geoms", []) if q.geom_type == "Polygon"]
res, log = {}, []
for _k in sorted(RAW):
    x = RAW[_k]; x = x.get("data", x); C = CR[x["uid"]]; G = {"peak": [], "hole": [], "foot": [], "slump": [], "rim": []}
    pass   # g19: every group is kept; only the peak outlines are used here
    for s in x["shapes"]:
        X, Y = fwd(C, s["lon"], s["lat"]); P = Polygon(np.column_stack([X, Y]))
        if not P.is_valid:
            fixed = [Polygon(q.exterior, [h for h in q.interiors if Polygon(h).area >= 0.05*Polygon(q.exterior).area]) for q in polys(make_valid(P))]
            log.append(f"{C['name'] or C['uid']}: self-intersecting {s['kind']} outline repaired ({len(fixed)} part(s))"); P = unary_union(fixed)
        G.setdefault(s["kind"], []).append(P)
    peak = unary_union(G["peak"]).difference(unary_union(G["hole"])) if G["peak"] else Polygon()
    U = lambda k: unary_union(G[k]) if G[k] else Polygon()
    def ll(g):
        outl = []
        for q in polys(g):
            if q.is_empty: continue
            rings = [np.array(q.exterior.coords)] + [np.array(h.coords) for h in q.interiors]
            outl.append([[[round(float(a), 5), round(float(b), 5)] for a, b in zip(*inv(C, r[:, 0], r[:, 1]))] for r in rings])
        return outl
    res[x["uid"]] = dict(uid=x["uid"], name=C["name"], group=C.get("group", "A"), status=x["status"], pass2=x.get("pass2", "todo"), note=x.get("note", ""), lat0=C["lat"], lon0=C["lon"],
                         peak_km2=round(peak.area, 1), slump_km2=round(U("slump").area, 1), peak=ll(peak), foot=ll(U("foot")), slump=ll(U("slump")), rim=ll(U("rim")), updated=x["updated"])
json.dump(dict(source="Hollow Host Terrain Tracer (artifact Hnvf7kvMZBqPpkaKvRX6So), collection hosts", geometry="[exterior, holes...] rings of [lon, lat]; lon continuous about the crater center",
               cleaning=log, craters=res), open(out, "w"), separators=(",", ":"))
print("\n".join(log)); print(len(res), "craters;", sum(v["slump_km2"] > 0 for v in res.values()), "with slump")
