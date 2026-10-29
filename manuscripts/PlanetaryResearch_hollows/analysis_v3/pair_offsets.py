import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Relative position error between overlapping NAC images: for each pair of images that both detect the
same hollows, match detections as mutual nearest neighbours within 3 km and take the median east/north
offset of the matches (a systematic shift between the two images). Done for recomputed (ISIS camera model)
and released positions. Pairs need >= 3 matches. Pixel scale class uses the coarser image of the pair."""
import os
import numpy as np, pandas as pd, json, itertools
from scipy.spatial import cKDTree
d = pd.read_csv(os.path.expanduser(_REPO + "")+"/code/rerun/out_v2/detections_qc.csv")
R = 2439.4
def en(lat, lon, lat0, lon0):
    return np.c_[np.radians(((lon - lon0 + 180) % 360) - 180) * np.cos(np.radians(lat0)) * R, np.radians(lat - lat0) * R]
px = d.groupby("image").px_m.first()
res = {}
for tag, (la, lo) in {"new": ("lat", "lon"), "rel": ("lat_release", "lon_release")}.items():
    rows = []
    # candidate image pairs: share a canonical hollow id (recomputed positions) -> overlapping images
    pairs = set()
    for hid, g in d.groupby("hollow_id"):
        ims = sorted(g.image.unique())
        for a, b in itertools.combinations(ims, 2): pairs.add((a, b))
    G = {k: g for k, g in d.groupby("image")}
    for a, b in pairs:
        A, B = G[a], G[b]
        lat0, lon0 = A[la].mean(), A[lo].mean()
        pa, pb = en(A[la].values, A[lo].values, lat0, lon0), en(B[la].values, B[lo].values, lat0, lon0)
        # for released positions the images may be offset by km; allow a coarse global shift first
        ta, tb = cKDTree(pa), cKDTree(pb)
        da, ia = tb.query(pa); db, ib = ta.query(pb)
        mutual = [(i, ia[i]) for i in range(len(pa)) if ib[ia[i]] == i and da[i] < 3.0]
        if len(mutual) < 3: continue
        off = np.array([pb[j] - pa[i] for i, j in mutual])
        med = np.median(off, axis=0)
        resid = np.linalg.norm(off - med, axis=1)
        rows.append(dict(a=a, b=b, n=len(mutual), shift_km=float(np.linalg.norm(med)), resid_med_km=float(np.median(resid)),
                         px_max=float(max(px[a], px[b]))))
    p = pd.DataFrame(rows); p.to_csv(f"pair_offsets_{tag}.csv", index=False)
    out = {"n_pairs": int(len(p))}
    for lo_, hi_ in [(0, 20), (20, 40), (40, 60), (0, 60)]:
        s = p[(p.px_max >= lo_) & (p.px_max < hi_)]
        out[f"{lo_}-{hi_}"] = dict(n_pairs=int(len(s)), shift_med_km=float(s.shift_km.median()), shift_med_px=float((s.shift_km*1000/s.px_max).median()),
                                   shift_p90_km=float(s.shift_km.quantile(.9)), resid_med_km=float(s.resid_med_km.median()),
                                   resid_med_px=float((s.resid_med_km*1000/s.px_max).median()))
    res[tag] = out
print(json.dumps(res, indent=1)); json.dump(res, open("pair_offsets_record.json", "w"), indent=1)
