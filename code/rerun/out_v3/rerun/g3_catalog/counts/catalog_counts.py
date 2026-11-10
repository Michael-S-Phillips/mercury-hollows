#!/usr/bin/env python3
"""g3_catalog (6 Nov 2026): catalog counts of Section 2.1, S1 and S2 for v2, v3 and v3b. New script; changes no existing file.
Hollows = canonical hollows not removed by the screen (v3: is_spurious_v3 False; v3b: also is_uncertain_v3 False), so the
64 v2 traces > 10 km stay in the count as in the paper's 6,731. Non-trace hollows = ~is_crater_trace (6,667 in v2).
Fields: (a) existing field_id with >= 1 remaining hollow (out_v3 README: ids unchanged, empty fields disappear);
(b) fields rebuilt by 20 km single linkage over the remaining hollows, as in build_hollow_products.py, because removing a
hollow can split a field. Pairs and detections from inputs/completeness_pairs_{v}.csv (S2)."""
import os, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
HE = os.path.expanduser("~/mnt/hollows-emissivity"); G3 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V3 = f"{HE}/code/rerun/out_v3"; RM = 2439.4
canon = {"v2": f"{HE}/code/rerun/out_v2/hollows_canonical.csv", "v3": f"{V3}/hollows_canonical_v3.csv", "v3b": f"{V3}/hollows_canonical_v3b.csv"}
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
def relink(d):
    T = cKDTree(xyz(d.lat.values, d.lon.values)); pr = np.array(list(T.query_pairs(2 * np.sin(20.0 / (2 * RM)))))
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(d), len(d))); return connected_components(g, directed=False)[0]
R = {}
for v, f in canon.items():
    c = pd.read_csv(f)
    rem = ~c.is_spurious_v3.astype(bool) if v != "v2" else np.ones(len(c), bool)
    if v == "v3b": rem &= ~c.is_uncertain_v3.astype(bool)
    h = c[rem]; nt = c[~c.is_crater_trace.astype(bool)]
    pairs = pd.read_csv(f"{G3}/inputs/completeness_pairs_{v}.csv")
    R[v] = dict(hollows=int(len(h)), nontrace_hollows=int(len(nt)), traces_gt10km_v2=int(c.is_crater_trace_v2.sum()) if v != "v2" else int(c.is_crater_trace.sum()),
                removed=int((~rem).sum()) if v != "v2" else 0,
                fields_existing_ids=int(h.field_id.nunique()), fields_relinked_20km=int(relink(h)),
                fields_with_nontrace=int(nt.field_id.nunique()),
                images_with_hollow_best_detection=int(h.image.nunique()),
                detections_of_nontrace_hollows=int(nt.n_det.sum()),
                pairs=int(len(pairs)), pair_detections=int(pairs.detected.sum()), pair_hollows=int(pairs.h.nunique()),
                nontrace_gt_0p3km=int((nt.D_eq_km > 0.3).sum()), nontrace_gt_1km=int((nt.D_eq_km > 1).sum()),
                nontrace_gt_3km=int((nt.D_eq_km > 3).sum()), nontrace_gt_5km=int((nt.D_eq_km > 5).sum()),
                median_D_nontrace_km=round(float(nt.D_eq_km.median()), 3))
json.dump(R, open(f"{G3}/counts/catalog_counts.json", "w"), indent=1); print(json.dumps(R, indent=1))
