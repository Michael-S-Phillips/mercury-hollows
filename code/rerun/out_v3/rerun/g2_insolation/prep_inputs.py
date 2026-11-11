#!/usr/bin/env python3
"""g2_insolation (6 Nov 2026): build the catalog-dependent inputs of the insolation / plains reruns for v2 (reproduction),
v3 and v3b. Writes inputs/<cat>/:
  hollows_canonical.csv      canonical hollows (v2: code/rerun/out_v2; v3/v3b: code/rerun/out_v3/hollows_canonical_<cat>.csv)
  hollows_canonical_inc.csv  rows used by latcorr5.py to map 60-75 deg detections to fields
                             (v2: all rows, as the original; v3/v3b: non-trace rows only, so a removed detection cannot
                             qualify a field)
  fields.csv                 v2: code/rerun/out_v2/fields.csv unchanged (1,228 fields, as used by every v2 field analysis);
                             v3/v3b: the same rows restricted to fields with >= 1 non-trace hollow, n_nontrace recounted.
                             Field ids and centroids are kept from v2 (README of out_v3: fields left with no hollows disappear).
  host_assignments.csv       v2: insolation_thermal/data/host_assignments.csv; v3/v3b: out_v3/host_assignments_<cat>.csv
  test2_craters_h18.csv      insolation_thermal/data/test2_craters_h18.csv with `hosts` reset from host_assignments
                             (= test2_h18.py with its host-assignment path changed)
  crater_features.csv        analysis_workshop/summary_model/crater_features.csv with hosts, h09, n09 recomputed from
                             host_assignments as in build_features.py; parent_hosts and ln_local (functions of hosts, not
                             used by these scripts) are dropped so they cannot be used stale
  canonical_stage_labels.csv code/rerun/latitude/incidence_null/canonical_stage_labels.csv restricted to the non-trace hollows
                             of the catalog (row-aligned with them, as stage_test.py asserts)
"""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
PR = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"; G = f"{ROOT}/manuscripts/insolation_thermal"
SRC = {"v2": dict(canon=f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", ha=f"{G}/data/host_assignments.csv"),
       "v3": dict(canon=f"{ROOT}/code/rerun/out_v3/hollows_canonical_v3.csv", ha=f"{ROOT}/code/rerun/out_v3/host_assignments_v3.csv"),
       "v3b": dict(canon=f"{ROOT}/code/rerun/out_v3/hollows_canonical_v3b.csv", ha=f"{ROOT}/code/rerun/out_v3/host_assignments_v3b.csv")}
F2 = pd.read_csv(f"{ROOT}/code/rerun/out_v2/fields.csv")
T2 = pd.read_csv(f"{G}/data/test2_craters_h18.csv")
CF = pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
L2 = pd.read_csv(f"{ROOT}/code/rerun/latitude/incidence_null/canonical_stage_labels.csv")
H2 = pd.read_csv(SRC["v2"]["canon"]); H2nt = H2[~H2.is_crater_trace.astype(bool)].reset_index(drop=True)
assert (H2nt.thumb.values == L2.thumb.values).all()
summary = {}
for cat, src in SRC.items():
    D = f"{HERE}/inputs/{cat}"; os.makedirs(D, exist_ok=True)
    H = pd.read_csv(src["canon"]); assert (H.hollow_id.values == H2.hollow_id.values).all()
    H.to_csv(f"{D}/hollows_canonical.csv", index=False)
    nt = ~H.is_crater_trace.astype(bool)
    (H if cat == "v2" else H[nt]).to_csv(f"{D}/hollows_canonical_inc.csv", index=False)
    cnt = H[nt].groupby("field_id").size()
    F = F2.copy()
    if cat != "v2":
        F["n_nontrace"] = F.field_id.map(cnt).fillna(0).astype(int); F = F[F.n_nontrace > 0]
    F.to_csv(f"{D}/fields.csv", index=False)
    ha = pd.read_csv(src["ha"]); ha.to_csv(f"{D}/host_assignments.csv", index=False)
    hid = set(ha.host_uid.dropna().astype(int))
    T = T2.copy(); T["hosts"] = T.unique_id.astype(int).isin(hid).astype(int); T.to_csv(f"{D}/test2_craters_h18.csv", index=False)
    h = ha[["host_uid", "rho"]].dropna().copy(); h["host_uid"] = h.host_uid.astype(int).astype(str)
    C = CF.drop(columns=["parent_hosts", "ln_local"]).copy()
    C["hosts"] = C.id.isin({str(x) for x in hid}).astype(int)
    C["n09"] = C.id.map(h[h.rho <= 0.9].host_uid.value_counts()).fillna(0).astype(int); C["h09"] = (C.n09 > 0).astype(int)
    if cat == "v2":   # reproduction check against the stored v2 table
        for k in ("hosts", "h09", "n09"): assert (C[k].values == CF[k].values).all(), k
        assert (T.hosts.values == T2.hosts.values).all()
    C.to_csv(f"{D}/crater_features.csv", index=False)
    keep = nt.values[~H2.is_crater_trace.astype(bool).values]       # v2 non-trace rows that stay non-trace
    L = L2[keep].reset_index(drop=True); assert (L.thumb.values == H[nt].thumb.values).all()
    L.to_csv(f"{D}/canonical_stage_labels.csv", index=False)
    summary[cat] = dict(nontrace_hollows=int(nt.sum()), fields=int(len(F)), hosts_any=int(len(hid)),
                        test2_hosts=int(T.hosts.sum()), cf_hosts=int(C.hosts.sum()), cf_h09=int(C.h09.sum()),
                        stage_labels_1_2=int(L.stage08.isin([1, 2]).sum()))
    print(cat, summary[cat])
json.dump(summary, open(f"{HERE}/inputs/prep_summary.json", "w"), indent=1)
