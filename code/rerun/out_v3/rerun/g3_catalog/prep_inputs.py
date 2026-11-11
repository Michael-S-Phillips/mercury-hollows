#!/usr/bin/env python3
"""g3_catalog (6 Nov 2026): derived inputs for the v3/v3b reruns. New script; it changes no existing file.
1. inputs/completeness_pairs_{v}.csv: code/rerun/out_v2/completeness_pairs.csv restricted to hollows that are non-trace in
   hollows_canonical_{v}.csv. In the pairs file, h is the row index of the v2 non-trace hollows (build_hollow_products.py:
   Hn = H.reset_index(drop=True)), not hollow_id; the map is checked against D_km. A column hollow_id is appended.
2. inputs/crater_dag_nodes_{v}.csv: analysis_v9/crater_dag_nodes.csv with n_hollows and hosts recomputed from
   host_assignments_{v}.csv (hosts = n_hollows > 0, which is how the v9 table defines it; checked on v2 below).
   Crater geometry, parents, completeness and controls are unchanged.
Also v2 copies (identity filter) so each copied script can be checked against its v2 record."""
import os, json, numpy as np, pandas as pd
HE = os.path.expanduser("~/mnt/hollows-emissivity"); OUT = os.path.dirname(os.path.abspath(__file__)) + "/inputs"
V3 = f"{HE}/code/rerun/out_v3"
canon = {"v2": f"{HE}/code/rerun/out_v2/hollows_canonical.csv", "v3": f"{V3}/hollows_canonical_v3.csv", "v3b": f"{V3}/hollows_canonical_v3b.csv"}
hosts = {"v2": f"{HE}/manuscripts/insolation_thermal/data/host_assignments.csv", "v3": f"{V3}/host_assignments_v3.csv", "v3b": f"{V3}/host_assignments_v3b.csv"}
v2 = pd.read_csv(canon["v2"]); hn = v2[~v2.is_crater_trace.astype(bool)].reset_index(drop=True)
CV = pd.read_csv(f"{HE}/code/rerun/out_v2/completeness_pairs.csv")
CV["hollow_id"] = hn.hollow_id.values[CV.h.values]
assert np.allclose(CV.D_km.values, hn.D_eq_km.values[CV.h.values]), "h -> hollow_id map failed"
nodes = pd.read_csv(f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_v9/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str})
cr = nodes.kind == "crater"
assert ((nodes.loc[cr, "hosts"] == 1) == (nodes.loc[cr, "n_hollows"] > 0)).all(), "hosts != (n_hollows > 0) in v9 nodes"
rec = {}
for v in ("v2", "v3", "v3b"):
    c = pd.read_csv(canon[v]); keep = set(c.hollow_id[~c.is_crater_trace.astype(bool)])
    cv = CV[CV.hollow_id.isin(keep)]; cv.to_csv(f"{OUT}/completeness_pairs_{v}.csv", index=False)
    h = pd.read_csv(hosts[v], usecols=["host_uid"]); k = h.host_uid.dropna().astype(int).astype(str).value_counts()
    n = nodes.copy(); n["n_hollows"] = n.id.map(k).fillna(0).astype(int)
    n.loc[n.kind == "crater", "hosts"] = (n.loc[n.kind == "crater", "n_hollows"] > 0).astype(int)
    if v == "v2":
        assert (n.n_hollows.values == nodes.n_hollows.values).all()
    n.to_csv(f"{OUT}/crater_dag_nodes_{v}.csv", index=False)
    rec[v] = dict(pairs=len(cv), pair_hollows=int(cv.h.nunique()), detections=int(cv.detected.sum()),
                  crater_hosts=int(n.loc[n.kind == "crater", "hosts"].sum()))
json.dump(rec, open(f"{OUT}/prep_inputs_record.json", "w"), indent=1); print(json.dumps(rec, indent=1))
