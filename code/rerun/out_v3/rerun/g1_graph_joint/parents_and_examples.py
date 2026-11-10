#!/usr/bin/env python3
"""g1 rerun (6 Nov 2026), new script for numbers that had no saved script:
(1) Supplement S14 'Parents without hollows': parent craters (kind crater) with at least one child in the tests (possible
    rings dropped) that hosts hollows (any assigned hollow), how many host no hollows themselves (any assigned hollow), and
    the completeness C1 at their centers (median; counts >= 0.2 and >= 0.5) against parents that host. The definition was
    chosen because it reproduces the v2 values (149, 75, 0.08, 0.21, 22, 4) exactly.
(2) Figure 3 caption: children (craters >= 5 km in the tests) of Praxiteles and Lermontov and how many host hollows within
    0.9 radii (h09), and the h09 hosting rate of all craters in the tests (1.9% in v2).
Usage: CATALOG=v3 python3 parents_and_examples.py  (reads <CATALOG>/crater_dag_nodes.csv; writes <CATALOG>/parents_and_examples.json)"""
import os, sys, json, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); from hostfiles import HOSTS, CAT, CRATERS
n = pd.read_csv(f"{H}/{CAT}/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str}, low_memory=False)
h = pd.read_csv(HOSTS, usecols=["host_uid", "rho"]).dropna(); h["host_uid"] = h.host_uid.astype(int).astype(str)
n["h09"] = n.id.map(h[h.rho <= 0.9].host_uid.value_counts()).fillna(0).gt(0).astype(int); n["hany"] = n.hosts.fillna(0).astype(int)
R = {}
s = n[(n.kind == "crater") & (~n.possible_ring.astype(bool))]
c = s[(s.parent >= 0) & (s.parent_kind == "crater")]
par = c[c.hany == 1].parent_id.unique(); P = n.set_index("id").loc[par]
no, yes = P[P.hany == 0], P[P.hany == 1]
R["parents_without_hollows"] = dict(parents_with_hosting_child=int(len(P)), host_no_hollows=int(len(no)), median_C1_no=float(no.C1.median()),
                                    median_C1_hosting_parents=float(yes.C1.median()), no_C1_ge_0p2=int((no.C1 >= 0.2).sum()), no_C1_ge_0p5=int((no.C1 >= 0.5).sum()))
t = pd.read_csv(CRATERS, usecols=["unique_id", "name"], dtype={"unique_id": str}); t["name"] = t.name.fillna("").str.strip()
for nm in ["Praxiteles", "Lermontov"]:
    pid = t[t.name == nm].unique_id.iloc[0]; ch = s[s.parent_id == pid]
    R[nm] = dict(id=pid, D=float(n.set_index("id").loc[pid, "D"]), children=int(len(ch)), children_h09=int(ch.h09.sum()), parent_h09=int(n.set_index("id").loc[pid, "h09"]))
R["h09_rate_all_tests"] = dict(hosts=int(s.h09.sum()), craters=int(len(s)), rate=float(s.h09.mean()))
json.dump(R, open(f"{H}/{CAT}/parents_and_examples.json", "w"), indent=1); print(CAT, json.dumps(R))
