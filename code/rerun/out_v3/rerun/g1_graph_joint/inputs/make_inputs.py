#!/usr/bin/env python3
"""Derived inputs for the g1 reruns (6 Nov 2026). The crater table GRL_insolation_hollows/data/test2_craters_h18.csv
carries the v2 hosting flag `hosts` (1 when any hollow is assigned to the crater, i.e. host_uid in host_assignments.csv;
707 hosts, checked to match exactly). For v3 and v3b, `hosts` is recomputed from code/rerun/out_v3/host_assignments_v3.csv
and _v3b.csv; every other column is copied unchanged. Writes inputs/test2_craters_h18_v3.csv and _v3b.csv."""
import os, pandas as pd
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
G = f"{ROOT}/manuscripts/GRL_insolation_hollows/data"
c = pd.read_csv(f"{G}/test2_craters_h18.csv")
h2 = set(pd.read_csv(f"{G}/host_assignments.csv", usecols=["host_uid"]).host_uid.dropna().astype(int))
assert set(c.unique_id[c.hosts == 1]) == h2, "v2 hosts flag does not equal the host_assignments host set"
for v in ["v3", "v3b"]:
    h = set(pd.read_csv(f"{ROOT}/code/rerun/out_v3/host_assignments_{v}.csv", usecols=["host_uid"]).host_uid.dropna().astype(int))
    o = c.copy(); o["hosts"] = o.unique_id.isin(h).astype(int)
    assert (o.drop(columns="hosts").equals(c.drop(columns="hosts")))
    o.to_csv(f"{HERE}/test2_craters_h18_{v}.csv", index=False)
    print(v, "hosts", int(o.hosts.sum()), "(v2 707); hosts lost", int(((c.hosts == 1) & (o.hosts == 0)).sum()), "gained", int(((c.hosts == 0) & (o.hosts == 1)).sum()))
