#!/usr/bin/env python3
"""Stage labels of Deutsch et al. (2025) (canonical_stage_labels.csv, stage08) for hollows inside host craters (h09 rule,
rho <= 0.9; host_assignments.csv) that have a Kinczyk et al. (2020) class. Counts of hollows and of host craters with at
least one Stage 1 hollow, for all labels and for labels from images at incidence < 70 deg (as in Section S8). 5 Nov 2026.
Writes stage_by_class.json."""
import os, json, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../..")
L = pd.read_csv(f"{ROOT}/code/rerun/latitude/incidence_null/canonical_stage_labels.csv")
h = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", low_memory=False).merge(L[["thumb", "stage08"]], on="thumb", how="left")
k = pd.read_csv(f"{HERE}/kinczyk_matched.csv", dtype={"id": str}, low_memory=False); k = k[k.id.notna()][["id", "kclass"]]
h = h[h.rho <= 0.9].copy(); h["host_uid"] = h.host_uid.astype("Int64").astype(str)
m = h.merge(k, left_on="host_uid", right_on="id", how="inner")
R = {}
for lab, x in [("all", m), ("inc_lt_70", m[m.inc < 70])]:
    R[lab] = {"hollows_by_class_stage": {int(c): {int(s): int(n) for s, n in g.stage08.value_counts().items()} for c, g in x.groupby("kclass")},
              "hosts_with_stage1_by_class": {int(c): int(n) for c, n in x[x.stage08 == 1].groupby("kclass").host_uid.nunique().items()}}
    R[lab]["hosts_with_stage1_class1to3"] = int(x[(x.stage08 == 1) & (x.kclass <= 3)].host_uid.nunique())
json.dump(R, open(f"{HERE}/stage_by_class.json", "w"), indent=1); print(json.dumps({k_: v["hosts_with_stage1_class1to3"] for k_, v in R.items()}))
