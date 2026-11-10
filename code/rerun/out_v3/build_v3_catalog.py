#!/usr/bin/env python3
"""Screened hollow catalogs v3 and v3b and their host assignments.

Written for the archive (November 2026). The catalog files were first made by commands run inline, without a saved script;
this script applies the rule recorded in out_v3/README.md and reproduces the archived files exactly:
  is_crater_trace (v3)  = v2 crater trace (D_eq > 10 km) OR screen class in {B, W, O}   (crater bowl, wall/rim/slope, other)
  is_crater_trace (v3b) = the v3 flag OR screen class U (uncertain)
  host assignments      = the v2 host assignments without the hollows flagged above (assignment is per hollow).
Inputs: code/rerun/out_v2/hollows_canonical.csv, manuscripts/GRL_insolation_hollows/data/host_assignments.csv,
code/rerun/out_v3/screen_classes.csv (visual classes; see out_v3/README.md and qc/).
Outputs (code/rerun/out_v3/): hollows_canonical_v3.csv, hollows_canonical_v3b.csv, host_assignments_v3.csv, host_assignments_v3b.csv.
Usage: python3 code/rerun/out_v3/build_v3_catalog.py [output folder, default: this folder]
"""
import os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
v2 = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv")
ha = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv")
sc = pd.read_csv(f"{HERE}/screen_classes.csv")[["hollow_id", "class"]]
cls = v2.hollow_id.map(dict(zip(sc.hollow_id, sc["class"])))
for tag, bad in (("v3", {"B", "W", "O"}), ("v3b", {"B", "W", "O", "U"})):
    h = v2.copy()
    h["screen_class"] = cls
    h["is_spurious_v3"] = cls.isin(["B", "W", "O"])
    h["is_uncertain_v3"] = cls.eq("U")
    h["is_crater_trace_v2"] = v2.is_crater_trace
    h["is_crater_trace"] = v2.is_crater_trace | cls.isin(bad)
    h.to_csv(f"{OUT}/hollows_canonical_{tag}.csv", index=False)
    drop = set(h.hollow_id[h.is_crater_trace])
    ha[~ha.hollow_id.isin(drop)].to_csv(f"{OUT}/host_assignments_{tag}.csv", index=False)
    print(tag, "non-trace hollows", int((~h.is_crater_trace).sum()), "host rows", int((~ha.hollow_id.isin(drop)).sum()))
