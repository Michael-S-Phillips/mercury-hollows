#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script. Copies code/rerun/latitude/aspect_robustness/aspect_robustness.py into aspect/ and changes only
# input/output paths (each asserted once): project root (the original hard-codes the Mac path), host assignments (G6_HA), catalog (G6_CANON),
# DEM (G6_DEM, in place of dempath.txt) and output folder (G6_OUT).
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); SRC = f"{HE}/code/rerun/latitude/aspect_robustness/aspect_robustness.py"
HERE = os.path.dirname(os.path.abspath(__file__)); t = open(SRC).read()
reps = [('import json, numpy as np, pandas as pd', 'import os, json, numpy as np, pandas as pd'),
        ('BASE = "/Users/phillipsm/Documents/Research/hollows-emissivity"', 'BASE = os.path.expanduser("~/mnt/hollows-emissivity")\nOUTD = os.environ["G6_OUT"]; os.makedirs(OUTD, exist_ok=True)'),
        ('HA = f"{BASE}/code/rerun/latitude/aspect/host_assignments.csv"', 'HA = os.environ.get("G6_HA", f"{BASE}/code/rerun/latitude/aspect/host_assignments.csv")'),
        ('CANON = f"{BASE}/hollow-networks/data_v2/hollows_canonical.csv"', 'CANON = os.environ.get("G6_CANON", f"{BASE}/hollow-networks/data_v2/hollows_canonical.csv")'),
        ('DEMPATH = open(f"{BASE}/code/rerun/latitude/aspect_robustness/dempath.txt").read().strip()', 'DEMPATH = os.environ["G6_DEM"]'),
        ('with open(f"{BASE}/code/rerun/latitude/aspect_robustness/aspect_robustness_record.json", "w") as f:', 'with open(f"{OUTD}/aspect_robustness_record.json", "w") as f:'),
        ('wall.to_csv(f"{BASE}/code/rerun/latitude/aspect_robustness/wall_with_dem_aspect.csv", index=False)', 'wall.to_csv(f"{OUTD}/wall_with_dem_aspect.csv", index=False)')]
for a, b in reps:
    assert t.count(a) == 1, a; t = t.replace(a, b)
head, t = t.split("\n", 1)
open(f"{HERE}/aspect/aspect_robustness.py", "w").write(head + "\n# g6_remaining copy of code/rerun/latitude/aspect_robustness/aspect_robustness.py (6 Nov 2026). Changed: input/output paths only (see make_copies_aspect.py).\n" + t)
print("aspect/aspect_robustness.py <- code/rerun/latitude/aspect_robustness/aspect_robustness.py (%d replacements)" % len(reps))
