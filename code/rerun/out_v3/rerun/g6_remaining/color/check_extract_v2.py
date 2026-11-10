#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script. v2 check of the copied extract_tiles.py on selected tiles.
"""Runs process_tile of the copied extract_tiles.py with the v2 mask (code/rerun/out_v2/hollows_canonical.csv, the default) for
the tiles named on the command line, writing to _scratch/tiles_v2check, and compares each output array with the original
analysis_workshop/color_mdr/tiles_qc/ files (exact equality), then empties the compared outputs (files cannot be deleted here). 'all' runs every tile not yet checked
(resumable). Appends the result to out_v2/extract_check.json."""
import os, sys, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); G6 = os.path.dirname(HERE)
os.environ["G6_TILES"] = f"{G6}/_scratch/tiles_v2check"; os.makedirs(os.environ["G6_TILES"], exist_ok=True)
assert "G6_MASK" not in os.environ
TILES = sys.argv[1:]; sys.argv = sys.argv[:1]
T0 = time.time(); BUDGET = 150.0
import extract_tiles as E
ORIG = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/color_mdr/tiles_qc")
F = E.build_features(); EX = E.build_exclusion()
import pandas as pd
fo = pd.read_csv(f"{ORIG}/features.csv"); same_F = bool(np.allclose(F.values.astype(float), fo.values.astype(float), equal_nan=True))
os.makedirs(f"{HERE}/out_v2", exist_ok=True); P = f"{HERE}/out_v2/extract_check.json"
R = json.load(open(P)) if os.path.exists(P) else {}
R["features_identical"] = same_F; R["n_exclusion_disks"] = int(len(EX[0]))
if TILES == ["all"]: TILES = [n for n in E.g.all_tiles() if n not in R]   # resumable: tiles already checked are skipped
for n in TILES:
    if time.time() - T0 > BUDGET: break
    t0 = time.time(); info = E.process_tile(n, F, EX); res = {"seconds": round(time.time() - t0, 1)}
    for kind in ["complete", "partial"]:
        a, b = f"{E.OUT}/{n}_{kind}.npz", f"{ORIG}/{n}_{kind}.npz"
        if not os.path.exists(b): res[kind] = "no original"; continue
        za, zb = np.load(a), np.load(b)
        res[kind] = {k: bool(np.array_equal(za[k], zb[k], equal_nan=(za[k].dtype.kind == "f"))) for k in zb.files}
    for kind in ["complete", "partial"]:   # the outputs are identical copies of the originals; empty them to save disk space
        a = f"{E.OUT}/{n}_{kind}.npz"
        if os.path.exists(a): open(a, "w").close()
    R[n] = res; print(n, res, flush=True); json.dump(R, open(P, "w"), indent=1)
json.dump(R, open(P, "w"), indent=1)
