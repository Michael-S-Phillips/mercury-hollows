#!/usr/bin/env python3
"""g2: compare every v2 output of the copied scripts with the stored record of the original script (reproduction check).
Writes out/v2_reproduction.json."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
PR = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"; W6 = f"{PR}/analysis_workshop/insolation_v6"; G = f"{ROOT}/manuscripts/GRL_insolation_hollows/records"
PAIRS = {"dT_record.json": f"{PR}/analysis_v3/dT_record.json", "dT_plains_record.json": f"{PR}/analysis_v3/dT_plains_record.json",
         "dT_hosting_completeness_record.json": f"{PR}/analysis_v7/dT_hosting_completeness_record.json",
         "transition_layer_hilat_record.json": f"{PR}/analysis_v3/transition_layer_hilat_record.json",
         "latcorr_density_record.json": f"{G}/latcorr_density_record.json", "latcorr_inc60_75_record.json": f"{G}/latcorr_inc60_75_record.json",
         "hotwarm_record.json": f"{ROOT}/code/rerun/latitude/hotwarm/hotwarm_record.json", "aspect_record.json": f"{ROOT}/code/rerun/latitude/aspect/aspect_record.json"}
for n in ("cos2lon_fit_lat15", "size_quirk", "band_crater_hotwarm", "hollow_rich_hosts", "fractions", "threshold", "tm_vs_dT_check",
          "stage_test", "stage_test_lowinc", "stage_fieldlevel"):
    PAIRS[f"{n}.json"] = f"{W6}/{n}.json"
def cmp(a, b, path=""):
    bad = []
    if isinstance(a, dict):
        for k in a:
            if k == "description": continue
            bad += cmp(a[k], b[k], f"{path}/{k}") if k in b else [f"{path}/{k} missing"]
    elif isinstance(a, list):
        bad += [f"{path} length {len(a)} vs {len(b)}"] if len(a) != len(b) else []
        for i, (x, y) in enumerate(zip(a, b)): bad += cmp(x, y, f"{path}[{i}]")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not (a == b or abs(a - b) <= 1e-6 * max(1.0, abs(a))): bad.append(f"{path}: record {a} vs rerun {b}")
    elif a != b: bad.append(f"{path}: record {a} vs rerun {b}")
    return bad
out = {}
for mine, rec in PAIRS.items():
    d = cmp(json.load(open(rec)), json.load(open(f"{HERE}/out/v2/{mine}")))
    out[mine] = {"record": os.path.relpath(rec, ROOT), "n_differences": len(d), "first_differences": d[:6]}
    print(f"{mine:38s} {'MATCH' if not d else str(len(d)) + ' differences: ' + '; '.join(d[:3])}")
json.dump(out, open(f"{HERE}/out/v2_reproduction.json", "w"), indent=1)
