#!/usr/bin/env python3
"""g4 (6 Nov 2026): compare every JSON / CSV record written by the copies here (out/<cat>/) with the g2 records
(g2_insolation/out/<cat>/). Usage: python3 check_repro.py v3 [v2 ...]. Writes out/repro_<cat>.json."""
import os, sys, json, filecmp
HERE = os.path.dirname(os.path.abspath(__file__)); G2 = os.path.abspath(os.path.join(HERE, "..", "..", "g2_insolation"))
def diff(a, b, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in set(a) | set(b):
            if k not in a or k not in b: out.append(f"{path}/{k}: missing")
            else: out += diff(a[k], b[k], f"{path}/{k}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): return [f"{path}: length {len(a)} vs {len(b)}"]
        return sum((diff(x, y, f"{path}[{i}]") for i, (x, y) in enumerate(zip(a, b))), [])
    if isinstance(a, float) and isinstance(b, float):
        return [] if (a == b or (a != a and b != b)) else [f"{path}: {a} vs {b}"]
    return [] if a == b else [f"{path}: {a} vs {b}"]
for cat in sys.argv[1:]:
    R = {}
    for f in sorted(os.listdir(f"{G2}/out/{cat}")):
        p1, p2 = f"{HERE}/out/{cat}/{f}", f"{G2}/out/{cat}/{f}"
        if not os.path.isfile(p2) or f.endswith((".png", ".pdf", ".txt")): continue
        if not os.path.exists(p1): R[f] = "missing"; continue
        if f.endswith(".json"):
            d = diff(json.load(open(p1)), json.load(open(p2))); R[f] = "identical" if not d else d[:10] + ([f"... {len(d)} differences"] if len(d) > 10 else [])
        else: R[f] = "identical" if filecmp.cmp(p1, p2, shallow=False) else "differs"
    json.dump(R, open(f"{HERE}/out/repro_{cat}.json", "w"), indent=1)
    print(cat, {k: v for k, v in R.items() if v != "identical"} or "all identical", len(R), "files")
