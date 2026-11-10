#!/usr/bin/env python3
"""g17 (7 Nov 2026), new script: one file per crater for the tracer page, pk/<uid>.json = {img, dem, col, ecol} as data URLs,
so that 194 more craters stay within the artifact's limit on published files. Reads assets/; writes pk/."""
import os, json, base64
H = os.path.dirname(os.path.abspath(__file__)); A = f"{H}/assets"; os.makedirs(f"{H}/pk", exist_ok=True)
SEL = json.load(open(f"{H}/selection.json")); U = [g["uid"] for k in ("group_B", "group_C", "group_D") for g in SEL[k]]; tot = 0
for u in U:
    d = {}
    for k, ext, mt in (("img", "jpg", "image/jpeg"), ("dem", "png", "image/png"), ("col", "jpg", "image/jpeg"), ("ecol", "jpg", "image/jpeg")):
        d[k] = f"data:{mt};base64," + base64.b64encode(open(f"{A}/{k}/{u}.{ext}", "rb").read()).decode()
    s = json.dumps(d, separators=(",", ":")); open(f"{H}/pk/{u}.json", "w").write(s); tot += len(s)
print("packs", len(U), "MB", round(tot / 1e6, 1))
