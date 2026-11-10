# g13 copy of g11_curvature/positions/merge_chunks.py (6 Nov 2026). Changed: folder g13_traced/positions; G11_RULES default traced (traced outlines where available, curvature v3 elsewhere via g13_traced/traced_classes.py; env CLS=v2 for curvature v2 elsewhere).
# g11 copy of g9_draft_reruns/positions/merge_chunks.py (6 Nov 2026). Changed: also joins the eta_b arrays (peak_etab, hollows_etab, jitter_etab, fields_etab; slope: peak_etab, hollows_etab, fields_etab, zbase_rel) and the slope hollows/fields reshape widths stay 6 and 7.
#!/usr/bin/env python3
"""g9 (6 Nov 2026), new script: joins the START-STOP chunks of positions_build.py (cache/H{H}_{a}_{b}.npz) into cache/H{H}_000_354.npz,
the single file positions_analyze.py reads, and the chunks of positions_slope_build.py (cache/_chunks/slope_H250_{a}_{b}.npz) into
cache/slope_H250.npz. Chunks must cover the host list contiguously from 0; this is asserted.
Usage: merge_chunks.py positions ODIR H | merge_chunks.py slope ODIR"""
import sys, glob, re, numpy as np, pandas as pd
kind, odir = sys.argv[1], sys.argv[2]
if kind == "positions":
    H = sys.argv[3]; fs = [f for f in glob.glob(f"{odir}/cache/H{H}_*.npz") if not f.endswith("_000_354.npz")]
    pat = re.compile(rf"H{H}_(\d+)_(\d+)\.npz$")
else:
    fs = glob.glob(f"{odir}/cache/_chunks/slope_H250_*.npz"); pat = re.compile(r"slope_H250_(\d+)_(\d+)\.npz$")
se = sorted((int(pat.search(f).group(1)), int(pat.search(f).group(2)), f) for f in fs)
assert se[0][0] == 0 and all(a[1] == b[0] for a, b in zip(se[:-1], se[1:])), [(a, b) for a, b, _ in se]
Z = [np.load(f, allow_pickle=True) for _, _, f in se]
if kind == "positions":
    cr = pd.concat([pd.DataFrame(z["craters"]) for z in Z], ignore_index=True)
    out = dict(craters=cr.to_records(index=False), peak=np.vstack([z["peak"] for z in Z]), floor=np.vstack([z["floor"] for z in Z]),
               hollows=np.vstack([z["hollows"] for z in Z]), jitter=np.concatenate([z["jitter"] for z in Z], axis=2),
               fields=np.vstack([z["fields"].reshape(-1, 13) for z in Z]),
               peak_etab=np.concatenate([z["peak_etab"] for z in Z]), hollows_etab=np.concatenate([z["hollows_etab"] for z in Z]),
               jitter_etab=np.concatenate([z["jitter_etab"] for z in Z], axis=1), fields_etab=np.concatenate([z["fields_etab"] for z in Z]))
    dst = f"{odir}/cache/H{H}_000_354.npz"
else:
    out = dict(craters=np.concatenate([z["craters"] for z in Z]), peak=np.vstack([z["peak"] for z in Z]),
               hollows=np.vstack([z["hollows"].reshape(-1, 6) for z in Z]), fields=np.vstack([z["fields"].reshape(-1, 7) for z in Z]),
               peak_etab=np.concatenate([z["peak_etab"] for z in Z]), hollows_etab=np.concatenate([z["hollows_etab"] for z in Z]),
               fields_etab=np.concatenate([z["fields_etab"] for z in Z]), zbase_rel=np.concatenate([z["zbase_rel"] for z in Z]))
    dst = f"{odir}/cache/slope_H250.npz"
np.savez_compressed(dst, **out); print("merged", [(a, b) for a, b, _ in se], "->", dst, {k: v.shape for k, v in out.items()})
