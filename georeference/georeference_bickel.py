#!/usr/bin/env python3
"""Recompute Bickel et al. (2024/2025) hollow-detection positions through the MDIS camera model with ISIS.

RECONSTRUCTED FOR THE ARCHIVE (29 Oct 2026). The positions in GIS/vector/bickel2024_hollows_corrected.gpkg were computed on
11 Sep 2026 with the procedure below (see GIS/vector/bickel2024_hollows_corrected_README.md), but the per-image driver script
was not saved. This script re-implements that procedure. Before it is cited as the code that produced the archived positions,
run it on a sample of images and check the result with compare_to_archive.py (median difference should be < 0.05 km).

Procedure for each MDIS narrow-angle image in the 0.7-confidence release:
  1. Download the EDR (PDS Imaging Node, MSGRMDS_1001/DATA/YYYY_DOY/EN<id>.IMG). Bickel IDs are CDR names (CN<id>);
     the EDR is EN<id>.
  2. mdis2isis -> spiceinit (reconstructed SPK/CK, system shape model = USGS MESSENGER DEM).
  3. campt at the detection bounding-box center and its four corners. The release's 'Image x/y coord min/max (pix)' columns
     are full-frame coordinates normalized to 0-1 in raw EDR orientation, so sample = x * samples + 0.5, line = y * lines + 0.5
     (ISIS pixel centers are at integer coordinates).
Requires ISIS 8 on PATH with ISISROOT and ISISDATA set (MESSENGER kernels installed), and network access to the PDS.

Usage: python3 georeference_bickel.py <Bickel_et_al_2024_MDIS_hollows_07_v1.0.csv> <work_dir> [--images ID ...] [--max N]
Writes <work_dir>/recomputed_positions.csv (one row per detection) and keeps the ISIS cubes in <work_dir>/cubes/.
"""
import argparse, csv, io, os, subprocess, sys, urllib.request
from datetime import datetime
import pandas as pd

EDR_BASE = "https://pds-imaging.jpl.nasa.gov/data/messenger/msgrmds_1001/DATA"

def read_release(path):
    raw = open(path, "rb").read().decode("latin-1")
    rows = list(csv.reader(io.StringIO(raw), delimiter=";"))
    hdr = [c.strip() for c in rows[0]]; ix = {c: i for i, c in enumerate(hdr)}
    out = []
    for k, r in enumerate(rows[1:]):
        if len(r) < len(hdr): continue
        g = lambda name: r[ix[name]]
        out.append(dict(row=k, mdis_id=g("MDIS ID"), earth_datetime=g("Earth datetime"),
                        xmin=float(g("Image x coord min (pix)")), xmax=float(g("Image x coord max (pix)")),
                        ymin=float(g("Image y coord min (pix)")), ymax=float(g("Image y coord max (pix)")),
                        lat_orig=float(g("Detection center LAT")), lon180_orig=float(g("Detection center LON 180d"))))
    return pd.DataFrame(out)

def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0: raise RuntimeError(f"{' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p.stdout

def get_edr(mdis_id, when, edr_dir):
    edr = "EN" + mdis_id[2:]
    dst = os.path.join(edr_dir, edr + ".IMG")
    if os.path.exists(dst): return dst
    t = datetime.fromisoformat(when.split(".")[0])
    for day in (t.timetuple().tm_yday, t.timetuple().tm_yday - 1, t.timetuple().tm_yday + 1):   # VERIFY: directory date convention
        url = f"{EDR_BASE}/{t.year}_{day:03d}/{edr}.IMG"
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(dst, "wb") as f: f.write(r.read())
            return dst
        except Exception:
            continue
    raise FileNotFoundError(f"EDR for {mdis_id} not found under {EDR_BASE}/{t.year}_DOY/")

def label_value(lbl, key):
    for line in lbl.splitlines():
        if line.strip().startswith(key):
            return line.split("=", 1)[1].strip()
    return None

def georeference_image(mdis_id, dets, work):
    cub = os.path.join(work, "cubes", "EN" + mdis_id[2:] + ".cub")
    if not os.path.exists(cub):
        edr = get_edr(mdis_id, dets.earth_datetime.iloc[0], os.path.join(work, "edr"))
        run(["mdis2isis", f"from={edr}", f"to={cub}"])
        run(["spiceinit", f"from={cub}", "shape=system"])
    dims = run(["catlab", f"from={cub}"])
    ns = int(label_value(dims, "Samples")); nl = int(label_value(dims, "Lines"))
    pts = []
    for _, d in dets.iterrows():
        cx, cy = (d.xmin + d.xmax) / 2, (d.ymin + d.ymax) / 2
        for tag, x, y in (("c", cx, cy), ("ul", d.xmin, d.ymin), ("ur", d.xmax, d.ymin), ("lr", d.xmax, d.ymax), ("ll", d.xmin, d.ymax)):
            pts.append((d.row, tag, x * ns + 0.5, y * nl + 0.5))
    lst = os.path.join(work, "tmp_coords.txt"); flat = os.path.join(work, "tmp_campt.csv")
    with open(lst, "w") as f:
        for _, _, s, l in pts: f.write(f"{s:.3f},{l:.3f}\n")
    if os.path.exists(flat): os.remove(flat)
    run(["campt", f"from={cub}", "usecoordlist=yes", f"coordlist={lst}", "coordtype=image", "format=flat", f"to={flat}", "allowoutside=yes"])
    c = pd.read_csv(flat)
    lat = c["PlanetocentricLatitude"].values; lon = c["PositiveEast360Longitude"].values
    out = {}
    for (row, tag, _, _), la, lo in zip(pts, lat, lon):
        out.setdefault(row, {})[f"lat_{tag}"] = la; out[row][f"lon180_{tag}"] = ((lo + 180) % 360) - 180
    res = pd.DataFrame.from_dict(out, orient="index"); res.index.name = "row"
    res["frame_samples"], res["frame_lines"] = ns, nl
    return res

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("release_csv"); ap.add_argument("work_dir")
    ap.add_argument("--images", nargs="*"); ap.add_argument("--max", type=int)
    a = ap.parse_args()
    for sub in ("cubes", "edr"): os.makedirs(os.path.join(a.work_dir, sub), exist_ok=True)
    rel = read_release(a.release_csv)
    ids = a.images or sorted(rel.mdis_id.unique())
    if a.max: ids = ids[:a.max]
    parts, failed = [], []
    for k, mid in enumerate(ids):
        try:
            parts.append(georeference_image(mid, rel[rel.mdis_id == mid], a.work_dir)); print(f"[{k+1}/{len(ids)}] {mid} ok")
        except Exception as e:
            failed.append(mid); print(f"[{k+1}/{len(ids)}] {mid} FAILED: {str(e)[:200]}", file=sys.stderr)
    res = rel.set_index("row").join(pd.concat(parts), how="inner") if parts else rel.iloc[:0]
    res.rename(columns={"lat_c": "lat", "lon180_c": "lon180"}).to_csv(os.path.join(a.work_dir, "recomputed_positions.csv"))
    print(f"done: {len(ids) - len(failed)} images, {len(failed)} failed")

if __name__ == "__main__":
    main()
