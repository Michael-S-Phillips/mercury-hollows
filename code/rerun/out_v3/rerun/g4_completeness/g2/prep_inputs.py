#!/usr/bin/env python3
"""g4 (6 Nov 2026): inputs of the g2 copies in this folder, built like g2_insolation/prep_inputs.py.
  inputs/v2, v3, v3b : copies of g2_insolation/inputs/<cat>/ (old completeness grid; used to check that the copies here
                       reproduce the g2 records)
  inputs/v3r, v3br   : the files of inputs/v3 (v3b) with only the completeness values replaced, from the grid built with the
                       refit detection model (../completeness_10ppd_v3.npz, _v3b.npz):
     hollows_canonical.csv, hollows_canonical_inc.csv, fields.csv : cov_*, C_D0.5km, C_D1km, C_D2km, C_D5km, C_own from
                       ../build_completeness.py <v3r|v3br> attach (verbatim code of build_hollow_products.py)
     crater_features.csv : C1 resampled at the crater centers exactly as crater_dag.py samples it
     crater_dag_nodes_C1.csv : id, C1 at the analysis_v9 graph nodes (read by dT_hosting_completeness.py instead of the
                       C1 column of analysis_v9/crater_dag_nodes.csv)
  Every other column and file is byte-identical to inputs/v3 (v3b). Checks: the same sampling on the old grid returns the
  stored C1 of crater_features.csv and of analysis_v9/crater_dag_nodes.csv exactly. Writes inputs/prep_g4_record.json."""
import os, json, shutil, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); G4 = os.path.dirname(HERE); RT = os.path.abspath(os.path.join(G4, *[".."] * 5))
G2 = f"{RT}/code/rerun/out_v3/rerun/g2_insolation"; PR = f"{RT}/manuscripts/PlanetaryResearch_hollows"
def C1grid(path):
    Cz = np.load(path); return Cz["C"][list(Cz["Dref"]).index(1.0)]
def samp_dag(C1, lat, lon):          # crater_dag.py: rr/cc as written there
    rr = ((90 - lat) * 10).astype(int).clip(0, 1799); cc = (lon * 10).astype(int).clip(0, 3599)
    return np.nan_to_num(C1[rr, cc])
OLD = C1grid(f"{RT}/code/rerun/out_v2/completeness_10ppd.npz")
R = {}
for cat in ("v2", "v3", "v3b"):
    D = f"{HERE}/inputs/{cat}"; os.makedirs(D, exist_ok=True)
    for f in os.listdir(f"{G2}/inputs/{cat}"): shutil.copyfile(f"{G2}/inputs/{cat}/{f}", f"{D}/{f}")
NV9 = pd.read_csv(f"{PR}/analysis_v9/crater_dag_nodes.csv", dtype={"id": str}, low_memory=False, usecols=["id", "lat", "lon", "C1"])
R["check_nodes_C1_old_grid_identical_float32"] = bool((samp_dag(OLD, NV9.lat.values, NV9.lon.values).astype(np.float32) == NV9.C1.values.astype(np.float32)).all())
for cat, base, grid in (("v3r", "v3", "completeness_10ppd_v3.npz"), ("v3br", "v3b", "completeness_10ppd_v3b.npz")):
    D = f"{HERE}/inputs/{cat}"; S = f"{HERE}/inputs/{base}"; A = f"{G4}/_scratch/attach_{cat}"; os.makedirs(D, exist_ok=True)
    NEW = C1grid(f"{G4}/{grid}"); r = {}
    for f in ("canonical_stage_labels.csv", "host_assignments.csv", "test2_craters_h18.csv"): shutil.copyfile(f"{S}/{f}", f"{D}/{f}")
    rd = lambda p: pd.read_csv(p, low_memory=False, float_precision="round_trip")
    # hollows and fields: completeness columns from the attach step
    for f in ("hollows_canonical", "fields"):
        b = rd(f"{S}/{f}.csv"); a = rd(f"{A}/{f}.csv"); key = "hollow_id" if f == "hollows_canonical" else "field_id"
        assert (a[key].values == b[key].values).all()
        cols = [c for c in b.columns if c.startswith("cov_") or c.startswith("C_")]
        assert cols == [c for c in a.columns if c.startswith("cov_") or c.startswith("C_")], (cols, a.columns)
        o = b.copy()
        for c in cols: o[c] = a[c].values
        o.to_csv(f"{D}/{f}.csv", index=False)
        chk = rd(f"{D}/{f}.csv"); assert all(chk[c].equals(b[c]) for c in b.columns if c not in cols)
        r[f] = {c: dict(median_ratio=float(np.nanmedian(o[c] / b[c])), max_abs_diff=float(np.nanmax(np.abs(o[c] - b[c])))) for c in cols if c.startswith("C_")}
        if f == "hollows_canonical":
            nt = ~o.is_crater_trace.astype(bool); o[nt].to_csv(f"{D}/hollows_canonical_inc.csv", index=False)
            bi = rd(f"{S}/hollows_canonical_inc.csv"); assert len(bi) == int(nt.sum()) and (bi.hollow_id.values == o[nt].hollow_id.values).all()
    # crater features: C1 resampled
    cf = rd(f"{S}/crater_features.csv"); cf_id = pd.read_csv(f"{S}/crater_features.csv", usecols=["id"], dtype={"id": str})
    old = samp_dag(OLD, cf.lat.values, cf.lon.values); r["check_crater_features_C1_old_grid_identical_float32"] = bool((old.astype(np.float32) == cf.C1.values.astype(np.float32)).all())
    assert r["check_crater_features_C1_old_grid_identical_float32"]
    new = samp_dag(NEW, cf.lat.values, cf.lon.values).astype(np.float64)
    raw = open(f"{S}/crater_features.csv").read().split("\n"); hdr = raw[0].split(","); k = hdr.index("C1")
    # write by replacing the C1 field only, so every other field keeps its exact text
    out = [raw[0]]
    for i, line in enumerate(raw[1:]):
        if not line: continue
        p = line.split(","); assert len(p) == len(hdr), i; p[k] = str(np.float32(new[i])); out.append(",".join(p))
    open(f"{D}/crater_features.csv", "w").write("\n".join(out) + "\n")
    chk = rd(f"{D}/crater_features.csv"); assert np.allclose(chk.C1.values, new, atol=1e-7) and all(chk[c].equals(cf[c]) for c in cf.columns if c != "C1")
    m = cf.C1.values > 0
    r["crater_features_C1"] = dict(n=int(len(cf)), mean_old=float(cf.C1.mean()), mean_new=float(new.mean()),
                                   ratio_p5_p50_p95=[float(x) for x in np.percentile(new[m] / cf.C1.values[m], [5, 50, 95])],
                                   h09_hosts_ratio_p5_p50_p95=[float(x) for x in np.percentile((new / np.where(cf.C1 > 0, cf.C1, 1))[m & (cf.h09.values == 1)], [5, 50, 95])])
    pd.DataFrame({"id": NV9.id, "C1": samp_dag(NEW, NV9.lat.values, NV9.lon.values)}).to_csv(f"{D}/crater_dag_nodes_C1.csv", index=False)
    R[cat] = r; print(cat, json.dumps(r["crater_features_C1"]), r["hollows_canonical"]["C_D1km"], r["hollows_canonical"]["C_own"])
json.dump(R, open(f"{HERE}/inputs/prep_g4_record.json", "w"), indent=1); print("nodes check", R["check_nodes_C1_old_grid_identical_float32"])
