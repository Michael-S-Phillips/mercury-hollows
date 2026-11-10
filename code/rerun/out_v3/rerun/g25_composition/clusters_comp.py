#!/usr/bin/env python3
"""g25 copy of g20_regional_clusters/clusters.py (9 Nov 2026), with the composition-group joint model (comp_model.py: LRM + MDIS color
PC1, PC2 and a no-color indicator; coefficients comp_coefs.csv). Otherwise unchanged. Original docstring: g20 (9 Nov 2026), for Michael: where is hosting clustered beyond the joint model? Uses the paper's joint model (thinned presence,
complementary log-log, seven physical groups, g4_completeness/g1/summary_model/v3r/out/improved_coefs.csv 'thinned_presence_all')
to give each of the 31,554 craters a predicted probability P of hosting (>= 1 hollow within 0.9 R). On a 2-degree grid, counts observed
hosts O and expected hosts E = sum P within 500 km (the radius of the regional hosting fraction), z = (O - E) / sqrt(sum P(1 - P)).
Family-wise threshold: 1,000 data sets simulated from P (Bernoulli), 95th percentile of the maximum z (and of the minimum) over nodes
with E >= 2. Significant nodes are grouped into clusters (8-connected, longitude wraps). Writes grid.npz, clusters.json, craters_P.csv."""
import os, sys, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.ndimage import label
HE = os.path.expanduser("~/mnt/hollows-emissivity"); SM = f"{HE}/code/rerun/out_v3/rerun/g4_completeness/g1/summary_model"
HERE = os.path.dirname(os.path.abspath(__file__)); SUF = "" if os.environ.get("G20_RAD", "500") == "500" else "_r" + os.environ["G20_RAD"]; RM = 2439.7; RAD = float(os.environ.get("G20_RAD", "500")); NSIM = 1000; EMIN = 2.0
os.environ["SM_DIR"] = f"{SM}/v3r"; sys.path.insert(0, SM)
import common
common.DATA = f"{SM}/v3r/crater_features_xrs.csv"          # the file step7_improved.py fitted (same 31,554 rows)
d = common.load(); d["lnA"] = 2 * np.log(d.D); d["lnC"] = np.log(d.C1 + 0.01)
cf = pd.read_csv(f"{HERE}/comp_coefs.csv")
cz = pd.read_csv(f"{HERE}/crater_color_z.csv", dtype={"id": str}); d["id"] = d.id.astype(str); d = d.merge(cz, on="id", how="left"); assert d.pc1_z.notna().all()
eta = d.lnA.values.copy()
for t, b in zip(cf.term, cf.coef): eta += b * (1.0 if t == "Intercept" else d[t].values.astype(float))
P = -np.expm1(-np.exp(np.clip(eta, -30, 5))); y = d.h09.values.astype(float)
print("craters", len(d), "hosts", int(y.sum()), "sum P", round(P.sum(), 1))
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon); return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
glat = np.arange(-89, 90, 2.0); glon = np.arange(1, 360, 2.0); LO, LA = np.meshgrid(glon, glat)
tree = cKDTree(xyz(d.lat.values, d.lon.values % 360)); chord = 2 * np.sin(RAD / RM / 2)
nb = tree.query_ball_point(xyz(LA.ravel(), LO.ravel()), chord)
rows = np.repeat(np.arange(len(nb)), [len(x) for x in nb]); cols = np.concatenate([np.asarray(x, int) for x in nb])
M = sparse.csr_matrix((np.ones(len(cols), np.float32), (rows, cols)), shape=(len(nb), len(d)))
O = M @ y; E = M @ P; V = M @ (P * (1 - P)); Z = (O - E) / np.sqrt(np.maximum(V, 1e-9)); ok = E >= EMIN
rng = np.random.default_rng(20261109); zmax, zmin = [], []
for _ in range(NSIM // 100):
    Ys = (rng.random((len(d), 100)) < P[:, None]).astype(np.float32); Zs = ((M @ Ys) - E[:, None]) / np.sqrt(np.maximum(V, 1e-9))[:, None]
    zmax += list(Zs[ok].max(0)); zmin += list(Zs[ok].min(0))
thr_hi, thr_lo = float(np.percentile(zmax, 95)), float(np.percentile(zmin, 5))
print("FWER thresholds z >=", round(thr_hi, 2), "and <=", round(thr_lo, 2), "| observed max", round(Z[ok].max(), 2), "min", round(Z[ok].min(), 2))
shape = LA.shape; Zg = np.where(ok, Z, np.nan).reshape(shape)
hc = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv", usecols=["unique_id", "name"]); hc["unique_id"] = hc.unique_id.astype(str)
d = d.merge(hc, left_on="id", right_on="unique_id", how="left"); d["P"] = P
def clusters(mask, sign):
    lab, n = label(mask, structure=np.ones((3, 3)))
    if n and mask[:, 0].any() and mask[:, -1].any():                   # wrap in longitude
        for i in range(shape[0]):
            a, b = lab[i, 0], lab[i, -1]
            if a and b and a != b: lab[lab == b] = a
    out = []
    for k in np.unique(lab[lab > 0]):
        idx = np.flatnonzero((lab == k).ravel()); j = idx[np.argmax(sign * Z[idx])]
        inside = np.unique(np.concatenate([nb[i] for i in idx])); peak = np.asarray(nb[j], int)
        dd = d.iloc[peak]; hosts = dd[dd.h09 == 1].sort_values("n09", ascending=False)
        out.append(dict(peak_lat=float(LA.ravel()[j]), peak_lon=float(LO.ravel()[j]), z_peak=float(Z[j]), O_peak=int(O[j]), E_peak=round(float(E[j]), 1),
                        OE_peak=round(float(O[j] / E[j]), 2), n_nodes=int(len(idx)), lat_range=[float(LA.ravel()[idx].min()), float(LA.ravel()[idx].max())],
                        lon_values=sorted(set(float(v) for v in LO.ravel()[idx])), O_union=int(y[inside].sum()), E_union=round(float(P[inside].sum()), 1),
                        top_hosts=[dict(id=r.id, name=(r["name"] if isinstance(r["name"], str) else ""), D=round(r.D, 1), lat=round(r.lat, 1), lon=round(r.lon % 360, 1),
                                        n09=int(r.n09), P=round(r.P, 3)) for _, r in hosts.head(8).iterrows()]))
    return sorted(out, key=lambda c: -abs(c["z_peak"]))
hot = clusters(np.nan_to_num(Zg, nan=0) >= thr_hi, 1); cold = clusters(np.nan_to_num(Zg, nan=0) <= thr_lo, -1)
for c in hot + cold: c["lon_range_text"] = f"{min(c['lon_values']):.0f}-{max(c['lon_values']):.0f}E" if max(c['lon_values']) - min(c['lon_values']) < 300 else "wraps 0E"
np.savez(f"{HERE}/grid{SUF}.npz", glat=glat, glon=glon, O=O.reshape(shape), E=E.reshape(shape), Z=Z.reshape(shape), ok=ok.reshape(shape))
in_hot = np.zeros(len(d), bool)
for c in hot: pass
mask_hot = (np.nan_to_num(Zg, nan=0) >= thr_hi).ravel(); hot_cr = np.unique(np.concatenate([nb[i] for i in np.flatnonzero(mask_hot)])) if mask_hot.any() else np.array([], int)
summary = dict(created="2026-11-09", note=__doc__, n_craters=int(len(d)), n_hosts=int(y.sum()), sum_P=round(float(P.sum()), 1), radius_km=RAD, grid_deg=2, E_min=EMIN, nsim=NSIM,
               z_threshold_excess=round(thr_hi, 2), z_threshold_deficit=round(thr_lo, 2),
               hosts_within_500km_of_excess_nodes=int(y[hot_cr].sum()) if len(hot_cr) else 0, expected_there=round(float(P[hot_cr].sum()), 1) if len(hot_cr) else 0,
               excess=[{k: v for k, v in c.items() if k != "lon_values"} for c in hot], deficit=[{k: v for k, v in c.items() if k != "lon_values"} for c in cold])
json.dump(summary, open(f"{HERE}/clusters{SUF}.json", "w"), indent=1)
d[["id", "name", "lat", "lon", "D", "C1", "h09", "n09", "P"]].to_csv(f"{HERE}/craters_P{SUF}.csv", index=False)
for lab_, cs in (("EXCESS", summary["excess"]), ("DEFICIT", summary["deficit"])):
    print(lab_)
    for c in cs:
        print(f"  peak {c['peak_lat']:+.0f} {c['peak_lon']:.0f}E z={c['z_peak']:.1f} O/E={c['O_peak']}/{c['E_peak']} lat {c['lat_range']} lon {c['lon_range_text']} nodes={c['n_nodes']} union O/E={c['O_union']}/{c['E_union']}")
        print("     ", "; ".join(f"{h['name'] or h['id']} {h['D']}km ({h['lat']},{h['lon']}) n={h['n09']} P={h['P']}" for h in c["top_hosts"][:6]))
print("hosts near excess nodes", summary["hosts_within_500km_of_excess_nodes"], "expected", summary["expected_there"])
