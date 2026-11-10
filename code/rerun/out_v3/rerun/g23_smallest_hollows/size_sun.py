#!/usr/bin/env python3
"""g23 step 9 (9 Nov): does sunlight set hollow size? Michael: hollows form with the crater and keep growing by sun-driven activity such as
scarp retreat, so poleward hollows are smaller, and smaller hollows are shallower. Catalog hollows (v3, non-trace), equivalent diameter
of the finest detection, with ln px of that detection controlled (detection and outlining depend on it). Insolation proxies: poleward
of 60N; cos(2 lon) (+1 at the hot longitudes 0/180, -1 at the warm 90/270). Repeated without the province (500 km of 65N 125E).
Writes size_sun.json."""
import json, numpy as np, pandas as pd
H = pd.read_csv("../g2_insolation/inputs/v3/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool) & (H.px_best <= 40)].copy()
def xyz(la, lo): la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
H["prov"] = (xyz(H.lat, H.lon) @ xyz(np.array([65.0]), np.array([125.0]))[0]) >= np.cos(500 / 2439.7)
H["lD"] = np.log(H.D_eq_km); H["lpx"] = np.log(H.px_best); H["pole"] = (H.lat >= 60).astype(float); H["hot"] = np.cos(np.radians(2 * H.lon)); H["north"] = (H.lat > 0).astype(float)
out = {}
def fit(d, cols, tag):
    X = np.column_stack([np.ones(len(d))] + [d[c] for c in cols]); y = d.lD; b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b; cov = r.var(ddof=X.shape[1]) * np.linalg.inv(X.T @ X); se = np.sqrt(np.diag(cov))
    res = {c: [round(float(np.exp(b[i + 1])), 3), round(float(np.exp(b[i + 1] - 1.96 * se[i + 1])), 3), round(float(np.exp(b[i + 1] + 1.96 * se[i + 1])), 3)] for i, c in enumerate(cols)}
    out[tag] = dict(n=len(d), size_ratio=res); print(tag, len(d), res)
fit(H, ["lpx", "pole", "hot"], "all")
fit(H[~H.prov], ["lpx", "pole", "hot"], "no_province")
eq = H[(H.lat.abs() < 45) & ~H.prov]; fit(eq, ["lpx", "hot"], "lat_lt45_no_province")
for nm, d in [("all", H), ("no_province", H[~H.prov])]:
    m = d.assign(pb=pd.cut(d.px_best, [0, 20, 40])).groupby(["pb", "pole"], observed=True).D_eq_km.median().unstack().round(3)
    print(nm, "median D by px bin and poleward:"); print(m.to_string()); out[f"median_{nm}"] = {str(k): {str(kk): vv for kk, vv in v.items()} for k, v in m.astype(float).to_dict().items()}
print("province share of poleward hollows", round(float(H[H.pole == 1].prov.mean()), 3))
json.dump(out, open("size_sun.json", "w"), indent=1, default=str)
