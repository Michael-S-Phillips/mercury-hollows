# g9 copy of g3_catalog/melt_source/knobs/t2_peak_hollows.py (6 Nov 2026). Unchanged; paths, catalog and grid come from common.py.
# g3_catalog copy of analysis_workshop/melt_source_tests/knobs/t2_peak_hollows.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Test 2b: completeness-weighted hollow density on a central-peak proxy for ALL craters 30-250 km (5 Nov 2026), unconditional on
hosting. Craters: Herrick et al. (2018) 30-250 km, plus Kinczyk-matched craters (kinczyk_matched.csv, 40-250 km) whose id is not a
Herrick crater in that range (their matched lat/lon/D are used). Regions: peak proxy rho <= 0.15 R; adjacent floor annulus
0.3-0.6 R; wall-and-terrace proxy 0.6-0.9 R; interior rho <= 0.9 R (the paper's host radius). Superposed craters (D >= 5 km, < 0.5 D_host, center within 1.05 R)
masked at 1.1 of their radii, with their hollows. Point grid spacing R/60, weighted by completeness for 1 km hollows (Dref 1.0).
Writes peak_hollows.csv. Usage: t2_peak_hollows.py START STOP"""
import sys, numpy as np, pandas as pd
from common import *
cr = load_craters(); h = load_hollows(); Cg = load_C(1.0)
kz = pd.read_csv(f"{AW}/kinczyk_degradation/kinczyk_matched.csv")
P = cr[(cr.diam >= 30) & (cr.diam <= 250)][["unique_id", "lat", "lon360", "diam", "cent_str", "preserved", "filled", "name"]].copy()
P["src"] = "herrick"
kx = kz[(kz.D >= 40) & (kz.D <= 250) & ~kz.id.isin(P.unique_id)]
P = pd.concat([P, pd.DataFrame(dict(unique_id=kx.id.values, lat=kx.lat.values, lon360=kx.lon.values % 360, diam=kx.D.values, cent_str="", preserved="", filled="", name="", src="kinczyk_only"))], ignore_index=True)
P = P.merge(kz[["id", "kclass", "kD"]].rename(columns={"id": "unique_id"}), on="unique_id", how="left").sort_values("unique_id").reset_index(drop=True)
a0, a1 = int(sys.argv[1]), int(sys.argv[2]); rows = []
for _, c in P.iloc[a0:a1].iterrows():
    R = c.diam / 2; la0, lo0 = c.lat, c.lon360; s = R / 60
    r, gla, glo = local_grid(la0, lo0, 0.9 * R, s); r = r / R; w = C_at(Cg, gla, glo)
    isup, _ = near(la0, lo0, cr.lat.values, cr.lon360.values, 1.05 * R)
    isup = isup[(cr.diam.values[isup] >= 5) & (cr.diam.values[isup] < 0.5 * c.diam)]
    ex = inside_any(gla, glo, cr.lat.values[isup], cr.lon360.values[isup], 1.1 * cr.diam.values[isup] / 2)
    ih, dh = near(la0, lo0, h.lat.values, h.lon360.values, 0.9 * R); dh = dh / R
    hex_ = inside_any(h.lat.values[ih], h.lon360.values[ih], cr.lat.values[isup], cr.lon360.values[isup], 1.1 * cr.diam.values[isup] / 2)
    row = dict(uid=int(c.unique_id), src=c.src, name=c["name"] if isinstance(c["name"], str) else "", D=c.diam, lat=la0, lon=lo0,
               cent_str=str(c.cent_str).strip() if isinstance(c.cent_str, str) else "", preserved=str(c.preserved).strip(), filled=str(c.filled).strip(),
               kclass=c.kclass, dcal_R=float(gcd(CAL["lat"], CAL["lon"], la0, lo0) / (CAL["D"] / 2)), n_sup=len(isup))
    for nm, (lo_, hi_) in {"pk": (0, 0.15), "fl": (0.3, 0.6), "wl": (0.6, 0.9), "int": (0, 0.9)}.items():
        g = (r >= lo_) & (r <= hi_) & ~ex; hh = (dh >= lo_) & (dh <= hi_) & ~hex_
        row[f"A_{nm}"] = float(g.sum() * s * s); row[f"Aw_{nm}"] = float(w[g].sum() * s * s); row[f"n_{nm}"] = int(hh.sum())
        f = h.field_id.values[ih][hh]; row[f"fields_{nm}"] = ";".join(str(int(x)) for x in sorted(set(f[pd.notna(f)])))
    rows.append(row)
os.makedirs(f"{HERE}/_chunks", exist_ok=True); out = f"{HERE}/_chunks/peak_hollows_{a0:05d}_{a1:05d}.csv"
pd.DataFrame(rows).to_csv(out, index=False); print(len(P), "craters total;", len(rows), "rows ->", out)
