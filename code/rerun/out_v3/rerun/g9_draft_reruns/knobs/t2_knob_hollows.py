# g9 copy of g3_catalog/melt_source/knobs/t2_knob_hollows.py (6 Nov 2026). Unchanged; paths, catalog and grid come from common.py.
# g3_catalog copy of analysis_workshop/melt_source_tests/knobs/t2_knob_hollows.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Test 2a: completeness-weighted hollow density on circum-Caloris knobs vs the adjacent matrix (5 Nov 2026).
All 2,011 Wright et al. (2020) knobs. Footprint: disk of radius a = long axis / 2 (3.5 km for unmeasured certain/uncertain knobs,
2.0 km for '<5 km' knobs); variant a + 0.5 km (as province_basins/knob_depths.py). Matrix: annulus 1.5a-3a excluding other knobs'
footprints and superposed Herrick craters 5-100 km (1.1 radii; a crater that contains the knob center is not excluded and the
knob is flagged in_crater). Areas from a local point grid (spacing a/15), each point weighted by catalog completeness for
1 km hollows (completeness_10ppd.npz, Dref = 1.0; row = (90 - lat) * 10, col = lon * 10, as peak_wall_test/peak_wall.py).
Hollows: hollows_canonical.csv without crater traces. Writes knob_hollows.csv (one row per knob) and, for regional context,
region_hollows.json: density over the Caloris exterior 1.0-2.0 basin radii (2 km grid) excluding craters 5-300 km and knobs."""
import json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from common import *
k = load_knobs(); cr = load_craters(); h = load_hollows(); Cg = load_C(1.0)
sup = cr[(cr.diam >= 5) & (cr.diam < 100)].reset_index(drop=True)
rows = []
for idx, q in k.iterrows():
    a = q.a_km; la0, lo0 = q.Latitude, q.lon360; s = max(a / 15, 0.12); rmax = 3 * a + 0.5
    r, gla, glo = local_grid(la0, lo0, rmax, s); w = C_at(Cg, gla, glo)
    isup, dsup = near(la0, lo0, sup.lat.values, sup.lon360.values, 56)
    host = isup[dsup <= 1.1 * sup.diam.values[isup] / 2]
    isup, _ = near(la0, lo0, sup.lat.values, sup.lon360.values, rmax + 56); isup = isup[~np.isin(isup, host)]
    ik, _ = near(la0, lo0, k.Latitude.values, k.lon360.values, rmax + 11); ik = ik[ik != idx]
    ih, dh = near(la0, lo0, h.lat.values, h.lon360.values, rmax)
    def ex_c(la, lo): return inside_any(la, lo, sup.lat.values[isup], sup.lon360.values[isup], 1.1 * sup.diam.values[isup] / 2)
    def ex_k(la, lo): return inside_any(la, lo, k.Latitude.values[ik], k.lon360.values[ik], k.a_km.values[ik])
    gc, gk = ex_c(gla, glo), ex_k(gla, glo)
    hc, hk = ex_c(h.lat.values[ih], h.lon360.values[ih]), ex_k(h.lat.values[ih], h.lon360.values[ih])
    reg = {"knob": (r <= a) & ~gc, "knobb": (r <= a + 0.5) & ~gc, "mat": (r >= 1.5 * a) & (r <= 3 * a) & ~gc & ~gk}
    hreg = {"knob": (dh <= a) & ~hc, "knobb": (dh <= a + 0.5) & ~hc, "mat": (dh >= 1.5 * a) & (dh <= 3 * a) & ~hc & ~hk}
    row = dict(knob_idx=int(idx), KnobID=q.KnobID, KnobType=q.KnobType, lat=la0, lon=lo0, long_axis_km=q.LongAxis_k, a_km=a, main=bool(q.main),
               dcal_R=q.dcal_R, in_crater=len(host) > 0, host_D=float(sup.diam.values[host].max()) if len(host) else 0.0)
    for nm in reg:
        row[f"A_{nm}"] = float(reg[nm].sum() * s * s); row[f"Aw_{nm}"] = float(w[reg[nm]].sum() * s * s)
        row[f"C_{nm}"] = float(w[reg[nm]].mean()) if reg[nm].any() else np.nan
        row[f"n_{nm}"] = int(hreg[nm].sum()); f = h.field_id.values[ih][hreg[nm]]
        row[f"fields_{nm}"] = ";".join(str(int(x)) for x in sorted(set(f[pd.notna(f)])))
    rows.append(row)
d = pd.DataFrame(rows); d.to_csv(f"{HERE}/knob_hollows.csv", index=False)
print(len(d), "knobs; hollows on knobs:", d.n_knob.sum(), "buffered:", d.n_knobb.sum(), "matrix:", d.n_mat.sum())
# regional context: Caloris exterior 1.0-2.0 R on a 2 km grid
RC = CAL["D"] / 2; s = 2.0
r, gla, glo = local_grid(CAL["lat"], CAL["lon"], 2 * RC, s); keep = r >= RC; r, gla, glo = r[keep], gla[keep], glo[keep]
P = xyz(gla, glo); tr = cKDTree(P); bad = np.zeros(len(P), bool)
cc = cr[(cr.diam >= 5) & (cr.diam < 300)]; cc = cc[gcd(CAL["lat"], CAL["lon"], cc.lat.values, cc.lon360.values) < 2 * RC + 170]
for lst in tr.query_ball_point(xyz(cc.lat.values, cc.lon360.values), r=2 * np.sin(1.1 * cc.diam.values / 2 / RM / 2)):
    bad[lst] = True
kk_ = np.zeros(len(P), bool)
for lst in tr.query_ball_point(xyz(k.Latitude.values, k.lon360.values), r=2 * np.sin(k.a_km.values / RM / 2)): kk_[lst] = True
w = C_at(Cg, gla, glo)
hd = gcd(CAL["lat"], CAL["lon"], h.lat.values, h.lon360.values); hm = (hd >= RC) & (hd <= 2 * RC)
hP = xyz(h.lat.values[hm], h.lon360.values[hm])
hbad = np.zeros(hm.sum(), bool); th = cKDTree(hP)
for lst in th.query_ball_point(xyz(cc.lat.values, cc.lon360.values), r=2 * np.sin(1.1 * cc.diam.values / 2 / RM / 2)): hbad[lst] = True
hk = np.zeros(hm.sum(), bool)
for lst in th.query_ball_point(xyz(k.Latitude.values, k.lon360.values), r=2 * np.sin(k.a_km.values / RM / 2)): hk[lst] = True
G = {}
for nm, (rlo, rhi) in {"1.0-2.0R": (1.0, 2.0), "1.0-1.5R": (1.0, 1.5), "1.5-2.0R": (1.5, 2.0)}.items():
    pm = (r >= rlo * RC) & (r <= rhi * RC); hr = (hd[hm] >= rlo * RC) & (hd[hm] <= rhi * RC)
    sel = pm & ~bad & ~kk_; hs = hr & ~hbad & ~hk
    G[nm] = dict(area_km2=float(sel.sum() * s * s), area_w_km2=float(w[sel].sum() * s * s), mean_C=float(w[sel].mean()), n_hollows=int(hs.sum()),
                 n_fields=int(pd.Series(h.field_id.values[hm][hs]).nunique()), density_per_1000km2_w=float(1000 * hs.sum() / (w[sel].sum() * s * s)),
                 all_terrain_incl_craters=dict(area_w_km2=float(w[pm & ~kk_].sum() * s * s), n_hollows=int((hr & ~hk).sum())))
json.dump(G, open(f"{HERE}/region_hollows.json", "w"), indent=1); print(json.dumps(G, indent=1))
