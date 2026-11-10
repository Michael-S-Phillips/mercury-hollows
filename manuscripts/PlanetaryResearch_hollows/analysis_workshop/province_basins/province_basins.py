#!/usr/bin/env python3
"""Is the hollow-rich province at 90-150E, 40-70N tied to the old basins b59 and b63 of Orgel et al. (2020), northwest of
Caloris? (5 Nov 2026, Michael's observation.) Crater-level hosting (h09, crater_features.csv) with the controls of the
paper's tests; hollows and fields from host_assignments.csv. Writes province_basins.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
B = pd.read_csv(f"{ROOT}/peak_wall_test/orgel_basins.csv"); B["lat"] = B.y_coord; B["lon"] = B.x_coord % 360
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
h = pd.read_csv(f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv", low_memory=False)
h = h[~h.is_crater_trace.astype(bool)]
def rho(lat, lon, b):
    u = xyz(np.atleast_1d(lat), np.atleast_1d(lon) % 360); v = xyz([b.lat], [b.lon])[0]
    return np.arccos(np.clip(u @ v, -1, 1)) * R_M / (b.Diam_km / 2)
R = {}
sel = {n: B[B.Basin_Name == n].iloc[0] for n in ["b59", "b63", "b60", "b79", "Caloris"]}
R["basins"] = {n: dict(D=round(float(b.Diam_km), 1), lat=round(float(b.lat), 2), lon=round(float(b.lon), 2), strat=str(b.Stratigrap), comment=str(b.Comment), confidence=str(b.Confidence)) for n, b in sel.items()}
prov_f = (f.lat.between(40, 70)) & (f.lon.between(90, 150)); prov_h = (h.lat.between(40, 70)) & ((h.lon % 360).between(90, 150))
for n in ["b59", "b63", "b60", "b79"]:
    f[f"r_{n}"] = rho(f.lat.values, f.lon.values, sel[n]); h[f"r_{n}"] = rho(h.lat.values, h.lon.values, sel[n])
f["in59"] = (f.r_b59 < 1).astype(int); f["in63"] = (f.r_b63 < 1).astype(int)
f["inB"] = ((f.r_b59 < 1) | (f.r_b63 < 1)).astype(int); f["ring"] = (((f.r_b59 >= 1) & (f.r_b59 < 1.5)) | ((f.r_b63 >= 1) & (f.r_b63 < 1.5))).astype(int) * (1 - f.inB)
f["prov"] = prov_f.astype(int); f["prov_out"] = (prov_f & (f.inB == 0)).astype(int)
h["inB"] = ((h.r_b59 < 1) | (h.r_b63 < 1)); h["prov"] = prov_h
# where are the province's hollows and fields?
hp = h[h.prov]
R["province_hollows"] = {"n": int(len(hp)), "fields": int(hp.field_id.nunique()),
    "in_b59": int((hp.r_b59 < 1).sum()), "in_b63": int((hp.r_b63 < 1).sum()), "in_either": int(hp.inB.sum()),
    "fields_in_either": int(hp[hp.inB].field_id.nunique()),
    "in_b59_ring_1_1.5": int(((hp.r_b59 >= 1) & (hp.r_b59 < 1.5)).sum()), "in_b60": int((hp.r_b60 < 1).sum())}
prov_area_frac = None
# area of the province inside either basin (Monte Carlo on the sphere)
rng = np.random.default_rng(1); la = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(40)), np.sin(np.radians(70)), 400000))); lo = rng.uniform(90, 150, 400000)
inside = (rho(la, lo, sel["b59"]) < 1) | (rho(la, lo, sel["b63"]) < 1); R["province_area_fraction_in_either_basin"] = round(float(inside.mean()), 3)
# hollows inside the basins but outside the province box
R["hollows_in_either_basin_total"] = int(h.inB.sum()); R["hollows_in_either_basin_outside_box"] = int((h.inB & ~h.prov).sum())
# crater-level hosting at 30-75N, outside nothing excluded; controls as in the paper's tests
d = f[f.lat.between(30, 75)].copy()
def OR(m, t):
    b, ci = m.params[t], m.conf_int().loc[t]; return [round(float(np.exp(b)), 2), round(float(np.exp(ci[0])), 2), round(float(np.exp(ci[1])), 2), float(m.pvalues[t])]
ctl = "logD + abslat + C1 + nsp + dT"
for lab, form, terms in [("province", f"h09 ~ prov + {ctl}", ["prov"]),
                         ("basins_and_rest_of_province", f"h09 ~ inB + prov_out + {ctl}", ["inB", "prov_out"]),
                         ("each_basin", f"h09 ~ in59 + in63 + prov_out + {ctl}", ["in59", "in63", "prov_out"]),
                         ("basins_ring_rest", f"h09 ~ inB + ring + prov_out + {ctl}", ["inB", "ring", "prov_out"])]:
    m = smf.logit(form, d).fit(disp=0); R[f"hosting_30_75N_{lab}"] = {t: OR(m, t) for t in terms}; R[f"hosting_30_75N_{lab}"]["aic"] = round(float(m.aic), 1)
R["counts_30_75N"] = {k: [int(len(g)), int(g.h09.sum())] for k, g in [("in_b59", d[d.in59 == 1]), ("in_b63", d[d.in63 == 1]), ("province_outside_basins", d[d.prov_out == 1]), ("elsewhere", d[(d.inB == 0) & (d.prov == 0)])]}
# do other old basins (Orgel, not Caloris, D >= 300 km) host more? craters >= 5 km anywhere, same controls
old = B[(B.Basin_Name != "Caloris") & (B.Diam_km >= 300)]
inside_any = np.zeros(len(f), bool)
for _, b in old.iterrows(): inside_any |= rho(f.lat.values, f.lon.values, b) < 1
f["in_old"] = inside_any.astype(int); f["in_old_other"] = (inside_any & (f.inB == 0)).astype(int)
m = smf.logit(f"h09 ~ inB + in_old_other + {ctl}", f).fit(disp=0)
R["hosting_all_lat"] = {t: OR(m, t) for t in ["inB", "in_old_other"]}; R["n_old_basins"] = int(len(old))
json.dump(R, open(f"{HERE}/province_basins.json", "w"), indent=1); print(json.dumps(R, indent=1))
