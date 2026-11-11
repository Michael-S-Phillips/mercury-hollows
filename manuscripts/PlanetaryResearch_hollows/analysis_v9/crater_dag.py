#!/usr/bin/env python3
"""Containment graph of Mercury craters and the depth of origin of material in hollow hosts.

Adapts the superposition graph of mars-crater-graph (Phillips, Zenodo 10.5281/zenodo.23675713) to Mercury. A crater whose
center lies inside a feature at least twice its diameter is younger than that feature (the larger impact would have erased
it). Edges point older -> younger; with containment only, each crater has one immediate parent (the smallest feature that
contains it), so the graph is a forest.

Depth of origin. An impact inside an earlier crater excavates and uplifts material that the earlier crater already brought
up. For each crater the script estimates the deepest pre-impact level it exposes, z_max = z_base + z_own, where
  z_own  = max(excavation depth 0.1 D_tc, stratigraphic uplift 0.086 D^1.03 for complex craters) (Croft 1985;
           Cintala and Grieve 1998), with D_tc = D_sc^0.15 D^0.85, D_sc = 12 km (D_tc = D/1.2 for simple craters);
  z_base = depth of origin of the parent's surface where the child formed: z_base(parent) + SU(parent) on the parent's
           peak (rho < 0.2), z_base(parent) + h_exc(parent) on its floor (0.2-0.6), z_base(parent) on its wall (>= 0.6),
           where rho is the child's center distance in parent radii. A child on a filled floor must first cut through
           the fill (t_fill from the depth deficit against fresh craters, Section S7); a child whose excavation does not
           exceed t_fill exposes only fill (fill_only).
These are first-order scalings; the tests below ask only whether the extra depth from nesting predicts hosting.

Inputs (relative to the hollows-emissivity root, passed as argv[1]): Herrick et al. (2018) craters with host flags
(manuscripts/insolation_thermal/data/test2_craters_h18.csv), Orgel et al. (2020) basins (peak_wall_test/orgel_basins.csv),
DEM morphometry (peak_wall_test/morph_all.csv), hollow host assignments (manuscripts/insolation_thermal/data/host_assignments.csv), completeness
grid (code/rerun/out_v2/completeness_10ppd.npz), plains mask, thermal grid.
Writes crater_dag_nodes.csv and crater_dag_record.json next to this script.
Usage: python3 crater_dag.py ~/mnt/hollows-emissivity

v8 (30 Oct 2026): one Caloris node. v6 kept both the Herrick et al. (2018) crater 200153 (Caloris, D = 1595 km,
31.93 N 164.10 E) and the Orgel et al. (2020) basin "Caloris" (D = 1556 km, 29.91 N 161.43 E): their centers are 130 km
apart, more than the 0.15-radius duplicate tolerance (117 km). v8 keeps one record with the Orgel center and diameter.
By default (CRATER_DAG_CALORIS=catalog) the kept record is the catalog crater 200153 (kind crater, its hollow-host flags,
which host_assignments.csv keys to the catalog id) with the Orgel geometry, and the Orgel row is dropped, as the duplicate
step does for every other cataloged basin. CRATER_DAG_CALORIS=orgel instead keeps the Orgel record (kind basin, no host
flags) and drops 200153. Either way the merged node takes as children every feature inside either footprint (Orgel or
catalog circle), so the children of both v6 nodes attach to it; rho is measured from the Orgel center and radius, so
children inside the catalog rim but beyond the Orgel rim have rho > 1 (zone wall). The Caloris fill stays 3.0 km.
Other options (environment): CRATER_DAG_SU = terrestrial (default, 0.086 D^1.03), lunar (0.022 D^1.45, Cintala and
Grieve 1998 eq. 12), or double (2 x 0.086 D^1.03); CRATER_DAG_OUT = output folder (default: this script's folder).

v9 (1 Nov 2026): basin relation. v8 applied the central-peak uplift laws to every complex crater, including peak-ring
craters and basins, far beyond the sizes they were fitted to (up to about 136 km for the lunar relation): 167 km
(terrestrial) or 935 km (lunar) for Caloris, and up to 213 or 1,324 km for the largest basin. v9 uses, for features with
peak-ring geometry, the deepest origin of peak-ring and central material u_s = k_ring x 0.1 D_tc, with k_ring = 1.05 (the
midpoint of the range 0.8-1.3 used in provenance_mc_v2_2.py): lunar peak rings come from near the maximum excavation
depth (Baker and Head 2015); Schrodinger's peak ring from <= 20-26 km in models against 19-24 km excavation depth and
observed lithologies up to about 30 km (Kring et al. 2016); Chicxulub's from the mid-crust, about 8-10 km
(Urrutia-Fucugauchi et al. 2019). Peak-ring geometry = Orgel basins, Herrick 'pr' or 'pb', or any crater larger than
D_pr = 120 km (the midpoint of the 100-140 km onset used in the provenance model) that is not marked 'cp' or 'mp'.
Central-peak craters keep the uplift laws. Options: CRATER_DAG_BASIN = v22 (default) or v8 (the v8 rule);
CRATER_DAG_KRING (default 1.05); CRATER_DAG_DPR (default 120). The record also counts craters and hosts whose z_max
exceeds 17, 35 and 53 km (the crust-base range 35 +- 18 km of Padovan et al. 2015), where the exposed material would be
mantle.
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/mnt/hollows-emissivity")
OUT = os.environ.get("CRATER_DAG_OUT", os.path.dirname(os.path.abspath(__file__))); os.makedirs(OUT, exist_ok=True)
CAL_MODE = os.environ.get("CRATER_DAG_CALORIS", "catalog"); SU_REL = os.environ.get("CRATER_DAG_SU", "terrestrial")
assert CAL_MODE in ("catalog", "orgel") and SU_REL in ("terrestrial", "lunar", "double")
BASIN_REL = os.environ.get("CRATER_DAG_BASIN", "v22"); K_RING = float(os.environ.get("CRATER_DAG_KRING", "1.05"))
D_PR = float(os.environ.get("CRATER_DAG_DPR", "120")); assert BASIN_REL in ("v22", "v8")
G = f"{ROOT}/manuscripts/insolation_thermal"
R_M, DSC = 2439.4, 12.0
sys.path.insert(0, f"{G}/scripts"); import thermal

# ---------------------------------------------------------------- nodes
c = pd.read_csv(f"{G}/data/test2_craters_h18.csv")
c = pd.DataFrame(dict(id=c.unique_id.astype(str), kind="crater", lat=c.latitude, lon=c.lon360 % 360, D=c.diameter,
                      hosts=c.hosts, filled=c.filled.fillna("").str.strip().str.lower(),
                      preserved=c.preserved.fillna("").str.strip(), cent=c.cent_str.fillna("").str.strip()))
o = pd.read_csv(f"{ROOT}/peak_wall_test/orgel_basins.csv")
o = o[o.Confidence.isin(["Certain", "Probable"])]
b = pd.DataFrame(dict(id="orgel_" + o.Basin_Name.astype(str), kind="basin", lat=o.y_coord, lon=o.x_coord % 360,
                      D=o.Diam_km, hosts=np.nan, filled="", preserved="", cent=""))
def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
# v8: merge the two Caloris nodes (see docstring)
ih = c.index[c.id == "200153"]; io = b.index[b.id == "orgel_Caloris"]
assert len(ih) == 1 and len(io) == 1
ih, io = ih[0], io[0]
CAL_ORGEL = tuple(float(x) for x in b.loc[io, ["lat", "lon", "D"]])
CAL_CAT = tuple(float(x) for x in c.loc[ih, ["lat", "lon", "D"]])
if CAL_MODE == "catalog":
    c.loc[ih, ["lat", "lon", "D"]] = list(CAL_ORGEL); b = b.drop(index=io); CAL_ID = "200153"
else:
    c = c.drop(index=ih); CAL_ID = "orgel_Caloris"
uc, ub = unit(c.lat, c.lon), unit(b.lat, b.lon)
# drop Orgel basins that duplicate a cataloged crater (near-concentric, similar size)
dup = []
for i in range(len(b)):
    ang = np.arccos(np.clip(uc @ ub[i], -1, 1)) * R_M
    rb = b.D.iloc[i] / 2
    m = (ang < 0.15 * rb) & (c.D.values / b.D.iloc[i] > 0.7) & (c.D.values / b.D.iloc[i] < 1.43)
    dup.append(bool(m.any()))
b = b[~np.array(dup)]
n = pd.concat([c, b], ignore_index=True).sort_values("D", ascending=False).reset_index(drop=True)
u = unit(n.lat, n.lon); D = n.D.values; r = D / 2

# ---------------------------------------------------------------- containment (center inside a feature >= 2x larger)
tree = cKDTree(u)
parent = np.full(len(n), -1); rho = np.full(len(n), np.nan); ncont = np.zeros(len(n), int)
p_cal = int(np.where(n.id.values == CAL_ID)[0][0])
u_alt = unit(CAL_CAT[0], CAL_CAT[1])[0]; r_alt = CAL_CAT[2] / 2      # the catalog (Herrick) Caloris footprint
for p in np.where(D >= 10)[0]:                       # D sorted descending: later (smaller) parents overwrite
    rq = r[p] if p != p_cal else max(r[p], np.arccos(np.clip(u_alt @ u[p], -1, 1)) * R_M + r_alt)
    hits = np.asarray(tree.query_ball_point(u[p], 2 * np.sin(min(rq / R_M, np.pi) / 2)))
    hits = hits[D[hits] <= D[p] / 2]
    if hits.size == 0: continue
    d = np.arccos(np.clip(u[hits] @ u[p], -1, 1)) * R_M
    inside = d < r[p]
    if p == p_cal:                                   # union of the Orgel and catalog footprints
        inside |= np.arccos(np.clip(u[hits] @ u_alt, -1, 1)) * R_M < r_alt
    hits, d = hits[inside], d[inside]
    ncont[hits] += 1
    parent[hits] = p; rho[hits] = d / r[p]
n["parent"] = parent; n["rho"] = rho; n["n_containers"] = ncont
# ring artifacts: large, near-concentric nested circles
n["possible_ring"] = (parent >= 0) & (n.rho < 0.2) & (n.D >= 20) & (n.D / np.where(parent >= 0, D[np.maximum(parent, 0)], 1) >= 0.3)

# ---------------------------------------------------------------- depth model
Dtc = np.where(D > DSC, DSC ** 0.15 * D ** 0.85, D / 1.2)
h_exc = 0.1 * Dtc
SU = np.where(D > DSC, {"terrestrial": 0.086 * D ** 1.03, "lunar": 0.022 * D ** 1.45, "double": 2 * 0.086 * D ** 1.03}[SU_REL], 0.0)
cent_ = n.cent.fillna("").astype(str).str.strip().str.lower().values
peak_ring = (n.kind.values == "basin") | np.isin(cent_, ["pr", "pb"]) | (~np.isin(cent_, ["cp", "mp"]) & (D > D_PR))
if BASIN_REL == "v22":
    SU = np.where(peak_ring & (D > DSC), K_RING * h_exc, SU)
n["peak_ring_geom"] = peak_ring
z_own = np.maximum(h_exc, SU)
m = pd.read_csv(f"{ROOT}/peak_wall_test/morph_all.csv")[["uid", "depth"]]
m["id"] = m.uid.astype(str)
n = n.merge(m[["id", "depth"]], on="id", how="left")
d_fresh = 0.241 * D ** 0.486                          # km, 90th percentile depth of isolated craters (Section S7)
t_fill = np.clip(d_fresh - n.depth.values / 1000, 0, None)     # NaN where no DEM depth (craters < 10 km, basins)
t_fill = np.where(n.id == "200001", np.nan, t_fill)
cal = n.index[(n.D > 1500) & (n.D < 1700) & (abs(n.lat - 31) < 3)]
t_fill[cal] = 3.0                                     # Caloris interior plains 2.5-3.5 km (Ernst et al. 2015)
n["h_exc"], n["SU"], n["z_own"], n["t_fill"] = h_exc, SU, z_own, t_fill
zone = np.where(n.parent < 0, "none", np.where(n.rho < 0.2, "peak", np.where(n.rho < 0.6, "floor", "wall")))
# a simple-crater parent has no peak: treat its interior as floor
simple_parent = (n.parent >= 0) & (D[np.maximum(n.parent, 0)] <= DSC)
zone = np.where(simple_parent & (zone == "peak"), "floor", zone)
n["zone"] = zone
z_base = np.zeros(len(n)); z_max = z_own.copy(); fill_only = np.zeros(len(n), bool); fill_known = np.ones(len(n), bool)
for i in range(len(n)):                               # parents precede children (sorted by D)
    p = n.parent.iat[i]
    if p < 0: continue
    zb_p = z_base[p] if not fill_only[p] else 0.0
    zz = n.zone.iat[i]
    loc = SU[p] if zz == "peak" else (h_exc[p] if zz == "floor" else 0.0)
    tf = t_fill[p] if zz in ("peak", "floor") else 0.0
    if np.isnan(tf): tf = 0.0; fill_known[i] = False
    reach = z_own[i] - tf
    if reach <= 0:
        fill_only[i] = True; z_base[i] = 0.0; z_max[i] = 0.0
    else:
        z_base[i] = zb_p + loc; z_max[i] = z_base[i] + reach
n["z_base"], n["z_max"], n["fill_only"], n["fill_known"] = z_base, z_max, fill_only, fill_known
# chain depth (number of ancestors)
depth_chain = np.zeros(len(n), int)
for i in range(len(n)):
    p = n.parent.iat[i]
    if p >= 0: depth_chain[i] = depth_chain[p] + 1
n["chain"] = depth_chain

# ---------------------------------------------------------------- covariates
Cz = np.load(f"{ROOT}/code/rerun/out_v2/completeness_10ppd.npz"); C1 = Cz["C"][list(Cz["Dref"]).index(1.0)]
nsp = np.load(f"{G}/data/nsp_mask_10ppd_r40.npy")
rr = ((90 - n.lat.values) * 10).astype(int).clip(0, 1799); cc = (n.lon.values * 10).astype(int).clip(0, 3599)
n["C1"] = np.nan_to_num(C1[rr, cc]); n["nsp"] = nsp[rr, cc].astype(int)
z = np.load(f"{G}/data/thermal_grid.npz")
n["Tm"] = thermal.lookup(z["la"], z["lo"], z["G"], n.lat.values, n.lon.values)
n["dT"] = n.Tm - thermal.lookup(z["la"], z["lo"], z["G"], n.lat.values, np.full(len(n), 90.0))
h = pd.read_csv(f"{G}/data/host_assignments.csv", usecols=["host_uid"])
k = h.host_uid.dropna().astype(int).astype(str).value_counts()
n["n_hollows"] = n.id.map(k).fillna(0).astype(int)
n["abslat"] = n.lat.abs(); n["logD"] = np.log10(n.D)
n["parent_id"] = np.where(n.parent >= 0, n.id.values[np.maximum(n.parent, 0)], "")
n["parent_D"] = np.where(n.parent >= 0, D[np.maximum(n.parent, 0)], np.nan)
n["parent_kind"] = np.where(n.parent >= 0, n.kind.values[np.maximum(n.parent, 0)], "")
n["parent_fresh_like"] = np.where(n.parent >= 0, n.preserved.values[np.maximum(n.parent, 0)], "")
n.drop(columns=["depth"]).to_csv(f"{OUT}/crater_dag_nodes.csv", index=False)

# ---------------------------------------------------------------- tests (cataloged craters >= 5 km)
R = {}
ch = n[n.parent == p_cal]
R["settings"] = dict(caloris_mode=CAL_MODE, caloris_kept=CAL_ID, caloris_geometry_lat_lon_D=CAL_ORGEL,
                     caloris_catalog_geometry_lat_lon_D=CAL_CAT, uplift=SU_REL,
                     caloris_children=int(len(ch)), caloris_children_beyond_orgel_rim=int((ch.rho >= 1).sum()),
                     basin_relation=BASIN_REL, k_ring=K_RING, D_pr_km=D_PR, n_peak_ring_geom=int(n.peak_ring_geom.sum()),
                     caloris_SU_km=float(n.SU.iat[p_cal]), caloris_z_own_km=float(n.z_own.iat[p_cal]),
                     caloris_children_hosts=int(ch.hosts.sum()), caloris_t_fill=float(n.t_fill.iat[p_cal]),
                     caloris_kind=str(n.kind.iat[p_cal]), caloris_hosts=None if pd.isna(n.hosts.iat[p_cal]) else int(n.hosts.iat[p_cal]))
s = n[(n.kind == "crater") & (~n.possible_ring)].copy()
s["nested"] = (s.parent >= 0).astype(int)
s["ln_zown"] = np.log(s.z_own); s["ln_zmax"] = np.log(np.maximum(s.z_max, 0.05))
R["counts"] = dict(nodes=int(len(n)), craters=int((n.kind == "crater").sum()), basins_kept=int((n.kind == "basin").sum()),
                   nested=int(s.nested.sum()), by_zone=s.zone.value_counts().to_dict(),
                   chain=s.chain.value_counts().sort_index().to_dict(), possible_ring=int(n.possible_ring.sum()),
                   fill_only=int(s.fill_only.sum()), fill_unknown_nested=int(((~s.fill_known) & (s.nested == 1)).sum()),
                   hosts=int(s.hosts.sum()), hosts_nested=int(s[s.nested == 1].hosts.sum()))
def fit(formula, d, name):
    mm = smf.logit(formula, data=d).fit(disp=0)
    R[name] = dict(n=int(mm.nobs), hosts=int(d.hosts.sum()), aic=float(mm.aic),
                   params={k: [float(v), float(mm.bse[k]), float(mm.pvalues[k])] for k, v in mm.params.items()})
    return mm
ctrl = "logD + abslat + C1 + nsp + dT"
fit(f"hosts ~ {ctrl}", s, "base")
fit(f"hosts ~ {ctrl} + nested", s, "nested")
s["zone_c"] = pd.Categorical(np.where(s.fill_only, "fill_only", s.zone), categories=["none", "peak", "floor", "wall", "fill_only"])
fit(f"hosts ~ {ctrl} + C(zone_c)", s, "zone")
fit(f"hosts ~ {ctrl} + z_base", s, "z_base")
fit("hosts ~ ln_zown + abslat + C1 + nsp + dT", s, "ln_zown")
fit("hosts ~ ln_zmax + abslat + C1 + nsp + dT", s, "ln_zmax")
fit("hosts ~ ln_zown + ln_zmax + abslat + C1 + nsp + dT", s, "ln_zown_zmax")
# same tests where the fill of every nested crater's parent is known, and with completeness >= 0.2
sk = s[(s.nested == 0) | (s.fill_known)]
fit(f"hosts ~ {ctrl} + z_base", sk, "z_base_fillknown")
fit("hosts ~ ln_zown + abslat + C1 + nsp + dT", sk, "ln_zown_fillknown")
fit("hosts ~ ln_zmax + abslat + C1 + nsp + dT", sk, "ln_zmax_fillknown")
sc = s[s.C1 >= 0.2]
fit(f"hosts ~ {ctrl} + C(zone_c)", sc, "zone_C02")
fit("hosts ~ ln_zown + abslat + C1 + nsp + dT", sc, "ln_zown_C02")
fit("hosts ~ ln_zmax + abslat + C1 + nsp + dT", sc, "ln_zmax_C02")
# raw rates by zone and size bin (completeness >= 0.2)
rows = []
for lo_, hi_ in [(5, 10), (10, 20), (20, 40), (40, 100), (100, 2000)]:
    t = sc[(sc.D >= lo_) & (sc.D < hi_)]
    row = dict(D=f"{lo_}-{hi_}")
    for zz in ["none", "peak", "floor", "wall", "fill_only"]:
        tt = t[t.zone_c == zz]; row[zz] = f"{int(tt.hosts.sum())}/{len(tt)}"
    rows.append(row)
R["rates_C02"] = rows
# chain depth at fixed size
fit(f"hosts ~ {ctrl} + C(chain.clip(upper=3))", s.assign(chain=s.chain), "chain")
# hollow-weighted: where do hollows sit?
hh = s[s.n_hollows > 0]
R["hollow_hosts"] = dict(n=int(len(hh)), hollows=int(hh.n_hollows.sum()),
                         frac_nested=float(hh.nested.mean()), frac_nested_all=float(s.nested.mean()),
                         z_base_median_nested_hosts=float(hh[hh.nested == 1].z_base.median()),
                         z_base_median_nested_all=float(s[s.nested == 1].z_base.median()))
R["crust_base"] = {grp: {f"z_max_gt_{T}": int((d_.z_max > T).sum()) for T in (17, 35, 53)} | {"n": int(len(d_)), "z_max_max": float(d_.z_max.max())}
                   for grp, d_ in [("all_craters", s), ("hosts", s[s.hosts == 1])]}
R["crust_base"]["note"] = "z_max beyond a crust base of 35 +- 18 km (Padovan et al. 2015) would expose mantle"
json.dump(R, open(f"{OUT}/crater_dag_record.json", "w"), indent=1)
print(json.dumps(R, indent=1))
