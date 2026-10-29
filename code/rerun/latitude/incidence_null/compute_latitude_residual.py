import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
"""
Incidence-controlled completeness null, STEP 3.
Rebuild the size-marginalized coverage x completeness null intensity
(net_pipeline.py logic, northern hemisphere) from the incidence-restricted
(60-75 deg) grids, and recompute the observed/null field-fraction profile per
10-deg latitude band using ONLY fields whose canonical hollow was detected in a
60-75 deg frame. Compare with the v2 panel-d profile and with the fixed-
resolution field densities. American English; every number with n and a
bootstrap CI.
"""
import os, json
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter

R = 2439.7
OUT = _REPO + "/code/rerun/latitude/incidence_null/"
Rd = OUT + "out_inc60_75/"
SEED = 20260911
rng = np.random.default_rng(SEED)

# reference profiles
V2_PANELD = [0.65, 0.70, 1.02, 0.98, 2.10, 2.59, 1.72, 0.83]      # obs/null, 0-80N
FIXEDRES_DENS = [42, 56, 73, 48, 68, 76, 47, 29]                   # per 1e6 km^2, <=40 m class
BANDS = [(b, b+10) for b in range(0, 80, 10)]

# ---------------------------------------------------------------- load restricted grids
fields = pd.read_csv(Rd+"fields.csv")
hollows = pd.read_csv(Rd+"hollows_canonical.csv")
det = pd.read_csv(Rd+"detections_qc.csv")
cov = np.load(Rd+"coverage_10ppd.npz"); nimg = cov["nimg"]
comp = np.load(Rd+"completeness_10ppd.npz"); C = comp["C"]; Dref = comp["Dref"]

ppd = 10; nr, nc = 180*ppd, 360*ppd
lat_c = 90 - (np.arange(nr)+0.5)/ppd
lon_c = (np.arange(nc)+0.5)/ppd
coslat = np.cos(np.radians(lat_c))
LAT = np.broadcast_to(lat_c[:, None], C[0].shape)
COSLAT = np.broadcast_to(coslat[:, None], C[0].shape)
covered = (nimg > 0)
deg_per_cell = 0.1

# ---------------------------------------------------------------- observed field set
# fields whose canonical hollow was detected in a 60-75 deg frame (any detection of a
# non-trace member hollow with inc in [60,75])
det_in = det[det.inc.between(60, 75)].hollow_id.unique()
hollows["in6075"] = hollows.hollow_id.isin(det_in)
nt = hollows[~hollows.is_crater_trace.astype(bool)]
qual_fields = nt[nt.in6075].field_id.unique()
fn_all = fields[fields.lat > 0].reset_index(drop=True)
fn = fn_all[fn_all.field_id.isin(qual_fields)].reset_index(drop=True)
N_obs = len(fn)

# ---------------------------------------------------------------- size-marginalized null
# non-trace field D_max distribution (size-marginalized detectability), matching net_pipeline
hnt = hollows[(hollows.lat > 0) & (~hollows.is_crater_trace.astype(bool))]
dmax_nt = hnt.groupby("field_id").D_eq_km.max()
fn = fn.merge(dmax_nt.rename("D_max_nt").reset_index(), on="field_id", how="left")
fn["D_max_nt"] = fn["D_max_nt"].fillna(fn["D_max_km"].clip(upper=9.71))
edges = np.concatenate([[0], (Dref[:-1]+Dref[1:])/2, [np.inf]])
w_k, _ = np.histogram(fn.D_max_nt.values, bins=edges); w_k = w_k/w_k.sum()
Cmix = np.tensordot(w_k, C, axes=(0, 0))
raw_mix = np.where((LAT > 0) & covered & (Cmix > 0), Cmix*COSLAT, 0.0)
lam_mix = gaussian_filter(raw_mix, sigma=3.0, mode=["reflect", "wrap"])
lam_mix = np.where((LAT > 0) & covered, lam_mix, 0.0)

def make_tables(lg):
    m = lg > 0; rr, cc = np.where(m); w = lg[rr, cc].astype(np.float64)
    return w/w.sum(), lat_c[rr], lon_c[cc]

def sample_null(n, p, cl, co):
    idx = rng.choice(len(p), size=n, p=p)
    la = cl[idx] + (rng.random(n)-0.5)*deg_per_cell
    lo = (co[idx] + (rng.random(n)-0.5)*deg_per_cell) % 360.0
    return la, lo

p, cl, co = make_tables(lam_mix)
la_pool, lo_pool = sample_null(200000, p, cl, co)

# ---------------------------------------------------------------- per-band residual + bootstrap
obs_lat = fn.lat.values
null_frac = np.array([((la_pool >= a) & (la_pool < b)).mean() for a, b in BANDS])
obs_frac = np.array([((obs_lat >= a) & (obs_lat < b)).mean() for a, b in BANDS])
ratio = obs_frac/np.where(null_frac == 0, np.nan, null_frac)
n_band = np.array([int(((obs_lat >= a) & (obs_lat < b)).sum()) for a, b in BANDS])

# bootstrap observed fields -> ratio CI (null held fixed)
NB = 2000
boot = np.empty((NB, len(BANDS)))
for t in range(NB):
    bl = obs_lat[rng.integers(0, N_obs, N_obs)]
    of = np.array([((bl >= a) & (bl < b)).mean() for a, b in BANDS])
    boot[t] = of/np.where(null_frac == 0, np.nan, null_frac)
ci_lo = np.nanpercentile(boot, 2.5, axis=0)
ci_hi = np.nanpercentile(boot, 97.5, axis=0)

# northern searched area per band (restricted coverage) and observed field DENSITY per 1e6 km2
lat_top = np.radians(np.minimum(lat_c+0.05, 90)); lat_bot = np.radians(np.maximum(lat_c-0.05, -90))
row_area = R*R*np.radians(0.1)*(np.sin(lat_top)-np.sin(lat_bot))
area_band = []
for a, b in BANDS:
    rows_b = np.where((lat_c >= a) & (lat_c < b))[0]
    area_band.append(float((covered[rows_b].sum(axis=1)*row_area[rows_b]).sum()))
area_band = np.array(area_band)
dens = n_band/(area_band/1e6)   # qualifying fields per 1e6 km^2 of restricted searched area

rows = []
for i, (a, b) in enumerate(BANDS):
    rows.append(dict(band=f"{a}-{b}", n_obs_fields=int(n_band[i]),
                     obs_frac=float(obs_frac[i]), null_frac=float(null_frac[i]),
                     ratio=float(ratio[i]), ratio_CI95=[float(ci_lo[i]), float(ci_hi[i])],
                     searched_area_km2=float(area_band[i]),
                     obs_field_density_per_1e6km2=float(dens[i]),
                     v2_paneld=V2_PANELD[i], fixedres_density=FIXEDRES_DENS[i]))
prof = pd.DataFrame(rows)
prof.to_csv(OUT+"latitude_residual_inc60_75.csv", index=False)

# correlations vs reference profiles (Spearman)
from scipy.stats import spearmanr
sp_v2 = spearmanr(ratio, V2_PANELD)
sp_fx = spearmanr(dens, FIXEDRES_DENS)

record = dict(
    N_obs_qualifying_fields=int(N_obs), N_north_fields_total=int(len(fn_all)),
    observed_set="north fields with >=1 non-trace member hollow detected in a 60-75 deg frame",
    null_definition="lambda_0 = smoothed[ cos(lat) * sum_k w_k C(Dref_k) ] on northern covered cells, "
                    "restricted (60-75 deg) coverage x completeness grids; w_k = qualifying-field D_max distribution",
    band_ratio=[float(x) for x in ratio],
    band_ratio_CI95=[[float(ci_lo[i]), float(ci_hi[i])] for i in range(len(BANDS))],
    band_n=[int(x) for x in n_band],
    v2_paneld=V2_PANELD, fixedres_density=FIXEDRES_DENS,
    obs_field_density_per_1e6km2=[float(x) for x in dens],
    spearman_ratio_vs_v2=dict(rho=float(sp_v2.correlation), p=float(sp_v2.pvalue)),
    spearman_density_vs_fixedres=dict(rho=float(sp_fx.correlation), p=float(sp_fx.pvalue)),
    nboot=NB, seed=SEED)
json.dump(record, open(OUT+"latitude_residual_inc60_75.json", "w"), indent=1)
np.savez_compressed(OUT+"out_inc60_75/net_null_intensity_inc60_75.npz",
                    lam=lam_mix.astype(np.float32), lat_c=lat_c, lon_c=lon_c, ppd=10,
                    size_weights=w_k, Dref=Dref)
print(prof.to_string(index=False))
print("Spearman ratio vs v2 panel-d: rho=%.3f p=%.3f" % (sp_v2.correlation, sp_v2.pvalue))
print("Spearman density vs fixed-res: rho=%.3f p=%.3f" % (sp_fx.correlation, sp_fx.pvalue))
