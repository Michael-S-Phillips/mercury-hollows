"""Shared pieces for the insolation_v6 workshop analyses (2 Nov 2026).

The null intensity of hollow fields is the one in code/rerun/latitude/hotwarm/hotwarm.py (size-marginalized
coverage x completeness, built per hemisphere, Caloris < 1,200 km and optionally the 90-150E / 40-70N province
excluded). The functions below are copied from hotwarm.py, which cannot be imported (its paths are hard-coded).
Here expected counts are sums of the null intensity over grid cells, not draws of 400,000 points per hemisphere.
"""
import os
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))          # hollows-emissivity
DD = f"{ROOT}/code/rerun/out_v2/"                                             # same files (md5) as hollow-networks/data_v2
PR = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"
R_MERCURY = 2439.7
CALORIS = (31.9, 161.9); CALORIS_KM = 1200.0
PROV_LON = (90.0, 150.0); PROV_LAT = (40.0, 70.0)


def grid_cell_centers(ppd=10):
    step = 1.0 / ppd
    return 90 - (np.arange(180 * ppd) + 0.5) * step, (np.arange(360 * ppd) + 0.5) * step


def gc_dist_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, (lat1, lon1, lat2, lon2))
    d = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    return R_MERCURY * np.arccos(np.clip(d, -1, 1))


def hot_mask(lon):
    lon = np.mod(lon, 360.0)
    return (np.minimum(np.abs(lon), np.abs(lon - 360)) <= 45.0) | (np.abs(lon - 180) <= 45.0)


def load():
    cov = np.load(DD + "coverage_10ppd.npz"); comp = np.load(DD + "completeness_10ppd.npz")
    fields = pd.read_csv(DD + "fields.csv"); hollows = pd.read_csv(DD + "hollows_canonical.csv")
    hnt = hollows[~hollows.is_crater_trace.astype(bool)]
    dmax_nt = hnt.groupby("field_id").D_eq_km.max()
    fields = fields.merge(dmax_nt.rename("D_max_nt").reset_index(), on="field_id", how="left")
    fields["D_max_nt"] = fields["D_max_nt"].fillna(fields["D_max_km"].clip(upper=9.71))
    return cov, comp, fields, hollows


def build_null(hemi, size_source_fields, comp, cov, lat_c, lon_c, exclude_caloris=False, exclude_province=False):
    Dref = comp["Dref"]; C = comp["C"]; nimg = cov["nimg"]
    coslat = np.cos(np.radians(lat_c))
    LAT = np.broadcast_to(lat_c[:, None], C[0].shape); LON = np.broadcast_to(lon_c[None, :], C[0].shape)
    COSLAT = np.broadcast_to(coslat[:, None], C[0].shape)
    covered = nimg > 0
    hemi_mask = (LAT > 0) if hemi == '+' else (LAT < 0)
    edges = np.concatenate([[0], (Dref[:-1] + Dref[1:]) / 2, [np.inf]])
    w_k, _ = np.histogram(size_source_fields["D_max_nt"].values, bins=edges); w_k = w_k / w_k.sum()
    Cmix = np.tensordot(w_k, C, axes=(0, 0))
    keep = hemi_mask & covered
    if exclude_caloris:
        keep = keep & (gc_dist_km(CALORIS[0], CALORIS[1], LAT, LON) > CALORIS_KM)
    if exclude_province:
        prov = (LON >= PROV_LON[0]) & (LON <= PROV_LON[1]) & (LAT >= PROV_LAT[0]) & (LAT <= PROV_LAT[1])
        keep = keep & (~prov)
    raw = np.where(keep & (Cmix > 0), Cmix * COSLAT, 0.0)
    lam = gaussian_filter(raw, sigma=3.0, mode=["reflect", "wrap"])
    return np.where(keep, lam, 0.0), keep, w_k, Cmix


def apply_point_exclusions(df, exclude_caloris, exclude_province):
    m = np.ones(len(df), bool)
    if exclude_caloris:
        m &= gc_dist_km(CALORIS[0], CALORIS[1], df.lat.values, df.lon.values) > CALORIS_KM
    if exclude_province:
        m &= ~((df.lon.values >= PROV_LON[0]) & (df.lon.values <= PROV_LON[1]) &
               (df.lat.values >= PROV_LAT[0]) & (df.lat.values <= PROV_LAT[1]))
    return df[m].reset_index(drop=True)


def wls_cos2(ratio, exp, x):
    g = np.isfinite(ratio) & (exp > 0)
    X = np.column_stack([np.ones(g.sum()), x[g]]); W = np.diag(exp[g])
    return np.linalg.solve(X.T @ W @ X, X.T @ W @ ratio[g])
