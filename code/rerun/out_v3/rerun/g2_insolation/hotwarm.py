# g2 copy of code/rerun/latitude/hotwarm/hotwarm.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""
Hot-pole vs warm-pole longitude preference of Mercury hollow fields, both
hemispheres, corrected (v2) positions.

Mercury's 3:2 spin-orbit resonance makes the hot-pole longitudes (0 and 180 E)
peak-hotter than the warm-pole longitudes (90 and 270 E) at low latitude. If
insolation drives hollow formation/exposure, hollow fields should be
over-represented at hot-pole longitudes relative to a coverage x completeness
null. Yesterday's northern-only result (net_pipeline size-marginalized null):
0-15N hot/warm 1.9, 15-30N 1.67, 45-60N 0.43 (the 90-150E province).

This script:
  (1) Builds the size-marginalized coverage null SEPARATELY for the southern and
      northern hemispheres, same recipe as net_pipeline.py
      (lam = smoothed[ cos(lat) * sum_k w_k C(Dref_k) ] on covered cells,
       w_k = observed non-trace field D_max distribution, gaussian sigma=3,
       mode reflect(lat)/wrap(lon)). South size-weight sensitivity vs the
      northern weights is reported.
  (2) Observed/null field fraction for hot-pole (0/180E +-45) vs warm-pole
      (90/270E +-45) longitudes in 15-deg latitude bands, both hemispheres, for
      (a) all fields, (b) Caloris < 1200 km excluded, (c) additionally the
      90-150E / 40-70N province excluded. Bootstrap CIs (resample fields) and a
      permutation p (draw fields from the null, hot/warm ratio distribution).
  (3) Continuous fit: observed/null vs cos(2*lon) within |lat|<30, both
      hemispheres, province + Caloris excluded.
  (4) The same hot/warm test on individual canonical hollows (not fields).

Mercury R = 2439.7 km. Grids 10 ppd, north at row 0, lon 0-360 from col 0.
American English throughout.
"""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_INSOLDATA = f"{G2_ROOT}/manuscripts/insolation_thermal/data"
# --- end g2 block ---
import os, json
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

R_MERCURY = 2439.7
DD = f"{G2_ROOT}/code/rerun/out_v2/"   # same files (md5) as hollow-networks/data_v2
OUT = G2_OUT + "/"
os.makedirs(OUT, exist_ok=True)
CALORIS = (31.9, 161.9)          # lat, lon E
CALORIS_KM = 1200.0
PROV_LON = (90.0, 150.0)         # 90-150 E
PROV_LAT = (40.0, 70.0)          # 40-70 N
SEED = 20260911
DEG_PER_CELL = 0.1


# ---------------------------------------------------------------- geometry
def grid_cell_centers(ppd=10):
    nr, nc = 180 * ppd, 360 * ppd
    step = 1.0 / ppd
    lat_c = 90 - (np.arange(nr) + 0.5) * step
    lon_c = (np.arange(nc) + 0.5) * step
    return lat_c, lon_c


def gc_dist_km(lat1, lon1, lat2, lon2):
    """Great-circle distance (km) from scalar point (lat1,lon1) to arrays."""
    la1, lo1 = np.radians(lat1), np.radians(lon1)
    la2, lo2 = np.radians(lat2), np.radians(lon2)
    d = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    return R_MERCURY * np.arccos(np.clip(d, -1, 1))


def hot_mask(lon):
    """Hot-pole = within +-45 deg of 0 E or 180 E."""
    lon = np.mod(lon, 360.0)
    d0 = np.minimum(np.abs(lon - 0), np.abs(lon - 360))
    d180 = np.abs(lon - 180)
    return (d0 <= 45.0) | (d180 <= 45.0)


# ---------------------------------------------------------------- null builder
def build_null(hemi, size_source_fields, comp, cov, lat_c, lon_c,
               exclude_caloris=False, exclude_province=False):
    """Size-marginalized coverage x completeness null intensity for one hemisphere."""
    Dref = comp["Dref"]; C = comp["C"]; nimg = cov["nimg"]
    coslat = np.cos(np.radians(lat_c))
    LAT = np.broadcast_to(lat_c[:, None], C[0].shape)
    LON = np.broadcast_to(lon_c[None, :], C[0].shape)
    COSLAT = np.broadcast_to(coslat[:, None], C[0].shape)
    covered = (nimg > 0)
    hemi_mask = (LAT > 0) if hemi == '+' else (LAT < 0)

    edges = np.concatenate([[0], (Dref[:-1] + Dref[1:]) / 2, [np.inf]])
    w_k, _ = np.histogram(size_source_fields["D_max_nt"].values, bins=edges)
    w_k = w_k / w_k.sum()
    Cmix = np.tensordot(w_k, C, axes=(0, 0))

    keep = hemi_mask & covered
    if exclude_caloris:
        dcal = gc_dist_km(CALORIS[0], CALORIS[1], LAT, LON)
        keep = keep & (dcal > CALORIS_KM)
    if exclude_province:
        prov = ((LON >= PROV_LON[0]) & (LON <= PROV_LON[1]) &
                (LAT >= PROV_LAT[0]) & (LAT <= PROV_LAT[1]))
        keep = keep & (~prov)

    raw = np.where(keep & (Cmix > 0), Cmix * COSLAT, 0.0)
    lam = gaussian_filter(raw, sigma=3.0, mode=["reflect", "wrap"])
    lam = np.where(keep, lam, 0.0)
    return lam, keep, w_k


def make_pool(lam, lat_c, lon_c, n, rng):
    m = lam > 0
    rr, cc = np.where(m)
    w = lam[rr, cc].astype(np.float64)
    p = w / w.sum()
    idx = rng.choice(len(p), size=n, p=p)
    la = lat_c[rr][idx] + (rng.random(n) - 0.5) * DEG_PER_CELL
    lo = (lon_c[cc][idx] + (rng.random(n) - 0.5) * DEG_PER_CELL) % 360.0
    return la, lo


# ---------------------------------------------------------------- hot/warm test
def apply_point_exclusions(df, exclude_caloris, exclude_province):
    m = np.ones(len(df), bool)
    if exclude_caloris:
        dcal = gc_dist_km(CALORIS[0], CALORIS[1], df.lat.values, df.lon.values)
        m &= (dcal > CALORIS_KM)
    if exclude_province:
        prov = ((df.lon.values >= PROV_LON[0]) & (df.lon.values <= PROV_LON[1]) &
                (df.lat.values >= PROV_LAT[0]) & (df.lat.values <= PROV_LAT[1]))
        m &= (~prov)
    return df[m].reset_index(drop=True)


def hotwarm_by_band(obs_df, pool_lat, pool_lon, hemi, rng, nboot=2000, nperm=2000):
    """Per 15-deg |lat| band: hot/warm over-representation, bootstrap CI, perm p."""
    obs_lat = obs_df.lat.values; obs_lon = obs_df.lon.values
    obs_hot = hot_mask(obs_lon)
    pool_hot = hot_mask(pool_lon)
    rows = []
    for b in range(0, 90, 15):
        if hemi == '+':
            m_obs = (obs_lat >= b) & (obs_lat < b + 15)
            m_pool = (pool_lat >= b) & (pool_lat < b + 15)
        else:
            m_obs = (obs_lat <= -b) & (obs_lat > -(b + 15))
            m_pool = (pool_lat <= -b) & (pool_lat > -(b + 15))
        n = int(m_obs.sum()); npool = int(m_pool.sum())
        if n == 0 or npool < 50:
            rows.append(dict(band=f"{b}-{b+15}", n=n, npool=npool,
                             null_hot_frac=None, over_hot=None, over_warm=None,
                             hot_over_warm=None, hw_ci=None, perm_p=None))
            continue
        f_null_hot = pool_hot[m_pool].mean()
        f_null_warm = 1.0 - f_null_hot
        oh = obs_hot[m_obs]
        f_obs_hot = oh.mean(); f_obs_warm = 1.0 - f_obs_hot
        over_hot = f_obs_hot / f_null_hot if f_null_hot > 0 else np.nan
        over_warm = f_obs_warm / f_null_warm if f_null_warm > 0 else np.nan
        hw = over_hot / over_warm if over_warm > 0 else np.nan
        # bootstrap: resample observed points in band
        boot = np.empty(nboot)
        for i in range(nboot):
            s = oh[rng.integers(0, n, n)]
            fh = s.mean(); fw = 1 - fh
            boot[i] = (fh / f_null_hot) / (fw / f_null_warm) if (fw > 0 and f_null_hot > 0) else np.nan
        hw_ci = np.nanpercentile(boot, [2.5, 97.5]).tolist()
        # permutation: draw n points from the null in this band, hot/warm ratio dist
        pool_hot_band = pool_hot[m_pool]
        perm = np.empty(nperm)
        for i in range(nperm):
            s = pool_hot_band[rng.integers(0, npool, n)]
            fh = s.mean(); fw = 1 - fh
            perm[i] = (fh / f_null_hot) / (fw / f_null_warm) if (fw > 0 and f_null_hot > 0) else np.nan
        perm = perm[np.isfinite(perm)]
        dev_obs = abs(np.log(hw)) if hw > 0 else np.inf
        perm_p = (1 + np.sum(np.abs(np.log(perm)) >= dev_obs)) / (1 + len(perm))
        rows.append(dict(band=f"{b}-{b+15}", n=n, npool=npool,
                         null_hot_frac=round(float(f_null_hot), 4),
                         over_hot=round(float(over_hot), 3),
                         over_warm=round(float(over_warm), 3),
                         hot_over_warm=round(float(hw), 3),
                         hw_ci=[round(x, 3) for x in hw_ci],
                         perm_p=round(float(perm_p), 4)))
    return rows


# ---------------------------------------------------------------- cos(2 lon) fit
def cos2lon_fit(obs_df, pool_lat, pool_lon, rng, nbins=12, latmax=30.0, nboot=2000):
    """Regress observed/null field-fraction on cos(2*lon) within |lat|<latmax."""
    obs = obs_df[np.abs(obs_df.lat.values) < latmax]
    m_pool = np.abs(pool_lat) < latmax
    plon = pool_lon[m_pool]
    edges = np.linspace(0, 360, nbins + 1)
    ctr = (edges[:-1] + edges[1:]) / 2
    obs_c, _ = np.histogram(np.mod(obs.lon.values, 360), bins=edges)
    null_c, _ = np.histogram(np.mod(plon, 360), bins=edges)
    null_frac = null_c / null_c.sum()
    exp = null_frac * obs_c.sum()
    ratio = obs_c / np.where(exp > 0, exp, np.nan)
    x = np.cos(np.radians(2 * ctr))
    good = np.isfinite(ratio) & (exp > 0)
    # weighted least squares, weight by expected counts
    w = exp[good]
    X = np.column_stack([np.ones(good.sum()), x[good]])
    W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ ratio[good])
    # bootstrap slope by resampling observed points
    n = len(obs); lon_o = np.mod(obs.lon.values, 360)
    slopes = np.empty(nboot)
    for i in range(nboot):
        idx = rng.integers(0, n, n)
        oc, _ = np.histogram(lon_o[idx], bins=edges)
        r = oc / np.where(exp > 0, exp, np.nan)
        gg = np.isfinite(r) & (exp > 0)
        Xi = np.column_stack([np.ones(gg.sum()), x[gg]])
        Wi = np.diag(exp[gg])
        try:
            bi = np.linalg.solve(Xi.T @ Wi @ Xi, Xi.T @ Wi @ r[gg])
            slopes[i] = bi[1]
        except Exception:
            slopes[i] = np.nan
    ci = np.nanpercentile(slopes, [2.5, 97.5]).tolist()
    p_boot = float(2 * min((slopes <= 0).mean(), (slopes >= 0).mean()))
    return dict(n=int(n), nbins=nbins, latmax=latmax,
                intercept=float(beta[0]), slope=float(beta[1]),
                slope_ci=[float(c) for c in ci], slope_boot_p=p_boot,
                lon_bin_centers=ctr.tolist(), obs_over_null_ratio=ratio.tolist(),
                cos2lon=x.tolist())


# ---------------------------------------------------------------- main
def main():
    lat_c, lon_c = grid_cell_centers(10)
    cov = np.load(DD + "coverage_10ppd.npz")
    comp = np.load(DD + "completeness_10ppd.npz")
    fields = pd.read_csv(IN + "/fields.csv")
    hollows = pd.read_csv(IN + "/hollows_canonical.csv")

    # non-trace field D_max (as in net_pipeline final null)
    hnt = hollows[~hollows.is_crater_trace.astype(bool)]
    dmax_nt = hnt.groupby("field_id").D_eq_km.max()
    fields = fields.merge(dmax_nt.rename("D_max_nt").reset_index(), on="field_id", how="left")
    fields["D_max_nt"] = fields["D_max_nt"].fillna(fields["D_max_km"].clip(upper=9.71))

    fn = fields[fields.lat > 0].reset_index(drop=True)
    fs = fields[fields.lat < 0].reset_index(drop=True)

    rng = np.random.default_rng(SEED)
    NPOOL = 400000
    record = {"description": "Hot-pole vs warm-pole longitude preference, both hemispheres, v2 corrected positions.",
              "definitions": {
                  "hot_pole_lon": "within +-45 deg of 0 or 180 E",
                  "warm_pole_lon": "within +-45 deg of 90 or 270 E",
                  "over_representation": "observed fraction in sector / null fraction in sector, within a latitude band",
                  "hot_over_warm": "over_rep(hot) / over_rep(warm)",
                  "null": "size-marginalized coverage x completeness (net_pipeline recipe), built per hemisphere",
                  "caloris_exclusion_km": CALORIS_KM,
                  "province": "90-150 E, 40-70 N"},
              "n_fields": {"north": int(len(fn)), "south": int(len(fs))},
              "variants": {}}

    variants = [("all", False, False),
                ("caloris_excl", True, False),
                ("caloris_and_province_excl", True, True)]

    for vname, exc_cal, exc_prov in variants:
        vrec = {}
        for hemi, fh, name in [('+', fn, "north"), ('-', fs, "south")]:
            # southern null uses SOUTHERN size weights; northern uses NORTHERN
            src = fh
            lam, keep, w_k = build_null(hemi, src, comp, cov, lat_c, lon_c,
                                        exclude_caloris=exc_cal, exclude_province=exc_prov)
            pool_lat, pool_lon = make_pool(lam, lat_c, lon_c, NPOOL, rng)
            obs = apply_point_exclusions(fh, exc_cal, exc_prov)
            rows = hotwarm_by_band(obs, pool_lat, pool_lon, hemi, rng)
            vrec[name] = {"n_after_excl": int(len(obs)),
                          "size_weights": [round(float(x), 4) for x in w_k],
                          "bands": rows}
        record["variants"][vname] = vrec

    # ---- south size-weight sensitivity: rebuild south null with NORTH weights
    lam_s_nw, _, w_north = build_null('-', fn, comp, cov, lat_c, lon_c)
    pl_lat, pl_lon = make_pool(lam_s_nw, lat_c, lon_c, NPOOL, rng)
    obs_s = fs
    rows_nw = hotwarm_by_band(obs_s, pl_lat, pl_lon, '-', rng)
    record["south_null_sensitivity_northweights"] = {
        "note": "Southern null intensity rebuilt using the NORTHERN non-trace field D_max weights instead of southern; south is far less complete so the size weights matter.",
        "north_weights": [round(float(x), 4) for x in w_north],
        "bands": rows_nw}

    # ---- (3) cos(2 lon) fit, |lat|<30, both hemispheres, province+Caloris excluded
    lam_n3, _, _ = build_null('+', fn, comp, cov, lat_c, lon_c, True, True)
    lam_s3, _, _ = build_null('-', fs, comp, cov, lat_c, lon_c, True, True)
    pl_n_lat, pl_n_lon = make_pool(lam_n3, lat_c, lon_c, NPOOL, rng)
    pl_s_lat, pl_s_lon = make_pool(lam_s3, lat_c, lon_c, NPOOL, rng)
    pool_lat3 = np.concatenate([pl_n_lat, pl_s_lat])
    pool_lon3 = np.concatenate([pl_n_lon, pl_s_lon])
    both3 = apply_point_exclusions(fields[fields.lat.abs() > 0].reset_index(drop=True), True, True)
    fit = cos2lon_fit(both3, pool_lat3, pool_lon3, rng)
    record["cos2lon_fit"] = fit

    # ---- (4) individual canonical hollows (non-trace), hot/warm test, all-fields variant
    hn = hollows[(hollows.lat > 0) & (~hollows.is_crater_trace.astype(bool))].reset_index(drop=True)
    hs = hollows[(hollows.lat < 0) & (~hollows.is_crater_trace.astype(bool))].reset_index(drop=True)
    # hollow-based size weights (D_eq)
    def hollow_weights(hdf):
        Dref = comp["Dref"]
        edges = np.concatenate([[0], (Dref[:-1] + Dref[1:]) / 2, [np.inf]])
        w, _ = np.histogram(hdf.D_eq_km.values, bins=edges)
        return w / w.sum()
    hollow_rec = {}
    for hemi, hh, name in [('+', hn, "north"), ('-', hs, "south")]:
        # reuse field-based lam but with hollow size weights: build null with a
        # temporary frame carrying D_max_nt = hollow D_eq
        tmp = pd.DataFrame({"D_max_nt": hh.D_eq_km.values})
        lam, keep, w_k = build_null(hemi, tmp, comp, cov, lat_c, lon_c)
        pool_lat, pool_lon = make_pool(lam, lat_c, lon_c, NPOOL, rng)
        rows = hotwarm_by_band(hh, pool_lat, pool_lon, hemi, rng)
        hollow_rec[name] = {"n": int(len(hh)),
                            "size_weights": [round(float(x), 4) for x in w_k],
                            "bands": rows}
    record["individual_hollows"] = {
        "note": "Robustness to field definition: same hot/warm test on individual non-trace canonical hollows, null size-weighted by hollow D_eq.",
        **hollow_rec}

    with open(OUT + "hotwarm_record.json", "w") as f:
        json.dump(record, f, indent=2)
    print("wrote", OUT + "hotwarm_record.json")
    return record


if __name__ == "__main__":
    rec = main()
    # brief console summary
    for v in ["all", "caloris_excl", "caloris_and_province_excl"]:
        print("===", v)
        for h in ["north", "south"]:
            print(" ", h, [(r["band"], r["hot_over_warm"], r["perm_p"]) for r in rec["variants"][v][h]["bands"]])
    print("cos2lon slope", rec["cos2lon_fit"]["slope"], rec["cos2lon_fit"]["slope_ci"],
          "p", rec["cos2lon_fit"]["slope_boot_p"])
