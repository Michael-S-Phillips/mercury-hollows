#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/insolation_v6/cos2lon_lat15.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Task (a), 2 Nov 2026: longitude profile of hollow fields within 15 deg of the equator (Figure 3a for paper v6).

Same method as cos2lon_fit in code/rerun/latitude/hotwarm/hotwarm.py (Figure 3a of v5: |lat| < 30, 642 fields):
observed / expected hollow fields in longitude bins, both hemispheres, fields within 1,200 km of the Caloris center
(and in the 90-150E, 40-70N province) excluded; weighted least squares fit of a + b cos 2 lon (weights = expected
counts); 95% interval on b from 2,000 bootstrap resamples of fields. Expected counts here are exact sums of the
per-hemisphere null intensity (each hemisphere normalized to the same total, as the equal 400,000-point pools of
hotwarm.py did). The |lat| < 30 case is recomputed as a check. Variants: hemisphere-stratified expectation, and
15 and 45 deg bins. Writes cos2lon_fit_lat15.json.
"""
import json, numpy as np
from common_null import *

SEED = 20261102
lat_c, lon_c = grid_cell_centers(10)
cov, comp, fields, hollows = load()
fn = fields[fields.lat > 0].reset_index(drop=True); fs = fields[fields.lat < 0].reset_index(drop=True)
lamN, *_ = build_null('+', fn, comp, cov, lat_c, lon_c, True, True)
lamS, *_ = build_null('-', fs, comp, cov, lat_c, lon_c, True, True)
both = apply_point_exclusions(fields[fields.lat.abs() > 0].reset_index(drop=True), True, True)
LON = np.broadcast_to(lon_c[None, :], lamN.shape); LAT = np.broadcast_to(lat_c[:, None], lamN.shape)


def expected_frac(latmax, edges, stratified_counts=None):
    """Null share per longitude bin within |lat| < latmax."""
    out = np.zeros(len(edges) - 1)
    for lam, sgn in ((lamN, 1), (lamS, -1)):
        P = lam / lam.sum()
        band = (np.abs(LAT) < latmax) & (np.sign(LAT) == sgn)
        h, _ = np.histogram(LON[band], bins=edges, weights=P[band])
        if stratified_counts is not None:            # expected = n_obs(hemisphere) x null share within hemisphere band
            h = h / h.sum() * stratified_counts[sgn]
        out += h
    return out / out.sum()


def fit(latmax, nbins=12, stratified=False, nboot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    obs = both[np.abs(both.lat.values) < latmax]
    edges = np.linspace(0, 360, nbins + 1); ctr = (edges[:-1] + edges[1:]) / 2; x = np.cos(np.radians(2 * ctr))
    sc = {1: int((obs.lat > 0).sum()), -1: int((obs.lat < 0).sum())} if stratified else None
    exp = expected_frac(latmax, edges, sc) * len(obs)
    lon_o = np.mod(obs.lon.values, 360)
    oc, _ = np.histogram(lon_o, bins=edges); ratio = oc / exp
    a, b = wls_cos2(ratio, exp, x)
    sl, ic = np.empty(nboot), np.empty(nboot)
    for i in range(nboot):
        o, _ = np.histogram(lon_o[rng.integers(0, len(obs), len(obs))], bins=edges)
        ic[i], sl[i] = wls_cos2(o / exp, exp, x)
    # sector hot/warm ratio of observed/expected (hot = within 45 deg of 0/180 E), with bootstrap
    e5 = np.array([0, 45, 135, 225, 315, 360.])
    f5 = expected_frac(latmax, e5, sc); eh = f5[[0, 2, 4]].sum()
    oh = hot_mask(lon_o)
    hw = (oh.mean() / eh) / ((1 - oh.mean()) / (1 - eh))
    bt = []
    for i in range(nboot):
        s = oh[rng.integers(0, len(oh), len(oh))].mean(); bt.append((s / eh) / ((1 - s) / (1 - eh)))
    return dict(latmax=latmax, nbins=nbins, stratified=stratified, n=int(len(obs)), n_north=int((obs.lat > 0).sum()),
                n_south=int((obs.lat < 0).sum()), intercept=float(a), slope=float(b),
                slope_ci=[float(v) for v in np.percentile(sl, [2.5, 97.5])],
                intercept_ci=[float(v) for v in np.percentile(ic, [2.5, 97.5])],
                slope_boot_p=float(2 * min((sl <= 0).mean(), (sl >= 0).mean())),
                fitted_hot=float(a + b), fitted_warm=float(a - b), fitted_hot_over_warm=float((a + b) / (a - b)),
                sector_null_hot_frac=float(eh), sector_hot_over_warm=float(hw),
                sector_hot_over_warm_ci=[float(v) for v in np.percentile(bt, [2.5, 97.5])],
                lon_bin_centers=ctr.tolist(), observed=oc.tolist(), expected=exp.tolist(),
                obs_over_null_ratio=ratio.tolist(), cos2lon=x.tolist())


rec = {"description": __doc__, "check_lat30_vs_hotwarm_record": fit(30.0),
       "lat15": fit(15.0), "lat15_stratified": fit(15.0, stratified=True),
       "lat15_bins15deg": fit(15.0, nbins=24), "lat15_bins45deg": fit(15.0, nbins=8),
       "lat30_stratified": fit(30.0, stratified=True)}
json.dump(rec, open(f"{G2_OUT}/cos2lon_fit_lat15.json", "w"), indent=1)
for k, v in rec.items():
    if k == "description": continue
    print(f"{k:32s} n={v['n']} (N {v['n_north']}, S {v['n_south']})  a={v['intercept']:.3f}  b={v['slope']:.3f} "
          f"({v['slope_ci'][0]:.3f}-{v['slope_ci'][1]:.3f}) p={v['slope_boot_p']:.4f}  hot/warm fit {v['fitted_hot_over_warm']:.2f}"
          f"  sector {v['sector_hot_over_warm']:.2f} ({v['sector_hot_over_warm_ci'][0]:.2f}-{v['sector_hot_over_warm_ci'][1]:.2f})")
print("lat15 ratios", np.round(rec["lat15"]["obs_over_null_ratio"], 2), rec["lat15"]["observed"])
