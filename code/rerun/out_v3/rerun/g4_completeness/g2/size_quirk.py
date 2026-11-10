# g4 copy of code/rerun/out_v3/rerun/g2_insolation/size_quirk.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g4.py).
#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/insolation_v6/size_quirk.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Task (b), 2 Nov 2026: is the hot-longitude excess of hollows within 15 deg of the equator a quirk of the crater
size distribution (more small craters at the equatorial hot longitudes)?

Craters: analysis_workshop/summary_model/crater_features.csv (31,554 craters >= 5 km of Herrick et al. 2018, possible
peak rings dropped; C1 = completeness for 1 km hollows at the crater center; h09 = >= 1 hollow within 0.9 radii).
Hot = within 45 deg of 0 or 180 E; warm = within 45 deg of 90 or 270 E. Band |lat| < 15.

(i)   Crater density by sector and size bin, per unit area and completeness-weighted (sum of C1 per unit area).
(ii)  Hot/warm ratio of observed/expected hollow fields, |lat| < 15, Caloris < 1,200 km excluded, with the expectation
      from (N0) the area null of hotwarm.py, (N1) completeness-weighted crater counts, (N2) completeness-weighted crater
      area (C1 x D^2), (N3) the hosting probability of each crater from a logistic model with no longitude term
      (h09 ~ logD + logD^2 + C1 + abslat, fitted to all craters within 60 deg). 2,000 bootstrap resamples of fields.
(iii) Crater-level logistic models of h09 on hot sector, logD and C1 (+ abslat) in each size bin; odds ratios.
Writes size_quirk.json.
"""
import json, numpy as np, pandas as pd
import statsmodels.formula.api as smf
from common_null import *

rng = np.random.default_rng(20261102)
c = pd.read_csv(f"{IN}/crater_features.csv", low_memory=False)
c["hot"] = hot_mask(c.lon.values).astype(int)
c["dcal"] = gc_dist_km(CALORIS[0], CALORIS[1], c.lat.values, c.lon.values)
c["logD2"] = c.logD ** 2
BINS = [(5, 10), (10, 20), (20, 50), (50, 1e5)]
lab = lambda lo, hi: f"{lo}-{hi}" if hi < 1e4 else f">={lo}"
OUT = {"description": __doc__}

# ---------- areas of the hot and warm sectors in the band, with and without the Caloris exclusion (0.1 deg cells)
lat_c, lon_c = grid_cell_centers(10)
LAT = np.broadcast_to(lat_c[:, None], (1800, 3600)); LON = np.broadcast_to(lon_c[None, :], (1800, 3600))
cellA = (np.radians(0.1) * R_MERCURY) ** 2 * np.cos(np.radians(LAT))
band = np.abs(LAT) < 15; hotc = hot_mask(LON); calc = gc_dist_km(CALORIS[0], CALORIS[1], LAT, LON) > CALORIS_KM
cov, comp, fields, hollows = load()
C1grid = comp["C"][list(comp["Dref"]).index(1.0)]
AREA = {}
for nm, m in (("all", band), ("noCaloris", band & calc)):
    AREA[nm] = dict(hot=float(cellA[m & hotc].sum()), warm=float(cellA[m & ~hotc].sum()),
                    hot_C1area=float((cellA * np.nan_to_num(C1grid))[m & hotc].sum()),
                    warm_C1area=float((cellA * np.nan_to_num(C1grid))[m & ~hotc].sum()))
OUT["areas_km2"] = AREA

# ---------- (i) crater densities
def dens(e, A):
    rows = {}
    for lo, hi in BINS + [(5, 1e5)]:
        s = e[(e.D >= lo) & (e.D < hi)]
        r = {}
        for sec, k in (("hot", 1), ("warm", 0)):
            t = s[s.hot == k]
            r[sec] = dict(n=int(len(t)), per_1e6km2=float(len(t) / A[sec] * 1e6), C1_sum=float(t.C1.sum()),
                          C1w_per_1e6km2=float(t.C1.sum() / A[sec] * 1e6), mean_C1=float(t.C1.mean()),
                          n_per_1e6km2_C1area=float(len(t) / A[sec + "_C1area"] * 1e6),
                          h09=int(t.h09.sum()), hosts11=int(t.hosts.sum()),
                          h09_per_C1sum=float(t.h09.sum() / t.C1.sum()))
        r["hot_over_warm_density"] = r["hot"]["per_1e6km2"] / r["warm"]["per_1e6km2"]
        r["hot_over_warm_C1w_density"] = r["hot"]["C1w_per_1e6km2"] / r["warm"]["C1w_per_1e6km2"]
        r["hot_over_warm_density_per_C1area"] = r["hot"]["n_per_1e6km2_C1area"] / r["warm"]["n_per_1e6km2_C1area"]
        r["hot_over_warm_h09_per_C1sum"] = r["hot"]["h09_per_C1sum"] / r["warm"]["h09_per_C1sum"]
        rows[lab(lo, hi)] = r
    return rows
E = c[c.abslat < 15]; Enc = E[E.dcal > CALORIS_KM]
OUT["i_density_all"] = dens(E, AREA["all"]); OUT["i_density_noCaloris"] = dens(Enc, AREA["noCaloris"])
# share of small craters in each sector
OUT["i_frac_5_10km"] = {s: float((E[E.hot == k].D < 10).mean()) for s, k in (("hot", 1), ("warm", 0))}

# ---------- (ii) field obs/exp hot/warm with crater-based nulls
fn = fields[fields.lat > 0].reset_index(drop=True); fs = fields[fields.lat < 0].reset_index(drop=True)
lamN, *_ = build_null('+', fn, comp, cov, lat_c, lon_c, True, True); lamS, *_ = build_null('-', fs, comp, cov, lat_c, lon_c, True, True)
obs = apply_point_exclusions(fields[fields.lat.abs() > 0].reset_index(drop=True), True, True)
obs = obs[obs.lat.abs() < 15]; oh = hot_mask(obs.lon.values)
def eh_area():
    num = den = 0.0
    for lam in (lamN, lamS):
        P = lam / lam.sum(); m = band
        num += P[m & hotc].sum(); den += P[m].sum()
    return num / den
c60 = c[c.abslat <= 60].copy()
base = smf.logit("h09 ~ logD + logD2 + C1 + abslat", data=c60).fit(disp=0)
Enc = Enc.copy(); Enc["p_base"] = base.predict(Enc)
nulls = {"N0_area_coverage_completeness": eh_area(),
         "N1_C1_weighted_crater_counts": float(Enc[Enc.hot == 1].C1.sum() / Enc.C1.sum()),
         "N1b_crater_counts_unweighted": float(Enc.hot.mean()),
         "N2_C1_weighted_crater_area": float((Enc.C1 * Enc.D ** 2)[Enc.hot == 1].sum() / (Enc.C1 * Enc.D ** 2).sum()),
         "N2b_C1_weighted_crater_area_D_lt_100": float((Enc.C1 * Enc.D ** 2)[(Enc.hot == 1) & (Enc.D < 100)].sum() / (Enc.C1 * Enc.D ** 2)[Enc.D < 100].sum()),
         "N3_predicted_hosting_no_longitude": float(Enc[Enc.hot == 1].p_base.sum() / Enc.p_base.sum())}
res = {}
fh = oh.mean()
for k, eh in nulls.items():
    hw = (fh / eh) / ((1 - fh) / (1 - eh))
    bt = []
    for i in range(2000):
        s = oh[rng.integers(0, len(oh), len(oh))].mean(); bt.append((s / eh) / ((1 - s) / (1 - eh)))
    # p: draws of the same number of fields from the null
    perm = rng.binomial(len(oh), eh, 4000) / len(oh); pr = (perm / eh) / ((1 - perm) / (1 - eh))
    p = (1 + np.sum(np.abs(np.log(pr)) >= abs(np.log(hw)))) / (1 + len(pr))
    res[k] = dict(null_hot_share=float(eh), hot_over_warm=float(hw), ci=[float(v) for v in np.percentile(bt, [2.5, 97.5])], p=float(p))
OUT["ii_fields"] = dict(n_fields=int(len(oh)), observed_hot_share=float(fh), base_model=str(base.model.formula),
                        base_model_params={k: float(v) for k, v in base.params.items()}, nulls=res)
# the same for h09 host craters (observed hosts vs expected) as a cross-check
hh = Enc[Enc.h09 == 1]; fhh = hh.hot.mean(); res2 = {}
for k in ("N1_C1_weighted_crater_counts", "N2_C1_weighted_crater_area", "N3_predicted_hosting_no_longitude"):
    eh = nulls[k]; hw = (fhh / eh) / ((1 - fhh) / (1 - eh)); hv = hh.hot.values
    bt = [(lambda s: (s / eh) / ((1 - s) / (1 - eh)))(hv[rng.integers(0, len(hv), len(hv))].mean()) for _ in range(2000)]
    res2[k] = dict(hot_over_warm=float(hw), ci=[float(v) for v in np.percentile(bt, [2.5, 97.5])])
OUT["ii_hosts_h09"] = dict(n_hosts=int(len(hh)), observed_hot_share=float(fhh), nulls=res2)

# ---------- (iii) crater-level logistic models within |lat| < 15
def orr(m, term="hot"):
    b, se = m.params[term], m.bse[term]
    return dict(OR=float(np.exp(b)), lo=float(np.exp(b - 1.96 * se)), hi=float(np.exp(b + 1.96 * se)), p=float(m.pvalues[term]))
iii = {}
for nm, d in (("all", E), ("noCaloris", E[E.dcal > CALORIS_KM])):
    rr = {}
    for lo, hi in BINS:
        s = d[(d.D >= lo) & (d.D < hi)]
        r = dict(n=int(len(s)), hosts_hot=int(s[s.hot == 1].h09.sum()), hosts_warm=int(s[s.hot == 0].h09.sum()))
        for f in ("h09 ~ hot + logD + C1", "h09 ~ hot + logD + C1 + abslat"):
            try: r[f] = orr(smf.logit(f, data=s).fit(disp=0))
            except Exception as ex: r[f] = str(ex)
        rr[lab(lo, hi)] = r
    m = smf.logit("h09 ~ hot + logD + C1 + abslat", data=d).fit(disp=0); rr["all_sizes"] = orr(m)
    mi = smf.logit("h09 ~ hot * logD + C1 + abslat", data=d).fit(disp=0)
    rr["hot_x_logD"] = dict(coef=float(mi.params["hot:logD"]), se=float(mi.bse["hot:logD"]), p=float(mi.pvalues["hot:logD"]))
    rr["all_sizes_hosts11"] = orr(smf.logit("hosts ~ hot + logD + C1 + abslat", data=d).fit(disp=0))
    iii[nm] = rr
OUT["iii_crater_logit_lat15"] = iii
# same, 15-60 N for contrast
N = c[(c.lat >= 15) & (c.lat <= 60)]
OUT["iii_contrast_15_60N"] = {**{lab(lo, hi): orr(smf.logit("h09 ~ hot + logD + C1 + abslat", data=N[(N.D >= lo) & (N.D < hi)]).fit(disp=0)) for lo, hi in BINS},
                              "all_sizes": orr(smf.logit("h09 ~ hot + logD + C1 + abslat", data=N).fit(disp=0)),
                              "n": int(len(N)), "hosts": int(N.h09.sum())}
json.dump(OUT, open(f"{G2_OUT}/size_quirk.json", "w"), indent=1)

# ---------- print
print("areas", {k: {kk: round(vv / 1e6, 3) for kk, vv in v.items()} for k, v in AREA.items()})
for nm in ("i_density_all", "i_density_noCaloris"):
    print("==", nm)
    for b, r in OUT[nm].items():
        print(f"  {b:7s} hot n={r['hot']['n']:5d} ({r['hot']['per_1e6km2']:.0f}/1e6km2, C1w {r['hot']['C1w_per_1e6km2']:.1f}, meanC1 {r['hot']['mean_C1']:.3f}, h09 {r['hot']['h09']})"
              f" | warm n={r['warm']['n']:5d} ({r['warm']['per_1e6km2']:.0f}, C1w {r['warm']['C1w_per_1e6km2']:.1f}, meanC1 {r['warm']['mean_C1']:.3f}, h09 {r['warm']['h09']})"
              f" | h/w dens {r['hot_over_warm_density']:.2f} C1w {r['hot_over_warm_C1w_density']:.2f} per-C1area {r['hot_over_warm_density_per_C1area']:.2f} hosts/C1sum {r['hot_over_warm_h09_per_C1sum']:.2f}")
print("frac 5-10 km", OUT["i_frac_5_10km"])
print("== ii fields n", OUT["ii_fields"]["n_fields"], "obs hot share", round(fh, 3))
for k, v in res.items(): print(f"  {k:45s} null hot {v['null_hot_share']:.3f}  hot/warm {v['hot_over_warm']:.2f} ({v['ci'][0]:.2f}-{v['ci'][1]:.2f}) p={v['p']:.4f}")
print("== ii hosts n", len(hh), "obs hot share", round(fhh, 3))
for k, v in res2.items(): print(f"  {k:45s} hot/warm {v['hot_over_warm']:.2f} ({v['ci'][0]:.2f}-{v['ci'][1]:.2f})")
for nm, rr in iii.items():
    print("== iii", nm)
    for k, v in rr.items(): print("  ", k, v)
print("== contrast 15-60N", OUT["iii_contrast_15_60N"])
