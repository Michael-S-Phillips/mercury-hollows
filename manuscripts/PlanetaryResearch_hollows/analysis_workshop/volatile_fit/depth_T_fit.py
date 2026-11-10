"""Task 1a. Fit ln h = a + b/T + c ln px + e ln tan(inc) to the 2,518 shadow-length depths of Blewett et al. (2016),
with standard errors clustered by image and 95% intervals from bootstrap resamples of images (2,000).
Apparent enthalpy: depth proportional to P_v^n with P_v ~ exp(-dH/RT) gives b = -n dH/R, so dH_app = -b R / n
(n = 1 for depth proportional to P_v; n = 1/2 for a diffusion-limited front, Schorghofer and Hsieh 2018, as in Section S9).
Temperatures: (i) modeled mean annual surface temperature T (column T of the data); (ii) peak (noon) surface temperature,
which the data do not contain, computed here for each hollow from its latitude and longitude as instantaneous radiative
equilibrium, T_peak = [S(r) cos(i) (1-A)/(eps sigma)]^(1/4), maximized over one Mercury solar day (3:2 resonance, zero obliquity,
e = 0.2056, a = 0.3871 AU, A = 0.07, eps = 0.95, flat ground, no conduction or shadowing). Writes depth_T_fit_record.json."""
import json, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "../../../GRL_insolation_hollows/data/test1_rows_inc.csv")
R = 8.314462618
d = pd.read_csv(D)
d["img"] = d.filename.astype("category").cat.codes
d["ltan"] = np.log(np.tan(np.radians(d.inc)))
# --- peak noon temperature from orbit geometry ------------------------------------------------------------------
a_AU, ecc, S1 = 0.387098, 0.205630, 1361.0
P_orb, P_rot = 87.9691, 58.6462
A_bond, eps, sig = 0.07, 0.95, 5.670374e-8
t = np.linspace(0, 2 * P_orb, 40001)                      # one solar day = two orbits
M = 2 * np.pi * t / P_orb
E = M.copy()
for _ in range(50): E = E - (E - ecc * np.sin(E) - M) / (1 - ecc * np.cos(E))
nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2), np.sqrt(1 - ecc) * np.cos(E / 2))
r = a_AU * (1 - ecc * np.cos(E))
lam_ss = nu - 2 * np.pi * t / P_rot                        # body-fixed subsolar longitude; 0 (hot pole) at perihelion
F = S1 / r ** 2
def tpeak(lat, lon):
    cosi = np.cos(np.radians(lat))[:, None] * np.cos(np.radians(lon)[:, None] - lam_ss[None, :])
    fl = (F[None, :] * np.clip(cosi, 0, None)).max(axis=1)
    return (fl * (1 - A_bond) / (eps * sig)) ** 0.25
d["Tpk"] = np.concatenate([tpeak(d.latitude.values[i:i + 200], d.longitude.values[i:i + 200]) for i in range(0, len(d), 200)])
rec = {"n": int(len(d)), "n_img": int(d.img.nunique()),
       "T_range": [float(d["T"].min()), float(d["T"].max())], "Tpk_range": [float(d.Tpk.min()), float(d.Tpk.max())],
       "Tpk_check_K": {"equator_hot_pole_0E": float(tpeak(np.array([0.]), np.array([0.]))[0]),
                       "equator_warm_pole_90E": float(tpeak(np.array([0.]), np.array([90.]))[0])},
       "corr_T_Tpk": float(np.corrcoef(d["T"], d.Tpk)[0, 1]),
       "corr_invT_invTpk": float(np.corrcoef(1 / d["T"], 1 / d.Tpk)[0, 1])}
# OLS of ln Tpk on ln T: if Tpk ~ T^k the coefficient on 1/Tpk ~ b * (Tpk/T)... report the median ratio too
rec["median_Tpk_over_T"] = float((d.Tpk / d["T"]).median())
def fitb(df, col):
    df = df.assign(x=1.0 / df[col])
    m = smf.ols("lnh ~ x + lnpx + ltan", data=df).fit(cov_type="cluster", cov_kwds={"groups": df.img})
    return m
rng = np.random.default_rng(20261102)
imgs = d.img.unique(); groups = {i: np.where(d.img.values == i)[0] for i in imgs}
NB = 2000
for col, lab in (("T", "mean_annual"), ("Tpk", "peak_noon")):
    m = fitb(d, col)
    b, se = float(m.params["x"]), float(m.bse["x"])
    bs = []
    for _ in range(NB):
        idx = np.concatenate([groups[i] for i in rng.choice(imgs, len(imgs))])
        s = d.iloc[idx]
        X = np.column_stack([np.ones(len(s)), 1.0 / s[col].values, s.lnpx.values, s.ltan.values])
        bs.append(np.linalg.lstsq(X, s.lnh.values, rcond=None)[0][1])
    bs = np.array(bs); lo, hi = np.percentile(bs, [2.5, 97.5])
    t0, t1 = float(d[col].min()), float(d[col].max()); dx = 1 / t0 - 1 / t1
    out = {"b_K": b, "se_cluster_K": se, "p_cluster": float(m.pvalues["x"]),
           "b_boot_95": [float(lo), float(hi)], "b_boot_se": float(bs.std(ddof=1)),
           "lnpx": float(m.params["lnpx"]), "ltan": float(m.params["ltan"]),
           "T_range": [t0, t1],
           "fold_change_over_range": [float(np.exp(-b * dx)), float(np.exp(-hi * dx)), float(np.exp(-lo * dx))],
           "dH_app_kJ_n1": [-b * R / 1e3, -hi * R / 1e3, -lo * R / 1e3],
           "dH_app_kJ_n0.5": [-2 * b * R / 1e3, -2 * hi * R / 1e3, -2 * lo * R / 1e3],
           "dH_app_kJ_n1_clusterCI": [-(b) * R / 1e3, -(b - 1.96 * se) * R / 1e3, -(b + 1.96 * se) * R / 1e3]}
    rec[lab] = out
# the paper's diffusion-limited sulfur front for comparison (B = -9,029 K, slope B/2)
dx = 1 / rec["T_range"][0] - 1 / rec["T_range"][1]
rec["paper_sulfur_front_fold"] = float(np.exp(9028.6679 / 2 * dx))
d[["filename", "latitude", "longitude", "T", "Tpk"]].to_csv(os.path.join(HERE, "hollow_peak_temperature.csv"), index=False)
json.dump(rec, open(os.path.join(HERE, "depth_T_fit_record.json"), "w"), indent=1)
print(json.dumps(rec, indent=1))
