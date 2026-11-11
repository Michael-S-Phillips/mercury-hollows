# g4 copy of code/rerun/out_v3/rerun/g2_insolation/dT_hosting.py (6 Nov 2026). Changed: G2_ROOT depth 5 -> 6 (see make_copies_g4.py).
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_v3/dT_hosting.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Hot-longitude test with the modeled temperature excess dT in place of longitude distance.
dT = mean annual surface temperature at the crater minus that at the warm longitude (90 or 270 E)
at the same latitude, from the 1-D thermal model grid (thermal.py, thermal_grid.npz)."""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 6))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_INSOLDATA = f"{G2_ROOT}/manuscripts/insolation_thermal/data"
# --- end g2 block ---
import os
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
G = os.path.expanduser("~/mnt/hollows-emissivity")+"/manuscripts/insolation_thermal"
sys.path.insert(0, f"{G}/scripts"); import thermal
z = np.load(f"{G}/data/thermal_grid.npz"); la, lo, T = z["la"], z["lo"], z["G"]
s = pd.read_csv(f"{IN}/test2_craters_h18.csv"); s = s[(s.diameter >= 5) & (s.abslat <= 60)].copy()
s["Tm"] = thermal.lookup(la, lo, T, s.latitude.values, s.longitude.values)
s["Tw"] = thermal.lookup(la, lo, T, s.latitude.values, np.full(len(s), 90.0))
s["dT"] = s.Tm - s.Tw
lo180 = s.longitude.values % 180.0; s["hotness"] = np.abs(np.where(lo180 > 90, 180 - lo180, lo180) - 90)
R = {}
R["dT_equator_hot"] = float(thermal.lookup(la, lo, T, 0, 0) - thermal.lookup(la, lo, T, 0, 90))
R["T_equator_hot"] = float(thermal.lookup(la, lo, T, 0, 0)); R["T_equator_warm"] = float(thermal.lookup(la, lo, T, 0, 90))
for a in (15, 30, 45, 60):
    R[f"dT_hot_lat{a}"] = float(thermal.lookup(la, lo, T, a, 0) - thermal.lookup(la, lo, T, a, 90))
R["n_craters"] = int(len(s)); R["n_hosts"] = int(s.hosts.sum()); R["Tmean_check_maxdiff"] = float(np.abs(s.Tm - s.Tmean).max())
m = smf.logit("hosts ~ dT * logD + abslat", data=s).fit(disp=0)
mh = smf.logit("hosts ~ hotness * logD + abslat", data=s).fit(disp=0)
b, bi = m.params["dT"], m.params["dT:logD"]
R["dT_model"] = dict(params={k: float(v) for k, v in m.params.items()}, bse={k: float(v) for k, v in m.bse.items()},
    p={k: float(v) for k, v in m.pvalues.items()}, aic=float(m.aic),
    OR10K_10km=float(np.exp(10 * (b + bi * 1))), OR10K_100km=float(np.exp(10 * (b + bi * 2))),
    ORfull_eq_10km=float(np.exp(R["dT_equator_hot"] * (b + bi * 1))), ORfull_eq_100km=float(np.exp(R["dT_equator_hot"] * (b + bi * 2))),
    zero_km=float(10 ** (-b / bi)))
# CI for OR per 10 K at 10 and 100 km
V = m.cov_params().loc[["dT", "dT:logD"], ["dT", "dT:logD"]].values
for D in (10, 100):
    x = np.array([1, np.log10(D)]); c = b + bi * np.log10(D); se = float(np.sqrt(x @ V @ x))
    R["dT_model"][f"OR10K_{D}km_ci"] = [float(np.exp(10 * (c - 1.96 * se))), float(np.exp(10 * (c + 1.96 * se)))]
R["hotness_model"] = dict(aic=float(mh.aic), inter=float(mh.params["hotness:logD"]), inter_se=float(mh.bse["hotness:logD"]), p=float(mh.pvalues["hotness:logD"]),
    OR_10km=float(np.exp(90 * (mh.params["hotness"] + mh.params["hotness:logD"]))), OR_100km=float(np.exp(90 * (mh.params["hotness"] + 2 * mh.params["hotness:logD"]))),
    zero_km=float(10 ** (-mh.params["hotness"] / mh.params["hotness:logD"])))
# size bins
bins = [5, 10, 20, 40, 80, 160, 800]; R["bins"] = []
for lo_, hi_ in zip(bins[:-1], bins[1:]):
    d2 = s[(s.diameter >= lo_) & (s.diameter < hi_)]
    mm = smf.logit("hosts ~ dT + abslat", data=d2).fit(disp=0)
    R["bins"].append(dict(lo=lo_, hi=hi_, n=int(len(d2)), hosts=int(d2.hosts.sum()), b=float(mm.params["dT"]), se=float(mm.bse["dT"]), p=float(mm.pvalues["dT"]),
                          OR10K=float(np.exp(10 * mm.params["dT"])), ci=[float(np.exp(10 * (mm.params["dT"] - 1.96 * mm.bse["dT"]))), float(np.exp(10 * (mm.params["dT"] + 1.96 * mm.bse["dT"])))]))
pass
json.dump(R, open(f"{G2_OUT}/dT_record.json", "w"), indent=1)
print(json.dumps(R, indent=1))
