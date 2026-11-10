# g3_catalog copy of analysis_workshop/completeness_sfd/sfd_numbers.py (6 Nov 2026): only input/output paths and env settings changed.
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")
CAT = _os.environ.get("G3_CAT", "v3"); THETA = _os.environ.get("G3_THETA", "fixed")
ODIR = f"{G3}/sfd/out_{CAT}_{THETA}"; _os.makedirs(ODIR, exist_ok=True)
THETA_JSON = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json")
              if THETA == "fixed" else f"{G3}/sfd/out_{CAT}_refit/theta_out_v2.json")
import numpy as np, json
from scipy.optimize import minimize_scalar
U = ODIR + "/"
Z = np.load(U + "sfd_grid.npz"); Dg, Ci, D = Z["Dg"], Z["Ci"].astype(float), Z["D"]; lDg = np.log(Dg); dl = lDg[1] - lDg[0]; lD = np.log(D)
def fitwin(D1, D2, sel=None):
    s = np.ones(len(D), bool) if sel is None else sel
    def nll(a):
        m = s & (D >= D1) & (D < D2); g = (Dg >= D1) & (Dg < D2); den = (Ci[m][:, g] * Dg[g] ** (-a)).sum(1) * dl
        return -np.sum(-a * lD[m] - np.log(den))
    r = minimize_scalar(nll, bounds=(-3, 7), method="bounded"); gr = np.linspace(r.x - 2, r.x + 2, 401); pr = np.array([nll(a) for a in gr]); ins = gr[pr - r.fun <= 1.92]
    return round(float(r.x), 2), round(float(ins.min()), 2), round(float(ins.max()), 2), int((s & (D >= D1) & (D < D2)).sum())
R = {}
R["slope_0.7_3"] = fitwin(0.7, 3.0); R["slope_3_12"] = fitwin(3.0, 12.0); R["slope_1_3"] = fitwin(1.0, 3.0); R["slope_0.3_0.7"] = fitwin(0.3, 0.7)
Cown = np.array([np.interp(np.log(x), lDg, r) for x, r in zip(D, Ci)]); w = 1 / np.clip(Cown, 1e-3, 1)
C03 = np.array([np.interp(np.log(0.3), lDg, r) for r in Ci]); C1 = np.array([np.interp(0.0, lDg, r) for r in Ci])
rng = np.random.default_rng(1); idx = [rng.integers(0, len(D), len(D)) for _ in range(500)]
a, aci = R["slope_0.7_3"][0], R["slope_0.7_3"][1:3]
def counts(Dthr):
    obs = int((D > Dthr).sum()); cor = w[D > Dthr].sum(); cb = [w[i][D[i] > Dthr].sum() for i in idx]
    return obs, round(float(cor)), [round(float(v)) for v in np.percentile(cb, [2.5, 97.5])]
# power law anchored on corrected count > 1 km, with the 0.7-3 km slope
N1 = w[D > 1].sum()
R["counts"] = {}
for Dt in (0.3, 0.5, 1.0, 3.0, 5.0):
    o, c, cci = counts(Dt); R["counts"][Dt] = dict(catalog=o, corrected=c, corrected_ci=cci, powerlaw=round(float(N1 * Dt ** (-a))), powerlaw_range=[round(float(N1 * Dt ** (-s))) for s in aci])
# well-imaged sites: detection probability of a 0.3 km hollow >= 0.5 at the site
for thr in (0.5, 0.3):
    s = C03 >= thr
    R[f"sites_C03_ge_{thr}"] = dict(n=int(s.sum()), n_ge03=int((s & (D >= 0.3)).sum()), median_C_own=float(np.median(Cown[s & (D >= 0.3)])),
        slope_03_1=fitwin(0.3, 1.0, s), slope_03_07=fitwin(0.3, 0.7, s), slope_1_3=fitwin(1.0, 3.0, s),
        obs_counts={Dt: int((s & (D > Dt)).sum()) for Dt in (0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0)})
R["median_C_at_size"] = {Dt: round(float(np.median([np.interp(np.log(Dt), lDg, r) for r in Ci])), 3) for Dt in (0.1, 0.3, 0.5, 1, 3, 5)}
print(json.dumps(R, indent=1)); json.dump(R, open(f"{ODIR}/sfd_numbers.json", "w"), indent=1)
