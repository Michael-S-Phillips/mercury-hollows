# g3_catalog copy of analysis_workshop/completeness_sfd/sfd_windows.py (6 Nov 2026): only input/output paths and env settings changed.
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")
CAT = _os.environ.get("G3_CAT", "v3"); THETA = _os.environ.get("G3_THETA", "fixed")
ODIR = f"{G3}/sfd/out_{CAT}_{THETA}"; _os.makedirs(ODIR, exist_ok=True)
THETA_JSON = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json")
              if THETA == "fixed" else f"{G3}/sfd/out_{CAT}_refit/theta_out_v2.json")
import numpy as np, json
from scipy.optimize import minimize, minimize_scalar
U = ODIR + "/"
Z = np.load(U + "sfd_grid.npz"); Dg, Ci, D = Z["Dg"], Z["Ci"].astype(float), Z["D"]; lDg = np.log(Dg); dl = lDg[1] - lDg[0]; lD = np.log(D)
def pw_nll(a, D1, D2):
    m = (D >= D1) & (D < D2); g = (Dg >= D1) & (Dg < D2); den = (Ci[m][:, g] * Dg[g] ** (-a)).sum(1) * dl
    return -np.sum(-a * lD[m] - np.log(den)), m.sum()
out = {}
for D1, D2 in [(0.2, 0.5), (0.3, 1.0), (0.5, 1.5), (1.0, 3.0), (2.0, 6.0), (3.0, 12.0)]:
    r = minimize_scalar(lambda a: pw_nll(a, D1, D2)[0], bounds=(-1, 6), method="bounded")
    # 95% interval from the likelihood profile
    f0 = r.fun; grid = np.linspace(r.x - 1.5, r.x + 1.5, 301); prof = np.array([pw_nll(a, D1, D2)[0] for a in grid]); ins = grid[prof - f0 <= 1.92]
    out[f"{D1}-{D2}"] = (float(round(r.x, 2)), float(round(ins.min(), 2)), float(round(ins.max(), 2)), int(pw_nll(r.x, D1, D2)[1]))
print("power-law slope by window (cumulative-equivalent a; 95% profile CI; n):"); [print("  ", k, v) for k, v in out.items()]
# D >= 0.3: power law vs truncated lognormal (2 params) vs curved power law (log-quadratic, 2 params)
m = D >= 0.3; g = Dg >= 0.3
def nll_f(logf_grid, logf_obs):
    den = (Ci[m][:, g] * np.exp(logf_grid)).sum(1) * dl; return -np.sum(logf_obs - np.log(den))
pl = minimize_scalar(lambda a: nll_f(-a * lDg[g], -a * lD[m]), bounds=(0, 5), method="bounded")
ln = minimize(lambda t: nll_f(-0.5 * ((lDg[g] - t[0]) / np.exp(t[1])) ** 2, -0.5 * ((lD[m] - t[0]) / np.exp(t[1])) ** 2), [np.log(0.03), np.log(1.5)], method="Nelder-Mead")
cq = minimize(lambda t: nll_f(-t[0] * lDg[g] - t[1] * lDg[g] ** 2, -t[0] * lD[m] - t[1] * lD[m] ** 2), [1.6, 0.0], method="Nelder-Mead")
print("D>=0.3 km, n =", m.sum())
print("  power law: a = %.3f, nll %.1f, AIC %.1f" % (pl.x, pl.fun, 2 * pl.fun + 2))
print("  lognormal: median %.3f km, sigma %.2f, nll %.1f, AIC %.1f" % (np.exp(ln.x[0]), np.exp(ln.x[1]), ln.fun, 2 * ln.fun + 4))
print("  curved power law: a = %.3f + %.3f ln D (local slope at 0.3/1/3 km: %.2f / %.2f / %.2f), nll %.1f, AIC %.1f" % (cq.x[0], 2 * cq.x[1], cq.x[0] + 2 * cq.x[1] * np.log(0.3), cq.x[0], cq.x[0] + 2 * cq.x[1] * np.log(3), cq.fun, 2 * cq.fun + 4))
json.dump({"windows": out, "pl": [float(pl.x), float(pl.fun)], "lognormal": [float(np.exp(ln.x[0])), float(np.exp(ln.x[1])), float(ln.fun)], "curved": [float(cq.x[0]), float(cq.x[1]), float(cq.fun)]}, open(f"{ODIR}/sfd_windows.json", "w"), indent=1)
