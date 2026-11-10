# g3_catalog copy of analysis_workshop/completeness_sfd/fit_sfd.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Size distribution of hollows corrected for detection (4 Nov 2026; Michael's question on a completeness measure that does not
depend on one reference size). For each cataloged hollow, the probability that a hollow of diameter D at its location is
detected in at least one image is C_i(D) = 1 - prod_k (1 - p(1000 D / px_k, inc_k)) over the searched images that cover it
(code/rerun/out_v2/completeness_pairs.csv), with p the per-image detection model (theta_out_v2.json) at each image's own pixel scale and incidence. The size distribution f(D) is fitted by
maximum likelihood conditional on location: L = prod_i f(D_i) C_i(D_i) / int f(D) C_i(D) dD. Forms: lognormal (all sizes) and
power law (differential, D >= Dmin, Dmin = 0.2, 0.3, 0.5 km). Writes sfd_fit.json and sfd_grid.npz (C_i on the D grid)."""
import os, json, numpy as np, pandas as pd
from scipy.special import erf
from scipy.optimize import minimize
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")
CAT = _os.environ.get("G3_CAT", "v3"); THETA = _os.environ.get("G3_THETA", "fixed")
ODIR = f"{G3}/sfd/out_{CAT}_{THETA}"; _os.makedirs(ODIR, exist_ok=True)
THETA_JSON = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json")
              if THETA == "fixed" else f"{G3}/sfd/out_{CAT}_refit/theta_out_v2.json")
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
PF = json.load(open(THETA_JSON))   # out_v2 detection model (refit_detection.py; matches the paper, S2)
def pdet(x, inc):
    Phi = lambda z: 0.5 * (1 + erf(z / np.sqrt(2)))
    p = PF["p_max"] * Phi((np.log(x) - np.log(PF["x50_px"])) / PF["sigma_lo"]) * (1 - Phi((np.log(x) - np.log(PF["x_hi_px"])) / PF["sigma_hi"]))
    return np.clip(p * np.exp(PF["b_inc"] * (inc - 60) / 25), 1e-6, 0.999)
CV = pd.read_csv(f"{G3}/inputs/completeness_pairs_{CAT}.csv")
Dg = np.exp(np.linspace(np.log(0.03), np.log(30), 121)); lDg = np.log(Dg); dl = lDg[1] - lDg[0]
H = CV.groupby("h").agg(D=("D_km", "first")).reset_index(); hix = {h: i for i, h in enumerate(H.h)}
L = np.zeros((len(H), len(Dg)))
ii = CV.h.map(hix).values
for j, Dj in enumerate(Dg):
    np.add.at(L[:, j], ii, np.log1p(-pdet(1000 * Dj / CV.px_m.values, CV.inc.values)))
Ci = 1 - np.exp(L)                     # (hollows x D grid)
lD = np.log(H.D.values)
def interp_rows(lx):                   # C_i at each hollow's own D, for diagnostics
    k = np.clip(((lx - lDg[0]) / dl).astype(int), 0, len(Dg) - 2); w = (lx - lDg[k]) / dl
    return Ci[np.arange(len(lx)), k] * (1 - w) + Ci[np.arange(len(lx)), k + 1] * w
def nll_lognormal(th):
    mu, ls = th; s = np.exp(ls)
    fg = np.exp(-0.5 * ((lDg - mu) / s) ** 2) / s           # density in ln D (unnormalized constant cancels)
    den = (Ci * fg).sum(1) * dl
    return -np.sum(-0.5 * ((lD - mu) / s) ** 2 - np.log(s) - np.log(den))
def nll_power(th, Dmin):
    a = th[0]; m = H.D.values >= Dmin; g = Dg >= Dmin
    fg = Dg[g] ** (-a)                                       # density in ln D: D^-a  (differential slope a+1 in D)
    den = (Ci[m][:, g] * fg).sum(1) * dl
    return -np.sum(-a * lD[m] - np.log(den))
R = {"n_hollows": int(len(H)), "D_grid_km": [float(Dg[0]), float(Dg[-1]), len(Dg)], "theta": PF}
r = minimize(nll_lognormal, [np.log(0.5), np.log(0.8)], method="Nelder-Mead")
mu, s = r.x[0], np.exp(r.x[1]); R["lognormal"] = {"median_km": float(np.exp(mu)), "sigma_lnD": float(s), "nll": float(r.fun)}
# observed-size check: predicted distribution of detected sizes vs observed (quantiles)
fg = np.exp(-0.5 * ((lDg - mu) / s) ** 2); pred = (Ci * fg) ; pred = pred / pred.sum(1, keepdims=True); predc = pred.mean(0).cumsum()
q_obs = np.quantile(H.D.values, [0.1, 0.25, 0.5, 0.75, 0.9]); q_pred = [float(np.interp(q, predc, Dg)) for q in [0.1, 0.25, 0.5, 0.75, 0.9]]
R["lognormal"]["detected_size_quantiles_obs_vs_pred"] = [list(np.round(q_obs, 3)), list(np.round(q_pred, 3))]
R["lognormal"]["fraction_of_population_below_0.2_0.3_0.5km"] = [float(0.5 * (1 + erf((np.log(x) - mu) / (s * np.sqrt(2))))) for x in (0.2, 0.3, 0.5)]
R["power"] = {}
for Dmin in (0.2, 0.3, 0.5):
    rp = minimize(lambda t: nll_power(t, Dmin), [1.5], method="Nelder-Mead"); a = rp.x[0]
    m = H.D.values >= Dmin; g = Dg >= Dmin; fgp = Dg[g] ** (-a); pp = Ci[m][:, g] * fgp; pp = pp / pp.sum(1, keepdims=True); pc = pp.mean(0).cumsum()
    qo = np.quantile(H.D.values[m], [0.1, 0.25, 0.5, 0.75, 0.9]); qp = [float(np.interp(q, pc, Dg[g])) for q in [0.1, 0.25, 0.5, 0.75, 0.9]]
    R["power"][str(Dmin)] = {"slope_cumulative": float(a), "n": int(m.sum()), "detected_size_quantiles_obs_vs_pred": [list(np.round(qo, 3)), list(np.round(qp, 3))]}
R["observed_D_quantiles_km"] = list(np.round(np.quantile(H.D.values, [0.01, 0.1, 0.5, 0.9, 0.99]), 3))
R["median_Ci_at_own_D"] = float(np.median(interp_rows(lD)))
json.dump(R, open(f"{ODIR}/sfd_fit.json", "w"), indent=1)
np.savez_compressed(f"{ODIR}/sfd_grid.npz", Dg=Dg, Ci=Ci.astype(np.float32), D=H.D.values)
print(json.dumps({k: v for k, v in R.items() if k != "theta"}, indent=1))
