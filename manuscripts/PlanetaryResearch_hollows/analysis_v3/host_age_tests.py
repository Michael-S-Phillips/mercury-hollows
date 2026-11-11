import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Host-age tests of the 'early burst + long-term activity' idea (27 Oct 2026, exploratory).
Herrick et al. (2018) preservation codes in column 'preserved': f, s, d (read here as fresh, intermediate,
degraded; confirm against the catalog documentation). Rays from column 'ejecta' containing 'r'.
(1) odds of hosting hollows by preservation class, at fixed log diameter and |latitude| (craters >= 10 km);
(2) warm-longitude preference (odds per 10 K of dT) within each class, and the dT x class interaction;
(3) corrected hollow depth (Blewett et al. 2016; pixel scale and incidence controlled) by host class,
    standard errors clustered by image."""
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
G = _REPO + "/manuscripts/insolation_thermal"
sys.path.insert(0, G + "/scripts"); import thermal
z = np.load(G + "/data/thermal_grid.npz"); la, lo, T = z["la"], z["lo"], z["G"]
R = {}
s = pd.read_csv(G + "/data/test2_craters_h18.csv")
s = s[s.preserved.isin(["f", "s", "d"])].copy(); s["rays"] = s.ejecta.astype(str).str.contains("r").astype(int)
s["pres"] = pd.Categorical(s.preserved, categories=["d", "s", "f"])
R["counts"] = {k: dict(n=int(len(g)), hosts=int(g.hosts.sum())) for k, g in s.groupby("preserved")}
m = smf.logit("hosts ~ C(pres) + logD + abslat", data=s).fit(disp=0)
for k in ["C(pres)[T.s]", "C(pres)[T.f]"]:
    b, se = m.params[k], m.bse[k]; R["OR_" + k[-2]] = [float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), float(m.pvalues[k])]
m = smf.logit("hosts ~ rays + logD + abslat", data=s).fit(disp=0); b, se = m.params["rays"], m.bse["rays"]
R["OR_rays"] = [float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), float(m.pvalues["rays"])]; R["n_rayed"] = int(s.rays.sum())
s2 = s[s.abslat <= 60].copy()
s2["dT"] = thermal.lookup(la, lo, T, s2.latitude.values, s2.longitude.values) - thermal.lookup(la, lo, T, s2.latitude.values, np.full(len(s2), 90.0))
for k, g in s2.groupby("preserved"):
    mm = smf.logit("hosts ~ dT + logD + abslat", data=g).fit(disp=0); b, se = mm.params["dT"], mm.bse["dT"]
    R["ORper10K_" + k] = [float(np.exp(10 * b)), float(np.exp(10 * (b - 1.96 * se))), float(np.exp(10 * (b + 1.96 * se))), int(g.hosts.sum())]
mi = smf.logit("hosts ~ dT * C(pres) + logD + abslat", data=s2).fit(disp=0)
R["interaction_p"] = {k: float(v) for k, v in mi.pvalues.items() if "dT:" in k}
# depth by host class (reuses the hollow-to-host assignment of depth_structure.py)
exec(open(_REPO + "/manuscripts/insolation_thermal/depthctl/depth_structure.py").read().split("R = {")[0])
h = d[d.host >= 0].copy(); h["pres"] = cr.preserved.values[h.host.values]
h = h[h.pres.isin(["f", "s", "d"])].copy(); h["pres"] = pd.Categorical(h.pres, categories=["d", "s", "f"]); h["abslat"] = np.abs(h.latitude)
md = smf.ols("res ~ C(pres) + abslat", data=h).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(h.filename)[0]})
for k in ["C(pres)[T.s]", "C(pres)[T.f]"]:
    b, se = md.params[k], md.bse[k]; R["depth_ratio_" + k[-2] + "_vs_d"] = [float(10 ** b), float(10 ** (b - 1.96 * se)), float(10 ** (b + 1.96 * se)), float(md.pvalues[k])]
R["depth_n"] = {k: dict(hollows=int(len(g)), craters=int(g.host.nunique())) for k, g in h.groupby("pres", observed=True)}
json.dump(R, open("host_age_record.json", "w"), indent=1); print(json.dumps(R, indent=1))
