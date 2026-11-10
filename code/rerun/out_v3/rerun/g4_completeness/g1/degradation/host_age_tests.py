# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/degradation/host_age_tests.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
"""Host-age tests of the 'early burst + long-term activity' idea (27 Oct 2026, exploratory).
Herrick et al. (2018) preservation codes in column 'preserved': f, s, d (read here as fresh, intermediate,
degraded; confirm against the catalog documentation). Rays from column 'ejecta' containing 'r'.
(1) odds of hosting hollows by preservation class, at fixed log diameter and |latitude| (craters >= 10 km);
(2) warm-longitude preference (odds per 10 K of dT) within each class, and the dT x class interaction;
(3) corrected hollow depth (Blewett et al. 2016; pixel scale and incidence controlled) by host class,
    standard errors clustered by image."""
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
# g1 rerun (6 Nov 2026): copy of analysis_v3/host_age_tests.py; only input/output paths changed; the depth block (3) needs the cloud copy of depth_structure.py and runs only with RUN_DEPTH=1 (depths do not depend on the hollow catalog).
import os as _os, sys as _sys; _G1 = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))); _sys.path.insert(0, _G1)
from hostfiles import HOSTS, CRATERS, CAT; OUTDIR = f"{_G1}/degradation/{CAT}"; _os.makedirs(OUTDIR, exist_ok=True)
FEATURES = f"{_G1}/summary_model/{CAT}/crater_features.csv"; KDIR = _os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/kinczyk_degradation")
G = _os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/GRL_insolation_hollows")
sys.path.insert(0, G + "/scripts"); import thermal
z = np.load(G + "/data/thermal_grid.npz"); la, lo, T = z["la"], z["lo"], z["G"]
R = {}
s = pd.read_csv(CRATERS)
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
if _os.environ.get("RUN_DEPTH") == "1":
  exec(open("/home/claude/grl/depthctl/depth_structure.py").read().split("R = {")[0])
  h = d[d.host >= 0].copy(); h["pres"] = cr.preserved.values[h.host.values]
  h = h[h.pres.isin(["f", "s", "d"])].copy(); h["pres"] = pd.Categorical(h.pres, categories=["d", "s", "f"]); h["abslat"] = np.abs(h.latitude)
  md = smf.ols("res ~ C(pres) + abslat", data=h).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(h.filename)[0]})
  for k in ["C(pres)[T.s]", "C(pres)[T.f]"]:
    b, se = md.params[k], md.bse[k]; R["depth_ratio_" + k[-2] + "_vs_d"] = [float(10 ** b), float(10 ** (b - 1.96 * se)), float(10 ** (b + 1.96 * se)), float(md.pvalues[k])]
  R["depth_n"] = {k: dict(hollows=int(len(g)), craters=int(g.host.nunique())) for k, g in h.groupby("pres", observed=True)}
json.dump(R, open(f"{OUTDIR}/host_age_record.json", "w"), indent=1); print(json.dumps(R, indent=1))
