#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Depth of hollows (Blewett et al. 2016 shadow-length depths) against mean annual temperature, latitude, and image geometry.
v2 corrections to the v1 analysis:
  (1) v1 held absolute latitude fixed, which removes most of the temperature variation; here temperature is not paired with latitude.
  (2) solar incidence of each image is added as a control. Shadow-length depths depend on incidence because the shortest
      measurable shadow sets a minimum measurable depth of about (shadow pixels x pixel scale) / tan(incidence).
Incidence comes from the MDIS CDR cumulative index (MSGRMDS_2001, CUMINDEX.TAB), joined on the product number (EN... -> CN...).
Writes data/test1_rows_inc.csv and records/depth_reanalysis_record.json. Standard errors are clustered by image."""
import json, os, pickle, numpy as np, pandas as pd, statsmodels.formula.api as smf
G = os.environ.get("INSOL_DIR", _REPO + "/manuscripts/insolation_thermal")
W = os.environ.get("INSOL_WORK", _REPO + "/manuscripts/insolation_thermal")
IDX = os.environ.get("CDR_INDEX", _REPO + "/external/mdis_index/cdr_index.pkl")
os.makedirs(f"{W}/records", exist_ok=True); os.makedirs(f"{W}/data", exist_ok=True)
d = pd.read_csv(f"{G}/data/test1_rows.csv")
x = pickle.load(open(IDX, "rb")); x["key"] = x.file.str.slice(1, 13).str.upper()
xi = x.drop_duplicates("key").set_index("key")
d["inc"] = d.filename.str.slice(1, 13).map(xi.inc); assert d.inc.notna().all()
d.to_csv(f"{W}/data/test1_rows_inc.csv", index=False)
d["img"] = d.filename.astype("category").cat.codes
d["x"] = 1.0 / d["T"]; d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["hi"] = (d.abslat >= 60).astype(int)
B = -9028.6679
R = {"n": int(len(d)), "n_img": int(d.img.nunique()), "inc_median": float(d.inc.median()), "px_median": float(d.px.median())}
def fit(f, sub, name):
    m = smf.ols(f, data=sub).fit(cov_type="cluster", cov_kwds={"groups": sub.img})
    R[name] = {k: [float(m.params[k]), float(m.bse[k]), float(m.pvalues[k])] for k in m.params.index}
    R[name]["n"] = int(len(sub)); R[name]["n_img"] = int(sub.img.nunique()); return m
fit("lnh ~ x + lnpx + abslat", d, "v1_with_lat")
m_px = fit("lnh ~ x + lnpx", d, "T_px")
m_pi = fit("lnh ~ x + lnpx + ltan", d, "T_px_inc")
lo = d[d.abslat < 60]; fit("lnh ~ x + lnpx", lo, "T_px_below60"); fit("lnh ~ x + lnpx + ltan", lo, "T_px_inc_below60")
s2 = d[(d.inc >= 60) & (d.inc < 75)]; fit("lnh ~ x + lnpx", s2, "T_px_inc60_75only")
m_hi = fit("lnh ~ hi + lnpx + ltan", d, "above60_px_inc"); fit("lnh ~ hi + lnpx", d, "above60_px")
fit("lnh ~ abslat + lnpx", d, "lat_px"); fit("lnh ~ abslat + lnpx + ltan", d, "lat_px_inc")
t0, t1 = d["T"].min(), d["T"].max(); dx = 1 / t0 - 1 / t1
R["T_range"] = [float(t0), float(t1)]; R["below60_T_range"] = [float(lo["T"].min()), float(lo["T"].max())]
for nm, m in [("T_px", m_px), ("T_px_inc", m_pi)]:
    s, e = m.params["x"], m.bse["x"]
    R[f"fold_{nm}"] = [float(np.exp(-s * dx)), float(np.exp(-(s + 1.96 * e) * dx)), float(np.exp(-(s - 1.96 * e) * dx))]
R["fold_pred_sulfur_front"] = float(np.exp(-B / 2 * dx))
b, e = m_hi.params["hi"], m_hi.bse["hi"]
R["above60_ratio_px_inc"] = [float(np.exp(b)), float(np.exp(b - 1.96 * e)), float(np.exp(b + 1.96 * e))]
# median depth by latitude band, corrected to median pixel scale and incidence; 95% intervals by bootstrap over images
c = smf.ols("lnh ~ lnpx + ltan", data=d).fit(); c0 = smf.ols("lnh ~ lnpx", data=d).fit()
d["hc"] = np.exp(d.lnh - c.params["lnpx"] * (d.lnpx - d.lnpx.median()) - c.params["ltan"] * (d.ltan - d.ltan.median()))
d["hp"] = np.exp(d.lnh - c0.params["lnpx"] * (d.lnpx - d.lnpx.median()))
rng = np.random.default_rng(20261025); R["bands"] = {}
for a, bb in [(0, 30), (30, 45), (45, 60), (60, 75), (75, 90)]:
    s_ = d[(d.abslat >= a) & (d.abslat < bb)]; imgs = s_.img.unique()
    out = dict(n=int(len(s_)), n_img=int(len(imgs)), T_med=float(s_["T"].median()), inc_med=float(s_.inc.median()), lat_med=float(s_.abslat.median()))
    for col in ("hc", "hp"):
        g = {i: s_[col].values[s_.img.values == i] for i in imgs}
        bt = [np.median(np.concatenate([g[i] for i in rng.choice(imgs, len(imgs))])) for _ in range(2000)]
        out[col] = [float(s_[col].median()), float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))]
    R["bands"][f"{a}-{bb}"] = out
R["ltan_coef_note"] = "ln depth falls with ln tan(incidence); consistent with a minimum measurable shadow length"
json.dump(R, open(f"{W}/records/depth_reanalysis_record.json", "w"), indent=1)
for k in ("v1_with_lat", "T_px", "T_px_inc", "T_px_below60", "T_px_inc_below60", "T_px_inc60_75only"):
    print(k, [round(v, 1) for v in R[k]["x"][:2]], "p=%.3f" % R[k]["x"][2], "n", R[k]["n"])
print("fold", {k: [round(v, 2) for v in R[k]] for k in ("fold_T_px", "fold_T_px_inc")}, "pred", round(R["fold_pred_sulfur_front"]))
print("above60 ratio", [round(v, 2) for v in R["above60_ratio_px_inc"]], "ltan", [round(v, 3) for v in R["T_px_inc"]["ltan"]])
for k, v in R["bands"].items(): print(k, v["n"], v["n_img"], "inc", round(v["inc_med"]), "hc", [round(z, 1) for z in v["hc"]], "hp", [round(z, 1) for z in v["hp"]])
