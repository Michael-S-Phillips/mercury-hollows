import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Shadow-length depths on sloping ground. Blewett et al. (2016) depths assume level ground: d = L tan(e), e = solar elevation.
If the ground under the shadow tilts by beta along the shadow direction (positive when the slope faces the Sun),
the true depth is L (tan e + tan beta), so d_true = d_meas * f with f = 1 + tan(beta) / tan(e).
tan(beta) = tan(slope) * cos(facing azimuth - solar azimuth). Slope and facing azimuth from the USGS 665 m DEM
(depth_points_slope.csv, made on the Mac); solar geometry at each hollow from the image's sub-solar point (CDR index).
Tests: (1) does ln(measured depth) follow -ln f (coefficient -1 if the DEM slope is the slope under the shadow)?
(2) do the temperature, poleward, and incidence results change when depths are corrected?
(3) does corrected depth vary with regional slope? Normal stress alone predicts d proportional to 1/cos^2(slope).
Standard errors clustered by image. Writes depth_slope_record.json."""
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf
d = pd.read_csv(_REPO + "/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")
s = pd.read_csv("depth_points_slope.csv").set_index("row"); d = d.join(s)
ss = pd.read_csv("subsolar_depth_images.csv").set_index("key")
k = d.filename.str.slice(1, 13).str.upper(); d["sslat"] = k.map(ss.sslat); d["sslon"] = k.map(ss.sslon)
la, lo = np.radians(d.latitude), np.radians(d.longitude % 360); las, los = np.radians(d.sslat), np.radians(d.sslon % 360)
dl = los - lo
cosi = np.sin(la) * np.sin(las) + np.cos(la) * np.cos(las) * np.cos(dl)
d["inc_h"] = np.degrees(np.arccos(np.clip(cosi, -1, 1)))
d["sun_az"] = np.degrees(np.arctan2(np.sin(dl) * np.cos(las), np.cos(la) * np.sin(las) - np.sin(la) * np.cos(las) * np.cos(dl))) % 360
R = {"n": int(len(d)), "inc_h_minus_inc_center": [float(x) for x in (d.inc_h - d.inc).describe()[["mean", "std", "min", "max"]]]}
d["img"] = d.filename.astype("category").cat.codes
d["x"] = 1.0 / d["T"]; d["ltan"] = np.log(np.tan(np.radians(d.inc))); d["hi"] = (d.abslat >= 60).astype(int)
d["eq_face"] = np.where(d.latitude >= 0, np.cos(np.radians(d.face_az_s1 - 180)), np.cos(np.radians(d.face_az_s1)))
for bl in ("s1", "s3"):
    tb = np.tan(np.radians(d[f"slope_{bl}"])) * np.cos(np.radians(d[f"face_az_{bl}"] - d.sun_az))
    d[f"tanb_{bl}"] = tb; d[f"f_{bl}"] = 1 + tb / np.tan(np.radians(90 - d.inc_h))
R["slope_quantiles_s1"] = {str(q): float(d.slope_s1.quantile(q)) for q in (.25, .5, .75, .9, .99)}
R["f_quantiles_s1"] = {str(q): float(d.f_s1.quantile(q)) for q in (.01, .05, .25, .5, .75, .95, .99)}
ok = (d.f_s1 > 0.5) & (d.f_s3 > 0.5); R["n_f_gt_0.5"] = int(ok.sum())
d = d[ok].copy()
for bl in ("s1", "s3"): d[f"lnf_{bl}"] = np.log(d[f"f_{bl}"]); d[f"lnh_{bl}"] = d.lnh + d[f"lnf_{bl}"]
def fit(f, sub, name):
    m = smf.ols(f, data=sub).fit(cov_type="cluster", cov_kwds={"groups": sub.img})
    R[name] = {kk: [float(m.params[kk]), float(m.bse[kk]), float(m.pvalues[kk])] for kk in m.params.index}; R[name]["n"] = int(len(sub)); return m
# (1) attenuation diagnostic
fit("lnh ~ lnf_s1 + x + lnpx + ltan", d, "diag_lnf_s1"); fit("lnh ~ lnf_s3 + x + lnpx + ltan", d, "diag_lnf_s3")
# (2) headline results, uncorrected (same subset) and corrected
B = -9028.6679; dx = 1 / d["T"].min() - 1 / d["T"].max()
for y in ("lnh", "lnh_s1", "lnh_s3"):
    m = fit(f"{y} ~ x + lnpx + ltan", d, f"T_px_inc_{y}")
    sl, e = m.params["x"], m.bse["x"]; R[f"fold_{y}"] = [float(np.exp(-sl * dx)), float(np.exp(-(sl + 1.96 * e) * dx)), float(np.exp(-(sl - 1.96 * e) * dx))]
    m = fit(f"{y} ~ hi + lnpx + ltan", d, f"above60_{y}")
    b, e = m.params["hi"], m.bse["hi"]; R[f"above60_ratio_{y}"] = [float(np.exp(b)), float(np.exp(b - 1.96 * e)), float(np.exp(b + 1.96 * e))]
# (3) depth against regional slope, corrected depths; normal-stress-only prediction is a coefficient of +1 on ln(1/cos^2)
d["lsec2"] = -2 * np.log(np.cos(np.radians(d.slope_s1)))
fit("lnh_s1 ~ slope_s1 + x + lnpx + ltan", d, "slope_deg_s1"); fit("lnh_s1 ~ lsec2 + x + lnpx + ltan", d, "lsec2_s1")
fit("lnh_s1 ~ slope_s1 + eq_face + x + lnpx + ltan", d, "slope_and_aspect_s1")
d["sbin"] = pd.cut(d.slope_s1, [0, 2, 4, 6, 10, 90], labels=["a0_2", "b2_4", "c4_6", "d6_10", "e10plus"], include_lowest=True).astype(str)
m = smf.ols("lnh_s1 ~ C(sbin) + x + lnpx + ltan", data=d).fit(cov_type="cluster", cov_kwds={"groups": d.img})
R["slope_bins_s1"] = {kk: [float(np.exp(m.params[kk])), float(np.exp(m.params[kk] - 1.96 * m.bse[kk])), float(np.exp(m.params[kk] + 1.96 * m.bse[kk]))] for kk in m.params.index if "sbin" in kk}
R["slope_bin_counts"] = d.sbin.value_counts().to_dict()
json.dump(R, open("depth_slope_record.json", "w"), indent=1)
print("inc_h - inc_center", [round(v, 2) for v in R["inc_h_minus_inc_center"]])
print("f quantiles", {k: round(v, 2) for k, v in R["f_quantiles_s1"].items()}, "kept", R["n_f_gt_0.5"])
for nm in ("diag_lnf_s1", "diag_lnf_s3"): print(nm, [round(v, 3) for v in R[nm][nm.split("_", 1)[1].replace("diag_", "") if False else ("lnf_s1" if "s1" in nm else "lnf_s3")]])
for y in ("lnh", "lnh_s1", "lnh_s3"):
    print(y, "T slope", [round(v, 1) for v in R[f"T_px_inc_{y}"]["x"][:2]], "fold", [round(v, 2) for v in R[f"fold_{y}"]], "above60", [round(v, 2) for v in R[f"above60_ratio_{y}"]], "ltan", [round(v, 3) for v in R[f"T_px_inc_{y}"]["ltan"][:2]])
print("slope_deg", [round(v, 4) for v in R["slope_deg_s1"]["slope_s1"]], "lsec2", [round(v, 2) for v in R["lsec2_s1"]["lsec2"]])
print("slope+aspect", {k: [round(v, 4) for v in R["slope_and_aspect_s1"][k]] for k in ("slope_s1", "eq_face")})
print("bins", {k[-9:]: [round(v, 2) for v in vv] for k, vv in R["slope_bins_s1"].items()}, R["slope_bin_counts"])
