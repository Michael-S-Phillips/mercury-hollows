#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Test 1: does hollow depth track mean annual temperature as the diffusion law predicts?

Diffusion-limited relief is h = sqrt(2 mu D Pv(T) f t / (R T rho)), and Pv = exp(A + B/T)
dominates, so ln h = const + B/(2T): an Arrhenius line with slope B/2 = -4514 K for sulfur.
A layer-limited depth instead gives slope 0 with a rolloff only at the cold margin.
"""
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from thermal import build_grid, lookup
import os

UP = _REPO + ""
d = pd.read_csv(_REPO + "/manuscripts/GRL_insolation_hollows/data/hollows_blewett2016_locations.csv")
d.columns = [c.strip() for c in d.columns]
d = d.rename(columns={"pixel scale (m)": "px", "depth (m)": "depth"})
d = d[["longitude", "latitude", "px", "depth", "filename"]].dropna()
d = d[(d.depth > 0) & (d.px > 0)]

g = np.load(_REPO + "/manuscripts/GRL_insolation_hollows/data/thermal_grid.npz")
d["T"] = lookup(g["la"], g["lo"], g["G"], d.latitude.values, d.longitude.values)
d["invT"] = 1.0 / d["T"]
d["lnh"] = np.log(d.depth)
d["lnpx"] = np.log(d.px)
d["abslat"] = d.latitude.abs()

print(f"n = {len(d)} hollows with shadow-length depths")
print(f"depth   median {d.depth.median():.1f} m, 5-95% {d.depth.quantile(.05):.1f}-{d.depth.quantile(.95):.1f}")
print(f"mean annual T spanned: {d['T'].min():.0f} - {d['T'].max():.0f} K")
print(f"pixel scale: median {d.px.median():.1f} m, range {d.px.min():.1f}-{d.px.max():.1f}\n")

B_S = -9028.6679
print("PREDICTION (sulfur, diffusion-limited): d ln(h) / d(1/T) = B/2 = "
      f"{B_S/2:.0f} K  ->  over the observed T range, depth should change by a factor of "
      f"{np.exp(B_S/2*(d.invT.max()-d.invT.min())):.3g}\n")

print("=" * 78)
print("A. Arrhenius regression, ln(depth) on 1/T")
for lab, f in [("raw", "lnh ~ invT"),
               ("+ pixel scale", "lnh ~ invT + lnpx"),
               ("+ pixel scale + |lat|", "lnh ~ invT + lnpx + abslat")]:
    m = smf.ols(f, data=d).fit()
    b, se = m.params["invT"], m.bse["invT"]
    print(f"  {lab:24s} slope {b:9.0f} +/- {se:6.0f} K   p={m.pvalues['invT']:.2g}   R2={m.rsquared:.3f}")
print(f"  predicted if diffusion-limited: {B_S/2:9.0f} K")
print(f"  predicted if layer-limited:             0 K")

print("=" * 78)
print("B. The clean design: longitude only, within latitude bands.")
print("   Pixel scale and viewing geometry track latitude, not longitude, so the")
print("   hot-pole/warm-pole contrast at fixed |lat| is the uncontaminated test.")
d["band"] = pd.cut(d.abslat, [0, 20, 40, 55, 70, 90])
for bnd, sub in d.groupby("band", observed=True):
    if len(sub) < 60:
        continue
    m = smf.ols("lnh ~ invT + lnpx", data=sub).fit()
    tr = sub["T"].max() - sub["T"].min()
    print(f"  |lat| {str(bnd):10s} n={len(sub):5d}  T range {tr:5.1f} K  "
          f"slope {m.params['invT']:9.0f} +/- {m.bse['invT']:6.0f}  p={m.pvalues['invT']:.2g}")

print("=" * 78)
print("C. Binned medians, pixel-scale corrected (residual depth at the median pixel scale)")
m0 = smf.ols("lnh ~ lnpx", data=d).fit()
d["lnh_corr"] = m0.resid + m0.params["Intercept"] + m0.params["lnpx"] * d.lnpx.median()
d["Tbin"] = pd.cut(d["T"], [180, 240, 270, 290, 305, 315, 325, 335])
tb = d.groupby("Tbin", observed=True).agg(n=("depth", "size"),
        T=("T", "median"), raw=("depth", "median"),
        corr=("lnh_corr", lambda s: float(np.exp(np.median(s)))))
print(tb.to_string(float_format=lambda v: f"{v:8.1f}"))
lo, hi = tb.iloc[0], tb.iloc[-1]
print(f"\n  observed corrected depth ratio, coldest to hottest bin: "
      f"{hi['corr']/lo['corr']:.2f}x   over {hi['T']-lo['T']:.0f} K")
print(f"  diffusion law predicts over that same interval:        "
      f"{np.exp(B_S/2*(1/hi['T']-1/lo['T'])):.3g}x")
d.to_csv(_REPO + "/manuscripts/GRL_insolation_hollows/data/test1_rows.csv", index=False)
