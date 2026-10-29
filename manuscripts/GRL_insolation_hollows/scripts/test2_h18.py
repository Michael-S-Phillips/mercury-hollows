#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Test 2 rebuilt on Herrick et al. (2018), the catalog the hollow host assignments were made in.
Hosts are flagged by Herrick unique_id (exact), not by position matching into another catalog."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf, json, os
from thermal import lookup
g = np.load(_REPO + "/manuscripts/GRL_insolation_hollows/data/thermal_grid.npz")
T = lambda la, lo: lookup(g["la"], g["lo"], g["G"], la, lo)
H = os.path.expanduser(_REPO + "/data/ME6M300TGT/herrick2018_all.csv")
HA = os.path.expanduser(_REPO + "/code/rerun/latitude/aspect/host_assignments.csv")
cr = pd.read_csv(H).rename(columns={"diam": "diameter", "lat": "latitude"})
cr["longitude"] = ((cr.lon + 180) % 360) - 180
cr["lon360"] = cr.lon % 360
cr["Tmean"] = T(cr.latitude.values, cr.lon360.values)
cr["logD"] = np.log10(cr.diameter); cr["abslat"] = cr.latitude.abs()
ha = pd.read_csv(HA)
hid = set(ha.host_uid.dropna().astype(int))
cr["hosts"] = cr.unique_id.astype(int).isin(hid).astype(int)
print(f"Herrick 2018: {len(cr)} craters, lon range {cr.lon.min():.1f}..{cr.lon.max():.1f}; hosts flagged {cr.hosts.sum()} of {len(hid)}")
sub = cr[cr.diameter >= 5].copy()
sub.to_csv(_REPO + "/manuscripts/GRL_insolation_hollows/data/test2_craters_h18.csv", index=False)
R = {"n_craters": int(len(sub)), "n_hosts": int(sub.hosts.sum())}
m = smf.logit("hosts ~ Tmean + logD + abslat", data=sub).fit(disp=0)
R["logD_all"] = [m.params["logD"], m.bse["logD"]]
print(f"A. all latitudes, hosts ~ Tmean + logD + abslat: logD {m.params['logD']:.2f} ± {m.bse['logD']:.2f} ({sub.hosts.sum()} hosts)")
s = sub[sub.abslat <= 60].copy()
lon = s.longitude.values % 180.0
s["hotness"] = np.abs(np.where(lon > 90, 180 - lon, lon) - 90)
print(f"\n|lat|<=60: {len(s)} craters, {int(s.hosts.sum())} hosts")
R["n60"] = int(len(s)); R["h60"] = int(s.hosts.sum())
print("B. cutoff sweep (hot-longitude coef, log-odds per degree)")
R["sweep"] = {}
for cut in [15, 20, 30, 40, 60, 80, 100, 150]:
    s["big"] = (s.diameter >= cut).astype(int); out = []
    for gg in (0, 1):
        d2 = s[s.big == gg]; mm = smf.logit("hosts ~ hotness + logD + abslat", data=d2).fit(disp=0)
        out.append((int(d2.hosts.sum()), mm.params["hotness"], mm.bse["hotness"], mm.pvalues["hotness"]))
    mi = smf.logit("hosts ~ hotness * big + logD + abslat", data=s).fit(disp=0)
    (h0, c0, e0, p0), (h1, c1, e1, p1) = out
    R["sweep"][cut] = dict(h_small=h0, c_small=c0, e_small=e0, p_small=p0, h_big=h1, c_big=c1, e_big=e1, p_big=p1,
                           inter=mi.params["hotness:big"], inter_se=mi.bse["hotness:big"], inter_p=mi.pvalues["hotness:big"])
    print(f"  {cut:4d} km: small {c0:+.4f}±{e0:.4f} p={p0:.1g} ({h0}) | big {c1:+.4f}±{e1:.4f} p={p1:.2g} ({h1}) | inter {mi.params['hotness:big']:+.4f}±{mi.bse['hotness:big']:.4f} p={mi.pvalues['hotness:big']:.2g}")
m = smf.logit("hosts ~ hotness * logD + abslat", data=s).fit(disp=0)
b0, b1 = m.params["hotness"], m.params["hotness:logD"]
R["continuous"] = dict(hot=b0, hot_se=m.bse["hotness"], inter=b1, inter_se=m.bse["hotness:logD"], inter_p=m.pvalues["hotness:logD"], zero_km=10**(-b0/b1))
print(f"C. continuous: hotness {b0:+.4f}±{m.bse['hotness']:.4f}; hotness:logD {b1:+.4f}±{m.bse['hotness:logD']:.4f} p={m.pvalues['hotness:logD']:.2g}; zero at {10**(-b0/b1):.0f} km")
print("D. size bins")
R["bins"] = {}
bins = [5, 10, 20, 40, 80, 160, 2000]
for lo_, hi in zip(bins[:-1], bins[1:]):
    d2 = s[(s.diameter >= lo_) & (s.diameter < hi)]
    mm = smf.logit("hosts ~ hotness + abslat", data=d2).fit(disp=0)
    R["bins"][f"{lo_}-{hi}"] = dict(n=len(d2), h=int(d2.hosts.sum()), c=mm.params["hotness"], e=mm.bse["hotness"], p=mm.pvalues["hotness"])
    print(f"  {lo_:4d}-{hi:<4d} n={len(d2):5d} hosts={int(d2.hosts.sum()):3d} coef {mm.params['hotness']:+.4f}±{mm.bse['hotness']:.4f} p={mm.pvalues['hotness']:.2g}")
json.dump(R, open(_REPO + "/manuscripts/GRL_insolation_hollows/records/test2_h18_record.json", "w"), indent=1, default=float)
