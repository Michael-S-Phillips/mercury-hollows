#!/usr/bin/env python3
"""g24 step 2 (9 Nov): Michael's blind hand depths against Blewett et al. (2016) on the same points (blewett_cmp_2026-11-09.json).
Point level: Michael's line for Blewett point k of box fid (marks_key.csv) against that point's Blewett depth. Box level: deepest of
Michael's lines against Blewett's deepest point. Ratios Blewett/Michael for depth and for shadow length (Blewett shadow length taken as
depth x tan(their incidence); Michael's as L_px x px x f), poleward (>= 60N) and equatorward, with 95% intervals from 2,000 bootstrap
resamples of frames. Then the paper's question on the same points: ln depth ~ poleward + ln px + ln tan i for Blewett's depths and for
Michael's, SEs clustered by frame. Writes cmp_summary.json and cmp_points.csv."""
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings; warnings.filterwarnings("ignore")
docs = {int(d["fid"]): d for d in json.load(open("blewett_cmp_2026-11-09.json"))}
man = json.load(open("manifest_cmp.json")); mk = pd.read_csv("marks_key.csv"); key = pd.read_csv("key_cmp.csv")
rows = []
for _, r in mk.iterrows():
    d = docs.get(int(r.fid), {}); l = (d.get("lines") or {}).get(str(int(r.k)))
    m = man[str(int(r.fid))]
    rows.append(dict(fid=int(r.fid), k=int(r.k), bid=int(r.bid), num=m["img"], b_depth=r.b_depth, b_inc=r.b_inc, b_px=r.b_px, lat=r.lat,
                     m_depth=l.get("depth_m") if l else np.nan, m_Lpx=l.get("L_px") if l else np.nan, skip=bool(l and l.get("skip")),
                     px=m["px"], f=m["f"], inc=m["inc"], cls=d.get("cls"), shadow=d.get("shadow")))
p = pd.DataFrame(rows); p["pole"] = (p.lat >= 60).astype(int)
p["m_L"] = p.m_Lpx * p.px * p.f; p["b_L"] = p.b_depth * np.tan(np.radians(p.b_inc))
p["r_depth"] = p.b_depth / p.m_depth; p["r_L"] = p.b_L / p.m_L; p["d_inc"] = p.b_inc - p.inc
p.to_csv("cmp_points.csv", index=False)
ok = p[p.m_depth > 0].copy()
out = dict(points=len(p), measured=int(len(ok)), skipped=int(p.skip.sum()), missing=int((p.m_depth.isna() & ~p.skip).sum()),
           boxes_all_points=int(sum(all((str(int(k)) in (docs[f].get("lines") or {})) for k in mk[mk.fid == f].k) for f in docs)))
rng = np.random.default_rng(9)
def boot_med(df, col, nb=2000):
    fr = df.num.unique(); g = {f: df[df.num == f][col].values for f in fr}; meds = []
    for _ in range(nb):
        s = rng.choice(fr, len(fr)); v = np.concatenate([g[f] for f in s]); meds.append(np.median(v))
    return [round(float(np.median(df[col])), 3), round(float(np.percentile(meds, 2.5)), 3), round(float(np.percentile(meds, 97.5)), 3)]
for lvl, d in [("point", ok), ("box", None)]:
    if d is None:
        g = ok.groupby("fid").agg(num=("num", "first"), pole=("pole", "first"), m_depth=("m_depth", "max"), b_depth=("b_depth", "max"), px=("px", "first"), inc=("inc", "first"), b_inc=("b_inc", "first"))
        d = g.reset_index(); d["r_depth"] = d.b_depth / d.m_depth; d["b_L"] = d.b_depth * np.tan(np.radians(d.b_inc)); d["m_L"] = d.m_depth * np.tan(np.radians(d.inc)); d["r_L"] = d.b_L / d.m_L
    for grp, s in [("equatorward", d[d.pole == 0]), ("poleward", d[d.pole == 1]), ("all", d)]:
        out[f"{lvl}_{grp}"] = dict(n=len(s), frames=int(s.num.nunique()), ratio_depth_B_over_M=boot_med(s, "r_depth"), ratio_shadow_B_over_M=boot_med(s, "r_L"),
                                   med_B=round(float(s.b_depth.median()), 1), med_M=round(float(s.m_depth.median()), 1),
                                   spearman=round(float(s[["b_depth", "m_depth"]].corr("spearman").iloc[0, 1]), 3))
        print(lvl, grp, out[f"{lvl}_{grp}"])
# does the ratio differ poleward? ln(B/M) ~ pole, clustered by frame
for lvl, d in [("point", ok)]:
    d = d.assign(lr=np.log(d.r_depth)); m = smf.ols("lr ~ pole", d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.num)[0]})
    out["pole_effect_on_B_over_M"] = [round(float(np.exp(m.params.pole)), 3), round(float(np.exp(m.params.pole - 1.96 * m.bse.pole)), 3), round(float(np.exp(m.params.pole + 1.96 * m.bse.pole)), 3), round(float(m.pvalues.pole), 4)]
    print("poleward factor on B/M", out["pole_effect_on_B_over_M"])
# the paper's poleward model on the same points, each measurement set
s = ok.copy(); s["lnpx"] = np.log(s.px); s["ltan"] = np.log(np.tan(np.radians(s.inc))); s["ltanb"] = np.log(np.tan(np.radians(s.b_inc))); s["lnpxb"] = np.log(s.b_px)
for nm, f in [("blewett", "np.log(b_depth) ~ pole + lnpxb + ltanb"), ("michael", "np.log(m_depth) ~ pole + lnpx + ltan")]:
    m = smf.ols(f, s).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(s.num)[0]})
    out[f"poleward_ratio_{nm}_same_points"] = [round(float(np.exp(m.params.pole)), 3), round(float(np.exp(m.params.pole - 1.96 * m.bse.pole)), 3), round(float(np.exp(m.params.pole + 1.96 * m.bse.pole)), 3)]
    print("poleward ratio on same points,", nm, out[f"poleward_ratio_{nm}_same_points"])
out["incidence_diff_B_minus_ours"] = dict(median=round(float(ok.d_inc.median()), 2), q10=round(float(ok.d_inc.quantile(.1)), 2), q90=round(float(ok.d_inc.quantile(.9)), 2))
out["auto_pilot_box"] = key.merge(ok.groupby("fid").m_depth.max().rename("m_max"), on="fid").assign(r=lambda x: x.auto_depth_pilot / x.m_max).groupby("pole").r.median().round(3).to_dict()
print("incidence diff", out["incidence_diff_B_minus_ours"], "| auto/M by pole (box)", out["auto_pilot_box"]); print({k: out[k] for k in ["points", "measured", "skipped", "missing", "boxes_all_points"]})
json.dump(out, open("cmp_summary.json", "w"), indent=1, default=str)
