#!/usr/bin/env python3
"""Add-on to stage_test.py, 2 Nov 2026: field-level check of the Stage 1 step near Delta T_w = 44-50 K. One row per hollow
field (labels from images with incidence < 70 deg, |lat| <= 60): share of Stage 1 among labeled hollows, field-median
incidence. Weighted GLM of the share (weights = labeled hollows), and an unweighted logistic model of 'all labeled hollows
Stage 1'. Lists the fields that carry most hollows above 50 K. Writes stage_fieldlevel.json."""
import json, numpy as np, pandas as pd
src = open(__file__.replace("stage_fieldlevel.py", "stage_test.py")).read().split("OUT[\"results\"] = R")[0]
exec(src)
d = S[(S.abslat <= 60) & (S.inc8_m < 70)]
F = d.groupby("field_id").agg(n=("s1", "size"), s1=("s1", "mean"), dT=("dT", "median"), inc10=("inc10", "median"), logpx=("logpx", "median"),
                              emi=("emi", "median"), abslat=("abslat", "median"), lat=("lat", "median"), lon=("lon", "median")).reset_index()
F["hi"] = (F.dT >= 50).astype(int); F["all1"] = (F.s1 == 1).astype(int)
out = {"n_fields": int(len(F)), "n_fields_dT_ge50": int(F.hi.sum())}
for t in (44, 50):
    F["st"] = (F.dT >= t).astype(float)
    m = smf.glm("s1 ~ st + inc10 + logpx + emi + abslat", data=F, family=sm.families.Binomial(), var_weights=F.n).fit(scale="X2")
    m2 = smf.glm("all1 ~ st + inc10 + logpx + emi + abslat", data=F, family=sm.families.Binomial()).fit()
    out[f"t{t}"] = dict(share_model_OR=float(np.exp(m.params["st"])), share_ci=[float(np.exp(m.params["st"] - 1.96 * m.bse["st"])), float(np.exp(m.params["st"] + 1.96 * m.bse["st"]))],
                        share_p=float(m.pvalues["st"]), all_s1_OR=float(np.exp(m2.params["st"])),
                        all_s1_ci=[float(np.exp(m2.params["st"] - 1.96 * m2.bse["st"])), float(np.exp(m2.params["st"] + 1.96 * m2.bse["st"]))], all_s1_p=float(m2.pvalues["st"]),
                        frac_fields_all_s1_above=float(F[F.st == 1].all1.mean()), frac_fields_all_s1_below=float(F[F.st == 0].all1.mean()),
                        mean_share_above=float(F[F.st == 1].s1.mean()), mean_share_below=float(F[F.st == 0].s1.mean()),
                        n_above=int((F.st == 1).sum()), n_below=int((F.st == 0).sum()))
    print(t, out[f"t{t}"])
top = F[F.hi == 1].sort_values("n", ascending=False).head(10)
out["top_fields_dT_ge50"] = top.round(3).to_dict("records"); print(top.round(2).to_string())
json.dump(out, open(f"{HERE}/stage_fieldlevel.json", "w"), indent=1, default=float)
