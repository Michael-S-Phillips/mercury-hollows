#!/usr/bin/env python3
"""g25 (11 Nov 2026): table of the fitted coefficients of the thinned occurrence model with all seven groups (supplement v29,
Section S16, Michael: the full equation belongs in the supplement). Same design and fit as comp_model.py; intervals are the
2.5th and 97.5th percentiles over the same 100 block resamples (seed 21) that give the rate ratios of the composition terms in
comp.json, so those three intervals are reproduced here as a check. Also checks that the coefficients predict a total of 439.9
host craters (S16, Location of the clustering). Run:
  SM_DIR=.../g4_completeness/g1/summary_model/v3r python3 coef_table.py
Writes coef_table.csv and coef_table.json here."""
import os, sys, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
SMD = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g4_completeness/g1/summary_model"); sys.path.insert(0, SMD)
import step7_improved as s7
col = pd.read_csv(os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g6_remaining/color/out_v3/crater_color.csv"), dtype={"id": str}, usecols=["id", "ok", "pc1", "pc2"])
d = s7.d.copy(); d["id"] = d.id.astype(str); d = d.merge(col, on="id", how="left")
ok = d.ok.fillna(False).astype(bool).values; d["nocolor"] = (~ok).astype(int)
for k in ["pc1", "pc2"]:
    m, s = d.loc[ok, k].mean(), d.loc[ok, k].std(); d[k + "_z"] = np.where(ok, (d[k] - m) / s, 0.0)
COMP = ["lrm", "pc1_z", "pc2_z", "nocolor"]; GR = {**s7.PG, "composition": COMP}
y, off, fit, _, _ = s7.make(d, GR, "cloglog")
cols = [c for g in GR for c in GR[g]]; names = ["Intercept", "lnC"] + cols; grp = {c: g for g in GR for c in GR[g]}
b, ll = fit(cols, np.ones(len(d)))
X = np.column_stack([np.ones(len(d)), d.lnC.values] + [d[c].values.astype(float) for c in cols])
P = -np.expm1(-np.exp(X @ b + off))
rng = np.random.default_rng(21); blocks = d.block.unique(); BB = []
for _ in range(100):
    cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); w = d.block.map(cnt).fillna(0).values.astype(float); BB.append(fit(cols, w)[0])
BB = np.array(BB); lo, hi = np.percentile(BB, 2.5, axis=0), np.percentile(BB, 97.5, axis=0)
rows = [dict(term=n, group=grp.get(n, ""), coef=float(b[i]), lo=float(lo[i]), hi=float(hi[i]), factor=float(np.exp(b[i])), factor_lo=float(np.exp(lo[i])),
             factor_hi=float(np.exp(hi[i])), n_ones=int((d[n] == 1).sum()) if n in cols and set(np.unique(d[n])) <= {0, 1} else None,
             hosts_ones=int(((d[n] == 1) & (y == 1)).sum()) if n in cols and set(np.unique(d[n])) <= {0, 1} else None,
             boot_abs_gt_10=int((np.abs(BB[:, i]) > 10).sum())) for i, n in enumerate(names)]
pd.DataFrame(rows).to_csv(f"{HERE}/coef_table.csv", index=False)
cm = json.load(open(f"{HERE}/comp.json"))
check = {c: dict(here=[rows[names.index(c)][k] for k in ("factor", "factor_lo", "factor_hi")], comp_json=cm["rate_ratios"][c]) for c in COMP}
cen = dict(logD_mean=float(d.logD.mean()), abslat_mean=float(d.abslat.mean()), ln_zmax_mean=float(d.ln_zmax.mean()), dT10_mean=float(d.dT10.mean()))
R = dict(created="2026-11-11", n=int(len(d)), events=int(y.sum()), loglik=float(ll), sum_P=float(P.sum()), centering=cen, check_vs_comp_json=check,
         n_nocolor=int(d.nocolor.sum()), nocolor_hosts=int(d[d.nocolor == 1].h09.sum()), coefs=rows)
json.dump(R, open(f"{HERE}/coef_table.json", "w"), indent=1)
print(json.dumps({k: R[k] for k in ["n", "events", "sum_P", "centering", "check_vs_comp_json", "n_nocolor", "nocolor_hosts"]}, indent=1))
print(pd.DataFrame(rows).round(4).to_string())
