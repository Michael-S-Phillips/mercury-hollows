# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/step8_figure_improved_v2.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Summary model, step 8: figure for the improved models (step 7). (a) Share of the explained variation in hosting by
physical group: thinned-occurrence model (detection as thinning, area as exposure; 95% block-bootstrap intervals), the
logistic model on well-imaged craters (C >= 0.2; detection group dropped and the rest rescaled to 100%), and the
negative-binomial count model. (b) Rate ratios from the thinned-occurrence model on the XRS subset, per standard
deviation of each element ratio, with 95% block-bootstrap intervals (200 draws of the full model).
Writes out/fig_improved.png and out/xrs_rate_ratios.csv."""
# v2 (5 Nov 2026): group labels renamed to match the main text (degradation state, crater size beyond area, vents, low-reflectance material); writes figures/figS8_joint_model_v2.png.
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/step8_figure_improved_v2.py; only input/output paths changed. SM_DIR (env) is the
# folder holding this catalog's crater_features*.csv and out/ (g1_graph_joint/summary_model/<CATALOG>).
import json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import step7_improved as s
from common import OUT
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"; S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
R = json.load(open(f"{OUT}/improved.json"))
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.facecolor": SURF, "figure.facecolor": SURF})
NAME = {"size": "Crater size (beyond area)", "insolation": "Insolation (latitude, ΔT$_w$)", "age": "Degradation state (fresh → filled)",
        "nesting": "Nesting (parent, position, depth of origin)", "cover": "Volcanic cover (plains, Caloris, fill)", "volcanism": "Vents",
        "composition": "Low-reflectance material"}
T = R["thinned_presence_all"]; G = sorted(T["point"]["share"], key=lambda g: T["point"]["share"][g])
W = R["logistic_C_ge_0.2"]["share"]; Wn = {g: W[g] / sum(W[k] for k in W if k != "detection") for g in W if g != "detection"}
Cn = R["thinned_counts_all"]["point"]["share"]
fig, ax = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.6})
a = ax[0]; yy = np.arange(len(G)); h = 0.26
v = np.array([T["point"]["share"][g] for g in G]); lo = np.array([T["ci"][g][0] for g in G]); hi = np.array([T["ci"][g][1] for g in G])
a.barh(yy + h, 100 * v, h, color=S1, label=f"Thinned occurrence model, presence (all {T['n']:,} craters)", zorder=3)
a.errorbar(100 * v, yy + h, xerr=[100 * (v - lo), 100 * (hi - v)], fmt="none", ecolor=INK2, lw=0.8, capsize=2, zorder=4)
a.barh(yy, [100 * Wn[g] for g in G], h, color=S2, label=f"Logistic model, well-imaged craters only (C ≥ 0.2; {R['logistic_C_ge_0.2']['n']:,} craters)", zorder=3)
a.barh(yy - h, [100 * Cn[g] for g in G], h, color=S3, label="Count model (negative binomial), hollows per crater", zorder=3)
a.set_yticks(yy); a.set_yticklabels([NAME[g] for g in G]); a.set_xlabel("Share of the variation explained by physical groups (%)")
a.grid(axis="x", color=GRID, lw=0.6, zorder=0); a.legend(frameon=False, loc="lower right", fontsize=6.8)
a.set_title(f"(a) Hosting with detection and crater area accounted for (R$^2$ beyond detection and area = {T['point']['full']:.2f})", loc="left", fontsize=8.5, color=INK)
for sp in ["top", "right"]: a.spines[sp].set_visible(False)
# (b) XRS rate ratios with block bootstrap
get, groups, fam = s.ANALYSES["thinned_presence_xrs_subset_lrm_xrs"]; dd = get(); y, off, fit, r2_all, dominance = s.make(dd, groups, fam)
cols = [c for g in groups for c in groups[g]]; b0, _ = fit(cols, np.ones(len(dd))); rng = np.random.default_rng(11); blocks = dd.block.unique(); BB = []
for _ in range(200):
    cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); BB.append(fit(cols, dd.block.map(cnt).fillna(0).values.astype(float))[0])
BB = np.array(BB); names = ["Intercept", "lnC"] + cols
terms = ["lrm", "s_si_z", "mg_si_z", "ca_si_z", "al_si_z"]; lab = {"lrm": "LRM fraction (0 to 1)", "s_si_z": "S/Si (per SD)", "mg_si_z": "Mg/Si (per SD)",
                                                                    "ca_si_z": "Ca/Si (per SD)", "al_si_z": "Al/Si (per SD)"}
rows = []
for t in terms:
    i = names.index(t); rows.append(dict(term=t, label=lab[t], rate_ratio=float(np.exp(b0[i])), lo=float(np.exp(np.percentile(BB[:, i], 2.5))), hi=float(np.exp(np.percentile(BB[:, i], 97.5)))))
X = pd.DataFrame(rows); X.to_csv(f"{OUT}/xrs_rate_ratios.csv", index=False)
b = ax[1]; yb = np.arange(len(X))[::-1]
b.axvline(1, color=MUTED, lw=0.8); b.errorbar(X.rate_ratio, yb, xerr=[X.rate_ratio - X.lo, X.hi - X.rate_ratio], fmt="o", ms=4, color=S1, ecolor=S1, lw=1.2)
b.set_xscale("log"); tk = [0.5, 0.7, 1, 1.4, 2, 3]; b.set_xticks(tk); b.set_xticklabels([f"{t:g}" for t in tk]); b.minorticks_off()
b.set_yticks(yb); b.set_yticklabels(X.label); b.grid(axis="x", color=GRID, lw=0.6, zorder=0)
P = R["thinned_presence_xrs_subset_lrm_only"]["point"]["region_add"]; Q = R["thinned_presence_xrs_subset_lrm_xrs"]["point"]["region_add"]
b.set_xlabel(f"Ratio of the density of hollows, all other groups held fixed\n(95% block-bootstrap interval); regional hosting fraction adds\nR$^2$ {P:.4f} without XRS and {Q:.4f} with XRS")
b.set_title(f"(b) Composition, {R['thinned_presence_xrs_subset_lrm_xrs']['n']:,} craters with XRS coverage", loc="left", fontsize=8.5, color=INK)
for sp in ["top", "right"]: b.spines[sp].set_visible(False)
fig.savefig(f"{OUT}/figS9_joint_model_{__import__('os').path.basename(__import__('os').environ['SM_DIR'])}.png", dpi=150, bbox_inches="tight"); print(X.round(3).to_string(index=False))
