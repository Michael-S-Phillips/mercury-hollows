# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/step9_color.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Summary model, step 9 (3 Nov 2026): MDIS 8-color composition (analysis_workshop/color_mdr, PDS MDR 64 ppd, brightness
removed by dividing each spectrum by its geometric mean; PC1 = spectral slope, PC2 = curvature with a low near 630 nm).
(1) Thinned-occurrence model (step 7) on craters with >= 20 valid hollow-masked interior pixels, with composition =
    LRM only, LRM + PC1 + PC2, PC1 + PC2 only, and LRM + interior-minus-surroundings contrasts of PC1 and PC2; dominance
    shares, the increase in R^2 from the regional hosting rate, and rate ratios per SD (200 block resamples of the full model).
    Shares and the regional increment for LRM + PC1 + PC2 get 60 block bootstraps, paired with LRM only.
(2) Depth: ln depth on the setting PC1 and PC2 within 10 km of each measured hollow (hollows masked), with ln pixel
    scale, ln tan(incidence), and |latitude| or a poleward flag; SEs clustered by image.
Writes out/color_record.json."""
import os, sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from multiprocessing import Pool
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/step9_color.py; only input/output paths changed. SM_DIR (env) is the
# folder holding this catalog's crater_features*.csv and out/ (g1_graph_joint/summary_model/<CATALOG>).
import step7_improved as s7
U = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts")
B = int(sys.argv[1]) if len(sys.argv) > 1 else 60
c = pd.read_csv(f"{U}/PlanetaryResearch_hollows/analysis_workshop/color_mdr/crater_color.csv", dtype={"id": str},
                usecols=["id", "ok", "hilat", "pc1", "pc2", "pc3", "ann_pc1", "ann_pc2", "albedo_gm", "slope", "curv"])
d0 = s7.d.copy(); d0["id"] = d0.id.astype(str)
d = d0.merge(c, on="id", how="left")
d = d[d.ok.fillna(False).astype(bool)].reset_index(drop=True)
for k in ["pc1", "pc2"]:
    d[k + "_z"] = (d[k] - d[k].mean()) / d[k].std()
    con = d[k] - d["ann_" + k]; d[k + "_con_z"] = ((con - con.mean()) / con.std()).fillna(0)
PG = s7.PG
VARIANTS = {"lrm_only": ["lrm"], "lrm_pc": ["lrm", "pc1_z", "pc2_z"], "pc_only": ["pc1_z", "pc2_z"], "lrm_contrast": ["lrm", "pc1_con_z", "pc2_con_z"]}
def mk(v): return s7.make(d, {**PG, "composition": VARIANTS[v]}, "cloglog")
def boot(args):
    v, seed = args; y, off, fit, r2_all, dom = mk(v)
    rng = np.random.default_rng(seed); blocks = d.block.unique()
    cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); w = d.block.map(cnt).fillna(0).values.astype(float)
    return v, dom(r2_all(w))
if __name__ == "__main__":
    R = {"n": int(len(d)), "events": int(d.h09.sum()), "B": B, "variants": {}}
    for v in VARIANTS:
        y, off, fit, r2_all, dom = mk(v); p = dom(r2_all(np.ones(len(d))))
        cols = [x for g in PG for x in PG[g]] + VARIANTS[v]; b0, _ = fit(cols, np.ones(len(d)))
        names = ["Intercept", "lnC"] + cols; rng = np.random.default_rng(21); blocks = d.block.unique(); BB = []
        for _ in range(200):
            cnt = pd.Series(rng.choice(blocks, len(blocks))).value_counts(); BB.append(fit(cols, d.block.map(cnt).fillna(0).values.astype(float))[0])
        BB = np.array(BB)
        rr = {t: [float(np.exp(b0[names.index(t)])), float(np.exp(np.percentile(BB[:, names.index(t)], 2.5))), float(np.exp(np.percentile(BB[:, names.index(t)], 97.5)))] for t in VARIANTS[v]}
        R["variants"][v] = {"point": p, "rate_ratios": rr}
        print(v, {k: round(x, 3) for k, x in p["share"].items()}, "full", round(p["full"], 4), "region add", round(p["region_add"], 4), {k: [round(z, 2) for z in x] for k, x in rr.items()}, flush=True)
    jobs = [(v, sd) for v in ["lrm_only", "lrm_pc"] for sd in range(1, B + 1)]
    with Pool(2) as pool: res = pool.map(boot, jobs, chunksize=3)
    a = [r for v, r in res if v == "lrm_only"]; b = [r for v, r in res if v == "lrm_pc"]
    q = lambda xs: [float(np.percentile(xs, 2.5)), float(np.percentile(xs, 97.5))]
    R["variants"]["lrm_pc"]["ci_share"] = {g: q([x["share"][g] for x in b]) for g in b[0]["share"]}
    R["variants"]["lrm_pc"]["ci_region_add"] = q([x["region_add"] for x in b])
    R["region_add_reduction_by_color"] = {"point": R["variants"]["lrm_only"]["point"]["region_add"] - R["variants"]["lrm_pc"]["point"]["region_add"],
                                          "ci": q([x["region_add"] - y["region_add"] for x, y in zip(a, b)])}
    R["composition_share_gain"] = {"point": R["variants"]["lrm_pc"]["point"]["share"]["composition"] - R["variants"]["lrm_only"]["point"]["share"]["composition"],
                                   "ci": q([y["share"]["composition"] - x["share"]["composition"] for x, y in zip(a, b)])}
    # (2) depth
    t = pd.read_csv(f"{U}/insolation_thermal/data/test1_rows_inc.csv"); h = pd.read_csv(f"{U}/PlanetaryResearch_hollows/analysis_workshop/color_mdr/hollow_color.csv")
    t["row"] = np.arange(len(t)); t = t.merge(h[["row", "set_pc1", "set_pc2", "all_pc1"]], on="row", how="left").dropna(subset=["set_pc1", "set_pc2"])
    t["ltan"] = np.log(np.tan(np.radians(t.inc))); t["pole"] = (t.abslat >= 60).astype(int); g = pd.factorize(t.filename)[0]
    for k in ["set_pc1", "set_pc2"]: t[k + "_z"] = (t[k] - t[k].mean()) / t[k].std()
    D = {}
    for name, f in [("pole", "lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan + pole"), ("abslat", "lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan + abslat"),
                    ("none", "lnh ~ set_pc1_z + set_pc2_z + lnpx + ltan")]:
        m = smf.ols(f, t).fit(cov_type="cluster", cov_kwds={"groups": g})
        D[name] = {k: [float(np.exp(m.params[k])), float(np.exp(m.params[k] - 1.96 * m.bse[k])), float(np.exp(m.params[k] + 1.96 * m.bse[k])), float(m.pvalues[k])] for k in ["set_pc1_z", "set_pc2_z"]}
        # 10th-90th percentile range in SD units
        D[name]["p10_p90_sd"] = {k: float((t[k].quantile(.9) - t[k].quantile(.1)) / t[k].std()) for k in ["set_pc1", "set_pc2"]}
    R["depth"] = {"n": int(len(t)), "fits": D}
    json.dump(R, open(f"{s7.OUT}/color_record.json", "w"), indent=1)
    print(json.dumps({k: R[k] for k in ["region_add_reduction_by_color", "composition_share_gain", "depth"]}, indent=1))
    print("lrm_pc CI", R["variants"]["lrm_pc"]["ci_share"], R["variants"]["lrm_pc"]["ci_region_add"])
