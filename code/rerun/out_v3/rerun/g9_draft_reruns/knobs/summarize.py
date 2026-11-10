# g9 copy of g3_catalog/melt_source/knobs/summarize.py (6 Nov 2026). Unchanged; paths, catalog and grid come from common.py.
# g3_catalog copy of analysis_workshop/melt_source_tests/knobs/summarize.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Summary of the knob tests (5 Nov 2026): bootstrap intervals, knob_tests_record.json, fig_knob_color.png, fig_knob_hollows.png.
Inputs: _chunks/knob_color_0000_0335.csv, _chunks/peak_color_0000_1163.csv, montes_color.csv (Test 1);
knob_hollows.csv, _chunks/peak_hollows_00000_04824.csv, region_hollows.json (Test 2). (Other files in _chunks/ are superseded
trial runs on the first 10/40/200 rows and are not read.) Intervals: 95% percentile bootstrap, B = 4000, resampling knobs,
craters or sectors; for knob/matrix ratios also a spatial block bootstrap (3 x 3 deg cells) and the exact conditional
(Clopper-Pearson) interval that treats hollows as independent."""
import json, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.ticker
from common import HERE, AW, ORIG
rng = np.random.default_rng(20261105); B = 4000
def r(x, n=4): return None if x is None or not np.isfinite(x) else round(float(x), n)
def boot_stat(v, f=np.median):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; idx = rng.integers(0, len(v), (B, len(v)))
    return [r(f(v)), r(np.percentile(f(v[idx], axis=1), 2.5)), r(np.percentile(f(v[idx], axis=1), 97.5))]
# ---------------- Test 1 ----------------
METS = {"dln_r750": ("r750", "lnratio"), "dln_r560": ("r560", "lnratio"), "d_ln_1000_480": ("ln_1000_480", "diff"),
        "d_slope_ln750_430": ("slope_750_430", "diff"), "d_pc1": ("pc1", "diff"), "d_pc2": ("pc2", "diff"),
        "d_lrm_old_any": ("lrm_old_any", "diff"), "d_lrm_conc_any": ("lrm_conc_any", "diff"), "d_lrm_old_128": ("lrm_old_128", "diff"), "d_lrm_conc_128": ("lrm_conc_128", "diff")}
def deltas(d, a, b):
    o = pd.DataFrame(index=d.index)
    for nm, (m, kind) in METS.items():
        o[nm] = np.log(d[f"{a}_{m}"] / d[f"{b}_{m}"]) if kind == "lnratio" else d[f"{a}_{m}"] - d[f"{b}_{m}"]
    return o
def summarize(o):
    o = o.dropna(subset=["dln_r750", "d_ln_1000_480"]); S = {"n": int(len(o))}
    for nm in METS:
        v = o[nm].dropna()
        if len(v) < 3: continue
        S[nm] = {"median_ci": boot_stat(v), "mean_ci": boot_stat(v, np.mean), "frac_negative": r((v < 0).mean(), 3),
                 "wilcoxon_p": r(stats.wilcoxon(v).pvalue if (v != 0).any() else 1.0, 6)}
    x, y = o.dln_r750, o.d_ln_1000_480
    S["quadrants_r750_vs_1000_480"] = {"darker_bluer_LRM_like": r(((x < 0) & (y < 0)).mean(), 3), "brighter_bluer": r(((x > 0) & (y < 0)).mean(), 3),
                                       "darker_redder": r(((x < 0) & (y > 0)).mean(), 3), "brighter_redder": r(((x > 0) & (y > 0)).mean(), 3)}
    return S, o
kc = pd.read_csv(f"{ORIG}/_chunks/knob_color_0000_0335.csv")
kc = kc[(kc.kin_npx >= 3) & (kc.kmat_npx >= 20)]
pc = pd.read_csv(f"{ORIG}/_chunks/peak_color_0000_1163.csv")
kz = pd.read_csv(f"{AW}/kinczyk_degradation/kinczyk_matched.csv")[["id", "kclass"]].rename(columns={"id": "uid"})
pc = pc.merge(kz, on="uid", how="left"); pc["relief"] = pc.pk_z_med - pc.fl_z_med
pc = pc[(pc.pk_npx >= 3) & (pc.fl_npx >= 20)]; pcx = pc[pc.dcal_R > 2]
mc = pd.read_csv(f"{ORIG}/montes_color.csv")
T1 = {}; DL = {}
groups = {
 "knob_interior_vs_matrix": (kc, "kin", "kmat"),
 "knob_interior_vs_matrix__not_in_crater": (kc[~kc.in_crater], "kin", "kmat"),
 "knob_interior_vs_matrix__long_axis_ge_8km": (kc[kc.long_axis_km >= 8], "kin", "kmat"),
 "knob_interior_vs_matrix__certain_only": (kc[kc.KnobType == "Knob (certain)"], "kin", "kmat"),
 "knob_whole_disk_vs_matrix": (kc, "kall", "kmat"),
 "pseudo_knob_control": (kc[kc.pseudo_ok], "pin", "pmat"),
 "peak_vs_floor__cp_mp_30_250km_outside_2R_Caloris": (pcx, "pk", "fl"),
 "peak_vs_floor__DEM_verified_relief_gt_100m": (pcx[pcx.relief > 100], "pk", "fl"),
 "peak_vs_floor__D30_50km": (pcx[pcx.D < 50], "pk", "fl"),
 "peak_vs_floor__D50_250km": (pcx[pcx.D >= 50], "pk", "fl"),
 "peak_vs_floor__Kinczyk_class1_2": (pcx[pcx.kclass.isin([1, 2])], "pk", "fl"),
 "peak_vs_floor__Kinczyk_class3_5": (pcx[pcx.kclass.isin([3, 4, 5])], "pk", "fl"),
 "peak_vs_floor__cp_only": (pcx[pcx.cent_str == "cp"], "pk", "fl"),
 "wall_vs_floor__cp_mp_30_250km_outside_2R_Caloris": (pcx[pcx.wl_npx >= 20], "wl", "fl"),
 "wall_vs_floor__Kinczyk_class1_2": (pcx[(pcx.wl_npx >= 20) & pcx.kclass.isin([1, 2])], "wl", "fl"),
 "peak_vs_wall__cp_mp_30_250km_outside_2R_Caloris": (pcx[pcx.wl_npx >= 20], "pk", "wl"),
 "montes_high20_vs_low50": (mc, "hi", "lo"),
 "montes_high20_vs_rest80": (mc, "hi", "rest"),
 "montes_high10_vs_low30": (mc, "hi10", "lo30"),
 "montes_high20_vs_low50__outside_rim_only": (mc, "hi_ext", "lo_ext")}
for g, (d, a, b) in groups.items():
    T1[g], DL[g] = summarize(deltas(d, a, b))
# knob minus its own pseudo-knob (paired)
kk = kc[kc.pseudo_ok]; dd = deltas(kk, "kin", "kmat") - deltas(kk, "pin", "pmat"); T1["knob_minus_pseudo_paired"], _ = summarize(dd)
# knob vs peak contrast: difference of medians
def diffmed(v1, v2):
    v1, v2 = np.asarray(v1.dropna()), np.asarray(v2.dropna()); i1 = rng.integers(0, len(v1), (B, len(v1))); i2 = rng.integers(0, len(v2), (B, len(v2)))
    dm = np.median(v1[i1], 1) - np.median(v2[i2], 1); return [r(np.median(v1) - np.median(v2)), r(np.percentile(dm, 2.5)), r(np.percentile(dm, 97.5))]
T1["knob_minus_peak_median_difference"] = {m: diffmed(DL["knob_interior_vs_matrix"][m], DL["peak_vs_floor__cp_mp_30_250km_outside_2R_Caloris"][m]) for m in ["dln_r750", "d_ln_1000_480", "d_pc1", "d_lrm_old_any"]}
T1["knob_minus_montes_median_difference"] = {m: diffmed(DL["knob_interior_vs_matrix"][m], DL["montes_high20_vs_low50"][m]) for m in ["dln_r750", "d_ln_1000_480", "d_pc1", "d_lrm_old_any"]}
T1["diagnostics"] = {"knobs_used": int(len(kc)), "knob_interior_pixels_median": r(kc.kin_npx.median(), 1), "matrix_pixels_median": r(kc.kmat_npx.median(), 1),
    "knob_interior_DEM_higher_than_matrix_frac": r((kc.kin_z_med > kc.kmat_z_med).mean(), 3),
    "knob_interior_minus_matrix_DEM_median_m": r((kc.kin_z_med - kc.kmat_z_med).median(), 1),
    "peaks_used": int(len(pcx)), "peak_DEM_relief_gt0_frac": r((pcx.relief > 0).mean(), 3), "peak_DEM_relief_gt0_frac_D_ge_50": r((pcx[pcx.D >= 50].relief > 0).mean(), 3),
    "wall_LRM_old_any_mean": r(pcx.wl_lrm_old_any.mean(), 3), "floor_LRM_old_any_mean": r(pcx.fl_lrm_old_any.mean(), 3), "peak_LRM_old_any_mean": r(pcx.pk_lrm_old_any.mean(), 3),
    "montes_sectors": int(len(mc)), "montes_hi_minus_lo_detrended_relief_m_median": r(mc.relief_hi_lo_m.median(), 1),
    "knob_interior_LRM_old_any_mean": r(kc.kin_lrm_old_any.mean(), 3), "matrix_LRM_old_any_mean": r(kc.kmat_lrm_old_any.mean(), 3),
    "knob_interior_LRM_conc_any_mean": r(kc.kin_lrm_conc_any.mean(), 3), "matrix_LRM_conc_any_mean": r(kc.kmat_lrm_conc_any.mean(), 3)}
# ---------------- Test 2 ----------------
kh = pd.read_csv(f"{HERE}/knob_hollows.csv"); ph = pd.read_csv(f"{HERE}/_chunks/peak_hollows_00000_04824.csv")
reg = json.load(open(f"{HERE}/region_hollows.json"))
def dens(n, aw): return 1000 * n / aw
def boot_ratio_sums(num, den, blocks=None):
    num, den = np.asarray(num, float), np.asarray(den, float)
    if blocks is not None:
        u, inv = np.unique(blocks, return_inverse=True); num = np.bincount(inv, num); den = np.bincount(inv, den)
    idx = rng.integers(0, len(num), (B, len(num))); return num[idx].sum(1) / den[idx].sum(1)
def fields(col):
    f = [x for s in col.dropna().astype(str) for x in s.split(";") if x]; return len(f), len(set(f))
def knob_group(d, buffered=False):
    kn = "knobb" if buffered else "knob"
    nk, nm, ak, am = d[f"n_{kn}"].sum(), d.n_mat.sum(), d[f"Aw_{kn}"].sum(), d.Aw_mat.sum()
    bk = boot_ratio_sums(d[f"n_{kn}"], d[f"Aw_{kn}"]) * 1000; bm = boot_ratio_sums(d.n_mat, d.Aw_mat) * 1000
    idx = rng.integers(0, len(d), (B, len(d))); a = [d[c].values[idx].sum(1) for c in [f"n_{kn}", f"Aw_{kn}", "n_mat", "Aw_mat"]]
    rb = (a[0] / a[1]) / (a[2] / a[3])
    blk = (np.floor(d.lat / 3).astype(int) * 1000 + np.floor(d.lon / 3).astype(int)).values
    u, inv = np.unique(blk, return_inverse=True); sb = [np.bincount(inv, d[c].values) for c in [f"n_{kn}", f"Aw_{kn}", "n_mat", "Aw_mat"]]
    ib = rng.integers(0, len(u), (B, len(u))); rbb = (sb[0][ib].sum(1) / sb[1][ib].sum(1)) / (sb[2][ib].sum(1) / sb[3][ib].sum(1))
    N = nk + nm; lo_p, hi_p = stats.beta.ppf(0.025, nk, N - nk + 1) if nk > 0 else 0.0, stats.beta.ppf(0.975, nk + 1, N - nk)
    fk_pairs, fk_u = fields(d[f"fields_{kn}"]); fm_pairs, fm_u = fields(d.fields_mat)
    return {"n_knobs": int(len(d)), "knobs_with_hollows": int((d[f"n_{kn}"] > 0).sum()),
            "hollows_on_knobs": int(nk), "hollow_fields_on_knobs": fk_u, "knob_field_pairs": fk_pairs,
            "area_knobs_km2": r(d[f"A_{kn}"].sum(), 0), "area_w_knobs_km2": r(ak, 0), "mean_C_knobs": r(ak / d[f"A_{kn}"].sum(), 3),
            "density_knobs_per1000km2_w": [r(dens(nk, ak), 3), r(np.percentile(bk, 2.5), 3), r(np.percentile(bk, 97.5), 3)],
            "hollows_matrix": int(nm), "hollow_fields_matrix": fm_u, "matrix_field_pairs": fm_pairs, "area_matrix_km2": r(d.A_mat.sum(), 0),
            "area_w_matrix_km2": r(am, 0), "mean_C_matrix": r(am / d.A_mat.sum(), 3),
            "density_matrix_per1000km2_w": [r(dens(nm, am), 3), r(np.percentile(bm, 2.5), 3), r(np.percentile(bm, 97.5), 3)],
            "ratio_knob_over_matrix": r((nk / ak) / (nm / am), 3), "ratio_CI_knob_bootstrap": [r(np.nanpercentile(rb, 2.5), 3), r(np.nanpercentile(rb, 97.5), 3)],
            "ratio_CI_block_bootstrap_3deg": [r(np.nanpercentile(rbb[np.isfinite(rbb)], 2.5), 3), r(np.nanpercentile(rbb[np.isfinite(rbb)], 97.5), 3)],
            "ratio_CI_exact_conditional": [r(lo_p / (1 - lo_p) * am / ak, 3), r(hi_p / (1 - hi_p) * am / ak, 3)],
            "field_pair_ratio_knob_over_matrix": r((fk_pairs / ak) / (fm_pairs / am), 3) if fm_pairs else None}
T2 = {"knobs": {}}
kg = {"main_certain_uncertain_long_axis_ge_5km": kh[kh.main], "all_2011": kh, "certain_uncertain_all": kh[kh.KnobType != "Knob (<5 km)"],
      "lt5km_knobs": kh[kh.KnobType == "Knob (<5 km)"], "all_not_in_crater": kh[~kh.in_crater]}
for g, d in kg.items(): T2["knobs"][g] = knob_group(d)
T2["knobs"]["all_2011__footprint_plus_0.5km"] = knob_group(kh, buffered=True)
T2["knobs"]["main__footprint_plus_0.5km"] = knob_group(kh[kh.main], buffered=True)
T2["caloris_exterior_region"] = reg
def crater_group(d):
    o = {"n_craters": int(len(d)), "craters_with_peak_hollows": int((d.n_pk > 0).sum())}
    for nm in ["pk", "fl", "wl", "int"]:
        b = boot_ratio_sums(d[f"n_{nm}"], d[f"Aw_{nm}"]) * 1000
        o[f"{nm}_hollows"] = int(d[f"n_{nm}"].sum()); o[f"{nm}_fields"] = fields(d[f"fields_{nm}"])[1]
        o[f"{nm}_area_w_km2"] = r(d[f"Aw_{nm}"].sum(), 0); o[f"{nm}_mean_C"] = r(d[f"Aw_{nm}"].sum() / d[f"A_{nm}"].sum(), 3)
        o[f"{nm}_density_per1000km2_w"] = [r(dens(d[f"n_{nm}"].sum(), d[f"Aw_{nm}"].sum()), 3), r(np.percentile(b, 2.5), 3), r(np.percentile(b, 97.5), 3)]
    idx = rng.integers(0, len(d), (B, len(d)))
    for a_, b_ in [("pk", "fl"), ("pk", "wl")]:
        s = [d[c].values[idx].sum(1) for c in [f"n_{a_}", f"Aw_{a_}", f"n_{b_}", f"Aw_{b_}"]]; rr = (s[0] / s[1]) / (s[2] / s[3])
        o[f"ratio_{a_}_over_{b_}"] = [r((d[f"n_{a_}"].sum() / d[f"Aw_{a_}"].sum()) / (d[f"n_{b_}"].sum() / d[f"Aw_{b_}"].sum()), 3), r(np.nanpercentile(rr[np.isfinite(rr)], 2.5), 3), r(np.nanpercentile(rr[np.isfinite(rr)], 97.5), 3)]
    return o
cg = {"all_30_250km": ph, "all_30_250km_abslat_le_62": ph[ph.lat.abs() <= 62], "cp_mp_coded": ph[ph.cent_str.isin(["cp", "mp"])],
      "Kinczyk_class1_2": ph[ph.kclass.isin([1, 2])], "Kinczyk_class1_2_abslat_le_62": ph[ph.kclass.isin([1, 2]) & (ph.lat.abs() <= 62)],
      "Kinczyk_class1": ph[ph.kclass == 1], "Kinczyk_class2": ph[ph.kclass == 2], "Kinczyk_class3": ph[ph.kclass == 3], "Kinczyk_class3_5": ph[ph.kclass.isin([3, 4, 5])],
      "hollow_hosting_only_(conditional)": ph[ph.n_int > 0]}
T2["crater_peak_proxy"] = {g: crater_group(d) for g, d in cg.items()}
def knob_over_peak(dk, dp):
    i1 = rng.integers(0, len(dk), (B, len(dk))); i2 = rng.integers(0, len(dp), (B, len(dp)))
    a = dk.n_knob.values[i1].sum(1) / dk.Aw_knob.values[i1].sum(1); b = dp.n_pk.values[i2].sum(1) / dp.Aw_pk.values[i2].sum(1)
    q = a / b; q = q[np.isfinite(q)]
    return [r((dk.n_knob.sum() / dk.Aw_knob.sum()) / (dp.n_pk.sum() / dp.Aw_pk.sum()), 3), r(np.percentile(q, 2.5), 3), r(np.percentile(q, 97.5), 3)]
T2["knob_density_over_peak_density"] = {f"{a}__over__{b}": knob_over_peak(kg[a], cg[b]) for a in ["main_certain_uncertain_long_axis_ge_5km", "all_2011"]
                                        for b in ["all_30_250km", "Kinczyk_class1_2", "Kinczyk_class3", "Kinczyk_class3_5"]}
T2["paper_reference"] = "Paper: peak/wall = 6.7 within hollow-hosting craters (peak_wall_test), conditional on hosting."
REC = {"date": "2026-11-05", "test1_color": T1, "test2_hollow_density": T2,
       "notes": {"units": "dln_* = ln(knob/matrix) (fractional difference); d_ln_1000_480 and d_slope = differences of ln band ratios; PC1/PC2 from color_mdr/pca_record.json (PC1 + = redder); LRM = area fraction with 8-bit index > 0 (any) or >= 128, georeferenced Klima et al. 2018 (conc) and Peplowski et al. 2016 (old) maps, 5.1 km pixels",
                 "density": "hollows per 1,000 km2 of completeness-weighted area (completeness for 1 km hollows, 10 ppd)"}}
json.dump(REC, open(f"{HERE}/knob_tests_record.json", "w"), indent=1)
# ---------------- figures ----------------
BLUE, ORANGE, AQUA, GRAY, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#b9b8b2", "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(11.5, 8.6)); gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1], hspace=0.38, wspace=0.55)
ax = fig.add_subplot(gs[0, :2])
sets = [("Central peaks − floors (craters 30–250 km)", DL["peak_vs_floor__cp_mp_30_250km_outside_2R_Caloris"], ORANGE, 6),
        ("Caloris Montes − low cells (per 5° sector)", DL["montes_high20_vs_low50"], AQUA, 14),
        ("Pseudo-knobs on the matrix (control)", DL["pseudo_knob_control"], GRAY, 7),
        ("Knob interiors − matrix (Wright et al. 2020)", DL["knob_interior_vs_matrix"], BLUE, 9)]
for lab, o, c, s in sets:
    ax.scatter(o.dln_r750, o.d_ln_1000_480, s=s, color=c, alpha=0.35, linewidths=0, label=f"{lab}, n = {len(o)}")
for lab, o, c, s in sets:
    mx = boot_stat(o.dln_r750); my = boot_stat(o.d_ln_1000_480)
    ax.errorbar(mx[0], my[0], xerr=[[mx[0] - mx[1]], [mx[2] - mx[0]]], yerr=[[my[0] - my[1]], [my[2] - my[0]]], fmt="o", ms=8, color=c, mec="white", mew=1.5, elinewidth=2, capsize=0, zorder=5)
ax.axhline(0, color="#b9b8b2", lw=0.8, zorder=0); ax.axvline(0, color="#b9b8b2", lw=0.8, zorder=0)
ax.set_xlim(-0.3, 0.3); ax.set_ylim(-0.12, 0.08)
ax.set_xlabel("Δ ln reflectance at 750 nm (high-standing − surroundings)"); ax.set_ylabel("Δ ln(R1000 / R480)  (negative = bluer)")
for tx, ty, t in [(-0.29, -0.115, "darker and bluer (LRM-like)"), (0.29, -0.115, "brighter and bluer (immature-like)"), (-0.29, 0.073, "darker and redder"), (0.29, 0.073, "brighter and redder")]:
    ax.text(tx, ty, t, ha="left" if tx < 0 else "right", va="bottom" if ty < 0 else "top", color=MUTED, fontsize=8.5, style="italic")
ax.legend(loc="upper left", bbox_to_anchor=(1.04, 1.0), frameon=False, fontsize=8.5, markerscale=2.2, title="Contrast (point = one knob, crater or sector)", title_fontsize=8.5, alignment="left")
ax.set_title("A. Color contrast; large markers = medians with 95% bootstrap intervals", loc="left", fontsize=10, color=INK)
rowsF = [("Knob interior", "knob_interior_vs_matrix", BLUE), ("Whole knob disk", "knob_whole_disk_vs_matrix", BLUE), ("Pseudo-knob control", "pseudo_knob_control", GRAY),
         ("Peaks, all", "peak_vs_floor__cp_mp_30_250km_outside_2R_Caloris", ORANGE), ("Peaks, DEM-verified", "peak_vs_floor__DEM_verified_relief_gt_100m", ORANGE),
         ("Peaks, Kinczyk 1–2", "peak_vs_floor__Kinczyk_class1_2", ORANGE), ("Walls − floors (slope ref.)", "wall_vs_floor__cp_mp_30_250km_outside_2R_Caloris", "#c9a48f"),
         ("Montes", "montes_high20_vs_low50", AQUA), ("Montes, outside rim", "montes_high20_vs_low50__outside_rim_only", AQUA)]
for j, (m, xl, st) in enumerate([("dln_r750", "Δ ln R750 (median)", "median_ci"), ("d_ln_1000_480", "Δ ln(R1000/R480) (median)", "median_ci"), ("d_lrm_old_any", "Δ LRM fraction, Peplowski map (mean)", "mean_ci")]):
    a2 = fig.add_subplot(gs[1, j])
    for i, (lab, g, c) in enumerate(rowsF):
        v = T1[g][m][st]; y = len(rowsF) - 1 - i
        a2.plot([v[1], v[2]], [y, y], color=c, lw=2.2, solid_capstyle="round"); a2.plot(v[0], y, "o", color=c, ms=6, mec="white", mew=1.2)
    a2.axvline(0, color="#b9b8b2", lw=0.8); a2.set_yticks(range(len(rowsF))); a2.set_yticklabels([x[0] for x in rowsF][::-1] if j == 0 else [])
    a2.set_xlabel(xl); a2.tick_params(axis="y", length=0)
    if j == 0: a2.set_title("B. Contrasts with 95% bootstrap intervals (LRM: means, since most medians are 0 at 5.1 km pixels)", loc="left", fontsize=10, color=INK)
fig.savefig(f"{HERE}/fig_knob_color.png", dpi=160, bbox_inches="tight"); plt.close(fig)
# figure 2
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 5.0), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.75})
K = T2["knobs"]; C_ = T2["crater_peak_proxy"]
rg = reg["1.0-2.0R"]; lo_p, hi_p = stats.chi2.ppf(0.025, 2 * rg["n_hollows"]) / 2, stats.chi2.ppf(0.975, 2 * rg["n_hollows"] + 2) / 2
rowsD = [("Knobs ≥5 km, certain + uncertain (335)", K["main_certain_uncertain_long_axis_ge_5km"]["density_knobs_per1000km2_w"], BLUE),
         ("Knobs, all 2,011", K["all_2011"]["density_knobs_per1000km2_w"], BLUE),
         ("Adjacent matrix of the 335 (1.5–3 radii)", K["main_certain_uncertain_long_axis_ge_5km"]["density_matrix_per1000km2_w"], GRAY),
         ("Adjacent matrix of all 2,011", K["all_2011"]["density_matrix_per1000km2_w"], GRAY),
         ("Caloris exterior 1–2 radii, no craters/knobs (Poisson)", [rg["density_per_1000km2_w"], 1000 * lo_p / rg["area_w_km2"], 1000 * hi_p / rg["area_w_km2"]], GRAY),
         ("Peak proxy, all craters 30–250 km", C_["all_30_250km"]["pk_density_per1000km2_w"], ORANGE),
         ("Peak proxy, Kinczyk Classes 1–2", C_["Kinczyk_class1_2"]["pk_density_per1000km2_w"], ORANGE),
         ("Peak proxy, Kinczyk Class 3", C_["Kinczyk_class3"]["pk_density_per1000km2_w"], ORANGE),
         ("Peak proxy, Kinczyk Classes 3–5", C_["Kinczyk_class3_5"]["pk_density_per1000km2_w"], ORANGE),
         ("Floor proxy, Kinczyk Classes 1–2", C_["Kinczyk_class1_2"]["fl_density_per1000km2_w"], "#c9a48f"),
         ("Wall proxy, Kinczyk Classes 1–2", C_["Kinczyk_class1_2"]["wl_density_per1000km2_w"], "#c9a48f")]
for i, (lab, v, c) in enumerate(rowsD):
    y = len(rowsD) - 1 - i; a1.plot([max(v[1], 0.05), v[2]], [y, y], color=c, lw=2.2, solid_capstyle="round"); a1.plot(v[0], y, "o", color=c, ms=6, mec="white", mew=1.2)
a1.set_xscale("log"); a1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{v:g}")); a1.set_yticks(range(len(rowsD))); a1.set_yticklabels([x[0] for x in rowsD][::-1]); a1.tick_params(axis="y", length=0)
a1.set_xlabel("hollows per 1,000 km² of completeness-weighted area"); a1.set_title("A. Completeness-weighted hollow density", loc="left", fontsize=10, color=INK)
mk = K["main_certain_uncertain_long_axis_ge_5km"]; ak_ = K["all_2011"]
rowsR = [("Knob / matrix, 335 knobs", [mk["ratio_knob_over_matrix"]] + mk["ratio_CI_knob_bootstrap"], BLUE),
         ("Knob / matrix, all 2,011", [ak_["ratio_knob_over_matrix"]] + ak_["ratio_CI_knob_bootstrap"], BLUE),
         ("Peak / floor, all craters", C_["all_30_250km"]["ratio_pk_over_fl"], ORANGE),
         ("Peak / floor, Kinczyk 1–2", C_["Kinczyk_class1_2"]["ratio_pk_over_fl"], ORANGE),
         ("Peak / floor, Kinczyk 3–5", C_["Kinczyk_class3_5"]["ratio_pk_over_fl"], ORANGE),
         ("Peak / wall, Kinczyk 1–2", C_["Kinczyk_class1_2"]["ratio_pk_over_wl"], ORANGE)]
for i, (lab, v, c) in enumerate(rowsR):
    y = len(rowsR) - 1 - i; a2.plot([v[1], v[2]], [y, y], color=c, lw=2.2, solid_capstyle="round"); a2.plot(v[0], y, "o", color=c, ms=6, mec="white", mew=1.2)
a2.axvline(1, color="#b9b8b2", lw=0.8); a2.set_xscale("log"); a2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{v:g}")); a2.set_yticks(range(len(rowsR))); a2.set_yticklabels([x[0] for x in rowsR][::-1]); a2.tick_params(axis="y", length=0)
a2.set_xlabel("density ratio (95% bootstrap interval)"); a2.set_title("B. High-standing / adjacent terrain", loc="left", fontsize=10, color=INK)
fig.savefig(f"{HERE}/fig_knob_hollows.png", dpi=160, bbox_inches="tight"); plt.close(fig)
print(json.dumps({k: T1[k] for k in ["knob_interior_vs_matrix", "pseudo_knob_control", "peak_vs_floor__cp_mp_30_250km_outside_2R_Caloris", "montes_high20_vs_low50"]}, indent=0)[:200])
