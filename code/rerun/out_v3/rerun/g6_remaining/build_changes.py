#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script. Writes changes_g6.csv from this folder's outputs (and, where stated, g1/g2/g3 outputs).
"""Columns as in the other groups' tables: section, quantity, text_value (as printed in v6.3 / supplement v13 or the workshop note),
v2_record (stored record, or this folder's v2 run where that reproduces the record exactly), v3, v3b, script, conclusion_change."""
import os, json, csv
H = os.path.dirname(os.path.abspath(__file__)); RR = os.path.dirname(H)
J = lambda p: json.load(open(p))
CATS = ["v2", "v3", "v3b"]
def ci(v, n=2): return f"{v[0]:.{n}f} ({v[1]:.{n}f}-{v[2]:.{n}f})"
rows = []
def add(section, quantity, text, vals, script, change="no"):
    rows.append(dict(section=section, quantity=quantity, text_value=text, v2_record=vals[0], v3=vals[1], v3b=vals[2], script=script, conclusion_change=change))
# ---------------- 1. Stage classifications against incidence (S8)
S = {c: J(f"{H}/stage/out_{c}/stage_incidence.json") for c in CATS}
f = lambda r: f"{r['OR_stage2_per10deg']:.2f} ({r['lo']:.2f}-{r['hi']:.2f}); n = {r['n']:,}"
add("S8", "Odds of a Stage 2 or 3 classification per 10 deg of solar incidence: all labeled detections of the 0.8 release (controls |lat|, log10 pixel scale, emission; the record's model)",
    "double (record: x2.0 per 10 deg)", [f(S[c]["A_release_all"]["fit"]) for c in CATS], "stage/stage_incidence.py (fit A; the release does not depend on the catalog)")
add("S8", "Same, labeled hollows of each catalog, Stage 1 vs Stage 2 (the S8 sample; SEs clustered by field)", "double",
    [f(S[c]["B1_hollows_stage1_vs_2"]) for c in CATS], "stage/stage_incidence.py (fit B1)")
add("S8", "Same, labeled hollows, Stage 1 vs Stage 2 or 3", "double", [f(S[c]["B2_hollows_stage1_vs_2or3"]) for c in CATS], "stage/stage_incidence.py (fit B2)")
add("S8", "Incidence term of the S8 Stage 1 model (Delta T_w, incidence, pixel scale, emission, |lat|), as odds of Stage 2 per 10 deg", "double",
    [f(S[c]["B3_S8_model_incidence_term"]) for c in CATS], "stage/stage_incidence.py (fit B3)")
L = {c: J(f"{RR}/g2_insolation/out/{c}/stage_test_lowinc.json")["inc_lt70"] for c in CATS}
def frac(d):
    s = d["s1_share_by_dT"]; lo = [x["s1"] for x in s if x["dT"] in ("(0, 20]", "(20, 30]", "(30, 40]", "(40, 50]")]; hi = [x["s1"] for x in s if x["dT"] not in ("(0, 20]", "(20, 30]", "(30, 40]", "(40, 50]")]
    return f"{100*min(lo):.1f}-{100*max(lo):.1f}% / {100*min(hi):.1f}-{100*max(hi):.1f}%; step {d['best_t']} K {ci([d['step'][k] for k in ('OR','lo','hi')])}; AIC linear - step {d['AIC']['linear'] - d['AIC']['step_best']:.1f}; n = {d['n']:,}"
add("3.5; S8", "Stage 1 fractions below ~45 K / above 50 K, 44 K step and AIC (incidence < 70 deg): confirmed from g2's outputs, not rerun here",
    "56-71% / 96-100%; 44 K, 3.0 (1.5-6.2); AIC 27; 1,647 hollows", [frac(L[c]) for c in CATS], "g2_insolation/stage_test_lowinc.py (g2 outputs)", "no (already in changes_g2.csv)")
# ---------------- 2. Color (S9, S16, main text 3.7 and 3.8)
P = {c: J(f"{H}/color/out_{c}/pca_record.json") for c in CATS}; D = {c: J(f"{H}/color/out_{c}/diagnostics_record.json") for c in CATS}
DL = {c: J(f"{H}/color/out_{c}/depth_color_lrm_record.json") for c in CATS}; CR = {c: J(f"{H}/summary_model/{c}/out/color_record.json") for c in CATS}
DS = {c: J(f"{H}/color/out_{c}/depth_sensitivity.json") for c in CATS}; BF = {c: J(f"{H}/summary_model/{c}/out/boot_failure_check.json") for c in CATS}
add("S9 (color)", "Cataloged hollows masked in the color map (1 km + min(D_eq, 5 km)/2)", "6,731", ["6,731", "6,363", "6,131"], "color/extract_tiles.py (mask: inputs/hollows_mask_<cat>.csv)")
add("S9; S16; Fig S7", "Craters with >= 20 valid hollow-masked interior pixels (share of all craters)", "31,373 (99.4%)",
    [f"{P[c]['n_pca']:,} ({100*P[c]['coverage']['frac_ge20_all']:.2f}%)" for c in CATS], "color/merge_spectra.py, pca_color.py")
add("S9; 3.7", "Variance explained by PC1 / PC2", "88% / 5%", [f"{100*P[c]['variance_ratio'][0]:.2f}% / {100*P[c]['variance_ratio'][1]:.2f}%" for c in CATS], "color/pca_color.py")
add("S9", "Correlation of PC1 with ln(R750/R430) (Pearson)", "0.93", [f"{P[c]['correlations']['PC1']['slope']['pearson']:.3f}" for c in CATS], "color/pca_color.py")
add("S9", "Spearman rho of PC2 / PC1 with the crater's concentrated-LRM fraction", "0.23 / -0.22",
    [f"{P[c]['correlations']['PC2']['lrm']['spearman']:.3f} / {P[c]['correlations']['PC1']['lrm']['spearman']:.3f}" for c in CATS], "color/pca_color.py")
add("S9", "Spearman rho of PC1 (slope) with |lat| across craters", "0.29", [f"{P[c]['correlations']['PC1']['abslat']['spearman']:.3f}" for c in CATS], "color/pca_color.py")
dep = lambda c, k, t: ci(CR[c]["depth"]["fits"][k][t][:3], 3)
add("S9; 3.7; Fig S7b", "Depth factor per SD of the slope (PC1) of the surroundings, |lat| control", "0.90 (0.86-0.95)", [dep(c, "abslat", "set_pc1_z") for c in CATS], "summary_model/step9_color.py (part 2)")
add("S9; 3.7", "Depth factor per SD of the curvature (PC2), |lat| control", "1.01 (0.97-1.06)", [dep(c, "abslat", "set_pc2_z") for c in CATS], "summary_model/step9_color.py")
add("S9", "Slope factor with the poleward indicator / without a latitude control", "0.92 (0.88-0.96) / 0.91 (0.87-0.94)",
    [dep(c, "pole", "set_pc1_z") + " / " + dep(c, "none", "set_pc1_z") for c in CATS], "summary_model/step9_color.py")
ds = lambda r: f"{r['factor']:.3f} ({r['lo']:.3f}-{r['hi']:.3f})"
add("S9", "Slope factor within 60 deg / with ln geometric-mean reflectance / with 10-deg latitude bands", "0.90 (0.86-0.94) / 0.89 / 0.91",
    [" / ".join(ds(DS[c][k]) for k in ("within_60", "with_ln_albedo", "lat_bands_10deg")) for c in CATS], "color/depth_sensitivity.py (new; no script was kept)")
add("S9", "10th-90th percentile range of the slope in SD; fitted depth change across it", "2.4 SD; 0.78-0.82",
    [f"{CR[c]['depth']['fits']['abslat']['p10_p90_sd']['set_pc1']:.2f} SD; " + f"{min(CR[c]['depth']['fits'][k]['set_pc1_z'][0] for k in ('pole','abslat','none'))**CR[c]['depth']['fits']['abslat']['p10_p90_sd']['set_pc1']:.2f}-{max(CR[c]['depth']['fits'][k]['set_pc1_z'][0] for k in ('pole','abslat','none'))**CR[c]['depth']['fits']['abslat']['p10_p90_sd']['set_pc1']:.2f}" for c in CATS], "summary_model/step9_color.py")
add("S9", "Spearman rho of the setting slope with |lat| / with the LRM fraction within 10 km", "0.65 / -0.46",
    [f"{D[c]['hollow_setting_pc1_vs_abslat']:.3f} / {DL[c]['rho_pc1_lrm']:.3f}" for c in CATS], "color/diagnostics.py; color/depth_color_lrm.py")
add("S9; 3.7", "With the LRM fraction and the poleward indicator: slope factor per SD; LRM factor none to complete", "0.93 (0.89-0.97); 1.08 (0.91-1.28)",
    [ci(DL[c]["pole"]["set_pc1_z"][:3], 3) + "; " + ci(DL[c]["pole"]["lrm_local"][:3], 3) for c in CATS], "color/depth_color_lrm.py")
add("S9", "Median share of pixels within 10 km of a measured hollow removed by the mask; measurements whose unmasked surroundings are bluer", "17%; 76%",
    [f"{100*P[c]['hollow_settings']['median_frac_masked']:.1f}%; {100*(1-D[c]['hollow_unmasked_minus_masked_10km']['pc1']['frac_positive']):.1f}%" for c in CATS], "color/pca_color.py; color/diagnostics.py")
add("S16; Fig S7c", "Color subset of the joint model: craters (hosts; hosts lost to the hollow mask)", "31,373 (563; 43)",
    [f"{CR[c]['n']:,} ({CR[c]['events']}; {[606, 446, 370][i] - CR[c]['events']})" for i, c in enumerate(CATS)], "summary_model/step9_color.py",
    "no (g1, with the v2 mask: 442 (4) / 366 (4) of 31,373)")
sh = lambda c, v: CR[c]["variants"][v]["point"]["share"]
add("S16; 3.8", "Low-reflectance share: LRM only -> LRM + PC1 + PC2 (95% block-bootstrap interval)", "1.5% -> 6.3% (2.9-9.8%); main text 6% (3-10%)",
    [f"{100*sh(c,'lrm_only')['composition']:.1f}% -> {100*sh(c,'lrm_pc')['composition']:.1f}% ({100*CR[c]['variants']['lrm_pc']['ci_share']['composition'][0]:.1f}-{100*CR[c]['variants']['lrm_pc']['ci_share']['composition'][1]:.1f}%)" for c in CATS],
    "summary_model/step9_color.py", "no (g1's v3b interval, -1.0-10.1%, came from the v2 mask; with the v3b mask it is 0.6-11.2%)")
add("S16", "Gain in the low-reflectance share, percentage points (paired interval)", "4.8 (2.0-7.8)",
    [f"{100*CR[c]['composition_share_gain']['point']:.1f} ({100*CR[c]['composition_share_gain']['ci'][0]:.1f}-{100*CR[c]['composition_share_gain']['ci'][1]:.1f})" for c in CATS], "summary_model/step9_color.py")
order = lambda c: ", ".join(f"{k} {100*v:.1f}" for k, v in sorted(sh(c, "lrm_pc").items(), key=lambda x: -x[1]) if k != "composition")
add("S16; 3.8", "Shares of the other groups with color added (order)", "degradation 32, insolation 22, size 15, nesting 15, cover 6.5, vents 3",
    [order(c) for c in CATS], "summary_model/step9_color.py", "yes (v3, v3b; as g1: nesting now ahead of size)")
rr = lambda c, t: ci(CR[c]["variants"]["lrm_pc"]["rate_ratios"][t])
add("S16", "Density ratio per SD of slope / per SD of curvature / LRM none to complete", "1.06 (0.87-1.25) / 0.90 (0.79-1.06) / 1.34 (1.12-1.66)",
    [" / ".join(rr(c, t) for t in ("pc1_z", "pc2_z", "lrm")) for c in CATS], "summary_model/step9_color.py")
add("S16; 3.8", "Regional R^2 gain without -> with color; reduction (interval)", "0.0307 -> 0.0297; 0.0010 (-0.0009 to 0.0029); main text 0.001-0.002 (with XRS)",
    [f"{CR[c]['variants']['lrm_only']['point']['region_add']:.4f} -> {CR[c]['variants']['lrm_pc']['point']['region_add']:.4f}; {CR[c]['region_add_reduction_by_color']['point']:.4f} ({CR[c]['region_add_reduction_by_color']['ci'][0]:.4f} to {CR[c]['region_add_reduction_by_color']['ci'][1]:.4f})" for c in CATS],
    "summary_model/step9_color.py")
add("S16 (diagnostic)", "Failed bootstrap draws (cloglog IRLS) in step9, LRM only / LRM + color, of 60", "not quoted",
    [f"{BF[c]['step9_color:lrm_only']['failed']} / {BF[c]['step9_color:lrm_pc']['failed']}" for c in CATS], "summary_model/boot_failure_check.py (g1's copy)", "n/a (widens intervals; see README)")
# ---------------- 3. Knob color tests (workshop)
K = {c: J(f"{H}/knobs/out_{c}/knob_tests_record.json")["test1_color"] for c in CATS}
kq = lambda t: f"dlnR750 {t['dln_r750']['median_ci'][0]:+.3f} [{t['dln_r750']['median_ci'][1]:.3f}, {t['dln_r750']['median_ci'][2]:.3f}]; dln(R1000/R480) {t['d_ln_1000_480']['median_ci'][0]:+.3f} [{t['d_ln_1000_480']['median_ci'][1]:.3f}, {t['d_ln_1000_480']['median_ci'][2]:.3f}]; n = {t['n']}"
add("workshop: knobs", "Knob interior - matrix (334 knobs >= 5 km): median contrasts; brighter and bluer / darker and bluer", "+0.061 [0.056, 0.067]; -0.033 [-0.037, -0.030]; 82% / 8%",
    [kq(K[c]["knob_interior_vs_matrix"]) + f"; {100*K[c]['knob_interior_vs_matrix']['quadrants_r750_vs_1000_480']['brighter_bluer']:.1f}% / {100*K[c]['knob_interior_vs_matrix']['quadrants_r750_vs_1000_480']['darker_bluer_LRM_like']:.1f}%" for c in CATS],
    "knobs/t1_knob_color.py + summarize.py")
add("workshop: knobs", "Pseudo-knob control (same geometry 5a away)", "0.000", [kq(K[c]["pseudo_knob_control"]) for c in CATS], "knobs/t1_knob_color.py + summarize.py")
for key, lab, txt in [("peak_vs_floor__cp_mp_30_250km_outside_2R_Caloris", "Central peaks - floors (cp/mp craters 30-250 km)", "+0.032 and -0.010"),
                      ("peak_vs_floor__Kinczyk_class1_2", "Age-matched peaks - floors (Kinczyk Class 1-2)", "+0.017 and -0.003"),
                      ("wall_vs_floor__cp_mp_30_250km_outside_2R_Caloris", "Walls - floors (slope reference)", "0.000 and -0.009"),
                      ("montes_high20_vs_low50", "Caloris Montes high 20% - low 50% (72 sectors)", "-0.008 and -0.006")]:
    add("workshop: knobs", lab, txt, [kq(K[c][key]) for c in CATS], "knobs/t1_peak_color.py / t1_montes_color.py + summarize.py")
add("workshop: knobs", "Test 2 (hollow density) in the same record", "see changes_g3.csv", ["original", "identical to g3", "point estimates identical to g3; bootstrap intervals differ slightly (shared random stream with Test 1); knobs >= 5 km block interval 0.77-28.6 (g3: 0.75-27.6)"],
    "knobs/summarize.py (reads g3's Test 2 outputs)", "no (g3's v3b 'yes' for the block interval holds)")
# ---------------- 4. DEM aspect check
A = {c: J(f"{H}/aspect/out_{c}/aspect_robustness_record.json") for c in CATS}
ca = lambda r: r["decompose_15_30N"]["by_caloris"]["ge1200km"]
af = lambda r: f"{r['frac_poleward']:.2f} ({r['frac_poleward_ci95'][0]:.2f}-{r['frac_poleward_ci95'][1]:.2f}); {r['per_crater_frac_poleward']:.2f} ({r['per_crater_ci95'][0]:.2f}-{r['per_crater_ci95'][1]:.2f}); {r['n']} hollows, {r['n_craters']} craters"
add("S8; Table S6", "Equator-facing fraction of northern wall hollows at 15-30N, Caloris excluded: per hollow; per host crater; n", "0.45 (0.41-0.48); 0.54 (0.42-0.66); 747 hollows, 53 craters",
    [af(ca(A[c])) for c in CATS], "aspect/aspect_robustness.py (agrees with g2 wall_fractions.py)", "no")
dp = lambda r: f"{r['dem_proxy_check']['agreement_rate_lt90deg']:.3f} / {r['dem_proxy_check']['agreement_rate_lt45deg']:.3f}; median {r['dem_proxy_check']['median_abs_diff_deg']:.1f} deg; n = {r['dem_proxy_check']['n']:,}"
add("not in v6.3 (aspect_robustness)", "DEM (8 ppd) slope azimuth vs position-based facing: agreement < 90 / < 45 deg", "not quoted", [dp(A[c]) for c in CATS], "aspect/aspect_robustness.py with _scratch/dem8_gis.npz")
# ---------------- 5. nested_test2
NS = J(f"{H}/nested/nested_test2_search.json")
add("none (peak_wall_test/nested_test2.py)", "Numbers of nested_test2_record.json found in v6.3 or supplement v13", "none expected (g1)",
    [f"{NS['n_values']} values searched; {len(NS['hits'])} string hits, all coincidental (no k/N count, sample size or coefficient of the record appears)", "not rerun", "not rerun"], "nested/search_nested_test2.py", "no")
# ---------------- 6. Reference size (S2), not rerun by g2/g3
RS = {c: J(f"{H}/refsize/out_{c}/ref_size_record.json") for c in CATS}; RI = {c: J(f"{H}/refsize/out_{c}/ref_size_insolation.json") for c in CATS}
HF = {c: J(f"{H}/refsize/out_{c}/hosting_fits.json") for c in CATS}; JT = {c: J(f"{H}/refsize/out_{c}/joint_thinned.json") for c in CATS}
rng = lambda xs, n=2: f"{min(xs):.{n}f}-{max(xs):.{n}f}"
def jl(c):
    o = {k[7:]: JT[c][k]["llf_full"] for k in JT[c] if k.startswith("lnCold_")}; b = max(o, key=o.get)
    return f"1 km {o['1.0']:,.1f}; others {max(v for k, v in o.items() if k != '1.0'):,.1f} to {min(v for k, v in o.items() if k != '1.0'):,.1f}; best {b} km"
add("S2", "Joint model log-likelihood with completeness for 1 km hollows vs 0.3-5 km (reason for choosing 1 km)", "-1,969 against -1,976 to -2,229",
    [jl(c) for c in CATS], "refsize/joint_thinned.py", "yes (v3, v3b: completeness for 0.7 km hollows now fits better than 1 km, by 4.0 and 6.0; 1 km is second)")
add("S2", "Containment-test log-likelihood with completeness for 1 vs 2 km hollows", "-2,100 and -2,093 (about equally well)",
    [f"{RS[c]['by_ref']['1.0']['llf']:,.1f} and {RS[c]['by_ref']['2.0']['llf']:,.1f}" for c in CATS], "refsize/ref_size_check.py", "qualified (v3, v3b: 1 km fits better by 10-12, which supports the choice)")
for lab, k, txt in [("craters inside older craters", "nested_OR", "1.85-1.91"), ("craters with volcanic vents", "vent_OR", "5.0-6.7"), ("fresh vs degraded craters >= 10 km", "fresh_vs_degraded_OR", "3.3-4.5")]:
    add("S2", f"Odds ratio with completeness for 0.3, 0.5 or 1 km hollows: {lab}", txt, [rng([RS[c]["by_ref"][r][k][0] for r in ("0.3", "0.5", "1.0")]) for c in CATS], "refsize/ref_size_check.py")
for lab, b, txt in [("hot vs warm, within 15 deg", "-15..15", "2.2-2.6"), ("hot vs warm, 15-30N", "15..30", "3.2-3.7"), ("hot vs warm, 30-45N", "30..45", "1.76-1.98")]:
    add("S2", f"Odds ratio with completeness for 0.3, 0.5 or 1 km hollows: {lab}", txt, [rng([RI[c][r][b]["hot_OR"][0] for r in ("0.3", "0.5", "1.0")]) for c in CATS], "refsize/ref_size_insolation.py",
        "yes (v3, v3b: the 1 km value is no longer significant, as in changes_g2.csv)" if b == "30..45" else "no")
add("S2", "Odds per 10 K of Delta T_w within 60 deg, completeness for 0.3, 0.5 or 1 km", "1.43-1.47", [rng([RI[c][r]["dT_per10K_lt60"][0] for r in ("0.3", "0.5", "1.0")]) for c in CATS], "refsize/ref_size_insolation.py")
def sig(c):
    o = []
    for r in ("0.3", "0.5", "1.0", "2.0", "5.0"):
        x = RI[c][r]["30..45"]["hot_OR"]; o.append(f"{r} km {x[0]:.2f} ({x[1]:.2f}-{x[2]:.2f}){'' if x[3] < 0.05 else ' n.s.'}")
    return "; ".join(o)
add("S2", "Hot-longitude excess at 30-45N by reference size (significance)", "loses significance only with 2 or 5 km (1.4, 0.9-2.2; 1.0, 0.6-1.6)", [sig(c) for c in CATS], "refsize/ref_size_insolation.py",
    "yes (v3: also n.s. at 1 km; v3b: significant only at 0.5 km)")
E = lambda c: HF[c]["lnE_1.6"]; C1 = lambda c: HF[c]["C1"]
add("S2", "Size-weighted completeness E (exponent 1.6): odds ratios for nesting, vents, fresh vs degraded (against C1)", "1.91 (1.59-2.29), 5.6 (3.4-9.2), 3.6 (2.4-5.4); against 1.85, 5.0, 3.3",
    [f"{ci(E(c)['nested_OR'])}, {ci(E(c)['vent_OR'],1)}, {ci(E(c)['fresh_vs_degraded_OR'],1)}; against {C1(c)['nested_OR'][0]:.2f}, {C1(c)['vent_OR'][0]:.2f}, {C1(c)['fresh_vs_degraded_OR'][0]:.2f}" for c in CATS], "refsize/hosting_fits.py")
add("S2", "E: hot vs warm odds ratios within 15 deg, 15-30N, 30-45N (against C1)", "2.0 (1.5-2.7), 2.8 (1.9-4.2), 1.70 (1.08-2.68); against 2.2, 3.2, 1.76",
    [f"{ci(E(c)['hot_OR_-15..15'][:3],1)}, {ci(E(c)['hot_OR_15..30'][:3],1)}, {ci(E(c)['hot_OR_30..45'][:3])}; against {C1(c)['hot_OR_-15..15'][0]:.2f}, {C1(c)['hot_OR_15..30'][0]:.2f}, {C1(c)['hot_OR_30..45'][0]:.2f}" for c in CATS],
    "refsize/hosting_fits.py", "yes (v3, v3b: 30-45N not significant with E either)")
add("S2", "Odds per 10 K of Delta T_w within 60 deg with E / with C1", "1.43 (1.32-1.54) with either measure",
    [f"{ci(E(c)['dT_OR_per10K_lt60'])} / {ci(C1(c)['dT_OR_per10K_lt60'])}" for c in CATS], "refsize/hosting_fits.py")
ks = ["nested_OR", "vent_OR", "fresh_vs_degraded_OR", "hot_OR_-15..15", "hot_OR_15..30", "hot_OR_30..45", "dT_OR_per10K_lt60"]
add("S2", "Largest change of these odds ratios with exponents 1.3 and 1.9 instead of 1.6", "less than 5%",
    [f"{100*max(abs(HF[c][a][k][0]/E(c)[k][0]-1) for a in ('lnE_1.3','lnE_1.9') for k in ks):.1f}%" for c in CATS], "refsize/hosting_fits.py", "qualified (v3b: 4.95%, at the limit)")
add("S2", "Joint model with E: largest change in group shares; log-likelihood (against C1)", "at most 1.2 percentage points; -2,022 (against -1,969)",
    [f"{100*max(abs(JT[c]['lnE_1.6']['share'][g]-JT[c]['lnC1']['share'][g]) for g in JT[c]['lnC1']['share']):.1f} pp; {JT[c]['lnE_1.6']['llf_full']:,.1f} (against {JT[c]['lnC1']['llf_full']:,.1f})" for c in CATS],
    "refsize/joint_thinned.py", "qualified (v3 1.7 pp, v3b 2.1 pp; E still fits worse)")
# ---------------- 6. other uncovered numbers checked here
V = {c: J(f"{RR}/g1_graph_joint/vents/{c}/test8_vents_record.json") for c in CATS}
add("3.4; S15", "Depth of hollows in host craters with vents vs other hosts (13 hosts, 405 hollows)", "1.10 (0.96-1.28)",
    [f"{ci([V[c]['depth_ratio_cluster_image'][k] for k in ('ratio','lo','hi')])}; {V[c]['depth_counts']['vent_hosts']} hosts, {V[c]['depth_counts']['in_vent_hosts']} hollows" for c in CATS],
    "g1_graph_joint/vents/test8_vents.py (g1 outputs; not in changes_g1.csv)")
W = {c: J(f"{RR}/g3_catalog/terrain/out_{c}/peak_wall_record.json")["rows_H250_r0.45"]["all"]["slope_stratified"] for c in CATS}
wa = lambda s: f"{100*s['s0']['A_wall']/(s['s0']['A_wall']+s['s1']['A_wall']+s['s2']['A_wall']):.1f}% (wall area < 5 deg in the terrain test)"
add("S8", "Share of host-crater wall area with regional slopes < 5 deg", "65%", [wa(W[c]) for c in CATS], "g3_catalog/terrain/peak_wall.py outputs (slope_stratified)",
    "not rerun (the script behind 65% was not found; the wall area < 5 deg of the 352 / 321 / 262 terrain-test hosts is 68.8 / 69.1 / 69.0%, so the screen barely moves it)")
for r in rows:   # g4_completeness, written at the same time, also reruns S2 (and with rebuilt completeness grids)
    if r["section"] == "S2": r["script"] += "; also rerun by g4_completeness (identical v3 values with the unchanged grids)"
with open(f"{H}/changes_g6.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["section", "quantity", "text_value", "v2_record", "v3", "v3b", "script", "conclusion_change"]); w.writeheader(); w.writerows(rows)
print(len(rows), "rows")
