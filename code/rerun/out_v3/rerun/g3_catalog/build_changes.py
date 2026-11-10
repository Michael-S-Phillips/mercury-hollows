#!/usr/bin/env python3
"""g3_catalog (6 Nov 2026): builds changes_g3.csv from the rerun outputs. text_value = value in Hollow_Formation_Model_Paper_v6_3.tex /
supplementary-materials_v13.tex (or 'not in v6.3' for workshop analyses); v2_record = value in the existing record JSON (each v2 rerun
here reproduces it unless noted); v3 / v3b = this rerun (size distribution with the out_v2 detection model unless the row says refit)."""
import os, json, csv
G = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog"); HE = os.path.expanduser("~/mnt/hollows-emissivity")
PR = f"{HE}/manuscripts/PlanetaryResearch_hollows"; J = lambda p: json.load(open(p))
CC = J(f"{G}/counts/catalog_counts.json")
SN = {v: J(f"{G}/sfd/out_{v}_fixed/sfd_numbers.json") for v in ("v3", "v3b")}; SN["rec"] = J(f"{PR}/analysis_workshop/completeness_sfd/sfd_numbers.json")
SF = {v: J(f"{G}/sfd/out_{v}_fixed/sfd_fit.json") for v in ("v3", "v3b")}; SF["rec"] = J(f"{PR}/analysis_workshop/completeness_sfd/sfd_fit.json")
DF, DR = J(f"{G}/sfd/sfd_derived_fixed.json"), J(f"{G}/sfd/sfd_derived_refit.json")
SL = {v: J(f"{G}/sfd/size_law_record_{v}.json")["size_" + v] for v in ("v3", "v3b")}; SL["rec"] = J(f"{HE}/code/rerun/closure/size_law_record.json")["size_v2"]
CA = {v: J(f"{G}/caloris/out_{v}/caloris_v9_record.json") for v in ("v3", "v3b")}; CA["rec"] = J(f"{PR}/analysis_workshop/caloris_v9/caloris_v9_record.json")
CT = {v: J(f"{G}/caloris/out_{v}/caloris_v9_controls_record.json") for v in ("v3", "v3b")}; CT["rec"] = J(f"{PR}/analysis_workshop/caloris_v9/caloris_v9_controls_record.json")
CO = {v: J(f"{G}/caloris/out_{v}/caloris_test_record.json") for v in ("v3", "v3b")}; CO["rec"] = J(f"{HE}/peak_wall_test/caloris_test_record.json")
PW = {v: J(f"{G}/terrain/out_{v}/peak_wall_record.json") for v in ("v3", "v3b")}; PW["rec"] = J(f"{HE}/peak_wall_test/peak_wall_record.json")
TT = {v: J(f"{G}/terrain/out_{v}/terrain_table_record.json") for v in ("v3", "v3b")}; TT["rec"] = J(f"{PR}/analysis_v3/terrain_table_record.json")
OE = {v: J(f"{G}/terrain/out_{v}/terrain_obs_exp.json") for v in ("v2", "v3", "v3b")}
PO = {v: J(f"{G}/melt_source/positions/out_{v}/positions_record.json") for v in ("v3", "v3b")}; PO["rec"] = J(f"{PR}/analysis_workshop/melt_source_tests/positions/positions_record.json")
KN = {v: J(f"{G}/melt_source/knobs/out_{v}/knob_tests_record.json")["test2_hollow_density"] for v in ("v3", "v3b")}; KN["rec"] = J(f"{PR}/analysis_workshop/melt_source_tests/knobs/knob_tests_record.json")["test2_hollow_density"]
V = ("rec", "v3", "v3b")
ci = lambda t: f"{t[0]:.2f} ({t[1]:.2f}-{t[2]:.2f})"
pct = lambda x: f"{100 * x:.1f}%"
def fmtp(p): return f"{p:.1e}" if p < 0.001 else f"{p:.3f}"
rows = []
def add(sec, q, text, vals, script, change):
    rows.append(dict(section=sec, quantity=q, text_value=text, v2_record=vals[0], v3=vals[1], v3b=vals[2], script=script, conclusion_change=change))
# ---------------- catalog counts
c = lambda k: [f"{CC[v][k]:,}" for v in ("v2", "v3", "v3b")]
add("2.1; S1", "hollows after merging repeat detections (incl. the 64 v2 traces > 10 km)", "6,731", [x + (" (recount)" if i == 0 else "") for i, x in enumerate(c("hollows"))], "counts/catalog_counts.py", "no")
add("2.1; S1", "hollow fields (20 km chains): existing field_id with >= 1 hollow [rebuilt chains]", "1,228",
    [f"{CC[v]['fields_existing_ids']:,} [{CC[v]['fields_relinked_20km']:,}]" for v in ("v2", "v3", "v3b")], "counts/catalog_counts.py", "no")
add("2.1; S1", "released detections / images (Bickel et al. release)", "19,110 / 4,915", ["19,110 / 4,915"] * 3, "none (release counts; unchanged)", "no")
add("2.1; S2", "hollow-image pairs in the detection model", "92,382", [f"{DF[v]['theta']['n_pairs']:,}" if v != "v3" and v != "v3b" else f"{CC[v]['pairs']:,}" for v in ("v2", "v3", "v3b")], "prep_inputs.py; counts/catalog_counts.py", "no")
add("S2", "detected pairs / non-trace hollows in the detection fit", "10,669 / 6,667", [f"{CC[v]['pair_detections']:,} / {CC[v]['pair_hollows']:,}" for v in ("v2", "v3", "v3b")], "counts/catalog_counts.py", "no")
th = lambda t: f"{t['p_max']:.3f}, {t['x50_px']:.1f}, {t['sigma_lo']:.3f}, {t['x_hi_px']:.1f}, {t['sigma_hi']:.3f}, {t['b_inc']:.3f}"
add("S2", "detection model p_max, x50, sigma_lo, x_hi, sigma_hi, b (refit on the v3/v3b pairs)", "0.243, 35.5, 0.590, 195.9, 0.267, -0.102",
    [th(SF["rec"]["theta"]), th(DR["v3"]["theta"]), th(DR["v3b"]["theta"])], "sfd/refit_detection.py (G3_THETA=refit)", "no (p_max rises 12-14%, so completeness at every site would rise if the grids were rebuilt)")
add("S9 (color)", "cataloged hollows masked in the color map", "6,731", c("hollows"), "not rerun (color_mdr; joint-model scope)", "not rerun")
add("2.2 (other group)", "host craters (any hollow) / within 0.9 R", "707 / 606 (graph) ", ["707 / 608", "534 / 448", "441 / 371"], "out_v3/README.md", "see g1/g2")
# ---------------- size distribution
s = lambda k: [ci(SN[v][k]) for v in V]
add("3.1; S2", "exponent of N(>D), 0.3-0.7 km", "1.57 (1.38-1.76)", s("slope_0.3_0.7"), "sfd/sfd_numbers.py", "no")
add("3.1; S2", "exponent of N(>D), 0.7-3 km", "1.38 (1.30-1.46)", s("slope_0.7_3"), "sfd/sfd_numbers.py", "no")
add("3.1; S2", "exponent of N(>D), 1-3 km", "1.43 (1.31-1.55)", s("slope_1_3"), "sfd/sfd_numbers.py", "no")
add("3.1; Fig 1", "exponent range 0.3-3 km (D^-a)", "1.4 to 1.6", [f"{DF[v]['rounded']['range'][0]}-{DF[v]['rounded']['range'][1]}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1; S2", "exponent above 3 km (3-12 km)", "2.7 (2.4-2.9); S2 2.68 (2.43-2.94)", s("slope_3_12"), "sfd/sfd_numbers.py", "no (steeper: the deficit of large hollows is stronger)")
W = {v: J(f"{G}/sfd/out_{v}_fixed/sfd_windows.json")["windows"]["2.0-6.0"] for v in ("v3", "v3b")}; W["rec"] = J(f"{PR}/analysis_workshop/completeness_sfd/sfd_windows.json")["windows"]["2.0-6.0"]
add("3.1 (support)", "exponent 2-6 km (where the steepening starts)", "not quoted", [ci(W[v]) for v in V], "sfd/sfd_windows.py", "no (break moves from ~3 km toward ~2 km)")
cnt = lambda k: [f"{SN[v]['counts'][k]['corrected']:,} ({SN[v]['counts'][k]['corrected_ci'][0]:,}-{SN[v]['counts'][k]['corrected_ci'][1]:,})" for v in V]
add("3.1", "corrected N(>5 km)", "214 (182-246)", cnt("5.0"), "sfd/sfd_numbers.py", "no")
add("3.1", "power-law N(>5 km) if the 0.3-3 km law continued", "430-590", [f"{DF[v]['rounded']['N_gt_5km_powerlaw'][0]:.0f}-{DF[v]['rounded']['N_gt_5km_powerlaw'][1]:.0f}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no (large hollows still far rarer than the power law: 0.24-0.28 of it in v3, 0.36-0.50 in v2)")
add("3.1", "median detection probability at 3 and 5 km", "82-86%", [f"{pct(SN[v]['median_C_at_size']['3'])}, {pct(SN[v]['median_C_at_size']['5'])}" for v in V], "sfd/sfd_numbers.py", "no")
add("3.1", "largest detection vs the one used, median excess, 3-5 and 5-10 km", "5-6%", [f"{pct(DF[v]['largest_detection_excess_median']['3-5'][0])}, {pct(DF[v]['largest_detection_excess_median']['5-10'][0])}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1", "largest detection vs the one used, median excess, 0.5-1, 1-2, 2-3 km", "7-10%", [", ".join(pct(DF[v]['largest_detection_excess_median'][k][0]) for k in ("0.5-1", "1-2", "2-3")) for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1; S2", "detections excluded as whole-crater outlines (> 10 km)", "64", ["64", "64 + 368 screened (B/W/O)", "64 + 600 screened (B/W/O/U)"], "out_v3/README.md", "no (text should say that screened detections are also excluded)")
add("2.1; 3.1", "median detection probability of a 1 km / 0.3 km hollow at cataloged sites", "47% / 4%", [f"{pct(SN[v]['median_C_at_size']['1'])} / {pct(SN[v]['median_C_at_size']['0.3'])}" for v in V], "sfd/sfd_numbers.py", "no (1 km is still where it reaches about one half)")
add("3.1", "share of covering images in which 0.3 km spans < 10 px; median pixel scale", "80%; 63 m", [f"{pct(DF[v]['pairs_0p3km_lt_10px'])}; {DF[v]['pairs_median_px_m']} m" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1; S2", "other covering images that detected a cataloged hollow", "5% (S2: 4.7%)", [pct(DF[v]['pairs_other_images_detected']) for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1; S2", "same, where the hollow spans 60-100 / 100-200 px (one detection dropped at random)", "14-15%", [" / ".join(pct(x) for x in DF[v]['pairs_detected_Dpx_60_100_100_200']) for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1", "cataloged hollows > 0.3 km", "6,044", [f"{SN[v]['counts']['0.3']['catalog']:,}" for v in V], "sfd/sfd_numbers.py", "no")
add("3.1", "corrected N(>0.3 km)", "24,100 (23,200-25,000)", cnt("0.3"), "sfd/sfd_numbers.py", "no")
add("3.1", "power-law N(>0.3 km) from 1 km", "30,000-38,000", [f"{DF[v]['rounded']['N_gt_0p3km_powerlaw'][0]:,.0f}-{DF[v]['rounded']['N_gt_0p3km_powerlaw'][1]:,.0f}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1", "catalog share of hollows > 0.3 km", "one in five to one in six", [f"1 in {1/DF[v]['rounded']['catalog_fraction_0p3'][0]:.1f} to 1 in {1/DF[v]['rounded']['catalog_fraction_0p3'][1]:.1f}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("3.1 (support)", "corrected N(>1 km) (anchor of the power-law lines)", "not quoted", [f"{SN[v]['counts']['1.0']['corrected']:,}" for v in V], "sfd/sfd_numbers.py", "no")
add("3.1; S2", "best-imaged sites (C(0.3 km) >= 0.3): hollows, hollows > 0.3 km, exponent 0.3-1 km", "193; 106; 2.4 (1.8-3.1)", [f"{SN[v]['sites_C03_ge_0.3']['n']}; {SN[v]['sites_C03_ge_0.3']['n_ge03']}; {ci(SN[v]['sites_C03_ge_0.3']['slope_03_1'])}" for v in V], "sfd/sfd_numbers.py", "no")
add("S2", "exponent above lower limits 0.2, 0.3, 0.5 km", "1.55, 1.63, 1.62", [", ".join(f"{SF[v]['power'][k]['slope_cumulative']:.2f}" for k in ("0.2", "0.3", "0.5")) for v in V], "sfd/fit_sfd.py", "no")
add("S2", "exponent used for the size-weighted completeness E (D^-a)", "1.6 (1.55-1.63)", ["1.6", "1.7 (1.61-1.72)", "1.7-1.8 (1.63-1.76)"], "sfd/fit_sfd.py", "no (S2 says exponents 1.3-1.9 change the hosting odds ratios by < 5%; E-based fits not rerun here)")
add("S2", "lognormal fit median (km), sigma", "not quoted", [f"{SF[v]['lognormal']['median_km']:.3g}, {SF[v]['lognormal']['sigma_lnD']:.2f}" for v in V], "sfd/fit_sfd.py", "no")
add("4.4", "HRIC: hollows > 0.1 km per hollow > 0.3 km", "4.7-5.8 times", [f"{DF[v]['rounded']['ratio_0p1_over_0p3'][0]:.1f}-{DF[v]['rounded']['ratio_0p1_over_0p3'][1]:.1f}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("4.4", "hollows > 0.1 km on the searched surface", "140,000-220,000", [f"{DF[v]['rounded']['N_gt_0p1km'][0]:,.0f}-{DF[v]['rounded']['N_gt_0p1km'][1]:,.0f}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no")
add("4.4", "hollows > 5 km: corrected total against the catalog", "about 210 against 185", [f"{DF[v]['N_gt_5km']['corrected']} against {DF[v]['N_gt_5km']['catalog']}" for v in ("v2", "v3", "v3b")], "sfd/sfd_derived.py", "no (prediction halves)")
for v in ("v3", "v3b"): pass
add("3.1 sensitivity", "REFIT detection model: exponents 0.3-0.7 / 0.7-3 / 1-3 km; above 3 km", "1.57 / 1.38 / 1.43; 2.68",
    ["1.57 / 1.38 / 1.43; 2.68"] + [f"{' / '.join(str(x) for x in DR[v]['exponents_0.3-0.7_0.7-3_1-3'])}; {DR[v]['exponent_3_12'][0]} ({DR[v]['exponent_3_12'][1]}-{DR[v]['exponent_3_12'][2]})" for v in ("v3", "v3b")], "sfd/refit_detection.py + fit_sfd.py + sfd_numbers.py (G3_THETA=refit)", "no")
add("3.1 sensitivity", "REFIT: corrected N(>5 km); N(>0.3 km); median C at 1 km", "214; 24,100; 47%",
    ["214; 24,071; 47.3%"] + [f"{DR[v]['N_gt_5km']['corrected']} ({DR[v]['N_gt_5km']['ci'][0]}-{DR[v]['N_gt_5km']['ci'][1]}); {DR[v]['N_gt_0p3km']['corrected']:,}; {pct(DR[v]['median_C_at_size']['1'])}" for v in ("v3", "v3b")], "same (G3_THETA=refit)", "no")
add("4.4 sensitivity", "REFIT: ratio > 0.1 / > 0.3 km; N(> 0.1 km)", "4.7-5.8; 140,000-220,000",
    ["4.7-5.8; 140,000-223,000"] + [f"{DR[v]['rounded']['ratio_0p1_over_0p3'][0]:.1f}-{DR[v]['rounded']['ratio_0p1_over_0p3'][1]:.1f}; {DR[v]['rounded']['N_gt_0p1km'][0]:,.0f}-{DR[v]['rounded']['N_gt_0p1km'][1]:,.0f}" for v in ("v3", "v3b")], "sfd/sfd_derived.py (G3_THETA=refit)", "no")
add("closure (not in v6.3)", "size_law.py differential exponent, north, 0.1-10 km (cumulative = alpha - 1)", "not in v6.3",
    [f"{SL[v]['alpha_point']:.3f} ({SL[v]['alpha_95CI'][0]:.3f}-{SL[v]['alpha_95CI'][1]:.3f}); n = {SL[v]['n_fitted']}" for v in V], "sfd/size_law.py", "no")
# ---------------- Caloris
add("3.2; S7; Fig 2c", "Caloris craters with excavation < 2.5 km hosting >= 3 hollows", "0 of 144", [CA[v]["a_step"]["shallow"] for v in V], "caloris/caloris_v9.py", "no")
add("3.2; S7; Fig 2c", "Caloris craters with excavation >= 2.5 km hosting >= 3 hollows; Fisher p", "14 of 35; 1.0e-11", [f"{CA[v]['a_step']['deep']}; {fmtp(CA[v]['a_step']['fisher_p'])}" for v in V], "caloris/caloris_v9.py", "no")
add("3.2; S7", "isolated craters 10-30 km, C >= 0.2, hosting >= 3", "2 of 762", [CA[v]["a_C1"]["iso_10_30_C02"] for v in V], "caloris/caloris_v9.py", "no")
add("S7", "Caloris craters 10-30 km, C >= 0.2: expected and observed with >= 3", "0.14 of 53; one (29 km)", [f"{CA[v]['a_C1']['expected_cal']} of {CA[v]['a_C1']['cal_10_30_C02'].split('/')[1]}; {CA[v]['a_C1']['cal_10_30_C02'].split('/')[0]}" for v in V], "caloris/caloris_v9.py", "no")
add("3.2; S7", "30-35 km: Caloris vs isolated with >= 3; p", "3 of 5 vs 1 of 52; p = 0.001", [f"{CA[v]['a_C1']['30-35']['cal']} vs {CA[v]['a_C1']['30-35']['iso']}; p = {fmtp(CA[v]['a_C1']['30-35']['fisher_p'])}" for v in V], "caloris/caloris_v9.py", "no")
add("S7", "35-40 km: Caloris vs isolated with >= 3", "3 of 6 and 5 of 44", [f"{CA[v]['a_C1']['35-40']['cal']} and {CA[v]['a_C1']['35-40']['iso']}" for v in V], "caloris/caloris_v9.py", "no")
add("S7", "20-40 km: outer 30% vs inner; p", "5 of 27 and 2 of 19; p = 0.68", [f"{CA[v]['a_edge_20_40']['outer_rho_ge0p7']} and {CA[v]['a_edge_20_40']['inner']}; p = {CA[v]['a_edge_20_40']['fisher_p']:.2f}" for v in V], "caloris/caloris_v9.py", "no")
sh = lambda v: CA[v]["a_shallow_with_1_2"]
add("3.2; S7", "shallower Caloris craters with one or two hollows (diameters; excavation depths)", "Seven (10.1-24.6 km; 0.8-2.2 km)",
    [f"{len(sh(v))} ({min(x['D'] for x in sh(v)):.1f}-{max(x['D'] for x in sh(v)):.1f} km; {min(x['d_ex'] for x in sh(v)):.1f}-{max(x['d_ex'] for x in sh(v)):.1f} km)" for v in V], "caloris/caloris_v9.py", "no")
add("S7", "small (excavation < 2.5 km) Caloris craters >= 5 km with any hollow: inside a larger crater vs basin parent; p", "5 of 12 vs 20 of 452; p = 0.0002",
    [f"{CA[v]['b']['any_all']['nested']} vs {CA[v]['b']['any_all']['direct']}; p = {fmtp(CA[v]['b']['any_all']['fisher_p'])}" for v in V], "caloris/caloris_v9.py", "no")
nh = lambda v: [x for x in CA[v]["b"]["nested_children"] if x["n_all"] > 0]
PD = {p["id"]: p["D"] for p in CA["rec"]["b"]["parents"]}
add("S7", "nested small hosts: number and parents", "five hosts in three parents 56-129 km",
    [f"{len(nh(v))} hosts in {len(set(x['parent_id'] for x in nh(v)))} parents {min(PD[x['parent_id']] for x in nh(v)):.0f}-{max(PD[x['parent_id']] for x in nh(v)):.0f} km" for v in V], "caloris/caloris_v9.py", "no")
add("S7", "basin-parent small hosts in the outer fifth of Caloris; p", "15 of the 20 (195 of 452 craters there); p = 0.005",
    [f"{CA[v]['b']['direct_hosts_outer_fifth']['hosts']} ({CA[v]['b']['direct_hosts_outer_fifth']['all_direct']}); p = {fmtp(CA[v]['b']['direct_hosts_outer_fifth']['fisher_p'])}" for v in V], "caloris/caloris_v9.py", "YES (p crosses 0.05)")
add("S7", "odds ratios with completeness and diameter controlled: nested; outer fifth (p)", "13 and 3.9",
    [f"{CT[v]['nested_OR']:.1f} (p = {fmtp(CT[v]['nested_logit']['nested'][2])}); {CT[v]['outer_OR']:.1f} (p = {fmtp(CT[v]['outer_logit_direct']['outer'][2])})" for v in V], "caloris/caloris_v9_controls.py", "YES (outer-fifth effect no longer significant)")
add("3.2", "'small craters with hollows inside Caloris lie mostly in the outer fifth of the basin ... or inside larger craters'", "qualitative",
    ["holds (15 of 20 outer, p = 0.005; 5 of 12 nested)", "nested part holds (3 of 12 vs 6 of 452); outer-fifth part not significant (4 of 6, p = 0.41)", "nested part holds (3 of 12 vs 4 of 452); outer-fifth part not significant (3 of 4, p = 0.32)"], "caloris/caloris_v9.py", "YES")
rb = lambda v: "; ".join(f"{x['D']} {x['iso']}" for x in CA[v]["a_C1"]["rates_C02"] if x["D"] in ("35-40", "50-70", "70-250"))
add("Fig 2c", "isolated-crater counts with >= 3 hollows that change (35-40, 50-70, 70-250 km bins)", "figure", [rb(v) for v in V], "caloris/caloris_v9.py", "no")
add("README only", "old peak_wall_test/caloris_test.py step counts", "not in v6.3", [str(CO[v]["caloris_counts"]) for v in V], "caloris/caloris_test.py", "no")
# ---------------- terrain
pw = lambda v, f, k="pooled_hollows", sub="all": PW[v][f][sub][k] if sub == "all" else PW[v][f][sub]["pooled_hollows"]
fr = lambda d: f"{d['ratio']:.2f} ({d['ci'][0]:.2f}-{d['ci'][1]:.2f})"
add("2.2; 3.2; S4", "host craters 30-250 km in the terrain test", "352", [str(PW[v]["rows_H250_r0.45"]["all"]["pooled_hollows"]["N"]) for v in V], "terrain/peak_wall.py", "no")
add("3.2; S4; Table S3", "peak/wall density, baseline (250 m, 0.45 R)", "6.7 (4.3-10.0)", [fr(pw(v, "rows_H250_r0.45")) for v in V], "terrain/peak_wall.py + analyze.py", "no")
add("3.2; S4; Table S3", "floor/wall density", "3.1 (2.1-4.4)", [fr(pw(v, "rows_H250_r0.45", "floor_vs_wall")) for v in V], "terrain/analyze.py", "no")
add("Table S3", "peak/wall, 150 m threshold", "5.7 (3.8-8.3)", [fr(pw(v, "rows_H150_r0.45")) for v in V], "terrain/analyze.py", "no")
add("Table S3", "peak/wall, 400 m threshold", "8.2 (5.5-12.2)", [fr(pw(v, "rows_H400_r0.45")) for v in V], "terrain/analyze.py", "no")
add("Table S3", "peak/wall, 0.35 R peak limit", "7.4 (4.9-10.9)", [fr(pw(v, "rows_H250_r0.35")) for v in V], "terrain/analyze.py", "no")
add("Table S3", "peak/wall, Mantel-Haenszel by slope class", "6.5", [f"{PW[v]['rows_H250_r0.45']['all']['slope_stratified']['MH_ratio']:.2f}" for v in V], "terrain/analyze.py", "no")
add("Table S3; S4", "peak/wall without the three craters richest in hollows", "5.4 (3.6-8.4)", [fr(PW[v]["rows_H250_r0.45"]["drop_top3"]) for v in V], "terrain/analyze.py", "no")
add("Table S3; S4", "peak/wall, hollow fields counted", "3.6 (2.9-4.6)", [fr(pw(v, "rows_H250_r0.45", "pooled_fields")) for v in V], "terrain/analyze.py", "no")
add("Table S3", "equator-facing / pole-facing walls", "1.52 (0.95-2.47)", [fr(pw(v, "rows_H250_r0.45", "eqwall_vs_polewall")) for v in V], "terrain/analyze.py", "no (interval still includes 1)")
rr = lambda v: [PW[v][f]["all"]["pooled_hollows"]["ratio"] for f in ("rows_H150_r0.45", "rows_H400_r0.45", "rows_H250_r0.35")] + [PW[v]["rows_H250_r0.45"]["drop_top3"]["ratio"]]
add("S4", "range of peak/wall under alternative definitions and without top 3", "5.4-8.2", [f"{min(rr(v)):.1f}-{max(rr(v)):.1f}" for v in V], "terrain/analyze.py", "no")
for sname, lab, txt in (("all", "All host craters", "352: 18.2 (13.3-23.4), 8.4 (6.3-10.9), 2.7 (2.0-3.8), 3.3 (2.1-5.0), 2.2 (1.5-2.9)"),
                        ("isolated", "Isolated", "211: 21.0 (14.7-27.7), 7.2 (4.9-10.1), 1.8 (1.2-2.6), 1.6 (1.0-2.5), 1.9 (1.1-3.0)"),
                        ("basin", "Basin-hosted", "73: 12.0 (5.4-19.7), 8.8 (4.0-14.7), 2.7 (1.4-4.5), 3.4 (1.4-6.3), 1.9 (1.2-2.8)"),
                        ("nested", "Nested", "68: 7.4 (3.0-14.2), 13.2 (8.3-20.5), 7.7 (3.4-13.9), 11.4 (4.0-22.4), 4.0 (2.3-6.7)")):
    add("Table S2; Fig 2d", f"{lab}: N; density per 1,000 km2 on peak, floor, wall, equator-facing wall, pole-facing wall", txt,
        [f"{TT[v][sname]['N']}: " + ", ".join(f"{TT[v][sname][t]['density']:.1f} ({TT[v][sname][t]['lo']:.1f}-{TT[v][sname][t]['hi']:.1f})" for t in ("peak", "floor", "wall", "wall_eq", "wall_pole")) for v in V],
        "terrain/terrain_table.py", "no")
oe = lambda v: ", ".join(f"{OE[v]['obs_over_exp'][t]['ratio']:.2f} ({OE[v]['obs_over_exp'][t]['lo']:.2f}-{OE[v]['obs_over_exp'][t]['hi']:.2f})" for t in ("peak", "floor", "wall"))
add("Table S2", "observed/expected hollows, peak, floor, wall (re-derived; v2 column is this script's v2 run)", "2.98 (2.45-3.60), 1.32 (1.17-1.48), 0.45 (0.34-0.58)", [oe("v2"), oe("v3"), oe("v3b")], "terrain/terrain_obs_exp.py", "no")
add("Table S2", "density over all terrain used for shading (re-derived)", "6.1", [f"{OE[v]['all_terrain_density_per_1000km2']:.1f}" for v in ("v2", "v3", "v3b")], "terrain/terrain_obs_exp.py", "no")
add("S4", "isolated craters: peak/wall", "12 times", [fr(PW[v]["rows_H250_r0.45"]["isolated"]["pooled_hollows"]) for v in V], "terrain/analyze.py", "no")
add("S4", "nested walls / isolated walls; densest terrain in nested craters", "4.4 times; floors", [f"{OE[v]['nested_wall_over_isolated_wall']:.1f}; {OE[v]['nested_densest']}" for v in ("v2", "v3", "v3b")], "terrain/terrain_obs_exp.py", "no")
# ---------------- composition
add("S9; S16", "composition/LRM numbers that use the catalog", "S9 depth-vs-LRM and XRS use Blewett depths only",
    ["composition_v6/depth_lrm.py reads no catalog file", "not applicable", "not applicable"], "none", "no (LRM share in the joint model belongs to the joint-model rerun)")
# ---------------- positions & knobs (workshop)
ph = lambda v, k, grp="hollows": PO[v][grp][k]
sd = lambda d: f"{d['n_in']} of {d['n_points']}; SDR {d['SDR'][0]:.2f} ({d['SDR'][1]:.2f}-{d['SDR'][2]:.2f})"
for k, lab in (("peak_eta_low_vs_upper", "peak hollows in the lowest third of eta"), ("peak_base_within2km", "peak hollows within 2 km of the base"),
               ("floor_margin_within3km_of_wall_or_peak", "floor hollows within 3 km of wall or peak"), ("floor_within3km_of_wall", "floor hollows within 3 km of the wall"),
               ("floor_within3km_of_peak", "floor hollows within 3 km of the peak")):
    add("workshop: melt-source positions", lab, "not in v6.3", [sd(ph(v, k)) for v in V], "melt_source/positions/positions_build.py + positions_analyze.py", "no")
add("workshop: melt-source positions", "floor FIELDS within 3 km of wall or peak", "not in v6.3", [sd(ph(v, "floor_margin_within3km_of_wall_or_peak", "fields")) for v in V], "melt_source/positions/positions_analyze.py", "yes (interval newly includes 1; workshop result only)")
add("workshop: melt-source positions", "mean height quantile q of peak hollows", "not in v6.3", [f"{PO[v]['mean_q_peak_hollows'][0]:.2f} ({PO[v]['mean_q_peak_hollows'][1]:.2f}-{PO[v]['mean_q_peak_hollows'][2]:.2f})" for v in V], "melt_source/positions/positions_analyze.py", "no")
kr = lambda d: f"{d['hollows_on_knobs']} on knobs, {d['hollows_matrix']} matrix; ratio {d['ratio_knob_over_matrix']:.2f}; knob boot {d['ratio_CI_knob_bootstrap'][0]:.2f}-{d['ratio_CI_knob_bootstrap'][1]:.1f}; block {d['ratio_CI_block_bootstrap_3deg'][0]:.2f}-{d['ratio_CI_block_bootstrap_3deg'][1]:.1f}; exact {d['ratio_CI_exact_conditional'][0]:.2f}-{d['ratio_CI_exact_conditional'][1]:.1f}"
add("workshop: knobs", "all 2,011 knobs: knob/matrix hollow density", "not in v6.3", [kr(KN[v]["knobs"]["all_2011"]) for v in V], "melt_source/knobs/t2_knob_hollows.py + summarize.py", "no")
add("workshop: knobs", "335 knobs >= 5 km: knob/matrix hollow density", "not in v6.3", [kr(KN[v]["knobs"]["main_certain_uncertain_long_axis_ge_5km"]) for v in V], "melt_source/knobs/t2_knob_hollows.py + summarize.py", "yes for v3b only (3-degree block-bootstrap interval newly includes 1; knob-bootstrap and exact intervals exclude 1)")
pp = lambda d: f"{d['pk_density_per1000km2_w'][0]:.2f} ({d['pk_density_per1000km2_w'][1]:.2f}-{d['pk_density_per1000km2_w'][2]:.2f}); peak/floor {d['ratio_pk_over_fl'][0]:.2f}; peak/wall {d['ratio_pk_over_wl'][0]:.2f}"
add("workshop: knobs", "peak proxy, all craters 30-250 km: density; peak/floor; peak/wall", "not in v6.3", [pp(KN[v]["crater_peak_proxy"]["all_30_250km"]) for v in V], "melt_source/knobs/t2_peak_hollows.py + summarize.py", "no")
with open(f"{G}/changes_g3.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["section", "quantity", "text_value", "v2_record", "v3", "v3b", "script", "conclusion_change"]); w.writeheader(); w.writerows(rows)
print(len(rows), "rows")
