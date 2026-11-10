#!/usr/bin/env python3
"""g9 (6 Nov 2026), new script: collects the numbers the draft quotes from every configuration's outputs into results.json.
Ratios are crater-matched (SDR) for positions, Mantel-Haenszel for slope, as in the draft; [estimate, lo, hi] with 95% intervals."""
import json, os
G9 = os.path.dirname(os.path.abspath(__file__)); WS = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/melt_source_tests")
J = lambda p: json.load(open(p))
rd = lambda v, k=3: None if v is None else round(v, k)
def tri(s, k="SDR"): return [rd(x) for x in s[k]]
def split(s): return dict(ratio=tri(s), pooled_DR=tri(s, "DR"), n=s["n_points"], n_in_band=s["n_in"], **({"craters_with_points": s["n_craters_with_points"]} if "n_craters_with_points" in s else {}))
POS = {"v2_record": ("positions/out_v2_old_gridv2", "v2 catalog, old classes, v2 grid: rerun of the workshop cache, identical to analysis_workshop/.../positions/positions_record.json"),
       "v3_old_classes_gridv2_g3": ("positions/out_v3_old_gridv2", "v3 catalog, old classes, v2 grid: rerun of the g3 cache, identical to g3_catalog/melt_source/positions/out_v3/positions_record.json"),
       "v3_old_classes_gridv3": ("positions/out_v3_old_gridv3", "v3 catalog, old classes, v3 refit grid (separates the grid change from the class fix)"),
       "v3r_fix": ("positions/out_v3_fix_gridv3", "v3 catalog, g7 classes, v3 refit grid (current configuration)")}
out = {"created": "2026-11-06", "folder": G9,
       "configuration_current": "v3r + terrain fix: hollows_canonical_v3.csv (non-trace), host_assignments_v3.csv, g4 completeness_10ppd_v3.npz (Dref 1 km), g7 terrain_classes.classify",
       "positions": {}, "slope": {}, "knobs": {}, "size_transition": {}}
for k, (d, note) in POS.items():
    r = J(f"{G9}/{d}/positions_record.json"); H, F, X = r["hollows"], r["fields"], r["g9_additions"]; fr = r["floor_density_by_radius"]
    out["positions"][k] = {"note": note, "sample": r["sample"],
        "peak_lowest_third_vs_upper_two_thirds": {"hollows": split(H["peak_eta_low_vs_upper"]), "fields": split(F["peak_eta_low_vs_upper"])},
        "peak_highest_third_vs_lower_two_thirds": {"hollows": split(H["peak_eta_high_vs_lower"]), "fields": split(F["peak_eta_high_vs_lower"])},
        "peak_within_2km_of_floor_vs_rest_of_peak": {"hollows": split(H["peak_base_within2km"]), "fields": split(F["peak_base_within2km"])},
        "mean_height_quantile_q_peak_hollows": [rd(x) for x in r["mean_q_peak_hollows"]],
        "floor_within_0.2R_of_center": {"OE_vs_whole_floor_as_quoted_in_draft": dict(OE=[rd(fr["OE"][0]), rd(fr["lo"][0]), rd(fr["hi"][0])], n_in_band=fr["n"][0], n=sum(fr["n"])),
                                        "hollows_vs_rest_of_floor": split(X["hollows"]["floor_within0.2R_of_center_vs_rest"]),
                                        "fields_vs_rest_of_floor": split(X["fields"]["floor_within0.2R_of_center_vs_rest"])},
        "floor_OE_by_radius": dict(edges=fr["edges"], n=fr["n"], OE=[rd(x) for x in fr["OE"]], lo=[rd(x) for x in fr["lo"]], hi=[rd(x) for x in fr["hi"]]),
        "floor_margin_3km_of_wall_or_peak": {"hollows": split(H["floor_margin_within3km_of_wall_or_peak"]), "fields": split(F["floor_margin_within3km_of_wall_or_peak"]),
                                             "hollows_radius_stratified": dict(ratio=tri(X["hollows"]["floor_within3km_of_wall_or_peak_radius_stratified"])),
                                             "fields_radius_stratified": dict(ratio=tri(X["fields"]["floor_within3km_of_wall_or_peak_radius_stratified"]))},
        "floor_margin_0.05R_of_wall_or_peak": {"hollows": split(X["hollows"]["floor_within0.05R_of_wall_or_peak"]), "fields": split(X["fields"]["floor_within0.05R_of_wall_or_peak"]),
                                               "hollows_radius_stratified": dict(ratio=tri(X["hollows"]["floor_within0.05R_of_wall_or_peak_radius_stratified"])),
                                               "fields_radius_stratified": dict(ratio=tri(X["fields"]["floor_within0.05R_of_wall_or_peak_radius_stratified"]))},
        "floor_margin_0.1R_of_wall_or_peak": {"hollows": split(X["hollows"]["floor_within0.1R_of_wall_or_peak"]), "fields": split(X["fields"]["floor_within0.1R_of_wall_or_peak"]),
                                              "hollows_radius_stratified": dict(ratio=tri(X["hollows"]["floor_within0.1R_of_wall_or_peak_radius_stratified"])),
                                              "fields_radius_stratified": dict(ratio=tri(X["fields"]["floor_within0.1R_of_wall_or_peak_radius_stratified"]))},
        "floor_within_3km_of_wall": split(H["floor_within3km_of_wall"]), "floor_within_3km_of_peak": split(H["floor_within3km_of_peak"])}
SL = {"v2_record": "positions/out_v2_old_gridv2", "v3_old_classes_gridv2": "positions/out_v3_old_gridv2", "v3_old_classes_gridv3": "positions/out_v3_old_gridv3", "v3r_fix": "positions/out_v3_fix_gridv3"}
mh = lambda m: [rd(m["RR"]), rd(m["lo"]), rd(m["hi"])]
for k, d in SL.items():
    r = J(f"{G9}/{d}/positions_slope_record.json"); o = {"screens": {}}
    for sc in ("loose_15cells_500m", "main_30cells_750m", "strict_60cells_1000m"):
        Hh, Ff = r["results"][f"s1|{sc}|hollows"], r["results"][f"s1|{sc}|fields"]
        o["screens"][sc] = dict(hosts_passing=r["screens"][sc]["craters"], hosts_with_peak_hollows=r["screens"][sc]["craters_with_peak_hollows"],
                                peak_hollows=r["screens"][sc]["peak_hollows"], peak_fields=r["screens"][sc]["peak_fields"],
                                hollows_lowest_third_vs_upper_two_thirds={"no_slope_control": mh(Hh["MH_lowest_vs_upper2_crater_only"]), "slope_bins": mh(Hh["MH_lowest_vs_upper2_by_sf"]), "slope_terciles": mh(Hh["MH_lowest_vs_upper2_by_st"])},
                                hollows_highest_third_vs_lower_two_thirds={"no_slope_control": mh(Hh["MH_upper_vs_lower2_crater_only"]), "slope_bins": mh(Hh["MH_upper_vs_lower2_by_sf"]), "slope_terciles": mh(Hh["MH_upper_vs_lower2_by_st"])},
                                fields_highest_third_vs_lower_two_thirds={"no_slope_control": mh(Ff["MH_upper_vs_lower2_crater_only"]), "slope_bins": mh(Ff["MH_upper_vs_lower2_by_sf"]), "slope_terciles": mh(Ff["MH_upper_vs_lower2_by_st"])})
    c = r["results"]["s1|main_30cells_750m|hollows"]["eta_slope_correlation"]
    o["eta_slope_correlation_main_screen_1cell"] = {kk: (rd(v) if isinstance(v, float) else v) for kk, v in c.items()}
    o["note"] = "draft quotes MH with slope bins (0-1, 1-2, 2-4, 4-8, > 8 deg), 1-cell slope; " + ("v3 old classes with the v2 grid was not run in g3; computed here" if k == "v3_old_classes_gridv2" else "")
    out["slope"][k] = o
KN = {"v2_record": f"{WS}/knobs/knob_tests_record.json", "v3_gridv2_g3": f"{G9}/knobs/out_v3_gridv2/knob_tests_record.json", "v3r": f"{G9}/knobs/out_v3_gridv3/knob_tests_record.json"}
for k, p in KN.items():
    t = J(p)["test2_hollow_density"]; a = t["knobs"]["all_2011"]; m = t["knobs"]["main_certain_uncertain_long_axis_ge_5km"]
    out["knobs"][k] = dict(knobs_with_hollows=a["knobs_with_hollows"], hollows_on_knobs=a["hollows_on_knobs"], hollow_fields_on_knobs=a["hollow_fields_on_knobs"],
        knob_density_per1000km2_w=a["density_knobs_per1000km2_w"], matrix_density_per1000km2_w=a["density_matrix_per1000km2_w"],
        knob_over_matrix_1p5_to_3_radii=dict(ratio=a["ratio_knob_over_matrix"], ci_knob_bootstrap=a["ratio_CI_knob_bootstrap"], ci_block_bootstrap_3deg=a["ratio_CI_block_bootstrap_3deg"], ci_exact_conditional=a["ratio_CI_exact_conditional"]),
        knob_over_kinczyk_class1_2_peak_proxy_0p15R=t["knob_density_over_peak_density"]["all_2011__over__Kinczyk_class1_2"],
        knobs_ge_5km=dict(knobs_with_hollows=m["knobs_with_hollows"], hollows=m["hollows_on_knobs"], ratio=m["ratio_knob_over_matrix"], ci_knob_bootstrap=m["ratio_CI_knob_bootstrap"],
                          ci_block_bootstrap_3deg=m["ratio_CI_block_bootstrap_3deg"], ci_exact_conditional=m["ratio_CI_exact_conditional"]),
        source=p)
ST = {"v2_record": "size_transition/out_v2", "v3_gridv2": "size_transition/out_v3", "v3r": "size_transition/out_v3r"}
for k, d in ST.items():
    r = J(f"{G9}/{d}/size_transition_record.json"); st2, st3, st4 = r["st2_breaks"], r["st3_nested"], r["st4_morphology"]; iso = st2["isolated|5-60|quad|C1"]; f12 = iso["fixed_12km"]["step"]
    A = st4["A_12-60_peak_vs_nopeak(f+s)"]; ck = J(f"{G9}/{d}/cache/isolated_small_check.json")
    out["size_transition"][k] = dict(
        isolated_5_60km=dict(craters=iso["profile"]["n"], hosts=iso["profile"]["hosts"], step_OR_at_12km_cluster95=[f12["OR"]] + f12["OR_CI95_cluster"], p_LR=f12["p_LR"],
                             best_step=iso["profile"]["best_step"], calibrated_p=r["summary"]["P1_isolated_rise_at_transition"]["calibrated_p_isolated"],
                             step_OR_variants={kk: [st2[kk]["fixed_12km"]["step"]["OR"]] + st2[kk]["fixed_12km"]["step"]["OR_CI95_cluster"] for kk in
                                               ("isolated|5-60|spline|C1", "isolated|5-60|quad|logitC", "isolated|5-100|quad|C1", "isolated|10-60|quad|C1", "isolated|5-60|quad|C1|C1>=0.1")},
                             blurred_transition_OR=r["st2c_blurred"]["isolated"]["OR_fully_complex_vs_simple_[est,lo,hi]_cluster"]),
        nested_step_OR_at_12km=[st2["nested|5-60|quad|C1"]["fixed_12km"]["step"]["OR"]] + st2["nested|5-60|quad|C1"]["fixed_12km"]["step"]["OR_CI95_cluster"],
        peak_vs_no_peak_unfilled_12_60km=dict(OR=A["OR_peak_[est,lo,hi]_cluster"], craters=A["n"], hosts=A["hosts"], peak_craters=A["n_peak"], peak_hosts=A["hosts_peak"]),
        complex_vs_simple_10_20km=st4["B_10-20_complex_vs_simple"]["OR_complex_[est,lo,hi]_cluster"],
        hosts_5_12km=dict(isolated=ck["hosts_5_12_isolated"], isolated_crater_sized_only=ck["hosts_5_12_isolated_crater_sized_only"], nested=ck["hosts_5_12_nested"],
                          nested_crater_sized_only=ck["hosts_5_12_nested_crater_sized_only"], MH_OR_nested_5_12=[st3["5-12"]["MH_C1_strata"]["OR"]] + st3["5-12"]["MH_C1_strata"]["CI95"],
                          MH_OR_nested_12_30=[st3["12-30"]["MH_C1_strata"]["OR"]] + st3["12-30"]["MH_C1_strata"]["CI95"]),
        features=r["inputs"]["features"])
out["relief_750m_rationale"] = ("Peak cells are defined as ground more than H = 250 m above the floor level, so eta = h/(z_peak95 - z_floor) cannot fall below "
    "250/(z_peak95 - z_floor). The lowest third of eta (eta < 1/3) exists only where z_peak95 - z_floor > 3 H = 750 m (positions README, caveats). The slope README adds "
    "that small peaks are the most smoothed by the 665 m DEM; the >= 30 peak-cell part of the screen addresses that. positions_slope.py itself gives no other reason "
    "(its docstring calls the screen a resolution screen). The g7 classes keep h > H for peak cells, so the bound still holds.")
out["reproduction_checks"] = {
    "positions_build": "v3/old classes/v2 grid, hosts 0-40: every array identical to g3 cache H250_000_354.npz for those craters (_scratch/test_chunks)",
    "positions_analyze": "on the g3 v3 cache: identical to g3 out_v3/positions_record.json; on the workshop v2 cache: identical to the workshop positions_record.json (0 differences, tolerance 1e-9)",
    "positions_slope_build + positions_slope": "v2/old classes/v2 grid: slope_H250.npz identical to the workshop cache; positions_slope_record.json identical (0 differences)",
    "knobs t2_knob_hollows + t2_peak_hollows + summarize": "v3/v2 grid: knob_hollows.csv, region_hollows.json and the peak chunk byte-identical to g3 out_v3; knob_tests_record.json identical",
    "size_transition st1-st6": "v2: size_transition_record.json identical to the workshop record except the fields created, folder and inputs",
    "isolated_small_check": "v2: 29 isolated hosts 5-12 km, 16 crater-sized only; 128 nested, 85 crater-sized only (as README_isolated_small.md)"}
out["notes"] = {
    "d_nonfloor": "kept as defined (distance to the nearest unmasked non-floor cell within 1.05 R). Under the g7 classes the floor reaches 1.0 R, so the radial cut it contained moves from 0.8 R to the rim (1.0 R): non-floor cells are peak, wall and the 1.0-1.05 R ring. It is not among the requested measures; the floor-margin measures use d_wall and d_peak only.",
    "radial_bins": "under the g7 classes the radius bins of floor_density_by_radius and of the radius-stratified tests gain a 0.8-1.0 R bin (old: 0-0.2-0.4-0.6-0.81)",
    "eta_in_annular_floor_craters": "eta is still measured from the inner floor level F0 (classify returns h = z - F0). In the 35 peak-ring craters whose annular floor is resolved, 91 of 411 peak hollows and 26% of peak cells lie beyond the ring crest, where the local base is the annular floor (_scratch/ring_outer_flank.py)",
    "knobs": "neither knob script uses terrain classes; t2_peak_hollows uses a geometric peak proxy (rho <= 0.15 R), so classify was not applied"}
json.dump(out, open(f"{G9}/results.json", "w"), indent=1)
print("wrote results.json")
