# g13 revision (7 Nov 2026): G13_REV=r2 (default) reads out_traced_r2 and writes results_r2.json; r1 keeps out_traced and results.json.
#!/usr/bin/env python3
"""g13_traced/positions (6 Nov 2026), new script: collects the position and slope results with traced classes (out_traced/) next to the
g11 curvature results (g11_curvature/positions/results.json). Reuses the extraction functions of g11 collect_results.py (executed up to
its data-reading line). Writes results.json."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); G11 = os.path.join(HERE, "..", "..", "g11_curvature", "positions")
src = open(os.path.join(G11, "collect_results.py")).read(); src = src[:src.index("rc, rf = J(")]
ns = {"__file__": os.path.join(G11, "collect_results.py")}; exec(src, ns)
J, positions, floors, slope = ns["J"], ns["positions"], ns["floors"], ns["slope"]
REV = os.environ.get("G13_REV", "r2"); OT = f"{HERE}/out_traced" + ("" if REV == "r1" else "_r2")
rt = J(f"{OT}/positions_record.json"); sl = {e: J(f"{OT}/positions_slope_record_{e}.json") for e in ("eta", "eta_b")}
NEW = ("main_30cells_norelief", "sens_30cells_500m", "sens_30cells_750m"); g11 = J(f"{G11}/results.json")
out = {"created": "2026-11-06", "configuration": "v3 catalog and grid as g11; terrain classes from Michael's traced outlines for 62 hosts (g13_traced/traces/traces_clean.json, status done), curvature v3 (g10 curv_classes_v3) elsewhere",
       "definitions": g11["definitions"],
       "positions": {"traced_eta": positions(rt, "eta"), "traced_eta_b": positions(rt, "eta_b"), "g11_curv_eta": g11["positions"]["curv_eta"], "g11_curv_eta_b": g11["positions"]["curv_eta_b"]},
       "floors": {"traced": floors(rt), "g11_curv": g11["floors"]["curv"]},
       "slope": {"traced_eta": slope(sl["eta"], NEW), "traced_eta_b": slope(sl["eta_b"], NEW), "g11_curv_eta": g11["slope"]["curv_eta"], "g11_curv_eta_b": g11["slope"]["curv_eta_b"]}}
json.dump(out, open(f"{HERE}/results" + ("" if REV == "r1" else "_r2") + ".json", "w"), indent=1)
f = lambda s: f"{s['ratio'][0]:.2f} ({s['ratio'][1]:.2f}-{s['ratio'][2]:.2f}); {s['n_in_band']}" if s.get("ratio") and s["ratio"][0] is not None else "-"
fm = lambda t: f"{t[0]:.2f} ({t[1]:.2f}-{t[2]:.2f})" if t and t[0] is not None else "-"
P = out["positions"]; F = out["floors"]; S = out["slope"]
print("counts traced", P["traced_eta"]["counts"]); print("counts g11", P["g11_curv_eta"]["counts"])
for k in ("peak_lowest_third_vs_upper_two_thirds", "peak_highest_third_vs_lower_two_thirds", "peak_within_2km_of_floor_vs_rest_of_peak"):
    for w in ("hollows", "fields"):
        print(f"{k} {w}: g11 eta {f(P['g11_curv_eta'][k][w])} | eta_b {f(P['g11_curv_eta_b'][k][w])} || traced eta {f(P['traced_eta'][k][w])} | eta_b {f(P['traced_eta_b'][k][w])}")
for m in ("central_peak", "peak_ring"):
    for w in ("hollows", "fields"):
        a = P["traced_eta"]["breakdown_by_morphology"][m][w]; b = P["g11_curv_eta"]["breakdown_by_morphology"][m][w]
        print(f"  {m} {w} low/high: g11 {f(b['peak_eta_low_vs_upper'])} / {f(b['peak_eta_high_vs_lower'])} || traced {f(a['peak_eta_low_vs_upper'])} / {f(a['peak_eta_high_vs_lower'])}")
print("mean q", P["g11_curv_eta"]["mean_height_quantile_q_peak_hollows"], "->", P["traced_eta"]["mean_height_quantile_q_peak_hollows"])
for k in ("floor_within_0.2R_of_center",):
    print(k, "hollows vs rest: g11", f(F["g11_curv"][k]["hollows_vs_rest_of_floor"]), "traced", f(F["traced"][k]["hollows_vs_rest_of_floor"]), "| fields g11", f(F["g11_curv"][k]["fields_vs_rest_of_floor"]), "traced", f(F["traced"][k]["fields_vs_rest_of_floor"]))
for k in ("floor_within_0.05R_of_peak_alone", "floor_within_0.1R_of_peak_alone", "floor_within_0.05R_of_wall_alone", "floor_within_0.1R_of_wall_alone", "floor_within_3km_of_peak_alone", "floor_within_3km_of_wall_alone", "floor_margin_0.1R_of_wall_or_peak"):
    print(k, "hollows g11", f(F["g11_curv"][k]["hollows"]), "| traced", f(F["traced"][k]["hollows"]), "|| fields g11", f(F["g11_curv"][k]["fields"]), "| traced", f(F["traced"][k]["fields"]))
print("floor OE by radius traced", F["traced"]["floor_OE_by_radius"]["edges"], F["traced"]["floor_OE_by_radius"]["OE"], F["traced"]["floor_OE_by_radius"]["n"])
print("floor OE by radius g11   ", F["g11_curv"]["floor_OE_by_radius"]["OE"], F["g11_curv"]["floor_OE_by_radius"]["n"])
for e in ("eta", "eta_b"):
    for sc in NEW:
        a = S[f"traced_{e}"]["screens"][sc]; b = S[f"g11_curv_{e}"]["screens"][sc]
        print(f"slope {e} {sc}: hosts {b['hosts_passing']}->{a['hosts_passing']}; low third slope-bins g11 {fm(b['hollows_lowest_third_vs_upper_two_thirds']['slope_bins'])} -> {fm(a['hollows_lowest_third_vs_upper_two_thirds']['slope_bins'])} (no control {fm(a['hollows_lowest_third_vs_upper_two_thirds']['no_slope_control'])}); "
              f"high third hollows g11 {fm(b['hollows_highest_third_vs_lower_two_thirds']['slope_bins'])} -> {fm(a['hollows_highest_third_vs_lower_two_thirds']['slope_bins'])}; fields high g11 {fm(b['fields_highest_third_vs_lower_two_thirds']['slope_bins'])} -> {fm(a['fields_highest_third_vs_lower_two_thirds']['slope_bins'])}; r {a['eta_slope_r_within_crater']}")
