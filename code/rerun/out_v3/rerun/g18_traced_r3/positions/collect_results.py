# g18 copy of g16_traced_only/positions/collect_results.py (7 Nov 2026). Changed: g18 caches (round-2 traces, 111 traced hosts) and the g18 traced-only rows; outputs in out_tracedonly_r3. In collect_results.py the comparison column ('g13r2' labels) is the g16 traced-only result (62 hosts).
# g16 copy of g13_traced/positions/collect_results.py (7 Nov 2026). Changed: reads g16_traced_only/positions/out_tracedonly_r2 (62 traced hosts only) and
# prints it next to g13 r2 (traced in 62 hosts + curvature v3 in the other 259); writes results_tracedonly_r2.json.
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); G11 = os.path.join(HERE, "..", "..", "g11_curvature", "positions"); G13 = os.path.join(HERE, "..", "..", "g16_traced_only", "positions")
src = open(os.path.join(G11, "collect_results.py")).read(); src = src[:src.index("rc, rf = J(")]
ns = {"__file__": os.path.join(G11, "collect_results.py")}; exec(src, ns)
J, positions, floors, slope = ns["J"], ns["positions"], ns["floors"], ns["slope"]
OT = f"{HERE}/out_tracedonly_r3"
rt = J(f"{OT}/positions_record.json"); sl = {e: J(f"{OT}/positions_slope_record_{e}.json") for e in ("eta", "eta_b")}
NEW = ("main_30cells_norelief", "sens_30cells_500m", "sens_30cells_750m"); g13 = J(f"{G13}/results_tracedonly_r2.json")
out = {"created": "2026-11-07", "configuration": "v3 catalog and grid; only the 62 hosts with Michael's traced outlines (r2: refined wall feet, slump as wall)",
       "positions": {"tracedonly_eta": positions(rt, "eta"), "tracedonly_eta_b": positions(rt, "eta_b"), "g13r2_eta": g13["positions"]["tracedonly_eta"], "g13r2_eta_b": g13["positions"]["tracedonly_eta_b"]},
       "floors": {"tracedonly": floors(rt), "g13r2": g13["floors"]["tracedonly"]},
       "slope": {"tracedonly_eta": slope(sl["eta"], NEW), "tracedonly_eta_b": slope(sl["eta_b"], NEW), "g13r2_eta": g13["slope"]["tracedonly_eta"], "g13r2_eta_b": g13["slope"]["tracedonly_eta_b"]}}
json.dump(out, open(f"{HERE}/results_tracedonly_r3.json", "w"), indent=1)
f = lambda s: f"{s['ratio'][0]:.2f} ({s['ratio'][1]:.2f}-{s['ratio'][2]:.2f}); {s['n_in_band']}" if s.get("ratio") and s["ratio"][0] is not None else "-"
fm = lambda t: f"{t[0]:.2f} ({t[1]:.2f}-{t[2]:.2f})" if t and t[0] is not None else "-"
P = out["positions"]; F = out["floors"]; S = out["slope"]
print("counts tracedonly", P["tracedonly_eta"]["counts"]); print("counts g13r2", P["g13r2_eta"]["counts"])
for k in ("peak_lowest_third_vs_upper_two_thirds", "peak_highest_third_vs_lower_two_thirds", "peak_within_2km_of_floor_vs_rest_of_peak"):
    for w in ("hollows", "fields"):
        print(f"{k} {w}: g13r2 eta {f(P['g13r2_eta'][k][w])} || tracedonly eta {f(P['tracedonly_eta'][k][w])} | eta_b {f(P['tracedonly_eta_b'][k][w])}")
for m in ("central_peak", "peak_ring"):
    for w in ("hollows", "fields"):
        a = P["tracedonly_eta"]["breakdown_by_morphology"][m][w]; b = P["g13r2_eta"]["breakdown_by_morphology"][m][w]
        print(f"  {m} {w} low/high: g13r2 {f(b['peak_eta_low_vs_upper'])} / {f(b['peak_eta_high_vs_lower'])} || tracedonly {f(a['peak_eta_low_vs_upper'])} / {f(a['peak_eta_high_vs_lower'])}")
for k in ("floor_within_0.2R_of_center",):
    print(k, "hollows: g13r2", f(F["g13r2"][k]["hollows_vs_rest_of_floor"]), "tracedonly", f(F["tracedonly"][k]["hollows_vs_rest_of_floor"]), "| fields g13r2", f(F["g13r2"][k]["fields_vs_rest_of_floor"]), "tracedonly", f(F["tracedonly"][k]["fields_vs_rest_of_floor"]))
for k in ("floor_within_0.05R_of_peak_alone", "floor_within_0.1R_of_peak_alone", "floor_within_0.05R_of_wall_alone", "floor_within_0.1R_of_wall_alone", "floor_within_3km_of_peak_alone", "floor_within_3km_of_wall_alone", "floor_margin_0.1R_of_wall_or_peak"):
    print(k, "hollows g13r2", f(F["g13r2"][k]["hollows"]), "| tracedonly", f(F["tracedonly"][k]["hollows"]), "|| fields g13r2", f(F["g13r2"][k]["fields"]), "| tracedonly", f(F["tracedonly"][k]["fields"]))
for e in ("eta", "eta_b"):
    for sc in NEW:
        a = S[f"tracedonly_{e}"]["screens"][sc]; b = S[f"g13r2_{e}"]["screens"][sc]
        print(f"slope {e} {sc}: hosts {b['hosts_passing']}->{a['hosts_passing']}; high third hollows {fm(b['hollows_highest_third_vs_lower_two_thirds']['slope_bins'])} -> {fm(a['hollows_highest_third_vs_lower_two_thirds']['slope_bins'])}; fields high {fm(b['fields_highest_third_vs_lower_two_thirds']['slope_bins'])} -> {fm(a['fields_highest_third_vs_lower_two_thirds']['slope_bins'])}; low third hollows {fm(b['hollows_lowest_third_vs_upper_two_thirds']['slope_bins'])} -> {fm(a['hollows_lowest_third_vs_upper_two_thirds']['slope_bins'])}")
