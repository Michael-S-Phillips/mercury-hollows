#!/usr/bin/env python3
"""g13 (6 Nov 2026), new script modeled on g11_curvature/terrain/build_summary.py: collects the terrain test with Michael's traced
outlines (62 hosts traced, curvature v3 elsewhere; terrain/out_v3/*traced_v3*) next to the g11 curvature v2 values
(g11_curvature/terrain/out_v3/*n24_h100_s0.03*), plus the decomposition runs (curvature v3 alone, traced + curvature v2), the
sensitivity without the 13 hosts whose tracing notes flag uncertainty, the curvature-vs-traced agreement for the traced hosts,
the positions and slope results (positions/results.json) and the knob comparison (knobs/class12_knob_ratio_traced.json).
Conclusion change = yes when a 95% interval moves across 1, a presence test crosses p = 0.05, or the class order changes.
Writes results.json and changes_traced.csv."""
import os, json, csv, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RR = os.path.dirname(HERE); O = f"{HERE}/terrain/out_v3"; G11 = f"{RR}/g11_curvature/terrain/out_v3"
J = lambda p: json.load(open(p))
rec = lambda p, k: J(p)[k]
F = rec(f"{G11}/peak_wall_record_rows_curv_n24_h100_s0.03.json", "rows_curv_n24_h100_s0.03"); C = rec(f"{O}/peak_wall_record_rows_traced_v3.json", "rows_traced_v3")
TF, TC = J(f"{G11}/terrain_table_record.json"), J(f"{O}/terrain_table_record.json"); EF, EC = J(f"{G11}/terrain_obs_exp.json"), J(f"{O}/terrain_obs_exp.json")
rf = pd.read_csv(f"{G11}/rows_curv_n24_h100_s0.03.csv"); rc = pd.read_csv(f"{O}/rows_traced_v3.csv"); rf = rf[(rf.D >= 30) & (rf.D <= 250)]; rc = rc[(rc.D >= 30) & (rc.D <= 250)]
fmt = lambda d: f"{d['ratio']:.2f} ({d['ci'][0]:.2f}-{d['ci'][1]:.2f})"
side = lambda lo, hi: "above 1" if lo > 1 else "below 1" if hi < 1 else "includes 1"
rows = []
def add(q, f, c, kind="ratio", note=""):
    if kind == "ratio":
        ch = side(*f["ci"]) != side(*c["ci"]); why = "; ".join(x for x in [f"interval {side(*f['ci'])} -> {side(*c['ci'])}" if ch else "", "" if f["ci"][0] <= c["ratio"] <= f["ci"][1] else "traced value outside the g11 interval"] if x)
        rows.append([q, fmt(f), fmt(c), "yes" if ch else "no", note or why])
    elif kind == "presence":
        fs, cs = f["presence_p"] < 0.05, c["presence_p"] < 0.05
        g = lambda d: f"{d['presence_obs']}/{d['presence_exp']} (z {d['presence_z']}, N {d['N']}); count obs/exp {d['count_obs_over_exp']} ({d['ci'][0]}-{d['ci'][1]})"
        rows.append([q, g(f), g(c), "yes" if fs != cs else "no", note])
    else: rows.append([q, f, c, kind, note])
A, B_ = F["all"], C["all"]
for q, k in [("peaks/walls, hollows", "pooled_hollows"), ("peaks/walls, hollow fields", "pooled_fields"), ("floors/walls, all floor", "floor_vs_wall"), ("floors/walls, inner floor (< 0.45 R)", "floor_in_vs_wall"),
             ("floors/walls, outer floor (>= 0.45 R)", "floor_out_vs_wall"), ("equator-/pole-facing walls", "eqwall_vs_polewall"), ("peaks/equator-facing walls", "peak_vs_eqwall"), ("peaks/walls, top 3 peak craters dropped", None)]:
    add(q, A[k] if k else F["drop_top3"], B_[k] if k else C["drop_top3"])
add("peaks/walls, slope-stratified (Mantel-Haenszel)", A["slope_stratified"]["MH_ratio"], B_["slope_stratified"]["MH_ratio"], "no", "")
GR = ["isolated", "nested", "central_peak", "ring", "unclassified", "filled_y", "filled_n", "fresh_s", "degraded_d", "D30_60", "D60_120", "D120_250", "isolated_central_peak"]
for g in GR:
    add(f"peaks/walls, {g} (N {C[g]['pooled_hollows']['N']})", F[g]["pooled_hollows"], C[g]["pooled_hollows"]); add(f"peaks/walls fields, {g}", F[g]["pooled_fields"], C[g]["pooled_fields"])
for q, k in [("crater-level presence, peak vs wall, hollows", "crater_level"), ("crater-level presence, peak vs wall, fields", "crater_level_fields"), ("crater-level presence, floor vs wall", "floor_vs_wall_crater"), ("crater-level presence, peak vs equator-facing wall", "peak_vs_eqwall_crater")]:
    add(q, A[k], B_[k], "presence")
for g in GR: add(f"crater-level presence, peak vs wall, {g}", F[g]["crater_level"], C[g]["crater_level"], "presence")
for s in ["all", "isolated", "basin", "nested"]:
    for t in ["peak", "floor", "wall"]:
        f, c = TF[s][t], TC[s][t]; rows.append([f"Table S2 density, {s} {t} (per 1,000 km2; n hollows)", f"{f['density']:.1f} ({f['lo']:.1f}-{f['hi']:.1f}); n {f['n']}", f"{c['density']:.1f} ({c['lo']:.1f}-{c['hi']:.1f}); n {c['n']}", "", ""])
    of = sorted(["peak", "floor", "wall"], key=lambda t: -TF[s][t]["density"]); oc = sorted(["peak", "floor", "wall"], key=lambda t: -TC[s][t]["density"])
    rows.append([f"Table S2 class order, {s}", " > ".join(of), " > ".join(oc), "yes" if of != oc else "no", ""])
for t in ["peak", "floor", "wall"]:
    f, c = EF["obs_over_exp"][t], EC["obs_over_exp"][t]; add(f"Table S2 observed/expected, {t}", dict(ratio=f["ratio"], ci=[f["lo"], f["hi"]]), dict(ratio=c["ratio"], ci=[c["lo"], c["hi"]]))
rows.append(["all-terrain density (per 1,000 km2)", f"{EF['all_terrain_density_per_1000km2']:.2f}", f"{EC['all_terrain_density_per_1000km2']:.2f}", "no", ""])
rows.append(["isolated peak/wall (S4)", f"{EF['isolated_peak_over_wall']:.2f}", f"{EC['isolated_peak_over_wall']:.2f}", "no", ""])
rows.append(["nested wall / isolated wall (S4)", f"{EF['nested_wall_over_isolated_wall']:.2f}", f"{EC['nested_wall_over_isolated_wall']:.2f}", "no", ""])
rows.append(["nested densest class", EF["nested_densest"], EC["nested_densest"], "no" if EF["nested_densest"] == EC["nested_densest"] else "yes", ""])
cnt = lambda r: "/".join(str(int(r[f"n_{k}"].sum())) for k in ["peak", "floor", "wall", "other", "masked"])
rows.append(["hollows per class, peak/floor/wall/other/masked", cnt(rf), cnt(rc), "", "3,761 hollows in 321 hosts"])
area = lambda r: "/".join(f"{r[f'A_{k}'].sum()/1000:.0f}" for k in ["peak", "floor", "wall"])
rows.append(["weighted area, peak/floor/wall (1,000 km2)", area(rf), area(rc), "", ""])
rows.append(["equator-/pole-facing wall hollows", f"{int(rf.n_wall_eq.sum())}/{int(rf.n_wall_pole.sum())}", f"{int(rc.n_wall_eq.sum())}/{int(rc.n_wall_pole.sum())}", "", ""])
# decomposition and sensitivity
DEC = {}
for cfg in ["curvonly_v2", "curvonly_v3", "traced_v3", "traced_v2", "traced_v3_certain"]:
    R = rec(f"{O}/peak_wall_record_rows_{cfg}.json", f"rows_{cfg}"); a = R["all"]
    DEC[cfg] = {k: dict(ratio=a[kk]["ratio"], ci=a[kk]["ci"], n_a=a[kk]["n_a"], n_b=a[kk]["n_b"]) for k, kk in [("peak_wall", "pooled_hollows"), ("fields_peak_wall", "pooled_fields"), ("floor_wall", "floor_vs_wall"), ("eq_pole_wall", "eqwall_vs_polewall")]}
    DEC[cfg]["N"] = a["pooled_hollows"]["N"]
lab = {"curvonly_v2": "curvature v2 everywhere (= g11)", "curvonly_v3": "curvature v3 everywhere", "traced_v3": "traced + curvature v3 (primary)", "traced_v2": "traced + curvature v2", "traced_v3_certain": "primary without the 13 flagged hosts"}
for cfg in ["curvonly_v3", "traced_v2", "traced_v3_certain"]:
    for k, l in [("peak_wall", "peaks/walls"), ("fields_peak_wall", "peaks/walls fields"), ("floor_wall", "floors/walls"), ("eq_pole_wall", "equator-/pole-facing walls")]:
        rows.append([f"sensitivity, {lab[cfg]}: {l}", fmt(DEC["traced_v3"][k]), fmt(DEC[cfg][k]), "yes" if side(*DEC["traced_v3"][k]["ci"]) != side(*DEC[cfg][k]["ci"]) else "no", "column 2 is the primary g13 value"])
# agreement of curvature classes with the tracings (traced hosts)
h = pd.read_csv(f"{O}/_chunks/hollow_classes_traced_v3_0000_0400.csv"); t = h[h.source == "traced"]
agree = {o: dict(agreement=round(float((t.cls == t[o]).mean()), 3), crosstab=pd.crosstab(t.cls, t[o]).to_dict()) for o in ("cls_v3", "cls_v2")}
ov = pd.read_csv(f"{O}/_chunks/area_overlap_traced_v3_0000_0400.csv"); AO = {}
for o in ("v3", "v2"):
    M = ov[ov.other == o].groupby(["traced", "cls"]).A.sum().unstack(); Tt, V, I = M.loc["peak"].sum(), M["peak"].sum(), M.loc["peak", "peak"]
    AO[o] = dict(traced_peak_km2=round(Tt), curv_peak_km2=round(V), overlap_km2=round(I), IoU=round(I / (Tt + V - I), 2), traced_inside_curv=round(I / Tt, 2), curv_inside_traced=round(I / V, 2), area_matrix_km2=(M.round(0)).to_dict())
P = J(f"{HERE}/positions/results.json"); K = J(f"{HERE}/knobs/class12_knob_ratio_traced.json")
tr3 = lambda v: dict(ratio=v[0], ci=[v[1], v[2]])
for e in ("eta", "eta_b"):
    a, b = P["positions"][f"g11_curv_{e}"], P["positions"][f"traced_{e}"]
    for k, l in [("peak_lowest_third_vs_upper_two_thirds", "lowest third of peak height"), ("peak_highest_third_vs_lower_two_thirds", "highest third of peak height"), ("peak_within_2km_of_floor_vs_rest_of_peak", "peak within 2 km of floor")]:
        for w in ("hollows", "fields"): add(f"positions ({e}), {l}, {w}", tr3(a[k][w]["ratio"]), tr3(b[k][w]["ratio"]))
Fa, Fb = P["floors"]["g11_curv"], P["floors"]["traced"]
for w in ("hollows", "fields"): add(f"floor within 0.2 R of center vs rest, {w}", tr3(Fa["floor_within_0.2R_of_center"][f"{w}_vs_rest_of_floor"]["ratio"]), tr3(Fb["floor_within_0.2R_of_center"][f"{w}_vs_rest_of_floor"]["ratio"]))
for k in ("floor_within_0.05R_of_peak_alone", "floor_within_0.1R_of_peak_alone", "floor_within_0.05R_of_wall_alone", "floor_within_0.1R_of_wall_alone"):
    for w in ("hollows", "fields"):
        add(f"{k.replace('_', ' ')}, {w}", tr3(Fa[k][w]["ratio"]), tr3(Fb[k][w]["ratio"]))
        add(f"{k.replace('_', ' ')}, {w}, radius-stratified", tr3(Fa[k][f"{w}_radius_stratified"]["ratio"]), tr3(Fb[k][f"{w}_radius_stratified"]["ratio"]))
for e in ("eta", "eta_b"):
    a, b = P["slope"][f"g11_curv_{e}"]["screens"]["main_30cells_norelief"], P["slope"][f"traced_{e}"]["screens"]["main_30cells_norelief"]
    for k in ("hollows_lowest_third_vs_upper_two_thirds", "hollows_highest_third_vs_lower_two_thirds", "fields_highest_third_vs_lower_two_thirds", "fields_lowest_third_vs_upper_two_thirds"):
        add(f"slope control ({e}, main screen), {k.replace('_', ' ')}, slope bins", tr3(a[k]["slope_bins"]), tr3(b[k]["slope_bins"]))
kc, kt = K["g11 curvature v2 (copied)"], K["traced where traced, curvature v3 elsewhere"]
add("knob density / Class 1-2 peak density", dict(ratio=kc["knob_over_peak"], ci=kc["ci"]), dict(ratio=kt["knob_over_peak"], ci=kt["ci"]), note=f"{K['n_traced']} of {K['n_craters']} Class 1-2 craters are traced")
res = dict(config=dict(catalog="v3", grid="g4_completeness/completeness_10ppd_v3.npz, Dref 1.0", hosts="30-250 km", N=int(len(rc)),
                       classes="Michael's traced outlines for 62 hosts (traces/traces_clean.json, status done); curvature v3 (g10 curv_classes_v3) for the other 259 hosts, including the 2 he marked can't tell"),
           terrain_test=dict(traced=C, g11_curvature=F), table_S2=dict(traced=TC, g11_curvature=TF), obs_exp=dict(traced=EC, g11_curvature=EF),
           decomposition=dict(labels=lab, values=DEC), uncertain_hosts=J(f"{O}/uncertain_traced_hosts.json"),
           curvature_vs_traced=dict(hollows_in_traced_hosts=int(len(t)), hollow_class_agreement=agree, peak_area=AO), positions=P, knob_comparison=K,
           changes=[dict(zip(["quantity", "g11_curvature_v2", "g13_traced", "conclusion_change", "note"], r)) for r in rows])
json.dump(res, open(f"{HERE}/results.json", "w"), indent=1, default=str)
with open(f"{HERE}/changes_traced.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["quantity", "g11 curvature v2", "g13 traced + curvature v3", "conclusion change", "note"]); w.writerows(rows)
for r in rows:
    if r[3] == "yes": print("CHANGE |", " | ".join(map(str, r)))
print(len(rows), "rows;", sum(r[3] == "yes" for r in rows), "conclusion changes")
