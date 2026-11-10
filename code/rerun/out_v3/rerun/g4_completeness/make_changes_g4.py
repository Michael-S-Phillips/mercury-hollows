#!/usr/bin/env python3
"""g4 (6 Nov 2026): writes changes_g4.csv. Columns as in the other changes tables, plus v3r / v3br (new completeness grid from the
refit detection model; v3 / v3b are the old-grid values of the earlier reruns) and conclusion_change_v2_to_v3r.
  conclusion_change          : does the conclusion change between old-grid v3 and new-grid v3r ('yes (v3r)'), or between v3b and
                               v3br ('yes (v3br)')? Same rules as the earlier tables (a 95% interval newly includes or excludes the
                               null value, p crosses 0.05, an AIC or dominance ranking flips). Hand-coded rows of the g1/g2 tables were
                               reviewed one by one (listed in OVERRIDE below).
  conclusion_change_v2_to_v3r: the earlier tables' v2 comparison redone with v3r / v3br in place of v3 / v3b.
Blocks: grid (compare_grids.py), g2 (g2/changes_g2_pass*.csv from g2/compile_changes.py), g1 (g1/changes_g1_pass*.csv from
g1/make_changes_g4.py), S2 reference size and E (s2/out), size distribution (g3 fixed and refit outputs; sfd/size_law.py).
Pass A = original mapping (reproduces the g1/g2 tables), B = v3r against v3, C = v3br against v3b, D = v3r / v3br against v2."""
import os, json, csv, re, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); RT = os.path.abspath(os.path.join(H, *[".."] * 5)); J = lambda p: json.load(open(p))
G3 = f"{RT}/code/rerun/out_v3/rerun/g3_catalog"; W = f"{RT}/manuscripts/PlanetaryResearch_hollows/analysis_workshop"
COLS = ["section", "quantity", "text_value", "v2_record", "v3", "v3b", "v3r", "v3br", "script", "conclusion_change", "conclusion_change_v2_to_v3r"]
rows = []
def add(section, quantity, text, v2, v3, v3b, v3r, v3br, script, cc, cc2):
    rows.append(dict(zip(COLS, [section, quantity, text, v2, v3, v3b, v3r, v3br, script, cc, cc2])))
ren = lambda s: re.sub(r"\bv3\b", "v3r", re.sub(r"\bv3b\b", "v3br", s))
# ------------------------------------------------------------------ grid
GC = J(f"{H}/grid_comparison.json"); PI = J(f"{H}/g2/inputs/prep_g4_record.json")
f3 = lambda q: f"{q[0]:.3f} / {q[1]:.3f} / {q[2]:.3f}"
for D in ("0.5", "1.0", "2.0", "5.0"):
    a, b = GC["v3"][D], GC["v3b"][D]
    add("S2 (grid)", f"completeness new/old over covered cells, {float(D):g} km: p5 / p50 / p95 (area-weighted in brackets)", "-", "1 (v2 grid)", "1 (v2 grid)", "1 (v2 grid)",
        f"{f3(a['p5_p50_p95'])} [{f3(a['area_weighted_p5_p50_p95'])}]", f"{f3(b['p5_p50_p95'])} [{f3(b['area_weighted_p5_p50_p95'])}]", "compare_grids.py", "n/a", "n/a")
add("S2 (grid)", "1 km: p5 and p95 relative to the median (rerun-depth criterion: within +-3%)", "-", "-", "-", "-",
    f"{100*GC['v3']['1.0']['rel_p5_p95'][0]:+.1f}% / {100*GC['v3']['1.0']['rel_p5_p95'][1]:+.1f}% (not uniform: full reruns)",
    f"{100*GC['v3b']['1.0']['rel_p5_p95'][0]:+.1f}% / {100*GC['v3b']['1.0']['rel_p5_p95'][1]:+.1f}%", "compare_grids.py", "n/a", "n/a")
for k, lab in (("ratio_p5_p50_p95", "all 31,554 craters with C1 > 0"), ("h09_hosts_ratio_p5_p50_p95", "h09 hosts")):
    add("S2 (grid)", f"C1 at crater centers, new/old p5 / p50 / p95: {lab}", "-", "-", "1", "1", f3(PI["v3r"]["crater_features_C1"][k]), f3(PI["v3br"]["crater_features_C1"][k]), "g2/prep_inputs.py", "n/a", "n/a")
add("S2 (grid)", "per-hollow C_own and C_D1km, median new/old", "-", "-", "1", "1",
    f"{PI['v3r']['hollows_canonical']['C_own']['median_ratio']:.3f}; {PI['v3r']['hollows_canonical']['C_D1km']['median_ratio']:.3f}",
    f"{PI['v3br']['hollows_canonical']['C_own']['median_ratio']:.3f}; {PI['v3br']['hollows_canonical']['C_D1km']['median_ratio']:.3f}", "build_completeness.py attach; g2/prep_inputs.py", "n/a", "n/a")
# ------------------------------------------------------------------ g2 and g1 tables
OVERRIDE = {  # (group, row index): (grid-effect flag, reason) for rows whose flag in pass B/C is hand-coded in the source table
 ("g2", 43): ("no", "p stays > 0.05 in both bands"), ("g2", 57): ("no", "still beyond the largest crater; interaction n.s."),
 ("g2", 92): ("no", "both p > 0.05"), ("g2", 93): ("no", "same pattern"),
 ("g2", 26): ("yes (v3br)", "bootstrap interval now excludes 1 (0.43-0.99) as in v3; p = 0.056 (v3b 0.075)"),
 ("g1", 102): ("no", "Class 4 -> 5 falls in both"),
 ("g1", 152): ("no", "but the v3r value is not usable: the cloglog fit fails to converge (negative cover share), as g1 found for v3b"),
}
for g, scr in (("g2", "g2/"), ("g1", "g1/")):
    A, B, C, D = (pd.read_csv(f"{H}/{g}/changes_{g}_pass{x}.csv", keep_default_na=False) for x in "ABCD")
    REF = pd.read_csv(f"{H}/../{'g2_insolation' if g == 'g2' else 'g1_graph_joint'}/changes_{g}.csv", keep_default_na=False)
    assert A.equals(REF), f"pass A does not reproduce changes_{g}.csv"
    for i in range(len(A)):
        v3r, v3br = B.v3[i], C.v3[i]; na_br = v3br.startswith("n/a")
        if na_br: v3br = "not run (summary model run for v3r only)"
        fl = []
        if v3r != A.v3[i] and B.conclusion_change[i] != "no": fl.append("v3r")
        if not na_br and C.v3[i] != A.v3b[i] and C.conclusion_change[i] != "no": fl.append("v3br")
        cc = "no" if not fl else f"yes ({', '.join(fl)})"
        if (g, i) in OVERRIDE: cc = OVERRIDE[(g, i)][0] + "; " + OVERRIDE[(g, i)][1]
        elif g == "g1" and i in (12, 107, 127) and "v3br" in fl: cc = "yes (v3br): p 0.043 -> 0.054 (degradation-state-controlled inheritance)"
        elif fl: cc += " - CHECK"
        d = D.conclusion_change[i]
        if na_br: d = re.sub(r",? ?\bv3b\b", "", d).replace("(, ", "(").replace("()", "").strip()
        d = ren(d) if not na_br else re.sub(r"\bv3\b", "v3r", d)
        if d in ("yes", "yes ("): d = "no"
        add(A.section[i], A.quantity[i], A.text_value[i], A.v2_record[i], A.v3[i], A.v3b[i], v3r, v3br, scr + A.script[i], cc, d)
# ------------------------------------------------------------------ S2: reference size and E
S = {c: dict(ins=J(f"{H}/s2/out/ref_size_insolation_{c}.json"), chk=J(f"{H}/s2/out/ref_size_record_{c}.json"), hf=J(f"{H}/s2/out/hosting_fits_{c}.json"),
             jt=J(f"{H}/s2/out/joint_thinned_{c}.json")) for c in ("v2", "v3", "v3b", "v3r", "v3br")}
S["v2"]["jt"] = {**S["v2"]["jt"], **J(f"{W}/completeness_sfd/joint_thinned.json")}   # stored record (identical where both exist) + lnE_1.7/1.8 from the copy
CS = ("v2", "v3", "v3b", "v3r", "v3br"); SCR = "s2/"
def rng(c, band):
    v = [S[c]["ins"][D][band]["hot_OR"] for D in ("0.3", "0.5", "1.0")]
    sig = "".join("*" if x[3] < 0.05 else "" for x in v[:0])
    return f"{min(x[0] for x in v):.2f}-{max(x[0] for x in v):.2f} (0.3/0.5/1 km p = " + ", ".join(f"{x[3]:.2g}" for x in v) + ")"
for band, lab, txt in (("-15..15", "within 15 deg", "2.2-2.6"), ("15..30", "15-30N", "3.2-3.7"), ("30..45", "30-45N", "1.76-1.98")):
    vals = [rng(c, band) for c in CS]
    p = {c: [S[c]["ins"][D][band]["hot_OR"][3] for D in ("0.3", "0.5", "1.0")] for c in CS}
    st = lambda c: tuple(x < 0.05 for x in p[c])
    cc = "no" if st("v3r") == st("v3") and st("v3br") == st("v3b") else "yes"
    cc2 = "no" if all(st(c) == st("v2") for c in ("v3r", "v3br")) else ("yes (v3r, v3br): 1 km p > 0.05 in both (old grid alike)" if band == "30..45" else "yes")
    add("S2 (reference size)", f"hot vs warm odds ratio of hosting with completeness for 0.3, 0.5 or 1 km hollows, {lab} (range)", txt, *vals, SCR + "ref_size_insolation.py", cc, cc2)
vals = [f"{min(S[c]['ins'][D]['dT_per10K_lt60'][0] for D in ('0.3', '0.5', '1.0')):.2f}-{max(S[c]['ins'][D]['dT_per10K_lt60'][0] for D in ('0.3', '0.5', '1.0')):.2f}" for c in CS]
add("S2 (reference size)", "odds per 10 K of Delta T_w within 60 deg, completeness for 0.3, 0.5 or 1 km (range)", "1.43-1.47", *vals, SCR + "ref_size_insolation.py", "no", "no")
fo = lambda x: f"{x[0]:.2f} ({x[1]:.2f}-{x[2]:.2f}; p = {x[3]:.2g})"
vals = [f"2 km {fo(S[c]['ins']['2.0']['30..45']['hot_OR'])}; 5 km {fo(S[c]['ins']['5.0']['30..45']['hot_OR'])}" for c in CS]
add("S2 (reference size)", "30-45N hot vs warm odds ratio with completeness for 2 and 5 km hollows", "1.4 (0.9-2.2); 1.0 (0.6-1.6)", *vals, SCR + "ref_size_insolation.py", "no", "no")
def sigsizes(c):
    return ", ".join(D.rstrip("0").rstrip(".") for D in ("0.3", "0.5", "1.0", "2.0", "5.0") if S[c]["ins"][D]["30..45"]["hot_OR"][3] < 0.05) or "none"
vals = [f"significant (p < 0.05) with: {sigsizes(c)} km" for c in CS]
add("S2 (reference size)", "'The hot-longitude excess at 30-45N loses significance only with completeness for 2 or 5 km hollows'", "only with 2 or 5 km", *vals,
    SCR + "ref_size_insolation.py", "no" if sigsizes("v3r") == sigsizes("v3") and sigsizes("v3br") == sigsizes("v3b") else "yes",
    "yes (v3r, v3br): already lost at 1 km (also in old-grid v3, v3b); kept only with 0.3 and 0.5 km" if sigsizes("v3r") == "0.3, 0.5" else "check")
vals = [f"<15 {fo(S[c]['ins']['2.0']['-15..15']['hot_OR'])}; 15-30 {fo(S[c]['ins']['2.0']['15..30']['hot_OR'])}" for c in CS]
add("S2 (reference size)", "hot vs warm odds ratio with completeness for 2 km hollows, within 15 deg and 15-30N (not quoted)", "not quoted", *vals, SCR + "ref_size_insolation.py",
    "no", "yes (v3r, v3br): both bands lose significance with 2 km completeness (old grid alike)")
vals = [f"1 km {S[c]['chk']['by_ref']['1.0']['llf']:,.1f}; 2 km {S[c]['chk']['by_ref']['2.0']['llf']:,.1f}" for c in CS]
add("S2 (reference size)", "log-likelihood of the containment (nesting) fit with completeness for 1 and 2 km hollows", "-2,100 and -2,093 (about equal)", *vals, SCR + "ref_size_check.py",
    "no", "yes (v3r, v3br): 1 km now fits better than 2 km, by " + " / ".join(f"{S[c]['chk']['by_ref']['1.0']['llf'] - S[c]['chk']['by_ref']['2.0']['llf']:.1f}" for c in ("v3r", "v3br")) + " (old grid v3 / v3b: " + " / ".join(f"{S[c]['chk']['by_ref']['1.0']['llf'] - S[c]['chk']['by_ref']['2.0']['llf']:.1f}" for c in ("v3", "v3b")) + ")")
for key, lab, txt in (("nested_OR", "craters inside older craters", "1.85-1.91"), ("vent_OR", "craters with volcanic vents", "5.0-6.7"), ("fresh_vs_degraded_OR", "fresh vs degraded craters >= 10 km", "3.3-4.5")):
    vals = [f"{min(S[c]['chk']['by_ref'][D][key][0] for D in ('0.3', '0.5', '1.0')):.2f}-{max(S[c]['chk']['by_ref'][D][key][0] for D in ('0.3', '0.5', '1.0')):.2f}" for c in CS]
    add("S2 (reference size)", f"odds ratio of hosting with completeness for 0.3, 0.5 or 1 km: {lab} (range; every interval excludes 1)", txt, *vals, SCR + "ref_size_check.py", "no", "no (values lower, as in old-grid v3)")
jt = lambda c, k: S[c]["jt"][k]["llf_full"]
vals = [f"1 km {jt(c, 'lnCold_1.0'):,.1f}; others {min(jt(c, 'lnCold_' + D) for D in ('0.3', '0.5', '0.7', '1.5', '2.0', '3.0', '5.0')):,.1f} to {max(jt(c, 'lnCold_' + D) for D in ('0.3', '0.5', '0.7', '1.5', '2.0', '3.0', '5.0')):,.1f} (best: {max(('0.3', '0.5', '0.7', '1.0', '1.5', '2.0', '3.0', '5.0'), key=lambda D: jt(c, 'lnCold_' + D))} km)" for c in CS]
add("S2 (reference size)", "joint model log-likelihood with completeness for 1 km against 0.3-5 km hollows", "-1,969 against -1,976 to -2,229", *vals, SCR + "joint_thinned.py",
    "no", "yes (v3r, v3br): 0.7 km now fits slightly better than 1 km (by 5.2 in v3r; old-grid v3: 4.0)")
for M, lab in (("lnE_1.6", "E with exponent 1.6 (as in S2)"), ("lnE_1.7", "E with exponent 1.7 (v3 value)")):
    h = lambda c, k: S[c]["hf"][M][k]; f = lambda x, d=2: f"{x[0]:.{d}f} ({x[1]:.{d}f}-{x[2]:.{d}f})"
    vals = [f"nested {f(h(c, 'nested_OR'))}; vents {f(h(c, 'vent_OR'), 1)}; fresh {f(h(c, 'fresh_vs_degraded_OR'), 1)}" for c in CS]
    add("S2 (E)", f"{lab}: odds ratios for nesting, vents, fresh vs degraded >= 10 km", "1.91 (1.59-2.29); 5.6 (3.4-9.2); 3.6 (2.4-5.4)" if M == "lnE_1.6" else "-", *vals, SCR + "hosting_fits.py", "no", "no (lower, all exclude 1)")
    vals = [f"{f(h(c, 'hot_OR_-15..15'))}; {f(h(c, 'hot_OR_15..30'))}; {f(h(c, 'hot_OR_30..45'))} p = {h(c, 'hot_OR_30..45')[3]:.2g}" for c in CS]
    add("S2 (E)", f"{lab}: hot vs warm odds ratio within 15 deg; 15-30N; 30-45N", "2.0 (1.5-2.7); 2.8 (1.9-4.2); 1.70 (1.08-2.68)" if M == "lnE_1.6" else "-", *vals, SCR + "hosting_fits.py",
        "no", "yes (v3r, v3br): 30-45N interval includes 1 (old grid alike)")
    vals = [f(h(c, "dT_OR_per10K_lt60")) for c in CS]
    add("S2 (E)", f"{lab}: odds per 10 K of Delta T_w within 60 deg", "1.43 (1.32-1.54)" if M == "lnE_1.6" else "-", *vals, SCR + "hosting_fits.py", "no", "no")
def maxdev(c, base="lnE_1.6"):
    keys = ["nested_OR", "vent_OR", "fresh_vs_degraded_OR", "hot_OR_-15..15", "hot_OR_15..30", "hot_OR_30..45", "dT_OR_per10K_lt60"]
    return max(abs(S[c]["hf"][m][k][0] / S[c]["hf"][base][k][0] - 1) for m in ("lnE_1.3", "lnE_1.9") for k in keys)
vals = [f"{100 * maxdev(c):.1f}% (from 1.7: {100 * maxdev(c, 'lnE_1.7'):.1f}%)" for c in CS]
add("S2 (E)", "largest change of these odds ratios for exponents 1.3 and 1.9 (from 1.6)", "< 5%", *vals, SCR + "hosting_fits.py", "no",
    "no" if maxdev("v3r") < .05 else "yes (v3r): just above 5% from the 1.7 base")
sh = lambda c, M: S[c]["jt"][M]["share"]
vals = [f"{100 * max(abs(sh(c, 'lnE_1.6')[g] - sh(c, 'lnC1')[g]) for g in sh(c, 'lnC1')):.1f} pp (exponent 1.7: {100 * max(abs(sh(c, 'lnE_1.7')[g] - sh(c, 'lnC1')[g]) for g in sh(c, 'lnC1')):.1f} pp); log-likelihood {jt(c, 'lnE_1.6'):,.0f} against {jt(c, 'lnC1'):,.0f}" for c in CS]
add("S2 (E)", "joint model with E: largest change of a group share against C1 (pp); log-likelihood E vs C1", "at most 1.2 pp; -2,022 (poorer than -1,969)", *vals, SCR + "joint_thinned.py",
    "no", "number changes (about 1.5-1.8 pp); E still fits worse")
# ------------------------------------------------------------------ size distribution: fixed (v3, v3b) and refit (v3r, v3br)
SN = {"rec": J(f"{W}/completeness_sfd/sfd_numbers.json"), "v3": J(f"{G3}/sfd/out_v3_fixed/sfd_numbers.json"), "v3b": J(f"{G3}/sfd/out_v3b_fixed/sfd_numbers.json"),
      "v3r": J(f"{G3}/sfd/out_v3_refit/sfd_numbers.json"), "v3br": J(f"{G3}/sfd/out_v3b_refit/sfd_numbers.json")}
SF = {"rec": J(f"{W}/completeness_sfd/sfd_fit.json"), "v3": J(f"{G3}/sfd/out_v3_fixed/sfd_fit.json"), "v3b": J(f"{G3}/sfd/out_v3b_fixed/sfd_fit.json"),
      "v3r": J(f"{G3}/sfd/out_v3_refit/sfd_fit.json"), "v3br": J(f"{G3}/sfd/out_v3b_refit/sfd_fit.json")}
WN = {"rec": J(f"{W}/completeness_sfd/sfd_windows.json"), "v3": J(f"{G3}/sfd/out_v3_fixed/sfd_windows.json"), "v3b": J(f"{G3}/sfd/out_v3b_fixed/sfd_windows.json"),
      "v3r": J(f"{G3}/sfd/out_v3_refit/sfd_windows.json"), "v3br": J(f"{G3}/sfd/out_v3b_refit/sfd_windows.json")}
DFx, DRx = J(f"{G3}/sfd/sfd_derived_fixed.json"), J(f"{G3}/sfd/sfd_derived_refit.json")
DD = {"rec": DFx["v2"], "v3": DFx["v3"], "v3b": DFx["v3b"], "v3r": DRx["v3"], "v3br": DRx["v3b"]}
SL = {c: J(f"{H}/sfd/size_law_record_{c}.json")["size_" + c] for c in ("v3", "v3b", "v3r", "v3br")}; SL["rec"] = J(f"{RT}/code/rerun/closure/size_law_record.json")["size_v2"]
V5 = ("rec", "v3", "v3b", "v3r", "v3br")
ci = lambda t: f"{t[0]:.2f} ({t[1]:.2f}-{t[2]:.2f})"; pct = lambda x: f"{100 * x:.1f}%"
def sadd(sec, q, text, fn, script, cc="no", cc2="no"):
    v = [fn(c) for c in V5]; add(sec, q, text, *v, script, cc, cc2)
th = lambda t: f"{t['p_max']:.3f}, {t['x50_px']:.1f}, {t['sigma_lo']:.3f}, {t['x_hi_px']:.1f}, {t['sigma_hi']:.3f}, {t['b_inc']:.3f}"
sadd("S2; 2.1", "detection model p_max, x50, sigma_lo, x_hi, sigma_hi, b (v3r/v3br: refit on the v3/v3b pairs, now primary)", "0.243, 35.5, 0.590, 195.9, 0.267, -0.102",
     lambda c: th(SF[c]["theta"]), "g3 sfd/refit_detection.py", "n/a (model change)", "values change; form unchanged")
for k, lab, txt in (("slope_0.3_0.7", "0.3-0.7 km", "1.57 (1.38-1.76)"), ("slope_0.7_3", "0.7-3 km", "1.38 (1.30-1.46)"), ("slope_1_3", "1-3 km", "1.43 (1.31-1.55)"), ("slope_3_12", "above 3 km (3-12 km)", "2.68 (2.43-2.94)")):
    sadd("3.1; S2", f"exponent of N(>D), {lab}", txt, lambda c, k=k: ci(SN[c][k]), "g3 sfd/sfd_numbers.py")
sadd("3.1; Fig 1", "exponent range 0.3-3 km (D^-a), rounded [exact]", "1.4 to 1.6", lambda c: f"{DD[c]['rounded']['range'][0]}-{DD[c]['rounded']['range'][1]} [{DD[c]['exact']['range'][0]}-{DD[c]['exact']['range'][1]}]", "g3 sfd/sfd_derived.py")
sadd("3.1 (support)", "exponent 2-6 km (where the steepening starts)", "not quoted", lambda c: ci(WN[c]["windows"]["2.0-6.0"]), "g3 sfd/sfd_windows.py")
cnt = lambda c, k: f"{SN[c]['counts'][k]['corrected']:,} ({SN[c]['counts'][k]['corrected_ci'][0]:,}-{SN[c]['counts'][k]['corrected_ci'][1]:,})"
sadd("3.1", "corrected N(>5 km)", "214 (182-246)", lambda c: cnt(c, "5.0"), "g3 sfd/sfd_numbers.py")
sadd("3.1", "power-law N(>5 km) if the 0.3-3 km law continued, rounded [exact]", "430-590", lambda c: f"{DD[c]['rounded']['N_gt_5km_powerlaw'][0]:.0f}-{DD[c]['rounded']['N_gt_5km_powerlaw'][1]:.0f} [{DD[c]['exact']['N_gt_5km_powerlaw'][0]:.0f}-{DD[c]['exact']['N_gt_5km_powerlaw'][1]:.0f}]", "g3 sfd/sfd_derived.py")
sadd("3.1", "median detection probability at 3 and 5 km", "82-86%", lambda c: f"{pct(SN[c]['median_C_at_size']['3'])}, {pct(SN[c]['median_C_at_size']['5'])}", "g3 sfd/sfd_numbers.py")
sadd("2.1; 3.1", "median detection probability of a 1 km / 0.3 km hollow at cataloged sites", "47% / 4%", lambda c: f"{pct(SN[c]['median_C_at_size']['1'])} / {pct(SN[c]['median_C_at_size']['0.3'])}", "g3 sfd/sfd_numbers.py")
sadd("3.1", "share of covering images in which 0.3 km spans < 10 px; median pixel scale", "80%; 63 m", lambda c: f"{pct(DD[c]['pairs_0p3km_lt_10px'])}; {DD[c]['pairs_median_px_m']} m", "g3 sfd/sfd_derived.py (no detection model)")
sadd("3.1; S2", "other covering images that detected a cataloged hollow", "5% (S2: 4.7%)", lambda c: pct(DD[c]["pairs_other_images_detected"]), "g3 sfd/sfd_derived.py (no detection model)")
sadd("3.1; S2", "same, where the hollow spans 60-100 / 100-200 px", "14-15%", lambda c: " / ".join(pct(x) for x in DD[c]["pairs_detected_Dpx_60_100_100_200"]), "g3 sfd/sfd_derived.py (no detection model)")
sadd("3.1", "largest detection vs the one used, median excess, 3-5 and 5-10 km", "5-6%", lambda c: f"{pct(DD[c]['largest_detection_excess_median']['3-5'][0])}, {pct(DD[c]['largest_detection_excess_median']['5-10'][0])}", "g3 sfd/sfd_derived.py (no detection model)")
sadd("3.1", "cataloged hollows > 0.3 km", "6,044", lambda c: f"{SN[c]['counts']['0.3']['catalog']:,}", "g3 sfd/sfd_numbers.py")
sadd("3.1", "corrected N(>0.3 km)", "24,100 (23,200-25,000)", lambda c: cnt(c, "0.3"), "g3 sfd/sfd_numbers.py")
sadd("3.1", "power-law N(>0.3 km) from 1 km, rounded [exact]", "30,000-38,000", lambda c: f"{DD[c]['rounded']['N_gt_0p3km_powerlaw'][0]:,.0f}-{DD[c]['rounded']['N_gt_0p3km_powerlaw'][1]:,.0f} [{DD[c]['exact']['N_gt_0p3km_powerlaw'][0]:,.0f}-{DD[c]['exact']['N_gt_0p3km_powerlaw'][1]:,.0f}]", "g3 sfd/sfd_derived.py")
sadd("3.1", "catalog share of hollows > 0.3 km", "one in five to one in six", lambda c: f"1 in {1/DD[c]['rounded']['catalog_fraction_0p3'][0]:.1f} to 1 in {1/DD[c]['rounded']['catalog_fraction_0p3'][1]:.1f}", "g3 sfd/sfd_derived.py")
sadd("3.1 (support)", "corrected N(>1 km) (anchor of the power-law lines)", "not quoted", lambda c: f"{SN[c]['counts']['1.0']['corrected']:,}", "g3 sfd/sfd_numbers.py")
sadd("3.1; S2", "best-imaged sites (C(0.3 km) >= 0.3): hollows, hollows > 0.3 km, exponent 0.3-1 km", "193; 106; 2.4 (1.8-3.1)", lambda c: f"{SN[c]['sites_C03_ge_0.3']['n']}; {SN[c]['sites_C03_ge_0.3']['n_ge03']}; {ci(SN[c]['sites_C03_ge_0.3']['slope_03_1'])}", "g3 sfd/sfd_numbers.py",
     "no (more sites qualify under the refit; exponent unchanged)", "no")
sadd("S2", "exponent above lower limits 0.2, 0.3, 0.5 km", "1.55, 1.63, 1.62", lambda c: ", ".join(f"{SF[c]['power'][k]['slope_cumulative']:.2f}" for k in ("0.2", "0.3", "0.5")), "g3 sfd/fit_sfd.py")
sadd("S2", "exponent for the size-weighted completeness E (D^-a)", "1.6 (1.55-1.63)", lambda c: {"rec": "1.6", "v3": "1.7 (1.61-1.72)", "v3b": "1.7-1.8 (1.63-1.76)",
     "v3r": "1.6-1.7 (1.57-1.68)", "v3br": "1.6-1.7 (1.57-1.69)"}[c], "g3 sfd/fit_sfd.py", "no", "no (E fits run with 1.6, 1.7 and 1.8; S2 (E) rows)")
sadd("S2", "lognormal fit median (km), sigma", "not quoted", lambda c: f"{SF[c]['lognormal']['median_km']:.3g}, {SF[c]['lognormal']['sigma_lnD']:.2f}", "g3 sfd/fit_sfd.py")
sadd("4.4", "HRIC: hollows > 0.1 km per hollow > 0.3 km, rounded [exact]", "4.7-5.8 times", lambda c: f"{DD[c]['rounded']['ratio_0p1_over_0p3'][0]:.1f}-{DD[c]['rounded']['ratio_0p1_over_0p3'][1]:.1f} [{DD[c]['exact']['ratio_0p1_over_0p3'][0]:.1f}-{DD[c]['exact']['ratio_0p1_over_0p3'][1]:.1f}]", "g3 sfd/sfd_derived.py")
sadd("4.4", "hollows > 0.1 km on the searched surface, rounded [exact]", "140,000-220,000", lambda c: f"{DD[c]['rounded']['N_gt_0p1km'][0]:,.0f}-{DD[c]['rounded']['N_gt_0p1km'][1]:,.0f} [{DD[c]['exact']['N_gt_0p1km'][0]:,.0f}-{DD[c]['exact']['N_gt_0p1km'][1]:,.0f}]", "g3 sfd/sfd_derived.py")
sadd("4.4", "hollows > 5 km: corrected total against the catalog", "about 210 against 185", lambda c: f"{DD[c]['N_gt_5km']['corrected']} against {DD[c]['N_gt_5km']['catalog']}", "g3 sfd/sfd_derived.py")
sadd("closure (not in v6.3)", "size_law.py differential exponent, north, 0.1-10 km (uses the 10 ppd grid and C_own)", "not in v6.3",
     lambda c: f"{SL[c]['alpha_point']:.3f} ({SL[c]['alpha_95CI'][0]:.3f}-{SL[c]['alpha_95CI'][1]:.3f}); n = {SL[c]['n_fitted']}", "sfd/size_law.py (g4 copy)")
# ------------------------------------------------------------------ write
with open(f"{H}/changes_g4.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(rows)
df = pd.DataFrame(rows); print(len(df), "rows; grid-effect conclusion changes:", int(df.conclusion_change.str.startswith("yes").sum()))
print(df[df.conclusion_change.str.startswith("yes") | df.conclusion_change.str.contains("CHECK")][["section", "quantity", "v3", "v3r", "v3b", "v3br", "conclusion_change"]].to_string()[:4000])
