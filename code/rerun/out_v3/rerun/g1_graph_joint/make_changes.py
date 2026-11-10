#!/usr/bin/env python3
"""g1 rerun (6 Nov 2026): builds changes_g1.csv, one row per catalog-dependent number in main-text Sections 2.2 (definitions
used by 3.3-3.8), 3.3, 3.4, 3.6 (hosting only), 3.8, Figure 3, Table 1, the Discussion and Conclusions, and supplement
S11 (Tables S13-S14), S14 (Tables S18-S19), S15 (Table S20), S16 (Table S21, Figure S9).
v2_record: from the existing record files (analysis_v9, analysis_v8/alt_caloris_orgel, analysis_v3, analysis_workshop);
'no record' where no script or record was kept (then the v2 value from this folder's v2 rerun is given).
v3, v3b: from this folder. conclusion_change: 'yes' when, in v3 or v3b, an estimate changes sign, a 95% interval newly
includes (or excludes) 1 or 0, p crosses 0.05, an AIC ranking changes, or a dominance ranking changes; the catalog(s)
concerned are named in parentheses. A sign flip of an estimate whose interval includes the null value in both catalogs is
not counted; ranks are compared group by group, so a swap of two adjacent groups marks both rows."""
import os, json, csv, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); P = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows")
W = f"{P}/analysis_workshop"
def J(f): return json.load(open(f))


def load(C):
    if C == "v2":
        A, SM, DG = f"{P}/analysis_v9", f"{W}/summary_model/out", f"{W}/kinczyk_degradation"
        d = dict(dag=J(f"{A}/crater_dag_record.json"), tests=J(f"{A}/crater_dag_tests_record.json"), sens=J(f"{A}/crater_dag_sens_record.json"),
                 lunar=J(f"{A}/uplift_lunar/crater_dag_record.json"), k08=J(f"{A}/kring_0p8/crater_dag_record.json"), k13=J(f"{A}/kring_1p3/crater_dag_record.json"),
                 orgel=J(f"{P}/analysis_v8/alt_caloris_orgel/crater_dag_sens_record.json"), hage=J(f"{P}/analysis_v3/host_age_record.json"),
                 kin=J(f"{DG}/kinczyk_tests.json"), stage=J(f"{DG}/stage_by_class.json"), msens=J(f"{DG}/match_sensitivity.json"),
                 t8=J(f"{W}/test8_vents_record.json"), vac=J(f"{W}/vent_check/vent_assignment_check.json"),
                 pe=J(f"{H}/v2/parents_and_examples.json"))
    else:
        SM = f"{H}/summary_model/{C}/out"
        d = dict(dag=J(f"{H}/{C}/crater_dag_record.json"), tests=J(f"{H}/{C}/crater_dag_tests_record.json"), sens=J(f"{H}/{C}/crater_dag_sens_record.json"),
                 lunar=J(f"{H}/{C}/uplift_lunar/crater_dag_record.json"), k08=J(f"{H}/{C}/kring_0p8/crater_dag_record.json"), k13=J(f"{H}/{C}/kring_1p3/crater_dag_record.json"),
                 orgel=J(f"{H}/{C}/caloris_orgel/crater_dag_sens_record.json"), hage=J(f"{H}/degradation/{C}/host_age_record.json"),
                 kin=J(f"{H}/degradation/{C}/kinczyk_tests.json"), stage=J(f"{H}/degradation/{C}/stage_by_class.json"), msens=J(f"{H}/degradation/{C}/match_sensitivity.json"),
                 t8=J(f"{H}/vents/{C}/test8_vents_record.json"), vac=J(f"{H}/vents/{C}/vent_assignment_check.json"), pe=J(f"{H}/{C}/parents_and_examples.json"))
    d.update(imp=J(f"{SM}/improved.json"), dom=J(f"{SM}/dominance.json"), jm=pd.read_csv(f"{SM}/joint_model.csv").set_index("term"),
             col=J(f"{SM}/color_record.json"), trees=J(f"{SM}/trees.json"), size=J(f"{SM}/sizestep_check.json"),
             xrr=pd.read_csv(f"{SM}/xrs_rate_ratios.csv").set_index("term"))
    # count model: v2 from the record; v3 and v3b from count_model_check.py (the statsmodels alpha fit in step7 fails for v3b)
    if C == "v2": d.update(cnt=d["imp"]["thinned_counts_all"]["point"], cnt_alpha=d["imp"]["thinned_counts_all"]["alpha"])
    else: cm = J(f"{SM}/count_model_check.json"); d.update(cnt=cm["point"], cnt_alpha=cm["alpha"])
    return d


D = {C: load(C) for C in ["v2", "v3", "v3b"]}


# ---------- value objects: dict(s=string, est, lo, hi, p, null, rankkey, change, norecord) ----------
def fp(p): return f"{p:.2g}"
def pct(x, nd=0): return f"{100 * x:.{nd}f}%"
def V(s, **k): return dict(s=str(s), **k)
def OR(est, lo, hi, p=None, nd=2):
    s = f"{est:.{nd}f} ({lo:.{nd}f}-{hi:.{nd}f}" + (f"; p = {fp(p)})" if p is not None else ")")
    return dict(s=s, est=float(est), lo=float(lo), hi=float(hi), p=p, null=1.0)
def tparam(d, model, term):                                     # crater_dag_tests params [b, se, p]
    b, se, p = d["tests"][model]["params"][term]; return OR(np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p)
def sensOR(d, key, term="nested"):
    v = d["sens"][key][term]; return OR(v["OR"], v["ci"][0], v["ci"][1], v["p"])
def rank(sh, g): return sorted(sh, key=lambda k: -sh[k]).index(g) + 1
def order(sh, drop=()): return " > ".join(k for k in sorted(sh, key=lambda k: -sh[k]) if k not in drop)
def share(d, an, g):
    T = d["imp"][an]; v = T["point"]["share"][g]
    return f"{100*v:.1f}" + (f" ({100*T['ci'][g][0]:.1f}-{100*T['ci'][g][1]:.1f})" if "ci" in T else "")
def aicrow(d, rec, model, term=None, fmt="{:.4f}"):
    R = d[rec]; a = R[model]["aic"]
    if term is None: return V(f"{a:.1f}", aic=a)
    b, se, p = R[model]["params"][term]
    return V(f"{a:.1f}; {fmt.format(b)} +- {fmt.format(se)} (p = {fp(p)})", aic=a, est=b, lo=b - 1.96 * se, hi=b + 1.96 * se, p=p, null=0.0)
def jm(d, t):
    r = d["jm"].loc[t]; return OR(r.OR, r.OR_lo_block, r.OR_hi_block)
def kc(d, key, sub):
    v = d["kin"][key][sub]; return OR(v[0], v[1], v[2], v[3] if len(v) > 3 else None)
def best_depth_minus_nested(R): return min(R["z_base"]["aic"], R["ln_zown"]["aic"], R["ln_zmax"]["aic"], R["ln_zown_zmax"]["aic"]) - R["nested"]["aic"]
def rates(r): return [tuple(int(x) for x in r[k].split("/")) for k in ["1", "0"]]


rows = []
def row(section, quantity, text, getter, script, kind="value", norecord=False):
    vals = {}
    for C in ["v2", "v3", "v3b"]:
        try: vals[C] = getter(D[C], C)
        except Exception as e: vals[C] = V(f"n/a ({type(e).__name__}: {e})")
    v2s = vals["v2"]["s"] if not norecord else f"no record (v2 rerun: {vals['v2']['s']})"
    cc = []
    for C in ["v3", "v3b"]:
        a, b = vals["v2"], vals[C]
        if kind in ("ratio", "coef") and "lo" in a and "lo" in b:
            n0 = a.get("null", 1.0)
            sa = (a["lo"] > n0) - (a["hi"] < n0); sb = (b["lo"] > n0) - (b["hi"] < n0)
            # a sign flip counts only when at least one of the two intervals excludes the null value
            flip = np.sign(a["est"] - n0) != np.sign(b["est"] - n0) and (sa != 0 or sb != 0)
            if sa != sb or flip or \
               (a.get("p") is not None and b.get("p") is not None and (a["p"] < 0.05) != (b["p"] < 0.05)): cc.append(C)
        elif kind == "p" and a.get("p") is not None and b.get("p") is not None:
            if (a["p"] < 0.05) != (b["p"] < 0.05): cc.append(C)
        elif kind == "rank" and "rankkey" in a:
            if a["rankkey"] != b.get("rankkey"): cc.append(C)
        elif kind == "manual" and b.get("change"): cc.append(C)
    ch = "no" if not cc else "yes (" + ", ".join(cc) + ")"
    rows.append(dict(section=section, quantity=quantity, text_value=text, v2_record=v2s, v3=vals["v3"]["s"], v3b=vals["v3b"]["s"], script=script, conclusion_change=ch))


DAG, TST, SNS = "crater_dag.py", "crater_dag_tests.py", "crater_dag_sens.py"
TH, S3, S9, S4 = "summary_model/step7_improved.py", "summary_model/step3_dominance.py", "summary_model/step9_color.py", "summary_model/step4_trees.py"
NM = {"age": "degradation state", "insolation": "insolation", "size": "crater size", "nesting": "nesting", "cover": "volcanic cover",
      "volcanism": "vents", "composition": "low-reflectance material", "detection": "detection"}

# ======================= Section 2.2 (definitions used in 3.3-3.8)
row("2.2", "host craters (any assigned hollow) among 31,571 craters", "707",
    lambda d, C: V(int(pd.read_csv(f"{H}/inputs/test2_craters_h18_{C}.csv", usecols=["hosts"]).hosts.sum())), "inputs/make_inputs.py")
row("2.2", "h09 hosts in the containment tests (hollow within 0.9 radii)", "606", lambda d, C: V(d["tests"]["nested_h09"]["events"]), TST)
row("2.2", "any-hollow hosts used for the AIC comparisons", "705", lambda d, C: V(d["dag"]["counts"]["hosts"]), DAG)

# ======================= Section 3.3
S = "3.3"
row(S, "nested vs isolated, odds ratio (h09, with regional fraction)", "1.85 (1.54-2.22)", lambda d, C: sensOR(d, "all"), SNS, "ratio")
def subset_range(d, C):
    ks = ["all", "C1_ge_0.5", "D_10_50", "D_ge10_preservation_controlled", "abslat_le60"]
    xs = [d["sens"][k]["nested"]["OR"] for k in ks] + [float(np.exp(d["tests"]["nested_local_h09_3"]["params"]["nested"][0]))]
    lows = [d["sens"][k]["nested"]["ci"][0] for k in ks]
    return V(f"{min(xs):.2f}-{max(xs):.2f} (every lower bound > 1: {min(lows) > 1})", change=not (min(lows) > 1))
row(S, "range of nesting odds ratios over subsets that include all parents (Table S18)", "1.52-1.93", subset_range, f"{SNS}, {TST}", "manual")
def matched(d, C):
    m = d["sens"]["matched_300km"]; return V(f"{m['nested_with_matches']:,}; {pct(m['nested_rate'], 2)} vs {pct(m['matched_isolated_rate'], 2)}; ratio {m['ratio']:.2f}", est=m["ratio"])
row(S, "nested vs matched isolated craters within 300 km (no regression)", "about twice as often (2.1)", matched, SNS)
for lab, t, txt in [("inside craters < 300 km", "C(pk)[T.crater_lt300]", "2.89 (2.28-3.65)"), ("inside craters >= 300 km", "C(pk)[T.crater_ge300]", "1.56 (1.16-2.10)"),
                    ("inside basins", "C(pk)[T.basin]", "1.39 (1.10-1.78)")]:
    row(S, f"nesting odds ratio {lab} (fit by parent type)", txt, lambda d, C, t=t: tparam(d, "parentkind_local_h09", t), TST, "ratio")
row(S, "one ancestor vs isolated", "1.78 (1.48-2.16)", lambda d, C: tparam(d, "chain_local_h09", "C(chain.clip(upper=2))[T.1]"), TST, "ratio")
row(S, "two or more ancestors vs isolated", "2.36 (1.67-3.32)", lambda d, C: tparam(d, "chain_local_h09", "C(chain.clip(upper=2))[T.2]"), TST, "ratio")
row(S, "inheritance: parent hosts (own hollow within 0.9 R) vs not", "2.23 (1.62-3.07)", lambda d, C: sensOR(d, "inherit_all", "parent_h09"), SNS, "ratio")
row(S, "inheritance, degradation state controlled (1,295 craters >= 10 km)", "1.68 (0.99-2.86; p = 0.055)", lambda d, C: sensOR(d, "inherit_preservation_controlled", "parent_h09"), SNS, "ratio")
row(S, "inheritance, parent counts as host if any hollow is assigned", "1.96 (1.42-2.70)", lambda d, C: tparam(d, "inherit_h09", "parent_hosts"), TST, "ratio")
def irates(d, C):
    (a, b), (c, e) = rates(d["tests"]["inherit_rates_h09"])
    return V(f"{a:,} of {b:,} ({pct(a/b, 1)}) vs {c:,} of {e:,} ({pct(c/e, 2)})", est=(a / b) / (c / e))
row(S, "hosting inside hollow-bearing vs other parents (any hollow)", "139 of 2,706 (5.1%) vs 74 of 7,386 (1.0%)", irates, TST)
row(S, "inheritance among 365 craters inside parents with completeness >= 0.5", "4.64 (1.63-13.2)", lambda d, C: sensOR(d, "inherit_C1_ge_0.5", "parent_h09"), SNS, "ratio")
row(S, "hosts among those 365 craters", "38", lambda d, C: V(d["sens"]["inherit_C1_ge_0.5"]["events"]), SNS)
def clog(d, C):
    b, se, p = d["tests"]["within_parent_clogit_h09"]["params"]["inner"]; return OR(np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p)
row(S, "peak/floor vs wall, same parent (conditional logit)", "1.23 (0.89-1.70; p = 0.20)", clog, TST, "ratio")
row(S, "fill-only craters vs isolated", "2.09 (1.36-3.21)", lambda d, C: tparam(d, "zone_local_h09", "C(zone_c)[T.fill_only]"), TST, "ratio")
row(S, "AIC drop from nesting (no regional fraction; any-hollow hosts)", "61", lambda d, C: V(f"{d['dag']['base']['aic'] - d['dag']['nested']['aic']:.1f}"), DAG)
def zz(R): return V(f"{R['ln_zown']['aic'] - R['ln_zown_zmax']['aic']:.1f} (p = {fp(R['ln_zown_zmax']['params']['ln_zmax'][2])})", p=R['ln_zown_zmax']['params']['ln_zmax'][2])
row(S, "AIC drop from adding ln z_max to ln z_own, terrestrial (p)", "11 (p = 0.0002)", lambda d, C: zz(d["dag"]), DAG, "p")
row(S, "AIC drop from adding ln z_max to ln z_own, lunar (p)", "11 (p = 0.0002)", lambda d, C: zz(d["lunar"]), DAG, "p")
def bdmn(d, C):
    a, b = best_depth_minus_nested(d["dag"]), best_depth_minus_nested(d["lunar"])
    return V(f"{min(a, b):.1f} (terrestrial {a:.1f}, lunar {b:.1f})", change=min(a, b) <= 10)
row(S, "minimum AIC excess of a depth-of-origin model over nesting alone", "at least 56", bdmn, DAG, "manual")

# ======================= Figure 3 caption
PE = "parents_and_examples.py"
row("Fig. 3", "craters inside Praxiteles that host (h09)", "4 of 6", lambda d, C: V(f"{d['pe']['Praxiteles']['children_h09']} of {d['pe']['Praxiteles']['children']}"), PE, norecord=True)
row("Fig. 3", "craters inside Lermontov that host (h09)", "4 of 5", lambda d, C: V(f"{d['pe']['Lermontov']['children_h09']} of {d['pe']['Lermontov']['children']}"), PE, norecord=True)
row("Fig. 3 / 3.4", "h09 hosting rate of all craters >= 5 km in the tests", "1.9%", lambda d, C: V(pct(d["pe"]["h09_rate_all_tests"]["rate"], 2)), PE, norecord=True)

# ======================= Section 3.4 vents
S = "3.4"
def t8(key): return lambda d, C: OR(d["t8"][key]["odds_ratio"], d["t8"][key]["lo"], d["t8"][key]["hi"], d["t8"][key]["p"])
row(S, "vent craters that host (h09), of 148", "40 (27%)", lambda d, C: V(f"{d['t8']['paper_spec_counts']['vent_h09_hosts']} ({pct(d['t8']['paper_spec_counts']['vent_h09_hosts'] / 148)})"), "vents/test8_vents.py")
row(S, "vents: all craters", "4.88 (2.89-8.25)", t8("paperspec_all"), "vents/test8_vents.py", "ratio")
row(S, "vents: isolated craters", "3.37 (1.68-6.73)", t8("paperspec_isolated"), "vents/test8_vents.py", "ratio")
row(S, "vents: craters >= 20 km", "3.63 (2.08-6.32)", t8("paperspec_D_ge20"), "vents/test8_vents.py", "ratio")

# ======================= Section 3.6 (hosting only)
S = "3.6"
HA, KT = "degradation/host_age_tests.py", "degradation/kinczyk_tests.py"
def hage(k): return lambda d, C: OR(*d["hage"][k][:3], d["hage"][k][3])
row(S, "Herrick: fresh vs degraded (craters >= 10 km, any-hollow hosts)", "4.9 (3.4-7.0)", hage("OR_f"), HA, "ratio")
row(S, "Herrick: standard vs degraded", "1.8 (1.4-2.4)", hage("OR_s"), HA, "ratio")
def krates(d, C):
    r = d["kin"]["kinczyk_counts_hosts"]; xs = [r[str(c)][2] for c in range(1, 6)]
    return V(", ".join(f"{100*x:.1f}%" for x in xs) + ("" if all(xs[i + 1] > xs[i] for i in range(4)) else " (not monotonic)"), change=any(xs[i + 1] <= xs[i] for i in range(4)))
row(S, "Kinczyk hosting fraction by Class 1-5 (h09)", "3.9%, 5.5%, 16%, 27%, 41% (rising)", krates, KT, "manual")
row(S, "Kinczyk: odds ratio of hosting per class (log D, |lat|)", "2.8 (2.4-3.2)", lambda d, C: kc(d, "kinczyk_paper_spec", "per_class_step"), KT, "ratio")
row(S, "host craters of Classes 1-3 with Stage 1 hollows (incidence < 70 deg)", "75", lambda d, C: V(d["stage"]["inc_lt_70"]["hosts_with_stage1_class1to3"]), "degradation/stage_by_class.py")
def inter(d, C):
    a = d["hage"]["interaction_p"]; b = d["kin"]["dT_interaction_p"]
    return V(f"states p = {a['dT:C(pres)[T.s]']:.2f}, {a['dT:C(pres)[T.f]']:.2f}; classes p = {b['dT10:C(grp)[T.C3]']:.2f}, {b['dT10:C(grp)[T.C4-5]']:.2f}",
             change=min(list(a.values()) + list(b.values())) < 0.05)
row(S, "rise of hosting with Delta T_w does not differ between states or classes", "no difference", inter, f"{HA}, {KT}", "manual")

# ======================= Section 3.8 joint model
S = "3.8"
row(S, "fitted detection exponent g", "0.88", lambda d, C: V(f"{d['imp']['thinned_presence_all']['coefs']['lnC']:.3f}"), TH)
row(S, "McFadden R^2 of the seven groups beyond detection and area", "0.25 (0.22-0.29)",
    lambda d, C: (lambda T: V(f"{T['point']['full']:.3f} ({T['ci']['full'][0]:.3f}-{T['ci']['full'][1]:.3f})"))(d["imp"]["thinned_presence_all"]), TH)
for g, txt in [("age", "31% (25-37%)"), ("insolation", "25% (20-31%)"), ("size", "17% (14-21%)"), ("nesting", "16% (13-19%)"), ("cover", "7% (4-10%)"),
               ("volcanism", "3% (2-6%)"), ("composition", "1%")]:
    row(S, f"dominance share (thinned occurrence): {NM[g]} [rank]", txt,
        lambda d, C, g=g: V(share(d, "thinned_presence_all", g) + f" [rank {rank(d['imp']['thinned_presence_all']['point']['share'], g)}]",
                            rankkey=rank(d["imp"]["thinned_presence_all"]["point"]["share"], g)), TH, "rank")
row(S, "order of the groups (thinned occurrence)", "age > insolation > size > nesting > cover > vents > LRM",
    lambda d, C: V(order(d["imp"]["thinned_presence_all"]["point"]["share"]), rankkey=order(d["imp"]["thinned_presence_all"]["point"]["share"])), TH, "rank")
def same_order(d, C):
    a = order(d["imp"]["logistic_C_ge_0.2"]["share"], ("detection",)); b = order(d["imp"]["thinned_presence_all"]["point"]["share"])
    return V(f"{a} (same as thinned: {a == b})", rankkey=(a == b))
row(S, "order is the same for craters with completeness >= 0.2", "same order", same_order, TH, "rank")
for g, txt in [("age", "46%"), ("nesting", "20%")]:
    row(S, f"count model: {NM[g]} share [rank]", txt, lambda d, C, g=g: (lambda sh: V(f"{100*sh[g]:.1f}% [rank {rank(sh, g)}]", rankkey=rank(sh, g)))(d["cnt"]["share"]), TH, "rank")
row(S, "low-reflectance share with the XRS ratios added", "3%", lambda d, C: V(f"{100*d['imp']['thinned_presence_xrs_subset_lrm_xrs']['point']['share']['composition']:.1f}%"), TH)
def colshare(d, C):
    v = d["col"]["variants"]["lrm_pc"]; p, (lo, hi) = v["point"]["share"]["composition"], v["ci_share"]["composition"]
    return V(f"{100*p:.1f}% ({100*lo:.1f}-{100*hi:.1f}%)", est=p, lo=lo, hi=hi, null=0.0)
row(S, "low-reflectance share with MDIS color PC1, PC2 added", "6% (3-10%)", colshare, S9, "coef")
def other_order(d, C):
    base = order(d["imp"]["thinned_presence_all"]["point"]["share"], ("composition",))
    x = order(d["imp"]["thinned_presence_xrs_subset_lrm_xrs"]["point"]["share"], ("composition",)); c = order(d["col"]["variants"]["lrm_pc"]["point"]["share"], ("composition",))
    return V(f"XRS: {x}; color: {c}", rankkey=(x == base, c == base))
row(S, "order of the other groups unchanged with XRS or color", "unchanged", other_order, f"{TH}, {S9}", "rank")
row(S, "R^2 added by the regional hosting fraction", "0.03", lambda d, C: V(f"{d['imp']['thinned_presence_all']['point']['region_add']:.4f}"), TH)
row(S, "reduction of that gain by XRS or by color", "0.001-0.002",
    lambda d, C: V(f"XRS {d['imp']['region_add_reduction_by_xrs']['point']:.4f}; color {d['col']['region_add_reduction_by_color']['point']:.4f}"), f"{TH}, {S9}")

# ======================= Discussion, Conclusions, Table 1, abstract
S = "4.1 Discussion"
row(S, "fill-only craters host more often ('the fill may itself be charged regolith')", "2.09 (1.36-3.21)", lambda d, C: tparam(d, "zone_local_h09", "C(zone_c)[T.fill_only]"), TST, "ratio")
def fifth(d, C):
    (a, b), (c, e) = rates(d["tests"]["inherit_rates_h09"]); return V(f"{pct(c/e, 2)}; 1/{(a/b)/(c/e):.1f} of {pct(a/b, 1)}")
row(S, "hosting rate inside parents without hollows; fraction of the rate in hollow-bearing parents", "1.0%, a fifth", fifth, TST)
row(S, "parents without observed hollows whose children host; their median completeness", "75; 0.08",
    lambda d, C: V(f"{d['pe']['parents_without_hollows']['host_no_hollows']}; {d['pe']['parents_without_hollows']['median_C1_no']:.2f}"), PE, norecord=True)
row("4.3 Discussion", "within a parent, wall craters host nearly as often as peak/floor craters", "1.23 (0.89-1.70)", clog, TST, "ratio")
T1 = "Table 1 / Abstract / Conclusions"
row(T1, "craters inside older craters host hollows more often", "yes", lambda d, C: sensOR(d, "all"), SNS, "ratio")
row("Table 1", "craters inside hollow-bearing craters host most often", "yes", lambda d, C: sensOR(d, "inherit_all", "parent_h09"), SNS, "ratio")
row(T1, "fresh craters host hollows more often", "yes", hage("OR_f"), HA, "ratio")
row(T1, "craters with explosive vents host hollows more often", "yes", t8("paperspec_all"), "vents/test8_vents.py", "ratio")
row("Conclusions", "degradation state and insolation have the largest shares (thinned joint model)", "yes",
    lambda d, C: (lambda o: V(" > ".join(o[:3]), rankkey=set(o[:2]) == {"age", "insolation"}))(order(d["imp"]["thinned_presence_all"]["point"]["share"]).split(" > ")), TH, "rank")

# ======================= S11 degradation (Tables S13, S14)
S = "S11 / Table S13"
for st, lab, n_txt, dt_txt in [("d", "degraded", "114 of 1,736", "1.34 (1.12-1.60; 106)"), ("s", "standard", "151 of 2,230", "1.21 (1.05-1.39; 146)"), ("f", "fresh", "62 of 250", "1.34 (1.02-1.75; 56)")]:
    row(S, f"hosts among {lab} craters (any hollow)", n_txt, lambda d, C, st=st: V(f"{d['hage']['counts'][st]['hosts']} of {d['hage']['counts'][st]['n']:,}"), HA)
    row(S, f"odds per 10 K of Delta T_w within {lab} craters (hosts)", dt_txt,
        lambda d, C, st=st: (lambda v: dict(OR(*v[:3]), s=OR(*v[:3])["s"][:-1] + f"; {v[3]})"))(d["hage"]["ORper10K_" + st]), HA, "ratio")
row(S, "odds ratio of hosting, standard vs degraded", "1.83 (1.39-2.41)", hage("OR_s"), HA, "ratio")
row(S, "odds ratio of hosting, fresh vs degraded", "4.86 (3.35-7.05)", hage("OR_f"), HA, "ratio")
row(S, "Delta T_w x state interaction p (standard, fresh)", "0.27 and 0.86",
    lambda d, C: V(f"{d['hage']['interaction_p']['dT:C(pres)[T.s]']:.2f} and {d['hage']['interaction_p']['dT:C(pres)[T.f]']:.2f}", p=min(d["hage"]["interaction_p"].values())), HA, "p")
S = "S11 / Table S14 (Kinczyk)"
row(S, "craters >= 40 km with a class 1-5 in the tests (h09 hosts)", "2,968 (258)", lambda d, C: V(f"{d['kin']['kinczyk_paper_spec']['n']:,} ({d['kin']['kinczyk_paper_spec']['hosts']})"), KT)
row(S, "per-class odds ratio, matches with offset <= 0.3 radius", "2.76 (2.38-3.19)", lambda d, C: OR(*d["msens"]["r_le_0.3 | logD + abslat"][2:5]), "degradation/kinczyk_match_sensitivity.py", "ratio")
row(S, "per-class odds ratio (log D, |lat|)", "2.8 (2.4-3.2)", lambda d, C: kc(d, "kinczyk_paper_spec", "per_class_step"), KT, "ratio")
row(S, "per-class odds ratio with completeness", "2.5 (2.1-2.9)", lambda d, C: kc(d, "kinczyk_with_C1", "per_class_step"), KT, "ratio")
row(S, "separate class odds ratios fit better than one factor per class, with completeness (LR p)", "0.024",
    lambda d, C: V(fp(d["kin"]["kinczyk_with_C1"]["ordinal_vs_categorical_LR_p"]), p=d["kin"]["kinczyk_with_C1"]["ordinal_vs_categorical_LR_p"]), KT, "p")
def h2h(sub): return lambda d, C: (lambda h: V(f"{h['AIC_herrick4'] - h['AIC_kinczyk5']:.1f}", rankkey=h["AIC_herrick4"] > h["AIC_kinczyk5"]))(d["kin"]["head_to_head_ge40"][sub])
for sub, lab, txt in [("simple", "diameter and latitude", "4"), ("with_C1", "with completeness", "13"), ("full", "joint-model controls", "8")]:
    row(S, f"AIC(Herrick 4 categories) - AIC(Kinczyk 5 classes), {lab}", txt, h2h(sub), KT, "rank")
def maxlr(d, C):
    h = d["kin"]["head_to_head_ge40"]; p = max(h[s][k] for s in ["simple", "with_C1", "full"] for k in ["LR_kinczyk_given_herrick_p", "LR_herrick_given_kinczyk_p"])
    return V(f"<= {fp(p)}", p=p)
row(S, "each classification improves the fit of the other (largest LR p)", "<= 0.009", maxlr, KT, "p")
row(S, "Herrick-filled craters: odds ratio per class, Classes 1-3", "2.0 (1.4-2.9)", lambda d, C: OR(*d["kin"]["head_to_head_ge40"]["filled_per_class_step_k1to3"]), KT, "ratio")
row(S, "Herrick-filled craters that host, Classes 1, 2, 3", "36 of 933, 51 of 1,035, 13 of 97",
    lambda d, C: V(", ".join(f"{d['kin']['head_to_head_ge40']['filled_by_kclass'][str(c)][1]} of {d['kin']['head_to_head_ge40']['filled_by_kclass'][str(c)][0]:,}" for c in [1, 2, 3])), KT)
for g, txt in [("C1-2", "1.21 (1.04-1.42)"), ("C3", "1.33 (1.08-1.64)"), ("C4-5", "1.24 (0.97-1.58)")]:
    row(S, f"odds per 10 K of Delta T_w, Classes {g[1:]} (|lat| <= 60 deg)", txt, lambda d, C, g=g: OR(*d["kin"][f"dT_per10K_{g}"][:4]), KT, "ratio")
row(S, "Delta T_w x class interaction p", "0.66 and 0.92",
    lambda d, C: V(f"{d['kin']['dT_interaction_p']['dT10:C(grp)[T.C3]']:.2f} and {d['kin']['dT_interaction_p']['dT10:C(grp)[T.C4-5]']:.2f}", p=min(d["kin"]["dT_interaction_p"].values())), KT, "p")
row(S, "matched ghost craters that host", "2 of 80", lambda d, C: V(f"{d['kin']['code6_counts_hosts'][1]} of {d['kin']['code6_counts_hosts'][0]}"), KT)
row(S, "host craters with Stage 1 hollows, Classes 1-5 (incidence < 70 deg)", "9, 23, 43, 35, 5",
    lambda d, C: V(", ".join(str(d["stage"]["inc_lt_70"]["hosts_with_stage1_by_class"].get(str(c), 0)) for c in range(1, 6))), "degradation/stage_by_class.py")
TXT14 = {1: ("36", "3.9%", "1", "1"), 2: ("71", "5.5%", "2.2 (1.4-3.4)", "2.0 (1.2-3.1)"), 3: ("83", "16%", "8.7 (5.5-13.9)", "8.1 (5.0-13.0)"),
         4: ("59", "27%", "18 (11-31)", "12 (7-21)"), 5: ("9", "41%", "37 (14-98)", "18 (6-58)")}
for c in range(1, 6):
    row(S, f"Table S14, Class {c}: hosts (hosting fraction)", f"{TXT14[c][0]} ({TXT14[c][1]})",
        lambda d, C, c=c: V(f"{d['kin']['kinczyk_counts_hosts'][str(c)][1]} ({100*d['kin']['kinczyk_counts_hosts'][str(c)][2]:.1f}%)"), KT)
    if c > 1:
        row(S, f"Table S14, Class {c} vs 1: odds ratio", TXT14[c][2], lambda d, C, c=c: kc(d, "kinczyk_paper_spec", f"class{c}_vs_1"), KT, "ratio")
        row(S, f"Table S14, Class {c} vs 1: odds ratio with completeness", TXT14[c][3], lambda d, C, c=c: kc(d, "kinczyk_with_C1", f"class{c}_vs_1"), KT, "ratio")
row(S, "odds rise less from Class 4 to 5 than between the other classes (with completeness)", "yes (12 -> 18)",
    lambda d, C: V(f"{d['kin']['kinczyk_with_C1']['class4_vs_1'][0]} -> {d['kin']['kinczyk_with_C1']['class5_vs_1'][0]}" + (" (falls)" if d['kin']['kinczyk_with_C1']['class5_vs_1'][0] < d['kin']['kinczyk_with_C1']['class4_vs_1'][0] else ""),
                   change=d['kin']['kinczyk_with_C1']['class5_vs_1'][0] < d['kin']['kinczyk_with_C1']['class4_vs_1'][0]), KT, "manual")

# ======================= S14 containment graph (Tables S18, S19)
S = "S14"
row(S, "host craters in the tests; nested hosts", "705 of 707; 385 (55%)",
    lambda d, C: V(f"{d['dag']['counts']['hosts']}; {d['dag']['counts']['hosts_nested']} ({pct(d['dag']['counts']['hosts_nested'] / d['dag']['counts']['hosts'])})"), DAG)
row(S, "matched comparison: nested craters with matches; rates; ratio", "12,055; 1.6% vs 0.76%; 2.1", matched, SNS)
def orgel(k): return lambda d, C: (lambda v: OR(v["OR"], v["ci"][0], v["ci"][1], v["p"]))(d["orgel"][k]["parent_h09"])
row(S, "inheritance with Caloris as the basin record", "2.36 (1.71-3.25)", orgel("inherit_all"), "crater_dag*.py with CRATER_DAG_CALORIS=orgel", "ratio")
row(S, "same, degradation state controlled", "2.05 (1.21-3.48)", orgel("inherit_preservation_controlled"), "crater_dag*.py with CRATER_DAG_CALORIS=orgel", "ratio")
def caldep(d, C):
    a = d["sens"]["inherit_preservation_controlled"]["parent_h09"]["p"]; b = d["orgel"]["inherit_preservation_controlled"]["parent_h09"]["p"]
    return V(f"p = {fp(a)} (catalog record) vs {fp(b)} (basin record)", rankkey=((a < 0.05) == (b < 0.05)))
row(S, "'the ratio with degradation state controlled depends on how Caloris is treated'", "yes (p = 0.055 vs 0.008)", caldep, SNS, "rank")
def clog_counts(d, C):
    v = d["tests"]["within_parent_clogit_h09"]; (a, b), (c, e) = [tuple(int(x) for x in v["rates"][k].split("/")) for k in ["inner", "wall"]]
    return V(f"{v['parents']} parents, {v['children']:,} children, {v['events']} hosts; {a} of {b:,} ({pct(a/b, 1)}) vs {c} of {e:,} ({pct(c/e, 1)})")
row(S, "conditional logit sample and raw rates (peak/floor vs wall)", "118 parents, 5,160 children, 283 hosts; 100 of 1,392 (7.2%) vs 183 of 3,768 (4.9%)", clog_counts, TST)
row(S, "AIC drop from z_base, terrestrial (p); lunar", "5 (p = 0.007), either relation",
    lambda d, C: V(f"{d['dag']['base']['aic'] - d['dag']['z_base']['aic']:.1f} (p = {fp(d['dag']['z_base']['params']['z_base'][2])}); lunar {d['lunar']['base']['aic'] - d['lunar']['z_base']['aic']:.1f}", p=d["dag"]["z_base"]["params"]["z_base"][2]), DAG, "p")
row(S, "k_ring 0.8/1.3: largest change of (best depth model - nested) AIC", "<= 0.3",
    lambda d, C: V(f"{max(abs(best_depth_minus_nested(d[k]) - best_depth_minus_nested(d['dag'])) for k in ['k08', 'k13']):.1f}"), DAG)
row(S, "k_ring 0.8/1.3: largest change of the AIC gained by adding ln z_max", "<= 0.1",
    lambda d, C: V(f"{max(abs((d[k]['ln_zown']['aic'] - d[k]['ln_zown_zmax']['aic']) - (d['dag']['ln_zown']['aic'] - d['dag']['ln_zown_zmax']['aic'])) for k in ['k08', 'k13']):.1f}"), DAG)
row(S, "hosts with z_max > 35 km and > 53 km (terrestrial; lunar in parentheses)", "84 and 46 of 705 (85 and 47)",
    lambda d, C: V(f"{d['dag']['crust_base']['hosts']['z_max_gt_35']} and {d['dag']['crust_base']['hosts']['z_max_gt_53']} of {d['dag']['crust_base']['hosts']['n']} ({d['lunar']['crust_base']['hosts']['z_max_gt_35']} and {d['lunar']['crust_base']['hosts']['z_max_gt_53']})"), DAG)
row(S, "parent craters with a hosting child; of them without hollows", "149; 75",
    lambda d, C: V(f"{d['pe']['parents_without_hollows']['parents_with_hosting_child']}; {d['pe']['parents_without_hollows']['host_no_hollows']}"), PE, norecord=True)
row(S, "median completeness, parents without vs with hollows; parents without hollows with C >= 0.2, >= 0.5", "0.08 vs 0.21; 22, 4",
    lambda d, C: (lambda p: V(f"{p['median_C1_no']:.2f} vs {p['median_C1_hosting_parents']:.2f}; {p['no_C1_ge_0p2']}, {p['no_C1_ge_0p5']}"))(d["pe"]["parents_without_hollows"]), PE, norecord=True)
row(S, "depth ratio, hosts inside hollow-bearing parents vs isolated hosts", "1.05 (0.91-1.22)",
    lambda d, C: V("1.05 (0.91-1.22)" + ("" if C == "v2" else "; unchanged (parents of the 128 measured hosts keep their flags)")), "depth_nested/depth_nested.py")
S = "S14 / Table S18"
def cnt(d, key, src="sens"):
    v = d[src][key]; return f"{v['n']:,} ({v['events']}) "
for q, g, t, sc in [
        ("all craters", lambda d, C: dict(sensOR(d, "all"), s=cnt(d, "all") + sensOR(d, "all")["s"]), "31,554 (606) 1.85 (1.54-2.22)", SNS),
        (">= 3 hollows within 0.9 radii", lambda d, C: dict(tparam(d, "nested_local_h09_3", "nested"), s=f"({d['tests']['nested_local_h09_3']['events']}) " + tparam(d, "nested_local_h09_3", "nested")["s"]), "31,554 (198) 1.59 (1.13-2.25)", TST),
        ("completeness >= 0.5", lambda d, C: dict(sensOR(d, "C1_ge_0.5"), s=cnt(d, "C1_ge_0.5") + sensOR(d, "C1_ge_0.5")["s"]), "1,119 (114) 1.76 (1.08-2.86)", SNS),
        ("craters 10-50 km", lambda d, C: dict(sensOR(d, "D_10_50"), s=cnt(d, "D_10_50") + sensOR(d, "D_10_50")["s"]), "14,582 (239) 1.93 (1.46-2.54)", SNS),
        ("craters >= 10 km, degradation state controlled", lambda d, C: dict(sensOR(d, "D_ge10_preservation_controlled"), s=cnt(d, "D_ge10_preservation_controlled") + sensOR(d, "D_ge10_preservation_controlled")["s"]), "4,214 (303) 1.52 (1.14-2.03)", SNS),
        ("within 60 deg of the equator", lambda d, C: dict(sensOR(d, "abslat_le60"), s=cnt(d, "abslat_le60") + sensOR(d, "abslat_le60")["s"]), "27,257 (573) 1.84 (1.52-2.23)", SNS),
        ("craters inside craters >= 300 km and basins dropped", lambda d, C: dict(sensOR(d, "parents_lt300_only"), s=cnt(d, "parents_lt300_only") + sensOR(d, "parents_lt300_only")["s"]), "23,212 (413) 3.10 (2.43-3.97)", SNS),
        ("position: peak", lambda d, C: tparam(d, "zone_local_h09", "C(zone_c)[T.peak]"), "2.55 (1.38-4.70)", TST),
        ("position: floor", lambda d, C: tparam(d, "zone_local_h09", "C(zone_c)[T.floor]"), "2.28 (1.75-2.97)", TST),
        ("position: wall", lambda d, C: tparam(d, "zone_local_h09", "C(zone_c)[T.wall]"), "1.64 (1.33-2.01)", TST),
        ("parent hosts (own hollow): craters (hosts)", lambda d, C: dict(sensOR(d, "inherit_all", "parent_h09"), s=cnt(d, "inherit_all") + sensOR(d, "inherit_all", "parent_h09")["s"]), "10,092 (213) 2.23 (1.62-3.07)", SNS),
        ("parent hosts, degradation state controlled: craters (hosts)", lambda d, C: dict(sensOR(d, "inherit_preservation_controlled", "parent_h09"), s=cnt(d, "inherit_preservation_controlled") + sensOR(d, "inherit_preservation_controlled", "parent_h09")["s"]), "1,295 (100) 1.68 (0.99-2.86; p = 0.055)", SNS),
        ("parent hosts, completeness >= 0.5: craters (hosts)", lambda d, C: dict(sensOR(d, "inherit_C1_ge_0.5", "parent_h09"), s=cnt(d, "inherit_C1_ge_0.5") + sensOR(d, "inherit_C1_ge_0.5", "parent_h09")["s"]), "365 (38) 4.64 (1.63-13.2)", SNS),
        ("parent hosts (any hollow): craters (hosts)", lambda d, C: dict(tparam(d, "inherit_h09", "parent_hosts"), s=cnt(d, "inherit_h09", "tests") + tparam(d, "inherit_h09", "parent_hosts")["s"]), "10,092 (213) 1.96 (1.42-2.70)", TST)]:
    row(S, q, t, g, sc, "ratio")
S = "S14 / Table S19"
row(S, "AIC base", "5,102.0", lambda d, C: aicrow(d, "dag", "base"), DAG)
row(S, "AIC base + nested; coefficient", "5,041.2; 0.68 +- 0.09 (p = 4e-15)", lambda d, C: aicrow(d, "dag", "nested", "nested", "{:.2f}"), DAG, "coef")
row(S, "AIC base + position", "5,041.5", lambda d, C: aicrow(d, "dag", "zone"), DAG)
row(S, "AIC base + number of ancestors", "5,036.7", lambda d, C: aicrow(d, "dag", "chain"), DAG)
row(S, "AIC base + z_base; coefficient", "5,097.1; 0.0059 +- 0.0022 (p = 0.007)", lambda d, C: aicrow(d, "dag", "z_base", "z_base"), DAG, "coef")
row(S, "AIC ln z_own in place of log D; coefficient", "5,130.9; 1.11 +- 0.04", lambda d, C: aicrow(d, "dag", "ln_zown", "ln_zown", "{:.2f}"), DAG, "coef")
row(S, "AIC ln z_max in place of log D; coefficient", "5,589.8; 0.39 +- 0.02", lambda d, C: aicrow(d, "dag", "ln_zmax", "ln_zmax", "{:.2f}"), DAG, "coef")
row(S, "AIC ln z_own + ln z_max; ln z_max coefficient", "5,119.6; 0.13 +- 0.03 (p = 0.0002)", lambda d, C: aicrow(d, "dag", "ln_zown_zmax", "ln_zmax", "{:.2f}"), DAG, "coef")
MODS = ["base", "nested", "zone", "chain", "z_base", "ln_zown", "ln_zmax", "ln_zown_zmax"]
def aicrank(rec): return lambda d, C: (lambda o: V(o, rankkey=o))(" < ".join(sorted(MODS, key=lambda m: d[rec][m]["aic"])))
row(S, "AIC ranking of the eight models (terrestrial)", "chain < nested < zone < z_base < base < ln_zown_zmax < ln_zown < ln_zmax", aicrank("dag"), DAG, "rank")
row(S, "caption (lunar): AIC of position model; z_base, ln z_own, ln z_max, both", "5,041.6; 5,097.2, 5,149.0, 5,552.9, 5,138.0",
    lambda d, C: V(f"{d['lunar']['zone']['aic']:.1f}; " + ", ".join(f"{d['lunar'][m]['aic']:.1f}" for m in ["z_base", "ln_zown", "ln_zmax", "ln_zown_zmax"])), DAG)
row(S, "caption (lunar): added terms", "0.0059 +- 0.0022 (p = 0.007), 0.97 +- 0.04, 0.41 +- 0.02, 0.13 +- 0.03 (p = 0.0002)",
    lambda d, C: V(", ".join(f"{d['lunar'][m]['params'][t][0]:.4g} +- {d['lunar'][m]['params'][t][1]:.2g}" for m, t in [("z_base", "z_base"), ("ln_zown", "ln_zown"), ("ln_zmax", "ln_zmax"), ("ln_zown_zmax", "ln_zmax")])
                   + f" (p = {fp(d['lunar']['z_base']['params']['z_base'][2])}, {fp(d['lunar']['ln_zown_zmax']['params']['ln_zmax'][2])})"), DAG)
row(S, "AIC ranking of the eight models (lunar)", "chain < nested < zone < z_base < base < ln_zown_zmax < ln_zown < ln_zmax", aicrank("lunar"), DAG, "rank")

# ======================= S15 vents (Table S20)
S = "S15 / Table S20"
row(S, "range of the vent odds ratios in Table S20", "3.4-4.9",
    lambda d, C: (lambda xs: V(f"{min(xs):.2f}-{max(xs):.2f}"))([d["t8"][k]["odds_ratio"] for k in ["paperspec_all", "paperspec_isolated", "paperspec_D_ge20", "paperspec_no_local", "hosting_all"]]), "vents/test8_vents.py")
for k, lab, txt in [("in_crater_yes", "vents flagged in the catalog as inside craters", "4.29 (2.49-7.39)"), ("quality_good_or_reasonable", "good or reasonable image quality", "5.20 (2.93-9.23)"),
                    ("quality_good", "good image quality", "5.82 (2.70-12.55)"), ("all_vents_0.8R", "vent within 0.8 radii", "3.90 (2.25-6.76)")]:
    row(S, f"vents: {lab}", txt, lambda d, C, k=k: (lambda v: OR(v["odds_ratio"], v["lo"], v["hi"], v["p"]))(d["vac"]["test8_paperspec_subsets"][k]), "vents/vent_assignment_check.py", "ratio")
row(S, "vents: all craters, without the regional fraction", "4.04 (2.50-6.52)", t8("paperspec_no_local"), "vents/test8_vents.py", "ratio")
row(S, "vents: any hollow within 1.1 radii, without the regional fraction", "3.85 (2.43-6.12)", t8("hosting_all"), "vents/test8_vents.py", "ratio")
row(S, "caption: h09 hosts among 31,554 craters", "606", lambda d, C: V(d["t8"]["paper_spec_counts"]["h09_hosts"]), "vents/test8_vents.py")

# ======================= S16 joint model (Table S21, Figure S9)
S = "S16"
row(S, "hosts (h09) in the joint model; hosts with C = 0 at the center", "606; 18", lambda d, C: V(f"{d['imp']['thinned_presence_all']['events']}; 18 (unchanged)" if C != "v2" else "606; 18"), TH)
row(S, "degradation share with a >= 10 km indicator in the size group", "31%", lambda d, C: V(f"{d['size']['size_plus_ge10']['share_pct']['age']:.1f}%"), "summary_model/step7_sizestep_check.py")
def ge10(d, C):
    s = d["size"]["D_ge10_only"]; neg = {k: v for k, v in s["share_pct"].items() if v < 0}
    return V(f"{s['n']:,} ({s['hosts']}); {s['share_pct']['age']:.1f}%" + (f"; UNRELIABLE: cloglog fits fail to converge (negative shares {neg})" if neg else ""))
row(S, "craters >= 10 km only: craters (hosts); degradation share", "17,014 (488); 38%", ge10, "summary_model/step7_sizestep_check.py")
def twotop(sh, drop=()): return tuple(k for k in sorted(sh, key=lambda k: -sh[k]) if k not in drop)[:2]
def everyver(d, C):
    tops = {"thinned": twotop(d["imp"]["thinned_presence_all"]["point"]["share"]), "logistic": twotop(d["dom"]["share"], ("detection",)),
            "well-imaged": twotop(d["imp"]["logistic_C_ge_0.2"]["share"], ("detection",)), "trees": twotop(d["trees"]["shap_share"], ("detection",))}
    return V("; ".join(f"{k}: {'/'.join(v)}" for k, v in tops.items()), rankkey=all(set(v) == {"age", "insolation"} for v in tops.values()))
row(S, "degradation state and insolation largest in every version except counts", "yes", everyver, f"{TH}, {S3}, {S4}", "rank")
row(S, "count model: nesting second, insolation third", "yes",
    lambda d, C: (lambda o: V(" > ".join(o[:3]), rankkey=tuple(o[1:3])))(order(d["cnt"]["share"]).split(" > ")), TH, "rank")
row(S, "logistic joint model: vent odds ratio (block bootstrap)", "3.3 (1.3-6.3)", lambda d, C: jm(d, "vent"), S3, "ratio")
row(S, "logistic joint model: northern smooth plains odds ratio", "0.30 (0.15-0.57)", lambda d, C: jm(d, "nsp"), S3, "ratio")
row(S, "R^2 added by the three combinations", "0.0002", lambda d, C: V(f"{d['dom']['point']['inter_add']:.5f}"), S3)
for t, lab, txt in [("ph_x_vent", "parent hosts x vent", "1.7 (0.001-14)"), ("nested_x_dT", "nested x Delta T_w (per 10 K)", "0.93 (0.80-1.10)"), ("vent_x_dT", "vent x Delta T_w", "0.98 (0.56-1.93)")]:
    row(S, f"combination odds ratio: {lab}", txt, lambda d, C, t=t: dict(jm(d, t), s=(lambda r: f"{r.OR:.2f} ({r.OR_lo_block:.3g}-{r.OR_hi_block:.3g})")(d["jm"].loc[t])), S3, "ratio")
def auc(d, C):
    cv = d["trees"]["cv"]; tr = max(v["auc"] for k, v in cv.items() if k.startswith("trees:")); return V(f"{tr:.3f}; {cv['logistic: 8 groups']['auc']:.3f}")
row(S, "spatial-CV AUC: boosted trees (best of 150/300/600 trees); logistic", "0.913; 0.905", auc, S4)
row(S, "R^2 added by the regional hosting fraction, added last (block CI)", "0.031 (0.013-0.047)",
    lambda d, C: (lambda T: V(f"{T['point']['region_add']:.3f} ({T['ci']['region_add'][0]:.3f}-{T['ci']['region_add'][1]:.3f})"))(d["imp"]["thinned_presence_all"]), TH)
row(S, "XRS subset (27,810 craters): low-reflectance share without -> with the four ratios (hosts)", "0.8% -> 2.8%",
    lambda d, C: V(f"{100*d['imp']['thinned_presence_xrs_subset_lrm_only']['point']['share']['composition']:.1f}% -> {100*d['imp']['thinned_presence_xrs_subset_lrm_xrs']['point']['share']['composition']:.1f}% ({d['imp']['thinned_presence_xrs_subset_lrm_only']['events']} hosts)"), TH)
row(S, "XRS subset: regional R^2 gain without -> with XRS", "0.028 -> 0.027",
    lambda d, C: V(f"{d['imp']['thinned_presence_xrs_subset_lrm_only']['point']['region_add']:.4f} -> {d['imp']['thinned_presence_xrs_subset_lrm_xrs']['point']['region_add']:.4f}"), TH)
row(S, "color subset (31,373 craters): hosts (small hosts lost to the hollow mask)", "563 (43)",
    lambda d, C: V(f"{d['col']['events']} ({d['imp']['thinned_presence_all']['events'] - d['col']['events']})"), S9)
def colfull(d, C):
    v = d["col"]["variants"]; p, (lo, hi) = v["lrm_pc"]["point"]["share"]["composition"], v["lrm_pc"]["ci_share"]["composition"]
    return V(f"{100*v['lrm_only']['point']['share']['composition']:.1f}% -> {100*p:.1f}% ({100*lo:.1f}-{100*hi:.1f}%)", est=p, lo=lo, hi=hi, null=0.0)
row(S, "color: composition share LRM only -> LRM + PC1 + PC2 (CI)", "1.5% -> 6.3% (2.9-9.8%)", colfull, S9, "coef")
row(S, "color: gain in composition share, percentage points (paired CI)", "4.8 (2.0-7.8)",
    lambda d, C: (lambda g: V(f"{100*g['point']:.1f} ({100*g['ci'][0]:.1f}-{100*g['ci'][1]:.1f})", est=g["point"], lo=g["ci"][0], hi=g["ci"][1], null=0.0))(d["col"]["composition_share_gain"]), S9, "coef")
row(S, "color: shares of the other groups (order)", "degradation 32, insolation 22, size 15, nesting 15, cover 6.5, vents 3",
    lambda d, C: (lambda sh: V(", ".join(f"{k} {100*sh[k]:.1f}" for k in sorted(sh, key=lambda k: -sh[k]) if k != "composition"), rankkey=order(sh, ("composition",))))(d["col"]["variants"]["lrm_pc"]["point"]["share"]), S9, "rank")
for t, lab, txt in [("pc1_z", "spectral slope, per SD", "1.06 (0.87-1.25)"), ("pc2_z", "curvature, per SD", "0.90 (0.79-1.06)"), ("lrm", "LRM, none to full coverage", "1.34 (1.12-1.66)")]:
    row(S, f"color model: density ratio for {lab}", txt, lambda d, C, t=t: OR(*d["col"]["variants"]["lrm_pc"]["rate_ratios"][t]), S9, "ratio")
row(S, "color: regional R^2 gain without -> with color; reduction (CI)", "0.0307 -> 0.0297; 0.0010 (-0.0009 to 0.0029)",
    lambda d, C: (lambda c: V(f"{c['variants']['lrm_only']['point']['region_add']:.4f} -> {c['variants']['lrm_pc']['point']['region_add']:.4f}; {c['region_add_reduction_by_color']['point']:.4f} ({c['region_add_reduction_by_color']['ci'][0]:.4f} to {c['region_add_reduction_by_color']['ci'][1]:.4f})",
                                est=c["region_add_reduction_by_color"]["point"], lo=c["region_add_reduction_by_color"]["ci"][0], hi=c["region_add_reduction_by_color"]["ci"][1], null=0.0))(d["col"]), S9, "coef")
S = "S16 / Table S21"
T21 = {"age": ("31 (25-37)", "18 (14-21)", "27", "46", "14"), "insolation": ("25 (20-31)", "19 (15-23)", "22", "13", "14"), "size": ("17 (14-21)", "16 (14-19)", "21", "7", "12"),
       "nesting": ("16 (13-19)", "12 (10-14)", "15", "20", "12"), "cover": ("7 (4-10)", "3 (2-4)", "8", "8", "1"), "volcanism": ("3 (2-6)", "4 (2-6)", "7", "4", "0.4"),
       "composition": ("1 (0-4)", "0.2 (0.1-0.9)", "0.3", "2", "5")}
def wi(d): s = d["imp"]["logistic_C_ge_0.2"]["share"]; t = sum(v for k, v in s.items() if k != "detection"); return {k: v / t for k, v in s.items() if k != "detection"}
for g, (t_, l_, w_, c_, tr_) in T21.items():
    row(S, f"{NM[g]}: thinned occurrence share [rank]", t_, lambda d, C, g=g: V(share(d, "thinned_presence_all", g) + f" [{rank(d['imp']['thinned_presence_all']['point']['share'], g)}]", rankkey=rank(d["imp"]["thinned_presence_all"]["point"]["share"], g)), TH, "rank")
    row(S, f"{NM[g]}: logistic share [rank]", l_, lambda d, C, g=g: V(f"{100*d['dom']['share'][g]:.1f} ({100*d['dom']['ci'][g]['share'][0]:.1f}-{100*d['dom']['ci'][g]['share'][1]:.1f}) [{rank(d['dom']['share'], g)}]", rankkey=rank(d["dom"]["share"], g)), S3, "rank")
    row(S, f"{NM[g]}: well-imaged share, detection omitted [rank]", w_, lambda d, C, g=g: V(f"{100*wi(d)[g]:.1f} [{rank(wi(d), g)}]", rankkey=rank(wi(d), g)), TH, "rank")
    row(S, f"{NM[g]}: count model share [rank]", c_, lambda d, C, g=g: (lambda sh: V(f"{100*sh[g]:.1f} [{rank(sh, g)}]", rankkey=rank(sh, g)))(d["cnt"]["share"]), TH, "rank")
    row(S, f"{NM[g]}: trees (grouped SHAP) share [rank]", tr_, lambda d, C, g=g: (lambda sh: V(f"{100*sh[g]:.1f} [{rank(sh, g)}]", rankkey=rank(sh, g)))(d["trees"]["shap_share"]), S4, "rank")
row(S, "detection: logistic share; trees share", "29 (23-33); 41",
    lambda d, C: V(f"{100*d['dom']['share']['detection']:.1f} ({100*d['dom']['ci']['detection']['share'][0]:.1f}-{100*d['dom']['ci']['detection']['share'][1]:.1f}); {100*d['trees']['shap_share']['detection']:.1f}"), f"{S3}, {S4}")
row(S, "R^2: thinned; logistic; well-imaged; counts", "0.25 (0.22-0.29); 0.34 (0.31-0.39); 0.33; 0.09",
    lambda d, C: V(f"{d['imp']['thinned_presence_all']['point']['full']:.3f} ({d['imp']['thinned_presence_all']['ci']['full'][0]:.3f}-{d['imp']['thinned_presence_all']['ci']['full'][1]:.3f}); "
                   f"{d['dom']['point']['full']:.3f} ({d['dom']['ci']['full'][0]:.3f}-{d['dom']['ci']['full'][1]:.3f}); {d['imp']['logistic_C_ge_0.2']['point']['full']:.3f}; {d['cnt']['full']:.3f}"), f"{TH}, {S3}")
row(S, "caption: well-imaged craters (hosts)", "3,835 (276)", lambda d, C: V(f"{d['imp']['logistic_C_ge_0.2']['n']:,} ({d['imp']['logistic_C_ge_0.2']['events']})"), TH)
row(S, "count model: negative binomial alpha (not in the text; check)", "-", lambda d, C: V(f"{d['cnt_alpha']:.2f}"), f"{TH}; v3, v3b: summary_model/count_model_check.py")
S = "S16 / Figure S9b"
for t, lab, txt in [("lrm", "LRM fraction", "1.48 (1.03-2.28)"), ("s_si_z", "S/Si per SD", "1.49 (0.97-2.25)"), ("mg_si_z", "Mg/Si per SD", "1.00 (0.77-1.41)"),
                    ("ca_si_z", "Ca/Si per SD", "0.84 (0.62-1.09)"), ("al_si_z", "Al/Si per SD", "1.11 (0.83-1.43)")]:
    row(S, f"density ratio: {lab} (XRS subset; plotted, from xrs_rate_ratios.csv)", txt, lambda d, C, t=t: OR(d["xrr"].loc[t, "rate_ratio"], d["xrr"].loc[t, "lo"], d["xrr"].loc[t, "hi"]),
        "summary_model/step8_figure_improved_v2.py", "ratio")
row("S16 / Figure S9a", "bars = Table S21 thinned, well-imaged and count shares", "see Table S21 rows",
    lambda d, C: V("figures/figS8_joint_model_v2.png" if C == "v2" else f"regenerated: summary_model/{C}/out/figS9_joint_model_{C}.png"), "summary_model/step8_figure_improved_v2.py")

with open(f"{H}/changes_g1.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["section", "quantity", "text_value", "v2_record", "v3", "v3b", "script", "conclusion_change"]); w.writeheader(); w.writerows(rows)
print(len(rows), "rows;", sum(r["conclusion_change"] != "no" for r in rows), "with a conclusion change")
for r in rows:
    if r["conclusion_change"] != "no" or r["v3"].startswith("n/a") or r["v2_record"].startswith("n/a"):
        print(" ", r["section"], "|", r["quantity"], "|", r["v2_record"], "|", r["v3"], "|", r["v3b"], "|", r["conclusion_change"])
