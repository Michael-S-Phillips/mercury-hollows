#!/usr/bin/env python3
"""Before/after table for the v8 Caloris merge (and the v8 variants). Reads the record JSONs of
v6_rerun_fullprecision/ (v6 scripts rerun unchanged except unrounded coefficients in the tests record; identical to the
analysis_v6 records to their rounding), analysis_v8/ (merged Caloris, catalog record), alt_caloris_orgel/ (merged
Caloris, Orgel record) and uplift_lunar/ (merged Caloris, lunar uplift; crater_dag.py only). Writes compare_v6_v8_record.json.
ORs are exp(b) with Wald 95% intervals exp(b +- 1.96 se)."""
import os, json, math
H = os.path.dirname(os.path.abspath(__file__))
runs = {"v6": f"{H}/v6_rerun_fullprecision", "v8": H, "v8_orgel": f"{H}/alt_caloris_orgel", "v8_lunarSU": f"{H}/uplift_lunar"}
def load(d, f):
    p = f"{d}/{f}"; return json.load(open(p)) if os.path.exists(p) else None
def orci(b, se): return [round(math.exp(b), 2), round(math.exp(b - 1.96 * se), 2), round(math.exp(b + 1.96 * se), 2)]
OUT = {}
for name, d in runs.items():
    A, T, S = load(d, "crater_dag_record.json"), load(d, "crater_dag_tests_record.json"), load(d, "crater_dag_sens_record.json")
    o = {}
    c = A["counts"]
    o["counts"] = dict(nodes=c["nodes"], basins_kept=c["basins_kept"], nested=c["nested"], **{f"zone_{k}": v for k, v in c["by_zone"].items()},
                       **{f"chain_{k}": v for k, v in c["chain"].items()}, fill_only=c["fill_only"],
                       fill_unknown_nested=c["fill_unknown_nested"], hosts=c["hosts"], hosts_nested=c["hosts_nested"])
    fits = {}
    for k in ["base", "nested", "zone", "chain", "z_base", "ln_zown", "ln_zmax", "ln_zown_zmax", "z_base_fillknown"]:
        f = A[k]; fits[k] = dict(n=f["n"], hosts=f["hosts"], aic=round(f["aic"], 1),
                                 coef={p: [round(v[0], 4), round(v[1], 4), float(f"{v[2]:.3g}")] for p, v in f["params"].items()
                                       if p not in ("Intercept", "logD", "abslat", "C1", "nsp", "dT")})
    fits["dAIC_vs_base"] = {k: round(A[k]["aic"] - A["base"]["aic"], 1) for k in ["nested", "zone", "chain", "z_base", "ln_zown", "ln_zmax", "ln_zown_zmax"]}
    fits["dAIC_vs_nested"] = {k: round(A[k]["aic"] - A["nested"]["aic"], 1) for k in ["zone", "chain", "z_base", "ln_zown", "ln_zmax", "ln_zown_zmax"]}
    o["fits_hosts"] = fits
    if T is not None:
        h = {}
        for k, v in S.items():
            if k == "matched_300km": continue
            if isinstance(v, str): h[k] = v; continue
            kk = [x for x in v if x not in ("n", "events")][0]
            h[k] = dict(n=v["n"], events=v["events"], OR=v[kk]["OR"], ci=v[kk]["ci"], p=float(f"{v[kk]['p']:.3g}"))
        def tk(rec, terms):
            return {t.split("[T.")[1].rstrip("]"): dict(zip(["OR", "lo", "hi"], orci(*rec["params"][t][:2])), p=float(f"{rec['params'][t][2]:.3g}")) for t in terms}
        pk = T["parentkind_local_h09"]; h["parent_kind"] = tk(pk, [t for t in pk["params"] if t.startswith("C(pk)")])
        ch = T["chain_local_h09"]; h["chain"] = tk(ch, [t for t in ch["params"] if t.startswith("C(chain")])
        zn = T["zone_local_h09"]; h["zone"] = tk(zn, [t for t in zn["params"] if t.startswith("C(zone_c)")])
        w = T["within_parent_clogit_h09"]; b, se, p = w["params"]["inner"]
        h["within_parent"] = dict(zip(["OR", "lo", "hi"], orci(b, se)), p=round(p, 3), parents=w["parents"], children=w["children"],
                                  events=w["events"], rates=w["rates"])
        ih = T["inherit_h09"]; b, se, p = ih["params"]["parent_hosts"]
        h["inherit_anyhollow_parent"] = dict(zip(["OR", "lo", "hi"], orci(b, se)), p=float(f"{p:.3g}"), n=ih["n"], events=ih["events"],
                                             rates=T["inherit_rates_h09"])
        h["nested_local_h09_tests"] = dict(zip(["OR", "lo", "hi"], orci(*T["nested_local_h09"]["params"]["nested"][:2])))
        o["h09"] = h
        o["matched_300km"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in S["matched_300km"].items()}
    OUT[name] = o
json.dump(OUT, open(f"{H}/compare_v6_v8_record.json", "w"), indent=1)
print(json.dumps(OUT, separators=(",", ":")))
