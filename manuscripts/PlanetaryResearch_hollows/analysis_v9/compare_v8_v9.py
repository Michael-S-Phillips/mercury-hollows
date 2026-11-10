#!/usr/bin/env python3
"""Table S17 (model comparison for nesting and depth of origin) from the v8 and v9 crater_dag records, terrestrial and
lunar uplift, plus the k_ring sensitivity of v9. Writes compare_v8_v9_record.json. M. S. Phillips (LPL)."""
import json, os
H = os.path.dirname(os.path.abspath(__file__)); P = os.path.dirname(H)
runs = {"v8_terrestrial": f"{P}/analysis_v8/crater_dag_record.json", "v8_lunar": f"{P}/analysis_v8/uplift_lunar/crater_dag_record.json",
        "v9_terrestrial": f"{H}/crater_dag_record.json", "v9_lunar": f"{H}/uplift_lunar/crater_dag_record.json",
        "v9_terrestrial_k0.8": f"{H}/kring_0p8/crater_dag_record.json", "v9_terrestrial_k1.3": f"{H}/kring_1p3/crater_dag_record.json"}
def term(R, model, name):
    v = R[model]["params"][name]; return dict(coef=v[0], se=v[1], p=v[2])
out = {}
for k, f in runs.items():
    R = json.load(open(f)); base = R["base"]["aic"]
    out[k] = dict(
        base_aic=base, nested_aic=R["nested"]["aic"], zone_aic=R["zone"]["aic"], chain_aic=R["chain"]["aic"],
        z_base_aic=R["z_base"]["aic"], z_base=term(R, "z_base", "z_base"),
        ln_zown_aic=R["ln_zown"]["aic"], ln_zown=term(R, "ln_zown", "ln_zown"),
        ln_zmax_aic=R["ln_zmax"]["aic"], ln_zmax=term(R, "ln_zmax", "ln_zmax"),
        ln_zown_zmax_aic=R["ln_zown_zmax"]["aic"], ln_zmax_added=term(R, "ln_zown_zmax", "ln_zmax"),
        dAIC_ln_zmax_added_to_ln_zown=R["ln_zown"]["aic"] - R["ln_zown_zmax"]["aic"],
        dAIC_best_depth_model_minus_nested=min(R["z_base"]["aic"], R["ln_zown"]["aic"], R["ln_zmax"]["aic"], R["ln_zown_zmax"]["aic"]) - R["nested"]["aic"],
        fillknown=dict(z_base_aic=R["z_base_fillknown"]["aic"], z_base=term(R, "z_base_fillknown", "z_base"),
                       ln_zown_aic=R["ln_zown_fillknown"]["aic"], ln_zown=term(R, "ln_zown_fillknown", "ln_zown"),
                       ln_zmax_aic=R["ln_zmax_fillknown"]["aic"], ln_zmax=term(R, "ln_zmax_fillknown", "ln_zmax")),
        fill_only=R["counts"]["fill_only"], hollow_hosts=R["hollow_hosts"],
        caloris=dict(SU=R["settings"].get("caloris_SU_km"), z_own=R["settings"].get("caloris_z_own_km")),
        crust_base=R.get("crust_base"))
json.dump(out, open(f"{H}/compare_v8_v9_record.json", "w"), indent=1)
for k, v in out.items():
    print(f"{k:22s} base {v['base_aic']:.1f} nested {v['nested_aic']:.1f} | z_base {v['z_base_aic']:.1f} ({v['z_base']['coef']:.4f}±{v['z_base']['se']:.4f}, p={v['z_base']['p']:.2g})"
          f" | ln_zown {v['ln_zown_aic']:.1f} | ln_zmax {v['ln_zmax_aic']:.1f} | +ln_zmax {v['ln_zown_zmax_aic']:.1f} (dAIC {v['dAIC_ln_zmax_added_to_ln_zown']:.1f};"
          f" b {v['ln_zmax_added']['coef']:.3f}±{v['ln_zmax_added']['se']:.3f}, p={v['ln_zmax_added']['p']:.2g}) | best-depth minus nested {v['dAIC_best_depth_model_minus_nested']:.1f} | fill_only {v['fill_only']}")
    fk = v["fillknown"]; print(f"{'':22s} fill-known: z_base {fk['z_base_aic']:.1f} ({fk['z_base']['coef']:.4f}±{fk['z_base']['se']:.4f}, p={fk['z_base']['p']:.2g}) ln_zown {fk['ln_zown_aic']:.1f} ({fk['ln_zown']['coef']:.2f}±{fk['ln_zown']['se']:.2f}) ln_zmax {fk['ln_zmax_aic']:.1f} ({fk['ln_zmax']['coef']:.2f}±{fk['ln_zmax']['se']:.2f})")
    if v["crust_base"]: print(f"{'':22s} crust base: {v['crust_base']['hosts']} ; all {v['crust_base']['all_craters']}; Caloris SU {v['caloris']['SU']}, z_own {v['caloris']['z_own']}")
