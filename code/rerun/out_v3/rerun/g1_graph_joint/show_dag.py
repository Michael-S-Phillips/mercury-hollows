import json, sys, numpy as np
for C in sys.argv[1:]:
    print("=====", C)
    R = json.load(open(f"{C}/crater_dag_record.json")); T = json.load(open(f"{C}/crater_dag_tests_record.json")); S = json.load(open(f"{C}/crater_dag_sens_record.json"))
    print("counts hosts", R["counts"]["hosts"], "nested hosts", R["counts"]["hosts_nested"], "AIC base/nested/zone/chain/zbase/lnzown/lnzmax/both",
          [round(R[k]["aic"], 1) for k in ["base", "nested", "zone", "chain", "z_base", "ln_zown", "ln_zmax", "ln_zown_zmax"]])
    print(" z_base", R["z_base"]["params"]["z_base"], "lnzmax added", R["ln_zown_zmax"]["params"]["ln_zmax"], "nested coef", R["nested"]["params"]["nested"])
    print(" crust", R["crust_base"]["hosts"], "hollow_hosts", R["hollow_hosts"])
    for k, v in T.items():
        if isinstance(v, dict) and "params" in v:
            ps = {kk: [round(float(np.exp(vv[0])), 2), round(float(np.exp(vv[0]-1.96*vv[1])), 2), round(float(np.exp(vv[0]+1.96*vv[1])), 2), float('%.2g' % vv[2])] for kk, vv in v["params"].items() if kk not in ('Intercept','abslat','C1','nsp','dT','logD','ln_local')}
            print(" ", k, v.get("n", v.get("children")), v["events"], round(v.get("aic", float("nan")), 1), ps, v.get("parents", ""), v.get("rates", ""))
        else: print(" ", k, v)
    for k, v in S.items(): print("  S", k, v)
