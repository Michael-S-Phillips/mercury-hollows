# g4 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/joint_thinned.py (6 Nov 2026). Changed: input/output paths by S2_CAT, added exponents/measures 1.7 and 1.8 (see make_copies_s2.py).
#!/usr/bin/env python3
"""Joint model of hosting (summary_model/step7_improved.py, thinned occurrence, complementary log-log with offset ln A) with each
completeness measure as the detection term in place of ln(C1 + 0.01). Point estimates: g (coefficient on the detection term; 1 =
pure thinning), log-likelihood of the full model, R^2 of the physical groups, shares, and the regional increment.
Usage: python3 joint_thinned.py <measure> ... ; appends to joint_thinned.json."""
import os, sys, json, numpy as np, pandas as pd
import os as _o, sys as _s; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__))); from s2cfg import *   # g4: inputs by S2_CAT
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SM_CODE); sys.argv_saved = sys.argv; sys.argv = [sys.argv[0]]
import step7_improved as s7
m = pd.read_csv(CCOMP, dtype={"id": str})
base = s7.d.copy(); out_p = f"{HERE}/out/joint_thinned_{CAT}.json"; R = json.load(open(out_p)) if os.path.exists(out_p) else {}
for M in sys.argv_saved[1:]:
    d = base.merge(m[["id", "C1_act", "E_1.3", "E_1.6", "E_1.7", "E_1.8", "E_1.9", "F_0.3"] + [c for c in m.columns if c.startswith("C_")]], on="id", how="left")
    if M.startswith("lnC_"): d["lnC"] = np.log(d[M[2:]] + 0.01)
    elif M.startswith("lnCold_"):          # the paper's completeness (60 deg incidence, band pixel scales) at another reference size
        Z = np.load(GRID); k = list(np.round(Z["Dref"], 3)).index(float(M[7:]))
        rr = ((90 - d.lat.values) * 10).astype(int).clip(0, 1799); cc = ((d.lon.values % 360) * 10).astype(int).clip(0, 3599)
        d["lnC"] = np.log(np.nan_to_num(Z["C"][k][rr, cc]) + 0.01)
    else: d["lnC"] = {"lnC1": np.log(d.C1 + 0.01), "lnC1_act": np.log(d.C1_act + 0.01), "lnE_1.6": np.log(d["E_1.6"] + 0.01),
                "lnE_1.3": np.log(d["E_1.3"] + 0.01), "lnE_1.9": np.log(d["E_1.9"] + 0.01), "lnE_1.7": np.log(d["E_1.7"] + 0.01), "lnE_1.8": np.log(d["E_1.8"] + 0.01), "lnF_0.3": np.log(d["F_0.3"] + 0.01)}[M]
    s7.d = d; get, groups, fam = s7.ANALYSES["thinned_presence_all"]; dd = get()
    y, off, fit, r2_all, dominance = s7.make(dd, groups, fam)
    cols = [c for g in groups for c in groups[g]]; b, llf = fit(cols, np.ones(len(dd))); b0, ll_det = fit([], np.ones(len(dd)))
    P = dominance(r2_all(np.ones(len(dd))))
    R[M] = {"g": float(b[1]), "llf_full": float(llf), "llf_detection_and_area_only": float(ll_det), "R2_physical": P["full"],
            "share": {k: round(v, 4) for k, v in P["share"].items()}, "region_add": P["region_add"]}
    print(M, json.dumps(R[M]), flush=True)
    json.dump(R, open(out_p, "w"), indent=1)
