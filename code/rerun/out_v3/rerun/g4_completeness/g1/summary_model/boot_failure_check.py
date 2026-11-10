# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/boot_failure_check.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""g1 rerun (6 Nov 2026), new diagnostic: the complementary log-log fits of step7_improved.py / step9_color.py (IRLS with
step halving, linear predictor clipped at 5) fail to converge in some bootstrap draws, which then show negative general
dominance (adding a group lowers R^2). This counts, per catalog and analysis, the cached bootstrap draws with any group's
general dominance < -0.002 or total R^2 < 0.1, and gives the 95% share intervals with and without those draws (the paper's
intervals are the 'all draws' ones, as the copied scripts compute them). Writes <CATALOG>/out/boot_failure_check.json."""
import os, sys, json, glob, pickle, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); R = {}
for C in sys.argv[1:]:
    R[C] = {}
    for script, pick in [("step7_improved", lambda r: r), ("step9_color", lambda r: r), ("step3_dominance", lambda r: (None, r))]:
        fs = glob.glob(f"{H}/{C}/cache/{script}/*.pkl")
        if not fs: continue
        res = [pick(pickle.load(open(f, "rb"))) for f in fs]
        for name in sorted(set(n for n, _ in res), key=str):
            bs = [r for n, r in res if n == name]; full = lambda b: b["full"]
            bad = [b for b in bs if min(b["general"].values()) < -0.002 or b["full"] < 0.1]; good = [b for b in bs if b not in bad]
            q = lambda xs: [round(float(np.percentile(xs, 2.5)) * 100, 1), round(float(np.percentile(xs, 97.5)) * 100, 1)] if len(xs) > 10 else None
            R[C][f"{script}:{name}"] = dict(draws=len(bs), failed=len(bad),
                ci_share_all={g: q([b["general"][g] / b["full"] for b in bs]) for g in bs[0]["general"]},
                ci_share_without_failed={g: q([b["general"][g] / b["full"] for b in good]) for g in bs[0]["general"]},
                ci_R2_all=q([b["full"] for b in bs]), ci_R2_without_failed=q([b["full"] for b in good]))
            print(C, script, name, "failed", len(bad), "of", len(bs), "| without failed:", R[C][f"{script}:{name}"]["ci_share_without_failed"])
    json.dump(R[C], open(f"{H}/{C}/out/boot_failure_check.json", "w"), indent=1)
