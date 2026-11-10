#!/usr/bin/env python3
"""Copies the size-distribution scripts into sfd/ and changes only input paths, output paths and env settings.
Env: G3_CAT = v2 | v3 | v3b (catalog); G3_THETA = fixed (out_v2 detection model, theta_out_v2.json) | refit (detection model
refit on the v3/v3b pairs by sfd/refit_detection.py). Outputs go to sfd/out_{G3_CAT}_{G3_THETA}/."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); SRC = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd"
DST = os.path.dirname(os.path.abspath(__file__)) + "/sfd"
HDR = ('import os as _os\nG3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")\n'
       'CAT = _os.environ.get("G3_CAT", "v3"); THETA = _os.environ.get("G3_THETA", "fixed")\n'
       'ODIR = f"{G3}/sfd/out_{CAT}_{THETA}"; _os.makedirs(ODIR, exist_ok=True)\n'
       'THETA_JSON = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json")\n'
       '              if THETA == "fixed" else f"{G3}/sfd/out_{CAT}_refit/theta_out_v2.json")\n')
def copy(name, reps, hdr_after=None):
    t = open(f"{SRC}/{name}").read()
    for a, b in reps:
        assert t.count(a) == 1, (name, a, t.count(a)); t = t.replace(a, b)
    if hdr_after:
        assert t.count(hdr_after) == 1; t = t.replace(hdr_after, hdr_after + HDR, 1)
    else:
        t = HDR + t
    t = f"# g3_catalog copy of analysis_workshop/completeness_sfd/{name} (6 Nov 2026): only input/output paths and env settings changed.\n" + t
    open(f"{DST}/{name}", "w").write(t)
copy("refit_detection.py", [('pd.read_csv(f"{HE}/code/rerun/out_v2/completeness_pairs.csv")', 'pd.read_csv(f"{G3}/inputs/completeness_pairs_{CAT}.csv")'),
                            ('open(f"{HERE}/theta_out_v2.json", "w")', 'open(f"{ODIR}/theta_out_v2.json", "w")')],
     hdr_after='from scipy.optimize import minimize\n')
copy("fit_sfd.py", [('json.load(open(f"{HERE}/theta_out_v2.json"))', 'json.load(open(THETA_JSON))'),
                    ('pd.read_csv(f"{HE}/code/rerun/out_v2/completeness_pairs.csv")', 'pd.read_csv(f"{G3}/inputs/completeness_pairs_{CAT}.csv")'),
                    ('open(f"{HERE}/sfd_fit.json", "w")', 'open(f"{ODIR}/sfd_fit.json", "w")'),
                    ('np.savez_compressed(f"{HERE}/sfd_grid.npz"', 'np.savez_compressed(f"{ODIR}/sfd_grid.npz"')],
     hdr_after='from scipy.optimize import minimize\n')
copy("sfd_windows.py", [('U = "/mnt/user-data/uploads/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/"', 'U = ODIR + "/"'),
                        ('open("sfd_windows.json", "w")', 'open(f"{ODIR}/sfd_windows.json", "w")')])
copy("sfd_numbers.py", [('U = "/mnt/user-data/uploads/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/"', 'U = ODIR + "/"'),
                        ('open("sfd_numbers.json", "w")', 'open(f"{ODIR}/sfd_numbers.json", "w")')])
print("copied")
