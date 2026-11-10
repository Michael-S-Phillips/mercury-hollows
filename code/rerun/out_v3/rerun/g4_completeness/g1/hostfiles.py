# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/hostfiles.py (6 Nov 2026). Changed: v3r/v3br catalogs map to the v3/v3b inputs; crater tables read from g1_graph_joint/inputs (see make_copies_g1.py).
"""Input selection for the g1 reruns (6 Nov 2026). CATALOG (env) = v3 or v3b selects the screened host assignments and the
crater table with the recomputed hosting flag; v2 points to the original inputs and is used only to check that the copied
scripts reproduce the v2 records."""
import os
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
CAT = os.environ.get("CATALOG", "v3"); assert CAT in ("v2", "v3", "v3b", "v3r", "v3br"), CAT
BASE = {"v3r": "v3", "v3br": "v3b"}.get(CAT, CAT)   # g4: v3r / v3br use the v3 / v3b catalog inputs
HOSTS = {"v2": f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv",
         "v3": f"{ROOT}/code/rerun/out_v3/host_assignments_v3.csv", "v3b": f"{ROOT}/code/rerun/out_v3/host_assignments_v3b.csv"}[BASE]
CRATERS = f"{ROOT}/code/rerun/out_v3/rerun/g1_graph_joint/inputs/test2_craters_h18_{BASE}.csv"     # v2: verbatim copy of GRL_insolation_hollows/data/test2_craters_h18.csv
