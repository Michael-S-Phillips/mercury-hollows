"""Input selection for the g1 reruns (6 Nov 2026). CATALOG (env) = v3 or v3b selects the screened host assignments and the
crater table with the recomputed hosting flag; v2 points to the original inputs and is used only to check that the copied
scripts reproduce the v2 records."""
import os
ROOT = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__))
CAT = os.environ.get("CATALOG", "v3"); assert CAT in ("v2", "v3", "v3b"), CAT
HOSTS = {"v2": f"{ROOT}/manuscripts/GRL_insolation_hollows/data/host_assignments.csv",
         "v3": f"{ROOT}/code/rerun/out_v3/host_assignments_v3.csv", "v3b": f"{ROOT}/code/rerun/out_v3/host_assignments_v3b.csv"}[CAT]
CRATERS = f"{HERE}/inputs/test2_craters_h18_{CAT}.csv"     # v2: verbatim copy of GRL_insolation_hollows/data/test2_craters_h18.csv
