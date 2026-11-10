# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/cached_pool.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Runs one of the copied summary-model scripts UNCHANGED within the ~180 s limit of each shell call on this machine
(background jobs do not survive a call). Usage: SM_DIR=... CATALOG=... python3 cached_pool.py <script.py> [script args]
multiprocessing.Pool is replaced by a pool whose map() returns results cached on disk (SM_DIR/cache/<script>/<job>.pkl) and
computes the missing jobs with a real fork pool of NPROC (env, default 4) workers, saving each result as it arrives. When the
time budget (BUDGET s, default 140) is spent, it exits with code 3; calling it again resumes. Each bootstrap job carries its
own seed, so the results do not depend on how the jobs are split across calls. When every job is cached, the script runs to
the end exactly as it would with its own Pool."""
import os, sys, time, pickle, runpy, multiprocessing as mp
T0 = time.time(); BUDGET = float(os.environ.get("BUDGET", "140")); NPROC = int(os.environ.get("NPROC", "4"))
SCRIPT = sys.argv[1]; CACHE = os.path.join(os.environ["SM_DIR"], "cache", os.path.basename(SCRIPT)[:-3]); os.makedirs(CACHE, exist_ok=True)
_RealPool = mp.get_context("fork").Pool; _FUNC = [None]
def _key(job): return "".join(c if c.isalnum() or c in "._-" else "_" for c in repr(job))
import types; _H = types.ModuleType("_cached_pool_helper"); sys.modules["_cached_pool_helper"] = _H   # picklable by reference under runpy
exec("def _call(job):\n    return job, _FUNC[0](job)\n", _H.__dict__); _H._FUNC = _FUNC; _call = _H._call
class CachedPool:
    def __init__(self, *a, **k): pass
    def __enter__(self): return self
    def __exit__(self, *a): return False
    def map(self, func, jobs, chunksize=1):
        jobs = list(jobs); todo = [j for j in jobs if not os.path.exists(os.path.join(CACHE, _key(j) + ".pkl"))]
        print(f"[cached_pool] {SCRIPT}: {len(jobs)} jobs, {len(todo)} to compute", file=sys.stderr, flush=True)
        if todo:
            _FUNC[0] = func; p = _RealPool(NPROC)
            try:
                for job, res in p.imap_unordered(_call, todo):
                    tmp = os.path.join(CACHE, _key(job) + ".tmp"); pickle.dump(res, open(tmp, "wb")); os.replace(tmp, tmp[:-4] + ".pkl")
                    if time.time() - T0 > BUDGET:
                        p.terminate(); n = sum(f.endswith(".pkl") for f in os.listdir(CACHE))
                        print(f"[cached_pool] budget spent; {n} results cached; call again", file=sys.stderr, flush=True); os._exit(3)
            finally:
                p.terminate()
        return [pickle.load(open(os.path.join(CACHE, _key(j) + ".pkl"), "rb")) for j in jobs]
mp.Pool = CachedPool
if len(sys.argv) > 3 and sys.argv[2] == "--fill":
    # fill mode: cached_pool.py <script.py> --fill <function> "<python list of jobs>": loads the script without running its
    # main block and computes the missing jobs of its own function into the same cache, so that the later full run only
    # repeats the serial part once.
    FN, JOBS = sys.argv[3], eval(sys.argv[4]); sys.argv = [SCRIPT]
    sys.path.insert(0, os.path.dirname(os.path.abspath(SCRIPT))); g = runpy.run_path(SCRIPT, run_name="g1_fill")
    CachedPool().map(g[FN], JOBS); print("[cached_pool] fill complete", file=sys.stderr)
else:
    sys.argv = sys.argv[1:]; sys.path.insert(0, os.path.dirname(os.path.abspath(SCRIPT)))
    runpy.run_path(SCRIPT, run_name="__main__")
