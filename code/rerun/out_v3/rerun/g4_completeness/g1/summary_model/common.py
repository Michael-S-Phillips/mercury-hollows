# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/common.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
"""Shared design for the summary model (3 Nov 2026). Reads crater_features.csv (build_features.py, on the Mac).
Groups of variables for the joint model of hosting (h09: >= 1 hollow within 0.9 radii, the paper's definition):
  detection    logit completeness at the crater center
  size         log10 D and its square
  insolation   |lat|, |lat|^2, Delta T_w, Delta T_w x log10 D (centered; the paper's size-dependent hot-longitude term)
  age          fresh, standard, degraded, filled (reference: unrated, i.e. < 10 km or no state), rayed
  nesting      nested, floor and peak position (reference: wall), parent hosts hollows, parent is a basin, ln depth of origin
  cover        on the northern smooth plains, inside Caloris shallower or deeper than the 2.5 km plains base, child excavates only fill
  volcanism    crater contains a vent, parent contains a vent
  composition  LRM fraction within one radius
Pre-specified cross-group combinations (added after the eight groups): parent hosts x vent, nested x Delta T_w, vent x Delta T_w.
Regional clustering: ln of the leave-one-out hosting rate within 500 km, added last."""
import os, numpy as np, pandas as pd
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/common.py; only input/output paths changed. SM_DIR (env) is the
# folder holding this catalog's crater_features*.csv and out/ (g1_graph_joint/summary_model/<CATALOG>).
HERE = os.path.dirname(os.path.abspath(__file__))
SM_DIR = os.environ["SM_DIR"]; _local = os.path.join(SM_DIR, "crater_features.csv")
DATA = os.environ.get("FEATURES", _local if os.path.exists(_local) else
                      "/mnt/user-data/uploads/hollows-emissivity/manuscripts/PlanetaryResearch_hollows/analysis_workshop/summary_model/crater_features.csv")
OUT = os.path.join(SM_DIR, "out"); os.makedirs(OUT, exist_ok=True)

def load():
    d = pd.read_csv(DATA, dtype={"id": str})
    d["lc"] = np.log((d.C1 + 0.01) / (1.01 - d.C1))
    d["logD2"] = (d.logD - d.logD.mean()) ** 2
    d["abslat2"] = (d.abslat - d.abslat.mean()) ** 2 / 100
    d["dT10"] = d.dT / 10
    d["dTxD"] = d.dT10 * (d.logD - d.logD.mean())
    for a in ["fresh", "standard", "degraded", "filled"]: d[a] = (d.age == a).astype(int)
    d["floor"] = (d.zone == "floor").astype(int); d["peak"] = (d.zone == "peak").astype(int)
    d["fill_only"] = d.fill_only.astype(int)
    d["ln_zmax_c"] = d.ln_zmax - d.ln_zmax.mean()
    d["ph_x_vent"] = d.parent_hosts * d.vent
    d["nested_x_dT"] = d.nested * (d.dT10 - d.dT10.mean())
    d["vent_x_dT"] = d.vent * (d.dT10 - d.dT10.mean())
    # spatial blocks ~640 km (15 deg lat bands, longitude cells of about equal area)
    band = np.floor((d.lat + 90) / 15).clip(0, 11).astype(int)
    latc = -90 + 15 * band + 7.5; ncell = np.maximum(1, np.round(24 * np.cos(np.radians(latc)))).astype(int)
    d["block"] = band * 100 + np.floor((d.lon % 360) / 360 * ncell).astype(int)
    return d

GROUPS = {
    "detection": ["lc"],
    "size": ["logD", "logD2"],
    "insolation": ["abslat", "abslat2", "dT10", "dTxD"],
    "age": ["fresh", "standard", "degraded", "filled", "rayed"],
    "nesting": ["nested", "floor", "peak", "parent_hosts", "basin_parent", "ln_zmax_c"],
    "cover": ["nsp", "cal_shallow", "cal_deep", "fill_only"],
    "volcanism": ["vent", "parent_vent"],
    "composition": ["lrm"],
}
INTER = ["ph_x_vent", "nested_x_dT", "vent_x_dT"]
REGION = ["ln_local"]
LABEL = {"lc": "completeness (logit)", "logD": "log D", "logD2": "(log D)^2", "abslat": "|latitude|", "abslat2": "|latitude|^2",
         "dT10": "Delta T_w (per 10 K)", "dTxD": "Delta T_w x log D", "fresh": "fresh", "standard": "standard", "degraded": "degraded",
         "filled": "filled", "rayed": "rayed", "nested": "nested", "floor": "on parent floor", "peak": "on parent peak",
         "parent_hosts": "parent hosts hollows", "basin_parent": "parent is a basin", "ln_zmax_c": "ln depth of origin",
         "nsp": "northern smooth plains", "cal_shallow": "Caloris, above plains base", "cal_deep": "Caloris, below plains base",
         "fill_only": "excavates only fill", "vent": "contains a vent", "parent_vent": "parent contains a vent", "lrm": "LRM fraction",
         "ph_x_vent": "parent hosts x vent", "nested_x_dT": "nested x Delta T_w", "vent_x_dT": "vent x Delta T_w", "ln_local": "regional hosting rate"}

def fit_logit(X, y, w=None, ridge=1e-4, iters=50):
    """Logistic regression by IRLS with a tiny ridge on the slopes (not the intercept) for sparse cells.
    Returns (beta, loglik, cov). X must include an intercept column first."""
    n, p = X.shape; w = np.ones(n) if w is None else w
    b = np.zeros(p); b[0] = np.log(np.average(y, weights=w) / (1 - np.average(y, weights=w)))
    R = ridge * np.eye(p); R[0, 0] = 0
    def pll(bb):
        e = X @ bb; return np.sum(w * (y * e - np.logaddexp(0, e))) - 0.5 * bb @ R @ bb
    cur = pll(b)
    for _ in range(iters):
        eta = X @ b; mu = 0.5 * (1 + np.tanh(eta / 2)); W = w * mu * (1 - mu)
        H = X.T @ (X * W[:, None]) + R; g = X.T @ (w * (y - mu)) - R @ b
        step = np.linalg.lstsq(H, g, rcond=None)[0]; t = 1.0
        while t > 1e-6:                                   # step halving keeps the penalized likelihood rising
            new = pll(b + t * step)
            if new >= cur - 1e-10: break
            t /= 2
        b = b + t * step; done = abs(new - cur) < 1e-9; cur = new
        if done: break
    eta = X @ b; mu = 0.5 * (1 + np.tanh(eta / 2)); W = w * mu * (1 - mu)
    H = X.T @ (X * W[:, None]) + R; ll = np.sum(w * (y * eta - np.logaddexp(0, eta)))
    return b, ll, np.linalg.pinv(H)

def design(d, cols):
    return np.column_stack([np.ones(len(d))] + [d[c].values.astype(float) for c in cols])

def ll_null(y, w=None):
    w = np.ones(len(y)) if w is None else w; p = np.average(y, weights=w)
    return np.sum(w * (y * np.log(p) + (1 - y) * np.log(1 - p)))
