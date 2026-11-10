# g3_catalog copy of code/rerun/closure/size_law.py (6 Nov 2026): only input/output paths and env settings changed (env G3_CAT).
"""Closure track: completeness-corrected size law, v1 reproduction + v2 re-fit.
Reproduces size_intrinsic_law.json (northern non-crater-trace hollows) on v1
inputs, then re-fits on data_v2. Also computes area-weighted catalog
completeness at D=0.5/1/2/5 km, north and south, with cell-bootstrap CIs.
"""
import json, os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from collections import defaultdict

BASE = os.path.expanduser("~/mnt/hollows-emissivity/hollow-networks/")
CAT = os.environ.get("G3_CAT", "v3")
CANON = {"v2": os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv"),
         "v3": os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_v3.csv"),
         "v3b": os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_v3b.csv")}[CAT]
OUT = os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog/sfd/")

Dlo, Dhi = 0.1, 10.0
G = np.geomspace(Dlo, Dhi, 300)
logG = np.log(G)
dD = np.gradient(G)


def cell_of(lat, lon):
    row = np.clip(((90.0 - lat) * 10).astype(int), 0, 1799)
    col = np.clip((np.mod(lon, 360) * 10).astype(int), 0, 3599)
    return row, col


def f_power(D, alpha):
    if abs(alpha - 1) < 1e-6:
        Z = np.log(Dhi / Dlo)
    else:
        Z = (Dhi ** (1 - alpha) - Dlo ** (1 - alpha)) / (1 - alpha)
    return D ** (-alpha) / Z


def f_expo(D, D0):
    Z = D0 * (np.exp(-Dlo / D0) - np.exp(-Dhi / D0))
    return np.exp(-D / D0) / Z


def f_logn(D, mu, sig):
    Z = norm.cdf((np.log(Dhi) - mu) / sig) - norm.cdf((np.log(Dlo) - mu) / sig)
    return np.exp(-0.5 * ((np.log(D) - mu) / sig) ** 2) / (D * sig * np.sqrt(2 * np.pi)) / Z


def raw_alpha(D, dlo, dhi):
    x = D[(D >= dlo) & (D <= dhi)]
    def nll(a):
        a = a[0]
        if abs(a - 1) < 1e-6:
            Z = np.log(dhi / dlo)
        else:
            Z = (dhi ** (1 - a) - dlo ** (1 - a)) / (1 - a)
        return -np.sum(-a * np.log(x) - np.log(Z))
    return float(minimize(nll, [2.0], method='Nelder-Mead').x[0]), int(len(x))


def fit_size(datadir, B=200):
    df = pd.read_csv(CANON)
    n = df[(df.lat > 0) & (~df.is_crater_trace.astype(bool))].copy()
    comp = np.load(datadir + "completeness_10ppd.npz")
    Dref = comp['Dref']; C = comp['C']; logDref = np.log(Dref)

    r, c = cell_of(n.lat.values, n.lon.values)
    Cprof = C[:, r, c].T
    N = len(n)
    Cgrid = np.empty((N, len(G)))
    for i in range(N):
        Cgrid[i] = np.interp(logG, logDref, Cprof[i])
    Dobs = n.D_eq_km.values
    Cown = n.C_own.values

    mask = (Dobs >= Dlo) & (Dobs <= Dhi) & (Cgrid.sum(1) > 0)
    n_f = int(mask.sum())
    Dobs_f = Dobs[mask]; Cown_f = Cown[mask]; Cgrid_f = Cgrid[mask]

    def Zsite(fg, Cg):
        return Cg @ (fg * dD)

    def negll(params, model, Dob, Cg, Cow):
        if model == 'power':
            a = params[0]; fobs = f_power(Dob, a); fg = f_power(G, a)
        elif model == 'expo':
            D0 = params[0]
            if D0 <= 0:
                return 1e12
            fobs = f_expo(Dob, D0); fg = f_expo(G, D0)
        elif model == 'logn':
            mu, sig = params
            if sig <= 0:
                return 1e12
            fobs = f_logn(Dob, mu, sig); fg = f_logn(G, mu, sig)
        Zi = Zsite(fg, Cg)
        if np.any(Zi <= 0) or np.any(fobs <= 0):
            return 1e12
        ll = np.log(fobs) + np.log(Cow) - np.log(Zi)
        if not np.all(np.isfinite(ll)):
            return 1e12
        return -np.sum(ll)

    rp = minimize(lambda p: negll(p, 'power', Dobs_f, Cgrid_f, Cown_f), [3.0],
                  method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-6, maxiter=3000))
    re_ = minimize(lambda p: negll(p, 'expo', Dobs_f, Cgrid_f, Cown_f), [1.0],
                   method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-6, maxiter=3000))
    rl = minimize(lambda p: negll(p, 'logn', Dobs_f, Cgrid_f, Cown_f), [0.5, 0.8],
                  method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-6, maxiter=3000))

    # Field-level block bootstrap (same construction/seed as v1 lineage)
    nf_df = n[mask].reset_index(drop=True)
    fid_f = nf_df.field_id.values
    fields = np.unique(fid_f)
    fidx = defaultdict(list)
    for i, f in enumerate(fid_f):
        fidx[f].append(i)
    fidx = {k: np.array(v) for k, v in fidx.items()}

    rng = np.random.default_rng(1)
    al_bs = []; mu_bs = []; sg_bs = []
    for b in range(B):
        samp = rng.choice(fields, size=len(fields), replace=True)
        idx = np.concatenate([fidx[f] for f in samp])
        Do = Dobs_f[idx]; Co = Cown_f[idx]; Cg = Cgrid_f[idx]
        def nllp(a):
            a = a[0]
            if abs(a - 1) < 1e-6:
                Zn = np.log(Dhi / Dlo)
            else:
                Zn = (Dhi ** (1 - a) - Dlo ** (1 - a)) / (1 - a)
            fg = G ** (-a) / Zn
            Zi = Cg @ (fg * dD)
            if np.any(Zi <= 0):
                return 1e12
            return -np.sum(-a * np.log(Do) - np.log(Zn) + np.log(Co) - np.log(Zi))
        al_bs.append(minimize(nllp, [2.7], method='Nelder-Mead').x[0])
        def nlll(p):
            mu, s = p
            if s <= 0:
                return 1e12
            fo = f_logn(Do, mu, s); fg = f_logn(G, mu, s)
            Zi = Cg @ (fg * dD)
            if np.any(Zi <= 0) or np.any(fo <= 0):
                return 1e12
            return -np.sum(np.log(fo) + np.log(Co) - np.log(Zi))
        rr = minimize(nlll, [-3.9, 1.5], method='Nelder-Mead')
        mu_bs.append(rr.x[0]); sg_bs.append(rr.x[1])
    al_bs = np.array(al_bs); mu_bs = np.array(mu_bs); sg_bs = np.array(sg_bs)

    # support-floor sensitivity (break diagnostic)
    floor = {}
    for dlo in [0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
        m2 = (Dobs_f >= dlo)
        Do = Dobs_f[m2]; Co = Cown_f[m2]; Cg = Cgrid_f[m2]
        def nllp(a):
            a = a[0]
            if abs(a - 1) < 1e-6:
                Zn = np.log(Dhi / dlo)
            else:
                Zn = (Dhi ** (1 - a) - dlo ** (1 - a)) / (1 - a)
            fg = np.where(G >= dlo, G ** (-a) / Zn, 0.0)
            Zi = Cg @ (fg * dD)
            if np.any(Zi <= 0):
                return 1e12
            return -np.sum(-a * np.log(Do) - np.log(Zn) + np.log(Co) - np.log(Zi))
        floor[str(dlo)] = float(minimize(nllp, [2.7], method='Nelder-Mead').x[0])

    ra_05_5 = raw_alpha(Dobs_f, 0.5, 5.0)
    ra_03_5 = raw_alpha(Dobs_f, 0.3, 5.0)

    return dict(
        N=int(N), n_fitted=n_f,
        alpha_point=float(rp.x[0]),
        alpha_bootstrap_median=float(np.median(al_bs)),
        alpha_95CI=[float(np.percentile(al_bs, 2.5)), float(np.percentile(al_bs, 97.5))],
        cumulative_slope=float(np.median(al_bs) - 1),
        nll_power=float(rp.fun),
        logn_mu=float(rl.x[0]), logn_sigma=float(rl.x[1]), nll_logn=float(rl.fun),
        expo_D0=float(re_.x[0]), nll_expo=float(re_.fun),
        alpha_by_Dlo=floor,
        raw_diff_0p5_5=ra_05_5[0], raw_diff_0p5_5_n=ra_05_5[1],
        raw_diff_0p3_5=ra_03_5[0], raw_diff_0p3_5_n=ra_03_5[1],
    )


def comp_at_size(datadir, sizes=(0.5, 1.0, 2.0, 5.0), B=300, seed=0):
    """Area-weighted catalog completeness at fixed D, north and south,
    reproducing the completeness_model formula, with cell-bootstrap CIs."""
    comp = np.load(datadir + "completeness_10ppd.npz")
    cov = np.load(datadir + "coverage_10ppd.npz")
    C = comp['C']; Dref = list(np.round(comp['Dref'], 6))
    nimg = cov['nimg']
    NY, NX = C.shape[1], C.shape[2]
    latc = 90 - (np.arange(NY) + 0.5) / 10.0
    W = np.cos(np.radians(latc))[:, None] * np.ones((1, NX))
    covm = nimg > 0
    latg = latc[:, None] * np.ones((1, NX))
    out = {}
    rng = np.random.default_rng(seed)
    for hemi, sel in [("north", (latg > 0) & covm), ("south", (latg < 0) & covm)]:
        wsel = W[sel]
        res = {}
        for D in sizes:
            di = Dref.index(round(float(D), 6))
            csel = C[di][sel]
            point = float(np.sum(wsel * csel) / np.sum(wsel))
            # cell bootstrap
            m = len(csel)
            bs = np.empty(B)
            for b in range(B):
                idx = rng.integers(0, m, m)
                bs[b] = np.sum(wsel[idx] * csel[idx]) / np.sum(wsel[idx])
            res[str(D)] = dict(point=point,
                               ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])
        out[hemi] = res
    return out


if __name__ == "__main__":
    rec = {}
    # g3: one catalog per run (env G3_CAT), on the unchanged v2 completeness grid (data_v2 = out_v2, identical bytes).
    # The v1 fit and comp_at_size (which does not use the catalog) are not rerun.
    print("=== %s fit ===" % CAT)
    rec['size_' + CAT] = fit_size(BASE + "data_v2/")
    print(json.dumps(rec['size_' + CAT], indent=1))
    json.dump(rec, open(OUT + "size_law_record_%s.json" % CAT, "w"), indent=1)
    print("wrote size_law_record.json")
