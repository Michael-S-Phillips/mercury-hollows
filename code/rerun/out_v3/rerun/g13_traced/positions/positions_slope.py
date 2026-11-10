# g13 revision (7 Nov 2026): with G11_RULES=traced and G13_REV=r2 (default) outputs go to out_traced_r2; G13_REV=r1 keeps out_traced.
# g13 copy of g11_curvature/positions/positions_slope.py (6 Nov 2026). Changed: folder g13_traced/positions; G11_RULES default traced (traced outlines where available, curvature v3 elsewhere via g13_traced/traced_classes.py; env CLS=v2 for curvature v2 elsewhere).
# g11 copy of g9_draft_reruns/positions/positions_slope.py (6 Nov 2026). Changed: env G11_RULES / G11_ETA (eta | eta_b: eta_b replaces eta for cells, hollows and fields) / G11_SCREENS (g11: main >= 30 peak cells with no relief screen, sensitivity >= 30 cells with relief >= 500 m and >= 750 m; g9: the g9 screens), the folder it reads and writes (out_<RULES>/), output names positions_slope_record_<ETA>[_g9screens].json and fig_positions_slope_<ETA>[_g9screens].png, the main screen key, and the eta-slope correlation computed for every screen in g11 mode.
# g9 copy of analysis_workshop/melt_source_tests/positions/positions_slope.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID and the folder it reads and writes (positions/out_<CAT>_<RULES>_grid<GRID>).
#!/usr/bin/env python3
"""Slope follow-up, step 2 (5 Nov 2026): does the excess of hollows high on peaks survive a control for DEM slope?
Reads cache/slope_H250.npz (positions_slope_build.py). Resolution screen on peak terrain: main >= 30 valid peak cells and
z_peak95 - z_floor >= 750 m; loose >= 15 / 500 m; strict >= 60 / 1,000 m.
Height bands: eta thirds (low < 1/3, middle 1/3-2/3, upper >= 2/3). Slope strata: within-crater completeness-weighted slope
terciles, fixed bins 0-1, 1-2, 2-4, 4-8, > 8 deg, and the bins 0-5, 5-10, 10-15, > 15 deg.
(a) Mantel-Haenszel rate ratio, upper third vs lower two thirds (and lowest third vs upper two thirds), strata = crater (no slope
    control) or crater x slope stratum; exposure = completeness-weighted area; 95% intervals from 2,000 crater bootstrap resamples.
(b) Poisson regression with crater fixed effects, fitted as the equivalent conditional Poisson on counts aggregated by
    crater x eta third x slope bin with log weighted area as offset; eta reference = middle third, slope reference = lowest bin.
    Intervals: cluster-robust (crater) Wald and 2,000 crater bootstrap resamples. One check fit with statsmodels GLM + dummies.
Hollows and hollow fields (one count per field at its centroid). Slope at 1-cell and 2-cell baselines.
Writes positions_slope_record.json and fig_positions_slope.png."""
import os, json, numpy as np, pandas as pd
import os as _os
import sys as _sys; _sys.dont_write_bytecode = True
G11 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced")
CAT = "v3"; RULES = _os.environ.get("G11_RULES", "traced"); GRID = "v3"; assert RULES in ("curv", "fix", "traced"), RULES
ETA = _os.environ.get("G11_ETA", "eta"); SCR = _os.environ.get("G11_SCREENS", "g11"); assert ETA in ("eta", "eta_b") and SCR in ("g11", "g9")
SUF = ETA + ("_g9screens" if SCR == "g9" else "")
HERE = f"{G11}/positions/out_{RULES}" + ("_r2" if RULES == "traced" and _os.environ.get("G13_REV", "r2") == "r2" else ""); NB = 2000; INF = np.inf
d = np.load(f"{HERE}/cache/slope_H250.npz")
cr = pd.DataFrame(d["craters"]); pk = pd.DataFrame(d["peak"], columns=["uid", "a", "eta", "s1", "s2", "h"])
ho = pd.DataFrame(d["hollows"], columns=["uid", "hid", "fid", "eta", "s1", "s2"]); fd = pd.DataFrame(d["fields"], columns=["uid", "fid", "n", "fb", "eta", "s1", "s2"])
if ETA == "eta_b": pk["eta"] = d["peak_etab"]; ho["eta"] = d["hollows_etab"]; fd["eta"] = d["fields_etab"]
cr["zbase_rel"] = d["zbase_rel"] if "zbase_rel" in d.files else np.nan
SCREENS = ({"loose_15cells_500m": (15, 500), "main_30cells_750m": (30, 750), "strict_60cells_1000m": (60, 1000), "none": (1, -INF)} if SCR == "g9" else
           {"main_30cells_norelief": (30, -INF), "sens_30cells_500m": (30, 500), "sens_30cells_750m": (30, 750)})
MAIN = "main_30cells_750m" if SCR == "g9" else "main_30cells_norelief"
E3 = [-INF, 1/3, 2/3, INF]; FIX = [0, 1, 2, 4, 8, INF]; FIX5 = [0, 5, 10, 15, INF]
REC = {"description": __doc__, "screens": {}}
for nm, (nc, rel) in SCREENS.items():
    ok = cr[(cr.n_peak_cells >= nc) & (cr.z95_rel >= rel)]; hh = ho[ho.uid.isin(ok.uid)]; ff = fd[fd.uid.isin(ok.uid)]
    REC["screens"][nm] = dict(min_peak_cells=nc, min_relief_m=None if not np.isfinite(rel) else rel, craters=len(ok), craters_with_peak_hollows=int(hh.uid.nunique()),
                              peak_hollows=len(hh), peak_fields=len(ff), median_relief_m=float(ok.z95_rel.median()), median_peak_cells=float(ok.n_peak_cells.median()))
REC["slope_distribution_deg"] = {k: [float(x) for x in np.nanpercentile(pk[k], [5, 25, 50, 75, 95])] for k in ("s1", "s2")}
REC["slope_distribution_deg"]["note"] = "5/25/50/75/95th percentiles over all valid peak cells (352 craters); s1 = 1-cell baseline, s2 = 2-cell baseline"

def prep(nc, rel, scol):
    ok = cr[(cr.n_peak_cells >= nc) & (cr.z95_rel >= rel)].reset_index(drop=True); idx = {u: i for i, u in enumerate(ok.uid)}
    c = pk[pk.uid.isin(ok.uid) & (pk.a > 0)].copy(); c["ci"] = c.uid.map(idx); c["s"] = c[scol]
    h = ho[ho.uid.isin(ok.uid)].copy(); h["ci"] = h.uid.map(idx); h["s"] = h[scol]
    f = fd[fd.uid.isin(ok.uid)].copy(); f["ci"] = f.uid.map(idx); f["s"] = f[scol]
    # within-crater weighted slope terciles
    ed = np.zeros((len(ok), 2))
    for i, g in c.groupby("ci"):
        o = np.argsort(g.s.values); cw = np.cumsum(g.a.values[o]); cw /= cw[-1]
        ed[i] = [g.s.values[o][np.searchsorted(cw, 1/3)], g.s.values[o][np.searchsorted(cw, 2/3)]]
    for t in (c, h, f):
        t["e3"] = np.digitize(t.eta.values, E3[1:-1]); t["st"] = (t.s.values > ed[t.ci.values, 0]).astype(int)+(t.s.values > ed[t.ci.values, 1]).astype(int)
        t["sf"] = np.digitize(t.s.values, FIX[1:-1]); t["sf5"] = np.digitize(t.s.values, FIX5[1:-1])
    return ok, c, h, f
def strata_tables(c, p, N, scol, ns):
    """A[c, s, e], n[c, s, e]: weighted area and counts by crater x slope stratum x eta third"""
    A = np.zeros((N, ns, 3)); n = np.zeros((N, ns, 3))
    sc = c[scol].values if scol else np.zeros(len(c), int); sp = p[scol].values if scol else np.zeros(len(p), int)
    np.add.at(A, (c.ci.values, sc, c.e3.values), c.a.values); np.add.at(n, (p.ci.values, sp, p.e3.values), 1)
    return A, n
def mh(A, n, inb, W):
    """MH rate ratio in-band vs out-of-band over crater x stratum; per-crater contributions so the bootstrap is a matrix product"""
    T1, T0 = A[:, :, inb].sum(2), A[:, :, ~inb].sum(2); a1, b0 = n[:, :, inb].sum(2), n[:, :, ~inb].sum(2); T = T1+T0
    with np.errstate(divide="ignore", invalid="ignore"):
        num = np.where(T > 0, a1*T0/T, 0); den = np.where(T > 0, b0*T1/T, 0)
    nc, dc = num.sum(1), den.sum(1); est = nc.sum()/dc.sum(); bs = (W@nc)/(W@dc); bs = bs[np.isfinite(bs)]
    per = []
    for s in range(A.shape[1]):
        ns_, ds_ = num[:, s].sum(), den[:, s].sum(); bss = (W@num[:, s])/(W@den[:, s]); bss = bss[np.isfinite(bss)]
        per.append(dict(RR=float(ns_/ds_) if ds_ > 0 else None, lo=float(np.percentile(bss, 2.5)) if len(bss) else None, hi=float(np.percentile(bss, 97.5)) if len(bss) else None,
                        n_in=int(a1[:, s].sum()), n_out=int(b0[:, s].sum()), area_frac_in=float(T1[:, s].sum()/T[:, s].sum()) if T[:, s].sum() > 0 else None))
    return dict(RR=float(est), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), per_stratum=per)
def bootW(N, seed=1):
    rng = np.random.default_rng(seed); return np.stack([np.bincount(rng.integers(0, N, N), minlength=N) for _ in range(NB)]).astype(float)
# ---- conditional Poisson (= Poisson with crater fixed effects) ----
def design(A, n, use_eta=True, use_slope=True):
    N, ns, _ = A.shape; rows = []
    for ci in range(N):
        for s in range(ns):
            for e in range(3):
                if A[ci, s, e] > 0: rows.append((ci, s, e, A[ci, s, e], n[ci, s, e]))
    r = np.array(rows); g = r[:, 0].astype(int); s = r[:, 1].astype(int); e = r[:, 2].astype(int)
    cols, names = [], []
    if use_eta: cols += [(e == 0).astype(float), (e == 2).astype(float)]; names += ["eta_low_vs_mid", "eta_upper_vs_mid"]
    if use_slope:
        for k in range(1, ns): cols.append((s == k).astype(float)); names.append(f"slope_bin{k}_vs_bin0")
    return g, np.column_stack(cols), r[:, 4], np.log(r[:, 3]), names
def cpois(g, X, y, off, w=None, it=50):
    N = g.max()+1; w = np.ones(N) if w is None else w; b = np.zeros(X.shape[1]); wr = w[g]
    for _ in range(it):
        eta = X@b+off; mx = np.zeros(N); np.maximum.at(mx, g, eta); ex = np.exp(eta-mx[g])
        Z = np.bincount(g, ex, N); p = ex/Z[g]; Y = np.bincount(g, y, N)
        Yr = Y[g]; mu = Yr*p
        grad = (wr*(y-mu))@X
        PX = np.stack([np.bincount(g, p*X[:, j], N) for j in range(X.shape[1])], 1)
        H = (X*(wr*mu)[:, None]).T@X-(PX*(w*Y)[:, None]).T@PX
        step = np.linalg.solve(H+1e-10*np.eye(len(b)), grad); b += step
        if np.max(np.abs(step)) < 1e-9: break
    return b, H, (y-mu)[:, None]*X - 0, PX, Y, p
def fit_cp(A, n, W, use_eta=True, use_slope=True):
    g, X, y, off, names = design(A, n, use_eta, use_slope)
    keep = np.bincount(g, y, g.max()+1)[g] > 0           # craters with no hollows carry no information
    gm = {v: i for i, v in enumerate(np.unique(g[keep]))}; g2 = np.array([gm[v] for v in g[keep]]); X, y, off = X[keep], y[keep], off[keep]
    b, H, _, PX, Y, p = cpois(g2, X, y, off)
    # cluster-robust sandwich (crater scores)
    N2 = g2.max()+1; S = np.stack([np.bincount(g2, (y-Y[g2]*p)*X[:, j], N2) for j in range(X.shape[1])], 1)
    Hi = np.linalg.inv(H); V = Hi@(S.T@S)@Hi*N2/(N2-1); se = np.sqrt(np.diag(V))
    Wk = W[:, np.unique(g[keep])]                        # bootstrap multiplicities of the informative craters
    bs = np.array([cpois(g2, X, y, off, w=Wk[i])[0] for i in range(len(Wk))])
    out = {}
    for j, nm in enumerate(names):
        q = bs[:, j][np.isfinite(bs[:, j])]
        out[nm] = dict(RR=float(np.exp(b[j])), robust_lo=float(np.exp(b[j]-1.96*se[j])), robust_hi=float(np.exp(b[j]+1.96*se[j])),
                       boot_lo=float(np.exp(np.percentile(q, 2.5))), boot_hi=float(np.exp(np.percentile(q, 97.5))))
    out["n_craters_informative"] = int(N2); out["n_counts"] = int(y.sum()); out["n_cells_strata"] = int(len(y))
    return out, (g2, X, y, off, b)
def analyse(nc, rel, scol, pts="hollows", full=True):
    ok, c, h, f = prep(nc, rel, scol); N = len(ok); W = bootW(N); p = h if pts == "hollows" else f
    R = dict(n_craters=N, n_points=len(p), n_craters_with_points=int(p.uid.nunique()))
    up = np.array([0, 0, 1], bool); lo = np.array([1, 0, 0], bool)
    A0, n0 = strata_tables(c, p, N, None, 1)
    R["MH_upper_vs_lower2_crater_only"] = mh(A0, n0, up, W)
    R["MH_lowest_vs_upper2_crater_only"] = mh(A0, n0, lo, W)
    for key, ns in (("st", 3), ("sf", len(FIX)-1)) + ((("sf5", len(FIX5)-1),) if full else ()):
        A, n = strata_tables(c, p, N, key, ns)
        R[f"MH_upper_vs_lower2_by_{key}"] = mh(A, n, up, W); R[f"MH_lowest_vs_upper2_by_{key}"] = mh(A, n, lo, W)
    # regression: crater x eta third x slope bin (fixed bins), and crater x eta third x within-crater slope tercile
    A, n = strata_tables(c, p, N, "sf", len(FIX)-1)
    R["poisson_eta_only"], _ = fit_cp(A, n, W, True, False)
    R["poisson_eta_plus_slopebins"], chk = fit_cp(A, n, W, True, True)
    R["poisson_slopebins_only"], _ = fit_cp(A, n, W, False, True)
    A, n = strata_tables(c, p, N, "st", 3); R["poisson_eta_plus_slopeterciles"], _ = fit_cp(A, n, W, True, True)
    # crude: hollow vs crater-matched mean slope and eta-slope correlation in this sample
    R["mean_slope_points"] = float(p.s.mean())
    Ac = c.groupby("ci").a.sum(); ms = (c.a*c.s).groupby(c.ci).sum()/Ac; cnt = p.groupby("ci").size()
    R["mean_slope_crater_matched_expected"] = float((ms.reindex(cnt.index)*cnt).sum()/cnt.sum())
    return R, chk, (ok, c, h, f)
def corr(c):
    w = c.a.values; x = c.eta.values; y = c.s.values
    def wc(x, y, w):
        mx, my = np.average(x, weights=w), np.average(y, weights=w)
        return float(np.average((x-mx)*(y-my), weights=w)/np.sqrt(np.average((x-mx)**2, weights=w)*np.average((y-my)**2, weights=w)))
    xd = x-pd.Series(x*w).groupby(c.ci.values).transform("sum").values/pd.Series(w).groupby(c.ci.values).transform("sum").values
    yd = y-pd.Series(y*w).groupby(c.ci.values).transform("sum").values/pd.Series(w).groupby(c.ci.values).transform("sum").values
    sp = [pd.Series(g.eta.values).rank().corr(pd.Series(g.s.values).rank()) for _, g in c.groupby("ci") if len(g) >= 10]
    rx = pd.Series(x).rank().values; ry = pd.Series(y).rank().values
    # mean slope by eta third (weighted)
    ms = [float(np.average(y[c.e3.values == k], weights=w[c.e3.values == k])) if (c.e3.values == k).any() else None for k in range(3)]
    return dict(pearson_weighted_pooled=wc(x, y, w), pearson_weighted_within_crater=wc(xd, yd, w), spearman_pooled=float(np.corrcoef(rx, ry)[0, 1]),
                spearman_per_crater_median=float(np.median(sp)), spearman_per_crater_iqr=[float(np.percentile(sp, 25)), float(np.percentile(sp, 75))],
                mean_slope_by_eta_third=ms)
RES = {}
for scol in ("s1", "s2"):
    for nm, (nc, rel) in SCREENS.items():
        if scol == "s2" and nm != MAIN: continue
        for pts in ("hollows", "fields"):
            full = (nm == MAIN)
            r, chk, (ok, c, h, f) = analyse(nc, rel, scol, pts, full)
            if full or SCR == "g11": r["eta_slope_correlation"] = corr(c)
            RES[f"{scol}|{nm}|{pts}"] = r
            if scol == "s1" and nm == MAIN and pts == "hollows": CHK = chk; FIGDATA = (ok, c, h, f)
REC["results"] = RES
# one check: statsmodels GLM Poisson with crater dummies and offset reproduces the conditional-Poisson coefficients
import statsmodels.api as sm
g2, X, y, off, b = CHK; Dm = np.eye(g2.max()+1)[g2]
m = sm.GLM(y, np.column_stack([X, Dm]), family=sm.families.Poisson(), offset=off).fit()
REC["check_glm_dummies_vs_conditional"] = dict(conditional=[float(x) for x in b], glm=[float(x) for x in m.params[:X.shape[1]]],
                                               pearson_chi2_over_df=float(m.pearson_chi2/m.df_resid))
REC["bins"] = dict(eta_thirds=["<1/3", "1/3-2/3", ">=2/3"], slope_fixed_deg=FIX[:-1], slope_fixed5_deg=FIX5[:-1], slope_terciles="within-crater, completeness-weighted")
REC["g11"] = dict(eta=ETA, screens_set=SCR, main_screen=MAIN, rules=RULES)
json.dump(REC, open(f"{HERE}/positions_slope_record_{SUF}.json", "w"), indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else o)
# figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
BLUE, ORANGE, INK, INK2, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d9d8d4"
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
ok, c, h, f = FIGDATA; M = RES[f"s1|{MAIN}|hollows"]; MF = RES[f"s1|{MAIN}|fields"]
fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.8), gridspec_kw=dict(width_ratios=[1, 1.05, 1.25]))
a = ax[0]; eb = np.linspace(0.2, 1.2, 11); xc = (eb[:-1]+eb[1:])/2
for lab, col in (("1-cell baseline", BLUE), ("2-cell baseline", ORANGE)):
    s = c.s1 if col == BLUE else c.s2; med, q1, q3 = [], [], []
    for lo_, hi_ in zip(eb[:-1], eb[1:]):
        k = (c.eta >= lo_) & (c.eta < hi_); med.append(np.median(s[k])); q1.append(np.percentile(s[k], 25)); q3.append(np.percentile(s[k], 75))
    a.plot(xc, med, color=col, marker="o", ms=4, lw=1.5, label=f"median, {lab}")
    if col == BLUE: a.fill_between(xc, q1, q3, color=col, alpha=0.15, lw=0, label="interquartile range, 1-cell")
r_ = M["eta_slope_correlation"]
a.set_xlabel("η = (z − z_floor)/(z_peak95 − z_floor)" if ETA == "eta" else "η_b = (z − z_base)/(z_peak95 − z_base)"); a.set_ylabel("DEM slope (°)"); a.yaxis.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
a.set_title(f"(a) Slope vs height in peak cells ({M['n_craters']} screened craters)", loc="left", fontsize=9.5)
a.text(0.02, 0.97, f"within-crater r = {r_['pearson_weighted_within_crater']:.2f}\nmedian per-crater Spearman = {r_['spearman_per_crater_median']:.2f}", transform=a.transAxes, va="top", color=INK2, fontsize=8)
a.legend(frameon=False, fontsize=7.5, loc="lower right")
a = ax[1]; per = M["MH_upper_vs_lower2_by_sf"]["per_stratum"]; labs = ["0–1°", "1–2°", "2–4°", "4–8°", ">8°"]
x = np.arange(len(per)); rr = np.array([p_["RR"] if p_["RR"] else np.nan for p_ in per], float)
lo_ = np.array([p_["lo"] if p_["lo"] else np.nan for p_ in per], float); hi_ = np.array([p_["hi"] if p_["hi"] else np.nan for p_ in per], float)
a.axhline(1, color=INK2, lw=0.9); a.errorbar(x, rr, yerr=[rr-lo_, hi_-rr], fmt="o", color=BLUE, ms=5, lw=1.2)

MHX = []
for j, (k, lab) in enumerate((("MH_upper_vs_lower2_crater_only", "crater only"), ("MH_upper_vs_lower2_by_sf", "crater × slope bin"), ("MH_upper_vs_lower2_by_st", "crater × slope tercile"))):
    v = M[k]; xx = len(per)+0.7+j*1.15; a.errorbar([xx], [v["RR"]], yerr=[[v["RR"]-v["lo"]], [v["hi"]-v["RR"]]], fmt="D", color=INK, ms=5, lw=1.2); MHX.append((xx, lab))
a.set_xticks(list(x)+[m_[0] for m_ in MHX]); a.set_xticklabels([f"{l_}\n{p_['n_in']+p_['n_out']}" for l_, p_ in zip(labs, per)]+["MH\nno\nslope", "MH\nslope\nbins", "MH\nslope\nterciles"], fontsize=7.5)
a.set_yscale("log"); a.set_yticks([0.5, 1, 2, 4, 8]); a.set_yticklabels(["0.5", "1", "2", "4", "8"]); a.set_ylim(0.4, 9); a.yaxis.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
a.axvline(len(per)-0.2, color=GRID, lw=1)
a.set_ylabel("rate ratio, upper third of η vs lower two thirds"); a.set_xlabel("DEM slope stratum, 1-cell baseline (hollows per stratum below)")
a.set_title("(b) Height effect within slope strata, hollows", loc="left", fontsize=9.5)
a = ax[2]; items = []
for lab, key, src in (("η upper vs middle, no slope", ("poisson_eta_only", "eta_upper_vs_mid"), M), ("η upper vs middle, + slope bins", ("poisson_eta_plus_slopebins", "eta_upper_vs_mid"), M),
                      ("η upper vs middle, + slope terciles", ("poisson_eta_plus_slopeterciles", "eta_upper_vs_mid"), M),
                      ("η low vs middle, no slope", ("poisson_eta_only", "eta_low_vs_mid"), M), ("η low vs middle, + slope bins", ("poisson_eta_plus_slopebins", "eta_low_vs_mid"), M),
                      ("slope 1–2° vs 0–1° (with η)", ("poisson_eta_plus_slopebins", "slope_bin1_vs_bin0"), M), ("slope 2–4° vs 0–1°", ("poisson_eta_plus_slopebins", "slope_bin2_vs_bin0"), M),
                      ("slope 4–8° vs 0–1°", ("poisson_eta_plus_slopebins", "slope_bin3_vs_bin0"), M), ("slope > 8° vs 0–1°", ("poisson_eta_plus_slopebins", "slope_bin4_vs_bin0"), M)):
    for off_, S, col in ((-0.15, M, BLUE), (0.15, MF, ORANGE)):
        v = S[key[0]][key[1]]; items.append((lab, off_, v, col))
labs = list(dict.fromkeys(i[0] for i in items)); ypos = {l: len(labs)-1-k for k, l in enumerate(labs)}
a.axvline(1, color=INK2, lw=0.9)
for lab, off_, v, col in items:
    a.errorbar([v["RR"]], [ypos[lab]+off_], xerr=[[v["RR"]-v["boot_lo"]], [v["boot_hi"]-v["RR"]]], fmt="o", color=col, ms=4, lw=1.2)
a.set_yticks(list(ypos.values())); a.set_yticklabels(list(ypos.keys()), fontsize=7.5); a.set_xscale("log"); a.set_xlabel("rate ratio (crater fixed effects, 95% crater bootstrap)")
a.set_xticks([0.25, 0.5, 1, 2, 4]); a.set_xticklabels(["0.25", "0.5", "1", "2", "4"]); a.minorticks_off(); a.set_xlim(0.18, 6.5)
a.xaxis.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
a.plot([], [], "o", color=BLUE, label=f"hollows (n = {M['n_points']})"); a.plot([], [], "o", color=ORANGE, label=f"fields (n = {MF['n_points']})"); a.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
a.set_title("(c) Poisson regression", loc="left", fontsize=9.5)
fig.tight_layout(); fig.savefig(f"{HERE}/fig_positions_slope_{SUF}.png", dpi=200, facecolor="white", bbox_inches="tight")
def P(k):
    r = RES[k]; s = f"{k}: craters {r['n_craters']} (with pts {r['n_craters_with_points']}), pts {r['n_points']}\n"
    for m_ in [x for x in r if x.startswith("MH")]: s += f"   {m_}: {r[m_]['RR']:.2f} ({r[m_]['lo']:.2f}-{r[m_]['hi']:.2f})\n"
    for m_ in [x for x in r if x.startswith("poisson")]:
        s += f"   {m_}: " + "; ".join(f"{kk} {vv['RR']:.2f} [rob {vv['robust_lo']:.2f}-{vv['robust_hi']:.2f}; boot {vv['boot_lo']:.2f}-{vv['boot_hi']:.2f}]" for kk, vv in r[m_].items() if isinstance(vv, dict)) + f" (inf. craters {r[m_]['n_craters_informative']})\n"
    s += f"   mean slope pts {r['mean_slope_points']:.2f} vs expected {r['mean_slope_crater_matched_expected']:.2f}\n"
    return s
for k in RES: print(P(k))
print(json.dumps(REC["screens"])); print(json.dumps(RES[f"s1|{MAIN}|hollows"]["eta_slope_correlation"])); print(json.dumps(RES[f"s2|{MAIN}|hollows"]["eta_slope_correlation"]))
print(json.dumps(REC["check_glm_dummies_vs_conditional"])); print(json.dumps(REC["slope_distribution_deg"]))
