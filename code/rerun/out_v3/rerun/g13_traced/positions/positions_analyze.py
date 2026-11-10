# g13 revision (7 Nov 2026): with G11_RULES=traced and G13_REV=r2 (default) outputs go to out_traced_r2; G13_REV=r1 keeps out_traced.
# g13 copy of g11_curvature/positions/positions_analyze.py (6 Nov 2026). Changed: folder g13_traced/positions; G11_RULES default traced (traced outlines where available, curvature v3 elsewhere via g13_traced/traced_classes.py; env CLS=v2 for curvature v2 elsewhere).
# g11 copy of g9_draft_reruns/positions/positions_analyze.py (6 Nov 2026). Changed: env G11_RULES (curv | fix) and paths (reads out_<RULES>/cache or G11_CACHE_DIR, writes out_<RULES>/); the H150/H400 threshold sensitivity runs only if those caches exist (the curvature classes have no relief threshold); new measures in REC['g11_additions']: the eta tests repeated with eta_b (height measured from the 5th percentile of peak-cell elevations), also by crater morphology, and floor near the wall alone and near the peak alone within 0.05 R and 0.1 R (plain and radius-stratified), for hollows and fields.
# g9 copy of g3_catalog/melt_source/positions/positions_analyze.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID and paths (reads <out>/cache or G9_CACHE_DIR), radial bins RB extended to 1.0 R under RULES=fix (the new floor class reaches 1.0 R), and new measures in REC['g9_additions']: floor margins within 0.05 R and 0.1 R of wall or peak cells (distance / host radius), their radius-stratified versions, and floor within 0.2 R of the center vs the rest of the floor as a split ratio.
# g3_catalog copy of analysis_workshop/melt_source_tests/positions/positions_analyze.py (6 Nov 2026): only input/output paths and env settings changed.
#!/usr/bin/env python3
"""Melt-source position test, step 2 (5 Nov 2026): where do hollows sit within peak and floor terrain?
Reads cache/H{150,250,400}_000_354.npz from positions_build.py. Each hollow (or hollow field, counted once at its centroid)
is compared with the completeness-weighted area of the same terrain class.
Two nulls: (i) pooled area -- density ratio DR = (n_in/A_in)/(n_out/A_out) summed over craters; (ii) crater-matched --
each crater's hollows are spread over that crater's own weighted-area distribution (expected E_in = sum_c n_c A_in,c/A_c),
giving the stratified ratio SDR = (O_in/E_in)/(O_out/E_out), which removes between-crater differences.
95% intervals: 2,000 bootstrap resamples of craters (seed 1). Position uncertainty: 25 jittered hollow positions
(0.5 km per axis) from the build step. Writes positions_record.json, hollow_positions.csv and fig_positions.png."""
import os, json, numpy as np, pandas as pd
import os as _os
import sys as _sys; _sys.dont_write_bytecode = True
G11 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced")
CAT = "v3"; RULES = _os.environ.get("G11_RULES", "traced"); GRID = "v3"; assert RULES in ("curv", "fix", "traced"), RULES
ODIR = f"{G11}/positions/out_{RULES}" + ("_r2" if RULES == "traced" and _os.environ.get("G13_REV", "r2") == "r2" else ""); _os.makedirs(ODIR + "/cache", exist_ok=True)
CACHE_DIR = _os.environ.get("G11_CACHE_DIR", ODIR + "/cache")
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))
HERE = ODIR; HE = os.path.expanduser("~/mnt/hollows-emissivity")
NB = 2000; R_M = 2439.7
HCOL = ["uid", "hid", "fid", "cls", "eta", "q", "d_base", "d_wall", "d_peak", "d_elev", "d_nonf", "rho"]
JCOL = ["cls", "eta", "q", "d_base", "d_wall", "d_peak", "d_elev", "d_nonf", "rho"]
FCOL = ["uid", "fid", "cls", "n", "fb", "eta", "q", "d_base", "d_wall", "d_peak", "d_elev", "d_nonf", "rho"]
def load(H):
    d = np.load(f"{CACHE_DIR}/H{H}_000_354.npz", allow_pickle=True)
    cr = pd.DataFrame(d["craters"]); cr = cr[(cr.D >= 30) & (cr.D <= 250)].reset_index(drop=True)
    pk = pd.DataFrame(d["peak"], columns=["uid", "a", "h", "eta", "q", "d_base", "rho"])
    fl = pd.DataFrame(d["floor"], columns=["uid", "a", "h", "rho", "d_wall", "d_peak", "d_elev", "d_nonf"])
    ho = pd.DataFrame(d["hollows"], columns=HCOL); fd = pd.DataFrame(d["fields"], columns=FCOL)
    fl["d_marg"] = np.minimum(fl.d_wall, fl.d_peak); ho["d_marg"] = np.minimum(ho.d_wall, ho.d_peak); fd["d_marg"] = np.minimum(fd.d_wall, fd.d_peak)
    Rk = dict(zip(cr.uid.astype(int), cr.D/2))
    for t in (fl, ho, fd): t["d_marg_R"] = t.d_marg/t.uid.astype(int).map(Rk)
    for t in (fl, ho, fd): t["d_wall_R"] = t.d_wall/t.uid.astype(int).map(Rk); t["d_peak_R"] = t.d_peak/t.uid.astype(int).map(Rk)
    if "peak_etab" in d.files: pk["eta_b"] = d["peak_etab"]; ho["eta_b"] = d["hollows_etab"]; fd["eta_b"] = d["fields_etab"]
    return cr, pk, fl, ho, fd, d["jitter"]
cr, pk, fl, ho, fd, jit = load(250)
N = len(cr); cidx = {u: i for i, u in enumerate(cr.uid.astype(int))}
for t in (pk, fl, ho, fd): t["ci"] = t.uid.astype(int).map(cidx)
# crater attributes: morphology, vent, size
cs = cr.cent_str.astype(str).str.strip()
cr["morph"] = np.where(cs.isin(["cp", "mp"]), "central_peak", np.where(cs.isin(["pr", "rpc", "pb"]), "peak_ring", "unclassified"))
v = pd.read_csv(f"{HE}/data/data/MER_Vent_Locations_Classifications.csv"); v.columns = [c.strip() for c in v.columns]
vlat, vlon = np.radians(v["Lon"].values), np.radians(v["Lat"].values % 360)   # header swap, as in test8_vents.py
def has_vent(la, lo, R):
    la, lo = np.radians(la), np.radians(lo)
    d = R_M*np.arccos(np.clip(np.sin(la)*np.sin(vlat)+np.cos(la)*np.cos(vlat)*np.cos(vlon-lo), -1, 1)); return bool((d <= R).any())
cr["vent"] = [has_vent(r.lat, r.lon, r.D/2) for r in cr.itertuples()]
cr["big"] = cr.D >= 100

def boot_weights(sel, seed=1):
    rng = np.random.default_rng(seed); ii = np.where(sel)[0]; W = np.zeros((NB, N))
    for b in range(NB): W[b] += np.bincount(ii[rng.integers(0, len(ii), len(ii))], minlength=N)
    return W
def mats(cells, pts, col, edges, cells_keep=None, pts_keep=None):
    """per-crater (N x nbins) weighted area and point counts"""
    c = cells if cells_keep is None else cells[cells_keep]; p = pts if pts_keep is None else pts[pts_keep]
    nb = len(edges)-1
    bc = np.clip(np.digitize(c[col].values, edges)-1, 0, nb-1); bp = np.clip(np.digitize(p[col].values, edges)-1, 0, nb-1)
    okc = np.isfinite(c[col].values); okp = np.isfinite(p[col].values)
    A = np.zeros((N, nb)); np.add.at(A, (c.ci.values[okc], bc[okc]), c.a.values[okc])
    n = np.zeros((N, nb)); np.add.at(n, (p.ci.values[okp], bp[okp]), 1)
    return A, n
def expected(A, n):
    At = A.sum(1, keepdims=True); ok = At[:, 0] > 0
    E = np.where(ok[:, None], n.sum(1, keepdims=True)*A/np.where(At > 0, At, 1), 0); O = np.where(ok[:, None], n, 0); return O, E
def split_stats(A, n, inb, sel, W):
    """inb: boolean over bins = 'in band'. Returns pooled DR and crater-matched SDR with bootstrap intervals"""
    w0 = sel.astype(float); O, E = expected(A, n)
    def f(w):
        Ai, Ao = w@A[:, inb].sum(1), w@A[:, ~inb].sum(1); ni, no = w@n[:, inb].sum(1), w@n[:, ~inb].sum(1)
        Oi, Oo = w@O[:, inb].sum(1), w@O[:, ~inb].sum(1); Ei, Eo = w@E[:, inb].sum(1), w@E[:, ~inb].sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            return dict(DR=(ni/Ai)/(no/Ao), SDR=(Oi/Ei)/(Oo/Eo), OE_in=Oi/Ei, f_h=ni/(ni+no), f_A=Ai/(Ai+Ao), f_E=Ei/(Ei+Eo))
    est = f(w0[None, :]); bs = f(W)
    out = dict(n_craters=int(sel.sum()), n_craters_with_area=int(((A.sum(1) > 0) & sel).sum()),
               n_craters_with_points=int(((n.sum(1) > 0) & sel).sum()), n_points=int(n[sel].sum()), n_in=int(n[sel][:, inb].sum()),
               frac_points_in=float(est["f_h"][0]), frac_area_in_pooled=float(est["f_A"][0]), frac_expected_in_crater_matched=float(est["f_E"][0]))
    for k in ("DR", "SDR", "OE_in"):
        b = bs[k][np.isfinite(bs[k])]
        out[k] = [float(est[k][0]), float(np.percentile(b, 2.5)) if len(b) else None, float(np.percentile(b, 97.5)) if len(b) else None]
    return out
def tercile_fracs(A, n, sel):
    O, E = expected(A, n); s = sel
    return dict(points=[float(x) for x in n[s].sum(0)/n[s].sum()], crater_matched_area=[float(x) for x in E[s].sum(0)/E[s].sum()],
                pooled_area=[float(x) for x in A[s].sum(0)/A[s].sum()], n_points=[int(x) for x in n[s].sum(0)])
INF = np.inf
E_ETA = [0, 1/3, 2/3, INF]; E_Q = [0, 1/3, 2/3, 1.0001]; E_DB = [0, 2, INF]; E_DM = [0, 3, INF]
def peak_floor_tests(pk, fl, pts_p, pts_f, sel, W, excl=False):
    """the four headline tests; pts_p/pts_f: point tables (hollows or fields) in peak / floor class"""
    R = {}
    kc = (pk.d_base > 0.7) if excl else None; kp = (pts_p.d_base > 0.7) if excl else None
    A, n = mats(pk, pts_p, "eta", E_ETA, kc, kp); R["peak_eta_low_vs_upper"] = split_stats(A, n, np.array([1, 0, 0], bool), sel, W)
    R["peak_eta_high_vs_lower"] = split_stats(A, n, np.array([0, 0, 1], bool), sel, W); R["peak_eta_terciles"] = tercile_fracs(A, n, sel)
    A, n = mats(pk, pts_p, "d_base", [0.7, 2, INF] if excl else E_DB, kc, kp); R["peak_base_within2km"] = split_stats(A, n, np.array([1, 0], bool), sel, W)
    A, n = mats(pk, pts_p, "q", E_Q, kc, kp); R["peak_q_lowest_third"] = split_stats(A, n, np.array([1, 0, 0], bool), sel, W)
    R["peak_q_terciles"] = tercile_fracs(A, n, sel)
    kc = (fl.d_marg > 0.7) if excl else None; kp = (pts_f.d_marg > 0.7) if excl else None
    A, n = mats(fl, pts_f, "d_marg", [0.7, 3, INF] if excl else E_DM, kc, kp); R["floor_margin_within3km_of_wall_or_peak"] = split_stats(A, n, np.array([1, 0], bool), sel, W)
    if not excl:
        for col, nm in (("d_wall", "floor_within3km_of_wall"), ("d_peak", "floor_within3km_of_peak"), ("d_elev", "floor_within3km_of_any_high_ground"),
                        ("d_nonf", "floor_within3km_of_any_nonfloor_cell")):
            A, n = mats(fl, pts_f, col, E_DM); R[nm] = split_stats(A, n, np.array([1, 0], bool), sel, W)
        # floor-peak-base test restricted to craters with peak terrain
        hp = sel & (cr.A_peak.values > 0); A, n = mats(fl, pts_f, "d_peak", E_DM)
        R["floor_within3km_of_peak_craters_with_peaks"] = split_stats(A, n, np.array([1, 0], bool), hp, boot_weights(hp))
    return R
REC = {"description": __doc__, "n_host_craters": N, "terrain_threshold_m": 250,
       "bands": {"peak_low": "eta < 1/3 (eta = (z - z_floor)/(z_peak95 - z_floor))", "peak_base": "distance to nearest floor-class cell <= 2 km",
                 "q": "completeness-weighted height quantile within the crater's peak terrain", "floor_margin": "min(distance to wall-class, to peak-class cell) <= 3 km"},
       "ratio_definitions": {"DR": "pooled density ratio in-band / out-of-band", "SDR": "crater-matched (stratified) density ratio",
                             "OE_in": "observed / crater-matched expected count in band"}}
ALL = np.ones(N, bool); W_ALL = boot_weights(ALL)
hp_, hf_ = ho[ho.cls == 1], ho[ho.cls == 3]; fp_, ff_ = fd[fd.cls == 1], fd[fd.cls == 3]
REC["sample"] = dict(craters=N, craters_with_peak_terrain=int((cr.A_peak > 0).sum()), peak_hollows=len(hp_), floor_hollows=len(hf_),
                     craters_with_peak_hollows=int(hp_.uid.nunique()), craters_with_floor_hollows=int(hf_.uid.nunique()),
                     peak_fields=len(fp_), floor_fields=len(ff_), field_centroid_fallbacks_peak=int(fp_.fb.sum()), field_centroid_fallbacks_floor=int(ff_.fb.sum()),
                     peak_hollows_with_eta_lt_1_3=int((hp_.eta < 1/3).sum()), craters_with_low_eta_area=int(pk[pk.eta < 1/3].uid.nunique()),
                     median_peak_relief_z95_m_hollow_hosts=float(cr[cr.uid.isin(hp_.uid)].z95_rel.median()))
REC["hollows"] = peak_floor_tests(pk, fl, hp_, hf_, ALL, W_ALL)
REC["fields"] = peak_floor_tests(pk, fl, fp_, ff_, ALL, W_ALL)
REC["mean_q_peak_hollows"] = float(hp_.q.mean()); 
_qb = [] ; _g = hp_.groupby("ci").q.agg(["sum", "size"]); _s = np.zeros(N); _c = np.zeros(N); _s[_g.index] = _g["sum"]; _c[_g.index] = _g["size"]
_b = (W_ALL@_s)/(W_ALL@_c); REC["mean_q_peak_hollows"] = [float(_s.sum()/_c.sum()), float(np.percentile(_b, 2.5)), float(np.percentile(_b, 97.5))]
REC["mean_q_peak_fields"] = float(fp_.q.mean())
# boundary exclusion (one DEM cell, 0.7 km)
REC["sensitivity_exclude_within_0.7km_of_boundary"] = peak_floor_tests(pk, fl, hp_, hf_, ALL, W_ALL, excl=True)
# breakdowns
BR = {}
for nm, s in [("central_peak", cr.morph == "central_peak"), ("peak_ring", cr.morph == "peak_ring"), ("unclassified_center", cr.morph == "unclassified"),
              ("vent_inside_rim", cr.vent), ("no_vent", ~cr.vent), ("D_ge_100km", cr.big), ("D_lt_100km", ~cr.big)]:
    s = s.values; W = boot_weights(s); BR[nm] = dict(hollows=peak_floor_tests(pk, fl, hp_, hf_, s, W), fields=peak_floor_tests(pk, fl, fp_, ff_, s, W))
REC["breakdowns"] = BR
# jitter (hollow position uncertainty)
JR = {k: [] for k in ("peak_eta_low_DR", "peak_eta_low_SDR", "peak_base2_DR", "peak_base2_SDR", "peak_q_low_frac", "floor_margin3_DR", "floor_margin3_SDR", "n_peak", "n_floor")}
K = jit.shape[1]
for k in range(K):
    j = pd.DataFrame(jit[:, k, :].T, columns=JCOL); j["uid"] = ho.uid.values; j["ci"] = ho.ci.values; j["d_marg"] = np.minimum(j.d_wall, j.d_peak)
    jp, jf = j[j.cls == 1], j[j.cls == 3]; one = np.ones((1, N))
    A, n = mats(pk, jp, "eta", E_ETA); s = split_stats(A, n, np.array([1, 0, 0], bool), ALL, one); JR["peak_eta_low_DR"].append(s["DR"][0]); JR["peak_eta_low_SDR"].append(s["SDR"][0])
    A, n = mats(pk, jp, "d_base", E_DB); s = split_stats(A, n, np.array([1, 0], bool), ALL, one); JR["peak_base2_DR"].append(s["DR"][0]); JR["peak_base2_SDR"].append(s["SDR"][0])
    A, n = mats(pk, jp, "q", E_Q); JR["peak_q_low_frac"].append(float(n[:, 0].sum()/n.sum()))
    A, n = mats(fl, jf, "d_marg", E_DM); s = split_stats(A, n, np.array([1, 0], bool), ALL, one); JR["floor_margin3_DR"].append(s["DR"][0]); JR["floor_margin3_SDR"].append(s["SDR"][0])
    JR["n_peak"].append(len(jp)); JR["n_floor"].append(len(jf))
REC["sensitivity_position_jitter_0.5km"] = {k: dict(median=float(np.median(v)), p2_5=float(np.percentile(v, 2.5)), p97_5=float(np.percentile(v, 97.5))) for k, v in JR.items()}
# threshold sensitivity
TH = {}
for H in [H_ for H_ in (150, 400) if os.path.exists(f"{CACHE_DIR}/H{H_}_000_354.npz")]:
    c2, p2, f2, h2, d2, _ = load(H); assert list(c2.uid) == list(cr.uid)
    for t in (p2, f2, h2, d2): t["ci"] = t.uid.astype(int).map(cidx)
    TH[f"H{H}"] = dict(n_peak_hollows=int((h2.cls == 1).sum()), n_floor_hollows=int((h2.cls == 3).sum()),
                       hollows=peak_floor_tests(p2, f2, h2[h2.cls == 1], h2[h2.cls == 3], ALL, W_ALL),
                       fields=peak_floor_tests(p2, f2, d2[d2.cls == 1], d2[d2.cls == 3], ALL, W_ALL))
REC["sensitivity_threshold"] = TH
# profiles for the figure: fine bins, crater-matched O/E with bootstrap intervals
def profile(cells, pts, col, edges):
    A, n = mats(cells, pts, col, edges); O, E = expected(A, n)
    oe = O.sum(0)/E.sum(0); bo = (W_ALL@O)/(W_ALL@E)
    return dict(edges=[float(x) if np.isfinite(x) else None for x in edges], frac_points=(n.sum(0)/n.sum()).tolist(),
                frac_crater_matched=(E.sum(0)/E.sum()).tolist(), frac_pooled_area=(A.sum(0)/A.sum()).tolist(), n=n.sum(0).astype(int).tolist(),
                OE=oe.tolist(), OE_lo=np.nanpercentile(bo, 2.5, axis=0).tolist(), OE_hi=np.nanpercentile(bo, 97.5, axis=0).tolist())
PR = dict(eta=profile(pk, hp_, "eta", [0, 1/6, 1/3, 1/2, 2/3, 5/6, 1.0, 7/6, INF]),
          d_base=profile(pk, hp_, "d_base", [0, 1, 2, 3, 4, 6, 8, 10, 13, 16, INF]),
          d_marg=profile(fl, hf_, "d_marg", [0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 20, INF]))
REC["profiles"] = PR
# ---- added checks: signed profile across the peak base; radius-stratified floor tests; vent-crater detail ----
def signed_tables(pk, fl, hp, hf):
    c = pd.concat([pd.DataFrame(dict(ci=pk.ci.values, a=pk.a.values, s=pk.d_base.values)),
                   pd.DataFrame(dict(ci=fl.ci.values, a=fl.a.values, s=-fl.d_peak.values))], ignore_index=True)
    p = pd.concat([pd.DataFrame(dict(ci=hp.ci.values, s=hp.d_base.values)), pd.DataFrame(dict(ci=hf.ci.values, s=-hf.d_peak.values))], ignore_index=True)
    return c[np.isfinite(c.s)], p[np.isfinite(p.s)]
SE = [-INF, -10, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 10, INF]
sc, sp = signed_tables(pk, fl, hp_, hf_); hpk = (cr.A_peak.values > 0)
REC["profiles"]["signed_peak_base"] = profile(sc, sp, "s", SE)
REC["profiles"]["signed_peak_base"]["note"] = "craters with peak terrain; negative = floor cells (distance to nearest peak cell), positive = peak cells (distance to nearest floor cell); crater-matched null over peak+floor"
sfc, sfp = signed_tables(pk, fl, fp_, ff_); REC["profiles"]["signed_peak_base_fields"] = profile(sfc, sfp, "s", SE)
# radius-stratified floor tests: strata = crater x rho bin
RB = [0, 0.2, 0.4, 0.6, 0.81] if RULES == "old" else [0, 0.2, 0.4, 0.6, 0.8, 1.01]; NRB = len(RB)-1
def strat_split(cells, pts, col, thr, sel, W):
    rc = np.clip(np.digitize(cells.rho.values, RB)-1, 0, NRB-1); rp = np.clip(np.digitize(pts.rho.values, RB)-1, 0, NRB-1)
    A = np.zeros((N*NRB, 2)); n = np.zeros((N*NRB, 2))
    ok = np.isfinite(cells[col].values); np.add.at(A, (cells.ci.values[ok]*NRB+rc[ok], (cells[col].values[ok] > thr).astype(int)), cells.a.values[ok])
    ok = np.isfinite(pts[col].values); np.add.at(n, (pts.ci.values[ok]*NRB+rp[ok], (pts[col].values[ok] > thr).astype(int)), 1)
    s = split_stats(A, n, np.array([1, 0], bool), np.repeat(sel, NRB), np.repeat(W, NRB, axis=1))
    return {k: v for k, v in s.items() if k not in ("n_craters", "n_craters_with_area", "n_craters_with_points")}
REC["floor_radius_stratified"] = {"note": "crater-matched null within crater x radial bins (0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8 R)",
    "within3km_of_wall_or_peak": strat_split(fl, hf_, "d_marg", 3, ALL, W_ALL), "within3km_of_wall": strat_split(fl, hf_, "d_wall", 3, ALL, W_ALL),
    "within3km_of_peak": strat_split(fl, hf_, "d_peak", 3, ALL, W_ALL),
    "fields_within3km_of_wall_or_peak": strat_split(fl, ff_, "d_marg", 3, ALL, W_ALL), "fields_within3km_of_peak": strat_split(fl, ff_, "d_peak", 3, ALL, W_ALL)}
# floor density by radius (crater-matched O/E per radial bin)
A, n = mats(fl.assign(r=fl.rho), hf_.assign(r=hf_.rho), "r", RB); O, E = expected(A, n); bo = (W_ALL@O)/(W_ALL@E)
REC["floor_density_by_radius"] = dict(edges=RB, n=n.sum(0).astype(int).tolist(), OE=(O.sum(0)/E.sum(0)).tolist(),
                                      lo=np.percentile(bo, 2.5, axis=0).tolist(), hi=np.percentile(bo, 97.5, axis=0).tolist())
# per-crater detail for the vent subset
A1, n1 = mats(pk, hp_, "eta", E_ETA); O1, E1 = expected(A1, n1); A2, n2 = mats(pk, hp_, "d_base", E_DB); O2, E2 = expected(A2, n2)
vd = []
for i in np.where(cr.vent.values & (n1.sum(1) > 0))[0]:
    vd.append(dict(uid=int(cr.uid[i]), D=float(cr.D[i]), lat=float(cr.lat[i]), lon=float(cr.lon[i]), morph=cr.morph[i], n_peak=int(n1[i].sum()),
                   n_eta_low=int(n1[i, 0]), E_eta_low=float(E1[i, 0]), n_base2=int(n2[i, 0]), E_base2=float(E2[i, 0])))
REC["vent_craters_with_peak_hollows"] = sorted(vd, key=lambda r: -r["n_peak"])
# leave-one-out: the low-band SDR in vent craters without the crater with the most peak hollows
if vd:
    top = cidx[REC["vent_craters_with_peak_hollows"][0]["uid"]]; s = cr.vent.values.copy(); s[top] = False
    REC["vent_low_band_without_top_crater"] = dict(dropped_uid=int(cr.uid[top]), eta_low=split_stats(A1, n1, np.array([1, 0, 0], bool), s, boot_weights(s)),
                                                   base2=split_stats(A2, n2, np.array([1, 0], bool), s, boot_weights(s)))
# ---- g9 additions (6 Nov 2026): radius-scaled floor margins (distance / host radius) and the central-floor split ----
G9X = {"note": "floor hollows/fields vs completeness-weighted floor area; d_marg_R = min(distance to wall-class, to peak-class cell) / host radius; "
               "radius_stratified = crater-matched null within crater x radial bins (RB)", "RB": RB}
for lab, pts_f in (("hollows", hf_), ("fields", ff_)):
    o = {}
    for thr in (0.05, 0.1):
        A, n = mats(fl, pts_f, "d_marg_R", [0, thr, INF]); o[f"floor_within{thr}R_of_wall_or_peak"] = split_stats(A, n, np.array([1, 0], bool), ALL, W_ALL)
        o[f"floor_within{thr}R_of_wall_or_peak_radius_stratified"] = strat_split(fl, pts_f, "d_marg_R", thr, ALL, W_ALL)
    o["floor_within3km_of_wall_or_peak_radius_stratified"] = strat_split(fl, pts_f, "d_marg", 3, ALL, W_ALL)
    A, n = mats(fl.assign(r=fl.rho), pts_f.assign(r=pts_f.rho), "r", [0, 0.2, INF]); o["floor_within0.2R_of_center_vs_rest"] = split_stats(A, n, np.array([1, 0], bool), ALL, W_ALL)
    G9X[lab] = o
REC["g9_additions"] = G9X
for lab in ("hollows", "fields"):
    for k, s in G9X[lab].items():
        print("G9", lab, k, {kk: (round(v[0], 3), round(v[1], 3), round(v[2], 3)) if isinstance(v, list) and len(v) == 3 and v[0] is not None and v[1] is not None else v for kk, v in s.items() if kk in ("n_points", "n_in", "n_craters_with_points", "SDR", "DR")})
# ---- g11 additions (6 Nov 2026): eta_b tests; floor near the wall alone and near the peak alone, radius-scaled ----
G11X = {"note": "eta_b = (z - z_base)/(z_peak95 - z_base), z_base = 5th percentile of the crater's valid peak-cell elevations; bands as for eta (< 1/3 low, >= 2/3 high; "
                "eta_b < 0 counts as low). Near-wall / near-peak: floor cells and floor hollows/fields with d_wall/R or d_peak/R below the threshold vs the rest of "
                "the floor; floor in craters without wall (peak) cells has infinite distance and drops out of that test, as in g9.", "RB": RB}
if "eta_b" in pk:
    E = {}
    for lab, pts in (("hollows", hp_), ("fields", fp_)):
        A, n = mats(pk, pts, "eta_b", E_ETA)
        E[lab] = dict(n_points_finite_eta_b=int(np.isfinite(pts.eta_b).sum()), n_points=len(pts), craters_with_points=int(pts[np.isfinite(pts.eta_b)].uid.nunique()),
                      peak_eta_b_low_vs_upper=split_stats(A, n, np.array([1, 0, 0], bool), ALL, W_ALL),
                      peak_eta_b_high_vs_lower=split_stats(A, n, np.array([0, 0, 1], bool), ALL, W_ALL), peak_eta_b_terciles=tercile_fracs(A, n, ALL))
    E["sample"] = dict(peak_hollows_with_eta_b_lt_1_3=int((hp_.eta_b < 1/3).sum()), craters_with_low_eta_b_area=int(pk[pk.eta_b < 1/3].uid.nunique()),
                       median_zbase_minus_zfloor_m_hollow_hosts=float(cr[cr.uid.isin(hp_.uid)].z5_rel.median()) if "z5_rel" in cr else None,
                       median_peak_relief_above_base_m_hollow_hosts=float((cr.z95_rel - cr.z5_rel)[cr.uid.isin(hp_.uid)].median()) if "z5_rel" in cr else None,
                       peak_hollows_eta_lt_1_3=int((hp_.eta < 1/3).sum()), peak_cells_frac_area_eta_lt_1_3=float(pk.a[pk.eta < 1/3].sum()/pk.a.sum()),
                       peak_cells_frac_area_eta_b_lt_1_3=float(pk.a[pk.eta_b < 1/3].sum()/pk.a.sum()))
    BRB = {}
    for nm, s_ in [("central_peak", cr.morph == "central_peak"), ("peak_ring", cr.morph == "peak_ring"), ("unclassified_center", cr.morph == "unclassified")]:
        s_ = s_.values; W_ = boot_weights(s_); BRB[nm] = {}
        for lab, pts in (("hollows", hp_), ("fields", fp_)):
            A, n = mats(pk, pts, "eta_b", E_ETA)
            BRB[nm][lab] = dict(peak_eta_b_low_vs_upper=split_stats(A, n, np.array([1, 0, 0], bool), s_, W_), peak_eta_b_high_vs_lower=split_stats(A, n, np.array([0, 0, 1], bool), s_, W_))
    E["breakdowns"] = BRB
    G11X["eta_b"] = E
for lab, pts_f in (("hollows", hf_), ("fields", ff_)):
    o = {}
    for col, nm in (("d_wall_R", "wall"), ("d_peak_R", "peak")):
        for thr in (0.05, 0.1):
            A, n = mats(fl, pts_f, col, [0, thr, INF]); o[f"floor_within{thr}R_of_{nm}"] = split_stats(A, n, np.array([1, 0], bool), ALL, W_ALL)
            o[f"floor_within{thr}R_of_{nm}_radius_stratified"] = strat_split(fl, pts_f, col, thr, ALL, W_ALL)
    for col, nm in (("d_wall", "wall"), ("d_peak", "peak")):
        o[f"floor_within3km_of_{nm}_radius_stratified"] = strat_split(fl, pts_f, col, 3, ALL, W_ALL)
    G11X[lab] = o
REC["g11_additions"] = G11X
json.dump(REC, open(f"{HERE}/positions_record.json", "w"), indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else o)
ho.drop(columns=["ci"]).to_csv(f"{HERE}/hollow_positions.csv", index=False)
def show(R, k): s = R[k]; return f"{k}: n={s['n_points']} in={s['n_in']} f_h={s['frac_points_in']:.3f} f_A={s['frac_area_in_pooled']:.3f} f_E={s['frac_expected_in_crater_matched']:.3f} DR={s['DR'][0]:.2f} ({s['DR'][1]:.2f}-{s['DR'][2]:.2f}) SDR={s['SDR'][0]:.2f} ({s['SDR'][1]:.2f}-{s['SDR'][2]:.2f})"
for lab, R in [("HOLLOWS", REC["hollows"]), ("FIELDS", REC["fields"]), ("EXCL0.7", REC["sensitivity_exclude_within_0.7km_of_boundary"])] + [(f"{k}", v["hollows"]) for k, v in TH.items()]:
    print("==", lab)
    for k in R:
        if "terciles" in k: print("  ", k, {kk: [round(x, 3) for x in vv] if isinstance(vv, list) and vv and isinstance(vv[0], float) else vv for kk, vv in R[k].items()})
        else: print("  ", show(R, k))
print(json.dumps(REC["sample"])); print("mean q", REC["mean_q_peak_hollows"], REC["mean_q_peak_fields"])
print(json.dumps(REC["sensitivity_position_jitter_0.5km"]))
