"""Automatic shadow-length measurement of one hollow detection in a raw MDIS NAC EDR frame.

Geometry (established on the frames in qc/explore_sun_az.png and by the dark/bright test in sun_convention.json):
the label SUB_SOLAR_AZIMUTH is measured clockwise from the +sample axis with lines increasing downward, and points TOWARD the
sun. Inside a depression the shadow therefore lies on the sunward side (low u below), cast by the sunward rim, and the
sun-facing (anti-sunward) wall is bright. For ~75% of MESSENGER NAC frames the azimuth is 179.2-179.5 deg, so the sun is
toward -sample (left) and shadows run along image rows.

Method (measure()):
1. Crop the Bickel box plus a margin (max(0.5 x box size, 12 px), clipped to the frame).
2. Resample the crop on a sun-aligned grid at 0.25 px: u along the direction the light travels (anti-sunward), w across it
   (bilinear, scipy.ndimage.map_coordinates).
3. Levels: B = median DN of the margin ring outside the box (the lit surroundings); S = 2nd percentile DN inside the box (the
   deepest shadow). Threshold T = S + FRAC*(B - S) (FRAC = 0.5: the half-contrast edge, the edge a person picks on the image).
   No shadow if B - S < MIN_CONTRAST*(B - dark), dark = label DARK_STRIP_MEAN.
4. Along every w-line, runs of DN < T. A run counts when it starts inside the box (within START_TOL px of its sunward edge),
   starts after a lit sample (not at the crop edge) and ends inside the box (within END_TOL px of its anti-sunward edge; a dark
   run that continues past the box belongs to something else, e.g. a crater-wall shadow or dark terrain, and is rejected). The per-line maximum run length is smoothed with
   a running median over 1.25 px across w (5 sub-lines) to suppress single-line noise; L_px = maximum of the smoothed profile.
5. Ground length L = L_px * px * f, with f = sqrt((cos D / cos e)^2 + sin^2 D) the foreshortening along the sun direction
   (D = sun azimuth - sub-spacecraft azimuth, e = emission); depth h = L / tan(i), i = incidence at the detection (sphere,
   label sub-solar point, camera-model position).
Iteration 2 (v2, after the first validation pass and the contact sheets qc/val_*_v1.png): (a) lines in the outer EDGE_FRAC
of the box across the sun direction are ignored, because the longest v1 runs often lay along the box's top or bottom edge on
dark terrain or a neighbouring crater-wall shadow; (b) a run must be followed within TAIL_PX by a surface at least as bright as
the surroundings B (the lit, sun-facing far wall of a depression); v1 only flagged this. Both are depression-geometry rules,
not fitted to the manual depths.
Flags: no_shadow, edge_sun (dark region continues past the sunward crop edge, i.e. the run is truncated or belongs to a larger
shadow), edge_far (run reaches the anti-sunward crop edge), multi (>= 2 separate shadow runs >= 1.5 px in the box on the same
or different lines, i.e. several hollows or a compound hollow), not_depression (the 1.5 px after the shadow tip are not brighter
than T, i.e. no lit far wall/floor), short (L_px < 2, at the resolution limit).
"""
import numpy as np
from scipy.ndimage import map_coordinates, median_filter, label as cc_label

STEP = 0.25
FRAC = 0.5
MIN_CONTRAST = 0.12
START_TOL = 2.0
END_TOL = 2.0
EDGE_FRAC = 0.10      # v2: ignore the outer 10% of the box across the sun direction (shadows of the box margins)
TAIL_PX = 3.0         # v2: a run must be followed within TAIL_PX by a surface at least as bright as B (lit far wall)
REQUIRE_TAIL = True
STAT = 'max'          # 'max' of the smoothed per-line lengths, or 'p90' (sensitivity)


def runs(mask):
    """Start/end indices (end exclusive) of True runs in a 1-D boolean array."""
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return np.where(d == 1)[0], np.where(d == -1)[0]


def measure(img, info, box, frac=FRAC, return_debug=False):
    """box = (x0, x1, y0, y1) in full-frame pixels (sample, line). Returns a dict."""
    H, W = img.shape
    x0, x1, y0, y1 = box
    bw, bh = x1 - x0, y1 - y0
    m = max(0.5 * max(bw, bh), 12.0)
    az = np.radians(info["subsol_az"])
    du = np.array([-np.cos(az), -np.sin(az)])        # light-travel direction in (sample, line)
    dw = np.array([-du[1], du[0]])
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    corners = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]) - [cx, cy]
    cu, cw = corners @ du, corners @ dw
    ub0, ub1, wb0, wb1 = cu.min(), cu.max(), cw.min(), cw.max()
    U = np.arange(ub0 - m, ub1 + m + 1e-6, STEP); Wv = np.arange(wb0 - m, wb1 + m + 1e-6, STEP)
    UU, WW = np.meshgrid(U, Wv)
    xs = cx + UU * du[0] + WW * dw[0]; ys = cy + UU * du[1] + WW * dw[1]
    valid = (xs >= 0) & (xs <= W - 1) & (ys >= 0) & (ys <= H - 1)
    A = map_coordinates(img, [ys, xs], order=1, mode="nearest")
    inbox = (UU >= ub0) & (UU <= ub1) & (WW >= wb0) & (WW <= wb1)
    ring = valid & ~inbox
    out = dict(B=np.nan, S=np.nan, T=np.nan, L_px=np.nan, flags="", n_runs=0, run_w=np.nan, run_u0=np.nan, run_u1=np.nan,
               along_px=float(ub1 - ub0), across_px=float(wb1 - wb0), frac_valid=float(valid.mean()))
    if ring.sum() < 20 or (valid & inbox).sum() < 20:
        out["flags"] = "frame_edge"; return out
    B = float(np.median(A[ring])); S = float(np.percentile(A[valid & inbox], 2)); dark = info.get("dark", 33.5)
    if not np.isfinite(dark): dark = 33.5
    T = S + frac * (B - S); out.update(B=B, S=S, T=T)
    flags = []
    if B - S < MIN_CONTRAST * max(B - dark, 1):
        out["flags"] = "no_shadow"; out["L_px"] = 0.0; return out
    sh = (A < T) & valid
    nW = len(Wv); best = np.zeros(nW); bu0 = np.full(nW, np.nan); bu1 = np.full(nW, np.nan)
    e_sun = np.zeros(nW, bool); e_far = np.zeros(nW, bool); e_over = np.zeros(nW, bool); allruns = []
    wi0 = wb0 + EDGE_FRAC * (wb1 - wb0); wi1 = wb1 - EDGE_FRAC * (wb1 - wb0)
    if wi1 - wi0 < 1.0: wi0, wi1 = (wb0 + wb1) / 2 - 0.5, (wb0 + wb1) / 2 + 0.5
    n_tailfail = 0
    for j in range(nW):
        if Wv[j] < wi0 or Wv[j] > wi1: continue
        s, e = runs(sh[j])
        for a, b in zip(s, e):
            ua, ub = U[a], U[b - 1] + STEP
            if ua < ub0 - START_TOL or ua > ub1:
                if a == 0 and U[b - 1] >= ub0: e_sun[j] = True
                continue
            if ub > ub1 + END_TOL:
                e_over[j] = True; continue
            if a == 0 or not valid[j, a - 1]: e_sun[j] = True; continue
            if b >= len(U) or not valid[j, min(b, len(U) - 1)]: e_far[j] = True; continue
            L = (b - a) * STEP
            if REQUIRE_TAIL:
                tail = A[j, b:b + int(TAIL_PX / STEP) + 1]
                if len(tail) == 0 or tail.max() < B:
                    if L >= 1.5: n_tailfail += 1
                    continue
            if L >= 1.5: allruns.append((j, ua, ub))
            if L > best[j]: best[j], bu0[j], bu1[j] = L, ua, ub
    sm = median_filter(best, size=5, mode="nearest")
    jbest = int(np.argmax(sm)); Lpx = float(sm[jbest])
    if STAT == "p90" and (sm > 0).sum() >= 3:
        Lpx = float(np.percentile(sm[sm > 0], 90)); jbest = int(np.argmin(np.abs(sm - Lpx)))
    out["n_tailfail"] = n_tailfail
    out["L_px"] = Lpx
    if Lpx > 0:
        # use the line whose raw run equals the smoothed maximum nearest to jbest
        js = np.arange(max(0, jbest - 2), min(nW, jbest + 3)); jj = js[np.argmin(np.abs(best[js] - Lpx))]
        out.update(run_w=float(Wv[jj]), run_u0=float(bu0[jj]), run_u1=float(bu1[jj]))
        # depression test: just beyond the shadow tip the surface is lit (above T)
        k1 = int(round((bu1[jj] - U[0]) / STEP)); tail = A[jj, k1:k1 + int(TAIL_PX / STEP) + 1]
        if len(tail) == 0 or np.max(tail) < B: flags.append("not_depression")
        if Lpx < 2: flags.append("short")
    else:
        flags.append("no_shadow")
    # independent shadow runs (components) inside the box
    lab, nlab = cc_label(sh & inbox)
    sizes = []
    for i in range(1, nlab + 1):
        jj_, kk_ = np.where(lab == i)
        if len(jj_) == 0: continue
        ext_u = (kk_.max() - kk_.min() + 1) * STEP; ext_w = (jj_.max() - jj_.min() + 1) * STEP
        if ext_u >= 1.5 and ext_w >= 1.0: sizes.append(ext_u)
    out["n_runs"] = len(sizes)
    if len(sizes) >= 2: flags.append("multi")
    if e_sun[max(0, jbest - 2):jbest + 3].any() or (e_sun[(Wv >= wi0) & (Wv <= wi1)].mean() > 0.3): flags.append("edge_sun")
    if e_far[max(0, jbest - 2):jbest + 3].any(): flags.append("edge_far")
    out["flags"] = ";".join(flags)
    if return_debug:
        out["_dbg"] = dict(A=A, U=U, W=Wv, sh=sh, inbox=inbox, du=du, dw=dw, cx=cx, cy=cy, best=best, sm=sm)
    return out


def foreshortening(info):
    D = np.radians(info["subsol_az"] - info["subsc_az"]); e = np.radians(info["emi"])
    return float(np.sqrt((np.cos(D) / np.cos(e)) ** 2 + np.sin(D) ** 2))
