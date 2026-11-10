"""Curvature-based terrain classes, version 3 (g10, 6 Nov 2026), from curv_classes.py (version 2, used for g11) with two fixes:
(1) ring radii are filled only across gaps of at most 2 sectors (30 deg); sectors in longer gaps have no ring;
(2) ring radii smoothed with a 3-sector circular median and ring presence interpolated between sector centers;
(3) peak and ring cells must stand at or above the elevation of their base knee (minus 25 m), so low floor between massifs and
pockets below the base are floor.
Version 2 description: No relief threshold except a 100 m minimum for declaring a
resolved peak. In each of 24 azimuth sectors: median elevation in radial bins of max(0.025 R, 1 pixel), filled across gaps,
Gaussian-smoothed with sigma max(sig_frac R, 1.5 pixels), sig_frac = 0.03 by default; curvature k = d2z/dr2 (concave up > 0).
- Rim: highest point of the profile at 0.8-1.15 R.
- Central-peak craters: summit = highest point within 0.45 R; floor low = lowest point between the summit and 0.8 R;
  peak base = curvature maximum (concave-up knee) between the summit and the floor low; a sector has a peak only if the
  summit stands >= 100 m above the base. Wall foot = curvature maximum between the floor low and the rim.
- Peak-ring craters (pr, rpc, pb), version 2: the crest is the highest local maximum of the profile at 0.25-0.75 R that stands
  >= 100 m above the lowest point inside it. Inner knee = curvature maximum between 0.05 R and the crest. Moat = lowest
  point between the crest and the rim; if the crest stands >= 50 m above it, outer knee = curvature maximum between the crest
  and the moat, otherwise the ring is taken as symmetric about its crest (outer edge = crest + (crest - inner knee)). Wall
  foot = curvature maximum between max(moat, outer edge) and the rim. Sectors without a crest take ring radii interpolated
  from neighbors if at least a quarter of the sectors have one; otherwise the crater has no ring class. The central-peak
  rule is also applied (rpc and pb have central peaks).
- All per-sector radii are smoothed with a 3-sector circular median before interpolation in azimuth.
Radii are interpolated linearly in azimuth between sector centers. peak = within the peak radii; wall = wall foot to 1.05 R;
floor = the rest within 1.0 R."""
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d
def _az(vals, az, nsec):
    w = 360.0 / nsec; cen = (np.arange(nsec) + 0.5) * w
    return np.interp(az, np.concatenate([cen - 360, cen, cen + 360]), np.tile(vals, 3))
def profiles(z, rho, az, R, nsec=24, sig_frac=0.03):
    w = max(0.025, 0.665 / R); edges = np.arange(0, 1.2 + 1e-9, w); rc = edges[:-1] + w / 2
    ok = np.isfinite(z) & (rho < 1.2); sec = np.minimum((az // (360 / nsec)).astype(int), nsec - 1)
    df = pd.DataFrame({"s": sec[ok], "b": np.digitize(rho[ok], edges) - 1, "z": z[ok]}); med = df.groupby(["s", "b"]).z.median()
    P = np.array([[med.get((s, i), np.nan) for i in range(len(rc))] for s in range(nsec)])
    for s in range(nsec):
        g = np.isfinite(P[s])
        if g.sum() > 3: P[s] = np.interp(rc, rc[g], P[s][g])
    Ps = gaussian_filter1d(P, max(sig_frac, 1.5 * 0.665 / R) / w, axis=1, mode="nearest")
    dr = w * R * 1000; K = np.gradient(np.gradient(Ps, dr, axis=1), dr, axis=1)
    return rc, P, Ps, K
def _argmax(a, idx): return idx[np.nanargmax(a[idx])] if len(idx) else None
def _argmin(a, idx): return idx[np.nanargmin(a[idx])] if len(idx) else None
def _med3(v):
    return np.array([np.median([v[(i - 1) % len(v)], v[i], v[(i + 1) % len(v)]]) for i in range(len(v))])
def _fill(v, ok):
    n = len(v); idx = np.arange(n); g = idx[ok]
    if len(g) == 0: return v
    return np.interp(idx, np.concatenate([g - n, g, g + n]), np.tile(v[ok], 3))
def _fill_short(v, ok, maxgap=2):
    """Fill NaN sectors by circular interpolation only where the run of missing sectors is <= maxgap; return values and validity."""
    n = len(v); ok = ok.copy(); out = v.copy()
    if ok.sum() == 0: return out, ok
    filled = _fill(v, ok); miss = ~ok; i = 0
    start = int(np.argmax(ok))   # walk the circle from a valid sector
    run = []
    for k in range(1, n + 1):
        j = (start + k) % n
        if miss[j]: run.append(j)
        else:
            if 0 < len(run) <= maxgap:
                for q in run: out[q] = filled[q]; ok[q] = True
            run = []
    return out, ok
def _sector(az, nsec): return np.minimum((az // (360.0 / nsec)).astype(int), nsec - 1)
def classify_curv(z, rho, az, R, cent, nsec=24, hmin=100.0, sig_frac=0.03, tol=25.0):
    rc, P, Ps, K = profiles(z, rho, az, R, nsec, sig_frac); ring = cent in ("pr", "rpc", "pb"); I = np.arange(len(rc))
    rpk = np.zeros(nsec); zpk = np.full(nsec, np.nan); rwall = np.full(nsec, np.nan)
    rin = np.full(nsec, np.nan); rout = np.full(nsec, np.nan); zr = np.full(nsec, np.nan); info = []
    for s in range(nsec):
        p, k = Ps[s], K[s]; irim = _argmax(p, I[(rc >= 0.8) & (rc <= 1.15)]); rec = dict(s=s)
        it = _argmax(p, I[rc <= 0.45]); il = _argmin(p, I[(rc >= rc[it]) & (rc <= 0.8)])
        ib = _argmax(k, I[(rc > rc[it]) & (rc <= rc[il])]) if il is not None and il > it else None
        if ib is not None and p[it] - p[ib] >= hmin and (not ring or rc[it] <= 0.2): rpk[s] = rc[ib]; zpk[s] = p[ib]
        rec["peak"] = float(rpk[s]); lowr = rc[il] if il is not None else 0.5
        if ring:
            loc = [i for i in I[(rc >= 0.25) & (rc <= 0.75)] if 0 < i < len(rc) - 1 and p[i] >= p[i - 1] and p[i] >= p[i + 1]]
            loc = [i for i in loc if p[i] - np.min(p[:i + 1]) >= hmin]
            if loc:
                ic = max(loc, key=lambda i: p[i]); ik = _argmax(k, I[(rc >= 0.05) & (rc < rc[ic])])
                im = _argmin(p, I[(rc > rc[ic]) & (rc < rc[irim])])
                if im is not None and p[ic] - p[im] >= 50: io = _argmax(k, I[(rc > rc[ic]) & (rc <= rc[im])]); ro = rc[io]; zo = p[io]
                else:
                    ro = min(rc[ic] + (rc[ic] - (rc[ik] if ik is not None else rc[ic] - 0.1)), 0.9); im = None
                    zo = np.interp(ro, rc, p)
                rin[s] = rc[ik] if ik is not None else max(rc[ic] - 0.1, 0.05); rout[s] = ro
                zr[s] = min(p[ik] if ik is not None else np.interp(rin[s], rc, p), zo)
                lowr = max(rc[im] if im is not None else ro, ro); rec.update(crest=float(rc[ic]), inner=float(rin[s]), outer=float(rout[s]), zbase=float(zr[s]))
        iw = _argmax(k, I[(rc >= lowr) & (rc < rc[irim])])
        rwall[s] = rc[iw] if iw is not None else 0.75; rec["wall"] = float(rwall[s]); info.append(rec)
    rwall = _med3(rwall); RW = _az(rwall, az, nsec); sec = _sector(az, nsec)
    okp = rpk > 0; rpk_s = _med3(rpk)
    zpk_f = _fill(np.where(okp, zpk, 0.0), okp) if okp.any() else np.zeros(nsec)
    peak = (rho <= _az(rpk_s, az, nsec)) & (z >= _az(zpk_f, az, nsec) - tol)
    okr = np.isfinite(rin); nring = int(okr.sum())
    if ring and nring >= nsec / 4:
        rin2, ok2 = _fill_short(rin, okr); rout2, _ = _fill_short(rout, okr); zr2, _ = _fill_short(zr, okr)
        rin2 = _med3(_fill(rin2, ok2)); rout2 = _med3(_fill(rout2, ok2)); zr2 = _med3(_fill(zr2, ok2))
        has = _az(ok2.astype(float), az, nsec) >= 0.5   # ring present where the interpolated sector flag is >= 0.5
        peak |= has & (rho >= _az(rin2, az, nsec)) & (rho <= _az(rout2, az, nsec)) & (z >= _az(zr2, az, nsec) - tol)
        nring = int(ok2.sum())
    wall = (rho >= RW) & (rho <= 1.05) & ~peak
    flo = (rho <= 1.0) & ~peak & ~wall & np.isfinite(z)
    return dict(peak=peak & np.isfinite(z), wall=wall & np.isfinite(z), flo=flo, info=info, rc=rc, Ps=Ps, n_ring_sectors=nring)
