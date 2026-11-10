"""Terrain classes inside a host crater (g7, 6 Nov 2026): a fix of the class rules in peak_wall_test/peak_wall.py.

Old rules (peak_wall.py): floor level F0 = 20th percentile of elevation at 0.15-0.6 R; h = z - F0;
  peak = rho <= 0.45 (0.6 for peak rings) and h > H; wall = 0.55 (0.7) <= rho <= 1.05 and h > H; floor = rho <= 0.8 and h <= H.
  Two gaps follow: ground more than H above F0 between the peak and wall limits (0.45-0.55 R, or 0.6-0.7 R in peak-ring
  craters) is in no class, and so is ground within H of F0 beyond 0.8 R. In peak-ring craters the annular floor between the
  ring and the wall stands 200-600 m above the inner floor, so it falls in the gap or is counted as wall.
New rules: every pixel inside the rim gets one class. In each of 24 azimuth sectors the radial profile of h (median in
  0.025 R bins, 3-bin running mean) sets a boundary radius b between the peak side and the wall side.
  - Peak-ring craters (cent_str pr, rpc, pb): the ring crest is the highest point of the profile at 0.3-0.7 R, and the
    annular floor is the lowest point between the crest and 0.9 R. If the crest stands >= 50 m above that low point and the
    profile rises >= 100 m beyond it, b is the low point's radius, and beyond the crest heights are measured from the
    annular floor level instead of F0. Otherwise b = 0.65 R and heights are measured from F0.
    The low point is searched only out to 0.85 R. If fewer than a quarter of the sectors show such a low point, the crater
    is treated as having no resolved annular floor; otherwise sectors without one take values interpolated from their
    neighbors. b, the crest radius and the annular floor level are interpolated linearly in azimuth between sector centers.
  - Other craters: b is the lowest point of the profile at 0.35-0.65 R, clipped to [RPK, 0.55] (the old gap), and also
    interpolated in azimuth.
  - peak = rho <= b and more than H above the local floor level (and above F0); wall = b < rho <= 1.05 and more than H
    above the local floor level; floor = everything else within 1.0 R.
"""
import numpy as np, pandas as pd

def _circ_fill(v, ok):
    """Fill v where ~ok by circular linear interpolation from the ok entries (sector centers)."""
    n = len(v); idx = np.arange(n); good = idx[ok]
    if len(good) == 0: return v
    x = np.concatenate([good - n, good, good + n]); y = np.tile(v[ok], 3)
    return np.interp(idx, x, y)

def _az_interp(vals, az, nsec):
    """Value at each pixel azimuth, linear in azimuth between sector centers (circular)."""
    w = 360.0 / nsec; cen = (np.arange(nsec) + 0.5) * w
    x = np.concatenate([cen - 360, cen, cen + 360]); y = np.tile(vals, 3)
    return np.interp(az, x, y)

def classify(z, rho, az, ring, HT=250.0, RPK=0.45, nsec=24):
    fz = z[(rho > 0.15) & (rho < 0.6)]; fz = fz[np.isfinite(fz)]
    F0 = np.percentile(fz, 20); h = z - F0
    sec = np.minimum((az // (360.0 / nsec)).astype(int), nsec - 1)
    edges = np.arange(0.2, 0.95 + 1e-9, 0.025); rc = edges[:-1] + 0.0125
    ok = np.isfinite(h) & (rho >= 0.2) & (rho < 0.95)
    df = pd.DataFrame({"s": sec[ok], "b": np.digitize(rho[ok], edges) - 1, "h": h[ok]})
    med = df.groupby(["s", "b"]).h.median()
    B = np.full(nsec, np.nan); CR = np.full(nsec, np.nan); ZA = np.full(nsec, np.nan); info = []
    for s in range(nsec):
        p = np.array([med.get((s, i), np.nan) for i in range(len(rc))])
        ps = pd.Series(p).rolling(3, center=True, min_periods=1).mean().values
        rec = dict(sector=s, moat=False)
        if ring:
            j = np.where((rc >= 0.3) & (rc <= 0.7) & np.isfinite(ps))[0]
            if len(j):
                ic = j[np.nanargmax(ps[j])]; k = np.where((rc > rc[ic]) & (rc <= 0.85) & np.isfinite(ps))[0]
                if len(k):
                    im = k[np.nanargmin(ps[k])]; beyond = ps[(rc > rc[im]) & np.isfinite(ps)]
                    if ps[ic] - ps[im] >= 50 and len(beyond) and np.nanmax(beyond) - ps[im] >= 100:
                        B[s], CR[s], ZA[s] = rc[im], rc[ic], ps[im]
                        rec.update(moat=True, crest=float(rc[ic]), b=float(rc[im]), annular_floor_m=float(ps[im]))
        else:
            j = np.where((rc >= 0.35) & (rc <= 0.65) & np.isfinite(ps))[0]
            if len(j): B[s] = float(np.clip(rc[j[np.nanargmin(ps[j])]], RPK, 0.55)); rec.update(b=B[s])
        info.append(rec)
    nmoat = int(np.isfinite(ZA).sum())
    if ring and nmoat >= nsec / 4:
        m = np.isfinite(ZA); B = _circ_fill(B, m); CR = _circ_fill(CR, m); ZA = _circ_fill(ZA, m)
        bnd = _az_interp(B, az, nsec); crest = _az_interp(CR, az, nsec); za = _az_interp(ZA, az, nsec)
        Floc = np.where(rho > crest, F0 + za, F0)
    else:
        if ring: bnd = np.full(z.shape, 0.65)
        else: bnd = _az_interp(_circ_fill(B, np.isfinite(B)) if np.isfinite(B).any() else np.full(nsec, 0.5), az, nsec)
        Floc = np.full(z.shape, F0)
    hl = z - Floc
    peak = (rho <= bnd) & (hl > HT) & (h > HT)
    wall = (rho > bnd) & (rho <= 1.05) & (hl > HT)
    flo = (rho <= 1.0) & ~peak & ~wall & np.isfinite(z)
    return dict(peak=peak, wall=wall, flo=flo, F0=F0, h=h, info=info, n_moat=nmoat, annular=bool(ring and nmoat >= nsec / 4))
