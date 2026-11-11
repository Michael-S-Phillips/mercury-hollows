# g6_remaining copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/color_mdr/extract_tiles.py (6 Nov 2026). Changed: hollow-mask input (G6_MASK) and tile output folder (G6_TILES) only.
"""Per-tile extraction of MDIS MDR 8-color spectra for craters and measured hollows (resumable).

Usage: python3 extract_tiles.py [time_budget_s]
Processes tiles that lack tiles/<TILE>.done, stopping before the time budget is spent.
Per tile it writes tiles_qc/<TILE>_complete.npz (features whose whole outer circle lies in this tile's
exclusive region: medians computed here) and tiles_qc/<TILE>_partial.npz (pixels of features that straddle
tile boundaries; merged later by merge_spectra.py).

Regions
  crater: inner  = d <= 0.9 R          (set1, hollow-masked)
          annulus = 1.2 R < d <= 2.0 R (set2, hollow-masked)
  hollow: d <= 10 km                   (set1 hollow-masked, set2 unmasked)
Hollow mask: pixel centers within r_excl of any cataloged hollow, r_excl = 1 km + min(D_eq, 5 km)/2
for hollows_canonical.csv (6,731 hollows, corrected positions) and 1 km for the 2,518 depth-measured
hollows (point positions). Large features are subsampled with a fixed per-feature stride so that the
outer radius spans <= 200 pixels; n_est = n_sampled * stride^2.
"""
import os, sys, time, glob
import numpy as np, pandas as pd
import mdr_georef as g

T0 = time.time()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 150.0
OUT = os.environ.get('G6_TILES', 'tiles_qc')
BASE = os.path.expanduser('~/mnt/hollows-emissivity')
CF = 'manuscripts/PlanetaryResearch_hollows/analysis_workshop/summary_model/crater_features.csv'
HF = 'manuscripts/insolation_thermal/data/test1_rows_inc.csv'
CAN = os.environ.get('G6_MASK', 'code/rerun/out_v2/hollows_canonical.csv')
SCALE0 = 0.665304   # km/px nominal
R_HOLLOW = 10.0


def build_features():
    cr = pd.read_csv(os.path.join(BASE, CF), usecols=['id', 'lat', 'lon', 'D'])
    ho = pd.read_csv(os.path.join(BASE, HF), usecols=['latitude', 'longitude'])
    R = cr.D.values / 2
    f = pd.DataFrame(dict(kind=np.r_[np.zeros(len(cr), int), np.ones(len(ho), int)],
                          fid=np.r_[np.arange(len(cr)), np.arange(len(ho))],
                          lat=np.r_[cr.lat.values, ho.latitude.values],
                          lon=np.r_[cr.lon.values % 360, ho.longitude.values % 360],
                          r1=np.r_[0.9 * R, np.full(len(ho), R_HOLLOW)],
                          r2a=np.r_[1.2 * R, np.full(len(ho), np.inf)],
                          rout=np.r_[2.0 * R, np.full(len(ho), R_HOLLOW)]))
    f['stride'] = np.maximum(1, np.ceil(f.rout / SCALE0 / 200.0)).astype(int)
    return f


def build_exclusion():
    can = pd.read_csv(os.path.join(BASE, CAN), usecols=['lat', 'lon', 'D_eq_km'])
    ho = pd.read_csv(os.path.join(BASE, HF), usecols=['latitude', 'longitude'])
    lat = np.r_[can.lat.values, ho.latitude.values]
    lon = np.r_[can.lon.values % 360, ho.longitude.values % 360]
    r = np.r_[1.0 + 0.5 * np.minimum(can.D_eq_km.fillna(0).values, 5.0), np.full(len(ho), 1.0)]
    return lat, lon, r


def cap_lon_halfwidth(lat0, dang_deg):
    """Half-width in longitude (deg) of a spherical cap of angular radius dang centered at lat0."""
    if abs(lat0) + dang_deg >= 89.9:
        return 180.0
    return float(np.degrees(np.arcsin(min(1.0, np.sin(np.radians(dang_deg)) / np.cos(np.radians(lat0))))))


def window(t, lat0, lon0, rout, stride=1, margin_px=2):
    """Return (ii, jj, complete) index vectors (strided) of the window covering the cap, clipped to the array."""
    dang = np.degrees(rout / g.R_KM)
    i0, j0 = t.ll2pix(lat0, lon0)
    i0 = float(i0); j0 = float(j0)
    if not t.polar:
        hi = rout / t.scale_km + margin_px
        dlon = cap_lon_halfwidth(lat0, dang)
        hj = (dlon * np.pi / 180 * t.R * t.coslat_c / t.scale_km + margin_px) if dlon < 180 else 1e9
        a0, a1 = int(np.floor(i0 - hi)), int(np.ceil(i0 + hi)) + 1
        b0, b1 = int(np.floor(j0 - hj)), int(np.ceil(j0 + hj)) + 1
        eps = 0.02
        dl_w = g.wrap180(lon0 - dlon - t.lon_w)
        complete = (dlon < 90 and lat0 + dang < t.lat_max - eps and lat0 - dang > t.lat_min + eps
                    and abs(lat0) + dang < g.POLAR_LIM - eps and dl_w > eps and dl_w + 2 * dlon < t.lon_span - eps)
    else:
        sgn = 1 if t.north else -1
        if sgn * lat0 + dang <= g.POLAR_LIM:      # cap never reaches |lat| > 65
            return None, None, False
        kmax = 2.0 / (1.0 + np.sin(np.radians(max(50.0, sgn * lat0 - dang))))
        h = rout * kmax / t.scale_km + margin_px
        a0, a1 = int(np.floor(i0 - h)), int(np.ceil(i0 + h)) + 1
        b0, b1 = int(np.floor(j0 - h)), int(np.ceil(j0 + h)) + 1
        complete = (sgn * lat0 - dang > g.POLAR_LIM + 0.02) and a0 >= 0 and b0 >= 0 and a1 <= t.H and b1 <= t.W
    a0c, a1c, b0c, b1c = max(a0, 0), min(a1, t.H), max(b0, 0), min(b1, t.W)
    if a0c >= a1c or b0c >= b1c:
        return None, None, False
    if not t.polar:
        complete = complete and a0 >= 0 and b0 >= 0 and a1 <= t.H and b1 <= t.W
    return np.arange(a0c, a1c, stride), np.arange(b0c, b1c, stride), complete


def stats(refl, cnt, stride):
    n = refl.shape[1]
    if n == 0:
        return n, 0, np.full(8, np.nan, np.float32), np.nan
    return n, n * stride * stride, np.median(refl, axis=1).astype(np.float32), float(np.median(cnt))


def process_tile(name, F, EX):
    t = g.Tile(name)
    other = g.Tile('2700_H15SP4') if name == 'H15SP4' else None
    refl, cnt, filled = g.load_refl(t, other)
    valid = np.isfinite(cnt) & np.isfinite(refl[0])
    if t.polar:
        I, J = np.mgrid[0:t.H, 0:t.W]
        LAT, LON = t.pix2ll(I, J); del I, J
        usable = valid & t.assigned(LAT, LON)
    else:
        lat1, _ = t.pix2ll(np.arange(t.H), np.zeros(t.H))
        _, lon1 = t.pix2ll(np.zeros(t.W), np.arange(t.W))
        usable = valid & t.assigned(lat1[:, None], lon1[None, :])
    def latlon(ii, jj):
        if t.polar:
            return LAT[np.ix_(ii, jj)], LON[np.ix_(ii, jj)]
        return np.broadcast_to(lat1[ii][:, None], (len(ii), len(jj))), np.broadcast_to(lon1[jj][None, :], (len(ii), len(jj)))

    # ---- hollow mask
    hm = np.zeros((t.H, t.W), bool)
    elat, elon, er = EX
    ei, ej = t.ll2pix(elat, elon)
    near = (ei > -20) & (ei < t.H + 20) & (ej > -60) & (ej < t.W + 60)
    if t.polar:
        near &= (elat > 55) if t.north else (elat < -55)
    for k in np.flatnonzero(near):
        ii, jj, _ = window(t, elat[k], elon[k], er[k], 1, margin_px=2) if not t.polar else \
            (np.arange(max(0, int(ei[k]) - 8), min(t.H, int(ei[k]) + 9)), np.arange(max(0, int(ej[k]) - 8), min(t.W, int(ej[k]) + 9)), False)
        if ii is None or len(ii) == 0 or len(jj) == 0:
            continue
        la, lo = latlon(ii, jj)
        d = g.gc_dist_km(elat[k], elon[k], la, lo)
        sub = hm[np.ix_(ii, jj)]
        hm[np.ix_(ii, jj)] = sub | (d <= er[k])

    # ---- features
    fi, fj = t.ll2pix(F.lat.values, F.lon.values)
    rpx = F.rout.values / t.scale_km * (1.2 if t.polar else 1.0)
    if t.polar:
        sgn = 1 if t.north else -1
        cand = (sgn * F.lat.values + np.degrees(F.rout.values / g.R_KM) > g.POLAR_LIM)
    else:
        cand = np.ones(len(F), bool)
    cand &= (fi + rpx > -2) & (fi - rpx < t.H + 2)
    comp_rows, part = [], dict(fidx=[], region=[], hm=[], refl=[], cnt=[])
    for k in np.flatnonzero(cand):
        r = F.iloc[k]
        ii, jj, complete = window(t, r.lat, r.lon, r.rout, int(r.stride))
        if ii is None or len(ii) == 0 or len(jj) == 0:
            continue
        la, lo = latlon(ii, jj)
        d = g.gc_dist_km(r.lat, r.lon, la, lo)
        u = usable[np.ix_(ii, jj)]
        reg = np.zeros(d.shape, np.int8)
        reg[d <= r.r1] = 1
        reg[(d > r.r2a) & (d <= r.rout)] = 2
        sel = u & (reg > 0)
        if not complete and not sel.any():
            continue
        ix = np.ix_(ii, jj)
        R8 = refl[:, ix[0], ix[1]][:, sel]
        C = cnt[ix][sel]; H_ = hm[ix][sel]; G_ = reg[sel]
        if complete:
            st = int(r.stride)
            if r.kind == 0:
                s1 = (G_ == 1) & ~H_; s2 = (G_ == 2) & ~H_
            else:
                s1 = (G_ == 1) & ~H_; s2 = (G_ == 1)
            n1, n1e, m1, c1 = stats(R8[:, s1], C[s1], st)
            n2, n2e, m2, c2 = stats(R8[:, s2], C[s2], st)
            comp_rows.append((k, n1, n1e, c1, n2, n2e, c2, int((G_ == 1).sum()), int(((G_ == 1) & H_).sum()), *m1, *m2))
        else:
            part['fidx'].append(np.full(sel.sum(), k, np.int32)); part['region'].append(G_)
            part['hm'].append(H_); part['refl'].append(R8.T.astype(np.float32)); part['cnt'].append(C.astype(np.float32))
    cols = ['fidx', 'n1', 'n1_est', 'cnt1', 'n2', 'n2_est', 'cnt2', 'n_reg1_all', 'n_reg1_masked'] + \
           ['s1_%d' % w for w in g.WAVES] + ['s2_%d' % w for w in g.WAVES]
    comp = np.array(comp_rows, dtype=np.float64).reshape(-1, len(cols))
    np.savez_compressed(os.path.join(OUT, name + '_complete.npz'), data=comp, cols=np.array(cols))
    if part['fidx']:
        np.savez_compressed(os.path.join(OUT, name + '_partial.npz'),
                            fidx=np.concatenate(part['fidx']), region=np.concatenate(part['region']),
                            hm=np.concatenate(part['hm']), refl=np.concatenate(part['refl']), cnt=np.concatenate(part['cnt']))
    info = dict(tile=name, n_complete=len(comp_rows), n_partial_feat=len(part['fidx']),
                n_partial_px=int(sum(len(a) for a in part['fidx'])), frac_valid=float(valid.mean()),
                frac_usable=float(usable.mean()), frac_hm=float(hm.mean()), n_filled_from_2700=int(filled.sum()), n_qc_rejected=g.LAST_QC.get('n_qc_rejected'))
    open(os.path.join(OUT, name + '.done'), 'w').write(repr(info))
    return info


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    F = build_features(); EX = build_exclusion()
    F.to_csv(os.path.join(OUT, 'features.csv'), index=False)
    todo = [n for n in g.all_tiles() if not os.path.exists(os.path.join(OUT, n + '.done'))]
    todo = sorted(todo, key=lambda n: ('P' in n[3:5]))  # polar tiles last
    for n in todo:
        est = 60 if n[3:5] in ('NP', 'SP') else 15
        if time.time() - T0 + est > BUDGET:
            break
        t1 = time.time()
        info = process_tile(n, F, EX)
        print('%-7s %5.1fs %s' % (n, time.time() - t1, info), flush=True)
    left = [n for n in g.all_tiles() if not os.path.exists(os.path.join(OUT, n + '.done'))]
    print('remaining tiles:', len(left))
