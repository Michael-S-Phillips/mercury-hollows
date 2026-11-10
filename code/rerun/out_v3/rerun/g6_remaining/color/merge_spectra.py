# g6_remaining copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/color_mdr/merge_spectra.py (6 Nov 2026). Changed: nothing (reads and writes tiles_qc/ in the working directory, a link to _scratch/tiles_<cat>).
"""Merge per-tile outputs into per-feature median spectra; compute brightness-normalized spectra and indices.

Writes tiles_qc/merged_raw.csv (one row per feature in tiles/features.csv).
Brightness removal: nl_k = ln r_k - mean_k ln r_k (division by the 8-band geometric mean);
albedo_gm = geometric mean reflectance (saved, not a composition variable).
Indices: slope = ln(r750/r430); curv = ln(r1000/r750) - ln(r750/r560) (= ln r1000 - 2 ln r750 + ln r560);
curv_wl = ln r1000 minus the linear-in-wavelength extrapolation of ln r through 560 and 750 nm.
"""
import os, glob, numpy as np, pandas as pd
import mdr_georef as g

W = g.WAVES
F = pd.read_csv('tiles_qc/features.csv')
nF = len(F)
tiles = g.all_tiles()

# ---- complete features
cols = None; rows = []
for n in tiles:
    z = np.load('tiles_qc/%s_complete.npz' % n)
    cols = list(z['cols']); d = pd.DataFrame(z['data'], columns=cols); d['tile'] = n; rows.append(d)
comp = pd.concat(rows, ignore_index=True)
comp['fidx'] = comp.fidx.astype(int)
dup_complete = int(comp.fidx.duplicated().sum())
comp = comp.drop_duplicates('fidx')

# ---- partial features, processed in chunks of feature index to bound memory
NCH = 6
part_rows = []
part_fidx_all = set()
for ch in range(NCH):
    acc = {k: [] for k in ['fidx', 'region', 'hm', 'refl', 'cnt']}
    for n in tiles:
        p = 'tiles_qc/%s_partial.npz' % n
        if not os.path.exists(p):
            continue
        z = np.load(p)
        m = (z['fidx'] % NCH) == ch
        for k in acc:
            acc[k].append(z[k][m])
    a = {k: np.concatenate(v) for k, v in acc.items()}
    o = np.argsort(a['fidx'], kind='stable')
    a = {k: v[o] for k, v in a.items()}
    u, starts = np.unique(a['fidx'], return_index=True)
    ends = np.r_[starts[1:], len(a['fidx'])]
    for k, s, e in zip(u, starts, ends):
        part_fidx_all.add(int(k))
        kind = F.kind.iat[k]; st = int(F.stride.iat[k])
        reg = a['region'][s:e]; hm = a['hm'][s:e]; R8 = a['refl'][s:e]; C = a['cnt'][s:e]
        s1 = (reg == 1) & ~hm
        s2 = ((reg == 2) & ~hm) if kind == 0 else (reg == 1)
        row = dict(fidx=int(k), tile='multi', n_reg1_all=int((reg == 1).sum()), n_reg1_masked=int(((reg == 1) & hm).sum()))
        for tag, m in [('1', s1), ('2', s2)]:
            nn = int(m.sum())
            row['n' + tag] = nn; row['n%s_est' % tag] = nn * st * st
            row['cnt' + tag] = float(np.median(C[m])) if nn else np.nan
            med = np.median(R8[m], axis=0) if nn else np.full(8, np.nan)
            for w, v in zip(W, med):
                row['s%s_%d' % (tag, w)] = v
        part_rows.append(row)
part = pd.DataFrame(part_rows)
both = sorted(set(comp.fidx) & part_fidx_all)
allr = pd.concat([comp[~comp.fidx.isin(both)], part], ignore_index=True).set_index('fidx').reindex(range(nF))
out = pd.concat([F, allr], axis=1)


def derive(df, tag, pre):
    r = df[['s%s_%d' % (tag, w) for w in W]].values.astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        lr = np.log(r)
        gm = lr.mean(axis=1)
        o = pd.DataFrame(index=df.index)
        for i, w in enumerate(W):
            o['%sr%d' % (pre, w)] = r[:, i]
        o[pre + 'albedo_gm'] = np.exp(gm)
        for i, w in enumerate(W):
            o['%snl%d' % (pre, w)] = lr[:, i] - gm
        i430, i560, i750, i1000 = [list(W).index(x) for x in (430, 560, 750, 1000)]
        o[pre + 'slope'] = lr[:, i750] - lr[:, i430]
        o[pre + 'curv'] = (lr[:, i1000] - lr[:, i750]) - (lr[:, i750] - lr[:, i560])
        o[pre + 'curv_wl'] = lr[:, i1000] - (lr[:, i750] + (lr[:, i750] - lr[:, i560]) * (1000 - 750) / (750 - 560))
    return o

d1 = derive(out, '1', 'S1_'); d2 = derive(out, '2', 'S2_')
out = pd.concat([out, d1, d2], axis=1)
out.to_csv('tiles_qc/merged_raw.csv', index=False)
print('features', nF, 'complete', len(comp), 'partial', len(part), 'dup complete', dup_complete, 'complete&partial', len(both))
print('craters with n1_est>=20:', int(((out.kind == 0) & (out.n1_est >= 20)).sum()), '/', int((out.kind == 0).sum()))
print('hollows with n1>=20:', int(((out.kind == 1) & (out.n1 >= 20)).sum()), '/', int((out.kind == 1).sum()))
print(out[out.kind == 0][['n1', 'n1_est', 'n2_est', 'S1_albedo_gm', 'S1_slope', 'S1_curv', 'cnt1']].describe().T)
