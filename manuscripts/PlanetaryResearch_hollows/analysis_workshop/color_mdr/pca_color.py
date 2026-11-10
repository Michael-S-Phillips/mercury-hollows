"""Assemble crater/hollow color tables, run PCA on brightness-normalized log spectra, compare hosts vs non-hosts.

Inputs : tiles_qc/merged_raw.csv (merge_spectra.py), crater_features.csv, test1_rows_inc.csv
Outputs: crater_color.csv, hollow_color.csv, pca_record.json, host_vs_nonhost.csv,
         fig_pca_loadings.png, fig_pc1_map.png, fig_pc_scatter.png
PCA: sklearn PCA on the 8 log-ratio values nl_k = ln r_k - mean(ln r) of crater interiors (<=0.9 R,
hollow-masked) for craters with n_pix >= 20; centered, not scaled. Rank <= 7 (rows sum to zero).
Sign conventions: PC1 + => loading(1000) > loading(430) (redder); PC2 + => loading(430)+loading(1000)
- 2 loading(630) > 0 (concave-up, relative low near 600-650 nm); PC3 + => loading(630) > loading(560).
"""
import json, os, numpy as np, pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

W = [430, 480, 560, 630, 750, 830, 900, 1000]
BASE = os.path.expanduser('~/mnt/hollows-emissivity')
CF = os.path.join(BASE, 'manuscripts/PlanetaryResearch_hollows/analysis_workshop/summary_model/crater_features.csv')
HF = os.path.join(BASE, 'manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv')
NMIN = 20
HILAT = 70.0

m = pd.read_csv('tiles_qc/merged_raw.csv')
cf = pd.read_csv(CF)
ho = pd.read_csv(HF)


def block(df, S, pre):
    o = pd.DataFrame(index=df.index)
    st2 = df.stride.values ** 2
    tag = S[1]
    o[pre + 'n_pix'] = df['n%s_est' % tag].values
    o[pre + 'n_sampled'] = df['n' + tag].values
    o[pre + 'img_count_med'] = df['cnt' + tag].values
    for w in W:
        o['%sr%d' % (pre, w)] = df['%s_r%d' % (S, w)].values
    o[pre + 'albedo_gm'] = df[S + '_albedo_gm'].values
    for w in W:
        o['%snl%d' % (pre, w)] = df['%s_nl%d' % (S, w)].values
    for k in ['slope', 'curv', 'curv_wl']:
        o[pre + k] = df['%s_%s' % (S, k)].values
    return o

# ---------------- crater table
mc = m[m.kind == 0].sort_values('fid').reset_index(drop=True)
assert (mc.fid.values == np.arange(len(cf))).all()
cr = pd.DataFrame(dict(id=cf.id.values, lat=cf.lat.values, lon=cf.lon.values, D=cf.D.values, stride=mc.stride.values))
inner = block(mc, 'S1', '')
ann = block(mc, 'S2', 'ann_')
cr = pd.concat([cr, inner, ann], axis=1)
cr['n_masked_hollow'] = (mc.n_reg1_masked.fillna(0).values * mc.stride.values ** 2)
cr['frac_masked_hollow'] = mc.n_reg1_masked / mc.n_reg1_all.replace(0, np.nan)
cr['ok'] = cr.n_pix >= NMIN
cr['hilat'] = cr.lat.abs() > HILAT

NL = ['nl%d' % w for w in W]
ok = cr.ok.values
X = cr.loc[ok, NL].values
pca = PCA().fit(X)
comp = pca.components_.copy()
i = {w: k for k, w in enumerate(W)}
sg = [np.sign(comp[0, i[1000]] - comp[0, i[430]]),
      np.sign(comp[1, i[430]] + comp[1, i[1000]] - 2 * comp[1, i[630]]),
      np.sign(comp[2, i[630]] - comp[2, i[560]])]
for k in range(3):
    comp[k] *= sg[k]
mu = pca.mean_


def project(df, pre=''):
    Y = df[[pre + c for c in NL]].values - mu
    return Y @ comp[:3].T

P = project(cr)
for k in range(3):
    cr['pc%d' % (k + 1)] = P[:, k]
    cr.loc[~ok, 'pc%d' % (k + 1)] = np.nan
Pa = project(cr, 'ann_')
for k in range(3):
    cr['ann_pc%d' % (k + 1)] = Pa[:, k]
    cr.loc[~(cr.ann_n_pix >= NMIN), 'ann_pc%d' % (k + 1)] = np.nan
cr.to_csv('crater_color.csv', index=False)

# ---------------- hollow table
mh = m[m.kind == 1].sort_values('fid').reset_index(drop=True)
hv = pd.DataFrame(dict(row=mh.fid.values, lat=ho.latitude.values, lon=ho.longitude.values, depth=ho.depth.values,
                       stride=mh.stride.values))
hs = block(mh, 'S1', 'set_')     # 10 km, hollow-masked (setting)
ha = block(mh, 'S2', 'all_')     # 10 km, unmasked
hv = pd.concat([hv, hs, ha], axis=1)
hv['set_n_masked_hollow'] = mh.n_reg1_masked.values
for pre in ['set_', 'all_']:
    Ph = project(hv, pre)
    for k in range(3):
        hv['%spc%d' % (pre, k + 1)] = np.where(hv[pre + 'n_pix'] >= NMIN, Ph[:, k], np.nan)
hv.to_csv('hollow_color.csv', index=False)

# ---------------- PCA record
c = cr[ok].join(cf[['h09', 'lrm', 'nsp', 'abslat', 'logD', 'age', 'rayed']].loc[ok])
def corr(a, b):
    v = np.isfinite(a) & np.isfinite(b)
    return dict(pearson=round(float(np.corrcoef(a[v], b[v])[0, 1]), 3), spearman=round(float(stats.spearmanr(a[v], b[v])[0]), 3), n=int(v.sum()))
cors = {('PC%d' % k): {v: corr(c['pc%d' % k].values, c[v].values.astype(float))
                        for v in ['slope', 'curv', 'curv_wl', 'albedo_gm', 'lrm', 'nsp', 'abslat', 'logD', 'img_count_med']}
        for k in (1, 2, 3)}
# sensitivity: PCA restricted to |lat| <= 65 (equirectangular tiles, better photometry), and on raw log reflectance
sub = (c.lat.abs() <= 65).values
p65 = PCA().fit(X[sub])
cos65 = [float(abs(np.dot(p65.components_[k], comp[k]))) for k in range(3)]
praw = PCA().fit(np.log(c[['r%d' % w for w in W]].values))
raw_pc1_vs_albedo = float(np.corrcoef(praw.transform(np.log(c[['r%d' % w for w in W]].values))[:, 0], np.log(c.albedo_gm))[0, 1])
# residual albedo dependence after removing latitude band means
def lrm_nsp():
    out = {}
    lb = pd.cut(c.lrm, [-0.001, 0.0, 0.25, 0.5, 1.0], labels=['0', '(0,0.25]', '(0.25,0.5]', '(0.5,1]'])
    out['by_lrm_fraction'] = {str(k): dict(n=int(len(g)), pc1=round(float(g.pc1.median()), 4), pc2=round(float(g.pc2.median()), 4),
                              pc3=round(float(g.pc3.median()), 4), slope=round(float(g.slope.median()), 4), albedo_gm=round(float(g.albedo_gm.median()), 4))
                              for k, g in c.groupby(lb, observed=True)}
    out['by_nsp'] = {str(k): dict(n=int(len(g)), pc1=round(float(g.pc1.median()), 4), pc2=round(float(g.pc2.median()), 4),
                     pc3=round(float(g.pc3.median()), 4), slope=round(float(g.slope.median()), 4), albedo_gm=round(float(g.albedo_gm.median()), 4))
                     for k, g in c.groupby('nsp')}
    # same comparisons restricted to |lat| <= 65 to avoid polar photometry
    cc = c[c.lat.abs() <= 65]
    out['lrm_ge_0.5_vs_0_within_65deg'] = dict(
        pc1=[round(float(cc[cc.lrm >= 0.5].pc1.median()), 4), round(float(cc[cc.lrm == 0].pc1.median()), 4)],
        pc2=[round(float(cc[cc.lrm >= 0.5].pc2.median()), 4), round(float(cc[cc.lrm == 0].pc2.median()), 4)],
        albedo=[round(float(cc[cc.lrm >= 0.5].albedo_gm.median()), 4), round(float(cc[cc.lrm == 0].albedo_gm.median()), 4)],
        auc_pc1=round(float(stats.mannwhitneyu(cc[cc.lrm >= 0.5].pc1, cc[cc.lrm == 0].pc1)[0] / ((cc.lrm >= 0.5).sum() * (cc.lrm == 0).sum())), 3),
        auc_pc2=round(float(stats.mannwhitneyu(cc[cc.lrm >= 0.5].pc2, cc[cc.lrm == 0].pc2)[0] / ((cc.lrm >= 0.5).sum() * (cc.lrm == 0).sum())), 3),
        auc_albedo=round(float(stats.mannwhitneyu(cc[cc.lrm >= 0.5].albedo_gm, cc[cc.lrm == 0].albedo_gm)[0] / ((cc.lrm >= 0.5).sum() * (cc.lrm == 0).sum())), 3))
    nn = cc[cc.lat > 0]
    out['nsp1_vs_nsp0_north_within_65deg'] = dict(
        pc1=[round(float(nn[nn.nsp == 1].pc1.median()), 4), round(float(nn[nn.nsp == 0].pc1.median()), 4)],
        pc2=[round(float(nn[nn.nsp == 1].pc2.median()), 4), round(float(nn[nn.nsp == 0].pc2.median()), 4)],
        albedo=[round(float(nn[nn.nsp == 1].albedo_gm.median()), 4), round(float(nn[nn.nsp == 0].albedo_gm.median()), 4)],
        auc_pc1=round(float(stats.mannwhitneyu(nn[nn.nsp == 1].pc1, nn[nn.nsp == 0].pc1)[0] / ((nn.nsp == 1).sum() * (nn.nsp == 0).sum())), 3))
    return out

# coverage
cf2 = cf.assign(ok=cr.ok.values, band=pd.cut(cf.abslat, [0, 30, 50, 65, 70, 90], include_lowest=True))
cov = {str(k): dict(n=int(len(g)), frac_ge20=round(float(g.ok.mean()), 4), hosts=int(g.h09.sum()),
                    hosts_frac_ge20=round(float(g[g.h09 == 1].ok.mean()), 4) if g.h09.sum() else None)
       for k, g in cf2.groupby('band', observed=True)}
Dcov = {str(k): round(float(g.ok.mean()), 4) for k, g in cf2.groupby(pd.cut(cf2.D, [0, 5, 7, 10, 20, 1e4]), observed=True)}

# host vs non-host by latitude band
rows = []
c['band'] = pd.cut(c.abslat, [0, 30, 50, 65, 90], include_lowest=True)
for b, g in list(c.groupby('band', observed=True)) + [('all', c), ('all |lat|<=65', c[c.abslat <= 65])]:
    h = g[g.h09 == 1]; nh = g[g.h09 == 0]
    for v in ['pc1', 'pc2', 'pc3', 'slope', 'curv', 'albedo_gm', 'ann_pc1', 'ann_pc2', 'ann_slope', 'ann_albedo_gm', 'logD', 'lrm']:
        a, bb = h[v].dropna(), nh[v].dropna()
        if len(a) < 3:
            continue
        u = stats.mannwhitneyu(a, bb)
        rows.append(dict(band=str(b), var=v, n_host=len(a), n_nonhost=len(bb), med_host=a.median(), med_nonhost=bb.median(),
                         diff=a.median() - bb.median(), diff_in_nonhost_sd=(a.median() - bb.median()) / bb.std(),
                         auc=u[0] / (len(a) * len(bb)), p_mw=u[1]))
hv_tab = pd.DataFrame(rows)
hv_tab.to_csv('host_vs_nonhost.csv', index=False)
# size-matched check (hosts are larger): compare hosts with non-hosts in same logD decile and band
c['Dq'] = pd.qcut(c.logD, 10, labels=False)
sm = {}
for v in ['pc1', 'pc2', 'pc3', 'slope', 'curv', 'albedo_gm']:
    diffs, wts = [], []
    for (b, q), g in c[c.abslat <= 65].groupby(['band', 'Dq'], observed=True):
        h = g[g.h09 == 1][v].dropna(); nh = g[g.h09 == 0][v].dropna()
        if len(h) >= 3 and len(nh) >= 10:
            diffs.append(h.median() - nh.median()); wts.append(len(h))
    sm[v] = round(float(np.average(diffs, weights=wts)), 4) if diffs else None

# hollow settings
hok = hv[hv.set_n_pix >= NMIN]
hollow_summary = dict(n=int(len(hv)), n_set_ge20=int((hv.set_n_pix >= NMIN).sum()),
                      median_set=dict(pc1=float(hok.set_pc1.median()), pc2=float(hok.set_pc2.median()), pc3=float(hok.set_pc3.median()),
                                      slope=float(hok.set_slope.median()), albedo_gm=float(hok.set_albedo_gm.median())),
                      median_all=dict(pc1=float(hok.all_pc1.median()), pc2=float(hok.all_pc2.median()), pc3=float(hok.all_pc3.median()),
                                      slope=float(hok.all_slope.median()), albedo_gm=float(hok.all_albedo_gm.median())),
                      median_frac_masked=float((hv.set_n_masked_hollow / hv.all_n_sampled).median()),
                      corr_depth=dict({v: corr(hok[v].values, hok.depth.values) for v in ['set_pc1', 'set_pc2', 'set_pc3', 'set_slope', 'set_curv', 'set_albedo_gm']}))
crater_ref = dict(pc1_median_all_craters=float(c.pc1.median()), pc1_median_hosts=float(c[c.h09 == 1].pc1.median()),
                  pc2_median_all_craters=float(c.pc2.median()), pc2_median_hosts=float(c[c.h09 == 1].pc2.median()))

rec = dict(
    method=dict(region='crater interior d <= 0.9 R, annulus 1.2-2.0 R; hollow-masked (r_excl = 1 km + min(D_eq,5 km)/2 around 6,731 canonical hollows; 1 km around 2,518 depth-measured hollows)',
                pixel_qc='all 8 bands > 0; r750 >= 0.035; 0.40 <= r430/r750 <= 0.85; 0.95 <= r1000/r750 <= 1.60',
                normalization='nl_k = ln r_k - mean_k ln r_k (divide by geometric mean); albedo_gm = exp(mean ln r)',
                indices=dict(slope='ln(r750/r430)', curv='ln(r1000/r750) - ln(r750/r560)', curv_wl='ln r1000 - [ln r750 + (ln r750 - ln r560)*250/190]'),
                pca='sklearn PCA, centered, unscaled, on nl430..nl1000 of craters with n_pix >= 20; signs fixed as in script docstring',
                n_min=NMIN, statistic='per-band median of valid pixels; spectra then normalized'),
    n_craters=int(len(cr)), n_pca=int(ok.sum()),
    variance_ratio=[round(float(v), 5) for v in pca.explained_variance_ratio_],
    eigenvalues=[float(v) for v in pca.explained_variance_],
    mean_nl=dict(zip(map(str, W), [round(float(v), 4) for v in mu])),
    loadings={('PC%d' % (k + 1)): dict(zip(map(str, W), [round(float(v), 4) for v in comp[k]])) for k in range(3)},
    pc_sd=[float(np.sqrt(pca.explained_variance_[k])) for k in range(3)],
    correlations=cors,
    sensitivity=dict(abs_cosine_loadings_all_vs_lat65=[round(v, 4) for v in cos65],
                     variance_ratio_lat65=[round(float(v), 5) for v in p65.explained_variance_ratio_[:4]],
                     raw_logR_pca_variance=[round(float(v), 5) for v in praw.explained_variance_ratio_[:3]],
                     raw_logR_pc1_corr_with_ln_albedo=round(abs(raw_pc1_vs_albedo), 4)),
    lrm_nsp=lrm_nsp(),
    coverage=dict(frac_ge20_all=round(float(cr.ok.mean()), 4), n_ge20=int(cr.ok.sum()), by_abslat=cov, by_D_km=Dcov,
                  frac_hosts_ge20=round(float(cr.ok[cf.h09.values == 1].mean()), 4)),
    host_size_matched_diff_lat65=sm,
    hosts_frac_masked_median=float(cr[(cf.h09 == 1).values & ok].frac_masked_hollow.median()),
    hollow_settings=hollow_summary, crater_ref=crater_ref)
json.dump(rec, open('pca_record.json', 'w'), indent=1, default=float)

# ---------------- figures
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False})
div = LinearSegmentedColormap.from_list('bluegrayred', ['#104281', '#3987e5', '#b7d3f6', '#f0efec', '#f4b9b2', '#e34948', '#8c1f20'])

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for k, col in enumerate([BLUE, ORANGE, AQUA]):
    ax[0].plot(W, comp[k], '-o', color=col, lw=2, ms=5, label='PC%d (%.1f%%)' % (k + 1, 100 * pca.explained_variance_ratio_[k]))
ax[0].axhline(0, color=GRID, lw=1, zorder=0)
ax[0].set_xlabel('wavelength (nm)'); ax[0].set_ylabel('loading on ln r$_k$ − mean ln r (unitless)')
ax[0].legend(frameon=False); ax[0].set_title('PCA loadings, crater interiors (n = %d)' % ok.sum(), fontsize=9, loc='left')
q = np.percentile(X, [5, 50, 95], axis=0)
ax[1].fill_between(W, q[0], q[2], color='#cde2fb', lw=0, label='5–95th percentile')
ax[1].plot(W, q[1], '-o', color=BLUE, lw=2, ms=5, label='median')
ax[1].set_xlabel('wavelength (nm)'); ax[1].set_ylabel('ln r$_k$ − mean ln r (unitless)')
ax[1].legend(frameon=False, loc='upper left'); ax[1].set_title('Brightness-normalized crater spectra', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig('fig_pca_loadings.png', dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 5.2))
cc = cr[ok].sort_values('D', ascending=False)
lim = np.nanpercentile(np.abs(cc.pc1), 98)
sz = np.clip(2 + 6 * np.log10(cc.D / 5), 1.5, 30)
sc = ax.scatter(cc.lon, cc.lat, c=cc.pc1, s=sz, cmap=div, norm=TwoSlopeNorm(0, -lim, lim), lw=0, rasterized=True)
ax.set_xlim(0, 360); ax.set_ylim(-90, 90); ax.set_aspect('equal')
ax.set_xticks(range(0, 361, 60)); ax.set_yticks(range(-90, 91, 30))
ax.set_xlabel('longitude (°E)'); ax.set_ylabel('latitude (°N)')
for y in (-HILAT, HILAT):
    ax.axhline(y, color=MUTED, lw=0.6, ls=':')
cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.01, extend='both')
cb.set_label('crater PC1 (spectral slope; + = redder), unitless')
ax.set_title('Crater-interior PC1 (n = %d; dot size ∝ log D; dotted lines = |lat| %d°, photometry less reliable poleward)' % (ok.sum(), HILAT), fontsize=9, loc='left')
fig.tight_layout(); fig.savefig('fig_pc1_map.png', dpi=150); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
cc = c
for a_, (x, y) in zip(ax, [('pc1', 'pc2'), ('pc1', 'pc3')]):
    nh = cc[cc.h09 == 0]; h = cc[cc.h09 == 1]
    a_.scatter(nh[x], nh[y], s=2, color='#b9b8b3', lw=0, alpha=0.6, label='non-host craters (n = %d)' % len(nh), rasterized=True)
    a_.scatter(h[x], h[y], s=14, color=ORANGE, edgecolor='white', lw=0.5, label='hollow hosts, h09 = 1 (n = %d)' % len(h))
    a_.set_xlabel('%s (unitless)' % x.upper()); a_.set_ylabel('%s (unitless)' % y.upper())
    lx = np.nanpercentile(cc[x], [0.2, 99.8]); ly = np.nanpercentile(cc[y], [0.2, 99.8])
    a_.set_xlim(lx[0] - 0.1 * np.ptp(lx), lx[1] + 0.1 * np.ptp(lx)); a_.set_ylim(ly[0] - 0.1 * np.ptp(ly), ly[1] + 0.1 * np.ptp(ly))
    a_.axhline(0, color=GRID, lw=1, zorder=0); a_.axvline(0, color=GRID, lw=1, zorder=0)
ax[0].legend(frameon=False, loc='upper left', markerscale=1.5)
ax[0].set_title('PC1 (slope) vs PC2 (curvature)', fontsize=9, loc='left'); ax[1].set_title('PC1 vs PC3 (560–630 nm step)', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig('fig_pc_scatter.png', dpi=150); plt.close(fig)
print(json.dumps({k: rec[k] for k in ['n_pca', 'variance_ratio', 'loadings', 'correlations', 'sensitivity', 'lrm_nsp', 'coverage', 'host_size_matched_diff_lat65', 'hosts_frac_masked_median', 'hollow_settings', 'crater_ref']}, indent=0, default=float))
