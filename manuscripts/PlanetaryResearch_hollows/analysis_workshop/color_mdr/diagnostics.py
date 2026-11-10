"""Extra diagnostics: residual albedo/latitude dependence, tile seams in PC1, hollow mask effect, size-matched host contrasts."""
import json, numpy as np, pandas as pd
from scipy import stats
cr = pd.read_csv('crater_color.csv'); cf = pd.read_csv('../summary_model/crater_features.csv')
hv = pd.read_csv('hollow_color.csv'); ho = pd.read_csv('../../../GRL_insolation_hollows/data/test1_rows_inc.csv')
c = cr.join(cf[['h09', 'lrm', 'nsp', 'abslat', 'logD']])
c = c[c.ok]
out = {}
def r(a, b):
    v = np.isfinite(a) & np.isfinite(b); return round(float(stats.spearmanr(a[v], b[v])[0]), 3)
out['spearman_pc1_albedo_by_band'] = {b: r(g.pc1.values, g.albedo_gm.values) for b, g in c.groupby(pd.cut(c.abslat, [0, 30, 50, 65, 90]).astype(str))}
out['spearman_pc2_albedo_by_band'] = {b: r(g.pc2.values, g.albedo_gm.values) for b, g in c.groupby(pd.cut(c.abslat, [0, 30, 50, 65, 90]).astype(str))}
s = c[c.lat < 0]
out['south_only_spearman_pc1_abslat'] = r(s.pc1.values, s.abslat.values)
out['north_nonNSP_spearman_pc1_abslat'] = r(c[(c.lat > 0) & (c.nsp == 0)].pc1.values, c[(c.lat > 0) & (c.nsp == 0)].abslat.values)
out['pc1_median_by_abslat_south'] = {b: round(float(g.pc1.median()), 4) for b, g in s.groupby(pd.cut(s.abslat, [0, 20, 40, 55, 65, 75, 90]).astype(str))}
# seams: craters whose 2R circle lies within 3 deg either side of a boundary
seam = {}
for name, lat0 in [('+65 polar/equirect', 65), ('-65 polar/equirect', -65), ('+43.75', 43.75), ('-43.75', -43.75), ('+22.5', 22.5), ('-22.5', -22.5), ('0', 0)]:
    a = c[(c.lat > lat0) & (c.lat < lat0 + 2) & (c.D < 30)]; b = c[(c.lat < lat0) & (c.lat > lat0 - 2) & (c.D < 30)]
    seam[name] = dict(n_above=len(a), n_below=len(b), pc1_above=round(float(a.pc1.median()), 4), pc1_below=round(float(b.pc1.median()), 4),
                      pc2_above=round(float(a.pc2.median()), 4), pc2_below=round(float(b.pc2.median()), 4),
                      albedo_above=round(float(a.albedo_gm.median()), 4), albedo_below=round(float(b.albedo_gm.median()), 4),
                      p_pc1=float(stats.mannwhitneyu(a.pc1, b.pc1)[1]))
out['seams_2deg_bands_D<30km'] = seam
out['pc_sd_all'] = [round(float(c['pc%d' % k].std()), 4) for k in (1, 2, 3)]
# size-matched host contrast per band (band x logD decile), in units of non-host SD
c['Dq'] = pd.qcut(c.logD, 10, labels=False)
sm = {}
for b, gb in c.groupby(pd.cut(c.abslat, [0, 30, 50, 65, 90]).astype(str)):
    row = {}
    for v in ['pc1', 'pc2', 'pc3', 'slope', 'curv', 'albedo_gm']:
        d, w = [], []
        for q, g in gb.groupby('Dq'):
            h = g[g.h09 == 1][v]; nh = g[g.h09 == 0][v]
            if len(h) >= 2 and len(nh) >= 10:
                d.append(h.median() - nh.median()); w.append(len(h))
        row[v] = round(float(np.average(d, weights=w) / gb[gb.h09 == 0][v].std()), 3) if d else None
    row['n_host'] = int((gb.h09 == 1).sum())
    sm[b] = row
out['host_minus_nonhost_size_matched_in_SD_by_band'] = sm
# hollow mask effect: paired (unmasked - masked) within 10 km of the measured hollows
h = hv[(hv.set_n_pix >= 20)]
pair = {}
for v in ['pc1', 'pc2', 'pc3', 'slope', 'curv', 'albedo_gm']:
    d = h['all_' + v] - h['set_' + v]
    pair[v] = dict(median=float(d.median()), frac_positive=round(float((d > 0).mean()), 3), p_wilcoxon=float(stats.wilcoxon(d)[1]))
out['hollow_unmasked_minus_masked_10km'] = pair
# hollow setting vs depth: raw and partial (on abslat) Spearman
def partial_spearman(x, y, z):
    rx = stats.rankdata(x); ry = stats.rankdata(y); rz = stats.rankdata(z)
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return round(float(np.corrcoef(ex, ey)[0, 1]), 3)
hh = h.join(ho[['abslat', 'inc', 'lnpx']])
out['hollow_depth_vs_setting'] = {v: dict(spearman=r(hh['set_' + v].values, hh.depth.values),
                                          partial_on_abslat=partial_spearman(hh['set_' + v].values, hh.depth.values, hh.abslat.values),
                                          partial_on_inc=partial_spearman(hh['set_' + v].values, hh.depth.values, hh.inc.values))
                                  for v in ['pc1', 'pc2', 'pc3', 'slope', 'curv', 'albedo_gm']}
out['hollow_setting_pc1_vs_abslat'] = r(hh.set_pc1.values, hh.abslat.values)
out['hollow_lat_quantiles'] = [round(float(x), 1) for x in hh.lat.quantile([.1, .5, .9])]
json.dump(out, open('diagnostics_record.json', 'w'), indent=1)
print(json.dumps(out, indent=0))

# ---- local contrast (interior minus 1.2-2.0 R annulus): cancels regional mosaic/photometric offsets
c['d_pc1'] = c.pc1 - c.ann_pc1; c['d_pc2'] = c.pc2 - c.ann_pc2; c['d_pc3'] = c.pc3 - c.ann_pc3
c['d_slope'] = c.slope - c.ann_slope; c['d_lnalb'] = np.log(c.albedo_gm / c.ann_albedo_gm)
loc = {}
for b, gb in list(c.groupby(pd.cut(c.abslat, [0, 30, 50, 65, 90]).astype(str))) + [('all', c)]:
    row = {}
    for v in ['d_pc1', 'd_pc2', 'd_pc3', 'd_slope', 'd_lnalb']:
        d, w = [], []
        for q, g in gb.groupby('Dq'):
            h = g[g.h09 == 1][v].dropna(); nh = g[g.h09 == 0][v].dropna()
            if len(h) >= 2 and len(nh) >= 10:
                d.append(h.median() - nh.median()); w.append(len(h))
        h = gb[gb.h09 == 1][v].dropna(); nh = gb[gb.h09 == 0][v].dropna()
        row[v] = dict(med_host=round(float(h.median()), 4), med_nonhost=round(float(nh.median()), 4),
                      size_matched_diff_SD=round(float(np.average(d, weights=w) / nh.std()), 3) if d else None,
                      p_mw_unmatched=float(stats.mannwhitneyu(h, nh)[1]))
    row['n_host'] = int((gb.h09 == 1).sum())
    loc[b] = row
out['local_contrast_interior_minus_annulus'] = loc
out['spearman_d_pc1_abslat'] = r(c.d_pc1.values, c.abslat.values)
out['spearman_d_pc1_lnalb'] = r(c.d_pc1.values, c.d_lnalb.values)
json.dump(out, open('diagnostics_record.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ['local_contrast_interior_minus_annulus', 'spearman_d_pc1_abslat', 'spearman_d_pc1_lnalb']}, indent=0))
