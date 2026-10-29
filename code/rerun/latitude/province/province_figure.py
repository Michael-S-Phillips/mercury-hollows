#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
"""Province figure: (a) map with LRM, fields, host craters; (b) province/rest
composition & crater-density ratios; (c) obs/null under null reweighting."""
import numpy as np, pandas as pd, json
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
# figure-style helpers: import from the skill's kernel.py when available,
# else expect them pre-injected into globals (Claude Science kernel context).
if 'apply_figure_style' not in dir():
    import importlib.util as _ilu, glob as _glob
    _c = _glob.glob('/Users/phillipsm/.claude-science/runtime/*/skills/figure-style/kernel.py')
    _spec = _ilu.spec_from_file_location('figure_style', sorted(_c)[-1])
    _fs = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_fs)
    apply_figure_style, panel_letter, META_GREY = _fs.apply_figure_style, _fs.panel_letter, _fs.META_GREY
apply_figure_style()

HERE = _REPO + '/code/rerun/latitude/province'
rec = json.load(open(HERE + '/province_record.json'))
fd = np.load(HERE + '/_figdata.npz')
PLAT, PLON = (40, 70), (90, 150)
FOCAL = '#c0392b'   # province
REST = '#5a6b7a'    # comparison

# ---- data ----
comp = rec['comparison']
ck = {c['key']: c for c in comp}
KL10 = fd['KL10']; latc = fd['latc']; lonc = fd['lonc']
flat, flon, fn = fd['field_lat'], fd['field_lon'], fd['field_n']
crlat, crlon, crD, crnh = fd['cr_lat'], fd['cr_lon'], fd['cr_D'], fd['cr_nh']

fig = plt.figure(figsize=(7.2, 8.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 1.0, 0.9],
                      hspace=0.42, wspace=0.34,
                      left=0.10, right=0.97, top=0.95, bottom=0.07)

# ============================================================ (a) map
axm = fig.add_subplot(gs[0, :])
mlon = (85, 172); mlat = (25, 75)
ci0 = int(mlon[0]*10); ci1 = int(mlon[1]*10)
ri0 = int((90-mlat[1])*10); ri1 = int((90-mlat[0])*10)
sub = KL10[ri0:ri1, ci0:ci1]
im = axm.imshow(sub, extent=[mlon[0], mlon[1], mlat[0], mlat[1]], origin='upper',
                cmap='cividis', vmin=0, vmax=60, aspect='auto', zorder=0)
# province box
axm.add_patch(Rectangle((PLON[0], PLAT[0]), PLON[1]-PLON[0], PLAT[1]-PLAT[0],
                        fill=False, ec=FOCAL, lw=1.8, zorder=5))
# fields inside map extent
inmap = (flon >= mlon[0]) & (flon <= mlon[1]) & (flat >= mlat[0]) & (flat <= mlat[1])
inprov = (flon >= PLON[0]) & (flon < PLON[1]) & (flat >= PLAT[0]) & (flat < PLAT[1])
sz = 4 + fn*0.9
axm.scatter(flon[inmap & ~inprov], flat[inmap & ~inprov], s=sz[inmap & ~inprov],
            facecolor='none', edgecolor='white', lw=0.5, alpha=0.7, zorder=6)
axm.scatter(flon[inprov], flat[inprov], s=sz[inprov], facecolor=FOCAL,
            edgecolor='white', lw=0.4, alpha=0.9, zorder=7)
# host craters (exclude Caloris giant + placeholder -999)
hc = (crlat > -900) & (crD < 300)
axm.scatter(crlon[hc], crlat[hc], s=(crD[hc]*1.4), facecolor='none',
            edgecolor='#f5f5f5', lw=0.9, ls='--', zorder=8)
# Caloris direction marker
axm.scatter([161.9], [31.9], marker='+', s=90, color='#ffd54f', lw=1.8, zorder=9)
axm.annotate('Caloris center', (161.9, 31.9), (150, 27.5), color='#ffd54f',
             fontsize=6, ha='center')
# label two largest fields + named craters
def lab(lo, la, txt, dx, dy, col='white'):
    axm.annotate(txt, (lo, la), (lo+dx, la+dy), color=col, fontsize=6,
                 arrowprops=dict(arrowstyle='-', color=col, lw=0.5), zorder=10)
lab(113.21, 55.02, 'n=124', -7, 4)
lab(132.53, 48.95, 'n=116', 4, 5)
for nm, la, lo in [('Kofi', 56.69, 118.2), ('Oskison', 60.41, 144.9), ('Pasch', 46.09, 134.76)]:
    axm.annotate(nm, (lo, la), (lo, la-3.2), color='#aee', fontsize=6, ha='center',
                 style='italic', zorder=10)
axm.set_xlim(*mlon); axm.set_ylim(*mlat)
axm.set_xlabel('Longitude (\u00b0E)'); axm.set_ylabel('Latitude (\u00b0N)')
axm.set_title('90\u2013150\u00b0E, 40\u201370\u00b0N province: hollow fields sit on LRM-poor terrane NE of Caloris',
              fontsize=8)
cb = fig.colorbar(im, ax=axm, fraction=0.025, pad=0.01)
cb.set_label('Klima LRM brightness', fontsize=6)
# legend
from matplotlib.lines import Line2D
leg = [Line2D([0],[0], marker='o', ls='', mfc=FOCAL, mec='white', ms=6, label='field in province'),
       Line2D([0],[0], marker='o', ls='', mfc='none', mec='white', ms=6, label='field elsewhere'),
       Line2D([0],[0], marker='o', ls='--', mfc='none', mec='#f5f5f5', ms=8, label='host crater')]
axm.legend(handles=leg, loc='upper left', fontsize=5.5, framealpha=0.25,
           facecolor='black', labelcolor='white')
panel_letter(axm, 'a')

# ============================================================ (b) ratios
axb = fig.add_subplot(gs[1, :])
order = ['klima', 'klima_flag', 'peplowski_flag', 'xrs_Mg_Si', 'xrs_S_Si',
         'xrs_Ca_Si', 'xrs_Al_Si', 'crater_dens_10_20', 'crater_dens_20_50', 'crater_dens_50_inf']
labels = ['Klima LRM bright.', 'Klima LRM flag frac', 'Peplowski LRM flag frac',
          'XRS Mg/Si', 'XRS S/Si', 'XRS Ca/Si', 'XRS Al/Si',
          'Crater dens 10\u201320 km', 'Crater dens 20\u201350 km', 'Crater dens >50 km']
y = np.arange(len(order))[::-1]
for yi, k in zip(y, order):
    c = ck[k]
    P, Ps = c['prov_val'], c['prov_sd']; Rr, Rs = c['rest_val'], c['rest_sd']
    ratio = P / Rr
    # delta-method sd of ratio
    sd = ratio * np.sqrt((Ps/P)**2 + (Rs/Rr)**2) if P > 0 else 0
    col = FOCAL if c['significant'] else META_GREY
    axb.plot([1, ratio], [yi, yi], color=col, lw=1.0, zorder=2)
    axb.errorbar(ratio, yi, xerr=sd, fmt='o', color=col, ms=5, capsize=2, zorder=3)
axb.axvline(1.0, color='k', lw=0.8, ls=':')
axb.set_yticks(y); axb.set_yticklabels(labels, fontsize=6.5)
axb.set_xlabel('Province \u00f7 rest of 40\u201370\u00b0N  (1 = equal; red = 95% CI excludes 1)')
axb.set_xlim(0, 1.35)
axb.set_title('Province is LRM-poor, volatile/Mg-depleted, and crater-poor', fontsize=8)
panel_letter(axb, 'b')

# ============================================================ (c) reweighting
axc = fig.add_subplot(gs[2, 0])
ex = rec['excess_after_reweighting']
keys = ['base', 'lrm_weighted', 'crater_weighted', 'both_weighted']
labs = ['coverage\nnull', '\u00d7 LRM', '\u00d7 crater\ninterior', '\u00d7 both']
vals = [ex[k]['ratio'] for k in keys]
n = ex['n_fields_province']
errs = [v / np.sqrt(n) for v in vals]  # Poisson (observed count) uncertainty
bars = axc.bar(range(4), vals, color=[REST, FOCAL, REST, FOCAL], width=0.7,
               yerr=errs, capsize=3, edgecolor='k', lw=0.5)
axc.axhline(1.0, color='k', lw=0.8, ls=':')
for i, v in enumerate(vals):
    axc.text(i, v + errs[i] + 0.5, f'{v:.1f}', ha='center', fontsize=6.5)
axc.set_xticks(range(4)); axc.set_xticklabels(labs, fontsize=6.5)
axc.set_ylabel('Province observed / null')
axc.set_title('Excess survives reweighting', fontsize=8)
axc.set_ylim(0, max(vals)*1.18)
panel_letter(axc, 'c')

# ============================================================ (c2) counts text
axd = fig.add_subplot(gs[2, 1]); axd.axis('off')
txt = (f"Province (90\u2013150\u00b0E, 40\u201370\u00b0N)\n"
       f"  fields: {rec['n_fields_province']} of {rec['n_fields_40_70N']} at 40\u201370\u00b0N\n"
       f"  hollows in those fields: {rec['n_hollows_province_fields']}\n"
       f"  host craters (\u22641.1R): {rec['n_host_craters_found']}/{rec['n_fields_province']}\n"
       f"  named craters: {', '.join(rec['named_craters'])}\n\n"
       f"obs/null: base {ex['base']['ratio']}, \u00d7LRM {ex['lrm_weighted']['ratio']},\n"
       f"  \u00d7crater {ex['crater_weighted']['ratio']}, \u00d7both {ex['both_weighted']['ratio']}\n"
       f"(N per bar = {n}; error = Poisson on obs count)")
axd.text(0.0, 0.95, txt, va='top', ha='left', fontsize=6.4, family='monospace')

fig.savefig(HERE + '/province.png', dpi=300)

# ---- §9.1 geometric check ----
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
ov = [(a.get_text(), b.get_text()) for i,(a,ba) in enumerate(texts)
      for b,bb in texts[i+1:] if ba.overlaps(bb)
      and a.axes is b.axes and a.axes not in (axm,)]  # ignore map overlays
print('text overlaps (non-map):', len(ov))
print('saved province.png')
