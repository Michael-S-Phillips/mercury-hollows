#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
"""
Province characterization: 90-150E, 40-70N vs rest of 40-70N (and whole northern
covered area). Part of the Mercury hollows latitude study.

Outputs (this folder): province_record.json, province_craters.csv
Figure written by province_figure.py.

Conventions: all grids north at row 0. Longitudes 0-360 from col 0 except XRS
(-180..180) which we roll to 0-360. Coverage/null/lat_grids are 10 ppd. Mercury
R = 2439.7 km.
"""
import numpy as np, pandas as pd, json, os, sys
from datetime import datetime, timezone

R_M = 2439.7  # km
BASE = _REPO + ''
HERE = BASE + '/code/rerun/latitude/province'
V2 = BASE + '/hollow-networks/data_v2'

# Province box
PLAT = (40.0, 70.0); PLON = (90.0, 150.0)
RNG = np.random.default_rng(20240611)
NBOOT = 2000

# ---------------------------------------------------------------- grid resolver
def resolve(fn):
    # artifact paths passed in via environment (written by driver)
    return ART[fn]

# ---------------------------------------------------------------- samplers
def make_sampler(grid, ppd, lon_origin='0'):
    """Return f(lat, lon0_360)->values via nearest-cell. north at row 0."""
    nrow, ncol = grid.shape
    def f(lat, lon):
        lat = np.asarray(lat, float); lon = np.asarray(lon, float) % 360.0
        row = np.floor((90.0 - lat) * ppd).astype(int)
        row = np.clip(row, 0, nrow - 1)
        col = np.floor(lon * ppd).astype(int) % ncol
        return grid[row, col]
    return f

def gc_dist_km(lat1, lon1, lat2, lon2):
    """great-circle distance on Mercury, degrees in."""
    p1 = np.radians(lat1); p2 = np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    d = np.arccos(np.clip(np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(dl), -1, 1))
    return R_M * d

# ---------------------------------------------------------------- load data
def main():
    global ART
    ART = json.load(open('/tmp/province_art.json'))
    fields = pd.read_csv(V2 + '/fields.csv')
    fields['lon'] = fields['lon'] % 360.0
    # grids
    lrm = np.load(ART['lrm_published_16ppd.npz'])
    K16, P16 = lrm['KLIMA'], lrm['PEPLOWSKI']
    grids8 = np.load(ART['mercury_grids_8ppd.npz'])
    DEM, BOUG = grids8['DEM'], grids8['BOUG']
    xrs = np.load(ART['xrs_ratios.npz'])
    # roll XRS -180..180 -> 0-360
    def roll360(a): return np.roll(a, -a.shape[1] // 2, axis=1)
    XRS = {k: roll360(xrs['map_' + k]) for k in ['Mg_Si', 'S_Si', 'Ca_Si', 'Al_Si']}
    cov = np.load(V2 + '/coverage_10ppd.npz')
    nimg = cov['nimg']  # 10 ppd
    latg = np.load(BASE + '/code/rerun/latitude/lat_grids.npz')
    LAM, KL10, CR = latg['lam'], latg['KL10'], latg['CR']

    # roughness: std of DEM in 1-deg windows (8 ppd -> 8x8 window)
    from scipy.ndimage import generic_filter, uniform_filter
    # fast std via uniform_filter of x and x^2
    m = uniform_filter(DEM.astype(float), size=8, mode='wrap')
    m2 = uniform_filter(DEM.astype(float)**2, size=8, mode='wrap')
    ROUGH = np.sqrt(np.maximum(m2 - m**2, 0.0))

    # samplers
    s_K = make_sampler(K16, 16); s_P = make_sampler(P16, 16)
    s_DEM = make_sampler(DEM, 8); s_BOUG = make_sampler(BOUG, 8); s_ROUGH = make_sampler(ROUGH, 8)
    s_xrs = {k: make_sampler(XRS[k], 4) for k in XRS}

    # ---- build covered 10ppd cell table (north 0-90 only, covered) ----
    ppd = 10
    nrow, ncol = nimg.shape
    latc = 90.0 - (np.arange(nrow) + 0.5) / ppd
    lonc = (np.arange(ncol) + 0.5) / ppd
    LON, LAT = np.meshgrid(lonc, latc)
    covered = nimg > 0
    north = LAT >= 0
    sel = covered & north
    cell_lat = LAT[sel]; cell_lon = LON[sel]
    # cell physical area (km2), cos-lat weighting
    dlat = 1.0 / ppd; dlon = 1.0 / ppd
    cell_area = (R_M**2) * np.radians(dlon) * (np.sin(np.radians(cell_lat + dlat/2)) -
                                               np.sin(np.radians(cell_lat - dlat/2)))
    # sample all continuous grids at covered northern cells
    samp = {
        'klima': s_K(cell_lat, cell_lon), 'peplowski': s_P(cell_lat, cell_lon),
        'dem_m': s_DEM(cell_lat, cell_lon), 'bouguer_mgal': s_BOUG(cell_lat, cell_lon),
        'roughness_m': s_ROUGH(cell_lat, cell_lon),
    }
    for k in XRS:
        samp['xrs_' + k] = s_xrs[k](cell_lat, cell_lon)
    cells = pd.DataFrame({'lat': cell_lat, 'lon': cell_lon, 'area_km2': cell_area, **samp})
    cells['klima_flag'] = (cells['klima'] > 0).astype(float)
    cells['peplowski_flag'] = (cells['peplowski'] > 0).astype(float)

    # region masks
    def in_prov(lat, lon):
        return (lat >= PLAT[0]) & (lat < PLAT[1]) & (lon >= PLON[0]) & (lon < PLON[1])
    cells['prov'] = in_prov(cells['lat'], cells['lon'])
    cells['band40_70'] = (cells['lat'] >= 40) & (cells['lat'] < 70)
    cells['rest'] = cells['band40_70'] & (~cells['prov'])
    # 5x5 block id
    cells['blk'] = ((np.floor(cells['lat'] / 5).astype(int)) * 1000 +
                    np.floor(cells['lon'] / 5).astype(int))

    # ---- Herrick craters ----
    h10 = pd.read_parquet('/tmp/h10.parquet'); h5 = pd.read_parquet('/tmp/h5.parquet')
    for h in (h10, h5):
        h['lon'] = h['lon'] % 360.0
    h10['name'] = h10['name'].astype(str).str.strip()
    h5['name'] = ''
    hall = pd.concat([h10[['lat','lon','diam','name']], h5[['lat','lon','diam','name']]],
                     ignore_index=True)

    # continuous grids to summarize
    CONT = ['klima', 'klima_flag', 'peplowski', 'peplowski_flag',
            'xrs_Mg_Si', 'xrs_S_Si', 'xrs_Ca_Si', 'xrs_Al_Si',
            'dem_m', 'roughness_m', 'bouguer_mgal']
    CBINS = [('crater_dens_10_20', 10, 20), ('crater_dens_20_50', 20, 50),
             ('crater_dens_50_inf', 50, 1e9)]

    # ---- per-block aggregates for bootstrap ----
    def block_table(cellmask, cropbox=None):
        """Return per-block dict of (sum_w, sum_w*x per cont, area, crater counts)."""
        sub = cells[cellmask]
        blocks = {}
        for blk, g in sub.groupby('blk'):
            w = g['area_km2'].values  # area-weighted (cos-lat already in area)
            rec = {'sum_w': w.sum(), 'area': g['area_km2'].sum()}
            for c in CONT:
                x = g[c].values.astype(float)
                good = np.isfinite(x)
                rec['w_' + c] = w[good].sum()
                rec['wx_' + c] = (w[good] * x[good]).sum()
            blocks[blk] = rec
        # crater counts per block: assign craters to 5x5 blocks within region lat/lon extent
        blkids = list(blocks.keys())
        for name, lo, hi in CBINS:
            for b in blkids:
                blocks[b][name] = 0
        # bounding region for craters = union of block boxes
        if len(blkids):
            b_lat = np.array([b // 1000 for b in blkids]) * 5
            b_lon = np.array([b % 1000 for b in blkids]) * 5
            for name, lo, hi in CBINS:
                hh = hall[(hall.diam >= lo) & (hall.diam < hi)]
                cblk = (np.floor(hh.lat.values / 5).astype(int) * 1000 +
                        np.floor(hh.lon.values / 5).astype(int))
                vc = pd.Series(cblk).value_counts()
                for b in blkids:
                    if b in vc.index:
                        blocks[b][name] += int(vc[b])
        return blocks

    def region_stats(blocks):
        """point estimate + block bootstrap CI for each quantity."""
        keys = list(blocks.keys())
        arr = {k: np.array([blocks[b][k] for b in keys]) for b in [keys[0]] for k in blocks[keys[0]]}
        def compute(idx):
            out = {}
            area = arr['area'][idx].sum()
            for c in CONT:
                wsum = arr['w_' + c][idx].sum()
                out[c] = arr['wx_' + c][idx].sum() / wsum if wsum > 0 else np.nan
            for name, lo, hi in CBINS:
                out[name] = arr[name][idx].sum() / area * 1e6  # per 1e6 km2
            return out
        pt = compute(np.arange(len(keys)))
        boots = {k: [] for k in pt}
        n = len(keys)
        for _ in range(NBOOT):
            idx = RNG.integers(0, n, n)
            c = compute(idx)
            for k in c: boots[k].append(c[k])
        res = {}
        for k in pt:
            b = np.array(boots[k])
            res[k] = {'value': float(pt[k]), 'sd': float(np.nanstd(b)),
                      'ci95': [float(np.nanpercentile(b, 2.5)), float(np.nanpercentile(b, 97.5))]}
        return res, {k: np.array(boots[k]) for k in pt}

    blk_prov = block_table(cells['prov'].values)
    blk_rest = block_table(cells['rest'].values)
    blk_north = block_table((cells['lat'] >= 0).values)
    st_prov, bt_prov = region_stats(blk_prov)
    st_rest, bt_rest = region_stats(blk_rest)
    st_north, bt_north = region_stats(blk_north)

    # difference significance (province - rest): bootstrap of difference
    def diff_sig(bt_a, bt_b, key):
        a = bt_a[key]; b = bt_b[key]
        n = min(len(a), len(b))
        d = a[:n] - b[:n]
        lo, hi = np.nanpercentile(d, [2.5, 97.5])
        return float(np.nanmean(d)), [float(lo), float(hi)], bool((lo > 0) or (hi < 0))

    comparison = []
    label = {'klima': 'Klima LRM brightness (0-255)', 'klima_flag': 'Klima LRM flag fraction',
             'peplowski': 'Peplowski LRM brightness (0-255)', 'peplowski_flag': 'Peplowski LRM flag fraction',
             'xrs_Mg_Si': 'XRS Mg/Si', 'xrs_S_Si': 'XRS S/Si', 'xrs_Ca_Si': 'XRS Ca/Si', 'xrs_Al_Si': 'XRS Al/Si',
             'dem_m': 'DEM elevation (m)', 'roughness_m': 'DEM roughness 1-deg std (m)',
             'bouguer_mgal': 'Bouguer anomaly (mGal)',
             'crater_dens_10_20': 'Crater density D 10-20 km (per 1e6 km2)',
             'crater_dens_20_50': 'Crater density D 20-50 km (per 1e6 km2)',
             'crater_dens_50_inf': 'Crater density D >50 km (per 1e6 km2)'}
    for k in CONT + [c[0] for c in CBINS]:
        dmean, dci, sig = diff_sig(bt_prov, bt_rest, k)
        def fmt(st): return f"{st[k]['value']:.3g} ± {st[k]['sd']:.2g}"
        comparison.append({'quantity': label[k], 'key': k,
                           'province': fmt(st_prov), 'rest_40_70N': fmt(st_rest),
                           'whole_north': fmt(st_north),
                           'diff_prov_minus_rest': f"{dmean:.3g} [{dci[0]:.3g},{dci[1]:.3g}]",
                           'significant': sig,
                           'prov_val': st_prov[k]['value'], 'prov_sd': st_prov[k]['sd'],
                           'prov_ci': st_prov[k]['ci95'],
                           'rest_val': st_rest[k]['value'], 'rest_sd': st_rest[k]['sd'],
                           'rest_ci': st_rest[k]['ci95'],
                           'north_val': st_north[k]['value']})
    pd.DataFrame([{'quantity': c['quantity'], 'province': c['province'],
                   'rest_40_70N': c['rest_40_70N'], 'whole_north': c['whole_north'],
                   'diff_prov_minus_rest': c['diff_prov_minus_rest'],
                   'significant': c['significant']} for c in comparison]).to_csv(
                       HERE + '/province_comparison.csv', index=False)

    # ---- field counts & obs/null with reweighting ----
    fld = fields.copy()
    fp = fld[in_prov(fld['lat'], fld['lon'])]
    n_prov = len(fp)
    n_4070 = int(((fld['lat'] >= 40) & (fld['lat'] < 70)).sum())

    # 10 ppd cell grids for null reweighting.
    # Use the network coverage null (net_null_intensity_v2) which reproduces the
    # record's 90-120E / 120-150E obs/null = 5.2 / 4.6. Reference domain = covered
    # northern hemisphere, N_ref = northern fields. Enrichment ratio is a fraction
    # so it is normalization-free; reweighting changes the null's SHAPE only.
    nn = np.load(ART['net_null_intensity_v2.npz'])
    LAMn = np.nan_to_num(nn['lam'])
    north10 = LAT >= 0
    ref_mask = covered & north10
    N_north = int((fld['lat'] >= 0).sum())
    def cellmask10(latm, lonm):
        return (LAT >= latm[0]) & (LAT < latm[1]) & (LON >= lonm[0]) & (LON < lonm[1])
    prov_cells = cellmask10(PLAT, PLON) & covered
    def obsnull(weight):
        w = LAMn * weight
        exp = N_north * (w[prov_cells].sum() / w[ref_mask].sum())
        return n_prov / exp, exp
    eps = 1.0
    crmask = np.where(CR > 0, 1.0, 0.0)
    ratio_base, exp_base = obsnull(np.ones_like(LAMn))
    ratio_lrm, exp_lrm = obsnull(KL10 + eps)
    ratio_cr, exp_cr = obsnull(crmask + eps)  # crater-interior emphasis
    ratio_both, exp_both = obsnull((KL10 + eps) * (crmask + eps))
    excess = {
        'n_fields_province': int(n_prov), 'n_fields_40_70N': n_4070, 'N_north_ref': N_north,
        'null_source': 'net_null_intensity_v2.npz (size-marginalized coverage null)',
        'reference_domain': 'covered northern hemisphere (0-90N)',
        'base': {'obs': int(n_prov), 'null_expected': round(float(exp_base), 2), 'ratio': round(float(ratio_base), 2)},
        'lrm_weighted': {'obs': int(n_prov), 'null_expected': round(float(exp_lrm), 2), 'ratio': round(float(ratio_lrm), 2)},
        'crater_weighted': {'obs': int(n_prov), 'null_expected': round(float(exp_cr), 2), 'ratio': round(float(ratio_cr), 2)},
        'both_weighted': {'obs': int(n_prov), 'null_expected': round(float(exp_both), 2), 'ratio': round(float(ratio_both), 2)},
        'eps': eps,
        'weights': 'LRM = Klima brightness KL10 + eps; crater = crater-interior mask (D>=10km) + eps; both = product',
        'note': 'enrichment ratio = n_prov / [N_north * (sum(w*lam)_prov / sum(w*lam)_north)]; reweighting redistributes the null shape'}

    # ---- host craters for province fields ----
    rows = []
    hlat = hall.lat.values; hlon = hall.lon.values; hd = hall.diam.values; hn = hall.name.values
    for _, f in fp.iterrows():
        d = gc_dist_km(f['lat'], f['lon'], hlat, hlon)
        rad = hd / 2.0
        contains = d <= 1.1 * rad
        if contains.any():
            idx = np.where(contains)[0]
            j = idx[np.argmin(hd[idx])]  # smallest containing crater
            rows.append({'field_id': int(f['field_id']), 'field_lat': round(float(f['lat']), 4),
                         'field_lon': round(float(f['lon']), 4), 'n_hollows': int(f['n']),
                         'n_nontrace': int(f['n_nontrace']),
                         'crater_lat': round(float(hlat[j]), 4), 'crater_lon': round(float(hlon[j]), 4),
                         'crater_D_km': round(float(hd[j]), 2),
                         'crater_name': (hn[j] if isinstance(hn[j], str) and hn[j] not in ('', 'nan', ' ') else '')})
        else:
            rows.append({'field_id': int(f['field_id']), 'field_lat': round(float(f['lat']), 4),
                         'field_lon': round(float(f['lon']), 4), 'n_hollows': int(f['n']),
                         'n_nontrace': int(f['n_nontrace']),
                         'crater_lat': None, 'crater_lon': None, 'crater_D_km': None, 'crater_name': ''})
    craters = pd.DataFrame(rows).sort_values('n_hollows', ascending=False)
    craters.to_csv(HERE + '/province_craters.csv', index=False)
    named = sorted({r['crater_name'] for r in rows if r['crater_name']})

    # ---- record ----
    record = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'province_box': {'lat': PLAT, 'lon_E': PLON},
        'n_fields_province': int(n_prov), 'n_fields_40_70N': n_4070,
        'n_hollows_province_fields': int(fp['n'].sum()),
        'comparison': comparison,
        'excess_after_reweighting': excess,
        'named_craters': named,
        'n_host_craters_found': int(craters['crater_D_km'].notna().sum()),
        'usgs_gazetteer': 'planetarynames.wr.usgs.gov not on sandbox allowlist; '
                          'names taken from the Herrick et al. 2018 >10 km catalogue name column '
                          '(IAU/USGS-sourced). 5-10 km craters have no names.',
        'lrm_flag_definition': 'flag = LRM brightness > 0 (0-255 scale); fraction over covered cells',
        'roughness_definition': 'std of 8 ppd DEM in 1-deg (8x8) windows',
        'bootstrap': f'{NBOOT} resamples of 5x5-deg blocks; area(cos-lat)-weighted means',
        'inputs': ART,
    }
    json.dump(record, open(HERE + '/province_record.json', 'w'), indent=1)
    # stash arrays for figure
    np.savez(HERE + '/_figdata.npz',
             lam=LAM, KL10=KL10, latc=latc, lonc=lonc,
             field_lat=fld['lat'].values, field_lon=fld['lon'].values, field_n=fld['n'].values,
             cr_lat=craters['crater_lat'].fillna(-999).values,
             cr_lon=craters['crater_lon'].fillna(-999).values,
             cr_D=craters['crater_D_km'].fillna(0).values,
             cr_nh=craters['n_hollows'].values,
             covered=covered.astype(np.int8))
    print('n_prov_fields', n_prov, 'of 40-70N', n_4070)
    print('host craters found', record['n_host_craters_found'], 'named', len(named))
    print('obs/null base %.2f lrm %.2f crater %.2f both %.2f' %
          (ratio_base, ratio_lrm, ratio_cr, ratio_both))
    print('named:', named)
    print('DONE')

if __name__ == '__main__':
    main()
