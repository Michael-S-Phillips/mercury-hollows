"""Georeferencing validation (a) corners/seams and (b) polar vs equirectangular overlap."""
import json, numpy as np
import mdr_georef as g

rec = {}
# ---- (a) corner check for every equirectangular tile + label/header facts
corner = []
for n in g.all_tiles(include_2700=True):
    t = g.Tile(n)
    L = t.lbl
    row = dict(tile=n, proj=L['PROJ'], lines=t.H, samples=t.W, record_bytes_eq_samples_x4=L['RECORD_BYTES'] == t.W * 4,
               image_ptr=L['IMAGE_PTR'], header_records=0)
    if not t.polar:
        la, lo = t.pix2ll(0, 0)
        hp = 0.5 / 64 * 1.0002
        row.update(first_pixel_lat=float(la), first_pixel_lon=float(lo),
                   dlat_from_maxlat_px=float((L['MAXIMUM_LATITUDE'] - la) * t.R * np.pi / 180 / t.scale_km),
                   dlon_from_west_px=float(g.wrap180(lo - L['WESTERNMOST_LONGITUDE']) * np.pi / 180 * t.R * t.coslat_c / t.scale_km))
        la2, lo2 = t.pix2ll(t.H - 1, t.W - 1)
        row.update(last_pixel_lat=float(la2), last_pixel_lon=float(lo2))
    else:
        la, lo = t.pix2ll(0, 0)
        row.update(corner_lat=float(la), label_min_or_max_lat=L['MINIMUM_LATITUDE'] if t.north else L['MAXIMUM_LATITUDE'])
    corner.append(row)
rec['corners'] = corner
eq = [c for c in corner if 'dlat_from_maxlat_px' in c]
rec['corner_summary'] = dict(n_equirect=len(eq),
    max_abs_dlat_px=max(abs(c['dlat_from_maxlat_px'] - 0.5) for c in eq),
    max_abs_dlon_px=max(abs(c['dlon_from_west_px'] - 0.5) for c in eq),
    note='line 1 / sample 1 centers sit 0.5 px inside MAXIMUM_LATITUDE / WESTERNMOST_LONGITUDE for every tile (values are deviation from 0.5 px)')

# ---- (a2) seam continuity: profiles just either side of tile boundaries, correct vs flipped
G = g.Globe(band_idx=4)
rng = np.random.default_rng(1)
seams = {}
# N-S seams at lat 0, +-22.5, +-43.75 ; E-W seams at lon 0, 36, 45, 90
eps = 0.6 / 64
for name, lat0 in [('lat+43.75', 43.75), ('lat+22.5', 22.5), ('lat0', 0.0), ('lat-22.5', -22.5), ('lat-43.75', -43.75)]:
    lon = rng.uniform(0, 360, 4000)
    a = G.sample(np.full_like(lon, lat0 + eps), lon); b = G.sample(np.full_like(lon, lat0 - eps), lon)
    ok = np.isfinite(a) & np.isfinite(b)
    seams[name] = dict(n=int(ok.sum()), r_across=float(np.corrcoef(a[ok], b[ok])[0, 1]),
                       med_abs_rel_diff=float(np.median(np.abs(a[ok] - b[ok]) / b[ok])))
for name, lon0 in [('lon0', 0.0), ('lon36', 36.0), ('lon45', 45.0), ('lon90', 90.0), ('lon180', 180.0)]:
    lat = rng.uniform(-64, 64, 4000)
    a = G.sample(lat, np.full_like(lat, lon0 + 0.6 / 64 / np.cos(np.radians(lat)))); b = G.sample(lat, np.full_like(lat, lon0 - 0.6 / 64 / np.cos(np.radians(lat))))
    ok = np.isfinite(a) & np.isfinite(b)
    seams[name] = dict(n=int(ok.sum()), r_across=float(np.corrcoef(a[ok], b[ok])[0, 1]),
                       med_abs_rel_diff=float(np.median(np.abs(a[ok] - b[ok]) / b[ok])))
rec['seams_750nm'] = seams

# ---- (b) polar vs equirectangular in the overlap band 55-64.9 deg
def polar_sample(t, lat, lon, band=4, tfill=None):
    i, j = t.ll2pix(lat, lon)
    i = np.rint(i).astype(int); j = np.rint(j).astype(int)
    ok = (i >= 0) & (i < t.H) & (j >= 0) & (j < t.W)
    v = np.full(lat.shape, np.nan, np.float32)
    v[ok] = t.mm[band, i[ok], j[ok]]
    if tfill is not None:
        bad = ok & ~(v > 0)
        v[bad] = tfill.mm[band, i[bad], j[bad]]
    v[~(v > 0)] = np.nan
    return v

over = {}
for pole, tname, sgn in [('north', 'H01NP4', 1), ('south', 'H15SP4', -1)]:
    t = g.Tile(tname); tf = g.Tile('2700_H15SP4') if pole == 'south' else None
    lat = sgn * rng.uniform(55, 64.9, 2000); lon = rng.uniform(0, 360, 2000)
    e = G.sample(lat, lon, source='equi')
    res = {}
    for lab, lon_t in [('as_derived', lon), ('mirrored_lon', -lon), ('rot180', lon + 180), ('rot90', lon + 90)]:
        p = polar_sample(t, lat, lon_t, tfill=tf)
        ok = np.isfinite(p) & np.isfinite(e)
        idx = np.flatnonzero(ok)[:200]
        res[lab] = dict(n=int(len(idx)), pearson_r=float(np.corrcoef(p[idx], e[idx])[0, 1]),
                        median_abs_diff=float(np.median(np.abs(p[idx] - e[idx]))),
                        median_abs_rel_diff=float(np.median(np.abs(p[idx] - e[idx]) / e[idx])),
                        median_ratio=float(np.median(p[idx] / e[idx])))
    # all 8 bands for the derived orientation, 3x3-insensitive check
    p8 = np.stack([polar_sample(t, lat, lon, band=b, tfill=tf) for b in range(8)])
    e8 = np.stack([G.sample(lat, lon, band=b, source='equi') for b in range(8)])
    ok = np.all(np.isfinite(p8) & np.isfinite(e8), axis=0); idx = np.flatnonzero(ok)[:200]
    res['all_bands_as_derived'] = dict(n=int(len(idx)),
        median_abs_rel_diff_by_band=[float(np.median(np.abs(p8[b, idx] - e8[b, idx]) / e8[b, idx])) for b in range(8)],
        median_ratio_by_band=[float(np.median(p8[b, idx] / e8[b, idx])) for b in range(8)])
    over[pole] = res
# south: 2700 tile alone vs H15SP4
t = g.Tile('H15SP4'); t2 = g.Tile('2700_H15SP4')
lat = -rng.uniform(65, 85, 4000); lon = rng.uniform(0, 360, 4000)
a = polar_sample(t, lat, lon); b = polar_sample(t2, lat, lon); ok = np.isfinite(a) & np.isfinite(b)
over['south_2700_vs_H15SP4'] = dict(n=int(ok.sum()), pearson_r=float(np.corrcoef(a[ok], b[ok])[0, 1]),
    median_abs_rel_diff=float(np.median(np.abs(a[ok] - b[ok]) / b[ok])), frac_H15SP4_missing=float(np.mean(~np.isfinite(a))),
    frac_missing_after_fill=float(np.mean(~np.isfinite(a) & ~np.isfinite(b))))
rec['polar_overlap_750nm'] = over
json.dump(rec, open('georef_validation.json', 'w'), indent=1)
print(json.dumps({k: rec[k] for k in ['corner_summary', 'seams_750nm', 'polar_overlap_750nm']}, indent=1))
