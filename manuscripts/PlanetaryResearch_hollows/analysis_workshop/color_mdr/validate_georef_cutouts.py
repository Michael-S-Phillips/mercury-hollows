"""Georeferencing validation (c): 750-nm cutouts of named craters with the cataloged rim drawn."""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdr_georef as g

CF = '../summary_model/crater_features.csv'
cr = pd.read_csv(CF)
named = [('Caloris', 31.9, 164.1, 1550), ('Rachmaninoff', 27.8, 57.3, 290), ('Tolstoj', -17.1, 195.3, 360),
         ('Raditladi', 27.3, 119.1, 263), ('N polar tile, host crater', 74.67, 141.57, 142), ('S polar tile crater', -76.05, 236.10, 198)]

def destination(lat0, lon0, az, dist):
    p0 = np.radians(lat0); c = dist / g.R_KM
    lat = np.arcsin(np.sin(p0) * np.cos(c) + np.cos(p0) * np.sin(c) * np.cos(az))
    lon = np.radians(lon0) + np.arctan2(np.sin(az) * np.sin(c) * np.cos(p0), np.cos(c) - np.sin(p0) * np.sin(lat))
    return np.degrees(lat), np.degrees(lon) % 360

G = g.Globe(band_idx=4)
fig, axs = plt.subplots(2, 3, figsize=(12, 8.4))
rec = []
for ax, (nm, la, lo, Dn) in zip(axs.ravel(), named):
    d = g.gc_dist_km(la, lo, cr.lat.values, cr.lon.values)
    cand = cr[(d < 150) & (cr.D > Dn / 2) & (cr.D < Dn * 2)]
    k = cand.index[np.argmin(d[cand.index])]
    c = cr.loc[k]
    R = c.D / 2
    half = 1.35 * R
    n = 500
    xs = np.linspace(-half, half, n)
    X, Y = np.meshgrid(xs, xs[::-1])
    rho = np.hypot(X, Y); az = np.arctan2(X, Y)
    lat, lon = destination(c.lat, c.lon, az, rho)
    img = G.sample(lat.ravel(), lon.ravel()).reshape(n, n)
    lo_, hi_ = np.nanpercentile(img, [1, 99])
    ax.imshow(img, cmap='gray', vmin=lo_, vmax=hi_, extent=[-half, half, -half, half])
    t = np.linspace(0, 2 * np.pi, 361)
    for f, ls, lw in [(1.0, '-', 1.4), (0.9, ':', 0.9)]:
        ax.plot(f * R * np.sin(t), f * R * np.cos(t), ls, color='#e8a33d', lw=lw)
    ax.plot(0, 0, '+', color='#e8a33d', ms=8)
    ax.set_title('%s  (id %d, D = %.0f km)\n%.2f°, %.2f°E' % (nm, c.id, c.D, c.lat, c.lon), fontsize=9)
    ax.set_xlabel('east (km)', fontsize=8); ax.set_ylabel('north (km)', fontsize=8); ax.tick_params(labelsize=7)
    rec.append(dict(name=nm, crater_id=int(c.id), lat=float(c.lat), lon=float(c.lon), D_km=float(c.D),
                    offset_from_named_km=float(d[k]), frac_valid=float(np.isfinite(img).mean())))
fig.suptitle('MDIS MDR 750-nm reflectance, azimuthal-equidistant cutouts; solid = cataloged rim (D/2), dotted = 0.9 R', fontsize=10)
fig.tight_layout()
fig.savefig('fig_georef_check.png', dpi=130)
json.dump(rec, open('georef_cutouts.json', 'w'), indent=1)
print(json.dumps(rec, indent=1))
