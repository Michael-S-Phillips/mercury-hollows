"""Figure S1 v4 (10 Nov 2026, Michael): the right-hand titles were cut off, so they read "Recomputed coordinates", and the count of
detections in other images moves into the lower right of each right-hand panel (white box, as the scale bar). Otherwise as v3.
Figure S1 v3 (4 Nov 2026, per Michael: show each frame's own detections, not detections from other images placed on blank
parts of the raster). Three MDIS NAC frames with many detections of their own. Left: each detection in the frame at its released
latitude and longitude (Bickel et al. 2025 release), joined by a line to its recomputed position. Right: the same detections at the
coordinates computed from their image pixel coordinates with the MDIS camera model in ISIS (GIS/vector/bickel2024_hollows_corrected.gpkg),
plus the recomputed positions of detections made in other images that fall on this frame.
Map projection: the CDR I/F frame (PDS MSGRMDS_2001, downloaded to data/mdis_frames/IMG) is resampled onto an equirectangular grid
with true scale at the frame-center latitude, using a projective transform from image pixel coordinates to map coordinates fitted
to the camera-model ground points of all detections in the frame (0.7-confidence release). Validation: for EN1035674191M and
EN1017816487M, projected the same way and compared with their ISIS cam2map GeoTIFFs (data/mdis_frames/geotiff_isis), the
phase-correlation offset is 0 pixels and the correlation is 0.92 and 0.99. Writes figures/figS1_positions_v4.png/.pdf and
analysis_v3/figS_positions_v4_record.json. Detections plotted are those of detections_qc.csv (confidence >= 0.75, incidence < 85 deg)."""
import os, json, sqlite3, numpy as np, pandas as pd, rasterio, matplotlib as mpl, matplotlib.pyplot as plt
HE = os.path.expanduser("~/mnt/hollows-emissivity"); PR = f"{HE}/manuscripts/PlanetaryResearch_hollows"; R = 2439700.0
FRAMES = ["CN0223615469M", "CN1002981281M", "CN1017701438M"]
def read_cdr(cn):
    raw = open(f"{HE}/data/mdis_frames/IMG/{cn}_IF_5.IMG", "rb").read()
    a = np.frombuffer(raw, dtype=">f4", offset=3 * 4096, count=1024 * 1024).reshape(1024, 1024).astype(float)
    a[(a < -1e30) | ~np.isfinite(a) | (a <= 0)] = np.nan; return a
def ground_points(cn):
    c = sqlite3.connect(f"{HE}/GIS/vector/bickel2024_hollows_corrected.gpkg")
    q = pd.read_sql("select img_x,img_y,frame_samples fs,frame_lines fl,lat,lon180 from hollows_points_corrected where mdis_id=?", c, params=(cn,))
    q["s"] = q.img_x * q.fs; q["l"] = q.img_y * q.fl; return q
class LocalMap:
    def __init__(self, lat0, lon0): self.lat0, self.lon0, self.k = lat0, lon0, np.cos(np.radians(lat0))
    def fwd(self, lat, lon):
        dl = ((np.asarray(lon, float) - self.lon0 + 180) % 360) - 180
        return R * np.radians(dl) * self.k, R * np.radians(np.asarray(lat, float) - self.lat0)
def fit_h(s, l, x, y):
    def norm(u, v):
        mu, mv = u.mean(), v.mean(); sc = np.sqrt(2) / np.mean(np.hypot(u - mu, v - mv))
        return np.array([[sc, 0, -sc * mu], [0, sc, -sc * mv], [0, 0, 1]])
    T1, T2 = norm(s, l), norm(x, y); p = T1 @ np.vstack([s, l, np.ones_like(s)]); q = T2 @ np.vstack([x, y, np.ones_like(x)]); A = []
    for i in range(p.shape[1]):
        X, Y, W = p[:, i]; u, v, w = q[:, i]
        A += [[0, 0, 0, -w * X, -w * Y, -w * W, v * X, v * Y, v * W], [w * X, w * Y, w * W, 0, 0, 0, -u * X, -u * Y, -u * W]]
    H = np.linalg.svd(np.array(A))[2][-1].reshape(3, 3); H = np.linalg.inv(T2) @ H @ T1; return H / H[2, 2]
def apply_h(H, u, v):
    u = np.atleast_1d(np.asarray(u, float)); v = np.atleast_1d(np.asarray(v, float)); p = H @ np.vstack([u, v, np.ones(len(u))]); return p[0] / p[2], p[1] / p[2]
def warp(img, H, x0, x1, y0, y1, res):
    xs = np.arange(x0 + res / 2, x1, res); ys = np.arange(y1 - res / 2, y0, -res); X, Y = np.meshgrid(xs, ys)
    s, l = apply_h(np.linalg.inv(H), X.ravel(), Y.ravel()); s = s.reshape(X.shape); l = l.reshape(X.shape)
    s0 = np.floor(s).astype(int); l0 = np.floor(l).astype(int); fs = s - s0; fl = l - l0
    ok = (s0 >= 0) & (s0 < 1023) & (l0 >= 0) & (l0 < 1023); out = np.full(X.shape, np.nan)
    a, b, c_, d = img[l0[ok], s0[ok]], img[l0[ok], s0[ok] + 1], img[l0[ok] + 1, s0[ok]], img[l0[ok] + 1, s0[ok] + 1]
    out[ok] = (a * (1 - fs[ok]) + b * fs[ok]) * (1 - fl[ok]) + (c_ * (1 - fs[ok]) + d * fs[ok]) * fl[ok]; return out
def validate():
    rec = {}
    for cn in ["CN1035674191M", "CN1017816487M"]:
        g = ground_points(cn); x = R * np.radians(((g.lon180.values + 180) % 360) - 180); y = R * np.radians(g.lat.values); H = fit_h(g.s.values, g.l.values, x, y)
        with rasterio.open(f"{HE}/data/mdis_frames/geotiff_isis/{cn.replace('CN', 'EN')}.tif") as ds:
            ref = ds.read(1).astype(float); b = ds.bounds; res = ds.res[0]
        ref[(ref <= -1e30) | ~np.isfinite(ref)] = np.nan
        mine = warp(read_cdr(cn), H, b.left, b.right, b.bottom, b.top, res)[:ref.shape[0], :ref.shape[1]]; m = np.isfinite(ref) & np.isfinite(mine)
        A = np.where(m, ref - ref[m].mean(), 0); B = np.where(m, mine - mine[m].mean(), 0); F = np.fft.fft2(A) * np.conj(np.fft.fft2(B))
        cc = np.fft.ifft2(F / (abs(F) + 1e-12)).real; iy, ix = np.unravel_index(cc.argmax(), cc.shape)
        iy = iy - cc.shape[0] if iy > cc.shape[0] // 2 else iy; ix = ix - cc.shape[1] if ix > cc.shape[1] // 2 else ix
        rec[cn] = dict(n_points=len(g), corr=float(np.corrcoef(ref[m], mine[m])[0, 1]), shift_px=[int(ix), int(iy)], grid_m=float(res))
    return rec
if __name__ == "__main__":
    d = pd.read_csv(f"{HE}/code/rerun/out_v2/detections_qc.csv"); fi = pd.read_csv(f"{HE}/hollow-networks/data/mdis_cdr_index_reduced.csv", usecols=["PRODUCT_ID", "INCIDENCE_ANGLE", "EMISSION_ANGLE", "HORIZONTAL_PIXEL_SCALE", "CENTER_LATITUDE", "CENTER_LONGITUDE"])
    fi = fi[fi.PRODUCT_ID.str.endswith("_IF_5")].assign(cn=lambda t: t.PRODUCT_ID.str[:13]).set_index("cn")
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": .6, "axes.edgecolor": "#5a5a57"})
    fig, axs = plt.subplots(3, 2, figsize=(5.9, 8.6)); REC = {"validation": validate(), "frames": {}}
    for row, cn in enumerate(FRAMES):
        g = ground_points(cn); lm = LocalMap(g.lat.median(), g.lon180.median()); gx, gy = lm.fwd(g.lat.values, g.lon180.values)
        H = fit_h(g.s.values, g.l.values, gx, gy); Hi = np.linalg.inv(H); px, py = apply_h(H, g.s.values, g.l.values); resid = np.hypot(px - gx, py - gy)
        cs, cl = apply_h(H, [0, 1023, 0, 1023, 511.5, 512.5], [0, 0, 1023, 1023, 511.5, 511.5]); res = float(np.hypot(cs[5] - cs[4], cl[5] - cl[4]))
        X0, X1, Y0, Y1 = cs[:4].min(), cs[:4].max(), cl[:4].min(), cl[:4].max(); img = warp(read_cdr(cn), H, X0, X1, Y0, Y1, res)
        own = d[d.image == cn]; oth = d[(d.image != cn) & (abs(d.lat - lm.lat0) < 2)]
        ox, oy = lm.fwd(own.lat.values, own.lon.values); rx, ry = lm.fwd(own.lat_release.values, own.lon_release.values)
        tx, ty = lm.fwd(oth.lat.values, oth.lon.values); ts, tl = apply_h(Hi, tx, ty); on_t = (ts >= 0) & (ts < 1024) & (tl >= 0) & (tl < 1024)
        rs, rl = apply_h(Hi, rx, ry); on_r = (rs >= 0) & (rs < 1024) & (rl >= 0) & (rl < 1024)
        allx = np.r_[ox, rx]; ally = np.r_[oy, ry]; m = 3000.0
        x0, x1, y0, y1 = allx.min() - m, allx.max() + m, ally.min() - m, ally.max() + m; w = max(x1 - x0, y1 - y0); cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        x0, x1, y0, y1 = cx - w / 2, cx + w / 2, cy - w / 2, cy + w / 2
        lo_, hi_ = np.nanpercentile(img, [1, 99.5])
        for col in (0, 1):
            ax = axs[row, col]; ax.set_facecolor("#e6e6e3")
            ax.imshow(img, cmap="gray", vmin=lo_, vmax=hi_, extent=[X0, X1, Y0, Y1], origin="upper", interpolation="nearest")
            if col == 0:
                for a, b_, c_, d_ in zip(ox, oy, rx, ry): ax.plot([a, c_], [b_, d_], color="#D4541E", lw=.6, alpha=.9)
                ax.plot(rx, ry, "o", ms=3.8, mfc="none", mec="#1a1a19", mew=1.9, ls="none"); ax.plot(rx, ry, "o", ms=3.8, mfc="none", mec="#D4541E", mew=.9, ls="none")
                ax.set_title("Released coordinates", fontsize=6.6, loc="left")
            else:
                ax.plot(tx[on_t], ty[on_t], "+", ms=3.6, mec="#1a1a19", mew=1.6, ls="none"); ax.plot(tx[on_t], ty[on_t], "+", ms=3.6, mec="#41c7e0", mew=.8, ls="none")
                ax.plot(ox, oy, "o", ms=3.8, mfc="none", mec="#1a1a19", mew=1.9, ls="none"); ax.plot(ox, oy, "o", ms=3.8, mfc="none", mec="#ffd84a", mew=.9, ls="none")
                ax.set_title("Recomputed coordinates", fontsize=6.6, loc="left")
                ax.text(.97, .04, f"+  {on_t.sum()} detections in {oth[on_t].image.nunique()} other images", transform=ax.transAxes, color="#1a1a19", ha="right", va="bottom", fontsize=6.2, zorder=7,
                        bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none", alpha=.85))   # v4: count inside the panel
            ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            sb = 5000.0; bx, by = x0 + .06 * w, y0 + .06 * w
            ax.plot([bx, bx + sb], [by, by], color="#1a1a19", lw=2.2, solid_capstyle="butt", zorder=6)
            ax.text(bx + sb / 2, by + .02 * w, "5 km", color="#1a1a19", ha="center", va="bottom", fontsize=6.4, zorder=6, bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=.85))
        r = fi.loc[cn]; lon = r.CENTER_LONGITUDE % 360
        axs[row, 0].text(-.03, .5, f"{cn}\n{r.HORIZONTAL_PIXEL_SCALE:.0f} m/pixel, incidence {r.INCIDENCE_ANGLE:.0f}°, emission {r.EMISSION_ANGLE:.0f}°\n{abs(r.CENTER_LATITUDE):.1f}°{'N' if r.CENTER_LATITUDE >= 0 else 'S'}, {lon:.1f}°E; {len(own)} detections, median difference {own.shift_km.median():.1f} km",
                         transform=axs[row, 0].transAxes, rotation=90, ha="right", va="center", fontsize=6.2)
        REC["frames"][cn] = dict(n_points=len(g), rms_m=float(np.sqrt((resid ** 2).mean())), max_m=float(resid.max()), px_m=res, n_own=len(own), own_released_on_frame=int(on_r.sum()),
                                 median_offset_km=float(own.shift_km.median()), n_other_on_frame=int(on_t.sum()), n_other_images=int(oth[on_t].image.nunique()))
        print(cn, REC["frames"][cn])
    for k, a in zip("abcdef", axs.ravel()): a.text(.02, .97, k, transform=a.transAxes, color="#1a1a19", fontsize=10, fontweight="bold", va="top", zorder=7, bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=.9))
    fig.subplots_adjust(left=.08, right=.99, top=.97, bottom=.01, wspace=.04, hspace=.12)
    fig.savefig(f"{PR}/figures/figS1_positions_v4.png", dpi=400); fig.savefig(f"{PR}/figures/figS1_positions_v4.pdf")
    json.dump(REC, open(f"{PR}/analysis_v3/figS_positions_v4_record.json", "w"), indent=1); print(REC["validation"])
