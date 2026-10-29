import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Supplementary figure: hollow detections of Bickel et al. (2025) at released and recomputed positions over three
map-projected MDIS NAC frames (ISIS 8.0 cam2map with reconstructed SPICE, uncontrolled; data/mdis_frames).
Detections are those of detections_qc.csv (confidence >= 0.75, incidence < 85 deg) from every source image."""
import os, numpy as np, pandas as pd, rasterio, matplotlib as mpl, matplotlib.pyplot as plt
HE = os.path.expanduser(_REPO + ""); OUT = f"{HE}/manuscripts/PlanetaryResearch_hollows/figures"
R = 2439700.0
d = pd.read_csv(f"{HE}/code/rerun/out_v2/detections_qc.csv")
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": .6, "axes.edgecolor": "#5a5a57",
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5})
frames = [("EN1017816487M", "CN1017816487M"), ("EN1035674191M", "CN1035674191M"), ("EN0259178827M", "CN0259178827M")]
fi = pd.read_csv(f"{HE}/data/mdis_frames/frames_index.csv").set_index("edr")
def xy(lat, lon):
    l = ((np.asarray(lon) + 180) % 360) - 180
    return R * np.radians(l), R * np.radians(np.asarray(lat))
fig, axs = plt.subplots(3, 2, figsize=(5.9, 8.4))
for row, (edr, cn) in enumerate(frames):
    with rasterio.open(f"{HE}/data/mdis_frames/geotiff_isis/{edr}.tif") as ds:
        img = ds.read(1).astype(float); b = ds.bounds
    img[~np.isfinite(img) | (img <= -1e30)] = np.nan
    lo_, hi_ = np.nanpercentile(img, [1, 99.5])
    xr, yr = xy(d.lat.values, d.lon.values); xo, yo = xy(d.lat_release.values, d.lon_release.values)
    inn = (xr > b.left) & (xr < b.right) & (yr > b.bottom) & (yr < b.top)
    ino = (xo > b.left) & (xo < b.right) & (yo > b.bottom) & (yo < b.top)
    sel = inn | ino
    s = d[sel]; nimg = s.image.nunique()
    # crop to the recomputed detections, with a 4 km margin, kept inside the frame
    m = 4000.0
    x0, x1 = max(b.left, xr[inn].min() - m), min(b.right, xr[inn].max() + m)
    y0, y1 = max(b.bottom, yr[inn].min() - m), min(b.top, yr[inn].max() + m)
    # square-ish window
    w = max(x1 - x0, y1 - y0); cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    x0, x1, y0, y1 = cx - w / 2, cx + w / 2, cy - w / 2, cy + w / 2
    for col, (X, Y, lab, mc) in enumerate([(xo[sel], yo[sel], "released positions", "#D4541E"), (xr[sel], yr[sel], "recomputed positions", "#ffd84a")]):
        ax = axs[row, col]
        ax.imshow(img, cmap="gray", vmin=lo_, vmax=hi_, extent=[b.left, b.right, b.bottom, b.top], origin="upper", interpolation="nearest")
        ax.plot(X, Y, "o", ms=3.8, mfc="none", mec="#1a1a19", mew=1.9, ls="none")
        ax.plot(X, Y, "o", ms=3.8, mfc="none", mec=mc, mew=.9, ls="none")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        k = int(((X > x0) & (X < x1) & (Y > y0) & (Y < y1)).sum())
        ax.set_title(f"{lab}: {k} of {len(s)} detections in view", fontsize=6.6, loc="left")
        # scale bar 5 km
        sb = 5000.0; bx, by = x0 + 0.06 * w, y0 + 0.06 * w
        ax.plot([bx, bx + sb], [by, by], color="#1a1a19", lw=2.2, solid_capstyle="butt", zorder=6)
        ax.text(bx + sb / 2, by + 0.02 * w, "5 km", color="#1a1a19", ha="center", va="bottom", fontsize=6.4, zorder=6,
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=.85))
    r = fi.loc[edr]
    axs[row, 0].text(-0.03, 0.5, f"{cn}\n{r.m_per_px:.0f} m/pixel, incidence {r.incidence:.0f}°\n{abs(r.center_lat):.1f}°{'N' if r.center_lat >= 0 else 'S'}, {r.center_lon_E:.1f}°E\n{len(s)} detections from {nimg} images",
                     transform=axs[row, 0].transAxes, rotation=90, ha="right", va="center", fontsize=6.3)
    print(edr, len(s), nimg)
for k, a in zip("abcdef", axs.ravel()): a.text(0.02, 0.97, k, transform=a.transAxes, color="#1a1a19", fontsize=10, fontweight="bold", va="top", zorder=7,
                                                    bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=.9))
fig.subplots_adjust(left=.08, right=.99, top=.97, bottom=.01, wspace=.04, hspace=.12)
fig.savefig(f"{OUT}/figS1_positions.png", dpi=400); fig.savefig(f"{OUT}/figS1_positions.pdf"); print("wrote")
