# g5 copy of manuscripts/PlanetaryResearch_hollows/analysis_v3/figS_terrain_examples.py (6 Nov 2026). Changed: host assignments = code/rerun/out_v3/host_assignments_v3.csv (was code/rerun/latitude/aspect/host_assignments.csv, v2); output g5_figures/figS2_terrain_examples_v3cat.png/.pdf.
"""Supplementary figure: terrain classes of Section 2.4 for four host craters, on a hillshade of the 665 m/pixel
USGS topography. Class definitions follow peak_wall_test/peak_wall.py exactly (250 m relief threshold, 0.45-radius peak
limit, 0.6 for peak rings; walls 0.55-1.05 radii, 0.7-1.05 for peak rings; floors within 0.8 radii and <= 250 m above
the floor elevation; superposed cataloged craters >= 5 km masked)."""
import os, numpy as np, pandas as pd, rasterio, matplotlib as mpl, matplotlib.pyplot as plt
from rasterio.windows import from_bounds
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
HE = os.path.expanduser("~/mnt/hollows-emissivity")
OUT = os.path.dirname(os.path.abspath(__file__))  # g5: output folder
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM = 2439.4; KD = np.pi * RM / 180
cr = pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"] = cr.lon % 360
ha = pd.read_csv(f"{HE}/code/rerun/out_v3/host_assignments_v3.csv")  # g5: v3 host assignments
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": .6, "axes.edgecolor": "#5a5a57",
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.labelsize": 7.5})
COL = {"peak": "#D4541E", "wall": "#1F6FB2", "floor": "#178f66", "masked": "#7a7a75"}
def gcd(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    return RM * np.arccos(np.clip(np.sin(p1)*np.sin(p2) + np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2 - lo1)), -1, 1))
def classes(c, ds, HT=250.0, RPK=0.45):
    R = c.diam / 2; la0, lo0 = c.lat, c.lon360; cl = np.cos(np.radians(la0))
    dla = 1.25 * R / KD; dlo = min(1.25 * R / (KD * max(cl, 0.05)), 179)
    xs = lambda lo: np.radians(lo - 180) * 2439400.0; ys = lambda la: np.radians(la) * 2439400.0
    win = from_bounds(xs(lo0 - dlo), ys(la0 - dla), xs(lo0 + dlo), ys(la0 + dla), ds.transform)
    z = ds.read(1, window=win, boundless=True, fill_value=-32768).astype(float); z[z == -32768] = np.nan; z *= 0.5
    T = ds.window_transform(win); nr, nc = z.shape
    xc = T.c + T.a * (np.arange(nc) + .5); yc = T.f + T.e * (np.arange(nr) + .5)
    LON = 180 + np.degrees(xc / 2439400.0)[None, :] * np.ones((nr, 1)); LAT = np.degrees(yc / 2439400.0)[:, None] * np.ones((1, nc))
    dx = (LON - lo0) * KD * cl; dy = (LAT - la0) * KD; rho = np.hypot(dx, dy) / R
    ring = str(c.cent_str).strip() in ("pr", "rpc", "pb")
    fz = z[(rho > 0.15) & (rho < 0.6)]; fz = fz[np.isfinite(fz)]; floor = np.percentile(fz, 20); h = z - floor
    peak = (rho <= (0.6 if ring else RPK)) & (h > HT)
    wall = (rho >= (0.7 if ring else 0.55)) & (rho <= 1.05) & (h > HT)
    flo = (rho <= 0.8) & (h <= HT)
    sub = cr[(cr.diam >= 5) & (cr.diam < 0.5 * c.diam) & (np.abs(cr.lat - la0) < dla)]
    sub = sub[gcd(la0, lo0, sub.lat.values, sub.lon360.values) < 1.05 * R]
    m = np.zeros_like(rho, bool)
    for s in sub.itertuples():
        m |= np.hypot((LON - s.lon360) * KD * np.cos(np.radians(s.lat)), (LAT - s.lat) * KD) < 1.1 * s.diam / 2
    lab = np.zeros(z.shape, int); lab[flo] = 1; lab[wall] = 2; lab[peak] = 3; lab[m & (rho <= 1.05)] = 4; lab[rho > 1.05] = 0
    return z, dx, dy, lab, R, ring
def hillshade(z, az=315, alt=35):
    gy, gx = np.gradient(np.nan_to_num(z, nan=np.nanmean(z)), 665.24, 665.24)
    sl = np.arctan(np.hypot(gx, gy)); asp = np.arctan2(-gx, gy)
    a, e = np.radians(az), np.radians(alt)
    return np.clip(np.sin(e) * np.cos(sl) + np.cos(e) * np.sin(sl) * np.cos(a - asp), 0, 1)
picks = [("Tyagaraja", None), ("Praxiteles", None), ("Scarlatti", None), ("Sander", None)]
fig, axs = plt.subplots(2, 2, figsize=(5.9, 6.0))
cmap = ListedColormap([(0, 0, 0, 0), COL["floor"], COL["wall"], COL["peak"], COL["masked"]])
with rasterio.open(DEM) as ds:
    for ax, (nm, _) in zip(axs.ravel(), picks):
        c = cr[cr.name == nm].iloc[0]
        z, dx, dy, lab, R, ring = classes(c, ds)
        ext = [dx.min(), dx.max(), dy.min(), dy.max()]
        hsd = hillshade(z); lo_, hi_ = np.percentile(hsd, [2, 98])
        ax.imshow(hsd, cmap="gray", extent=ext, origin="upper", vmin=lo_, vmax=hi_, interpolation="nearest")
        ax.imshow(np.ma.masked_equal(lab, 0), cmap=cmap, vmin=0, vmax=4, alpha=.30, extent=ext, origin="upper", interpolation="nearest")
        th = np.linspace(0, 2 * np.pi, 361)
        ax.plot(R * np.cos(th), R * np.sin(th), color="white", lw=.7, ls=(0, (3, 2)))
        hs = ha[ha.host_uid == c.unique_id]
        hx = ((hs.lon.values % 360 - c.lon360 + 180) % 360 - 180) * KD * np.cos(np.radians(c.lat)); hy = (hs.lat.values - c.lat) * KD
        ax.plot(hx, hy, "o", ms=2.3, mfc="#ffd84a", mec="#1a1a19", mew=.35, ls="none")
        ax.set_xlim(-1.2 * R, 1.2 * R); ax.set_ylim(-1.2 * R, 1.2 * R); ax.set_aspect("equal")
        ax.set_xlabel("east (km)"); ax.set_ylabel("north (km)")
        kind = "peak-ring limits" if ring else "central-peak limits"
        ax.set_title(f"{nm}, {c.diam:.0f} km, {abs(c.lat):.1f}°{'N' if c.lat >= 0 else 'S'}, {c.lon360:.1f}°E\n{kind}; {len(hs)} hollows assigned", fontsize=6.6, loc="left")
        ax.tick_params(length=2)
h = [Patch(fc=COL[k], alpha=.6, label=l) for k, l in [("peak", "peak"), ("wall", "wall"), ("floor", "floor"), ("masked", "superposed crater (masked)")]]
h += [Line2D([], [], marker="o", ms=3.5, mfc="#ffd84a", mec="#1a1a19", mew=.4, ls="none", label="hollow assigned to host"),
      Line2D([], [], color="#7a7a75", lw=.8, ls=(0, (3, 2)), label="cataloged rim"), Patch(fc="white", ec="#9a9a95", lw=.5, label="unclassified")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=6.6, bbox_to_anchor=(0.5, 0.0))
fig.subplots_adjust(left=.09, right=.99, top=.93, bottom=.14, wspace=.30, hspace=.38)
for k, a in zip("abcd", axs.ravel()): a.text(-0.2, 1.10, k, transform=a.transAxes, fontsize=10, fontweight="bold")
fig.savefig(f"{OUT}/figS2_terrain_examples_v3cat.png", dpi=400); fig.savefig(f"{OUT}/figS2_terrain_examples_v3cat.pdf")  # g5: output names
print("wrote")
