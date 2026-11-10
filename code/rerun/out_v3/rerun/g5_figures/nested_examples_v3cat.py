# g5 new script (6 Nov 2026): no script for figures/fig_nested_examples.png (29 Oct 2026) was saved, so this script redraws its layout from the data: same two parents, same 166 m/pixel MDIS mosaic, same symbols and colors (sampled from the v2 PNG), on the v3 catalog (children's hosting from host_assignments_v3.csv, hollows = non-trace hollows of hollows_canonical_v3.csv). Arrows end at each child's rim (the v2 arrows ended at varying distances short of the children).
"""Figure 3 (fig_nested_examples) on the v3 catalog. Children = craters >= 5 km in the containment tests (kind crater, possible
peak rings dropped) whose parent in g1_graph_joint/v3/crater_dag_nodes.csv is Praxiteles (200309) or Lermontov (200362), the
definition g1's parents_and_examples.py uses for the caption counts. Host = at least one hollow assigned within 0.9 radii (h09).
Writes g5_figures/fig_nested_examples_v3cat.png (400 dpi) and .pdf, and nested_examples_v3cat_record.json."""
import os, json, numpy as np, pandas as pd, rasterio, matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
HE = os.path.expanduser("~/mnt/hollows-emissivity"); G5 = os.path.dirname(os.path.abspath(__file__))
G1 = f"{HE}/code/rerun/out_v3/rerun/g1_graph_joint"
BDR = f"{HE}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif"
RM = 2439.4; KD = np.pi * RM / 180; R0 = 2439400.0
GOLD, BLUE, HOL, DARK, INK = "#E9A91F", "#4F97EA", "#FFE24D", "#1a1a19", "#000000"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": .6})
lon180 = lambda x: ((x + 180) % 360) - 180
def crop(lat0, lon0, half_km, npx):   # same windowed read as analysis_workshop/caloris_v9/figS_caloris_map_v9.py (no antimeridian case needed here)
    lo0 = lon180(lon0); dla = half_km / KD; dlo = half_km / (KD * np.cos(np.radians(lat0)))
    with rasterio.open(BDR) as ds:
        win = from_bounds(np.radians(lo0 - dlo) * R0, np.radians(lat0 - dla) * R0, np.radians(lo0 + dlo) * R0, np.radians(lat0 + dla) * R0, ds.transform)
        a = ds.read(1, window=win, out_shape=(npx, npx), resampling=Resampling.average, boundless=True).astype(float)
    p0, p1 = np.percentile(a[a > 0], [1, 99.5]); return np.clip((a - p0) / (p1 - p0), 0, 1)
def local(lat, lon, lat0, lon0): return lon180(np.asarray(lon) - lon0) * KD * np.cos(np.radians(lat0)), (np.asarray(lat) - lat0) * KD
def ring(ax, x, y, r, **kw): t = np.linspace(0, 2 * np.pi, 240); ax.plot(x + r * np.cos(t), y + r * np.sin(t), **kw)

n = pd.read_csv(f"{G1}/v3/crater_dag_nodes.csv", dtype={"id": str, "parent_id": str}, low_memory=False)
ha = pd.read_csv(f"{HE}/code/rerun/out_v3/host_assignments_v3.csv", usecols=["host_uid", "rho"]).dropna(); ha["host_uid"] = ha.host_uid.astype(int).astype(str)
n["n09"] = n.id.map(ha[ha.rho <= 0.9].host_uid.value_counts()).fillna(0).astype(int); n["h09"] = (n.n09 > 0).astype(int)
s = n[(n.kind == "crater") & (~n.possible_ring.astype(bool))]
hc = pd.read_csv(f"{HE}/code/rerun/out_v3/hollows_canonical_v3.csv"); hc = hc[~hc.is_crater_trace.astype(bool)]   # v3 non-trace hollows (6,299)
PAR = [("Praxiteles", "200309"), ("Lermontov", "200362")]
fig = plt.figure(figsize=(7.2, 3.75)); REC = {}
for k, (name, pid) in enumerate(PAR):
    p = n.set_index("id").loc[pid]; R = p.D / 2; half = R + 30
    ax = fig.add_axes([0.034 + 0.5 * k, 0.090, 0.4325, 0.829])
    ax.imshow(crop(p.lat, p.lon, half, 1250), cmap="gray", extent=[-half, half, -half, half], origin="upper", interpolation="lanczos")
    ring(ax, 0, 0, R, color="white", lw=1.4, ls=(0, (5, 3)), zorder=3)
    ch = s[s.parent_id == pid].sort_values("D", ascending=False)
    cx, cy = local(ch.lat.values, ch.lon.values, p.lat, p.lon)
    for (x, y), r, h in zip(zip(cx, cy), ch.D.values / 2, ch.h09.values):
        d = np.hypot(x, y); stop = max(d - r - 1.0, 0.0)
        ax.annotate("", xy=(x * stop / d, y * stop / d), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color="white", lw=.8, mutation_scale=7, shrinkA=0, shrinkB=0), zorder=4)
        ring(ax, x, y, r, color=GOLD if h else BLUE, lw=1.6, zorder=5)
    v = hc[(np.abs(hc.lat - p.lat) < half / KD) & (np.abs(lon180(hc.lon - p.lon)) < half / (KD * np.cos(np.radians(p.lat))))]
    hx, hy = local(v.lat.values, v.lon.values, p.lat, p.lon)
    ax.scatter(hx, hy, s=5.5, c=HOL, edgecolors=DARK, linewidths=.35, zorder=6)
    ax.set_xlim(-half, half); ax.set_ylim(-half, half); ax.set_xticks([]); ax.set_yticks([])
    nh = int(ch.h09.sum()); ax.set_title(f"{name}, {p.D:.0f} km: {nh} of {len(ch)} craters inside have hollows", loc="left", fontsize=7.5, color=INK)
    ax.text(0.025, 0.965, f"{p.lat:.1f}°N, {p.lon % 360:.1f}°E", transform=ax.transAxes, color="white", fontsize=7, va="top")
    ax.plot([0.95 * half - 50, 0.95 * half], [-0.9 * half] * 2, color="white", lw=2, solid_capstyle="butt", zorder=7)
    ax.text(0.95 * half - 25, -0.86 * half, "50 km", color="white", ha="center", va="bottom", fontsize=7.5, zorder=7)
    REC[name] = dict(id=pid, D=float(p.D), lat=float(p.lat), lon=float(p.lon), half_km=float(half), hollows_in_view=int(len(v)), parent_n09=int(p.n09),
                     children=[dict(id=i, D=round(float(D), 2), lat=float(la), lon=float(lo), n09=int(c9), h09=int(h)) for i, D, la, lo, c9, h in
                               zip(ch.id, ch.D, ch.lat, ch.lon, ch.n09, ch.h09)], children_h09=nh)
REC["h09_rate_all_tests"] = dict(hosts=int(s.h09.sum()), craters=int(len(s)), rate=float(s.h09.mean()))
H = [Line2D([], [], color="#555555", lw=1.4, ls=(0, (5, 3)), label="parent rim (white on map)"),
     Line2D([], [], marker="o", ms=9, mfc="none", mec=GOLD, mew=1.6, ls="none", label="crater inside, with hollows"),
     Line2D([], [], marker="o", ms=9, mfc="none", mec=BLUE, mew=1.6, ls="none", label="crater inside, no hollows"),
     Line2D([], [], marker="o", ms=3.6, mfc=HOL, mec=DARK, mew=.4, ls="none", label="hollow")]
fig.legend(handles=H, loc="center", bbox_to_anchor=(0.5, 0.037), ncol=4, frameon=False, fontsize=7.5, handlelength=2.0, columnspacing=2.4)
fig.savefig(f"{G5}/fig_nested_examples_v3cat.png", dpi=400); fig.savefig(f"{G5}/fig_nested_examples_v3cat.pdf")
json.dump(REC, open(f"{G5}/nested_examples_v3cat_record.json", "w"), indent=1)
print(json.dumps({k: (v["children_h09"], len(v["children"]), v["hollows_in_view"]) if isinstance(v, dict) and "children" in v else v for k, v in REC.items()}))
