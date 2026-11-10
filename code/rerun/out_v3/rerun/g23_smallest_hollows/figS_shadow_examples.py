#!/usr/bin/env python3
"""g23 (9 Nov 2026), for Michael: supplement figure of automatically measured shadows (S9).
Top two rows: six hollows from the hand-measured set (<= 7 m/pixel; bench_depths.csv), one per log-size bin from 40 m to 1 km, each the
hollow whose |ln(hand/automatic)| is closest to the median over the set (typical agreement, not the best); measurements marked
uncertain, and automatic runs flagged at the edge of the crop, are excluded (the automatic set excludes them too). Bottom two rows: six accepted automatic measurements at 10-20 m/pixel (depth_all.csv; the range that holds most
automatic depths), drawn at random (seed 7) within size bins 150-300, 300-500 (2), 500-800 (2), and 800-1,500 m.
Each panel: the raw EDR crop at native resolution (detection box plus max(0.5 box size, 15 px) on each side, 0.5-99.5% stretch),
the detection box (yellow), the automatic shadow (red; shadow_measure.py v2 run, the longest run after the running median),
Michael's hand line where there is one (cyan, dashed), the direction the light travels, and a scale bar.
Writes figS_shadow_examples_v1.png (400 dpi) and .pdf, 15 cm wide, and figS_shadow_examples_v1.json (the hollows and depths)."""
import sys, os, json, numpy as np, pandas as pd, matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import small_shadow as G, depth_vs_size as DV
C, V, M = G.C, G.V, G.M
YEL, RED, CYAN, INK = "#f2c230", "#ff3b2f", "#38c6f4", "#1a1a19"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.4})
# ---- selection
b = pd.read_csv(f"{HERE}/bench_depths.csv"); b = b[(b.cls == "hollow") & (b.shadow == "measured") & b.auto_h.notna() & ~b.uncertain.astype(bool) & ~b.auto_flags.fillna("").str.contains("edge_sun|edge_far|frame_edge")].copy()
b["alr"] = np.log(b.depth_m / b.auto_h).abs(); med = b.alr.median()
edges = np.exp(np.linspace(np.log(40), np.log(1000), 7)); top = []
for lo, hi in zip(edges[:-1], edges[1:]):
    s = b[(b.D_m >= lo) & (b.D_m < hi)]
    if len(s): top.append(int(s.iloc[(s.alr - med).abs().argsort().iloc[0]].hollow_id))
a = pd.read_csv(f"{HERE}/depth_all.csv"); fl = a["flags"].fillna("")
a = a[(a.depth_m > 0) & ~fl.str.contains("no_shadow|edge_sun|frame_edge|no_edr") & a.px.between(10, 20)].copy()
rng = np.random.default_rng(7); bot = []
for lo, hi, k in [(150, 300, 1), (300, 500, 2), (500, 800, 2), (800, 1500, 1)]:
    s = a[(a.D_can_m > lo) & (a.D_can_m <= hi)].sort_values("hollow_id")
    bot += [int(x) for x in rng.choice(s.hollow_id.values, k, replace=False)]
# ---- detection boxes and image times
bk = C.load_bickel(with_corners=False)[["thumb", "num", "datetime", "x0", "x1", "y0", "y1", "lat", "lon180"]].drop_duplicates("thumb")
sel = DV.select(); sel = sel[(sel.tier == "le7")]
rows = []
for hid in top:
    r = sel[sel.hollow_id == hid].iloc[0]; rows.append(dict(hid=hid, thumb=r.thumb, set="hand"))
for hid in bot:
    r = a[a.hollow_id == hid].iloc[0]; rows.append(dict(hid=hid, thumb=r.thumb, set="auto"))
R = pd.DataFrame(rows).merge(bk, on="thumb", how="left"); R["num"] = R.num.astype("int64").astype(str).str.zfill(10)
hand = {int(d["hollow_id"]): d for d in json.load(open(f"{HERE}/bench_shadows_2026-11-09.json"))}
man = json.load(open(f"{HERE}/manifest_meas.json"))
V.fetch_many([(n, d) for n, d in R.drop_duplicates("num")[["num", "datetime"]].values])
def stretch(x):
    lo, hi = np.percentile(x, [0.5, 99.5]); return np.clip((x - lo) / (hi - lo + 1e-6), 0, 1)
fig, axs = plt.subplots(4, 3, figsize=(5.9, 8.2)); out = []
for ax, (_, r) in zip(axs.flat, R.iterrows()):
    p = C.edr_path(r.num); img, info = G.read_any(p); H, W = img.shape; fs = M.foreshortening(info)
    x0, x1, y0, y1 = r.x0 * W, r.x1 * W, r.y0 * H, r.y1 * H; bw = max(x1 - x0, y1 - y0); cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    half = bw / 2 + max(0.5 * bw, 15); X0, X1, Y0, Y1 = int(max(0, cx - half)), int(min(W, cx + half)), int(max(0, cy - half)), int(min(H, cy + half))
    ax.imshow(stretch(img[Y0:Y1, X0:X1].astype(float)), cmap="gray", interpolation="nearest", extent=(X0 - .5, X1 - .5, Y1 - .5, Y0 - .5))
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=YEL, lw=0.8))
    o = M.measure(img, info, (x0, x1, y0, y1), return_debug=True); dbg = o.pop("_dbg", None)
    inc = float(C.local_incidence(r.lat, r.lon180, info["subsol_lat"], info["subsol_lon"])); ti = np.tan(np.radians(inc))
    h_auto = o["L_px"] * info["px"] * fs / ti if np.isfinite(o["L_px"]) else np.nan
    if dbg is not None and np.isfinite(o.get("run_u0", np.nan)):
        du, dw = dbg["du"], dbg["dw"]; c0 = np.array([dbg["cx"], dbg["cy"]])
        pa = c0 + o["run_u0"] * du + o["run_w"] * dw; pb = c0 + o["run_u1"] * du + o["run_w"] * dw
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=RED, lw=1.3, solid_capstyle="butt", zorder=5)
    h_hand = np.nan
    if r.set == "hand" and r.hid in hand and all(k in hand[r.hid] for k in ("ax", "ay", "bx", "by")):
        m = man[str(r.hid)]; hh = hand[r.hid]; h_hand = hh.get("depth_m", np.nan)
        ax.plot([hh["ax"] + m["X0"], hh["bx"] + m["X0"]], [hh["ay"] + m["Y0"], hh["by"] + m["Y0"]], color=CYAN, lw=1.0, ls=(0, (2, 1.4)), zorder=6)
    az = np.radians(info["subsol_az"]); ux, uy = -np.cos(az), -np.sin(az)     # light-travel direction, image coordinates (lines down)
    span = X1 - X0; sx, sy = X0 + 0.12 * span, Y0 + 0.12 * (Y1 - Y0); L = 0.13 * span
    ax.annotate("", xy=(sx + ux * L / 2, sy + uy * L / 2), xytext=(sx - ux * L / 2, sy - uy * L / 2), arrowprops=dict(arrowstyle="-|>", color="white", lw=0.8, mutation_scale=6))
    wm = span * info["px"]; nice = [5, 10, 20, 50, 100, 200, 500]; sb = max([n for n in nice if n <= 0.3 * wm] or [nice[0]]); sbp = sb / info["px"]
    bx0, by0 = X1 - 0.06 * span - sbp, Y1 - 0.07 * (Y1 - Y0)
    ax.plot([bx0, bx0 + sbp], [by0, by0], color="white", lw=1.6, solid_capstyle="butt"); ax.text(bx0 + sbp / 2, by0 - 0.03 * (Y1 - Y0), f"{sb} m", color="white", fontsize=5.4, ha="center", va="bottom")
    Dm = man[str(r.hid)]["D"] if r.set == "hand" else float(a[a.hollow_id == r.hid].D_can_m.iloc[0])
    ax.set_title(f"{Dm:.0f} m across, {info['px']:.1f} m/pixel, i = {inc:.0f}$^\\circ$", fontsize=5.9, pad=2)
    lab = f"automatic {h_auto:.1f} m" + (f"   hand {h_hand:.1f} m" if np.isfinite(h_hand) else "")
    ax.text(0.03, 0.03, lab, transform=ax.transAxes, color="white", fontsize=5.6, va="bottom", bbox=dict(boxstyle="round,pad=0.2", fc="black", ec="none", alpha=0.55))
    ax.set_xlim(X0 - .5, X1 - .5); ax.set_ylim(Y1 - .5, Y0 - .5); ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_linewidth(0.4)
    out.append(dict(hollow_id=int(r.hid), set=r.set, image=f"EN{r.num}M", D_m=Dm, px=round(info["px"], 2), inc=round(inc, 1), depth_auto=round(h_auto, 2), depth_hand=None if not np.isfinite(h_hand) else round(h_hand, 2), flags=o["flags"]))
    os.remove(p) if os.path.exists(p) else None
for k, lab in [(0, "a"), (2, "b")]:
    axs[k, 0].text(-0.08, 1.06, lab, transform=axs[k, 0].transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")
fig.subplots_adjust(left=0.045, right=0.99, top=0.97, bottom=0.01, wspace=0.06, hspace=0.16)
fig.savefig(f"{HERE}/figS_shadow_examples_v1.png", dpi=400); fig.savefig(f"{HERE}/figS_shadow_examples_v1.pdf")
json.dump(dict(median_abs_ln_ratio=round(float(med), 3), panels=out), open(f"{HERE}/figS_shadow_examples_v1.json", "w"), indent=1)
print(json.dumps(out))
