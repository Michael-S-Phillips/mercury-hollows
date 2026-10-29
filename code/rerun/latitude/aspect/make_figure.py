#!/usr/bin/env python3
"""Aspect figure: rose diagrams (facing aspect) by hemisphere + fraction-poleward
by |lat| band. Reads aspect_wall_hollows.csv and aspect_record.json."""
import json, numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

# Publication style (equivalent to figure-style apply_figure_style; standalone-safe)
def apply_figure_style(sizes=(9, 8, 7)):
    base, mid, small = sizes
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300,
        "font.family": "DejaVu Sans", "font.size": mid,
        "axes.titlesize": base, "axes.labelsize": base,
        "xtick.labelsize": small, "ytick.labelsize": small,
        "legend.fontsize": small, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "figure.facecolor": "white", "savefig.facecolor": "white",
    })
def panel_letter(ax, letter, case="lower"):
    s = letter.upper() if case == "upper" else letter.lower()
    ax.annotate(s, xy=(0, 1.0), xytext=(-24, 10), xycoords="axes fraction",
                textcoords="offset points", fontsize=11, fontweight="bold", va="bottom")

apply_figure_style(sizes=(9, 8, 7))

w = pd.read_csv("aspect_wall_hollows.csv")
rec = json.load(open("aspect_record.json"))

fig = plt.figure(figsize=(11, 4.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.35], wspace=0.35)

def rose(ax, ang, eq_dir, title, thomas=None):
    nb = 16
    edges = np.linspace(0, 360, nb+1)
    counts, _ = np.histogram(ang % 360, bins=edges)
    frac = counts / max(counts.sum(), 1)
    theta = np.radians(edges[:-1] + np.diff(edges)/2)
    width = np.radians(360/nb)
    ax.bar(theta, frac, width=width, bottom=0.0, color="#4C72B0",
           edgecolor="white", linewidth=0.5, align="center")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_xticks(np.radians([0, 90, 180, 270]))
    ax.set_xticklabels(["N (0\u00b0)", "E (90\u00b0)", "S (180\u00b0)", "W (270\u00b0)"])
    ax.set_yticklabels([])
    rmax = frac.max()*1.15
    ax.set_ylim(0, rmax)
    # equator-facing direction (green) and Thomas (orange) as radial rays
    ax.plot([np.radians(eq_dir)]*2, [0, rmax], color="#2CA02C", lw=2.2, zorder=5,
            label="equator-facing")
    if thomas is not None:
        ax.plot([np.radians(thomas)]*2, [0, rmax], color="#D55E00", lw=2.0,
                ls="--", zorder=6, label="Thomas 2014 (187\u00b0)")
    ax.set_title(title, pad=14)

# Panel A: North roses
wn = w[w.lat >= 0]; ws_ = w[w.lat < 0]
axA = fig.add_subplot(gs[0, 0], projection="polar")
rose(axA, wn["facing_aspect_deg"].values, 180, f"Northern wall hollows\n(n={len(wn)}, facing aspect)", thomas=187)
axB = fig.add_subplot(gs[0, 1], projection="polar")
rose(axB, ws_["facing_aspect_deg"].values, 0, f"Southern wall hollows\n(n={len(ws_)}, facing aspect)")

# Panel C: fraction poleward by band (n>=20 only; small-n bands uninformative)
NMIN = 20
axC = fig.add_subplot(gs[0, 2])
bands = [f"{e}-{e+15}" for e in range(0, 90, 15)]
xc = np.arange(len(bands))
for h, col, off, mk in [("N", "#4C72B0", -0.12, "o"), ("S", "#C44E52", 0.12, "s")]:
    fp, lo, hi, ns = [], [], [], []
    for b in bands:
        d = rec["wall_by_band"][h].get(b, {})
        if d.get("n", 0) >= NMIN:
            fp.append(d["frac_poleward"]); lo.append(d["frac_poleward_ci95"][0])
            hi.append(d["frac_poleward_ci95"][1]); ns.append(d["n"])
        else:
            fp.append(np.nan); lo.append(np.nan); hi.append(np.nan); ns.append(0)
    fp = np.array(fp); lo = np.array(lo); hi = np.array(hi)
    yerr = np.vstack([fp-lo, hi-fp])
    m = ~np.isnan(fp)
    axC.errorbar(xc[m]+off, fp[m], yerr=yerr[:, m], fmt=mk, color=col, ms=6,
                 capsize=3, lw=1.4, label=f"{h} hemisphere")
    for xi, yi, nn in zip(xc[m]+off, lo[m], np.array(ns)[m]):
        axC.annotate(f"n={nn}", (xi, yi-0.012), ha="center", va="top",
                     fontsize=6, color=col)
axC.axhline(0.5, color="0.4", ls=":", lw=1.2)
axC.annotate("null 0.5", (5.35, 0.505), va="bottom", ha="right", fontsize=6.5, color="0.4")
axC.set_xticks(xc); axC.set_xticklabels(bands, rotation=0)
axC.set_xlabel("|latitude| band (deg)")
axC.set_ylabel("fraction on poleward side\n(= equator-facing wall)")
axC.set_ylim(0.28, 0.86); axC.margins(x=0.06)
axC.set_title("Equator-facing preference strengthens poleward (N)", pad=18)
axC.legend(frameon=False, loc="lower right", fontsize=7)
axC.annotate("> 0.5 = equator-facing wall", (0.03, 0.98), xycoords="axes fraction",
             fontsize=6.5, color="0.3", va="top")

panel_letter(axA, "a"); panel_letter(axB, "b")
axC.annotate("c", xy=(0, 1.0), xytext=(-46, 20), xycoords="axes fraction",
             textcoords="offset points", fontsize=11, fontweight="bold", va="bottom")
# shared rose legend
axA.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.32), fontsize=6)

fig.suptitle("Scarp-aspect proxy (hollow position on host-crater wall; 665 m DEM cannot resolve hollow-scale slopes)",
             fontsize=8, y=1.02, color="0.25")
fig.savefig("aspect.png", dpi=300, bbox_inches="tight")

# 9.1 bbox overlap check
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
ov = [(a.get_text(), b.get_text()) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)]
print("text overlaps:", len(ov))
for o in ov[:12]: print("  ", o)
print("saved aspect.png")
