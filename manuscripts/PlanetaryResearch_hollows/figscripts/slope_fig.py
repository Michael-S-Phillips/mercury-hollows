import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Figure 3 of the Planetary Research paper: hollow occurrence against regional slope within host craters.
Reads occ_analysis_record.json (slope_test/occ_analysis.py on the Mac; 1.33 km baseline, 's1')."""
import json, os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
REC = os.environ.get("OCC_REC", _REPO + "/slope_test/occ_analysis_record.json")
OUT = os.environ.get("PR_FIGS", _REPO + "/manuscripts/PlanetaryResearch_hollows/figures")
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.linewidth": .6, "axes.edgecolor": "#5a5a57",
    "xtick.color": "#5a5a57", "ytick.color": "#5a5a57", "xtick.major.width": .6, "ytick.major.width": .6, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.labelcolor": "#1a1a19", "text.color": "#1a1a19", "axes.labelsize": 7.5})
INK, MUTED, GRID = "#1a1a19", "#6b6b67", "#dcdcd8"
BLUE, ORANGE, GREEN = "#1F6FB2", "#D4541E", "#178f66"
R = json.load(open(REC))["s1"]
def tidy(a, letter):
    a.grid(True, axis="y", color=GRID, lw=.5, zorder=0); a.set_axisbelow(True)
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
    a.text(-0.02, 1.03, letter, transform=a.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(5.9, 2.75))
bins = [0, 5, 10, 15, 20]; xc = np.array(bins) + 2.5
def ser(d, keys): return [d[str(k)] for k in keys]
for key, mk, col, mfc, off, lab in [("all_hollows", "o", INK, INK, -0.5, "hollows"), ("all_fields", "s", MUTED, "white", 0.5, "hollow fields")]:
    v = ser(R[key], bins); y = np.array([x["ratio"] for x in v]); lo = np.array([x["lo"] for x in v]); hi = np.array([x["hi"] for x in v])
    a1.errorbar(xc + off, y, yerr=[y - lo, hi - y], fmt=mk, ms=4.4, color=col, mfc=mfc, mew=1.1, lw=.9, capsize=2, zorder=5, label=lab)
for x, v in zip(xc, ser(R["all_hollows"], bins)):
    a1.text(x, -0.13, f"{v['A_km2'] / 1e3:.0f}" if v["A_km2"] >= 1e4 else f"{v['A_km2'] / 1e3:.1f}", ha="center", va="bottom", fontsize=6, color=MUTED)
a1.text(-0.6, -0.04, "area", ha="left", va="bottom", fontsize=6, color=MUTED)
a1.axhline(1, color=MUTED, lw=.9, ls=(0, (3, 2)), zorder=1)
a1.set_xlim(-1, 25.5); a1.set_ylim(-0.16, 1.85); a1.set_xticks([0, 5, 10, 15, 20, 25])
a1.set_xlabel("regional slope (°)"); a1.set_ylabel("observed / expected")
a1.legend(frameon=False, fontsize=6.4, loc="upper right", handletextpad=.3, borderpad=.2)
tidy(a1, "a")
bb = [5, 10, 15, 20]; xb = np.array(bb) + 2.5
for asp, mk, col, off, lab in [("eq", "o", ORANGE, -0.7, "equator-facing"), ("ew", "s", GREEN, 0.0, "east- or west-facing"), ("pole", "D", BLUE, 0.7, "pole-facing")]:
    v = [R["by_aspect"][f"('{asp}', {b})"] for b in bb]
    y = np.array([x["ratio"] for x in v]); lo = np.array([x["lo"] for x in v]); hi = np.array([x["hi"] for x in v])
    a2.errorbar(xb + off, y, yerr=[y - lo, hi - y], fmt=mk, ms=4.4, color=col, mfc=col, mew=1, lw=.9, capsize=2, zorder=5, label=lab)
a2.axhline(1, color=MUTED, lw=.9, ls=(0, (3, 2)), zorder=1)
a2.set_xlim(4, 25.5); a2.set_ylim(-0.05, 2.3); a2.set_xticks([5, 10, 15, 20, 25])
a2.set_xlabel("regional slope (°)"); a2.set_ylabel("observed / expected")
a2.legend(frameon=False, fontsize=6.4, loc="upper left", handletextpad=.3, borderpad=.2)
tidy(a2, "b")
fig.subplots_adjust(left=.085, right=.985, top=.9, bottom=.17, wspace=.32)
import figure_data; figure_data.dump(fig, "fig3_slope", f"{OUT}/data")  # archive: table of plotted values
fig.savefig(f"{OUT}/fig3_slope.png", dpi=400); fig.savefig(f"{OUT}/fig3_slope.pdf"); print("wrote fig3_slope")
