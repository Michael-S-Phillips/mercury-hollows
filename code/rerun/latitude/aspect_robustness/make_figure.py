#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
"""Aspect robustness figure: 3 panels.
A: frac_poleward vs latitude band, all-hosts vs D20-100 hosts (per-crater + per-hollow).
B: 15-30N decomposition (Caloris / host D / longitude), per-hollow vs per-crater.
C: DEM slope-azimuth vs position facing agreement.
"""
import json, numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
# figure-style helpers (apply_figure_style, panel_letter) are provided by the
# figure-style skill kernel plugin. When present in globals (run in-kernel) they
# are used; otherwise fall back to no-op styling so the script still runs.
if "apply_figure_style" not in globals():
    try:
        from figure_style_kernel import apply_figure_style, panel_letter
    except Exception:
        def apply_figure_style(*a, **k):
            import matplotlib.pyplot as _p; _p.rcParams.update({"font.size": 8})
        def panel_letter(ax, L, **k):
            ax.text(-0.12, 1.05, L, transform=ax.transAxes, fontweight="bold", fontsize=11)

BASE = _REPO + ""
FOLD = f"{BASE}/code/rerun/latitude/aspect_robustness"
rec = json.load(open(f"{FOLD}/aspect_robustness_record.json"))
wall = pd.read_csv(f"{FOLD}/wall_with_dem_aspect.csv")

apply_figure_style(sizes=(8, 7, 6))

# recorded (all-hosts) N-band per-hollow & per-crater from aspect_record.json
ar = json.load(open(f"{BASE}/code/rerun/latitude/aspect/aspect_record.json"))
bands = ["0-15", "15-30", "30-45", "45-60", "60-75"]
x = np.arange(len(bands))

def get(src, b):
    return src[b]

all_ph = [ar["wall_by_band"]["N"][b]["frac_poleward"] for b in bands]
all_pc = [ar["wall_by_band"]["N"][b]["per_crater_frac_poleward"] for b in bands]
all_pc_ci = np.array([ar["wall_by_band"]["N"][b]["per_crater_ci95"] for b in bands]).T
all_n = [ar["wall_by_band"]["N"][b]["n"] for b in bands]

re_ph = [rec["reassign_D20_100"]["N_bands"][b]["frac_poleward"] for b in bands]
re_pc = [rec["reassign_D20_100"]["N_bands"][b]["per_crater_frac_poleward"] for b in bands]
re_pc_ci = np.array([rec["reassign_D20_100"]["N_bands"][b]["per_crater_ci95"] for b in bands]).T
re_n = [rec["reassign_D20_100"]["N_bands"][b]["n"] for b in bands]

C_ALL = "#1f77b4"   # all hosts (D>=15, incl basins)
C_RE = "#d62728"    # D 20-100 km hosts

fig = plt.figure(figsize=(13, 4.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.35, 1.0], wspace=0.32)
axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1]); axC = fig.add_subplot(gs[2])

# ---------- Panel A ----------
axA.axhline(0.5, color="0.5", lw=0.8, ls="--", zorder=1)
# per-crater with CI (the robust weighting)
axA.errorbar(x-0.08, all_pc, yerr=[np.array(all_pc)-all_pc_ci[0], all_pc_ci[1]-np.array(all_pc)],
             fmt="o-", color=C_ALL, capsize=2, lw=1.6, ms=5, label="all hosts (D\u226515 km)", zorder=4)
axA.errorbar(x+0.08, re_pc, yerr=[np.array(re_pc)-re_pc_ci[0], re_pc_ci[1]-np.array(re_pc)],
             fmt="s-", color=C_RE, capsize=2, lw=1.6, ms=5, label="hosts D 20\u2013100 km", zorder=5)
# per-hollow as faint background points
axA.plot(x-0.08, all_ph, "o", color=C_ALL, alpha=0.35, ms=4, zorder=2)
axA.plot(x+0.08, re_ph, "s", color=C_RE, alpha=0.35, ms=4, zorder=2)
axA.set_xticks(x); axA.set_xticklabels([b+"\u00b0N" for b in bands])
axA.set_ylabel("frac. on equator-facing wall")
axA.set_ylim(0.20, 0.90)
axA.set_title("Poleward-wall preference vs latitude\n(per-crater \u00b195% CI; faint = per-hollow)")
axA.legend(frameon=False, fontsize=6, loc="upper left")
axA.annotate("reversal\n(all hosts)", xy=(1-0.08, all_ph[1]), xytext=(0.35, 0.30),
             fontsize=6, ha="center", color=C_ALL,
             arrowprops=dict(arrowstyle="->", color=C_ALL, lw=0.8))

# ---------- Panel B: 15-30N decomposition ----------
d = rec["decompose_15_30N"]
splits = [("all", d["overall"]),
          ("Cal<1200", d["by_caloris"]["lt1200km"]),
          ("Cal\u22651200", d["by_caloris"]["ge1200km"]),
          ("D 15-30", d["by_hostD"]["D15-30"]),
          ("D 30-100", d["by_hostD"]["D30-100"]),
          ("D >100", d["by_hostD"]["D>100"]),
          ("hot lon", d["by_hotwarm"]["hot"]),
          ("warm lon", d["by_hotwarm"]["warm"])]
xb = np.arange(len(splits))
ph = [s["frac_poleward"] for _, s in splits]
pc = [s["per_crater_frac_poleward"] for _, s in splits]
pc_ci = np.array([s["per_crater_ci95"] for _, s in splits]).T
nn = [s["n"] for _, s in splits]
axB.axhline(0.5, color="0.5", lw=0.8, ls="--", zorder=1)
axB.bar(xb-0.19, ph, width=0.36, color="#8c9fb0", label="per hollow", zorder=2)
axB.errorbar(xb+0.19, pc, yerr=[np.array(pc)-pc_ci[0], pc_ci[1]-np.array(pc)],
             fmt="D", color="#b5651d", capsize=2, ms=4, lw=1.2, label="per crater \u00b1CI", zorder=4)
for xi, n in zip(xb, nn):
    axB.text(xi-0.19, 0.02, f"n={n}", rotation=90, va="bottom", ha="center", fontsize=5, color="0.15")
axB.set_xticks(xb); axB.set_xticklabels([s[0] for s in splits], fontsize=6, rotation=30, ha="right")
axB.set_ylabel("frac. on equator-facing wall")
axB.set_ylim(0, 0.90)
axB.set_title("15\u201330\u00b0N reversal decomposed\n(per-hollow deficit vanishes per-crater)")
axB.legend(frameon=False, fontsize=6, loc="upper right")

# ---------- Panel C: DEM slope azimuth vs position facing ----------
diff = wall["dem_pos_diff_deg"].values
axC.hist(diff, bins=np.arange(0, 181, 15), color="#4c78a8", edgecolor="white")
axC.axvline(90, color="0.5", ls="--", lw=0.8)
p = rec["dem_proxy_check"]
axC.set_xlabel("|DEM slope azimuth \u2212 position facing| (deg)")
axC.set_ylabel("wall hollows")
axC.set_xlim(0, 180); axC.set_xticks([0, 45, 90, 135, 180])
axC.set_title("Aspect-proxy sanity check\n(8 ppd DEM vs position)")
axC.text(0.97, 0.95, f"agree <90\u00b0: {p['agreement_rate_lt90deg']*100:.0f}%\n"
                     f"agree <45\u00b0: {p['agreement_rate_lt45deg']*100:.0f}%\n"
                     f"median |\u0394|: {p['median_abs_diff_deg']:.0f}\u00b0\nn={p['n']}",
         transform=axC.transAxes, ha="right", va="top", fontsize=6.5,
         bbox=dict(boxstyle="round", fc="white", ec="0.7", lw=0.6))

for ax, L in [(axA, "a"), (axB, "b"), (axC, "c")]:
    panel_letter(ax, L)

fig.savefig(f"{FOLD}/aspect_robustness.png", dpi=300, bbox_inches="tight")
print("wrote aspect_robustness.png")

# §9.1 geometric check
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
tickset = {ax: set(ax.get_xticklabels() + ax.get_yticklabels()) for ax in fig.axes}
ov = [(a.get_text(), b.get_text()) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:]
      if ba.overlaps(bb)]
print("text overlaps:", ov[:10], "count", len(ov))
