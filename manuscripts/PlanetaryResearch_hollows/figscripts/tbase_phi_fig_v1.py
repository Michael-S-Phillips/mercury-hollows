"""Supplement figure: hollow depth for complete collapse of the charged layer, h = phi P_S(T)/(rho g), against the peak
temperature at the base of the layer (T_base) and the volume fraction lost on collapse (phi).
P_S(T): Clausius-Clapeyron relation anchored at the normal boiling point of sulfur (717.8 K, 0.101 MPa; Section S13), with the
mean of the two slopes of Section S13 (7,750 and 8,420 K) in ln P, i.e. B = 8,085 K. rho = 2,500 kg m-3, g = 3.70 m s-2.
Writes figures/figS_tbase_phi_v1.png and .pdf. The two slopes change T_base at fixed h by about +/-10-20 K over this range,
and rho = 2,000 or 3,000 kg m-3 by about -/+20 K (printed below)."""
import os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "figures")
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5,
    "axes.linewidth": .6, "axes.edgecolor": "#5a5a57", "xtick.color": "#5a5a57", "ytick.color": "#5a5a57",
    "xtick.major.width": .6, "ytick.major.width": .6, "xtick.minor.width": .4, "ytick.minor.width": .4,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.labelcolor": "#1a1a19", "text.color": "#1a1a19", "axes.labelsize": 7.5,
    "figure.facecolor": "white", "axes.facecolor": "white"})
INK, MUTED, ORANGE = "#1a1a19", "#6b6b67", "#D4541E"
Tb, Pb, g, rho = 717.8, 101325.0, 3.70, 2500.0
B = 0.5 * (7750.0 + 8420.0)
def P_S(T, B=B): return Pb * np.exp(B * (1 / Tb - 1 / T))
def T_of_P(P, B=B): return 1 / (1 / Tb - np.log(P / Pb) / B)
T = np.linspace(700, 1200, 501); phi = np.geomspace(0.05, 0.5, 400)
TT, FF = np.meshgrid(T, phi)
H = FF * P_S(TT) / (rho * g)
# sequential single hue (house blue), light -> dark
cmap = LinearSegmentedColormap.from_list("blue_seq", ["#f3f7fb", "#cfe0f1", "#8db8de", "#4a8cc6", "#1F6FB2", "#14508a", "#0b2f52"])
norm = LogNorm(vmin=0.4, vmax=500)
fig = plt.figure(figsize=(12 / 2.54, 8.6 / 2.54))
ax = fig.add_axes([0.105, 0.15, 0.70, 0.66])
cax = fig.add_axes([0.86, 0.15, 0.025, 0.66])
ax.pcolormesh(TT, FF, H, cmap=cmap, norm=norm, shading="auto", rasterized=True)
# observed band, 24-28 m
ax.contourf(TT, FF, H, levels=[24, 28], colors=[ORANGE], alpha=0.9, zorder=3)
levels = [1, 2, 5, 10, 50, 100, 200]          # 25 m is shown by the observed band
cs = ax.contour(TT, FF, H, levels=levels, colors=[INK if L <= 10 else "white" for L in levels], linewidths=0.55, zorder=4)
# label each contour at phi ~ 0.075 (right side) or at the top edge, chosen to avoid the band
def lab_pos(L):
    for ph in (0.065,):
        Tl = T_of_P(L / ph * rho * g)
        if 715 < Tl < 1185: return (Tl, ph)
    Tl = T_of_P(L / 0.42 * rho * g); return (Tl, 0.42)
ax.clabel(cs, levels, inline=True, inline_spacing=3, fontsize=6.3, fmt=lambda v: f"{v:g} m", manual=[lab_pos(L) for L in levels])
ax.set_yscale("log")
ax.set_xlim(700, 1200); ax.set_ylim(0.05, 0.5)
ax.set_yticks([0.05, 0.1, 0.2, 0.3, 0.4, 0.5]); ax.set_yticklabels(["0.05", "0.1", "0.2", "0.3", "0.4", "0.5"])
ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.set_xlabel("peak temperature at the base of the charged layer, $T_{\\mathrm{base}}$ (K)")
ax.set_ylabel("volume fraction lost on collapse, $\\phi$")
# top axis: thickness of the charged layer d = P_S(T)/(rho g)
top = ax.twiny(); top.set_xlim(ax.get_xlim())
dt = [10, 20, 50, 100, 200, 500, 1000]
top.set_xticks([T_of_P(d * rho * g) for d in dt]); top.set_xticklabels([f"{d:g}" for d in dt])
top.set_xlabel("thickness of the charged layer, $d$ (m)", labelpad=4)
for a in (ax, top):
    for sp in a.spines.values(): sp.set_color("#5a5a57")
# band label: white text parallel to the band, on its hot side
fig.canvas.draw()
ph0, ph1 = 0.13, 0.20
p0 = ax.transData.transform((T_of_P(26 / ph0 * rho * g), ph0)); p1 = ax.transData.transform((T_of_P(26 / ph1 * rho * g), ph1))
ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
ang = ang - 180 if ang > 90 else (ang + 180 if ang < -90 else ang)   # keep text upright
phm = 0.16; Tm = T_of_P(26 / phm * rho * g) + 27
ax.text(Tm, phm, "observed median depth, 24–28 m", rotation=ang, rotation_mode="anchor", ha="center", va="center",
        fontsize=6.5, color="white", zorder=6)
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.set_ticks([1, 2, 5, 10, 25, 50, 100, 200, 500]); cb.set_ticklabels(["1", "2", "5", "10", "25", "50", "100", "200", "500"])
cb.ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor("#5a5a57")
cb.set_label("hollow depth, $h = \\phi d$ (m)")
cb.ax.axhspan(24, 28, color=ORANGE, lw=0)
for ext in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"figS_tbase_phi_v1.{ext}"), dpi=400)
# numbers for the caption
for phi_ in (0.1, 0.2, 0.4):
    Tm = [T_of_P(h / phi_ * rho * g) for h in (24, 28)]
    Tsl = [T_of_P(h / phi_ * rho * g, b) for h in (24, 28) for b in (7750, 8420)]
    Tr = [T_of_P(h / phi_ * r * g) for h in (24, 28) for r in (2000, 3000)]
    print(f"phi={phi_}: T_base {Tm[0]:.0f}-{Tm[1]:.0f} K (rho 2500, B 8085); slopes {min(Tsl):.0f}-{max(Tsl):.0f}; rho 2000-3000 {min(Tr):.0f}-{max(Tr):.0f}")
print("d ticks at T:", {d: round(float(T_of_P(d * rho * g))) for d in dt})
print("h range in panel:", H.min(), H.max())
