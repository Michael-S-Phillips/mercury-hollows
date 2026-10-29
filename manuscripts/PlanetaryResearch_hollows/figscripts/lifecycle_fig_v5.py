#!/usr/bin/env python3
"""Figure 5 (manuscript v5): life cycle of a hollow in the zoned model (sulfidized breccia at depth, native
sulfur in cooler breccia above).
Panels 1-6 are schematic steps (source, charging, storage, activity, dormancy, fossilization); panel 7 is
data: daily maximum and minimum temperature with depth at the hot and warm equator from the thermal model
(thermal.py physics, depth step one-tenth of the diurnal skin depth), with the base of the layer that
crosses 368.5 K every Mercury solar day.
Usage: python3 lifecycle_fig_v4.py <dir with thermal.py> <out.png>"""
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, Patch, Circle
from matplotlib.lines import Line2D
sys.path.insert(0, sys.argv[1]); import thermal as th
OUT = sys.argv[2]
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5, "axes.linewidth": 0.6, "figure.facecolor": "white"})
INK, MUTED, BLUE, ORANGE = "#1a1a19", "#6b6b67", "#1F6FB2", "#D4541E"
SKY, BREC, BRECT, PLAINS = "#f4f2ed", "#d9c9b2", "#a8967d", "#b9c3cc"
SULF, SULFT, SULFIDE = "#E3B93A", "#9c7a12", "#4a3f33"
HOT, LAG, HALO, VAP, FRAC, LRM = "#e7a07c", "#3f3d3b", "#cfe3f1", "#2E8FC4", "#7a5a3a", "#b9a68c"
SULFZ, SULFZT = "#c3cba6", "#3F7A46"   # sulfidized breccia: Ca-Mg sulfide + silica (Renggli et al. 2022)

# ---------- thermal profiles (panel 7) ----------
def profiles(lat, lon, nz=140, nt=900, ndays=45):
    lat = np.atleast_1d(lat).astype(float); lon = np.atleast_1d(lon).astype(float)
    t, F = th._flux(lat, lon, nt); dt, dz = t[1]-t[0], th.DSKIN/10
    c = th.KAPPA*dt/dz**2; ns = lat.size; T = np.full((ns, nz), 300.0)
    for day in range(ndays):
        net = np.zeros(ns); Tx = np.full_like(T, -1.0); Tn = np.full_like(T, 1e4)
        for n in range(nt):
            T[:, 1:-1] += c*(T[:, 2:]-2*T[:, 1:-1]+T[:, :-2]); T[:, -1] = T[:, -2]
            Ts, T1 = T[:, 0].copy(), T[:, 1]
            for _ in range(8):
                f = th.EMIS*th.SIG*Ts**4 - F[:, n] + th.K*(Ts-T1)/dz
                Ts -= f/(4*th.EMIS*th.SIG*Ts**3 + th.K/dz)
            T[:, 0] = Ts; Tx = np.maximum(Tx, T); Tn = np.minimum(Tn, T)
            net += F[:, n] - th.EMIS*th.SIG*Ts**4
        net /= nt
        T += np.clip(net/(4*th.EMIS*th.SIG*np.maximum(T[:, -1], 50)**3), -40, 40)[:, None]
    return np.arange(nz)*dz*100, Tx, Tn
z, Tx, Tn = profiles([0, 0], [0, 90])
def crossing(tx, thr=368.5):
    j = np.where(tx > thr)[0].max(); return z[j] + (z[1]-z[0])*(tx[j]-thr)/(tx[j]-tx[j+1])
L = [crossing(Tx[k]) for k in range(2)]

# ---------- layout ----------
fig = plt.figure(figsize=(5.9, 4.15), dpi=100)
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.16], left=0.01, right=0.735, top=0.945, bottom=0.02, wspace=0.07, hspace=0.20)
W, ZT, ZB, ZBASE = 100.0, 14.0, -100.0, -46.0
rng0 = np.random.default_rng(3)
CLASTS = [(rng0.uniform(1, W-1), rng0.uniform(ZB+2, -2), rng0.uniform(0.6, 1.6)) for _ in range(150)]
FR = [(18, [(18, -100), (22, -80), (17, -62), (21, -46), (19, -30)]),
      (50, [(50, -100), (46, -84), (52, -66), (48, -50), (51, -34)]),
      (80, [(80, -100), (84, -82), (79, -64), (83, -48), (81, -32)])]

def block(ax, n, title, plains=False):
    ax.set_xlim(0, W); ax.set_ylim(ZB, ZT); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_color(MUTED); sp.set_linewidth(0.5)
    ax.add_patch(Rectangle((0, 0), W, ZT, fc=SKY, ec="none", zorder=0))
    ax.add_patch(Rectangle((0, ZB), W, -ZB, fc=BREC, ec="none", zorder=0))
    for x, y, r in CLASTS:
        ax.add_patch(Circle((x, y), r, fc=BRECT, ec="none", alpha=0.35, zorder=1))
    ax.text(0.0, 1.02, f"{n}", transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK, va="bottom")
    ax.text(0.09, 1.02, title, transform=ax.transAxes, fontsize=7, color=INK, va="bottom")

def bodies(ax, cut=None, alpha=1.0):
    rng = np.random.default_rng(11)
    for x0, _ in FR:
        for k in range(24):
            x = x0 + rng.normal(0, 9); y = rng.uniform(ZBASE+3, -6)
            if cut is not None and cut(x, y): continue
            rx, ry = rng.uniform(2.0, 4.6), rng.uniform(1.2, 2.8)
            th_ = np.linspace(0, 2*np.pi, 14); rr = 1 + 0.25*rng.normal(size=14)
            ax.add_patch(Polygon(np.c_[x + rx*rr*np.cos(th_), y + ry*rr*np.sin(th_)], closed=True,
                                 fc=SULF, ec=SULFT, lw=0.3, alpha=alpha, zorder=3))

def sulfides(ax, n=120):
    rng = np.random.default_rng(5)
    cen = rng.uniform(8, 92, 6), rng.uniform(-88, -42, 6)
    for i in range(n):
        k = i % 6; x = cen[0][k] + rng.normal(0, 7); y = cen[1][k] + rng.normal(0, 7)
        if 1 < x < W - 1 and ZB + 1 < y < -33: ax.add_patch(Circle((x, y), 0.9, fc=SULFIDE, ec="none", zorder=2))

def sulfidized(ax, top=-52.0, bot=-77.0):
    ax.add_patch(Rectangle((0, bot), W, top - bot, fc=SULFZ, ec="none", zorder=0.9))
    ax.add_patch(Rectangle((0, bot), W, top - bot, fc="none", ec=SULFZT, lw=0, hatch="....", alpha=0.35, zorder=0.95))

def fractures(ax, ls="-", lw=0.8):
    for _, pts in FR:
        p = np.array(pts); ax.plot(p[:, 0], p[:, 1], color=FRAC, lw=lw, ls=ls, zorder=2)

def base_line(ax, x0=0, x1=W):
    ax.plot([x0, x1], [ZBASE, ZBASE], color=INK, lw=0.7, ls=(0, (4, 2)), zorder=6)

def arrow(ax, x0, y0, x1, y1, color=VAP, lw=0.9, ms=5):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms, color=color, lw=lw, zorder=7))

def hollow(ax, xs, floor_lag=True, soft=False, lagscarp=False, halo=True):
    if soft:
        xx = np.linspace(0, W, 200); top = -5.0; bot = ZBASE + 10
        prof = bot + (top - bot) / (1 + np.exp(-(xx - xs) / 4.0))
        ax.add_patch(Polygon(np.r_[np.c_[xx, prof], [[W, ZT], [0, ZT]]], closed=True, fc=SKY, ec="none", zorder=5))
        ax.plot(xx, prof, color=INK, lw=0.6, zorder=6)
        ax.fill_between(xx, prof - 2.5, prof, color=LAG, alpha=0.55, lw=0, zorder=5.5)
        return
    pol = np.array([[0, ZT], [0, ZBASE], [xs - 6, ZBASE], [xs, -5], [xs, ZT]])
    ax.add_patch(Polygon(pol, closed=True, fc=SKY, ec="none", zorder=5))
    ax.plot([0, xs - 6, xs, W], [ZBASE, ZBASE, -5, -5], color=INK, lw=0.6, zorder=6)
    if floor_lag: ax.add_patch(Rectangle((0, ZBASE - 3), xs - 6, 3, fc=FRAC, ec="none", alpha=0.45, zorder=4))
    ax.add_patch(Rectangle((xs, -5), W - xs, 5, fc=LAG, ec="none", zorder=4))
    if lagscarp: ax.add_patch(Polygon([[xs - 6, ZBASE], [xs, -5], [xs + 0.01, -5], [xs - 3.5, ZBASE]], closed=True, fc=LAG, ec="none", zorder=6))
    if halo: ax.add_patch(Rectangle((xs, -1.2), 24, 1.6, fc=HALO, ec="#8fb6d4", lw=0.3, zorder=6))

# 1 source
ax = fig.add_subplot(gs[0, 0]); block(ax, "1", "Source")
ax.add_patch(Rectangle((0, -30), W, 30, fc=PLAINS, ec="none", zorder=1.5))
ax.add_patch(Rectangle((0, ZB), W, 70, fc=LRM, ec="none", alpha=0.55, zorder=1.4))
sulfides(ax)
ax.text(50, -15, "younger plains", ha="center", va="center", fontsize=6, color=INK, zorder=8)
ax.text(50, -88, "older crust with sulfides", ha="center", va="center", fontsize=6, color=INK, zorder=8,
        bbox=dict(fc=LRM, ec="none", pad=0.5, alpha=0.9))
# 2 charging
ax = fig.add_subplot(gs[0, 1]); block(ax, "2", "Charging"); sulfidized(ax)
ax.add_patch(Rectangle((0, ZB), W, 22, fc=HOT, ec="none", alpha=0.85, zorder=1.5))
ax.text(66, ZB + 11, "> 1200 K", ha="center", va="center", fontsize=6, color=INK, zorder=8)
fractures(ax); bodies(ax); base_line(ax)
for x0, pts in FR:
    p = np.array(pts)
    for i in (0, 2): arrow(ax, p[i, 0] + 1.5, p[i, 1] + 3, p[i + 1, 0] + 1.5, p[i + 1, 1] - 1)
ax.text(66, ZBASE - 2.0, "d", ha="center", va="top", fontsize=6.5, style="italic", color=INK, zorder=8,
        bbox=dict(fc=BREC, ec="none", pad=0.6))
# 3 storage
ax = fig.add_subplot(gs[0, 2]); block(ax, "3", "Storage"); sulfidized(ax)
fractures(ax, ls=(0, (2, 1.5)), lw=0.6); bodies(ax); base_line(ax)
ax.add_patch(Rectangle((0, -5), W, 5, fc=LAG, ec="none", zorder=4))
for x in (15, 40, 65, 90): arrow(ax, x, -4, x, 4, lw=0.6, ms=3.5)
# 4 activity
xs = 58.0
ax = fig.add_subplot(gs[1, 0]); block(ax, "4", "Activity"); sulfidized(ax)
fractures(ax, ls=(0, (2, 1.5)), lw=0.6); bodies(ax, cut=lambda x, y: x < xs - 2 + 0.25*(y - ZBASE))
hollow(ax, xs)
for (x, y) in [(xs - 3, -14), (xs - 5, -26), (xs - 6.5, -38)]: arrow(ax, x, y, x - 9, y + 6, lw=0.8, ms=4)
base_line(ax, xs + 1, W)
ax.text(xs + 2, 2.2, "halo", fontsize=6, color=INK, zorder=8)
# 5 dormancy
xs2 = 72.0
ax = fig.add_subplot(gs[1, 1]); block(ax, "5", "Dormancy"); sulfidized(ax)
fractures(ax, ls=(0, (2, 1.5)), lw=0.6); bodies(ax, cut=lambda x, y: x < xs2 - 2 + 0.25*(y - ZBASE))
hollow(ax, xs2, lagscarp=True, halo=False)
base_line(ax, xs2 + 1, W)
# 6 fossilization
ax = fig.add_subplot(gs[1, 2]); block(ax, "6", "Fossilization"); sulfidized(ax)
fractures(ax, ls=(0, (2, 1.5)), lw=0.6); bodies(ax, cut=lambda x, y: x < 72, alpha=0.8)
hollow(ax, 72, soft=True)
# 7 thermal profiles
ax = fig.add_axes([0.80, 0.18, 0.15, 0.765])
for k, col in [(0, ORANGE), (1, BLUE)]:
    ax.plot(Tx[k], z, color=col, lw=1.1, zorder=3); ax.plot(Tn[k], z, color=col, lw=1.1, ls=(0, (3, 1.5)), zorder=3)
    ax.plot([368.5, 760], [L[k], L[k]], color=col, lw=0.7, zorder=2)
    ax.text(745, L[k] + 0.6, f"{L[k]:.0f} cm", ha="right", va="top", fontsize=6, color=col)
ax.axvline(368.5, color=INK, lw=0.7, ls=(0, (1, 1.5)), zorder=2)
ax.text(382, 31.0, "368 K", fontsize=6, color=INK, va="bottom", rotation=90)
ax.set_xlim(0, 760); ax.set_ylim(36, 0)
ax.set_xlabel("Temperature (K)", labelpad=1); ax.set_ylabel("Depth (cm)", labelpad=1)
ax.yaxis.set_label_position("right"); ax.yaxis.tick_right(); ax.tick_params(length=2, pad=1, labelsize=6)
ax.text(0.0, 1.01, "7", transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK, va="bottom")
# legend
lg = fig.add_subplot(gs[2, 0:3]); lg.axis("off")
lg.legend(handles=[Patch(fc=SULFIDE, ec="none", label="sulfide"), Patch(fc=SULF, ec=SULFT, lw=0.3, label="native sulfur"),
                   Patch(fc=BREC, ec="none", label="breccia"), Patch(fc=SULFZ, ec=SULFZT, lw=0.3, hatch="....", label="sulfidized breccia"),
                   Patch(fc=PLAINS, ec="none", label="plains"),
                   Patch(fc=HOT, ec="none", label="hot zone"), Patch(fc=LAG, ec="none", label="sulfur-free lag"),
                   Patch(fc=HALO, ec="#8fb6d4", lw=0.3, label="halo"), Line2D([], [], color=FRAC, lw=0.8, label="fracture"),
                   Line2D([], [], color=VAP, lw=0.9, marker=">", ms=3, label="S$_2$ vapor"),
                   Line2D([], [], color=INK, lw=0.7, ls=(0, (4, 2)), label="base of charged layer")],
          loc="center", ncol=6, fontsize=5.4, frameon=False, handlelength=1.5, columnspacing=0.8, handletextpad=0.35, borderaxespad=0)
import os, figure_data; figure_data.dump(fig, os.path.splitext(os.path.basename(OUT))[0], os.path.join(os.path.dirname(OUT), "data"))  # archive: table of plotted values
fig.savefig(OUT, dpi=450, facecolor="white")
print("wrote", OUT, [round(v, 1) for v in L])
