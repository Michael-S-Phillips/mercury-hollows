#!/usr/bin/env python3
"""Figure 7 (figure v11, 9 Nov 2026; Michael's v6.10 comment 83): the thermal panel b becomes inset i of step 4. A dashed box
labeled i on the scarp of step 4 marks where the inset applies; the inset keeps its place in the grid, gets the dashed frame of the
other inset, loses its panel letter, and shades the zone above the 368 K lines as the zone of activity. Otherwise as v10.
Figure 5 (figure v10): lifecycle_fig_v9.py with the impact added. An unnumbered context pane, "Impact (context)", sits
between steps 1 and 2: a half cross-section of a complex crater cut through the younger plains into the older crust,
with ejecta of older crust on the plains, a central uplift that raises deeper crust, impact melt on the floor, and
impact-heated rock beneath the floor; a box marks the part of the floor that steps 2-6 enlarge, and an arrow leads from
it to step 2. A small inset beneath it, framed with a dashed line, shows a possible second source of heat or sulfur for
isolated host craters: magma intruded as a sill beneath a crater floor, which heats the floor, supplies sulfur, and feeds an explosive
vent. The panes run in a 4 x 2 grid (1, impact, 2, 3 / 4, 5, 6, b); the step panes are 1.30 in wide instead of 1.36 in.
Figure 5 (manuscript v4, figure v9): identical to lifecycle_fig_v8.py except that step 2 is titled "Charging and collapse".
Figure 5 (manuscript v8): life cycle of a hollow in the zoned model (sulfidized regolith at depth, native
sulfur in cooler regolith above, compacted sulfur-free residue at the base of the charged layer).
Panel a holds the six schematic steps, numbered 1-6 (source, charging and collapse, storage, activity, dormancy, fossilization).
Panel b is data: daily maximum and minimum temperature with depth at the hot and warm equator from the thermal model
(thermal.py physics, depth step one-tenth of the diurnal skin depth), with the base of the layer that crosses 368.5 K
every Mercury solar day.
Changes from v5: panels lettered a (schematic) and b (thermal); "breccia" relabeled "regolith"; hot zone > 1,300 K;
new unit "compacted residue" (regolith that has lost its sulfur and compacted); shallow surface sag over the charged
layer from step 2 on; residue band above the base d in steps 2-3; in steps 4-6 the hollow floor lies at depth h,
above the base d, with residue beneath it, and d runs across the whole panel.
Usage: python3 lifecycle_fig_v10.py <dir with thermal.py> <out.png>   (also writes <out>.pdf)"""
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch, Patch, Circle, ConnectionPatch
from matplotlib.lines import Line2D
from matplotlib.path import Path
sys.path.insert(0, sys.argv[1]); import thermal as th
OUT = sys.argv[2]; OUT_PDF = os.path.splitext(OUT)[0] + ".pdf"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5, "axes.linewidth": 0.6, "figure.facecolor": "white",
                     "pdf.fonttype": 42})
INK, MUTED, BLUE, ORANGE = "#1a1a19", "#6b6b67", "#1F6FB2", "#D4541E"
SKY, REG, REGT, PLAINS = "#f4f2ed", "#d9c9b2", "#a8967d", "#b9c3cc"
SULF, SULFT, SULFIDE = "#E3B93A", "#9c7a12", "#4a3f33"
HOT, LAG, HALO, VAP, FRAC, LRM = "#e7a07c", "#3f3d3b", "#cfe3f1", "#2E8FC4", "#7a5a3a", "#b9a68c"
SULFZ, SULFZT = "#c3cba6", "#3F7A46"   # sulfidized regolith: Ca-Mg sulfide + silica (Renggli et al. 2022)
RESID, RESIDT = "#c2baaf", "#6f685f"   # compacted residue: regolith that has lost its sulfur and compacted
HATCH_LW = 0.35

# ---------- thermal profiles (panel b) ----------
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

# ---------- layout (inches; 5.9 in = 15 cm) ----------
FW, FH = 5.9, 4.30
fig = plt.figure(figsize=(FW, FH), dpi=100)
def axes_in(x, y, w, h): return fig.add_axes([x/FW, y/FH, w/FW, h/FH])
PX0, PW, PGX, PH = 0.06, 1.30, 0.10, 1.55          # schematic panel left edge, width, gap, height (4 columns)
ROWY = (2.37, 0.58)                                  # bottom of schematic rows 1 and 2
def panel(row, col): return axes_in(PX0 + col*(PW + PGX), ROWY[row], PW, PH)

# ---------- schematic geometry (data units; surface before charging at 0) ----------
W, ZT, ZB, ZBASE = 100.0, 14.0, -100.0, -46.0      # ZBASE: base of the charged layer (depth d)
SAG2, SAG = 3.5, 5.0                                 # surface sag over the charged layer: step 2, steps 3-6
LT = 3.5                                             # lag thickness
RES2, RES3 = ZBASE + 5.0, -30.0                      # top of the compacted residue in steps 2 and 3
YF = -25.0                                           # hollow floor in steps 4-6 (h = SAG + YF below the surroundings)
UNERODED_RESIDUE_TOP = None   # None: native sulfur down to d beside the hollow (steps 4-6); a y value adds residue below it
rng0 = np.random.default_rng(3)
CLASTS = [(rng0.uniform(1, W-1), rng0.uniform(ZB+2, -2), rng0.uniform(0.6, 1.6)) for _ in range(150)]
FR = [(18, [(18, -100), (22, -80), (17, -62), (21, -46), (19, -30)]),
      (50, [(50, -100), (46, -84), (52, -66), (48, -50), (51, -34)]),
      (80, [(80, -100), (84, -82), (79, -64), (83, -48), (81, -32)])]
_rb = np.random.default_rng(11)
BODY = [(x0 + _rb.normal(0, 9), (k + _rb.uniform(0, 1))/24, _rb.uniform(2.0, 4.6), _rb.uniform(1.2, 2.8), 1 + 0.25*_rb.normal(size=14))
        for x0, _ in FR for k in range(24)]   # depth fraction stratified so bodies fill the layer evenly

def sagw(x, a=3.0, b=20.0):
    """Sag weight: 0 at the panel edges, 1 over the charged body, cosine ramps between."""
    x = np.asarray(x, float)
    l = np.clip((x - a)/(b - a), 0, 1); r = np.clip(((W - a) - x)/(b - a), 0, 1)
    return 0.25*(1 - np.cos(np.pi*l))*(1 - np.cos(np.pi*r))
def surf(x, sag=SAG): return -sag*sagw(x)
XX = np.linspace(0, W, 301)

def block(ax, n, title):
    ax.set_xlim(0, W); ax.set_ylim(ZB, ZT); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_color(MUTED); sp.set_linewidth(0.5)
    ax.add_patch(Rectangle((0, 0), W, ZT, fc=SKY, ec="none", zorder=0))
    ax.add_patch(Rectangle((0, ZB), W, -ZB, fc=REG, ec="none", zorder=0))
    for x, y, r in CLASTS:
        ax.add_patch(Circle((x, y), r, fc=REGT, ec="none", alpha=0.35, zorder=1))
    ax.text(0.0, 1.02, f"{n}", transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK, va="bottom")
    ax.text(0.09, 1.02, title, transform=ax.transAxes, fontsize=7, color=INK, va="bottom")

def bodies(ax, ylo, yhi, region=None, alpha=1.0, ys=1.0):
    """Native-sulfur bodies between ylo and yhi; with region (data-space polygon), only bodies centered inside it,
    clipped to it."""
    clip = None
    if region is not None:
        region = np.asarray(region, float); inside = Path(region).contains_point
        clip = Polygon(region, closed=True, transform=ax.transData)
    th_ = np.linspace(0, 2*np.pi, 14)
    for x, u, rx, ry, rr in BODY:
        y = ylo + u*(yhi - ylo)
        if clip is not None and not inside((x, y)): continue
        p = Polygon(np.c_[x + rx*rr*np.cos(th_), y + ys*ry*rr*np.sin(th_)], closed=True,
                    fc=SULF, ec=SULFT, lw=0.3, alpha=alpha, zorder=3)
        ax.add_patch(p)
        if clip is not None: p.set_clip_path(clip)

def sulfides(ax, n=120):
    rng = np.random.default_rng(5)
    cen = rng.uniform(8, 92, 6), rng.uniform(-88, -42, 6)
    for i in range(n):
        k = i % 6; x = cen[0][k] + rng.normal(0, 7); y = cen[1][k] + rng.normal(0, 7)
        if 1 < x < W - 1 and ZB + 1 < y < -33: ax.add_patch(Circle((x, y), 0.9, fc=SULFIDE, ec="none", zorder=2))

def sulfidized(ax, top=-52.0, bot=-77.0):
    ax.add_patch(Rectangle((0, bot), W, top - bot, fc=SULFZ, ec="none", zorder=0.9))
    ax.add_patch(Rectangle((0, bot), W, top - bot, fc="none", ec=SULFZT, lw=0, hatch="....", alpha=0.35, zorder=0.95))

def residue(ax, pts):
    """Compacted residue: darker, desaturated regolith with a fine diagonal hatch."""
    ax.add_patch(Polygon(pts, closed=True, fc=RESID, ec="none", zorder=0.92))
    hp = Polygon(pts, closed=True, fc="none", ec=RESIDT, lw=0, hatch="//////", alpha=0.55, zorder=0.96)
    hp.set_hatch_linewidth(HATCH_LW); ax.add_patch(hp)

def ground(ax, sag, lag=False):
    """Sky above the sagged surface, and the sulfur-free lag that follows it."""
    s = surf(XX, sag)
    ax.add_patch(Polygon(np.r_[np.c_[XX, s], [[W, ZT], [0, ZT]]], closed=True, fc=SKY, ec="none", zorder=4.5))
    if lag: ax.add_patch(Polygon(np.r_[np.c_[XX, s], np.c_[XX[::-1], s[::-1] - LT]], closed=True, fc=LAG, ec="none", zorder=4))
    else: ax.plot(XX, s, color=INK, lw=0.5, zorder=6)

def fractures(ax, ls="-", lw=0.8):
    for _, pts in FR:
        p = np.array(pts); ax.plot(p[:, 0], p[:, 1], color=FRAC, lw=lw, ls=ls, zorder=2)

def base_line(ax, x0=0, x1=W):
    ax.plot([x0, x1], [ZBASE, ZBASE], color=INK, lw=0.7, ls=(0, (4, 2)), zorder=6)

def dlabel(ax, x=66):
    ax.text(x, ZBASE - 2.0, "d", ha="center", va="top", fontsize=6.5, style="italic", color=INK, zorder=8,
            bbox=dict(fc=REG, ec="none", pad=0.6))

def arrow(ax, x0, y0, x1, y1, color=VAP, lw=0.9, ms=5):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms, color=color, lw=lw, zorder=7))

def uneroded(ax, region, alpha=1.0):
    """Charged layer beside a hollow: native sulfur from beneath the lag down to d (optionally residue at the base)."""
    ylo = ZBASE + 3
    if UNERODED_RESIDUE_TOP is not None:
        clip = Polygon(region, closed=True, transform=ax.transData)
        pts = [[0, ZBASE], [W, ZBASE], [W, UNERODED_RESIDUE_TOP], [0, UNERODED_RESIDUE_TOP]]
        f = Polygon(pts, closed=True, fc=RESID, ec="none", zorder=0.92); ax.add_patch(f); f.set_clip_path(clip)
        hp = Polygon(pts, closed=True, fc="none", ec=RESIDT, lw=0, hatch="//////", alpha=0.55, zorder=0.96)
        hp.set_hatch_linewidth(HATCH_LW); ax.add_patch(hp); hp.set_clip_path(clip)
        ylo = UNERODED_RESIDUE_TOP + 3
    bodies(ax, ylo, -SAG - LT - 4, region=region, alpha=alpha)

def scarp_hollow(ax, xs, active):
    """Steps 4-5: hollow floor at YF (depth h below the surroundings), steep scarp at xs, residue beneath the floor
    down to the base d, which does not rise."""
    S0 = -SAG; top = S0 - LT; xt = xs - 4.0; xb = xt - 3.0
    xr = XX[XX >= xs]; sr = surf(xr)
    residue(ax, [[0, YF], [xt, YF], [xb, ZBASE], [0, ZBASE]])
    uneroded(ax, np.r_[[[xb, ZBASE], [xt, YF], [xs, top]], np.c_[xr, sr - LT], [[W, ZBASE]]])
    prof = np.r_[[[0, YF], [xt, YF], [xs, top], [xs, S0]], np.c_[xr, sr]]
    ax.add_patch(Polygon(np.r_[prof, [[W, ZT], [0, ZT]]], closed=True, fc=SKY, ec="none", zorder=5))
    ax.add_patch(Polygon(np.r_[np.c_[xr, sr], np.c_[xr[::-1], sr[::-1] - LT]], closed=True, fc=LAG, ec="none", zorder=4))
    if not active:   # lag mantles the floor and the scarp
        ax.add_patch(Polygon([[0, YF], [xt, YF], [xs, top], [xs + 0.7, top], [xt + 2.8, YF - 2.2], [0, YF - 2.2]],
                             closed=True, fc=LAG, ec="none", zorder=5.8))
    ax.plot([0, xt, xs], [YF, YF, top], color=INK, lw=0.6, zorder=6.2)
    base_line(ax)
    return S0, top, xt

# ---------- panel a: steps 1-6 ----------
# 1 source
ax = panel(0, 0); block(ax, "1", "Source")
ax.add_patch(Rectangle((0, -30), W, 30, fc=PLAINS, ec="none", zorder=1.5))
ax.add_patch(Rectangle((0, ZB), W, 70, fc=LRM, ec="none", alpha=0.55, zorder=1.4))
sulfides(ax)
ax.text(50, -15, "younger plains", ha="center", va="center", fontsize=6, color=INK, zorder=8)
ax.text(50, -88, "older crust with sulfides", ha="center", va="center", fontsize=6, color=INK, zorder=8,
        bbox=dict(fc=LRM, ec="none", pad=0.5, alpha=0.9))
# impact (context; not a numbered step): half cross-section of a complex crater at crater scale, crater axis at x = 0.
# The impact excavates through the younger plains into the older crust and spreads ejecta of older crust over the plains;
# the central uplift raises deeper crust into the central peak; impact melt lines the floor, and impact-heated rock lies
# beneath it. The box marks the part of the floor that steps 2-6 enlarge.
MELT, DEEP, CRUST = "#b9532f", "#9b8a72", "#c8b79e"   # impact melt and magma; deeper crust; older crust (REG under LRM)
IX0, IX1, IYB, IYT = -1.55, 0.52, -1.0, 0.33
RRIM, RFE, RPK, YFL, YPK, YRIM, PLB = 1.0, 0.60, 0.20, -0.30, -0.04, 0.11, -0.16
HZW, HZD = 0.72, 0.62                                   # heated rock: half-width and depth below the floor
BOX = (-0.56, -0.58, 0.28, 0.35)                        # x, y, width, height of the box enlarged in steps 2-6
XI = np.linspace(IX0, IX1, 801)
def c_up(x): return 0.065*np.exp(-np.abs(np.abs(x) - RRIM)/0.25)    # structural uplift of the plains at the rim
def c_surf(x):
    r = np.abs(np.asarray(x, float)); y = np.full_like(r, YFL)
    pk = r < RPK; y[pk] = YFL + (YPK - YFL)*0.5*(1 + np.cos(np.pi*r[pk]/RPK))
    wl = (r >= RFE) & (r < RRIM); t = (r[wl] - RFE)/(RRIM - RFE)
    y[wl] = YFL + (YRIM - YFL)*(t - 0.8*np.sin(6*np.pi*t)/(6*np.pi))    # terraced wall
    out = r >= RRIM; y[out] = c_up(r[out]) + 0.045*np.exp(-(r[out] - RRIM)/0.35)   # rim uplift plus ejecta
    return y
def c_deep(x): return -0.70 + 0.62*np.exp(-(np.asarray(x, float)/0.30)**2)   # top of the deeper crust
def c_heat(x): return YFL - HZD*np.sqrt(np.clip(1 - (np.asarray(x, float)/HZW)**2, 0, 1))
SXI = c_surf(XI)
IMP_BOTTOM, IMP_INSET_H, IMP_GAP = ROWY[0], 0.47, 0.14
axI = axes_in(PX0 + (PW + PGX), IMP_BOTTOM + IMP_INSET_H + IMP_GAP, PW, PH - IMP_INSET_H - IMP_GAP)
ax = axI
ax.set_xlim(IX0, IX1); ax.set_ylim(IYB, IYT); ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values(): sp.set_color(MUTED); sp.set_linewidth(0.5)
ax.text(0.0, 1.02 + 0.0, "Impact (context)", transform=ax.transAxes, fontsize=7, color=INK, va="bottom")
ax.add_patch(Rectangle((IX0, IYB), IX1 - IX0, IYT - IYB, fc=SKY, ec="none", zorder=0))
gpoly = np.r_[np.c_[XI, SXI], [[IX1, IYB], [IX0, IYB]]]
ax.add_patch(Polygon(gpoly, closed=True, fc=CRUST, ec="none", zorder=0.5))
def clipped(patch):
    ax.add_patch(patch); patch.set_clip_path(Polygon(gpoly, closed=True, transform=ax.transData)); return patch
clipped(Polygon(np.r_[np.c_[XI, c_deep(XI)], [[IX1, IYB], [IX0, IYB]]], closed=True, fc=DEEP, ec="none", zorder=0.6))
rng = np.random.default_rng(7)
for _ in range(90):
    x, y = rng.uniform(IX0 + 0.03, IX1 - 0.03), rng.uniform(IYB + 0.03, -0.22)
    if y > np.interp(x, XI, SXI) - 0.04 or y > PLB + c_up(x) - 0.03: continue
    if -1.55 < x < -0.70 and -0.56 < y < -0.36: continue        # keep the "older crust" label clear
    if -1.55 < x < -0.60 and y < -0.76: continue                 # keep the "deeper crust" label clear
    clipped(Circle((x, y), 0.011, fc=SULFIDE, ec="none", zorder=0.7))
pt, pb = c_up(XI), PLB + c_up(XI)
clipped(Polygon(np.r_[np.c_[XI, pt], np.c_[XI[::-1], pb[::-1]]], closed=True, fc=PLAINS, ec="none", zorder=1.0))
m = np.abs(XI) >= RRIM
for sgn in (-1, 1):
    mm = m & (np.sign(XI) == sgn)
    if mm.sum() > 1:
        xe = XI[mm]; ax.add_patch(Polygon(np.r_[np.c_[xe, c_surf(xe)], np.c_[xe[::-1], c_up(xe)[::-1]]], closed=True,
                                          fc=CRUST, ec="none", zorder=1.1))
xh = XI[np.abs(XI) < HZW]
clipped(Polygon(np.r_[np.c_[xh, c_heat(xh)], [[xh[-1], IYT], [xh[0], IYT]]], closed=True, fc=HOT, ec="none", alpha=0.8, zorder=1.5))
ln, = ax.plot(XI, c_deep(XI), color=MUTED, lw=0.4, ls=(0, (2, 1.5)), zorder=1.6)
ln.set_clip_path(Polygon(gpoly, closed=True, transform=ax.transData))
xm = XI[(np.abs(XI) > 0.6*RPK) & (np.abs(XI) < RFE + 0.04)]
for sgn in (-1, 1):
    xs_ = xm[np.sign(xm) == sgn]
    if xs_.size > 1:
        clipped(Polygon(np.r_[np.c_[xs_, c_surf(xs_) + 0.01], np.c_[xs_[::-1], c_surf(xs_[::-1]) - 0.04]], closed=True,
                        fc=MELT, ec="none", zorder=1.7))
for pts in ([(-0.47, -0.86), (-0.43, -0.70), (-0.48, -0.55), (-0.44, -0.42), (-0.46, -0.32)],
            [(-0.66, -0.62), (-0.62, -0.50), (-0.66, -0.38), (-0.63, -0.30)]):
    p = np.array(pts); l2, = ax.plot(p[:, 0], p[:, 1], color=FRAC, lw=0.5, zorder=2)
    l2.set_clip_path(Polygon(gpoly, closed=True, transform=ax.transData))
ax.plot(XI, SXI, color=INK, lw=0.5, zorder=3)
ax.add_patch(Rectangle(BOX[:2], BOX[2], BOX[3], fc="none", ec=INK, lw=0.7, zorder=6))
ax.add_patch(FancyArrowPatch((0.0, -0.56), (0.0, -0.13), arrowstyle="-|>", mutation_scale=4.5, color=INK, lw=0.7, zorder=7))
ax.text(0.05, -0.38, "uplift", ha="left", va="center", fontsize=6, color=INK, zorder=8)
ax.text(-1.30, -0.065, "plains", ha="center", va="center", fontsize=6, color=INK, zorder=8)
ax.text(-1.52, -0.46, "older crust", ha="left", va="center", fontsize=6, color=INK, zorder=8)
ax.text(-1.52, -0.87, "deeper crust", ha="left", va="center", fontsize=6, color=INK, zorder=8)
ax.text(-0.04, -0.71, "heated rock", ha="center", va="center", fontsize=6, color=INK, zorder=8)
ax.annotate("melt", xy=(0.36, -0.30), xytext=(0.33, -0.13), fontsize=6, color=INK, ha="center", va="center", zorder=8,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.4, shrinkA=0.5, shrinkB=0))

# possible second source of heat or sulfur (inset): magma intruded beneath a crater floor as a sill heats the floor and supplies sulfur, and
# feeds an explosive vent on the floor
axV = axes_in(PX0 + (PW + PGX), IMP_BOTTOM, PW, IMP_INSET_H)
ax = axV
JX0, JX1, JYB, JYT = -1.0, 1.0, -1.0, 0.42
ax.set_xlim(JX0, JX1); ax.set_ylim(JYB, JYT); ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values(): sp.set_color(MUTED); sp.set_linewidth(0.5); sp.set_linestyle((0, (3, 2)))
ax.text(0.0, 1.03, "Second source: intrusion", transform=ax.transAxes, fontsize=6, style="italic", color=INK, va="bottom")
XV = np.linspace(JX0, JX1, 801)
def v_surf(x):
    x = np.asarray(x, float); r = np.abs(x); y = np.where(r < 0.42, -0.22, 0.0)
    wl = (r >= 0.42) & (r < 0.78); t = (r[wl] - 0.42)/0.36; y[wl] = -0.22 + 0.32*0.5*(1 - np.cos(np.pi*t))
    out = r >= 0.78; y[out] = 0.10*np.exp(-(r[out] - 0.78)/0.12)
    return y - 0.09*np.exp(-((x - 0.22)/0.045)**2)       # vent pit on the floor
SV = v_surf(XV)
vpoly = np.r_[np.c_[XV, SV], [[JX1, JYB], [JX0, JYB]]]
ax.add_patch(Rectangle((JX0, JYB), JX1 - JX0, JYT - JYB, fc=SKY, ec="none", zorder=0))
ax.add_patch(Polygon(vpoly, closed=True, fc=CRUST, ec="none", zorder=0.5))
def vclip(patch):
    ax.add_patch(patch); patch.set_clip_path(Polygon(vpoly, closed=True, transform=ax.transData)); return patch
vu = 0.04*np.exp(-np.abs(np.abs(XV) - 0.78)/0.15)
vclip(Polygon(np.r_[np.c_[XV, vu], np.c_[XV[::-1], (vu - 0.13)[::-1]]], closed=True, fc=PLAINS, ec="none", zorder=1.0))
th_ = np.linspace(0, 2*np.pi, 80)
vclip(Polygon(np.c_[0.66*np.cos(th_), -0.60 + 0.26*np.sin(th_)], closed=True, fc=HOT, ec="none", alpha=0.8, zorder=1.5))
xs_ = np.linspace(-0.48, 0.48, 60); hw = 0.10*np.sqrt(1 - (xs_/0.48)**2)
ax.add_patch(Polygon(np.r_[np.c_[xs_, -0.60 + hw], np.c_[xs_[::-1], (-0.60 - hw)[::-1]]], closed=True, fc=MELT, ec="none", zorder=2))
ax.plot([-0.30, -0.22], [JYB, -0.66], color=MELT, lw=1.4, solid_capstyle="butt", zorder=2)
ax.plot([0.24, 0.22], [-0.56, -0.29], color=MELT, lw=0.9, solid_capstyle="butt", zorder=2)
for x in (-0.28, -0.04):
    arrow(ax, x, -0.47, x, -0.27, lw=0.7, ms=4)
xd = np.linspace(0.02, 0.42, 60); td = 0.07*np.exp(-((xd - 0.22)/0.09)**2)   # pyroclastic deposit around the vent
ax.add_patch(Polygon(np.r_[np.c_[xd, v_surf(xd) + td], np.c_[xd[::-1], v_surf(xd[::-1])]], closed=True, fc=MELT, ec="none",
                     alpha=0.75, zorder=2.5))
ax.plot(XV, SV, color=INK, lw=0.5, zorder=3)
ax.text(-0.94, -0.60, "sill", ha="left", va="center", fontsize=6, color=INK, zorder=8)
ax.text(0.22, -0.03, "vent", ha="center", va="center", fontsize=6, color=INK, zorder=8)
# 2 charging: hot zone at depth, S2 rising along fractures, sulfidized regolith near the hot zone, native sulfur in the
# cooler regolith above d; the lower part of the charged layer already loses its sulfur and compacts, and the surface sags
ax = panel(0, 2); block(ax, "2", "Charging and collapse"); sulfidized(ax)
ax.add_patch(Rectangle((0, ZB), W, 22, fc=HOT, ec="none", alpha=0.85, zorder=1.5))
ax.text(66.5, ZB + 17.5, "> 1,300 K", ha="center", va="center", fontsize=6, color=INK, zorder=8)
residue(ax, [[0, ZBASE], [W, ZBASE], [W, RES2], [0, RES2]])
fractures(ax); bodies(ax, RES2 + 3, -SAG2 - 5.5); ground(ax, SAG2); base_line(ax)
for x0, pts in FR:
    p = np.array(pts)
    for i in (0, 2): arrow(ax, p[i, 0] + 1.5, p[i, 1] + 3, p[i + 1, 0] + 1.5, p[i + 1, 1] - 1)
dlabel(ax)
# arrow from the box in the impact pane to step 2, which enlarges it
fig.add_artist(ConnectionPatch(xyA=(BOX[0] + BOX[2], BOX[1] + BOX[3]), coordsA=axI.transData, xyB=(0.0, 0.94),
                               coordsB=ax.transAxes, arrowstyle="-|>", mutation_scale=6, color=INK, lw=0.7,
                               connectionstyle="arc3,rad=-0.4", shrinkA=1, shrinkB=1, zorder=20))
# 3 storage: lag on top following the sag; native sulfur in the upper part, thicker residue below it down to d
ax = panel(0, 3); block(ax, "3", "Storage"); sulfidized(ax)
residue(ax, [[0, ZBASE], [W, ZBASE], [W, RES3], [0, RES3]])
fractures(ax, ls=(0, (2, 1.5)), lw=0.6); bodies(ax, RES3 + 3, -SAG - LT - 4, ys=0.85); ground(ax, SAG, lag=True); base_line(ax)
for x in (15, 40, 65, 90): arrow(ax, x, surf(x) - 1, x, surf(x) + 7, lw=0.6, ms=3.5)
# 4 activity
xs = 58.0
ax = panel(1, 0); block(ax, "4", "Activity"); sulfidized(ax)
fractures(ax, ls=(0, (2, 1.5)), lw=0.6)
S0, top, xt = scarp_hollow(ax, xs, active=True)
for y in (-12.0, -17.0, -22.0):
    xf = xt + (xs - xt)*(y - YF)/(top - YF)
    arrow(ax, xf - 1.0, y, xf - 10.0, y + 4.0, lw=0.8, ms=4)
ax.add_patch(Rectangle((xs, S0 - 1.2), 24, 1.6, fc=HALO, ec="#8fb6d4", lw=0.3, zorder=6))
ax.add_patch(Rectangle((xs - 7.0, S0 - 5.5), 9.0, 7.0, fc="none", ec=INK, lw=0.6, ls=(0, (2, 1.2)), zorder=9))   # where inset i applies
ax.text(xs - 8.0, S0 + 2.2, "i", fontsize=6.5, fontweight="bold", color=INK, ha="right", va="bottom", zorder=9)
ax.text(xs + 2, S0 + 2.2, "halo", fontsize=6, color=INK, zorder=8)
xh = 28.0                                                  # floor depth h below the surrounding surface
ax.plot([xh - 3, xs], [S0, S0], color=MUTED, lw=0.5, ls=(0, (1, 1.2)), zorder=6)
ax.add_patch(FancyArrowPatch((xh, S0), (xh, YF), arrowstyle="<|-|>", mutation_scale=3.5, color=INK, lw=0.6,
                             shrinkA=0, shrinkB=0, zorder=7))
ax.text(xh - 2.0, 0.5*(S0 + YF), "h", ha="right", va="center", fontsize=6.5, style="italic", color=INK, zorder=8)
dlabel(ax)
# 5 dormancy
xs2 = 72.0
ax = panel(1, 1); block(ax, "5", "Dormancy"); sulfidized(ax)
fractures(ax, ls=(0, (2, 1.5)), lw=0.6)
scarp_hollow(ax, xs2, active=False)
# 6 fossilization: softened scarp; residue beneath the floor and the base d retained
xc = 72.0
ax = panel(1, 2); block(ax, "6", "Fossilization"); sulfidized(ax)
fractures(ax, ls=(0, (2, 1.5)), lw=0.6)
prof = YF + (surf(XX) - YF)/(1 + np.exp(-(XX - xc)/4.0))
xbt, xbb = 66.0, 63.0; ybt = np.interp(xbt, XX, prof); m = XX < xbt
residue(ax, np.r_[[[0, ZBASE]], np.c_[XX[m], prof[m]], [[xbt, ybt], [xbb, ZBASE]]])
uneroded(ax, np.r_[[[xbb, ZBASE], [xbt, ybt - 2.5]], np.c_[XX[~m], prof[~m] - 2.5], [[W, ZBASE]]], alpha=0.8)
ax.add_patch(Polygon(np.r_[np.c_[XX, prof], [[W, ZT], [0, ZT]]], closed=True, fc=SKY, ec="none", zorder=5))
ax.plot(XX, prof, color=INK, lw=0.6, zorder=6)
ax.fill_between(XX, prof - 2.5, prof, color=LAG, alpha=0.55, lw=0, zorder=5.5)
base_line(ax)

# ---------- inset i (was panel b): thermal profiles at the exposed scarp of step 4 ----------
BX = PX0 + 3*(PW + PGX) + 0.12
ax = axes_in(BX, ROWY[1], PW - 0.22, PH)
for s_ in ax.spines.values(): s_.set_linestyle((0, (2, 1.2))); s_.set_linewidth(0.6)
ax.axhspan(0, L[0], color=ORANGE, alpha=0.10, lw=0, zorder=1); ax.axhspan(0, L[1], color=BLUE, alpha=0.10, lw=0, zorder=1)
ax.text(290, 4.8, "zone of\nactivity", fontsize=6, color=INK, ha="center", va="center", linespacing=1.0, zorder=6, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
ax.text(0.03, 1.02, "i", transform=ax.transAxes, fontsize=6.5, fontweight="bold", color=INK, ha="left", va="bottom")
ax.text(0.12, 1.02, "step 4, top of the scarp", transform=ax.transAxes, fontsize=6, color=INK, ha="left", va="bottom")
for k, col in [(0, ORANGE), (1, BLUE)]:
    ax.plot(Tx[k], z, color=col, lw=1.1, zorder=3); ax.plot(Tn[k], z, color=col, lw=1.1, ls=(0, (3, 1.5)), zorder=3)
    ax.plot([368.5, 760], [L[k], L[k]], color=col, lw=0.7, zorder=2)
    ax.text(745, L[k] + 0.6, f"{L[k]:.0f} cm", ha="right", va="top", fontsize=6, color=col)
ax.axvline(368.5, color=INK, lw=0.7, ls=(0, (1, 1.5)), zorder=2)
ax.text(382, 31.0, "368 K", fontsize=6, color=INK, va="bottom", rotation=90)
ax.set_xlim(0, 760); ax.set_ylim(36, 0)
ax.set_xlabel("Temperature (K)", labelpad=1); ax.set_ylabel("Depth (cm)", labelpad=1)
ax.yaxis.set_label_position("right"); ax.yaxis.tick_right(); ax.tick_params(length=2, pad=1, labelsize=6)

# ---------- panel letters ----------
# v11: no panel letters (one panel with insets)


# ---------- legend (column by column) ----------
lg = axes_in(PX0, 0.03, 3*PW + 2*PGX, 0.47); lg.axis("off")
h_sz = Patch(fc=SULFZ, ec=SULFZT, lw=0.3, hatch="....", label="sulfidized regolith")
h_res = Patch(fc=RESID, ec=RESIDT, lw=0.3, hatch="//////", label="compacted residue")
h_res.set_hatch_linewidth(HATCH_LW)
leg = lg.legend(handles=[Patch(fc=SULFIDE, ec="none", label="sulfide"), Patch(fc=SULF, ec=SULFT, lw=0.3, label="native sulfur"), h_sz,
                         Patch(fc=REG, ec="none", label="regolith"), h_res, Patch(fc=LAG, ec="none", label="sulfur-free lag"),
                         Patch(fc=PLAINS, ec="none", label="plains"), Patch(fc=HOT, ec="none", label="hot zone"),
                         Patch(fc=HALO, ec="#8fb6d4", lw=0.3, label="halo"),
                         Line2D([], [], color=FRAC, lw=0.8, label="fracture"),
                         Line2D([], [], color=VAP, lw=0.9, marker=">", ms=3, label="S$_2$ vapor"),
                         Line2D([], [], color=INK, lw=0.7, ls=(0, (4, 2)), label="base of charged layer")],
                loc="center", ncol=4, fontsize=6, frameon=False, handlelength=1.8, handleheight=0.9, columnspacing=1.4,
                handletextpad=0.45, labelspacing=0.55, borderaxespad=0)
for hp in leg.get_patches():
    if hp.get_hatch(): hp.set_hatch_linewidth(HATCH_LW)
fig.savefig(OUT, dpi=300, facecolor="white")
fig.savefig(OUT_PDF, facecolor="white")
print("wrote", OUT, OUT_PDF, [round(v, 1) for v in L])
