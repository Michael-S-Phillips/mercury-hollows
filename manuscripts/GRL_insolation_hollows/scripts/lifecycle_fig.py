#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Figure 4: schematic life cycle of a hollow (draft v2). Not to scale; vertical exaggeration differs between rows.
Top row: crater scale (km). Bottom row: hollow scale (tens of meters)."""
import os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Rectangle, ConnectionPatch

OUT = os.environ.get("GRL_FIGS", _REPO + "/manuscripts/GRL_insolation_hollows/figs")
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.0, "text.color": "#1a1a19"})
INK, MUTED = "#1a1a19", "#6b6b67"
PLAINS, PLAINS_EDGE = "#dedbd3", "#a9a59b"
CRUST, CRUST_DARK = "#8f8a80", "#6f6a61"
VOL = "#D4541E"; HOLLOW = "#fbf8ee"; HALO = "#f3e9c6"; LAG = "#b8ad97"; HEAT = "#f2c9b4"
rng = np.random.default_rng(7)

def speckle(ax, poly_xy, n, color=VOL, s=1.6, z=3):
    """Scatter n points inside polygon poly_xy (rejection sampling)."""
    path = mpl.path.Path(poly_xy); x0, y0 = poly_xy.min(0); x1, y1 = poly_xy.max(0)
    pts = np.column_stack([rng.uniform(x0, x1, n * 4), rng.uniform(y0, y1, n * 4)])
    pts = pts[path.contains_points(pts)][:n]
    ax.scatter(pts[:, 0], pts[:, 1], s=s, color=color, lw=0, zorder=z)

TIMECOL = "#8a3d1c"
def panel(ax, title, sub, time=None, tag=None):
    ax.set_axis_off()
    ax.text(0.0, 1.02, title, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")
    t = ax.text(0.0, -0.03, sub, transform=ax.transAxes, fontsize=6.4, color=MUTED, va="top", ha="left", linespacing=1.25)
    last = t
    if time:
        last = ax.annotate("Time: " + time, xy=(0, 0), xycoords=t, xytext=(0, -2), textcoords="offset points", fontsize=6.2, color=TIMECOL,
                           va="top", ha="left", linespacing=1.2, annotation_clip=False)
    if tag:
        ax.annotate(tag, xy=(0, 0), xycoords=last, xytext=(0, -2), textcoords="offset points", fontsize=6.0, color=INK,
                    va="top", ha="left", annotation_clip=False, fontstyle="italic")

def sun(ax, x, y, r, n=8, rays=True):
    ax.add_patch(Circle((x, y), r, color="#f0a43a", zorder=6, lw=0))
    if rays:
        for a in np.linspace(0, 2 * np.pi, n, endpoint=False):
            ax.plot([x + 1.35 * r * np.cos(a), x + 1.9 * r * np.cos(a)], [y + 1.35 * r * np.sin(a), y + 1.9 * r * np.sin(a)], color="#f0a43a", lw=.9, zorder=6)

fig = plt.figure(figsize=(7.0, 5.6))
# ------------------------------------------------------------------ top row (crater scale)
A1 = fig.add_axes([0.02, 0.625, 0.27, 0.30]); A2 = fig.add_axes([0.335, 0.625, 0.645, 0.30])

# Stage 1: source
ax = A1; ax.set_xlim(0, 10); ax.set_ylim(0, 10)
top = 8.0; base_pl = 5.6
ax.add_patch(Polygon([[0, top], [10, top], [10, base_pl], [0, base_pl]], closed=True, fc=PLAINS, ec=PLAINS_EDGE, lw=.6, zorder=1))
crust = np.array([[0, base_pl], [10, base_pl], [10, 0.4], [0, 0.4]])
ax.add_patch(Polygon(crust, closed=True, fc=CRUST, ec=CRUST_DARK, lw=.6, zorder=1)); speckle(ax, crust, 260)
ax.text(5, (top + base_pl) / 2, "younger volcanic plains\n(2.5–3.5 km thick in Caloris)", ha="center", va="center", fontsize=6.3, linespacing=1.15)
ax.text(5, 3.0, "older low-reflectance crust\nwith a volatile-bearing phase", ha="center", va="center", fontsize=6.3, color="white", linespacing=1.15,
        bbox=dict(boxstyle="round,pad=0.25", fc=CRUST_DARK, ec="none", alpha=.85), zorder=5)
panel(ax, "1  Source", "Volatile-bearing crust lies\nbeneath younger plains.", time="plains emplaced ≥ 2.5 Gyr ago")

# Stage 2: exposure (small crater confined to plains; complex crater cutting through them)
ax = A2; ax.set_xlim(0.3, 25.2); ax.set_ylim(0, 12.0)
top, base_pl = 8.0, 5.6
# heat glow beneath the large crater
# plains and crust with the craters cut out: build surface profile
xs = np.linspace(0, 25.2, 1009)
FLAP_H, FLAP_W = 2.1, 2.9
def surf(x):
    y = np.full_like(x, top)
    # small simple crater at x=3.7 (radius 1.6, depth 1.3) confined to plains
    m = np.abs(x - 3.0) < 1.6; y[m] = top - 1.3 * (1 - ((x[m] - 3.0) / 1.6) ** 2)
    rim = (np.abs(x - 3.0) >= 1.6) & (np.abs(x - 3.0) < 2.4); y[rim] += .25 * (1 - (np.abs(x[rim] - 3.0) - 1.6) / .8)
    # complex crater centered at 15.5, radius 7: rim, terraced walls, flat floor at 3.6, central peak
    c, R = 15.5, 7.0; d = np.abs(x - c)
    m = d < R
    wall = m & (d > 4.2)
    y[wall] = 3.6 + (top + FLAP_H - 3.6) * ((d[wall] - 4.2) / (R - 4.2)) ** 1.15
    y[wall] -= .25 * ((np.sin((d[wall] - 4.2) / (R - 4.2) * 3 * np.pi)) > .6)       # terraces
    fl = m & (d <= 4.2); y[fl] = 3.6
    pk = d < 1.5; y[pk] = 3.6 + 2.6 * np.cos(d[pk] / 1.5 * np.pi / 2) ** 1.4
    rim = (d >= R) & (d < R + FLAP_W); y[rim] = top + FLAP_H * (1 - (d[rim] - R) / FLAP_W) ** 1.3
    return y
ys = surf(xs)
sub = Polygon(np.vstack([np.column_stack([xs, ys]), [[25.2, 0.4], [0, 0.4]]]), closed=True, fc="none", ec="none"); ax.add_patch(sub)
for rr_, al in [(5.2, .22), (4.0, .3), (2.8, .4)]:
    e_ = mpl.patches.Ellipse((15.5, 3.0), 2 * rr_ * 1.5, 2 * rr_ * .6, color=HEAT, alpha=al, lw=0, zorder=2.5); ax.add_patch(e_); e_.set_clip_path(sub)
# crust region: below min(surface, base_pl) everywhere except where the peak uplifts crust above the floor
crust_top = np.minimum(ys, base_pl)
d = np.abs(xs - 15.5); uplift = d < 2.4
crust_top[uplift] = np.minimum(ys[uplift], base_pl + 0.0)
crust_top[d < 1.5] = ys[d < 1.5]                                 # central peak is uplifted older crust
crust_poly = np.vstack([np.column_stack([xs, crust_top]), [[25.2, 0.4], [0, 0.4]]])
plains_poly = np.vstack([np.column_stack([xs, ys]), np.column_stack([xs[::-1], crust_top[::-1]])])
ax.add_patch(Polygon(plains_poly, closed=True, fc=PLAINS, ec="none", zorder=1))
ax.add_patch(Polygon(crust_poly, closed=True, fc=CRUST, ec="none", zorder=1.5))
speckle(ax, crust_poly, 900)
ax.plot(xs, ys, color=INK, lw=.8, zorder=4)
# overturned flap at both rims: inverted plains below, older crust (from depth) on top
for side in (-1, 1):
    dd_ = side * (xs - 15.5)
    sel = (dd_ >= 7.0 - 0.35) & (dd_ < 7.0 + FLAP_W) & (ys > top + 1e-6)
    xf, yf = xs[sel], ys[sel]
    split = top + 0.42 * (yf - top)
    lower = np.vstack([np.column_stack([xf, split]), np.column_stack([xf[::-1], np.full(len(xf), top)])])
    upper = np.vstack([np.column_stack([xf, yf]), np.column_stack([xf[::-1], split[::-1]])])
    ax.add_patch(Polygon(lower, closed=True, fc=PLAINS, ec="none", zorder=3.2))
    ax.add_patch(Polygon(upper, closed=True, fc=CRUST, ec="none", zorder=3.3))
    speckle(ax, upper, 70, z=3.4)
    hx = 15.5 + side * 6.75; ax.annotate("", xy=(hx + side * 1.15, top + 1.0), xytext=(hx - side * 0.05, top + 0.25),
        arrowprops=dict(arrowstyle="-|>", color=INK, lw=.6, mutation_scale=5, connectionstyle=f"arc3,rad={-0.5 * side}"), zorder=6)
    ax.plot(xf, split, color=CRUST_DARK, lw=.5, zorder=3.5)
    ax.plot([xf.min(), xf.max()], [top, top], color=PLAINS_EDGE, lw=.5, ls=(0, (2, 1.5)), zorder=3.5)
# melt sheet on the crater floor
fl = np.abs(xs - 15.5) <= 4.2; pk = np.abs(xs - 15.5) < 1.5
for side in (-1, 1):
    sel = fl & ~pk & (np.sign(xs - 15.5) == side)
    ax.fill_between(xs[sel], ys[sel] - .45, ys[sel], color="#5d5852", lw=0, zorder=3)
# hollows: bright marks on the peak flanks, the lower walls, and the floor near the peak
def mark(x, w=.38, h=.13):
    y = np.interp(x, xs, ys); ax.add_patch(Rectangle((x - w / 2, y - h), w, h + .02, fc=HOLLOW, ec="#c9b98a", lw=.4, zorder=5))
for x in (14.55, 14.95, 16.05, 16.45, 11.0, 20.0, 12.6, 18.5, 8.3, 22.7):
    mark(x)
# fracture / vapor arrows from the heated zone
for x0 in (13.2, 15.5, 17.8):
    ax.add_patch(FancyArrowPatch((x0, 1.4), (x0 + (15.5 - x0) * .25, 3.2), arrowstyle="-|>", mutation_scale=6, color="#b5532a", lw=.7, ls=(0, (2, 1.5)), zorder=4))
ax.text(3.0, 8.75, "small crater stays\nwithin the plains:\nno hollows", ha="center", va="bottom", fontsize=6.2, linespacing=1.1)
ax.text(13.6, 10.55, "larger crater cuts through the plains: hollows on the\nuplifted peak, lower walls, and overturned rim", ha="center", va="bottom", fontsize=6.2, linespacing=1.1)
ax.text(25.0, 3.0, "overturned flap:\nolder crust above\ninverted plains", ha="right", va="center", fontsize=6.0, color="white", linespacing=1.1,
        bbox=dict(boxstyle="round,pad=0.25", fc=CRUST_DARK, ec="none", alpha=.9), zorder=8)
ax.annotate("", xy=(23.3, 8.75), xytext=(23.3, 4.25), arrowprops=dict(arrowstyle="-|>", color=INK, lw=.7, mutation_scale=7), zorder=8)
ax.text(6.0, 2.3, "impact heat drives the volatile\nout of the breccia; the vapor\nrecondenses as a pore-filling\ncement in cold rock near the surface", ha="center", va="center", fontsize=6.0, color="white", linespacing=1.1,
        bbox=dict(boxstyle="round,pad=0.25", fc="#8a3d1c", ec="none", alpha=.85), zorder=6)
ax.text(17.15, 5.35, "3–6", ha="left", va="bottom", fontsize=6.0, color=INK, zorder=7)
# zoom box on the peak flank
zb = Rectangle((16.0, 4.15), 1.0, 1.25, fill=False, ec=INK, lw=.7, ls=(0, (2, 1)), zorder=7); ax.add_patch(zb)
panel(ax, "2  Exposure", "Impacts that penetrate the plains expose older crust on walls, on central peaks (the deepest material),\nand on rims, where the overturned flap places deep material on top.", time="impact in minutes; heat lasts 10³–10⁶ yr; hosts of active hollows are < 0.3 to > 1.7 Gyr old")

# ------------------------------------------------------------------ bottom row (hollow scale): collapse model on the peak flank
# The dashed box in panel 2 sits on the sloping flank of the central peak, so panels 3-6 show the same sloping ground.
# Geology is drawn in "flat" coordinates and sheared onto the slope; vertical lines (H, compaction) stay vertical.
from matplotlib.transforms import Affine2D
from matplotlib.patches import Ellipse
B = [fig.add_axes([0.02 + i * 0.245, 0.185, 0.225, 0.32]) for i in range(4)]
YTOP = 12.8                                    # data units; 1 unit = 0.14 in vertically, 0.1575 in horizontally
BRECCIA, CLOSED, ROOF = "#bcb6aa", "#7c766c", CRUST
SLOPE = 0.4                                    # ground falls 0.4 units per unit x (about 20 degrees on the page)
SURF, ROOF_B, H_B = 7.7, 7.05, 2.9             # surface, base of cemented roof, base of void-bearing column (depth H; not to scale)
SKY_TXT = 12.4
x = np.linspace(0, 10, 801)

def shear(ax):
    return Affine2D(np.array([[1, 0, 0], [-SLOPE, 1, 5 * SLOPE], [0, 0, 1]])) + ax.transData

def ground(xx, yy):
    """Page position (data coords) of a point given in flat coords."""
    return xx, yy - SLOPE * (xx - 5)

def rings(ax, T, poly_xy, n, s=3.2, z=1.4):
    path = mpl.path.Path(poly_xy); x0, y0 = poly_xy.min(0); x1, y1 = poly_xy.max(0)
    pts = np.column_stack([rng.uniform(x0, x1, n * 4), rng.uniform(y0, y1, n * 4)])
    pts = pts[path.contains_points(pts)][:n]
    ax.scatter(pts[:, 0], pts[:, 1], s=s, facecolor="#f7f4ec", edgecolor=VOL, linewidths=.35, zorder=z, transform=T)

def dots(ax, T, poly_xy, n, s=1.2, z=1.6):
    path = mpl.path.Path(poly_xy); x0, y0 = poly_xy.min(0); x1, y1 = poly_xy.max(0)
    pts = np.column_stack([rng.uniform(x0, x1, n * 4), rng.uniform(y0, y1, n * 4)])
    pts = pts[path.contains_points(pts)][:n]
    ax.scatter(pts[:, 0], pts[:, 1], s=s, color=VOL, lw=0, zorder=z, transform=T)

def runs(mask):
    idx = np.where(mask)[0]
    return [r for r in (np.split(idx, np.where(np.diff(idx) > 1)[0] + 1) if len(idx) else []) if len(r) >= 3]

G_ROOF = 1.9                                   # depth over which the cement fraction falls from full to zero
CEMENT_EDGE = "#9c3a12"
def mix(c1, c2, f):
    a1, a2 = np.array(mpl.colors.to_rgb(c1)), np.array(mpl.colors.to_rgb(c2)); return tuple((1 - f) * a1 + f * a2)

def column(ax, prof, compacted, roof_mask=None, strip_w=0.0):
    """Void-bearing breccia down to H_B whose pores near the surface are filled with recondensed volatile (the cemented roof),
    grading to open voids with depth; void-free breccia below H_B and under the hollow."""
    T = shear(ax)
    ax.add_patch(Rectangle((-2, -6), 14, H_B + 6, fc=CLOSED, ec="none", zorder=1, transform=T))
    ok = ~compacted if roof_mask is None else roof_mask
    ax.fill_between(x, H_B, prof, where=ok, color=BRECCIA, lw=0, zorder=1.2, transform=T, interpolate=True)
    ax.fill_between(x, H_B, prof, where=~ok, color=CLOSED, lw=0, zorder=1.2, transform=T, interpolate=True)
    # background tint: cemented (crust color) at the surface, fading to uncemented breccia by depth G_ROOF
    NG = 24
    for k in range(NG):
        d0, d1 = k * G_ROOF / NG, (k + 1) * G_ROOF / NG
        f = (k + .5) / NG
        ax.fill_between(x, prof - d1, prof - d0, where=ok, color=mix(ROOF, BRECCIA, f ** .8), lw=0, zorder=1.25, transform=T,
                        interpolate=True)
    # pores: filled with cement near the surface, open below; the filled fraction falls linearly with depth
    for r_ in runs(ok):
        xa, xb = x[r_].min() + .08, x[r_].max() - .08
        if xb - xa < .2: continue
        n = int(62 * (xb - xa))
        pts = np.column_stack([rng.uniform(xa, xb, n), rng.uniform(H_B + .12, SURF - .1, n)])
        top = np.interp(pts[:, 0], x, prof) - .1
        pts = pts[pts[:, 1] < top]; depth = np.interp(pts[:, 0], x, prof) - pts[:, 1]
        pfill = np.clip(1 - depth / G_ROOF, 0, 1)
        if strip_w > 0 and (~ok).any():                     # sunlight has removed the cement next to the scarps
            xe = x[~ok]; dist = np.abs(pts[:, [0]] - xe[None, :]).min(1)
            pfill = pfill * np.clip(dist / strip_w, 0, 1)
        filled = rng.random(len(pts)) < pfill
        ax.scatter(pts[filled, 0], pts[filled, 1], s=3.2, facecolor=VOL, edgecolor=CEMENT_EDGE, linewidths=.3, zorder=1.4, transform=T)
        ax.scatter(pts[~filled, 0], pts[~filled, 1], s=3.2, facecolor="#f7f4ec", edgecolor=VOL, linewidths=.35, zorder=1.4, transform=T)
    ax.plot(x, prof, color=INK, lw=.8, zorder=4, transform=T)
    return T

def hdim(ax, T, xx=9.45):
    ax.annotate("", xy=(xx, H_B), xycoords=T, xytext=(xx, SURF), textcoords=T,
                arrowprops=dict(arrowstyle="<->", color=INK, lw=.6, mutation_scale=5, shrinkA=0, shrinkB=0), zorder=6)
    ax.text(xx - .18, (H_B + ROOF_B) / 2, "H", ha="right", va="center", fontsize=6.4, color=INK, zorder=6, transform=T,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=.85))

def halo(ax, T, xc, half):
    ax.plot([xc - half, xc + half], [SURF + .05, SURF + .05], color="#e6c76e", lw=2.2, solid_capstyle="butt", zorder=5, transform=T)

def meteor(ax, xh, yh, size=1.0, ang=38, ux=0.1575, uy=0.14):
    """Small impactor: head at (xh, yh) in data coords, trail toward the upper right; ux, uy = inches per data unit."""
    r = .055 * size                                              # head radius, inches
    ax.add_patch(Ellipse((xh, yh), 2 * r / ux, 2 * r / uy, color="#8a3d1c", lw=0, zorder=7))
    a = np.deg2rad(ang); dx, dy = np.cos(a), np.sin(a); px, py = -dy, dx
    for off, L, lw, al in [(0, .24, .9, .9), (.028, .17, .6, .7), (-.028, .17, .6, .7)]:
        x0 = (xh * ux + off * size * px + r * dx) / ux; y0 = (yh * uy + off * size * py + r * dy) / uy
        x1 = x0 + L * size * dx / ux; y1 = y0 + L * size * dy / uy
        ax.plot([x0, x1], [y0, y1], color="#8a3d1c", lw=lw, alpha=al, solid_capstyle="round", zorder=7)

def rot_label(ax, T, xx, yy, s, **kw):
    return ax.text(xx, yy, s, transform=T, transform_rotates_text=True, rotation=0, rotation_mode="anchor", **kw)

def setup(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, YTOP)

# Phase 3: initiation (sunlight strips the roof cement; the roof fails over a small area)
ax = B[0]; setup(ax)
prof = np.full_like(x, SURF); m = np.abs(x - 6.0) < .7
prof[m] = SURF - 0.9 * np.clip((0.7 - np.abs(x[m] - 6.0)) / .25, 0, 1)
T = column(ax, prof, m, strip_w=1.0)
for xx in (4.75, 7.25): halo(ax, T, xx, .4)
sun(ax, 8.7, 11.4, .42)
for k in range(3):
    ax.annotate("", xy=(6.6 - k * .55, SURF + .2), xycoords=T, xytext=(8.2 - k * .6, 10.7), textcoords="data",
                arrowprops=dict(arrowstyle="-|>", color="#e39a2d", lw=.7, mutation_scale=6))
ax.text(.2, SKY_TXT, "sunlight sublimates the\ncement near the scarp;\nthe roof fails", fontsize=6.2, va="top", linespacing=1.1)
rot_label(ax, T, .25, SURF - .38, "cemented roof", fontsize=5.6, color="white", va="center", zorder=7,
          bbox=dict(boxstyle="round,pad=0.12", fc=ROOF, ec="none", alpha=.75))
rot_label(ax, T, .25, 5.45, "devolatilized breccia\nwith open voids", fontsize=5.6, color=INK, va="center", linespacing=1.05, zorder=7,
          bbox=dict(boxstyle="round,pad=0.15", fc="#e4e0d8", ec="none", alpha=.92))
rot_label(ax, T, .25, H_B - 1.05, "voids crushed: overburden\nexceeds framework strength", fontsize=5.6, color="white", va="center",
          linespacing=1.05, zorder=7)
hdim(ax, T)
panel(ax, "3  Initiation", "Loss of roof cement begins where\nsunlight is strongest.", time="can begin > 10⁸ yr after\nthe impact", tag="Deutsch et al. Stage 1")

# Phase 4: growth (the roof fails at the scarps; the column beneath compacts; the floor drops by phi*H)
ax = B[1]; setup(ax)
prof = np.full_like(x, SURF); m = np.abs(x - 5) < 2.9
prof[m] = SURF - 0.9 * np.clip((2.9 - np.abs(x[m] - 5)) / .35, 0, 1)
T = column(ax, prof, m, strip_w=1.0)
ax.fill_between(x[m], prof[m], prof[m] + .1, color=HOLLOW, zorder=3, lw=0, transform=T)
for xx in (1.15, 8.85): halo(ax, T, xx, .8)
for sgn in (-1, 1):
    ax.annotate("", xy=(5 + sgn * 3.9, SURF - .33), xycoords=T, xytext=(5 + sgn * 2.95, SURF - .33), textcoords=T,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=.8, mutation_scale=7), zorder=6)
for xx in (3.9, 5.0, 6.1):
    ax.annotate("", xy=(xx, H_B + .5), xycoords=T, xytext=(xx, SURF - 1.4), textcoords=T,
                arrowprops=dict(arrowstyle="-|>", color="white", lw=.7, mutation_scale=6), zorder=6)
ax.annotate("floor drops by d = φH", xy=(5.6, SURF - .85), xycoords=T, xytext=(6.2, 10.0), textcoords="data", ha="center", va="bottom",
            fontsize=6.0, color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=.5, shrinkA=1, shrinkB=1), zorder=6)
ax.text(.2, SKY_TXT, "bright scarps retreat as\nthe roof fails", fontsize=6.2, va="top", linespacing=1.1)
hdim(ax, T)
panel(ax, "4  Growth", "The column compacts; overburden\nsets H, the same everywhere.", time="retreat ≥ 1 cm per 10⁴ yr;\nfloors ≤ ~10⁵ yr old", tag="Deutsch et al. Stage 1")

# Phase 5: dormancy (lag mantles floor and scarps; halo fades; a small impact can restart retreat)
ax = B[2]; setup(ax)
prof = np.full_like(x, SURF); m = np.abs(x - 5) < 3.3
prof[m] = SURF - 0.9 * np.clip((3.3 - np.abs(x[m] - 5)) / .75, 0, 1)
T = column(ax, prof, m)
lag_top = np.minimum(prof + .22, SURF)
ax.fill_between(x[m], prof[m], lag_top[m], color=LAG, zorder=3.5, lw=0, hatch="....", edgecolor="#9a8f78", transform=T)
ax.plot(x[m], lag_top[m], color="#7d735f", lw=.5, zorder=4, transform=T)
mx, my = ground(7.6, SURF); meteor(ax, mx + .05, my + 1.05, size=1.0)
ax.text(.2, SKY_TXT, "lag mantles floor and scarps;\na small impact can breach\nthe roof again", fontsize=6.2, va="top", linespacing=1.1)
hdim(ax, T)
panel(ax, "5  Dormancy", "Retreat stops where the roof holds\nor a lag forms; the halo fades.", time="~0.7 × the growth time; halos\nfade within ~3 × 10⁸ yr", tag="Deutsch et al. Stage 2")

# Phase 6: fossilization (softened, darkened, cratered)
ax = B[3]; setup(ax)
prof = np.full_like(x, SURF); m = np.abs(x - 5) < 3.7
prof[m] = SURF - 0.75 * np.cos((x[m] - 5) / 3.7 * np.pi / 2) ** 2
for cxx, rr_ in [(2.0, .35), (6.3, .45), (8.4, .3)]:
    mm_ = np.abs(x - cxx) < rr_; prof[mm_] -= .3 * (1 - ((x[mm_] - cxx) / rr_) ** 2)
T = column(ax, prof, m)
for (gx, dy_, sz) in [(5.6, 1.5, .75), (7.6, 2.3, .65), (8.9, 1.2, .8)]:
    mx, my = ground(gx, SURF); meteor(ax, mx, my + dy_, size=sz)
ax.text(.2, SKY_TXT, "softened scarps,\nsuperposed craters,\nno reflectance contrast", fontsize=6.2, va="top", linespacing=1.1)
hdim(ax, T)
panel(ax, "6  Fossilization", "Impact gardening and space\nweathering erase the hollow.", time="diffusion alone needs > 4 Gyr\nto erase 20–30 m of relief", tag="Deutsch et al. Stage 3")

# ------------------------------------------------------------------ key to the symbols (two rows of three)
K = fig.add_axes([0.02, 0.0, 0.96, 0.058]); K.set_axis_off(); K.set_xlim(0, 100); K.set_ylim(0, 10)
kx = [0.5, 34.0, 67.5]; ky = [7.2, 2.4]; tx = 3.6
def key_text(i, j, s): K.text(kx[i] + tx, ky[j], s, fontsize=6.2, va="center", ha="left")
# row 1
K.scatter([kx[0] + .9, kx[0] + 1.8, kx[0] + 2.5, kx[0] + 1.4], [ky[0] + .7, ky[0] - .6, ky[0] + .5, ky[0] + .1], s=3, color=VOL, lw=0)
key_text(0, 0, "volatile-bearing phase")
K.scatter([kx[1] + .9, kx[1] + 1.9, kx[1] + 1.4], [ky[1 - 1] + .6, ky[0] + .4, ky[0] - .7], s=6, facecolor="#f7f4ec", edgecolor=VOL, linewidths=.45)
key_text(1, 0, "void where the volatile was lost")
K.scatter([kx[2] + .9, kx[2] + 1.9, kx[2] + 1.4], [ky[0] + .6, ky[0] + .4, ky[0] - .7], s=6, facecolor=VOL, edgecolor=CEMENT_EDGE, linewidths=.4)
key_text(2, 0, "pore refilled by recondensed volatile")
# row 2
K.add_patch(Rectangle((kx[0] + .3, ky[1] - 1.6), 2.8, 3.2, fc=CLOSED, ec="none"))
key_text(0, 1, "voids closed by compaction or overburden")
K.plot([kx[1] + .3, kx[1] + 3.1], [ky[1], ky[1]], color="#e6c76e", lw=2.4, solid_capstyle="butt")
key_text(1, 1, "bright halo")
meteor(K, kx[2] + 1.0, ky[1] - .5, size=.55, ux=0.96 * 7.0 / 100, uy=0.058 * 5.6 / 10)
key_text(2, 1, "small impact")

fig.savefig(f"{OUT}/fig4.png", dpi=400); fig.savefig(f"{OUT}/fig4.pdf"); print("wrote fig4")
