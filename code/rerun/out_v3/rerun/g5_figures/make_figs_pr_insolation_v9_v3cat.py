# g5 copy of code/rerun/out_v3/rerun/g2_insolation/make_figs_pr_insolation_v9.py (itself the g2 copy of manuscripts/PlanetaryResearch_hollows/figscripts/make_figs_pr_insolation_v9.py), 6 Nov 2026. Changed: inputs still read from the g2 folder (G2_CAT=v3); output to g5_figures/figS5_insolation_v9_v3cat.png/.pdf (the main-text Figure 3 crop is no longer written); figure height 5.3 -> 7.95 in with panels a-d unchanged in size and spacing; new panel e (hosting vs Delta T_w: binned points, linear fit, best single step), model code from g2_insolation/fig_step_vs_linear.py, record in figS5_panel_e_v3cat_record.json.
#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/figscripts/make_figs_pr_insolation_v9.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""v9 (2 Nov 2026), copied from make_figs_pr_insolation_v8.py. One change: panel a shows hollow fields within 15 deg of the
equator (was 30 deg), to match the values quoted in Section 3.4, read from
analysis_workshop/insolation_v6/cos2lon_fit_lat15.json (made by cos2lon_lat15.py there: same null, exclusions, 30-deg bins,
weighted fit of a + b cos 2 lon, and 2,000-resample bootstrap interval on b as hotwarm.py). Panels b-d are unchanged.
Writes figures/figS5_insolation_v9.png/.pdf and figures/fig3_insolation_v9.png/.pdf. Run without arguments.

v8 docstring follows.
v8 (30 Oct 2026), copied from make_figs_pr.py, the script that made figures/fig2_insolation_v3.png (Figure S5; its top half
is main-text Figure 3, fig3_insolation_v6.png). Changes, and nothing else: (1) the two labels that read "ΔT" in panel b now
read "Δ$T_\mathrm{w}$" (mean annual surface temperature minus that at the warm longitude of the same latitude); (2) default
data paths point to the hollows-emissivity tree on this machine; (3) the script stops after the insolation figure and writes
figures/figS5_insolation_v8.png/.pdf (all four panels) and figures/fig3_insolation_v8.png/.pdf (panels a-b: the top half of
the figure, as fig3_insolation_v6.py cropped it). Run without arguments.

Original docstring follows.
Planetary Research version (15 cm wide). Figures 1, 2, and 4 of the paper, from the GRL figure script on the setting and depth of hollows on Mercury (draft v2).
Figure 4 (life-cycle schematic) is drawn by lifecycle_fig.py.

Fig 1  insolation: (a) longitude profile of hollow fields, (b) hot/warm odds of hosting vs crater size,
                   (c) equator-facing fraction of wall hollows vs latitude, (d) noon temperature of flat ground and walls
Fig 2  exposure of older crust: (a) north polar map with the smooth-plains mask, (b) field density by latitude in and out of
       the plains, (c) hollow density on peaks, floors, and walls by crater setting, (d) Caloris interior craters vs
       isolated craters, against crater diameter, with the diameter at which excavation reaches the base of the plains
Fig 3  depth: (a) median depth by latitude band against mean annual temperature, (b) median depth per host crater vs size

Inputs (environment variables give the folders):
  GRL_DIR      data/test2_craters_h18.csv, data/aspect_wall_hollows.csv, records/hotwarm_record.json
  HOLLOWS_DIR  code/rerun/out_v2/coverage_10ppd.npz, fields.csv; peak_wall_test/peak_wall_by_setting.json,
               peak_wall_test/caloris_test_record.json
  GRL_SCRIPTS  nsp_laea.npz, latcorr_density_record.json, latcorr_inc60_75_record.json
  GRL_WORK     data/test1_rows_inc.csv and records/depth_reanalysis_record.json (from depth_reanalysis.py)
Writes figs/fig1-3 (.png, .pdf) and figs/figure_numbers_v2.json.
"""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g2_insolation");  # g5: inputs stay in the g2 folder
G5 = _g2_os.path.dirname(_g2_os.path.abspath(__file__));  # g5: output folder
G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_GRLDATA = f"{G2_ROOT}/manuscripts/GRL_insolation_hollows/data"
# --- end g2 block ---
import json, os
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import pyproj
from scipy.stats import beta
from scipy.spatial import cKDTree

_ROOT = G2_ROOT   # v8: hollows-emissivity
G = os.environ.get("GRL_DIR", f"{_ROOT}/manuscripts/GRL_insolation_hollows")
U = os.environ.get("HOLLOWS_DIR", _ROOT)
S = os.environ.get("GRL_SCRIPTS", "/home/claude/grl/scripts")
W = os.environ.get("GRL_WORK", "/home/claude/grl")
OUT = G5   # g5: output folder (was GRL_FIGS or G2_OUT)
PW = f"{U}/peak_wall_test"

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5,
    "axes.linewidth": .6, "axes.edgecolor": "#5a5a57", "xtick.color": "#5a5a57", "ytick.color": "#5a5a57",
    "xtick.major.width": .6, "ytick.major.width": .6, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.labelcolor": "#1a1a19", "text.color": "#1a1a19", "axes.labelsize": 7.5,
    "figure.facecolor": "white", "axes.facecolor": "white"})
INK, MUTED, GRID = "#1a1a19", "#6b6b67", "#dcdcd8"
BLUE, ORANGE, GREEN = "#1F6FB2", "#D4541E", "#178f66"     # validated together (all pairs, light surface)
HOTTXT, WARMTXT = "#a8431a", "#4a6b8a"
NUM = {}

def tidy(a, letter, grid="both"):
    if grid: a.grid(True, axis=grid if grid != "both" else "both", color=GRID, lw=.5, zorder=0)
    a.set_axisbelow(True)
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
    a.text(-0.02, 1.03, letter, transform=a.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=400); fig.savefig(f"{OUT}/{name}.pdf"); plt.close(fig); print("wrote", name)

def cp(k, n, a=.05):
    lo = beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi

R_M = 2439.7
# ===================================================================================== Figure 1: insolation
# g5: one row added below for panel e; the figure grows from 5.3 to 7.95 in, and panels a-d keep their size and spacing in inches
FH = 5.3 + 2.65; _p = lambda x, y, w, h: [x, (y * 5.3 + 2.65) / FH, w, h * 5.3 / FH]
fig = plt.figure(figsize=(5.9, FH))
a1 = fig.add_axes(_p(0.075, 0.585, 0.40, 0.36)); a2 = fig.add_axes(_p(0.585, 0.585, 0.40, 0.36))
b1 = fig.add_axes(_p(0.075, 0.085, 0.40, 0.36)); b2 = fig.add_axes(_p(0.585, 0.085, 0.40, 0.36))
e1 = fig.add_axes([0.075, 0.085 * 5.3 / FH, 0.40, 0.36 * 5.3 / FH])

# (a) longitude profile
hw = json.load(open(f"{G2_OUT}/cos2lon_fit_lat15.json"))["lat15"]   # v9
lon = np.array(hw["lon_bin_centers"]); r = np.array(hw["obs_over_null_ratio"]); xs = np.linspace(0, 360, 361)
for c in (0, 180, 360): a1.axvspan(c - 45, c + 45, color=ORANGE, alpha=.07, lw=0, zorder=0)
for c in (90, 270): a1.axvspan(c - 45, c + 45, color="#6f8fae", alpha=.09, lw=0, zorder=0)
lo_c, hi_c = hw["slope_ci"]
a1.fill_between(xs, hw["intercept"] + lo_c * np.cos(np.radians(2 * xs)), hw["intercept"] + hi_c * np.cos(np.radians(2 * xs)), color=BLUE, alpha=.18, lw=0, zorder=2)
a1.plot(xs, hw["intercept"] + hw["slope"] * np.cos(np.radians(2 * xs)), color=BLUE, lw=1.5, zorder=3)
a1.plot(lon, r, "o", ms=4.2, color=INK, mfc="white", mew=1.1, zorder=5)
a1.axhline(1, color=MUTED, lw=.9, ls=(0, (3, 2)), zorder=1)
a1.set_xlim(0, 360); a1.set_xticks([0, 90, 180, 270, 360]); a1.set_ylim(0, 2.6)
a1.set_xlabel("longitude (°E)"); a1.set_ylabel("observed / expected hollow fields")
for c, t, col in [(0, "hot", HOTTXT), (90, "warm", WARMTXT), (180, "hot", HOTTXT), (270, "warm", WARMTXT)]:
    a1.text(c if c else 4, 2.5, t, ha="center" if c else "left", va="top", fontsize=6.6, color=col)
a1.text(355, 0.12, f"|latitude| ≤ {hw['latmax']:.0f}°, n = {hw['n']} fields", ha="right", fontsize=6.4, color=MUTED)
tidy(a1, "a")

# (b) odds of hosting per 10 K of modeled temperature excess dT, against crater diameter
import sys; sys.path.insert(0, f"{G}/scripts"); import thermal
_z = np.load(f"{G}/data/thermal_grid.npz"); _la, _lo, _T = _z["la"], _z["lo"], _z["G"]
s = pd.read_csv(f"{IN}/test2_craters_h18.csv"); s = s[(s.diameter >= 5) & (s.abslat <= 60)].copy()
s["dT"] = thermal.lookup(_la, _lo, _T, s.latitude.values, s.longitude.values) - thermal.lookup(_la, _lo, _T, s.latitude.values, np.full(len(s), 90.0))
m = smf.logit("hosts ~ dT * logD + abslat", data=s).fit(disp=0)
V = m.cov_params().loc[["dT", "dT:logD"], ["dT", "dT:logD"]].values
Dg = np.logspace(np.log10(5), np.log10(800), 120); X = np.column_stack([np.ones_like(Dg), np.log10(Dg)])
c = m.params["dT"] + m.params["dT:logD"] * np.log10(Dg); se = np.sqrt(np.einsum("ij,jk,ik->i", X, V, X))
a2.axhspan(1, 1e3, color=ORANGE, alpha=.07, lw=0, zorder=0); a2.axhspan(1e-3, 1, color="#6f8fae", alpha=.09, lw=0, zorder=0)
a2.fill_between(Dg, np.exp(10 * (c - 1.96 * se)), np.exp(10 * (c + 1.96 * se)), color=BLUE, alpha=.18, lw=0, zorder=2)
a2.plot(Dg, np.exp(10 * c), color=BLUE, lw=1.5, zorder=4)
bins = [5, 10, 20, 40, 80, 160, 800]
for lo, hi in zip(bins[:-1], bins[1:]):
    d2 = s[(s.diameter >= lo) & (s.diameter < hi)]
    mm = smf.logit("hosts ~ dT + abslat", data=d2).fit(disp=0)
    cb, eb = mm.params["dT"], mm.bse["dT"]; xc = np.sqrt(lo * hi)
    a2.errorbar([xc], [np.exp(10 * cb)], yerr=[[np.exp(10 * cb) - np.exp(10 * (cb - 1.96 * eb))], [np.exp(10 * (cb + 1.96 * eb)) - np.exp(10 * cb)]],
                fmt="o", ms=4.2, color=INK, mfc="white", mew=1.1, lw=.8, capsize=1.8, zorder=5)
    a2.text(xc, np.exp(10 * (cb + 1.96 * eb)) * 1.025, str(int(d2.hosts.sum())), fontsize=6, color=MUTED, ha="center", va="bottom")
a2.axhline(1, color=MUTED, lw=.9, ls=(0, (3, 2)), zorder=3)
a2.set_xscale("log"); a2.set_xlim(3.5, 900); a2.set_ylim(.75, 2.3)
a2.set_xticks([5, 10, 20, 50, 100, 300]); a2.set_xticklabels(["5", "10", "20", "50", "100", "300"]); a2.minorticks_off()
a2.set_yticks([.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2])
a2.set_xlabel("crater diameter (km)"); a2.set_ylabel(r"odds ratio of hosting per 10 K of Δ$T_\mathrm{w}$")
a2.text(3.8, 2.22, r"↑ hosting rises with Δ$T_\mathrm{w}$", fontsize=6.4, color=HOTTXT, va="center")
tidy(a2, "b")
NUM["fig2b_odds"] = dict(n_craters=int(len(s)), n_hosts=int(s.hosts.sum()), inter=float(m.params["dT:logD"]), inter_se=float(m.bse["dT:logD"]),
                    inter_p=float(m.pvalues["dT:logD"]), OR10K_10km=float(np.exp(10 * (m.params["dT"] + m.params["dT:logD"]))),
                    OR10K_100km=float(np.exp(10 * (m.params["dT"] + 2 * m.params["dT:logD"]))),
                    zero_km=float(10 ** (-m.params["dT"] / m.params["dT:logD"])))

# (c) wall aspect
a = pd.read_csv(f"{G2_OUT}/aspect_wall_hollows.csv")
w = a[(a.rho > .55) & (a.rho < 1.05) & (a.host_D_km >= 15)].copy()
cal = (31.9, 161.9)
def gc(lat1, lon1, lat2, lon2):
    p1, p2, l1, l2 = map(np.radians, (lat1, lat2, lon1, lon2))
    return R_M * np.arccos(np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(l1 - l2), -1, 1))
w["dcal"] = gc(w.lat, w.lon, *cal)
rng = np.random.default_rng(20260911)
def stats(sub, nb=4000):
    k = sub.poleward.astype(bool); f = k.mean()
    bt = [rng.choice(k.values, len(k)).mean() for _ in range(nb)]
    pc = sub.groupby("host_uid").poleward.mean().values
    bc = [rng.choice(pc, len(pc)).mean() for _ in range(nb)]
    return dict(n=int(len(sub)), nc=int(len(pc)), f=float(f), lo=float(np.percentile(bt, 2.5)), hi=float(np.percentile(bt, 97.5)),
                fc=float(pc.mean()), clo=float(np.percentile(bc, 2.5)), chi=float(np.percentile(bc, 97.5)))
N = w[w.hemi == "N"]
bands = [(0, 15), (15, 30), (30, 45), (45, 60), (60, 75)]
ST = {f"{lo}-{hi}": stats(N[(N.abslat >= lo) & (N.abslat < hi)]) for lo, hi in bands}
ST["15-30_noCaloris"] = stats(N[(N.abslat >= 15) & (N.abslat < 30) & (N.dcal >= 1200)])
ST["south_all"] = stats(w[w.hemi == "S"])
xc = np.array([7.5, 22.5, 37.5, 52.5, 67.5]); fh = [ST[f"{lo}-{hi}"] for lo, hi in bands]
b1.axhline(.5, color=MUTED, lw=.9, ls=(0, (3, 2)), zorder=1)
b1.errorbar(xc - 1.2, [x["f"] for x in fh], yerr=[[x["f"] - x["lo"] for x in fh], [x["hi"] - x["f"] for x in fh]],
            fmt="o", ms=4.6, color=BLUE, mfc=BLUE, lw=1, capsize=2, zorder=5, label="per hollow")
b1.errorbar(xc + 1.2, [x["fc"] for x in fh], yerr=[[x["fc"] - x["clo"] for x in fh], [x["chi"] - x["fc"] for x in fh]],
            fmt="s", ms=4.2, color=INK, mfc="white", mew=1.1, lw=.8, capsize=2, zorder=5, label="per host crater")
nc = ST["15-30_noCaloris"]
b1.errorbar([22.5 - 1.2], [nc["f"]], yerr=[[nc["f"] - nc["lo"]], [nc["hi"] - nc["f"]]], fmt="D", ms=4.2, color=ORANGE, mfc="white", mew=1.1, lw=1,
            capsize=2, zorder=6, label="per hollow, Caloris excluded")
for x, v in zip(xc, fh): b1.text(x, .17, f"{v['n']}\n{v['nc']}", ha="center", va="bottom", fontsize=5.8, color=MUTED, linespacing=1.1)
b1.text(-9.5, .17, "hollows\ncraters", ha="left", va="bottom", fontsize=5.8, color=MUTED, linespacing=1.1)
b1.set_xlim(-10, 76); b1.set_ylim(.15, 1.0); b1.set_xticks([0, 15, 30, 45, 60, 75])
b1.set_xlabel("latitude (°N)"); b1.set_ylabel("fraction on equator-facing walls")
b1.legend(frameon=False, fontsize=6.6, loc="upper left", handletextpad=.3, borderpad=.2)
tidy(b1, "c")

# (d) noon radiative-equilibrium temperature
sig, A, eps = 5.670e-8, 0.09, 0.95
Sp, Sa = 1361 / 0.3075 ** 2, 1361 / 0.4667 ** 2
def Tn(Sf, inc_deg): return ((1 - A) * Sf * np.clip(np.cos(np.radians(inc_deg)), 0, None) / (eps * sig)) ** .25
lat = np.linspace(0, 85, 171)
for Sf, ls in [(Sp, "-"), (Sa, (0, (4, 2)))]:
    b2.plot(lat, Tn(Sf, np.abs(lat - 30)), ls=ls, color=ORANGE, lw=1.4)
    b2.plot(lat, Tn(Sf, lat), ls=ls, color=INK, lw=1.4)
    tp = Tn(Sf, lat + 30); tp[lat + 30 >= 89.5] = np.nan
    b2.plot(lat, tp, ls=ls, color=BLUE, lw=1.4)
b2.text(60.5, 158, "pole-facing\nwall shadowed\nat noon", fontsize=6.0, color=BLUE, ha="left", va="bottom", linespacing=1.1)
H = [mpl.lines.Line2D([], [], color=ORANGE, lw=1.4, label="30° equator-facing wall"), mpl.lines.Line2D([], [], color=INK, lw=1.4, label="flat ground"),
     mpl.lines.Line2D([], [], color=BLUE, lw=1.4, label="30° pole-facing wall"),
     mpl.lines.Line2D([], [], color=MUTED, lw=1.2, label="hot longitude, perihelion"), mpl.lines.Line2D([], [], color=MUTED, lw=1.2, ls=(0, (4, 2)), label="warm longitude, aphelion")]
b2.legend(handles=H, frameon=False, fontsize=6.1, loc="lower left", handlelength=2.2, borderpad=.2, labelspacing=.25)
b2.set_xlim(0, 85); b2.set_ylim(150, 760); b2.set_xlabel("latitude (°)"); b2.set_ylabel("noon surface temperature (K)")
tidy(b2, "d")
# g5 (6 Nov 2026), new panel e: hosting against Delta T_w, model code and sample of g2_insolation/fig_step_vs_linear.py (top-middle
# panel, v3): craters >= 5 km within 60 deg of the equator, h09, controls logD + abslat + C1; probabilities for the median crater
# from a model with Delta T_w in bins (points, 95% intervals), the linear model and the best single step (threshold.json
# lat60_main; 34 K in v3). The step-plus-linear line of that figure is omitted.
import patsy
CTRL = "logD + abslat + C1"; BINS = [0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 63]
fit = lambda f, d: smf.logit(f, data=d).fit(disp=0, maxiter=200)
def prob_curve(m, ref):
    X = np.asarray(patsy.build_design_matrices([m.model.data.model_spec], ref)[0])
    eta = X @ m.params.values; se = np.sqrt(np.einsum("ij,jk,ik->i", X, m.cov_params().values, X))
    ex = lambda z: 100 / (1 + np.exp(-z)); return ex(eta), ex(eta - 1.96 * se), ex(eta + 1.96 * se)
_c = pd.read_csv(f"{IN}/crater_features.csv", low_memory=False, usecols=["id", "dT", "abslat", "logD", "C1", "h09", "lat", "lon"])
de = _c[_c.abslat <= 60].copy(); TR = json.load(open(f"{G2_OUT}/threshold.json"))["lat60_main"]; tstep = TR["step_dT_t"]
de["dT10"] = de.dT / 10; de["st"] = (de.dT >= tstep).astype(float); de["bin"] = pd.cut(de.dT, BINS, include_lowest=True).astype(str)
med = {k: float(de[k].median()) for k in ("logD", "abslat", "C1")}
mb = fit(f"h09 ~ {CTRL} + C(bin)", de); ml = fit(f"h09 ~ {CTRL} + dT10", de); ms = fit(f"h09 ~ {CTRL} + st", de)
cats = sorted(de["bin"].unique(), key=lambda q: float(q.strip("([").split(",")[0]))
pe, plo, phi = prob_curve(mb, pd.DataFrame({**{k: [v] * len(cats) for k, v in med.items()}, "bin": cats}))
mid = np.array([(float(q.strip("([").split(",")[0]) + float(q.strip("])").split(",")[1])) / 2 for q in cats])
gg = np.linspace(0, 63, 631); rg = pd.DataFrame({**{k: [v] * len(gg) for k, v in med.items()}, "dT10": gg / 10, "st": (gg >= tstep).astype(float)})
e1.errorbar(mid, pe, yerr=[pe - plo, phi - pe], fmt="o", ms=4.2, color=INK, mfc="white", mew=1.1, lw=.8, capsize=1.8, zorder=5, label=r"Δ$T_\mathrm{w}$ bins (95% intervals)")
e1.plot(gg, prob_curve(ml, rg)[0], color=BLUE, lw=1.5, zorder=4, label="steady rise (linear)")
e1.plot(gg, prob_curve(ms, rg)[0], color=ORANGE, lw=1.5, ls=(0, (4, 2)), zorder=4, label=f"single step at {tstep:.0f} K")
e1.set_xlim(0, 65); e1.set_ylim(0, 1.08 * phi.max()); e1.set_xticks([0, 10, 20, 30, 40, 50, 60])
e1.set_xlabel(r"Δ$T_\mathrm{w}$ (K)"); e1.set_ylabel("hosting probability (%)")
tidy(e1, "e")
_aic = TR["AIC"]; dA = _aic["step_dT"] - _aic["linear_dT"]
_h, _l = e1.get_legend_handles_labels()
_o = [_l.index(x) for x in sorted(_l, key=lambda q: 0 if "bins" in q else 1 if "linear" in q else 2)]
fig.legend([_h[i] for i in _o], [_l[i] for i in _o], frameon=False, fontsize=6.6, loc="upper left", bbox_to_anchor=(0.585, (0.085 * 5.3 + 1.908) / FH), handletextpad=.4, borderpad=.2)
fig.text(0.59, (0.085 * 5.3 + 1.42) / FH, f"{TR['hosts']} hosts among {TR['n']:,} craters ≥ 5 km\nwithin 60° of the equator; median crater\n"
         f"(D = {10 ** med['logD']:.1f} km, |latitude| {med['abslat']:.0f}°, C = {med['C1']:.3f})\n"
         f"AIC of the step {abs(dA):.1f} {'lower' if dA < 0 else 'higher'} than linear", fontsize=6.4, color=MUTED, va="top", ha="left", linespacing=1.25)
json.dump(dict(source="g2_insolation/fig_step_vs_linear.py (v3, top-middle panel), step-plus-linear omitted", catalog=G2_CAT, t_step=tstep, hosts=TR["hosts"], n=TR["n"],
               dAIC_step_minus_linear=round(dA, 1), median_controls=med, D_median_km=10 ** med["logD"],
               bins=[dict(bin=q, p=float(x), lo=float(y), hi=float(z)) for q, x, y, z in zip(cats, pe, plo, phi)],
               linear_OR_per10K=float(np.exp(ml.params["dT10"])), step_OR=float(np.exp(ms.params["st"]))),
          open(f"{OUT}/figS5_panel_e_v3cat_record.json", "w"), indent=1)
# g5: write Figure S5 only (main-text Figure 3 is g2's fig3_insolation_v9.png, copied unchanged)
_names = [f"{OUT}/{n}" for n in ("figS5_insolation_v9_v3cat.png", "figS5_insolation_v9_v3cat.pdf")]
assert _g2_os.environ.get("G5_REDRAW") == "1" or not any(os.path.exists(n) for n in _names), "outputs exist; set G5_REDRAW=1 to redraw this script's own outputs"
fig.savefig(_names[0], dpi=400); fig.savefig(_names[1])
plt.close(fig); print("wrote", *_names, "panel e: t", tstep, "dAIC", round(dA, 1), "hosts", TR["hosts"])
raise SystemExit(0)

# ===================================================================================== Figure 2: exposure of older crust
fig = plt.figure(figsize=(5.9, 6.5))
c1 = fig.add_axes([0.035, 0.55, 0.42, 0.43]); c2 = fig.add_axes([0.575, 0.585, 0.41, 0.36])
d2 = fig.add_axes([0.075, 0.075, 0.41, 0.36]); d1 = fig.add_axes([0.575, 0.075, 0.41, 0.36])   # c = Caloris, d = terrain
# (a) polar map
cov = np.load(f"{U}/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"].astype(float)
f = pd.read_csv(f"{U}/code/rerun/out_v2/fields.csv")
nm = np.load(f"{S}/nsp_laea.npz"); msk = nm["m40"]; pix = float(nm["pix"]); half = float(nm["half"])
R = 2439400.0
ps = pyproj.CRS.from_proj4(f"+proj=laea +lat_0=90 +lon_0=0 +R={R} +no_defs"); geo = pyproj.CRS.from_proj4(f"+proj=longlat +R={R} +no_defs")
n = msk.shape[0]; xx = (np.arange(n) + .5) * pix - half; Xg, Yg = np.meshgrid(xx, -xx)
lo_, la_ = pyproj.Transformer.from_crs(ps, geo, always_xy=True).transform(Xg, Yg); lo_ = np.mod(lo_, 360)
rr = ((90 - la_) * 10).astype(int).clip(0, 1799); cc = (lo_ * 10).astype(int).clip(0, 3599)
bp = cov[rr, cc]; good = (bp <= 40) & (la_ >= 30)
img = np.where(la_ < 30, np.nan, np.where(good, 1.0, np.where(np.isfinite(bp), 0.55, 0.0)))
c1.imshow(img, extent=[-half / 1e3, half / 1e3, -half / 1e3, half / 1e3], cmap=mpl.colors.ListedColormap(["#ffffff", "#e6e6e3", "#bdbdb8"]),
          norm=mpl.colors.BoundaryNorm([-.1, .3, .8, 1.1], 3), interpolation="nearest", zorder=0)
c1.contour(Xg / 1e3, Yg / 1e3, msk.astype(float), levels=[.5], colors=[ORANGE], linewidths=1.3, zorder=3)
t = pyproj.Transformer.from_crs(geo, ps, always_xy=True); fx, fy = t.transform(f.lon.values, f.lat.values)
sel = f.lat >= 30
c1.scatter(fx[sel] / 1e3, fy[sel] / 1e3, s=np.clip(f.n[sel].values, 2, 50) * .5 + 2, color=BLUE, edgecolors="white", linewidths=.3, zorder=4)
th = np.linspace(0, 2 * np.pi, 400)
for lc in (30, 45, 60, 75):
    rad = 2 * R * np.sin(np.radians(90 - lc) / 2) / 1e3
    c1.plot(rad * np.cos(th), rad * np.sin(th), color=MUTED, lw=.4, ls=(0, (2, 2)), zorder=2)
    c1.text(rad * np.cos(np.radians(-58)), rad * np.sin(np.radians(-58)), f"{lc}°", fontsize=5.8, color=MUTED, ha="left", va="top")
for lg_ in (0, 90, 180, 270):
    ang = np.radians(lg_ - 90); rad = 2 * R * np.sin(np.radians(60) / 2) / 1e3
    c1.plot([0, rad * np.cos(ang)], [0, rad * np.sin(ang)], color=MUTED, lw=.4, ls=(0, (2, 2)), zorder=2)
    c1.text(rad * (1.06 if lg_ in (0, 180) else 1.14) * np.cos(ang), rad * (1.06 if lg_ in (0, 180) else 1.14) * np.sin(ang), f"{lg_}°E", fontsize=5.8, color=MUTED, ha="center", va="center")
c1.set_xlim(-2900, 2900); c1.set_ylim(-2650, 2650); c1.set_aspect("equal"); c1.axis("off")
c1.text(-2850, 2600, "a", fontsize=10, fontweight="bold", va="top")
c1.legend(handles=[mpl.patches.Patch(color="#bdbdb8", label="imaged at ≤ 40 m/pixel"), mpl.patches.Patch(color="#e6e6e3", label="imaged at 40–150 m/pixel"),
                   mpl.lines.Line2D([], [], color=ORANGE, lw=1.3, label="northern smooth plains"),
                   mpl.lines.Line2D([], [], marker="o", ls="", color=BLUE, ms=4, label="hollow field")],
          loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=False, fontsize=6.2, handlelength=1.4, columnspacing=1.0)
# (b) latitude density in and out of the plains
D1 = json.load(open(f"{S}/latcorr_density_record.json")); D2 = json.load(open(f"{S}/latcorr_inc60_75_record.json"))
xb = {"0-30": 15, "30-40": 35, "40-60": 50, "60-80": 70}
def series(get, keys):
    out = [(xb[k],) + get(k) for k in keys]; return [o for o in out if o[1] is not None]
def g_all(zone, noprov=False):
    def g(k):
        if k == "0-30": v = D1[f"N_le40_{zone}_0-30"]
        elif k == "30-40": v = D1[f"N_le40_{zone}{'_noprov' if noprov else ''}"]["30-40"] if zone == "out" else {"d": None, "lo": None, "hi": None}
        else: v = D1[f"N_le40_{zone}_{k}{'_noprov' if noprov else ''}"]
        return (v["d"], v["lo"], v["hi"])
    return g
def g_inc(key):
    def g(k): v = D2[key][k]; return (v["d"], v["lo"], v["hi"])
    return g
keys = ["0-30", "30-40", "40-60", "60-80"]
for data, mk, col, mfc, off, lab in [
        (series(g_all("out"), keys), "o", BLUE, BLUE, -1.6, "outside plains"),
        (series(g_inc("le40_out"), keys), "s", BLUE, "white", 0.0, "outside, 60–75° incidence only"),
        (series(g_all("out", True), keys), "^", GREEN, "white", 1.6, "outside, 90–150°E province removed"),
        (series(g_all("in"), ["40-60", "60-80"]), "o", ORANGE, ORANGE, 0.8, "inside plains")]:
    x = np.array([d[0] for d in data]) + off; y = np.array([d[1] for d in data], float)
    lo = np.array([d[2] for d in data], float); hi = np.array([d[3] for d in data], float)
    c2.errorbar(x, y, yerr=[y - lo, hi - y], fmt=mk, ms=4.4, color=col, mfc=mfc, mew=1.1, lw=.9, capsize=2, zorder=5, label=lab)
for e in (30, 40, 60): c2.axvline(e, color=GRID, lw=.6, zorder=0)
c2.set_xlim(0, 80); c2.set_ylim(0, 152); c2.set_xticks([0, 30, 40, 60, 80])
c2.set_xlabel("latitude (°N)"); c2.set_ylabel("hollow fields per 10$^6$ km$^2$ imaged")
c2.legend(frameon=False, fontsize=6.1, loc="upper left", handletextpad=.3, borderpad=.2)
tidy(c2, "b", grid="y")

# (c) density on peaks, floors, walls by crater setting
PS = json.load(open(f"{PW}/peak_wall_by_setting.json"))
groups = [("isolated", "isolated"), ("inside_Orgel_basin_only", "basin-hosted"), ("inside_larger_Herrick_crater", "nested in a\nlarger crater")]
cls = [("peak", "central peak or ring", ORANGE, "o"), ("floor", "floor", GREEN, "s"), ("wall", "wall", BLUE, "D")]
for gi, (key, lab) in enumerate(groups):
    for ci, (k, nm_, col, mk) in enumerate(cls):
        v = PS[key][k]; xg = gi + (ci - 1) * 0.22
        d1.errorbar([xg], [v["d"]], yerr=[[v["d"] - v["ci"][0]], [v["ci"][1] - v["d"]]], fmt=mk, ms=5, color=col, mfc=col, mew=1, lw=1.1, capsize=2.2, zorder=5,
                    label=nm_ if gi == 0 else None)
    d1.text(gi, 0.62, f"{PS[key]['N']} craters", ha="center", va="bottom", fontsize=6, color=MUTED)
d1.set_yscale("log"); d1.set_ylim(.55, 40); d1.set_xlim(-.5, 2.5)
d1.set_yticks([1, 2, 5, 10, 20]); d1.set_yticklabels(["1", "2", "5", "10", "20"]); d1.minorticks_off()
d1.set_xticks([0, 1, 2]); d1.set_xticklabels([g[1] for g in groups], fontsize=6.8)
d1.set_ylabel("hollows per 1000 km$^2$ imaged")
d1.legend(frameon=False, fontsize=6.4, loc="upper right", handletextpad=.3, borderpad=.2, numpoints=1)
tidy(d1, "d", grid="y")
NUM["fig1d_terrain"] = PS

# (d) Caloris interior vs isolated craters
CR = json.load(open(f"{PW}/caloris_test_record.json"))
rows = CR["rates_Cmin0.2"]
edges = {"10-15": (10, 15), "15-20": (15, 20), "20-25": (20, 25), "25-30": (25, 30), "30-35": (30, 35), "35-40": (35, 40),
         "40-50": (40, 50), "50-70": (50, 70), "70-250": (70, 140)}
Dlo, Dhi = CR["D_for_excavation_depth"]["2.5"], CR["D_for_excavation_depth"]["3.5"]
d2.axvspan(Dlo, Dhi, color="#d8d3e8", alpha=.6, lw=0, zorder=0)
d2.text(np.sqrt(Dlo * Dhi) * 1.25, 1.27, "excavation reaches\nbase of plains", ha="center", va="bottom", fontsize=6.0, color="#4b3f72", linespacing=1.1)
for setk, col, mk, mfc, off, lab in [("iso_field", INK, "s", "white", 0.97, "isolated"), ("caloris_field", ORANGE, "o", ORANGE, 1.03, "inside Caloris")]:
    xs_, ys_, lo_, hi_ = [], [], [], []
    for rw in rows:
        k, n_ = map(int, rw[setk].split("/"))
        if n_ == 0: continue
        e0, e1 = edges[rw["D"]]; xg = np.sqrt(e0 * e1) * off
        l, h = cp(k, n_); xs_.append(xg); ys_.append(k / n_); lo_.append(l); hi_.append(h)
        if setk == "caloris_field": d2.text(xg, h + .02, f"{k}/{n_}", ha="center", va="bottom", fontsize=5.6, color=HOTTXT, rotation=90)
    ys_, lo_, hi_ = map(np.array, (ys_, lo_, hi_))
    d2.errorbar(xs_, ys_, yerr=[ys_ - lo_, hi_ - ys_], fmt=mk, ms=4.4, color=col, mfc=mfc, mew=1.1, lw=.9, capsize=1.8, zorder=5, label=lab)
d2.set_xscale("log"); d2.set_xlim(9, 150); d2.set_ylim(-.03, 1.47)
d2.set_xticks([10, 20, 30, 50, 100]); d2.set_xticklabels(["10", "20", "30", "50", "100"]); d2.minorticks_off()
d2.set_xlabel("crater diameter (km)"); d2.set_ylabel("fraction of craters hosting ≥ 3 hollows"); d2.set_yticks([0, .2, .4, .6, .8, 1.0])
d2.legend(frameon=False, fontsize=6.4, loc="upper left", handletextpad=.3, borderpad=.2, bbox_to_anchor=(0.0, 1.0))
tidy(d2, "c", grid="y")
save(fig, "fig1_exposure"); NUM["fig1c_caloris"] = dict(rows=rows, D_band=[Dlo, Dhi], counts=CR["caloris_counts"])

# ===================================================================================== Figure 3: depth
DR = json.load(open(f"{W}/records/depth_reanalysis_record.json"))
d = pd.read_csv(f"{W}/data/test1_rows_inc.csv"); d["ltan"] = np.log(np.tan(np.radians(d.inc)))
B = -9028.6679
fig, (e1, e2) = plt.subplots(1, 2, figsize=(5.9, 2.75))
bk = ["0-30", "30-45", "45-60", "60-75", "75-90"]
T_ = np.array([DR["bands"][k]["T_med"] for k in bk])
for col, mk, mfc, off, lab in [("hp", "o", "white", -0.8, "pixel scale corrected"), ("hc", "o", INK, 0.8, "pixel scale and incidence corrected")]:
    y = np.array([DR["bands"][k][col][0] for k in bk]); lo = np.array([DR["bands"][k][col][1] for k in bk]); hi = np.array([DR["bands"][k][col][2] for k in bk])
    e1.errorbar(T_ + off, y, yerr=[y - lo, hi - y], fmt=mk, ms=4.6, color=INK if col == "hc" else MUTED, mfc=mfc, mew=1.1, lw=.9, capsize=2, zorder=5, label=lab)
yhc = np.array([DR["bands"][k]["hc"][0] for k in bk])
hihc = np.array([max(DR["bands"][k]["hc"][2], DR["bands"][k]["hp"][2]) for k in bk])
for k, x_, y_ in zip(bk, T_, hihc):
    e1.text(x_, y_ * 1.12, k.replace("-", "–") + "°", fontsize=5.9, color=MUTED, ha="center", va="bottom", rotation=90)
Tg = np.linspace(236, 316, 100); T0 = DR["bands"]["0-30"]["T_med"]; h0 = DR["bands"]["0-30"]["hc"][0]
e1.plot(Tg, h0 * np.exp(B / 2 * (1 / Tg - 1 / T0)), ls=(0, (5, 2)), color=ORANGE, lw=1.5, zorder=3)
e1.text(258, 2.0, "diffusion-limited sulfur front\nanchored at 0–30°", fontsize=6.4, color=HOTTXT, ha="left", va="top", linespacing=1.15)
e1.set_yscale("log"); e1.set_ylim(.3, 120); e1.set_xlim(234, 318)
e1.set_yticks([.5, 1, 2, 5, 10, 20, 50, 100]); e1.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50", "100"]); e1.minorticks_off()
e1.set_xlabel("mean annual surface temperature (K)"); e1.set_ylabel("median hollow depth (m)")
e1.legend(frameon=False, fontsize=6.3, loc="lower right", handletextpad=.3, borderpad=.2)
tidy(e1, "a")

# (b) depth vs host crater diameter (smallest Herrick crater >= 5 km whose 1.1 R contains the hollow)
cr = pd.read_csv(f"{G}/data/test2_craters_h18.csv"); cr = cr[cr.diameter >= 5]
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
tree = cKDTree(xyz(cr.latitude.values, cr.lon360.values)); P = xyz(d.latitude.values, d.longitude.values % 360)
Dkm = cr.diameter.values; host = np.full(len(d), -1)
for i, cl in enumerate(tree.query_ball_point(P, r=0.36)):
    if not cl: continue
    cl = np.array(cl); ang = np.arccos(np.clip(P[i] @ tree.data[cl].T, -1, 1)) * R_M
    ok = cl[ang <= 1.1 * Dkm[cl] / 2]
    if len(ok): host[i] = ok[np.argmin(Dkm[ok])]
d["host"] = host; d["hostD"] = np.where(host >= 0, Dkm[np.maximum(host, 0)], np.nan)
dd = d[d.host >= 0].copy(); dd["logdep"] = np.log10(dd.depth); dd["logpx"] = np.log10(dd.px); dd["li"] = dd.ltan / np.log(10)
pc = dd.groupby("host").agg(n=("depth", "size"), D=("hostD", "first"), ld=("logdep", "median"), lp=("logpx", "median"), li=("li", "median"), al=("abslat", "median"))
pc = pc[pc.n >= 3].copy(); pc["logD"] = np.log10(pc.D)
mm = smf.ols("ld ~ logD + lp + li + al", data=pc).fit()
pc["ldc"] = pc.ld - mm.params["lp"] * (pc.lp - pc.lp.median()) - mm.params["li"] * (pc.li - pc.li.median()) - mm.params["al"] * (pc.al - pc.al.median())
e2.scatter(pc.D, 10 ** pc.ldc, s=np.clip(pc.n, 3, 80) * .6 + 4, facecolor="white", edgecolor=INK, lw=.8, zorder=4)
Dg = np.logspace(np.log10(pc.D.min() * .9), np.log10(pc.D.max() * 1.1), 50)
e2.plot(Dg, 10 ** (mm.params["Intercept"] + mm.params["logD"] * np.log10(Dg) + mm.params["lp"] * pc.lp.median() + mm.params["li"] * pc.li.median()
                   + mm.params["al"] * pc.al.median()), color=INK, lw=1.3, zorder=5)
e2.set_xscale("log"); e2.set_yscale("log"); e2.set_ylim(3, 150)
e2.set_xticks([5, 10, 30, 100, 300, 1000]); e2.set_xticklabels(["5", "10", "30", "100", "300", "1000"]); e2.minorticks_off()
e2.set_yticks([5, 10, 20, 50, 100]); e2.set_yticklabels(["5", "10", "20", "50", "100"])
e2.set_xlabel("host crater diameter (km)"); e2.set_ylabel("median depth per crater, corrected (m)")
sl_ = mm.params['logD']; sl_ = 0.0 if abs(sl_) < .005 else sl_
e2.text(.97, .95, f"slope {sl_:.2f} ± {mm.bse['logD']:.2f}".replace("-", "−") + f"\n{len(pc)} craters, {int(pc.n.sum())} hollows", transform=e2.transAxes,
        ha="right", va="top", fontsize=6.6, linespacing=1.2)
tidy(e2, "b")
fig.subplots_adjust(left=.09, right=.985, top=.9, bottom=.17, wspace=.34); save(fig, "fig4_depth")
# variance among host craters (one-way ANOVA on log depth, craters with >= 3 hollows)
from scipy.stats import f_oneway
cc_ = smf.ols("logdep ~ logpx + li", data=dd).fit(); dd["res"] = cc_.resid
grp = [g.res.values for _, g in dd[dd.host.isin(pc.index)].groupby("host")]
F, pF = f_oneway(*grp)
NUM["fig3b"] = dict(n_crater=int(len(pc)), n_hollows=int(pc.n.sum()), slope=float(mm.params["logD"]), se=float(mm.bse["logD"]), p=float(mm.pvalues["logD"]),
                    maxD=float(pc.D.max()), minD=float(pc.D.min()), anova_F=float(F), anova_p=float(pF))
json.dump(NUM, open(f"{OUT}/../figure_numbers_pr.json", "w"), indent=1)
print(json.dumps({k: NUM[k] for k in ("fig2b_odds", "fig3b")}, indent=1))
