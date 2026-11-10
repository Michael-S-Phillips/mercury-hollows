# g4 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/crater_completeness.py (6 Nov 2026). Changed: input/output paths by S2_CAT, added exponents/measures 1.7 and 1.8 (see make_copies_s2.py).
#!/usr/bin/env python3
"""Completeness at crater centers, old and size-distribution-weighted (4 Nov 2026).
Images and footprints exactly as code/rerun/build_hollow_products.py (searched set: MDIS NAC CDR, orbital, px <= 150 m,
inc < 85; emission-corrected square footprints on the 10 ppd grid); a crater is covered by an image when its 10 ppd cell lies
in the image's cell range, as in the grid. Measures per crater:
  C1_old : 1 km hollows, each image at 60 deg incidence and its band's representative pixel scale (the paper's C1; checked).
  C1_act : 1 km hollows, each image at its own incidence and pixel scale.
  C(D)   : on a D grid 0.03-30 km, own incidence and pixel scale.
  E_a    : expected number of hollows detected per hollow >= 1 km present, for a power-law size distribution with cumulative
           slope a (N(>D) ~ D^-a; a = 1.6 fitted in fit_sfd.py (1.55-1.63 for cutoffs 0.2-0.5 km), 1.3 and 1.9 as checks): E_a = a * int D^-a C(D) dlnD.
           Cutoff-free: C(D) -> 0 at small D, so the integral converges; a constant factor drops out of the regressions.
  F_0.3  : fraction of hollows >= 0.3 km detected (same distribution, a = 1.6), a probability-like version with a cutoff.
Writes crater_completeness.csv (id + measures)."""
import os, json, numpy as np, pandas as pd
import os as _o, sys as _s; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__))); from s2cfg import *   # g4: inputs by S2_CAT
from scipy.special import erf
HE = os.path.expanduser("~/mnt/hollows-emissivity"); HERE = os.path.dirname(os.path.abspath(__file__)); RM_ = 2439.7; KM_PER_DEG = RM_ * np.pi / 180
PF = json.load(open(THETA))   # out_v2 detection model (refit_detection.py; matches the paper, S2)
def pdet(x, inc):
    Phi = lambda z: 0.5 * (1 + erf(z / np.sqrt(2)))
    p = PF["p_max"] * Phi((np.log(x) - np.log(PF["x50_px"])) / PF["sigma_lo"]) * (1 - Phi((np.log(x) - np.log(PF["x_hi_px"])) / PF["sigma_hi"]))
    return np.clip(p * np.exp(PF["b_inc"] * (inc - 60) / 25), 1e-6, 0.999)
IDX = pd.read_csv(f"{HE}/hollow-networks/data/mdis_cdr_index_reduced.csv", low_memory=False, usecols=["PRODUCT_ID", "INSTRUMENT_ID", "MISSION_PHASE_NAME", "CENTER_LATITUDE", "CENTER_LONGITUDE", "HORIZONTAL_PIXEL_SCALE", "INCIDENCE_ANGLE", "EMISSION_ANGLE"])
for c in [c for c in IDX.columns if c not in ('PRODUCT_ID', 'INSTRUMENT_ID', 'START_TIME', 'MISSION_PHASE_NAME', 'DATA_QUALITY_ID')]:
    IDX[c] = pd.to_numeric(IDX[c], errors='coerce')
IDX['cam'] = np.where(IDX.INSTRUMENT_ID.str.contains('NAC'), 'N', 'W'); IDX['pid'] = IDX.PRODUCT_ID.str.split('_').str[0].str.strip()
ok = IDX.HORIZONTAL_PIXEL_SCALE.gt(0) & IDX.CENTER_LATITUDE.between(-90, 90) & IDX.INCIDENCE_ANGLE.between(0, 180)
IDX = IDX[ok & IDX.MISSION_PHASE_NAME.str.contains('ORBIT', na=False)].drop_duplicates('pid').reset_index(drop=True)
S = IDX[(IDX.cam == 'N') & (IDX.HORIZONTAL_PIXEL_SCALE <= 150) & (IDX.INCIDENCE_ANGLE < 85)].reset_index(drop=True)
PPD = 10; NY, NX = 180 * PPD, 360 * PPD
clat = S.CENTER_LATITUDE.values; clon = S.CENTER_LONGITUDE.values % 360; px = S.HORIZONTAL_PIXEL_SCALE.values; inc = S.INCIDENCE_ANGLE.values; emi = S.EMISSION_ANGLE.values
grow = (1 / np.cos(np.radians(np.clip(emi, 0, 80)))) ** 0.8; half_km = 1024 * px / 1000 / 2 * grow
half_lat = half_km / KM_PER_DEG; half_lon = half_lat / np.maximum(np.cos(np.radians(clat)), 0.05)
BANDS = [(0, 10), (10, 20), (20, 35), (35, 60), (60, 100), (100, 150.01)]; band_of = np.digitize(px, [b[1] for b in BANDS[:-1]])
bpx = np.array([np.sqrt(max(a, 4) * b) for a, b in BANDS])
i0 = np.clip(((90 - (clat + half_lat)) * PPD).astype(int), 0, NY - 1); i1 = np.clip(((90 - (clat - half_lat)) * PPD).astype(int), 0, NY - 1)
j0 = ((clon - half_lon) * PPD).astype(int); j1 = ((clon + half_lon) * PPD).astype(int)
d = pd.read_csv(FEAT, dtype={"id": str}, usecols=["id", "lat", "lon", "C1"])
rr = ((90 - d.lat.values) * 10).astype(int).clip(0, 1799); cc = ((d.lon.values % 360) * 10).astype(int).clip(0, 3599)
# crater-image pairs: craters whose 10 ppd cell lies in the image's cell range (rows i0..i1, columns j0..j1 with wrap)
order = np.argsort(rr, kind="stable"); rs = rr[order]; PI, PK = [], []
for k in range(len(S)):
    a_, b_ = np.searchsorted(rs, i0[k], "left"), np.searchsorted(rs, i1[k], "right")
    if b_ <= a_: continue
    cand_ = order[a_:b_]; sel = cand_[((cc[cand_] - j0[k]) % NX) <= (j1[k] - j0[k])]
    if len(sel): PI.append(sel); PK.append(np.full(len(sel), k))
PI = np.concatenate(PI); PK = np.concatenate(PK); so = np.argsort(PI, kind="stable"); PI, PK = PI[so], PK[so]
starts = np.searchsorted(PI, np.arange(len(d) + 1)); print("images", len(S), "crater-image pairs", len(PI), flush=True)
Dg = np.exp(np.linspace(np.log(0.03), np.log(30), 61)); dl = np.log(Dg[1] / Dg[0])
out = np.zeros((len(d), 2)); CD = np.zeros((len(d), len(Dg)), np.float32); nimg = np.zeros(len(d), int)
p_old = pdet(1000 / bpx, 60.0)
for n in range(len(d)):
    k = PK[starts[n]:starts[n + 1]]
    nimg[n] = len(k)
    if not len(k): continue
    out[n, 0] = 1 - np.prod(1 - p_old[band_of[k]])
    out[n, 1] = 1 - np.prod(1 - pdet(1000 / px[k], inc[k]))
    L = np.log1p(-pdet(1000 * Dg[:, None] / px[k][None, :], inc[k][None, :])).sum(1); CD[n] = 1 - np.exp(L)
d["nimg"] = nimg; d["C1_old"] = out[:, 0]; d["C1_act"] = out[:, 1]
for a in (1.3, 1.6, 1.7, 1.8, 1.9):
    d[f"E_{a}"] = a * (CD * Dg[None, :] ** (-a)).sum(1) * dl
g = Dg >= 0.3; w = Dg[g] ** (-1.6); d["F_0.3"] = (CD[:, g] * w).sum(1) / w.sum()
for D0 in (0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0): d[f"C_{D0}_act"] = [np.interp(np.log(D0), np.log(Dg), r) for r in CD]
print("check C1_old vs paper C1: max abs diff", float(np.abs(d.C1_old - d.C1).max()), " median", float(np.median(np.abs(d.C1_old - d.C1))))
print(d[["C1", "C1_old", "C1_act", "E_1.6", "F_0.3", "nimg"]].describe().round(3).to_string())
print("Spearman with C1:", d[["C1", "C1_act", "E_1.3", "E_1.6", "E_1.9", "F_0.3", "C_0.3_act", "C_2.0_act"]].corr("spearman").iloc[0].round(3).to_dict())
d.drop(columns=["lat", "lon"]).to_csv(CCOMP, index=False)
