# g4 copy of code/rerun/build_hollow_products.py, blocks "Build coverage grid", "Build completeness grid" and "Attach coverage / completeness" (6 Nov 2026). Changed: detection-model parameters (theta read from a JSON file instead of the in-script fit), input catalog (hollows / fields of the selected catalog instead of the freshly built v2 files), and input/output paths (/Users/... -> ~/mnt/hollows-emissivity; outputs in g4_completeness/). Code otherwise verbatim.
"""Usage: python3 build_completeness.py <tag> [coverage|grid|attach|all]
  tag v2chk : out_v2 theta (analysis_workshop/completeness_sfd/theta_out_v2.json) on the v2 catalog -> reproduction check
  tag v3r   : g3 refit theta for v3 (g3_catalog/sfd/out_v3_refit/theta_out_v2.json) on hollows_canonical_v3.csv
  tag v3br  : g3 refit theta for v3b on hollows_canonical_v3b.csv
  tag v3chk : attach only, stored out_v2 grid on the v3 catalog (check of the attach code against hollows_canonical_v3.csv)
coverage: rebuilds coverage_10ppd.npz from the MDIS index (same code) into _scratch/ (tag-independent).
grid    : completeness grid -> completeness_10ppd_<v3|v3b|v2chk>.npz in this folder.
attach  : per-hollow and per-field cov_* and C_* columns -> inputs_<tag>/hollows_canonical.csv, fields.csv."""
import numpy as np
import pandas as pd
import json
import os
import csv
import io
import re
import sys
from scipy.special import erf

RM_ = 2439.7
KM_PER_DEG = RM_ * np.pi / 180
PJ = os.path.expanduser("~/mnt/hollows-emissivity/")
HN = PJ + "hollow-networks/"
G4 = os.path.dirname(os.path.abspath(__file__)) + "/"
TAG = sys.argv[1]; STEP = sys.argv[2] if len(sys.argv) > 2 else "all"
CFG = {"v2chk": dict(theta=PJ + "manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json",
                     hol=PJ + "code/rerun/out_v2/hollows_canonical.csv", fld=PJ + "code/rerun/out_v2/fields.csv", grid="_scratch/completeness_10ppd_v2chk.npz"),
       "v3r": dict(theta=PJ + "code/rerun/out_v3/rerun/g3_catalog/sfd/out_v3_refit/theta_out_v2.json",
                   hol=PJ + "code/rerun/out_v3/hollows_canonical_v3.csv", fld=PJ + "code/rerun/out_v3/rerun/g2_insolation/inputs/v3/fields.csv", grid="completeness_10ppd_v3.npz"),
       "v3br": dict(theta=PJ + "code/rerun/out_v3/rerun/g3_catalog/sfd/out_v3b_refit/theta_out_v2.json",
                    hol=PJ + "code/rerun/out_v3/hollows_canonical_v3b.csv", fld=PJ + "code/rerun/out_v3/rerun/g2_insolation/inputs/v3b/fields.csv", grid="completeness_10ppd_v3b.npz"),
       # attach-only check: the stored v2 grid on the v3 catalog must give back the C columns already in hollows_canonical_v3.csv
       "v3chk": dict(theta=PJ + "manuscripts/PlanetaryResearch_hollows/analysis_workshop/completeness_sfd/theta_out_v2.json",
                     hol=PJ + "code/rerun/out_v3/hollows_canonical_v3.csv", fld=PJ + "code/rerun/out_v3/rerun/g2_insolation/inputs/v3/fields.csv", grid=PJ + "code/rerun/out_v2/completeness_10ppd.npz")}[TAG]
GRIDF = CFG["grid"] if CFG["grid"].startswith("/") else G4 + CFG["grid"]
COVF = G4 + "_scratch/coverage_10ppd_rebuilt.npz"
th = np.array(json.load(open(CFG["theta"]))["theta_raw"])
PPD = 10
NY, NX = 180 * PPD, 360 * PPD
BANDS = [(0, 10), (10, 20), (20, 35), (35, 60), (60, 100), (100, 150.01)]

if STEP in ("coverage", "all"):
    # ---------------- Build coverage grid (verbatim; the index is already on disk, so no download)
    WANT = ['PRODUCT_ID', 'INSTRUMENT_ID', 'START_TIME', 'FILTER_NUMBER', 'CENTER_LATITUDE',
            'CENTER_LONGITUDE', 'HORIZONTAL_PIXEL_SCALE', 'RETICLE_POINT_LATITUDE_1',
            'RETICLE_POINT_LATITUDE_2', 'RETICLE_POINT_LATITUDE_3', 'RETICLE_POINT_LATITUDE_4',
            'RETICLE_POINT_LONGITUDE_1', 'RETICLE_POINT_LONGITUDE_2', 'RETICLE_POINT_LONGITUDE_3',
            'RETICLE_POINT_LONGITUDE_4', 'INCIDENCE_ANGLE', 'EMISSION_ANGLE', 'PHASE_ANGLE',
            'EXPOSURE_DURATION', 'MISSION_PHASE_NAME', 'DATA_QUALITY_ID']
    idx_path = HN + "data/mdis_cdr_index_reduced.csv"
    IDX = pd.read_csv(idx_path, low_memory=False)
    num = [c for c in IDX.columns if c not in ('PRODUCT_ID', 'INSTRUMENT_ID', 'START_TIME',
                                                 'MISSION_PHASE_NAME', 'DATA_QUALITY_ID')]
    for c in num:
        IDX[c] = pd.to_numeric(IDX[c], errors='coerce')
    IDX['cam'] = np.where(IDX.INSTRUMENT_ID.str.contains('NAC'), 'N', 'W')
    IDX['pid'] = IDX.PRODUCT_ID.str.split('_').str[0].str.strip()
    ok = (IDX.HORIZONTAL_PIXEL_SCALE.gt(0) & IDX.CENTER_LATITUDE.between(-90, 90) &
          IDX.INCIDENCE_ANGLE.between(0, 180))
    IDX = IDX[ok & IDX.MISSION_PHASE_NAME.str.contains('ORBIT', na=False)].drop_duplicates('pid').reset_index(drop=True)
    S = IDX[(IDX.cam == 'N') & (IDX.HORIZONTAL_PIXEL_SCALE <= 150) & (IDX.INCIDENCE_ANGLE < 85)].reset_index(drop=True)
    clat = S.CENTER_LATITUDE.values
    clon = S.CENTER_LONGITUDE.values % 360
    px = S.HORIZONTAL_PIXEL_SCALE.values
    inc = S.INCIDENCE_ANGLE.values
    emi = S.EMISSION_ANGLE.values
    grow = (1 / np.cos(np.radians(np.clip(emi, 0, 80)))) ** 0.8
    half_km = 1024 * px / 1000 / 2 * grow
    half_lat = half_km / KM_PER_DEG
    half_lon = half_lat / np.maximum(np.cos(np.radians(clat)), 0.05)
    band_of = np.digitize(px, [b[1] for b in BANDS[:-1]])
    nimg = np.zeros((NY, NX), np.int32)
    bestpx = np.full((NY, NX), np.inf, np.float32)
    nband = np.zeros((len(BANDS), NY, NX), np.int16)
    i0 = np.clip(((90 - (clat + half_lat)) * PPD).astype(int), 0, NY - 1)
    i1 = np.clip(((90 - (clat - half_lat)) * PPD).astype(int), 0, NY - 1)
    j0 = ((clon - half_lon) * PPD).astype(int)
    j1 = ((clon + half_lon) * PPD).astype(int)
    for k in np.argsort(-px):
        a, b = i0[k], i1[k] + 1
        js = np.arange(j0[k], j1[k] + 1) % NX
        nimg[a:b][:, js] += 1
        nband[band_of[k], a:b][:, js] += 1
        sub = bestpx[a:b][:, js]
        np.minimum(sub, px[k], out=sub)
        bestpx[a:b][:, js] = sub
    bestpx[~np.isfinite(bestpx)] = np.nan
    np.savez_compressed(COVF, nimg=nimg, bestpx=bestpx, nband=nband, bands=np.array(BANDS), ppd=PPD,
        note="Searched set: MDIS NAC CDR, orbital, px<=150 m, inc<85. Footprint corrected for emission.")
    o = np.load(PJ + "code/rerun/out_v2/coverage_10ppd.npz")
    chk = {k: bool(np.array_equal(o[k], np.load(COVF)[k], equal_nan=(k == "bestpx"))) for k in ("nimg", "bestpx", "nband")}
    print("images", len(S), "coverage identical to out_v2:", chk)
    json.dump(dict(images=int(len(S)), identical_to_out_v2=chk), open(G4 + "_scratch/coverage_check.json", "w"))

def Phi(z):
    return 0.5 * (1 + erf(z / np.sqrt(2)))

def pdet(x, inc_v, th):
    pmax, lx50, ls1, lxhi, ls2, binc = th
    p = pmax * Phi((np.log(x) - lx50) / np.exp(ls1)) * (1 - Phi((np.log(x) - lxhi) / np.exp(ls2)))
    return np.clip(p * np.exp(binc * (inc_v - 60) / 25), 1e-6, 0.999)

Dref = np.array([0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
if STEP in ("grid", "all"):
    # ---------------- Build completeness grid (verbatim)
    nband = np.load(COVF)["nband"]
    bpx = np.array([np.sqrt(max(a, 4) * b) for a, b in BANDS])
    C = np.ones((len(Dref), NY, NX), np.float32)
    for bi in range(len(BANDS)):
        for di, D_r in enumerate(Dref):
            pb = pdet(D_r * 1000 / bpx[bi], 60.0, th)
            C[di] *= (1 - pb) ** nband[bi].astype(np.float32)
    C = 1 - C
    assert TAG != "v3chk"
    np.savez_compressed(GRIDF, C=C, Dref=Dref, ppd=PPD,
        note=f"Catalog completeness C(D,cell). Axis 0 = Dref (km). 10 ppd grid. g4 ({TAG}): theta from {CFG['theta'].replace(PJ, '')}")
    print("wrote", CFG["grid"], "theta", np.round(th, 3).tolist())

if STEP in ("attach", "all"):
    # ---------------- Attach coverage / completeness (verbatim, applied to the selected catalog's hollows and fields)
    HOLa = pd.read_csv(CFG["hol"])
    FLD = pd.read_csv(CFG["fld"])
    Z = np.load(GRIDF)
    Cg = Z['C']
    cv = np.load(COVF)
    def samp(grid, lat, lon):
        i = np.clip(((90 - lat) * PPD).astype(int), 0, NY - 1)
        j = ((lon % 360) * PPD).astype(int) % NX
        return grid[i, j]
    for df in (HOLa, FLD):
        df['cov_nframes'] = samp(cv['nimg'], df.lat.values, df.lon.values)
        df['cov_bestpx_m'] = samp(cv['bestpx'], df.lat.values, df.lon.values)
        for D_r in [0.5, 1.0, 2.0, 5.0]:
            df[f'C_D{D_r:g}km'] = samp(Cg[list(Dref).index(D_r)], df.lat.values, df.lon.values)
    di = np.clip(np.searchsorted(Dref, HOLa.D_eq_km.values), 1, len(Dref) - 1)
    i_idx = np.clip(((90 - HOLa.lat.values) * PPD).astype(int), 0, NY - 1)
    j_idx = ((HOLa.lon.values % 360) * PPD).astype(int) % NX
    lo_ = Cg[di - 1, i_idx, j_idx]
    hi_ = Cg[di, i_idx, j_idx]
    f = ((np.log(HOLa.D_eq_km.values) - np.log(Dref[di - 1])) /
         (np.log(Dref[di]) - np.log(Dref[di - 1])))
    HOLa['C_own'] = np.clip(lo_ + (hi_ - lo_) * np.clip(f, 0, 1), 1e-4, 1)
    D = G4 + f"_scratch/attach_{TAG}"; os.makedirs(D, exist_ok=True)
    HOLa.to_csv(D + "/hollows_canonical.csv", index=False)
    FLD.to_csv(D + "/fields.csv", index=False)
    print("attached", len(HOLa), "hollows", len(FLD), "fields ->", D)
