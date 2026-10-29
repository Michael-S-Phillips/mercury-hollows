import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
import numpy as np
import pandas as pd
import json
import os
import csv
import io
import re
import glob
import time
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.optimize import minimize
from scipy.special import erf

RM_ = 2439.7
KM_PER_DEG = RM_ * np.pi / 180
import sys
HN = _REPO + "/hollow-networks/"
USE_CORRECTED = sys.argv[1] == "corrected"
OUTDIR = sys.argv[2].rstrip("/") + "/"
GPKG = _REPO + "/GIS/vector/bickel2024_hollows_corrected.gpkg"
os.makedirs(OUTDIR, exist_ok=True)
PJ = _REPO + "/"

os.makedirs(HN + "data", exist_ok=True)
os.makedirs(HN + "results", exist_ok=True)

# Parse Bickel catalog
def sph_poly_area_km2(lon, lat):
    la0 = np.mean(lat); lo0 = np.mean(lon)
    x = np.radians(((lon - lo0 + 180) % 360) - 180) * np.cos(np.radians(la0)) * RM_
    y = np.radians(lat - la0) * RM_
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

def parse(fpath):
    raw = open(fpath, 'rb').read().decode('latin-1')
    rows = list(csv.reader(io.StringIO(raw), delimiter=';'))
    hdr = [c.strip() for c in rows[0]]
    ix = {c: i for i, c in enumerate(hdr)}
    out = []
    for r in rows[1:]:
        if len(r) < len(hdr): continue
        try:
            g = r[ix['Polgyon geometry (° 180d)']]
            pts = np.array([[float(a), float(b)] for a, b in
                            re.findall(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', g)])
            if len(pts) < 3: continue
            A = sph_poly_area_km2(pts[:, 0], pts[:, 1])
            xs = np.radians(((pts[:, 0] - pts[:, 0].mean() + 180) % 360) - 180) * np.cos(np.radians(pts[:, 1].mean())) * RM_
            ys = np.radians(pts[:, 1] - pts[:, 1].mean()) * RM_
            ext = max(np.ptp(xs), np.ptp(ys))
            out.append(dict(
                image=r[ix['MDIS ID']], lat=float(r[ix['Detection center LAT']]),
                lon=float(r[ix['Detection center LON 180d']]) % 360.0,
                px_m=float(r[ix['Pixel scale (m)']]), inc=float(r[ix['Incidence (°)']]),
                emi=float(r[ix['Emission (°)']]), phase=float(r[ix['Phase (°)']]),
                conf=float(r[ix['Confidence']]), stage=int(float(r[ix['Stage']])),
                area_km2=A, D_eq_km=2 * np.sqrt(A / np.pi), max_ext_km=ext, nvert=len(pts),
                bx=float(r[ix['Image x coord max (pix)']]) - float(r[ix['Image x coord min (pix)']]),
                by=float(r[ix['Image y coord max (pix)']]) - float(r[ix['Image y coord min (pix)']]),
                thumb=r[ix['Thumbnail tag']].strip())
            )
        except (ValueError, IndexError):
            continue
    df = pd.DataFrame(out)
    if USE_CORRECTED:
        import sqlite3
        con = sqlite3.connect(GPKG)
        C = pd.read_sql("select thumbnail, lat as lat_c, lon180 as lon_c, shift_km from hollows_points_corrected", con)
        C['thumbnail'] = C.thumbnail.str.strip()
        assert C.thumbnail.is_unique
        df = df.merge(C, left_on='thumb', right_on='thumbnail', how='left')
        assert df.lat_c.notna().all(), f"{df.lat_c.isna().sum()} detections without corrected coordinates"
        df['lat_release'] = df.lat; df['lon_release'] = df.lon
        df['lat'] = df.lat_c; df['lon'] = df.lon_c % 360.0
        df = df.drop(columns=['thumbnail', 'lat_c', 'lon_c'])
    return df

D = PJ + "data/bickel-2024/"
D7 = parse(glob.glob(D + "Bickel_et_al_2024_MDIS_hollows_07_v1.0.csv")[0])
qc = (D7.conf >= 0.75) & (D7.inc < 85)
DT = D7[qc].reset_index(drop=True)

# Deduplicate
def xyz(lat, lon):
    return np.c_[np.cos(np.radians(lat)) * np.cos(np.radians(lon)),
                 np.cos(np.radians(lat)) * np.sin(np.radians(lon)),
                 np.sin(np.radians(lat))]

XYZ = xyz(DT.lat.values, DT.lon.values)
tree = cKDTree(XYZ)
R = (DT.D_eq_km.values / 2)
pairs = tree.query_pairs(2 * np.sin(3.0 / (2 * RM_)))
same = []
img = DT.image.values
for i, j in pairs:
    d = 2 * RM_ * np.arcsin(np.linalg.norm(XYZ[i] - XYZ[j]) / 2)
    if d < 0.5 * (R[i] + R[j]) and img[i] != img[j]:
        same.append((i, j))

if same:
    a = np.array(same)
    g = coo_matrix((np.ones(len(a)), (a[:, 0], a[:, 1])), shape=(len(DT), len(DT)))
    ncomp, lab = connected_components(g, directed=False)
    DT['hollow_id'] = lab
else:
    DT['hollow_id'] = np.arange(len(DT))

idx = DT.sort_values('px_m').groupby('hollow_id').head(1).index
HOL = DT.loc[idx].copy().sort_values('hollow_id').reset_index(drop=True)
grp = DT.groupby('hollow_id')
HOL['n_det'] = grp.size().reindex(HOL.hollow_id).values
HOL['D_eq_sd_frac'] = (grp.D_eq_km.std() / grp.D_eq_km.mean()).reindex(HOL.hollow_id).fillna(0).values
HOL['px_best'] = HOL.px_m
HOL['is_crater_trace'] = HOL.D_eq_km > 10.0

DT.to_csv(OUTDIR + "detections_qc.csv", index=False)
HOL.to_csv(OUTDIR + "hollows_canonical.csv", index=False)

# Build coverage grid
WANT = ['PRODUCT_ID', 'INSTRUMENT_ID', 'START_TIME', 'FILTER_NUMBER', 'CENTER_LATITUDE',
        'CENTER_LONGITUDE', 'HORIZONTAL_PIXEL_SCALE', 'RETICLE_POINT_LATITUDE_1',
        'RETICLE_POINT_LATITUDE_2', 'RETICLE_POINT_LATITUDE_3', 'RETICLE_POINT_LATITUDE_4',
        'RETICLE_POINT_LONGITUDE_1', 'RETICLE_POINT_LONGITUDE_2', 'RETICLE_POINT_LONGITUDE_3',
        'RETICLE_POINT_LONGITUDE_4', 'INCIDENCE_ANGLE', 'EMISSION_ANGLE', 'PHASE_ANGLE',
        'EXPOSURE_DURATION', 'MISSION_PHASE_NAME', 'DATA_QUALITY_ID']

import urllib.request

lbl_path = HN + "data/CUMINDEX.LBL"
if not os.path.exists(lbl_path):
    ROOT = "https://pds-imaging.jpl.nasa.gov/data/messenger/"
    with urllib.request.urlopen(urllib.request.Request(ROOT + "msgrmds_2001/INDEX/CUMINDEX.LBL"), timeout=90) as r:
        lbl = r.read().decode('utf8', 'ignore')
    open(lbl_path, "w").write(lbl)
else:
    lbl = open(lbl_path).read()

blocks = re.findall(r'OBJECT\s*=\s*COLUMN(.*?)END_OBJECT\s*=\s*COLUMN', lbl, re.S)
COLS = []
for b in blocks:
    nm = re.search(r'\bNAME\s*=\s*"?([A-Z_0-9]+)"?', b).group(1)
    cn = re.search(r'COLUMN_NUMBER\s*=\s*(\d+)', b)
    COLS.append((int(cn.group(1)) if cn else None, nm))
COLS.sort()
pos = {nm: cn - 1 for cn, nm in COLS}
keep = [pos[w] for w in WANT]

idx_path = HN + "data/mdis_cdr_index_reduced.csv"
if not os.path.exists(idx_path):
    URL = "https://pds-imaging.jpl.nasa.gov/data/messenger/msgrmds_2001/INDEX/CUMINDEX.TAB"
    out = open(idx_path, "w", newline="")
    w = csv.writer(out)
    w.writerow(WANT)
    n = 0
    buf = b""
    with urllib.request.urlopen(urllib.request.Request(URL), timeout=300) as r:
        while True:
            ch = r.read(1 << 22)
            if not ch: break
            buf += ch
            lines = buf.split(b"\n")
            buf = lines.pop()
            for ln in lines:
                s = ln.decode('latin-1').strip()
                if not s: continue
                p = next(csv.reader([s]))
                if len(p) <= max(keep): continue
                row = [p[i].strip() for i in keep]
                w.writerow(row)
                n += 1
    out.close()

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

S_full = IDX[(IDX.cam == 'N') & (IDX.HORIZONTAL_PIXEL_SCALE <= 150) & (IDX.INCIDENCE_ANGLE < 85)].reset_index(drop=True)
# INCIDENCE-CONTROLLED NULL: only frames with 60 <= incidence <= 75 deg count as searching.
INC_LO, INC_HI = 60.0, 75.0
S = S_full[S_full.INCIDENCE_ANGLE.between(INC_LO, INC_HI)].reset_index(drop=True)
print("searched frames: full=%d  inc60_75=%d" % (len(S_full), len(S)))

PPD = 10
NY, NX = 180 * PPD, 360 * PPD
latc = 90 - (np.arange(NY) + 0.5) / PPD
W = np.cos(np.radians(latc))[:, None] * np.ones((1, NX))

clat = S.CENTER_LATITUDE.values
clon = S.CENTER_LONGITUDE.values % 360
px = S.HORIZONTAL_PIXEL_SCALE.values
inc = S.INCIDENCE_ANGLE.values
emi = S.EMISSION_ANGLE.values

grow = (1 / np.cos(np.radians(np.clip(emi, 0, 80)))) ** 0.8
half_km = 1024 * px / 1000 / 2 * grow
half_lat = half_km / KM_PER_DEG
half_lon = half_lat / np.maximum(np.cos(np.radians(clat)), 0.05)

BANDS = [(0, 10), (10, 20), (20, 35), (35, 60), (60, 100), (100, 150.01)]
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
cov = nimg > 0

np.savez_compressed(OUTDIR + "coverage_10ppd.npz", nimg=nimg, bestpx=bestpx, nband=nband,
    bands=np.array(BANDS), ppd=PPD,
    note="Searched set: MDIS NAC CDR, orbital, px<=150 m, inc<85. Footprint corrected for emission.")

# Build completeness pairs
HOLa = pd.read_csv(OUTDIR + "hollows_canonical.csv")
H = HOLa[~HOLa.is_crater_trace].copy()
hx = xyz(H.lat.values, H.lon.values)
Hn = H.reset_index(drop=True)

det_by_hollow = DT.groupby('hollow_id').image.apply(lambda s: set(x.strip() for x in s)).to_dict()

TI = cKDTree(xyz(clat, clon))
rmax = 2 * np.sin((half_km.max() * 1.5) / (2 * RM_))
hx2 = xyz(Hn.lat.values, Hn.lon.values)
nb = TI.query_ball_point(hx2, rmax)
pidS = S.pid.values
pid_ix = {p: i for i, p in enumerate(pidS)}

rows = []
for i in range(len(Hn)):
    hl, ho = Hn.lat.iloc[i], Hn.lon.iloc[i]
    D_h = Hn.D_eq_km.iloc[i]
    dets = det_by_hollow.get(Hn.hollow_id.iloc[i], set())
    seen = set()
    for k in nb[i]:
        if (abs(hl - clat[k]) <= half_lat[k] and
                abs(((ho - clon[k] + 180) % 360) - 180) <= half_lon[k]):
            rows.append((i, D_h, px[k], inc[k], D_h * 1000 / px[k], pidS[k] in dets))
            seen.add(k)
    for p_ in dets:
        k = pid_ix.get(p_)
        if k is not None and k not in seen:
            rows.append((i, D_h, px[k], inc[k], D_h * 1000 / px[k], True))

CV = pd.DataFrame(rows, columns=['h', 'D_km', 'px_m', 'inc', 'D_px', 'detected']).sort_values('h').reset_index(drop=True)
CV.to_csv(OUTDIR + "completeness_pairs.csv", index=False)

# Fit per-frame detection function
def Phi(z):
    return 0.5 * (1 + erf(z / np.sqrt(2)))

def pdet(x, inc_v, th):
    pmax, lx50, ls1, lxhi, ls2, binc = th
    p = pmax * Phi((np.log(x) - lx50) / np.exp(ls1)) * (1 - Phi((np.log(x) - lxhi) / np.exp(ls2)))
    return np.clip(p * np.exp(binc * (inc_v - 60) / 25), 1e-6, 0.999)

# INCIDENCE-CONTROLLED NULL: the truncated binomial conditions on >=1 detection per hollow.
# Under the restricted frame set many hollows have zero detections in 60-75 deg frames
# (their detecting frames fell outside the band); those violate the truncation and must be
# excluded from the fit. Keep the full restricted pairs on disk (saved above) but fit theta
# only on hollows with >=1 detection in the restricted frames.
_det_any = CV.groupby('h').detected.transform('any').values
CV_fit = CV[_det_any].copy().sort_values('h').reset_index(drop=True)
print("fit conditioning: hollows in pairs=%d  with>=1 restricted detection=%d"
      % (CV.h.nunique(), CV_fit.h.nunique()))
CV = CV_fit
hid = CV.h.values
X = CV.D_px.values
I = CV.inc.values
Y = CV.detected.values.astype(float)
uh, starts = np.unique(hid, return_index=True)
seg = np.r_[starts, len(hid)]
nfr = np.diff(seg)

def nll_trunc(th):
    if not (0.005 < th[0] < 0.95) or th[1] > th[3]: return 1e18
    p = pdet(X, I, th)
    lp = Y * np.log(p) + (1 - Y) * np.log(1 - p)
    l1p = np.log(1 - p)
    S1 = np.add.reduceat(lp, starts)
    S0 = np.add.reduceat(l1p, starts)
    return -np.sum(S1 - np.log(1 - np.exp(S0) + 1e-300))

best = None
for p0 in [0.05, 0.1, 0.2]:
    for x0 in [12, 18, 25]:
        r = minimize(nll_trunc, [p0, np.log(x0), np.log(0.5), np.log(170), np.log(0.35), 0.0],
                     method='Nelder-Mead', options=dict(maxiter=10000, xatol=1e-5, fatol=1e-3))
        if best is None or r.fun < best.fun: best = r

th = best.x
pmax, lx50, ls1, lxhi, ls2, binc = th

PF = dict(
    theta=dict(p_max=float(pmax), x50_px=float(np.exp(lx50)), sigma_lo=float(np.exp(ls1)),
               x_hi_px=float(np.exp(lxhi)), sigma_hi=float(np.exp(ls2)), b_inc=float(binc)),
    form="p(x,inc)=p_max*Phi((ln x-ln x50)/sigma_lo)*(1-Phi((ln x-ln x_hi)/sigma_hi))*exp(b_inc*(inc-60)/25); x=D_eq/px",
    likelihood="binomial over hollow-frame pairs conditional on >=1 detection per hollow (truncated)",
    catalog_completeness="C(D,site)=1-prod_f(1-p(D/px_f,inc_f)) over searched frames covering the site",
    n_pairs=int(len(CV)), n_hollows=int(len(uh)),
    caveats=["footprint model contains 83% of detecting frames",
             "hollows never detected are unobservable"])
json.dump(PF, open(OUTDIR + "per_frame_detection_function_inc60_75.json", "w"), indent=1)

# Build completeness grid
Dref = np.array([0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
bpx = np.array([np.sqrt(max(a, 4) * b) for a, b in BANDS])

C = np.ones((len(Dref), NY, NX), np.float32)
for bi in range(len(BANDS)):
    for di, D_r in enumerate(Dref):
        pb = pdet(D_r * 1000 / bpx[bi], 60.0, th)
        C[di] *= (1 - pb) ** nband[bi].astype(np.float32)
C = 1 - C

np.savez_compressed(OUTDIR + "completeness_10ppd.npz", C=C, Dref=Dref, ppd=PPD,
    note="Catalog completeness C(D,cell). Axis 0 = Dref (km). 10 ppd grid.")

# Decluster canonical hollows into fields
HOLa = pd.read_csv(OUTDIR + "hollows_canonical.csv")
hx_all = xyz(HOLa.lat.values, HOLa.lon.values)
T = cKDTree(hx_all)
pr = np.array(list(T.query_pairs(2 * np.sin(20.0 / (2 * RM_)))))
g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(HOLa), len(HOLa)))
nf, lab = connected_components(g, directed=False)
HOLa['field_id'] = lab

def cen(gr):
    v = xyz(gr.lat.values, gr.lon.values).mean(0)
    v /= np.linalg.norm(v)
    return pd.Series(dict(lat=np.degrees(np.arcsin(v[2])),
                          lon=np.degrees(np.arctan2(v[1], v[0])) % 360))

FC = HOLa.groupby('field_id').apply(cen)
FLD = HOLa.groupby('field_id').agg(
    n=('hollow_id', 'size'),
    n_nontrace=('is_crater_trace', lambda s: int((~s).sum())),
    D_max_km=('D_eq_km', 'max'),
    D_med_km=('D_eq_km', 'median'),
    area_sum_km2=('area_km2', 'sum'),
    px_med_m=('px_best', 'median'),
    px_min_m=('px_best', 'min'),
    inc_med=('inc', 'median'),
    conf_med=('conf', 'median')
).join(FC)

def extent(gr):
    if len(gr) < 2: return 0.0
    c = FC.loc[gr.name]
    d = RM_ * np.arccos(np.clip(
        np.sin(np.radians(c.lat)) * np.sin(np.radians(gr.lat)) +
        np.cos(np.radians(c.lat)) * np.cos(np.radians(gr.lat)) * np.cos(np.radians(gr.lon - c.lon)),
        -1, 1))
    return float(d.max())

FLD['extent_km'] = HOLa.groupby('field_id').apply(extent)

# Attach coverage / completeness
Z = np.load(OUTDIR + "completeness_10ppd.npz")
Cg = Z['C']
cv = np.load(OUTDIR + "coverage_10ppd.npz")

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

HOLa.to_csv(OUTDIR + "hollows_canonical.csv", index=False)
FLD.reset_index().to_csv(OUTDIR + "fields.csv", index=False)
print("canonical", len(HOLa), "fields", len(FLD), "detections_qc", len(DT), "completeness pairs", len(CV), "theta", np.round(th,3).tolist())

# ---------------------------------------------------------------- per-band frame + area report
# Northern searched area per 10-deg band, restricted (S) vs full (S_full).
def band_report(Sset, tag_):
    clat_ = Sset.CENTER_LATITUDE.values
    clon_ = Sset.CENTER_LONGITUDE.values % 360
    px_ = Sset.HORIZONTAL_PIXEL_SCALE.values
    emi_ = Sset.EMISSION_ANGLE.values
    grow_ = (1 / np.cos(np.radians(np.clip(emi_, 0, 80)))) ** 0.8
    half_km_ = 1024 * px_ / 1000 / 2 * grow_
    half_lat_ = half_km_ / KM_PER_DEG
    half_lon_ = half_lat_ / np.maximum(np.cos(np.radians(clat_)), 0.05)
    cov_ = np.zeros((NY, NX), bool)
    i0_ = np.clip(((90 - (clat_ + half_lat_)) * PPD).astype(int), 0, NY - 1)
    i1_ = np.clip(((90 - (clat_ - half_lat_)) * PPD).astype(int), 0, NY - 1)
    j0_ = ((clon_ - half_lon_) * PPD).astype(int)
    j1_ = ((clon_ + half_lon_) * PPD).astype(int)
    for k in range(len(Sset)):
        js = np.arange(j0_[k], j1_[k] + 1) % NX
        cov_[i0_[k]:i1_[k] + 1][:, js] = True
    # cell area (km^2) as function of row
    lat_top = np.radians(np.minimum(latc + 0.5 / PPD, 90))
    lat_bot = np.radians(np.maximum(latc - 0.5 / PPD, -90))
    row_area = RM_ * RM_ * np.radians(1.0 / PPD) * (np.sin(lat_top) - np.sin(lat_bot))
    out_ = []
    for b in range(0, 80, 10):
        rows_b = np.where((latc >= b) & (latc < b + 10))[0]
        area_cov = float((cov_[rows_b].sum(axis=1) * row_area[rows_b]).sum())
        nfr = int(((clat_ >= b) & (clat_ < b + 10)).sum())
        out_.append(dict(band="%d-%d" % (b, b + 10), n_frames_center=nfr,
                         covered_area_km2=area_cov))
    return out_

rep_restricted = band_report(S, "inc60_75")
rep_full = band_report(S_full, "full")
band_frames_area = dict(
    n_frames_full=int(len(S_full)), n_frames_inc60_75=int(len(S)),
    inc_lo=INC_LO, inc_hi=INC_HI,
    per_band_restricted=rep_restricted, per_band_full=rep_full)
json.dump(band_frames_area, open(OUTDIR + "band_frames_area.json", "w"), indent=1)
print("BAND REPORT (restricted): ", json.dumps(rep_restricted))
print("theta_inc60_75 b_inc =", float(binc))
