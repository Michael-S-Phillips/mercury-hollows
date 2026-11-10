# g9 copy of g3_catalog/melt_source/knobs/common.py (6 Nov 2026). Changed: env G9_CAT/G9_GRID, output folder knobs/out_<CAT>_grid<GRID>, and the completeness grid in load_C (v3 = g4 completeness_10ppd_v3.npz).
# g3_catalog copy of analysis_workshop/melt_source_tests/knobs/common.py (6 Nov 2026): only input/output paths and env settings changed.
"""Shared helpers for the circum-Caloris knob tests (5 Nov 2026): paths, catalogs, local grids, and samplers for
MDIS MDR color (64 ppd), the georeferenced LRM maps, the 1 km hollow completeness grid and the USGS 665 m DEM.
Workshop analysis only; reads existing files and writes nothing outside this folder."""
import os, sys, numpy as np, pandas as pd, rasterio, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.expanduser("~/mnt/hollows-emissivity")
AW = f"{ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_workshop"
G9CAT = os.environ.get("G9_CAT", "v3"); G9GRID = os.environ.get("G9_GRID", "v3")
HERE = f"{ROOT}/code/rerun/out_v3/rerun/g9_draft_reruns/knobs/out_{G9CAT}_grid{G9GRID}"; os.makedirs(HERE, exist_ok=True)
ORIG = f"{AW}/melt_source_tests/knobs"
CANON = (f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv" if G9CAT == "v2"
         else f"{ROOT}/code/rerun/out_v3/hollows_canonical_" + G9CAT + ".csv")
sys.path.insert(0, f"{AW}/color_mdr"); sys.path.insert(0, f"{AW}/province_basins/wright2020_knobs")
import mdr_georef as M
from read_dbf import read_dbf
RM = 2439.4; KD = np.pi * RM / 180.0
CAL = dict(lat=29.91, lon=161.43, D=1556.0)          # Orgel et al. (2020)
WAVES = [430, 480, 560, 630, 750, 830, 900, 1000]     # MDR bands 0-7 (label BAND_NAME); band 8+ = stdev/count
DEM = f"{ROOT}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
LRM_FILES = {"conc": f"{ROOT}/provenance_pilot/lrm/Concentrated_LRM_560_params_BW_georef.tif",   # Klima et al. 2018 def.
             "old": f"{ROOT}/provenance_pilot/lrm/Old_LRM_newmos_georef.tif"}                     # Peplowski et al. 2016 def.

def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]

def gcd(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2); dl = np.radians(np.asarray(lo2) - np.asarray(lo1))
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * RM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

def inv_aeqd(la0, lo0, dx, dy):
    """Azimuthal equidistant inverse: local east/north offsets (km) -> lat, lon (0-360)."""
    p0 = np.radians(la0); rho = np.hypot(dx, dy) / RM; az = np.arctan2(dx, dy)
    lat = np.arcsin(np.clip(np.sin(p0) * np.cos(rho) + np.cos(p0) * np.sin(rho) * np.cos(az), -1, 1))
    lon = np.radians(lo0) + np.arctan2(np.sin(az) * np.sin(rho) * np.cos(p0), np.cos(rho) - np.sin(p0) * np.sin(lat))
    return np.degrees(lat), np.degrees(lon) % 360.0

def local_grid(la0, lo0, rmax, s):
    """Square grid (spacing s km) clipped to a disk of radius rmax around (la0, lo0). Returns r_km, lat, lon."""
    n = int(np.ceil(rmax / s)); v = (np.arange(-n, n + 1)) * s
    dx, dy = np.meshgrid(v, v); dx = dx.ravel(); dy = dy.ravel(); r = np.hypot(dx, dy); k = r <= rmax
    la, lo = inv_aeqd(la0, lo0, dx[k], dy[k]); return r[k], la, lo

def near(la0, lo0, la, lo, rad):
    """Indices of catalog points (la, lo) whose great-circle distance to (la0, lo0) is <= rad (km)."""
    dla = rad / KD + 0.01
    pre = np.where(np.abs(la - la0) <= dla)[0]
    d = gcd(la0, lo0, la[pre], lo[pre]); return pre[d <= rad], d[d <= rad]

def inside_any(la, lo, fla, flo, frad):
    """Boolean: point (la, lo) lies within frad of any feature (fla, flo)."""
    m = np.zeros(len(la), bool)
    for a, b, r in zip(fla, flo, frad):
        m |= gcd(a, b, la, lo) <= r
    return m

# ---------------- catalogs ----------------
def load_knobs():
    k = read_dbf(f"{AW}/province_basins/wright2020_knobs/CircumCalorisKnobs_SupplementaryShapefile.dbf")
    k["lon360"] = k.Longitude % 360
    k["a_km"] = np.where(k.LongAxis_k > 0, k.LongAxis_k / 2, np.where(k.KnobType == "Knob (<5 km)", 2.0, 3.5))
    k["main"] = k.KnobType.isin(["Knob (certain)", "Knob (uncertain)"]) & (k.LongAxis_k >= 5)
    k["dcal_R"] = gcd(CAL["lat"], CAL["lon"], k.Latitude.values, k.lon360.values) / (CAL["D"] / 2)
    return k.reset_index(drop=True)

def load_craters():
    c = pd.read_csv(f"{ROOT}/data/ME6M300TGT/herrick2018_all.csv"); c["lon360"] = c.lon % 360
    c["cent_str"] = c.cent_str.fillna("").astype(str).str.strip(); return c

def load_hollows():
    h = pd.read_csv(CANON, usecols=["lat", "lon", "hollow_id", "field_id", "is_crater_trace", "D_eq_km"])
    h = h[~h.is_crater_trace.astype(bool)].reset_index(drop=True); h["lon360"] = h.lon % 360; return h

def load_C(dref=1.0):
    C = np.load(f"{ROOT}/code/rerun/out_v2/completeness_10ppd.npz" if G9GRID == "v2" else f"{ROOT}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz")
    return C["C"][list(C["Dref"]).index(dref)]          # rows from 90N (row = (90-lat)*10), cols from 0E (col = lon*10)

def C_at(Cg, la, lo):
    ci = np.clip(((90 - np.asarray(la)) * 10).astype(int), 0, 1799); cj = ((np.asarray(lo) % 360) * 10).astype(int) % 3600
    w = Cg[ci, cj]; return np.where(np.isfinite(w), w, 0.0)

# ---------------- LRM (georeferenced copies; 0.12 deg = 5.1 km pixels; 8-bit index, 0 = not LRM) ----------------
class LRM:
    def __init__(self):
        self.a = {}
        for k, f in LRM_FILES.items():
            with rasterio.open(f) as d:
                self.a[k] = d.read(1); self.T = d.transform
        self.ny, self.nx = self.a["conc"].shape
    def at(self, la, lo, key):
        x = np.radians(((np.asarray(lo) + 180) % 360) - 180) * RM * 1000; y = np.radians(np.asarray(la)) * RM * 1000
        c = ((x - self.T.c) / self.T.a).astype(int) % self.nx; r = np.clip(((y - self.T.f) / self.T.e).astype(int), 0, self.ny - 1)
        return self.a[key][r, c]

# ---------------- MDR color ----------------
class MDR:
    """Pixel-level access to the MDR 64 ppd tiles (nearest pixel, each tile's exclusive region, as mdr_georef.Globe)."""
    def __init__(self):
        self.G = M.Globe(); self.tiles = self.G.tiles
    def keys(self, la, lo):
        """Unique MDR pixels touched by points (la, lo): returns tile index, i, j arrays."""
        la = np.atleast_1d(la).astype(float); lo = np.atleast_1d(lo).astype(float) % 360
        T, I, J = [], [], []
        for ti, t in enumerate(self.tiles):
            m = t.assigned(la, lo)
            if not m.any(): continue
            i, j = t.ll2pix(la[m], lo[m]); i = np.rint(i).astype(np.int64); j = np.rint(j).astype(np.int64)
            ok = (i >= 0) & (i < t.H) & (j >= 0) & (j < t.W)
            T.append(np.full(ok.sum(), ti)); I.append(i[ok]); J.append(j[ok])
        if not T: return np.zeros(0, int), np.zeros(0, int), np.zeros(0, int)
        T, I, J = np.concatenate(T), np.concatenate(I), np.concatenate(J)
        u = np.unique(np.c_[T, I, J], axis=0); return u[:, 0], u[:, 1], u[:, 2]
    def center_ll(self, T, I, J):
        la = np.empty(len(T)); lo = np.empty(len(T))
        for ti in np.unique(T):
            m = T == ti; la[m], lo[m] = self.tiles[ti].pix2ll(I[m], J[m])
        return la, lo
    def spectra(self, T, I, J):
        """(n, 8) reflectance with the color_mdr pixel QC applied (failing pixels -> NaN rows)."""
        S = np.full((len(T), 8), np.nan, np.float32)
        for ti in np.unique(T):
            m = np.where(T == ti)[0]; t = self.tiles[ti]
            o = np.argsort(I[m] * 100000 + J[m]); mm = m[o]
            v = np.asarray(t.mm[:8, I[mm], J[mm]]).T
            if t.name == "H15SP4":
                bad = ~np.all(v > 0, axis=1)
                if bad.any(): v[bad] = np.asarray(self.G.t2700.mm[:8, I[mm][bad], J[mm][bad]]).T
            S[mm] = v
        with np.errstate(invalid="ignore", divide="ignore"):
            ok = np.all(np.isfinite(S) & (S > 0) & (S < 2), axis=1)
            ok &= (S[:, 4] >= M.QC_R750_MIN) & (S[:, 0] / S[:, 4] >= M.QC_430_750[0]) & (S[:, 0] / S[:, 4] <= M.QC_430_750[1])
            ok &= (S[:, 7] / S[:, 4] >= M.QC_1000_750[0]) & (S[:, 7] / S[:, 4] <= M.QC_1000_750[1])
        S[~ok] = np.nan; return S

# PCA of brightness-normalized log spectra from color_mdr/pca_record.json (crater interiors, global)
import json as _json
_P = _json.load(open(f"{AW}/color_mdr/pca_record.json"))
MEAN_NL = np.array([_P["mean_nl"][str(w)] for w in WAVES]); L1 = np.array([_P["loadings"]["PC1"][str(w)] for w in WAVES])
L2 = np.array([_P["loadings"]["PC2"][str(w)] for w in WAVES])

def indices(spec):
    """Color indices of one region from its (n, 8) pixel spectra: per-band median (as color_mdr), then indices."""
    S = spec[np.all(np.isfinite(spec), axis=1)]
    out = {"npx": int(len(S))}
    if len(S) < 3:
        for k in ["r560", "r750", "albedo_gm", "ln_1000_480", "slope_750_430", "pc1", "pc2"]: out[k] = np.nan
        return out
    r = np.median(S, axis=0); lr = np.log(r); nl = lr - lr.mean()
    out.update(r560=float(r[2]), r750=float(r[4]), albedo_gm=float(np.exp(lr.mean())), ln_1000_480=float(lr[7] - lr[1]),
               slope_750_430=float(lr[4] - lr[0]), pc1=float((nl - MEAN_NL) @ L1), pc2=float((nl - MEAN_NL) @ L2))
    return out

# ---------------- DEM (USGS 665 m; x = radians(lon-180) R, y = radians(lat) R; int16 x 0.5 m as in peak_wall.py) ----------------
class DEMWin:
    def __init__(self, lat0, lat1, lon0, lon1):
        from rasterio.windows import from_bounds
        xs = lambda lo: np.radians(lo - 180) * 2439400.0; ys = lambda la: np.radians(la) * 2439400.0
        with rasterio.open(DEM) as ds:
            win = from_bounds(xs(lon0), ys(lat0), xs(lon1), ys(lat1), ds.transform)
            z = ds.read(1, window=win, boundless=True, fill_value=-32768).astype(np.float32); self.T = ds.window_transform(win)
        z[z == -32768] = np.nan; self.z = z * 0.5; self.ny, self.nx = z.shape
    def rc(self, la, lo):
        x = np.radians(np.asarray(lo) % 360 - 180) * 2439400.0; y = np.radians(np.asarray(la)) * 2439400.0
        c = np.floor((x - self.T.c) / self.T.a).astype(int); r = np.floor((y - self.T.f) / self.T.e).astype(int)
        return r, c
    def at(self, la, lo):
        r, c = self.rc(la, lo); ok = (r >= 0) & (r < self.ny) & (c >= 0) & (c < self.nx)
        out = np.full(len(r), np.nan, np.float32); out[ok] = self.z[r[ok], c[ok]]; return out
