"""Shared loaders for the hollow shadow-depth pilot (depth_expansion_pilot).

Inputs (ROOT = hollows-emissivity root, default ~/mnt/hollows-emissivity):
- Blewett et al. (2016) shadow-length depths: ROOT/manuscripts/insolation_thermal/data/test1_rows_inc.csv
- Bickel et al. (2024) release CSV (0.7 confidence, 19,110 detections): ROOT/data/bickel-2024/Bickel_et_al_2024_MDIS_hollows_07_v1.0.csv
  (row i of the release = fid i+1 of the GeoPackage; checked on MDIS ID and box centre)
- Camera-model positions and box corners: ROOT/GIS/vector/bickel2024_hollows_corrected.gpkg
  (hollows_points_corrected; hollows_bbox_cameramodel, corner order (x0,y0),(x1,y0),(x1,y1),(x0,y1), checked against a
  frame-level fit)
- Herrick et al. (2018) craters: ROOT/manuscripts/insolation_thermal/data/test2_craters_h18.csv
- MDIS NAC EDRs (PDS MSGRMDS_1001) from the PDS Imaging Node mirror
  https://planetarydata.jpl.nasa.gov/img/data/messenger/MSGRMDS_1001/DATA/YYYY_DOY/EN<num>M.IMG, cached outside the
  project (EDR_CACHE, default ~/edr_cache on the VM); the directory is the UTC date of the Bickel 'Earth datetime'
  (the neighbouring days are tried when the file is not there).
"""
import os, re, struct, sqlite3, subprocess, datetime as dt
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

ROOT = os.path.expanduser(os.environ.get("HOLLOWS_ROOT", "~/mnt/hollows-emissivity"))
HERE = os.path.dirname(os.path.abspath(__file__))
F_BLEWETT = f"{ROOT}/manuscripts/insolation_thermal/data/test1_rows_inc.csv"
F_BICKEL = f"{ROOT}/data/bickel-2024/Bickel_et_al_2024_MDIS_hollows_07_v1.0.csv"
F_GPKG = f"{ROOT}/GIS/vector/bickel2024_hollows_corrected.gpkg"
F_CRATERS = f"{ROOT}/manuscripts/insolation_thermal/data/test2_craters_h18.csv"
EDR_URL = "https://planetarydata.jpl.nasa.gov/img/data/messenger/MSGRMDS_1001/DATA"
EDR_CACHE = os.path.expanduser(os.environ.get("EDR_CACHE", "~/edr_cache"))
R_M = 2439.7  # km, as in the paper's host assignment


def load_blewett():
    b = pd.read_csv(F_BLEWETT)
    b["num"] = b.filename.str[2:12]
    b["bid"] = np.arange(len(b))
    return b


def _gp_coords(gb):
    flags = gb[3]; env = (flags >> 1) & 7; n = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]; w = gb[8 + n:]
    bo = "<" if w[0] == 1 else ">"; npt = struct.unpack(bo + "I", w[9:13])[0]
    return np.array(struct.unpack(bo + "%dd" % (2 * npt), w[13:13 + 16 * npt])).reshape(-1, 2)


def _poly_eqdiam_km(wkt):
    """Equivalent-circle diameter (km) of the release polygon (lon/lat degrees), local equirectangular area."""
    try:
        nums = np.array(re.findall(r"-?\d+\.?\d*(?:[eE]-?\d+)?", wkt), float).reshape(-1, 2)
    except Exception:
        return np.nan
    lon, lat = nums[:, 0], nums[:, 1]; la0 = np.radians(lat.mean())
    x = R_M * np.radians(lon - lon.mean()) * np.cos(la0); y = R_M * np.radians(lat - lat.mean())
    a = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return 2 * np.sqrt(a / np.pi)


def load_bickel(with_corners=True):
    k = pd.read_csv(F_BICKEL, sep=";", encoding="latin1")
    k.columns = ["mdis_id", "datetime", "emission_rel", "incidence_rel", "phase_rel", "px_rel", "subsol_az_rel", "subsol_lon",
                 "subsol_lat", "sol_dist", "exposure", "saturated", "lon_rel", "lat_rel", "x0", "x1", "y0", "y1", "conf_raw",
                 "stage", "thumb", "polygon"]
    k["confidence"] = pd.to_numeric(k.conf_raw.astype(str).str.split().str[-1], errors="coerce")  # merged fields
    k["fid"] = np.arange(1, len(k) + 1)
    c = sqlite3.connect(F_GPKG)
    g = pd.read_sql("select fid, mdis_id as mid2, lat, lon180, frame_samples, frame_lines, img_x, img_y from hollows_points_corrected", c)
    k = k.merge(g, on="fid"); assert (k.mdis_id == k.mid2).all() and np.allclose(k.img_x, (k.x0 + k.x1) / 2)
    k["num"] = k.mdis_id.str[2:12]
    k["poly_eqd_km"] = [_poly_eqdiam_km(w) for w in k.polygon.astype(str)]
    if with_corners:
        cor = {r[0]: _gp_coords(r[1])[:4] for r in c.execute("select fid, geom from hollows_bbox_cameramodel")}
        k["corners"] = k.fid.map(cor)
    return k


class FrameMap:
    """Affine map lon/lat (local tangent plane, km) -> raw EDR pixel (sample, line) for one frame, fitted to all camera-model
    box corners and centres of the Bickel detections in that frame (rms is reported)."""
    def __init__(self, q):
        N = q.frame_samples.values[0]; L = q.frame_lines.values[0]
        self.lat0, self.lon0 = q.lat.mean(), q.lon180.mean()
        E, Pp = [], []
        for _, r in q.iterrows():
            E.append([r.lon180, r.lat]); Pp.append([r.img_x * N, r.img_y * L])
            if isinstance(r.corners, np.ndarray):
                for (lo, la), (u, v) in zip(r.corners, [(r.x0, r.y0), (r.x1, r.y0), (r.x1, r.y1), (r.x0, r.y1)]):
                    E.append([lo, la]); Pp.append([u * N, v * L])
        E = np.array(E); Pp = np.array(Pp); A = np.column_stack([self.en(E[:, 0], E[:, 1]), np.ones(len(E))])
        self.coef = np.linalg.lstsq(A, Pp, rcond=None)[0]
        self.rms = float(np.sqrt(((A @ self.coef - Pp) ** 2).sum(1).mean()))
        self.km_per_px = float(np.sqrt(abs(np.linalg.det(np.linalg.inv(self.coef[:2])))))

    def en(self, lon, lat):
        dlon = (np.asarray(lon) - self.lon0 + 180) % 360 - 180
        return np.column_stack([R_M * np.radians(dlon) * np.cos(np.radians(self.lat0)), R_M * np.radians(np.asarray(lat) - self.lat0)])

    def to_px(self, lon, lat):
        E = self.en(lon, lat); return np.column_stack([E, np.ones(len(E))]) @ self.coef


# ------------------------------------------------------------------ EDR access
def edr_path(num):
    return os.path.join(EDR_CACHE, f"EN{num}M.IMG")


def fetch_edr(num, datetime_str):
    os.makedirs(EDR_CACHE, exist_ok=True); p = edr_path(num)
    if os.path.exists(p) and os.path.getsize(p) > 200000:
        return p
    t = dt.datetime.fromisoformat(datetime_str[:19])
    for dd in (0, -1, 1):
        d = t + dt.timedelta(days=dd); url = f"{EDR_URL}/{d.year}_{d.timetuple().tm_yday:03d}/EN{num}M.IMG"
        r = subprocess.run(["curl", "-s", "-f", "-m", "60", "-o", p + ".part", url])
        if r.returncode == 0 and os.path.getsize(p + ".part") > 200000:
            os.replace(p + ".part", p); return p
    if os.path.exists(p + ".part"):
        os.remove(p + ".part")
    return None


def read_edr(path):
    """Return (image float32 [lines, samples], label dict). Raw EDR orientation, 8-bit DN (companded)."""
    raw = open(path, "rb").read()
    end = raw.find(b"\r\nEND\r\n"); end = raw.find(b"\nEND\n") if end < 0 else end
    lab = raw[:end].decode("latin1")
    kw = {}
    for m in re.finditer(r"^\s*([A-Z0-9_:\^]+)\s*=\s*(.+?)\s*$", lab, re.M):
        kw.setdefault(m.group(1), m.group(2))   # first occurrence (outer level before IMAGE object)
    def num(key):
        v = kw.get(key); 
        if v is None: return np.nan
        try: return float(re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", v)[0])
        except Exception: return np.nan
    rb = int(num("RECORD_BYTES")); ptr = int(num("^IMAGE"))
    obj = lab[lab.find("OBJECT = IMAGE"):] if "OBJECT = IMAGE" in lab else lab
    L = int(re.search(r"LINES\s*=\s*(\d+)", obj).group(1)); S = int(re.search(r"LINE_SAMPLES\s*=\s*(\d+)", obj).group(1))
    bits = int(re.search(r"SAMPLE_BITS\s*=\s*(\d+)", obj).group(1)); assert bits == 8
    dark = re.search(r"DARK_STRIP_MEAN\s*=\s*([-\d.]+)", obj)
    img = np.frombuffer(raw, np.uint8, count=L * S, offset=(ptr - 1) * rb).reshape(L, S).astype(np.float32)
    info = dict(lines=L, samples=S, subsol_az=num("SUB_SOLAR_AZIMUTH"), subsc_az=num("SUB_SPACECRAFT_AZIMUTH"),
                north_az=num("NORTH_AZIMUTH"), inc=num("INCIDENCE_ANGLE"), emi=num("EMISSION_ANGLE"), phase=num("PHASE_ANGLE"),
                px=num("HORIZONTAL_PIXEL_SCALE"), subsol_lat=num("SUB_SOLAR_LATITUDE"), subsol_lon=num("SUB_SOLAR_LONGITUDE"),
                center_lat=num("CENTER_LATITUDE"), center_lon=num("CENTER_LONGITUDE"), dark=float(dark.group(1)) if dark else np.nan,
                pixelbin=num("MESS:PIXELBIN"), comp=num("MESS:COMP12_8"), start=kw.get("START_TIME", ""))
    return img, info


def local_incidence(lat, lon, sslat, sslon):
    la, lo, sa, so = map(np.radians, (lat, lon, sslat, sslon))
    return np.degrees(np.arccos(np.clip(np.sin(la) * np.sin(sa) + np.cos(la) * np.cos(sa) * np.cos(lo - so), -1, 1)))


# ------------------------------------------------------------------ host craters (depth_nested.py / depth_structure.py)
def assign_hosts(lat, lon360):
    cr = pd.read_csv(F_CRATERS); cr = cr[cr.diameter >= 5].reset_index(drop=True)
    def xyz(la, lo):
        la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    tree = cKDTree(xyz(cr.latitude.values, cr.lon360.values)); Pt = xyz(np.asarray(lat), np.asarray(lon360) % 360)
    Dkm = cr.diameter.values; host = np.full(len(Pt), -1); rho = np.full(len(Pt), np.nan)
    for i, cl in enumerate(tree.query_ball_point(Pt, r=0.36)):
        if not cl: continue
        cl = np.array(cl); ang = np.arccos(np.clip(Pt[i] @ tree.data[cl].T, -1, 1)) * R_M
        ok = ang <= 1.1 * Dkm[cl] / 2
        if ok.any():
            j = np.argmin(np.where(ok, Dkm[cl], np.inf)); host[i] = cl[j]; rho[i] = ang[j] / (Dkm[cl[j]] / 2)
    hid = np.where(host >= 0, cr.unique_id.values[np.maximum(host, 0)].astype(str), "")
    hD = np.where(host >= 0, Dkm[np.maximum(host, 0)], np.nan)
    return hid, hD, rho
