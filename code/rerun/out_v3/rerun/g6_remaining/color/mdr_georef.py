# g6_remaining copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/color_mdr/mdr_georef.py (6 Nov 2026). Changed: nothing (imported module).
"""Georeferencing and I/O for MESSENGER MDIS MDR 64-ppd tiles (PDS3, detached labels).

Pixel convention: PDS 1-based line/sample coordinates in which the center of pixel
(line L, sample S) is at coordinate (L, S); array index i = L-1, j = S-1.
Equirectangular (standard parallel = CENTER_LATITUDE):
    lat = (LINE_PROJECTION_OFFSET - L) * MAP_SCALE / R
    lon = CENTER_LONGITUDE + (S - SAMPLE_PROJECTION_OFFSET) * MAP_SCALE / (R cos(CENTER_LATITUDE))
Polar stereographic (true scale at the pole):
    north: rho = 2R tan(pi/4 - lat/2); x = rho sin(lon-lon0); y = -rho cos(lon-lon0)
    south: rho = 2R tan(pi/4 + lat/2); x = rho sin(lon-lon0); y = +rho cos(lon-lon0)
    S = SAMPLE_PROJECTION_OFFSET + x/MAP_SCALE ; L = LINE_PROJECTION_OFFSET - y/MAP_SCALE
Longitudes are returned 0-360 E (planetocentric, east-positive).
"""
import os, re, glob
import numpy as np

MDR_DIR = os.path.expanduser('~/mnt/hollows-emissivity/data/mdis_mdr')
R_KM = 2439.4
WAVES = np.array([430, 480, 560, 630, 750, 830, 900, 1000])
MISSING = -3.4028226550889045e+38
POLAR_LIM = 65.0
QC_R750_MIN = 0.035          # darkest equatorial terrain is ~0.057 at 750 nm
QC_430_750 = (0.40, 0.85)    # equatorial 0.01-99.99 pct: 0.57-0.73
QC_1000_750 = (0.95, 1.60)   # equatorial 0.01-99.99 pct: 1.03-1.26
LAST_QC = {}   # |lat| > 65 taken from polar tiles, <= 65 from equirectangular tiles


def parse_label(path):
    t = open(path).read()
    def g(k, num=True):
        m = re.search(r'^\s*' + k + r'\s*=\s*"?([^"\n<]+)', t, re.M)
        if m is None:
            return None
        v = m.group(1).strip()
        return float(v) if num else v
    d = {k: g(k) for k in ['RECORD_BYTES', 'FILE_RECORDS', 'LINES', 'LINE_SAMPLES', 'BANDS',
                           'CENTER_LATITUDE', 'CENTER_LONGITUDE', 'MAXIMUM_LATITUDE', 'MINIMUM_LATITUDE',
                           'WESTERNMOST_LONGITUDE', 'EASTERNMOST_LONGITUDE', 'LINE_PROJECTION_OFFSET',
                           'SAMPLE_PROJECTION_OFFSET', 'MAP_SCALE', 'MAP_RESOLUTION', 'A_AXIS_RADIUS']}
    d['PROJ'] = g('MAP_PROJECTION_TYPE', False)
    d['IMAGE_PTR'] = g(r'\^IMAGE', False)
    d['SAMPLE_TYPE'] = g('SAMPLE_TYPE', False)
    d['BAND_STORAGE_TYPE'] = g('BAND_STORAGE_TYPE', False)
    for k in ['LINES', 'LINE_SAMPLES', 'BANDS', 'RECORD_BYTES', 'FILE_RECORDS']:
        d[k] = int(d[k])
    return d


def wrap180(x):
    return (np.asarray(x) + 180.0) % 360.0 - 180.0


class Tile:
    def __init__(self, name):
        self.name = name
        self.lbl_path = os.path.join(MDR_DIR, 'MDIS_MDR_064PPD_%s.LBL' % name)
        self.img_path = os.path.join(MDR_DIR, 'MDIS_MDR_064PPD_%s.IMG' % name)
        L = parse_label(self.lbl_path)
        self.lbl = L
        self.polar = L['PROJ'].upper().startswith('POLAR')
        self.H, self.W = L['LINES'], L['LINE_SAMPLES']
        self.lpo, self.spo = L['LINE_PROJECTION_OFFSET'], L['SAMPLE_PROJECTION_OFFSET']
        self.scale_km = L['MAP_SCALE'] / 1000.0
        self.R = L['A_AXIS_RADIUS'] or R_KM
        self.clat, self.clon = L['CENTER_LATITUDE'], L['CENTER_LONGITUDE']
        self.north = self.polar and self.clat > 0
        # header check: detached label, ^IMAGE is a bare filename -> image starts at byte 0
        assert L['IMAGE_PTR'].endswith('.IMG') and ',' not in L['IMAGE_PTR']
        assert L['RECORD_BYTES'] == self.W * 4
        assert os.path.getsize(self.img_path) == self.H * self.W * 4 * L['BANDS']
        if not self.polar:
            W_ = L['WESTERNMOST_LONGITUDE']; E_ = L['EASTERNMOST_LONGITUDE']
            span = (E_ - W_) % 360.0
            self.lon_w = W_ % 360.0                 # nominal west edge, 0-360
            self.lon_span = float(np.round(span))   # nominal width (36 or 45 deg)
            self.lat_max = L['MAXIMUM_LATITUDE']
            self.lat_min = float(np.round(L['MINIMUM_LATITUDE'] * 4) / 4)  # nominal (e.g. 43.75)
            self.coslat_c = np.cos(np.radians(self.clat))
        self._mm = None

    # ---------------- I/O ----------------
    @property
    def mm(self):
        if self._mm is None:
            self._mm = np.memmap(self.img_path, dtype='<f4', mode='r',
                                 shape=(self.lbl['BANDS'], self.H, self.W))
        return self._mm

    # ---------------- projection ----------------
    def ll2pix(self, lat, lon):
        """lat, lon (deg; lon any convention) -> fractional array indices (i, j), 0-based, pixel centers at integers."""
        lat = np.asarray(lat, float); lon = np.asarray(lon, float)
        if not self.polar:
            y = np.radians(lat) * self.R
            x = np.radians(wrap180(lon - self.clon)) * self.R * self.coslat_c
        else:
            phi = np.radians(lat); dl = np.radians(lon - self.clon)
            if self.north:
                rho = 2 * self.R * np.tan(np.pi / 4 - phi / 2)
                x = rho * np.sin(dl); y = -rho * np.cos(dl)
            else:
                rho = 2 * self.R * np.tan(np.pi / 4 + phi / 2)
                x = rho * np.sin(dl); y = rho * np.cos(dl)
        S = self.spo + x / self.scale_km
        L = self.lpo - y / self.scale_km
        return L - 1.0, S - 1.0

    def pix2ll(self, i, j):
        """0-based array indices (pixel centers) -> lat, lon (deg, lon 0-360)."""
        L = np.asarray(i, float) + 1.0; S = np.asarray(j, float) + 1.0
        x = (S - self.spo) * self.scale_km
        y = (self.lpo - L) * self.scale_km
        if not self.polar:
            lat = np.degrees(y / self.R)
            lon = self.clon + np.degrees(x / (self.R * self.coslat_c))
        else:
            rho = np.hypot(x, y)
            c = 2 * np.arctan(rho / (2 * self.R))
            if self.north:
                lat = 90.0 - np.degrees(c)
                lon = self.clon + np.degrees(np.arctan2(x, -y))
            else:
                lat = -90.0 + np.degrees(c)
                lon = self.clon + np.degrees(np.arctan2(x, y))
        return lat, lon % 360.0

    def assigned(self, lat, lon):
        """True where (lat, lon) belongs to this tile's exclusive region (no double counting)."""
        lat = np.asarray(lat); lon = np.asarray(lon) % 360.0
        if self.polar:
            return lat > POLAR_LIM if self.north else lat < -POLAR_LIM
        dl = (lon - self.lon_w) % 360.0
        return (dl < self.lon_span) & (lat <= self.lat_max) & (lat > self.lat_min) & (np.abs(lat) <= POLAR_LIM)


def all_tiles(include_2700=False):
    names = sorted(os.path.basename(p)[16:-4] for p in glob.glob(os.path.join(MDR_DIR, 'MDIS_MDR_064PPD_*.LBL')))
    if not include_2700:
        names = [n for n in names if not n.startswith('2700')]
    return names


def load_refl(tile, other=None):
    """Return refl (8,H,W) float32 with NaN for invalid, count (H,W). If `other` (Tile on the same grid)
    is given, fill pixels invalid in `tile` from `other` (used for H15SP4 <- 2700_H15SP4)."""
    a = np.array(tile.mm[:9], dtype=np.float32)
    refl, cnt = a[:8], a[8]
    bad = ~np.all(np.isfinite(refl) & (refl > 0) & (refl < 2), axis=0)
    filled = np.zeros(bad.shape, bool)
    if other is not None:
        b = np.array(other.mm[:9], dtype=np.float32)
        ok2 = np.all(np.isfinite(b[:8]) & (b[:8] > 0) & (b[:8] < 2), axis=0)
        fill = bad & ok2
        refl[:, fill] = b[:8][:, fill]; cnt[fill] = b[8][fill]
        bad = bad & ~ok2; filled = fill
        del b
    # pixel QC (physically implausible values; mostly shadowed / extrapolated pixels poleward of ~70 deg)
    with np.errstate(invalid='ignore', divide='ignore'):
        q430 = refl[0] / refl[4]; q1000 = refl[7] / refl[4]
        qc_bad = ~bad & ~((refl[4] >= QC_R750_MIN) & (q430 >= QC_430_750[0]) & (q430 <= QC_430_750[1])
                          & (q1000 >= QC_1000_750[0]) & (q1000 <= QC_1000_750[1]))
    bad = bad | qc_bad
    refl[:, bad] = np.nan
    cnt = np.where(bad, np.nan, cnt)
    LAST_QC['n_qc_rejected'] = int(qc_bad.sum())
    return refl, cnt, filled


def gc_dist_km(lat1, lon1, lat2, lon2, R=R_KM):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


class Globe:
    """Nearest-pixel sampler over all tiles using each tile's exclusive region."""
    def __init__(self, band_idx=4):
        self.tiles = [Tile(n) for n in all_tiles()]
        self.t2700 = Tile('2700_H15SP4')
        self.band = band_idx

    def sample(self, lat, lon, band=None, source='auto'):
        band = self.band if band is None else band
        lat = np.atleast_1d(np.asarray(lat, float)); lon = np.atleast_1d(np.asarray(lon, float)) % 360
        out = np.full(lat.shape, np.nan, np.float32)
        for t in self.tiles:
            if source == 'auto':
                m = t.assigned(lat, lon)
            elif source == 'polar':
                m = np.ones(lat.shape, bool) if t.polar else np.zeros(lat.shape, bool)
                m &= (lat > 0) == t.north if t.polar else m
            else:  # equirect only, ignoring the 65-deg limit
                if t.polar:
                    continue
                dl = (lon - t.lon_w) % 360.0
                m = (dl < t.lon_span) & (lat <= t.lat_max) & (lat > t.lat_min)
            if not m.any():
                continue
            i, j = t.ll2pix(lat[m], lon[m])
            i = np.rint(i).astype(int); j = np.rint(j).astype(int)
            ok = (i >= 0) & (i < t.H) & (j >= 0) & (j < t.W)
            v = np.full(i.shape, np.nan, np.float32)
            v[ok] = t.mm[band, i[ok], j[ok]]
            if t.name == 'H15SP4':
                bad = ~(v > 0)
                if bad.any():
                    vv = np.full(bad.sum(), np.nan, np.float32)
                    okb = ok[bad]
                    vv[okb] = self.t2700.mm[band, i[bad][okb], j[bad][okb]]
                    v[bad] = vv
            v[~(v > 0)] = np.nan
            out[m] = v
        return out
