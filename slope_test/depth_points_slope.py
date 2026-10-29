import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Regional slope and facing azimuth at each Blewett et al. (2016) depth measurement, from the USGS 665 m DEM.
s1: central differences on the native grid (1.33 km baseline); s3: after a 3 x 3 mean filter (about 2.7 km baseline).
East-west spacing scaled by cos(latitude). Writes depth_points_slope.csv (row order matches test1_rows.csv)."""
import os, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
from scipy.ndimage import uniform_filter
HE = os.path.expanduser(_REPO + "")
DEM = f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
d = pd.read_csv(f"{HE}/manuscripts/GRL_insolation_hollows/data/test1_rows.csv")
PX = 665.24; out = []
with rasterio.open(DEM) as ds:
    for r in d.itertuples():
        lon = r.longitude % 360; lat = r.latitude
        x = np.radians(lon - 180) * 2439400.0; y = np.radians(lat) * 2439400.0
        col, row = ~ds.transform * (x, y); col, row = int(np.floor(col)), int(np.floor(row))
        z = ds.read(1, window=Window(col - 4, row - 4, 9, 9), boundless=True, fill_value=-32768).astype(float)
        z[z == -32768] = np.nan; z *= 0.5
        res = dict(row=r.Index)
        for bl, zz in (("s1", z), ("s3", uniform_filter(np.where(np.isfinite(z), z, np.nanmean(z)), 3, mode="nearest"))):
            gr, gc = np.gradient(zz)
            dzdn = -gr[4, 4] / PX; dzde = gc[4, 4] / (PX * np.cos(np.radians(lat)))
            res[f"slope_{bl}"] = float(np.degrees(np.arctan(np.hypot(dzde, dzdn))))
            res[f"face_az_{bl}"] = float(np.degrees(np.arctan2(-dzde, -dzdn)) % 360)
        res["z_m"] = float(z[4, 4]); out.append(res)
o = pd.DataFrame(out); o.to_csv(f"{HE}/slope_test/depth_points_slope.csv", index=False)
print(len(o), o[["slope_s1", "slope_s3"]].describe().round(2))
