import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Sub-solar point for each image in the Blewett et al. (2016) depth set, from the MDIS CDR cumulative index (fixed-width rows,
1084 bytes; FILE_NAME bytes 33-54, SUB_SOLAR_LATITUDE 930-938, SUB_SOLAR_LONGITUDE 940-948, INCIDENCE_ANGLE 950-958; 1-based).
Writes subsolar_depth_images.csv (key, subsolar lat/lon, image-center incidence)."""
import pandas as pd
d = pd.read_csv(_REPO + "/manuscripts/GRL_insolation_hollows/data/test1_rows_inc.csv")
keys = set(d.filename.str.slice(1, 13).str.upper())
rows = {}
with open(_REPO + "/external/mdis_index/CDR_CUMINDEX.TAB", "rb") as f:
    for line in f:
        fn = line[32:54].decode().strip().strip('"')
        k = fn[1:13].upper()
        if k in keys and k not in rows:
            rows[k] = dict(key=k, file=fn, sslat=float(line[929:938]), sslon=float(line[939:948]), inc_center=float(line[949:958]))
out = pd.DataFrame(rows.values()); out.to_csv("subsolar_depth_images.csv", index=False)
print(len(keys), "images needed;", len(out), "found"); print(out.head())
