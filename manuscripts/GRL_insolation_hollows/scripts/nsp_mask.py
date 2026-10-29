import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
import geopandas as gpd, numpy as np, pyproj
from scipy import ndimage
U=_REPO + ""
g=gpd.read_file(f"{U}/data/data/NSP_map_data/MLA_North_Smooth_Plains_lin_V0.shp")
src=pyproj.CRS.from_wkt(open(f"{U}/data/data/NSP_map_data/MLA_North_Smooth_Plains_lin_V0.prj").read())
geo=pyproj.CRS.from_proj4("+proj=longlat +R=2439400 +no_defs")
tr=pyproj.Transformer.from_crs(src,geo,always_xy=True)
# densify each line to ~2 km spacing and transform
xs=[];ys=[]
for geom in g.geometry:
    L=geom.length; n=max(2,int(L/2000)+1)
    pts=[geom.interpolate(t,normalized=True) for t in np.linspace(0,1,n)]
    xs+= [p.x for p in pts]; ys+=[p.y for p in pts]
lon,lat=tr.transform(np.array(xs),np.array(ys))
lon=np.mod(lon,360)
print("lineament points",len(lon),"lat range",lat.min(),lat.max())
ppd=10; H=180*ppd; W=360*ppd
r=((90-lat)*ppd).astype(int).clip(0,H-1); c=(lon*ppd).astype(int).clip(0,W-1)
cnt=np.zeros((H,W),np.int32); np.add.at(cnt,(r,c),1)
np.savez_compressed("lineament_count_10ppd.npz",cnt=cnt)
