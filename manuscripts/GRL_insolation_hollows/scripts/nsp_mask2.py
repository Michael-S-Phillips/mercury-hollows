import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
import numpy as np, pyproj
from scipy import ndimage
import geopandas as gpd
U=_REPO + ""
R=2439400.0
# re-derive points in a north polar stereographic frame (true scale at pole) for morphology
g=gpd.read_file(f"{U}/data/data/NSP_map_data/MLA_North_Smooth_Plains_lin_V0.shp")
src=pyproj.CRS.from_wkt(open(f"{U}/data/data/NSP_map_data/MLA_North_Smooth_Plains_lin_V0.prj").read())
ps=pyproj.CRS.from_proj4(f"+proj=laea +lat_0=90 +lon_0=0 +R={R} +no_defs")   # equal-area polar
tr=pyproj.Transformer.from_crs(src,ps,always_xy=True)
xs=[];ys=[]
for geom in g.geometry:
    n=max(2,int(geom.length/2000)+1)
    for t in np.linspace(0,1,n):
        p=geom.interpolate(t,normalized=True); xs.append(p.x); ys.append(p.y)
X,Y=tr.transform(np.array(xs),np.array(ys))
pix=5000.0; half=3_200_000; n=int(2*half/pix)
ci=((X+half)/pix).astype(int); ri=((half-Y)/pix).astype(int)
occ=np.zeros((n,n),bool); occ[ri.clip(0,n-1),ci.clip(0,n-1)]=True
out={}
for rad_km in [25,40,60]:
    k=int(rad_km*1000/pix); yy,xx=np.mgrid[-k:k+1,-k:k+1]; st=(xx**2+yy**2)<=k*k
    m=ndimage.binary_dilation(occ,structure=st)
    m=ndimage.binary_fill_holes(m)
    m=ndimage.binary_erosion(m,structure=st)   # closing: dilate, fill, erode back
    lab,nl=ndimage.label(m); sizes=ndimage.sum(m,lab,range(1,nl+1))
    keep=np.isin(lab,1+np.where(sizes*pix*pix/1e6>2e4)[0])   # drop specks <20,000 km2
    out[rad_km]=keep
    print(rad_km,"km closing: area %.2f e6 km2, %d parts kept of %d"%(keep.sum()*pix*pix/1e12, len(np.unique(lab[keep])), nl))
np.savez_compressed("nsp_laea.npz",**{f"m{k}":v for k,v in out.items()},occ=occ,pix=pix,half=half)
# map to 10 ppd lat/lon grid
ppd=10; H=180*ppd; W=360*ppd
lat=90-(np.arange(H)+0.5)/ppd; lon=(np.arange(W)+0.5)/ppd
LON,LAT=np.meshgrid(lon,lat)
geo=pyproj.CRS.from_proj4(f"+proj=longlat +R={R} +no_defs")
t2=pyproj.Transformer.from_crs(geo,ps,always_xy=True)
north=LAT>20
Xg,Yg=t2.transform(LON[north],LAT[north])
cg=((Xg+half)/pix).astype(int); rg=((half-Yg)/pix).astype(int); ok=(cg>=0)&(cg<n)&(rg>=0)&(rg<n)
for k,m in out.items():
    G=np.zeros((H,W),bool); v=np.zeros(north.sum(),bool); v[ok]=m[rg[ok],cg[ok]]; G[north]=v
    a=np.cos(np.radians(LAT))*(np.pi*2439.4/180/ppd)**2
    print(k,"10ppd mask area %.2f e6 km2, lat range %.1f-%.1f"%((a*G).sum()/1e6, LAT[G].min(), LAT[G].max()))
    np.save(f"nsp_mask_10ppd_r{k}.npy",G)
