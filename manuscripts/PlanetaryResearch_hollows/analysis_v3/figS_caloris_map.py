import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
import os, json, numpy as np, pandas as pd, rasterio, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
HE=os.path.expanduser(_REPO + "")
BDR=f"{HE}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif"
COL=f"{HE}/GIS/raster/Mercury_MESSENGER_MDIS_Basemap_EnhancedColor_Mosaic_Global_665m_cog.tif"
RM=2439.4; KD=np.pi*RM/180; R0=2439400.0
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.5,"axes.linewidth":.6})
hc=pd.read_csv(f"{HE}/code/rerun/out_v2/hollows_canonical.csv")
cal=pd.read_csv(f"{HE}/peak_wall_test/caloris_children.csv"); m=pd.read_csv(f"{HE}/peak_wall_test/morph_all_nested.csv"); ch=pd.read_csv(f"{HE}/peak_wall_test/nested_children_classified.csv")
ha=pd.read_csv(f"{HE}/code/rerun/latitude/aspect/host_assignments.csv"); g=ha.groupby("host_uid").agg(n_in=("rho",lambda r:int((r<=0.9).sum())))
for d in (m,ch):
    d2=d.merge(g,left_on="uid",right_index=True,how="left"); d["n_in"]=d2.n_in.fillna(0).values; d["field"]=(d.n_in>=3).astype(int)
lon180=lambda x: ((x+180)%360)-180
def _read(ds,lo_a,lo_b,la_a,la_b,h,w,bands):
    win=from_bounds(np.radians(lo_a)*R0,np.radians(la_a)*R0,np.radians(lo_b)*R0,np.radians(la_b)*R0,ds.transform)
    return ds.read(bands,window=win,out_shape=(len(bands),h,w),resampling=Resampling.average,boundless=True)
def crop(src,lat0,lon0,half_km,maxpx=900,bands=None):
    lo0=lon180(lon0); cl=np.cos(np.radians(lat0)); dla=half_km/KD; dlo=half_km/(KD*cl)
    with rasterio.open(src) as ds:
        bands=bands or list(range(1,ds.count+1)); px_deg=ds.res[0]/R0*180/np.pi
        H=int(min(maxpx,2*dla/px_deg)); W=int(H*dlo/dla)
        a0,a1=lo0-dlo,lo0+dlo
        if a0<-180 or a1>180:   # wrap across the antimeridian
            if a1>180: parts=[(a0,180),(-180,a1-360)]
            else: parts=[(a0+360,180),(-180,a1)]
            ws=[int(round(W*(p[1]-p[0])/(2*dlo))) for p in parts]; ws[-1]=W-ws[0]
            a=np.concatenate([_read(ds,p[0],p[1],lat0-dla,lat0+dla,H,wi,bands) for p,wi in zip(parts,ws)],axis=2)
        else: a=_read(ds,a0,a1,lat0-dla,lat0+dla,H,W,bands)
    return a,[-half_km,half_km,-half_km,half_km]
def local(lat,lon,lat0,lon0): return (lon180(lon-lon0))*KD*np.cos(np.radians(lat0)), (lat-lat0)*KD
def gray(a):
    a=a[0].astype(float); p0,p1=np.percentile(a[a>0],[1,99.5]) if (a>0).any() else (0,1); return np.clip((a-p0)/(p1-p0),0,1)
def rgb(a):
    a=np.moveaxis(a.astype(float),0,-1); out=np.zeros_like(a)
    for k in range(3):
        b=a[...,k]; p0,p1=np.percentile(b[b>0],[1,99.5]) if (b>0).any() else (0,1); out[...,k]=np.clip((b-p0)/(p1-p0),0,1)
    return out
def circle(ax,x,y,r,**kw):
    t=np.linspace(0,2*np.pi,200); ax.plot(x+r*np.cos(t),y+r*np.sin(t),**kw)
# ---------------- Figure A: Caloris natural experiment

# ---------------- Supplementary figure: craters inside Caloris (map only; the counts are in Figure 1c of the main text)
OUT=f"{HE}/manuscripts/PlanetaryResearch_hollows/figures"
fig=plt.figure(figsize=(5.9,5.9)); ax=fig.add_axes([0.01,0.01,0.98,0.98])
lat0,lon0=31.93,164.1; half=950
a,ext=crop(BDR,lat0,lon0,half,maxpx=1600,bands=[1]); ax.imshow(gray(a),cmap="gray",extent=ext,origin="upper")
circle(ax,0,0,1595/2,color="w",lw=.9,ls=(0,(4,3)))
x,y=local(cal.lat.values,cal.lon.values,lat0,lon0)
below=cal.through==0; fld=cal.n_in>=3
ax.scatter(x[below],y[below],s=cal.D[below]*1.2,facecolors="none",edgecolors="#9ad0ff",lw=.7,label=f"excavation < 2.5 km ({int(below.sum())} craters)")
ax.scatter(x[~below&~fld],y[~below&~fld],s=cal.D[~below&~fld]*1.2,facecolors="none",edgecolors="#ffd166",lw=1.1,label=f"excavation ≥ 2.5 km, < 3 hollows ({int((~below&~fld).sum())})")
ax.scatter(x[fld],y[fld],s=cal.D[fld]*1.2,facecolors="#ff6b35",edgecolors="w",lw=.7,label=f"≥ 3 hollows within 0.9 radii ({int(fld.sum())})")
for nm in ["Sander","Nawahi","Apollodorus","Cunningham","Balanchine","Munch","Poe","Atget"]:
    r=cal[cal.name==nm]
    if len(r): xx,yy=local(r.lat.values[0],r.lon.values[0],lat0,lon0); ax.text(xx+20,yy+12,nm,color="w",fontsize=7)
ax.set_xlim(-half,half); ax.set_ylim(-half,half); ax.set_xticks([]); ax.set_yticks([])
ax.legend(loc="lower left",fontsize=7,frameon=True,facecolor="black",labelcolor="w",edgecolor="none",framealpha=.6,markerscale=.8)
ax.plot([half-325,half-25],[-half+40,-half+40],color="w",lw=2); ax.text(half-175,-half+60,"300 km",color="w",ha="center",fontsize=7)
fig.savefig(f"{OUT}/figS3_caloris_map.png",dpi=400); fig.savefig(f"{OUT}/figS3_caloris_map.pdf"); plt.close(fig)
print("wrote", len(cal), int(below.sum()), int((~below&~fld).sum()), int(fld.sum()), int((below&fld).sum()))
