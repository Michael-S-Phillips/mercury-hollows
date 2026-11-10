# g5 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/caloris_v9/figS_caloris_map_v9.py (6 Nov 2026). Changed: crater table = g5_figures/_scratch/caloris_v3/caloris_children_v9.csv (v3 hosting, from caloris_test_record_v9_v3cat.py); output g5_figures/figS3_caloris_map_v9_v3cat.png/.jpg/.pdf.
"""v9 copy of analysis_v3/figS_caloris_map.py (Figure S3). Craters >= 10 km inside Caloris in the v9 crater containment graph
(179 craters, caloris_children_v9.csv from caloris_test_record_v9.py), on the 166 m/pixel MDIS mosaic. v9 membership uses the
union of the Orgel et al. (2020) outline (dashed; the merged Caloris node's geometry) and the Herrick et al. (2018) catalog
outline (dotted, the circle drawn in the v5 figure). Writes figures/figS3_caloris_map_v9.png/.jpg/.pdf without overwriting."""
import os, numpy as np, pandas as pd, rasterio, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
from PIL import Image
HE=os.path.expanduser("~/mnt/hollows-emissivity")
BDR=f"{HE}/data/Mercury_MESSENGER_MDIS_Basemap_BDR_Mosaic_Global_166m.tif"
RM=2439.4; KD=np.pi*RM/180; R0=2439400.0
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.5,"axes.linewidth":.6})
cal=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"_scratch/caloris_v3/caloris_children_v9.csv"),dtype={"id":str})  # g5: v3 (caloris_test_record_v9_v3cat.py)
cal["name"]=cal.name.fillna("").astype(str).str.strip()
assert len(cal)==179
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
        if a0<-180 or a1>180:   # wrap across the antimeridian (as in the v5 script)
            if a1>180: parts=[(a0,180),(-180,a1-360)]
            else: parts=[(a0+360,180),(-180,a1)]
            ws=[int(round(W*(p[1]-p[0])/(2*dlo))) for p in parts]; ws[-1]=W-ws[0]
            a=np.concatenate([_read(ds,p[0],p[1],lat0-dla,lat0+dla,H,wi,bands) for p,wi in zip(parts,ws)],axis=2)
        else: a=_read(ds,a0,a1,lat0-dla,lat0+dla,H,W,bands)
    return a,[-half_km,half_km,-half_km,half_km]
def local(lat,lon,lat0,lon0): return (lon180(lon-lon0))*KD*np.cos(np.radians(lat0)), (lat-lat0)*KD
def gray(a):
    a=a[0].astype(float); p0,p1=np.percentile(a[a>0],[1,99.5]) if (a>0).any() else (0,1); return np.clip((a-p0)/(p1-p0),0,1)
def circle(ax,x,y,r,**kw):
    t=np.linspace(0,2*np.pi,200); ax.plot(x+r*np.cos(t),y+r*np.sin(t),**kw)
OUT=os.path.dirname(os.path.abspath(__file__)); NAME="figS3_caloris_map_v9_v3cat"  # g5: output folder and name
for ext in ("png","jpg","pdf"): assert os.environ.get("REDRAW_V9")=="1" or not os.path.exists(f"{OUT}/{NAME}.{ext}"), f"refusing to overwrite {NAME}.{ext}"   # REDRAW_V9=1 redraws only this script's own _v9 outputs
fig=plt.figure(figsize=(5.9,5.9)); ax=fig.add_axes([0.01,0.01,0.98,0.98])
lat0,lon0=31.93,164.1; half=950                                # map center = catalog center, as in v5
a,ext=crop(BDR,lat0,lon0,half,maxpx=1600,bands=[1]); ax.imshow(gray(a),cmap="gray",extent=ext,origin="upper")
ox,oy=local(29.9087500241,161.431849032,lat0,lon0)            # Orgel et al. (2020) Caloris: 1,556 km
circle(ax,ox,oy,1555.91479183/2,color="w",lw=.9,ls=(0,(4,3)))
circle(ax,0,0,1595.019469/2,color="#e6e6e6",lw=.8,ls=(0,(1,1.6)))  # Herrick et al. (2018) catalog outline
x,y=local(cal.lat.values,cal.lon.values,lat0,lon0)
below=(cal.through==0).values; fld=(cal.n_in>=3).values
ax.scatter(x[below],y[below],s=cal.D[below]*1.2,facecolors="none",edgecolors="#9ad0ff",lw=.7,label=f"excavation < 2.5 km ({int(below.sum())} craters)")
ax.scatter(x[~below&~fld],y[~below&~fld],s=cal.D[~below&~fld]*1.2,facecolors="none",edgecolors="#ffd166",lw=1.1,label=f"excavation ≥ 2.5 km, < 3 hollows ({int((~below&~fld).sum())})")
ax.scatter(x[fld],y[fld],s=cal.D[fld]*1.2,facecolors="#ff6b35",edgecolors="w",lw=.7,label=f"≥ 3 hollows within 0.9 radii ({int(fld.sum())})")
for nm in ["Sander","Nawahi","Apollodorus","Cunningham","Balanchine","Munch","Poe","Atget"]:
    r=cal[cal.name==nm]
    if len(r): xx,yy=local(r.lat.values[0],r.lon.values[0],lat0,lon0); ax.text(xx+20,yy+12,nm,color="w",fontsize=7)
ax.set_xlim(-half,half); ax.set_ylim(-half,half); ax.set_xticks([]); ax.set_yticks([])
h,l=ax.get_legend_handles_labels()
h+=[Line2D([],[],color="w",lw=.9,ls=(0,(4,3))),Line2D([],[],color="#e6e6e6",lw=.8,ls=(0,(1,1.6)))]; l+=["basin (Orgel et al. 2020)","catalog outline (Herrick et al. 2018)"]
ax.legend(h,l,loc="upper left",fontsize=7,frameon=True,facecolor="black",labelcolor="w",edgecolor="none",framealpha=.6,markerscale=.8)
ax.plot([half-325,half-25],[-half+40,-half+40],color="w",lw=2); ax.text(half-175,-half+60,"300 km",color="w",ha="center",fontsize=7)
fig.savefig(f"{OUT}/{NAME}.png",dpi=400); fig.savefig(f"{OUT}/{NAME}.pdf"); plt.close(fig)
Image.open(f"{OUT}/{NAME}.png").convert("RGB").save(f"{OUT}/{NAME}.jpg",quality=90,dpi=(400,400))
print("wrote", len(cal), int(below.sum()), int((~below&~fld).sum()), int(fld.sum()), int((below&fld).sum()), "outside Orgel rim:", int((cal.rho_cal>=1).sum()))
