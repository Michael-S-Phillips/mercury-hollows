import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
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
cal=pd.read_csv("caloris_children.csv"); m=pd.read_csv("morph_all_nested.csv"); ch=pd.read_csv("nested_children_classified.csv")
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
fig=plt.figure(figsize=(7.2,3.6)); ax=fig.add_axes([0.0,0.06,0.5,0.88]); bx=fig.add_axes([0.6,0.15,0.38,0.78])
lat0,lon0=31.93,164.1; half=950
a,ext=crop(BDR,lat0,lon0,half,maxpx=1100,bands=[1]); ax.imshow(gray(a),cmap="gray",extent=ext,origin="upper")
circle(ax,0,0,1595/2,color="w",lw=.8,ls=(0,(4,3)))
x,y=local(cal.lat.values,cal.lon.values,lat0,lon0)
below=cal.through==0; fld=cal.n_in>=3
ax.scatter(x[below],y[below],s=cal.D[below]*0.6,facecolors="none",edgecolors="#9ad0ff",lw=.6,label="excavates < 2.5 km (no hollow fields)")
ax.scatter(x[~below&~fld],y[~below&~fld],s=cal.D[~below&~fld]*0.6,facecolors="none",edgecolors="#ffd166",lw=1.0,label="excavates ≥ 2.5 km, no hollow field")
ax.scatter(x[fld],y[fld],s=cal.D[fld]*0.6,facecolors="#ff6b35",edgecolors="w",lw=.6,label="hollow field (≥ 3 hollows)")
for nm in ["Sander","Nawahi","Apollodorus","Cunningham","Balanchine","Munch","Poe","Atget"]:
    r=cal[cal.name==nm]
    if len(r): xx,yy=local(r.lat.values[0],r.lon.values[0],lat0,lon0); ax.text(xx+18,yy+10,nm,color="w",fontsize=6)
ax.set_xlim(-half,half); ax.set_ylim(-half,half); ax.set_xticks([]); ax.set_yticks([])
ax.legend(loc="lower left",fontsize=5.8,frameon=True,facecolor="black",labelcolor="w",edgecolor="none",framealpha=.6,markerscale=.8)
ax.text(-half+25,half-60,"a",color="w",fontsize=9,fontweight="bold")
ax.plot([half-325,half-25],[-half+40,-half+40],color="w",lw=2); ax.text(half-175,-half+60,"300 km",color="w",ha="center",fontsize=6.5)
iso=m[(m.parent<0)&~m.in_orgel&(m.C_mean>=0.2)]
bins=np.array([10,15,20,25,30,35,40,50,70,130]); xc=np.sqrt(bins[:-1]*bins[1:])
def frac(df):
    out=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        s=df[(df.D>=lo)&(df.D<hi)]; n=len(s); k=int(s.field.sum()) if "field" in s else int((s.n_in>=3).sum())
        from scipy.stats import beta
        lo_,hi_=(beta.ppf(.025,k,n-k+1) if k>0 else 0, beta.ppf(.975,k+1,n-k) if n>k else 1) if n else (np.nan,np.nan)
        out.append((k/n if n else np.nan,lo_,hi_,n,k))
    return np.array(out)
calf=cal.copy(); calf["field"]=(calf.n_in>=3).astype(int)
for df,c,lab,mk,off in [(iso,"#6b6b67","isolated craters","s",0.97),(calf,"#ff6b35","inside Caloris","o",1.03)]:
    f=frac(df); ok=f[:,3]>0
    bx.errorbar(xc[ok]*off,f[ok,0],yerr=[f[ok,0]-f[ok,1],f[ok,2]-f[ok,0]],fmt=mk,color=c,ms=4,lw=.8,capsize=1.5,label=lab)
    if c=="#ff6b35":
        for xx,(fr,lo_,hi_,n,k) in zip(xc[ok]*off,f[ok]): bx.text(xx*1.04,hi_+0.02,f"{int(k)}/{int(n)}",color=c,fontsize=5.4,ha="center",rotation=90,va="bottom")
d1,d2=(2.5/(0.1*12**0.15))**(1/0.85),(3.5/(0.1*12**0.15))**(1/0.85)
bx.axvspan(d1,d2,color="#7B4EA8",alpha=.15,lw=0); bx.text(np.sqrt(d1*d2),1.24,"excavation depth =\nplains thickness\n(2.5–3.5 km)",ha="center",va="top",fontsize=5.4,color="#5a3a80")
bx.set_xscale("log"); bx.set_xlim(9,140); bx.set_ylim(-0.02,1.25); bx.set_xticks([10,20,30,50,100]); bx.set_xticklabels(["10","20","30","50","100"])
bx.set_xlabel("crater diameter (km)"); bx.set_ylabel("fraction hosting a hollow field"); bx.legend(fontsize=6,frameon=False,loc="upper left",bbox_to_anchor=(0,0.9))
for sp in ("top","right"): bx.spines[sp].set_visible(False)
bx.text(-0.16,1.02,"b",transform=bx.transAxes,fontsize=9,fontweight="bold")
fig.savefig("fig_nested_A_caloris.png",dpi=220); plt.close(fig)
# ---------------- example panels
def example_row(axg,axc,uid,title,half_mult=1.6,hollow_all=True):
    r=m[m.uid==uid].iloc[0]; half=max(12,half_mult*r.D/2)
    a,ext=crop(BDR,r.lat,r.lon,half,maxpx=700,bands=[1]); axg.imshow(gray(a),cmap="gray",extent=ext)
    c,_=crop(COL,r.lat,r.lon,half,maxpx=400); axc.imshow(rgb(c),extent=ext)
    hh=hc[(np.abs(hc.lat-r.lat)<half/KD)&(np.abs(lon180(hc.lon-r.lon))<half/(KD*np.cos(np.radians(r.lat))))]
    hx,hy=local(hh.lat.values,hh.lon.values,r.lat,r.lon)
    for ax in (axg,axc):
        circle(ax,0,0,r.D/2,color="#9ad0ff",lw=.7,ls=(0,(3,2))); ax.set_xlim(-half,half); ax.set_ylim(-half,half); ax.set_xticks([]); ax.set_yticks([])
    axg.scatter(hx,hy,s=5,c="#ff6b35",edgecolors="none")
    L=10 if half>20 else 5; axg.plot([half*0.9-L,half*0.9],[-half*0.9]*2,color="w",lw=2); axg.text(half*0.9-L/2,-half*0.82,f"{L} km",color="w",ha="center",fontsize=6)
    axg.set_title(title,fontsize=6.2,loc="left")
    return r
cases=[(112324,"a  {D:.0f} km, excavates ~{dx:.1f} km\nwithin the plains: no hollows"),
       (110684,"b  Cunningham, {D:.0f} km, ~{dx:.1f} km\njust through: {n} hollows"),
       (109122,"c  Sander, {D:.0f} km, ~{dx:.1f} km\nwell through: {n} hollows"),
       (110242,"d  {D:.0f} km, ~{dx:.1f} km, degraded\nthrough, well imaged: no hollows")]
fig,axs=plt.subplots(2,4,figsize=(7.2,3.9))
for j,(uid,t) in enumerate(cases):
    r=m[m.uid==uid].iloc[0]; dx=0.1*12**0.15*r.D**0.85
    example_row(axs[0,j],axs[1,j],uid,t.format(D=r.D,dx=dx,n=int(r.n_in)))
axs[0,0].text(-0.08,0.5,"166 m/px mosaic\n+ hollows (orange)",transform=axs[0,0].transAxes,rotation=90,va="center",ha="right",fontsize=6.5)
axs[1,0].text(-0.08,0.5,"enhanced color\n(LRM is dark blue)",transform=axs[1,0].transAxes,rotation=90,va="center",ha="right",fontsize=6.5)
fig.subplots_adjust(left=0.05,right=0.995,top=0.9,bottom=0.01,wspace=0.04,hspace=0.04); fig.savefig("fig_nested_B_caloris_examples.png",dpi=220); plt.close(fig)
# ---------------- Mendelssohn: one parent, two children
p=m[m.uid==200199].iloc[0]; A=m[m.uid==112519].iloc[0]; B=m[m.uid==114342].iloc[0]
fig=plt.figure(figsize=(7.2,3.5))
ax=fig.add_axes([0.0,0.02,0.36,0.9]); half=p.D/2*1.08
a,ext=crop(BDR,p.lat,p.lon,half,maxpx=900,bands=[1]); ax.imshow(gray(a),cmap="gray",extent=ext); circle(ax,0,0,p.D/2,color="w",lw=.7,ls=(0,(4,3)))
for r,c in [(A,"#ff6b35"),(B,"#9ad0ff")]:
    xx,yy=local(r.lat,r.lon,p.lat,p.lon); circle(ax,xx,yy,r.D/2,color=c,lw=1.0)
hh=hc[(np.abs(hc.lat-p.lat)<half/KD)&(np.abs(lon180(hc.lon-p.lon))<half/(KD*np.cos(np.radians(p.lat))))]; hx,hy=local(hh.lat.values,hh.lon.values,p.lat,p.lon)
ax.scatter(hx,hy,s=4,c="#ff6b35",edgecolors="none"); ax.set_xlim(-half,half); ax.set_ylim(-half,half); ax.set_xticks([]); ax.set_yticks([])
tf=0.0
ax.set_title(f"Mendelssohn ({p.D:.0f} km, {p.lat:.0f}°N): floor fill ~{ch[ch.uid==A.uid].iloc[0].tf:.1f} km (depth deficit)",fontsize=6.6,loc="left")
ax.plot([half*0.9-50,half*0.9],[-half*0.9]*2,color="w",lw=2); ax.text(half*0.9-25,-half*0.83,"50 km",color="w",ha="center",fontsize=6)
spots=[(A,"#ff6b35",f"{A.D:.0f} km child: excavates ~{0.1*12**0.15*A.D**0.85:.1f} km,\nthrough the fill: {int(A.n_in)} hollows"),(B,"#9ad0ff",f"{B.D:.0f} km child: excavates ~{0.1*12**0.15*B.D**0.85:.1f} km,\nwithin the fill: no hollows")]
for k,(r,c,t) in enumerate(spots):
    axg=fig.add_axes([0.38+k*0.31,0.47,0.3,0.42]); axc=fig.add_axes([0.38+k*0.31,0.02,0.3,0.42])
    example_row(axg,axc,int(r.uid),t,half_mult=1.7)
fig.savefig("fig_nested_C_mendelssohn.png",dpi=220); plt.close(fig)
print("figs ok")
