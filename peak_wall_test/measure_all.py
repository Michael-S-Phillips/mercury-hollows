import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Per-crater DEM morphometry and terrain-resolved hollow counts for all Herrick craters D >= 15 km (hosts and non-hosts).
depth = median over 36 azimuth sectors of the max elevation at 0.85-1.2 R, minus the 20th-percentile elevation at 0.15-0.6 R.
Terrain classes as in peak_wall.py (250 m relief threshold). Usage: python3 measure_all.py START STOP"""
import sys, os, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
HE=os.path.expanduser(_REPO + ""); DEM=f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM=2439.4; KD=np.pi*RM/180; H=250.0
cr=pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"]=cr.lon%360
ha=pd.read_csv(f"{HE}/code/rerun/latitude/aspect/host_assignments.csv")
C=np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz"); Cg=C["C"][list(C["Dref"]).index(1.0)]
cov=np.load(f"{HE}/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"]
DMIN=float(sys.argv[3]) if len(sys.argv)>3 else 15.0; DMAX=float(sys.argv[4]) if len(sys.argv)>4 else 1e9
todo=cr[(cr.diam>=DMIN)&(cr.diam<DMAX)].sort_values("unique_id").reset_index(drop=True)
a0,a1=int(sys.argv[1]),int(sys.argv[2])
def gcd(la1,lo1,la2,lo2):
    p1,p2=np.radians(la1),np.radians(la2); return RM*np.arccos(np.clip(np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2-lo1)),-1,1))
xs=lambda lo: np.radians(lo-180)*2439400.0; ys=lambda la: np.radians(la)*2439400.0
rows=[]
with rasterio.open(DEM) as ds:
  for _,c in todo.iloc[a0:a1].iterrows():
    R=c.diam/2; la0,lo0=c.lat,c.lon360; cl=np.cos(np.radians(la0))
    dla=1.25*R/KD; dlo=1.25*R/(KD*max(cl,0.05))
    if lo0-dlo<0 or lo0+dlo>360 or abs(la0)+dla>88: continue
    step=max(1,int(R/150))      # decimate large basins
    win=from_bounds(xs(lo0-dlo),ys(la0-dla),xs(lo0+dlo),ys(la0+dla),ds.transform)
    z=ds.read(1,window=win,boundless=True,fill_value=-32768)[::step,::step].astype(float); z[z==-32768]=np.nan; z*=0.5
    T=ds.window_transform(win); nr,nc=z.shape
    xc=T.c+T.a*step*(np.arange(nc)+.5); yc=T.f+T.e*step*(np.arange(nr)+.5)
    LON=180+np.degrees(xc/2439400.0)[None,:]*np.ones((nr,1)); LAT=np.degrees(yc/2439400.0)[:,None]*np.ones((1,nc))
    dx=(LON-lo0)*KD*cl; dy=(LAT-la0)*KD; rho=np.hypot(dx,dy)/R; az=np.degrees(np.arctan2(dx,dy))%360
    fz=z[(rho>0.15)&(rho<0.6)]; fz=fz[np.isfinite(fz)]
    if len(fz)<10: continue
    floor=np.percentile(fz,20)
    rimz=[]
    for k in range(36):
        s=z[(rho>0.85)&(rho<1.2)&(az>=k*10)&(az<k*10+10)]; s=s[np.isfinite(s)]
        if len(s): rimz.append(s.max())
    rim=np.median(rimz) if len(rimz)>=18 else np.nan
    sur=z[(rho>1.3)&(rho<1.6)] if False else np.array([])
    h=z-floor; ring=str(c.cent_str).strip() in ("pr","rpc","pb")
    peak=(rho<=(0.6 if ring else 0.45))&(h>H); wall=(rho>=(0.7 if ring else 0.55))&(rho<=1.05)&(h>H); flo=(rho<=0.8)&(h<=H)
    sub=cr[(cr.diam>=5)&(cr.diam<0.5*c.diam)&(np.abs(cr.lat-la0)<dla)]; sub=sub[gcd(la0,lo0,sub.lat.values,sub.lon360.values)<1.05*R]
    m=np.zeros_like(rho,bool)
    for s in sub.itertuples(): m|=np.hypot((LON-s.lon360)*KD*np.cos(np.radians(s.lat)),(LAT-s.lat)*KD)<1.1*s.diam/2
    valid=np.isfinite(z)&~m&(rho<=1.05)
    ci=np.clip(((90-LAT)*10).astype(int),0,1799); cj=np.clip((LON*10).astype(int),0,3599)
    w=np.nan_to_num(Cg[ci,cj]); px=cov[ci,cj]
    pxa=(0.66524*step)**2*np.cos(np.radians(LAT)); A=lambda msk: float(np.sum(pxa*w*(msk&valid)))
    hs=ha[ha.host_uid==c.unique_id]
    if len(hs):
        hi=np.clip(np.round((T.f-ys(hs.lat.values))/(-T.e*step)-.5).astype(int),0,nr-1); hj=np.clip(np.round((xs(hs.lon.values%360)-T.c)/(T.a*step)-.5).astype(int),0,nc-1)
        cls=np.where(peak[hi,hj],"peak",np.where(wall[hi,hj],"wall",np.where(flo[hi,hj],"floor","other")))
    else: cls=np.array([])
    inside=valid&(rho<=1.0)
    rows.append(dict(uid=int(c.unique_id),name=c["name"] if isinstance(c["name"],str) else "",D=c.diam,lat=la0,lon=lo0,cent_str=str(c.cent_str).strip(),filled=c.filled,preserved=c.preserved,
       rim=rim,floor=floor,depth=rim-floor,A_peak=A(peak),A_wall=A(wall),A_floor=A(flo),A_all=A(inside),C_mean=float(np.mean(w[inside])) if inside.any() else 0.0,
       frac_covered_100m=float(np.mean(np.isfinite(px[inside])&(px[inside]<=100))) if inside.any() else 0.0,
       n_all=len(hs),n_peak=int((cls=="peak").sum()),n_wall=int((cls=="wall").sum()),n_floor=int((cls=="floor").sum())))
TAG=("_"+sys.argv[3]) if len(sys.argv)>3 else ""
pd.DataFrame(rows).to_csv(f"{HE}/peak_wall_test/morph{TAG}_{a0:05d}_{a1:05d}.csv",index=False); print(len(rows))
