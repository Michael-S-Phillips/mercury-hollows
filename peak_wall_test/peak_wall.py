import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Peak vs wall vs floor hollow density inside hollow-hosting craters (D >= 30 km), from the 665 m USGS DEM.
Per host crater: terrain classes from radial position and height above the floor; areas weighted by catalog completeness for 1 km
hollows and with superposed Herrick craters (D >= 5 km) masked, because hollows inside them are assigned to them; hollows counted
only if assigned to this host. Writes per-crater rows. Usage: python3 peak_wall.py START STOP"""
import sys, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
B="/Users/x"  # placeholder
import os
HE=os.path.expanduser(_REPO + "")
DEM=f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM=2439.4; KD=np.pi*RM/180
cr=pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"]=cr.lon%360
ha=pd.read_csv(f"{HE}/code/rerun/latitude/aspect/host_assignments.csv")
hc=pd.read_csv(f"{HE}/code/rerun/out_v2/hollows_canonical.csv")[["hollow_id","field_id"]]
ha=ha.merge(hc,on="hollow_id",how="left",suffixes=("","_c"))
C=np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz"); Cg=C["C"][list(C["Dref"]).index(1.0)]
org=pd.read_csv(f"{HE}/peak_wall_test/orgel_basins.csv"); org=org[org.Confidence.isin(["Certain","Probable"])]
hosts=cr[cr.unique_id.isin(set(ha.host_uid.dropna().astype(int)))&(cr.diam>=30)].sort_values("unique_id").reset_index(drop=True)
a0,a1=int(sys.argv[1]),int(sys.argv[2]); HT=float(sys.argv[3]) if len(sys.argv)>3 else 250.0; RPK=float(sys.argv[4]) if len(sys.argv)>4 else 0.45
def gcd(la1,lo1,la2,lo2):
    p1,p2=np.radians(la1),np.radians(la2); return RM*np.arccos(np.clip(np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2-lo1)),-1,1))
rows=[]
with rasterio.open(DEM) as ds:
  for _,c in hosts.iloc[a0:a1].iterrows():
    R=c.diam/2; la0,lo0=c.lat,c.lon360; cl=np.cos(np.radians(la0))
    dla=1.25*R/KD; dlo=min(1.25*R/(KD*max(cl,0.05)),179)
    xs=lambda lo: np.radians(lo-180)*2439400.0; ys=lambda la: np.radians(la)*2439400.0
    lo_a,lo_b=lo0-dlo,lo0+dlo
    if lo_a<0 or lo_b>360 or la0+dla>89 or la0-dla<-89: continue   # skip seam / pole crossings
    win=from_bounds(xs(lo_a),ys(la0-dla),xs(lo_b),ys(la0+dla),ds.transform)
    z=ds.read(1,window=win,boundless=True,fill_value=-32768).astype(float); z[z==-32768]=np.nan; z*=0.5
    T=ds.window_transform(win); nr,nc=z.shape
    xc=T.c+T.a*(np.arange(nc)+.5); yc=T.f+T.e*(np.arange(nr)+.5)
    LON=180+np.degrees(xc/2439400.0)[None,:]*np.ones((nr,1)); LAT=np.degrees(yc/2439400.0)[:,None]*np.ones((1,nc))
    dx=(LON-lo0)*KD*cl; dy=(LAT-la0)*KD; rho=np.hypot(dx,dy)/R
    pole=np.sign(la0 if la0!=0 else 1)*dy>0          # poleward side of center: wall faces the equator
    ring=c.cent_str in ("pr","rpc","pb")
    fz=z[(rho>0.15)&(rho<0.6)]; fz=fz[np.isfinite(fz)]
    if len(fz)<20: continue
    floor=np.percentile(fz,20); h=z-floor
    H=HT
    peak=(rho<=(0.6 if ring else RPK))&(h>H)
    wall=(rho>=(0.7 if ring else 0.55))&(rho<=1.05)&(h>H)
    flo=(rho<=0.8)&(h<=H)
    # mask superposed craters (D>=5 km, smaller than half the host, center inside the host rim)
    sub=cr[(cr.diam>=5)&(cr.diam<0.5*c.diam)&(np.abs(cr.lat-la0)<dla)]
    sub=sub[gcd(la0,lo0,sub.lat.values,sub.lon360.values)<1.05*R]
    m=np.zeros_like(rho,bool)
    for s in sub.itertuples():
        m|=np.hypot((LON-s.lon360)*KD*np.cos(np.radians(s.lat)),(LAT-s.lat)*KD)<1.1*s.diam/2
    valid=np.isfinite(z)&~m&(rho<=1.05)
    ci=np.clip(((90-LAT)*10).astype(int),0,1799); cj=np.clip((LON*10).astype(int),0,3599)
    w=Cg[ci,cj]; w=np.where(np.isfinite(w),w,0)
    pxa=(0.66524**2)*np.cos(np.radians(LAT))
    A=lambda msk: float(np.sum(pxa*w*(msk&valid)))
    hs=ha[ha.host_uid==c.unique_id]
    hi=np.clip(np.round((T.f-ys(hs.lat.values))/(-T.e)-.5).astype(int),0,nr-1)
    hj=np.clip(np.round((xs(hs.lon.values%360)-T.c)/T.a-.5).astype(int),0,nc-1)
    cls=np.where(peak[hi,hj],"peak",np.where(wall[hi,hj],"wall",np.where(flo[hi,hj],"floor","other")))
    cls=np.where(m[hi,hj],"masked",cls)
    eqw=wall&pole; pw=wall&~pole
    hpole=pole[hi,hj]
    gy,gx=np.gradient(z,665.24,665.24); slope=np.degrees(np.arctan(np.hypot(gx,gy)))
    sb=np.digitize(slope,[5,15]); hsb=sb[hi,hj]
    extra={}
    fin=flo&(rho<0.45); fout=flo&(rho>=0.45)
    extra["A_floor_in"]=A(fin); extra["A_floor_out"]=A(fout)
    extra["n_floor_in"]=int(((cls=="floor")&fin[hi,hj]).sum()); extra["n_floor_out"]=int(((cls=="floor")&fout[hi,hj]).sum())
    for k,(nm,msk) in enumerate([("peak",peak),("wall",wall)]):
        for b in range(3):
            extra[f"A_{nm}_s{b}"]=A(msk&(sb==b)); extra[f"n_{nm}_s{b}"]=int(((cls==nm)&(hsb==b)).sum())
    # parents: Herrick craters >= 2D containing the center, and Orgel certain/probable basins
    par=cr[(cr.diam>=2*c.diam)&(np.abs(cr.lat-la0)<20)]; par=par[gcd(la0,lo0,par.lat.values,par.lon360.values)<par.diam.values/2]
    ob=org[gcd(la0,lo0,org.y_coord.values.astype(float),org.x_coord.values.astype(float)%360)<org.Diam_km.values.astype(float)/2]
    rows.append(dict(uid=int(c.unique_id),name=c["name"] if isinstance(c["name"],str) else "",D=c.diam,lat=la0,lon=lo0,cent_str=str(c.cent_str).strip(),filled=c.filled,preserved=c.preserved,
        floor_m=floor,peak_relief_m=float(np.nanmax(np.where(rho<=0.45,h,np.nan))),n_sub_masked=len(sub),
        A_peak=A(peak),A_wall=A(wall),A_wall_eq=A(eqw),A_wall_pole=A(pw),A_floor=A(flo),
        n_peak=int((cls=="peak").sum()),n_wall=int((cls=="wall").sum()),n_wall_eq=int(((cls=="wall")&hpole).sum()),n_wall_pole=int(((cls=="wall")&~hpole).sum()),
        n_floor=int((cls=="floor").sum()),n_other=int((cls=="other").sum()),n_masked=int((cls=="masked").sum()),
        f_peak=int(hs[cls=="peak"].field_id.nunique()),f_wall=int(hs[cls=="wall"].field_id.nunique()),f_floor=int(hs[cls=="floor"].field_id.nunique()),
        parent_D=float(par.diam.max()) if len(par) else 0.0,parent_name=(par.sort_values("diam").iloc[-1]["name"] if len(par) and isinstance(par.sort_values("diam").iloc[-1]["name"],str) else ""),
        in_orgel_basin=len(ob)>0,orgel=";".join(ob.Basin_Name.astype(str)),**extra))
out=f"{HE}/peak_wall_test/rows_H{int(HT)}_r{RPK}.csv"; pd.DataFrame(rows).to_csv(out,index=False); print(len(hosts),"hosts total;",len(rows),"rows written",out)
