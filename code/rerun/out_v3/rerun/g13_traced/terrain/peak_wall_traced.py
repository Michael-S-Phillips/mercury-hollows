# g13 revision (7 Nov 2026): env G13_REV (r1 | r2, default r2): r2 reads traces/traces_clean_r2.json (Michael's refined wall feet and slump outlines, 7 Nov) and counts slump as wall (slump inside the wall foot becomes wall, peak still wins); r1 reproduces the 6 Nov results.
# g13 copy of g11_curvature/terrain/peak_wall_curv.py (6 Nov 2026). Changed: for hosts Michael traced (g13_traced/traces/traces_clean.json, status done) the classes come from his outlines (peak = peak outlines minus floor-inside-ring outlines; wall = outside the wall foot out to 1.05 R, not peak; floor = inside the wall foot, not peak); every other host uses curvature classes (env CLS: v3 = g10 curv_classes_v3 [default], v2 = g10 curv_classes as in g11). Columns class_source, A/n by source unchanged otherwise. Also writes out_<CAT>/_chunks/hollow_classes_<CLS>_START_STOP.csv with each hollow's class under the primary rule, curvature v3 and curvature v2, and out_<CAT>/_chunks/area_overlap_<CLS>_START_STOP.csv with weighted peak/wall/floor areas under each rule for traced hosts. Env NOTRACE=1 ignores the traces (curvature classes everywhere; rows_curvonly_<CLS>.csv). Usage: peak_wall_traced.py START STOP (curvature defaults nsec 24, hmin 100 m, sig_frac 0.03).
# g11 copy of g10_threshold/peak_wall_scaled.py (6 Nov 2026). Changed: classes from curv_classes.classify_curv (g10_threshold, imported, not modified) in place of terrain_classes.classify with a relief threshold; usage peak_wall_curv.py START STOP [NSEC HMIN SIG_FRAC] (defaults 24 100 0.03); hosts limited to D <= 250 km (rows above 250 km were dropped by analyze.py anyway); output g11_curvature/terrain/out_<CAT>/_chunks/rows_curv_n<NSEC>_h<HMIN>_s<SIG>_START_STOP.csv; floor_m and peak_relief_m kept from the old floor level (20th percentile at 0.15-0.6 R) for reference only; n_moat_sectors replaced by n_ring_sectors; HT_used replaced by nsec, hmin, sig_frac. Everything else (superposed-crater mask, completeness weights, hollow counting, slope bins, equator/pole walls, parents, Orgel basins) unchanged.
# g10 copy of g7_terrain_fix/peak_wall_fix.py (6 Nov 2026). Changed: output folder g10_threshold/out_<CAT>_<G10_MODE>; with G10_MODE=depth<f>
# the relief threshold of each crater is f x its depth (median over 36 sectors of the max elevation at 0.85-1.2 R minus the floor level),
# in place of the fixed HT (at least 100 m); columns depth_m and HT_used added. terrain_classes is imported from g7_terrain_fix.
# g7 copy of g3_catalog/terrain/peak_wall.py (itself a copy of peak_wall_test/peak_wall.py), 6 Nov 2026. Changed: the class rules
# (terrain_classes.classify, which closes the unclassified gaps), output folder g7_terrain_fix/out_<CAT>, and env G7_GRID
# (v2 = out_v2 completeness grid, v3 = g4_completeness/completeness_10ppd_v3.npz).
"""Peak vs wall vs floor hollow density inside hollow-hosting craters (D >= 30 km), from the 665 m USGS DEM.
Per host crater: terrain classes from radial position and height above the floor; areas weighted by catalog completeness for 1 km
hollows and with superposed Herrick craters (D >= 5 km) masked, because hollows inside them are assigned to them; hollows counted
only if assigned to this host. Writes per-crater rows. Usage: python3 peak_wall.py START STOP"""
import sys, numpy as np, pandas as pd, rasterio
sys.path.insert(0, __import__("os").path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g10_threshold")); from curv_classes import classify_curv as classify_v2; from curv_classes_v3 import classify_curv as classify_v3
import json as _json; from matplotlib.path import Path as _Path
from rasterio.windows import from_bounds
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced/terrain"); CLS = _os.environ.get("CLS", "v3")
REV = _os.environ.get("G13_REV", "r2"); TR = {int(k): v for k, v in _json.load(open(_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced/traces/" + ("traces_clean.json" if REV == "r1" else "traces_clean_r2.json"))))["craters"].items() if v["status"] == "done" and not _os.environ.get("NOTRACE")}
GRID = _os.environ.get("G7_GRID", "v3")
CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/out_{CAT}"; _os.makedirs(ODIR, exist_ok=True)
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))
B="/Users/x"  # placeholder
import os
HE=os.path.expanduser("~/mnt/hollows-emissivity")
DEM=f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM=2439.4; KD=np.pi*RM/180
cr=pd.read_csv(f"{HE}/data/ME6M300TGT/herrick2018_all.csv"); cr["lon360"]=cr.lon%360
ha=pd.read_csv(HOSTS)
hc=pd.read_csv(CANON)[["hollow_id","field_id"]]
ha=ha.merge(hc,on="hollow_id",how="left",suffixes=("","_c"))
C=np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz" if GRID=="v2" else f"{HE}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_{CAT if CAT in ('v3','v3b') else 'v3'}.npz"); Cg=C["C"][list(C["Dref"]).index(1.0)]
org=pd.read_csv(f"{HE}/peak_wall_test/orgel_basins.csv"); org=org[org.Confidence.isin(["Certain","Probable"])]
hosts=cr[cr.unique_id.isin(set(ha.host_uid.dropna().astype(int)))&(cr.diam>=30)&(cr.diam<=250)].sort_values("unique_id").reset_index(drop=True)
a0,a1=int(sys.argv[1]),int(sys.argv[2]); NSEC=int(sys.argv[3]) if len(sys.argv)>3 else 24; HMIN=float(sys.argv[4]) if len(sys.argv)>4 else 100.0; SIG=float(sys.argv[5]) if len(sys.argv)>5 else 0.03
CFG=(f"curvonly_{CLS}" if _os.environ.get("NOTRACE") else (f"traced_{CLS}" if REV == "r1" else f"traced_{REV}_{CLS}"))
def aeqd(lon,lat,la0,lo0):
    p,p0,dl=np.radians(lat),np.radians(la0),np.radians(lon-lo0); cc=np.clip(np.sin(p0)*np.sin(p)+np.cos(p0)*np.cos(p)*np.cos(dl),-1,1); c=np.arccos(cc)
    k=np.where(c>1e-12,c/np.sin(np.maximum(c,1e-12)),1.0); return RM*k*np.cos(p)*np.sin(dl), RM*k*(np.cos(p0)*np.sin(p)-np.sin(p0)*np.cos(p)*np.cos(dl))
def inside(polys,X,Y,t):
    pts=np.column_stack([X.ravel(),Y.ravel()]); m=np.zeros(len(pts),bool)
    for rings in polys:
        rings=[r for r in rings if len(r)>=3]
        if not rings: continue
        P=[np.column_stack(aeqd(np.array(r)[:,0],np.array(r)[:,1],t["lat0"],t["lon0"])) for r in rings]
        q=_Path(P[0]).contains_points(pts)
        for h in P[1:]: q&=~_Path(h).contains_points(pts)
        m|=q
    return m.reshape(X.shape)
def traced_classes(t,LON,LAT,rho,z):
    X,Y=aeqd(LON,LAT,t["lat0"],t["lon0"]); pk=inside(t["peak"],X,Y,t)&np.isfinite(z); ft=inside(t["foot"],X,Y,t)
    sl=inside(t.get("slump",[]),X,Y,t) if REV!="r1" else np.zeros_like(ft)
    wall=((~ft)|sl)&(rho<=1.05)&~pk&np.isfinite(z); flo=ft&~sl&~pk&np.isfinite(z); return dict(peak=pk,wall=wall,flo=flo,n_ring_sectors=-1)
hrows=[]; arows=[]
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
    az=np.degrees(np.arctan2(dy,dx))%360
    F0_=np.percentile(fz,20); rim_=[np.nanmax(np.where((rho>=0.85)&(rho<=1.2)&(az>=k*10)&(az<k*10+10),z,np.nan)) for k in range(36)]
    depth_=float(np.nanmedian(rim_)-F0_)
    cent=str(c.cent_str).strip() if isinstance(c.cent_str,str) else ""   # g11: Herrick cent_str stripped ("" if missing)
    K3=classify_v3(z,rho,az,R,cent,nsec=NSEC,hmin=HMIN,sig_frac=SIG); K2=classify_v2(z,rho,az,R,cent,nsec=NSEC,hmin=HMIN,sig_frac=SIG)
    uid_=int(c.unique_id); src_="traced" if uid_ in TR else f"curv_{CLS}"
    K=traced_classes(TR[uid_],LON,LAT,rho,z) if uid_ in TR else (K3 if CLS=="v3" else K2)   # g13: traced outlines where available
    floor=F0_; h=z-F0_; peak,wall,flo=K["peak"],K["wall"],K["flo"]
    nmoat=K["n_ring_sectors"]
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
    def _cls(KK):
        q=np.where(KK["peak"][hi,hj],"peak",np.where(KK["wall"][hi,hj],"wall",np.where(KK["flo"][hi,hj],"floor","other"))); return np.where(m[hi,hj],"masked",q)
    c3,c2=_cls(K3),_cls(K2)
    for _i in range(len(hs)): hrows.append(dict(hollow_id=hs.hollow_id.values[_i],host_uid=uid_,source=src_,cls=cls[_i],cls_v3=c3[_i],cls_v2=c2[_i]))
    if uid_ in TR:
        for a_ in ("peak","wall","flo"):
            for b_ in ("peak","wall","flo"): arows.append(dict(uid=uid_,traced=a_,other="v3",cls=b_,A=A(K[a_]&K3[b_]))); arows.append(dict(uid=uid_,traced=a_,other="v2",cls=b_,A=A(K[a_]&K2[b_])))
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
        in_orgel_basin=len(ob)>0,n_ring_sectors=nmoat,class_source=src_,depth_m=depth_,nsec=NSEC,hmin=HMIN,sig_frac=SIG,orgel=";".join(ob.Basin_Name.astype(str)),**extra))
_os.makedirs(f"{ODIR}/_chunks", exist_ok=True); out=f"{ODIR}/_chunks/rows_{CFG}_{a0:04d}_{a1:04d}.csv"; pd.DataFrame(rows).to_csv(out,index=False); pd.DataFrame(hrows).to_csv(f"{ODIR}/_chunks/hollow_classes_{CFG}_{a0:04d}_{a1:04d}.csv",index=False); pd.DataFrame(arows).to_csv(f"{ODIR}/_chunks/area_overlap_{CFG}_{a0:04d}_{a1:04d}.csv",index=False); print(len(hosts),"hosts total;",len(rows),"rows written",out)
