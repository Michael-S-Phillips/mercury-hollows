import re, glob, os, itertools, sqlite3, numpy as np, pandas as pd
HE=os.path.expanduser("~/mnt/hollows-emissivity"); R=2439.7
c=sqlite3.connect(f"{HE}/GIS/vector/bickel2024_hollows_corrected.gpkg")
P=pd.read_sql("select mdis_id,img_x,img_y,lat,lon180,lat_orig,lon180_orig from hollows_points_corrected",c)
def local(lat,lon,lat0,lon0): return R*np.radians(((np.asarray(lon)-lon0+180)%360)-180)*np.cos(np.radians(lat0)), R*np.radians(np.asarray(lat)-lat0)
def bbox_map(la,lo,u,v,swap,fu,fv):
    a,b=(v,u) if swap else (u,v)
    a=1-a if fu else a; b=1-b if fv else b
    return la.max()-b*(la.max()-la.min()), lo.min()+a*(lo.max()-lo.min())
def best_bbox(la,lo,g):
    lat0,lon0=g.lat_orig.median(),g.lon180_orig.median(); B=np.c_[local(g.lat_orig.values,g.lon180_orig.values,lat0,lon0)]; out=[]
    for sw,fu,fv in itertools.product([0,1],[0,1],[0,1]):
        pl,pn=bbox_map(la,lo,g.img_x.values,g.img_y.values,sw,fu,fv); Q=np.c_[local(pl,pn,lat0,lon0)]
        out.append((np.median(np.hypot(*(Q-B).T)),(sw,fu,fv)))
    return min(out)
print("== label reticle corners (6 local frames) ==")
for f in sorted(glob.glob(f"{HE}/data/mdis_frames/IMG/CN*_IF_5.IMG")):
    im=os.path.basename(f)[:13]; g=P[P.mdis_id==im]
    if len(g)<3: continue
    t=open(f,'rb').read(4096*3).decode('latin1')
    vec=lambda k: np.array([float(x) for x in re.findall(r'(-?\d+\.\d+)\s*<DEG>',re.search(k+r'\s*=\s*\(([^)]*)\)',t).group(1))])
    la=vec('RETICLE_POINT_LATITUDE'); lo=((vec('RETICLE_POINT_LONGITUDE')+180)%360)-180
    e,k=best_bbox(la,lo,g); print(f"  {im}: bbox-of-reticles mapping vs released, median {e*1000:.0f} m  (swap,flipx,flipy)={k}")
print("== camera-model corners (all images with >=6 detections) ==")
def H_fit(u,v,x,y):
    A=[]
    for uu,vv,xx,yy in zip(u,v,x,y): A+=[[uu,vv,1,0,0,0,-xx*uu,-xx*vv,-xx],[0,0,0,uu,vv,1,-yy*uu,-yy*vv,-yy]]
    return np.linalg.svd(np.array(A))[2][-1].reshape(3,3)
res=[]; ks=[]
for im,g in P.groupby("mdis_id"):
    if len(g)<6: continue
    lat0,lon0=g.lat.median(),g.lon180.median(); x,y=local(g.lat.values,g.lon180.values,lat0,lon0); H=H_fit(g.img_x.values,g.img_y.values,x,y)
    q=H@np.array([[0,1,0,1],[0,0,1,1],[1,1,1,1.]]); X,Y=q[:2]/q[2]
    la=lat0+np.degrees(Y/R); lo=lon0+np.degrees(X/(R*np.cos(np.radians(lat0))))
    e,k=best_bbox(la,lo,g); res.append(e); ks.append(k)
res=np.array(res); print("  images:",len(res)," median error (km) quantiles 25/50/75/90:", np.round(np.percentile(res,[25,50,75,90]),3))
from collections import Counter; print("  orientation counts:", Counter(ks).most_common(4))
print("  fraction with median error < 0.3 km:", (res<0.3).mean().round(3), " < 1 km:", (res<1).mean().round(3))
