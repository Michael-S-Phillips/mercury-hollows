"""Which way the release orients each image, and what a correct use of the same label corners would give.
(1) For 933 images (>= 6 detections), predict the axis directions from the camera-model corners: column 0 is west unless the
left corners lie east of the right corners; row 0 is north unless the top corners lie south of the bottom corners. Compare
with the orientation that best reproduces the released coordinates (bbox_test.py logic).
(2) For the 6 frames with PDS labels, interpolate the label corner points bilinearly in their own order (corner 1..4 at
image (0,0),(1,0),(0,1),(1,1)) and compare with the camera-model (recomputed) coordinates."""
import re, glob, os, itertools, sqlite3, numpy as np, pandas as pd
HE=os.path.expanduser("~/mnt/hollows-emissivity"); R=2439.7
c=sqlite3.connect(f"{HE}/GIS/vector/bickel2024_hollows_corrected.gpkg")
P=pd.read_sql("select mdis_id,img_x,img_y,lat,lon180,lat_orig,lon180_orig from hollows_points_corrected",c)
def local(lat,lon,lat0,lon0): return R*np.radians(((np.asarray(lon)-lon0+180)%360)-180)*np.cos(np.radians(lat0)), R*np.radians(np.asarray(lat)-lat0)
def bbox_map(la,lo,u,v,fu,fv):
    a=1-u if fu else u; b=1-v if fv else v
    return la.max()-b*(la.max()-la.min()), lo.min()+a*(lo.max()-lo.min())
def H_fit(u,v,x,y):
    A=[]
    for uu,vv,xx,yy in zip(u,v,x,y): A+=[[uu,vv,1,0,0,0,-xx*uu,-xx*vv,-xx],[0,0,0,uu,vv,1,-yy*uu,-yy*vv,-yy]]
    return np.linalg.svd(np.array(A))[2][-1].reshape(3,3)
agree=0; n=0; errs=[]
for im,g in P.groupby("mdis_id"):
    if len(g)<6: continue
    lat0,lon0=g.lat.median(),g.lon180.median(); x,y=local(g.lat.values,g.lon180.values,lat0,lon0); H=H_fit(g.img_x.values,g.img_y.values,x,y)
    q=H@np.array([[0,1,0,1],[0,0,1,1],[1,1,1,1.]]); X,Y=q[:2]/q[2]   # corners (0,0),(1,0),(0,1),(1,1)
    la=lat0+np.degrees(Y/R); lo=lon0+np.degrees(X/(R*np.cos(np.radians(lat0))))
    B=np.c_[local(g.lat_orig.values,g.lon180_orig.values,lat0,lon0)]
    fits={(fu,fv):np.median(np.hypot(*(np.c_[local(*bbox_map(la,lo,g.img_x.values,g.img_y.values,fu,fv),lat0,lon0)]-B).T)) for fu,fv in itertools.product([0,1],[0,1])}
    best=min(fits,key=fits.get)
    pred=(int((lo[0]+lo[2])/2>(lo[1]+lo[3])/2), int((la[0]+la[1])/2<(la[2]+la[3])/2))
    n+=1; agree+=(pred==best); errs.append(fits[pred])
print(f"(1) orientation predicted by comparing corner coordinates per axis matches the best-fitting orientation in {agree} of {n} images")
print("    median error of the predicted reconstruction vs released (km), quantiles 50/90:", np.round(np.percentile(errs,[50,90]),3), "(corners from the camera model, not the label)")
print("(2) label corners interpolated within the tilted frame vs camera-model coordinates:")
for f in sorted(glob.glob(f"{HE}/data/mdis_frames/IMG/CN*_IF_5.IMG")):
    im=os.path.basename(f)[:13]; g=P[P.mdis_id==im]
    if len(g)<3: continue
    t=open(f,'rb').read(4096*3).decode('latin1')
    vec=lambda k: np.array([float(x) for x in re.findall(r'(-?\d+\.\d+)\s*<DEG>',re.search(k+r'\s*=\s*\(([^)]*)\)',t).group(1))])
    la=vec('RETICLE_POINT_LATITUDE'); lo=((vec('RETICLE_POINT_LONGITUDE')+180)%360)-180
    lat0,lon0=g.lat.median(),g.lon180.median(); C=np.c_[local(la,lo,lat0,lon0)]; u,v=g.img_x.values,g.img_y.values
    Q=((1-u)*(1-v))[:,None]*C[0]+(u*(1-v))[:,None]*C[1]+((1-u)*v)[:,None]*C[2]+(u*v)[:,None]*C[3]
    A=np.c_[local(g.lat.values,g.lon180.values,lat0,lon0)]; Bm=np.c_[local(g.lat_orig.values,g.lon180_orig.values,lat0,lon0)]
    print(f"    {im}: n={len(g)}, median distance to camera-model coordinates {np.median(np.hypot(*(Q-A).T)):.2f} km (released: {np.median(np.hypot(*(Bm-A).T)):.2f} km)")
