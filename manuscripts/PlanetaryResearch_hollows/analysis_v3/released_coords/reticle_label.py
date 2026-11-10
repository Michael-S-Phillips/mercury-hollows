import re, glob, os, itertools, sqlite3, numpy as np, pandas as pd
HE=os.path.expanduser("~/mnt/hollows-emissivity"); R=2439.7
c=sqlite3.connect(f"{HE}/GIS/vector/bickel2024_hollows_corrected.gpkg")
P=pd.read_sql("select mdis_id,img_x,img_y,lat,lon180,lat_orig,lon180_orig from hollows_points_corrected",c)
def lab(f):
    t=open(f,'rb').read(4096*3).decode('latin1')
    def vec(k):
        m=re.search(k+r'\s*=\s*\(([^)]*)\)',t); return [float(x) for x in re.findall(r'(-?\d+\.\d+)\s*<DEG>',m.group(1))]
    return np.array(vec('RETICLE_POINT_LATITUDE')), np.array(vec('RETICLE_POINT_LONGITUDE')), float(re.search(r'TWIST_ANGLE\s*=\s*([-\d.]+)',t).group(1)), float(re.search(r'NORTH_AZIMUTH\s*=\s*([-\d.]+)',t).group(1))
def local(lat,lon,lat0,lon0): return R*np.radians(((np.asarray(lon)-lon0+180)%360)-180)*np.cos(np.radians(lat0)), R*np.radians(np.asarray(lat)-lat0)
def bilin(C,u,v): return ((1-u)*(1-v))[:,None]*C[0]+(u*(1-v))[:,None]*C[1]+((1-u)*v)[:,None]*C[2]+(u*v)[:,None]*C[3]
for f in sorted(glob.glob(f"{HE}/data/mdis_frames/IMG/CN*_IF_5.IMG")):
    im=os.path.basename(f)[:13]; g=P[P.mdis_id==im]
    if len(g)<3: continue
    la,lo,tw,na=lab(f); lo180=((lo+180)%360)-180
    lat0,lon0=g.lat.median(),g.lon180.median()
    C=np.c_[local(la,lo180,lat0,lon0)]; B=np.c_[local(g.lat_orig.values,g.lon180_orig.values,lat0,lon0)]; A=np.c_[local(g.lat.values,g.lon180.values,lat0,lon0)]
    best=min(((np.median(np.hypot(*(bilin(C[list(p)],g.img_x.values,g.img_y.values)-B).T)),p) for p in itertools.permutations(range(4))))
    # also: bilinear with the best order vs recomputed (camera model)
    pr=bilin(C[list(best[1])],g.img_x.values,g.img_y.values)
    # camera-model corners: fit homography on recomputed points, evaluate at the frame corners
    import numpy.linalg as L
    u,v=g.img_x.values,g.img_y.values; Xs=[];
    A2=[]
    for (uu,vv),(x,y) in zip(zip(u,v),A): A2+= [[uu,vv,1,0,0,0,-x*uu,-x*vv,-x],[0,0,0,uu,vv,1,-y*uu,-y*vv,-y]]
    h=L.svd(np.array(A2))[2][-1].reshape(3,3)
    cc=np.array([[0,0],[1,0],[0,1],[1,1]]); q=(h@np.c_[cc,np.ones(4)].T); cam=(q[:2]/q[2]).T
    Cb=C[list(best[1])]
    print(f"{im} n={len(g)} reticle-bilinear vs released: median {best[0]*1000:.1f} m, corner order {best[1]};  twist {tw:.1f}, north az {na:.1f}")
    print("    reticle corners vs camera-model corners (km):", np.round(np.hypot(*(Cb-cam).T),2), "  frame width (camera) km:", round(np.hypot(*(cam[1]-cam[0])),1), " reticle width km:", round(np.hypot(*(Cb[1]-Cb[0])),1))
