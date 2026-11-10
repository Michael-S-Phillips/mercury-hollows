import sqlite3, numpy as np, pandas as pd, os
HE=os.path.expanduser("~/mnt/hollows-emissivity"); R=2439.7
c=sqlite3.connect(f"{HE}/GIS/vector/bickel2024_hollows_corrected.gpkg")
P=pd.read_sql("select mdis_id,img_x,img_y,lat,lon180,lat_orig,lon180_orig from hollows_points_corrected",c)
def local(lat,lon,lat0,lon0):
    return R*np.radians(((lon-lon0+180)%360)-180)*np.cos(np.radians(lat0)), R*np.radians(lat-lat0)
def umeyama(A,B,scale=True):
    ma,mb=A.mean(0),B.mean(0); a,b=A-ma,B-mb; U,S,Vt=np.linalg.svd(b.T@a/len(A)); D=np.eye(2)
    if np.linalg.det(U@Vt)<0: D[1,1]=-1
    Rm=U@D@Vt; s=(S*np.diag(D)).sum()/a.var(0).sum() if scale else 1.0
    return s,Rm,mb-s*Rm@ma
rows=[]
for im,g in P.groupby("mdis_id"):
    if len(g)<6: continue
    lat0,lon0=g.lat.median(),g.lon180.median()
    A=np.c_[local(g.lat.values,g.lon180.values,lat0,lon0)]; B=np.c_[local(g.lat_orig.values,g.lon180_orig.values,lat0,lon0)]
    off=np.hypot(*(B-A).T); rt=np.hypot(*((B-A)-(B-A).mean(0)).T)
    s,Rm,t=umeyama(A,B,False); rr=np.hypot(*(B-(A@Rm.T*s+t)).T); ang=np.degrees(np.arctan2(Rm[1,0],Rm[0,0]))
    s2,R2,t2=umeyama(A,B,True); rs=np.hypot(*(B-(A@R2.T*s2+t2)).T)
    X=np.c_[A,np.ones(len(A))]; M,*_=np.linalg.lstsq(X,B,rcond=None); ra=np.hypot(*(B-X@M).T); det=np.linalg.det(M[:2].T)
    rows.append(dict(image=im,n=len(g),offset=np.median(off),res_transl=np.sqrt((rt**2).mean()),res_rigid=np.sqrt((rr**2).mean()),rot_deg=ang,
                     res_sim=np.sqrt((rs**2).mean()),scale=s2,res_affine=np.sqrt((ra**2).mean()),det=det))
T=pd.DataFrame(rows).set_index("image")
pd.set_option("display.width",200)
print(T.loc[["CN0223615469M","CN1002981281M","CN1017701438M"]].round(3).to_string())
print("\nall images with >=6 detections:", len(T))
for k in ["offset","res_transl","res_rigid","res_sim","res_affine"]: print(f"  median {k}: {T[k].median():.3f} km")
print("  |rotation| deg quantiles:", np.round(np.percentile(abs(T.rot_deg),[10,25,50,75,90]),1))
print("  scale quantiles:", np.round(np.percentile(T.scale,[10,50,90]),3), " det<0 fraction:", (T.det<0).mean().round(3))
T.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"released_error_per_image.csv"))
