import sqlite3, numpy as np, pandas as pd, os
HE=os.path.expanduser("~/mnt/hollows-emissivity"); R=2439.7
c=sqlite3.connect(f"{HE}/GIS/vector/bickel2024_hollows_corrected.gpkg")
P=pd.read_sql("select mdis_id,img_x,img_y,lat,lon180,lat_orig,lon180_orig,emission,px_scale_m from hollows_points_corrected",c)
T=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"released_error_per_image.csv")).set_index("image")
def local(lat,lon,lat0,lon0): return R*np.radians(((np.asarray(lon)-lon0+180)%360)-180)*np.cos(np.radians(lat0)), R*np.radians(np.asarray(lat)-lat0)
def H_fit(u,v,x,y):
    A=[]
    for uu,vv,xx,yy in zip(u,v,x,y): A+=[[uu,vv,1,0,0,0,-xx*uu,-xx*vv,-xx],[0,0,0,uu,vv,1,-yy*uu,-yy*vv,-yy]]
    return np.linalg.svd(np.array(A))[2][-1].reshape(3,3)
out=[]
for im,g in P.groupby("mdis_id"):
    if im not in T.index: continue
    lat0,lon0=g.lat.median(),g.lon180.median(); x,y=local(g.lat.values,g.lon180.values,lat0,lon0)
    H=H_fit(g.img_x.values,g.img_y.values,x,y); Hi=np.linalg.inv(H)
    def toimg(X,Y): q=Hi@np.array([X,Y,1.0]); return q[:2]/q[2]
    q0=H@np.array([.5,.5,1]); cx,cy=q0[:2]/q0[2]
    a=toimg(cx,cy); b=toimg(cx,cy+0.5)  # 0.5 km north of center
    d=(b-a)*1024  # pixels: d[0] along samples (+right), d[1] along lines (+down)
    na=np.degrees(np.arctan2(d[1],d[0]))%360   # clockwise from +sample (screen coords, y down)
    out.append(dict(image=im,north_az=na,rot=T.loc[im,"rot_deg"],scale=T.loc[im,"scale"],n=len(g),emi=g.emission.median(),lat=lat0,px=g.px_scale_m.median(),res_sim=T.loc[im,"res_sim"]))
O=pd.DataFrame(out).set_index("image")
print(O.loc[["CN0223615469M","CN1002981281M","CN1017701438M"],["north_az","rot","scale"]].round(2))
pred=((90-O.north_az+180)%360)-180
err=((O.rot-pred+180)%360)-180
print("images:",len(O)," |rot - (90 - north_az)| quantiles (deg):", np.round(np.percentile(abs(err),[25,50,75,90]),2))
print("corr(scale, 1/cos(lat)):", np.corrcoef(O.scale,1/np.cos(np.radians(O.lat)))[0,1].round(3), " corr(scale,emi):",np.corrcoef(O.scale,O.emi)[0,1].round(3), " corr(scale, |rot|):",np.corrcoef(O.scale,abs(O.rot))[0,1].round(3))
# try scale vs 1/cos(rot) and 1/cos(lat) etc
for name,f in [("1/cos(lat)",1/np.cos(np.radians(O.lat))),("1+|sin(rot)|",1+abs(np.sin(np.radians(O.rot)))),("1/cos(emi)",1/np.cos(np.radians(O.emi)))]:
    print(f"  median scale/{name}: {np.median(O.scale/f):.3f}  IQR {np.percentile(O.scale/f,25):.3f}-{np.percentile(O.scale/f,75):.3f}")
O.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"rot_vs_north.csv"))
