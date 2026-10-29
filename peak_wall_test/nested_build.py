import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Nested craters: parent assignment, parent fill estimate from depth deficit, child position on the parent floor or raised terrain,
child excavation depth. Writes nested_children.csv and fresh_envelope.json."""
import os, json, numpy as np, pandas as pd, rasterio
from rasterio.windows import from_bounds
HE=os.path.expanduser(_REPO + ""); DEM=f"{HE}/data/DEMs/Mercury_Messenger_USGS_DEM_Global_665m_v2_cog.tif"
RM=2439.4; KD=np.pi*RM/180
m=pd.read_csv("morph_all.csv"); m["lon360"]=m.lon%360
org=pd.read_csv("orgel_basins.csv"); org=org[org.Confidence.isin(["Certain","Probable"])]
def gcd(la1,lo1,la2,lo2):
    p1,p2=np.radians(la1),np.radians(la2); return RM*np.arccos(np.clip(np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(np.radians(lo2-lo1)),-1,1))
big=m[m.D>=60].reset_index(drop=True)
par=[];prho=[]
for r in m.itertuples():
    c=big[(big.D>=2.5*r.D)&(np.abs(big.lat-r.lat)<big.D/2/KD+0.1)]
    if len(c):
        d=gcd(r.lat,r.lon360,c.lat.values,c.lon360.values); ok=d<0.95*c.D.values/2
        if ok.any():
            j=np.argmin(np.where(ok,c.D.values,np.inf)); par.append(int(c.uid.values[j])); prho.append(d[j]/(c.D.values[j]/2)); continue
    par.append(-1); prho.append(np.nan)
m["parent"]=par; m["rho_p"]=prho
m["in_orgel"]=[bool((gcd(r.lat,r.lon360,org.y_coord.values.astype(float),org.x_coord.values.astype(float)%360)<org.Diam_km.values.astype(float)/2).any()) for r in m.itertuples()]
# fresh depth envelope from isolated craters: 90th percentile in log-D bins, power-law fit
iso=m[(m.parent<0)&~m.in_orgel&(m.depth>0)]
edges=np.logspace(np.log10(15),np.log10(400),12); xb=[];yb=[];nb=[]
for a,b in zip(edges[:-1],edges[1:]):
    s=iso[(iso.D>=a)&(iso.D<b)].depth
    if len(s)>=8: xb.append(np.sqrt(a*b)); yb.append(np.percentile(s,90)/1000); nb.append(len(s))
B,A=np.polyfit(np.log10(xb),np.log10(yb),1); fresh=lambda D: 10**A*D**B
json.dump(dict(coef=10**A,exponent=B,bins_D=xb,bins_d90_km=yb,n=nb,note="d_fresh(km)=coef*D^exponent, 90th percentile of DEM depth of isolated craters"),open("fresh_envelope.json","w"),indent=1)
print("fresh envelope d = %.3f D^%.3f  (D=50: %.2f km, D=100: %.2f, D=300: %.2f, D=1500: %.2f)"%(10**A,B,fresh(50),fresh(100),fresh(300),fresh(1500)))
P=m.set_index("uid")
ch=m[m.parent>=0].copy()
ch["D_p"]=P.loc[ch.parent,"D"].values; ch["depth_p"]=P.loc[ch.parent,"depth"].values/1000; ch["floor_p"]=P.loc[ch.parent,"floor"].values
ch["name_p"]=P.loc[ch.parent,"name"].values; ch["filled_p"]=P.loc[ch.parent,"filled"].values
ch["t_fill"]=np.clip(fresh(ch.D_p)-ch.depth_p,0,None)
ch["d_ex"]=0.1*12**0.15*ch.D**0.85          # excavation depth, 0.1 x transient diameter (Croft 1985 scaling, Dt = 12 km)
ch["su"]=0.086*ch.D**1.03                    # stratigraphic uplift of the central peak (terrestrial scaling)
# surrounding pre-impact surface of the child relative to the parent floor
xs=lambda lo: np.radians(lo-180)*2439400.0; ys=lambda la: np.radians(la)*2439400.0
hs=[]
with rasterio.open(DEM) as ds:
    for r in ch.itertuples():
        R=r.D/2; cl=np.cos(np.radians(r.lat)); dla=1.9*R/KD; dlo=1.9*R/(KD*max(cl,0.05))
        if r.lon360-dlo<0 or r.lon360+dlo>360: hs.append(np.nan); continue
        win=from_bounds(xs(r.lon360-dlo),ys(r.lat-dla),xs(r.lon360+dlo),ys(r.lat+dla),ds.transform)
        z=ds.read(1,window=win,boundless=True,fill_value=-32768).astype(float); z[z==-32768]=np.nan; z*=0.5
        T=ds.window_transform(win); nr,nc=z.shape
        LON=180+np.degrees((T.c+T.a*(np.arange(nc)+.5))/2439400.0)[None,:]; LAT=np.degrees((T.f+T.e*(np.arange(nr)+.5))/2439400.0)[:,None]
        rho=np.hypot((LON-r.lon360)*KD*cl,(LAT-r.lat)*KD)/R; s=z[(rho>1.3)&(rho<1.8)]; s=s[np.isfinite(s)]
        hs.append(np.median(s)-r.floor_p if len(s) else np.nan)
ch["h_pos"]=hs
ch.to_csv("nested_children.csv",index=False); m.to_csv("morph_all_nested.csv",index=False)
print(len(ch),"nested children;",(ch.n_all>0).sum(),"hosts")
print(ch[["D","D_p","t_fill","d_ex","h_pos","rho_p"]].describe().round(2).to_string())
