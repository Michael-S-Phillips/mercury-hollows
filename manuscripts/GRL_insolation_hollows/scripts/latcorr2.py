import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
import numpy as np, pandas as pd
from scipy import stats
T=pd.read_csv("density_by_band.csv"); f=pd.read_csv("fields_flagged.csv")
def pois_ci(k):
    return (stats.chi2.ppf(0.025,2*k)/2 if k>0 else 0.0, stats.chi2.ppf(0.975,2*k+2)/2)
t=T[(T.hemi=="N")&(T.cls=="<=40")&(T.zone=="outNSP")].copy()
t[["lo","hi"]]=[ [x/a for x in pois_ci(k)] if a>0 else [np.nan,np.nan] for k,a in zip(t.fields,t.area_1e6km2)]
print("North, <=40 m/px, outside NSP, with 95% Poisson CI:"); print(t[["band","area_1e6km2","fields","field_density","lo","hi","hollows","hollow_density"]].round(1).to_string(index=False))
# exclude the 90-150E, 40-70N province (area and fields)
U=_REPO + ""; ppd=10
cov=np.load(f"{U}/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"].astype(float); nsp=np.load("nsp_mask_10ppd_r40.npy")
H,W=cov.shape; latc=90-(np.arange(H)+0.5)/ppd; lonc=(np.arange(W)+0.5)/ppd
cell=np.cos(np.radians(latc))[:,None]*(np.pi*2439.4/180/ppd)**2*np.ones((1,W))
prov=((latc>=40)&(latc<70))[:,None]&((lonc>=90)&(lonc<150))[None,:]
f["prov"]=(f.lat>=40)&(f.lat<70)&(np.mod(f.lon,360)>=90)&(np.mod(f.lon,360)<150)
print("\nNorth, <=40 m/px, outside NSP and outside the 90-150E/40-70N province:")
for lo in range(0,80,10):
    gm=((latc>=lo)&(latc<lo+10))[:,None]&(cov<=40)&~nsp&~prov
    A=(cell*gm).sum()/1e6; k=int(((f.lat>=lo)&(f.lat<lo+10)&(f.gpx<=40)&~f.nsp&~f.prov).sum()); ci=pois_ci(k)
    print(f"  {lo:2d}-{lo+10:<2d}N area {A:5.2f}  fields {k:4d}  density {k/A:6.1f}  [{ci[0]/A:5.1f}, {ci[1]/A:5.1f}]" if A>0 else f"  {lo}-{lo+10} no area")
print("\nSouth by resolution class (all zones):")
print(T[(T.hemi=="S")&(T.zone=="all")&(T.cls!="all<=150")].pivot(index="band",columns="cls",values=["area_1e6km2","fields","field_density"]).round(2).to_string())
