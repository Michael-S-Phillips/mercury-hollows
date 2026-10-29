import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Incidence-controlled version: searched area = cells covered by 60-75 deg incidence frames; fields counted only if
at least one member hollow was detected in a 60-75 deg frame."""
import numpy as np, pandas as pd, json
from scipy import stats
U=_REPO + "/code/rerun"; ppd=10
cov=np.load(f"{U}/latitude/incidence_null/out_inc60_75/coverage_10ppd.npz")["bestpx"].astype(float); nsp=np.load(_REPO + "/manuscripts/GRL_insolation_hollows/data/nsp_mask_10ppd_r40.npy")
det=pd.read_csv(f"{U}/out_v2/detections_qc.csv"); h=pd.read_csv(f"{U}/out_v2/hollows_canonical.csv"); f=pd.read_csv(f"{U}/out_v2/fields.csv")
ok_h=set(det.loc[(det.inc>=60)&(det.inc<=75),"hollow_id"].dropna().astype(int))
ok_f=set(h.loc[h.hollow_id.isin(ok_h),"field_id"].dropna().astype(int))
print("hollows with a 60-75 deg detection:",len(ok_h),"| fields:",len(ok_f))
f=f[f.field_id.isin(ok_f)].copy()
H,W=cov.shape; latc=90-(np.arange(H)+0.5)/ppd; lonc=(np.arange(W)+0.5)/ppd
cell=np.cos(np.radians(latc))[:,None]*(np.pi*2439.4/180/ppd)**2*np.ones((1,W))
r=((90-f.lat)*ppd).astype(int).clip(0,H-1); c=(np.mod(f.lon,360)*ppd).astype(int).clip(0,W-1)
f["gpx"]=cov[r,c]; f["nsp"]=nsp[r,c]; f["prov"]=(f.lat>=40)&(f.lat<70)&(np.mod(f.lon,360)>=90)&(np.mod(f.lon,360)<150)
prov=((latc>=40)&(latc<70))[:,None]&((lonc>=90)&(lonc<150))[None,:]
def pci(k): return (stats.chi2.ppf(0.025,2*k)/2 if k>0 else 0.0, stats.chi2.ppf(0.975,2*k+2)/2)
R={"n_hollows_inc60_75":len(ok_h),"n_fields_inc60_75":len(ok_f)}
for px in [40,80]:
  for zone in ["out","in"]:
    for ep in [False,True]:
      key=f"le{px}_{zone}{'_noprov' if ep else ''}"; R[key]={}
      for lo,hi in [(0,30),(30,40),(40,60),(60,80)]:
        gm=((latc>=lo)&(latc<hi))[:,None]&(cov<=px)&(nsp if zone=="in" else ~nsp)
        fm=(f.lat>=lo)&(f.lat<hi)&(f.gpx<=px)&(f.nsp if zone=="in" else ~f.nsp)
        if ep: gm&=~prov; fm&=~f.prov
        A=(cell*gm).sum()/1e6; k=int(fm.sum()); ci=pci(k); ok=A>0.02
        R[key][f"{lo}-{hi}"]=dict(A=round(A,3),k=k,d=round(k/A,1) if ok else None,lo=round(ci[0]/A,1) if ok else None,hi=round(ci[1]/A,1) if ok else None)
      print(key, {b:(v['k'],v['A'],v['d'],v['lo'],v['hi']) for b,v in R[key].items()})
for key in ["le40_out","le40_out_noprov","le80_out","le80_out_noprov"]:
    a=R[key]["60-80"]; b=R[key]["40-60"]; n=a["k"]+b["k"]; p0=a["A"]/(a["A"]+b["A"]); bt=stats.binomtest(a["k"],n,p0); ci=bt.proportion_ci()
    R[key+"_rr"]=dict(rr=round((a["k"]/a["A"])/(b["k"]/b["A"]),2),ci=[round(ci.low/(1-ci.low)*b["A"]/a["A"],2),round(ci.high/(1-ci.high)*b["A"]/a["A"],2)],p=round(float(bt.pvalue),3))
    print(key,"60-80 vs 40-60:",R[key+"_rr"])
json.dump(R,open(_REPO + "/manuscripts/GRL_insolation_hollows/records/latcorr_inc60_75_record.json","w"),indent=1)
