# g2 copy of manuscripts/insolation_thermal/scripts/latcorr3.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_INSOLDATA = f"{G2_ROOT}/manuscripts/insolation_thermal/data"
# --- end g2 block ---
import numpy as np, pandas as pd, json
from scipy import stats
import statsmodels.formula.api as smf
U=G2_ROOT; ppd=10
T=pd.read_csv(f"{G2_OUT}/density_by_band.csv"); f=pd.read_csv(f"{G2_OUT}/fields_flagged.csv")
cov=np.load(f"{U}/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"].astype(float); nsp=np.load(f"{G2_INSOLDATA}/nsp_mask_10ppd_r40.npy")
H,W=cov.shape; latc=90-(np.arange(H)+0.5)/ppd; lonc=(np.arange(W)+0.5)/ppd
cell=np.cos(np.radians(latc))[:,None]*(np.pi*2439.4/180/ppd)**2*np.ones((1,W))
def pci(k): return (stats.chi2.ppf(0.025,2*k)/2 if k>0 else 0.0, stats.chi2.ppf(0.975,2*k+2)/2)
prov=((latc>=40)&(latc<70))[:,None]&((lonc>=90)&(lonc<150))[None,:]
f["prov"]=(f.lat>=40)&(f.lat<70)&(np.mod(f.lon,360)>=90)&(np.mod(f.lon,360)<150)
R={}
def dens(sgn,lo,hi,pxmax,zone,excl_prov=False):
    gm=((latc*sgn>=lo)&(latc*sgn<hi))[:,None]&(cov<=pxmax)
    fm=(f.lat*sgn>=lo)&(f.lat*sgn<hi)&(f.gpx<=pxmax)
    if zone=="out": gm&=~nsp; fm&=~f.nsp
    if zone=="in": gm&=nsp; fm&=f.nsp
    if excl_prov: gm&=~prov; fm&=~f.prov
    A=(cell*gm).sum()/1e6; k=int(fm.sum())
    if A<0.02: return dict(A=A,k=k,d=None,lo=None,hi=None)
    c=pci(k); return dict(A=round(A,3),k=k,d=round(k/A,1),lo=round(c[0]/A,1),hi=round(c[1]/A,1))
for name,args in [("N_le40_out",(1,40,"out",False)),("N_le40_in",(1,40,"in",False)),("N_le40_out_noprov",(1,40,"out",True)),
                  ("N_le150_out",(1,150,"out",False)),("S_le40",(-1,40,"all",False)),("S_le80",(-1,80,"all",False)),("S_le150",(-1,150,"all",False))]:
    sgn,px,zone,ep=args; R[name]={f"{lo}-{lo+10}":dens(sgn,lo,lo+10,px,zone,ep) for lo in range(0,80,10)}
    print(name); print("  "+" | ".join(f"{b}: {v['k']}/{v['A']}={v['d']} [{v['lo']},{v['hi']}]" for b,v in R[name].items()))
# coarser bands for the high north
for name,(lo,hi) in {"40-60":(40,60),"60-80":(60,80),"0-30":(0,30)}.items():
    for z in ["out","in"]:
        R[f"N_le40_{z}_{name}"]=dens(1,lo,hi,40,z); R[f"N_le40_{z}_{name}_noprov"]=dens(1,lo,hi,40,z,True)
        print(name,z,R[f"N_le40_{z}_{name}"],"noprov",R[f"N_le40_{z}_{name}_noprov"])
# rate ratio 60-80 vs 40-60 outside NSP (<=40 m), exact conditional binomial test
for tag in ["","_noprov"]:
    a=R[f"N_le40_out_60-80{tag}"]; b=R[f"N_le40_out_40-60{tag}"]
    n=a["k"]+b["k"]; p0=a["A"]/(a["A"]+b["A"]); bt=stats.binomtest(a["k"],n,p0)
    rr=(a["k"]/a["A"])/(b["k"]/b["A"]); ci=bt.proportion_ci()
    rlo=(ci.low/(1-ci.low))*(b["A"]/a["A"]); rhi=(ci.high/(1-ci.high))*(b["A"]/a["A"])
    R[f"rr_60-80_vs_40-60_out{tag}"]=dict(rr=round(rr,2),ci=[round(rlo,2),round(rhi,2)],p=float(bt.pvalue))
    print("rate ratio 60-80 / 40-60 outside NSP",tag, R[f"rr_60-80_vs_40-60_out{tag}"])
json.dump(R,open(f"{G2_OUT}/latcorr_density_record.json","w"),indent=1)
