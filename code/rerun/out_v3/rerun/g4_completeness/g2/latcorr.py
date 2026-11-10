# g4 copy of code/rerun/out_v3/rerun/g2_insolation/latcorr.py (6 Nov 2026). Changed: G2_ROOT depth 5 -> 6 (see make_copies_g4.py).
# g2 copy of manuscripts/GRL_insolation_hollows/scripts/latcorr.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Hollow-field density and per-crater host fraction by latitude, corrected for NAC image coverage and the northern smooth plains (NSP)."""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 6))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_GRLDATA = f"{G2_ROOT}/manuscripts/GRL_insolation_hollows/data"
# --- end g2 block ---
import numpy as np, pandas as pd, json
from scipy import stats
import statsmodels.formula.api as smf
U=G2_ROOT
ppd=10; H=180*ppd; W=360*ppd
cov=np.load(f"{U}/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"].astype(float)     # NaN = not searched
nsp=np.load(f"{G2_GRLDATA}/nsp_mask_10ppd_r40.npy")
latc=90-(np.arange(H)+0.5)/ppd
cell_km2=np.cos(np.radians(latc))[:,None]*(np.pi*2439.4/180/ppd)**2*np.ones((1,W))
def rc(lat,lon):
    return ((90-np.asarray(lat))*ppd).astype(int).clip(0,H-1), ((np.mod(lon,360))*ppd).astype(int).clip(0,W-1)
f=pd.read_csv(f"{IN}/fields.csv"); h=pd.read_csv(f"{IN}/hollows_canonical.csv")
for d in (f,h):
    r,c=rc(d.lat,d.lon); d["gpx"]=cov[r,c]; d["nsp"]=nsp[r,c]
print("fields",len(f),"in NSP",int(f.nsp.sum()),"| hollows",len(h),"in NSP",int(h.nsp.sum()), "| fields at unsearched cells", int(f.gpx.isna().sum()))
classes={"all<=150":(0,150.1),"<=40":(0,40),"40-80":(40,80),"80-150":(80,150.1)}
def pois_ci(k):
    lo=stats.chi2.ppf(0.025,2*k)/2 if k>0 else 0.0; hi=stats.chi2.ppf(0.975,2*k+2)/2; return lo,hi
rows=[]
bands=[(a,a+10) for a in range(0,90,10)]
for hemi,sgn in [("N",1),("S",-1)]:
    for lo,hi in bands:
        bl=(latc*sgn>=lo)&(latc*sgn<hi)
        for cname,(p0,p1) in classes.items():
            gm=bl[:,None]&(cov>p0)&(cov<=p1)
            for zone,zm in [("all",np.ones_like(nsp)),("outNSP",~nsp),("inNSP",nsp)]:
                A=(cell_km2*(gm&zm)).sum()
                sel=lambda d:(d.lat*sgn>=lo)&(d.lat*sgn<hi)&(d.gpx>p0)&(d.gpx<=p1)&({"all":True,"outNSP":~d.nsp,"inNSP":d.nsp}[zone])
                kf=int(sel(f).sum()); kh=int(sel(h).sum())
                dens=kf/A*1e6 if A>0 else np.nan; ci=[x/A*1e6 for x in pois_ci(kf)] if A>0 else [np.nan,np.nan]
                rows.append(dict(hemi=hemi,band=f"{lo}-{hi}",cls=cname,zone=zone,area_1e6km2=A/1e6,fields=kf,hollows=kh,
                                 field_density=dens,ci_lo=ci[0],ci_hi=ci[1],hollow_density=kh/A*1e6 if A>0 else np.nan))
T=pd.DataFrame(rows); T.to_csv(f"{G2_OUT}/density_by_band.csv",index=False)
pd.set_option("display.width",250)
for cname in ["all<=150","<=40","40-80","80-150"]:
    t=T[(T.hemi=="N")&(T.cls==cname)].pivot(index="band",columns="zone",values=["area_1e6km2","fields","field_density"]).round(2)
    print(f"\n=== North, searched class {cname} m/px: area (1e6 km2), fields, fields per 1e6 km2 ==="); print(t.to_string())
t=T[(T.hemi=="S")&(T.cls=="all<=150")&(T.zone=="all")][["band","area_1e6km2","fields","field_density","ci_lo","ci_hi"]].round(2)
print("\n=== South, all searched ==="); print(t.to_string(index=False))
f.to_csv(f"{G2_OUT}/fields_flagged.csv",index=False)
