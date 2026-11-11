# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_v3/dT_plains.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Plains effect on crater hosting at 40-80N, craters >= 10 km whose center cell was imaged at <= 80 m/pixel
(reproduces percrater_record.json nsp_model: n = 1,179, 107 hosts), with dT in place of longitude distance."""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_INSOLDATA = f"{G2_ROOT}/manuscripts/insolation_thermal/data"
# --- end g2 block ---
import os
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
G=os.path.expanduser("~/mnt/hollows-emissivity")+"/manuscripts/insolation_thermal"; sys.path.insert(0,f"{G}/scripts"); import thermal
s = pd.read_csv(f"{IN}/test2_craters_h18.csv")
z=np.load(f"{G}/data/thermal_grid.npz"); la,lo,T=z["la"],z["lo"],z["G"]
s["dT"]=thermal.lookup(la,lo,T,s.latitude.values,s.longitude.values)-thermal.lookup(la,lo,T,s.latitude.values,np.full(len(s),90.0))
lo180 = s.longitude.values % 180.0; s["hot"] = np.abs(np.where(lo180 > 90, 180 - lo180, lo180) - 90)
nsp=np.load(f"{G2_ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_v3/nsp_mask_10ppd_r40.npy"); ppd=10
cov=np.load(os.path.expanduser("~/mnt/hollows-emissivity")+"/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"].astype(float)
r=((90-s.latitude.values)*ppd).astype(int).clip(0,1799); c=((np.mod(s.lon360.values,360))*ppd).astype(int).clip(0,3599)
s["nsp"]=nsp[r,c].astype(int); s["lat"]=s.latitude; s["gpx"]=cov[r,c]
d=s[(s.diameter>=10)&(s.latitude>=40)&(s.latitude<80)&(s.gpx<=80)].copy()
out={}
for name,f in [("hot","hosts ~ logD + nsp + hot + lat"),("dT","hosts ~ logD + nsp + dT + lat")]:
    m=smf.logit(f,data=d).fit(disp=0)
    out[name]=dict(n=int(len(d)),hosts=int(d.hosts.sum()),n_nsp=int(d.nsp.sum()),hosts_nsp=int(d[d.nsp==1].hosts.sum()),OR=float(np.exp(m.params.nsp)),ci=[float(np.exp(m.params.nsp-1.96*m.bse.nsp)),float(np.exp(m.params.nsp+1.96*m.bse.nsp))],p=float(m.pvalues.nsp),aic=float(m.aic))
print(json.dumps(out,indent=1)); json.dump(out,open(f"{G2_OUT}/dT_plains_record.json","w"),indent=1)
