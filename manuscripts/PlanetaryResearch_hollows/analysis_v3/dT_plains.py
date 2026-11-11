import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Plains effect on crater hosting at 40-80N, craters >= 10 km whose center cell was imaged at <= 80 m/pixel
(reproduces percrater_record.json nsp_model: n = 1,179, 107 hosts), with dT in place of longitude distance."""
import os
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
G=os.path.expanduser(_REPO + "")+"/manuscripts/insolation_thermal"; sys.path.insert(0,f"{G}/scripts"); import thermal
s = pd.read_csv(f"{G}/data/test2_craters_h18.csv")
z=np.load(f"{G}/data/thermal_grid.npz"); la,lo,T=z["la"],z["lo"],z["G"]
s["dT"]=thermal.lookup(la,lo,T,s.latitude.values,s.longitude.values)-thermal.lookup(la,lo,T,s.latitude.values,np.full(len(s),90.0))
lo180 = s.longitude.values % 180.0; s["hot"] = np.abs(np.where(lo180 > 90, 180 - lo180, lo180) - 90)
nsp=np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)),"nsp_mask_10ppd_r40.npy")); ppd=10
cov=np.load(os.path.expanduser(_REPO + "")+"/code/rerun/out_v2/coverage_10ppd.npz")["bestpx"].astype(float)
r=((90-s.latitude.values)*ppd).astype(int).clip(0,1799); c=((np.mod(s.lon360.values,360))*ppd).astype(int).clip(0,3599)
s["nsp"]=nsp[r,c].astype(int); s["lat"]=s.latitude; s["gpx"]=cov[r,c]
d=s[(s.diameter>=10)&(s.latitude>=40)&(s.latitude<80)&(s.gpx<=80)].copy()
out={}
for name,f in [("hot","hosts ~ logD + nsp + hot + lat"),("dT","hosts ~ logD + nsp + dT + lat")]:
    m=smf.logit(f,data=d).fit(disp=0)
    out[name]=dict(n=int(len(d)),hosts=int(d.hosts.sum()),n_nsp=int(d.nsp.sum()),hosts_nsp=int(d[d.nsp==1].hosts.sum()),OR=float(np.exp(m.params.nsp)),ci=[float(np.exp(m.params.nsp-1.96*m.bse.nsp)),float(np.exp(m.params.nsp+1.96*m.bse.nsp))],p=float(m.pvalues.nsp),aic=float(m.aic))
print(json.dumps(out,indent=1)); json.dump(out,open("dT_plains_record.json","w"),indent=1)
