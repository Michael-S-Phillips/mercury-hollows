"""Layer-thickness (L368) vs dT as predictors of hollow hosting, including latitudes above 60 deg.
Usage: python3 transition_layer_hilat.py <GRL scripts dir with thermal.py; the grid is read from ../data/thermal_grid.npz> <test2_craters_h18.csv> <transition_layer_grid.npz>
Writes transition_layer_hilat_record.json in the current directory."""
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
TH, CR, GR = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, TH); import thermal as th
z=np.load(f"{TH}/../data/thermal_grid.npz"); la,lo,T=z["la"],z["lo"],z["G"]  # archive: thermal.py in scripts/, grid in data/
g=np.load(GR)
s=pd.read_csv(CR); s=s[s.diameter>=5].copy()
s["dT"]=th.lookup(la,lo,T,s.latitude.values,s.longitude.values)-th.lookup(la,lo,T,s.latitude.values,np.full(len(s),90.0))
s["L368"]=100*th.lookup(g["la"],g["lo"],g["L368"],s.latitude.values,s.longitude.values)
out={}
for amax in (60,75,90):
    d=s[s.abslat<=amax]
    r={}
    for name,f in [("dT","hosts ~ dT * logD + abslat"),("L368","hosts ~ L368 * logD + abslat"),("both","hosts ~ dT * logD + L368 * logD + abslat")]:
        m=smf.logit(f,data=d).fit(disp=0); r[name]=dict(aic=round(m.aic,1),p={k:float(f"{v:.2g}") for k,v in m.pvalues.items()})
    out[f"abslat_le{amax}"]=dict(n=int(len(d)),hosts=int(d.hosts.sum()),corr=round(float(np.corrcoef(d.L368,d.dT)[0,1]),3),fits=r)
    # hosts and craters where L368 < 3 cm
    t=d[d.L368<3]; out[f"abslat_le{amax}"]["thin_layer_lt3cm"]=dict(n=int(len(t)),hosts=int(t.hosts.sum()))
print(json.dumps(out,indent=1))
json.dump(out,open("transition_layer_hilat_record.json","w"),indent=1)
bands={}
for hemi,cond in [("N",s.latitude>0),("S",s.latitude<0)]:
  for a0,a1 in [(60,75),(75,90)]:
    d=s[cond&(s.abslat>a0)&(s.abslat<=a1)]
    for lab,c in [("lt3cm",d.L368<3),("3to6cm",(d.L368>=3)&(d.L368<6)),("ge6cm",d.L368>=6)]:
      x=d[c]; bands[f"{hemi}{a0}-{a1}_{lab}"]=dict(n=int(len(x)),hosts=int(x.hosts.sum()))
out["highlat_bands"]=bands
out["note"]="dT model vs L368 model: logit hosts ~ X*logD + abslat; craters >=5 km from test2_craters_h18.csv; L368 from transition_layer_grid.npz"
json.dump(out,open("transition_layer_hilat_record.json","w"),indent=1)
