# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_v3/transition_layer_hilat.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Layer-thickness (L368) vs dT as predictors of hollow hosting, including latitudes above 60 deg.
Usage: python3 transition_layer_hilat.py <dir with thermal.py and thermal_grid.npz> <test2_craters_h18.csv> <transition_layer_grid.npz>
Writes transition_layer_hilat_record.json in the current directory."""
# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---
import os as _g2_os
G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))
G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"
_g2_os.makedirs(G2_OUT, exist_ok=True)
G2_GRLDATA = f"{G2_ROOT}/manuscripts/GRL_insolation_hollows/data"
# --- end g2 block ---
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
TH, CR, GR = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, TH); import thermal as th
z=np.load(f"{G2_GRLDATA}/thermal_grid.npz"); la,lo,T=z["la"],z["lo"],z["G"]
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
json.dump(out,open(f"{G2_OUT}/transition_layer_hilat_record.json","w"),indent=1)
bands={}
for hemi,cond in [("N",s.latitude>0),("S",s.latitude<0)]:
  for a0,a1 in [(60,75),(75,90)]:
    d=s[cond&(s.abslat>a0)&(s.abslat<=a1)]
    for lab,c in [("lt3cm",d.L368<3),("3to6cm",(d.L368>=3)&(d.L368<6)),("ge6cm",d.L368>=6)]:
      x=d[c]; bands[f"{hemi}{a0}-{a1}_{lab}"]=dict(n=int(len(x)),hosts=int(x.hosts.sum()))
out["highlat_bands"]=bands
out["note"]="dT model vs L368 model: logit hosts ~ X*logD + abslat; craters >=5 km from test2_craters_h18.csv; L368 from transition_layer_grid.npz"
json.dump(out,open(f"{G2_OUT}/transition_layer_hilat_record.json","w"),indent=1)
