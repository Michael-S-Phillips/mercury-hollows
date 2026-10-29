import numpy as np, pandas as pd, json, warnings; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
m=pd.read_csv("morph_all_nested.csv"); ch=pd.read_csv("nested_children.csv"); E=json.load(open("fresh_envelope.json"))
bx,by=np.log10(E["bins_D"]),np.log10(E["bins_d90_km"])
fresh=lambda D: 10**np.interp(np.log10(D),bx,by)          # flat beyond the last bin (~2.7 km above ~250 km)
for d in (m,ch):
    d["abslat"]=d.lat.abs(); d["logD"]=np.log10(d.D); d["host"]=(d.n_all>0).astype(int); d["wallh"]=(d.n_wall>0).astype(int)
ch["tf"]=np.clip(fresh(ch.D_p)-ch.depth_p,0,None); ch["margin"]=ch.d_ex-ch.tf
onf=(ch.h_pos<400)&(ch.rho_p<0.8); rais=(ch.h_pos>=800)|(ch.rho_p>=0.85)
ch["pos"]=np.where(onf,"floor",np.where(rais,"raised","ambiguous"))
ch["group"]=np.where(onf&(ch.margin>=0),"floor_through",np.where(onf&(ch.margin<0),"floor_fill_only",np.where(rais,"raised","ambiguous")))
ch.to_csv("nested_children_classified.csv",index=False)
R={}
def summ(s): return dict(N=int(len(s)),hosts=int(s.host.sum()),rate=round(float(s.host.mean()),3) if len(s) else None,wall_hosts=int(s.wallh.sum()),D_med=round(float(s.D.median()),1) if len(s) else None)
def fit(formula,data):
    try:
        mo=smf.logit(formula,data=data).fit(disp=0,maxiter=200); return {k:[round(float(mo.params[k]),3),round(float(mo.bse[k]),3),round(float(mo.pvalues[k]),4)] for k in mo.params.index}
    except Exception as e: return str(e)
for nm,sel in [("parents_60_300",(ch.D_p>=60)&(ch.D_p<=300)),("parents_60_300_C0.2",(ch.D_p>=60)&(ch.D_p<=300)&(ch.C_mean>=0.2))]:
    c=ch[sel]; r={g:summ(c[c.group==g]) for g in ["floor_fill_only","floor_through","raised","ambiguous"]}
    tab=[]
    for a,b in zip([10,15,20,30,50],[15,20,30,50,250]):
        s=c[(c.D>=a)&(c.D<b)]; tab.append({"D":f"{a}-{b}",**{g:f"{int(s[s.group==g].host.sum())}/{int((s.group==g).sum())}" for g in ["floor_fill_only","floor_through","raised"]}})
    r["by_D"]=tab
    r["logit_groups"]=fit("host ~ C(group, Treatment('floor_fill_only')) + logD + C_mean + abslat",c[c.group!="ambiguous"])
    r["logit_margin_floor"]=fit("host ~ margin + logD + C_mean + abslat",c[c.pos=="floor"])
    R[nm]=r
# size-slope interaction: isolated vs nested on parent floor vs nested on raised terrain (all parent sizes, incl. basins)
iso=m[(m.parent<0)&~m.in_orgel&(m.D<=250)].assign(setting="isolated")
nf=ch[ch.pos=="floor"].assign(setting="nested_floor"); nr=ch[ch.pos=="raised"].assign(setting="nested_raised")
A=pd.concat([iso,nf,nr])[["host","logD","C_mean","abslat","setting","D"]]
R["slope_interaction"]=fit("host ~ logD*C(setting, Treatment('isolated')) + C_mean + abslat",A)
R["rates_by_setting_D"]=[{"D":f"{a}-{b}",**{s:f"{int(A[(A.setting==s)&(A.D>=a)&(A.D<b)].host.sum())}/{int(((A.setting==s)&(A.D>=a)&(A.D<b)).sum())}" for s in ["isolated","nested_floor","nested_raised"]}}
                         for a,b in zip([10,15,20,30,50],[15,20,30,50,250])]
json.dump(R,open("nested_test2_record.json","w"),indent=1)
for k,v in R.items(): print("=====",k); print(json.dumps(v,indent=0)[:3000])
