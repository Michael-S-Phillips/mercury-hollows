import numpy as np, pandas as pd, json, warnings; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
m=pd.read_csv("morph_all_nested.csv"); ch=pd.read_csv("nested_children.csv")
m["abslat"]=m.lat.abs(); ch["abslat"]=ch.lat.abs(); ch["logD"]=np.log10(ch.D); m["logD"]=np.log10(m.D)
for d in (m,ch):
    d["host"]=(d.n_all>0).astype(int); d["wallh"]=(d.n_wall>0).astype(int); d["peakh"]=(d.n_peak>0).astype(int)
def classify(ch,hf=400,hr=800,cap=None,scale=1.0):
    c=ch.copy(); tf=c.t_fill*scale
    if cap: 
        fresh=lambda D: 0.245*np.minimum(D,cap)**0.481
        tf=np.clip(fresh(c.D_p)-c.depth_p,0,None)*scale
    c["tf"]=tf; c["margin"]=c.d_ex-tf
    onf=(c.h_pos<hf)&(c.rho_p<0.8); rais=(c.h_pos>=hr)|(c.rho_p>=0.85)
    c["group"]=np.where(onf&(c.margin>=0),"floor_through",np.where(onf&(c.margin<0),"floor_fill_only",np.where(rais,"raised","ambiguous")))
    return c
R={}
iso=m[(m.parent<0)&~m.in_orgel&(m.D<=250)]
def rates(s): return dict(N=len(s),host=round(s.host.mean(),3),wall=round(s.wallh.mean(),3),peak=round(s.peakh.mean(),3),n_hosts=int(s.host.sum()),
                         wall_dens=round(1000*s.n_wall.sum()/max(s.A_wall.sum(),1),2),D_med=round(s.D.median(),1),C_med=round(s.C_mean.median(),2))
for nm,kw,sub in [("main_Dp<=400",{},lambda c:c.D_p<=400),("hf300_hr600",dict(hf=300,hr=600),lambda c:c.D_p<=400),
                  ("fill_x0.5",dict(scale=0.5),lambda c:c.D_p<=400),("fill_x1.5",dict(scale=1.5),lambda c:c.D_p<=400),
                  ("basins_Dp>400_capped400",dict(cap=400),lambda c:c.D_p>400),("all_capped400",dict(cap=400),lambda c:c.D_p>0)]:
    c=classify(ch,**kw); c=c[sub(c)&(c.C_mean>=0.2)]
    r={"groups":{g:rates(c[c.group==g]) for g in ["floor_fill_only","floor_through","raised","ambiguous"]}}
    # D-matched comparison: within D bins, compare host rates
    bins=[15,20,30,50,250]; tab=[]
    for a,b in zip(bins[:-1],bins[1:]):
        s=c[(c.D>=a)&(c.D<b)]
        tab.append({"D":f"{a}-{b}",**{g:f"{int(s[s.group==g].host.sum())}/{(s.group==g).sum()}" for g in ["floor_fill_only","floor_through","raised"]}})
    r["by_D"]=tab
    f=c[c.group!="ambiguous"].copy()
    try:
        mod=smf.logit("host ~ C(group, Treatment('floor_fill_only')) + logD + C_mean + abslat",data=f).fit(disp=0)
        r["logit_group"]={k:[round(float(mod.params[k]),3),round(float(mod.bse[k]),3),float(mod.pvalues[k])] for k in mod.params.index}
    except Exception as e: r["logit_group"]=str(e)
    fl=c[c.group.str.startswith("floor")].copy()
    try:
        mod2=smf.logit("host ~ margin + logD + C_mean + abslat",data=fl).fit(disp=0)
        r["logit_margin_floor"]={k:[round(float(mod2.params[k]),3),round(float(mod2.bse[k]),3),float(mod2.pvalues[k])] for k in mod2.params.index}
    except Exception as e: r["logit_margin_floor"]=str(e)
    R[nm]=r
isoC=iso[iso.C_mean>=0.2]; R["isolated_baseline"]=rates(isoC)
bins=[15,20,30,50,250]; R["isolated_by_D"]=[{"D":f"{a}-{b}","rate":f"{int(isoC[(isoC.D>=a)&(isoC.D<b)].host.sum())}/{((isoC.D>=a)&(isoC.D<b)).sum()}"} for a,b in zip(bins[:-1],bins[1:])]
json.dump(R,open("nested_test_record.json","w"),indent=1,default=float)
for k,v in R.items(): print("=====",k); print(json.dumps(v,indent=0,default=float)[:2500])
