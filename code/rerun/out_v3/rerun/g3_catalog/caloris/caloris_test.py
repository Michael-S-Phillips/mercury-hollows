# g3_catalog copy of peak_wall_test/caloris_test.py (6 Nov 2026): only input/output paths and env settings changed.
import numpy as np, pandas as pd, json, warnings; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")
CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/caloris/out_{CAT}"; _os.makedirs(ODIR, exist_ok=True)
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/insolation_thermal/data/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
ch=pd.read_csv(_os.path.expanduser("~/mnt/hollows-emissivity/peak_wall_test/nested_children_classified.csv")); m=pd.read_csv(_os.path.expanduser("~/mnt/hollows-emissivity/peak_wall_test/morph_all_nested.csv"))
ha=pd.read_csv(HOSTS)
g=ha.groupby("host_uid").agg(n_in=("rho",lambda r:int((r<=0.9).sum())))
for d in (ch,m):
    d2=d.merge(g,left_on="uid",right_index=True,how="left"); d["n_in"]=d2.n_in.fillna(0).values
    d["field"]=(d.n_in>=3).astype(int); d["any"]=(d.n_in>=1).astype(int); d["logD"]=np.log10(d.D); d["abslat"]=d.lat.abs()
d_ex=lambda D: 0.1*12**0.15*D**0.85
Dthr={t:(t/(0.1*12**0.15))**(1/0.85) for t in (2.5,3.0,3.5)}
cal=ch[ch.name_p=="Caloris"].copy(); cal["set"]="caloris"
iso=m[(m.parent<0)&~m.in_orgel&(m.D<=250)].copy(); iso["set"]="isolated"
R={"Ernst2015_plains_km":[2.5,3.5],"D_for_excavation_depth":{str(k):round(v,1) for k,v in Dthr.items()}}
bins=[10,15,20,25,30,35,40,50,70,250]
for cmin in (0.0,0.2):
    rows=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        a=cal[(cal.D>=lo)&(cal.D<hi)&(cal.C_mean>=cmin)]; b=iso[(iso.D>=lo)&(iso.D<hi)&(iso.C_mean>=cmin)]
        rows.append(dict(D=f"{lo}-{hi}",caloris_field=f"{int(a.field.sum())}/{len(a)}",caloris_any=f"{int(a['any'].sum())}/{len(a)}",
                         iso_field=f"{int(b.field.sum())}/{len(b)}",iso_any=f"{int(b['any'].sum())}/{len(b)}",
                         C_cal=round(a.C_mean.median(),2) if len(a) else None,C_iso=round(b.C_mean.median(),2) if len(b) else None))
    R[f"rates_Cmin{cmin}"]=rows
A=pd.concat([cal,iso])[["field","any","logD","C_mean","abslat","set","D"]]
A=A[(A.D>=10)&(A.D<=100)]
for out in ("field","any"):
    mo=smf.logit(f"{out} ~ logD*C(set, Treatment('isolated')) + C_mean + abslat",data=A).fit(disp=0,maxiter=300)
    R[f"logit_{out}"]={k:[round(float(mo.params[k]),3),round(float(mo.bse[k]),3),round(float(mo.pvalues[k]),5)] for k in mo.params.index}
    # predicted probability at median C for D = 15, 25, 35, 45 km
    pr={}
    for s in ("caloris","isolated"):
        X=pd.DataFrame({"logD":np.log10([15,25,35,45]),"C_mean":0.3,"abslat":32,"set":s}); pr[s]=[round(float(p),4) for p in mo.predict(X)]
    R[f"pred_{out}_C0.3"]=pr
# step test within Caloris: excavation depth above vs below 2.5 km (D ~ 28 km)
cal["through"]=(d_ex(cal.D)>=2.5).astype(int)
mo=smf.logit("field ~ through + C_mean",data=cal).fit(disp=0); R["caloris_step_field"]={k:[round(float(mo.params[k]),3),round(float(mo.bse[k]),3),round(float(mo.pvalues[k]),6)] for k in mo.params.index}
R["caloris_counts"]={"d_ex<2.5":f"{int(cal[cal.through==0].field.sum())}/{int((cal.through==0).sum())}","d_ex>=2.5":f"{int(cal[cal.through==1].field.sum())}/{int((cal.through==1).sum())}"}
# edge vs center (Ernst: thinner fill toward the edge) for 20-40 km
for lo,hi in [(20,40)]:
    s=cal[(cal.D>=lo)&(cal.D<hi)]
    R["caloris_20_40_edge_vs_center"]={"rho_p>=0.7":f"{int(s[s.rho_p>=0.7].field.sum())}/{int((s.rho_p>=0.7).sum())}","rho_p<0.7":f"{int(s[s.rho_p<0.7].field.sum())}/{int((s.rho_p<0.7).sum())}"}
cal.to_csv(f"{ODIR}/caloris_children.csv",index=False)
json.dump(R,open(f"{ODIR}/caloris_test_record.json","w"),indent=1)
print(json.dumps(R,indent=0))
