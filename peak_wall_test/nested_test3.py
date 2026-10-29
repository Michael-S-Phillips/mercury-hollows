import numpy as np, pandas as pd, json, warnings; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from scipy import stats
ch=pd.read_csv("nested_children_classified.csv"); m=pd.read_csv("morph_all_nested.csv")
ha=pd.read_csv("../code/rerun/latitude/aspect/host_assignments.csv")
g=ha.groupby("host_uid").agg(n_in=("rho",lambda r:int((r<=0.9).sum())))
for d in (ch,m):
    d2=d.merge(g,left_on="uid",right_index=True,how="left"); d["n_in"]=d2.n_in.fillna(0).values
    d["field"]=(d.n_in>=3).astype(int); d["any"]=(d.n_in>=1).astype(int); d["logD"]=np.log10(d.D); d["abslat"]=d.lat.abs()
R={}
def fit(f,data):
    try:
        mo=smf.logit(f,data=data).fit(disp=0,maxiter=300); return {k:[round(float(mo.params[k]),3),round(float(mo.bse[k]),3),round(float(mo.pvalues[k]),4)] for k in mo.params.index}
    except Exception as e: return str(e)
P=ch[(ch.D_p>=60)&(ch.D_p<=300)]
for out in ["field","any"]:
    r={}
    for grp in ["floor_fill_only","floor_through","raised","ambiguous"]:
        s=P[P.group==grp]; r[grp]=dict(N=len(s),k=int(s[out].sum()),rate=round(s[out].mean(),4),hollows=int(s.n_in.sum()),D_med=round(s.D.median(),1),C_med=round(s.C_mean.median(),2))
    # exact comparison through vs fill-only at matched D (10-20 km, where both groups exist)
    a=P[(P.group=="floor_fill_only")&(P.D<20)]; b=P[(P.group=="floor_through")&(P.D<20)]
    r["D10_20_fill_vs_through"]=dict(fill=f"{int(a[out].sum())}/{len(a)}",through=f"{int(b[out].sum())}/{len(b)}",
        fisher_p=round(float(stats.fisher_exact([[a[out].sum(),len(a)-a[out].sum()],[b[out].sum(),len(b)-b[out].sum()]])[1]),4),
        C_fill=round(a.C_mean.median(),2),C_through=round(b.C_mean.median(),2))
    r["logit"]=fit(f"{out} ~ C(group, Treatment('floor_fill_only')) + logD + C_mean + abslat",P[P.group!="ambiguous"])
    r["logit_margin"]=fit(f"{out} ~ margin + logD + C_mean + abslat",P[P.pos=="floor"])
    R[f"parents_60_300_{out}"]=r
# all settings vs isolated, field outcome, D-binned
iso=m[(m.parent<0)&~m.in_orgel&(m.D<=250)]
tab=[]
for lo,hi in zip([10,15,20,30,50],[15,20,30,50,250]):
    row={"D":f"{lo}-{hi}"}
    row["isolated"]=f"{int(iso[(iso.D>=lo)&(iso.D<hi)].field.sum())}/{int(((iso.D>=lo)&(iso.D<hi)).sum())}"
    for grp in ["floor_fill_only","floor_through","raised"]:
        s=P[(P.group==grp)&(P.D>=lo)&(P.D<hi)]; row[grp]=f"{int(s.field.sum())}/{len(s)}"
    tab.append(row)
R["field_rates_by_D_parents60_300"]=tab
json.dump(R,open("nested_test3_record.json","w"),indent=1,default=float)
for k,v in R.items(): print("=====",k); print(json.dumps(v,indent=0,default=float)[:2600])
