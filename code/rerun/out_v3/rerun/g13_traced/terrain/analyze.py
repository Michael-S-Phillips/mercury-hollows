# g13 copy of g11_curvature/terrain/analyze.py (6 Nov 2026). Changed: folder g13_traced/terrain; rows file rows_traced_<CLS>.csv (env CURV_CFG, default traced_v3).
# g11 copy of g10_threshold/analyze.py (6 Nov 2026). Changed: folder g13_traced/terrain/out_<CAT>; G10_MODE removed; the record of each input file is written to peak_wall_record_<file stem>.json (same structure, one key) so configurations can be run in separate calls.
# g10 copy of g7_terrain_fix/analyze.py: output folder g10_threshold/out_<CAT>_<G10_MODE> only.
# g7 copy of g3_catalog/terrain/analyze.py (6 Nov 2026): output folder g7_terrain_fix/out_<CAT>_grid<G7_GRID> only.
# g3_catalog copy of peak_wall_test/analyze.py (6 Nov 2026): only input/output paths and env settings changed.
import pandas as pd, numpy as np, json, sys
from scipy import stats
import os as _os
G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g13_traced/terrain")
CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/out_{CAT}"; _os.makedirs(ODIR, exist_ok=True)
HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))
CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"
         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))
rng=np.random.default_rng(1)
def load(f):
    d=pd.read_csv(f); d=d[(d.D>=30)&(d.D<=250)].copy()
    d["nested"]=(d.parent_D>0)|d.in_orgel_basin
    d["struct"]=np.where(d.cent_str.isin(["cp","mp"]),"peak",np.where(d.cent_str.isin(["pr","rpc","pb"]),"ring","unclassified"))
    return d
def pooled(d,a="peak",b="wall",n="n",nb=2000):
    na,nb_,Aa,Ab=d[f"{n}_{a}"].values,d[f"{n}_{b}"].values,d[f"A_{a}"].values,d[f"A_{b}"].values
    r=lambda i:(na[i].sum()/Aa[i].sum())/(nb_[i].sum()/Ab[i].sum()) if nb_[i].sum()>0 and Aa[i].sum()>0 else np.nan
    idx=np.arange(len(d)); bs=[r(rng.choice(idx,len(idx))) for _ in range(nb)]
    return dict(ratio=round(r(idx),2),ci=[round(np.nanpercentile(bs,2.5),2),round(np.nanpercentile(bs,97.5),2)],
                dens_a=round(1000*na.sum()/Aa.sum(),2),dens_b=round(1000*nb_.sum()/Ab.sum(),2),N=len(d),n_a=int(na.sum()),n_b=int(nb_.sum()))
def crater_level(d,a="peak",b="wall",n="n",minA=(5,20),nb=2000):
    e=d[(d[f"A_{a}"]>=minA[0])&(d[f"A_{b}"]>=minA[1])&((d[f"{n}_{a}"]+d[f"{n}_{b}"])>0)].copy()
    tot=(e[f"{n}_{a}"]+e[f"{n}_{b}"]).values; p=(e[f"A_{a}"]/(e[f"A_{a}"]+e[f"A_{b}"])).values; k=e[f"{n}_{a}"].values
    # presence: craters with >=1 hollow in class a, vs expectation if each crater's hollows fell by area (Poisson-binomial)
    pe=1-(1-p)**tot; obs_pres=int((k>0).sum()); exp_pres=float(pe.sum()); var=float((pe*(1-pe)).sum())
    z=(obs_pres-exp_pres)/np.sqrt(var)
    # count ratio obs/exp with crater bootstrap
    idx=np.arange(len(e)); r=lambda i: k[i].sum()/(tot[i]*p[i]).sum()
    bs=[r(rng.choice(idx,len(idx))) for _ in range(nb)]
    return dict(N=len(e),presence_obs=obs_pres,presence_exp=round(exp_pres,1),presence_z=round(float(z),1),presence_p=float(stats.norm.sf(z)),
                count_obs_over_exp=round(float(r(idx)),2),ci=[round(float(np.percentile(bs,2.5)),2),round(float(np.percentile(bs,97.5)),2)],
                median_area_frac=round(float(np.median(p)),2))
def slope_strat(d):
    num=den=0; out={}
    for b in range(3):
        np_,nw,Ap,Aw=d[f"n_peak_s{b}"].sum(),d[f"n_wall_s{b}"].sum(),d[f"A_peak_s{b}"].sum(),d[f"A_wall_s{b}"].sum()
        out[f"s{b}"]=dict(dens_peak=round(1000*np_/Ap,2) if Ap>0 else None,dens_wall=round(1000*nw/Aw,2) if Aw>0 else None,n_peak=int(np_),n_wall=int(nw),A_peak=round(Ap),A_wall=round(Aw))
        if Ap>0 and Aw>0: T=Ap+Aw; num+=np_*Aw/T; den+=nw*Ap/T
    out["MH_ratio"]=round(num/den,2) if den>0 else None; return out
R={}
for f in sys.argv[1:]:
    d=load(f); key=f.split("/")[-1].replace(".csv","")
    r={"all":dict(pooled_hollows=pooled(d),pooled_fields=pooled(d,n="f"),floor_vs_wall=pooled(d,"floor","wall"),crater_level=crater_level(d),crater_level_fields=crater_level(d,n="f"),
                  eqwall_vs_polewall=pooled(d,"wall_eq","wall_pole"),floor_in_vs_wall=pooled(d,"floor_in","wall"),floor_out_vs_wall=pooled(d,"floor_out","wall"),floor_vs_wall_crater=crater_level(d,"floor","wall"),peak_vs_eqwall=pooled(d,"peak","wall_eq"),peak_vs_eqwall_crater=crater_level(d,"peak","wall_eq",minA=(5,10)),slope_stratified=slope_strat(d))}
    for nm,sub in [("isolated",d[~d.nested]),("nested",d[d.nested]),("central_peak",d[d.struct=="peak"]),("ring",d[d.struct=="ring"]),("unclassified",d[d.struct=="unclassified"]),
                   ("filled_y",d[d.filled=="y"]),("filled_n",d[d.filled=="n"]),("fresh_s",d[d.preserved=="s"]),("degraded_d",d[d.preserved=="d"]),
                   ("D30_60",d[d.D<60]),("D60_120",d[(d.D>=60)&(d.D<120)]),("D120_250",d[d.D>=120]),("isolated_central_peak",d[(~d.nested)&(d.struct=="peak")])]:
        r[nm]=dict(pooled_hollows=pooled(sub),pooled_fields=pooled(sub,n="f"),crater_level=crater_level(sub))
    # leave-one-out: drop the crater with the most peak hollows
    top=d.sort_values("n_peak",ascending=False).head(3)[["uid","name","D","n_peak","n_wall","A_peak","A_wall"]]
    r["top_peak_craters"]=top.to_dict("records"); r["drop_top3"]=pooled(d[~d.uid.isin(top.uid)])
    R[key]=r
for _k in R: json.dump({_k:R[_k]},open(f"{ODIR}/peak_wall_record_{_k}.json","w"),indent=1,default=str)   # g11: one record per input file
k=list(R)[0]
for key,r in R.items():
    print("=====",key); a=r["all"]
    for kk in ["pooled_hollows","pooled_fields","floor_vs_wall","floor_in_vs_wall","floor_out_vs_wall","floor_vs_wall_crater","crater_level","crater_level_fields","eqwall_vs_polewall","peak_vs_eqwall","peak_vs_eqwall_crater"]: print(" ",kk,a[kk])
    print("  slope",a["slope_stratified"])
    for nm in ["isolated","nested","central_peak","ring","unclassified","filled_y","filled_n","fresh_s","degraded_d","D30_60","D60_120","D120_250","isolated_central_peak"]:
        print(f"  {nm:22s} pooled {r[nm]['pooled_hollows']}  fields {r[nm]['pooled_fields']['ratio']} {r[nm]['pooled_fields']['ci']}  crater {r[nm]['crater_level']}")
    print("  top",r["top_peak_craters"]); print("  drop_top3",r["drop_top3"])
