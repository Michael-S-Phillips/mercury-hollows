import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../../.."))
"""
Incidence-controlled completeness null, STEP 2.
Fit the per-frame detection function separately on hollows whose canonical
detection is Stage 1 vs Stage >= 2 in the Bickel 0.8 release (degraded /
relief-only), and report the incidence term b_inc for each with a bootstrap CI.

b_inc enters pdet as exp(b_inc*(inc-60)/25): positive => detection RISES with
incidence (relief-only hollows need long shadows), which the v2 completeness
model (b_inc = -0.10, fitted on halo-bearing repeat detections) does not capture.

Match canonical detection to 0.8 Stage by thumbnail tag, else by position
within 300 m (0.8 positions are UNCORRECTED -> matched against canonical
release coordinates lat_release/lon_release). American English throughout.
"""
import os, json, csv, io
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.special import erf
from scipy.optimize import minimize

RM = 2439.7
PJ = _REPO + "/"
OUT = PJ + "code/rerun/latitude/incidence_null/"
SEED = 20260911
rng = np.random.default_rng(SEED)

def xyz(lat, lon):
    return np.c_[np.cos(np.radians(lat))*np.cos(np.radians(lon)),
                 np.cos(np.radians(lat))*np.sin(np.radians(lon)),
                 np.sin(np.radians(lat))]

# ---------------------------------------------------------------- load
HOL = pd.read_csv(PJ+"hollow-networks/data_v2/hollows_canonical.csv")
H = HOL[~HOL.is_crater_trace.astype(bool)].reset_index(drop=True)  # h index space
CV = pd.read_csv(PJ+"hollow-networks/data_v2/completeness_pairs.csv")

# 0.8 release stages
raw = open(PJ+"data/bickel-2024/Bickel_et_al_2024_MDIS_hollows_08_v1.0.csv",'rb').read().decode('latin-1')
rows = list(csv.reader(io.StringIO(raw), delimiter=';'))
hdr = [c.strip() for c in rows[0]]; ix={c:i for i,c in enumerate(hdr)}
rec=[]
for r in rows[1:]:
    if len(r)<len(hdr): continue
    try:
        rec.append(dict(thumb=r[ix['Thumbnail tag']].strip(), stage=int(float(r[ix['Stage']])),
                        lat=float(r[ix['Detection center LAT']]),
                        lon=float(r[ix['Detection center LON 180d']])%360.0))
    except (ValueError,IndexError): continue
B08=pd.DataFrame(rec)

# ---------------------------------------------------------------- match
tag2stage=dict(zip(B08.thumb,B08.stage))
H['stage08']=H.thumb.map(tag2stage)
n_tag=int(H.stage08.notna().sum())
un=H.stage08.isna().values
tree=cKDTree(xyz(B08.lat.values,B08.lon.values))
d,i=tree.query(xyz(H.lat_release.values[un],H.lon_release.values[un]),k=1)
gc=2*RM*np.arcsin(np.clip(d/2,0,1))
m=gc<=0.3
stfb=np.where(m,B08.stage.values[i],np.nan)
tmp=H.loc[un,'stage08'].copy(); tmp.iloc[:]=stfb; H.loc[un,'stage08']=tmp.values
n_pos=int(m.sum())
H['stage_grp']=np.where(H.stage08>=2,'ge2',np.where(H.stage08==1,'s1',None))
CV['grp']=CV.h.map(H['stage_grp'].to_dict())

# ---------------------------------------------------------------- model
def Phi(z): return 0.5*(1+erf(z/np.sqrt(2)))
def pdet(x,inc_v,th):
    pmax,lx50,ls1,lxhi,ls2,binc=th
    p=pmax*Phi((np.log(x)-lx50)/np.exp(ls1))*(1-Phi((np.log(x)-lxhi)/np.exp(ls2)))
    return np.clip(p*np.exp(binc*(inc_v-60)/25),1e-6,0.999)

def make_blocks(sub):
    """Return list of (X,I,Y) arrays, one per hollow (contiguous blocks)."""
    sub=sub.sort_values('h')
    blocks=[]
    for _,g in sub.groupby('h'):
        blocks.append((g.D_px.values, g.inc.values, g.detected.values.astype(float)))
    return blocks

def nll_from_blocks(th, X, I, Y, starts):
    if not (0.005<th[0]<0.95) or th[1]>th[3]: return 1e18
    p=pdet(X,I,th)
    lp=Y*np.log(p)+(1-Y)*np.log(1-p); l1p=np.log(1-p)
    S1=np.add.reduceat(lp,starts); S0=np.add.reduceat(l1p,starts)
    return -np.sum(S1-np.log(1-np.exp(S0)+1e-300))

def assemble(blocks):
    X=np.concatenate([b[0] for b in blocks]); I=np.concatenate([b[1] for b in blocks])
    Y=np.concatenate([b[2] for b in blocks])
    lens=np.array([len(b[0]) for b in blocks]); starts=np.r_[0,np.cumsum(lens)[:-1]]
    return X,I,Y,starts

def fit(blocks, multistart=True, th0=None):
    X,I,Y,starts=assemble(blocks)
    if multistart:
        grid=[[p0,np.log(x0),np.log(0.5),np.log(170),np.log(0.35),0.0]
              for p0 in [0.05,0.1,0.2] for x0 in [12,18,25]]
    else:
        grid=[th0]
    best=None
    for g0 in grid:
        r=minimize(nll_from_blocks,g0,args=(X,I,Y,starts),method='Nelder-Mead',
                   options=dict(maxiter=10000,xatol=1e-5,fatol=1e-3))
        if best is None or r.fun<best.fun: best=r
    return best.x

NBOOT=600
res={}
blocks_by_grp={}
for g in ['s1','ge2']:
    blk=make_blocks(CV[CV.grp==g]); blocks_by_grp[g]=blk
    th=fit(blk, multistart=True)
    n_h=len(blk); n_pairs=sum(len(b[0]) for b in blk); n_det=int(sum(b[2].sum() for b in blk))
    # bootstrap over hollows
    bs=[]
    for _ in range(NBOOT):
        idx=rng.integers(0,n_h,n_h)
        thb=fit([blk[k] for k in idx], multistart=False, th0=th)
        bs.append(thb[5])
    bs=np.array(bs)
    ci=np.percentile(bs,[2.5,97.5])
    res[g]=dict(theta=dict(p_max=float(th[0]),x50_px=float(np.exp(th[1])),sigma_lo=float(np.exp(th[2])),
                           x_hi_px=float(np.exp(th[3])),sigma_hi=float(np.exp(th[4])),b_inc=float(th[5])),
                b_inc=float(th[5]), b_inc_boot_mean=float(bs.mean()), b_inc_boot_sd=float(bs.std(ddof=1)),
                b_inc_CI95=[float(ci[0]),float(ci[1])],
                n_hollows=int(n_h), n_pairs=int(n_pairs), n_detected_pairs=n_det,
                det_inc_median=float(CV[(CV.grp==g)&CV.detected].inc.median()))
    print(f"{g}: b_inc={th[5]:.3f} CI95=[{ci[0]:.3f},{ci[1]:.3f}] n_hollows={n_h} det_inc_med={res[g]['det_inc_median']:.0f}")

# difference in b_inc with bootstrap (paired resample would need same draws; report independent)
rng2=np.random.default_rng(SEED+1)
diff=[]
for _ in range(NBOOT):
    i1=rng2.integers(0,len(blocks_by_grp['ge2']),len(blocks_by_grp['ge2']))
    i0=rng2.integers(0,len(blocks_by_grp['s1']),len(blocks_by_grp['s1']))
    b1=fit([blocks_by_grp['ge2'][k] for k in i1],multistart=False,th0=[res['ge2']['theta'][q] for q in ['p_max']]+[np.log(res['ge2']['theta']['x50_px']),np.log(res['ge2']['theta']['sigma_lo']),np.log(res['ge2']['theta']['x_hi_px']),np.log(res['ge2']['theta']['sigma_hi']),res['ge2']['theta']['b_inc']])
    b0=fit([blocks_by_grp['s1'][k] for k in i0],multistart=False,th0=[res['s1']['theta'][q] for q in ['p_max']]+[np.log(res['s1']['theta']['x50_px']),np.log(res['s1']['theta']['sigma_lo']),np.log(res['s1']['theta']['x_hi_px']),np.log(res['s1']['theta']['sigma_hi']),res['s1']['theta']['b_inc']])
    diff.append(b1[5]-b0[5])
diff=np.array(diff)
dci=np.percentile(diff,[2.5,97.5])
p_ge=float((diff<=0).mean())
res['b_inc_diff_ge2_minus_s1']=dict(point=float(res['ge2']['b_inc']-res['s1']['b_inc']),
    boot_mean=float(diff.mean()),CI95=[float(dci[0]),float(dci[1])],
    frac_bootstrap_le0=p_ge)
print("diff ge2-s1:", round(res['ge2']['b_inc']-res['s1']['b_inc'],3),
      "CI95=[%.3f,%.3f]"%(dci[0],dci[1]), "frac<=0=%.3f"%p_ge)

out=dict(n_labeled=int(H.stage_grp.notna().sum()), n_tag_matches=n_tag, n_pos_matches=n_pos,
         n_total_canonical=int(len(H)), nboot=NBOOT, seed=SEED,
         match="thumbnail tag or position <=300 m (release coords)",
         v2_b_inc=-0.10, groups=res)
json.dump(out, open(OUT+"binc_by_stage.json","w"), indent=1)
H[['thumb','lat','lon','stage08','stage_grp']].to_csv(OUT+"canonical_stage_labels.csv", index=False)
print("saved binc_by_stage.json ; labeled", int(H.stage_grp.notna().sum()))
