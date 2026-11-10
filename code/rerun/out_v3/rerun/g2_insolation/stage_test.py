#!/usr/bin/env python3
# g2 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/insolation_v6/stage_test.py (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).
"""Task (d), 2 Nov 2026: do Stage 1 hollows (bright, sharp; Deutsch et al. 2025) favor the hot longitudes and
equator-facing walls more than Stage 2 hollows?

Stage labels: the Deutsch et al. (2025) stages carried in the Bickel et al. (2024) v1.0 confidence >= 0.8 release
(data/bickel-2024/Bickel_et_al_2024_MDIS_hollows_08_v1.0.csv; 3,268 detections: Stage 1 2,129, Stage 2 1,129, Stage 3 10),
matched to the project's canonical (deduplicated, re-georeferenced) non-trace hollows as in
code/rerun/latitude/incidence_null/fit_binc_by_stage.py (thumbnail tag, else nearest 0.8 detection within 300 m of the
release position); canonical_stage_labels.csv there holds the result (1,500 Stage 1, 968 Stage 2, 3 Stage 3). The
incidence and pixel scale of the labeled 0.8 detection are carried along, because stage labels depend on incidence
(stage_vs_incidence.json: odds of Stage >= 2 double per 10 deg of incidence; repeat-image labels agree 65%).

Models: GLM binomial for P(Stage 1 | Stage 1 or 2), standard errors clustered by hollow field (field_id).
Controls: incidence (per 10 deg), log10 pixel scale, emission, |lat|. Insolation terms: Delta T_w (per 10 K), hot sector,
equator-facing wall (north, 0.55-1.05 radii of hosts >= 15 km, as in Figure S5c). Writes stage_test.json.
"""
import json, csv, io, sys, numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy.spatial import cKDTree
from common_null import *
sys.path.insert(0, f"{ROOT}/manuscripts/GRL_insolation_hollows/scripts"); import thermal
G = f"{ROOT}/manuscripts/GRL_insolation_hollows/data"

H = pd.read_csv(IN + "/hollows_canonical.csv"); H = H[~H.is_crater_trace.astype(bool)].reset_index(drop=True)
L = pd.read_csv(f"{IN}/canonical_stage_labels.csv")
assert (H.thumb.values == L.thumb.values).all(); H["stage08"] = L.stage08.values

# 0.8 release rows, to carry the incidence / pixel scale / emission of the labeled detection
raw = open(f"{ROOT}/data/bickel-2024/Bickel_et_al_2024_MDIS_hollows_08_v1.0.csv", "rb").read().decode("latin-1")
rows = list(csv.reader(io.StringIO(raw), delimiter=';')); hdr = [c.strip() for c in rows[0]]; ix = {c: i for i, c in enumerate(hdr)}
rec = []
for r in rows[1:]:
    if len(r) < len(hdr): continue
    try:
        rec.append(dict(thumb=r[ix['Thumbnail tag']].strip(), stage=int(float(r[ix['Stage']])), lat=float(r[ix['Detection center LAT']]),
                        lon=float(r[ix['Detection center LON 180d']]) % 360.0, inc8=float(r[ix['Incidence (°)']]),
                        emi8=float(r[ix['Emission (°)']]), px8=float(r[ix['Pixel scale (m)']])))
    except (ValueError, IndexError): continue
B = pd.DataFrame(rec)
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
bt = B.drop_duplicates("thumb").set_index("thumb")
j = H.thumb.map(lambda t: bt.index.get_loc(t) if t in bt.index else -1).values
tree = cKDTree(xyz(B.lat.values, B.lon.values)); d, i = tree.query(xyz(H.lat_release.values, H.lon_release.values), k=1)
gc = 2 * R_MERCURY * np.arcsin(np.clip(d / 2, 0, 1))
for col in ("stage", "inc8", "emi8", "px8"):
    v_tag = np.where(j >= 0, bt[col].values[np.maximum(j, 0)], np.nan)
    v_pos = np.where(gc <= 0.3, B[col].values[i], np.nan)
    H[col + "_m"] = np.where(np.isfinite(v_tag), v_tag, v_pos)
chk = H.stage08.notna(); agree = float((H.loc[chk, "stage_m"] == H.loc[chk, "stage08"]).mean())
S = H[H.stage08.isin([1, 2])].copy(); S["s1"] = (S.stage08 == 1).astype(int)
S["inc10"] = S.inc8_m / 10; S["logpx"] = np.log10(S.px8_m); S["emi"] = S.emi8_m; S["abslat"] = S.lat.abs()
_z = np.load(f"{G}/thermal_grid.npz")
S["dT"] = thermal.lookup(_z["la"], _z["lo"], _z["G"], S.lat.values, S.lon.values) - thermal.lookup(_z["la"], _z["lo"], _z["G"], S.lat.values, np.full(len(S), 90.0))
S["dT10"] = S.dT / 10; S["hot"] = hot_mask(S.lon.values).astype(int)
S["dcal"] = gc_dist_km(CALORIS[0], CALORIS[1], S.lat.values, S.lon.values)
A = pd.read_csv(f"{G2_OUT}/aspect_wall_hollows.csv")[["hollow_id", "host_uid", "host_D_km", "rho", "hemi", "poleward"]]
S = S.merge(A, on="hollow_id", how="left")
OUT = {"description": __doc__, "n_labeled_stage1_2": int(len(S)), "n_stage1": int(S.s1.sum()), "n_stage2": int((1 - S.s1).sum()),
       "label_rematch_agreement": agree, "n_fields": int(S.field_id.nunique())}

def glm(f, d, term):
    m = smf.glm(f, data=d, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": d.field_id.values})
    b, se = m.params[term], m.bse[term]
    return dict(n=int(len(d)), n_s1=int(d.s1.sum()), n_fields=int(d.field_id.nunique()), OR=float(np.exp(b)),
                lo=float(np.exp(b - 1.96 * se)), hi=float(np.exp(b + 1.96 * se)), p=float(m.pvalues[term]))
CTRL = "inc10 + logpx + emi + abslat"
R = {}
subsets = {"all": S, "lat_lt15": S[S.abslat < 15], "15_75N": S[(S.lat >= 15) & (S.lat < 75)],
           "lat_lt60": S[S.abslat <= 60], "lat_lt60_noCaloris": S[(S.abslat <= 60) & (S.dcal > CALORIS_KM)]}
for nm, d in subsets.items():
    R[nm] = dict(dT10=glm(f"s1 ~ dT10 + {CTRL}", d, "dT10"), hot=glm(f"s1 ~ hot + {CTRL}", d, "hot"),
                 hot_nocontrols=glm("s1 ~ hot", d, "hot"),
                 raw=dict(s1_hot=int(d[d.hot == 1].s1.sum()), s2_hot=int((1 - d[d.hot == 1].s1).sum()),
                          s1_warm=int(d[d.hot == 0].s1.sum()), s2_warm=int((1 - d[d.hot == 0].s1).sum()),
                          mean_inc_hot=float(d[d.hot == 1].inc8_m.mean()), mean_inc_warm=float(d[d.hot == 0].inc8_m.mean())))
# Delta T_w threshold look: Stage 1 share by Delta T_w bin (raw), and step vs linear (AIC, ordinary GLM)
d60 = S[S.abslat <= 60].copy(); e = pd.cut(d60.dT, [0, 20, 30, 40, 50, 55, 60, 63])
R["s1_share_by_dT"] = d60.groupby(e, observed=True).agg(n=("s1", "size"), s1=("s1", "mean"), inc=("inc8_m", "median")).reset_index().astype({"dT": str}).to_dict("records")
aic = {"base": smf.glm(f"s1 ~ {CTRL}", data=d60, family=sm.families.Binomial()).fit().aic,
       "linear": smf.glm(f"s1 ~ dT10 + {CTRL}", data=d60, family=sm.families.Binomial()).fit().aic}
st = {}
for t in range(20, 61, 2):
    d60["st"] = (d60.dT >= t).astype(float); st[t] = smf.glm(f"s1 ~ st + {CTRL}", data=d60, family=sm.families.Binomial()).fit().aic
tb = min(st, key=st.get); aic["step_best"] = st[tb] + 2; R["s1_step_vs_linear_lat60"] = dict(AIC={k: round(v, 1) for k, v in aic.items()}, best_t=tb)
d60["st"] = (d60.dT >= tb).astype(float); dnc = d60[d60.dcal > CALORIS_KM].copy()
R["s1_step_vs_linear_lat60"]["step_OR_clustered"] = glm(f"s1 ~ st + {CTRL}", d60, "st")
R["s1_step_vs_linear_lat60"]["step_OR_clustered_noCaloris"] = glm(f"s1 ~ st + {CTRL}", dnc, "st")
R["s1_step_vs_linear_lat60"]["step_OR_clustered_inc_lt60"] = glm(f"s1 ~ st + {CTRL}", d60[d60.inc8_m < 60], "st")
R["s1_step_vs_linear_lat60"]["n_fields_above_t"] = int(d60[d60.st == 1].field_id.nunique())
# equator-facing walls (north, hosts >= 15 km, 0.55-1.05 radii)
w = S[(S.rho > .55) & (S.rho < 1.05) & (S.host_D_km >= 15) & (S.hemi == "N")].copy(); w["eqf"] = w.poleward.astype(int)
R["walls_north"] = {"all": glm(f"s1 ~ eqf + {CTRL}", w, "eqf"), "lat_ge30": glm(f"s1 ~ eqf + {CTRL}", w[w.lat >= 30], "eqf"),
                    "lat_ge45": glm(f"s1 ~ eqf + {CTRL}", w[w.lat >= 45], "eqf"), "lat_lt30": glm(f"s1 ~ eqf + {CTRL}", w[w.lat < 30], "eqf"),
                    "eqfacing_frac_s1": float(w[w.s1 == 1].eqf.mean()), "eqfacing_frac_s2": float(w[w.s1 == 0].eqf.mean()),
                    "n_s1": int(w.s1.sum()), "n_s2": int((1 - w.s1).sum()),
                    "eqfacing_frac_by_band": {f"{lo}-{hi}": dict(s1=float(w[(w.lat >= lo) & (w.lat < hi) & (w.s1 == 1)].eqf.mean()),
                                                                 s2=float(w[(w.lat >= lo) & (w.lat < hi) & (w.s1 == 0)].eqf.mean()),
                                                                 n1=int(w[(w.lat >= lo) & (w.lat < hi)].s1.sum()), n2=int((1 - w[(w.lat >= lo) & (w.lat < hi)].s1).sum()))
                                              for lo, hi in ((0, 30), (30, 45), (45, 60), (60, 75))}}
OUT["results"] = R
json.dump(OUT, open(f"{G2_OUT}/stage_test.json", "w"), indent=1, default=float)
print({k: OUT[k] for k in ("n_labeled_stage1_2", "n_stage1", "n_stage2", "label_rematch_agreement", "n_fields")})
for nm in subsets:
    r = R[nm]; print(f"== {nm}: dT10 OR {r['dT10']['OR']:.2f} ({r['dT10']['lo']:.2f}-{r['dT10']['hi']:.2f}) p={r['dT10']['p']:.2g} n={r['dT10']['n']}"
                     f" | hot OR {r['hot']['OR']:.2f} ({r['hot']['lo']:.2f}-{r['hot']['hi']:.2f}) p={r['hot']['p']:.2g} | no controls {r['hot_nocontrols']['OR']:.2f}"
                     f" ({r['hot_nocontrols']['lo']:.2f}-{r['hot_nocontrols']['hi']:.2f}) | raw {r['raw']}")
print("s1 share by dT:", [(x['dT'], x['n'], round(x['s1'], 2), round(x['inc'], 1)) for x in R["s1_share_by_dT"]])
print("step vs linear:", R["s1_step_vs_linear_lat60"])
for k, v in R["walls_north"].items(): print("walls", k, v)
