#!/usr/bin/env python3
"""Hosting around Caloris by azimuth and distance (5 Nov 2026): is the hollow-rich province NW of Caloris part of a
ring of enhanced hosting in the Caloris exterior (Van Eyck lineated terrain, up to ~1,000 km beyond the rim), or is it
specific to the NW sector (where b59 lies)? Craters >= 5 km outside Caloris (rho 1-3 Caloris radii, Orgel et al. 2020
center 29.9N 161.4E, D 1,556 km), logistic hosting (h09) with the paper's controls, against craters beyond 3 radii.
Writes caloris_sectors.json."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(f"{HERE}/../../../.."); R_M = 2439.4
f = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})
la0, lo0, Rc = np.radians(29.91), np.radians(161.43), 1555.9 / 2
la, lo = np.radians(f.lat.values), np.radians(f.lon.values % 360)
c = np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0)
f["rc"] = np.arccos(np.clip(c, -1, 1)) * R_M / Rc
az = np.degrees(np.arctan2(np.sin(lo - lo0) * np.cos(la), np.cos(la0) * np.sin(la) - np.sin(la0) * np.cos(la) * np.cos(lo - lo0))) % 360
f["az"] = az
sect = pd.cut((f.az + 22.5) % 360, np.arange(0, 361, 45), labels=["N", "NE", "E", "SE", "S", "SW", "W", "NW"], right=False)
f["sect"] = sect.astype(str)
f["zone"] = np.where(f.rc < 1, "inside", np.where(f.rc < 2.3, "ring", np.where(f.rc < 3, "mid", "far")))
f["grp"] = np.where(f.zone == "ring", "ring_" + f.sect, np.where(f.zone == "far", "far", "other"))
d = f[f.zone.isin(["ring", "far"])].copy(); d["grp"] = pd.Categorical(d.grp, categories=["far"] + [f"ring_{s}" for s in ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]])
m = smf.logit("h09 ~ C(grp) + logD + abslat + C1 + nsp + dT", d).fit(disp=0)
R = {"ring_def": "1.0-2.3 Caloris radii (rim to ~1,000 km beyond it)", "reference": "craters beyond 3 Caloris radii", "sectors": {}}
for s in ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]:
    t = f"C(grp)[T.ring_{s}]"; g = d[d.grp == f"ring_{s}"]
    R["sectors"][s] = {"craters": int(len(g)), "hosts": int(g.h09.sum()), "OR": [round(float(np.exp(m.params[t])), 2), round(float(np.exp(m.conf_int().loc[t][0])), 2), round(float(np.exp(m.conf_int().loc[t][1])), 2)]}
# b59 position relative to Caloris
b_la, b_lo = np.radians(50.05), np.radians(119.26)
R["b59_center_rc"] = round(float(np.arccos(np.sin(la0)*np.sin(b_la)+np.cos(la0)*np.cos(b_la)*np.cos(b_lo-lo0)) * R_M / Rc), 2)
R["b59_center_az"] = round(float(np.degrees(np.arctan2(np.sin(b_lo-lo0)*np.cos(b_la), np.cos(la0)*np.sin(b_la)-np.sin(la0)*np.cos(b_la)*np.cos(b_lo-lo0))) % 360), 1)
# NW sector: does the excess sit inside/near b59 or along the whole sector?
nw = d[d.grp == "ring_NW"].copy()
bc = np.arccos(np.clip(np.sin(b_la)*np.sin(np.radians(nw.lat))+np.cos(b_la)*np.cos(np.radians(nw.lat))*np.cos(np.radians(nw.lon % 360)-b_lo), -1, 1)) * R_M / 350.3
nw["near59"] = bc < 1.5
R["NW_ring_near_b59"] = [int(nw.near59.sum()), int(nw[nw.near59].h09.sum())]; R["NW_ring_away_b59"] = [int((~nw.near59).sum()), int(nw[~nw.near59].h09.sum())]
json.dump(R, open(f"{HERE}/caloris_sectors.json", "w"), indent=1); print(json.dumps(R, indent=1))
