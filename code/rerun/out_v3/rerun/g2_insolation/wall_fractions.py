#!/usr/bin/env python3
"""g2 (6 Nov 2026), new script: the wall-orientation numbers of Table S6 / Figure S5c. The per-hollow and per-host-crater
fractions of northern wall hollows on equator-facing walls, with 95% bootstrap intervals, computed exactly as panel (c) of
figscripts/make_figs_pr_insolation_v9.py (same wall definition, Caloris center, seed 20260911, 4,000 resamples, same order
of calls), which draws the panel but does not save the numbers. Reads out/<cat>/aspect_wall_hollows.csv (aspect_stats.py).
Writes out/<cat>/wall_fractions.json."""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CAT = os.environ.get("G2_CAT", "v3"); OUTD = f"{HERE}/out/{CAT}"
R_M = 2439.4      # as in the figure script (R = 2439400 m)
a = pd.read_csv(f"{OUTD}/aspect_wall_hollows.csv")
w = a[(a.rho > .55) & (a.rho < 1.05) & (a.host_D_km >= 15)].copy()
cal = (31.9, 161.9)
def gc(lat1, lon1, lat2, lon2):
    p1, p2, l1, l2 = map(np.radians, (lat1, lat2, lon1, lon2))
    return R_M * np.arccos(np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(l1 - l2), -1, 1))
w["dcal"] = gc(w.lat, w.lon, *cal)
rng = np.random.default_rng(20260911)
def stats(sub, nb=4000):
    k = sub.poleward.astype(bool); f = k.mean()
    bt = [rng.choice(k.values, len(k)).mean() for _ in range(nb)]
    pc = sub.groupby("host_uid").poleward.mean().values
    bc = [rng.choice(pc, len(pc)).mean() for _ in range(nb)]
    return dict(n=int(len(sub)), nc=int(len(pc)), f=float(f), lo=float(np.percentile(bt, 2.5)), hi=float(np.percentile(bt, 97.5)),
                fc=float(pc.mean()), clo=float(np.percentile(bc, 2.5)), chi=float(np.percentile(bc, 97.5)))
N = w[w.hemi == "N"]
bands = [(0, 15), (15, 30), (30, 45), (45, 60), (60, 75)]
ST = {f"{lo}-{hi}": stats(N[(N.abslat >= lo) & (N.abslat < hi)]) for lo, hi in bands}
ST["15-30_noCaloris"] = stats(N[(N.abslat >= 15) & (N.abslat < 30) & (N.dcal >= 1200)])
ST["south_all"] = stats(w[w.hemi == "S"])
ST["n_wall_hollows_all"] = int(len(w))
json.dump(ST, open(f"{OUTD}/wall_fractions.json", "w"), indent=1)
for k, v in ST.items(): print(k, v)
