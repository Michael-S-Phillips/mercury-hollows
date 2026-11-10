#!/usr/bin/env python3
"""g17 (7 Nov 2026), new script: craters for the second tracing round (Michael, 7 Nov: trace the next 50 hosts and the peaks of
degraded craters).
Group B, hosts round 2: the 50 untraced terrain-test hosts (g13 rows_traced_r2_v3.csv, class_source curv_v3, 30-250 km) with the most
  hollows (n_peak + n_wall + n_floor + n_other), excluding the two hosts Michael marked "can't tell" in round 1 (110709, 111168);
  ties broken by diameter, larger first. Full task: peaks, wall foot, slump.
Group C, degraded hosts: every Kinczyk et al. (2020) Class 1-2 crater with assigned hollows (g13 knobs/class12_peaks_traced_r2.csv,
  n_assigned > 0) that is not traced yet and not in group B. Peaks only. Needed so that every hollow on a Class 1-2 peak is counted.
Group D, degraded sample: 60 Class 1-2 craters without assigned hollows, drawn without replacement with probability proportional to
  the completeness-weighted area within 0.15 R of the center (A_p015; systematic PPS on a random order, seed 7; craters whose
  probability would reach 1 are taken with certainty). Peaks only. Their
  traced peak area, divided by the inclusion probability, estimates the peak area of all 2,052 Class 1-2 craters without hollows
  (Horvitz-Thompson); with curvature areas as a stand-in, the relative standard error of that total is 0.23 (5 certainty craters), and
  craters with A_p015 = 0 (no sampling chance) hold 4% of the curvature peak area of craters without hollows.
Writes selection.json."""
import os, json, numpy as np, pandas as pd
RR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = pd.read_csv(f"{RR}/g13_traced/terrain/out_v3/rows_traced_r2_v3.csv"); rows["n_tot"] = rows[["n_peak", "n_wall", "n_floor", "n_other"]].sum(axis=1)
cant = {110709, 111168}
b = rows[(rows.class_source != "traced") & ~rows.uid.isin(cant)].sort_values(["n_tot", "D"], ascending=[False, False]).head(50)
k = pd.read_csv(f"{RR}/g13_traced/knobs/class12_peaks_traced_r2.csv")
c = k[(k.n_assigned > 0) & ~k.traced & ~k.uid.isin(b.uid)].sort_values("n_assigned", ascending=False)
nh = k[k.n_assigned == 0].reset_index(drop=True); x = nh.A_p015.clip(lower=0).values; n = 60
rng = np.random.default_rng(7); perm = rng.permutation(len(nh)); xp = x[perm]
cert = np.zeros(len(xp), bool)                  # certainty units: craters whose share would give an inclusion probability >= 1
while True:
    m = n - cert.sum(); X = xp[~cert].sum(); new = (~cert) & (m * xp / X >= 1)
    if not new.any(): break
    cert |= new
pi = np.where(cert, 1.0, m * xp / X)
cum = np.cumsum(np.where(cert, 0, pi)); start = rng.uniform(0, 1); picks = np.searchsorted(cum, start + np.arange(m))
idx = np.concatenate([np.where(cert)[0], picks]); sel = nh.iloc[perm[idx]].copy(); sel["pi"] = pi[idx]
assert sel.uid.nunique() == n, sel.uid.nunique()
out = dict(created="2026-11-07", note=__doc__,
    group_B=[dict(uid=int(r.uid), name=str(r["name"]).strip() if isinstance(r["name"], str) else "", D=round(float(r.D), 2), n_tot=int(r.n_tot)) for _, r in b.iterrows()],
    group_C=[dict(uid=int(r.uid), D=round(float(r.D), 2), kclass=int(r.kclass), n_assigned=int(r.n_assigned)) for _, r in c.iterrows()],
    group_D=[dict(uid=int(r.uid), D=round(float(r.D), 2), kclass=int(r.kclass), A_p015=round(float(r.A_p015), 3), pi=round(float(r.pi), 6)) for _, r in sel.iterrows()],
    design=dict(N_nonhost=int(len(nh)), n=n, n_certainty=int(cert.sum()), X_noncertainty=float(X), class12_in_B=int(k[(k.n_assigned > 0) & ~k.traced & k.uid.isin(b.uid)].shape[0])))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "selection.json"), "w"), indent=1)
cov = rows[rows.class_source == "traced"].n_tot.sum() + b.n_tot.sum()
print("B", len(b), "hollows", int(b.n_tot.sum()), "min n_tot", int(b.n_tot.min()), "| terrain coverage with B", int(cov), "of", int(rows.n_tot.sum()), round(cov / rows.n_tot.sum(), 3))
print("C", len(c), "hollows", int(c.n_assigned.sum()), "| Class 1-2 hosts in B", out["design"]["class12_in_B"])
print("D", len(sel), "D range", round(sel.D.min()), round(sel.D.max()), "pi range", sel.pi.min().round(4), sel.pi.max().round(4), "zero-x", int((sel.A_p015 <= 0).sum()))
print("ties at the B cut:", int((rows[(rows.class_source != 'traced') & ~rows.uid.isin(cant)].n_tot == b.n_tot.min()).sum()), "hosts have", int(b.n_tot.min()), "hollows;", int((b.n_tot == b.n_tot.min()).sum()), "taken")
