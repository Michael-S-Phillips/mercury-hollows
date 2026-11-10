"""Validation scatter: automatic (v2) vs Blewett et al. (2016) depth per matched detection box (emission <= 30 deg, clean
measurements), coloured by |latitude| >= 60 deg. Writes qc/validation_scatter.png."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
v = pd.read_csv("validation.csv"); v = v[v.emi_le30 & v.clean]; R = json.load(open("validation_record.json"))["emi_le30_clean"]["vs_max"]
fig, ax = plt.subplots(figsize=(5.2, 5)); pole = v.lat.abs() >= 60
ax.scatter(v.blewett_max_m[~pole], v.depth_m[~pole], s=12, c="#3b6ea5", alpha=0.6, lw=0, label="|lat| < 60°")
ax.scatter(v.blewett_max_m[pole], v.depth_m[pole], s=16, c="#c8553d", alpha=0.9, lw=0, label="|lat| ≥ 60°")
lim = [1.5, 300]; ax.plot(lim, lim, "k-", lw=0.8, label="1:1"); ax.plot(lim, np.array(lim) * R["median_ratio"], "k--", lw=0.8, label=f"median ratio {R['median_ratio']:.2f}")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("Blewett et al. (2016) depth, deepest in box (m)"); ax.set_ylabel("Automatic shadow depth (m)")
ax.set_title(f"n = {R['n']} boxes, {R['n_images']} images; r(ln) = {R['pearson_r_log']:.2f}", fontsize=10)
ax.legend(fontsize=8, loc="upper left", frameon=False); fig.tight_layout(); fig.savefig("qc/validation_scatter.png", dpi=150)
