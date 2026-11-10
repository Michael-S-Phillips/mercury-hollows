# g4 copy of code/rerun/out_v3/rerun/g1_graph_joint/summary_model/add_xrs.py (6 Nov 2026). Changed: none (verbatim) (see make_copies_g1.py).
#!/usr/bin/env python3
"""Summary model, step 7 input (3 Nov 2026): add MESSENGER XRS element ratios (Nittler et al. RDR maps, 2015-04-24, in
GRL_insolation_hollows/data/xrs_maps) to crater_features.csv. Maps are 4 pixels per degree, simple cylindrical,
-180..180 E, 90N at the top, 8-bit, 0 = missing. The PDS labels list LINES = 1440 and LINE_SAMPLES = 720, but the files
are 720 lines x 1440 samples (1,036,800 bytes; read this way, the high-Mg region and the Caloris high-Al/Si interior fall
in the right places). Values are sampled at the crater center: S/Si, Mg/Si, Ca/Si, Al/Si, Fe/Si (NaN where missing), and
the footprint size of the S/Si samples (km). Writes crater_features_xrs.csv."""
import os, numpy as np, pandas as pd
# g1 rerun (6 Nov 2026): copy of analysis_workshop/summary_model/add_xrs.py; only input/output paths changed. SM_DIR (env) is the
# folder holding this catalog's crater_features*.csv and out/ (g1_graph_joint/summary_model/<CATALOG>).
HERE = os.environ["SM_DIR"]
X = os.path.expanduser("~/mnt/hollows-emissivity/manuscripts/GRL_insolation_hollows/data/xrs_maps")
SCALE = {"s__si": 0.147, "mg_si": None, "ca_si": None, "al_si": None, "fe_si": None}
def scale_from_label(fn):
    t = open(fn).read(); import re
    return float(re.search(r"SCALING_FACTOR\s*=\s*([0-9.Ee+-]+)", t).group(1))
d = pd.read_csv(f"{HERE}/crater_features.csv", dtype={"id": str})
row = np.clip(((90 - d.lat.values) * 4).astype(int), 0, 719); col = np.clip(((((d.lon.values + 180) % 360)) * 4).astype(int), 0, 1439)
for el, name in [("s__si", "s_si"), ("mg_si", "mg_si"), ("ca_si", "ca_si"), ("al_si", "al_si"), ("fe_si", "fe_si")]:
    m = np.fromfile(f"{X}/xrs_map_{el}_20150424.img", dtype=np.uint8).reshape(720, 1440)
    sf = scale_from_label(f"{X}/xrs_map_{el}_20150424.lbl"); v = m[row, col].astype(float)
    d[name] = np.where(v > 0, v * sf, np.nan)
r = np.fromfile(f"{X}/xrs_res_s__si_20150424.img", dtype=np.uint8).reshape(720, 1440)
d["xrs_res_km"] = np.where(r[row, col] > 0, r[row, col] * scale_from_label(f"{X}/xrs_res_s__si_20150424.lbl"), np.nan)
d.to_csv(f"{HERE}/crater_features_xrs.csv", index=False)
ok = d[["s_si", "mg_si", "ca_si", "al_si"]].notna().all(axis=1)
print("craters with S, Mg, Ca, Al:", int(ok.sum()), "of", len(d), "; h09 hosts:", int(d[ok].h09.sum()), "of", int(d.h09.sum()))
print(d.loc[ok, ["s_si", "mg_si", "ca_si", "al_si", "fe_si", "xrs_res_km"]].describe().round(3).to_string())
print("mean ratios, hosts vs non-hosts:"); print(d[ok].groupby("h09")[["s_si", "mg_si", "ca_si", "al_si"]].mean().round(4))
