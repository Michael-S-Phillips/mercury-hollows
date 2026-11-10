#!/usr/bin/env python3
"""g2_insolation (6 Nov 2026): copy each original analysis script into this folder and change only its input paths,
output paths and environment settings. Every replacement is listed below and must match the original text exactly
(the count is asserted). A block that defines IN (inputs/<G2_CAT>) and G2_OUT (out/<G2_CAT>) is inserted before the
first import; G2_CAT (v2, v3 or v3b) comes from the environment. The originals are not modified."""
import os, re
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
PR = "manuscripts/PlanetaryResearch_hollows"; W6 = f"{PR}/analysis_workshop/insolation_v6"
BLOCK = ('# --- g2 rerun (6 Nov 2026): catalog inputs and output folder; paths and env only ---\n'
         'import os as _g2_os\n'
         'G2_HERE = _g2_os.path.dirname(_g2_os.path.abspath(__file__)); G2_ROOT = _g2_os.path.abspath(_g2_os.path.join(G2_HERE, *[".."] * 5))\n'
         'G2_CAT = _g2_os.environ.get("G2_CAT", "v3"); IN = f"{G2_HERE}/inputs/{G2_CAT}"; G2_OUT = f"{G2_HERE}/out/{G2_CAT}"\n'
         '_g2_os.makedirs(G2_OUT, exist_ok=True)\n'
         'G2_GRLDATA = f"{G2_ROOT}/manuscripts/GRL_insolation_hollows/data"\n'
         '# --- end g2 block ---\n')
SPEC = {
 f"{PR}/analysis_v3/dT_hosting.py": [
   ('s = pd.read_csv(f"{G}/data/test2_craters_h18.csv")', 's = pd.read_csv(f"{IN}/test2_craters_h18.csv")', 1),
   ('open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"dT_record.json"), "w")', 'open(f"{G2_OUT}/dT_record.json", "w")', 1)],
 f"{PR}/analysis_v3/dT_plains.py": [
   ('s = pd.read_csv(f"{G}/data/test2_craters_h18.csv")', 's = pd.read_csv(f"{IN}/test2_craters_h18.csv")', 1),
   ('np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)),"nsp_mask_10ppd_r40.npy"))', 'np.load(f"{G2_ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_v3/nsp_mask_10ppd_r40.npy")', 1),
   ('open("dT_plains_record.json","w")', 'open(f"{G2_OUT}/dT_plains_record.json","w")', 1)],
 f"{PR}/analysis_v7/dT_hosting_completeness.py": [
   ('s = pd.read_csv(f"{G}/data/test2_craters_h18.csv")', 's = pd.read_csv(f"{IN}/test2_craters_h18.csv")', 1),
   ('pd.read_csv(f"{MCG}/outputs/crater_dag_nodes.csv.gz", dtype={"id": str}, low_memory=False)',
    'pd.read_csv(f"{G2_ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_v9/crater_dag_nodes.csv", dtype={"id": str}, low_memory=False)', 1),
   ('open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dT_hosting_completeness_record.json"), "w")', 'open(f"{G2_OUT}/dT_hosting_completeness_record.json", "w")', 1)],
 f"{PR}/analysis_v3/transition_layer_hilat.py": [
   ('z=np.load(f"{TH}/thermal_grid.npz")', 'z=np.load(f"{G2_GRLDATA}/thermal_grid.npz")', 1),
   ('open("transition_layer_hilat_record.json","w")', 'open(f"{G2_OUT}/transition_layer_hilat_record.json","w")', 2)],
 "manuscripts/GRL_insolation_hollows/scripts/latcorr.py": [
   ('U="/mnt/user-data/uploads/hollows-emissivity"', 'U=G2_ROOT', 1),
   ('nsp=np.load("nsp_mask_10ppd_r40.npy")', 'nsp=np.load(f"{G2_GRLDATA}/nsp_mask_10ppd_r40.npy")', 1),
   ('f=pd.read_csv(f"{U}/code/rerun/out_v2/fields.csv"); h=pd.read_csv(f"{U}/code/rerun/out_v2/hollows_canonical.csv")',
    'f=pd.read_csv(f"{IN}/fields.csv"); h=pd.read_csv(f"{IN}/hollows_canonical.csv")', 1),
   ('T.to_csv("density_by_band.csv",index=False)', 'T.to_csv(f"{G2_OUT}/density_by_band.csv",index=False)', 1),
   ('f.to_csv("fields_flagged.csv",index=False)', 'f.to_csv(f"{G2_OUT}/fields_flagged.csv",index=False)', 1)],
 "manuscripts/GRL_insolation_hollows/scripts/latcorr3.py": [
   ('U="/mnt/user-data/uploads/hollows-emissivity"', 'U=G2_ROOT', 1),
   ('T=pd.read_csv("density_by_band.csv"); f=pd.read_csv("fields_flagged.csv")', 'T=pd.read_csv(f"{G2_OUT}/density_by_band.csv"); f=pd.read_csv(f"{G2_OUT}/fields_flagged.csv")', 1),
   ('nsp=np.load("nsp_mask_10ppd_r40.npy")', 'nsp=np.load(f"{G2_GRLDATA}/nsp_mask_10ppd_r40.npy")', 1),
   ('open("latcorr_density_record.json","w")', 'open(f"{G2_OUT}/latcorr_density_record.json","w")', 1)],
 "manuscripts/GRL_insolation_hollows/scripts/latcorr5.py": [
   ('U="/mnt/user-data/uploads/hollows-emissivity/code/rerun"', 'U=f"{G2_ROOT}/code/rerun"', 1),
   ('nsp=np.load("nsp_mask_10ppd_r40.npy")', 'nsp=np.load(f"{G2_GRLDATA}/nsp_mask_10ppd_r40.npy")', 1),
   ('h=pd.read_csv(f"{U}/out_v2/hollows_canonical.csv"); f=pd.read_csv(f"{U}/out_v2/fields.csv")',
    'h=pd.read_csv(f"{IN}/hollows_canonical_inc.csv"); f=pd.read_csv(f"{IN}/fields.csv")', 1),
   ('open("latcorr_inc60_75_record.json","w")', 'open(f"{G2_OUT}/latcorr_inc60_75_record.json","w")', 1)],
 "code/rerun/latitude/hotwarm/hotwarm.py": [
   ('DD = "/Users/phillipsm/Documents/Research/hollows-emissivity/hollow-networks/data_v2/"', 'DD = f"{G2_ROOT}/code/rerun/out_v2/"   # same files (md5) as hollow-networks/data_v2', 1),
   ('OUT = "/Users/phillipsm/Documents/Research/hollows-emissivity/code/rerun/latitude/hotwarm/"', 'OUT = G2_OUT + "/"', 1),
   ('fields = pd.read_csv(DD + "fields.csv")', 'fields = pd.read_csv(IN + "/fields.csv")', 1),
   ('hollows = pd.read_csv(DD + "hollows_canonical.csv")', 'hollows = pd.read_csv(IN + "/hollows_canonical.csv")', 1)],
 "code/rerun/latitude/aspect/aspect_stats.py": [
   ('df = pd.read_csv("host_assignments.csv")', 'df = pd.read_csv(f"{IN}/host_assignments.csv")', 1),
   ('wall.to_csv("aspect_wall_hollows.csv", index=False)', 'wall.to_csv(f"{G2_OUT}/aspect_wall_hollows.csv", index=False)', 1),
   ('with open("aspect_record.json", "w") as f:', 'with open(f"{G2_OUT}/aspect_record.json", "w") as f:', 1)],
 f"{W6}/common_null.py": [
   ('ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))', 'ROOT = G2_ROOT', 1),
   ('fields = pd.read_csv(DD + "fields.csv"); hollows = pd.read_csv(DD + "hollows_canonical.csv")',
    'fields = pd.read_csv(IN + "/fields.csv"); hollows = pd.read_csv(IN + "/hollows_canonical.csv")', 1)],
 f"{W6}/cos2lon_lat15.py": [('open(f"{HERE}/cos2lon_fit_lat15.json", "w")', 'open(f"{G2_OUT}/cos2lon_fit_lat15.json", "w")', 1)],
 f"{W6}/size_quirk.py": [
   ('pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv"', 'pd.read_csv(f"{IN}/crater_features.csv"', 1),
   ('open(f"{HERE}/size_quirk.json", "w")', 'open(f"{G2_OUT}/size_quirk.json", "w")', 1)],
 f"{W6}/band_crater_hotwarm.py": [
   ('pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv"', 'pd.read_csv(f"{IN}/crater_features.csv"', 1),
   ('open(f"{HERE}/band_crater_hotwarm.json", "w")', 'open(f"{G2_OUT}/band_crater_hotwarm.json", "w")', 1)],
 f"{W6}/hollow_rich_hosts.py": [
   ('pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv"', 'pd.read_csv(f"{IN}/crater_features.csv"', 1),
   ('open(f"{HERE}/hollow_rich_hosts.json", "w")', 'open(f"{G2_OUT}/hollow_rich_hosts.json", "w")', 1)],
 f"{W6}/fractions.py": [
   ('pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv"', 'pd.read_csv(f"{IN}/crater_features.csv"', 1),
   ('open(f"{HERE}/fractions.json", "w")', 'open(f"{G2_OUT}/fractions.json", "w")', 1)],
 f"{W6}/threshold.py": [
   ('pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv"', 'pd.read_csv(f"{IN}/crater_features.csv"', 1),
   ('open(f"{HERE}/threshold.json", "w")', 'open(f"{G2_OUT}/threshold.json", "w")', 1),
   ('fig.savefig(f"{HERE}/threshold_shape.png", dpi=200)', 'fig.savefig(f"{G2_OUT}/threshold_shape.png", dpi=200)', 1)],
 f"{W6}/tm_vs_dT_check.py": [
   ('pd.read_csv(f"{PR}/analysis_workshop/summary_model/crater_features.csv"', 'pd.read_csv(f"{IN}/crater_features.csv"', 1),
   ('open(f"{HERE}/tm_vs_dT_check.json", "w")', 'open(f"{G2_OUT}/tm_vs_dT_check.json", "w")', 1)],
 f"{W6}/stage_test.py": [
   ('H = pd.read_csv(DD + "hollows_canonical.csv")', 'H = pd.read_csv(IN + "/hollows_canonical.csv")', 1),
   ('pd.read_csv(f"{ROOT}/code/rerun/latitude/incidence_null/canonical_stage_labels.csv")', 'pd.read_csv(f"{IN}/canonical_stage_labels.csv")', 1),
   ('A = pd.read_csv(f"{G}/aspect_wall_hollows.csv")', 'A = pd.read_csv(f"{G2_OUT}/aspect_wall_hollows.csv")', 1),
   ('open(f"{HERE}/stage_test.json", "w")', 'open(f"{G2_OUT}/stage_test.json", "w")', 1)],
 f"{W6}/stage_test_lowinc.py": [('open(f"{HERE}/stage_test_lowinc.json", "w")', 'open(f"{G2_OUT}/stage_test_lowinc.json", "w")', 1)],
 f"{PR}/figscripts/make_figs_pr_insolation_v9.py": [
   ('_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))', '_ROOT = G2_ROOT', 1),
   ('OUT = os.environ.get("GRL_FIGS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))', 'OUT = os.environ.get("GRL_FIGS", G2_OUT)', 1),
   ('json.load(open(f"{_ROOT}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/insolation_v6/cos2lon_fit_lat15.json"))', 'json.load(open(f"{G2_OUT}/cos2lon_fit_lat15.json"))', 1),
   ('s = pd.read_csv(f"{G}/data/test2_craters_h18.csv")', 's = pd.read_csv(f"{IN}/test2_craters_h18.csv")', 1),
   ('a = pd.read_csv(f"{G}/data/aspect_wall_hollows.csv")', 'a = pd.read_csv(f"{G2_OUT}/aspect_wall_hollows.csv")', 1)],
 f"{W6}/stage_fieldlevel.py": [('open(f"{HERE}/stage_fieldlevel.json", "w")', 'open(f"{G2_OUT}/stage_fieldlevel.json", "w")', 1)],
}
NO_BLOCK = {f"{W6}/cos2lon_lat15.py", f"{W6}/size_quirk.py", f"{W6}/band_crater_hotwarm.py", f"{W6}/hollow_rich_hosts.py",
            f"{W6}/fractions.py", f"{W6}/threshold.py", f"{W6}/tm_vs_dT_check.py", f"{W6}/stage_test.py",
            f"{W6}/stage_test_lowinc.py", f"{W6}/stage_fieldlevel.py"}      # these get IN / G2_OUT from common_null import *
RENAME = {"fractions.py": "fractions_g2.py"}   # a file named fractions.py here would shadow the standard-library module
log = []
for rel, reps in SPEC.items():
    src = open(f"{ROOT}/{rel}").read(); txt = src
    for old, new, n in reps:
        assert txt.count(old) == n, (rel, old, txt.count(old)); txt = txt.replace(old, new)
    if rel not in NO_BLOCK:
        lines = txt.split("\n"); i = next(k for k, l in enumerate(lines) if re.match(r"^(import|from) ", l))
        txt = "\n".join(lines[:i]) + ("\n" if i else "") + BLOCK + "\n".join(lines[i:])
    hdr = f"# g2 copy of {rel} (6 Nov 2026). Changed: input/output paths and env only (see make_copies.py).\n"
    if txt.startswith("#!"): a, b = txt.split("\n", 1); txt = a + "\n" + hdr + b
    else: txt = hdr + txt
    name = RENAME.get(os.path.basename(rel), os.path.basename(rel)); open(f"{HERE}/{name}", "w").write(txt); log.append(f"{name} <- {rel} ({len(reps)} replacements)")
open(f"{HERE}/make_copies.log", "w").write("\n".join(log) + "\n"); print("\n".join(log))
