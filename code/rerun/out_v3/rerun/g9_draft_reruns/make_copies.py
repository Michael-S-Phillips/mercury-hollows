#!/usr/bin/env python3
"""g9_draft_reruns (6 Nov 2026): copies the scripts behind four draft analyses into this folder and changes only inputs,
outputs, class rules and the requested new measures. Every change is an exact string replacement asserted to occur once,
so the changes can be read here. Sources are never modified.
Env settings of the copies: G9_CAT = v2 | v3 (catalog), G9_RULES = fix | old (terrain classes: g7 terrain_classes.classify or
the old peak_wall.py lines), G9_GRID = v3 | v2 (completeness grid: g4 completeness_10ppd_v3.npz or out_v2/completeness_10ppd.npz),
G9_ST_CAT = v2 | v3 | v3r (size-transition crater table)."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); RR = f"{HE}/code/rerun/out_v3/rerun"; G9 = f"{RR}/g9_draft_reruns"
WS = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/melt_source_tests"
G3P = f"{RR}/g3_catalog/melt_source/positions"; G3K = f"{RR}/g3_catalog/melt_source/knobs"
def copy(src, dst, first, reps=()):
    t = open(src).read()
    for a, b in reps:
        assert t.count(a) == 1, (src, a[:80], t.count(a)); t = t.replace(a, b)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    open(dst, "w").write(first.rstrip("\n") + "\n" + t); print("wrote", dst)
CFG = ('G9 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g9_draft_reruns")\n'
       'CAT = _os.environ.get("G9_CAT", "v3"); RULES = _os.environ.get("G9_RULES", "fix"); GRID = _os.environ.get("G9_GRID", "v3")\n')
G3HDR = ('G3 = _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g3_catalog")\n'
         'CAT = _os.environ.get("G3_CAT", "v3"); ODIR = f"{G3}/melt_source/positions/out_{CAT}"; _os.makedirs(ODIR + "/cache", exist_ok=True)\n')
GRIDLINE_OLD = 'C = np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz"); Cg = C["C"][list(C["Dref"]).index(1.0)]'
GRIDLINE_NEW = ('C = np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz" if GRID == "v2" else '
                'f"{HE}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz"); Cg = C["C"][list(C["Dref"]).index(1.0)]')
CLS_IMPORT = ('sys.path.insert(0, _os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v3/rerun/g7_terrain_fix")); '
              'from terrain_classes import classify\n')

# ---------------- 1. positions (from the g3 v3 copies) ----------------
copy(f"{G3P}/positions_build.py", f"{G9}/positions/positions_build.py",
     "# g9 copy of g3_catalog/melt_source/positions/positions_build.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID, output folder "
     "positions/out_<CAT>_<RULES>_grid<GRID>/cache, completeness grid (v3 = g4 completeness_10ppd_v3.npz), and the terrain-class lines "
     "(RULES=fix: g7 terrain_classes.classify, az = degrees(arctan2(dy, dx)) % 360). d_nonfloor keeps its definition (nearest non-floor "
     "cell within 1.05 R); under the new classes the floor reaches 1.0 R, so its radial cut moves from 0.8 R to 1.0 R.",
     [(G3HDR, CFG + 'ODIR = f"{G9}/positions/out_{CAT}_{RULES}_grid{GRID}"; _os.makedirs(ODIR + "/cache", exist_ok=True)\n' + CLS_IMPORT),
      (GRIDLINE_OLD, GRIDLINE_NEW),
      ("    floor = np.percentile(fz, 20); h = z-floor; H = HT\n"
       "    peak = (rho <= (0.6 if ring else RPK)) & (h > H)\n"
       "    wall = (rho >= (0.7 if ring else 0.55)) & (rho <= 1.05) & (h > H)\n"
       "    flo = (rho <= 0.8) & (h <= H)\n",
       "    floor = np.percentile(fz, 20); h = z-floor; H = HT\n"
       "    if RULES == \"fix\":\n"
       "        az = np.degrees(np.arctan2(dy, dx)) % 360\n"
       "        KC = classify(z, rho, az, ring, HT, RPK); assert KC[\"F0\"] == floor\n"
       "        peak, wall, flo = KC[\"peak\"], KC[\"wall\"], KC[\"flo\"]\n"
       "    else:\n"
       "        peak = (rho <= (0.6 if ring else RPK)) & (h > H)\n"
       "        wall = (rho >= (0.7 if ring else 0.55)) & (rho <= 1.05) & (h > H)\n"
       "        flo = (rho <= 0.8) & (h <= H)\n")])
ADD = '''# ---- g9 additions (6 Nov 2026): radius-scaled floor margins (distance / host radius) and the central-floor split ----
G9X = {"note": "floor hollows/fields vs completeness-weighted floor area; d_marg_R = min(distance to wall-class, to peak-class cell) / host radius; "
               "radius_stratified = crater-matched null within crater x radial bins (RB)", "RB": RB}
for lab, pts_f in (("hollows", hf_), ("fields", ff_)):
    o = {}
    for thr in (0.05, 0.1):
        A, n = mats(fl, pts_f, "d_marg_R", [0, thr, INF]); o[f"floor_within{thr}R_of_wall_or_peak"] = split_stats(A, n, np.array([1, 0], bool), ALL, W_ALL)
        o[f"floor_within{thr}R_of_wall_or_peak_radius_stratified"] = strat_split(fl, pts_f, "d_marg_R", thr, ALL, W_ALL)
    o["floor_within3km_of_wall_or_peak_radius_stratified"] = strat_split(fl, pts_f, "d_marg", 3, ALL, W_ALL)
    A, n = mats(fl.assign(r=fl.rho), pts_f.assign(r=pts_f.rho), "r", [0, 0.2, INF]); o["floor_within0.2R_of_center_vs_rest"] = split_stats(A, n, np.array([1, 0], bool), ALL, W_ALL)
    G9X[lab] = o
REC["g9_additions"] = G9X
for lab in ("hollows", "fields"):
    for k, s in G9X[lab].items():
        print("G9", lab, k, {kk: (round(v[0], 3), round(v[1], 3), round(v[2], 3)) if isinstance(v, list) and len(v) == 3 and v[0] is not None and v[1] is not None else v for kk, v in s.items() if kk in ("n_points", "n_in", "n_craters_with_points", "SDR", "DR")})
'''
copy(f"{G3P}/positions_analyze.py", f"{G9}/positions/positions_analyze.py",
     "# g9 copy of g3_catalog/melt_source/positions/positions_analyze.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID and paths "
     "(reads <out>/cache or G9_CACHE_DIR), radial bins RB extended to 1.0 R under RULES=fix (the new floor class reaches 1.0 R), and new "
     "measures in REC['g9_additions']: floor margins within 0.05 R and 0.1 R of wall or peak cells (distance / host radius), their "
     "radius-stratified versions, and floor within 0.2 R of the center vs the rest of the floor as a split ratio.",
     [(G3HDR, CFG + 'ODIR = f"{G9}/positions/out_{CAT}_{RULES}_grid{GRID}"; _os.makedirs(ODIR + "/cache", exist_ok=True)\n'
              'CACHE_DIR = _os.environ.get("G9_CACHE_DIR", ODIR + "/cache")\n'),
      ('d = np.load(f"{HERE}/cache/H{H}_000_354.npz", allow_pickle=True)', 'd = np.load(f"{CACHE_DIR}/H{H}_000_354.npz", allow_pickle=True)'),
      ('    fl["d_marg"] = np.minimum(fl.d_wall, fl.d_peak); ho["d_marg"] = np.minimum(ho.d_wall, ho.d_peak); fd["d_marg"] = np.minimum(fd.d_wall, fd.d_peak)\n',
       '    fl["d_marg"] = np.minimum(fl.d_wall, fl.d_peak); ho["d_marg"] = np.minimum(ho.d_wall, ho.d_peak); fd["d_marg"] = np.minimum(fd.d_wall, fd.d_peak)\n'
       '    Rk = dict(zip(cr.uid.astype(int), cr.D/2))\n'
       '    for t in (fl, ho, fd): t["d_marg_R"] = t.d_marg/t.uid.astype(int).map(Rk)\n'),
      ('RB = [0, 0.2, 0.4, 0.6, 0.81]; NRB = len(RB)-1', 'RB = [0, 0.2, 0.4, 0.6, 0.81] if RULES == "old" else [0, 0.2, 0.4, 0.6, 0.8, 1.01]; NRB = len(RB)-1'),
      ('json.dump(REC, open(f"{HERE}/positions_record.json", "w")', ADD + 'json.dump(REC, open(f"{HERE}/positions_record.json", "w")')])

# ---------------- 2. slope follow-up (from the workshop originals; no g3 copy exists) ----------------
SLHDR = ('import os as _os, sys\n' + CFG + 'ODIR = f"{G9}/positions/out_{CAT}_{RULES}_grid{GRID}"; _os.makedirs(ODIR + "/cache/_chunks", exist_ok=True)\n'
         'HOSTS = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/latitude/aspect/host_assignments.csv") if CAT == "v2"\n'
         '         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/host_assignments_{CAT}.csv"))\n'
         'CANON = (_os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/hollows_canonical.csv") if CAT == "v2"\n'
         '         else _os.path.expanduser(f"~/mnt/hollows-emissivity/code/rerun/out_v3/hollows_canonical_{CAT}.csv"))\n')
copy(f"{WS}/positions/positions_slope_build.py", f"{G9}/positions/positions_slope_build.py",
     "# g9 copy of analysis_workshop/melt_source_tests/positions/positions_slope_build.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID, "
     "host and catalog paths, completeness grid, output folder positions/out_<CAT>_<RULES>_grid<GRID>/cache, the peak-class line (RULES=fix: "
     "g7 terrain_classes.classify), and optional START STOP arguments (chunks in cache/_chunks, joined by merge_chunks.py) for the 180 s limit.",
     [('import os, numpy as np, pandas as pd, rasterio\n', 'import os, numpy as np, pandas as pd, rasterio\n' + SLHDR + CLS_IMPORT),
      ('HE = os.path.expanduser("~/mnt/hollows-emissivity"); OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")',
       'HE = os.path.expanduser("~/mnt/hollows-emissivity"); OUT = ODIR + "/cache"'),
      ('ha = pd.read_csv(f"{HE}/code/rerun/latitude/aspect/host_assignments.csv")', 'ha = pd.read_csv(HOSTS)'),
      ('hc = pd.read_csv(f"{HE}/code/rerun/out_v2/hollows_canonical.csv")[["hollow_id", "field_id", "is_crater_trace"]]',
       'hc = pd.read_csv(CANON)[["hollow_id", "field_id", "is_crater_trace"]]'),
      (GRIDLINE_OLD, GRIDLINE_NEW),
      ('  for _, c in hosts.iterrows():', '  A0 = int(sys.argv[1]) if len(sys.argv) > 2 else 0; A1 = int(sys.argv[2]) if len(sys.argv) > 2 else len(hosts)\n  for _, c in hosts.iloc[A0:A1].iterrows():'),
      ('    peak = (rho <= (0.6 if ring else RPK)) & (h > H)\n',
       '    if RULES == "fix":\n'
       '        az = np.degrees(np.arctan2(dy, dx)) % 360\n'
       '        KC = classify(z, rho, az, ring, HT, RPK); assert KC["F0"] == floor; peak = KC["peak"]\n'
       '    else:\n'
       '        peak = (rho <= (0.6 if ring else RPK)) & (h > H)\n'),
      ('np.savez_compressed(f"{OUT}/slope_H250.npz", craters=craters, peak=np.vstack(PK).astype(np.float32), hollows=np.vstack(HO), fields=np.array(FD, float))',
       'np.savez_compressed(f"{OUT}/slope_H250.npz" if (A0, A1) == (0, len(hosts)) else f"{OUT}/_chunks/slope_H250_{A0:03d}_{A1:03d}.npz", craters=craters, '
       'peak=np.vstack(PK).astype(np.float32) if PK else np.zeros((0, 6), np.float32), hollows=np.vstack(HO) if HO else np.zeros((0, 6)), '
       'fields=np.array(FD, float) if FD else np.zeros((0, 7)))')])
copy(f"{WS}/positions/positions_slope.py", f"{G9}/positions/positions_slope.py",
     "# g9 copy of analysis_workshop/melt_source_tests/positions/positions_slope.py (6 Nov 2026). Changed: env G9_CAT/G9_RULES/G9_GRID and the "
     "folder it reads and writes (positions/out_<CAT>_<RULES>_grid<GRID>).",
     [('HERE = os.path.dirname(os.path.abspath(__file__)); NB = 2000; INF = np.inf',
       'import os as _os\n' + CFG + 'HERE = f"{G9}/positions/out_{CAT}_{RULES}_grid{GRID}"; NB = 2000; INF = np.inf')])

# ---------------- 3. knobs (from the g3 v3 copies; no terrain classes are used) ----------------
copy(f"{G3K}/common.py", f"{G9}/knobs/common.py",
     "# g9 copy of g3_catalog/melt_source/knobs/common.py (6 Nov 2026). Changed: env G9_CAT/G9_GRID, output folder knobs/out_<CAT>_grid<GRID>, "
     "and the completeness grid in load_C (v3 = g4 completeness_10ppd_v3.npz).",
     [('HERE = f"{ROOT}/code/rerun/out_v3/rerun/g3_catalog/melt_source/knobs/out_" + os.environ.get("G3_CAT", "v3"); os.makedirs(HERE, exist_ok=True)',
       'G9CAT = os.environ.get("G9_CAT", "v3"); G9GRID = os.environ.get("G9_GRID", "v3")\n'
       'HERE = f"{ROOT}/code/rerun/out_v3/rerun/g9_draft_reruns/knobs/out_{G9CAT}_grid{G9GRID}"; os.makedirs(HERE, exist_ok=True)'),
      ('CANON = (f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv" if os.environ.get("G3_CAT", "v3") == "v2"\n'
       '         else f"{ROOT}/code/rerun/out_v3/hollows_canonical_" + os.environ.get("G3_CAT", "v3") + ".csv")',
       'CANON = (f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv" if G9CAT == "v2"\n'
       '         else f"{ROOT}/code/rerun/out_v3/hollows_canonical_" + G9CAT + ".csv")'),
      ('    C = np.load(f"{ROOT}/code/rerun/out_v2/completeness_10ppd.npz")',
       '    C = np.load(f"{ROOT}/code/rerun/out_v2/completeness_10ppd.npz" if G9GRID == "v2" else f"{ROOT}/code/rerun/out_v3/rerun/g4_completeness/completeness_10ppd_v3.npz")')])
for nm in ("t2_knob_hollows.py", "t2_peak_hollows.py", "summarize.py"):
    copy(f"{G3K}/{nm}", f"{G9}/knobs/{nm}", f"# g9 copy of g3_catalog/melt_source/knobs/{nm} (6 Nov 2026). Unchanged; paths, catalog and grid come from common.py.")

# ---------------- 4. size transition (from the workshop originals) ----------------
copy(f"{WS}/size_transition/st_common.py", f"{G9}/size_transition/st_common.py",
     "# g9 copy of analysis_workshop/melt_source_tests/size_transition/st_common.py (6 Nov 2026). Changed: env G9_ST_CAT = v2 | v3 | v3r selects "
     "the crater table (v2: summary_model/crater_features.csv; v3: g1_graph_joint/summary_model/v3; v3r: g4_completeness/g1/summary_model/v3r, "
     "v3 hosts with C1 from the v3 refit grid) and the hollow catalog (v3, v3r: out_v3/hollows_canonical_v3.csv); output folder size_transition/out_<CAT>.",
     [('HERE = os.path.dirname(os.path.abspath(__file__))',
       'STCAT = os.environ.get("G9_ST_CAT", "v3r")\n'
       'HERE = f"{ROOT}/code/rerun/out_v3/rerun/g9_draft_reruns/size_transition/out_{STCAT}"; os.makedirs(HERE, exist_ok=True)'),
      ('FEAT = f"{PR}/analysis_workshop/summary_model/crater_features.csv"',
       'FEAT = {"v2": f"{PR}/analysis_workshop/summary_model/crater_features.csv",\n'
       '        "v3": f"{ROOT}/code/rerun/out_v3/rerun/g1_graph_joint/summary_model/v3/crater_features.csv",\n'
       '        "v3r": f"{ROOT}/code/rerun/out_v3/rerun/g4_completeness/g1/summary_model/v3r/crater_features.csv"}[STCAT]'),
      ('CANON = f"{ROOT}/hollow-networks/data_v2/hollows_canonical.csv"',
       'CANON = f"{ROOT}/hollow-networks/data_v2/hollows_canonical.csv" if STCAT == "v2" else f"{ROOT}/code/rerun/out_v3/hollows_canonical_v3.csv"')])
for nm in ("st1_bins.py", "st2_breaks.py", "st2b_null.py", "st2c_blurred.py", "st3_nested.py", "st4_morph.py", "st5_caloris.py", "st6_record.py"):
    copy(f"{WS}/size_transition/{nm}", f"{G9}/size_transition/{nm}",
         f"# g9 copy of analysis_workshop/melt_source_tests/size_transition/{nm} (6 Nov 2026). Unchanged; inputs and outputs come from st_common.py.")
print("done")
