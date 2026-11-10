#!/usr/bin/env python3
"""Copies the knob hollow-density scripts (melt_source_tests/knobs: common.py, t2_knob_hollows.py, t2_peak_hollows.py, summarize.py)
into melt_source/knobs/ and changes only input paths, output paths and env settings. Env G3_CAT = v2 | v3 | v3b; outputs in
melt_source/knobs/out_{G3_CAT}/. The color tests (t1_*) mask hollows only in the matrix annulus and read the MDR maps; they were
not rerun, and summarize.py reads their v2 outputs from the original folder."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); SRC = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/melt_source_tests/knobs"
DST = os.path.dirname(os.path.abspath(__file__)) + "/melt_source/knobs"; os.makedirs(DST, exist_ok=True)
def copy(name, reps):
    t = open(f"{SRC}/{name}").read()
    for a, b in reps:
        assert t.count(a) == 1, (name, a, t.count(a)); t = t.replace(a, b)
    open(f"{DST}/{name}", "w").write(f"# g3_catalog copy of analysis_workshop/melt_source_tests/knobs/{name} (6 Nov 2026): only input/output paths and env settings changed.\n" + t)
copy("common.py", [('HERE = f"{AW}/melt_source_tests/knobs"',
                    'HERE = f"{ROOT}/code/rerun/out_v3/rerun/g3_catalog/melt_source/knobs/out_" + os.environ.get("G3_CAT", "v3"); os.makedirs(HERE, exist_ok=True)\n'
                    'ORIG = f"{AW}/melt_source_tests/knobs"\n'
                    'CANON = (f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv" if os.environ.get("G3_CAT", "v3") == "v2"\n'
                    '         else f"{ROOT}/code/rerun/out_v3/hollows_canonical_" + os.environ.get("G3_CAT", "v3") + ".csv")'),
                   ('h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", usecols=', 'h = pd.read_csv(CANON, usecols=')])
copy("t2_knob_hollows.py", [])
copy("t2_peak_hollows.py", [])
copy("summarize.py", [('from common import HERE, AW', 'from common import HERE, AW, ORIG'),
                      ('pd.read_csv(f"{HERE}/_chunks/knob_color_0000_0335.csv")', 'pd.read_csv(f"{ORIG}/_chunks/knob_color_0000_0335.csv")'),
                      ('pd.read_csv(f"{HERE}/_chunks/peak_color_0000_1163.csv")', 'pd.read_csv(f"{ORIG}/_chunks/peak_color_0000_1163.csv")'),
                      ('pd.read_csv(f"{HERE}/montes_color.csv")', 'pd.read_csv(f"{ORIG}/montes_color.csv")')])
print("copied")
