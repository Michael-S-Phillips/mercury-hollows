#!/usr/bin/env python3
# g6_remaining (6 Nov 2026): new script. Copies the knob color tests and changes only input/output paths (each asserted once).
"""Copies analysis_workshop/melt_source_tests/knobs/{common, t1_knob_color, t1_peak_color, t1_montes_color, summarize}.py into
knobs/. Env G6_CAT = v2 | v3 | v3b selects the hollow catalog (v2: code/rerun/out_v2/hollows_canonical.csv; v3/v3b:
code/rerun/out_v3/hollows_canonical_<cat>.csv; non-trace hollows in both, as in the original) and the PCA loadings for the PC indices
(v2: color_mdr/pca_record.json; v3/v3b: this folder's color/out_<cat>/pca_record.json). Outputs go to knobs/out_<cat>/.
summarize.py reads the Test 2 (hollow density) inputs from g3's rerun (g3_catalog/melt_source/knobs/out_<cat>) for v3/v3b and from
the original folder for v2, so its record holds both tests for each catalog."""
import os
HE = os.path.expanduser("~/mnt/hollows-emissivity"); SRC = f"{HE}/manuscripts/PlanetaryResearch_hollows/analysis_workshop/melt_source_tests/knobs"
HERE = os.path.dirname(os.path.abspath(__file__)); DST = f"{HERE}/knobs"; LOG = []
def copy(name, reps, tag):
    t = open(f"{SRC}/{name}").read()
    for a, b in reps:
        assert t.count(a) == 1, (name, a, t.count(a)); t = t.replace(a, b)
    head = ""
    if t.startswith("#!"): head, t = t.split("\n", 1); head += "\n"
    open(f"{DST}/{name}", "w").write(head + f"# g6_remaining copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/melt_source_tests/knobs/{name} (6 Nov 2026). Changed: {tag}\n" + t)
    LOG.append(f"knobs/{name} <- analysis_workshop/melt_source_tests/knobs/{name} ({len(reps)} replacements)")
copy("common.py", [
    ('HERE = f"{AW}/melt_source_tests/knobs"',
     'G6_CAT = os.environ.get("G6_CAT", "v3"); G6 = f"{ROOT}/code/rerun/out_v3/rerun/g6_remaining"\n'
     'HERE = f"{G6}/knobs/out_{G6_CAT}"; os.makedirs(HERE, exist_ok=True)\n'
     'ORIG = f"{AW}/melt_source_tests/knobs"\n'
     'T2 = ORIG if G6_CAT == "v2" else f"{ROOT}/code/rerun/out_v3/rerun/g3_catalog/melt_source/knobs/out_{G6_CAT}"\n'
     'CANON = f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv" if G6_CAT == "v2" else f"{ROOT}/code/rerun/out_v3/hollows_canonical_{G6_CAT}.csv"\n'
     'PCAREC = f"{AW}/color_mdr/pca_record.json" if G6_CAT == "v2" else f"{G6}/color/out_{G6_CAT}/pca_record.json"'),
    ('h = pd.read_csv(f"{ROOT}/code/rerun/out_v2/hollows_canonical.csv", usecols=', 'h = pd.read_csv(CANON, usecols='),
    ('_P = _json.load(open(f"{AW}/color_mdr/pca_record.json"))', '_P = _json.load(open(PCAREC))')],
    "output folder (knobs/out_<G6_CAT>), hollow catalog, PCA-record and Test 2 input paths only.")
for n in ["t1_knob_color.py", "t1_peak_color.py", "t1_montes_color.py"]:
    copy(n, [], "nothing (paths come from common.py).")
copy("summarize.py", [('from common import HERE, AW', 'from common import HERE, AW, T2'),
                      ('pd.read_csv(f"{HERE}/knob_hollows.csv")', 'pd.read_csv(f"{T2}/knob_hollows.csv")'),
                      ('pd.read_csv(f"{HERE}/_chunks/peak_hollows_00000_04824.csv")', 'pd.read_csv(f"{T2}/_chunks/peak_hollows_00000_04824.csv")'),
                      ('json.load(open(f"{HERE}/region_hollows.json"))', 'json.load(open(f"{T2}/region_hollows.json"))')],
     "Test 2 input paths (T2) only; Test 1 inputs and all outputs in knobs/out_<G6_CAT>.")
open(f"{HERE}/make_copies_knobs.log", "w").write("\n".join(LOG) + "\n"); print("\n".join(LOG))
