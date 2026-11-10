#!/usr/bin/env python3
"""g4 (6 Nov 2026): copies of the S2 completeness-group scripts (manuscripts/PlanetaryResearch_hollows/analysis_workshop/
completeness_ref/ref_size_check.py, ref_size_insolation.py; completeness_sfd/crater_completeness.py, hosting_fits.py,
joint_thinned.py). Changes, each asserted to match exactly once: input paths (crater table, completeness grid, detection model
and crater_completeness table chosen by the env S2_CAT), output paths (this folder, file names tagged by S2_CAT), the size-
distribution exponents of E (1.7, the v3 value of g3, and 1.8 added to 1.3, 1.6, 1.9), and the measure lists of hosting_fits.py
and joint_thinned.py (lnE_1.7, lnE_1.8 added). The imports of common_null.hot_mask and step7_improved come from the g4 copies
(g4/g2/common_null.py is not imported; hot_mask is reproduced verbatim in s2cfg.py) so that nothing is written under manuscripts/.
S2_CAT: v2 (original inputs; reproduction), v3 (v3 crater table, old grid, v2 detection model), v3r (v3 crater table, grid and
detection model of the refit), v3br (the same for v3b)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.abspath(os.path.join(HERE, *[".."] * 6, "manuscripts/PlanetaryResearch_hollows/analysis_workshop"))
HDR = "import os as _o, sys as _s; _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__))); from s2cfg import *   # g4: inputs by S2_CAT\n"
SPEC = {
 "completeness_ref/ref_size_check.py": [
  ('d = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", dtype={"id": str})', 'd = pd.read_csv(FEAT, dtype={"id": str})'),
  ('Cz = np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz")', 'Cz = np.load(GRID)'),
  ('json.dump(R, open(f"{HERE}/ref_size_record.json", "w"), indent=1)', 'json.dump(R, open(f"{HERE}/out/ref_size_record_{CAT}.json", "w"), indent=1)')],
 "completeness_ref/ref_size_insolation.py": [
  ('sys.path.insert(0, f"{HERE}/../insolation_v6")\nfrom common_null import hot_mask', 'pass'),
  ('c = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})', 'c = pd.read_csv(FEAT, low_memory=False, dtype={"id": str})'),
  ('Cz = np.load(f"{HE}/code/rerun/out_v2/completeness_10ppd.npz")', 'Cz = np.load(GRID)'),
  ('json.dump(R, open(f"{HERE}/ref_size_insolation.json", "w"), indent=1)', 'json.dump(R, open(f"{HERE}/out/ref_size_insolation_{CAT}.json", "w"), indent=1)')],
 "completeness_sfd/crater_completeness.py": [
  ('PF = json.load(open(f"{HERE}/theta_out_v2.json"))', 'PF = json.load(open(THETA))'),
  ('d = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", dtype={"id": str}, usecols=["id", "lat", "lon", "C1"])', 'd = pd.read_csv(FEAT, dtype={"id": str}, usecols=["id", "lat", "lon", "C1"])'),
  ('for a in (1.3, 1.6, 1.9):', 'for a in (1.3, 1.6, 1.7, 1.8, 1.9):'),
  ('d.drop(columns=["lat", "lon"]).to_csv(f"{HERE}/crater_completeness.csv", index=False)', 'd.drop(columns=["lat", "lon"]).to_csv(CCOMP, index=False)')],
 "completeness_sfd/hosting_fits.py": [
  ('sys.path.insert(0, f"{HERE}/../insolation_v6")\nfrom common_null import hot_mask', 'pass'),
  ('d = pd.read_csv(f"{HERE}/../summary_model/crater_features.csv", low_memory=False, dtype={"id": str})', 'd = pd.read_csv(FEAT, low_memory=False, dtype={"id": str})'),
  ('m = pd.read_csv(f"{HERE}/crater_completeness.csv", dtype={"id": str})', 'm = pd.read_csv(CCOMP, dtype={"id": str})'),
  ('for a in ("1.3", "1.6", "1.9"):', 'for a in ("1.3", "1.6", "1.7", "1.8", "1.9"):'),
  ('MEAS = ["C1", "lnC1", "C1_act", "lnE_1.6", "lnE_1.3", "lnE_1.9", "F_0.3"]', 'MEAS = ["C1", "lnC1", "C1_act", "lnE_1.6", "lnE_1.3", "lnE_1.9", "F_0.3", "lnE_1.7", "lnE_1.8"]'),
  ('json.dump(R, open(f"{HERE}/hosting_fits.json", "w"), indent=1)', 'json.dump(R, open(f"{HERE}/out/hosting_fits_{CAT}.json", "w"), indent=1)')],
 "completeness_sfd/joint_thinned.py": [
  ('sys.path.insert(0, f"{HERE}/../summary_model")', 'sys.path.insert(0, SM_CODE)'),
  ('m = pd.read_csv(f"{HERE}/crater_completeness.csv", dtype={"id": str})', 'm = pd.read_csv(CCOMP, dtype={"id": str})'),
  ('out_p = f"{HERE}/joint_thinned.json"', 'out_p = f"{HERE}/out/joint_thinned_{CAT}.json"'),
  ('m[["id", "C1_act", "E_1.3", "E_1.6", "E_1.9", "F_0.3"]', 'm[["id", "C1_act", "E_1.3", "E_1.6", "E_1.7", "E_1.8", "E_1.9", "F_0.3"]'),
  ('Z = np.load(os.path.expanduser("~/mnt/hollows-emissivity/code/rerun/out_v2/completeness_10ppd.npz"))', 'Z = np.load(GRID)'),
  ('"lnE_1.3": np.log(d["E_1.3"] + 0.01), "lnE_1.9": np.log(d["E_1.9"] + 0.01)', '"lnE_1.3": np.log(d["E_1.3"] + 0.01), "lnE_1.9": np.log(d["E_1.9"] + 0.01), "lnE_1.7": np.log(d["E_1.7"] + 0.01), "lnE_1.8": np.log(d["E_1.8"] + 0.01)')],
}
for rel, reps in SPEC.items():
    t = open(f"{W}/{rel}").read()
    for o, n in reps:
        assert t.count(o) == 1, (rel, o, t.count(o)); t = t.replace(o, n)
    lines = t.split("\n"); i = next(k for k, l in enumerate(lines) if l.startswith("import ") or l.startswith("from "))
    # the config import goes after the first import line, so HERE / HE of the original are defined by the script itself
    lines.insert(i + 1, HDR.rstrip("\n"))
    hdr = f"# g4 copy of manuscripts/PlanetaryResearch_hollows/analysis_workshop/{rel} (6 Nov 2026). Changed: input/output paths by S2_CAT, added exponents/measures 1.7 and 1.8 (see make_copies_s2.py)."
    open(f"{HERE}/{os.path.basename(rel)}", "w").write(hdr + "\n" + "\n".join(lines)); print("copied", rel, len(reps), "replacements")
