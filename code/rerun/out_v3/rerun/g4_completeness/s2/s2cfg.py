"""g4 (6 Nov 2026): input selection for the S2 copies (env S2_CAT = v2 | v3 | v3r | v3br)."""
import os as _os, numpy as _np
_RT = _os.path.expanduser("~/mnt/hollows-emissivity"); _G = f"{_RT}/code/rerun/out_v3/rerun"; _H = _os.path.dirname(_os.path.abspath(__file__))
_W = f"{_RT}/manuscripts/PlanetaryResearch_hollows/analysis_workshop"
CAT = _os.environ.get("S2_CAT", "v3r"); assert CAT in ("v2", "v3", "v3b", "v3r", "v3br"), CAT
FEAT = {"v2": f"{_W}/summary_model/crater_features.csv", "v3": f"{_G}/g1_graph_joint/summary_model/v3/crater_features.csv",
        "v3b": f"{_G}/g1_graph_joint/summary_model/v3b/crater_features.csv",
        "v3r": f"{_G}/g4_completeness/g1/summary_model/v3r/crater_features.csv", "v3br": f"{_G}/g4_completeness/g1/summary_model/v3br/crater_features.csv"}[CAT]
GRID = {"v3r": f"{_G}/g4_completeness/completeness_10ppd_v3.npz", "v3br": f"{_G}/g4_completeness/completeness_10ppd_v3b.npz"}.get(CAT, f"{_RT}/code/rerun/out_v2/completeness_10ppd.npz")
THETA = {"v3r": f"{_G}/g3_catalog/sfd/out_v3_refit/theta_out_v2.json", "v3br": f"{_G}/g3_catalog/sfd/out_v3b_refit/theta_out_v2.json"}.get(CAT, f"{_W}/completeness_sfd/theta_out_v2.json")
GTAG = {"v3r": "v3refit", "v3br": "v3brefit"}.get(CAT, "v2model")        # crater_completeness depends only on the detection model
CCOMP = f"{_H}/out/crater_completeness_{GTAG}.csv"
SM_CODE = f"{_G}/g4_completeness/g1/summary_model"                          # g4 copy of step7_improved.py (verbatim)
_os.environ.setdefault("SM_DIR", {"v2": f"{_G}/g1_graph_joint/summary_model/v2", "v3": f"{_G}/g1_graph_joint/summary_model/v3", "v3b": f"{_G}/g1_graph_joint/summary_model/v3b",
                                  "v3r": f"{_G}/g4_completeness/g1/summary_model/v3r", "v3br": f"{_G}/g4_completeness/g1/summary_model/v3br"}[CAT])
_os.makedirs(f"{_H}/out", exist_ok=True)
def hot_mask(lon):   # verbatim from analysis_workshop/insolation_v6/common_null.py
    lon = _np.mod(lon, 360.0)
    return (_np.minimum(_np.abs(lon), _np.abs(lon - 360)) <= 45.0) | (_np.abs(lon - 180) <= 45.0)
