#!/usr/bin/env python3
"""Charged-layer thickness and charging temperature when a hollow forms by collapse of the charged layer (29 Oct 2026).
Hollow depth h = phi * d, where d = P_S(T) / (rho g) is the depth at which the vapor pressure of sulfur at the charging
temperature T equals the overburden of the whole column, and phi is the volume fraction of the layer lost on collapse
(native sulfur plus the pore space that closes). Vapor pressure as in analysis_v3/sulfur_charge.py (Section S16):
ln(P/p_b) = B (1/T_b - 1/T), T_b = 717.8 K, p_b = 101,325 Pa, B = 7,750-8,420 K; rho = 2,000-3,000 kg/m3; g = 3.70 m/s2.
Writes collapse_depth_record.json next to this script."""
import json, os, numpy as np
Tb, pb, g = 717.8, 101325.0, 3.70
RHO, BS = (2000.0, 3000.0), (7750.0, 8420.0)
def T_for(P, B): return 1.0 / (1.0 / Tb - np.log(P / pb) / B)
def P_at(T, B): return pb * np.exp(B * (1.0 / Tb - 1.0 / T))
R = {"assumptions": __doc__.split("\n")[1:5], "depth_targets_m": {"raw_median": 20.8, "corrected_medians_lt60": [24, 28]}, "by_phi": {}}
for phi in (1.0, 0.4, 0.3, 0.2, 0.1, 0.05):
    row = {}
    for name, D in (("20.8m", [20.8]), ("24-28m", [24.0, 28.0])):
        d = [x / phi for x in D]
        Ts = [T_for(rho * g * dd, B) for dd in d for rho in RHO for B in BS]
        row[name] = dict(d_m=[round(min(d), 1), round(max(d), 1)], T_K=[round(min(Ts)), round(max(Ts))],
                         P_MPa=[round(RHO[0] * g * min(d) / 1e6, 3), round(RHO[1] * g * max(d) / 1e6, 3)])
    Tmid = 0.5 * sum(row["24-28m"]["T_K"])
    row["K_drop_for_depth_ratio_0.77_at_Tmid"] = [round(Tmid - 1.0 / (1.0 / Tmid + np.log(1 / 0.77) / B), 1) for B in BS]
    row["depth_ratio_for_20K_drop_at_Tmid"] = [round(P_at(Tmid - 20, B) / P_at(Tmid, B), 3) for B in BS]
    row["rock_volume_heated_above_1200K_per_layer_volume_10volpctS"] = [8.0, 8.6]
    row["rock_column_m_for_24m_hollow_if_10volpctS"] = [round(8.0 * 24 / phi), round(8.6 * 24 / phi)]
    R["by_phi"][str(phi)] = row
R["critical_point_check"] = {"T_K": 1314, "P_measured_MPa": 20.7, "P_relation_MPa": [round(P_at(1314, B) / 1e6, 1) for B in BS]}
json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "collapse_depth_record.json"), "w"), indent=1)
for k, v in R["by_phi"].items(): print(k, v["24-28m"], v["K_drop_for_depth_ratio_0.77_at_Tmid"], v["rock_column_m_for_24m_hollow_if_10volpctS"])
print(R["critical_point_check"])
