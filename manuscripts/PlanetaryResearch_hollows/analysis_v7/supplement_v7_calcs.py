#!/usr/bin/env python3
"""Numbers for Sections S10, S14, and S16 of supplement v7 (30 Oct 2026).

Writes supplement_v7_calcs_record.json next to this script. Inputs:
  - NIST-JANAF table Fe-025, troilite FeS(cr,l), delta-f G (kJ/mol), https://janaf.nist.gov/tables/Fe-025.txt.
    JANAF takes S2(g) as the sulfur reference state above 882 K, so for Fe(s) + 1/2 S2 = FeS,
    p(S2) [bar] = exp(2 dfG / RT).
  - collapse_depth_record.json (this folder): charging temperatures and depth ratios for h = phi d.
  - depth_reanalysis_record.json (mercury-hollows, GRL_insolation_hollows/records): band medians of mean annual
    temperature and counts of measured hollows, and the depth regression on 1/T.
  - sulfur_charge_record.json (analysis_v3): vapor-pressure anchor and slopes B; iron budget.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
R = 8.314462618

# --- 1. FeS -> Fe + 1/2 S2 from NIST-JANAF Fe-025 -------------------------------------------------------------
JANAF_FeS_dfG_kJ = {900: -103.437, 1000: -97.670, 1100: -91.766, 1200: -85.786, 1300: -79.802,
                    1400: -73.918, 1500: -68.943, 1600: -65.382}


def pS2_Pa(T):
    return math.exp(2 * JANAF_FeS_dfG_kJ[T] * 1e3 / (R * T)) * 1e5


p_table = {T: pS2_Pa(T) for T in JANAF_FeS_dfG_kJ}


def T_at(p_target):
    """Interpolate ln p linearly in 1/T between tabulated temperatures."""
    Ts = sorted(p_table)
    for a, b in zip(Ts, Ts[1:]):
        if p_table[a] <= p_target <= p_table[b]:
            x = (math.log(p_target) - math.log(p_table[a])) / (math.log(p_table[b]) - math.log(p_table[a]))
            return 1 / (1 / a + x * (1 / b - 1 / a))
    return None


fes = {"pS2_Pa": {str(T): float(f"{p:.3g}") for T, p in p_table.items()},
       "T_at_0.1Pa_K": round(T_at(0.1)), "T_at_1Pa_K": round(T_at(1.0)),
       "T_298K_estimate_K": 1186,
       "note": "FeS melts near 1,463 K; JANAF rows at 1,500-1,600 K are for the liquid."}

# --- 2. Condensation temperature of sulfur vapor at 0.1 and 1 Pa -----------------------------------------------
# 1 Pa at 375 K is the lowest tabulated vapor pressure of sulfur (CRC Handbook). Below it, extrapolate with the
# least-squares slope of the six tabulated values (9,029 K) and with the slope of the two lowest values (10,684 K).
cond = {}
for B in (9029, 10684):
    cond[str(B)] = round(1 / (1 / 375 + math.log(10) / B), 1)
condensation = {"T_at_1Pa_K": 375, "T_at_0.1Pa_K_by_slope": cond}

# --- 3. Sunlight and depth ----------------------------------------------------------------------------------------
dr = json.load(open(os.path.join(HERE, "inputs", "depth_reanalysis_record.json")))
b = dr["bands"]
eq = [b[k] for k in ("0-30", "30-45", "45-60")]
po = [b[k] for k in ("60-75", "75-90")]
T_eq = sum(x["n"] * x["T_med"] for x in eq) / sum(x["n"] for x in eq)
T_po = sum(x["n"] * x["T_med"] for x in po) / sum(x["n"] for x in po)
coef = dr["T_px_inc"]["x"]

cd = json.load(open(os.path.join(HERE, "collapse_depth_record.json")))
Bs = (7750, 8420)
T_s = 290.0  # n-weighted band-median mean annual temperature of the measured hollows, rounded
pred = {}
for phi in ("0.4", "0.2", "0.1"):
    Tlo, Thi = cd["by_phi"][phi]["24-28m"]["T_K"]
    Tmid = 0.5 * (Tlo + Thi)
    # If the charging temperature includes the local mean surface temperature one for one, ln h = const - B/Tc with
    # Tc = Tc0 + Ts, so d ln h / d(1/Ts) = -B Ts^2 / Tc^2.
    pred[phi] = {"Tc_mid_K": round(Tmid), "coef_K": [round(-B * T_s ** 2 / Tmid ** 2) for B in Bs],
                 "K_drop_for_ratio_0.77": cd["by_phi"][phi]["K_drop_for_depth_ratio_0.77_at_Tmid"]}
# Depth ratio if the charging temperature falls by the full surface-temperature difference, h ~ P_S(T) with
# ln(P/p_b) = B (1/T_b - 1/T), T_b = 717.8 K.
dTs = T_eq - T_po
for phi in ("0.4", "0.2", "0.1"):
    Tmid = pred[phi]["Tc_mid_K"]
    pred[phi]["depth_ratio_for_full_dTs"] = [round(math.exp(-B * (1 / (Tmid - dTs) - 1 / Tmid)), 3) for B in Bs]
sunlight = {"T_mean_equatorward_K": round(T_eq, 2), "T_mean_poleward_K": round(T_po, 2),
            "difference_K": round(T_eq - T_po, 1),
            "measured_coef_K": [round(coef[0]), round(coef[1])],
            "predicted_coef_if_Tc_includes_Ts": pred,
            "misfit_in_SE": [round((min(min(v["coef_K"]) for v in pred.values()) - coef[0]) / coef[1], 1),
                             round((max(max(v["coef_K"]) for v in pred.values()) - coef[0]) / coef[1], 1)]}

# --- 4. Iron budget ----------------------------------------------------------------------------------------------
sc = json.load(open(os.path.join(HERE, "inputs", "sulfur_charge_record.json")))
Fe = sc["iron_budget"]["Fe_wt_pct"]
S_max = Fe * 32.06 / 55.845            # wt% S if all Fe is FeS
rho_S, rho_rock = 2.07, (2.8, 3.0)      # g/cm3: native sulfur; rock
ratio = [0.10 * rho_S / (S_max / 100 * r) for r in rho_rock]
iron = {"Fe_wt_pct": Fe, "S_max_wt_pct": round(S_max, 2), "rho_native_S": rho_S, "rho_rock": rho_rock,
        "rock_volume_per_layer_volume_10volpct_S": sorted(round(x, 1) for x in ratio),
        "layer_m_for_24m_hollow_phi0.1": 240,
        "rock_column_m": sorted(round(240 * x) for x in ratio)}

# --- 5. Floor lowering -------------------------------------------------------------------------------------------
day_yr = sc["solar_day_yr"]
lower = {"m_per_1e5yr": 2.0, "um_per_solar_day": round(2.0 / 1e5 * day_yr * 1e6, 1),
         "yr_for_24m": 24 / 2.0 * 1e5, "cycles_per_Gyr": sc["cycles_per_Gyr"]}

rec = {"FeS_JANAF": fes, "sulfur_condensation": condensation, "sunlight_and_depth": sunlight,
       "iron_budget": iron, "floor_lowering": lower}
json.dump(rec, open(os.path.join(HERE, "supplement_v7_calcs_record.json"), "w"), indent=1)
print(json.dumps(rec, indent=1))
