"""Numbers behind the revised fumarole (sulfur charging and cycling) model in Section 4 (manuscript v4).
1. Sulfur vapor pressure near and above its normal boiling point, anchored at the normal boiling point
   (717.8 K, 101.325 kPa; Lide 2007) with a Clausius-Clapeyron slope B between 7,750 K (fit to the CRC
   tabulation at 591-717 K) and 8,420 K (chord from the boiling point to the critical point, 1,313 K and
   20.7 MPa; Baker 1971). The ln p = 24.40978 - 9028.67/T fit used for sublimation fronts at 200-330 K
   overestimates p by ~1.5x above 750 K, so it is not used here.
2. Depth at which that pressure equals the overburden, d = p/(rho g), g = 3.70 m/s2, rho = 2000-3000 kg/m3.
3. S2 release by sulfide breakdown: temperature at which pS2 = 1e-6 bar (sulfide_s2.py, standard-state data).
4. Iron budget: rock with 1.5 wt% Fe (Weider et al. 2014), all as FeS, releases at most 0.86 wt% S.
5. Daily cycles per Gyr and the lowering per Mercury solar day needed for ~2 m per 1e5 yr.
Writes sulfur_charge_record.json."""
import json, numpy as np
Tb, pb = 717.8, 101325.0
g, rhos = 3.70, (2000.0, 3000.0)
R = {"anchor": {"Tb_K": Tb, "pb_Pa": pb, "B_K": [7750, 8420]}}
def p(T, B): return pb*np.exp(B*(1/Tb - 1/T))
R["vapor_pressure_MPa"] = {f"{T}K": [round(p(T, B)/1e6, 3) for B in (7750, 8420)] for T in (717.8, 750, 775, 800, 850)}
R["depth_m"] = {f"{T}K": [round(min(p(T, B) for B in (7750, 8420))/(g*rhos[1]), 1),
                          round(max(p(T, B) for B in (7750, 8420))/(g*rhos[0]), 1)] for T in (717.8, 750, 775, 800, 850)}
R["fit_used_for_fronts_MPa"] = {f"{T}K": round(np.exp(24.40978 - 9028.66789/T)/1e6, 3) for T in (750, 800)}
R["S2_release_T_K_at_1e-6bar"] = {"FeS": 1186, "MgS_congruent_Mg_gas": 1471, "MgS_with_SiO2": [1693, 1930],
                                  "CaS_congruent_Ca_gas": 1911, "CaS_with_SiO2": [2282, 2590], "source": "sulfide_s2_record.json"}
Fe = 1.5; S = Fe*32.06/55.85
R["iron_budget"] = {"Fe_wt_pct": Fe, "max_S_wt_pct": round(S, 2),
                    "rock_volume_per_layer_volume_10volpct_S": [round(0.10*2070/(r*S/100), 1) for r in (3000, 2800)]}
day_yr = 175.94/365.25
R["solar_day_yr"] = round(day_yr, 4); R["cycles_per_Gyr"] = round(1e9/day_yr, -6)
R["lowering_um_per_day_for_2m_per_1e5yr"] = round(2.0/1e5*day_yr*1e6, 1)

# charging temperature implied by a depth d: p(T) = rho g d, over rho = 2000-3000 and B = 7750-8420
def T_for(d):
    out = []
    for r in rhos:
        for B in (7750, 8420):
            P = r*g*d; out.append(1/(1/Tb - np.log(P/pb)/B))
    return [round(min(out)), round(max(out))]
R["charging_T_K_for_depth"] = {"20.8m_raw_median": T_for(20.8), "24m": T_for(24.0), "28m": T_for(28.0)}
# depth ratio for a 20 K lower charging temperature near 770 K
R["depth_ratio_minus20K_at_770K"] = [round(np.exp(B*(1/770 - 1/750)), 2) for B in (7750, 8420)]
# S2 partial pressure over FeS + Fe at 1190 K (from the 1e-6 bar definition in sulfide_s2.py)
R["FeS_pS2_at_1186K_Pa"] = 0.1
json.dump(R, open("sulfur_charge_record.json", "w"), indent=1); print(json.dumps(R, indent=1))
