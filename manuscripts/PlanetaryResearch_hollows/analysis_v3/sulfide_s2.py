"""Temperature at which each sulfide reaction reaches an S2 pressure of 1e-6 bar (0.1 Pa).
298 K standard enthalpies (kJ/mol) and entropies (J/mol/K), CRC/NIST values; heat-capacity changes neglected
(uncertainty ~+-100 K); pure phases at unit activity except where a_Si is given. Writes sulfide_s2_record.json."""
import json, numpy as np
from scipy.optimize import brentq
R = 8.314
H = {'FeS': -100.0, 'MgS': -346.0, 'CaS': -482.4, 'S2g': 128.6, 'Mgg': 147.1, 'Cag': 177.8, 'MgO': -601.6,
     'CaO': -634.9, 'SiO2': -910.7, 'Si': 0, 'Fe': 0}
S = {'FeS': 60.3, 'MgS': 50.3, 'CaS': 56.5, 'S2g': 228.2, 'Mgg': 148.6, 'Cag': 154.9, 'MgO': 26.9, 'CaO': 38.1,
     'SiO2': 41.5, 'Si': 18.8, 'Fe': 27.3}
def dG(prod, reac, T):
    h = sum(n*H[k] for k, n in prod.items()) - sum(n*H[k] for k, n in reac.items())
    s = sum(n*S[k] for k, n in prod.items()) - sum(n*S[k] for k, n in reac.items())
    return h*1e3 - T*s
p = 1e-6
cases = {
 'FeS -> Fe(s) + 1/2 S2': lambda T: dG({'Fe': 1, 'S2g': .5}, {'FeS': 1}, T) + R*T*np.log(p**.5),
 'MgS -> Mg(g) + 1/2 S2 (congruent)': lambda T: dG({'Mgg': 1, 'S2g': .5}, {'MgS': 1}, T) + R*T*np.log(2*p*p**.5),
 'CaS -> Ca(g) + 1/2 S2 (congruent)': lambda T: dG({'Cag': 1, 'S2g': .5}, {'CaS': 1}, T) + R*T*np.log(2*p*p**.5),
 'MgS + 1/2 SiO2 -> MgO + 1/2 Si + 1/2 S2': lambda T: dG({'MgO': 1, 'Si': .5, 'S2g': .5}, {'MgS': 1, 'SiO2': .5}, T) + R*T*np.log(p**.5),
 'MgS + 1/2 SiO2 -> MgO + 1/2 Si + 1/2 S2, a_Si=0.01': lambda T: dG({'MgO': 1, 'Si': .5, 'S2g': .5}, {'MgS': 1, 'SiO2': .5}, T) + R*T*np.log(p**.5*0.01**.5),
 'CaS + 1/2 SiO2 -> CaO + 1/2 Si + 1/2 S2': lambda T: dG({'CaO': 1, 'Si': .5, 'S2g': .5}, {'CaS': 1, 'SiO2': .5}, T) + R*T*np.log(p**.5),
 'CaS + 1/2 SiO2 -> CaO + 1/2 Si + 1/2 S2, a_Si=0.01': lambda T: dG({'CaO': 1, 'Si': .5, 'S2g': .5}, {'CaS': 1, 'SiO2': .5}, T) + R*T*np.log(p**.5*0.01**.5),
}
rec = {"pS2_bar": p, "T_K": {k: round(brentq(f, 300, 6000)) for k, f in cases.items()}}
# iron budget: S released per unit rock if all Fe is FeS; volume of source rock per volume of a 10 vol% native-S layer
Fe_wt = 1.5; S_wt = Fe_wt*32.06/55.85; rho_rock, rho_S = (2800, 3000), 2070
rec["iron_budget"] = {"Fe_wt_pct": Fe_wt, "max_S_wt_pct": round(S_wt, 2),
    "source_to_layer_volume_ratio_10volpct_S": [round(0.10*rho_S/(r*S_wt/100), 1) for r in rho_rock]}
# depth where S2 vapor pressure equals overburden, g = 3.7 m/s2; ln p(Pa) = 24.40978 - 9028.66789/T
rec["vapor_overburden_depth_m"] = {f"{T}K_rho{r}": round(np.exp(24.40978 - 9028.66789/T)/(r*3.7), 1)
                                   for T in (718, 750, 800) for r in (2000, 3000)}
rec["mercury_day_lowering"] = {"m_per_1e5yr": 2.0, "um_per_solar_day": round(2.0/1e5*175.94/365.25*1e6, 1)}
json.dump(rec, open("sulfide_s2_record.json", "w"), indent=1); print(json.dumps(rec, indent=1))
