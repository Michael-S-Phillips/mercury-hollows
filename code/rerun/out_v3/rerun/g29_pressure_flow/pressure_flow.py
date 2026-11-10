#!/usr/bin/env python3
"""g29 (10 Nov 2026): numbers for the "Pressure-driven flow" paragraph of supplement S13 (v28) and its row in Table S14, and the
2 kPa / 2,150 K statement in main text 4.1 and 4.3 (v6.16). Michael's comment on v27 p. 38: sulfur cannot outrun heat by diffusion, so
explain pressure-driven flow through fractures. Inputs are the constants of S13: S2 at 1,000 K; 100 m of regolith with 5 vol% native
sulfur (1.0e4 kg m-2, the "Sulfur delivery" paragraph); heat crossing 100 m in L^2/kappa = 317 yr (kappa = 1e-6 m2 s-1); FeS S2 pressure
5.4 Pa at 1,600 K and an enthalpy of 314 kJ mol-1 (Table S16); MgS + SiO2 -> MgO + SiO(g) + 1/2 S2 at 350 Pa (2,000 K) and 7,000 Pa
(2,300 K) (Table S15); native sulfur ln(P/101,325 Pa) = B (1/717.8 - 1/T), B = 7,750-8,420 K (S13). Permeabilities 1e-12 to 1e-14 m2
from Alexander et al. (2024, doi:10.1029/2023EA003383) for the Chicxulub structure; open fractures 1 mm wide at 1% of the volume with a
tortuosity of 3 (cubic law, k = eps a^2 / (12 tau)), as in S13. Viscosity of S2 from kinetic theory (hard spheres, d = 3.6 A).
Viscous (Darcy) flux of an ideal gas from a source at p into p ~ 0: J = k M p^2 / (2 mu R T L). Delivery time = mass / J, which scales
as L^2 like the heat time, so the condition p >= p_min does not depend on L. Writes pressure_flow_record.json."""
import json, numpy as np
kB, R, NA = 1.380649e-23, 8.314, 6.02214e23
M, d = 0.064, 3.6e-10
m = M / NA
T = 1000.0
mu = 5 / 16 * np.sqrt(np.pi * m * kB * T) / (np.pi * d ** 2)
lam = lambda p, T=T: kB * T / (np.sqrt(2) * np.pi * d ** 2 * p)
L, mass, kappa = 100.0, 1.0e4, 1e-6
t_heat = L ** 2 / kappa
def p_min(k):
    return float(np.sqrt(mass * 2 * mu * R * T * L / (k * M * t_heat)))
k_frac = 0.01 * (1e-3) ** 2 / (12 * 3)
T_FeS = lambda P: float(1 / (1 / 1600 - np.log(P / 5.4) * R / 314000))
_lnP, _iT = np.log([350.0, 7000.0]), np.array([1 / 2000, 1 / 2300])
T_MgS = lambda P: float(1 / (_iT[0] + (np.log(P) - _lnP[0]) / (_lnP[1] - _lnP[0]) * (_iT[1] - _iT[0])))
T_S = lambda P, B: float(1 / (1 / 717.8 - np.log(P / 101325) / B))
rec = dict(created="2026-11-10", note=__doc__,
           mu_S2_1000K_Pa_s=float(mu), heat_time_100m_yr=t_heat / 3.156e7,
           mean_free_path_1000K_m={"1 Pa": lam(1.0), "5.4 Pa": lam(5.4)},
           p_min_kPa={f"{k:.0e}": p_min(k) / 1e3 for k in (1e-12, 1e-13, 1e-14, 1e-15)},
           k_open_fractures_m2=k_frac, p_min_open_fractures_kPa=p_min(k_frac) / 1e3,
           mean_free_path_at_half_pmin_um={"1e-12": lam(p_min(1e-12) / 2) * 1e6, "1e-14": lam(p_min(1e-14) / 2) * 1e6, "open fractures": lam(p_min(k_frac) / 2) * 1e6},
           T_reach_K={"FeS (extrapolated), open fractures": T_FeS(p_min(k_frac)), "MgS + SiO2, open fractures": T_MgS(p_min(k_frac)),
                      "native S, open fractures, B 7,750 / 8,420": [T_S(p_min(k_frac), 7750), T_S(p_min(k_frac), 8420)],
                      "native S, k 1e-12, B 7,750 / 8,420": [T_S(p_min(1e-12), 7750), T_S(p_min(1e-12), 8420)],
                      "native S, k 1e-14, B 7,750 / 8,420": [T_S(p_min(1e-14), 7750), T_S(p_min(1e-14), 8420)],
                      "FeS (extrapolated), k 1e-12": T_FeS(p_min(1e-12))},
           mean_impact_velocity_km_s={"Mercury": 42.5, "Earth": 20.4, "source": "Le Feuvre and Wieczorek (2008), Table 2"})
json.dump(rec, open("pressure_flow_record.json", "w"), indent=1)
print(json.dumps({k: v for k, v in rec.items() if k != "note"}, indent=1))
