#!/usr/bin/env python3
"""Electron sinks for making native sulfur from sulfide (5 Nov 2026). Equilibrium S2 pressures over pure phases from the
NASA polynomial data shipped with Cantera (nasa_gas.yaml, nasa_condensed.yaml; pip install cantera). Unit activities, so
order of magnitude only. Check: FeS -> Fe + S2/2 gives 0.041 Pa at 1,300 K and 1.1 Pa at 1,471 K (JANAF in S13: 1 Pa at 1,471 K).
(1) Fe sink: FeS -> Fe + 1/2 S2. (2) Si sink: MS + SiO2 -> MO + SiO(g) + 1/2 S2 for M = Fe, Ca, Mg, with P_SiO = 2 P_S2.
(3) Carbon: COS -> CO + 1/2 S2 (K = 1 near 1,180 K); C + S2 -> CS2 (K 4-18 at 600-2,000 K); CS2(g) -> C + 2 S(cond),
favorable below ~750 K."""
# 9 Nov 2026: the printed pressures are also written to redox_check_record.json (Pa, two significant figures as printed); computation unchanged.
import cantera as ct, numpy as np, re, json
REC = {"source": "Cantera nasa_gas.yaml / nasa_condensed.yaml, cantera " + ct.__version__, "units": "Pa (S2); kJ/mol (dG)", "S2_Pa": {}, "CS2_to_C_2S_dG_kJ": {}}
sp = {s.name: s for s in ct.Species.list_from_file('nasa_gas.yaml')}
sc = {s.name: s for s in ct.Species.list_from_file('nasa_condensed.yaml')}
R = 8.314
def gas(s, T): th = sp[s].thermo; return (th.h(T) - T * th.s(T)) / 1000
def cond(base, T):
    best = None
    for n, s in sc.items():
        if re.fullmatch(re.escape(base) + r"\((?:[a-z]+\d*|[A-Za-z]{1,4})\)", n):
            th = s.thermo
            if th.min_temp <= T <= th.max_temp:
                v = (th.h(T) - T * th.s(T)) / 1000; best = v if best is None else min(best, v)
    return best
for T in (1300, 1471, 1600, 1800, 2000, 2300, 2600):
    p = np.exp(-2*(cond('Fe',T)+0.5*gas('S2',T)-cond('FeS',T))/(R*T))*1e5; REC["S2_Pa"][str(T)] = {"FeS->Fe+S2/2": float(f"{p:.2g}")}
    out = [f"{T} K | FeS->Fe+S2/2: {p:.2g} Pa"]
    for sulf, ox in (("FeS", "FeO"), ("CaS", "CaO"), ("MgS", "MgO")):
        dG = cond(ox, T) + gas('SiO', T) + 0.5 * gas('S2', T) - cond(sulf, T) - cond('SiO2', T)
        q = (np.exp(-dG/(R*T))/2)**(2/3)*1e5; REC["S2_Pa"][str(T)][f"{sulf}+SiO2->{ox}+SiO+S2/2"] = float(f"{q:.2g}")
        out.append(f"{sulf}+SiO2->{ox}+SiO+S2/2: {q:.2g} Pa")
    print(" | ".join(out))
for T in (400, 600, 700, 800):
    print(T, "CS2 -> C + 2S(cond) dG kJ/mol", round((cond('C', T) + 2 * cond('S', T) - gas('CS2', T)) / 1000, 1))
    REC["CS2_to_C_2S_dG_kJ"][str(T)] = round((cond('C', T) + 2 * cond('S', T) - gas('CS2', T)) / 1000, 1)
json.dump(REC, open("redox_check_record.json", "w"), indent=1)
