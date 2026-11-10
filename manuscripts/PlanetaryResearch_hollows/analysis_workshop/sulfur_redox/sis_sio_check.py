# Check (6 Nov 2026): SiS(g) route vs SiO(g) route for S2 release from CaS/MgS + SiO2; MgSiO3 product; SiO and S2
# condensation pressures. NASA polynomials shipped with Cantera (same data as redox_check.py). Pure phases, unit activities.
import cantera as ct, numpy as np, re
sp = {s.name: s for s in ct.Species.list_from_file('nasa_gas.yaml')}
sc = {s.name: s for s in ct.Species.list_from_file('nasa_condensed.yaml')}
R=8.314
def gas(s,T): th=sp[s].thermo; return (th.h(T)-T*th.s(T))/1000
def cond(base,T):
    best=None
    for n,s in sc.items():
        if re.fullmatch(re.escape(base)+r"\((?:[a-z]+\d*|[A-Za-z]{1,4}|I+)\)", n):
            th=s.thermo
            if th.min_temp<=T<=th.max_temp:
                v=(th.h(T)-T*th.s(T))/1000; best=v if best is None else min(best,v)
    return best
P=lambda dG,T: (np.exp(-dG/(R*T))/2)**(2/3)*1e5   # gas pair with p1 = 2 pS2: K = 2 pS2^1.5 (bar)
for T in (1600,1800,2000,2300):
    out=[f"{T} K"]
    for sulf,ox in (("CaS","CaO"),("MgS","MgO")):
        dG=cond(ox,T)+gas('SiO',T)+0.5*gas('S2',T)-cond(sulf,T)-cond('SiO2',T)
        dG2=2*cond(ox,T)+gas('SiS',T)+0.5*gas('S2',T)-2*cond(sulf,T)-cond('SiO2',T)
        out.append(f"{sulf}: via SiO {P(dG,T):.2g} Pa, via SiS {P(dG2,T):.2g} Pa")
    m=cond('MgSiO3',T)
    if m is not None:
        dG3=m+gas('SiO',T)+0.5*gas('S2',T)-cond('MgS',T)-2*cond('SiO2',T)
        out.append(f"MgS+2SiO2->MgSiO3+SiO+S2/2: {P(dG3,T):.2g} Pa")
    print(" | ".join(out))
for T in (1000,1100,1200,1300,1400,1500,1600):
    dG=2*gas('SiO',T)-cond('SiO2',T)-cond('Si',T); print(T,"K pSiO over SiO2+Si", f"{np.sqrt(np.exp(-dG/(R*T)))*1e5:.2g} Pa")
for T in (350,400,450,500,600,700):
    c=cond('S',T)
    if c is None: print(T,"no condensed S"); continue
    dG=0.5*gas('S2',T)-c; print(T,"K pS2 over condensed S", f"{np.exp(-dG/(R*T))**2*1e5:.2g} Pa")
