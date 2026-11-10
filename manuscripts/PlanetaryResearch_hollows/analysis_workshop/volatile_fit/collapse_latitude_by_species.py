"""Supplement to volatile_candidates.py. In the collapse model, ln h = const + ln P_v(T_base). If the base temperature rises one
for one with the local mean annual surface temperature T_s, the coefficient of ln h on 1/T_s is -B_eff T_s^2 / T_base^2, where
B_eff = -d ln P_v / d(1/T) at T_base; with half coupling (contact temperature, equal thermal inertias) it is half that.
The poleward depth ratio for a 21.5 K lower T_s is exp(-B_eff c dT_s / T_base^2), c = 1 or 0.5. Measured: -291 +/- 99 K; 0.77
(0.69-0.86). T_s = 290 K. Uses the vapor-pressure functions of volatile_candidates.py at phi = 0.4, 0.2, 0.1. Writes
collapse_latitude_by_species_record.json."""
import json, os, math, importlib.util, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("vc", os.path.join(HERE, "volatile_candidates.py"))
vc = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(vc)
Ts, dTs = 290.0, 21.5
out = {}
for nm, fs in vc.Pmeth.items():
    out[nm] = {}
    for phi in vc.phis:
        coef, rat = [], []
        for f in fs:
            for r in vc.rhos:
                for h in vc.hs:
                    P = r * vc.g * h / phi; Ta, Tb = f(P), f(1.01 * P)
                    B = math.log(1.01) / (1 / Ta - 1 / Tb)
                    for c in (1.0, 0.5):
                        coef.append((c, -c * B * Ts ** 2 / Ta ** 2)); rat.append((c, math.exp(-B * c * dTs / Ta ** 2)))
        out[nm][str(phi)] = {f"coef_c{c}": [round(min(v for cc, v in coef if cc == c)), round(max(v for cc, v in coef if cc == c))] for c in (1.0, 0.5)}
        out[nm][str(phi)].update({f"ratio_c{c}": [round(min(v for cc, v in rat if cc == c), 3), round(max(v for cc, v in rat if cc == c), 3)] for c in (1.0, 0.5)})
json.dump(out, open(os.path.join(HERE, "collapse_latitude_by_species_record.json"), "w"), indent=1)
for k, v in out.items(): print(k, v)
