"""Tasks 1b and 1c. Candidate volatiles: (b) enthalpy that sets the temperature dependence of each vapor pressure, and the
depth ratio a sunlight-driven front of each would predict across 219-324 K; (c) base temperature at which each vapor pressure
equals the overburden rho g d of the charged layer beneath a 24-28 m hollow, d = h/phi (complete collapse, f = 0).
Data sources (NIST Chemistry WebBook pages reproduce NIST-JANAF; Chase 1998, J. Phys. Chem. Ref. Data Monograph 9), all
accessed 2 November 2026; values are quoted in SOURCES below. Writes volatile_candidates_record.json and volatile_candidates.csv."""
import json, os, math, numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
R = 8.314462618
fit = json.load(open(os.path.join(HERE, "depth_T_fit_record.json")))
T0, T1 = fit["T_range"]; dx = 1 / T0 - 1 / T1
Tp0, Tp1 = fit["peak_noon"]["T_range"]; dxp = 1 / Tp0 - 1 / Tp1
g, rhos, hs, phis = 3.70, (2000.0, 3000.0), (24.0, 28.0), (0.4, 0.2, 0.1)

SOURCES = {
 "S8(g) dfH298 = 100.416 kJ/mol": "NIST-JANAF table S-021, S8(g) (janaf.nist.gov/tables/S-021.html)",
 "S2(g) dfH298 = 128.60 kJ/mol, S298 = 228.19; Shomate 298-6000 K": "NIST WebBook (C23550450), Chase 1998",
 "S(l) dfH298 = 1.85 kJ/mol": "NIST WebBook (C7704349), Chase 1998",
 "S liquid vapor pressure fit ln p = 24.41 - 9029/T (375-717 K); two lowest points 10,684 K; above 718 K B = 7,750-8,420 K": "Section S9/S13 of the supplement (lide2007 = CRC Handbook, 88th ed.)",
 "NaCl(cr) dfH298 = -411.12, S = 72.11; NaCl(l) -385.92, 95.06; NaCl(g) -181.42, 229.79": "NIST WebBook (C7647145), Chase 1998",
 "NaCl Antoine 1138-1738 K: A 5.07184, B 8388.497, C -82.638 (bar, K)": "NIST WebBook (C7647145), Stull 1947",
 "KCl(cr) dfH298 = -436.68, S = 82.56; KCl(l) -421.79, 86.66; KCl(g) -214.68, 239.09": "NIST WebBook (C7447407), Chase 1998",
 "KCl Antoine 1094-1680 K: A 4.78236, B 7440.691, C -122.709": "NIST WebBook (C7447407), Stull 1947",
 "Na2S(cr) dfH298 = -366.10, Shomate 298-1445 K": "NIST WebBook (C1313822), Chase 1998",
 "K2S(cr) dfH298 = -376.56, Shomate 298-1221 K": "NIST WebBook (C1312738), Chase 1998",
 "Na(g) dfH298 = 107.30, S = 153.65; Shomate 1170-6000 K": "NIST WebBook (C7440235), Chase 1998",
 "K(g) dfH298 = 89.00, S = 160.34; Shomate 1039-1800, 1800-6000 K": "NIST WebBook (C7440097), Chase 1998",
 "MgS(cr) dfH298 = -345.72, S = 50.30, Shomate 298-3000 K; MgS(g) 145.23, 225.49": "NIST WebBook (C12032369), Chase 1998",
 "CaS(cr) dfH298 = -473.21, S = 56.54, Shomate 298-3000 K; CaS(g) 123.59, 232.59": "NIST WebBook (C20548543), Chase 1998",
 "Mg(g) dfH298 = 147.10, S = 148.65; Shomate 1366-2200, 2200-6000 K": "NIST WebBook (C7439954), Chase 1998",
 "Ca(g) dfH298 = 177.80, S = 154.89; Shomate 1774-6000 K": "NIST WebBook (C7440702), Chase 1998",
 "FeS dfG (S2 reference) 900-1600 K": "NIST-JANAF table Fe-025 (janaf_fe025), as tabulated in analysis_v7/supplement_v7_calcs.py",
 "C(g) dfH298 = 716.669, S = 158.1; Shomate 298-6000 K": "NIST WebBook (C7440440), Chase 1998; graphite S298 = 5.74 J/mol/K (standard value, not fetched)",
 "H2O(g) dfH298 = -241.83, H2O(l) -285.83 kJ/mol": "NIST WebBook (C7732185), Chase 1998; fusion enthalpy 6.01 kJ/mol (standard value, not fetched)",
 "H2O Antoine 379-573 K: A 3.55959, B 643.748, C -198.043": "NIST WebBook (C7732185), Liu and Lindsay 1970",
}
# ---- Shomate helpers: H_abs (kJ/mol, on the dfH298 scale) and S (J/mol/K) ------------------------------------
def shom(p, T):
    A, B, C, D, E, F, G, H = p; t = T / 1000.0
    Habs = A * t + B * t * t / 2 + C * t ** 3 / 3 + D * t ** 4 / 4 - E / t + F
    S = A * math.log(t) + B * t + C * t * t / 2 + D * t ** 3 / 3 - E / (2 * t * t) + G
    return Habs, S
def Gabs(ranges, T):
    """ranges: list of (Tlo, Thi, params). Outside the ranges, extend the nearest range with its Cp held constant."""
    for lo, hi, p in ranges:
        if lo - 1e-6 <= T <= hi + 1e-6: H, S = shom(p, T); return H * 1e3 - T * S
    lo, hi, p = ranges[-1] if T > ranges[-1][1] else ranges[0]
    Te = hi if T > hi else lo
    H, S = shom(p, Te); t = Te / 1000
    A, B, C, D, E = p[:5]; cp = A + B * t + C * t * t + D * t ** 3 + E / (t * t)
    H += cp * (T - Te) / 1e3; S += cp * math.log(T / Te)
    return H * 1e3 - T * S
SH = {
 "S2g":  [(298, 6000, (33.51313, 5.065360, -1.059670, 0.089905, -0.211911, 117.6855, 266.0919, 128.6003))],
 "Nag":  [(1170.525, 6000, (20.80573, 0.277206, -0.392086, 0.119634, -0.008879, 101.0386, 178.7095, 107.2999))],
 "Kg":   [(1039.54, 1800, (20.66122, 0.391869, -0.417344, 0.145582, 0.003764, 82.83860, 185.2650, 88.99996)),
          (1800, 6000, (58.70570, -27.38277, 6.730509, -0.420844, -25.87921, 32.37931, 198.3777, 88.99996))],
 "Mgg":  [(1366.104, 2200, (20.77306, 0.035592, -0.031917, 0.009109, 0.000461, 140.9071, 173.7799, 147.1002)),
          (2200, 6000, (47.60848, -15.40875, 2.875965, -0.120806, -27.01764, 97.40017, 177.2305, 147.1002))],
 "Cag":  [(1774, 6000, (121.5470, -74.95390, 19.17230, -1.400821, -64.51340, 42.23540, 217.4470, 177.8000))],
 "MgSs": [(298, 3000, (48.60009, 3.685472, 0.021388, 0.015067, -0.369222, -361.6286, 105.9347, -345.7243))],
 "MgSg": [(298, 800, (51.11480, -139.9390, 360.4070, -223.8690, -0.078348, 133.2030, 314.5860, 145.2270)),
          (800, 2000, (199.0350, -190.3110, 84.05370, -13.15290, -24.82510, 25.63380, 421.0220, 145.2270)),
          (2000, 6000, (34.52960, 1.196611, 0.015877, -0.004516, 25.05680, 168.0760, 291.2810, 145.2270))],
 "CaSs": [(298, 3000, (48.60260, 7.280161, -1.933724, 0.329739, -0.289224, -488.9966, 111.6450, -473.2104))],
 "CaSg": [(298, 1400, (33.08887, 15.89464, -21.61998, 10.84635, -0.123820, 112.7768, 268.0601, 123.5949)),
          (1400, 3000, (-247.8309, 258.9883, -74.23546, 6.965901, 106.5113, 395.8616, 104.5720, 123.5949)),
          (3000, 6000, (254.7052, -86.85356, 13.05132, -0.689578, -287.7881, -317.4526, 244.2067, 123.5949))],
 "Na2Ss":[(298, 1000, (78.38780, 13.94820, 2.753281, -2.347990, 0.008301, -390.0830, 186.8820, -366.1000)),
          (1000, 1276, (58023.20, -104689.0, 70037.20, -16374.80, -6903.890, -32139.00, 71878.20, -366.1000)),
          (1276, 1445, (-542020.0, 789111.0, -431902.0, 84244.40, 103416.0, 373259.0, -549542.0, -366.1000))],
 "K2Ss": [(298, 800, (69.17220, 16.75400, 13.43840, -6.573700, -0.045088, -398.1860, 192.9800, -376.5600)),
          (800, 1050, (81853.30, -157961.0, 113333.0, -28439.40, -8628.280, -42482.30, 106679.0, -376.5600)),
          (1050, 1221, (-219.1790, 607.2850, -381.9350, 85.26800, 51.00700, -238.6090, -198.2940, -376.5600))],
 "Cg":   [(298, 6000, (21.17510, -0.812428, 0.448537, -0.043256, -0.013103, 710.3470, 183.8734, 716.6690))],
}
S298 = {k: shom(v[0][2], 298.15)[1] for k, v in SH.items() if v[0][0] <= 298.2}
# ---- (b) enthalpy per mole of vapor that sets d ln P_total / d(1/T) ----------------------------------------
# JANAF FeS dfG (kJ/mol, Fe + 1/2 S2(g) = FeS): linear fit 900-1400 K (solid) gives dH, dS
FeS = {900: -103.437, 1000: -97.670, 1100: -91.766, 1200: -85.786, 1300: -79.802, 1400: -73.918, 1500: -68.943, 1600: -65.382}
Ts = np.array([900, 1000, 1100, 1200, 1300, 1400.]); Gs = np.array([FeS[int(x)] for x in Ts])
sl, ic = np.polyfit(Ts, Gs, 1)          # dfG = dH - T dS  ->  dH = ic, dS = -sl
dH_FeS_form = ic                        # kJ/mol (negative)
cand = [
 # name, dH (kJ per mol of vapor, controls slope of total P), basis
 ("H2O ice -> H2O(g) (reference)", (-241.83 + 285.83) + 6.01, "JANAF dfH(g) - dfH(l) at 298 K + fusion"),
 ("S liquid -> total vapor (CRC fit 375-717 K)", 9028.6679 * R / 1e3, "slope 9,029 K of the CRC fit used in Section S9"),
 ("S8: 8 S(rhombic) -> S8(g)", 100.416, "JANAF dfH298 of S8(g)"),
 ("S8: 8 S(liquid) -> S8(g)", 100.416 - 8 * 1.85, "JANAF dfH298 of S8(g) and S(l)"),
 ("S2: 2 S(rhombic) -> S2(g)", 128.60, "JANAF dfH298 of S2(g)"),
 ("S2: 2 S(liquid) -> S2(g)", 128.60 - 2 * 1.85, "JANAF dfH298 of S2(g) and S(l)"),
 ("NaCl(cr) -> NaCl(g)", -181.42 + 411.12, "JANAF dfH298"),
 ("KCl(cr) -> KCl(g)", -214.68 + 436.68, "JANAF dfH298"),
 ("Na2S(cr) -> 2 Na(g) + 1/2 S2(g)", (2 * 107.30 + 0.5 * 128.60 + 366.10) / 2.5, "JANAF dfH298; congruent, per mol of gas (2.5)"),
 ("K2S(cr) -> 2 K(g) + 1/2 S2(g)", (2 * 89.00 + 0.5 * 128.60 + 376.56) / 2.5, "JANAF dfH298; congruent, per mol of gas (2.5)"),
 ("MgS(cr) -> Mg(g) + 1/2 S2(g)", (147.10 + 0.5 * 128.60 + 345.72) / 1.5, "JANAF dfH298; congruent, per mol of gas (1.5)"),
 ("CaS(cr) -> Ca(g) + 1/2 S2(g)", (177.80 + 0.5 * 128.60 + 473.21) / 1.5, "JANAF dfH298; congruent, per mol of gas (1.5)"),
 ("FeS(cr) -> Fe(s) + 1/2 S2(g)", -2 * dH_FeS_form, "JANAF Fe-025 dfG at 900-1400 K; P_S2 ~ K^2, per mol S2"),
 ("C(graphite) -> C(g)", 716.669, "JANAF dfH298 of C(g)"),
]
rows = []
for nm, dH, basis in cand:
    lr_half = 0.5 * dH * 1e3 / R * dx / math.log(10)     # log10 of predicted depth ratio, front ~ P^(1/2)
    lr_one = dH * 1e3 / R * dx / math.log(10)
    lr_pk = 0.5 * dH * 1e3 / R * dxp / math.log(10)
    rows.append(dict(candidate=nm, dH_kJ_per_mol=round(dH, 1), basis=basis,
                     log10_ratio_219_324K_sqrtP=round(lr_half, 2), log10_ratio_219_324K_P=round(lr_one, 2),
                     log10_ratio_peakT_sqrtP=round(lr_pk, 2), dH_over_dHapp_meanT_sqrtP=round(dH / fit["mean_annual"]["dH_app_kJ_n0.5"][0], 1)))
# ---- (c) temperature at which P_total equals rho g h / phi ---------------------------------------------------
Pb, Tb = 101325.0, 717.8
def T_S(P, B): return 1 / (1 / Tb - math.log(P / Pb) / B)
def T_antoine(P, A, B, C): return B / (A - math.log10(P / 1e5)) - C
def T_cc_from_antoine(P, A, B, C, Tn):   # Clausius-Clapeyron with the Antoine slope at the normal boiling point Tn
    slope = math.log(10) * B * Tn ** 2 / (Tn + C) ** 2       # K
    return 1 / (1 / Tn - math.log(P / Pb) / slope)
def Ptot_congruent(T, cond, metal, nmetal, mol=None):
    """total pressure (Pa) over cond -> nmetal M(g) + 1/2 S2(g), congruent (P_M = 2 nmetal P_S2), plus molecular MS(g) if given"""
    dG = nmetal * Gabs(SH[metal], T) + 0.5 * Gabs(SH["S2g"], T) - Gabs(SH[cond], T)
    lnK = -dG / (R * T)                       # K = P_M^n P_S2^0.5 (bar)
    k = 2 * nmetal                            # P_M = k x, P_S2 = x
    x = math.exp((lnK - nmetal * math.log(k)) / (nmetal + 0.5))
    P = (k + 1) * x
    if mol:
        dGm = Gabs(SH[mol], T) - Gabs(SH[cond], T); P += math.exp(-dGm / (R * T))
    return P * 1e5
def T_for(Pfun, P):
    return brentq(lambda T: math.log(Pfun(T)) - math.log(P), 400, 7000)
def FeS_P(T):    # P_S2 over FeS + Fe, JANAF dfG; linear extrapolation in T beyond 1600 K from the 1500-1600 K (liquid) rows
    if T <= 1600: G = np.interp(T, sorted(FeS), [FeS[k] for k in sorted(FeS)])
    else: G = FeS[1600] + (FeS[1600] - FeS[1500]) / 100 * (T - 1600)
    return math.exp(2 * G * 1e3 / (R * T)) * 1e5
def C_P(T): return math.exp(-(Gabs(SH["Cg"], T) - (0 - T * 5.74)) / (R * T)) * 1e5 if False else math.exp(-(716669 - T * (158.1 - 5.74)) / (R * T)) * 1e5
Pmeth = {
 "H2O liquid (Antoine, Liu and Lindsay 1970)": [lambda P: T_antoine(P, 3.55959, 643.748, -198.043)],
 "Native sulfur (Section S13 curve)": [lambda P: T_S(P, 7750.0), lambda P: T_S(P, 8420.0)],
 "NaCl liquid, total vapor (Stull 1947; Antoine and CC extrapolated above 1,738 K)": [lambda P: T_antoine(P, 5.07184, 8388.497, -82.638), lambda P: T_cc_from_antoine(P, 5.07184, 8388.497, -82.638, 1738.0)],
 "KCl liquid, total vapor (Stull 1947; extrapolated above 1,680 K)": [lambda P: T_antoine(P, 4.78236, 7440.691, -122.709), lambda P: T_cc_from_antoine(P, 4.78236, 7440.691, -122.709, 1680.0)],
 "Na2S -> 2Na + 1/2 S2 (JANAF; solid extrapolated)": [lambda P: T_for(lambda T: Ptot_congruent(T, "Na2Ss", "Nag", 2), P)],
 "K2S -> 2K + 1/2 S2 (JANAF; solid extrapolated)": [lambda P: T_for(lambda T: Ptot_congruent(T, "K2Ss", "Kg", 2), P)],
 "MgS -> Mg + 1/2 S2, + MgS(g) (JANAF)": [lambda P: T_for(lambda T: Ptot_congruent(T, "MgSs", "Mgg", 1, "MgSg"), P)],
 "CaS -> Ca + 1/2 S2, + CaS(g) (JANAF)": [lambda P: T_for(lambda T: Ptot_congruent(T, "CaSs", "Cag", 1, "CaSg"), P)],
 "FeS -> Fe + 1/2 S2 (JANAF Fe-025; extrapolated above 1,600 K)": [lambda P: T_for(FeS_P, P)],
 "Graphite -> C(g) (JANAF, monatomic only; 298 K dH, dS)": [lambda P: T_for(C_P, P)],
}
tb = {}
for nm, fs in Pmeth.items():
    tb[nm] = {}
    for phi in (1.0,) + phis:
        vals = [f(r * g * h / phi) for f in fs for r in rhos for h in hs]
        tb[nm][str(phi)] = [round(min(vals)), round(max(vals))]
    tb[nm]["T_at_1bar"] = [round(min(f(1e5) for f in fs)), round(max(f(1e5) for f in fs))]
# check the 298 K-only estimates for the congruent sulfides (dCp = 0), as in analysis_v3/sulfide_s2.py
def T298(dH, dS, nmetal, P):
    k = 2 * nmetal; x = P / 1e5 / (k + 1); lnK = nmetal * math.log(k * x) + 0.5 * math.log(x)
    return dH * 1e3 / (dS - R * lnK)
chk = {}
for nm, dH, dS, n in (("Na2S", 2 * 107.30 + 64.30 + 366.10, 2 * 153.65 + 0.5 * 228.19 - S298["Na2Ss"], 2),
                      ("K2S", 2 * 89.00 + 64.30 + 376.56, 2 * 160.34 + 0.5 * 228.19 - S298["K2Ss"], 2),
                      ("MgS", 147.10 + 64.30 + 345.72, 148.65 + 0.5 * 228.19 - 50.30, 1),
                      ("CaS", 177.80 + 64.30 + 473.21, 154.89 + 0.5 * 228.19 - 56.54, 1)):
    chk[nm] = {str(phi): [round(min(T298(dH, dS, n, r * g * h / phi) for r in rhos for h in hs)),
                          round(max(T298(dH, dS, n, r * g * h / phi) for r in rhos for h in hs))] for phi in phis}
rec = {"pressures_MPa": {str(phi): [round(rhos[0] * g * hs[0] / phi / 1e6, 3), round(rhos[1] * g * hs[1] / phi / 1e6, 3)] for phi in (1.0,) + phis},
       "dHapp_kJ": {"meanT_n1": fit["mean_annual"]["dH_app_kJ_n1"], "meanT_n0.5": fit["mean_annual"]["dH_app_kJ_n0.5"],
                    "peakT_n1": fit["peak_noon"]["dH_app_kJ_n1"], "peakT_n0.5": fit["peak_noon"]["dH_app_kJ_n0.5"]},
       "candidates_dH": rows, "T_base_K": tb, "T_base_298K_estimates_sulfides": chk,
       "S298_from_Shomate": {k: round(v, 2) for k, v in S298.items()},
       "FeS_formation_fit_900_1400K": {"dH_kJ": round(ic, 2), "dS_J_K": round(-sl * 1e3, 2)},
       "sources": SOURCES}
json.dump(rec, open(os.path.join(HERE, "volatile_candidates_record.json"), "w"), indent=1)
import csv
with open(os.path.join(HERE, "volatile_candidates.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow(x) for x in rows]
print(json.dumps({k: rec[k] for k in ("pressures_MPa", "dHapp_kJ", "T_base_K", "T_base_298K_estimates_sulfides", "S298_from_Shomate", "FeS_formation_fit_900_1400K")}, indent=1))
for x in rows: print(x)
