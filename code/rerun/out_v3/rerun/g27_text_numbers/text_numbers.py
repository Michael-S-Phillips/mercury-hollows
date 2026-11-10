#!/usr/bin/env python3
"""g27 (9 Nov 2026): numbers that the paper computes from stated inputs, recorded so the archive can point to them.
(1) S13 "Steady conduction above the hot rock": diffusion time H^2/kappa for H = 60-280 m; gradients (T_base - T_s)/h for the
    base temperatures and depths of Section 4.2 (855 K at 65 m for phi = 0.4; 1,000 K at 260 m for phi = 0.1; T_s = 300 K); the
    exponent of base depth on gradient; the range of H allowed by a twofold range of depth; heat flux k * gradient for
    k = 1-2 W/m/K; and the rise in mean surface temperature q / (4 eps sigma T_s^3) for eps = 0.95.
(2) S8 "Noon temperature of walls": radiative equilibrium at noon at perihelion (0.3075 AU), albedo 0.09, emissivity 0.95, as drawn
    in Figure S6b (g8_figS5/make_figS5_v10.py): a 30 deg equator-facing wall at 70 N (incidence 40 deg) against flat ground at 70 N.
(3) S8 "Wall orientation": share of the wall area of the 111 hand-traced host craters with 665 m regional slope < 5 deg, from the
    traced terrain record (g18_traced_r3/terrain/out_v3/peak_wall_record_rows_tracedonly_r3_v3.json, rows_tracedonly_r3_v3, all,
    slope_stratified; s0 = slope < 5 deg, s1 = 5-15 deg, s2 >= 15 deg in peak_wall.py). 10 Nov 2026: was the automated relief-threshold
    classes (g3_catalog/terrain/out_v3/peak_wall_record.json, rows_H250_r0.45; 69.1% over 321 hosts); the paper uses traced classes only.
(4) Main text: conduction times L^2/kappa (Section 4.1.1); excavation depth of a 5 km crater, 0.1 D_t with D_t = D / 1.2 for simple
    craters (Section S7 scaling); pixels across a 0.1 km hollow at 6 m/pixel (Section 4.4).
(5) S2: 1st-99th percentiles of hollow diameter in pixels and of solar incidence over the 85,602 hollow-image pairs of the v3 catalog
    (g3_catalog/inputs/completeness_pairs_v3.csv), the pairs to which the detection model was refit.
Writes text_numbers.json."""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RR = os.path.dirname(HERE)
yr = 3.15576e7; kappa = 1e-6; sig = 5.670e-8; R = {}
# (1)
Ts = 300.0; cases = {"phi0.4": (855.0, 65.0), "phi0.1": (1000.0, 260.0)}
grad = {k: (Tb - Ts) / h for k, (Tb, h) in cases.items()}
expo = np.log(cases["phi0.1"][1] / cases["phi0.4"][1]) / np.log(grad["phi0.1"] / grad["phi0.4"])
q = [1.0 * min(grad.values()), 2.0 * max(grad.values())]
R["S13_steady_conduction"] = dict(time_yr_H60_H280=[60.0 ** 2 / kappa / yr, 280.0 ** 2 / kappa / yr], gradient_K_per_m=grad,
                                  exponent_depth_on_gradient=expo, H_range_for_twofold_depth=2 ** (1 / abs(expo)),
                                  heat_flux_W_m2=q, surface_warming_K=[x / (4 * 0.95 * sig * Ts ** 3) for x in q])
# (2)
A, eps, Sp = 0.09, 0.95, 1361 / 0.3075 ** 2
Tn = lambda inc: ((1 - A) * Sp * np.cos(np.radians(inc)) / (eps * sig)) ** .25
R["S8_noon_wall"] = dict(T_wall_70N_K=Tn(40.0), T_flat_70N_K=Tn(70.0), excess_K=Tn(40.0) - Tn(70.0))
# (3)
s = json.load(open(f"{RR}/g18_traced_r3/terrain/out_v3/peak_wall_record_rows_tracedonly_r3_v3.json"))["rows_tracedonly_r3_v3"]["all"]["slope_stratified"]
aw = [s[k]["A_wall"] for k in ("s0", "s1", "s2")]
R["S8_wall_area_slope_lt5"] = dict(classes="hand-traced, 111 hosts", A_wall_km2_lt5_5to15_ge15=aw, fraction_lt5=aw[0] / sum(aw))
# (4)
R["main_conduction_time"] = {f"L{L}m": dict(s=L ** 2 / kappa, yr=L ** 2 / kappa / yr) for L in (200.0, 5000.0)}
R["main_excavation_5km_km"] = 0.1 * 5.0 / 1.2
R["main_pixels_0.1km_at_6m"] = 100.0 / 6.0
# (5)
p = pd.read_csv(f"{RR}/g3_catalog/inputs/completeness_pairs_v3.csv")
R["S2_pair_ranges_v3"] = dict(n_pairs=int(len(p)), n_detected=int(p.detected.sum()), n_hollows=int(p.hollow_id.nunique()),
                              D_px_p1_p99=[float(x) for x in np.percentile(p.D_px, [1, 99])], inc_p1_p99=[float(x) for x in np.percentile(p.inc, [1, 99])])
json.dump(R, open(f"{HERE}/text_numbers.json", "w"), indent=1)
print(json.dumps(R, indent=1))
