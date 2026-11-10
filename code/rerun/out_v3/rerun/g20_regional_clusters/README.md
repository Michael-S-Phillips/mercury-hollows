# g20_regional_clusters: where hosting is clustered beyond the joint model (9 Nov 2026)

Michael asked to see the regional clustering behind the joint-model result that the regional hosting fraction raises R^2 by 0.03.
Not in the paper.

- `clusters.py`: predicted hosting probability P for each of the 31,554 craters from the paper's joint model (thinned presence,
  complementary log-log, seven physical groups; coefficients from g4_completeness/g1/summary_model/v3r/out/improved_coefs.csv,
  'thinned_presence_all'; sum P = 438.8 against 446 hosts). On a 2-degree grid, observed hosts O and expected hosts E within 500 km,
  z = (O - E) / sqrt(sum P(1 - P)). Family-wise thresholds from 1,000 simulated data sets: z >= 5.10 (excess), z <= -3.22 (deficit),
  over nodes with E >= 2. Writes grid.npz, clusters.json, craters_P.csv.
- `fig_clusters.py` -> fig_regional_clusters_v2.png (v1 had a clipped caption line).

Result: three regions of excess, no significant deficit (lowest z = -3.19).
1. Northern province, 45-75N, 95-145E: 23 hosts within 500 km of 65N 125E against 4.9 expected (z = 9.4); 55 against 21.5 over the
   whole region. Many hosts there are small craters the model gives low P (e.g. 25 km at 71.7N 104.6E with 14 hollows, P = 0.03).
2. 5-25N, 325-340E (Abu Nuwas, Asvaghosa, Sinan, Yeats, Ts'ai Wen-chi, Geddes): 23 against 8.0 near 13N 329E (z = 6.1).
3. Around Tyagaraja, 0-10N, 195-210E: 19 against 6.2 near 7N 203E (z = 6.0); Tyagaraja itself is expected (P = 0.999), the excess is
   its neighbors, including the 112 km ghost crater at 17.8N 207.3E.
Hosts within 500 km of a significant node: 118 of 446, against 51.5 expected.

## Radius check and paper figure (9 Nov 2026)
- clusters.py takes G20_RAD (km; default 500). G20_RAD=300 -> clusters_r300.json (threshold 5.08): the province (as two parts, 51-65N and 39-45N near Pasch), 9-19N 327-335E, and -1-3N 199-201E. G20_RAD=750 -> clusters_r750.json (threshold 4.83): the province, 7-37N 311-341E (toward Praxiteles and Velazquez), 1N 209-211E, and a deficit near 19-21N 35-39E (8 vs 22.8; threshold -3.36).
- fig_regional_paper.py -> fig6_regional_excess_v1.png/.pdf, copied to manuscripts/.../figures (main-text Figure 6 of v6.9).

## What sets the three regions apart (9 Nov 2026; not in the paper)
- region_profiles.py -> region_profiles.json: craters within 500 km of each peak against all craters.
  - Province (65N 125E): the excess is in small craters (5-10 km: 6 hosts vs 0.8; 10-30 km: 9 vs 0.9); hollows small (median 0.51 km, 49% < 0.5 km) and detected at ~15 m/pixel; LRM, S/Si, Mg/Si, Ca/Si below average (-0.6 to -1.5 SD); 3 vents (3.2 expected by area); 1,776 km from the Caloris center.
  - 13N 329E: excess in craters >= 30 km (15 vs 5.1); hollows large (median 1.6 km); S/Si, Mg/Si, Ca/Si, Al/Si above average (+0.5 to +0.65 SD; hosts S/Si +0.8); bluer, low PC2; no vents within 500 km; more standard craters (27% of craters >= 10 km vs 13%).
  - Tyagaraja (7N 203E): excess in craters >= 10 km; hollows large (median 1.7 km); LRM -0.45 SD, low PC2; 3 vents (as expected by area); 1,925 km from the Caloris center.
- region_checks.py -> region_checks.json: 17 of the province's 23 hosts owe their status only to hollows < 1 km (10 only to < 0.5 km), at a median 15 m/pixel; none in the other two regions. Observed/predicted hosts by distance from the Caloris center is 2.14 at 1,600-2,000 km (49 vs 22.9) and 0.67-1.11 in the other rings; in the 1,500-2,200 km ring the excess lies in two opposite sectors (azimuth 315-360: 27 vs 6.5; 90-135: 17 vs 5.6).
- region_tests.py -> region_tests.json: refits with completeness for 0.3 or 0.5 km hollows (or 0.3 and 1 km) leave every region in excess (province z 9.4-11.2). Within 30 deg of the equator, observed/predicted hosts shows no rise at high Delta T_w (0.76-1.16 above 45 K), so a threshold response to insolation does not explain regions 2 and 3.
- Not done: the province's sub-km detections were never screened (the v3 screen covered detections >= 1 km), so false positives there remain untested.
