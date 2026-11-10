# analysis_v9: crater containment graph with a peak-ring basin relation (1 November 2026)

**Change from v8.** v8 applied the central-peak uplift relations of Cintala and Grieve (1998) to every complex crater.
Extrapolated to Caloris they gave 167 km (terrestrial) and 935 km (lunar) of uplift. v9 sets the uplift beneath features
with peak-ring geometry to u_s = k_ring x 0.1 D_tc with k_ring = 1.05 (79 km for Caloris). Peak-ring geometry means an
Orgel et al. (2020) basin, a Herrick et al. (2018) crater coded pr or pb, or a crater > 120 km not coded cp or mp.
Sources: Baker and Head (2015), Kring et al. (2016), Urrutia-Fucugauchi et al. (2019); Rumil tag mercury-provenance.

**Files.** crater_dag.py (v9 patch; CRATER_DAG_BASIN=v8 reproduces v8 exactly, checked in check_v8_rule/),
crater_dag_tests.py and crater_dag_sens.py (unchanged; records identical to v8), uplift_lunar/, kring_0p8/,
kring_1p3/, compare_v8_v9.py -> compare_v8_v9_record.json. Each step runs in a few seconds in the foreground
(run_v9.sh does not survive as a background job on this machine).

**Results used in Hollow_Formation_Model_Paper_v4 and supplementary-materials_v9 (Section S14, Table S16).**
- Nesting, inheritance, position, ancestors: unchanged from v8 (nested 5,041.2 vs base 5,102.0).
- z_base lowers the AIC by 5 (p = 0.007) with either uplift relation (v8 lunar: p = 0.51).
- ln z_max added to ln z_own lowers the AIC by 11 (p = 0.0002) with either relation.
- Best depth-of-origin model minus nested: 55.9 (terrestrial), 56.0 (lunar); 55.7 and 56.2 for k_ring = 0.8 and 1.3.
- Hosts with z_max > 35 km (mean crust, Padovan et al. 2015): 84 of 705 (85 lunar); > 53 km: 46 (47 lunar).

**Moved out of the paper.** The full terrestrial/lunar comparison, the fill-known subset, and the k_ring runs are kept
here for the provenance (depth-of-origin) paper; see ../cut_material_v4.tex.
