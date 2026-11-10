# Host-crater degradation with the Kinczyk et al. (2020) classes (5 Nov 2026)

Workshop analysis only; the manuscript (v6.2) and supplement (v12) are unchanged.

## Data
- data/crater_degradation_database_MJK_2019.txt: the supplement of Kinczyk et al. (2020, Icarus 341, 113637), from Michael's Downloads.
  - Columns: FID, x_coord (longitude, -180..180), y_coord (latitude), Diameter_km, Class.
  - Contents: 3,253 craters >= 40 km from Fassett et al. (2011).
  - Class 1 is the most degraded and Class 5 the freshest (rayed; Kuiperian). Class 4 is Mansurian (no rays).
  - Classes 1-5 have 1,048 / 1,300 / 521 / 219 / 22 craters, matching their Table 2, N(40).
  - Checked against the published paper (Downloads/1-s2.0-S0019103519303598-main.pdf). For classes 1-5, the file reproduces every N(D) column of their Table 2 (40-100, 40-200, >=40, >=100, >=200, >=400 km).
  - Code 6 (140 craters) is their "ghost" class (Section 2.2: "craters entirely buried by subsequent lava flows"; Section 3.1: craters embayed by the smooth plains were "identified as ghost craters" rather than given C1 or C2). The paper doesn't number it, and it is not in Table 2.
    - Their median latitude is 64.5 N.
    - Of the 80 matched, all are Herrick "filled" and 56 (70%) lie on the northern smooth plains, against 2-6% for classes 1-5.
  - Code 0 (3 craters) is not described in the paper and is not counted in Table 2, so it is most likely unclassified entries (0 being the default of an empty integer field).
    - Two lie at the 40 km threshold (40.0 and 41.0 km). The third is 118 km, at 179.5 E on the longitude seam.
    - Herrick rates two as standard and one as filled; none hosts hollows.
  - Both codes are left out of the class fits; including them would not change the results.
- data/craterstats_files.zip: their crater-size-frequency files by class (not used).

## Scripts and records
- match_kinczyk.py → kinczyk_matched.csv, match_summary.json.
  - Each Kinczyk crater is matched one-to-one to a crater in summary_model/crater_features.csv. The centers must be within one Kinczyk radius and the diameter ratio within 0.67-1.5.
  - 3,051 of 3,253 matched. Median offset is 0.12 radius (95th percentile 0.28), and the diameter ratio is 0.90-1.12 (5th-95th percentile).
  - Unmatched: 113 Class 1, 60 code 6 and 29 other craters.
  - Restricting to offsets <= 0.3 radius leaves the results unchanged (match_sensitivity.json).
- kinczyk_tests.py → kinczyk_tests.json: hosting (h09) by class, the head-to-head with Herrick, Delta T_w by class, and depth.
- kinczyk_depth_check.py → kinczyk_depth_check.json: depth by Kinczyk and Herrick class on the same hollows, clustered by image or host, with and without Herrick-"filled" hosts.

## Results
- **Coverage.** Kinczyk classifies 2,968 of the matched craters >= 40 km (258 hosts).
  - Herrick gives a preservation state to only 908 craters >= 40 km (160 hosts), because he leaves filled craters unrated (2,150 of the matched ones).
  - Herrick also covers 10-40 km, which Kinczyk does not: 4,214 craters >= 10 km with a state, 303 hosts.
- **Hosting by class** (logit with log D and |lat|, as in S11):
  - Hosting rate by class: 3.9%, 5.5%, 16%, 27% and 41% for classes 1-5.
  - Odds ratios against Class 1: 2.2 (1.4-3.4), 8.7 (5.5-13.9), 18 (11-31) and 37 (14-98) for classes 2-5.
  - The odds rise 2.8 (2.4-3.2) times per class step.
  - With completeness (C1) added, the odds ratios are 2.0, 8.1, 12 and 18, and the per-step factor is 2.45 (2.1-2.9).
  - With the full control set (adding the northern smooth plains, Delta T_w and the local hosting rate), the per-step factor is 2.8 (2.4-3.3).
- **Herrick in the same data (h09).**
  - Craters >= 10 km: fresh vs degraded 5.2 (3.6-7.6) and standard vs degraded 1.8 (1.4-2.4). With C1, fresh vs degraded is 3.3.
  - Craters >= 40 km: fresh vs degraded 6.2 (3.5-10.8).
- **Head to head** (2,966 craters >= 40 km with both classifications; Herrick as filled/degraded/standard/fresh):
  - Kinczyk fits better in every specification: AIC 1,386.9 vs 1,390.6 with log D and |lat|, 1,229.2 vs 1,242.1 with C1, and 1,084.9 vs 1,093.1 with the full controls.
  - Each classification still adds information to the other (likelihood-ratio p = 1e-5 to 0.009).
  - Within Herrick-filled craters, Kinczyk's class separates hosting: 36 of 933, 51 of 1,035 and 13 of 97 craters host in classes 1-3, and the odds rise 2.0 (1.4-2.9) times per class.
- **Insolation by class** (|lat| <= 60):
  - Odds per 10 K of Delta T_w: 1.21 (1.04-1.42) in classes 1-2, 1.33 (1.08-1.64) in Class 3 and 1.24 (0.97-1.58) in classes 4-5.
  - The interaction p values are 0.66 and 0.92, so the dependence on Delta T_w does not differ by class (as with Herrick).
- **Depth** (1,484 hollows in 80 hosts >= 40 km; ln depth with ln pixel scale, ln tan i and |lat|):
  - Hollows in Class 3 hosts are 1.32 (1.14-1.51) times as deep as hollows in Class 1-2 hosts, and hollows in Class 4-5 hosts 1.27 (1.12-1.44) times. Clustering by host gives 1.32 (1.13-1.53) and 1.27 (1.05-1.54).
  - Without Herrick-filled hosts (41 hosts), the ratios are 1.25 (1.04-1.51) and 1.23 (1.00-1.51).
  - Herrick on the same hollows: fresh vs degraded is about 1.06, consistent with S11's 1.02.
  - So the depth signal is in the oldest hosts, mostly filled craters, which Herrick does not grade. Only 2 Class 5 hosts have depths.
- **Ghost craters** (code 6): 2 of the 80 matched craters host hollows (2.5%), about the rate of Class 1.

## Craters superposed on ghost craters (5 Nov 2026; Michael's idea)
Question: do craters that dig through a ghost crater's volcanic fill host hollows, like the Caloris craters that excavate below the plains (0/144 vs 14/35)?
- **Method** (ghost_superposed.py → ghost_superposed.json/.csv; ghost_expected.py → ghost_expected.json):
  - Fill over the ghost floor, lower bound: the fresh rim-to-floor depth 0.241 D^0.486, capped at 3 km for ghost basins. Variants add 0.5 and 1 km of plains above the buried rim.
  - Reach: excavation depth 0.1 D_tc, or the containment graph's z_own (which includes uplift).
- **Sample.** 582 craters >= 5 km are superposed on 83 ghost craters, with 7 hosts; 216 of them are in ghost interiors (rho < 0.6), with 3 hosts.
  - Fills are 1.8-3 km thick, so only craters of roughly 20-45 km or more reach below them: 17 in ghost interiors.
- **Result.** In ghost interiors, 1 of 17 piercing craters and 2 of 199 non-piercing craters host hollows (Fisher p = 0.22). At fixed size the odds ratio is 1.6 (0.25-9.8).
  - A hosting model fitted to craters off the plains, with size, latitude and completeness, expects 0.8 hosts among the 17 piercing craters. The ghost craters lie mostly at high northern latitude, where hosting is rare.
  - So the test cannot tell the two outcomes apart: observed counts match expected in both groups.
  - Craters superposed on ghost craters host at 0.86 (0.35-1.77) times the expected rate. Other northern-smooth-plains craters host at 0.47 (0.28-0.75) times.
- **Hollows inside ghost craters** (ghost_hollows.py): 30 of 5,409 hollows, in 6 ghost craters.
  - 18 of the 30 are in one 109 km ghost crater at 17.8 N, 207.4 E (Herrick 200653, filled, 0.42 km deep on the DEM, so the fill is about 2.0 km).
  - The host craters there are 6-12 km across and excavate 0.5-1.2 km, short of that fill. These hollows therefore come from the fill, or the fill is thinner than estimated. This crater is worth a look in the images.
- **General version** (fill_only_check.json): all filled parents with a DEM depth, outside Caloris (3,236 children on parent peaks or floors).
  - Children whose excavation stays within the fill host hollows at 0.37 (0.22-0.64) times the odds of children that dig through it, at fixed size, latitude and completeness. Adding the northern smooth plains gives 0.39 (0.22-0.67).
  - Raw rates: 24 of 1,891 (1.3%) vs 61 of 1,345 (4.5%).
  - Within the same parent (conditional logit; 23 parents, 31 hosts): 0.61 (0.12-3.2).
