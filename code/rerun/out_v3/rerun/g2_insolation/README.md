# g2_insolation: insolation and plains results rerun on the v3 and v3b catalogs (6 Nov 2026)

This folder reruns the analyses behind main-text Section 3.5, the plains part of Section 3.2, Supplement Sections S3, S5, S6
and S8 (Tables S5 and S6), Figures 1b and 3, and the matching statements in Table 1 and the Discussion of
`Hollow_Formation_Model_Paper_v6_3.tex` / `supplementary-materials_v13.tex`, after the spurious detections were removed
(`code/rerun/out_v3/README.md`). No existing file outside this folder was modified.

## Files
- `changes_g2.csv`: one row per catalog-dependent number (94 rows): `text_value` as printed in the paper, `v2_record` from the
  stored record of the original script, then `v3` and `v3b`. `conclusion_change` is "yes" when, relative to v2, a ratio or
  odds ratio crosses 1 or its 95% interval newly includes or excludes 1, a p value crosses 0.05, or an AIC ranking flips. The
  bracket names the catalog(s). Rows marked "plotted" are figure values that the text does not quote.
- `prep_inputs.py`: builds `inputs/<v2|v3|v3b>/` (hollows, fields, host assignments, `test2_craters_h18.csv` and
  `crater_features.csv` with hosts, h09 and n09 recomputed from the host assignments, and stage labels).
- `make_copies.py`: copies the 21 original scripts here and applies only path and env replacements; each replacement is listed
  in the script and asserted to match once. `make_copies.log` maps each copy to its original. `fractions.py` is saved as
  `fractions_g2.py`, because a file named `fractions.py` here shadows the Python standard-library module.
- `run_g2.sh`: runs the copies for one catalog (`G2_CAT=v3 bash run_g2.sh`); outputs go to `out/<cat>/`, logs to
  `out/<cat>/logs/`, and exit codes to `out/<cat>/run_status.txt` (all 0).
- New scripts: `wall_fractions.py` (the Table S6 wall fractions, computed with the code and seed of panel c of
  `figscripts/make_figs_pr_insolation_v9.py`, which draws but does not save them), `plains_hosts.py` (the S6 counts of plains
  hosts), `check_v2_repro.py`, and `compile_changes.py` (writes `changes_g2.csv`; it also derives the Table S5
  likelihood-ratio p values from the AICs).
- `inputs/coverage_10ppd_v2.npz` is a relative link to the v2 coverage grid left from setup (files cannot be deleted here); the
  scripts do not use it.
- `out/<cat>/fig3_insolation_v9.png` and `figS5_insolation_v9.png`: Figure 3 and Figure S5 redrawn per catalog. The axis
  limits are those of the original script, so in the v3b figure one point of panel a (2.8 at 315E) is clipped.

## Inputs and one rule change for fields
- Host craters: `host_assignments_v3.csv` / `_v3b.csv`. Any-hollow hosts within 60 deg drop from 657 to 485 (v3) and 392
  (v3b). h09 hosts in the containment-graph table drop from 606 to 446 and 370 (the README's 448 and 371 include two hosts that
  the table drops as possible peak rings).
- Fields: field ids and centroids are kept from v2, and a field is kept only if it has at least one non-trace hollow. Fields
  drop from 1,228 to 935 (v3) and 757 (v3b). Within 15 deg of the equator (Caloris excluded), fields drop from 379 to 217 and 128.
  The v2 field analyses counted all 1,228 fields, including 33 fields made only of the 64 v2 crater traces; the v2 column
  keeps that rule, so the v3 change includes dropping those 33 fields.
- For the 60-75 deg incidence densities (`latcorr5.py`), a field qualifies in v3/v3b only through a non-trace hollow.

## v2 reproduction
Run on the v2 inputs, the copies reproduce every stored record exactly (`out/v2_reproduction.json`), with two exceptions that do
not affect any quoted value: the completeness-controlled fit (AIC 4,658.03 against 4,657.96; same values at the precision
printed), because the crater-graph completeness table it reads has changed slightly since the record was written; and the
circular-median intervals of `aspect_record.json`, which the paper does not quote.

## What changes
1. **Size dependence of the insolation effect.** Removing spurious detections removes most hosts 5-10 km across (130 to 46 and
   42 within 60 deg), and the Delta T_w x log D interaction is no longer significant: p = 1.5e-5 becomes 0.061 (v3) and 0.10
   (v3b); with completeness, 0.0003 becomes 0.18 and 0.23. The odds ratio per 10 K is 1.38 (1.21-1.57) for 10 km craters and
   1.19 (1.09-1.30) for 100 km craters in v3 (v2: 1.60 and 1.18). The interaction is not significant at 15-30N or 30-45N
   (p = 0.88 and 0.18 in v3), so the statements that the rise weakens with crater size and that this comes from 15-45N
   (Section 3.5, S3, S8, Figure 3b, Discussion) no longer hold.
2. **Hot longitudes to 45N.** The crater-level odds ratio stays above 1 within 15 deg (2.19, 1.48-3.24 in v3; 2.06,
   1.30-3.27 in v3b) and at 15-30N (2.64, 1.55-4.48; 2.43, 1.38-4.27), but at 30-45N it becomes 1.49 (0.93-2.40) and 1.42
   (0.88-2.31). "Near the equator and to 45N" (Section 3.5, Table 1, Discussion) becomes "near the equator and to 30N".
3. **Field-level hot-longitude excess.** Within 15 deg it falls from 1.95 (1.59-2.47) to 1.63 (1.25-2.13) in v3 and 1.29
   (0.91-1.79) in v3b; the cos 2 lon amplitude falls from 0.67 to 0.49 (0.33-0.65) and 0.39 (0.19-0.59). In v3b the sector
   ratio and the crater-population comparisons (1.71/1.84/2.08 become 1.14/1.23/1.41) no longer exclude 1. At 15-30N the
   field ratio is 1.20 (0.87-1.70) in v3. The reversal at higher latitudes strengthens: 30-45N becomes 0.69 (0.49-0.97,
   p = 0.038), and 45-60N without the province becomes 0.59 (0.37-0.90, p = 0.021), so "the remaining reversal is marginal"
   no longer holds.
4. **Shape of the rise with Delta T_w (S8).** The linear odds ratio holds (1.37, 1.25-1.50 in v3), but the best single step,
   now at 34 K, fits slightly better than the linear term (AIC lower by 1.6 in v3 and 3.0 in v3b, against 15.6 higher in v2),
   and with both terms the step is significant (p = 0.004 and 0.007) while the linear term weakens (p = 0.035 and 0.12). The
   statement that a step at any single value fits worse than a steady rise is reversed.
5. **368 K layer against Delta T_w (S5, Table S5, Discussion).** Delta T_w still has the lower AIC in every latitude range,
   but the contrast in the text reverses: within 60 deg Delta T_w is now better by 8.5 (v3) or 7.2 (v3b) instead of "about
   equally", and across all latitudes the difference shrinks from 14.3 to 5.3 and 2.8 (Discussion: "AIC 5,790 against 5,804"
   becomes 4,423 against 4,428 in v3). Adding the layer to Delta T_w now improves the all-latitude fit (p = 0.017 and 0.022,
   against 0.58).
6. **Plains (3.2, S6).** The plains scarcity holds: density inside / outside becomes 0.35 (v3) and 0.32 (v3b) at 40-60N
   (v2: 0.41) and stays 0.13 at 60-80N; the hosting odds ratio becomes 0.40 (0.23-0.70) and 0.38 (0.21-0.67); 19 of 20 (v3)
   and 18 of 19 (v3b) plains hosts are >= 10 km. Field density outside the plains at 0-30N falls from 61 to 52 and 44, so it
   still does not increase toward the equator.

Unchanged within rounding: the Stage 1 step at 44 K (3.05, 1.50-6.20 in v3; AIC lower by 27.5), the Stage 1 fractions
(56-71% and 96-100%), the equator-facing wall fractions at 60-75N (0.75 per hollow, 0.71 per crater), the 45-60N field
reversal (0.39 in v3), the hollow-rich-host counts (8 of the 10 richest hosts near the warm longitudes, holding 50% of the
hollows), and the cos 2 lon comparison (AIC higher by 20).

## Not rerun
- Equator-facing / pole-facing wall density 1.52 (0.95-2.47): from the terrain analysis of S4 / Table S2.
- The reference-size checks of S2 for the hot-longitude odds ratios (2.2-2.6, 3.2-3.7, 1.76-1.98): completeness-group
  scripts. Note that S2 says the 30-45N excess loses significance "only with completeness for 2 or 5 km hollows", which v3
  contradicts (it already does so at 1 km).
- "The odds of a Stage 2 or 3 classification double per 10 deg of solar incidence" (stage_vs_incidence.json): not
  rerun; it depends on the catalog only through the 2,468 labeled hollows, which become 2,405 and 2,372.
- Figure 1 was not redrawn, because panels c and d use other groups' records. Its panel b values are in `changes_g2.csv`.
- The DEM aspect check in `code/rerun/latitude/aspect_robustness/` needs a DEM path that does not exist on this machine. The
  only number of that script quoted in the paper (15-30N with Caloris excluded) is reproduced by `wall_fractions.py`.
- The joint-model share of insolation (25%) quoted in the Discussion belongs to the joint-model group.

- stage_by_band.py (8 Nov 2026; G2_CAT=v3 -> out/v3/stage_by_band.json): the Stage 1 vs Stage 2 hot-sector test of stage_test.py by 15-degree latitude band, labels from images at incidence < 70 deg and all labels. Not in the paper. With < 70 deg labels, Stage 1 share hot/warm: |lat|<15 0.99/0.85, 15-30N 0.97/0.66, 30-45N 0.97/0.63, 45-60N 0.82/0.53, 60-75N 0.73/0.38 (n = 19); controlled OR 30-60N 2.16 (0.99-4.71), 30-75N 2.35 (1.13-4.89); at low latitudes the hot-sector labels come from much lower incidence (28 vs 52 deg), so the controlled OR is uninformative there.
