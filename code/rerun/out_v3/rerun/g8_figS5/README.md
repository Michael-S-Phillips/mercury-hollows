# g8_figS5: Figure S5 v10 (6 Nov 2026)

`make_figS5_v10.py` (a copy of `g5_figures/make_figs_pr_insolation_v9_v3cat.py`) writes `figS5_insolation_v10_v3cat.png/.pdf`.
- Old panels a (field ratio against longitude) and b (odds ratio per 10 K against crater diameter) are dropped, because
  they are the same panels as main-text Figure 4a-b. The script still draws them, on a figure it never saves.
- Old panels c (equator-facing wall fraction) and d (noon temperatures) become a and b, unchanged.
- Old panel e becomes c: hosting against Delta T_w for the v3 catalog with the refit completeness grid (v3r inputs from
  `g4_completeness/g2/inputs/v3r`), without the single-step model. The fit is a logistic regression with Delta T_w
  entering linearly in the log odds (odds x1.37, 1.25-1.50, per 10 K), so it curves on a linear probability axis; the
  probability axis is now logarithmic, on which the fit is close to straight (at these probabilities log odds and log
  probability differ by less than 0.5%).
Record: `figS5_panel_c_v10_record.json`.
