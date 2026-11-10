# g26: measurement uncertainty of the shadow depths (9 Nov 2026)
Michael: is the uncertainty of the automatic shadow depths, relative to his hand depths, propagated into the results? And his
disagreement with Blewett et al. measures person-to-person variability. me.py -> me.json.
- Already in the intervals: random error of each depth enters the regression residuals; SEs clustered by image and image bootstraps
  carry it, including errors shared within an image. Unbiased random error in depth does not bias the coefficients.
- (B) Three-cornered hat on the 87 boxes measured three ways (Blewett max, Michael max, pilot automatic), ln depth: per-measurement error
  SD Michael 0.16 (robust 0.12; 95% 0-0.25), Blewett 0.28 (0.23; 0.20-0.35), automatic 0.29 (0.20; 0.19-0.39). The automatic method
  scatters about as much as a second person. Systematic person offset: Blewett/Michael 0.83 (g24), not in any interval; affects absolute
  depths only (reference depth, base temperature, layer thickness), not ratios within one set; no latitude dependence.
- (A) Differential error vs hand (bench, 111 hollows, 54 images): ln(hand/auto) does not depend on ln D (-0.07 +- 0.10), ln px, ln tan i,
  poleward (-0.01 +- 0.08) or multi (without uncertain or edge-flagged, 102 hollows: -0.05 +- 0.10 and -0.03 +- 0.08). Size exponent on the same hollows: hand 0.60 +- 0.10, auto 0.67 +- 0.07.
  On the Blewett boxes (coarser images, pilot run): ln(auto/hand) -0.27 +- 0.09 per ln px (absorbed by the px control) and +0.18 +- 0.09
  per ln D (could inflate the automatic size exponent by up to ~0.2); poleward 0.01 +- 0.07.
- (C) Error in hollow diameter (a regressor; not in any interval before): repeat detections at <= 20 m/px scatter by 0.30 in ln D
  (393 hollows, 941 detections) against 0.54 in the sample at fixed px and incidence: reliability 0.68. Method-of-moments correction
  (full / half error variance): size exponent 0.61 -> 0.90 (0.82-1.00) / 0.73 (0.65-0.80); size-controlled poleward ratio 0.98 -> 1.01 /
  0.99; Kinczyk 4-5 vs 1-2 0.97 -> 0.89 (0.80-0.97) / 0.93 (0.85-1.01); 3 vs 1-2 0.99 -> 0.95 / 0.98; nested 1.00 -> 1.03 / 1.02;
  324 K to 219 K 0.88 -> 0.89 / 0.88. A quadratic SIMEX was unstable at this reliability and was replaced.
- (D) Censoring: 320 boxes (17%) have no shadow component and were dropped. As upper limits (Tobit): size exponent 0.68 (0.59-0.74)
  vs 0.61; size-controlled poleward ratio 0.91 (0.81-1.02) vs 0.98. Some no-shadow boxes are not depressions, so this bounds the effect.
- Reading: the paper's conclusions hold, but the reported depth-size exponent (0.57-0.61) is a lower bound (0.7-0.9 after the diameter
  correction), the size-controlled poleward ratio should carry the censored range (0.91-0.98), and absolute depths carry a ~20%
  person-to-person scale uncertainty. In v6.12 and supplement v24 (Sections 3.7, 4.2, Figure 8; S9 and Figure S8).
- (E) Are the no-shadow boxes shallow hollows? (no_shadow_check.py -> no_shadow_check.json; Michael asked why the 'upper limits' sit at
  the bottom of the cloud.) They plot at each image's minimum measurable depth by construction. Hollows following the fitted depth-size
  relation (residual SD 0.54, or 30% more) would fall below the threshold in 26-70 of the boxes, against 320 observed; at 0.8-1.5 km,
  1-3 against 52. Most no-shadow boxes are not shallow hollows (non-depressions, offset boxes, low contrast), so the Tobit fit in (D)
  overstates the censoring effect and is an outer bound; the measured-only estimates (0.61, 0.98) are probably closer. In v6.12/v24.
