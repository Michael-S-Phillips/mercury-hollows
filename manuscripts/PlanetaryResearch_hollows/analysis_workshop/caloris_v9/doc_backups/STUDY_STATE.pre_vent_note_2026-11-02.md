# Study state — candidate hollow-forming volatile and its observational test

Updated 2026-08-28. Supersedes earlier phase summaries where they conflict.

## Where the study stands

Phases 1–3 complete. 56 candidate phases inventoried across six classes that the 2021 list of 57
did not cover, screened on thermal and athermal loss routes as **rate-versus-latitude curves**,
then on lag evolution, morphology, and geologic plausibility.

**Leading answer: elemental sulfur, as a concentrated and replenished reservoir.** It is the only
candidate whose loss rate exceeds the requirement imposed by the observed hollow age, and it scores
9/10 on the geologic rubric. But its weakest axis is *concentration*, and the open
question is whether that is met by local enrichment or by fracture-fed access to depth (below).

## The three constraints that did the discriminating

**1. The rate window, reset by the hollow age.** Wang et al. 2020 (10.1029/2020JE006559) give a
maximum model age of 103 (+200, −96) kyr for the global hollow population. With observed depths this
requires **6.9×10⁴ – 1.6×10⁷ m/Gyr**, about 3.2 orders above the top of the 3.3–130 m/Gyr window
used in the 2021 paper. Consequences:
  - Every athermal route fails. VUV photolysis of NaCl, even at a quantum yield ~50× the measured
    PSD anchor, reaches 56 m/Gyr — short by 3.1 orders. Native Na PSD short by 0.8; cusp sputtering
    by 5.8. **This eliminates photon- and particle-driven mechanisms as a class, on timescale
    grounds alone, independent of the unmeasured yields.**
  - C₆₀ and C₇₀ fall out (2.6×10³ and 3.1×10² m/Gyr at the cold limit).
  - Elemental S *exceeds* the requirement by 4.7 orders — the supply-limited regime.
  Caveat: the 103 kyr rests on a single probable impact crater across all sub-10 m/px NAC coverage
  of hollow floors. It is an upper limit on a thin statistical basis, and the authors say so. All
  screening is reported at **both** windows.

**2. Conservation of protons kills the carbon route.** Wang propose C + solar-wind H → CH₄. That
reaction cannot consume more carbon than there are protons to react with. At the cusp-peak proton
flux with *every* arriving proton forming methane, graphite removal tops out at **7.2×10⁻³ m/Gyr —
7.0 orders short**. Independently, their own global budget (1,266 km³ in ≤103 kyr) needs 5.7×10⁴¹
protons against 1.9×10³⁴ delivered: **7.5 orders short**, or 3,062 Gyr at perfect efficiency. Two
independent routes agreeing to 0.5 dex, with no free parameters. Diluting the carbon helps only
linearly (5.7 dex short at Mercury's measured 1.4 wt%); even a nonsensical 1:1 H:C stoichiometry
leaves 5.1 dex.
  **What this does and does not overturn.** Wang's *observations* stand — the graphite absorption is
  absent from hollow floors and present in background terrain, and hollows concentrate in LRM. What
  fails is the loss mechanism offered to explain them. Their spectral evidence is equally consistent
  with graphite as the **passive host removed along with excavated material**, which is the role
  assigned independently in Phase 1. Their own paper notes the CH₄ route needs an age near 12 kyr,
  making the demand 8.6× worse still.

**3. A concentration trade-off that binds our own answer.** Wang's eq 5 gives the integrated
volatile column removed; with their own halo parameters it reduces to h_volatile ~ 0.95 x hollow
depth at every depth (verified — it reproduces their published >22 m). So the volatile column is set
by the hollow depth, and the source volatile fraction is f = h_volatile / h_source, where the source
thickness is a **free parameter**. The budget fixes only the *product*:

| source thickness | required volatile fraction (median / deepest hollow) |
|---|---|
| 1x hollow depth | 95 / 99 vol% |
| 2x | 47 / 50 vol% |
| 5x | 19 / 20 vol% |
| 10x | 9.5 / 9.9 vol% |
| 20x | 4.7 / 5.0 vol% |

**The honest statement:** to source the volatile at Mercury's global-average sulfur (4 wt% = 5.4
vol%), the reservoir must extend **~18x deeper than the hollow, about 365 m below a median 20.8 m
hollow**. Either the volatile is locally enriched well above the planetary average, or hollows tap a
reservoir extending some hundreds of meters below them. The budget forbids both being false; it
cannot choose between them.

Two earlier framings are **withdrawn**: that ~20 vol% is a binding *floor* (it is the value for one
arbitrary geometry — a median hollow tapping a source as deep as the deepest hollow), and that the
exposed-layer model is *impossible* for median hollows (an artifact of holding the 22 m column fixed
while varying depth, which Wang's own equation couples). The claim that this supersedes the
pore-space floor is withdrawn too: the pore-space threshold (x_v* ~ 4 vol%, bracket 0–27%) is an
independent *floor* — below it, removal empties pores without lowering the surface — and both
constraints apply simultaneously.

Which branch holds is testable. Local enrichment implies a concentrating process (magmatic sulfide
segregation, exhalative deposition). A deep reservoir implies fracture-mediated access — the
cooling-fracture structural control Wang document — and this is independently favored by the
spatial statistics: hollows show a **significant crater association** (z = 3.4 at 5 km, z = 4.0 at
10 km) but **no near-field vent association** at any radius from 25–200 km.

## Reviewer 1's abundance objection, revisited

R1 objected in 2021 that natural sulfur abundances are too dilute to build near-pure layers in
impact-mixed regolith. The pore-space reply (4 vol%) settled it comfortably, because at 4 wt% S
Mercury already sits at 5.4 vol% — exceeding that floor with no enrichment at all.

The budget trade-off reopens the objection, though less sharply than I first claimed. It does not
demand a fixed 3.7x enrichment; it demands **either** local enrichment **or** a reservoir extending
~18x deeper than the hollow. R1's concern applies squarely to the first branch. The second sidesteps
it, but requires fracture-mediated access to depth — a different geologic claim needing its own
support, and one the crater/vent statistics happen to favor.

## Geologic plausibility (rubric: source, abundance, concentration, emplacement, survival; /10)

| class / entity | S | A | C | E | Sv | total | role |
|---|---|---|---|---|---|---|---|
| Elemental sulfur (α-S) | 2 | 2 | 1 | 2 | 2 | **9** | FORMER — only survivor under the Wang age |
| Crustal sulfides (oldhamite etc.) | 2 | 2 | 1 | 2 | 2 | **9** | RESERVOIR / geochemical marker |
| Graphite | 2 | 2 | 0 | 2 | 2 | **8** | HOST — CH₄ route rebutted |
| Alkali chlorides (halite, KCl) | 1 | 1 | 1 | 1 | 2 | 6 | eliminated — no route survives |
| Water ice | 2 | 2 | 2 | 0 | 0 | 6 | eliminated on latitude/thermal |
| Alkali sulfides | 1 | 0 | 1 | 1 | 1 | 4 | eliminated on abundance |
| Divalent/transition-metal chlorides | 1 | 1 | 0 | 1 | 1 | 4 | eliminated — residue armouring at 11–31 nm |
| Fluorides | 1 | 0 | 0 | 1 | 2 | 4 | eliminated — bandgaps miss Lyman-α |
| Fullerenes | 0 | 0 | 0 | 1 | 1 | 2 | eliminated — source unobserved |
| Organics | 1 | 1 | 0 | 0 | 0 | 2 | eliminated — radiolysis + location |
| Native alkali metals | 0 | 0 | 0 | 1 | 0 | 1 | eliminated — no bulk reservoir |

Both leading entities score **1/2 on concentration** — the axis the budget constrains. These scores were set against the withdrawn 20 vol% figure and are being re-scored against the trade-off; the ranking is not expected to change (concentration is a demand on sulfur under either branch), but the axis reasoning does.

## Spatial associations, and a null-model caution

Hollows cluster in fields, so a null randomising each hollow independently treats correlated points
as independent, understates the variance, and manufactures significance. The conservative null is a
rigid longitude rotation, preserving clustering and the latitude distribution.

- **Crater association: significant, and the strongest spatial result.** Under the rotation null,
  hollows sit inside a >=10 km crater far more often than chance. Two containment definitions were
  computed and reconciled: *nearest-crater* (inside only if the closest crater contains it) gives
  0.489 observed vs 0.177 +/- 0.047 null, **z = +6.7**; *any-crater* (inside if any crater contains
  it) gives 0.753 vs 0.454 +/- 0.082, **z = +3.7**. The difference is the definition, not the size
  cut — the any-crater fraction is stable at 0.75/0.73/0.68/0.62 across >=10/20/50/100 km cuts.
  Significant under both; the definition must be stated for either to be reproducible. The nearest
  crater is also larger than chance (39.8 km mean vs 21.9 +/- 4.7, z = +3.8).
- **No near-field vent association** at any radius from 25–200 km (at 25 km the observed fraction,
  0.24%, sits *below* the 0.86% null mean). A fumarolic point source predicts the concentration that
  is absent.
- **LRM: a result obtained and then withdrawn.** The naive null gave z = -4.2; under rotation the SD
  inflates 7x and z collapses to -0.6 (p = 0.28). Honest statement: *not supported by this coarse
  global raster* (0.12 deg/px), **not** "hollows are not associated with LRM" — the published
  association rests on higher-resolution local mapping, and Wang's local evidence is cited in
  preference to this null.

Together the crater and vent results point at **impact-fracture access** to a volatile reservoir
rather than fumarolic emplacement — independently the fracture-fed branch of the concentration
trade-off above.

**A methodological note for the PI.** `von_mises_fit.ipynb` uses independence-based longitude fits.
The vent and LRM-raster signals near-vanish only once clustering is preserved, so this affects the
*significance* of those two associations — not the fits themselves, and not the crater association,
which survives both nulls and both containment definitions.

## Resolution bias in the depth data (PI-flagged)

Measured depth correlates with image pixel scale at **ρ = +0.52** (p ≈ 10⁻¹⁷⁸), three times stronger
than the depth–latitude trend. Poleward of 70°N, **56% of cataloged hollows are shallower than one
pixel**. The *trend* survives control (partial ρ = −0.157) and holds within every resolution stratum,
but the *magnitude* is estimate-dependent: **1.9–4.0×** (the raw 3.98× is withdrawn). Mechanism
ranking is unchanged — thermal predicts 615×, cusp 0.04× (wrong sign), photolysis 4.46× — but
photolysis now sits at or above the top of the observed range rather than landing on it.

## Poleward limit — settled from the Bickel CNN catalog

**78.60°N**, superseding the 71.4° used throughout the 2021 paper and every brief in this study.
Derived from the Bickel et al. MDIS catalog (n = 19,110 at confidence ≥ 0.7 — matching the published
count exactly), supplied by the PI after the Bern repository moved to a network allowlist.

The raw catalog maximum is 88.49°, and it is **not** the limit: it is a single detection at incidence
90.17° — the Sun at or below the local horizon, so it cannot be a sunlit surface feature — whose
thumbnail shows diffuse texture with no bright-floored depression. The next two extremes (85.09,
84.75) are also grazing-incidence, low-confidence detections whose thumbnails are small bright blobs
consistent with sunlit slope facets. Only 0.22% of detections have incidence ≥ 85°, yet they supply
every |lat| above 79°: the poleward tail *is* the shadow/hollow confusion population.

After quality control (incidence < 85°, confidence ≥ 0.75, n = 11,166) the maximum is 78.60°,
corroborated by the ≥5-detection threshold (78.56°), the ≥10-detection threshold (78.12°) and the
99.9th percentile (78.01°). The value is stable at every confidence cut from 0.75 to 0.85.

**Convergence worth noting:** the human-mapped Blewett 2016 catalog (n = 2,518) gives 78.53° — within
0.07° of an independently constructed CNN catalog 7.6× larger.

**No verdicts flip.** Nothing has a thermal crossing latitude in (71.4°, 78.60°]. But the chloride
exclusion margins widen substantially — KCl 0.8° → 8.0°, FeCl₂ 2.2° → 9.4°, NaCl 4.5° → 11.7° — so
the halite exclusion that motivated the athermal extension is reinforced, not overturned. It also
narrows the headroom of the permissive survivors (C₆₀ at 83.0° retains 4.4°, C₇₀ at 81.7° retains
3.1°).

### The resolution bias is resolved, and it runs the other way

From the Blewett catalog, measured depth correlates with pixel scale at ρ = +0.52 — stronger than the
depth-latitude trend — raising the worry that poleward hollows were imaged more coarsely, inflating
their depths and masking a real decline. **The Bickel catalog shows the opposite:** Spearman(|lat|,
pixel scale) = −0.65 per detection and **−0.68 per unique MDIS image** (n = 4,915), with median
per-image pixel scale falling from 76.0 m at 0–15° to 13.9 m at 60–80°.

Poleward hollows are measured on *better* imaging. The observed poleward shallowing therefore cannot
be a coarse-resolution artifact, and the bias runs against the observation — finer images resolve
shallower features, so the true equator-to-pole depth contrast is if anything larger than measured.
This removes the principal caveat on the depth-latitude discriminant.

*Caveat:* computed on images that contain detections; images with no hollows are absent from the
catalog, so this is the resolution of hollow-bearing imaging rather than MDIS coverage in general.
The per-image and per-detection statistics agree closely, which argues against a detection-selection
effect, but full coverage modeling would need the MDIS index.

### Occurrence per imaged scene

Detections per unique imaged scene decline only **2.1×** equator-to-pole, and **non-monotonically** —
the 20–30° band is the global maximum at 2.06× the equatorial value. This is an *occurrence
probability*, not a rate, and must not be read as another vote for a loss mechanism merely because
2.1× sits near the cos-latitude photolytic 3.1×: for a supply-limited former, occurrence tracks where
the *reservoir* is accessible, not where the loss rate is high. The non-monotonic shape is what one
expects if hollows track specific geologic settings (LRM-bearing crater walls) rather than any smooth
latitude function of a driver.

### A data defect to report upstream

187 rows of the 0.7 catalog (33 of the 0.8) carry two space-separated tokens in the Confidence field,
e.g. `9 0.8466`, with all 22 fields otherwise present. I took the second token (all 187 land in
[0.70, 0.89], consistent with the catalog's own cut) and blanked the preceding y-coordinate field
rather than guessing. **The cause is not established:** a plausible digit-truncation hypothesis was
tested and found *vacuous* (appending digits to a value already in [0,1] trivially stays in range, so
the test had no power), and decimal-place counts in affected rows are indistinguishable from clean
rows. Worth reporting to the catalog authors.

## Non-volatile alternatives and the null hypothesis

All fail against the measured relief. Combined bare-surface erosion (sputtering + micrometeoroid
vaporisation) tops out at **0.27 m/Gyr** — 78 Gyr to excavate a median hollow, a shortfall of 5.9
orders against the median-depth/best-age case and 7.8 orders in the extreme bracket (deepest hollow,
shortest age). Impact-melt collapse fails the youth constraint. The optical-only null — hollows as a
reflectance phenomenon with no material removed — is refuted by the measured 20.8 m relief itself and
by Wang's >1,266 km3 loss budget.

## The observational test (answering question 2)

Of 110 hypothesis x observable cells: **33 diagnostic, 61 degenerate, 16 null**. Ranked
discriminating power: **S abundance (5) > local exosphere (4) > Cl/Na (3) >> thermal IR (1) >
reflectance (0)**.

The decisive measurement is **S/Si over hollow floors versus adjacent host terrain**, because the two
surviving hypotheses predict opposite signs:

| observable | sulfur predicts | Wang carbon predicts | instrument |
|---|---|---|---|
| S/Si over hollow vs host | elevated, in the top ~10 um | **no S enrichment** | MIXS (S Ka 2.308 keV) |
| C abundance | decrease if graphite co-removed | decrease (Wang's key evidence) | GRS / NS |
| local exosphere | S/S2 enhanced, interannual variability | CH4/C2, flux tiny | PHEBUS, Strofio |
| 500–600 nm graphite band | featureless brightening | loss of darkening feature | MASCS VIRS |

Note the asymmetry: carbon predicts *no* sulfur enrichment while sulfur predicts a clear one, so a
single MIXS measurement discriminates rather than merely checks consistency. The carbon-side
observables are **degenerate** — graphite-as-passive-host reproduces Wang's spectral evidence without
his loss mechanism.

**This is a negative result for the study's original premise.** 7–14 um emissivity is degenerate and
reflectance carries no discriminating power: **MERTIS cannot settle the hollow-former question.** The
instruments that can are MIXS and the exosphere package.

## Referee objections from the 2021 paper: 10 answered, 5 partial, 2 open

Open: Nathair Facula shows S and C *depleted* by explosive eruption, which bears on the supply term
of a replenished reservoir; and the S budget considered only explosive pyroclastics, neglecting
larger-volume effusive volcanism (would enlarge the reservoir, favorable but unquantified).

## Outstanding

1. **Poleward limit — CLOSED.** Bickel catalog obtained and analyzed; limit is 78.60°N with the
   depth-resolution caveat resolved in the favorable direction. Depths still come from Blewett
   (n=2518), since the Bickel product is a detection catalog with no depth column — a proper
   completeness-corrected *depth* model would need morphometry this catalog does not carry.
2. **The enrichment-versus-depth branch** — the sharpest open question. Either a concentrating
   process raising sulfur well above the planetary average, or fracture-mediated access to a
   reservoir ~18x deeper than the hollow. The spatial statistics favor the second; it needs
   quantitative development (fracture aperture, connectivity, transport rate).
3. **Phase 4 (detectability) not yet run**, but now well-posed: the target is the S/Si contrast
   between hollow floors and adjacent host terrain at MIXS spatial resolution, plus the predicted
   local exospheric S signal. MDIS image products remain unreachable (server-side refusal, not a
   proxy block).
4. **Two referee objections still open**, both on sulfur *supply* rather than loss physics: Nathair
   Facula shows S and C depleted by explosive eruption (unfavorable for juvenile outgassing as the
   recharge term), and the 2021 budget counted only explosive pyroclastics while neglecting
   larger-volume effusive volcanism (would enlarge the reservoir, unquantified).
5. **Unmeasured yields remain bracketed**, but this now matters less: under the Wang window the whole
   athermal class fails on timescale regardless of where in each bracket the true value sits.
6. **Ten of 16 scored candidate entities rest on a process not documented on Mercury** (flagged in
   the CSV). For the surviving candidate the load-bearing unmeasured quantity is the sulfur
   concentration step.
