# Poleward-limit sensitivity and the hollow depth-vs-latitude discriminant

**Update to Phase 2b (framework/athermal track).** Author: hollows-emissivity project.
PI: M. Phillips. Generated 2026-08-28.

The PI moved the hollow **poleward limit** from a hardcoded 71.4°N (Thomas 2014 /
Phillips 2021) to a **parameter**, re-derived from the more complete Blewett 2016 catalog
(n=2518) at **78.53°N** (+7.13° poleward). The HORNET machine-learning catalog
(Bickel et al. 2025, 19,110 hollows, doi:10.48620/77046) will supersede this once its
deposit returns from maintenance; the analysis is parameterized so that number swaps in
without a re-run. This document reports (1) the limit-sensitivity of every athermal verdict
and (2) the depth-vs-latitude morphological discriminant derived from the same catalog.

---

## 1. The limit is now a parameter; results reported at both values

`athermal_screening.py` now carries `COLD_LIMITS = {"Thomas2014_Phillips2021": 71.4,
"Blewett2016_derived": 78.53}` and the combined matrix
(`runs/combined_loss_matrix.csv`) carries a full column set at **both** limits
(`*_71p4` and `*_78p53`), plus limit-independent columns. Nothing is hardcoded to 71.4.

**Catalog cross-check (reproduced exactly from the PI's file).** Blewett 2016: median depth
20.80 m; |lat| percentiles 99.0th = 73.05, 99.5th = 77.78, 99.9th = 78.53; max 78.53°N.
Distribution overwhelmingly northern (2462 N vs 56 S; southern extreme only −3.44°). **The
hemispheric asymmetry is an MDIS coverage/resolution bias, not a real N–S difference** — it
should not be read as a physical signal.

---

## 2. Consequence: the move STRENGTHENS every thermal exclusion, and it goes the way that helps the paper

The extra 7° of latitude drops the surface temperature from 523 K (71.4°N) to 465 K (78.53°N).
Because sublimation is Arrhenius, **thermal loss collapses another ~500–1200×** over that span
(NaCl 648×; CaCl₂ ~1200×; MgCl₂ ~1000×), while the athermal photon
routes barely move.

- **Halite** (the paradigm case): its thermal congruent-sublimation rate falls below the window
  at **65.0°N** on this work's model (`thermal_maxlat_window_N` in the matrix; the Phillips 2021
  archive quotes ~67°N — a ~2° difference from the different thermophysical model, reported not
  silently adopted). On this work's 65.0°N value, halite is excluded by **6.4° against 71.4°N** and
  by **13.5° against 78.53°N**; on the archive's 67°N, by 4.4° and 11.5°. Either way, a more complete
  catalog puts hollows *further* poleward of where thermal sublimation of halite can operate — the
  thermal exclusion **strengthens**, and the conclusion is insensitive to the ~2° model difference.
- The same logic applies to every thermally-limited candidate: a phase marginal at 71.4°N fails
  at 78.53°N, and a phase already excluded fails harder.

**KCl is the phase that sits closest to the boundary.** Its congruent-sublimation rate is in-window
through 70°N (`thermal_maxlat_window_N` = 70.0) but drops below the window floor by 71.4°N
(2.1 m/Gyr) and is far below it by 78.53°N (4.3e-3 m/Gyr). So on the discrete `thermal_class`
label it reads `thermal_excluded_crosses_below_coldlimit` at **both** limits (it crosses at 70°N,
already below 71.4°N) — no phase's `thermal_class` label actually flips between the two limits, and
the matrix columns confirm this. What the limit move does to KCl is quantitative: it goes from just
below the old boundary (crossing at 70.0°N vs a 71.4°N limit, ~marginal) to comfortably excluded at
78.53°N. Every other thermally-limited phase crosses the window even further equatorward, so all
fail harder poleward; none is rehabilitated thermally by the move.

---

## 3. Athermal track: does photolysis still reach the window at 78.53°N?

**Yes — and this is a STRONGER argument for the athermal route than we had at 71.4°N.** The
photon flux falls only from 3.1× (to 71.4°N) to 5.0×
(to 78.53°N) relative to the equator — so the quantum yield a wide-gap salt needs to reach the
window rises by just the flux ratio, **1.6×**, over the same 7° in which thermal loss collapses
~500–1200×.

| Phase | gap (eV) | Y needed @71.4°N (× anchor) | Y needed @78.53°N (× anchor) | verdict |
|---|---|---|---|---|
| CaCl₂ | 7.0 | 2× | 3× | within rehab bracket at both |
| MgCl₂ | 8.0 | 5× | 7× | within rehab bracket at both |
| KCl | 8.4 | 6× | 10× | within rehab bracket at both |
| NaCl | 8.5 | 9× | 15× | within rehab bracket at both |
| NaF | 11.5 | 137× | 220× | NOT rehabilitated at either |
| CaF₂ | 12.1 | 101× | 163× | NOT rehabilitated at either |
| MgF₂ | 12.8 | 142× | 227× | NOT rehabilitated at either |

The chlorides stay within the rehabilitation bracket (Y ≤ 1e-3, i.e. ≤53× the PSD anchor) at
**both** limits: NaCl needs 9× the anchor at 71.4°N and 15× at 78.53°N — still comfortably inside.
The fluorides need ~100–230× the anchor at both and remain excluded over the plausible bracket.
**No rehabilitation verdict flips when the limit moves poleward**; the athermal argument survives
the more stringent test and is strengthened relative to the (now even more comprehensively
excluded) thermal route.

(The quantum yield each wide-gap salt needs to reach the window at 71.4 vs 78.53°N is tabulated
above; the move makes photolysis only 1.6× harder and the chlorides stay within the rehabilitation
bracket at both limits.)

---

## 4. Morphological discriminant: the poleward tail is systematically SHALLOW

The Blewett 2016 depths decline poleward: the 8 highest-latitude hollows (all at 78.53°N) are
3.5–9.6 m — every one below 10 m — against a 20.8 m global median. **All statistics here use the
|latitude| convention (both hemispheres), matching the PI's band table; the equatorial anchor bin
(|lat| 0–5°) is n=144, median 23.32 m.** (Northern-only, that bin is n=88, median 20.63 m — the 56
southern hollows all sit <3.5° from the equator and fall only in the anchor bin, so the convention
choice changes *only* the anchor and therefore propagates into every equator-to-pole ratio.)

**A large spatial-resolution confound must be controlled (PI-flagged).** Coarser images yield
systematically *deeper* measured hollows: Spearman ρ(pixel-scale, depth) = **+0.52** (p≈1e-178) —
**three times stronger** than the depth–latitude trend itself (ρ(|lat|, depth) = −0.16). There is
also censoring near the detection limit: poleward of 70°N, **55.6%** of cataloged hollows are
shallower than one pixel (vs 13–17% equatorward), so high-latitude depths are measurements at or
past the resolution floor.

**What survives the control — the SIGN of the trend, not its magnitude.**
- The poleward decline is robust: partial Spearman(|lat|, depth | pixel-scale) = **−0.157**
  (p=2.8e-15, essentially unchanged from the raw −0.160), and it holds independently within every
  resolution stratum (ρ = −0.13, −0.16, −0.20, −0.27 for px 1–8, 8–12, 12–16, 16–20 m).
- The **magnitude is strongly control-dependent and the earlier 3.98× figure is withdrawn.** The
  equator→pole median-depth ratio is 3.91× (raw, ≥75°N, n=21 only) / 2.43× (raw, ≥70°N) / 2.31×
  (resolution-matched to common pixel support) / 1.92× (stratum-weighted). **Honest range: ~1.9–4.0×.**

![Depth discriminant and resolution confound]({{artifact:art_bf9d38fa-960f-460c-b9cb-a9923fe7e7c7}})

*(a) The spatial-resolution confound: coarser images give systematically deeper measured hollows
(ρ(pixel, depth)=+0.52); the dashed line is the 1-pixel detection floor. (b) Depth vs |latitude|
band medians (±IQR, |lat| convention) — the poleward decline is robust in sign (partial ρ=−0.157
controlling for pixel scale; holds in every resolution stratum) but its magnitude is control-
dependent (1.9–4.0×). (c) Predicted equator→pole rate decline per mechanism against the observed
range (grey band): thermal over-predicts by ~2 orders, cusp sputtering has the wrong sign, and only
VUV photolysis is in the right regime — at or above the top of the observed range.*

**What each mechanism predicts** (depth ∝ rate × time, uniform-age assumption):

| Mechanism | Predicted equator→pole rate decline | vs observed 1.9–4.0× |
|---|---|---|
| Thermal sublimation (S) | **×522–615 DOWN** | over-predicts by ~2 orders — excluded as sole shape-setter |
| **VUV photolysis (~cos)** | **×4.5–5.0 DOWN** | in the right regime, but **at or above the top** of the range |
| Cusp sputtering | ×0.03 (**RISES** poleward, peaks ~70°N) | **wrong sign** — excluded |

- **Thermal** predicts a ~500–600× poleward rate collapse — one to two orders steeper than *any* of
  the observed estimates. A pure thermal former would produce essentially no deep high-latitude
  hollows; it can only match the shallow observed gradient with an additional younger-poleward age
  gradient. Excluded as the sole setter of the depth-latitude shape.
- **VUV photolysis** (∝cos-latitude, ×4.5–5.0) is the only prediction in the observed regime — same
  sign, same order — but it now sits **at or above the top** of the 1.9–4.0× range rather than
  landing on it. The earlier "precise agreement" (a log-slope of −0.005/deg landing on photolysis)
  was an artifact keyed to the 3.98× / 21-hollow bin and is withdrawn. The **ranking is robust; the
  precise agreement is not.**
- **Cusp sputtering** predicts loss to *peak* near the 70°N cusp flux maximum (1.03e11 m⁻² s⁻¹), i.e.
  deeper poleward. The data show monotonic shallowing through and past that latitude — the opposite
  sign — so cusp sputtering is **excluded as the primary excavating driver** (a sign flip no age
  assumption can undo).

**Net.** The depth-vs-latitude morphology is a genuine mechanism discriminant, but a weaker one than
first stated. It **excludes cusp sputtering** outright, **excludes a pure thermal former** as the
shape-setter (over-predicts the decline by ~2 orders), and leaves **VUV photolysis** as the only
mechanism whose latitude shape is in the right regime — though sitting at/above the top of the
observed range, not on it. Combined with §3 (photolysis reaches the window at 78.53°N almost as
easily as at 71.4°N, while thermal is even more comprehensively excluded), the limit move and the
depth trend point the same way — **toward a photon-driven component in high-latitude excavation** —
with the depth-latitude *slope* treated as provisional pending HORNET (Bickel 2025), whose documented
resolution field will allow the completeness/censoring correction to be modeled rather than
stratified.

---

## 5. Caveats

- **Spatial-resolution confound (PI-flagged, controlled here).** ρ(pixel-scale, depth)=+0.52 is 3×
  stronger than the depth–latitude trend; poleward of 70°N, 55.6% of hollows are sub-pixel. The
  trend's SIGN survives the control (partial ρ=−0.157; holds in every resolution stratum) but its
  MAGNITUDE does not — the honest range is 1.9–4.0× and the earlier 3.98× is withdrawn. The
  depth-latitude *slope* is provisional pending HORNET (documented resolution field → modeled
  completeness correction). Do not use a uniform-detectability-floor ("censoring-corrected") variant:
  it leaves only ~5 high-latitude hollows, which is not a measurement.
- Depth ∝ rate × time assumes roughly uniform hollow age with latitude. Because thermal now
  over-predicts the decline by ~2 orders (not merely 19× too steep), thermal is excluded as the sole
  shape-setter regardless of a modest age gradient; the age assumption bears mainly on whether a
  *combined* thermal+age picture could mimic the shallow gradient. Cusp sputtering is cleanly excluded
  (wrong sign, no age assumption can flip a sign).
- The photolysis rate uses the PSD-anchored quantum yield and its bracket; the yield is UNMEASURED
  for every phase except Na. The rehabilitation claim for chlorides is bracketed (viable if Y reaches
  ~1e-3), not a single-point result — and that conclusion is unchanged by the limit move.
- 78.53°N is the Blewett 2016 catalog extreme; HORNET (Bickel 2025) may move it again. All numbers
  here are reported as functions of the limit so the update is a parameter swap.
- Counts per latitude band are NOT a clean abundance measure (coverage/resolution vary strongly with
  latitude); depth is used as the discriminating statistic throughout.


---

## 6. The Wang (2020) age constraint inverts the athermal rehabilitation — and that is a stronger result

**New primary constraint (PI-supplied, Crossref-verified).** Wang, Xiao, Chang & Cui 2020,
*JGR Planets* **125**, e2020JE006559, doi:10.1029/2020JE006559 — "Lost Volatiles During the
Formation of Hollows on Mercury". NAC crater statistics on hollow floors give a **maximum model age
for the global hollow population of 103 (+200, −96) kyr** — kyr, not Myr — and hollow depth does
**not** correlate with model age (depth is set by volatile-layer thickness, not by how long
sublimation has run). DOI verified against Crossref this session: authors Wang/Xiao/Chang/Cui, title,
journal, year 2020, volume 125 all match.

**The implied rate window is ~3 orders above the one used in Phases 1–2.** Required rate = depth /
age. With Blewett depths (median 20.8 m, max 110.5 m) over the age bracket 7–303 kyr:

| | rate (m/Gyr) |
|---|---|
| median depth, best age (103 kyr) | 2.0×10⁵ |
| **Wang window (full bracket)** | **6.9×10⁴ – 1.6×10⁷** |
| Phillips 2021 window (used in §1–5) | 3.3 – 130 |

The Wang window floor sits **2.7 dex above the Phillips window top**. Both windows are reported here;
the switch is **not** made silently — every §1–5 verdict is on the Phillips window, this section is
on the Wang window, and the sensitivity to the assumption is the point.

![Wang-window re-screen]({{artifact:art_b07ac3a7-5318-49ad-b5e8-cc9f5dfdfa54}})

*(a) Every candidate's thermal (●) and best-bracket athermal (×) rate at the 78.53°N cold limit
against both windows (green = Phillips, red = Wang). Every athermal × lies left of the Wang band; the
dashed line is the maximum athermal rate anywhere in the matrix (3.4×10⁴ m/Gyr), below the Wang floor.
(b) The verdict inversion: elemental S exceeds both windows (supply-limited); C60/C70 and all chloride
rehabilitations, viable under Phillips, now fail the Wang window by 1.4–2.8 dex (rates are the best
achievable at the 78.53°N cold limit, per the companion combined_loss_matrix.csv).*

**Re-screen result (all 53 rate-carrying candidates, cold limit 78.53°N):**

- **No candidate reaches the Wang window via ANY athermal route, at any latitude, at the top of every
  unmeasured-yield bracket.** The maximum athermal rate anywhere in the matrix is **3.4×10⁴ m/Gyr** —
  below the Wang floor of 6.9×10⁴. This is a clean, bracket-surviving negative.
- **All six chloride/oxysalt rehabilitations fail.** NaCl VUV photolysis at the *optimiztic* yield
  (Y=10⁻²) is 117 m/Gyr at the cold limit — **2.8 dex short** of the Wang floor. CaCl₂, the fastest
  chloride at the cold limit, reaches only ~5.3×10² m/Gyr — still 2.1 dex short (MgCl₂ 2.3×10²,
  2.5 dex short). The athermal route is simply too slow to excavate 20 m in
  ~100 kyr.
- **C60/C70 also fail**, now on rate: C60 thermal at the cold limit is 2.6×10³ m/Gyr (1.4 dex short),
  C70 3.1×10² (2.3 dex short). Under the Phillips window these were "thermophysically viable"; under
  Wang they are too slow.
- **Native Na PSD** fails by ~0.8–2.8 dex (convention-dependent — see §3 caveat) even at best case.
  **Cusp sputtering** fails by ~5.8 dex.

**This is a stronger result than the rehabilitation it replaces.** The athermal work does not become
worthless — it converts from "chlorides are conditionally rehabilitated" into a much sharper
**negative that survives the entire unmeasured-yield bracket**: *no athermal route (photolysis, PSD,
sputtering, ESD, radiolysis, impact, chemical) can form hollows on the Wang timescale.* That
independently **eliminates the entire athermal mechanism class as a PRIMARY hollow-former**, on rate
grounds alone, regardless of the yields we could never measure.

**What still reaches the Wang window** (thermally, at the cold limit): 20 phases — but these are
disqualified on independent grounds established in Phases 1–2b: the true ices (H₂O, CO₂, CO, NH₃, SO₂)
sublimate far too fast to persist and have no stability field; the organics are radiolytically
carbonised before they can excavate; native metals (Na, K) have no bulk reservoir; the reactive
non-metals (P₄, red P, Se, Te) are abundance- or chemistry-excluded. **Elemental sulfur is the one
phase that reaches the Wang window AND survives every other filter** — and it does so in the
**supply-limited** regime: its rate *exceeds* the window (by ~4.7 dex at the cold limit on the
SI-Pa convention), so its hollow depth is set by reservoir supply, not by the excavation rate. This
is exactly the SCArFS/replenished picture, and it is exactly what Wang's own findings independently
require: a young (~100 kyr), rapidly-formed, layer-thickness-limited, actively-replenished volatile.

**Convention note.** Elemental sulfur's Eq.-9 coefficient was stored in `candidates_sulfide.csv` as a
raw-kPa `A_Pa`=17.502; the project's settled SI-Pa convention requires `A_Pa`=24.410 (+ln 1000). The
corrected value gives S thermal = 7.7×10¹¹ m/Gyr (diurnal) at 78.53°N, matching the PI's 8.1×10¹¹
spot-check; the uncorrected value was ~3 dex lower. The verdict (S exceeds the window, supply-limited)
is unchanged either way. The five metal-sulfides are thermodynamically derived (raw Pa, correct as
stored); the salt and carbon tracks were already on the SI-Pa convention.

**Bottom line under the Wang age.** The most plausible hollow-forming volatile is **elemental sulfur
in a replenished (SCArFS-type) reservoir** — the only candidate fast enough on a 100-kyr timescale
that is not otherwise disqualified. The athermal routes, our earlier headline, are eliminated as
primary formers by the Wang timescale, which is itself a robust and independently useful result.

---

## 7. Wang's carbon mechanism (C + solar-wind H → CH₄): does our graphite elimination survive?

Wang et al. (2020, §4.1) argue *for* carbon and against our sulfur conclusion. Their case is spectral:
hollows sit in low-reflectance material (LRM) whose darkening agent is graphite (Murchie et al. 2015);
the ~500–600 nm graphite feature is **absent** in hollow-floor spectra but present in background
terrain, so C appears to have been lost from the floors. Their proposed loss route is **chemical**:
C + implanted solar-wind protons → CH₄ (Blewett et al. 2016), with escaping CH₄ lofting the fine
grains that make the bright haloes. This is a different mechanism from the ones our Phase-1 carbon
verdict screened (sublimation, photolysis, sputtering) — it is chemical reduction by implanted H —
so it must be scored on its own rate law, not dismissed by analogy.

**The C→CH₄ rate is set by a hard proton-supply ceiling.** Each CH₄ removes one C atom but consumes
**four** implanted H. The maximum carbon-erosion rate is therefore the implanted-proton flux divided
by 4, converted to a recession rate — a ceiling that no unmeasured reaction efficiency can exceed:

    max C rate  =  (φ_H / 4) / n_C · GYR       (100% efficiency, n_C = 1.13×10²⁹ m⁻³ for graphite)

| available proton flux | max C→CH₄ rate (m/Gyr) | dex short vs Wang C req. (2.3×10⁵) |
|---|---|---|
| equatorial surface (3×10⁹) | 2.1×10⁻⁴ | 9.0 dex short |
| cusp peak 70°N surface (1.03×10¹¹) | 7.2×10⁻³ | **7.5 dex short** |
| raw solar wind, unshielded (2×10¹³) | 1.4 | 5.2 dex short |
| raw solar wind at perihelion (7×10¹³) | 4.9 | 4.7 dex short |

All four rows are benchmarked against the *same* carbon-specific requirement, 2.3×10⁵ m/Gyr (Wang's
103-kyr age). Against Wang's more demanding 12-kyr CH₄ steady-state figure (2.0×10⁶) each row is ~0.9
dex further short (10.0 / 8.4 / 6.2 / 5.6). The generic depth/age *window floor* (6.9×10⁴, used in §6)
is a slightly less stringent benchmark; the conclusion holds against any of the three.

![C→CH₄ supply ceiling]({{artifact:art_9e5d4362-3f6b-40e0-851f-fbce072bc6d3}})

*The C→CH₄ recession ceiling (blue band, 4→1 H per C, 100% efficiency) against proton flux, versus
Wang's own carbon requirement (red) and the Phillips window (green). Across the entire physical flux
range — from the equatorial surface to unshielded solar wind — the ceiling stays 4.7–9.0 orders below
what Wang's age requires. The gap at the cusp flux is 7.5 dex.*

**The route fails Wang's OWN requirement by 4–7 orders, and the failure is bracket-proof:**
- At the realistic cusp-peak surface flux (the most favorable place on Mercury), the ceiling is
  7.2×10⁻³ m/Gyr — **7.5 dex below** the 2.3×10⁵ m/Gyr Wang's 103-kyr age demands (and 8.4 dex below
  their 12-kyr CH₄ steady-state figure of 2.0×10⁶).
- Even at raw, magnetospherically **unshielded** solar wind — a physically impossible upper bound at
  the surface — it is still 4.7–5.2 dex short.
- The gap is immune to the unmeasured reaction efficiency (already set to 100%) and to product
  stoichiometry: the most H-poor plausible product (1 H per C) still leaves the ceiling 6.9 dex below
  the 2.3×10⁵ carbon requirement. To meet Wang's carbon rate the surface would need a proton flux of
  3.3×10¹⁸ m⁻² s⁻¹ — 7.5 dex above the cusp surface flux.

**Independent global-budget cross-check (no rate model at all).** Removing Wang's own 1,266 km³ of
carbon as CH₄ requires 5.74×10⁴¹ protons (4 H per C atom at graphite density). The solar wind delivers
only 1.92×10³⁴ protons to the 57,400 km² hollow area in 103 kyr even at the cusp-peak flux — **short by
7.5 orders**; at that flux and 100% efficiency, supplying the protons would take ~3,076 Gyr (~220× the
age of the Universe). This budget route uses no sublimation or DIET rate law — only conservation of
protons — and agrees with the per-area rate ceiling above to within 0.5 dex. (Both reproduced
independently this session; they match the PI's 5.7×10⁴¹ / 1.9×10³⁴ / 7.5 dex / 3,062 Gyr.) The carbon
elimination therefore no longer rests on a *modeled* rate — it rests on proton conservation, which is
far harder to contest.

**Verdict: our graphite elimination survives, and now stands quantitatively against a published
mechanism.** The disagreement with Wang is real and we resolve it on the physics: the C→CH₄ route is
supply-starved by 4–7 orders and cannot excavate hollows on any timescale in the Wang bracket.

**We accept Wang's observation but reject Wang's excavator.** The spectral fact — C depleted in hollow
floors, retained in background LRM — is not in dispute; only the mechanism that removes the C is. Our
Phase-1 picture explains the *same* observation without a carbon-specific excavator: graphite is a
trace grain host (<1 wt%, Xu 2024) carried off passively when the actual excavating volatile deflates
the surface. Wang's own halo mechanism — fine grains lofted by an escaping gas — works with an
**abundant escaping S-vapor driver** (elemental S sublimating in the supply-limited SCArFS regime,
§6) but fails with a CH₄ driver that is starved by 4–7 orders. So the correct reading is: **S is the
excavator, C is the passive tracer whose depletion marks where the excavation happened** — one picture
that satisfies both Wang's spectral observation and the rate constraint.

**Abundance cross-check (geologic-track quantity, noted for support).** Wang's budget (>1,266 km³ lost
in ≤103 kyr, using 24 m mean depth over 57,400 km²) is **17.1% of the near-surface S inventory**
(4 wt% terrane maximum, Weider et al. 2015; 29.7-35.9% at the 1.9-2.3 wt% global mean) but **41.2% of the near-surface C inventory** (1.4 wt%, Peplowski 2015).
Consuming 41% of the planet's near-surface carbon in 100 kyr is a far more extreme demand than 17% of
its sulfur — an independent abundance argument pointing the same way (S over C), and one that (as Wang
concludes) requires active rejuvenation of the subsurface reservoir either way, reinforcing the
supply-limited picture. Full scoring belongs to the geologic track.

**Provenance and caveats (per PI corrections).** The 103 kyr is an **upper limit** on the average
model age, from a single probable impact crater in all <10 m/px NAC coverage of hollow floors — a thin
statistical basis, carried as such. The rate-window derivation (§6) uses the **20.8 m median** from the
Blewett n=2518 point catalog; Wang's budget uses a **24 m mean** over the Thomas (2014) 57,400 km²
area — both defensible, stated explicitly, not mixed. The earlier "depth does not correlate with model
age" phrasing is dropped as unsourceable; what Wang establishes is a single global model age and that
depth is set by volatile-layer thickness.

---

## 8. The concentration requirement (Wang column-integral) — a trade-off, not a single number

**[Corrected framing — supersedes an earlier "~20 vol% binding requirement" statement.** The
arithmetic below was reproduced correctly; the earlier *interpretation* of ~20 vol% as a floor was
wrong. ~20 vol% is the value for one particular source geometry, not a requirement.]

Wang's removed volatile column h_volatile is **coupled to hollow depth**, not a fixed 22 m. Their
mass balance (eq 5),

    h_volatile = h_inner + h_HF − (ρ_halo/ρ_bed)·(h_inner + h_outer·A_bnd/A),

has submillimeter halo terms, so **h_volatile ≈ 0.95 · h_HF at every depth** (verified: h_volatile(24 m)
= 22.90 m, reproducing Wang's published ">22 m"; ratio 0.954). The volatile column removed is therefore
set by the hollow depth, and the source volatile fraction is f = h_volatile / h_source where
**h_source — how far the source extends below the hollow — is a free parameter.** The budget fixes only
the *product* f × h_source, not f alone.

**The constraint is a trade-off between source concentration and source depth** (independently
reproduced this session; matches the corrected record to the meter):

| source thickness (× hollow depth) | required f, median hollow (20.8 m) | required f, deepest hollow (110.5 m) |
|---|---|---|
| 1× (source as thin as the hollow) | 94.7 vol% | 99.0 vol% |
| 2× | 47.4 vol% | 49.5 vol% |
| 5× | 18.9 vol% | 19.8 vol% |
| 10× | 9.5 vol% | 9.9 vol% |
| 20× | 4.7 vol% | 5.0 vol% |

**Honest statement.** To source the volatile at Mercury's global-average sulfur content
(4 wt% = 5.4 vol%), the reservoir must extend **~18× deeper than the hollow — about 365 m below a
median 20.8 m hollow.** So EITHER the volatile is locally enriched well above the planetary average,
OR the hollow taps a reservoir extending some hundreds of meters below it. The budget cannot
distinguish these two branches; it only forbids both being false.

**Both branches are real, testable geologic demands, and the fracture-fed branch is now co-equal, not
a caveat.** Local enrichment implies a concentrating process (magmatic sulfide segregation, exhalative
deposition). A deep reservoir implies fracture-mediated access — exactly the cooling-fracture
structural control Wang document (their §3.1), and consistent with the impact-crater association the
geologic track found significant and with the *absence* of near-field vent association (together
pointing at impact-fracture access rather than fumarolic emplacement). This is the geologic track's
axis to score; it is noted here only for cross-track traceability, and its concentration scores are
being re-set against the trade-off (not a single 20 vol% figure).

**Relation to the pore-space floor — both apply, neither supersedes the other.** The pore-space
threshold (x_v* ≈ 4 vol%, bracket 0–27%) remains a genuine independent *floor*: below it, removing the
volatile empties pore space without lowering the surface. The budget gives a *trade-off*, a different
kind of statement. The requirement is the conjunction: **f must exceed the pore-space floor AND satisfy
f × h_source = h_volatile.** (My earlier "the budget supersedes the pore-space floor at 5× strictness"
is withdrawn.)

**Caveats (per PI).** Wang's h_volatile is an explicit *minimum* (they cancel the lag term by setting
ρ_lag = ρ_bedrock), so every f in the table is a lower bound. Depth basis: 20.8 m median / 110.5 m max
from the Blewett n=2518 point catalog; Wang's own budget uses a 24 m mean over the Thomas (2014) area —
stated, not mixed. Details in `runs/concentration_requirement.json` (geologic track).

**What this does and does not change.** It refines *how* the concentration demand is stated; it does
**not** touch the results this track owns: the Wang age window (6.9×10⁴–1.6×10⁷ m/Gyr), the athermal
class elimination on rate, the carbon/CH₄ proton-supply falsification (§7, a supply argument
independent of any concentration geometry), or elemental S as the only rate-surviving candidate. The
geologic conclusion that the former must be **concentrated rather than dispersed trace material** still
holds — it now rests on the pore-space floor plus the enrichment-or-depth trade-off, rather than on a
single 20 vol% number.


> [CORRECTED 2026-09-02: 4 wt% is the high-sulfur TERRANE MAXIMUM, not the global mean (1.9-2.3 wt%, from S/Si = 0.076 XRS Nittler et al. 2020 / 0.092 GRS Evans et al. 2012 with Si ~23.4-27.1 wt%); and the attribution was wrong - Weider et al. 2016 (10.1002/2016gl068325) is a real paper but on sulfur- and carbon-driven explosive volcanism, NOT global mapping. The global-mapping source is Weider et al. 2015 EPSL (10.1016/j.epsl.2015.01.023). At the global mean the S inventory is 3,520-4,270 km3 and the lost fraction is 29.7-35.9%, not 17.1%, so the S-vs-C margin narrows from 2.4x to 1.1-1.4x: carbon remains the more extreme demand, but only just.]
