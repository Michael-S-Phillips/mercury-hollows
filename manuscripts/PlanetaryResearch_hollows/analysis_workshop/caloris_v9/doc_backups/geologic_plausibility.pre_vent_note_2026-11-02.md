# Geologic plausibility of candidate hollow-forming volatiles (Phase 3a)

**Question.** Physics narrowed the field of candidate volatiles by loss-rate vs latitude. Geology
must now do the discriminating: a phase is only a credible hollow-former if Mercury can **make** it,
**concentrate** it into a losable body, and **emplace** it where hollows are actually seen, and if
that reservoir can **survive** to the present. This is the logic of Phillips et al. 2021
([10.1016/j.icarus.2021.114306](https://doi.org/10.1016/j.icarus.2021.114306)).

Every DOI cited here is Crossref-verified in `docs/verified_bibliography.md` (all 26 checked this
session against that file; none newly recalled). Scores are **geologic** and are stated plainly
where they disagree with the physics ranking.

---

## 1. Spatial associations — what the observations actually support

The single most important methodological point is inherited from the null-model work already run on
this project's data and is **non-negotiable**: hollows **cluster in fields**, so a null that
randomises each hollow independently (latitude-matched random longitude) treats 2518 correlated
points as independent, understates the variance, and manufactures significance. **All significance
below uses the conservative rigid-longitude-rotation null**, which rotates the whole catalog by a
random longitude and so preserves both the clustering and the latitude distribution. Data:
`runs/spatial_associations.json`.

### 1.1 Impact craters — the positive discriminant (this session)

Test: hollow centroid → nearest crater, against the Herrick/Fassett-lineage global crater catalogs
(`2017je005516r_*`; n = 31 612 ≥ 5 km, 17 031 ≥ 10 km), rotation null N = 1000.

| statistic | observed | rotation null | z | verdict |
|---|---|---|---|---|
| frac. **inside** a ≥10 km crater (≤1 radius, *nearest-crater* def.) | **0.489** | 0.177 ± 0.047 | **+6.7** | highly significant |
| frac. within 2 crater-radii | 0.565 | 0.364 ± 0.064 | +3.1 | significant |
| frac. within 10 km of any ≥5 km crater | 0.228 | 0.114 ± 0.029 | +4.0 | significant |
| mean diameter of nearest crater (km) | 39.8 | 21.9 ± 4.7 | +3.8 | significant |
| frac. within 25 / 50 / 100 km | 0.53 / 0.93 / 1.00 | 0.52 / 0.93 / 1.00 | ≈0 | n.s. (craters saturate) |

**Result: ~49% of hollows lie inside the rim of a large (≥10 km) crater versus ~18% expected under
the conservative null (z = +6.7), and the crater a hollow sits near is systematically larger than a
random meridian would place it (39.8 vs 21.9 km, z = +3.8).** The 25–100 km radii are not
significant only because ≥5 km craters saturate those scales everywhere. This is the strongest,
most robust spatial association in the global data, and it is exactly the association predicted by a
source that requires **impact excavation / exposure of the deeper crust** (walls, rims, central
peaks of complex craters).

**Containment definition (state it explicitly — two defensible criteria give different fractions).**
The "inside a ≥10 km crater" figure depends on how containment is counted, and the two give 0.49 vs
0.75, so a reader must know which is used. This table uses the **nearest-crater** definition: a
hollow is inside only if the *nearest* crater center contains it within its rim radius (0.489). An
independent recomputation using an **any-crater** definition — inside if *any* ≥10 km crater contains
it — gives 0.753 (rotation null 0.454 ± 0.082, z = +3.7); that criterion is more permissive because
large, overlapping and basin-scale craters cover enormous area. The size cut is not the source of the
difference (the any-crater fraction is stable at 0.75/0.73/0.68/0.62 across ≥10/20/50/100 km cuts).
The nearest-crater criterion is the stricter and better-behaved one — a stricter criterion has a lower
null baseline, which is why it yields the larger z (+6.7) despite the smaller observed fraction. **The
association is strongly significant under either definition** (reconciliation recorded in
`runs/crater_association_reconciliation.json`); the conclusion is unaffected and doubly supported.

### 1.2 Pyroclastic vents — clean negative (inherited)

Hollow → nearest vent (n = 309): median 346.7 km. Near-field enrichment against the rotation null is
**not significant at any radius 25–200 km**; at 25 km the observed fraction (0.24%) sits *below* the
null mean (0.86%). A fumarolic/exhalative point-source origin predicts near-field concentration and
the data show its absence. **This argues against vent-fed volatiles as the general hollow-forming
mechanism** — and specifically against the fumarole-anchored concentration route (SCArFS) as the
*general* explanation, though it does not exclude it for individual vent-adjacent fields.

### 1.3 Concentrated-LRM raster — not demonstrable in this product (inherited)

Mean value of `concentrated_lrm_560_params_bw.tif` at hollows: observed 4.57 vs rotation null
6.82 ± 3.64 (z = −0.6, n.s.); fraction on mapped LRM observed 0.186 vs null 0.176 ± 0.051 (z = +0.2).
The apparent association (z = −4.2) under the naive null is an artifact of ignoring clustering (the
null SD inflates 7×). **This coarse (~0.12°/px ≈ 5 km) global "concentrated-LRM" product cannot
support an LRM-association claim.** This does **not** refute the published LRM association (Blewett
et al. 2011 [10.1126/science.1211681](https://doi.org/10.1126/science.1211681); Thomas et al. 2016
[10.1016/j.icarus.2016.05.036](https://doi.org/10.1016/j.icarus.2016.05.036)), which rests on
higher-resolution local reflectance mapping we do not currently hold. The crater association (§1.1)
is the geologically informative one at this resolution, and note that impact craters are precisely
where LRM is exposed (Klima et al. 2018 [10.1002/2018GL077544](https://doi.org/10.1002/2018GL077544)),
so the two associations are physically the same story sampled at different resolutions.

### 1.4 Reconciliation with the PI's own fits

The PI's `von_mises_fit.ipynb` fits **1-D longitude** distributions (bimodal von Mises) of hollow
counts/area against crater and vent longitude structure, finding hollow and crater longitude peaks
near the warm meridians (~98°/287°). That is a complementary, orbital-forcing-oriented analysis. The
tests here are **2-D nearest-neighbor + containment** tests with an explicit clustering-preserving
null. They agree on the crater association and add a quantified significance to it. **One gentle,
quantitative flag for the PI:** the von Mises count fits, and any association read from them, are
implicitly independence-based; the near-vanishing of the vent and LRM-raster signals here happens
*only* when clustering is preserved in the null. To be precise about scope: this affects the
**significance** attached to those two associations, **not the fits themselves**, and it does **not**
touch the crater association — that one is significant under both the independence and rotation nulls
*and* under both the nearest-crater and any-crater containment definitions (§1.1). Only the vent and
LRM-raster significances are sensitive to the null choice; worth revisiting those on the rotation
null.

### 1.5 The Wang et al. 2020 age constraint — and the LRM association resolved (new primary constraint)

Wang, Xiao, Chang & Cui 2020 (JGR Planets 125, e2020JE006559,
[10.1029/2020JE006559](https://doi.org/10.1029/2020JE006559), Crossref-verified this session) is a
**primary constraint** that changes the rate window and sharpens this track. From crater statistics
on hollow floors in <10 m/px MDIS NAC images they derive a **maximum model age for the global hollow
population of 103 (+200, −96) kyr** — an *upper limit* — and, from a mass-conservation budget over
all visible hollows, a **minimum lost-volatile volume ≥ 1,266 km³**. They conclude hollows are
actively forming and predict **ongoing rejuvenation of C- and S-bearing volatiles in the shallow
crust**, with hollows developed in **LRM that is richer in Mg/Ca/S/C than the global average and
exposed from depth** by impact excavation, mass wasting and cooling fractures.

Two consequences for this track:

1. **The LRM association is real; my null result was resolution-limited.** Wang's *local*
   high-resolution mapping places hollows in S/C-rich LRM excavated from depth — precisely the
   association my coarse (~0.12°/px ≈ 5 km) global "concentrated-LRM" raster could not resolve
   (§1.3). **Cite Wang for the LRM association, not my null.** My null result is not evidence against
   the association; it is evidence that *this particular global product cannot test it*. Because
   impact craters are where deep LRM is exposed (Klima et al. 2018), the strong crater association
   (§1.1) and Wang's LRM association are the same physical story at two resolutions.

2. **The concentration axis is sharpened decisively.** A ~100-kyr formation time with hollow depth
   set by **the thickness of the volatile-bearing layer** means the reservoir must be a
   **concentrated, recently exposed body** — a crustal volatile-bearing layer unroofed by impact, not
   a slowly-accumulated dispersed trace and not a fumarolic point source (consistent with the vent
   null, §1.2). This is exactly the concentration requirement the pore-space result imposes, now with
   an independent geologic timescale behind it.

Three caveats on how Wang's numbers are used here, per the PI's own reading of the full text: (i) the
103 kyr is explicitly an **upper limit** on the average model age, resting on a single probable impact
crater found across all <10 m/px NAC coverage of hollow floors — a thin statistical basis the authors
themselves flag; it is carried as an upper limit, not a measured age. (ii) Wang's area and depth come
from Thomas et al. 2014 (57,400 km²) and Blewett et al. 2016 (**24 m mean**); our own n=2518 point
catalog gives a **20.8 m median** — a mean-vs-median / different-subset difference. Both are stated
and never silently mixed: budget figures below use Wang's 24 m, the depth distribution elsewhere uses
our 20.8 m median. (iii) I do **not** claim "hollow depth does not correlate with model age" — that
is not sourceable to this paper; what Wang establishes is one global model age and that depth reflects
volatile-layer thickness.

The rate window this implies (required rate = depth / age; Blewett depths median 20.8 m, max 110.5 m;
Wang age bracket 7–303 kyr) is **6.9×10⁴ – 1.6×10⁷ m/Gyr**, about **3.2 orders of magnitude above the
top** of the Phillips 2021 Table-5 window (3.3–130 m/Gyr) used in Phases 1–2. **Both windows are
reported here so the sensitivity to this assumption is explicit; the switch is not silent.** The
loss-mechanism consequence — screened across the combined matrix — is stark:

| candidate | rate at 78.5N cold limit (m/Gyr) | vs Wang window (6.9e4–1.6e7) |
|---|---|---|
| Elemental S | 8.1×10¹¹ | **exceeds** by ~4.7 dex → supply-limited (the regime already identified) |
| C60 / C70 | 2.6×10³ / 3.1×10² | too slow by 1.4 / 2.3 dex — *was "thermophysically viable"* |
| NaCl / KCl (thermal) | 3.8×10⁻⁴ / 4.3×10⁻³ | too slow by 8.3 / 7.2 dex |
| Graphite | 5.5×10⁻⁶⁶ | too slow by ~70 dex |
| NaCl VUV photolysis (optimiztic yield) | 56 | **too slow by ~3.1 dex** |
| native Na PSD (best case) / cusp sputtering | — | too slow by ~0.8 / ~5.8 dex |

**Elemental sulfur is the only candidate whose rate exceeds the requirement**, and it exceeds it —
the supply-limited regime, exactly consistent with Wang's active-rejuvenation prediction. **The
athermal-chloride rehabilitation inverts**: under the young age it fails by ~3 orders even at the
most generous unmeasured quantum yield. That does not make the athermal work worthless — it converts
"chlorides are rehabilitated" into a stronger, independent negative: **no athermal route can form
hollows on the observed ~100-kyr timescale**, which eliminates an entire class of mechanisms. Alkali
chlorides therefore keep a 6/10 *geologic* score but their loss-mechanism role flips to ELIMINATED —
a clean illustration that geology and physics can disagree and that the loss-rate verdict, not the
geologic score, is decisive.

### 1.6 Wang argues for carbon — we engage it quantitatively and reject the *rate*, not the chemistry

Wang et al. 2020 (their §4.1) advance **carbon**, not sulfur, as the hollow-forming volatile: hollows
sit in graphite-darkened LRM; the graphite 500–600 nm feature is absent on hollow floors and present
in background terrain, so C appears lost from the floors; and they propose the loss route is **C +
implanted solar-wind protons → CH₄** (following Blewett et al. 2016), with escaping CH₄ lofting the
fine grains that make the bright haloes. This is in direct tension with our carbon verdict — but the
difference is **mechanism, not chemistry**: Phases 1–2 screened graphite on sublimation, photolysis
and sputtering (all ≤0.7 m/Gyr); Wang invokes a *reactive* route with a different rate law. So we
re-derived that route explicitly.

**The CH₄ route is falsified by conservation of protons — no free parameter — and two independent
routes agree to <0.5 dex.** The reaction C(s) + 4H → CH₄(g) cannot consume more carbon than there are
protons to react with, so the solar-wind proton flux is a *hard ceiling* on the rate, independent of
any assumed efficiency. I reproduced the PI's arithmetic from scratch and confirm it:

- **Per-area ceiling.** At the highest proton flux anywhere on Mercury — the northern cusp peak,
  ~1.03×10¹¹ m⁻²s⁻¹ near 70° N (Raines et al. 2022; Diego et al. 2020), which happens to coincide with
  the hollow latitudes — with every arriving proton forming CH₄ at 4 H per C, graphite removal is
  **7.5×10⁻³ m/Gyr** (PI: 7.2×10⁻³; the small difference is graphite density). Against the Wang window
  (6.9×10⁴–1.6×10⁷ m/Gyr) that is **short by 7.0–9.3 orders**.
- **Global-budget cross-check.** Removing 1,266 km³ of carbon as CH₄ needs 5.5×10⁴¹ protons (PI:
  5.7×10⁴¹); the solar wind delivers ~1.9×10³⁴ to the 57,400 km² hollow area in 103 kyr — **short by
  7.5 orders**. At cusp-peak flux and 100% efficiency the removal would take ~2,900 Gyr (PI: 3,062 Gyr;
  ~220× the age of the universe).

Two independent routes agreeing to within 0.5 dex. Rescue attempts fail: diluting the carbon helps
only linearly (at the measured 1.4 wt% the gap is still ~5.7 dex), and even a chemically nonsensical
1:1 H:C stoichiometry leaves ~5.1 dex. **The carbon elimination therefore no longer rests on a
modeled sublimation rate but on conservation of protons, which is far harder to contest.** Framed
precisely: Wang's *observations* stand (the graphite 500–600 nm feature absent on hollow floors,
present in background terrain; hollows concentrated in LRM); what fails is the *loss mechanism*
proposed to explain them. Their spectral evidence is equally consistent with graphite as the passive
host removed along with excavated material — exactly the Phase-1 assignment. The carbon is a
**passenger** of hollow formation, not its driver.

**The volume budget compounds the carbon problem and modestly favors sulfur.** Wang's own
mass-conservation budget: ≥1,266 km³ of volatiles lost, from a >22 m volatile layer over 57,400 km²
of mapped hollows, in ≤103 kyr. Against the near-surface inventory of the top 24 m (S at 4 wt% →
7,420 km³; C at 1.4 wt% → 3,070 km³):

| | lost / available |
|---|---|
| Sulfur (4 wt% terrane max, Weider et al. 2015) | **17.1%** |
| Sulfur (1.9-2.3 wt% global mean) | **29.7-35.9%** |
| Carbon (1.4 wt%, Peplowski 2015) | **41.2%** |

Consuming **41% of the planet's entire near-surface carbon inventory** in a geological instant, from
only the ~0.5% of LRM that hosts hollows, is a severe demand — and it lands on top of the 7.5-dex CH₄
rate shortfall. The sulfur demand (17%) is large but an order less extreme, and elemental S is the one
phase whose loss rate *exceeds* the requirement. Both classes require Wang's inferred **active
rejuvenation** of subsurface volatiles — which is independent support for the supply-limited /
replenished picture already reached from the depth distribution, and on the abundance/concentration
axes it favors S over C.

### 1.7 A concentration trade-off that binds every candidate — including sulfur

This constrains our own preferred answer, so it is the sharpest open question, not a point in
sulfur's favor — and it is a **trade-off, not a single threshold** (an earlier framing of it as a
"~20 vol% binding requirement" was a mis-framing: the arithmetic was right, the interpretation was
not, and it is corrected here). Wang's eq. 5 **couples** the volatile column to hollow depth: because
the halo is sub-millimeter-thin, h_volatile ≈ 0.95 × h_HF at every depth (verified independently:
h_volatile(24 m) = 22.9 m, reproducing Wang's published ">22 m"). The source volume-fraction is
therefore *f* = h_volatile / h_source, where **h_source — how far the source extends below the hollow
— is a free parameter**. The budget fixes only the *product* *f* × h_source, i.e. a trade-off between
how enriched the source is and how deep it reaches:

| source thickness (× hollow depth) | required *f*, median hollow | required *f*, deepest hollow |
|---|---|---|
| 1× | 94.7 vol% | 99.0 vol% |
| 2× | 47.4 vol% | 49.5 vol% |
| 5× | 18.9 vol% | 19.8 vol% |
| 10× | 9.5 vol% | 9.9 vol% |
| 20× | 4.7 vol% | 5.0 vol% |

**The single honest statement:** to source the volatile at Mercury's global-average sulfur content
(4 wt% ≈ 5.4 vol%), the reservoir must extend ~18× deeper than the hollow — about 365 m below a median
20.8 m hollow. So **either** the volatile is locally enriched well above the planetary average,
**or** the hollow taps a reservoir extending some hundreds of meters below it. The budget cannot
distinguish these; it only forbids both being false. These two branches are **co-equal hypotheses**,
not a requirement plus a caveat.

**The fracture-fed branch is favored by our own spatial results.** A deeper reservoir implies
fracture-mediated access — exactly the cooling-fracture structural control Wang document (their §3.1)
— and it is consistent with the significant crater association (§1.1: z = +4.0 at 10 km, z = +3.4 at
5 km) and with the *absence* of near-field vent association (§1.2), which together point at
impact-fracture access rather than fumarolic emplacement.

**This does not supersede the pore-space floor; both apply.** The pore-space threshold (x_v\* ≈ 4 vol%,
bracket 0–27%) remains a genuine, independent *floor* — below it, removal empties pores without
lowering the surface. The budget is a different kind of statement (a trade-off). The requirement is
that *f* exceed the pore-space floor **and** satisfy *f* × h_source = h_volatile. Note also that
Wang's ≥22 m is an explicit *minimum* (they cancel the lag term by setting ρ_lag = ρ_bedrock), so every
*f* in the table is a lower bound. This applies symmetrically to every candidate; sulfur passes the
proton-conservation test that kills the CH₄ route (it needs no external reactant) but faces this
concentration trade-off on equal terms.

---

## 2. Scoring rubric

Five axes, each 0/1/2 (higher = stronger geologic case for being the hollow-**forming** volatile):

| axis | 2 | 1 | 0 |
|---|---|---|---|
| **Source** | Hermean process documented | inferred from analogues/experiment | none/implausible |
| **Abundance** | MEASURED sufficient Hermean abundance | cosmochemical estimate, plausible | known-low/trace/absent |
| **Concentration** | satisfies the volume-budget trade-off (§1.7) by local enrichment **or** a documented deeper (fracture-fed) reservoir, as a **concentrated body whose loss lowers the surface** at 10–100 m, and clears the pore-space floor | mechanism exists but unproven/marginal, sits near the pore-space floor, or loss armours | no losable concentrated body (below pore-space floor, refractory-inert, or empties pores) |
| **Emplacement** | predicted association matches observation (crater z=+6.7; LRM lit.; equator-facing aspects) | partial/untested | contradicts observation |
| **Survival** | survives 4 Gyr gardening OR documented replenishment | marginal/needs unproven supply | destroyed |

**Concentration is the decisive axis**, and it is governed by two *distinct* statements (§1.7), both
of which must hold: the pore-space **floor** (x_v\* ≈ 4 vol% — below this, removal merely empties pores)
and the Wang 2020 volume-budget **trade-off** (source enrichment vs reservoir depth; *f* × h_source =
h_volatile). A phase that cannot be delivered as a *losable* concentrated body — either locally
enriched or from a fracture-fed deeper reservoir, and above the pore-space floor — fails as a former no
matter how good its physics. "†" marks any class whose case rests on a process **not documented on
Mercury**.

---

## 3. Scores by class

Full table with per-cell reasoning and DOIs: `runs/geologic_plausibility.csv`. Figure:
`geologic_plausibility.png`.

| class / entity | S | A | C | E | Su | /10 | role |
|---|:-:|:-:|:-:|:-:|:-:|:-:|---|
| Elemental sulfur (α-S) | 2 | 2 | 1 | 2 | 2 | **9†** | **FORMER (only survivor under Wang age)** |
| Crustal sulfides (CaS/MgS/MnS/FeS) | 2 | 2 | 1 | 2 | 2 | **9** | **RESERVOIR** |
| Graphite (elemental C) | 2 | 2 | 0 | 2 | 2 | 8 | HOST |
| Alkali chlorides (NaCl/KCl) | 1 | 1 | 1 | 1 | 2 | 6† | ELIMINATED under Wang age |
| Water ice | 2 | 2 | 2 | 0 | 0 | 6 | ELIMINATED (latitude/thermal) |
| Alkali sulfides (Na₂S/K₂S) | 1 | 0 | 1 | 1 | 1 | 4† | ELIMINATED (abundance) |
| Divalent chlorides (CaCl₂/MgCl₂/FeCl₂) | 1 | 1 | 0 | 1 | 1 | 4† | ELIMINATED (armour) |
| Fluorides (NaF/CaF₂/MgF₂) | 1 | 0 | 0 | 1 | 2 | 4† | ELIMINATED (both routes) |
| Phosphorus / phosphides | 1 | 0 | 0 | 1 | 2 | 4† | ELIMINATED |
| Other ices (CO₂/NH₃/SO₂/CO) | 1 | 1 | 2 | 0 | 0 | 4 | ELIMINATED (latitude/thermal) |
| Chalcogens (Se/Te) | 1 | 0 | 0 | 1 | 1 | 3† | ELIMINATED (abundance) |
| Nitrides (Si₃N₄/TiN/AlN) | 0 | 0 | 0 | 1 | 2 | 3† | ELIMINATED |
| Fullerenes (C₆₀/C₇₀) | 0 | 0 | 0 | 1 | 1 | 2† | ELIMINATED (source unobserved) |
| Organics (acids/PAHs/IOM) | 1 | 1 | 0 | 0 | 0 | 2 | ELIMINATED (radiolysis + location) |
| Native alkaline-earth metals (Mg/Ca) | 0 | 0 | 0 | 1 | 1 | 2† | ELIMINATED |
| Native alkali metals (Na/K) | 0 | 0 | 0 | 1 | 0 | 1† | ELIMINATED |

### The three verdicts that matter

1. **Elemental sulfur + crustal sulfides tie at 9/10 and are the same geologic system split by
   physics.** Sulfides are the highest-scoring class geologically — a documented magmatic source in
   a reduced, FeO-poor, S-rich crust (McCubbin 2012
   [10.1029/2012GL051711](https://doi.org/10.1029/2012GL051711); Vander Kaaden & McCubbin 2015
   [10.1002/2014je004733](https://doi.org/10.1002/2014je004733)), the highest MEASURED volatile
   abundance of any candidate (S up to ~4 wt%; Nittler et al. 2020
   [10.1016/j.icarus.2020.113716](https://doi.org/10.1016/j.icarus.2020.113716); Weider et al. 2014
   [10.1016/j.icarus.2014.03.002](https://doi.org/10.1016/j.icarus.2014.03.002)), a demonstrated
   occurrence control (Barraud et al. 2023
   [10.1126/sciadv.add6452](https://doi.org/10.1126/sciadv.add6452)), and the observed crater/LRM
   emplacement (deep-crust exposure). But **physics makes the sulfides the reservoir/marker, not the
   excavating phase**: they are 7–15 orders below the hollow-rate window and their incongruent loss
   armours. The phase that actually leaves is the **elemental sulfur** they release, cycled and
   replenished (Renggli et al. 2022 [10.1016/j.epsl.2022.117647](https://doi.org/10.1016/j.epsl.2022.117647);
   Verkercke et al. 2025 [10.1029/2025gl118112](https://doi.org/10.1029/2025gl118112)). **The
   reservoir↔former split is set by physics, not geology.** Elemental S carries a "†" on
   concentration/survival because a *static* surface reservoir is removed instantly — it is viable
   only as a **dynamically replenished** near-surface reservoir, and the specific concentration route
   remains the least-proven step (the vent test rules out the SCArFS fumarolic-redeposition version
   as the *general* mechanism).

2. **Graphite is the geologically ubiquitous phase that is physically inert — the LRM host, not a
   former, and Wang's CH₄ route is rebutted on rate (§1.6).** It scores 8/10 on source, abundance,
   emplacement and survival (flotation cumulate from a graphite-saturated magma ocean, Vander Kaaden &
   McCubbin 2015; measured <1 wt%, Xu et al. 2024
   [10.1038/s41550-023-02169-5](https://doi.org/10.1038/s41550-023-02169-5); the darkening agent of
   LRM, Klima et al. 2018) but **0 on concentration**: its only reactive loss channel, the C + H → CH₄
   route Wang et al. 2020 propose, caps at a proton-supply ceiling of 7.5×10⁻³ m/Gyr — 7.5 orders below
   Wang's own required rate — and the budget would demand 41% of the near-surface C inventory. So its
   loss cannot lower a surface on the observed timescale. It is the passive LRM host/residue whose
   depletion on hollow floors *marks* where hollows form without being the volatile that forms them.

3. **The clearest disagreements between physics and geology, stated plainly.**
   - **Fullerenes** are thermophysically among the *better* candidates (viable to 81–83° N, not
     athermally disqualified) yet score 2/10 — no Hermean process is documented that makes them at
     any abundance. *Physically fine, geologically impossible because the source is unobserved.*
   - **Water ice** has a real, concentrated, MESSENGER-confirmed reservoir (Paige et al. 1992
     [10.1126/science.258.5082.643](https://doi.org/10.1126/science.258.5082.643); Chabot et al. 2018
     [10.1017/9781316650684.014](https://doi.org/10.1017/9781316650684.014)) and scores 6/10 on
     source/abundance/concentration — but 0 on emplacement and survival: hollows are equator-facing,
     warm, sunlit, to 71–78° N, the **opposite** of permanently shadowed cold traps. *A perfect
     volatile in exactly the wrong place.*
   - **Organics** (incl. stearic acid, the 2021 "viable" organic) are radiolytically carbonised in
     125–251 Myr (central) before a hollow forms, and the free-organic reservoirs that do exist are
     polar cold-trap ones — again the wrong latitude.

---

## 4. Candidates whose case rests on an unobserved process (flagged)

Eleven of sixteen entries carry "†". The load-bearing ones for the final verdict:
- **Elemental S** — the *concentration* step (how a losable ~4% body is built and maintained) is the
  unproven link; the source, abundance and emplacement are all documented.
- **Alkali chlorides** — depend on **both** an unmeasured VUV quantum yield (~1e−3, ~50× the PSD
  anchor; Yakshinskiy & Madey 1999 [10.1038/23204](https://doi.org/10.1038/23204)) **and** an
  unobserved chloride-condensate concentration route (Llado et al. 2023
  [10.7185/gold2023.17838](https://doi.org/10.7185/gold2023.17838)) whose spatial prediction (vent
  proximity) is contradicted by §1.2.
- **Fullerenes, native metals** — source unobserved outright.

A conclusion that survives regardless of where these unmeasured quantities fall is a result; one
that depends on them is flagged as such throughout.

---

## 5. Bottom line for Phase 3a

The geology **confirms and sharpens** the physics ranking rather than overturning it, and the Wang
et al. 2020 young age (§1.5) sharpens it further. The sulfur–sulfide system is both the most
physically viable and the most geologically supported: a documented reduced-crust source, the
highest measured abundance, the one spatial association that survives a conservative null (impact
craters, z = +6.7), and — under the ~100-kyr formation timescale — **the only candidate whose loss
rate exceeds the required excavation rate** (supply-limited, matching Wang's active-rejuvenation
prediction). The Wang age simultaneously eliminates the athermal-chloride route by ~3 orders,
converting the earlier "chlorides rehabilitated" result into a stronger negative: no athermal route
can excavate 20 m in 100 kyr.

Wang's carbon hypothesis is met head-on: the C + H → CH₄ route is falsified by conservation of
protons — 7.0–7.5 orders short on two independent routes, with no free parameter (§1.6) — a mechanism
disagreement, not a chemistry one, and the budget would demand 41% of the near-surface carbon
inventory. Carbon is depleted from hollow floors as a passenger of formation, not its driver; the
budget's 17% sulfur demand and elemental S's rate excess both point the same way.

The open question is not *which element* — it is **how the sulfur is delivered as a concentrated
body**: either locally enriched well above the planetary average, or tapped from a reservoir extending
some hundreds of meters below the hollow through fractures (§1.7). These two branches are co-equal;
the volume budget forbids only both being false, and both sit above the independent pore-space floor.
The fracture-fed branch is favored by our own spatial results — the significant crater association
and the absence of near-field vent association together point at impact-fracture access rather than
fumarolic emplacement. That is the target the next-phase observational test must address: not
detecting sulfur *per se*, but distinguishing an actively-replenished, locally-enriched or
fracture-fed near-surface S reservoir in freshly exposed crater-wall LRM from a spent or dispersed
one.

**Hand-off to Phase 4 (observational test).** Wang themselves flag that the CH₄-vs-age tension is
testable with BepiColombo high-resolution imaging, and the carbon and sulfur hypotheses make
divergent, testable predictions: (i) a **sulfur** former predicts high S in fresh hollow-floor LRM
(SERENA/Strofio neutral S, MIXS S Kα), tracking the crater/deep-crust association; (ii) Wang's
**carbon** former predicts CH₄ (or its fragments) in the exosphere over active hollows and the
graphite 500–600 nm depletion — but must also survive the 7.5-dex rate objection, which it does not on
current fluxes. The budget's active-rejuvenation requirement is itself a target: a detectable,
present-day S resupply flux (consistent with Verkercke et al. 2025) is the signature that distinguishes
a replenished sulfur reservoir from a spent one.

**Deliverables:** `scripts/geologic_plausibility.py`, `runs/geologic_plausibility.csv`,
`runs/spatial_associations.json`, figure `geologic_plausibility.png`, this document.


> [CORRECTED 2026-09-02: 4 wt% is the high-sulfur TERRANE MAXIMUM, not the global mean (1.9-2.3 wt%, from S/Si = 0.076 XRS Nittler et al. 2020 / 0.092 GRS Evans et al. 2012 with Si ~23.4-27.1 wt%); and the attribution was wrong - Weider et al. 2016 (10.1002/2016gl068325) is a real paper but on sulfur- and carbon-driven explosive volcanism, NOT global mapping. The global-mapping source is Weider et al. 2015 EPSL (10.1016/j.epsl.2015.01.023). At the global mean the S inventory is 3,520-4,270 km3 and the lost fraction is 29.7-35.9%, not 17.1%, so the S-vs-C margin narrows from 2.4x to 1.1-1.4x: carbon remains the more extreme demand, but only just.]
