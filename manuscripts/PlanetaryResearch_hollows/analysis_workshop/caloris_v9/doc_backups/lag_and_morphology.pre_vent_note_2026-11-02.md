# Lag-deposit evolution, self-limitation, and hollow morphology (Phase 2c)

*Track: lag & morphology. PI: M. Phillips. Parent paper: Phillips, Moersch, Viviano &
Emery (2021), Icarus 359, 114306, doi:10.1016/j.icarus.2021.114306.*

A phase that sublimates fast at a bare surface but armours itself after a few microns cannot
excavate a 20 m hollow. Phase 1 screened **surface** rates; this track asks whether each
surviving candidate can **keep going** (self-limitation), and whether what it would produce
actually **looks like a hollow** (morphology). All numbers produced here are **MODEL OUTPUT**;
measured inputs are labeled at point of use.

Surviving-candidate set carried forward from Phase 1: elemental sulfur (S8, reference/SCArFS),
fullerene C60, alkali chlorides NaCl/KCl, divalent chlorides CaCl2/MgCl2/FeCl2, and the fast-enough
alkali sulfides Na2S/K2S (eliminated on abundance, retained here to test the lag physics).

---

## 0. Validation gate (reproduce the paper before trusting new numbers)

The lag machinery (`scripts/sublimation_halite.py`, Phillips 2021 Eqs. 14–15, Knudsen diffusion
through a growing lag) has a single free parameter — the effective pore radius — pinned to sulfur's
**published** hot-limit lag. On the paper's own vapor-pressure convention:

| quantity | this work | published (Phillips 2021) |
|---|---|---|
| S hot-limit lag (0°N, threshold 6.4 m/Gyr) | **15.00 m** | 15 m |
| implied original-layer purity | **63.1%** | ~62–63% |
| S cold-limit purity (70°N, warm meridian) | **99.34%** | ~99.5% |
| calibrated pore radius | **0.464 µm** | (implicit) |

Both gates pass. The calibrated pore radius (0.46 µm) is physically reasonable for a fine
regolith and is consistent with the archive grain-spacing range (40–130 µm inter-grain, µm-scale
throats). We use **this physical pore radius** for the candidate screen.

**A convention finding (the "we can do better now" number).** The paper's finite 15 m self-limiting
lag for sulfur holds on the tabulated (kPa-as-Pa) vapor-pressure convention. On the project's
**settled SI-Pa convention** (A_Pa = A + ln 1000; established Phase 1), sulfur's true vapor pressure
is ~1000× higher, and reproducing a 15 m lag would demand a sub-nm (~0.46 nm) pore radius —
unphysical. With true vapor pressures **and** a physical pore radius, a sublimator as volatile as
sulfur at the hot equator **cannot be throttled by any diffusive lag**: it is **supply-limited**, not
self-limited. This does not contradict the paper — it refines it. It is exactly the SCArFS picture
(sublimation cycling around fumaroles requires **replenishment**; a static reservoir is removed at
1.9×10¹¹ m/Gyr, Phase 1). We report on the SI-Pa convention and flag the difference, per the PI's
standing instruction.

---

## 1. Self-limitation: three modes, and which candidates stall

The operative self-limitation mechanism depends on the loss route:

- **Congruent thermal sublimators** (S8, C60): no residue; the only brake is a Knudsen-diffusion
  lag of *inert host* regolith (paper Eqs. 14–15).
- **Athermal photolytic / incongruent phases**: loss removes one sublattice (Cl by VUV photolysis;
  S2 by decomposition) and leaves the **counter-species metal** behind. Self-limitation is then a
  **chemical-armouring** question: does the metal residue accumulate into a coherent barrier, or is
  it removed as fast as it is exposed?

The primary switch is **whether the residue is removed** (not its instantaneous volume fraction).
A refractory metal residue that is *not* removed accumulates and armours; for VUV photolysis a
continuous metal film is **VUV-opaque** (metal skin depth ~5–15 nm at 100–150 nm), so it blocks
the photolysis of the chloride beneath it. The parent thickness consumed before that opaque film
percolates sets the **depth-to-stall**.

| candidate | formula | route | self-limit mode | stalls? | depth-to-stall | residue (mass%) | lowers surface |
|---|---|---|---|---|---|---|---|
| orthorhombic sulfur | S8 | thermal sublimation | congruent no residue | no | — | — | True |
| fullerene C60 | C60 | thermal sublimation | congruent no residue | no | — | — | True |
| halite | NaCl | vuv photolysis | self sustaining | no | — | 39% | True |
| sylvite | KCl | vuv photolysis | self sustaining | no | — | 52% | True |
| CaCl2 | CaCl2 | vuv photolysis | armouring | YES | 11 nm | 36% | conditional |
| MgCl2 | MgCl2 | vuv photolysis | armouring | YES | 16 nm | 26% | conditional |
| FeCl2 | FeCl2 | vuv photolysis | armouring | YES | 31 nm | 44% | conditional |
| Na2S | Na2S | incongruent | self sustaining | no | — | 59% | True |
| K2S | K2S | incongruent | self sustaining | no | — | 71% | True |

**Result.** Three groups emerge:

1. **Armouring → STALLS (divalent chlorides CaCl2, MgCl2, FeCl2).** The refractory Ca/Mg/Fe metal
   residue is not PSD-mobile and accumulates into a VUV-opaque film after only
   **11, 16, and 31 nm** of parent consumed
   (CaCl2/MgCl2/FeCl2). They cannot excavate more than a monolayer-scale film — **eliminated as
   bulk hollow-formers by self-limitation alone**, independent of their (already conditional)
   Phase-1 yield brackets. FeCl2's huge nominal photolysis rate (narrow 2.5 eV gap) is
   **self-defeating**: it plates out its own opaque Fe barrier fastest of the three.

2. **Self-sustaining residue (alkali chlorides NaCl/KCl; alkali sulfides Na2S/K2S).** The Na/K
   counter-species is PSD-released ~10³–10⁵× faster than the anion is liberated (Phase 1), so no
   coherent barrier accumulates. These do **not** stall on residue grounds. (NaCl/KCl remain
   yield-limited from Phase 1: viable only at quantum yield ≥ 10⁻³; the sulfides remain eliminated
   on **abundance**, not self-limitation.)

3. **Congruent, no residue (S8, C60).** No armouring possible.
   - **C60** self-limits diffusively with a thin lag (0.26 m at 0°N, 0.024 m at 70°N) but demands a
     near-pure source: **99.0% purity at the hot limit and 99.8% at the cold limit.**
     Consistent with Phase 1's source-limited verdict — a >99% pure C60 lens has no known Hermean
     source.
   - **S8** cannot self-limit diffusively at all (§0): supply-limited, the SCArFS regime.

---

## 2. Morphology consistency

### 2.1 Depth distribution: no self-limiting ceiling

A self-limiting phase predicts a **hard depth ceiling** — hollows should pile up at a maximum and
the distribution should truncate. The measured Blewett distribution (n=2518) shows the opposite:

- percentiles (MEASURED, m): 5th 5.08, 25th 12.21, **median 20.8**, 75th 32.1, 95th 53.49, max 110.47
- the excess-over-median tail is **consistent with a memoryless exponential** (scale 15.0 m; KS
  D=0.038, p=0.05) — the signature of a **supply/age-limited** process, not a bounded one
- max/95th ratio = **2.07** (a ceiling-limited distribution would be ~1.0–1.1); only 2 of 2518
  hollows lie within 5 m of the maximum — **no pileup**

**This conclusion is resolution-robust.** The depth measurements carry a strong image-resolution
bias (§4.3: Spearman(pixel scale, depth) = +0.52, coarse images inflate depths and preferentially
populate the deep tail), so the tail fit was re-run *within* each pixel-scale stratum. The memoryless
exponential tail and the absence of a ceiling **survive in every stratum** (px 1–8/8–12/12–16/16–20 m:
KS p = 0.86/0.09/0.51/0.12, none rejecting exponential; max/95th = 2.44/2.16/1.75/1.77, no pileup).
The overall median does inflate with pixel scale (11.9 → 29.5 m across the strata), confirming the
bias, but the tail *shape* — the load-bearing feature — does not depend on it.

**Implication.** The observed morphology is inconsistent with a phase that self-limits at a fixed
depth (the divalent chlorides, had they survived §1; or a shallow-lag former). It **is** consistent
with a supply-limited former whose depth is set by how much volatile was locally available and how
long it has been eroding — i.e. **sulfur in the SCArFS/replenished regime**, or a source-limited lens
former where the lens thickness varies. This demotes any residual self-limiting candidate and
**favors the supply-limited sulfur picture**.

### 2.2 Slope aspect: a discriminant already in the archive

Hollows are systematically **equator-facing** (archive n=31: **24/31 within 45° of due south**,
median aspect 187.5° — MEASURED). Modeled diurnally-integrated forcing on a 25° slope at 40°N
peaks at exactly **180° (due south)**, with a south/north ratio of:

| mechanism | scaling | south/north forcing ratio | predicts equator-facing? |
|---|---|---|---|
| thermal sublimation | ∝ exp(B/T) | **≤ 19×** (upper bound; T_S=702 K vs T_N=571 K) | strongly YES |
| VUV photolysis / PSD | ∝ photon dose (linear) | **2.3×** | weakly yes |
| cusp ion sputtering | magnetospheric geometry | **1.0× (flat)** | NO |

Under a formation **threshold** (a hollow forms only where the local rate clears the excavation
window), the observed 24/31 = 0.77 equator-facing fraction is reproduced by a **thermal** former at a
modest threshold (≈0.5–0.7 of the maximum rate) but requires an **extreme** threshold (>0.9) for a
linear-in-dose photolytic former:

| threshold (fraction of max rate) | thermal frac-within-45°-S | photolytic frac-within-45°-S |
|---|---|---|
| 0.1 | 0.35 | 0.26 |
| 0.3 | 0.51 | 0.26 |
| 0.5 | 0.66 | 0.33 |
| 0.7 | 0.90 | 0.49 |
| 0.9 | 1.00 | 0.90 |

**Discriminant verdict.** The equator-facing aspect preference **decisively excludes an
aspect-neutral cusp-sputtering former** (predicts a flat aspect distribution; observed strongly
peaked). It **favors a thermally-driven former** (exponential temperature amplification reproduces
the observed concentration at a modest threshold) over a purely photolytic one, though VUV
photolysis is not excluded because photon dose tracks the same equator-facing geometry (weaker, but
same sign). This is a **potential discriminant sitting in data the PI already has** (RJThomas aspect
table, n=31); the strongest test would be to enlarge that aspect sample and measure the
concentration quantitatively against the threshold curves above.

**The aspect discriminant is free of the depth-measurement resolution bias** (§4.3). Aspect is a
slope *azimuth*, not a depth, so it does not carry the pixel-scale confound structurally, and it does
not: matching the 31 aspect hollows to the catalog by location, the recovered pixel scale is
uncorrelated with the deviation from due south (Spearman = 0.11, p = 0.55), and the equator-facing
preference holds equally in the fine- and coarse-resolution halves (12/15 and 12/16 within 45° of S).
All 31 aspect hollows lie at 4–60°N, equatorward of the ≥70° censoring regime, so the detection-limit
censoring that afflicts the high-latitude depths does not touch this sample either. The aspect
verdict therefore stands unchanged after the confound is controlled.

### 2.3 Pore-space objection: surface lowering requires the volatile to be a major solid constituent

The referee's objection (quoted): *"if S is added in pore spaces, then how does its removal lead to
lowering of the ground surface? Why doesn't this just result in empty pore space?"* — is correct for
a **dispersed** volatile, and we answer it quantitatively.

Model the bulk as (refractory host grains) + (volatile) + (pores), total solid volume fraction
φ_total. Let x_v be the volatile's fraction **of the solids**. After the volatile is removed, the
refractory framework occupies φ_total·(1−x_v) of the bulk. That framework **self-supports** (holds
the surface up → empty pores, the referee's outcome) if it exceeds the granular rigidity /
loose-jamming percolation threshold φ_c (measured granular physics: high-friction loose jamming
φ_c≈0.44, loose random packing ≈0.55, up to ≈0.60). The surface **collapses** (lowers) only if

  φ_total·(1 − x_v) < φ_c  ⇔  **x_v > x_v\* = 1 − φ_c/φ_total**

| φ_c (rigidity threshold) | φ_total = 0.50 | φ_total = 0.60 |
|---|---|---|
| 0.44 | x_v* = 0.12 | x_v* = 0.27 |
| 0.55 | x_v* = 0.00 | x_v* = 0.08 |
| 0.60 | x_v* = 0.00 | x_v* = 0.00 |

**Answer.** The crossover spans **x_v\* ≈ 0–27% of the solid volume** across the φ_c=0.44–0.60,
φ_total=0.50–0.60 bracket (six computed values 0, 0, 0, 8%, 12%, 27%; median ~4%). The upper end
(~27%) applies to a densely-packed, high-rigidity regolith (φ_total=0.60, φ_c=0.44); at the lower
end (loose packing already near the rigidity threshold) even a **few percent** volatile suffices, and
where φ_total ≤ φ_c the framework never self-supports (x_v\*=0). Below x_v\*, removing the volatile
leaves empty pores and the surface does **not** lower — the referee is right. Above it, the residual
framework cannot self-support, compacts, and the surface lowers. Surface lowering therefore requires
the hollow-forming volatile to be either a **non-trivial solid constituent** (a few percent up to
~27% of the solids, depending on packing) or a **near-pure buried lens** (x_v→1, full collapse,
surface lowers by ~the lens thickness). This is a real constraint on **every** candidate:

- **S8 as a near-pure lens / SCArFS-deposited layer**: x_v→1, surface lowers fully. ✔
- **S8 dispersed at trace level in LRM pore space**: x_v ≪ x_v*, surface would **not** lower —
  removal just empties pores. So the SCArFS/replenished picture, which concentrates S into
  fumarole-fed layers, is **required** not just for the rate (§0) but also for surface lowering.
- **Divalent chlorides**: moot — they stall at nm scale (§1) before removing enough to matter.
- **C60**: needs >99% purity anyway (§1.3); a pure lens lowers the surface, but no such source exists.

### 2.4 The Wang budget gives an enrichment-or-depth TRADE-OFF (not a ~20 vol% floor), on top of the pore-space floor

Wang 2020's global volatile budget adds a second geologic constraint that binds **every** candidate
including sulfur — but it is a **trade-off, not a single number**. (An earlier version of this
section mis-stated it as a "~20 vol% floor"; the arithmetic was right, the framing was wrong, and it
is corrected here.) Wang eq 5 gives the integrated volatile column removed,

  h_volatile = h_inner + h_HF − (ρ_halo/ρ_bed)·(h_inner + h_outer·A_bnd/A),

and with their own parameters (ρ_halo=1100, ρ_bed=2800 kg/m³, h_inner=75 µm, h_outer=2 m,
A_bnd/A=1.40) the halo terms are sub-millimeter, so **h_volatile ≈ 0.95·h_HF at every depth**
(verified: h_volatile(24 m) = 22.9 m, reproducing Wang's published ">22 m"). The volatile *column*
is therefore set by the hollow depth. The source volatile fraction is *f* = h_volatile / h_source,
and **h_source — how far the source extends below the hollow — is a free parameter**: the budget
fixes only the *product* f·h_source, i.e. a trade-off between enrichment and reservoir depth:

| source extends below hollow | required *f* (median 20.8 m) | required *f* (deepest 110.5 m) |
|---|---|---|
| 1× hollow depth | 94.7 vol% | 99.0 vol% |
| 2× | 47.4 vol% | 49.5 vol% |
| 5× | 18.9 vol% | 19.8 vol% |
| 10× | 9.5 vol% | 9.9 vol% |
| 20× | 4.7 vol% | 5.0 vol% |

The "~20 vol%" figure is simply the 5×-depth row — one geometry among many, not a floor. **The honest
statement:** to source the volatile at Mercury's global-average sulfur content (4 wt% = 5.4 vol%),
the reservoir must extend **~18× deeper than the hollow — about 365 m below a median 20.8 m hollow**.
So **either** the volatile is locally enriched well above the planetary average, **or** the hollow
taps a reservoir extending some hundreds of meters below it. The budget cannot distinguish these; it
only forbids **both** being false — a dispersed trace volatile in a shallow exposed layer.

These two branches — **local enrichment** vs **fracture-fed depth** — are **co-equal hypotheses**,
not a requirement plus a caveat. Local enrichment implies a concentrating process (magmatic sulfide
segregation, exhalative deposition); a deep reservoir implies fracture-mediated access — exactly the
cooling-fracture structural control Wang document (their §3.1), and consistent with the significant
crater association the geologic track found (5 km z=3.4, 10 km z=4.0) and the *absence* of near-field
vent association, which together point at impact-fracture access rather than fumarolic emplacement.

**Both constraints apply; the budget does not supersede the pore-space floor.** They are different
kinds of statement: the pore-space threshold (x_v\* ≈ 4 vol%, bracket 0–27%; §2.3) is a genuine
independent *floor* — below it, removal empties pores without lowering the surface — while the budget
is the trade-off above. The requirement is *f* > pore-space floor **and** f·h_source = h_volatile.
The net geologic conclusion is unchanged and robust: the former must be **concentrated rather than
dispersed trace material**, now resting on the pore-space floor plus the enrichment-or-depth
trade-off rather than on a single 20 vol% figure. This is the sharpest open question for sulfur and
is handed to the geologic track.

### 2.5 Independent check of the carbon/CH₄ falsification (proton-supply ceiling)

The PI's Result 1 falsifies Wang's proposed carbon loss route (C + solar-wind H → CH₄) on proton
conservation, and asked for an independent arithmetic check. **Confirmed to 0.5 dex.** The reaction
cannot consume more carbon than there are protons to react with, so the proton flux is a hard
ceiling with no free parameters. At the **cusp-peak** proton flux (1.03×10¹¹ m⁻² s⁻¹, the highest
anywhere on Mercury) with 4 H per C, graphite removal ≤ **7.4×10⁻³ m/Gyr** — short of the Wang window
floor (6.9×10⁴) by **7.0 dex**. Independent global-budget cross-check: removing Wang's 1266 km³ of C
as CH₄ needs 5.6×10⁴¹ protons, while even cusp-peak flux over the 57 400 km² hollow area in 103 kyr
delivers 1.9×10³⁴ — short by **7.5 dex**. The two routes agree to 0.5 dex, and the falsification is
*conservative* because it assumes the maximum proton supply on the planet (disk-averaged solar wind,
~3×10⁹ m⁻² s⁻¹, is ~1.5 dex worse). This supersedes any modeled-sublimation-rate argument for
carbon: the elimination now rests on conservation of protons, which is far harder to contest.
Wang's *observations* stand (the graphite absorption feature is absent from hollow floors and present
in background terrain; hollows concentrate in LRM); what fails is the *loss mechanism* proposed to
explain them — and those observations are equally consistent with graphite as the **passive host**
removed along with excavated material, exactly the Phase-1 assignment. (Rate law and full numbers:
`concentration_requirement.json`; helper `proton_supply_ceiling()` in `lag_morphology.py`.)

---

## 3. Synthesis — what survives the lag & morphology filter

| candidate | §1 self-limitation | §2.1 depth ceiling | §2.2 aspect | §2.4 enrichment-or-depth trade-off | verdict |
|---|---|---|---|---|---|
| **elemental S (SCArFS/replenished)** | supply-limited (no diffusive stall) | matches heavy tail | matches (thermal) | local enrichment OR ~365 m fracture-fed reservoir | **survives — favored (concentration is the open question)** |
| C60 | self-limits, needs >99% purity | matches | matches (thermal) | pure lens satisfies it, but no source | survives thermophysically; **source-limited; too slow vs Wang window (§5)** |
| NaCl / KCl | self-sustaining residue | matches | weaker (photolytic) | major constituent | survives §2c filter; **eliminated by Wang window (§5) — athermal too slow** |
| Na2S / K2S | self-sustaining residue | — | — | — | survive §2c; **eliminated on abundance (Phase 1) + Wang window** |
| **CaCl2 / MgCl2 / FeCl2** | **ARMOUR → stall at 11/16/31 nm** | would predict ceiling (not seen) | — | fail (stall first) | **ELIMINATED by self-limitation** |
| graphite (Wang's CH₄ route) | — | — | — | — | **CH₄ loss route FALSIFIED on proton supply (§2.5); passive host** |

**Headline results of this track:**

1. **The divalent chlorides are eliminated by self-limitation** (stall at 11–31 nm), a filter that
   pure surface-rate arguments miss. This sharpens Phase 1, where they were only *conditionally*
   rehabilitated.
2. **The depth distribution has no self-limiting ceiling** (memoryless exponential tail to 110 m),
   so the hollow-former is **supply/age-limited**, not depth-self-limited — favoring the
   replenished-sulfur (SCArFS) picture over any fixed-lag former.
3. **The equator-facing aspect preference** (24/31 within 45° of S) **favors a thermal former and
   excludes an aspect-neutral cusp-sputtering former** — a discriminant already in the PI's archive.
4. **The hollow-former must be concentrated, not dispersed trace material.** Two constraints of
   different kinds, **both** applying (the budget does *not* supersede the floor): (a) the pore-space
   floor for surface lowering to occur at all, x_v\* ≈ 4 vol% median (§2.3); (b) the Wang budget
   **enrichment-or-depth trade-off** (§2.4), f·h_source = h_volatile — *either* local enrichment above
   the planetary average *or* a reservoir extending ~18× deeper than the hollow (~365 m below a median
   hollow to source at global-average S). These two branches are **co-equal hypotheses**; the budget
   forbids only both being false (a dispersed trace volatile in a shallow exposed layer). The
   fracture-fed branch aligns with the geologic track's significant crater association and the absent
   vent association. This is the sharpest open question for sulfur, handed to the geologic track.
   **Wang's carbon/CH₄ loss route is falsified on proton supply** (§2.5, independently confirmed to
   0.5 dex: graphite removal ≤ 7.4×10⁻³ m/Gyr vs a 6.9×10⁴ requirement, short by 7 orders) — carbon
   is the passive LRM host, not the excavated volatile, as Phase 1 concluded.
5. **The poleward limit is now a parameter** (71.4 → 78.53°N; §4). Every thermally-limited candidate
   fails harder at the new limit (halite's thermal exclusion widens 4.4° → 11.5°). Against the
   *Phillips 2021* window VUV photolysis still reaches at 78.53°N (Y ≥ 3.2×10⁻⁴); against the *Wang*
   window (§5) it does not. The depth-vs-latitude discriminant (§4.3), once the image-resolution
   confound is controlled, **excludes cusp sputtering by the wrong decline sign** (robust); its
   magnitude is provisional and, in the supply-limited case, not an age diagnostic.
6. **The Wang 2020 age (~103 kyr) raises the required rate window by ~3 orders and inverts the
   verdicts** (§5). Against 6.9×10⁴–1.6×10⁷ m/Gyr, **only elemental sulfur survives** — with the
   supply-limited sign this track already identified — while C60/C70 flip to too-slow and **every
   athermal route fails by 0.8–8.3 dex**. Wang's independent finding that depth is set by
   volatile-layer thickness (not duration) is direct support for the §2.1 no-ceiling/supply-limited
   conclusion. Both rate windows are reported; the earlier 189 Myr–5.5 Gyr excavation bracket is
   withdrawn as the operative timescale.

All four independently point the same way: **elemental sulfur in a replenished/concentrated
(SCArFS) regime** is the candidate most consistent with lag physics and hollow morphology; the
divalent chlorides are removed; C60 and the alkali chlorides/sulfides survive the morphology filter
but remain limited by source abundance or quantum yield as established in Phase 1.

---

## 4. The poleward limit is a PARAMETER (71.4 → 78.53°N), and it sharpens every latitude test

The 2021 paper's 71.4°N cold limit came from the Thomas 2014 database. The PI re-derived it from
the more complete Blewett 2016 catalog (n=2518, lat/lon/depth supplied locally) — **78.53°N**
(|lat| 99.9th pct; 99.0th = 73.05°, 99.5th = 77.78°), a **+7.13° poleward shift** — and will swap in
the ML HORNET catalog (Bickel et al. 2025, 19,110 hollows, PI-supplied doi:10.48620/77046 — a
University of Bern BORIS *data* deposit indexed on DataCite, not Crossref-resolvable; unverified
here, pending its return from maintenance) when its deposit returns. The catalog is overwhelmingly northern (**2462 N vs 56 S**, southern extreme only −3.44°);
this is an MDIS coverage/resolution bias, **not** a real hemispheric asymmetry. The limit is now a
parameter throughout (`POLEWARD_LIMITS` in `lag_morphology.py`); every latitude-sensitive result
below is reported at **both 71.4 and 78.53°N** so HORNET drops in without a re-run.

### 4.1 Thermally-limited candidates fail harder at 78.53°N

| candidate | route | surf rate 71.4°N | reaches window? | surf rate 78.53°N | reaches window? |
|---|---|---|---|---|---|
| elemental S | thermal | 6.8×10¹² | yes (supply-limited) | 7.7×10¹¹ | yes (supply-limited) |
| C60 | thermal | 5.7×10⁵ | yes (needs ~100% purity) | 2.5×10³ | yes (needs ~100% purity) |
| Na2S | incongruent | 4.4×10⁻² | **no** | 2.6×10⁻⁵ | **no (harder)** |
| K2S | incongruent | 2.2×10⁻¹ | **no** | 1.6×10⁻⁴ | **no (harder)** |

(rates m/Gyr; window floor 3.3 m/Gyr). The alkali sulfides were already below the window at 71.4°N
on rate; at 78.53°N they fall a further 3–4 orders — so their Phase-1 elimination (on abundance)
is reinforced by rate too. **Halite** (the athermal paradigm): its thermal cold limit is ~67°N, so
it was thermally excluded by **4.4°** at 71.4°N and is now excluded by **11.5°** at 78.53°N — the
thermal exclusion **strengthens**, exactly as the PI noted, and the same holds for every
thermally-limited refractory.

### 4.2 Athermal route: photolysis still reaches the window at 78.53°N — a STRONGER argument

The poleward shift widens the gap the athermal route must cover, because photolysis falls only
~cos(lat) (**×5** equator→78.5°N) while thermal sublimation collapses (**×662**). Quantitatively,
thermal is suppressed relative to photolysis by **22×** at 71.4°N but **132×** at 78.53°N — a **6×**
widening of the thermal/athermal gap. Crucially, **VUV photolysis of NaCl still reaches the hollow
window at 78.53°N**: at the mid yield (Y=10⁻³) the rate is **10.5 m/Gyr** (window floor 3.3), and the
threshold yield to reach the window is **Y ≥ 3.2×10⁻⁴** (17× the PSD anchor, still below mid). So a
more complete catalog puts hollows *further* poleward of where thermal sublimation of any refractory
can operate, yet leaves the photolytic route intact — a **stronger** case for the athermal route
than we had at 71.4°N, because the thermal alternative is now even more comprehensively excluded.

### 4.3 A morphological discriminant: the poleward depth decline — with the image-resolution confound controlled

The poleward tail is systematically **shallow** (MEASURED, Blewett 2016, |lat| convention): the 8
highest-latitude hollows are all <10 m against a 20.8 m global median, and band medians fall from
~20–26 m through 0–60° to 11.8 m (60–65°), 9.9 m (65–70°), 6.0 m (75–80°).

**But the depth measurements carry a strong image-resolution bias** (PI-flagged confound). The
catalog's per-hollow pixel scale correlates with measured depth at **Spearman = +0.52 (p = 4×10⁻¹⁷⁸)**
— coarser images give systematically deeper hollows, a bias **3× stronger** than the raw
depth–latitude trend itself (−0.16). There is also **censoring**: poleward of 70°, **55.6%** of
cataloged hollows are shallower than one pixel (vs 13–17% equatorward), so high-latitude depths sit
at or past the detection limit. This must be controlled before the trend can carry any mechanism
weight. What survives:

- **The trend (sign) is ROBUST.** Partial Spearman(|lat|, depth | pixel scale) = **−0.157
  (p = 3×10⁻¹⁵)**, essentially unchanged from the raw −0.16, and it holds independently *within* every
  resolution stratum (px 1–8/8–12/12–16/16–20 m: ρ = −0.130/−0.161/−0.202/−0.270, all p<0.01).
- **The MAGNITUDE is NOT robust.** The equator→pole decline factor is cut- and control-dependent
  over **~1.9–4.0×** (raw ≥75° 3.98× on n=21; raw ≥70° 2.48×; resolution-matched 2.34×;
  stratum-weighted 1.92×). The precise **3.98× is withdrawn** — it was an artifact of a 21-hollow bin.
  (A censoring-corrected variant with a uniform detectability floor leaves only 5 high-latitude
  hollows and is not usable.)

The discriminant compares each mechanism's predicted equator→pole **rate** decline against this
observed **depth** decline (valid, under depth ∝ cumulative rate, only for a *near-threshold /
rate-limited* former):

| mechanism | predicted equator→pole rate decline | vs observed 1.9–4.0× | verdict |
|---|---|---|---|
| cusp ion sputtering | **0.04×** (i.e. ~25× **UP** poleward, flux peaks at 70°) | wrong sign | **EXCLUDED (robust)** — predicts *deeper* poleward; data show the opposite |
| VUV photolysis | **4.46×** down | at/above the top of the range | **only prediction in the right regime** |
| thermal sublimation, *rate-limited* | **615×** down | ~2 orders too steep | a rate-limited thermal former is **excluded by magnitude** |

**Corrected reading (this supersedes the earlier "only thermal tracks the decline" claim).** The
robust part is the **sign**: the decline excludes cusp sputtering outright (it predicts the opposite
sign) and this is the sharpest single discriminant in the morphology data. On **magnitude**, the
only prediction in the observed regime is **VUV photolysis (4.46×)**, sitting at or just above the
top of the 1.9–4.0× range; a *rate-limited* thermal former (615×) is ~2 orders too steep.

**This does not demote sulfur, because sulfur is not rate-limited.** Its sublimation rate exceeds the
hollow window by ~11 orders at *every* hollow latitude (7.7×10¹¹ m/Gyr even at 78.5°; §4.1), so its
depth is set by reservoir **supply**, not by the rate collapse — the depth ∝ rate mapping, and hence
the 615× magnitude test, simply does not apply to it. For sulfur only the decline *sign* is
constraining, and it is satisfied. The magnitude test is diagnostic for a **near-threshold** former,
i.e. the marginal photolytic chloride, and there photolysis is the mechanism that lands in the
observed range.

**Reconciling §4.2 and §4.3.** No tension: (a) cusp sputtering is excluded by the decline sign; (b)
the *existence* of hollows out to 78.5°, beyond any refractory's thermal reach, needs an athermal
contributor for thermally-excluded phases, and photolysis both reaches the window there (§4.2) and
matches the depth-decline magnitude for a near-threshold former; (c) elemental sulfur, being
supply-limited, is unconstrained by the magnitude and remains consistent via the sign. The
depth-latitude *slope* is **provisional** pending the HORNET catalog (Bickel 2025, 19,110 hollows),
whose documented resolution field will let the completeness/censoring correction be modeled rather
than stratified around.

---

## 5. The Wang 2020 age constraint raises the rate window by ~3 orders — and only sulfur survives it

**New primary constraint** (PI): Wang, Xiao, Chang & Cui 2020 (*JGR Planets* 125, e2020JE006559,
doi:10.1029/2020JE006559 — **Crossref-verified**: authors Wang/Xiao/Chang/Cui, JGR Planets, 2020).
Crater statistics on hollow floors in <10 m/px NAC images give a **maximum model age of 103 (+200,
−96) kyr** for the global hollow population. This is explicitly an **upper limit on the average model
age**, derived from a **single** probable impact crater found in *all* <10 m/px NAC coverage of
hollow floors (summed area >18 km²) — a deliberately thin statistical basis the authors flag as
such. It is carried here as an upper-limit estimate, not a measured age; the qualitative conclusions
below hold across the whole 7–303 kyr bracket, not just the best-fit value. Two findings bear on this
track:

- **Hollow depth is set by the thickness of the volatile-bearing LAYER**, not by how long
  sublimation has run (Wang's implied vertical growth rate is >0.23 mm/yr = 2.3×10⁵ m/Gyr, >1000×
  prior stratigraphic-age estimates). *(The stronger phrasing "depth does not correlate with model
  age" is not sourceable to this paper — they establish one global model age from one candidate
  crater, so there is no depth-vs-age regression to cite; only the layer-thickness statement is
  theirs.)*
- Hollows are concentrated in **impact-crater settings associated with LRM excavated from depth**,
  consistent with a crustal reservoir exposed by impacts.

**The rate window moves by ~3 orders.** Required rate = depth / age. Two depth statistics are in play
and must not be mixed: Wang's own budget uses the **24 m average** depth (Blewett 2016) over the
Thomas 2014 hollow area; this track's depth-distribution work uses the **20.8 m median** from the
Blewett n=2518 point catalog. They differ (mean vs median, different subsets) but agree to ~15%; the
window below uses the 20.8 m median for continuity with §2–§4, and the difference is immaterial
against the 3-order window shift. With Blewett depths (median 20.8 m, max 110.5 m) and the age
bracket (7–303 kyr):

| | rate window (m/Gyr) | source |
|---|---|---|
| Phillips 2021 Table 5 (used in Phases 1–2) | **3.3 – 130** | latitude-dependent effective-erosion floors |
| **Wang 2020-implied** | **6.9×10⁴ – 1.6×10⁷** | depth/age; median 20.8 m @ 103 kyr = 2.0×10⁵ |

The Wang requirement is **~2.7 dex above the top** of the window we used throughout. Both are
reported here; the switch is not silent. Re-screening the surviving candidates at the 78.53°N cold
limit against **both** windows (`runs/wang_window_screen.json`):

| candidate / route | rate (m/Gyr) | source | vs Phillips 2021 | vs Wang 2020 |
|---|---|---|---|---|
| elemental S (thermal) | 7.7×10¹¹ | this work | above (supply-limited) | **above +4.7 dex (supply-limited — SURVIVES)** |
| C60 (thermal) | 2.5×10³ | this work | above +1.3 dex | **TOO SLOW −1.4 dex** |
| C70 (thermal) | 3.1×10² | PI-supplied | above +0.4 dex | **TOO SLOW −2.3 dex** |
| NaCl (thermal) | 3.8×10⁻⁴ | PI-supplied | below | too slow −8.3 dex |
| KCl (thermal) | 4.3×10⁻³ | PI-supplied | below | too slow −7.2 dex |
| NaCl VUV photolysis (optimiztic Y) | 56 | PI-supplied (athermal track) | in-window | **TOO SLOW −3.1 dex** |
| native Na PSD (best case Y=10⁴) | 1.1×10⁴ | PI-supplied (athermal track) | above | too slow −0.8 dex |
| cusp sputtering | 0.1 | PI-supplied (athermal track) | below | too slow −5.8 dex |
| graphite (all routes) | 5.5×10⁻⁶⁶ | PI-supplied (athermal track) | below | too slow −70 dex |

*Provenance: the elemental-S and C60 thermal rates are independently recomputed in this track on the
corrected SI-Pa convention (7.70×10¹¹ and 2.54×10³ m/Gyr) and reproduce the PI's preliminary-check
values (8.1×10¹¹, 2.6×10³) to within 5% and 2% respectively — a genuine cross-check. C70, the
chloride thermal rates, and all athermal-route rates are the PI's supplied values from the athermal
track (this track does not re-derive them); the dex verdicts are unchanged by the small S/C60
differences.*

**Only elemental sulfur exceeds the Wang requirement, and it does so with the right sign** — the
supply-limited regime this track already identified (sulfur cannot self-limit diffusively; §1, §3).
C60/C70 **flip** from "thermophysically viable" (old window) to too slow. Every **athermal** route
fails by 0.8–8.3 dex: the athermal rehabilitation of the chlorides — the earlier headline
methodological result — **does not survive the Wang age.** As the PI notes, this is not a loss but a
stronger negative: *no athermal route can excavate a 20 m hollow in ~100 kyr*, which independently
eliminates an entire mechanism class.

**Direct support for §2.1 (no self-limiting ceiling → supply/age-limited).** Wang's finding that
hollow depth is set by the **volatile-bearing layer thickness** (not by sublimation duration) is
the same supply-limited signature §2.1 inferred independently from the depth *distribution*
(memoryless exponential tail, no ceiling). The two arguments are independent (our tail-shape
statistic on the Blewett n=2518 depths; Wang's layer-thickness reasoning from the young model age
and halo/rim-height budget) and converge: **hollow depth reflects how much volatile was locally
available, not how long erosion ran.** In this supply-limited picture the depth-latitude *slope* is
not an age/rate diagnostic — consistent with §4.3, where its magnitude was already withdrawn as
provisional.

**Excavation-timescale revision.** My earlier excavation bracket (189 Myr – 5.5 Gyr) was built on the
Phillips 2021 window and is correct *for that window*, but it **does not apply** under Wang: formation
is ~100 kyr and depth is layer-thickness-limited, so the 189 Myr–5.5 Gyr figure is 3.3–4.7 dex too
long and is withdrawn as the operative excavation time. The operative picture is a **recently exposed,
concentrated crustal reservoir** (impact-excavated LRM layer) sublimating sulfur essentially at the
supply rate — not a slow multi-Gyr erosion of a dispersed phase.

---

*Deliverables: `scripts/lag_morphology.py` (poleward limit parameterized via `POLEWARD_LIMITS`;
both rate windows via `WINDOW_PHILLIPS2021` / `WINDOW_WANG2020` and `screen_against_windows()`),
`runs/lag_screening_results.csv` (lag/purity/window-reach at both 71.4 and 78.53°N + both-window
verdicts), `runs/poleward_limit_discriminant.json` (depth-latitude discriminant + mechanism
profiles), `runs/wang_window_screen.json` (all-routes screen against both windows),
`runs/lag_morphology_figure.png`, `runs/poleward_limit_figure.png`. Validation: sulfur hot-limit lag
15.00 m / 63.1% and cold-limit 99.34% reproduce the published values. Convention: SI-Pa vapor
pressures (settled Phase-1 standard); physical pore radius 0.464 µm from the sulfur gate. Poleward
limit: 78.53°N (Blewett 2016) with 71.4°N (Thomas 2014) reported alongside; HORNET (Bickel 2025)
pending. Rate windows: BOTH Phillips 2021 (3.3–130 m/Gyr) and Wang 2020 (6.9×10⁴–1.6×10⁷ m/Gyr,
doi:10.1029/2020JE006559) reported; only elemental sulfur survives the Wang window.*
