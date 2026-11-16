# TH-5 — absorbing boundary condition: step history

The §7 row and entry carry the status and the Done-when; this file carries
the step narrative (`OPS-46`). Earlier steps (1a, attempts 1–2) are recorded
in the entry itself.

## Step 1b-0 — the dipole / sphere-reflection anchor, numpy only (§9 item 107, 2026-09-30 12:00 slot)

**Why a sphere, not the box (10:30 review, journal
`docs/planning/reviews/2026-09-30-daily-1030.md` §4).** The entry's original
step-1b clause cannot produce a band, for two reasons:

1. **R(θ) is planar.** Jin (9.48) is the reflection of a plane wave off a flat
   first-order ABC. Its rms over a cube, weighted by the dipole's power
   pattern sin²θ_z, is ρ ≈ 0.113 whatever the box size (measured here:
   0.11309, max \|R\| 0.26795 in the corner direction). No box size makes
   that number smaller, so it cannot set a band that tightens as the box
   grows.
2. **Curvature is missing.** On a curved boundary the first-order condition
   also reflects through the 1/(2(k₀R)²) term: 5.1 % at 0.5 λ₀. A box has
   no closed-form reflection. A sphere has one.

**The anchor.** Take a z-directed Hertzian dipole with e^{+jωt}, x = k₀r,
C = 1, f = 1 − j/x − 1/x² and g = 1 − j/x:

- E_r = 2C cosθ (1/x² − j/x³) e^{−jx}
- E_θ = C sinθ (j f/x) e^{−jx}
- η₀H_φ = C sinθ (j g/x) e^{−jx}

These are the standard small-dipole fields (Balanis, *Antenna Theory*,
eqs. 4-8 / 4-10). That citation is from memory and **not verified in-repo**;
the converted Jin text does not carry the dipole. For that reason (n1) and
(n2) verify the field independently.

On r = R the first-order ABC n × ∇ × E + jk₀ n × n × E = 0 (Jin §9.2 (9.44)
at s = 1, §9.3 (9.54)) reduces to E_θ = η₀H_φ. Write the field as
E = A·E_out + B·conj(E_out); the incoming wave's field is H_in = −conj(H_out).
For a general surface impedance E_θ = z η₀ H_φ this gives

  Γ(z) = B/A = (f − z g) e^{−2jx} / (f̄ + z ḡ)

and four special cases:

- z = 1 (the ABC): Γ_sph = −e^{−2jx}/(2x² − 1 + 2jx), with \|Γ_sph\| = 1/√(1 + 4x⁴).
- z = f/g (the exact impedance): Γ = 0.
- z = 0 (PEC): Γ = (f/f̄) e^{−2jx}, with \|Γ\| = 1.
- z = −1 (γ → −γ): \|Γ\| = √(1 + 4x⁴).

The ABC reflects a fraction 1/(2x²) of the wave for x ≫ 1, and fully as x → 0.

**Verified (`tests/validation/test_th5_dipole_sphere.py`,
`docs/testing/logs/20260930T171635Z_TH-5-step1b0-W1.log`, 46 s, 17 passed,
`-n 2`, complex image; numpy only):**

- **(n1)** ∇×∇×E − k₀²E by nested central differences, at 20 off-axis points
  in 0.2 ≤ r ≤ 0.5:
  - δ = 1e-3: 3.614e-4. δ = 2e-3: 1.445e-3. The ratio is 3.998, i.e. pure δ²
    truncation.
  - δ = 2.5e-4: **2.259e-5 ≤ 1e-4** (`:88`).
  - **Step-size finding.** At the item's δ = 1e-3 the 1e-4 bound is red by
    3.6×. The item's trap clause ("a red by a factor ≲ 10 with the residual
    falling as δ² is a step-size finding") applies. The bound was not
    touched: it is asserted at 2.5e-4, and the δ² ratio is asserted in
    [3.8, 4.2]. The first window is `20260930T171523Z_TH-5-step1b0-W1.log:88`.
  - (H) The analytic η₀H is (j/k₀)∇×E to 3.902e-6 at δ = 1e-4 (`:89`).
- **(n2)** The radiated power at r = 0.3 equals 4πC²/(3η₀k₀²) to 1.33e-15
  (64-point Gauss–Legendre, exact for this polynomial integrand).
- **(n3)** On the standing wave E_out + Γ_sph·conj(E_out), with H from a
  central-difference curl of E at δ = 1e-4, max \|E_θ − η₀H_φ\|/\|E_θ\| over
  20 polar angles:
  - x = π: 2.427e-7.
  - x = 2π: 4.817e-8.
  - At δ = 1e-3 both are 100× larger, the δ² signature.
  - Control: Γ = 0 misses by 0.106275, equal to \|f − g\|/\|f\|.
- **(n4)** Reflection magnitudes:
  - \|Γ_sph(π)\| = 0.0505957066, \|Γ_sph(2π)\| = 0.0126641323, ratio 3.9952.
  - \|Γ_sph(1e-3) − 1\| = 1.3e-9.
  - Γ(z = f/g) = 2.9e-17.
  - Γ(z = 1) matches Γ_sph to 3.5e-18.
  - The PEC case equals (f/f̄)e^{−2jx} to 1.6e-16.
  - \|Γ(z = −1)\| = 19.7645.
- **(n5)** ρ = 0.11309 against the predicted 0.113, and max \|R\| = 0.26795.

**Consequence.** The derivation holds, so the entry's conditional step-1b
clause (the FEM gate on a sphere at k₀R = π, §9 item 114) is effective. The
bar it may not freeze above is \|Γ_sph(π)\|/3 = 1.69e-2. `TH-5` stays 🟡.

## Step 1b — the dipole gate on the sphere, FEM (§9 item 114, 2026-09-30 13:30 slot)

**Fixture.** `MeshGenerator.ball_domain(radius=0.5)`, λ₀ = 1 m (k₀R = π),
degree 1, vacuum. Every dof of the cells with r ≤ 0.15 is pinned to E_out
(r clamped ≥ 0.05). The shell carries A E_out + B conj(E_out), and Γ = B/A
is read by a volume least-squares projection over 0.22 < r < 0.42 (χ is DG0
by midpoints, quadrature degree 8). E_out is written as
e^{−jx}[2cosθ(1/x² − j/x³) r̂ + (jf/x)(cosθ r̂ − ẑ)], in numpy for the drive
and in UFL for the projection, so nothing divides by ρ. New code sits in
`tests/validation/test_th5_dipole_sphere.py` (tests only); `_assemble_complex`
is imported from the step-1a module.

**W1** (`20260930T184708Z_TH-5-step1b-W1.log`, `BAND = None`, 19 passed,
214 s, `-n 2`):

- **(P)** A planted Γ = 0.01 − 0.02j is recovered to 8.4e-15 (`:108`).
- **Ladder** (`:112–124`):

  | | λ₀/20 (19 553 cells) | λ₀/40 (151 498 cells) |
  |---|---|---|
  | (E) \|Γ_E\| | 5.0600e-3 | 1.7695e-3 |
  | (A) \|Γ_A − Γ_sph\| | 5.2789e-3 | 1.8336e-3 |

  The predictions were 3.8e-3 and 9.5e-4 (× / 3). λ₀/20 is inside that
  range. λ₀/40 is 1.9× above its prediction, still inside × 3, and the
  observed rate over the pair is 1.52 against the TEM line's 2.0. Each fine
  solve took 37–45 s, above the 20–31 s the item composed.
- **Band candidate** `:125`: 3 × 1.7695e-3 = 5.31e-3, rounded up to **6e-3**.
  The ceiling is 1.69e-2, so it was frozen.
- **Controls** (printed, `:126–127`):
  - PEC: Γ = 0.782555 − 0.623784j (predicted 0.777 − 0.629j).
  - \|Γ_PEC\|/\|Γ_A\| = 21.61. The item predicted 18–22 and gave a ceiling of
    19.76. The excess over the ceiling is the discretisation error on \|Γ_A\|
    at λ₀/20, since that ceiling holds only for the exact Γ_sph.
  - The miss ratio is 198.6 (predicted 275, range 90–830).
  - The conjugated γ gives \|Γ\| = 21.65.
- **Printed only** (`:128–130`):
  - k₀R = 2π at λ/20: \|Γ_A\| 0.009915 against 0.012664, difference 4.64e-3.
    This is the λ/20 error at a smaller reflection, not a gate.
  - Degree 2 at res 0.1 (2 502 cells, R/h = 5): the (E) / (A) errors are
    3.18e-2 / 3.34e-2. That rung is too coarse to read.
  - The box leg was not run: it is printed-only and unpriced in-slot.

**W2** (`20260930T185052Z_TH-5-step1b-W2.log`, `BAND = 6e-3`, the only
edit): 19 passed in 183 s, with the same figures to the printed digit. It
asserts \|Γ_E\| ≤ B, \|Γ_A − Γ_sph\| ≤ B at λ₀/40, B ≤ 1.69e-2, and that
PEC misses Γ_sph by > B (≈ 1.05).

**Reading.** The first-order ABC term imposes E_θ = η₀H_φ on a curved,
faceted boundary as well as on the plane, to 1.8e-3 at λ₀/40, which is
27.6× below \|Γ_sph\|. The exact-impedance premise (E) reads the same error
floor, so the residual is discretisation and faceting, not the ABC. `TH-5`
stays 🟡: ✅ stays step 2 and the §2.1 line, and the weekly ratifies the
sphere re-pose.

## §7 row moved verbatim 2026-10-01 (rotate_plan_archive.py chunks)

 Absorbing boundary condition (ABC) — **feature ladder B1** (operator directive 2026-09-04, §9 item 5; HFSS *Radiation* boundary; entry below) | 🟡 *(**Step 1a DONE 2026-09-30 (04:30 take-next slot, §9 item 94, attempt 2 of 4) ⇒ `TH-5` ⬜ → 🟡: read by volume projection, the first-order ABC's TEM-line reflection falls as h² to 5.9e-5.** New `src/fem_em_solver/core/boundary_terms.py` (additive: `impedance_boundary_term`, `first_order_abc_term`, `abc_first_order_reflection_te`; no existing `src/` file or `core/__init__.py` edited) and `tests/validation/test_th5_abc_reflection.py`, both from the parked attempt; `_box`, `_solve`, `_fit_gamma` and `boundary_terms.py` byte-identical to it, the test gains the projection. **(P)** planted Γ = 0.01 − 0.02j recovered to 6.2e-15 (band 1e-10); ∫χ dV = 0.60·(2h)² to ≤ 2.7e-13 at all four rungs (band 1e-12; `20260930T094044Z_TH-5-step1a2-W1.log:94–101`). **(N), degree 1, by projection over x ∈ [0.15, 0.75] m: \|Γ\| 3.814e-3 / 9.500e-4 / 2.373e-4 / 5.932e-5 at λ₀/20 / 40 / 80 / 160** — ≤ 1e-3 at /80 and /160, strictly decreasing, **LSQ rate 2.0021** against ≥ 1.5 (`:105–112`); predicted 2.3e-3 / 5.7e-4 / 1.4e-4 / 3.5e-5 within ×3 (read ×1.66–1.69 above, same h² law). The same four solves through attempt 1's 41-point fit reproduce its ladder to four figures (4.916e-3 / 2.018e-3 / 1.531e-3 / 1.296e-3, rate 0.6169), so the solves did not change and **the degree-1 floor was the point sampling**, by the item's pre-registered test (≤ 4.3e-4 / ≤ 1.1e-4 at /80 / /160). (O1) 1.161e-4 and (O2) \|Γ + 1/3\| 2.621e-4 by the 41-point fit, inside 3e-3 (`:124–125`); numpy anchors green (`:80`). Printed only: (O2) by projection 2.815e-4 (predicted 2.6e-4); degree 2 at λ₀/40 (N) 1.369e-7 by projection against 2.319e-5 by points, (O2) 1.070e-6 against 4.122e-6; controls at λ₀/80 by projection — PEC Γ = −1.000000 + 0.000538j, \|1 + Γ(0)\| 1.1751, PMC +1.000000 − 0.000471j, 1.6183, factor 4 214 (predicted ≈ 7 000, range 2 300–21 000); (O) PEC \|1 + Γ(0)\| 2.0000, factor 2 544 (`:114–129`). **Two predictions missed low, recorded not gated (rule (e)):** Im Γ_PEC reads 5.38e-4 against a predicted 3.2e-3 (÷× 3) — the point fit's 0.0272 was the in-cell bias, and what remains is below the dispersion estimate; k_fit − k₀ reads +5.59e-4 rad/m on both the ABC and PEC solves against a predicted +1.8e-3, i.e. the fixture's numerical dispersion is ≈ 3× smaller than the `TH-6`-scaled figure; \|Γ\| at k_fit 2.126e-4 (`:117–119`). W1 16 passed, 82 s, `-n 2`, complex, first try, no edit after it. Not claimed: TM / second order / PML, the dipole gate (1b), D∞ (step 2), any open-boundary birdcage figure. ✅ stays steps 1b and 2.)* *(step 1a attempt 1 BLOCKED 2026-09-30, §9 item 90: numpy anchors, (O1) 1.161e-4 and (O2) \|Γ + 1/3\| 2.621e-4 green; (N) degree-1 TEM ladder \|Γ\| 4.916e-3 / 2.018e-3 / 1.531e-3 / 1.296e-3 at λ₀/20…/160 against ≤ 1e-3 at /80 and /160, LSQ rate 0.617 (`20260930T022458Z_TH-5-step1a-W1.log:95–102,119–121`); degree 2 at λ₀/40 prints 2.3e-5. Unattributed; the executor's hypothesis is fit-point spacing off the cell lattice (0.01575 m) leaking the degree-1 in-cell error into B — (O) samples at 2h and reaches 1e-4. Code on `attempt/TH-5-step1a-20260930T022458Z`. **Ruled 2026-09-30 03:00 review (journal §4):** no band moves; (N)'s Γ is read by a volume least-squares projection, with this 41-point ladder printed beside it — attempt 2 of 4, §9 item 94, entry amended.)* *(**Step 1b DONE 2026-09-30 (13:30 take-next slot, §9 item 114, attempt 2 of 4 since 1a's gate) — the first-order ABC on a curved boundary: the dipole on the ball k₀R = π reflects Γ_sph(π) to 1.83e-3, band 6e-3 frozen from W1; 🟡 held (✅ stays step 2 and the §2.1 line; the weekly ratifies the sphere re-pose).** (P) 8.4e-15; at λ₀/40 (151 498 cells) \|Γ_E\| 1.7695e-3 and \|Γ_A − Γ_sph\| 1.8336e-3, both against 6e-3 ≤ \|Γ_sph\|/3; PEC misses Γ_sph by ≈ 1.05 > B. `20260930T184708Z_TH-5-step1b-W1.log:108–131` (band None), `…185052Z_TH-5-step1b-W2.log` (asserted, 19 passed, 183 s). Narrative: `docs/planning/chunks/TH-5.md`.)* *(**Step 1b-0 DONE 2026-09-30 (12:00 take-next slot, §9 item 107):** the dipole field and Γ_sph(x) = −e^{−2jx}/(2x² − 1 + 2jx) are verified in numpy — \|Γ_sph(π)\| 0.0505957066, the first-order condition holds on the standing wave to 2.4e-7, and the radiated power matches to 1.3e-15. The (n1) curl-curl residual at the item's δ = 1e-3 is a δ² step-size finding (3.6e-4, ratio 4.00), asserted at δ = 2.5e-4 (2.3e-5 ≤ 1e-4). Item 114's gate is now effective. Narrative: `docs/planning/chunks/TH-5.md`; `20260930T171635Z_TH-5-step1b0-W1.log:88–104`. 🟡 held. **Confirmed 2026-09-30 18:00 review:** the δ change is the item's own pre-registered step-size branch (a truncation bound, not a physics band; red ≲ 10× with a measured δ² ratio 4.00, first window `20260930T171523Z_TH-5-step1b0-W1.log:88` kept); no bound moved.)* 

### §9 item 277 — 2026-10-06 (lane 2, TH-5 step 2a)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §9 rule (5) lane 2, the weekly's enumerated target (PROJECT_PLAN §10, "Lane 1 and lane 2, targets for the daily": "`TH-5` 2a / 2b as re-posed 09-30"); the re-pose is in `docs/planning/plan-archive.md` (grep `step 2 — re-posed in two parts`). **No band, tolerance, Done-when or validation target is set or moved here.** The [0.3, 0.7] band and the sign clause are the weekly's; every other number is imported by name or labelled predicted.

**`TH-5` step 2a — the first-order ABC on `PORT-1`'s two-torus box: a banked negative.** Three solves at `PORT-1`'s 3b-xi paddings W = 0.08 / 0.10 / 0.12 m, with the PEC wall replaced by the first-order ABC (Jin 3e §9.2 (9.44) at s = 1; plane-wave reflection (9.48); the curved-boundary form §9.1.2 (9.21); this file's 1b-0 Γ_sph). At 10 MHz the box sits at k₀·half-width ≈ 0.022–0.035 (k₀ = 0.20958 rad/m), where \|Γ_sph\| = 1/√(1 + 4x⁴) ≈ 1 − 2e-6. There the ABC is a reflector, equal to the natural (magnetic) wall to O(k₀R). The sphere model's static reaction fields (+2 under PEC, −1 under a magnetic wall; entry step-2 note) predict that the box term **flips sign and halves**. The item measures that, and banks the entry's original "≥ 5× fall" as unreachable.

- **Executor:** implementer. **Tier:** heavy (one window ≈ 5–6 min, over standard's 180 s), complex build, `-n 2` (`PORT-1` / `GEO-4`'s width). **Class:** lane-2 solve. **Cost:** ≈ 7 measured (the weekly's lane median) / **≈ 40 predicted** min: 15 fixed + ≈ 12 code + W1 ≈ 6 + W2 ≈ 1 + ≈ 6 spare. **If W1 is not green by minute 40, park** (rule (d)).
- **Depends on:** nothing. **Attempt:** 1 of 4. **No `src/` edit:** the ABC term exists (`fem_em_solver.core.boundary_terms.first_order_abc_term`, step 1a), and the solver takes it through `TimeHarmonicSolver.solve(extra_bilinear_terms=…)` (`core/time_harmonic.py:400, 556–559`, `PORT-9`'s hook). No re-gate under rule (c).

**Settled from the code (read before designing).**
- **The box term, as the entry defines it.** The §7 entry's step 2 says "the box term `D∞`", and its step-2 note says "the box term would flip sign and about halve". Both mean B(W) = deficit(W) − D∞, where D∞ is the free-exponent extrapolated box-free deficit. **Read here (for the review to ratify):** the weekly's "\|D∞^ABC\| / \|D∞^PEC\|" is the ratio of *box terms*, r(W) = \|B^ABC(W)\| / \|B^PEC(W)\|, and the sphere model's 0.5 is that ratio. Under the same model the extrapolated *box-free* deficits would agree (ratio ≈ 1), so the literal reading contradicts its own 0.5.
- **PEC leg: no new solve beyond the premise.** `PORT1_PADDING_RUNGS_M` (0.08, 0.10, 0.12) / `PORT1_PADDING_DEFICITS_PP` (−8.0324, −5.0256, −3.2733) and `_free_exponent_fit` (`src/fem_em_solver/io/mesh.py:436–440`) give D∞ = +1.6934077797 pp, p = 1.6574 (`GEO-4` (A0), bitwise the probe, `docs/planning/chunks/GEO-4.md:11`). Hence B^PEC = −9.7258 / −6.7190 / −4.9667 pp at the three rungs. The common reference is D∞^PEC: the mesh, drive and discretisation are shared per rung, so the discretisation part of the deficit cancels in B.
- **Fixture.** `MeshGenerator.two_torus_domain` returns `(msh, cell_tags, facet_tags)`, and **facet tag 1 is the outer boundary** (`io/mesh.py:1353–1354`). Deficit = 100·(Im Z₂₁ / ωM₁₂ − 1) on one driven column (tag 1), as in `GEO-4`'s `gate_deficits_pp` (`tests/validation/test_geo4_box_padding_heuristic.py:236–249`). The reaction is `Z = −∫_tag E·J / (I_d I_t)` (`test_port_reaction_impedance.py:227–239`).
- **`_solve_reaction_z` hard-codes `boundary_condition="pec_zero_tangential_a"` and takes no surface terms** (`test_port_reaction_impedance.py:261–317`). Do **not** add a kwarg to it: that edits `PORT-1`'s gated module and buys its re-gate. Instead, copy its mesh / current / solve / `_reaction` lines into the new module as `_solve_z_column(padding, wall)` with `wall ∈ {"pec", "natural", "abc"}`, where `"abc"` = `boundary_condition="natural"` **plus** `extra_bilinear_terms=[first_order_abc_term(msh, facet_tags, 1, OMEGA)]`. Import the constants (`MAJOR_RADIUS`, `MINOR_RADIUS`, `SEPARATION`, `H_FAR`, `H_WIRE`, `CURRENT_A`, `TORUS_TAGS`, `OMEGA`, `FREQUENCY_HZ`), `_azimuthal_current_density`, `_tag_volume`, `_reaction` and `mutual_inductance`. Skip the step-2f energy bookkeeping.
- **The volume form's scaling matches the ABC term's.** The solver assembles (1/μᵣ)∫∇×u·∇×v̄ − k₀²∫ε_c u·v̄ with load −jωμ₀∫J′·v̄ (`time_harmonic.py:520–549`). That is the convention `boundary_terms.py:3–10` assumes (γ = +jk₀).

**New test module `tests/validation/test_th5_abc_two_torus.py`** (tests only). One module-scoped fixture runs five solves: PEC at 0.08 m (premise), ABC at 0.08 / 0.10 / 0.12 m, and natural at 0.08 m (control). `@complex_only` on every solve node.

**Anchors (asserted).**
- **(P0) Premise.** The copied solver at W = 0.08 m under PEC reproduces `GEO-4`'s 0.11-image rung −8.0404 pp within `P0_BAND_PP` (0.05 pp, imported from `test_geo4_box_padding_heuristic.py`). This proves the copy is the fixture.
- **(S) Sign (the weekly's clause).** B^ABC(W) = deficit^ABC(W) − D∞^PEC **> 0** at all three rungs, against B^PEC(W) < 0. Predicted ≈ +4.86 / +3.36 / +2.48 pp.
- **(R) Ratio (the weekly's band).** \|B^ABC(W)\| / \|B^PEC(W)\| ∈ [0.3, 0.7] at all three rungs; the sphere model gives 0.5. Predicted ABC deficits (r = 0.5): +6.556 / +5.053 / +4.177 pp. The band maps to [4.611, 8.501] / [3.709, 6.397] / [3.183, 5.170] pp.
- **(E) The ABC term is in the matrix (sign asserted).** Under ABC, Re Z₁₁ > 0 (the boundary absorbs). Under PEC, \|Re Z₂₁\| ≤ 1e-9·\|Z₂₁\| (`PORT-1`'s structural-zero test, `IDENTITY_TOLERANCE`).
  - Printed beside it, not asserted: the Galerkin power identity Re Z₁₁·I₁² = ∮_{tag 1}\|n × E_h\|² ds / η₀. It is the imaginary part of the weak form at v = E_h; ∫_tag E·J = ∫E·J′ still holds under the ABC because n × ∇q = 0 for q ∈ CG1 ∩ H¹₀. Relative residual **predicted ≤ 1e-6**.
  - Re Z₁₁ itself is predicted at ≈ 1e-4–5e-3 Ω and printed: the free-space loop R_rad of 9.7e-7 Ω times the near-field factor ≈ 1 + 1/(k₀R)² ≈ 800–2 000. Here the ABC acts as a resistive η₀ wall on the near field.

**Negative controls (rule (e)).**
- **(c1) PEC in the ABC's place — asserted.** A broken ABC run that leaves the Dirichlet set on, or drops the term into a PEC solve, returns the PEC deficit. (S) separates the two by \|B^ABC − B^PEC\| ≈ 14.6 pp at 0.08 m (predicted), against a rung precision of 5e-5 pp. Asserted through (S) and (P0) on the same mesh.
- **(c2) Conjugated γ — asserted via (E).** γ = −jk₀ makes the boundary a source and gives Re Z₁₁ < 0 (the 1b trap: \|Γ(z = −1)\| = 19.76). (E)'s sign catches it.
- **(c3) Natural wall at 0.08 m — predicted, printed only.** \|deficit^ABC − deficit^NAT\| is predicted ≤ 0.05·\|B^PEC(0.08)\| ≈ 0.5 pp (expected ≈ (k₀R)²·\|B\| ≈ 0.01 pp). It shows that the ABC is the magnetic wall here. It is **not** step 2b's combination and gates nothing.

**Printed (not asserted).**
- Cell counts per rung (0.08 m: 119 533 on the 0.11 image, `GEO-4.md:13`), per-solve seconds, and the complex Z₂₁ and Z₁₁.
- The free-exponent fit through the three ABC deficits (`_free_exponent_fit`): D∞^ABC, p^ABC, and \|D∞^ABC − D∞^PEC\|, predicted ≈ 0 under the sphere model. The fit may fail if the ABC ladder is not monotone; if so, print the failure.
- The entry's original bar, \|B^ABC\| ≤ \|B^PEC\|/5. (R)'s lower bound 0.3 > 0.2 makes it false by construction: **that is the banked negative**.

**Windows.** Foreground; each harness call bare and whole; Bash `timeout: 660000`; `-s`; node IDs only, never `-k`; `tests/environment` first; no pipe or redirect (the window is under 10 min, so no durable capture).
- **W1** (heavy). Predicted ≈ 330 s, bracket 250–450 s: environment ≈ 45 s plus five mesh-and-solves. The rung costs come from `GEO-4` at 43 / ≈ 52 / ≈ 70 s for the 0.08 / 0.10–0.11 / 0.12–0.14 m rungs (`20261003T035018Z_GEO-4-W1.log`: 226 s for three solves).
  `scripts/testing/run_and_log.sh TH-5-step2a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_th5_abc_two_torus.py -v -s --tb=short'"`
- **W2** (real build: collection and skips, ≈ 15 s).
  `scripts/testing/run_and_log.sh TH-5-step2a-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_th5_abc_two_torus.py -v -s --tb=short'"`

**Traps already paid for.**
- **Signed deficits and signed box terms** (`GEO-4`). Never compose relative percents. C is negative in the commissioning convention. (S) is on signed B; (R) is on magnitudes.
- **`project_source` stays at its default, `True`** (CG1 ∩ H¹₀, pinned exterior): the projection every PEC rung was measured with. `"matched"` would unpin the exterior under the natural mode and change the drive between the walls.
- **A conjugated γ is a source** (`boundary_terms.py:10–11`). Use `first_order_abc_term` as it is; never re-derive γ.
- **The ABC needs `boundary_condition="natural"`.** Under PEC the Dirichlet rows overwrite it silently (c1).
- **Rank-safety.** Allreduce every `assemble_scalar` (the reaction and the ∮\|n × E\|² integral). Take global counts via `index_map(…).size_global`. Print on rank 0.
- **Complex build.** No `ufl.max_value` and no ordering on complex operands. The boundary integral is `ufl.inner(ufl.cross(n, E), ufl.cross(n, E))`, with the real part taken after the allreduce. Pin `metadata={"quadrature_degree": …}` on any `SpatialCoordinate`-bearing facet integral (`POST-5`).
- **A kill at 590 s is not a result.** Run `pgrep -c python3` in the container first, and sweep orphan `.c` files only when it reads 0. Never lengthen the timeout or change `-n`. Status 124 is a finding.
- **Rule (i):** no edit to the new module after its last green window. Write ≈, never `~`, in tracked Markdown.

**Scope.**
- **Covered:** one fixture (`PORT-1`'s two-torus, d = 0.04 m, 10 MHz, degree 1, `-n 2`), with the first-order ABC on the whole outer box at three paddings.
- **Not this item:** step 2b (the ⅓ PEC + ⅔ natural mixed wall, its weights and its ≤ \|B^PEC\| / 5 gate); the §2.1 boundary-model line; any second-order ABC or PML; any birdcage or open-boundary figure; any statement about larger boxes (k₀R ≳ π belongs to 1b). `PORT-1`'s ✅ and `GEO-4`'s heuristic are untouched.

**Status it can move.** No glyph: `TH-5` stays 🟡. Its Done-when as re-posed 09-30 is ✅ on 1a, 1b, 2a, 2b and the §2.1 line. On green, the §7 row gains a step-2a record sentence (the B^ABC signs and r(W) at three rungs, the fit's D∞^ABC, elapsed time). The §10 lane-2 entry marks 2a DONE, and 2b becomes the lane's next step with the measured ratio beside the sphere model's.

**Negative result.** In every case: report the measurement, park on `attempt/TH-5-step2a-*`, mark item 277 BLOCKED in the same commit (rule (d)), and journal it in attempts.md. No band moves.
- **(P0) red:** the copied solver is not the fixture. Stop before reading (S) / (R).
- **(E) red** (Re Z₁₁ ≤ 0 under ABC): the term did not enter `a`, or entered with the wrong sign. Stop: this is a code defect, not a physics reading.
- **(S) red** (B^ABC < 0 at any rung) with (E) green: the box does not behave as a magnetic wall. That contradicts 1b-0's closed form at x → 0; a finding for the weekly.
- **(R) red with (S) green:** the cube departs from the sphere model, recorded as that (the weekly's own disposition in the plan-archive re-pose). The measured r(W) goes into the §7 row as the 2a record. The weekly decides whether 2a stands closed on the record, and re-derives 2b's weights from it.
- **W1 killed or over 450 s:** BLOCKED without a result. The review shrinks the case (drop (c3), or split the rungs over two windows).

No `PREDICTED RED:` line. The asserted clauses (P0, S, R, E) are all predicted green. The banked negative is the original ≥ 5× bar, which is printed and never asserted.

**Outcome (2026-10-07 00:00 slot, attempt 1 of 4): 🚫 BLOCKED — (R) red with (S) green.** (P0) −8.0404 pp (\|Δ\| 0.0080); (E) Re Z₁₁ > 0 at all rungs; (S) B^ABC +4.7133 / +2.3965 / +1.1256 pp; (R) r(W) = 0.4846 / 0.3567 / 0.2266 against [0.3, 0.7], red at 0.12 m. `20261007T052330Z_TH-5-step2a-W1.log` (313.7 s) on `attempt/TH-5-step2a-20261007T052330Z`. Recorded in the §7 row; the weekly rules on 2a and on 2b's weights.

### §9 item 290 — 2026-10-06 (lane 2)

*Drafted by a general-purpose agent for the 2026-10-06 18:00 daily review; the review reviews it.* Licence: §9 rule (5) lane 2, the weekly's enumerated target (PROJECT_PLAN §10, "Lane 1 and lane 2, targets for the daily": "`TH-5` 2a / 2b as re-posed 09-30"); the re-pose text is in `docs/planning/plan-archive.md` (grep `step 2 — re-posed in two parts`): "**2b (the mixed wall):** the ⅓ PEC + ⅔ natural-wall combination, weights from the sphere model's +2 / −1 reaction fields and stated before the run — asserted, its box term ≤ \|D∞^PEC\| / 5 at the smallest padding; red ⇒ the cube's weights differ, a record." **No band, tolerance, Done-when or validation target is set or moved here.** The ⅓ / ⅔ weights and the ÷ 5 bar are the weekly's; every other number is imported by name or labelled predicted.

**`TH-5` step 2b — the mixed wall on `PORT-1`'s two-torus box.** Two solves at `PORT-1`'s smallest 3b-xi padding W = 0.08 m (PEC, the premise; natural, the measurement), and two printed natural rungs at 0.10 / 0.12 m. The mixed impedance is the post-solve combination Z^mix = ⅓ Z^PEC + ⅔ Z^NAT (the §7 entry's "mixed average", PROJECT_PLAN §7 `TH-5` entry, step-2 note ≈ line 2441). It is not a boundary condition and not a face split.

- **Executor:** implementer. **Tier:** heavy (W1 predicted ≈ 230 s, above standard's 180 s), complex build, `-n 2` (`PORT-1` / `GEO-4`'s width). **Class:** lane-2 solve. **Cost:** ≈ 7 measured (lane-2 median) / **≈ 35 predicted** min: 15 fixed + ≈ 10 code + W1 ≈ 4 + W2 ≈ 1 + ≈ 5 spare. **If W1 is not green by minute 40, park** (rule (d)).
- **Depends on:** nothing. Step 2b carries its own premise (P0) and its own natural solve; it does **not** need item 277's ABC term or its result. **Code source, deterministic:** if `tests/validation/test_th5_abc_two_torus.py` is on `main` with `_solve_z_column(padding, wall)`, import it (no edit to that module, rule (i)); otherwise write the copy item 277's brief specifies (this file, §9 item 277, "Settled from the code", the `_solve_z_column` bullet), restricted to `wall ∈ {"pec", "natural"}`, into the new module. Either path is proven to be the fixture by (P0). Its effect on the glyph needs 277 (see "Status it can move"). **Attempt:** 1 of 4. **No `src/` edit**; no re-gate under rule (c).

**The anchor — the sphere model's reaction weights (stated before the run).** Take a magnetic dipole m at the centre of a sphere of radius R, quasi-static (k₀R ≈ 0.03–0.04 here; 1b-0's \|Γ_sph\| = 1/√(1 + 4x⁴) → 1). The uniform interior reaction field that cancels B_r at r = R (PEC) is −2m/(4πR³); the one that cancels H_θ (magnetic / natural wall) is +m/(4πR³). For a general l-pole the two coefficients are (l + 1)/l and −1 (interior potential C r^l P_l against the source A r^{−(l+1)} P_l, with ∂_r Φ = 0 vs Φ = 0 at R). Hence on the l = 1 term **B^NAT = −½ B^PEC** (the §7 entry's "+2/R³ under PEC and −1/R³ under a magnetic wall", ≈ line 2435), and the weights w_PEC = ⅓, w_NAT = ⅔ (2w_PEC − w_NAT = 0, w_PEC + w_NAT = 1) cancel it. The residual for l ≥ 2 is (1 − l)/(3l) of each pole's PEC term. The loops sit ±d/2 = ±0.02 m off-centre and a/R ≈ 0.36, so a residual of a few percent of \|B^PEC\| is expected, plus the cube's departure from a sphere, which has no closed form. Dipole-field form: Jackson, *Classical Electrodynamics* 3e §5.6 — **from memory, not held in the repo**; the derivation above is two lines, and the executor may check it in numpy beside the run (printed only).

**Settled from the code and the records (read before designing).**
- **Deficit and box term, as item 277 defines them** (the 03:00 review's ruling on item 277: the weekly's "D∞" ratio is a ratio of box terms). d^X(W) = 100·(Im Z₂₁^X / ωM₁₂ − 1) pp on one driven column (tag 1), `GEO-4`'s `gate_deficits_pp` form (`tests/validation/test_geo4_box_padding_heuristic.py:236–249`). B^X(W) = d^X(W) − D∞, with D∞ = D∞^PEC = +1.6934077797 pp from `_free_exponent_fit(PORT1_PADDING_RUNGS_M, PORT1_PADDING_DEFICITS_PP)` (`src/fem_em_solver/io/mesh.py:436–440`; `GEO-4` (A0)). Because the weights sum to 1, **B^mix = ⅓ B^PEC + ⅔ B^NAT = d^mix − D∞**, and d^mix is the deficit of Z^mix (the same ωM₁₂ denominator, so linear).
- **The bar, read as 277 reads D∞:** \|B^mix(0.08)\| ≤ \|B^PEC(0.08)\| / 5, with both from this window's solves. Predicted B^PEC(0.08) = −8.0404 − 1.6934 = **−9.7338 pp** (`20261003T035018Z_GEO-4-W1.log:294`; `GEO-4` (C1) measured 9.7338), so the bar is **≈ 1.947 pp**. Algebra: B^mix = B^PEC(⅓ − ⅔r) with r = \|B^NAT\| / \|B^PEC\|, so the bar is met iff **r ∈ [0.2, 0.8]**; item 277's (R) band [0.3, 0.7] lies inside it. **Literal reading, printed only:** \|D∞^PEC\| / 5 = 0.3387 pp, met iff r ∈ [0.448, 0.552] (a question for the review).
- **Fixture.** `MeshGenerator.two_torus_domain` returns `(msh, cell_tags, facet_tags)`; facet tag 1 is the outer boundary (`io/mesh.py:1353–1354`). The natural wall is `boundary_condition="natural"` with no surface term: no Dirichlet rows, no `extra_bilinear_terms`. Constants and helpers come from `tests/validation/test_port_reaction_impedance.py` (`MAJOR_RADIUS` 0.04, `SEPARATION` 0.04, `H_FAR` 0.03, `FREQUENCY_HZ` 10 MHz, `OMEGA`, `TORUS_TAGS`, `CURRENT_A`, `_azimuthal_current_density`, `_reaction`, `mutual_inductance`, `IDENTITY_TOLERANCE` = 1e-9 at `:192`). `_solve_reaction_z` hard-codes PEC (`:316`): **do not add a kwarg to it** (that edits `PORT-1`'s gated module).
- **One mesh per padding, both walls on it.** At 0.08 m build the mesh once and solve PEC and natural on it, so the discretisation part of the deficit is shared and cancels in B.

**New test module `tests/validation/test_th5_mixed_wall_two_torus.py`** (tests only). Constants at module top, before any solve: `W_PEC = 1.0 / 3.0`, `W_NAT = 2.0 / 3.0`, `MIX_BAR_DIVISOR = 5.0` (the weekly's). One module-scoped fixture runs four solves: PEC and natural at 0.08 m, and natural at 0.10 and 0.12 m. `@complex_only` on every solve node.

**Anchors (asserted).**
- **(P0) Premise.** The PEC solve at 0.08 m reproduces `P0_RUNG_PP` (−8.0324 pp) within `P0_BAND_PP` (0.05 pp), both imported from `test_geo4_box_padding_heuristic.py`. This is exactly `GEO-4`'s premise (measured there: −8.0404, \|Δ\| 0.0080).
- **(E) Lossless identity.** Under both walls \|Re Z₂₁\| ≤ `IDENTITY_TOLERANCE`·\|Z₂₁\|. There is no loss term, so the system is real up to the j of the load; `GEO-4` printed Re Z = 0.000000e+00 under PEC (`:290`).
- **(S) Sign flip of the natural wall.** B^NAT(0.08) > 0 against B^PEC(0.08) < 0. Predicted d^NAT(0.08) ≈ +6.56 pp at r = 0.5 (range +3.64 … +9.48 for r ∈ [0.2, 0.8]).
- **(M) The weekly's gate.** \|B^mix(0.08)\| ≤ \|B^PEC(0.08)\| / `MIX_BAR_DIVISOR`. Predicted ≈ 0 (sphere model), and ≤ ≈ 1.3 pp across 277's band; the bar is ≈ 1.947 pp.

**Negative controls (rule (e)).**
- **(c1) Wall not switched — asserted via (S).** A run that leaves the Dirichlet set on returns d^NAT = d^PEC bitwise and B^mix = B^PEC, missing (M) by exactly 5× (the ceiling, by construction). (S) separates the two by \|B^NAT − B^PEC\| ≈ 14.6 pp (predicted) against a rung precision of ≈ 5e-5 pp. Also assert that d^NAT(0.08) and d^PEC(0.08) differ by more than `P0_BAND_PP`.
- **(c2) Swapped weights — asserted, conditional on (M) green.** ⅔ PEC + ⅓ NAT gives B^swap = B^PEC(⅔ − ⅓r); assert \|B^swap(0.08)\| > the bar. **Ceiling computed:** \|B^swap\| / bar = 5(⅔ − ⅓r) = **2.5** at r = 0.5, and 3.0 … 2.0 across r ∈ [0.2, 0.8]. So any r that passes (M) makes (c2) miss by ≥ 2.0×: the weights discriminate, not any mixture passes. Equal weights (½, ½) miss by only 1.25× at r = 0.5: printed, not asserted.

**Printed (not asserted).**
- Cell counts per rung (0.08 m: 119 533 on the 0.11 image), per-solve seconds, complex Z₂₁ and Z₁₁ per wall, and the reaction/energy residual per solve.
- r_NAT(W) and B^mix(W) at 0.10 / 0.12 m against B^PEC from the fit (−6.7190 / −4.9667 pp; labelled fit-referenced, since no PEC solve runs there).
- The weights that would zero B^mix on this cube, w_PEC* = r/(1 + r), beside ⅓; and the literal-reading margin, \|B^mix\| against 0.3387 pp.
- If item 277 has landed: its (c3) natural deficit at 0.08 m beside this one (it should agree to print precision: same code, same mesh recipe).

**Windows.** Foreground; each harness call bare and whole; Bash `timeout: 660000`; `-s`; node IDs only, never `-k`; `tests/environment` first; no pipe or redirect (the window is under 10 min, so no durable capture).
- **W1** (heavy). Predicted ≈ 230 s, bracket 180–340 s. Priced from `20261003T035018Z_GEO-4-W1.log`: environment and overhead ≈ 53 s (226 s total, less three rungs); the 0.08 m mesh 18.4 s + solve 19.1 s (`:288`), so PEC + natural ≈ 57 s; 0.10 m ≈ 52 s and 0.12 m ≈ 65 s, interpolated between the 0.08 / 0.11 / 0.14 m rungs (`:288, 474, 671`). A natural-wall solve has no measured price on this fixture (no Dirichlet rows; MUMPS on an indefinite system). If W1 passes 340 s, drop the two printed rungs, never the asserted pair.
  `scripts/testing/run_and_log.sh TH-5-step2b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_th5_mixed_wall_two_torus.py -v -s --tb=short'"`
- **W2** (real build: collection and skips, ≈ 15 s).
  `scripts/testing/run_and_log.sh TH-5-step2b-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_th5_mixed_wall_two_torus.py -v -s --tb=short'"`

**Traps already paid for.**
- **Weights are constants, never fitted.** Solving for the w that zeroes B^mix makes (M) a tautology; w_PEC* is printed only.
- **Signed deficits and signed box terms** (`GEO-4`). Combine Z, or signed d on the shared ωM₁₂; never compose relative percents with different denominators. (S) is on signed B; (M) and (c2) are on magnitudes.
- **`project_source` stays at its default, `True`,** under both walls (item 277's trap: `"matched"` unpins the exterior under the natural mode and changes the drive between the walls). The tori sit ≥ 0.08 m from the wall, so no boundary-node gradient overlaps J′.
- **The natural wall needs `boundary_condition="natural"`.** Under PEC the Dirichlet rows win silently (c1).
- **Rank-safety.** Allreduce every `assemble_scalar` (the reaction). Take global counts via `index_map(…).size_global`. Print on rank 0.
- **Complex build.** No `ufl.max_value`, no ordering on complex operands, no `sqrt` inside a comparison; take real and imaginary parts after the allreduce.
- **A kill at 590 s is not a result.** Run `pgrep -c python3` in the container first, and sweep orphan `.c` stems only when it reads 0 (the rubric's one-call sweep). Never lengthen the timeout or change `-n`. Status 124 is a finding.
- **Rule (i):** no edit to the new module after its last green window, and never an edit to item 277's module. Write ≈, never `~`, in tracked Markdown.

**Scope.**
- **Covered:** one fixture (`PORT-1`'s two-torus, d = 0.04 m, 10 MHz, degree 1, `-n 2`), and the post-solve ⅓ / ⅔ combination at one asserted padding.
- **Not this item:** a mixed wall as a boundary condition any user can set; any second-order ABC or PML; any birdcage or open-boundary figure; a padding heuristic for the mixed combination; the §2.1 line itself (the review's). `PORT-1`'s ✅ and `GEO-4`'s heuristic are untouched.

**Status it can move.** `TH-5`'s Done-when as re-posed 09-30 is ✅ on 1a, 1b, 2a, 2b and the §2.1 line "valid for k₀R ≳ π; a reflector for small boxes". On green, the §7 row gains a step-2b record sentence (B^PEC, B^NAT and B^mix at 0.08 m, r_NAT(W) at three rungs, (c2)'s miss, elapsed time), and the §10 lane-2 entry marks 2b DONE. **If item 277 is DONE**, the row then lacks only the §2.1 line: the next review writes it and flips 🟡 → ✅ after an `auditor` PASS. **If 277 is not DONE**, 🟡 holds and 2b stands as a record awaiting 2a. The implementer holds 🟡 in either case.

**Negative result.** In every case: report the measurement, park on `attempt/TH-5-step2b-*`, mark item 290 BLOCKED in the same commit (rule (d)), and journal it in attempts.md. No band, weight or divisor moves.
- **(P0) red:** the copied solver is not the fixture. Stop before reading (S) / (M).
- **(E) red:** a loss term entered a lossless solve, which is a code defect. Stop.
- **(S) red** (B^NAT ≤ 0): first check whether d^NAT == d^PEC bitwise (c1's signature: the wall did not switch). If not, the box does not act as a magnetic wall at k₀R ≈ 0.03, contradicting 1b-0's closed form as x → 0: a finding for the weekly.
- **(M) red with (S) green:** r_NAT(0.08) ∉ [0.2, 0.8]. **The cube's weights differ**, which is the weekly's own disposition: a record. r_NAT, w_PEC* and B^mix go into the §7 row as the 2b record, and the weekly decides whether 2b closes on it.
- **W1 killed or over 450 s:** BLOCKED without a result. The review shrinks the case (drop the 0.10 / 0.12 m rungs).

**Outcome (2026-10-07 04:30 slot, attempt 1 of 4): ✅ DONE — every asserted clause green, first try.** (P0) −8.0404 pp (\|Δ\| 0.0080); (E) \|Re Z₂₁\|/\|Z₂₁\| = 0 under both walls; (S) B^PEC −9.7338, B^NAT +4.4889 pp, (c1) separation 14.2227 pp; (M) B^mix −0.2520 pp ≤ bar 1.9468 pp, r_NAT(0.08) 0.4612, w_PEC* 0.3156; (c2) 2.5647×. Printed: r_NAT 0.3359 / 0.2080 at 0.10 / 0.12 m (fit-referenced). Cells 119 533 / 135 685 / 153 952; solves 21.4 + 20.6 / 30.4 / 35.6 s. `20261007T093129Z_TH-5-step2b-W1.log:621–636` (240 s), W2 `20261007T093540Z_TH-5-step2b-W2.log` (65 s). **Correction to the anchor paragraph (printed numpy check, W1 `:79–81`):** the l ≥ 2 residual (1 − l)/(3l) is in the source normalisation; as a fraction of each pole's PEC term it is (1 − l)/(3(l + 1)) (−1/9 at l = 2, not −1/6). Gates nothing. 🟡 held: 2a is BLOCKED.

No `PREDICTED RED:` line: (P0), (E), (S), (M) and (c2) are all predicted green.


## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 Absorbing boundary condition (ABC) | 🟡 Steps 1a (first-order ABC reflection on the TEM line falling at LSQ rate 2.0021 by volume projection), 1b-0 and 1b (the dipole on a curved ball boundary) all landed 2026-09-30. ✅ still needs step 2 and the §2.1 line; the sphere re-pose was ratified by the 2026-09-30 weekly. *(**Step 2a record 2026-10-07 (§9 item 277, attempt 1, BLOCKED — (R) red with (S) green):** first-order ABC on `PORT-1`'s two-torus box, degree 1, `-n 2`, 10 MHz. (P0) PEC 0.08 m deficit −8.0404 pp (\|Δ\| 0.0080 ≤ `P0_BAND_PP`); (E) Re Z₁₁ +3.784e-3 / +2.620e-3 / +1.924e-3 Ω > 0, PEC \|Re Z₂₁\|/\|Z₂₁\| = 0, power identity Re Z₁₁·I₁² = ∮\|n × E\|²/η₀ to ≤ 1.2e-11 (printed); (S) green, B^ABC = +4.7133 / +2.3965 / +1.1256 pp against B^PEC −9.7258 / −6.7190 / −4.9667 pp; **(R) r(W) = 0.4846 / 0.3567 / 0.2266 at 0.08 / 0.10 / 0.12 m against [0.3, 0.7] — red at 0.12 m**, falling with W rather than holding at the sphere model's 0.5. ABC deficits +6.4067 / +4.0899 / +2.8190 pp; the free-exponent fit through them fails (printed); (c3) \|deficit^ABC − deficit^NAT\| at 0.08 m 0.2245 pp (predicted ≤ 0.486, printed); original ≥ 5× bar false at all rungs (banked negative). Executor's hypothesis, untested: the ABC ladder approaches a box-free limit below D∞^PEC, so B^ABC measured against D∞^PEC carries that offset and r(W) drifts. `20261007T052330Z_TH-5-step2a-W1.log:952–974` (on `attempt/TH-5-step2a-20261007T052330Z`), 313.7 s. The weekly decides whether 2a stands closed on this record and re-derives 2b's weights.)* *(**Step 2b DONE 2026-10-07 (04:30 slot, §9 item 290, attempt 1 of 4) — the ⅓ PEC + ⅔ natural mixed wall meets the weekly's ≤ \|B^PEC\|/5 bar at 0.08 m; 🟡 held (2a is BLOCKED, so 2b stands as a record awaiting the weekly's 2a ruling).** Post-solve Z^mix = ⅓ Z^PEC + ⅔ Z^NAT on one mesh, degree 1, `-n 2`, 10 MHz, new `tests/validation/test_th5_mixed_wall_two_torus.py` (tests only, own `_solve_z_column` copy). (P0) PEC 0.08 m −8.0404 pp (\|Δ\| 0.0080 ≤ 0.05); (E) \|Re Z₂₁\|/\|Z₂₁\| = 0 under both walls; (S) B^PEC −9.7338 pp, B^NAT +4.4889 pp (d^NAT +6.1823; (c1) \|d^NAT − d^PEC\| 14.2227 pp); **(M) B^mix −0.2520 pp against the bar 1.9468 pp**, r_NAT(0.08) 0.4612, w_PEC* 0.3156 against ⅓; (c2) swapped weights miss by 2.5647× (predicted 2.5); literal-reading bar 0.3387 pp also met (printed). Printed, fit-referenced: r_NAT 0.3359 / 0.2080 and \|B^mix\|/\|B^PEC\| 0.1094 / 0.1947 at 0.10 / 0.12 m — the cancellation weakens as the box grows, the same fall as 2a's r(W). `20261007T093129Z_TH-5-step2b-W1.log:621–636` (16 passed, 240 s, first try), W2 real build `20261007T093540Z_TH-5-step2b-W2.log` (9 passed, 7 skipped, 65 s).)* *History: `docs/planning/chunks/TH-5.md`.* 

### §9 item 320 — 2026-10-08 (lane 2 — TH-5 step 2c)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §9 rule (5) lane 2, the 2026-10-07 weekly's ruling (`PROJECT_PLAN.md:4454–4468`, "`TH-5` step 2a — stands on its record; (R) retired as posed"), enumerated as "**Lane 2: `TH-5` 2c** … ⇒ `TH-5` ✅" (`:4777–4778`) and flagged for this review (`:4892, 4899`). **No band, tolerance, Done-when or target is set or moved.** Every gated number below is imported by name or from a log; the §2.1 sentence is the weekly's, copied verbatim.

**`TH-5` step 2c — land 2a's module with (R) retired, re-run 2b, write the §2.1 line.** Three windows, no new physics. The change to 2a's module is one deletion: the `RATIO_BAND` assertion loop at the end of `test_s_r_box_term_sign_and_ratio`. The weekly ruled that (R) was the sphere model's 0.5 applied to a cube ladder whose reference D∞^PEC is itself a fit, and that 2b's natural wall shows the same fall (r_NAT 0.46 → 0.34 → 0.21). That is a changed quantity, not a loosened band (the T21′ pattern).

- **Executor:** implementer. **Tier:** heavy (W1 ≈ 280 s and W2 ≈ 200 s, each over standard's 180 s), complex build, `-n 2` (2a's and 2b's width). **Class:** lane-2 land-and-close. **Cost:** ≈ 10 measured (2b's slot: 04:30 CDT start → journal 09:40Z) / **≈ 45 predicted** min = 15 fixed + ≈ 5 code + ≈ 10 compute (W1 + W2 + W3) + ≈ 9 docs (row, §2.1, §10, this file) + ≈ 6 spare. **If W1 is not green by minute 35, park** (rule (d)).
- **Depends on:** nothing. Both modules' fixtures and imports are on `main`; 2b's module (`6b27711d`) is independent of 2a's (its own `_solve_z_column` copy). **Attempt:** 1 of 4. **No `src/` edit; no edit to 2b's module, `PORT-1`'s or `GEO-4`'s** (rule (c) not triggered).

**Settled from the branch (read before editing).** `git diff --name-status main...attempt/TH-5-step2a-20261007T052330Z` = `A tests/validation/test_th5_abc_two_torus.py` (261 lines), `A docs/testing/logs/20261007T052330Z_TH-5-step2a-W1.log`, `M docs/testing/test-results.md` (one row). Branch tip `0a90dee7`.
- **Take the code by path, only the code:** `git checkout attempt/TH-5-step2a-20261007T052330Z -- tests/validation/test_th5_abc_two_torus.py`. **Never path-checkout `docs/testing/test-results.md` or any `docs/testing/logs/` file from the branch** (rubric item 4, `8d4cf58` dropped nine rows). 2c's own windows append their own rows. 2a's W1 log stays cited "on the attempt branch", as the §7 row already says.
- **Edits to the module, and nothing else:** (1) delete `for r in ratios: assert RATIO_BAND[0] <= r <= RATIO_BAND[1]`; (2) delete the `RATIO_BAND` constant, and in the S/R print replace `(band {RATIO_BAND})` with `(record; (R) retired, weekly 2026-10-07; sphere model 0.5)`; (3) in the docstring, move (R) from "Asserted" to "Printed (record)" and add "§9 item 320". Keep the test names and every other line byte-identical, so the log lines read like 2a's. Before W1, show `git diff attempt/TH-5-step2a-20261007T052330Z -- tests/validation/test_th5_abc_two_torus.py` in the journal: it must contain only those three hunks.

**Anchors (asserted; all inherited unchanged from item 277's brief, this file `:165–171`).**
- **(P0)** PEC 0.08 m deficit within `P0_BAND_PP` of `P0_RUNG_PP` (`test_geo4_box_padding_heuristic.py`). 2a record −8.0404 pp, \|Δ\| 0.0080.
- **(E)** Re Z₁₁^ABC > 0 at all three rungs (record +3.784e-3 / +2.620e-3 / +1.924e-3 Ω), and PEC \|Re Z₂₁\| ≤ `IDENTITY_TOLERANCE`·\|Z₂₁\| (record 0).
- **(S)** B^PEC < 0 and B^ABC > 0 at 0.08 / 0.10 / 0.12 m (record B^ABC +4.7133 / +2.3965 / +1.1256 pp).
- **2b re-run:** its own (P0), (E), (S), (M), (c2) as on `main`. Record: (M) B^mix −0.2520 pp against bar 1.9468 pp; (c2) 2.5647×.
- **Printed record, not asserted:** r(W) = 0.4846 / 0.3567 / 0.2266; the power-identity residual (≤ 1.2e-11); (c3) 0.2245 pp; the failed free-exponent fit; the original ≥ 5× bar false. Same code, image and `-n 2`, so everything is **predicted to reproduce 2a's / 2b's digits as printed**.

**Negative controls (rule (e)).**
- **(c1) PEC in the ABC's place — asserted via (S) + (P0).** Separation \|B^ABC − B^PEC\| 14.4391 pp at 0.08 m (2a record) against a rung precision of ≈ 5e-5 pp.
- **(c2) Conjugated γ — asserted via (E)'s sign.** It makes the boundary a source and gives Re Z₁₁ < 0.
- **2b's (c1) / (c2) — asserted in 2b's module.** Swapped weights miss by 2.5647× (ceiling 2.0–3.0 for r ∈ [0.2, 0.8], item 290).
- **(L) the landed code is 2a's code — journal check, not a test.** The three-hunk diff above. A module rewritten instead of lifted voids 2a's record as the comparand.

**Windows.** Foreground; each harness call bare and whole; Bash `timeout: 660000`; `-s`; node IDs only, never `-k`; `test_complex_mode.py` first. Each window is under 10 min, so no durable capture. Priced from `20261007T052330Z_TH-5-step2a-W1.log` (313.68 s with all of `tests/environment`, footer 315 s) and `20261007T093129Z_TH-5-step2b-W1.log` (240 s) / `…093540Z_…-W2.log` (65 s). Dropping the other three environment files saves ≈ 35–45 s.
- **W1** (heavy; 2a's module). Predicted ≈ 280 s, bracket 230–450 s.
  `scripts/testing/run_and_log.sh TH-5-step2c-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_th5_abc_two_torus.py -v -s --tb=short'"`
- **W2** (heavy; 2b re-run). Predicted ≈ 200 s, bracket 160–340 s.
  `scripts/testing/run_and_log.sh TH-5-step2c-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_th5_mixed_wall_two_torus.py -v -s --tb=short'"`
- **W3** (real build: collection and skips for both modules, ≈ 30–65 s).
  `scripts/testing/run_and_log.sh TH-5-step2c-W3 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_th5_abc_two_torus.py tests/validation/test_th5_mixed_wall_two_torus.py -v -s --tb=short'"`
- W1 and W2 run serially, never concurrently, and never beside another heavy window. Combined into one window they would total ≈ 505 s against 590 s, which is too thin.

**Traps already paid for.**
- **Path-checkout of records** (above). After the commit is staged, `git diff --cached --stat` must show only insertions on `docs/testing/test-results.md`.
- **(R) is retired, not relaxed.** No replacement band and no new constant that reads as a bound. "Sphere model 0.5" is a printed comparand. Never re-assert r(W) at a wider interval.
- **Rule (i):** the module is edited once, before W1, and never after its last green window. 2b's module is not touched at all.
- **`project_source` stays `True`; the ABC needs `boundary_condition="natural"`; a conjugated γ is a source.** All three are already in the lifted code (item 277 traps, `:191–193`), so do not touch them.
- **Rank-safety and complex build:** already in the lifted code (allreduce before the real part; the `ds` integral has no `SpatialCoordinate`).
- **A kill at 590 s is not a result:** run `pgrep -c python3` in the container first; use the rubric's one-call orphan-`.c` sweep only when it reads 0; status 124 is a finding. Never lengthen the timeout or change `-n`.
- **Branch deletion is not this item's job** (flag 6, `PROJECT_PLAN.md:4909–4910`, is the review's). Write ≈, never `~`.

**Scope.**
- **Covered:** `PORT-1`'s two-torus at 10 MHz, degree 1, `-n 2`, at three paddings (k₀·half-width ≈ 0.02–0.04), plus 1a / 1b's closed forms at k₀R ≳ π.
- **Not this item:** the band k₀R ∈ (0.04, π) between them, which includes the 64 / 128 MHz birdcage box (§7 entry step 3's question stays open); a second-order ABC or PML; a user-settable mixed wall; any open-boundary birdcage, SAR or HFSS-*Radiation* parity figure; any reading of *why* r(W) falls (cube vs sphere, or the D∞^PEC reference — 2a's untested hypothesis stays untested).

**Status it can move.** On W1–W3 green, with the printed records reproducing 2a's / 2b's digits:
- the §7 row (`PROJECT_PLAN.md:2136`) 🟡 → ✅, with a step-2c record sentence (passed counts, elapsed, logs);
- §2.1 gains the weekly's line verbatim (`:4464–4467`) as a new bullet after "Layered impedance boundary" (`:360`);
- the §10 lane-2 entry marks 2c DONE;
- an outcome line goes under this brief.

The next review owes an `auditor` run ("audit TH-5 closed at <sha>") and the flag-6 branch deletion.

**Negative result.** In every case: report the measurement, park on `attempt/TH-5-step2c-*`, mark item 320 BLOCKED in the same commit (rule (d)), and journal it in attempts.md. 🟡 holds, the §2.1 line is not written, and no band moves.
- **(P0) / (E) / (S) red in W1, or any 2b assertion red in W2:** something on `main` changed under a fixture that was green on 10-07 (image, an imported helper, `boundary_terms.py`). That is a finding for the review, not for re-posing here.
- **Every assertion green, but a printed record drifts** (r(W) at two figures, \|B^mix\| at 0.25 pp, or (P0) beyond its fourth decimal): do **not** write the §2.1 line, because its numbers would no longer be what the logs read. Record both readings in the row; the review rules.
- **W1 or W2 killed, or over 450 s:** BLOCKED without a result. The review shrinks the case (split 2a's rungs).

No `PREDICTED RED:` line: every asserted clause is predicted green, on a record of the same code.


**Review rulings (2026-10-08 03:00 daily review, item 320).** (1) **The 2a W1 log lands on `main` by path** (`git checkout attempt/TH-5-step2a-20261007T052330Z -- docs/testing/logs/<that log>` — a new file, nothing overwritten) and its `test-results.md` row **by hand**, in the item's commit, so the §7 citation survives the branch's deletion; never a path checkout of `test-results.md`. (2) The implementer flips 🟡 → ✅ and writes the weekly's §2.1 line verbatim, as drafted; the next review audits. (3) The drift stop is confirmed: all asserts green but a printed record that does not reproduce 2a / 2b to their printed digits ⇒ no §2.1 line, hold 🟡, report. (4) The line goes where drafted (after "Layered impedance boundary"). (5) Yes: the slot adds a rule (j) annotation to the stale §7 entry-body text (`:2382`, `:2399–2400`, `:2434–2445` at drafting) in the same commit; no Done-when is moved.

**Outcome (2026-10-08 13:30 slot, attempt 1 of 4): ✅ DONE — every asserted clause green, every printed record reproduces 2a / 2b to its printed digit.** (L) diff against the branch: the three edits, as four hunks (docstring; `RATIO_BAND` constant and print text; assert loop). W1 7 passed, 225 s: (P0) −8.0404 pp (\|Δ\| 0.0080); (E) Re Z₁₁^ABC +3.784e-3 / +2.620e-3 / +1.924e-3 Ω, PEC \|Re Z₂₁\|/\|Z₂₁\| 0, power residual ≤ 1.216e-11; (S) B^ABC +4.7133 / +2.3965 / +1.1256 pp against B^PEC −9.7258 / −6.7190 / −4.9667; r(W) 0.4846 / 0.3567 / 0.2266 (record); (c1) 14.4391 pp; (c3) 0.224456 pp; ABC fit fails (printed) — `20261008T183105Z_TH-5-step2c-W1.log:928–943`. W2 8 passed, 169 s: B^PEC −9.7338, B^NAT +4.4889, (c1) 14.2227 pp, (M) B^mix −0.2520 pp vs 1.9468, (c2) 2.5647×, w_PEC* 0.3156 — `20261008T183455Z_TH-5-step2c-W2.log:597–613`. W3 real build 2 passed, 9 skipped, 2 s — `20261008T183751Z_TH-5-step2c-W3.log`. 2a's W1 log and its row now on `main` (ruling (1)). §7 row ✅, §2.1 line written verbatim after "Layered impedance boundary", §10 lane 2 marked, rule (j) annotations added at the three stale entry-body spans. Owed by the next review: the `auditor` run and the 2a branch deletion.

### §9 item 596 — 2026-10-22 (`TH-5` step 3, lane 2)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** §10, 2026-10-21 weekly: "**`TH-5` lane 2 — step 3 opened (the PML),** the one Tier B member left: a polynomial-graded PML on the TEM line, its reflection against the graded-PML closed form (Jin 3e ch. 9 — the daily cites the equation from `INDEX.md` when it briefs), ≤ 1e-4 at λ₀/40 where the first-order ABC reads 9.5e-4 (step 1a's record); the ABC's own ladder is the negative control. Additive `boundary_terms.py`; no birdcage use." (`PROJECT_PLAN.md:4826–4831`). The step is also listed as "Lane 2: `TH-5` step 3 (the PML)" (`:4999–5000`), priced in the weekly's supply sheet at the chain-solve median 16.7 (`:5014`), and queued last in the flagged order (`:5114`, `:5122`). The §7 entry's own step 3 is "(optional, review's call)" (`:2464–2466`). **No band, tolerance, Done-when or validation target is set or moved.** The 1e-4 bar is the weekly's. The PML design constants below (thickness, grading order, the two design reflections) are the item's fixture parameters, stated before the run. **The read the bar is applied to is a question for the review (For the review, 1), because at step 1a's own read the bar sits under that read's floor.**

**The anchor in Jin (cited from the converted text; the equation bodies are not in it).** `INDEX.md` lists only section titles, and states that equations are placeholders `[eq (n)]` numbered as in print (`docs/references/jin-fem-3e/INDEX.md:3–6`). The PML is §9.6, `09_absorbing-boundary-conditions.md:761–1209`:
- **§9.6.1 (9.133)** (`:877–881`): the metal-backed PML with a constant stretch s = s′ − js″. Its reflection is "reduced systematically by increasing the value of s″L", with \|R\| = exp(−2κ₂ s″ L cosθ) (the form the prose describes, assuming κ₂ is real).
- **§9.6.4 (9.165)** (`:1039–1043`): the polynomial profile of s_z across N layers, "where m can be 0, 1, 2, or 3", with δ_max "the maximum loss tangent" (`:1047`).
- **(9.166) / (9.167)** (`:1045–1049`): Jin's *discrete* R^TM / R^TE of the N-layer PML, computed by the 1D FEM of §3.4. They are not closed forms.
- **§9.6.2 (9.152)–(9.154)** (`:961–973`): the anisotropic absorber Λ.
- **§9.6.3 (9.157)** (`:993–995`): the FE vector wave equation with Λ.
- `:1027` is the item's premise: "Although the PML interface is reflectionless in the continuum space, this is not the case in the discretized space [41]."

**The graded closed form is (9.133) with s″L replaced by its integral over (9.165)'s profile.** For s(ξ) = 1 − j s″_max (ξ/d)^m on 0 < ξ < d, PEC-backed, e^{+jωt}, θ = 0:

  Γ_c = −e^{−2jk₀d} · exp(−2k₀ s″_max d / (m + 1)),

or in σ form (s″ = σ/(ωε₀)), \|R\| = exp(−2η₀ σ_max d cosθ / (m + 1)). **It is not a numbered equation in the converted text.** That (9.165) is s_z = 1 − jδ_max(z/L)^m, i.e. δ_max ≡ s″_max, is the drafter's reading of the prose, unverified against print (For the review, 2).

**What exists (read 2026-10-22).**
- **`src/fem_em_solver/core/boundary_terms.py` (61 lines, `ac65e2fe`):**
  - the convention docstring (`:3–16`: e^{+jωt}; the weak form ∫∇×u·∇×v̄ − k₀²∫u·v̄ + …; "Nothing here is exported from `core/__init__.py`");
  - `first_order_abc_term` (`:33–41`), `_tangential_mass_term` (`:44–51`), `abc_first_order_reflection_te` (`:54–61`).
  - Importers: the three `TH-5` modules (`test_th5_abc_reflection.py`, `test_th5_dipole_sphere.py`, `test_th5_abc_two_torus.py`) and `study.py:2276` (`impedance_boundary_term` only). No generated API page names the module.
- **The TEM-line fixture, `tests/validation/test_th5_abc_reflection.py` (393 lines):**
  - λ₀ = 1 m, `K0`, `OMEGA` (`:58–61`); tags (`:63`); `L_N` = 0.9 (`:66`); `PROJ_N` = (0.15, 0.75) (`:70`);
  - `_box` (`:106–134`, tags `TAG_TERM` at x = length, cells allreduced);
  - `_solve` (`:137–168`; terminations `abc` / `exact` / `pec` / `pmc` only, **no extra-term argument**);
  - `_assemble_complex` (`:185–187`);
  - `_Projection` (`:190–240`): a volume LSQ onto e^{∓jβx} over a DG0 window set in numpy. `k_fit` (`:233–240`) is a 21-point residual scan at ±5e-3 relative with a parabolic vertex, and returns `(k, Γ(k), interior)`;
  - `_tem_box` (W = 2h, 2 × 2 cells across; `:247–250`) and `_tem(comm, h, degree, termination, scan)` (`:253–265`);
  - the ABC ladder test (`:317–364`) and the TE10 test (`:367–393`).
- **Step 1a's record** (`docs/testing/logs/20260930T094044Z_TH-5-step1a2-W1.log.gz`; 16 passed, 82 s, `-n 2`, complex):
  - degree 1, projection at k₀: \|Γ\| 3.814078e-3 / **9.499917e-4** / 2.373045e-4 / 5.931735e-5 at λ₀/20 / 40 / 80 / 160 (`:105`, **`:107`**, `:109`, `:111`), rate 2.0021 (`:112`);
  - degree 2 at λ₀/40: 1.368981e-7 (`:114`);
  - PEC control at λ₀/80: Γ = −1.000000 + 0.000538j (`:115`);
  - **k_fit at λ₀/80:** ABC k_fit − k₀ = +5.5894e-4 rad/m, \|Γ\| at k_fit **2.126188e-4** (`:118`); PEC +5.5789e-4 (`:119`).
- **Solver hook (not used here; the term has its shape):** `TimeHarmonicSolver.solve(extra_bilinear_terms=…)` (`src/fem_em_solver/core/time_harmonic.py:400`, `:556–559`).

**The floor of step 1a's read at λ₀/40 (rubric element 2: the ceiling, here a floor, is computed before the bar is applied).** On the TEM line at normal incidence the first-order ABC is *exact* in the continuum: R(0) = 0 (Jin (9.48); `abc_first_order_reflection_te(1.0) == 0`, `test_th5_abc_reflection.py:87`). So the ABC's 9.5e-4 at λ₀/40 is discretisation only, and the comparison the weekly asks for is discrete against discrete.
- The fixture's discrete wavenumber is k_h ≠ k₀. The record measures k_fit − k₀ = +5.59e-4 rad/m at λ₀/80 on both ABC and PEC solves (`:118–119`); h² scaling gives ≈ +2.24e-3 rad/m at λ₀/40.
- **A perfectly reflectionless forward wave e^{−jk_h x}, projected at k₀ over `PROJ_N`, reads \|Γ\| = 1.18e-4–1.20e-4 at λ₀/40.** This is a drafter's stdlib computation (piecewise-linear interpolant, six-point Gauss per cell, 1a's window), unlogged, a reading and not a window. It is the read's own dispersion leak. It does not depend on the termination.
- **Cross-check against the record:** the same computation gives 2.97e-5 at λ₀/80, against the record's measured gap between the k₀ and k_fit reads of the ABC, 2.373045e-4 − 2.126188e-4 = 2.47e-5 (`:109`, `:118`; a vector leak bounds the magnitude gap from above).
- **So at the k₀ read the weekly's ≤ 1e-4 at λ₀/40 is below the read's floor (≈ 1.2e-4), whatever terminates the line.** Including a PML's own discrete reflection, the predicted band at the k₀ read is [0.9e-4, 1.6e-4]. The bar sits inside it, so a green there would be cancellation luck, not absorption, and a red would not be a finding about the PML. This band straddles the bar, so it cannot be written as a `PREDICTED RED` range either.
- **At the k_fit read** (1a's own `_Projection.k_fit`, recorded at λ₀/80) the leak is removed. The predicted PML reading is **≈ 5e-6**, bracket [2e-6, 5e-5]. **The bar is reachable with ≈ 2–50× margin (point ≈ 19×).**

**The discretisation estimate (drafter's stdlib 1D model, unlogged; predictions only).**
- **The model.** On the TEM line, E = E_z(x) ẑ, and Λ reduces to −(u′/s)′ − k₀² s u = 0: the yy entry of Λ⁻¹ and the zz entry of Λ, since ∇×(E_z ẑ) = −E_z′ ŷ. The model discretises this with 1D P1 elements, the mass a blend α · consistent + (1 − α) · lumped.
- **Calibration of α:** α = 0.674 reproduces the record's dispersion (k_h − k₀ 5.61e-4 against 5.59e-4 at λ₀/80).
  - It reproduces the record's **ABC k_fit read at λ₀/80 to 0.15 %** (2.123e-4 against 2.126188e-4, `:118`).
  - It sits 1.25× below the record's k₀ reads (7.57e-4 against 9.50e-4 at λ₀/40), because the leak's phase differs.
  - α = 0.5 (no dispersion) brackets it from the other side (ABC 1.03e-3 at λ₀/40).
- **Design scan, PEC-backed, m ∈ {1, 2, 3}, d ∈ {0.5, 1, 1.5, 2} λ₀, design R ∈ {1e-3 … 1e-6}:**
  - d = 0.5 λ₀ (20 cells) floors at ≥ 4e-4 at λ₀/40, whatever s″_max is.
  - d = 1.0 λ₀ reaches 3.3e-5–4.6e-5 at the k_fit read (R = 1e-6, m = 2), only 2–3× under the bar.
  - **d = 1.5 λ₀, m = 2 reaches 5.2e-6 (α 0.674) / 1.4e-5 (α 0.5).** These are the item's design constants.

**`TH-5` step 3 — a polynomial-graded PML (Jin's anisotropic absorber) on step 1a's TEM line, PEC-backed, read by step 1a's projection; its reflection against the graded closed form at λ₀/40 at the weekly's 1e-4, with the first-order ABC on the same line and the same read as the negative control.**
- **Executor:** implementer.
- **Tier:** standard (W1 predicted ≈ 100 s), complex build, `-n 2` (1a's width; the read's assemblies are allreduced).
- **Class:** lane-2 feature step. **Cost:** ≈ 10.6 measured (= 46 predicted × 0.231). The 10-21 weekly's supply sheet priced it at 16.7 (`PROJECT_PLAN.md:5014`).
  - The ≈ 46 predicted minutes are: 15 fixed, ≈ 15 code (≈ 60 lines in `boundary_terms.py` and a ≈ 220-line module), W1 ≈ 2, W2 ≈ 1, docs (§7 row sentence, §10, outcome line) ≈ 8, and ≈ 5 spare.
  - The window prices come from 1a's W1: 82 s for all of `tests/environment` plus the 1a module at `-n 2` (`20260930T094044Z_TH-5-step1a2-W1.log.gz`, footer). Dropping the other three environment files saves ≈ 35–45 s (2c's measurement, this file `:315`).
  - The PML meshes are 1 152 / 2 304 / 4 608 / 9 216 cells at λ₀/20 … /160 (length 2.4 m × 24 cells per x-step), each solve well under 1 s. The cost is the one complex JIT of the PML form.
- **Depends on:** nothing. **Daytime:** no. **Slot-first:** no. **Attempt:** 1 of 4.

**Code.** No edit to any existing function in `boundary_terms.py`, to `core/__init__.py`, `time_harmonic.py`, `study.py`, or any `TH-5` test module (1a, 1b, 2a, 2b).
- **`src/fem_em_solver/core/boundary_terms.py`, additive only (appended below `abc_first_order_reflection_te`; the module docstring gains one line naming step 3):**
  - `pml_slab_term(msh, x_start, thickness, s_max, order, omega_rad_s)` returns `f(trial, test)`, the *correction* that turns the vacuum form into Jin's anisotropic-absorber PML (§9.6.2 (9.152)–(9.154), §9.6.3 (9.157)) in the slab x_start < x < x_start + thickness. The slab's normal is x̂, outgoing +x.
    - The correction is ∫ χ [(Λ⁻¹ − I)∇×u·∇×v̄ − k₀²(Λ − I)u·v̄] dx.
    - Λ = diag(1/s, s, s) and Λ⁻¹ = diag(s, 1/s, 1/s), built with `ufl.as_matrix`.
    - s = 1 − j s_max χ ξ^m, ξ = (x[0] − x_start)/thickness (the (9.165) profile).
    - χ is a DG0 indicator **set in numpy at cell midpoints** (1a's pattern, `test_th5_abc_reflection.py:201–203`), never a UFL comparison.
    - The measure has `metadata={"quadrature_degree": 6}` pinned (1a's `_Projection` degree, `:209`).
    - The tensor sits on the trial side, and the test function is conjugated by `ufl.inner` (the docstring convention, `:3–10`).
    - It is vacuum only, like `first_order_abc_term`. The docstring states that `s_max < 0` is the conjugated, amplifying stretch (the γ trap's analogue, `:10–11`), kept callable for the negative control.
  - `pml_reflection_pec_backed(s_max, thickness, order, k0, cos_theta=1.0)` returns the complex Γ_c = −exp(−2jk₀cosθ d) · exp(−2k₀cosθ s_max d/(m + 1)) above. Its docstring cites (9.133) and (9.165), says the graded integral is the drafter's derivation and not a numbered equation, and gives the σ form.
  - `pml_s_max_for_reflection(r, thickness, order, k0, cos_theta=1.0)` = −ln r · (m + 1)/(2k₀cosθ d).
- **New module `tests/validation/test_th5_pml_reflection.py`** (tests only).
  - It imports from `tests.validation.test_th5_abc_reflection`: `L_N`, `PROJ_N`, `K0`, `OMEGA`, `LAMBDA0`, the tags, `_box`, `_Projection`, `_assemble_complex` and `_tem` (1b's precedent, this file `:96–98`), plus `IDENTITY_TOLERANCE` (`tests/validation/test_port_reaction_impedance.py:192`, 1e-9).
  - **Constants at module top, before any solve:**
    - `PML_REFLECTION_BAR = 1e-4`, with a comment quoting the weekly (`PROJECT_PLAN.md:4829`);
    - `PML_D = 1.5` (λ₀; 60 cells at λ₀/40), `PML_M = 2`;
    - `R_ABSORB = 1e-6` (leg A, s″_max = 2.19880) and `R_CLOSED = 1e-3` (leg C, s″_max = 1.09940), both through `pml_s_max_for_reflection`.
    - d = 1.5 λ₀ makes e^{−2jk₀d} = 1, so Γ_c = −R.
  - **Own solve, copied from 1a's `_solve` (`:139–168`, the `pec` path):** `_solve_terms(msh, tags, degree, drive, pec_tags, extra)` adds `sum(t(u, v) for t in extra)` to `a`, with PEC on `TAG_Z` and `TAG_TERM`.
  - **`_pml_tem(comm, h, s_max, scan=True)`:**
    - builds `_box(comm, L_N + PML_D, 2h, 2h, round((L_N + PML_D)/h), 2, 2)`, with node planes at x = 0.9 and 2.4 at every rung (0.9/h = 18 / 36 / 72 / 144);
    - solves with `pml_slab_term(msh, L_N, PML_D, s_max, PML_M, OMEGA)`;
    - reads `_Projection(e_h[2], msh, *PROJ_N, K0)` at k₀ (`.gamma(K0, L_N)`) and at `.k_fit(K0, L_N)`;
    - returns both, plus the cell count and seconds.
  - `test_pml_closed_form_numpy`: ungated, numpy only.
  - `test_pml_tem_line`: `@complex_only` (1a's marker, copied), `@pytest.mark.integration`; one solve set, every assertion and print.

**Anchors (asserted).**
- **(Z) numpy identities, each ≤ `IDENTITY_TOLERANCE`:**
  - \|Γ_c\| = R for both designs;
  - `pml_s_max_for_reflection` inverts `pml_reflection_pec_backed`;
  - the σ form equals the s″ form;
  - s_max = 0 gives \|Γ_c\| = 1;
  - \|Γ_c(θ)\| = \|Γ_c(0)\|^{cosθ}.
  - Predicted ≤ 1e-15.
- **(V) premises at every rung:**
  - ∫χ_slab dV = PML_D·(2h)² and 1a's window ∫χ dV = 0.60·(2h)², each relative ≤ `IDENTITY_TOLERANCE`. 1a read ≤ 2.7e-13 on the window (`:94–101`); predicted ≤ 1e-12.
  - The cell count is 24 · round(2.4/h).
- **(K) every k_fit read that an assertion uses has `interior` True.** k_fit − k₀ is predicted ≈ +2.24e-3 rad/m at λ₀/40, inside the scan's central step (3.1e-2 rad/m half-width, 3.1e-3 rad/m step).
- **(A) absorption leg, λ₀/40, k_fit read:** \|Γ − Γ_c\| ≤ `PML_REFLECTION_BAR` at R_ABSORB. Predicted 5.2e-6 (α 0.674) / 1.4e-5 (α 0.5), bracket [2e-6, 5e-5].
- **(C) closed-form leg, λ₀/40, k_fit read:** \|Γ − Γ_c\| ≤ `PML_REFLECTION_BAR` at R_CLOSED, where Γ_c = −1e-3, so the bar is 10 % of the closed form. Predicted \|Γ\| ≈ 1.005e-3 and \|Γ − Γ_c\| ≈ 5.5e-6 (both models), bracket [2e-6, 3e-5]. **This leg is what makes the comparison a test of the closed form.** At R_ABSORB, Γ_c is invisible under the discretisation term.
- **The 1a module runs green in the same window** (its three tests, unedited). This is the regression for the additive `src/` edit; rule (c) is satisfied in-slot.

**Negative controls (rule (e)).**
- **(NC-ABC) the ABC on the same line and read, λ₀/40: asserted.** The run uses `_tem(comm, LAMBDA0/40, 1, "abc", scan=True)` from 1a.
  - **At the k₀ read, asserted > `PML_REFLECTION_BAR`, backed:** the record 9.499917e-4 (`:107`, the same code and fixture). The ceiling is 9.5× the bar.
  - **At the k_fit read, asserted > `PML_REFLECTION_BAR`; the floor is the weekly's bar and the factor is predicted:** ≈ 8.5e-4 (model 8.498e-4), backed in kind by the λ₀/80 record 2.126188e-4 > 1e-4 (`:118`). The ceiling is ≈ 8.5×.
  - **Printed:** the k₀ read against 9.499917e-4 (predicted to reproduce to every printed digit: same code, image and `-n 2`), and \|Γ_ABC\| / \|Γ_PML\| at k_fit (predicted ≈ 160, 74–180 across the models).
  - **"The ABC's own ladder"** is 1a's ladder test, re-run green in the same window: ≤ 1e-4 is reached only at λ₀/160 (5.93e-5, `:111`), where the PML is asserted at λ₀/40.
- **(NC-0) s″_max = 0 (the term enters but does not absorb; also catches a dropped term): asserted.** \|Γ\| > `PML_REFLECTION_BAR` at λ₀/40, with the floor imported and the factor predicted. Predicted \|Γ\| = 1.000000: vacuum to a PEC at 2.4 m, the 1a PEC control's analogue, 1.000000 at λ₀/80 (`:115`).
- **Printed only (new comparisons):**
  - **(NC-conj)** s_max = −1.09940 (the conjugated stretch): \|Γ\| predicted ≈ 9.9e2–1.0e3, i.e. 1/R_CLOSED: the layer amplifies.
  - **(NC-grade)** leg C against the un-graded (9.133) exponent (R = e^{−2k₀s″d} = 1e-9): miss predicted 1.0e-3, ceiling exactly \|Γ_c\|/bar = 10.
  - Leg C against the m = 3 formula applied to the m = 2 profile (R = 5.62e-3): miss predicted 4.6e-3.
  - Together these show the anchor resolves the 1/(m + 1) grading factor.

**Predictions (written first; printed, never asserted).** The k_fit read is the α 0.674 model, the k₀ read the same model's:
- **Leg A ladder at λ₀/20 / 40 / 80 / 160, k_fit read:** 1.9e-5 / 5.2e-6 / 1.7e-6 / 1.0e-6. It flattens at /160 onto the design's continuum 1e-6; the LSQ rate over /20–/80 is ≈ 1.7–1.9, printed. At the k₀ read: 4.8e-4 / 1.16e-4 / 2.8e-5 / 6.5e-6.
- **Leg C at λ₀/40, k₀ read:** \|Γ − Γ_c\| ≈ 1.16e-4, the floor above.
- **ABC at k_fit:** k_fit − k₀ ≈ +2.24e-3 rad/m at λ₀/40, and \|Γ\| 3.41e-3 / 8.50e-4 / 2.12e-4 / 5.31e-5 at /20 … /160 (only /40 is run here).
- **The 3D/1D level ratio is unmeasured for the PML.** The ×3 brackets above cover it, and the α 0.5 model is the upper bracket.

**Windows.** Each runs in the foreground, with Bash `timeout: 660000`, `-s`, node paths only (never `-k`), and no pipe or redirect. Both are under 10 min, so no durable capture. **The harness call is the whole Bash command.**
- **W1** (standard, `-n 2`; predicted ≈ 100 s, bracket 60–170 s):
  `scripts/testing/run_and_log.sh TH-5-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_th5_abc_reflection.py tests/validation/test_th5_pml_reflection.py -v -s --tb=short'"`
- **W2** (real build: collection, skips and the numpy test; ≈ 15–30 s):
  `scripts/testing/run_and_log.sh TH-5-step3-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_th5_abc_reflection.py tests/validation/test_th5_pml_reflection.py -v -s --tb=short'"`
- Before W1, put `git diff --stat src/fem_em_solver/core/boundary_terms.py` in the journal: insertions only.

**Traps already paid for.**
- **Never assert at the k₀ read.** Its floor (≈ 1.2e-4 at λ₀/40) is above the bar. A k₀-read green comes from leak cancellation: the model's d = 1.0, m = 3, R = 1e-4 cell reads 1.8e-5 at k₀ by phase luck. **Design constants are fixed at module top before W1 and never tuned between windows.**
- **No UFL comparison anywhere in the slab.** χ is DG0 and set in numpy (no `ufl.conditional`, `ufl.gt`, `ufl.max_value`, or `sqrt` in a comparison: `OPS-22`, `WF-6` 4c). ξ^m with integer m on a negative ξ outside the slab is harmless only because χ = 0 multiplies it. Keep χ inside s, not outside the form.
- **Pin the quadrature degree.** The integrand carries the rational 1/s and a `SpatialCoordinate` (`POST-5`'s compile trap, by analogy).
- **Sign convention:** s = 1 − js″ absorbs under e^{+jωt}, and the conjugate amplifies (the γ trap, `boundary_terms.py:10–11`). Never re-derive it; (NC-conj) prints it.
- **Node planes:** x = 0.9 and 2.4 must be mesh planes at every rung, or χ mis-assigns a straddling cell. Do not change `PML_D` to a value that breaks 0.9/h, 2.4/h ∈ ℤ.
- **`_solve` cannot take the term.** Copy its lines (`:139–168`). Never add a kwarg to 1a's module: rule (i), and it is a gated module.
- **Rank-safety:** `_Projection` and `_box` allreduce already. Take any new scalar (∫χ_slab dV) through `_assemble_complex`. Print on rank 0.
- **A kill at 180 s is not a result.** Run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first, and the rubric's one-call orphan-`.c` sweep only if it reads 0. Never lengthen the timeout or change `-n`.
- **Rule (i):** no edit to either file after W1 starts. Write ≈, never `~`, in tracked Markdown.

**Scope.** Covered: one fixture (1a's TEM line, normal incidence, degree 1, vacuum, `-n 2`), one PML design family (d = 1.5 λ₀, m = 2, PEC-backed), and the k_fit read. **Not claimed:**
- oblique incidence: the TE10-60° guide, where the ABC's continuum −1/3 is what a PML would actually beat (`:125`);
- degree 2: the ABC already reads 1.4e-7 at λ₀/40 (`:114`);
- an ABC-backed or complementary PML (Jin §9.6.5–9.6.6);
- curved or corner PMLs;
- any use through `TimeHarmonicSolver`, `Study` or YAML;
- the 64 / 128 MHz birdcage box (the entry's step-3 question, `PROJECT_PLAN.md:2464–2466`);
- any HFSS *Radiation* or PML parity statement;
- the §2.1 line's "What this is not: … PML" clause (`:396–399`), which stays until the review re-words it.

**Negative result.** In every case: report the measurement, park on `attempt/TH-5-step3-*`, mark item 596 BLOCKED in the same commit (rule (d)), and journal it in attempts.md. No band, design constant or read moves in-slot.
- **(Z) or (V) red:** a code defect in the closed form or the slab indicator. Stop.
- **(K) red:** the k_fit vertex sits at the scan edge, and that read is invalid. Report k_fit − k₀; never widen `half_width` in-slot.
- **(NC-0) red** (\|Γ\| ≤ bar with s″ = 0): the term absorbs when it should not, a code defect. Stop.
- **(NC-ABC) red at the k₀ read:** 1a's record did not reproduce, so something on `main` changed under the fixture. This is a finding for the review.
- **(A) or (C) red with every control green:** the degree-1 PML's discrete reflection at λ₀/40 exceeds the bar on this design. The measured ladder goes into the §7 row as the step-3 record, and the weekly decides (a thicker layer is a new design, not an in-slot change).
- **(C) red with (A) green:** the grading or the closed form is mis-implemented. Read the printed (NC-grade) misses and report.

No `PREDICTED RED:` line: with the k_fit read every asserted clause is predicted green. At the k₀ read no valid range exists, because the floor's band [0.9e-4, 1.6e-4] straddles the bar.

**Status it can move.** No glyph. `TH-5` is ✅, and its entry Done-when ("Steps 1–2 executed and asserted; §2.1 gains the boundary-model line") was met at 2c (`PROJECT_PLAN.md:2467–2473`, row `:2216`); step 3 is the entry's optional step, so rule (j) binds nothing. On W1–W2 green:
- the §7 row gains a step-3 record sentence (legs A / C, the ladder, the ABC control at both reads, elapsed time, logs);
- the §10 lane-2 entry marks step 3 DONE, and **Tier B is complete** (ladder row B1, `:3804`, gains the record);
- an outcome line goes under this brief.

The §2.1 boundary line's "a second-order ABC or PML" exclusion is the review's to re-word, not the implementer's.

**For the review.**
1. **The read (the item's crux).** At step 1a's k₀ read the weekly's ≤ 1e-4 at λ₀/40 is **below that read's floor (≈ 1.2e-4)**, whatever the termination. Same lesson as step 2: on 2026-09-30 the 10:30 review kept step 2 out of the queue because its ≥ 5× bar sat above the sphere model's ≈ 2× ceiling (`PROJECT_PLAN.md:2524–2535`). 2a then measured the banked negative (r(W) 0.4846 → 0.2266, this file `:214`), and 2c closed on the re-pose (10-08). Here the analogue is a floor, and the fix is in the reading, not the physics. The drafter's default asserts at 1a's existing k_fit read and prints the k₀ read beside it.
   - **Whether that is the item's measurement definition (the review's) or a change to the weekly's target (back to the 10-24 weekly) is the ruling asked for.**
   - Alternatives if the k₀ read binds: λ₀/80, where the floor is ≈ 3e-5 but the PML's k₀ read is predicted at ≈ 2.8e-5, so marginal again; or degree 2. Both re-pose the weekly's "at λ₀/40".
2. **Citation.** The converted Jin carries the equation *numbers* but not their bodies (`INDEX.md:3–6`).
   - The graded closed form is (9.133)'s exponent integrated over (9.165)'s profile. It is a two-line derivation and **not a numbered equation**.
   - "δ_max ≡ s″_max with s_z = 1 − jδ_max(z/L)^m" is read from the prose.
   - A check against the print copy is owed (the operator's, if the review wants it before W1). The item does not depend on Jin's exact notation, because (Z) asserts the derived form's own identities.
3. **Design constants:** d = 1.5 λ₀, m = 2, legs at R = 1e-6 and 1e-3. d = 1.0 λ₀ has only 2–3× margin in the model.
4. **Leg C.** The drafter reads "its reflection against the graded-PML closed form, ≤ 1e-4" as \|Γ − Γ_c\| ≤ 1e-4, and adds the R = 1e-3 leg so that the closed form is visible. Leg A alone only shows "small". Ratify, or drop leg C.
5. **What the fixture can and cannot show.** At normal incidence the first-order ABC is exact in the continuum, so this compares discretisation errors (a PML at λ₀/40 against an ABC that needs λ₀/160). The PML's continuum advantage, oblique incidence, is untested. A TE10-60° PML leg on 1a's guide would be the natural next step; that is the weekly's to open.
6. **Price.** By the lane-2 rule (× 0.231) this is ≈ 10.6 measured. The weekly's supply sheet used 16.7.

**Execution order.** Read the brief and any review ruling under it → code (`boundary_terms.py` additions, new module) → journal the `git diff --stat` → W1 → W2 → commit (both files, both harness logs, the §7 row sentence, the §10 lane-2 mark, the outcome line here) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 596 — `TH-5` step 3 (§10 lane 2, the PML) — per `docs/planning/chunks/TH-5.md` §9 item 596. Read the brief and any review ruling beneath it first; the rulings hold, in particular on which read the 1e-4 bar is asserted at. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run W1 then W2, in the foreground, with `timeout: 660000` on every Bash call that runs a window, pytest with -s, container `timeout -k 30`, never `-k` on pytest. Add only the three functions named in the brief to `src/fem_em_solver/core/boundary_terms.py`, appended, with no existing line changed and nothing exported from `core/__init__.py`; create only `tests/validation/test_th5_pml_reflection.py`. Do not edit `tests/validation/test_th5_abc_reflection.py`, any other `TH-5` module, `time_harmonic.py` or `study.py`. Fix the design constants (`PML_D` = 1.5, `PML_M` = 2, `R_ABSORB` = 1e-6, `R_CLOSED` = 1e-3, `PML_REFLECTION_BAR` = 1e-4) before W1 and never change them between windows. Never assert at the k₀ read unless the review's ruling says so. No edit to either file after W1 starts. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. After any killed window, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first. On a red, park on `attempt/TH-5-step3-*`, mark item 596 BLOCKED in the same commit, and report."

### Review ruling 2026-10-22 03:00 (item 596, step 3)

Binds over the brief. Queued with the reflection **asserted at step 1a's own k_fit read** (the residual-minimising read the module already carries), with the k₀ read printed beside it: the drafter's floor arithmetic (k₀ read ≈ 1.2e-4 for a reflectionless wave at λ₀/40, checked against the λ₀/80 record, 2.97e-5 predicted vs 2.47e-5) shows the k₀ read cannot resolve 1e-4. The bar (≤ 1e-4 at λ₀/40), the fixture and the ABC control are the weekly's, unmoved; **the read is a provisional execution reading the 10-24 weekly ratifies or reverts** (rule (h): the same physics, the PML's reflection, with the fixture's dispersion leak removed). (1) The Jin citation — §9.6.1 eq. (9.133) integrated over the §9.6.4 (9.165) profile — is a derivation; the print-copy check is the operator's (the converted text's equations are image placeholders). (2) Design constants (d = 1.5 λ₀, m = 2, legs at 1e-6 / 1e-3) stand under rule (f), provisional. (3) Normal incidence only: no oblique-incidence claim.

**Outcome 2026-10-22 (20:00 slot, attempt 1 of 4): DONE, every assertion green.** W1 `20261022T203122Z_TH-5-step3-W1.log` (10 passed, 47 s, `-n 2`, complex); W2 `20261022T203218Z_TH-5-step3-W2.log` (3 passed / 7 skipped, real, `-n 1`). `boundary_terms.py` +63 insertions, 0 deletions. (Z) worst identity 1.332e-15; (V) slab / window volumes ≤ 2.93e-13 rel; (K) k_fit − k₀ +2.2451e-3 rad/m at λ₀/40, interior. (A) R 1e-6: k_fit read \|Γ − Γ_c\| 1.6680e-5 / **1.3859e-6** / 1.6651e-7 / 2.4307e-8 at λ₀/20 / 40 / 80 / 160, \|Γ\| flattening at ≈ 1.0e-6 (the design value); LSQ rate /20–/80 3.323 (predicted 1.7–1.9: faster than the 1D model, recorded, not asserted). (C) R 1e-3: Γ = −9.985820e-4 + 7.597845e-6j, \|Γ − Γ_c\| **7.7290e-6** (predicted ≈ 5.5e-6). NC-grade misses 9.9861e-4 / 4.6248e-3 (predicted 1.0e-3 / 4.6e-3). NC-ABC k₀ read 9.499917e-4 = record to every digit, k_fit read 8.5100e-4 (predicted 8.50e-4), ratio 784.5 (predicted ≈ 160: the PML reads ≈ 5× below the model). NC-0 \|Γ\| 1.000000; NC-conj 9.1059e+2. PML k₀ reads at λ₀/40: leg A 1.1821e-4, leg C 1.2434e-4 — the read's floor as predicted (≈ 1.2e-4), printed only.

### §9 item 628 — 2026-10-25 (§10 lane 2, `TH-5` step 3b: the PML at oblique incidence)

*Written by the 2026-10-25 03:00 daily review from the 10-24 weekly's lane-2 target (§9 rule (5); `PROJECT_PLAN.md` § 10, "Lane 1 and lane 2"). Two execution readings below go to the 10-26 weekly.*

**Question.** Does step 3's graded, PEC-backed PML reflect a TE10 mode at 60° by its closed form R(θ) = R(0)^{cos θ} (magnitude), within step 3's absolute band at λ₀/40?

**Execution readings (this review's):**
1. **Degree 2, not 1.** On step 1a's TE10 guide at degree 1, the matched termination's own point-fit floor already reads 1.16e-4 at λ₀/80 (`EX-143`'s (O1), `20261024T003217Z_EX-143.log:55`). The ABC's discretisation miss at λ₀/40 is 1.03e-3 at degree 1 against **4.12e-6 at degree 2** (`20260930T094044Z_TH-5-step1a2-W1.log.gz:127–128`). A 1e-4 bar at degree 1 sits at or below the read's floor, which would make the target unassertable. Degree 2 is the production default (order rule, 2026-10-01), and the weekly's target names no degree. So the gate is asserted **at degree 2**, with degree 1 printed beside it.
2. **The read.** Step 3 asserted at 1a's `_Projection.k_fit` read (ratified 10-24). The TE10 helper `_te10` has no k_fit read (`test_th5_abc_reflection.py:269`); the k_fit scan exists only in `_tem` (`:254`, `_Projection.k_fit` at `:233–240`). The new module adds the same scan on the TE10 projection (β scanned ±5e-3 relative about k₀ cos θ, 21 points, weighted by sin(πy/a)), imported if `_Projection` takes the weight. If it does not, the scan is copied with the source lines cited. The k₀ cos θ read is printed beside it.

**Anchors.**
- **(B)** \|Γ_h − Γ_c(θ)\| ≤ `PML_REFLECTION_BAR` (1e-4, `tests/validation/test_th5_pml_reflection.py:55`) at λ₀/40, degree 2, θ = 60°, R(0) = `R_ABSORB` = 1e-6 (s_max from `pml_s_max_for_reflection` at normal incidence, `S_ABSORB` 2.19881).
- Γ_c(θ) comes from `pml_reflection_pec_backed(…, cos_theta=0.5)` (`src/fem_em_solver/core/boundary_terms.py:107`): −exp(−2jk₀ cos θ · d) · exp(−2k₀ cos θ · s_max · d/(m + 1)). At d = `PML_D` = 1.5 λ₀ the phase factor is e^{−3πj} = −1, so **Γ_c = +1.000e-3**, with magnitude R(0)^{cos θ} = 1e-3. Print both.
- **(Z)** \|Γ_c(θ)\| = \|Γ_c(0)\|^{cos θ}, already asserted in the numpy test (`:86`); re-run, not re-written.
- **Closed form:** Jin 3e §9.6.1, eq. (9.133) (the PEC-backed reflection, attenuation exp(κ₂ s″ z cos θ), `docs/references/jin-fem-3e/09_absorbing-boundary-conditions.md:877–879`) over the polynomial profile (9.165) (§9.6.4, `:1041`). The graded form R(0)^{cos θ} is (9.133) integrated over (9.165), not a numbered equation, and the print-copy check is the operator's (dashboard item 0).

**Control (asserted, backed).** **(Cabc)** On the same guide, mesh and degree, the first-order ABC termination's \|Γ\| is ≈ 1/3 (Jin (9.48); `EX-143`'s (O2) −3.332454e-01). Assert \|Γ_ABC\| / \|Γ_PML\| > 100. Ceiling: (1/3) / 1e-3 = 333, so 100 is reachable with margin 3.3×, and the ratio is printed. **(NC-0)** s_max = 0 (the PML switched off, PEC at the end) reads \|Γ\| ≈ 1, as step 3's NC-0 does (`20261022T203122Z_TH-5-step3-W1.log:129`). Assert \|Γ\| > 0.9.

**Route.** New module `tests/validation/test_th5_pml_oblique.py` (the gated `test_th5_abc_reflection.py` and `test_th5_pml_reflection.py` are not edited).
- **Mesh:** step 1a's `_te10` geometry (`A_O` = λ₀/(2 sin 60°), `L_O` = 1.5, PEC on `TAG_Y` / `TAG_Z`, the drive sin(πy/a)·ẑ; `:75–80`, `:269–277`), extended by `PML_D` in x as `_pml_tem` extends the TEM line (`test_th5_pml_reflection.py:128–135`), PEC backing at the far end. The PML term is `pml_slab_term` (`boundary_terms.py:68`), `PML_M` = 2, χ DG0, quadrature degree 6.
- **Solve:** step 3's `_solve_terms` (`test_th5_pml_reflection.py:98`), imported.
- **Fit window:** `FIT_O` = (0.25, 1.25) (`:79`), inside the physical guide.

**Windows.**
1. W0: the numpy (Z) node, ≈ 5 s.
2. W1 `TH-5-step3b-W1`: degree 2 at λ₀/40 with PML, the ABC control, NC-0, and degree 1 printed; `-n 2`, complex, standard, `timeout -k 30 180`.

**Price:** step 3's TEM line at λ₀/40 solved in 1.0 s (2 304 cells). This guide at λ₀/40 is ≈ 120 × 23 × 2 hexes ≈ 33 k tets; degree 2 is ≈ 2 × 10⁵ dofs, a few seconds per solve. With compiles, W1 ≈ 40–90 s.

**Tier:** standard, `-n 2`. **Class:** lane-2 feature. **Cost:** ≈ 27.9 measured / 40 predicted slot-minutes. **Executor:** implementer. **Depends on:** nothing. **Family:** `TH-5` step 3b, attempt 1 of 4.

**Traps already paid for.**
- `ufl` comparisons on complex operands. The χ indicator is DG0 from coordinates, built in numpy: never `ufl.conditional` on a complex `SpatialCoordinate` expression, and no `sqrt` in a comparison (rubric 4).
- Pin `quadrature_degree` on every facet / volume form that carries the PML profile.
- The fit is collective: reduce before the LSQ.
- Point-sample bytes depend on the batch: compare call by call (item 537).
- `-k 30`; never background a window; orphan `.c` sweep only with `pgrep -c python3` reading 0.

**Scope boundary.** One angle, one thickness, one profile. Not closed: other angles or a sweep in θ, evanescent modes, the PML at the birdcage's box (k₀R below π, unchanged), a TH-5 glyph change. `TH-5` stays ✅ on step 2c; step 3b is a record on the optional step, like step 3.

**A negative result means:**
- (B) red at degree 2 with (Cabc) and NC-0 green: a measured negative on the oblique PML. Record it in `TH-5.md` with the degree-1 / degree-2 pair and the k₀ / k_fit reads, and park. Never loosen `PML_REFLECTION_BAR`.
- (Cabc) not separating: the read is broken. Stop.

**Status it can move.** §10 lane 2 `TH-5` 3b → DONE. `TH-5`'s row gains "step 3b: oblique incidence (TE10 at 60°), degree 2" as a record. Dashboard item 4's open `TH-5` question closes.

**Outcome 2026-10-25 (06:00 slot, attempt 1 of 4): DONE, every assertion green.** New module `tests/validation/test_th5_pml_oblique.py` only (no edit to `boundary_terms.py`, `test_th5_abc_reflection.py` or `test_th5_pml_reflection.py`; `_Projection` takes the weight, so the scan is imported, not copied). W0 `20261025T111701Z_TH-5-step3b-W0.log` (step 3's numpy (Z) node, 1 passed, worst identity 1.332e-15, `:50`, 2 s). W1 `20261025T111709Z_TH-5-step3b-W1.log` (1 passed, 131 s, `-n 2`, complex; predicted 40–90 s, the degree-2 solve took 66.4 s). Γ_c(60°) = 1.000000e-03 − 3.2e-18j, \|Γ_c\|/\|Γ_c(0)\|^{cos θ} − 1 = 1.110e-15 (`:52`). **(B)** degree 2, λ₀/40 (33 120 cells): k_fit read Γ = 9.999865e-04 + 1.094756e-07j, **\|Γ − Γ_c\| 1.1031e-07** (bar 1e-4), k_fit − β −2.5009e-06 rad/m, interior; k₀cos θ read \|Γ − Γ_c\| 2.8512e-07 (`:53`). Degree 1 (printed): k_fit read \|Γ − Γ_c\| 6.6733e-05, k₀cos θ read 1.8789e-04, k_fit − β +1.5765e-03 (`:54`) — the degree-1 k₀ read would have missed the bar, as the review's reading 1 predicted. **(Cabc)** degree-2 ABC Γ = −3.333334e-01 (Jin (9.48): −1/3), \|Γ_ABC\|/\|Γ_PML\| **333.3** (bar 100, ceiling 333) (`:59`). **(NC-0)** s_max = 0: \|Γ\| 1.000000 (`:60`). The ABC control ran on 1a's own L_O guide (same cross-section, h and degree; no PML extension).

## §7 row moved verbatim 2026-10-24 (rotate_plan_archive.py chunks)

 Absorbing boundary condition (ABC) | ✅ **Step 3b DONE 2026-10-25 (§9 item 628, attempt 1 of 4; record on the optional step): oblique incidence (TE10 at 60°), degree 2** — step 3's PML (m = 2, d = 1.5 λ₀, s_max from R(0) = 1e-6) on 1a's TE10 guide, read at 1a's weighted projection's k_fit: (B) \|Γ − Γ_c(60°)\| 1.1031e-07 against Γ_c = +1e-3 = R(0)^{cos θ} (bar 1e-4); degree 1 printed 6.6733e-05 (k₀cos θ read 1.8789e-04); (Cabc) ABC Γ = −0.3333334, ratio 333.3; NC-0 \|Γ\| 1.000000. W1 1 passed, 131 s (`20261025T111709Z_TH-5-step3b-W1.log`). One angle, one thickness, one profile. **Step 3 DONE 2026-10-22 (§9 item 596, attempt 1 of 4; the PML, record on the optional step):** a polynomial-graded (m = 2, d = 1.5 λ₀) PEC-backed PML (`pml_slab_term`, additive in `boundary_terms.py`) on 1a's TEM line, degree 1, `-n 2`, asserted at 1a's k_fit read (03:00 ruling, provisional; the 10-24 weekly ratifies): (A) R = 1e-6, \|Γ − Γ_c\| 1.3859e-6 at λ₀/40 (bar 1e-4; ladder 1.6680e-5 / 1.3859e-6 / 1.6651e-7 / 2.4307e-8 at /20…/160, rate /20–/80 3.323); (C) R = 1e-3, Γ = −9.985820e-4 + 7.597845e-6j against Γ_c = −1e-3, \|Γ − Γ_c\| 7.7290e-6; controls: ABC k₀ read 9.499917e-4 (record, every digit), k_fit read 8.5100e-4, \|Γ_ABC\|/\|Γ_PML\| 784.5; s″ = 0 \|Γ\| 1.000000; conjugated stretch \|Γ\| 9.1059e+2; NC-grade misses 9.9861e-4 (un-graded) / 4.6248e-3 (m = 3); PML k₀ read at λ₀/40 1.1821e-4 (the read's floor, printed only). W1 10 passed, 47 s (`20261022T203122Z_TH-5-step3-W1.log`); W2 real build 3 passed / 7 skipped (`20261022T203218Z_TH-5-step3-W2.log`). Normal incidence only. **Step 2c DONE 2026-10-08 (§9 item 320, attempt 1 of 4):** 2a's module landed on `main` with (R) retired (r(W) printed, never asserted); W1 7 passed, 225 s — (P0) −8.0404 pp (\|Δ\| 0.0080), (E) Re Z₁₁^ABC +3.784e-3 / +2.620e-3 / +1.924e-3 Ω, PEC \|Re Z₂₁\|/\|Z₂₁\| 0, (S) B^ABC +4.7133 / +2.3965 / +1.1256 pp, r(W) 0.4846 / 0.3567 / 0.2266 (record) — every digit 2a's (`20261008T183105Z_TH-5-step2c-W1.log:928–943`); W2 2b re-run 8 passed, 169 s, (M) B^mix −0.2520 pp against 1.9468 pp, (c2) 2.5647× (`20261008T183455Z_TH-5-step2c-W2.log:604–607`); W3 real build 2 passed, 9 skipped, 2 s (`20261008T183751Z_TH-5-step2c-W3.log`); 2a's W1 log now on `main` (`20261007T052330Z_TH-5-step2a-W1.log`). §2.1 line written (weekly 2026-10-07, verbatim). Prior state: steps 1a / 1b-0 / 1b landed 2026-09-30; step 2b (2026-10-07, §9 item 290) meets the mixed-wall bar at 0.08 m (\|B^mix\| 0.2520 pp against 1.9468 pp) and step 2a (item 277) is a banked negative — (S) green at every padding, (R) red at 0.12 m because r(W) 0.4846 / 0.3567 / 0.2266 falls with W as the natural wall's does. Ruled 2026-10-07 (weekly, §10 lane 2): 2a stands on its record, (R) retired as a sphere-model band on a cube; ✅ on step 2c (the module landed with (S) asserted, (R) printed) and the §2.1 line. *History: `docs/planning/chunks/TH-5.md`.* 

### §9 item 644 — 2026-10-25 (§9 rule (5) lane 2, `TH-5` step 3c: the graded PML's angular law at 45°, on the closed leg)

*Written by the 2026-10-25 18:00 daily review (fact sheet: the step 3b module, the Γ_c function and the angle table). Lane 2's anchor is a closed form (Jin 3e, cited below) at an imported band; no band, tolerance, record or Done-when is set or moved. Step 3b's own scope said "one angle … Not closed: other angles or a sweep in θ" (`:593` of this file).*

**Why 45° and why the closed leg.** At step 3b's absorbing design (R(0) = `R_ABSORB` = 1e-6), |Γ_c(θ)| = R(0)^{cos θ} sits *below* `PML_REFLECTION_BAR` (1e-4) at every angle under 60° (5.7e-5 at 45°, 6.4e-6 at 30°), so a run returning Γ = 0 would pass (the trap step 3 named: "at R_ABSORB, Γ_c is invisible under the discretisation term", `:465`). On step 3's **closed leg** (`R_CLOSED` = 1e-3, `S_CLOSED`, `test_th5_pml_reflection.py:59`, `:61`), |Γ_c(45°)| = (1e-3)^{0.7071} ≈ 7.6e-3 = 76× the bar, and its phase −2k₀ cos θ d (d = 1.5 λ₀) is not a multiple of π, so Γ_c is complex. **30° is excluded:** TE20 is exactly at cut-off on that guide (k_c2 = 2k₀ sin 30° = k₀). **75° is excluded:** (Cabc)'s ceiling there is 21, under `CABC_RATIO_BAR`, and λ_g = 3.9 m is 4× the 1 m fit window. 45° is clean: a = λ₀/√2 ≈ 0.707 m, β = k₀/√2 ≈ 4.443 rad/m, λ_g ≈ 1.414 m, TE20 evanescent (k_c2 = 2√2 π > k₀).

**The anchor.** `pml_reflection_pec_backed(S_CLOSED, PML_D, PML_M, K0, cos_theta=cos 45°)` (`src/fem_em_solver/core/boundary_terms.py:107–119`): Γ_c = −e^{−2jk₀ cos θ d}·exp(−2k₀ cos θ s_max d/(m + 1)), Jin 3e §9.6.1 eq. (9.133) integrated over the §9.6.4 profile (9.165) — the same citation as step 3b (the operator's print-copy check, dashboard item 0, still stands). Asserted:
- **(Gc45)** `| |Γ_c(45°)| / |Γ_c(0)|^{cos 45°} − 1 | ≤ IDENTITY_TOLERANCE` (imported from `test_port_reaction_impedance`, as step 3b does), and `| |Γ_c(45°)| − R_CLOSED^{cos 45°} | ≤ 1e-12·R_CLOSED^{cos 45°}`.
- **(B45)** degree 2, λ₀/40, PML at `S_CLOSED`, read at k_fit (step 3b's read, β scanned ±5e-3 relative, 21 points, `scan_interior` asserted): `|Γ_h − Γ_c(45°)| ≤ PML_REFLECTION_BAR`.

**The negative control (asserted, ceiling computed).** **(Law)** the normal-incidence law — the reflection a PML that ignored the angle would give on this leg, Γ_c(cos θ = 1) = −e^{−2jk₀d}·R_CLOSED = −1e-3 (d = 1.5 λ₀) — misses: `|Γ_h − Γ_c(cos θ = 1)| > 10 × PML_REFLECTION_BAR`. Ceiling: |Γ_c(45°) − Γ_c(0)| ≥ 7.6e-3 − 1e-3 = 6.6e-3 = 66× the bar, so the 10× bar has ≥ 6.6× margin. Also printed (not asserted): the 60° law on this leg, |Γ_h − Γ_c(cos 60°)| (|Γ_c(60°)| = 3.16e-2). **(Cabc) is printed only:** at 45° on the closed leg its ceiling is |Γ_ABC|/|Γ_c| = 0.1716/7.6e-3 ≈ 22.6, under `CABC_RATIO_BAR` by construction — the step-3b control does not transfer to this leg, and the brief does not pretend it does. **Degree 1 printed** (its λ₀/40 floor ≈ 1.2e-4 is at the bar, step 3b's reading).

**Route — a new module, no gate module edited.** `tests/validation/test_th5_pml_angle.py` (knob-free; `@complex_only`, `@pytest.mark.integration`). Write `_guide_at(comm, theta, n_per_lambda, degree, s_max)` modelled on step 3b's `_guide` (`test_th5_pml_oblique.py:76–98`), computing **locally** a = λ₀/(2 sin θ), β = k₀ cos θ, `ny = int(a / h)`, a local drive ẑ sin(πy/a), and the `_Projection` weight sin(πy/a). Import by name: from `test_th5_abc_reflection` — `FIT_O`, `K0`, `L_O`, `LAMBDA0`, `OMEGA`, `TAG_Y`, `TAG_Z`, `_box`, `_Projection`, `_solve`, `complex_only`; from `test_th5_pml_reflection` — `PML_D`, `PML_M`, `PML_REFLECTION_BAR`, `R_CLOSED`, `S_CLOSED`, `_solve_terms`; from `test_th5_pml_oblique` — `N_GATE`, `CABC_RATIO_BAR` (printed beside the ratio); `pml_reflection_pec_backed`, `pml_slab_term` from `fem_em_solver.core.boundary_terms`. **Never import `A_O`, `BETA_O`, `THETA_O` or `COS_O`** — they are 60°'s. The ABC guide for the printed (Cabc) uses `_solve(..., "abc")` on the 45° box. Γ_ABC(45°) closed form (Jin (9.48)): (cos θ − 1)/(cos θ + 1) = −0.1716, printed beside it.

**Windows.**
- **W0** (smoke, `-n 1`, complex): `--collect-only` plus the numpy node if you write one — the (Gc45) identity and the (Law) ceiling, no mesh.
- **W1** (heavy, `-n 2`, complex, `timeout -k 30 300`, after `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0): `scripts/testing/run_and_log.sh TH-5-step3c-W1 "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_th5_pml_angle.py -v -s --tb=short'"`. Then re-run step 3b's module unedited in the same window or a second one (`test_th5_pml_oblique.py`, 131 s), as the regression for the shared imports.
- **Prediction:** the 45° box has ≈ 1.2× step 3b's cells (ny 28 against 23; ≈ 40 000 at degree 2); degree-2 PML ≈ 80 s, degree 1 ≈ 14 s, ABC ≈ 19 s, mesh and scan ≈ 10 s ⇒ ≈ 125–160 s. Over 300 s is a kill: park, never raise the timeout.

**The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**

**Executor:** implementer. **Tier:** heavy, `-n 2`, complex. **Class:** lane-2 feature. **Cost:** ≈ 7.3 measured / 40 predicted slot-minutes. **Depends on:** nothing. **Family:** step 3c, attempt 1 of 4. **Daytime:** no.

**Traps already paid for.**
- Γ_c(45°) is complex: never port step 3b's `abs(gc - 1e-3)` literal (`test_th5_pml_oblique.py:132`); assert the magnitude identity and compare Γ_h to the complex Γ_c.
- Importing the module-level 60° constants silently re-builds the 60° guide; compute a, β, cos θ locally.
- Degree 1 at λ₀/40 sits at the bar (step 3b's execution reading) — print it, do not gate it.
- `ufl.max_value` and ordering comparisons on complex operands do not compile (`OPS-22`); a `SpatialCoordinate` facet integral needs a pinned `quadrature_degree`.
- `assemble_scalar` is rank-local; pytest captures prints without `-s`; durable capture with the trailing `; exit $rc` if you redirect; never pipe the harness command.

**Scope.** One more angle on the closed leg — the angular law R(0)^{cos θ} read at a second θ where the comparand is above the bar. **Not claimed:** the absorbing leg at 45° (Γ_c under the bar), a sweep in θ, grazing incidence, evanescent modes, the birdcage's box, any change to step 3b's gates.

**Negative result.** (B45) red with (Law) green: the PML misses its closed form at 45° — record both distances, park on `attempt/TH-5-step3c-*`, item 644 BLOCKED, known-issues entry if the miss reproduces. (Law) red: the solve cannot tell the angular law from the normal-incidence one at this fixture — report, stop. Never loosen `PML_REFLECTION_BAR`. No `PREDICTED RED:` line.

**Status it can move.** §7 `TH-5` stays ✅; its state line gains "step 3c: the angular law at 45° on the closed leg" as a registered record (green), or a known-issues entry (red). The dashboard's PML line ("works at 60°") becomes "at 45° and 60°".

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 644 — `TH-5` step 3c — per `docs/planning/chunks/TH-5.md` "§9 item 644". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Write a new module; never edit `test_th5_pml_oblique.py`, `test_th5_pml_reflection.py` or `test_th5_abc_reflection.py`, and never import their 60° constants. Spawned executors run every window foreground with Bash `timeout: 660000`."

### Review ruling 2026-10-26 03:00 (item 644)

Binds over the brief above. It rests on a `log-pathologist` ruling on `20261026T005627Z_TH-5-step3c-W1.log` (on `attempt/TH-5-step3c-20261026T0104Z`).
- **The window is UNCOUNTABLE** (Status 124, Elapsed 300, `W1:105–106`). (B45) 8.8931e-08 and (Law) 6.8743e-03 are prints from a killed run (`W1:93–95`). They are not passes and are not banked.
- **The journal's "rank 1 passed both nodes" is overruled.** The stream holds 26 `PASSED` tokens, 13 per rank. So neither rank finished `test_pml_angle_45_closed_leg`. The `PASSED` at the end of `W1:76` is rank 0's last migration verdict, landing on rank 1's node line. `_box` builds on `COMM_WORLD` (`test_th5_abc_reflection.py:107–134`), so rank 1 could not have finished the solve alone. The real cost is both ranks running the migration sweep. The rank-0 lag was a few seconds. The attempts entry (2026-10-26T01:03Z) and the known-issues 2026-10-26 Symptom row carry this correction.
- **The cost reading fits but is not shown.** The log has no per-phase times for `tests/environment`, and it cannot exclude a hang inside the (Cabc) solve. **The re-posed window is the discriminator.** It prints `[TH-5 s3c Cabc] … N s` with its own duration (module `:151`). Two outcomes:
  - It finishes in ≈ 140–175 s: the cost reading is confirmed and a hang is excluded.
  - It dies after (Law) again, well under 300 s: report it as a hang. Do not retry.
- **Re-posed as §9 item 649** (attempt 2 of 4) under the environment-first ruling (journal `2026-10-26-daily.md` §4). The parked module runs unedited.

### §9 item 649 — 2026-10-26 (§9 rule (5) lane 2, `TH-5` step 3c re-posed: the parked 45° module run unedited, `tests/environment` in its own window)

*Written by the 2026-10-26 03:00 daily review. This re-poses item 644 (attempt 2 of 4). The anchor (Jin 3e §9.6.1 (9.133) over (9.165)), the (Gc45) / (B45) gates, the (Law) control and its ceiling, the printed (Cabc) and degree 1, the scope, the negative-result disposition and the status it can move are item 644's above, unchanged. Nothing here sets or moves a band, tolerance, record, Done-when or target.*

**Land the parked work first (code and logs by path, record rows by hand).**
1. `git checkout attempt/TH-5-step3c-20261026T0104Z -- tests/validation/test_th5_pml_angle.py docs/testing/logs/20261026T005619Z_TH-5-step3c-W0.log docs/testing/logs/20261026T005627Z_TH-5-step3c-W1.log`.
2. Append the branch's two `test-results.md` rows **by hand** (`git diff main...attempt/TH-5-step3c-20261026T0104Z -- docs/testing/test-results.md`). **Never path-checkout `test-results.md`.** `git diff --stat` on it must show insertions only.
3. Do not edit the module.

**Pre-check (the pathologist's first step).** The 300 s kill may have landed inside an FFCx compile. When `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0, run the rubric's orphan-stem sweep as one bare harness call, labelled `TH-5-step3c-cache-sweep`. Report any orphan it deletes. A red in W1b that follows a deleted orphan is read against the stub signature before it is counted.

**Windows** (foreground, Bash `timeout: 660000`, `-s`; the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; `pgrep -c python3` reads 0 before each).
- **W1b** (heavy, `-n 2`, complex, `timeout -k 30 300`): `scripts/testing/run_and_log.sh TH-5-step3c-W1b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_th5_pml_angle.py -v -s --tb=short'"`. **Prediction ≈ 150 s (range ≈ 140–175 s; pathologist's estimate):** startup ≈ 3 s, degree-2 PML 103.4 s and degree 1 16.3 s (measured in the killed W1), (Cabc) ≈ 19–25 s (scaled from 3b's 15.9 s at ny 23 → 28). Every node must be PASSED on both ranks: count `PASSED` tokens per test ID, never by line adjacency.
- **Wreg** (heavy, `-n 2`, complex, `timeout -k 30 300`): step 3b's module unedited, as the regression for the shared imports: `scripts/testing/run_and_log.sh TH-5-step3c-Wreg "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_th5_pml_oblique.py -v -s --tb=short'"`. **Prediction ≈ 131 s** (`20261025T111709Z_TH-5-step3b-W1.log`).
- **W-env** (heavy by measurement, `-n 1`, `timeout -k 30 300`): `scripts/testing/run_and_log.sh TH-5-step3c-Wenv "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 300 python3 -m pytest tests/environment -v --tb=short --durations=5'"`. **Prediction ≈ 150–165 s.**

**Executor:** implementer. **Tier:** heavy, `-n 2` (W1b, Wreg) / `-n 1` (W-env), complex. **Class:** lane-2 feature. **Cost:** ≈ 7.3 measured / 25 predicted slot-minutes. **Depends on:** nothing. **Family:** step 3c, attempt 2 of 4. **Daytime:** no.

**Traps already paid for** (item 644's, plus): the rank-interleave trap above; a W1b that dies after (Law) **well under** 300 s is a hang, so report it and do not retry. A kill **at** 300 s is a price miss: park, never raise the timeout, and report the last phase line each rank printed.

**Status it can move:** as item 644. In the closing commit, delete `attempt/TH-5-step3c-20261026T0104Z` once its module, logs and rows are on `main`.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 649 — `TH-5` step 3c re-posed — per `docs/planning/chunks/TH-5.md` "§9 item 649". Read the brief, item 644's brief above it and the 2026-10-26 03:00 review ruling first; the ruling holds. Land the parked module and logs by path and the test-results rows by hand; never edit the module, `test_th5_pml_oblique.py`, `test_th5_pml_reflection.py` or `test_th5_abc_reflection.py`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-26, 09:46 UTC / 04:30 CDT slot, attempt 2 of 4) — DONE, green; `TH-5` stays ✅ and gains the step-3c record.** Run inline by the slot. Landed by path from `attempt/TH-5-step3c-20261026T0104Z`: the module and its two logs. The branch's two `test-results.md` rows were appended by hand (2 insertions, 0 deletions). The module is unedited.
- **Cache sweep** `20261026T094629Z_TH-5-step3c-cache-sweep.log` (164 s, Status 0) ran after `pgrep -c python3` read 0. It found **no orphan**, so the 300 s kill had not left a half-written compile.
- **W1b** `20261026T094918Z_TH-5-step3c-W1b.log` (`-n 2`, complex, **109 s** against a predicted ≈ 150, Status 0). Both rank summaries read "6 passed", with 12 `PASSED` tokens, two per test ID.
  - (Gc45) identity 6.661e-16.
  - Γ_c(45°) = −5.469868e-3 + 5.222615e-3j.
  - **(B45)** degree 2, 40 320 cells (ny 28), 73.9 s: \|Γ_h − Γ_c(45°)\| = **8.8931e-08** ≤ `PML_REFLECTION_BAR` 1e-4.
  - **(Law)** 6.8743e-03, 68.7× the bar (asserted > 10×).
  - Printed only:
    - degree 1: 3.8115e-05 (12.5 s);
    - the 60° law: 3.7459e-02;
    - **(Cabc)**: 19.5 s, Γ_ABC −1.715730e-01 against Jin (9.48) −1.715729e-01, ratio 22.69, under `CABC_RATIO_BAR` by construction, as the brief said.
  - **The ruling's discriminator:** (Cabc) printed its own duration, and the window finished well inside the ceiling. That excludes a hang. Item 644's kill was the migration sweep's cost: (B45) and (Law) match the killed run's prints to the digit.
- **Wreg** `20261026T095117Z_TH-5-step3c-Wreg.log` (`-n 2`, 121 s, Status 0). Step 3b's module ran unedited: "5 passed" per rank, 10 tokens. (B) 1.1031e-07 / degree 1 6.6733e-05 equal `20261025T111709Z_TH-5-step3b-W1.log:61` to the digit.
- **W-env** `20261026T095326Z_TH-5-step3c-Wenv.log` (`-n 1`, 105 s, Status 0): 12 passed.
- The attempt branch is deleted in this slot. The dashboard's PML line ("at 45° and 60°") is the daily review's to update.
