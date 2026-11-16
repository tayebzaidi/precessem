# ANS-9 — Coil-driven heating of the gel phantom at 10 MHz: HFSS SAR as the heat source of a transient thermal solve (one-way, gel, no perfusion)

**Commissioned 2026-10-03 (weekly planning review; PROJECT_PLAN §10 chain P4,
operator directive 2026-09-30).** Phase 8's first external check: every
thermal gate so far is a closed form on a sphere or cylinder, or an identity
on our own coil-driven solve (`MAT-15` ✅ 2026-10-02). The trigger written
into the placement — `MAT-15` step 1 landing — fired on 2026-10-01.

**Why this case, stated as the gap it closes.** `MAT-15` asserts that our
thermal solve conserves energy, starts with the calorimetric slope
`dT/dt|₀ = SAR/c` and respects C4 — identities a solve wrong by a constant
factor in the EM → thermal hand-off (a time-average factor of two, a
per-volume against per-mass source, a power normalisation) passes exactly.
This case is the first independent check of the hand-off, and it is cheap
for the operator because **the EM half is byte-for-byte `ANS-2`'s project**
(adjudicated AGREE 2026-09-18) plus one thermal design.

This is step 1 of `ANS-9`'s progression (one-way, gel). Steps 2 (two-way
σ(T)) and 3 (perfusion) are not specified here.

## Status

- [x] Runnable half implemented (2026-10-04, §9 item 218) (§7 chunk `ANS-9` step 1a, to be scoped by
      the daily review: a script in this directory dispatching through
      `./run_examples.sh -e ans:9`, importing every constant and band from
      `tests/validation/test_mat15_coil_driven_heating_phantom.py` —
      nothing restated — writing `metrics.json`, `COMPARISON.md` with our
      columns filled and the AED columns blank, and combined XDMF)
- [x] Operator replication in AED (HFSS + the Ansys transient thermal
      solver the operator's licence carries — Icepak or Mechanical thermal;
      **state which in the results file**) — **done 2026-10-06** (results
      in the gitignored `aed_results/`, recorded here 2026-10-08 by an
      interactive session because the 10-07 weekly still listed the run as
      owed). Thermal solver: AEDT 2026 R1 IcepakFEA (PyAEDT `Mechanical`,
      Transient Thermal), First Order EM design. **Read before
      adjudicating:** (i) the raw harvest (the private results JSON in the
      gitignored `aed_results/`) is superseded for rows 2 / 3 / 6 by the
      operator's offline correction pass (its corrected-rows JSON, beside it
      in `aed_results/`; script: the operator's offline processing script in
      the gitignored `aed/`; it rescales the field-calculator
      volume / surface integrals from the thermal mesh's own V / A to the
      exact cylinder's — the raw row 3 is not usable without it);
      (ii) AED's time stepping departs from the SPEC: one setup per report
      time, each a single load step with AEDT's automatic substepping, not
      40 × 31.35 s steps — state it in the adjudication; (iii) row 4 is a
      point value (monitor points). An interactive first pass read rows 2,
      3, 6 and A–C inside their bands, B / C / D ordering preserved, and
      probe D outside its 10 % band — the weekly's adjudication decides.
- [x] Adjudication (weekly; numbers in the gitignored
      `docs/private/ans9-adjudication-<date>.md`) — **done 2026-10-10**
      (`docs/private/t13-t29-ans9-adjudication-2026-10-10.md`): AGREE on
      rows 1 / 2 / 3 / 6 and probes A–C, DISAGREE at probe D (ours 11 % below
      AED, band 10 %); verdict table in `COMPARISON.md`; known-issues
      2026-10-10; step 1b (PROJECT_PLAN §10 P21) diagnoses D.
      *Superseded 2026-10-14 (next item).*
- [x] Re-adjudication on the degree-2 column (weekly, **2026-10-14**;
      numbers in the gitignored
      `docs/private/ans9-readjudication-2026-10-14.md`; no new AED run) —
      on step 1d's primary column (N1curl degree 2, exact rectangles, κ off,
      h_p 0.0025 m, the production convention) every row **AGREES**: probe
      D within 2.5 % and probes A–C within 2.3 % (band 10 %), B > D > C
      held in both codes; ⟨θ⟩ within 0.8 % at every report time, the
      initial slope within 0.5 %, the surface heat within 1.3 % on the
      Simpson figure (bands 5 %); rows 5 / 7 records. The 2026-10-10
      DISAGREE at D was our discretisation — 13 % phantom resolution and
      8 % element order (steps 1b / 1c). `ANS-9` ✅ (step 1); known-issues
      2026-10-10 retired; verdict table in `COMPARISON.md`.

## EM half — `ANS-2`, unchanged

Geometry, materials, boundary condition, the four lumped ports (50 Ω, exact
0.007 × 0.008 rectangles), solution frequency **10 MHz**, mesh and
convergence settings: **exactly `ANS-2`'s SPEC.md** (itself `ANS-4`'s
geometry with the phantom's mass density). If this file and `ANS-2`'s ever
disagree on the EM half, `ANS-2`'s is authoritative. Use the **First Order**
design (operator order ruling 2026-10-01); a Zero Order column is optional.

**Excitation for the thermal source:** port **P1 driven alone**, the other
three terminated in 50 Ω, source scaled so that the **power absorbed in the
phantom is exactly 1 W** (`∫_phantom ½ σ |E|² dV = 1 W`; in HFSS: edit
sources so that the volume-loss integral over the phantom reads 1 W, and
record the incident power that gives it). This is our `G = 1 W / P₁`
normalisation; it removes the port-feed systematic from the thermal
comparison, so a thermal disagreement cannot be blamed on the feed.

## Thermal half

| Item | Value |
|---|---|
| Thermal domain | **the phantom only** (cylinder radius 0.030 m, `z ∈ [−0.040, +0.040]` m) — no air, no coil |
| Thermal conductivity k | **0.6 W/(m·K)**, uniform, temperature-independent |
| Volumetric heat capacity ρc | **4.18e6 J/(m³·K)** (ρ = 1000 kg/m³, c = 4180 J/(kg·K)) |
| Perfusion / flow | none (gel) |
| Heat source | the EM loss density `½ σ |E|²` of the P1 drive above, mapped onto the thermal mesh (one-way; σ does not depend on T) |
| Boundary, all three phantom faces | convection (Robin): **h = 10 W/(m²·K)** to ambient **T∞ = 20 °C**; no radiation |
| Initial condition | uniform **20 °C** |
| Duration | **1254 s** (= 0.2·a²ρc/k), source on throughout |
| Time step | ≤ 31.35 s (40 steps or finer); state the scheme and step used |

If the thermal solver cannot take a convection coefficient on a bare solid
face without a fluid region, state what was used instead; do not substitute
a fixed-temperature wall silently.

## Quantities to export (rise θ = T − 20 °C, in K)

| Row | Quantity |
|---|---|
| 1 | Total power deposited in the phantom by the mapped source, W (should read 1.000; the mapping's conservation) |
| 2 | Volume-mean rise ⟨θ⟩ at t = 313.5, 627, 940.5, 1254 s |
| 3 | Initial slope of the volume-mean rise, K/s, from the first time step (calorimetry: equals 1 W / (ρc·V)) |
| 4 | Point rise at t = 1254 s at five probes, metres: A (0, 0, 0); B (0.020, 0, 0); C (0, 0.020, 0); D (−0.020, 0, 0); E (0, 0, 0.030) |
| 5 | Maximum rise in the phantom at t = 1254 s and its location |
| 6 | Total heat leaving through the Robin faces integrated over 0–1254 s, J (energy balance: row 1 × 1254 s = ρc·V·⟨θ⟩(1254) + row 6) |
| 7 | Incident power at P1 that gives 1 W in the phantom, W |

**Probe volume for row 4 (added 2026-10-05 by the weekly review, after
`MAT-19` step 2).** Row 4 may be exported either as the field value at the
probe point or as the mean over a ball of radius ≤ 5 mm centred on it — state
which. On our side the 5 mm ball mean sits within 1.9 % of the point value at
A–E (−1.24 % at B, the largest) and moves 0.44 % when B is displaced by
0.1 mm, both inside row 4's 10 % readout, so neither choice changes the
adjudication. No new row and no extra run is asked for.

Also export the solver, version, element order of both solves, thermal mesh
size, and wall time.

## Readout (pre-registered)

- Rows 2, 3 and 6 are integral quantities: agreement within **5 %** (the
  level `ANS-2`'s phantom power integral was adjudicated at) is AGREE.
- Row 4 is pointwise on a source with a hot region near the driven leg:
  within **10 %** at probes A–D, with B / C / D ordered the same way in both
  codes, is AGREE; probe E is a record.
- Row 5 is a record (a maximum on a mapped field is mesh-dependent).
- Row 7 against ours measures the feed, not the thermal solve; it is
  `ANS-2`'s question and is recorded, not adjudicated here.
- A disagreement on rows 2 / 3 with row 1 at 1 W on both sides is a finding
  about one code's thermal solve or source mapping; a known-issues entry and
  a diagnosis chunk follow.

AED / Ansys values are private (CLAUDE.md hard rule): they go in the
gitignored `aed_results/` and `COMPARISON_private.md` only.
