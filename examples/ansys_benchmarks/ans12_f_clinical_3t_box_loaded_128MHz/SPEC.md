# ANS-12 — The ASTM-F2182-box-loaded F-clinical-3T birdcage under the two-port quadrature drive at its in-model Larmor root: `|B₁⁺|`, the located 10 g hotspot and the feed quantities at 128 MHz

**Commissioned 2026-11-09 by the weekly planning review (PROJECT_PLAN §10
chain C — the 3 T `ANS` case posed by `xl` entry 58's decision rule,
`docs/testing/xl-pending.md` entry 58).** The 3 T twin of `ANS-11`: the
production fixture's 3 T variant (operator directive 2026-10-14, "then the 3 T
twin") loaded with the standard's gel box, no implant. Our degree-1 column
exists (`xl` entry 58, `20261109T080008Z_GEO-38-C6-3T.log`, a record); **the
adjudication column is the 3 T degree-2 set, `xl-pending.md` entry 66**, not
yet run. Every value below is taken from the code at the cited line or from
`ANS-11`'s SPEC where the two fixtures share it.

## Open before the operator starts

1. **The node — TBD (SPEC v2).** The adjudication column is at degree 2 (the
   10-12 ruling; `ANS-8` v4's and `ANS-11` v2's precedent), so the capacitor
   is **C_in(F-clinical-3T, d2) at 0.25 pF** — `xl` entry 64 (Sunday
   2026-11-15, the three-node degree-2 grid) and its refinement (entry 56's
   two-node straddle, the first open night after the 11-16 weekly reads 64).
   The weekly that reads the refinement writes SPEC v2 with C, R_p ‖ C and
   `Z_th`. **The 44.0525 pF of entry 58's record column is not the node to
   enter.** Until v2: build the geometry, materials, sheets and ports;
   do not solve.
2. **The `|B₁⁺|` mean — ruled (as `ANS-11`): the comparison row is the
   18-point lattice mean**; the 400-centroid mean is a second record (row 4).
3. **Readout bands — `ANS-8` step 2's** (its SPEC `:315–329`), with `ANS-11`
   v2.1's row-3 tie clause; adjudicated on entry 66 only.
4. **The gel at 128 MHz — provisional.** σ = 0.5 S/m, ε_r = 78 are the
   2026-10-21 weekly's provisional ruling (dashboard W3), not the F2182
   clause's nominal; our column carries them
   (`test_geo38_f_clinical_3t_esr_root.py:426`). If the operator's copy of
   the standard gives 128 MHz values, they go into the repo under
   `docs/references/` first and both columns are re-posed on them.

## Status

- [x] SPEC v1 (2026-11-09, skeleton; the node TBD).
- [ ] SPEC v2 (the node — after `xl` entry 64 and its refinement).
- [ ] Runnable half (`ans:12`): imports C21's module
      (`tests/validation/test_geo38_f_clinical_3t_box_loaded.py`), never
      re-implements it; `COMPARISON.md` with our rows, AED columns blank
      (§10 C48 / C49).
- [x] Our degree-1 column (`xl` entry 58, 2026-11-09) — **a record** (below).
- [ ] Our degree-2 column (entry 66, `xxl`) — **the adjudication column**.
- [ ] Operator replication in AED (HFSS + Circuit), First Order.
- [ ] Adjudication (weekly, on entry 66).

## Geometry (SI units, metres)

The fixture is `MeshGenerator.birdcage_port_domain(**F_CLINICAL_3T_KWARGS)`
(`tests/validation/test_geo38_f_clinical_3t.py:44–87`: the 1.5 T table with
the 3 T overrides) with the phantom cylinder replaced by the box
(`test_geo38_f_clinical_box.py:54–65`, the same box as `ANS-11`),
`conductor_resolution = 1.6e-3` (rung r1). The nominal 3 T table is
`RUNG_DIAMETER_3T_M = 0.704`, `SHIELD_DIAMETER_3T_M = 0.743`,
`SHIELD_LENGTH_3T_M = 0.70`, `RUNG_LENGTH_3T_M = 0.58`
(`test_geo38_f_clinical_3t.py:45–48`); everything not listed is `ANS-11`'s.

| Part | Definition |
|---|---|
| RF shield / outer boundary | Closed **PEC can**, radius **R_s = 0.3715** (Ø 0.743), **z ∈ ±0.35** (length 0.70: `air_padding` is computed so the can spans ±0.35, `:78–86`). Lateral face and both end caps PEC. |
| Legs (rungs) | **32** solid cylinders, radius **0.006** (leg_width 0.012), axis ∥ ẑ at **R = 0.352** (rung Ø 0.704), azimuths **θ_n = n·11.25°**, n = 0…31, spanning **z ∈ ±0.29** (`coil_length` = 0.58). Unbroken: high-pass. |
| End rings | Two tori, major radius 0.352, minor radius **0.004**, centred at **z = ±0.2578571429** (`leg_spacing` = 0.58 − 2 × 0.0321428571 = 0.5157142858, the 1.5 T leg overhang kept, `:73–81`). Conductor = union of legs and rings. |
| Ring gaps | On **both** rings, at the 32 mid-azimuths **φ_j = (j + ½)·11.25°**, remove the arc `\|φ − φ_j\| ≤ α`, **α = g/(2R) = 0.0171428571 / 0.704 = 0.0243506493 rad**, radial cut planes: **64 gaps**, each leaving two planar circular terminal faces of radius 0.004. |
| Gap sheets (exact rectangles, κ off) | In each gap one planar rectangle in the plane **u = R**, centred at (R cos φ_j, R sin φ_j, ±0.2578571429): **chord 2R·tan α ≈ 0.0171462** along φ̂ (review arithmetic; the module prints the exact area) × **0.010** along ẑ. Current line along +φ̂ at radius R. |
| Gel box (ASTM F2182) | **0.42 (x) × 0.09 (y) × 0.65 (z)**, centred on the origin — `ANS-11`'s box unchanged (`test_geo38_f_clinical_box.py:60–62`; `TABLE_OFFSET_M = 0.0`). |

**Port numbering (ours):** as `ANS-11`'s — gap j is **P(33 + j) on the bottom
ring (z = −0.2578571429)** and **P(65 + j) on the top ring**. **Feeds: P33
(φ = 5.625°) and P41 (φ = 95.625°)**, both bottom ring, 8 slots = 90° apart
(`test_geo38_f_clinical_3t_esr_root.py:289`).

*In our mesh but not modelled in AED, by design:* (i) 32 floating leg-port
boxes 4 × 4 × 21.43 mm at the mid-leg azimuths (inert vacuum; the 3 T leg-box
corner sits at r = 0.363814913104, `test_geo38_f_clinical_3t.py:90`); (ii) a
rod-shaped sub-volume (tag 403, radius 1.0e-3, length 0.030 on ẑ, centred at
(0.105, 0, 0)) **carrying the gel's material** (entry 58's log `:12011`:
14 848 cells, gel-mapped) — a refinement artefact of C20's rod-in-box mesh
(1 128 430 cells, U₁ 1 508 256), not an implant. The AED gel is one
homogeneous box.

## Materials

| Region | σ [S/m] | εᵣ | μᵣ | ρ [kg/m³] |
|---|---|---|---|---|
| Air inside the can | 0 | 1 | 1 | — |
| Coil (legs ∪ rings) | **copper 5.8e7**: ours Leontovich `Z_s = (1 + j)R_s` on the conductor surface, R_s = 2.952e-3 Ω at 128 MHz (√2 × `ANS-11`'s 2.087e-3 at 64 MHz; `src/fem_em_solver/core/cavity.py:543–546`); **AED: Finite Conductivity σ = 5.8e7, Solve Inside off** (`ANS-6`'s treatment) | 1 | 1 | — |
| Gel box | **0.5** (provisional, item 4 above) | **78** (provisional) | 1 | **1000** |

## Excitation — all 64 gaps capacitors with ESR, two feeds in quadrature

- **Capacitor on every one of the 64 gaps: C = TBD (SPEC v2)** with a
  **0.2 Ω series ESR** (binding ruling 2026-10-06). In HFSS as the parallel
  equivalent at 128 MHz, high-Q form as `ANS-8` v4 / `ANS-11` v2: R_p =
  1/(ω²C²r) ‖ C, written by v2. *(The degree-1 record column ran at
  44.0525 pF: `Z_RC = 0.2 − 28.225j Ω` at 128 MHz, R_p ≈ 3 983 Ω ‖
  44.0525 pF, Q_C ≈ 141 — review arithmetic, for the record column only;
  not to enter.)*
- **62 non-feed gaps:** Lumped RLC sheets (R_p ‖ C).
- **Two feed gaps P33, P41:** 50 Ω lumped ports **with the same R_p ‖ C
  present across each feed gap in the solved field** — ours drives each feed
  through `Z_th = 50 ‖ Z_RC` with V_th = 1 V (the module's form, as
  `ANS-11`'s `:144–151`); v2 writes the number. State how the feed-gap
  R_p ‖ C is realised.
- **Drive (ii), primary: ccw = (V_33, V_41) = (1, −1j)**; **cw = (1, +1j)**
  printed beside it at c\* only — the conventions are `ANS-11`'s.
- **Normaliser: P_acc** — the sum over the two feeds of HFSS's accepted
  power at the lumped port (the 2026-10-12 ruling; `ANS-8` SPEC `:291–304`);
  ours is the lumped line cut (`test_geo38_f_clinical_box_loaded.py:376–381`,
  run unchanged under C21's patch set). P_box and P_src printed beside it.
- **No re-tune on the AED side**; state so in the results file. Im Y_in at
  the node is a record, not a root check: the root was found on C12's
  cylinder stand-in and the box detunes the coil (entry 58: the two feeds'
  Im Y_in have opposite signs under both senses, `:12232–12233`).

## Frequency, order and mesh

- **128 MHz only** (`FREQUENCY_3T_HZ`, `test_geo38_f_clinical_3t.py:49`; entry
  58's log `:138`). No sweep.
- **First Order** (operator order ruling 2026-10-01) — our counterpart is the
  degree-2 column (entry 66). A Zero Order column (our degree 1, entry 58) is
  optional and a record.
- Adaptive at 128 MHz to **ΔS ≤ 0.002** on the 2×2, ≥ 2 consecutive converged
  passes, max 10, ≥ 2 elements across each gap along φ̂; report if the cap is
  hit. Mixed Order forbidden. Our mesh: 1 128 430 tets (C20's rod-in-box mesh
  with the rod gel-mapped), conductor 1.6 mm, global 15 mm. Resources: 20
  cores / 1 TB (`ANS-8` SPEC `:158`).

## Quantities to export (all printed digits; per √W or per W of P_acc where stated)

`ANS-11`'s seven rows unchanged: (1) `S` 2×2 (P33, P41) at 50 Ω, 128 MHz, and
`Y_in` at each feed under drive (ii); (2) `|B₁⁺|` at the 18 lattice points —
x ∈ {−0.105, 0, +0.105}, y = 0, z ∈ {−0.25, −0.15, −0.05, +0.05, +0.15,
+0.25} (the same lattice; C21 runs C6a's `_lattice` unchanged) — per √W of
P_acc; (3) SAR10g at the 18 points (ball radius 0.013365; a C95.3 cube is
admissible — state which), c\* = the largest with its (r, dz, az), per W of
P_acc; the box-wide 10 g maximum and location (record); (4) the gel-volume
mean `|B₁⁺|`, per √W of P_acc; (5) power partition per P_acc: gel / copper /
capacitor-ESR fractions; (6) mean `|B₁⁺|` / mean `|B₁⁻|` over the gel; (7)
solve metadata (tets, passes, final ΔS, converged yes / no, order, unknowns,
wall time; whether AED was re-tuned — must be no).

## Readout (pre-registered, adjudicated on entry 66 only)

`ANS-11`'s rows and bands with the v2.1 tie clause: row 1 Im Y_in of the same
sign at both feeds ⇒ record, a sign flip is a detuning finding (our degree-1
column reads opposite signs at the degree-1 root — the degree-2 column at the
degree-2 root is the one adjudicated); row 2 within **10 %** per point with
the z-ordering at x = 0 preserved ⇒ AGREE; row 3 c\* within **15 %** and
within one ball radius of any lattice ball in the tie set (SAR10g within 1 %
of our maximum) ⇒ AGREE, the box-wide maximum a record; row 5 gel fraction
within **5 %**, copper and ESR within **10 %** ⇒ AGREE; row 6 a ratio ≫ 1 in
both codes (ours 1.67 at degree 1) is the sense check. No absolute tuning,
compliance, homogeneity or C95.3 claim follows.

## Our column (degree 1, `xl` entry 58 — a record; raw values at V_th = 1 V per feed)

`20261109T080008Z_GEO-38-C6-3T.log`, `-n 8`, 780 s, `memory.peak` 33.48 GiB;
every gate green ((R1) 1.893e-13, (P) 5.120e-13 / 1.810e-12 with NEG
1.2375e-02 / 1.2827e-02, (P_RC) 5.120e-13, coverage (i) 2.975e-14 / 3.220e-15;
build 448.3 s, columns 194.1 + 13.4 s, readouts 117.1 s, `:12261`).

| Quantity | Value | Log line |
|---|---|---|
| S₁₁, S₁₂ = S₂₁, S₂₂ | 3.6878613057e-01 − 6.7206972791e-02j; −1.5727100565e-02 + 4.2880439168e-02j; 2.0011568377e-01 + 6.0905743336e-02j | `:12228` |
| Y_in ccw, P33 / P41 [S] | +8.297670657e-03 + 1.031943733e-03j / +1.451253979e-02 − 1.347374887e-03j (Im of opposite signs) | `:12232` |
| Y_in cw, P33 / P41 [S] | +1.005048930e-02 + 1.879628016e-03j / +1.205845746e-02 − 1.976458913e-03j | `:12233` |
| 18-point lattice mean `\|B₁⁺\|`_ccw (the comparison row) | 1.174932617e-08 T | `:12235` |
| 400-centroid gel mean `\|B₁⁺\|`_ccw / mean `\|B₁⁻\|` | 1.562098621e-08 T / 9.328405808e-09 T; (Sense) 1.674561; CV 0.359276 | `:12237` |
| c\* | x +0.105, z +0.050 (r 0.105, az 0°): SAR10g 1.748477085e-04 W/kg (ccw); cw ball 1.427632775e-04 | `:12256` |
| P_box / P_src / P_acc [W] | 3.659141247e-03 / 5.347172304e-02 / 1.874513035e-02 | `:12259` |
| SAR10g(c\*) per P_acc; centroid mean `\|B₁⁺\|`/√P_acc | 9.327634e-03 W/kg/W; 1.140944e-07 T/√W (lattice mean: 8.581609e-08 T/√W, C33) | `:12259` |
| 3 T / 1.5 T per accepted watt (C33 reader, record) | lattice `\|B₁⁺\|`/√P_acc 1.410682; centroid 1.344035; SAR10g(c\*) 4.101440; P_box / P_acc 3.693500 | `20261109T103557Z_GEO-38-C33-W1-r2.log:195–198` |

Lattice `|B₁⁺|`_ccw [T] / SAR10g [W/kg] (`:12238–12255`; the `ans:12` page
(C48) carries all digits and the per-P_acc columns):

| x | dz −0.25 | −0.15 | −0.05 | +0.05 | +0.15 | +0.25 |
|---|---|---|---|---|---|---|
| −0.105 | 9.020981e-09 / 1.014257e-04 | 1.474242e-08 / 1.247755e-04 | 1.712664e-08 / 1.714032e-04 | 1.747062e-08 / 1.674858e-04 | 1.460517e-08 / 1.249719e-04 | 8.644831e-09 / 9.714595e-05 |
| 0 | 7.447523e-09 / 1.313789e-04 | 7.651445e-09 / 3.791187e-05 | 7.420153e-09 / 6.935344e-06 | 7.351443e-09 / 7.045244e-06 | 7.356692e-09 / 3.746219e-05 | 7.355732e-09 / 1.257708e-04 |
| +0.105 | 8.849241e-09 / 1.009853e-04 | 1.534804e-08 / 1.270816e-04 | 1.819490e-08 / 1.734193e-04 | 1.821289e-08 / 1.748477e-04 | 1.552658e-08 / 1.263618e-04 | 9.162570e-09 / 9.877796e-05 |

Row 3's tie set on this column: x+105 mm dz+50 mm (1.748477e-04), x+105 mm
dz−50 mm (1.734193e-04, −0.82 %), x−105 mm dz−50 mm (1.714032e-04, −1.97 %,
outside 1 %). The 1.5 T column's near-tie was between mirror balls 210 mm
apart (`ANS-11` v2.1); here the tie is between the two z-mirror balls on the
+x side.

## Our side

The adjudicated column is **entry 66** — C21's module at degree 2 on the same
mesh (`xxl`; node per item 1) — read through the same rows; the degree-1 table
above is a labelled record, as `ANS-11`'s degree-1 rows are. `COMPARISON.md`
carries our rows with the AED columns blank (C49).

**AED / Ansys values are private (CLAUDE.md hard rule): they go in the
gitignored `aed_results/` and `COMPARISON_private.md` only; the qualitative
verdict and a relative agreement level are all that reaches a tracked file.**

**Out of scope:** no implant (`GEO-36` step 11 / entry 61 is separate); no
thermal half (chain P, entry 60); no matching network; no 1.5 T row
(`ANS-11`); no homogeneity or compliance statement; no claim at degree 1.

## Revision history

- **v1.1 (2026-11-11, weekly planning review; citations only):** four
  log-line citations into entry 58's log were off by one (S `:12228`, Y_in
  ccw `:12232`, cw `:12233`, price `:12261`; found by the 2026-11-10 03:00
  daily review); no value, row or band changes.
- **v1 (2026-11-09, weekly planning review; skeleton):** fixture, materials,
  excitation, rows and the entry-58 record column as above; items 2 and 3
  ruled as `ANS-11`'s; item 1 (the node) is v2's, from `xl` entry 64 and its
  refinement; item 4 (the gel at 128 MHz) stays provisional (W3).
