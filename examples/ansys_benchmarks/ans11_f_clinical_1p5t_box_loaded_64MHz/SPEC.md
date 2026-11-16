# ANS-11 — The ASTM-F2182-box-loaded F-clinical-1.5T birdcage under the two-port quadrature drive at its in-model Larmor root: `|B₁⁺|`, the located 10 g hotspot and the feed quantities at 64 MHz

**Commissioned 2026-11-04 by the weekly planning review (PROJECT_PLAN §10
chain C, step C7 — the `ANS` case posed on F-clinical by `xl` entry 52's
decision rule, `docs/testing/xl-pending.md` entry 52).** The first external
check at clinical scale: the production fixture (operator directive
2026-10-14) loaded with the standard's gel box, no implant. Our degree-1
column exists (`xl` entry 52, `20261103T080009Z_GEO-38-C6.log`, a record);
**the adjudication column is the degree-2 one, C19** (`xl-pending.md:2335`),
not yet run. Every value below is taken from the code at the cited line.

## Open before the operator starts — all three ruled (2026-11-04 and 2026-11-07 weeklies); nothing is left to decide

1. **The node — ruled 2026-11-07 (weekly, SPEC v2): C = 185.7629 pF.**
   The adjudication column is C19 at degree 2 (the 10-12 ruling; `ANS-8` v4's
   precedent), so the capacitor is **C_in(F-clinical-1.5T, d2) at 0.25 pF
   from `xl` entry 56 (2026-11-05)**: one − → + crossing of Im Y_in on P33
   (1, +1j), 185.762904 pF by the bracketing pair
   (`20261105T080008Z_GEO-38-C5b-prime.log:12814`), −4.05 % from the degree-1
   root; written to four decimals as the degree-1 root was (193.5996 from
   193.599624). The 193.5996 pF of entry 52's record column **is not the node
   to enter**; the R_p ‖ C figures in "Excitation" are recomputed below.
   **Nothing else in this SPEC changes with the node.** Our column at this
   node is `xxl` entry 62 (Saturday 2026-11-14, READY).
2. **The `|B₁⁺|` mean — ruled: the comparison row is the 18-point lattice
   mean** (the readout chain C registered); the 400-centroid mean is a second
   record, kept beside it in row 4 with its definition. Both are quoted in
   `COMPARISON.md`.
3. **Readout bands — ratified** as `ANS-8` step 2's (its SPEC `:315–329`),
   adjudicated on C19 only; the degree-1 column is never adjudicated.

## Status

- [x] SPEC v1 (2026-11-04, skeleton; TBDs 2 and 3 ruled that night).
- [x] SPEC v2 (2026-11-07, the node entered — ready to build and replicate).
- [ ] Runnable half (`ans:11`): imports C6's module
      (`tests/validation/test_geo38_f_clinical_box_loaded.py`), never
      re-implements it; `COMPARISON.md` with our rows, AED columns blank.
- [x] Our degree-1 column (`xl` entry 52, 2026-11-03) — **a record** (below).
- [ ] Our degree-2 column (C19, `xxl`) — **the adjudication column**.
- [ ] Operator replication in AED (HFSS + Circuit), First Order.
- [ ] Adjudication (weekly, on C19).

## Geometry (SI units, metres)

The fixture is `MeshGenerator.birdcage_port_domain(**F_CLINICAL_15T_KWARGS)`
(`tests/validation/test_geo38_f_clinical_build.py:44–64`) with the phantom
cylinder replaced by the box (`tests/validation/test_geo38_f_clinical_box.py:54–65`),
`conductor_resolution = 1.6e-3` (rung r1, `:71`). Values evaluated from
those lines; closed forms from `src/fem_em_solver/io/mesh.py`.

| Part | Definition |
|---|---|
| RF shield / outer boundary | Closed **PEC can**, radius **R_s = 0.3725** (`test_geo38_f_clinical_build.py:59`; Ø 0.745), `z ∈ ±0.4771428571` (L_z/2 + air_padding 0.1521428571, `test_geo38_f_clinical_build.py:58`; `mesh.py:5700`, `:4587`); length 0.954285714 (`GEO-38.md:450`). Lateral face and both end caps PEC (`mesh.py:4316–4320`). |
| Legs (rungs) | **32** solid cylinders, radius **0.006** (leg_width 0.012, `test_geo38_f_clinical_build.py:49`; `mesh.py:4770`, `:5431`), axis ∥ ẑ at **R = 0.355** (`test_geo38_f_clinical_build.py:46`; rung Ø 0.710), azimuths **θ_n = n·11.25°**, n = 0…31, spanning **z ∈ [−0.2571428571, +0.2571428571]** (coil_length = 0.45 + 2·0.0321428571, `test_geo38_f_clinical_build.py:48`). Unbroken: high-pass. |
| End rings | Two tori, major radius 0.355, minor radius **0.004** (`test_geo38_f_clinical_build.py:50`), centred at **z = ±0.225** (leg_spacing 0.45, `test_geo38_f_clinical_build.py:47`). Conductor = union of legs and rings. |
| Ring gaps | On **both** rings, at the 32 mid-azimuths **φ_j = (j + ½)·11.25°**, j = 0…31, remove the arc `\|φ − φ_j\| ≤ α`, **α = g/(2R) = 0.0241448692 rad** (g = 0.0171428571, `test_geo38_f_clinical_build.py:51`; `mesh.py:4442–4448`), radial cut planes: **64 gaps**, each leaving two planar circular terminal faces of radius 0.004. |
| Gap sheets (exact rectangles, κ off) | In each gap one planar rectangle in the plane **u = R** (normal r̂(φ_j)), centred at (R cos φ_j, R sin φ_j, ±0.225): **chord 2R·tan α = 0.0171461892** along φ̂ × **0.010** along ẑ (2·0.004 + 2·0.001 port clearance, `test_geo38_f_clinical_matrix.py:63`; `mesh.py:4461–4470`); area 1.714618916e-04 (log `:143`). Current line along +φ̂ at radius R. |
| Gel box (ASTM F2182) | Axis-aligned box **0.42 (x) × 0.09 (y) × 0.65 (z)**, centred on the origin: `x ∈ ±0.21, y ∈ ±0.045, z ∈ ±0.325` (A, B, C = 0.325, 0.21, 0.045, `tests/validation/test_mat21_box_closed_form.py:34`; size (2B, 2C, 2A), `test_geo38_f_clinical_box.py:60–62`; `TABLE_OFFSET_M = 0.0`, `test_geo38_f_clinical_box.py:51`; `mesh.py:4300–4303`). The long side is the coil axis; the 0.09 depth is vertical (y). |

**Port numbering (ours):** gap j is **P(33 + j) on the bottom ring (z = −0.225)** and **P(65 + j) on the top ring** (`test_geo38_f_clinical_matrix.py:197`, `:267`; log `:12060` z sign −1). **Feeds: P33 (φ = 5.625°) and P41 (φ = 95.625°)**, both bottom ring, 8 slots = 90° apart (`test_geo38_f_clinical_esr_root.py:65`, `:74–75`; log `:12060`).

*In our mesh but not modelled in AED, by design:* (i) 32 floating leg-port boxes 4 × 4 × 21.43 mm at the mid-leg azimuths (`test_geo38_f_clinical_build.py:55`; inert vacuum); (ii) a rod-shaped sub-volume (tag 403), radius 1.0e-3, cylindrical length 0.030 on ẑ, centred at **(0.105, 0, 0)** (`tests/validation/test_geo36_rod_f_clinical.py:88–90`; `tests/mesh/test_birdcage_implant_rod.py:93–94`), **carrying the gel's material** (`test_geo36_rod_f_clinical_drive.py:248–256`; log `:12064`) — a refinement artefact of entry 45's mesh, not an implant. The AED gel is one homogeneous box.

## Materials

| Region | σ [S/m] | εᵣ | μᵣ | ρ [kg/m³] |
|---|---|---|---|---|
| Air inside the can | 0 | 1 | 1 | — |
| Coil (legs ∪ rings) | **copper 5.8e7** (`test_geo38_f_clinical_matrix.py:58`): ours Leontovich `Z_s = (1 + j)R_s` on the conductor surface, R_s = 2.087e-3 Ω at 64 MHz (`src/fem_em_solver/core/cavity.py:543–546`; `test_th14_birdcage_copper.py:160–171`); **AED: Finite Conductivity σ = 5.8e7, Solve Inside off** (`ANS-6`'s treatment) | 1 | 1 | — |
| Gel box | **0.5** | **78** | 1 | **1000** (`test_lossy_sphere_fullwave.py:92–93`; `test_birdcage_b1_plus_map.py:139`; log `:12052`) |

## Excitation — all 64 gaps capacitors with ESR, two feeds in quadrature

- **Capacitor on every one of the 64 gaps: C = 185.7629 pF exactly** (SPEC v2 — `xl` entry 56's degree-2 root; the degree-1 record column below ran at 193.5996 pF, log `:12063`) **with a 0.2 Ω series ESR** (`test_wf10_f_human_copper_esr_root.py:53`, binding ruling 2026-10-06), `Z_RC = 0.2 − 13.386935347j Ω` at 64 MHz (review arithmetic, ω = 2π · 64 MHz; the module's form `:188–191`). In HFSS as the parallel equivalent at 64 MHz, high-Q form as `ANS-8` v4: **R_p = 1/(ω²C²r) = 896.1 Ω ‖ 185.7629 pF**, Q_C = 66.9 (exact C_p = 185.7215 pF, 0.041 pF below). Enter C with all its digits and its unit, `185.7629pF`. *(v1's figures at the degree-1 node, for the record column only — not to enter: R_p = 825.0 Ω ‖ 193.5996 pF, Q_C = 64.2, C_p = 193.5527 pF, `Z_RC = 0.2 − 12.845047127j Ω`.)*
- **62 non-feed gaps:** Lumped RLC sheets (R_p ‖ C).
- **Two feed gaps P33, P41:** 50 Ω lumped ports **with the same R_p ‖ C present across each feed gap in the solved field** — ours drives each feed through `Z_th = 50 ‖ Z_RC = 3.505604063 − 12.39875443j Ω` at the v2 node with V_th = 1 V (review arithmetic from the module's form `:144–151`; the degree-1 record column's value was `3.259458821 − 11.95984969j Ω`, log `:12063`). How the feed-gap R_p ‖ C is realised is the operator's choice; state it in the results file.
- **Drive (ii), primary: ccw = (V_33, V_41) = (1, −1j)** — P41 lags P33 by 90°; the sense that co-rotates with `B₁⁺` (`quadrature_phase_weights`, `test_birdcage_b1_quadrature.py:175–188`; `test_wf10_f_human_tuned_b1_sar.py:230–247`). **cw = (1, +1j)** is printed beside it at c\* only.
- **Normaliser: P_acc** = sum over the two feeds of HFSS's accepted power at the lumped port (incident − reflected at the 50 Ω plane), the 2026-10-12 ruling (`ANS-8` SPEC `:291–304`); ours is the lumped line cut `|V_s|²/(8·50)·(Σ|a|² − Σ|b|²)`, a = (1, −1j), b = S a (`test_geo38_f_clinical_box_loaded.py:376–381`). P_box (gel loss) and P_src (field partition) are printed beside it.
- **No re-tune on the AED side**; state so in the results file. Im Y_in at the node is a record, not a root check: the root was found on C4's cylinder stand-in and the box detunes the coil (`xl-pending.md` entry 52, Readout).

## Frequency, order and mesh

- **64 MHz only** (`test_geo38_f_clinical_esr_root.py:70`; log `:12055`). No sweep.
- **First Order** (operator order ruling 2026-10-01) — our counterpart is C19 at degree 2. A Zero Order column (our degree 1, entry 52) is optional and a record.
- Adaptive at 64 MHz to **ΔS ≤ 0.002** on the 2×2, ≥ 2 consecutive converged passes, max 10, ≥ 2 elements across each gap along φ̂; report if the cap is hit. Mixed Order forbidden. Our mesh: 1 217 025 tets, conductor 1.6 mm, global 15 mm (log `:12055`, `:144`). Resources: 20 cores / 1 TB (`ANS-8` SPEC `:158`).

## Quantities to export (all printed digits; per √W or per W of P_acc where stated)

| Row | Quantity |
|---|---|
| 1 | `S` 2×2 (P33, P41) at 50 Ω, 64 MHz, and `Y_in` at each feed under drive (ii) (sign and value of Im) |
| 2 | **`\|B₁⁺\|` at the 18 lattice points** (metres): x ∈ {−0.105, 0, +0.105}, y = 0, z ∈ {−0.25, −0.15, −0.05, +0.05, +0.15, +0.25} (`test_geo38_f_clinical_box_loaded.py:53`, `:116–127`); `\|B₁⁺\| = \|B_x + jB_y\|/2` (`test_birdcage_b1_quadrature.py:214–228`); point values per √W of P_acc |
| 3 | **SAR10g at the 18 points**: ball radius **0.013365** (10 g at ρ = 1000, `src/fem_em_solver/post/sar.py:214`; a C95.3 cube is admissible — state which), **c\*** = the largest with its (r, dz, az), per W of P_acc; the box-wide 10 g maximum and location (record) |
| 4 | Gel-volume mean `\|B₁⁺\|` over the box (AED: volume mean), per √W of P_acc — read against TBD 2 |
| 5 | Power partition per P_acc: gel / copper / capacitor-ESR (62 + 2 feeds) fractions |
| 6 | mean `\|B₁⁺\|` / mean `\|B₁⁻\|` over the gel (the sense check) |
| 7 | Solve metadata: tets, passes, final ΔS, converged yes / no, order, unknowns, wall time; whether AED was re-tuned (must be no) |

## Readout (pre-registered, adjudicated on C19 only; TBD 3)

- Row 1: Im Y_in of the **same sign** at both feeds ⇒ record (C19 vs entry 52
  gives our order sensitivity first); a sign flip is a detuning finding.
- Row 2: within **10 %** per point with the z-ordering at x = 0 preserved
  (ours: `|dz| = 0.05 > 0.15 > 0.25`) ⇒ AGREE.
- Row 3: c\* within **15 %** and within one ball radius ⇒ AGREE; the box-wide
  maximum is a record. **Tie clause (v2.1, 2026-11-09 weekly):** our 18-point
  SAR10g column has near-ties — on the degree-1 column the top two balls are
  mirror balls 210 mm apart with SAR within ≈ 6.3e-04 relative (`:12291`,
  `:12303`), and the 3 T twin's c\* landed on the other member of its pair
  (`xl` entry 58). So the location half is met if AED's c\* lies within one
  ball radius of **any** lattice ball whose SAR10g is within **1 %** of our
  maximum on the adjudication column (the tie set, listed in
  `COMPARISON.md`), and the value half compares AED's SAR10g(c\*) with that
  ball's. A c\* outside the tie set is the location finding.
- Row 5: gel fraction within **5 %**, copper and ESR within **10 %** ⇒ AGREE.
- Row 6: a ratio ≫ 1 in both codes (ours 2.96) is the sense check.
- No absolute tuning, compliance, homogeneity or C95.3 claim follows.

## Our column (degree 1, `xl` entry 52 — a record; raw values at V_th = 1 V per feed)

`20261103T080009Z_GEO-38-C6.log`, `-n 8`, 868 s, 40.95 GiB; every gate green
((R1) 1.578e-14, (P) 1.802e-12 / 1.425e-12, coverage (i) 1.976e-14 / 4.663e-15;
build 460.8 s, column 219.0 s, readouts 169.7 s, `:12312`).

| Quantity | Value | Log line |
|---|---|---|
| S₁₁, S₁₂ = S₂₁, S₂₂ | −4.8167140712e-01 + 7.9527084190e-04j; 4.6496927594e-03 + 1.0179009021e-01j; −5.0211295667e-01 − 2.1994222327e-02j | `:12281` |
| Y_in ccw, P33 / P41 [S] | 4.450128950e-02 + 4.009154872e-04j / 8.079211370e-02 + 4.413545175e-03j | `:12285` |
| mean `\|B₁⁺\|`_ccw over 400 gel-cell centroids / mean `\|B₁⁻\|` | 2.090267723e-08 T / 7.053675108e-09 T; ratio 2.963374; CV 0.392744 (cw Y_in: `:12286`) | `:12288` |
| 18-point lattice mean `\|B₁⁺\|`_ccw | 1.49792e-08 T (arithmetic on `:12289–12306`, not printed) | — |
| c\* | x −0.105, z −0.050 (r 0.105, az 180°): SAR10g 1.378896513e-04 W/kg (ccw); cw ball 1.059925022e-04 | `:12307` |
| P_box / P_src / P_acc [W] | 3.204417752e-03 / 8.266643116e-02 / 6.063124201e-02 | `:12310` |
| SAR10g(c\*) per P_acc; mean `\|B₁⁺\|`/√P_acc | 2.274234e-03 W/kg/W; 8.488944e-08 T/√W (lattice mean: 6.0833e-08, derived) | `:12310` |

Lattice `|B₁⁺|`_ccw [T] / SAR10g [W/kg] (`:12289–12306`; six digits here, full digits in the log):

| x | dz −0.25 | −0.15 | −0.05 | +0.05 | +0.15 | +0.25 |
|---|---|---|---|---|---|---|
| −0.105 | 5.319736e-09 / 7.896193e-05 | 1.664794e-08 / 1.153728e-04 | 2.230737e-08 / 1.378897e-04 | 2.255457e-08 / 1.368278e-04 | 1.691890e-08 / 1.118693e-04 | 4.981575e-09 / 7.397515e-05 |
| 0 | 5.345957e-09 / 1.106314e-04 | 1.526970e-08 / 7.066740e-05 | 2.084487e-08 / 1.150657e-05 | 2.082999e-08 / 9.807226e-06 | 1.528202e-08 / 6.505141e-05 | 5.246063e-09 / 1.038151e-04 |
| +0.105 | 5.140520e-09 / 8.277083e-05 | 1.895999e-08 / 1.160654e-04 | 2.519919e-08 / 1.378033e-04 | 2.486936e-08 / 1.349824e-04 | 1.842140e-08 / 1.106522e-04 | 5.485955e-09 / 7.651317e-05 |

## Our side

The adjudicated column is **C19** — C6's module at degree 2 on the same mesh
(`xxl`, `xl-pending.md:2182`, `:2335`; node per TBD 1) — read through the same
rows; the degree-1 table above is a labelled record, as `ANS-8`'s degree-1
rows are. `COMPARISON.md` carries our rows with the AED columns blank.

**AED / Ansys values are private (CLAUDE.md hard rule): they go in the
gitignored `aed_results/` and `COMPARISON_private.md` only; the qualitative
verdict and a relative agreement level are all that reaches a tracked file.**

**Out of scope:** no implant (`GEO-36` step 8 / `IMP-6` step 4 are separate
cases); no thermal half (chain P, entry 59); no matching network; no 3 T twin
(entry 58); no homogeneity or compliance statement; no claim at degree 1.

## Revision history

- **v1 (2026-11-04, weekly planning review; skeleton):** fixture, materials,
  excitation, rows and the entry-52 record as above; TBDs 2 and 3 ruled the
  same night (the lattice mean is the comparison row; `ANS-8`'s bands);
  TBD 1 (the node) is v2's, from `xl` entry 56.
- **v2.1 (2026-11-09, weekly planning review; readout only):** row 3 gains
  the tie clause (the location half against the tie set of lattice balls
  within 1 % of our maximum); nothing the operator builds changes.
- **v2 (2026-11-07, weekly planning review; the node):** C = 185.7629 pF
  (entry 56's degree-2 root at 0.25 pF, first branch), R_p = 896.1 Ω ‖ C at
  64 MHz (Q_C = 66.9, C_p = 185.7215 pF), `Z_RC` and `Z_th` recomputed;
  everything else unchanged. Ready to build; our adjudication column is
  `xxl` entry 62 (2026-11-14).
