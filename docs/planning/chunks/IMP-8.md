# IMP-8 — the heated rod in gel: the finite line source in closed form

*History file for the §7 `IMP-8` row (`OPS-46`). Opened by the 2026-11-02 10:30 daily review under §9 rule (5) lane 3 (a Phase 7 chunk whose validation target is a cited closed form); **for the weekly to ratify**. Not a chain step, so it takes no count on chain I.*

**Why this chunk.** `IMP-7` is the implant-heating canonical on a sphere, and I36 and entry 55 read a rod in the clinical box. **Nothing in the repo gives a closed form for the temperature around a rod,** which is the F2182 geometry: `erfc` and `exp1` appear in no validation module and not in `src/`. In gelled saline there is no perfusion, so the infinite-medium solution is exact until the walls are felt.

**The closed form.** Carslaw & Jaeger, *Conduction of Heat in Solids*, 2nd ed., ch. X, the continuous point source and the continuous line source. **The text is not in the repo, so the chapter citation is from memory; the derivation below is what the module rests on.** A point source of power P switched on at t = 0 in an infinite medium gives θ(R, t) = P / (4πkR) · erfc(R / (2√(κt))), the time integral of the instantaneous source. So a uniform line of power q′ per length on |z′| ≤ L/2 gives

θ(r, z, t) = q′ / (4πk) · Φ, with Φ(r, z, t) = ∫ erfc(R / (2√(κt))) / R dz′ over |z′| ≤ L/2, and R² = r² + (z − z′)².

**Its limits.** As t → ∞, Φ∞ = asinh((z + L/2) / r) − asinh((z − L/2) / r). As L → ∞ at z = 0, Φ → E₁(r² / (4κt)).

**Done-when.** Steps 0 and 1 green at their gates, each module asserting its limits and controls. The FEM half is the weekly's to pose.

### §9 item 759 — 2026-11-02 (lane 3, closed form; `IMP-8` step 0)

**Outcome (2026-11-02 13:30 slot): DONE, attempt 1, `IMP-8` ⬜ → 🟡.** `20261102T193113Z_IMP-8-item759-W1.log`: Status 0, 7 passed, 4 s, `-n 1`, the image's default build (no solve, nothing scalar-type dependent). (Inf) 1.110e-16; (St) 4.586e-08 / 1.945e-07 / 5.806e-07 / 2.758e-07 against 1.329e-07 / 5.602e-07 / 1.190e-06 / 3.628e-07; (En) 2.220e-16 / 0; (Far) 4.162e-04 against 4.167e-04; (Mono) 0.3761431 … 4.689659, Φ∞ 5.424931. Controls, measured beside predicted: (c-κ) 0.9544 / 0.954, (c-fin) 9.538e-04 / 9.5e-04, (c-L) 0.2546 / 0.254. Printed at 900 s, (2, 0) and (2, 15) mm: Φ 4.006102 / 2.134933, Φ∞ 5.424931 / 3.402307, K per watt 17.71086 / 9.438477 (steady 23.98346 / 15.0415). `quad`, split at z′ = z, `limit=800`, reached every bar on the first window, so the review's unlogged arithmetic stands confirmed.

**Anchor.** New module `tests/validation/test_imp8_rod_line_source_closed_form.py`, numpy and `scipy.integrate.quad`, with L = `ROD_LENGTH` (`tests/mesh/test_birdcage_implant_rod.py:94`, 30 mm), `K_GEL` (`test_mat15_coil_driven_heating_phantom`, 0.6), `RHO_C` (`test_mat12_transient_sar_heating_sphere`, 4.18e6), κ = k / ρc, and `B1_BAND` from `test_imp4_pec_cylinder_closed_form` (1e-12). All imported. `Φ` is defined in the module; **no `src/` edit**.
- **(Inf)** r = 2 mm, z = 0, t = 10 s: |Φ / E₁(r² / 4κt) − 1| ≤ `B1_BAND`. Review's arithmetic ≈ 1e-16.
- **(St)** at t = 1e9 s, at (r, z) = (2, 0), (2, 16), (5, 30), (20, 0) mm: the residual Φ∞ − Φ equals its first term L / √(πκt) to a relative bound **x² / 3** at x = R_max / (2√(κt)), the series' own next term, computed in the module. Review's arithmetic: rel 4.6e-08 / 1.9e-07 / 5.8e-07 / 2.8e-07 against bounds 1.3e-07 / 5.6e-07 / 1.2e-06 / 3.6e-07.
- **(En)** the energy identities behind both sources: ∫₀^∞ E₁(s) ds = 1 (the line's q′t) and ∫₀^∞ u erfc(u) du = 1/4 (the point's Pt), each ≤ `B1_BAND`.
- **(Far)** steady, r = 10 L, z = 0: |Φ∞ / (L / r) − 1| ≤ (L / r)² / 24, the alternating series' first term. Review's arithmetic ≈ 4.16e-04 against 4.17e-04.
- **(Mono)** Φ at (2, 0) mm strictly rises on t = 10, 30, 100, 300, 900, 3600 s and stays below Φ∞.

**Negative controls.**
- **(c-κ)** κ → κ / 4 inside the erfc misses (Inf) by > `B1_BAND`. Ceiling, review's arithmetic: 0.954.
- **(c-fin)** the (Inf) comparison at t = 120 s misses by > `B1_BAND`, 9.5e-04 by the review's arithmetic, so the finite length is visible and (Inf) is not a tautology.
- **(c-L)** L → 2L in Φ∞ misses the t = 1e9 Φ at (2, 0) mm by > 0.1. Review's arithmetic ≈ 0.254.

**Printed (no band):** Φ and Φ∞ at t = 900 s at (2, 0) mm and (2, 15) mm, and the same as K per watt of total power, Φ / (4πkL). Review's arithmetic for Φ: 4.006 / 5.425 and 2.135 / 3.402. **Print 7 digits or fewer.**

**The review ran this arithmetic with stdlib `math.erfc` and numpy Gauss–Legendre, unlogged, for ceilings only.** The module's own values are the record. **`quad` was not run** (no scipy on the host), so if `quad` cannot reach a bar at (Inf) or (St), say so with the value. Fix the quadrature first (split at z′ = z, `limit=800`), never the bar.

**Tier, ranks, wall clock.** Standard, `-n 1`, `timeout -k 30 170`; state the build. Predicted seconds. **Class:** lane 3. **Cost:** ≈ 2.7 measured / 20 predicted min. **Depends on:** nothing. **Family:** `IMP-8` step 0, attempt 1 of 4.

**Traps already paid for.** `-k a or b` splits inside a quoted container command. pytest captures prints without `-s`. The imports pull dolfinx, so the module runs only in the container. Every point here has r > 0, and `quad` must be split at z′ = z when |z| < L/2.

**Scope.** A line of zero radius in an infinite homogeneous medium, the rod's own heat capacity ignored, uniform power. **Not claimed:** any F2182 figure, a rod's SAR shape, a wall effect, or a Robin box. **The rod's radius is 1 mm** (`ROD_RADIUS`), so r = 2 mm is 1 mm off its surface; the line is the far-field comparand.

**Negative result.** A red at a quadrature floor: fix the quadrature and say so. A red on the closed form itself is a wrong derivation in this brief: stop, keep the log, BLOCKED; never move a bound. No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` ⬜ → 🟡 on step 0.

### §9 item 760 — 2026-11-02 (lane 3, closed form; `IMP-8` step 1)

**Outcome (2026-11-02 13:30 slot): DONE, attempt 1, `IMP-8` 🟡 → ✅ on the closed form.** `20261102T193249Z_IMP-8-item760-W1.log`: Status 0, 12 passed (step 0's seven unchanged, five new), 3 s, `-n 1`. (Rate) rel 8.929e-03 / 2.215e-03 / 5.527e-04 / 1.381e-04, slope 2.0047 in the provisional [1.95, 2.05], **for the weekly to ratify**; (Sym) 0 and 0; (Far), named `(Far2)` in the module to keep it apart from step 0's, miss 1.871e-06 against 1e-4 with Φ∞ the larger; (c-one) 0.4994 against predicted 0.5. Printed at 900 s: (2, 15) mm Φ₂ 6.787982, Φ 2.134933, ratio 3.179482; (2, 0) mm Φ₂ 0.6868772, Φ 4.006102, ratio 0.1714577. The FEM half stays the weekly's to pose.

**Why.** A rod's SAR concentrates at its tips, not uniformly along it. Two tip sources and the uniform line are the two clean shapes with the same total power, and the ratio between them at a probe is what a lumped reading of a rod cannot give.

**Anchor.** Same module, `Φ` from step 0, with the tip pair Φ₂(r, z, t) = (L / 2) Σ_± erfc(R_± / (2√(κt))) / R_±, R_±² = r² + (z ∓ L/2)², the same total power.
- **(Rate)** the N-point midpoint sum of point sources against Φ at (20, 0) mm, t = 900 s, N = 4, 8, 16, 32: LSQ slope of log error in [1.95, 2.05]. **That band is this review's, provisional, for the weekly to ratify.** Review's arithmetic: rel 8.93e-03 / 2.22e-03 / 5.53e-04 / 1.38e-04, slope ≈ 2.00.
- **(Sym)** Φ and Φ₂ even in z at (2, ±10) mm, ≤ `B1_BAND`.
- **(Far)** steady, a = L / r = 0.1, z = 0: | |Φ₂∞ − Φ∞| / a − a² / 12 | ≤ a⁴, the two expansions' difference and its next-order bound. Φ∞ is the larger.

**Negative control.** **(c-one)** dropping one of the pair's two terms misses (Far) by ≈ 0.5 against a bound of 1e-4, since the pair is symmetric at z = 0.

**Printed (no band):** Φ₂ and Φ at t = 900 s at (2, 15) mm and (2, 0) mm and their ratio. Review's arithmetic: 6.788 / 2.135 = 3.18 and 0.6869 / 4.006.

**Tier, ranks, wall clock.** Standard, `-n 1`, `timeout -k 30 170`. **Class:** lane 3. **Cost:** ≈ 2.7 measured / 15 predicted min. **Depends on:** item 759 landing; skip at no cost if it has not. **Family:** `IMP-8` step 1, attempt 1 of 4.

**Traps already paid for.** As item 759.

**Scope.** Two shapes of a source. **Not claimed:** which shape a real rod has, a tip Q, or anything about I36's η_rod. Entry 55's tip θ is the weekly's to read, and no comparison with it is printed here.

**Negative result.** As item 759. No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` 🟡 → ✅ on the closed form, the FEM half left for the weekly to pose.

**Review note 2026-11-03 03:00 (from the 11-02 weekly's rulings, §10).** Ratified, with **(Rate) [1.95, 2.05]** (the N-point midpoint sum of point sources is a second-order quadrature of a smooth integrand; 2.0047 is its slope). The Carslaw & Jaeger chapter citation is from memory and stays marked so until the operator's print copy makes it firm (dashboard item 0). **Done-when narrowed:** "each with a control that fires" now reads "on each gate that has content at finite t" — Φ₂ at finite t has no closed-form gate and (Sym) is near-tautological, so Φ₂'s first gate is the FEM (I56). Widening it back was declined. The FEM half as posed: step 2 = I55 (the finite line source in FEM on a ≥ 0.10 m half-width gel box, (I) at `INTERIOR_BAND`, (E), (S0), (St); controls L → 2L and a perturbed source), step 3 = I56 (two balls at z = ±L/2, Φ₂'s first gate; one ball dropped the control), step 4 = I59 (Duhamel duty cycle, closed form), then I60 (I59's FEM twin, serial on I55 + I59), I61 (an `imp:` example, serial on I55) and I58 (the lumped reading of entry 48's rod power, numpy). The ratios 3.179482 / 0.1714577 stay printed, no band.

### §9 item 778 — 2026-11-03 (§10 I55, chain I — `IMP-8` step 2, the finite line source in FEM; chain step (solve))

**Outcome (2026-11-03, slot started 10:30 UTC): BLOCKED, attempt 1. W1 was killed at its tier ceiling, so no gate was read. This is not the predicted red.** Module parked on `attempt/IMP-8-I55-20261103T1111Z`.
- **W0** (`20261103T105052Z_IMP-8-I55-W0.log`) was killed at 300 s inside gmsh. That is not a result.
- **W0b, a cost probe** (`20261103T105626Z_IMP-8-I55-W0b-probe.log`, HXT mesher by an env knob, 67 s) measured **4 101 108 cells**, red on the 800 000 cap. The cause is in the fixture: beyond its DistMax a gmsh `Threshold` returns its SizeMax unless `StopAtDistMax` is set. So the inner field returned 2 mm everywhere, and the `Min` capped the whole box at 2 mm.
- **W0c** (`20261103T105748Z_IMP-8-I55-W0c.log`, gmsh default mesher, `StopAtDistMax = 1` on the inner field, 23 s, Status 0) measured **485 291 cells** and 6 495 rod cells. Selected volume 9.465847e-08 m³, so the effective radius is **1.0022 mm** against 1 mm. Mesh size h within 3 mm of each point: (10, 0) 2.131 mm, (20, 0) 3.245, (10, 15) 2.136, (20, 15) 3.330, (2, 0) 1.167.
- **W1** (`20261103T105819Z_IMP-8-I55-W1.log`, `timeout -k 30 600`) was **killed at 600 s**. Only the W0 build test printed. The CN solve (degree 2, N = 40) never returned, and the two steady solves were still ahead of it. After the kill `pgrep -c python3` read 0.
- **Unblock:** re-price W1 on the measured mesh. At ≈ 485 k P2 cells, one transient and two steady solves do not fit 600 s at -n 2. The options are to split the steady (St) / (c-E) solves into their own window, pose a coarser inner grading, or move the run to an `xl` window. That is the review's call. Shrinking the mesh in-slot would re-pose the fixture.

**Executes** §10 I55, with one fixture premise replaced. **`ThermalBoxGeometry` (`study.py:1678`) is a uniform box.** At 0.10 m half-width it cannot resolve a 1 mm rod: `MAT-21` step 1's finest rung is h ≈ 4.5 mm on a 0.045 m box. The execution detail is a **graded gmsh box**; the target and the gates stay as the weekly wrote them.

**Fixture.** New `tests/validation/test_imp8_rod_line_source_fem.py`, real build.
- **The builder.** A module-level `_graded_box(comm, *, half_width=0.10, axis_half=L/2 + 2e-3)`, which I56 and I60 import. It is a gmsh OCC box [−0.10, 0.10]³ m: one volume (cell tag 1), outer facets tag 1 (Robin), **no embedded rod**.
- **The grading.** A gmsh `Distance` field on ≥ 60 points sampled along the axis segment |z| ≤ axis_half, x = y = 0, feeds two `Threshold`s under a `Min`. Inner: 0.4 mm within 1.5 mm of the segment, rising to 2.0 mm at 25 mm. Outer: 2.0 mm rising to 12 mm at 100 mm.
- **W0 prints the cell count and stops if it exceeds 800 000** (estimate 150 000–400 000).
- **The rod source.** A DG0 indicator of cells whose midpoint has r ≤ `ROD_RADIUS` 1e-3 and |z| ≤ L/2 (the `_interior_chi` pattern, `test_mat21_box_fem.py:111–121`). It is **normalised so the owned-volume sum gives ∫q dV = P exactly** (`test_mat21_box_dg0_source_fem.py:152–180`), with P = 1 W (results per watt) and q′ = P/L. Print the effective radius implied by the selected cells' volume beside 1 mm.

**Solve.** `solve_transient_conduction` (`src/fem_em_solver/thermal/transient.py:94`):
- k = `K_GEL` 0.6, ρc = `RHO_C` 4.18e6, h = `H_GEL` 10, t_inf 0, θ₀ = 0.
- **Degree 2**, θ = 0.5, `startup_half_steps=2`, Δt = t\*/40 (`N_MAIN`), t\* = `T_PRINT` 900 s.
- Sample points with `evaluate_vector_field_parallel` (`src/fem_em_solver/post/evaluation.py:224`).

**Gates (imported, unmoved).**
- **(E):** `max_balance_residual` ≤ `E_BAND` 1e-9.
- **(S0):** the source's ∫q = P within `B_BAND` 1e-10, with the selected-cell count > 0 on every run.
- **(I):** θ_FEM against q′/(4πk)·`_phi(r, z, t*)` (`test_imp8_rod_line_source_closed_form.py:63`) at (r, z) = (10, 0), (20, 0), (10, 15) and (20, 15) mm, rel ≤ `INTERIOR_BAND` 1e-2. Print (a/r)² beside each point; (2, 0) mm is printed only.
- **(St):** at t = `T_STEADY` 1e9 s, against `_phi_steady` at the same points. **This one is pre-registered red.** The comparand is the infinite medium, but the box is Robin-cooled at 0.10 m. Review arithmetic puts the box's steady offset at ≈ −0.9 P (K per W), against θ∞ ≈ 10.6 P at (10, 0) mm and 6.1 P at (20, 0) mm: a miss of ≈ 0.09–0.15. At t\* the walls are ≈ 9 diffusion lengths away (√(κt\*) ≈ 11.4 mm), so (I) is unaffected.

PREDICTED RED: (St) range [0.03, 0.30]

**Controls.**
- **(c-2L):** L → 2L in the comparand misses (I) by > 0.1 at (10, 15) mm.
- **(c-E):** I32a's `PERTURBED_SOURCE` pattern (K10-es). Solve with q × (1 + 1e-3), and balance against the nominal P. The balance misses by ≈ 1e-3, ≫ `E_BAND`. The scheme's own residual cannot see a consistent perturbation, so this control checks the bookkeeping.

**Windows.** W0: `-n 2`, build only (cells, h at the four points, the selected-cell count), `timeout -k 30 300`. W1: `-n 2`, the gates, `timeout -k 30 600`. Window form: `20261101T151011Z_MAT-21-P54-W1.log:12` (real, `-n 2`, 194 s), with the module changed. Pre-check `pgrep -c python3` reads 0 before each window.

**Tier.** Heavy, `-n 2`, real. **Class:** chain step (solve). **Cost:** ≈ 6.1 measured / 45 predicted min.

**Traps already paid for.**
- MAT-21's box helpers assume `create_box`. A gmsh mesh goes through `gmshio.model_to_mesh` with the cell and facet tags, as `two_layer_ball_domain` does (`src/fem_em_solver/io/mesh.py`).
- `.eval` bytes depend on the batch, so no bitwise comparison here.
- A `SpatialCoordinate` facet integral without a pinned quadrature degree can hang FFCx: pin it.
- A killed compile leaves an orphan `.c`. Sweep only when `pgrep` reads 0.

**Scope.** One source shape, one t\*, one mesh. **Not claimed:** a rate, the two-tip twin (I56), a duty cycle (I60), anything about entry 48's rod.

**Negative result.**
- (I) red: stop, keep the log, BLOCKED. The mesh near the points is the first suspect; print h there.
- (St) outside its predicted range: the item stops, as on any red.

**Status it can move.** `IMP-8` gains step 2, a 🟡 record line on its ✅ row and the first FEM gate on Φ. §10 I55 closes, and I56 / I60 unblock. **Depends on:** nothing.

### §9 item 786 — 2026-11-03 (§10 I58, chain I — `IMP-8` step 6, the lumped reading of entry 48's rod; lane 3 → chain step)

**Executes** §10 I58, with two premises made exact.
1. **The tip 1 g ball centre lies on the line's axis, at its end.** Entry 48's top and bottom balls sit at the rod centre ± L/2 (`_probe_points`, `test_geo36_rod_f_clinical_drive.py:273–292`). There Φ and Φ₂ are singular: 1/R is integrable in 3D, but not at a point. So the lumped θ at a tip is read in two ways, both printed:
   - **as the 1 g ball average**, with radius `averaging_ball_radius(mass_kg=1e-3, rho=…)` at the gel ρ (`src/fem_em_solver/post/sar.py:214`);
   - **at the rod surface point** (r = `ROD_RADIUS`, z = L/2).
2. **(En) and (Far) do not take P** (`test_imp8_rod_line_source_closed_form.py:118–136`). They are re-run as imported, as the kernel's own gates. This item's own gates are on the new quadrature.

**Inputs, printed and taken from the logs.**
- P_rod = 9.495544579e-08 W on (1, +1j), arm M (`20261102T080008Z_GEO-36-step9-I36.log:12494`).
- The rod centre (0.105, 0, 0) m, its tips at z = ±0.015 m and its apexes at ±0.016 m (`20261031T003859Z_GEO-36-I36c-W0b.log:90`).
- P37's six points (`PROBE_PTS`, `test_mat21_box_probe_points_fem.py:72–79`), mapped into the clinical box by `_map_box_point` (`test_geo36_rod_f_clinical_heating.py:145`). Read the box centre from the module; never re-derive it.

New module `tests/validation/test_imp8_entry48_lumped.py`, numpy / scipy.

**Gates.**
- **(Bq):** the ball average of Φ by a product Gauss rule (radial × angular) agrees at n and 2n points to ≤ 1e-6 relative. Beside it, the ball average of the point kernel 1/R over a ball centred on the source equals its closed form 3/(2a) to ≤ 1e-6. That is the singular-integrand check that gives (Bq) its meaning.
- **(Mono):** θ(t\*) < θ∞ at every point (a monotone rise).
- (En) and (Far), imported.

**Control.** Moving the 1 g ball one ball radius off the axis changes the ball-averaged value by > 1e-2 against the on-axis one. The readout is not flat, so the gate set can see position.

**Printed, as a record.** θ_line and θ_tips at t\* = 900 s in K:
- at the tip ball (average and surface point), for the top and bottom tips;
- at the six box points, where it is ≈ 0 (√(κt\*) ≈ 11 mm against ≥ 0.1 m) and printed as such;
- the ratio θ_tips / θ_line at the tip ball;
- all per watt and at P_rod.

It is a record for entry 55 to be read against. 🟡.

**Tier.** Smoke, `-n 1`, `timeout -k 30 60`. **Class:** lane 3 → chain step. **Cost:** ≈ 2.1 / 15.

**Traps already paid for.**
- `_phi`'s `quad` needs its split at z′ = z for points near the axis.
- **The leak guard:** this item prints a log-derived power. Print P_rod at the log's own `.9e` form, never a reformatted digit string.
- Never open `docs/private/`.

**Scope.** An infinite-medium lumped shape at entry 48's power. **Not claimed:** the box walls, the rod's own conduction, the FEM tip θ (that is entry 55).

**Negative result.** A gate red: stop, BLOCKED.

**Status it can move.** §10 I58 (a record for reading entry 55) and the `IMP-8` row's I58 clause. **Depends on:** nothing.

**Outcome (2026-11-03, 15:00 slot): DONE, attempt 1, a 🟡 record.** `20261103T211525Z_IMP-8-item786-W1.log`: Status 0, 9 s, `-n 1`, the image's default (real) build, 8 passed / 3 skipped (the complex-only environment tests). New module `tests/validation/test_imp8_entry48_lumped.py`, numpy / scipy; `_phi`, `_phi_steady`, `_phi_tips` imported from step 0's module; no `src/` edit.
- **The rule.** Product Gauss–Legendre in (ρ, μ, φ) about the ball centre, ρ = a·u² clustered at the centre and 1 ∓ μ = s⁴ at each pole (the axis line's log singularity), n = 16 against 2n = 32. 1 g ball radius a = 6.203505e-03 m at ρ = 1000 kg/m³ (`PHANTOM_RHO_KG_PER_M3`).
- **(Bq):** line, rel 8.078e-10 at both tip balls; tips, 2.109e-15 / 2.331e-15; band 1e-6. **(Bq-k):** 1/R over a ball centred on its source against 3/(2a) = 2.417987931e+02, rel 8.882e-16 (n) / 1.776e-15 (2n).
- **(Mono):** θ(t\*) < θ∞ at both tip balls, both surface points and the six box points, both shapes. Tip-ball steady Φ 2.909416 (line) / 4.126982 (tips).
- **(En) / (Far), imported:** 2.220e-16 / 0, and 4.162e-04 against 4.167e-04, as step 0.
- **Control (c-off):** the ball moved one radius off the axis: 1.652975 → 1.049444, change 3.651e-01 > 1e-2. Off-axis rule at n = 8, 16 φ nodes on [0, π].
- **Inputs, asserted against the logs:** P_rod printed as the log string 9.495544579e-08 W; centre (0.105, 0, 0), tips ±0.015, apexes ±0.016 m from `_probe_points(i35._rod_centre())`.
- **Printed (t\* = 900 s), top = bottom by symmetry:**
  - tip 1 g ball: Φ 1.652975, Φ₂ 2.926267; θ_line 7.307753 K/W = 6.939109e-07 K at P_rod; θ_tips 12.93694 K/W = 1.228433e-06 K; θ_tips / θ_line **1.770303**;
  - rod-surface point (r = 1 mm, z = ±L/2): θ_line 12.48967 K/W = 1.185962e-06 K; θ_tips 63.16158 K/W = 5.997536e-06 K; ratio 5.057107;
  - P37's six points (`_p37_fixture(centre[1])`): ≈ 0, largest θ_line 7.093e-11 K/W (6.735e-18 K) at P0 and P1, both 0.105 m off the axis; steady 0.39–1.26 K/W.
- **Disclosed:** `quad` raises an `IntegrationWarning` (roundoff) inside `_phi` at ball nodes nearest the axis; (Bq) at 8.1e-10 shows it does not reach the reading.

### §9 item 787 — 2026-11-03 (§10 I59, chain I — `IMP-8` step 4, the line source under a duty cycle, closed form; lane 3 → chain step)

**Executes** §10 I59. New `tests/validation/test_imp8_rod_line_source_duty_closed_form.py`, numpy / scipy only.
- **The step response:** θ_step(r, z, t) = q′/(4πk)·`_phi(r, z, t)`, with P = 1 W.
- **The duty sum:** MAT-17's Duhamel rule (`_duhamel`, `test_mat17_duty_cycle_closed_form.py:92–108`) re-implemented on θ_step. The source switches on at jP and off at jP + DP, with period P = 0.5 t\* = 450 s and D ∈ {1, 0.5, 0.25}.
- **Where and when:** at (2, 15) and (10, 0) mm, on a 41-point time grid to t\* = 900 s.

**Gates (imported).**
- **(D1):** the D = 1 sum equals θ_step within `D1_BAND` 1e-13 relative.
- **(Duh):** as `test_mat21_wall_duty_closed_form.py:119` defines it, within `DUH_BAND` 1e-10. Read that test body first and apply the same identity to Φ.
- **(Q-d) and (Win-d):** as `test_mat21_pre_equilibration_duty_dose_closed_form.py:89–117` does them (`_q_check_g`, `_win_check_g`). `cem43` (`src/fem_em_solver/thermal/dose.py:45`) against its quadrature within `DOSE_BAR` 1e-10; `running_average` (`:108`, window 360 s) within `AVG_BAR` 1e-12.
- **The dose's temperature scale:** use the baseline and per-watt scaling that `test_mat21_box_dose_closed_form.py` uses (P39′'s convention; read it), and print it.

**Control.** Substituting the continuous form (D = 1) for D = 0.5 misses the D = 0.5 value at t\* by > 0.1.

**Printed, as a record (🟡):** θ(t\*) per watt, the six-minute-window maximum and the dose at each duty and point.

**Tier.** Smoke, `-n 1`, `timeout -k 30 60`. The `IMP-8` W1 windows took 3–4 s; `MAT-21` step 32 took 14–18 s. **Class:** lane 3 → chain step. **Cost:** ≈ 2.1 / 15.

**Traps already paid for.** `_phi` is a `quad` split at z′ = z; keep the split at r = 2 mm. Read `cem43`'s units in its docstring.

**Scope.** Closed form only. **Not claimed:** the FEM (that is I60), or a recommended duty.

**Negative result.** A gate red: stop, BLOCKED.

**Status it can move.** `IMP-8` step 4 (a 🟡 record line) and §10 I59. With item 778 landed, I60 becomes queueable. **Depends on:** nothing.

**Outcome (2026-11-03, 15:00 slot): DONE, attempt 1, a 🟡 record.** `20261103T211919Z_IMP-8-item787-W1.log`: Status 0, 3 s, `-n 1`, the image's default (real) build, 6 passed / 3 skipped (the complex-only environment tests). New module `tests/validation/test_imp8_rod_line_source_duty_closed_form.py`, numpy / scipy; `_phi`, `_phi_steady`, `KAPPA`, `QUAD`, `L` imported from step 0's module; `D1_BAND` / `DUH_BAND` from `test_mat21_wall_duty_closed_form.py`, `DOSE_BAR` / `AVG_BAR` / `W_B2` from `test_mat21_box_dose_closed_form.py`, `_q_check_g` / `_win_check_g` from `test_mat21_box_dt_duty125_dose_closed_form.py`; no `src/` edit.
- **Schedule (asserted on the grid):** P = 450 s, Δt = 22.5 s; switches D = 0.5 at 0+, 225−, 450+, 675−, 900+; D = 0.25 at 0+, 112.5−, 450+, 562.5−, 900+. q′/(4πk) at 1 W = 4.420970641 K per unit Φ.
- **(D1):** 0.000e+00 at both points (band 1e-13). Disclosed: exact by construction, since each off edge at (j + 1)P and the next on edge evaluate the identical `_phi` term; the gate checks the edge bookkeeping, nothing numerical.
- **(Duh), how it was read.** Φ has no modal duty form, so `MAT-21` step 21a's "from rest = Duhamel step superposition" is applied as: the step sum (one `quad` in z′ per edge) against the schedule convolved with the impulse kernel, ∂Φ/∂t = e^{−r²/4κt} / (2t) · [erf((z + L/2)/2√(κt)) − erf((z − L/2)/2√(κt))] (the z′ integral of ∂/∂t of the erfc integrand, done in closed form), one `quad` in τ per on-interval. Worst 2.952e-16 of θ_ss over the 41 points, three duties and two points (band 1e-10).
- **(Q-d) / (Win-d):** worst rel 4.042e-16 (bar 1e-10) / 5.720e-16 (bar 1e-12); six-minute maxima inside their window's sample range, t_max in [360, 900] s. Temperature scale: T = T_INF (20 °C, MAT-15's, `test_mat21_box_closed_form.T_INF`) + θ at 1 W; max T 29.438477 °C < 43, asserted.
- **Control (c-cont):** θ_D1(t\*) / θ_D0.5(t\*) − 1 = 3.5028 at (2, 15) mm and 1.0722 at (10, 0) mm, > 0.1.
- **Printed (t\* = 900 s, the end of an off phase; per watt; dose at T = 20 °C + θ(1 W)):**

  | D | point (mm) | θ(t\*) K/W | six-minute max K/W (t_max s) | CEM43 s |
  |---|---|---|---|---|
  | 1 | (2, 15) | 9.438477 | 8.948814 (900) | 1.571787e-06 |
  | 1 | (10, 0) | 4.645356 | 4.079805 (900) | 1.84337e-09 |
  | 0.5 | (2, 15) | 2.096145 | 5.135368 (832.5) | 3.548044e-08 |
  | 0.5 | (10, 0) | 2.241722 | 2.329229 (900) | 1.969098e-10 |
  | 0.25 | (2, 15) | 0.8986021 | 2.656731 (810) | 1.816366e-09 |
  | 0.25 | (10, 0) | 1.00437 | 1.231445 (900) | 4.96861e-11 |

  θ_ss per watt: 15.0415 at (2, 15) mm, 10.56403 at (10, 0) mm. At (10, 0) mm the D = 0.5 off-phase θ(t\*) still exceeds half of D = 1 (lag: √(κt) ≈ 11 mm against r = 10 mm).
- **For I60:** the module's `_duhamel(r, z, t, duty)` is the comparand item 799 names (θ at t\* and t₂ = 675 s), per watt.

### §9 item 796 — 2026-11-03 (§10 I56, chain I — `IMP-8` step 3, the two-tip twin in FEM; chain step (solve))

**Executes** §10 I56 on item 778's `_graded_box`, imported unchanged. Its `axis_half` already covers the balls out to ±(L/2 + 1 mm).

**Source.** Two DG0 indicators: the cells whose midpoint lies within `ROD_RADIUS` of (0, 0, ±L/2), each normalised to exactly P/2.

**Gates.**
- **(I):** against q′/(4πk)·`_phi_tips(r, z, t*)` (`:89`, i.e. (L/2)·Σ± erfc(R±/2√(κt))/R±), at item 778's four points, within `INTERIOR_BAND`.
- (E) and (S0) as item 778. There is no (St) here.

**Printed, as a record with no band:** the FEM's Φ₂/Φ at (2, 15) and (2, 0) mm, taken from item 778's and this run's θ at those points, beside step 1's 3.179482 / 0.1714577.

**Control.** Dropping one ball (P/2 in total) misses (I) by ≈ 0.5.

**Windows.** W1 only (item 778's W0 already measured the mesh), `timeout -k 30 600`.

**Tier.** Heavy, `-n 2`, real. **Class:** chain step (solve). **Cost:** ≈ 6.1 / 30.

**Traps, scope and negative result:** as item 778.

**Status it can move.** `IMP-8` step 3: the first gate on Φ₂ at finite t, the gate the narrowed Done-when names. §10 I56 closes, and I61 (the `imp:` example) becomes queueable. **Depends on:** item 778 landing; skip at no cost if it has not.

**Outcome (2026-11-04, 22:30 CST slot, started 04:30 UTC): DONE, a 🟡 record — Φ₂'s first gate at finite t.** New module `tests/validation/test_imp8_rod_tips_fem.py`; it imports item 803's `box` fixture (`_graded_box`), `_need`, `_theta_pts`, `GATE_PTS`, `QP`, `DEGREE` unchanged — the landed module is not edited — and runs under `FEM_EM_IMP8_I55_PHASE=transient`. Width `-n 8` per the 10:30 ruling 6 and the 18:00 ruling 5. Pre-check `pgrep -c python3` read 0. A phase-unset smoke run (`20261104T043212Z_IMP-8-I56-W0-smoke.log`, `-n 2`, 23 s, 4 skipped per rank) checked the import only.
- **W1** (`-n 8`, real; `20261104T043245Z_IMP-8-I56-W1-item796.log`, Status 0, 461 s): cells 485 291 = W0c, t_build 18.4 s; balls 294 / 300 cells (equivalent radii 1.0085 / 0.9963 mm), each ∫q/(P/2) − 1 ≤ 4.4e-16; CN t_call 439.6 s, `n_factorizations` 2.
  - (E) 1.754e-14 ≤ `E_BAND`; (S0) 4.441e-16.
  - **(I) green:** (10, 0) 2.288e-03, (20, 0) 2.267e-03, (10, 15) 5.994e-04, (20, 15) 1.268e-03 against 1e-2. Printed only: (2, 0) 2.491e-03, (2, 15) 2.181e-03 (θ_FEM 30.07493 against the closed form's 30.00947 K/W, the same value item 808 measured).
  - **(c-one) asserted** at (I)'s band, factor 1: the one-ball comparand (z = −L/2 kept) misses θ_FEM by |θ/cmp − 1| = 1.0046 / 1.0045 / 34.37 / 15.48 at the four points, min 1.0045 > 1e-2. The predicted ≈ 0.5, as |cmp/θ − 1|, reads 0.5011 at z = 0 and 0.972 / 0.939 at z = 15 mm — printed, not asserted.
  - **Record (no band), Φ₂/Φ:** (2, 0) mm 0.171766 — this run's θ over item 803's line θ_FEM 17.72315 K/W (`20261104T030134Z_IMP-8-I55-W2-item803.log:145`), FEM over FEM, beside step 1's 0.1714577. (2, 15) mm 3.186418 — item 803 did not evaluate this point, so the ratio is this run's θ over the line's **closed form** 9.438477 K/W, not FEM over FEM, beside step 1's 3.179482. No second solve.
  - Per rank, counted by token: 4 passed on all 8 ranks.
- **Price for 799:** ≈ 440 s of solve, ≈ 461 s of window at `-n 8`, the same as item 803.

### §9 item 799 — 2026-11-03 (§10 I60, chain I — `IMP-8` step 5, the duty-cycled line source in FEM; chain step (solve))

**Executes** §10 I60 on item 778's module (its `_graded_box` and rod source).
- **The schedule:** `source_schedule` (`solve_transient_conduction`, `transient.py:94`) at D = 0.5 with period 450 s: on for [0, 225) s, off for [225, 450) s, and so on.
- **The step:** Δt = 22.5 s, 10 steps per phase. Assert that every switch time lies on the step grid.

**Gates.**
- **(I):** θ_FEM against item 787's Duhamel closed form at item 778's four points, within `INTERIOR_BAND`. Read it at t\* = 900 s (the end of an off phase) and t₂ = 675 s (the end of an on phase).
- (E) and (S0) as item 778.

**Control.** The weekly's: the D = 1 run reproduces item 778's θ at (10, 0) mm to `.6e`.

**Tier.** Heavy, `-n 2`, real, `timeout -k 30 600`. **Class:** chain step (solve). **Cost:** ≈ 6.1 / 30.

**Traps, scope and negative result:** as item 778, with no (St).

**Status it can move.** `IMP-8` step 5 and §10 I60. **Depends on:** items 778 and 787 both landing; skip at no cost if either has not.

**Outcome (2026-11-04, 22:30 CST slot, started 04:30 UTC, take-next after item 796): DONE, a 🟡 record — the duty-cycled line source in FEM.** New module `tests/validation/test_imp8_rod_line_source_duty_fem.py`; it imports item 803's `box` fixture (`_graded_box`), `_source_and_s0` (the rod source), `_theta_pts`, `GATE_PTS`, `DEGREE`, `T_INF`, `W0C_CELLS` unchanged, and item 787's `_duhamel` / `PERIOD` as the comparand; no edit to either landed module, no `src/` edit. Windows by `FEM_EM_IMP8_I60_DUTY` = `0.5` (W1) / `1` (W2). t₂ is the solver's step-30 snapshot (`snapshot_steps`), so one solve gives both read times. Width `-n 8` per the 10:30 ruling 6 and the 18:00 ruling 5. Pre-check `pgrep -c python3` read 0 before each window. A duty-unset smoke run (`20261104T044356Z_IMP-8-I60-W0-smoke.log`, `-n 2`, 22 s, 1 passed / 5 skipped per rank) checked the imports and the grid test.
- **Schedule (asserted):** every switch (0, 225, 450, 675, 900 s) is an integer multiple of Δt = 22.5 s; 10 steps per phase; the solver's 41 (half-)step multipliers are exactly [1, 1] + 9 × 1 + 10 × 0 + 10 × 1 + 10 × 0 at D = 0.5 (sum 21) and all 1 at D = 1.
- **W1** (D = 0.5, `-n 8`, real; `20261104T044429Z_IMP-8-I60-W1-item799.log`, Status 0, 418 s): cells 485 291 = W0c, t_build 17.2 s, 6 495 rod cells; CN t_call 396.9 s, `n_factorizations` 2.
  - (E) 1.435e-14 ≤ `E_BAND`; (S0) 4.441e-16.
  - **(I) green** (θ_FEM against `_duhamel`, K/W, band 1e-2): at t\* = 900 s — (10, 0) 2.240893 vs 2.241722, rel 3.700e-04; (20, 0) 1.637e-03; (10, 15) 1.078e-03; (20, 15) 2.752e-03. At t₂ = 675 s — (10, 0) 2.409751 vs 2.403634, rel 2.545e-03; (20, 0) 2.190e-03; (10, 15) 4.381e-03; (20, 15) 3.896e-03. Max 4.381e-03.
  - Per rank, counted by token: 5 passed, 1 skipped (the D = 1 control) on all 8 ranks.
- **W2, the control** (D = 1, `-n 8`, real; `20261104T045139Z_IMP-8-I60-W2-item799.log`, Status 0, 444 s; CN t_call 421.7 s): θ(10, 0) mm at t\* this run `4.651377e+00`, item 803 `4.651377e+00` (`20261104T030134Z_IMP-8-I55-W2-item803.log:141`) — reproduced to `.6e`, asserted as string equality. The other three t\* values also match item 803's printed digits; (I) at t\* and t₂ max 4.210e-03, (E) 2.021e-14. Per rank: 6 passed on all 8 ranks.
- **Disclosed:** CN runs through the switches with no restart at each edge (the brief asks for none); the t₂ (on-phase end) residuals are the larger ones, ≤ 4.4e-03, against ≤ 2.8e-03 at the off-phase end.

### Review ruling 2026-11-03 10:30 (item 778 BLOCKED — the re-pose as item 800)

**Diagnosis, from the four logs and the parked module (`attempt/IMP-8-I55-20261103T1111Z`, `062a56dc`).** Nothing failed. The item was un-priced. The fixture is ready: W0c measured 485 291 cells and an effective rod radius of 1.0022 mm in 23 s (`20261103T105748Z_IMP-8-I55-W0c.log`). W1 (`20261103T105819Z_IMP-8-I55-W1.log`) then asked one 600 s window at `-n 2` for three solves on that mesh at degree 2: the CN transient, the (St) steady solve and the (c-E) perturbed transient.

`solve_transient_conduction` factors with MUMPS once per (Δt, θ) (`src/fem_em_solver/thermal/transient.py:20`, `:308–313`). `startup_half_steps=2` therefore costs two factorizations, and the (c-E) run repeats them. At ≈ 485 k tets, P2 has ≈ 0.6–0.7 M dofs (review arithmetic: vertices ≈ cells / 5.5, edges ≈ 7 × vertices). The brief's estimate of 150 000–400 000 cells never priced a factorization at that size. W1 never returned from the first transient.

**Ruling.**
1. **Attempt 2 of the step-2 family is a split-window re-pose** (§9 rule (6): "split windows"). The fixture, the degree, the points, the gates, the bands and the `PREDICTED RED` stay as item 778 posed them.
2. **The (St) solve goes first, alone, and is the price.**
3. **(c-E) runs on one full step**, as two backward-Euler half-steps. The balance residual is per (half-)step (`transient.py`'s docstring), so a 1e-3 bookkeeping mismatch shows on the first half-step. The control keeps its meaning.
4. **The transient's width follows the measured price.** `-n 8` is the heavy tier's width, inside the 12-rank ceiling. The rank-local reductions (∫q normalisation, the balance integrals) are seen at `-n 2` in W0c and W1s.
5. **No coarser grading in-slot.** A coarser grading would be a new fixture, and that is the weekly's.
6. Items 796 and 799 inherit item 800's measured width and stop rule. They stay serial on it.

### §9 item 800 — 2026-11-03 (§10 I55, chain I — `IMP-8` step 2, the finite line source in FEM, attempt 2: split windows; chain step (solve))

**Executes** §10 I55 under the ruling above. Item 778's brief binds wherever this one is silent: fixture, source normalisation, points, gates, controls, traps, scope, negative result.

**Code.** Take `tests/validation/test_imp8_rod_line_source_fem.py` from `attempt/IMP-8-I55-20261103T1111Z` by path checkout (code only; append record rows by hand). Keep W0c's grading: `StopAtDistMax = 1` on the inner `Threshold`, the default mesher, no HXT knob. Split the module's solve test into three knob-selected tests (`FEM_EM_IMP8_I55_PHASE` = `steady` / `transient` / `ce`). Each one builds the mesh, prints `t_build`, and prints the solve call's wall time and `n_factorizations` (`TransientConductionResult`, `transient.py:89`).

**Windows.** Before each window, the pre-check `pgrep -c python3` must read 0. Every window is a bare harness call in the foreground, `timeout -k 30 600` container-side, Bash timeout 660000 ms. Window form: `20261103T105819Z_IMP-8-I55-W1.log:12`, with the knob added.

- **W1s** (`-n 2`, real): the build and the (St) solve.
  - Asserts (S0) and (St).
  - Prints the cell count beside W0c's 485 291. A different count stops the item: the grading moved.
  - Prints `t_build` and `t_call` (the (St) call's wall time).
  - `PREDICTED RED: (St) range [0.03, 0.30]`. A red inside this range, with (S0) green, lets the item go on to W2.
  - If W1s is killed at 600 s, re-run it **once** at `-n 8`. If that is killed too, the item is BLOCKED (below).
- **W2** (the transient: (E), (S0), (I), (c-2L)). Projection: T = `t_build` + 3·`t_call`. The transient makes two factorizations plus 42 back-solves, and 3× covers that with margin.
  - Run at `-n 2` if T ≤ 480 s.
  - Otherwise run at `-n 8`. If W1s already ran at `-n 8`, use `-n 8` again.
  - State the projection and the width in the outcome.
- **W3** ((c-E), at W2's width): `n_steps = 1`, `startup_half_steps = 2`, q × (1 + 1e-3), balanced against the nominal P.
  - The miss is ≈ 1e-3 on each half-step, ≫ `E_BAND`.
  - A third test asserts that the nominal one-step run's own balance is ≤ `E_BAND`. That shows the miss is the perturbation's, not the short run's.

PREDICTED RED: (St) range [0.03, 0.30]

**Stop rules.**
- W2 killed at `-n 8` ⇒ BLOCKED with the measured `t_build` / `t_call` / width. The item becomes an `xl` cost probe for the 11-04 weekly to commission. It is never a coarser mesh and never a lower degree in-slot.
- Start no new window after minute 45. Unrun windows stay listed, and the item lands partial.
- Any unpredicted red: stop, as on item 778.

**Tier.** Heavy, `-n 2` (W2 / W3 at `-n 8` by the projection), real. **Class:** chain step (solve). **Cost:** ≈ 6.1 measured / 45 predicted min. Three windows, each ≤ 600 s.

**Traps already paid for.** All of item 778's traps, and:
- A gmsh `Threshold` without `StopAtDistMax` caps the whole box at its SizeMax (W0b, 4 101 108 cells).
- W0 was killed at 300 s inside gmsh. W0c's settings built the mesh in 23 s.
- A window that needs more than ≈ 600 s is unrunnable; never background it.

**Status it can move.** As item 778: `IMP-8` step 2 (🟡 record line, the first FEM gate on Φ), §10 I55, and unblocking items 796 / 799. On a W2 kill: §10 I55 re-priced as an `xl` probe for the weekly. **Depends on:** nothing.

**Outcome 2026-11-03 (12:00 slot, take-next after item 790): BLOCKED, attempt 2. The slot stopped on an unpredicted red: the (c-2L) control in W2. Every gate is green or the predicted red.** The module (`FEM_EM_IMP8_I55_PHASE` = `steady` / `transient` / `ce`, W0c's grading, 485 291 cells asserted) is parked on `attempt/IMP-8-I55-20261103T1835Z` (`4ade96c3`) with both logs.
- **W1s** (`-n 2`, real, 307 s, rc 1; `20261103T182153Z_IMP-8-I55-W1s.log:55–98`):
  - t_build 20.9 s; cells 485 291 = W0c; 6 495 rod cells; effective radius 1.0022 mm.
  - (S0) 1.110e-16, green.
  - **(St), predicted red, lands inside [0.03, 0.30]:** (10, 0) 5.890e-02, (20, 0) 1.020e-01, (10, 15) 7.691e-02, (20, 15) 1.181e-01. The steady box offset is −0.62 K/W at every point, against the review's ≈ −0.9 estimate.
  - t_call for the steady solve: 284.0 s.
  - Per rank: 3 passed, 1 failed (St), 6 skipped.
- **W2 projection and width:** 20.9 + 3 × 284.0 = 873 s > 480 s, so W2 ran at **`-n 8`**.
- **W2** (`-n 8`, real, 415 s, rc 1; `20261103T182723Z_IMP-8-I55-W2-n8.log:109–146`, `:198`):
  - t_build 17.3 s; CN t_call 395.6 s; `n_factorizations` 2.
  - (E) 1.898e-14 ≤ `E_BAND`; (S0) 4.441e-16.
  - **(I) green at all four points:** (10, 0) 1.296e-03, (20, 0) 1.943e-03, (10, 15) 2.867e-03, (20, 15) 3.216e-03, against `INTERIOR_BAND` 1e-2. (2, 0) is printed only, at 6.942e-04.
  - **(c-2L) red:** the L → 2L comparand at (10, 15) mm misses by **7.517e-02**, against the brief's "> 0.1".
  - Per rank: 5 passed, 1 failed (c-2L), 4 skipped, on all 8 ranks.
- **Not run:** W3 ((c-E) and the nominal one-step balance). The slot stops on an unpredicted red.
- **Reading, for the review (not acted on).** (c-2L) was posed in item 778 as a *predicted* factor ("Controls (predicted)") but asserted in the module. Rule (e) says a predicted factor is printed, never asserted. Its 0.1 threshold was never measured on this comparison. θ_FEM matches the L comparand to 2.9e-03 at that point, so the 7.5e-02 miss is the closed-form gap between L and 2L at (10, 15) mm, t\* = 900 s. With √(κt\*) ≈ 11.4 mm, the extra length at |z| > L/2 = 15 mm is barely felt, and halving q′ is mostly compensated. The control still separates: 7.5e-02 is ≈ 26× the (I) residual there. **Unblock:** the review re-registers (c-2L)'s threshold from the closed form's own L-vs-2L ratio at that point (rule (e): asserted only on a measured comparison, and this log is that measurement), or moves the control to a point where it is larger. Then the next attempt re-runs W2 at `-n 8` (≈ 7 min) and W3. W1s's (St) record stands, unless the module edit forces a re-run under rule (i).
- **Price for items 796 / 799 (ruling 6):** at `-n 8`, a 40-step CN transient on this mesh at degree 2 costs ≈ 396 s of solve and ≈ 415 s of window.

### Review ruling 2026-11-03 18:00 (item 800 BLOCKED — (c-2L) re-registered under rule (e); the re-pose as item 803)

**Diagnosis, from the two logs and the parked module (`attempt/IMP-8-I55-20261103T1835Z`, `4ade96c3`).** Nothing in the solve failed. Every gate W1s and W2 evaluated is green or the predicted red: (S0) 1.110e-16, (St) inside [0.03, 0.30] (`20261103T182153Z_IMP-8-I55-W1s.log:55–98`); (E) 1.898e-14, (I) ≤ 3.216e-03 at all four points (`20261103T182723Z_IMP-8-I55-W2-n8.log:109–146`). The one red is the (c-2L) **control**, 7.517e-02 at (10, 15) mm against "> 0.1". Item 778 wrote that threshold without a measurement (its brief: "Controls."; no log line backs 0.1), so under implementer-run.md rule (e) it is a *predicted* factor and should have been printed, not asserted. The slot was right to stop: the label was missing, and rule (e) says an executor meeting an unlabelled failing pre-registered assertion stops.

**Ruling.**
1. **(c-2L) is re-registered as an asserted control at the gate's own band:** the L → 2L comparand misses θ_FEM at (10, 15) mm by **> `INTERIOR_BAND`** (1e-2). That is what a control on (I) must show — the gate at its band would reject the wrong length — and the W2 log is the prior measurement on the same comparison and fixture that rule (e) asks for (7.517e-02, `20261103T182723Z_IMP-8-I55-W2-n8.log`, the (c-2L) line). The old 0.1 is printed beside the measured miss as "predicted 0.1 (item 778), not asserted". No band moves: `INTERIOR_BAND` is (I)'s own, imported.
2. **W3's third test is dropped as redundant, not loosened.** The nominal one-step balance it asked for is the first two half-steps of W2's own run at the same Δt and `startup_half_steps`, and W2's (E) is the maximum over every half-step (1.898e-14 ≤ `E_BAND`). W3 therefore runs the (c-E) perturbed solve alone, and its comparison cites W2's (E) line as the nominal reference.
3. **W1s's (St) record stands without a re-run** (rule (i): "or the row says it was not"). The module edit is confined to the (c-2L) assertion inside the `transient` phase and the `ce` phase; the `steady` phase's code is untouched. The outcome says so, and quotes `git diff` of the module against `4ade96c3` restricted to the steady test (empty).
4. **Family count:** attempt 3 of 4 in the step-2 family. A further red that is not a control-label question freezes the family at the next attempt (§9 rule (2)).
5. Items 796 and 799 depend on item 803 now (read "item 800" in their briefs as "item 803"), at `-n 8` and ≈ 415 s per 40-step window (ruling 6 above).

### §9 item 803 — 2026-11-03 (§10 I55, chain I — `IMP-8` step 2, the finite line source in FEM, attempt 3: (c-2L) re-registered; chain step (solve))

**Executes** §10 I55 under the ruling above. Item 800's brief binds wherever this one is silent (fixture, phases, windows, stop rules, traps); item 778's beneath it.

**Code.** Take `tests/validation/test_imp8_rod_line_source_fem.py` from `attempt/IMP-8-I55-20261103T1835Z` (`4ade96c3`) by path checkout — code only. Copy its two logs by path as the records of attempt 2 and **append their two test-results rows by hand**; `git diff --stat` on `test-results.md` shows insertions only. Edits, and only these:
- In the `transient` phase, (c-2L): assert the L → 2L miss at (10, 15) mm **> `INTERIOR_BAND`**, and print beside it "predicted 0.1 (item 778), not asserted". Label it *asserted* in the docstring, citing `20261103T182723Z_IMP-8-I55-W2-n8.log`.
- In the `ce` phase: run only the perturbed one-step solve (q × (1 + 1e-3), `n_steps = 1`, `startup_half_steps = 2`), asserting its balance against the nominal P misses by > `E_BAND` on the first half-step — factor 1, the gate's own band, the weakest separation a control can assert — with the predicted ≈ 1e-3 printed beside it, labelled *predicted* (rule (e): no measurement of this comparison on this fixture exists yet; this window is it). Remove the nominal one-step test; the docstring cites W2's (E) line as the nominal reference.

**Windows.** Pre-check `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 before each. Bare, foreground, `timeout -k 30 600`, Bash timeout 660000 ms. Window form `20261103T182723Z_IMP-8-I55-W2-n8.log:12`, with the phase knob.
- **W2** (`FEM_EM_IMP8_I55_PHASE=transient`, `-n 8`, real): (E), (S0), (I), (c-2L). Measured ≈ 415 s.
- **W3** (`FEM_EM_IMP8_I55_PHASE=ce`, `-n 8`, real): (c-E). Two factorizations and two back-solves: ≈ 150–330 s (W2's 395.6 s CN call made two factorizations and 42 back-solves).
- No W1s: its record stands (ruling 3). State `git diff 4ade96c3 -- <module>` restricted to the steady test in the outcome.

**Gates.** All of item 778's, imported and unmoved: (E) `E_BAND`, (S0) `B_BAND`, (I) `INTERIOR_BAND` at the four points, (St) as W1s's record. Nothing is pre-registered red in W2 / W3: **no `PREDICTED RED:` line.** (St)'s predicted red was read in attempt 2 and is a record.

**Controls.** (c-2L) asserted > `INTERIOR_BAND` (measured 7.517e-02 at attempt 2); (c-E) asserted > `E_BAND` (factor 1; the ceiling is the 1e-3 perturbation itself, ≈ 1e6 × `E_BAND`, so it is reachable by construction), the predicted ≈ 1e-3 printed.

**Tier.** Heavy, `-n 8` (by attempt 2's measured projection), real. **Class:** chain step (solve). **Cost:** ≈ 6.1 measured / 25 predicted min (two windows, ≈ 10–13 min of compute).

**Traps already paid for.** Item 800's, and: a predicted factor is printed, never asserted (rule (e), this ruling); path checkouts take code only; count pytest verdicts per rank by token (8 ranks here).

**Scope.** As item 778. **Not claimed:** a rate, I56, I60, entry 48's rod.

**Negative result.** Any red: stop, park on `attempt/IMP-8-I55-*`, BLOCKED. A W2 or W3 kill at `-n 8` ⇒ an `xl` cost probe for the weekly (item 800's stop rule).

**Status it can move.** `IMP-8` step 2 (a 🟡 record line on its ✅ row, the first FEM gate on Φ), §10 I55 DONE, and items 796 / 799 unblock. **Depends on:** nothing.

**Outcome (2026-11-04, 21:00 slot, started 03:00 UTC): DONE, attempt 3, a 🟡 record — the first FEM gate on Φ.** Module `tests/validation/test_imp8_rod_line_source_fem.py` taken by path from `4ade96c3`, with attempt 2's two logs and their two test-results rows appended by hand (insertions only). Edits only as the ruling lists: (c-2L) asserted > `INTERIOR_BAND` with "predicted 0.1 (item 778), not asserted" printed; the `ce` phase runs the perturbed one-step solve alone, asserted > `E_BAND` on the first half-step with the ≈ 1e-3 printed as predicted; the nominal one-step test removed, its docstring citing W2's (E) line. `git diff 4ade96c3 -- tests/validation/test_imp8_rod_line_source_fem.py` has no hunk in the steady block (fixture `steady_run` and its two tests, lines 256–299 of the parked module), so W1s's (St) record stands without a re-run (ruling 3). No `src/` edit. Pre-check `pgrep -c python3` read 0 before each window.
- **W2** (`FEM_EM_IMP8_I55_PHASE=transient`, `-n 8`, real; `20261104T030134Z_IMP-8-I55-W2-item803.log`, Status 0, 461 s): cells 485 291 = W0c, t_build 17.1 s, 6 495 rod cells, effective radius 1.0022 mm; CN t_call 439.8 s, `n_factorizations` 2.
  - (E) 2.098e-14 ≤ `E_BAND` (attempt 2 read 1.898e-14 on the same command; both are round-off, disclosed); (S0) 4.441e-16.
  - **(I) green:** (10, 0) 1.296e-03, (20, 0) 1.943e-03, (10, 15) 2.867e-03, (20, 15) 3.216e-03 against 1e-2 — the same digits as attempt 2. (2, 0) printed only, 6.942e-04.
  - **(c-2L) asserted:** 7.5169e-02 > `INTERIOR_BAND`.
  - Per rank, counted by token: 6 passed, 3 skipped on all 8 ranks.
- **W3** (`FEM_EM_IMP8_I55_PHASE=ce`, `-n 8`, real; `20261104T030930Z_IMP-8-I55-W3-item803.log`, Status 0, 210 s): t_build 17.7 s; the perturbed one-step solve t_call 189.0 s, `n_factorizations` **1** (both startup half-steps share one backward-Euler operator at Δt/2; the brief's "two factorizations" was a prediction), 2 residuals.
  - **(c-E) asserted:** own first-half-step residual 1.264e-14; ∫Qe/P − 1 = 1.000000e-03; first-half-step miss against the nominal P in [1.000000e-03, 1.000000e-03] > `E_BAND` 1e-9 (predicted ≈ 1e-3 — met). (S0) 4.441e-16.
  - Per rank: 3 passed, 6 skipped on all 8 ranks.
- **Price for 796 / 799:** a 40-step CN transient at `-n 8` costs ≈ 440 s of solve, ≈ 461 s of window — up from attempt 2's 396 / 416 s; still inside a 600 s window.

### §9 item 808 — 2026-11-03 (§10 I62, chain I — daily, 2026-11-03 18:00, provisional under rule (6); `IMP-8` step 7, the two-tip source under a duty cycle, closed form; lane 3 → chain step)

*Written by the 2026-11-03 18:00 daily review under rule (6): chain I's next step toward Phase 7's stated subgoal (the temperature around an implant in gel, `IMP-8`'s closed forms and their FEM twins). It is I59's twin on step 1's Φ₂, so that I56's tip source has the duty-cycled comparand I60 has for the line. Provisional until the 11-04 weekly ratifies or reverts it. Chain I provisional count since the 11-02 weekly: 1 of 6. Family I62, attempt 1 of 4.*

**Executes** I59 (item 787) with the source shape swapped. New `tests/validation/test_imp8_rod_tips_duty_closed_form.py`, numpy / scipy. Import item 787's module (`test_imp8_rod_line_source_duty_closed_form.py`) for its schedule, Duhamel step sum and impulse-kernel convolution, and step 1's `_phi_tips` (`test_imp8_rod_line_source_closed_form.py:89`) as the kernel. Never edit either landed module; if item 787's helpers are bound to `_phi`, wrap them with the kernel as an argument in the new module.

**Anchors (imported, unmoved).** At D = 1, 0.5, 0.25, period 0.5 t\*, at (2, 15) and (10, 0) mm, exactly item 787's grid:
- **(D1):** at D = 1 the duty sum equals the continuous Φ₂ form exactly (item 787's (D1) construction).
- **(Duh):** the step sum against the impulse-kernel convolution of Φ₂ at `MAT-17`'s bar (1e-10 of θ_ss, as item 787).
- **(Q-d) / (Win-d):** `cem43` against the quadrature at `DOSE_BAR`, the six-minute window at `AVG_BAR`, as item 787.

**Control.** The continuous Φ₂ form misses the D = 0.5 history at t\* (an off-phase end) at both points: **asserted > `INTERIOR_BAND`** (1e-2 relative, factor 1 — rule (e): item 787's measured ≥ 1.072 is on Φ, not on Φ₂, so it backs no factor here), with the *predicted* 0.1 printed beside the measured miss.

**Printed (🟡 record, no band).** θ(t\*) per watt at both points and three duties, beside item 787's line-source values (9.438477 / 4.645356, 2.096145 / 2.241722, 0.8986021 / 1.00437 K/W), and the ratio tips / line per duty.

**Tier.** Smoke, `-n 1`, real, `timeout -k 30 60` (item 787: 3 s). **Class:** lane 3 closed form → chain step. **Cost:** ≈ 2.1 measured / 15 predicted min.

**Traps already paid for.** Item 787's: every switch time on the step grid; `MAT-17`'s bars by import; `_phi_tips` takes (r, z, t) with the balls at z = ±L/2 — at (2, 15) mm the point is 1 mm from a ball's surface, so print (a / R₊)² beside it.

**Scope.** Closed form only. **Not claimed:** any FEM (I56 / I60 are those), entry 48's rod.

**Negative result.** A red: stop, BLOCKED; no `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 7 (a 🟡 record line on its ✅ row) and §10 I62 (provisional) DONE; with I56 landed it is the comparand an FEM tip-duty twin would need — the weekly's to enumerate. **Depends on:** nothing.

**Outcome (2026-11-04, 21:00 CST slot, take-next at minute 29): BLOCKED, attempt 1. An unpredicted red on an inherited precondition, not on a gate band.** New `tests/validation/test_imp8_rod_tips_duty_closed_form.py` (item 787's `_edges`, grid, duties, points and `SCALE` imported; its step sum and convolution wrapped with the kernel as an argument; ∂Φ₂/∂t in closed form, (L/2)·Σ± e^{−R±²/4κt}/(2√(πκ) t^{3/2})), parked with its log on `attempt/IMP-8-I62-20261104T0331Z` (`735ecb96`). `20261104T033016Z_IMP-8-item808-W1.log`: `-n 1`, real, Status 1, 3 s; 5 passed, 1 failed, 3 skipped (complex-only environment tests).
- **Green:** (D1) 0.000e+00 / 3.018e-17 (band 1e-13); (Duh) worst 2.951e-16 (band 1e-10); control (c-cont), asserted > `INTERIOR_BAND`: 1.1127e+01 at (2, 15) mm and 6.9531e-01 at (10, 0) mm, beside the predicted 0.1. Printed: closed ∂Φ₂/∂t against a central difference, worst rel 4.424e-05 (h = 1e-3 u, printed only); (a/R₊)² at (2, 15) mm 2.500e-01.
- **Red:** `test_q_d_win_d_and_printed_record` stops at its first case on item 787's `max T < 43 °C` assertion, which the brief inherits ("as item 787"): at D = 1, (2, 15) mm, θ_tips(t\*) = 30.00946732 K/W, so T = 20 °C + θ(1 W) = 50.009 °C. Item 787's own docstring gives the reason for that assertion: the quadrature reference (`_q_check_g`) integrates R_below only. So (Q-d) / (Win-d) and the printed record never ran. θ(t\*) per watt from the (Duh) lines: (2, 15) / (10, 0) mm 30.00947 / 1.927920 (D = 1), 2.474674 / 1.137210 (D = 0.5), 0.9807666 / 0.5680058 (D = 0.25) K/W.
- **Not done in-slot:** P is not rescaled, the 43 °C assertion is not dropped, and no other point is chosen. Each of those would change the item's target or its inherited gate. (2, 15) mm sits 1 mm from the top ball's surface, where 1 W on a 1 mm ball is ≈ 3× the line's 9.44 K/W.
- **Unblock (the review's):** pose the dose's temperature scale for Φ₂: a power such that max T < 43 °C at every point (for example P scaled from the measured 30.009 K/W peak), or a `_q_check_g` path that covers R_above. Then re-run the one window (≈ 3 s).

### Review ruling 2026-11-04 03:00 (item 808 BLOCKED — the dose's temperature scale on Φ₂; I62 re-posed as item 811)

**Reading, from `20261104T033016Z_IMP-8-item808-W1.log` (on `attempt/IMP-8-I62-20261104T0331Z`, `735ecb96`).** Every anchor that ran is green: (D1) ≤ 3.018e-17, (Duh) 2.951e-16, and the control (c-cont) at 1.1127e+01 / 6.9531e-01 against `INTERIOR_BAND`. The red is item 787's inherited precondition `max T < 43 °C`, at D = 1 and (2, 15) mm. There θ_tips(t\*) = 30.00947 K/W, so at item 787's 1 W, T = 50.009 °C. Item 787's docstring says why the precondition exists: `_q_check_g`'s quadrature integrates R_below only. Nothing is wrong with the closed form. The 1 W scale suits the line (peak 9.44 K/W, max T 29.44 °C) but not two 1 mm balls carrying the same watt, which run ≈ 3 × hotter 1 mm from a ball's surface. The slot was right not to rescale in-slot.

**Ruling.**
1. **The power that sets the dose's temperature scale is part of the fixture, not the target.** Item 787's brief makes it a convention ("the baseline and per-watt scaling that `test_mat21_box_dose_closed_form.py` uses … print it"). The targets are (D1), (Duh), (Q-d) and (Win-d) at item 787's bars. They are relative bars, and they do not change with P. The precondition that keeps `_q_check_g` valid stays, unmoved.
2. **The dose for Φ₂ is evaluated at P_DOSE = 0.5 W:** T = T_INF + 0.5 · θ(1 W). The measured peak is 30.00947 K/W at D = 1, (2, 15) mm, t\*, so max T = 35.005 °C < 43 °C at every duty and point. That margin is the reason for 0.5 W rather than a value tuned to just under 43 °C. θ per watt and the tips / line ratio are printed per watt, unchanged. **The doses are not comparable in absolute terms with item 787's 1 W doses,** and the print says which scale each one is on.
3. **Not chosen: a `_q_check_g` path that covers R_above.** That would be a change to an imported reference, a larger item, and nothing in chain I needs a dose above 43 °C.
4. **Family count:** I62 attempt 2 of 4, rule (6) provisional (chain I, 1 of 6, unchanged; a re-pose adds no step).

### §9 item 811 — 2026-11-04 (§10 I62, chain I — daily, provisional; `IMP-8` step 7, the two-tip source under a duty cycle, closed form, attempt 2: the dose at 0.5 W; lane 3 → chain step)

**Executes** §10 I62 under the ruling above. Item 808's brief binds wherever this one is silent: anchors, control, printed record, traps, scope.

**Code.** Take `tests/validation/test_imp8_rod_tips_duty_closed_form.py` from `attempt/IMP-8-I62-20261104T0331Z` (`735ecb96`) by path checkout, code only. Copy `20261104T033016Z_IMP-8-item808-W1.log` by path as attempt 1's record, and **append its test-results row by hand** (insertions only). Make these edits and no others:
- Add `P_DOSE_W = 0.5` with a comment that cites this ruling.
- In `test_q_d_win_d_and_printed_record`, set T = `T_INF_0` + `P_DOSE_W` · θ. Keep the `max T < 43 °C` assertion unchanged.
- Print the scale on the (Q-d) line as "dose at P = 0.5 W (item 787's line: 1 W)".
- Leave θ per watt, the six-minute window per watt and the tips / line ratio as they are.
- The (D1), (Duh) and control tests are untouched.
- Update the module docstring's "Temperature scale" sentence to say the same.

**Window.** W1 `IMP-8-item811-W1`: `-n 1`, real, `timeout -k 30 60`, bare and foreground. Window form: `20261104T033016Z_IMP-8-item808-W1.log:12`. Measured: 3 s.

**Anchors (imported, unmoved).** (D1) at `D1_BAND`; (Duh) at `DUH_BAND`; (Q-d) at `DOSE_BAR`; (Win-d) at `AVG_BAR`. Precondition: max T < 43 °C, asserted. **No `PREDICTED RED:` line.**

**Control (asserted).** The continuous Φ₂ form misses the D = 0.5 history at t\* by > `INTERIOR_BAND` at both points. Measured: 1.1127e+01 and 6.9531e-01.

**Printed (🟡 record, no band).**
- θ(t\*) per watt at both points and three duties, beside item 787's line values.
- The tips / line ratio per duty.
- The six-minute maximum per watt.
- CEM43 at 0.5 W.
- The predicted max T, 35.005 °C, beside the measured one.

**Tier.** Smoke, `-n 1`, real. **Class:** lane 3 closed form → chain step. **Cost:** ≈ 4.1 measured / 10 predicted min.

**Traps already paid for.** Item 808's, and:
- The 43 °C precondition guards `_q_check_g`'s R_below-only quadrature. It is never dropped.
- Path checkouts take code only.

**Scope.** Closed form only. **Not claimed:** any FEM (that is I63), a recommended duty, or a dose at any other power.

**Negative result.** A red: stop, BLOCKED.

**Status it can move.**
- `IMP-8` gains step 7, a 🟡 record line on its ✅ row.
- §10 I62 (provisional) is DONE.
- Item 814 (I63) gets its comparand.

**Depends on:** nothing.

**Outcome (2026-11-04, 04:30 CST slot): DONE, attempt 2 — a 🟡 record line on `IMP-8`'s ✅ row (step 7).** Attempt 1's module taken from `attempt/IMP-8-I62-20261104T0331Z` (`735ecb96`) by path, code only, with the ruling's edits and no others: `P_DOSE_W = 0.5` (with `MAX_T_PREDICTED = 35.005`, printed only), T = `T_INF_0` + `P_DOSE_W` · θ, the scale printed on the (Q-d) line, the docstring's temperature-scale sentence; (D1), (Duh) and the control untouched. Attempt 1's log copied by path with its test-results row appended by hand. `20261104T103134Z_IMP-8-item811-W1.log`: `-n 1`, real, Status 0, 6 s; 6 passed, 3 skipped (complex-only environment tests).
- **Asserted, green:** (D1) 0.000e+00 / 3.018e-17 (band 1e-13); (Duh) worst 2.951e-16 (band 1e-10); (Q-d) worst rel 2.052e-16 (bar 1e-10); (Win-d) worst rel 3.987e-16 (bar 1e-12), maxima inside their brackets; precondition max T < 43 °C: measured 35.004734 °C (D = 1, (2, 15) mm) beside the predicted 35.005 °C; control (c-cont) 1.1127e+01 / 6.9531e-01 > `INTERIOR_BAND` (predicted 0.1, not asserted).
- **Printed (🟡 record):** θ_tips(t\*) per watt at (2, 15) / (10, 0) mm 30.00947 / 1.92792 (D = 1), 2.474674 / 1.13721 (D = 0.5), 0.9807666 / 0.5680058 (D = 0.25) K/W; tips / line 3.179482 / 0.4150209, 1.180583 / 0.5072929, 1.091436 / 0.5655344; six-minute maximum per watt 29.54154 / 1.531514, 17.81965 / 0.8924947, 9.05116 / 0.4840151 K/W; CEM43 at 0.5 W 5.531217e-03 / 2.533198e-11, 2.971548e-04 / 1.928274e-11, 1.651446e-05 / 1.598229e-11 s (not comparable in absolute terms with item 787's 1 W line doses).

### §9 item 814 — 2026-11-04 (§10 I63, chain I — daily, 2026-11-04 03:00, provisional under rule (6); `IMP-8` step 8, the two-tip source under a duty cycle in FEM; chain step (solve))

*Written by the 2026-11-04 03:00 daily review under rule (6). It is chain I's next step toward Phase 7's stated subgoal, the temperature around an implant in gel (`IMP-8`'s closed forms and their FEM twins). It is I60's construction with I56's source and I62's comparand, so the duty-cycled tip source gets the same FEM gate that the line source has. Provisional until the 11-04 weekly ratifies or reverts it. Chain I's provisional count since the 11-02 weekly is now 2 of 6 (I62, I63). Family I63, attempt 1 of 4.*

**Executes** §10 I63. New module `tests/validation/test_imp8_rod_tips_duty_fem.py`, real build. Import by name, and never edit a landed module:
- item 803's `box` fixture (`_graded_box`), `_theta_pts`, `GATE_PTS`, `DEGREE`, `T_INF` and `W0C_CELLS` (`test_imp8_rod_line_source_fem.py`), the way item 799's module imports them;
- item 796's `_tips_source` (`test_imp8_rod_tips_fem.py:74`), which puts P/2 on each of the two balls;
- item 799's `_schedule`, `_switches`, `DT`, `N_STEPS` and `STEP_T2` (`test_imp8_rod_line_source_duty_fem.py:54–79`);
- item 811's `_duhamel(kernel, r, z, t, duty)` with `_phi_tips` as the kernel, as the comparand (`test_imp8_rod_tips_duty_closed_form.py:73`), per watt.

**Run.** One CN solve at D = 0.5, period 450 s, Δt = 22.5 s, 40 steps to t\* = 900 s, `snapshot_steps=(STEP_T2,)` for t₂ = 675 s. This is item 799's `duty_run` body with the rod source swapped for `_tips_source`.

**Gates (imported, unmoved).**
- **(I):** θ_FEM against item 811's Φ₂ Duhamel sum, at item 803's four gate points ((10, 0), (20, 0), (10, 15) and (20, 15) mm), at t\* and at t₂, within `INTERIOR_BAND`.
- **(E):** ≤ `E_BAND`.
- **(S0):** each ball's ∫q / (P/2) − 1 at `B_BAND`, as item 796.
- The switches lie on the Δt grid, and the solver's multipliers are as item 799 asserted.

**Control (asserted).** The continuous Φ₂ comparand (D = 1, `_phi_tips` step response) misses θ_FEM at (10, 0) mm at t\* by > `INTERIOR_BAND`, factor 1. Item 811's closed-form control measured 6.9531e-01 for this comparison at this point, and θ_FEM there is within (I) of the D = 0.5 sum. **Printed only:** the same miss at the other three points.

**Printed (records, no band).**
- θ at (2, 15) and (2, 0) mm at t\* and t₂, beside item 811's closed form at those points. 1 mm from a ball the closed form is not the comparand: print (a / R₊)².
- The tips / line ratio, FEM over FEM, at (10, 0) mm, over item 799's θ_FEM at t\*: `2.240893` K/W (`20261104T044429Z_IMP-8-I60-W1-item799.log`).
- t_build, t_call and `n_factorizations`.

**Windows.**
- Pre-check: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
- W0 smoke: `-n 2`, `timeout -k 30 120`. Imports and the grid test only, with no solve, as item 799's W0.
- W1: `-n 8`, real, bare and foreground, `timeout -k 30 600`, Bash timeout 660000 ms. Window form: `20261104T044429Z_IMP-8-I60-W1-item799.log:12`.
- **Price:** item 799's W1 took 418 s and item 796's 461 s on the same mesh, steps and width, so ≈ 420–470 s.

**Tier.** Heavy, `-n 8` (items 796 / 799's measured width), real. **Class:** chain step (solve). **Cost:** ≈ 11.2 measured / 30 predicted min.

**Traps already paid for.** Items 778–803, 796 and 799:
- Rank-local integrals are all-reduced before any assertion.
- Count verdicts per rank, by token, across 8 ranks.
- One fixed batch per `.eval`.
- CN runs through the switches with no restart. Item 799 measured the on-phase-end residuals as the larger ones; that is a reading, not a defect.
- A predicted factor is printed, never asserted (rule (e)).
- Pytest fixtures cannot be imported across modules except by the `noqa: F401` pattern item 799 uses for `box`.

**Scope.** One mesh, one Δt, D = 0.5. **Not claimed:** a rate, entry 48's rod, or the dose in FEM.

**Negative result.** Any red: stop, park on `attempt/IMP-8-I63-*`, BLOCKED. No `PREDICTED RED:` line. A W1 kill at 600 s parks the item with its measured times; never raise the timeout.

**Status it can move.**
- `IMP-8` gains step 8, a 🟡 record line on its ✅ row: the duty-cycled tip source's first FEM gate.
- §10 I63 (provisional) is DONE.

**Depends on:** item 811 landing, since its module is the comparand. If 811 is BLOCKED, skip this item at no cost.

**Outcome (2026-11-04, 04:30 CST slot, take-next): DONE, attempt 1 — a 🟡 record line on `IMP-8`'s ✅ row (step 8).** New module `tests/validation/test_imp8_rod_tips_duty_fem.py`, every named object imported, no landed module edited. W0 `20261104T105756Z_IMP-8-I63-W0-item814.log`: `-n 2`, Status 0, 38 s, 1 passed / 5 skipped per rank (grid + comparand precondition; no solve). W1 `20261104T105843Z_IMP-8-I63-W1-item814.log`: `-n 8`, real, Status 0, 544 s (t_call 511.3 s, `n_factorizations` 2, t_build 0.4 s); "6 passed" on all 8 ranks, no FAILED token. Above the brief's 420–470 s price, inside the 600 s ceiling.
- **Asserted, green:** (I) at (10, 0) / (20, 0) / (10, 15) / (20, 15) mm, t\* = 900 s: 2.105e-03 / 1.558e-03 / 1.829e-03 / 1.123e-03; t₂ = 675 s: 2.169e-03 / 3.477e-03 / 2.216e-03 / 1.308e-03 (band 1e-02); (E) 1.678e-14 (band 1e-09); (S0) 4.441e-16 total, per ball 2.220e-16 / 4.441e-16 (294 / 300 cells); multipliers 41, item 799's pattern; control (c-cont) at (10, 0) mm, t\*: |θ_FEM / cont − 1| = 4.0889e-01 > `INTERIOR_BAND`. **Disclosed:** item 811's 6.9531e-01 is the same comparison in the other direction (|cont / D=0.5 sum − 1| = 1.927920 / 1.137210 − 1); this module uses (I)'s convention, |θ_FEM / comparand − 1|, as item 796 did.
- **Printed (records, no band):** θ_FEM at (2, 15) mm 2.286806 (t\*) / 27.84130 (t₂) K/W beside the closed form 2.474674 / 27.53479 (rel 7.592e-02 / 1.113e-02, (a/R₊)² 0.25 — not the comparand); at (2, 0) mm 1.693052 / 1.350644 beside 1.689945 / 1.346719 (rel 1.839e-03 / 2.914e-03); tips / line at (10, 0) mm, t\*, FEM / FEM 1.139603 / 2.240893 = 0.508549 (item 811's closed form 0.5072929); (c-cont) misses at the other three points 3.8725e-01 / 5.6978e-01 / 3.9383e-01.

## §7 row moved verbatim 2026-11-04 (rotate_plan_archive.py chunks)

 The heated rod in gel: the finite line source in closed form (Phase 7; opened 2026-11-02 10:30 daily review, rule (5) lane 3; the weekly ratifies) | ✅ on the closed form, 2026-11-02, **ratified by the 11-02 weekly with its (Rate) band [1.95, 2.05]** (slope 2.0047): steps 0 and 1 green (§9 items 759, 760; `20261102T193249Z_IMP-8-item760-W1.log`, 12 passed, 3 s), controls (c-κ) 0.954, (c-fin) 9.5e-4, (c-L) 0.255, (c-one) 0.499, ratios 3.179482 / 0.1714577 printed. The FEM half is posed, each a 🟡 record line: step 2 = I55 (the line source in FEM, heavy, real; **🟡 record 2026-11-04**, §9 item 803, `-n 8`, `20261104T030134Z_IMP-8-I55-W2-item803.log` / `20261104T030930Z_IMP-8-I55-W3-item803.log`: on a graded gmsh box (485 291 cells, degree 2, CN N = 40) (I) θ_FEM against q′/(4πk)·Φ at t\* = 900 s ≤ 3.216e-03 at four points, (E) 2.098e-14, (S0) 4.441e-16; controls (c-2L) 7.517e-02 > `INTERIOR_BAND`, (c-E) 1.000e-03 > `E_BAND`; (St) predicted red inside [0.03, 0.30], 5.9e-02 … 1.18e-01, a record from `20261103T182153Z_IMP-8-I55-W1s.log`), step 3 = I56 (the two-tip twin, Φ₂'s first gate at finite t; **🟡 record 2026-11-04**, §9 item 796, `-n 8`, `20261104T043245Z_IMP-8-I56-W1-item796.log`, 461 s: on item 803's box (485 291 cells, degree 2, CN N = 40), two DG0 balls of P/2 each (294 / 300 cells) — (I) θ_FEM against q′/(4πk)·Φ₂ at t\* = 900 s ≤ 2.288e-03 at four points, (E) 1.754e-14, (S0) 4.441e-16; control (c-one) min 1.0045 > `INTERIOR_BAND` (one ball short by 0.501 at z = 0, predicted ≈ 0.5); Φ₂/Φ printed 0.171766 at (2, 0) mm (FEM / FEM) and 3.186418 at (2, 15) mm (FEM over the line's closed form) beside 0.1714577 / 3.179482), step 4 = I59 (the duty cycle, closed form; **🟡 record 2026-11-03**, §9 item 787, `20261103T211919Z_IMP-8-item787-W1.log`: P = 450 s, θ(t\*) per watt at (2, 15) / (10, 0) mm 9.438477 / 4.645356 (D = 1), 2.096145 / 2.241722 (D = 0.5), 0.8986021 / 1.00437 (D = 0.25) K/W; (Duh) 2.95e-16), step 5 = I60 (I59's FEM twin; **🟡 record 2026-11-04**, §9 item 799, `-n 8`, `20261104T044429Z_IMP-8-I60-W1-item799.log`, 418 s: on item 803's box and rod source, `source_schedule` D = 0.5, P = 450 s, Δt = 22.5 s (switches on the grid, asserted) — (I) θ_FEM against item 787's Duhamel sum ≤ 4.381e-03 at four points at t\* = 900 s (off-phase end) and t₂ = 675 s (on-phase end), (E) 1.435e-14, (S0) 4.441e-16; control `20261104T045139Z_IMP-8-I60-W2-item799.log`, 444 s: the D = 1 run reproduces item 803's θ(10, 0) mm 4.651377e+00 to `.6e`), step 7 = I62 (I59 on Φ₂, closed form; **🟡 record 2026-11-04**, §9 item 811, `20261104T103134Z_IMP-8-item811-W1.log`, 6 s: (D1) ≤ 3.018e-17, (Duh) 2.951e-16, (Q-d) ≤ 2.052e-16, (Win-d) ≤ 3.987e-16 at item 787's bars; dose at P = 0.5 W (review ruling 2026-11-04 03:00), max T 35.004734 °C < 43; θ_tips(t\*) per watt at (2, 15) / (10, 0) mm 30.00947 / 1.92792 (D = 1), 2.474674 / 1.13721 (D = 0.5), 0.9807666 / 0.5680058 (D = 0.25) K/W, tips / line 3.179482 / 0.4150209, 1.180583 / 0.5072929, 1.091436 / 0.5655344; control (c-cont) 1.1127e+01 / 6.9531e-01 > `INTERIOR_BAND`), step 8 = I63 (I62's FEM twin; **🟡 record 2026-11-04**, §9 item 814, `-n 8`, `20261104T105843Z_IMP-8-I63-W1-item814.log`, 544 s: on item 803's box with item 796's two balls (294 / 300 cells) under item 799's schedule (D = 0.5, P = 450 s, Δt = 22.5 s) — (I) θ_FEM against item 811's Φ₂ Duhamel sum ≤ 3.477e-03 at four points at t\* and t₂, (E) 1.678e-14, (S0) 4.441e-16; control (c-cont) at (10, 0) mm 4.0889e-01 > `INTERIOR_BAND`; tips / line at (10, 0) mm, t\*, FEM / FEM 0.508549), then I61 (an `imp:` example) and I58 (the lumped reading of entry 48's rod; **🟡 record 2026-11-03**, §9 item 786, `20261103T211525Z_IMP-8-item786-W1.log`: at the tip 1 g ball θ_line 7.307753 / θ_tips 12.93694 K/W at t\* = 900 s, ratio 1.770303, at P_rod 6.939109e-07 / 1.228433e-06 K; (Bq) n vs 2n ≤ 8.1e-10). *History and brief: `docs/planning/chunks/IMP-8.md`.* **Done-when:** steps 0 and 1 green: the steady, infinite-line and far-field limits, the energy identities and the midpoint rate in [1.95, 2.05], with a control that fires on each gate that has content at finite t (narrowed by the 11-02 weekly: Φ₂ at finite t and (Sym) have no closed-form gate; Φ₂'s first gate is the FEM, I56). 

### Review ruling 2026-11-05 03:00 — chain I's I64 / I65 as execution detail

Writing the briefs for the 11-04 weekly's I64 and I65, this review fixed the values the step text leaves open and corrected two premises. Each target is kept. **All of this is named in the dashboard for the 11-07 weekly to ratify.**

**The fixture values.**
- a = `ROD_RADIUS` = 1 mm.
- **b = 5a = 5 mm.** This is sized so that I65's cylinder, of length 10 b = 50 mm, fits a `-n 2` window.
- **t_c = a²/κ_gel = 6.966667 s.** This is the weekly's control time.
- **t\* = 5 t_c = 34.83333 s.** Chain I's 900 s is 5.2 b²/κ_gel here, and steady to ≈ 1e-13, so a read there would repeat the steady check.
- **`ROD_K` = 6.7 W/(m·K) and `ROD_RHO_C` = 2.33018e6 J/(m³·K)** are Ti-6Al-4V room-temperature values, quoted from memory. No source is in the repo, so they are registered as `# nominal`.
- **Carslaw & Jaeger ch. XIII** is cited from memory, and marked so.

**I64's (Lim).** `MAT-9` is the Robin sphere, not a cylinder. The repo's only homogeneous cylinder is `MAT-15`'s Robin-cooled `robin_finite_cylinder_*`, which takes a whole-volume source. (Lim)'s comparand is therefore the **Dirichlet J0 series with the source confined to the shell**, written in the module. Its roots are `jn_zeros(0, N)/b`, independent of the composite bisection, so (Lim) is not true by construction. It plays the role of `MAT-18`'s (R4).

**I65's control, as written, cannot fire.**
- The weekly's control is "k_r = `K_GEL` misses (I) by the closed form's own difference, asserted at factor 1 on the band". At the weekly's points that difference is below `INTERIOR_BAND`:
  - 0 at the steady limit (the core is uniform at θ(a) whatever k_r is);
  - ≈ 4e-4 at (a, t\*);
  - ≤ 7.4e-3 at (a, t_c).
- The control is therefore read **on the axis at (0, t_c)**, where the predicted difference is ≈ 0.133 (unlogged review arithmetic, printed as predicted).
- **(I0), θ_FEM(0, t_c) against the composite at `INTERIOR_BAND`, is added as a fifth asserted (I) read**, so that the control has a matched gate. No band is set or moved; the weekly ratifies or reverts the added read.

**I65's Dirichlet wall and geometry.**
- `solve_transient_conduction` takes Robin only, so θ = 0 at r = b is a Robin wall at `H_DIR` = 1e7 W/(m²·K). The predicted offset is 2.8e-5 of θ_ss(a), printed.
- Affine cells at r = b would move θ_ss(a) by ≈ −1.7 %, over the band. So the mesh is **curved P2** (`setOrder(2)`, the `IMP-7` step 3 route).

### §9 item 824 — 2026-11-05 (§10 I64, chain I — `IMP-8` step 9, the composite infinite cylinder in closed form; lane 3 → chain step)

**Outcome (07:30 slot, 2026-11-05): BLOCKED, attempt 1 of 4, parked on `attempt/IMP-8-I64-20261105T134551Z` (`c2393f89`).**
- Every gate is green at `B1_BAND`, and so is the control (0.1507 / 0.1396).
- The only red is the precondition: "|F| ≤ 1e-12 of F's pre-cancellation scale", a scale the brief leaves undefined. The executor's |A J0| + |B Y0| collapses to |F| on the gel core, reading 1.0 there.
- Against the envelope √(Z0² + Z1²)(b), max |F| ≤ 2.25e-13 (`20261105T134551Z_IMP-8-item824-W1.log`, `20261105T134627Z_IMP-8-item824-diag.log`).
- **Unblock:** a review defines the scale.

**Executes** §10 I64 under the ruling above. The target and gates are the weekly's. Family I64, attempt 1 of 4.

**Code.** New `tests/validation/test_imp8_composite_cylinder_closed_form.py`, numpy / scipy only. No `src/` edit, and no landed module edited.
- **Imports:**
  - `ROD_RADIUS`;
  - `K_GEL` (`test_mat15_coil_driven_heating_phantom`);
  - `RHO_C` (`test_mat12_transient_sar_heating_sphere`);
  - `B1_BAND` 1e-12 (`test_imp4_pec_cylinder_closed_form`);
  - `INTERIOR_BAND` 1e-2 (`test_mat13_transient_perfusion_sphere`).
- **Template:** `MAT-18`'s two-layer sphere (`analytical.py:2913–3174`): the same root scan and the same Green-form projection, with cylinder functions in place of sin/r.
- **Source:** q = `Q_SHELL` = 1/(π(b² − a²)) W/m³ in a < r < b, which is 1 W per metre of rod. Results are in K per W/m.
- **Steady state:** θ_ss = q(b² − r²)/(4k_g) − q a² ln(b/r)/(2k_g) on a ≤ r ≤ b, and θ_ss(a) for r < a. It does not depend on (k_r, ρc_r). Review arithmetic: θ_ss(a) = 0.1148409 and θ_ss(3 mm) = 0.0827736 K per W/m.
- **Modes**, with μᵢ = √(s/κᵢ):
  - core: J0(μ₁r);
  - shell: Z = A J0(μ₂r) + B Y0(μ₂r), with θ and k∂ᵣθ continuous at a;
  - A and B come from the Wronskian (Jin 3e App. E, (E.18)–(E.19));
  - the eigencondition is F(√s) = Z(b) = 0.
- **Roots:**
  - scan √s at a step of π/(32Λ), with Λ = a/√κ_r + (b − a)/√κ_g;
  - refine each sign change by bisection to 4ε relative;
  - start the scan at x > 0.
- **Series:**
  - Nₙ by Lommel's integral;
  - cₙ = q ∫ₐᵇ Z r dr / (sₙNₙ);
  - θ = θ_ss − Σ cₙψₙe^{−sₙt}, with `N_MODES` 400.

**Gates (asserted, bands imported).**
- **(St)**
  - At r ∈ {0, a, (a+b)/2}, the series at t = 50 b²/κ_gel equals θ_ss within `B1_BAND`.
  - θ_ss's residuals on the scale qb²/k_g: θ_ss(b) = 0, −k_g θ_ss′(a⁺) = 0, and k_g(θ″ + θ′/r) + q = 0 at three shell radii.
- **(En)**
  - For each mode, sₙ∫ρcψₙ r dr = −k_g b ψₙ′(b) (`MAT-18`'s (R3)), relative to the sum of the magnitudes, ≤ `B1_BAND`.
  - The steady balance −k_g b θ_ss′(b) = q(b² − a²)/2.
- **(Lim)** At k_r → `K_GEL` and ρc_r → `RHO_C`:
  - the composite equals the J0 series at r ∈ {0, a, (a+b)/2} and t ∈ {t_c, t\*}, ≤ `B1_BAND`;
  - the first 50 composite roots equal `jn_zeros(0, 50)/b`, ≤ `B1_BAND`.
- **(Bq)** N = 200 against N = 400 at the same six points, ≤ `B1_BAND`.
- **Precondition:**
  - |F| at every root ≤ 1e-12 of F's pre-cancellation scale;
  - the root count is N;
  - the Gram matrix's largest off-diagonal at N = 20 is printed.

**Control (asserted, factor 1).** ρc_r → `RHO_C` alone (k_r = `ROD_K`) misses the composite at t_c by > `INTERIOR_BAND` at r = 0 and at r = a; the minimum of the two is asserted.
- **Predicted ≈ 0.151 / 0.140** (unlogged review arithmetic), printed beside the measurement.
- At (a+b)/2 the miss is ≈ 5e-3, under the band, so it is printed only.

**Printed (🟡 record, no band, ≤ 7 digits).**
- θ in K per W/m at r ∈ {0, a, (a+b)/2} and t ∈ {t_c, t\*, steady}, for three cores:
  - **Ti** (the composite);
  - **gel core**, which is (Lim), the kept-cell rod thermally;
  - **adiabatic hole**, shell only with ψ′(a) = 0, which is the PEC-hole rod.
- The three cores' ordering at r = a.
- The t = 0 tail at N = 400 / 800. It converges as 1/N², so it is never a gate.

**Window.** W1, bare and foreground, Bash timeout 660000 ms. Window form: `20261104T103134Z_IMP-8-item811-W1.log:12`.
`scripts/testing/run_and_log.sh IMP-8-item824-W1 "docker compose -f docker/docker-compose.yml exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 60 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp8_composite_cylinder_closed_form.py -v -s --tb=short'"`

**Tier.** Smoke, `-n 1`, real (the image default), no solve. Predicted < 20 s; `IMP-8`'s closed-form windows took 3–9 s. **Class:** lane 3 → chain step. **Cost:** ≈ 4.1 measured / 20 predicted min.

**Traps already paid for.**
- The imports pull in dolfinx, so the module runs in the container only.
- `-s` is needed, or prints are lost.
- Never evaluate the series at t = 0 (1/N² convergence).
- Composite roots are not one per π interval: scan for them.
- Y0 is singular at 0.
- A predicted factor is printed, never asserted (rule (e)).
- Print ≤ 7 digits.

**Scope.** An infinite cylinder with a Dirichlet wall, a uniform shell source, constant properties and perfect contact. **Not claimed:** a finite rod (I66), a SAR shape, a Ti property (the pair is nominal), or entry 55.

**Negative result.**
- A gate red is a derivation or root-finding error: stop, BLOCKED, never move `B1_BAND`.
- A precondition red is a missed or duplicated root: print the brackets.
- A control under the band means the core's capacity barely matters at t_c: a fixture finding, reported.
- No `PREDICTED RED:` line.

**Status it can move.**
- `IMP-8` gains step 9, a 🟡 record line.
- §10 I64 is DONE.
- Item 837 gets its comparand.

**Depends on:** nothing.

### §9 item 837 — 2026-11-05 (§10 I65, chain I — `IMP-8` step 10, the composite cylinder in FEM; chain step (solve))

**Outcome (16:30 slot, 2026-11-05): DONE, attempt 1 of 4. `IMP-8` step 10 is a 🟡 record; §10 I65 is DONE; I66 is the next review's to write.**
- New `tests/validation/test_imp8_composite_cylinder_fem.py`; no `src/` edit, no landed module edited, not edited after W1.
- **W0** `20261105T225535Z_IMP-8-item837-W0.log` was red on the executor's own facet classifier: the end-face test `|z_max − z_min| < 1e-9` let OCC's tolerance-padded bounding boxes of the end annuli through, so they were tagged as the wall (a tagged geometry node at |r − b| = 4 mm, a tagged midpoint at |z| = L/2). The test is now `< 1e-3 L`; re-run as **W0r** `20261105T225605Z_IMP-8-item837-W0r.log` (8 s, green). This is the fragment-branch trap the brief names, met through a different door.
- **Preconditions (asserted, W0r and W1):** coordinate degree 2; 53 504 cells (cap 250 000; W0b not needed); cell tags {1, 2}; every geometry node of the 918 owned tag-1 facets at r = b to 8.7e-19 m, max |z| of their midpoints 24.55 mm < L/2; P2 wiring max error 2.6e-18. Printed: V_core,h/(πa²L) − 1 = −1.25e-5, V_shell,h − 1 = −3.11e-5.
- **W1** `20261105T225620Z_IMP-8-item837-W1.log`: heavy label, real, `-n 2`, 37 s, Status 0, 7 passed per rank. t_build 4.5 s, transient t_call 22.7 s (2 factorizations), steady 4.0 s. Pre-check `pgrep -c python3` read 0 before each window.
- **Gates (asserted, bands imported):**
  - (I), band 1e-2: (a, t\*) 2.225e-5, ((a+b)/2, t\*) 6.017e-4, (a, steady) 1.324e-5, ((a+b)/2, steady) 2.373e-4;
  - (I0) (0, t_c) 1.657e-3;
  - (E) 3.382e-13 (`E_BAND` 1e-9); (S0) 2.776e-15 (`B_BAND` 1e-10), ∫q dx(1) = 0 exactly; no NaN left.
- **Control (asserted, factor 1 on `INTERIOR_BAND`):** k_r = `K_GEL` misses θ_FEM(0, t_c) by **0.1521**, predicted ≈ 0.133. Printed at the weekly's points: (a, t_c) 9.429e-3 (predicted ≤ 7.4e-3 — over the prediction, still under the band), (a, t\*) 3.286e-4 (predicted ≈ 4e-4), (a, steady) 1.32e-5 (predicted 0; this is the FEM's own steady error, since the swap leaves θ_ss unchanged).
- **Printed:** z-invariance at t\*: θ(a, 0) 7.723997e-2, θ(a, ±L/4) 7.724034e-2 / 7.724128e-2 K per W/m. Robin offset 2.7717e-5 of θ_ss(a) (predicted 2.8e-5). Off-gate t_c reads at (a, 0) 1.669e-3 and ((a+b)/2, 0) 3.739e-3.
- **Surprise:** the control's axis miss (0.152) and its (a, t_c) miss (9.4e-3) both read above the review's predictions (0.133, ≤ 7.4e-3). Both are comparand-only arithmetic (the FEM matches the Ti composite to 1.7e-3 at both points), so the predictions, not the solve, are the likely source; not diagnosed.

**Executes** §10 I65 under the ruling above. Item 824's constants and `_theta` are imported as the comparand. Family I65, attempt 1 of 4.

**Code.** New `tests/validation/test_imp8_composite_cylinder_fem.py`, real build, degree 2, with knob `FEM_EM_IMP8_I65_PHASE` = `build` / `solve`.
- **`_composite_cylinder(comm, *, h_a=0.30e-3, h_b=2.0e-3)`:**
  - gmsh OCC `occ.fragment` of the radius-a and radius-b cylinders, length 10b = 50 mm, centred on z = 0;
  - cell tags: core 1, shell 2, by bbox r_max;
  - facets from `getBoundary` of both volumes, classified **z first**;
  - end faces untagged (Neumann), and the r = b lateral face as `ROBIN_TAG` 1;
  - **do not use `MeshGenerator.cylindrical_domain`'s fragment branch**: it tags the end annuli as the outer wall (`mesh.py:1175–1185`);
  - sizing: a `Distance` field on the core's lateral surface into one `Threshold` (h_a out to a, up to h_b at b − a);
  - `setOrder(2)`, then `_model_to_mesh`, on the pattern of `_graded_box` (`test_imp8_rod_line_source_fem.py:87`).
- **Properties:** k and ρc as DG0, via `_dg0(msh, ct, ROD_K, K_GEL)` / `_dg0(msh, ct, ROD_RHO_C, RHO_C)`, and q = `_dg0(msh, ct, 0.0, Q_SHELL)` (`test_mat18_two_layer_sphere_fem.py:77`). Assert that no NaN is left.
- **Sampling:** `_theta_pts` (`test_imp8_rod_line_source_fem.py:180`) at (r, z) mm = (1, 0), (3, 0), (0, 0).
- **Solves:**
  - one CN run: `theta=0.5`, `startup_half_steps=2`, N = 80, Δt = t\*/80, and `snapshot_steps=(16,)` for t_c;
  - the steady state: `solve_steady_conduction` with DG0 k and h = `H_DIR`.

**Gates (imported, unmoved).**
- **(I):** |θ_FEM/θ_cf − 1| ≤ `INTERIOR_BAND` at (a, 0) and ((a+b)/2, 0), at t\* and at the steady limit.
- **(I0):** the same at (0, 0), t_c (the ruling above).
- **(E):** `max_balance_residual` ≤ `E_BAND`.
- **(S0):** ∫q dx = `Q_SHELL`·V_shell,h within `B_BAND` (all-reduced), and ∫q dx(1) = 0.
- **Preconditions (W0):**
  - the coordinate element is degree 2;
  - cells ≤ `CELL_CAP` 250 000;
  - the tag-1 facet midpoints lie at r = b;
  - f = x + 2y + 3z, in P2, evaluates at the three read points to ≤ 1e-12 (a wiring check of point evaluation on curved cells).

**Control (asserted, factor 1 on `INTERIOR_BAND`).** θ_cf with k_r = `K_GEL` misses θ_FEM at (0, t_c) by more than the band.
- **Predicted ≈ 0.133**, printed beside the measurement.
- The same miss is printed at the weekly's points, predicted ≤ 7.4e-3, ≈ 4e-4 and 0.
- It is a comparand swap and costs no solve.

**Printed.**
- Cells.
- V_core,h/(πa²L) − 1 and V_shell,h − 1.
- θ_FEM(a, ±L/4) beside (a, 0), a z-invariance witness.
- The Robin offset (predicted 2.8e-5).
- `t_build`, `t_call`, `n_factorizations`.

**Windows.** Pre-check: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0. Every window is a bare foreground harness call, Bash timeout 660000 ms, in the durable-capture form of `20261103T182153Z_IMP-8-I55-W1s.log:12`, including the trailing `; exit $rc`.
- **W0** (`build`, `-n 2`, `timeout -k 30 300`): the preconditions and the census. If the count is over 250 000, run W0b once at (0.40, 2.5) mm; if that is still over, BLOCKED.
- **W1** (`solve`, `-n 2`, `timeout -k 30 600`): the gates.
- **Price:**
  - build ≈ 30–90 s (the `IMP-7` step 3 curved build took 93.3 s, `20261026T003804Z_IMP-7-step3-W2.log`);
  - each factorization ≈ 15–60 s (item 800's steady solve took 284.0 s on 485 291 P2 cells at `-n 2`, scaled);
  - W1 ≈ 120–300 s in total.

**Tier.** Heavy, real, `-n 2`. **Class:** chain step (solve). **Cost:** ≈ 11.2 measured / 40 predicted min.

**Traps already paid for.**
- The ruling's three premises.
- A gmsh `Min` of `Threshold`s needs `StopAtDistMax` (`20261103T105626Z_IMP-8-I55-W0b-probe.log`); this build uses one threshold.
- Pin the quadrature degree on any facet integral that uses `SpatialCoordinate`.
- An orphan `.c` sweep runs only when `pgrep` reads 0.
- Reduce rank-local integrals before asserting.
- Count pytest verdicts per rank, by token.
- `.eval` bytes depend on the batch.
- Keep Rannacher's start at the near-Dirichlet wall.
- A W1 kill is never answered by raising the timeout.

**Scope.** One mesh, one Δt, an infinite-cylinder comparand. **Not claimed:** a rate, the finite rod (I66), entry 55, or a Ti property.

**Negative result.**
- (I) red with (E) / (S0) green: print the volume errors, stop, BLOCKED.
- (I0) red alone: the early transient on the axis is under-resolved. Report it and stop.
- A control under the band is a finding: report it.
- No `PREDICTED RED:` line.

**Status it can move.**
- `IMP-8` gains step 10, a 🟡 record line: the first FEM gate with an implant's thermal properties in it.
- §10 I65 is DONE.
- I66 becomes the next review's to write.

**Depends on:** item 824 landing. If 824 is BLOCKED, skip at no cost. *(10:30 review, 2026-11-05: 824 is BLOCKED and re-posed as item 843; read "item 843 landing".)*

### Review ruling 2026-11-05 10:30 — item 824 (I64) re-posed as item 843 under rule (6)

The red is on a precondition the 03:00 brief worded and did not define: "|F| at every root ≤ 1e-12 of F's pre-cancellation scale". The executor's reading, |A J0(μ₂b)| + |B Y0(μ₂b)|, is |F| itself wherever B = 0, so on the gel core (Lim's homogeneous case) the ratio is identically 1.0 and the precondition cannot pass for any root finder. **Defined here (the band 1e-12 unchanged):** the scale is the outer cylinder function's local amplitude at b,

  env(√s) = √( Z₀(μ₂b)² + Z₁(μ₂b)² ),  Z_ν = A J_ν + B Y_ν with the same A, B as F,

which does not vanish at a root of Z₀ (Z₀ and Z₁ have no common zero), so |F| / env measures how close the refined root is to the zero, which is what the precondition is for. Attempt 1's diagnostic read max |F| ≤ 2.25e-13 against exactly this envelope on all four cores (`20261105T134627Z_IMP-8-item824-diag.log`), so a green is expected; that diagnostic is the precondition's measurement, not a gate. No gate, band or target moves, and this is not a new §10 step.

### §9 item 843 — 2026-11-05 (§10 I64, chain I — re-pose of item 824 under rule (6), attempt 2 of 4; `IMP-8` step 9, the composite infinite cylinder in closed form; lane 3 → chain step)

**Outcome (12:00 slot, 2026-11-05): DONE, attempt 2 of 4. `IMP-8` step 9 is a 🟡 record; §10 I64 is DONE; item 837 has its comparand.**
- W1 `20261105T181438Z_IMP-8-item843-W1.log`: smoke, `-n 1`, 3 s, Status 0, 2 passed / 3 skipped (complex-mode checks skip in the real build).
- The only code change from attempt 1 is the precondition's scale, now env(√s) per the 10:30 ruling (Z₁ with the same A, B; B·Y₁ only where B ≠ 0).
- **Precondition (asserted):** max |F|/env = 2.037e-13 (Ti), 2.245e-13 (gel), 1.987e-13 (hole), 1.950e-13 (ctrl), all ≤ 1e-12. Every core has 400 roots. The N = 20 Gram off-diagonal is 1.850e-15 (Ti) and 2.243e-15 (hole).
- **Precondition control (printed):** the attempt-1 scale |A J0| + |B Y0| reads 7.282e-12 / 1.000 / 2.801e-13 / 1.097e-11. The gel core's 1.0 is the definitional failure the envelope removes.
- **Gates (asserted, `B1_BAND` 1e-12):**
  - (St) series 0 and residuals ≤ 1.607e-17;
  - (En) modes 3.151e-16 and steady 0;
  - (Lim) θ 3.625e-16 and roots 2.464e-16;
  - (Bq) 0.
  These reproduce attempt 1's digits.
- **Control (asserted, factor 1 on `INTERIOR_BAND`):** 0.1507 / 0.1396 at r = 0 / a, against the predicted 0.151 / 0.140. The (a+b)/2 value, 0.00494, is printed only.
- Attempt 1's two logs and their test-results rows were carried onto `main` by hand. `scripts/diagnostics/imp8_i64_pre_diag.py` was carried as a record and left unedited. The attempt branch is kept for the review to dispose of, per the slot's prompt.

**Executes** §10 I64 under the 10:30 ruling above. Item 824's brief binds wherever this one is silent — every gate, band, control, print and scope line. Family I64, attempt 2 of 4.

**Code.**
- Take `tests/validation/test_imp8_composite_cylinder_closed_form.py` (and `scripts/diagnostics/imp8_i64_pre_diag.py`, as a record) from `attempt/IMP-8-I64-20261105T134551Z` (`c2393f89`) by path — new files, so a path checkout is safe. Change only the precondition's scale to env(√s) as defined above, printing max |F| / env per core. Nothing else.
- Carry attempt 1's two logs (`20261105T134551Z_IMP-8-item824-W1.log`, `20261105T134627Z_IMP-8-item824-diag.log`) onto `main` by path and append their test-results rows by hand (insertions only).

**Window.** Item 824's W1 form, new tag: `scripts/testing/run_and_log.sh IMP-8-item843-W1 "docker compose -f docker/docker-compose.yml exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 60 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp8_composite_cylinder_closed_form.py -v -s --tb=short'"` — bare, foreground, Bash timeout 660000.

**Anchors.** Item 824's: (St), (En), (Lim) (composite = homogeneous Dirichlet J0 series; 50 roots = `jn_zeros(0, 50)/b`), (Bq), all ≤ `B1_BAND`; the precondition (|F| / env ≤ 1e-12, root count N, the N = 20 Gram off-diagonal printed). Attempt 1 read every gate green (`20261105T134551Z_IMP-8-item824-W1.log`), and those digits are the reproduction expectation.

**Control (asserted, factor 1).** Item 824's: ρc_r → `RHO_C` alone misses the composite at t_c by > `INTERIOR_BAND` at r = 0 and r = a (attempt 1: 0.1507 / 0.1396, predicted 0.151 / 0.140). **The precondition's own control (printed, the reason for the ruling):** the executor's attempt-1 scale |A J0| + |B Y0| on the gel core prints 1.0 — the definitional failure the envelope removes.

**Tier.** Smoke, `-n 1`, real, no solve. Attempt 1's W1 took 3 s. **Class:** lane 3 → chain step. **Cost:** ≈ 4.1 measured / 10 predicted min.

**Traps already paid for.** Item 824's list. Z₁ uses the same A, B as Z₀ (never re-normalised per core); Y₁ is singular at 0 exactly as Y₀ — the gel core has B = 0, so evaluate B·Y₁ only where B ≠ 0.

**Scope.** As item 824. **Not claimed:** a finite rod (I66), a SAR shape, a Ti property, or entry 55.

**Negative result.** A gate red: stop, BLOCKED, never move `B1_BAND`. A precondition red against the envelope is a missed or duplicated root: print the brackets and park on a new `attempt/IMP-8-I64-*`. No `PREDICTED RED:` line.

**Status it can move.**
- `IMP-8` gains step 9, a 🟡 record line; §10 I64 is DONE.
- Item 837 (I65) gets its comparand.

**Depends on:** nothing. Delete `attempt/IMP-8-I64-20261105T134551Z` in the landing slot once its files and logs are on `main`.

*(18:00 review, 2026-11-05: the landing slot left the branch for the review; every file on it — the module, `scripts/diagnostics/imp8_i64_pre_diag.py`, both attempt-1 logs and their two test-results rows — is on `main` (`git diff --stat main...branch` lists those five paths only), so the review deleted it.)*

### Review ruling 2026-11-05 18:00 — chain I's I66 as execution detail, and three provisional steps

Writing the brief for the 11-04 weekly's I66, this review fixed what the step text leaves open. The target (a record, (E), (S0)) is kept. **Named in the dashboard for the 11-07 weekly to ratify or revert.**

- **The rod cells and the source overlap.** The step says "the rod cells (r ≤ `ROD_RADIUS`, |z| ≤ L/2) given (`ROD_K`, `ROD_RHO_C`) … and I56's two-tip source in the gel around them". On I55's box both sets are midpoint selections, and I56's two balls (midpoint within `ROD_RADIUS` of (0, 0, ±L/2)) overlap the rod cylinder's end caps. **Ruled:** a cell in both sets keeps the gel's properties and carries the source. The rod is the cylinder minus the two ball sets. I56's source is then byte-identical (same `_tips_source` call), it lies wholly in gel-property cells as the step says, and the only change from I56 is the material of the remaining rod cells. The overlap's cell count and volume are printed.
- **The window split.** One CN run on item 803's box costs 439.6 s of `t_call` at `-n 8` (`20261104T043245Z_IMP-8-I56-W1-item796.log:104`), so one arm per window: **W1** the Ti arm, **W2** the gel arm as the DG0 path's reproduction control (it must return I56's six printed θ).
- **The negative control.** A DG0 path that ignored its coefficients would return I56's θ unchanged on the Ti arm. So the Ti arm's θ must differ from I56's at (2, 0) mm by more than 1e-9 relative (the round-off floor; I56's (E) is 1.754e-14). The size of the move is predicted, not asserted: the rod's axial conductance k_r·πa²/L ≈ 7.0e-4 W/K against a 1 mm ball's spreading conductance 4πk_g·a ≈ 7.5e-3 W/K (unlogged review arithmetic) puts it at the per-cent level near the rod.
- **Three provisional steps under rule (6)** (chain I at 3 of 6 after them; each keeps (7e)'s target "implant-driven heating in closed form, then FEM", sets no band): **I69** the near-adiabatic rod arm on I66's module (I64's "adiabatic hole" core, the bracket); **I70** the composite cylinder with the source in the core, closed form (the implant as the heat source); **I71** I70's FEM twin on I65's curved-P2 cylinder. **I72**, an `imp:` example of I64 / I65, is a fourth (`EX-184`).

### §9 item 846 — 2026-11-05 (§10 I66, chain I — `IMP-8` step 11, I56's rod as a conducting body; chain step (solve))

**Outcome (2026-11-05, 19:30 slot, implementer): DONE.** Attempt 1 of 4. `IMP-8` step 11 is a 🟡 record; §10 I66 is DONE; item 851 is unblocked.
- New `tests/validation/test_imp8_rod_conducting_body_fem.py`, imports as the brief lists them. No `src/` change, no landed module edited, and the module was not edited after W0.
- W0 `20261106T013504Z_IMP-8-item846-W0.log` (`-n 2`, 5 s, 1 passed). W1 `20261106T013519Z_IMP-8-item846-W1.log` (arm `ti`, `-n 8`, 433 s, 5 passed / 1 skipped per rank). W2 `20261106T014244Z_IMP-8-item846-W2.log` (arm `gel`, `-n 8`, 420 s, 5 passed / 1 skipped per rank). `pgrep -c python3` read 0 before each heavy window.
- **Census:** 485 291 cells (= `W0C_CELLS`). The cylinder has 6 495 cells. The rod has 6 186 cells, 9.036954e-08 m³ (0.958851 of πa²L). The overlap has 309 cells, 4.288928e-09 m³, with gel properties and the source. Balls 294 / 300 cells as I56.
- **Gates:** (E) 1.554e-14 (Ti) / 1.910e-14 (gel) against 1e-9; (S0) 4.441e-16 against 1e-10 on both arms; per-ball ∫q/(P/2) − 1 ≤ 4.4e-16; non-finite coefficient dofs 0.
- **Reproduction (W2, asserted):** gel-arm |θ − θ₅₆| = 8.6e-8, 1.4e-8, 4.7e-7, 3.5e-9, 1.4e-6, 2.9e-7 at (10, 0), (20, 0), (10, 15), (20, 15), (2, 15), (2, 0) mm, each under its sixth-digit unit.
- **Record (W1, Ti, K/W):** θ = 2.365053, 0.7194612, 3.506584, 0.7502096, 25.25271, 4.735218 in that order. θ/θ₅₆ = 1.223938, 1.128767, 0.961888, 0.994523, 0.839660, 1.555474. Against q′/(4πk)·Φ₂ (not a gate: the comparand has no rod) the ratios are 1.226738, 1.131326, 0.962465, 0.995784, 0.841491, 1.559349. The tip-point ratio θ(2, 15)/θ(2, 0) is 5.332956 (Ti) against 9.879328 (gel arm) / 9.879329 (I56). t_call 412.1 s (Ti) / 399.0 s (gel); n_factorizations 2.
- **Control (asserted):** θ_ti(2, 0)/θ₅₆ − 1 = +5.555e-01, against the 1e-9 floor.
- **Findings, not asserted:**
  - The (2, 0) mm move, +55.5 %, is ≈ 6× the review's predicted ≈ 0.09 (unlogged axial-conductance arithmetic), printed only.
  - The gel arm reads θ(2, 0)/θ₅₆ − 1 = 9.4e-8, so it would itself clear the brief's "1e-9 round-off floor". The cause is the comparand: θ₅₆ is I56's 7-digit print, which carries up to ≈ 1.6e-7 relative rounding at 3.044228. The 1e-9 floor is therefore below what that comparand resolves. This window does not measure the DG0-vs-`fem.Constant` path difference. The Ti move is 6.7 decades above both numbers, so the control's verdict stands.

**Executes** §10 I66, the 11-04 weekly's enumerated step, under the ruling above. Family I66, attempt 1 of 4.

**Code.** New `tests/validation/test_imp8_rod_conducting_body_fem.py`, real build, degree 2. No `src/` edit, no landed module edited.
- **Imports, unchanged:** the `box` fixture, `GATE_PTS`, `P`, `QP`, `T_INF`, `W0C_CELLS`, `DEGREE`, `_need`, `_theta_pts` (`test_imp8_rod_line_source_fem`); `_tips_source`, `RATIO_PTS`, `SIGNS` (`test_imp8_rod_tips_fem`); `L`, `MM`, `T_PRINT`, `_phi_tips` (`test_imp8_rod_line_source_closed_form`); `ROD_K`, `ROD_RHO_C` (`test_imp8_composite_cylinder_closed_form`, both `# nominal`); `ROD_RADIUS`; `K_GEL`, `H_GEL`, `N_MAIN`; `RHO_C`, `E_BAND`; `B_BAND`, `ROBIN_TAG`, `_allsum`, `_say`; `INTERIOR_BAND`.
- **Knob** `FEM_EM_IMP8_I66_ARM` = `ti` (default) / `gel` / `hole` (`hole` is item 851's; this item refuses it with a message naming item 851 until 851's brief is executed). The phase knob `FEM_EM_IMP8_I55_PHASE=transient` is set on the command, as I56's was.
- **`_rod_cells(msh)`**: midpoints over owned + ghost cells (as `_tips_source` builds them); rod = r ≤ `ROD_RADIUS` ∧ |z| ≤ L/2; ball = within `ROD_RADIUS` of (0, 0, ±L/2); rod cells = rod ∧ ¬ball. Return the index arrays; print owned counts allreduced (rod, overlap) and the rod-cell volume against πa²L.
- **Coefficients:** DG0 `k` = `K_GEL` everywhere then `ROD_K` on rod cells; DG0 `rho_c` = `RHO_C` then `ROD_RHO_C` on rod cells (arm `gel`: no overwrite — the same DG0 objects, gel values). `scatter_forward`; assert no NaN.
- **Solve:** exactly I56's call (`solve_transient_conduction(msh, ft, k=k, rho_c=rho_c, q=q, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN_TAG, dt=T_PRINT / N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, degree=DEGREE)`), q from `_tips_source(msh)`.
- **Reads:** `_theta_pts` at `GATE_PTS + RATIO_PTS` in that order, one batch (the batch rule).

**Gates (asserted, bands imported).**
- **(E)** `max_balance_residual` ≤ `E_BAND`.
- **(S0)** ∫q = P and each ball's ∫q = P/2 within `B_BAND`, cell counts > 0 (I56's test, same numbers: 294 / 300 cells).
- **Preconditions:** cells = `W0C_CELLS` (485 291); rod cells > 0 and overlap cells > 0 (both allreduced).
- **W2 only (arm `gel`), the reproduction control:** each of the six θ equals I56's printed value within one unit of its sixth significant digit, |θ − θ₅₆| ≤ 10^(⌊log₁₀|θ₅₆|⌋ − 5). I56's values (`20261104T043245Z_IMP-8-I56-W1-item796.log:106–115`, K/W): (10, 0) 1.932331e+00, (20, 0) 6.373868e-01, (10, 15) 3.645521e+00, (20, 15) 7.543411e-01, (2, 15) 3.007493e+01, (2, 0) 3.044228e+00.
- **W1 only (arm `ti`), the negative control:** |θ_ti(2, 0) / 3.044228 − 1| > 1e-9.

**Printed (🟡 record, no band, ≤ 7 digits).** On W1: θ at the six points beside I56's, the ratio θ_ti / θ₅₆ per point; θ_ti against q′/(4πk)·Φ₂ (I56's comparand) per point, labelled "not a gate: the comparand has no rod"; the tip-point ratio θ(2, 15) / θ(2, 0) on both arms; rod and overlap cell counts and volumes; `t_call`, `n_factorizations`. The axial-conductance prediction (≈ 0.09, unlogged review arithmetic) beside the (2, 0) mm move.

**Windows.** Pre-check `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 before each. Bare foreground harness calls, Bash timeout 660000 ms, durable capture as `20261103T182153Z_IMP-8-I55-W1s.log:12` including the trailing `; exit $rc`, container `timeout -k 30 600`, `mpiexec -n 8`, env `FEM_EM_IMP8_I55_PHASE=transient` plus the arm knob.
- **W0** (`-n 2`, `timeout -k 30 120`): `-k k0` only — a knob test that the arm parser accepts `ti` / `gel` and refuses `hole` and junk (no mesh).
- **W1** arm `ti`, **W2** arm `gel`. Each ≈ 460–480 s (I56: 461 s total, 439.6 s `t_call`, 18.4 s build). A kill parks the item; never raise the ceiling.

**Tier.** Heavy, real, `-n 8`. **Class:** chain step (solve). **Cost:** ≈ 5.5 measured / 25 predicted min (two 8-minute windows).

**Traps already paid for.**
- `_tips_source` builds `one * w * ufl.dx`, which the complex build rejects (`ArityMismatch`, `20261105T193522Z_EX-182-figure.log:104`): this module runs in the **real** build only, as I56 did. Do not source the complex mode.
- Two heavy windows never at once; `pgrep` reads 0 first.
- Reduce every rank-local count and integral before asserting; count verdicts per rank by token.
- `.eval` bytes depend on the batch: read the six points as one batch, in I56's order, or W2's reproduction can miss in the last digit for a reason that is not the physics.
- `-s` on pytest, or the prints are lost.

**Scope.** One mesh, one Δt, nominal Ti properties, perfect thermal contact, a rod resolved only as I55's grading resolves it. **Not claimed:** a gate on the rod's effect (no closed form has the finite rod), entry 55, a Ti property, I69's bracket.

**Negative result.** (E) or (S0) red: stop, BLOCKED. W2's reproduction red: the DG0 path is not I56's path — print every θ, stop, BLOCKED (W1's record is not banked without it). W1's control red (θ unmoved): the coefficient never reached the form — BLOCKED. No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 11, a 🟡 record line (the first number for what the rod's own conduction does to the tip temperature); §10 I66 DONE; item 851 unblocked.

**Depends on:** nothing (I56's records and I65 are on `main`).

### §9 item 847 — 2026-11-05 (§10 I70, chain I — daily, 2026-11-05 18:00, provisional under rule (6); `IMP-8` step 12, the composite cylinder with the source in the core, closed form; lane 3 → chain step)

**Outcome (2026-11-05, 19:30 slot, implementer): DONE.** Attempt 1 of 4. `IMP-8` step 12 is a 🟡 record; §10 I70 is DONE; item 852 has its comparand.
- New `tests/validation/test_imp8_composite_cylinder_core_source_closed_form.py`, imports as the brief lists them plus item 843's `_theta_ss` (aliased `_theta_ss_shell`, used only for the printed shell-vs-core comparison). No `src/` change, no landed module edited, not edited after W1.
- W1 `20261106T015337Z_IMP-8-item847-W1.log` (`-n 1`, real, smoke, 3 s, 2 passed / 3 skipped).
- **Gates:** (St) series at `T_ST` vs θ_ss 0 at all three radii; residuals θ(b) 0, θ continuity at a 4.478e-14, flux continuity 0, PDE ≤ 1.8e-16 at 0.3 / 0.7 / 2 / 4 mm (scale q_c a²/k_g = 0.5305165, PDE on q_c, flux on q_c a/2). (En) −k_g b θ′(b) = 0.1591549 = q_c a²/2, rel 0. (Lim) ≤ 2.093e-16. (Bq) 0.
- **Control (asserted > 1e-2):** k_r → `K_GEL` alone, steady θ(0) (series at `T_ST`) 0.5595458 against 0.4387939, miss 0.2752 (predicted 0.2752). At r = a the miss is 0 (printed).
- **Record (K per W/m), θ at r = 0, a, (a+b)/2:** Ti t_c 0.1935047, 0.1827126, 0.009287464; t\* 0.3551997, 0.3435803, 0.08512798; steady 0.4387939, 0.4269167, 0.1355007. Gel core t_c 0.2558521, 0.1510219, 0.006445654; t\* 0.4552613, 0.3286419, 0.07908581; steady 0.5595458, 0.4269167, 0.1355007. The ρc-only swap at t_c on the axis misses by 0.1129. θ(a, steady) is 0.4269167 with the source in the core and 0.1148409 with item 843's shell source (ratio 3.717461).

**Executes** §10 I70 under the ruling above. Family I70, attempt 1 of 4.

**Code.** New `tests/validation/test_imp8_composite_cylinder_core_source_closed_form.py`, numpy / scipy. No edit to item 843's module; import from it `A_ROD`, `B_WALL`, `KAPPA_GEL`, `T_C`, `T_STAR`, `T_ST`, `ROD_K`, `ROD_RHO_C`, `N_MODES`, `R_PTS`, `_modes`, `_psi`, `_gram`.
- **Source:** q_c = `Q_CORE` = 1/(πa²) W/m³ in r < a (1 W per metre of rod), zero in the shell, θ = 0 at r = b, from θ = 0.
- **Steady state** (write it; review arithmetic, unlogged): core θ_ss = q_c(a² − r²)/(4k_r) + q_c a² ln(b/a)/(2k_g); shell θ_ss = q_c a² ln(b/r)/(2k_g). It depends on k_r in the core only. Predicted θ_ss(0) = 0.43879 (Ti) / 0.55955 (gel core) and θ_ss(a) = 0.42692 K per W/m.
- **Series:** `_modes(ROD_K, ROD_RHO_C, N_MODES)` already returns `int_core` = ∫₀ᵃ J0(μ₁r) r dr and `norm`; cₙ = q_c·int_core/(sₙ·normₙ); θ = θ_ss − Σ cₙψₙe^{−sₙt}.

**Gates (asserted, bands imported: `B1_BAND` 1e-12 from `test_imp4_pec_cylinder_closed_form`, `INTERIOR_BAND` 1e-2).**
- **(St)** the series at t = `T_ST` equals θ_ss at r ∈ `R_PTS` within `B1_BAND`; θ_ss's residuals on the scale q_c a²/k_g: θ_ss(b) = 0, flux continuity k_r θ′(a⁻) = k_g θ′(a⁺), θ continuous at a, and k(θ″ + θ′/r) + q = 0 at two core and two shell radii.
- **(En)** −k_g b θ_ss′(b) = q_c a²/2 relative ≤ `B1_BAND`; the modal identity is item 843's (same modes, not re-asserted).
- **(Lim)** k_r → `K_GEL`, ρc_r → `RHO_C`: the composite equals the homogeneous Dirichlet J0 series with the source confined to r < a — λₙ = `jn_zeros(0, N)/b`, cₙ = q_c·a J1(λₙa)/λₙ / (κ_gel λₙ² · RHO_C · ½b²J1(λₙb)²) — at r ∈ `R_PTS`, t ∈ {`T_C`, `T_STAR`}, ≤ `B1_BAND`.
- **(Bq)** N = 200 vs 400 at the same six points ≤ `B1_BAND`.

**Control (asserted, factor 1 on `INTERIOR_BAND`).** k_r → `K_GEL` alone (ρc_r = `ROD_RHO_C`) misses the composite's steady θ(0) by > `INTERIOR_BAND`. Predicted 0.55955 / 0.43879 − 1 = 0.2752 (unlogged), printed beside. At r = a the steady miss is 0 by construction (printed, the reason the control is read on the axis).

**Printed (record, ≤ 7 digits).** θ at r ∈ `R_PTS`, t ∈ {`T_C`, `T_STAR`, steady} for Ti and gel cores; the ρc-only swap's miss at t_c on the axis; the source-in-core vs source-in-shell (item 843) θ(a, steady) side by side per W/m.

**Window.** W1, smoke, `-n 1`, real, `timeout -k 30 60`, the form of `20261105T181438Z_IMP-8-item843-W1.log:12` with this module. Predicted < 10 s.

**Tier.** Smoke. **Class:** lane 3 → chain step. **Cost:** ≈ 4.1 measured / 15 predicted min.

**Traps already paid for.** Item 843's list (never t = 0; Y0 singular at 0; `-s`; a predicted factor printed, never asserted beyond the band; ≤ 7 digits). `_modes` is `lru_cache`d on its arguments: pass floats, as `_theta` does.

**Scope.** Infinite cylinder, Dirichlet wall, uniform core source, constant properties. **Not claimed:** a realistic implant source (Joule heating in a Ti rod is a SAR shape no step has), a finite rod, a Ti property.

**Negative result.** A gate red is a derivation error: stop, BLOCKED, `B1_BAND` unmoved. A control under the band contradicts the closed-form arithmetic above: report it. No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 12, a 🟡 record line; §10 I70 DONE; item 852 gets its comparand.

**Depends on:** nothing.

### §9 item 851 — 2026-11-05 (§10 I69, chain I — daily, 2026-11-05 18:00, provisional under rule (6); `IMP-8` step 13, the near-adiabatic rod on I66's box; chain step (solve))

**Outcome (2026-11-06, 21:00 slot, implementer): DONE.** Attempt 1 of 4. `IMP-8` step 13 is a 🟡 record; §10 I69 is DONE.
- **Code:** item 846's `tests/validation/test_imp8_rod_conducting_body_fem.py`, arm `hole` enabled (`ARMS`, the k0 knob test no longer refuses it). On the rod cells: k = 1e-6·`K_GEL` and ρc = 1e-6·`RHO_C` (`EPS_HOLE`). There is a new `test_w3_hole_control_moves`, item 846's Ti θ is recorded as `TI_THETA` (cited to its W1 log) for the printed bracket, and the module docstring is updated. No `src/` change. No other gate changed. Not edited after W0 / W3.
- W0 `20261106T031512Z_IMP-8-item851-W0.log` (`-n 2`, k0, 1 passed). W3 `logs/leak-held/20261106T031528Z_IMP-8-item851-W3.log` (gitignored, uncommitted: the commit hook's leak guard matched its line 117, one of the six printed θ lines of this closed-form-free gel-box solve, with no AED file opened in this slot, so it reads as a coincidental digit match; moved unchanged per the standing rule, for the operator to rule) (arm `hole`, `-n 8`, real, 438 s, 5 passed / 2 skipped on each of 8 ranks, `[capture] rc=0`). `pgrep -c python3` read 0 before W3.
- **Gates:** (E) 2.021e-14 against 1e-9; (S0) 4.441e-16, balls 294 / 300 cells, per-ball ≤ 4.4e-16; cells 485 291; rod 6 186 / overlap 309 cells as item 846; non-finite dofs 0. **Control (asserted):** θ_hole(2, 0)/θ₅₆ − 1 = −8.4602e-02, against the 1e-9 floor. ε = 1e-6 did not degrade (E): 2.0e-14 vs the Ti arm's 1.6e-14.
- **Record (K/W):** θ_hole = 1.829011, 0.6164280, 3.672656, 0.7566369, 32.36154, 2.786680 at (10, 0), (20, 0), (10, 15), (20, 15), (2, 15), (2, 0) mm. hole/gel = 0.946531, 0.967118, 1.007443, 1.003043, 1.076030, 0.915398. Tip-point ratio θ(2, 15)/θ(2, 0) = 11.61 (hole) against 9.879 (gel) and 5.333 (Ti). t_call 416.7 s, n_factorizations 2.
- **Ordering (printed, a finding, not a gate):**
  - At z = 0 (r = 2, 10, 20 mm) the FEM orders **Ti > gel > hole**. At z = 15 mm it orders hole > gel > Ti.
  - I64's closed form at r = a reads hole > Ti > gel at t_c and t\*. That case is an infinite cylinder with a shell source, so it is not this geometry and the comparison is printed only.
  - Here the source sits at the tips. The rod's axial conductance carries tip heat to the mid-plane, raising θ at z = 0 for Ti and lowering it for the hole. The hole also confines heat near the tips, which makes the tip-point ratio larger.
- **Disclosed:** the shared `[imp8-i66] (2, 0) mm move` print line still says "asserted on arm ti". On arm `hole` the same quantity is asserted by `test_w3_hole_control_moves`. The text is a print only and was left as committed (rule (i)). `_rod_cells` still builds `one * w * ufl.dx` (real-only, line 125 before this edit). That is outside `OPS-82`'s two helpers and is not changed here.

**Executes** §10 I69 under the ruling above. Family I69, attempt 1 of 4.

**Code.** Item 846's module, arm `hole` enabled: on the rod cells, k = ε·`K_GEL` and ρc = ε·`RHO_C` with **ε = 1e-6** (I64's adiabatic-hole core in the FEM's terms — the rod takes neither heat nor flux; the source cells keep gel properties, as the ruling above). Only the arm parser and the coefficient line change; item 846's gates stand.

**Gates.** (E), (S0) and the preconditions as item 846. **Negative control:** |θ_hole(2, 0) / θ_gel(2, 0) − 1| > 1e-9, with θ_gel = I56's 3.044228e+00.

**Printed (record).** θ_hole at the six points beside θ_ti (item 846's W1 log) and I56's θ_gel; the ordering at each point; I64's closed-form ordering at r = a (Ti / gel / hole, `20261105T181438Z_IMP-8-item843-W1.log`) beside it, labelled "infinite cylinder, shell source — not the same geometry".

**Window.** W3, arm `hole`, as item 846's W1 (`-n 8`, `timeout -k 30 600`, real). Predicted ≈ 460–480 s; ε = 1e-6 makes the matrix less well conditioned, so if the direct solve's (E) comes back above `E_BAND`, that is the finding (not a ceiling to raise).

**Tier.** Heavy, real, `-n 8`. **Class:** chain step (solve). **Cost:** ≈ 5.5 measured / 15 predicted min.

**Traps.** Item 846's list.

**Scope.** A bracket by FEM against FEM on one mesh. **Not claimed:** the PEC-hole mesh (the rod is not a mesh hole here), a gate on the bracket.

**Negative result.** (E) red at ε = 1e-6: report the residual, stop, BLOCKED — a review may re-pose ε, never the band. No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 13, a 🟡 record line; §10 I69 DONE.

**Depends on:** item 846 landing (its module and W1 log). If 846 is open or BLOCKED, skip at no cost.

### §9 item 852 — 2026-11-05 (§10 I71, chain I — daily, 2026-11-05 18:00, provisional under rule (6); `IMP-8` step 14, the core-source composite cylinder in FEM; chain step (solve))

**Outcome (2026-11-06, 21:00 slot, implementer): DONE.** Attempt 1 of 4. `IMP-8` step 14 is a 🟡 record; §10 I71 is DONE.
- **Code:** a new module `tests/validation/test_imp8_composite_cylinder_core_source_fem.py`, real build, degree 2. It imports `_composite_cylinder`, `H_DIR`, `CELL_CAP`, the CN settings, `READ_PTS` / `WIT_PTS` and `SNAP` from item 837's module, `_dg0` from `test_mat18`, and `_theta` / `_theta_ss` / `Q_CORE` / `PRED_CTRL` from item 847's module, all unmoved. The only physics change is q = `_dg0(msh, ct, Q_CORE, 0.0)`. Item 837's W0 preconditions are re-run inside the solve fixture (this module has no build-only phase). Knob `FEM_EM_IMP8_I71_PHASE` = `solve`. No `src/` edit and no landed module edited.
- **Window:** W1 `20261106T032719Z_IMP-8-item852-W1.log` (`-n 2`, real, 32 s, 7 passed per rank, `[capture] rc=0`). `pgrep` read 0 first.
- **Preconditions:** 53 504 cells, coordinate degree 2, tag-1 geometry nodes within 8.7e-19 m of r = b, max |z| of tagged facets 24.55 mm < L/2, P2 wiring 2.6e-18.
- **Gates:**
  - (I) relative error at (1, 0) / (3, 0) mm: t\* 5.426e-04 / 1.713e-03; steady 3.979e-04 / 1.423e-03.
  - (I0) at (0, 0): t\* 5.458e-04; steady 4.044e-04.
  - (E) 3.009e-13; (S0) 3.553e-15, with ∫q dx(2) = 0.
- **Control (asserted):** the k_r = `K_GEL` steady θ(0) comparand 0.5595458 misses θ_FEM 0.4386164 by 0.2161. The predicted 0.2752 is the comparand-over-composite ratio; the gate reads θ_FEM/comparand − 1, whose predicted value is 1/1.2752 − 1 = −0.2158.
- **Printed:** θ_FEM at (0, 0) is 0.1933441 at t_c, 0.3550059 at t\* and 0.4386164 at steady, K per W/m. V_core,h/(πa²L) − 1 = −1.25e-05. Robin offset 7.456e-06. z-invariance at (a, ±L/4) holds to 8.5e-05 relative. t_build 4.3 s; t_call 19.6 s transient, 3.0 s steady.
- **Disclosed:**
  1. The ungated t_c read at (3, 0) mm is 1.143e-02, above `INTERIOR_BAND`. The brief gates t\* and steady only, as item 837 did. The shell's mid-radius at t_c is barely heated by a core source (θ 9.2e-03 against 0.19 on the axis), so this is a small-number relative error, printed only.
  2. A dolfinx `cmap` DeprecationWarning, the same call item 837 makes.

**Executes** §10 I71 under the ruling above. Family I71, attempt 1 of 4.

**Code.** New `tests/validation/test_imp8_composite_cylinder_core_source_fem.py`, real, degree 2: import `_composite_cylinder`, `H_DIR`, `CELL_CAP` and item 837's helpers (`_dg0` route, `_theta_pts`, the CN settings and `snapshot_steps`) from `test_imp8_composite_cylinder_fem`; the comparand is item 847's θ. The only physics change from item 837: q = `_dg0(msh, ct, Q_CORE, 0.0)`.

**Gates (imported, unmoved).** (I) |θ_FEM/θ_cf − 1| ≤ `INTERIOR_BAND` at (a, 0) and ((a+b)/2, 0) at t\* and steady, and (I0) on the axis at t\* and steady (the core is the hot side now); (E) ≤ `E_BAND`; (S0) ∫q dx = `Q_CORE`·V_core,h within `B_BAND`, ∫q dx(2) = 0. Preconditions: item 837's (degree-2 coordinates, 53 504 cells, tag-1 facets at r = b, the P2 wiring check).

**Control (asserted, factor 1).** The k_r = `K_GEL` comparand misses θ_FEM(0, steady) by > `INTERIOR_BAND` (predicted ≈ 0.275 from item 847's arithmetic, printed).

**Printed.** The Robin offset (`H_DIR`), V_core,h/(πa²L) − 1, z-invariance at (a, ±L/4), `t_call`.

**Window.** W1 (`solve`, `-n 2`, `timeout -k 30 600`), item 837's form and knob (`FEM_EM_IMP8_I65_PHASE=solve` is item 837's; this module gets its own `FEM_EM_IMP8_I71_PHASE`). Predicted ≈ 40 s (item 837's W1: 37 s).

**Tier.** Heavy label, real, `-n 2`. **Class:** chain step (solve). **Cost:** ≈ 5.5 measured / 20 predicted min.

**Traps.** Item 837's list, plus its W0 finding: classify end faces with `|z_max − z_min| < 1e-3 L`, not 1e-9 (OCC pads bounding boxes); import the builder rather than copying it.

**Scope.** One mesh, one Δt. **Not claimed:** a rate, a finite rod in gel, an implant SAR.

**Negative result.** (I0) red with (I) green: the core's steep profile is under-resolved at h_a = 0.30 mm — report, stop, BLOCKED. No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 14, a 🟡 record line; §10 I71 DONE.

**Depends on:** item 847 landing. If 847 is open or BLOCKED, skip at no cost.

### Review ruling 2026-11-06 03:00 — chain I's last two provisional steps (I73, I74)

The 03:00 review found the queue drained (the 22:30 and 00:00 slots stopped with nothing takeable) and wrote chain I's fifth and sixth provisional steps under rule (6), toward (7e) "implant-driven heating in closed form, then FEM". **Chain I is at 6 of 6** after them; no further provisional chain I step can stand before the 11-07 weekly. No band, tolerance or target is set or moved. **Named in the dashboard for the 11-07 weekly to ratify or revert.**

- **I73 (§9 item 856) — the self-heated rod.** I70 / I71 put the source in the core of an infinite cylinder. On I55's box the source already lies in the rod: `_rod_source` is the DG0 indicator of the cylinder cells (r ≤ `ROD_RADIUS`, |z| ≤ L/2, 6 495 cells). I73 gives exactly those cells (`ROD_K`, `ROD_RHO_C`). This is the finite-rod analogue of I70 / I71 — the implant as both the heat source and a conductor. The conducting body here is the whole cylinder: unlike I66 there is no ball source to keep in gel, so nothing is subtracted.
- **Its anchor is I55's.** The gel arm through the DG0-coefficient path must return I55's five printed θ (`20261104T030134Z_IMP-8-I55-W2-item803.log:142–146`) within one sixth-digit unit, and I55's (I) closed-form gate at `INTERIOR_BAND` is re-asserted on it. The Ti arm is the record. Its negative control is I66's: θ_ti(2, 0) must move off I55's 1.772315e+01 by more than 1e-9 relative.
- **I74 (§9 item 857, `EX-185`, `imp:25`) — an example of I70 / I71**, as I72 was of I64 / I65: the implant as the heat source, closed form beside FEM.

### §9 item 856 — 2026-11-06 (§10 I73, chain I — daily, 2026-11-06 03:00, provisional under rule (6); `IMP-8` step 15, the self-heated rod on I55's box; chain step (solve))

**Outcome (2026-11-06, 04:30 slot, implementer): DONE.** Attempt 1 of 4. `IMP-8` step 15 is a 🟡 record; §10 I73 is DONE.
- **Code:** new `tests/validation/test_imp8_rod_self_heated_fem.py`, real, degree 2. Imports as the brief (`_allsum` not needed and not imported). Body cells are the **positive DG0 dofs of q** from `_source_and_s0` (the simpler route), owned count allreduced. Knob `FEM_EM_IMP8_I73_ARM`. No `src/` change, no landed module edited, not edited after W0 / W1 / W2.
- **Windows:** W0 `20261106T103136Z_IMP-8-item856-W0.log` (`-n 2`, `-k k0`, 1 passed, 5 s). W1 `20261106T103149Z_IMP-8-item856-W1.log` (arm `ti`, `-n 8`, real, 430 s, 5 passed / 2 skipped per rank, `[capture] rc=0`). W2 `20261106T103928Z_IMP-8-item856-W2.log` (arm `gel`, 429 s, 6 passed / 1 skipped per rank, `[capture] rc=0`). `pgrep -c python3` read 0 before each.
- **Gates:** (E) 2.698e-14 (ti) / 2.032e-14 (gel) against 1e-9; (S0) 4.441e-16 on both; cells 485 291; body cells 6 495 on both arms; V_body/(πa²L) 1.004358; non-finite dofs 0.
- **Anchor (arm gel, asserted):** I55's five θ reproduced, |diff| 3.501e-07, 3.421e-07, 4.462e-07, 2.442e-08, 4.345e-06 against sixth-digit units 1e-05, 1e-05, 1e-05, 1e-06, 1e-04. I55's (I) at the four `GATE_PTS`: 1.296e-03, 1.943e-03, 2.867e-03, 3.216e-03 against 1e-2.
- **Control (arm ti, asserted):** θ_ti(2, 0)/θ₅₅ − 1 = −2.0707e-02 against the 1e-9 floor.
- **Record (K/W, t\* = 900 s):** θ_ti = 4.583697, 1.098572, 2.774057, 0.7112581, 17.35615, 10.41661 at (10, 0), (20, 0), (10, 15), (20, 15), (2, 0), (2, 15) mm; ti/I55 = 0.9854495, 0.9937559, 1.0151136, 1.0095081, 0.9792928. New gel point θ(2, 15) = 9.482226. Axial ratio θ(2, 15)/θ(2, 0): gel 0.535019, Ti 0.600168. Ti against QP·Φ (not a gate: no rod in the comparand): 1.3e-02 … 1.0e-01. t_call 408.4 s (ti) / 407.6 s (gel), n_factorizations 2.
- **Prediction (printed, not asserted):** both signs came out as predicted — θ(2, 0) falls and the axial ratio moves toward 1 on the Ti arm. The rod carries heat from mid-length toward its ends: θ(2, 15) is 9.9 % above the gel arm's, and the z = 15 mm gel points rise 1.0–1.5 %.

**Executes** §10 I73 under the ruling above. Family I73, attempt 1 of 4.

**Code.** New `tests/validation/test_imp8_rod_self_heated_fem.py`, real build, degree 2. No `src/` edit and no landed module edited.
- **Imports, unchanged:**
  - from `test_imp8_rod_line_source_fem`: the `box` fixture, `GATE_PTS`, `PRINT_PTS`, `P`, `QP`, `T_INF`, `W0C_CELLS`, `DEGREE`, `_need`, `_theta_pts`, `_source_and_s0`, `_rod_source`;
  - I66's constants and the solve arguments, as item 846's module imports them (`ROD_K`, `ROD_RHO_C` (both `# nominal`), `ROD_RADIUS`, `K_GEL`, `H_GEL`, `N_MAIN`, `RHO_C`, `E_BAND`, `B_BAND`, `ROBIN_TAG`, `_allsum`, `_say`, `INTERIOR_BAND`, `L`, `MM`, `T_PRINT`);
  - `_phi` (`test_imp8_rod_line_source_closed_form`).
  Read every signature before use.
- **Knob** `FEM_EM_IMP8_I73_ARM` = `ti` (default) / `gel`; anything else is refused. The phase knob `FEM_EM_IMP8_I55_PHASE=transient` goes on the command, as I66's did.
- **The body cells** are exactly `_rod_source`'s selection (midpoint r ≤ `ROD_RADIUS` ∧ |z| ≤ L/2, owned + ghost). Rebuild that boolean the same way, or take it from the positive DG0 dofs of the returned q (the simpler route; say which). Print the allreduced owned count, which must be 6 495 as I55's W0c.
- **Coefficients:** as item 846. DG0 `k` = `K_GEL` everywhere, then `ROD_K` on the body cells; DG0 `rho_c` = `RHO_C`, then `ROD_RHO_C`. Arm `gel` makes no overwrite. Call `scatter_forward`, then assert there is no NaN.
- **Solve:** exactly I55's W2 call (`test_imp8_rod_line_source_fem.py:322–325`), with the DG0 `k`, `rho_c` in place of the constants and q from `_source_and_s0(msh)`.
- **Reads:** `_theta_pts` at `GATE_PTS + PRINT_PTS + ((2.0, 15.0),)`, one batch, on both arms.

**Gates (asserted, bands imported).**
- **(E)** `max_balance_residual` ≤ `E_BAND`. **(S0)** |∫q − P|/P ≤ `B_BAND`. **Preconditions:** cells = `W0C_CELLS` (485 291); body cells = 6 495.
- **Arm `gel` only (W2): the anchor and the reproduction control.**
  - Each of I55's five θ is reproduced within one unit of its sixth significant digit, |θ − θ₅₅| ≤ 10^(⌊log₁₀|θ₅₅|⌋ − 5). I55's values (K/W, t\* = 900 s, `20261104T030134Z_IMP-8-I55-W2-item803.log:142–146`) are (10, 0) 4.651377e+00, (20, 0) 1.105475e+00, (10, 15) 2.732755e+00, (20, 15) 7.045591e-01 and (2, 0) 1.772315e+01.
  - I55's (I): |θ/(QP·Φ) − 1| ≤ `INTERIOR_BAND` at the four `GATE_PTS`, with Φ = `_phi(r·MM, z·MM, T_PRINT)`.
- **Arm `ti` only (W1): the negative control.** |θ_ti(2, 0)/1.772315e+01 − 1| > 1e-9. A DG0 path that ignored its coefficients would return I55's θ. I66's caveat applies: the 1e-9 floor sits below what a seven-digit comparand resolves, so the verdict rests on the size of the move, which is printed.

**Printed (🟡 record, no band, ≤ 7 digits).**
- On W1: θ_ti at the six points beside I55's, and θ_ti/θ₅₅ per point. Also θ_ti against QP·Φ, labelled "not a gate: the comparand has no rod".
- The axial ratio θ(2, 15)/θ(2, 0) on both arms. The arm `gel` value at (2, 15) is new (I55 did not read it).
- The body-cell count and volume against πa²L, `t_call` and `n_factorizations`.
- **Predicted, not asserted** (unlogged review reading): a conducting rod evens out θ along its own length. So θ(2, 15)/θ(2, 0) should move toward 1 on the Ti arm, and θ(2, 0) should fall. The sign is a prediction. If it comes out the other way, that is the finding.

**Windows.** Pre-check `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before each. Every window is a bare foreground harness call, Bash timeout 660000 ms. Use the durable capture of `20261104T030134Z_IMP-8-I55-W2-item803.log:12`, including the trailing `; exit $rc`, container `timeout -k 30 600`, `mpiexec -n 8`, and the arm knob added to the env.
- **W0** (`-n 2`, `timeout -k 30 120`): `-k k0` only. This is a knob test that the parser accepts `ti` / `gel` and refuses junk (no mesh).
- **W1** arm `ti`; **W2** arm `gel`. Each should take ≈ 440–480 s (I55's W2: 461 s, `t_call` 439.8 s; I66's W1: 433 s). A kill parks the item. Never raise the ceiling.

**Tier.** Heavy, real, `-n 8`. **Class:** chain step (solve). **Cost:** ≈ 5.5 measured / 25 predicted min (two 8-minute windows).

**Traps already paid for.** Item 846's list:
- real build only, since I55's record is real;
- two heavy windows never at once;
- reduce before asserting, and count verdicts per rank by token;
- one batch per read;
- `-s`.
Also: a leak-guard hit on a printed θ line (item 851's W3) means the log moves unchanged to `logs/leak-held/`, and the slot reports it. The item's verdict is still read from that log, and the gate values are quoted in this file.

**Scope.** One mesh, one Δt, nominal Ti properties, perfect thermal contact, and a uniform source in the rod cells (I55's). **Not claimed:** a Joule-heating SAR shape inside a metal rod, a gate on the rod's effect (no closed form has the finite conducting rod), entry 55, or a Ti property.

**Negative result.**
- (E) or (S0) red: stop, BLOCKED.
- W2's reproduction red, or I55's (I) red on the gel arm: the DG0 path is not I55's. Print every θ, stop, BLOCKED; W1's record is not banked without it.
- W1's control red (θ unmoved): the coefficient never reached the form. BLOCKED.
- No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 15, a 🟡 record line: the first number for what the rod's own conduction does when the rod is the heat source. §10 I73 DONE.

**Depends on:** nothing (I55's record and I66's module are on `main`).

## §7 row moved verbatim 2026-11-07 (rotate_plan_archive.py chunks)

 The heated rod in gel: the finite line source in closed form (Phase 7; opened 2026-11-02 10:30 daily review, rule (5) lane 3; the weekly ratifies) | ✅ on the closed form (ratified by the 11-02 weekly); the FEM half landed as 🟡 records — I55, I56, I60, I62, I63 and the lumped reading I58 — I61 is sized as `EX-182` (`imp:23`, ≤ 120 000 cells, §9 item 829), and the rod's own thermal properties are posed as steps 9–10: I64 the composite cylinder in closed form (step 9, 🟡 record, item 843 landed 2026-11-05: (St) / (En) / (Lim) / (Bq) ≤ 3.6e-16 against `B1_BAND`), I65 its FEM twin (step 10, 🟡 record, item 837 landed 2026-11-05: (I) ≤ 6.0e-4 and (I0) 1.66e-3 against `INTERIOR_BAND` — the first FEM gate with an implant's thermal properties in it), then I66 I56's rod as a conducting body (step 11, 🟡 record, item 846 landed 2026-11-05: the Ti rod moves θ(2, 0) mm by +55.5 % and the tip-point ratio 9.879 → 5.333 against I56, the gel arm reproducing I56's six θ, (E) ≤ 1.9e-14 — the first number for what the rod's own conduction does to the tip temperature), and I70 the composite cylinder with the source in the core (step 12, 🟡 record, item 847 landed 2026-11-05: (St) / (En) / (Lim) / (Bq) ≤ 4.5e-14 against `B1_BAND`; steady θ(0) 0.4388 (Ti) / 0.5595 (gel core) K per W/m — the implant as the heat source), and I69 the near-adiabatic rod (step 13, 🟡 record, item 851 landed 2026-11-06: rod cells at 1e-6 × gel k and ρc move θ(2, 0) mm by −8.46 % and the tip-point ratio 9.879 → 11.61 against I56, (E) 2.0e-14 — the bracket's third arm; near the rod at z = 0 the FEM orders Ti > gel > hole, opposite to I64's shell-source cylinder at r = a), and I71 the core-source cylinder in FEM (step 14, 🟡 record, item 852 landed 2026-11-06: (I) ≤ 1.7e-3 and (I0) ≤ 5.5e-4 at t\* and steady against item 847's series at `INTERIOR_BAND` — the implant as the heat source, closed form to FEM), and I73 the self-heated rod on I55's box (step 15, 🟡 record, item 856 landed 2026-11-06: the gel arm through DG0 coefficients reproduces I55's five θ and its (I) ≤ 3.2e-3; the Ti source cells move θ(2, 0) mm by −2.07 % and the axial ratio θ(2, 15)/θ(2, 0) 0.535 → 0.600, (E) ≤ 2.7e-14 — what the rod's own conduction does when the rod is the heat source). *History: `docs/planning/chunks/IMP-8.md`.* 

### §9 item 874 — 2026-11-08 (§10 I78, chain I — IMP-8 step 16, the lumped F2182-shaped prediction; lane 3 → chain step (closed form), smoke, numpy)

**Outcome (07:30, 2026-11-08 slot, implementer): DONE.** Attempt 1 of 4. `IMP-8` step 16 is a 🟡 record; §10 I78 is DONE. New `tests/validation/test_imp8_probe_f2182_closed_form.py` (imports as the brief, no `src/` or landed-module edit, not edited after W1). W1 `20261108T133157Z_IMP-8-item874-W1.log` (`-n 1`, real, 5 passed, 5 s harness, `[capture] rc=0`; `pgrep -c python3` read 0 before). (Lim) 1.110e-16, (Bq) 4.030e-14 (t = 900 s) / 2.520e-14 (t = ∞), (Dut) 2.486e-16 line / 5.839e-17 twin, all ≤ band. Control 1 − θ(Z_C)/θ(Z_P) at 900 s: line 0.2941583, twin 0.5253344 (steady-form prediction 0.1926 / 0.4854, printed). Record, the F2182-shaped probe (placement and time quoted from memory; W3 confirms or corrects), K per watt per rod at (0, ±16 mm): line θ(900 s) 9.648662 continuous, 1.978384 at D = 0.5, θ_ss 15.18156; twin 63.14 / 2.498341 / 68.45374; twin/line 6.543912 continuous, 1.262819 at D = 0.5; ±Z_P symmetry 0 exactly.

**Executes** §10 I78, the 11-07 weekly's enumerated step (PROJECT_PLAN §10 chain I). The target, anchor and control are the weekly's; what follows is execution detail. Family I78, attempt 1 of 4. **Executor:** the implementer.

**Code.** New `tests/validation/test_imp8_probe_f2182_closed_form.py`, numpy / scipy, the image's default (real) build. No `src/` edit and no landed module edited.
- **Imports, unchanged** (read every signature first):
  - step 0 (`test_imp8_rod_line_source_closed_form`): `_phi`, `_phi_tips`, `_point`, `KAPPA`, `L`, `MM`, `T_PRINT`, `QUAD`;
  - I58 (`test_imp8_entry48_lumped`, item 786): `_tips_steady`, `_gl01`, `BQ_BAND`;
  - I59 (`test_imp8_rod_line_source_duty_closed_form`, item 787): `_phi_dot`, `PERIOD`, `SCALE`, `T_STAR`;
  - I62 (`test_imp8_rod_tips_duty_closed_form`, item 811): `_duhamel`, `_convolution`, `_phi_tips_dot`. These are item 787's `MAT-17` Duhamel rule made kernel-generic: `_duhamel(kernel, r, z, t, duty)` is the step sum Σⱼ ±θ_step(t − tⱼ) over `_edges(duty)`, and `_convolution(kernel_dot, r, z, t, duty)` is the from-rest convolution;
  - `B1_BAND` (`test_imp4_pec_cylinder_closed_form`, 1e-12) and `INTERIOR_BAND` (`test_mat13_transient_perfusion_sphere`, 1e-2), by name.
- **The probe:** on the axis, 1 mm beyond each tip, (r, z) = (0, ±`Z_P`) with `Z_P` = L/2 + 1 mm = 16 mm. The control point is `Z_C` = L/2 + 2 mm. Both are module constants.
- **The duty form** is D = 0.5 at `PERIOD` = 450 s (on [0, 225) s, off [225, 450) s, and so on). That is the period I59–I63 share. I56 itself has no period. The weekly's "period as I56's" is read as the two-tip duty runs' period (I62 / I63 import the same `PERIOD`), and the brief says so here. **Read at t = 900 s, the end of an off phase:** the edge at 900 s is excluded by both routes (`t > te`, `j·PERIOD < t`).
- **Per watt:** 1 W per rod, through `SCALE` = 1/(4π`K_GEL`L). Both shapes carry the same total power.

**Gates (asserted, bands imported).**
- **(Lim)** t → ∞ equals the steady line-source form. `_phi(0.0, Z_P, np.inf)` (erfc → 1, a `quad` of 1/|z − z′|) against the on-axis steady form ln((Z_P + L/2)/(Z_P − L/2)), relative ≤ `B1_BAND`. This axis form is the r → 0 limit of `_phi_steady` and is written in the module, because `_phi_steady` and I58's `_line(r, z, np.inf)` divide by r (see Traps). Printed beside, not asserted: `_phi_tips(0.0, Z_P, np.inf)` against `_tips_steady(0.0, Z_P)`, which is exact by construction; and `_phi_steady(1e-9, Z_P)` against the axis form, whose expected O(r²/(Z_P − L/2)²) ≈ 1e-12 sits at the band, so it is printed only.
- **(Bq)** N vs 2N on the line kernel. A Gauss–Legendre rule in z′ on [−L/2, L/2] (I58's `_gl01(n)`, mapped) applied to `_point(|Z_P − z′|, t)` at n = 64 and 2n = 128, at t = `T_PRINT` and t = ∞: |GL_n/GL_2n − 1| ≤ `BQ_BAND` (I58's (Bq) band, imported). |GL_2n/`_phi` − 1| is printed beside. Predicted round-off, unlogged review arithmetic: the integrand's nearest singularity is 1 mm past the endpoint, Bernstein ρ ≈ 1.44, so ρ^(−2n) ≈ 1e-20 at n = 64.
- **(Dut)** the duty-cycled read at t = 900 s equals Duhamel's form. The step sum `_duhamel(_phi, 0.0, Z_P, T_STAR, 0.5)` against `_convolution(_phi_dot, 0.0, Z_P, T_STAR, 0.5)`, and the same pair for the twin (`_phi_tips` / `_phi_tips_dot`). Each |difference| is normalised by θ_ss at the probe, as I59's (Duh) normalises (`test_imp8_rod_line_source_duty_closed_form.py:162`), and must be ≤ `B1_BAND`. Precedents at the looser `DUH_BAND` (1e-10): 2.952e-16 (`20261103T211919Z_IMP-8-item787-W1.log:60`) and 2.951e-16 (`20261104T103134Z_IMP-8-item811-W1.log:61`).

**Control (the weekly's, asserted at factor 1 on `INTERIOR_BAND`).** The probe moved to 2 mm reads lower: 1 − θ(0, Z_C)/θ(0, Z_P) > `INTERIOR_BAND` at t = 900 s under the continuous source, for both shapes, with θ(Z_C) < θ(Z_P). The ceiling is 1, since θ > 0. Predicted from the steady forms (unlogged review arithmetic, printed beside the measurement): line 1 − ln 16/ln 31 ≈ 0.19, twin ≈ 0.49. A flat readout would mean the probe cannot see its own placement.

**Printed (🟡 record, no band, ≤ 7 digits, K per watt per rod).** At (0, ±Z_P), for the line and the twin:
- θ(900 s) continuous, θ(900 s) at D = 0.5 and θ_ss;
- duty/continuous and twin/line;
- the same reads at Z_C;
- |θ(+Z_P) − θ(−Z_P)|, which is round-off by symmetry.
Label the block "the F2182-shaped probe (placement and time quoted from memory; W3 confirms or corrects)". These are the numbers I75 and I76 are read beside.

**Window.** W1, smoke, `-n 1`, real, container `timeout -k 30 60`, label `IMP-8-item874-W1`. Make one bare foreground harness call with Bash timeout 660000 ms, using the durable capture of `20261106T103149Z_IMP-8-item856-W1.log:12` verbatim including the trailing `; exit $rc`, with this module, `-n 1`, no knob env and raw log `/workspace/logs/imp8-item874-w1-raw.log`. Pre-check `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Predicted < 15 s: I58 took 9 s (`20261103T211525Z_IMP-8-item786-W1.log:91`), I59 3 s (`20261103T211919Z_IMP-8-item787-W1.log:89`) and I62 6 s (`20261104T103134Z_IMP-8-item811-W1.log:90`).

**Tier.** Smoke, numpy, `-n 1`. **Class:** lane 3 → chain step (closed form). **Cost:** ≈ 3.6 measured / 15 predicted min.

**Traps already paid for.**
- **`_phi_steady` and I58's `_line(·, ·, np.inf)` divide by r.** At r = 0 that is a `ZeroDivisionError` on Python floats, or inf − inf = nan on numpy. Never call them on the axis; use the axis form above.
- At Z_P > L/2, `_phi` takes one piece (no split at z′ = z). The integrand peaks at the endpoint z′ = L/2, 1 mm from the probe. Use the `QUAD` settings unchanged. I58 disclosed a `quad` `IntegrationWarning` near the axis, and (Bq) is what shows whether it reaches the reading.
- `_step` guards t ≤ 0. Never evaluate a kernel at age 0.
- Use `-s`, or the prints are lost. Print ≤ 7 digits.
- **The leak guard:** if a printed θ line is refused by the hook, the log moves unchanged to `logs/leak-held/` and the slot reports it, as for item 851's W3. Never `--no-verify`.

**Scope.** An infinite gel medium, two idealised source shapes (a zero-radius line and two point tips), 1 W per rod, one probe placement. **Not claimed:** the box walls, the rod's own conduction (I73 / I75), a ball-averaged read, any F2182 figure or published value. The weekly's caveat stands: the standard's probe placement (1 mm beyond the tip, on the axis) and time (900 s) are quoted from memory and marked so (§10 I32b), and W3 confirms or corrects them. The twin's point source sits 1 mm from the probe, so it reads a point's 1/R there. The FEM's 1 mm balls (I56) put their surface at the probe, so the twin is not the FEM's source shape at that distance. That is printed, never compared.

**Negative result.**
- (Lim), (Bq) or (Dut) red is a derivation or quadrature error: stop, BLOCKED, band unmoved, every value printed.
- A control under the band contradicts the steady arithmetic above: report the measurement and stop.
- No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 16, a 🟡 record line: the lumped prediction at the F2182-shaped probe. §10 I78 is DONE, and I75 / I76 have their closed form beside them.

**Depends on:** nothing (items 786, 787 and 811 are on `main`).

### §9 item 881 — 2026-11-08 (§10 I75, chain I — IMP-8 step 17, the F2182-shaped readout on the bracket fixture; chain step (solve), heavy, real, -n 8)

**Executes** §10 I75, the 11-07 weekly's enumerated step. The target, anchor and control are the weekly's, and what follows is execution detail. Family I75, attempt 1 of 4. **Executor:** the implementer.

**The control as the weekly wrote it is degenerate, so it is re-posed here for the review to ratify.** The step says "the probe read at t = 900 s exceeds the read at t\* on every arm (monotone heating under a continuous source; asserted at factor 1)". On I73's fixture t\* *is* 900 s: `T_PRINT` = 900.0 (`test_imp8_rod_line_source_closed_form.py:48`), I73's CN run is `N_MAIN` steps of `T_PRINT / N_MAIN` (`test_imp8_rod_self_heated_fem.py:131–136`), and its log prints `t*=900 s` (`20261106T103149Z_IMP-8-item856-W1.log:113`). The weekly's own words are "monotone heating under a continuous source", so this brief asserts that monotonicity against the step-30 snapshot, t₂ = 675 s (I60's `T2` / `STEP_T2`, `test_imp8_rod_line_source_duty_fem.py:56–57`):

θ_probe(900 s) / θ_probe(675 s) > 1 (factor 1) at both probes on every arm, under the continuous source.

The size of the rise is predicted, not asserted: per-cent level (unlogged review arithmetic, 1 mm from a heated end). **Flagged in the slot report for the review.**

**The duty read is a superposition, not a second solve, and that choice is flagged for the review.** "Under I56's duty form" is read as D = 0.5 at `PERIOD` = 450 s, the two-tip duty runs' period (I56 has none of its own; I62 / I63 import I59's `PERIOD`). It is taken from the continuous run by linearity: the problem is linear and time-invariant from rest, with θ₀ = `T_INF` = 0. The read is θ_D(900) = θ(900) − θ(675) + θ(450) − θ(225), from `snapshot_steps=(10, 20, 30)`.

This costs no extra solve. It is not bit-identical to a switched solve: each shifted step response restarts with the two backward-Euler half-steps, while I60's direct run goes through the switches in CN with no restart (I60's disclosure, item 799). The gel arm prints the difference. **If the review reads the step as a direct switched solve**, that is three more windows of ≈ 7.3 min each (I60's D = 0.5 window: 418 s, `20261104T044429Z_IMP-8-I60-W1-item799.log:181`). Together with W1–W3 that does not fit one slot, so it is a second item, not this one.

**Code.** New `tests/validation/test_imp8_rod_probe_f2182_fem.py`, real build, degree 2. No `src/` edit and no landed module edited.
- **Imports, unchanged** (read every signature first):
  - I73 (`test_imp8_rod_self_heated_fem`): `BODY_CELLS`, `I55_THETA`, `READ_PTS`, `CONTROL_FLOOR`;
  - I55 (`test_imp8_rod_line_source_fem`): the `box` fixture, `GATE_PTS`, `QP`, `T_INF`, `W0C_CELLS`, `DEGREE`, `_need`, `_theta_pts`, `_source_and_s0`;
  - I66 / I69 (`test_imp8_rod_conducting_body_fem`): `_coefficients(msh, rod_idx, arm)`, which handles `ti` / `gel` / `hole` with `EPS_HOLE`, and `EPS_HOLE` itself, printed;
  - `ROD_K`, `ROD_RHO_C` (both `# nominal`) and `ROD_RADIUS` as I73 imports them;
  - step 0's `L`, `MM`, `T_PRINT`, `_phi` and I59's `_duhamel` (`test_imp8_rod_line_source_duty_closed_form`; the line's, `_duhamel(r, z, t, duty)`), for the printed closed forms;
  - I60's `T2`, `STEP_T2`;
  - `K_GEL`, `H_GEL`, `N_MAIN`, `RHO_C`, `E_BAND`, `B_BAND`, `ROBIN_TAG`, `_say`, `INTERIOR_BAND` from the modules I73 takes them from.
- **Knob** `FEM_EM_IMP8_I75_ARM` = `ti` (default) / `gel` / `hole`. Anything else is refused. The phase knob `FEM_EM_IMP8_I55_PHASE=transient` goes on the command, as I73's did.
- **Body cells are I73's**, the positive DG0 dofs of q from `_source_and_s0(msh)`, converted to **cell** indices through the DG0 dofmap (owned + ghost). Then call I66's `_coefficients(msh, body_idx, arm)`. Print the allreduced owned count, which must equal `BODY_CELLS` (6 495). **Do not call I73's `_coefficients` with `hole`:** it overwrites only on `ti`, so any other arm silently returns gel (`test_imp8_rod_self_heated_fem.py:109–111`).
- **Solve:** I73's call exactly (`test_imp8_rod_self_heated_fem.py:133–136`) plus `snapshot_steps=(10, 20, 30)`. Assert that snapshot step j is at t = j·`T_PRINT`/`N_MAIN` (225, 450, 675 s).
- **Reads**, one batch per call (the batch rule):
  1. `_theta_pts(res.temperature, READ_PTS)`, I73's batch unchanged, so the reproductions compare call for call;
  2. `_theta_pts(·, PROBE_PTS)` with `PROBE_PTS` = ((0.0, 16.0), (0.0, −16.0)) mm on the final field and on each snapshot;
  3. on arm `gel` only, `_theta_pts(·, GATE_PTS)` on each snapshot, for the printed duty comparison.

**Gates (asserted, bands imported).**
- **Every arm:**
  - (E) `max_balance_residual` ≤ `E_BAND`;
  - (S0) |∫q − P|/P ≤ `B_BAND`;
  - cells = `W0C_CELLS` (485 291), body cells = `BODY_CELLS`, no non-finite coefficient dofs.
- **Arm `gel` — I73's anchor and its (I).**
  - Each of `I55_THETA`'s five θ is reproduced within one unit of its sixth significant digit.
  - The new gel point is reproduced the same way: θ(2, 15) = 9.482226e+00 (`20261106T103928Z_IMP-8-item856-W2.log:118`).
  - I55's (I) holds: |θ/(QP·`_phi`) − 1| ≤ `INTERIOR_BAND` at the four `GATE_PTS`.
- **Arm `ti` — the route control.** I73's six Ti θ are reproduced within one sixth-digit unit through I66's `_coefficients` route: 4.583697e+00, 1.098572e+00, 2.774057e+00, 7.112581e-01, 1.735615e+01, 1.041661e+01 at `READ_PTS` (`20261106T103149Z_IMP-8-item856-W1.log:113–118`). A route that missed the body cells would return I55's θ.
- **Arm `hole` — the coefficient control.** |θ_hole(2, 0)/1.772315e+01 − 1| > `CONTROL_FLOOR`, as item 851's on I66's box. A path that ignored ε would return the gel arm's θ.
- **Every arm — the re-posed weekly control:** θ_probe(900 s) > θ_probe(675 s) at both probes (factor 1, strict).

**Printed (🟡 record, no band, ≤ 7 digits, K/W, 1 W in the rod).**
- **Per arm:**
  - θ_probe(900 s) at +16 and −16 mm, their mean and |difference| (mesh asymmetry);
  - the superposed D = 0.5 θ_probe(900 s);
  - θ_probe at 225 / 450 / 675 s and the ratio θ(900)/θ(675);
  - `t_call`, `n_factorizations` and t_build.
- **Beside them, labelled "infinite medium, line source, no rod":** `QP·_phi(0.0, ±Z_P, T_PRINT)` and the line's `_duhamel(0.0, Z_P, T_PRINT, 0.5)`. These are the numbers I78 gates (item 874). They are computed inline here, so this item does not depend on 874 landing.
- **Arm `gel`:** the superposed duty θ at `GATE_PTS` beside I59's `_duhamel` and beside I60's direct FEM reads 2.240893e+00, 6.683594e-01, 1.366131e+00, 4.289727e-01 (`20261104T044429Z_IMP-8-I60-W1-item799.log:112–115`). It also prints the restart-vs-no-restart difference. **Rule (e):** printed, never asserted.
- **Across the three windows (slot report):** the probe ordering of Ti / gel / hole. Item 856 measured a −2.07 % move at (2, 0) mm and the axial ratio 0.535 → 0.600 toward the tips (`…item856-W1.log:119–120`, `…item856-W2.log:119`), so the Ti probe is *predicted* to read above the gel's. That is a sign only, printed; if it comes out the other way, that is the finding.

**Windows.** Pre-check `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before each.
- Every window is a bare foreground harness call, Bash timeout 660000 ms, using the durable capture of `20261106T103149Z_IMP-8-item856-W1.log:12` verbatim including the trailing `; exit $rc`.
- Use this module, raw log `/workspace/logs/imp8-item881-w<N>-raw.log`, env `FEM_EM_IMP8_I55_PHASE=transient FEM_EM_IMP8_I75_ARM=<arm>`, and container `timeout -k 30 600`.
- **Each arm is its own window**, since one CN run on this box is ≈ 408 s of `t_call` (`…item856-W1.log:112`, `…item856-W2.log:112`).

The windows:
- **W0** `IMP-8-item881-W0` (`-n 2`, `timeout -k 30 120`): `-k k0` only. This is the knob test that `ti` / `gel` / `hole` parse and junk is refused (no mesh).
- **W1** `IMP-8-item881-W1`, arm `ti`, `-n 8`. Predicted ≈ 430–445 s; I73's W1 took 430 s (`…item856-W1.log:190`). The three snapshot copies add only copies.
- **W2** `IMP-8-item881-W2`, arm `gel`, `-n 8`. Predicted ≈ 430–445 s; I73's W2 took 429 s (`20261106T103928Z_IMP-8-item856-W2.log:190`).
- **W3** `IMP-8-item881-W3`, arm `hole`, `-n 8`. Predicted ≈ 430–450 s; item 851's hole window took 438 s per its Outcome above (that log is leak-held, so none of its digits is quoted here).

A kill parks the item. Never raise the ceiling, and never run two heavy windows at once.

**Tier.** Heavy, real, `-n 8`. **Class:** chain step (solve). **Cost:** ≈ 6.8 measured / 35 predicted min (W0 plus three ≈ 7.3-minute windows).

**Traps already paid for.**
- **Items 846 / 856's list:**
  - real build only, since I55's and I73's records are real and every reproduction here compares real-build digits;
  - reduce every rank-local count and integral before asserting, and count verdicts per rank by token;
  - one batch per read, with I73's `READ_PTS` batch unchanged;
  - `-s` on pytest.
- **The hole arm with the source in the hole is new** (I69's source sat in gel). With k = ε·`K_GEL`, the rod core reaches ≈ q′/(4π·ε·`K_GEL`) ≈ 4e6 K/W within a few seconds (κ is unchanged by ε; unlogged review arithmetic). **Print no rod-interior or max θ.** The gel reads are what is recorded. If (E) comes back above `E_BAND`, that is the finding (item 851's rule): report the residual and stop. A review may re-pose ε, never the band.
- **`_phi_steady` divides by r.** The printed closed forms on the axis use `_phi` at finite t only.
- **Leak guard:** a hook refusal on a printed θ line moves the log unchanged to `logs/leak-held/` and the slot reports it (item 851's W3). The gate values are quoted in this file. Never `--no-verify`.

**Scope.** One mesh, one Δt, nominal Ti, perfect thermal contact, and I55's uniform source in the rod cells. The probe placement and time are the weekly's, quoted from memory of the standard and marked so (§10 I32b); W3 confirms or corrects them. **Not claimed:**
- a published F2182 value or an F2182 measurement (W3);
- a direct switched-source duty solve (see above);
- the clinical probe (entry 63 / I76);
- a Joule SAR shape in a metal rod;
- a gate on the rod's effect at the probe (no closed form has the finite conducting rod).

**Negative result.**
- (E), (S0) or a precondition red: stop, BLOCKED.
- The gel arm's anchor or (I) red, or the Ti arm's reproduction red: the body route is not I73's. Print every θ, stop, BLOCKED; no arm's record is banked without both.
- The hole control red (θ unmoved): ε never reached the form. BLOCKED.
- The re-posed control red, meaning θ falls under a continuous source from rest: that contradicts the physics of a heating from rest. Print every snapshot read, stop, BLOCKED, and loosen nothing.
- No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 17, a 🟡 record line: the F2182-shaped probe on the gel / Ti / hole bracket, continuous and D = 0.5, per watt. §10 I75 is DONE, and I76 gets its bracket arms.

**Depends on:** nothing. I73's module and logs, I66 / I69's module and I60's log are on `main`; I78 (item 874) is computed inline, not imported.

**Outcome (2026-11-08 13:30 slot, take-next second item, implementer, not delegated; attempt 1 of 4): DONE.** New `tests/validation/test_imp8_rod_probe_f2182_fem.py`; no `src/` or landed module edited. `pgrep -c python3` read 0 before each window.
- **W0** `20261108T194310Z_IMP-8-item881-W0.log` (`-k k0`, `-n 2`): 1 passed per rank, 4 s.
- **W1** `20261108T194323Z_IMP-8-item881-W1.log` (arm `ti`, `-n 8`): 6 passed / 3 skipped per rank, 470 s (t_call 449.6 s), Status 0.
  - The six Ti θ are reproduced, with the largest |diff| 2.619e-06 at (2, 0) against the 1e-04 unit (`:113–118`).
  - (E) 2.509e-14, (S0) 4.441e-16, body cells 6495 (`:111`, `:131`).
  - Probe +16 / −16 mm: 1.064117e+01 / 1.049960e+01 K/W; θ(900)/θ(675) 1.059614 / 1.060263; D = 0.5 superposed 2.109553e+00 / 2.099602e+00 (`:120–122`).
- **W2** `20261108T195128Z_IMP-8-item881-W2.log` (arm `gel`): crashed after the solve, before any gate value. `dict(I55_THETA, **GEL_NEW)` raised `TypeError` (tuple keys), in non-gate code. Fixed to `{**I55_THETA, **GEL_NEW}` and re-run once at the same `-t` (the standing ruling); first log kept.
- **W2b** `20261108T195853Z_IMP-8-item881-W2b.log` (arm `gel`): 7 passed / 2 skipped per rank, 441 s, Status 0.
  - The six gel θ are reproduced, including (2, 15) 9.482226e+00 (`:113–118`).
  - (I) worst 2.867e-03 at `GATE_PTS`; (E) 1.898e-14 (`:129`).
  - Probe 9.220827e+00 / 9.102961e+00; θ(900)/θ(675) 1.066831 / 1.067542; D = 0.5 superposed 1.979610e+00 / 1.970586e+00 (`:120–121`).
  - Gel duty at `GATE_PTS`, superposed vs I60's direct run (restart vs no restart): rel +5.255e-04 / −1.045e-03 / +2.824e-04 / −8.501e-04. Against I59's `_duhamel`: rel +1.554e-04 … +1.899e-03, printed (`:125–128`).
- **W3** `20261108T200636Z_IMP-8-item881-W3.log` (arm `hole`): 6 passed / 3 skipped per rank, 459 s, Status 0.
  - (2, 0) move 4.3944e-03 > `CONTROL_FLOOR` (`:119`); (E) 2.920e-14 (`:125`).
  - Probe 8.815610e+00 / 8.750725e+00; θ(900)/θ(675) 1.068047 / 1.068419; D = 0.5 superposed 1.895292e+00 / 1.888582e+00 (`:120–121`). No rod-interior θ printed.
- **Bracket (slot report).** Probe mean per watt at 900 s: Ti 10.57 > gel 9.162 > hole 8.783 K. The predicted sign (Ti above gel) holds. The line closed form at the probe is 9.648662 K continuous and 1.978384 K at D = 0.5. Mesh asymmetry between the two probes is 1.4 / 1.2 / 0.65 % of the mean.
- **Disclosures for the review:** (1) the weekly's control is re-posed as θ(900) > θ(675) (this brief); (2) D = 0.5 is read by superposition, not a switched solve; (3) rule (i): W1 ran on the module before the W2 fix. The fix touches only the `gel` branch of a conditional and the gel-only test, so the `ti` path is byte-equivalent, but W1 was not re-run.
- **Status moved:** `IMP-8` step 17 is a 🟡 record and §10 I75 is DONE; I76 gets its bracket arms.

### §9 item 894 — 2026-11-08 (§10 I80, chain I — daily, 2026-11-08, provisional, rule (6); IMP-8 step 18, the lumped probe-placement ladder; lane-3 closed form → chain step, smoke, numpy)

**Outcome (22:30, 2026-11-08 slot, implementer): DONE.** Attempt 1 of 4. `IMP-8` step 18 is a 🟡 record; §10 I80 is DONE. New `tests/validation/test_imp8_probe_ladder_closed_form.py` (imports as the brief, no `src/` or landed-module edit, not edited after W1). W1 `20261109T043829Z_IMP-8-item894-W1.log` (`-n 2`, real, 12 passed per rank — 7 new + I78's 5 — 4 s harness, pytest 2.46 s; `pgrep -c python3` read 0 before). Worst per rung (16 / 17 / 18 / 20 mm): (Lim) 1.110e-16 / 0 / 2.220e-16 / 2.220e-16; (Bq) 4.030e-14 / 2.554e-14 / 1.965e-14 / 1.499e-14 (falling outward, as predicted); (Dut) 2.486e-16 / 1.387e-16 / 7.366e-17 / 2.051e-16. (Rep) all 12 `.7g` strings match I78's record at rungs 1–2. Control drops θ(900 s) / θ_ss: line 0.2941583 / 0.1926037, 0.2302515 / 0.1351421, 0.3462264 / 0.1884924; twin 0.5253344 / 0.4853516, 0.3689553 / 0.315508, 0.4670608 / 0.3714286 — the θ_ss drops equal the review's predictions to its four digits. Record, K per watt per rod, continuous θ(900 s) line 9.648662 / 6.810428 / 5.242317 / 3.427288, twin 63.14 / 29.97039 / 18.91265 / 10.07929; D = 0.5 line 1.978384 / 1.822237 / 1.668552 / 1.374611; ±z symmetry 0 exactly at every rung.

*Written by the 2026-11-08 18:00 daily review under §9 rule (6). It works toward the chain's stated readout: I76 sets the clinical probe "beside I75's bracket arms and I78's closed form" (`PROJECT_PLAN.md:5068–5071`). Target, bands and controls are imported, and none is set. Family I80, attempt 1 of 4. **Executor:** the implementer. Depends on nothing.*

**Why.** The F2182 placement is quoted from memory (W3 confirms or corrects), and the two probes in this chain are not the same point relative to the source:
- I75's bracket probe is 1 mm beyond a flat-ended source, at z = 16 mm = I78's `Z_P`;
- entry 63's clinical probe is 1 mm beyond a capsule apex, at z = `ROD_LENGTH`/2 + `ROD_RADIUS` + 1 mm = 17 mm = I78's `Z_C` (`test_geo36_rod_f_clinical_drive.py:288–289`, `test_geo36_rod_f_clinical_heating.py:293`, `:582–584`; item 893's (P)).

I78 gated the closed form at 16 mm only. Its 17 mm read is a control and a print. A gated closed-form ladder of placements makes the lumped column of I76's table hold at the clinical placement, and keeps it valid if W3 moves the offset to any rung.

**Code.** New `tests/validation/test_imp8_probe_ladder_closed_form.py`, numpy / scipy, the image's default (real) build. No `src/` edit and no landed module edited.
- **Imports, unchanged** (read every signature first):
  - I78 (`test_imp8_probe_f2182_closed_form`): `SHAPES`, `_axis_steady`, `_gl_line`, `_reads`, `Z_P`, `Z_C`, `DUTY`, `N_GL`, `_p`;
  - step 0 (`test_imp8_rod_line_source_closed_form`): `L`, `MM`, `T_PRINT`, `_phi`;
  - I59 (`test_imp8_rod_line_source_duty_closed_form`): `SCALE`, `T_STAR`;
  - I62 (`test_imp8_rod_tips_duty_closed_form`): `_duhamel`, `_convolution`, which are the kernel-generic forms I78 uses;
  - `B1_BAND` (`test_imp4_pec_cylinder_closed_form`), `BQ_BAND` (`test_imp8_entry48_lumped`) and `INTERIOR_BAND` (`test_mat13_transient_perfusion_sphere`), by name.
- **The ladder.** `LADDER_MM = (1.0, 2.0, 3.0, 5.0)` mm beyond L/2, so z_k = 16, 17, 18, 20 mm on the axis, both signs. Rung 1 is I75's probe and rung 2 the clinical one; label each on its print.

**Gates (asserted, bands imported).** At every rung:
- **(Lim)** |`_phi(0, z_k, ∞)` / `_axis_steady(z_k)` − 1| ≤ `B1_BAND`. Unlogged review arithmetic: the steady forms are ln 31, ln 16, ln 11, ln 7.
- **(Bq)** |`_gl_line(z_k, t, N_GL)` / `_gl_line(z_k, t, 2·N_GL)` − 1| ≤ `BQ_BAND` at t = `T_PRINT` and ∞. Rung 1 is I78's worst: 4.030e-14 (`20261108T133157Z_IMP-8-item874-W1.log:47`). Farther rungs move the pole away from the endpoint, so they converge faster.
- **(Dut)** For line and twin, |`_duhamel(kern, 0, z_k, T_STAR, DUTY)` − `_convolution(kern_dot, 0, z_k, T_STAR, DUTY)`| / (`SCALE`·steady(z_k)) ≤ `B1_BAND`. I78 measured 2.486e-16 / 5.839e-17 at rung 1 (`:50–51`).
- **(Rep)** At rungs 1 and 2, `_reads` reproduces I78's printed record by `.7g` string (`:58–59`, `:61–62`):
  - line: `9.648662` / `1.978384` / `15.18156` and `6.810428` / `1.822237` / `12.25753`;
  - twin: `63.14` / `2.498341` / `68.45374` and `29.97039` / `2.451249` / `35.22961`.

**Control (I78's, extended along the ladder; asserted at factor 1 on `INTERIOR_BAND`).** For both shapes, continuous θ(900 s) and θ_ss strictly decrease from rung to rung, and each successive drop 1 − θ(z_{k+1})/θ(z_k) > `INTERIOR_BAND`. The ceiling is 1 (θ > 0).
- Predicted from the steady forms (unlogged review arithmetic, printed beside): line 0.1926 / 0.1352 / 0.1885 and twin 0.4854 / 0.3155 / 0.3714. The 18 → 20 mm step is 2 mm.
- At 900 s I78 measured the first drop at 0.2941583 (line) and 0.5253344 (twin) (`:54–55`).
- A flat ladder would mean the probe cannot see its own placement.
- **The D = 0.5 reads are not asserted monotone.** At 900 s the source has been off for 225 s and the heat is spreading outward, so they are printed only.

**Printed (🟡 record, no band, ≤ 7 digits, K per watt per rod, infinite gel).** For each rung and shape:
- θ(900 s) continuous, θ(900 s) at D = 0.5, θ_ss, duty/continuous;
- θ(z_k)/θ(16 mm);
- twin/line;
- |θ(+z) − θ(−z)|.

Label rung 1 "the bracket's probe (I75: 1 mm beyond a flat end)" and rung 2 "the clinical probe (entry 63: 1 mm beyond the capsule apex, 2 mm beyond the line's end)". Head the block "placement and time quoted from memory; W3 confirms or corrects".

**Window.** Before it, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. **W1** `IMP-8-item894-W1`, smoke, real, `-n 2`. Make one bare foreground harness call, Bash timeout 660000 ms; I78's module rides along as its own regression:

`scripts/testing/run_and_log.sh IMP-8-item894-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/validation/test_imp8_probe_ladder_closed_form.py tests/validation/test_imp8_probe_f2182_closed_form.py -v -s --tb=short'"`

**Price.** I78's module took 5 s at `-n 1` (`20261108T133157Z_IMP-8-item874-W1.log:74`; pytest 2.73 s, `:69`) for two placements. Four rungs, plus I78 re-run, predict ≈ 10–20 s.

**Tier.** Smoke, numpy, `-n 2`. **Class:** lane-3 closed form → chain step. **Cost:** ≈ 4.1 measured / 15 predicted min.

**Traps already paid for** (item 874's list).
- `_phi_steady` and I58's `_line(·, ·, ∞)` divide by r. Never call them on the axis; use `_axis_steady`.
- At |z| > L/2, `_phi` takes one piece. Use the `QUAD` settings unchanged. A `quad` `IntegrationWarning` is not a red; (Bq) shows whether it reaches the reading.
- `_step` guards t ≤ 0. Never evaluate a kernel at age 0. At t = 900 s, the edge at 900 s is excluded by both routes.
- Prints go through I78's rank-0 `_p`. Use `-s`, and print ≤ 7 digits.
- **The leak guard:** a hook refusal on a printed θ line moves the log unchanged to `logs/leak-held/`. Never `--no-verify`.

**Scope.** Infinite gel, the zero-radius line and the two point tips, 1 W per rod, four placements on the axis. **Not claimed:**
- the box walls or the rod's own conduction (I75 / I81);
- a ball-averaged read;
- the capsule's caps as a source;
- any F2182 figure or published value.

The rungs are placements, not a convergence ladder. No rate is read. The twin's point source sits 1 mm from rung 1, a point's 1/R there, so it is printed beside the line, never compared.

**Negative result.**
- (Lim), (Bq) or (Dut) red at any rung: a derivation or quadrature error. Stop, BLOCKED, band unmoved, every value printed.
- (Rep) red: I78's module no longer reproduces its record. BLOCKED.
- A control red contradicts the steady arithmetic above. Report the measurement and stop.
- No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 18, a 🟡 record line: the lumped probe-placement ladder, gated at all four rungs. §10 I80 is DONE. I76's lumped column is gated at the clinical placement, and I81 has a closed form at each FEM rung.

**Depends on:** nothing (items 874, 786, 787 and 811 are on `main`).

### §9 item 895 — 2026-11-08 (§10 I81, chain I — daily, 2026-11-08, provisional, rule (6); IMP-8 step 19, the bracket at the clinical probe placement; chain step (solve), heavy, real, -n 8)

*Written by the 2026-11-08 18:00 daily review under §9 rule (6), toward I76's stated readout (`PROJECT_PLAN.md:5068–5071`). Bands, controls and anchors are imported, and none is set. Family I81, attempt 1 of 4. **Executor:** the implementer. Depends on nothing; item 894 is printed beside it if landed, never imported.*

**Why.** I76 will set entry 63's clinical probe beside I75's bracket arms. The two probes differ:
- the clinical probe is 1 mm beyond a capsule apex, at z = 17 mm (item 893's (P));
- I75's is 1 mm beyond a flat-ended source, at z = 16 mm (`test_imp8_rod_probe_f2182_fem.py:74`).

On I78's line source that 1 mm alone moves θ(900 s) by 29 % (`20261108T133157Z_IMP-8-item874-W1.log:54`). This step reads I75's three arms **at the clinical placement** and on a short ladder past it, on I75's fixture and solve unchanged. I76's bracket rows can then be read at the same distance from the source as the clinical row. If W3 moves the offset, the ladder still has a rung for it.

**Code.** New `tests/validation/test_imp8_rod_probe_ladder_fem.py`, real build, degree 2. No `src/` edit and no landed module edited. (I75's fixture keeps no field, so its module cannot be extended by import.)
- **Imports, unchanged** (read every signature first):
  - I75 (`test_imp8_rod_probe_f2182_fem`): `ARM_ENV`, `ARMS`, `PROBE_PTS`, `SNAP_STEPS`, `DUTY`, `GEL_NEW`, `TI_THETA`, `_arm`, `_body_cells`, `_sixth`, `_superpose`. **Not** its `probe_run` fixture: importing it would make a second solve available to collection;
  - everything I75 imports for its fixture body (`:39–70`), from the same modules, including the `box` fixture from I55 (`# noqa: F401`, as I75 does).
- **The knob is I75's** `FEM_EM_IMP8_I75_ARM` (`ti` / `gel` / `hole`, junk refused), read through `_arm()`. Put `FEM_EM_IMP8_I55_PHASE=transient` on the command.
- **The fixture `ladder_run(box)`** is I75's `probe_run` body (`:133–208`) with the solve call byte-identical (`:148–151`, `snapshot_steps=SNAP_STEPS`) and three reads, **one batch per call** (the batch rule):
  1. `_theta_pts(res.temperature, READ_PTS)`: I73's batch, unchanged;
  2. `_theta_pts(·, PROBE_PTS)` on the final field and each snapshot: I75's batch, unchanged, so the reproduction below compares call for call;
  3. **new:** `_theta_pts(·, LADDER_PTS)` on the final field and each snapshot, with `LADDER_PTS` = ((0, 17), (0, −17), (0, 18), (0, −18), (0, 20), (0, −20)) mm, i.e. L/2 + 2, 3 and 5 mm. Assert those z values from `L` and `MM` in a no-mesh `test_k0_ladder_points`.

**Gates (asserted, bands imported).**
- **Every arm, I75's:**
  - (E) ≤ `E_BAND`; (S0) ≤ `B_BAND`;
  - cells = `W0C_CELLS` (485 291), body cells = `BODY_CELLS` (6 495), no non-finite coefficient dofs;
  - snapshots at 225 / 450 / 675 s;
  - θ_probe(900 s) > θ_probe(675 s) at both ±16 mm probes.
- **(R75) every arm, the anchor.** I75's printed probe reads are reproduced within one unit of the sixth significant digit (`_sixth`), at +16 and −16 mm, for θ(900 s) and for the D = 0.5 superposed θ(900 s):
  - ti `1.064117e+01` / `1.049960e+01`, D `2.109553e+00` / `2.099602e+00` (`20261108T194323Z_IMP-8-item881-W1.log:120–121`);
  - gel `9.220827e+00` / `9.102961e+00`, `1.979610e+00` / `1.970586e+00` (`20261108T195853Z_IMP-8-item881-W2b.log:120–121`);
  - hole `8.815610e+00` / `8.750725e+00`, `1.895292e+00` / `1.888582e+00` (`20261108T200636Z_IMP-8-item881-W3.log:120–121`).

  I75 reproduced I73's θ across sessions at `-n 8` to ≤ 2.619e-06 against the 1e-04 unit (`…item881-W1.log:113–118`).
- **Arm `gel`:** I73's six θ (`{**I55_THETA, **GEL_NEW}`) and I55's (I) ≤ `INTERIOR_BAND` at `GATE_PTS`.
- **Arm `ti`:** `TI_THETA`.
- **Arm `hole`:** |θ(2, 0)/`I55_THETA`[(2, 0)] − 1| > `CONTROL_FLOOR`. These three are I75's arm gates, unchanged.

**Control (I78's, on the FEM; asserted at factor 1 on `INTERIOR_BAND`).** At 900 s under the continuous source, at each sign of z and on every arm:
- θ(16) > θ(17) > θ(18) > θ(20) strictly;
- 1 − θ(17)/θ(16) > `INTERIOR_BAND`.

The ceiling is 1 (θ > 0). The infinite-medium line predicts 0.2941583 for the first drop (`…item874-W1.log:54`). The FEM's source has a 1 mm radius and a body, which is predicted to flatten that, as a sign only. A ladder that does not fall would mean the read cannot see its own placement.

**Printed (🟡 record, no band, ≤ 7 digits, K/W, 1 W in the rod).** Per arm, at each ladder point:
- θ(900 s); θ at 225 / 450 / 675 s; θ(900)/θ(675);
- the D = 0.5 superposed θ(900 s), by I75's `_superpose`, a superposition and not a switched solve (I75's disclosure 2);
- θ(z)/θ(16) at the same sign; the ± mean and |difference| (mesh asymmetry).

Beside each, labelled "infinite medium, line source, no rod": `QP·_phi(0, z·MM, T_PRINT)` and the line's `_duhamel(0, z·MM, T_PRINT, DUTY)` computed inline, as I75 does at `:191–194`, with FEM / closed form. The slot report gives, per arm, the clinical-placement (17 mm) mean.

**Windows.** Before each, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Every window is a bare foreground harness call, Bash timeout 660000 ms. The heavy windows use the durable capture of `20261108T194323Z_IMP-8-item881-W1.log:12` verbatim, including the trailing `; exit $rc`. Its logged Command line must read as that one does, with this module and the label changed.
- **W0** `IMP-8-item895-W0`, `-n 2`, no mesh:
  `scripts/testing/run_and_log.sh IMP-8-item895-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/validation/test_imp8_rod_probe_ladder_fem.py -k k0 -v -s --tb=short'"`
  I75's W0 took 4 s (`20261108T194310Z_IMP-8-item881-W0.log:60`).
- **W1** `IMP-8-item895-W1`, arm `ti`, `-n 8`:
  `scripts/testing/run_and_log.sh IMP-8-item895-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src FEM_EM_IMP8_I55_PHASE=transient FEM_EM_IMP8_I75_ARM=ti timeout -k 30 600 mpiexec -n 8 python3 -m pytest tests/validation/test_imp8_rod_probe_ladder_fem.py -v -s --tb=short > /workspace/logs/imp8-item895-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp8-item895-w1-raw.log; cat /workspace/logs/imp8-item895-w1-raw.log; exit \$rc'"`
  I75's W1 took 470 s (t_call 449.6 s; `…item881-W1.log:112`, `:210`).
- **W2** `IMP-8-item895-W2`, arm `gel`: W1's command with `gel` and `w2`. I75's W2b took 441 s (`…item881-W2b.log:214`).
- **W3** `IMP-8-item895-W3`, arm `hole`: with `hole` and `w3`. I75's W3 took 459 s (`…item881-W3.log:210`).

The ladder batch adds four `_theta_pts` calls on a solved field, so predict ≈ 445–485 s per arm window against the 600 s ceiling. **Why `-n 8` and not `-n 2`:** one CN run on this 485 291-cell box is ≈ 420–450 s of `t_call` at `-n 8`, and at `-n 2` it would not return inside 600 s. A kill parks the item. Never raise the ceiling, and never run two heavy windows at once.

**Tier.** Heavy, real, `-n 8`. **Class:** chain step (solve). **Cost:** ≈ 8.2 measured / 35 predicted min (W0 plus three ≈ 7.6-minute windows).

**Traps already paid for.**
- **I75's W2 crash:** `dict(I55_THETA, **GEL_NEW)` raises `TypeError` on tuple keys (`20261108T195128Z_IMP-8-item881-W2.log`). Use `{**I55_THETA, **GEL_NEW}`.
- **I75's rule (i):** all three arm windows run on one committed module. If a fix lands after W1, re-run W1.
- **The batch rule:** the PROBE batch is I75's, unchanged and separate from the ladder batch. (R75) compares call for call; never fold the ladder into the probe batch.
- **The hole arm** (item 881's trap): the rod core reaches ≈ 4e6 K/W at ε = `EPS_HOLE`. Print no rod-interior or max θ. Every ladder point is ≥ 2 mm beyond the source cells' |z| ≤ 15 mm, and I75's 16 mm hole read was 8.8 K/W.
- **Do not call I73's `_coefficients` with `hole`:** it overwrites only on `ti` (`test_imp8_rod_self_heated_fem.py:109–111`). Use I66's, as I75 does.
- Real build only. Reduce every rank-local count before asserting, and count verdicts per rank by token. Use `-s`.
- **The leak guard:** a hook refusal moves the log unchanged to `logs/leak-held/`, and the slot reports it. Never `--no-verify`.

**Scope.** I75's mesh, Δt, nominal Ti, perfect contact and uniform rod source, read at three more placements on the axis. **Not claimed:**
- the clinical geometry: a capsule in the ASTM box, heated by gel SAR with the rod unheated (entry 63);
- a like-for-like clinical comparison (I76's to label);
- a switched-source duty solve;
- a gate on the rod's effect at the probe (no closed form has the finite conducting rod);
- any F2182 value (W3).

**Negative result.**
- (E), (S0), a precondition, or an arm gate red: stop, BLOCKED (I75's dispositions).
- **(R75) red:** this route is not I75's. Print every θ, stop, BLOCKED; bank no ladder read without it.
- **The control red** (θ not falling outward along the axis beyond a heated rod end in a gel box, from rest under a continuous source) contradicts the physics. Print every ladder read, stop, BLOCKED, and loosen nothing.
- No `PREDICTED RED:` line.

**Status it can move.** `IMP-8` gains step 19, a 🟡 record line: the gel / Ti / hole bracket at the clinical probe placement and on a 1 / 2 / 3 / 5 mm ladder, continuous and D = 0.5, per watt. §10 I81 is DONE. I76's bracket rows can be read at the clinical row's distance from the source.

**Depends on:** nothing (I75's module and logs, I73's, I66's and I55's modules are on `main`).

**Outcome (00:00, 2026-11-09 slot, implementer, delegated, foreground): DONE.** Attempt 1 of 4. `IMP-8` step 19 is a 🟡 record; §10 I81 is DONE. New `tests/validation/test_imp8_rod_probe_ladder_fem.py` (imports as the brief, no `src/` or landed-module edit; not edited after W0, so all three arm windows ran on one module text). `pgrep -c python3` read 0 before each heavy window. W0 `20261109T060323Z_IMP-8-item895-W0.log` (`-n 2`, 1 passed per rank, 5 s). W1 ti `20261109T060343Z_IMP-8-item895-W1.log` (8 passed / 3 skipped × 8 ranks, 424 s, t_call 400.7 s); W2 gel `20261109T061107Z_IMP-8-item895-W2.log` (9 / 2 × 8, 417 s, t_call 395.4 s); W3 hole `20261109T061832Z_IMP-8-item895-W3.log` (8 / 3 × 8, 448 s, t_call 425.5 s); every `[capture] rc=0`, verdicts counted per rank by token. (E) 2.776e-14 / 1.954e-14 / 2.975e-14 ≤ `E_BAND`; (S0) 4.441e-16; cells 485 291, body cells 6 495, non-finite 0; arm gates green (gel I73 worst |Δ| 4.345e-06, (I) ≤ 3.216e-03; hole θ(2, 0) off I55 by 4.3944e-03). (R75) every arm inside one sixth-digit unit (worst ti +16 mm θ(900 s) 1.993e-06 against 1e-04). Control green on every arm and sign: θ(16) > θ(17) > θ(18) > θ(20), first drop ti 0.295119 / 0.294595, gel 0.273386 / 0.272288, hole 0.262420 / 0.263417 (+z / −z) > `INTERIOR_BAND`; the ti drop is steeper than the infinite line's 0.2941583 by ≈ 0.001, gel and hole flatter — printed, not asserted. Record at the clinical placement (|z| = 17 mm, ± mean, K per watt, 900 s): ti 7.453618 / gel 6.662160 / hole 6.473925 (± difference 9.429e-02 / 7.565e-02 / 5.659e-02); D = 0.5 by superposition 1.903297 / 1.818740 / 1.765061; FEM / line QP·Φ 1.094442 / 0.978229 / 0.950590, FEM / line `_duhamel` 1.044484 / 0.998081 / 0.968624 (`W1:127`, `W2:127`, `W3:127`). Disclosed: item 894's line record is printed beside from a constant citing `20261109T043829Z_IMP-8-item894-W1.log:98–107`, not imported; `ARM_ENV` is imported per the brief's list but unused.

### §9 item 911 — 2026-11-09 (§10 I83, chain I — daily, 2026-11-09, provisional, rule (6); IMP-8 step 20, the rod's finite radius at the probe, closed form; lane-3 closed form → chain step, smoke, numpy)

*Written by the 2026-11-09 03:00 daily review under §9 rule (6), toward I76's stated readout (entry 63's probe "beside I75's bracket arms and I78's closed form", `PROJECT_PLAN.md:5077–5080`). Bands are imported and none is set or moved. **One new assertion is flagged for the weekly to ratify or revert:** (G) below carries I55's (I) comparison (FEM against closed form at `INTERIOR_BAND`) to the probe points, with the comparand's radius set to the source's. Family I83, attempt 1 of 4. **Executor:** the implementer. Depends on nothing.*

**Why.** Item 895's Outcome flags a sign: ti's first FEM drop (16 → 17 mm) is steeper than the infinite line's 0.2941583, where the brief predicted flatter. The logged first drops are ti 0.295119 / 0.294595 (`20261109T060343Z_IMP-8-item895-W1.log:149–150`), gel 0.273386 / 0.272288 (`20261109T061107Z_IMP-8-item895-W2.log:134–135`) and hole 0.262420 / 0.263417 (`20261109T061832Z_IMP-8-item895-W3.log:134–135`), all +z / −z.

Two readings are open:
- **(a) discretisation.** The ti-vs-line excess is +0.000961 / +0.000437. That is the same size as each arm's +z/−z spread of the first drop (gel 0.001098, ti 0.000524, hole 0.000997).
- **(b) physics.** The brief's "flatter" was a prediction about the source's 1 mm radius, and the conducting rod is a second effect that I81 could not separate from it.

I75 / I81 printed every arm beside the **zero-radius** line. I55's source, however, is the DG0 cylinder r ≤ `ROD_RADIUS`, |z| ≤ L/2 (`test_imp8_rod_line_source_fem.py:152–177`; effective radius 1.0022 mm, `20261103T182153Z_IMP-8-I55-W1s.log:56`), and 1 mm from its end face the line is the wrong comparand.

The uniform finite cylinder has a closed form on its axis. The disc integral of the point kernel is exact:

2π∫₀ᵃ erfc(√(r′² + d²)/s)/√(r′² + d²) r′ dr′ = 2πs[ierfc(d/s) − ierfc(√(a² + d²)/s)], with s = 2√(κt).

Reading the gel arm against that form separates the two readings:
- if the gel arm sits on it, the gel arm's flattening is the radius, not the mesh;
- ti − gel is then the rod's effect, on a gated baseline, on the same mesh and source cells.

**Unlogged drafter arithmetic (numpy Gauss–Legendre, host; printed beside by the module, never asserted as values).** SCALE·Φ_a(900 s) at 16 / 17 / 18 / 20 mm is 9.197899 / 6.682151 / 5.184508 / 3.407212 K/W. The same code reproduces I80's line 9.648662 / 6.810428 / 5.242317 / 3.427288 (`20261109T043829Z_IMP-8-item894-W1.log:98`, `:101`). Steady values are 14.72925 / 12.12776 / 10.54176 / 8.581375. The cylinder's first drop is 0.273513, against the line's 0.2941583.

Against I81's ± means, (G) is predicted at:
- continuous: −3.914e-03 / −2.992e-03 / −2.420e-03 / −1.921e-03;
- D = 0.5: +6.4e-04 / +3.6e-04 / +1.9e-04 / +3.7e-05.

The decomposition the record will carry, on the ± means, is:

| Term | First-drop change |
|---|---|
| radius (cylinder − line) | −0.020645 |
| mesh (gel − cylinder) | −0.000673 |
| rod (ti − gel) | +0.022019 |

So reading (b) is predicted: the radius flattens, the rod steepens by slightly more, and ti's net +0.0007 over the line is the near-cancellation of two effects ≈ 30× larger.

**Code.** New `tests/validation/test_imp8_probe_cylinder_source_closed_form.py`, numpy / scipy, the image's default (real) build. No `src/` edit and no landed module edited.
- **Imports, unchanged** (read every signature first):
  - I80 (`test_imp8_probe_ladder_closed_form`): `RUNGS`, `LADDER_MM`, `_lab`;
  - I78 (`test_imp8_probe_f2182_closed_form`): `DUTY`, `N_GL`, `_axis_steady`, `_p`;
  - I58 (`test_imp8_entry48_lumped`): `_gl01`, `BQ_BAND`;
  - step 0 (`test_imp8_rod_line_source_closed_form`): `L`, `MM`, `KAPPA`, `QUAD`, `T_PRINT`, `_phi`;
  - I59 (`test_imp8_rod_line_source_duty_closed_form`): `SCALE`, `T_STAR`, `PERIOD`;
  - I62 (`test_imp8_rod_tips_duty_closed_form`): `_duhamel`, `_convolution`, the kernel-generic forms;
  - `ROD_RADIUS` (`tests.mesh.test_birdcage_implant_rod`);
  - `B1_BAND` (`test_imp4_pec_cylinder_closed_form`) and `INTERIOR_BAND` (`test_mat13_transient_perfusion_sphere`), by name.
- **The kernel (per unit SCALE, so θ = SCALE·Φ_a; Φ_a → `_phi` as a → 0):**
  - `_ierfc(u)` = exp(−u²)/√π − u·erfc(u).
  - `_phi_cyl(r, z, t, a=ROD_RADIUS, length=L)`. Assert r == 0.0 and |z| > L/2: the axis beyond the end only, one piece. With d = |z − z′| and D = hypot(a, d), it is (2s/a²)·∫ [`_ierfc`(d/s) − `_ierfc`(D/s)] dz′ over [−L/2, L/2] by `quad` with `QUAD` unchanged. At t = `np.inf` the kernel is the steady 2/(D + d).
  - `_phi_cyl_dot(r, z, u)` is closed, with no quadrature: (s²/(2a²u))·(−expm1(−a²/s²))·[erf((|z| + L/2)/s) − erf((|z| − L/2)/s)], s² = 4κu, and 0 for u ≤ 0. Derivation, for the module docstring: d/dt[s·ierfc(x/s)] = (s/2t)·e^{−x²/s²}/√π, and e^{−d²/s²} − e^{−D²/s²} = e^{−d²/s²}(1 − e^{−a²/s²}). The a → 0 limit is the line's (1/2t)[erf − erf]. The drafter checked ∫₀⁹⁰⁰ Φ̇_a du = Φ_a(900) to 1.0e-13 at 16 / 17 / −16 mm (unlogged).
  - `_cyl_steady(z)` = g(|z| + L/2) − g(|z| − L/2), with g(u) = u/(√(a² + u²) + u) + asinh(u/a). This is the cancellation-free antiderivative of 2∫(√(a² + u²) − u) du / a².
  - `_gl_cyl(z, t, n)` is I78's `_gl_line` pattern: `_gl01(n)` mapped onto [−L/2, L/2], on the same integrand.
- **The FEM reads, custody by string (I79's (K-b) pattern).** Module literals for I81's reads, each with its log:line:
  - the ± mean θ(900 s) and D = 0.5 mean per arm and rung, from the `ladder |z|=NN mm … arm X: +/- mean theta(900 s) V K/W … D=0.5 mean V` lines: gel W2 `:124 / :127 / :130 / :133`, ti W1 `:124 / :127 / :133 / :145`, hole W3 `:124 / :127 / :130 / :133`;
  - the per-sign θ at 16 / 17 / 18 / 20 mm, from the `control sign ±1 arm X: theta 16/17/18/20 …` lines: W1 `:149–150`, W2 / W3 `:134–135`.

  `test_k_reads_from_logs` opens the three tracked logs under `/workspace/docs/testing/logs/`, finds each line by `re.search` on its content (not `startswith`; see Traps), and asserts every literal appears there as the same `.6e` / `.7g` string. It also asserts each arm's log has exactly four mean lines and two control lines (the drafter counted 4 / 2 in each). The gel literals, for the brief's own check: 9.161894e+00 / 6.662160e+00 / 5.171959e+00 / 3.400665e+00, and D means 1.975098e+00 / 1.818740e+00 / 1.665104e+00 / 1.371623e+00. The ti means are 1.057039e+01 / 7.453618e+00 / 5.672034e+00 / 3.653661e+00, and the hole means 8.783168e+00 / 6.473925e+00 / 5.068436e+00 / 3.361168e+00.

**Gates (asserted, bands imported).** At every rung of `RUNGS`, both signs where the form is evaluated:
- **(Lim)** |`_phi_cyl`(0, z, ∞) / `_cyl_steady`(z) − 1| ≤ `B1_BAND`. This is the z′ quadrature against the analytic antiderivative. Predicted round-off: the drafter's 400-point GL read 7.6e-15 … 9.2e-14.
- **(Bq)** |`_gl_cyl`(z, t, `N_GL`) / `_gl_cyl`(z, t, 2·`N_GL`) − 1| ≤ `BQ_BAND` at t = `T_PRINT` and ∞. The integrand is entire in d except the branch points of √(a² + d²) at d = ±ia, farther from the endpoint than the line's real pole 1 mm away. Predicted at or below I80's 1.499e-14 … 4.030e-14 (`20261109T043829Z_IMP-8-item894-W1.log`, I80's Outcome).
- **(Dut)** |`_duhamel`(`_phi_cyl`, 0, z, `T_STAR`, `DUTY`) − `_convolution`(`_phi_cyl_dot`, 0, z, `T_STAR`, `DUTY`)| / (`SCALE`·`_cyl_steady`(z)) ≤ `B1_BAND`. The step sum runs on the quadrature route and the convolution on the closed-form Φ̇_a, so the two routes are independent. I80's line precedents are ≤ 2.486e-16.
- **(G), the anchor; flagged for ratification.** I81's gel arm, z-even part: |θ̄_gel(z) / (`SCALE`·`_phi_cyl`(0, z, `T_PRINT`)) − 1| ≤ `INTERIOR_BAND`, and |θ̄ᴰ_gel(z) / `_duhamel`(`_phi_cyl`, 0, z, `T_STAR`, `DUTY`) − 1| ≤ `INTERIOR_BAND`, at all four rungs. θ̄ is the ± mean.
  - **Why the mean and not each sign:** the source and box are even in z, so the odd part of the FEM read is mesh asymmetry. I75 / I81 record it as such: |+ − −| is 1.179e-01 on the 16 mm mean of 9.161894, i.e. 1.29 % (W2 `:124`), larger than the band by itself.
  - **Per sign is printed only.** It is predicted at +2.493e-03 / −1.032e-02 at ±16 mm, so a per-sign gate would be red at −16 mm by I81's own recorded asymmetry. This step does not pose it, and says so on the print line.

**Controls (asserted at factor 1 on `INTERIOR_BAND`; ceilings computed).**
- **(C1) the blind comparand.** The zero-radius line, which I75 / I81 printed beside, misses the gel mean by more than the band at 16 and 17 mm: 1 − θ̄_gel/(`SCALE`·`_phi`) = 0.050449 / 0.021771 (from `mean FEM/QP*Phi 0.949551 / 0.978229`, W2 `:124`, `:127`).
  - **Ceiling:** on an exact FEM the line's miss is 1 − Φ_a/Φ = 0.04672 / 0.01883 / 0.01103 / 0.00586 (drafter arithmetic). It clears the band only at 16 and 17 mm, so (C1) is asserted there and printed at 18 / 20 mm. The separation from (G) is ≈ 12.9× at 16 mm and ≈ 7.3× at 17 mm.
- **(C2) the gate sees the arm.** The ti and hole means miss `SCALE`·Φ_a by more than the band at 16 and 17 mm. Predicted: ti 0.14922 / 0.11545, hole 0.04509 / 0.03116; printed at 18 / 20 mm (hole 0.0135 at 20 mm is near the band). A reader that took the wrong arm, or a gate that any arm passes, fails here.

**Printed (🟡 record, no band, ≤ 7 digits, K per watt, 1 W in the rod, infinite gel).** Head the block "placement and time quoted from memory; W3 confirms or corrects".
- **Per rung:** `SCALE`·Φ_a at 900 s continuous, at D = 0.5, and steady; Φ_a/Φ (the radius factor); the same at a = 1.0022 mm (I55's effective radius, `20261103T182153Z_IMP-8-I55-W1s.log:56`; the drafter's 16 mm value is 9.196238, a 1.8e-04 shift).
- **Per arm and rung:** FEM ± mean / Φ_a and per sign / Φ_a; ti/gel − 1 and hole/gel − 1, the rod's effect on the gated baseline (predicted ti +0.1537 / +0.1188 / +0.0967 / +0.0744).
- **The first-drop decomposition, on the means and per sign:** line, cylinder, gel, ti and hole drops; radius = cylinder − line; mesh = gel − cylinder; rod = ti − gel; and each arm's +z/−z spread of the drop.
- **The a → 0 check, printed:** Φ_a/Φ − 1 at a = `ROD_RADIUS`/10, /20 and /40, with the observed order log₂ of successive deficits (expected ≈ 2). It is not asserted, since no rate band is imported for it.

**Window.** Before it, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. **W1** `IMP-8-item911-W1`, smoke, real, `-n 2`. Make one bare foreground harness call, Bash timeout 660000 ms; I80's module rides along as its own regression:

`scripts/testing/run_and_log.sh IMP-8-item911-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/validation/test_imp8_probe_cylinder_source_closed_form.py tests/validation/test_imp8_probe_ladder_closed_form.py -v -s --tb=short'"`

**Price.** I80's window ran 12 tests in 2.46 s (`20261109T043829Z_IMP-8-item894-W1.log:142`), 4 s harness (`:149`). The new module adds a `quad` per Φ_a call, four Φ_a calls per (Dut) step sum, a 1-D `quad` of a closed form per (Dut) convolution, and three log reads. Predict ≤ 30 s against the 120 s ceiling.

**Tier.** Smoke, numpy, `-n 2`. **Class:** lane-3 closed form → chain step. **Cost:** ≈ 4.1 measured / 15 predicted min.

**Traps already paid for** (items 874 / 894's list, plus this module's own).
- `_phi_steady` and I58's `_line(·, ·, ∞)` divide by r. Never call them on the axis; use `_axis_steady` for the line.
- **Cancellation.** Write the steady kernel as 2/(D + d), never (D − d)·2/a², and the antiderivative as g(u) above, never u√(a² + u²) − u². The naive forms lose ≈ 3 digits at u = 31 mm, which is at `B1_BAND`. Write Φ̇_a with `expm1`.
- `_step` (I62) guards t ≤ 0, and `_phi_cyl_dot` returns 0 for u ≤ 0. Never evaluate a kernel at age 0. At t = 900 s the edge at 900 s is excluded by both routes.
- **The log lines are interleaved with pytest progress.** The ti `control sign +1` line begins `tests/validation/…::test_w_probe_rises_675_to_900 [imp8-i81] …` (W1 `:149`). Match with `re.search`, not `startswith`, and never re-sample a FEM value.
- **If a cited log has been gzipped by rotation** before the slot, read `<name>.log.gz` with `gzip.open`; the line content is unchanged. If both are absent, (K) is red: BLOCKED, and never retype a value from this brief.
- Prints go through I78's rank-0 `_p`. Use `-s`, and print ≤ 7 digits.
- **The leak guard:** a hook refusal on a printed θ line moves the log unchanged to `logs/leak-held/`, and the slot reports it. Never `--no-verify`.

**Scope.** An infinite gel medium, the uniform finite cylinder (nominal `ROD_RADIUS`, length L) on its axis beyond the end, 1 W, four placements. I81's gel / ti / hole FEM reads are taken by custody from the logs, with no solve. **Not claimed:**
- a per-sign gate (the mesh asymmetry is recorded, not gated);
- an h-statement: the mesh is not varied, and "mesh = gel − cylinder" is a residual on one mesh, not a convergence;
- the box walls: I55's (I) holds at 10–20 mm from the axis at 900 s, and nothing here is read nearer the walls;
- any closed form for the conducting (ti) or insulating (hole) rod: ti/gel and hole/gel are records against the gated gel baseline, not gates;
- the clinical capsule, whose source is gel SAR with the rod unheated (entry 63);
- any F2182 figure or published value (W3).

The rungs are placements, not a convergence ladder. **I76's comparand set is not changed:** whether I76 prints Φ_a beside I78's line is the weekly's to decide.

**Negative result.**
- (Lim), (Bq) or (Dut) red at any rung: a derivation or quadrature error. Stop, BLOCKED, band unmoved, every value printed.
- (K) red: the logs do not carry the literals. BLOCKED; never retype.
- **(G) red at any rung:** the gel arm at the probe is not the finite-radius source within the band. Report every rel, both signs, and stop. I81's record stays as it stands, and reading (b)'s radius half is refuted. Never switch the comparand to a = 1.0022 mm or to a single sign to rescue it, and never loosen.
- (C1) red contradicts I81's logged 0.949551 / 0.978229: a reader error, BLOCKED.
- (C2) red: the reader took the wrong arm, or the arms are indistinguishable at the probe. Report and stop.
- No `PREDICTED RED:` line.

**Status it can move.**
- `IMP-8` gains step 20, a 🟡 record line: the finite-radius rod source in closed form on the axis, gated by (Lim) / (Bq) / (Dut), with I81's gel-arm probe reads gated against it on their z-even part at 16 / 17 / 18 / 20 mm, continuous and D = 0.5. §10 I83 is DONE.
- Item 895's disclosed sign flag is resolved into the decomposition: radius and rod, each ≈ 0.02 on the first drop, with ti's +0.0007 over the line as their residue.
- I76's bracket rows can quote the rod's effect at the clinical placement as ti/gel − 1 on a gated baseline instead of against the zero-radius line.
- A red (G) is equally a status: the gel row at the probe stays unanchored, with the measured miss recorded on step 20.

**Depends on:** nothing. Items 874, 894, 895, 786, 787 and 811 and their logs are on `main`.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 911 — `IMP-8` step 20 (§10 I83) — per `docs/planning/chunks/IMP-8.md` "§9 item 911". Read the brief and any review ruling beneath it first; the rulings hold. Create only `tests/validation/test_imp8_probe_cylinder_source_closed_form.py`; edit no landed module and no `src/` file; import module objects only, never `test_*`. No gmsh, no solve: the FEM reads come from the three tracked I81 logs by string, never retyped and never re-sampled. Before W1, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; run pytest with -s; container `timeout -k 30`. Never move a band, and never change (G)'s comparand (nominal `ROD_RADIUS`, ± mean) to rescue a red. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-8-step20-*`, mark item 911 BLOCKED in the same commit, and report."

**Outcome (2026-11-09, 04:30 slot, fifth item): DONE.** New `tests/validation/test_imp8_probe_cylinder_source_closed_form.py`; no `src/` edit, no landed module edited. W1 `20261109T110006Z_IMP-8-item911-W1.log` — smoke, real, `-n 2`, container `timeout -k 30 120`: 14 passed on each rank (28 `PASSED` tokens, 0 `FAILED`), 2.15 / 2.21 s pytest, 4 s harness; I80's module rode along green. `pgrep -c python3` read 0 before the window.
- **(Lim)** worst 2.220e-16 (`:59`), both signs at every rung, against `B1_BAND`. **(Bq)** worst 3.542e-14 (`:80`) against `BQ_BAND`. **(Dut)** worst 8.971e-16 (`:90`) against `B1_BAND`. **(K)** every literal found as the same string; 4 mean / 2 control lines per arm (`:92–94`).
- **(G), flagged for the weekly to ratify or revert (a new read at the imported `INTERIOR_BAND`).** Gel ± mean / SCALE·Φ_a − 1 = −3.914e-03 / −2.992e-03 / −2.420e-03 / −1.921e-03 at 16 / 17 / 18 / 20 mm; D = 0.5 +6.380e-04 / +3.610e-04 / +1.922e-04 / +3.656e-05 (`:97–110`); worst 3.914e-03 (`:111`), all as predicted. Per sign, printed only: +2.493e-03 / −1.032e-02 at 16 mm (`:99`) — the −16 mm read alone is outside the band, as the brief predicted; not gated.
- **Controls (asserted at factor 1 on `INTERIOR_BAND` at 16 / 17 mm).** (C1) 1 − gel/(SCALE·Φ_line) 0.050449 / 0.021771, against exact-FEM ceilings 0.04672 / 0.01884; (C2) ti 0.14922 / 0.11545, hole 0.04509 / 0.03116 (`:113–114`). Printed at 18 / 20 mm: (C1) 0.013421 / 0.007768, hole 0.02239 / 0.01351 (`:115–116`).
- **Record (🟡).** SCALE·Φ_a(900 s) 9.197899 / 6.682151 / 5.184508 / 3.407212 K/W, steady 14.72925 / 12.12776 / 10.54176 / 8.581375, Φ_a/Φ 0.9532823 / 0.9811647 / 0.9889726 / 0.9941423; at a = 1.0022 mm, 9.196238 at 16 mm (`:119–139`). ti/gel − 1 +0.1537 / +0.1188 / +0.0967 / +0.0744; hole/gel − 1 −0.0413 / −0.0283 / −0.0200 / −0.0116. First drop on the means: line 0.294158, cylinder 0.273513, gel 0.272840, ti 0.294859, hole 0.262917; **radius −0.020645, mesh −0.000673, rod +0.022018, ti − line +0.000700** (`:141–142`) — reading (b): item 895's sign flag is the near-cancellation of the radius and the rod. a → 0: deficits −5.6912e-04 / −1.4255e-04 / −3.5653e-05, observed order 1.9973 / 1.9993 (`:146`), printed only.
- **Disclosure.** The bracket ierfc(d/s) − ierfc(D/s) is evaluated as the identical ∫ erfc(v) dv over [d/s, D/s] by 16-point Gauss–Legendre, to avoid the two-ierfc cancellation; `_ierfc` is kept, and the two forms are printed side by side and agree to 3e-17 absolute at one test point (`:147`). The steady kernel is 2/(D + d) and Φ̇_a uses `expm1`, as briefed. (G) and the record were evaluated at both ±z (the kernel is even); per-sign FEM reads are printed against +z's form.

### §9 item 934 — 2026-11-10 (§10 I86, chain I — rule (4), the 11-09 weekly's step; IMP-8 step 21, the finite-radius two-tip twin at the probe, closed form; lane-3 closed form → chain step, smoke, numpy)

*Enumerated by the 11-09 weekly (`PROJECT_PLAN.md:4964–4974`). Written by the 2026-11-10 03:00 daily review as execution detail. The bands are imported, and none is set or moved. **One execution detail is flagged for the weekly:** the step's "short cylinders" get a length here (below). Family I86, attempt 1 of 4. **Executor:** the implementer. Depends on nothing; no queued window's module is touched (numpy only, new module).*

**Why.** I83 measured how much the uniform source's radius moves the 16 mm read: 1 − Φ_a/Φ = 0.0467 (`20261109T110006Z_IMP-8-item911-W1.log:119`). The point twin reads 6.543912 × the line there (`20261109T043829Z_IMP-8-item894-W1.log:100`). The twin's point sources sit 1 mm from the probe, where a 1/R kernel is most sensitive to the source's size. So it is open whether the radius or the tip weighting dominates the twin's probe read.

**The fixture (execution detail, flagged).** Each tip is a uniform cylinder of radius a = `ROD_RADIUS` and length h = 2a, centred on z = ±L/2 and carrying P/2 (I58's split, `test_imp8_rod_line_source_closed_form.py:89–91`).
- It circumscribes I56's tip ball, radius `ROD_RADIUS` about ±L/2 (this file, item 846's `_rod_cells`).
- As a → 0 it collapses to I58's points, which is what (Lim-a) needs.
- At rung 1 (16 mm) the probe sits on the near cylinder's end face, as I56's balls put their surface there (`test_imp8_probe_f2182_closed_form.py:31–32`). On the axis d ≥ 0 and 2/(D + d) ≤ 2/a, so the kernel is finite.
- A second length is printed beside as a sensitivity: the cylinder [L/2 − a, L/2] (h = a, inside the rod's end).

**The kernel**, per unit `SCALE` (θ = `SCALE`·Φ₂,a):

Φ₂,a(0, z, t) = (L/2)·[M(|z − L/2|, t) + M(|z + L/2|, t)], with:
- M(ζ, t) = ½∫₋₁¹ B(|ζ − a u|, t) du, the segment mean in the unit variable u. Assert ζ ≥ a: the probe is on or beyond the face.
- B(d, t) = (2s/a²)[ierfc(d/s) − ierfc(D/s)], with D = hypot(a, d) and s = 2√(κt). At t = ∞, B = 2/(D + d). This is I83's disc integral (`test_imp8_probe_cylinder_source_closed_form.py:5–13`).
- Steady: M(ζ, ∞) = [g(ζ + a) − g(ζ − a)]/(2a), with I83's `_g` (`:122–123`).
- Rate: per tip, (L/(2h))·`_phi_cyl_dot`(0, ζ, u, a=a, length=h), I83's closed form (`:112–119`). Use it at a = `ROD_RADIUS` only.

**Code.** New `tests/validation/test_imp8_probe_twin_cylinder_closed_form.py`, numpy / scipy, the image's default (real) build. No `src/` edit and no landed module edited.
- **Imports, unchanged** (read every signature first):
  - I83 (`test_imp8_probe_cylinder_source_closed_form`): `_g`, `_phi_cyl`, `_phi_cyl_dot`, `_GL_B`;
  - I80 (`test_imp8_probe_ladder_closed_form`): `RUNGS`, `LADDER_MM`, `_lab`;
  - I78 (`test_imp8_probe_f2182_closed_form`): `DUTY`, `N_GL`, `_axis_steady`, `_p`;
  - I58 (`test_imp8_entry48_lumped`): `_gl01`, `_tips_steady`, `BQ_BAND`;
  - step 0 (`test_imp8_rod_line_source_closed_form`): `L`, `MM`, `KAPPA`, `QUAD`, `T_PRINT`, `_phi`, `_phi_tips`;
  - I59 (`test_imp8_rod_line_source_duty_closed_form`): `SCALE`, `T_STAR`, `PERIOD`;
  - I62 (`test_imp8_rod_tips_duty_closed_form`): `_duhamel`, `_convolution`, `_phi_tips_dot`;
  - `ROD_RADIUS` (`tests.mesh.test_birdcage_implant_rod`);
  - `B1_BAND` (`test_imp4_pec_cylinder_closed_form`) and `INTERIOR_BAND` (`test_mat13_transient_perfusion_sphere`), by name.
- **`_bracket_cf(d, t, a)`.** I83's `_bracket` with y − x written as a²/((D + d)s), never D/s − d/s. The GL16 nodes go on [x, x + (y − x)]. I83's `_bracket` takes the difference (`:100–103`), and at a = 1e-10 m that loses every digit.
- **`_seg_mean(zeta, t, a, n=None)`.** `quad` over u ∈ [−1, 1] with `QUAD` when n is None; otherwise `_gl01(n)` mapped onto [−1, 1].
- **The twin functions.**
  - `_twin_a(r, z, t, a=ROD_RADIUS)` (assert r == 0.0) builds the sum above.
  - `_twin_a_dot(r, z, u, a=ROD_RADIUS)` and `_twin_a_steady(z, a)` are the matching closed forms.

**Gates (asserted, bands imported).** At every rung of `RUNGS`, on both signs:
- **(Lim)** |`_twin_a`(0, z, ∞) / `_twin_a_steady`(z) − 1| ≤ `B1_BAND`, at a = `ROD_RADIUS`: the u quadrature against the analytic antiderivative. Predicted at the round-off floor; I83's worst was 2.220e-16 (item 911 Outcome).
- **(Bq)** |GL_n / GL_2n − 1| ≤ `BQ_BAND`, n = `N_GL`, at t = `T_PRINT` and ∞. The integrand is analytic on [−1, 1], and the nearest branch point of √(a² + d²) sits a unit off the end face in u. Predicted ≤ 1e-14, against I80's worst 4.030e-14 (`…item894-W1.log:68`).
- **(Dut)** |`_duhamel`(`_twin_a`, 0, z, `T_STAR`, `DUTY`) − `_convolution`(`_twin_a_dot`, 0, z, `T_STAR`, `DUTY`)| / (`SCALE`·`_twin_a_steady`(z)) ≤ `B1_BAND`. The step sum runs on the quadrature route and the convolution on the closed-form Φ̇, so the two routes are independent. Precedents: I80's point twin 5.839e-17 (`…item894-W1.log:71`) and I83's cylinder worst 8.971e-16.
- **(Lim-a)** At a = `ROD_RADIUS`·1e-7 (1e-10 m), `_twin_a` against I80's two-tip line `_phi_tips`(0, z, t) at t = `T_PRINT` and ∞, and `_duhamel` of each at D = 0.5: |ratio − 1| ≤ `B1_BAND`.
  - The small-a deficit is O(a²), ≈ +(a/1 mm)²/12 on the near tip, so ≈ 8e-16 at 1e-10 m (drafter arithmetic). The `QUAD` floor (1e-13) sets the read.
  - Printed beside: the deficit at a = `ROD_RADIUS`/10, /20 and /40, and the observed order (expected ≈ 2; I83's line read 1.9973 / 1.9993, `…item911-W1.log:146`).

**Control (the weekly's; asserted at factor 1 on `INTERIOR_BAND`).** At 16 mm, continuous θ(900 s): 1 − `_twin_a`/`_phi_tips` > `INTERIOR_BAND`.
- **Ceiling:** 1 (θ > 0).
- **Predicted ≈ 0.044 at 900 s, ≈ 4.4× the band** (unlogged drafter arithmetic, against I80's logged point twin):
  - I80's point twin is 63.14 K/W (`…item894-W1.log:99`). That is `SCALE`·(L/2)·[erfc(1 mm/s)/1 mm + erfc(31 mm/s)/31 mm], with `SCALE` 4.420971 and κ 1.435407e-07 (`:112`), so s = 22.732 mm.
  - The cylinder's near mean is [g(2a) − g(0)]/(2a) − 2/(s√π) ≈ 0.95789 − 0.04964 per mm, against the point's 0.95039. That gives ≈ 60.35 K/W.
- **Steady:** 65.661 against 68.45374, predicted 0.0408, printed.
- **D = 0.5 is not asserted.** With the source off for the last 225 s, the 1/R parts cancel between the on and off edges. The radius term left is O(a²/s³), predicted ≈ 2e-3 … 5e-3, below the band. I83's uniform-source analogue is 1 − 1.973839/1.978384 = 2.3e-3 (`…item911-W1.log:119`, `…item894-W1.log:98`). Print it, labelled "below the band; not the weekly's control".
- **(Lim-a)'s own negative.** At a = `ROD_RADIUS`/10 the deficit is > `B1_BAND` (predicted ≈ 8e-4). Asserted as > `B1_BAND` only (rule (e)). It shows (Lim-a) is not comparing a kernel with itself.

**Printed (🟡 record, no band, ≤ 7 digits, K per watt per rod, infinite gel).** Head the block "placement and time quoted from memory; W3 confirms or corrects".
- **Per rung:** `SCALE`·Φ₂,a at 900 s continuous, at D = 0.5, and steady.
- **Per rung, the radius factor on the twin, Φ₂,a/Φ₂.** Beside it, I83's Φ_a/Φ on the line: 0.9532823 / 0.9811647 / 0.9889726 / 0.9941423 (`…item911-W1.log:119`, `:125`, `:130`, `:135`).
- **Per rung, the tip weighting at finite radius, twin_a / line_a**, with I83's `_phi_cyl` beyond the line's end. Beside it, I80's point twin / line: 6.543912 / 4.400661 / 3.60769 / 2.940895 (`…item894-W1.log:100`, `:103`, `:106`, `:109`).
- **The answer line:** radius factor against tip weighting at 16 mm (predicted ≈ 4 % against × 6.5, so the tip weighting dominates), printed from the computed values.
- **Sensitivity, h = a inside the rod's end:** Φ₂,a/Φ₂, predicted ≈ 0.60 at 900 s and ≈ 0.63 steady. Label it "the cylinder's length is an execution detail; the weekly's to fix".
- |θ(+z) − θ(−z)| (round-off, by symmetry).

**Window.** Before it, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. **W1** `IMP-8-item934-W1`, smoke, real, `-n 2`. Make one bare foreground harness call, Bash timeout 660000 ms; I83's and I80's modules ride along as the comparands' regressions:

`scripts/testing/run_and_log.sh IMP-8-item934-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/validation/test_imp8_probe_twin_cylinder_closed_form.py tests/validation/test_imp8_probe_cylinder_source_closed_form.py tests/validation/test_imp8_probe_ladder_closed_form.py -v -s --tb=short'"`

**Price.** I83's window, with I80 riding along, ran 14 tests in 2.15 / 2.21 s, 4 s harness (item 911 Outcome, `20261109T110006Z_IMP-8-item911-W1.log`). The new module adds two segment `quad`s per kernel call, ≈ 4 kernel calls per step sum, and the (Lim-a) ladder (≈ 7 radii × 4 rungs × 3 reads). Predicted ≤ 30 s against the 120 s ceiling.

**Tier.** Smoke, numpy, `-n 2`. **Class:** lane-3 closed form → chain step. **Cost:** ≈ 4.8 measured / 15 predicted min.

**Traps already paid for** (items 874 / 894 / 911, plus this module's own).
- **Small-a cancellations, this module's own:**
  - D/s − d/s: use a²/((D + d)s);
  - g(ζ + a) − g(ζ − a) for a ≪ ζ: use the antiderivative at a = `ROD_RADIUS` only; (Lim-a)'s steady read goes through the quadrature;
  - the erf difference in `_phi_cyl_dot` for h ≪ s: use Φ̇ at a = `ROD_RADIUS` only; (Lim-a)'s D = 0.5 read goes through `_duhamel`, which needs only the kernel.
- **`quad`'s absolute tolerance.** `QUAD`'s `epsabs` of 1e-13, against a z′ integral of size h·B ≈ 2e-10 m · 1e3 m⁻¹ at a = 1e-10 m, is a 5e-7 relative error. Integrate the segment **mean** in u, never in z′.
- I83's `_phi_cyl` asserts |z| > L/2 strictly, so it cannot read the end face. Import it for twin_a / line_a beyond the line's end only.
- `_phi_steady` and I58's `_line(·, ·, ∞)` divide by r. Never call them on the axis.
- `_step` guards t ≤ 0, and `_phi_cyl_dot` returns 0 for u ≤ 0. At t = 900 s, both routes exclude the edge at 900 s.
- Prints go through I78's rank-0 `_p`. Use `-s`, and print ≤ 7 digits.
- **The leak guard:** a hook refusal on a printed θ line moves the log unchanged to `logs/leak-held/`. Never `--no-verify`.

**Scope.** Infinite gel, two uniform tip cylinders (the length fixed above), 1 W per rod, four placements on the axis. **Not claimed:**
- I56's FEM balls (a ball is not a cylinder);
- the box walls or the rod's own conduction (I75 / I81);
- a ball-averaged read;
- any F2182 figure (W3).

The rungs are placements, not a convergence ladder.

**Negative result.**
- (Lim), (Bq), (Dut) or (Lim-a) red: a derivation or quadrature error. Stop, BLOCKED, band unmoved, every value printed.
- **The control red** (the cylinder twin within 1 % of the point twin at 16 mm, continuous) contradicts the arithmetic above. Report both values and stop. Never re-fix the cylinder's length to rescue it.
- No `PREDICTED RED:` line.

**Status it can move.**
- `IMP-8` gains step 21, a 🟡 record line: the finite-radius two-tip twin on the axis, gated by (Lim) / (Bq) / (Dut) / (Lim-a), with the radius-against-tip-weighting answer at the probe.
- §10 I86 is DONE.
- I76's lumped column can carry Φ₂,a beside I78's point twin; whether it does is the weekly's call.

**Depends on:** nothing. Items 874, 894 and 911 and their logs are on `main`.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 934 — `IMP-8` step 21 (§10 I86) — per `docs/planning/chunks/IMP-8.md` "§9 item 934". Create only `tests/validation/test_imp8_probe_twin_cylinder_closed_form.py`; edit no landed module and no `src/` file; import module objects only, never `test_*`. Numpy / scipy, no mesh, no solve. Before W1, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; run pytest with -s; container `timeout -k 30`. Never move a band, and never change the cylinder's length (h = 2a, centred on ±L/2) to rescue a red. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-8-step21-*`, mark item 934 BLOCKED in the same commit, and report."

**Outcome (2026-11-10, 06:00 slot): DONE.** New `tests/validation/test_imp8_probe_twin_cylinder_closed_form.py`; no `src/` edit, no landed module edited. W1 `20261110T120247Z_IMP-8-item934-W1.log` — smoke, real, `-n 2`, container `timeout -k 30 120`: 20 passed on each rank (40 `PASSED` tokens, 0 `FAILED`; the new module 6, I83's and I80's modules riding along green), 3.33 / 3.35 s pytest, 5 s harness. `pgrep -c python3` read 0 before the window.
- **(Lim)** worst 4.441e-16 (`:59`), both signs at every rung, against `B1_BAND`. **(Bq)** worst 2.442e-15 (`:81`) against `BQ_BAND`. **(Dut)** worst 4.058e-16 (`:92`) against `B1_BAND`. **(Lim-a)** at a = 1e-10 m, t = 900 s / ∞ / D = 0.5, worst 4.441e-15 (`:120`) against `B1_BAND`.
- **(Lim-a) ladder (printed) and own negative (asserted > `B1_BAND`).** Φ₂,a/Φ₂ − 1 at 16 mm, 900 s, at a/10, /20, /40: +8.5713e-04 / +2.1774e-04 / +5.4651e-05, observed order 1.9769 / 1.9943 (`:121`); the a/10 deficit 8.571e-04 > `B1_BAND` (`:122`), as predicted (≈ 8e-4). The sign is an excess, the cylinder's quadrupole term +a²/(12ζ³) in the small-a regime.
- **Control (asserted at factor 1 on `INTERIOR_BAND`).** At 16 mm, continuous θ(900 s): cylinder twin 60.34935 against point twin 63.14 K/W, 1 − ratio 4.4198e-02 (`:124`), predicted ≈ 0.044. Steady 65.66113 against 68.45374, 4.0796e-02 (predicted 0.0408, `:125`). D = 0.5, printed below the band: 2.488589 against 2.498341, 3.9031e-03 (`:126`).
- **Record (🟡).** SCALE·Φ₂,a at 900 s / D = 0.5 / steady: 60.34935 / 2.488589 / 65.66113 at 16 mm, 30.29385 / 2.441764 / 35.55115 at 17 mm, 19.07042 / 2.374175 / 24.27028 at 18 mm, 10.12165 / 2.18383 / 15.19819 at 20 mm; |θ(+z) − θ(−z)| 0 everywhere (`:129–142`). Radius factor Φ₂,a/Φ₂ at 900 s 0.9558023 / 1.010793 / 1.008342 / 1.004202 (I83's line 0.9532823 / 0.9811647 / 0.9889726 / 0.9941423). Tip weighting twin_a/line_a at 900 s 6.561211 / 4.533547 / 3.678348 / 2.970653 (point twin/line 6.543912 / 4.400661 / 3.60769 / 2.940895, I80's logged strings reproduced).
- **Answer (`:143`).** At 16 mm the radius moves the twin's read by 4.42 % (factor 0.9558023), against a tip weighting of × 6.561211: **the tip weighting dominates**, as predicted. Not predicted, and recorded: from 17 mm on, the radius factor is above 1 (+1.08 % / +0.83 % / +0.42 %): once the probe is a radius clear of the cylinder's face, the h = 2a cylinder's quadrupole excess outweighs the near face's distance deficit. The line's I83 factor stays below 1 at every rung.
- **Sensitivity (the cylinder's length is an execution detail; the weekly's to fix).** h = a inside the rod's end: Φ₂,a/Φ₂ at 900 s / steady 0.6013158 / 0.6325601 at 16 mm (predicted ≈ 0.60 / 0.63), 0.7560154 / 0.7929387, 0.8203262 / 0.859705, 0.8755296 / 0.9178775 at 17 / 18 / 20 mm (`:145–148`). The twin's read at the probe is therefore much more sensitive to the tip's length than to its radius.


## §7 row moved verbatim 2026-11-11 (rotate_plan_archive.py chunks)

 The heated rod in gel: the finite line source in closed form (Phase 7; opened 2026-11-02 10:30 daily review, rule (5) lane 3; the weekly ratifies) | ✅ on the closed form (ratified by the 11-02 weekly); the FEM half and the rod's own thermal properties are 🟡 records in the chunk file — I55, I56, I58, I60, I62, I63, I64, I65, I66, I69, I70, I71 and I73 (steps 9–15, every gate at `B1_BAND` / `INTERIOR_BAND`); step 16 (I78, §9 item 874, 2026-11-08) is a 🟡 record — the lumped F2182-shaped prediction on the axis 1 mm beyond each tip at 900 s, (Lim) / (Bq) / (Dut) green at `B1_BAND` / `BQ_BAND`, θ per watt 9.648662 K line / 63.14 K twin continuous, 1.978384 / 2.498341 K at D = 0.5 (`20261108T133157Z_IMP-8-item874-W1.log`); step 17 (I75, §9 item 881, 2026-11-08) is a 🟡 record — the F2182-shaped probe on the gel / Ti / hole bracket, θ_probe(900 s) per watt Ti 10.57 / gel 9.162 / hole 8.783 K, D = 0.5 by superposition (`20261108T200636Z_IMP-8-item881-W3.log`); step 18 (I80, §9 item 894, 2026-11-08) is a 🟡 record — the lumped probe-placement ladder at 16 / 17 / 18 / 20 mm on the axis, (Lim) / (Bq) / (Dut) green at every rung, I78's record reproduced at rungs 1–2, θ(900 s) per watt line 9.648662 / 6.810428 / 5.242317 / 3.427288 K (`20261109T043829Z_IMP-8-item894-W1.log`); step 19 (I81, §9 item 895, 2026-11-09) is a 🟡 record — I75's bracket read at the clinical placement and on a 17 / 18 / 20 mm ladder, (R75) and I75's gates green on every arm, θ(900 s) per watt at 17 mm Ti 7.454 / gel 6.662 / hole 6.474 K (`20261109T061832Z_IMP-8-item895-W3.log`); step 21 (I86, §9 item 934, 2026-11-10) is a 🟡 record — the finite-radius two-tip twin on the axis (tips as uniform cylinders a = h/2 = `ROD_RADIUS` on ±L/2), (Lim) / (Bq) / (Dut) / (Lim-a) green at every rung, θ(900 s) per watt 60.34935 K at 16 mm against the point twin's 63.14 (control 1 − ratio 0.0442 > `INTERIOR_BAND`); at 16 mm the radius factor 0.9558 against the tip weighting × 6.561 — the tip weighting dominates (`20261110T120247Z_IMP-8-item934-W1.log`); with I61 sized as `EX-182` (`imp:23`). *History: `docs/planning/chunks/IMP-8.md`.* 
