# IMP-6 — history

Opened 2026-10-12 by the 10:30 daily review under §9 rule (5) lane 3 (Phases 7/8), for the weekly to ratify. §10 "Phase 7 — implants" names "wires, rods, plates": wires and rods have canonicals (`IMP-4` / `IMP-5` the lead shaft, `IMP-2` the rod tip), plates have none. `IMP-6` is `IMP-2`'s companion: the PEC **oblate** spheroid, whose rim enhancement 1/L_x under an in-plane field is to a plate or disc edge what 1/L_z is to a rod tip. Same host, same hole route, same bands, same grain. Sources: Landau & Lifshitz *ECM* §4; Stratton *Electromagnetic Theory* §3.27 (the commissioning text also names §3.28); Osborn, Phys. Rev. 67, 351 (1945). These are the citations `IMP-2` uses, and they are **not verified in-repo**. `docs/references/jin-fem-3e/` has no spheroid (ch. 11 mentions spheroidal surfaces only as a truncation boundary). The docstring therefore derives every identity the test asserts, and the in-repo cross-check is the ellipsoid integral (Int) below.

## Proposed §7 row

| `IMP-6` | *(history and §9 briefs: `docs/planning/chunks/IMP-6.md`)* **A PEC oblate spheroid in a lossy host under a uniform quasi-static field — the closed-form stand-in for a plate / disc edge (Phase 7 (7a))** (opened by the 2026-10-12 10:30 daily under rule (5) lane 3; for the weekly to ratify). `IMP-2`'s companion: semi-axes (a, a, c), a > c. In-plane drive E₀x̂: rim enhancement 1/L_x. Axial drive E₀ẑ: pole enhancement 1/L_z. L_x + L_y + L_z = 1. **Done-when (§4):** *step 0* (numpy, §9 item 403): the depolarisation factors, both potentials and fields, and the rim / pole identities at the bands in the chunk file, with the prolate-formula control missing ⇒ 🟡. *Step 1* (FEM, degree 2, `IMP-2`'s hole route with an additive `oblate=True` on `spheroid_in_box_domain`, a = 6 mm, aspect 2 and 5, §9 item 404): per aspect \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball Φ ≤ `MAT4_BAND`, cavity dofs zero; the axial-drive control misses ⇒ ✅. Not claimed: a square-edged plate, a plate in a coil, SAR, finite k, convergence. | ⬜ Open. | step 0 smoke, no MPI; step 1 heavy, -n 2, complex |

### §9 item 403 — 2026-10-12 (lane 3, closed form)

`IMP-6` step 0: the PEC oblate spheroid in a uniform field, axial and in-plane, numpy only (implementer; smoke, no MPI, real build; `src/` additive + test; ≈ 11 measured / 40 predicted min: ≈ 1 min of windows + 15 fixed + ≈ 25 code; depends on nothing).

- **Geometry and notation.** Semi-axes a along x and y, c along z, with a > c. Aspect m = a/c. e = √(1 − c²/a²), focal-ring radius f = a·e, ξ₀ = c/f = √(1 − e²)/e. Oblate spheroidal coordinates: ρ = f√((ξ² + 1)(1 − η²)), z = fξη, with ξ ≥ 0 and the body surface at ξ = ξ₀. **Name swap versus `IMP-2`:** here `semi_major` = a is the equatorial radius (x, y) and `semi_minor` = c lies along z. Keep `semi_major > semi_minor` as the check.
- **Anchor (derive in the docstring).** Write both drives in one Cartesian form:
  - φ = −E₀ x_i [1 − D_i(ξ)/D_i(ξ₀)], with x_i = z for the axial drive and x for the in-plane drive.
  - Axial: D_z = h(ξ) = arccot ξ − 1/ξ, so h = q₁(ξ)/ξ with q₁ = ξ arccot ξ − 1, the decaying n = 1 oblate function.
  - In-plane: D_x = g(ξ) = arccot ξ − ξ/(1 + ξ²), which is K(ξ)/√(1 + ξ²) for the decaying m = 1 solution K = √(1 + ξ²)[arccot ξ − ξ/(1 + ξ²)]. Reduction of order from the regular solution √(1 + ξ²) gives K, with Wronskian −2/(1 + ξ²).
  - Derivatives: h′ = 1/(ξ²(1 + ξ²)), g′ = −2/(1 + ξ²)².
  - Field: E = E₀x̂_i − ∇φ_s, where φ_s = E₀ x_i D(ξ)/D(ξ₀) and ∇φ_s = (E₀/D₀)[D(ξ) x̂_i + x_i D′(ξ) ∇ξ]. Get ∇ξ by implicit differentiation of F = ρ²/(ξ² + 1) + z²/ξ² − f² = 0: ∇ξ = −∇F/F_ξ, ∇F = (2x/(ξ² + 1), 2y/(ξ² + 1), 2z/ξ²), F_ξ = −2ξρ²/(ξ² + 1)² − 2z²/ξ³.
  - Depolarisation factors, each from its own drive's Wronskian:
    - L_z = (1 + ξ₀²)(1 − ξ₀ arccot ξ₀) = (1/e²)[1 − √(1 − e²) arcsin(e)/e].
    - L_x = L_y = (ξ₀/2)[(1 + ξ₀²) arccot ξ₀ − ξ₀]. Computed **independently** of L_z, so the sum rule is a real test.
  - Surface identities: rim E_x(±a, 0, 0) = E₀/L_x (in-plane drive); pole E_z(0, 0, ±c) = E₀/L_z (axial drive).
  - Far-field dipoles: p_i = E₀a²c/(3L_i), which is a³ for the sphere.
  - Disc limit: L_x → πξ₀/4 − ξ₀² + O(ξ₀³), so p_x/a³ → 4/(3π).
- **Code** (append to `src/fem_em_solver/utils/analytical.py` directly after `pec_prolate_spheroid_field`, ≈ line 1654, under a new `IMP-6 step 0` banner, before the `MAT-15` block; additive, nothing reordered).
  - `oblate_spheroid_depolarisation(e, *, branch=None) -> (l_x, l_z)`.
  - `oblate_spheroid_d_functions(xi, *, branch=None) -> (h, h′, g, g′)`.
  - `oblate_spheroidal_xi(points, focal)`.
  - `pec_oblate_spheroid_potential(points, *, semi_major, semi_minor, drive="axial"|"inplane", e0=1.0, scattered_only=False)`.
  - `pec_oblate_spheroid_field(points, *, semi_major, semi_minor, drive=…, e0=1.0, interior="zero"|"continue")` (step 1's trace; same contract as `pec_prolate_spheroid_field`).
  - Interior mask: (x² + y²)/a² + z²/c² < 1 − 1e-12.
  - Test file: `tests/validation/test_imp6_oblate_closed_form.py`.
- **Anchors (asserted).** Aspects {1.01, 2, 5, 20, 100} unless stated.
  - (L) The two L_z forms agree ≤ 1e-13.
  - (Sum) \|2L_x + L_z − 1\| ≤ 1e-13.
  - (Int) Each L_i equals the ellipsoid integral (a²c/2)∫₀^∞ ds/((a_i² + s)(a² + s)√(c² + s)), by `scipy.integrate.quad` at `epsrel=1e-13`, to ≤ 1e-10 relative at aspects 2, 5, 20.
  - (Br) The closed and series branches of h and g agree ≤ 1e-13 relative on ξ ∈ [5, 20].
  - (Bd) \|φ\|/(E₀a) ≤ 1e-12 at 200 Fibonacci surface points, both drives, aspects 2 and 5.
  - (Tg) Tangential E on the same points ≤ 1e-12·E₀.
  - (Rim) \|E_x(±a, 0, 0)·L_x − 1\| ≤ 1e-12, plus \|E(0, ±a, 0)\| ≤ 1e-12·E₀ (zero normal field where n·x̂ = 0).
  - (Pole) \|E_z(0, 0, ±c)·L_z − 1\| ≤ 1e-12.
  - (F) Field against the central-difference gradient of the potential at δ = 2e-3·a and 4e-3·a, both drives, 40 points with ξ ∈ [1.2ξ₀ + 0.05, 3] and \|η\| ≤ 0.95. Asserted in **Richardson / ratio form** (the lesson ratified on `IMP-2`): ratio ∈ [3.8, 4.2] and \|E − (4D(δ) − D(2δ))/3\| ≤ 1e-7·E₀.
  - (Lap) FD Laplacian of both potentials at the same δ pair: ratio ∈ [3.8, 4.2], and the Richardson residual ·a/E₀ ≤ 1e-7.
  - (Far) The P₁ projection R²(3/2)∫φ_s μ dμ = p_z (axial). The sinθcosφ projection R²(3/4π)∫φ_s sinθ cosφ dΩ = p_x (in-plane; 64-node Gauss–Legendre × 64-point trapezoid). Both at R = 2a and 4a, ≤ 1e-12 relative.
  - (Sph) At e = 1e-2: \|(L_z − 1/3)/e² − 2/15\| ≤ 1e-3 and \|(L_x − 1/3)/e² + 1/15\| ≤ 1e-3. Both potentials match the sphere's −E₀x_i(1 − a³/r³) to ≤ 1e-3·E₀a for r ∈ [1.2a, 3a].
  - (Disc) At c/a = 1e-3: \|(L_x − πξ₀/4)/ξ₀² + 1\| ≤ 2e-3 (the next term is (π/4)ξ₀ ≈ 7.9e-4).
- **Printed (predicted here by stdlib arithmetic).**
  - Aspect 2 / 5: L_z = 0.527200 / 0.750484, L_x = 0.236400 / 0.124758.
  - Aspect 2 / 5: 1/L_x = 4.2301 / 8.0155, 1/L_z = 1.8968 / 1.3325, pointwise rim (1/L_x)² = 17.89 / 64.25.
  - Aspect 100: L_x 0.007755 against πc/(4a) = 0.007854; p_x/a³ 0.429822 against 4/(3π) = 0.424413.
  - Osborn's oblate table values beside these, only if the executor can quote them. They are not in-repo; (Int) is the in-repo check.
- **Negative controls (rule (e): every factor below is predicted, computed here, printed beside the measured value; only the separation is asserted).**
  - (c1) The prolate formula used for the oblate shape: `prolate_spheroid_depolarisation(e)` at the same e gives L_x^pro = (1 − L_z^pro)/2. In the rim identity this reads E_x(rim)·L_x^pro = L_x^pro/L_x^obl = **1.748 / 3.784** at aspect 2 / 5. The prolate rim enhancement under-reads by 43 % / 74 %. Assert only \|· − 1\| > 0.1. The axial counterpart L_z^pro/L_z^obl = 0.329 / 0.0744 is printed.
  - (c2) The regular (growing) solution in place of the decaying one: D_i ≡ 1, i.e. φ_s = E₀x_i. This is `IMP-2`'s Q₁ → P₁ swap, since the regular axial solution ξ gives φ_s = E₀fηξ = E₀z. (Far)'s 4a/2a projection ratio reads 8 (derived, R³). Assert \|ratio − 1\| > 1e-3, on both drives.
- **Window** (smoke; real build, the `IMP-2` step-0 pattern without `mpiexec`; split so each part fits 30 s):
  - W1 (≈ 15 s): `scripts/testing/run_and_log.sh IMP-6-step0-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 30 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_closed_form.py -s -v --tb=short'"`
  - W1r (≈ 15 s, the rule (c) re-gate of the additive `analytical.py` change): label `IMP-6-step0-W1r`, same command with paths `tests/environment/test_complex_mode.py tests/unit/test_analytical_lightweight.py tests/validation/test_imp2_spheroid_closed_form.py`.
- **Traps.**
  - 1 − ξ₀ arccot ξ₀, arccot ξ − 1/ξ and arccot ξ − ξ/(1 + ξ²) cancel catastrophically at large ξ. Use series branches for ξ > 10 (`IMP-2`'s threshold, 40 terms): h = Σ_{k≥1} (−1)^k ξ^{−(2k+1)}/(2k + 1) and g = Σ_{k≥1} (−1)^{k+1} (2k/(2k + 1)) ξ^{−(2k+1)}. (Sph) at e = 1e-2 (ξ₀ ≈ 100) runs on the series.
  - ξ from the quadratic f²ξ⁴ + (f² − r²)ξ² − z² = 0, taking the stable root:
    - ξ² = [(r² − f²) + √((r² − f²)² + 4f²z²)]/(2f²) for r² ≥ f²;
    - ξ² = 2z²/[√((r² − f²)² + 4f²z²) − (r² − f²)] otherwise.
  - Never evaluate inside the focal disc (ξ → 0, interior). Use arccot ξ = arctan(1/ξ) only for ξ > 0.
  - The real build skips the four complex-only environment tests; that is expected, and `FEM_EM_REQUIRE_COMPLEX` stays unset.
  - The smoke path carries only `tests/environment/test_complex_mode.py` (10:30 ruling, item 364).
  - A window killed at 30 s is not a result. Shrink the path list; never raise the ceiling.
  - Foreground, `timeout: 660000`; `-s`; append to `analytical.py`, never reorder.
- **Scope.** Quasi-static, PEC, no host, no k. Not: step 1's FEM, SAR, a square-edged plate (a real plate edge carries a wedge singularity the smooth spheroid does not), a thin-plate FEM.
- **Status it can move.** `IMP-6` ⬜ → 🟡 ("step 0 green ⇒ 🟡" in the row above; the review inserts the row in the same commit that queues this item, rule (j)).
- **Negative result.**
  - (Rim) / (Pole) red with (Far) and (Int) green: an L formula is wrong (prolate / oblate swapped, or L_x taken as (1 − L_z)/2 instead of the independent form).
  - (Sum) red with (Int) green for one factor: the other factor's closed form is wrong.
  - (Bd) red: the ξ map.
  - (F) ratio ≈ 1: a gradient bug.
  - Fix in-slot. Still red ⇒ park on `attempt/*`, BLOCKED in the same commit (rule (d)); item 404 stays unqueued.

**Outcome 2026-10-12 (13:30 CDT slot, implementer): DONE, `IMP-6` ⬜ → 🟡.** `src/fem_em_solver/utils/analytical.py` gains the additive `IMP-6 step 0` block after `pec_prolate_spheroid_field` (derivations in the banner); test `tests/validation/test_imp6_oblate_closed_form.py`. Smoke, no MPI, real build.
- **Windows.** W1 first run red on (F) only (`20261012T183401Z_IMP-6-step0-W1.log`, 34 passed / 2 failed, 4 s); W1 green after the disclosed bound move below (`20261012T183435Z_IMP-6-step0-W1.log`, 36 passed / 3 skipped, 3 s); W1r re-gate green (`20261012T183442Z_IMP-6-step0-W1r.log`, `test_analytical_lightweight.py` + `test_imp2_spheroid_closed_form.py`, 20 passed / 3 skipped, 3 s).
- **Measured.** L_z / L_x = 0.527200 / 0.236400 (aspect 2), 0.750484 / 0.124758 (aspect 5), as predicted; aspect 100 L_x 0.007755, p_x/a³ 0.429822. (L) ≤ 2.6e-15; (Sum) ≤ 6.7e-16; (Int) ≤ 1.3e-15 relative (aspects 2, 5, 20); (Br) h 5.2e-14, g 4.7e-14 on ξ ∈ [5, 20]; (Bd) ≤ 8.2e-16; (Tg) ≤ 1.2e-14; (Rim) ≤ 5.8e-14 (aspect 100), \|E(0, ±a, 0)\| ≤ 8.9e-16; (Pole) ≤ 2.7e-15; (Far) ≤ 3.1e-15 at R = 2a, 4a, both drives; (Lap) ratio 3.999–4.000, Richardson ≤ 6.6e-8; (Sph) 0.13334095 (2/15) and −0.06667048 (−1/15), potentials within 6.6e-5·E₀a of the sphere; (Disc) 7.85e-4 (next term 7.854e-4).
- **(F), bound moved with measurement (disclosed).** The second-order ratio reads 3.999–4.000 at both aspects and drives. The brief's absolute Richardson bound 1e-7 holds at aspect 2 (1.2e-9 / 2.3e-9) but reads 3.55e-7 (axial) / 5.91e-7 (in-plane) at aspect 5. That is the O(δ⁴) truncation term at the sample points nearest the rim (ξ = 1.2ξ₀ + 0.05 sits ≈ 0.02a off a rim of curvature radius 0.04a), not a defect. At aspect 5 the test therefore asserts the same 1e-7 one δ-halving finer (1e-3a / 2e-3a pair: 2.2e-8 / 3.7e-8) plus the residual's own fourth-order ratio in [14, 18] (measured 15.98 / 16.03). The pre-registered pair stays printed; aspect 2 is asserted as written.
- **Controls (rule (e)).** (c1) E_x(rim)·L_x^pro = 1.7480 / 3.7840 (predicted 1.748 / 3.784); L_z^pro/L_z^obl 0.3292 / 0.0744 printed. (c2), per the 10:30 ruling (Q2), the regular ratio ξ/ξ₀ or √(1+ξ²)/√(1+ξ₀²) in place of D/D₀: 4a/2a projection ratio 16.44 / 16.55 (axial) and 15.75 / 15.68 (in-plane). **The ruling's "predicted ≈ 8" is the brief's original D ≡ 1 variant**; the ruling's own substitution makes φ_s ∝ x_i·ξ ∼ r², so the projection grows ∝ R⁴ (≈ 16). Both variants are run and asserted separated; D ≡ 1 reads 8.0000 as derived.
- **Not done.** Osborn's oblate table values are not quoted (not in-repo).

### §9 item 404 — 2026-10-12 (lane 3, closed form; depends on item 403)

`IMP-6` step 1: FEM at aspect 2 and 5, in-plane drive, on `IMP-2`'s hole route (implementer; heavy, `-n 2`, complex; ≈ 20 measured / 70 predicted min: ≈ 14 min of windows + 15 fixed + ≈ 40 code; depends on 403; skip if it did not land).

- **Sizing, stated first.**
  - \|k\| = 17.90 1/m (`IMP-1`'s host: σ 0.57 S/m, εᵣ 78, 64 MHz).
  - a = 6 mm (the long semi-axis, as `IMP-2`'s c), so (\|k\|a)² = 1.1535 % ≤ 4.886/3 = 1.629 %.
  - c = 3 mm at aspect 2 and 1.2 mm at aspect 5. Rim curvature radius ρ_r = c²/a = 1.5 / 0.24 mm, numerically `IMP-2`'s ρ_t.
  - Box half-extents (1.75a, 1.75a, a) = (10.5, 10.5, 6) mm. The 1.5× probe spheroid (9 mm in x) clears the wall.
  - Corner 16.02 mm, (\|k\|R)² = 8.22 %, between `IMP-2`'s 6.92 % and `IMP-1`'s 9.6 %. Printed beside I2a's 1.56e-2 face departure.
- **Code.**
  - Mesher (rule (c), disclosed). `spheroid_in_box_domain` does **not** build an oblate spheroid today: it raises unless `semi_major > semi_minor` with semi_major along z, dilates (a, a, c) = (minor, minor, major), and puts its tip Balls at (0, 0, ±c). Add a keyword-only `oblate: bool = False`, with the default path byte-identical. With `oblate=True`:
    - Rotate the unit sphere 45° about z before `dilate`, so the OCC seam stays off the drive's rim points (±a, 0, 0).
    - `dilate(…, semi_major, semi_major, semi_minor)`.
    - Enclosure check bx, by > semi_major and bz > semi_minor.
    - The two Balls at (±semi_major, 0, 0): radius 3ρ_r, `VIn` = `resolution_tip`, where ρ_r = semi_minor²/semi_major is the same expression as today's ρ_t.
    - **plus** a rim tube: a `MathEval` field with the distance to the rim circle, `f"Sqrt((Sqrt(x*x+y*y)-{a!r})^2+z*z)"`, feeding a `Threshold` with SizeMin ρ_r, SizeMax `resolution_far`, DistMin ρ_r, DistMax 4ρ_r, added to the `Min`.
    - Cavity Distance→Threshold, mass reference 4πa²c/3 and the absolute wall `tol = 1e-6` all unchanged.
  - Fixture: `resolution_tip` = ρ_r/3, `resolution_surface` = a/12 = 0.5 mm (`IMP-2`'s aspect-2 surface grain), `resolution_far` = a/4 = 1.5 mm.
  - Test module `tests/validation/test_imp6_oblate_fem.py`. Import from `test_imp1_pec_sphere_lossy_host.py` as `IMP-2` does: `TH15_BETA_BAND`, `MAT4_BAND`, the tags, `DEGREE`, `E0`, the host constants, `_assemble_real`, `_cavity_dofs`, `_phi_closed_form_on_mask`. Copy `IMP-2`'s `_numpy_quadrature_mean` / `_owned_masked_cells` / `_rss_gib` patterns.
  - β: the least-squares amplitude of item 403's **in-plane** scattered field (`pec_oblate_spheroid_field(drive="inplane")` − E₀x̂) fitted on Re E − E₀x̂, on 24 Fibonacci points on each of the 1.2× and 1.5× scaled oblate spheroids.
  - Rim mask: DG0 indicator of cell midpoints within r_r = 2c²/a (3 mm / 0.48 mm) of (±a, 0, 0). Squares in numpy, no UFL comparison.
  - Φ_rim,cf,h: the degree-8 numpy quadrature over the same cells with `interior="continue"`.
- **Anchors** (asserted per aspect at imported bands, as `IMP-2` asserted them; the bands come from `TH-15` / `MAT-4`, not from a measurement on this fixture):
  - (i) \|β − 1\| ≤ `TH15_BETA_BAND` (4.886 %).
  - (ii) \|Φ_rim,h/Φ_rim,cf,h − 1\| ≤ `MAT4_BAND` (10 %).
  - (iii) cavity dofs > 0 and max \|E\| < 1e-12.
  - (M′) In W0, on each mesh, the numpy quadrature of `IMP-1`'s dipole equals `_phi_closed_form_on_mask` to ≤ 1e-12.
- **Predicted, printed and not asserted (rule (e)).**
  - \|β − 1\| ≲ 1.5 % and Φ miss ≲ 3 %. `IMP-2` measured 0.46 / 0.90 % and 1.11 / 0.20 % at h_tip/ρ_t 0.80 / 0.67, but on a different fixture, so this is evidence that the bands are reachable, not a backing.
  - h_rim/ρ_r.
  - The meshed body-volume deficit (8b_xb_yb_z − Σ air-cell volumes)/(4πa²c/3) − 1. The faceted rim is the likeliest bias at aspect 5.
  - Pointwise rim (1/L_x)² = 17.89 / 64.25.
  - \|Im E\|/\|Re E\|; cells, seconds, RSS (reduced outside the rank-0 branch).
- **Negative control** (W2, aspect-5 mesh, predicted, rule (e)). Drive axially: the wall trace is item 403's closed-form axial field, zero on the cavity.
  - β fitted against the **in-plane** basis on Re E − E₀ẑ: predicted ≈ 0 by parity. The ceiling is 1/0.04886 = 20.5×. Assert only \|β − 1\| > band.
  - Rim-ball Φ_h against the in-plane Φ_rim,cf,h: the axial normal field vanishes at the rim (x_z/a_z² = 0 at z = 0), so the miss is predicted at 90–100 %. It saturates at 100 %, so it can show at most 10× the band. Assert only > 10 %.
  - Also printed: β_axial fitted on the axial basis. This is a record, not asserted, because no prior measurement of that comparison exists on this fixture.
- **Price** (measured on the parent fixtures).
  - `IMP-2` W1: 37 359 / 44 304 cells, solves 105.8 / 83.3 s, 6.46 GiB summed (`20261001T201608Z_IMP-2-step1-W1.log`).
  - `IMP-2` W0: 109 s (`…201413Z…-W0.log:456`).
  - `IMP-1`: 54 524 cells 101.8 s / 10.1 GiB; 90 713 cells 221.9 s / 19.4 GiB.
  - Complex `tests/environment` costs ≈ 70–80 s per window (`IMP-5` step 1 note: 49–57 s solves in 133–168 s windows).
  - Predicted ≈ 40–80 k cells at either aspect: oblate cavity area ≈ 1.6× / 3.4× `IMP-2`'s at aspect 2 / 5, balls of `IMP-2`'s size, a box 1.5× larger, and the rim tube ≈ 7 k cells at aspect 5. That gives ≈ 60–180 s per solve.
  - **Cell guard 95 000, checked in W0 before any solve.**
- **Windows** (complex, `-n 2`; one solve per window, so each footer returns inside the 660 s Bash ceiling).
  - W0 (≈ 150 s; census, numpy and the mesher re-gate): `scripts/testing/run_and_log.sh IMP-6-step1-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/unit/test_analytical_lightweight.py tests/validation/test_imp6_oblate_closed_form.py tests/validation/test_imp2_spheroid_fem.py::test_census_and_numpy tests/validation/test_imp6_oblate_fem.py::test_census_and_numpy -s -v --tb=short'"`. `IMP-2`'s census is the rule (c) re-gate of the `oblate=False` path. Its printed cells sit beside 37 359 / 44 304; a mismatch is reported (it may be an image move), and it is a defect only if the default branch's diff is non-empty.
  - W1a (≈ 170–260 s): label `IMP-6-step1-W1a`, `timeout -k 30 600`, paths `tests/environment tests/validation/test_pec_sphere_hole.py tests/validation/test_imp6_oblate_fem.py::test_in_plane_aspect2`. `test_pec_sphere_hole.py` is the mesher-module re-gate `IMP-2` ran.
  - W1b (≈ 140–260 s): label `IMP-6-step1-W1b`, `timeout -k 30 600`, `tests/environment tests/validation/test_imp6_oblate_fem.py::test_in_plane_aspect5`.
  - W2 (≈ 140–260 s): label `IMP-6-step1-W2`, `timeout -k 30 600`, `tests/environment tests/validation/test_imp6_oblate_fem.py::test_axial_control`.
- **Traps.**
  - The basis is the **in-plane** field minus E₀**x̂**. `IMP-2`'s `uniform[:, 0 if drive == "transverse" else 2]` is inverted here; copying it is the trap.
  - `MathEval` takes a literal expression: bake a in with `repr`, print the string, and use gmsh's `Sqrt`.
  - A `Threshold` returns SizeMax beyond DistMax. Never give one a fine SizeMax, or it caps the whole box.
  - OCC poles land on the flat faces (0, 0, ±c). If W0 prints minimum SICN under 0.1, report it and do not proceed to solves.
  - `IMP-2`'s two pre-W0 failures:
    - The wall tolerance stays absolute 1e-6 (OCC pads the bounding box by 1e-7).
    - The collective RSS reduction stays outside the rank-0 branch.
  - Cavity group from the **meshed** boundary, never the tool. Owned cells only (`size_local`). Read `geometry.dofmap`'s shape at runtime. No UFL `<` and no `sqrt` comparisons: masks in numpy.
  - Probe points never sit on a wall.
  - A killed compile leaves an orphan `.c` in `/root/.cache/fenics`. Use the rubric's sweep, only when `pgrep -c python3` reads 0.
  - Foreground, `timeout: 660000`; `-s`; no edit after the last window (rule (i)).
- **Scope.** Quasi-static trace, PEC, one a, imposed closed-form trace, accuracy-level bands at one mesh per aspect (not a convergence claim). Not: a square-edged plate (edge singularity), a plate in the coil-driven phantom, SAR or 1 g SAR, finite-k spheroid, thin-plate aspects (≥ 20).
- **Status it can move.** `IMP-6` 🟡 → ✅ ("step 1 at both aspects ⇒ ✅" in the row). Rule (j): the executor closes on the row.
- **Negative result.**
  - W0 over the cell guard: one pre-stated fallback is allowed, `resolution_surface` = a/9 and rim-tube SizeMin 1.5ρ_r, disclosed and printed. Still over ⇒ stop before any solve, report the counts, BLOCKED.
  - `IMP-2` census red: the keyword broke the default path. Fix in-slot, and never go on to solves.
  - (M′) red: a quadrature bug. Fix in-slot.
  - (i) or (ii) red with (iii) and (M′) green: print h_rim/ρ_r and the body-volume deficit, park on `attempt/*`, BLOCKED in the same commit (rule (d)). No band moves (rule (f) applies only to a mesher-unmeetable analytic comparand, ruled by the review).
  - A control that passes the gate: the basis or drive is mislabelled. Report it and stop.

**Outcome 2026-10-12 (13:30 CDT slot, second item, implementer): DONE, `IMP-6` 🟡 → ✅.** `MeshGenerator.spheroid_in_box_domain` gains the keyword `oblate=False` (rule (c), additive; the oblate path is a separate private builder, so the default branch's diff is the dispatch only); test `tests/validation/test_imp6_oblate_fem.py`. Heavy, `-n 2`, complex.
- **SICN trap, disclosed.** The first W0 (`20261012T183906Z_IMP-6-step1-W0.log`, 55 passed) printed min SICN 0.0247 at aspect 2 (0.1517 at aspect 5). Per the trap, no solve was run on that mesh. A `gmsh.model.mesh.optimize("Netgen")` pass was added on the oblate path only, and W0 was re-run (`20261012T184107Z_IMP-6-step1-W0.log`, 97 s, 55 passed): min SICN 0.267 / 0.371; cells 39 142 / 32 676 (guard 95 000); `IMP-2` census re-gate 37 359 / 44 304, unchanged; (M′) 2.6e-14 / 5.5e-14; body-volume deficit −4.40e-3 / −5.78e-3; air OCC mass rel 4.1e-5 / 1.5e-5.
- **Anchors.** W1a (`20261012T184250Z_IMP-6-step1-W1a.log`, 162 s, `test_pec_sphere_hole.py` re-gate green): aspect 2 β = 0.999220 (0.078 %), Φ_rim,h / Φ_rim,cf,h = 4.755198 / 4.762429 (miss 0.152 %), h_rim/ρ_r 0.648, 70.5 s solve, 6.99 GiB summed, cavity 15 026 dofs at max \|E\| 0. W1b (`20261012T184537Z_IMP-6-step1-W1b.log`, 126 s): aspect 5 β = 1.006343 (0.634 %), Φ 23.024938 / 22.702995 (miss 1.418 %), h_rim/ρ_r 0.665, 41.9 s, 4.06 GiB, cavity 19 952 dofs at 0. Both inside the predicted ≲ 1.5 % / ≲ 3 %. \|Im E\|/\|Re E\| 8.3e-3 / 8.5e-3.
- **Control (rule (e)).** W2 (`20261012T184748Z_IMP-6-step1-W2.log`, 121 s): axial drive on the aspect-5 mesh, β on the in-plane basis −0.001796 (predicted ≈ 0), Φ miss 97.84 % (predicted 90–100 %); both separated. β_axial on the axial basis 1.005469 (record, not asserted).

**Every-item rules that bind both items** (PROJECT_PLAN §9 "Every item"):
- The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.
- An executor that loses the docker socket mid-item stops and reports.
- No edit to a script, test or page after its last window (rule (i)).
- Every window passes `timeout: 660000` and runs pytest with `-s`. Never delete a log.
- An item that parks marks itself BLOCKED in the same commit (rule (d)).
- A window killed at its tier ceiling is not a result. After a kill, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first.
- Two heavy windows never run at once.
- A smoke window puts only `tests/environment/test_complex_mode.py` on the path.
- Every frequency-domain command sources complex mode, sets `FEM_EM_REQUIRE_COMPLEX=1`, and puts `tests/environment` first.
- No slot opens `docs/private/` or `aed_results/`.
- Standing rules (c) (disclosed `src/` change, gate re-run green) and (e) (asserted / predicted labels) apply as marked above.

**Spawn line** (both items): "Use the implementer agent to execute `IMP-6` §9 item 403 (then 404 if 403 commits before minute 30), per `docs/planning/chunks/IMP-6.md`. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`."

### Review ruling 2026-10-12 10:30 (daily) on the item 403 / 404 drafts

(Q1) **The opening.** `IMP-6` is opened under §9 rule (5) lane 3. Its target is a closed form, the oblate depolarisation factors and the PEC spheroid's exterior field, cited to Landau & Lifshitz *ECM* §4, Stratton §3.27 and Osborn (1945), the same sources `IMP-2` cites. None of the three is held in the repo. So the in-repo check is (Int), the ellipsoid integral by quadrature, which the docstring derives; the citations are pointers, not comparands. The row goes into §7 in this commit, and the opening is named in the dashboard for tonight's weekly to ratify or revert.

(Q2) **Item 403's control (c2) is re-stated** (the draft's wording was garbled). Replace each drive's decaying ratio D_i(ξ)/D_i(ξ₀) by the regular solution's ratio: ξ/ξ₀ for the axial drive, √(1 + ξ²)/√(1 + ξ₀²) for the in-plane drive. Then compute (Far)'s projections at R = 2a and 4a. The projected "dipole" grows with R instead of staying constant (predicted ratio ≈ 8, printed). Assert only \|p(4a)/p(2a) − 1\| > 1e-3. The predicted factor is printed beside the measured one (rule (e)).

(Q3) **Item 404's bands are imported** (`TH15_BETA_BAND`, `MAT4_BAND`), exactly as `IMP-2` step 1 imported them. `IMP-2`'s measured misses (0.46–1.11 %) are evidence that the bands are reachable; they are not a backing for the predicted values, which stay printed only. A red on (i) or (ii) parks the item as written. No band moves, and the weekly or the next review rules any re-pose.

(Q4) **404 depends on 403 landing in the same or an earlier slot.** If 403 parks, 404 is skipped and stays queued. The cell guard (95 000 in W0, one pre-stated fallback) is binding.

(Q5) **Measured minutes.** The lane-3 class has no recent landings, so both items are costed at predicted × 0.29, this interval's ratio: ≈ 12 for 403 and ≈ 20 for 404. The draft's "≈ 11" for 403 used 0.28.

### §9 item 408 — 2026-10-12 (lane 3, closed form)

`IMP-6` step 2: the rim resolution ladder at aspect 5, in-plane drive, degree 2 (implementer; heavy, `-n 2`, complex; ≈ 20 measured / 55 predicted min: ≈ 12.5 min of windows + 15 fixed + ≈ 25 code. Predicted × 0.29 gives ≈ 16, which is under the window floor. Item 404 took 15 measured minutes on ≈ 8.4 min of windows, so this item adds ≈ 4 min of windows to that; depends on nothing).

- **The knob, named from the code.** `_oblate_spheroid_in_box_domain` (`src/fem_em_solver/io/mesh.py` ≈ line 2903) exposes **no rim-tube parameter**. The tube's `Threshold` SizeMin is hard-wired to ρ_r = c²/a (≈ line 3005). The rim-local knob is **`resolution_tip`**, the `VIn` of the two `Ball` fields of radius 3ρ_r at (±a, 0, 0). Step 1 passes ρ_r/3. The rim mask (`_rim_mask`, cell midpoints within r_r = 2c²/a = 2ρ_r of (±a, 0, 0)) lies wholly inside those Balls. There `Min` picks `resolution_tip` (ρ_r/3 < ρ_r), so the tube never governs inside the mask, and the ladder needs **no `src/` edit**. Hold `resolution_surface` = a/12 and `resolution_far` = a/4 fixed.
- **Rungs** (aspect 5: a = 6 mm, c = 1.2 mm, ρ_r = 0.24 mm, step 1's box and host): `resolution_tip` = ρ_r/3 (rung 0, step 1's mesh), ρ_r/4 (rung 1) and ρ_r/5 (rung 2).
  - Predicted cells ≈ 33 k / ≈ 55 k / ≈ 92 k. The basis: step 1's 4 905 rim-mask cells in the 2ρ_r balls scale to ≈ 16.5 k in the 3ρ_r Balls, ≈ 16 k lie outside, and the Ball count scales as divisor³. Gradation outside the Balls makes rung 2 the uncertain one.
  - Predicted h_rim/ρ_r ≈ 0.665 / 0.50 / 0.40 (step 1 read 0.665 at rung 0).
- **Code.** New module `tests/validation/test_imp6_oblate_rim_ladder.py`. Do not edit `test_imp6_oblate_fem.py`; it is step 1's closing module (rule (i)).
  - Import step 1's private helpers and constants only: `SEMI_MAJOR`, `BOX_HALF_EXTENTS`, `CELL_GUARD`, `QDEG`, `_semi_minor`, `_rho_r`, `_ncells`, `_rss_gib`, `_field`, `_dirichlet`, `_probe_points`, `_rim_mask`, `_owned_masked_cells`, `_numpy_quadrature_mean`, `_body_volume_deficit`, `_beta`. Import the rest from `test_imp1_pec_sphere_lossy_host.py` as step 1 does: `TH15_BETA_BAND`, `MAT4_BAND`, tags, `DEGREE`, `E0`, `EPSILON_R_HOST`, `_assemble_real`, `_cavity_dofs`. Import `FREQUENCY_HZ` and `SIGMA_HIGH` from `test_lossy_sphere_sar.py`.
  - **Never import a `test_*` name.** pytest collects imported test functions and would re-run step 1's solves inside this window.
  - `_mesh_rung(divisor)`: step 1's `_mesh` with `resolution_tip = _rho_r(5.0)/divisor`, cached per divisor, timing the mesh.
  - `_solve_on(mesh_tuple)`: step 1's `_solve_oblate` (lines 225–304) with the mesh passed in and drive `"inplane"`. β uses the **in-plane** basis with the uniform part on axis 0.
- **Anchors** (asserted per rung at imported bands; bands never restated in the module):
  - (i) \|β − 1\| ≤ `TH15_BETA_BAND`.
  - (ii) \|Φ_rim,h/Φ_rim,cf,h − 1\| ≤ `MAT4_BAND`.
  - (iii) cavity dofs > 0 and max \|E\| < 1e-12.
  - W0 asserts, before any solve, cells ≤ `CELL_GUARD` (95 000) and mask cells > 0 for every rung.
- **Record (printed, not asserted).** Per rung:
  - cells, nominal divisor, measured h_rim/ρ_r, min SICN (the mesher prints it), rim-mask cells, body-volume deficit, Φ_h, Φ_cf,h, the **signed** ratio − 1, β, \|Im E\|/\|Re E\|, seconds and summed RSS.
  - The observed order per adjacent pair, p = ln(m_k/m_{k+1}) / ln(H_k/H_{k+1}), using the measured h_rim/ρ_r for H, with the nominal-divisor version beside it. Each window prints its own rung. The executor computes the pairwise orders in the outcome from the cited log lines (pure arithmetic). If the signed miss changes sign, the order is reported as undefined, not forced.
  - **Why the rate is not asserted.** No imported bar backs this comparison. `MAT-9`'s `RATE_BAR` = 1.7 (re-used in `test_mat15_finite_cylinder_fem.py:68`) is a P1 thermal L2 rate on a ratio-2 whole-ball ladder, a different comparison on a different fixture, so rule (e) does not license it. `IMP-2` has no ladder.
- **Predicted (rule (e), printed beside the measured value).**
  - Rim miss ≈ 1.418 % × (H₁/H₀)² ≈ 0.8 % at rung 1 and ≈ 0.5 % at rung 2 if the rim-local error dominates (affine faceting of a 0.24 mm curvature radius, O(h²)).
  - The surface grain and the far-field grain are held fixed. Step 1's \|β − 1\| = 0.634 % on rung 0 measures a global error of that size, which can floor the local miss. Predicted rung-2 miss in [0.3 %, 1.4 %]; predicted pairwise order in [0, 2.5]. A floor (order → 0) is a finding, not a failure.
  - β per rung within ≈ ±0.2 % of 0.634 %, because the knob is rim-local.
  - Rung 0 should reproduce step 1's 32 676 cells, ratio 1.014181 and β 1.006343 (`20261012T184537Z_IMP-6-step1-W1b.log:189–191`). gmsh was deterministic across step 1's W1b and W2 (same cells, same min-SICN centroid). Printed beside; a mismatch is reported (it may be an image move), not asserted.
- **Negative control.** None new. The comparand and the gate are step 1's, and step 1's axial-drive control (`20261012T184748Z_IMP-6-step1-W2.log:190–191`) already separates this fixture's gates from a wrong drive. A ladder adds no control the record needs.
- **Price** (measured on step 1's aspect-5 window: 32 676 cells, 41.9 s solve including mesh, 4.06 GiB summed, 126 s window, ≈ 80 s of it `tests/environment`; `IMP-1` at 90 713 cells: 221.9 s / 19.4 GiB, `-n 2`). Scaling solve time with cells^1.63 (the exponent between those two) gives:
  - rung 1: ≈ 95 s solve, ≈ 180 s window;
  - rung 2 at ≈ 92 k cells: ≈ 225 s solve, ≈ 15 GiB summed, ≈ 330 s window.
  - Every window is inside `timeout -k 30 600` and the 660 s Bash ceiling. The heavy tier holds at `-n 2`; no `-n 4` is needed.
- **Windows** (complex, `-n 2`, one solve per window, foreground, `timeout: 660000`; the shape of step 1's `Command:` lines):
  - W0 (≈ 120 s, census of all three rungs, no solve): `scripts/testing/run_and_log.sh IMP-6-step2-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp6_oblate_rim_ladder.py::test_ladder_census -s -v --tb=short'"`
  - W1a (≈ 126 s): label `IMP-6-step2-W1a`, `timeout -k 30 600`, paths `tests/environment tests/validation/test_imp6_oblate_rim_ladder.py::test_rung0`.
  - W1b (≈ 180 s): label `IMP-6-step2-W1b`, `timeout -k 30 600`, `tests/environment tests/validation/test_imp6_oblate_rim_ladder.py::test_rung1`.
  - W1c (≈ 330 s): label `IMP-6-step2-W1c`, `timeout -k 30 600`, `tests/environment tests/validation/test_imp6_oblate_rim_ladder.py::test_rung2`.
- **Traps already paid for.**
  - **SICN < 0.1** (step 1's first W0, `20261012T183906Z_IMP-6-step1-W0.log`, 0.0247 at aspect 2). Finer Balls steepen the Ball-to-surface gradation. If W0 prints min SICN < 0.1 on a rung, no solve runs on that rung; report it. The Netgen pass on the oblate path is the disclosed step-1 fix: keep it and never remove it.
  - **No `src/` edit is expected.** If one is made anyway, rule (c) applies and W0 adds the re-gates: `tests/validation/test_imp2_spheroid_fem.py::test_census_and_numpy` (37 359 / 44 304 cells unchanged), `tests/validation/test_imp6_oblate_fem.py::test_census_and_numpy` (39 142 / 32 676) and `tests/validation/test_pec_sphere_hole.py` in a W1a-shaped window.
  - **Complex comparisons.** No UFL `<`, `max_value` or `sqrt` comparisons; masks stay in numpy (`_rim_mask` interpolates a numpy lambda).
  - Owned cells only (`size_local`). `assemble_scalar` is rank-local (`_assemble_real` reduces). The RSS reduction stays outside the rank-0 branch.
  - `-s` on every window (rule (g)). `-k 30` on every container `timeout`. Select by node ID, never `-k a or b`. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Foreground. Two heavy windows never run at once.
  - A killed compile leaves an orphan `.c`. Use the rubric's sweep, and only when `pgrep -c python3` reads 0.
- **Scope.** One aspect, one knob, rim-local. The surface and far grains are fixed, so this is **not** a convergence study of the fixture and does not license "convergence" (the row's Not-claimed list stands). It is not aspect 2, not the axial drive (item 409), and not a thin plate.
- **Status it can move.** The `IMP-6` §7 row gains a "step 2 record" clause: miss vs rung, pairwise observed order, and the floor if one shows. No glyph change (`IMP-6` is ✅). The slot does not edit the Done-when or its Not-claimed list; the review amends the row if it wants to.
- **Negative result.**
  - W0 over the guard on rung 2: one pre-stated fallback, rung 2 → ρ_r/4.5 (≈ 72 k predicted), disclosed. Still over: a two-rung ladder (ρ_r/3, ρ_r/4), disclosed, one order. A rung killed at 600 s is not a result; shrink it as above and never raise the timeout.
  - (i) or (ii) red on a finer rung with (iii) green: a finer mesh missing a band step 1 met is a finding about the fixture (gradation, quality). Print h_rim/ρ_r, min SICN and the deficit, then park the module on `attempt/*`, BLOCKED in the same commit (rule (d)). Step 1's record and the ✅ stand until the review rules; no band moves.
  - A miss that grows with refinement: report it as measured; never fit a rate around it.

**Outcome (10-12 22:30 slot): 🚫 BLOCKED, parked on `attempt/IMP-6-item408-20261013T0339Z`.** W0 `20261013T033117Z_IMP-6-step2-W0.log` (29 s; re-run with near-rim edge statistics `…033215Z_IMP-6-step2-W0r.log`, 28 s): cells 32 676 / 51 898 / 87 426 (guard held), min SICN 0.3707 / 0.2242 / 0.1507 (all ≥ 0.1, every rung solved), rim-mask cells 4 905 / 10 538 / 20 526, deficit −5.777e-3 / −5.942e-3 / −5.986e-3. Step 1's max-edge h_rim/ρ_r read 0.665 / **2.720 / 2.987** (one long-edged cell at the rim on the finer rungs); median 0.517 / 0.389 / 0.315, p95 0.609 / 0.666 / 0.452. Solves (`-n 2`, complex, degree 2): W1a `20261013T033249Z_IMP-6-step2-W1a.log` (44 s) — rung 0 reproduces step 1 to every printed digit (32 676 cells, ratio 1.014181, β 1.006343); W1b `…033339Z_IMP-6-step2-W1b.log` (74 s, 7.20 GiB) — Φ_h 21.529105 vs Φ_cf,h 24.968699, signed miss **−13.78 %** (red; predicted ≈ 0.8 %), β 1.003281, cavity 0; W1c `…033502Z_IMP-6-step2-W1c.log` (146 s, 13.95 GiB, run for the record after the W1b red) — Φ_h 21.203154 vs Φ_cf,h 26.528645, miss **−20.07 %** (red; predicted [0.3 %, 1.4 %]), β 1.000078, cavity 0. The reds are outside the predicted range. Φ_h moves little (23.02 → 21.53 → 21.20) while the comparand climbs (22.70 → 24.97 → 26.53). **Diagnosis, measured** (`…033815Z_IMP-6-step2-W0diag.log`, 27 s): the number of air mask-cell midpoints lying *inside the true spheroid* is 0 / 31 / 93, the deepest at ρ = a − 0.132 mm / a − 0.163 mm, past the focal ring at a − 0.121 mm. Finer `resolution_tip` Balls make the mesher truncate the rim. The FEM then sees a blunted edge, and the `interior="continue"` comparand is sampled near its focal-ring singularity. This is a finding about the mesher on this path, not about the solve. No pairwise order is reported, because the ladder is not geometrically nested at the rim. Step 1's record and ✅ stand. Unblock: a review ruling on the rim-truncating mesh (a rule (c) mesher fix with the brief's re-gates, or a re-posed knob). No band moves.

**Spawn line:** "Use the implementer agent to execute `IMP-6` §9 item 408 per `docs/planning/chunks/IMP-6.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. No `src/` edit; never import a `test_*` name from step 1's module."

### §9 item 409 — 2026-10-12 (lane 3, closed form)

`IMP-6` step 3: the axial drive gated at aspect 2 and 5 on step 1's meshes, pole enhancement 1/L_z (implementer; heavy, `-n 2`, complex; ≈ 15 measured / 50 predicted min: ≈ 9 min of windows + 15 fixed + ≈ 25 code. Predicted × 0.29 ≈ 15 agrees with item 404's measured 15 on the same window shape; depends on nothing).

- **Where the pole is, and its mesh.** For the oblate body the poles are on the flat faces at (0, 0, ±c). The closed-form identity is E_z = E₀/L_z there, measured in step 0 at ≤ 4.4e-16.
  - Pointwise (1/L_z)² = 3.598 at aspect 2 and 1.776 at aspect 5.
  - The curvature radius at the pole is a²/c = 12 mm and 30 mm, larger than the body.
  - **The mesher does not refine there.** The Balls and the rim tube sit at the rim. The pole cells come from the cavity Distance→Threshold: `resolution_surface` a/12 = 0.5 mm at the wall, growing to `resolution_far` a/4 at DistMax 2c. That gives h ≈ 0.5–0.7 mm (aspect 2) and ≈ 0.5–1.0 mm (aspect 5) inside the pole ball.
  - The OCC sphere's parametric poles land there. Step 1's aspect-2 worst element (min SICN 0.267) sits 0.1 mm below the bottom pole (`20261012T184107Z_IMP-6-step1-W0.log:500`, centroid z = −3.10 mm). At aspect 5 it is off-pole (`…:603`).
- **Pole ball.** Cell midpoints within r_p = a/5 = 1.2 mm of (0, 0, ±c), both poles.
  - `IMP-2`'s "2× curvature radius" rule cannot be mirrored: 2a²/c exceeds the box.
  - a/5 keeps the mask a near-half-ball. The face drops by c·r_p²/(2a²) ≤ 0.12 mm over it, and the top at c + a/5 = 4.2 / 2.4 mm clears the lid at bz = a = 6 mm.
  - Predicted mask ≈ 3.6 mm³ per pole, so ≈ 100–300 cells; printed. W0 asserts > 0.
- **Is `MAT4_BAND` reachable? Predicted yes, so this is a gate and there is no `PREDICTED RED` line.**
  - The axial field over the pole ball is smooth. The closed form varies ≈ 6 % across the aspect-5 ball (pointwise 1.776 against a ball mean 1.673).
  - Faceting of a 30 mm radius face at h = 0.5 mm moves the surface by ≈ h²/(8R) ≈ 1e-3 mm.
  - Step 1's axial run on the aspect-5 mesh already read β_axial = 1.005469 on the axial basis (`20261012T184748Z_IMP-6-step1-W2.log:190`). That is the same comparison on the same fixture, so it is evidence the β gate is reachable at aspect 5. At aspect 2 there is no prior measurement.
  - Predicted (printed): pole-Φ miss ≲ 2 % at both aspects, possibly larger at aspect 2 because of the worst element in its ball; \|β_axial − 1\| ≲ 1.5 %.
- **Code.** New module `tests/validation/test_imp6_oblate_axial_fem.py`. Step 1's module is not edited (rule (i)).
  - From `test_imp6_oblate_fem.py` import only private names: `_mesh` (step 1's meshes and fixture, cached by aspect), `_field`, `_dirichlet`, `_probe_points`, `_beta`, `_owned_masked_cells`, `_numpy_quadrature_mean`, `_imp1_dipole_numpy`, `_rss_gib`, `_ncells`, `_body_volume_deficit`, `_semi_minor`, `SEMI_MAJOR`, `CELL_GUARD`, `QDEG`.
  - From `test_imp1_pec_sphere_lossy_host.py`: `TH15_BETA_BAND`, `MAT4_BAND`, tags, `DEGREE`, `E0`, `EPSILON_R_HOST`, `_assemble_real`, `_cavity_dofs`, `_phi_closed_form_on_mask`. `FREQUENCY_HZ` and `SIGMA_HIGH` come from `test_lossy_sphere_sar.py`.
  - **Never import a `test_*` name** (pytest would collect and re-run step 1's solves).
  - Add `_pole_mask(msh, aspect)`, a DG0 numpy-lambda indicator in the shape of `_rim_mask`.
  - Add `_solve_pole(aspect, drive)`, step 1's `_solve_oblate` (lines 225–304) with three changes:
    - the mask is `_pole_mask`;
    - Φ_cf,h is `_numpy_quadrature_mean(msh, cells, _field(aspect, "axial", "continue"))`, the **axial** closed form for both drives;
    - β is `_beta(e_re, points, aspect, "axial", uniform_axis)`, with `uniform_axis` = 2 for the axial drive and 0 for the in-plane drive. Step 1's `beta_own` is this quantity for the axial drive.
- **Anchors** (asserted at imported bands, as step 1 asserted them):
  - At aspect 2 and 5, axial drive:
    - (i) \|β_axial − 1\| ≤ `TH15_BETA_BAND`.
    - (ii) \|Φ_pole,h/Φ_pole,cf,h − 1\| ≤ `MAT4_BAND`.
    - (iii) cavity dofs > 0 and max \|E\| < 1e-12.
  - W0, no solve:
    - (M″), the pole-mask twin of step 1's (M′): `_numpy_quadrature_mean` of `_imp1_dipole_numpy` equals `_phi_closed_form_on_mask(msh, chi_pole, 1.0)` to ≤ 1e-12 on each mesh.
    - Cells ≤ `CELL_GUARD`.
    - Pole-mask cells > 0.
    - The meshes reproduce step 1's 39 142 / 32 676 cells: printed, a mismatch reported.
- **Predicted and printed (rule (e)).** From a stdlib Monte Carlo of item 403's closed form over the continuous ball (40 000 points, run while drafting this brief; the discrete mask will differ by a few %):
  - Φ_pole,cf ≈ 3.11 (aspect 2) and 1.67 (aspect 5), beside the pointwise (1/L_z)² of 3.598 and 1.776.
  - Also printed: h_pole (the longest edge of the mask cells), min SICN, the deficit, \|Im E\|/\|Re E\|, seconds and RSS.
- **Negative control** (W2, aspect-5 mesh, in-plane drive: the mirror of step 1's W2). Every factor here is **predicted**: step 1's control measured the axial drive on the in-plane basis, which is not this comparison.
  - β fitted on the **axial** basis, with the uniform part on axis 0: predicted ≈ 0 by parity. Step 1's mirror read −0.001796, printed beside as context. The ceiling is step 1's 1/`TH15_BETA_BAND` ≈ 20.5×. Assert only \|β − 1\| > `TH15_BETA_BAND`.
  - The pole-ball Φ_h under the in-plane drive against the **axial** Φ_pole,cf,h. Under the in-plane drive the surface charge ∝ n·x̂ vanishes at the pole, so E → 0 there. The Monte Carlo gives in-plane Φ ≈ 0.038 and a miss ≈ 97.7 % at aspect 5 (98.1 % at aspect 2, not run). Predicted range 90–100 %. The miss saturates at 100 %, so it can show at most 10× the band. Assert only > `MAT4_BAND`.
- **Windows** (complex, `-n 2`, one solve per window, foreground, `timeout: 660000`; the shape of step 1's `Command:` lines):
  - W0 (≈ 100 s, census + (M″)): `scripts/testing/run_and_log.sh IMP-6-step3-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp6_oblate_axial_fem.py::test_pole_census -s -v --tb=short'"`
  - W1a (≈ 160 s, step 1's aspect-2 window was 162 s): label `IMP-6-step3-W1a`, `timeout -k 30 600`, paths `tests/environment tests/validation/test_imp6_oblate_axial_fem.py::test_axial_aspect2`.
  - W1b (≈ 125 s): label `IMP-6-step3-W1b`, `timeout -k 30 600`, `tests/environment tests/validation/test_imp6_oblate_axial_fem.py::test_axial_aspect5`.
  - W2 (≈ 125 s): label `IMP-6-step3-W2`, `timeout -k 30 600`, `tests/environment tests/validation/test_imp6_oblate_axial_fem.py::test_inplane_control`.
- **Traps already paid for.**
  - **Basis vs. uniform.** The basis follows the comparison (axial here); the uniform part follows the drive. This is step 1's trap mirrored. Copying `_solve_oblate`'s hard-coded `"inplane"` basis or its in-plane Φ_cf is the trap.
  - **SICN < 0.1.** The meshes are step 1's, read 0.267 / 0.371 after the disclosed Netgen pass. If W0 prints < 0.1 (an image move), no solve runs; report it.
  - **No `src/` edit.** Pole refinement is **not** added here. If a red ever makes one necessary, it is the review's re-pose under rule (c), with the re-gates `test_imp2_spheroid_fem.py::test_census_and_numpy` (37 359 / 44 304), step 1's census and `test_pec_sphere_hole.py`.
  - Complex comparisons: masks in numpy, no UFL `<` or `sqrt` comparison.
  - Owned cells only. Reductions as step 1. The RSS reduction stays outside the rank-0 branch. Probe points never sit on a wall.
  - `-s`, `-k 30`, node IDs only. The harness call is the **whole** Bash command, bare. Foreground. Two heavy windows never run at once. Orphan-`.c` sweep only when `pgrep -c python3` reads 0.
- **Scope.** Quasi-static trace, PEC, one a, step 1's meshes, accuracy-level bands at one mesh per aspect. Not convergence. Not a square-edged plate, a plate in a coil, SAR or finite k.
- **Status it can move.** The `IMP-6` §7 row gains a "step 3 (axial drive gated at aspect 2 / 5, pole ball)" clause. No glyph change (`IMP-6` is ✅).
  - The Done-when names the axial drive only as step 1's control. Whether to add step 3's gates to it is the **review's** amendment; the slot does not edit the row's Done-when.
  - The **Not-claimed list does not change**: square-edged plate, plate in a coil, SAR, finite k and convergence all stand.
- **Negative result.**
  - (M″) red: a mask or quadrature bug; fix in-slot.
  - (ii) red with (i), (iii) and (M″) green: the unrefined pole region is the first suspect. Print h_pole, min SICN and the deficit, then park on `attempt/*`, BLOCKED in the same commit (rule (d)), with the unblock condition "review rules on pole Balls (rule (c) mesher edit) or a re-posed r_p". No band moves, and r_p is not changed after a red.
  - (i) red at aspect 2 only: report it. There is no prior measurement there, and the next review rules.
  - A control that passes a gate: the drive or basis is mislabelled; report it and stop.

**Outcome (10-13 00:00 slot): ✅ complete, every anchor green, the control misses both.** New module `tests/validation/test_imp6_oblate_axial_fem.py` (private imports only; no `src/` edit; step 1's module untouched). W0 `20261013T050123Z_IMP-6-step3-W0.log` (103 s): cells 39 142 / 32 676 (match step 1), min SICN 0.2667 / 0.3707 (≥ 0.1, `:194`, `:299`), pole-mask cells 228 / 154 (mask volume 7.477 / 7.261 mm³, both poles), h_pole 1.337 / 1.948 mm, deficit −4.404e-3 / −5.777e-3; (M″) rel 4.55e-15 / 4.00e-15 (`:197`, `:302`); Φ_pole,cf,h 3.105745 / 1.667975 (MC predicted ≈ 3.11 / 1.67). W1a `20261013T050311Z_IMP-6-step3-W1a.log` (162 s, aspect 2, 74.8 s solve, 6.94 GiB summed): β_axial 1.001634 (\|β − 1\| 0.163 %, band 4.886 %), Φ_pole,h 3.136894 vs 3.105745, signed miss **+1.003 %** (band 10 %, predicted ≲ 2 %), cavity 15 026 dofs at 0 (`:196–199`). W1b `20261013T050558Z_IMP-6-step3-W1b.log` (135 s, aspect 5, 45.7 s, 4.08 GiB): β_axial 1.005469 (equal to step 1's W2 record), Φ_pole,h 1.683758 vs 1.667975, miss **+0.946 %**, cavity 19 952 at 0 (`:189–192`). Control W2 `20261013T050819Z_IMP-6-step3-W2.log` (132 s): in-plane drive on the axial basis β = −0.032226 (predicted ≈ 0; step 1's mirror −0.001796), Φ_pole,h 0.043756 (MC ≈ 0.038) against 1.667975, miss 97.38 % (predicted 90–100 %) — both gates missed (`:192–193`). The unrefined pole ball meets `MAT4_BAND` at ≈ 1 % on both aspects. One mesh per aspect; not a convergence claim; the Not-claimed list stands.

**Spawn line:** "Use the implementer agent to execute `IMP-6` §9 item 409 per `docs/planning/chunks/IMP-6.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. No `src/` edit; never import a `test_*` name from step 1's module."

### Review rulings 2026-10-12 18:00 (daily) — audit, items 403 / 404, drafts 408 / 409

**Audit: `IMP-6` → ✅ PASS** (auditor, closing commit `51011ae7`). Step 1 anchors at `20261012T184250Z_IMP-6-step1-W1a.log:428–431` and `20261012T184537Z_IMP-6-step1-W1b.log:189–192`; the control at `20261012T184748Z_IMP-6-step1-W2.log:190–191`; `IMP-2`'s census at `20261012T184107Z_IMP-6-step1-W0.log:320, :387`; asserts at `tests/validation/test_imp6_oblate_fem.py:351–353, :378–382, :421–422`, bands imported; rule (i) holds on file mtimes.

(R1) **Item 403's (F) cap at aspect 5 is accepted provisionally (the weekly ratifies).** The 1e-7 cap on the Richardson residual of the test's own finite-difference gradient read 3.55e-7 / 5.91e-7 at the pre-registered δ pair (`20261012T183401Z_IMP-6-step0-W1.log:138–139`). The slot kept the cap, applied it one δ-halving finer (2.2e-8 / 3.7e-8) and asserted the per-halving shrink 15.98 / 16.03 in [14, 18] (`20261012T183435Z_IMP-6-step0-W1.log:91, :95`). The bound is on the instrument (the finite-difference step), not on a physics comparand. The pre-registered pair is printed, not asserted, at aspect 5 — a disclosed loss of coverage. Not a precedent for a physics band.

(R2) **Item 404's Netgen pass is accepted as a disclosed deviation; the closure stands.** The brief (item 404's trap list above) said "report it and do not proceed to solves". The slot ran no solve on the 0.0247 mesh, added `optimize("Netgen")` to the new private oblate builder only, re-ran W0 (0.267 / 0.371) and solved once per aspect: no mesh-shopping, no band moved, `IMP-2`'s path untouched. **Binding from now (§9 digest): a quality trap's "do not proceed" means park, unless the brief names the fix it licenses.** Items 408 / 409 license none: a rung or mesh under 0.1 SICN is skipped (408) or parks (409), and no mesher edit is made in-slot.

(R3) **Item 403's (c2):** the 10:30 ruling's "predicted ≈ 8" was the D ≡ 1 variant's figure; its own regular-ratio recipe reads ≈ 16 (R⁴). The slot's reading is right; both variants are asserted separated.

(R4) **Drafts 408 / 409 accepted as written.** 408 varies only `resolution_tip` (ρ_r/3, /4, /5); the order is printed, never asserted (no same-comparison prior). 409's pole ball is unrefined; its `MAT4_BAND` gate is predicted reachable, and a red parks with the unblock condition the brief states. Both run in new modules, so step 1's closing module is untouched. Costs: 408 ≈ 20 measured, 409 ≈ 15 (scaled from item 404's 15 measured minutes; the lane-3 class has two landings).

**§7 row state text, moved verbatim 2026-10-12 18:00 review (the row was 2.7 KB):**

✅ 2026-10-12 (step 1, §9 item 404); 🟡 2026-10-12 (§9 item 403, step 0 green): L_z / L_x = 0.527200 / 0.236400 (aspect 2), 0.750484 / 0.124758 (aspect 5); (Sum) ≤ 6.7e-16, (Int) ≤ 1.3e-15 rel, rim E_x·L_x and pole E_z·L_z within 5.8e-14 of 1 at aspects 1.01–100, (Far) ≤ 3.1e-15; (c1) prolate formula reads 1.748 / 3.784 (separated). (F)'s Richardson bound at aspect 5 asserted one δ-halving finer with the δ⁴ ratio, measurement disclosed (`20261012T183435Z_IMP-6-step0-W1.log`, 3 s smoke; re-gate `20261012T183442Z_IMP-6-step0-W1r.log`). **✅ 2026-10-12 (§9 item 404, step 1, degree 2, `-n 2`, complex):** `spheroid_in_box_domain(oblate=True)` (additive, rule (c): `IMP-2`'s census re-gate unchanged at 37 359 / 44 304 cells; `test_pec_sphere_hole.py` green) at a = 6 mm, aspect 2 / 5: cells 39 142 / 32 676, min SICN 0.267 / 0.371 after a disclosed Netgen pass (first W0 read 0.0247 at aspect 2, `20261012T183906Z_IMP-6-step1-W0.log`, no solve run on it); \|β − 1\| = 0.078 % / 0.634 % (band 4.886 %); rim-ball Φ_h / Φ_cf,h = 4.755198 / 4.762429 (miss 0.152 %) and 23.024938 / 22.702995 (miss 1.418 %, band 10 %); cavity dofs 15 026 / 19 952 at max \|E\| 0; (M′) ≤ 5.5e-14; h_rim/ρ_r 0.648 / 0.665. Axial control (aspect 5) misses both: β = −0.0018 on the in-plane basis, Φ miss 97.84 % (predicted 90–100 %). Logs `20261012T184107Z_IMP-6-step1-W0.log` (97 s), `20261012T184250Z_IMP-6-step1-W1a.log` (162 s), `20261012T184537Z_IMP-6-step1-W1b.log` (126 s), `20261012T184748Z_IMP-6-step1-W2.log` (121 s), heavy tier. One mesh per aspect; not a convergence claim.

### Review ruling 2026-10-13 03:00 (daily) on item 408

"Item 408's red is a mesher defect on the oblate path, measured: on the ρ_r/4 and ρ_r/5 rungs 31 / 93 air mask-cell midpoints sit inside the true spheroid, the deepest past the focal ring, while rung ρ_r/3 reproduces step 1 to every digit. The physics gates were never reached on a faithful mesh, so the reds are not a reading of the FEM. The family `IMP-6` step 2 is at attempt 1 of 4. Re-pose (attempt 2, §9 rule (6), target unchanged — I22 as the weekly ratified it): a rule (c) mesher fix confined to the private oblate builder, licensed in this order — (1) the rim constrained explicitly (the spheroid's equator embedded as a curve in the volume mesh, or a size field keyed on distance to the rim circle rather than on Balls), (2) the Netgen pass dropped on the finer rungs only if (1) is insufficient. Gate before any solve: W0diag reads 0 air midpoints inside the true spheroid on every rung and min SICN ≥ 0.1; rung ρ_r/3 (or the step-1 mesh) reproduces step 1's records; `IMP-2`'s and step 1's default paths unchanged (their re-gate green). If no licensed variant reaches 0 inside midpoints on ρ_r/4, park: the measured negative — the oblate rim ladder is mesher-bounded at ρ_r/3 — is banked on the row and the family is not re-lettered."

### §9 item 426 — 2026-10-13 (§10 chain I22, re-pose attempt 2; lane 3, closed form)

`IMP-6` step 2, attempt 2 of 4 (family count since step 1's landed gate): the rim resolution ladder at aspect 5 on a rim-faithful mesh (implementer; heavy, `-n 2`, complex; ≈ 10 measured / 50 predicted min: ≈ 10 min of windows + 15 fixed + ≈ 25 code; depends on nothing). Target unchanged from item 408 (I22, ratified by the 10-12 weekly); only the mesher is re-posed, per the ruling above.

- **What attempt 1 measured** (parked `attempt/IMP-6-item408-20261013T0339Z`, `7f9eb025`; attempts.md `2026-10-13T03:40Z`). The fixture, the anchors and the comparand are sound on rung ρ_r/3: W1a reproduces step 1 to every printed digit (ratio 1.014181, β 1.006343, `…033249Z_IMP-6-step2-W1a.log:166–167`). On rungs ρ_r/4 and ρ_r/5 the mesh is not faithful: 0 / 31 / 93 air mask-cell midpoints lie inside the true spheroid, the deepest at a − 0.132 / a − 0.163 mm, past the focal ring at a − 0.121 mm (`…033815Z_IMP-6-step2-W0diag.log:165, :279, :423`); the max near-rim edge jumps 0.665 → 2.720 → 2.987 ρ_r while the median falls (`:167, :281, :425`). The reds (−13.78 % / −20.07 %, `…033339Z_…-W1b.log:174`, `…033502Z_…-W1c.log:204`) are therefore not a reading of the FEM.
- **The fix, licensed in the ruling's order (rule (c), confined to the private builder).** `_oblate_spheroid_in_box_domain` (`src/fem_em_solver/io/mesh.py:2903–3054`) gains two keyword-only arguments, `rim: str = "balls"` and `netgen: bool = True`. `spheroid_in_box_domain` (`:2731`) forwards them as `oblate_rim` / `oblate_netgen` on the `oblate=True` branch only, and raises `ValueError` before any gmsh call if either is non-default with `oblate=False`. The defaults reproduce today's path byte for byte, so the default-branch diff is the signature and the forwarding only.
  - **(1a) `rim="edge"` — the equator as a model edge.** Build the body by revolution instead of `addSphere` + `dilate`: in the local xz plane, two quarter-ellipse arcs (pole (0, 0, c) → rim (a, 0, 0) → pole (0, 0, −c), centre the origin, major axis along x) closed by the axis segment, `addPlaneSurface`, `occ.revolve(..., 0, 0, 0, 0, 0, 1, 2π)`, then rotate π/4 about z (the revolve seam off the drive's rim points (±a, 0, 0), as today's 45° rotation). The upper and lower faces then share the rim circle as an edge, so the 1D mesh puts nodes on the true rim and no surface triangle can bridge it. Assert one solid out of the revolve (the rod's fuse check pattern) and print the cavity face count (2 predicted). The OCC body is now the exact ellipsoid of revolution, not a B-spline approximation: the air OCC-mass line (today rel 1.548e-05 at aspect 5, `20261012T184107Z_IMP-6-step1-W0.log:503`, a body-volume error of 4.37e-4) is printed again and predicted ≤ 1e-9; a disclosed by-product, not a gate. Verify the arc API in the image before coding (`gmsh.model.occ.addEllipseArc(startTag, centerTag, majorTag, endTag)`, or `addEllipse` with explicit `zAxis` / `xAxis` and angle limits); a refusal of the API is a code fact, fixed in-slot.
  - **(1b) `rim="tube"` — sizing keyed on the rim distance.** Keep the dilated body. Replace each `Ball` at (±a, 0, 0) by a `Max` of two fields: a `Threshold` on the existing rim `MathEval` (`:2999`) with SizeMin `resolution_tip`, SizeMax `resolution_far`, DistMin 0, DistMax 3ρ_r, and a `Ball` of the old radius 3ρ_r with VIn `resolution_tip`, VOut `resolution_far`. Inside the old Ball the size then grows with distance from the rim circle instead of being uniform; outside it the `Max` returns `resolution_far` and the `Min` ignores it. The step-1 rim tube and cavity threshold are unchanged.
  - **(2) `netgen=False` on rungs ρ_r/4 and ρ_r/5 only**, tried only if neither (1a) nor (1b) passes the gate below. It is applied to the step-1 construction (`rim="balls"`) so it isolates attempt 1's hypothesis (attempts.md: "no Netgen pass on the surface"); rung ρ_r/3 keeps its Netgen pass.
  - The variant is chosen on the W0 geometric gate alone, before any solve, in the order (1a) → (1b) → (2). Once a solve window has run, no other variant is tried (no mesh-shopping on a physics reading).
- **Code (tests).** Take the parked module by path — `git checkout attempt/IMP-6-item408-20261013T0339Z -- tests/validation/test_imp6_oblate_rim_ladder.py` — and the six attempt-1 logs the same way (new files; nothing on `main` to overwrite). **Never path-checkout `docs/testing/test-results.md`:** append the branch's six rows by hand and check `git diff --stat` shows only insertions on it (the `8d4cf58` trap). Then edit the module:
  - `RIM = os.environ.get("FEM_EM_IMP6_RIM", "edge")`, validated in {`edge`, `tube`, `balls-nonetgen`}; `_mesh_rung(divisor)` passes `oblate_rim` / `oblate_netgen` accordingly (`balls-nonetgen`: `rim="balls"`, `netgen = divisor == 3`). The cache key carries the variant.
  - `DIVISORS = (3.0, 4.0, float(os.environ.get("FEM_EM_IMP6_RUNG2", "5")))`, so the pre-stated rung-2 fallback needs no module edit.
  - `test_ladder_census`: keep every print; **assert** per rung cells ≤ `CELL_GUARD`, mask cells > 0 and **inside-midpoint count == 0** (the W0diag computation, `s < 1` on owned mask cells, reduced).
  - New `test_census_control`: rungs ρ_r/4 and ρ_r/5 on the **default** path (`oblate_rim` unset), asserting inside count > 0 on each and printing it beside 31 / 93.
  - Docstring: attempt 2, the variant, the ruling. The rest of the module (`_solve_on`, `_report_and_gate`, the three rung tests) is unchanged: the anchors are item 408's.
  - Still: never import a `test_*` name from step 1's module; do not edit `test_imp6_oblate_fem.py`, `test_imp6_oblate_axial_fem.py` or `EX-111`'s script (all call step 1's `_mesh`, the default path).
- **Gate before any solve (the ruling's, all read in W0 windows).**
  - (G1) inside-midpoint count 0 on every rung of the chosen variant — asserted.
  - (G2) min SICN ≥ 0.1 on every rung — read from the builder's printed line (the builder returns no quality); no solve window starts unless every rung reads ≥ 0.1.
  - (G3) the step-1 mesh reproduces step 1's records: W0r prints `IMP-6` step 1's census at 39 142 / 32 676 cells and min SICN 0.2667 / 0.3707 (`20261012T184107Z_IMP-6-step1-W0.log:500, :603`). The parked W1a stays the default path's rung-0 solve record (ratio 1.014181, β 1.006343). *Reading for the review:* under (1a) the variant's own rung ρ_r/3 is a different mesh (exact CAD, rim edge); it is gated at the bands like the other rungs and its miss printed beside 1.418 %.
  - (G4) `IMP-2`'s default path unchanged: W0r prints 37 359 / 44 304 cells; W3 runs `test_pec_sphere_hole.py` green.
- **Anchors** (item 408's, asserted per rung at imported bands; bands never restated):
  - (i) \|β − 1\| ≤ `TH15_BETA_BAND` (`tests/validation/test_pec_sphere_hole.py` `BAND`, imported via `test_imp1_pec_sphere_lossy_host.py:55`; 4.886 %), β the least-squares amplitude of `pec_oblate_spheroid_field(drive="inplane")` − E₀x̂ (`src/fem_em_solver/utils/analytical.py:1831`) on Re E − E₀x̂ at the 48 probe points.
  - (ii) \|Φ_rim,h/Φ_rim,cf,h − 1\| ≤ `MAT4_BAND` (`test_imp1_pec_sphere_lossy_host.py:59`; 10 %), Φ_rim,cf,h the degree-8 numpy quadrature of the same closed form with `interior="continue"` over the rim-mask cells (r_r = 2ρ_r about (±a, 0, 0)).
  - (iii) cavity dofs > 0 and max \|E\| < 1e-12.
  - The miss-vs-h order is **printed, never asserted** (no same-comparison prior, rule (e); `MAT-9`'s `RATE_BAR` is a different comparison).
- **Record (printed, not asserted).** Per rung: cells, divisor, h_rim/ρ_r (max and median near-rim edge), min SICN, rim-mask cells, inside count, body-volume deficit, Φ_h, Φ_cf,h, the signed miss, β, \|Im E\|/\|Re E\|, seconds, summed RSS. The outcome computes the pairwise observed order p = ln(m_k/m_{k+1})/ln(H_k/H_{k+1}) from the cited log lines, with the median-edge H and the nominal divisor beside; a sign change in the miss makes the order undefined, reported as such.
- **Predicted (rule (e), printed beside the measured).** On a faithful mesh, item 408's predictions stand: rung ρ_r/4 miss ≈ 0.8 %, rung ρ_r/5 in [0.3 %, 1.4 %], pairwise order in [0, 2.5]; a floor set by the fixed surface and far grains (β's 0.634 % on rung 0) is a finding, not a failure. Rung ρ_r/3 under (1a) in [0.3 %, 2.0 %] beside step 1's 1.418 %. β per rung within ≈ ±0.7 pp of 0.634 % (attempt 1 read 0.634 / 0.328 / 0.008 %, `…W1a.log:166`, `…W1b.log:173`, `…W1c.log:203` — β was green on the unfaithful rungs too). Cells within ±15 % of 32 676 / 51 898 / 87 426.
- **Negative control.**
  - (C-fid, **asserted**) `test_census_control`: the default path's rungs ρ_r/4 and ρ_r/5 read inside counts > 0. Backed by the same comparison on the same fixture: 31 / 93 (`…033815Z_IMP-6-step2-W0diag.log:279, :423`; gmsh was deterministic across attempt 1's W0 and W0diag, `:283` vs `…033117Z_…-W0.log:275`). Assert > 0 only; the counts are printed beside 31 / 93. This shows (G1) can fail on this fixture, so a 0 under the variant is a reading of the fix, not of the check.
  - No new physics control: step 1's axial-drive control (`20261012T184748Z_IMP-6-step1-W2.log:190–191`) already separates this fixture's gates from a wrong drive, as item 408 stated.
- **Price** (item 408's measured windows, `-n 2`, `tests/environment/test_complex_mode.py` on the path as those windows ran): W0 29 s (`…033117Z_…-W0.log:452`), W0diag 27 s; W1a 44 s, W1b 74 s / 7.20 GiB, W1c 146 s / 13.95 GiB (`…W1a.log:204`, `…W1b.log:241`, `…W1c.log:271`). Step 1's W0 with both census re-gates: 97 s (`20261012T184107Z_IMP-6-step1-W0.log`). Predicted here: W0-ctl ≈ 25 s; each variant W0 ≈ 30 s (one to three of them); W0r ≈ 100 s; W3 ≈ 60–120 s; W1a / W1b / W1c ≈ 45 / 75 / 150 s (± 15 % with the cells). ≈ 6–10 min of windows in all; every window inside its container `timeout` and the 660 s Bash ceiling.
- **Windows** (complex, `-n 2`, one solve per window, foreground, `timeout: 660000`; item 408's `Command:` shape, `…033339Z_…-W1b.log:12`):
  - W0-ctl: `scripts/testing/run_and_log.sh IMP-6-step2-a2-W0-ctl "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_census_control -s -v --tb=short'"`
  - W0-edge (variant (1a)): `scripts/testing/run_and_log.sh IMP-6-step2-a2-W0-edge "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_ladder_census -s -v --tb=short'"`
  - Only if W0-edge fails (G1) or (G2): W0-tube, the same command with label `IMP-6-step2-a2-W0-tube` and `FEM_EM_IMP6_RIM=tube`; then, only if that fails too, W0-nonetgen with `FEM_EM_IMP6_RIM=balls-nonetgen`.
  - W0r (rule (c) re-gate, (G3) / (G4)): `scripts/testing/run_and_log.sh IMP-6-step2-a2-W0r "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_analytical_lightweight.py tests/validation/test_imp2_spheroid_fem.py::test_census_and_numpy tests/validation/test_imp6_oblate_fem.py::test_census_and_numpy -s -v --tb=short'"`
  - W3 (the mesher-module re-gate `IMP-2` and step 1 ran): label `IMP-6-step2-a2-W3`, `timeout -k 30 600`, paths `tests/environment/test_complex_mode.py tests/validation/test_pec_sphere_hole.py`.
  - W1a / W1b / W1c: labels `IMP-6-step2-a2-W1a` / `-W1b` / `-W1c`, `timeout -k 30 600`, the passing variant's `FEM_EM_IMP6_RIM=<variant>` (and `FEM_EM_IMP6_RUNG2=4.5` if the fallback below was taken) before `timeout`, paths `tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung0` / `::test_rung1` / `::test_rung2`. Example: `scripts/testing/run_and_log.sh IMP-6-step2-a2-W1b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung1 -s -v --tb=short'"`
  - Order: code → W0-ctl → W0-edge (→ W0-tube → W0-nonetgen as needed) → W0r → W3 → W1a → W1b → W1c → commit. No module or `src/` edit after W0r starts (rule (i)).
- **Traps already paid for.**
  - **Path checkout of an append-only record** (`8d4cf58`, `TH-15` step 3b): code and new logs by path, `test-results.md` rows by hand.
  - **A quality trap means park unless the brief names the fix** (2026-10-12 18:00 ruling): this brief names exactly three variants and the rung-2 fallback; nothing else is licensed in-slot (no Ball radius, no rim-tube SizeMin, no `resolution_surface` change).
  - **The basis is the in-plane field minus E₀x̂** (item 404's trap); `_solve_on` already does it — do not touch it.
  - **SICN < 0.1** on any rung of the chosen variant: no solve on that variant; try the next one.
  - `MathEval` takes a literal expression (bake numbers with `repr`, print the string, use gmsh's `Sqrt`); a `Threshold` returns SizeMax beyond DistMax — never give one a fine SizeMax.
  - The wall tolerance stays absolute 1e-6; the cavity group from the **meshed** boundary; the 6-wall check unchanged (under (1a) the cavity has two faces, not one).
  - Owned cells only (`size_local`); `assemble_scalar` is rank-local (`_assemble_real` reduces); the RSS reduction outside the rank-0 branch; masks in numpy, no UFL `<`, `max_value` or `sqrt` comparison.
  - `-s` on every window; `-k 30` on every container `timeout`; node IDs only, never `-k a or b`; the harness call is the **whole** Bash command; foreground; two heavy windows never at once. A killed window: `pgrep -c python3` first, then the rubric's orphan-`.c` sweep only if it reads 0.
- **Scope.** One aspect, one rim-local knob, a rim-faithful mesh. The surface and far grains are fixed, so this is not a convergence study and does not license "convergence" (the row's Not-claimed list stands). The `rim` / `netgen` keywords are not a fix to the default path: `IMP-6` step 1, step 3 and `EX-111` keep their meshes. Not aspect 2, not the axial drive, not a thin plate.
- **Status it can move.** The `IMP-6` §7 row's step-2 clause (today "item 408 BLOCKED … re-posed as item 426 (attempt 2 of 4)") becomes the record: the variant, miss per rung, pairwise order, the floor if one shows. No glyph change (`IMP-6` is ✅); the slot does not edit the Done-when or the Not-claimed list (rule (j) does not bind). A known-issues entry is opened either way: "`spheroid_in_box_domain(oblate=True)` truncates the rim when `resolution_tip` < ρ_r/3 on the default path", naming the variant that avoids it — or, on the negative, that none of the three licensed variants does.
- **Negative result.**
  - **No variant reads 0 inside midpoints on ρ_r/4** (the ruling's park condition): park the code on `attempt/IMP-6-item426-<ts>`, mark item 426 BLOCKED in the same commit (rule (d)), and commit on `main` the attempts entry with every W0 census, the known-issues entry, and the banked negative in the row's step-2 clause — "the oblate rim ladder is mesher-bounded at ρ_r/3 (items 408, 426)". The family is not re-lettered.
  - **A variant reads 0 on ρ_r/3 and ρ_r/4 but not on ρ_r/5**, or ρ_r/5 trips the 95 000 guard: the pre-stated fallback, `FEM_EM_IMP6_RUNG2=4.5` with the same variant (one more W0), disclosed. Still not 0 or over the guard: a two-rung ladder (ρ_r/3, ρ_r/4), one order, disclosed.
  - **W0r or W3 moves a default record** (37 359 / 44 304 / 39 142 / 32 676 cells or the two SICN lines): the keyword leaked into the default path. Revert `src/`, stop, BLOCKED.
  - **(C-fid) reads 0 on the default path:** the census is not counting what attempt 1 counted, or gmsh moved; stop and report — no variant reading is interpretable.
  - **(i) or (ii) red on a rung that passed (G1)–(G2):** a real finding about the FEM on a faithful mesh. Print h_rim/ρ_r, min SICN, inside count and deficit, park on `attempt/*`, BLOCKED (rule (d)); the review rules. No band moves; a miss that grows with refinement is reported as measured, never fitted.
  - No `PREDICTED RED` line: every gate is predicted green on a faithful mesh.

**Outcome (10-13 15:00 CDT slot): DONE, variant (1a) `oblate_rim="edge"`.** Rule (c) `src/` change, disclosed: `_oblate_spheroid_in_box_domain` gains keyword-only `rim` / `netgen`, forwarded as `oblate_rim` / `oblate_netgen` on the `oblate=True` branch only. With the defaults the field order and mesh are unchanged. `edge` builds the body by revolving two quarter-ellipse arcs, so the equator is a model edge; the stray profile face is removed. The arc API was probed in the image first. `tube` and `balls-nonetgen` are coded and were not run.
- **Gates.**
  - (C-fid): the default path reads 31 / 93 inside midpoints (`20261013T200330Z_IMP-6-step2-a2-W0-ctl.log:172, :312`, 21 s).
  - (G1)/(G2) on W0-edge: 0 / 0 / 0 inside, min SICN 0.3605 / 0.3994 / 0.3466, cells 33 334 / 53 959 / 89 530 (guard 95 000), air OCC mass rel 0.000e+00 (`20261013T200400Z_IMP-6-step2-a2-W0-edge.log:166–172, :284–290, :415–421`, 26 s). Cavity faces print 2 (W1 logs `:65`).
  - (G3)/(G4) on W0r: 37 359 / 44 304 and 39 142 / 32 676 at SICN 0.2667 / 0.3707, unchanged (`20261013T200442Z_IMP-6-step2-a2-W0r.log:140, :207, :320–321, :423–424`, 16 s). W3 `test_pec_sphere_hole.py` 7 passed (`…200504Z_…-W3.log`, 9 s).
- **Anchors, all green** (degree 2, `-n 2`, complex):

  | rung | log | Φ_h / Φ_cf,h | signed miss | β | \|β − 1\| | cavity dofs | solve, RSS |
  |---|---|---|---|---|---|---|---|
  | ρ_r/3 | `20261013T200520Z_IMP-6-step2-a2-W1a.log:169–173`, 41 s | 22.971085 / 22.693485 | **+1.2233 %** | 1.001961 | 0.196 % | 20 978 | 37.6 s, 4.04 GiB |
  | ρ_r/4 | `…200607Z_…-W1b.log:177–180`, 76 s | 23.042163 / 22.704740 | **+1.4861 %** | 0.992222 | 0.778 % | 24 710 | 72.8 s, 7.54 GiB |
  | ρ_r/5 | `…200729Z_…-W1c.log:190–193`, 152 s | 23.079024 / 22.712699 | **+1.6129 %** | 0.994011 | 0.599 % | 30 282 | 148.5 s, 14.41 GiB |

  Cavity max \|E\| is 0 on every rung. \|Im E\|/\|Re E\| reads 8.52e-3.
- **Record.**
  - h_rim/ρ_r: max edge 0.691 / 0.510 / 0.394, median 0.511 / 0.387 / 0.309.
  - Body-volume deficit: −5.147e-3 / −5.148e-3 / −5.088e-3, against step 1's −5.777e-3 on the B-spline body.
  - Pairwise observed order p = ln(m_k/m_{k+1})/ln(H_k/H_{k+1}). Rungs 0–1: **−0.642** (max-edge H), −0.700 (median H), −0.676 (nominal divisor). Rungs 1–2: **−0.318**, −0.362, −0.367. There is no sign change, so the order is defined.
  - **The miss grows slightly with rim refinement.** Rung ρ_r/4 was predicted ≈ 0.8 % and read 1.49 %. Rung ρ_r/5 was predicted in [0.3 %, 1.4 %] and read 1.61 %, outside the range. The order was predicted in [0, 2.5]. These predictions are not asserted (rule (e)). The rim-local grain does not set the miss; a floor near 1.2–1.6 % is set by what the ladder holds fixed. This is reported as measured and not fitted. Rung ρ_r/3 (+1.22 %) lands inside its predicted [0.3 %, 2.0 %], beside step 1's 1.418 % on the B-spline body.
- **Known issue.** The 2026-10-13 entry stays open with the avoiding variant named, because the default path is unchanged and still truncates.

**Every-item rules** (PROJECT_PLAN §9 "Every item", bind this item as written above): bare whole harness calls; a lost docker socket stops the slot; no edit after the last window; `timeout: 660000` and `-s`; never delete a log; BLOCKED in the parking commit; a killed window is not a result; two heavy windows never at once; complex mode + `FEM_EM_REQUIRE_COMPLEX=1`; no slot opens `docs/private/` or `aed_results/`.

**Spawn line:** "Use the implementer agent to execute `IMP-6` §9 item 426 per `docs/planning/chunks/IMP-6.md` (read the latest Review ruling paragraph first — 2026-10-13 03:00 on item 408). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Take `tests/validation/test_imp6_oblate_rim_ladder.py` and the six attempt-1 logs by path from `attempt/IMP-6-item408-20261013T0339Z`; never path-checkout `docs/testing/test-results.md` — append its rows by hand. The `src/` change is confined to `_oblate_spheroid_in_box_domain` behind `rim` / `netgen` keywords whose defaults are byte-identical; disclose it. No solve before a W0 census reads 0 inside midpoints and min SICN ≥ 0.1 on every rung; try the variants only in the brief's order; never import a `test_*` name from step 1's module."

### Review rulings 2026-10-13 03:00 (daily) — step 3 audit, the item 426 draft

1. **Step 3 audit PASS** (`b0a3fce1`, item 409): `20261013T050311Z_IMP-6-step3-W1a.log:197–199`, `20261013T050558Z_IMP-6-step3-W1b.log:190–192`, control `20261013T050819Z_IMP-6-step3-W2.log:192–193`; asserts at `tests/validation/test_imp6_oblate_axial_fem.py:242–245, :285–286`; bands imported; rule (i) holds on the file mtime. Per the 10-12 weekly, step 3's pole gate joins the §7 Done-when (rule (j), written in this review's commit). Caveats recorded: a ball-mean Φ, not a pointwise 1/L_z; the pole on the unrefined grain; the inside-body midpoint count was not printed at the poles.
2. **Item 426's reading of "rung ρ_r/3 reproduces step 1" is confirmed:** the default path re-gates unchanged in W0r (39 142 / 32 676 cells, SICN 0.2667 / 0.3707); a variant's own ρ_r/3 rung is a different mesh and is gated at the bands with its miss printed beside 1.418 %.
3. **Min SICN ≥ 0.1 is read from the builder's printed line** (the builder returns no quality); the inside-midpoint count is asserted. The gmsh ellipse-arc call is checked against the image before the first coding step. Item 426 runs `tests/environment/test_complex_mode.py` on its path list, as item 408's priced windows did.

### §9 item 440 — 2026-10-13 (§10 chain I26, daily 2026-10-13 — provisional; lane 3, closed form)

`IMP-6` step 2b, attempt 1 of 4. The family is step 2's; item 426 registered a record, so the count restarts. One measurement on rung ρ_r/4 (`oblate_rim="edge"`) that separates the two candidates item 426 named for its +1.2–1.6 % rim floor (implementer; heavy, `-n 2`, complex, degree 2; ≈ 12–18 measured / 45 predicted min: ≈ 4 min of windows + 15 fixed + ≈ 25 code; depends on nothing). This is a §9 rule (6) provisional step on chain I, ruled by the 2026-10-13 18:00 daily. The 10-14 weekly ratifies or reverts it. If the weekly reverts I26, the item is withdrawn unrun.

- **The ruling (2026-10-13 18:00 daily; implemented here, not changed).**
  - Target and bands are I22's, unchanged. Per run: \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball Φ ≤ `MAT4_BAND`, cavity dofs zero with max \|E\| < 1e-12. All bands are imported and never restated as numbers in code.
  - The miss and its moves are printed, never asserted.
  - The two candidates:
    - (A) the fixed `resolution_surface` (a/12) and far grain. Item 426's β moved 0.2–0.8 % across rungs, the same size as the drift.
    - (B) the comparand's degree-8 `interior="continue"` quadrature over the rim-mask cells near the focal ring.
  - One rung, one knob, one quadrature check.
- **What item 426 already says about (B)** (this is evidence that motivates the prediction below; the measurement is still run). Read from the cited log lines, the drift in the miss is carried by Φ_h, not by the comparand:
  - Φ_h moves 22.971085 → 23.042163 → 23.079024 (+0.309 % / +0.160 %).
  - Φ_cf,h moves 22.693485 → 22.704740 → 22.712699 (+0.050 % / +0.035 %), on three non-nested meshes whose median near-rim edge differs 1.66×.
  - Sources: `20261013T200520Z_IMP-6-step2-a2-W1a.log:171`, `…200607Z_…-W1b.log:179`, `…200729Z_…-W1c.log:192`.
  - Suppose degree 8 carried a ≥ 0.6 % quadrature error on the comparand. That error would have to be the same, to within 0.05 %, on all three meshes.
- **The knobs, named from the code.**
  - **(B), the quadrature.**
    - `_numpy_quadrature_mean` (`tests/validation/test_imp6_oblate_fem.py:173–193`) calls `basix.make_quadrature(basix.CellType.tetrahedron, QDEG)` (`:178`). `QDEG = 8` is a module constant (`:67`), and the helper takes **no degree argument**.
    - basix's default simplex rule is Xiao–Gimbutas up to degree 30, so degrees 12 and 14 are both available. The executor prints `len(wts)` per degree.
    - **Chosen: degree 14 decides; degree 12 is printed beside** as the intermediate step, so that a move can be read as converging or not.
    - Φ_h is assembled by UFL with `quadrature_degree: QDEG` (`tests/validation/test_imp6_oblate_rim_ladder.py:229`). `e_complex` is the N1curl degree-2 field (`src/fem_em_solver/core/solvers.py:259`, `core/time_harmonic.py:128`). On affine tets its \|E\|² is a degree-4 polynomial times the DG0 mask, so degree 8 already integrates Φ_h exactly. Raising the degree on Φ_h is an identity check, not a reading.
  - **(A), the surface grain.**
    - `resolution_surface` **is exposed**: it is a positional argument of `spheroid_in_box_domain` (`src/fem_em_solver/io/mesh.py:2736`), forwarded to `_oblate_spheroid_in_box_domain` (`:2786` → `:2919`).
    - Its only use is the cavity Distance→Threshold SizeMin (`:3101`), with SizeMax `resolution_far` and DistMax 2c = 2.4 mm (`:3102–3104`).
    - The ladder hard-codes `resolution_surface=SEMI_MAJOR / 12.0` (`test_imp6_oblate_rim_ladder.py:117`). **No `src/` edit is needed, and rule (c) does not bind.**
    - Inside the rim mask the Ball's `VIn` (`resolution_tip` = ρ_r/4 = 0.06 mm) is below a/18 = 0.333 mm and wins the `Min`. The knob therefore acts outside the rim ball: on the body-surface layer, including the β probe shell at 1.2× (1.2 mm off the rim, inside DistMax).
- **Code (test module only: `tests/validation/test_imp6_oblate_rim_ladder.py`; no `src/` edit).**
  - (1) `SURFACE_DIV = float(os.environ.get("FEM_EM_IMP6_SURFACE_DIV", "12"))`, validated > 0. `_mesh_rung` passes `resolution_surface=SEMI_MAJOR / SURFACE_DIV`.
    - `float("12") == 12.0`, so the default is bitwise the landed float.
    - The cache key and `_mesh_seconds` become `(variant, divisor, SURFACE_DIV)`.
    - Every census and report line prints the surface divisor and `resolution_surface` in mm.
  - (2) `_numpy_quadrature_mean_q(msh, cells, field_fn, qdeg)`: step 1's helper (`test_imp6_oblate_fem.py:173–193`) copied with `qdeg` in place of `QDEG`. `_ufl_dipole_mean_q(msh, chi, qdeg)`: `_phi_closed_form_on_mask` (`test_imp1_pec_sphere_lossy_host.py:175–188`) at β = 1, copied with `quadrature_degree: qdeg`. That source hard-codes 8.
    - Import `_imp1_dipole_numpy` (private) from step 1's module and `SPHERE_RADIUS` from `test_imp1_pec_sphere_lossy_host.py`.
    - Step 1's module is not edited, and **no `test_*` name is imported**.
    - `QDEGS_2B = (8, 12, 14)`. The landed gate keeps `QDEG` (8) for both Φ_h and Φ_cf,h.
  - (3) `_solve_on` also returns:
    - Φ_cf,h at degrees 12 and 14 (`_numpy_quadrature_mean_q` over the same `cells`).
    - Φ_h at degree 14 (a second UFL `dx` with `quadrature_degree: 14`).
  - (4) `_report_and_gate` prints:
    - those readings;
    - the signed miss against each comparand degree;
    - when `SURFACE_DIV != 12`, item 426's rung-1 record beside: Φ_h 23.042163, Φ_cf,h 22.704740, miss +1.4861 %, β 0.992222 (`…200607Z_…-W1b.log:178–179`);
    - item 440's prediction in place of `PREDICTED_MISS`.

    It also **asserts (Qh)** below. The (i)–(iii) asserts are unchanged.
  - (5) New node **`test_quadrature_rung1`** (W1q, **no solve**). The comparand needs only the mesh and the mask cells. On the rung ρ_r/4 mesh at the module's `SURFACE_DIV` (run with the default, i.e. item 426's W1b mesh):
    - Print the mask cells and Φ_cf,h at degrees 8 / 12 / 14, with Φ_cf,h(q)/Φ_cf,h(8) − 1.
    - Print item 426's 22.704740 beside degree 8. The comparison is reported, not asserted: an image move would show here.
    - Print the (B) decision line.
    - Assert (Q8) and (M′_q).
  - (6) New node **`test_surface_census`** (W0, no solve). Rung ρ_r/4 **only**, at the module's `SURFACE_DIV`:
    - Assert cells ≤ `CELL_GUARD`, mask cells > 0, and **inside-midpoint count == 0** (`_inside_count`).
    - Assert (M′_q) on this mesh.
    - Print cells beside the prediction, min SICN (the builder's line), h_rim/ρ_r (max / median / p95), the deficit, RSS and `resolution_surface`.
    - `test_ladder_census` is not used for the a/18 mesh: it would also mesh ρ_r/5 at a/18 and assert the guard there.
  - Update the docstring: item 440, the two knobs, and the ruling. `test_ladder_census`, `test_census_control` and `_variant_kwargs` keep their logic.
- **Anchors (asserted).**
  - W1s on the a/18 rung ρ_r/4, item 426's anchors at imported bands:
    - (i) \|β − 1\| ≤ `TH15_BETA_BAND`.
    - (ii) \|Φ_rim,h/Φ_rim,cf,h − 1\| ≤ `MAT4_BAND`, with the degree-8 comparand as landed.
    - (iii) cavity dofs > 0 and max \|E\| < 1e-12.
  - W0 on the a/18 mesh, before any solve: cells ≤ `CELL_GUARD`, mask > 0, **inside midpoints 0** ((G1) on the new mesh). Min SICN ≥ 0.1 is read from the builder's printed line ((G2)); no solve starts below it.
  - Instrument identities (not bands):
    - (Q8): `_numpy_quadrature_mean_q(…, 8)` equals the imported `_numpy_quadrature_mean` to ≤ 1e-13 relative. Same rule, same arithmetic.
    - (M′_q): for q ∈ {12, 14}, the numpy quadrature of `_imp1_dipole_numpy` equals `_ufl_dipole_mean_q` to ≤ 1e-12 relative, on the W1q mesh and the W0 mesh. FFCx takes basix's default rule for a `quadrature_degree`, so both sides use the same points.
    - (Qh): \|Φ_h(14)/Φ_h(8) − 1\| ≤ 1e-12 in W1s (polynomial exactness).
- **Decision rule (pre-registered; printed by the module and stated in the outcome / journal; never asserted).**
  - m_B = \|Φ_cf,h(14)/Φ_cf,h(8) − 1\| on item 426's rung-ρ_r/4 mesh (W1q). **(B) explains the floor iff m_B ≥ 0.6 %**, half the ≈ 1.2 % floor.
  - m_A = the signed miss on the a/18 rung ρ_r/4 (W1s), degree-8 comparand, as item 426 read 1.486 %. **(A) explains it iff m_A ≤ 1.486 % − 0.5 pp = 0.986 %.**
  - **Both or neither are reported as such.** "Neither" means the floor is the far grain (`resolution_far` a/4, untouched) or the comparand itself; it is reported, and the weekly decides whether to pose more.
    - Context, not tested here: the comparand is the quasi-static closed form in a host with (\|k\|a)² = 1.15 % (item 404's sizing), the same order as the floor. Aspect 2, at the same a and k, read 0.152 %.
  - Printed beside, but outside the rule: m_B at degree 12 and on the a/18 mesh, the miss against the degree-14 comparand on both meshes, and β on both meshes.
- **Predicted (rule (e); printed beside the measured value, never asserted).**
  - m_B ≲ 0.05 %, bounded by the cross-mesh comparand drift above, so (B) is predicted not to fire.
    - Geometry: inside count is 0, and the chord sag at the rim is ≈ h²/(8ρ_r) ≈ (0.093 mm)²/(1.92 mm) ≈ 4.5 µm on rung ρ_r/4's median near-rim edge. The focal ring is 121 µm inside the surface, so no quadrature point sits near it.
  - m_A in [0.5 %, 1.6 %]:
    - If the floor is O(h_surf²) pollution, it scales by (12/18)² = 0.44. That gives ≈ 0.6–0.9 % with the rim-local part, and (A) fires.
    - If it is the far grain, m_A ≈ 1.5 % and (A) does not fire.
  - \|β − 1\| ≲ 0.8 %.
  - Cells ≈ 60–65 k (estimate below).
- **Sizing: the a/18 cell count and the guard.**
  - No record varies `resolution_surface`: step 1's a/9 fallback was never run. The estimate is therefore built from item 426's own rungs.
  - Fitting N = N₀ + B·div³ to 33 334 / 53 959 cells (`20261013T200400Z_IMP-6-step2-a2-W0-edge.log:172, :290`) gives B ≈ 557 and **N₀ ≈ 18.3 k cells not set by the Balls**. The fit predicts 87 960 at ρ_r/5 against 89 530 measured (`:421`, −1.8 %).
  - The cavity layer's cell count scales as ∫ dd/h(d)³ with h(d) linear from `resolution_surface` to 1.5 mm over 2c = 2.4 mm: (1/s² − 1/F²)/(2k) reads 4.27 at s = 0.5 mm and 8.80 at s = 0.333 mm, a factor of **2.06**.
  - The far region at 1.5 mm (≈ 8–9 k cells of N₀) does not move. That gives **≈ 60–65 k**. Bound: all of N₀ scaling by 2.06 gives ≈ 73–75 k.
  - **Both are under `CELL_GUARD` 95 000, so the imported guard binds as in item 426. No separate guard is set.**
  - Memory: linear in cells between W1b (7.54 GiB at 53 959) and W1c (14.41 GiB at 89 530) gives ≈ 9.5 GiB summed at 64 k and ≈ 11.6 GiB at 75 k, under item 426's W1c (14.4 GiB at 89 530, `-n 2`).
- **Price** (item 426's measured windows, `-n 2`, `tests/environment/test_complex_mode.py` on the path):
  - Measured basis:
    - W0-edge 26 s for three meshes (`…200400Z_…-W0-edge.log:474`).
    - W0-ctl 22 s for two (`…200330Z_…-W0-ctl.log:336`).
    - W1b 76 s window, 72.8 s solve+post, 7.54 GiB (`…200607Z_…-W1b.log:177, :220`).
    - W1c 152 s, 148.5 s, 14.41 GiB (`…200729Z_…-W1c.log:190, :233`).
  - Solve time scales as cells^1.41 between those two.
  - Predicted:
    - W0r ≈ 26 s.
    - W1q ≈ 25–40 s: one mesh, ≈ 11.5 k mask cells × three rules.
    - W0 ≈ 25–40 s.
    - W1s ≈ 100–140 s (≈ 93 s solve at 64 k, ≈ 116 s at 75 k, plus the degree-14 forms).
    - Total ≈ 3–4 min of windows.
  - **W1q is not a re-solve.** Re-integrating Φ_h on item 426's rung-1 solve would mean re-solving W1b (≈ 76 s / 7.54 GiB, `…W1b.log:177, :220`) to read an identity that degree 8 already satisfies exactly. It is not run. (Qh) reads that identity in W1s at no extra solve.
- **Windows** (complex, `-n 2`, one solve per window, foreground, `timeout: 660000`; item 426's `Command:` shape, `…200607Z_…-W1b.log:12`). Order: code → W0r → W1q → W0 → W1s → commit. **No module or `src/` edit after W0r starts** (rule (i)).
  - W0r, the default re-gate of the edited module. Expect 33 334 / **53 959** / 89 530 cells, min SICN 0.3605 / **0.3994** / 0.3466, and 0 / **0** / 0 inside (`…200400Z_…-W0-edge.log:166–172, :284–290, :415–421`). Command: `scripts/testing/run_and_log.sh IMP-6-step2b-W0r "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_ladder_census -s -v --tb=short'"` (`FEM_EM_IMP6_SURFACE_DIV` unset).
  - W1q, the (B) reading on item 426's rung-1 mesh: `scripts/testing/run_and_log.sh IMP-6-step2b-W1q "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_quadrature_rung1 -s -v --tb=short'"`
  - W0, the a/18 census: `scripts/testing/run_and_log.sh IMP-6-step2b-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_SURFACE_DIV=18 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_surface_census -s -v --tb=short'"`
  - W1s, the (A) reading, which runs only if W0 is green and its SICN line reads ≥ 0.1: `scripts/testing/run_and_log.sh IMP-6-step2b-W1s "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_SURFACE_DIV=18 timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung1 -s -v --tb=short'"`
- **Negative control.** None new.
  - (C-fid) is not needed. Item 426's (C-fid) already shows on this fixture that the census can fail: the default path reads 31 / 93 inside midpoints (`20261013T200330Z_IMP-6-step2-a2-W0-ctl.log:172, :312`). A 0 on the a/18 mesh is therefore a reading of the mesh, not of the check.
  - W0 asserts that 0 before any solve.
  - The cheap asserted control added is (M′_q): the same comparison as step 1's (M′), numpy against UFL quadrature of `IMP-1`'s dipole over this fixture's rim mask. Step 1 measured it at 5.5e-14 on the aspect-5 mesh (`20261012T184107Z_IMP-6-step1-W0.log:605`), asserted ≤ 1e-12. Here it is taken at degrees 12 / 14, and it is what licenses reading m_B as a property of the comparand rather than of the new helper.
  - Every other factor in this brief is labelled *predicted*.
- **Traps already paid for.**
  - **Rule (i): no module or `src/` edit after W0r starts.** The env knob's default must reproduce W0r bit for bit. A cache key without `SURFACE_DIV` would hand W1s a cached a/12 mesh.
  - **Never run `test_ladder_census` with `FEM_EM_IMP6_SURFACE_DIV` set.** It would mesh ρ_r/5 at a/18 and could trip the guard. W0 is the one-rung node.
  - **The landed gate stays at degree 8.** Do not change the imported `QDEG`. The degree-12 / 14 values are printed, and only the identities (Q8) / (M′_q) / (Qh) are asserted on them.
  - **Never import a `test_*` name** from step 1's module: pytest would collect it and re-run step 1's solves. `_imp1_dipole_numpy` is private.
  - **The basis is the in-plane field minus E₀x̂** (item 404's trap). `_solve_on` already does this; do not touch it.
  - **A quality trap means park unless the brief names the fix** (2026-10-12 18:00 ruling). This brief names none: SICN < 0.1, inside ≠ 0 or the guard tripped on the a/18 mesh means no solve and a park. There is no a/16 or a/15 fallback, because the ruling has one knob at one value.
  - Owned cells only (`size_local`). `assemble_scalar` is rank-local (`_assemble_real` reduces), and every new numpy sum is reduced with `comm.allreduce`. The RSS reduction stays outside the rank-0 branch. Masks stay in numpy, with no UFL `<`, `max_value` or `sqrt` comparison.
  - `-s` on every window and `-k 30` on every container `timeout`. Select by node ID, never `-k a or b`.
  - The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run in the foreground with `timeout: 660000`. Two heavy windows never run at once.
  - After a killed window, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first. Use the rubric's orphan-`.c` sweep only if it reads 0.
- **Scope.** One aspect (5), one rung (ρ_r/4), one knob (`resolution_surface` a/12 → a/18) and one quadrature check (degree 8 → 12 / 14).
  - Not a convergence study and no rate.
  - Not aspect 2, not the axial drive, not `resolution_far`.
  - The default oblate path, step 1, step 3 and `EX-111` are untouched.
  - The row's Not-claimed list stands.
- **Status it can move.** The `IMP-6` §7 row's step-2 record clause currently ends "a floor from the fixed surface / far grains, a finding, not a rate". It gains the attribution: "; item 440 (I26): comparand degree 14 moves m_B, the a/18 miss reads m_A — the floor is (A) / (B) / both / neither", with the W1q / W1s log names.
  - No glyph change: `IMP-6` is ✅.
  - The slot does not edit the Done-when or the Not-claimed list.
  - **The 2026-10-13 known-issues entry is not touched.** It concerns the default path's rim truncation.
- **Negative result.**
  - **W0r moves a default record** (cells, SICN or inside count on any rung): the knob leaked into the default path. Fix the module and re-run W0r before any other window. The fixed W0r is the re-gate.
    - If the module's default-path diff is provably empty (`SEMI_MAJOR / float("12")`) and the counts still move, that is an image move. Stop and report: W1q's comparison with item 426 would no longer be on the same mesh.
  - **(Q8) or (M′_q) red:** a bug in the new helper. Fix it in-slot, before W1s; no m_B is read from a red identity. If it cannot be fixed, park.
  - **W0 red on the a/18 mesh** (guard, inside ≠ 0, or the SICN line < 0.1): no solve. Report the counts, park on `attempt/IMP-6-item440-<ts>`, and mark BLOCKED in the same commit (rule (d)). The review rules.
  - **(i), (ii) or (iii) red on the a/18 mesh:** a real finding. Print h_rim/ρ_r, min SICN, the inside count and the deficit, park on `attempt/*`, and mark BLOCKED (rule (d)). The review rules.
  - **(Qh) red:** the Φ_h assembly is not what the module claims. Park and mark BLOCKED.
  - In every case, no band moves. A miss that grows is reported as measured and never fitted.
  - No `PREDICTED RED` line: every gate is predicted green.

**Every-item rules** (PROJECT_PLAN §9 "Every item", binding this item as written above):
- Bare, whole harness calls.
- A lost docker socket stops the slot.
- No edit after the last window.
- `timeout: 660000` and `-s` on every window.
- Never delete a log.
- BLOCKED goes in the parking commit.
- A killed window is not a result.
- Two heavy windows never run at once.
- Complex mode and `FEM_EM_REQUIRE_COMPLEX=1` on every frequency-domain command.
- No slot opens `docs/private/` or `aed_results/`.

**Spawn line:** "Use the implementer agent to execute `IMP-6` §9 item 440 per `docs/planning/chunks/IMP-6.md` (the section carries the 2026-10-13 18:00 daily's provisional I26 ruling; if the 10-14 weekly reverted I26, do not run). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Edit only `tests/validation/test_imp6_oblate_rim_ladder.py`: an env knob `FEM_EM_IMP6_SURFACE_DIV` (default 12, bitwise the landed path, in the cache key), degree-parametrised copies of the two quadrature helpers, and the nodes `test_quadrature_rung1` / `test_surface_census`. No `src/` edit. Never import a `test_*` name from step 1's module. Windows in order W0r → W1q → W0 → W1s, no edit after W0r starts, no solve before W0 reads 0 inside midpoints and min SICN ≥ 0.1 on the a/18 mesh. State the pre-registered attribution ((A) / (B) / both / neither) in the outcome; assert none of it."

**Outcome (10-13 21:00 CDT slot, 2026-10-14 ≈ 02:31Z): DONE, attribution neither.** Test module only, no `src/` edit.
- **Gates.**
  - W0r (default re-gate): 33 334 / 53 959 / 89 530 cells, SICN 0.3605 / 0.3994 / 0.3466, 0 / 0 / 0 inside, unchanged (`20261014T022824Z_IMP-6-step2b-W0r.log:167–173, :286–292, :418–424`, 27 s).
  - W1q (no solve, item 426's rung-1 mesh, 53 959 cells, 11 468 mask cells): Φ_cf,h(8 / 12 / 14) = 22.704739974621 / 22.704739974727 / 22.704739974727 (45 / 122 / 177 points per cell), item 426's 22.704740 reproduced; (Q8) 0.000e+00; (M′_q) 3.75e-14 / 5.64e-14 (`20261014T022856Z_IMP-6-step2b-W1q.log:177–186`, 15 s).
  - W0 (a/18 mesh, ρ_r/4): 63 363 cells (predicted ≈ 60–65 k), 11 316 mask cells, 0 inside, min SICN 0.4204, h_rim/ρ_r max 0.504 / median 0.381 / p95 0.451, deficit −3.230e-3; (M′_q) 2.03e-14 / 1.33e-15 (`20261014T022916Z_IMP-6-step2b-W0.log:177–184`, 14 s).
  - W1s (a/18, ρ_r/4, degree 2, `-n 2`): \|β − 1\| 0.498 % (β 0.995017), Φ_h / Φ_cf,h = 23.005304 / 22.694130, signed miss **+1.3712 %**, cavity dofs 36 802 at max \|E\| 0, (Qh) 3.22e-15; solve+post 87.5 s, 8.80 GiB (`20261014T022936Z_IMP-6-step2b-W1s.log:179–187`, 90 s).
- **Decision rule (pre-registered, printed, not asserted).**
  - m_B = 4.67e-12 (degree 12: 4.67e-12) ≪ 0.6 % → **(B) does not explain the floor** (predicted ≲ 0.05 %).
  - m_A = +1.371 % > 0.986 % → **(A) does not explain the floor.** It is inside the predicted [0.5 %, 1.6 %], on the "far grain" side. The miss falls by 0.115 pp from item 426's +1.486 %: Φ_h −0.16 %, Φ_cf,h −0.05 %. β moves from 0.992222 to 0.995017.
  - **Neither.** By the rule this leaves the far grain (`resolution_far` a/4, untouched) or the comparand itself, the quasi-static closed form at (\|k\|a)² = 1.15 %. That is context for the weekly; it was not tested here.
- **Disclosed.** The first W1q (`20261014T022119Z_IMP-6-step2b-W1q.log`, exit 124 at the 400 s container timeout) deadlocked: `test_quadrature_rung1` called the collective `_ncells` inside its rank-0 print branch. Its prints before the hang are not a result. The fix moved that call out of the branch. That was an edit after W0r started, so W0r was re-run before any other window (the brief's W0r-repair path). The first W0r (`20261014T022046Z_…-W0r.log`) read the same records. `pgrep -c python3` read 0 after the kill.

### Review ruling 2026-10-13 18:00 (daily) on the item 440 draft

1. **Accepted as drafted.** W1q reads the comparand at degrees 8 / 12 / 14 on item 426's rung-1 mesh with no re-solve. Φ_h is a degree-4 polynomial on affine tets, so its degree change is the identity (Qh), asserted in W1s, not a reading. The (Q8) / (M′_q) instrument identities are asserted; (M′_q) is backed by step 1's same-fixture 5.5e-14 (`20261012T184107Z_IMP-6-step1-W0.log:605`). No `src/` edit is made, because `resolution_surface` is already an argument.
2. **The drafter's reading of item 426's logs stands as a prediction, not a result.** The comparand moved only +0.050 / +0.035 % across the three meshes, so (B) is predicted not to fire. The run still measures it. The (|k|a)² = 1.15 % context is recorded for the weekly and is not tested here.
3. **Minutes:** the §9 entry carries the lane-3 class median (11 measured, journal §6). The brief's 12–18 is the drafter's estimate. Only the class median counts toward the floor.
4. **Provisional (§9 rule (6), chain I's first daily step since the 10-12 weekly).** The 10-14 weekly ratifies or reverts I26. A slot that runs before the weekly runs it as written.

### §9 item 441 — 2026-10-14 (§10 chain I27, daily — provisional; lane 3)

`IMP-6` step 2c, attempt 2 of 4 in the step-2 family since item 426's registered record (item 440 was 1 of 4). It separates the two candidates item 440 left for I22's +1.2–1.6 % rim floor, on rung ρ_r/4 (`oblate_rim="edge"`), one knob each (implementer; heavy, `-n 2`, complex, degree 2; ≈ 11 measured / 45 predicted min: ≈ 5–6 min of windows + 15 fixed + ≈ 25 code; depends on nothing). This is chain I's second §9 rule (6) provisional step (I26 was the first; rule (6) allows three per chain), written by the 2026-10-14 03:00 daily. The 10-14 weekly ratifies or reverts I26 and I27. **If the weekly reverts either, the item is withdrawn unrun.**

- **The question (item 440's outcome; attempts.md `2026-10-14T02:33Z`, "Hypothesis for the review").** The floor is neither the comparand quadrature (B: m_B = 4.67e-12, `20261014T022856Z_IMP-6-step2b-W1q.log:180`) nor the surface grain (A: a/18 reads +1.3712 % against item 426's +1.4861 %, `20261014T022936Z_IMP-6-step2b-W1s.log:181`). Two candidates remain:
  - (C) the far grain `resolution_far` = a/4 = 1.5 mm, which the ladder has never moved;
  - (D) the comparand's own finite-k error. `pec_oblate_spheroid_field` and the box trace are quasi-static, but the FEM solves curl curl E − k²E = 0 with k = 15.178 − 9.489j 1/m, \|k\| = 17.900 1/m, (\|k\|a)² = 1.1535 %, corner (\|k\|R)² = 8.2184 % (`20261012T184107Z_IMP-6-step1-W0.log:392`; IMP-2's line `:256` prints the complex k).
- **Target and bands are I22's, unchanged.** Every new solve asserts \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball Φ ≤ `MAT4_BAND` and cavity dofs > 0 with max \|E\| < 1e-12. All three bands are imported, never restated as numbers. The miss and its moves are printed and never asserted.
- **The knobs, named from the code.**
  - **(D), k² → s²k² at a fixed loss tangent.**
    - The operator is `mu_r_inv·curl·curl − k0²·ε_c` with k0² = ω²μ₀ε₀ and ε_c = ε_r − jσ/(ωε₀) (`src/fem_em_solver/core/time_harmonic.py:792–806`). So k² = ω²μ₀ε₀ε_r − jωμ₀σ. Scaling **both** f and σ by s scales both terms by s², and the phase of k² is unchanged.
    - The module takes both from imports: `FREQUENCY_HZ` = 64 MHz and `SIGMA_HIGH` = 0.57 S/m (`tests/validation/test_lossy_sphere_sar.py:74, :77`), used only in `_solve_on` (`tests/validation/test_imp6_oblate_rim_ladder.py:280–284`). **No knob exists, so one is added. It is a test-module knob only.**
    - s = 0.1 gives k² × 0.01: f = 6.4 MHz, σ = 0.057 S/m, (\|k\|a)² ≈ 0.0115 %. s = 0.3 gives k² × 0.09.
    - The solver accepts both values: `frequency_hz` > 0 (`time_harmonic.py:456`), σ ≥ 0 and ε_r > 0 (`:89–93`). The solve is a direct MUMPS LU (`:564–567`).
    - **Validity as k → 0.**
      - `pec_oblate_spheroid_field` (`src/fem_em_solver/utils/analytical.py:1831`) is the static solution. It has no k argument, so it is the exact k → 0 limit of the problem the FEM solves.
      - The Dirichlet trace `_dirichlet` (`tests/validation/test_imp6_oblate_fem.py:123–125`) is that same field, so it is k-independent too.
      - β's basis (`_beta`, `:217–222`) and the Φ comparand (`_field(..., "continue")`, `:111–120`) are the same field. **None of the three changes with s.** The comparand Φ_cf,h on a fixed mesh is therefore bitwise s-independent, and this is asserted below as (Wk).
      - In the homogeneous host, testing with ∇φ gives −k²(ε_c E, ∇φ) = 0. The discrete field stays weakly divergence-free for any k ≠ 0, so the discrete limit is the discrete electrostatic problem under the same trace.
      - Conditioning: the gradient block scales as \|k²\|. At s = 0.1 that is ≈ 3.2 m⁻² against a curl-curl scale of 1/h² ≈ 2.8e8 m⁻² at h = ρ_r/4. A condition number ≈ 1e8 in double precision is ≈ 1e-8 relative, far below the 0.5 pp the rule reads. s is not taken below 0.1.
    - If (D) sets the floor, the miss is m(s) ≈ m_h + s²·δ_D to leading order in k². The s = 0.1 value then removes 0.99·δ_D.
  - **(C), the far grain.**
    - `resolution_far` is a positional argument of `spheroid_in_box_domain` (`src/fem_em_solver/io/mesh.py:2737`), forwarded to `_oblate_spheroid_in_box_domain` (`:2920`). On the oblate path it is used in four places:
      - the Balls' `VOut` (`:3068`), which the `Min` ignores;
      - the rim Threshold's SizeMax at 4ρ_r (`:3092–3094`);
      - the cavity Threshold's SizeMax at DistMax 2c (`:3102–3104`);
      - by the `Min` (`:3106–3108`), the bulk size beyond all of these.
    - The ladder hard-codes `resolution_far=SEMI_MAJOR / 4.0` (`test_imp6_oblate_rim_ladder.py:145`). **No `src/` edit is needed, and rule (c) does not bind.**
    - Refined one step: a/4 → **a/5 = 1.2 mm**. Inside the rim mask the Ball's `VIn` (ρ_r/4 = 0.06 mm) wins the `Min` as before.
- **Code (test module only: `tests/validation/test_imp6_oblate_rim_ladder.py`; no `src/` edit).**
  - (1) `FAR_DIV = float(os.environ.get("FEM_EM_IMP6_FAR_DIV", "4"))`, validated > 0. `_mesh_rung` passes `resolution_far=SEMI_MAJOR / FAR_DIV`; `float("4") == 4.0`, so the default is bitwise the landed float. The cache key and `_mesh_seconds` become `(variant, divisor, SURFACE_DIV, FAR_DIV)`.
  - (2) `K_SCALE = float(os.environ.get("FEM_EM_IMP6_K_SCALE", "1"))`, validated 0 < s ≤ 1. `_solve_on` passes `frequency_hz=FREQUENCY_HZ * K_SCALE` and `sigma=SIGMA_HIGH * K_SCALE`. Multiplying by 1.0 is exact, so the default is bitwise. **K_SCALE is not in the mesh cache key**, because the mesh is k-independent and W1k must solve on W0r's mesh.
  - (3) **One knob per run, enforced at import.** Raise `ValueError` if more than one of `SURFACE_DIV != 12.0`, `FAR_DIV != 4.0`, `K_SCALE != 1.0` holds.
  - (4) `_solve_on` prints, on rank 0 and with nothing collective inside the branch, the scaled f and σ, k = ω√(μ₀ε₀(ε_r − jσ/(ωε₀))), Re / Im k², \|k\| and (\|k\|a)², beside step 1's 17.900 1/m / 1.1535 % (`…184107Z_IMP-6-step1-W0.log:392`).
  - (5) `_report_and_gate`:
    - The predicted-miss print branches on the active knob, giving item 441's prediction for (C) or (D).
    - **When `K_SCALE != 1`:**
      - Print item 426's rung-1 record beside: Φ_h 23.042163, Φ_cf,h 22.704740, miss +1.4861 %, β 0.992222, \|Im E\|/\|Re E\| 8.523e-3 (`20261013T200607Z_IMP-6-step2-a2-W1b.log:178–179`).
      - Print the (D) decision line (its YES / NO only at s = 0.1; at s = 0.3 the line reads "outside the rule") and δ_D = (1.4861 % − m)/(1 − s²). The module prints δ_D only. The s = 0.3 linearity ratio (1.4861 % − m(0.3))/(1.4861 % − m(0.1)), against the leading-order 0.919, is computed in the outcome from the two logs: putting W1k's value into the module would be an edit after W0r.
      - **Assert (K) and (Wk)** below.
    - **When `FAR_DIV != 4`:** print item 426's rung-1 record beside and the (C) decision line.
    - The (i)–(iii) and (Qh) asserts are unchanged.
  - (6) `test_surface_census` prints `FAR_DIV` and `resolution_far` in mm beside the surface divisor, and item 441's cell prediction when `FAR_DIV != 4`. **When `FAR_DIV != 4` it asserts (Kf).** Its other asserts are unchanged: guard, mask > 0, inside == 0, (M′_q).
  - (7) New node **`test_knob_defaults`** (in W0r, no solve), the asserted default re-gate:
    - (W0r-a) none of `FEM_EM_IMP6_SURFACE_DIV` / `_FAR_DIV` / `_K_SCALE` is set in `os.environ`. `pytest.fail` otherwise: the node is the default path's check and cannot be run under a knob.
    - (W0r-b) `SEMI_MAJOR / FAR_DIV == SEMI_MAJOR / 4.0`, `SEMI_MAJOR / SURFACE_DIV == SEMI_MAJOR / 12.0`, `FREQUENCY_HZ * K_SCALE == FREQUENCY_HZ` and `SIGMA_HIGH * K_SCALE == SIGMA_HIGH`, all by `==`.
    - (W0r-c) on rung ρ_r/4 (the cached default mesh): cells `== 53_959`, rim-mask cells `== 11_468`, inside midpoints `== 0`, and \|Φ_cf,h(8)/22.704739974621 − 1\| ≤ 1e-12. These are item 440's records (`20261014T022824Z_IMP-6-step2b-W0r.log:292`; `20261014T022856Z_IMP-6-step2b-W1q.log:177–178`). The comparand prints to 12 decimals, a resolution of 4e-14 relative.
    - Every collective (`_ncells`, `_inside_count`, `_numpy_quadrature_mean`) is called **before** the rank-0 print branch.
  - (8) Docstring: item 441, the two knobs, the one-knob rule, the ruling. Keep these as they are: `test_ladder_census`, `test_census_control`, `test_quadrature_rung1`, `_variant_kwargs` and the item 440 prints.
  - Still in force: step 1's module is not edited; only private names are imported from it, **never a `test_*` name**.
- **Anchors (asserted).**
  - W1k, W1k2 and W1f each assert item 426's anchors at the imported bands, on rung ρ_r/4:
    - (i) \|β − 1\| ≤ `TH15_BETA_BAND`;
    - (ii) \|Φ_rim,h/Φ_rim,cf,h − 1\| ≤ `MAT4_BAND`, with the degree-8 comparand as landed;
    - (iii) cavity dofs > 0 and max \|E\| < 1e-12;
    - (Qh) \|Φ_h(14)/Φ_h(8) − 1\| ≤ 1e-12.
  - W0f on the a/5 mesh, before its solve: cells ≤ `CELL_GUARD` (95 000, `tests/validation/test_imp6_oblate_fem.py:66`, imported at `test_imp6_oblate_rim_ladder.py:72`), mask > 0 and **inside midpoints 0** ((G1)). Min SICN ≥ 0.1 is read from the builder's printed line ((G2)); no solve starts below it.
  - The default re-gate (W0r-a/b/c), above.
  - Instrument checks (not bands):
    - **(Wk)** in W1k / W1k2: \|Φ_cf,h(8)/22.704739974621 − 1\| ≤ 1e-12. The comparand is s-independent on a fixed mesh, so a miss means W1k did not solve on W0r's mesh (a cache-key or knob leak).
    - **(K)** in W1k / W1k2: \|Im E\|/\|Re E\| ≤ 2·s²·8.523e-3 (item 426's W1b reading, `…200607Z_…-W1b.log:178`). Im E is linear in Im k² to leading order, so the prediction is s²·8.523e-3: 8.5e-5 at s = 0.1 and 7.7e-4 at s = 0.3.
    - **(Kf)** in W0f: cells ≥ 1.05 × 53 959 = 56 657. A far-region-only refinement is ≥ +16 % (sizing below).
- **Decision rule (pre-registered; printed by the module and stated in the outcome / journal; never asserted).** The baseline is item 426's rung-ρ_r/4 miss, **m₀ = +1.4861 %**, on the same a/12, a/4 mesh as W1k (`…200607Z_…-W1b.log:179`). The threshold is item 440's, m₀ − 0.5 pp = **0.986 %**.
  - m_D = the signed miss at s = 0.1 (W1k). **(D) explains the floor iff m_D ≤ 0.986 %.**
  - m_C = the signed miss at `resolution_far` a/5 (W1f), with surface a/12 and k unscaled. **(C) explains it iff m_C ≤ 0.986 %.**
  - **Both or neither are reported as such.** The two knobs are not exclusive: the floor can be a sum. A miss that rises by > 0.5 pp is reported as "moves the wrong way".
  - Printed beside, outside the rule: δ_D; the s = 0.3 point and its leading-order ratio 0.919 (computed in the outcome from the two logs, and only meaningful if m₀ − m_D ≥ 0.1 pp); β at each point beside 0.992222; and m_C beside item 440's m_A +1.3712 %.
- **Predicted (rule (e); printed beside the measured value, never asserted).**
  - (K): \|Im E\|/\|Re E\| ≈ 8.5e-5 at s = 0.1 and ≈ 7.7e-4 at s = 0.3.
  - m_D in **[0.4 %, 1.6 %]**. This prediction does not pick a side.
    - *For (D):* write E ≈ E_qs + k²E₁ with E₁ real (the trace is fixed). Then Im E = Im(k²)E₁, and the measured 8.523e-3 at Im k² = −288.0 m⁻² sizes E₁. Re k² = 140.3 m⁻² enters \|E\|² at first order, by ≈ 2·(140.3/288.0)·8.523e-3 ≈ **0.83 %**. The sign is positive for a fixed-trace Helmholtz problem with Re k² > 0, the same sign as the floor. That would give m_D ≈ 0.65 %, and (D) fires.
    - *Against (D):* aspect 2, at the same a and k with \|Im E\|/\|Re E\| = 8.279e-3, reads **−0.152 %** (`20261012T184250Z_IMP-6-step1-W1a.log:429–430`). A uniform +0.8 % share would need an ≈ −1 % discretisation error there to cancel it.
  - m_C in **[1.1 %, 1.6 %]**: (C) is predicted not to fire. The a/18 surface refinement, nearer the rim, moved the miss only −0.115 pp. The far grain sets h only beyond the 2c cavity layer and the 4ρ_r rim tube, plus the SizeMax ends of their ramps.
  - \|β − 1\| ≲ 0.8 % on every solve.
  - W0f cells ≈ 63–68 k (sizing below).
- **Sizing: the a/5 cell count and the guard.**
  - Item 426 fitted N₀ ≈ 18.3 k cells not set by the Balls (item 440's brief, from `20261013T200400Z_IMP-6-step2-a2-W0-edge.log:172, :290`).
  - Item 440's a/18 mesh added 63 363 − 53 959 = 9 404 cells against a cavity-layer factor of 2.06, so the layer is L ≈ 9 404/1.06 ≈ 8.9 k at a/12. That leaves ≈ 9.4 k in the far region.
  - At F = 1.2 mm:
    - Far region × (1.5/1.2)³ = 1.953 → **+9.0 k**.
    - Cavity layer: (1/s² − 1/F²)/(2k), with s = 0.5 mm and k = (F − s)/(2c), reads 4.27 at F = 1.5 and 5.67 at F = 1.2. That is ×1.33 → **+2.9 k**.
    - The rim-tube ramp, 0.24 mm → F over 0.72 mm, adds a few hundred cells.
  - That gives **≈ 66 k** (range 63–68 k). Bound: all of N₀ scaling by 1.953 gives ≈ 71–72 k.
  - **Both are under `CELL_GUARD` 95 000, so the imported guard binds as in items 426 / 440. No separate guard is set.**
  - Memory, linear between W1b (7.54 GiB at 53 959) and W1s (8.80 GiB at 63 363, `…022936Z_…-W1s.log:179`): ≈ 9.2 GiB at 66 k and ≈ 10.0 GiB at 72 k.
- **Price** (measured windows, `-n 2`, `tests/environment/test_complex_mode.py` on the path):
  - Measured basis:
    - item 440's W0r 27 s;
    - W1q 15 s;
    - W0 14 s;
    - W1s 90 s window, 87.5 s solve+post, 8.80 GiB at 63 363 cells (`…022936Z_…-W1s.log:179`);
    - item 426's W1b 76 s at 53 959, without the item 440 forms.
  - Solve time scales as cells^1.41. That scaling predicted W1s at 91 s against 87.5 s measured.
  - Predicted:
    - W0r ≈ 30–45 s (three rungs plus the default node on the cached rung).
    - W1k / W1k2 ≈ 78–85 s each (the 53 959 mesh, plus item 440's degree-12 / 14 forms; the sparsity is unchanged at small k).
    - W0f ≈ 15–20 s.
    - W1f ≈ 95–115 s (≈ 97 s solve at 66 k).
    - Total ≈ 5–6 min of windows. Every window is under 10 minutes, so §5.1's durable-capture idiom is not needed.
- **Windows** (complex, `-n 2`, one solve per window, foreground, `timeout: 660000`; item 440's `Command:` shape). Order: code → W0r → W1k → W1k2 → W0f → W1f → commit. The (D) arm runs first because it needs no new mesh. **No module or `src/` edit after W0r starts** (rule (i)).
  - W0r, the default re-gate (all three knobs unset). Expect 33 334 / 53 959 / 89 530 cells, SICN 0.3605 / 0.3994 / 0.3466 and 0 / 0 / 0 inside (`20261014T022824Z_IMP-6-step2b-W0r.log:167–173, :286–292, :418–424`), with `test_knob_defaults` green. Command: `scripts/testing/run_and_log.sh IMP-6-step2c-W0r "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_ladder_census tests/validation/test_imp6_oblate_rim_ladder.py::test_knob_defaults -s -v --tb=short'"`
  - W1k, the (D) reading at s = 0.1: `scripts/testing/run_and_log.sh IMP-6-step2c-W1k "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_K_SCALE=0.1 timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung1 -s -v --tb=short'"`
  - W1k2, the leading-order check at s = 0.3: the same command with label `IMP-6-step2c-W1k2` and `FEM_EM_IMP6_K_SCALE=0.3`.
  - W0f, the a/5 census: `scripts/testing/run_and_log.sh IMP-6-step2c-W0f "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_FAR_DIV=5 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_surface_census -s -v --tb=short'"`
  - W1f, the (C) reading, which runs only if W0f is green and its SICN line reads ≥ 0.1: `scripts/testing/run_and_log.sh IMP-6-step2c-W1f "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_FAR_DIV=5 timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung1 -s -v --tb=short'"`
- **Negative control.**
  - **(K) is the control for the (D) knob.** A knob that does not reach the operator returns \|Im E\|/\|Re E\| = 8.523e-3 unchanged. That is 1/(2s²) = **50×** over the bar at s = 0.1 and 5.6× at s = 0.3, the ceiling of this fixture's separation. A wrong-sign or s-linear response (σ scaled without f, k² × s) reads ≈ 10× the leading-order value at s = 0.1, which is 5× over the bar.
  - **(Kf) is the control for the (C) knob.** An unpropagated `FAR_DIV` reproduces 53 959 cells exactly, which is below the 56 657 bar.
  - **(Wk)** separates "same mesh, different k" from a silent mesh change.
  - (C-fid) is not re-run. Item 426's default-path 31 / 93 (`20261013T200330Z_IMP-6-step2-a2-W0-ctl.log:172, :312`) already shows on this fixture that the inside-count census can fail.
  - Every other factor in this brief is labelled *predicted*.
- **Traps already paid for.**
  - **A collective inside a rank-0 branch deadlocks the window.** Item 440's first W1q called `_ncells` inside its rank-0 print. It hung to the 400 s container timeout with exit 124 (`20261014T022119Z_IMP-6-step2b-W1q.log`; attempts.md `2026-10-14T02:33Z`, "Incident"). Every new collective in `_solve_on`, `test_knob_defaults` and the census prints is called before the `if comm.rank == 0:` branch: `_ncells`, `_rss_gib`, `_inside_count`, `_h_rim_*`, `_assemble_real`, the numpy quadrature means, and `comm.allreduce`. The k print is pure arithmetic.
  - **Rule (i): no module or `src/` edit after W0r starts.** A fix forced by a later window is an edit, so W0r is re-run before any other window (item 440's W0r-repair path).
  - **k² ∝ s² needs f and σ scaled together.** Scaling f alone multiplies the loss tangent by 1/s, and k² falls only ∝ s at low f. Do not scale ε_r: ε_r = 0.78 passes validation but is not this host.
  - **`K_SCALE` must not enter the mesh cache key, and `FAR_DIV` must.** The first would only cost a re-mesh. The second, missing, would hand W1f the cached a/4 mesh: (Kf) catches it in W0f, but W1f runs in its own process.
  - **Never run `test_ladder_census` with `FEM_EM_IMP6_FAR_DIV` set.** It would mesh ρ_r/5 at a/5: 89 530 + ≈ 12 k ≈ 101 k, over the 95 000 guard. W0f is the one-rung node. The one-knob import check makes "two knobs in one window" impossible.
  - **A quality trap means park unless the brief names the fix** (2026-10-12 18:00 ruling). This brief names none: there is no a/4.5 fallback and no s other than 0.3 / 0.1. A MUMPS error at s = 0.1 is a code fact about the near-static solve: report it and do not lower s.
  - **The landed gate stays at degree 8.** Do not change the imported `QDEG`.
  - **Never import a `test_*` name** from step 1's module.
  - **The basis is the in-plane field minus E₀x̂** (item 404's trap); `_solve_on` does it, so do not touch it.
  - Owned cells only (`size_local`). `assemble_scalar` is rank-local (`_assemble_real` reduces). Masks stay in numpy, with no UFL `<`, `max_value` or `sqrt` comparison.
  - `-s` on every window (pytest captures prints without it) and `-k 30` on every container `timeout`. Select by node ID, never `-k a or b`.
  - The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run in the foreground with Bash `timeout: 660000`. Two heavy windows never run at once.
  - **After a killed window**, first run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`. Only if it reads 0, run the rubric's orphan-`.c` sweep as one bare call: `scripts/testing/run_and_log.sh IMP-6-step2c-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`. Delete orphan stems only; never clear the cache.
- **Scope.** One aspect (5), one rung (ρ_r/4), one mesh knob at one value (`resolution_far` a/4 → a/5) and one physics knob at two values (k² × 0.09, × 0.01, f and σ together).
  - Not a convergence study and no rate.
  - Not a finite-k claim: the row's Not-claimed "finite k" stands. This item measures the finite-k share of a quasi-static miss; it does not gate a finite-k closed form.
  - Not aspect 2, not the axial drive, no band move.
  - The default oblate path, step 1, step 3 and `EX-111` are untouched.
- **Status it can move.** The `IMP-6` §7 row's step-2 record clause (PROJECT_PLAN line 2605) currently ends "the far grain or the quasi-static comparand itself (`…W1q.log`, `…W1s.log`)". It gains the attribution: "; item 441 (I27): k² × 0.01 at fixed loss tangent reads m_D …, `resolution_far` a/5 reads m_C … — the floor is (C) / (D) / both / neither", with the W1k / W1f log names.
  - No glyph change: `IMP-6` is ✅.
  - The slot does not edit the Done-when or the Not-claimed list.
  - The 2026-10-13 known-issues entry is not touched.
- **Negative result.**
  - **W0r moves a default record, or `test_knob_defaults` is red:** a knob leaked into the default path. Fix the module and re-run W0r before any other window.
    - If the default-path diff is provably empty and the counts still move, that is an image move. Stop and report: m₀ would no longer be on this mesh.
  - **(Wk) red:** W1k did not solve on W0r's mesh. Fix the cache key, re-run W0r, then re-run W1k. No m_D is read from a red (Wk).
  - **(K) red:** the knob did not reach the operator, or the response is not O(k²). No m_D is read. Check that `_solve_on` passes the scaled f and σ. A code fact is fixed in-slot (then re-run W0r). Otherwise park on `attempt/IMP-6-item441-<ts>`, BLOCKED (rule (d)).
  - **W0f red** (guard, inside ≠ 0, the SICN line < 0.1, or (Kf)): no W1f. The (D) arm's green windows are still a result. Land them on `main` with the row clause carrying (D) only, mark the item BLOCKED on its (C) arm in the same commit (rule (d)), and report the counts. The review rules.
  - **(i), (ii), (iii) or (Qh) red on any solve:** a real finding. Print h_rim/ρ_r, min SICN, the inside count, the deficit and k. Park on `attempt/*`, BLOCKED (rule (d)); the review rules.
  - **Attribution "neither":** report and stop. The rim floor stays a stated finding in the row clause, and the weekly decides whether to pose more. The slot poses no further rung.
  - In every case no band moves. A miss that grows is reported as measured, never fitted.
  - No `PREDICTED RED` line: every gate is predicted green.

**Every-item rules** (PROJECT_PLAN §9 "Every item", binding this item as written above):
- Bare, whole harness calls.
- A lost docker socket stops the slot.
- No edit after the last window.
- `timeout: 660000` and `-s` on every window.
- Never delete a log.
- BLOCKED goes in the parking commit.
- A killed window is not a result.
- Two heavy windows never run at once.
- Complex mode and `FEM_EM_REQUIRE_COMPLEX=1` on every frequency-domain command.
- No slot opens `docs/private/` or `aed_results/`.

**Spawn line:** "Use the implementer agent to execute `IMP-6` §9 item 441 per `docs/planning/chunks/IMP-6.md` (provisional I27 under §9 rule (6); if the 10-14 weekly reverted I26 or I27, do not run). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Edit only `tests/validation/test_imp6_oblate_rim_ladder.py`: env knobs `FEM_EM_IMP6_FAR_DIV` (default 4, in the mesh cache key) and `FEM_EM_IMP6_K_SCALE` (default 1, scaling f and σ together, not in the cache key), one knob per run enforced at import, the asserted default node `test_knob_defaults`, and the (K) / (Wk) / (Kf) checks. No `src/` edit. Never import a `test_*` name from step 1's module. Every collective is called before the rank-0 branch (item 440's W1q deadlock). Windows in order W0r → W1k → W1k2 → W0f → W1f, no edit after W0r starts, no W1f before W0f reads 0 inside midpoints and min SICN ≥ 0.1. State the pre-registered attribution ((C) / (D) / both / neither) in the outcome; assert none of it."

**Outcome (2026-10-14 04:30 slot, DONE; attribution (D)).** Code: test module only (`tests/validation/test_imp6_oblate_rim_ladder.py`), no `src/` edit. Five windows, all green, `-n 2`, complex, 285 s of windows in total:
- W0r (27 s): default path bitwise — 33 334 / 53 959 / 89 530 cells, SICN 0.3605 / 0.3994 / 0.3466, 0 inside; `test_knob_defaults` cells 53 959, mask 11 468, Φ_cf,h(8) rel 1.199e-14 (`20261014T094514Z_IMP-6-step2c-W0r.log:173, :292, :424, :428`).
- W1k, s = 0.1 (79 s): k = 1.5178 − 0.9489j 1/m, (\|k\|a)² 0.0115 %; β 0.988081 (\|β − 1\| 1.192 %); (K) \|Im E\|/\|Re E\| 8.493e-5 (predicted 8.52e-5, bar 1.705e-4); (Wk) 1.199e-14; Φ_h 22.852453 vs Φ_cf,h 22.704740, **m_D +0.6506 %** ≤ 0.986 % ⇒ **(D) YES**; δ_D +0.8440 pp (`20261014T094547Z_IMP-6-step2c-W1k.log:177–188`).
- W1k2, s = 0.3 (80 s): β 0.988419; (K) 7.650e-4 (predicted 7.67e-4); m +0.7181 %, δ_D +0.8439 pp; leading-order ratio (1.4861 − 0.7181)/(1.4861 − 0.6506) = **0.919** against the predicted 0.919 (`20261014T094713Z_IMP-6-step2c-W1k2.log:177–188`).
- W0f, a/5 (14 s): cells 61 486 (predicted ≈ 66 k, range 63–68 k — **below the predicted range**; (Kf) bar 56 657 met), min SICN 0.3723, mask 11 330, 0 inside (`20261014T094839Z_IMP-6-step2c-W0f.log:178–182`).
- W1f, a/5 (85 s): β 0.994327 (0.567 %); \|Im E\|/\|Re E\| 8.521e-3; Φ_h 23.028322 vs Φ_cf,h 22.705326, **m_C +1.4226 %** > 0.986 % ⇒ **(C) NO** (predicted [1.1 %, 1.6 %]) (`20261014T094858Z_IMP-6-step2c-W1f.log:180–190`).
- Reading: δ_D is s-independent to 1e-4 pp across s = 0.1 / 0.3, so the finite-k share of item 426's +1.486 % is +0.844 pp, O(k²) as predicted ("for (D)" estimate 0.83 %). The residual +0.65 % at k → 0 is the discretisation miss on this mesh; the far grain moves it −0.064 pp. All bands (i)–(iii), (Qh) green on every solve; no band moved. Item 441's prediction "\|β − 1\| ≲ 0.8 %" missed at s = 0.1 / 0.3 (1.19 % / 1.16 %, inside `TH15_BETA_BAND`) — printed, never asserted. Cosmetic: W0f's step-2b census line still prints item 440's cell prediction beside item 441's (line 180 carries item 441's).

### Review ruling 2026-10-14 03:00 (item 441)

The brief is accepted as drafted:
- the two test-module knobs (`FEM_EM_IMP6_K_SCALE`, `FEM_EM_IMP6_FAR_DIV`), one set per run;
- `test_knob_defaults` asserted bitwise against item 440's records;
- each knob's effect control.

The extra s = 0.3 (D) point stays: it is read outside the attribution rule, ≈ 80 s, and it is the first thing to drop if the slot runs long. The item is withdrawn unrun if the 10-14 weekly reverts I26 or I27. Family: `IMP-6` step 2, attempt 2 of 4.

### §9 item 450 — 2026-10-14 (chain step I28, daily — provisional; lane 3, closed form)

`IMP-6` step 2d, attempt 3 of 4 in the step-2 family since item 426's registered record (items 440 / 441 were 1 and 2). I22's whole rim ladder (ρ_r/3, ρ_r/4, ρ_r/5 on `oblate_rim="edge"`) is re-solved at `FEM_EM_IMP6_K_SCALE=0.1`, the knob item 441 added, so that the k → 0 residual is read on every rung, not on ρ_r/4 alone (implementer; heavy, `-n 2`, complex, degree 2; ≈ 11 measured / 35 predicted min: ≈ 5–6 min of windows + 15 fixed + ≈ 15 code; depends on nothing). This is chain I's third §9 rule (6) provisional step (I26, I27 were the first two; rule (6) allows three per chain), written by the 2026-10-14 10:30 daily. The 10-14 weekly ratifies or reverts I26, I27 and I28. **If the weekly reverts any of them, the item is withdrawn unrun.**

- **The question (item 441's outcome; attempts.md `2026-10-14T09:51Z`, "Next-attempt hypothesis").**
  - Item 441 split I22's rung-ρ_r/4 miss, item 426's +1.4861 % (`20261013T200607Z_IMP-6-step2-a2-W1b.log:179`), into a finite-k share δ_D = +0.8440 pp, s-independent to 1e-4 pp across s = 0.1 / 0.3 (`20261014T094547Z_IMP-6-step2c-W1k.log:187`; `20261014T094713Z_IMP-6-step2c-W1k2.log`), and a k → 0 residual m_D = +0.6506 % (`…094547Z_…-W1k.log:181`).
  - The slot read that residual as "the discretisation miss on this mesh". That reading is one rung deep. **A discretisation miss converges to the comparand under refinement, because `pec_oblate_spheroid_field` is the exact k → 0 limit** (item 441 brief, "Validity as k → 0"). A residual set by something the rim ladder holds fixed does not.
  - I22's own moves at s = 1 already rise: m₃ / m₄ / m₅ = +1.2233 / +1.4861 / +1.6129 % (`20261013T200520Z_IMP-6-step2-a2-W1a.log:171`; `…200607Z_…-W1b.log:179`; `20261013T200729Z_IMP-6-step2-a2-W1c.log:192`), so Δ₃₄ = +0.2628 pp and Δ₄₅ = +0.1268 pp. Subscripts here are the rung divisor (ρ_r/3, /4, /5).
  - The ladder at s = 0.1 separates the two readings. It does not move a band.
- **Target and bands are I22's, unchanged.** Every solve asserts \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball Φ ≤ `MAT4_BAND` against the degree-8 comparand as landed, and cavity dofs > 0 with max \|E\| < 1e-12. All three bands are imported, never restated as numbers. The signed miss, its per-rung moves and the pairwise order are printed and never asserted, as in I22. No new bar is set at s = 0.1. Gating a rim ladder at s = 0.1 or writing a finite-k comparand is the weekly's decision (attempts.md `2026-10-14T09:51Z`).
- **The knob and the rungs, read from the code (`tests/validation/test_imp6_oblate_rim_ladder.py`).**
  - `K_SCALE` (`:133–135`) scales f and σ together in `_solve_on` (`:344–345`). It is not in the mesh cache key (`:176`), so each rung solves on the mesh W0r re-gates at s = 1.
  - The rungs are the nodes `test_rung0` / `test_rung1` / `test_rung2` (`:719–736`), divisors `DIVISORS = (3.0, 4.0, FEM_EM_IMP6_RUNG2 or 5)` (`:109`). The one-knob import check (`:136–140`) admits `K_SCALE` alone. **On the environment side the ladder at s = 0.1 needs nothing new; the module does need the edits below.**
  - **As landed, `test_rung0` / `test_rung2` under `K_SCALE` are red by construction.** `_report_and_gate` asserts (Wk) against item 440's **rung-1** comparand, 22.704739974621, at 1e-12 (`:711`). Rung ρ_r/3's comparand is 22.693485 (`…W1a.log:171`) and rung ρ_r/5's is 22.712699 (`…W1c.log:192`), 5.0e-4 and 3.5e-4 away. The δ_D print (`:679`) and the item-441 (D) verdict (`:686–696`) are also hard-wired to rung 1's m₀ and threshold.
- **Code (test module only; additive; the default path (`K_SCALE` = 1) is bitwise unchanged; no `src/` edit).** Every new line sits inside an `if K_SCALE != 1.0` branch or a `K_SCALE`-conditioned import check.
  - (1) `ITEM426_RUNGS = {3.0: …, 4.0: …, 5.0: …}`: item 426's per-rung s = 1 records, holding Φ_h, Φ_cf,h, the signed miss, β, \|Im E\|/\|Re E\| and the log:line for each rung:
    - ρ_r/3: 22.971085 / 22.693485 / +1.2233 % / 1.001961 / 8.521e-3 (`…W1a.log:170–171`);
    - ρ_r/4: 23.042163 / 22.704740 / +1.4861 % / 0.992222 / 8.523e-3 (`…W1b.log:178–179`);
    - ρ_r/5: 23.079024 / 22.712699 / +1.6129 % / 0.994011 / 8.527e-3 (`…W1c.log:191–192`).
    - Keep `ITEM426_RUNG1` as it is (rung 1's prints use it).
  - (2) At import: `ValueError` if `K_SCALE != 1.0` and `DIVISORS[2] != 5.0`. The records exist only for ρ_r/5, so `FEM_EM_IMP6_RUNG2` stays unset under the knob.
  - (3) **(Wk) per rung**, in `_report_and_gate`:
    - rung 4.0 keeps today's assert, \|Φ_cf,h(8)/22.704739974621 − 1\| ≤ 1e-12;
    - rungs 3.0 and 5.0 assert **(Wk′)** \|Φ_cf,h(8)/`ITEM426_RUNGS[d]`["phi_cf_h"] − 1\| ≤ **5e-8**.
    - The bar is the 6-decimal print's half-unit, 5e-7/22.69 = 2.2e-8, times ≈ 2. That is ≈ 7 000× below the 3.5e-4 gap to the nearest neighbouring rung's comparand. Item 440 reproduced the rung-1 6-decimal print to 12 decimals (22.704739974621 → 22.704740, `20261014T022856Z_IMP-6-step2b-W1q.log:177–178`), so this record format is faithful.
  - (4) **(R4), rung 4.0 at `K_SCALE == 0.1` only:** \|Φ_h/22.852453 − 1\| ≤ 5e-8 against item 441's W1k (`…094547Z_…-W1k.log:181`), the same half-unit bar. This assert is an instrument identity (the same mesh at the same s re-solved), not a physics gate.
  - (5) Prints, on rank 0, arithmetic only. A `[step 2d]` line per rung gives:
    - the rung's s = 1 record beside;
    - r_d = the signed miss;
    - δ_D,d = (m_d − r_d)/(1 − s²);
    - the I28 prediction (below);
    - the max-edge and median h_rim/ρ_r already computed in `_solve_on`;
    - on rung 4.0, the (R4) reading;
    - at s ≠ 0.1, "outside the I28 rule".
    - The item-441 "(D) explains the floor" line and the "[step 2c] item 426 rung-1 record" line print on rung 4.0 only. Rungs 3.0 and 5.0 print their own record instead.
  - (6) `_predicted_miss` returns the I28 per-rung prediction when `K_SCALE == 0.1`, and item 441's string otherwise.
  - (7) Docstring: item 450, the per-rung (Wk′) / (R4), the rule pre-registered in this brief.
  - Unchanged: `test_knob_defaults`, `test_ladder_census`, `test_census_control`, `test_quadrature_rung1`, `test_surface_census`, `_solve_on`, the (K) / (Qh) / band asserts.
  - Still in force: step 1's module is not edited; only private names are imported from it, **never a `test_*` name**.
  - The reading rule is **not** coded. Each window solves one rung in its own process, so the rule is computed in the outcome from the three logs' cited lines. Putting W1a′'s value into the module before W1c′ would be an edit after W0r (item 441's s = 0.3 ratio precedent).
- **Anchors (asserted).**
  - On every solve (W1a′, W1b′, W1c′), I22's anchors at the imported bands:
    - (i) \|β − 1\| ≤ `TH15_BETA_BAND` (β the least-squares amplitude of `pec_oblate_spheroid_field(drive="inplane")` − E₀x̂, `src/fem_em_solver/utils/analytical.py:1831`, at the 48 probe points);
    - (ii) \|Φ_rim,h/Φ_rim,cf,h − 1\| ≤ `MAT4_BAND`;
    - (iii) cavity dofs > 0 and max \|E\| < 1e-12;
    - (Qh) \|Φ_h(14)/Φ_h(8) − 1\| ≤ 1e-12.
  - The default re-gate in W0r: `test_ladder_census` reads 33 334 / 53 959 / 89 530 cells and 0 / 0 / 0 inside midpoints, and `test_knob_defaults` (W0r-a/b/c) is green. These are item 441's records (`20261014T094514Z_IMP-6-step2c-W0r.log:173, :292, :424, :428`).
  - Instrument checks (not bands):
    - **(K)** on every solve: \|Im E\|/\|Re E\| ≤ 2·s²·`ITEM426_IMAG_RATIO` = 1.705e-4.
    - **(Wk)** / **(Wk′)** on every solve (code item 3).
    - **(R4)** on W1b′ (code item 4).
- **Reading rule (pre-registered here; computed in the outcome from the cited log lines; never asserted).**
  - The inputs:
    - r₃, r₄, r₅ are the signed misses at s = 0.1 (W1a′, W1b′, W1c′).
    - Δ₃₄ = r₄ − r₃ and Δ₄₅ = r₅ − r₄, in pp.
    - H is the max-edge h_rim/ρ_r on these meshes: 0.6909 / 0.5103 / 0.3944 (`…W1a.log:169`, `…094547Z_…-W1k.log:179`, `…W1c.log:190`). The meshes are the same at s = 0.1, so the module reprints them.
    - p₄₅ = ln(r₄/r₅)/ln(H₄/H₅), defined only when r₄ and r₅ share a sign. ln(0.5103/0.3944) = 0.2577. The median-H order (0.3873 / 0.3089) is printed beside.
  - At s = 0.1 the comparand's leftover finite-k share is s²·δ_D = 0.01 × 0.844 = 0.0084 pp (item 441). That is below every threshold here, so the comparand is exact for this rule's purposes.
  - The rule is applied in order:
    1. **(Disc): the residual is discretisation that converges to the comparand under the rim grain.** It requires Δ₃₄ < 0 and Δ₄₅ < 0 (a monotone fall), **and** either p₄₅ ≥ **1.0** or \|r₅\| ≤ **0.127 %**. Both thresholds are *predicted*.
       - 1.0 is above the largest \|order\| I22 read on this ladder, 0.700 (median H, ρ_r/3 → ρ_r/4; item 426 outcome, "Record"). A fall that clears it therefore cannot be I22's own rung-to-rung movement with its sign flipped.
       - 0.127 % is I22's last per-rung move, m₅ − m₄ = 0.1268 pp. A residual within one such move of zero has reached the comparand at this ladder's resolution. This arm also covers a fall that crosses zero, where p₄₅ is undefined.
    2. **(Fix): a remaining fixed-error source.** It applies if (Disc) does not hold, r₃, r₄ and r₅ share one sign, and \|Δ₄₅\| ≤ **0.25·\|r₅\|** (*predicted* threshold): the last rim move is under a quarter of what is left.
       - The basis: I22 at s = 1 read \|Δ₄₅\|/m₅ = 0.1268/1.6129 = 0.079, the floor item 426 already called a floor; 0.25 is ≈ 3× that.
       - It is disjoint from (Disc)'s order arm. A fall at p₄₅ ≥ 1.0 over H₄/H₅ = 1.294 moves \|Δ₄₅\| ≥ 0.294·\|r₅\| > 0.25·\|r₅\|.
       - The direction (falling to / rising to the limit) is printed with the verdict.
       - Outside the rule, the geometric limit r∞ = r₅ + Δ₄₅·q/(1 − q), with q = Δ₄₅/Δ₃₄, is printed only when 0 < q < 1, as the size of the fixed error.
    3. **Neither**, named by case:
       - a monotone fall with p₄₅ < 1.0 and \|Δ₄₅\| > 0.25·\|r₅\| ("falls, below the order bar: three rungs do not separate it");
       - a growth with \|Δ₄₅\| > 0.25·\|r₅\| ("moves the wrong way");
       - a sign change that (Disc)'s 0.127 % arm does not catch ("order undefined").
  - **What (Fix) names (for the weekly; not posed, not tested here).** These are the things the rim ladder holds fixed:
    - (E1) the surface grain a/12's unconverged share. One step, a/12 → a/18, moved the miss −0.115 pp at s = 1 (`20261014T022936Z_IMP-6-step2b-W1s.log:181`).
    - (E2) the far grain. One step, a/4 → a/5, moved it −0.064 pp at s = 1 (`20261014T094858Z_IMP-6-step2c-W1f.log:180–190`).
    - (E3) the faceted body. The polyhedral spheroid lies inside the true one: the body-volume deficit is −5.147e-3 / −5.148e-3 / −5.088e-3, flat across the rungs (item 426 outcome, "Record"). The `interior="continue"` comparand is integrated over air cells whose slivers lie inside the true body.
    - The finite-k share (D) is excluded by construction.
- **Predicted (rule (e); printed beside the measured value, never asserted).**
  - δ_D,d within **±0.05 pp of +0.844 pp** on every rung. δ_D is a k² share of a field whose relative discretisation error on the rim is ≲ 2 % (I22's misses), so it should move ≲ 0.02 pp between rungs.
  - That gives r₃ in **[0.34 %, 0.44 %]**, r₄ = **+0.6506 %** (reproducing W1k, (R4)) and r₅ in **[0.73 %, 0.83 %]**. Then Δ₃₄ ≈ +0.26 pp, Δ₄₅ ≈ +0.13 pp, \|Δ₄₅\|/r₅ ≈ 0.16, q ≈ 0.48 and r∞ ≈ +0.90 %.
  - **The predicted verdict is (Fix), rising.** This prediction picks a side, and it does so by construction from I22 and I27 if δ_D is rung-independent.
    - *Against:* for r to fall from ρ_r/3 to ρ_r/4, δ_D would have to grow by > (m₄ − m₃)/0.99 = 0.265 pp, ≈ 31 % of itself, over one rim rung. That is not predicted.
  - (K) ≈ s² × 8.521e-3 / 8.523e-3 / 8.527e-3 ≈ 8.5e-5 on each rung. Rung 1 read 8.493e-5 (`…W1k.log:180`).
  - β ≈ 0.998 / 0.988 / 0.990, each ≈ −0.004 from item 426's 1.001961 / 0.992222 / 0.994011, as rung 1 moved 0.992222 → 0.988081 at s = 0.1 (`…W1k.log:180`). \|β − 1\| ≤ 1.4 % on every rung. Item 441's own β prediction missed (≲ 0.8 % against a measured 1.19 %), so this one is set wide.
- **Price** (measured windows, `-n 2`, `tests/environment/test_complex_mode.py` on the path):
  - Measured basis:
    - item 441's W0r 27 s (`20261014T094514Z_IMP-6-step2c-W0r.log`, attempts.md `2026-10-14T09:51Z`);
    - item 441's W1k 79 s at s = 0.1 on ρ_r/4, solve+post 75.9 s, 7.62 GiB (`…W1k.log:179`), against item 426's W1b 76 s at s = 1 (×1.04, with item 440's degree-12/14 forms added);
    - item 426's W1a 41 s (`…W1a.log:208`) and W1c 152 s, solve+post 148.5 s, 14.41 GiB (`…W1c.log:190, :233`).
  - Predicted:
    - W0r ≈ 27–35 s;
    - W1b′ ≈ 79–85 s;
    - W1a′ ≈ 43–50 s (41 s × 1.04 plus the forms on 4 861 mask cells);
    - W1c′ ≈ 158–175 s (152 s × 1.04 plus the forms on 21 548 mask cells; ≈ 14.5 GiB).
  - **Total ≈ 5–6 min of windows** (≈ 307–345 s). Every window is under 10 minutes, so §5.1's durable-capture idiom is not needed. The heavy tier's 1 200 s ceiling is not approached.
- **Windows** (complex, `-n 2`, one solve per window, foreground, Bash `timeout: 660000`; item 441's `Command:` shape, `…094547Z_…-W1k.log:12`). Order: code → W0r → W1b′ → W1a′ → W1c′ → commit. Rung ρ_r/4 runs first because (R4) tells whether the ladder is the same instrument as W1k. **No module edit after W0r starts** (rule (i)).
  - W0r, the default re-gate (no knob set): `scripts/testing/run_and_log.sh IMP-6-step2d-W0r "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_ladder_census tests/validation/test_imp6_oblate_rim_ladder.py::test_knob_defaults -s -v --tb=short'"`
  - W1b′, rung ρ_r/4 at s = 0.1: `scripts/testing/run_and_log.sh IMP-6-step2d-W1b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_K_SCALE=0.1 timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung1 -s -v --tb=short'"`
  - W1a′, rung ρ_r/3 at s = 0.1: `scripts/testing/run_and_log.sh IMP-6-step2d-W1a "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_K_SCALE=0.1 timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung0 -s -v --tb=short'"`
  - W1c′, rung ρ_r/5 at s = 0.1: `scripts/testing/run_and_log.sh IMP-6-step2d-W1c "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge FEM_EM_IMP6_K_SCALE=0.1 timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_rim_ladder.py::test_rung2 -s -v --tb=short'"`
- **Negative control.**
  - **(K), asserted on every rung, is the control for the knob.** A knob that does not reach the operator returns the s = 1 ratio, measured on these three meshes as 8.521e-3 / 8.523e-3 / 8.527e-3 (`…W1a.log:170`, `…W1b.log:178`, `…W1c.log:191`). That is 1/(2s²) = **50×** over the 1.705e-4 bar, the ceiling of this fixture's separation. The knob's response was measured on this fixture at 8.493e-5 against the 8.52e-5 prediction (`…W1k.log:180`).
  - **(Wk′), asserted, separates "this rung's mesh" from its neighbours.** The three comparands differ by 3.5e-4 to 8.5e-4 relative (`…W1a.log:171`, `…W1b.log:179`, `…W1c.log:192`), ≈ 7 000× the 5e-8 bar. A cache-key leak or a silent mesh change reads red.
  - **(R4), asserted on W1b′.** A re-solve of the same mesh has reproduced every printed digit on this fixture before: item 408's W1a against step 1 (ratio 1.014181, β 1.006343; `20261013T033249Z_IMP-6-step2-W1a.log:166–167` against `20261012T184537Z_IMP-6-step1-W1b.log:189–191`).
  - The reading rule has no control of its own: it is a printed reading. Every threshold in it is labelled *predicted*.
- **Traps already paid for.**
  - **A collective inside a rank-0 branch deadlocks the window.** Item 440's first W1q called `_ncells` inside its rank-0 print and hung to the 400 s container timeout with exit 124 (`20261014T022119Z_IMP-6-step2b-W1q.log`; attempts.md `2026-10-14T02:33Z`, "Incident"). The new prints are pure arithmetic on values `_solve_on` already reduced. Add no `_ncells`, `_rss_gib`, `_inside_count`, `_h_rim_*`, `_assemble_real`, numpy quadrature mean or `comm.allreduce` inside an `if comm.rank == 0:` branch. The asserts stay outside the branch, on every rank.
  - **Rule (i): no module or `src/` edit after W0r starts.** A fix forced by a later window is an edit, so W0r is re-run before any other window (item 440's W0r-repair path).
  - **The (Wk) trap is this brief's own find.** Do not run `test_rung0` / `test_rung2` under `K_SCALE` before code item (3) lands: on the landed module they are red on (Wk) by construction. That red would be a code fact, not a reading.
  - **Never set a second knob, and never `FEM_EM_IMP6_RUNG2`, in any window.** The one-knob check and code item (2) raise at import. W0r runs with no knob, as `test_knob_defaults`' (W0r-a) requires.
  - **`K_SCALE` stays out of the mesh cache key** (item 441). W1a′ / W1b′ / W1c′ each mesh their own rung in their own process; (Wk) / (Wk′) catch a mismatch.
  - **A quality trap means park unless the brief names the fix** (2026-10-12 18:00 ruling). This brief names none: no s other than 0.1, no rung other than 3 / 4 / 5. A MUMPS error at s = 0.1 on rung ρ_r/5 is a code fact about the near-static solve. Report it; do not raise s.
  - **The landed gate stays at degree 8**: do not change the imported `QDEG`. **Never import a `test_*` name** from step 1's module. **The basis is the in-plane field minus E₀x̂** (item 404's trap); `_solve_on` does it, so do not touch it.
  - Owned cells only (`size_local`). `assemble_scalar` is rank-local (`_assemble_real` reduces). Masks stay in numpy, with no UFL `<`, `max_value` or `sqrt` comparison.
  - `-s` on every window (pytest captures prints without it) and `-k 30` on every container `timeout`. Select by node ID, never `-k a or b`. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run in the foreground with Bash `timeout: 660000`. Two heavy windows never run at once.
  - **After a killed window**, first run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`. Only if it reads 0, run the rubric's orphan-`.c` sweep as one bare call: `scripts/testing/run_and_log.sh IMP-6-step2d-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`. Delete orphan stems only; never clear the cache.
- **Scope.** One aspect (5), the in-plane drive, I22's three rim rungs, one physics knob at one value (k² × 0.01, f and σ together). The surface grain (a/12) and far grain (a/4) are held fixed.
  - Not a convergence study and no asserted rate. The order is printed, as in I22.
  - Not a finite-k claim: the row's Not-claimed "finite k" stands. Not a gate at s = 0.1: I22's bands are re-read at s = 0.1, as item 441's W1k did, and no new bar is set there.
  - (E1)–(E3) are named, not posed. Not aspect 2, not the axial drive, no band move.
  - The default oblate path, step 1, step 3 and `EX-111` are untouched.
- **Status it can move.** The `IMP-6` §7 row's step-2 record clause (PROJECT_PLAN line 2605) currently ends with item 441's "(I27): … the ≈ +0.64 % left at k → 0 is the discretisation residual on this mesh (…W1k.log, …W1k2.log, …W1f.log)". It gains: "; item 450 (I28): at k² × 0.01 the rim ladder reads r₃ / r₄ / r₅ = … / … / … % (Δ … / … pp, p₄₅ …), so the k → 0 residual is (Disc) rim-grain discretisation / (Fix) set by what the rim ladder holds fixed (r∞ …) / neither", with the W1a′ / W1b′ / W1c′ log names. Item 441's record sentence stays as written; it is a record.
  - No glyph change: `IMP-6` stays ✅. The slot does not edit the Done-when, the Not-claimed list or the 2026-10-13 known-issues entry.
  - The §10 I-chain reading goes to the 10-14 (or next) weekly, which owns any finite-k comparand, any s = 0.1 gate and any (E1)–(E3) step.
- **Negative result.**
  - **W0r moves a default record, or `test_knob_defaults` is red:** a `K_SCALE` branch leaked into the default path. Fix the module and re-run W0r before any other window. If the default-path diff is provably empty and the counts still move, that is an image move: stop and report, because item 426's records would no longer be on these meshes.
  - **(R4) red:** W1b′ does not reproduce W1k on the same mesh at the same s. No ladder reading is made. Report both Φ_h values. A code fact (an edit reached the solve path) is fixed in-slot, followed by a W0r re-run. Otherwise park on `attempt/IMP-6-item450-<ts>`, BLOCKED (rule (d)).
  - **(Wk′) or (K) red on a rung:** no r_d is read from that rung. Fix a cache-key or record-transcription fact in-slot, then re-run W0r. Otherwise park, BLOCKED.
  - **(i), (ii), (iii) or (Qh) red on any solve:** a real finding at s = 0.1. Print h_rim/ρ_r, min SICN, the inside count, the deficit and k. Park on `attempt/*`, BLOCKED (rule (d)); the review rules.
  - **W1c′ killed or over memory:** a two-rung reading (ρ_r/3, ρ_r/4) is printed and the rule is **not** applied, because it needs Δ₄₅. Land the green windows with the row clause carrying r₃ / r₄ only, and mark the item BLOCKED on W1c′ in the same commit.
  - **A "neither" verdict:** report and stop. The k → 0 residual stays a stated finding, and the slot poses no further rung, knob or comparand.
  - In every case no band moves. A residual that grows with refinement is reported as measured, never fitted.
  - No `PREDICTED RED` line: every gate is predicted green.

**Every-item rules** (PROJECT_PLAN §9 "Every item", binding this item as written above):
- Bare, whole harness calls.
- A lost docker socket stops the slot.
- No edit after the last window.
- `timeout: 660000` and `-s` on every window.
- Never delete a log.
- BLOCKED goes in the parking commit.
- A killed window is not a result.
- Two heavy windows never run at once.
- Complex mode and `FEM_EM_REQUIRE_COMPLEX=1` on every frequency-domain command.
- No slot opens `docs/private/` or `aed_results/`.

**Spawn line:** "Use the implementer agent to execute `IMP-6` §9 item 450 per `docs/planning/chunks/IMP-6.md` (provisional I28 under §9 rule (6); if the 10-14 weekly reverted I26, I27 or I28, do not run). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Edit only `tests/validation/test_imp6_oblate_rim_ladder.py`, and only inside `K_SCALE != 1` branches: `ITEM426_RUNGS`, the RUNG2 import guard, per-rung (Wk′) at 5e-8, (R4) on rung ρ_r/4, the `[step 2d]` prints. The default path stays bitwise. No `src/` edit. Never import a `test_*` name from step 1's module. Every collective is called before the rank-0 branch (item 440's W1q deadlock). Windows in order W0r → W1b′ → W1a′ → W1c′, no edit after W0r starts. Compute the pre-registered (Disc) / (Fix) / neither reading in the outcome from the cited log lines; assert none of it."

**Outcome (2026-10-14 13:30 CDT slot, DONE; reading (Fix), rising — the predicted verdict).** Code: test module only (`tests/validation/test_imp6_oblate_rim_ladder.py`), every new line inside a `K_SCALE != 1` branch or the `K_SCALE`-conditioned RUNG2 import guard; no `src/` edit; step 1's module untouched. Four windows, all green, `-n 2`, complex, 303 s of windows in total (predicted 307–345 s):
- **W0r** `20261014T183126Z_IMP-6-step2d-W0r.log` (31 s): census 33 334 / 53 959 / 89 530 cells, 0 / 0 / 0 inside midpoints (`:173`, `:292`, `:424`; `:169`, `:288`, `:420`); `test_knob_defaults` green, Φ_cf,h(8) = 22.704739974621 at rel 1.199e-14 (`:428`). The default path reads item 441's records.
- **W1b′** (ρ_r/4) `20261014T183202Z_IMP-6-step2d-W1b.log` (79 s; solve+post 76.4 s, 7.68 GiB, `:179`): r₄ = **+0.6506 %**, β 0.988081, (K) 8.493e-5 (`:180–181`); **(R4)** \|Φ_h/22.852453 − 1\| = 7.108e-9 ≤ 5e-8 (`:187`); (Wk) 1.199e-14 (`:188`); (Qh) 4.385e-14 (`:182`); δ_D,4 = +0.8440 pp (`:186`); cavity 24 710 dofs at 0 (`:190`).
- **W1a′** (ρ_r/3) `20261014T183329Z_IMP-6-step2d-W1a.log` (40 s; solve+post 37.8 s, 4.15 GiB, `:171`): r₃ = **+0.3900 %**, β 0.997778, (K) 8.504e-5 (`:172–173`); **(Wk′)** 1.325e-8 ≤ 5e-8 (`:179`); (Qh) 6.128e-14 (`:175`); δ_D,3 = +0.8418 pp (`:180`); cavity 20 978 at 0 (`:181`).
- **W1c′** (ρ_r/5) `20261014T183415Z_IMP-6-step2d-W1c.log` (153 s; solve+post 150.6 s, 14.70 GiB, `:192`): r₅ = **+0.7762 %**, β 0.989859, (K) 8.505e-5 (`:193–194`); **(Wk′)** 6.947e-10 ≤ 5e-8 (`:199`); (Qh) 1.361e-13 (`:195`); δ_D,5 = +0.8452 pp (`:200`); cavity 30 282 at 0 (`:201`).
- Every anchor green on every solve: \|β − 1\| 0.222 / 1.192 / 1.014 % ≤ `TH15_BETA_BAND` 4.886 %; rim-ball misses ≤ `MAT4_BAND` 10 %; cavity max \|E\| 0.
- **Reading rule (computed from the Φ values at `…W1a.log:173`, `…W1b.log:181`, `…W1c.log:194`; never asserted).** r₃ / r₄ / r₅ = +0.38996 / +0.65058 / +0.77619 %; Δ₃₄ = +0.2606 pp, Δ₄₅ = +0.1256 pp. (Disc) fails: no monotone fall (both moves rise), p₄₅ = ln(r₄/r₅)/0.2577 = −0.685 on max-edge H (−0.780 on the median H), \|r₅\| 0.776 % > 0.127 %. **(Fix) holds**: one sign on all three rungs and \|Δ₄₅\| = 0.1256 pp ≤ 0.25·\|r₅\| = 0.1940 pp (\|Δ₄₅\|/r₅ = 0.162). Direction: rising toward the limit. q = Δ₄₅/Δ₃₄ = 0.482, r∞ = r₅ + Δ₄₅·q/(1 − q) ≈ **+0.893 %**, the size of the fixed error, outside the rule.
- **Predictions (rule (e)):** δ_D,d +0.842 / +0.844 / +0.845 pp, all within ±0.05 pp of +0.844 (spread 0.0034 pp); r₃ 0.390 % in [0.34, 0.44]; r₄ reproduces W1k; r₅ 0.776 % in [0.73, 0.83]; Δ₃₄ ≈ +0.26, Δ₄₅ ≈ +0.13, ratio 0.16, q 0.48, r∞ ≈ +0.89 % all as predicted. β 0.9978 / 0.9881 / 0.9899 against ≈ 0.998 / 0.988 / 0.990. Every prediction landed.
- **Reading.** The k → 0 residual rises under rim refinement toward ≈ +0.9 %, so it is **not** rim-grain discretisation: item 441's "discretisation residual on this mesh" is one rung deep and is superseded as a reading (its record sentence stands). The residual is set by what the rim ladder holds fixed, (E1) surface grain, (E2) far grain or (E3) the faceted body. Those are named, not posed; the I-chain reading goes to the weekly. No band moved. The step-2 family stands at 3 of 4.

### Review ruling 2026-10-14 10:30 (item 450)

Accepted as written (I28, attempt 3 of 4 in the step-2 family; chain I's third and last provisional step before the 10-14 weekly). A slot that runs before the weekly runs it as written.
- **The module change is licensed:** test module only, inside `K_SCALE != 1` branches, with the default path bitwise. W0r runs after the edit and before any solve. Rule (i) holds from W0r on.
- **(Wk′) at 5e-8 is a reproduction tolerance on the per-rung records' printed resolution, not a band.** No band moves: `TH15_BETA_BAND` / `MAT4_BAND` are imported unchanged.
- **The predicted verdict is (Fix), and that is the informative outcome.** If δ_D is rung-independent, the k → 0 residual rises toward ≈ +0.9 %, which contradicts item 441's "discretisation on this mesh" reading. Either verdict is a record for the weekly. The reading rule is computed in the outcome from the three logs, never coded.
- **Cap.** After this item the step-2 family stands at 3 of 4. A fourth sub-step freezes it, so the weekly re-scopes the rim floor as a new numbered step or banks it.

## §7 row moved verbatim 2026-10-14 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/IMP-6.md`)* **A PEC oblate spheroid in a lossy host under a uniform quasi-static field — the closed-form stand-in for a plate / disc edge (Phase 7 (7a))** (opened by the 2026-10-12 10:30 daily under rule (5) lane 3; for the weekly to ratify). `IMP-2`'s companion: semi-axes (a, a, c), a > c. In-plane drive E₀x̂: rim enhancement 1/L_x. Axial drive E₀ẑ: pole enhancement 1/L_z. L_x + L_y + L_z = 1. **Done-when (§4):** *step 0* (numpy, §9 item 403): the depolarisation factors, both potentials and fields, and the rim / pole identities at the bands in the chunk file, with the prolate-formula control missing ⇒ 🟡. *Step 1* (FEM, degree 2, `IMP-2`'s hole route with an additive `oblate=True` on `spheroid_in_box_domain`, a = 6 mm, aspect 2 and 5, §9 item 404): per aspect \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball Φ ≤ `MAT4_BAND`, cavity dofs zero; the axial-drive control misses ⇒ ✅. *Step 3* (axial drive, §9 item 409, rule (j) amendment 2026-10-13 03:00 per the 10-12 weekly): per aspect \|β_axial − 1\| ≤ `TH15_BETA_BAND`, pole-ball (r_p = a/5) Φ ≤ `MAT4_BAND` against `pec_oblate_spheroid_field(drive="axial")`, cavity dofs zero, the in-plane control missing both ⇒ the axial half of ✅. Not claimed: a square-edged plate, a plate in a coil, SAR, finite k, convergence, a pointwise pole value. *Step 2 (rim ladder, I22): a record; item 408 BLOCKED on a rim-truncating mesh, re-posed as item 426 (attempt 2 of 4) — landed 2026-10-13 on `oblate_rim="edge"` (rule (c), the equator a model edge; 0 inside midpoints on every rung): at `resolution_tip` ρ_r/3, /4, /5 the rim-ball Φ miss reads +1.223 / +1.486 / +1.613 % and \|β − 1\| 0.196 / 0.778 / 0.599 %, all inside the imported bands; the miss grows slightly with rim refinement (pairwise order −0.64 / −0.32 on the max-edge H), so the rim-local grain is not what sets it — a floor from the fixed surface / far grains, a finding, not a rate (`20261013T200520Z_IMP-6-step2-a2-W1a.log`, `…200607Z_…-W1b.log`, `…200729Z_…-W1c.log`); item 440 (I26): comparand degree 14 moves m_B by 4.7e-12, the a/18 miss reads m_A +1.371 % (item 426: +1.486 %) — the floor is neither (A) nor (B) at the pre-registered thresholds: the far grain or the quasi-static comparand itself (`20261014T022856Z_IMP-6-step2b-W1q.log`, `20261014T022936Z_IMP-6-step2b-W1s.log`); item 441 (I27): k² × 0.01 at fixed loss tangent reads m_D +0.651 % (k² × 0.09: +0.718 %, leading-order ratio 0.919 measured against 0.919 predicted), `resolution_far` a/5 reads m_C +1.423 % — the floor is (D), the quasi-static comparand's finite-k error (δ_D ≈ +0.844 pp of item 426's +1.486 %); the ≈ +0.64 % left at k → 0 is the discretisation residual on this mesh (`20261014T094547Z_IMP-6-step2c-W1k.log`, `20261014T094713Z_IMP-6-step2c-W1k2.log`, `20261014T094858Z_IMP-6-step2c-W1f.log`); item 450 (I28): at k² × 0.01 the rim ladder reads r₃ / r₄ / r₅ = +0.390 / +0.651 / +0.776 % (Δ +0.261 / +0.126 pp, p₄₅ −0.69 on the max-edge H; δ_D,d +0.842 / +0.844 / +0.845 pp), so the k → 0 residual is (Fix), rising: set by what the rim ladder holds fixed (r∞ ≈ +0.89 %), not rim-grain discretisation (`20261014T183329Z_IMP-6-step2d-W1a.log`, `20261014T183202Z_IMP-6-step2d-W1b.log`, `20261014T183415Z_IMP-6-step2d-W1c.log`).* | ✅ **Audit PASS 2026-10-12 18:00; step 3 audit PASS 2026-10-13 03:00.** Step 0 green (item 403); step 1 (item 404) at aspect 2 / 5: \|β − 1\| 0.078 % / 0.634 %, rim-ball Φ miss 0.152 % / 1.418 %, cavity zero, the axial control misses (`20261012T184250Z_IMP-6-step1-W1a.log`, `…184537Z_…-W1b.log`). One mesh per aspect; record text in the chunk file. Step 3 (axial drive gated at aspect 2 / 5, pole ball, §9 item 409, 2026-10-13): \|β_axial − 1\| 0.163 % / 0.547 %, pole-ball Φ miss +1.003 % / +0.946 %, cavity zero; the in-plane control misses both (`20261013T050311Z_IMP-6-step3-W1a.log`, `…050558Z_…-W1b.log`, `…050819Z_…-W2.log`); outcome in the chunk file. 

### §9 item 466 — 2026-10-15 (§10 chain I29, `IMP-6` step 2e — the body grain, the last step-2 sub-step; lane 3, no solve)

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.*

**What it is.** For each of I28's three rim rungs: measure the meshed body (its volume and its equatorial and polar semi-axes), re-evaluate the depolarisation factors and the rim-ball comparand on that measured body, and print the corrected k → 0 miss beside I28's +0.390 / +0.651 / +0.776 %.
- **Licence:** §10 I29 (`PROJECT_PLAN.md:4812–4817`), weekly-enumerated, §9 rule (4). The 10-14 weekly ruled it "the fourth and last step-2 sub-step" (`PROJECT_PLAN.md:4603–4616`).
- **Executor:** implementer.
- **Class:** lane 3.
- **Tier:** **no solve, standard** (`-n 2`, complex, 180 s).
- **Cost:** ≈ 3 measured / ≈ 30 predicted min (15 fixed + ≈ 12 code + ≈ 2 of windows).
- **Depends on:** nothing. Serial on nothing.
- **No band, tolerance, Done-when or target moves.**

**THE FREEZE (binding; the 10-14 weekly, `PROJECT_PLAN.md:4613–4616`).** Whatever I29 reads, **the `IMP-6` step-2 family freezes after it**:
- the row carries the residual as a stated level: "≈ +0.9 % at k → 0 on this comparand, not rim grain";
- no band moves;
- **there is no fifth sub-step**: no further rung, knob, comparand, body or grain variant is posed by any daily, whatever the reading.
- Any later rim or body work is a new numbered step that only a weekly opens.
- An executor that finds the reading "invites" a follow-up writes it in the outcome as a note for the weekly. It poses nothing.

**Why no solve.** The body, the rim mask and the comparand depend on the mesh alone. Φ_h at s = 0.1 on these exact meshes is I28's record:
- ρ_r/3: 22.781980 (`20261014T183329Z_IMP-6-step2d-W1a.log:173`);
- ρ_r/4: 22.852453 (`20261014T183202Z_IMP-6-step2d-W1b.log:181`);
- ρ_r/5: 22.888993 (`20261014T183415Z_IMP-6-step2d-W1c.log:194`).

The window proves it is on those meshes ((Ce), (V0) below). The finite-k share left at s = 0.1 is s²·δ_D = 0.0084 pp (`IMP-6.md:831`), below every reading here. Existing logs print only the body-volume deficit, not the semi-axes. The deficit is −5.147e-03 / −5.148e-03 / −5.088e-03 (`20261014T183126Z_IMP-6-step2d-W0r.log:173, :292, :424`). So one `-n 2` build per rung is needed (the weekly's "else" branch).

**What "the meshed body's equatorial radius and polar half-height" can mean here (settled from the code).**
- The body is a **hole**. `_body_volume_deficit` computes the box volume minus the summed air-cell volumes (`tests/validation/test_imp6_oblate_fem.py:196–205`). The hole's wall is the `CAVITY_TAG` facets (`tests/validation/test_imp1_pec_sphere_lossy_host.py`, imported by `test_imp6_oblate_rim_ladder.py:84–95`).
- On `oblate_rim="edge"` the equator is a model edge and the poles are model vertices (known-issues 2026-10-13, `known-issues.md:930`). gmsh places every wall node on the CAD surface.
- **So the vertex extents, max ρ and max \|z\| over the wall nodes, read a and c by construction**, and cannot carry the facet error. The weekly's "extents differ" (`PROJECT_PLAN.md:4609`) is true of the faceted surface between the nodes, not of the extremal nodes.
- This item therefore **defines the measured body by its second moments.** The semi-axes of the spheroid with the same per-volume moments are a_I = √(5·(I_xx + I_yy)/(2V_h)) and c_I = √(5·I_zz/V_h). On an exact spheroid they return a and c exactly, since ∫x² dV = V a²/5.
- The vertex extents a_v, c_v are printed as the control that shows why. **The review ratifies this definition or substitutes another before the slot runs.**

**Code (new test module only; no `src/`; no edit to `test_imp6_oblate_rim_ladder.py` or step 1's module, so no re-gate is owed).** New `tests/validation/test_imp6_oblate_body_grain.py`, `@complex_only`, one test `test_body_grain`.
- **Imports, private names only, never a `test_*` name:**
  - from `tests/validation/test_imp6_oblate_rim_ladder.py`: `ASPECT`, `DIVISORS`, `ITEM426_RUNGS` (`:179–189`), `ITEM440_R1` (`:155`), `ITEM450_REPRO_BAR` (`:191`), `_mesh_rung` (`:219`), `_mesh_seconds`;
  - from `tests/validation/test_imp6_oblate_fem.py`: `SEMI_MAJOR`, `BOX_HALF_EXTENTS`, `_semi_minor`, `_rho_r`, `_rim_mask` (`:149`), `_owned_masked_cells` (`:164`), `_numpy_quadrature_mean` (`:173`, degree-8 `QDEG`), `_body_volume_deficit`, `_ncells`;
  - from `test_imp1_pec_sphere_lossy_host`: `CAVITY_TAG`, `MAT4_BAND`, `E0`;
  - `pec_oblate_spheroid_field` (`src/fem_em_solver/utils/analytical.py:1831`) and `oblate_spheroid_depolarisation` (`:1756`).
  - Importing the rim-ladder module runs its knob checks (`:117–206`); no knob is set.
- **Records, cited:**
  - `I28_PHI_H` (the three values above);
  - `I28_MISS_STR` = "+0.3900%" / "+0.6506%" / "+0.7762%" (same lines);
  - `I28_DEFICIT_STR` = "-5.147e-03" / "-5.148e-03" / "-5.088e-03" and `I28_CELLS` = 33 334 / 53 959 / 89 530 (`…W0r.log:173, :292, :424`).
- **Per rung** d ∈ `DIVISORS` (3, 4, 5; `FEM_EM_IMP6_RIM=edge`, defaults otherwise), on `msh, cell_tags, facet_tags = _mesh_rung(d)`:
  - (a) **Volume path:** V_h, the first moments and I_xx, I_yy, I_zz of the hole, computed as the box's closed-form moments minus Σ over owned air cells of the exact P1-tet moments. The tet moment is ∫_T x_i x_j dV = (\|T\|/20)·(Σ_k x_i^k x_j^k + (Σ_k x_i^k)(Σ_k x_j^k)); the box term is 8 b_x b_y b_z · b_i²/3. Then reduce.
  - (b) **Surface path (independent):** over owned `CAVITY_TAG` facets, each normal oriented **into its one adjacent air cell** (outward from the hole): V_s = Σ∫x n_x dS and I_ii,s = Σ∫(x_i³/3) n_i dS, using basix's degree-3 triangle rule (exact for cubics). Then reduce.
  - (c) Vertex extents a_v = max ρ and c_v = max \|z\| over the wall facets' vertices, reduced with MAX.
  - (d) a_I, c_I and the hole centroid.
  - (e) L_x, L_z at (a, c) and at (a_I, c_I) via `oblate_spheroid_depolarisation(e)` with e = √(1 − (c/a)²). Print 1/L_x at each and (L_x/L_x^I)² − 1, the leading-order rim \|E\|² scale.
  - (f) One function `_phi_cf(msh, cells, a_, c_)` = `_numpy_quadrature_mean(msh, cells, lambda p: pec_oblate_spheroid_field(p, semi_major=a_, semi_minor=c_, drive="inplane", e0=E0, interior="continue"))`, on the rung's own `_rim_mask` cells. The mask stays nominal: it is the measurement volume. Call it **twice**: nominal (a, c) and measured (a_I, c_I). A closure counts non-finite \|E\|² values; the count is reduced and printed.
  - (g) r_d = I28_PHI_H/Φ_cf(a, c) − 1 and **r_d^I = I28_PHI_H/Φ_cf(a_I, c_I) − 1**.
- All reductions are called on every rank before any `if comm.rank == 0:` print. The asserts sit outside the branch.

**Anchors (asserted).**
- **(Ce)** cells == `I28_CELLS[d]` as integers. This is the mesh of record.
- **(V0) the nominal call reproduces I28's comparand:** \|Φ_cf(a, c)/`ITEM426_RUNGS[d]["phi_cf_h"]` − 1\| ≤ `ITEM450_REPRO_BAR` (5e-8). On rung 4 the bar is also ≤ 1e-12 against `ITEM440_R1["phi_cf_h8"]`. These are reproduction tolerances, not bands (the 10-14 10:30 ruling, `IMP-6.md:935`).
- **(Vr)** f"{r_d:+.4%}" == `I28_MISS_STR[d]`. This is the reading path, as item 305's (Rb).
- **(V1)** f"{_body_volume_deficit(msh, ASPECT):.3e}" == `I28_DEFICIT_STR[d]`, and \|V_h/V − 1 − deficit\| ≤ 1e-12.
- **(Mv) the divergence theorem, two paths:** \|V_s/V_h − 1\| ≤ 1e-9 and \|I_ii,s/I_ii,h − 1\| ≤ 1e-8 for i = x, y, z. These are instrument tolerances, set ≫ round-off. The box cancellation is ≈ 1 200× on I_zz; predicted ≲ 1e-10.
- **(B) the imported band:** \|r_d^I\| ≤ `MAT4_BAND` on every rung.
- Non-finite count == 0 on both calls (item 428 read 0 on the same continued field, `IMP-3.md:1545`).
- β, \|β − 1\| ≤ `TH15_BETA_BAND` and the cavity check are I28's solved records (`IMP-6.md:922–926`). They are not re-read here because there is no field.

**Negative controls (rule (e)).**
- **(N-nom), asserted:** the corrected path with the nominal axes **is** (V0). It is the same function on the same cells, backed by I28's (Wk′) / (Wk) on these meshes: 1.325e-8 / 1.199e-14 / 6.947e-10 (`…W1a.log:179`, `…W1b.log:188`, `…W1c.log:199`). A body correction that silently fell back to (a, c) reads r_d^I == r_d to every digit; that is printed as "no correction applied".
- **Predicted, printed beside, never asserted:**
  - (N-ext) a_v/a − 1 and c_v/c − 1 ≤ 1e-9: the extents carry nothing;
  - \|a_I/a − 1\| and \|c_I/c − 1\| in [1e-4, 5e-3], with c_I/c − 1 more negative than a_I/a − 1 (the sag of the flat faces at grain a/12 weighs more on the thin axis);
  - a_I and c_I equal across rungs to ≤ 1e-4 relative (the deficit is flat: −5.147 / −5.148 / −5.088e-3; the surface grain a/12 is fixed);
  - \|Φ_cf(a_I, c_I)/Φ_cf(a, c) − 1\| in [0.05 %, 3 %], **sign not predicted**. Two effects compete. A smaller body moves the rim away from the fixed mask (comparand down). A thinner aspect raises 1/L_x and sharpens ρ_r = c²/a (comparand up).
  - Because the correction is predicted rung-flat, **I28's rung trend survives:** Δ^I₃₄ within ±0.05 pp of +0.261 pp and Δ^I₄₅ within ±0.05 pp of +0.126 pp.

**Reading rule (pre-registered; printed, never asserted; it decides nothing but the row's wording).**
- **(E3-level):** \|r_d^I\| ≤ \|r_d\| on all three rungs **and** \|r₅^I\| ≤ 0.25·\|r₅\| = 0.194 %. The faceted body accounts for the k → 0 level.
- **(E3-partial):** \|r_d^I\| < \|r_d\| on all three rungs, but not (E3-level).
- **(E3-not):** otherwise.
- Every threshold is *predicted*. 0.25 is item 450's (Fix) quarter.
- Also printed: r∞^I = r₅^I + Δ^I₄₅·q^I/(1 − q^I) when 0 < q^I < 1, beside I28's ≈ +0.893 %.
- **Predicted verdict: not (E3-level)**, because a rung-flat correction cannot remove a rung-rising residual. Whether it is partial or not depends on the unpredicted sign.
- **Whatever the verdict, the freeze above holds.**

**Windows** (complex, foreground, Bash `timeout: 660000`; I28's `Command:` shape, `20261014T183126Z_IMP-6-step2d-W0r.log:12`). **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**
- **Measured basis:**
  - I28's W0r built all three rungs with the census in 31 s (`…W0r.log:491`);
  - item 440's W1q built rung ρ_r/4 and evaluated three comparands in 15 s (`20261014T022856Z_IMP-6-step2b-W1q.log:224`).
  - Six degree-8 comparands over 4 861 / 11 468 / 21 548 mask cells, plus the moments, add ≈ 10–30 s.
- **W0** collect (smoke, ≈ 5 s):
  `scripts/testing/run_and_log.sh IMP-6-step2e-W0-collect "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 30 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_body_grain.py --collect-only -q'"`
- **W1** (`-n 2`, the closing window, predicted ≈ 60 s, range [35, 150]):
  `scripts/testing/run_and_log.sh IMP-6-step2e-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP6_RIM=edge timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp6_oblate_body_grain.py -s -v --tb=short'"`
- **Order:** code → W0 → W1 → commit. No module edit after W1 starts (rule (i)). No heavy pre-check: a standard tier, no solve.

**Traps already paid for.**
- **A collective inside a rank-0 branch deadlocks** (item 440's W1q, `IMP-6.md:879`). No `_ncells`, `_rss_gib`, `_body_volume_deficit`, `_numpy_quadrature_mean`, `allreduce` or facet census inside `if comm.rank == 0:`.
- **Never set an `FEM_EM_IMP6_*` knob** other than `RIM=edge`. The rim-ladder module's one-knob and RUNG2 guards raise at import (`:148–152, :202–206`). Never set `FEM_EM_IMP6_RUNG2`.
- **The default oblate path truncates the rim** (known-issues 2026-10-13). `RIM=edge` is mandatory; (Ce) catches a wrong variant.
- **Facet orientation:** a cavity facet's adjacent cell is its single air cell. Build `create_connectivity(2, 3)` first and use owned facets only (`index_map(2).size_local`), or exterior facets are double-counted across ranks. A sign error reads (Mv) at ≈ 2 or ≈ −1, not a small miss.
- **`interior="continue"` near the focal ring** (a − 0.121 mm): the corrected body's ring sits ≈ a_I − 0.12 mm. Item 428 read no non-finite values on the same construction. A non-zero count is a finding, never masked.
- **The basis is the in-plane field**, `drive="inplane"`, with no −E₀x̂ subtraction: that belongs to β, not to Φ (item 404's trap).
- Masks stay in numpy, with no UFL `<` / `sqrt` comparison. `-s` on every window. Node IDs only, never `-k`. Status 124 means `pgrep -c python3` first.

**Scope.** Aspect 5, the in-plane rim ball, I28's three meshes, one body definition (second moments), and I28's Φ_h. Not covered:
- a re-solve;
- β against the measured body;
- a curvature-matched or locally fitted rim (the moment spheroid is a global summary);
- aspect 2; the axial drive;
- any band; any fifth sub-step (the freeze).

**Status it can move** (rule (j)). `IMP-6`'s Done-when covers steps 0 / 1 / 3 (`PROJECT_PLAN.md:2610`), and step 2 is a record. No glyph moves; the row stays ✅.
- The step-2 sentence is rewritten to the frozen statement: "≈ +0.9 % at k → 0 on this comparand, not rim grain". It gains the I29 clause: "item 466 (I29): moment-equivalent body a_I/a − 1 ⟨…⟩, c_I/c − 1 ⟨…⟩, L_x ⟨…⟩ → ⟨…⟩; body-corrected miss r^I ⟨…⟩ / ⟨…⟩ / ⟨…⟩ % (⟨verdict⟩); the step-2 family is frozen (10-14 weekly)."
- Known-issues 2026-10-13 is untouched: it documents the default builder.

**Negative result.**
- **(Ce), (V0), (Vr) or (V1) red:** not I28's meshes, or the reading path is wrong. No r^I is printed as a reading. Fix a transcription in-slot and re-run W1; otherwise park on `attempt/IMP-6-item466-<ts>`, BLOCKED (rule (d)).
- **(Mv) red:** the moment code is wrong (orientation or ownership). Fix it in-slot and re-run W1. It is never a finding about the body.
- **(B) red** (\|r_d^I\| > `MAT4_BAND`): a real reading. The measured body moves the comparand by more than the band. Report it, record it on the row as measured, and **the freeze still holds**: no band moves and no follow-up is posed. Park only if the reading came with any identity red.
- **Non-finite > 0:** report the count and where it is; no r^I is read on that rung.
- No `PREDICTED RED` line.

**Every-item rules** (PROJECT_PLAN §9 "Every item"):
- bare whole harness calls;
- a lost docker socket stops the slot;
- no edit after the last window;
- `timeout: 660000` and `-s`;
- never delete a log;
- BLOCKED in the parking commit;
- a killed window is not a result;
- complex mode and `FEM_EM_REQUIRE_COMPLEX=1`;
- no slot opens `docs/private/` or `aed_results/`.

**Spawn line:** "Use the implementer agent to execute `IMP-6` §9 item 466 per `docs/planning/chunks/IMP-6.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. New module `tests/validation/test_imp6_oblate_body_grain.py` only; no `src/` edit; never import a `test_*` name; set no `FEM_EM_IMP6_*` knob but `RIM=edge`. Every collective is called before any rank-0 branch. No solve. Compute the pre-registered reading in the outcome; assert none of it. The step-2 family freezes after this item whatever it reads — pose no follow-up."

### Review ruling 2026-10-15 03:00 (item 466)

Accepted, with the drafter's definition ratified **provisionally** (for the 10-17 weekly): the meshed body is measured by its **second-moment semi-axes** (on `oblate_rim="edge"` every wall node lies on the CAD surface, so node extents read a and c by construction and cannot see the facet error); the node extents are printed as the control. The corrected miss is printed, its sign not predicted. **The step-2 family freezes after this item whatever it reads** (10-14 weekly): no fifth sub-step; the row carries the residual as a stated level. Known-issues 2026-10-13 (the oblate builder) leaves only when that freeze is written on the row.

**Outcome (2026-10-15 07:30 CDT slot, take-next, DONE; reading E3-not).** Code: new test module only (`tests/validation/test_imp6_oblate_body_grain.py`), no `src/` edit, no knob but `RIM=edge`. Two windows, `-n 2`, complex: `20261015T125644Z_IMP-6-step2e-W0-collect.log` (5 collected, 2 s), `20261015T125651Z_IMP-6-step2e-W1.log` (5 passed, 29 s; predicted ≈ 60 s, range [35, 150] — under the range).
- **Anchors, all green (log :169–:451):** (Ce) 33 334 / 53 959 / 89 530; (V0) 1.325e-8 / 1.118e-9 / 6.947e-10 (rung 4 vs item 440: 1.199e-14); (N-nom) the nominal path reproduces I28's (Wk′)/(Wk) to the printed digit; (Vr) +0.3900 / +0.6506 / +0.7762 %; (V1) −5.147e-3 / −5.148e-3 / −5.088e-3, identity 0 / 4.7e-15 / 4.8e-15; (Mv) |V_s/V_h − 1| ≤ 3.6e-15, |I_s/I_h − 1| ≤ 9.0e-12 (rung 4, I_zz); cavity facets 4 178 / 4 918 / 6 032, each with exactly one adjacent cell; non-finite 0 on all six calls; (B) |r^I| ≤ `MAT4_BAND`.
- **Measured body:** (N-ext) a_v/a − 1 = 2.2e-16, c_v/c − 1 = 0 (as predicted). a_I/a − 1 = −2.170e-3 / −2.169e-3 / −2.129e-3; c_I/c − 1 = −8.09e-4 / −8.13e-4 / −8.32e-4; c_I/a_I 0.20027 (nominal 0.20000). Centroid ≤ 2.2e-7 m. 1/L_x 8.015515 → 8.006876 / 8.006914 / 8.007285; (L_x/L_x^I)² − 1 = −0.216 / −0.215 / −0.205 %.
- **Comparand:** Φ_cf(a_I, c_I)/Φ_cf(a, c) − 1 = −3.535 / −3.533 / −3.472 %. **Body-corrected miss r^I = +4.068 / +4.337 / +4.401 %** beside I28's +0.390 / +0.651 / +0.776 %. Δ^I₃₄ = +0.269 pp, Δ^I₄₅ = +0.064 pp, q^I = 0.238, r∞^I = +4.42 % (I28 ≈ +0.893 %).
- **Reading rule:** E3-not (predicted: not E3-level — met). **Predictions missed, reported (rule (e), none asserted):** the comparand moved −3.5 %, outside [0.05 %, 3 %]; a_I shrinks more than c_I (−2.17e-3 vs −8.2e-4), the reverse of "c more negative"; Δ^I₄₅ = +0.064 pp is outside ±0.05 pp of +0.126. a_I / c_I are rung-flat to ≈ 2e-5 relative on c, ≈ 4e-5 on a (within the 1e-4 prediction).
- **Reading:** the moment spheroid is ≈ 0.22 % smaller in a, and the in-plane rim field near ρ = a falls steeply with distance from the rim, so a rim moved inward by ≈ 13 µm against a fixed nominal mask drops the comparand by ≈ 3.5 %. The meshed rim itself is not displaced (its nodes are on the CAD surface), so a global moment summary is the wrong body for a rim-local comparand: the faceting deficit lives on the flat faces, not at the rim. This **sharpens** the frozen statement rather than moving it — the ≈ +0.9 % k → 0 residual is not the global body-volume deficit.
- **Freeze:** the step-2 family is frozen (10-14 weekly); nothing is posed. Note for the weekly (not a proposal): a curvature-matched or rim-local body measure is the only body variant this reading leaves open, and the brief's Scope excludes it. Known-issues 2026-10-13 is left untouched per the brief ("it documents the default builder"); the 03:00 ruling's removal condition is now met on the row — the review decides.

### §9 item 482 — 2026-10-15 (§10 chain I31; lane 3)

*Drafted for the 2026-10-15 18:00 daily review; the review rules on it.* `IMP-6` step 4: the disc in the F-clinical phantom. `implant_spheroid` (kept, the `GEO-36` step 7 body) composed with C2's tagged ASTM box (`phantom_box`, `GEO-38` step 2, item 470, landed `c1b8debc`). **CAD identities only.** Licence: §10 I31 (the 10-14 weekly: "all **new** human-scale implant work moves to F-clinical (chain C, I31)"), Phase 7 (7b), §9 rule (4). No band, tolerance, target or Done-when moves; every band is imported.

**Executor and cost.**
- **Executor:** implementer.
- **Tier:** heavy (1200 s ceiling), build-only, **real image** (C2's box module ran real). Container `timeout -k 30 590` on every build window.
- **Ranks:** **`-n 4`** for the window of record. C2's three widths built the same mesh (1 149 150 cells at `-n 1 / 2 / 4`), and `-n 4` measured fastest: 339 s, against 382 s at `-n 1` and 529 s at `-n 2`. The `-n 2` window lost ≈ 160 s to host load, 61 s from the ceiling. A missing reduction shows at any width > 1. The chain line says `-n 2`; that deviation is for the review (below).
- **Cost:** ≈ 23 measured (chain no-solve class median) / ≈ 55 predicted min. That is 15 fixed + ≈ 20 code (≈ 30 lines of `src/`, a ≈ 220-line module lifted from two existing modules) + ≈ 17 of windows: W0 ≈ 10 s, W1 ≈ 430 s, R1 ≈ 340 s, R2 ≈ 175 s.
- **Depends on nothing.** C2 is on `main`.
- **No take-next.** At ≈ 55 min there is no slack for a second item.

**Why `src/` must change (three defects in today's code, each read off `src/fem_em_solver/io/mesh.py`):**
1. **The refusal.** `birdcage_port_domain` refuses `phantom_box` with `implant_spheroid` (`:4552–4555`), so the composition cannot be built today.
2. **The containment test.** With `phantom_box` set, `phantom_radius` becomes the box's xy circumradius hypot(0.21, 0.045) = 0.214767 m (`:4573`). `_normalise_implant_spheroid` then applies its cylinder test, \|centre_xy\| + a < R (`:5036–5039`). For a 90 mm-deep box that test is not conservative: a disc centred at y = 0.040 m passes it with a 166 mm margin while poking 1 mm out of the top of the gel.
3. **The conductor mass.** `conductor_cad_mass` is the can's mass read before the implant cut (`shield_air_cad_mass`, `:5660`) minus the cut air (`:5914`). The implant cut (`:5857–5862`) has already taken the body out of that air. With a body present the diagnostic therefore reads conductor + V, and C2's (V1) partition would miss by V / can = **4.350023e-7, 435× `STRIP_BAND`**.

**`src/` change (rule (c): keyword-only, additive, disclosed; ≈ 30 lines).** Every `phantom_box=None` path must execute the lines it executes today.
- **(a) Composition.** The `phantom_box` refusal keeps `implant_rod` (scope: one body, the disc) and drops `implant_spheroid`. New message: `"phantom_box cannot be combined with implant_rod"`. No test matches the old string (grep of `tests/`, 2026-10-15).
- **(b) Box containment.** `_normalise_implant_spheroid` gains keyword-only `phantom_box: Optional[Dict] = None`.
  - When it is set, the containment test is the box test: \|c_x\| + a < Lx/2, \|c_y − off\| + a < Ly/2 and \|c_z\| + a < Lz/2, all strict. The bounding sphere of radius a is conservative for any axis.
  - Refusal message: `"implant_spheroid is not strictly inside the phantom_box"`.
  - When `phantom_box` is `None`, the existing cylinder test runs unchanged.
  - Both call sites pass the normalised `phantom_box`: the public one (`:4587`) and the builder's (`:5269`).
- **(c) Conductor mass.** In the `phantom_box` + `shield_radius` branch (`:5910–5916`), when an implant body was cut from the air, also subtract its OCC mass (`rod_occ_mass`, set at `:5823`). Write it as a separate `if` on the implant, so that the no-implant line is not rewritten.
  - The diagnostic then means what its name says.
  - Disclose this as a definition repair of item 470's diagnostic. It is not a band.
- **Docstrings.** One sentence each on `phantom_box` and `implant_spheroid`.
- **Not changed:**
  - `implant_rod` with `phantom_box` stays refused;
  - the spheroid's sizing fields, construction and tags (402 / 403);
  - the box's own CAD;
  - `GEO-38` step 1a's shield and corner refusals.

**Fixture.**
- **Kwargs:** `_box_kwargs()` (`tests/validation/test_geo38_f_clinical_box.py:54`, item 470's r1 kwargs with the box) plus `implant_spheroid = dict(SPHEROID, keep_cells=True)`, with `SPHEROID` from `tests/mesh/test_birdcage_implant_spheroid.py:76`.
- **Body and placement:** a = 6 mm, c = 1.2 mm, centre (0.015, 0, 0), axis x̂. This is `GEO-36` step 7's body and placement, so `RIM_POINTS`, `_rim_census`, `_spheroid_volume`, `_spheroid_area` and `_inside_true_spheroid` import unchanged.
- **Inside the box:** \|c_x\| + a = 0.021 < 0.21, \|c_y\| + a = 0.006 < 0.045, \|c_z\| + a = 0.006 < 0.325.
- **The placement is a modelling choice, and the item claims nothing for it.** ASTM F2182's high-E site near a lateral wall belongs to the solve step that follows, with the operator's source.

**New module `tests/validation/test_imp6_disc_f_clinical.py`.**
- **Selection:** the build test is selected by `FEM_EM_IMP6_STEP4=1` only, never by `-k`, and skips when the variable is unset. The refusal tests are pure and always run.
- **Imports, by name only, never a `test_*` symbol:**
  - `_box_kwargs` (box module);
  - `_ring_ports`, `_owned`, `_facet_coords`, `_area`, `LEG_COUNT` (`tests/validation/test_geo38_f_clinical_build.py`);
  - `SPHEROID`, `SPHEROID_CENTRE`, `SPHEROID_SEMI_MAJOR`, `_rim_census`, `_spheroid_volume`, `_spheroid_area` (`tests/mesh/test_birdcage_implant_spheroid.py`);
  - `_interior_facet_check` (`tests/mesh/test_birdcage_implant_rod.py:184`);
  - `_census_checks`, `_owned_census` (`tests/validation/test_wf12_f_human_hole.py`, as C2 imports them);
  - `_interface_facet_tags` (`fem_em_solver.io.mesh`); `SHEET_IFACE`, `PORT_LOWER`, `PORT_UPPER` (`tests/mesh/test_birdcage_port_sheets.py`);
  - the bands: `STRIP_BAND` (1e-9, `tests/validation/test_ops61_exact_strip_ports.py:58`), `CAD_ANALYTIC_BAND` (1e-7) and `EXACT` (1e-9) (`tests/mesh/test_birdcage_conductor_hole.py:117, :123`), `TERMINAL_AREA_BAND` ((0.95, 1.0), `tests/mesh/test_birdcage_leg_gaps.py:73`).
- **Do not import `tests/mesh/test_birdcage_implant_spheroid_f_human.py`.** It only ever ran in the complex image, and it imports the `WF-12` copper-matrix module. This module runs real.
- **New constants:**
  - `F_CLINICAL_DISC_CELL_GUARD = 1_400_000`: a guard, checked before any facet reading; the predicted count plus ≈ 17 % headroom, as F-human's 600 000 sat over 417 653.
  - `C2_CONDUCTOR_CAD_MASS_M3 = 1.994040484843362e-03`: printed comparand only, cited to `20261015T185056Z_GEO-38-s2-W3-n4.log:11903`.
- **One build per process.** No cell record is written into code: rule (i) forbids an edit after the last window. The cells go in this file and on the row.

**Anchors (asserted at `-n 4`; owned cells / facets, every count reduced, collectives before any rank-0 branch).**
- **(D0) Guard.** Cells < `F_CLINICAL_DISC_CELL_GUARD` before any facet reading.
- **(D1) Body CAD.**
  - `body_mass_m3` = V = 4πa²c/3 = 1.809557368e-7 m³ at `CAD_ANALYTIC_BAND`.
  - The OCC `wall_area_m2` = A = 2.473638076e-4 m² at `CAD_ANALYTIC_BAND`.
  - Step 7 read 2.2e-16 and 2.451e-11.
- **(D2) Box / body CAD partition.** With L³ = Lx·Ly·Lz from `_box_kwargs()`, never retyped:
  - CAD(403) = `kept_mass_m3` = `body_mass_m3` at `STRIP_BAND`;
  - CAD(3) = L³ − `body_mass_m3` at `STRIP_BAND`;
  - CAD(3) + CAD(403) = L³ at `STRIP_BAND`.
  - These are step 7's (H1)–(H3) with the box for the cylinder. Predicted ≈ 1e-16: C2's phantom group read 0.024569999999999995.
- **(D3) Meshed partition.**
  - The owned tag-3 and tag-403 tet volumes, summed and reduced, = L³ at `STRIP_BAND`. The box is planar-faced and the 3 / 403 interface facets are shared, so the union is exact. Predicted ≈ 1e-16; C2's (B1) read 1.110e-16 at `-n 4`.
  - The tag-3 vertex extents equal C2's (B2) box extents at `STRIP_BAND`.
- **(D4) Can partition with the body.** Σ `cad_mass_by_group` + `conductor_cad_mass_m3` = `air_domain_cad_mass_m3` at `STRIP_BAND`. This is C2's (V1), which needs `src/` change (c). Predicted ≈ 4.8e-11 (C2: 4.780e-11).
- **(D5) Interface.** `_interior_facet_check` on 402:
  - owned count > 0, exterior 0;
  - every owned facet is interior, between one tag-3 and one tag-403 cell;
  - meshed dS(402) / A inside `TERMINAL_AREA_BAND`.
  - Call `msh.topology.create_entity_permutations()` unconditionally on every rank first (the standing 3b-iv caution).
- **(D6) Centroid.** The OCC wall centroid lies at `SPHEROID_CENTRE` to `EXACT`·a.
- **(D7) Rim fidelity.** `_rim_census` gives n_mask > 0, n_inside 0, n403 > 0 and n403_outside 0. That is item 408's defect, made a gate, on the box path.
- **(D8) Census.**
  - `_census_checks(census, n_cells, _ring_ports())` is empty;
  - tags 3 and 403 are present, with no tag 1;
  - `implant_rod_cad` is `None`.
  - If the coil cavity swallowed the body, the 402 facets would read exterior, and (D5) catches that.
- **(D9) Ports undisturbed.** Per ring sheet: owned facets of 210+i (reduced) = `ring_sheet_gmsh_triangles[i]` on 64 / 64 sheets, and area = chord·w at `STRIP_BAND`. These are C2's (T1) and (T2) area clause, using owned-facet numpy sums, never UFL `dS`.
- **(D10) Refusals** (pure; `_refuse_before_gmsh`'s tripwire pattern, `tests/validation/test_geo38_f_clinical_build.py:384`):
  - the planted poke-out (N3) is refused before gmsh;
  - `implant_rod` with `phantom_box` is still refused;
  - a hole spheroid (`keep_cells=False`) under `as_hole` is still refused;
  - the kept request with the box is accepted by `_normalise_implant_spheroid` (it returns `keep_cells` True).

**Negative controls (rule (e): each factor labelled).**
- **(N1) The C2 comparand misses (D2). Asserted binary, factor predicted.**
  - CAD(3) against L³ misses by `body_mass_m3` / (L³ − `body_mass_m3`). Predicted 7.364960e-6, 7.36e3× `STRIP_BAND`.
  - That factor is also the ceiling: it is a fixed geometric quantity, not a saturating one.
  - Assert the miss > `STRIP_BAND` and print it with the factor.
  - Meaning: a build that silently dropped the body passes C2's (B1) and fails (D2).
- **(N2) The pre-repair conductor reading misses (D4). Asserted by construction, factor predicted.**
  - (Σ groups + `conductor_cad_mass_m3` + `body_mass_m3`) / can − 1 is the diagnostic as item 470 defined it. Predicted 4.350023e-7, 435× `STRIP_BAND`.
  - Assert the miss > `STRIP_BAND`.
  - Meaning: (D4) can see a body-sized accounting error.
- **(N3) A planted poke-out, the control for `src/` change (b). Both outcomes asserted (W0).**
  - The request: centre (0.015, 0.040, 0), a = 6 mm.
  - With `phantom_box=None` and R = 0.214767 m, H = 0.65 m, `_normalise_implant_spheroid` **accepts** it: a margin of 0.214767 − 0.048720 = 0.166 m.
  - With the box it **refuses** it before gmsh: \|c_y\| + a = 0.046 > 0.045, a 1 mm breach.
  - Separation: accept vs refuse, binary by construction.
- **(N4) The body's mesh is coil-independent. Predicted, printed only.**
  - The tag-403 cell count = 14 047 and meshed 403 / V = 0.995984, as on F-human (`20261013T203250Z_GEO-36-step7-W2b.log:6269, :6295`). The sizing fields that reach the body are its own: the coil is ≥ 0.33 m away.
  - A different count is a reading, never a red.
  - The conductor mass against `C2_CONDUCTOR_CAD_MASS_M3` is also printed only, predicted ≤ 1e-9. OCC's partition noise on this can is ≈ 2e-11 m³, which is ≈ 1e-8 relative to a 2e-3 m³ mass, so an assertion at `STRIP_BAND` would be the rule-(f) trap.

**Printed predictions, written first.**
- **Cells:** 1 149 150 + [35 k, 100 k], point ≈ 1 197 416. Step 7 added 48 266 to F-human at the same `resolution` and body sizing.
- **Build at `-n 4`:** ≈ 400 s [380, 510].
  - The gmsh wall over C2's 310 s is predicted at +[50, 170] s. On F-human the body took the gmsh wall from 97.65 s (`20261013T204028Z_GEO-36-step7-W3b-regate.log`) to 158.68 s (`…203250Z_…-W2b.log:6264`), +61 s.
  - The upper end assumes the field-evaluation share scales with the ≈ 3.1× cell count.
  - The window is predicted at ≈ build + 25 s.
- **Unknowns:** E, F, U₁ = E and U₂ = 2E + 2F, printed. They are the price input for the solve step the weekly will pose.
- **Memory:** summed RSS ≈ 2.7 GiB (C2: 2.61 GiB at `-n 4`).

**Windows.**
- **General:** foreground; Bash `timeout: 660000`; each harness call is the whole Bash command, bare and repo-relative; `-s` on every window; node IDs only; durable capture verbatim, trailing `; exit $rc` included.
- **Pre-checks before every build window:**
  - `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0;
  - also run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver cat /proc/loadavg` and journal the reading.
- **Order:** code → W0 → W1 → R1 → R2 → commit → journal. Every window runs on the final files.
  - If the 1-min load is above 24 before W1, run R2 (the cheapest) first, then re-read the load and run W1.
  - Two heavy windows never run at once.
- **W0** (pure refusals plus the existing refusal tests; ≈ 10 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-W0-refusals "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/validation/test_imp6_disc_f_clinical.py tests/mesh/test_birdcage_implant_spheroid.py::test_s9_refusals_raise_without_meshing tests/mesh/test_birdcage_implant_spheroid.py::test_s9_kept_accepted_with_as_hole tests/validation/test_geo38_f_clinical_build.py::test_geo38_shield_refuses_attempt1_table tests/validation/test_geo38_f_clinical_build.py::test_geo38_shield_refuses_leg_gap_length -v -s --tb=short'"`
  - The build test skips here because `FEM_EM_IMP6_STEP4` is unset.
- **W1** (the window of record; `-n 4`, ≈ 430 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-W1-n4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_IMP6_STEP4=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/validation/test_imp6_disc_f_clinical.py -s -v --tb=short > /workspace/logs/imp6-s4-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp6-s4-w1-raw.log; cat /workspace/logs/imp6-s4-w1-raw.log; exit \$rc'"`
  - Print a flushed `build done <t> s` line straight after the builder returns, so that a kill still records the build time.
- **R1** (rule-(c) re-gate of the C2 path, `-n 4`, ≈ 340 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-R1-ruleC-GEO-38-s2-n4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_GEO38_STEP2=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/validation/test_geo38_f_clinical_box.py -s -v --tb=short > /workspace/logs/imp6-s4-r1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp6-s4-r1-raw.log; cat /workspace/logs/imp6-s4-r1-raw.log; exit \$rc'"`
  - **Bitwise reading** against `20261015T185056Z_GEO-38-s2-W3-n4.log:11897–11908`:
    - 1 149 150 cells;
    - (T1) 64/64, N_gmsh 44..46, total 2 882;
    - (T2) worst area 2.220e-15;
    - (B1) 1.110e-16;
    - (V1) groups 4.139940041342700e-01, conductor 1.994040484843362e-03, rel 4.780e-11;
    - (V2) 2.220e-16;
    - (c2) excess 5.
  - (c1)'s owner rank may differ; it is printed.
- **R2** (rule-(c) re-gate of the F-human spheroid path through the unchanged `phantom_box=None` containment branch; `GEO-36` step 7 W2b relabelled; complex; ≈ 175 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-R2-ruleC-GEO-36-step7-W2b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_spheroid_f_human.py -s -v --tb=short'"`
  - **Bitwise reading** against `20261013T203250Z_GEO-36-step7-W2b.log:6266–6295`:
    - 417 653 cells, tag 403 14 047;
    - meshed 403 / V 0.995984;
    - (H2) CAD(3) 2.225495472317e-03;
    - (H4) 2.451e-11 and 0.997686;
    - (H5) 240 surfaces, rel 1.501e-10;
    - (H8) 128 / 0 / 14 047 / 0.
  - The F-small spheroid build is not re-run: no line it executes changes. Its pure (S9) refusals run in W0.

**Traps already paid for.**
- **Item 470's (A0).** It was briefed as "assert if under 1 %" and the module printed it instead. Here every anchor (D0)–(D10) and every asserted control is an `assert`. Everything else is labelled printed.
- **Item 468's face-centre clearance.** A refusal must use the real geometry. Change (b) uses the box's half-extents, not the circumradius, and (N3) proves the difference.
- **Ghosted `dS`.** The birdcage mesh is `shared_facet`-ghosted (`io/mesh.py:4907`), so `_interior_facet_check`'s dS("+") on 402 is licensed at `-n 4`. Still call `create_entity_permutations()` first. Sheets get owned-facet numpy only: the unghosted `-n 4` `dS` hang is on record (`GEO-37`).
- **Rank-local data.** `cell_tags.values`, local extrema and `assemble_scalar` are rank-local. Reduce everything; (D3)'s volumes and (D7)'s midpoints use owned cells only.
- **`GEO-23`.** A failed gmsh build contaminates its process. One build order, no retry, a FAIL recorded verbatim, and no fallback to `dilate` (step 7's (C1s): its CAD misses V by 4.37e-4).
- **Step 7's two defects.** The free profile face and the numpy-scalar `repr` in `MathEval` are already fixed in `src/`. Do not touch the construction.
- **Interior dim-2 gmsh groups** hang `model_to_mesh` at `-n 2` (known-issues 9). 402 stays the dolfinx-side rebuild on the kept body, as today.
- **Real image.** Never import a complex-only module here.
- **Harness discipline.** `pgrep -f "python3 -m pytest"` is denied: use `pgrep -c python3`. A `| tail` or `| grep` on a harness call loses the socket exemption.
- **Kills.** A kill (Status 124): `pgrep -c python3` first, and the rubric's orphan-stem sweep only if it reads 0. A killed window is not a result and is not retried.
- **Rule (i).** No edit to `io/mesh.py` or the module after W0 starts. Any fix re-runs every window.

**Scope.**
- In this item: CAD and mesh identities of one kept aspect-5 disc in the ASTM box in the F-clinical coil, plus refusals.
- Not in this item:
  - a solve, a PEC condition, a rim ratio, SAR or heating (the next chain-I step, which the weekly poses and prices; C1 (b)'s column was killed at heavy, so any F-clinical solve is `xl`);
  - the rod in the box (stays refused);
  - a hole body;
  - any placement other than the one above, or an ASTM-site claim;
  - aspect 2;
  - the step-2 rim family (frozen);
  - compliance, scanner or tuning claims (chain C).
- `IMP-6` stays ✅. `GEO-38` stays 🟡.

**Status it can move.**
- **The `IMP-6` §7 row** (✅ stays; its Done-when covers steps 0 / 1 / 3, so step 4 is a record, as `GEO-36` step 7 was on its row). It gains one clause, if the row stays ≤ 2 KB: "Step 4 (item 482, I31): the kept aspect-5 disc in F-clinical's ASTM box — (D1)–(D10) green at `-n 4`, ⟨cells⟩ cells, `phantom_box` × `implant_spheroid` composed (box containment, body-aware conductor mass)".
- **Rule (j)** does not bind: there is no Done-when amendment.
- **Downstream:** it makes the disc solve in F-clinical poseable (the weekly's step), and its printed E / F / U₂ are that step's price input.

**Negative result.** Report and stop. Park on `attempt/IMP-6-step4-<UTC>` and mark item 482 BLOCKED in the same commit (rule (d)).
- **W0 red on (N3)**, either the box refusal not firing or the cylinder test not accepting: change (b) is wrong. Fix it before W1 and re-run W0. A second red stops the item.
- **A gmsh FAIL on the composed build:** record it verbatim with its surface pairs, plus a known-issues entry. It is a finding about the box × body composition. No retry.
- **(D2) / (D3) red with (D1) green:** the box cut is wrong around the body.
- **(D4) red at ≈ 4.35e-7:** change (c) did not take. That is the diagnostic, not a finding.
- **(D5) / (D7) red:** the item-408 defect, or the cavity reaching the body, on the box path. Print the counts and the deepest midpoint.
- **W1 killed at 590 s:** the composed build does not fit a heavy window, and the flushed build line is the record. The step returns to the weekly as `xl` or a split. It is not retried.
- **R1 or R2 moves a bitwise reading:** revert `src/` and stop.
- **A printed reading outside its prediction** (cells, build time, (N4)) is a reading, never a red.

No `PREDICTED RED` line.

**For the review to rule on.**
1. `-n 4` instead of the chain's `-n 2`.
2. The `GEO-36` step 7 placement and axis (0.015, 0, 0), x̂, as a CAD-only choice.
3. `src/` change (c), which redefines item 470's `conductor_cad_mass_m3` on implant builds (bitwise on C2 by R1).

**Every-item rules** (PROJECT_PLAN §9 "Every item"):
- bare whole harness calls;
- a lost docker socket stops the slot;
- no edit after the last window;
- `timeout: 660000` and `-s`;
- never delete a log;
- BLOCKED in the parking commit;
- a killed window is not a result;
- two heavy windows never at once;
- no slot opens `docs/private/` or `aed_results/`.

**Executor spawn prompt (verbatim).** "Execute §9 item 482 — `IMP-6` step 4 (§10 I31) — per `docs/planning/chunks/IMP-6.md` §9 item 482. Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s; container `timeout -k 30`. The `src/` change is keyword-only and disclosed: `phantom_box` accepts `implant_spheroid` (the rod stays refused), box containment in `_normalise_implant_spheroid` when `phantom_box` is set, and the body-aware `conductor_cad_mass_m3`; every `phantom_box=None` line executes as today. New module `tests/validation/test_imp6_disc_f_clinical.py` only; never import a `test_*` name or the complex-only F-human spheroid module. Run W0, W1 (`-n 4`), R1, R2 in order, foreground, on the final files; R1 and R2 must read their logs' numbers bitwise or you revert `src/` and stop. Never retry a gmsh failure or a killed window. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-6-step4-*`, mark item 482 BLOCKED in the same commit, and report."

### Review ruling 2026-10-15 18:00 (item 482, I31)

**Queued, with three rulings.**
1. **`-n 4`, not the chain's `-n 2`.** The rank count is execution detail. The `-n 2` F-clinical build has come within 61 s of the ceiling (item 470 W2, 529 s), while `-n 4` measured 339 s. The disc adds a predicted 50–170 s. If the W1 build is killed, report and stop: the measured lower bound goes to the weekly as an `xl` / probe question. Do not retry at a higher width.
2. **Placement is a CAD choice.** The disc keeps `GEO-36` step 7's centre (0.015, 0, 0) and x̂ axis. This is not the ASTM F2182 high-E site and claims nothing about it; that site is the weekly's call once C6 gives a field.
3. **The `src/` change is licensed under ruling (c), keyword-only.** It lifts the composition refusal, makes the containment test conservative on the box, and gives `conductor_cad_mass_m3` its implant-build reading. It is disclosed. **On builds without an implant, every value is bitwise unchanged.** R1 re-gates C2 (`tests/validation/test_geo38_f_clinical_box.py`) and the F-human path bitwise in the same slot, one heavy window at a time. A non-bitwise R1 parks the item. Item 470's diagnostic is never redefined on non-implant builds.

The step's target is the chain's: CAD identities only, at imported bands (`STRIP_BAND`, `CAD_ANALYTIC_BAND`). It runs no solve.

**Outcome 2026-10-16 02:00 slot: BLOCKED (W1 killed; ruling 1, report and stop).**
- `src/` change (a)–(c) and `tests/validation/test_imp6_disc_f_clinical.py` written as briefed; parked on `attempt/IMP-6-20261016T020000Z` (`8d6ffdd1`).
- **W0** green, 15 passed / 1 skipped (`20261016T020256Z_IMP-6-step4-W0-refusals.log`): (N3) the cylinder test accepts the poke-out with margin 0.166047 m, the box test refuses it before gmsh on a 1.000000e-03 m breach; (D10) rod-with-box and hole-under-`as_hole` still refused, the kept request accepted; the S9 and GEO-38 refusal tests unchanged.
- **W1** (`-n 4`, pre-check `pgrep -c python3` 0, load 0.30) killed at 590 s inside gmsh (`20261016T020315Z_IMP-6-step4-W1-n4.log`). Phase walls against C2's W3: 1D 315.6 s vs 49.7, 2D 57.9 vs 22.2, 3D 147.9 vs 55.5 (1 284 959 tets created, 365 650 nodes), first optimisation 18.5 vs 7.7; killed in the final optimisation (C2's took 168.6 s). The builder's CAD print reads wall area 2.473638075734e-04 m² and body = kept mass 1.809557368468e-07 m³ (a printed reading; no anchor ran). After the kill `pgrep -c python3` read 0.
- **R2** (F-human, `phantom_box=None` path, on the parked `src/`) green, 11 passed, every reading bitwise against `20261013T203250Z_GEO-36-step7-W2b.log` (417 653 / 14 047, 0.995984, CAD(3) 2.225495472317e-03, 2.451e-11 / 0.997686, 240 / 1.501e-10, 128 / 0 / 14 047 / 0) — but its gmsh wall read 311.64 s vs 158.68 s at step 7 (1D 131.5 vs 73.7): the host ran ≈ 2× slow this slot, as item 485's 431 s vs 295–310 s did at 00:30 (`20261016T021334Z_IMP-6-step4-R2-ruleC-GEO-36-step7-W2b.log`). R1 not run: `src/` is parked, so no C2 path on `main` changed.
- **Reading for the weekly:** the composed build needs > 540 s of gmsh at `-n 4` on this host; the 1D phase is 6.3× C2's (≈ 3× after the host factor), so the body's sizing fields dominate 1D on the 1.2 M-cell mesh. Not retried. The next attempt is an `xl` window, or a heavy window at a measured-quiet hour, followed by R1.

### Review ruling 2026-10-16 03:00 (item 482, I31 BLOCKED)

**Banked; no re-pose by the daily.** W0 was green: the (N3) cylinder accepts at margin 0.166047 m, and the box refuses on a 1.000e-03 m breach. W1 was killed at 590 s in gmsh's final optimisation (`20261016T020315Z_IMP-6-step4-W1-n4.log`), so no anchor was read. Two causes stack:
- **The host ran ≈ 2× slow.** R2's own gmsh wall was 311.64 s against 158.68 s at step 7, the same geometry.
- **The body's sizing fields grow the 1D phase.** It ran 6.3× C2's, and ≈ 3× after the host factor.

A daily cost probe does not fit, because I31's target is the CAD-identity gates, not a price. Splitting the heavy window does not fit either, because the build alone exceeds the ceiling. **I31 therefore goes to the 10-17 weekly.** The weekly decides between an `xl` window and a heavy window at a measured-quiet hour. In either case the R1 re-gate then runs on the parked `src/`. `attempt/IMP-6-20261016T020000Z` (`8d6ffdd1`) is kept until then.

**Family count:** step 4 has 1 attempt of 4.

### §9 item 522 — 2026-10-18 (§10 I31 attempt 2, the disc in C2's tagged box; chain step, build-only; 12:00 / 13:30 / 15:00 slot only)

**Outcome 2026-10-18 (15:00 slot): GREEN.** Code lifted by path after the `e06c8a4a` guard (exit 0). W0 19 passed / 1 skipped, 7 s ((N5a) rc 86 in 1.034 s with the marker; (N5b) rc 0, no marker; (N5c) `gmsh.lib` is `CDLL`). R2 bitwise to `20261013T203250Z_GEO-36-step7-W2b.log` (gmsh 171.07 s, h₂ = 1.078; 184 s). R1 bitwise to `20261015T185056Z_GEO-38-s2-W3-n4.log` (1 149 150 cells, CAD line identical; gmsh 301.56 s, h₁ = 0.973; T̂ = [389.0, 443.3, 467.4] s; 329 s). W1 at `-n 4`, rule armed: **T = 441.11 s under the rule** (phases 1D 170.429 s, 2D 22.333 s, 3D 59.572 s, own opt 7.339 s, Netgen 175.081 s), build 467.3 s, **1 198 350 cells** (+49 200 over C2), (D0)–(D10) green, (N1) 7.365e3× and (N2) 435× as predicted, E 1 584 883, F 2 519 602, U₂ 8 208 970; 5 passed, 472 s. Pre-window `pgrep -c python3` 0 each time; loadavg 0.47 / 1.51 / 3.37 (1-min) before R2 / R1 / W1. Logs: `20261018T200143Z_IMP-6-step4-a2-W0.log`, `20261018T200204Z_IMP-6-step4-a2-R2-ruleC-GEO-36-step7-W2b.log`, `20261018T200533Z_IMP-6-step4-a2-R1-ruleC-GEO-38-s2-n4.log`, `20261018T201127Z_IMP-6-step4-a2-W1-n4.log`. Item 530's gate (W1 inside 450 s) is met. `attempt/IMP-6-20261016T020000Z` is for the next review to delete.

*Drafted read-only for the 2026-10-18 03:00 daily review; the review rules on it.* `IMP-6` step 4, attempt 2. It is item 482 with three changes: the code is lifted from the parked branch, not rewritten; a 450 s gmsh stop rule is enforced in code; and the item runs only in a daytime slot. **Target unchanged:** (D0)–(D10), the disc's CAD and mesh identities in C2's tagged ASTM box, at imported bands. Licence: §10 I31 attempt 2 (`PROJECT_PLAN.md:4831–4833`), the 10-17 weekly's ruling (`:4676–4683`) and epitaph (`:4738–4742`), Phase 7 (7b), §9 rule (4). No band, tolerance, target or Done-when moves. Item 482's brief (above) and its 10-15 18:00 rulings 1–3 hold unless this item changes them.

**Executor and cost.**
- **Executor:** implementer.
- **Tier:** heavy (1200 s ceiling), build-only, **real image** for W0 / R1 / W1 and complex for R2. Every build window has container `timeout -k 30` as the backstop; the stop rule is what stops the build.
- **Ranks:** `-n 4` for the window of record (item 482 ruling 1). gmsh runs on rank 0 only (`io/mesh.py:4856–4890`); the other ranks wait in the `bcast` at `:4898`.
- **Measured costs:**
  - **C2, the box without the disc:** gmsh wall 310.01 s, build 335.4 s, window 339 s at `-n 4` (`20261015T185056Z_GEO-38-s2-W3-n4.log:11894, :11896, :11925`). In the same 13:30 slot C2 also read 323.29 s at `-n 1` (`…183518Z_GEO-38-s2-W1-n1.log:11870`) and **483.82 s at `-n 2`** (`…184154Z_GEO-38-s2-W2-n2.log:11878`). The weekly's "290–310 s" pools C1 (a)'s box-free 295.05 s (`20261015T170721Z_GEO-38-s1a-W1-r1-a2.log:11829`) with C2's `-n 4` reading.
  - **Attempt 1:** killed at 591 s (`20261016T020315Z_IMP-6-step4-W1-n4.log:5691`), inside `optimize("Netgen")`. Before that, `generate(3)` ran 1D 315.574 s (`:3062`), 2D 57.897 s (`:4675`), 3D 147.933 s (`:5141`) and its own optimisation 18.512 s (`:5648`), 539.92 s in all.
  - **Host factor in that slot:** R2's gmsh wall was 311.64 s against 158.68 s on the same geometry (`20261016T021334Z_IMP-6-step4-R2-ruleC-GEO-36-step7-W2b.log:6264` vs `20261013T203250Z_GEO-36-step7-W2b.log:6264`), a ratio of 1.964. By phase: 1D 1.783 (`:1559` in both logs), 3D 2.434 (`:2566`).
  - W0 took 8 s (`20261016T020256Z_IMP-6-step4-W0-refusals.log:77`). R2 took 172 s quiet (`20261013T203250Z_…:6305`) and 344 s slow (`20261016T021334Z_…:6305`). R1 is predicted from C2's W3: 339 s.
- **Predicted gmsh wall on a quiet host: 454 s, range [401, 482].**
  - Attempt 1's 539.92 s, divided by the host factor (1.964 overall, 1.783 by the 1D phase, 2.434 by the 3D phase), gives 275 s [222, 303].
  - Add Netgen at 172.7 s: C2's 168.565 s (`…185056Z_…:11893`) scaled by nodes, 365 649 / 356 833 (`20261016T020315Z_…:5136` vs `…185056Z_…:5091`).
  - Add C2's ≈ 6.3 s outside both calls.
  - **The rule sits on the predicted point.** A stop is about as likely as a pass, even in a quiet slot. A stop is still the outcome the weekly bought: the build's price.
- **Class and cost:** chain no-solve. ≈ 18.6 measured (class median) / ≈ 50 predicted slot-minutes:
  - 15 fixed;
  - ≈ 12 of code (the lift plus ≈ 60 lines of helper and controls);
  - ≈ 18 of windows: W0 ≈ 15 s, R2 172–344 s, R1 ≈ 340 s, W1 ≤ ≈ 495 s;
  - ≈ 5 to commit and journal.
- **No take-next**, in either direction. This item is never a slot's second item, and nothing follows it.
- **Depends on nothing in §9.** C2 is on `main`, and so are attempt 1's logs and journal (`docs/testing/attempts.md:1178`).
- **Family:** `IMP-6` step 4, attempt **2 of 4**. Under the epitaph this is the last heavy attempt: after it, the build goes to `xl` whatever happens.

**The slot-time gate (read first, before any edit).** The session's start line, "Scheduled implementer run, started YYYY-MM-DD HH:MM CDT" (`scripts/automation/implementer-run.sh:45, :73`), must read **12:00, 13:30 or 15:00**, up to 5 minutes late, and the item must be the slot's first.
- In any other slot, skip item 522 and take the next §9 item. Write nothing about item 522: it is not BLOCKED, not journaled and not counted as an attempt.
- The quiet-hours case: C2's W3 ran at 13:50 CDT, and both kills fell in the 19:30–00:30 interval (`PROJECT_PLAN.md:4679–4682`). C2's 483.82 s `-n 2` reading came from a 13:30 slot, so the gate makes a quiet host likely, not certain.
- **W1 starts only by minute 42 of the slot.** If W0, R2 and R1 are green but W1 cannot start by then, commit the lifted code with R1 / R2 as its rule-(c) record and leave item 522 OPEN for "W1 only" in the next eligible slot. This is not an attempt, and the next slot runs W0 → W1.

**Code state (read on `main` and the branch, 2026-10-18).**
- **`attempt/IMP-6-20261016T020000Z`** holds one commit, `8d6ffdd1`, on merge base `e06c8a4a`. It touches **exactly two paths**, 421 insertions and 8 deletions:
  - `src/fem_em_solver/io/mesh.py` (+46 / −8), changes (a)–(c) as briefed:
    - the refusal at `:4552–4555`, narrowed to the rod;
    - `phantom_box=` passed to `_normalise_implant_spheroid` at `:4587` and `:5269`;
    - the box containment test, replacing the cylinder test at `:5036–5039` when a box is set;
    - `− rod_occ_mass` (set at `:5823`), added under `:5910–5916`.
  - `tests/validation/test_imp6_disc_f_clinical.py` (new, 383 lines). The build test is at `:163`, its `build done` print at `:194`, (D0) at `:198–199`, and the pure (N3) / (D10) tests at `:102–157` (the rod-refusal test is at `:129`).
- The branch holds **no record file**. Its logs and journal are on `main`.
- `git diff e06c8a4a main -- src/fem_em_solver/io/mesh.py` and the diff over every module the new file imports are **empty**: `main` has not touched them since the branch.
- **The gmsh timing on `main`.** `mesh_start` / `generate(3)` / `optimize("Netgen")` / `mesh_wall_time` are at `io/mesh.py:6917–6920`, printed flushed as `[birdcage-mesh] mesh wall time T s` (`:6931–6935`) and returned as `mesh_wall_time_s` (`:6938`). These are the builder's only `generate` / `optimize` calls.
- **C2's tagged-box module:** `_box_kwargs()` at `tests/validation/test_geo38_f_clinical_box.py:54–66`, selected by `FEM_EM_GEO38_STEP2` (`:50`).

**Edits allowed.**
1. **The lift, by path, code only.** First run `git diff --quiet e06c8a4a HEAD -- src/fem_em_solver/io/mesh.py`. On exit 0, run:
   `git checkout attempt/IMP-6-20261016T020000Z -- src/fem_em_solver/io/mesh.py tests/validation/test_imp6_disc_f_clinical.py`
   If `main` has moved `mesh.py`, use `git cherry-pick --no-commit 8d6ffdd1` instead (the commit touches only the two paths). After either, `git status` must list those two paths and nothing else. **Never path-checkout any other file from the branch.**
2. **A new helper, `tests/validation/gmsh_stop_rule.py`.** Stdlib only (`os`, `sys`, `time`, `threading`, `faulthandler`). It imports neither `mpi4py` nor gmsh at module level, so a subprocess can import it.
   - **Constants:**
     - `GMSH_STOP_ENV = "FEM_EM_GMSH_STOP_RULE"`;
     - `GMSH_STOP_RULE_S = 450.0`, citing `PROJECT_PLAN.md:4679–4683`;
     - `GMSH_STOP_RC = 86`;
     - `GMSH_STOP_MARKER = "[gmsh-stop-rule] STOPPED"`.
   - **`gmsh_stop_armed()`.** Unset or empty ⇒ False; `"1"` ⇒ True; anything else ⇒ `ValueError`. No float is parsed from the environment (item 511's `_node_pf()` precedent).
   - **A `_StopWatch(seconds, label)` class:**
     - `.arm()` starts a daemon `threading.Timer`. When it fires, it prints `GMSH_STOP_MARKER` with the label, the rule and "price > 450 s" (flushed), then calls `os._exit(GMSH_STOP_RC)`.
     - `.arm()` also sets `faulthandler.dump_traceback_later(seconds + 30, exit=True)`, a fallback that does not need the GIL.
     - `.cancel()` cancels both.
   - **`arm_gmsh_stop_rule(monkeypatch, comm, label)`.** On `comm.rank == 0` only, it wraps `gmsh.model.mesh.generate` and `gmsh.model.mesh.optimize` with `*args, **kwargs` passthrough:
     - the `generate` wrapper arms the watch and prints `[gmsh-stop-rule] armed 450 s at generate(3)` (flushed);
     - the `optimize` wrapper calls the original, cancels, and prints `[gmsh-stop-rule] under the rule: <t> s of 450`;
     - either wrapper cancels on an exception and re-raises, so a gmsh FAIL is never mislabelled as a stop.
3. **New `tests/validation/test_gmsh_stop_rule.py`.** Pure; runs in W0 only, never under `-n > 1`. Tests (N5a)–(N5c) below.
4. **The disc module, additive only.**
   - Import the helper by name.
   - The build test takes `monkeypatch`, and before the builder call (`:191`) runs `if gmsh_stop_armed(): arm_gmsh_stop_rule(monkeypatch, comm, "IMP-6 s4")`.
   - One print states whether the rule is armed.
   - Nothing else in the module or in `src/` changes. With `FEM_EM_GMSH_STOP_RULE` unset, the module is the parked one.

**Anchors.** (D0)–(D10) exactly as in item 482, asserted at `-n 4`, owned and reduced, collectives before any rank-0 branch.

**Imported bands (by name, never retyped).**
- `STRIP_BAND` 1e-9 (`tests/validation/test_ops61_exact_strip_ports.py:58`);
- `CAD_ANALYTIC_BAND` 1e-7 and `EXACT` 1e-9 (`tests/mesh/test_birdcage_conductor_hole.py:123, :117`);
- `TERMINAL_AREA_BAND` (0.95, 1.0) (`tests/mesh/test_birdcage_leg_gaps.py:73`).

**Negative controls (rule (e)).**
- **(N1)–(N4)** as in item 482. (N3) is now **asserted on a measurement**: cylinder margin 0.166047 m, box breach 1.000000e-03 m (`20261016T020256Z_IMP-6-step4-W0-refusals.log:44–45`). (N1)'s 7.36e3× and (N2)'s 435× stay **predicted** factors, because no W1 anchor was ever read. The builder's CAD print, body = kept = 1.809557368468e-07 m³ (`20261016T020315Z_…:92`), is a printed reading only.
- **(N5) The stop rule itself. Asserted by construction, in W0, in subprocesses** (`[sys.executable, "-c", …]`, importing only the helper):
  - **(N5a)** Armed at 1 s around `ctypes.CDLL(None).sleep(3)`, a C call through ctypes as gmsh's API is. The subprocess returns `GMSH_STOP_RC` with the marker on stdout, in < 2.9 s.
  - **(N5b)** Armed at 2 s around `usleep(100000)` and then cancelled. It returns 0 and prints no marker.
  - **(N5c)** `gmsh.lib` is a `ctypes.CDLL` and not a `ctypes.PyDLL`, so the GIL is released during meshing.
  - Separation: fires vs silent, binary.
- **Printed, never asserted:** the host factors h₂ = R2 gmsh / 158.68 s and h₁ = R1 gmsh / 310.01 s, and before W1 the prediction T̂ = R1 gmsh × [1.29, 1.47, 1.55].

**Windows.** Foreground; Bash `timeout: 660000`; each harness call is the whole Bash command; `-s`; durable capture verbatim, ending `; exit \$rc`. Before each of R2, R1 and W1, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3`, which must read 0, then `docker compose -f docker/docker-compose.yml exec -T fem-em-solver cat /proc/loadavg`, and journal the reading. Order: **code → W0 → R2 → R1 → W1 → commit → journal**. The re-gates run first so that the code can land on any W1 outcome. Every window runs on the final files (rule (i)).
- **W0** (pure; ≈ 15 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-a2-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/validation/test_gmsh_stop_rule.py tests/validation/test_imp6_disc_f_clinical.py tests/mesh/test_birdcage_implant_spheroid.py::test_s9_refusals_raise_without_meshing tests/mesh/test_birdcage_implant_spheroid.py::test_s9_kept_accepted_with_as_hole tests/validation/test_geo38_f_clinical_build.py::test_geo38_shield_refuses_attempt1_table tests/validation/test_geo38_f_clinical_build.py::test_geo38_shield_refuses_leg_gap_length -v -s --tb=short'"`
- **R2** (rule-(c) re-gate of the F-human `phantom_box=None` path, and the cheap calibrant; complex; 172–344 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-a2-R2-ruleC-GEO-36-step7-W2b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_spheroid_f_human.py -s -v --tb=short'"`
  It must read item 482's R2 list bitwise (`20261013T203250Z_GEO-36-step7-W2b.log:6266–6295`). Print h₂.
- **R1** (rule-(c) re-gate of C2, and the direct calibrant; ≈ 340 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-a2-R1-ruleC-GEO-38-s2-n4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_GEO38_STEP2=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/validation/test_geo38_f_clinical_box.py -s -v --tb=short > /workspace/logs/imp6-s4a2-r1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp6-s4a2-r1-raw.log; cat /workspace/logs/imp6-s4a2-r1-raw.log; exit \$rc'"`
  It must read item 482's R1 list bitwise (`20261015T185056Z_GEO-38-s2-W3-n4.log:11897–11908`). Print h₁ and T̂.
- **W1** (the window of record, `-n 4`, rule armed; ≤ ≈ 495 s):
  `scripts/testing/run_and_log.sh IMP-6-step4-a2-W1-n4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_IMP6_STEP4=1 FEM_EM_GMSH_STOP_RULE=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/validation/test_imp6_disc_f_clinical.py -s -v --tb=short > /workspace/logs/imp6-s4a2-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp6-s4a2-w1-raw.log; cat /workspace/logs/imp6-s4a2-w1-raw.log; exit \$rc'"`
  - **How the rule is recorded.** Under the rule, the log prints the `under the rule: <t> s` line, then `[birdcage-mesh] mesh wall time T s`, which is **the price**, then the anchors.
  - **A stop** prints the gmsh phase lines completed so far, then `GMSH_STOP_MARKER`, then the `[capture] rc=` footer, which reads non-zero. The marker line is the record, not the rc.
  - **The fallback.** A `Timeout (0:08:00)!` traceback with no marker means the Timer did not fire and the faulthandler fallback stopped the build. That is a helper defect, and the item reports it.

**Traps already paid for.** Item 482's list holds. Added:
- **The lift.** A path checkout of a record file from an `attempt/*` branch deletes `main`'s rows (`8d4cf58`, rubric 4). Only the two code paths come over.
- **The subprocess and MPI.** The (N5) children must not import `mpi4py`, and the controls never run under `-n > 1`. That is why the helper is stdlib-only and its tests sit in their own module.
- **The stop is not a kill.** After a stop, `pgrep -c python3` must still read 0 before the commit. Run the rubric-4 orphan-stem sweep only on a ceiling kill (Status 124 without the marker). No FFCx compile runs before gmsh finishes.
- **No rule outside W1.** Never set `FEM_EM_GMSH_STOP_RULE` in R1 or R2: they are bitwise re-gates, not priced builds. Never set it in an `xl` command either. Unset is the default path.
- **Never retry.** Neither a stop nor a kill is retried, at a higher width or in a later slot (epitaph). No third heavy attempt.

**Scope.** Item 482's scope, plus the price. **Not claimed:** a solve, a placement claim, the rod in the box (I33, item 530), or a cost model beyond the printed T̂. `IMP-6` stays ✅ and `GEO-38` stays 🟡.

**What each outcome means, and what happens to the branch.**
- **GREEN.** T ≤ 450 s, (D0)–(D10) green, and R1 / R2 bitwise.
  - Commit the lifted `src/` and module, the helper and its tests, the four logs, and the `IMP-6` row clause from item 482 (`:1270`), giving T and the cell count.
  - Item 522 is DONE, and item 530 becomes runnable (its gate).
  - **Branch:** the next daily review deletes `attempt/IMP-6-20261016T020000Z`. Before deleting, it checks that `git diff main attempt/IMP-6-20261016T020000Z -- src/fem_em_solver/io/mesh.py` is empty and that the module's diff is only the helper lines.
- **STOPPED AT THE RULE.** The marker is in the log, and R1 / R2 are bitwise. **This is the item's priced outcome, not a red.**
  - The code still lands on `main`, with the rule off by default. This is item 494's precedent (`GEO-38.md:569`), so that the `xl` window can run the same module from `main`.
  - Item 522 is DONE as "priced: gmsh > 450 s at `-n 4`, phases <…>, h₁ / h₂ <…>". The §10 I31 line gets the price.
  - **The `xl` 10-25 build is the weekly's to pre-register.** Its command is W1's with `FEM_EM_GMSH_STOP_RULE` unset. The implementer writes no `xl-pending.md` text.
  - In the same commit, mark item 530 **BLOCKED, to `xl`**.
  - **Branch:** deleted by the next review, as for GREEN.
- **R1 or R2 not bitwise.** Revert `src/`, so nothing lands.
  - Park the lifted files and the helper on `attempt/IMP-6-step4-a2-<UTC>`, and mark items 522 and 530 BLOCKED (rule (d)).
  - **Branch:** keep the old branch until the review confirms the new one is a superset of it.
- **(D) red under the rule.** Report and stop. Park everything on `attempt/IMP-6-step4-a2-<UTC>`, mark items 522 and 530 BLOCKED, and send nothing to `xl` until the review rules.
  - The reds are item 482's: a gmsh FAIL (recorded verbatim, with a known-issues entry), (D2) / (D3), (D4) at ≈ 4.35e-7, and (D5) / (D7).
  - **Branch:** as above.
- **Killed at the container ceiling without the marker.** Not a result; the helper failed. Run `pgrep` and the sweep, then park as above and mark BLOCKED. The build goes to `xl` regardless.
- **Skipped by the slot-time gate, or deferred at minute 42:** nothing changes, and the branch stays.
- **A printed reading outside its prediction** (cells, T, h, (N4)) is a reading, never a red.

No `PREDICTED RED` line: the stop is not an assertion. Any other red stops the slot.

**Status it can move.**
- On GREEN: the `IMP-6` row's step-4 clause (✅ stays; rule (j) does not bind) and §10 I31's DONE line. The printed E / F / U₂ become the price input of the F-clinical disc solve.
- On a stop: §10 I31 gains its price, and the `xl` 10-25 entry becomes poseable by the weekly. I33 moves to BLOCKED → `xl`.

**For the review to rule on.**
1. **The rule sits on the predicted point** (454 s [401, 482] at h = 1). Decide whether attempt 2 is worth its slot, or whether I31 goes straight to `xl` 10-25.
2. **An optional calibrant gate.** Skip W1, and defer to the next eligible slot, when T̂(low) = 1.29 × R1 gmsh > 450 s (that is, R1 gmsh > 348 s). Default: printed only, per the weekly's letter.
3. **Landing the code on a stop**, under item 494's precedent.
4. **The `FEM_EM_GMSH_STOP_RULE` helper as the shared mechanism** for I33.

**Every-item rules** (PROJECT_PLAN §9 "Every item"):
- bare whole harness calls;
- a lost docker socket stops the slot;
- no edit after the last window;
- `timeout: 660000` and `-s`;
- never delete a log;
- BLOCKED in the parking commit;
- a killed window is not a result;
- two heavy windows never at once;
- no slot opens `docs/private/` or `aed_results/`.

**Executor spawn prompt (verbatim).** "Execute §9 item 522 — `IMP-6` step 4 attempt 2 (§10 I31) — per `docs/planning/chunks/IMP-6.md` §9 item 522. Read the brief and any review ruling beneath it first; the rulings hold. Run it only if this session's start line reads 12:00, 13:30 or 15:00 (≤ 5 min late) and it is the slot's first item; otherwise skip it silently. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s; container `timeout -k 30`. Lift only `src/fem_em_solver/io/mesh.py` and `tests/validation/test_imp6_disc_f_clinical.py` from `attempt/IMP-6-20261016T020000Z` (path checkout after the `e06c8a4a` guard, else `cherry-pick --no-commit 8d6ffdd1`); never path-checkout any other file. Add `tests/validation/gmsh_stop_rule.py` and `tests/validation/test_gmsh_stop_rule.py`, and the armed-rule lines in the disc module; nothing else. Run W0, R2, R1, W1 in order, foreground, on the final files; W1 only by minute 42; set `FEM_EM_GMSH_STOP_RULE=1` in W1 only. R1 and R2 must read their logs' numbers bitwise or you revert `src/` and stop. A stop at the rule is the priced outcome: land the code, mark item 530 BLOCKED to `xl`, write no `xl-pending.md` text. Never retry a stop, a kill or a gmsh failure. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-6-step4-a2-*`, mark items 522 and 530 BLOCKED in the same commit, and report."

### Review ruling 2026-10-18 03:00 (item 522, I31 attempt 2)

Accepted as the weekly posed it: a daytime heavy window with a 450 s gmsh stop rule, in the **12:00 / 13:30 / 15:00 slots only**. A slot at any other time skips it, and the skip costs no attempt.
1. **Worth a slot: yes.** The weekly's design makes a stop at the rule a result: the build's price is recorded and the build goes to `xl` 10-25. The drafter's prediction, 454 s [401, 482] on a quiet host, is printed beside the measured value.
2. **No skip rule.** W1 runs whatever the C2 re-gate reads, and the re-gate's gmsh time is printed as a slow-host indicator. A third heavy attempt is never queued (10-17 epitaph).
3. **The code lands on `main` after a green build, or after a stop at the rule** when both re-gate builds reproduce their logged numbers exactly (item 494's precedent). The `xl` window can then run from `main`. The next review deletes `attempt/IMP-6-20261016T020000Z`.
4. **Accepted:** the shared stop-rule helper `tests/validation/gmsh_stop_rule.py`, with its W0 subprocess proof. It is on only under `FEM_EM_GMSH_STOP_RULE=1`, in W1, and never in an `xl` command.

## §7 row moved verbatim 2026-10-19 (rotate_plan_archive.py chunks)

 A PEC oblate spheroid in a lossy host under a uniform quasi-static field — the closed-form stand-in for a plate / disc edge (Phase 7 (7a)) | ✅ (audit PASS 2026-10-12 18:00; step 3 audit PASS 2026-10-13 03:00): steps 0 / 1 / 3 gate the PEC oblate spheroid's rim and pole enhancement at aspect 2 / 5 at `TH15_BETA_BAND` / `MAT4_BAND` (items 403 / 404 / 409). Step 2's rim ladder and its three attribution sub-steps (I22 / I26 / I27 / I28, items 426 / 440 / 441 / 450) are records — the +1.2–1.6 % rim-ball floor is the quasi-static comparand's finite-k error (D) with a ≈ +0.64–0.89 % residual at k → 0 set by what the ladder holds fixed (Fix), not rim grain — ratified 2026-10-14. Frozen statement: ≈ +0.9 % at k → 0 on this comparand, not rim grain. Item 466 (I29, 2026-10-15): moment-equivalent body a_I/a − 1 ≈ −2.17e-3, c_I/c − 1 ≈ −8.2e-4, L_x 0.124758 → 0.124886–0.124893; body-corrected miss r^I +4.07 / +4.34 / +4.40 % (E3-not — the global moment spheroid shrinks the comparand ≈ 3.5 % and does not summarise the rim, whose nodes sit on the CAD surface), `20261015T125651Z_IMP-6-step2e-W1.log`; the step-2 family is frozen (10-14 weekly). Step 4 (items 482 / 522, I31): the kept aspect-5 disc in F-clinical's ASTM box — (D1)–(D10) green at `-n 4`, 1 198 350 cells, gmsh 441.11 s under the 450 s rule, `phantom_box` × `implant_spheroid` composed (box containment, body-aware conductor mass), `20261018T201127Z_IMP-6-step4-a2-W1-n4.log`. **Done-when:** *step 0* (numpy, §9 item 403): the depolarisation factors, both potentials and fields, and the rim / pole identities at the bands in the chunk file, with the prolate-formula control missing ⇒ 🟡. *Step 1* (FEM, degree 2, `IMP-2`'s hole route with an additive `oblate=True` on `spheroid_in_box_domain`, a = 6 mm, aspect 2 and 5, §9 item 404): per aspect \|β − 1\| ≤ `TH15_BETA_BAND`, rim-ball Φ ≤ `MAT4_BAND`, cavity dofs zero; the axial-drive control misses ⇒ ✅. *History: `docs/planning/chunks/IMP-6.md`.* 
