# IMP-4 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record. Opened by the 2026-10-03 03:00 daily review under §9 rule (5), lane 3, inside Phase 7 subgoal (7a) (canonical conductors in a lossy host against closed forms); the weekly ratifies or reverts the opening.

**Citation, corrected.** The 2026-10-02 18:00 review named this candidate "the `IMP-1`-style finite cylinder (Jin §11.3)" and the eqs (11.56)–(11.63). In the repo's copy those equations sit in **§11.2.1, Two-Dimensional Scattering** (`docs/references/jin-fem-3e/11_fem-eigenfunction-expansion-methods.md:355–419`); §11.3 (`:575`) is the unimoment method.
- What §11.2.1 contributes: eq (11.56), the outgoing expansion in H⁽²⁾ₙ(kρ)e^{jnφ} (stated in prose at `:383`), and **Exercise 11.4** (`:461`): plane-wave scattering by a circular conducting cylinder, for both polarisations, compared with Harrington.
- The Bessel identities are Jin App. E: Bessel's equation (E.1) and the Wronskian (E.19) (`appE_special-functions.md:15, 69–71`).
- The equation bodies in the repo's markdown are image placeholders (`[eq (11.56)]`), so only their numbers can be cited. The series itself is Harrington's (*Time-Harmonic Electromagnetic Fields*, §5-8), which the repo does not hold. That is the same footing as Carslaw & Jaeger for `MAT-9`…`MAT-18`: the docstring derives every form.
- "Finite-k" means a finite wavenumber, not the quasi-static limit. The cylinder is **infinite** and the problem is 2D-invariant. There is no closed form for a finite-length rod; that is `IMP-2`'s spheroid and `IMP-3`'s rod.

### §9 item 202 — 2026-10-03 (lane 3)

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.*

**`IMP-4` step 0: the PEC circular cylinder in `IMP-1`'s lossy host at finite k, in closed form.** It is numpy only, verified by residuals and by the current identity before any FEM.
- **Executor:** implementer. **Tier:** standard, `-n 2`, real build. **Class:** lane-3 closed form (measured median 8 min).
- **Cost:** ≈ 8 measured / ≈ 40 predicted min (15 fixed + ≈ 20 code + one ≈ 1-min window + spare).
- **Depends on:** nothing. Shared file: item 199 also appends to `utils/analytical.py`. This item is **append-only** after item 199's `MAT-18` block; rebase if 199 lands first.

**Not already covered.** `IMP-1` (the sphere, quasi-static, then the exact l = 1 field) and `IMP-2` (the spheroid, quasi-static) have no 2D canonical conductor and no lead-current quantity. A 2D PEC wire under an axial field has **no quasi-static limit**: the n = 0 coefficient carries ln(ka). Finite k is intrinsic here, not a correction.

**Fixture.**
- **Host:** `IMP-1`'s — σ 0.57 S/m (`SIGMA_HIGH`), εᵣ 78, 64 MHz, k = `K_FIXTURE` (`test_imp1_step3_exact_comparand.py:51–53`). Principal root ≈ 15.18 − 9.49j 1/m, \|k\| = 17.90 1/m (`IMP-2`'s W0 print, `docs/planning/chunks/IMP-2.md`).
- **Cylinder:** radius a = 1 mm (`GEO-36`'s rod radius). Outer radius for step 1, c = 10 mm, so \|k\|c = 0.179, `MAT-4`'s gated reach and `IMP-1`'s box half-width.
- Convention e^{+jωt}, with H⁽²⁾ as outgoing (Im k < 0, so it decays outward).

**Closed forms** (single-mode drives; the docstring derives each):
- **TM_z, n = 0 (axial E):** E_z = E₀[J₀(kρ) − J₀(ka)H₀⁽²⁾(kρ)/H₀⁽²⁾(ka)], with H_φ = ∂_ρE_z/(jωμ₀) from Faraday.
- **TE_z, n = 1 (transverse E):** H_z = H₀[J₁(kρ) − J₁′(ka)H₁⁽²⁾(kρ)/H₁⁽²⁾′(ka)] cos φ, with H₀ = −2jωε_cE₀/k so that the incident field on the axis is E₀ŷ. E = ∇ × (H_zẑ)/(jωε_c). The PEC condition is E_φ(a) = 0 ⇔ ∂_ρH_z(a) = 0.

**Code** (additive in `utils/analytical.py`; scipy imported inside the functions):
- `pec_cylinder_tm0_field(points, *, k, radius, e0=1.0, interior="zero")` and `pec_cylinder_te1_field(points, *, k, radius, omega, e0=1.0, interior="zero")`.
  - Both take points of shape (N, 3) and return complex (N, 3) Cartesian E.
  - `interior="zero"` returns 0 for ρ < a, as `IMP-1`'s trace does for the chord points. `"continue"` evaluates the formula there, as `pec_prolate_spheroid_field` does (`analytical.py:1576`).
- `pec_cylinder_tm0_current(*, k, radius, omega, e0=1.0)`, which returns 4E₀/(ωμ₀H₀⁽²⁾(ka)).
- The test, `tests/validation/test_imp4_pec_cylinder_closed_form.py`.
- `scipy.special.jv`, `hankel2`, `jvp` and `h2vp` all accept complex arguments; check that in the test rather than assume it.

**Anchors (asserted).**
- **(B1) Bessel's equation** (Jin (E.1)) for both radial functions, with f′ and f″ from scipy's recurrence derivatives (`jvp`/`h2vp`, n = 1, 2), at 64 radii in [a, c]. The residual is bounded by ≤ 1e-12 of the **sum of the term magnitudes** \|f″\| + \|f′/ρ\| + \|(k² − n²/ρ²)f\|, the scale before cancellation (item 190's lesson).
- **(B2) Wronskian** (Jin (E.19), written for H⁽²⁾): J_nH_n⁽²⁾′ − J_n′H_n⁽²⁾ = −2j/(πz) at z ∈ {ka, kc} for n ∈ {0, 1, 2}, ≤ 1e-13 relative.
- **(P) PEC traces** on the Cartesian output at 64 points on ρ = a: the tangential components φ̂·E and ẑ·E are ≤ 1e-14·E₀ for both drives. This catches a bad Cartesian conversion, not just the radial algebra.
- **(I) Wire current:** 2πa·H_φ(a), from the field's analytic ρ-derivative, equals `pec_cylinder_tm0_current` = 4E₀/(ωμ₀H₀⁽²⁾(ka)) to ≤ 1e-12 relative. This is the Wronskian applied to the axial-drive current. Predicted \|I\| ≈ 2.9e-3 A per V/m, from the small-argument H₀⁽²⁾(ka) ≈ 0.64 + 2.64j.
- **(Q) TE static limit:** at k → s·k with s = 1e-6 (arg k kept), on 200 points with a ≤ ρ ≤ c, the field equals the 2D PEC-cylinder static field E_ρ = E₀(1 + a²/ρ²) sin φ, E_φ = E₀(1 − a²/ρ²) cos φ (potential −E₀(ρ − a²/ρ) sin φ) to ≤ 1e-9·E₀. Predicted ≲ 1e-11: the correction goes as (s\|k\|c)² ln.
- **(L) Layout:** a dolfinx-style adapter on x of shape (3, N) returns the (N, 3) function's output transposed, bitwise. Step 1 passes it as `dirichlet_e_field`; this is `IMP-1` step 3(b)'s (L) lesson.

**Printed (predicted, never asserted).**
- k, \|k\|a, \|k\|c.
- The TM n = 0 scattering coefficient J₀(ka)/H₀⁽²⁾(ka) (predicted \|·\| ≈ 0.37).
- \|E_z\|/E₀ at ρ = 1.2a, 1.5a, 2a: the axial field the wire suppresses.
- The TE surface field \|E_ρ(a, π/2)\|/E₀, predicted ≈ 2.00, the static doubling; so the side-of-wire SAR is ≈ 4× the far-field SAR.
- \|I\| and arg I.

**Negative controls.** Rule (e): all are predicted, since nothing has been measured on this fixture. Each is asserted only as missing its anchor's band, with the measured value printed beside the prediction.
- **(c1) H⁽¹⁾ in place of H⁽²⁾** in the TM field, the wrong outgoing convention for e^{+jωt}, against (I). The current becomes −4E₀/(ωμ₀H₀⁽¹⁾(ka)), so the miss is predicted ≈ 0.7 from the small-argument forms. Asserted > 1e-3. (B1) and (P) still pass, which is why (I) is the discriminating anchor.
- **(c2) The Dirichlet condition in TE** (H_z(a) = 0 instead of ∂_ρH_z(a) = 0), against (P). The tangential E_φ(a) is predicted O(E₀). Asserted > 1e-3·E₀.
- **(c3) σ dropped** (real k = k₀√εᵣ) in the current, against (I) evaluated at `K_FIXTURE`. Predicted miss ≈ 0.15. Asserted > 1e-3.

**Window.** Foreground, Bash `timeout: 660000`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Predicted 45–60 s: `IMP-1` step 3(a) took 45 s (`20261001T142052Z_IMP-1-step3a-W1.log:263`) and `IMP-2` step 0 W1b 45 s (`20261001T142617Z_IMP-2-step0-W1b.log:308`), both numpy windows in the real build.

`scripts/testing/run_and_log.sh IMP-4-step0-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/environment tests/unit/test_analytical_lightweight.py tests/validation/test_imp4_pec_cylinder_closed_form.py -s -v --tb=short'"`

The lightweight test is the rule (c) re-gate of the additive `analytical.py` change.

**Traps.**
- **Convention:** e^{+jωt} throughout, H⁽²⁾ outgoing. `K_FIXTURE` is the principal root (Re > 0, Im < 0, `:53`); never `np.sqrt` of a negated ε_c.
- **Derivatives:** use analytic or recurrence derivatives. An FD check on these functions would need δ ≳ 1e-2·a and the Richardson or ratio form (`IMP-2`'s (F) ruling); none is pre-registered here.
- **The ρ = 0 axis** is outside [a, c]. The TE Cartesian conversion still divides by ρ, so guard ρ = 0 and never sample it.
- **Normalise (B1) by the parts before cancellation,** never by \|f\|: J₁ and H₁ are small or large at small kρ.
- **Run hygiene:** print on rank 0; the real build skips the four complex-only environment tests, as expected; `-s`; rule (i) — no edit after the last window.
- **Rule (f)** applies to a parameter (s, the radii), never to a bound.
- **Shared file:** append at the end of `analytical.py`, and rebase if item 199 lands first.

**Scope.** Numpy only, two single-mode drives on an infinite cylinder. **Not included:**
- the FEM (step 1);
- the plane-wave multi-mode series;
- a finite rod or its tip (`IMP-2`, `IMP-3`);
- a coil, heating, or any published-data comparison.

**Status it can move.** The proposed §7 `IMP-4` row moves ⬜ → 🟡 on its step-0 Done-when, which this item matches (rule (j)). The review adds the row in the commit that queues this item, and names it in the dashboard for the weekly to ratify (rule (5) lane 3). Step 1 then becomes draftable.

**Negative result.**
- **(B1), (B2), (P) or (L) red:** a code bug. Fix it in-slot, run W1b, and disclose.
- **(I) red with (B2) green:** a sign or convention error in H_φ. Print both sides and their ratio, then stop.
- **(Q) red:** print the miss at s ∈ {1e-4, 1e-5, 1e-6}. If it scales as s² ln s, s is a parameter: re-register it under rule (f), disclosed; otherwise stop.
- **On a stop:** park on `attempt/IMP-4-step0-<UTC>` and mark BLOCKED in the same commit (rule (d)); the row stays ⬜.
- There is no PREDICTED RED line.

**Outcome — DONE 2026-10-03 (12:00 CDT slot, second item, implementer).** W1 `20261003T172120Z_IMP-4-step0-W1.log` passed first time: `-n 2`, real build, 55 s, 25 passed, 4 skipped (complex-only environment tests). No W1b.

Anchors:
- **(B1)** residual over term sum 1.8e-16 – 2.6e-16 (J₀, H₀, J₁, H₁, both radial functions).
- **(B2)** 4.467e-16. **(P)** TM φ̂ 0 / ẑ 2.235e-16; TE φ̂ 1.142e-16 / ẑ 0 (·E₀).
- **(I)** 2.601e-16; I = 6.92952e-4 − 2.83487e-3j A per V/m.
- **(Q)** 9.302e-15·E₀ at s = 1e-6 (predicted ≲ 1e-11). **(L)** bitwise, both drives.

Controls: (c1) 0.6749 (predicted ≈ 0.7; the H⁽¹⁾ field still passes B1 2.3e-16 and P 2.3e-16); (c2) 1.999·E₀ (predicted O(1)); (c3) 0.1443 (predicted ≈ 0.15).

Printed: k = 15.178 − 9.48858j 1/m, \|k\|a 0.0179, \|k\|c 0.179; J₀(ka)/H₀⁽²⁾(ka) magnitude 0.3687 (H₀⁽²⁾(ka) = 0.644055 + 2.63483j); \|E_z\|/E₀ 0.0428 / 0.0952 / 0.1627 at 1.2a / 1.5a / 2a; TE surface \|E_ρ(a, π/2)\|/E₀ 2.00022; \|I\| 2.91833e-3, arg −1.33106 rad.

Disclosed:
- Private control keywords (`_hankel_kind`, `_condition` on the public functions; `hankel_kind` / `condition` on the radial helpers); defaults are the brief's forms.
- (P) evaluates with `interior="continue"`, because `hypot` can return ρ a hair below a, and `"zero"` would then pass vacuously. A guard asserts that the normal TE component on that ring is > E₀.
- The TE field is independent of `omega` (μ = μ₀ host; −H₀/(jωε_c) = 2E₀/k); the ρ = 0 axis is guarded.
- An extra `test_scipy_accepts_complex_arguments` checks `jv` against the power series for J₀ and J₁ to 1e-14. That bound is the executor's choice and was not pre-registered.
- Fixture constants are imported from `test_imp1_step3_exact_comparand.py`.

Row ⬜ → 🟡; step 1 (item 206) is unblocked.

**What step 1 would be (not briefed; for the row).**
- **Mesh:** `TH-15`'s hole route, degree 2. An additive `as_hole=True` on `MeshGenerator.cylindrical_domain` (`io/mesh.py:833`, today a fragment): the inner cylinder is cut, the cavity wall gets facet tag 2 from the **meshed** boundary (`TH-15`'s lesson, `:2335–2340`), the outer lateral wall tag 1, the end annuli a new tag 3, with graded cavity sizing as in `IMP-2`. `False` stays bit-identical, re-gated by `tests/mesh/test_cylindrical_domain.py` and `tests/solver/test_cylinder.py`.
- **Boundary conditions:** `TimeHarmonicProblem(boundary_condition="pec_zero_tangential_a", dirichlet_e_field=…, pec_facet_tags=…)`, `IMP-1`'s `_solve` pattern (`test_imp1_pec_sphere_lossy_host.py:198–257`).
  - TM constrains tags {1, 2, 3}: the exact z-invariant field has zero tangential E on the end planes, so PEC there is exact.
  - TE constrains {1, 2} and leaves the ends natural, which is PMC (`core/time_harmonic.py:345–354`). That is exact because n × curl E = 0 there.
- **Gates at each drive:**
  - the complex LSQ amplitude of the closed-form scattered field on probe circles at ρ = 1.2a, 1.5a, 2a in the mid-plane: \|β − 1\| ≤ `TH15_BETA_BAND`;
  - the annulus-mean \|E\|² (a < ρ < 1.5a) against the closed form on the same discrete mask: ≤ `MAT4_BAND`;
  - the cavity dofs are zero;
  - control: the natural cavity misses.
- **Price:** `IMP-1` at degree 2 ran 54 524 cells in 101.8 s (`20261001T033100Z_IMP-1-land-W1.log:446`), so ≈ 20–40k cells per solve here and a heavy window.

### §9 item 206 — 2026-10-03 (lane 3)

*Drafted by a general-purpose agent for the 2026-10-03 daily review; the review reviews it.*

**`IMP-4` step 1: the FEM half. The PEC cylinder is cut from a coaxial host annulus, the outer wall carries step 0's exact finite-k trace, and the end planes carry the condition that makes the 3D solve z-invariant for each drive.**
- **Executor:** implementer. **Tier:** heavy, `-n 2`; complex build for W1–W2, real build for W0. **Class:** lane-3 closed form (measured median 8 min).
- **Cost:** ≈ 8 measured / ≈ 65 predicted min (15 fixed + ≈ 25 code + W0 ≈ 4 + W1 ≈ 6 + W2 ≈ 5 + ≈ 10 spare).
- **Depends on item 202 (step 0's closed form); if it did not land, skip.** The check is `grep -n 'def pec_cylinder_tm0_field' src/fem_em_solver/utils/analytical.py` non-empty, plus 202's W1 green in `test-results.md`. This item **does not edit `utils/analytical.py`**: it imports 202's three functions, and the incident-mode split is test-local. Should an edit become necessary, the file is append-only after items 199 / 201 / 202's blocks; rebase if any lands first.

**Why this fixture is exact, not approximate.** The z-invariant reduction is Jin §1.4.2 (eqs (1.44), (1.45), `01_basic-electromagnetic-theory.md:265–277`) and §11.2.1 eq (11.53) (`11_fem-eigenfunction-expansion-methods.md:361`): ϕ = E_z for E_z-polarisation, ϕ = H_z for H_z-polarisation. A finite slab of the infinite problem reproduces the infinite solution exactly when every boundary condition is one that the exact z-invariant field already satisfies. With the solve well posed (lossy host), the FEM's continuous solution is then the closed form, and only discretisation and faceting remain.
- **TM_z, n = 0** (E = E_zẑ, H = H_φφ̂). On a plane z = const, n̂ = ±ẑ:
  - n̂ × E = 0 identically, so the **PEC** end condition (Jin §1.5.2 eq (1.52), `01:315–321`; the Dirichlet "electric wall" n̂ × E = 0, `08_vector-finite-elements.md:1533`) is exact.
  - n̂ × H = ∓H_φ ẑ × φ̂ ≠ 0, so the natural end (the Neumann "magnetic wall" n̂ × (∇ × E) = 0, `08:1533`; `core/time_harmonic.py:345–348` "NATURAL … is PMC") is **wrong**: it forces H_φ = 0 on the ends and breaks z-invariance.
  - The ends are therefore tag 3 **constrained**: pec tags (1, 2, 3).
- **TE_z, n = 1** (H = H_zẑ, E transverse). On the end planes:
  - n̂ × E = (E_x, E_y) rotated ≠ 0, so a PEC end is **wrong**.
  - n̂ × H = H_z ẑ × ẑ = 0, so the natural (PMC) end is exact.
  - The ends are tag 3 **left out**: pec tags (1, 2).
- **This is the step-0 brief's assignment** (`:92–93`), and the §7 row's ("end planes PEC (TM) or natural/PMC (TE)"). The opposite pairing (PMC ends for E_z, PEC for TE) is the trap, not the fixture.
- **Truncation.** The outer wall ρ = c carries the exact total field's tangential trace from step 0, so it is exact at finite k. Jin's own §11.2.1 route, the global third-kind condition (11.61)/(11.63) from the outgoing expansion (11.56), is **not** implemented or tested here. Imposing the exact trace is the `IMP-1` step 3(b) idiom, and it removes truncation from the error budget.
- **The sizing rule does not bind.** The rule "(\|k\|c)² under a third of the band" (`IMP-2.md:9`) bounded a quasi-static comparand. Here comparand and trace are exact at finite k. c = 10 mm keeps \|k\|c = 0.179, `IMP-1`'s and `MAT-4`'s reach, for continuity only.

**Fixture.**
- **Host:** step 0's, which is `IMP-1`'s. Import `FREQUENCY_HZ`, `SIGMA_HIGH`, `EPSILON_R_HOST`, `E0`, `DEGREE` = 2, `OUTER_BOUNDARY_TAG` = 1, `CAVITY_TAG` = 2, `TH15_BETA_BAND`, `MAT4_BAND`, `_assemble_real` and `_cavity_dofs` from `test_imp1_pec_sphere_lossy_host.py`. Import `K_FIXTURE` from `test_imp1_step3_exact_comparand.py:53`. Import `_owned_masked_cells` and `_numpy_quadrature_mean` from `test_imp2_spheroid_fem.py:172–201`. Never re-derive these.
- **Geometry:** a = 1 mm, c = 10 mm, length L = 4 mm, centred on z = 0.
- **Permittivity:** ε_c = K_FIXTURE²/(ω²μ₀) wherever the TE normalisation needs it. That keeps the incident split consistent with the k step 0 used; never rebuild ε_c from σ and εᵣ separately.

**Code.**
- **Mesher (additive; rule (c), disclosed).** Add keyword-only `*, as_hole=False, resolution_wall=None, resolution_far=None` to `MeshGenerator.cylindrical_domain` (`io/mesh.py:833`). `as_hole=False` takes the existing branch untouched; passing either resolution without `as_hole` raises, so nothing is silently ignored. Under `as_hole=True`, modelled on `sphere_in_box_domain`'s branch (`:2403–2436`):
  - `occ.cut` the inner cylinder from the outer; exactly one volume remains (raise otherwise). Volume tag 2.
  - Facet groups come from `getBoundary` of the **meshed** air volume (`TH-15`'s rule, `:2338–2341`), classified **z first**. Faces with z-extent < tol at \|z\| = L/2 are tag 3. Then bbox r_max within tol of c is tag 1. The rest is the cavity, tag 2. Raise if any group is empty, and print the face count per group.
  - Use tol = `_WALL_TOL_FRACTION`·(c − a) = 9e-5 m, ≥ 10× OCC's 1e-7 bbox padding (`IMP-2`'s W0 failure `20261001T200341Z`).
  - Sizing: `Distance` from the cavity faces → `Threshold` (SizeMin `resolution_wall` = a/6 at DistMin 0, SizeMax `resolution_far` = c/6 at DistMax c − a), set as the background mesh. This is h(ρ) = ρ/6 exactly, the radial grading the ln and 1/ρ² fields want. Turn off `MeshSizeExtendFromBoundary`, `MeshSizeFromPoints` and `MeshSizeFromCurvature` as `sphere_in_box_domain` does. `generate(3)` only, with no Netgen optimise on this branch. Print the minimum SICN.
- **Module** `tests/validation/test_imp4_pec_cylinder_fem.py`:
  - **Trace:** `lambda x: np.ascontiguousarray(f(x.T).T).astype(PETSc.ScalarType)` around `pec_cylinder_tm0_field` / `pec_cylinder_te1_field` with `interior="zero"`. This is `IMP-2`'s `_dirichlet_axial` pattern and step 0's (L).
  - **Incident split (test-local, numpy):**
    - TM: E_inc = E₀J₀(kρ)ẑ.
    - TE: H_inc = H₀J₁(kρ) cos φ with H₀ = −2jωε_cE₀/k, and E_inc = ∇ × (H_incẑ)/(jωε_c) in Cartesian form via `jv` / `jvp`.
    - The scattered basis is E_s = E_cf − E_inc.
  - **Probes:** 24 points on each of ρ = 1.2a, 1.5a and 2a at z = 0, angles (i + ½)·2π/24, evaluated with `evaluate_vector_field_parallel` on `e_real` and `e_imag`; raise unless every point is valid. β is the **complex** LSQ over all components, β = Σ conj(E_s)·(E_h − E_inc) / Σ \|E_s\|². Do not use `IMP-1`'s Re-only fit: arg k = −32°, and the 18:00 pathologist found the Re-only fit absorbs h² error (`IMP-1.md:33`).
  - **Annulus mask:** a DG0 indicator of cell midpoints with x² + y² < (1.5a)², numpy only. Φ_h = ⟨\|E_h\|²⟩/E₀² by UFL at `quadrature_degree` 8. Φ_cf,h comes from `_numpy_quadrature_mean` over the same owned cells, using the closed form with **`interior="continue"`**: quadrature points in the chord sliver have ρ < a, and the "zero" mask would read them as zero (`IMP-1`'s trap, item 146).
  - **Solve:** a module-local `_solve_cylinder(drive, pec_tags)` following `test_imp2_spheroid_fem.py:213–295`, using `boundary_condition="pec_zero_tangential_a"`, `dirichlet_e_field`, `pec_facet_tags=pec_tags`. It returns β, Φ_h, Φ_cf,h, the cavity dof count and max \|E\|, cells, seconds and summed RSS.
  - **Tests:** `test_census_and_numpy` (W0, **not** `complex_only`: mesh, numpy and real UFL only), `test_both_drives` (W1, `complex_only`), `test_natural_cavity_controls` (W2, `complex_only`).

**Anchors (asserted).**
- **W0 (census, before any solve):**
  - (Guard) global cells ≤ **60 000**. Predicted 20–45 k: h = ρ/6 gives 6·216·2πL(1/a − 1/c) ≈ 29 k, calibrated on `IMP-1` a/6's ≈ 6.1 tets per h³ (13 467 cells, `20261001T033100Z_IMP-1-land-W1.log:370`).
  - (Tag) reduced `ds(t)` areas (pinned `quadrature_degree` 2): tag 1 against 2πcL, tag 2 against 2πaL, tag 3 against 2π(c² − a²), each within 2 %. The faceting deficit is predicted ≤ 0.5 % ((Δφ)²/6 at Δφ = 1/6). A mis-tag is O(1): annuli tagged 1 would read +250 %.
  - (M′) on the annulus mask, the numpy quadrature of the 2D static PEC-cylinder TE field E = (2E₀a²xy/ρ⁴, E₀ + E₀a²(y² − x²)/ρ⁴, 0) equals the same field's UFL cell integral at degree 8 to ≤ 1e-12 relative. The field is a polynomial ratio, so it needs no `sqrt` and no comparison. The continuous value 1 + ⟨a⁴/ρ⁴⟩ = 1.4444 is printed beside it.
  - (S) the incident split: TM (E_cf,z − E₀J₀(kρ))/H₀⁽²⁾(kρ) is constant and equals −E₀J₀(ka)/H₀⁽²⁾(ka); TE E_cf − E_inc equals E_inc's formula with J₁ → H₁⁽²⁾ (and J₁′ → H₁⁽²⁾′) times −J₁′(ka)/H₁⁽²⁾′(ka). Each holds to ≤ 1e-12 relative on 64 points in [1.2a, 2a].
- **W1, per drive** (TM pec tags (1, 2, 3); TE pec tags (1, 2)). These are the §7 row's step-1 Done-when, verbatim in substance:
  - (i) \|β − 1\| ≤ `TH15_BETA_BAND` (4.886 %);
  - (ii) \|Φ_h/Φ_cf,h − 1\| ≤ `MAT4_BAND` (10 %);
  - (iii) cavity dofs located > 0 (reduced), with max \|E\| < 1e-12 (reduced with `MPI.MAX`).

**Printed (predicted, never asserted; rule (e)).**
- **TM, faceting-led.** Faceting puts the effective wall at a(1 − (Δφ)²/12) ≈ a(1 − 0.23 %).
  - Φ_cf,h is ≈ 3.6e-3: \|E_z\|/E₀ ≈ 0.235 ln(ρ/a) by the small-argument forms, and ⟨ln²⟩ = 0.0661 over the annulus. Note that this is **small**, so it is the fragile gate.
  - The miss is predicted **positive, +1 to +3 %**: 2·0.23 %·⟨ln⟩/⟨ln²⟩ ≈ 1.7 %, plus discretisation.
  - \|β − 1\| is predicted ≤ 1.5 %; faceting alone gives ≈ 0.235·0.23 % ≈ 0.05 %.
- **TE.**
  - β is predicted ≈ 0.995, since the 2D dipole goes as a_eff². \|β − 1\| ≤ 1.5 %.
  - Φ_cf,h ≈ 1.444, with a miss of ≲ 1 %.
- **For each drive:** Φ_cf,h beside the continuous annulus value (a numpy radial / angular quadrature of step 0's field), which isolates faceting; cells, degree-2 dofs, seconds and summed RSS; h_wall/a; min SICN.
- **z-invariance:** max\|E_h(z = ±L/4) − E_h(0)\| / max\|E_h(0)\| on the 1.5a circle, predicted ≲ 1e-3. A large value means an end condition is wrong.
- **TM current:** β·`pec_cylinder_tm0_current` (\|I\| ≈ 2.9e-3 A per V/m from step 0), the induced current the fitted amplitude implies.

**Negative controls** (rule (e): nothing has been measured on this fixture, so each is predicted and asserted only as missing the band).
- **(c1) TM, natural cavity** (pec tags (1, 3)), which makes the wire PMC: H_φ(a) = 0.
  - **Prediction.** The scattered coefficient J₁(ka)/H₁⁽²⁾(ka) is O((ka)²) ≈ 1e-3. The annulus therefore carries the regular mode pinned to the trace, E_z ≈ C₁J₀(kρ) with C₁ = E_cf(c)/J₀(kc) ≈ 0.525 + 0.128j (\|C₁\| 0.54). The LSQ β is then ≈ 0.53 − 0.13j, so \|β − 1\| ≈ **0.49 (range [0.3, 0.7])**, and Φ_h ≈ \|C₁\|² ≈ 0.29 against Φ_cf,h ≈ 3.6e-3 gives a miss ≈ +7 900 % (range [+3 000, +15 000] %).
  - **Ceilings.** On β, 0.49/0.04886 ≈ **10×** the band. On Φ there is no saturation, because the comparand is small.
  - **Asserted:** \|β − 1\| > `TH15_BETA_BAND` and \|Φ_h/Φ_cf,h − 1\| > `MAT4_BAND`.
- **(c2) TE, natural cavity** (pec tags (1,)), which makes H_z(a) = 0, the void cylinder.
  - **Prediction.** In the static annulus, ψ = A(ρ + a²/ρ) sin φ with ∂_ρψ(a) = 0 and ψ(c) on the PEC trace gives A = −E₀(1 − (a/c)²)/(1 + (a/c)²). Hence β ≈ **−0.980**, and \|β − 1\| ≈ 1.98, with a ceiling of 1.98/0.04886 ≈ **40×**. The uniform residual (A + E₀) is orthogonal to the dipole on full circles, so it does not leak into β.
  - **Calibration of the method.** The same annulus argument for `IMP-1`'s 3D control, A = −(1 − (a/b)³)/(1 + (a/b)³/2), gives −0.41 against the measured −0.433 / −0.439 (`IMP-1.md:94`; `20261001T033747Z_IMP-1-land-W2.log:152`), on a box rather than a sphere.
  - **Asserted:** \|β − 1\| > `TH15_BETA_BAND` **only**.
  - **Φ cannot discriminate in TE.** A 2D line dipole of either sign carries the same annulus-mean \|E\|², A²(1 + ⟨q²⟩) against 1 + ⟨q²⟩. The predicted miss is 1 − A² ≈ **3.9 %**, inside the band. It is printed beside its prediction, and not asserted either way. **For the review:** the row's "natural-cavity control misses" is met on β for TE, and on β and Φ for TM.

**Windows** (foreground, Bash `timeout: 660000`; each harness call is the whole Bash command, bare and repo-relative; `-s`; no `-k`; no pipes). `PYTHONPATH` carries `:/workspace` because the module imports `tests.validation.*`, as `IMP-1` / `IMP-2`'s windows do.
- **W0** (real build). Predicted 2–5 min; the two re-gate modules are unmeasured as a pair, the rest is `IMP-2` W0's 109 s (`20261001T201413Z_IMP-2-step1-W0.log`). This window is the rule (c) re-gate of `as_hole=False` (bit-identical branch), step 0's module re-run as the dependency check, and the census:
  `scripts/testing/run_and_log.sh IMP-4-step1-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_cylindrical_domain.py tests/solver/test_cylinder.py tests/validation/test_imp4_pec_cylinder_closed_form.py tests/validation/test_imp4_pec_cylinder_fem.py::test_census_and_numpy -v -s --tb=short'"`
- **W1** (complex). Predicted 250–400 s: two solves at 20–45 k cells, against `IMP-2`'s 105.8 s / 83.3 s per solve at 37 359 / 44 304 cells and 6.46 GiB summed (`20261001T201608Z_IMP-2-step1-W1.log:381,449`), where that W1 took 307 s.
  `scripts/testing/run_and_log.sh IMP-4-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp4_pec_cylinder_fem.py::test_both_drives -v -s --tb=short'"`
- **W2** (complex). Predicted 250–350 s: two solves, against `IMP-2` W2's 179 s for one at 76.8 s (`20261001T202123Z_IMP-2-step1-W2.log:151`).
  `scripts/testing/run_and_log.sh IMP-4-step1-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp4_pec_cylinder_fem.py::test_natural_cavity_controls -v -s --tb=short'"`

**Traps (all paid for).**
- **End conditions are polarisation-specific** (above). Swapping them is an O(1) wrong field that the probes would read as a β miss, not a small error.
- **Mesh classification.** `cylindrical_domain`'s existing bbox-r_max classification would put the end annuli (r_max = c) into tag 1, so classify z first. The cavity group comes from the meshed boundary, never the cut tool (`TH-15`).
- **Masks and fields.** Use the "zero" interior for the trace and "continue" for Φ_cf,h. Guard N = 0 on a rank inside every callable. Never sample ρ = 0 (the TE conversion divides by ρ).
- **Complex build:** no UFL comparison and no `ufl.max_value`. Real operands do not save a comparison containing `sqrt` (`OPS-22`, `WF-6` 4c). Every mask is numpy.
- **Quadrature:** every `SpatialCoordinate`-bearing form carries `metadata={"quadrature_degree": …}` (`POST-5` step 1); the `ds` area forms are pinned too.
- **Rank safety.**
  - `assemble_scalar`, `cell_tags.values` and local max/min are rank-local; reduce before asserting.
  - Do the RSS `allreduce` **outside** any rank-0 branch: `IMP-2`'s W0 failure `20261001T200855Z`.
  - Owned cells only (`size_local`). Read the `geometry.dofmap` shape at runtime.
- **JIT cache.** If an earlier window today was killed, sweep orphan `.c` stems first, and only when `pgrep -c python3` in the container reads 0, with the rubric's one bare harness call (`<ID>` = `IMP-4-step1`). Never clear the cache wholesale.
- **Run hygiene.** Foreground; container `timeout -k 30`; never pipe or filter inside the harness command. Print on rank 0. Rule (i): no edit to a file after its last window.
- **Rule (f)** reaches a parameter (a resolution, L), never a bound.

**Scope.** One mesh, degree 2, two single-mode drives on an infinite PEC cylinder with an exact trace. These are accuracy-level bands, **not a convergence claim**. **Not included:**
- Jin's (11.63) global boundary condition;
- the plane-wave multi-mode series of Exercise 11.4 proper;
- modes n ≥ 2;
- a finite rod or tip (`IMP-2`, `IMP-3`);
- a coil, heating, 1 g SAR, or published data.

**Status it can move.** The `IMP-4` §7 row, whose Done-when reads: "*Step 1* (FEM, degree 2, `TH-15` hole route): additive `as_hole=True` on `cylindrical_domain`; exact trace on the outer wall; end planes PEC (TM) or natural/PMC (TE); per drive \|β − 1\| ≤ `TH15_BETA_BAND`, annulus ⟨\|E\|²⟩ on the same mask ≤ `MAT4_BAND`, cavity dofs zero, natural-cavity control misses ⇒ ✅".
- With step 0 green (🟡), this item green moves **🟡 → ✅**, auditable. The opening itself is the weekly's to ratify (rule (5) lane 3), as `IMP-1`'s was.
- The review may then add `IMP-4` to §2's closed-form-lane sentence: imposed trace, no coil.
- **Flag for the review (rule (j)):** the control clause is read as β for TE, because Φ cannot miss there (see (c2)).

**Negative result.**
- **Over the cell guard in W0:** stop before any solve, report the counts, BLOCKED.
- **(Tag), (S) or (M′) red:** a mesher or helper bug. Fix it in-slot, run W0b, and disclose.
- **A re-gate red on `as_hole=False`:** the branch is not bit-identical. Fix it; if it was red at `HEAD` already, cite known-issues and continue.
- **(i) or (ii) red with (iii), (M′) and (Tag) green:** print h_wall/a, Φ_cf,h against the continuous value and the z-invariance number. Park on `attempt/IMP-4-step1-<UTC>` and mark BLOCKED in the same commit (rule (d)); no band moves. If TM (ii) is the red, positive, and near the faceting prediction, the review may re-pose the wall resolution (a/10) under rule (f).
- **A control inside its band:** the cavity condition is not reaching the probes, which points at a tag or solver bug. Stop and report.
- There is no PREDICTED RED line.

**Outcome (2026-10-03 15:00 slot): DONE, row 🟡 → ✅.** Executor `implementer`, foreground; W0 / W1 / W2 green on first run, rule (i) holds (both files predate W0).
- **Windows:** `20261003T200309Z_IMP-4-step1-W0.log` (real, 22 passed / 4 skipped, 96 s; re-gates `test_cylindrical_domain.py`, `test_cylinder.py` and step 0 green); `20261003T200454Z_IMP-4-step1-W1.log` (complex, 13 passed, 123 s); `20261003T200703Z_IMP-4-step1-W2.log` (complex, 13 passed, 117 s).
- **W0:** 25 568 cells (guard 60 000), faces 1 / 1 / 2, min SICN 0.2988, h_wall/a 0.396. (Tag) −0.088 % / −0.088 % / −0.455 % (band 2 %). (M′) 2.0e-15 (numpy 1.443486, continuous 1.4444). (S) TM 3.0e-16, TE 6.0e-16.
- **W1:** TM β 0.999655 − 0.000093j (\|β − 1\| 0.036 %), Φ_h 3.742e-3 vs Φ_cf,h 3.705e-3 (+1.01 %; Φ_cf,h vs continuous +1.77 %, the predicted faceting offset); TE β 0.992258 (0.774 %), Φ −0.23 %. Cavity 10 704 dofs, max\|E\| 3.7e-20 / 2.1e-34. 174 112 dofs, ≈ 37 / 29 s per solve, 4.28 GiB summed. TM current β·I \|I\| 2.917e-3 vs 2.918e-3 A per V/m.
- **W2 controls:** c1 β 0.527 − 0.126j (\|β − 1\| 0.489, predicted 0.49), Φ +7 790 % (predicted ≈ +7 900 %); c2 β −0.972 (\|β − 1\| 1.97, predicted 1.98), Φ −4.0 % (printed, predicted ≈ 3.9 %).
- **Off prediction (printed only):** z-invariance 7.3e-3 (TM) / 1.3e-2 (TE) against ≲ 1e-3. Every gate is tight, so this is not read as a wrong end condition; the executor attributes it to unstructured discretisation differing between z-planes — untested.
- **Disclosed `src/` change (rule (c)):** `cylindrical_domain` gains keyword-only `as_hole` / `resolution_wall` / `resolution_far`; the `as_hole=False` body moved verbatim into `_cylindrical_fragment_domain`; the new branch is the module-level `_cylindrical_hole_domain`.

### Audit — 2026-10-03 18:00 review

**PASS** (`auditor`, closing commit `9b8fefc7`). Evidence: `20261003T200454Z_IMP-4-step1-W1.log:134–135` (TM \|β − 1\| 0.0357 %, Φ +1.0110 %), `:140–141` (TE 0.7742 %, −0.2345 %), `:136,142` (cavity max\|E\| ≤ 3.7e-20), `:294–296` (Status 0, 123 s); controls `…200703Z…W2.log:133–143` (c1 \|β − 1\| 48.92 %, Φ +7 789.79 %; c2 \|β − 1\| 197.23 %); W0 `…200309Z…W0.log:465–470`. Caveats (no demotion): the z-invariance prints 7.3e-3 (TM) / 1.3e-2 (TE) against a pre-registered ≲ 1e-3 — a 7–13× miss on an unasserted print whose plane-to-plane explanation is untested (for the weekly); one mesh at degree 2, no convergence claim; the TE Φ control cannot see a wrong cavity condition (the 10:30 rule (j) ruling), while the β control does; the opening awaits the weekly's ratification. Journal: `docs/planning/reviews/2026-10-03-daily-1800.md` §3.
