# IMP-5 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

**Opened** by the 2026-10-03 18:00 daily review under §9 rule (5), lane 3, on a closed-form licence: the anchor is the exact Bessel/Hankel series of a dielectric-coated PEC circular cylinder. This is Phase 7 subgoal (7a), canonical conductors in a lossy host against closed forms, and the natural successor of `IMP-4` (✅ 2026-10-03, the bare PEC cylinder on the hole route). The weekly ratifies or reverts the opening.

**What it adds over `IMP-4`.** The insulated lead shaft. A lossless dielectric sleeve a < ρ < b around the PEC conductor puts the project's first **conforming dielectric interface** (n̂ × E, n̂ × H and n̂·D continuous) under a closed form at finite k, in the lossy host. In the TE drive the contrast across that interface is \|ε_c,host/ε_c,coat\| ≈ 59.

**Done-when.**
- *Step 0* (numpy closed form, §9 item 210) green ⇒ ⬜ → 🟡.
- *Step 1* (FEM on the `TH-15` hole route with a coating subdomain, §9 item 211) green ⇒ 🟡 → ✅.

**Citation.**
- **Jin 3e does not print the coated-cylinder series.** It uses it without stating it. §10.2.3, Fig. 10.11 (`10_fem-boundary-integral-methods.md:625–631`) compares an FE–BI solution for a coated circular conducting cylinder with "the exact eigenfunction expansion solution". §11.1, Fig. 11.4 (`11_fem-eigenfunction-expansion-methods.md:127–133`) is a dielectric-coated conducting post in a waveguide.
- **What Jin does supply, cited for each piece the docstring derives:**
  - §1.4.2 eqs (1.44)/(1.45) (`01_basic-electromagnetic-theory.md:269–277`): the z-invariant scalar reductions with piecewise εᵣ, μᵣ. They give the scalar interface conditions: TM, E_z and μᵣ⁻¹∂_ρE_z continuous; TE, H_z and ε_c⁻¹∂_ρH_z continuous.
  - §1.5.1 eqs (1.46)–(1.49) (`01:285–299`): the four field-continuity conditions, of which only two are independent.
  - §1.5.2 eq (1.52) (`01:319–321`): PEC, n̂ × E = 0.
  - §11.2.1 eq (11.53) with (u_r, v_r) per polarisation, and eq (11.56), the outgoing H⁽²⁾ₙ expansion (`11:355–383`).
  - App. E: Bessel's equation (E.1) and the Wronskian (E.19) (`appE_special-functions.md:15, 69–71`).
- **The series itself** is the standard layered-cylinder eigenfunction solution. Cite Harrington, *Time-Harmonic Electromagnetic Fields* (1961), ch. 5 (§5-8, scattering by cylinders; the coated case is an extension of that section's method), and Balanis, *Advanced Engineering Electromagnetics*, ch. 11 (scattering by circular cylinders). **Neither book is held in the repo, and neither section was re-verified against its text.** That is the same footing as `IMP-4` (Harrington §5-8) and `MAT-9`…`MAT-18` (Carslaw & Jaeger): the docstring derives every coefficient from the Jin equations above.
- The quasi-static anchor (Q) is derived independently, by 2D potential theory for a coated PEC cylinder in a uniform field. It shares no Bessel algebra with the series.

**Fixture.**
- **Host:** `IMP-4`'s, which is `IMP-1`'s: σ 0.57 S/m, εᵣ 78, 64 MHz, k_host = `K_FIXTURE` = 15.178 − 9.4886j 1/m.
- **Conductor and coating:** a = 1 mm, b = 1.5 mm, c = 10 mm.
  - A real lead conductor (≈ 0.3–0.6 mm radius, ≈ 0.1–0.3 mm insulation) is too fine for `IMP-4`'s mesh. So a keeps `IMP-4`'s radius, and b/a = 1.5 puts ≈ 2–3 degree-2 cells across the sleeve under `IMP-4`'s h(ρ) = ρ/6 sizing (measured h_wall/a 0.396).
  - The thin-insulation limit is carried only by the analytic b → a anchor.
- **Coating:** εᵣ 3, σ 0 (lossless), μ = μ₀. This is order-of-magnitude silicone or polyurethane, not a measured material.
  - k_coat = K0·√3 ≈ 2.3217 1/m (real), with `K0` the same free-space wavenumber `K_FIXTURE` is built from.
  - The ratio ε_c,coat/ε_c,host = k_coat²/k_host² ≈ 7.39e-3 + 1.515e-2j (\|·\| ≈ 0.0169).

### §9 item 210 — 2026-10-03 (lane 3)

*Drafted by a general-purpose agent for the 2026-10-03 18:00 daily review; the review reviews it.*

**`IMP-5` step 0: the dielectric-coated PEC circular cylinder in `IMP-4`'s lossy host at finite k, in closed form.** It is numpy only, verified by residuals, interface identities, two exact limits onto `IMP-4` and an independent quasi-static form, before any FEM.
- **Executor:** implementer. **Tier:** standard, `-n 2`, real build. **Class:** lane-3 closed form (measured median ≈ 12 min).
- **Cost:** ≈ 12 measured / ≈ 45 predicted min (15 fixed + ≈ 25 code + one ≈ 1-min window + spare).
- **Depends on:** nothing. Shared file: `utils/analytical.py` is **append-only** after `IMP-4`'s block (`analytical.py:2230–2376`). Rebase if another item appends first.

**Closed forms** (single-mode drives, e^{+jωt}, H⁽²⁾ outgoing; the docstring derives each). Regions: coating a ≤ ρ < b (k₁ = k_coat), host ρ ≥ b (k₂ = k_host).
- **TM_z, n = 0 (axial E).**
  - Coating: E_z = E₀·c·g(ρ), with g = J₀(k₁ρ)Y₀(k₁a) − Y₀(k₁ρ)J₀(k₁a), so g(a) = 0.
  - Host: E_z = E₀[J₀(k₂ρ) + s H₀⁽²⁾(k₂ρ)].
  - At ρ = b, E_z and ∂_ρE_z are continuous (μ equal, so H_φ = ∂_ρE_z/(jωμ₀) is continuous). That is a 2 × 2 solve for (c, s).
  - Wire current I = 2πaH_φ(a) = −4E₀c/(jωμ₀), because g′(a) = −2/(πa) by the J/Y Wronskian (E.19).
- **TE_z, n = 1 (transverse E).**
  - Coating: H_z = H₀·c·h(ρ) cos φ, with h = J₁(k₁ρ)Y₁′(k₁a) − Y₁(k₁ρ)J₁′(k₁a), so h′(a) = 0 ⇔ E_φ(a) = 0.
  - Host: H_z = H₀[J₁(k₂ρ) + t H₁⁽²⁾(k₂ρ)] cos φ, with `IMP-4`'s H₀ = −2jωε_c,host E₀/k₂ (incident E₀ŷ on the axis).
  - At ρ = b, H_z and ε_c⁻¹∂_ρH_z are continuous, with ε_c,coat/ε_c,host = k₁²/k₂². That is a 2 × 2 solve for (c, t).
  - E = ∇ × (H_zẑ)/(jωε_c) uses the region's own ε_c, so in the coating E carries the factor k₂²/k₁². ω cancels out of E, as in `IMP-4`.
- **Quasi-static TE (k → 0, ratio r = k₁²/k₂² kept),** derived independently. Let q = a²/b² and Δ = (1 − q) + r(1 + q).
  - Host: E_ρ = (E₀ + B/ρ²) sin φ and E_φ = (E₀ − B/ρ²) cos φ, with B = E₀b²[r(1 + q) − (1 − q)]/Δ.
  - Coating: E_ρ = −C(1 + a²/ρ²) sin φ and E_φ = −C(1 − a²/ρ²) cos φ, with C = −2E₀/Δ.
  - At r = 1 this is `IMP-4`'s static field (B = E₀a²).

**Code** (additive in `utils/analytical.py`, appended after `pec_cylinder_tm0_current`; scipy imported inside the functions):
- **`coated_pec_cylinder_tm0_field(points, *, k_host, k_coat, radius, coat_radius, e0=1.0, interior="zero", region=None)`** and **`coated_pec_cylinder_te1_field(…, omega, …)`**.
  - Points have shape (N, 3); the output is complex (N, 3) Cartesian E.
  - `region=None` selects by ρ: below a it applies `interior` ("zero" / "continue"), as `IMP-4` does; a ≤ ρ < b is the coating, ρ ≥ b the host.
  - `region="coat"` / `"host"` forces that region's formula at every point. Step 1 needs this for per-cell-tag quadrature, and (IF) needs it to evaluate both sides at exactly ρ = b.
  - Validate 0 < radius ≤ coat_radius. **coat_radius == radius is legal**: the solve stays non-singular because g′(a) and h(a) are non-zero by the Wronskian.
- **`coated_pec_cylinder_tm0_current(*, k_host, k_coat, radius, coat_radius, omega, e0=1.0)`** returns −4E₀c/(jωμ₀).
- **Private helpers** `_coated_cylinder_coefficients(drive, …, condition="exact")` (the 2 × 2 solve) and `_coated_cylinder_radial(rho, drive, region, …)` (f, f′, f″ from `jv`/`yv`/`jvp`/`yvp`/`hankel2`/`h2vp`).
  - A private `_condition="swapped"` on both public functions applies the *other* polarisation's flux condition; it exists for the controls only. This is `IMP-4`'s private-keyword pattern (`_hankel_kind`, `_condition`).
- **The test:** `tests/validation/test_imp5_coated_cylinder_closed_form.py`.
  - Fixture constants are imported from `test_imp1_step3_exact_comparand.py` (`K_FIXTURE`, `K0`) and `test_imp1_pec_sphere_lossy_host.py` (`FREQUENCY_HZ`, `E0`, `TH15_BETA_BAND`).
  - `IMP-4`'s three public functions are imported as the limit comparands.

**Anchors (asserted).**
- **(Z) scipy takes complex arguments for Y.** `yv`, `yvp` at n = 0 and 1 agree with the ascending series (docstring-derived) at z = k_host·a to ≤ 1e-13 relative. The (Lε) limit evaluates Y at a complex argument, so this is checked rather than assumed (`IMP-4`'s `test_scipy_accepts_complex_arguments`, extended to Y).
- **(B1) Bessel's equation (Jin (E.1)) in each region with that region's k**, for both radial functions of both drives, at 64 radii in [a, b] (k₁) and 64 in [b, c] (k₂). The bound is ≤ 1e-12 of the term sum \|f″\| + \|f′/ρ\| + \|(k² − n²/ρ²)f\|. `IMP-4` measured 1.8e-16 – 2.6e-16.
- **(B2) Wronskian (Jin (E.19)).**
  - J_nY_n′ − J_n′Y_n = 2/(πz) at z ∈ {k₁a, k₁b}.
  - J_nH_n⁽²⁾′ − J_n′H_n⁽²⁾ = −2j/(πz) at z ∈ {k₂b, k₂c}.
  - Both for n ∈ {0, 1, 2}, ≤ 1e-13 relative.
- **(P) PEC trace at ρ = a** (`region="coat"`, 64 points). On the Cartesian output, φ̂·E and ẑ·E are ≤ 1e-14·E₀ for both drives. `region="coat"` replaces `IMP-4`'s `interior="continue"` fix for `hypot` landing a hair below a.
- **(IF) Interface at ρ = b, 64 points, `region="coat"` against `region="host"`.** Each residual is ≤ 1e-12 of the sum of the two sides' magnitudes (the scale before cancellation):
  - tangential E (TM ẑ·E; TE φ̂·E, ẑ·E) — Jin (1.46);
  - tangential H (TM H_φ from the radial f′; TE H_z from f) — Jin (1.47);
  - TE normal D, k₁²·E_ρ,coat = k₂²·E_ρ,host — Jin (1.48), dependent but checked on the Cartesian output, where it catches a wrong ε factor in the conversion.
- **(I) Wire current.** 2πa·H_φ(a), from the analytic coating-region derivative, equals `coated_pec_cylinder_tm0_current` to ≤ 1e-12 relative.
- **(Lε) Coating → host.** With k_coat = k_host (complex), on 200 points with a ≤ ρ ≤ c and `region=None`, both coated fields equal `pec_cylinder_tm0_field` / `pec_cylinder_te1_field`, and the current equals `pec_cylinder_tm0_current`, to ≤ 1e-12 relative to max\|E\| (or \|I\|). Predicted ≈ 3e-16: the coefficient misses were computed at 3.4e-16 (TM) and 2.9e-16 (TE) by an independent pure-Python series.
- **(Lb) Thickness → 0.** With coat_radius = radius exactly, the same three equalities hold to ≤ 1e-12.
- **(Q) TE quasi-static limit.** With k_host → s·k_host and k_coat → s·k_coat at s = 1e-6 (r kept), on 200 points with a ≤ ρ ≤ c spanning both regions, the field equals the quasi-static form above to ≤ 1e-9·E₀ (`IMP-4`'s (Q) band). Predicted ≲ 1e-11.
- **(L) Layout.** A dolfinx-style adapter on x of shape (3, N) returns the (N, 3) output transposed, bitwise, for both drives (`IMP-4` (L)).

**Printed (predicted, never asserted; rule (e)).** All predictions come from an independent pure-Python series computed while drafting. That series reproduces `IMP-4`'s measured H₀⁽²⁾(ka) = 0.644055 + 2.63483j.
- k_host, k_coat, \|k_coat\|b ≈ 3.48e-3, r = k₁²/k₂² ≈ 7.39e-3 + 1.515e-2j.
- **TM: the coating is invisible at leading order.** \|s/s₀ − 1\| ≈ **3.1e-6** and \|I/I_IMP-4 − 1\| ≈ 4.3e-5, with I ≈ 6.9305e-4 − 2.8348e-3j A per V/m.
  - Why: in the quasi-static sleeve E_z = C ln(ρ/a), and the matching reproduces `IMP-4`'s H₀⁽²⁾(ka) denominator exactly. With μ equal, ε enters only through k²ρ² terms.
- **TE: the coating reverses the dipole.** t/t₀ ≈ −2.158 + 0.166j, against the quasi-static B/(E₀a²) = −2.159 + 0.170j.
  - In the sleeve C/E₀ ≈ −3.53 + 0.14j.
  - At ρ = b: \|E_ρ(b⁻)\|/\|E_ρ(b⁺)\| ≈ 59.4, \|E_ρ(b⁺)\| ≈ 0.086·E₀ and \|E_φ(b)\| ≈ 1.960·E₀. The tissue next to the insulation sees the normal field suppressed and the tangential field doubled, a void-like response. That is the reverse of the bare wire's normal doubling (`IMP-4`: 2.00).
- **Continuous annulus means ⟨\|E\|²⟩/E₀²** (step 1's comparands):
  - TE: coating (a, b) 17.99; host (b, 1.5b) 1.411, against the bare wire's 1.088 on the same host annulus, i.e. near-lead SAR +30 % in TE.
  - TM: coating 3.640e-3 (= `IMP-4`'s TM annulus); host 2.296e-2.

**Negative controls** (rule (e): predicted, not measured on this fixture; each is asserted only as missing its anchor's band, with the measured value printed beside the prediction).
- **(c1) TE: the uncoated `IMP-4` comparand against the coated β.**
  - **Measurement.** The complex LSQ amplitude, on step 1's probe circles (ρ = 1.2b, 1.5b, 2b; 24 points each; z = 0), of the coated scattered field fitted to `pec_cylinder_te1_field` − E_inc. The mode shapes are identical in the host, so β = t₀/t.
  - **Prediction:** β ≈ −0.461 − 0.035j, so \|β − 1\| ≈ **1.461**. The ceiling is 1.461/0.04886 ≈ 30× the band.
  - **Asserted:** \|β − 1\| > `TH15_BETA_BAND`.
- **TM has no "uncoated misses" control.** The same fit gives \|β − 1\| ≈ 3.1e-6, far inside the band. It is printed as a record, and it is the reason (c2) exists.
- **(c2) TM, swapped flux condition** (`_condition="swapped"`: ε-weighted ∂_ρE_z continuity, TE's condition, applied to E_z).
  - Tested against (IF)'s H_φ residual. The derivative ratio becomes r, so the predicted residual is \|r − 1\|/(\|r\| + 1) ≈ **0.98**. Asserted > 1e-3.
  - Its β separation from the correct coated TM is printed: predicted \|s_w/s − 1\| ≈ 0.104.
- **(c3) TE, swapped flux condition** (unweighted ∂_ρH_z continuity: the ε jump is dropped).
  - Tested against (IF)'s φ̂·E residual. The predicted residual is ≈ **0.98**. Asserted > 1e-3.
  - Its β separation is printed: predicted ≈ 1.461, the same as (c1) to 1e-4, because at \|k\|b ≪ 1 dropping the ε weight *is* omitting the coating.
- (B1) and (P) still pass under (c2) and (c3). That is why (IF) is the discriminating anchor.

**Window.** Foreground, Bash `timeout: 660000`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.
- Predicted 55–75 s: `IMP-4` step 0 W1 took 55 s with the same environment and lightweight modules (`20261003T172120Z_IMP-4-step0-W1.log`); this module adds ≈ 10–20 s of numpy.
- `IMP-4`'s closed-form module and `test_analytical_lightweight.py` are the rule (c) re-gate of the additive `analytical.py` change.
- `PYTHONPATH` carries `:/workspace` because the module imports `tests.validation.*`.

`scripts/testing/run_and_log.sh IMP-5-step0-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/unit/test_analytical_lightweight.py tests/validation/test_imp4_pec_cylinder_closed_form.py tests/validation/test_imp5_coated_cylinder_closed_form.py -s -v --tb=short'"`

**Traps (paid for).**
- **Convention.** e^{+jωt}, H⁽²⁾ outgoing in the host. `K_FIXTURE` is the principal root (Re > 0, Im < 0); never `np.sqrt` of a negated ε_c.
  - The sleeve uses J/Y with real k₁ (Y is regular there: the sleeve excludes ρ = 0). (Lε) puts a complex k into Y, hence (Z).
- **The ε ratio is k₁²/k₂²** (μ = μ₀ in both). Never rebuild ε_c from σ and εᵣ with a separately typed ε₀: k_coat = `K0`·√3, with `K0` imported.
- **Derivatives** are analytic or recurrence derivatives only. No finite-difference check is pre-registered (`IMP-2`'s (F) ruling).
- **Normalise residuals by the parts before cancellation,** never by \|f\| (item 190's lesson). Y₁ ≈ −270 and J₁ ≈ 1e-3 in the sleeve.
- **The ρ = 0 axis** is outside every sample set. Guard the TE conversion's 1/ρ.
- **Run hygiene.** Print on rank 0. The real build skips the four complex-only environment tests, as expected. Use `-s`. Rule (i): no edit after the last window. Rule (f) reaches a parameter (s, b), never a bound.
- **Shared file:** append only; rebase if another item appends first.

**Scope.** Numpy only: two single-mode drives on an infinite coated PEC cylinder, with one lossless coating at one thickness. **Not included:**
- the FEM (step 1);
- the plane-wave multi-mode series;
- **k_z ≠ 0 propagation along an insulated lead** (the insulated-antenna / transmission-line regime). That is where insulation actually controls lead heating, so TM's coating-invisibility is a k_z = 0 statement and must not be read as "insulation does not matter";
- a lossy or thin real-lead coating;
- a finite lead or tip (`IMP-2`, `IMP-3`);
- a coil, heating, or published data.

**Status it can move.** The proposed §7 `IMP-5` row moves ⬜ → 🟡 on its step-0 Done-when, which this item matches (rule (j)). The review adds the row in the commit that queues this item and names it in the dashboard for the weekly to ratify (rule (5) lane 3).

**Negative result.**
- **(Z), (B1), (B2), (P) or (L) red:** a code bug. Fix it in-slot, run W1b, and disclose.
- **(IF) red with (B1) green:** a wrong ε factor or region selection. Print both sides and their ratio, then stop.
- **(Lε) or (Lb) red:** the coated algebra disagrees with `IMP-4`. Print the coefficient pairs, then stop.
- **(Q) red:** print the miss at s ∈ {1e-4, 1e-5, 1e-6}. If it scales as s² ln s, s is a parameter: re-register it under rule (f), disclosed; otherwise stop.
- **A control inside its band:** the anchor does not discriminate. Stop and report.
- **On a stop:** park on `attempt/IMP-5-step0-<UTC>` and mark BLOCKED in the same commit (rule (d)); the row stays ⬜.
- There is no PREDICTED RED line.

**Outcome (2026-10-04 00:30 UTC, 19:30 CDT slot; implementer): ✅ step 0, row ⬜ → 🟡.** Windows `20261004T003434Z_IMP-5-step0-W1.log` (39 passed, 4 complex-only skips, 63 s) and `20261004T003556Z_IMP-5-step0-W1b.log` (same, Status 0, 59 s; the closing window). **Disclosed:** W1 was all green, but its (Q) area-uniform draw put 1 of 200 points in the sleeve; the test now draws 100 coat + 100 host points and asserts ≥ 90 per region (no band moved), and W1b re-ran the same command under rule (i). Anchors: (Z) 2.3e-16; (B1) ≤ 3.2e-16; (B2) ≤ 6.7e-16; (P) ≤ 8.9e-16; (IF) ≤ 9.0e-15 (TE normal D); (I) 4.5e-16; (Lε) ≤ 6.8e-16 (predicted ≈ 3e-16); (Lb) ≤ 5.6e-16; (Q) 4.3e-14 at s = 1e-6 (W1 1.75e-14 on the one-coat-point draw); (L) bitwise. Readings vs prediction: TM \|s/s₀ − 1\| 3.084e-6 (3.1e-6), \|I/I_IMP-4 − 1\| 4.27e-5 (4.3e-5), I = 6.9305e-4 − 2.8348e-3j; TE t/t₀ −2.1584 + 0.16582j vs quasi-static B/(E₀a²) −2.1587 + 0.17033j; C/E₀ −3.527 + 0.136j; \|E_ρ(b⁻)\|/\|E_ρ(b⁺)\| 59.36, \|E_ρ(b⁺)\| 0.0859·E₀, \|E_φ(b)\| 1.960·E₀. Annulus ⟨\|E\|²⟩/E₀²: TE coat 17.987, TE host 1.411 (bare 1.0878, ×1.297), TM coat 3.6404e-3, TM host 2.2963e-2. Controls: (c1) \|β − 1\| 1.461 (> 0.04886); (c2) H_φ residual 0.9763, \|s_w/s − 1\| 0.104; (c3) E_φ residual 0.9763, \|β − 1\| 1.461. **Brief wording miss (not a gate):** under (c3) the printed, unasserted (P) TE φ reads 1.005e-14, a hair over the 1e-14 band the brief says "still passes" — the swapped coefficient enlarges the sleeve field ≈ 11×, so rounding scales with it; (IF) still discriminates. Code: additive block after `pec_cylinder_tm0_current` in `utils/analytical.py` (three public functions, four private helpers, private `_condition`); test `tests/validation/test_imp5_coated_cylinder_closed_form.py`. `IMP-4`'s closed-form module and `test_analytical_lightweight.py` re-gated green in the same window (rule (c)).

### §9 item 211 — 2026-10-03 (lane 3)

*Drafted by a general-purpose agent for the 2026-10-03 18:00 daily review; the review reviews it.*

**`IMP-5` step 1: the FEM half.** It is `IMP-4` step 1's fixture with a conforming εᵣ 3 sleeve a < ρ < b meshed as its own cell tag and assigned through `material_map`. The outer wall carries step 0's exact coated trace, and the end planes carry `IMP-4`'s polarisation-specific conditions.
- **Executor:** implementer. **Tier:** standard, `-n 2`; real build for W0, complex for W1a–W2b. **Class:** lane-3 closed form (measured median ≈ 12 min).
- **Cost:** ≈ 12 measured / ≈ 70 predicted min (15 fixed + ≈ 30 code + five windows of ≈ 1.5–2.5 min + ≈ 10 spare).
- **Depends on item 210; if it did not land, skip.** The check is `grep -n 'def coated_pec_cylinder_tm0_field' src/fem_em_solver/utils/analytical.py` non-empty, plus 210's W1 green in `test-results.md`.
  - This item **does not edit `utils/analytical.py`**. Should an edit become necessary, the file is append-only; rebase if another item lands first.

**Why the fixture is still exact.** `IMP-4` step 1's argument (`IMP-4.md`, item 206) carries region by region.
- **TM:** E = E_zẑ in both regions, so n̂ × E = 0 on z = const and the **PEC** ends are exact: pec tags (1, 2, 3).
- **TE:** H = H_zẑ in both regions, so n̂ × H = 0 there and the **natural (PMC)** ends are exact: pec tags (1, 2).
- **The interface needs no condition.** The curl-curl weak form imposes tangential E and tangential μ⁻¹ curl E continuity naturally. Nédélec elements carry tangential continuity and leave the normal component free to jump by 59×.
- **Truncation:** the outer wall carries the exact trace. Only discretisation and faceting remain.

**Fixture.**
- **Imports from `test_imp1_pec_sphere_lossy_host.py`:** `FREQUENCY_HZ`, `SIGMA_HIGH`, `EPSILON_R_HOST`, `E0`, `DEGREE` = 2, `OUTER_BOUNDARY_TAG` = 1, `CAVITY_TAG` = 2, `TH15_BETA_BAND`, `MAT4_BAND`, `_assemble_real`, `_cavity_dofs`.
- **Imports from `test_imp1_step3_exact_comparand.py`:** `K_FIXTURE`, `K0`.
- **Imports from `test_imp2_spheroid_fem.py`:** `QDEG`, `_numpy_quadrature_mean`, `_owned_masked_cells`, `_rss_gib`.
- **Imports from `test_imp4_pec_cylinder_fem.py`:** `A`, `C`, `LENGTH`, `END_TAG`, `RESOLUTION_WALL`, `RESOLUTION_FAR`, `CELL_GUARD`, `TAG_BAND`, `M_BAND`, `S_BAND`, `N_PROBE`, `PEC_TAGS`, `_say`, `_ncells`, `_circle`, `_eval_complex`, `_te_mode`, `_incident`, `_static_te_numpy`, `_h_wall_over_a`.
- **Never re-derive these.** New constants: `B = 1.5 * A`, `COAT_TAG = 4`, `EPS_R_COAT = 3.0`, `K_COAT = K0 * sqrt(3)`.

**Code.**
- **Mesher (one additive keyword; rule (c), disclosed).** Add keyword-only `coat_radius: Optional[float] = None` to `MeshGenerator.cylindrical_domain` (`io/mesh.py:941`), passed through to `_cylindrical_hole_domain` (`io/mesh.py:70`) with default `None`.
  - Passing it without `as_hole=True` raises, as `resolution_wall` does. It must satisfy a < coat_radius < c, or raise.
  - `None` takes the existing hole branch **untouched**.
  - When it is set:
    - Build the sleeve as `cut(cyl_b, cyl_a)` and the host as `cut(cyl_c, cyl_b)`, then `occ.fragment` them so the ρ = b face is **shared** (conforming). Raise unless exactly two volumes remain.
    - Identify the volumes by bbox r_max (b against c). The host is volume tag **2** (unchanged); the sleeve is the new cell tag **4** "coating".
    - Facet groups come from `getBoundary` of **both** volumes combined, so the shared face drops out. They are classified z first exactly as now: tag 3 ends (now four annular faces), tag 1 outer wall, tag 2 cavity.
    - No interior-interface facet group is emitted, because the material jump is carried by cell tags.
    - The sizing field is unchanged: a Distance from the cavity faces gives h(ρ) = ρ/6.
    - Print the face count per group, the cells per volume tag and the min SICN.
- **Module** `tests/validation/test_imp5_coated_cylinder_fem.py`:
  - **Material:** `HomogeneousMaterial(sigma=SIGMA_HIGH, epsilon_r=EPSILON_R_HOST)` as the default, plus `material_map={COAT_TAG: HomogeneousMaterial(sigma=0.0, epsilon_r=EPS_R_COAT)}`, and `cell_tags`.
  - **Trace:** step 0's fields with `region=None, interior="zero"`, through `IMP-4`'s `(3, N)` adapter.
  - **Incident split:** `IMP-4`'s `_incident(drive)`, which is host k and host ε_c. The scattered basis is E_s = E_cf(`region="host"`) − E_inc.
  - **Probes:** 24 points on each of ρ = 1.2b, 1.5b and 2b at z = 0, all in the host. β is `IMP-4`'s complex LSQ over all components; raise unless every point evaluates.
  - **Two masks, both built from cell tags (numpy DG0, never a radius test at the interface):**
    - **coat** = cell tag 4;
    - **host annulus** = cell tag 2 with midpoint x² + y² < (1.5b)².
    - Φ_h,m = ⟨\|E_h\|²⟩/E₀² by UFL at `quadrature_degree` `QDEG`.
    - Φ_cf,h,m is `_numpy_quadrature_mean` over the same owned cells, using the closed form with **`region` forced to the mask's region** ("coat" on tag 4, "host" on tag 2). Quadrature points in the chord slivers lie on the wrong side of the true circles, at ρ < a or ρ < b, and region-by-radius would hand a host point the coating formula: a 59× E_ρ error.
    - Φ_cont,m is a numpy radial / angular quadrature of step 0's field over the true annulus, printed beside Φ_cf,h,m to isolate faceting.
  - **Solve:** a module-local `_solve_coated(drive, pec_tags, coated=True)` following `IMP-4`'s `_solve_cylinder`. `coated=False` drops the `material_map`, which is control (c1). It returns β, both Φ triples, the cavity dof count and max\|E\|, z-invariance, cells, seconds and summed RSS.
  - **Tests:**
    - `test_census_and_numpy` (W0, **not** `complex_only`);
    - `test_tm_drive` (W1a), `test_te_drive` (W1b);
    - `test_control_te_coating_omitted` (W2a), `test_control_tm_natural_cavity` (W2b).
    - All four solve tests are `complex_only`, with one solve per window to hold the standard tier.

**Anchors (asserted).**
- **W0 (census, before any solve):**
  - **(Guard)** global cells ≤ `CELL_GUARD` (60 000). Predicted 26–33 k: `IMP-4`'s 25 568, plus the conforming ρ = b surface, which the ρ/6 field already resolves.
  - **(R)** `cylindrical_domain(as_hole=True, coat_radius=None)` reproduces `IMP-4` W0's **25 568** cells exactly (`20261003T200309Z_IMP-4-step1-W0.log:465`). This is the `None` branch's bit-identity, alongside the re-gates below.
  - **(Tag)** reduced `ds(t)` areas (pinned `quadrature_degree` 2): tag 1 against 2πcL, tag 2 against 2πaL, tag 3 against 2π(c² − a²), each within `TAG_BAND` (2 %). `IMP-4` measured −0.088 / −0.088 / −0.455 %.
  - **(Vol)** reduced `dx(4)` and `dx(2)` volumes against π(b² − a²)L and π(c² − b²)L within `TAG_BAND`. Predicted ≈ −0.5 % (polygon deficit). Swapped tags read O(1). This is the check that `material_map` targets the sleeve.
  - **(M′)** On each of the two masks, the numpy quadrature of `IMP-4`'s polynomial static TE field equals its UFL cell integral to ≤ `M_BAND` (1e-12).
  - **(S)** The incident split, on 64 points in [1.2b, 2b]: TM (E_cf,z − E₀J₀(k₂ρ))/H₀⁽²⁾(k₂ρ) is constant = s; TE E_cf − E_inc = t·(the H-mode of `_te_mode`). Each holds to ≤ `S_BAND` (1e-12).
- **W1a / W1b, per drive** (TM pec tags (1, 2, 3); TE pec tags (1, 2)):
  - (i) \|β − 1\| ≤ `TH15_BETA_BAND` (4.886 %);
  - (ii) \|Φ_h/Φ_cf,h − 1\| ≤ `MAT4_BAND` (10 %) on the **coat** mask;
  - (iii) the same on the **host-annulus** mask;
  - (iv) cavity dofs located > 0 (reduced), with max\|E\| < 1e-12 (reduced with `MPI.MAX`).

**Printed (predicted, never asserted; rule (e)).**
- **TM.**
  - \|β − 1\| ≲ 0.2 % (`IMP-4` measured 0.036 %; the coating moves s by 3e-6).
  - Coat Φ miss +1 to +3 %. This is `IMP-4`'s TM annulus a < ρ < 1.5a, measured +1.01 % there, wall-faceting-led; Φ_cont 3.640e-3.
  - Host Φ miss +0.3 to +2 %, with Φ_cont 2.296e-2.
  - TM current β·I against step 0's I.
- **TE.**
  - \|β − 1\| ≲ 2 %: the dipole goes as b_eff², the faceting deficit is ≈ −0.46 %, and `IMP-4` TE measured 0.774 %.
  - Coat Φ \|miss\| ≲ 2 %: the sleeve thins by ≈ 0.23 % under faceting, giving ≈ +0.46 % in \|E\|². Φ_cont 17.99.
  - Host Φ \|miss\| ≲ 1.5 %, with Φ_cont 1.411.
- **Both drives:** cells per tag, degree-2 dofs, seconds, summed RSS, h_wall/a, the sleeve's max edge/(b − a) (predicted 0.5–0.8), min SICN, and z-invariance (`IMP-4` measured 7.3e-3 / 1.3e-2).

**Negative controls** (rule (e): predicted by the drafting series with the exact coated trace at c, not measured; each is asserted only as missing its band).
- **(c1) TE, coating omitted** (W2a; `coated=False`: host material in tag 4, same mesh, pec tags (1, 2)). The material map must reach the sleeve, or TE cannot tell.
  - **Prediction.** β ≈ −0.476 − 0.036j, \|β − 1\| ≈ **1.476** (range [1.3, 1.65]), with a ceiling of 1.476/0.04886 ≈ **30×**. Coat Φ miss ≈ **−91.4 %** (range [−95, −85] %). Host Φ miss ≈ **−17.8 %** (range [−22, −14] %, 1.8× the band).
  - **Asserted:** \|β − 1\| > `TH15_BETA_BAND`, and both \|Φ misses\| > `MAT4_BAND`.
- **(c2) TM, natural cavity** (W2b; PMC wire, pec tags (1, 3), coated).
  - **Prediction.** β ≈ 0.590 − 0.123j, \|β − 1\| ≈ **0.428** (range [0.25, 0.6]), with a ceiling of ≈ **8.8×**. Coat Φ miss ≈ **+7 940 %** (range [+3 000, +15 000] %). Host Φ miss ≈ **+1 175 %** (range [+500, +2 500] %).
  - **Asserted:** \|β − 1\| > `TH15_BETA_BAND`, and both \|Φ misses\| > `MAT4_BAND`.
  - **The calibration of the method:** the same trace-pinned annulus argument predicted `IMP-4`'s c1 at 0.49, measured 0.489.
- **Not run:** TM with the coating omitted. It is predicted at \|β − 1\| ≈ 2.4e-6, the step-0 invisibility. A solve would only restate it, so no window is spent there.

**Windows** (foreground, Bash `timeout: 660000`; each harness call is the whole Bash command, bare and repo-relative; `-s`; node IDs, never `-k`; no pipes).
- **W0** (real build). Predicted 100–140 s: `IMP-4` W0 took 96 s with the same two re-gates plus a closed-form module (`20261003T200309Z_IMP-4-step1-W0.log:620`), and this window adds one census (mesh 2.1 s there).
  - It is the rule (c) re-gate of `coat_radius=None`: `test_cylindrical_domain.py`, `test_cylinder.py`, and `IMP-4`'s census, whose tag areas must match its record.
  - Step 0's module is re-run as the dependency check.

  `scripts/testing/run_and_log.sh IMP-5-step1-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_cylindrical_domain.py tests/solver/test_cylinder.py tests/validation/test_imp4_pec_cylinder_fem.py::test_census_and_numpy tests/validation/test_imp5_coated_cylinder_closed_form.py tests/validation/test_imp5_coated_cylinder_fem.py::test_census_and_numpy -v -s --tb=short'"`
- **W1a / W1b / W2a / W2b** (complex, one solve each). Predicted 70–110 s each.
  - `IMP-4` W1 was 123 s for two solves of 36.2 s and 28.2 s (`20261003T200454Z_IMP-4-step1-W1.log:133,139,296`), i.e. ≈ 57 s of fixed overhead plus one solve of ≈ 30–50 s at ≤ 33 k cells.
  - Two solves per window would sit at ≈ 125–160 s against the 180 s ceiling, which is why each window holds one.

  `scripts/testing/run_and_log.sh IMP-5-step1-W1a "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp5_coated_cylinder_fem.py::test_tm_drive -v -s --tb=short'"`

  `scripts/testing/run_and_log.sh IMP-5-step1-W1b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp5_coated_cylinder_fem.py::test_te_drive -v -s --tb=short'"`

  `scripts/testing/run_and_log.sh IMP-5-step1-W2a "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp5_coated_cylinder_fem.py::test_control_te_coating_omitted -v -s --tb=short'"`

  `scripts/testing/run_and_log.sh IMP-5-step1-W2b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp5_coated_cylinder_fem.py::test_control_tm_natural_cavity -v -s --tb=short'"`
- **An overrun** of any window is killed and shrunk (CLAUDE.md), never re-timed. Check `pgrep -c python3` in the container for orphaned ranks first.

**Traps (all paid for).**
- **End conditions are polarisation-specific** (`IMP-4`): TM ends PEC, TE ends natural. Swapping them is an O(1) β miss, not a small error.
- **Region by cell tag, never by radius, at either circle.** The TE normal E jumps 59× at ρ = b, and the chord slivers put quadrature points of host cells at ρ < b and of sleeve cells at ρ < a. This generalises `IMP-1` / `IMP-4`'s chord-sliver trap (item 146). Use the "zero" interior for the trace only.
- **Mesh classification.** Classify z first (`IMP-4`). Facet groups come from the **meshed** boundary of the combined volumes (`TH-15`). The fragment must leave the ρ = b face shared, or the sleeve is a disconnected solid and the weak form never couples the regions.
- **ε_c consistency.** The solver builds ε_c from (σ, εᵣ) with its own ε₀ (`core/time_harmonic.py:782–783`). The comparand uses k₂ = `K_FIXTURE` and k₁ = `K0`·√3, which are the same constants. Never mix in a separately typed ε₀.
- **Complex build.** No UFL comparison, no `ufl.max_value`, and no comparison containing `sqrt`, even on real operands (`OPS-22`, `WF-6` 4c). Every mask is numpy, from cell tags.
- **Quadrature.** Every `SpatialCoordinate`-bearing form, and the `ds` / `dx` area and volume forms, carries `metadata={"quadrature_degree": …}` (`POST-5` step 1).
- **Rank safety.**
  - `cell_tags.values`, `assemble_scalar` and local max/min are rank-local; reduce before asserting.
  - `material_map` tag validation is rank-reduced (`core/time_harmonic.py:150–180`); do not pre-check tags rank-locally.
  - Do the RSS `allreduce` **outside** any rank-0 branch (`IMP-2` W0 `20261001T200855Z`). Owned cells only.
- **JIT cache.** If an earlier window today was killed, sweep orphan `.c` stems first, and only when `pgrep -c python3` in the container reads 0, using the rubric's one bare harness call with `<ID>` = `IMP-5-step1`. Never clear the cache wholesale.
- **Run hygiene.**
  - Foreground; container `timeout -k 30`; never pipe or filter inside the harness command.
  - `PYTHONPATH` carries `:/workspace`. Print on rank 0. Rule (i): no edit to a file after its last window.
  - Rule (f) reaches a parameter (b, a resolution), never a bound. Rule (h) applies: no in-slot comparand substitution.
- **Reading TM.** A green TM gate is **not** a coating test, because TM cannot see the sleeve (step 0). The coating is tested by TE's gates and by (c1).

**Scope.** One mesh, degree 2, two single-mode drives on an infinite coated PEC cylinder with an exact trace and one lossless coating. These are accuracy-level bands, **not a convergence claim**. **Not included:**
- k_z ≠ 0 guided waves along an insulated lead (the regime where insulation sets lead heating);
- a lossy or real-thickness insulation;
- Jin's (11.63) global boundary condition;
- modes n ≥ 2;
- a finite lead or tip (`IMP-2`, `IMP-3`);
- a coil, heating, 1 g SAR, or published data.

**Status it can move.** The `IMP-5` §7 row, whose Done-when reads: "*Step 1* (FEM, degree 2, `TH-15` hole route): additive `coat_radius=` on `cylindrical_domain(as_hole=True)` (sleeve cell tag 4; `None` bit-identical); sleeve via `material_map`; exact trace; ends as `IMP-4`; per drive \|β − 1\| ≤ `TH15_BETA_BAND`, sleeve and near-lead host annulus ⟨\|E\|²⟩ each ≤ `MAT4_BAND`, cavity dofs zero; TE coating-omitted and TM natural-cavity controls miss ⇒ ✅".
- With step 0 green (🟡), this item green moves **🟡 → ✅**, auditable. The opening is the weekly's to ratify.
- The review may then add `IMP-5` to §2's closed-form-lane sentence: imposed trace, no coil, k_z = 0.

**Negative result.**
- **Over the cell guard in W0:** stop before any solve, report the counts, BLOCKED.
- **(R) red:** the `None` branch is not bit-identical. Fix it in-slot. If the image changed since `IMP-4` W0, stop instead: that is a record re-registration for `record-reconciler`, not a code bug.
- **(Tag), (Vol), (M′) or (S) red:** a mesher or helper bug. Fix it in-slot, run W0b, and disclose.
- **A re-gate red:** fix it; if it was red at `HEAD` already, cite known-issues and continue.
- **(i), (ii) or (iii) red with (iv), (Vol) and (M′) green:** print h_wall/a, the sleeve's edge/(b − a), each Φ_cf,h against Φ_cont, and z-invariance. Park on `attempt/IMP-5-step1-<UTC>` and mark BLOCKED in the same commit (rule (d)); no band moves.
  - If TE (ii) is the red and it tracks the sleeve's resolution, the review may re-pose `RESOLUTION_WALL` (a/10) under rule (f).
- **A control inside its band:** the sleeve's material or the cavity condition is not reaching the probes, which points at a tag or `material_map` bug. Stop and report.
- There is no PREDICTED RED line.

**Outcome (2026-10-04 22:30 CDT slot, implementer): ✅ complete — five windows green on the first run, no W0b.**
- **Windows** (`-n 2`): W0 real `20261004T034835Z_IMP-5-step1-W0.log` 26 passed / 4 complex-only skipped, 104 s; W1a `20261004T035027Z_IMP-5-step1-W1a.log` 133 s; W1b `20261004T035247Z_IMP-5-step1-W1b.log` 168 s; W2a `20261004T035543Z_IMP-5-step1-W2a.log` 162 s; W2b `20261004T035832Z_IMP-5-step1-W2b.log` 152 s (13 passed each). Complex windows ran over the predicted 70–110 s (solves 49–57 s; the rest is the environment tests); W1b sits 12 s under the 180 s ceiling, so one solve per window is the right size.
- **W0:** 25 878 cells (host 16 255, sleeve 9 623) ≤ 60 000; (R) `coat_radius=None` 25 568 cells = `IMP-4` W0, `IMP-4` census tag areas −0.0882 / −0.0878 / −0.4550 % as recorded; (Tag) same three; (Vol) sleeve −0.3368 %, host −0.3743 %; (M′) 6.44e-15 / 7.11e-15; (S) TM 3.04e-16, TE 4.82e-16, t/t₀ −2.1584 + 0.16582j. Printed: h_wall/a 0.383, min SICN 0.2965, sleeve max edge/(b − a) 0.858 (predicted 0.5–0.8 — over, printed only).
- **Gates:** TM β 0.999703 − 0.000089j (0.0310 %), sleeve +1.0438 %, host +0.3631 %, cavity 10 704 dofs max 0, z-inv 4.87e-3, β·I 6.92592e-4 − 2.83402e-3j vs 6.93050e-4 − 2.83480e-3j; TE β 0.994548 − 0.000062j (0.5452 %), sleeve −0.1656 %, host −0.1850 %, cavity max 7.3e-34, z-inv 2.36e-2. Φ_cont 3.640e-3 / 2.296e-2 (TM), 17.987 / 1.4110 (TE), as predicted.
- **Controls:** c1 (TE, coating omitted) β −0.473107 − 0.035938j, \|β − 1\| 1.4735, sleeve −91.46 %, host −18.01 %; c2 (TM, natural cavity) β 0.590535 − 0.122462j, \|β − 1\| 0.4274, sleeve +8 024 %, host +1 166 % — all inside the drafting series' predicted ranges.
- **Rule (c), disclosed:** keyword-only `coat_radius=None` on `MeshGenerator.cylindrical_domain` → `_cylindrical_hole_domain` in `io/mesh.py`; raises without `as_hole=True` or outside a < b < c; the `None` geometry block is unchanged inside an `if coat_radius is None:` branch; shared facet classification now collects a list of volumes and passes `combined=True` to `getBoundary` explicitly (gmsh's default; (R) and the tag areas confirm the `None` path). Coated path: `cut` / `cut` / `fragment`, exactly two volumes, host tag 2 / sleeve tag 4 by bbox radius, 1 outer / 1 cavity / 4 end faces.
- **Deviations:** the test module imports the private `_coated_cylinder_coefficients` from `utils/analytical.py` (read-only) for the (S) targets, which the brief left unsourced; new module-local helpers `_tag_mask`, `_continuous_mean`, `_sleeve_edge_over_t`, `_probe_points` and `IMP4_W0_CELLS`. Rule (i): `io/mesh.py` last edited 03:47:00 UTC, the module 03:48:27, both before W0 (03:48:35).

### §7 row text moved here verbatim — 2026-10-04 03:00 daily review (row over 2 KB)

> `IMP-4`'s host and a = 1 mm; lossless εᵣ 3 sleeve to b = 1.5 mm; TM_z n = 0 and TE_z n = 1. Series: Harrington ch. 5 / Balanis ch. 11 (not held, not re-verified); Jin 3e uses it (§10.2.3 Fig. 10.11) without printing it, so the docstring derives it from (1.44)–(1.52), (11.56), (E.1), (E.19). **Done-when:** *step 0* (numpy): Bessel residuals per region ≤ 1e-12 (scale before cancellation); Wronskians ≤ 1e-13; PEC traces at a ≤ 1e-14·E₀; tangential E, H and TE normal D continuous at b ≤ 1e-12; 2πaH_φ(a) = closed-form current ≤ 1e-12; coating ε → host and b → a each reproduce `IMP-4` ≤ 1e-12; TE quasi-static ≤ 1e-9·E₀; three controls miss ⇒ 🟡.

Audited PASS 2026-10-04 03:00 review (journal `docs/planning/reviews/2026-10-04-daily.md` §3).
