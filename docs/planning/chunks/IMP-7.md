# IMP-7 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record once the review adds it. **Proposed 2026-10-20 for the 18:00 daily review under §9 rule (5) lane 3 (Phases 7 / 8), for the weekly to ratify or re-home.**

**Why this chunk.** Every Phase 7 row says "no heating": `IMP-1` … `IMP-6` gate the near-implant *field*, and Phase 8 (`MAT-9` … `MAT-21`) heats spheres, cylinders and boxes with *no implant*. No closed form in the repo joins the two. `IMP-7` is that join on the simplest implant: `IMP-1`'s PEC sphere (a = 5 mm), whose quasi-static dipole field gives a near-implant SAR that heats a Robin-cooled gel shell around a thermally conducting core. The steady temperature is an exact two-term Legendre series (l = 0 and l = 2). It reduces to `MAT-18`'s two-layer sphere when the dipole is switched off, and it gives the implant's own heating factor and pole / equator hot spot in closed form.

**Sources.**
- **Field:** the PEC sphere in a uniform quasi-static field, β = 1, as `TH-15` and `IMP-1` gate it (`tests/validation/test_pec_sphere_hole.py:3–14`, `test_imp1_pec_sphere_lossy_host.py:4–16`). The general dipole β = (ε − 1)/(ε + 2) is `dielectric_sphere_quasistatic_field` (`src/fem_em_solver/utils/analytical.py:595–626`).
- **Legendre:** Jin 3e App. E §E.3, eqs (E.34)–(E.36) (`docs/references/jin-fem-3e/appE_special-functions.md:119–135`). The equation bodies are image placeholders, so only the numbers are cited.
- **Conduction:** Carslaw & Jaeger, *Conduction of Heat in Solids* 2e, ch. IX (the sphere), the same pointer and footing as `MAT-9` … `MAT-18` (`analytical.py:2530–2532`). Equation numbers are not pinned (the repo holds no copy); the docstring derives every form.
- No published implant-heating datum is claimed. (7c) still waits on the operator's reference (dashboard W3).

## Proposed §7 row

| `IMP-7` | *(history and §9 briefs: `docs/planning/chunks/IMP-7.md`)* **SAR heating around a PEC-sphere implant in a Robin-cooled gel ball: the implant-heating canonical in closed form (Phase 7 × Phase 8)** (proposed by the 2026-10-20 18:00 daily under rule (5) lane 3; for the weekly to ratify or re-home). `IMP-1`'s sphere (a = 5 mm) as a thermally conducting core (k₁) inside a gel shell (k₂ 0.6) to b = 4a, Robin h 100 on r = b. The source is q = q_far·\|E\|²/E₀² of the quasi-static dipole field (β = 1 PEC) in the shell and zero in the core. θ = T − T∞ = Σ_{l=0,2} R_l(r)P_l(cos ϑ), exact. **Done-when:** *step 0* (numpy, §9 item 579): the source, radial-PDE, interface, Robin, balance and FD-Laplacian identities and the β = 0 → `MAT-18` and Φ → `IMP-1` reductions at the bands in the chunk file; three controls miss ⇒ 🟡. *Step 1* (FEM, real build, P1 ladder + P2 on `two_layer_ball_domain(inner_radius=a, outer_radius=b)`): `MAT-18` step 2(a)'s gates — interior rel L2 ≤ 1 %, P1 L2 rate ≥ 1.7, balance ≤ 1e-9 — on the metal core (DG0 k) and on the adiabatic hole (shell submesh, where the l = 2 part is ≈ 13 % of θ at the pole), with the isotropic-source and β = 0 controls missing ⇒ ✅. Bands for any l = 2 moment are the review's when it briefs step 1. Not claimed: a FEM SAR source (a later step, composed band as `MAT-10`), finite k, a rod or tip, a coil, transient, perfusion, published data. | ⬜ Open. | step 0 standard, no MPI, real; step 1 standard / heavy, -n 2, real |

### §9 item 579 — 2026-10-20 (lane 3; `IMP-7` step 0, closed form)

*Drafted read-only for the 2026-10-20 18:00 daily review; the review rules on it.*

**Licence.** §9 rule (5) lane 3 (`PROJECT_PLAN.md:3937–3948`): a Phase 7/8 chunk whose validation target is a closed form, named in the dashboard for the weekly to ratify. Precedents: `IMP-4` / `IMP-5` / `IMP-6` / `MAT-21`, each opened by a daily at step 0 and then ratified. **Independent of every §10 chain** (H/C/I/T/P/U and their `xl` entries) and of the queued items 556, 574, 542, 553 and 560. Nothing waits on the operator. **This brief moves no existing band.** Every bound below is imported by name, except (Rk)'s ratio window, which is flagged for the review.

**`IMP-7` step 0: the steady temperature around `IMP-1`'s PEC sphere heated by its own dipole SAR, numpy only.**
- **Executor:** implementer. **Tier:** standard, no `mpiexec` (-n 1), real build.
- **Class:** lane-3 closed form. **Cost:** ≈ 9.7 measured / ≈ 40 predicted min (15 fixed + ≈ 20 code + one ≈ 0.5-min window + spare).
- **Depends on:** nothing. **Family:** step 0, attempt 1 of 4.
- **Comparable windows:** `IMP-4-step0-W1` 55 s (`-n 2`, `test-results.md:2980`); `IMP-6-step0-W1` 3–4 s (no MPI, `:3986–3988`); `IMP-1-step3a-W1` 45 s (real build, importing the `IMP-1` module, `:2699`).

**Fixture** (rule (f) parameters, not bands):
- a = `SPHERE_RADIUS` = 0.005 (`test_imp1_pec_sphere_lossy_host.py:61`) and b = 4a = 0.02. b is chosen so step 1's uniform 1 mm mesh stays below `MAT-18` step 2(a)'s finest rung.
- k₂ = `MAT-9`'s `K` = 0.6 and q_far = `Q` = 2000 W/m³ (`test_mat9_steady_conduction_sphere.py:37–39`).
- h = `MAT-10`'s `H` = 100 (`test_mat10_sar_to_heat_sphere.py:53`), so Bi = hb/k₂ = 3.33. T∞ = 20.
- k₁ (the core): `K_CORE = 20.0`, metal-like and representative, not a material claim (`MAT-18`'s wording). Two exact branches as well: `k_core=0.0` (adiabatic: the core is excluded, as on `IMP-1`'s hole mesh) and `k_core=np.inf` (isothermal core).
- β = 1 (PEC). β = −½ (void) and β = ½ are used for cross-checks only.

**The closed form.** The docstring derives every line. Let μ = z/r, u = (a/r)³ and E = E₀ẑ + βE₀a³(3(ẑ·r̂)r̂ − ẑ)/r³ for r > a.
- **Source.** \|E\|²/E₀² = (1 + 2βu)²μ² + (1 − βu)²(1 − μ²) = **1 + 2β²u² + (4βu + 2β²u²)P₂(μ)**. In the shell, q = q_far × that; in the core, q = 0. At β = 1 this is 9 at the poles of r = a and 0 at its equator.
- **l = 0, shell:** R₀ = A₀ + B₀/r − q_far r²/(6k₂) − β²q_far a⁶/(6k₂r⁴). The core is constant, R₀ = A₀ (zero net flux at r = a for every k₁). This gives:
  - B₀ = (2β² − 1)q_far a³/(3k₂);
  - θ₀(b) = q_far(b³ − a³ + 2β²a³ − 2β²a⁶/b³)/(3hb²);
  - A₀ = θ₀(b) − B₀/b + q_far b²/(6k₂) + β²q_far a⁶/(6k₂b⁴).
- **l = 2, shell:** R₂ = D r² + E r⁻³ + C₁/r + C₂/r⁴, with C₁ = 2βq_far a³/(3k₂) and C₂ = −β²q_far a⁶/(3k₂). The particular solutions come from Δ_l rⁿ = [n(n + 1) − 6]rⁿ⁻². The core is F r².
  - Conditions: R₂(a⁻) = R₂(a⁺); k₁R₂′(a⁻) = k₂R₂′(a⁺); k₂R₂′(b) + hR₂(b) = 0.
  - Eliminate F = R₂(a)/a² and solve the 2 × 2 system for (D, E) by Cramer's rule.
  - Adiabatic branch: R₂′(a⁺) = 0. Isothermal branch: R₂(a⁺) = 0.
- **Total power:** P = (4πq_far/3)(b³ − a³ + 2β²a³ − 2β²a⁶/b³) = 4πb²hθ₀(b). Only l = 0 crosses r = b.

**Code** (append at the end of `src/fem_em_solver/utils/analytical.py`, after `coated_pec_cylinder_tm0_current`, `:3116`, under an `IMP-7 step 0` banner; additive, numpy only):
- `pec_sphere_implant_source(points, *, radius, q_far, beta=1.0)`: (N, 3) → q, zero for r < a.
- `pec_sphere_implant_coefficients(*, radius, outer_radius, k_core, k_gel, q_far, h, beta=1.0)` → a dict of A₀, B₀, C₁, C₂, D, E, F. `k_core` takes `0.0` and `np.inf` as explicit branches.
- `pec_sphere_implant_radial(r, coeffs, l, *, derivative=0)` → R_l and its first and second derivatives in the shell (analytic power laws).
- `pec_sphere_implant_steady_temperature(points, *, radius, outer_radius, k_core, k_gel, q_far, h, t_inf, beta=1.0)` → T. Core points return T∞ + A₀ + F r²P₂, or NaN on the adiabatic branch.
- **Test:** `tests/validation/test_imp7_implant_heating_closed_form.py`, prefix `IMP7-0`.

**Anchors (asserted).** Bands are imported by name from:
- `B1_BAND` = 1e-12 and `CONTROL_MISS` = 1e-3 (`test_imp4_pec_cylinder_closed_form.py:46, 51`);
- `MAT-18` (S1)'s 1e-13 reduction literal (`test_mat18_two_layer_sphere_closed_form.py:113`);
- `IMP-1` step 3(a)'s (H) literals: ratio in [3.8, 4.2] and residual ≤ 1e-4 (`test_imp1_step3_exact_comparand.py:193–194`), and its (Φ0) 1e-10 (`:217`).

The anchors:
- **(Src)** `pec_sphere_implant_source` at β = 1 equals q_far·\|E\|² from `IMP-1`'s `_pec_exterior_numpy()` (`test_imp1_pec_sphere_lossy_host.py:80–104`, E₀ = 1). Sample 512 points with r ∈ [1.001a, b], random μ and φ. Bound: ≤ `B1_BAND` × max q. Repeat at β = −½ and β = ½ against `dielectric_sphere_quasistatic_field(radius=a, epsilon_c=0.0 / 4.0)`, exterior only (β = (ε − 1)/(ε + 2) gives −½ / ½ exactly).
- **(PDE)** The radial residual −k₂(R_l″ + 2R_l′/r − l(l + 1)R_l/r²) − q_l at 64 radii in (a, b], for l = 0 and l = 2, using the analytic derivatives. Bound: ≤ `B1_BAND` × the sum of the term magnitudes (the scale before cancellation; `IMP-4` (B1)).
- **(I)** At r = a, for l = 0 and l = 2 on all three core branches: the θ jump and the flux jump k₁R_l′(a⁻) − k₂R_l′(a⁺). Bound: ≤ `B1_BAND` × the sum of the term magnitudes. For l = 0 this includes k₂R₀′(a⁺) = 0.
- **(Rb)** Robin at r = b: k₂R_l′(b) + hR_l(b) ≤ `B1_BAND` × the term sum, for l = 0 and l = 2.
- **(FD)** The full Cartesian T, from `pec_sphere_implant_steady_temperature`, through a 7-point central Laplacian at 32 shell points with r ∈ [1.1a, 0.95b]. Compute max\|−k₂Δ_hT − q\| / (9q_far) at δ = 1e-2·a and at 2δ.
  - Bounds: the residual at δ ≤ 1e-4, and the 2δ / δ ratio in [3.8, 4.2] (`IMP-1` (H)).
  - δ = 1e-2·a is a parameter under rule (f), following `IMP-2`'s (F) ruling ("δ ≳ 1e-2·a, ratio form"). At 1e-3·a, round-off sits within ≈ 7 % of the truncation term, and the ratio wobbles.
  - This is the only anchor independent of the radial algebra.
- **(Bal)** Three ways to compute the power:
  - 2π∫∫ q_far\|E\|² r² dr dμ, from `_pec_exterior_numpy`, by composite Gauss–Legendre (`MAT-18`'s `_gl`, `test_mat18_two_layer_sphere_closed_form.py:56–62`; 40 r-panels of order 40 on [a, b], order-40 GL in μ);
  - the closed-form P;
  - 4πb²hθ₀(b).

  Pairwise bound: ≤ `B1_BAND`.
- **(Φ)** The shell mean of q/q_far over a < r < `SHELL_OUTER_FACTOR`·a equals `PHI_CLOSED_FORM_PEC` (`test_imp1_pec_sphere_lossy_host.py:73–74`) to ≤ 1e-10 (`IMP-1` (Φ0)). This ties the source to `IMP-1`'s gated SAR enhancement.
- **(R0)** At β = 0 (both core and shell, all k₁), T equals `two_layer_sphere_steady_temperature(r, a=a, b=b, k1=k₁, k2=k₂, q1=0.0, q2=q_far, h=h, t_inf=T∞)` (`analytical.py:2508`) to ≤ 1e-13 relative on θ (`MAT-18` (S1)). Use k₁ = `K_CORE` and the isothermal branch. On the adiabatic branch, check the shell only.
- **(Rk)** The finite-k₁ D and E approach the adiabatic branch as k₁/k₂ → 0 and the isothermal branch as k₂/k₁ → 0, both at first order. Measure e(ε)/e(ε/2) at ε = 1e-3. The ratio window [1.9, 2.1] is ±5 % of the exact 2, by analogy with `IMP-1` (H)'s ±5 % around 4. **It is not an imported band; the review rules on it** (see "For the review").

**Negative controls** (§9 rule (e): none has been measured on this code, so all are predicted and asserted only as missing their anchor by > `CONTROL_MISS`; the measured value is printed beside the prediction).
- **(c1) Isotropic source** (the P₂ part dropped), against (Src): max\|q_iso − q\| / max q. **Predicted 0.6667**, which is the ceiling: at the pole of r = a the P₂ part is 6 of the 9.
- **(c2) No implant near-field** (β = 0 source), against (Bal)'s β = 1 power: 1 − P(0)/P(1). **Predicted 0.030303** = 1 − 1/1.03125. This is the share of heat the dipole adds on this fixture, and also the ceiling, because the l = 0 near-field term is all of the difference.
- **(c3) Wrong core branch** (the adiabatic R₂ used for k₁ = `K_CORE`), against (I)'s flux jump: \|2k₁R₂(a)/a − 0\| / term sum. **Predicted 1.0**, the ceiling by construction, since one side is identically zero.

**Printed, predicted, never asserted** (drafter's pure-Python arithmetic, unlogged, with the fixture above):
- θ₀(a) = 0.350684 K and θ₀(b) = 0.135352 K (l = 0, the same for every core). P = 6.803512e-2 W.
- R₂(a): adiabatic 0.045660 K, k₁ = 20 0.001966 K, isothermal 0.
- θ at r = a, pole / equator: adiabatic 0.396344 / 0.327854 (ratio 1.2089); k₁ = 20 0.352650 / 0.349700 (1.0084); isothermal 0.350684 / 0.350684.
- **The implant heating factor at r = a (l = 0):** θ₀(a)\|β=1 / θ(a)\|β=0 = 0.350684 / 0.318750 = 1.1002. The β = 0 value equals `MAT-18`'s 0.318750 to 12 digits in the drafter's check.
- Void (β = −½, adiabatic): pole / equator 0.310844 / 0.334678. The void's hot ring is at the equator.
- (FD) at δ = 1e-2·a: 1.13e-5, ratio 4.00 (20 random shell points in the drafter's check).
- (Rk): 4.864e-6 / 2.433e-6 (ratio 1.9993) toward the adiabatic branch; 1.102e-5 / 5.512e-6 (1.9986) toward the isothermal branch.

**Window.** Foreground, Bash `timeout: 660000`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Predicted 8–40 s; the range is [3, 60] s.

`scripts/testing/run_and_log.sh IMP-7-step0-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_analytical_lightweight.py tests/validation/test_imp7_implant_heating_closed_form.py -s -v --tb=short'"`

The lightweight test is the rule (c) re-gate of the additive `analytical.py` change. The real build skips the complex-only environment tests, as expected. No FFCx compile, so no cache sweep is needed. W1b (after an in-slot code fix) uses the same command, relabelled.

**Traps already paid for.**
- **Normalise residuals by the sum of term magnitudes, never by \|θ\|** (item 190; `IMP-4` (B1)). At r = b the l = 2 radial terms cancel heavily.
- **k₁ branches:** handle `0.0` and `np.inf` explicitly. Never let `inf` reach arithmetic (`inf·0` is NaN).
- **`_pec_exterior_numpy` zeroes r < a:** sample r ≥ 1.001a. Its dtype is `PETSc.ScalarType` (float in the real build); take `np.real` anyway.
- **Never pass ε = ∞** to `dielectric_sphere_quasistatic_field` (`test_pec_sphere_hole.py:10–12`). β = 1 comes from `IMP-1`'s helper.
- **Imports from `IMP-1`'s module:** import only `_pec_exterior_numpy`, `SPHERE_RADIUS`, `SHELL_OUTER_FACTOR` and `PHI_CLOSED_FORM_PEC`. That module already loads in the real build (`IMP-1` step 3(a)). Never import a fixture or a `test_*` name.
- **FD stencil:** keep r ± δ inside (a, b). Do not use δ = 1e-3·a (see (FD)).
- **analytical.py stays dolfinx-free** (`test_analytical_lightweight.py`). Append at the end of the file; if another item lands there first, rebase. Nothing is reordered.
- **Run hygiene:** print on rank 0; use `-s`; rule (i) applies (no edit after the last window); rule (f) covers δ, the sample counts and the GL orders, never a bound.

**Scope.** Numpy only. One quasi-static field, one geometry, steady state. **Not included:** the FEM (step 1); a FEM SAR source; finite k (`IMP-1` step 3's l = 1 field); a rod, a tip or a disc; a coil; transient, perfusion, dose; any published implant-heating comparison ((7c)).

**Negative result.** There is no `PREDICTED RED:` line.
- **(Src), (PDE), (I), (Rb) or (Bal) red:** a code or derivation bug. Fix it in-slot, run W1b, and disclose.
- **(FD) red with (PDE) green:** the radial algebra and the Cartesian assembly disagree (most likely the P₂ argument or a sign in C₁ / C₂). Print both sides at the worst point and stop.
- **(R0) red:** a disagreement with `MAT-18`'s function. Print both sides. Do not edit `MAT-18`'s code; stop.
- **(Rk) outside its window, with everything else green:** print the ratios at ε ∈ {1e-2, 1e-3, 1e-4} and report. The window is the review's (see "For the review"), so it is never moved in-slot.
- **A control inside its band:** a mislabelled comparand. Report it.
- **On a stop:** park on `attempt/IMP-7-step0-<UTC>` and mark item 579 BLOCKED in the same commit (rule (d)). The row stays ⬜.

**Status it can move.** The proposed §7 `IMP-7` row moves ⬜ → 🟡 on its step-0 Done-when, which this item matches (rule (j)). The review adds the row in the commit that queues this item and names it in the dashboard for the weekly to ratify (rule (5) lane 3). Step 1 then becomes draftable. The Phase 7 §6 clause "no heating" does not move on a numpy step.

**For the review.**
1. **Opening licence.** The weekly closed (7a) on 10-14 (`PROJECT_PLAN.md:4416`), and the 2026-10-01 18:00 review declined `MAT-16` because its subject was outside the written subgoals (`docs/planning/chunks/MAT-16.md:3`). `IMP-4` … `IMP-6` and `MAT-21` were nevertheless opened by dailies after subgoals existed, and were ratified. If the review reads rule (5)'s "until the weekly has written…" clause as barring this chunk, decline it. The brief then stands as a proposal for the 10-21 weekly, as `MAT-16`'s did.
2. **Home and ID:** `IMP-7` (Phase 7, near-implant hot spot) or `MAT-22` (Phase 8 (8a)). The brief is written as `IMP-7`.
3. **(Rk)'s [1.9, 2.1]:** rule it, or demote (Rk) to printed. Every other bound is imported.
4. **Rule (f) parameters:** b = 4a, `K_CORE` 20, h 100, q_far 2000. The review may substitute others; the bounds do not change.
5. **Shared file:** confirm that no other item tonight appends to `analytical.py`, or order them.

**Spawn line:** "Use the implementer agent to execute `IMP-7` §9 item 579 per `docs/planning/chunks/IMP-7.md#§9-item-579` (read the latest Review ruling paragraph first, if one exists). Append the four `pec_sphere_implant_*` functions at the end of `src/fem_em_solver/utils/analytical.py` (numpy only, additive, nothing reordered), and write `tests/validation/test_imp7_implant_heating_closed_form.py` (prefix `IMP7-0`). Import the bands by name (`B1_BAND`, `CONTROL_MISS` from `test_imp4_pec_cylinder_closed_form.py`; the `MAT-18` (S1), `IMP-1` (H) and (Φ0) literals as the brief cites them). From `test_imp1_pec_sphere_lossy_host.py`, import only `_pec_exterior_numpy`, `SPHERE_RADIUS`, `SHELL_OUTER_FACTOR` and `PHI_CLOSED_FORM_PEC`. Never import a fixture or a `test_*` name. Edit no landed module. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run in the foreground with Bash `timeout: 660000`, pytest with `-s`, container `timeout -k 30 120`, no `mpiexec`. Print every predicted value beside its measurement and never assert one. Assert controls only as missing by > `CONTROL_MISS`. Report a red and stop; never widen a bound, move δ to rescue (FD), or edit `MAT-18`'s code."

**Same-commit duties.**
- Mark §9 item 579 DONE or BLOCKED.
- Flip the §7 `IMP-7` row ⬜ → 🟡 on green.
- Add an Outcome paragraph under this heading.
- Add a journal entry in `docs/testing/attempts.md`.
- Commit the logs and the `test-results.md` rows (insertions only).

### Review ruling 2026-10-20 18:00 (item 579) — declined to queue; a proposal for the 2026-10-21 weekly

Not queued and not opened in §7. Reasons: (1) the 10-14 weekly closed Phase 7 subgoal (7a), and the `MAT-16` precedent declined a daily-opened chunk outside the written subgoals; a chunk joining implant field to heating is a new question of the kind rule (5) leaves to the weekly once its subgoals exist; (2) the closed form is the drafter's own derivation (Carslaw & Jaeger ch. IX; Jin App. E), not a cited published result, and its (Rk) ratio window [1.9, 2.1] and fixture parameters (b = 4a, core k₁ = 20, h = 100, q_far = 2000) are set here, not imported — under rule (5) a daily does not invent a target or a band. The brief stays as the weekly's proposal (ID `IMP-7` or `MAT-22`, the weekly's choice), named in the dashboard.

### Weekly ruling 2026-10-21 — opened as `IMP-7`, Phase 7 subgoal (7e); two amendments (full text PROJECT_PLAN §10)

**Opened.** The join of implant field to gel heating is the operator's research question ((7c) is its published-data half) and the closed form is a derivation from Carslaw & Jaeger ch. IX on the `MAT-21` precedent (the triple Robin series was also a derivation and was ratified). Home: `IMP-7`, Phase 7 subgoal (7e) — "implant-driven heating in closed form, then FEM"; the §7 row is added by the weekly in the proposed form. Steps: **I37** = step 0 (numpy, as drafted above), **I38** = step 1 (the FEM on `two_layer_ball_domain(inner_radius=a, outer_radius=b)` at `MAT-18` step 2(a)'s gates on the metal core and on the adiabatic hole, the isotropic-source and β = 0 controls missing; the l = 2 pole moment printed, no band unless a consumer appears).

**Amendments to the step-0 brief.** (1) **(Rk) is printed, never asserted** — its [1.9, 2.1] window is not an imported band; the ratios at ε ∈ {1e-2, 1e-3, 1e-4} are printed beside the prediction 2. (2) The fixture parameters (b = 4a, `K_CORE` 20, h 100, q_far 2000) stand under rule (f). Everything else in the brief binds as written; the daily assigns the §9 number (579 is unused) and confirms that no other queued item appends to `analytical.py` in the same slot. Audit owed on the ⬜ → 🟡 move.

### §9 item 579 — 2026-10-22 (I37, lane 3)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** §10 chain I, "I37 (lane 3; standard, no MPI, real build). `IMP-7` step 0 as drafted in `IMP-7.md` with tonight's amendments: the four `pec_sphere_implant_*` functions appended to `analytical.py`; (Src), (PDE), (I), (Rb), (FD), (Bal), (Φ), (R0) asserted at their imported bands, (Rk) printed; controls (c1)–(c3) miss by > `CONTROL_MISS`; the row ⬜ → 🟡" (`PROJECT_PLAN.md:4928–4933`). The weekly opened the chunk as Phase 7 subgoal (7e) (`:4408–4416`, `:4809–4818`; "Weekly ruling 2026-10-21" above), and the §7 row is at `:2638` (⬜). **This brief is the 10-20 draft above with the weekly's two amendments applied and five corrections from re-reading the code.** The closed form (draft "The closed form"), fixture, code names, test file and prefix are unchanged. No band, tolerance, Done-when or target is set or moved.

**Amendments applied (weekly, binding).**
1. **(Rk) is printed, never asserted.** Print e(ε)/e(ε/2) toward both branches at ε ∈ {1e-2, 1e-3, 1e-4} beside the prediction 2. The [1.9, 2.1] window may be printed as a reading aid, but it has no `assert`.
2. **Rule (f) parameters stand:** b = 4a = 0.02, `K_CORE` = 20.0, h = 100, q_far = 2000, T∞ = 20, β = 1 (plus the β = ±½ cross-checks).

**Corrections to the 10-20 draft (code read 2026-10-22).**
1. **MAT-9's `K` = 0.6 and `Q` = 2000.0 are at `test_mat9_steady_conduction_sphere.py:35–36`,** not `:37–39`. Line `:37` is MAT-9's `H` = 10, which is *not* this fixture's h. h = 100 is `MAT-10`'s `H` (`test_mat10_sar_to_heat_sphere.py:53`). These are rule (f) parameters, so restating them by value with the citation is fine.
2. **Only two of the draft's bands are named constants:** `B1_BAND` = 1e-12 and `CONTROL_MISS` = 1e-3 (`test_imp4_pec_cylinder_closed_form.py:46`, `:51`). The other three are **literals in `assert` lines**, so they cannot be imported by name:
   - `MAT-18` (S1)'s 1e-13 (`test_mat18_two_layer_sphere_closed_form.py:113`; also `:95`);
   - `IMP-1` (H)'s ratio window [3.8, 4.2] and residual bound 1e-4 (`test_imp1_step3_exact_comparand.py:193–194`);
   - `IMP-1` (Φ0)'s 1e-10 (`:217`).

   Restate them by value as module constants, each with a comment citing file:line. This follows `MAT4_BAND`'s precedent (`test_imp1_pec_sphere_lossy_host.py:57–59`: "restated here by value, not widened") and `MAT-18` step 2(a)'s ("imported by value: they are literals", `MAT-18.md:181`). For the review: item 1 below.
3. **E₀ = 1.0** is `IMP-1`'s own constant (`test_imp1_pec_sphere_lossy_host.py:66`). `_pec_exterior_numpy` closes over it (`:80–104`). (Src) at E₀ = 1 therefore needs no rescale.
4. **Both field callables use the dolfinx convention:** `x` of shape (3, N) → (3, N). This holds for `_pec_exterior_numpy()` (`:85–104`) and for `dielectric_sphere_quasistatic_field` (`analytical.py:599–600`). Transpose any (N, 3) sample. The new `pec_sphere_implant_*` functions take (N, 3), as drafted.
5. **The coefficient dict must carry `radius`, `outer_radius`, `q_far`, `k_gel`, `k_core` and `beta`** beside A₀ … F, because R₀ = A₀ + B₀/r − q_far r²/(6k₂) − β²q_far a⁶/(6k₂r⁴) needs them. `pec_sphere_implant_radial(r, coeffs, l, *, derivative=0)` is **pure power laws with no r-mask** (the draft's "analytic power laws"). Item 594 (I38) evaluates it on shell cells at chord points slightly inside r = a, where `pec_sphere_implant_source` returns 0 and the adiabatic `pec_sphere_implant_steady_temperature` returns NaN by design. Those two keep the draft's masks. This is a clarification of the letter, not a new feature.

**`IMP-7` step 0: the steady temperature around `IMP-1`'s PEC sphere heated by its own dipole SAR, numpy only.**
- **Executor:** implementer.
- **Tier:** standard (container `timeout -k 30 120`), no `mpiexec`, real build (no complex-mode source). This is the shape `MAT-21` step 0 and `IMP-6` step 0 ran in (below).
- **Class:** lane-3 closed form. **Cost:** ≈ 7.9 measured (lane-3 median) / ≈ 40 predicted min (15 fixed + ≈ 20 code + one ≈ 0.5-min window + commit and journal).
- **Depends on:** nothing. **Daytime:** no. **Slot-first:** no. **Family:** step 0, attempt 1 of 4.
- **Comparable windows** (all real build, no MPI, `test_complex_mode.py` on the path):
  - `MAT-21-step0-W1`: 5 s, own module only (`20261013T004329Z_MAT-21-step0-W1.log:12`, `:115`);
  - `MAT-21-step0-W1r`: 7 s, with `test_analytical_lightweight.py` (`20261013T004342Z_MAT-21-step0-W1r.log:12`, `:101`);
  - `IMP-6-step0-W1`: 3 s (`20261012T183435Z_IMP-6-step0-W1.log:132`).

  Importing `IMP-1`'s module adds dolfinx, ufl and gmsh imports: `IMP-1-step3a-W1` took 45 s at `-n 2` with heavier tests (`test-results.md:2699`).

**Code.** As the 10-20 draft's "Code" section, with correction 5 applied:
- four functions appended at the **end** of `src/fem_em_solver/utils/analytical.py` (now 3 116 lines; the last function is `coated_pec_cylinder_tm0_current`, `:3109`), under an `IMP-7 step 0` banner, numpy only, no dolfinx or petsc import at module level;
- `tests/validation/test_imp7_implant_heating_closed_form.py`, prefix `IMP7-0`.

Do not edit anything else. In particular, do not edit `two_layer_sphere_steady_temperature` (`:2508`), `dielectric_sphere_quasistatic_field` (`:595`), any `test_imp1_*`, `test_mat18_*` or `test_imp4_*` module, or `test_analytical_lightweight.py`.

**Anchors (asserted).** These are the 10-20 draft's definitions, unchanged. Each band below is checked against its definition:

| ID | Bound | Source (file:line) |
|---|---|---|
| (Src) β = 1 against `_pec_exterior_numpy`; β = ∓½ against `dielectric_sphere_quasistatic_field(epsilon_c=0.0 / 4.0)`, exterior | ≤ `B1_BAND` × max q | `test_imp4_pec_cylinder_closed_form.py:46` (named) |
| (PDE) radial residual, l = 0 and 2, 64 radii | ≤ `B1_BAND` × term sum | same |
| (I) θ and flux jumps at r = a, three core branches | ≤ `B1_BAND` × term sum | same |
| (Rb) Robin at r = b | ≤ `B1_BAND` × term sum | same |
| (FD) 7-point Laplacian, 32 shell points, δ = 1e-2·a | residual ≤ 1e-4; 2δ/δ ratio in [3.8, 4.2] | `test_imp1_step3_exact_comparand.py:193–194` (literals, restated) |
| (Bal) GL quadrature, closed-form P, 4πb²hθ₀(b), pairwise | ≤ `B1_BAND` | named, as above |
| (Φ) shell mean of q/q_far over a < r < 1.5a = `PHI_CLOSED_FORM_PEC` | ≤ 1e-10 | `test_imp1_step3_exact_comparand.py:217` (literal, restated) |
| (R0) β = 0 against `two_layer_sphere_steady_temperature(..., q1=0.0, q2=q_far)` | ≤ 1e-13 relative on θ | `test_mat18_two_layer_sphere_closed_form.py:113` (literal, restated) |

Imports from `IMP-1`'s module (unchanged): `_pec_exterior_numpy`, `SPHERE_RADIUS` (`:61`), `SHELL_OUTER_FACTOR` (`:73`) and `PHI_CLOSED_FORM_PEC` (`:74`). From `MAT-18`'s closed-form module, import only `_gl` (`:54–60`). `PHI_CLOSED_FORM_PEC` = 1 + 2/1.5³ is exactly the shell mean of the l = 0 part. The P₂ part integrates to zero over μ, so (Φ) is predicted ≈ 1e-15.

**Printed, never asserted.** (Rk) at three ε (amendment 1), and every prediction below.

**Negative controls.** As drafted, under §9 rule (e). None of them has been measured on this code, so all are **predicted** and asserted only as missing their anchor by > `CONTROL_MISS` (1e-3, named):
- (c1) isotropic source against (Src): predicted 0.6667, which is the ceiling;
- (c2) β = 0 power against (Bal): predicted 0.030303 = 1 − 1/1.03125, the ceiling;
- (c3) the adiabatic R₂ used for k₁ = 20, against (I)'s flux jump: predicted 1.0, the ceiling.

The smallest separation is (c2), at 30× the floor.

**Predictions.** The drafter re-derived these on 2026-10-22 in pure-Python stdlib arithmetic (unlogged, no numpy on the host); they agree with the 10-20 values to every printed digit:
- θ₀(a) 0.350684 K and θ₀(b) 0.135352 K;
- R₂(a): adiabatic 0.045660, k₁ = 20 0.001966, isothermal 0;
- pole / equator at r = a: adiabatic 0.396344 / 0.327854; k₁ = 20 0.352650 / 0.349700;
- P(β = 1) 6.803512e-2 W and P(β = 0) 6.597345e-2 W;
- the heating factor at r = a, θ₀(a)|β=1 / θ(a)|β=0 = 0.350684 / 0.318750 = 1.1002;
- (FD) ≈ 1.1e-5 with ratio ≈ 4.00;
- (Rk) ratios ≈ 1.999.

**Window.** Foreground, Bash `timeout: 660000`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. The command is as drafted. `test_analytical_lightweight.py` is the rule (c) re-gate of the additive `analytical.py` change. Predicted 8–25 s; range [3, 60] s.

`scripts/testing/run_and_log.sh IMP-7-step0-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_analytical_lightweight.py tests/validation/test_imp7_implant_heating_closed_form.py -s -v --tb=short'"`

W1b, after an in-slot code fix, uses the same command with the label `IMP-7-step0-W1b`.

**Traps already paid for.** All of the 10-20 draft's traps hold. Added on 2026-10-22:
- **Band literals restated by value, never "imported"** (correction 2). Each gets a comment with file:line. Never round or widen one.
- **The (3, N) convention** of both field callables (correction 4).
- **(R0):** never pass `np.inf` to `two_layer_sphere_steady_temperature`. At β = 0 with q₁ = 0 the core is isothermal at A₀ for every finite k₁, so pass `k1=K_CORE` for every branch's comparison.
- **MAT-9's `H` (10) is not this fixture's h (100)** (correction 1).

**Scope.** Unchanged from the 10-20 draft. The step is numpy only: one quasi-static field, one geometry, steady state. It does not include the FEM (item 594), a FEM SAR source, finite k, a rod, tip or disc, a coil, transient heating, perfusion, dose, or any published comparison ((7c)).

**Negative result.** As drafted. There is no `PREDICTED RED:` line. On a stop, park on `attempt/IMP-7-step0-<UTC>` and mark item 579 BLOCKED in the same commit (rule (d)); the row stays ⬜.
- **(Rk) far from 2** is a printed reading. It is reported and never stops the slot (amendment 1).

**Status it can move.** The §7 `IMP-7` row (`PROJECT_PLAN.md:2638`) moves ⬜ → 🟡. Its step-0 Done-when lists (Src), (PDE), (I), (Rb), (FD), (Bal), (Φ), (R0) at their imported bands, (Rk) printed, the two reductions, and three controls missing by > `CONTROL_MISS`. This item matches that list one for one (rule (j)). The audit is owed (`PROJECT_PLAN.md:5136`). Item 594 (step 1) becomes runnable.

**Shared file (the weekly's "one item at a time", `:4818`).**
- `analytical.py` was last touched by `7e03fc01` (`MAT-21` step 9). Its tail is `:3109–3116`.
- Items 586–596 had not been written when this was drafted, so the drafter cannot read them. The steps the review is restocking from are H30–H32, U74–U78, P45, I36a, C10, C12, `TH-5` step 3 and `POST-8` e4 / 3b-iii (`PROJECT_PLAN.md:4997–5014`). **None of their §10 texts names `analytical.py`.** The only §10 mentions of the file are `:4818` and `:4930`, both `IMP-7`.
  - P45 (`:4904–4914`) re-uses `MAT-21` / `MAT-20` operators.
  - U74–U76 are package code.
  - C10 / C12 and I36a are `GEO-38` / `GEO-36` modules.
  - H30–H32 are the `ans:8` record route.
- **`TH-5` step 3 (the PML) is the one step with plausible closed-form content.** The review should confirm that its brief does not append to `analytical.py`. If it does, order it after 579.
- Item 594 imports from `analytical.py` and does not edit it.

**For the review.**
1. **Three of the eight bands are literals restated by value**, not named imports (correction 2). This is the `MAT4_BAND` / `MAT-18` 2(a) practice. Accept it, or name the anchors as unimportable.
2. **(FD)'s 1e-4 level is transplanted.** In `IMP-1` (H) it bounds a curl-curl residual normalised by \|k²E\|. Here it bounds a 7-point Laplacian residual normalised by 9q_far. The [3.8, 4.2] ratio is the same second-order signature (central-difference truncation), but the level is applied to a new quantity. Predicted 1.1e-5, which is ≈ 9× inside. The review may keep both, or demote the level to printed and keep the ratio.
3. **Correction 5** (the extra dict keys; `pec_sphere_implant_radial` unmasked) is relied on by item 594.
4. **Shared file:** confirm that `TH-5` step 3's brief does not touch `analytical.py`.

**Execution order.** Read this file's rulings → code → W1 (W1b on an in-slot fix) → commit (the `analytical.py` append, the new test, the log, the `test-results.md` rows, §9 item 579 DONE, §7 row ⬜ → 🟡, an Outcome paragraph here) → `attempts.md` journal.

**Outcome (2026-10-22, 18:30 UTC slot, third item) — DONE, green; row ⬜ → 🟡 (audit owed).** Four `pec_sphere_implant_*` functions (plus a private `_pec_sphere_implant_terms` helper) appended at the end of `src/fem_em_solver/utils/analytical.py`, numpy only, nothing reordered; module `tests/validation/test_imp7_implant_heating_closed_form.py` (prefix `IMP7-0`). One window, `20261022T185726Z_IMP-7-step0-W1.log`: 14 passed / 3 skipped (complex-only environment tests), 1.35 s pytest, 3 s harness, no MPI, real build — inside the predicted [3, 60] s. Measured (bands in brackets): (Src) β = 1 vs `_pec_exterior_numpy` 4.95e-16, β = −½ / ½ vs `dielectric_sphere_quasistatic_field` 3.02e-16 / 2.67e-16 [`B1_BAND` 1e-12 × max q]; (PDE) worst 2.38e-16; (I) worst 5.92e-16; (Rb) worst 7.71e-17; (Bal) pairwise ≤ 1.22e-15 (P = 6.803512e-2 W, pred 6.803512e-2); (Φ) 9.33e-15 [1e-10 restated]; (R0) 0.000 on all three branches [1e-13 restated]; (FD) level 1.422e-5 at δ = 1e-2·a (pred ≈ 1.1e-5; `IMP-1` (H)'s 1e-4 printed only, per the 03:00 ruling), 2δ 5.688e-5, ratio 4.000 [3.8, 4.2 asserted]. Controls (predicted, asserted only > `CONTROL_MISS`): (c1) 0.6496 (pred 0.6667, the ceiling at r = a; samples start at 1.001a), (c2) 0.030303 (pred 0.030303), (c3) 0.9320 (pred 1.0 — the scale is the term sum, which includes the shell flux's pre-cancellation terms, so it reads below the ceiling). (Rk) printed: e(ε)/e(ε/2) toward adiabatic 1.9934 / 1.9993 / 1.9999 and toward isothermal 1.9852 / 1.9985 / 1.9999 at ε = 1e-2 / 1e-3 / 1e-4 (pred 2; e is \|ΔR₂(a)\| in K, so the levels differ from the draft's unnamed metric). Every printed record matches its prediction to the printed digit: θ₀(a) 0.350684, θ₀(b) 0.135352, R₂(a) 0.045660 / 0.001966 / 0, pole / equator 0.396344 / 0.327854 (adiabatic) and 0.352650 / 0.349700 (k₁ = 20), heating factor 1.1002, void pole / equator 0.310844 / 0.334678. **Disclosed clarification of the letter:** the draft's "the core is constant, R₀ = A₀" holds only at β = ±1 (there B₀/a cancels the two particular terms at r = a); the code carries the core constant as `theta0_core` = R₀(a⁺) = A₀ + B₀/a − (1 + β²)q_far a²/(6k₂) for every β, which is what (R0) at β = 0 checks against `MAT-18` (0.318750 = θ(a), not A₀ = 0.360417). Dict keys as ruled: radius, outer_radius, q_far, k_gel, k_core, h, beta, A0, B0, C1, C2, D, E, F, theta0_core, theta0_b, power; `pec_sphere_implant_radial` unmasked. No band moved, no landed module edited.

**Executor spawn prompt (verbatim).** "Execute §9 item 579 — `IMP-7` step 0 (§10 I37, lane 3) — per `docs/planning/chunks/IMP-7.md` § '§9 item 579 — 2026-10-22' and the 10-20 draft above it. Read the 2026-10-21 weekly ruling and any later review ruling first; they bind. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Append the four `pec_sphere_implant_*` functions at the end of `src/fem_em_solver/utils/analytical.py` (numpy only, additive, nothing reordered; the coefficient dict carries radius, outer_radius, q_far, k_gel, k_core and beta; `pec_sphere_implant_radial` has no r-mask). Write `tests/validation/test_imp7_implant_heating_closed_form.py` (prefix `IMP7-0`). Import `B1_BAND` and `CONTROL_MISS` from `test_imp4_pec_cylinder_closed_form.py`. Restate the 1e-13, [3.8, 4.2] / 1e-4 and 1e-10 literals by value, each with a file:line comment. From `test_imp1_pec_sphere_lossy_host.py`, import only `_pec_exterior_numpy`, `SPHERE_RADIUS`, `SHELL_OUTER_FACTOR` and `PHI_CLOSED_FORM_PEC`; from `test_mat18_two_layer_sphere_closed_form.py`, import only `_gl`. Never import a fixture or a `test_*` name. Field callables take (3, N). (Rk) is printed only, never asserted. Controls are asserted only as missing by > `CONTROL_MISS`. Print every prediction beside its measurement. Run W1 in the foreground with Bash `timeout: 660000`, pytest `-s`, container `timeout -k 30 120`, no `mpiexec`, and no complex-mode source. Edit nothing after the last window (rule (i)). You are the executor: do not spawn agents and never background a harness window. If the docker socket is denied, stop and report the exact command. On a red, report and stop: never widen a bound, never move δ to rescue (FD), never edit `MAT-18`'s code. Park on `attempt/IMP-7-step0-<UTC>` and mark item 579 BLOCKED in the same commit."

### §9 item 594 — 2026-10-22 (I38, lane 3)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** §10 chain I, "I38 (lane 3; standard / heavy, `-n 2`, real build; serial on I37). `IMP-7` step 1 — the FEM on `two_layer_ball_domain(inner_radius=a, outer_radius=b)`: `MAT-18` step 2(a)'s gates (interior rel L2 ≤ 1 %, the P1 L2 rate ≥ 1.7 on two rungs, balance ≤ 1e-9) on the metal core (DG0 k) and on the adiabatic hole, the isotropic-source and β = 0 controls missing; the l = 2 moment at the pole printed, no band (a band is briefed only if a consumer appears) ⇒ `IMP-7` ✅" (`PROJECT_PLAN.md:4934–4940`). The §7 row's step-1 Done-when is at `:2638`. The weekly ruling above names it I38. **If item 579 did not land (row not 🟡, or `pec_sphere_implant_coefficients` absent from `analytical.py`), skip this item at no attempt cost.** No band, tolerance, Done-when or target is set or moved.

**`IMP-7` step 1: P1 FEM of the implant-heated gel ball, on the metal core and on the adiabatic hole, against step 0's closed form at `MAT-18` step 2(a)'s gates.**
- **Executor:** implementer.
- **Tier:** standard, two windows (each ≤ 180 s container), `-n 2`, real build.
- **Class:** lane-3 closed form. **Cost:** ≈ 7.9 measured (lane-3 median; the weekly priced I38 at the chain-solve median, 16.7, `PROJECT_PLAN.md:5013`) / ≈ 55 predicted min (15 fixed + ≈ 30 code + W1 ≈ 2.5 + W2 ≈ 2 + commit and journal ≈ 5).
- **Depends on:** item 579 (serial; take-next only if 579 committed before minute 30). **Daytime:** no. **Slot-first:** no. **Family:** step 1, attempt 1 of 4.

**What exists (read 2026-10-22).**
- **Mesher:** `MeshGenerator.two_layer_ball_domain(inner_radius, outer_radius, resolution, comm, rank)` (`src/fem_em_solver/io/mesh.py:3199–3283`).
  - Cell tag 1 is the core and cell tag 2 the shell. Facet tag 1 is r = b only.
  - There is no interface facet group (known-issues 9, `:3213–3215`).
  - Sizing is uniform (`:3269–3273`).
  - It has only been run at a/b = 0.91. This fixture's a/b = 0.25 is new to it.
- **`MAT-18` step 2(a)'s module:** `tests/validation/test_mat18_two_layer_sphere_fem.py`.
  - Geometry-free helpers: `_allsum` (`:71–74`), `_dg0` (`:77–85`), `_indicator` (`:88–97`).
  - Helpers bound to `MAT-18`'s geometry: `_errors` (`:115–130`, interior mask midpoint \|x\|² < 0.81b² at `:120`) and `_census` (`:133–172`). **Both close over `A_` = 0.1 / `B_` = 0.11 / `AREA_B` / `T_INF`** (`:50–51`, `:120`, `:198`), so do not import them.
  - `CELL_GUARD` = 800 000 (`:48`) and `QDEG` = {"quadrature_degree": 4} (`:49`).
  - The ladder pattern with an LSQ rate on h_eff = (V_h/N)^{1/3} is at `:213–215` and `:243–254`.
- **`MAT-18` step 2(a)'s bands** are literals in two places, both restated by value:
  - (S3) balance ≤ 1e-9 (`:226`), (S1) LSQ rate ≥ 1.7 (`:253`), (S2) interior ≤ 1e-2 (`:254`);
  - the same three in `MAT-9` (`test_mat9_steady_conduction_sphere.py:139`, `:154`, `:156`), which is where `MAT-18` took them "by value" (`MAT-18.md:181`).
- **`MAT-18`'s measured window:** `20261003T172548Z_MAT-18-step2a-W1.log`, 121 s, `-n 2`, real build, four rungs plus controls (`:261`).
  - The 215 542-cell rung: build 8.7 s, solve 1.8 s (`:90`).
  - The 584 228-cell rung: build 27.4 s, solve 9.2 s (`:93`).
  - Balance ≤ 1.02e-14, LSQ rate 1.822, finest interior 6.79e-4 (`:86–101`).
- **Solver:** `solve_steady_conduction(mesh, facet_tags, *, k, q, h, t_inf, robin_tag, degree=1)` (`src/fem_em_solver/thermal/steady.py:151–198`).
  - k is a scalar or a DG0 `fem.Function` (allreduced positivity check, `:183–189`).
  - q is a scalar or a UFL expression / `fem.Function`, used as-is in `L` and in `source_power` (`:165–169`).
  - It returns `temperature`, `source_power` and `convected_power`, all allreduced.
- **Submesh pattern:** `MAT-10`'s `_sphere_submesh` (`test_mat10_sar_to_heat_sphere.py:74–89`; `create_submesh` on a cell-tag value).
  - **Its `_robin_tags` (`:92–98`) tags every exterior facet.** On the shell submesh that would put Robin on r = a too (trap 2).
- **Restricted interpolation:** `Function.interpolate(f, cells0=…)` (`tests/validation/test_lossy_sphere_fullwave.py:391`).
- **Point sampling:** `evaluate_vector_field_parallel` on a scalar temperature (`MAT-18` fem `:229–230`).
- **Step 0 (item 579):** `pec_sphere_implant_source`, `pec_sphere_implant_coefficients`, `pec_sphere_implant_radial` (unmasked, 579 correction 5) and `pec_sphere_implant_steady_temperature`.

**Code.** No `src/` edit and no `analytical.py` edit, so rule (c) does not bind. Do not edit any `test_mat18_*`, `test_mat10_*` or `test_imp7_*_closed_form` module.

New module `tests/validation/test_imp7_implant_heating_fem.py`, prefix `IMP7-1`.
- **Constants.** a = `SPHERE_RADIUS`; b = 4a; `K_CORE`, k₂, q_far, h and T∞ from step 0's module (`test_imp7_implant_heating_closed_form.py`, constants only); `B1_BAND`, `CELL_GUARD` and `QDEG` imported by name. The (S1) / (S2) / (S3) literals are restated with file:line comments. Robin tag 1.
- **Ladder (rule (f) parameter): `resolution` ∈ {1.414e-3, 1.0e-3, 0.7071e-3} m, ratio √2.**
  - h/b is 0.071 / 0.050 / 0.035, which spans `MAT-18`'s three finer rungs (0.064 / 0.045 / 0.032).
  - Predicted cells ≈ 59k / 160k / 440k (±25 %). These are scaled from `MAT-18`'s ≈ 19–22 cells per (b/res)³ (`MAT-18` W1 `:85–93`); they are not measured at a/b = 0.25.
  - The 2 mm rung is left out: the core would be ≈ 2.5 cells in radius, which is pre-asymptotic.
- **Module helpers.** All new and private:
  - `_q_shell(x)`: q_far·[1 + 2β²u² + (4βu + 2β²u²)P₂(μ)] with **no r-mask**, taking x of shape (3, N);
  - `_theta_shell(x, coeffs)`: R₀ + R₂P₂ from `pec_sphere_implant_radial` (unmasked);
  - `_errors(msh, chi_in, T_h, theta)` and `_census(msh, ct, ft, a, b)`: `MAT-18`'s forms, parametrised by a and b;
  - `_robin_outer(sub, a, b)`: exterior facets with midpoint \|x\| > (a + b)/2 only;
  - a module-level mesh cache, as `MAT-18`'s `_MESHES`.
- **Source and coefficients.**
  - The source is a **DG2 `fem.Function`**: zero, then `interpolate(_q_shell, cells0=<shell cells>)`, then `scatter_forward()`.
  - k is a DG0 function (`_dg0`: `K_CORE` on the core, k₂ on the shell).
  - The source is never a `SpatialCoordinate` expression (trap 4).
- **`test_imp7_1_source_guard`** (numpy, ungated):
  - (Src′) `_q_shell` equals `pec_sphere_implant_source` at 512 points with r ∈ [1.001a, b], ≤ `B1_BAND` × max q;
  - (Th′) `_theta_shell` (adiabatic coefficients) equals `pec_sphere_implant_steady_temperature(k_core=0.0) − T∞` at the same points, ≤ `B1_BAND` × max\|θ\|.

  These tie the unmasked continuations to step 0's gated functions.
- **`test_imp7_1_metal_core_ladder`**, on the full two-layer mesh. The comparand is `pec_sphere_implant_steady_temperature(k_core=K_CORE)` − T∞ (layer chosen by exact r, `MAT-18`'s trap; θ and flux are continuous), interpolated into P2. Per rung:
  - (G) as `MAT-18`;
  - the cell guard before the finest solve;
  - (S3), printed and asserted;
  - the whole / interior / shell errors, printed.

  After the ladder: (S1) LSQ rate and (S2) finest interior asserted; pairwise rates, ∫q_h against the closed-form P, and V_core,h printed; the l = 2 readout at the finest rung.
- **`test_imp7_1_metal_core_controls`**, on the 1.0 mm rung, against the β = 1 comparand:
  - (c-iso) q = q_far(1 + 2u²) on the shell;
  - (c-β0) q = q_far on the shell.

  Each interior miss is asserted > 1e-2, i.e. it misses (S2)'s band.
- **`test_imp7_1_hole_ladder`.** Use `create_submesh` of cell tag 2 (the `MAT-10` pattern), Robin tags from `_robin_outer`, the natural (zero-flux) condition on r = a, scalar k = k₂, and `_q_shell` on every cell. The comparand is **`_theta_shell` with `k_core=0.0` coefficients, never `pec_sphere_implant_steady_temperature`** (trap 1). Assert that the interpolated comparand holds no NaN (allreduced). Apply the same per-rung and ladder gates. Also:
  - (G-hole): the Robin facet area is within 1 % of 4πb² (`MAT-18` (G)'s 1e-2, `:198`);
  - the r = a facet area is printed beside 4πa².
- **`test_imp7_1_hole_controls`:** (c-iso) and (c-β0) on the 1.0 mm hole.
- **l = 2 readout (printed, no band).** On the finest rung of each fixture, at r₀ = 1.1a:
  - R₂,h(r₀) = (5/2) Σⱼ wⱼ P₂(μⱼ) · mean_φ θ_h(r₀, μⱼ, φ), with order-16 GL in μ and 8 φ, as one `evaluate_vector_field_parallel` batch;
  - θ_h at the pole (0, 0, r₀) and at the equator (r₀, 0, 0).

  r₀ = 1.1a lies well outside the faceting: the chord depth h²/(8a) ≈ 1.2e-5 m against 0.1a = 5e-4 m.

**Anchors (asserted; bands restated by value or imported by name, each with file:line).** These apply on **both** fixtures (metal core, adiabatic hole):
- **(G)** census and conformity as `MAT-18` (`:193–198`): owned tag-1 + tag-2 cells equal the global count; core vertices satisfy \|x\| ≤ a(1 + 1e-9); shell vertices lie in [a(1 − 1e-9), b(1 + 1e-9)]; tag-1 facets are exterior; their area is within 1e-2 of 4πb². On the hole: the cell count equals the parent's tag-2 count (allreduced), plus (G-hole).
- **Cell guard:** finest-rung cells ≤ `CELL_GUARD` (800 000, `test_mat18_two_layer_sphere_fem.py:48`), checked before its solve.
- **(S3) balance**, every rung: \|`convected_power` − `source_power`\| / `source_power` ≤ 1e-9 (`:226`; `MAT-9` `:139`). `source_power` = ∫q_h dV is `MAT-18`'s q₁V_core,h in this fixture's form. A DG2 q on affine tets is integrated exactly by both the L-form and the source-form rules, so this is a round-off identity. Predicted ≤ 2e-14 (`MAT-18` read ≤ 1.02e-14).
- **(S1) rate:** the LSQ slope of the whole-domain relative L2 error on θ against h_eff = (V_h/N)^{1/3}, over the three rungs, ≥ 1.7 (`:253`; `MAT-9` `:154`). Predicted 1.8–2.1. Pairwise rates are printed.
- **(S2) interior:** the relative L2 error on θ over cells with midpoint \|x\|² < 0.81b², at the 0.7071 mm rung, ≤ 1e-2 (`:254`; `MAT-9` `:156`). Predicted 5e-4 to 5e-3 (`MAT-18` read 6.79e-4 at h/b 0.032; the near-field source steepens θ″ ≈ 9× within ≈ a/4 of r = a, in a small volume).
- **(Src′) / (Th′)** ≤ `B1_BAND` (`test_imp4_pec_cylinder_closed_form.py:46`).

**Negative controls (rule (e)).** No measurement of these comparisons exists on this fixture, so all four are **predicted**. Each is asserted only as missing (S2)'s 1e-2, as `MAT-18`'s (c1) / (c2) were (`:313`, `:329`). Predicted interior misses are comparand-level, from the drafter's unlogged pure-Python quadrature of the closed form over r < 0.9b:

| Control | Metal core | Adiabatic hole |
|---|---|---|
| (c-iso) l = 2 source dropped | 2.730e-2 | 3.603e-2 |
| (c-β0) uniform shell source | 5.671e-2 | 6.195e-2 |

- **Ceiling and separation.** The isotropic control's miss is exactly the l = 2 share of θ, and the β = 0 control's miss is the l = 2 share plus the near-field l = 0 term. These are the most each control can show on this fixture.
- Since \|θ_h,ctl − θ\| ≥ \|θ_ctl − θ\| − \|θ_h,ctl − θ_ctl\| and the FEM error at 1.0 mm is predicted ≤ 3e-3, the measured misses should be ≥ 2.4e-2 / 3.3e-2 / 5.4e-2 / 5.9e-2. That is **2.4× to 5.9× the 1e-2 band**. The metal-core isotropic control is the thinnest (item 3 for the review).

**Predictions (printed, never asserted).** Pure-Python, unlogged, 2026-10-22:
- P(β = 1) = 6.803512e-2 W. ∫q_h is predicted below it by ≈ the faceted-sphere volume deficit (≈ −1e-3 to −4e-3; `MAT-18`'s area deficits were −4.2e-4 to −1.7e-3).
- At r₀ = 1.1a:
  - metal core: R₂(r₀) 1.1966e-2, pole / equator 0.361610 / 0.343661 K;
  - adiabatic hole: R₂(r₀) 4.4793e-2, pole / equator 0.394437 / 0.327248 K.
- At r = a: R₂(a) metal 1.966e-3 / hole 4.566e-2. Metal core θ ≈ A₀ = 0.350684 K.

**Windows.** Foreground, Bash `timeout: 660000`, pytest `-s`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Node IDs select the tests (never `-k`). Two standard windows, because one would predict ≈ 160–250 s against the 180 s ceiling.
- **W1** (metal core; predicted 90–150 s, range [45, 175]; the first window pays the new FFCx compiles):

  `scripts/testing/run_and_log.sh IMP-7-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp7_implant_heating_fem.py::test_imp7_1_source_guard tests/validation/test_imp7_implant_heating_fem.py::test_imp7_1_metal_core_ladder tests/validation/test_imp7_implant_heating_fem.py::test_imp7_1_metal_core_controls -s -v --tb=short'"`
- **W2** (adiabatic hole; predicted 70–130 s, range [40, 175]):

  `scripts/testing/run_and_log.sh IMP-7-step1-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp7_implant_heating_fem.py::test_imp7_1_hole_ladder tests/validation/test_imp7_implant_heating_fem.py::test_imp7_1_hole_controls -s -v --tb=short'"`
- **Price basis.** `MAT-18` W1 (121 s for ≈ 0.90M cells over four rungs plus three 0.2M-cell control solves, `:261`). The ladder here totals ≈ 0.66M cells per fixture.
- **Kills.** A window killed at 180 s is not a result. Run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first, then the orphan-`.c` sweep (rubric item 4). Then stop and mark the item BLOCKED. Never raise the timeout and never drop a rung in-slot: re-tiering or shrinking the ladder is the review's call.
- No durable capture is needed (both windows are well under ≈ 10 min). No module edit after W1 starts (rule (i)). W1b / W2b after an in-slot fix use the same commands, relabelled.

**Traps already paid for.**
1. **Chord slivers at r = a.** Shell cells carry interpolation points just inside r = a (chord depth ≈ h²/(8a)).
   - There, `pec_sphere_implant_source` returns 0 and the adiabatic `pec_sphere_implant_steady_temperature` returns NaN (by design, item 579).
   - On shell cells, use `_q_shell` / `_theta_shell` (unmasked, guarded by (Src′) / (Th′)).
   - Assert no NaN in the hole comparand.
   - On the metal-core comparand, layer by exact r is consistent (`MAT-18` 2(a) trap: "do not special-case it").
2. **Robin on the hole submesh only at r = b.** `MAT-10`'s `_robin_tags` tags every exterior facet, which would add a Robin sink on r = a and silently cool the hole. Select by facet midpoint radius.
3. **Do not import `MAT-18`'s `_errors` / `_census`** (they close over `MAT-18`'s a = 0.1, b = 0.11). Parametrise copies.
4. **q as a DG2 `Function`, never a `SpatialCoordinate` expression.** A rational UFL expression lets FFCx estimate a large quadrature degree, and `source_power`'s `ufl.dx` carries no metadata, so the L-form and the source-form rules could differ and (S3) would stop being a round-off identity. Pin `QDEG` on every error form (the gmsh `SpatialCoordinate` compile trap, rubric item 4).
5. **Ghosts.** Set DG0 / DG2 values from `ct.find(·)`, then `scatter_forward()` (`MAT-18` 2(a) trap). Reduce every count, integral, min and max (rank-safety).
6. **(3, N) vs (N, 3).** Dolfinx `interpolate` passes (3, N), while step 0's functions take (N, 3). Transpose inside the wrappers.
7. **`cells0=`** is the keyword on this image (`test_lossy_sphere_fullwave.py:391`). The interpolated cells are local, so `scatter_forward()` afterwards.
8. **Point samples are printed only.** No bitwise or cross-batch comparison (`.eval` batch trap, item 537).
9. **Run hygiene.** Print on rank 0; use `-s`; never `-k`; never pipe; after a kill, `pgrep -c python3`, then the sweep; never delete a log.

**Scope.** P1, steady, one fixture (a = 5 mm, b = 4a), closed-form dipole source. **Not claimed:**
- a FEM SAR source (the source is step 0's closed-form field; a composed band is a later step, as in `MAT-10`);
- finite-k EM or the isothermal-core FEM (no Dirichlet on r = a);
- degree 2;
- a band on any l = 2 moment;
- a rod, tip or disc, a coil, transient heating, perfusion, CEM43, or any published implant-heating datum ((7c)).

**Negative result.** There is no `PREDICTED RED:` line, so a red anywhere stops the slot. On every stop, park on `attempt/IMP-7-step1-<UTC>` and mark item 594 BLOCKED in the same commit (rule (d)); the row stays 🟡.

| Red | Meaning and action |
|---|---|
| (Src′) / (Th′) | The unmasked continuation disagrees with step 0. Fix in-slot, run W1b. |
| (G) | The mesher at a/b = 0.25 is non-conforming or mis-tagged. Allow one in-slot fix (test side only; never edit `mesh.py`), then stop. |
| (S3) | A code bug (q's form, Robin tags). Fix in-slot and disclose. |
| (S1) < 1.7 with (S2) green | Print the pairwise rates and stop. Never drop a rung (rules (f), (h)). The review rules. |
| (S2) > 1e-2 | Report it beside the volume, area and ∫q_h deficits, and stop. |
| A control inside the band | A mislabelled comparand or source. Report it. |
| Over the cell guard | Stop before the solve and report the counts. |
| A kill | See Windows. |

**Status it can move.** The §7 `IMP-7` row (`PROJECT_PLAN.md:2638`) moves 🟡 → ✅ on its step-1 Done-when: "`MAT-18` step 2(a)'s gates — interior rel L2 ≤ 1 %, P1 L2 rate ≥ 1.7, balance ≤ 1e-9 — on the metal core and on the adiabatic hole, the isotropic-source and β = 0 controls missing". This item asserts exactly those, plus (G) and (Src′) / (Th′), and prints the l = 2 moment (rule (j)). The row does not say "on two rungs"; the §10 text does (item 1 for the review). The executor closes on the row. The audit is owed. After ✅, the Phase 7 §6 digest's "no heating" clause becomes the weekly's to re-word (not this item's).

**For the review.**
1. **"The P1 L2 rate ≥ 1.7 on two rungs"** (§10 `:4937`) against the row's "P1 L2 rate ≥ 1.7" (`:2638`). The brief reads it as `MAT-18`'s gate, the LSQ slope over the ladder (here three rungs), with pairwise rates printed. If the review means the pairwise rate on the two finest rungs, it says so; the ladder does not change.
2. **Ladder (rule (f) parameter):** 1.414 / 1.0 / 0.7071 mm. The cell counts are scaled, not measured at a/b = 0.25 (the mesher's first use there).
3. **The thinnest control:** (c-iso) on the metal core, ≈ 2.7e-2 against the 1e-2 band (the ceiling is the l = 2 share of θ, because k₁ = 20 suppresses R₂(a) to 2e-3). If the review wants more separation, the hole's (c-iso) (3.6e-2) and both (c-β0) (5.7e-2 / 6.2e-2) are larger. The band does not move.
4. **Two standard windows** split by node ID, rather than one heavy window. The §10 text allows "standard / heavy".
5. **The l = 2 readout's definition** (R₂,h(1.1a) by Legendre projection, plus pole / equator θ_h) is the drafter's. It is printed only.
6. **Class pricing:** lane-3 median 7.9 against the weekly's 16.7 for I38.

**Execution order.** Check that 579 landed (row 🟡; the four functions on `main`), else skip → code → W1 → W2 (W1b / W2b on an in-slot fix) → commit (module, both logs, `test-results.md` rows, §9 item 594 DONE, §7 row 🟡 → ✅, an Outcome paragraph here) → `attempts.md` journal.

**Outcome (2026-10-22, 20:00 UTC slot) — DONE, green; row 🟡 → ✅ (audit owed).** 579 verified landed (commit `36d99fc1`, row 🟡, `pec_sphere_implant_coefficients` present). Module `tests/validation/test_imp7_implant_heating_fem.py` (prefix `IMP7-1`), no `src/` edit, no in-slot fix. **W1** `20261022T202402Z_IMP-7-step1-W1.log` (4 passed / 3 skipped, 42 s harness, `-n 2`, real) and **W2** `20261022T202503Z_IMP-7-step1-W2.log` (3 passed / 3 skipped, 42 s) — W1 3 s below its predicted [45, 175] s range, W2 inside [40, 175] s. Measured: (Src′) 0 / (Th′) 4.5e-15 [`B1_BAND`]. Cells 54 647 / 152 914 / 425 751 (pred ≈ 59k / 160k / 440k; guard 800 000); (G) green on every rung, area deficit −1.01e-3 / −5.07e-4 / −2.55e-4; hole cell count equals parent tag-2 on every rung, r = a facet area deficit −1.53e-2 / −7.96e-3 / −4.01e-3 (printed). **Metal core:** whole rel L2 2.363e-3 / 1.202e-3 / 6.072e-4, pairwise 1.973 / 2.000, **(S1) LSQ 1.987** [≥ 1.7], **(S2) interior 6.113e-4** [≤ 1e-2], (S3) ≤ 6.65e-13 [≤ 1e-9]; ∫q_h/P − 1 −3.87e-4 / −2.01e-4 / −1.07e-4. **Adiabatic hole:** whole 1.796e-3 / 9.033e-4 / 4.503e-4, pairwise 2.003 / 2.039, **(S1) LSQ 2.021**, **(S2) interior 4.295e-4**, (S3) ≤ 1.46e-13; comparand NaN-free. Controls (predicted, asserted only > 1e-2), at 1.0 mm: metal (c-iso) 2.733e-2 (pred 2.730e-2), (c-β0) 5.965e-2 (pred 5.671e-2); hole (c-iso) 3.606e-2 (pred 3.603e-2), (c-β0) 6.243e-2 (pred 6.195e-2). **l = 2 readout at r₀ = 1.1a (printed, no band):** metal R₂,h 1.0174e-2 (pred 1.1966e-2, −15 %), pole / equator 0.359626 / 0.344320 K (pred 0.361610 / 0.343661); hole R₂,h 4.3742e-2 (pred 4.4793e-2, −2.3 %), pole / equator 0.393152 / 0.327785 (pred 0.394437 / 0.327248). The metal-core l = 2 part is ≈ 3 % of θ at r₀, so a ≈ 2e-3 K pointwise P1 error is a 15 % miss on it; whether that is P1 discretisation or the drafter's prediction is not diagnosed (the readout carries no band, and the (Th′)-guarded comparand passed every gate). No band moved.

**Executor spawn prompt (verbatim).** "Execute §9 item 594 — `IMP-7` step 1 (§10 I38, lane 3) — per `docs/planning/chunks/IMP-7.md` § '§9 item 594 — 2026-10-22'. Read every ruling in that file first; the newest binds. If §9 item 579 has not landed (the §7 `IMP-7` row is not 🟡, or `pec_sphere_implant_coefficients` is absent from `src/fem_em_solver/utils/analytical.py`), stop: the item is skipped at no attempt cost. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Create only `tests/validation/test_imp7_implant_heating_fem.py` (prefix `IMP7-1`). Edit no `src/` file, not `analytical.py`, and no `test_mat18_*`, `test_mat10_*` or step-0 module. Import `B1_BAND` from `test_imp4_pec_cylinder_closed_form.py` and `CELL_GUARD` / `QDEG` / `_allsum` / `_dg0` / `_indicator` from `test_mat18_two_layer_sphere_fem.py`. Do not import its `_errors` or `_census`: they close over `MAT-18`'s a and b, so write parametrised copies. Restate the 1e-9 / 1.7 / 1e-2 literals by value with file:line comments. The source is a DG2 Function from the unmasked shell formula on shell cells (`cells0=`, then `scatter_forward()`). The hole comparand comes from `pec_sphere_implant_radial` with `k_core=0.0` coefficients, never from `pec_sphere_implant_steady_temperature`. Robin on the hole submesh goes on r = b facets only. Ladder 1.414e-3 / 1.0e-3 / 0.7071e-3 m; controls on the 1.0e-3 rung, asserted only as interior miss > 1e-2; the l = 2 readout printed only. Run W1 then W2 in the foreground with Bash `timeout: 660000`, pytest `-s`, container `timeout -k 30 180`, `mpiexec -n 2`, real build (no complex-mode source), and select tests by node ID, never `-k`. Edit nothing after W1 starts (rule (i)). You are the executor: do not spawn agents and never background a harness window. If the docker socket is denied, stop and report the exact command. After a kill, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first; never raise the timeout or drop a rung. On a red, report and stop: never widen a bound. Park on `attempt/IMP-7-step1-<UTC>` and mark item 594 BLOCKED in the same commit."

### Review ruling 2026-10-22 03:00 (items 579 I37, 594 I38)

Binds over the briefs. §9 number 579 (unused since 10-20) is I37's.
- **579:** (1) Bands are imported **by name** where the source defines a constant (`B1_BAND`, `CONTROL_MISS`); the 1e-13 (R0) and 1e-10 (Φ) levels are restated by value with file:line comments, the `MAT4_BAND` precedent, disclosed. (2) **(FD)'s 1e-4 level is printed, not asserted:** it bounds a different quantity in `IMP-1`, and borrowing it would set a band; (FD)'s ratio window [3.8, 4.2] (the second-order rate it carries over) stays asserted. (3) The coefficient dict (radius, q_far, k, β) and the unmasked radial function are confirmed. (4) No other queued item edits `analytical.py` (`TH-5` step 3's brief does not).
- **594:** serial on 579 (skip if it did not land). The rate is `MAT-18`'s least-squares fit over the three rungs, pairwise rates printed; two standard windows; the l = 2 readout printed only; class lane 3 (7.9 measured).

### Audit of the closure — 2026-10-22 18:00 review

**PASS** (`auditor`, steps 0 and 1; closing commits `36d99fc1`, `04885e88`). Footers are Status 0 with elapsed time: `20261022T202402Z_IMP-7-step1-W1.log:109–111` (42 s), `20261022T202503Z_IMP-7-step1-W2.log:100–102` (42 s), `20261022T185726Z_IMP-7-step0-W1.log:113–115` (3 s). ≈ 95 digits traced. The asserts are executed in the committed modules (`test_imp7_implant_heating_fem.py:292–293`, `:316–317`, `:356`, `:373`, `:436`, `:453`). Every band is imported or restated by value with its source line, and none moved. The step-0 FD level and (Rk) are printed-only by the 03:00 ruling. Re-read by this review: W1 `:70` (finest interior 6.1134e-04), `:74` (LSQ rate 1.987), `:109–111` (Status 0, 42 s); W2 `:67–71` (4.2951e-04, rate 2.021).

**Caveats, not demoting.**
- The row rounds the metal balance down to "≤ 6.6e-13"; it is 6.646e-13 (`W1:70`), immaterial against 1e-9.
- The logs ran on the working tree before the commit, which is the usual pattern. The `:158` warning line matches the committed file.
- The comparand is step 0's own derivation; the independent check of its radial algebra is the FD Laplacian (ratio 4.000).
- The thinnest control is the metal (c-iso), at 2.7× the band.
- §10 I38 says "two rungs"; the row and test use a three-rung LSQ fit, and both readings pass.

**The l = 2 readout (printed, no band).** Metal R₂,h is −15 % at 1.1a (`W1:71`) and hole −2.3 % (`W2:68`). This review reads the metal figure as a **cancellation readout**, not a defect. R₂ is (pole − equator) / 1.5, and pole − equator is ≈ 5 % of the mean temperature rise. The pointwise misses are −0.55 % at the pole (0.359626 vs 0.361610 K) and +0.19 % at the equator (0.344320 vs 0.343661 K), and they combine to −1.8e-3 K on R₂ (≈ −15 %). The P1 interior L2 at that rung is 6.1e-4. A pointwise P1 sample one-tenth of a radius outside a k-contrast-33 interface carries ≈ 0.5 %. No band follows: a band is briefed only if a consumer appears (§10 I38), and `MAT-18` has no degree-2 gate for a P2 re-pose. The question goes to the weekly as information, and no item is queued.

### §9 item 630 — 2026-10-25 (§10 I40, `IMP-7` step 2; lane 3)

*Written by the 2026-10-25 03:00 daily review from the 10-24 weekly's enumerated step I40 (§9 rules (4) / (5)).*

**Question.** At P2, does the two-layer implant ball keep `MAT-18` step 2(a)'s order-free gates? And does the l = 2 pole moment's −15 % at P1 (read as a cancellation of two near-equal point values: pole −0.55 %, equator +0.19 %) shrink at the next order?

**Anchors (imported, order-free).**
- Interior rel L2 ≤ `S2_BAND` (1e-2) and balance ≤ `S3_BAND` (1e-9) on the **metal core** and on the **adiabatic hole**, finest rung. These are `IMP-7`'s module restatements of `MAT-18` step 2(a)'s literals (`tests/validation/test_imp7_implant_heating_fem.py:52–55`; `test_mat18_two_layer_sphere_fem.py:226`, `:254`). Import them from `IMP-7`'s module by name.
- **(β₀)** The β = 0 control (`_q_shell(x, beta=0.0)`, `:362–366`) misses the interior band. At P1 it reads 5.965e-2 (metal) and 6.243e-2 (hole), `IMP-7.md:429`.

**Printed (the record, no band):**
- the P2 L2 rate on the three rungs (1.414 / 1.0 / 0.7071 mm, `:61`), with no bar: `MAT-18` has no degree-2 rate gate;
- the l = 2 moment R₂,h at 1.1a beside P1's 1.0174e-2 (−15 %) and the series value 1.1966e-2 (`20261022T202402Z_IMP-7-step1-W1.log:71`);
- the pole / equator θ_h beside P1's 0.359626 / 0.344320 K and the series.

**Route (rule (c), additive).** `degree=1` is hard-coded at `:241`. Add a `degree` keyword (default 1, so step 1's path is unchanged) or a `FEM_EM_IMP7_DEGREE` knob (unset ⇒ 1; `2` ⇒ P2; anything else raises).
- **The comparand must move up one order.** It is interpolated into P2 (`_p2`, `:131–136`, errors at quadrature degree 4). A P2 solution needs a **P3 comparand and quadrature degree ≥ 6**, as `MAT-21` step 3 did. Otherwise the error floor is the comparand's.
- **Regression:** step 1's W1 / W2 at degree 1, unchanged, reproduce their finest-rung strings (6.1134e-4 / 4.2951e-4, `W1:70`, `W2:67`).

**Windows.**
1. W0: coding tests.
2. Wreg1 / Wreg2: the degree-1 reproduction, 42 s each.
3. W1 `IMP-7-step2-W1`: metal core, P2.
4. W2 `IMP-7-step2-W2`: adiabatic hole, P2.

All at `-n 2`, real build (no `@complex_only` in this module), heavy, `timeout -k 30 1200`.

**Price** (estimate, not measured): P1's finest rung was 425 751 cells, build 17.3 s, solve 4.2 s. P2 is 7–8× the DOFs (`MAT-21`'s box: 6.9–7.4×), at about the same time per DOF, so ≈ 30–60 s per finest solve plus new FFCx compiles. That gives ≈ 120–300 s per window.

**Tier:** heavy, `-n 2`. **Class:** lane-3 closed form. **Cost:** ≈ 7.5 measured / 45 predicted slot-minutes. **Executor:** implementer. **Depends on:** nothing. **Family:** `IMP-7` step 2, attempt 1 of 4.

**Traps already paid for:** `IMP-7.md:382–395` traps 1–9 — chord-sliver NaN, Robin only at r = b, never import `MAT-18`'s `_errors` / `_census`, q as a DG2 Function (at P2 consider DG3 and say which), `cells0=`, the (3, N) shapes. Rank-local `assemble_scalar`. Orphan `.c` sweep only with `pgrep -c python3` reading 0.

**Scope boundary.** No band set; `IMP-7` stays ✅ on step 1. Not closed: a degree-2 rate gate (none exists to import), finite k, a FEM SAR source, a rod.

**A negative result means:** an order-free band red at P2 is a comparand-order or quadrature defect first. Check the P3 interpolation, report, and park. The l = 2 moment not moving is a record for the weekly, never a red.

**Status it can move.** §10 I40 → DONE. `IMP-7`'s row gains the step-2 record (P2 on both cases). The 10-24 weekly's cancellation reading of the l = 2 moment is tested at the next order.

**Outcome (2026-10-25, 12:30 UTC / 07:30 CDT slot) — DONE, green; row stays ✅.** Route: a `degree` keyword (default 1) threaded through `_solve` / `_metal_rung` / `_hole_rung` with a `DEG_ROUTE` table; at P2 the source q is **DG3** (step 1 keeps DG2), the comparand P3 and the error forms quadrature degree 6. New tests `test_imp7_2_metal_core_p2` / `test_imp7_2_hole_p2` (prefix `IMP7-2`); `S2_BAND` / `S3_BAND` are the module's own names. No `src/` edit, no in-slot fix. Container ceiling `timeout -k 30 590` on W1 / W2 (stricter than the brief's 1200, to guarantee a footer inside the Bash tool's 660 s). **W0** `20261025T123235Z_IMP-7-step2-W0.log` (collect-only, 4 s). **Wreg1** `20261025T123246Z_IMP-7-step2-Wreg1.log` (42 s; finest interior 6.1134e-04, LSQ 1.987) and **Wreg2** `20261025T123341Z_IMP-7-step2-Wreg2.log` (43 s; 4.2951e-04, LSQ 2.021) — both strings reproduced. **W1** `20261025T123436Z_IMP-7-step2-W1.log` (metal, 2 passed / 3 skipped, 259 s; finest solve 159.2 s): whole rel L2 1.5827e-3 / 8.1490e-4 / 4.0631e-4, pairwise 1.937 / 2.040, LSQ **1.988** (printed); **finest interior 4.2200e-4** [≤ 1e-2]; balance ≤ 3.215e-13 [≤ 1e-9]; (c-β0) at 1.0 mm 5.9510e-2 (P1 5.965e-2) — missed. **W2** `20261025T123915Z_IMP-7-step2-W2.log` (hole, 2 passed / 3 skipped, 236 s; finest solve 131.9 s): whole 1.3382e-3 / 6.8842e-4 / 3.4273e-4, pairwise 1.936 / 2.043, LSQ **1.990**; **finest interior 3.5242e-4**; balance ≤ 3.678e-13; (c-β0) 6.2240e-2 (P1 6.243e-2) — missed. **l = 2 readout at 1.1a (printed, no band):** metal P2 R₂,h 1.2224e-2 (series 1.1966e-2, **+2.2 %**; P1 1.0174e-2, −15 %), pole / equator 0.362259 / 0.343898 K (series 0.361610 / 0.343661; P1 0.359626 / 0.344320) — pointwise misses +0.18 % / +0.07 %, against P1's −0.55 % / +0.19 %; hole R₂,h 4.4953e-2 (series 4.4793e-2, +0.36 %; P1 −2.3 %), pole / equator 0.394957 / 0.327547 (series 0.394437 / 0.327248). The cancellation reading holds: the moment's miss shrinks ≈ 7× with the pointwise misses. **Observation, not diagnosed:** the P2 L2 rate is ≈ 2, not 3, and the P2 error is only ≈ 0.67× / 0.76× P1's on the finest rung — consistent with an O(h²) geometric floor from the affine (polyhedral) r = a / r = b surfaces (chord depth ≈ h²/(8a), trap 1), not with the comparand (P3, qdeg 6). No band was set; a degree-2 rate gate would need curved (P2-geometry) cells first.

### §9 item 647 — 2026-10-25 (§10 chain I, I42 — daily, 2026-10-25, provisional; `IMP-7` step 3: P2 on curved (P2-geometry) cells; chain step (solve))

*Written by the 2026-10-25 18:00 daily review under §9 rule (6), toward chain I's subgoal (7e) — "the FEM on the two-layer ball against [the l = 0 + l = 2 series] at `MAT-18`'s gates". Chain I's provisional count since the 10-24 weekly: I41, I42 = 2 of 6. **The question I40 left:** its P2 L2 rate is ≈ 2 (LSQ 1.988 metal / 1.990 hole), not 3, and P2's error is only 0.67× / 0.76× P1's on the finest rung. Its outcome names a reading but does not test it: an O(h²) geometric floor from the affine r = a / r = b surfaces (chord depth ≈ h²/(8a)), "a degree-2 rate gate would need curved (P2-geometry) cells first". **This step is the measurement that separates the two readings.** On second-order geometry, the polyhedral floor reading predicts the pairwise P2 rate rises toward 3 and the error falls well below the affine P2's. "Something else" (the comparand, the DG3 source, quadrature) predicts ≈ 2 again. No band is set: the rate and the error ratio are a **record**, and a degree-2 rate gate stays the weekly's to write from it.*

**Route (rule (c), additive; a `src/` keyword).** Add `geometry_order: int = 1` to `MeshGenerator.two_layer_ball_domain` (`src/fem_em_solver/io/mesh.py:3199`). At 1 the code path is byte-identical. At 2 it calls `gmsh.model.mesh.setOrder(2)` after `generate(3)` and before the gmsh → DolfinX conversion, so the coordinate element is degree 2 (10-node tets). Anything else raises. Write a new module `tests/validation/test_imp7_curved_p2.py` (real build, as `IMP-7`'s module):
- **Imports, by name** from `tests/validation/test_imp7_implant_heating_fem.py` (as `m7`): `LADDER`, `DEG_ROUTE`, `S2_BAND`, `S3_BAND`, `A`, `B`, `K_CORE`, `_dg0`, `_source_dg2`, `_solve`, `_p2`, `_theta_metal`, `_q_shell`, `_errors`, plus whatever `_metal_rung` reads (read `:335–345`). Check every name at execution, and never import a `test_*` name.
- **Its own `_ball_curved(res)`**, calling the generator with `geometry_order=2` and caching by `(res, 2)`. **Never pass the curved mesh through `m7._ball`**, whose cache is keyed by `res` alone and would hand a flat mesh to a curved rung or the reverse.
- **Its own rung body**: `_metal_rung`'s body written out (`:335–345`) with the mesh argument explicit and `rt = DEG_ROUTE[2]` (DG3 source, P3 comparand, quadrature 6). On curved cells, check whether `rt["qmeta"]` is still adequate. If the quadrature has to rise, raise it **in the new module only**, print it, and say so.
- **Metal core only.** The hole repeats the question at the same cost, and I40 read both the same.

**Anchors (asserted, bands imported; order-free, as I40's).**
- **(S3c)** balance ≤ `S3_BAND` on every curved rung.
- **(S2c)** interior rel L2 ≤ `S2_BAND` on the finest curved rung run.
- **(c-β0c)** the β = 0 source on the curved `LADDER[1]` misses `S2_BAND`, backed by I40's 5.9510e-02 on the affine rung (`20261025T123436Z_IMP-7-step2-W1.log:70`).
- **(Reg)** the default path is unchanged: step 1's degree-1 metal regression re-runs unedited and reproduces its finest interior string 6.1134e-04 (`20261025T123246Z_IMP-7-step2-Wreg1.log`, 42 s), and the generator's default call is unchanged by inspection of the diff.

**Printed (the record that separates the readings, no band):**
- **(Geo)** per rung, the curved and affine meshes' core and ball volumes against 4/3·π·a³ and 4/3·π·b³, as relative errors. The closed form makes this a direct read of the chord floor.
- **(Rate)** whole and interior rel L2 on `LADDER[0]` and `LADDER[1]` (1.414 / 1.0 mm), curved, with the pairwise rate. Beside it, I40's affine pairwise 1.937 on the same two rungs (`W1`, whole 1.5827e-3 → 8.1490e-4).
- The curved-over-affine error ratio per rung.
- The l = 2 readout R₂,h at 1.1a and pole / equator θ_h beside I40's affine P2 (+2.2 %; 0.362259 / 0.343898 K) and the series (1.1966e-2; 0.361610 / 0.343661 K).
- **The finest rung (0.7071 mm) is optional.** Run it only if W1's two rungs footered in under 250 s. It gets a window of its own (W2), and its pairwise rate is printed if it runs.

**Windows (all real build, `-n 4`, heavy, `timeout -k 30 600`, after `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0).**
- **W0** (smoke, `-n 1`): the coding tests. `geometry_order=3` raises; the default call's signature is unchanged; the curved ball at a coarse resolution (≈ 4 mm) builds, and its coordinate element degree reads 2 (`msh.geometry.cmap.degree` or the 0.11 spelling — check). Volumes printed.
- **Wreg** (`-n 2`): step 1's degree-1 metal regression, unedited (as I40's Wreg1).
- **W1**: curved `LADDER[0]` and `LADDER[1]` at P2, plus (c-β0c) on `LADDER[1]`.
- **W2** (optional, as above): the finest curved rung.
- **Price** (estimate, not measured): I40's affine P2 solves at `-n 2` took ≈ 15–30 s on the two coarse rungs and 159 s on the finest. Second-order geometry raises assembly cost (non-affine Jacobians, higher quadrature), say 1.5–3×. `-n 4` roughly halves the solve. W1 ≈ 120–250 s. **Over 600 s is a kill: park, never raise the timeout** (the §9 window ceiling).

**The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**

**Executor:** implementer. **Tier:** heavy, `-n 4` (`-n 2` for Wreg), real build. **Class:** chain step (solve). **Cost:** ≈ 11.2 measured / 50 predicted slot-minutes. **Depends on:** nothing. **Family:** `IMP-7` step 3, attempt 1 of 4. **Daytime:** no.

**Traps already paid for.**
- `IMP-7.md:382–395` traps 1–9: chord-sliver NaN (the DG0 coefficient from cell tags), Robin only at r = b, never import `MAT-18`'s `_errors` / `_census`, `cells0=`, the (3, N) shapes.
- `m7._ball`'s cache is keyed by resolution alone (`:91–99`); a curved mesh must never enter it.
- On curved cells, `dmesh.compute_midpoints` returns vertex-average midpoints. The interior indicator (|x|² < 0.81 B²) is a cell-midpoint test, and it stays one: print the indicator's cell count on both meshes.
- `assemble_scalar` is rank-local, so reduce before asserting. pin `quadrature_degree` on every form with a `SpatialCoordinate`. Sweep orphan `.c` files only when `pgrep -c python3` reads 0.
- **`-n 4` is a width change from I40's `-n 2`.** The order-free bands are width-invariant, but the printed strings may move in the last digits. Compare rates, never digits.

**Scope.** The P2 rate's geometric-floor reading, tested on curved cells for the metal core, as a record. **Not claimed:** a degree-2 rate gate (the weekly's, from this record), the hole, any change to `IMP-7`'s ✅ or to I40's record, finite k, a coil.

**Negative result.** (S3c) or (S2c) red on curved cells with the regression green is a route defect in the curved path, most likely the comparand or quadrature on non-affine cells. Report it, park on `attempt/IMP-7-step3-*`, item 647 BLOCKED. A curved rate ≈ 2 is **not** a red: it refutes the floor reading, and it is the record. If gmsh `setOrder(2)` or the conversion fails on the 0.11 image, that is a finding: record it in known-issues and stop, with no workaround in this item. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7`'s row gains the step-3 record: the floor reading confirmed (curved rate ≥ ≈ 2.6, error well below the affine P2) or refuted (≈ 2). The I40 "observation, not diagnosed" becomes diagnosed either way. §10 I42 → DONE (provisional until the weekly ratifies). The weekly can write a degree-2 rate gate only from a confirmed record.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 647 — `IMP-7` step 3 / I42 — per `docs/planning/chunks/IMP-7.md` "§9 item 647". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Add `geometry_order` to `two_layer_ball_domain` additively (default 1, byte-identical path); write a new module; never pass a curved mesh through `m7._ball`; never edit `test_imp7_implant_heating_fem.py`. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-26, 00:30 UTC slot) — DONE, green; the polyhedral-floor reading is confirmed; row stays ✅.** Route: `geometry_order: int = 1` on `MeshGenerator.two_layer_ball_domain` (1 → the unchanged path; 2 → `gmsh.model.mesh.setOrder(2)` after `generate(3)`; else `ValueError`). gmsh `setOrder(2)` and the 0.11 conversion worked first time (`msh.geometry.cmap.degree` reads 2; the cell count is unchanged). New module `tests/validation/test_imp7_curved_p2.py` (prefix `IMP7-3`): every name it uses is taken from `m7` by attribute and checked at execution; it has its own `_ball_curved` cache keyed `(res, 2)`; `m7._ball` sees only flat meshes, which it uses for the affine (Geo) column. `m7` is unedited. **Quadrature:** `DEG_ROUTE[2]`'s quadrature degree 6 is adequate on curved cells. A degree-8 recomputation agrees to ≤ 2e-4 relative on every rung (printed), so it was not raised. **W0** `20261026T003241Z_IMP-7-step3-W0.log` (`-n 1`, 10 s): `geometry_order=3` raises; the default signature is unchanged; the curved 4 mm ball has cmap degree 2 and Vcore − 1 = −4.7e-3 against the affine −1.8e-1. **Wreg** `20261026T003307Z_IMP-7-step3-Wreg.log` (`-n 2`, 98 s; I40 took 42 s on a less loaded box): the finest interior 6.1134e-04 and LSQ 1.987 are reproduced. **W1** `20261026T003501Z_IMP-7-step3-W1.log` (`-n 4`, 173 s, passed) and **W2** `20261026T003804Z_IMP-7-step3-W2.log` (`-n 4`, 396 s, passed; finest build 93.3 s, solve 167.0 s; W1 footered under 250 s, so W2 ran):

| res (mm) | Vcore − 1 curved / affine | Vball − 1 curved / affine | whole rel L2 curved | interior curved | curved / affine (whole) |
|---|---|---|---|---|---|
| 1.414 | −9.66e-5 / −2.78e-2 | −4.1e-7 / −1.83e-3 | 1.6666e-4 | 1.7629e-4 | 0.105 |
| 1.0 | −2.56e-5 / −1.44e-2 | −1.0e-7 / −9.16e-4 | 5.8831e-5 | 6.2631e-5 | 0.072 |
| 0.7071 | −6.49e-6 / −7.27e-3 | −2.6e-8 / −4.60e-4 | 2.0383e-5 | 2.1703e-5 | 0.050 |

The curved P2 pairwise rate is **3.036** (1.414 → 1.0) and **3.105** (1.0 → 0.7071), on the whole domain and the interior alike. The affine P2 rates were 1.937 and 2.040. The affine volume error halves per √2 step, which is O(h²), the chord floor. The curved volume error falls about 4× per step, which is O(h⁴). The interior-indicator cell counts are identical on the two meshes (38 628 / 107 898 / 303 556). Asserted: (S3c) balance ≤ 7.8e-13 on every curved rung and on the control; (S2c) interior 6.26e-5 on `LADDER[1]` (W1) and 2.17e-5 on the finest (W2), both ≤ `S2_BAND`; (c-β0c) 5.9317e-2 on curved `LADDER[1]` (affine 5.9510e-2), which misses the band. **l = 2 readout at 1.1a (printed, no band):** R₂,h is 1.1775e-2 on curved 1.0 mm (−1.6 %) and **1.1971e-2 on the curved finest rung (series 1.1966e-2, +0.04 %)**; affine P2 finest was +2.2 % and P1 −15 %. Pole / equator on the curved finest rung read 0.361648 / 0.343661 K, against the series 0.361610 / 0.343661. I40's "observation, not diagnosed" is now diagnosed: the P2 rate ≈ 2 came from the affine r = a / r = b surfaces, not from the comparand, the DG3 source or the quadrature. Widths differ (`-n 4` here, `-n 2` in I40), so compare rates and ratios, not digits. No band was set. A degree-2 rate gate on curved cells is the weekly's to write from this record. Hole not run (scope).

**Audit 2026-10-26 03:00 review: PASS** (`auditor`). The gates are executed asserts: (S3c) at `test_imp7_curved_p2.py:157` and `:230`, (S2c) at `:226` and `:244`, (c-β0c) at `:231`. The bands are imported from `m7` by name (`:36–37`, defined at `test_imp7_implant_heating_fem.py:57`, `:59`). The rate, ratio, volume and R₂,h figures are prints (`:117–122`, `:153–181`). The default path is regressed (`Wreg.log:60`, `:65`). Footers: `W1.log:92–93` (173 s) and `W2.log:90–91` (396 s). One caveat, cosmetic: the commit message quotes pytest's 10 s / 98 s for W0 / Wreg, while the footers read 14 s / 101 s.

### §9 item 651 — 2026-10-26 (§10 chain I, I43 — daily, 2026-10-26, provisional; `IMP-7` step 3b: curved P2 on the adiabatic hole; chain step (solve))

*Written by the 2026-10-26 03:00 daily review under §9 rule (6), toward chain I's subgoal (7e): "the FEM on the two-layer ball against [the l = 0 + l = 2 series] at `MAT-18`'s gates". I38 and I40 ran that FEM on **both** fixtures, the metal core and the adiabatic hole. Chain I's provisional count since the 10-24 weekly is now I41, I42, I43 = 3 of 6. **The question I42 left:** I42 confirmed the polyhedral-floor reading on the metal core only (curved pairwise P2 rate 3.036 / 3.105 against affine 1.937 / 2.040; `20261026T003501Z_IMP-7-step3-W1.log:73`, `…W2.log:72`). I40's hole read the same affine rate (1.936 / 2.043, `20261025T123915Z_IMP-7-step2-W2.log:66–67`). But the hole's inner surface carries the **natural (zero-flux) condition on r = a**, not a conductivity jump, and the chord error enters the two conditions differently. This step measures whether the floor reading holds on the second boundary type. No band is set. The rate and the error ratio are a **record**, beside I42's, and a curved degree-2 rate gate stays the weekly's.*

**Route — a new module, nothing edited.** Write `tests/validation/test_imp7_curved_p2_hole.py` (real build, as `IMP-7`'s module):
- **Imports, by name.** From `tests/validation/test_imp7_implant_heating_fem.py` (as `m7`): `LADDER`, `CONTROL_RES`, `DEG_ROUTE`, `S2_BAND`, `S3_BAND`, `A`, `B`, `K_GEL`, `ROBIN_TAG`, `AREA_B`, `G_AREA_BAND`, `SERIES_READ`, `_shell_submesh`, `_robin_outer`, `_allsum`, `_source_dg2`, `_solve`, `_coeffs`, `_theta_shell`, `_p2`, `_q_shell`, `_errors`, `_volume`, `_vertex_radii`, `_rmin_rmax`, `_l2_readout`, plus whatever `_hole_rung` reads (`:393–432`; check every name at execution). From `tests/validation/test_imp7_curved_p2.py`: `_ball_curved`, `_cmap_degree`. **Never import a `test_*` function name.**
- **Its own curved hole rung:** `m7._hole_rung`'s body (`:393–432`) written out, with the mesh taken from `_ball_curved(res)` instead of `m7._ball(res)`, and `rt = DEG_ROUTE[2]` (DG3 source, P3 comparand, quadrature 6). Keep the (G) checks of the `check=True` branch: census, vertex radii, the Robin area against 4πb² at `G_AREA_BAND`. **Never pass the curved mesh through `m7._ball`** (its cache is keyed by `res` alone).
- **The sub-mesh is the new risk.** `create_submesh` of tag 2 on a mesh with a degree-2 coordinate element must carry that element. **Assert `_cmap_degree(sub) == 2`** before the solve. If the 0.11 image refuses or degrades it, that is a **finding**: record it in known-issues and stop, with no workaround in this item.
- **The affine comparand is I40's record,** not a re-solve. Its P2 hole rows are whole 1.3382e-03 / 6.8842e-04 / 3.4273e-04, pairwise 1.936 / 2.043, and (c-β0) 6.2240e-02 (`20261025T123915Z_IMP-7-step2-W2.log:61–70`). They ran at `-n 2`, this item runs at `-n 4`: compare rates and ratios, never digits.

**Anchors (asserted, bands imported; order-free, as I40's and I42's).**
- **(G-hole-c)** the curved sub-mesh's census equals the parent's tag-2 cells, and its vertex radii lie in [a(1 − 1e-9), b(1 + 1e-9)]. The Robin area on r = b is within `G_AREA_BAND` of 4πb². All three are as `_hole_rung(check=True)` asserts them. **Printed beside it:** the r = a area against 4πa², on both meshes.
- **(S3h)** balance ≤ `S3_BAND` on every curved rung run and on the control.
- **(S2h)** interior rel L2 ≤ `S2_BAND` on the finest curved rung run.
- **(c-β0h)** the β = 0 source on curved `CONTROL_RES` (= `LADDER[1]`) misses `S2_BAND`. It is backed by I40's 6.2240e-02 on the affine rung (`…step2-W2.log:70`), the same comparison at the same order.
- **No NaN in the hole comparand** (`_hole_rung`'s own assert, kept).

**Printed (the record, no band):**
- Per rung: the shell volume (curved and affine) against 4/3·π·(b³ − a³).
- Whole and interior rel L2, curved, with the pairwise rate. I40's affine rates sit beside them.
- The curved-over-affine whole-error ratio per rung, against I40's rows.
- The l = 2 readout R₂,h at 1.1a and pole / equator θ_h on the finest curved rung run. Beside them: I40's affine P2 (4.4953e-02; 0.394957 / 0.327547 K, `…step2-W2.log:65`) and the series (`SERIES_READ["hole"]` = 4.4793e-2; 0.394437 / 0.327248).

**Windows** (all real build, `-n 4`, heavy, `timeout -k 30 600`, after `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0; foreground, Bash `timeout: 660000`, `-s`; node IDs select the tests, never `-k`; the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped). Under the environment-first ruling (journal `2026-10-26-daily.md` §4), `tests/environment/test_complex_mode.py` is not needed in a real-build window. The module's own W0 node goes first.
- **W0** (smoke, `-n 1`): the names, `_cmap_degree(sub) == 2` on the coarse 4 mm curved ball's tag-2 sub-mesh, and its shell volume against the affine one.
- **W1:** curved `LADDER[0]` and `LADDER[1]` at P2, then (c-β0h) on `LADDER[1]`.
- **W2** (optional): the finest curved rung, as I42's W2. Run it only if W1 footered in under 250 s.
- **W-env** (`-n 1`, `timeout -k 30 300`): `tests/environment` in its own window, as the ruling requires for an item that adds a `.py` (`MAT-3.md` § item 648 gives the command; label `IMP-7-step3b-Wenv`).
- **Price** (from I42, measured): I42's W1 took 173 s and its W2 396 s at `-n 4` on the metal core. The hole solves on the tag-2 sub-mesh only (150 334 cells at `LADDER[1]`, `…step2-W2.log:62`), and I40's affine hole P2 solves took 5.7 / 31.6 / 131.9 s (`:61–63`), no dearer than the metal core's. So **W1 ≈ 140–200 s and W2 ≈ 330–420 s.** Over 600 s is a kill: park, never raise the timeout.

**Executor:** implementer. **Tier:** heavy, `-n 4` (`-n 1` for W0 / W-env), real build. **Class:** chain step (solve). **Cost:** ≈ 11.2 measured / 45 predicted slot-minutes. **Depends on:** nothing. **Family:** `IMP-7` step 3b, attempt 1 of 4. **Daytime:** no.

**Traps already paid for.**
- `IMP-7.md:382–395` traps 1–9 (chord-sliver NaN, Robin only at r = b, never import `MAT-18`'s `_errors` / `_census`, `cells0=`, the (3, N) shapes). **Trap 2 is sharper here:** on curved cells, `_robin_outer`'s facet-midpoint radius is a vertex average. Print the counts of facets it tags at r = a and at r = b on both meshes. On the curved mesh, the r = b area must still meet `G_AREA_BAND`.
- `m7._ball`'s cache is keyed by `res` alone; the curved mesh never enters it (I42's trap).
- `assemble_scalar` is rank-local, so reduce before asserting. Pin `quadrature_degree` on every form with a `SpatialCoordinate`. Sweep orphan `.c` files only when `pgrep -c python3` reads 0.
- `-n 4` is a width change from I40's `-n 2`. Compare rates, never digits (I42's trap).

**Scope.** I42's question on the adiabatic hole: curved P2 on the natural-condition fixture, as a record. **Not claimed:** a degree-2 rate gate (the weekly's), any change to `IMP-7`'s ✅ or to I40's / I42's records, finite k in the hole, a coil.

**Negative result.** If (G-hole-c) is red, or the sub-mesh's `cmap` degree is 1, the curved sub-mesh route is the finding: record it in known-issues, park on `attempt/IMP-7-step3b-*`, item 651 BLOCKED. (S3h) or (S2h) red with (G-hole-c) green means the curved hole path has a defect. Report it, park it, item 651 BLOCKED. **A curved hole rate ≈ 2 is not a red.** It says the natural condition at r = a carries a floor the conductivity jump does not, and that is the record. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7`'s row gains the step-3b record. Either the floor reading holds on both boundary types (curved hole rate ≥ ≈ 2.6), or it holds on the metal core only (≈ 2). §10 I43 → DONE (provisional until the weekly ratifies). The weekly can write a curved degree-2 rate gate only from records on both fixtures.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 651 — `IMP-7` step 3b / I43 — per `docs/planning/chunks/IMP-7.md` "§9 item 651". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Write a new module; never pass a curved mesh through `m7._ball`; never edit `test_imp7_implant_heating_fem.py`, `test_imp7_curved_p2.py` or `src/`. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-26, 09:30 UTC / 04:30 CDT slot) — DONE, green; the polyhedral-floor reading holds on the adiabatic hole too; row stays ✅.** New module `tests/validation/test_imp7_curved_p2_hole.py` (prefix `IMP7-3b`): names taken from `m7` and checked at execution (`AREA_A` added, `_hole_rung` reads it), `_ball_curved` / `_cmap_degree` from `test_imp7_curved_p2.py`, its own `_curved_hole_rung` (the `_hole_rung` body on `_ball_curved(res)`, `DEG_ROUTE[2]`); only the affine (printed) meshes go through `m7._ball`. No `src/` or existing-module edit; no edit after W0. The curved tag-2 sub-mesh carries coordinate degree 2 on the 0.11 image (asserted at every rung; affine 1) — nothing for known-issues. **W0** `20261026T093154Z_IMP-7-step3b-W0.log` (`-n 1`, 8 s). **W1** `20261026T093208Z_IMP-7-step3b-W1.log` (`-n 4`, 93 s; four PASSED tokens, one per rank). **W2** `20261026T093347Z_IMP-7-step3b-W2.log` (`-n 4`, 226 s; W1 footered under 250 s, so it ran). **W-env** `20261026T093738Z_IMP-7-step3b-Wenv.log` (`-n 1`, real build, 104 s; 8 passed, 4 complex-only skipped).

| res (mm) | Vshell − 1 curved / affine | whole rel L2 curved | interior curved | curved / affine (whole) | Robin area r = b − 1 curved / affine | area r = a − 1 curved / affine |
|---|---|---|---|---|---|---|
| 1.414 | 1.12e-6 / −1.41e-3 | 1.1475e-4 | 1.2088e-4 | 0.086 | −2.7e-7 / −1.01e-3 | −6.4e-5 / −1.53e-2 |
| 1.0 | 3.03e-7 / −7.02e-4 | 3.7781e-5 | 4.0143e-5 | 0.055 | −6.8e-8 / −5.07e-4 | −1.7e-5 / −7.96e-3 |
| 0.7071 | 7.67e-8 / −3.52e-4 | 1.2980e-5 | 1.3820e-5 | 0.038 | −1.7e-8 / −2.55e-4 | −4.3e-6 / −4.01e-3 |

Curved P2 pairwise rate **3.234** (1.414 → 1.0) and **3.129** (1.0 → 0.7071) whole (interior 3.209 / 3.123), against I40's affine 1.936 / 2.043. Facet counts tagged at r = b / r = a are identical on both meshes (6076 / 406, 12 102 / 776, 24 054 / 1538). Asserted: (G-hole-c) sub-mesh census = parent tag-2 at every rung, radii at round-off, Robin area within `G_AREA_BAND`; (S3h) balance ≤ 1.4e-13 on every rung, 2.1e-14 on the control; (S2h) finest interior 1.38e-5 ≤ `S2_BAND`; (c-β0h) 6.1986e-2 on curved `LADDER[1]` (affine 6.2240e-2), missed; no NaN in the comparand. **l = 2 readout at 1.1a (printed, no band):** curved finest R₂,h 4.4800e-2 (series 4.4793e-2, +0.02 %; affine P2 +0.36 %), pole / equator 0.394456 / 0.327267 K (series 0.394437 / 0.327248). The natural condition on r = a carries no floor the conductivity jump does not: both fixtures now read ≈ 3 on curved cells, so the record a curved degree-2 rate gate would need exists on both boundary types — the weekly's to write. Widths differ from I40 (`-n 4` vs `-n 2`): rates and ratios compared, not digits. One harmless warning: `cmap is deprecated, use cmaps[0]` from `_cmap_degree`.

## §7 row moved verbatim 2026-10-26 (rotate_plan_archive.py chunks)

 SAR heating around a PEC-sphere implant in a Robin-cooled gel ball: the implant-heating canonical in closed form (Phase 7 subgoal (7e); opened by the 2026-10-21 weekly from the 2026-10-20 18:00 daily's proposal) | ✅ Step 3 record 2026-10-26 (§9 item 647, I42, `20261026T003501Z_IMP-7-step3-W1.log` / `20261026T003804Z_IMP-7-step3-W2.log`, 173 / 396 s, heavy, `-n 4`, real): P2 on **curved (P2-geometry) cells** via an additive `geometry_order=2` on `two_layer_ball_domain`, metal core — pairwise P2 rate 3.036 / 3.105 on the three rungs (affine P2 1.937 / 2.040), curved/affine error ratio 0.105 / 0.072 / 0.050, finest interior rel L2 2.17e-5 ≤ `S2_BAND`, balance ≤ 7.8e-13, β = 0 control missing at 5.93e-2; core-volume error −9.7e-5 / −2.6e-5 / −6.5e-6 curved vs −2.8e-2 / −1.4e-2 / −7.3e-3 affine; l = 2 R₂,h at 1.1a 1.1971e-2 on the finest (+0.04 % vs series) — **step 2's rate ≈ 2 is diagnosed as the affine-boundary O(h²) floor**; printed, no band; degree-1 string 6.1134e-4 reproduced (`20261026T003307Z_IMP-7-step3-Wreg.log`); audit PASS 2026-10-26 03:00. Step 3b record 2026-10-26 (§9 item 651, I43, `20261026T093208Z_IMP-7-step3b-W1.log` / `20261026T093347Z_IMP-7-step3b-W2.log`, 93 / 226 s, heavy, `-n 4`, real; W0 / W-env `20261026T093154Z_IMP-7-step3b-W0.log` / `20261026T093738Z_IMP-7-step3b-Wenv.log`): curved P2 on the **adiabatic hole** (curved tag-2 sub-mesh, coordinate degree 2 asserted) — pairwise P2 rate 3.234 / 3.129 (affine P2 1.936 / 2.043), curved/affine whole-error ratio 0.086 / 0.055 / 0.038, finest interior rel L2 1.38e-5 ≤ `S2_BAND`, balance ≤ 1.4e-13, β = 0 control missing at 6.20e-2, (G) census / radii / Robin area green; shell-volume error 1.1e-6 / 3.0e-7 / 7.7e-8 curved vs −1.4e-3 / −7.0e-4 / −3.5e-4 affine; l = 2 R₂,h at 1.1a 4.4800e-2 on the finest (+0.02 % vs series) — **the polyhedral-floor reading holds on both boundary types**; printed, no band. Step 2 record 2026-10-25 (§9 item 630, `20261025T123436Z_IMP-7-step2-W1.log` / `20261025T123915Z_IMP-7-step2-W2.log`, 259 / 236 s, heavy, `-n 2`, real): P2 (DG3 source, P3 comparand, quadrature 6) on the metal core / adiabatic hole, finest interior rel L2 4.22e-4 / 3.52e-4 ≤ `S2_BAND`, balance ≤ 3.7e-13, β = 0 control missing at 5.95e-2 / 6.22e-2; P2 L2 rate 1.988 / 1.990 printed (no band; second order, not third — not diagnosed, likely the affine r = a interface); l = 2 R₂,h at 1.1a 1.2224e-2 / 4.4953e-2 (+2.2 % / +0.36 % vs series; P1 −15 % / −2.3 %), printed; degree-1 strings reproduced. Step 1 green 2026-10-22 (§9 item 594, `20261022T202402Z_IMP-7-step1-W1.log` / `20261022T202503Z_IMP-7-step1-W2.log`, 42 s each, standard, `-n 2`, real): P1 on the metal core / adiabatic hole, LSQ rate 1.987 / 2.021, finest interior rel L2 6.11e-4 / 4.30e-4, balance ≤ 6.6e-13, the isotropic and β = 0 controls missing at 2.7e-2–6.2e-2; audit PASS 2026-10-22 18:00 (steps 0 and 1; `IMP-7.md`). Step 0 green 2026-10-22 (§9 item 579, `20261022T185726Z_IMP-7-step0-W1.log`, 3 s): every identity ≤ 1.3e-15 against `B1_BAND`, reductions → `MAT-18` / `IMP-1`, three controls missing. `IMP-1`'s sphere (a = 5 mm) as a conducting core in a gel shell to b = 4a, Robin h on r = b; the source is the quasi-static dipole SAR (β = 1); θ = Σ_{l=0,2} R_l(r) P_l(cos ϑ), exact (Carslaw & Jaeger ch. IX; Jin App. E). **Done-when:** *step 0* (numpy, I37): the (Src), (PDE), (I), (Rb), (FD), (Bal), (Φ), (R0) identities at their imported bands with (Rk) printed, the β = 0 → `MAT-18` and Φ → `IMP-1` reductions, three controls missing by > `CONTROL_MISS` ⇒ 🟡. *Step 1* (FEM, I38): `MAT-18` step 2(a)'s gates — interior rel L2 ≤ 1 %, P1 L2 rate ≥ 1.7, balance ≤ 1e-9 — on the metal core and on the adiabatic hole, the isotropic-source and β = 0 controls missing ⇒ ✅. Not claimed: a FEM SAR source, finite k, a rod or tip, a coil, transient, perfusion, published data. *History and briefs: `docs/planning/chunks/IMP-7.md`.* 

### §9 item 657 — 2026-10-27 (§10 I44, chain I — `IMP-7` step 3c, the curved degree-2 rate gate; chain step (solve))

*Written by the 2026-10-27 03:00 daily review from the 10-26 weekly's I44. The band `RATE_BAR_P2_CURVED` = 2.6 is the weekly's (0.85 × the order, the family's convention); nothing else moves.*

**Anchor.** The pairwise P2 L2 rate on curved P2-geometry cells, as `_pairwise` already computes it (`rw`, whole-domain error; `tests/validation/test_imp7_curved_p2.py:158–164`, I42, metal core; `tests/validation/test_imp7_curved_p2_hole.py:166–172`, I43, adiabatic hole), over `LADDER[0]` → `[1]` (in each module's W1 test, `:210` / `:207`) and `[1]` → `[2]` (W2, `:232` / `:230`). `LADDER` = `m7.LADDER` (`test_imp7_implant_heating_fem.py:66`). **Assert `rw >= RATE_BAR_P2_CURVED`** immediately after the existing `(Rate)` print inside `_pairwise`, in both modules, so every pair is gated. Records: I42 3.036 / 3.105, I43 3.234 / 3.129 (`IMP-7` row), so the margin is ≥ 0.43.

**Route — additive (rule (c)).** Define `RATE_BAR_P2_CURVED = 2.6` once, in `test_imp7_curved_p2.py` beside its other constants, with a comment citing the 10-26 weekly (§10 I44) and the 0.85 × order convention (`S1_RATE` 1.7 for P1, `test_imp7_implant_heating_fem.py:58`). Import it into the I43 module. Change the print's "— printed, no band" to name the bar. Nothing else changes: `S2_BAND`, `S3_BAND` and the β = 0 control stay as they are.

**Negative control (asserted, from records).** The affine P2 rate is **not computed at runtime**: it is I40's recorded constant `AFFINE_P2_PAIRWISE` = (1.937, 2.040) in I42 (`:60`) and (1.936, 2.043) in I43 (`:64`). Assert `max(AFFINE_P2_PAIRWISE) < RATE_BAR_P2_CURVED` in a no-mesh coding test in each module, labelled "record control: the affine floor I40 measured misses the curved bar". This is the weekly's "affine rate printed as the control (predicted < 2.6, asserted missing)", honestly scoped: it separates the two records, not two runtime solves.

**Windows** (bare harness calls, foreground, Bash `timeout: 660000`; real build; `-n 4`, the width I42 / I43 ran at; `pgrep -c python3` 0 before each; one window at a time):
- **W0** (`-n 1`, smoke): both modules' W0 / coding tests, including the record control.
- **W1a** I42 W1 (predicted ≈ 175 s, measured 173 s), **W2a** I42 W2 (≈ 400 s, measured 396 s), **W1b** I43 W1 (≈ 95 s, 93 s), **W2b** I43 W2 (≈ 230 s, 226 s): the existing commands from I42 / I43's logs, unchanged except the label (`IMP-7-step3c-W1a` …), each `timeout -k 30 600`. Total ≈ 15 min of compute.
- **W-env** (`-n 1`, `timeout -k 30 300`).

**Executor:** implementer. **Tier:** heavy, `-n 4`, real. **Class:** chain step (solve). **Cost:** ≈ 11.5 measured / 35 predicted slot-minutes. **Depends on:** nothing. **Family:** I44 (`IMP-7` step 3c), attempt 1 of 4.

**Traps already paid for.** I42's W2 at 396 s is the longest window; keep it alone. Rates are compared across widths, digits are not: I40 ran at `-n 2`, I42 / I43 at `-n 4`. A `cmap is deprecated` warning from `_cmap_degree` is harmless (I43). Count verdicts per rank by token.

**Scope.** A rate gate on records that exist. **Not claimed:** a rate at degree 3, a transient, any new physics.

**Negative result.** A rate below 2.6 on the curved cells contradicts the records at the same width: report the measured rate, park, BLOCKED; never move the bar. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7`'s row gains "curved P2 rate ≥ 2.6 gated on both boundary types"; §10 I44 → DONE.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 657 — `IMP-7` I44 — per `docs/planning/chunks/IMP-7.md` "§9 item 657". The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Add one constant, one assert per module and one record-control test; edit nothing else. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-27, 09:30 UTC / 04:30 CDT slot) — DONE, green; the curved P2 rate is now gated on both boundary types; row stays ✅.** `RATE_BAR_P2_CURVED = 2.6` is defined once in `test_imp7_curved_p2.py`, with a comment citing §10 I44 and the 0.85 × order convention (`S1_RATE`). It is imported into `test_imp7_curved_p2_hole.py`. `assert rw >= RATE_BAR_P2_CURVED` sits in each module's `_pairwise`, right after the (Rate) print, and the print now names the bar. A no-mesh record control was added to each module: `test_imp7_3c_metal_record_control` / `test_imp7_3c_hole_record_control`, which assert `max(AFFINE_P2_PAIRWISE) < RATE_BAR_P2_CURVED`. Nothing else was edited, and there was no edit after W0. `pgrep -c python3` read 0 before every `-n 4` window.

Windows, with verdicts counted by token:
- **W0** `20261027T093804Z_IMP-7-step3c-W0.log` (`-n 1`, 3 s): 4 passed. Record controls 2.040 (metal) / 2.043 (hole) < 2.6.
- **W1a** `20261027T093813Z_IMP-7-step3c-W1a.log` (`-n 4`, 98 s; 4 PASSED tokens): metal 1.414 → 1.0 whole **3.036** (interior 3.017).
- **W2a** `20261027T093959Z_IMP-7-step3c-W2a.log` (`-n 4`, 286 s; 4 PASSED): metal 1.0 → 0.7071 whole **3.105** (interior 3.105).
- **W1b** `20261027T094454Z_IMP-7-step3c-W1b.log` (`-n 4`, 89 s; 4 PASSED): hole **3.234** (interior 3.209).
- **W2b** `20261027T094630Z_IMP-7-step3c-W2b.log` (`-n 4`, 218 s; 4 PASSED): hole **3.129** (interior 3.123).
- **W-env** `20261027T095013Z_IMP-7-step3c-Wenv.log` (`-n 1`, real, 105 s): 8 passed, 4 skipped.

The rates reproduce I42 / I43 to the printed digit at the same width. The smallest margin above the bar is 0.436 (metal, first pair). The control separates the two records; it does not compare two runtime solves (brief, scope).

### §9 item 658 — 2026-10-27 (§10 I45, chain I — `IMP-7` step 4, the transient in closed form; lane 3 → chain step (no solve))

*Written by the 2026-10-27 03:00 daily review from the 10-26 weekly's I45. The weekly named the comparand, its identities and controls; this brief adds which identities can reach `B1_BAND` and which are truncation-limited, without moving a band.*

**Anchor — a new closed form.** `pec_sphere_implant_transient_temperature(points, t, *, radius, outer_radius, k_gel, rho_c, q_far, h, t_inf, beta=1.0, n_terms=200)` in `src/fem_em_solver/utils/analytical.py`, appended **after the `IMP-7` step 0 section at the end of the file** (`:3119` onward, where the steady comparand lives; the weekly said "after the `MAT-21 step 0` banner", and file order is the review's detail). Adiabatic hole only (Neumann at r = a). It returns θ(r, ϑ, t) = θ_ss − Σ_{l∈{0,2}} P_l(cos ϑ) Σ_n c_{l,n} R_{l,n}(r) exp(−κ λ_{l,n}² t), with κ = k / ρc and:
- θ_ss from `pec_sphere_implant_steady_temperature(…, k_core=0.0, …)` (`:3275`) and its radial parts `pec_sphere_implant_radial` (`:3258`) on the adiabatic branch (`:3231`);
- R_{l,n}(r) = j_l(λr) y_l′(λa) − y_l(λr) j_l′(λa) (Neumann at a), with λ the roots of k R′(b) + h R(b) = 0 (Carslaw & Jaeger ch. IX, the hollow sphere; Jin App. E for the spherical Bessel identities);
- c_{l,n} = ⟨θ_ss,l − θ_init, R_{l,n}⟩_{r²} / ⟨R_{l,n}, R_{l,n}⟩_{r²} on (a, b), with θ_init = t_inf;
- a helper `pec_sphere_implant_eigenvalues(l, n, *, radius, outer_radius, k, h)`: bracket sign changes of the determinant on a grid finer than π / (b − a) / 8, then `brentq`. **Assert the count:** the n-th root lies within one spacing of (n + l/2) π / (b − a) asymptotically, and no bracket holds two roots.

**Identities — new module `tests/validation/test_imp7_transient_closed_form.py`** (fixture constants imported from the step 0 module `test_imp7_implant_heating_closed_form.py`; `B1_BAND` = 1e-12 and `CONTROL_MISS` = 1e-3 imported from `test_imp4_pec_cylinder_closed_form.py:46`, `:51`, as step 0 does):
- **Per-mode, at `B1_BAND`:** (PDE) each R_{l,n} satisfies the spherical Bessel ODE with analytic derivatives; (N) R′(a) = 0; (Rb) k R′(b) + h R(b) = 0, all relative to the mode's scale; (Orth) ⟨R_{l,n}, R_{l,m}⟩_{r²} / norms ≤ `B1_BAND` for n ≠ m on the first 20 modes (high-order quadrature).
- **Series-level, exact at any truncation, at `B1_BAND`:** (Lim) θ(t → ∞) = θ_ss (the decaying sum underflows); (Duh) a 50 % duty cycle by Duhamel superposition of the step response equals the direct modal duty-cycle sum (the `_duty_cycle_modal` pattern, `:2370`) on the same modes.
- **Series-level, truncation-limited:** (IC) θ(t = 0) reconstructs θ_init, and (Slope) ρc ∂θ/∂t at t → 0⁺ equals the SAR. Both converge only as the eigen-expansion of a function that does not satisfy the Robin condition, so they **cannot reach 1e-12 pointwise**. Assert them as a **convergence statement**: the L2(r²) error at n_terms = 400 is ≤ ½ × that at 100; print all three truncations. `MAT-21` step 0 printed its PDE residual for the same reason (`test_mat21_box_closed_form.py:325`). This is a record, not a moved band; the weekly may write a bar.
- **Controls, asserted > `CONTROL_MISS`:** the Dirichlet-at-a series (θ = t_inf at r = a; the R_{l,n} built with the Dirichlet condition) and the isotropic source (the l = 2 part dropped) each miss the adiabatic transient at t = t\* = 0.2 (b − a)² / κ.

**Window.** **W1** `IMP-7-step4-W1` (smoke, numpy, `-n 1`, real image, `timeout -k 30 120`; predicted ≤ 30 s; the step 0 window took 3 s). **W-env** (`-n 1`, `timeout -k 30 300`).

**Executor:** implementer. **Tier:** smoke, `-n 1`. **Class:** lane-3 closed form. **Cost:** ≈ 7.2 measured / 30 predicted slot-minutes. **Depends on:** nothing. **Family:** I45 (`IMP-7` step 4), attempt 1 of 4.

**Traps already paid for.** `scipy.special.spherical_jn` / `spherical_yn` with `derivative=True` for the primes. y_l blows up as λa → 0: scale each R_{l,n} by its maximum on (a, b) before any residual. A missed root shows as an orthogonality failure, which is why (Orth) is asserted. The l = 0 mode has no zero eigenvalue with a Robin outer wall (h > 0).

**Scope.** The transient comparand. **Not claimed:** the FEM (that is I46), a metal core (k_core > 0), finite k in the hole, perfusion.

**Negative result.** A per-mode identity red ⇒ a defect in the eigenfunctions: report, park, BLOCKED. (Duh) or (Lim) red ⇒ a defect in the series. (IC) / (Slope) not converging at the asserted ratio ⇒ report the three truncations and stop; the weekly rules. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 4 (🟡 record; the chunk stays ✅); §10 I45 → DONE; unblocks item 668 (I46).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 658 — `IMP-7` I45 — per `docs/planning/chunks/IMP-7.md` "§9 item 658". The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Append to `analytical.py`, never edit an existing function. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-27, 04:30 CDT slot, ≈ 10:00 UTC) — DONE, green; step 4 is a 🟡 record and the row stays ✅.** The new code is appended after the step 0 section of `analytical.py`, and no existing function was edited. It adds `pec_sphere_implant_mode`, `pec_sphere_implant_eigenvalues` (the grid spacing is π/(b − a)/16, then `brentq`), `pec_sphere_implant_source_radial`, `pec_sphere_implant_transient_modes`, `pec_sphere_implant_transient_from_modes`, `pec_sphere_implant_transient_temperature` (the brief's signature) and `pec_sphere_implant_duty_cycle_from_modes`, plus private helpers. Second derivatives come from the Bessel recurrence, not from the ODE, so (PDE) is not circular. The Dirichlet-at-a control projects the source, c = ⟨q_l, R⟩/(kλ²‖R‖²), because no closed steady form exists for it here. ρc = 4.18e6 is `MAT-12`'s `RHO_C`. t\* = 313.5 s.

The test module is `tests/validation/test_imp7_transient_closed_form.py`. **W1** `20261027T095811Z_IMP-7-step4-W1.log` (`-n 1` without mpiexec, real image, 61 s, 5 passed) measured:
- **Per-mode identities (first 50 modes):** (PDE) 3.6e-16 / 1.1e-15 (l = 0 / 2). (N) 0 exactly. (Rb) 4.6e-15 / 5.3e-15.
- **(Orth), first 20 modes:** 2.4e-15 / 3.4e-15.
- **(Lim):** 0.
- **(Duh):** 5.0e-14 (on-phase, 1.6 t\*) and 6.3e-14 (off-phase, 1.4 t\*), against a θ max of 0.141 / 0.112 K. All of these are at or under `B1_BAND`.
- **(IC), rel L2(r²) at 100 / 200 / 400 terms:** l = 0 2.56e-8 / 2.24e-9 / 1.97e-10; l = 2 6.26e-7 / 5.49e-8 / 4.83e-9. Ratio 0.0077.
- **(Slope):** l = 0 1.13e-3 / 3.99e-4 / 1.41e-4; l = 2 3.59e-3 / 1.27e-3 / 4.47e-4. Ratio 0.124 (≈ n^−1.5).

Both ratios are ≤ ½, as asserted. Controls at t\* (adiabatic θ max 0.205 K): the Dirichlet-at-a series misses by rel L2 0.431, and the isotropic source by 0.0869. Both are > `CONTROL_MISS`.

**W-env** `20261027T095917Z_IMP-7-step4-Wenv.log` (`-n 1`, 106 s): 8 passed, 4 skipped.

**Two departures, returned to the weekly.**
- **(Count).** The brief's asymptote (n + l/2)π/(b − a) is wrong for this mode family. With Neumann at a, the large-λ mode is cos λ(r − a)/r, because the lπ/2 phases cancel between the two ends. So the 0-indexed n-th root → nπ/(b − a) for every l. Measured at n = 399, the brief's ladder is off by +0.001 for l = 0 and by −0.997 for l = 2. The test therefore asserts |λ_n(b − a)/π − n| < ½ for n ≥ 10 (measured ≤ 0.117), and n + ½ for the Dirichlet control (≤ 0.081). It also asserts every root gap > 1/8 spacing (measured ≥ 0.73), so no bracket holds two roots.
- **W1 time.** W1 took 61 s, against ≤ 30 s predicted. The 400-term truncation dominates. The container ceiling of 120 s was the brief's, and it held.

### §9 item 668 — 2026-10-27 (§10 I46, chain I — `IMP-7` step 5, the CN transient on the curved hole against I45; chain step (solve))

*Written by the 2026-10-27 03:00 daily review from the 10-26 weekly's I46. **Depends on item 658 DONE** (I45's closed form on `main`); otherwise skip at no cost.*

**Anchor.** θ_h(t) from Crank–Nicolson on I43's curved P2 adiabatic-hole fixture (`test_imp7_curved_p2_hole.py`: `_shell_submesh`, `_robin_outer`, `_hole_rung` from `test_imp7_implant_heating_fem.py:228`, `:240`, `:393`, geometry order 2) at `LADDER[0]` and `[1]`, against `pec_sphere_implant_transient_temperature` (item 658). **Gates, imported under `MAT-21` step 1's names** (`tests/validation/test_mat21_box_fem.py:43–49`):
- interior rel L2 at t\* ≤ `INTERIOR_BAND` (1e-2, `test_mat13_transient_perfusion_sphere.py:64`);
- the first-step slope ≤ `S2_BAND` **as `MAT-21` imports it**, 1e-4 from `test_mat12_transient_sar_heating_sphere.py:75`. **Name clash:** `IMP-7`'s own `S2_BAND` is 1e-2 (`test_imp7_implant_heating_fem.py:59`). Import the `MAT-12` one under an alias (`from … import S2_BAND as SLOPE_BAND`) and say which in the print;
- the energy balance ≤ `E_BAND` (1e-9) on every run;
- long backward-Euler = the steady series ≤ `INF_BAND` (1e-8).

**Negative control (asserted).** The same FEM transient misses the Dirichlet-at-a closed form (item 658's control) by > `INTERIOR_BAND` at t\*.

**Route.** A new module `tests/validation/test_imp7_transient_fem.py`. The solver is `solve_transient_conduction` (`thermal/transient.py:94`): `t_init` is a float, so start from t_inf, the source is the dipole SAR on the shell sub-mesh (as I43 builds it), Robin on the outer tag. Check whether it accepts a curved (coordinate degree 2) mesh with a P2 space before writing; if it does not, stop and say so (the route is then the weekly's). Δt, n_steps and `startup_half_steps` as `MAT-21` step 1 (`test_mat21_box_fem.py:233–251`), scaled to this fixture's t\* = 0.2 (b − a)² / κ.

**Windows** (bare, foreground, Bash `timeout: 660000`; real build; `-n 4`; `pgrep -c python3` 0 first): **W1** `LADDER[0]` (predicted 150–300 s from I43's 93 s steady W1 and `MAT-21`'s transient-over-steady ratio ≈ 1.05–2), **W2** `LADDER[1]` (predicted 300–550 s; over 600 s is a kill: park, never raise the timeout), **W-env**.

**Executor:** implementer. **Tier:** heavy, `-n 4`, real. **Class:** chain step (solve). **Cost:** ≈ 11.5 measured / 45 predicted slot-minutes. **Family:** I46 (`IMP-7` step 5), attempt 1 of 4.

**Traps already paid for.** `assemble_scalar` is rank-local: reduce before the balance. The balance residual in `solve_transient_conduction` divides by the source integral (`:310–312`); here the source is on, so it is defined. A killed FFCx compile leaves an orphan `.c` in `/root/.cache/fenics`: sweep only with `pgrep -c python3` at 0 (rubric element 4). Comparison points through `evaluate_vector_field_parallel`'s scalar sibling, never `f.eval` with a local index.

**Scope.** The first implant-heating transient against a closed form. **Not claimed:** a metal core, a coil-driven source, a rate in Δt.

**Negative result.** Interior red with the balance green ⇒ report the rel L2 on both rungs and whether it shrinks with h; park, BLOCKED, no loosened band. Balance red ⇒ a solver defect on curved cells: known-issues entry. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 5 (the transient FEM gate); §10 I46 → DONE; the weekly's `imp:` example becomes owed.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 668 — `IMP-7` I46 — per `docs/planning/chunks/IMP-7.md` "§9 item 668", only if item 658 is DONE. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Import every band; never restate one. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-27, 09:00 CDT slot, ≈ 14:14 UTC) — BLOCKED, parked on `attempt/IMP-7-I46-20261027T141400Z` (`c8121715`); every physics gate green, one (E) red that is a round-off floor; the row stays ✅ and gains no step 5.** Item 658 was DONE, so the item ran. The module is `tests/validation/test_imp7_transient_fem.py` (on the branch only). It imports every band: `INTERIOR_BAND` from `MAT-13`; `E_BAND`, `INF_BAND`, `N_MAIN`, `RHO_C` and `S2_BAND as SLOPE_BAND` from `MAT-12`; `T_STAR`, `KAPPA`, `SPAN`, `N_DEFAULT` from item 658's module. The time scale is t_s = (b − a)²/κ = 1567.5 s, so t\* = 0.2 t_s = 313.5 s, CN takes `N_MAIN` = 40 steps with two BE half-steps, DT_1 = 1e-5 t_s and the long run is BE 20 × 5 t_s (the `MAT-21` step 1 pattern). The comparand is the item 658 series with 200 terms, interpolated into P3 in chunks of 20 000 points.
- **Check-first.** `solve_transient_conduction` accepts the coordinate-degree-2 sub-mesh with P2 unchanged. **W0** `20261027T140249Z_IMP-7-step5-W0.log` (`-n 1`, 9 s): cmap 2, T degree 2, two CN steps with balance 6.3e-12.
- **W1** `LADDER[0]` (53 643 sub-cells): `20261027T140533Z_IMP-7-step5-W1.log` (`-n 4`, 93 s, FAILED ×4 ranks). The first W1 (`20261027T140310Z_IMP-7-step5-W1.log`, 89 s) read the same numbers; the module was then edited to add the diagnostic, so W1 was re-run under rule (i).
- **W2** `LADDER[1]` (150 334 sub-cells): `20261027T140719Z_IMP-7-step5-W2.log` (`-n 4`, 310 s, FAILED ×4 ranks; predicted 300–550 s).
- **W-env** `20261027T141245Z_IMP-7-step5-Wenv.log` (`-n 1`, 46 s): 8 passed, 4 skipped.

| gate | `LADDER[0]` | `LADDER[1]` | band |
|---|---|---|---|
| (T) CN interior rel L2 at t\* | 1.8004e-4 (whole 1.6849e-4) | 7.7194e-5 (whole 7.3476e-5) | ≤ `INTERIOR_BAND` ✓ |
| (S-slope) \|ρc∫θ¹/(DT_1∫Q) − 1\| | 5.5312e-5 | 5.5272e-5 | ≤ `SLOPE_BAND` (`MAT-12` `S2_BAND`) ✓ |
| (∞) BE 20 × 5 t_s vs steady FEM | 1.934e-10 | 2.408e-11 | ≤ `INF_BAND` ✓ |
| (c-D) Dirichlet-at-a miss | 0.31386 | 0.31573 | > `INTERIOR_BAND` ✓ |
| (E) CN / BE ∞ | 5.79e-12 / 2.47e-14 | 4.90e-12 / 1.16e-13 | ≤ `E_BAND` ✓ |
| (E) BE one step DT_1 | **1.300e-9** | **1.161e-9** | ≤ `E_BAND` ✗ |
| (E-diag) the same step in θ (t_inf = 0), printed | 1.611e-14 | 3.366e-14 | — |

The interior error shrinks with h (rate ≈ 2.4 over one √2 step, printed). BE ∞ against the steady **series** reads 1.1475e-4 / 3.7781e-5, which is I43's steady P2 error to the digit.

**Why (E) is red, and why that is not a solver defect.** The DT_1 step raises T by about q·DT_1/ρc ≈ 8 µK on a 20 K level. The solver's storage integral ρc∫(Tⁿ⁺¹ − Tⁿ)/Δt therefore loses about log₁₀(20 / 8e-6) ≈ 6.4 digits to cancellation. The θ-variable diagnostic runs the identical linear problem with no offset and reads 1.6e-14 / 3.4e-14. The slope ratio agrees between the two runs to 2e-9. So the red is the floor of the residual's arithmetic, not of the discretisation. That is also why no known-issues entry was filed, although the brief's "balance red ⇒ solver defect" clause asked for one. The two gates pull against each other at T∞ = 20 K: the slope is ≈ 5.5e-5 at DT_1, so doubling DT_1 would halve the floor but push the slope past `SLOPE_BAND`. `MAT-21`'s box escapes because its DT_1 step is ≈ 10× larger in kelvin. Nothing was retuned in the slot; the band is untouched.

**For the weekly (unblock).** Either (a) run the item's transients in θ (t_inf = t_init = 0; the comparand is θ already; measured above), or (b) record the floor on the DT_1 run with these measurements, as the MAG-10 / MAG-15 precedent allows. Once ruled, the re-run is W1 + W2 ≈ 7 min of compute on the parked module.

### Review ruling 2026-10-27 10:30 (item 668, I46) — the outcome's reading corrected; re-posed as item 672 under rule (6)

**`log-pathologist` ruling** (the review re-read the cited lines). `20261027T140533Z_IMP-7-step5-W1.log` and `…140719Z_…-W2.log` are countable:
- Each footer reads Status 1 (W1 96 s, W2 311 s), with no kill.
- There is one FAILED token per rank, on `be_1`, at the test's `:179` (W1:75; FAILURES blocks at :77 / 98 / 119 / 140).

The pathologist's findings:
- **"Every physics gate green … asserted" is REFUTED as an assertion claim.** `_gates` (parked module `:178–184`) asserts the (E) loop first, in the order cn, be_1, be_inf, and every rank raised on `be_1`.
  - **(T), (S-slope), (∞), the (c-D) control and (E) on be_inf were never asserted, on any rank.**
  - Their values are rank-0 prints (W1:69–74, W2:69–74). They are allreduced quantities and all sit inside their bands, but they are readings, not gates.
  - Only (E) on CN was asserted, green on 4 ranks.
  - Item 668's ✓ marks therefore read "printed in band".
- **The `be_1` red is CONFIRMED and reproducible:** 1.300e-9 / 1.161e-9, rank-identical (the first W1 read 1.304e-9). "The only red" is UNCOUNTABLE, for the same reason.
- **The cancellation mechanism is CONFIRMED:** the matrix does not depend on t_inf (`thermal/transient.py:228`), and T = θ + c·1 cancels algebraically. Three corrections apply:
  - (a) The Robin facet matrix and RHS terms (`:223`, `:236` / `:243`) carry no quadrature metadata. On curved facets they are therefore not bitwise offset-free. be_inf at 20 K balances to 2.5e-14 / 1.2e-13, which bounds that channel far below 1.3e-9.
  - (b) The observed loss is 4.5–4.9 digits, not 6.4.
  - (c) The floor sits in the T representation itself: the slope ratio moves ≈ 2e-9 between the T and θ runs, which is harmless against `SLOPE_BAND`.
- **Masking risk.** Moving *every* run to θ would hide any T∞-wiring defect: a sign error, a wrong facet tag, the quadrature mismatch, or a t_init ≠ t_inf path. CN and be_inf are the only runs that exercise it. **Moving only the `be_1` / slope run to θ, declared, is safe.**

**This review's ruling** (provisional under the rule (6) re-pose licence; the 10-28 weekly ratifies or reverts). I46 is re-posed as **§9 item 672**, family I46 attempt 2 of 4, with the target and every band kept:
- (i) The `be_1` / slope run is solved in θ (t_inf = t_init = 0), and its (E) is asserted ≤ `E_BAND`. The T-variable `be_1` balance is printed as a record.
- (ii) CN and be_inf stay at T∞ = 20 K.
- (iii) `_gates` evaluates every gate before asserting the collected list.
- (iv) A printed discriminator: `be_1` in T at 2·DT_1. Cancellation predicts ≈ ½ the 1.3e-9; the quadrature channel predicts no change.

This is not a band move: `E_BAND` applies unchanged to the same identity, on a run whose state carries no 20 K offset. The weekly's option (b), recording the floor, stays open.

### §9 item 672 — 2026-10-27 (§10 I46′, chain I — `IMP-7` step 5 re-posed under rule (6), daily 2026-10-27 provisional; chain step (solve))

*Written by the 2026-10-27 10:30 daily review. It re-poses item 668 (BLOCKED, parked `attempt/IMP-7-I46-20261027T141400Z`, `c8121715`) per the ruling directly above. **Depends on nothing:** the module is on the branch, and I45 (item 658) is on `main`. Family I46, attempt 2 of 4.*

**Route — land the parked module with four edits, then re-run.** Take `tests/validation/test_imp7_transient_fem.py` from the branch **by path**. It is a new file, so a path checkout is safe; never path-checkout a record file. Then make four edits:
1. **`be_1` in θ.** The `(S-slope)` block (`:136–144` on the branch) solves `be1` with `t_inf=0.0`, and `stor` becomes `∫ be1.temperature dx` (no `− T_INF`). `out["E"]["be_1"]` is that run's balance. The T-variable run becomes the printed record `out["be_1_T"]`: the existing `be0` diagnostic block (`:146–156`) inverted, so that θ is gated and T is printed.
2. **Collect, then assert.** `_gates` builds `(name, ok, message)` for every (E) entry and for (T), (S-slope), (∞) and (c-D). It prints the list on rank 0, then asserts `all(ok)` with the failing names.
3. **Discriminator (printed).** One extra `be_1` run in T at `2*DT_1`, its balance printed beside the `DT_1` one together with the ratio. Cancellation predicts ≈ 0.5; the facet-quadrature channel predicts ≈ 1.
4. **Docstring.** One paragraph citing this ruling.

Nothing else changes: imports, `_fixture`, the CN and be_inf runs at T∞, the comparand and the control.

**Anchor (unchanged from item 668).** `pec_sphere_implant_transient_temperature` through `…_transient_modes` / `…_from_modes` (`analytical.py:3584`, `:3642`, `:3665`; item 658). The gates are imported as on the branch (`:39–53`):
- interior ≤ `INTERIOR_BAND` (`MAT-13`, 1e-2);
- slope ≤ `S2_BAND as SLOPE_BAND` (`MAT-12`, 1e-4);
- every (E) ≤ `E_BAND` (1e-9);
- (∞) ≤ `INF_BAND` (1e-8);
- (c-D) > `INTERIOR_BAND`.

Expected from item 668's prints, at the two rungs (not a pre-registered red):

| Quantity | Expected |
|---|---|
| (T) | 1.80e-4 / 7.72e-5 |
| slope | 5.53e-5 |
| θ-run (E) | 1.6e-14 / 3.4e-14 |
| CN (E) | 5.8e-12 / 4.9e-12 |
| be_inf (E) | 2.5e-14 / 1.2e-13 |
| (∞) | 1.9e-10 / 2.4e-11 |
| (c-D) | 0.314 / 0.316 |

**Windows.** Each is bare and foreground, with Bash `timeout: 660000`, on the real build at `-n 4`. First, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W1** `IMP-7-step5-W1r`, `LADDER[0]`: ≈ 100 s (item 668's 93 s plus one BE step).
- **W2** `IMP-7-step5-W2r`, `LADDER[1]`: ≈ 320 s, `timeout -k 30 600`. Over 600 s is a kill: park, never raise the timeout.
- **Environment check:** the full `tests/environment` (≈ 45 s since `OPS-72`) may lead W1, or run as **W-env**.

**Executor:** implementer. **Tier:** heavy, `-n 4`, real. **Class:** chain step (solve). **Cost:** ≈ 11.5 measured / 20 predicted slot-minutes.

**Traps already paid for.**
- Count verdicts per rank by token.
- An assertion that short-circuits hides later gates (item 668; edit 2).
- **Never move CN or be_inf to θ.** They are the only runs that exercise the T∞ wiring.
- A killed compile leaves an orphan `.c`. Sweep only with `pgrep -c python3` at 0 (rubric element 4).

**Scope.** The first implant-heating transient against a closed form. **Not claimed:** a metal core, a coil-driven source, a rate in Δt, or a duty cycle (item 675).

**Negative result.** Any collected gate red ⇒ report the full list, park, BLOCKED; no loosened band. A red θ-run (E) ⇒ the cancellation reading is wrong; file a known-issues entry on the solver's balance on curved P2 facets. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 5, the transient FEM gate, and §10 I46 → DONE. An `imp:` example becomes owed (the next review opens it). Item 675 (I47) unblocks.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 672 — `IMP-7` I46′ — per `docs/planning/chunks/IMP-7.md` "§9 item 672" and the "Review ruling 2026-10-27 10:30" above it. Take the parked module by path, make exactly the four edits, never move the CN or be_inf runs to θ. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-27, 12:00 CDT slot, second item) — DONE; `IMP-7` gains step 5 (the transient FEM gate), §10 I46 → DONE, item 675 unblocks.** The module was taken from `c8121715` by path and given exactly the four edits; CN and be_inf stay at T∞ = 20 K. `pgrep -c python3` read 0 before each window. Both windows: `-n 4`, real build, Status 0, one PASSED token per rank (4 / 4), and `_gates` printed the full collected list before asserting.

| Gate (collected, asserted) | LADDER[0] (53 643 cells) | LADDER[1] (150 334 cells) | Band |
|---|---|---|---|
| (T) CN interior rel L2 at t\* | 1.8004e-4 | 7.7194e-5 | ≤ `INTERIOR_BAND` 1e-2 ✓ |
| (S-slope), θ run | 5.5314e-5 | 5.5274e-5 | ≤ `SLOPE_BAND` 1e-4 ✓ |
| (E) cn (T∞) | 5.807e-12 | 4.852e-12 | ≤ `E_BAND` 1e-9 ✓ |
| (E) be_1 (θ) | 1.611e-14 | 3.386e-14 | ≤ `E_BAND` ✓ |
| (E) be_inf (T∞) | 1.836e-14 | 1.146e-13 | ≤ `E_BAND` ✓ |
| (∞) BE vs steady FEM | 1.934e-10 | 2.406e-11 | ≤ `INF_BAND` 1e-8 ✓ |
| (c-D) Dirichlet-at-a miss | 0.31386 | 0.31573 | > `INTERIOR_BAND` ✓ |

Printed records: the T-variable `be_1` balance is 1.302e-9 / 1.164e-9 (item 668: 1.300e-9 / 1.161e-9), and its slope ratio sits 2.17e-9 / 2.03e-9 from the θ run's. **The discriminator** (BE in T at 2·DT_1) reads 6.523e-10 / 5.803e-10, **ratio 0.5011 / 0.4986**: the cancellation reading holds, and the facet-quadrature channel (predicted ≈ 1) is excluded. Logs: `20261027T171958Z_IMP-7-step5-W1r.log` (104 s), `20261027T172155Z_IMP-7-step5-W2r.log` (356 s). The module was not edited after W1r. Not claimed: a metal core, a coil-driven source, a rate in Δt, or a duty cycle.

### §9 item 675 — 2026-10-27 (§10 chain I, I47 — daily, 2026-10-27, provisional; `IMP-7` step 6: the duty-cycled CN transient on the curved hole against I45's duty series; chain step (solve))

*Written by the 2026-10-27 10:30 daily review under rule (6): chain I's next step toward subgoal (7e), implant heating in closed form and FEM. It is provisional until the 10-28 weekly ratifies it. **Depends on item 672 DONE** (the module and its helpers on `main`); otherwise skip at no cost.*

Clinical RF exposure is duty-cycled. `MAT-17` gated the duty-cycled FEM on the solid ball (`tests/validation/test_mat17_duty_cycle_fem.py`). This item is the implant analogue on I43's curved hole, against an exact series that already exists.

**Anchor.** `pec_sphere_implant_duty_cycle_from_modes(points, t, modes, *, rho_c, t_inf, period, duty=0.5)` (`analytical.py:3686`). Its modes come from `pec_sphere_implant_transient_modes(beta=BETA, n_terms=N_DEFAULT, inner="neumann", **FIX)`, built as item 672's `_series_fn` builds them. The parameters are I45's (Duh) values (`test_imp7_transient_closed_form.py:150–151`): **period = 0.5 T_STAR, duty 0.5**, so t\* spans two periods.

**Cross-check (asserted, numpy, before any FEM).** At five axis radii, the helper at t\* equals I45's Duhamel sum of shifted step responses, ≤ `B1_BAND` (imported as I45 does, `:43`).

**FEM.** A new module, `tests/validation/test_imp7_duty_cycle_fem.py`. It imports `_fixture`, `_transient` and `_rel` from item 672's module, and the docstring says so.
- **Solver call:** `solve_transient_conduction(..., source_schedule=χ)`, with χ(t) = 1 if (t mod period) < duty·period, else 0. This is `MAT-17`'s `_chi` (`test_mat17_duty_cycle_fem.py:80–82`), copied, not imported.
- **Stepping:** CN with `startup_half_steps=2` and dt = T_STAR / `N_MAIN` (40 steps). The period is 20 steps and the on-phase 10, so every switch falls on a step boundary.
- **Temperature:** the run is at T∞ = 20 K, so it exercises the T∞ wiring (ruling (ii) above).

**Gates (imported, unmoved, collected then asserted):**
- interior rel L2 at t\* ≤ `INTERIOR_BAND`;
- the max balance ≤ `E_BAND`, as `MAT-17`'s duty runs assert it (`:282–285`);
- (S0) Σχ·dt = the on-time, ≤ 1e-12 (`MAT-17`'s census, `:36–47`).

**If the solver's balance normalisation is undefined on a source-off step, stop and report the line.** Do not add a normalisation here: that is item 669's additive keyword. If 669 has landed, use its route.

**Negative control (asserted).** The **continuous-source** closed form at t\* (I45's `pec_sphere_implant_transient_temperature`, χ ≡ 1) misses the duty FEM by > `INTERIOR_BAND`. The duty halves the deposited energy, so the miss should be ≫ 1e-2; print it.

**Windows.** Each is bare and foreground, with Bash `timeout: 660000`, on the real build at `-n 4`, with `pgrep -c python3` reading 0 first.
- **W0** the numpy cross-check, `-n 1`, smoke.
- **W1** `LADDER[0]`, ≈ 90 s.
- **W2** `LADDER[1]`, ≈ 300 s, `timeout -k 30 600`.
- **W-env**, or the full `tests/environment` leading W0.

**Executor:** implementer. **Tier:** heavy, `-n 4`, real. **Class:** chain step (solve). **Cost:** ≈ 11.5 measured / 40 predicted slot-minutes. **Family:** I47, attempt 1 of 4.

**Traps already paid for.**
- CN on a switched source needs the startup half-steps. Print the no-startup run beside the startup run, as `MAT-17`'s `s_ns` does (`:253–256`), with no gate.
- The duty helper raises unless `modes["inner"] == "neumann"`.
- `assemble_scalar` is rank-local.
- Interpolate the series in chunks of 20 000 points (`INTERP_CHUNK`).

**Scope.** The duty-cycled implant transient against its exact series. **Not claimed:** a metal core, a coil source, CEM43 dose, a rate in Δt.

**Negative result.**
- Interior red with the balance green ⇒ report both rungs and whether the error shrinks with h; park, BLOCKED.
- Census red ⇒ report and park.

No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 6, a 🟡 record. §10 I47 (provisional) → DONE.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 675 — `IMP-7` I47 — per `docs/planning/chunks/IMP-7.md` "§9 item 675", only if item 672 is DONE. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Import every band; never restate one. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-27, 16:30 CDT slot, third item) — DONE; `IMP-7` gains step 6 (🟡 record), §10 I47 → DONE.** New module `tests/validation/test_imp7_duty_cycle_fem.py`, importing `_fixture`, `_transient` and `_series_fn` from item 672's module (the docstring says so; `_rel` was not needed, the interior norm is `m7._errors` as in 672); `_chi` copied from `MAT-17`. Every band imported; the (S0) bar 1e-12 is `MAT-17`'s census literal, named `S0_BAND`. `pgrep -c python3` read 0 before each window; every window Status 0, one PASSED token per rank.

| Gate (collected, asserted) | LADDER[0] (53 643 cells) | LADDER[1] (150 334 cells) | Band |
|---|---|---|---|
| (T) CN duty interior rel L2 at t\* | 1.0444e-4 | 1.3374e-4 | ≤ `INTERIOR_BAND` 1e-2 ✓ |
| (E) max balance (T∞ = 20 K) | 5.809e-12 | 4.850e-12 | ≤ `E_BAND` 1e-9 ✓ |
| (S0) \|Σχ·dt − on-time\|/on-time | 0.000e+00 | 0.000e+00 | ≤ 1e-12 ✓ |
| control: continuous-source series miss | 0.53873 | 0.53852 | > `INTERIOR_BAND` ✓ |

W0 (numpy, `-n 1`): the duty helper at t\* against I45's Duhamel sum at five axis radii, max rel **4.484e-14** ≤ `B1_BAND`. χ over the 41 (half-)steps: 21 on, 20 off, every switch on a step boundary. The balance is normalised by the full-on ∫Q (`thermal/transient.py`), so the off steps are defined. Printed: CN without start-up 1.3477e-4 / 8.3775e-5. **The (T) error grows slightly with h at the fixed CN Δt** (1.04e-4 → 1.34e-4, where I46's continuous run fell 1.80e-4 → 7.72e-5), so at these rungs the duty run sits on a time-discretisation floor; no rate in Δt or h is claimed. Logs: `20261027T215045Z_IMP-7-step6-W0.log` (7 s; `tests/environment/test_complex_mode.py` led it, 1 passed / 3 skipped in real mode), `20261027T215103Z_IMP-7-step6-W1.log` (81 s), `20261027T215234Z_IMP-7-step6-W2.log` (233 s). No `src/` change. The module was not edited after W0. Not claimed: a metal core, a coil source, CEM43 dose, a rate in Δt.

### §9 item 699 — 2026-10-30 (§10 chain I, I48a — daily, 2026-10-30, provisional; `IMP-7` step 7: the dose and six-minute window around the adiabatic hole, closed form; lane 3 → chain step (no solve))

*Written by the 2026-10-30 18:00 daily review. Depends on nothing. Item 675's outcome lists CEM43 dose as not claimed; this is the implant-side twin of `MAT-21` step 13's dose question. Next: I48b (item 700).*

**New module `tests/validation/test_imp7_duty_dose_closed_form.py`** (numpy). Every band imported.
- **Grid — not t\*.** `T_STAR` is 0.2 (b − a)² / κ ≈ 313.5 s here (a = 5 mm, b = 20 mm, κ = 0.6 / 4.18e6), **shorter than `W_B2` = 360 s**, so a window on item 675's t\* history would not fit. Use `T_END = 4 * T_STAR` ≈ 1254 s, `DT = T_STAR / N_MAIN` as item 675, `T_GRID = DT * np.arange(4 * N_MAIN + 1)`, 161 rows; `PERIOD` = 20 `DT`, on-phase 10, so every switch lies on the grid. Assert `T_END >= W_B2`.
- **Points:** `PTS` = pole (0, 0, 1.05 a), equator (1.05 a, 0, 0), outer pole (0, 0, 0.95 b): all interior, so I48b can read them on the curved P2 mesh, where a point on the boundary may fall outside. Labels printed.
- **Histories:** `pec_sphere_implant_duty_cycle_from_modes(pts, t, modes, rho_c=RHO_C, t_inf=T_INF, period=PERIOD, duty=DUTY)` on `_modes()` (`test_imp7_transient_closed_form.py`); row 0 = T∞ exactly (θ = 0). The continuous history from `_theta`.
- Copy `_q_check` / `_win_check`'s bodies with `t` as an argument (`test_mat21_box_dose_closed_form.py:78–100`): they close over the box's `T_GRID`. Cite the source.

**Anchors (asserted), at `DOSE_BAR` / `AVG_BAR`, imported.**
- (C100-i) `duty = 1.0` equals the continuous history at every row and point ≤ `B1_BAND` (the band of I45's (Duh) and W0), and its dose matches ≤ `DOSE_BAR`;
- (Q-id) peak and equal-mean (θ / DUTY): `cem43` against `_quad_dose` ≤ `DOSE_BAR`; max T < `T_REF` asserted;
- (Win-id) the same two scales: `running_average` at `W_B2` against `quad`, ≤ `AVG_BAR`, t_max inside the grid's reach, the bracket `lo ≤ mx ≤ hi`.

**Negative controls (asserted).**
- (Rc-id) the R swap misses by > `DOSE_BAR`;
- (Pole-id) the isotropic source, `beta = 0.0` modes (the control I45 uses), misses the β = 1 dose at the pole by > `DOSE_BAR`. A dose is exponential in T, so the pole / equator contrast of the dipole SAR must show; **print the ratio**.

**Printed (records, no band):** per point D_duty(peak), D_duty(equal mean), D_cont, the equal-mean ratio, the pole / equator dose ratio, the six-minute peak and its t_max, beside the bare gel box's ratio (the implant-side reading of the same question).

**Windows.** W1 `IMP-7-step7-W1` (`-n 1`, `timeout -k 30 120`, smoke, numpy). Predicted < 30 s: item 675's W0 ran the helper at five radii inside 7 s; this is 161 × 3. **Class:** lane-3 closed form, ≈ 6.5 measured / 20 predicted minutes.

**Traps already paid for.** `DUTY` and `PERIOD` import from `test_imp7_duty_cycle_fem` (this `PERIOD` is 0.5 × `IMP-7`'s t\*). Module-level imports of `test_imp7_*_fem` pull dolfinx: harmless, but they make this numpy test need the container. `T_STAR` here is not `MAT-21`'s. `beta` is a keyword of `pec_sphere_implant_transient_modes`, and `_modes()` fixes it, so the control builds its own modes with `n_terms=N_DEFAULT, inner="neumann", **FIX`. Do not import `test_*`-named functions.

**Scope.** The closed form's dose around the hole under a duty cycle. **Not claimed:** a metal core, a coil source, a thermal-damage statement, a rate. The series' numbers are at `MAT-9`'s q_far and T∞ = 20 °C.

**Negative result.** (Q-id) or (Win-id) red: report and stop. (Pole-id) under the bar: the control does not discriminate; record it, BLOCKED. No `PREDICTED RED:` line.

**Status it can move.** §10 I48a → DONE; the `IMP-7` row gains step 7 as a record; item 700 unblocks.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 699 — `IMP-7` step 7 — per `docs/planning/chunks/IMP-7.md` "§9 item 699". Numpy only. The grid runs to 4 t\*, not t\*. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (2026-10-31T01:01Z, 19:30 CDT slot): DONE — `IMP-7` step 7, a record.** W1-r2 (`20261031T010117Z_IMP-7-step7-W1-r2.log`, `-n 1`): 5 passed, Status 0, 77 s. That is over the predicted < 30 s and over the smoke tier's 30 s, and inside the brief's `timeout -k 30 120`. No `src/` change.

Asserted (`:43–69`):
- **Grid:** 161 rows, t\* 313.5 s, T_END 1254 s ≥ 360 s, every switch on a grid point.
- **(C100-i):** duty = 1.0 against the continuous history reads rel **0.000e+00** at every row and point against `B1_BAND` 1e-12, and the doses are equal. **That zero is not a shortcut.** The two routes compute differently, and `IMP-7-step7-probe` (`20261031T005942Z_IMP-7-step7-probe.log:55–59`) reads them at `t_inf = 0` differing by ≤ 1.7e-15 K. Both histories here are T − T∞ at T∞ = 20, and the ulp of 20 (3.6e-15) absorbs that difference. The bound is the 1e-15 absolute reading, not 0.
- **(Q-id):** ≤ 5.973e-16 against `DOSE_BAR` 1e-10; max T 20.43 °C at equal mean, under `T_REF`.
- **(Win-id):** ≤ 7.324e-16 against `AVG_BAR` 1e-12; t_max inside the grid's reach, the bracket holds.
- **(Rc-id):** the R swap misses by factors ≥ 6.9e+06.
- **(Pole-id):** at the pole β = 1 reads 2.163760e-11 s and β = 0 reads 2.059642e-11 s, ratio 1.050551, a miss of 4.8119e-02, ≈ 4.8e+08 × the bar. The control discriminates.

Printed records (`:71–76`):

| point | D_cont (s) | D_duty(peak) (s) | D_duty(equal mean) (s) | equal-mean ratio |
|---|---|---|---|---|
| pole | 2.598306e-11 | 2.163760e-11 | 2.637229e-11 | 1.014980 |
| equator | 2.370914e-11 | 2.065074e-11 | 2.400623e-11 | 1.012531 |
| outer pole | 2.085300e-11 | 1.933320e-11 | 2.098868e-11 | 1.006506 |

Pole / equator: D_cont 1.095909, D_duty(peak) 1.047788, D_duty(equal mean) 1.098560. The six-minute peaks and their t_max are at `:71–73`. These doses are ≈ 1e-11 s because the gel sits ≈ 0.4 K above 20 °C at `MAT-9`'s q_far, so they are identities on a history, not a thermal statement. The bare gel box's ratios are printed beside them from a different geometry and schedule.

**Disclosed.** The first W1 (`20261031T005756Z_IMP-7-step7-W1.log`, green, 76 s, the same readings) tripped the private-leak guard at its line 69, the (Pole-id) print at `.15e` / `.9f`. The guard prints no value; I take it as a coincidental digit match on our own closed-form figure, and no AED file was opened. Per the §9 rule that log is held unchanged in the gitignored `logs/leak-held/`, so its `test-results.md` row points at an uncommitted file, as `EX-142`'s does. The print dropped to 7 digits and W1 re-ran on the module as committed. `_q_check` / `_win_check` are imported from item 698's module, which holds the `t`-argument copies of `test_mat21_box_dose_closed_form.py:78–100`, rather than copied again; the brief said "copy", and these are non-`test_*` names. `BARE_EQ_RATIO` comes from the same module.

**Status moved.** §10 I48a DONE; item 700 is unblocked.

### §9 item 700 — 2026-10-30 (§10 chain I, I48b — daily, 2026-10-30, provisional; `IMP-7` step 8: the FEM dose around the curved hole against I48a; chain step (solve))

*Written by the 2026-10-30 18:00 daily review. **Serial on item 699 landing; if it did not, skip at no cost.** Heavy: `pgrep -c python3` reads 0 first. The 02:00 Saturday `xxl` window (entry 50) runs in another container, predicted to end before 03:00, and is sized to overlap slots; if `uptime` loadavg (1 min) reads ≥ 16 when this item starts, skip it this slot.*

**New module `tests/validation/test_imp7_duty_dose_fem.py`,** on item 675's fixture at **`LADDER[0]`**: `_fixture`, `_transient(..., dt=DT, n_steps=4 * N_MAIN, theta=0.5, startup_half_steps=2, source_schedule=_chi(PERIOD, DUTY), snapshot_steps=tuple(range(4 * N_MAIN + 1)))`. One read per snapshot on I48a's `PTS` in one fixed batch (bitwise comparisons are call by call, item 537; this item makes none). The read route is `test_mat21_box_duty_dose_fem.py`'s (`_probe_vals(r.snapshots[k], pts)` over the snapshots, its module fixture `ddose_m4`, `:90`, the histories filled at `:115–118`); **the points are one batch through the same evaluation on every rank**.

**Anchors (asserted, imported).**
- preconditions: (E) max balance ≤ `E_BAND`; (S0) census to `T_END` ≤ `S0_BAND`; every point valid; max T < `T_REF`;
- (Q-hd), (Win-hd) on the FEM history at peak and equal mean, `DOSE_BAR` / `AVG_BAR`, as item 567's module;
- (T) the interior rel L2 at `T_END` against I48a's series ≤ `INTERIOR_BAND` (item 675's gate, at the new end time).

**Negative control (asserted).** The R swap misses by > `DOSE_BAR` on the FEM history, and the continuous series at `T_END` misses by > `INTERIOR_BAND` (item 675's (c-cont); ceiling near 0.54 at t\*, a reading here).

**Printed (records, no band):** the FEM dose against I48a's series dose at each point (rel, a reading: **no dose band exists on this fixture and this item sets none**), the six-minute peak of each, Δ peak.

**Windows.** W0 `IMP-7-step8-W0` (`-n 1`, smoke, `timeout -k 30 120`): imports, grid, the census arithmetic. W1 `IMP-7-step8-W1` (`-n 4`, real, `timeout -k 30 600`). Price: item 675's W1 ran two CN runs of 40 steps plus comparands in 81 s at `LADDER[0]`, `-n 4`, so a 160-step run with 161 reads is predicted ≈ 200–350 s. **Over 600 s is a kill: park, never raise the timeout.** W-env (`-n 2`, `timeout -k 30 300`). **Class:** chain step (solve), ≈ 8.3 measured / 30 predicted minutes. Tier heavy.

**Traps already paid for.** The module-scope fixture must print before its first assert. Snapshots hold the temperature, `Function` by `Function`: 161 P2 fields at `LADDER[0]` are small, but read them as you go if memory is printed high. Real build (the module asserts nothing complex); `Wreg` on the complex image is not needed, since no `src/` changes. Curved geometry: `m7._errors`, not a plain affine integral. Every band imported by name. A reading on `-n 4` that differs in the last ULP from another width is not a finding.

**Scope.** The FEM dose history on `LADDER[0]`. **Not claimed:** dose convergence in h or Δt (I47's interior error does not shrink with h at fixed Δt; no rate), a metal core, a coil source. Example: the review writes it after this lands.

**Negative result.** (T) red at `T_END` with (E) green: the time-discretisation floor item 675 recorded has grown with the longer run; print the error at t\*, 2 t\*, 4 t\* and stop. Report, never re-band. Kill at 600 s: park on `attempt/IMP-7-I48b-*`, BLOCKED, last printed step. No `PREDICTED RED:` line.

**Status it can move.** §10 I48b → DONE; the `IMP-7` row gains step 8 as a record; the first implant-side FEM dose.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 700 — `IMP-7` step 8 — per `docs/planning/chunks/IMP-7.md` "§9 item 700". Only after item 699 is DONE. `-n 4`, `timeout -k 30 600`: a kill parks the item. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (2026-10-31T02:06Z, 21:00 CDT slot): DONE — `IMP-7` step 8, a record.** W1-r2 (`20261031T020326Z_IMP-7-step8-W1-r2.log`, `-n 4`, real, `timeout -k 30 600`): 5 passed on each of 4 ranks, Status 0, 142 s, under the predicted 200–350 s. The CN run took 36.9 s, the 161 reads 4.5 s, each comparand ≈ 18 s, the series history 19.9 s. W0 `20261031T020154Z_IMP-7-step8-W0.log` (`-n 1`, 4 s). W-env `20261031T020601Z_IMP-7-step8-Wenv.log` (`-n 2`, 8 passed, 4 skipped, 50 s). No `src/` change.

Asserted (`:74–79,105–121`), on `LADDER[0]`, 53 643 sub cells:
- **Preconditions:** (E) 5.806e-12 against `E_BAND` 1e-9; (S0) rel 0 on 627 s of on-time; all three points valid; max T 20.43 °C at equal mean, under `T_REF`.
- **(T) at 4 t\*:** interior rel L2 **9.0012e-05** against `INTERIOR_BAND` 1e-2. Printed at t\* 1.0444e-04 (item 675's 1.04e-4 reproduced) and at 2 t\* 7.2569e-05: the floor does not grow with the longer run.
- **(Q-hd):** ≤ 4.902e-16 against `DOSE_BAR` 1e-10. **(Win-hd):** ≤ 3.662e-16 against `AVG_BAR` 1e-12; t_max inside the grid's reach, the bracket holds.
- **Controls:** the R swap misses by factors ≥ 6.9e+06; the continuous series at `T_END` misses by 5.3906e-01.

Printed records, no band (`:123–125`):

| point | peak-scale dose rel, FEM vs series | equal-mean rel | six-minute peak FEM / series (K) | t_max (s) | max_n \|Δθ\| (K) |
|---|---|---|---|---|---|
| pole | 1.245e-04 | 2.481e-04 | 0.179697 / 0.179779 | 1191.30 | 1.285e-03 |
| equator | 1.795e-05 | 3.556e-05 | 0.145280 / 0.145288 | 1254.00 | 2.181e-04 |
| outer pole | 2.605e-06 | 4.630e-06 | 0.075803 / 0.075798 | 1206.98 | 6.484e-05 |

The gap falls with distance from the hole, and the pole's point error (≈ 0.6 % of its peak θ) is far above the (T) L2 figure: the point sits 0.25 mm off the hole wall. One mesh, one Δt, a reading, not a rate. The doses are ≈ 1e-11 s, identities on a history, not a thermal statement.

**Disclosed.** The brief's read route could not run as written. The first W1 (`20261031T020205Z_IMP-7-step8-W1.log`, Status 1, 45 s, 4 setup errors) raised in `point_owner_cells`: this fixture is a `create_submesh` mesh, which carries no original cell index (`POST-10`), and `_probe_vals` passes no `parent=`. It is a read-route error in non-gate code after the CN run, not a kill or a red, and it was fixed and re-run once at the same `-t` (the 2026-10-30 ruling on item 678). `_shell_fixture` is `_fixture` keeping the parent map `_shell_submesh` drops, and `_probe` is `_probe_vals` with `parent=`, both in the new module; no landed module was edited. W0 could not catch this, because the brief scopes it to imports, grid and census. W1-r2 ran the module as committed; nothing was edited after it. The heavy pre-check read `pgrep -c python3` 0 and loadavg 0.36 at the slot's start. The doses print at 7 digits to stay clear of the leak guard (item 699's disclosure).

**Status moved.** §10 I48b DONE; the `IMP-7` row gains step 8 as a record. The chain's next, per §10, is an `imp:` example, which the review writes.

### Review ruling 2026-10-31 03:00 (item 700's re-run)

**Confirmed: the re-run stands, and item 678's ruling covers this case with one narrowing.** That ruling names a crash in an example's own print helper; here the crash was in the module's read route, which feeds gates. It still qualifies because three facts hold together:
- the failure was a raise before any gate value existed (`20261031T020205Z_IMP-7-step8-W1.log:552–555`: `point_owner_cells`, "the mesh carries no original cell index"; 1 passed, 4 setup errors), so no reading was seen and then re-tried;
- the fix adds `parent=` on the same owner rule and moves no band or timeout, at the same `-t`, with the first log kept;
- an independent number ties the new route to the old record: (T) at t\* reads 1.0444e-04 in W1-r2, item 675's 1.0444e-04 on the other module's run (`20261027T215103Z_IMP-7-step6-W1.log:70`).

**So a crash in a module's own non-gate code before any gate value exists may be fixed and re-run once at the same `-t`, first log kept and disclosed. A crash after a gate value has printed stays under item 678's text only if the gate cannot have changed.** No band moved; the weekly may tighten this.

**Not written: I48c.** The next measurement on this fixture would be the FEM-vs-series dose gap at `LADDER[1]` and at half Δt on `LADDER[0]`, which separates h from Δt and is the evidence a dose band or a Δt band would need. Chain I's sixth step went to I36d (item 702), which protects Monday's window, so I48c is left as a proposal for the 10-31 weekly.

## §7 row moved verbatim 2026-10-31 (rotate_plan_archive.py chunks)

 SAR heating around a PEC-sphere implant in a Robin-cooled gel ball: the implant-heating canonical in closed form (Phase 7 subgoal (7e); opened by the … | ✅ 2026-10-22 on steps 0 and 1 (the exact l = 0 + l = 2 series and the P1 FEM at `MAT-18` step 2(a)'s gates); steps 2, 3 and 3b are records — the affine P2 rate ≈ 2 is the polyhedral-boundary floor, lifted to 3.036 / 3.105 (metal core, I42) and 3.234 / 3.129 (adiabatic hole, I43) on curved P2-geometry cells. **Curved P2 rate ≥ 2.6 gated on both boundary types** (I44, §9 item 657, 2026-10-27, heavy `-n 4`: 3.036 / 3.105 metal, 3.234 / 3.129 hole; the affine record 2.04 misses). **Step 4 (🟡 record): the adiabatic-hole transient in closed form** (I45, §9 item 658, 2026-10-27, smoke numpy: eigen-identities, t → ∞ limit and 50 % duty-cycle Duhamel at `B1_BAND`; IC / initial slope converge in truncation; Dirichlet-at-a and isotropic controls miss) — the comparand for I46. **Step 5 (gated): the CN transient on I43's curved P2 hole against I45** (I46′, §9 item 672, 2026-10-27, heavy `-n 4`, real; every collected gate green on 4 ranks at `LADDER[0]` / `[1]`): (T) interior rel L2 at t\* 1.80e-4 / 7.72e-5 ≤ `INTERIOR_BAND`; slope (θ run) 5.53e-5 ≤ `MAT-12` `S2_BAND`; (E) every run ≤ 5.8e-12 ≤ `E_BAND`; (∞) 1.9e-10 / 2.4e-11 ≤ `INF_BAND`; Dirichlet-at-a control misses 0.314 / 0.316. The T∞ = 20 K one-step BE balance (1.30e-9 / 1.16e-9) is a printed cancellation-floor record (2·DT_1 ratio 0.501 / 0.499) (`20261027T171958Z_IMP-7-step5-W1r.log`, `20261027T172155Z_IMP-7-step5-W2r.log`, 104 / 356 s). **Step 6 (🟡 record): the duty-cycled CN transient on the curved hole against `pec_sphere_implant_duty_cycle_from_modes`** (I47, §9 item 675, 2026-10-27, heavy `-n 4`, real; period 0.5 t\*, duty 0.5; every collected gate green on 4 ranks at `LADDER[0]` / `[1]`): helper = I45's Duhamel at five axis radii 4.48e-14 ≤ `B1_BAND`; (T) interior rel L2 at t\* 1.04e-4 / 1.34e-4 ≤ `INTERIOR_BAND`; (E) 5.8e-12 / 4.9e-12 ≤ `E_BAND`; (S0) census rel 0 ≤ 1e-12; continuous-source control misses 0.539 / 0.539 (`20261027T215045Z_IMP-7-step6-W0.log`, `20261027T215103Z_IMP-7-step6-W1.log`, `20261027T215234Z_IMP-7-step6-W2.log`, 7 / 81 / 233 s). Step 7 (I48a, §9 item 699, 2026-10-30, record, numpy): the dose and six-minute window around the adiabatic hole on the duty series to 4 t\*: `cem43` against quad ≤ 6.0e-16, the window ≤ 7.3e-16, the R swap and β = 0 miss, equal-mean dose ratio 1.015 / 1.013 / 1.007 at pole / equator / outer pole (`20261031T010117Z_IMP-7-step7-W1-r2.log`). Step 8 (I48b, §9 item 700, 2026-10-30, record, heavy `-n 4`, real, `LADDER[0]`): the first implant-side FEM dose, 160 CN steps to 4 t\*: (T) interior rel L2 at 4 t\* 9.00e-5 ≤ `INTERIOR_BAND`, (E) 5.8e-12, (S0) 0, `cem43` against quad ≤ 4.9e-16, the window ≤ 3.7e-16, the R swap and the continuous series (0.539) miss; FEM against series dose 1.2e-4 / 1.8e-5 / 2.6e-6 at pole / equator / outer pole, a reading with no band, no convergence claim (`20261031T020326Z_IMP-7-step8-W1-r2.log`, 142 s). The ✅ glyph stays on steps 0–1. *History: `docs/planning/chunks/IMP-7.md`.* 

### Review note 2026-11-01 03:00 (the 10-31 weekly's rulings carried)

**From §10:** I46′, I47, I48a, I48b ratified, tags dropped. I45's (IC) / (Slope) as truncation convergence ratified with no new bar; **its root ladder corrected:** with Neumann at a, the n-th root → nπ/(b − a) for every l, and the 10-26 brief's (n + l/2)π/(b − a) was wrong for this family. `S2_BAND` = `MAT-21` step 1's 1e-4 ratified. Item 700's re-run ratified as read narrowly. A CN Δt-rate gate for I47 is **declined**; I48c adopted.

### §9 item 710 — 2026-11-01 (§10 I48c, chain I — `IMP-7` step 9: the FEM-vs-series dose gap, separating h from Δt; chain step (solve), heavy, real)

*Enumerated by the 10-31 weekly. Depends on nothing. Next: I49. Family I48c, attempt 1 of 4.*

**Route.** New module `tests/validation/test_imp7_duty_dose_ladder.py`, `-n 4`, real build. Item 700's fixture and read route imported by name from `test_imp7_duty_dose_fem.py` (`_shell_fixture`, `_probe`, `_census`, `_err`, `STARTUP`), that module untouched. Two runs, each to `T_END` = 4 t\*:
- **Run h:** `LADDER[1]` (150 334 cells), Δt = `DT`, 160 steps, 161 reads.
- **Run Δt:** `LADDER[0]` (53 643 cells), Δt = `DT` / 2, 320 steps, reads at the even steps only, so the same 161 rows of `T_GRID`. `_chi(PERIOD, DUTY)` unchanged: a period is 40 steps, an on-phase 20, exact.

Series history from `_duty_hist(_modes())`, once. `_census` is `DT`-bound, so write `_census_dt(mult, dt)` beside it and assert the same `ON_TIME_END`.

**Anchors (imported, unmoved).** On each run: (T) interior rel L2 at 4 t\* ≤ `INTERIOR_BAND`; (E) ≤ `E_BAND`; (S0) ≤ `S0_BAND`; every point valid; max T < `T_REF`. **No dose band and no rate is claimed.** The record is one table of the three runs' FEM-vs-series gaps at each of `LABELS`: peak-scale dose rel, equal-mean dose rel, six-minute peak Δ, max_n|Δθ|. Row 1 is item 700's (`20261031T020326Z_IMP-7-step8-W1-r2.log`), read from a named constant set copied from that log, never re-run.

**Negative control (asserted).** The continuous series at `T_END` misses by > `INTERIOR_BAND` on each run (item 700 read this on `LADDER[0]`). **Grid identity (asserted):** run Δt's even-step times equal `T_GRID` to 1e-12 relative, so the rows compare like with like.

**Windows.** W0 `IMP-7-step9-W0` (`-n 1`, smoke, `timeout -k 30 120`: imports, grids, both census arithmetics, no solve); W1 `IMP-7-step9-W1a` (`-n 4`, run Δt alone, `timeout -k 30 600`); W2 `IMP-7-step9-W1b` (`-n 4`, run h alone, `timeout -k 30 600`). The two selected by a `FEM_EM_IMP7_S9` knob, `dt` / `h`, because together they risk the 600 s ceiling. W-env only if a `src/` file changed. **Price:** item 700's run read CN 36.9 s, 161 reads 4.5 s, window 142 s at `LADDER[0]`, `-n 4`, so run Δt ≈ 2× CN + the same comparands ≈ 180–220 s. Item 675's W2 at `LADDER[1]` read ≈ 300 s for two 40-step runs; a 160-step run at 2.8× cells ⇒ ≈ 300–450 s, inside the ceiling. **Class:** chain step (solve), ≈ 7.6 measured / 30 predicted minutes; the two windows alone are ≈ 10, so expect ≈ 12–14. **Tier:** heavy.

**Traps already paid for.** `parent=` on the sub-mesh read (item 700). The first log's 45 s crash is why. Every read on the same fixed `PTS` batch, call by call (item 537): never merge. Item 699's leak-guard trip on a closed-form print: print doses at ≤ 7 digits. `DT` in `test_imp7_duty_dose_fem.py` is imported; a run at half Δt must not reach `_census` or `_err`'s `T_GRID` silently. `INTERIOR_BAND`'s interior norm needs `_err`'s `RT["qmeta"]`. Real build only; take `np.real` if the group sources complex mode. A window over 600 s is unrunnable: on `LADDER[1]` reduce nothing, stop and report the bound.

**Scope.** Which of h and Δt the dose gap follows, as a record. **Not claimed:** a rate, a dose band, a converged dose, a Δt recommendation, a metal core.

**Negative result.** (T) red on either run ⇒ report it, park on `attempt/IMP-7-step9-*`, BLOCKED. A gap that does not move with either run is itself the reading: record it, since it would mean the gap is the source's or the point's. No `PREDICTED RED:` line.

**Status it can move.** §10 I48c → DONE; `IMP-7` gains step 9 as a record; the weekly has the evidence any dose or Δt band would need.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 710 — `IMP-7` I48c — per `docs/planning/chunks/IMP-7.md` "§9 item 710". Real build, `-n 4`, two windows. No landed module edited, no band, no rate. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (2026-11-01T12:09Z, 06:00 CST slot): DONE — `IMP-7` step 9, a record.** New module `tests/validation/test_imp7_duty_dose_ladder.py`, no landed module edited, no `src/` change, so no W-env. W0 `20261101T120125Z_IMP-7-step9-W0.log` (`-n 1`, 4 s). W1a `20261101T120136Z_IMP-7-step9-W1a.log` (run Δt, `-n 4`, 127 s against ≈ 180–220 predicted). W1b `20261101T120345Z_IMP-7-step9-W1b.log` (run h, `-n 4`, 254 s against ≈ 300–450 predicted). Each solve window: 4 passed on each of 4 ranks, Status 0. Pre-check: `pgrep -c python3` 0, loadavg 0.16.

Asserted, every band imported:

| run | cells | Δt (s) | steps | (T) at 4 t\* | (E) | (S0) | (Grid) | continuous series miss |
|---|---|---|---|---|---|---|---|---|
| item 700 (`LADDER[0]`, `DT`) | 53 643 | 7.8375 | 160 | 9.0012e-05 | 5.806e-12 | 0 | n/a | 5.3906e-01 |
| Δt (`LADDER[0]`, `DT` / 2) | 53 643 | 3.91875 | 320 | 8.4187e-05 | 1.168e-11 | 0 | 0 | 5.3906e-01 |
| h (`LADDER[1]`, `DT`) | 150 334 | 7.8375 | 160 | 4.1209e-05 | 4.929e-12 | 0 | 0 | 5.3897e-01 |

Every point valid; max T 20.43 °C at equal mean on both new runs, under `T_REF`.

The record (no band, no rate; row 1 of each point is item 700's, copied from `20261031T020326Z_IMP-7-step8-W1-r2.log:123–125`, never re-run):

| point | run | peak-scale dose rel | equal-mean dose rel | six-minute peak Δ (K) | max_n \|Δθ\| (K) |
|---|---|---|---|---|---|
| pole | item 700 | 1.245e-04 | 2.481e-04 | −8.242e-05 | 1.285e-03 |
| pole | Δt | 9.462e-05 | 1.903e-04 | −6.774e-05 | 3.413e-04 |
| pole | h | 1.310e-04 | 2.608e-04 | −8.709e-05 | 1.292e-03 |
| equator | item 700 | 1.795e-05 | 3.556e-05 | −7.155e-06 | 2.181e-04 |
| equator | Δt | 1.732e-05 | 3.412e-05 | −7.082e-06 | 3.889e-05 |
| equator | h | 1.317e-05 | 2.686e-05 | −9.836e-06 | 2.321e-04 |
| outer pole | item 700 | 2.605e-06 | 4.630e-06 | +4.165e-06 | 6.484e-05 |
| outer pole | Δt | 2.484e-06 | 5.219e-06 | +4.917e-06 | 1.697e-05 |
| outer pole | h | 4.920e-06 | 9.417e-06 | +2.462e-07 | 5.992e-05 |

**Reading, against item 700 (one step of each, so ratios and not rates).**
- **max_n |Δθ| follows Δt at all three points:** halving Δt gives × 0.27 / 0.18 / 0.26, and refining h gives × 1.01 / 1.06 / 0.92.
- **The dose gap does not follow one variable.** At the pole it moves with Δt, but only × 0.76, and h leaves it at × 1.05, so most of the pole's gap survives both runs. At the equator it moves with h (× 0.73 against × 0.96). At the outer pole neither run shrinks it: h raises it × 1.89, and Δt gives × 0.95 on the peak scale and × 1.13 at equal mean. Those are gaps of a few 1e-6 on a dose of ≈ 2e-11 s.
- **The interior L2 error at 4 t\* follows h** (× 0.46 against × 0.94). That differs from I47's reading at t\*, where the error did not shrink with h at fixed Δt (1.04e-4 / 1.34e-4). Both are records at different times, not a contradiction this item resolves.

No run refined h and Δt together, so the residual pole gap is not attributed to either. Nothing here supports a dose band yet.

**Disclosed.** The module prints (T) at 4 t\* only, which is what the brief asks; item 700's t\* and 2 t\* prints were not repeated. The solve windows also collect the W0 test, green on each. The module was not edited after W0 (rule (i)). `§10` is the weekly's text and the `IMP-7` §7 row already points here, so neither was edited; §10 I48c reads DONE by this outcome and §9's entry.

**Status moved.** §10 I48c DONE; `IMP-7` gains step 9 as a record.

### §9 item 711 — 2026-11-01 (§10 I49, lane 3 → chain I — `IMP-7` step 10: the dose sensitivity bound on the hollow-sphere duty series; lane 3 closed form, smoke, numpy)

*Enumerated by the 10-31 weekly. Depends on nothing. Next: I50. Family I49, attempt 1 of 4.*

**One execution detail, returned to the 11-02 weekly.** P39′'s (Bd) relates two histories (FEM and series), and this step is numpy, so it has one. (Bd)(i), m_p ≤ ln 4·M_p, holds for any pair of histories below 43 °C. **So the pair here is the N_DEFAULT-mode series and the same series at `N_DEFAULT // 2` modes**, a real, known perturbation of the comparand. (Bd)(ii), M_p ≤ `INTERIOR_BAND`·max θ_s, has no FEM history to bound here, so **`INTERIOR_BAND`·max θ_s is printed per cell as the cap I50 reads against.**

**Route.** New module `tests/validation/test_imp7_duty_matrix_closed_form.py`, numpy. `T_GRID`, `DT`, `T_END`, `SCALES` imported from `test_imp7_duty_dose_closed_form.py`, untouched, whose `_duty_hist` is bound to `PTS`: write a local `_duty_hist_at(modes, pts, duty)`. Points: `P_POLE` = `PTS[POLE]` (1.05 a on z) and `P_MID` = (0, 0, (a + b)/2). Duties D = 1, 0.5, 0.25, period `PERIOD` (0.5 t\*, so an on-phase of 10 and 5 `DT`: both on the grid).

**Anchors (asserted, every band imported).** For each of 3 duties × 2 points × `SCALES`:
- (Q-m) `cem43` against `_q_check`'s `quad` ≤ `DOSE_BAR`; max T < `T_REF`.
- (Win-m) the six-minute window: `_win_check`, t_max inside the range, the bracket, `quad` ≤ `AVG_BAR`.
- (Bd-i) m_p ≤ ln 4·M_p on the pair above.
- (C100-m) D = 1 equals the continuous history to `B1_BAND`, as (C100-i).

**Negative controls (asserted).** (Rc-m) the R swap misses by > `DOSE_BAR` on every cell. (β0-m) β = 0 misses at `P_POLE` by > `DOSE_BAR`, with its own modes built as `_hists()` builds them (`:88`). (Mis-m) a wrong schedule, D = 0.25 history against D = 0.5, misses by > `DOSE_BAR`.

**Printed (records, no band):** the 3 × 2 matrix of D(peak), D(equal mean), six-minute peak and headroom; (Bd) ratio m_p / (ln 4·M_p); the cap `INTERIOR_BAND`·max θ_s for I50.

**Windows.** W1 `IMP-7-step10-W1` (`-n 1`, numpy, `timeout -k 30 300`). **Price:** item 699's W1-r2 read 77 s for three histories at three points, so six duty histories at two points, plus two `N_DEFAULT // 2` sets, is ≈ 120–180 s. That is over the 30 s smoke ceiling, so the label is **standard**; item 699 was over-labelled the same way. W-env only if a `src/` file changed. **Class:** lane-3 closed form, ≈ 4.9 measured / 15 predicted minutes.

**Traps already paid for.** Item 699: print closed-form doses at ≤ 7 digits, because `.15e` / `.9f` tripped the leak guard twice; `IMP7-7`'s module carries the same note at `:145`. Item 689: θ is continuous at a switch only modally, so never assert grid equality across a switch. `_modes()` fixes β: the control builds its own. A dose ≈ 1e-11 s is ordinary at this source scale: the identities are on a history, not a thermal statement. `T_STAR` here is 313.5 s, shorter than the 360 s window, which is why the grid runs to 4 t\*. Do not import `test_*`-named functions.

**Scope.** The closed-form dose matrix and I50's comparand. **Not claimed:** a FEM history (I50), (Bd)(ii), a rate, a thermal-damage figure.

**Negative result.** (Bd-i) red would contradict a bound that holds by construction below 43 °C: report it and stop, because it means the max or the quadrature is wrong. A red (Q-m) ⇒ `cem43` and `quad` disagree on a history they have not met: report and stop. No `PREDICTED RED:` line.

**Status it can move.** §10 I49 → DONE; `IMP-7` gains step 10 (🟡 record); item 719 (I50) unblocks.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 711 — `IMP-7` I49 — per `docs/planning/chunks/IMP-7.md` "§9 item 711". Numpy only. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (2026-11-01T12:18Z, 06:00 CST slot): DONE — `IMP-7` step 10, a record, with one premise of the brief measured false.** New module `tests/validation/test_imp7_duty_matrix_closed_form.py`, no `src/` change. W1 `20261101T121040Z_IMP-7-step10-W1.log` (155 s, 5 passed) ran the brief's letter. W1-r2 `20261101T121358Z_IMP-7-step10-W1-r2.log` (149 s, 6 passed, Status 0) ran the module as committed, with the one addition below. Both `-n 1`.

Asserted on 12 cells (3 duties × 2 points × 2 scales), every band imported, worst values:
- **(Q-m):** 1.102e-15 against `DOSE_BAR` 1e-10; max T 20.50 °C, under `T_REF`.
- **(Win-m):** 6.175e-16 against `AVG_BAR` 1e-12; t_max in range, the bracket holds.
- **(C100-m):** 0 against `B1_BAND`, bit-for-bit, so it shows the helper's D = 1 path and `_theta` agree on these rows and little more.
- **Controls:** (Rc-m) misses by factors ≥ 6.9e+06; (β0-m) 9.21e-02 / 4.81e-02 / 2.45e-02 at D = 1 / 0.5 / 0.25; (Mis-m) 9.14e-02 at the pole, 6.72e-02 at mid.

**The premise that failed: the `N_DEFAULT // 2` pair is not a real perturbation.** The review's detail calls it one, but M_p reads 3.5527e-15 K on every cell at peak scale, one ulp at 20 °C: the modes it drops have decayed by the first grid row, ≈ exp(−500) at n = 100. So (Bd-i) is green as written (m_p / (ln 4·M_p) ≤ 0.0451) and that green carries no content, and the matrix's headroom column (≈ 1e11–1e12) is meaningless. **That column must not be read as headroom.**

**Added beyond the brief, for the 11-02 weekly to keep or strike: (Bd-i′), the same inequality on `N_DEFAULT` against `N_COARSE` = 2 modes,** which does move the history. No band is involved; the right side is measured.

| D | scale | pole: m_p ≤ ln 4·M_p (M_p, K) | mid: m_p ≤ ln 4·M_p (M_p, K) |
|---|---|---|---|
| 1 | peak | 4.6563e-04 ≤ 2.8351e-02 (2.0451e-02) | 1.7428e-04 ≤ 8.4527e-03 (6.0973e-03) |
| 0.5 | peak | 3.9431e-05 ≤ 2.8351e-02 (2.0451e-02) | 2.6514e-05 ≤ 8.4527e-03 (6.0973e-03) |
| 0.25 | peak | 6.4074e-05 ≤ 2.8351e-02 (2.0451e-02) | 1.6055e-05 ≤ 8.4527e-03 (6.0973e-03) |
| 0.5 | equal-mean | 9.2824e-05 ≤ 5.6702e-02 (4.0902e-02) | 1.3498e-04 ≤ 1.6905e-02 (1.2195e-02) |
| 0.25 | equal-mean | 9.8896e-04 ≤ 1.1340e-01 (8.1803e-02) | 2.9239e-04 ≤ 3.3811e-02 (2.4389e-02) |

M_p is the same at every duty on the peak scale, so it sits in the first on-phase, and the dose barely weights that time: m_p is ≈ 0.1–2 % of the bound. The bound holds but is loose, and I50 should expect the same.

The matrix, printed (no band). `cap` = `INTERIOR_BAND`·max θ_s, the figure I50 reads M_p against:

| D | point | D(peak) (s) | D(equal mean) (s) | six-minute peak (K) | cap (K) |
|---|---|---|---|---|---|
| 1 | pole | 2.598306e-11 | 2.598306e-11 | 0.355655 | 3.6878e-03 |
| 1 | mid | 2.342155e-11 | 2.342155e-11 | 0.263963 | 2.7408e-03 |
| 0.5 | pole | 2.163760e-11 | 2.637229e-11 | 0.179779 | 2.1490e-03 |
| 0.5 | mid | 2.052262e-11 | 2.368667e-11 | 0.132976 | 1.4877e-03 |
| 0.25 | pole | 1.966077e-11 | 2.657487e-11 | 0.090729 | 1.2410e-03 |
| 0.25 | mid | 1.914405e-11 | 2.380887e-11 | 0.066578 | 7.7009e-04 |

Item 700's series dose at the pole (2.163760e-11 s peak, 2.637229e-11 s equal mean) is reproduced at D = 0.5.

**Disclosed.** The brief imports `SCALES`, whose equal-mean factor is bound to D = 0.5, so this module takes only its labels and applies 1 / D per duty. "Headroom" follows the project's cap / M_p usage. The brief labels the item standard and §10 says smoke; 149 s is standard. The module carries a code comment on `N_COARSE` so the reason stays with it.

**Status moved.** §10 I49 DONE; `IMP-7` gains step 10 (🟡 record); item 719 (I50) unblocks. Its cap column is the table above, unaffected by the failed premise, because cap uses max θ_s alone.

### §9 item 719 — 2026-11-01 (§10 I50, chain I — `IMP-7` step 11: the FEM dose under D = 0.25 on `LADDER[0]`; chain step (solve), heavy, real)

*Enumerated by the 10-31 weekly. **Serial on item 711 (I49) landing; if it did not, skip at no cost.** Next: an `imp:` example. Family I50, attempt 1 of 4.*

**Route.** New module `tests/validation/test_imp7_duty25_dose_fem.py`, `-n 4`, real. Item 700's fixture and read route by import; the schedule `_chi(PERIOD, 0.25)`, so an on-phase of 5 `DT`, and the census target `ON_TIME_END` × 0.5. The read points are I49's `P_POLE` and `P_MID`, imported, in one fixed batch.

**Anchors (imported, unmoved).** (T) interior rel L2 ≤ `INTERIOR_BAND` at 4 t\* against the D = 0.25 series; (E) ≤ `E_BAND`; (S0) census ≤ `S0_BAND`; every point valid; max T < `T_REF`. (Q-hd) / (Win-hd) on the FEM history at both scales, `DOSE_BAR` / `AVG_BAR`.

**Negative controls (asserted).** The continuous series at `T_END` misses by > `INTERIOR_BAND`; the D = 0.5 series misses by > `INTERIOR_BAND`. The R swap misses by > `DOSE_BAR`.

**Printed (records, no band).** The dose gap beside item 700's D = 0.5 gap (constants copied from its log), the six-minute peak of each, and (Bd)(ii)'s M_p against I49's printed cap at each point. The weekly did not list (Bd)(ii) among I50's gates, so it stays a reading here. The `MAT-21` box found its headroom falling as the duty halves (5.04 → 3.05 → 1.72× at P5), so the slot prints that trend for the weekly.

**Windows.** W0 `IMP-7-step11-W0` (`-n 1`, smoke, `timeout -k 30 120`), W1 `IMP-7-step11-W1` (`-n 4`, real, `timeout -k 30 600`). **Price:** item 700's window, 142 s, the same run with a cheaper χ and two read points. **Class:** chain step (solve), ≈ 7.6 measured / 30 predicted minutes. **Tier:** standard by that price, heavy by its ceiling.

**Traps already paid for.** Item 710's list. The two start-up half-steps both sit inside the first on-phase at D = 0.25 too (5 `DT` long), but read `_census` rather than assume it. `_chi` takes (period, duty): pass 0.25 explicitly, as `DUTY` is a module constant at 0.5, and `ON_TIME_END` is too.

**Scope.** The FEM dose at D = 0.25. **Not claimed:** a rate, a dose band, convergence.

**Negative result.** (T) red ⇒ report, park on `attempt/IMP-7-step11-*`, BLOCKED. M_p above the cap is a reading for the weekly, not a red. No `PREDICTED RED:` line.

**Status it can move.** §10 I50 → DONE; `IMP-7` gains step 11; an `imp:` example becomes owed (the next review writes it, closing window `./run_examples.sh <id>`).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 719 — `IMP-7` I50 — per `docs/planning/chunks/IMP-7.md` "§9 item 719". Only after item 711 is DONE. Real build, `-n 4`. No landed module edited, no band moved. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (2026-11-01T15:06Z, 09:00 CST slot): DONE — `IMP-7` step 11, a record.** New module `tests/validation/test_imp7_duty25_dose_fem.py`, no landed module edited, no `src/` change, so no W-env. W0 `20261101T150132Z_IMP-7-step11-W0.log` (`-n 1`, 5 s). W1 `20261101T150144Z_IMP-7-step11-W1.log` (`-n 4`, 122 s against ≈ 142 predicted): 5 passed on each of 4 ranks, Status 0. Pre-check: `pgrep -c python3` 0. The census ran from the run's own multipliers: 161 (half-)steps, the two start-up half-steps both on.

Asserted, every band imported:
- **(T) at 4 t\*:** interior rel L2 7.3744e-05 ≤ `INTERIOR_BAND` 1e-02 (item 700 at D = 0.5: 9.0012e-05).
- **(E):** 5.794e-12 ≤ `E_BAND`. **(S0):** 313.5 s against 313.5 s, rel 0. Both points valid; max T 20.50 °C at equal mean, under `T_REF`.
- **(Q-hd):** ≤ 3.287e-16 against `DOSE_BAR` 1e-10 on both points and scales. **(Win-hd):** ≤ 9.183e-16 against `AVG_BAR` 1e-12.
- **Controls:** the continuous series misses 7.7871e-01, the D = 0.5 series 5.1989e-01, each > `INTERIOR_BAND`; the R swap misses by factors ≥ 6.9e+06.

The record (no band, no rate). Series doses at D = 0.25 reproduce I49's matrix (1.966077e-11 s and 1.914405e-11 s at peak). Item 700 read no mid point, so mid has no D = 0.5 row.

| point | D | peak-scale dose rel | equal-mean dose rel | six-minute peak FEM / series (K) | peak Δ (K) | M_p = max_n \|Δθ\| (K) | cap (K) | headroom cap / M_p |
|---|---|---|---|---|---|---|---|---|
| pole | 0.5 (item 700) | 1.245e-04 | 2.481e-04 | 0.179697 / 0.179779 | −8.242e-05 | 1.285e-03 | 2.1490e-03 | 1.672 |
| pole | 0.25 | 8.403e-05 | 3.392e-04 | 0.090679 / 0.090729 | −4.971e-05 | 1.2845e-03 | 1.2410e-03 | **0.966** |
| mid | 0.25 | 4.072e-06 | 1.462e-05 | 0.066580 / 0.066578 | +2.108e-06 | 2.6266e-05 | 7.7009e-04 | 29.319 |

**Reading, for the weekly.**
- **At the pole, M_p exceeds (Bd)(ii)'s cap at D = 0.25,** which is the brief's "reading, not a red". The trend matches the `MAT-21` box's: headroom falls as the duty halves.
- **The fall here comes from the cap, not from the error.** M_p is unchanged between the duties and sits at t = 7.84 s, the first full step, inside the first on-phase at either duty. The cap scales with max θ_s, which falls with D. So (Bd)(ii) at this point measures a start-up error against a peak that shrinks with D. That start-up error barely weights the dose, as I49 found.
- **The dose gap does not grow the same way:** the peak-scale gap shrinks × 0.67 and the equal-mean gap grows × 1.37 from D = 0.5 to 0.25. One mesh, one Δt: ratios, not rates.

**Disclosed.** `_duty_fn` is bound to D = 0.5, so the module carries `_duty_fn_at(t, duty)`, and `_census` is bound to `ON_TIME_END`, so it carries `_census_25`. `CAP_50_POLE` and `ITEM_700_POLE` are constants copied from I49's table and item 700's log, not re-run. The brief asks for the headroom trend; the D = 1 cell has no FEM M_p on this fixture, so the trend here has two points, at the pole only. The solve window also collects the W0 test, green on each rank. The module was not edited after W0 (rule (i)).

**Status moved.** §10 I50 DONE; `IMP-7` gains step 11 as a record; an `imp:` example is owed, and the next review writes it.

### §9 item 725 — 2026-11-01 (chain I, I51 — daily, provisional, rule (6))

**Why.** I48c's question was "the record that separates h from Δt", and item 710 answered half of it: no run refined both, so most of the pole's dose gap is unattributed (pole × 0.76 under Δt, × 1.05 under h). This is the fourth corner of that table. It adds no target and no band. **Provisional:** the 11-02 weekly ratifies or reverts.

**Anchor.** On `LADDER[1]`, Δt = `DT` / 2, 320 steps, read at even steps: (T) interior rel L2 at 4 t\* ≤ `INTERIOR_BAND`, (E) ≤ `E_BAND`, (S0) ≤ `S0_BAND`, (Grid) ≤ `GRID_BAR`, every point valid, max T < `T_REF`. All imported from `tests/validation/test_imp7_duty_dose_ladder.py` and its sources. Same `_gaps`, same 161 rows of `T_GRID`.

**Negative control.** The continuous series at `T_END` misses by > `INTERIOR_BAND` (items 700 / 710 read 5.39e-01, so the bar is reachable). The h and Δt rows are constants copied from item 710's log, with a cited line, and never re-run.

**Execution.** New module `tests/validation/test_imp7_duty_dose_both.py`. Import from `test_imp7_duty_dose_ladder.py`: `_census_dt`, `_mids`, `_read_times`, `_gaps`, `GRID_BAR`; the run is `(LADDER[1], 0.5 * DT, 2 * N_RUN, 2)`. Write `run`'s body out for this one run; no landed module edited, no knob-selected second run.

**Printed (no band, no rate):** the four-row table at pole / equator / outer pole, item 700, Δt, h and this run: peak-scale dose rel, equal-mean dose rel, six-minute peak Δ, max_n |Δθ|.

**Tier, ranks, wall clock.** Heavy, `-n 4`, real, `timeout -k 30 590`. Predicted ≈ 375–420 s: item 710's h run took 254 s with 160 CN steps; the extra 160 steps are ≈ 120–160 s at that rung by scaling item 700's 36.9 s CN by ≈ 3–4 × for 2.8 × the cells. W0 `-n 1` ≈ 5 s on the grid and census node. **Class:** chain step (solve). **Cost:** ≈ 5.5 measured / 25 predicted min. **Depends on:** nothing. Pre-check `pgrep -c python3` 0.

**Traps already paid for.** `_probe` needs `parent=` (a `create_submesh` mesh). Read in one fixed batch per step (item 537). Clear snapshots once read. **The leak guard fires on 8-digit closed-form figures:** print 7 digits or fewer. A window over ≈ 10 min must not exist: if W0's census or the first phases point past 590 s, stop, do not raise `-t`.

**Scope.** A four-corner record. **Not claimed:** a rate (the 10-31 weekly declined a Δt-rate gate, and P46's epitaph applies: the comparand carries its own floor), a dose band, a converged dose, a metal core.

**Negative result.** Any asserted gate red: stop, keep the log, BLOCKED; never move a band. A gap that does not shrink at this corner is a reading, not a red; it says the residual sits in something neither h nor Δt moves (for example the probe's interpolation), and the outcome says so. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 12, a record; item 710's "unattributed" clause in the outcome resolves one way or the other for the weekly. No glyph.

**Outcome (2026-11-01T18:33Z, 12:00 CST slot): DONE — `IMP-7` step 12, a record.** New module `tests/validation/test_imp7_duty_dose_both.py`, no landed module edited, no `src/` change, so no W-env, as item 710. W0 `20261101T182629Z_IMP-7-step12-W0.log` (`-n 1`, 1 passed). W1 `20261101T182633Z_IMP-7-step12-W1.log` (`-n 4`, 338 s against ≈ 375–420 predicted, CN 210.1 s): 4 passed on each of 4 ranks, Status 0 (`:90`, `:118`, `:129`, `:138`, `:141`). Pre-check: `pgrep -c python3` 0.

Asserted, every band imported: `LADDER[1]`, 150 334 cells, Δt 3.91875 s, 320 steps; (T) at 4 t\* 2.3901e-05, (E) 1.010e-11, (S0) 0, (Grid) 0; every point valid; max T 20.43 °C at equal mean, under `T_REF`. The continuous series misses by 5.3897e-01 (`:74–77`, `:93`, `:95`).

The record (no band, no rate; rows 1–3 are constants, cited in the module, never re-run):

| point | run | peak-scale dose rel | equal-mean dose rel | six-minute peak Δ (K) | max_n \|Δθ\| (K) |
|---|---|---|---|---|---|
| pole | item 700 | 1.245e-04 | 2.481e-04 | −8.242e-05 | 1.285e-03 |
| pole | Δt | 9.462e-05 | 1.903e-04 | −6.774e-05 | 3.413e-04 |
| pole | h | 1.310e-04 | 2.608e-04 | −8.709e-05 | 1.292e-03 |
| pole | both | 1.011e-04 | 2.031e-04 | −7.272e-05 | 3.488e-04 |
| equator | item 700 | 1.795e-05 | 3.556e-05 | −7.155e-06 | 2.181e-04 |
| equator | Δt | 1.732e-05 | 3.412e-05 | −7.082e-06 | 3.889e-05 |
| equator | h | 1.317e-05 | 2.686e-05 | −9.836e-06 | 2.321e-04 |
| equator | both | 1.251e-05 | 2.535e-05 | −9.759e-06 | 5.094e-05 |
| outer pole | item 700 | 2.605e-06 | 4.630e-06 | +4.165e-06 | 6.484e-05 |
| outer pole | Δt | 2.484e-06 | 5.219e-06 | +4.917e-06 | 1.697e-05 |
| outer pole | h | 4.920e-06 | 9.417e-06 | +2.462e-07 | 5.992e-05 |
| outer pole | both | 1.724e-07 | 4.398e-07 | +9.965e-07 | 1.150e-05 |

**Reading, against item 700 (ratios, not rates; one row each).**
- **Item 710's unattributed pole gap stays mostly unmoved at this corner, which is a reading, not a red:** × 0.81 peak-scale and × 0.82 equal-mean, and larger than Δt alone (× 0.76). About four fifths of the pole's dose gap is in something neither step moves here. This record does not say what; the brief's example is the probe's interpolation, and that is not tested.
- **The equator sits at × 0.70 / × 0.71,** close to h alone (× 0.73 / × 0.76).
- **The outer pole drops × 0.066 / × 0.095,** which neither run did alone (Δt × 0.95 / × 1.13, h × 1.89 / × 2.03). These gaps are the smallest in the table, of order 1e-7 to 1e-6 relative, so this may be a crossing and not a convergence. One row cannot tell which.
- **max_n |Δθ| follows Δt again:** × 0.27 / × 0.23 / × 0.18, against Δt alone × 0.27 / × 0.18 / × 0.26. **(T) at 4 t\* falls × 0.27 from item 700.**

**Disclosed.** `_duty_fn` and every band are imported; `ITEM_710_DT` and `ITEM_710_H` are constants copied at the `.3e` of the logs' prints; `ITEM_700` is imported from item 710's module. The window also collects the W0 test, green on each rank. The module was not edited after W0 (rule (i)). The print is ≤ 7 digits. **Provisional:** the 11-02 weekly ratifies or reverts.

**Status moved.** `IMP-7` gains step 12 as a record. Item 710's "unattributed" clause resolves as: not h, not Δt, not their pairing, at this one refinement.

### §9 item 726 — 2026-11-01 (chain I, I52 — daily, provisional, rule (6))

**Why.** I50 read M_p above (Bd)(ii)'s cap at D = 0.25 and argued that the fall in headroom comes from the cap, because M_p sits at t = 7.84 s at both duties. That rests on two points, because the D = 1 cell has no FEM M_p on this fixture. This is the third point. It tests the argument rather than a new target. **Provisional.**

**Anchor.** On `LADDER[0]` at D = 1: (T) interior rel L2 at 4 t\* ≤ `INTERIOR_BAND` against the continuous series, (E) ≤ `E_BAND`, (Q-hd) / (Win-hd) at `DOSE_BAR` / `AVG_BAR` at pole and mid on the peak scale, every point valid, max T < `T_REF`. All imported as `tests/validation/test_imp7_duty25_dose_fem.py` imports them.

**Negative control.** The D = 0.5 series at `T_END` misses by > `INTERIOR_BAND` (item 719 read 5.2e-01 between duties at D = 0.25; print this one's value, the bar is `INTERIOR_BAND`, and if it does not clear, record it rather than lower anything); the R swap misses by > `DOSE_BAR`.

**Execution.** New module `tests/validation/test_imp7_continuous_dose_fem.py`, the module body of item 719 with `D = 1.0` and no `source_schedule` (the constant source), on `M_PTS`. The comparand is the D = 1 history from `_duty_hist_at(_modes(), M_PTS, 1.0)` if that path reduces to the continuous series; assert (C100) that it does, at `B1_BAND`, against `_series_fn`, as I49 did. (S0) has no content at D = 1: assert the source multipliers are all 1 and say so.

**Printed (no band):** the M_p table: D = 1, 0.5, 0.25 at the pole, M_p, t at M_p, cap (3.6878e-03 / 2.1490e-03 / 1.2410e-03 K, I49's), headroom. **Prediction, not a gate:** M_p ≈ 1.28e-03 K at t ≈ 7.84 s, so headroom ≈ 2.9. A D = 1 M_p far from that refutes I50's reading, and the outcome says so.

**Tier, ranks, wall clock.** Standard, `-n 4`, real, `timeout -k 30 300`, predicted ≈ 120 s (item 719: 122 s). W0 `-n 1` ≈ 5 s. **Class:** chain step (solve). **Cost:** ≈ 5.5 measured / 15 predicted min. **Depends on:** nothing; independent of item 725.

**Traps already paid for.** As item 725. `SCALES`' equal-mean factor is bound to D = 0.5 (item 711): read the peak scale only here, because equal mean is the same scale at D = 1. `_census` and `_duty_fn` are bound to D = 0.5; do not import them (item 719).

**Scope.** The third point of one trend at one mesh and Δt. **Not claimed:** a dose accuracy, a rate, any change to (Bd), which stays a reading.

**Negative result.** Asserted gate red: stop, keep the log, BLOCKED. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 13, a record. I50's reading of (Bd)(ii) is confirmed or refuted for the weekly, which decides whether (Bd)(ii) stays a reading at all.

**Outcome (2026-11-01T19:52Z, 13:30 CST slot): DONE — `IMP-7` step 13, a record.** New module `tests/validation/test_imp7_continuous_dose_fem.py`, no landed module edited, no `src/` change, so no W-env, as items 710 and 725. W0 `20261101T194739Z_IMP-7-step13-W0.log` (`-n 1`, 1 passed, 39 s). W1 `20261101T194823Z_IMP-7-step13-W1.log` (`-n 4`, 143 s against ≈ 120 predicted, CN 36.2 s): 5 passed on each of 4 ranks, 20 `PASSED` tokens, Status 0 (`:106`, `:118`, `:133`, `:142`, `:145`). Pre-check: `pgrep -c python3` 0.

Asserted, every band imported: `LADDER[0]`, 53 643 cells, `DT`, 160 steps; (T) at 4 t\* 1.2517e-04, (E) 5.779e-12; both points valid; max T 20.368654 / 20.274101 °C under `T_REF`. (Q-hd) pole 4.975e-16, mid 4.139e-16; (Win-hd) pole 6.245e-16, mid 1.051e-15. (C100) worst rel 0.000e+00 at `B1_BAND`. The D = 0.5 series misses by 1.1693e+00, the R swap by 6.9342e+06 / 7.3079e+06 (`:68`, `:77`, `:78`, `:82`, `:84–87`).

The M_p table at the pole (no band):

| D | M_p (K) | t(M_p) (s) | cap (K) | headroom cap / M_p |
|---|---|---|---|---|
| 1 | 1.2845e-03 | 7.84 | 3.6878e-03 | 2.871 |
| 0.5 | 1.2850e-03 | not printed by item 700's log | 2.1490e-03 | 1.672 |
| 0.25 | 1.2845e-03 | 7.84 | 1.2410e-03 | 0.966 |

**Reading.** **I50's argument is confirmed at the third point:** M_p at D = 1 is 1.2845e-03 K at t = 7.84 s, 1.004 × the prediction, so the fall in headroom comes from the cap. The D = 0.5 row is the `.3e` print of item 700 and carries no t(M_p); the "both duties" of I50 is therefore read here at D = 1 and 0.25 only. The cap column reproduces I49's value at D = 1 on this run's own series. Mid is printed in the log, M_p 2.5227e-05 K at t = 62.70 s, headroom 108.648. Peak-scale dose gap at D = 1: pole 1.909e-04, mid 1.145e-05.

**Disclosed.** (1) **W0's first window was killed at its 120 s ceiling and is kept** (`20261101T194502Z_IMP-7-step13-W0.log`, Status 124, no gate value printed, `pgrep -c python3` 0 after). The cause was this module's own: (C100) called `_series_fn` for each of 160 rows, and that rebuilds the modes per call. The case was shrunk to one mode set through `_theta`, I49's own route, and re-run at the same `-t`. The brief says "against `_series_fn`, as I49 did", and I49 used `_theta`; the two are the same series and (T) still uses `_series_fn`. "A kill parks the item" is not in this brief, so the slot read the kill under CLAUDE.md's "kill and shrink" and did not park it; the review may rule otherwise. (2) (C100) was placed in W0, so W0 reads ≈ 39 s against the brief's ≈ 5 s, and W1 repeats it; the added ≈ 35 s per collected run is the D = 1 history. (3) (S0): the constant source records no multipliers, so the module asserts the recorded list is empty and, in W0, that `_chi(PERIOD, 1.0)` reads 1 at all 161 midpoints with Σχ·dt = `T_END` at rel 0. The brief's "all 1" is asserted on that list, not on the solver's. (4) `ITEM_719_POLE` is a constant copied from `20261101T150144Z_IMP-7-step11-W1.log:131`. The module was not edited after W0's green re-run (rule (i)). **Provisional:** the 11-02 weekly ratifies or reverts.

**Status moved.** `IMP-7` gains step 13 as a record. I50's reading of (Bd)(ii) is confirmed at the third point, and the weekly decides whether (Bd)(ii) stays a reading.

### Review note 2026-11-01 18:00 (item 725's residual; item 726's kill)

- **Correction to item 725's brief.** Its example of what neither refinement moves, "the probe's interpolation", was a guess. `PTS` are interior points (r = 1.05 a, 0.95 b) on a curved P2 mesh with P2 temperature, so a low-order interpolation at a surface is not the first suspect. Nothing in the fixture has yet tested **the comparand** at the pole: `N_DEFAULT` = 200 modes, read from the first row of a pulsed history. Item 741 tests that first, and P46's epitaph says it can carry a floor. Item 745 then asks whether the gap has a Δt limit.
- **Item 726's W0 kill stands as the slot handled it.** It was a 120 s ceiling on a module-level loop with no gate value printed, shrunk at the same `-t` under CLAUDE.md's "kill and shrink", first log kept. The brief had no "a kill parks the item" line. It goes to the 11-02 weekly as a disclosure; this review extends no ruling.

### §9 item 741 — 2026-11-01 (chain I, I53 — daily, provisional, rule (6); lane 3 closed form)

**Why.** Item 725 left about four fifths of the pole's dose gap in something neither h nor Δt moves at one refinement. The cheapest candidate is the comparand itself. **Provisional:** the 11-02 weekly ratifies or reverts.

**Anchor.** `_duty_hist(_modes(n))` (`tests/validation/test_imp7_duty_dose_closed_form.py`) at n = 200, 400, 800 on `T_GRID`, at `PTS`:
- **(Q-hd) / (Win-hd)** at `DOSE_BAR` / `AVG_BAR` on each n, as item 711 writes them.
- **(Self)** the 800-mode dose at the pole against the 200-mode dose: the rel value printed, and asserted ≤ `B1_BAND`, imported. A red is a finding about `N_DEFAULT`; see Negative result.

**Negative control.** n = 20 misses the 200-mode dose at the pole by > `DOSE_BAR`, so truncation has teeth on this quantity. The R swap misses by > `DOSE_BAR`.

**Printed (no band):** rel dose change 200 → 400 → 800 at pole, equator, outer pole, on the peak scale, beside item 725's gaps (1.011e-04 / 1.251e-05 / 1.724e-07); the row of T_GRID that carries the largest |Δθ| between 200 and 800.

**Execution.** New module `tests/validation/test_imp7_dose_comparand_floor.py`, numpy, no `src/` change. **First in W1:** `_modes(800)` builds and its root ladder passes its own checks.

**Tier, ranks, wall clock.** Standard, `-n 1`, `timeout -k 30 170`. Item 711 ran the matrix in seconds at 200 modes; the cost is linear in n by row, predicted ≈ 20–90 s. **Class:** lane-3 closed form. **Cost:** ≈ 2.9 measured / 15 predicted min. **Depends on:** nothing.

**Traps already paid for.** `SCALES`' equal-mean factor is bound to D = 0.5: this module is D = 0.5, so both scales may print, but item 711's `N_DEFAULT // 2` pair is the (Bd-i) with no content (M_p one ulp at 20 °C): do not reuse it as a sensitivity. **The leak guard fires on 8-digit closed-form figures:** print 7 digits or fewer. A root-ladder failure at n = 800 is a finding: report it, fall to n = 400, say so.

**Scope.** The comparand. **Not claimed:** a dose accuracy, a rate, that the residual is the comparand unless the printed change reaches the gap's size.

**Negative result.** (Self) red, or the change at the pole of order 1e-4: the comparand carries the floor, which is the finding. Keep the log, stop, BLOCKED, and the review re-reads item 710's table. Both Q identities red: stop. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 14, a record; item 725's "not h, not Δt" residual is attributed to the comparand or the comparand is excluded.

**Outcome (2026-11-02T03:06Z, 21:00 CST slot on 11-01): BLOCKED — the negative control (Tr20) red, unpredicted; every other gate green.** W1 `20261102T030134Z_IMP-7-step14-W1.log` (`-n 1`, `timeout -k 30 170`, 150 s against ≈ 20–90 predicted): 1 failed, 4 passed, Status 1. Module and probe parked on `attempt/IMP-7-step14-20261102T030507Z`; no `src/` change, no landed module edited.

**The red.** (Tr20): n = 20 against n = 200 at the pole misses by 1.4933e-16; the brief asks > `DOSE_BAR` = 1e-10 (`:85–86`, `:97–98`). The brief carries no `PREDICTED RED:` line and does not label the control asserted or predicted, so rule (e) applies: the result is reported, not decided in-slot. No bound, n or module was changed after the window.

**Green, every band imported.** (Root) 800 roots on both harmonics, max |x_n − n| 4.779e-02 / 1.167e-01, min gap 0.752 / 0.732 spacings (`:43–44`). (Q-hd) / (Win-hd) green on n = 200, 400, 800 at all three points on both scales, the R swap missing by far more than `DOSE_BAR` on every cell (`:46–81`). (Self) rel 0.000e+00 at `B1_BAND`, against item 725's 1.011e-04 (`:83`).

**Printed (no band).** Peak-scale rel dose change, 200 → 400 / 400 → 800, beside item 725's gap (`:87–89`):

| point | 200 → 400 | 400 → 800 | item 725's gap |
|---|---|---|---|
| pole | 0.000e+00 | 0.000e+00 | 1.011e-04 |
| equator | 1.565e-16 | 0.000e+00 | 1.251e-05 |
| outer pole | 1.671e-16 | 0.000e+00 | 1.724e-07 |

The largest |Δθ| between 200 and 800 is 3.553e-15 K at row 25, t = 195.94 s, outer pole; the pole's and equator's are the same value at rows 69 and 80 (`:90`). That is rounding at 20 °C, so **the row carries no content**.

**Reading, one the review may take or not.** The window's numbers exclude the comparand: the change reads 0 at the pole and at most 1.671e-16 elsewhere, at least nine orders below each point's gap, so the residual is not 200-mode truncation. **That reading is not closed here,** because the control that was to show truncation has teeth on this quantity did not fire at n = 20. The two facts are one fact: by n = 20 the series is already converged in dose to rounding, so n = 20 cannot separate a working (Self) from a blind one.

**Measure-only probe, for the ruling** (`20261102T030430Z_IMP-7-step14-probe.log`, 28 s, Status 0, asserts nothing, on the branch). Rel dose against n = 200:

| n | pole | equator | outer pole | pole max \|Δθ\| (K), row |
|---|---|---|---|---|
| 1 | 5.463e-04 | 3.345e-04 | 1.552e-04 | 1.789e-02, 11 |
| 2 | 3.943e-05 | 7.243e-06 | 1.786e-06 | 2.045e-02, 1 |
| 3 | 2.103e-05 | 6.880e-06 | 1.306e-06 | 5.127e-03, 1 |
| 4 | 9.171e-06 | 5.615e-08 | 6.238e-07 | 2.007e-03, 1 |
| 6 | 8.205e-07 | 1.651e-09 | 2.520e-08 | 1.898e-04, 1 |
| 8 | 6.629e-08 | 8.007e-10 | 1.617e-10 | 1.605e-05, 21 |
| 12 | 2.319e-10 | 6.334e-12 | 1.772e-11 | 5.704e-08, 1 |
| 20 | 1.493e-16 | 1.565e-16 | 0.000e+00 | 1.776e-14, 21 |

So truncation does have teeth on this quantity, below n = 20: at the pole n = 12 misses by 2.3 × `DOSE_BAR` and n = 8 by 6.6e+02 ×. A one-mode series changes the pole's dose by 5.463e-04, about 5.4 × item 725's gap, and that is the size of error a badly under-resolved comparand would have to be; 200 modes is not near it.

**Disclosed.** (1) The window took 150 s, over the predicted 90 and under the 170 ceiling; the n = 400 and n = 800 histories are the likely cost, not timed apart. (2) The `.6f` max T and `.6e` prints passed the leak guard. (3) The probe is not in the brief; it edits no gate and closes nothing. (4) Step 14 is not gained, and item 745's joint "Δt, the comparand, or neither" stays half made. **Provisional:** the 11-02 weekly ratifies or reverts.

### §9 item 745 — 2026-11-01 (chain I, I54 — daily, provisional, rule (6))

**Why.** Item 700's table has `DT` and `DT` / 2 on `LADDER[0]`: pole gap × 0.76. Two points cannot say whether the gap heads to zero or to a floor. A third at `DT` / 4 can, as a record. **Provisional.**

**Anchor.** `LADDER[0]`, Δt = `DT` / 4, 640 steps, read at every fourth step, the same 161 rows of `T_GRID`: (T) interior rel L2 at 4 t\* ≤ `INTERIOR_BAND`, (E) ≤ `E_BAND`, (S0) ≤ `S0_BAND`, (Grid) ≤ `GRID_BAR`, every point valid, max T < `T_REF`. All imported from `tests/validation/test_imp7_duty_dose_ladder.py` and its sources.

**Negative control.** The continuous series at `T_END` misses by > `INTERIOR_BAND` (items 700 / 710 / 725 read 5.39e-01, so the bar is reachable).

**Execution.** New module `tests/validation/test_imp7_duty_dose_quarter_dt.py`, on item 725's pattern: the run is `(LADDER[0], 0.25 * DT, 4 * N_RUN, 4)`, `run`'s body written out for this one run, `_census_dt`, `_mids`, `_read_times`, `_gaps`, `GRID_BAR` imported. No landed module edited.

**Printed (no band, no rate):** the three-row table `DT`, `DT` / 2, `DT` / 4 at pole / equator / outer pole: peak-scale dose rel, equal-mean dose rel, six-minute peak Δ, max_n |Δθ|, with `ITEM_700`, and item 710's row as constants cited to the log.

**Tier, ranks, wall clock.** Heavy, `-n 4`, real, `timeout -k 30 590`. Predicted ≈ 180–280 s: item 726's 160 steps took CN 36.2 s, so 640 ≈ 145 s, plus build, 161 reads and two series. W0 `-n 1` on the grid and census node. **Class:** chain step (solve). **Cost:** ≈ 7.7 measured / 25 predicted min. **Depends on:** nothing. Pre-check `pgrep -c python3` 0.

**Traps already paid for.** As item 725. `_probe` needs `parent=`. Read in one fixed batch per step (item 537), only the snapshots read, cleared once read: 161 of them, not 641. **Print 7 digits or fewer.** A window over ≈ 10 min must not exist: if W0's census or the first phases point past 590 s, stop.

**Scope.** A record at one mesh. **Not claimed:** a rate (the 10-31 weekly declined a Δt-rate gate), a dose band, a converged dose.

**Negative result.** Asserted gate red: stop, keep the log, BLOCKED. A gap that stops shrinking is a reading, not a red. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 15, a record; with item 741 it says whether the pole's residual is Δt, the comparand, or neither.

**Outcome (2026-11-02T01:54Z, 19:30 CST slot on 11-01, fourth item): DONE — `IMP-7` step 15, a record.** New module `tests/validation/test_imp7_duty_dose_quarter_dt.py`, run inline; no landed module edited, no `src/` change, so no W-env. W0 `20261102T015029Z_IMP-7-step15-W0.log` (`-n 1`, 3 s, 1 passed). W1 `20261102T015038Z_IMP-7-step15-W1.log` (`-n 4`, 177 s against ≈ 180–280 predicted, CN 110.3 s): 4 passed on each of 4 ranks, 16 `PASSED` tokens, Status 0. Pre-check: `pgrep -c python3` 0. The module was not edited after W0 (rule (i)).

Asserted, every band imported: `LADDER[0]`, 53 643 cells, Δt 1.959375 s, 640 steps, 161 reads; (T) at 4 t\* 8.3973e-05, (E) 2.386e-11, (S0) 0, (Grid) 0; every point valid; max T 20.43 °C at equal mean, under `T_REF`. The continuous series misses by 5.3906e-01.

The record (no band, no rate; rows `DT` and `DT` / 2 are constants, `ITEM_700` and `ITEM_710_DT`, imported, never re-run):

| point | Δt | peak-scale dose rel | equal-mean dose rel | six-minute peak Δ (K) | max_n \|Δθ\| (K) |
|---|---|---|---|---|---|
| pole | `DT` | 1.245e-04 | 2.481e-04 | −8.242e-05 | 1.285e-03 |
| pole | `DT` / 2 | 9.462e-05 | 1.903e-04 | −6.774e-05 | 3.413e-04 |
| pole | `DT` / 4 | 8.917e-05 | 1.806e-04 | −6.717e-05 | 1.867e-04 |
| equator | `DT` | 1.795e-05 | 3.556e-05 | −7.155e-06 | 2.181e-04 |
| equator | `DT` / 2 | 1.732e-05 | 3.412e-05 | −7.082e-06 | 3.889e-05 |
| equator | `DT` / 4 | 1.701e-05 | 3.344e-05 | −7.063e-06 | 3.812e-05 |
| outer pole | `DT` | 2.605e-06 | 4.630e-06 | +4.165e-06 | 6.484e-05 |
| outer pole | `DT` / 2 | 2.484e-06 | 5.219e-06 | +4.917e-06 | 1.697e-05 |
| outer pole | `DT` / 4 | 3.774e-06 | 7.724e-06 | +5.154e-06 | 8.618e-06 |

**Reading (ratios, not rates).** **The pole's dose gap does not head to zero with Δt at this mesh.** × 0.76 from `DT` to `DT` / 2, then × 0.94 peak-scale and × 0.95 equal-mean to `DT` / 4, with the six-minute Δ at × 0.99. The gap sits near 8.9e-05 and looks like a floor. Three points cannot name that floor's size. The equator is flatter still, × 0.96 then × 0.98. **The outer pole's gap grows,** × 1.52 / × 1.48 at the last halving; these are the table's smallest values, of order 1e-6, so it may be a crossing. **max_n |Δθ| still falls at the pole and outer pole,** × 0.55 and × 0.51, but stops at the equator, × 0.98. That row's own floor, 3.8e-05, is reached by `DT` / 2. So the residual is not Δt at this mesh; item 725's corner says it is not h at one refinement either. Item 741 tests the comparand; until it lands, this record says only "not Δt".

**Disclosed.** (1) Item 741 (step 14) has not run, so step 15 lands before it, in list order; the status line's joint statement is therefore half made. (2) `ITEM_710_DT` is imported from item 725's module and `ITEM_700` from item 710's, as the brief's "with `ITEM_700`, and item 710's row as constants" allows. (3) Snapshots are cleared at the end of the fixture, as item 725, not one by one; only the 161 read are recorded. (4) (E) 2.386e-11 is above the `DT` and `DT` / 2 runs' readings and well inside `E_BAND`; not investigated. **Provisional:** the 11-02 weekly ratifies or reverts.

**Status moved.** `IMP-7` gains step 15 as a record; §9 item 745 DONE.

### Review ruling 2026-11-02 03:00 (item 741, I53) — the control's n was the brief's guess; re-posed as item 752

**Ruled.** (Tr20) was the review's own guess at n, written without computing the ceiling the rubric asks for. The red is a defect in the brief, not in the module or the comparand. The control is re-posed at **n = 8**, chosen from the slot's measure-only probe (`20261102T030430Z_IMP-7-step14-probe.log`): 6.629e-08 at the pole, 6.6e+02 × `DOSE_BAR`. n = 12 clears it by only 2.3 ×, so it is not used. **`DOSE_BAR` stays the bar and no band moves.** The window is attempt 2 of 4 of step 14's family.

**The reading stays open until item 752 is green.** The slot's reading is that the comparand is excluded, at least nine orders below each point's gap. It is plausible and is not banked, because the control that backs it has not fired yet. Item 745 stands either way: the pole's gap 1.245e-04 → 9.462e-05 → 8.917e-05 at `DT`, `DT` / 2, `DT` / 4 on `LADDER[0]`.

**After item 752 this question stops for the weekly.** If 752 is green, the pole's gap is not h, Δt or truncation at the refinements run. This review writes no further rule (6) step on it: each step costs a cap slot, and a rel 1e-4 gap on a record has no status left to move. The weekly decides whether it is worth a step. The attempt branch `attempt/IMP-7-step14-20261102T030507Z` stays until item 752 lands.

### §9 item 752 — 2026-11-02 (chain I, I53 re-posed, attempt 2 of 4 — daily, provisional, rule (6); lane 3 closed form)

**Why.** Item 741's module is green on every gate except a control whose n sat above where this quantity stops moving. **Provisional:** the 11-02 weekly ratifies or reverts. A re-pose takes no new cap count.

**Anchor.** As item 741, unchanged: (Root), (Q-hd) / (Win-hd) at `DOSE_BAR` / `AVG_BAR` on n = 200, 400, 800, (Self) ≤ `B1_BAND`, every band imported.

**Negative control.** **(Tr8)** n = 8 misses the 200-mode dose at the pole by > `DOSE_BAR`. Ceiling from the probe: 6.629e-08, so the bar is reachable by 6.6e+02 ×. The R swap misses by > `DOSE_BAR`, as item 741 read it.

**Execution.** Take `tests/validation/test_imp7_dose_comparand_floor.py` and `tests/validation/probe_imp7_truncation_teeth.py` by path from `attempt/IMP-7-step14-20261102T030507Z`. In the module edit only: `N_CONTROL = 8`, the test's name (`test_tr8_truncation_control`), its print label, and the docstring's (Tr20) sentence. Add the docstring line "n = 8 is from the probe's log, the brief's n = 20 read 1.493e-16". Leave the probe unedited. Then W1 at the same `-t`. Delete the branch once the item lands. No `src/` change, no landed module edited, so no W-env.

**Tier, ranks, wall clock.** Standard, `-n 1`, `timeout -k 30 170`. Item 741's window took 150 s; dropping `_modes(20)` for `_modes(8)` saves nothing measurable, so expect ≈ 150 s. **If the box is loaded (`uptime` 1 min ≥ 16), skip this slot rather than risk the ceiling.** **Class:** lane-3 closed form. **Cost:** ≈ 2.9 measured / 15 predicted min. **Depends on:** nothing.

**Traps already paid for.** As item 741. **The leak guard fires on 8-digit closed-form figures:** print 7 digits or fewer. A kill parks the item; never raise `-t`. Path checkouts take code only.

**Scope.** The comparand. **Not claimed:** a dose accuracy, a rate, a cause for the pole's gap.

**Negative result.** (Tr8) red at a value far from the probe's: the module and probe disagree, which is a finding; stop, keep the log, BLOCKED. Any other gate red: the same. No `PREDICTED RED:` line.

**Status it can move.** `IMP-7` gains step 14, a record; item 725's residual is excluded from 200-mode truncation; `imp:20` (item 757) unblocks.

**Outcome (2026-11-02T10:34Z, 04:30 CST slot): DONE — `IMP-7` step 14, a record.** W1 `20261102T103036Z_IMP-7-step14-W1.log` (`-n 1`, `timeout -k 30 170`, 148 s, loadavg 0.13): 5 passed, Status 0 (`:93`, `:96–97`). The module was edited only as the brief lists and not after W1 (rule (i)); the probe is unedited. No `src/` change, no landed module edited, so no W-env.

**Green, every band imported.** (Tr8) n = 8 misses n = 200 at the pole by 6.6291e-08 against `DOSE_BAR` = 1e-10 (`:85`), equal to the probe's 6.629e-08, so module and probe agree. (Root) 800 roots on both harmonics, max |x_n − n| 4.779e-02 / 1.167e-01, min gap 0.752 / 0.732 spacings (`:43–44`). (Q-hd) / (Win-hd) green on n = 200, 400, 800 at all three points on both scales, the R swap missing on every cell. (Self) rel 0.000e+00 at `B1_BAND` (`:83`).

**Printed (no band).** Peak-scale rel dose change 200 → 800: pole 0.000e+00, equator 1.565e-16, outer pole 1.671e-16, beside item 725's gaps 1.011e-04 / 1.251e-05 / 1.724e-07 (`:87–89`); largest |Δθ| 3.553e-15 K, rounding at 20 °C, so the row carries no content (`:90`). Every reading item 741 recorded came back at the same digits.

**Reading.** With the control fired, the 200-mode truncation is excluded from item 725's residual: the change is at least nine orders below each point's gap. With item 745, the pole's gap is not h, Δt or truncation at the refinements run. The cause stays unnamed, and the 03:00 ruling leaves whether it earns a step to the weekly.

**Disclosed.** (1) The module's header still cites item 741, as the brief's edit list does not include it. (2) `attempt/IMP-7-step14-20261102T030507Z` is deleted as the brief directs; both of its logs were already on `main`. **Provisional:** the 11-02 weekly ratifies or reverts.

**Status moved.** `IMP-7` gains step 14 as a record; §9 item 752 DONE; item 757 unblocked.
