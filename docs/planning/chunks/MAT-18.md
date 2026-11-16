# MAT-18 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record. Opened by the 2026-10-02 18:00 daily review under §9 rule (5), lane 3, inside Phase 8 subgoal (8a) (the closed-form lane on the sphere); the weekly ratifies or reverts the opening.

### §9 item 190 — 2026-10-02 (lane 3)

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.*

**`MAT-18` step 1: the two-layer sphere in closed form** — a gel core, a container-wall shell and a Robin outer surface, verified by residuals before any FEM.
- **Executor:** implementer. **Tier:** standard, `-n 2`, real build.
- **Cost:** ≈ 8 measured / ≈ 45 predicted min (≈ 1 min window + 15 fixed + ≈ 30 code).
- **Depends on:** nothing.

Every thermal gate so far uses one homogeneous medium. This is the first material interface: θ and k∂ᵣθ continuous across it.

**Fixture.**
- **Core,** r < a = 0.1, gel: k₁ 0.6, ρc₁ 4.18e6, q₁ 2000.
- **Shell,** a < r < b = 0.11: k₂ 0.19, ρc₂ 1.7e6, q₂ 0. These are representative acrylic-wall parameters, not a material claim.
- **Outside:** h 10 on r = b; T∞ = T(0) = 20; t_R = a²ρc₁/k₁ = 69 667 s.

**Anchor (closed form).**
- **Steady, general q₂:**
  - θ(b) = [q₁a³ + q₂(b³ − a³)]/(3hb²);
  - shell: θ = θ(b) + (q₁ − q₂)a³(1/r − 1/b)/(3k₂) + q₂(b² − r²)/(6k₂);
  - core: θ = θ(a) + q₁(a² − r²)/(6k₁).

  This is the spherical-wall resistance (1/r₁ − 1/r₂)/(4πk) of Incropera ch. 3 (radial systems), in series with the convective film.
- **Transient,** with u = rψ and μᵢ = √(s/αᵢ):
  - core: u = sin μ₁r;
  - shell: u₂ = u₂(a)cos μ₂(r − a) + u₂′(a) sin μ₂(r − a)/μ₂, with u₂(a) = sin μ₁a and u₂′(a) = [k₁μ₁cos μ₁a − (k₁ − k₂) sin μ₁a/a]/k₂ (flux continuity);
  - eigencondition: F(s) = k₂u₂′(b) − (k₂/b − h)u₂(b) = 0. F is entire in s, with no poles;
  - solution: θ = θ_ss − Σ cₙψₙe^{−sₙt}, with cₙ = ∫qψₙr²dr/(sₙNₙ) and Nₙ = ∫ρc uₙ²dr.
- **Source:** Carslaw & Jaeger 2e, ch. IX (the sphere; the composite-sphere case). Equation numbers are unpinned (no copy in the repo); the docstring derives every form.

**Code** (additive in `utils/analytical.py`, scipy imported inside the functions):
- `two_layer_sphere_steady_temperature(r, *, a, b, k1, k2, q1, q2=0.0, h, t_inf)`;
- `two_layer_sphere_eigenvalues(*, a, b, k1, k2, rho_c1, rho_c2, h, n)`;
- `two_layer_sphere_transient_temperature(r, t, *, …, t_init, n_terms=400)`;
- the test, `tests/validation/test_mat18_two_layer_sphere_closed_form.py`.

**Anchors (asserted).**
- **(S1) Steady.**
  - θ and k∂ᵣθ continuous at a, ≤ 1e-13 relative.
  - Robin residual at b, ≤ 1e-13.
  - 4πb²hθ(b) = q₁·4πa³/3 to ≤ 1e-13.
  - With k₂ = k₁ and q₂ = q₁, it equals `robin_sphere_steady_temperature(radius=b)` to ≤ 1e-13.
  - With b = a(1 + 1e-9), it equals `…(radius=a)` to ≤ 1e-7 (predicted ≈ 4e-9).
- **(R0) Roots.**
  - |F(sₙ)| ≤ 1e-12 of scale.
  - For n ≤ 60, ψₙ has exactly n − 1 sign changes in (0, b): Sturm–Liouville oscillation, so no root is missed or doubled.
  - Nₙ analytic against Gauss–Legendre split at r = a: ≤ 1e-12 in |ΔN|/max N.
  - Weighted Gram: |⟨ψₘ, ψₙ⟩_{ρc r²}|/√(NₘNₙ) ≤ 1e-10 for m ≠ n ≤ 60.
- **(R1)** cₙ in the Green form against the direct ∫ρcθ_ssψₙr²dr/Nₙ by Gauss–Legendre: ≤ 1e-10 in |Δc|/max|c|.
- **(R2)** Transient interface residuals (θ and k∂ᵣθ) and the Robin residual, from analytic derivatives at Fo 0.05, 0.2 and 1: ≤ 1e-10·q₁a²/k₁.
- **(R3)** The balance d/dt∫ρcθ dV + 4πb²hθ(b) = q₁V_core, analytic per mode at Fo 0.05, 0.2 and 1: ≤ 1e-10 relative.
- **(R4) Limits.**
  - Fo 50 gives θ_ss to ≤ 1e-12.
  - At equal properties with q₂ = q₁, it matches `robin_sphere_transient_temperature(radius=b)` to ≤ 1e-10·θ_ss(0) at Fo 0.05, 0.2 and 1.
  - The root difference against `robin_sphere_eigenvalues(hb/k)` is printed.

**Printed (predicted by hand, never asserted).**
- θ(b) 5.5096 K, θ(a) 8.6994 K, θ(0) 14.2550 K, against `MAT-9`'s 12.2222 K with no wall: the wall adds 14.3 % of the composite value.
- Effective Bi at r = a ≈ 1.28, so λ₁ ≈ 1.73 and τ₁ ≈ 2.3e4 s.
- Root spacing ≈ π/Λ in √s, with Λ = a/√α₁ + (b − a)/√α₂.
- The t = 0 tail at N = 100, 200 and 400.

**Negative controls.** Rule (e): all are predicted, since nothing has been measured on this fixture. Each is asserted only as exceeding its band, with the measured value printed beside the prediction.
- **(c1) Flux jump omitted** (u₂′(a) = u₁′(a)): the (R2) flux residual is predicted O(|1 − k₂/k₁|) ≈ 0.7. Asserted > 1e-3.
- **(c2) Homogeneous-ball roots** `robin_sphere_eigenvalues(hb/k₁)` with the composite eigenfunctions: (R0)'s |F| is predicted O(1). Asserted > 1e-3.
- **(c3) Gram with a uniform weight** ρc₁: off-diagonals predicted 1e-2 to 1e-1. Asserted > 1e-6.

**Window.** Foreground, Bash `timeout: 660000`, the harness call whole and bare. Predicted 50–70 s, against 45–52 s for the analogues.

`scripts/testing/run_and_log.sh MAT-18-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat18_two_layer_sphere_closed_form.py tests/unit/test_analytical_lightweight.py -v -s --tb=short'"`

**Traps.**
- **Composite roots do not fall one per π interval,** so the bracket pattern of `robin_sphere_eigenvalues` does not carry over. Scan √s at a step of π/(32Λ), `brentq` each sign change, then apply the Sturm count.
- **Never divide by F, and never use cot.**
- **Split every radial quadrature at r = a.** One Gauss–Legendre rule across the kink loses digits.
- **ρc is piecewise** in Nₙ and in the Gram weight.
- **Never assert the series near t = 0:** the tail goes as 1/N² (the `MAT-9` lesson).
- **Use analytic derivatives.** An FD check on a summed series needs δ ≳ 1e-2 of the length scale (the 10-01 lesson).
- **Run hygiene.** Print on rank 0. The real build skips the complex-only environment tests. No edit after the last window. Rule (f) re-registration applies to parameters only, never to bounds.
- **Shared file:** item 189 also appends to `analytical.py`. Append at the end, and rebase if 189 lands first.

**Scope.** Numpy only, with perfect thermal contact and a uniform core source. **Not included:**
- the FEM (step 2): it needs a two-region ball generator, because `ball_domain` is homogeneous and `sphere_in_box_domain` is a cube, and it needs a k / ρc that may be a Function in `thermal/`;
- contact resistance;
- the cylinder;
- perfusion.

**Status it can move.** §7 `MAT-18` ⬜ → 🟡 on the row's step-1 Done-when, which this item matches (rule (j)). Step 2 then becomes draftable.

**Negative result.**
- **(S1), (R0)'s |F|, (R1) or (R2) red:** a code bug. Fix it in-slot and run W1b, disclosed.
- **The Sturm count red:** a root was missed. Print the scan around it and stop.
- **(R4) red with (R0)–(R3) green:** the composite form and the homogeneous form disagree. Print both at three points and stop.
- **On either stop:** park on `attempt/MAT-18-step1-<UTC>`, mark BLOCKED in the same commit (rule (d)); the row stays ⬜.
- There is no PREDICTED RED line.

**Outcome — BLOCKED 2026-10-03 (21:00 CDT slot, executor `implementer`).** Red on (R0)'s |F| alone after the one permitted fix window; parked on `attempt/MAT-18-step1-20261003T021930Z`, row stays ⬜.
- **Windows (on the branch):** `20261003T021446Z_MAT-18-step1-W1.log` (`-n 2`, 63 s, 1 failed / 19 passed / 4 skipped; R0 7.770e-11 with the post-cancellation scale), `20261003T021651Z_MAT-18-step1-W1b.log` (63 s, same counts; R0 **1.823e-11 at n = 364** vs 1e-12), `20261003T021825Z_MAT-18-step1-probe.log` (diagnostic, 1 s; script `scratch_mat18_probe.py` on the branch).
- **W1b fix (disclosed, band unchanged):** F scaled by the size of its parts before cancellation (W1's scale kept as `want="scale_post"`, printed 6.829e-11); brentq root polished to the |F|-minimising double within ±16 ulp.
- **Why it stays red:** at n = 364 (x = 3.886) one ulp of x moves F/scale by 2.74e-11, so the best float64 root sits ≈ 1.4e-11 from zero; 1e-12 is unreachable in double precision above n ≈ 300. For n ≤ 60 the maximum is 2.09e-13; the median over 400 roots is 4.7e-14. Note the §7 row's step-1 Done-when names no |F| band — R0's 1e-12 is the brief's.
- **Green (W1b):** S1 interface / Robin / balance ≤ 2.13e-16, equal properties 2.53e-16, thin-wall limit 1.78e-9 (band 1e-7); Sturm count exact n ≤ 60; Nₙ 1.23e-15; Gram 8.96e-15; R1 1.26e-12; R2 ≤ 7.1e-17; R3 ≤ 2.1e-16; R4 2.49e-16 and equal-property transient ≤ 5.1e-16; equal-property roots vs `robin_sphere_eigenvalues` 4.5e-13.
- **Controls:** c1 0.153, c2 max 0.999 / median 0.918, c3 0.2455 — all three fire.
- **Records:** θ(b) / θ(a) / θ(0) = 5.5096 / 8.6994 / 14.2550 K; the wall adds 14.26 % over `MAT-9` (predicted 14.3 %); λ₁ = 1.7297, Bi_eff 1.2772.
- **Proposed re-pose (executor's, for the review):** 1e-12 on n ≤ 60, and for all 400 roots a sign change of F within ±2 ulp with the stored root the |F|-minimiser among its neighbouring doubles.

**Ruling, 2026-10-03 03:00 daily review (R0 re-posed; the 10-03 weekly ratifies or reverts, together with the chunk itself).** R0's |F(sₙ)| ≤ 1e-12·scale was the item-190 brief's own bound; the §7 step-1 Done-when names no |F| band. The W1b probe shows it is unreachable in float64 above n ≈ 300: at n = 364 one ulp of x moves F/scale by 2.74e-11 (`20261003T021825Z_MAT-18-step1-probe.log`, on the branch). A bound below one ulp's worth of F tests the arithmetic, not the root. **Re-posed R0 (attempt 2 of 4, same target, same code):** (R0a) |F(sₙ)| ≤ 1e-12·scale for n ≤ 60, the range the Done-when's Sturm and Gram criteria already use (measured max 2.09e-13); (R0b) for all 400 roots, F changes sign between the stored root's neighbouring doubles within ±2 ulp (`np.nextafter`), and the stored root is the |F|-minimiser among the doubles within ±2 ulp. (R0b) certifies each root to float64 resolution, which is the most any double-precision root can carry. Every other anchor, band and control of item 190 is unchanged.

### §9 item 199 — 2026-10-03 (lane 3, attempt 2 of 4)

*Written by the 2026-10-03 03:00 daily review. Licence: rule (5) lane 3, re-pose of item 190 under the ruling above. Item 190 is superseded by this one.*

**Outcome — DONE 2026-10-03 (09:00 CDT slot, executor `implementer`).** Pre-check empty; `analytical.py` and the test landed from the branch by path, the three branch logs landed and their rows appended by hand; the test's R0 alone edited to (R0a) + (R0b) with the +8-ulp control. W2 `20261003T140149Z_MAT-18-step1-W2.log` (`-n 2`, 64 s, standard): 20 passed / 4 skipped. (R0a) 2.093e-13 at n = 30 (band 1e-12); (R0b) 0 / 400 without a ±2-ulp sign change and 0 / 400 non-minimisers; +8-ulp control 400 / 400 fail the minimiser check (pred ≥ 1, asserted > 0); all-400 max \|F\|/scale 1.823e-11 at n = 364 (record). Every other anchor and control reproduces W1b to the printed digits (S1 ≤ 2.13e-16, thin wall 1.78e-9, Sturm exact, Nₙ 1.23e-15, Gram 8.96e-15, R1 1.26e-12, R2 ≤ 7.1e-17, R3 ≤ 2.1e-16, R4 ≤ 5.1e-16; c1 0.153, c2 0.999, c3 0.2455). Row ⬜ → 🟡; branch deleted.

- **Executor:** implementer. **Tier:** standard, `-n 2`. **Class:** lane-3 closed form. **Cost:** ≈ 8 measured / ≈ 20 predicted min (W1b 63 s). **Depends on:** nothing.

**Do.**
1. **Pre-check:** `git log --oneline fcdade64..HEAD -- src/fem_em_solver/utils/analytical.py` must be empty (it was at `341fe583`). If it is empty, take `src/fem_em_solver/utils/analytical.py` and `tests/validation/test_mat18_two_layer_sphere_closed_form.py` from `attempt/MAT-18-step1-20261003T021930Z` by path. If it is **not** empty (another item appended to `analytical.py` first), apply the branch's diff instead: `git diff fcdade64 attempt/MAT-18-step1-20261003T021930Z -- src/fem_em_solver/utils/analytical.py` applied with `git apply`; never path-checkout over another item's append.
2. Land the branch's three logs by path. Append their three `test-results.md` rows **by hand** (never path-checkout `test-results.md`; `git diff --stat` must show insertions only on it). Do not land `scratch_mat18_probe.py`.
3. Edit the test module's R0 only, to (R0a) + (R0b) as ruled above. Print, for all 400 roots, max |F|/scale and the n at which it occurs (a record), and the count of roots whose stored double is not the ±2-ulp |F|-minimiser (must be 0). Change nothing else in the module.
4. Run W2. Delete the branch after the landing commit (`git branch -D`).

**Anchor.** Item 190's anchors (S1), (R0a), (R0b), (R1)–(R4), its Sturm count for n ≤ 60, Nₙ and the weighted Gram, all at item 190's bands (above, unchanged).

**Negative controls.** Item 190's c1, c2 and c3, unchanged. The measured values from W1b (c1 0.153, c2 max 0.999, c3 0.2455) are **asserted** labels now for their exceedance bands (same comparison, same fixture: W1b). Added control for (R0b), **predicted**: the stored roots perturbed by +8 ulp must fail (R0b)'s minimiser check for at least one root (printed count beside the prediction "≥ 1"); asserted as > 0.

**Window.** Foreground, Bash `timeout: 660000`, the harness call whole and bare. Predicted 55–75 s.

`scripts/testing/run_and_log.sh MAT-18-step1-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat18_two_layer_sphere_closed_form.py tests/unit/test_analytical_lightweight.py -v -s --tb=short'"`

**Traps.** Item 190's list, plus: `np.nextafter` on a numpy float64 scalar, not a Python float cast through float32; evaluate F at the neighbours with the same scaled form W1b used (`want="scale"`), never the post-cancellation scale; no edit after W2 (rule (i)); landing a parked branch by path checkout of an append-only record deletes rows (the 2026-09-10 lesson).

**Scope.** Numpy only, step 1. Not covered: the FEM (step 2, item 201), contact resistance, the cylinder, perfusion.

**Negative result.** (R0b) red ⇒ a root is not float64-correct: print the ±4-ulp neighbourhood of the worst root and stop; park and mark BLOCKED (rule (d)). Any other red ⇒ the landed code is not the branch's: stop and report. There is no PREDICTED RED line.

**Status it can move.** §7 `MAT-18` ⬜ → 🟡 on the row's step-1 Done-when, which this item matches (rule (j)); §9 item 190 → superseded; item 201 (step 2) becomes runnable.

### §9 item 201 — 2026-10-03 (lane 3)

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.*

**`MAT-18` step 2(a): the conforming two-region ball and piecewise k in the steady solver, gated against step 1's steady closed form.**
- **Executor:** implementer. **Tier:** standard, `-n 2`; W1, W2a and W2b real build, W3a and W3b complex. **Class:** lane-3 closed form (measured median 8 min).
- **Cost:** ≈ 8 measured / ≈ 50 predicted min (15 fixed + ≈ 25 code + five ≈ 2-min windows).
- **Depends on:** item 199. Step 1's `two_layer_sphere_steady_temperature` must be on `main`; **if item 199 did not land, skip this item.**
- `utils/analytical.py` is only imported here, not edited.

**Why split.** The row's step 2 bundles a new mesher, piecewise k and ρc in both solvers, a steady ladder, a CN run and a six-module re-gate. That is ≈ 70 predicted minutes, more than one slot. So:
- **2(a), this item:** the mesher, piecewise k in `solve_steady_conduction`, the steady anchors and the uniform-k control, plus the re-gate of `steady.py`'s users;
- **2(b), not briefed here:** see the end of this brief.

The review amends the row's step-2 Done-when into 2(a) / 2(b) in the commit that queues this item (rule (j)).

**Code (additive; rule (c), disclosed).**
- **`io/mesh.py`:** `MeshGenerator.two_layer_ball_domain(inner_radius=0.1, outer_radius=0.11, resolution=0.01, comm=MPI.COMM_WORLD, rank=0)`, modelled on `ball_domain` (`io/mesh.py:2649–2701`).
  - OCC spheres of radius a and b, then `occ.fragment`, following the `sphere_in_box_domain` precedent (`:2320`). Identify the two volumes by `occ.getMass(3, ·)` (the core is 4πa³/3), never by tag order.
  - Cell tag 1 is the core and cell tag 2 the shell. Facet tag 1 is the r = b surface only.
  - Uniform `Mesh.MeshSizeMin` = `Mesh.MeshSizeMax` = `resolution`, with point, boundary and curvature sizing off, exactly as `ball_domain` does.
  - **No dim-2 physical group on the r = a interface.** Interior facet groups hang `model_to_mesh` at `-n 2` (known-issues 9, `io/mesh.py:124–130`). If a later step needs the interface, use `_interface_facet_tags(msh, cell_tags, {2: (1, 2)}, existing=facet_tags)` (`:108`).
- **`thermal/steady.py`:** `k` in `solve_steady_conduction` and `_solve` may be a UFL expression or a `fem.Function` (DG0).
  - A scalar still goes through `fem.Constant` exactly as now (`steady.py:96`), so the scalar path's forms are byte-identical.
  - The positivity check (`:174`) applies to scalars. A Function is checked by an allreduced minimum of its array > 0.
  - `solve_steady_conduction_picard` is untouched.
- **New test** `tests/validation/test_mat18_two_layer_sphere_fem.py`.
  - Constants come from step 1's module (`A_, B_, K1, Q1, K2, Q2, H, T_INF`; `test_mat18_two_layer_sphere_closed_form.py:37–41` on the parked branch, landed by item 199).
  - DG0 k and q are set from `cell_tags.find(1)` / `find(2)`, then `x.scatter_forward()`: ghost cells are not carried by `cell_tags` (`io/mesh.py:132–136`).
  - Comparand: `two_layer_sphere_steady_temperature`, interpolated into P2. It is cheap: there is no series.

**Fixture.**
- a = 0.1 and b = 0.11.
- Gel core: k₁ 0.6, q₁ 2000. Shell: k₂ 0.19, q₂ 0.
- h 10 on r = b, T∞ 20. P1.
- Rungs: `resolution` ∈ {0.01, 0.007071, 0.005, 0.003536}, ratio √2.
- Predicted cells ≈ 26k / 73k / 0.20M / 0.57M: `MAT-9`'s R = 0.1 ball read 19 463 / 151 548 cells at 0.01 / 0.005 (`20260929T034247Z_MAT-9-step1.log:87–88`), scaled by (b/R)³ = 1.331.
- **Cell guard: 800 000 at the finest rung, checked before its solve.** Over it, stop.

**Anchors (asserted).** The bands are `MAT-9` step 1's, imported by value: they are literals at `test_mat9_steady_conduction_sphere.py:139, 154, 156`.
- **(G) Conformity and census,** on every rung, rank-reduced, owned cells only (`size_local`):
  - owned tag-1 cells plus owned tag-2 cells equal the global cell count;
  - every vertex of a core cell has \|x\| ≤ a(1 + 1e-9);
  - every vertex of a shell cell has a(1 − 1e-9) ≤ \|x\| ≤ b(1 + 1e-9), so no tet straddles r = a;
  - tag-1 facets are exterior, and their area is within 1 % of 4πb² (predicted −1e-3 to −4e-3).
- **(S1) Rate:** the whole-ball relative L2 error on θ, fitted as the LSQ slope of log e against log h_eff with h_eff = (V_h/N)^{1/3} (`MAT-15`'s form), is ≥ 1.7. Predicted ≈ 2.0.
- **(S2) Interior:** the interior relative L2 error ≤ 1 % at 0.003536. Interior means cells with midpoint \|x\|² < 0.81 b², `MAT-9`'s mask scaled to the outer radius. Predicted 4e-4 to 1.5e-3.
- **(S3) Balance,** every rung: |∮hθ dS − q₁V_core,h| / (q₁V_core,h) ≤ 1e-9, with V_core,h the allreduced tag-1 volume. Predicted ≤ 2e-14.
- **(D) Coefficient path:** DG0 k ≡ K1 everywhere against scalar k = K1, on the 0.01 rung with the same piecewise q. Whole-ball relative difference ≤ 1e-13; predicted ≤ 1e-15.
- **(E) Equal properties:** DG0 k ≡ K1 and q ≡ Q1 in both regions, against `robin_sphere_steady_temperature(radius=b, k=K1, q=Q1, h=H)` at 0.005. Interior ≤ 1 %; predicted ≈ 8e-4 (`…MAT-9-step1.log:88`).

**Printed.**
- Per rung: cells; h_eff; build and solve seconds; the whole, interior and shell-only errors; V_core,h/(4πa³/3) − 1; the outer-area deficit.
- The pairwise rates.
- θ_h at r = 0, a and b, beside step 1's records 14.2550 / 8.6994 / 5.5096 K (this file, item 190's outcome).

**Negative controls.** Rule (e): both are predicted. Each is asserted only as missing its band.
- **(c1) Uniform k:** scalar k = K1 everywhere, with the piecewise q, against the two-layer comparand at 0.005. The core shifts uniformly by θ(a) − θ_uniform(a) = 8.6994 − 6.5197 = 2.180 K, so the interior miss is predicted 0.197. Asserted > 1e-2.
- **(c2) Source in both regions:** q₁ applied to both regions, with the balance taken against q₁V_core,h. Predicted V_ball,h/V_core,h − 1 ≈ 0.331. Asserted > 1e-3.

**Windows.** Foreground, Bash `timeout: 660000`. Each harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.
- **W1** (real build; predicted 90–150 s). The environment tests take ≈ 30 s. `MAT-9`'s 0.005 m rung was built and solved in 6.4 s at 152k cells (`…MAT-9-step1.log:88`), so the four rungs and the controls take ≈ 60–100 s.

  `scripts/testing/run_and_log.sh MAT-18-step2a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat18_two_layer_sphere_fem.py -v -s --tb=short'"`
- **W2a** (real build; these three modules measured 123 s together, `20260930T020340Z_MAT-5-step1.log:300`).

  `scripts/testing/run_and_log.sh MAT-18-step2a-W2a-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat9_steady_conduction_sphere.py tests/validation/test_mat9_transient_conduction_sphere.py tests/validation/test_mat5_linear_sigma_sphere.py -v -s --tb=short'"`
- **W2b** (real build; predicted 110–140 s). `MAT-11` measured 62 s (`20260930T095245Z_MAT-11-W1.log:291`) and `MAT-13` 75 s (`20260930T172102Z_MAT-13-W1.log:245`), each including the environment, which they share here.

  `scripts/testing/run_and_log.sh MAT-18-step2a-W2b-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat11_perfusion_sphere.py tests/validation/test_mat13_transient_perfusion_sphere.py -v -s --tb=short'"`
- **W3a** (complex build; predicted 100–120 s). `MAT-10` measured 70 s (`20260930T021802Z_MAT-10-step1-W1.log:302`) and `MAT-12` 77 s (`20260930T110230Z_MAT-12-W1.log:333`), each including ≈ 44 s of environment and collection, which they share here.

  `scripts/testing/run_and_log.sh MAT-18-step2a-W3a-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat10_sar_to_heat_sphere.py tests/validation/test_mat12_transient_sar_heating_sphere.py -s -v --tb=short'"`
- **W3b** (complex build; `MAT-14` measured 79 s, `20261001T004928Z_MAT-14-W1.log:309`). It is kept separate so that no complex window crowds the 180 s ceiling.

  `scripts/testing/run_and_log.sh MAT-18-step2a-W3b-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat14_sar_pennes_transient_sphere.py -s -v --tb=short'"`
- **Not re-run, disclosed:** `MAT-5` step 2 (`test_mat5_em_coupled_sphere.py`, 178 s, complex) and `MAT-15`'s two modules also call `solve_steady_conduction`. They are heavy and outside the row's list. The evidence for them is the byte-identical scalar path.

**Traps.**
- **The interface physical group** hangs `model_to_mesh` at `-n 2` (above).
- **Ghost cells and the DG0 coefficients:** set the values from the tags, then `scatter_forward` (above).
- **Rank-local data:** `cell_tags.find` and `assemble_scalar` are rank-local, so reduce every count and integral.
- **Geometry:** read `geometry.dofmap`'s shape at runtime (the `IMP-2` lesson).
- **The outer facet group** is the surface whose `occ.getMass(2, s)` ≈ 4πb², or the shell's boundary minus the shared interface. Never the core's boundary.
- **At 0.01 the shell is ≈ 1 cell thick.** That is expected. The rung stays: never drop a rung in-slot (rules (f), (h)).
- **The comparand chooses its layer by exact r.** Points between the faceted interface and r = a get the other layer's formula. That is consistent, because θ and the flux are continuous there; do not special-case it.
- **Masks are built in numpy** (a DG0 indicator from cell midpoints), never with UFL comparisons.
- **Run hygiene:** print on rank 0; run with `-s`; never `-k` inside the quoted command. Rule (i): no edit after a file's last window.
- **After a kill:** `pgrep -c python3` in the container first, then the any-`.c`-without-`.so` sweep (daily-review.md rubric item 4).

**Scope.** Steady state only, P1, perfect contact, a uniform core source. **Not included:**
- the transient (2(b)) and a piecewise ρc;
- contact resistance;
- a real container's properties;
- the cylinder (`MAT-15`);
- perfusion;
- degree 2.

**Status it can move.** No glyph. After item 199 the row reads 🟡, and ✅ needs 2(b). The §7 `MAT-18` row registers step 2(a) green as a record, on the Done-when as the queueing review amends it, and 2(b) becomes draftable. If item 199 left the row ⬜, this item is skipped.

**Negative result.**
- **(G) red:** the mesher is non-conforming or mis-tagged. One fix in-slot with a W1b; if it is still red, park and mark BLOCKED.
- **(D) or (S3) red:** a code bug. Fix it in-slot, run W1b, and disclose.
- **(S1) below 1.7 with (S2) green:** print the pairwise rates. If only the 0.01 pair is low, do not drop that rung (rules (f), (h)). Record it, park, BLOCKED; the review rules.
- **(S2) above 1 %, or (E) red:** report the reading beside the volume and area deficits, then stop. Park, BLOCKED.
- **Any W2a / W2b / W3a / W3b red:** the scalar path is not byte-identical, so stop. Park, BLOCKED.
- **Over the cell guard:** stop before the solve, report the counts, BLOCKED.
- On every stop, park on `attempt/MAT-18-step2a-<UTC>` and mark the item BLOCKED in the same commit (rule (d)).
- There is no PREDICTED RED line.

**Outcome — DONE 2026-10-03 (12:00 CDT slot, third item, implementer).** All five windows were green on the first run (`-n 2`, standard): no kill, no in-slot fix, no W1b.

Windows:
- W1 `20261003T172548Z_MAT-18-step2a-W1.log`: 13 passed / 4 skipped, 121 s.
- W2a `…172759Z…-W2a-regate.log`: 20 / 4, 114 s.
- W2b `…173000Z…-W2b-regate.log`: 24 / 4, 80 s.
- W3a `…173127Z…-W3a-regate.log` (complex): 14 passed, 93 s.
- W3b `…173307Z…-W3b-regate.log` (complex): 14 passed, 83 s.

| res | cells (core / shell) | h_eff | whole | interior | shell | V_core,h/V − 1 | area deficit | balance |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 29 788 (19 470 / 10 318) | 5.714e-3 | 4.131e-3 | 4.219e-3 | 3.422e-3 | −3.648e-3 | −1.662e-3 | 2.34e-15 |
| 0.007071 | 74 830 (54 351 / 20 479) | 4.206e-3 | 2.931e-3 | 2.541e-3 | 4.920e-3 | −1.820e-3 | −8.35e-4 | 1.91e-15 |
| 0.005 | 215 542 (151 410 / 64 132) | 2.957e-3 | 1.121e-3 | 1.121e-3 | 1.121e-3 | −9.15e-4 | −4.19e-4 | 1.02e-14 |
| 0.003536 | 584 228 (423 509 / 160 719) | 2.121e-3 | 7.603e-4 | 6.791e-4 | 1.189e-3 | −4.61e-4 | −2.11e-4 | 5.52e-15 |

Anchors:
- **(G)** green on every rung. Vertex radii are within 2.2e-16 of a and b; every tag-1 facet is exterior.
- The 0.01 area deficit is −1.66e-3, inside the 1 % band; the predicted −1e-3 to −4e-3 held from 0.005 down. The cell guard read 584 228, under 800 000.
- **(S1)** LSQ 1.822, against a predicted ≈ 2.0. The pairwise rates are uneven, **1.120 / 2.728 / 1.168**. This is recorded for the review; the assertion is on the LSQ rate, per the brief.
- **(S2)** 6.791e-4. **(S3)** ≤ 1.02e-14. **(D)** 6.165e-17. **(E)** interior 7.250e-4 at 0.005.
- θ_h at 0.003536: 14.2498 / 8.6897 K at r = 0 / a, and 5.5392 K at r = 0.999b. The step-1 records are 14.2550 / 8.6994 / 5.5096 K (the last at b).

Controls: (c1) 0.1986 (predicted 0.197); (c2) 0.3312 (predicted 0.331).

Disclosed:
- `io/mesh.py` gains `two_layer_ball_domain`. The core is found by `getMass(3,·)`, and the outer facet group by `getMass(2,·)` ≈ 4πb². There is no interface group.
- `thermal/steady.py`: `k` goes through `_as_source`, so a scalar still becomes the same `fem.Constant`. Positivity is checked as before for scalars, by an allreduced minimum for a `fem.Function`, and only h > 0 is checked for other UFL expressions.
- `steady.py`'s docstring and type hint still say homogeneous `k: float`. They were left unedited so that rule (i) holds.
- The r = b point is evaluated at 0.999b through `evaluate_vector_field_parallel`, because a point on the faceted surface can fall outside the mesh.
- `MAT-5` step 2 and `MAT-15` were not re-run; their evidence is the byte-identical scalar path.

Row stays 🟡 (a record); 2(b), item 205, is runnable now that items 200 and 201 have both landed.

**What 2(b) would be (not briefed).**
- Piecewise k and ρc in `solve_transient_conduction`, with the scalar path byte-identical. It is rebased on item 200's `source_schedule` keyword if that has landed, since both edit `transient.py`.
- CN with the Rannacher start, N 40, to Fo 0.2 (t = 0.2·a²ρc₁/k₁) at 0.005, against `two_layer_sphere_transient_temperature`: interior ≤ 1 %.
- Every run's per-step balance ≤ 1e-9, normalised by q₁V_core,h.
- Backward Euler 20 × 5 t_R equals 2(a)'s steady solve to ≤ 1e-8 (`MAT-9`'s (iv)).
- A uniform-ρc control, printed.
- `MAT-9`…`MAT-14` re-run green ⇒ ✅.

### §9 item 205 — 2026-10-03 (lane 3)

*Drafted by a general-purpose agent for a 2026-10-03 daily review; the review reviewed it.*

**`MAT-18` step 2(b): piecewise k and ρc in the transient solver, gated against step 1's composite series at Fo 0.2 on 2(a)'s two-region ball.**
- **Executor:** implementer. **Tier:** standard, `-n 2`; W1, W2a, W2b and W2c real build, W3a and W3b complex. No width check: the `MAT-9`…`MAT-17` thermal idiom has none, and `-n 2` is the width where a missing reduction shows. **Class:** lane-3 closed form (measured median 8 min).
- **Cost:** ≈ 8 measured / ≈ 50 predicted min (15 fixed + ≈ 20 code + six ≈ 2-min windows).
- **Depends on:** item 201 (2(a): `MeshGenerator.two_layer_ball_domain` and DG0 k in `solve_steady_conduction`) **and** item 200 (`MAT-17` step 2: `source_schedule` in `thermal/transient.py`, the same file). **If either did not land, skip this item.** At drafting (2026-10-03, after `657d4721`) neither had: `transient.py` has no `source_schedule`, and `tests/validation/test_mat18_two_layer_sphere_fem.py` does not exist.
- `utils/analytical.py` and `io/mesh.py` are only imported here, not edited.

**Target — the §7 row's 2(b) Done-when, verbatim:** "*2(b)* is the transient: CN at Fo 0.2 ≤ 1 % with piecewise ρc, and `MAT-9`…`MAT-14` re-run green. 2(a) + 2(b) ⇒ ✅." The row names no rung and no Fo convention. This item pins them from item 201's "What 2(b) would be" notes: the 0.005 rung, with Fo = t/t_R and t_R = a²ρc₁/k₁, item 190's core time.

**Code (additive; rule (c), disclosed).**
- **`thermal/transient.py`:** `k` and `rho_c` in `solve_transient_conduction` and `_solve_transient` may be a UFL expression or a `fem.Function` (DG0). Rebase on item 200's `source_schedule` edit; the two compose and do not interact.
  - A scalar still goes through `fem.Constant` exactly as now (`transient.py:158–159` before item 200), so the scalar path's forms are byte-identical and hit the same JIT cache entries. Reuse 2(a)'s scalar-or-Function helper from `steady.py` if it exposes one (as `_as_source` is reused at `:55`, `:160`); otherwise mirror it.
  - The positivity check (`:143`) applies to scalars. A Function is checked by an allreduced minimum of `np.real(f.x.array)` > 0.
  - `rc_c` enters three forms: `a_form` (`:180`), the `L` storage term (`:182`) and `storage_form` (`:190`). All three take the same object. A uniform ρc in `storage_form` alone would leave the solve right and the balance wrong by control (c2)'s ≈ 0.07.
  - The balance stays `|stor + conv + perf − src| / src` with `src = ∫q dx`. Under 2(a)'s piecewise q (q₂ = 0) that is already q₁V_core,h, as the row's normalisation requires.
- **New test** `tests/validation/test_mat18_two_layer_sphere_transient_fem.py`. It is new so that 2(a)'s module stays as committed (rule (i)).
  - Constants come from step 1's module (`A_, B_, K1, RC1, Q1, K2, RC2, Q2, H, T_INF, T_R`; `test_mat18_two_layer_sphere_closed_form.py:37–41`). Helpers from 2(a)'s module (mesh cache, DG0 builders from `cell_tags`, interior indicator) are imported read-only where it exposes them; anything else is re-built here. `_say`, `_allsum` and `_rel_diff`'s form come from `test_mat9_transient_conduction_sphere.py`, read-only.
  - DG0 k, ρc and q are set from `cell_tags.find(1)` / `find(2)`, then `x.scatter_forward()` (2(a)'s ghost-cell rule).
  - Comparand: `two_layer_sphere_transient_temperature(r, t, a=A_, b=B_, k1=K1, k2=K2, rho_c1=RC1, rho_c2=RC2, q1=Q1, q2=Q2, h=H, t_inf=T_INF, t_init=T_INF)` (`analytical.py:2189`, default `n_terms=400`). It is evaluated in numpy at the P2 dof coordinates in chunks of 20 000 points, as `MAT-9`'s `_exact_theta` does (`test_mat9_transient_conduction_sphere.py:93–106`).

**Fixture.** 2(a)'s ball and coefficients (a 0.1, b 0.11; core k₁ 0.6, q₁ 2000; shell k₂ 0.19, q₂ 0; h 10, T∞ 20), plus ρc₁ 4.18e6 and ρc₂ 1.7e6. T(0) = T∞, with the source on for t > 0. P1.
- t* = 0.2·T_R = 0.2·a²ρc₁/k₁ = 13 933 s. That is numerically `MAT-9`'s `T_END`, because a = R and the core constants are `MAT-9`'s. So N 40 gives Δt 348.33 s, `MAT-9`'s (i) step.
- Rungs: 0.005 for (T1), (E) and (c1); 0.01 for (D), (T3) and the printed CN ladder.
- **Cell guard: 400 000 at 0.005, checked before the first solve on it.** 2(a)'s W1 log has the measured count; item 201 predicted ≈ 0.20M. Over it, stop.

**Anchors (asserted).** The bands are `MAT-9` step 2's, imported by value: they are literals at `test_mat9_transient_conduction_sphere.py:274–275, 301, 321–322`, with no module constant to import. (D)'s band is 2(a)'s (D) literal.
- **(T1) Interior, the row's criterion:** CN with the Rannacher start (`theta=0.5`, `startup_half_steps=2`), N 40 to t*, on the 0.005 rung with DG0 k and DG0 ρc. Interior relative L2 on θ against the comparand ≤ 1 % (`:275`). Interior means 2(a)'s mask, cell midpoints with \|x\|² < 0.81 b²; that is r < 0.099, inside the core. Predicted 5e-4 to 2e-3: `MAT-9`'s CN N 40 read 4.6327e-4 on its 0.005 ball (`20260929T170037Z_MAT-9-step2b.log:114`), and 2(a)'s steady interior is predicted 4e-4 to 1.5e-3.
- **(T2) Balance:** every run's `max_balance_residual` ≤ 1e-9 (`:274`, `:301`), normalised by q₁V_core,h. Predicted ≤ 5e-14 (`MAT-9` read ≤ 3.086e-14, `…step2b.log:138`).
- **(T3) Long BE = steady:** BE, 20 steps of 5·T_R, piecewise k, ρc and q, on the 0.01 rung, against `solve_steady_conduction` with 2(a)'s DG0 k and the same q on the same mesh. Whole-ball relative difference ≤ 1e-8 (`:322`). Predicted ≈ 1e-14 (`MAT-9`'s (iv) read 4.301e-15, `…step2b.log:141`; the slowest decay is item 190's τ₁ ≈ 2.3e4 s, so the BE factor (1 + 5T_R/τ₁)⁻²⁰ ≈ 16⁻²⁰).
- **(D) Coefficient path:** DG0 k ≡ K1 and DG0 ρc ≡ RC1 against scalar k = K1, rho_c = RC1, with the same piecewise q, on the 0.01 rung, BE N 10 to t*. Whole-ball relative difference ≤ 1e-13 (2(a)'s (D) band); predicted ≤ 1e-15.
- **(E) Equal properties:** DG0 k ≡ K1, ρc ≡ RC1 and q ≡ Q1 in both regions against `robin_sphere_transient_temperature(r, t*, radius=B_, k=K1, rho_c=RC1, q=Q1, h=H, t_inf=T_INF, t_init=T_INF)`, CN N 40 + Rannacher at 0.005. Interior ≤ 1 % (`:275`); predicted ≈ 5e-4 to 2e-3. This separates a mesh or comparand fault from a ρc fault if (T1) goes red.
- **(R) Re-gate, the row's second criterion:** W2a–W3b green. W2a's `MAT-9` transient lines reproduce `…step2b.log:114, 124, 130` digit for digit (interior 4.6327e-04, BE LSQ 0.984, CN 2.001). W2b's `MAT-13` (P_w) reads 1.6426e-03 (`20260930T172102Z_MAT-13-W1.log:97`). W3a's `MAT-12` (P_t) reads 1.6845e-3 (`20260930T110230Z_MAT-12-W1.log:161`). W3b's `MAT-14` (P_w′) reads 1.6089e-03 (`20261001T004928Z_MAT-14-W1.log:162`). W2c's `MAT-17` (i) digits equal item 200's own W1 log. Journal the comparison. A changed digit on the scalar path is a stop.

**Printed (predicted, never asserted).**
- Per run: cells; V_core,h/(4πa³/3) − 1 and V_shell,h/(4π(b³ − a³)/3) − 1; build, solve and comparand seconds.
- (T1)'s shell-only error. θ_h beside the comparand at r = 0, a and b at t*. A two-node lumped estimate puts the core mean at ≈ 4.8 K and the shell at ≈ 3.0 K. It is crude, so it is printed only.
- The CN + Rannacher self-difference on the 0.01 rung, N 20 / 40 / 80, and its pairwise rate. Predicted ≈ 2.0 (`MAT-9` read 2.001, `…step2b.log:130`).
- CN without the Rannacher start at N 40 on 0.005.

**Negative controls.** Rule (e): both are **predicted**, since this comparison has never been measured on this fixture. Each is asserted only as missing its anchor's band, with the measured value printed beside the prediction.
- **(c1) Uniform ρc:** DG0 k piecewise, scalar `rho_c = RC1` everywhere (the shell gets 4.18e6 instead of 1.7e6, so its capacity is 5 796 instead of 2 357 J/K), against the two-layer comparand. Same run as (T1) otherwise.
  - Predicted interior miss ≈ 0.05 (range 0.02–0.10). This comes from a two-node lumped model (core C₁ 17 509 J/K; shell node at the middle of the wall resistance; R_core 0.265, R_wall 0.381, R_film 0.658 K/W; its slow time constant 2.42e4 s reproduces item 190's τ₁ ≈ 2.3e4 s). The model takes the core mean at t* from 4.80 K to 4.56 K.
  - **Ceiling:** an infinitely capacitive wall in the same model holds the core mean to 3.15 K, a 0.34 miss. A bar is reachable below ≈ 0.3; the asserted 1e-2 sits ≈ 5× under the central prediction.
  - Asserted > 1e-2.
- **(c2) Uniform ρc in the balance:** on the (T1) run, take snapshots 39 and 40. Compute |ρc₁∫(T₄₀ − T₃₉)/Δt + ∮h(½T₄₀ + ½T₃₉ − T∞) − ∫q| / ∫q, i.e. the storage term weighted by ρc₁ everywhere instead of the DG0 ρc. This is the form of `MAT-9`'s control (b) (`test_mat9_transient_conduction_sphere.py:285–292`).
  - Predicted ≈ 0.07 (range 0.03–0.15): the shell's mis-weighted storage rate, (ρc₁ − ρc₂)·V_shell·dθ_shell/dt ≈ 3 439 J/K × 1.70e-4 K/s ≈ 0.58 W, over P = q₁V_core = 8.38 W.
  - Asserted > 1e-3.

**Windows.** Foreground, Bash `timeout: 660000`. Each harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run W1 first; the re-gate windows follow on the final code.
- **W1** (real build; predicted 80–130 s). The environment tests take ≈ 30 s. At 0.005 there are four CN N 40 runs ((T1), (c1), (E), no-startup) at ≈ 5 s each: `MAT-9` took 3.7 s at 152k cells (`…step2b.log:114`), scaled to ≈ 0.2M. Two comparand evaluations take ≈ 5–10 s; each chunk's call recomputes the 400 roots (≈ 0.5 s predicted, printed). Add 2(a)'s ≈ 10–15 s build at 0.005 and the 0.01-rung runs at ≈ 10 s.

  `scripts/testing/run_and_log.sh MAT-18-step2b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat18_two_layer_sphere_transient_fem.py -v -s --tb=short'"`
- **W2a** (real build). These two modules measured 115 s (`20260929T170037Z_MAT-9-step2b.log`, test-results row of 2026-09-29 17:02) and 105 s as `MAT-12`'s rule-(c) re-gate (`20260930T110355Z_MAT-12-W2-ruleC.log`).

  `scripts/testing/run_and_log.sh MAT-18-step2b-W2a-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat9_steady_conduction_sphere.py tests/validation/test_mat9_transient_conduction_sphere.py -v -s --tb=short'"`
- **W2b** (real build; predicted 110–140 s). `MAT-11` measured 62 s (`20260930T095245Z_MAT-11-W1.log`) and `MAT-13` 75 s (`20260930T172102Z_MAT-13-W1.log`), each including the shared environment. `MAT-11` is steady-only and does not import `transient.py`; it runs because the row names `MAT-9`…`MAT-14`, at ≈ 30 s marginal cost.

  `scripts/testing/run_and_log.sh MAT-18-step2b-W2b-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat11_perfusion_sphere.py tests/validation/test_mat13_transient_perfusion_sphere.py -v -s --tb=short'"`
- **W2c** (real build; rule (c), not the row's). `MAT-17`'s FEM module is item 200's gate on the same file. Its cost is item 200's W1 wall clock; item 200 predicted 80–140 s.

  `scripts/testing/run_and_log.sh MAT-18-step2b-W2c-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat17_duty_cycle_fem.py -v -s --tb=short'"`
- **W3a** (complex build; predicted 100–120 s). `MAT-10` measured 70 s (`20260930T021802Z_MAT-10-step1-W1.log`) and `MAT-12` 77 s (`20260930T110230Z_MAT-12-W1.log`), each including ≈ 44 s of environment and collection. `MAT-10` is steady-only and runs for the row's letter.

  `scripts/testing/run_and_log.sh MAT-18-step2b-W3a-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat10_sar_to_heat_sphere.py tests/validation/test_mat12_transient_sar_heating_sphere.py -s -v --tb=short'"`
- **W3b** (complex build; `MAT-14` measured 79 s, `20261001T004928Z_MAT-14-W1.log`).

  `scripts/testing/run_and_log.sh MAT-18-step2b-W3b-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat14_sar_pennes_transient_sphere.py -s -v --tb=short'"`
- **Not re-run, disclosed:** `MAT-15`'s two FEM modules (`test_mat15_finite_cylinder_fem.py`, `test_mat15_coil_driven_heating_phantom.py`; heavy, 252–273 s), 2(a)'s module, and `examples/materials/04_sar_heats_gel_ball.py`, `07_sar_heats_perfused_ball.py` and `08_gel_ball_settling_time.py` also call `solve_transient_conduction` or the steady solver with scalar coefficients. The evidence for them is the byte-identical scalar path, shown by W2a–W3b's digits.

**Traps.**
- **The positivity check on a Function:** `k <= 0.0` at `transient.py:143` is a truthiness error on a `fem.Function`, and an ordering on complex values in the complex build. Branch scalar or Function; use the real part and allreduce.
- **One ρc object in all three forms** (above). A mismatch reads as a (T2) red of ≈ 0.07, not as a (T1) red.
- **Ghost cells:** set DG0 values from the tags, then `scatter_forward`. `cell_tags.find` and `assemble_scalar` are rank-local, so reduce every count and integral (`MAT-9`'s `_allsum`).
- **Time truncation is not the lever.** `MAT-9`'s CN self-difference at N 40→80 was 4.3898e-5 whole-ball (`…step2b.log:130`), ≈ 200× under the band. A (T1) red is spatial or a coefficient fault: never raise N in-slot.
- **The shell's stiff modes:** at the 0.005 rung, h²/α₂ = 0.005²/(0.19/1.7e6) ≈ 224 s < Δt 348 s, and the source jumps at r = a. CN leaves those modes undamped without the Rannacher start, so (T1) and (E) use `startup_half_steps=2` (even). The no-startup run is printed only.
- **The comparand:** never evaluate it at t = 0 (the series tail goes as 1/N²); t* is Fo 0.2. Its points × 400 arrays are chunked at 20 000 points. It chooses its layer by exact r; points between the faceted interface and r = a get the other layer's formula, which is consistent because θ and the flux are continuous there. Do not special-case it.
- **Masks are built in numpy** (a DG0 indicator from cell midpoints), never with UFL comparisons.
- **The new DG0-coefficient forms compile on first use.** A killed window leaves an orphan `.c` (a 0-byte or partial file). Before any re-run, `pgrep -c python3` in the container must read 0; then run the rubric's sweep as one bare call:

  `scripts/testing/run_and_log.sh MAT-18-step2b-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **Shared files:** `analytical.py` is not edited. If a helper must be added, append it after the last block (item 199's `MAT-18` block; item 202's `IMP-4` block if it has landed) and rebase. `io/mesh.py` is not edited either: a defect in 2(a)'s generator is a stop, never an in-slot edit or a path checkout from a branch.
- **Records:** append `test-results.md` rows by hand; `git diff --stat` must show insertions only on it.
- **Builds:** real-build windows skip the four complex-only environment tests; that is expected.
- **Run hygiene:** print on rank 0; run with `-s`; never `-k` inside the quoted command; no pipes. Rule (i): no edit to `transient.py` or the new test after its last window. A W1b fix that touches `transient.py` re-runs W2a–W3b after it.

**Scope.** The transient with piecewise k and ρc, P1, perfect contact, a uniform core source, from rest, to Fo 0.2. **Not included:**
- contact resistance;
- a real container's properties (k₂ and ρc₂ are representative acrylic, not a material claim);
- the cylinder (`MAT-15`);
- perfusion;
- a SAR or coil-driven source in the two-layer ball;
- degree 2;
- the new module in the complex build (the forms use `ufl.inner(·, v)` and should compile there, but W3a/W3b gate only the scalar path);
- interface-flux extraction.

**Status it can move.** §7 `MAT-18` 🟡 → ✅ on the row's Done-when: "*2(a)* … *2(b)* is the transient: CN at Fo 0.2 ≤ 1 % with piecewise ρc, and `MAT-9`…`MAT-14` re-run green. 2(a) + 2(b) ⇒ ✅." This item matches 2(b) term for term (rule (j)): (T1) is the CN criterion, and W2a, W2b, W3a and W3b cover `MAT-9`…`MAT-14`. ✅ needs 2(a) recorded green by item 201. If 2(a) is on the row only as BLOCKED or parked, the row is held 🟡 whatever this item reads. The chunk's opening and step 1's R0 re-pose are the weekly's to ratify, so the ✅ stands subject to that ratification.

**Negative result.**
- **(D), (T2) or (T3) red:** a code bug. Fix it in-slot, run W1b, re-run W2a–W3b if `transient.py` moved, and disclose.
- **(E) red:** the two-region mesh or the transient comparand path is at fault, not ρc. Print (E) beside 2(a)'s steady (E) reading and the volume deficits, then stop.
- **(T1) red with (D), (T2), (T3) and (E) green:** the piecewise-ρc transient disagrees with the composite series. Print θ_h against the comparand at r = 0, a and b, the shell-only error, and one extra CN N 80 run at 0.005 (≈ 5 s) to separate time from space. Then stop. Step 1's R2 / R3 are green, so suspect the FEM first, but rule nothing in-slot. No band moves (rule (f) never applies to a bound).
- **(c1) ≤ 1e-2 or (c2) ≤ 1e-3:** the shell's ρc does not reach the reading, and the anchor does not discriminate. Report the values and stop.
- **Any W2a / W2b / W2c / W3a / W3b red, or a changed digit:** the scalar path is not byte-identical, so stop.
- **Over the cell guard:** stop before the solve, report the counts.
- On every stop, park on `attempt/MAT-18-step2b-<UTC>` and mark the item BLOCKED in the same commit (rule (d)). The row stays 🟡.
- There is no PREDICTED RED line.

**Outcome (2026-10-03 13:30 slot, `implementer`): DONE, row 🟡 → ✅** (subject to the weekly's ratification of the opening and the R0 re-pose).
- **W1** `20261003T184348Z_MAT-18-step2b-W1.log`: status 124 at the 180 s ceiling, in the comparand (`f.interpolate(callable)` at ≈ 1.1M points per rank). `pgrep` 0; sweep `20261003T184659Z_MAT-18-step2b-cache-sweep.log` found no orphans; probe `20261003T184829Z_MAT-18-step2b-probe-comparand.log` (2.23 s per 20 000-point call). Fix in the new test only: evaluate at `tabulate_dof_coordinates()` as the brief specifies (same P2 interpolant). No band, N or timeout moved.
- **W1b** `20261003T184850Z_MAT-18-step2b-W1b.log` (13 passed / 4 skipped, 119 s): cells 215 542 at 0.005 (guard 400 000), 29 788 at 0.01; volume deficits core −9.15e-4 / shell −2.79e-4 at 0.005. (T1) 8.2066e-04 (pred 5e-4–2e-3); (T2) ≤ 2.599e-14; (T3) 2.703e-15; (D) 5.112e-16; (E) 4.6771e-04 (just under its predicted low end); (c1) 5.3140e-02 (pred ≈ 0.05); (c2) 6.7332e-02 (pred ≈ 0.07). Printed: shell-only 2.22e-3; θ_h vs comparand at r = 0 / a / 0.999b 5.9481 / 5.9448, 3.7877 / 3.7940, 2.3703 / 2.3708 K; no-Rannacher interior 7.67e-4; CN self-difference rate 2.000.
- **Re-gates (digits identical to the cited records):** W2a `20261003T185100Z` (`MAT-9` 4.6327e-04 / 0.984 / 2.001), W2b `20261003T185304Z` (`MAT-13` 1.6426e-03), W2c `20261003T185430Z` (`MAT-17` (i) 4.3870e-04 / 1.0683e-03), W3a `20261003T185646Z` (`MAT-12` 1.6845e-03), W3b `20261003T185824Z` (`MAT-14` 1.6089e-03). Only round-off-level readings moved (`MAT-9` (iv) 4.302e-15 vs 4.301e-15).
- **Disclosed (rule (c)):** `transient.py` routes k / ρc through `steady._as_source` (scalar → same `fem.Constant`), plus `_is_scalar` and `_check_positive_coefficient` (allreduced real minimum for a Function); one ρc object in all three forms. Rule (i): `transient.py` 18:42:42Z precedes W1; the test 18:48:47Z precedes W1b.

### Audit — 2026-10-03 18:00 review

**PASS** (`auditor`, closing commit `c371c5a4`). Evidence: `20261003T184850Z_MAT-18-step2b-W1b.log:88` ((T1) 8.2066e-4), `:98` ((E) 4.6771e-4; balance 2.599e-14), `:102` ((D) 5.112e-16), `:106` ((T3) 2.703e-15), `:93–94` (c2 / c1 miss), `:243–244` (Status 0, 119 s); re-gates W2a–W3b Status 0, each under 180 s. Caveats (no demotion): the steady piece rests on 2(a)'s record (LSQ rate 1.822, pairwise 1.12 / 2.73 / 1.17); (T1) / (E) gate the core interior only, and the shell-only error 2.22e-3 is printed; (E) landed just under its predicted 5e-4 … 2e-3 (print only); ✅ stays subject to the weekly's ratification of the opening and the R0 re-pose. Journal: `docs/planning/reviews/2026-10-03-daily-1800.md` §3.

## §7 row moved verbatim 2026-10-03 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/MAT-18.md`)* **The gel ball in a container wall: conduction through a two-layer sphere (gel core, low-k shell, Robin outside) — Phase 8 (8a), gel, no perfusion** (opened by the 2026-10-02 18:00 daily under rule (5) lane 3; the weekly ratifies) | ✅ **Audit PASS 2026-10-03 18:00 review.** 2026-10-03 (§9 item 205, step 2(b); subject to the weekly's ratification of the opening and the R0 re-pose): DG0 k and ρc in `solve_transient_conduction`, CN + Rannacher N 40 to Fo 0.2 on the 0.005 rung, interior 8.21e-4 ≤ 1 %, balance ≤ 2.6e-14, long BE = steady 2.7e-15, the uniform-ρc controls miss (0.053 / 0.067), `MAT-9`…`MAT-14` and `MAT-17` re-gated with their recorded digits (`20261003T184850Z_MAT-18-step2b-W1b.log`). Step 1 green 2026-10-03 (§9 item 199, attempt 2 of 4, R0 re-posed by the 03:00 review; the weekly ratifies): two-layer Robin sphere closed form with (R0a) 2.09e-13 for n ≤ 60, (R0b) all 400 roots float64-certified, every other step-1 band and all three controls green (`20261003T140149Z_MAT-18-step1-W2.log`). Step 2(a) green 2026-10-03 (§9 item 201, record, no glyph): conforming two-region ball + DG0 k in `solve_steady_conduction`, LSQ rate 1.822, interior 6.79e-4 at 0.003536, balance ≤ 1.02e-14, coefficient path 6.2e-17, equal-properties 7.25e-4, uniform-k control misses 0.199; `MAT-5`/`9`/`10`/`11`/`12`/`13`/`14` re-gated green (`20261003T172548Z_MAT-18-step2a-W1.log`, W2a / W2b / W3a / W3b same stem). Step 2(b) (§9 item 205) is runnable; chunk file. The first thermal material interface (θ and k∂ᵣθ continuous); the wall's series resistance in the gel phantom. **Done-when:** *step 1* (numpy) `two_layer_sphere_steady_temperature` / `two_layer_sphere_transient_temperature`: interface and Robin residuals ≤ 1e-13 (steady) / 1e-10 (transient); global and per-time balance ≤ 1e-10; the n-th eigenfunction has n − 1 interior zeros for n ≤ 60; weighted Gram off-diagonals ≤ 1e-10; analytic vs quadrature norms ≤ 1e-12; equal properties reproduce `robin_sphere_{steady,transient}_temperature(radius=b)` ≤ 1e-13 / 1e-10; three controls miss ⇒ 🟡. *Step 2* (FEM; split 2026-10-03 03:00 review into 2(a) / 2(b) for slot size, bands and criteria unchanged): P1 on a conforming two-region ball with piecewise k and ρc (additive in `thermal/` and `io/mesh.py`). *2(a)* is the mesh and the steady solve: LSQ rate ≥ 1.7 and interior ≤ 1 % at the finest rung, balance ≤ 1e-9 (`MAT-9`'s bands), and the uniform-k control misses. *2(b)* is the transient: CN at Fo 0.2 ≤ 1 % with piecewise ρc, and `MAT-9`…`MAT-14` re-run green. 2(a) + 2(b) ⇒ ✅. Not claimed: contact resistance, a real container's properties, the cylinder (`MAT-15`), perfusion. 

### §9 item 270 — 2026-10-06 (chain P11)

*Drafted for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §10 P11 (quoted, not changed): "`MAT-18` step 2 — the container wall on the coil-driven phantom. `MAT-18`'s wall model (thickness, k_wall, outer Robin) on `MAT-15` P1's fixture. Asserted: the scheme's energy identity with the wall at `MAT-15`'s `E_BAND`, the wall → 0 limit reproducing P1's ⟨θ⟩(t\*) at `B_BAND`, the C4 identity; printed, the wall's effect on ⟨θ⟩(t\*) and on probe B. Target: identities and a limit. Serial on nothing." **Numbering:** `MAT-18` step 2 (2(a) / 2(b)) is already green, so this is filed as **step 3**. A record on a ✅ row: no glyph moves.

`MAT-18` step 3 puts `MAT-18`'s container wall (thickness b − a, k₂, outer Robin h) on `MAT-15` step 1's coil-driven phantom. That fixture is `ANS-2`'s F-small: four single-leg drives at 10 MHz, the 0.0025 m phantom rung, 58 866 tag-3 cells, G = 1 W / P₁ (`test_mat15_coil_driven_heating_phantom.py:180–301`, `_run_heated`).
- **Executor:** implementer. One new module `tests/validation/test_mat18_wall_coil_driven.py`, `@complex_only`, one module-scoped fixture. No `src/` change. **Do not edit** `test_mat15_coil_driven_heating_phantom.py`: `ans:9` imports `_run_heated`, and an edit owes its windows under rule (a).
- **Tier:** heavy (step 1's own window is 218 s at `-n 2`), `-n 2`, complex. **Time:** ≈ 17 measured / ≈ 30 predicted min.
- **Depends on nothing.**

**The wall model on this fixture — read first.** The P1 thermal domain is a submesh of the EM mesh. It has no wall cells, and building a conforming wall layer would be a new mesh (out of scope). So the wall enters as **its series thermal resistance in the Robin coefficient**. That is the steady limit of `MAT-18`'s conforming shell, with no wall heat capacity.
- **h_eff = 1 / (1/h + t_w/k₂), planar.** Import `MAT-18`'s constants (`A_`, `B_`, `K1`, `K2`, `Q1`, `H`) from `test_mat18_two_layer_sphere_closed_form.py` (or from wherever that module imports them): t_w = `B_ − A_` = 0.01 m and k₂ = 0.19, with step 1's `H_GEL` = 10. Predicted h_eff ≈ 6.552 W/m²/K.
- `solve_transient_conduction` takes a scalar `h` (`transient.py:101`), so one h_eff serves the lateral faces and the caps. The lateral cylindrical value, 1 / (a·ln(1 + t_w/a)/k₂ + a/((a + t_w)·h)) at a = `A` = 0.03 m, is **printed beside it** and never used.
- **(R) precondition, asserted — the tie to `MAT-18`.** On the sphere the series reduction is exact in steady state. `two_layer_sphere_steady_temperature(r, a=A_, b=B_, k1=K1, k2=K2, q1=Q1, q2=0, h=H, t_inf=T_INF)` at 20 points r ∈ [0, A_] must equal `robin_sphere_steady_temperature(r, radius=A_, k=K1, q=Q1, h=h_sph)` with h_sph = 1 / (A_²(1/A_ − 1/B_)/K2 + A_²/(B_²·H)), relative ≤ 1e-12. This is numpy round-off; `MAT-18` step 1's steady residual band is 1e-13. It shows the resistance algebra is right before the algebra touches the coil fixture.

**Code.**
- The fixture calls `_run_heated(return_fields=True)` once. That gives the four EM drives, P1's CN runs, the ball centres, and `q_h` and `T_star` keyed by k.
- For each k, run CN with h_eff through a local `_transient_h(sub, ft, q, h, …)`: `_transient`'s body (`:145–151`) with `h` as an argument.
- Use a local `_e_prime_h`: `:164–171` with `H_GEL` replaced by h. Step 1's `_e_prime` hard-codes `H_GEL`, so using it with the wall would fail (E′) on a correct field.
- Ball means at t\* reuse step 1's compiled ball forms: copy `:231–261` into the new module. Do not import `test_*` names.

**Anchors** (asserted; every band imported by name):
- **(R)** above, ≤ 1e-12.
- **(E′) energy identity with the wall:** ρc∫(θᴺ − θᴺ⁻¹)/Δt + ∮h_eff θ^(n+½) = G·P_k, from the EM side, ≤ `E_BAND` (1e-9, `test_mat12_transient_sar_heating_sphere.py:76`) for all four drives. (E): every wall run's `max_balance_residual` ≤ `E_BAND`.
- **(L0) wall → 0 limit:** t_w = 0 (h_eff = `H_GEL`) reproduces P1's volume-mean ⟨θ⟩(t\*) per drive, relative ≤ `B_BAND` (1e-10, `test_mat10_sar_to_heat_sphere.py:61`). The comparand is the mean of `_run_heated`'s `T_star[k]`, the same run that printed 1.153 K per deposited watt (`20261001T123309Z_MAT-15-step1-W1.log:1929`).
- **(C4) with the wall:** the four drives' 10 g ball means of θ at t\*, in the cyclic pairs step 1 forms, each ≤ `C4_COVARIANCE_BAND` (5 %, `test_birdcage_b1_plus_map.py:120`). Step 1 measured 0.020–0.049 % (`…MAT-15-step1-W1.log:1924–1954`). The prediction is the same order with the wall, because a uniform h does not break C4.

**Negative controls.**
- **(d) doubled source, asserted** (prior same-fixture measurement): (E′) with 2·q_h against G·P_k reads 1.000000 on this fixture (`20261001T123309Z_MAT-15-step1-W1.log:1948`). Assert a miss > `D_CONTROL_FLOOR`, as step 1 does.
- **(W-ctl) the wall is felt, predicted (rule (e)).** The wall-on ⟨θ⟩(t\*) misses P1's by more than `B_BAND`.
  - Estimate: the adiabatic ⟨θ⟩(t\*) = t\*/(ρc·V_sub) ≈ 1 254/(4.18e6 × 2.262e-4) ≈ 1.326 K per watt, against P1's 1.153 K. So ≈ 13 % is lost through the surface by t\*.
  - h → 0.655 h cuts that loss by ≈ 20–35 %, so the predicted **wall effect is +2 % … +6 %** on ⟨θ⟩(t\*).
  - Asserted only as a miss > 1e-6. Ceiling: the effect cannot exceed the 15 % adiabatic gap, so 1e-6 sits ≥ 10⁴× under any non-null reading.

**Printed (records):** h_eff (planar) and the lateral cylindrical value; per-drive ⟨θ⟩(t\*) with the wall and its ratio to P1's; the peak 10 g θ̄ with the wall and its ratio; the C4 pairs; max|Im θ|. **Probe B is not printed.** The probe centres live in `MAT-19`'s module on the tuned copper fixture, not on P1's. Deferred, disclosed.

**Window** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first).
- **Price:** step 1's window was 218 s at `-n 2` (`20261001T123309Z_MAT-15-step1-W1.log`, §7 row), and its CN at 0.0025 took 2.5 s (`:1929`). Four wall CN runs, the (L0) run, (d) and the ball forms add ≈ 20–40 s. Predicted 240–300 s.
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0. Then run the any-orphan cache sweep (rubric item 4) as `MAT-18-step3-cache-sweep`.
- **W1:**

  `scripts/testing/run_and_log.sh MAT-18-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat18_wall_coil_driven.py -v -s --tb=short'"`
- No edit after W1 (rule (i)).

**Traps already paid for.**
- Orphan `.c` of any size (`MAT-15` item 160, 44 KB): sweep only at `pgrep` 0.
- Step 1's `_e_prime` and `_transient` hard-code `H_GEL` (above).
- `cell_tags.values` and `assemble_scalar` are rank-local: reduce through `_allsum`.
- No point evaluation on a `create_submesh` mesh (known-issues 2026-10-04, `ANS-9`). Every θ reading is an assembled integral.
- Complex build: put `ufl.real` inside `ufl.lt` for the ball indicator, and keep `sqrt` out of every comparison (`WF-6` 4c). Apply `np.real` before every `float`.
- No `-k` in the quoted command. A killed window is not a result.

**Scope.** Covered: the series-resistance wall (no wall heat capacity) on P1's untuned 10 MHz, σ = 800 S/m fixture; P1 thermal; gel constants; no perfusion. **Not claimed:** a conforming wall on the cylinder, the wall's transient storage, contact resistance, a real container's properties, any tuned-coil, copper, 64 MHz or absolute temperature figure, `ANS-9`.

**Negative result.**
- **(R) red:** the resistance algebra or the import is wrong. Stop before any solve, fix, re-run.
- **(E′) or (E) red:** the h_eff path in the local helper. Print stor / conv / deposited. BLOCKED; `MAT-15` and `MAT-18` ✅ stand.
- **(L0) red:** the local helper is not step 1's run (args, start-up half-steps, snapshots). Fix it in the new module only.
- **(C4) red:** print the q_h pairs. If those are red too, the fixture moved (known-issues entry). Otherwise BLOCKED.
- **(W-ctl) does not miss:** h_eff never reached the solver. Stop.
- No band moves.

**Status it can move.** None: `MAT-18` stays ✅. On green, the §7 row gains a "step 3 record" sentence (h_eff, the wall effect, the identities) and §10 P11 closes. Rule (j): the reviewing daily adds "step 3 (§10 P11): identities green ⇒ a record, no glyph" to the Done-when.

**Outcome — 2026-10-06 15:00 slot (second item): GREEN, a record; `MAT-18` stays ✅.** New module `tests/validation/test_mat18_wall_coil_driven.py` (implementer, `-n 2`, heavy). `20261006T200649Z_MAT-18-step3-W1.log`: Status 0, 266 s (predicted 240–300 s), 31 passed per rank. (R) rel 4.985e-16 (`:80`, h_sph 7.663333). h_eff planar 6.551724138 W/m²/K, used; lateral cylindrical 8.304028064, printed (`:81`). P1 fixture reproduced (P₁ rel 2.598e-14, 58 866 cells; `:1928–1929`). Per drive k = 0…3 (`:1959–1973`): (E′) 6.661e-15 / 2.267e-14 / 1.133e-14 / 9.337e-15; (E) worst 4.554e-14; (L0) 0 / 6.661e-16 / 2.220e-16 / 2.220e-16; ⟨θ⟩(t\*) with the wall / P1 = 1.046314–1.046316 (wall effect +4.631 %, inside the predicted +2 … +6 %; (W-ctl) miss ≫ 1e-6); peak 10 g θ̄ 2.664050 K on k = 0, ratio 1.029608. (d) 1.000000 (`:1974`). C4 θ pairs 0.0402 / 0.0210 / 0.0292 / 0.0484 %, q_h pairs 0.0133 / 0.0084 / 0.0466 / 0.0683 % (`:1975–1978`). max|Im θ| 0. Caveat for the review: (L0) is exact to round-off because h_eff(t_w = 0) returns exactly `H_GEL`, so it proves the local helper is step 1's run, as the brief asks — not an independent limit. Probe B not printed (deferred, disclosed). Sweep: `20261006T200520Z_MAT-18-step3-cache-sweep.log` (Status 0, 84 s, no orphans); a first attempt `20261006T200510Z_…` failed host-side on unescaped `$` (Status 1, never reached the container, nothing deleted).
