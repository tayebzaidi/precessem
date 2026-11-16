# MAT-17 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record. Opened by the 2026-10-02 18:00 daily review under §9 rule (5), lane 3, inside Phase 8 subgoal (8a) (the closed-form lane on the sphere); the weekly ratifies or reverts the opening.

### §9 item 189 — 2026-10-02 (lane 3)

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.*

**`MAT-17` step 1: the duty-cycled Robin gel ball in closed form.** It is verified by residuals and by two independent routes before any FEM.
- **Executor:** implementer. **Tier:** standard, `-n 2`, real build.
- **Cost:** ≈ 8 measured / ≈ 40 predicted min (≈ 1 min window + 15 fixed + ≈ 25 code).
- **Depends on:** nothing.

**Not already covered.** `MAT-9`…`MAT-14` and `MAT-5` use one homogeneous ball with a constant or once-switched source. `solve_transient_conduction` cannot switch a source off ("step source", `thermal/transient.py:20`). No row mentions duty cycles.

**Fixture.**
- `MAT-9`'s ball: R 0.1, k 0.6, ρc 4.18e6, q 2000, h 10 (Bi 5/3), T∞ = T(0) = 20, t_R = R²ρc/k = 69 667 s.
- The source is on while 0 ≤ (t mod P) < DP, from t = 0.
- **Case A:** P = 1200 s, D = 0.5. **Case B:** P = 0.2 t_R, D = 0.25.

**Anchor (closed form).**
- **Modal form:** θ = Σ Bₙ Xₙ(r) gₙ(t). λₙ come from `robin_sphere_eigenvalues`; Bₙ and Xₙ from `robin_sphere_transient_coefficients` (θ_ss = Σ Bₙ Xₙ); sₙ = λₙ²/t_R.
- **Cycle j, t = jP + τ:**
  - start value g_{j,0} = Gₙ(1 − e^{−j sₙ P}), with Gₙ = e^{−sₙ(1−D)P}(1 − e^{−sₙDP})/(1 − e^{−sₙP});
  - on phase: θ = θ_ss(r) − Σ Bₙ Xₙ (1 − g_{j,0}) e^{−sₙτ}, with θ_ss in closed form;
  - off phase: θ = Σ Bₙ Xₙ g_end e^{−sₙ(τ − DP)}.
- **Fourier route** (uses neither λₙ nor Bₙ):
  - θ_p = D·θ_ss + 2 Re Σ_{m≥1} ûₘ Θ(r; ωₘ) e^{iωₘt}, with ωₘ = 2πm/P and ûₘ = (1 − e^{−iωₘDP})/(iωₘP);
  - Θ = (q/w)[1 − Bi·R·S(r)/Dₙ] at w = iωρc, γ = √(w/k). This is `MAT-11`'s perfused steady state continued to complex w.
- **Duhamel:** θ(t) = Σⱼ (−1)ʲ θ_step(t − tⱼ), summed over the switch instants tⱼ.
- **Sources:** Carslaw & Jaeger, *Conduction of Heat in Solids* 2e, ch. I (Duhamel's theorem) and ch. IX (the sphere with heat production and a radiation boundary); Incropera §5.6 for the λₙ equation, as `robin_sphere_eigenvalues` already cites. Equation numbers are unpinned because the repo holds no copy, the same footing as `MAT-9`…`MAT-15`. The docstring derives every form.

**Code** (additive in `utils/analytical.py`; scipy imported inside the functions; nothing exported from `thermal/__init__`, which imports dolfinx):
- `robin_sphere_duty_cycle_temperature(r, t, *, radius, k, rho_c, q, h, t_inf, period, duty, n_terms=400, periodic=False)`;
- `_robin_sphere_complex_w(r, w, *, radius, k, q, h)`, the complex Θ, with `robin_sphere_periodic_harmonic(r, omega, *, radius, k, rho_c, q, h)` on top of it;
- the test, `tests/validation/test_mat17_duty_cycle_closed_form.py`.

**Evaluation set.** r ∈ {0, R/2, R}; eight times per phase, at fractions 0.05 … 0.95 of the phase. Never within 0.05 of a switch.

**Anchors (asserted).**
- **(H0)** `_robin_sphere_complex_w` at a real w > 0 equals `robin_sphere_perfusion_temperature` (t_inf 0) to ≤ 1e-13 relative.
- **(H1) Robin residuals.**
  - −kΘ′(R) − hΘ(R) at m ∈ {1, 10, 1000}: ≤ 1e-12 relative to h|Θ(R)|.
  - The periodic modal form's residual, from analytic derivatives: ≤ 1e-10·qR/h.
- **(P1)** The periodic cycle mean (each phase's time integral taken analytically per mode) equals D·θ_ss(r) to ≤ 1e-12·θ_ss(0), in both cases.
- **(P2)** The periodic modal form matches the Fourier route at M = 2¹⁶ to ≤ 1e-9·θ_ss(0), in both cases (predicted ≲ 5e-11).
- **(P3)** The from-rest form matches the Duhamel superposition of `robin_sphere_transient_temperature` (t_init = t_inf) over cycles 0–5 to ≤ 1e-10·θ_ss(0).
- **(P4)** duty = 1 equals `robin_sphere_transient_temperature` to ≤ 1e-13. duty ∉ (0, 1] raises `ValueError`.
- **(P5)** From rest, at cycle j* = ⌈30τ₁/P⌉, the form equals the periodic one to ≤ 1e-12·θ_ss(0).

**Printed (predicted, never asserted).**
- θ_ss(0) = 12.2222 K.
- λ₁ ≈ 1.91, so τ₁ = t_R/λ₁² ≈ 1.91e4 s.
- **Case A:**
  - centre mean 6.1111 K (exact by (P1));
  - centre ripple peak to peak ≈ q(1−D)DP/ρc = 0.1435 K, so peak/(D·θ_ss) − 1 ≈ +1.2 %;
  - ≈ 110 cycles to come within 1e-3 of the periodic state.
- **Case B:** peak, trough and mean (no hand prediction).
- **The Fourier ladder** M = 2¹⁰…2¹⁶ on the evaluation set and at the switch instants, where the convergence is O(1/M).

**Negative controls.** Rule (e): all are predicted, since nothing has been measured on this fixture. Each is asserted only as exceeding its anchor's band, with the measured factor printed beside the prediction.
- **(c1)** Dirichlet roots λₙ = nπ (Bₙ re-projected) in the periodic modal form, against (P2). Predicted ≈ 0.27·θ_ss(0) at the centre in case A. Asserted > 1e-3.
- **(c2)** The geometric sum dropped (g_{j,0} = 0 every cycle), against (P5). Predicted ≈ 0.49. Asserted > 1e-3.
- **(c3)** The e^{−iωt} convention (Θ at −ωₘ, same ûₘ), against (P2). Predicted ≈ 1e-2 in case A and ≈ 1e-1 in case B. Asserted > 1e-6.

**Window.** Foreground, Bash `timeout: 660000`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Predicted 50–70 s; the three lane-3 numpy windows on 10-01 measured 45–52 s (closest: `20261001T143248Z_MAT-15-step2a-W1b.log`).

`scripts/testing/run_and_log.sh MAT-17-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat17_duty_cycle_closed_form.py tests/unit/test_analytical_lightweight.py -v -s --tb=short'"`

The lightweight test is the rule (c) re-gate of the additive `analytical.py` change.

**Traps.**
- **Series convergence.** Σ Bₙ Xₙ alone converges only as 1/N²: at r = 0 the tail is ≈ 2e-5 K at N = 400. Use the closed-form θ_ss in the on phase and keep the decaying factor in the off phase. s_N ≈ 22.7 /s, so never assert near a switch.
- **Cancellation.** Compute 1 − e^{−sP} with `-np.expm1`.
- **Overflow.** Θ overflows once |γR| reaches a few hundred. Use the scaled form, as `_robin_sphere_perfusion_scaled` does:
  - S(r) = e^{−γ(R−r)}(1 − e^{−2γr})/(2r);
  - S(0) = γe^{−γR};
  - Dₙ = [γR(1 + e^{−2γR}) + (Bi − 1)(1 − e^{−2γR})]/2.

  Check that `np.expm1` accepts complex input, or write the expression out.
- **m = 0** is handled separately as D·θ_ss: q/(iω) is singular there.
- **Sign convention.** Use e^{+jωt} throughout (the repo convention), with ûₘ taken against e^{−iωₘt}.
- **Memory.** Sum the harmonics in blocks of ≤ 8192 to bound memory on two ranks.
- **Output.** Print on rank 0. The real build skips the four complex-only environment tests; that is expected.
- **Rule (i):** no edit after the last window.
- **Rule (f):** a W1 red on a pre-registered *parameter* (δ, a time fraction), never on a bound, may be re-registered by measurement. Disclose it, run W1b, and the review ratifies (the 10-01 lesson, items 136–137).
- **Shared file.** Item 190 also appends to `analytical.py`: append at the end, and rebase if 190 lands first.

**Scope.** Numpy only: a uniform source on `MAT-9`'s ball. Not included:
- the FEM (step 2) and the solver's source schedule;
- an SAR- or coil-driven duty cycle;
- perfusion;
- any pulse-sequence or 6-min-averaging claim.

**Status it can move.** §7 `MAT-17` ⬜ → 🟡 on the row's step-1 Done-when, which this item matches (rule (j)). Step 2 then becomes draftable.

**Negative result.**
- **(H0), (H1), (P1), (P4) or (P5) red:** a code bug. Fix it in-slot, run W1b and disclose.
- **(P2) red with (P3) green, or the reverse:** one route is wrong. Print both routes and the M ladder at the failing points, then stop. Park on `attempt/MAT-17-step1-<UTC>`, mark BLOCKED in the same commit (rule (d)); the row stays ⬜.
- There is no PREDICTED RED line.

**Outcome — DONE 2026-10-03 (19:30 CDT slot, implementer).** W1 `20261003T010113Z_MAT-17-step1-W1.log` passed first time: `-n 2`, 66 s, 26 passed, 4 skipped. No W1b was needed.

Anchors, case A / case B:
- **(H0)** ≤ 1.6e-16.
- **(H1)** The complex Robin residual is ≤ 3.3e-16, and the periodic modal residual is 1.1e-15 / 1.3e-15 of qR/h.
- **(P1)** 1.5e-16 / 0.
- **(P2)** 2.1e-12 / 8.1e-11. Case B is above the predicted ≲ 5e-11 but inside the 1e-9 band.
- **(P3)** 1.2e-15 / 1.0e-15.
- **(P4)** 5.8e-16, and a duty outside (0, 1] raises.
- **(P5)** 4.9e-14 at j* = 479 / 1.1e-14 at j* = 42.

Controls, case A / case B:
- **(c1)** 0.278 / 0.442 (predicted ≈ 0.27 in case A).
- **(c2)** 0.494 / 0.195 (predicted ≈ 0.49, which was case A's figure).
- **(c3)** 1.06e-2 / 9.7e-2 (predicted ≈ 1e-2 / 1e-1).

Records:
- θ_ss(0) = 12.2222 K, λ₁ = 1.90709, τ₁ = 1.9155e4 s.
- **Case A, centre:** mean 6.1111 K, ripple 0.14354 K, peak +1.17 % above D·θ_ss. It takes 111 cycles to come within 1e-3 of D·θ_ss(0), the reading the executor chose because the brief left the reference unstated.
- **Case B, centre:** peak 3.6801 K, trough 2.4042 K, mean 3.0556 K.
- **Fourier ladder:** O(1/M) at the switch instants.

Disclosed:
- The complex 1 − e^{−z} uses a 13-term series below \|z\| = 0.1 rather than `np.expm1`.
- The case B c1 print line reuses case A's prediction text. It is a print only.

Row ⬜ → 🟡; step 2 is draftable.

### §9 item 200 — 2026-10-03 (lane 3)

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.*

**`MAT-17` step 2: a per-step source schedule in `solve_transient_conduction`, gated against step 1's closed form on `MAT-9`'s ball.**
- **Executor:** implementer. **Tier:** standard, `-n 2`; W1 and W2 real build, W3 complex. **Class:** lane-3 closed form (measured median 8 min).
- **Cost:** ≈ 8 measured / ≈ 45 predicted min (15 fixed + ≈ 25 code + three ≈ 2-min windows).
- **Depends on:** nothing. `utils/analytical.py` is not edited by this item; if the executor finds it must add a helper there, append-only after item 199's `MAT-18` block, and rebase.

**Target — the §7 row's step-2 Done-when, verbatim:** "an additive per-step source schedule in `solve_transient_conduction` reproduces step 1 on `MAT-9`'s ball at `MAT-9` step 2's imported bands (interior ≤ 1 %, BE rate ∈ [0.9, 1.2], per-step balance ≤ 1e-9), with `MAT-9`…`MAT-14`'s transient gates re-run green ⇒ ✅". "Default digit-identical" is the item's letter (below).

**Code (additive; rule (c), disclosed).**
- **`thermal/transient.py`:** keyword `source_schedule=None` on `solve_transient_conduction` and `_solve_transient` (`transient.py:73–91`, `:122–140`). It is a callable `χ(t_mid) -> float`, evaluated at the midpoint of every (half-)step, that multiplies the `ufl.inner(q_c, v) * ufl.dx` term of `L` (`:184`) through a new `fem.Constant` `s_c` set before each solve.
  - Step times come from the integer index, never accumulated: half-step i spans [i·Δt/2, (i + 1)·Δt/2]; full step i after the start spans [(m/2 + i)·Δt, (m/2 + i + 1)·Δt].
  - Build the scaled `L` **only when a schedule is given** — the `perfused` pattern (`:171–189`). With `None` every form is byte-identical, so it hits the same JIT cache entry and gives the same digits.
  - Balance under a schedule: `|stor + conv + perf − s_n·src| / |src|`, normalised by the full-on `source_power` (`:199–201`, `:247`), never by `s_n·src` (zero on an off step). The `None` path is unchanged.
  - Additive result field `source_multipliers`: one entry per (half-)step, filled only under a schedule.
- **New test** `tests/validation/test_mat17_duty_cycle_fem.py`.
  - Imported read-only from `test_mat9_transient_conduction_sphere.py`: `R, K, Q, H, T_INF, RHO_C, T_R, ROBIN_TAG, QDEG, _ball, _interior_error, _rel_diff, _say`. From step 1's module: `CASES` (`test_mat17_duty_cycle_closed_form.py:44`).
  - Comparand: `robin_sphere_duty_cycle_temperature(…, periodic=False)` (`analytical.py:1909`), evaluated in numpy at the P2 dof coordinates in chunks of 20 000 points (`MAT-9`'s `_exact_theta`, `:93–106`).

**Fixture.** `MAT-9`'s ball (R 0.1, k 0.6, ρc 4.18e6, q 2000, h 10, T∞ = T(0) = 20), P1.
- **Case A** (P 1200 s, D 0.5): t_end = 11P + 900 s = 14 100 s (Fo 0.2024), the middle of cycle 11's off phase. Every Δt used divides both DP = 600 s and (1 − D)P/2 = 300 s. So every switch is a step boundary, every step midpoint is ≥ Δt/2 from a switch, and χ is exactly 0 or 1.
- **Case B** (P = 0.2 t_R, D 0.25): t_end = 2.5·DP = 8 708.3 s (Fo 0.125), the middle of cycle 0's off phase; Δt = DP/8, so 20 steps.

**Anchors (asserted).** The bands are `MAT-9` step 2's, imported by value: they are literals at `test_mat9_transient_conduction_sphere.py:274–275, 301–302`, with no module constant to import.
- **(S0) Schedule census, every run:** Σ s_n·Δt_n equals the on-time in [0, t_end] to ≤ 1e-12 relative. That is 12 × 600 = 7 200 s for case A, and DP for case B.
- **(D) Default path:** `source_schedule=lambda t: 1.0` against the default call, on the 0.01 m mesh, BE, Δt 300, 47 steps. Whole-ball ‖ΔT‖/‖θ‖ ≤ 1e-13 (predicted ≤ 1e-15).
- **(L) Linearity:** θ[χ] + θ[1 − χ] = θ[1], BE, Δt 150, 0.01 m. Whole-ball relative ≤ 1e-12 (predicted ≈ 1e-15). It is exact because the θ-scheme is linear in the source when T(0) = T∞.
- **(i) Interior:** CN with the Rannacher start (`startup_half_steps=2`) on the 0.005 m mesh — case A at Δt 150 (N 94), case B at Δt DP/8 (N 20). Interior relative L2 on θ (`MAT-9`'s mask \|x\|² < 0.81R²) ≤ 1 % in each case. Predicted 5e-4 – 2e-3; `MAT-9`'s CN N 40 read 4.6327e-4 on this mesh (`20260929T170037Z_MAT-9-step2b.log:114`).
- **(ii) BE rate:** 0.01 m mesh, case A, Δt ∈ {300, 150, 75, 37.5, 18.75} (N 47 … 752). The whole-ball self-difference LSQ rate is in [0.9, 1.2] (predicted ≈ 1.0; `MAT-9`'s BE read 0.984, `…step2b.log:124`).
- **(iii) Balance:** every run's `max_balance_residual` ≤ 1e-9 (predicted ≤ 5e-14).
- **(R) Re-gate:** W2 and W3 run the `MAT-9`, `MAT-12`, `MAT-13` and `MAT-14` transient modules green. The slot also checks that W2's `MAT-9` lines reproduce `…step2b.log:114, 124, 130` digit for digit (interior 4.6327e-04, BE LSQ 0.984, CN 2.001) and journals the comparison. A changed digit on the `None` path is a stop.

**Printed (predicted, never asserted).**
- The CN ladder rate on case A (predicted ≈ 2.0), and the interior error of CN without the Rannacher start.
- The BE N 752 interior error against the series. Predicted ≈ 2.8e-3, the 0.01 m spatial floor: `MAT-9`'s BE N 160 read 2.8063e-3 (`20260929T112146Z_MAT-9-step2.log:122`).
- θ_h(0) beside the closed form at t_end.
- On/off step counts (case A at N 94: 48 on, 46 off).
- The wall clock of every run.

**Negative controls.** Rule (e): all are predicted, since this comparison has never been measured. Each is asserted only as missing its anchor's band, with the measured value printed beside the prediction.
- **(c1) Schedule ignored** — the default call, source always on — against case A's duty comparand at t_end, on (i)'s mesh and Δt. Predicted ≈ 1.0 (range 0.8–1.2): the step response is about twice the D = 0.5 history. Asserted > 1e-2.
- **(c2) Schedule dropped from the balance.** On case A's CN run, take snapshots 93 and 94. Step 93→94 is an off step (midpoint 14 025 s, τ = 825 s). Compute |ρc∫(T₉₄ − T₉₃)/Δt + ∮h(½T₉₄ + ½T₉₃ − T∞) − ∫Q|/∫Q, the form of `MAT-9`'s control (b) (`:285–292`). Predicted 1.000; asserted > 1e-3.
- **(c3) Half-phase shift, printed only:** χ(t − 300 s) against the true comparand. Interior miss predicted 2–4 % (the half ripple, 0.072 K, over θ ≈ 2.5 K). It is not asserted because its ceiling is too close to the band.

**Windows.** Foreground, Bash `timeout: 660000`. Each harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.
- **W1** (real build; predicted 80–140 s). The environment tests take ≈ 30 s. `MAT-9`'s 0.005 m CN N 40 solve took 3.7 s (`…step2b.log:114`) and its five-run BE ladder 3.4 s (`:124`); this module runs ≈ 1 500 ladder steps and ≈ 115 CN steps plus the controls.

  `scripts/testing/run_and_log.sh MAT-17-step2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat17_duty_cycle_fem.py -v -s --tb=short'"`
- **W2** (real build; predicted 110–150 s). `MAT-9`'s two modules measured 115 s with the environment (`…step2b.log:275`), of which steady ≈ 25 s (`20260929T034247Z_MAT-9-step1.log:234`, 54 s with the environment). `MAT-13` measured 75 s (`20260930T172102Z_MAT-13-W1.log:245`).

  `scripts/testing/run_and_log.sh MAT-17-step2-W2-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat9_transient_conduction_sphere.py tests/validation/test_mat13_transient_perfusion_sphere.py -v -s --tb=short'"`
- **W3** (complex build; predicted 110–150 s). `MAT-12` measured 77 s (`20260930T110230Z_MAT-12-W1.log:333`) and `MAT-14` 79 s (`20261001T004928Z_MAT-14-W1.log:309`), each including ≈ 44 s of environment and collection, which the two share here.

  `scripts/testing/run_and_log.sh MAT-17-step2-W3-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat12_transient_sar_heating_sphere.py tests/validation/test_mat14_sar_pennes_transient_sphere.py -s -v --tb=short'"`
- **Not re-run, disclosed:**
  - `MAT-10` and `MAT-11` are steady-only: neither imports `thermal/transient.py`.
  - `MAT-15`'s two modules (`test_mat15_finite_cylinder_fem.py`, `test_mat15_coil_driven_heating_phantom.py`) also call `solve_transient_conduction`. They are heavy (252–273 s) and outside the row's list. The evidence for them is the byte-identical `None` path.

**Traps.**
- **Division by zero on an off step** if the balance is normalised by the scheduled source (above).
- **Never evaluate χ at a step endpoint.** At a switch instant it is ambiguous and shifts the source by one step. That is an O(Δt) error the BE rate would hide.
- **Time from the integer step index.** Case B's P is not binary-exact, but every midpoint is Δt/2 from a switch, so a float `mod` is safe there.
- **The closed form:** never evaluate it within 0.05 of a phase from a switch; both t_end values above are mid-phase. Its (points × 400) arrays are chunked at 20 000 points.
- **The new `L` form compiles on first use.** A killed window leaves an orphan `.c`. Before any re-run, `pgrep -c python3` in the container must read 0; then run the any-`.c`-without-`.so` sweep (daily-review.md rubric item 4), never a 0-byte-only sweep.
- **Builds:** W1 and W2 run in the real build, so the four complex-only environment tests skip; that is expected. W3 needs the complex build because `MAT-12` and `MAT-14` solve the EM field.
- **Run hygiene:** print on rank 0; allreduce every integral (`MAT-9`'s `_allsum`); run with `-s`; never `-k` inside the quoted command.
- **Read-only imports:** do not edit `test_mat9_transient_conduction_sphere.py` or any band in it.
- **Rule (i):** no edit to `transient.py` or the test after its last window.

**Scope.** A uniform source switched on and off on `MAT-9`'s ball from rest, up to Fo ≈ 0.2. **Not included:**
- the periodic steady state in FEM (≈ 110 cycles away);
- SAR- or coil-driven duty cycles;
- perfusion under a schedule (the code path composes, but it is not gated);
- switches off step boundaries (a fractional χ, not gated);
- any pulse-sequence, 6-min-averaging or regulatory statement.

**Status it can move.** §7 `MAT-17` 🟡 → ✅ on the row's step-2 Done-when. This item matches it term for term (rule (j)):
- the schedule is additive and the default is digit-identical ((D) and (R));
- interior ≤ 1 % ((i));
- the BE rate is in band ((ii));
- the per-step balance holds ((iii));
- the `MAT-9`…`MAT-14` transient gates re-run green (W2, W3).

**Negative result.**
- **(S0), (D), (L) or (iii) red:** a code bug. Fix it in-slot, re-run W1 as W1b, and disclose.
- **(R) red** — a changed digit on the `None` path, or any W2 / W3 failure: the `None` path is not byte-identical, so stop. Park on `attempt/MAT-17-step2-<UTC>` and mark BLOCKED in the same commit (rule (d)).
- **(i) red with (L) and (iii) green:** print the interior error at every snapshot, the closed form at r = 0, R/2 and R, and the per-step χ list around the last switch. Park, BLOCKED; no band moves.
- **(ii) outside [0.9, 1.2]:** print the pairwise rates; park, BLOCKED; the review rules. Rule (f) never applies to a bound.
- There is no PREDICTED RED line.

**Outcome — DONE 2026-10-03 (12:00 CDT slot, implementer).** All anchors green; row 🟡 → ✅.

Windows (all `-n 2`, standard, container `timeout -k 30 180`):
- `20261003T170330Z_MAT-17-step2-W1.log` — killed at the 180 s ceiling (exit 124) in the printed-only CN ladder; no assertion had failed. Cause, from the probe's timings: the closed form interpolated per P2 point on the 0.005 m mesh (≈ 37 s per call).
- `20261003T170658Z_MAT-17-step2-cache-sweep.log` (`pgrep -c python3` 0 first; no orphans), `20261003T170842Z_MAT-17-step2-probe.log` (test module only, 3 passed, 135 s).
- Fix (test only): the comparand is evaluated once per distinct radius (`np.unique(…, return_inverse=True)`), ≈ 9.5 s per call; the printed values are digit-identical between W1, the probe and W1b. The timeout was not raised.
- `20261003T171113Z_MAT-17-step2-W1b.log` — 11 passed, 4 skipped, 124 s.
- `20261003T171324Z_MAT-17-step2-W2-regate.log` — 22 passed, 4 skipped, 112 s.
- `20261003T171523Z_MAT-17-step2-W3-regate.log` (complex) — 15 passed, 105 s.

Anchors (W1b):
- **(S0)** relative 0 on every scheduled run (7 200 s / 3 483.333333 s).
- **(D)** 5.700e-16 (default call also leaves `source_multipliers` empty). **(L)** 4.012e-15.
- **(i)** case A 4.3870e-04 (151 548 cells, N 94), case B 1.0683e-03 (N 20).
- **(ii)** BE d = 2.6051e-3 … 3.2832e-4, LSQ 0.996, pairwise 0.993 / 0.997 / 0.998.
- **(iii)** worst 2.358e-13 (printed CN ladder N 752); asserted 0.005 m runs 2.823e-14.
- **(R)** W2's `MAT-9` lines equal `20260929T170037Z_MAT-9-step2b.log:114, 124, 130` digit for digit: interior 4.6327e-04 (balance 1.358e-14), BE d-list and LSQ 0.984 (0.972 / 0.986 / 0.993), CN d-list and LSQ 2.001 (2.002 / 2.000 / 2.000); only wall times differ.

Controls: (c1) 0.97087 (predicted ≈ 1.0); (c2) 1.0000 (predicted 1.000); (c3, printed) 1.4357e-02 — **below its predicted 2–4 %**, recorded as a missed prediction, not asserted.

Printed: CN ladder LSQ 2.104 (2.298 / 2.034 / 2.005); CN without Rannacher 4.2243e-04; BE N 752 vs series 1.8819e-03 (predicted ≈ 2.8e-3); θ_h(0) 3.053003 vs 3.050426 K (A), 1.537399 vs 1.535862 K (B); case A N 94 multipliers 49 on / 46 off over 95 entries (the two Rannacher half-steps counted separately = 48 / 46 full-step).

Disclosed: `thermal/transient.py` additive change (rule (c)); the test also imports `_allsum` from the `MAT-9` module; `MAT-15`'s two modules not re-run (evidence: the byte-identical `None` path).

### Audit — 2026-10-03 18:00 review

**PASS** (`auditor`, closing commit `f320a1f6`). Evidence: `20261003T171113Z_MAT-17-step2-W1b.log:98,101` (interior 4.3870e-4 / 1.0683e-3), `:125` (BE rate 0.996), `:141,144` (balance ≤ 2.358e-13), `:104–105` (controls), `:277–278` (Status 0, 124 s); the `MAT-9` 2b re-gate digits at `…171324Z…W2-regate.log:95,105,111` equal `20260929T170037Z_MAT-9-step2b.log:114,124,130`. Caveats (no demotion): "byte-identical default path" is round-off-identical on the gated quantities, with last-digit balance variation run to run; `MAT-15` was not re-run (the `None` path is unchanged, disclosed); c3 printed 1.44e-2 against a predicted 2–4 %. Journal: `docs/planning/reviews/2026-10-03-daily-1800.md` §3.

### §9 item 271 — 2026-10-06 (chain P12)

*Drafted for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §10 P12 (quoted, not changed): "`MAT-17` step 2 — duty cycle on the coil-driven phantom. `MAT-17`'s on / off schedule on P1's source. Asserted: the Duhamel superposition identity — the duty-cycled response equals the sum of shifted step responses at `MAT-17`'s imported band — and the energy identity; printed, ⟨θ⟩(t\*) against the continuous-source P1 scaled by the duty factor. Target: identities. Serial on nothing." **Numbering:** `MAT-17` step 2 (the FEM `source_schedule`) is already green, so this is filed as **step 3**. A record on a ✅ row: no glyph moves.

`MAT-17` step 3: `source_schedule` (`transient.py:112`; χ evaluated at each (half-)step midpoint, `:294`) on `MAT-15` step 1's coil-driven source — `ANS-2`'s F-small, drive k = 0 (P1), 10 MHz, 58 866 tag-3 cells, G = 1 W / P₁ (`test_mat15_coil_driven_heating_phantom.py:180–301`, `_run_heated`).
- **Executor:** implementer. One new module `tests/validation/test_mat17_duty_cycle_coil_driven.py`, `@complex_only`, one module-scoped fixture. No `src/` change; **do not edit** step 1's module (`ans:9` imports `_run_heated`; an edit owes its windows under rule (a)).
- **Tier:** heavy (step 1's window 218 s at `-n 2`), `-n 2`, complex. **Time:** ≈ 17 measured / ≈ 30 predicted min.
- **Depends on nothing.**

**Schedule (pinned here).** Step 1's CN grid: t\* = 0.2 t_R = 1 254 s, N = `N_MAIN` = 40, Δt = t\*/40, two BE start-up half-steps (`_cn`, `:154–157`). Period P = t\*/5 = 8 steps, D = 0.5: χ(t) = 1 on [jP, jP + P/2), 0 on [jP + P/2, (j + 1)P). Every switch t_m = m·P/2 (m = 0 … 9) is a grid point, so every midpoint χ is exactly 0 or 1. The fixture calls `_run_heated(return_fields=True)` once, then runs through a local `_transient_s(sub, ft, q, schedule)` — `_transient`'s body (`:145–151`) plus `source_schedule=` — with step 1's `startup_half_steps=2` and `snapshot_steps=(N_MAIN − 1, N_MAIN)`. Step 1's `_transient` does not pass a schedule.

**Anchors** (asserted; bands imported by name):
- **(Δ) Duhamel superposition:** θ[χ] − T∞ = Σ_{m=0}^{9} (−1)^m (θ[H_m] − T∞), where θ[H_m] is a **full run from t = 0** with χ_m(t) = 1 for t > t_m. Same scheme, same start-up, so the identity is exact linear algebra, not a shifted-restart approximation. Whole-submesh relative L2 ≤ 1e-12: `MAT-17`'s (L) bar (`test_mat17_duty_cycle_fem.py:201`), pre-stated, measured ≈ 1e-15 there. Ten alternating terms amplify round-off ≈ 10–20×, so the prediction is ≤ 1e-13.
- **(D′) duty 1 is P1:** χ ≡ 1 reproduces `_run_heated`'s `T_star[0]`, relative ≤ 1e-13 (`MAT-17`'s (D) bar, `test_mat17_duty_cycle_fem.py:183`).
- **(E) / (E′) energy identity:** every run's `max_balance_residual` ≤ `E_BAND` (1e-9, `test_mat12_transient_sar_heating_sphere.py:76`; `MAT-17`'s `BALANCE_BAND`, `test_mat17_duty_cycle_fem.py:69`, has the same value). From the EM side on the duty run: |ρc∫(θᴺ − θᴺ⁻¹)/Δt + ∮hθ^(n+½) − χ_N·G·P₀| / (G·P₀) ≤ `E_BAND`. It is normalised by the **unscheduled** deposited power, exactly as `transient.py:310` does. The last step is an off step (χ_N = 0), so a form relative to χ_N·G·P₀ would divide by zero.

**Negative controls.**
- **(d) doubled source, asserted** (prior same-fixture measurement): step 1's (E′) with 2·q_h reads 1.000000 (`20261001T123309Z_MAT-15-step1-W1.log:1948`). On the duty run, assert the doubled-source (E′) misses by > `D_CONTROL_FLOOR`.
- **(Δ-ctl) one term dropped, predicted (rule (e)):** omit m = 9 (the last switch-off). At t\*, θ[H₉] − T∞ carries ≈ (P/2)/(ρc·V_sub) ≈ 125/945 ≈ 0.13 K of mean, against ⟨θ_duty⟩ ≈ 0.55–0.6 K. Predicted miss 0.15–0.3; asserted only as a miss > 1e-12 (ceiling ≈ 0.3, ≥ 10¹¹× the band).

**Printed (records):**
- ⟨θ⟩(t\*) of the duty run against D·⟨θ⟩_P1(t\*), where P1's continuous value is 1.153 K per deposited watt (`…MAT-15-step1-W1.log:1929`). Predicted ratio 0.97–1.00: the deposited energy is the same, t\* ends an off phase, and t_R/P ≈ 25 keeps the ripple small.
- The 10 g peak θ̄ ratio, the ⟨θ⟩ ripple over the last on/off pair, and max|Im θ|.

**Window** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first). Price: step 1's 218 s at `-n 2`. CN at 0.0025 is 2.5 s per run (`:1929`), and this item adds 13 runs (duty, ten H_m, χ ≡ 1, doubled), ≈ 35 s. Predicted 250–300 s.
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0. Then run the any-orphan cache sweep (rubric item 4) as `MAT-17-step3-cache-sweep`.
- **W1:**

  `scripts/testing/run_and_log.sh MAT-17-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat17_duty_cycle_coil_driven.py -v -s --tb=short'"`
- No edit after W1 (rule (i)).

**Traps already paid for.**
- Orphan `.c` of any size (`MAT-15` item 160, 44 KB): sweep only at `pgrep` 0.
- `_transient` does not pass `source_schedule`: use the local helper, never edit step 1's module.
- Keep every switch on a grid point. Off-grid switches make a midpoint χ fractional; (Δ) would survive, but the schedule would no longer be the pinned one.
- Balance normalisation on off steps (above).
- `assemble_scalar` and `cell_tags.values` are rank-local: reduce through `_allsum`.
- No point evaluation on a `create_submesh` mesh (known-issues 2026-10-04, `ANS-9`). `np.real` before `float`. No `-k` in the quoted command. A killed window is not a result.

**Scope.** Covered: a 50 % duty, P = t\*/5 schedule on P1's untuned 10 MHz drive k = 0, P1 thermal, gel constants, Robin h 10, no perfusion. **Not claimed:** any pulse sequence, SAR-averaging or regulatory statement, a closed form for the coil-driven θ, the tuned copper fixture, 64 MHz, `ANS-9`.

**Negative result.**
- **(D′) red:** the local helper is not step 1's run. Fix it in the new module only.
- **(Δ) red with (D′) green:** print the per-term norms and `result.source_multipliers`. A χ that is not 0/1 is a schedule-grid bug in the module. Otherwise open a known-issues entry against `source_schedule` and park BLOCKED; `MAT-17` ✅ stands.
- **(E) or (E′) red:** print stor / conv / χ_N·G·P₀. BLOCKED.
- **A control does not miss:** stop.
- No band moves.

**Status it can move.** None: `MAT-17` stays ✅. On green, the §7 row gains a "step 3 record" sentence and §10 P12 closes. Rule (j): the reviewing daily adds "step 3 (§10 P12): identities green ⇒ a record, no glyph" to the Done-when.

**Outcome — 2026-10-06 15:00 slot (third item): GREEN, a record; `MAT-17` stays ✅.** New module `tests/validation/test_mat17_duty_cycle_coil_driven.py` (implementer, `-n 2`, heavy). `20261006T201547Z_MAT-17-step3-W1.log`: Status 0, 276 s (predicted 250–300 s), 19 passed per rank. Schedule string `11111000011110000111100001111000011110000`, on-time 627 s = D·t\* (`:1956`). (Δ) whole-submesh rel L2 3.421e-14 (`:1970`, bar 1e-12, predicted ≤ 1e-13); (D′) 2.119e-15 (`:1959`); (E′) at N (off) 3.317e-15 normalised by the unscheduled P₀ (`:1957`), at n = 36 (on) 1.399e-14 (`:1958`); (E) worst of 12 runs 3.743e-14. (Δ-ctl) m = 9 dropped: 0.298210 (`:1971`, predicted 0.15–0.3). Records: duty ⟨θ⟩(t\*) 0.568836 K against D·P1 0.576657 K, ratio 0.986437 (`:1973`, predicted 0.97–1.00); last on/off pair 0.466250 → 0.584439 → 0.568836 K (`:1974`); 10 g ball at c₀ duty / (D·P1) 0.967065 (`:1975`, no prediction stated); max|Im θ| 0.

**Departure for the review to rule on — (d).** The brief pins (d) to the duty run's (E′), which is evaluated at N, an off step. There stored + surface does not depend on the source amplitude at that step (linear response, χ_N = 0), so a doubled source cannot be seen at N: the executor measured it there at 2.942e-15 and printed it as blind. Before W1 it moved the asserted (d) to the last on-step n = 36 (χ = 1), against the same `D_CONTROL_FLOOR` and the same prior reading 1.000000 (`…MAT-15-step1-W1.log:1948`), and added the on-step (E′) beside the pinned step-N one. No band moved and nothing was edited after W1. Under rule (h) it is a re-placement of the same comparison, not a substitute true by construction. The review ratifies or reverts it. MAT-17's (L)/(D) bars are restated as module literals with citations, because that module asserts them inline and has no importable names.

## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/MAT-17.md`)* **Duty-cycled RF heating of the Robin gel ball: the on/off source's history and periodic steady state — Phase 8 (8a), gel, no perfusion** (opened by the 2026-10-02 18:00 daily under rule (5) lane 3; ratified by the 2026-10-03 weekly) | ✅ **Audit PASS 2026-10-03 18:00 review.** 2026-10-03 (step 2, §9 item 200): an additive `source_schedule` in `solve_transient_conduction` reproduces the duty-cycle closed form on `MAT-9`'s ball — interior 4.39e-4 / 1.07e-3 (cases A / B, CN + Rannacher, 0.005 m), BE LSQ rate 0.996, balance ≤ 2.4e-13, default path 5.7e-16 and `MAT-9` 2b digits identical, `MAT-9`/`12`/`13`/`14` re-gated green (`20261003T171113Z_MAT-17-step2-W1b.log`, `…171324Z…-W2-regate.log`, `…171523Z…-W3-regate.log`). Step 1 ✅ 2026-10-03 (§9 item 189, `20261003T010113Z_MAT-17-step1-W1.log`). **Step 3 record (2026-10-06, §9 item 271, §10 P12; `20261006T201547Z_MAT-17-step3-W1.log`, 276 s):** the 50 % duty schedule on `MAT-15` P1's coil-driven source keeps Duhamel superposition at 3.4e-14, duty 1 = P1 at 2.1e-15 and the energy identity ≤ 3.7e-14; duty ⟨θ⟩(t\*) / (D·P1) = 0.986 — a record, no glyph; (d) sits on the last on-step (disclosed, the review rules). `MAT-9`'s ball with its uniform source on for D·P and off for (1 − D)·P from T∞; nothing in `thermal/` switches a source off yet. **Done-when:** *step 1* (numpy) `robin_sphere_duty_cycle_temperature` (from rest and periodic) matches its Fourier route (`MAT-11`'s sphere continued to w = iωρc) ≤ 1e-9·θ_ss(0) and the superposition of `robin_sphere_transient_temperature` ≤ 1e-10; the periodic cycle mean equals D·θ_ss ≤ 1e-12; duty 1 reproduces the step response ≤ 1e-13; Robin residuals ≤ 1e-10·qR/h; three controls miss ⇒ 🟡. *Step 2* (FEM): an additive per-step source schedule in `solve_transient_conduction` reproduces step 1 on `MAT-9`'s ball at `MAT-9` step 2's imported bands (interior ≤ 1 %, BE rate ∈ [0.9, 1.2], per-step balance ≤ 1e-9), with `MAT-9`…`MAT-14`'s transient gates re-run green ⇒ ✅. Not claimed: SAR- or coil-driven duty cycles, perfusion, any pulse-sequence or regulatory statement. 
