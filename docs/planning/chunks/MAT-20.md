# MAT-20 — history

Item briefs and execution history. The §7 row (`PROJECT_PLAN.md:2538`) is the status of record. The 2026-10-07 weekly opened the chunk as §10 chain P13 / P14 (`PROJECT_PLAN.md:4657–4668`).

### §9 item 318 — 2026-10-08 (chain P13 — CEM43 dose against its closed form)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 P13 (`PROJECT_PLAN.md:4657–4662`), a weekly-enumerated step (§9 rule (4)). It is serial on nothing. **This brief sets and moves no band, tolerance, Done-when or target.** The 1e-10 below is the §7 row's own bar.

**`MAT-20` step 1 — `thermal/dose.py`, the CEM43 post-processor.**
- **Executor:** implementer.
- **Tier:** smoke. Numpy only: no mesh, no solve, no form compile.
- **Cost:** ≈ 0 measured / ≈ 35 predicted min (15 fixed + 15 code + ≈ 1 compute + spare). No prior measurement exists. The MAT-17 step 1 closed-form module is the comparable case: numpy series, smoke.
- **Ranks:** `-n 2` (rank-safety house rule), even though nothing reduces.
- **New files:** `src/fem_em_solver/thermal/dose.py` and `tests/validation/test_mat20_cem43_closed_form.py`.
- **Command:** `scripts/testing/run_and_log.sh MAT-20-step1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 30 mpiexec -n 2 python3 -m pytest tests/validation/test_mat20_cem43_closed_form.py -v -s --tb=short'"`

**Operator.** `cem43(times, temps, *, r_above=0.5, r_below=0.25, t_ref=43.0)` takes a 1-D series, or a (n_t, n_cells) array (the "one DG0 series": one value per cell per time). Between samples it treats T as **piecewise linear in t** and integrates R^(t_ref − T) exactly on each segment:
- ∫ R^(43−T) dt = (R^(43−T₀) − R^(43−T₁)) / (b ln R), with b = (T₁ − T₀)/Δt.
- A segment that crosses 43 °C is split at t₍₄₃₎ = t₀ + (43 − T₀)/b. R = 0.5 applies on the side T ≥ 43 and R = 0.25 on the side T < 43. The integrand is continuous at 43, where both sides equal 1.
- The b → 0 limit uses the `expm1` form, Δt·R^(43−T₀)·(−expm1(−ΔT ln R))/(ΔT ln R), never the raw difference quotient, which cancels catastrophically when ΔT ≈ 0.
- The §7 row's Σ R^(43−T_i) Δt_i is the piecewise-constant special case. It is printed as a record, never asserted: a rectangle rule cannot reach 1e-10 on a ramp.

**Anchors (asserted at 1e-10 relative — the §7 row's bar).**
- (A1) **Three linear ramps straddling 43 °C** against the analytic piecewise integral above, written out in the test as its own two-line function. The test does **not** call the operator for this; that independence is the point. The ramps, one per branch of the operator:
  - 41 → 45 °C over 600 s (rising);
  - 42.5 → 43.5 °C over 60 s (a shallow straddle);
  - 45 → 40 °C over 900 s (falling, b < 0).

  For a straight line the operator is exact to rounding on any sampling, so each ramp is also run at 2, 7 and 1 000 samples. That catches a split-point bug that only one sampling would hide.
- (A2) **`MAT-9`'s closed-form step response.** The source is `utils.analytical.robin_sphere_transient_temperature` (`analytical.py:1266`) at r = 0, with `MAT-9`'s constants (R = 0.1, k = 0.6, Q = 2000, h = 10; `test_mat9_transient_conduction_sphere.py:41–44`; ρc as `test_mat17_duty_cycle_closed_form.py:41`).
  - The temperatures are re-based to t_inf = t_init = 37 °C. Only the additive offset changes, so θ is unchanged. The centre's steady rise of ≈ 12.2 K means T crosses 43 °C. At MAT-9's 20 °C the response never reaches 43, so the dose would test only the R = 0.25 branch.
  - Reference: `scipy.integrate.quad` of R^(43−T(t)) with `points=[t₄₃]`, where t₄₃ comes from `brentq` on the same closed form, at `epsabs` / `epsrel` ≤ 1e-13. Quad's own error estimate is printed.
  - Operator: run on N uniform samples of the same closed form. **The executor sizes N from a printed N, 2N, 4N sequence**, whose rate should be ≈ 2 for a piecewise-linear interpolant, and asserts at the first N whose predicted error is ≤ 1e-11.
  - Sizing estimate: N ≈ 10⁵ gives an error near 1e-11. That is 4 × 10⁷ series operations at 400 terms, under a second.

**Negative controls (rule (e)).**
- **(Rc) asserted — the planted R swap misses.** With `r_above=0.25, r_below=0.5` the operator must miss (A1) by more than 1e-10 on every ramp. Predicted misses, computed analytically while drafting:
  - 41 → 45: 1.378;
  - 42.5 → 43.5: 0.194;
  - 45 → 40: 1.398.

  The shallow ramp has the smallest separation (≈ 2 × 10⁹ × the bar) and is the binding case.
- **(Sc) predicted, printed — no split at 43.** Using one R across a crossing segment misses (A1). The size is printed and not asserted.

**Traps.**
- **Real build only.** Nothing here needs complex DolfinX, and a `np.real` belongs nowhere in the operator. If a DG0 series from a complex run is ever passed in, the operator raises on a non-zero imaginary part rather than dropping it.
- Pytest captures prints unless `-s` is passed. Never pipe the harness command.
- A `brentq` bracket at t = 0 needs T(0) < 43 < T(t_end). Assert the bracket; do not let it be assumed.

**Public API (rule).** If `dose.py`'s names enter `thermal/__init__.py`'s `__all__` (`thermal/__init__.py:13`), the same slot re-runs `scripts/docs/generate_api_reference.py` and `tests/package/test_wf17_api_reference.py`. The same holds if the slot touches `tests/package/validation_envelope.toml`. That costs one more smoke window. Keeping the module unexported in step 1 avoids it. Step 2 adds to the same module, so either choice is fine.

**Scope.** This is a post-processor against a closed form. It makes no claim to a dose limit, a damage threshold, compliance, tissue or perfusion, and does not apply it to any FEM field.

**Negative result.** A red on (A1) is an operator defect: report it, park the slot on `attempt/MAT-20-step1`, and journal it in attempts.md. A red on (A2) with (A1) green is a sampling or quad-accuracy matter. Print the N-sequence, do not raise N past the smoke ceiling, and stop. Either way the bar stays as it is.

**Status it can move.** The §7 `MAT-20` row goes ⬜ → 🟡 under step 1's Done-when (rule (j)). ✅ waits on step 2.

**Outcome — ✅ DONE 2026-10-08 (12:00 slot).** `MAT-20-step1`, smoke, `-n 2`, 12 s, 8 passed (`docs/testing/logs/20261008T171434Z_MAT-20-step1.log`). `src/fem_em_solver/thermal/dose.py` (`cem43`, not exported from `thermal/__init__.py`) and `tests/validation/test_mat20_cem43_closed_form.py`.
- (A1) the three ramps at 2 / 7 / 1000 samples: relative residuals 0 … 6.3e-16 against the independent two-piece closed form; a (1000, 3) DG0-shaped call ≤ 2.4e-15. The 7-sample rise puts 43 °C exactly on a sample (the touching-end branch). The rectangle sum Σ R^(43−T_i)Δt_i, printed only, is 37.5 / 572.4 / 749.5 against 750.65 on the rise (2 / 7 / 1000 samples).
- (A2) `MAT-9`'s centre response re-based to 37 °C, t_end = 36 000 s (T rises to 47.24 °C; brentq t₄₃ = 14 087.0956 s); quad 1.708516120835631e5 with error estimate 1.1e-14 relative. N = 2¹⁰ … 2¹⁸: relative errors 5.885e-7, 1.471e-7, 3.678e-8, 9.196e-9, 2.299e-9, 5.747e-10, 1.437e-10, 3.592e-11, 8.980e-12, measured rate 2.000 throughout; the Richardson prediction equals the measured error to the printed digits. Asserted at the first N with predicted ≤ 1e-11: N = 262 144, 8.980e-12 ≤ 1e-10. The truncated 400-term series gives T(0) = 37.000000028 (the same function feeds quad and the operator, so it cancels).
- (Rc) asserted: the R swap misses by 1.3784 / 0.1937 / 1.3982 on every sampling — the drafted predictions (1.378 / 0.194 / 1.398) to the printed digits.
- (Sc) printed: one R (the start side) across the crossing segment misses by 1.2973 / 0.1292 / 0.1096.

### §9 item 319 — 2026-10-08 (chain P14 — the six-minute running average)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 P14 (`PROJECT_PLAN.md:4663–4668`), a weekly-enumerated step. **It is serial on item 318, the same module: if 318 did not land, skip.** **This brief sets and moves no band, tolerance, Done-when or target.**

**`MAT-20` step 2 — `running_average(times, values, *, window=360.0)` in `thermal/dose.py`.**
- **Executor:** implementer.
- **Tier:** smoke. Numpy plus the closed form, with no new solve.
- **Cost:** ≈ 0 measured / ≈ 35 predicted min.
- **Ranks:** `-n 2`.
- **New file:** `tests/validation/test_mat20_running_average.py`.
- **Command:** as item 318's, with ID `MAT-20-step2` and that path.

**Operator.** It integrates the piecewise-linear interpolant over [t − W, t] exactly, including fractional end segments. It returns the series of window means, the maximum, and the position t_max. The window must not need to be an integer number of samples: on the FEM grid, 360 s / 31.35 s = 11.48 samples.

**Anchors.**
- (B1) **Integer number of periods.** The input is `MAT-17`'s periodic closed-form history: `robin_sphere_duty_cycle_temperature(..., periodic=True)` (`analytical.py:1943`), case A (P = 1200 s, D = 0.5; `test_mat17_duty_cycle_closed_form.py:44`), at r = 0, with switch times on the sampling grid.
  - The window is W = kP for k = 1 and 2.
  - The window mean must equal D·θ_ss(0) at **1e-12**, the §7 row's bar. θ_ss(0) is `robin_sphere_steady_temperature` and appears as `TH_SS0`, `:48`.
  - MAT-17 P1 reached 1e-12 only by an **analytic** cycle integral (`:146–163`). A sampled operator relies on the trapezoid kink terms cancelling (the jumps in θ′ at on and off are equal and opposite), and that is not guaranteed below O(Δt²). **See question 2.**
- (B2) **The phase of the window maximum.** For W = 360 s on case A, the maximum of the window mean lies where θ(t) = θ(t − W) on the falling side. That point is pre-computed by `brentq` on the closed form, independently of the operator. The operator's t_max must fall within one sample of it, and the root and the offset are printed.
- (B3) **The FEM history at `D_BAND`.** `D_BAND = 1.0e-13` is imported by name from `tests/validation/test_mat17_duty_cycle_coil_driven.py:65` (restated there from `test_mat17_duty_cycle_fem.py:183`). **See question 3 before queueing. As written, this assertion is unreachable on the log's records.**

**Negative controls (rule (e)).**
- **(Wc) asserted — a window shifted by half a period** misses (B1). For W = kP the window mean is phase-invariant, so a shift alone does not move it. The control is therefore **W = kP + P/2**. Its miss from D·θ_ss is O(ripple). The predicted size comes from the closed form's peak and trough, which the MAT-17 P5 records print at `:241`. The miss is asserted > the bar.
- **(Pc) predicted, printed** — the maximum's position under a window one sample short.

**Traps.** These are item 318's, plus three:
- `MAT-17` case B's P = 0.2 t_R is not a multiple of any round Δt. Pin the grid to P/M.
- **The log values are text.** Parse them with a regex anchored on the `MAT17s3   H_m` prefix (`20261006T201547Z_MAT-17-step3-W1.log:1960–1969`). Never hand-copy them into constants without citing the line.
- The log prints every line twice (two ranks' worth of test IDs). Take the first match.

**Scope.** This is a time-average operator against closed forms. It makes no claim to the six-minute SAR/temperature limit itself, compliance or tissue, and says nothing about the coil-driven fixture beyond reproducing a record.

**Negative result.** If (B1) misses 1e-12 on a sampled history, report the M / 2M sequence and stop. Never loosen the bar. Park the slot on `attempt/MAT-20-step2` and record a §7 annotation. A (B3) red is answered per question 3.

**Status it can move.** The §7 `MAT-20` row goes 🟡 → ✅ under step 2's Done-when, provided step 1 is 🟡-landed (rule (j)).

**Questions for the review.**
1. Item 318: the §7 row writes CEM43 as Σ R^(43−T_i)Δt_i but asks for 1e-10 against the analytic ramp integral. The brief resolves this by making the operator exact on piecewise-linear segments, with the Σ form printed. Ratify, or rule otherwise.
2. Item 319 (B1): is 1e-12 meant for a *sampled* operator, or for the analytic cycle integral that MAT-17 P1 already holds? The brief asserts it on samples, with the cancellation argument unproven. A dense M (2¹⁶ per period, as MAT-17 P2 used) may or may not reach it.
3. Item 319 (B3) has a **ceiling**:
   - The W1 log records only end-time values: ⟨θ⟩(t\*) per H_m to 13 significant figures, the duty ⟨θ⟩ at n = 32 / 36 / 40 to 9 decimals, and no per-step series.
   - Measured while drafting: Σ(−1)^m ⟨θ[H_m]⟩(t\*) from `:1960–1969` = 0.5688355135425 K, against the duty record 0.568835514 K (`:1973`). That is a relative difference of 8.0e-10. The print-precision ceiling is 5e-10 / 0.5688 ≈ 8.8e-10.
   - `D_BAND` = 1e-13 is therefore unreachable on the records by ≈ 4 orders of magnitude. Nor can a 360 s window be formed on three duty samples.

   The options are:
   - (a) re-run `_run_duty(return_fields=True)` with every step snapshotted. That is standard, ≈ 276 s for W1's whole window, and the §7 row says "no new solve";
   - (b) the weekly re-poses (B3);
   - (c) leave (B3) out, which holds the row at 🟡.

   The daily may not choose a band. If (a) is chosen, the anchor becomes: on the per-step ⟨θ⟩ series, the operator's window mean of the duty series equals Σ(−1)^m of the operator applied to the H_m series, at `D_BAND`. That identity is linearity, and it is exact to rounding.

**Review rulings (2026-10-08 03:00 daily review, items 318 / 319).** (318-1) **Ratified:** the operator integrates each linear segment of θ(t) exactly — the exact CEM43 of the piecewise-linear interpolant — and the plain sum Σ R^(43−T_i)Δt_i is printed beside it; the Done-when's anchor (1e-10 against the analytic ramp integral) is unchanged. (319-2) No range can be pre-registered for the sampled 1e-12 clause: a red there stops the slot and is reported. (319-3) **The FEM clause is not queued.** The Done-when says "no new solve" and the W1 log holds end-time records only, so `D_BAND` cannot be met from it (the 9-decimal print limits agreement to ≈ 8.8e-10). Option (a) is a new solve the row forbids; the choice between (a) re-posing the clause on a re-run with the per-step series kept, and (b) another comparand, is **the weekly's** (dashboard). **Item 319 executes the closed-form half only** — the integer-period mean = D × θ_ss at 1e-12, the window-maximum phase, and the kP + P/2 control — and records it in the row; `MAT-20` stays 🟡 until the weekly rules on the FEM clause.

**Outcome — ✅ DONE 2026-10-08 (12:00 slot), closed-form half only.** `MAT-20-step2`, smoke, `-n 2`, 19 s, 11 passed (8 step-1 tests re-run against the edited module + 3 new; `docs/testing/logs/20261008T171843Z_MAT-20-step2.log`). `running_average(times, values, *, window=360.0)` added to `src/fem_em_solver/thermal/dose.py` (not exported from `thermal/__init__.py`); test `tests/validation/test_mat20_running_average.py`. Grid: M samples per period, switch times on the grid, `MAT-17` case A periodic closed form at r = 0.
- (B1) max |mean − D·θ_ss(0)| / (D·θ_ss(0)) over every window: k = 1: 5.96e-15 (M = 2¹²), 8.14e-15 (M = 2¹³); k = 2: 1.03e-14, 9.01e-15 — rounding level, no M-trend (the trapezoid kink terms cancel, as the brief argued). Bar 1e-12 unchanged.
- (B2) W = 360 s (1228.8 / 2457.6 samples — fractional): brentq root of θ(t) = θ(t − W) on the falling side 779.999309 s; t_max 779.882812 s (offset −0.398 sample, M = 2¹²) and 780.029297 s (+0.205 sample, M = 2¹³).
- (Wc) asserted: W = kP + P/2 misses by 1.957e-3 (k = 1) and 1.174e-3 (k = 2), inside the drafted O(ripple) ceiling (3.9e-3 / 2.3e-3).
- (Pc) printed: a window one sample short moves t_max by 0 (M = 2¹²) and −1 sample (M = 2¹³) — within one sample, as the half-sample prediction implies.
- A linear ramp (fractional W = 360, 120π, 100 s; 87 samples) is reproduced to ≤ 1.2e-15; the 2-D (DG0-shaped) call equals the column-wise calls exactly.
- (B3) not executed (ruling 319-3): it awaits the weekly. `MAT-20` stays 🟡.

## §7 row moved verbatim 2026-10-10 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/MAT-20.md`)* **Thermal dose and the six-minute window: the regulatory post-processors against closed forms — Phase 8 (8a), gel** (opened by the 2026-10-07 weekly; §10 chain P13 / P14). A `thermal/dose.py` module: *step 1* CEM43 thermal dose on a θ(t) series, `CEM43 = Σ R^(43 − T_i) Δt_i` with R = 0.5 above 43 °C and 0.25 below; *step 2* a running time-average of a series over a fixed window (the six-minute average the exposure limits are written on), returning the window maximum and its position. **Done-when (§4):** *step 1* — for a linear ramp T = a + bt the dose integral is analytic piecewise ((R^(43−a) − R^(43−a−bτ)) / (b ln R) on each side of 43 °C); the operator matches it to 1e-10 relative on three ramps that straddle 43 °C and on `MAT-9`'s closed-form step response, a planted R swap misses ⇒ 🟡. *Step 2* — on `MAT-17`'s duty-cycled closed-form history the running mean over an integer number of periods equals D × the continuous-source value at 1e-12 and the window maximum sits at the pre-computable phase; on `MAT-17` step 3's FEM history (`20261006T201547Z_MAT-17-step3-W1.log`'s records) the operator reproduces the Duhamel sum at `D_BAND`; a shifted-window control misses ⇒ ✅. numpy + one DG0 series; no new solve. Not claimed: a limit value, a compliance statement, tissue, perfusion. **State 2026-10-08: step 1 landed — `cem43` exact on the piecewise-linear interpolant matches the three ramps to ≤ 6.3e-16 and `MAT-9`'s re-based step response to 8.98e-12 at N = 2¹⁸ (rate 2.000), the R swap misses by ≥ 0.194 (`20261008T171434Z_MAT-20-step1.log`, smoke, 12 s). Step 2 closed-form half landed 2026-10-08 (item 319): `running_average` (exact on the interpolant, fractional windows) gives the integer-period mean = D × θ_ss(0) to ≤ 1.03e-14 for W = P, 2P at M = 2¹² / 2¹³ per period, the W = 360 s maximum within 0.40 sample of the brentq root (779.9993 s), the kP + P/2 control misses by ≥ 1.17e-3 (`20261008T171843Z_MAT-20-step2.log`, smoke, 19 s, step 1 re-run green). The FEM clause (`D_BAND` on `MAT-17` step 3's history) is not executed — it awaits the weekly per ruling 319-3; the row stays 🟡.** | 🟡 

### §9 item 370 — 2026-10-11 (chain P19 — `MAT-20` step 2's FEM clause; chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain P19 ("Chain P — Phase 8, gel first", the 2026-10-10 weekly), a weekly-enumerated step (§9 rule (4)). **This brief sets and moves no band, tolerance, Done-when or target.** `D_BAND`, `L_BAND` and the step-2 controls are imported by name.

**Weekly ruling being executed (2026-10-10 weekly, §10 P19), with the drafter's finding.** The weekly wrote: "ruling 319-3, ratified tonight as written: `running_average` on `MAT-17` step 3's FEM history reproduces its Duhamel sum at `D_BAND`, the shifted-window control misses ⇒ `MAT-20` ✅. Smoke, no solve." **The "smoke, no solve" half cannot be executed as written.** The FEM history is not stored anywhere:
- `logs/` on the host has no `MAT-17` or `mat20` JSON / npz (checked 2026-10-11).
- `_run_duty` (`tests/validation/test_mat17_duty_cycle_coil_driven.py:107`) snapshots only `SNAPS` = steps 32, 35, 36, 39, 40 (`:73`). The ten Heaviside runs keep only their end state (`:167–175`).
- The W1 log prints end-time ⟨θ⟩(t\*) per H_m to 13 significant figures and the duty ⟨θ⟩ at n = 32 / 36 / 40 to 9 decimals (`20261006T201547Z_MAT-17-step3-W1.log:1960–1974`). That is three samples of the duty history, and no 360 s window fits in three samples. Ruling 319-3's measured print ceiling (≈ 8.8e-10) stands against `D_BAND` = 1e-13.

Ruling 319-3, as written, names exactly one route to the clause: **option (a)**, a re-run of the step-3 fixture with the per-step series kept, under the anchor "the operator's window mean of the duty series equals Σ(−1)^m of the operator applied to the H_m series, at `D_BAND`; that identity is linearity" (`MAT-20.md`, item 319 question 3). This brief executes option (a) and **re-tiers the item honestly: heavy, `-n 2`, one window ≈ 4–5 min**. It needs the §7 row amendment given under "Status it can move" (rule (j)). **If the review will not amend the row, the item is "write the anchor" for the 10-12 weekly and is not queued.**

**`MAT-20` step 2, FEM clause — the six-minute operator on `MAT-17` step 3's per-step history.**
- **Executor:** implementer.
- **Tier:** heavy, `-n 2` (the step-3 record's width), complex build. **Class:** chain step (solve).
- **Cost:** M measured / ≈ 40 predicted min (15 fixed + 12 code + one window ≈ 4–5 min + one re-run ≈ 5 + spare).
- **Depends on:** nothing open.
- **Gitignored inputs:** none. The fixture is `MAT-15` step 1's mesh, built in-run.

**Code.** One new module, `tests/validation/test_mat20_fem_history.py`: `@complex_only`, one module-scoped fixture, no knob (the `MAT-17` step 3 module's precedent). There is **no `src/` change and no edit to any existing module**.
- **Imports only.**
  - From `tests/validation/test_mat17_duty_cycle_coil_driven.py`: `_transient_s` (`:86`, it already takes `snapshot_steps=`), `_chi_duty` (`:76`), `_heaviside` (`:81`), `N_SWITCH`, `HALF`, `DT`, `STEPS_PER_HALF`, `D_BAND`, `L_BAND` (`:64–71`).
  - From `tests/validation/test_mat15_coil_driven_heating_phantom.py`: `_run_heated` (`:180`), `T_INF`, `N_MAIN`, `T_STAR`.
  - From `tests/validation/test_mat10_sar_to_heat_sphere.py`: `_allsum`, `_say`.
  - From `src/fem_em_solver/thermal/dose.py`: `running_average`, which is unexported, so the import is direct.
- **Fixture body.**
  1. `base = _run_heated(return_fields=True)`. Take `sub`, `ft`, `v_sub` and `q0 = base["q_h"][0]`.
  2. `steps = tuple(range(N_MAIN + 1))`. Run the duty schedule `_transient_s(sub, ft, q0, _chi_duty, snapshot_steps=steps)` and the ten `_transient_s(sub, ft, q0, _heaviside(m), snapshot_steps=steps)`.
  3. At every snapshot, reduce ⟨θ⟩ₙ = ∫(Tₙ − T_INF)/V_sub. Drop each run's snapshots once its series is taken (memory).
  4. Times are tₙ = n·`DT`. Snapshot 1 lands after the two Rannacher half-steps, so the grid is uniform (`transient.py:125–126`).
- **Persist (rank 0):** write `logs/mat20/mat17s3-series-n2.json` (gitignored, `.gitignore:82`). It holds times, the duty series, the ten H_m series, `DT`, `HALF`, `N_SWITCH`, and the rank count. P20 does not need it, but a later example (`mat:18`) can then avoid this solve.

**Anchors (asserted; bands imported).**
- **(H0) The series obeys Duhamel before any averaging:** max over n of |⟨θ_duty⟩ₙ − Σ_m (−1)^m ⟨θ_Hm⟩ₙ| / maxₙ|⟨θ_duty⟩ₙ| ≤ `L_BAND` (1e-12), `MAT-17` step 3's own (Δ) band.
  - Backed on the same fixture: (Δ) read 3.421e-14 in whole-submesh L2 at t\* (`…MAT-17-step3-W1.log:1970`).
  - Predicted ≤ 5e-14.
- **(FEM) The six-minute operator reproduces the Duhamel sum (the row's clause):** with `W = 360.0`, `a = running_average(t, duty, window=W)[0]` and `b = Σ_m (−1)^m running_average(t, H_m, window=W)[0]`, max over valid entries of |a − b| / max|a| ≤ **`D_BAND`** (1e-13).
  - The operator is linear in `values` on a fixed grid, so this is (H0) carried through the operator, plus rounding.
  - W / `DT` = 11.48 samples, so the fractional-segment path is exercised.
  - Predicted ≤ 5e-14. The margin is ≈ 2–3×, which is thinner than the other anchors (see Negative result).
- **(FEM-2d)** The same identity through the operator's 2-D call: an (n_t, 11) array of [duty, H_0 … H_9], with the signed column sum taken afterwards, also ≤ `D_BAND`. Step 2 showed the 2-D call equals the column-wise calls exactly (item 319 outcome).

**Negative control (rule (e); the Done-when's "shifted-window control misses").**
- **(Wc-FEM)** Compare a at tᵢ with the Duhamel sum b at tᵢ₋₄, a shift of P/2 = `STEPS_PER_HALF` samples (item 319's half-period (Wc) analogue).
  - The miss is asserted only > `D_BAND`, the row's form.
  - The factor is *predicted*, printed: ≈ 1e-2–1e-1 relative. The duty ⟨θ⟩ moves 0.466 → 0.584 K over one on-phase (`…:1974`), and a 360 s mean shifted by 125 s moves by a comparable fraction.
  - No same-fixture prior exists.
- **(Dc) Printed only:** b with the m = 9 term dropped, the windowed analogue of `MAT-17`'s (Δ-ctl), which read 0.298 at t\* (`:1971`).

**Printed records.**
- The reproduction of the step-3 log's end-time records, which shows this is "`MAT-17` step 3's history":
  - ⟨θ_Hm⟩(t\*) for m = 0 … 9 against `:1960–1969`;
  - the duty ⟨θ⟩ at n = 32 / 36 / 40 against `:1974` (0.466250330 / 0.584438780 / 0.568835514 K).

  Predicted equal at the printed digits (same code, same width). **Printed, not asserted.** A mismatch beyond the last printed digit stops the item, and it is reported.
- t_max and the maximum six-minute mean of the duty series.
- The maximum six-minute mean of H_0 (P1 alone, rising throughout; t_max predicted at t\*).

**Window.** Foreground; the harness call is the whole Bash command, bare, repo-relative; Bash `timeout: 660000`; `-s`; `tests/environment` first.
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W1** (heavy, `-n 2`, complex):
  - **Predicted harness 230–300 s.** The prior run was `20261006T201547Z_MAT-17-step3-W1.log`, 276 s, with the whole step-3 module. This window has the same `_run_heated` (EM 137 s, total 161 s at `-n 2`, `20261004T142653Z_ANS-9-step1a-W1.log:1895, :1925`), plus 11 transients at ≈ 3.4 s each (`MAT-17` step 3, `:1956`) and ≈ 450 scalar reductions.
  - `scripts/testing/run_and_log.sh MAT-20-step2-fem-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/mat20 && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat20_fem_history.py -v -s --tb=short'"`
- The window is under 10 min, so plain capture suffices (the step-3 precedent at 590 s).

**Traps already paid for.**
- Snapshot index 1 is after **two** BE half-steps (`transient.py:325`). Do not insert a t = Δt/2 sample.
- **Rule (h):** never form b from a itself, or the Duhamel sum from the duty run. The ten H_m runs are the comparand.
- Every reduction is `np.real(_allsum(...))` on every rank. No collective inside `rank == 0`. Only the JSON write is rank-0.
- Pytest captures prints without `-s`. Never pipe the harness command.
- After a kill: `pgrep -c python3`, then the any-orphan cache sweep (daily-review.md rubric item 4) as one bare call.
- **Rule (i):** no edit to the module after W1. An edit re-runs W1.
- No `src/` edit: `dose.py` is not touched, so neither `MAT-20` step-1 nor step-2 module owes a re-run.

**Scope.** A post-processor identity on one FEM history: the F-small coil-driven gel phantom, degree-1 EM, 50 % duty cycle at 1 W. **Not claimed:**
- a six-minute SAR or temperature limit, or compliance;
- tissue, perfusion, or an absolute temperature;
- any history other than `MAT-17` step 3's.

**Negative result.**
- **(H0) red:** the re-run is not step 3's Duhamel history. Stop, park on `attempt/MAT-20-step2-fem-<ts>`, and mark BLOCKED (rule (d)).
- **(FEM) red with (H0) green:** a rounding-ceiling finding on the window means. Print the per-window residual series and stop; `D_BAND` is never widened (ruling 319-2's form). `MAT-20` stays 🟡 and the 10-12 weekly reads it.
- **(Wc-FEM) does not miss:** a wiring error in the shift. Fix it once and re-run W1.
- No `PREDICTED RED` line: every assertion is predicted green.

**Status it can move.** §7 `MAT-20` 🟡 → ✅ under step 2's Done-when, step 1 being landed (rule (j)), **only with this row amendment, made in the review commit that queues the item.** In the row's Done-when, replace:

> on `MAT-17` step 3's FEM history (`20261006T201547Z_MAT-17-step3-W1.log`'s records) the operator reproduces the Duhamel sum at `D_BAND`

with

> on `MAT-17` step 3's FEM history — its fixture re-run at `-n 2` with every step's ⟨θ⟩ kept (ruling 319-3 option (a); the log's end-time records reproduced, printed) — the operator applied to the duty series reproduces Σ(−1)^m of it applied to the ten Heaviside series at `D_BAND`

Also replace "numpy + one DG0 series; no new solve" in the archived row text with "numpy, plus one re-run of `MAT-17` step 3's fixture (no new physics)". The tier stays smoke except for this window, which is heavy.

**Spawn prompt (verbatim sentences).** "Execute §9 item 370 per docs/planning/chunks/MAT-20.md#§9-item-370. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Every Bash call that runs a window passes `timeout: 660000`; run pytest with `-s`; never background a window. No edit to the module after W1. If (H0) or (FEM) is red, report the measurement and stop; never widen a band."

**Same-commit duties.** Mark the §9 item DONE, or BLOCKED with its unblock condition. Update the `MAT-20` row's state line (✅ only if every anchor and the control are green). Add a test-results row through the harness, an attempts.md journal entry, and an outcome paragraph here.

**Outcome — ✅ DONE 2026-10-11 (12:00 slot).** `MAT-20-step2-fem-W1`, heavy, `-n 2`, complex, 270 s harness, 13 passed (`docs/testing/logs/20261011T172810Z_MAT-20-step2-fem-W1.log`). New module `tests/validation/test_mat20_fem_history.py`, imports only; no `src/` or existing-module edit. `_run_heated` 166.8 s (EM 143.3 s); 11 transients + 451 reductions 28.2 s. Series persisted to the gitignored `logs/mat20/mat17s3-series-n2.json`.
- (H0) max_n |⟨θ_duty⟩ − Σ(−1)^m ⟨θ_Hm⟩| / max|⟨θ_duty⟩| = 5.689e-14 ≤ `L_BAND` 1e-12 (predicted ≤ 5e-14: a reading, slightly over the prediction, not a red).
- (FEM) W = 360 s (11.4833 samples), 29 valid windows: max residual 3.434e-14 ≤ `D_BAND` 1e-13 (per-window residuals 4.7e-15 … 3.4e-14, rising with t, printed). (FEM-2d) 3.444e-14; the 2-D duty column equals the 1-D call exactly.
- (Wc-FEM) a(tᵢ) vs b(tᵢ₋₄): miss 1.326e-1 over 25 windows (asserted > `D_BAND`; predicted 1e-2–1e-1, so it landed just above the predicted range). (Dc) m = 9 dropped, printed: 0.043296.
- Records: ⟨θ_Hm⟩(t\*) for m = 0 … 9 equal `20261006T201547Z_MAT-17-step3-W1.log:1960–1969` at all 13 printed digits except H_8 (2.564498589888e-01 vs …887, one unit in the last printed digit: a rounding-boundary print, not a mismatch beyond it). Duty ⟨θ⟩ at n = 32 / 36 / 40 = 0.466250330 / 0.584438780 / 0.568835514 K, equal to `:1974`. Duty six-minute maximum 5.270483664e-01 K at t_max = 1254 s (= t\*); H_0's 1.004390808 K at t_max = t\*, as predicted.
- `MAT-20` 🟡 → ✅ under the amended Done-when (the 03:00 review's row amendment was already in place).

### §9 item 380 — 2026-10-11 (chain P20 — `MAT-20` step 3, dose and the six-minute window on the tuned F-human's heating; chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain P20 ("Chain P — Phase 8, gel first", the 2026-10-10 weekly). **Serial on item 370 (P19), as the weekly wrote it. If 370 did not land, or was not queued, skip.** No file is shared: the seriality is the weekly's ordering, not a code dependency. **This brief sets and moves no band, tolerance, Done-when or target.**

**What the weekly asked, and what exists.** P20 reads "`MAT-15` step 5's probe and ball-mean series, re-read from its JSON". **No such JSON exists, and no θ(t) series was ever computed:**
- `tests/validation/test_mat15_f_human_tuned_heating.py` writes nothing outside the harness log (docstring `:32–34`).
- Its CN runs go through `_cn`, which snapshots only steps N−1 and N (`test_mat15_coil_driven_heating_phantom.py:154–157`).
- Its ball means are taken at t\* only (`:380–480`).
- The W1 log `20261011T003520Z_MAT-15-step5-W1.log` prints t\* values only (`:6311–6314`, `:6384–6385`).
- On the host, `logs/` has no `mat15_s5/` directory; only `logs/mat15-step5-w1-raw.log` exists.

**The step therefore first adds a persistence knob to the step-5 module and regenerates the series. That is a heavy `-n 8` window, ≈ 3 min.** The MAT-20 operators then run on the JSON in a smoke window. Rule (a)/(c): this is an additive edit to an existing gate module, disclosed. Its own window (W1) re-runs the whole module green in the same slot.

**Which "accepted watt".** The 19:30 slot asked this (attempts.md 2026-10-11T00:40Z): P_acc as ruling 361-1 defines it (P_src − P_sheet(P17) − P_sheet(P21)) differs from ½Re Σ V_port I\* by 8.75e-3 relative (`…MAT-15-step5-W1.log:6297`). The weekly then ratified P17 "as written" (§10, the 2026-10-10 weekly's "Ratified" bullet).
- **This step normalises by ruling 361-1's P_acc,ccw**, the module's own G = 1 / P_acc,ccw (`:245`). That is the definition the ratified record carries.
- Every per-watt figure is **also printed per W of ½Re Σ V_port I\***, scaled by P_acc / P_VI, and per W of P_src.
- Both P_acc,ccw and P_VI,ccw are written to the JSON, so the question stays answerable without a re-run.

**Executor and cost.**
- **Executor:** implementer.
- **Tier:** W0 smoke `-n 1`; W1 heavy `-n 8`, complex; W2 smoke `-n 2`, real build. **Class:** chain step (solve).
- **Cost:** M measured / ≈ 45 predicted min (15 fixed + 15 code + W0 ≈ 1 + W1 ≈ 3–4 + W2 ≈ 1 + spare).
- **Gitignored inputs, both checked on the host 2026-10-11:**
  - `logs/wf10/h9p-node-118.2437.json` (5 468 bytes, rewritten by item 361's builder);
  - `logs/wf12/h12c-cu-s32-a2.npz`.

  **If either is absent, stop and mark BLOCKED. Never re-run 358 in-slot.**
- **If W1 has not started by minute 40, park** (rule (d)).

**Code.**
1. **`tests/validation/test_mat15_f_human_tuned_heating.py`: additive knob `FEM_EM_MAT15_STEP5_SERIES`.**
   - Unset ⇒ every print and number is byte-identical and nothing is written. `"1"` ⇒ on, honoured only when `FEM_EM_MAT15_STEP5=1`. Anything else raises `ValueError`, using the module's `_enabled` pattern (`:158–164`).
   - Extend the ungated `test_mat15_s5_coding` to the new parser.
   - Under the knob, after step 5's ball block, when c\*_θ is known:
     - run `_transient(sub, ft, q["ccw"], dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, snapshot_steps=tuple(range(N_MAIN + 1)))`. That is `_cn`'s schedule with every step kept, `ans:9` step 1a's precedent (`ANS-9.md`, item 218, computation step 2). Import `_transient` from step 1's module (`:145`).
     - At each snapshot, compute three series:
       - the **volume mean** ⟨θ⟩ₙ;
       - the **10 g ball mean** θ̄ₙ(c\*_θ), through the module's own `f_th` / `f_vol` forms with `centre_c` set to c\*_θ (`:388–398`);
       - the **probe**: `probe_series(snapshots, T_STAR / N_MAIN, [ProbeBall(c*_θ, A_PROBE)])` from `src/fem_em_solver/thermal/probe.py:174`, with `A_PROBE` = 5e-3 imported from `tests/validation/test_mat19_probe_coil_phantom.py:81`.
     - Rank 0 writes `logs/mat15_s5/theta-series-ccw.json` with:
       - times nΔt and `T_INF`;
       - `G`, `p_acc` (ccw, cw), `p_vi` (ccw, cw), `p_src`, `esr`, `p_phantom` (ccw);
       - c\*_θ, `a10`, `A_PROBE`, the three series, `cells`, `n_sub` and `comm.size`.
   - Print the new block under a `MAT15s5 SERIES` prefix. Its n = N entries must equal the module's own t\* values (mean θ 3.553024444e-02, `:6311`; θ̄(c\*) 4.926111e-02, `:6385`). They are printed, predicted equal to rounding.
2. **New `tests/validation/test_mat20_tuned_heating_dose.py`** (numpy only, real build).
   - It reads the JSON. If the JSON is absent it skips with a reason, which keeps CI green on a fresh clone. **In this slot, W2 must report passed, not skipped.**
   - **Imports:** `cem43`, `running_average` (`thermal/dose.py`); `BAR` as `DOSE_BAR` (1e-10, `tests/validation/test_mat20_cem43_closed_form.py:35`); `BAR` as `AVG_BAR` (1e-12, `tests/validation/test_mat20_running_average.py:34`); `W_B2` (360.0, `:36`).

**Anchors (asserted; bars imported).** All are on the 10 g ball-mean series (the weekly's "ball-mean"), at the stated power P★ = 1 W accepted (361-1), T(t) = `T_INF` + θ̄(t).
- **(Q) Dose = quadrature of the stored series.**
  - Compute `cem43(t, T)` against an independent per-segment `scipy.integrate.quad` of R^(43 − T_lin(t)) on [tᵢ, tᵢ₊₁], summed, with T_lin linear between samples and R = 0.25 (T < 43 throughout; assert that).
  - Relative difference ≤ `DOSE_BAR`. Predicted ≤ 1e-14.
- **(Win) The six-minute maximum sits inside its window.** For `means, mx, t_mx = running_average(t, θ̄, window=W_B2)`:
  - (i) t₀ + W ≤ t_mx ≤ t_N;
  - (ii) the mean-value sandwich: min ≤ mx ≤ max of the interpolant over [t_mx − W, t_mx], evaluated at the fractional start and at every sample inside. This is exact.
  - (iii) mx equals an independent `quad` of `np.interp` over [t_mx − W, t_mx] with `points=` the samples, ≤ `AVG_BAR`.
  - Predicted: t_mx = t\* = 1 254 s (source on throughout; θ̄ rising); (iii) ≤ 1e-15.

**Negative controls (rule (e)).**
- **(Rc′) R swap.** `cem43(..., r_above=0.25, r_below=0.5)` misses (Q)'s reference.
  - Asserted only > `DOSE_BAR`, step 1's (Rc) form.
  - The factor is *predicted*, printed: 2^(43 − T̄) ≈ 8 × 10⁶ relative. No same-series prior exists.
- **(Sh) Printed:** the window one sample short (W − Δt). It gives t_mx and mx.
- **(Rect) Printed:** the rectangle sum Σ R^(43−Tᵢ)Δtᵢ beside (Q), item 318's record form.

**Printed readings (no limit value, no compliance statement).**
- **The ball-mean series:**
  - the six-minute peak θ̄ per accepted W (361-1);
  - the same per W of ½Re Σ V I\* (× P_acc / P_VI, predicted ≈ 0.991);
  - the same per W of P_src;
  - CEM43 at P★ = 1 W under each of the three normalisations, in s and in min.
- **Predicted dose:** ≈ 1.8e-11 s (≈ 3e-13 min). A gel at 20 °C rising by ≈ 0.05 K per accepted W accrues a negligible dose. Printed as such; it is not a safety statement.
- **The same two figures on the volume-mean and probe series**, printed only. A 5 mm probe on a 3 306-cell phantom reads the mesh (P10's red, 2026-10-07 weekly), so nothing is asserted on it.

**Windows.** Foreground; the harness call is the whole Bash command, bare, repo-relative; Bash `timeout: 660000`; `-s`.
- **First:** `pgrep -c python3` must read 0. Then run the any-orphan cache sweep as one bare call, ID `MAT-20-step3-cache-sweep` (daily-review.md rubric item 4).
- **W0** (smoke, `-n 1`): the coding test with both knobs' parsers. The gated tests skip.

  `scripts/testing/run_and_log.sh MAT-20-step3-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -v -s --tb=short'"`
- **W1** (heavy, `-n 8`). Predicted harness 165–200 s: item 361's W1 was 162 s with builder 152.7 s and thermal 4.3 s (`20261011T003520Z_MAT-15-step5-W1.log`). This window adds one CN run (≈ 0.4 s, `:6311`) and 41 × 2 reductions plus one probe geometry. §5.1's durable-capture idiom is copied verbatim, trailing `; exit $rc` included, as item 361 did:

  `scripts/testing/run_and_log.sh MAT-20-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf10 /workspace/logs/mat15_s5 && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_MAT15_STEP5=1 FEM_EM_MAT15_STEP5_SERIES=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -s -v --tb=short > /workspace/logs/mat20-step3-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/mat20-step3-w1-raw.log; cat /workspace/logs/mat20-step3-w1-raw.log; exit \$rc'"`
- **W2** (smoke, `-n 2`, real build; the `MAT-20` step 1/2 precedent, 12–19 s):

  `scripts/testing/run_and_log.sh MAT-20-step3-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 30 mpiexec -n 2 python3 -m pytest tests/validation/test_mat20_tuned_heating_dose.py -v -s --tb=short'"`

**Traps already paid for.**
- **(S_rep)** needs 358's JSON read *before* the builder rewrites it. The module already does this (`:205–214`); do not move the series block above it.
- `quad`'s default `epsabs` (1.49e-8) is ≈ 10⁶ × the integrand (≈ 1.4e-14). Pass `epsabs=0`, `epsrel=1e-13`, and integrate per segment, or the reference is noise.
- `probe_series` is collective: call it on every rank with identical `ProbeBall`s. Snapshots are absolute T (`transient.py:126`), so subtract `T_INF` for θ.
- Snapshot index 1 is after two BE half-steps; times are nΔt.
- **Never apply G twice.** `q["ccw"]` already carries G (ruling 361-4, `:311`).
- `cem43` takes °C. Pass T = `T_INF` + θ, never θ.
- **Rule (i):** no edit to either module after its window; an edit re-runs that window. Rule (a)/(c): W1 *is* the step-5 gate's re-run. All 25 of its prior tests must stay green, with their readings equal to item 361's at the printed digits.

**Scope.** These are operator identities on FEM data from one fixture: the tuned copper F-human, degree-1 EM, the cylindrical 0.015 m gel phantom at 20 °C, t\* = 1 254 s. **Not claimed:**
- a dose or six-minute limit, or compliance;
- tissue, perfusion, or σ(T);
- the ASTM phantom, or degree 2;
- an absolute temperature;
- `ANS-9`.

**Negative result.**
- **W1 red on any prior step-5 anchor:** the knob leaked. Revert it, report, and mark BLOCKED.
- **(Q) or (Win)(iii) red:** an operator defect on real data (step 1/2's closed forms were green). Park on `attempt/MAT-20-step3-<ts>`, file a known-issues entry, and mark BLOCKED.
- **(Win)(i)/(ii) red:** an operator bug. Same disposition.
- No `PREDICTED RED` line.

**Status it can move.** A record on the `MAT-20` row with no glyph, so rule (j) does not bind. The row has no step 3. The review adds to the row's Done-when, in the same commit that queues this item:

> *Step 3 (record, §10 P20):* on `MAT-15` step 5's 10 g ball-mean θ(t) series (regenerated with every step kept), `cem43` equals an independent quadrature at the step-1 bar and the six-minute maximum sits inside its window (step-2 bar); dose and six-minute peak per accepted watt printed — no glyph.

The result feeds the `mat:18` example the weekly asked for after `MAT-20` ✅ (§10, "an example after `MAT-20` ✅ (P20's dose figures, `mat:18`)"). `MAT-15` holds ✅.

**Spawn prompt (verbatim sentences).** "Execute §9 item 380 per docs/planning/chunks/MAT-20.md#§9-item-380; skip it if item 370 is not DONE. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Every Bash call that runs a window passes `timeout: 660000`; never background a window; run pytest with `-s`. If an input under logs/wf10 or logs/wf12 is missing, stop and mark the item BLOCKED. Report a red and stop; never widen a band."

**Same-commit duties.** Mark the §9 item DONE or BLOCKED. Add a state-line clause on the `MAT-20` row (step 3 record) and an outcome paragraph here. Add a journal entry that names which accepted-watt definition the printed figures use.

**Review ruling 2026-10-11 (03:00 daily) on items 370 / 380.** (1) Item 370 is re-tiered heavy (`-n 2`): `MAT-17` step 3's per-step history is not stored and the log's printed digits cannot reach `D_BAND`, so the fixture is re-run once with every step kept (ruling 319-3 option (a), ratified by the 10-10 weekly); the row's Done-when is amended under rule (j) to name the data source — the identity and `D_BAND` are unchanged. (2) Item 380 keeps the weekly's serial order on 370. It normalises per accepted watt by ruling 361-1's definition (P_src minus the two feed sheets) and prints every figure per W of ½Re ΣV·I\* and per W of P_src beside it (journal ruling 2: the definition is the 10-12 weekly's to confirm). The step-3 clause is a record — no glyph.

**Outcome — ✅ DONE 2026-10-11 (13:30 slot).** Inputs present; `pgrep` 0; cache sweep `20261011T183056Z_MAT-20-step3-cache-sweep.log`. Code: additive knob `FEM_EM_MAT15_STEP5_SERIES` in the step-5 module (parser `_series_enabled`, coding test extended; the series block sits after the PRICE line so no prior print moves) and new `tests/validation/test_mat20_tuned_heating_dose.py`; no `src/` edit.
- W0 `20261011T183406Z_MAT-20-step3-W0.log`: 5 passed, 20 skipped, 5 s.
- W1 `20261011T183420Z_MAT-20-step3-W1.log` (heavy, `-n 8`, 160 s): all 25 step-5 tests passed. 72 of 85 `MAT15s5` lines are byte-equal to item 361's log; the other 13 differ only in residual figures at rounding ((P_w) 4.919e-13 vs 4.965e-13, (E) worst 3.310e-10 vs 3.269e-10, (S2) |ratio−1| last digit). Series block 0.7 s. The n = N entries equal the module's t\* values (mean 3.553024444e-02 K, rel 4.5e-14; ball 4.926110850e-02 K, rel 7.4e-14); the probe reads 4.987727e-02 K (vol 5.386e-7 vs nominal 5.236e-7 m³).
- W2 `20261011T183717Z_MAT-20-step3-W2.log` (smoke, `-n 2`, real build, 3 s, 3 passed), on the 10 g ball-mean series at c\*_θ, per accepted W (ruling 361-1):
  - (Q) `cem43` 1.846214193739e-11 s vs per-segment `quad` rel 1.750e-16 (≤ `DOSE_BAR` 1e-10; predicted ≤ 1e-14); max T 20.0493 °C < 43.
  - (Win) W = 360 s (11.48 samples): max six-minute mean 4.277005e-02 K at t_max = t\* = 1 254 s (predicted); sandwich [3.613e-02, 4.926e-02] held; (iii) vs `quad` 3.245e-16 (≤ `AVG_BAR` 1e-12).
  - (Rc′) R swap misses by 8.241e6 (> 1e-10; predicted ≈ 8e6). (Rect) −8.53e-4 vs (Q). (Sh) W − Δt: 4.335e-02 K at t\*.
  - Readings, printed (no limit, no safety statement): ball CEM43 1.846e-11 s = 3.08e-13 min per accepted W; per W of ½Re ΣV·I\* (× 0.991324) peak 4.240e-02 K, 1.846e-11 s; per W of P_src (× 0.335728) peak 1.436e-02 K, 1.803e-11 s. Volume mean: peak 3.070e-02 K, 1.828e-11 s; 5 mm probe: 4.330e-02 K, 1.847e-11 s (per accepted W; not asserted).
- The dose barely depends on the normalisation (≤ 2.4 % across the three), since R^(43 − T) is set by T_INF at a 0.05 K rise; the six-minute peak scales 1:1 with it, which is what the 10-12 weekly's confirmation of 361-1 bears on.

## §7 row moved verbatim 2026-10-12 (rotate_plan_archive.py chunks)

 Thermal dose and the six-minute window: the regulatory post-processors against closed forms — Phase 8 (8a), gel | ✅ 2026-10-11 (§9 item 370) — step 1 landed 2026-10-08 (`cem43` exact on the piecewise-linear interpolant matches three ramps to ≤ 6.3e-16 and `MAT-9`'s step response to 8.98e-12); step 2's closed-form half landed 2026-10-08 (integer-period mean = D × θ_ss(0) to ≤ 1.03e-14). Step 2's FEM clause landed 2026-10-11: `MAT-17` step 3's fixture re-run at `-n 2` with all 41 steps' ⟨θ⟩ kept; (H0) Duhamel on the series 5.689e-14 (≤ `L_BAND`), `running_average` (W = 360 s, 11.48 samples) on the duty series vs Σ(−1)^m on the ten H_m series 3.434e-14 (2-D call 3.444e-14) ≤ `D_BAND`; the P/2-shifted control misses by 1.326e-1; the step-3 log's end-time records reproduced at the printed digits (H_8 differs by one unit in the 13th digit) (`20261011T172810Z_MAT-20-step2-fem-W1.log`, heavy, 270 s). Step 3 (record, §9 item 380) landed 2026-10-11: `MAT-15` step 5's module gained the additive knob `FEM_EM_MAT15_STEP5_SERIES` and re-ran green at `-n 8` (25 passed, 160 s heavy; prior readings equal item 361's at the printed digits, residual figures at rounding), writing the ccw CN θ(t) series (41 steps) to the gitignored `logs/mat15_s5/theta-series-ccw.json`; on the 10 g ball-mean series at c\*_θ, per accepted W (ruling 361-1), `cem43` equals per-segment `quad` to 1.750e-16 (≤ 1e-10), the 360 s maximum 4.277e-2 K sits at t\* = 1 254 s inside its window, equal to `quad` to 3.245e-16 (≤ 1e-12), the R swap misses by 8.241e6; CEM43 1.846e-11 s (printed; 1.846e-11 / 1.803e-11 s per W of ½Re ΣV·I\* / P_src) — a record, no glyph, no safety statement (`20261011T183420Z_MAT-20-step3-W1.log`, `20261011T183717Z_MAT-20-step3-W2.log`, smoke 3 s). **Done-when:** *step 1* — for a linear ramp T = a + bt the dose integral is analytic piecewise ((R^(43−a) − R^(43−a−bτ)) / (b ln R) on each side of 43 °C); the operator matches it to 1e-10 relative on three ramps that straddle 43 °C and on `MAT-9`'s closed-form step response, a planted R swap misses ⇒ 🟡. *Step 2* — on `MAT-17`'s duty-cycled closed-form history the running mean over an integer number of periods equals D × the continuous-source value at 1e-12 and the window maximum sits at the pre-computable phase; on `MAT-17` step 3's FEM history — its fixture re-run at `-n 2` with every step's ⟨θ⟩ kept (ruling 319-3 option (a), ratified 2026-10-10; the log's end-time records reproduced, printed; amended 2026-10-11 under rule (j) because the log's printed digits cannot reach `D_BAND`) — the operator applied to the duty series reproduces Σ(−1)^m of it applied to the ten Heaviside series at `D_BAND`; a shifted-window control misses ⇒ ✅. *Step 3 (record, §10 P20):* on `MAT-15` step 5's 10 g ball-mean θ(t) series (regenerated with every step kept), `cem43` equals an independent quadrature at the step-1 bar and the six-minute maximum sits inside its window (step-2 bar); dose and six-minute peak per accepted watt printed — no glyph. *History: `docs/planning/chunks/MAT-20.md`.* 

### §9 item 557 — 2026-10-20 (§10 P40, `MAT-20` step 4; lane 1, chain no-solve)

*Drafted read-only for the 2026-10-20 03:00 daily review; the review reviews it.* Licence: §10 chain P40 ("Chain P — Phase 8, gel first", the 2026-10-19 weekly, `PROJECT_PLAN.md` ≈ :4841–4848), under §9 rule (4); the weekly's entry-42 ruling ("first branch fires … the dose at degree 2 is **P40** (lane 1, no solve) on the saved series", rulings bullet 1) is its reason. **Serial on nothing.** **This brief sets and moves no band, tolerance, Done-when or target.** Template: item 380 above (step 3, P20) — this is the same operator pair on the same fixture's series at the degree-2 root.

**Why `MAT-20` and not `MAT-15`.** The step applies `MAT-20`'s operators through step 3's module pattern (its helpers are imported); the data is `MAT-15` step 6's (entry 42). Item 380 filed the degree-1 twin here, so the degree-2 twin is step 4 of this chunk; the `MAT-15` row gains a pointer (below).

**The input and its provenance (checked on the host 2026-10-20 03:00).**
- `logs/mat15_s5/theta-series-ccw-d2.json` (gitignored; 4 500 bytes, mtime 2026-10-19 02:09) — written by `xl` entry 42 (`docs/testing/xl-pending.md` § 42; `docs/testing/xl-ledger.md` row 2026-10-19), `docs/testing/logs/20261019T070008Z_MAT-15-step6-P22-d2root.log:6887–6929` (`MAT15s5 SERIES` block; `:6929` "wrote logs/mat15_s5/theta-series-ccw-d2.json (rank 0)"). Keys: `degree` 2, `node_knob` "d2root", `c_in_pf` 115.5209, `ranks` 16, 41 samples at Δt 31.35 s to t\* = 1 254 s, `G_per_w` 36.6406, P_acc,ccw 2.72921e-2 W, P_src 8.41662e-2 W.
- The degree-1 twin `logs/mat15_s5/theta-series-ccw.json` (item 380; 4 422 bytes, mtime 2026-10-11), its records in this file's item-380 Outcome and `20261011T183717Z_MAT-20-step3-W2.log`.
- **If the degree-2 file is absent, the step records that and stops — entry 42 is not re-run** (the weekly's words). The module skips with a reason on a fresh clone; **in this slot the window must report passed, not skipped**; a skip ⇒ mark the item BLOCKED with "series absent" and stop.

**Executor and cost.**
- **Executor:** implementer.
- **Tier:** one smoke window, `-n 2`, **real build** (numpy on a saved series; item 380's W2 ran 3 s, `20261011T183717Z_MAT-20-step3-W2.log`). The weekly's heading says "lane 1 → heavy, no solve": the heavy part is the `xl` window that wrote the series; no window here needs more than smoke.
- **Class:** chain no-solve. **Cost:** ≈ 6.6 measured / ≈ 25 predicted min (module ≈ 12, window ≈ 1, journal and row ≈ 10).
- **Depends on:** nothing. No FFCx form is compiled, so no cache sweep is owed.

**Code: one new module, `tests/validation/test_mat20_tuned_heating_dose_d2.py`** (numpy only; prefix `[MAT-20 s4]`). Edit no landed module (not step 3's, not `MAT-15` step 5's).
- **Imports:** `cem43`, `running_average` (`src/fem_em_solver/thermal/dose.py:45`, `:108`); `BAR` as `DOSE_BAR` (1e-10, `tests/validation/test_mat20_cem43_closed_form.py:35`); `BAR` as `AVG_BAR` (1e-12, `tests/validation/test_mat20_running_average.py:34`) and `W_B2` (360.0, `:36`); from step 3's module `tests/validation/test_mat20_tuned_heating_dose.py`: `SERIES_JSON` (the degree-1 path, `:47`), `R_BELOW`, `T_REF` (`:48–49`), `_quad_dose` (`:78–89`), `_rect` (`:92–93`). Never import its `series` fixture or a `test_*` name.
- `SERIES_D2 = SERIES_JSON.with_name("theta-series-ccw-d2.json")`; a module-scoped fixture loads both files (d2 absent ⇒ `pytest.skip`; d1 absent ⇒ the d1 → d2 prints are skipped with a reason, the d2 identities still run).
- **Provenance strings** (module constant `E42_LOG`, citing `20261019T070008Z_MAT-15-step6-P22-d2root.log:6889` and `:6928`): n = 1 — mean `9.408837938992e-04`, ball `1.293205450941e-03`, probe `1.287310587223e-03`; n = N — mean `3.507742275686e-02`, ball `4.916453675956e-02`, probe `4.979676417684e-02` (each the log's `:.12e`).
- A helper `_is_d2_entry42(d)` returns True iff `d.get("degree") == 2`, `d.get("node_knob") == "d2root"`, `d.get("c_in_pf") == 115.5209`, 41 samples, and the six `f"{x:.12e}"` strings above equal `E42_LOG`.

**Anchors (asserted; every bar imported).** The weekly's target is "the operator identities at `MAT-20`'s imported bands on the degree-2 series". Step 3 asserted them on the 10 g ball-mean series; the weekly's −1.27 % is the **volume-mean** series (entry 42's "ball-mean θ(t\*)" is the phantom volume mean — label trap below), so both are asserted here, at the same bars (an execution choice; the review may strike the volume-mean rows to printed). T(t) = `t_inf_c` + θ, per accepted watt (ruling 361-1; the JSON already carries G).
- **(P0-d2) Provenance.** `_is_d2_entry42(d2)` is True; `d2["times_s"] == d1["times_s"]` and `t_inf_c` equal (list equality) when d1 is present.
- **(Q-d2)** On each asserted series: `cem43(t, T)` against `_quad_dose(t, T, ln R_BELOW)` ≤ `DOSE_BAR`; max T < `T_REF` asserted. Predicted ≤ 1e-14 (step 3: 1.750e-16).
- **(Win-d2)** On each asserted series, for `means, mx, t_mx = running_average(t, θ, window=W_B2)`: (i) t₀ + W ≤ t_mx ≤ t_N; (ii) the mean-value sandwich over the interpolant on [t_mx − W, t_mx] (step 3's `:123–126` construction, exact); (iii) mx against an independent `quad` of `np.interp` on the window with `points=` the samples, ≤ `AVG_BAR`. Predicted: t_mx = t\* = 1 254 s on both (both series rise monotonically, entry-42 log `:6888–6928`); (iii) ≤ 1e-15 (step 3: 3.245e-16).

**Negative controls (rule (e)).**
- **(Rc′-d2) asserted:** the R swap `cem43(..., r_above=0.25, r_below=0.5)` misses (Q-d2)'s reference by > `DOSE_BAR` on each asserted series. Factor *predicted*, printed: ≈ 8.2e6 (step 3 measured 8.241e6 on the degree-1 ball series; T moves by < 1e-3 K between orders).
- **(P0-ctl) asserted:** `_is_d2_entry42(d1)` is False (the degree-1 file has no `degree` key and its n = N mean prints `3.553024444118e-02`) — the provenance check discriminates. Skipped with a reason only if d1 is absent.

**Printed readings (a record; no limit, no compliance statement).** For each of the three series (volume mean, 10 g ball mean at c\*_θ, 5 mm probe), per accepted watt, d2 beside d1:
- CEM43 at P★ = 1 W accepted (s and min), the **d1 → d2 dose move** D₂/D₁ − 1, and beside it θ(t\*)'s move (volume mean −1.27 %, 10 g ball −0.196 %, the weekly's figures) and the six-minute peak's move.
- **The dose excess over the baseline** D − R^(43 − T∞)·t\* and its d1 → d2 move. Reason (predicted, and the reason this print exists): at 1 W accepted the gel rises ≈ 0.05 K from 20 °C, so ≈ 96.5 % of the dose is the baseline R^23·t\* = 1.782041e-11 s and the dose itself moves ≈ 1e-5 – 3e-4 while θ moves 1.27 %; the excess carries the move. Not a target.
- Per W of ½Re ΣV·I\* (× P_acc / P_VI) and per W of P_src (× P_acc / P_src) — the second labelled "not the pre-registered normaliser" (the −4.65 % label trap, xl-ledger 2026-10-19).
- The (Rect) rectangle sum beside (Q-d2); (Sh) the window one sample short.

**Predicted values** (stdlib, exact per-segment integral on the two JSON files, 2026-10-20; d1 reproduces item 380's records 1.846214e-11 s / 4.277005e-02 K):

| series | CEM43 d1 → d2 (s) | dose move | excess move | six-minute peak d1 → d2 (K) | peak move | θ(t\*) move |
|---|---|---|---|---|---|---|
| volume mean | 1.827699e-11 → 1.827100e-11 | −3.28e-4 | −1.313 % | 3.070252e-02 → 3.030951e-02 | −1.280 % | −1.2745 % |
| 10 g ball (c\*_θ) | 1.846214e-11 → 1.846197e-11 | −9.2e-6 | −0.026 % | 4.277005e-02 → 4.270632e-02 | −0.149 % | −0.196 % |
| 5 mm probe | 1.846939e-11 → 1.846950e-11 | +6.0e-6 | +0.017 % | 4.329679e-02 → 4.324857e-02 | −0.111 % | −0.161 % |

**Window.** Foreground; **the harness call is the whole Bash command** — bare, repo-relative, nothing chained or piped; Bash `timeout: 660000`; pytest `-s`. It runs step 3's module beside the new one, so the degree-1 records print in the same log (rule (a): a re-run of a landed module, unedited):

`scripts/testing/run_and_log.sh MAT-20-step4-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 mpiexec -n 2 python3 -m pytest tests/validation/test_mat20_tuned_heating_dose.py tests/validation/test_mat20_tuned_heating_dose_d2.py -v -s --tb=short'"`

Predicted 3–10 s; step 3's three tests must pass with item 380's printed digits.

**Traps already paid for.**
- **Label trap — "ball-mean".** Entry 42's readout and the weekly call the phantom **volume mean** (`theta_mean_k`) "ball-mean θ(t\*)"; step 3's asserted series is the **10 g ball mean** (`theta_ball10g_k`). Print the JSON key beside every figure.
- **Label trap — per P_src.** The θ move per P_src is −4.65 % (the matching term); it is not the pre-registered normaliser and is printed under its own label.
- The d2 JSON's `source` string still reads "(MAT-20 step 3, item 380)" — the knob's text, unchanged by the degree knob. Provenance is (P0-d2), never the `source` string.
- **Never apply G again:** the series are already per accepted W (ruling 361-4; `G_per_w` differs between the files, 36.64 vs 34.42, because P_acc moved).
- `cem43` takes °C: pass `t_inf_c` + θ, never θ.
- `quad`'s default `epsabs` (1.49e-8) is ≈ 10⁶ × the integrand: `_quad_dose` already passes `epsabs=0`, `epsrel=1e-13` per segment — import it, do not re-write it.
- The probe series' n = 0 entry is 3.55e-15 K, not 0 (entry-42 log `:6888`): never assert θ(0) = 0 on the probe; the probe is printed only (P10's red — a 5 mm probe on a 3 306-cell phantom reads the mesh).
- `ranks` is 16 in the d2 JSON and 8 in d1: the files are rank-independent records (rank 0 wrote them); nothing here reduces across ranks; print on rank 0 only (`_p`).
- pytest captures prints without `-s`; a pipe after pytest inside the harness command records the pipe's status — none here.
- Rule (i): no edit to the module after the window starts; an edit re-runs it.

**Scope.** Operator identities and a printed order move on one fixture's saved FEM series: the tuned copper F-human at its degree-2 root (115.5209 pF), the cylindrical 0.015 m gel phantom at 20 °C, t\* = 1 254 s, per accepted watt. **Not claimed:** a dose or six-minute limit, compliance or a safety statement; tissue, perfusion, σ(T); the ASTM phantom; an absolute temperature; a convergence statement (two element orders on one mesh); any external check.

**Negative result.** No `PREDICTED RED:` line.
- **Series absent / skipped:** record "series absent" in the journal and the item line, mark BLOCKED, stop. Entry 42 is not re-run in-slot or queued by the implementer.
- **(P0-d2) red:** the file on disk is not entry 42's series (overwritten or partial). Print every key and string that differs, mark BLOCKED, stop; the review decides (a known-issues entry if the file was overwritten by a later run).
- **(Q-d2) / (Win-d2) red:** an operator finding on real data (step 3 was green on the same fixture at degree 1). Park on `attempt/MAT-20-step4-<ts>`, write a known-issues entry, mark BLOCKED.
- **Step 3's module red in this window:** a regression on `main` unrelated to this item — report it with the failing node ID, mark BLOCKED.
- **A move that differs from the table:** printed, never asserted; the weekly reads it.

**Status it can move.** A record on the ✅ `MAT-20` row — no glyph, so rule (j) does not bind. The row's state line gains one clause: "Step 4 (record, §10 P40): on entry 42's degree-2 series (gitignored `theta-series-ccw-d2.json`), (Q) / (Win) green at `DOSE_BAR` / `AVG_BAR` on the 10 g ball-mean and volume-mean series; the d1 → d2 dose move per accepted watt printed beside θ's −1.27 %." The `MAT-15` row's state line gains "dose at the degree-2 root: `MAT-20` step 4 (P40)", closing the entry-42 readout's `cem43` / six-minute row (UNCOUNTABLE from entry 42's module, the 10-18 ruling). `MAT-15` and `MAT-20` hold ✅.

**Done-when (this item).** The window is green with (P0-d2), (Q-d2), (Win-d2), (Rc′-d2) and (P0-ctl) asserted, not skipped, and the move table printed; the two row clauses written; a journal entry naming the log and the measured moves.

**Spawn line:** "Use the implementer agent to execute `MAT-20` §9 item 557 per `docs/planning/chunks/MAT-20.md#§9-item-557` (read the latest Review ruling paragraph first, if one exists). Write one new module, `tests/validation/test_mat20_tuned_heating_dose_d2.py`; edit no landed module and nothing in `src/`; never import a `test_*` name or step 3's `series` fixture. If `logs/mat15_s5/theta-series-ccw-d2.json` is absent, record that and stop (BLOCKED) — never re-run entry 42. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run the window in the foreground with Bash `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Report a red and stop; never widen a band."

**Same-commit duties.** Mark §9 item 557 DONE or BLOCKED; the two row clauses; an Outcome paragraph under this heading; a journal entry in `docs/testing/attempts.md`; the harness log and its `test-results.md` row (insertions only).

**Outcome — DONE 2026-10-20 (06:00 slot), first try.** Input present (`theta-series-ccw-d2.json`, 4 500 bytes, mtime 2026-10-19 02:09). Code: new `tests/validation/test_mat20_tuned_heating_dose_d2.py` only (imports step 3's `SERIES_JSON`, `R_BELOW`, `T_REF`, `_quad_dose`, `_rect`; no `test_*` name, no fixture); no landed module and nothing in `src/` edited. Window `20261020T110139Z_MAT-20-step4-W1.log` (smoke, `-n 2`, real build, harness 5 s, 10 passed — step 3's three with item 380's digits unchanged, plus seven new, none skipped).
- (P0-d2) `_is_d2_entry42(d2)` True; `times_s` and `t_inf_c` equal d1's. (P0-ctl) `_is_d2_entry42(d1)` False (no `degree` / `node_knob` / `c_in_pf`; all six log strings differ).
- (Q-d2) 10 g ball [`theta_ball10g_k`]: CEM43 1.846197272499e-11 s vs per-segment `quad` rel 0.000e+00; volume mean [`theta_mean_k`]: 1.827099791245e-11 s, rel 1.768e-16 (≤ `DOSE_BAR` 1e-10; predicted ≤ 1e-14); max T 20.0492 / 20.0351 °C < 43.
- (Win-d2) W = 360 s: ball 4.270632e-02 K, volume mean 3.030951e-02 K, both at t_max = t\* = 1 254 s (predicted); sandwiches held; (iii) vs `quad` 1.625e-16 / 3.434e-16 (≤ `AVG_BAR` 1e-12; predicted ≤ 1e-15).
- (Rc′-d2) R swap misses by 8.2412e6 (ball) / 8.2843e6 (volume mean) (> 1e-10; predicted ≈ 8.2e6). (Rect) −8.508e-4 / −6.074e-4 vs (Q). (Sh) W − Δt: 4.328e-02 / 3.073e-02 K at t\*.
- Move table, per accepted W (ruling 361-1), d1 → d2, printed (baseline R^23·t\* = 1.782041e-11 s):

| series [key] | CEM43 d1 → d2 (s) | dose move | excess move | six-minute peak d1 → d2 (K) | peak move | θ(t\*) move |
|---|---|---|---|---|---|---|
| volume mean [`theta_mean_k`] | 1.827699e-11 → 1.827100e-11 | −3.279e-4 | −1.3127 % | 3.070252e-02 → 3.030951e-02 | −1.2801 % | −1.2745 % |
| 10 g ball [`theta_ball10g_k`] | 1.846214e-11 → 1.846197e-11 | −9.165e-6 | −0.0264 % | 4.277005e-02 → 4.270632e-02 | −0.1490 % | −0.1960 % |
| 5 mm probe [`theta_probe_k`] | 1.846939e-11 → 1.846950e-11 | +6.031e-6 | +0.0172 % | 4.329679e-02 → 4.324857e-02 | −0.1114 % | −0.1614 % |

  Every move equals the brief's predicted table at its printed digits. Per W of ½Re ΣV·I\* (× 0.977588 d2, × 0.991324 d1) the volume-mean θ(t\*) moves −2.64 %, the ball −1.58 %; per W of P_src (not the pre-registered normaliser; × 0.324265 / × 0.335728) −4.65 % / −3.60 %. The dose itself barely moves because ≈ 96.5 % of it is the 20 °C baseline; the excess carries θ's move. A record: no limit, no compliance or safety statement, no convergence claim (two orders on one mesh).

### Review ruling 2026-10-20 03:00 (item 557)

Binds over the brief. P40 is filed under `MAT-20` (step 4) because item 380, its degree-1 twin, lives here and the module reuses step 3's helpers; the weekly tied it to the `MAT-15` row, whose state line gains "dose at the degree-2 root: `MAT-20` step 4" when it lands. Asserting (Q) / (Win) on the volume-mean series as well as the 10 g ball series is accepted (an execution choice at unmoved bars). The dose move is expected to be small (≈ 96.5 % of the dose is the 20 °C baseline); the brief prints the excess over baseline beside it — a record, not a target.

### §9 item 638 — 2026-10-25 (§10 lane 1, R3: `MAT-16`'s kink and `expm1` readings against `MAT-20`'s `cem43`)

*Written by the 2026-10-25 10:30 daily review from the 10-24 weekly's lane-1 entry R3 ("optional: `MAT-16`'s (K) kink and `expm1` reading against `MAT-20`'s implementation, printed"). The weekly wrote "printed"; this brief adds no band — every assertion below is `MAT-20`'s own `DOSE_BAR`, imported.*

**Question.** `MAT-16`'s draft (never opened; `docs/planning/chunks/MAT-16.md`) named two hazards of a CEM43 integrator: the 43 °C kink, and cancellation in `(exp(x) − 1)/x` at small slope. `MAT-20` landed `cem43` without a test at either: no test sits at a ramp of half-width ε ≪ 1 K across 43 °C, none at a near-zero slope, and none on a constant segment (the `_phi` `x == 0` branch). Does `cem43` hold `DOSE_BAR` there?

**Facts (read 2026-10-25).**
- `src/fem_em_solver/thermal/dose.py`: `_phi(x)` = `expm1(x)/x`, 1.0 at exactly `x == 0.0` (`:32–37`); `_piece(t_a, t_b, dur, ln_r, t_ref)` = `dur·exp(ln_r·(t_ref − t_a))·_phi(−ln_r·(t_b − t_a))` (`:40–42`); `cem43(times, temps, *, r_above=0.5, r_below=0.25, t_ref=43.0)` (`:45`), output in the units of `times`. Crossing is strict, `(T0 − 43)(T1 − 43) < 0` (`:89`); a non-crossing segment's side is `0.5(T0 + T1) >= 43` (`:96`).
- `tests/validation/test_mat20_cem43_closed_form.py`: `BAR` = 1e-10 (`:35`, imported elsewhere as `DOSE_BAR`); `_side` / `ramp_dose` (`:50–60`) is the **raw difference quotient** `(r^(43−ta) − r^(43−tb)) / (b·ln r)` — itself the naive form `MAT-16` warned of, accurate at ε = 1e-3 (cancellation ≈ 1e-13) but not at ε ≈ 1e-9.

**Route.** New module `tests/validation/test_mat20_kink_expm1.py`, numpy + `math` only, real build, prefix `[MAT-20 R3]`. Import `cem43` (`dose.py:45`), `BAR as DOSE_BAR` and `ramp_dose` from the closed-form module (helpers, never a `test_*` name). No `src/` edit, no edit to any landed module.

**Anchors (asserted at `DOSE_BAR`, imported).**
- **(K) the kink at ε = 1e-3.** `cem43([0, τ], [43 − ε, 43 + ε])` against `ramp_dose(43 − ε, 43 + ε, τ)` (the split two-piece closed form) ≤ `DOSE_BAR`, τ = 600 s; and the same on the descending ramp. **Printed beside it:** mean = dose/τ, its miss from `MAT-16`'s 1 − (ln 2/4)ε, and the predicted remainder 5(ln 2)²ε²/12 ≈ 2.0e-7 (this review's expansion, not asserted — `MAT-16`'s 1e-6 bar is not adopted).
- **(E) the small-slope limit.** On one non-crossing segment T₀ = 41 → T₀ + δ over Δt = 600 s, δ ∈ {1e-6, 1e-9, 6e-10}: `cem43` against the second-order Taylor form Δt·R^(43 − T₀)·(1 + x/2 + x²/6), x = −δ·ln R (truncation ≈ x³/24 < 1e-27) ≤ `DOSE_BAR`; the same at T₀ = 45 (R = 0.5).
- **(Z) the constant segment.** `cem43([0, Δt], [T, T])` = Δt·R^(43 − T) at T = 41 and 45, ≤ `DOSE_BAR` (the `x == 0` branch).
- **`MAT-20` step 1's module re-run unedited in the same window** (its A1 / A2 / R-swap assertions stay green).

**Controls (rule (e)).**
- **(Kc) R swapped**, *predicted*, printed: `cem43(…, r_above=0.25, r_below=0.5)` on the ε-ramp moves the mean by ≈ +(ln 2/2)ε ≈ 3.5e-4 relative. Printed beside the measured miss; not asserted (no prior measurement on this fixture).
- **(Ec) the naive quotient**, *predicted*, printed: `ramp_dose(T₀, T₀ + δ, Δt)` against the Taylor form at δ = 6e-10 — predicted ≈ 2e-7 relative (stdlib arithmetic, ≈ 2.3e-7 at Δt = 600 s), i.e. ≈ 2 000× `DOSE_BAR`. Printed with the ratio; if it does **not** miss, say so — that would mean the trap is not live on this platform, a finding.
- The asserted, backed control is step 1's `test_rc_r_swap_misses` (re-run, unedited).

**Window (bare harness call, foreground, Bash `timeout: 660000`).** One smoke window, `-n 1`, real build, no `mpiexec`:
`scripts/testing/run_and_log.sh MAT-20-R3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 30 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat20_cem43_closed_form.py tests/validation/test_mat20_kink_expm1.py -v -s --tb=short'"`
(the smoke rule: only `test_complex_mode.py` from `tests/environment`). No FFCx form compiles, so no cache sweep is owed.

**Tier:** smoke, `-n 1`. **Class:** lane-1 revalidation. **Cost:** ≈ 7.4 measured / 20 predicted slot-minutes. **Executor:** implementer. **Depends on:** nothing. **Family:** R3, attempt 1 of 4.

**Traps.** Inputs are absolute °C, not θ. Do not import from `thermal/__init__.py`'s exports for anything but `cem43`'s module path. `_piece` / `_phi` are private: call `cem43` only (the test is of the public function). Write `≈`, never a tilde, in the journal. Never pipe the pytest call.

**Scope boundary.** A numerical-hazard read of one post-processor. Not claimed: Arrhenius damage, any threshold, any biological or regulatory statement, `MAT-16`'s other draft items (C3 exponential, C4 trapezoid convergence — `cem43` is exact per linear piece, so they do not apply).

**A negative result means:** a red (K), (E) or (Z) is a defect in `cem43` at that input: stop, report the measured miss, park the module on `attempt/*`, add a known-issues entry naming the input; never loosen `DOSE_BAR`.

**Status it can move.** §10 R3 → DONE. The `MAT-20` row's state line gains "R3: the 43 °C kink at ε = 1e-3, the small-slope limit and the constant segment at `DOSE_BAR` (§9 item 638)". `MAT-16.md` gains a one-line closing note that its (K) and `expm1` hazards were read against `MAT-20`. `MAT-20` stays ✅.

**Outcome (2026-10-25, 15:00 CDT slot, DONE).** `tests/validation/test_mat20_kink_expm1.py`; window `20261025T200112Z_MAT-20-R3-W1.log` (smoke, `-n 1`, real, 14 s harness-wall, 12.69 s pytest; 21 passed, 3 skipped — the complex-mode environment tests on the real build; step 1's module re-run unedited, `test_rc_r_swap_misses` green). **(K)** rise / fall rel 2.255e-14 / 2.274e-14 against `ramp_dose` (bar 1e-10); mean 0.999826913345066, its miss from 1 − (ln 2/4)ε = 2.0014e-7 against the predicted remainder 5(ln 2)²ε²/12 = 2.0019e-7. **(E)** at δ = 1e-6 / 1e-9 / 6e-10, T₀ = 41 and 45: rel 0, 0 / 1.895e-16, 0 / 1.895e-16, 0 against the Taylor form (truncation ≤ 1.1e-19). **(Z)** rel 0.000e+00 at T = 41 and 45. **(Kc)** *predicted, printed*: the R swap moves the mean by +3.466337e-4 against the predicted +3.465736e-4 (ratio 1.000173), both directions. **(Ec)** *predicted, printed*: the naive quotient misses the Taylor form at δ = 6e-10 by 2.218e-9 at T₀ = 41 (22× `DOSE_BAR`; predicted ≈ 2e-7 — the miss is ≈ 100× smaller than predicted, a rounding-luck matter of where 41 + 6e-10 lands) and 2.645e-7 at T₀ = 45 (2 645× the bar); `cem43` itself misses by ≤ 1.9e-16 at both. **Finding:** the trap is live on this platform at both T₀, but the T₀ = 41 magnitude was over-predicted by ≈ 100×. No defect in `cem43`.
