# PORT-24 — capacitor ESR (circuit layer + in-model series-RC sheet); §7 row: `grep -n '^| \`PORT-24\`' PROJECT_PLAN.md` (`PROJECT_PLAN.md:2716`)

### §9 item 308 — 2026-10-08 (chain T25 — `PORT-24` step 1, circuit-layer ESR, no solve)

**Outcome (2026-10-08 04:30 implementer slot): BLOCKED as a record — the pre-registered `PREDICTED RED: (Im-bitwise)` fired inside its range [1e-18, 1e-13].** Run `docs/testing/logs/20261008T094657Z_PORT-24-step1.log` (6 s, smoke, `-n 1`). Code is parked on `attempt/PORT-24-step1-20261008T0948Z` (`src/fem_em_solver/ports/ring_modes.py` + `tests/unit/test_port24_esr_circuit.py`). Measured: (Q2) local DFT helper vs `k_blocks` 0.000e+00 (band 1e-14). (B) every `ring_mode_rows` / `_v2` / `terminated_currents` output on the stored F-human 32×32, at `SPEC_C_F` and `C2_RECORD_F`, is bitwise equal to pinned `b2174873`; the F-small `tuned_input` r = 0 path is bitwise on 12 (record, termination) pairs; `test_ans8_private_comparison.py` is green. (S)-Re: max \|ΔRe λ − r\|/r 1.724e-13 over 33 (record, k), band 1e-12, green (F-human k = 2 the worst; F-small ≤ 5.3e-15). **(S)-Im bitwise: red. Max \|ΔIm\|/\|Im\| is 2.889e-15 and \|ΔC_mode\|/C_mode is 4.571e-16, both on F-human (`np.linalg.eig` of the shifted 2×2 k-blocks). Every F-small value is exactly 0 (the DFT diagonal is a direct sum).** (P) planted: min \|ΔRe − r\|/r 3.125e-02 (F-human, ≈ 0.20625 Ω, as predicted) and 0.25 on F-small, against the 1e-12 bar, green. (T) printed: off by 1.0. Unrelated: three `test_wf17_api_reference.py` tests fail on clean `main` (stale `load_study` / `run_study_file` pages, known-issues 2026-10-08). `ring_modes` has no API page, so no regeneration was due.

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 T25 (`PROJECT_PLAN.md:4579–4585`), weekly-enumerated (§9 rule (4)); binding ruling §9 `PROJECT_PLAN.md:3735–3754`, adoption §10 `:4411–4418, 4513–4522`. **No band, tolerance, Done-when or target is set or moved.** Target (weekly's): closed-form identity.

- **Executor:** implementer. **Tier:** smoke, `-n 1` (pure NumPy; no MPI reduction exists to hide). **Class:** chain step (no solve). **Cost:** ≈ 30 predicted min (≈ 20 code, < 1 compute, rest fixed); no prior measurement of this module's runtime. **Depends on:** nothing.

**Settled from the code.**
- `ring_mode_rows` (`src/fem_em_solver/ports/ring_modes.py:130`), `ring_mode_rows_v2` (`:186`, composes v1, passes `esr_ohm` through to every inner call incl. controls `:206`), `terminated_currents` (`:116`, `y_load = jωC` at `:121`; feed rows `yin` `:163`; feed 2×2 `jwc·I` `:172`). λ_k is the antisymmetric eigenvalue of the k-block of `Ẑ = W Z Wᴴ` (`k_blocks` `:74`, `antisymmetric_eigenpair` `:88`); `C_mode = 1/(ω Im λ₁)` `:147, 180`.
- 4-port circuit terminations (F-small): `_capacitor_impedance` / `tuned_input` (`tests/validation/test_port_circuit_layer_field.py:355, 365`, via `series_rlc_impedance(f, c_f=C)` — already takes `r_ohm`, `src/fem_em_solver/ports/lumped.py:138`) and the `ans:7` SPEC circuit layer `z_cap` / `add_c1` / `route_b_s` / `route_b_z` (`examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py:167–185`). Each gains `esr_ohm=0.0` (keyword, additive).
- Stored records: F-small 4×4 `rungs["10 MHz"|"64 MHz"|"128 MHz"]["s_matrix"]` in `examples/ansys_benchmarks/ans4_…/metrics.json` and `ans6_…/metrics.json` (tracked); F-human degree-1 32×32 `logs/ans8/s32.npz` (gitignored, present on this box) with its tracked rows `ans8_…/metrics.json` keys `rows`, `rows_v2`.
- New module `tests/unit/test_port24_esr_circuit.py`; edits only the functions above.

**Anchors (asserted).**
- (B) default bitwise: with `esr_ohm` unset, every output of the edited functions on each stored record equals the **pre-change** function's output bitwise (the pre-change source loaded from `git show <sha>^:path` pinned, never `HEAD:`), and `tests/unit/test_ans8_private_comparison.py` stays green (its tracked-file byte-identity check).
- (S) at r = 0.2 Ω: Re λ_k(r) − Re λ_k(0) = r for every k, at 1e-12 relative; Im λ_k and `C_mode` bitwise unchanged — on the F-human 32×32 (k = 0…8) and on the F-small 4×4 at 10 / 64 / 128 MHz (λ_k there is the N = 4 DFT diagonal through a local helper, per ruling Q2 below).
- (P) planted non-uniform r (one gap at 0.4 Ω, the rest 0.2): the uniform shift fails on every mode.

**Negative controls (rule (e)).**
- (P) **predicted, asserted as the Done-when's control:** a bottom-gap plant shifts the bottom ring's k-domain diagonal by mean(r) = 0.2125 Ω, the top's by 0.2, so the antisymmetric eigenvalue moves ≈ 0.206 Ω, ≈ 3 % of r off the identity — ≈ 1e10 above the 1e-12 bar (ceiling: the full 0.2 Ω plant excess). No prior measurement; predicted.
- (T) **predicted:** an implementation that wires r into the termination only (`y_load`) leaves λ_k unchanged and fails (S) by exactly r (100 %). This is the reason r must enter the loop impedance (Z + diag r), not only the load.

**Traps already paid for.** Bitwise default needs a branch, not arithmetic: `1/(0 + 1/(jωC))` is not bitwise `jωC`, so the r = 0 path must keep the old expression; pin the pre-change commit (`test_orphan_guard.sh`, `OPS-43`); pytest needs `-s` for the printed rows; never `-k a or b` inside the quoted command; harness foreground.

**Harness (W1, after the complex-mode check):**
`scripts/testing/run_and_log.sh PORT-24-step1 "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_port24_esr_circuit.py tests/unit/test_ring_mode_rows.py tests/unit/test_ans8_private_comparison.py tests/unit/test_port_circuit_layer.py -v -s --tb=short'"`

**Scope boundary.** Closes step 1 only ⇒ `PORT-24` 🟡. Not claimed: any solve, any in-model sheet (item 309), a measured ESR, frequency dependence, any AED re-run; lossless records are never re-read with r.

**Negative result.** Report and stop; park on `attempt/PORT-24-step1`, journal in attempts.md. A red (S)-Im-bitwise is a §7 annotation, never a loosened comparison (see Q1).

**Status it can move.** `PORT-24` ⬜ → 🟡 (the row's step-1 Done-when, `PROJECT_PLAN.md:2716`); unblocks T26 / `PORT-24` step 2.

PREDICTED RED: (Im-bitwise) range [1e-18, 1e-13]

(Max relative |Im λ_k(r) − Im λ_k(0)| / |Im λ_k(0)| over k and records, and the same for `C_mode`. A red inside the range with every other assertion green parks the item as a record for the weekly. A red above 1e-13 stops the slot. The default-bitwise clause (B) is a hard assertion, not covered by this line.)

**Questions for the review (answered below).** (Q1) The Done-when requires Im λ_k *bitwise* while Re shifts at 1e-12. If r enters as `Z + diag(r)` and λ is re-computed by `np.linalg.eig`, Im is not guaranteed bitwise (LAPACK is not shift-invariant to the last bit); if λ(r) is formed as λ(0) + r, Im is bitwise but (S) is tautological for uniform r and only (P) exercises the eig path. The review rules which reading of "λ_k(r)" the step implements (network λ unchanged plus a new loop eigen-impedance key, or Z + diag r), or pre-registers a `PREDICTED RED` for Im-bitwise; the implementer does not choose. (Q2) On the F-small 4×4, λ_k is undefined in the code today (`ring_mode_rows` is 32-port only): the review names it (4-point DFT diagonal of Z, or `eigvals` of Z + rI sorted by k).

**Review rulings (2026-10-08 03:00 daily review).** (Q1) λ_k(r) is computed through the production path on the terminated matrix (Z + diag(r) through the new `esr_ohm` keyword), never formed as λ(0) + r. "Im λ_k and `C_mode` bitwise unchanged" is asserted as written, with the `PREDICTED RED: (Im-bitwise)` line above. The default (r = 0 / unset) bitwise-vs-stored-records clause stays a hard assertion via the r = 0 branch, compared against a pinned pre-change commit, never `HEAD:`. (Q2) On the F-small 4×4, λ_k is the N = 4 DFT diagonal v_kᴴ Z v_k with the normalisation `ring_mode_rows` uses at N = 32. `ring_mode_rows` cannot take N = 4: `N_SLOTS = 16` and two rings are fixed at module level (`ring_modes.py:38–39`), and `_canonical_permutation` (`:56`) rejects anything but 32 ports. So the test module carries a **local helper**, first asserted equal to `ring_mode_rows`' construction (`_dft_matrix` `:68` / `k_blocks` `:74`) on the F-human 32×32 to 1e-14. (Q3) Settled in item 309 below.

### §9 item 309 — 2026-10-08 (chain T26 prerequisite — `PORT-24` step 2, in-model series-RC sheet)

**Outcome (2026-10-08 12:00 implementer slot): GREEN on all three arms ⇒ `PORT-24` ⬜ → 🟡 on step 2 alone (step 1, item 308, stays BLOCKED).** Module `tests/validation/test_port24_series_rc_sheet.py` (new; helpers imported by name, no edit to any other file). Fixture: F-small exact-rectangle hole mesh (`_exact_geo`, `GATED_WIDTH_FRACTION`), 80 798 cells, 64 MHz, degree 1, κ off (`width_correction_kappa=None` asserted on every spec), copper Leontovich column, C = `C_F` = 6.629e-11 F, Z = r − j37.51389328 Ω. Construction per arm: the 50 Ω 4×4 (4 solves), then for each sheet i the field-solved 3×3 with sheet i carrying its series-RC against `reduce_terminated_ports(S₄, 50, {i: series_rlc_impedance(f, r_ohm=r_ref, l_h=0, c_f=C)})` — `PORT-14`'s identity on every sheet — and one P1 drive (P2–P4 in the arm's series-RC) for the power identity.
- **r02** (`20261008T170256Z_PORT-24-step2-r02.log`, 94 s): (I) residuals P1–P4 3.049e-5 / 3.036e-5 / 3.041e-5 / 3.004e-5 (≤ `REDUCTION_BAND` 1e-3, asserted); (P) power identity 1.011e-12 (≤ `DISCRETE_IDENTITY_RTOL` 1e-6, asserted); 50 Ω 4×4 reciprocity 7.4e-16.
- **r0** (`20261008T170618Z_PORT-24-step2-r0.log`, 140 s): (I) 3.097e-5 / 3.084e-5 / 3.088e-5 / 3.051e-5; (Z) r = 0 series-RC sheet vs the capacitor-only sheet (`ans:7` `z_cap`) per-sheet 3×3: 5.2e-16 / 7.6e-16 / 2.5e-16 / 2.8e-16 (≤ `LEG_D0_REPRODUCTION_BAND` 1e-9, asserted; not bitwise 0 — printed as measured); (P) 1.017e-12; RC sheets dissipate exactly 0 W.
- **plant** (`20261008T170439Z_PORT-24-step2-plant.log`, 91 s): (Pl) the planted P2 (0.4 Ω) against the 0.2 Ω reduction **1.095e-2 = 10.95× `REDUCTION_BAND`** (asserted > band); unplanted max 3.049e-5 (asserted ≤ band); (P) 1.012e-12. **The brief's predicted ceiling (≈ 5.3e-3, the terminal-impedance level) is exceeded by ≈ 2×** — the S-level residual is larger than the impedance-level excess (|Γ| ≈ 1 termination, the (I − S_bb Γ)⁻¹ factor). The asserted direction holds; the ceiling was printed, never asserted — recorded as a missed prediction, not a failure.
- **(Dz), printed only:** with the RC sheets' dissipation omitted the identity residual reads 7.031e-3 of P_src at r = 0.2 (1.0e-12 at r = 0). The RC sheets take 70.8 % of the non-port dissipation (P_RC 2.545e-5 W vs P_vol + P_surf 1.051e-5 W; P1's 50 Ω sheet dominates P_src). The brief's "≈ Q_Cu/(Q_Cu + Q_C)" is a share of the dissipated power excluding the 50 Ω load, not of P_src, so the two printed numbers do not compare directly.
- **Q(r)/Q(0), printed only:** at fixed stored energy, P_rest/(P_rest + P_RC) = 0.2923; with Q_C = 1/(ωCr) = 187.57 this implies Q_rest = 454.1 for the copper + phantom losses. P_rest itself moves 0.7 % between r = 0 (1.058e-5 W) and r = 0.2 (1.051e-5 W) on the same drive.
- Price: 23 s mesh, ≈ 5.7 s per 50 Ω 4×4, ≈ 14 s per 3×3, 4.5 s per field solve at `-n 2` — every arm far under the 590 s heavy ceiling; the weekly may re-label the §7 tier column from this.

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* **If item 308 did not land, skip.** Licence: §7 row step 2 (`PROJECT_PLAN.md:2716`), §10 T26 (`:4586–4592`). **No band, tolerance, Done-when or target is set or moved.** Target (row's): identities on one fixture.

- **Executor:** implementer. **Tier:** standard, `-n 2`, one arm per window (r = 0.2; r = 0; planted), container `timeout -k 30 180`. **Class:** chain step (solve). **Depends on:** item 308.
- **Cost:** ≈ 60 predicted min (15 fixed + 20 code + ≈ 12 compute + spare). Prior measurement: `ANS-7-item195-W2-routeA`, 548 s (`docs/testing/test-results.md:2961`). That window is the in-model exact-strip Route A on F-small, but it covers several sweep nodes. A single 64 MHz point should be well under it, and T26 calls the class "257 s" (`PROJECT_PLAN.md:4592`). **If a single-arm window overruns 180 s, kill and report; never raise the timeout.** The review may re-tier to heavy, since the §7 row's "standard" was not priced against this.
- **Settled so far:** the sheet impedance already accepts r — `series_rlc_impedance(f, r_ohm=r, l_h=0, c_f=C)` (`src/fem_em_solver/ports/lumped.py:138`); capacitor sheets are placed via `LumpedPortSheet` with `z_term = series_rlc_impedance(…)` (`tests/validation/test_port_lumped_rlc_termination.py:247`). The F-small exact-rectangle fixture with in-model C is `tests/validation/test_ans7_exact_strip.py` (`_exact_geo` `:409`, `route_a` `:425`, `C_F` `:91`); κ off (`width_correction_kappa=None`, §9 binding ruling `:3733`).
- **Anchors (asserted; imported by name).** The sheet (r, C) is compared with `series_rlc_impedance(0.2, 0, C)` at `REDUCTION_BAND` (`tests/validation/test_port_lumped_rlc_termination.py:88`). Its source comment (`:85`) names it "the `PORT-14` §7 entry's pre-stated band for the reduction identity", ‖S′_meas − S′_pred‖_F / ‖S′_pred‖_F (`:15`), so it is kept. The r = 0 vs capacitor-only 4×4 at `LEG_D0_REPRODUCTION_BAND` (`tests/validation/test_port_birdcage_four_port.py:165`); exact power identity with the sheet's own dissipation at `DISCRETE_IDENTITY_RTOL` (`tests/validation/test_birdcage_power_identity.py:227`). Printed: Q(r)/Q(0) vs 1/Q = 1/Q_Cu + 1/Q_C, Q_C = 1/(ωCr).
- **(Pl) planted non-uniform r — predicted, asserted to fail on the planted sheet.** One sheet is at 0.4 Ω, the other three at 0.2. That sheet's identity against `series_rlc_impedance(0.2, 0, C)` must fail. At `C_F` = 6.629e-11 F (`ans7 metrics.json` `c_f`; `test_ans7_exact_strip.py:91`) and 64 MHz, |1/(ωC)| ≈ 37.5 Ω, so the 0.2 Ω excess is ≈ 5.3e-3 relative at the terminal-impedance level. That is ≈ 5× `REDUCTION_BAND`, with a ceiling of ≈ 5.3e-3. The S-level residual is not measured and may be smaller, so the margin is predicted. A planted residual inside `REDUCTION_BAND` means the identity cannot see r at this C: report and stop, never a pass. The three unplanted sheets must still pass.
- **(Dz) dissipation omitted — predicted, printed, never asserted.** With the sheet's own dissipation left out, the power-identity residual should read ≈ the ESR share of loss, Q_Cu/(Q_Cu + Q_C), with Q_C = 1/(ωC r) ≈ 188 at this C. No prior measurement on the same fixture, so it is not asserted.
- **Traps:** complex build first (`tests/environment/test_complex_mode.py`, `FEM_EM_REQUIRE_COMPLEX=1`); exact rectangles only, never a sub-sheet out of a box; no `ufl.max_value`/ordering on complex operands; pin facet-integral quadrature degree; container `timeout -k 30 ≤ 590`.
- **Harness (per arm, `FEM_EM_PORT24_ARM` ∈ {`r02`, `r0`, `plant`}; one test dispatches on it, never `-k`):** `scripts/testing/run_and_log.sh PORT-24-step2-<arm> "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_PORT24_ARM=<arm> PYTHONPATH=/workspace/src timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_port24_series_rc_sheet.py -v -s --tb=short'"`. New module only; helpers imported from `test_ans7_exact_strip.py` / `test_port_lumped_rlc_termination.py` by name, with no edit to either file.
- **Scope / status:** Done-when step 2 ⇒ `PORT-24` ✅ (with step 1 green). Not claimed: measured ESR, frequency dependence, AED. Negative result: report and stop, `attempt/PORT-24-step2`.

**Review rulings, second pass (2026-10-08 03:00 daily review, items 308 / 309).** (309-tier) **Heavy, not standard**: no 64 MHz degree-1 exact four-port window at `-n 2` has been measured on this fixture (the nearest, `ANS-7-item195-W2-routeA`, took 548 s over several points), so each of the three arm windows runs at container `timeout -k 30 590`; a kill parks the item and is the price record. The §7 row's tier column is the weekly's to re-label from the measurement. (308-scope) **Item 308 edits `src/` and adds one new test module only.** It does not edit `examples/ansys_benchmarks/…/07_tuned_copper_birdcage_driven_sweep.py` (that is T27's, `example-runner`, serial on T26) nor `tuned_input` in `test_port_circuit_layer_field.py`: the new module re-implements those terminations with `esr_ohm`, and asserts its r = 0 path bitwise equal to the existing helpers by import. Any `src/` public signature change re-runs `scripts/docs/generate_api_reference.py` and `tests/package/test_wf17_api_reference.py` in the same slot.

**Review rulings (2026-10-08 10:30 daily review, after item 308 parked).** (308-disp) Item 308 stays BLOCKED for the 10-10 weekly. Its red is the Done-when's own clause ("Im λ_k and `C_mode` bitwise unchanged"), measured at 2.889e-15 / 4.571e-16 on F-human through `np.linalg.eig` of the shifted 2×2 k-blocks and exactly 0 on F-small (`20261008T094657Z_PORT-24-step1.log`). Landing the branch as written would put a red test on `main`, and changing the clause to a tolerance is moving a Done-when — neither is a daily's to do. The weekly re-rules the clause (e.g. "bitwise on the F-small DFT diagonal; at the eig round-off floor on F-human") or lands the branch as a record. (309-dep) **Item 309 no longer depends on item 308.** The 03:00 dependency was an ordering, not a code need: 309's sheet uses `series_rlc_impedance(f, r_ohm=r, l_h=0, c_f=C)` (`src/fem_em_solver/ports/lumped.py:138`), which is already on `main`, and none of 308's `esr_ohm` keywords. §10 makes T26 (item 311) serial on T25 **and** `PORT-24` step 2, but it does not make step 2 serial on T25. The §7 row's ✅ needs both steps, so 309 green holds `PORT-24` at ⬜ → 🟡 on step 2 alone, not ✅. Item 311 stays blocked on 308.

**Review ruling 2026-10-11 (recording the 2026-10-10 weekly, §10 rulings "T25 (`PORT-24` step 1, item 308)" and "Ratified … the `PORT-24` 'Q_Cu' naming").** (a) The step-1 clause "Im λ_k and `C_mode` bitwise unchanged" is **re-stated, not loosened**: `np.linalg.eig(A + rI)` has no bitwise contract with `eig(A)`, and the identity is Im λ(r) = Im λ(0) in exact arithmetic. The Done-when now reads "within 1e-12 relative", the Re clause's own band. The planted-control and Re clauses stand. (b) **T25′** lands the parked branch with that one edit. The 2.889e-15 reading is asserted under the re-stated clause and the `PREDICTED RED` line is removed ⇒ `PORT-24` ✅, since step 2 is green (item 309). (c) The step-2 formula's "Q_Cu" is renamed **Q_0**, the lossless-capacitor coil's loaded Q, because T26's formula reads 1/Q = 1/Q_0 + 1/Q_C. Item 309's landed text above keeps its words as history; the row and the live briefs carry Q_0. (d) Item 308 is superseded by item 363 and keeps its BLOCKED record. `attempt/PORT-24-step1-20261008T0948Z` stays until item 363 lands; the next review then deletes it (weekly flag 6). The slot does not delete it.

### §9 item 363 — 2026-10-11 (T25′ — `PORT-24` step 1 lands under the re-stated clause; lane-1 revalidation)

**Outcome (2026-10-11 04:30 implementer slot): GREEN ⇒ `PORT-24` ✅.** Code taken by path from `attempt/PORT-24-step1-20261008T0948Z` (`ring_modes.py` byte-identical to the branch; the test module carries the one re-statement edit plus the printed (P-Im) line). W1 `docs/testing/logs/20261011T093040Z_PORT-24-step1-T25p.log`, smoke, `-n 1`, 43 passed, elapsed 10 s (incl. `test_wf17_api_reference.py`, green). Readings identical to item 308: (Q2) 0.000e+00; (B) bitwise vs `b2174873` on F-human and 12 F-small pairs; (S)-Re 1.724e-13; **(S)-Im 2.889e-15, `C_mode` 4.571e-16, asserted ≤ 1e-12**; (P) 3.125e-02. **(P-Im), printed:** F-human plant max_k \|ΔIm\|/\|Im\| 8.122e-07, \|C_mode(plant)/C_mode(0) − 1\| 2.035e-09 — both ≫ 1e-12, so the re-stated Im comparator does see a non-uniform r on F-human (prediction held). (T) printed 1.0. The parked branch is left for the review to delete (weekly flag 6).

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain T, step T25′ (grep "T25′. `PORT-24` step 1 — land the parked branch"), quoted, not changed: "`esr_ohm` default bitwise against every stored record; at r = 0.2 Ω, Re λ_k(r) − Re λ_k(0) = r at 1e-12 relative, Im λ_k and `C_mode` within 1e-12 relative, the planted 0.4 Ω gap breaks the shift ⇒ `PORT-24` ✅. Lane 1, smoke, no solve. Serial on nothing." **No band, tolerance, Done-when or target is set or moved here.** The Done-when re-statement is the weekly's (ruling above), and the 1e-12 is the existing `SHIFT_RTOL`.

- **Executor:** implementer. **Class:** lane-1 revalidation (no solve). **Tier:** smoke, `-n 1` (pure NumPy; no MPI reduction to hide, as item 308). **Cost:** ≈ 6 s compute, measured: item 308's W1 ran the same module set in 6 s (`20261008T094657Z_PORT-24-step1.log`, footer `Elapsed (s): 6`; `docs/testing/test-results.md` row "PORT-24-step1"). Predicted ≈ 20 min: 10 fixed, 5 for the land and edit, 1 window, 4 for the records. **Depends on nothing.**
- **Inputs, checked on the host 2026-10-11:**
  - `logs/ans8/s32.npz` (gitignored, 19 802 B, present);
  - the tracked `examples/ansys_benchmarks/ans4_birdcage_four_port_10_64_128MHz/metrics.json` and `ans6_copper_birdcage_four_port_10_64_128MHz/metrics.json`;
  - the parked branch `attempt/PORT-24-step1-20261008T0948Z` at `6c959e41`, a single commit on merge-base `b2174873`.
  - If `s32.npz` is absent, stop and mark BLOCKED. Never re-solve it.

**What the branch carries (read 2026-10-11).** `git show --stat 6c959e41` and `git diff --stat main...attempt/PORT-24-step1-20261008T0948Z` both list exactly two files and 363 insertions / 9 deletions:
- `src/fem_em_solver/ports/ring_modes.py`: +69 / −9, the `esr_ohm` keyword on `terminated_currents`, `ring_mode_rows` and `ring_mode_rows_v2`, plus the helpers `_esr_canonical` and `_gap_admittance`;
- `tests/unit/test_port24_esr_circuit.py`: new, 303 lines.

`git log b2174873..main -- src/fem_em_solver/ports/ring_modes.py` is **empty**: `main` has not touched the file since the branch point, so a path checkout of it is exact. The branch carries **no** record file. Its log and two test-results rows were committed to `main` on 10-08 (`docs/testing/test-results.md`, rows "PORT-24-step1" and "PORT-24-step1-wf17-baseline").

**Steps.**
1. **Take code by path, never merge or check out a record:** `git checkout attempt/PORT-24-step1-20261008T0948Z -- src/fem_em_solver/ports/ring_modes.py tests/unit/test_port24_esr_circuit.py`. Nothing else comes from the branch.
2. **The one edit, to `tests/unit/test_port24_esr_circuit.py` only. `ring_modes.py` stays byte-identical to the branch.**
   - Module docstring, the (S) bullet (lines 17–20 on the branch): replace "and ``Im λ_k``, ``C_mode`` bitwise unchanged, carrying the brief's ``PREDICTED RED: (Im-bitwise) range [1e-18, 1e-13]``" with "and ``Im λ_k``, ``C_mode`` unchanged within 1e-12 relative (``SHIFT_RTOL``). The Done-when was re-stated from 'bitwise' by the 2026-10-10 weekly (T25): ``eig(A + rI)`` has no bitwise contract with ``eig(A)``. Item 308 measured 2.889e-15 / 4.571e-16."
   - Rename `test_uniform_shift_im_and_c_mode_bitwise` (`:267`) to `test_uniform_shift_im_and_c_mode_unchanged`, and delete its `# PREDICTED RED:` comment (`:268`).
   - Its assertion becomes `assert im_worst <= SHIFT_RTOL and c_worst <= SHIFT_RTOL, (im_worst, c_worst)`.
   - Its summary print ends "(ASSERTED <= 1e-12 relative, SHIFT_RTOL; item 308 read 2.889e-15 / 4.571e-16)" in place of "(asserted bitwise: 0)".
   - Keep the per-record `[PORT-24 S-Im]` lines as they are.
   - **(P-Im), printed only (see Negative controls).** In `test_planted_nonuniform_breaks_shift`, after the F-human loop, print `max_k |Im λ_k(plant) − Im λ_k(0)| / |Im λ_k(0)|` and `|C_mode(plant)/C_mode(0) − 1|` on F-human, from `r0` / `rp` (already computed). Add no assertion.
   - `SHIFT_RTOL` (`:60`, "the Done-when's 1e-12 relative") is reused, not re-declared. Add no new constant.
3. **API pages.** `ring_modes` has no page in the WF-17 API reference: item 308's outcome says so, and `scripts/docs/generate_api_reference.py` names no `ring_mode*`. So nothing is regenerated. `tests/package/test_wf17_api_reference.py` runs in W1 as the regression. It was red on 10-08 for an unrelated stale page, and that entry is retired: known-issues, "`WF-17` API pages, §9 item 331".
4. **W1:** item 308's W1 command, with the label changed and `tests/package/test_wf17_api_reference.py` appended, as in the 10-08 log's command line. The Bash call's `timeout` is 660000. The harness call is the whole Bash command: bare and repo-relative, with nothing chained or piped.

   `scripts/testing/run_and_log.sh PORT-24-step1-T25p "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_port24_esr_circuit.py tests/unit/test_ring_mode_rows.py tests/unit/test_ans8_private_comparison.py tests/unit/test_port_circuit_layer.py tests/package/test_wf17_api_reference.py -v -s --tb=short'"`

   The `timeout -k 30 120` is item 308's reviewed ceiling. The measured 6 s is 20× under it, so a kill means a wedged container (recover per known-issues), not a price.
5. **Commit** the two code files, the W1 log, its test-results row, the §7 row edit and the §9 marking together. Before committing, run `git diff --stat --cached`. It must show `test-results.md` with insertions only, and no `docs/testing/logs/*PORT-24-step1.log` (10-08) modification.

**Anchors (asserted; every band imported or already in the module, none new):**
- **(B) default bitwise.** With `esr_ohm` unset, `ring_mode_rows` / `_v2` / `terminated_currents` on F-human's stored 32×32 equal the pinned pre-change source (`PRE_CHANGE_SHA = "b2174873"`, loaded via `git show`, never `HEAD:`) at `SPEC_C_F` and `C2_RECORD_F`. The F-small `tuned_input` r = 0 path is bitwise on 12 (record, termination) pairs. `test_ans8_private_comparison.py` is green.
- **(Q2)** The local DFT helper matches `k_blocks` to ≤ `HELPER_RTOL` (1e-14). Item 308 read 0.000e+00.
- **(S)-Re** max \|ΔRe λ − r\|/r ≤ `SHIFT_RTOL` over 33 (record, k). Item 308 read 1.724e-13 (`20261008T094657Z_PORT-24-step1.log:96`).
- **(S)-Im, re-stated.** Max \|ΔIm λ_k\|/\|Im λ_k(0)\| ≤ 1e-12 and \|ΔC_mode\|/C_mode ≤ 1e-12. Item 308 read 2.889e-15 / 4.571e-16 (`20261008T094657Z_PORT-24-step1.log:107`); F-small is exactly 0 (`:101–106`). **Predicted:** the same digits, since the module path is unchanged; 346× and 2 188× under the bar.
- **(P) planted.** Min \|ΔRe − r\|/r over every mode > `SHIFT_RTOL`. Item 308 read 3.125e-02 (`:142`).

**Negative controls (rule (e)).**
- **(P): asserted.** It is backed by item 308's reading on the same fixture and comparison: 3.125e-02 on F-human, 0.25 on F-small (`20261008T094657Z_PORT-24-step1.log:109–142`), ≈ 3e10 above the bar. Ceiling: the full plant excess, 0.2 Ω.
- **(P-Im): predicted, printed, never asserted.** It asks whether the re-stated Im comparator can see anything at all.
  - On F-small, a one-gap plant shifts every N = 4 DFT diagonal entry by the same real mean(r). Im is therefore predicted **exactly 0**, so it cannot serve as a control there.
  - On F-human, the plant at P17 shifts only the bottom ring's k-diagonal (by 0.2/16 = 0.0125 Ω beyond the uniform 0.2 Ω). The 2×2 k-block becomes unequal on its diagonal, and its antisymmetric eigenvalue's Im moves at second order. That is predicted ≫ 1e-12, with magnitude unknown (no prior measurement).
  - A reading ≤ 1e-12 on F-human would mean the Im clause is blind to non-uniform r on this fixture. That is a §7 annotation for the review, not a red.
- **(T) termination-only wiring: predicted, printed** as on the branch, off by exactly 1.0.

**Traps already paid for.**
- **Path checkout of append-only records.** `git checkout attempt/… -- docs/testing/test-results.md` replaces the file and silently deletes every row `main` appended since the branch. `8d4cf58` (`TH-15` step 3b, 2026-09-09) dropped nine rows across four chunks. Take code by path (step 1 names exactly two files), append record rows by hand, and check that `git diff --stat` shows only insertions on record files.
- The pinned pre-change commit, never `HEAD:` (`test_orphan_guard.sh` / `OPS-43`). The control would go stale the moment this change lands.
- The r = 0 path is a branch, not arithmetic: `1/(0 + 1/(jωC))` is not bitwise `jωC`. It is already on the branch; do not "simplify" it.
- pytest `-s`. Never `-k a or b` inside the quoted command. Harness foreground. Never pipe pytest inside the harness command.
- Rule (i): no edit to the module after W1 starts. If one is needed, re-run W1.

**Scope boundary.** Closes `PORT-24` step 1, and with item 309 the row ⇒ ✅. **Not claimed:** any solve; any edit to `tuned_input`, `ans:7`'s terminations or `route_b_z` (those are T26's test-side and T27's example work, the 308-scope ruling); a measured ESR; frequency dependence; any AED re-run. Lossless records are never re-read with r.

**What a negative result means.**
- **(B), (Q2), (S)-Re or (P) red:** this contradicts item 308's green on an unchanged `src/` file. Report the reading and stop. Do not commit the land; leave the branch where it is, and mark item 363 BLOCKED in the same commit as the attempts.md journal (rule (d)).
- **(S)-Im above 1e-12:** the measurement moved by ≥ 346× on unchanged code. Report and stop, BLOCKED; never widen.
- **`test_wf17_api_reference.py` red:** a pre-existing failure unrelated to this item. Check known-issues first. If it is not listed, land item 363 on its own green tests and open a known-issues entry; do not edit the API pages in this slot.
- No `PREDICTED RED:` line.

**Status it can move.** `PORT-24` 🟡 → ✅, by rule (j) against the row's Done-when **as amended by the review commit that queues this item**. The amended row carries "Im λ_k and `C_mode` are unchanged within 1e-12 relative" and "Q_0". If the row still reads "bitwise unchanged" when the slot opens, the slot closes on the row: it holds 🟡, lands the code as a record and flags the mismatch (rule (j)). In the landing commit, the implementer:
- flips the glyph 🟡 → ✅;
- replaces the row's last bold sentence ("**Step 1: T25′ (§9 item 363) lands the parked branch under the re-stated clause.**") with "**✅ 2026-10-11 (item 363, T25′): step 1 landed; (S)-Im ≤ `<read>` / `C_mode` ≤ `<read>` against 1e-12; `<W1 log>`.**";
- marks §9 item 363 DONE, and item 308 "superseded by 363".

The closure then goes to the daily's auditor. T26 (item 311) becomes runnable.

### Audit 2026-10-11 10:30 review — ✅ PASS

`auditor` on `9b90d7b3` (step 1) and `91254c5a` (step 2): PASS.

**Step 1** (`20261011T093040Z_PORT-24-step1-T25p.log`):
- Status 0 after 10 s (`:261–263`), 43 passed (`:259`).
- (S)-Im 2.889e-15 and `C_mode` 4.571e-16 (`:100`, `:107`) against 1e-12, asserted at `tests/unit/test_port24_esr_circuit.py:277`.
- Re λ shift = r (`:65`, `:96`; assert `:255`).
- Planted 0.4 Ω control (`:110–117`, `:143`; assert `:312`).
- Unset bitwise (B) at `:185–210`.

**Step 2** (`20261008T170256Z_PORT-24-step2-r02.log:940–947`, `…170439Z_…plant.log:942–949`, `…170618Z_…r0.log:941–945`):
- Status 0 at 94 / 91 / 140 s.
- Asserts at `tests/validation/test_port24_series_rc_sheet.py:261–272`.

The red baseline is item 308's bitwise clause (`20261008T094657Z_PORT-24-step1.log:107–108`). The code is unchanged; the bar was re-stated to 1e-12 by the 2026-10-10 weekly (T25), with a measured basis (346× / 2 188× above the readings).

**Caveats (no demotion):**
- **"Circuit-layer terminations" are test-local on F-small** (the 308-scope ruling). `esr_ohm` exists in `src/` only on `ring_mode_rows`, `ring_mode_rows_v2` and `terminated_currents`. F-small's `tuned_input_esr` / `capacitor_impedance_esr` live in the test module, and `tuned_input` and `ans:7`'s `route_b_z` are deliberately unedited (T26 / T27).
- F-small's Im clause is exact by construction (the N = 4 DFT diagonal of Z + diag(r) shifts only Re); only F-human exercises the `eig` path.
- Step 2's "r = 0 vs capacitor-only" is per sheet (3×3, ≤ 7.6e-16), not the 4×4.
- Step 2's Q(r)/Q_0 0.2923 is printed at fixed stored energy on a P1-driven solve, not a resonance Q. Item 311's resonance-Q ratio is 0.1694 against the closed form 0.1676.
- The planted P2 reads 1.095e-2, above the brief's predicted 5.3e-3 ceiling, which was never asserted.

The parked branch `attempt/PORT-24-step1-20261008T0948Z` (`6c959e41`) is deleted: `ring_modes.py` on main is byte-identical to it.
