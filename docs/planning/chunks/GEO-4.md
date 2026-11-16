# GEO-4 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 185 — 2026-10-02 (lane 1)

**Outcome (2026-10-03 22:30 CDT slot): DONE — W1 29 passed, W2 20 passed / 9 skipped; `GEO-4` 🧪 → ✅.**

- **Code (rule (c), additive, disclosed):** `src/fem_em_solver/io/mesh.py` +113 lines — `PORT1_PADDING_RUNGS_M` / `_DEFICITS_PP`, `_free_exponent_fit` (verbatim from the probe), `_padding_for_box_term` (the ceil with the 1e-12 guard, shared by C2), and `MeshGenerator.two_torus_box_padding`. Its pre-existing gate `tests/mesh/test_domain_sizing_heuristics.py` is green in both windows.
- **Windows:** smoke `20261003T034956Z_GEO-4-smoke.log` (real, Status 1, before W1; see the disclosure below), W1 `20261003T035018Z_GEO-4-W1.log` (complex `-n 2`, 226 s, predicted ≈ 190 s), W2 `20261003T035410Z_GEO-4-W2.log` (real `-n 1`, 55 s). No edit after W1.
- **(A0):** D∞ 1.6934077797495957 pp, C −0.1478719225310089, p 1.6574194246323086, `==` the probe. **(A1):** τ = 6 → 0.11 m (unrounded 0.1070673, \|B\| 5.737202; \|B(0.10)\| 6.719 > 6); τ = 4 → 0.14 m (0.136742, 3.846889; \|B(0.13)\| 4.350 > 4).
- **Disclosure, A1 approximate values:** the executor first pinned the brief's "≈ 0.10706 / B 3.8467" to ½ unit of the last printed digit; the smoke window measured 0.1070673 and 3.846889 (the brief's digits were truncated / rounded), and these two informational checks were re-pinned to 1e-4 relative before W1. The brief's asserted A1 clauses (predicted(padding) ≤ τ < predicted(padding − q); predicted(unrounded) = τ to 1e-12) are unchanged and green.
- **(P0):** −8.0404 pp at 0.08 m vs 3b-xi −8.0324 (\|Δ\| 0.0080 ≤ 0.05; 119 533 cells vs 119 738).
- **(G):** −4.0256 pp at 0.11 m (145 101 cells), −2.2071 pp at 0.14 m (175 228 cells); \|deficit − D∞\| 5.7190 / 3.9005 pp, margins +0.2810 / +0.0995 pp (fit predicted +0.2628 / +0.1531). The p = 3 fit predicted −3.9912 / −2.6718 pp and would fail τ = 4 (4.365 > 4); the free exponent is 0.054 pp off at 0.14 m.
- **(C1):** box term at 0.08 m 9.7258 pp (fit) / 9.7338 pp (measured) = 1.62× τ = 6, 2.43× τ = 4. **(C2):** p = 3 → 0.10 m for τ = 4, recorded 6.7190 pp > 4 (asserted); τ = 6 → 0.09 m at 8.001 pp (printed).
- **Record for `PORT-1`:** −2.2071 pp at 0.14 m, a fourth 3b-xi rung (factor 1.75).

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.*

**`GEO-4` step 2 — the box-term padding heuristic on `PORT-1`'s fit.**
- **Executor:** implementer. **Tier:** heavy, complex, `-n 2`.
- **Cost:** ≈ 4 measured / ≈ 40 predicted min.
- **Depends on:** item 182 having landed or been marked BLOCKED (file order only: 182 takes `io/mesh.py` by path from its branch, and this item adds to it). **Attempt:** 1 of 4.

**What changed since the rejection.** The 09-29 18:00 review rejected `GEO-4`'s 0.35·extent rule as a tautological formula. The weekly now names a physical target: `PORT-1`'s decision-(4) fit.
- The fit is deficit(W) = D∞ + C·W^(−p), with D∞ = +1.6934 pp, C = −1.478719e-01 and p = 1.6574 (`20260812T213337Z_PORT-1-dec4-fit.log`, `[fit]` block).
- It was fitted on 3b-xi's rungs: −8.0324 / −5.0256 / −3.2733 pp at W = 0.08 / 0.10 / 0.12 m (`scripts/probes/port1_padding_fit.py:72–73`).
- The box term is B(W) = C·W^(−p).

**Code (additive, rule (c)), in `src/fem_em_solver/io/mesh.py` beside `coil_phantom_domain_sizing_diagnostics`.**
- **Constants:** `PORT1_PADDING_RUNGS_M` and `PORT1_PADDING_DEFICITS_PP`, with provenance.
- **Fit helper:** `_free_exponent_fit`, copied **verbatim** from the probe's `power_law_extrapolation` (`:81–127`).
- **`MeshGenerator.two_torus_box_padding(tolerance_pp, *, quantum_m=0.01)`.**
  - It returns `padding_m` = ceil((|C|/τ)^(1/p) / q)·q, with a 1e-12 guard, plus `unrounded_padding_m`, `predicted_box_term_pp`, `d_inf_pp`, `c_pp` and `p`.
  - It raises `ValueError` for τ ≤ 0 or when there is no fit.
- **Do not hard-code the printed C and p.** Rounded, they reproduce the 0.08 m rung only to 5e-4 pp. Recompute the fit from the rungs.

**Test: new `tests/validation/test_geo4_box_padding_heuristic.py`.**
- **Imports.**
  - The probe, loaded by `importlib` (as `test_post8_owner_rule_ans2.py:107` does).
  - From `test_port_reaction_impedance`: `_solve_reaction_z`, `MAJOR_RADIUS`, `OMEGA`, `SEPARATION`, `TORUS_TAGS` and `mutual_inductance`.
  - From `test_port_box_padding_sweep.py:87–88`: `STEP2F_DEFICIT_AT_008` / `_TOLERANCE`.
- **(A0)** The rung constants `==` the probe's. The fit `==` the probe's output bitwise, and lies within ½ ulp of the printed D∞ / C / p.
- **(A1)** τ = 6 gives 0.11 m (unrounded ≈ 0.10706, predicted B 5.737). τ = 4 gives 0.14 m (≈ 0.13674, B 3.8467).
  - Asserted: predicted(padding) ≤ τ < predicted(padding − q).
  - Asserted: predicted(unrounded) = τ to 1e-12.
- **(P0) Premise: one solve at 0.08 m.**
  - Within the step-2f pin.
  - |deficit − (−8.0324)| ≤ 0.05 pp, a third of the smaller margin.
  - This band is **predicted**: the fixture has no 0.11-image run, and only `20260807T124038Z_PORT-1-step3bxi-gate.log` prints these deficits.
- **(G) The gate: one solve each at 0.11 and 0.14 m** (one driven column, as 3b-xi `:109–117`).
  - Every deficit < 0.
  - |deficit| strictly decreasing over 0.08 → 0.11 → 0.14.
  - |deficit_pp(W_τ) − D∞| ≤ τ, with margins 0.263 / 0.153 pp under the fit.
  - Printed beside the readings: the p = 3 fit's predictions, −3.99 / −2.67 pp. At 0.14 m that fit gives |B| = 4.37 > 4, so **0.14 m tells the two exponents apart**: they differ by 0.52 pp against a rung precision of 5e-5.

**Negative controls.**
- **(C1) Un-sized padding (asserted).** The production padding 0.08 m has box term ≈ 9.726 pp, which is 1.62× τ = 6 and 2.43× τ = 4. Backing: the 3b-xi log `:438` measured the same fixture and comparison, and P0 re-measures it.
- **(C2) The dipolar-exponent heuristic** (p = 3, C = −3.410151e-3, from the probe's diagnostic). No extra solve.
  - **Asserted** for τ = 4: it returns 0.0948 → 0.10 m, whose recorded box term is 6.719 pp > 4 (`:828`).
  - **Predicted, printed** for τ = 6: 0.0828 → 0.09 m, at ≈ 8.0 pp.

**Windows.** Foreground, Bash `timeout: 660000`, each harness call bare and whole.
- **W1** (≈ 190 s, on 3b-xi timings: mesh + solve 43 / ≈ 52 / ≈ 70 s; ≈ 145 k and ≈ 177 k cells, under the 250 k MUMPS line):

  `scripts/testing/run_and_log.sh GEO-4-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_domain_sizing_heuristics.py tests/validation/test_geo4_box_padding_heuristic.py -v -s --tb=short'"`
- **W2** (real build, ≈ 15 s; the solve tests skip):

  `scripts/testing/run_and_log.sh GEO-4-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment tests/mesh/test_domain_sizing_heuristics.py tests/validation/test_geo4_box_padding_heuristic.py -v -s --tb=short'"`

**Traps:**
- **Signed deficits.** The deficits are signed pp; never compose relative percents.
- **Sign of C.** C is negative in the commissioning convention; the gate uses |B|.
- **Solve tests.** Mark them `complex_only` and share them through one module-scoped fixture.
- **A kill at 590 s is not a result.** Run `pgrep -c python3` first, and never lengthen the timeout.
- **Run discipline.** Run with `-s`. Make no edit after W1.

**Scope:**
- The heuristic is two-torus only, an effective-range fit, not a converged box-free deficit.
- The coil+phantom 0.35 rule stays as it is.
- `PORT-1`'s ✅ is untouched. The 0.14 m reading is a fourth-rung **record** for `PORT-1` (factor 1.75, short of the named factor ≥ 2).

**Status it can move:** `GEO-4` 🧪 → ✅ on the Done-when this review wrote into its row (rule (j)). The τ = 4 leg may go red honestly: 0.14 m is outside the fit's range, and the p = 3 fit predicts a miss there. That would be an informative negative.

**Negative result:**
- **P0 off by more than 0.05 pp:** the fit's data moved on 0.11. Stop before the gate; re-registering the fit is the review's job under (1*).
- **(G) red at τ = 4 with τ = 6 green:** the fit does not extrapolate to 0.14 m. Record the rung under `PORT-1` 3b-xi in §7, park, and mark BLOCKED; the weekly re-poses.
- **τ = 6 red:** stop.
- There is no PREDICTED RED line.

**Audit, 2026-10-03 03:00 daily review (§4): PASS.** The anchor is `assert excess <= tau` (`test_geo4_box_padding_heuristic.py:305`), which passed at τ = 6 and at τ = 4 (`20261003T035018Z_GEO-4-W1.log:684–689`). The window ran 226 s with Status 0 (W1:823–824).
- The fit constants trace to `20260812T213337Z_PORT-1-dec4-fit.log.gz:52–59`.
- The smoke red, and the A1 re-pin to 1e-4 relative before W1, are disclosed (`…034956Z_GEO-4-smoke.log:68–80`). The binding A1 asserts are unchanged.
- Caveats: the gate is internally referenced, since D∞ is the fit's own extrapolation. The τ = 4 margin is thin (+0.0995 pp), and 0.14 m lies outside the fit's 0.08–0.12 m range. C2's asserted leg is arithmetic on recorded rungs, not a fresh solve. The row's "3.901" rounds the log's 3.9005.
- Review spot-check: W1:689 re-read, and it matches.
