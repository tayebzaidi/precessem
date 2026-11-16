# TH-2 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 153 — 2026-10-01 (lane 1 revalidation; queued by the 10:30 review)

`TH-2` — the ⚠️ row's functional surface gated through the harness (implementer; standard, `-n 2`, complex; ≈ 3 measured / ≈ 23 predicted min; depends on nothing).

- **Anchor.** Nothing new is written. The row's one functional surface, the per-tag `material_map` DG0 assignment (`core/time_harmonic.py:184–255`), is already gated by quantities in two modules: `OPS-13` — tagged volume = 1/162, ∫σ = 200·V, ∫μ_r and ∫ε_r, all to 1e-12 and reduced on every rank (`tests/materials/test_material_map_rank_safety.py:164–184`); `MAT-1` (P1) DG0 σ / ε_r / μ_r bitwise on 74 020 owned dofs, (P2) ‖A_ph − A_map‖∞ = 0.0, (P3) 3.536 % < 10 %, presets SAR 3.523 / 3.548 / 3.557 %, Im/Re 0.404 / 0.835 / 1.499 % (`20260930T003229Z_MAT-1-reanchor-W1.log:485–499`). **The three raise-before-solve tests (`tests/solver/test_time_harmonic_smoke.py:625–696`) assert no quantity** — they check a `ValueError` and its message; API coverage, no anchor.
- **Negative controls (asserted, already in the modules).** `OPS-13`: tag 4242 rejected on every rank (`:187–213`), and the phantom-tag twins (`:334–380`). `MAT-1` (c1): a 1-ULP σ reads 1.11e-16 on 19 496 dofs; (c2): a 1e-9 σ reads 2.406e-14 relative (`…:491–493`).
- **Windows.**
  - W1 (≈ 140 s; the re-anchor window ran 123 s for env + phantom, 17 passed, plus ≈ 15 s for the cube and three cylinder meshes): `scripts/testing/run_and_log.sh TH-2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/materials/test_material_map_rank_safety.py tests/materials/test_phantom_material_model.py tests/solver/test_time_harmonic_smoke.py::test_time_harmonic_solver_rejects_non_hz_frequency_unit_before_solve tests/solver/test_time_harmonic_smoke.py::test_time_harmonic_solver_rejects_material_map_without_cell_tags_before_solve tests/solver/test_time_harmonic_smoke.py::test_time_harmonic_solver_rejects_unknown_material_map_tag_before_solve -s -v --tb=short'"` — predicted 26 passed, 0 skipped, the `…reanchor-W1.log:485–499` readings reproduced (the solve-level line printed only; it drifts ≤ 4.45e-11, `OPS-43`).
  - W2 (record, real build, ≈ 40 s): label `TH-2-W2-real`, the same paths without `source …complex-mode` and without `FEM_EM_REQUIRE_COMPLEX`, `timeout -k 30 200` — predicted: the two `@complex_only` phantom tests skip and the three raise tests pass (validation runs before `require_complex_mode()`, `:483–485`).
- **Traps.** Select tests by node IDs, never `-k a or b`. `PYTHONPATH` must include `/workspace` (the phantom module imports `tests.validation`). `cell_tags.values` is rank-local — add no local reads. Edit no module (rule (i)). Harness bare, foreground, `timeout: 660000`.
- **Scope.** No new test, no `src/`. Not claimed: the non-finite, ≤ 0 and > 1e12 frequency guards and the phantom-frequency mismatch guard (`:447–466`) have no test.
- **Status it can move.** `TH-2` ⚠️ → ✅ on the Done-when written by the 10:30 review (the italic "retire or merge" note replaced).
- **Negative result.** `OPS-13` red: the reduction in `_validate_material_map_tags` regressed — known-issues entry, stop, ⚠️ stays. `MAT-1` red with (P1)/(P2) green: report the numbers and stop. W2 red only: known-issues entry; the ✅ stands on W1.

**Outcome (2026-10-02 07:30 slot, second item, run by the slot itself):** DONE, ⚠️ → ✅.
- **W1** (`20261002T124713Z_TH-2-W1.log`): complex `-n 2`, 26 passed, 0 skipped as predicted, 118 s. `OPS-13`'s five tests are green. `MAT-1` (P1) σ / ε_r / μ_r read 0.0 with 0 differing on 74 020 owned dofs, and (P2) ‖A_ph − A_map‖∞ = 0.0. (c1) read 1.1102e-16 on 19 496 dofs and (c2) 2.406e-14 relative. (P3) read 3.536 %, presets SAR 3.523 / 3.548 / 3.557 % and Im/Re 0.404 / 0.835 / 1.499 % (`:499–512`). The solve-level record read 8.620e-12, inside the `OPS-43` drift floor.
- **W2** (`…124920Z_TH-2-W2-real.log`): real build, 20 passed and 6 skipped, 50 s. The brief predicted two skips, the `@complex_only` phantom tests. The other four skips are the `tests/environment` complex-only checks, which skip on a real build by design. The three raise tests passed.


### Audit — 2026-10-02 10:30 daily review (`auditor`, closure `54f9ce89`)

**PASS.**
- **Digits:** 25 / 25 traced. The review re-read W1:509 ((P3) 3.536 %, band 10 %).
- **W2 skips:** all six are complex-only by design: the two phantom tests and four environment tests.
- **Done-when:** asks only for the complex window, which ran 26 passed with 0 skipped (`20261002T124713Z_TH-2-W1.log:1039`).

Caveats:
- W1's container ceiling was `timeout -k 30 400`, above the standard tier. It measured 118 s.
- The frequency and phantom-frequency guards stay untested, as the row says.
