# OPS-82 — the `IMP-8` DG0 source helpers build in both builds

*History file for the §7 `OPS-82` row (`OPS-46`). Opened by the 2026-11-05 18:00 daily review, lane 1, from `EX-182`'s finding (§9 item 829, attempts.md 2026-11-05T19:54Z): `_rod_source` (`tests/validation/test_imp8_rod_line_source_fem.py:152`) and `_tips_source` (`tests/validation/test_imp8_rod_tips_fem.py:74`) assemble `one * w * ufl.dx` with an unconjugated test function, which the complex build rejects (`ArityMismatch`, `20261105T193522Z_EX-182-figure.log:104`). Every gate that calls them runs in the real build, so nothing on `main` is red; `imp:23` carries disclosed local `ufl.conj(w)` copies. Known-issues 2026-11-05 ("the `IMP-8` source helpers are real-build only") files it, and this chunk is its "Resolves with".*

**Done-when.** Both helpers assemble `one * ufl.conj(w) * ufl.dx`; a new unit module gates them in both builds — (V) the owned DG0 volume sum equals the box volume, (S0) ∫q = P and each ball P/2 within `B_BAND` — with the old form's complex-build refusal asserted as the control and the real-build vector bitwise unchanged; no other line of either module edited, no `src/` edit; the known-issues entry retired.

### §9 item 849 — 2026-11-05 (lane 1 / ops; `OPS-82` step 1)

**Edit.** In each helper, `fem.form(one * w * ufl.dx)` → `fem.form(one * ufl.conj(w) * ufl.dx)`. Nothing else in either module changes (`git diff` shows one line per file). `conj` is the identity on a real operand, so every real-build gate that imports them (I55, I56, I60, I63, I66's module once it lands) is unchanged by construction; the bitwise check below makes that a measurement.

**Test.** New `tests/unit/test_ops82_dg0_source_helpers_both_builds.py`:
- Mesh: `dolfinx.mesh.create_box` on [−0.05, 0.05]² × [−0.03, 0.03] m, 12 × 12 × 8 cells (no gmsh), so the module runs in either build in seconds.
- Call `_rod_source(msh, power=1.0, length=0.03, radius=0.012)` and `_tips_source(msh, power=1.0, radius=0.012, length=0.03)` (radii large enough that each selection is non-empty at this resolution; assert the counts > 0, allreduced).
- **(V)** the helper's own `bvol` route — re-built in the test from the same form — sums over owned cells to the box volume 0.1 × 0.1 × 0.06 = 6e-4 m³ within 1e-12 relative (allreduced).
- **(S0)** ∫q dx = P (rod) and each ball's ∫q dx = P/2, `np.real` of `_allsum`, within `B_BAND` (imported from `test_mat10_sar_to_heat_sphere`).
- **Control, complex build:** a local copy of the pre-fix form, `fem.form(one * w * ufl.dx)`, raises — `pytest.raises(Exception)` with `ArityMismatch` in the message; skipped in the real build.
- **Control, real build:** the pre-fix form's assembled vector equals the fixed form's **bitwise** (`np.array_equal`) on the same mesh; skipped in the complex build.

**Windows.** Bare, foreground, Bash timeout 660000 ms.
- W1 real: `scripts/testing/run_and_log.sh OPS-82-W1-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/unit/test_ops82_dg0_source_helpers_both_builds.py -v -s --tb=short'"`
- W2 complex: the same with `source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1` and `tests/environment/test_complex_mode.py` first (`20261105T182800Z_OPS-81-W0.log:12`'s form).
- No I55 / I56 re-run is bought: the real-build bitwise control is the evidence that their numbers cannot move.

**Tier.** Smoke, `-n 2`, both builds. Predicted < 30 s each. **Class:** ops / docs. **Cost:** ≈ 3.4 measured / 15 predicted min.

**Traps already paid for.** In the complex build a raising rank can leave its peer hung to the ceiling (`imp:23`'s W1, 540 s, Status 124): keep the refusal inside `pytest.raises` on every rank (the form is built collectively), and keep the 120 s ceiling. Count verdicts per rank by token. `-s` for prints. Reduce rank-local sums before asserting. Do not edit `examples/implants/23_rod_line_source_fem_twin.py` (its local copies are disclosed and its closure audited; dropping them is a later example edit, not this item).

**Scope.** Two helper lines and one unit module. **Not claimed:** any `IMP-8` record, `imp:23`'s copies, other real-only helpers elsewhere.

**Negative result.** (V) or (S0) red in either build: stop, BLOCKED, revert nothing on `main` (the edit is not committed). The real-build bitwise control red: the fix changed real-build numbers — report, BLOCKED. No `PREDICTED RED:` line.

**Status it can move.** `OPS-82` ✅; known-issues 2026-11-05 (the helpers) retired.

**Depends on:** nothing.

**Outcome (2026-11-06T03:05Z, 21:00 slot): DONE.** The two helper lines now read `one * ufl.conj(w) * ufl.dx` (`git diff` one line per file); `tests/unit/test_ops82_dg0_source_helpers_both_builds.py` is green in both builds at `-n 2` on the 6 912-cell box: (V) owned volume sum 6.000000000000003e-04 vs 6.000000000000001e-04 m³, rel 3.614e-16; (S0) rod n_sel 144, ∫q |d| 7.772e-16; balls n_sel 76 / 76, ∫q |d| 2.220e-16 / 3.331e-16, total 2.220e-16 (`B_BAND` 1e-10); real control: pre-fix and fixed vectors bitwise equal on every rank; complex control: the pre-fix form refused with `ArityMismatch: Failure to conjugate test function in complex Form`, the fixed form passing the same check. Closing windows `20261106T030402Z_OPS-82-W2-complex-r2.log` (8 passed, 1 skipped per rank, 3 s) and `20261106T030411Z_OPS-82-W1-real-r2.log` (4 passed, 1 skipped per rank, 3 s). **Disclosed:** the first complex window `20261106T030126Z_OPS-82-W2-complex.log` ran the control as `pytest.raises(Exception)` around `fem.form(...)`. That control went red and the window hung to its 120 s ceiling (Status 124; `pgrep -c python3` read 0 afterwards). Two causes: `ArityMismatch` derives from `BaseException` (MRO checked in the container), so `pytest.raises(Exception)` does not catch it; and dolfinx's `mpi_jit` compiles on rank 0 only, so the uncaught raise there left rank 1 waiting. The refusal itself was measured in that log (line 119). The control now calls `ufl.algorithms.check_arities.check_form_arity(..., complex_mode=True)` directly — the function that raised at `fem.form`'s line 118 — which is rank-local and needs no JIT. The module was edited after the first windows, so both builds were re-run (rule (i)). **Trap for future briefs:** "keep the refusal inside `pytest.raises` on every rank" does not hold for a UFL arity refusal through `fem.form`; catch `ArityMismatch` by name, at the UFL check.

### Audit 2026-11-06 03:00 daily review — PASS

- **Evidence (`auditor`, re-cited):** `20261106T030402Z_OPS-82-W2-complex-r2.log:63–81` ((V) rel 3.614e-16; rod |d| 7.772e-16; balls ≤ 3.331e-16; control `:75–76`; 8 passed / 1 skipped per rank `:81`, `:84`; Status 0, 3 s `:86–88`) and `20261106T030411Z_OPS-82-W1-real-r2.log:51–75` (bitwise control `:65`; 4 passed / 1 skipped per rank). The first complex window `20261106T030126Z_OPS-82-W2-complex.log` (Status 124 at 120 s) is kept and disclosed.
- **Ruling on the edit and re-run:** compliant with rule (i) and §4. The red was test mechanics (`ArityMismatch` is a `BaseException`), not a numeric bound. The rewritten control asserts the same refusal by name, adds a positive arm, and both builds were re-run.
- **Caveat for readers:** the control calls UFL's `check_form_arity(..., complex_mode=True)` directly. It is not a `fem.form` refusal test; the kept red log is the record of that route.
