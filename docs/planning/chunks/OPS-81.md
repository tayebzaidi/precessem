# OPS-81 — the two-tag power-identity helper states its scope and refuses a lossy tag outside (1, 3)

*History file for the §7 `OPS-81` row (`OPS-46`). Opened by the 2026-11-05 03:00 daily review under §9 rule (5), lane 1, on the 2026-11-04 weekly's ruling (PROJECT_PLAN §10, "The two-tag power-identity helper … not widened"). The known-issues entry "OPEN 2026-11-03 (`GEO-36` I57 …)" retires on this chunk's commit.*

### Review ruling 2026-11-05 03:00 — the guard's branch, and an acknowledgment for the (K) wrappers

**The guard branch applies.** The case object the helper receives carries its material map. `build_four_port_sweep` returns `"problem"`, a `TimeHarmonicProblem` with `.material` and `.material_map`, in both its return dicts (`tests/validation/test_port_birdcage_four_port.py:469–487`, `:521`, `:625`).

**One premise in the weekly's text is corrected, as execution detail; the target is kept.** Every kept-cell (K) caller in `OPS-80`'s table writes a lossy tag 403 into that same map, then reaches `_exact_shares` through `_rb._arm` (`test_imp3_rod_birdcage.py:329`):
- `_arm_rod._wrapped` sets 403 to σ_wire (`test_imp3_tip_q_entry48_def.py:111`). Its all-cells (P) gate reads `_exact_shares`' `p_src` and `sheet_field_total`, and its asserted control reads the two-tag `p_vol`.
- `_arm_with_403._wrapped` sets 403 to the gel (`test_imp3_rod_heating_matched.py:115`). `_gel_seen` and `examples/implants/15_rod_heating_matched_phantom.py` both go through it.

An unconditional raise would therefore turn `test_imp3_tip_q_entry48_def.py`, I41 and `imp:15` red. That contradicts the target's "the existing callers' gates unchanged". So:
- **The guard raises by default.**
- The case dict may carry an explicit acknowledgment key. Exactly the two `_wrapped` functions set it, each with its reason.

The ring-quadrature route builds a sweep with no `"problem"` key (`test_port_drive_superposition.py:1326`). The guard skips it, and the docstring says so. **The acknowledgment is named in the dashboard for the 11-07 weekly to ratify or revert.**

### §9 item 828 — 2026-11-05 (lane 1, known-issues "Resolves with"; `OPS-81` the two-tag helper's docstring and guard)

**Executes** the weekly's lane-1 target under the ruling above. Family `OPS-81`, attempt 1 of 4.

**Code.**
- **`tests/validation/test_port_drive_superposition.py`:**
  - Add the constant `TWO_TAG_ACK_KEY = "p_vol_two_tag_acknowledged"`.
  - Add `_assert_two_tag_scope(sweep)`.
    - `problem = sweep.get("problem")`. If it is `None`, return; this is the ring route (`:1326–1331`), documented as unguarded.
    - Otherwise collect three things:
      - every `material_map` key not in `(CONDUCTOR_CELL_TAG, PHANTOM_CELL_TAG)` whose `float(np.real(getattr(m, "sigma", nan)))` is not `== 0.0` (a missing or non-numeric σ counts as lossy, so the guard never passes silently);
      - `problem.material`, the default, if its σ ≠ 0;
      - `problem.phantom_material`, if it is set and `problem.phantom_tag` is not in (1, 3).
    - If that set is non-empty and `not sweep.get(TWO_TAG_ACK_KEY)`, raise `ValueError`. The message names the tags and says: "`_loss_power_w` sums tags 3 and 1 only; gate (P) on the all-cells integral (entry 48's form, `test_geo36_rod_f_clinical_drive.py:381`; I57's `_all_cells_resid`, `test_imp3_tip_q_entry48_def.py:204`), or set `case[TWO_TAG_ACK_KEY]` with the reason".
  - Call it as the **first statement** of `_loss_power_w` (`:205`).
  - Rewrite `_loss_power_w`'s docstring to state: the two-tag scope only, the guard, the ack key, the unguarded no-`problem` route, and the all-cells form.
- **`tests/validation/test_birdcage_power_identity.py:313`:** `_exact_shares`' docstring states the same scope and points to `_loss_power_w`. No code change; it reaches the guard at `:330`.
- **The two `_wrapped` functions** each add one line after the map write, e.g. `case[TWO_TAG_ACK_KEY] = "(P) gated over all cells by _all_cells_resid"` and `"I41 arm M: shares computed, p_vol never read"`. Import the constant by name.
- **New `tests/validation/test_ops81_two_tag_guard.py`,** with no mesh.
  - Cases are `{"problem": TimeHarmonicProblem(mesh=None, frequency_hz=64e6, material=HomogeneousMaterial(0, 1, 1), material_map={1: σ_wire, 3: gel, 403: …})}`.
  - `TimeHarmonicProblem` is a frozen dataclass with no `__post_init__` (`src/fem_em_solver/core/time_harmonic.py:101–121`).
- **No edit** to bands, gates, `src/`, `examples/ports/13_birdcage_asymmetric_drive.py`'s local copy (`:187`, a (U) caller, out of scope) or `examples/implants/15_*` (it inherits the key).

**Windows.** Bare, foreground, Bash timeout 660000. Pre-check: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
1. `scripts/testing/run_and_log.sh OPS-81-grep-pre "timeout -k 30 30 grep -rnwI -e _exact_shares -e _loss_power_w tests src examples scripts"`. This should print 46 lines, `OPS-80`'s count at `42dd675d`. Any other count is a new caller: classify it before editing.
2. Code.
3. W0 (smoke): `scripts/testing/run_and_log.sh OPS-81-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ops81_two_tag_guard.py -v --tb=short'"`. Predicted ≈ 30 s.
4. W1: `scripts/testing/run_and_log.sh OPS-81-W1-entry48 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_I57=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 600 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_tip_q_entry48_def.py -v -s --tb=short'"`. Measured 353 s (`20261104T104820Z_GEO-36-I57-W1d.log`).
5. W2: the same form, with `FEM_EM_IMP3_I41=1` and `tests/validation/test_imp3_rod_heating_matched.py`, label `OPS-81-W2-i41`. Measured 290 s (`20261025T201613Z_IMP-3-I41-W1b.log`).
6. W3: the same form, with `timeout -k 30 300` on `tests/environment tests/validation/test_birdcage_power_identity.py`, label `OPS-81-W3-port16`. Measured 92 s (`20260912T123429Z_PORT-14-step2d-10mhz.log`, 17 passed).

**Anchor.** The guard's contract on a synthetic case, asserted in W0.
- **Positive control:** `_loss_power_w(case_403, None, None)` raises `ValueError`, and the message contains `403` and "all-cells". Calling through `_loss_power_w` proves the guard is wired first.

**Control (asserted).**
- **Negative controls in W0, none of which raises:**
  - the same map with 403 at σ = 0 (the unmapped / vacuum case, `test_port_birdcage_four_port.py:472`);
  - 403 lossy with the ack key set;
  - a sweep without `problem`;
  - the stock tags-1/3 map.
- **Further positive controls in W0, which raise:**
  - a default material with σ = 0.5 and no 403;
  - a material without `sigma`.
- **The existing callers, unchanged:**
  - W1: 14 passed / 3 skipped (as at `:2028`). The two-tag control `test_control_two_tag_p1_reproduces_omitted_term` reads 2.8818e-06 (`:2010`).
  - W2: 16 passed (as at `20261025T201613Z_IMP-3-I41-W1b.log:3791`).
  - W3: 0 failed, counted per rank by token.
  - A last-digit drift in W1 or W2 is the `OPS-43` MUMPS observation, not a red. The verdicts are what is asserted.

**Tier.** Heavy: W1 and W2 exceed 180 s. `-n 2`, every window ≤ 600 s container-side. **Class:** ops / docs. **Cost:** ≈ 3.2 measured (the class median) / 40 predicted min. The compute is ≈ 12.5 min.

**Traps already paid for.**
- No module-level import of `test_birdcage_power_identity` in the superposition module: that is the import cycle at `:224–238`.
- `pgrep -c python3` is the orphan check. `pgrep -f "python3 -m pytest"` is denied by `bash_guard.py`.
- Count verdicts per rank, by token (item 644).
- `grep -w` hides wrappers. The case is mutated inside `_wrapped`, not at the call site (`OPS-80.md`'s traps).
- The ack key goes on the **case dict**, not on `problem`: the frozen dataclass is replaced in `_wrapped`. Set the key after the `replace`.
- The guard runs before any assembly.
- Never loosen I57's bands.

**Scope.** The helper's docstrings, the guard and two ack lines. **Not claimed:**
- widening the helper (the weekly ruled it is not widened);
- any (P) re-read;
- `examples/ports/13`'s local copy;
- the unguarded ring route.

**Negative result.**
- W1 or W2 red: the ack wiring is wrong. Stop, park on `attempt/OPS-81-*`, and the entry stays open.
- A (U) caller raises in W3: the stock map holds a lossy tag outside (1, 3). That is a finding. Name it with file:line in the known-issues Scope, and change nothing.
- There is no `PREDICTED RED:` line.

**Status it can move.**
- `OPS-81` ⬜ → ✅.
- Known-issues 2026-11-03 → RETIRED in the same commit:
  - Its heading becomes `### ✅ RETIRED <date> (`OPS-81`, §9 item 828, <slot> implementer slot — the scope docstrings and the guard landed; the positive control raised and the callers' gates held) — was 🟡 OPEN 2026-11-03 (…unchanged…)`.
  - A last row is added: `| **Retired by** | `OPS-81` (`<sha>`): the two-tag scope stated, the all-cells form named; `_assert_two_tag_scope` raises on a lossy tag outside (1, 3) unless the case carries `TWO_TAG_ACK_KEY` (positive control `<W0 log>:<line>`); the acknowledging wrappers and the (U) caller stay green (`<W1>`, `<W2>`, `<W3>`); the helper is not widened (weekly 2026-11-04). |`

**Depends on:** nothing.

**Outcome (2026-11-05, 12:00 slot, attempt 1 of 4) — ✅.**
- Pre-grep: 46 matching lines, `OPS-80`'s count; no new caller (`20261105T182655Z_OPS-81-grep-pre.log`).
- Code as briefed: `TWO_TAG_ACK_KEY`, `_sigma_is_zero`, `_assert_two_tag_scope` (first statement of `_loss_power_w`), both docstrings, one ack line in each `_wrapped`. No band, gate, `src/` or example edit.
- W0 (smoke, 3 s): 12 passed per rank. The positive control's message names `tag 403` and "all-cells" (`20261105T182800Z_OPS-81-W0.log:63`); default σ = 0.5 and a σ-less material raise; σ = 0, acknowledged, no-`problem` and stock-1/3 cases pass the guard.
- W1 (302 s): 14 passed / 3 skipped per rank, two-tag control 2.8818e-06 (`20261105T182816Z_OPS-81-W1-entry48.log:2019`, `:2028`).
- W2 (310 s): 16 passed per rank (`20261105T183330Z_OPS-81-W2-i41.log:3791`).
- W3 (134 s): 18 passed per rank, 0 failed, no refusal (`20261105T183847Z_OPS-81-W3-port16.log:2007`, `:2010`). The brief's 17 was the 09-12 count; `tests/environment` has gained a test since. No (U) caller holds a lossy tag outside (1, 3).
- Orphan check `pgrep -c python3` read 0 before every heavy window and after the last.
