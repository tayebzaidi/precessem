# PORT-1 — history and briefs

*Created 2026-10-25 by the 03:00 daily review (OPS-46 split: the §7 row keeps its state line; briefs live here). The chunk itself closed 2026-08-15 (`‖S − Sᵀ‖/‖S‖ = 2.5494e-05` against the 1e-3 gate, two-torus fixture, step 4, `20260813T183606Z`). This file starts with the lane-1 item the 10-24 weekly gave the open padding-test entry.*

### §9 item 635 — 2026-10-25 (§10 lane 1, R1: the padding test under the owner rule)

*Written by the 2026-10-25 03:00 daily review from the 10-24 weekly's ruling (`PROJECT_PLAN.md` § 10, "Known-issues 2026-09-09 (`PORT-1`'s padding test, red on `main`, no xfail) — given an owner: lane 1, R1").*

**The red.** `tests/validation/test_port_gap_voltage_padding.py::test_box_enlargement_discriminates_between_the_two_routes` (`:163`) asserts `worst <= CONVERGENCE_THRESHOLD` (0.025, `:77`, assert `:224`) at `DISCRIMINATOR_AIR_PADDING` = 0.10 (`:71`). It reads **3.111508e-02** (known-issues 2026-09-09, `known-issues.md:93`).
- The voltage is the **test-side** path voltage: `_solve_gap_ports`, imported from `tests/validation/test_port_gap_voltage_impedance.py` (`:975`, used `:1267`), → `_path_voltage` (`:747`) → `_segment_voltage` → `evaluate_vector_field_parallel(e_field, points, comm)` (`:736`).
- It is called **with no `owner=`**, i.e. the legacy point evaluation, the family `POST-8` ruled on (✅ 2026-10-22).

**The weekly's question.** Under the owner rule (`owner="lowest_original_cell"`, `src/fem_em_solver/post/evaluation.py:224–249`), does the deviation fall under 2.5 %?
- If yes, the test is green by consequence and the entry retires.
- If not, the test is marked `xfail(strict=True)` with the entry as its reason, never loosened.

**Route (rule (c), additive; nothing else in the impedance module moves).**
- Add one keyword, `owner=None`, threaded through `_solve_gap_ports` → `_path_voltage` → `_segment_voltage` → `evaluate_vector_field_parallel(…, owner=owner)`. None keeps today's legacy call byte for byte. The impedance module's own tests (its 0.08 gates) call it unchanged.
- The padding module's fixture `gap_ports_padded` (`:95`) passes `owner="lowest_original_cell"` for **both** routes it compares, the estimator and the control. A one-route change would compare two evaluation rules.
- `src/fem_em_solver/ports/gap_voltage.py:284` is also legacy, and it is **not** edited here. Name it in the commit as the next owner-rule site.

**Anchor.** The test's own assertion, unedited: `worst <= CONVERGENCE_THRESHOLD`. This is the pre-decided 2.5 % bar on the two-torus estimator / control pair. If green, it is a quantitative closure by consequence of the owner rule.

**Controls.**
- **(Leg)** The same window, with `owner=None`, reproduces the red, 3.111508e-02 at `.6e`. This proves the change is the rule and not the fixture. It runs as its own window only if time allows (see Windows): its 2026-09-10 `-n 4` run already printed the value (`20260910T033937Z_TH-15.log.gz:1146`).
- **(Imp)** The impedance module's node that uses `_solve_gap_ports` re-runs green with the keyword at its default. Name the node from the module, and say which.

**Windows (bare harness calls, foreground; `pgrep -c python3` reading 0 before each).**
1. W0: `--collect-only` on both modules plus a keyword-default unit check, no solve.
2. W1 `PORT-1-R1-W1`: the padding module at **`-n 2`**, heavy, container `timeout -k 30 1150`.
   - **The price is the risk.** No `-n 2` run of this module has ever returned a footer: it timed out at 400, 500 and 590 s (`known-issues.md:96–97`; `20260909T184758Z_TH-15.log.gz:12`, `:785`). At `-n 4` it took 537 s (1 failed, 12 passed; `20260910T033937Z_TH-15.log.gz:1146`, `:1149`).
   - The weekly asked for `-n 2`. The heavy tier's 1200 s ceiling is the room it has, so the timeout is 1150 s. **If W1 is killed, do not re-run at a higher timeout or another width**: park and report the last flushed line (the shrink is the weekly's).
3. W2 `PORT-1-R1-W2` (only if W1 returned and the slot clock allows): (Imp).

**Tier:** heavy, `-n 2`. **Class:** lane-1 revalidation. **Cost:** ≈ 11.4 measured / 50 predicted slot-minutes. **Park (rule (d)) if W1 has not started by minute 30.** **Executor:** implementer. **Depends on:** nothing. **Family:** R1, attempt 1 of 4. **Daytime:** no. Never run beside another heavy window.

**Traps.**
- `mpiexec` needs `timeout -k 30`: a plain TERM leaves ranks burning (MAT-6 step 10). After any kill, `pgrep -c python3` before anything else.
- The owner rule is width-invariant by construction. It is the reason `-n 2` is the meaningful width here; do not "save time" at `-n 4`.
- The module is module-scoped (`gap_ports_padded`), so one node runs the whole fixture. Select by node ID, never `-k`.
- Never pipe the test run inside the harness (the footer would record the pipe's status).

**Scope boundary.** One test's evaluation rule. Not closed: `ports/gap_voltage.py:284` (legacy, named), the impedance module's own gates (unchanged), any band.

**A negative result means:**
- Deviation ≥ 2.5 % under the owner rule ⇒ mark the test `@pytest.mark.xfail(strict=True, reason="known-issues 2026-09-09: …")`, the house pattern (`tests/validation/test_port_package_sparameters.py:289–296`). Re-run W1 to show XFAIL, and keep the entry open with the owner-rule reading.
- W1 killed ⇒ park with the last line, the entry unchanged.

**Status it can move.** Known-issues 2026-09-09: green ⇒ **retired** with this commit; red ⇒ the test becomes strict-xfail with the entry as reason, and `main` is no longer red without a marker. §10 R1 → DONE either way. `PORT-1` stays ✅.

**Outcome 2026-10-25 (09:00 CDT slot, take-next spawn at minute 28) — BLOCKED, W1 not run.**
- Route done and parked on `attempt/PORT-1-R1-20261025T142815Z` (`e2df2a6e`): `owner=None` threaded `_solve_gap_ports` → `_path_voltage` and the retiling sum → `_segment_voltage` → `evaluate_vector_field_parallel(…, owner=owner)`; `gap_ports_padded` passes `owner="lowest_original_cell"`. No assertion touched; `ports/gap_voltage.py:284` not edited (next owner-rule site).
- W0 (`--collect-only`, both modules, `-n 2`, 3 s): 22 collected per rank, `20261025T142758Z_PORT-1-R1-W0.log`. The keyword-default unit check was not run.
- W1 not started. Two reasons, both structural: (1) the brief's own park clause (W1 not started by minute 30; the spawn arrived at minute 28 with edits and W0 still owed); (2) a spawned implementer's Bash call caps at 660 s and may not background a harness window, so it cannot host the brief's 1150 s container timeout. No `-n 2` run of this module has footered under 590 s, so a 590 s window would most likely be killed, which is not a result.
- Unblock: take 635 as a slot's **first** item, with an executor that can hold one foreground harness call ≥ 1180 s; or the weekly's shrink. Cherry-pick `e2df2a6e`, then W1 as briefed.

### Review ruling 2026-10-25 10:30 (item 635 → item 637; rule (6) re-pose, provisional)

*Binds over item 635's brief.* **The unblock condition cannot be met.** The Bash tool's foreground ceiling is 660 s in every session (slot and spawned executor alike), and no harness call may be backgrounded (rubric trap, three slots 2026-08-10/11). The 1150 s container timeout that the 03:00 review wrote into item 635 was unrunnable by construction — a briefing error, not the slot's. **Re-posed under §9 rule (6) as item 637, attempt 2 of 4 in the R1 family:** the same target (deviation ≤ `CONVERGENCE_THRESHOLD` under the owner rule, else strict xfail), on the same fixture, at **`-n 4`** with container `timeout -k 30 600`.
- **Why `-n 4` keeps the question.** The owner rule is width-invariant by construction (`POST-8` ✅ 2026-10-22), so its reading at `-n 4` is its reading at `-n 2`. The red it tests — 3.111508e-02 — was itself measured at `-n 4` (`20260910T033937Z_TH-15.log.gz:1146`), so `-n 4` is the like-for-like width for (Leg), and that log serves as (Leg) without a re-run. Item 635's trap "do not save time at `-n 4`" is withdrawn.
- **The price is the risk, still.** The legacy module took 537.36 s at `-n 4` (`…033937Z_TH-15.log.gz:1149`), so a 600 s window has ≈ 10 % headroom. If the owner-rule evaluation adds more than that, the window is killed — and that is a park, not a result. The 10-26 weekly owns the shrink then.
- **The 10-26 weekly ratifies or reverts this re-pose** (it changes the weekly's "re-run the module at `-n 2`" wording, not its target).

### §9 item 637 — 2026-10-25 (§10 lane 1, R1 re-posed at `-n 4`; rule (6), provisional; attempt 2 of 4)

*Written by the 2026-10-25 10:30 daily review. Item 635's brief stands except where this block differs: anchor, route, scope, negative-result disposition and status it can move are item 635's, unchanged.*

**Slot-first.** Take 637 only as a slot's **first** item (W1 must start by minute 15; the window alone is ≈ 10 min).

**Route.** `git cherry-pick e2df2a6e` from `attempt/PORT-1-R1-20261025T142815Z` (the `owner=` threading and the `gap_ports_padded` owner keyword; 29 insertions, 10 deletions across the two modules). Read the diff before the first window: `owner=None` must leave the impedance module's legacy call byte for byte.

**Windows (bare harness calls, foreground, Bash `timeout: 660000`; `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reading 0 before each).**
1. W0 `PORT-1-R1-W0`: `--collect-only` on both modules at `-n 2`, smoke (item 635's W0 command; it already read 22 collected).
2. W1 `PORT-1-R1-W1`, heavy, `-n 4`, complex, durable capture (§5.1 idiom verbatim, `; exit $rc` included):
   `scripts/testing/run_and_log.sh PORT-1-R1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 600 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_port_gap_voltage_padding.py -v -s --tb=short > /workspace/logs/port1-r1-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/port1-r1-w1-raw.log; cat /workspace/logs/port1-r1-w1-raw.log; exit \$rc'"`
3. W2 (Imp) only if W1 returned and the slot clock allows: the impedance module node that calls `_solve_gap_ports` with the keyword at its default, `-n 2`, standard (name the node in the journal).
4. If the test is red under the owner rule: add the strict xfail (item 635's negative branch), then **re-run W1 once** to show XFAIL — rule (i). That second window also needs ≈ 10 min, so it runs only if it can start by minute 40; otherwise commit the red reading with the marker added but the XFAIL window owed, and say so (the row records it as not re-run).

**Tier:** heavy, `-n 4`. **Class:** lane-1 revalidation. **Cost:** ≈ 7.4 measured / 25 predicted slot-minutes. **Executor:** implementer (a spawned executor gets the bare-harness sentence verbatim). **Depends on:** nothing. **Family:** R1, attempt 2 of 4.

**Traps (added to item 635's).** A W1 killed at 600 s is not a result: park on `attempt/*` with the last flushed line of `logs/port1-r1-w1-raw.log`; check `pgrep -c python3` before anything else; **do not** retry at a higher timeout or at `-n 6` in-slot. Never pipe the pytest call (the footer reads the pipe).

**A negative result means:** as item 635 (strict xfail with the entry as reason, entry kept open with the owner-rule reading). Killed ⇒ BLOCKED with the last line; the shrink is the 10-26 weekly's.

**Status it can move.** As item 635: known-issues 2026-09-09 retired (green) or the test strict-xfail (red); §10 R1 → DONE either way; `PORT-1` stays ✅.

**Outcome 2026-10-25 (12:00 CDT slot, first item, W1 started at minute 1) — BLOCKED, W1 killed at its ceiling, not a result.**
- Route: `e2df2a6e` cherry-picked as briefed; the impedance diff was read and `owner=None` is the default of `evaluate_vector_field_parallel` (keyword-only, `post/evaluation.py:229`), so the legacy call is unchanged. The cherry-pick is parked again on `attempt/PORT-1-R1-20261025T171120Z` (`9b8ed86f`); `main` carries no code change.
- W0 (`--collect-only`, both modules, `-n 2`): 22 collected per rank, `20261025T170042Z_PORT-1-R1-W0.log`. pgrep 0 before W1.
- W1 (`-n 4`, container `timeout -k 30 600`, durable capture): `[capture] rc=124`, elapsed 601 s, `20261025T170055Z_PORT-1-R1-W1.log`. It was killed inside the module-scoped `gap_ports_padded` fixture, so neither the padding assertion nor any other node in the module reported. Last flushed line (`:844`): `step 3b-xiv gapped(sigma = 0.000e+00): … Im Z_gap/I' = -3.976149386e+02 Ohm`, then the PETSc SIGTERM. pgrep 0 after the kill. The windows ran in that order; no XFAIL window, no W2.
- Reading: the legacy module took 537 s at `-n 4`, so the owner rule plus this run's load add more than the ≈ 10 % headroom the 10:30 ruling priced. The per-port and per-σ point evaluations (path orders up to 1025 points per segment, retiling, the σ ladders) all go through the owner rule now. Following the trap, no retry was made at a higher timeout or at `-n 6`.
- Unblock (the 10-26 weekly's): shrink so that one ≤ 600 s window reaches the assert. For example, a fixture that computes only what `test_box_enlargement_discriminates_between_the_two_routes` reads, or a width the weekly licenses. R1 family: attempt 2 of 4 spent.

### Review note 2026-10-25 18:00 (items 635 / 637, R1) — not re-posed; facts for the 10-26 weekly's shrink

R1 stays BLOCKED at attempt 2 of 4. Two attempts remain, and a third spent on an extrapolated width would leave one, so this review does not re-pose. These facts come from an `Explore` fact sheet, re-checked by this review on the cited lines:
- **The 637 kill was at 9 of 10 solves.** The kill line `20261025T170055Z_PORT-1-R1-W1.log:844` corresponds to legacy `20260910T033937Z_TH-15.log.gz:831`. After it, the legacy run did one more solve (26.9 s), the summary and the two tests.
- **Per-solve times match the legacy run** (W1 26.3–30.5 s; legacy 24.8–27.2 s). The extra ≈ 60–90 s is in the untimed parts: the mesh build, the owner-rule path-voltage sweep (8 quadrature orders up to 4 097 points × 2 ports × 2 columns, `impedance.py:1265–1270`) and the retiling sum (`:1284–1293`), or machine load.
- **Correction to item 637's outcome hypothesis:** the σ ladders do **not** go through the owner rule. `_closure_decomposition` (`:768–788`) still calls `_segment_voltage` without `owner`.
- **The padding assertion reads 3 of the fixture's 10 solves:** the two port-driven ones and the σ = 0 control (`test_port_gap_voltage_padding.py:107–125`; `CONVERGENCE_THRESHOLD` = 0.025 at `:77`, asserted at `:224`). The sibling structural test reads the same three. The ladders, the σ sweep, the closure and the retiling feed neither.
- **Width evidence:** `-n 2` → `-n 4` cut a σ×2 solve from 39.2 to 26.4 s (`20260909T184758Z_TH-15.log.gz:786` against `20260910T033937Z_TH-15.log.gz:803`). No `-n 8` run exists. Extrapolated, `-n 8` would be ≈ 380–450 s, but that is not a measurement.
- **Options for the weekly:**
  - (a) An additive keyword on `_solve_gap_ports` that skips what the padding tests do not read. This is ≈ 3 solves, ≈ 150 s at `-n 4`. The impedance module's default path is unchanged by construction; its own full re-run needs > 537 s at `-n 4`, so it is not a regression window.
  - (b) `-n 8` on the full fixture, extrapolated as above.

### §9 item 666 — 2026-10-27 (§10 lane 1, R1′ — the weekly's shrink, 2026-10-26; attempt 3 of 4)

*Written by the 2026-10-27 03:00 daily review from the 10-26 weekly's ruling (`PROJECT_PLAN.md` §10, "R1 … the shrink, ruled: R1′ (lane 1)"). Item 635's anchor, scope and negative-result disposition stand; this block changes the route and the window only. Item 637 is overtaken by it, not reverted.*

**Anchor (unchanged, imported).** `test_box_enlargement_discriminates_between_the_two_routes` (`tests/validation/test_port_gap_voltage_padding.py:163`): the worst per-column estimator/control deviation ≤ `CONVERGENCE_THRESHOLD` = 0.025 (`:77`, asserted at `:224`) on the padding-0.10 box, under the owner rule. The red under test is 3.111508e-02 at `-n 4` without the owner rule (`known-issues.md` 2026-09-09, `:119–120`; `20260910T033937Z_TH-15.log.gz:1146`). **Negative control (unchanged):** `WEDGE_ONLY_ESTIMATOR` 0.493653 against `WEDGE_ONLY_FLOOR` 0.20 (`:91–92`), as the module already asserts.

**Route — one additive keyword (rule (c)).**
1. `git cherry-pick 9b8ed86f` from `attempt/PORT-1-R1-20261025T171120Z` (the `owner=` threading; 2 files, +29 / −10). Read the diff first: with `owner=None` every legacy call is byte-for-byte what it was (`evaluate_vector_field_parallel`'s `owner` defaults to `None`, `post/evaluation.py:229`).
2. Add a keyword-only `padding_only: bool = False` to `_solve_gap_ports` (`tests/validation/test_port_gap_voltage_impedance.py:975`). When `True`, compute only what the two padding tests read: the mesh and its `cells` / `air_padding`, the two port-driven solves (`currents`, with their gap `voltages` / `driven`), the σ = 0 `control` (its `z21`), and `omega_m`. **Skip** the σ sweep, the control and production σ ladders, the closure decomposition, the retiling sum and the 8-order path-convergence sweep, and return `None` for their keys (not absent: a reader that touches one fails loudly). With `False` nothing changes: the default path is the code it was.
3. `gap_ports_padded` (`test_port_gap_voltage_padding.py:96`) passes `padding_only=True` alongside the owner keyword from step 1. Correct its docstring's "five solves" to the count it now runs.
4. Before W1, read every key the two padding tests access (`:109–112`, `:139–145`, `:180` and anything `test_box_enlargement…` reads further down) and confirm that each is computed under `padding_only=True`. Write that list in the outcome.

**Windows** (bare harness calls, foreground, Bash `timeout: 660000`; `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reading 0 before W1; complex build; the environment-first ruling holds, so only `test_complex_mode.py` goes first):
- **W0** `PORT-1-R1p-W0` (smoke, `-n 2`, `timeout -k 30 60`): `--collect-only` on both modules. 22 collected, as item 637's W0.
- **W1** `PORT-1-R1p-W1` (heavy, `-n 4`, durable capture, §5.1 idiom verbatim with `; exit $rc`):
  `scripts/testing/run_and_log.sh PORT-1-R1p-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 600 mpiexec -n 4 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_port_gap_voltage_padding.py -v -s --tb=short > /workspace/logs/port1-r1p-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/port1-r1p-w1-raw.log; cat /workspace/logs/port1-r1p-w1-raw.log; exit \$rc'"`
  **Predicted ≈ 150 s** (3 of the legacy 10 solves at 26–30 s each plus the mesh and the owner-rule evaluation of the two gap voltages; the weekly's figure). Over 600 s is a kill: park it, never raise the timeout.
- **W-env** (`-n 1`, `timeout -k 30 300`; the item edits `.py`): label `PORT-1-R1p-Wenv`, full `tests/environment`, as in `MAT-3.md` § item 648. Measured 102–120 s on the last five items.
- **XFAIL re-run** (only on a red): rule (i), the same W1 command under label `PORT-1-R1p-W1x`.

**Regression.** The impedance module's own default path is not re-run: it needs > 537 s at `-n 4`. The proof that it is unchanged is the diff: `padding_only` defaults to `False`, and every new branch sits under `if padding_only`. Paste the `git diff` hunk headers for `test_port_gap_voltage_impedance.py` into the outcome.

**Executor:** implementer. **Tier:** heavy (W1), `-n 4`, complex. **Class:** lane-1 revalidation. **Cost:** ≈ 6.0 measured / 20 predicted slot-minutes. **Depends on:** nothing. **Slot-first:** no (W1 is predicted under 3 min). **Family:** R1, **attempt 3 of 4**: a fourth may not be lettered, and a kill here leaves one.

**Traps already paid for.**
- Item 637's kill was **inside the module-scoped fixture**, so nothing reported. If W1 is killed again, the last flushed line of `logs/port1-r1p-w1-raw.log` is the whole result. Check `pgrep -c python3` before anything else.
- Never pipe the pytest call (the footer then reads the pipe's status). Count verdicts per rank by token: in a 4-rank stream a `PASSED` can be another rank's.
- `ports/gap_voltage.py:284` is the next owner-rule site (named in `9b8ed86f`'s message). It is **not** on this fixture's path; do not edit it here.
- `cell_tags.values` and `assemble_scalar` are rank-local. `omega_m` and `z21` must come out of the same reductions as before; the cherry-pick does not touch them.

**Scope.** The padding discriminator under the owner rule. **Not claimed:** the impedance module's ladders under the owner rule (the σ ladders still call `_segment_voltage` without `owner`, `_closure_decomposition`, `impedance.py:768–788`), `gap_voltage.py:284`, or any `PORT-1` record.

**Negative result.** As item 635: under 2.5 % ⇒ green by consequence and known-issues 2026-09-09 retires in the same commit. Over ⇒ add a strict `xfail` naming the entry as its reason, re-run W1 once to show XFAIL (rule (i)), and keep the entry open with the owner-rule reading as its new measurement. **Never loosen `CONVERGENCE_THRESHOLD`.** Killed ⇒ BLOCKED with the last flushed line, and the family is at 3 of 4: the next review may not re-letter it. No `PREDICTED RED:` line.

**Status it can move.** Known-issues 2026-09-09 retired (green) or re-measured under the owner rule (strict xfail); §10 R1′ → DONE either way; `PORT-1` stays ✅. Delete `attempt/PORT-1-R1-20261025T171120Z` once its commit is on `main`.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 666 — `PORT-1` R1′ — per `docs/planning/chunks/PORT-1.md` "§9 item 666". Read it and any review ruling beneath it first. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Cherry-pick `9b8ed86f`, add the one keyword, and never edit `CONVERGENCE_THRESHOLD` or any assert. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome 2026-10-27 (07:30 slot) — DONE: red under the owner rule, strict xfail.**
- `9b8ed86f` cherry-picked onto `main` (new hash; `attempt/PORT-1-R1-20261025T171120Z` left for the daily review to dispose of).
- **Keys the two padding tests read**, all computed under `padding_only=True`: `air_padding`, `cells`, `omega_m`, `control["z21"]`, `control["current_prime"]`, `control["current_prescribed"]`, `control["imag_ratio"]`, `currents[tag]["voltages"]` (the gated fine order, same `_path_voltage(..., owner=owner)` call as the sweep), `currents[tag]["driven"]`. `None` under `padding_only`: `sigma_sweep`, `control_sigma_ladder`, `production_sigma_ladder`, and per column `closure`, `retiling`, `path_voltages_by_order`.
- **Diff hunk headers, `test_port_gap_voltage_impedance.py`** (every new branch under `padding_only`, default `False`): `@@ -990 +990,6` (signature), `@@ -999,0 +1005,7` (docstring), `@@ -1281,3 +1293,8`, `@@ -1286,6 +1303,12` (gate-order-only path V), `@@ -1294 +1317` (closure), `@@ -1301 +1324` (retiling), `@@ -1344 +1367` (sweep print), `@@ -1409 +1432,4` (closure print), `@@ -1568 +1594,2` (control σ ladder), `@@ -1683 +1710,2` (σ sweep), `@@ -1759 +1787,2` (production ladder), `@@ -1909 +1938`, `@@ -1911,2 +1940,4` (return `None`s).
- **W0** `20261027T125429Z_PORT-1-R1p-W0.log`: 22 collected, 2 s. **W1** `20261027T125439Z_PORT-1-R1p-W1.log` (`-n 4`): **1 failed, 5 passed in 218.41 s**, elapsed 220 s (prediction ≈ 150 s); FAILED ×4 ranks. Worst deviation **3.101935e-02** (`:728`; −3.0515e-02 / −3.1019e-02, `:726–727`) vs 2.5 %; without the owner rule it was 3.111508e-02, so rank ownership moves it ≈ 1e-4 and is excluded as the cause. Negative control: wedge-only ratio 0.5181 (deviation −0.4819 > 0.20). 195 837 cells, mesh 31.9 s, solves 33.3 / 27.3 s, control 33.3 s.
- Strict `xfail` naming known-issues 2026-09-09 added; **W1x** `20261027T125835Z_PORT-1-R1p-W1x.log`: XFAIL ×4, **5 passed, 1 xfailed in 207.73 s**, Status 0, elapsed 209 s. **W-env** `20261027T130209Z_PORT-1-R1p-Wenv.log` (`-n 1`): 8 passed / 4 skipped in 43.87 s, elapsed 45 s.
- `CONVERGENCE_THRESHOLD` and every assert unmoved. The impedance module's default path was not re-run (its proof is the diff above).


### Review note 2026-11-01 03:00 (R1 banked by the 10-31 weekly; the branch disposed of)

**R1 is banked at attempt 3 of 4 and the family is closed** (§10, 10-31 rulings and epitaph; known-issues 2026-09-09 re-headed OBSERVATION with its strict xfail as the record). The entry retires when the padding module is deleted or a consumer reappears. `attempt/PORT-1-R1-20261025T171120Z` (tip `9b8ed86f`) is deleted, as the weekly allowed: its one commit was cherry-picked onto `main` by item 666, and the 3.1 % reading is on `main` in `20261027T125439Z_PORT-1-R1p-W1.log`.
