# OPS-72 — history

Opened 2026-10-26 by the weekly review (lane 1; owner of known-issues 2026-10-26). The §7 row carries the Done-when; this file carries the briefs and outcomes.

### §9 item 665 — 2026-10-27 (§10 lane 1, `OPS-72`: the migration sweep's cost; lane-1 revalidation)

*Written by the 2026-10-27 03:00 daily review from the §7 row's Done-when (`PROJECT_PLAN.md` `OPS-72`), unchanged. No band, gate or scope moves: the checker's roots (`src`, `tests`, `examples`, `scripts`) and every assertion in `tests/environment/test_dolfinx_api_migration.py` stay as they are.*

**The cost, measured (the before).** `tests/environment` at `-n 1`: 3 passed in 144.42 s (`20261026T021500Z_MAT-3-step1b-envprice.log`; known-issues 2026-10-26, `:31–39`). Per test: `test_negative_control_detects_reverted_migrations` 78.76 s (a baseline sweep of `src` + `tests` plus a sweep of the temp copy), `test_no_unmigrated_call_sites` 39.01 s (a third sweep of `src` + `tests`), `test_filed_survivors_outside_the_gated_roots_are_unchanged` 26.63 s (`examples` + `scripts`, then `scripts/probes`). The checker is loaded by `importlib` in a module fixture (`test_dolfinx_api_migration.py:79–85`), and nothing in it or the checker mentions MPI, so **every rank runs all of it**. The last five `W-env` windows measured 102–120 s (`test-results.md:5462–5493`).

**Anchor (the Done-when's measurements, asserted by the windows' own verdicts).**
- (V) Every `tests/environment` verdict unchanged: the same test IDs, all passing, at `-n 1` and `-n 2`, counted per rank by token.
- (N) The planted negative control still detects: `status=pass` and `applied=6` (`:121–124`, unedited).
- (R) The checker's CLI `RESULT:` line byte-identical on a clean tree before and after the change: `python3 scripts/testing/check_dolfinx_api_migration.py` (default `--roots src tests`, `:778`) run once on the parent commit and once on the change, each line pasted into the outcome, compared by `cmp` or string equality in the outcome note.
- (T) The full `tests/environment` ≤ 60 s by the harness footer at `-n 1` **and** at `-n 2`.

**Negative control.** The before run is the control: 144.42 s at `-n 1`, ≈ 2.4× the 60 s bar. A change that only moves the sweep to rank 0 cannot meet (T) at `-n 1`, because one rank still does three `src` + `tests` sweeps there. So the route has to remove duplicated work, not just ranks.

**Route (in order; stop at the first stage that meets (T)).**
1. **Profile, no edit** (W0): one `checker.sweep(["src", "tests"], REPO_ROOT)` under `cProfile` at `-n 1`, top 15 cumulative entries printed. That shows whether the time is parsing, signature resolution or I/O. Paste the top 5 lines into the outcome.
2. **One sweep per root set per session, on rank 0.** Add a module-scoped fixture that computes the `src` + `tests` report once on rank 0 and broadcasts it (`MPI.COMM_WORLD.bcast`; `SweepReport` / `Finding` are dataclasses, `check_dolfinx_api_migration.py:90–103`; check they pickle). Then `test_no_unmigrated_call_sites` reads it, and `negative_control` takes it as its baseline through a **new keyword-only argument** `baseline: SweepReport | None = None` (default `None` keeps the CLI path byte-identical). Rank 0 runs the negative control and the survivors sweep and broadcasts `(status, captured_output)` / the reports. Each rank asserts on the same broadcast values, so the verdicts stay per rank.
3. **Only if 1–2 do not reach 60 s: memoise the expensive resolution** that the profile names (e.g. one signature lookup per distinct API rather than per call site), as an internal cache inside the checker. The `RESULT:` line (R) is the proof that it changed nothing.

**Do not** narrow what any sweep reaches (the Done-when's "skips nothing it checked before"): the negative control still sweeps the whole temp copy, and the file / call-site / API floors stay. If stage 3 is still over 60 s, stop: report the measurement and leave the row ⬜ with the stage-by-stage times. Narrowing the control's scope is the weekly's call.

**Windows** (bare harness calls, foreground, Bash `timeout: 660000`; real build is fine, since `tests/environment` runs either way, but use the complex image as every `W-env` has: `source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1`):
- **W0** `OPS-72-W0-profile` (`-n 1`, `timeout -k 30 300`): `python3 -c` running the profile of route step 1 inside the container, with `PYTHONPATH=/workspace/src:/workspace`. It is a measurement, not a pytest run, so it does not trip the pytest-routing guard; keep it one bare harness call.
- **W-rb** `OPS-72-W-rb` (`-n 1`, `timeout -k 30 120`): the CLI on the **parent** commit's checker (`git show HEAD:scripts/testing/check_dolfinx_api_migration.py > /workspace/logs/ops72-checker-before.py` from the host first; then run that copy with `--roots src tests` in the container). Its `RESULT:` line is (R)'s before.
- **W1** `OPS-72-W1-n1` (`-n 1`, `timeout -k 30 300`): `python3 -m pytest tests/environment -v --tb=short --durations=5`, the `W-env` command verbatim (`20261027T003839Z_WF-15-s8-Wenv.log:12`).
- **W2** `OPS-72-W2-n2` (`-n 2`, `timeout -k 30 300`): `mpiexec -n 2 python3 -m pytest tests/environment -v --tb=short --durations=5`.
- **W-ra** `OPS-72-W-ra` (`-n 1`, `timeout -k 30 120`): the CLI on the changed checker, `--roots src tests`. Its `RESULT:` line must equal W-rb's.

**Executor:** implementer. **Tier:** smoke (target), `-n 1` / `-n 2`, the before run is standard. **Class:** lane-1 revalidation (the weekly placed `OPS-72` on lane 1). **Cost:** ≈ 6.0 measured / 30 predicted slot-minutes. **Depends on:** nothing. **Family:** `OPS-72` step 1, attempt 1 of 4.

**Traps already paid for.**
- `bash_guard.py` denies pytest runs that bypass the harness and any harness call that is not a bare whole command. W-rb's host-side `git show … >` is its own Bash call, never chained to the harness call.
- A string `pgrep -f "python3 -m pytest"` is denied by the guard; `pgrep -c python3` in the container is the check.
- Count pytest verdicts per rank by token. In a 2-rank stream a trailing `PASSED` can be the other rank's.
- `capsys` captures only the calling rank's stdout. Under route step 2 the non-root ranks must assert on the broadcast output string, not on `capsys`.
- A collective call (`bcast`) inside a test that one rank skips deadlocks the window. Every rank must reach every `bcast`; keep them in fixtures, not in branches.

**Scope.** The cost of `tests/environment`. **Not claimed:** any change to what the sweep checks, its roots, its floors or the negative control's six reversions; any other environment test.

**Negative result.** (V), (N) or (R) red ⇒ the change altered the gate: revert it, park on `attempt/OPS-72-*`, BLOCKED. (T) over 60 s with (V)/(N)/(R) green ⇒ land the stages that are green as an improvement, record each window's footer, keep the row ⬜ and known-issues 2026-10-26 open with the new measurement, and the environment-first ruling keeps binding. No `PREDICTED RED:` line.

**Status it can move.** `OPS-72` ⬜ → ✅ when (V), (N), (R) and (T) are all green; known-issues 2026-10-26 retires in the same commit; the environment-first ruling is then lifted by the next daily review (not by the slot). Otherwise the row records the measured time.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 665 — `OPS-72` — per `docs/planning/chunks/OPS-72.md` "§9 item 665". The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Never edit an assertion, floor or reversion in `tests/environment/test_dolfinx_api_migration.py`; the checker's CLI default stays byte-identical. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-27, 07:30 CDT slot, implementer, no spawn) — ✅, all four green; route steps 2 and 3 both landed.**
- **W0 profile** (`20261027T124634Z_OPS-72-W0-profile.log:63–70`, one `src` + `tests` sweep under cProfile, 108.10 s profiled): `sweep` 108.097 cum; `sweep_file` (582 calls) 99.120; `ast.walk` (22 716 788 calls) 87.526; `iter_child_nodes` 59.279; **`enclosing_bindings` (1 998 calls) 52.423** — a full re-walk of every scope's subtree per call site. Parsing and signature resolution are not the cost; tree walking is.
- **Change.** Checker (internal only): `scope_map` stores each scope's set of contained node ids, built by one walk per scope, and `enclosing_bindings` tests membership, giving the same set union. `sweep_file` walks its tree once and reuses the list for its four passes (the tree is never mutated, so the order is the same). `negative_control` gains keyword-only `baseline: SweepReport | None = None`; `None` (the CLI) sweeps as before. Gate module: three module fixtures compute the `src` + `tests` report, the control's `(status, stdout)` and the `examples` + `scripts` / `scripts/probes` reports **on rank 0 and `bcast`** them. Every rank asserts on the broadcast values, so `capsys` is gone from the control test. No assertion, floor, root or reversion is edited.
- **(R)** W-rb (parent checker, `20261027T124558Z_OPS-72-W-rb.log:56`) and W-ra (changed, `…124919Z_OPS-72-W-ra.log`) print the identical line `RESULT: files=582 calls=1471 apis=39 methods=17 violations=0 uncheckable=0 shadowed=0` (string-equal; CLI 28 s → 16 s).
- **(V)** W1 `-n 1` (the W-env command verbatim, no `-s`): 12 passed. The sorted `ID PASSED` list `cmp`-identical to `20261027T003839Z_WF-15-s8-Wenv.log`'s. W2 `-n 2` (with `-s`, rule (g)): 24 `PASSED` tokens = 12 per rank, 0 `FAILED` / `ERROR`, two `12 passed` summaries.
- **(N)** `test_negative_control_detects_reverted_migrations` PASSED on every rank at both widths, so `status == 0`, `"status=pass" in out` and `"applied=6" in out` held on the broadcast output. The control still sweeps the whole temp copy.
- **(T)** Footers: W1 **44.90 s** (`12 passed, 22 warnings in 44.90s`, harness elapsed 46 s); W2 **44.73 s / 44.66 s** per rank (harness elapsed 46 s). The before was 144.42 s (`MAT-3` envprice) / 119.34 s (last W-env). Setup times at `-n 1`: control 18.27 s, `src` + `tests` sweep 15.77 s, ungated 10.14 s.
- Not done here: lifting the environment-first ruling (the next daily's).

### Audit of the closure — 2026-10-27 10:30 review

**PASS on §4, with the tier re-declared from smoke to standard by measurement; no demotion** (`auditor`, re-cited by the review). The harness logs exist with Status 0 and an elapsed time. The verification was executed by the slot, and the asserted gates are unchanged in `tests/environment/test_dolfinx_api_migration.py:131–195`.
- **Timing.** `20261027T124941Z_OPS-72-W1-n1.log:114` reads `12 passed … in 44.90s`, with harness elapsed 46 s at `:118`. `20261027T125038Z_OPS-72-W2-n2.log:145` and `:154` read 12 passed per rank (44.66 s / 44.73 s), 24 `PASSED` tokens.
- **(R).** `…124558Z_OPS-72-W-rb.log:56` and `…124919Z_OPS-72-W-ra.log:56` print string-equal `RESULT:` lines (re-read with `cmp`).
- **Tier.** The row read "smoke", but smoke's ceiling is 30 s (§5.1). The Done-when's own ≤ 60 s bar and the measured 46 s are standard-tier. The row's tier cell is corrected, as `MAT-3`'s was (10-26 weekly).
- **Caveats, recorded and not gating:**
  - The "before" of 144.42 s (`20261026T021500Z_MAT-3-step1b-envprice.log:119`) ran the migration module's 3 tests only. The last full 12-test `tests/environment` before the change was 119.34 s (`20261027T003839Z_WF-15-s8-Wenv.log:124`). The comparison is conservative, and the ≤ 60 s bar is met either way.
  - The windows ran on the parent commit `e38ee1b0` with the change uncommitted in the tree. The fixture-setup shape (`W1-n1.log:109–111`) and W-ra's 16 s against W-rb's 28 s show the change was present.
  - The CLI default `negative_control(baseline=None)` is not exercised by any logged run.

**The environment-first ruling (03:00 daily, 2026-10-26) is lifted by this review** (§9 digest).
