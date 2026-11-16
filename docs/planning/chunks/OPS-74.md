# OPS-74 — the four example closures' owed windows, and the first real `./run_examples.sh` through the harness

*History file for the §7 `OPS-74` row (`OPS-46`). Opened by the 2026-11-01 18:00 daily review under the 10-31 weekly's closing-window rule ("a closure whose last window is anything else is audited with the caveat and the daily re-queues the one window").*

**The defect.**
- `EX-166` … `EX-169` are the first closures since the rule, and all four closed on the dry-run's emitted line with `-T` and `FEM_EM_SETUP_FIGURES=1` added. All four audited PASS with that caveat (2026-11-01 18:00 journal §3).
- **The rule has never been observed to run.** No plain harness log carries a non-dry `run_examples` command line (grep of `docs/testing/logs/*.log`, 2026-11-01 18:00; `.gz` logs not swept).
- **Three standing texts say the opposite of the rule:** `scripts/run_examples.sh:30–43`, `docs/automation/implementer-run.md:133–150` and `.claude/agents/example-runner.md:57–67` all make emit-then-harness the supported path. The rule's premise is that the harness is sandbox-exempt, so `./run_examples.sh` inside it reaches the socket. Those texts' failures were all direct, unwrapped calls. That is a reading, not an observation.
- **`scripts/run_examples.sh:285` has no `-T`,** and a harness window has no TTY, so the real form may stop on "the input device is not a TTY". Not observed.

**Done-when.** Each of `EX-166` … `EX-169` has a Status 0 harness log whose Command line is `./scripts/run_examples.sh …`, cited in its §7 row beside the existing log; or the window's failure is recorded here and the rule returned to the weekly with the measurement.

### §9 item 738 — 2026-11-01 (ops; the 10-31 weekly's closing-window rule)

**Why.** The weekly told the daily to re-queue these, and each brief below this one inherits the rule. This item is also the rule's first test, so it goes first.

**Anchor.** Every gate in each example, unchanged, asserted by the run itself: exit 0 and the script's own final line. Per example, the asserted records equal the §7 row's at the printed digits, or the difference is named.

**Negative control.** `scripts/testing/run_and_log.sh OPS-74-Wc "./scripts/run_examples.sh -e mat:999 -n 2 -t 60"` reads a non-zero Status and runs no example: the footer tells a bad selector from a pass. Run it first.

**Execution, run inline, not through `example-runner`** (its definition says the opposite). Each is a bare harness call, cheapest first:
1. `scripts/testing/run_and_log.sh EX-169-close "./scripts/run_examples.sh -e mat:33 -n 2 -t 300"` (74 s)
2. `scripts/testing/run_and_log.sh EX-166-close "./scripts/run_examples.sh -e mat:32 -n 2 -t 480"` (117 s)
3. `scripts/testing/run_and_log.sh EX-168-close "./scripts/run_examples.sh -e imp:19 -n 4 -t 480"` (190 s)
4. `scripts/testing/run_and_log.sh EX-167-close "./scripts/run_examples.sh -e mri:37 -n 2 -t 560"` (211 s)

**Licensed, and only this:**
- **If step 1 stops on the TTY message before any example code ran:** add `-T` to `scripts/run_examples.sh:285`, and to the echoed `[dry-run]` line at `:283` so the two stay equal. Re-run step 1 under the same label. Keep the first log. Nothing else in that script changes.
- **`examples/materials/33_gel_box_pre_equilibration_bias.md:43` reads "All equal the gate's real-build records at the printed digits", which its own bullet above contradicts** (audit). Change that sentence to say (I), (c-ref) and (c-Dz) equal them, and (E) and (Sym′) differ at round-off on the complex build, with the new log's digits.
- **Each §7 row** gains the new log, tagged "closing window", and keeps the old one.

**Tier, ranks, wall clock.** The examples' own tiers, ≈ 10 min of compute in four calls, none over 600 s. Bash `timeout: 660000`, foreground. Pre-check `pgrep -c python3` 0 before each. **Class:** ops / docs. **Cost:** ≈ 5.6 measured (the class median; expect ≈ 15) / 25 predicted min. **Depends on:** nothing.

**Traps already paid for.** `run_examples.sh` is host-side: never wrap it in `docker compose exec` (Status 127, three times). A direct, unwrapped call from a sandboxed session gets `permission denied` from the docker API. `set -e` means one red example aborts. The harness call is the whole Bash command, nothing chained or piped. No `FEM_EM_SETUP_FIGURES`: the figures exist, and setting it would be a variation on the rule's command. A kill parks the item; never raise a `-t`.

**Scope.** Four windows and one flag. **Not claimed:** any new gate. **Not to edit:** the header at `scripts/run_examples.sh:30–43` or `implementer-run.md`; the review amends those from this item's result, and `.claude/agents/example-runner.md` is read-only in `main` (the operator's worktree session).

**Negative result.** Step 1 fails with `-T` in place and no example code run, for any reason: stop, BLOCKED, keep every log, skip steps 2–4. That is the measurement the weekly needs: the rule as written cannot run from a slot, and it decides between amending the rule and adding the sandbox exclusion (the operator's). An example gate red on this form that was green on the emitted form is a divergence: stop and report. No `PREDICTED RED:` line.

**Status it can move.** The §7 rows `EX-166` … `EX-169` gain a closing-window log; the audit caveats retire; the 10-31 rule gets its first observation, either way. `OPS-74` ✅ on four green windows.

**Outcome (2026-11-01, 19:30 slot): complete, `OPS-74` ✅.** The rule's form runs from a slot. Negative control `20261102T013018Z_OPS-74-Wc.log`: Status 2, "No mat example found for number: 999", no example run. Windows, each Status 0 with its script's own final line: `mat:33` 74 s, `mat:32` 112 s, `imp:19` 183 s, `mri:37` 203 s (`…013021Z_EX-169-close.log`, `…013146Z_EX-166-close.log`, `…013342Z_EX-168-close.log`, `…013649Z_EX-167-close.log`). `pgrep -c python3` read 0 before each.
- **No TTY stop:** compose v5.6.1's `docker compose exec` ran without `-T` under the harness, so the licensed edit to `scripts/run_examples.sh` was not made and the script is unchanged. The socket reached fine, which observes the rule's premise.
- **Records against the §7 rows, by `diff` of the old and new logs:** every asserted value equal at the printed digits except round-off terms. `mat:32`: (E) 2.218e-13 against 2.248e-13, the (Q-wf) / (Win-wf) rels ≤ 4.1e-15. `imp:19`: (E) 5.740e-12 against 5.910e-12, (Q-hd) / (Win-hd) rels ≤ 1.1e-15. `mat:33`: (E) and (Sym′), now on the guide. `mri:37`: (S) 6.598e-14 against 1.295e-15, (i) 1.767e-15 against 9.935e-16, cells, σ_max, (c-off) equal.
- **(S) on `mri:37` moved by a factor ≈ 51 between two runs of one script** and now clears `IDENTITY_BAND` 1e-12 by ≈ 15 times, not ≈ 770. Green, not diagnosed; the two runs differ in `FEM_EM_SETUP_FIGURES` and `-T` only.
- **Guide:** `33_gel_box_pre_equilibration_bias.md:43` corrected as licensed, after the window; the script does not read it. `OPS-74-Wdoc` `dead=0 guide=0` (Status 2, staleness only) and `OPS-74-Wpin` 20 passed ran after the edit.
- **Left to the review:** `scripts/run_examples.sh:30–43` and `implementer-run.md:133–150` still say the opposite of the rule, and `.claude/agents/example-runner.md:57–67` is the operator's worktree session.
