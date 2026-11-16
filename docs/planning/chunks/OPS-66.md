# OPS-66 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 139 — 2026-10-01 (ops; queued by the 2026-10-01 03:00 review; moved verbatim from §9 by the 10:30 review)

139. **`OPS-66` (new row) — `rotate_plan_archive.py known-issues` finds the four undated RETIRED entries** (implementer; independent; smoke; no MPI; real; tooling + test; **≈ 36 slot-min**: ≈ 0.5 min of windows + 15 fixed + ≈ 20 code). Discharges the 2026-09-30 weekly's flag 8; row §7 `OPS-66` (opened this review).
   - **The defect.** Not `_ki_closed` (`scripts/maintenance/rotate_plan_archive.py:171–177`), which already accepts all four, but the local `entry_re = re.compile(r"^(### |## \d{4}-\d{2}-\d{2})")` at `:182`. The four headings (`docs/testing/known-issues.md`): `## RETIRED 2026-09-20 (`OPS-59` (b)…`, `## ✅ RETIRED 2026-09-21…`, `## ✅ RETIRED 2026-09-20 (`PORT-20`…`, `## ✅ RETIRED 2026-09-27…`. The weekly's runs read "10 RETIRED entries" (`20261001T020036Z_WEEKLY-ROTATE.log:35`) and left these four.
   - **Code.** Lift to `KI_ENTRY_RE = re.compile(r"^(### |## (?:✅\s*)?(?:(?:RETIRED|RESOLVED|FIXED|CLOSED)\s+)?\d{4}-\d{2}-\d{2})")`, read by `rotate_known_issues`; docstring `:51–59` to match. New `tests/unit/test_rotate_known_issues_entries.py`, loading the tool via `importlib.util.spec_from_file_location` (`test_doc_reference_exit_codes.py:40`'s pattern).
   - **Anchors (asserted).** (A1) the four literal headings match `KI_ENTRY_RE` and `_ki_closed`. (A2) over every non-fenced heading of the live file and the archive: old ⊆ new and new − old ⊆ {headings `_ki_closed` accepts}; on the live file new − old is exactly 4. (A3) `## How to check…`, `## Failing tests`, `## Non-test issues`, `## Recording a new entry` and `KI_INDEX_HEADING` are not entries. The two struck-through "🟡 OPEN" clauses inside two of the headings stay closed (the OPEN test is anchored at the heading start): assert it.
   - **Negative control (predicted, binary).** The pinned pre-change tool (`git show e5af2cae111b35c95a5b17022ae299b7e5d1139e:scripts/maintenance/rotate_plan_archive.py`, the file's last commit — never `HEAD:`) into scratch under the gitignored `/workspace/logs/ops66-negative-control/` (with `git -c safe.directory=…`, `test_log_label_literals.py:53–55`'s idiom); `KNOWN_ISSUES` monkeypatched on both modules to a scratch copy of the live file; `rotate_known_issues(date, dry_run=True)`: old reports 0, new reports 4 — asserted as the binary.
   - **Windows.** W1 (≈ 10 s): `scripts/testing/run_and_log.sh OPS-66-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 30 python3 -m pytest tests/unit/test_rotate_known_issues_entries.py -s -v --tb=short'"`. W2 (≈ 2 s, writes nothing): `scripts/testing/run_and_log.sh OPS-66-W2-dry "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 30 python3 scripts/maintenance/rotate_plan_archive.py known-issues 2026-10-01 --dry-run'"` — predicted `4 RETIRED entries`, all four index lines listed.
   - **Traps.** Never run without `--dry-run`: the move is the weekly's (10-03). No hard-coded live line numbers in the test. Pin `<sha>`, never `HEAD:`.
   - **Scope.** The entry pattern only. Not: `_ki_closed`, the `attempts` / `plan` / `history` modes, `weekly-review.md`.
   - **Status it can move.** `OPS-66` ⬜ → ✅; the 10-03 weekly's rotation then moves the four entries.
   - **Negative result.** A section header becomes an entry, or old ⊄ new: the pattern is wrong — fix once, then park. W2 lists ≠ 4: report the list and stop.
   - **Outcome (2026-10-01 12:00 slot, `implementer`): DONE.** W1 `20261001T170610Z_OPS-66-W1.log` Status 0, 2 s, 8 passed; W2 `20261001T170616Z_OPS-66-W2-dry.log` Status 0 — `4 RETIRED entries, 94 lines`, exactly the four. (A1) 4/4 match `KI_ENTRY_RE` and `_ki_closed` (old pattern 0/4); (A2) live 39 headings, old 29 / new 33, new − old = the four; archive 81, old 77 / new 77; (A3) no section header or `KI_INDEX_HEADING` admitted, the two struck "🟡 OPEN" headings stay closed. Control (pinned `e5af2cae`): old 0, new 4. **Deviation (disclosed):** so the test survives the 10-03 rotation, it reads the working `known-issues.md` while all four headings are present and otherwise falls back to the file at `32fabbec`; A2 also accepts the four in the archive.

### §9 item 384 — 2026-10-11 (weekly flag 8, `OPS-66` step 2; ops/docs, lane 1)

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.*

**Licence.** The 2026-10-10 weekly's flag 8 (PROJECT_PLAN.md §10, "Flagged for the 2026-10-11 03:00 daily review", item 8): "`rotate_plan_archive.py chunks`: a state line beginning with the glyph is doubled — `OPS-66`-class tool fix (strip a leading glyph), lane 1."

**Home.** The `chunks` subcommand was built by `OPS-46` (✅ 2026-09-12; `docs/planning/chunks/OPS-46.md`) and extended by `OPS-47` (`--narratives`). The weekly named the class `OPS-66` (rotation-tool defects the weekly finds), so the brief lives here as `OPS-66` step 2. The review may instead open a new row (next free `OPS-71`); the brief is unchanged either way. No band, Done-when or target is set or moved.

**The defect (read on `510f41f9`).** `_replacement` builds the status cell as `f" {parsed['title']} | {parsed['glyph']} {state} {dw}"` (`scripts/maintenance/rotate_plan_archive.py:428–431`). `parsed['glyph']` is the row's current glyph (`_parse_row`, `:401–405`). A spec state line written from the row's last ruling usually opens with that glyph ("✅ 2026-10-05 …"), so the cell becomes "✅ ✅ 2026-10-05 …". The 2026-10-10 weekly's run hit it on all twelve rows (`docs/testing/logs/20261011T020900Z_WEEKLY-ROTATE-CHUNKS.log:35–46`) and de-duplicated them by hand in `397b69d0` (its commit message: "the tool's prepended glyph de-duplicated by hand on all twelve"). The spec file itself was not kept. The script's last commit is `15bdce61ca16199c6b1f3e947b6f13ad020e6ed8` (2026-10-01).

**The fix (≈ 10 lines; the docstring at `:29–47` gains one sentence).**
- In `_replacement`, or in a helper `_strip_leading_glyph(cid, state, glyph)` called from `rotate_chunks` before `_check_state`, `state.lstrip()` is tested against `GLYPHS` (`:351`):
  - **(i)** If it opens with the row's own glyph, remove that glyph and the whitespace after it.
  - **(ii)** If it opens with a **different** member of `GLYPHS`, raise `RowError(f"{cid}: state line opens with {g} but the row's status is {glyph}; the rotation never changes a status")`.
  - **(iii)** Otherwise, unchanged.
- **(ii) goes beyond the weekly's letter ("strip a leading glyph") and is the review's to accept.** The reason: a strip-any rule would silently drop a status the reviewer typed, and the tool's contract is a verbatim move ("never edits a line it moves", `:23`). If the review declines (ii), the leading glyph is stripped whatever it is, and test (G3) asserts the row's glyph survives.
- `_check_state`'s digit and sentence checks run on the stripped state, so their behaviour is unchanged (a glyph carries no digit and no sentence stop).

**Test: new `tests/unit/test_rotate_chunks_glyph.py`.** It loads the tool with `_load` / `importlib.util.spec_from_file_location`, as `tests/unit/test_rotate_known_issues_entries.py:68–75` does, and uses `git -c safe.directory=…` (`:78–89`). Every fixture is synthetic and written under `tmp_path`: a minimal plan with a `## 7.` heading, a `| ID | Title | Status | Tier |` table, and two rows `ZZ-1` (✅) and `ZZ-2` (🟡) whose history carries the digit runs the state lines use, closed by a `## 8.` heading. The chunks dir is empty under `tmp_path`, and every spec is written by the test.

**Anchors (asserted, exact).**
- **(G1) Glyph-led state.** Spec `=== ZZ-1` / `✅ 2026-10-10 state.` through `rotate_chunks(spec, plan, chunks_dir, dry_run=False, date="2026-10-11")`: rc 0, and the row's new status cell begins `| ✅ 2026-10-10 state.` with **exactly one** `✅` between the title/status boundary and `*History:`.
- **(G2) No-glyph state, unchanged behaviour, byte for byte.** Spec `=== ZZ-2` / `state without a glyph 2026.` (digits present in the row), with and without a `Done-when: ` line, run once through the **new** tool and once through the **pinned pre-change tool** on identical fixture copies: plan bytes, chunk-file bytes and `--dry-run` stdout are all equal. The same equality holds on a `--census` run.
- **(G3) Mismatched glyph** (as ruled). Spec `=== ZZ-1` / `🟡 2026-10-10 state.` on the ✅ row: rc 1, `REFUSED: ZZ-1: state line opens with 🟡 …` printed, plan bytes unchanged, no chunk file created. If the review declines (ii): rc 0 and the cell opens `✅ 2026-10-10`.
- **(G4) The variation-selector glyph.** A `⚠️` row (U+26A0 U+FE0F, as in `GLYPHS`) with state `⚠️ 2026-10-10 x.` → exactly one `⚠️`.

**Negative control (rule (e)): the pinned pre-change tool, never `HEAD:`.**
- `git show 15bdce61ca16199c6b1f3e947b6f13ad020e6ed8:scripts/maintenance/rotate_plan_archive.py` goes into the gitignored `/workspace/logs/ops66-step2-control/` (OPS-66 step 1's idiom). That is the file's last commit before the fix, i.e. the fix commit's `^` version.
- Run on the (G1) fixture, it writes `✅ ✅ 2026-10-10` (two glyphs). **Asserted as a binary:** old = 2, new = 1.
- It is backed by the same comparison on real data: the weekly's twelve doubled rows (`…WEEKLY-ROTATE-CHUNKS.log`, `397b69d0`).
- A control that copied `HEAD:` would go stale the moment the fix lands. That trap cost `OPS-43` (a) a window (rubric item 4).

**Printed only (real-data replay, no assertion).**
- Inputs: `git show 38aa4afd:PROJECT_PLAN.md` (`397b69d0^`, the pre-rotation plan) and a spec reconstructed from `397b69d0:PROJECT_PLAN.md`. The reconstruction takes each row's state = the row glyph + the text up to ` **Done-when:** ` or ` *History:`, and the Done-when line from the same row. It covers the nine rows whose titles were not hand-condensed: `OPS-70`, `GEO-36`, `MAT-15`, `WF-1`, `WF-10`, `WF-12`, `WF-15`, `WF-16` and `WF-17` (`MAT-20`, `WF-8` and `EX-89` were, per the commit message).
- Run into a `tmp_path` plan and chunks dir. Print, per row, whether the new tool's row equals `397b69d0`'s row byte for byte, and the old tool's doubled-glyph count.
- Predicted: 9 / 9 equal, and old doubled = 9. A mismatch is a reading (another hand edit in `397b69d0`), never a red.

**Window.** Smoke. Foreground; the harness call is the whole Bash command, bare; Bash `timeout: 660000`; `-s`. No MPI and no complex build (OPS-66 step 1's W1 ran this way, 8 passed, 2 s; outcome above).

`scripts/testing/run_and_log.sh OPS-66-s2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 30 python3 -m pytest tests/unit/test_rotate_chunks_glyph.py tests/unit/test_rotate_known_issues_entries.py -s -v --tb=short'"`

The OPS-66 step 1 module rides along: it loads the same tool, and standing rule (c) asks that a pre-existing gate on an edited file re-runs green. **Never run the tool on the live `PROJECT_PLAN.md`:** the rotation is the weekly's, and every fixture here is a copy.

**Cost.** ≈ M measured / ≈ 25 predicted min (15 fixed + ≈ 8 code + W1 ≈ 0.1 + journal). Class: ops/docs (median 3, carried; this one carries code). Depends on nothing.

**Traps already paid for.**
- Pin `<sha>`, never `HEAD:`.
- `git` in the container needs `-c safe.directory=/workspace`.
- No pipe or filter in the harness command; `-s` always.
- No hard-coded live line numbers in the test.
- No edit after W1 (rule (i)).
- Never path-checkout a record file.
- The scratch copy goes under gitignored `/workspace/logs/`, never `tests/`.

**Scope.** Covered: the `chunks` replacement's leading glyph. Not covered: `--narratives`, `--census` output (beyond (G2)'s equality), the `attempts` / `known-issues` / `history` / `plan` modes, the title derivation (the three hand-condensed titles are a separate question for the weekly), and `weekly-review.md`.

**Negative result.** Report and stop. Park on `attempt/OPS-66-s2-*` and mark item 384 BLOCKED in the same commit (rule (d)).
- (G2) unequal → the fix touched the no-glyph path: report the first differing byte.
- The control does not double → the pinned tool is not the defective one: report the sha.
- The step 1 module red → a regression in the shared loader: report it.

No `PREDICTED RED:` line.

**Status it can move.** `OPS-66` stays ✅: its Done-when is the known-issues entry pattern, met 2026-10-01. The row gains one sentence: "Step 2 (§9 item 384, 2026-10-11): `chunks` strips a state line's leading row glyph (refuses a different one); the pinned pre-change tool doubles it, the new one does not; no-glyph output byte-identical (`…_OPS-66-s2-W1.log`)." The weekly's flag 8 is discharged, and its next `chunks` rotation needs no hand de-duplication; the weekly confirms on that run. **No rule (j) amendment**, as no glyph moves. If the review opens `OPS-71` instead, that row's Done-when is (G1), (G2) and the control asserted ⇒ ✅.

**Executor spawn prompt (verbatim).** "Execute §9 item 384 — `OPS-66` step 2 (weekly flag 8) — per `docs/planning/chunks/OPS-66.md` §9 item 384. Read the brief and any review rulings beneath it first. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Pin the control to `15bdce61ca16199c6b1f3e947b6f13ad020e6ed8`, never `HEAD:`. Never run the tool against the live `PROJECT_PLAN.md` or `docs/planning/chunks/`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/OPS-66-s2-*`, mark item 384 BLOCKED in the same commit, and report."

**Outcome (2026-10-12, 22:30 CDT slot fourth item, `implementer`): DONE.** Fix as drafted, (ii) included (no review ruling declined it): `_strip_leading_glyph` called from `rotate_chunks` before `_check_state`; docstring gains one sentence. W1 `20261012T035822Z_OPS-66-s2-W1.log` (Status 1, 2 s): the new module 7 passed — (G1) new cell ` ✅ 2026-10-10 state.`, one ✅; control (pinned `15bdce61`) ` ✅ ✅ 2026-10-10 state.` — binary old 2 / new 1; (G2) plan / chunk / `--dry-run` / `--census` byte-identical without and with `Done-when:` (430 / 112 / 201 / 455 B and 472 / 112 / 244 / 455 B); (G3) rc 1, `REFUSED: ZZ-1: state line opens with 🟡 but the row's status is ✅; the rotation never changes a status`, plan unchanged, no chunk file; (G4) ⚠️ old 2 / new 1. Replay (printed only): new == `397b69d0` on 9/9 rows (predicted 9/9); the old tool's status cell held exactly two glyphs on 4/9 (`WF-1`, `WF-10`, `WF-12`, `WF-15`; predicted 9) — the criterion was "exactly 2", and why the other five read otherwise is not diagnosed (a reading, not a red). **The one red is not this item's:** step 1's `test_a2_old_subset_of_new_and_difference_is_exactly_the_four` fails identically on clean `main` `fc95e41b` with the fix stashed (`20261012T035840Z_OPS-66-s2-W0-main.log`, 1 failed / 7 passed) — the archive now carries three later undated RETIRED headings (`ANS-7` 2026-10-06, `WF-17` 2026-10-08, `OPS-69` 2026-10-08) that A2's `archive ⊆ FOUR` never anticipated. Filed in known-issues, not fixed in passing (deviation from the brief's "step 1 red ⇒ park": that branch presumed a loader regression, which the clean-`main` run rules out).

**Review ruling 2026-10-12 03:00 (daily) on item 384.** Accepted as landed; the deviation (file, do not park) is ratified — the clean-`main` W0 is the evidence that the red is step 1's test data, not the tool. The one-line fix is queued as item 393 below.

### §9 item 393 — 2026-10-12 (`OPS-66` step 3: A2's archive clause re-stated to its own docstring; ops / docs, lane 1)

*Written by the 2026-10-12 03:00 daily review.*

**What.** `tests/unit/test_rotate_known_issues_entries.py::test_a2_old_subset_of_new_and_difference_is_exactly_the_four` is red on `main` (known-issues "🟡 OPEN 2026-10-12 (found by `OPS-66` step 2 …)"). Its last clause, `assert diffs["archive"] <= set(FOUR)` (`:139–141`), asserts a *snapshot* of the archive — true only until the next rotation moved another undated RETIRED heading there, which the 10-06 / 10-08 rotations did three times. The module docstring (`:9–11`) states A2 as "new − old ⊆ {closed headings}; new − old is exactly the four" — the second half for the **live** file (the four's home, or `KI_PINNED_SHA` once they left). The archive half of that statement is already asserted at `:131–132` (`unclosed` empty, both labels). So the fix restores the clause to its pre-registered letter; it does not loosen it.

**Edit (tests only; no tool change).**
1. Replace `:139–141` with a durable archive clause: `assert set(FOUR) <= diffs["archive"] | diffs["live"]` (the four are still found, live or archived) — plus a print of the archive's new − old count and headings (60-char cut). Keep `:131–132` (the closed-only check) unchanged for both labels; keep the live `== set(FOUR)` clause unchanged.
2. Update the comment at `:139–140` to say why the snapshot clause was replaced (this item, the known-issues entry).

**Anchors.**
- (A2-live) unchanged: live new − old == the four (pinned live text where the four have left).
- (A2-closed) unchanged: every new − old heading in live and archive is closed (`tool._ki_closed`).
- (A2-four) new: the four ⊆ live ∪ archive new − old.
- The whole module green: 8 passed (A1, A2, A3, NC and the rest), plus `tests/unit/test_rotate_chunks_glyph.py` 7 passed (standing rule (c): the same tool is loaded; it is not edited).

**Negative control (asserted; backed by `20261012T035840Z_OPS-66-s2-W0-main.log`, the same comparison on the same archive).** In the same test, before the new clause, compute `diffs["archive"] - set(FOUR)` and assert it is **non-empty** (it read three headings on `fc95e41b`) — i.e. the old snapshot clause would still fail on today's archive, so the replacement is what turns the module green, not a change in the data. Print the three (or more) headings. Second control (asserted): `(A2-closed)` on a copy of the archive text with one planted heading `## RETIRED-NOT 2026-10-12 planted` appended is checked only if `tool.KI_ENTRY_RE` admits it and `tool._ki_closed` rejects it — if no planted form is admitted-but-unclosed, say so and print the attempt (the closed check then guards nothing new, a reading).

**Tier, cost.** Smoke, no MPI, ≈ 2 s window. ≈ 3 measured (ops / docs median) / 15 predicted min. Depends on nothing.

**Window (bare, whole, foreground, Bash `timeout: 660000`):**

`scripts/testing/run_and_log.sh OPS-66-s3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 30 python3 -m pytest tests/unit/test_rotate_known_issues_entries.py tests/unit/test_rotate_chunks_glyph.py -s -v --tb=short'"`

**Traps already paid for.** Never run the tool against the live `PROJECT_PLAN.md`, `known-issues.md` or `docs/planning/chunks/` — every fixture is a copy; the pinned control commits stay pinned (`15bdce61`, `KI_PINNED_SHA`), never `HEAD:`; `-k a or b` splits inside the quoted container command (do not use `-k`); pytest prints need `-s`.

**Scope.** Test data only. Not claimed: anything about future rotations beyond "closed headings only, and the four stay findable".

**Status it can move.** The known-issues entry "🟡 OPEN 2026-10-12 (found by `OPS-66` step 2 …)" is **retired in the same commit** (✅ RETIRED with the W1 log). `OPS-66` stays ✅; its row gains nothing (≤ 2 KB) — the step-3 sentence goes in this file's outcome. No rule (j) amendment.

**What a negative result means.** Any other red (A1, A3, NC, the glyph module) is a tool regression since `864b90a5`: stop, revert the test edit, mark item 393 BLOCKED with the log, leave the known-issues entry open. Never edit the tool in this item. No `PREDICTED RED` line.

**Executor spawn sentence.** "The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped."

**Outcome (2026-10-12, 09:00 CDT slot second item, run by the slot itself, not delegated): DONE.** Edit as briefed, tests only, and the tool is not touched. A small `_new_minus_old` helper now carries the per-file heading scan. W1 `20261012T141002Z_OPS-66-s3-W1.log` (Status 0, 2 s, 14 passed):
- (A2-live): live 4 = the four, read on `known-issues.md@32fabbec`.
- (A2-closed): holds on both files; the archive's new − old is 7.
- (A2-four): holds.
- First control: asserted non-empty, 3 beyond the four (the `ANS-7` 10-06, `WF-17` 10-08 and `OPS-69` 10-08 headings, as on `fc95e41b`).
- Second control: the brief's `## RETIRED-NOT 2026-10-12 planted` is **not admitted** by `KI_ENTRY_RE`, which needs `RETIRED\s+` (printed: admitted=False). The next form tried, `## RETIRED 2026-10-12 PARTIAL planted`, is admitted and not closed, and the closed check catches it on the archive copy (asserted).

**Count disclosed:** the total is 14 (8 + 6), not the 8 + 7 the anchors predicted. `test_rotate_chunks_glyph.py` has six nodes (g2 is parametrized ×2), and item 384's own W1 also totalled 14 (13 passed, 1 failed), so its "7 passed" was a miscount. The known-issues entry is retired in the same commit.

## §7 row moved verbatim 2026-10-12 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/OPS-66.md`)* **The known-issues rotation finds undated RETIRED headings** — opened 2026-10-01 03:00 review from the 2026-09-30 weekly's flag 8. `rotate_known_issues` takes an entry to be a `###` heading or a `## YYYY-MM-DD` one (`scripts/maintenance/rotate_plan_archive.py:182`, the local `entry_re` — not `_ki_closed`, which already accepts them); four retired entries open `## RETIRED <date>` / `## ✅ RETIRED <date>` and never rotate. **Done-when (§4):** the entry pattern becomes a module constant `KI_ENTRY_RE` that also admits `## [✅ ]RETIRED\|RESOLVED\|FIXED\|CLOSED <date>`; a unit test shows (a) the four headings, as verbatim literals, are entries and are closed; (b) on the live file and `known-issues-archive.md` every heading the old pattern matched is still matched, new − old is exactly the closed undated-retired headings, and no section header becomes an entry; (c) the pinned pre-change tool's dry-run on a copy of the live file reports none of the four and the new one all four; no rotation is run (the weekly's) ⇒ ✅ (§9 item 139). **✅ 2026-10-01 (12:00 slot):** (a) 4/4; (b) live 29 → 33 (new − old = the four), archive 77 → 77, no section header admitted; (c) pinned `e5af2cae` reports 0, new 4; dry run lists exactly the four (`20261001T170610Z_OPS-66-W1.log`, 8 passed, 2 s; `20261001T170616Z_OPS-66-W2-dry.log`). The test falls back to `known-issues.md` at `32fabbec` once the weekly has moved the four (disclosed deviation). *Audited PASS 2026-10-01 18:00 review (`…170610Z_OPS-66-W1.log:36–39,47–48,53–57,60,63–64`; `…170616Z_OPS-66-W2-dry.log:35–39,42–43`). Caveat: the regex also admits `## ✅ <date>` with no keyword; harmless today, since new − old is exactly the four.* Step 2 (§9 item 384, 2026-10-11): `chunks` strips a state line's leading row glyph (refuses a different one); the pinned pre-change tool doubles it, the new one does not; no-glyph output byte-identical (`20261012T035822Z_OPS-66-s2-W1.log`). | ✅ 
