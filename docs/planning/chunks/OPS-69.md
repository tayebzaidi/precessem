# OPS-69 — docrefs dead=2: the ans9 SPEC names private results files

*Opened 2026-10-08 by the 18:00 daily review, under §9 rule (5) lane 1 (an open known-issues entry whose own "Resolves it" a slot can execute). Documentation only: no solve, no `src/` change, no checker change.*

**Proposed §7 row.** `OPS-69` | docrefs `dead=2`: the `ans:9` SPEC names gitignored private results files | ⬜ The known-issues entry "2026-10-08 — 🟡 OPEN (found by `EX-89` …) — the examples docrefs census reads `dead=2`, `exit=1` on clean `main`" (`known-issues.md` ≈ :816) traces both dead references to `examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md` :36 / :38, which the operator's commit `4d4d93c5` added. Those lines name two gitignored JSON files in `aed_results/` as filenames. **Done-when:** those SPEC lines (and the script name on :39, the same class) are reworded so that they form no filename token. The meaning is kept, no number is added, and the change is disclosed in the SPEC, the commit and the journal. Then the census reads `dead=0` and `exit != 1` (W1, after W0 read `dead=2 … exit=1`), the checker's module goes green with its planted-dead fixtures still exiting 1 (W2), and the entry is re-marked RETIRED in place (W3 dry-run lists it closed) ⇒ ✅. | smoke

### §9 item 336 — 2026-10-08 (ops / docs)

**Executor:** implementer. The slot may execute it directly. **Tier:** smoke, no MPI. There are four bare harness windows, each a few seconds: the census measured `Elapsed (s): 1` (`20261008T173603Z_EX-89-item332-post-docrefs.log`, footer), and the checker's unit module took 12.25 s for 24 tests together with a second module (`20261007T003811Z_ANS-7-item288-pin.log:488`). **Cost:** ≈ 4 measured / 20 predicted min. In the ops / docs class, `OPS-68` measured 5 min against 25 min predicted.

**Anchor.** The census's `RESULT:` line, from `scripts/testing/check_example_doc_references.py` (exit-code contract `OPS-19`: 0 clean, 1 hard violation, 2 stale-only).
- **Before:** `20261008T173603Z_EX-89-item332-post-docrefs.log:98` reads `RESULT: dead=2 guide=0 stale=56 stale_severity=report exit=1`. The two dead references are listed at `:36–38`: `ans9_aed_results.json` [`SPEC.md:36`] and `ans9_corrected_rows.json` [`SPEC.md:38`]. The last clean reading is `20261008T004224Z_EX-102-post-docrefs.log:96` (`dead=0 … exit=2`).
- **After:** `dead=0`, `guide=0`, `exit != 1`. Expect `exit=2`, because `stale=56` is staleness only and is not this item's to clear. `stale` must not move because of this edit: compare W0 and W1 within the same hour.
- **Second identity (W2):** `tests/unit/test_doc_reference_exit_codes.py::test_current_tree_is_free_of_hard_violations_and_obeys_the_contract` asserts `counts["dead"] == 0` on the committed tree (`:124`). That test is red on `main` by construction now and must go green.

**What counts as a reference.** The checker scans every `*.md` under `examples/` except `*_private.md` (`:446–448`). On each raw line, including HTML comments, it matches `REFERENCE_RE` (`:120–122`), a filename-shaped token ending in `.py` / `.xdmf` / `.h5` / `.bp` / `.csv` / `.json` / `.png` / `.msh`. A `.py` resolves if any file with that basename exists under `src/ examples/ scripts/ tests/` **on disk** (`repo_basenames`, `:266–272`). An artifact resolves if git tracks it or a fresh run wrote it. A directory name such as `` `aed_results/` `` or the word "JSON" is not a reference.

**The edit.** `SPEC.md` is an operator-authored benchmark spec, so change the fewest words possible and disclose the change.
1. **:36–39.** Replace the two JSON filenames and the script path with descriptions, for example: "(i) the raw harvest (the private results JSON in the gitignored `aed_results/`) is superseded for rows 2 / 3 / 6 by the operator's offline correction pass (its corrected-rows JSON, beside it in `aed_results/`; script: the operator's offline processing script in the gitignored `aed/`; it rescales …". Leave the rest of the parenthetical, and every other line, byte-identical.
   - **Why :39 is in scope:** `aed/ans9_private_process_offline.py` resolves on this box only because the operator's untracked `aed/` directory is on disk. It is excluded by `.git/info/exclude:9`, not committed. On a fresh clone it is a third dead reference that the census here cannot see.
2. **One disclosure line** directly under the bullet, in the SPEC's own indentation: `*(Wording of (i) changed 2026-10-xx by `OPS-69`: the private file names were replaced by descriptions so the examples docrefs census resolves this page; the files and their meaning are unchanged.)*` This line must contain no filename token.
3. **Check before W1:** `git diff` on `SPEC.md` shows only those lines. The new text has no `.json` / `.py` token and contains no digit that was not already on the line.

**The windows.** Each harness call is the whole Bash command: bare, repo-relative, nothing chained. Use Bash `timeout: 660000`.
- `OPS-69-W0`, **before** any edit, is the exact command from `20261008T173603Z_EX-89-item332-post-docrefs.log:12` with the ID changed:
  `scripts/testing/run_and_log.sh OPS-69-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
- `OPS-69-W1`, **after** the SPEC edit, is the same command under `OPS-69-W1`.
- `OPS-69-W2` runs the checker's module, which contains both the tree anchor and the planted-dead controls. The command is the form from `20260929T094710Z_ANS-7-item69-W3-pin-rerun.log:12`:
  `scripts/testing/run_and_log.sh OPS-69-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 60 python3 -m pytest tests/unit/test_doc_reference_exit_codes.py -s -v --tb=short'"`
- `OPS-69-W3` runs **after** the heading re-mark. It is the exact command from `20261008T005620Z_OPS-68-W1.log:12`, with the date moved to the slot's date:
  `scripts/testing/run_and_log.sh OPS-69-W3 "timeout -k 30 30 python3 scripts/maintenance/rotate_plan_archive.py known-issues <slot date> --dry-run"`. The entry must appear in its closed list. **Never run the rotation itself** (`OPS-66`: rotation belongs to the weekly).

PREDICTED RED: OPS-69-W0 census `dead` range [2, 3] (`exit=1`; 3 only if the tree lacks the operator's `aed/` script)

**Negative control.**
- **W0 is the control for the fix.** On the unedited tree it must read `dead=2 … exit=1`, naming exactly the two SPEC sites. If W0 already reads `dead=0`, the premise is wrong (someone fixed it): stop and report.
- **The checker's own control** is `test_dead_reference_still_exits_one_after_the_split[dead-artifact|dead-script]` (`test_doc_reference_exit_codes.py:161–178`). Each case plants one dead reference in a temp-dir fixture and asserts `dead == 1` and exit 1. In W2 both cases must pass beside the now-green tree anchor. That separates "the references became resolvable or disappeared" from "the checker stopped seeing dead references". No checker change is allowed, so these cases guard against an accidental one.

**Traps.**
- **Never open `aed_results/`, `COMPARISON_private.md` or `docs/private/`.** The edit needs only the tracked SPEC. No AED value goes anywhere.
- Do not add the names to `ALLOWLIST` (`:126`). The review chose the reword: an allowlist entry would publish private file names inside a tracked checker, and an `.json` allowlist entry still leaves the `.py` on a fresh clone.
- An HTML comment does not hide a filename from the census, because it scans raw lines. The disclosure line must not contain the old names.
- `stale` can drift between windows as artifacts age. Run W0 and W1 within the same few minutes and read only `dead` / `guide` / `exit` as the anchor.
- Write `≈`, never `~`. Bare harness calls only. `python3` on the host runs only through the harness, as in item 307's W3 form.

**Scope.** This item changes only `SPEC.md` :36–39 plus one disclosure line, the known-issues heading, and the `OPS-69` row. It does not touch the 56 stale artifacts (their examples' re-runs belong to other items), the checker, or the `ANS-9` adjudication: the weekly reads the private files under their real names. It does not move `ANS-9`'s status.

**Negative result.** If W1 still reads `dead >= 1`, report the remaining reference and its site from the W1 log and stop. Do not widen the reword beyond the named lines without a review ruling. If W2's tree anchor goes green but a planted-dead case goes red, the checker regressed independently: file a known-issues entry and do not retire this one. If the dry-run tool errors in W3, retire anyway by heading (W1 + W2 carry the evidence) and record the tool error in the journal.

**Retirement mechanics** follow `known-issues.md` ≈ :728–734. Re-mark the heading **in place** as `✅ RETIRED <date> (OPS-69, §9 item 336) — was 🟡 OPEN 2026-10-08 (found by EX-89 …)`, keep the body verbatim, and add one retired-by blockquote line citing W1's `RESULT:` line and W2's footer. The weekly's `rotate_plan_archive.py known-issues` moves the body verbatim to `known-issues-archive.md` and leaves the index line. Commit the SPEC edit, the four logs, the heading re-mark, the §7 row and the journal entry together.

**Status it can move.** `OPS-69` ⬜ → ✅. The known-issues entry OPEN → RETIRED, ready for the weekly's rotation. It also makes the "docrefs `exit != 1`" done-when clause of every `EX-*` item meetable again (entry, "Consequence").

**Review rulings (2026-10-08 18:00 review) — these outrank the brief above.** (Q1) Retirement follows `known-issues.md`'s own convention: re-mark the heading RETIRED in place; the weekly's rotation moves the body. (Q2) SPEC :39 (`aed/ans9_private_process_offline.py`) is in scope — a fresh clone would read it dead. (Q3) No disclosure line inside the operator's SPEC; the change is disclosed in the commit message and the slot journal only. (Q4) No allowlist route: it would put private filenames in a tracked checker. The reword names no file and adds no number.

**Outcome: ✅ DONE 2026-10-09 02:15Z** (2026-10-08 21:00 CDT slot, take-next after item 345, attempt 1). Four windows, first run each, no checker change. W0 control (`20261009T021416Z_OPS-69-W0.log:98`, 1 s): `RESULT: dead=2 guide=0 stale=56 stale_severity=report exit=1`, naming exactly `SPEC.md:36` / `:38`. Edit: SPEC (i) :36–39 → "the raw harvest (the private results JSON in the gitignored `aed_results/`) … (its corrected-rows JSON, beside it in `aed_results/`; script: the operator's offline processing script in the gitignored `aed/`; it rescales …" — no `.json` / `.py` token, no new digit, rest of the file byte-identical, no disclosure line in the SPEC (ruling Q3). W1 (`20261009T021429Z_OPS-69-W1.log:96`, 2 s): `RESULT: dead=0 guide=0 stale=56 stale_severity=report exit=2` (`stale` unmoved). W2 (`20261009T021433Z_OPS-69-W2.log`, 13 s): 20 passed in 11.49 s — tree anchor green (`dead=0 … exit=2`), `[dead-artifact]` / `[dead-script]` each `dead=1 … exit=1`, passed. Known-issues entry re-marked RETIRED in place; W3 dry-run (`20261009T021502Z_OPS-69-W3.log`, Status 0) lists it as the 1 RETIRED entry. Rotation not run (weekly's).
