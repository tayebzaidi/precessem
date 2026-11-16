# OPS-76 — `femem run` on a swept file with no `tune.at_hz` ends in an uncaught `ValueError`

*History file for the §7 `OPS-76` row (`OPS-46`). Opened by the 2026-11-02 03:00 daily review, lane 1, from item 744's disclosure (`WF-17` U114) and the known-issues entry of the same date.*

**The defect.** `_run_load_result` (`src/fem_em_solver/study_yaml.py:2034`) refuses a bad `tune.at_hz` by field, but not a missing one. A multi-frequency `load_result:` with no `tune.at_hz` reaches `Study.tune` (`study.py:1232`), which raises a plain `ValueError`. `cli.py:209` catches only `StudyConfigError`, so `femem run` ends in a traceback instead of the rc 2 by-field refusal that item 730's three refusals give. The page says the API raises `ValueError`; it says nothing false about the CLI, but the CLI's route is the one the user meets.

**Done-when.** `femem run` on that file exits 2 and prints the field `tune.at_hz`, on every rank's return code, with no output written and no solve; the API's `ValueError` naming `at_hz` is unchanged; the known-issues entry is retired in the same commit.

### §9 item 754 — 2026-11-02 (lane 1, known-issues; item 744's disclosure)

**Why.** Item 744's executor found it and correctly left it: "A one-line `src` fix is not licensed here and is the review's call." This is that call. It is a defect in a landed refusal, not a new step, so it does not count on chain U's cap.

**Anchor.** New node(s) in `tests/package/test_wf1_yaml_load_result_sweep.py` if that module exists, otherwise the module that item 730 landed (grep `tune.at_hz` under `tests/package/`). On item 730's fixture YAML with the `at_hz` line removed:
- **(Rc)** `cli.main(["run", path])` returns 2 at `-n 2` on both ranks, and rank 0's stderr names `tune.at_hz`.
- **(S0)** 0 mesh builds and 0 solves, by the raising wrapper those nodes already use.
- **(O)** no summary file is written.
- **(API=)** `Study.load_result` then `tune(1)` still raises `ValueError` matching `at_hz`, so the fix moved no behaviour of the API. `test_wf17_loaded_tune_tutorial.py:434` is the imported pin and stays unedited.

**Negative control.** A single-frequency `load_result:` with no `tune.at_hz` still tunes, bit-equal to the record that node already pins: the check must not fire on the common case.

**Execution.** In `_run_load_result`, after `loaded = …`: when `len(loaded.frequencies_hz) >= 2` and `"at_hz"` is not in `outputs["tune"]`, raise `StudyConfigError(path, "tune.at_hz", "a multi-frequency load_result needs at_hz")`. Every rank has read the file by then, as the comment at `:2039` says. Add one sentence to the saved-tune page's refusal paragraph and to `study_yaml.py:117–120`'s docstring. The page sentence is pinned at `test_wf17_loaded_tune_tutorial.py:102–105`, so extend that tuple in the same commit (**licensed here, one line**) and name it in the commit message. No other `src/` edit.

**Tier, ranks, wall clock.** Standard, `-n 2`, complex, `timeout -k 30 170`; the nodes it adds are no-solve, and the neighbouring module took 8 s (`20261102T015645Z_WF-17-s28-W1.log`). W-env, since `src/` moves. **Class:** ops / docs. **Cost:** ≈ 5.5 measured / 15 predicted min. **Depends on:** nothing.

**Traps already paid for.** The refusal must fire before anything writes, and every rank must raise, or `-n 2` hangs in a collective. `StudyConfigError` subclasses `ValueError` (`study_yaml.py:203`), so an existing `pytest.raises(ValueError)` on the YAML route stays green, and that is not evidence the rc changed: assert the return code. `cli.main` returns and never exits. `-k a or b` splits inside a quoted container command.

**Scope.** One refusal on one route. **Not claimed:** any change to `tune`, to what `at_hz` selects, or to the other families.

**Negative result.** A landed node that pinned the old traceback goes red: stop, keep the log, BLOCKED, and name it. No `PREDICTED RED:` line.

**Status it can move.** `OPS-76` ✅; the known-issues entry of 2026-11-02 retired.

**Outcome (2026-11-02T10:42Z, 04:30 CST slot, third item): BLOCKED — a landed node pinned the old traceback, the brief's named negative result; every new gate green.** W1 `20261102T103825Z_OPS-76-W1.log` (`-n 2`, complex, `timeout -k 30 170`, 58 s, the full `tests/environment` leading): 1 failed, 23 passed, 3 skipped per rank, Status 1 (`:223`, `:285`, `:288–289`). Fix, page sentence, the licensed `SWEPT_REFUSAL_SENTENCES` line and the new node are parked on `attempt/OPS-76-20261102T104100Z`. `main` carries no `src/` change, and the known-issues entry stays open.

**The red, named.** `tests/package/test_wf17_loaded_tune_yaml_tutorial.py::test_wf17_s28_stale_lines_corrected`, `:292–293`: `pytest.raises(ValueError, match="at_hz")` around `cli.main(["run", …])` reads `DID NOT RAISE` on both ranks (`:218–222`, `:228`, `:284`), because `cli.main` now returns 2, which is the defect fixed. The node's `:304` print reads `err_cli` too. The `FAILED` on the `test_wf17_s20b_blocks_and_refusals` line (`:159`) is the other rank's `s28` verdict: that node's own tokens read `PASSED`, and both summaries name `s28` alone (by token, the 2026-10-26 rule).

**Green on the branch, new node `test_ops76_no_at_hz_refused_by_field`** (`:120–124`). (Rc) rc (2, 2), rank 0's stderr names `tune.at_hz`. (S0) 0 calls. (O) no output directory. (API=) `tune(1)` raises a plain `ValueError` naming `at_hz`, type asserted. Control: a one-frequency file with no `at_hz` tunes, `C_mode` hex-equal to the API, differing []. The page's new sentence passed its pin. `test_wf16_load_result_sweep.py:182`'s run-route `ValueError` still holds, as `StudyConfigError` subclasses it.

**Ruling asked.** License editing `test_wf17_loaded_tune_yaml_tutorial.py:292–293` and `:304` to assert the return code, rc (2, 2) with the field in stderr, as the node now pins the behaviour this chunk removes. Then take the branch and re-run the same window at the same `-t`. `docs/guide/tuning_a_saved_result_from_yaml.md:72–73` says the run stops with a `ValueError` naming `at_hz`: still true of `run_study_file`, and the CLI's status 2 is unsaid there; the same ruling may want a sentence, which the brief did not license.

**Not done, deliberately.** The landed node was not edited: the brief licenses one line elsewhere and names this outcome as a stop. The comment at `test_wf16_load_result_sweep.py:180` ("tune's own refusal") goes stale with the fix, and it is not edited either.

**Disclosed.** The three skips are `test_wf17_s23_yaml_loaded_tune`, `test_wf17_s20b_loaded_tune_block` and `test_wf16_s17_on_the_sweep_file`, none of them on this route's refusal; the reasons were not read. `test_wf1_yaml_load_result.py`, where the brief's control record lives, was not in the window; the control compares against the API instead.

**Status moved.** None; §9 item 754 BLOCKED, unblock on that ruling. Attempt 1.

### Review ruling 2026-11-02 10:30 (item 754, attempt 1) — the pin moves; re-posed as item 758

**The red is the fix working, and the brief named this outcome.** `test_wf17_s28_stale_lines_corrected:292–293` asserted that `cli.main` raises `ValueError`, which is the uncaught traceback this chunk exists to remove. A test that pins a defect is not evidence against the fix. The slot stopped correctly and did not edit it.

**Licensed, in item 758 only:**
1. `tests/package/test_wf17_loaded_tune_yaml_tutorial.py:292–293`: `cli.main` returns; assert its return code is 2 on both ranks and rank 0's stderr names `tune.at_hz`. **This is stricter than the old pin,** which any `ValueError` passed. No band is involved.
2. `:304`'s print reads the return code in place of `err_cli`.
3. `docs/guide/tuning_a_saved_result_from_yaml.md:72–73`, one sentence, as below. A search of `tests/` finds no pin on the old text.
4. `tests/package/test_wf16_load_result_sweep.py:180`, comment only, since "tune's own refusal" is no longer the route. The `pytest.raises(ValueError, …)` under it stays, because `StudyConfigError` subclasses `ValueError`.

`:290–291`'s `run_study_file` assertion stays unedited, and so does the API pin at `test_wf17_loaded_tune_tutorial.py:434`.

### §9 item 758 — 2026-11-02 (lane 1, known-issues; `OPS-76` attempt 2 of 4; item 754 re-posed)

**Why.** A defect in a landed refusal with its fix already written and green on a branch. No new cap count on chain U.

**Anchor.** Unchanged from item 754: (Rc) rc (2, 2) and rank 0's stderr names `tune.at_hz`, (S0) 0 calls, (O) no output directory, (API=) `tune(1)` raises a plain `ValueError` naming `at_hz`, type asserted. New: (Pin) `test_wf17_s28_stale_lines_corrected` green with the return-code assertion.

**Negative control.** The one-frequency file with no `at_hz` still tunes, hex-equal to the API (the branch's node). Item 754's W1 is the other direction: the old assertion reads `DID NOT RAISE` on the fixed route.

**Execution.**
1. Take the branch's code by path, not by merge: `git checkout attempt/OPS-76-20261102T104100Z -- docs/guide/tuning_a_saved_result.md src/fem_em_solver/study_yaml.py tests/package/test_wf17_loaded_tune_tutorial.py tests/package/test_wf1_yaml_load_result_sweep.py`. No commit since the branch point touches `tests/package`, `docs/guide` or `src` (read 2026-11-02 10:30), so this drops nothing. Then delete the branch.
2. Replace `:292–293`, at the same indent inside the `with`, by three lines: `capsys.readouterr()`, then `rc_cli = cli.main(["run", str(no_yml), "--output-dir", str(d / "no_cli")])`, then `err_text = capsys.readouterr().err`. Beside `n_calls`: `rcs = tuple(comm.allgather(rc_cli))` and `named = comm.bcast("tune.at_hz" in err_text, root=0)`. Assert `rcs == (2, 2)` and `named` with the node's other asserts. `:304` prints `rcs`.
3. Guide `:72–73`: replace "without it the run stops with a `ValueError` naming `at_hz`." with "without it `run_study_file` stops with a `ValueError` naming `at_hz`, and `femem run` refuses the field `tune.at_hz` with status 2." Re-wrap to the page's width. The test module's `:34` docstring stays true.
4. Comment at `test_wf16_load_result_sweep.py:180`: "the by-field refusal, a `ValueError` subclass (`OPS-76`)".
5. Retire the known-issues entry of 2026-11-02 in the same commit and flip the §7 row.

**Tier, ranks, wall clock.** Standard, `-n 2`, complex, `timeout -k 30 170`: item 754's W1 command with `tests/package/test_wf1_yaml_load_result.py` added as the last path, the control-record module the first window left out. W1 took 58 s. If the added module takes the window past 170 s that is a size kill: split it into a second window at the same `-t`, never raise it. **Class:** ops / docs. **Cost:** ≈ 5.5 measured / 15 predicted min. **Depends on:** nothing. **Family:** attempt 2 of 4.

**Traps already paid for.** Item 754's traps stand. Assert the return code, not an exception. `capsys.readouterr()` before the call, or the module's earlier prints land in the text. Rank 0 alone may print the refusal, so `bcast` rank 0's reading. Path checkouts take code only. Count verdicts per rank by token.

**Scope.** One refusal on one route. **Not claimed:** any change to `tune` or to what `at_hz` selects.

**Negative result.** Any other red: stop, keep the log, BLOCKED, and name it. No `PREDICTED RED:` line.

**Status it can move.** `OPS-76` ✅; the known-issues entry of 2026-11-02 retired.

**Outcome (2026-11-02T18:02Z, 12:00 CST slot, first item): complete.** `20261102T180047Z_OPS-76-item758-W1.log` (`-n 2`, complex, `timeout -k 30 170`, the full `tests/environment` leading, `test_wf1_yaml_load_result.py` last): Status 0, 87 s, 26 passed, 3 skipped per rank (`:1189–1196`). (Rc) rc (2, 2), rank 0's stderr names `tune.at_hz` (`:120`). (S0) 0 calls (`:121`). (O) no output directory (`:122`). (API=) `tune(1)` a plain `ValueError` naming `at_hz` (`:123`). Control: the one-frequency file hex-equal to the API, `0x1.7ae3f5a5c44f3p-35`, differing [] (`:124`). (Pin) `test_wf17_s28_stale_lines_corrected` green: `cli.main rcs (2, 2)` (`:158`).

**Disclosed.** The branch's four paths were taken as the brief says, and the branch is deleted. The three skips are the same three as item 754's. The `s28` node's assertion sits one line after `selected == FREQS[1]`, and the print's words around `rcs` are the old sentence's. `run_study_file`'s `ValueError` at `s28 M` is now a `StudyConfigError` by type (`StudyConfigError True`, `:158`): still a `ValueError`, the page's sentence still true, and the node's `:290–291` is unedited.
