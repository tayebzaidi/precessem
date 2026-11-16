# OPS-67 — known-issues: retire three live entries whose resolutions have already landed

*Opened 2026-10-07 by the 10:30 daily review, under §9 rule (5) lane 1 ("open known-issues entries whose own 'Resolves with' a slot can execute"). Found by the review's supply sweep. Documentation only: no solve and no `src/` change.*

**Proposed §7 row.** `OPS-67` | known-issues: retire three live entries whose resolutions already landed | ⬜ The `PORT-13` transverse-ring-sheet entry (`known-issues.md` ≈ :503), the `_DummyComm` / `allgather` entry (≈ :1224) and the `test_coil_loading_degree2.py` runtime entry (≈ :730) are still in the live file, although `GEO-26` ✅, `OPS-28` ✅ and `TH-13` step 3a‴ (archive ≈ :3280) appear to have resolved them. **Done-when:** for each entry, either (i) its own "Resolves with" condition is shown met by a cited landed log (Status 0, the named identity), its heading is re-marked `✅ RETIRED <date>` with a retired-by line naming that commit and log, and the rotation dry-run lists it as closed; or (ii) the condition is shown *not* met, and the entry gets a dated annotation saying what is still owed and stays OPEN. No rotation is run (that is the weekly's, `OPS-66`) ⇒ ✅. | smoke

### §9 item 303 — 2026-10-07 (lane 1, known-issues housekeeping)

**Executor:** implementer. **Tier:** smoke, one harness window (`-n 1`, no MPI). **Cost:** ≈ 5 measured / 25 predicted min. This uses the ops / docs class (uncalibrated, n = 4): 25 predicted × the 03:00 → 10:30 interval's measured / predicted ratio of 0.18 (52 / 285).

**Anchor.** For each of the three entries, the entry's own **"Resolves with"** text is the condition. Its evidence is a harness log that already exists on `main`, read at its footer (Status 0) and at the line carrying the named identity:
1. **`PORT-13` ring sheets** (2026-09-03, ≈ :503). It retires with "`GEO-26` step 2's 16-leg record (the control `PORT-13` step 1 re-opens on), which must point its retired-by row here". Find `GEO-26` step 2's 16-leg record log in `docs/planning/chunks/GEO-26.md`, and check that `PORT-13`'s row (§7, ≈ :2702) is ✅ on longitudinal sheets.
2. **`_DummyComm` / `allgather`** (2026-08-17, ≈ :1224). Read its own Resolves-with. `OPS-28` (§7 ≈ :1335, ✅ 2026-08-28) added `allgather`, **but its row says the S-matrix assertion then reached its check "and is red there"**. So check whether that red was itself resolved and by what, and whether the archive's retired entry (`known-issues-archive.md` ≈ :588) is **the same entry**. If the two are the same, the live copy is a stale duplicate. If they differ, the live entry is judged on its own condition. If its condition is not met, use branch (ii).
3. **`test_coil_loading_degree2.py` runtime** (≈ :730). The archive carries "✅ RETIRED 2026-09-01 … by `TH-13` step 3a‴" (≈ :3280). Decide whether the live text is a stale copy of that entry or a distinct one, with a `diff` of the two bodies printed in the window.

**The window** is a bare harness call, repo-relative, with nothing chained:
`scripts/testing/run_and_log.sh OPS-67-W1 "timeout -k 30 30 python3 scripts/maintenance/rotate_plan_archive.py known-issues 2026-10-07 --dry-run"`
Check the exact argv form in `rotate_plan_archive.py` (≈ :810–820) before running. Run it **after** the heading edits. Its listing of closed entries must contain exactly the entries you marked RETIRED (no more, no fewer) plus whatever it listed before. Also run it once **before** the edits (`OPS-67-W0`): the difference W1 − W0 is the asserted set. The comparison of the two lists is the item's quantitative identity: an exact set equality, printed in the journal.

**Negative control.** W0 is the control. Before the edits, none of the three appears in the dry-run's closed list (their headings carry no RETIRED marker). If W0 already lists one, the premise is wrong: stop and report.

**Traps.**
- The harness call is the **whole** Bash command: bare, repo-relative, no `cd`, nothing chained or piped. Pass `timeout: 660000`. Never delete a log.
- `python3` on the host is denied except through the harness. This window needs no container (if the harness denies a host `python3`, use `docker compose exec -T fem-em-solver bash -lc 'cd /workspace && …'` as the whole harness argument).
- **Never run the rotation itself** (no `--dry-run` omitted); `OPS-66` made rotation the weekly's.
- Known-issues discipline: an entry leaves only with the commit that fixes it. Here, the fixes landed earlier. This item **re-marks** the headings (in place, the text kept verbatim) and adds retired-by lines. It deletes a live body only when it is shown byte-identical to an archived one, and says so.
- Write `≈`, never `~`.

**Scope.** No code, no test, no solve. No §7 glyph other than `OPS-67`'s moves. If an entry's condition is not met, nothing moves for it except the dated annotation.

**Negative result.** If none of the three conditions is met, the item still closes ✅ on branch (ii): three dated annotations, W0 = W1. That is a finding (the sweep's premise was wrong), stated in the journal. If the dry-run tool errors, park on `attempt/OPS-67-…`, mark §9 BLOCKED, and report.

**Status it can move.** `OPS-67` ⬜ → ✅. Up to three known-issues entries change OPEN → RETIRED and become ready for the weekly's rotation. `PORT-13`'s and `GEO-26`'s rows gain the retired-by pointer the entry asked for, if it is missing.

**Outcome (2026-10-07, 13:30 CDT slot, third item, executed by the slot directly): ✅ DONE.** Control W0 (`20261007T185729Z_OPS-67-W0.log`, Status 0) listed one closed entry (the 2026-09-19 `ANS-4` entry) and none of the three — premise holds. All three conditions met (branch (i)); each heading re-marked `✅ RETIRED 2026-10-07 (OPS-67)` in place with a retired-by line, bodies verbatim:
1. **`PORT-13` ring sheets:** `GEO-26` step 1 `41cd5c84`, step 2's 16-leg record `1ad8ba3c` (`20260903T170351Z_GEO-26.log`, the terminal-class STOP), step 3's band green `2445b9ed` (`20260904T034441Z_GEO-26.log`, 4 passed, Status 0); `GEO-26` ✅ and `PORT-13` ✅.
2. **`_DummyComm` / `allgather`:** the archive's ≈ :588 entry is the same defect under a different heading (found again by `OPS-26`'s census) — so the live copy was a stale duplicate. Its own condition was met by `OPS-28` `fd1479ca` (double repaired, entry 3 re-dated in the same commit), bracketed by `20260828T033037Z_OPS-28-red-baseline.log` / `20260828T033055Z_OPS-28-gate.log` (no `AttributeError`; its Status 1 is entry 3's two reds alone). The module was later rewritten without the double (`PORT-8`, `5418f66b`).
3. **`test_coil_loading_degree2.py` runtime:** diff in `20261007T185802Z_OPS-67-diff.log` — the live text is the original open body, the archive (≈ :3278) holds the retirement summary; not byte-identical, so nothing deleted. Resolved by `TH-13` step 3a‴ `f494dc75` (three footered windows `20260901T033335Z` / `T033434Z` / `T034059Z_TH-13-step3a3-*.log`).

W1 (`20261007T185849Z_OPS-67-W1.log`, Status 0) lists 4 closed: W1 − W0 = exactly the three re-marked headings (set equality). No rotation run. The `PORT-13` / `GEO-26` §7 rows were not given a retired-by pointer (both ✅ with history files; the pointer lives in the known-issues retired-by line).

**Audit (2026-10-07 18:00 daily review): PASS** (`auditor`; W0 `20261007T185729Z_OPS-67-W0.log:35–36`, W1 `20261007T185849Z_OPS-67-W1.log:35–39`, set identity W1 − W0 = the three re-marked headings; every retired-by log exists with a footer, e.g. `20260904T034441Z_GEO-26.log.gz:26418,26421`). **Two bookkeeping notes for the weekly's rotation (`OPS-66`), no item queued:** (1) the `PORT-13` entry's Resolves-with asks the `GEO-26` record to point its retired-by row at the entry; that back-pointer sits in the entry's retired-by line instead of the §7 rows. (2) The `test_coil_loading_degree2.py` retired-by line says all three `TH-13` 3a‴ windows carried the identity red; the probe window (`20260901T033335Z_…-original-probe.log.gz:476,500`) is green — only the loaded and free windows carry it.
