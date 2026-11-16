# OPS-68 — known-issues: retire four more live entries whose resolutions have already landed

*Opened 2026-10-07 by the 18:00 daily review, under §9 rule (5) lane 1 ("open known-issues entries whose own 'Resolves with' a slot can execute"), the same kind of work as `OPS-67` (✅ 2026-10-07, audited PASS). Found by the review's `Explore` sweep of the 24 live entries. Documentation only: no solve and no `src/` change.*

**Proposed §7 row.** `OPS-68` | known-issues: retire four more live entries whose resolutions already landed | ⬜ The `MAG-18` anchor (ii) entry (`known-issues.md` ≈ :413), the `OPS-19` docrefs exit-code contract (≈ :514), the four `OPS-17` step 2 defects (≈ :910) and the `EX-49` control-floor observation (≈ :1247) each carry their own resolution in the body but keep live headings. **Done-when:** for each entry, either (i) its own resolution condition is shown met by a cited landed log (footer Status 0, the named identity) and a ✅ §7 row, its heading is re-marked `✅ RETIRED <date> (OPS-68)` with a retired-by line naming that commit and log, and the rotation dry-run lists it as closed; or (ii) the condition is shown *not* met, and the entry gets a dated annotation saying what is still owed and stays OPEN. W1 − W0 equals the re-marked set exactly. No rotation is run (`OPS-66`: the weekly's) ⇒ ✅. | smoke

### §9 item 307 — 2026-10-07 (lane 1, known-issues housekeeping)

**Outcome (2026-10-08 00:57Z, ✅):** executed directly by the 19:30 slot. W0 (`20261008T005516Z_OPS-68-W0.log`, Status 0) listed 4 closed entries (the three `OPS-67` retirements and the 2026-09-19 `ANS-4` entry), none of the four. Every condition met, all on branch (i): (1) `MAG-18` ✅, (ii) re-registered at ≤ 1e-6, both MAG-18 record logs Status 0; (2) `OPS-19` ✅, rerun 8 passed, and `exit != 1` still in §9 standing rules (PROJECT_PLAN.md :3684); (3) `GEO-17` / `MAG-17` / `POST-5` ✅, the four fix logs 3 / 6 / 12 / 8 passed, Status 0; (4) `EX-49` ✅, closed observation, log Status 0. Headings re-marked in place with one retired-by blockquote each, bodies verbatim. W1 (`20261008T005620Z_OPS-68-W1.log`, Status 0) lists 8: W1 − W0 = {`MAG-18` (ii), docrefs exit codes, `OPS-17` four defects, `EX-49` floor}, exactly the re-marked set. No rotation run.

**Executor:** implementer (the slot may execute it directly, as the 13:30 slot did item 303). **Tier:** smoke, three harness windows of ≈ 0 s (`-n 1`, no MPI). **Cost:** ≈ 5 measured / 25 predicted min. Ops / docs class: item 303 measured 3 min against 25 predicted; this interval's measured / predicted ratio is 0.24 (67 / 280).

**Anchor.** For each entry, its own resolution text is the condition, and the evidence is a harness log already on `main`, read at its footer and at the line carrying the named identity. The sweep's citations (verify each; do not trust them):
1. **`MAG-18` anchor (ii)** (≈ :413). Condition: "Re-registering (ii) at the measured floor … is a review decision; `MAG-18` stays 🟡 until it is made." Body ≈ :459: "✅ RESOLVED 2026-08-23 (03:00 review) — (ii) re-registered at ≤ 1e-6". `MAG-18` ✅ (§7 ≈ :1666). Log `20260823T003327Z_MAG-18-record-probe.log.gz`.
2. **`check_example_doc_references.py` exit codes** (≈ :514). The body opens "Resolved 2026-08-16 (`OPS-19` step 1)"; the heading calls it "contract, not an issue". `OPS-19` ✅ (§7 ≈ :1327). Log `20260816T213312Z_OPS-19-step1-rerun.log.gz`. **The contract itself stays stated in §9's standing rules** (`exit != 1`), so retiring the entry loses no instruction; check that sentence is present before re-marking.
3. **Four defects from `OPS-17` step 2** (≈ :910). Each defect names its own fix, and the body marks each struck or RESOLVED: (1) `GEO-17` ✅ (≈ :1841), `20260820T110549Z_GEO-17-step1-final.log.gz`; (2) `MAG-17` ✅ (≈ :1665), `20260820T123823Z_MAG-17-step1-final2.log.gz`; (3) `POST-5` step 4 ✅ (≈ :2584), `20260819T201005Z_POST-5-step4-smoke-final.log.gz`; (4) `POST-5` step 1, `20260818T215117Z_POST-5-step1-negcontrol.log.gz`. **All four** must be met for branch (i); if any one is not, use branch (ii) for the entry and name the open defect.
4. **`EX-49` control floor** (≈ :1247). Its body: "Nothing further needed; this is a closed observation, not an open chunk." The 2026-09-05 18:00 ruling sits in the body (≈ :1257). `EX-49` ✅ (≈ :2988). Log `20260905T201151Z_EX-49.log.gz`.

Logs are gzipped: read them with `zcat` (not through the harness; reading is not a window).

**The windows** (bare harness calls, repo-relative, nothing chained; Bash `timeout: 660000`):
- `OPS-68-W0`, **before** any edit: `scripts/testing/run_and_log.sh OPS-68-W0 "timeout -k 30 30 python3 scripts/maintenance/rotate_plan_archive.py known-issues 2026-10-07 --dry-run"` — the exact command of `20261007T185849Z_OPS-67-W1.log:12`, with the ID changed.
- `OPS-68-W1`, **after** the heading edits: the same command under `OPS-68-W1`.
The comparison of the two closed lists is the item's quantitative identity: W1 − W0 equals exactly the entries re-marked RETIRED, an exact set equality printed in the journal.

**Negative control.** W0. Before the edits, none of the four appears in the dry-run's closed list (W0 should list the three `OPS-67` retirements and the 2026-09-19 `ANS-4` entry, four in all, per `OPS-67-W1.log:35–39`). If W0 already lists one of the four, the premise is wrong for it: stop and report.

**Traps.**
- The harness call is the **whole** Bash command: bare, repo-relative, no `cd`, nothing chained or piped. Never delete a log.
- `python3` on the host is denied except through the harness (item 303's windows ran this exact form, Status 0).
- **Never run the rotation itself** (no `--dry-run` omitted); `OPS-66` made rotation the weekly's.
- Re-mark headings **in place**; keep every body verbatim; add one retired-by blockquote line per entry, as `OPS-67` did (`known-issues.md` ≈ :503–505). Delete nothing.
- A retired-by line states only what its log shows (the `OPS-67` audit caught one line that claimed a red on a green window). Cite a line number for each claim.
- Write `≈`, never `~`.

**Scope.** No code, no test, no solve. No §7 glyph other than `OPS-68`'s moves. The entries the sweep found ambiguous (≈ :51, :69, :122, :549, :646, :707, :819, :836, :846, :1275) are **not** touched: their conditions wait on a review or weekly decision, and the duplicate reading of ≈ :122 against ≈ :51 is the weekly's rotation question.

**Negative result.** If no condition is met, the item still closes ✅ on branch (ii): four dated annotations, W0 = W1, and the journal says the sweep's premise was wrong. If the dry-run tool errors, park on `attempt/OPS-68-…`, mark the §9 item BLOCKED in the same commit (rule (d)), and report.

**Status it can move.** `OPS-68` ⬜ → ✅. Up to four known-issues entries change OPEN → RETIRED and become ready for the weekly's rotation.
