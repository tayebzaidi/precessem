# POST-4 — history file

*Created by the 2026-10-23 10:30 daily review for a §9 brief (`OPS-46`). The §7 `POST-4` row (✅, closed 2026-08-12) carries the chunk's status; its step 4 measurement is the known-issues entry "Every XDMF/VTX export of a Nédélec or DG field ships a Lagrange-P1 interpolant …" (`docs/testing/known-issues.md`, filed 2026-08-12).*

### §9 item 607 — 2026-10-23 (ops / docs; lane 1, a known-issues entry's own scope)

**Why.** The entry's "What this affects" table lists **11** P1 interpolation sites in `examples/` as of 2026-08-12, ten of them fed by a non-conforming source (N1curl / DG). The census is now 138 examples, and `grep -rn '"Lagrange", 1\|("CG", 1\|("P", 1' examples --include=*.py` finds ≈ 82 constructions (this review's count; the executor re-counts). The entry's one *evaluated* site, `mag:1`, stopped evaluating from its interpolant on 2026-10-21 (`a49b4e16`; noted in the entry by this review). The entry's open question — export DG1, or keep P1 with the caveat — is a review decision that needs the current table, not the August one. **Nothing here decides it.**

**The anchor (a census identity, measured).** A probe script `scripts/probes/post4_p1_site_census.py` (new, stdlib only, no dolfinx import) walks `examples/**/*.py`, finds every `functionspace(…, ("Lagrange"|"CG"|"P", 1 …))` construction, and prints one row per site: file:line, and the **name** of the function that is interpolated into it on the following `.interpolate(` call when one exists in the same function body. The executor then classifies each row **by reading the code**, into exactly one of:
- **(N)** source is N1curl / RT / any `ufl` curl of one;
- **(D)** source is DG of any degree;
- **(C)** source is conforming (Lagrange of any degree, a numpy/analytic expression) — safe by construction;
- **(U)** cannot tell from the file (say why).

and, orthogonally, **(x)** export only, or **(e)** the interpolant is evaluated, sampled, integrated or asserted on.

**Asserted by the script (window 1):** (i) the 08-12 table's ten file paths each still exist or are named as removed / moved, with their current line; (ii) the row count equals the independent `grep -c` total printed in the same window (the negative control — a parser that misses a spelling shows here). Window: smoke, `scripts/testing/run_and_log.sh POST-4-census "timeout -k 30 30 python3 scripts/probes/post4_p1_site_census.py"` (host `python3` through the harness, as `PLAN-INDEX` runs).

**Deliverable.** A dated sub-section appended to the known-issues entry, **"Re-census 2026-10-23 (`POST-4` §9 item 607)"**: the counts per class (N / D / C / U × x / e), the full row list or a pointer to the probe's log, and — the part that matters — **every (N)- or (D)-sourced site in class (e)**, by file:line, with what it prints. The August table stays as history above it (do not edit or strike it). Commit the probe script, the log, the entry and a one-line attempts journal.

**Executor:** implementer (no compute; one smoke harness window). **Tier:** smoke, no MPI.
- **Class:** ops / docs. **Cost:** ≈ 3.0 measured (the class median; this item is larger than the class's typical guide-text fix) / ≈ 40 predicted slot-minutes.
- **Depends on:** nothing. **Family:** attempt 1 of 4. **Daytime:** no.

**Traps already paid for.** The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. No `docker` call is needed. Never delete a log. Write `≈`, never a tilde, in the entry. Do not touch any example: a site found in class (e) with an (N)/(D) source is **reported, not fixed** — fixing it is a chunk the weekly opens.

**Scope.** Not claimed: any export change; the DG1-vs-P1 decision; any statement about ParaView renderings beyond the entry's existing one; any gate. The entry stays open.

**Negative result.** If the parser's count and the grep count differ, report both and the missing spellings, commit the probe and log, and stop without the sub-section (rule (d): mark item 607 BLOCKED). If an (N)/(D)-sourced (e) site is found whose printout a guide quotes as a number, that is the finding: list it in the sub-section and in the journal's "Hypothesis for next", nothing more.

No `PREDICTED RED:` line.

**Status it can move.** The known-issues entry's scope: its affected-site list re-dated to 2026-10-23 with counts by class, and a named list of evaluated-from-P1 sites (zero, or each one by file:line) — the input the weekly needs for the entry's open DG1-vs-P1 decision. No §7 glyph moves.

**Outcome (2026-10-23 17:00 slot, implementer).** Done: `20261023T171007Z_POST-4-census.log` (smoke, 2 s, exit 0) — 82 P1 constructions in 151 example files, AST count = line-regex control (82 = 82; the brief's grep reads 83, the extra line is `mri/01:45` docstring prose); all ten August paths present. Classes: (N) 54 [52 x / 2 e], (D) 6 [2 x / 4 e], (C) 22 [21 x / 1 e], (U) 0. The six (N)/(D)×(e) sites — `ans2:336`, `th:5:129`, `mag:6:192`, `mag:5:218`, `mag:1:251`, `mag:2:315` (the last two a VTX file-fidelity identity) — are listed with their printouts in the known-issues entry's "Re-census 2026-10-23" sub-section; none fixed; whether a guide quotes any of those printouts was not audited.

### §9 item 636 — 2026-10-25 (§10 lane 1, R2: item 607's six P1-export sites against the guides)

*Written by the 2026-10-25 03:00 daily review from the 10-24 weekly's ruling (`PROJECT_PLAN.md` § 10, "Known-issues 2026-09-09 … R2"): audit whether any guide quotes those printouts as a figure, and where one does, re-quote it from the native field under the owner rule, with the P1 value kept as a labelled record. The picture exports stay P1. No band, no gate edit.*

**The audit (done by this review, read-only; `docs/guide/` quotes none of them).**

| Site | Guide quote | Disposition |
|---|---|---|
| `ans2:336` (`examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/02_birdcage_coil_driven_sar_10MHz.py:336`) | `02_…md:40–51` quotes the CG1 row as a "labelled record, never published" | **none**: already labelled as asked |
| `th:5:129` (`examples/time_harmonic/05_resonance_guard_sweep.py:129`) | `05_resonance_guard_sweep.md:47–48`: "Exported peak \|E\| … 6.0531e+03 V/m", quiet 6.1951e+02, factor 9.77 | **re-quote**: a figure from the P1 export |
| `mag:6:192` (`examples/magnetostatics/06_h_convergence_rate.py:192`) | `06_h_convergence_rate.md:45`, `:130`: "Error of the *exported* CG1 field 17.1451%" | **none**: the page's point is the export's own error, labelled "exported", and it is owner-sampled since `POST-8` (`:172–178`) |
| `mag:5:218` (`examples/magnetostatics/05_gauge_cross_check.py:218`) | `05_gauge_cross_check.md:39`, `:109`: 0.0033% | **re-quote**, if the printout is the P1-export cross-check (W0 decides from the script) |
| `mag:1:251`, `mag:2:315` | `01_straight_wire.md:80–81`, `02_circular_loop.md:71–72`: in-memory = read-back max\|B\| | **none**: a VTX file-fidelity identity, which is about the file by design |

**The item: two examples, `th:5` and `mag:5`.**
1. **W0 (no solve):** read both scripts. Confirm that each quoted figure is computed from the P1 export (the known-issues "Re-census 2026-10-23" table, `known-issues.md:749–758`). If a figure turns out to be computed from the native field already, mark it "none" and say so.
2. **For each confirmed site,** add one printed line computing the same quantity from the native field through `evaluate_vector_field_parallel(…, owner="lowest_original_cell")` (`src/fem_em_solver/post/evaluation.py:224–249`). The P1 line stays, relabelled "P1 export (picture), record". No assert is added or changed: rule (a) — the examples' existing asserts stay byte-identical.
3. **Run each example** through `./run_examples.sh` as its own harness window, then the docrefs and setup-figure census windows, as every example item does (rule (i)).
4. **Re-quote** each guide line from the native-field print, with the P1 value beside it labelled as the export's.

**Anchor.** Each example's existing asserts re-run green unedited; the window is the gate (rule (a)). The re-quoted guide number equals the printed native-field string, and the P1 string equals the old guide's.

**Control.** The old P1 line reproduces the old guide number to its printed digits, so the export path did not move.

**Tier:** standard (both examples are small). **Class:** lane-1 revalidation. **Cost:** ≈ 11.4 measured / 35 predicted slot-minutes. **Executor:** `example-runner` (example scripts plus same-stem guides, census before and after). Spawn it foreground, with the no-spawn and no-edit-after-last-run rules (dashboard item 9). **Depends on:** nothing. **Family:** R2, attempt 1 of 4.

**Traps.**
- `.eval` bytes depend on the batch (item 537): evaluate the native field in one call per quantity.
- The owner rule's keyword is `owner="lowest_original_cell"`, never the legacy path.
- An example edited after its last run breaks rule (i): edit, then run, then census.
- The leak guard: a printed closed-form column can trip it (`EX-142`). If it does, move the log to `logs/leak-held/` and report.

**Scope boundary.** Two guide re-quotes. Not closed: the picture exports' P1-vs-DG1 question (the operator's, dashboard item 13), the other four sites (dispositioned above), any band. `POST-4` stays ✅.

**A negative result means:** an existing assert red is an example defect, so stop and report. A native-field number differing materially from P1 is the expected reading (that is why the re-quote is owed): print both.

**Outcome (2026-10-25 13:30 CDT slot, `example-runner`) — BLOCKED, attempt 1 of 4.** W0 confirmed both quoted figures come from the P1 export. Native prints were added; the examples' existing asserts are green and unedited (`20261025T183230Z_POST-4.log` th:5 32 s, exit 0; `20261025T183305Z_POST-4.log` mag:5 6 s, exit 0; `20261025T183149Z_POST-4.log` a 35 s exit-1 helper crash, `index_map` is a method, fixed before the green run). Readings: th:5 native vertex peak \|E\| 5.8778e+03 / 5.9984e+02 V/m (factor 9.80), P1 6.0350e+03 / 6.1719e+02 (9.78), against the old guide's 6.0531e+03 / 6.1951e+02 (9.77). mag:5 native DG volume L2 0.0033 %, P1 0.0040 %, against the old guide's 0.0033 %. **The Control missed at both sites:** the P1 line does not reproduce the old guide digits, so either the export path moved or the August guide figures predate the current image. The control carries no asserted / predicted label, so under rule (e) the slot did not decide it. Censuses: docrefs `dead=0 exit=2` both sides, stale 147 → 146; setup figures `broken=0` both sides. The code, guides and known-issues disposition are parked as `5b55cd7c` on `attempt/POST-4-R2-20261025T183400Z`. Deviation: the pre-census ran after the two script edits and before any run.

**Status it can move.** Item 607's open clause ("whether a guide quotes any of those printouts was not audited") closes: the audit table above plus two re-quotes. The known-issues 2026-10-23 re-census sub-section gains the disposition per site. §10 R2 → DONE.

**Outcome (2026-10-25 13:30 slot, implementer).** Done in part: both W0 reads confirmed the guide figures were P1-export quantities. `th:5` (`20261025T183230Z_POST-4.log`, 32 s, exit 0; a first window `20261025T183149Z_POST-4.log` failed on `index_map` being a method in 0.11 after the table and all anchors had printed, fixed and re-run): native vertex peaks 5.8778e+03 / 5.9984e+02 V/m, factor 9.80, against the export's 6.0350e+03 / 6.1719e+02, 9.78. `mag:5` (`20261025T183305Z_POST-4.log`, 6 s, exit 0): native DG volume L2 difference 0.0033%, export 0.0040%. Existing asserts green and untouched. **The pre-registered control missed at both sites** (the export line does not reproduce the old guide's 6.0531e+03 / 0.0033%): the old guide figures predate the current image. Guides re-quoted from the native prints with the export's values labelled beside them. Census: docrefs pre `dead=0 guide=0 stale=147 exit=2`, predicted unchanged, measured `stale=146` (prediction missed by one stale entry; dead/guide/exit as predicted); setup figures `broken=0` before and after. Deviations: the pre-census ran after the script edits (guides untouched at that point); `mag:5` native number is an assembled volume integral, not `evaluate_vector_field_parallel`; `mag:5` is not selectable as `mag:5` in `run_examples.sh` (selector `-e 5`, as the runner's group is the bare number).

**Review ruling (2026-10-25 18:00 daily review) — the control holds on the current image; `5b55cd7c` landed, item 636 DONE, R2 attempt 1 of 4 closes it.** The control's purpose was "the export path did not move". Its comparand, the guides' digits, dates from 2026-08-10, before the 2026-08-22 image change (`OPS-18`). Against the current image's own record, the P1 lines reproduce exactly. **th:5:** 6.0531e+03 / 6.1951e+02 (factor 9.77) at `20260810T033313Z_EX-8-gate.log.gz:66`; 6.0350e+03 / 6.1719e+02 (9.78) at `20260825T050232Z_EX-30-th-run-5to6.log.gz:66`, which is today's print digit for digit (`20261025T183230Z_POST-4.log`). **mag:5:** the export reads 0.0033 % on 08-10 and 08-20 (`20260810T110311Z_EX-10-run1.log.gz:442`, `20260820T003532Z_EX-22-mag-56.log.gz:442`), then 0.0040 % from 08-25 (`20260825T140737Z_EX-30-root-run-mag2to5.log.gz:2227`, `20260916T095433Z_EX-57-mag5-control.log.gz:229`, `20261020T214139Z_POST-8-step3b-ii-e2d-e1-real-n2.log:662`). Today's 0.0040 % (`20261025T183305Z_POST-4.log:231`) is the same figure. So the export path has not moved since the image change, and the guides were stale records. The native 0.0033 % equalling the old export figure is a coincidence of rounding at four digits. No pathologist ruling was needed: the readings are direct greps of footered logs, and no assert or band is involved. Landed as its own commit, with main's §9 and test-results rows kept (the logs and rows were already on main from `e8eac1eb`).
