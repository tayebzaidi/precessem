# MAG-14 — Helmholtz magnitude comparison in the test suite (history)

### §9 item 281 — 2026-10-06 (lane 1, (1*) re-record)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: the (1*) licence, under the 10-03 weekly's ruling in PROJECT_PLAN.md §10 ("`MAG-14`'s centre-field figure (0.671 % against the pre-0.11 0.728 %) is a (1*) re-record candidate — `record-reconciler`, lane 1"; grep `re-record candidate`). **No band, tolerance, gate constant or Done-when is set or moved.** `MAGNITUDE_REL_TOL` 0.05 and `UNIFORMITY_CV_MAX` stay as they are (`tests/validation/test_helmholtz_magnitude.py:67–74`).

**Executor: `record-reconciler`** (`.claude/agents/record-reconciler.md`). Spawn it in the foreground. Its windows run foreground too (rubric 4).

**What the record is.** The §7 row `^| \`MAG-14\` |` (PROJECT_PLAN.md:1661) says "0.728% vs closed form (1.731% before `GEO-8`); 11 s, in CI". 0.728 % was measured on the pre-0.11 image: `GEO-8`'s fix (2026-08-01, `20260801T003528Z_GEO-8-after.log.gz`), then again at `OPS-22`'s real baseline, with centre 0.728 %, mean 0.644 % and CV 0.1602 % (2026-08-19, `20260819T093105Z`, PROJECT_PLAN.md:1531). On the 0.11 image (dolfinx 0.11.0.post0, real build) the same test node prints the following, the same in all six `WF-13` installed-package windows (most recent: `20261006T040812Z_WF-15-s2-W4-wf13.log:257–260`):
`cells=52844`; `centre B_z: FEM=3.554758e-09 T  closed form=3.531057e-09 T  rel err=0.671%`; mean 0.572 %; central CV 0.2056 %; node 10.91 s.
Those six windows ran the **installed** package from `/tmp/wf13/venv`, with `PYTHONPATH` unset. No in-repo window of this module on the 0.11 image has printed the figure. W1 below is that window.

**Sites (from the readings above, not from grep):**
- (s1) PROJECT_PLAN.md:1661, the §7 row's figure: **move**.
- (s2) `examples/magnetostatics/04_helmholtz_analytic_comparison.md:81`, the guide copy of the same figure ("**0.728%** vs closed form"): **move in the same commit** (the GEO-16 rule).
- **Not sites:** the dated `OPS-22` annotation (PROJECT_PLAN.md:1531) and the `GEO-8` narrative (:1911). Both are dated measurements on the old image and stand as written. The test module's docstring (`:29–43`, 1.731 % / 53 941 cells, 2026-07-29) is pre-`GEO-8` history, and no constant in the module carries 0.728. **The test module is not edited**, so there is no rule (c) re-gate and `WF-13`'s sha256 copy check is unaffected.

**Anchors (asserted).**
- **(A) The module itself, in-repo, on the 0.11 image:** `test_helmholtz_centre_field_magnitude` green. Inside it: helper vs closed form (4/5)^1.5·μ₀I/R < 1e-12 (`AnalyticalSolutions.helmholtz_coil_field_on_axis`, `:156–157`); centre and mean rel err < `MAGNITUDE_REL_TOL`; CV < `UNIFORMITY_CV_MAX`; Im B_z ≤ 1e-12·max (`:139`).
- **(Rep) Reproduction identity, read by the executor:** W1 must print the installed windows' line character for character: `cells=52844`, `FEM=3.554758e-09`, `rel err=0.671%`, mean `0.572%`, CV `0.2056%`. Only a footered log can supply the new value (contract rule 2). If W1 and the installed runs disagree, that is the finding (see Negative result).
- **(D) The (1*) drift gate, computed on B_z and not on the error %.** Read the old FEM B_z from `20260819T093105Z`'s log (or from the `GEO-8`-after log; both are `.gz`, so use a reader the sandbox allows). If the old FEM value was **above** the closed form, as the new one is, the field moved by (1.006712 − 1.00728) / 1.00728 ≈ **−0.056 %**. That is inside the mesher class (0.032–0.233 %), so re-record. If it was **below**, the move is ≈ 1.40 %, above the contract's ≈ 0.5 % stop: do not re-record, and file it. Also print the cell-count drift if the old log carries `cells=`. Note that the error % itself moves 7.8 % relative; that is not the drift statistic.

**Negative controls (rule (e)).**
- **(n1) Asserted, by the executor on the readings:** the drift gate (D). This is the control that keeps a physics move from being re-recorded.
- **(n2) Not re-bought:** the band's detection power is already recorded in the module docstring (1R padding → 7.43 %, `:31–33`). A blind solver (B_z ≡ 0) reads 100 %, a ceiling of 20× the 5 % band. Nothing new is run.

**Window.** Real build (magnetostatics; **no** `dolfinx-complex-mode`, **no** `FEM_EM_REQUIRE_COMPLEX`). Standard tier, `-n 2` (the row's and `OPS-22`'s width). Predicted ≈ 60–80 s: `tests/environment` ≈ 45 s plus the node at 10.91 s (`…WF-15-s2-W4-wf13.log:263`). Foreground, bare, whole-command harness call; Bash `timeout: 660000`.
`scripts/testing/run_and_log.sh MAG-14-rerecord-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_helmholtz_magnitude.py -v -s --tb=short'"`
- Collected-count identity: W1 collects the module's node count plus `tests/environment`. Print it. No census found a red here, so this is a record of the count, not a comparison.
- Cost: class lane 1, ≈ 7 measured min. **≈ 20 predicted**: 15 fixed + W1 ≈ 1.5 + two-line edit + commit.

**Traps already paid for.**
- **`-s` or nothing prints.** pytest captures prints without it.
- **No pipe or redirect** inside or after the harness call (the footer records the pipe's status; `bash_guard` rule 3).
- **Real build only.** A complex-mode source is not needed here. Since `OPS-22` the module also runs in complex, but that is not this record's build.
- **The installed and in-repo paths are different imports** (`WF-13`). Never take the installed log's digits as the in-repo record without W1.
- **`.gz` logs:** a shell decompressor may be denied in the sandbox. If the old FEM B_z cannot be read, (D) is unresolved: stop as INCONCLUSIVE and do not re-record.
- **Write ≈, never `~`**; keep the old value in-text with provenance (contract rule 1). The row form: "0.671 % vs closed form on the 0.11 image (`<W1 log>`); 0.728 % on the pre-0.11 image (`20260819T093105Z`, 2026-08-19); 1.731 % before `GEO-8`; 11 s, in CI".

**Scope.** Covered: the centre-field figure in the §7 row and its guide copy. Not claimed: any band move; the mean / CV figures as new gates (W1 prints them); the docstring's history; complex-mode digits; any other magnetostatic record (`MAG-1` / `MAG-4` loop relL2 is a separate (1*) question if it drifted).

**Negative result.** Report and stop. Never edit the module's bands.
- **(A) red:** the 0.11 image breaks the gate itself. File a known-issues entry and mark item 281 BLOCKED.
- **(Rep) mismatch between in-repo W1 and the installed windows:** the two import paths solve differently. Record both readings in known-issues and **do not re-record**. The review owns the finding.
- **(D) above ≈ 0.5 %, or unreadable:** a physics-sized move. File it in known-issues and do not re-record (contract rule 3).
- A killed window is a finding, not a retry. Check `pgrep -c python3` first.

No `PREDICTED RED:` line.

**Status it can move.** No glyph: `MAG-14` stays ✅. On green: the §7 row's figure becomes the re-recorded 0.671 % (0.11 image) with 0.728 % kept in-text and its provenance; the guide copy moves in the same commit; the §10 ruling clause "a (1*) re-record candidate" is marked done by the review. Commit W1's log, the two edits and the test-results row together.

**Outcome (2026-10-07 06:00 slot, third item, `record-reconciler`): ✅ re-recorded, first try.** W1 `20261007T112319Z_MAG-14-rerecord-W1.log` (77 s, exit 0; 13 collected, 9 passed, 4 complex-only skips). (A) `test_helmholtz_centre_field_magnitude` green. (Rep) `:238–241` prints `cells=52844`, `FEM=3.554758e-09 T`, closed form `3.531057e-09 T`, `rel err=0.671%`, mean 0.572 %, CV 0.2056 % — character-equal to the installed `WF-13` windows. (D) old FEM B_z 3.556767e-09 T (above the closed form, as now) from `20260801T003528Z_GEO-8-after.log.gz:1919–1922`, old cells 53 365: B_z drift −0.0565 %, inside the mesher class 0.032–0.233 %; cell drift −0.98 % (printed). The `20260819T093105Z` log is not on disk; its 0.728 % / 0.644 % / 0.1602 % in the `OPS-22` annotation agree with the `GEO-8`-after log. Sites (s1) §7 row and (s2) `04_helmholtz_analytic_comparison.md:81` moved to 0.671 % with 0.728 % and its provenance kept in-text; the test module is untouched, no band moved.
