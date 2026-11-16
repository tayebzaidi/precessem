# GEO-1 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 186 — 2026-10-02 (lane 1)

**Outcome (2026-10-03 22:30 CDT slot): DONE — `GEO-1` 🧪 → ✅; known-issues 2026-08-27 retired.**

- **Edit:** `tests/mesh/test_birdcage_port_tags.py` only, +92 lines, nothing removed: `GRADED_CONTROL_RESOLUTION = 4.8e-3` (0.11-tagged), the partition test on it with function-level (A0) `BASELINE_CONTROL_RESOLUTION ==` check, and `test_ungraded_build_still_hits_the_documented_floor` last. No `src/` change.
- **Windows (real):** W1 `-n 2` `20261003T035654Z_GEO-1-W1.log` (102 s vs predicted ≈ 70 s; 14 passed / 4 skipped), W2 `-n 1` `20261003T035842Z_GEO-1-W2.log` (62 s; 12 passed / 4 skipped). No edit after W1.
- **Readings (identical at both widths):** 33 185 cells (= prediction); V_mesh / V_box and Σ tags / V_mesh 1.000000000000; ports 8.000000e-07 each; conductor 0.8104 (= 08-26's 0.810427); phantom 0.9705 (predicted ≈ 0.97).
- **Control:** rank 0 `Invalid boundary mesh (overlapping facets) on surface 59 surface 79` (the old docstring quoted 3 / 49 — a note, the substring check holds); rank 1 `birdcage_port_domain geometry generation failed on rank 0`; allreduce caught = size at both widths.

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.*

**`GEO-1`: the partition test onto the graded control.**
- **Executor:** implementer. **Tier:** standard, real build, `-n 2` then `-n 1`.
- **Cost:** ≈ 4 measured / ≈ 15 predicted min.
- **Depends on:** nothing. **Attempt:** 1 of 4.

**Why this is queueable.**
- **The red test.** `test_birdcage_volumes_partition_the_box` (`test_birdcage_port_tags.py:145–258`) is the 🔴 known-issues entry of 2026-08-27. It calls `birdcage_port_domain` ungraded, and that call aborts with "overlapping facets" on 0.11.
- **The floor is ruled.** The 08-30 weekly ruled that floor a documented, uncommissioned limitation (2026-08-25 entry): "every gate that reads this generator is green on a graded control (h_c = 4.8e-3, ruling (b) of 08-26)".
- **This is the one call site never moved.** `GEO-21` step 2 and `GEO-23` step 2b each moved a call site under that ruling, not a band. This test is the exception (`GEO-23` row, "a review's call").
- **The weekly named the fix.** The 09-30 weekly named it as `GEO-1`'s known-issues fix.

**Review ruling (2026-10-02 18:00).** This review extends the 08-26 ruling (b) to this last call site. The move is the same one `GEO-21` and `GEO-23` made: a test-side sizing move, with every identity and band unchanged and the ungraded floor asserted as a control. It is not a generator fix (uncommissioned, 08-30 weekly), and not a loosening: the identities stay at 1e-9. The weekly ratifies the ruling.

**Edit: `tests/mesh/test_birdcage_port_tags.py` only.** The edit is additive and removes no name, since ten modules import from this file.
- **The constant.** Add `GRADED_CONTROL_RESOLUTION = 4.8e-3`, version-tagged (0.11 image), citing the `GEO-21` ladder (`test_birdcage_conductor_sizing.py:65–101`).
- **The partition test.**
  - Pass the constant as `conductor_resolution`.
  - Its docstring cites the 08-26 ruling (b), the 08-30 weekly, both known-issues entries and this ruling.
  - Print the global cell count.
- **(A0)** Inside the test function, import `BASELINE_CONTROL_RESOLUTION` and assert it equals the new constant. Import it at function level: a module-level import would be circular.
- **The negative control.** Add `test_ungraded_build_still_hits_the_documented_floor` as the last test in the module.
  - Call with `conductor_resolution=None` inside `pytest.raises(Exception)`.
  - Rank 0's message contains "overlapping facets". Every other rank's contains "birdcage_port_domain geometry generation failed on rank 0" (`mesh.py:4028–4040`).
  - Assert `comm.allreduce(caught) == comm.size`.
  - **Asserted, with backing:** `20260828T140313Z_GEO-23-step1b-birdcagepart-n1.log` and the `-n 2` log `20260828T140326Z…` measured this exact call.

**Anchors (identities unchanged, at 1e-9).**
- **Identities:** V_mesh / V_box; Σ tags / V_mesh; the four port boxes = 8.0e-7; the tag set exactly.
- **Region bands, unchanged:**
  - **Conductor**, in (0.65, 1)·the analytic sum: **asserted**, since this exact call measured 0.810427 (`20260826T093202Z_GEO-21-step2-gate.log:11075`).
  - **Phantom**, in (0.90, 1): **predicted** ≈ 0.97 (0.9734 on 0.7.2, `20260803T200358Z_GEO-9-step2b-bands.log`).
- **Cells:** predicted 33 185 (`:11075`). Printed, not asserted.

**Windows.** Bare, foreground, `timeout: 660000`.
- **W1** (≈ 70 s):

  `scripts/testing/run_and_log.sh GEO-1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_conductor_sizing.py tests/mesh/test_birdcage_port_tags.py -v -s --tb=short'"`
- **W2** (≈ 25 s): the same with `-n 1` and without the conductor-sizing module, labelled `GEO-1-W2`.

**Traps:**
- **Contamination.** A failed gmsh build contaminates later builds in the same process (`GEO-23`, the `IndexError … size 0` signature). Keep the raise-control **last** and the importer module **first**, as ordered above.
- **Never call the ungraded build outside `pytest.raises`.**
- **Do not touch** `test_birdcage_conductor_sizing.py`.
- **Run hygiene:** use `-s`, and make no edit after W2.
- **File sharing with item 184.** Item 184 imports from this module but does not edit it; the two items are independent.

**Scope:**
- No `src/` change, no band move, no re-record.
- The 2026-08-25 entry (the generator floor) stays OPEN. Its Retire-when gains: "the item-186 control goes red".

**Status it can move:**
- `GEO-1` 🧪 → ✅ on the row this review amended (rule (j)).
- Known-issues 2026-08-27 **retires**: its test-red is gone, and the limitation lives on in 2026-08-25.

**Negative result:**
- **An identity is red at 4.8e-3:** a regression since 08-26. Annotate the 08-27 entry and stop.
- **The phantom band misses:** report and stop; never widen.
- **The control does not raise:** the floor moved. Record the cell count and stop. Retiring 08-25 is the review's decision, and only after it re-measures the ladder.
- **Any stop:** park and mark BLOCKED (rule (d)).
- There is no PREDICTED RED line.

**Audit, 2026-10-03 03:00 daily review (§4): PASS, and the known-issues 2026-08-27 retirement is justified.**
- **Anchors** are asserted at `test_birdcage_port_tags.py:266, :271, :283, :287, :292`. The control is `n_caught == comm.size`.
- **Readings:** 33 185 cells (`20261003T035654Z_GEO-1-W1.log:5846`), identities and regions at W1:5849–5857, and the control raising on every rank at W1:6156–6160. The window ran 102 s with Status 0 (W1:6293–6294). At `-n 1` the window ran 62 s (`…035842Z_GEO-1-W2.log:1936–1937`).
- **Red baseline:** `20260828T140313Z_GEO-23-step1b-birdcagepart-n1.log.gz:58, :358`. It shows the same "overlapping facets on surface 59 surface 79" string.
- **No assertion was edited:** the diff is +92 lines, −0. No `src/` file changed.
- **Caveats:**
  - The ruling (b) extension is still unratified; the weekly ratifies it.
  - The `-n 2` red log `20260828T140326Z…-n2-ports` is missing from `docs/testing/logs/`, a retention gap. It gates nothing, because W1 and W2 re-measure the control.
- **Review spot-check:** W1:5846 was re-read, and it matches.
