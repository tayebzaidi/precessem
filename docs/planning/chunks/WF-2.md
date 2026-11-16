# WF-2 — reproducible output bundle manifest

### §9 item 346 — 2026-10-08 (lane 1, rule (j) re-pose at L0_BAND, provisional)

**Outcome (2026-10-09, 00:00 CDT slot, take-next): ✅ complete, provisional for the 10-10 weekly.** New module `tests/package/test_wf2_bundle_reproducible.py`; no `src/` edit; `control_sigma_ulp.yaml` already existed. W1 sphere `20261009T053015Z_WF-2-W1-sphere.log` 6 passed (48 s): (D) no deterministic diff, `sha256` = recomputed; (F) file = writer's memory bitwise (each rank's own summary differs only in `elapsed_s`, printed); (R) SAR / power 1.453e-12, `volume_m3` bitwise equal as predicted; (X) 6.222e-11 / 6.076e-11 vs `MAT4_DEG1_RECORD`; (c1) flags `{degree}`, miss 3.3757 % (against `run_a`, not the brief's record-to-record 3.49 %); (c3) exact; (c4) flags `{sha256, source}`, SAR rel 7.566e-12 (≤ 1e-10 predicted, float layer blind); (c5) 9.999e-09. W2 four-port `20261009T053110Z_WF-2-W2-four-port.log` 5 passed (77 s): (D) / (F) clean; (R) 5.010e-16; (X) 9.557e-11 both; (T) 5 header lines bitwise bar `generated_utc`, no `PLACEHOLDER`, parsed rel 5.521e-14; Z 5.628e-16 at cond(I − S_a) 19.28 (printed); (c1) 0.0495 (against `run_a`), (c2) 0.2492; (c5) 1.000e-08. **Disclosed:** the first W1 (`20261009T052650Z_WF-2-W1-sphere.log`, Status 124) was a test defect — (F) compared each rank's own summary (rank-local `elapsed_s`) to rank 0's file, rank 1's assert failed (1 failed / 5 passed in 31.16 s) and rank 0 hung in the next collective; the slot fixed the comparator (broadcast the writer's summary) and re-ran the verbatim command at the same ceiling rather than the W1b split, which is for a cost overrun. Orphan check 0. Rule (i): module 05:30:09Z < W1 05:30:15Z, untouched before W2 05:31:10Z.

**Executor:** `implementer`. **Lane 1** (§9 rule (5), revalidation of a ⚠️ chunk). **Provisional rule-(j) re-pose by the 2026-10-08 18:00 review, under the 10-07 weekly's licence** (PROJECT_PLAN.md, the 10-07 weekly's lane-1 targets: "`WF-2` (as 09-30; its 1e-12 target sits below the `-n 2` drift — the daily may re-pose it at `L0_BAND` under rule (j), the weekly ratifies)"). **The 10-10 weekly ratifies or reverts.**

**Why this item exists, and what changes.**
- The 09-30 weekly's target (`docs/planning/plan-archive.md:36207`): "a re-run at fixed width reproduces the manifest's S and SAR to 1e-12 and a planted change is detected."
- It was declined three times (2026-09-29 18:00 §6; 2026-10-02 18:00, journal lines 85–88; 2026-10-06 18:00). There were two reasons:
  - (a) The only manifest `WF-2` ever wrote (`examples/mri/01_coil_phantom_fields.py:209–267`, commit `74ee3d6e`) records neither S nor SAR. It holds existence flags, parameters and E / B means, and no test covers it.
  - (b) 1e-12 sits below the measured `-n 2` run-to-run drift. On this very fixture, `MAT-4` step 1's lossy sphere (74 020 cells), mean SAR drifts by up to **4.45e-11** relative between repeated solves in one process (known-issues "Second fixture (2026-09-29)" row; `20260929T183957Z_MAT-1-probe-noise.log:105–106`).
- (a) is dissolved by work that has since landed, and no new feature is needed. Since `WF-1` step 1 (§9 item 256, ✅ 2026-10-06), `run_study_file` (`src/fem_em_solver/study_yaml.py:496–557`) writes the package's output bundle. That bundle is `summary.json` (schema `femem`, the `sha256` of the study file, `family`, `degree`, `frequency_hz`, `ranks`, `cells`, the `mean_sar` dict or the `s_matrix` / `z_matrix`, the `touchstone` name and `elapsed_s`), plus `<stem>.s4p` for a port study. This is the only manifest in the repository that carries S and SAR.
- **The re-pose keeps the weekly's question and changes only two things:**
  - **(1) The comparator:** 1e-12 becomes `L0_BAND` = 1e-9 (`tests/validation/test_mat5_em_coupled_sphere.py:104`) on floating fields, and bitwise equality on deterministic fields.
  - **(2) The object:** "the manifest" is read as `femem run`'s bundle, because the `mri:1` manifest records neither quantity the question names.
  - Both are this review's reading. The weekly ratifies both.
- **The ⚠️ no-new-features rule holds.** No `src/` edit is made. The comparator is test-side, in the new test module only.

**Anchor (what is asserted).** Module `tests/package/test_wf2_bundle_reproducible.py`, complex build, `-n 2` (the width of both records: `ANS4_EXACT_STRIP_RECORD_RANK_WIDTH`, as `WF-1` asserts). Each fixture is run twice through `cli.main(["run", <file>, "--output-dir", <dir>])` in one process, into fresh directories `logs/wf2/<stem>/run_a` and `run_b`. The fixture files are `tests/package/wf1/sphere_mat4.yaml` and `tests/package/wf1/f_small_four_port.yaml`, used as-is and not copied or edited. The test-side comparator then reads both `summary.json` files. The labels are:
- **(D) Deterministic fields, bitwise.**
  - `femem`, `source`, `sha256`, `family`, `degree`, `frequency_hz` (compared by `float.hex`), `ranks` and `cells` are equal between `run_a` and `run_b`, and so is `touchstone` for the port study. The two key sets are also identical.
  - `sha256` also equals `hashlib.sha256(<file>.read_bytes()).hexdigest()`, recomputed in the test.
  - `elapsed_s` is the only key excluded, because it is wall clock; it is printed for both runs.
  - Backing: `cells` equal across two meshings in one process, `WF-1` (S) (`20261006T003504Z_WF-1-W1-parse-sphere.log`, `…003650Z_WF-1-W2-f-small.log`). The rest are deterministic by construction.
- **(F) File equals memory, bitwise.** Each `summary.json` loaded with `json.load` equals the `StudyRun.summary` the same `run_study_file` call returned, compared by `float.hex` on every float. JSON float `repr` round-trips exactly.
  - This needs `run_study_file` (the same function `cli.main` calls) for the in-memory side. Call it directly for (F), or take the `StudyRun` from a direct call in place of the CLI on `run_b`. Either way the comparison is file against memory from one call.
- **(R) Floating fields reproduce at `L0_BAND`.**
  - Sphere: |SAR_b − SAR_a| / SAR_a ≤ `L0_BAND` on `mean_sar_w_per_kg`, and the same comparison on `dissipated_power_w` and `volume_m3`.
  - Four-port: ‖S_b − S_a‖_F / ‖S_a‖_F ≤ `L0_BAND`, using `_rel` from `tests/package/test_wf16_study_four_port.py`.
  - Backing (same fixture, same width):
    - Sphere: in-process drift ≤ 4.45e-11 (cited above).
    - Four-port: in-process (I) 9.110e-16 (`20261005T202623Z_WF-16-W1-deg1.log`).
- **(X) Cross-process re-run: each bundle against the record written by an earlier process.**
  - Sphere: `run_a` and `run_b` SAR against `MAT4_DEG1_RECORD`, imported from `tests/package/test_wf15_study_lossy_sphere.py`, ≤ `L0_BAND`. Backing: `WF-1` W1 (R1) 6.222e-11.
  - Four-port: `run_a` and `run_b` S against `_record(1)` (`P1["64"]`), ≤ `L0_BAND`. Backing: 9.557e-11, `WF-16` W1 and `WF-1` W2 (R1).
  - **(X) is the "re-run at fixed width" in its cross-process sense.** (R) is the same reading in-process.
- **(T) Touchstone, port study only.** `run_a/f_small_four_port.s4p` and `run_b/f_small_four_port.s4p`:
  - Every `!` / `#` header line except the `! generated_utc:` line (`src/fem_em_solver/ports/touchstone.py:177`, wall clock by design) is bitwise equal.
  - Neither file contains `PLACEHOLDER`.
  - The parsed matrices (`load_touchstone`) agree at `L0_BAND` (`_rel`).
- **Printed, not asserted (records):**
  - `z_matrix` rel(b, a), with cond(I − S_a) printed beside it. Z is not in the question, and its run-to-run drift through (I − S)⁻¹ is unmeasured.
  - Both `elapsed_s`.
  - Whether `volume_m3` is bitwise equal between runs (prediction: yes, because it comes from the mesh only).
  - Every ratio of a reading to `L0_BAND`.

**Negative controls (rule (e) labels).**
- **(c1) Degree plant, asserted.** This needs no solve. Take a copy of `run_b`'s summary dict, set `degree` to 2, and replace the floating payload with the degree-2 record on the same fixture: `WF15_DEG2_RECORD` 8.0116577245e-08 (import it from `tests/package/test_wf1_study_yaml.py`) for the sphere, and `_record(2)` (`P2_X1`) for the four-port. The comparator must then:
  - flag exactly `{degree}` in the bitwise layer;
  - read a floating miss above `L0_BAND`.
  - Backed: sphere 3.49 % (|`MAT4_DEG1_RECORD` / `WF15_DEG2_RECORD` − 1|, `WF-1` W1 (c2), asserted > `L0_BAND` there). Four-port **0.0491** (`WF-16` W1 (c1), `WF-1` W2 (c3)). Same comparison, same fixture.
- **(c2) Port-swap plant, asserted.** `run_b`'s S with P1↔P2 permuted (`perm = [1, 0, 2, 3]`) misses `run_a` by more than `L0_BAND`. Backed: **0.2492** (`WF-16` W1 (c2), `WF-1` W2 (c4)), same comparator and fixture.
- **(c3) Key plants, asserted.** These are structural and true by construction. In a copy of `run_b`'s summary:
  - `cells` + 1 is flagged as exactly `{cells}`;
  - one hex digit of `sha256` changed is flagged as exactly `{sha256}`;
  - a deleted `ranks` key is flagged as a key-set mismatch.
- **(c4) One-ULP config, sphere only, one extra solve.**
  - Run `tests/package/wf1/control_sigma_ulp.yaml` through the CLI into `logs/wf2/control_sigma_ulp/run_c`. That file is `sphere_mat4.yaml` with σ = `nextafter(0.57, 1)` and a different comment line.
  - **Asserted:** its `sha256` differs from `run_a`'s, and the bitwise layer flags exactly `{source, sha256}`.
  - **Predicted, printed only:** the float layer is blind here. SAR rel(c, a) is expected at ≤ 1e-10 (2⁻⁵³ in σ plus drift). `WF-1` (c1) measured the DG0 σ difference at exactly 2⁻⁵³, but no one has measured this SAR comparison.
  - This control is the reason the manifest carries a hash: it catches a change that the floats cannot see.
- **(c5) Resolution plant, predicted, printed only.** `run_b`'s SAR and S × (1 + 1e-8) read a comparator value of ≈ 1e-8 (≈ 10× `L0_BAND`). This is arithmetic on the test's own numbers, but there is no prior same-fixture measurement, so it is not asserted.

**Tier, ranks, measured cost.** Complex, `-n 2`. Two windows; never run them at once.
- **W1 — the sphere and the comparator nodes.** Standard, 180 s ceiling.
  - Contents: three degree-1 sphere CLI runs (a, b, the (c4) control) plus the solve-free (c1) / (c3) plants.
  - Predicted 70–130 s. Basis: `WF-1` W1 measured **100 s** for two degree-1 sphere solves plus parse nodes and DG0 builds (`20261006T003504Z_WF-1-W1-parse-sphere.log`). A sphere solve is ≈ 14.2 s (`WF-15` W1).
- **W2 — the F-small four-port.** Heavy by ceiling, 400 s, the ceiling `WF-1` W2 and `WF-16` W1 had on the same fixture at the same width.
  - Contents: two degree-1 CLI runs.
  - Predicted 110–170 s. Basis: `WF-1` W2 measured **142 s** for two F-small solves at `-n 2` (`20261006T003650Z_WF-1-W2-f-small.log`).
- **Item tier:** heavy by ceiling (W2), `-n 2`. Total ≈ 4 min of compute.

**Harness commands** (each one whole, bare Bash call; Bash `timeout: 660000`; `tests/environment/test_complex_mode.py` first in each):

```
scripts/testing/run_and_log.sh WF-2-W1-sphere "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf2_bundle_reproducible.py::test_wf2_comparator_flags_planted_changes tests/package/test_wf2_bundle_reproducible.py::test_wf2_sphere_bundle_reproduces -s -v --tb=short'"
```

```
scripts/testing/run_and_log.sh WF-2-W2-four-port "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/package/test_wf2_bundle_reproducible.py::test_wf2_four_port_bundle_reproduces -s -v --tb=short'"
```

- **Node split.** `test_wf2_comparator_flags_planted_changes` holds the solve-free (c3) plants on synthetic dicts. The (c1) / (c2) plants need the solved `run_b` bundles, so they live inside the two solve nodes.
- **Fallback if W1 is killed at 180 s.** Move (c4) to its own node and window, `WF-2-W1b-ulp`, with the same command shape and ceiling. Never lengthen the timeout.

**Traps already paid for.**
- The harness call is the whole Bash command: bare, repo-relative, with no `cd … &&` prefix, no `| tail` and no redirect. `bash_guard` rule 3 denies anything else. Pass `timeout: 660000` on every window and run foreground. Never background a harness run and end the turn.
- `-s` on every pytest (rule (g)). `timeout -k 30` inside the container: a plain TERM does not stop `mpiexec` (MAT-6 step 10). After a kill, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` before anything else.
- Use the complex build: `source /usr/local/bin/dolfinx-complex-mode` plus `FEM_EM_REQUIRE_COMPLEX=1`, with `tests/environment/test_complex_mode.py` collected first. In the real build `TimeHarmonicSolver.solve` raises.
- **Stale bundles pass vacuously.** `logs/wf1/<stem>/` already holds `WF-1`'s bundles, and an old `summary.json` in a reused directory would compare equal to itself.
  - Write only under `logs/wf2/…`.
  - Call `_fresh_dir` (import it from `test_wf1_study_yaml.py` or copy its three lines) on every output directory before its run.
  - Assert that `run_a` and `run_b` are distinct paths and that each `summary.json`'s mtime is later than the test start.
- Rank 0 writes and then `Barrier`, so read `summary.json` only after `cli.main` returns. Reduce `rc` with `MPI.MIN` / `MPI.MAX` (`_rc_all`). The whole comparison is done on every rank from the same file, with no rank-local assertions on solve data.
- **No solve output is compared bitwise.** `-n 2` MUMPS is not bitwise run to run (known-issues "Second fixture (2026-09-29)" row; `WF-15`'s and `WF-1`'s docstrings). Bitwise applies only to the deterministic keys and to file against memory.
- The Touchstone `! generated_utc:` line differs by design. Exclude it by its exact prefix, never by "all comment lines".
- Import bands and records; never restate them: `L0_BAND`, `MAT4_DEG1_RECORD`, `_record`, `_rel`, `WF15_DEG2_RECORD`, `ANS4_EXACT_STRIP_RECORD_RANK_WIDTH`. `LEG_D0_REPRODUCTION_BAND` equals `L0_BAND` by value, but the licence names `L0_BAND`, so use `L0_BAND`.
- No edit to the module after its last green window (rule (i)). No edit under `src/` (⚠️ no-new-features rule). If the comparator wants to become a package API, that is Phase 9 work for the weekly, not this item.

**Scope boundary — not claimed:**
- Reproducibility across widths (`-n 1` against `-n 2`), degree 2, or the thermal summary (`WF-1` step 3).
- The `mri:1` example manifest (`_write_output_manifest`). Its fields stay an example artifact covered only by the census. Whether it is retired, or rewired to `run_study_file`, is the weekly's call (question below).
- Z reproduction (printed only), wall-clock fields, or any bundle-comparison API.
- The physics of the records reproduced, which belongs to `MAT-4` / `ANS-4` / `PORT-11`.

**What a negative result means.**
- **A floating miss in (R) or (X)** contradicts two same-fixture measurements (4.45e-11; 9.110e-16 / 9.557e-11).
  - Report the reading and stop.
  - Add a known-issues entry beside the "Second fixture (2026-09-29)" determinism row, with both bundles kept under `logs/wf2/`.
  - `WF-2` stays ⚠️ and no band moves.
- **A deterministic miss in (D)** (a `cells` or `sha256` difference between two runs of one file) is a defect: meshing or hashing non-determinism. Open a known-issues entry and stop; `WF-2` stays ⚠️.
- **A Z-only drift above `L0_BAND` with S green** is a printed record, not a red. Name it in the journal with cond(I − S).
- **A control that does not fire** ((c1)–(c3)) is a defect in the test-side comparator. Fix it in the same slot (it is test code) and re-run the window. The (c4) and (c5) readings are records either way.
- Park any incomplete work on `attempt/WF-2-…`, journal it in `docs/testing/attempts.md`, and mark the §9 item BLOCKED in the same commit (rule (d)).

**Status this result can move.**
- **Green ⇒ `WF-2` ⚠️ → ✅** on the amended Done-when below.
- `WF-2` leaves §2's ⚠️ list (PROJECT_PLAN.md ≈ line 780, leaving `TH-2`) and §9 rule (5)'s lane-1 list (≈ line 3811).
- The note under the WF table ("`WF-2`/`WF-3` produce structurally valid manifests … physically meaningless numbers", ≈ line 2935) is struck for `WF-2`. `WF-3` is already ✅.
- The weekly ratifies the reading. If it reverts the reading, the result stands as a record and the glyph returns to ⚠️.

**Proposed §7 Done-when** (the row has none today; the review writes it under rule (j), provisional, for the weekly to ratify):

> **Done-when (§4) (written 2026-10-08 18:00 review, rule (j), provisional — the 10-07 weekly's `L0_BAND` licence; the weekly ratifies; §9 item 346):** two `femem run`s of the same study file at `-n 2`, complex, degree 1 — `WF-1`'s lossy sphere and F-small four-port files — write output bundles (`summary.json` + `.s4p`, `run_study_file`'s manifest, read as "the manifest" because the `mri:1` manifest records neither S nor SAR) whose deterministic fields (`femem`, `source`, `sha256` = the file's recomputed hash, `family`, `degree`, `frequency_hz`, `ranks`, `cells`, `touchstone`, Touchstone header bar `generated_utc`) are bitwise equal and equal file-to-memory, and whose S and SAR agree with each other and with `MAT4_DEG1_RECORD` / `P1["64"]` within `L0_BAND` (1e-9, the measured `-n 2` drift being ≤ 4.45e-11 — the old 1e-12 sat below it); a planted degree change, a P1↔P2 swap and a one-ULP σ config are each detected (the last by `sha256` alone); harness, elapsed recorded ⇒ ✅. Not claimed: cross-width reproduction, degree 2, the `mri:1` example manifest, Z.

**Review rulings (2026-10-08 18:00 review) on item 346 — these outrank the brief above.** Written as a **provisional** rule-(j) re-pose under the 10-07 weekly's `L0_BAND` licence; the weekly ratifies or reverts, including the reading of "the manifest" as `run_study_file`'s bundle (Q1) and the fate of `mri:1`'s `_write_output_manifest` (Q2). (3) Two calls in one process, with the cross-process half carried by the registered records, are accepted. (4) Z printed beside cond(I − S), not asserted. (5) The §7 tier becomes heavy (by ceiling) in this review's row edit. The §7 row gains the drafted Done-when, tagged provisional.
