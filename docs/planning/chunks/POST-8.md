# POST-8 — Pointwise samples of a degree-1 N1curl field depend on the MPI partition: `ans:2`'s Row 1 moves −19 % … +112 % in SAR between `-n 2` and `-n 4`, while every …

 **Pointwise samples of a degree-1 N1curl field depend on the MPI partition: `ans:2`'s Row 1 moves −19 % … +112 % in SAR between `-n 2` and `-n 4`, while every integral agrees to ≈ 12 figures** (opened 2026-09-27 by the 03:00 daily review from the 09-26 09:00 slot's disclosure and a `log-pathologist` ruling; known-issues 2026-09-27) | 🟡 **Step 1 DONE 2026-09-27 (11:20 take-next; heavy-budgeted, 213 s / 174 s):** `scripts/probes/post8_pointwise_partition.py` on `0b00091`. Anchors green: P1 row reproduces its own width's logged row to 8.94e-11 (`-n 2`) / 4.95e-10 (`-n 4`); P1 phantom power within 4.80e-12 / 6.62e-14 of `metrics.json:100`. Mechanism reading: each named point has **one** incident cell; the raw N1curl `E` in it is width-invariant (16 point×drive SARs agree `-n 2` vs `-n 4` to ≤ 1.0e-9, print resolution), while the CG1 sample differs from it at both widths (+x, P1: N1curl 5.514140460e-07; CG1 4.858355042e-07 at `-n 2`, 5.878191741e-07 at `-n 4`, W/kg) and the cell's CG1 nodal |E|² sits at 0.83–0.98 of its own N1curl vertex value ⇒ the N1curl → CG1 nodal interpolation carries the dependence. The pre-registered "bracket" readout does not apply (one cell, no range). Census printed by text heuristic (18 candidate files; magnetostatic hits likely false positives — weekly to classify). `20260927T112115Z_POST-8-step1-n2.log`, `20260927T112501Z_POST-8-step1-n4.log`. Sampling rule: the weekly's. Original step text — **Step 1 (§9 item 54, standard):** a same-code control and a mechanism reading. Two probe windows of `ans:2`'s case at `-n 2` and `-n 4` on one commit, with the tracked-artifact writes switched off. **Asserted:** each width reproduces its own logged Row 1 to 1e-9 relative, and the whole-phantom power is width-invariant to 1e-9. **Printed:** at each named point, the raw N1curl `E` in every incident cell next to both CG1 values. The Nédélec kernel is exact for this fixture's straight-sided tets, so the per-cell values are the partition-free ground truth. Also a census of which examples and gates *sample points* from an N1curl → CG1 interpolant. **Done-when (§4):** a sampling rule whose output is width-invariant to 1e-9 relative, asserted at `-n 2` and `-n 4` on `ans:2`'s points, and adopted by every point-sampling site the census found. The rule is chosen by the weekly from step 1's reading (e.g. the owning cell by a deterministic global rule, or a DG0 cell average). Then `ans:2`'s Row 1 is re-issued and the 09-26 `ANS-2` pointwise relative level re-read. **No band moves; no integral quantity is in question.** 

### §9 item 177 — 2026-10-02 (lane 1; the 09-30 weekly's sampling rule)

**Outcome (2026-10-02 15:00 slot): DONE — the rule implemented opt-in and gated at `-n 1 / 2 / 4`.** All four windows Status 0, one tree.

- **Code:** `post.evaluation.point_owner_cells` (owned candidates, int64 `Allreduce(MIN)` on `original_cell_index`, `Allreduce(SUM)` claim guard) and the keyword-only `owner=None | "lowest_original_cell"` on `evaluate_vector_field_parallel`; `owner` rides the OPS-40/55 allgather key; the legacy path is untouched (guard module 6 passed / 1 skipped at every width; C1a reproduces both logged CG1 rows).
- **Windows:** W1 `-n 2` `20261002T200331Z_POST-8-step2-W1.log` (195 s; 23 passed / 1 skipped); W2 `-n 4` `20261002T200654Z_POST-8-step2-W2.log` (191 s); W3 `-n 1` `20261002T201013Z_POST-8-step2-W3.log` (228 s); W4 real build `20261002T201406Z_POST-8-step2-W4-real.log` (49 s; 17 passed / 5 skipped). `git status --porcelain examples/` empty after W4.
- **K0:** `original_cell_index` is an int64 ndarray covering owned + ghost cells in 0.11 (shape (223,) at `-n 2`, 191 owned + 32 ghosts); owned part unique, 384 = `size_global`. **K1 / K2:** CG1 ≤ 8.9e-16, N1curl ≤ 2.2e-15 (≤ 1e-12). **K3:** oracle incident sets 1 / 2 / 24 (vertex prediction 24 met); the rule returns the minimum 217 / 223 / 126 at all three widths, imaginary part 0. **K4, K5, K6** green (K5 raises on every rank at `-n 2 / 4`, does not at `-n 1`). **CK2:** the legacy DG0 value at the vertex is in the incident set and is 126 (= min) at `-n 2`, 128 at `-n 4`, 131 at `-n 1`; at the facet 232 ≠ min at `-n 4` — the legacy path's partition dependence shown on the known answer.
- **A0:** 199 920 cells at every width; P1 phantom power rel 4.91e-12 / 8.93e-13 / 1.07e-11 (`-n 2 / 4 / 1`). **A1:** one owned containing cell per point at every width. **A2 (the gate):** max rel vs `STEP1_N1CURL_SAR` **1.37e-10 / 6.55e-11 / 3.03e-10** (≤ 1e-9). **A3:** owner-rank `eval` bitwise and all ranks bitwise at every width.
- **Cross-width record (not asserted):** pairwise max relative difference of the 16 SARs: `-n 2` vs `-n 4` 1.33e-10, `-n 2` vs `-n 1` 1.66e-10, `-n 4` vs `-n 1` **2.98e-10**, all at P3 / −y (the smallest-magnitude entry class, 7.48e-9 W/kg). **The printed prediction ≤ 1e-11 is missed** (a negative result on a prediction, rule (e)); the asserted bound ≤ 2e-9 holds. The owner cell is the same (original index 150135 / 162464 / 142550 / 169932 at every width), so this is the solve's own width-dependent round-off at ≈ 1e-10 (P1 phantom power also moves 1e-11 across widths), not sampling.
- **C1a:** 8.89e-11 (`-n 2`) / 5.89e-10 (`-n 4`) ≤ 1e-9. **C1b:** min |rel| vs the other width 0.1293 / 0.1485 ≥ 0.10. **C1c (`-n 1`, printed):** P1 CG1 row −3.43 / +87.8 / −4.30 / +29.8 % from the `-n 2` row and −20.2 / +8.82 / −16.7 / +12.5 % from the `-n 4` row (prediction ≥ 1 % from each: met). **C1d:** P1 +x CG1 / raw 0.881 (`-n 2`), 1.066 (`-n 4`; step 1's reading), 0.851 (`-n 1`). **C2:** owner `local_to_global` 85735 / 68816 / 194905 at +x for `-n 2 / 4 / 1` while `original_cell_index` stays 150135 — **even `-n 1` renumbers**, so `local_to_global` is not the input order at any width.

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 10:30 daily review.*

`POST-8` step 2 puts the 09-30 weekly's sampling rule into the point-evaluation path, as an **opt-in** keyword on `post.evaluation.evaluate_vector_field_parallel`. The rule is gated on `ans:2`'s named points at `-n 1`, `-n 2` and `-n 4`. It also gets a known answer on a shared vertex and facet, and the CG1 route serves as the negative control.
- **Executor:** implementer. **Tier:** heavy, since each `ans:2` window is ≈ 2–5 min, above the 180 s standard tier.
- **Builds:** complex, plus one cheap real-build window.
- **Cost:** ≈ 4 measured (lane-1 median) / ≈ 50 predicted min.
- **Depends on:** nothing; step 1 landed on `0b00091` (item 54).

**The rule as ruled (§10, 09-30 weekly):** "Raw point evaluation takes a deterministic owner (lowest global cell index among the cells containing the point), gated width-invariant to 1e-9 at `-n 1 / 2 / 4`; public comparison rows quote cell or mass averages."

**Which global cell index: `mesh.topology.original_cell_index`, never `index_map(tdim).local_to_global`.**
- The original index is the cell's position in the input array given to `create_mesh`. Every birdcage mesh goes through `_model_to_mesh` (`src/fem_em_solver/io/mesh.py:12`), so that position is set by the gmsh mesh and not by the partition.
- The index-map global number depends on the partition. Step 1 measured it: the same +x cell, carrying the same `E` digits, is global cell 85735 at `-n 2` and 68816 at `-n 4` (`20260927T112115Z_POST-8-step1-n2.log:1918`, `20260927T112501Z_POST-8-step1-n4.log:1920`). The probe used `local_to_global` (`scripts/probes/post8_pointwise_partition.py:102`); do not copy that line.

**Code (a `src/` change; rule (c) disclosure; additive only).** In `src/fem_em_solver/post/evaluation.py`:
- **New public `point_owner_cells(mesh, points, comm)`.** It returns `owner_original` (int64[N], identical on every rank, −1 for a point outside the mesh) and `owner_local_cell` (int32[N]; the local cell on the owning rank, −1 elsewhere).
  1. `bb_tree` → `compute_collisions_points` → `compute_colliding_cells` (the existing calls, `:99–101`).
  2. Keep only **owned** candidates (`c < size_local`).
  3. `local_min` = the minimum `original_cell_index` over those, with sentinel `np.iinfo(np.int64).max`.
  4. Reduce with the buffer-form `comm.Allreduce(local_min, global_min, op=MPI.MIN)`.
  5. The rank whose `local_min == global_min != sentinel` owns the point.
  6. Defensive: an `Allreduce(SUM)` of claim counts above 1 raises on **every** rank.
- **New keyword-only `owner=None` on `evaluate_vector_field_parallel`.**
  - The accepted values are `None` (the legacy path, byte for byte) and `"lowest_original_cell"`. Anything else raises `ValueError` before any collective.
  - Under the rule, each owning rank calls `function.eval` only on its owned points, in the owner cell. The existing gather / scatter / `bcast` (`:120–135`) is unchanged.
- **`owner` joins the OPS-40/55 key.** The `(n_points, digest)` allgather becomes `(n_points, digest, owner)`, so a rank-mismatched `owner` raises on every rank instead of hanging.
- **The default stays the legacy path.** Every existing caller therefore stays digit-identical by construction: ≈ 30 test modules, `core/solvers.py:763`, `post/sar.py:354`, `ports/gap_voltage.py:289`, and the examples. Flipping the default is adoption, which is step 3.
- Export `point_owner_cells` from `src/fem_em_solver/post/__init__.py` (additive).

**Module 1: `tests/post/test_point_owner_rule.py`** — a known answer, cheap, run in both builds.
- **Fixture:** `create_unit_cube(comm, 4, 4, 4)`, with three points: interior (0.13, 0.27, 0.61), the hex-face plane (0.5, 0.3, 0.6), and the interior vertex (0.5, 0.5, 0.5).
- **Oracle,** independent of DolfinX collision: allgather each owned cell's original index and its four vertices. Rank 0 takes the incident set at barycentrics ≥ −1e-12 and broadcasts it.
- **Asserted:**
  - **(K0)** `original_cell_index` (the owned part, allgathered) is unique, and its count equals `size_global`. Print its type and shape: the 0.11 API may return a per-cell-type list, in which case take `[0]`.
  - **(K1)** CG1 `(x, 2y, 3z)` under the rule meets its closed form to ≤ 1e-12 at all three points (the guard module's bound, `tests/post/test_evaluation_collective_guard.py:83`).
  - **(K2)** N1curl `a + b × x`, with a = (1, −2, 0.5) and b = (0.3, 0.1, −0.2), meets its closed form to ≤ 1e-12.
  - **(K3)** The tie-break. A DG0 scalar holds `original_cell_index` in each cell (owned and ghost, then `scatter_forward`). Under the rule, each point returns exactly the minimum over the oracle's incident set, with imaginary part 0 in complex. The oracle count is 1 at the interior point and ≥ 2 at the facet and vertex points.
  - **(K4)** Every rank holds bitwise-equal values.
  - **(K5)** If rank 0 passes `owner=None` while the others pass the rule, every rank raises. At `-n 1` the same construction does not raise (the control).
  - **(K6)** An unknown `owner` string raises `ValueError`.
- **Controls:**
  - **(CK1, asserted** by construction once K3's count is ≥ 2 with distinct indices**):** the incident set's maximum ≠ its minimum at the vertex and facet points.
  - **(CK2):** the legacy value at the vertex is a member of the incident set (asserted). That it differs from the minimum at some width is predicted and printed only.
- **Printed:** the incident-set sizes (predicted 24 at the vertex, not asserted) and the rule's three answers. The slot compares those answers across W1–W3.

**Module 2: `tests/validation/test_post8_owner_rule_ans2.py`** — `complex_only`, one module-scoped fixture, width taken from `comm.size`.
- **Fixture:**
  - `sweep = build_four_port_sweep(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG)`, imported directly from `tests.validation.test_port_birdcage_four_port` and `tests.validation.test_birdcage_sar_1g_rung`.
  - `_solve_driven(sweep, pid)` for P1…P4, with the pids as `_build_mass_averaged` derives them (`tests/validation/test_birdcage_sar_mass_averaged.py:232, 239–249`), but without its 21 ball integrals.
  - Load `ans:2`'s script by importlib, as the probe does (`post8_pointwise_partition.py:77–82`). Assert `not ans2.RESOLUTION_OVERRIDDEN`.
  - Use only `ans2._pointwise_sar` and `ans2._phantom_power`; **never `main()`**.
  - Import `LOGGED_P1_ROW`, `LOGGED_P1_SOURCE`, `TRACKED_P1_POWER_W`, `ROW_RTOL` and `POWER_RTOL` from the probe; never restate them.
  - Named points: (±`CENTRE_RADIUS_M`, 0, 0) and (0, ±`CENTRE_RADIUS_M`, 0), in `ans:2`'s order (`02_birdcage_coil_driven_sar_10MHz.py:747–755`).
- **Anchors (asserted):**
  - **(A0)** Cells == `ONE_GRAM_RUNG_CELL_RECORD` 199 920 exactly (n2 log :1901, n4 log :1903). P1 phantom power is within `POWER_RTOL` 1e-9 of `TRACKED_P1_POWER_W`; measured 4.80e-12 / 6.62e-14 (:1913 / :1915).
  - **(A1)** Each named point has exactly one owned containing cell globally (Allreduce SUM == 1; n2 log :1916 / 1930 / 1944 / 1958).
  - **(A2) The gate.** The 16 raw N1curl SARs σ|E|²/(2ρ) (with `SALINE_SIGMA` and `PHANTOM_RHO_KG_PER_M3`, as `_pointwise_sar` uses them), taken from `evaluate_vector_field_parallel(e_complex, points, comm, owner="lowest_original_cell")`, each lie within **1e-9 relative** of `STEP1_N1CURL_SAR`. That record is restated with its log citation:
    - Each entry is the mean of the `-n 2` print (n2 log :1917–1965) and the `-n 4` print (n4 log :1919–1967). The two prints agree in 8 of 16 entries and differ by one unit in the 10th figure in the other 8.
    - Comparand budget ≤ 5e-10. Predicted measured ≤ 5e-10 at every width.
    - The `-n 2` values (+x / +y / −x / −y):
      - P1: 5.514140460e-07, 8.056036637e-09, 6.147781056e-08, 1.000625322e-08;
      - P2: 5.413193615e-09, 5.365427637e-07, 7.538407567e-09, 6.062393526e-08;
      - P3: 6.887651145e-08, 7.937373888e-09, 5.038331027e-07, 7.475028492e-09;
      - P4: 9.595073748e-09, 6.651644441e-08, 9.294666807e-09, 4.980486020e-07.
    - Print with `.17e`: the 16 SARs and the E components, each point's owner `original_cell_index` and its `local_to_global`.
    - After W3 the slot records the pairwise maximum relative difference across the three logs, predicted ≤ 1e-11, as a record that is not asserted. Each window's assertion bounds that pairwise difference at ≤ 2e-9 through the common comparand; the docstring says so.
  - **(A3)** On the owning rank the rule's value equals `e_complex.eval(x, [owner_local_cell])` bitwise, and every rank holds the same 16 values bitwise.
- **Negative control, the CG1 route** (`ans2._pointwise_sar` on the same solves, which interpolates to CG1 and then uses the default `owner=None`):
  - **(C1a, asserted, at `-n 2` / `-n 4` only):** the P1 CG1 row reproduces `LOGGED_P1_ROW[width]` to `ROW_RTOL` 1e-9; measured 8.94e-11 / 4.95e-10 (n2 :1907, n4 :1909). This is also the check that the default path did not move.
  - **(C1b, asserted, at `-n 2` / `-n 4`):** the same row against the *other* width's logged row has minimum |rel| ≥ 0.10. Measured −17.35 / −42.04 / −12.93 / −13.36 % (n2 :1908) and +20.99 / +72.54 / +14.85 / +15.42 % (n4 :1910). The minimum is 0.1293, i.e. 1.29e8× the 1e-9 the rule is held to.
  - **(C1c, predicted, printed, `-n 1`):** no logged `-n 1` CG1 row exists. Print it beside both logged rows; predicted ≥ 1 % from each.
  - **(C1d, printed):** the CG1 / raw ratio per (point, drive). Step 1 read 0.881 / 1.066 at P1 +x.
  - **(C2, predicted, printed):** the owner's `local_to_global` differs across widths, while its `original_cell_index` does not.

**Windows.** Foreground, Bash `timeout: 660000`, each harness call bare and whole, pytest `-s`, no pytest `-k` selection, no pipes. Run W1 → W4 in that order.
- **Basis:** step 1's 213 s (`-n 2`) and 174 s (`-n 4`) windows (n2 :2016, n4 :2018). Those included the 21 ball integrals this module skips. The driven solves take ≈ 17–21 s each and the mesh 36 s.
- **W1 (`-n 2`),** predicted ≈ 170 s:

  `scripts/testing/run_and_log.sh POST-8-step2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py tests/validation/test_post8_owner_rule_ans2.py -v -s --tb=short'"`
- **W2 (`-n 4`),** predicted ≈ 140 s: the same command with `-n 4`, labelled `POST-8-step2-W2`.
- **W3 (`-n 1`),** unmeasured, predicted ≈ 260 s: the same command with `-n 1`, labelled `POST-8-step2-W3`.
- **W4,** the real-build re-gate under rule (c), predicted ≈ 30 s:

  `scripts/testing/run_and_log.sh POST-8-step2-W4-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py -v -s --tb=short'"`
- `tests/post/test_evaluation_collective_guard.py` is `evaluation.py`'s pre-existing gate (OPS-40/55; 6 passed / 1 skipped complex, `20260921T111802Z_OPS-55-callers-complex.log:260–272`).
- After W4, run `git status --porcelain examples/` as its own call. It must be empty.

**Traps:**
- Use the buffer-form `Allreduce` on int64 with `MPI.MIN`. The lowercase `allreduce` compares objects. Never decide ownership rank-locally.
- Filter to owned cells: `create_unit_cube` ghosts by default, and `bb_tree` may include ghost cells.
- Read `original_cell_index`'s shape in-container (K0 prints it).
- The rule is collective, which is why `owner` rides the digest allgather. Raise on every rank (OPS-40).
- Loading `ans:2` also loads EX-34's script and patches `build_four_port_sweep` inside the mass-averaged module. Import the builder directly, and never call `ans2.main()`: it writes `metrics.json`, `COMPARISON.md` and the private comparison.
- `cell_tags.values` and `assemble_scalar` are rank-local. `_phantom_power` → `mean_sar` already reduces, so do not re-sum.
- Use the complex build. The DG0 index field is complex there: compare `.real` and assert `.imag == 0`.
- Rule (i): no edit to either module after its last window. Rule (g): `-s`.
- A W3 kill at 590 s is not a result. Run `pgrep -c python3` first, and never lengthen the timeout.
- If the oracle's incident set disagrees with DolfinX's at the vertex, print both: that is a collision-tolerance finding. Never tune the 1e-12.

**Scope (step 3, not queued here):**
- adoption at the census sites (`ans:2`'s `_pointwise_sar`, `post/sar.py:354`, `phantom_fields.py:163` and the 18 heuristic candidates);
- flipping the default;
- classifying the census;
- re-issuing `ans:2`'s Row 1 as a cell or mass average;
- re-reading the 09-26 `ANS-2` "within 6.3 %", which keeps its `-n 4`-only caveat.

No band moves, and no integral is touched.

**Status it can move:**
- §7 `POST-8` stays 🟡. Its state line gains "step 2 done: the rule implemented opt-in and gated at `-n 1 / 2 / 4`", and step 3 becomes queueable.
- The 2026-09-27 known-issues entry (`known-issues.md:573–576`) stays OPEN; its "Verified at" row gains step 2's reading.
- The review that queues step 3 (which writes "→ ✅") amends the Done-when under rule (j) to name `-n 1 / 2 / 4` and the averaged public Row 1.

**Negative result: report and stop, never widen a bound.**
- **(A2) misses at `-n 2` / `-n 4`** with A0 / A1 green: step 1's reading is contradicted. Annotate the known-issues Cause row, park on `attempt/*`, and mark the item BLOCKED in the same commit (rule (d)).
- **(A2) misses at `-n 1` only:** record it and mark BLOCKED with W1 / W2 green; the weekly re-poses the `-n 1` clause.
- **(A0) or (A1) fails:** the mesh depends on the width. Open a new known-issues entry and stop.
- **K1–K3 red:** revert `src/` and park.
- **The guard module or W4 red:** the legacy path moved. Revert `src/` and stop.
- **C1a fails:** the fixture drifted. Report it and stop.
- **W3 killed:** W1 / W2 stand. Mark BLOCKED, with the unblock condition "re-brief `-n 1` on P1 only".
- There is no PREDICTED RED line.

### §9 item 188 — 2026-10-02 (lane 1)

**Outcome (2026-10-03 22:30 CDT slot): DONE — 3a gated at `-n 1 / 2 / 4`, 19 passed in each window, no `src/` change.**

- **Windows:** W1 `-n 2` `20261003T033225Z_POST-8-step3a-W1.log` (272 s), W2 `-n 4` `20261003T033703Z_POST-8-step3a-W2.log` (231 s), W3 `-n 1` `20261003T034100Z_POST-8-step3a-W3.log` (338 s; predicted ≈ 290 s). `git status --porcelain examples/` empty.
- **(A0):** 199 920 cells; P1 power rel 4.22e-12 / 6.45e-12 / 1.07e-11. **(A2):** 1.37e-10 / 1.30e-10 / 3.03e-10; owner `original_cell_index` 150135 / 162464 / 142550 / 169932 at every width. **(R):** 2.1e-15 / 2.1e-15 / 2.3e-15, route (ii) Im 0. **(M):** max \|c\| + a₁g 2.120e-2 m < 0.03; a₁g 6.20e-3 m < 0.04. **(G):** exactly the planted entry in each family at every width.
- **(W)** max rel vs `-n 2`, `-n 4` / `-n 1`: raw 2.21e-10 / 1.66e-10, cell 2.15e-10 / 1.84e-10, 1 g 1.07e-10 / 8.59e-11 (≤ 1e-9). Predicted 1 g ≤ 1e-10 **missed at `-n 4`** (printed, rule (e)).
- **C1b:** 0.1293 (`-n 2`) / 0.1485 (`-n 4`); `-n 1` printed 0.0343 / 0.0882.
- **Ratios (identical to six digits across widths):** cell/raw 0.639 … 1.596 — **predicted [0.8, 1.25] missed** on 6 of 16 entries, all at off-leg points; 1 g/raw 1.000 … 2.993 (P2 +x). Both bear on the weekly's cell vs 1 g choice for the published row.

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.* Step 3 as scoped above (`POST-8.md:115–120`) does not fit one slot, and three of its parts are blocked:
- re-issuing the tracked Row 1 needs `ans:2`'s `main()`, which does not honour `FEM_EM_ANS_SKIP_PRIVATE` (the 10:30 ruling: do not run it);
- the census classification is the weekly's;
- the 09-26 re-read touches private data.

This item takes **3a**: the averaged Row 1, computed inside a gate module and proven width-invariant. **3b** is the census-site adoption and the tracked Row 1.
- **Executor:** implementer. **Tier:** heavy, complex. **Cost:** ≈ 4 measured / ≈ 45 predicted min.
- **Depends on:** nothing (step 2 is done, item 177). No `src/` change.

**Module `tests/validation/test_post8_step3a_averaged_row1.py`** (`complex_only`).
- **Fixture:** built as item 177's `case` fixture (`test_post8_owner_rule_ans2.py:114–148`). `ans:2` is loaded by importlib only for `_pointwise_sar`, `_phantom_power` and its constants. **Never `main()`.**
- **Imported, never restated:**
  - `STEP1_N1CURL_SAR`, `A2_RTOL` and `C1B_MIN_SEPARATION` (step 2's module);
  - `LOGGED_P1_ROW`, `TRACKED_P1_POWER_W` and `POWER_RTOL` (the probe);
  - `ONE_GRAM_KG` and `QUADRATURE_DEGREE` (`test_mass_averaged_sar_standard_masses`);
  - `PHANTOM_RADIUS` and `PHANTOM_HEIGHT` (`test_birdcage_port_tags`).
- **Readings: three families, each 4 × 4** (drives P1–P4 × points +x, +y, −x, −y):
  - **raw:** the owner-rule SAR.
  - **cell:** the owner cell's average of σ|E|²/(2ρ), by two routes, both exact because |E|² is quadratic for degree-1 N1curl:
    - (i) `e_complex` at the 4-point degree-2 rule (a = 0.5854101966249685, b = 0.1381966011250105), evaluated on the owning rank and shared by buffer `Allreduce` SUM;
    - (ii) σ·inner(E,E)·χ/(2ρ) assembled over a DG0 indicator χ of the owner cell and divided by ∫χ, with `metadata={"quadrature_degree": 2}` and both integrals reduced.
  - **1 g:** `mass_averaged_sar` centred on each point, radius `averaging_ball_radius(ONE_GRAM_KG, ρ)`, `QUADRATURE_DEGREE`, with ρ from `build_density_field` as at `test_birdcage_sar_mass_averaged.py:254–260`.
- **npz guard:**
  - At `-n 2`, rank 0 writes `/workspace/logs/post8/step3a-n2.npz` (gitignored) with `raw`, `cell`, `one_gram`, `owner_original` and `cells`, then `Barrier`.
  - At `-n 4` / `-n 1` the module loads it. If it is missing, `pytest.fail`; never recompute.

**Asserted.**
- **(A0)** 199 920 cells, and P1 phantom power within `POWER_RTOL` of `TRACKED_P1_POWER_W` (step 2 measured 4.91e-12 / 8.93e-13 / 1.07e-11).
- **(A2)** raw within `A2_RTOL` of `STEP1_N1CURL_SAR` (step 2 measured ≤ 3.03e-10).
- **(R)** The two cell-average routes agree to ≤ 1e-12 relative, and route (ii)'s |Im| ≤ 1e-14·Re.
- **(M)** Every 1 g ball lies inside the phantom: |c| + a₁g < `PHANTOM_RADIUS` and a₁g < `PHANTOM_HEIGHT`/2.
- **(W)** At `-n 4` and `-n 1`:
  - `owner_original` and `cells` equal the npz exactly;
  - each of the 48 values is within 1e-9 relative of its `-n 2` entry (the weekly's 1e-9).

  *Predicted, printed beside:* raw and cell ≤ 5e-10, 1 g ≤ 1e-10. This adopts step 2's executor hypothesis (attempts.md 2026-10-02T20:15Z): the same 16 raw SARs measured 2.98e-10 across widths, so the 1e-11 prediction is retired.
- **(G) Planted control on the comparator:** one entry per family multiplied by (1 + 1e-8) is flagged, and only that entry.

**Negative control (rule (e)).**
- *Asserted* at `-n 2` and `-n 4`: C1b — `ans:2`'s CG1 P1 row against the other width's `LOGGED_P1_ROW` has min |rel| ≥ `C1B_MIN_SEPARATION` (item 177 measured 0.1293 / 0.1485 on this fixture).
- *Printed only:* C1b at `-n 1`; the cell/raw and 1 g/raw ratios (cell/raw predicted in [0.8, 1.25]).

**Windows.** Bare, foreground, Bash `timeout: 660000`, in the order W1 → W2 → W3. Basis: step 2's windows took 195 / 191 / 228 s; this adds 16 ball integrals and 32 cell integrals, and W1 pays the JIT for four new ball centres.
- **W1 (`-n 2`, ≈ 260 s):**

  `scripts/testing/run_and_log.sh POST-8-step3a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_post8_step3a_averaged_row1.py -v -s --tb=short'"`
- **W2 (`-n 4`, ≈ 240 s):** the same command with `-n 4`, label `POST-8-step3a-W2`.
- **W3 (`-n 1`, ≈ 290 s):** the same command with `-n 1`, label `POST-8-step3a-W3`.
- **Then, as its own call:** `git status --porcelain examples/` must print nothing.

**Traps:**
- **Private data:** never `main()`; never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`; no AED number anywhere.
- **Ownership** is decided only by buffer `Allreduce` MIN, never rank-locally.
- **Reductions:** `assemble_scalar` is rank-local, but `mean_sar` / `mass_averaged_sar` already reduce — do not sum them again.
- **The ball centre** must be plain floats (`sar.py:257`).
- **A killed compile** leaves an orphan `.c`: run the rubric's cache sweep before a re-run.
- **Rule (i):** a re-run of W1 means W2 and W3 re-run too. `-s`, no `-k`, no pipes.

**Scope.** Not 3b: adoption at the census sites, flipping the default, switching `ans:2`'s `_pointwise_sar`, the tracked Row 1, and the 09-26 re-read. No tracked example file changes.

**Negative result.** Report and stop; never loosen the 1e-9.
- **(A0), (A2) or C1b red:** fixture drift. Annotate known-issues 2026-09-27, park, and mark BLOCKED.
- **(R) red:** a cell-average code bug. Fix the code, never the 1e-12, and re-run W1.
- **(W) red:** record the per-entry readings, annotate known-issues, and mark BLOCKED. The weekly re-poses the bound.
- **A killed window:** the passed windows stand; unblock by re-running the missing width.
- There is no `PREDICTED RED` line.

**Status it can move.**
- The §7 state line gains "step 3a done"; the row stays 🟡 (Done-when amended by this review under rule (j)).
- Known-issues 2026-09-27's "Verified at" gains 3a.
- The weekly is handed the cell vs 1 g choice for the published row, plus 3b's two blockers: `ans:1`–`ans:3` lack the skip-private knob (ops work), and the census classification is the weekly's.

### §9 item 526 — 2026-10-18 (`POST-8` step 3b, owner-cell rule adopted; lane 1)

*Drafted by a read-only drafter for the 2026-10-18 03:00 daily review. The review checks every command against its source before queueing.*

**Licence.** The 2026-10-17 weekly classified 3b (`PROJECT_PLAN.md` ≈ 4711–4717). The published `ans:2` Row 1 is the **1 g mass average**. Every other point-sampling site of the census (`post/sar.py`'s point sampler, `phantom_fields.py:163`, and the 18 heuristic candidates) adopts the **owner-cell rule** as its default. 3b is lane 1 and "queueable now". **No band, tolerance, Done-when or target is set or moved.** Templates: item 177 (step 2's known-answer module) and item 188 (3a's npz width guard).

**Split: 3b does not fit one slot, so this item is 3b-i only.**
- **3b-i (this item):** the two `src/` sites. `point_sar` defaults to the rule. `phantom_fields.py:163` is asserted owner-compliant by construction.
- **3b-ii (next):** the 17 non-`ans:2` heuristic candidates, i.e. 12 examples and 5 tests (census, `20260927T112115Z_POST-8-step1-n2.log:1975–2013`). It splits by what each site samples:
  - **(α) a DG1 split of the N1curl field** through the legacy path. Pass the rule: `time_harmonic/03:227–228`, `06:157–158`, `07:170–171`, `09:161–162`; `test_lossy_sphere_degree2.py:176–177`; `test_lossy_sphere_fullwave.py:218–219`.
  - **(β) a Lagrange interpolant of a discontinuous field.** The rule cannot remove the interpolant's own partition dependence (step 1's mechanism). These sites need a switch to the source field, as `POST-4` step 3 did at `mri/01_coil_phantom_fields.py:478–487`: `magnetostatics/01:273` (`B_lag`), `test_time_harmonic_smoke.py:206`, `test_coil_phantom_bfield_metrics.py:149, :239`.
  - **(γ) source fields:** the remainder (`magnetostatics/02:269`, `04:200`, `05:123`, `06:212`; `ports/01:303`, `12:503–506`; `time_harmonic/04:227–228`; `test_straight_wire.py:289`; `test_coil_phantom_bfield_metrics.py:134, :204`).
  - The site lists are the drafter's code reading; 3b-ii's drafter verifies each field type. The tests can run as an implementer item. The examples are `example-runner` work, with guide re-sync.
- **3b-iii (blocked as posed):** re-issue the tracked Row 1. Its prerequisite is `ans:2` honouring `FEM_EM_ANS_SKIP_PRIVATE`.
  - Today it does not (`grep` finds the knob only in `ans:4` (`04_…py:196`), `ans:6`, `ans:7` and `ans:8`).
  - With `aed_results/` on disk, `main()` rewrites `COMPARISON_private.md` (`02_birdcage_coil_driven_sar_10MHz.py:1056–1057`), which `EX-57.md:40` forbids.
  - 3b-iii adds the knob on `ans:4`'s pattern, makes Row 1 the 1 g average (3a's family), regenerates `metrics.json` / `COMPARISON.md` / the guide, and retires known-issues 2026-09-27. `POST-8` → ✅ needs 3b-ii and 3b-iii both.

**`POST-8` step 3b-i: `point_sar` defaults to the owner rule and `phantom_fields.py:163` is asserted owner-compliant. Both are gated width-invariant at `-n 1 / 2 / 4` on step 2's unit-cube fixture, and `point_sar`'s five gate modules are re-run green (rule (c)).**
- **Executor:** implementer.
- **Tier:**
  - W1–W4 standard (`-n 2 / 4 / 1`, plus real build `-n 2`). Basis: step 2's W4-real ran `tests/environment` + guard + `test_point_owner_rule.py` in 49 s (`20261002T201406Z_POST-8-step2-W4-real.log`). Predicted ≈ 40 s each.
  - W5a heavy by ceiling (`-n 2`). Measured 53 s for the two mass-averaged modules plus `test_lossy_sphere_sar.py` (`20260907T003812Z_MAT-4-step4-src-regression.log:69`) and 81 s for `test_pec_sphere_hole.py` + `test_imp1_pec_sphere_lossy_host.py` (`20261001T034025Z_IMP-1-land-W3.log:518`). Predicted ≈ 140 s.
  - W5b heavy (`-n 2`). Measured 194 s for `test_birdcage_sar_integral.py` + `test_birdcage_sar_map.py` (`20260903T003309Z_WF-6-step3h.log:9712`) and 117 s for `test_birdcage_sar_fine_phantom.py` (`20260903T020607Z_WF-6-step3f-prime.log:4970`). Predicted ≈ 320 s.
- **Class:** lane 1. **Cost:** ≈ 3.5 measured (lane-1 median) / ≈ 35 predicted slot-minutes. Windows ≈ 10.5 min, code ≈ 15, commit and journal ≈ 6. **If W5a has not started by minute 40, park** (rule (d)).
- **Depends on:** nothing in §9. Steps 2 and 3a are done (items 177, 188). On disk: none. This item writes its own npz and never reads `logs/post8/step3a-n2.npz`.
- **Family:** `POST-8` step 3b-i; attempt 1.

**What the code is today.**
- **`post.sar.point_sar`** (`src/fem_em_solver/post/sar.py:332–365`) calls `evaluate_vector_field_parallel(e_real, …)` / `(e_imag, …)` with no `owner` (`:354–355`), so it takes the legacy path. That path is `evaluation.py:379–388`: `links[0]` on every rank that finds the point, and the last gathered rank wins (`:413–416`). The rule branch is `:389–396` → `point_owner_cells` (`:133–221`). `OWNER_RULES` is at `:14`.
- **Callers that import `point_sar`** (the drafter's `grep`):
  - five gate modules: `test_mass_averaged_sar.py`, `test_mass_averaged_sar_standard_masses.py`, `test_imp1_pec_sphere_lossy_host.py`, `test_birdcage_sar_map.py`, `test_birdcage_sar_fine_phantom.py`;
  - three examples: `materials/02:171` and `materials/03:104` (cell centroids, `02_lossy_sphere_sar.py:57`, `03_…py:100–104`) and `mri/02:325–331` (origin, imposed uniform field).
  - None uses `create_submesh` (`grep`), so the rule's no-index raise (`evaluation.py:99–106`) cannot fire.
- **`phantom_fields._evaluate_on_cells`** (`phantom_fields.py:149–202`; batch `.eval` at `:163`, fallback at `:181`). It is called with `cells = _tagged_cells(...)`, which are owned cells only (`:73–77`), and `points = _cell_centroids(mesh, cells)` (`:80–92`), via `:294–303` and `:373–381`. A centroid lies strictly inside one cell, so its incident set is that cell alone and the rule's owner is the sampled cell. **The site is compliant by construction; this item asserts that and changes no code there.** Its production callers pass source fields (`mri/01:394`, `:415`; `test_quicklook_report.py:178–179`).
- **The time-harmonic split** (`core/time_harmonic.py:745–754`): `e_real` / `e_imag` are DG1 interpolants of the degree-1 N1curl phasor. That is exact, and it is discontinuous at every facet, so a facet or vertex sample depends on the cell chosen.

**Code. Edit only `src/fem_em_solver/post/sar.py` (`point_sar` only) and add `tests/post/test_post8_step3b_point_sar_owner.py`.** No edit to `evaluation.py`, `phantom_fields.py`, any example, the `ans:2` script, or any existing test module. This is a `src/` change under rule (c): disclose it, and re-run its pre-existing gates green in the same slot (W4, W5a, W5b).
- `point_sar` gains a keyword-only `owner: str | None = "lowest_original_cell"`, forwarded as `owner=owner` to both calls at `:354–355`.
  - The docstring names `POST-8` step 3b and says that `owner=None` is the legacy path, byte for byte.
  - An unknown value raises in `evaluate_vector_field_parallel` (`:288–291`); add no check of your own.
- Nothing else changes in `point_sar`: signature order, `comm` default, the raise on unevaluated points, and the return.

**Module `tests/post/test_post8_step3b_point_sar_owner.py`** (both builds, apart from (W)).
- **Imported, never restated:**
  - `POINTS`, `LABELS`, `CLOSED_FORM_TOL`, `_mesh`, `_original`, `_oracle` and `_bitwise_equal_across_ranks` from `tests/post/test_point_owner_rule.py` (`:35`, `:43`, `:46`, `:51`, `:55`, `:62`, `:84`);
  - `WIDTH_RTOL` from `tests/validation/test_post8_step3a_averaged_row1.py:88`, **imported inside the (W) test only**. That module's import graph pulls in the birdcage builders and is `complex_only`.
  - `_cell_centroids` / `_evaluate_on_cells` from `fem_em_solver.post.phantom_fields`.
- **Fixture.** `create_unit_cube(comm, 4, 4, 4)`, built with `_mesh`. Two real fields are interpolated into N1curl degree 1 and then into DG1, exactly as `time_harmonic.py:745–754` splits:
  - `f_r = (sin 3y + z, x·cos 2z, x·y)`;
  - `f_i = (y², sin 2x, cos y + z)`;
  - σ = 0.5, ρ = 1000.
  - Both are non-polynomial, so the interpolant's normal component jumps across facets.
  - The DG0 original-index field of item 177's K3 (owned + ghost, `scatter_forward`) is rebuilt here for (N1).
- **Asserted:**
  - **(D0) The default is the rule, and the legacy path is reachable.** `point_sar(e_r, e_i, POINTS, sigma, rho)` is bitwise equal to the same call with `owner="lowest_original_cell"`. The call with `owner=None` is bitwise equal to σ|E|²/(2ρ) built by hand from `evaluate_vector_field_parallel(..., owner=None)` on both parts. Every rank holds the same bytes (`_bitwise_equal_across_ranks`).
  - **(K) The known answer.** At the three points, `point_sar` equals σ|E_c|²/(2ρ) to ≤ `CLOSED_FORM_TOL`. E_c is evaluated in the oracle's minimum cell (`_oracle`): the rank whose owned `_original(msh)[:size_local]` holds that index calls `.eval` in it, and the result is shared by buffer `Allreduce` SUM. That route is independent of `point_owner_cells`' collision path.
  - **(P) `phantom_fields.py:163` is compliant.**
    - Take all owned cells. Allgather their `_cell_centroids` and original indices; that is 384 points at every width.
    - `point_owner_cells(msh, centroids)` returns each centroid's own original index, for all 384.
    - On each rank, `_evaluate_on_cells(e_r, local_centroids, local_cells)` equals `evaluate_vector_field_parallel(e_r, centroids, owner="lowest_original_cell")` bitwise on that rank's rows.
  - **(W) Width invariance** (`complex_only`).
    - At `-n 2`, rank 0 writes `/workspace/logs/post8/step3b-n2.npz` (gitignored) with the 3 SARs and the 3 owner original indices, then `Barrier`. **Never touch `step3a-n2.npz`.**
    - At `-n 4` / `-n 1` the module loads it; if it is missing, `pytest.fail`, never recompute.
    - The owner indices must equal the npz exactly, and each SAR must lie within `WIDTH_RTOL` (1e-9) relative of it.
    - *Predicted, printed beside:* ≤ 1e-15. There is no solve, and interpolation is cell-local.
- **Printed:** the 3 SARs to `.17e`; the owner index per point (item 177 measured 217 / 223 / 126); the max/min of SAR over the vertex's 24 incident cells (predicted > 1.01, not asserted).

**Negative control.**
- **(N1) Asserted (rule (e); the same comparison on the same fixture and code path):** the legacy cell, read through the DG0 index field with `owner=None`, differs from the owner.
  - At the vertex at `-n 4` (128 vs 126, `20261002T200654Z_POST-8-step2-W2.log:166`) and at `-n 1` (131 vs 126, `20261002T201013Z_POST-8-step2-W3.log:79`).
  - At the facet at `-n 4` (232 vs 223, W2 `:165`).
  - At `-n 2` the legacy cell equals the owner (W1 `:106–108`). That is printed and not asserted, so W1 carries no (N1) assertion.
- **(N2) Predicted, printed:** at the widths where (N1) fires, the legacy `point_sar` at the vertex differs from the rule's value by ≥ 1 %.
- **The −19 % … +112 % is not cited as asserted.** It is `ans:2`'s CG1 route, and step 2's (A1) shows each `ans:2` named point has exactly one owned containing cell. On that fixture, flipping `point_sar` is therefore a no-op; the dependence there is the CG1 interpolation at `02_…py:327–333`, which 3b-iii removes.

**Windows.** Foreground, Bash `timeout: 660000`, `-s`, no `-k`, no pipes. **The harness call is the whole Bash command.** Run W1 → W2 → W3 → W4 → `pgrep` → W5a → `pgrep` → W5b.
- **W1** (`-n 2`, complex, writes the npz):
  `scripts/testing/run_and_log.sh POST-8-step3b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py tests/post/test_post8_step3b_point_sar_owner.py -v -s --tb=short'"`
- **W2** (`-n 4`): the same command with `-n 4`, label `POST-8-step3b-W2`.
- **W3** (`-n 1`): the same command with `-n 1`, label `POST-8-step3b-W3`.
- **W4** (real build, `-n 2`; rule (c) for the real path; (W) skips):
  `scripts/testing/run_and_log.sh POST-8-step3b-W4-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py tests/post/test_post8_step3b_point_sar_owner.py -v -s --tb=short'"`
- **Before W5a and before W5b:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W5a** (`point_sar`'s sphere gates, `-n 2`):
  `scripts/testing/run_and_log.sh POST-8-step3b-W5a-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mass_averaged_sar.py tests/validation/test_mass_averaged_sar_standard_masses.py tests/validation/test_pec_sphere_hole.py tests/validation/test_imp1_pec_sphere_lossy_host.py -v -s --tb=short'"`
- **W5b** (`point_sar`'s birdcage gates, `-n 2`):
  `scripts/testing/run_and_log.sh POST-8-step3b-W5b-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_birdcage_sar_map.py tests/validation/test_birdcage_sar_fine_phantom.py -v -s --tb=short'"`
- Then, as its own call, `git status --porcelain examples/` must print nothing.

**Traps already paid for.**
- **Do not flip `evaluate_vector_field_parallel`'s own default.** ≈ 30 test modules, `core/solvers.py:766`, `ports/gap_voltage.py:289` and `study.py:1975` ride it. They are not in the weekly's census, and changing them is out of scope.
- **The rule is collective.** Every `point_sar` caller already passes one identical point list (OPS-40/55); `(P)` must allgather the centroids before calling the rule. Ownership is decided only by buffer `Allreduce` MIN.
- **Complex build:** the DG0 index field is complex; compare `.real` and assert `.imag == 0`.
- **The import graph:** importing `test_post8_step3a_averaged_row1` at module level would drag the birdcage builders into the real-build W4. Import `WIDTH_RTOL` inside the `complex_only` (W) test.
- **The npz:** rank 0 writes and every rank waits at `Barrier`. Use a new file name; 3b-iii may still need `step3a-n2.npz`.
- **Rule (i):** an edit after W1 starts means W1–W5b again.
- **A killed window:** run `pgrep -c python3` first; at 0, run the rubric's cache sweep as `POST-8-step3b-cache-sweep`. A kill is not a result; never lengthen a timeout.
- **W5a / W5b last ran 09-07 … 10-01** (the cited logs). A red there may predate this edit (see Negative result).
- Never load or run `ans:2`; never open `docs/private/` or `aed_results/`; no AED number anywhere.

**Scope.** `point_sar`'s default and an assertion that `phantom_fields.py:163` complies. **Not claimed:**
- 3b-ii's 17 candidates;
- the (β) source-field switches;
- the three examples that import `point_sar` (predicted digit-identical: they sample centroids or an imposed uniform field; 3b-ii's example window re-runs them);
- the `src/` legacy callers outside the census;
- `ans:2`'s knob and the tracked Row 1;
- the 09-26 `ANS-2` re-read (private);
- any band.

**Negative result.** Report and stop. Park on `attempt/POST-8-3b-i-*`, mark item 526 BLOCKED in the same commit (rule (d)), and never loosen 1e-9.
- **(D0) or (K) red:** the `src/` edit is wrong. Revert it on the branch.
- **(W) red:** record the per-point readings and annotate known-issues 2026-09-27; the review re-poses the bound.
- **(N1) red** (legacy equals owner at `-n 4` / `-n 1`): the partition has moved since step 2. The control has lost its separation on this fixture, which is not a defect of the rule; the review re-poses it.
- **(P) red:** `phantom_fields.py:163` is not compliant. Open a new known-issues entry; 3b-ii is unaffected.
- **W4 / W5a / W5b red:** the change is not admitted. `git stash` the `sar.py` edit and re-run only the red module, with the same command shape and label `…-baseline`. If the baseline is also red, the failure predates this item: open a known-issues entry. Either way the item stays BLOCKED.
- **Status 124:** not a result; park with the measured time.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** `POST-8` stays 🟡.
- Its Done-when (`PROJECT_PLAN.md:2670`) needs the rule adopted by every classified census site *and* the tracked Row 1 re-issued; this item delivers two sites of twenty and no Row 1.
- On green, the §7 state line gains: "**Step 3b-i done** <date> (§9 item 526): `point_sar` defaults to the owner rule (`owner=None` keeps the legacy path); `phantom_fields.py:163` owner-compliant by construction, asserted on 384 centroids at `-n 1 / 2 / 4`; width-invariant to <max> (1e-9 asserted); five gate modules re-run green (<W5a log>, <W5b log>)."
- The stale clause "waits on the weekly's census classification" is replaced by "classified 10-17; 3b-ii (17 candidates) and 3b-iii (`ans:2` skip-private knob, 1 g Row 1) open".
- Known-issues 2026-09-27's "Verified at" gains 3b-i, and the entry stays OPEN until 3b-iii.
- The review that queues 3b-iii amends the Done-when's "cell or mass average (… the weekly picks …)" to "the 1 g mass average" (rule (j)).

**Execution order.** Code → W1 → W2 → W3 → W4 → `pgrep` → W5a → `pgrep` → W5b → `git status --porcelain examples/` → commit (`sar.py`, the new module, harness logs, the §7 sentence, the known-issues line) → journal.

**Spawn line (if delegated):** "Use the implementer agent to execute `POST-8` §9 item 526 (step 3b-i) per `docs/planning/chunks/POST-8.md`. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W1, W2, W3, W4, W5a and W5b in order, in the foreground, with `timeout: 660000`, pytest with `-s`, and container `timeout -k 30`. Before W5a and W5b, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. Edit only `point_sar` in `src/fem_em_solver/post/sar.py`, and add `tests/post/test_post8_step3b_point_sar_owner.py`. Make no edit to `evaluation.py`, `phantom_fields.py`, any example, the `ans:2` script or any existing test. Never load `ans:2`, and never open `docs/private/` or `aed_results/`. You are the executor: do not spawn agents, and never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/POST-8-3b-i-*`, mark item 526 BLOCKED in the same commit, and report."

### Review ruling 2026-10-18 03:00 (item 526, step 3b-i)

Accepted as 3b-i only. Spot-check: `src/fem_em_solver/post/sar.py:332` is `def point_sar(` (re-read). The weekly's "`post/sar.py`'s `_pointwise_sar`" means `point_sar`; `_pointwise_sar` is `ans:2`'s own function.
1. **Yes.** At `phantom_fields.py:163`, "owner-compliant by construction, asserted" counts as adopting the rule. The code is not edited.
2. **`core/solvers.py:766`, `ports/gap_voltage.py:289` and `study.py:1975` are out of 3b-i.** They are not in the weekly's classified census and are named for the weekly.
3. **The (β) interpolant sites** are 3b-ii's question and the weekly's.
4. **No Done-when amendment.** The item holds `POST-8` 🟡 (rule (j)). The drafter's "1 g mass average" wording is the weekly's to adopt.
5. **3b-iii stays unqueued.** Running `ans:2`'s `main()` would rewrite `COMPARISON_private.md`, and the 10-02 10:30 ruling forbids running it. Adding a skip-private knob is named for the weekly.
6. Noted for 3b-iii: `metrics.json["pointwise"]` is read by `test_ans9_far_side_sar.py:88, :220–228`.

**Outcome 2026-10-18 (07:30 CDT slot, attempt 1): done, green first try.** `point_sar` gained keyword `owner="lowest_original_cell"` forwarded to both evaluations; new `tests/post/test_post8_step3b_point_sar_owner.py`. (D0) bitwise at every width; (K) rel 0.0; (P) 384 / 384 at `-n 1 / 2 / 4`; (W) 3.611e-16 / 4.087e-16; (N1) fired at `-n 4` (facet 232 vs 223, vertex 128 vs 126) and `-n 1` (vertex 131 vs 126), `-n 2` legacy = owner (printed); (N2) 11.9 % / 12.2 %; vertex incident-set SAR max/min 1.6557. W5a 14 passed (60 s), W5b 123 passed (156 s) (W5b log held out of git in the gitignored `logs/leak-held/`: `check_private_leak.py` flagged line 3948 — our own `test_birdcage_sar_fine_phantom.py` best-approximation print, a coincidental digit-prefix collision; hook not bypassed, log unedited, EX-112 precedent; ask: may it land?). Logs `20261018T125312Z_POST-8-step3b-W1.log`, `20261018T125322Z_POST-8-step3b-W2.log`, `20261018T125326Z_POST-8-step3b-W3.log`, `20261018T125330Z_POST-8-step3b-W4-real.log`, `20261018T125346Z_POST-8-step3b-W5a-regate.log`, `20261018T125455Z_POST-8-step3b-W5b-regate.log`. The (P) equality is checked as bytes of the rank's own rows.

### Review ruling 2026-10-18 10:30 (item 526's held W5b log)

- **3b-i's DONE stands.** Its anchors, (D0) / (K) / (P) / (W) and the (N1) / (N2) controls, are asserted in W1–W4, and those logs are on `main`. W5a, also on `main`, re-gates the step-2 / 3a modules. W5b is the rule (c) re-gate of `point_sar`'s other gate modules. Its harness row is in `test-results.md:4682` (Status 0, 156 s). Its log sits unchanged in the gitignored `logs/leak-held/`, read at 10:30 (299 629 B). No status rests on W5b alone: `POST-8` stays 🟡 and no glyph moved.
- **Landing the log is the operator's call**, as for `EX-112`'s two held logs (2026-10-14). The options are a guard allowlist entry for that line, or a ruling that the collision may land. The log is never edited, trimmed or committed past the hook. Named on the dashboard beside `EX-112`.
- **For 3b-ii:** the flagged print is in `test_birdcage_sar_fine_phantom.py`. Any later window that re-runs that module will hit the guard again. Briefs that re-gate it carry the hold-and-report instruction.

### §9 item 536 — 2026-10-18 (`POST-8` step 3b-ii-t1; lane 1)

*Drafted read-only for the 2026-10-18 10:30 daily review; the review reviews it (ruling below).*

**Outcome 2026-10-18 (16:30 CDT slot, attempt 1): done, green first try.** Four calls on the owner rule; the smoke test samples `fields.e_imag` directly (CG1 interpolant deleted), prints at 6 significant digits. New `tests/post/test_post8_step3b_ii_th_sites.py`. (X) W1 25 passed in 72 s (all five modules at their bands). (S) 25/25 probes single-incidence at degree 1 and 2 at every width, rule = legacy bitwise. (O) owner indices (25 probes + V) equal the `-n 2` npz at `-n 4` / `-n 1`. (W) max spread d1 1.326e-12 / d2 8.532e-12 at `-n 4`, 1.877e-12 / 1.108e-11 at `-n 1` (predicted ≤ 5e-10). (N0) V has 26 incident cells, max/min `|E|` 1.007059 (asserted > 1 + 1e-6; the predicted > 1.01 missed). (N1, printed) legacy cell at V differs from owner 4 at every width (50 / 580 / 491 at `-n 2 / 4 / 1`). (P-smoke) 1.527e-14 (predicted ≤ 1e-6); `-n 2` writes the full-precision values to the gitignored `logs/post8/step3b-ii-t1-smoke-n2.npz` for that print. No leak-guard hit. Logs `20261018T213254Z_POST-8-step3b-ii-t1-W1.log` (74 s), `20261018T213419Z_POST-8-step3b-ii-t1-W2.log` (8 s), `20261018T213438Z_POST-8-step3b-ii-t1-W3.log` (13 s).

**Licence.** The 10-17 weekly's 3b classification, plus the 526 brief's 3b-ii split (`POST-8.md:229–233`). No band, tolerance, Done-when or target is set or moved. Template: item 526.

**Why 3b-ii's test half is split.** Three of the 526 brief's γ sites are mislabelled. It calls them "source fields" needing no change. But `test_coil_phantom_bfield_metrics.py:134`, `:204` sample `b_dg0` (a DG0 interpolant), and `test_straight_wire.py:289` samples `compute_b_field()`, DG1 of curl A (`solvers.py:644–654`). Both fields are discontinuous, so the owner rule applies. A 1e-9 width bound is out of reach on magnetostatic solves: MAG-18 measured the LU's cross-width floor at 7.28e-8 (`test_straight_wire.py:113–122`). The straight-wire module also ran 371 s at `-n 2` (`20260828T050256Z_MAG-20-step1-anchor-module.log.gz`). So this item is **3b-ii-t1**, the time-harmonic sites in the complex build. 3b-ii-t2 (magnetostatic, real build) is described at the end and needs a ruling first.

**`POST-8` step 3b-ii-t1: the three time-harmonic test sites sample under the owner rule. The rule is asserted to be a bitwise no-op at single-incidence probes and width-invariant at `-n 1 / 2 / 4`, and the edited modules plus their downstream gates are re-run green.**
- **Executor:** implementer.
- **Tier:**
  - W1: heavy by ceiling (`-n 2`, 400 s), predicted ≈ 110 s. Basis:
    - `test_lossy_sphere_fullwave` + `test_lossy_sphere_degree2`: 27 s (`20260825T033221Z_OPS-25-gate-green.log.gz`);
    - `test_geometry_floor_discriminator`: 28 s (`20260813T213156Z_GEO-14-step1-discriminator.log.gz`);
    - `test_degree2_energy_mechanism`: 10 s (`20260819T183425Z_TH-12-step3-warm.log.gz`);
    - `test_time_harmonic_smoke`: 8 s (`20260819T201005Z_POST-5-step4-smoke-final.log.gz`), not re-measured since the module grew.
  - W2: standard (`-n 4`), predicted ≈ 30 s.
  - W3: heavy by ceiling (`-n 1`, 300 s), predicted ≈ 90 s. Basis: the smoke node plus the sigma hook ran 65 s at `-n 1` (`20261005T050254Z_MAT-5-step3-W0-hook.log.gz`).
- **Class:** lane 1. **Cost:** ≈ 6.1 measured / ≈ 35 predicted slot-minutes. **If W1 has not started by minute 35, park.**
- **Depends on:** nothing (item 526 landed). **Family:** `POST-8` 3b-ii-t1; attempt 1.

**The sites as they are today.** Re-read 10-18; no line drift since the census at `20260927T112115Z_POST-8-step1-n2.log.gz:1973–2013`.
- **(α) `test_lossy_sphere_degree2.py:176–177`.** Samples `fields.e_real` / `e_imag`, the DG(degree) splits of the N1curl phasor (`time_harmonic.py:745–754`), through the legacy path. Its `field_error` is gated at 1 % downstream, in `test_degree2_energy_mechanism.py:321–326` (`_run_at_degree`, imported at `:107–112`).
- **(α) `test_lossy_sphere_fullwave.py:218–219`.** The same field, inside `_solve`. `test_geometry_floor_discriminator.py:89` imports `_solve` and gates on it.
- **(β) `test_time_harmonic_smoke.py:206`.** Samples a CG1 interpolant of `fields.e_imag` (`:190–193`). The samples are print-only (`:264–265`); only `valid_mask` is asserted (`:208`).

**Code.**
- In the two α sites, add `owner="lowest_original_cell"` to all four calls.
- In the smoke test:
  - delete `:190–193`;
  - sample `fields.e_imag` with the owner rule instead, as `POST-4` step 3 did at `examples/mri/01_coil_phantom_fields.py:478–487`;
  - print the five `|E_imag|` values at ≤ 6 significant digits.
- Add `tests/post/test_post8_step3b_ii_th_sites.py` (`complex_only`).
- No other edit: not `evaluation.py`, not `src/`, not an example, not a t2 module.

**New module.** Imported, never restated:

| Name | Source |
|---|---|
| `_run_at_degree` | `test_lossy_sphere_degree2` |
| `_interior_probe_points` (25 points) | `test_lossy_sphere_fullwave` |
| `_original`, `_bitwise_equal_across_ranks`, `CLOSED_FORM_TOL` | `tests/post/test_point_owner_rule.py:55`, `:84`, `:46` |
| `_dg0_index` | `tests/post/test_post8_step3b_point_sar_owner.py:76` |
| `WIDTH_RTOL` (imported inside the (W) test only, as 3b-i does at `:228`) | `tests.validation.test_post8_step3a_averaged_row1` |

The fixture is `{d: _run_at_degree(d) for d in (1, 2)}`, with `msh` taken from `row["fields"].e_real.function_space.mesh`. `n_incident` per point is `compute_colliding_cells` over owned cells, then `Allreduce` SUM.

**Anchors (asserted).**
- **(X)** Every pre-existing assertion passes at its existing band in W1, untouched. That covers the three edited modules plus `test_degree2_energy_mechanism.py` and `test_geometry_floor_discriminator.py`.
- **(S)** Per degree, on the same solved field:
  - at every probe with `n_incident == 1`, the owner-rule sample is bitwise equal to the `owner=None` sample on every rank (`_bitwise_equal_across_ranks`);
  - at any probe with `n_incident > 1`, it equals `.eval` in the owner cell (owning rank, then `Allreduce` SUM) to ≤ `CLOSED_FORM_TOL`.
  
  *Predicted, printed:* 25/25 probes are single-incidence.
- **(O)** The owner original indices, at all probes and at the vertex (V), equal the `-n 2` npz exactly at `-n 4` and `-n 1`.
- **(W)** At `-n 2`, rank 0 writes `/workspace/logs/post8/step3b-ii-t1-n2.npz` (gitignored), then `Barrier`. At `-n 4` / `-n 1` the module loads it, with `pytest.fail` if it is missing. Each sampled E satisfies ‖E_w − E_2‖/‖E_2‖ ≤ `WIDTH_RTOL` (1e-9). *Predicted, printed:* ≤ 5e-10, from 3a's complex MUMPS cross-width floor of 2.21e-10.
- **(V)** The lexicographically smallest vertex of probe 0's owner cell at degree 1, with its coordinates shared by `Allreduce` SUM.

**Negative controls.**
- **(N0), asserted (rule (e)):** at V, `n_incident ≥ 4` and `max/min |E|` over the incident cells > 1 + 1e-6. *Predicted, printed:* > 1.01.
- **(N1), predicted, printed:** the legacy cell at V (`_dg0_index`, `owner=None`) differs from the owner at `-n 4` or `-n 1`. There is no measurement on this mesh.
- **(P-smoke), printed only:** the five smoke `|E_imag|` values at `-n 2` against `-n 1`, as a `.3e` relative difference. *Predicted* ≤ 1e-6.

**Windows.** Each runs in the foreground with Bash `timeout: 660000` and `-s`, and each harness call is the whole bare Bash command. Run W1 → W2 → W3, with `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reading 0 before each.
- **W1** (`-n 2`, writes the npz):
  `scripts/testing/run_and_log.sh POST-8-step3b-ii-t1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_th_sites.py tests/validation/test_lossy_sphere_fullwave.py tests/validation/test_lossy_sphere_degree2.py tests/solver/test_time_harmonic_smoke.py tests/validation/test_degree2_energy_mechanism.py tests/validation/test_geometry_floor_discriminator.py -v -s --tb=short'"`
- **W2** (`-n 4`):
  `scripts/testing/run_and_log.sh POST-8-step3b-ii-t1-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 4 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_th_sites.py -v -s --tb=short'"`
- **W3** (`-n 1`): W2's command with these changes:
  - `-n 1` and `timeout -k 30 300`;
  - label `POST-8-step3b-ii-t1-W3`;
  - `tests/solver/test_time_harmonic_smoke.py::test_time_harmonic_smoke_solve_conserves_real_power` appended after the new module, by node ID, never `-k`.

**Traps already paid for.**
- **Rule (i):** an edit after W1 starts means W1–W3 again.
- Use `-s` and container `timeout -k 30`. A kill is not a result: run `pgrep -c python3`, then the rubric cache sweep as `POST-8-step3b-ii-t1-cache-sweep`.
- **Complex build:**
  - `tests/environment/test_complex_mode.py` goes first on the path;
  - the DG0 index field is complex, so compare `.real` and assert `.imag == 0`.
- The rule is collective: pass one identical point list on every rank.
- Use the new npz name. Never touch `step3a-n2.npz` or `step3b-n2.npz`. Never flip `evaluate_vector_field_parallel`'s default.
- **The bash hook** denies non-harness commands containing "pytest", and `pgrep -f "python3 -m pytest"`.
- **Leak guard:** if `check_private_leak.py` fires, move the log unchanged to `logs/leak-held/`, never edit or bypass it, and report (item 526 precedent).
- No AED number; never open `docs/private/` or `aed_results/`.

**Scope.** Four call sites in three test modules. **Not claimed:**
- t2 (`test_coil_phantom_bfield_metrics.py:134, :149, :204, :239`; `test_straight_wire.py:289`);
- the 12 example sites;
- 3b-iii;
- `POST-8` ✅.

**Negative result.** Report and stop. Park on `attempt/POST-8-3b-ii-t1-*`, mark item 536 BLOCKED in the same commit, and never loosen 1e-9 or any band.
- **(X) red:** stash the edits and re-run the red module as `…-baseline`. A red baseline predates this item: open a known-issues entry.
- **(S) red:** a rule-path defect. New known-issues entry.
- **(O) red:** the indices are not partition-free. New known-issues entry.
- **(W) red with (O) green:** the solve floor is above 1e-9. Record the readings and annotate known-issues 2026-09-27; the review re-poses the bound.
- **(N0) red:** the item parks for a re-pose.
- **Status 124:** not a result.

There is no `PREDICTED RED:` line.

**Status it can move.** `POST-8` stays 🟡.
- The §7 state line gains "Step 3b-ii-t1 done <date> (§9 item 536)", with the single-incidence count, the max width-spread and the W1 log.
- "3b-ii (17 candidates)" becomes "t1 done; t2 (5 magnetostatic sites) and 12 examples open".
- Known-issues 2026-09-27's "Verified at" gains t1 and stays OPEN.

**Spawn line:** "Use the implementer agent to execute `POST-8` §9 item 536 (step 3b-ii-t1) per `docs/planning/chunks/POST-8.md` and the 2026-10-18 10:30 ruling beneath it.
- The harness call is the **whole** bare Bash command.
- Run W1, W2, W3 in order, in the foreground, with `timeout: 660000`, `-s` and container `timeout -k 30`. `pgrep -c python3` must read 0 before each.
- Edit only the four calls in `test_lossy_sphere_degree2.py`, `test_lossy_sphere_fullwave.py` and `test_time_harmonic_smoke.py`, and add `tests/post/test_post8_step3b_ii_th_sites.py`.
- Do not spawn agents, and never background a window.
- If the leak guard fires, move the log to `logs/leak-held/` and report.
- On a red, park on `attempt/POST-8-3b-ii-t1-*` and mark the item BLOCKED."

**Not drafted, for a later review: 3b-ii-t2 (magnetostatic, real build).**
- **Pass the owner rule** at `test_coil_phantom_bfield_metrics.py:134`, `:204` (`b_dg0`) and `test_straight_wire.py:289` (DG1 `b_field`).
- **The drafter recommends no edit at `:149` / `:239`.** These are CG1 "retired path, never gated" prints that the docstring keeps as its 3.03× partition-dependence record.
- **The width anchor would be exact owner-index equality,** plus a spread printed against the 7.28e-8 floor.
- **Records at risk if a sample sits on a facet:**
  - `NPOINTS_CONTROL_BY_VERSION`, at its 1e-4 band (`test_straight_wire.py:150–175`);
  - the coil-phantom record 0.311170 / 0.311166.
  
  The centerline's 0.341 % cross-width scatter suggests multi-incidence axis points, so t2 measures `n_incident` first.
- **Cost:**
  - straight wire: 371 s at `-n 2`;
  - coil-phantom: 9 s at `-n 4` and 20 s at `-n 1` (`20260808T170334Z_MAG-6-step3-gate-n4.log.gz`, `…170126Z_…redbaseline-n1.log.gz`).

### Review ruling 2026-10-18 10:30 (item 536)

- **The split is accepted.** t1 is queued. t2 waits for a review to rule on `:149` / `:239` and on its width anchor. The re-classification of the three γ sites as discontinuous is recorded here, not yet ruled.
- **Spot-checked by the review (re-read):**
  - the four α calls at `test_lossy_sphere_degree2.py:176–177` and `test_lossy_sphere_fullwave.py:218–219`;
  - the smoke call at `test_time_harmonic_smoke.py:206`;
  - `_dg0_index` at `test_post8_step3b_point_sar_owner.py:76`;
  - `_original` / `_bitwise_equal_across_ranks` / `CLOSED_FORM_TOL` at `test_point_owner_rule.py:55`, `:84`, `:46`;
  - `WIDTH_RTOL`, imported from `test_post8_step3a_averaged_row1` at `:228`.
- **The smoke prints are cut to ≤ 6 significant digits.** The drafter proposed `.17e`, which raises the leak-guard odds for no gated purpose.
- **(N0) stays asserted,** on the drafter's reasoning (the same code path as 3b-i's vertex control, max/min 1.6557 on the cube). The mesh differs, so a red parks the item for a re-pose; it never loosens.

### §9 item 537 — 2026-10-18 (`POST-8` step 3b-ii-t2; lane 1)

*Drafted for the 2026-10-18 18:00 daily review; the review rules on it (ruling below).*

**Licence.** The 10-17 weekly's 3b classification, the 526 brief's 3b-ii split (`POST-8.md:229–233`), item 536's re-classification of the three γ sites as discontinuous (`:378`), and the 18:00 review's t2 ruling, which binds this brief:
- the owner rule at the three DG sites;
- no edit at `:149` / `:239`;
- exact owner-index equality across widths, plus (S), as the width anchor, with the field spread printed only;
- records imported, never restated, and never re-recorded in the slot.

No band, tolerance, Done-when or target is set or moved. Template: item 536.

**`POST-8` step 3b-ii-t2: the three magnetostatic DG test sites sample under the owner rule. The rule is asserted to be a bitwise no-op at single-incidence probes, owner indices are asserted exact across `-n 1 / 2 / 4` on the coil-phantom fixture, and both edited modules re-run green at their existing bands.**
- **Executor:** implementer.
- **Tier** (real build throughout; magnetostatics needs no complex-mode source):
  - W0 smoke (collect only, `-n 1`, 30 s).
  - W1 standard (`-n 2`, 180 s), predicted ≈ 40 s. It runs two coil-phantom solves (the new module's fixture, plus the module's own node), each ≈ 14 s plus ≈ 6 s of mesh. Basis: `test_coil_phantom_bfield_metrics` took 14.30 s per call at `-n 2` (`20260821T033534Z_OPS-17-step3j-sar-padding-group.log.gz:95`), and the 0.010 m mesh takes 6.4 s (module comment `:60–63`).
  - W2 standard (`-n 4`, 180 s), predicted ≈ 30 s. Basis: 10 s at `-n 4` on the 0.010 m rung (`MAG-6-step3-gatefinal-n4`, `test-results.md:367`).
  - W3 heavy by ceiling (`-n 1`, 400 s), predicted ≈ 160 s, with **one** solve: the new module only, because its fixture executes the original test body (see "New module").
    - Basis: **143 s** at `-n 1` on the current 0.010 m rung (`20260808T170549Z_MAG-6-step3-gatefinal-n1.log.gz`, `test-results.md:369`). That log was pruned from `main` by `505dd2bd` and is readable with `git show 505dd2bd^:docs/testing/logs/…`.
    - **The 20 s figure the ruling quotes is the 0.015 m red baseline** (`20260808T170126Z_MAG-6-step3-redbaseline-n1.log.gz`, Status 1). It is not this fixture's cost.
    - The 9 s `-n 4` figure (`…170334Z_MAG-6-step3-gate-n4.log.gz`) is also a Status 1 run. W2's basis is therefore the green `gatefinal-n4`.
  - W4 heavy (`-n 2`, container `timeout -k 30 560`), predicted ≈ 430 s.
    - The straight-wire module ran 7 passed in 369.95 s, harness 371 s (`20260828T050256Z_MAG-20-step1-anchor-module.log.gz:3641`, `test-results.md:1107`).
    - The new module's straight-wire test adds one h = 0.0025 solve and one h = 0.004 solve. Splitting the 371 s by cells across the module's 11 solves gives an estimate of ≈ 38 s + ≈ 10 s; this is an estimate, not a measurement.
    - Worst case 590 s, inside the Bash tool's 660 s.
- **Scope of the windows (ruling).** The width windows (`-n 4`, `-n 1`) run the coil-phantom sites only. The straight-wire site gets (X) + (S) + (R) at `-n 2` only, in the one heavy window W4.
- **Class:** lane 1. **Cost:** ≈ 6.1 measured / ≈ 40 predicted slot-minutes. Windows ≈ 11 min, code ≈ 15, `pgrep` checks, commit and journal ≈ 10.
- **Park rules:** if W1 has not started by minute 25, park. If W4 has not started by minute 40, park (rule (d)).
- **Depends on:** nothing in §9. Item 536 landed (`b71e7c8a`), and this module imports from its test file.
- **Family:** `POST-8` 3b-ii-t2; attempt 1.

**The sites as they are today** (re-read 10-18 18:00):
- **`test_coil_phantom_bfield_metrics.py:134–138`** samples `b_dg0` on the 9 centerline points (`:120–123`; x = y = 0, on the axis). This is the DG0 interpolant of `compute_b_field()` (`:116–118`, name `"B_dg0"`), through the legacy path. Gated at the 0.60 jump-ratio bound (`:175`).
- **`:204`** samples `b_dg0` on the 18 ±x symmetry points (`:195–203`; the middle row is z = 0). Gated at `PHANTOM_SYMMETRY_REL_TOL` 0.35 (`:229`).
  - The record, docstring `:30–31` / message `:234`, reads 0.311170 / 0.311166 at `-n 2 / 4` (and 0.311226 at `-n 1`, comment `:101`).
  - **It dates from `MAG-6` step 5 (2026-08-09), before the 0.11 image.** No log on `main` prints it today (`zgrep` for `RANKSPREAD_INPUT max_rel_diff` over `docs/testing/logs/` finds none).
- **`:149`, `:239`: not edited (ruling).** These are the CG1 `b_lagrange` (name `"B"`, `:111`) print-only paths, the 3.03× partition record of the docstring `:12–22`.
- **`test_straight_wire.py:289`**, inside `_sample_radial` (`:282–295`), samples `compute_b_field()` (DG1 of curl A) along +x at y = z = 0. Its consumers, and the gates at risk:
  - `test_straight_wire_b_field` (`:307`, reported only);
  - `test_analytic_bc_improves_on_natural_bc` (`:379`, the 0.85 ratio at h = 0.004);
  - `test_straight_wire_convergence` (`:418`, `RATE_MIN < rate < RATE_MAX` at n = 8; the n = 10 margin is 0.0246, per `:450–452`);
  - `test_domain_l2_record` (`:515`): `NPOINTS_CONTROL_BY_VERSION` (`:150–153`) at `NPOINTS_CONTROL_BAND` 1e-4 (`:175`), read through `_npoints_control()` (`:156–172`).

**Code.**
- Add `owner="lowest_original_cell"` to the calls at `test_coil_phantom_bfield_metrics.py:134–138` and `:204`, and at `test_straight_wire.py:289`. That is three calls; nothing else in either module changes.
- Add `tests/post/test_post8_step3b_ii_mag_sites.py` (real build; no `complex_only`).
- No other edit: not `:149` / `:239`, not `evaluation.py`, not `src/`, not an example, not a t1 module. Never flip `evaluate_vector_field_parallel`'s default.

**New module.** Imported, never restated:

| Name | Source |
|---|---|
| `evaluate_vector_field_parallel`, `point_owner_cells` | `fem_em_solver.post.evaluation` (`:224`, `:133`) |
| `_n_incident`, `_owned_incident_cells`, `N0_MIN_INCIDENT`, `N0_SPREAD_FLOOR` | `tests/post/test_post8_step3b_ii_th_sites.py:61`, `:77`, `:46`, `:47` |
| `CLOSED_FORM_TOL`, `_bitwise_equal_across_ranks` | `tests/post/test_point_owner_rule.py:46`, `:84` |
| `_dg0_index` | `tests/post/test_post8_step3b_point_sar_owner.py:76` |
| the coil-phantom module, as `cpm` | `import tests.validation.test_coil_phantom_bfield_metrics as cpm` |
| the straight-wire module, as `sw` (`_solve_straight_wire` `:248`, `_sample_radial` `:282`, `_npoints_control` `:156`, `RESOLUTION` `:65`, `R_MIN` / `R_MAX` / `R_MAX_BC` `:70–77`, `E_OMEGA_LADDER` `:101`) | `import tests.validation.test_straight_wire as sw` |

- **Import the two modules, never their test names.** Pytest collects any `test_*` function or `Test*` class in the new module's namespace, which would double every solve.
- If W0 shows that the t1 import fails in the real build, copy `_n_incident` / `_owned_incident_cells` verbatim with a source-line comment instead. Do this before W1; it is not a rule (i) event.

**Fixtures (the spy).**
- `cp` (module scope). It wraps `evaluate_vector_field_parallel` in a spy that records `(function, points copy, owner kwarg)` and forwards every positional and keyword argument unchanged.
  - Apply it with `pytest.MonkeyPatch()`, `setattr(cpm, "evaluate_vector_field_parallel", spy)`, inside `try` / `finally: undo()`.
  - Then call `cpm.test_coil_phantom_bfield_metrics_are_finite_smooth_and_symmetric()`.
  - That one call is the coil-phantom solve, its own assertions included, so W3 needs no second solve.
- `test_sw_s_r` builds its own spy the same way on `sw`. It calls:
  - `sw._solve_straight_wire(sw.RESOLUTION, comm)`, then `_sample_radial` at each n of `_npoints_control()[1]` over `(R_MIN, R_MAX_BC)` (n = 10 is also `test_straight_wire_b_field`'s set), and at n = 8 over `(R_MIN, R_MAX)`;
  - `sw._solve_straight_wire(sw.E_OMEGA_LADDER[0], comm)`, the 0.004 rung that `:373` and `:412` use as literals, then n = 10 over `(R_MIN, R_MAX_BC)` and n = 8 over `(R_MIN, R_MAX)`.
  - The natural-BC field of `:378` shares that mesh and those points, so its incidence is the same. It is not solved again.

**Anchors (asserted).**
- **(X)** Every pre-existing assertion passes at its existing band, untouched:
  - the coil-phantom module in W1 / W2 (its own node) and inside `cp` at every width;
  - the straight-wire module in W4. `test_domain_l2_record` is the 1e-4 `NPOINTS_CONTROL_BY_VERSION` gate.
- **(E) The edit is the one ruled.**
  - `cp`'s spy sees exactly four calls, in source order, with names `"B_dg0"`, `"B"`, `"B_dg0"`, `"B"` and owners `RULE`, `None`, `RULE`, `None`. So `:149` / `:239` are untouched and on the legacy path.
  - Every `sw` spy call carries `owner == RULE`.
- **(S)** For each captured DG point set (centerline 9, symmetry 18; straight wire 46 points at h = 0.0025 and 18 at h = 0.004), on the captured field:
  - `n_incident = _n_incident(msh, pts, comm)` is measured **first** and printed per probe;
  - the owner-rule sample is `_bitwise_equal_across_ranks`;
  - at `n_incident == 1`, owner bytes equal `owner=None` bytes;
  - at `n_incident > 1`, the owner sample equals `.eval` in the owner cell to ≤ `CLOSED_FORM_TOL` relative. The owning rank evaluates from `point_owner_cells`' `owner_local`, and the result is shared by buffer `Allreduce` SUM. This is a ≈ 10-line real-field analogue of t1's `_e_in_cell`, which is complex-specific.
  - Print the single-incidence count per set, and the legacy cell (`_dg0_index`, `owner=None`) beside the owner at each multi-incidence probe.
- **(R) The records do not move.**
  - At the record widths (coil-phantom `-n 2` and `-n 4`; straight wire `-n 2`), every probe's owner sample is bytewise equal to its legacy sample on the same solve. At `-n 1`, coil-phantom (R) is printed only.
  - This is an identity on the same solve, not a band, and it is the only valid test of "the record moved": the docstring digits predate the 0.11 image (see "The sites"). The edited module's printed metric is set beside the docstring record for information only.
  - (R) can be red only at a multi-incidence probe whose owner cell differs from the legacy cell at that width.
- **(O)** Coil-phantom only.
  - At `-n 2`, rank 0 writes `/workspace/logs/post8/step3b-ii-t2-n2.npz` (gitignored) with the 27 + 1 owner original indices (the probes plus V), `n_incident`, and the 27 owner-path `b_dg0` samples, then `Barrier`.
  - At `-n 4` / `-n 1` the module loads it, with `pytest.fail` if it is missing, never a recompute. The owner indices must equal the npz **exactly**.
  - **Printed only (ruling):** max over probes of ‖B_w − B_2‖/‖B_2‖ (`.3e`), beside the `MAG-18` LU cross-width floor 7.28e-8 (`test_straight_wire.py:113–122`, a different fixture; quoted from the comment, not imported). *Predicted:* ≤ 1e-6. **No band is set on it.**
- **(V)** The lexicographically smallest vertex of centerline probe 0's owner cell (the pattern of t1's `_vertex_v`, on `b_dg0`'s mesh), with its coordinates shared by `Allreduce` SUM.

**Negative controls.**
- **(N0), asserted (rule (e)):** at V, `n_incident ≥ N0_MIN_INCIDENT` (4), and `max/min |B_dg0|` over V's incident cells > 1 + `N0_SPREAD_FLOOR`. The values come from each rank's `_owned_incident_cells` evaluated in-cell, then `allgather`.
  - This shows a cell choice at a shared vertex moves the DG0 sample, so (O) has teeth.
  - *Predicted, printed:* > 1.01. Degree-1 curl A is cell-wise constant, and the module docstring's 3.03× CG1 swing is that jump averaged.
- **(N1), printed:** the legacy cell at V against the owner, at every width.

**Prediction for (S) / (R), printed, not asserted.** The centerline jump ratio scatters 0.341 % across widths, and the mirror metric 0.022 % (`:99–102`). Both are ≥ 10⁴ × the LU floor. The DG0 values themselves are cell-local, so that scatter points at **multi-incidence centerline and/or z = 0 probes**. A red (R) is therefore a live outcome; see "Negative result".

**Windows.** Each runs in the foreground with Bash `timeout: 660000` and `-s`, and each harness call is the whole bare Bash command. Before each window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Run W0 → W1 → W2 → W3 → W4.
- **W0** (collect, real build):
  `scripts/testing/run_and_log.sh POST-8-step3b-ii-t2-W0-collect "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 30 python3 -m pytest tests/post/test_post8_step3b_ii_mag_sites.py --collect-only -q -s'"`
- **W1** (`-n 2`, writes the npz):
  `scripts/testing/run_and_log.sh POST-8-step3b-ii-t2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_mag_sites.py::test_cp_s_r_n0 tests/post/test_post8_step3b_ii_mag_sites.py::test_cp_o_width tests/validation/test_coil_phantom_bfield_metrics.py -v -s --tb=short'"`
- **W2** (`-n 4`): W1's command with `-n 4`, label `POST-8-step3b-ii-t2-W2`.
- **W3** (`-n 1`):
  `scripts/testing/run_and_log.sh POST-8-step3b-ii-t2-W3 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 1 python3 -m pytest tests/post/test_post8_step3b_ii_mag_sites.py::test_cp_s_r_n0 tests/post/test_post8_step3b_ii_mag_sites.py::test_cp_o_width -v -s --tb=short'"`
- **W4** (`-n 2`, straight wire, heavy):
  `scripts/testing/run_and_log.sh POST-8-step3b-ii-t2-W4-wire "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 560 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_mag_sites.py::test_sw_s_r tests/validation/test_straight_wire.py -v -s --tb=short'"`
- The form is the `MAG-20` anchor-module window's (`20260828T050256Z_…log.gz`, header `Command:`) plus `/workspace` on `PYTHONPATH` for the `tests.*` imports (3b-i's W4-real form, `:303`). There is no `dolfinx-complex-mode`, no `FEM_EM_REQUIRE_COMPLEX`, and no `tests/environment/test_complex_mode.py`.
- Node IDs select tests; never `-k`.

**Traps already paid for.**
- **Rule (i):** an edit after W1 starts means W1–W4 again.
- Use `-s` and container `timeout -k 30`. A kill is not a result: run `pgrep -c python3` first, then at 0 run the rubric cache sweep as `POST-8-step3b-ii-t2-cache-sweep`. Never lengthen a timeout.
- **W4 sits under the ≈ 10-minute line** (ceiling 590 s with `-k`). It needs no durable capture, but if it is killed, check for orphaned ranks before anything else.
- **The bash hook** denies non-harness commands containing "pytest", and `pgrep -f "python3 -m pytest"`. `pgrep -c python3` is the check that runs.
- **Spy hygiene:**
  - a module-scoped fixture cannot take the function-scoped `monkeypatch`, so use `pytest.MonkeyPatch()` and `undo()` in `finally`;
  - patch the attribute on `cpm` / `sw`, which is what their call sites resolve, never on `fem_em_solver.post.evaluation`;
  - copy the points, because `_sample_radial` builds them locally.
- The rule is collective: every point set is a deterministic array, identical on every rank.
- `assemble_scalar`, `n_incident` and `.eval` are rank-local; reduce before asserting.
- The npz name is new. Never touch `step3a-n2.npz`, `step3b-n2.npz` or `step3b-ii-t1-*.npz`.
- **Leak guard:** if `check_private_leak.py` fires, move the log unchanged to `logs/leak-held/`, never edit or bypass it, and report (item 526 precedent). Print field values at ≤ 6 significant digits; ratios and spreads at `.3e`.
- No AED number; never open `docs/private/` or `aed_results/`.

**Scope.** Three call sites in two test modules, plus the new module. **Not claimed:**
- `:149` / `:239` (kept by ruling);
- the 12 example sites (`example-runner` work);
- 3b-iii;
- `POST-8` ✅;
- any band, and any width bound on magnetostatic field values;
- straight-wire width invariance (not run at `-n 4` / `-n 1`, on cost).

**Negative result.** Report and stop. Park on `attempt/POST-8-3b-ii-t2-*`, mark item 537 BLOCKED in the same commit (rule (d)), and never loosen a band or re-record in the slot.
- **(X) red:** stash the three edits and re-run only the red module, with the same command shape and label `…-baseline`. A red baseline predates this item: the straight-wire module last ran green on 08-28, and the coil-phantom module on 08-21 (complex build). Open a known-issues entry.
- **(E) red:** the edit is not the one ruled. Park.
- **(S) red:** a rule-path defect. New known-issues entry.
- **(R) red (the ruling's reported negative):** a record-bearing probe is multi-incidence, and the rule picks a different cell than the legacy path. Journal, per probe, `n_incident`, the owner and legacy original indices, the sample delta, and the edited module's printed metric beside its record. Never re-record. The review re-poses the site, for example a `(1*)` re-record by `record-reconciler`, or a source-field switch.
- **(O) red:** the indices are not partition-free. New known-issues entry.
- **(N0) red:** the item parks for a re-pose.
- **Status 124:** not a result; park with the measured time.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** `POST-8` stays 🟡. Its Done-when (`PROJECT_PLAN.md:2670`) still needs the 12 examples and 3b-iii.
- On green, the §7 state line's "3b-ii t1 done; t2 (5 magnetostatic sites) and 12 examples open" becomes "3b-ii t1 and t2 done; 12 examples open".
- It gains: "**Step 3b-ii-t2 done** <date> (§9 item 537): the three magnetostatic DG sites (`test_coil_phantom_bfield_metrics.py` centerline / symmetry, `test_straight_wire.py::_sample_radial`) use the owner rule; `:149` / `:239` kept as the CG1 3.03× record (ruling); <k>/<N> probes single-incidence, the rule bytewise = legacy at every record-width probe; coil-phantom owner indices exact at `-n 1 / 2 / 4`, field spread <max> printed beside the 7.28e-8 LU floor; both modules green (<W1 log>, <W4 log>)."
- Known-issues 2026-09-27's "Verified at" gains t2, and the entry stays OPEN until 3b-iii.

**Execution order.** Code → W0 → W1 → W2 → W3 → W4 → `git status --porcelain examples/` (must print nothing) → commit (the three edits, the new module, harness logs, the §7 sentence, the known-issues line) → journal.

**Spawn line:** "Use the implementer agent to execute `POST-8` §9 item 537 (step 3b-ii-t2) per `docs/planning/chunks/POST-8.md` and the 2026-10-18 18:00 ruling beneath it.
- The harness call is the **whole** bare Bash command.
- Run W0, W1, W2, W3, W4 in order, in the foreground, with `timeout: 660000`, `-s` and container `timeout -k 30`. `pgrep -c python3` must read 0 before each.
- Edit only the three calls at `test_coil_phantom_bfield_metrics.py:134` / `:204` and `test_straight_wire.py:289`, and add `tests/post/test_post8_step3b_ii_mag_sites.py`. Make no edit at `:149` / `:239`.
- Do not spawn agents, and never background a window.
- If the leak guard fires, move the log to `logs/leak-held/` and report.
- On a red, park on `attempt/POST-8-3b-ii-t2-*` and mark the item BLOCKED."

### Review ruling 2026-10-18 18:00 (item 537)

- **The t2 design is ruled.** No edit at `test_coil_phantom_bfield_metrics.py:149` / `:239` (the CG1 prints stay as the 3.03× partition record). The width anchor is **exact owner original-index equality at `-n 1 / 2 / 4`** on the coil-phantom sites, plus (S), owner = legacy bitwise at single-incidence probes. Both are identities. The field-value width spread is **printed only**, beside the 7.28e-8 LU floor. No band is set or moved.
- **The drafter's cost correction is accepted.** The cited 20 s / 9 s runs were Status 1 baselines; the current fixture's `-n 1` window is 143 s (`test-results.md:369`). W3 is heavy by ceiling and runs one solve.
- **The module-scoped `pytest.MonkeyPatch` spy on `cpm.evaluate_vector_field_parallel` is accepted.** It captures the field and points from the unedited test body. The fallback (extracting a helper inside the coil-phantom module) is **not** licensed; if the spy cannot capture, park.
- **(R) is printed only.** An owner-vs-legacy byte check at a multi-incidence probe differs by construction, so it is not a gate. The records are gated by (X): the edited modules' own assertions at their existing bands, on the owner rule. Print per probe `n_incident`, the owner and legacy indices and the sample delta, and the metric beside the 0.311170 / 0.311166 docstring record. An (X) red on a record is the reported negative of the brief's `:655` disposition: journal, park, never re-record.
- **W0 (collect-only) is kept** to check that the t1 helpers import in the real build. If they do not, copy the two helpers verbatim with a source citation; nothing else.
- Spot-checked by the review: `test_coil_phantom_bfield_metrics.py:30`, `:234`, `:273` (the record text); `tests/validation/test_straight_wire.py:150–175`.

**Outcome 2026-10-19 (21:00 CDT slot, attempt 1): done, after one rule-(i) re-run.** Three calls on the owner rule; new `tests/post/test_post8_step3b_ii_mag_sites.py`; no edit at `:149` / `:239`. First pass: W0–W3 green, W4 red on a check **outside the anchors** that the module itself added — re-sampling the four h = 0.0025 batches as one 46-point batch differed from the spy-captured values in the last ULP at some probes (`.eval` bytes depend on the batch; the anchors (S) / (E) / (X) were green, all 7 straight-wire tests passed) — `20261019T020649Z_POST-8-step3b-ii-t2-W4-wire.log`, Status 1. The check was changed to compare call by call on the same batch, with that measurement in a code comment, and W1–W4 re-ran (rule (i)). Re-run: (E) four `cpm` calls `B_dg0` / `B` / `B_dg0` / `B`, owners rule / None / rule / None; six `sw` calls, all on the rule. (S) single-incidence 9/9 centerline and 18/18 symmetry at `-n 1 / 2 / 4`, 46/46 (h = 0.0025) and 18/18 (h = 0.004) straight wire at `-n 2`; rule bitwise = legacy everywhere; no multi-incidence probe, so (R) has no delta to print. (O) owner indices (27 probes + V) exact at `-n 4` / `-n 1`; field spread (printed) 5.614e-05 / 9.647e-04 beside the 7.28e-8 LU floor — the predicted ≤ 1e-6 missed (a run-to-run spread of the solve itself: the same-width `cpm` mirror metric reads 0.311031 in the fixture and 0.311034 in the module's own node). (N0) V has 26 incident cells, max/min `|B_dg0|` 4.225523 / 4.225331 / 4.229599 at `-n 2 / 4 / 1`. (N1, printed) legacy cell at V 3762 / 3762 / 4073 vs owner 3013. Mirror metric owner path 0.311031 / 0.311037 / 0.310984 at `-n 2 / 4 / 1` beside the docstring record 0.311170 / 0.311166 (record not re-recorded). (X) `cpm` 3 passed at `-n 2` / `-n 4` and the fixture at `-n 1`; straight wire 8 passed (rate 0.7900, E_Omega 1.0617170107e-01, n_points rows ≤ 3.390e-06 relative of their records, inside the 1e-4 band). No leak-guard hit. Logs `20261019T020246Z_POST-8-step3b-ii-t2-W0-collect.log` (4 s), `20261019T021325Z_POST-8-step3b-ii-t2-W1.log` (24 s), `20261019T021402Z_POST-8-step3b-ii-t2-W2.log` (22 s), `20261019T021435Z_POST-8-step3b-ii-t2-W3.log` (152 s), `20261019T021728Z_POST-8-step3b-ii-t2-W4-wire.log` (356 s); first-pass `…020259Z_…W1.log`, `…020333Z_…W2.log`, `…020407Z_…W3.log` green.

### §9 item 547 — 2026-10-19 (lane 1, POST-8 3b-ii-e1)

**Executor: `example-runner`.** *Drafted read-only for the 2026-10-19 03:00 daily review; the review rules on it. It sets and moves no band, tolerance, Done-when or target.*

**Licence.**
- The 10-17 weekly's 3b classification (`PROJECT_PLAN.md` ≈ 4753–4758): every point-sampling site of the census adopts the owner-cell rule.
- Item 526's 3b-ii split: "the examples are `example-runner` work, with guide re-sync" (`POST-8.md:233`).
- The 2026-10-18 18:00 ruling on item 537, applied here to the example half:
  - the width anchor is exact owner original-index equality at `-n 1 / 2 / 4`, plus (S);
  - the field-value spread is printed only, and so is (R);
  - records are gated by (X), the edited files' own assertions at their existing bands, and are **never re-recorded in the slot**.
- Templates: items 536 / 537 (the spy, the (S) / (O) / (N0) forms) and item 544 (example-runner windows and censuses).

**The twelve sites, re-read 2026-10-19.** The field type was verified at each call, as item 526 asked 3b-ii's drafter to do (`POST-8.md:233`). This table supersedes 526's α / β / γ labels.

| Tag | Call site (today) | Sampled field | Family | What consumes it | Record width | Item |
|---|---|---|---|---|---|---|
| `mag:1` | `01_straight_wire.py:273` | `B_lag` (`:250–252`), the CG1 interpolant of `compute_b_field()` | CG1 of a DG field (526 β) | the guide table only (no assert); `B_lag` also feeds the export / VTX at `:386`, `:425`, `:444` | 2 | 547 |
| `mag:5` | `05_gauge_cross_check.py:123` (`solve_with`) | `compute_b_field()` (`:116`) | DG | `AGREEMENT_RTOL` 0.05 at `:285` and the gauge asserts at `:207–306` | 2 | 547 |
| `th:3` | `03_dielectric_sphere_in_uniform_field.py:227–228` (`_probe`, called at `:344`, `:423`) | `fields.e_real` / `e_imag` | DG split of N1curl | `:369–470` (interior / exterior / jump / blind) | 2 | 547 |
| `th:4` | `04_evanescent_waveguide_decay.py:227–228` (`_probe_ey`) | the DG split at `:333`; the CG1 export `e_re` / `e_im` at `:375`, `:380` | DG + CG1 | `GAMMA_RTOL`, `FIELD_L2_MAX`, `IMAG_RATIO_MAX`, `CG1_VS_NEDELEC_MAX`, `PROFILE_RMS_MAX` (`:356–479`) | 2 | 547 |
| `th:6` | `06_larmor_lossy_sphere.py:157–158` (`_interior_errors`, `:234`) | DG split | DG | `_check_record` at the 1 % `REPRODUCTION_BAND` (`:202`), `:253–426` | 2 | 547 |
| `th:7` | `07_element_order_lossy_sphere.py:170–171` | DG split, degree 1 and 2 | DG | `:259–425`, including exact cell and DOF counts | 2 | 547 |
| `th:9` | `09_pec_sphere_as_hole.py:161–162` (`_solve`, called at `:270`, `:348`) | DG split | DG | `BAND` and the controls (`:319–379`) | 2 | 547 |
| `mag:2` | `02_circular_loop.py:269` (no `comm` argument) | `compute_b_field()` (`:255`) | DG | the guide table only (no physics assert) | 2 | 548 |
| `mag:4` | `04_helmholtz_analytic_comparison.py:200` (`run_case`) | `compute_b_field()` (`:147`) | DG | the guide table only (0 asserts) | 2 | 548 |
| `mag:6` | `06_h_convergence_rate.py:212–214` | `b_num_cg` (`:192–194`), a CG1 interpolant **on purpose**: "the exported field is the asserted field" (`:205–210`) | CG1 | `export_error < errors[0]` (`:272`), monotone decay (`:285`) | 2 | 548 |
| `ports:1` | `01_two_torus_port_pair.py:303` (`_path_voltage`, orders 2049 / 4097) | `fields.e_complex` (`:545`) | raw N1curl | `PATH_QUADRATURE_TOLERANCE` (`:585`); 1e-6 reproduction of four 0.11-image records (`:739`); mutual / blind / S gates (`:750–771`) | 2 | 548 |
| `ports:12` | `12_birdcage_ring_quarter_turn.py:503–508` | `col_p17` / `col_p21` `["fields"].e_complex` | raw N1curl | an ungated 64-point ring print (`:529–543`) | 4 | 548 |

**Corrections to item 526's labels.**
- `mag:2` / `mag:4` / `mag:5` are listed by 526 as γ "source fields", but they sample `compute_b_field()`, which is discontinuous. Item 536 found the same for the test twins (`POST-8.md:378`), so the rule applies.
- `mag:6`, also γ in 526, samples a deliberate CG1 interpolant (ruling question R2 below).
- `ports:1`, `ports:12` and `th:4` (γ in 526) sample the raw N1curl phasor or its DG split. Both are discontinuous.

**Out of the census's 18, named for the weekly, not edited:** `mri:1` `01_coil_phantom_fields.py:486–487` samples the source fields `e_field` / `b_field` on the legacy path. The census row was "-", not a candidate.

**Why two items.** The split follows run cost and width.
- **547 (e1):** the seven sites whose example costs ≤ ≈ 30 s at `-n 2`, so the full `-n 1 / 2 / 4` ladder fits a slot. It also carries the three `point_sar` example re-gates owed by item 526's scope (`POST-8.md:325`).
- **548 (e2):** the five examples of ≥ 110 s each, run at record width only.

**`POST-8` step 3b-ii-e1: seven example sites sample under the owner rule. Each edited example runs as `__main__` under a spy: its own assertions hold at their imported bands, the rule is bitwise equal to legacy at single-incidence probes (call by call, same batch), and the owner original indices are exact across `-n 1 / 2 / 4`. The three `point_sar` examples re-run green.**
- **Executor:** `example-runner`.
- **Tier:** standard, two windows heavy by ceiling; ranks `-n 1 / 2 / 4`. Measured run-of-record basis (`-n 2`):

  | Tag(s) | Seconds | Source |
  |---|---|---|
  | `mag:1` | 9 | `20260826T170155Z_EX-30-root2-run-mag1.log` |
  | `mag:5` | 8 | `05_gauge_cross_check.md:62–63` (`20260810T110311Z_EX-10-run1.log`) |
  | `th:3` + `th:4` | 14 together | `20260825T050214Z_EX-30-th-run-3to4.log` |
  | `th:6` | 27 (with `th:5`: 55) | `20260813T200415Z_EX-19-example-n2.log`; `20260825T050232Z_EX-30-th-run-5to6.log` |
  | `th:7` | 16 | `20260819T140334Z_EX-25-example-n2.log` |
  | `th:9` | 4 | `20260906T020307Z_EX-50.log` |
  | `mat:3` | ≈ 65 | `03_gelled_saline_presets.md:48` |
  | `mri:2` | 14 | `20260808T020414Z_EX-3-gate.log` |
  | `mat:2` | ≈ 40, unmeasured alone (predicted from the gate's 40 s, `02_lossy_sphere_sar.md:74–76`) | — |

  Predicted windows: the seven `-n 2` spied runs ≈ 130 s in total (example time plus ≈ 5 s of pytest import each); `-n 4` groups ≈ 25 s (mag) and ≈ 80 s (th); `-n 1` groups ≈ 40 s (mag) and ≈ 160 s (th, heavy by ceiling); the three re-gates ≈ 120 s; collects, dry run and censuses ≈ 30 s. **Total ≈ 10 min of windows.**
- **Class:** lane 1, example. **Cost:** ≈ 8.4 measured (example-class median) / ≈ 50 predicted slot-minutes: fixed ≈ 5, code ≈ 15, windows ≈ 10, plus per-call overhead ≈ 5, guides ≈ 8, commit and journal ≈ 6.
- **Slot-first only** (never take-next). **Park (rule (d)):**
  - the first `-n 2` spied window not started by minute 25;
  - the `-n 1` th group not started by minute 40;
  - the re-gates not started by minute 46.
- **Depends on:** nothing in §9. Items 526 / 536 / 537 landed (`859fa8af`), and the new module imports from their test files. The queueing review rules R1 first (below).
- **Family:** `POST-8` 3b-ii-e1; attempt 1.

**Code.** Each edit adds `owner="lowest_original_cell"` to the call, the same keyword the t1 / t2 sites use (`evaluation.py:224`, `OWNER_RULES` `:14`). Nothing else changes in those files.
- `th:3` `:227–228`; `th:4` `:227–228`; `th:6` `:157–158`; `th:7` `:170–171`; `th:9` `:161–162`: two calls each.
- `mag:5` `:123`: one call.
- **`mag:1` `:273`, under ruling R1 (the drafter's recommendation):**
  - sample the DG source `B` (`:245`) instead of `B_lag`, with the rule, as `POST-4` step 3 did at `mri/01:478–487` and t1 did in the smoke test;
  - **keep `B_lag`** (`:250–252`), because the export and the VTX round trip use it;
  - rewrite the comment at `:247–248`: "DG functions … need proper cell indices for eval" is now false, so it says "for export and the VTX round trip".
  - **If R1 is declined,** add only the keyword on `B_lag` and leave the guide table alone.
- No edit to `src/`, `evaluation.py`'s default, any t1 / t2 module, `mri:1`, the `ans:2` script, or a figure.

**New module `tests/post/test_post8_step3b_ii_example_sites.py`** (licensed here; the executor's only `tests/` file).
- **Gates:**
  - module-level `pytestmark = [pytest.mark.slow, skipif(os.environ.get("FEM_EM_POST8_EXAMPLES") != "1")]`, so no CI job ever runs whole examples;
  - the complex tests are `complex_only`; the mag tests skip when `IS_COMPLEX` (`tests/complex_mode.py:23`).
  - Every window sets `FEM_EM_POST8_EXAMPLES=1` and must read **N passed, 0 skipped**.
- **One test per example:** `test_mag1`, `test_mag5`, `test_th3`, `test_th4`, `test_th6`, `test_th7`, `test_th9`. A dict maps each tag to its path and record width.
- **Imported, never restated:**

  | Names | Source |
  |---|---|
  | `point_owner_cells`, `evaluate_vector_field_parallel` (captured as `real` at import) | `fem_em_solver.post.evaluation:133`, `:224` |
  | `_n_incident`, `_owned_incident_cells`, `N0_MIN_INCIDENT`, `N0_SPREAD_FLOOR` | `tests/post/test_post8_step3b_ii_th_sites.py:61`, `:77`, `:46`, `:47` (imports in the real build: t2's W0) |
  | `CLOSED_FORM_TOL`, `_bitwise_equal_across_ranks` | `tests/post/test_point_owner_rule.py:46`, `:84` |
  | `_dg0_index` | `tests/post/test_post8_step3b_point_sar_owner.py:76` |
  | `FIELD_SCALE_FLOOR` | `tests.tolerances` |
  | `complex_only`, `IS_COMPLEX` | `tests/complex_mode.py` |

  **Do not import t2's module.** Its `_s_and_r` / `_b_in_cell` (`test_post8_step3b_ii_mag_sites.py:89–154`) cast to float64, and importing it pulls in `cpm` / `sw`. Write dtype-preserving analogues, buffer dtype `f.x.array.dtype`, with a source-line comment.
- **The spy runs the example as its run of record.** `_run_spied(tag)`:
  1. Take a `pytest.MonkeyPatch()`.
  2. `setattr` the spy on **both** `fem_em_solver.post.evaluation` and `fem_em_solver.post`. `mag:1` imports the name from the package (`01_straight_wire.py:18`); the others import from the module. Both bind when `run_path` executes the example's imports.
  3. Set `sys.argv = [str(path)]`.
  4. Call `runpy.run_path(str(path), run_name="__main__")`, catching `SystemExit`: the code must be 0 or None, otherwise re-raise.
  5. `undo()` in `finally`.

  The example's own assertions run unmodified inside step 4: that is (X).
- **What the spy records.** It records a call only when `sys._getframe(1).f_code.co_filename` resolves to the example's path. `src` callers bound the real function at their own import, and anything else is forwarded unrecorded. For a recorded call, on every rank (the example's call is collective, OPS-40/55):
  1. Forward the call unchanged; this is what the example receives.
  2. Measure `n_inc = _n_incident(msh, pts, comm)` **first**.
  3. `own = real(f, pts, comm, owner=RULE)` and `leg = real(f, pts, comm, owner=None)`, on **the identical `pts` array of this call**.
  4. `owner_orig, owner_local = point_owner_cells(msh, pts, comm)`.
  5. The legacy cell, from `_dg0_index(msh)` evaluated with `owner=None` (t2's `_legacy_cells` pattern; complex build: `.real`, assert `.imag == 0`).
  6. At each multi-incidence probe, an in-cell evaluation on the owning rank, shared by buffer `Allreduce` SUM.

  Store the caller line, `f.name`, the element family, the `owner` kwarg, and copies of everything above. `comm` is `kwargs.get("comm")`, else the first positional argument after `points`, else `MPI.COMM_WORLD`, because `mag:2` passes none.

**Anchors (asserted).**
- **(X)** The example's own assertions pass at their imported bands, untouched, in every spied window (table above), and `run_path` exits 0. Each window's log carries the example's full stdout (`-s`).
- **(E)** Every recorded call carries `owner == RULE`, and at least one call is recorded. The returned values are bytewise equal to `own` (the same call, the same batch).
  - *Predicted, printed:* call counts `mag:1` 1, `mag:5` 2, `th:3` 4, `th:4` 6, `th:6` 4, `th:7` 4, `th:9` 4.
  - The call count equals the `-n 2` npz's at the other widths.
- **(S)** Per recorded call:
  - `own` is `_bitwise_equal_across_ranks`;
  - at every probe with `n_inc == 1`, `own` bytes equal `leg` bytes, **call by call on that call's own point array**;
  - at `n_inc > 1`, `own` equals the in-cell evaluation to ≤ `CLOSED_FORM_TOL` relative, with the denominator floored at `FIELD_SCALE_FLOOR`.
- **(O)** At `-n 2`, rank 0 writes `/workspace/logs/post8/step3b-ii-e1-<tag>-n2.npz` (gitignored) and then `Barrier`. The file holds, per call, `owner_orig`, `n_inc`, the caller line and `own`. At `-n 4` / `-n 1` the test loads it, with `pytest.fail` if it is missing (never recompute), and asserts the call count, every `n_inc` array and every `owner_orig` array **exactly** equal.
  - *Printed only, no band:* the max over calls and probes of ‖own_w − own_2‖ / ‖own_2‖ (`.3e`).
  - *Predicted:* th ≤ 5e-10 (t1 measured 1.108e-11 on the same sphere family); mag ≤ 1e-3 (t2 measured 9.647e-04, a run-to-run spread of the magnetostatic solve).
- **The re-gates (rule (c) for item 526):** the unedited `mat:2`, `mat:3` and `mri:2` run bare at `-n 2`, Status 0 with their own assertions.
  - *Predicted:* every printed figure digit-identical to its guide's record. They sample cell centroids or an imposed uniform field, all single-incidence (526 Scope, `POST-8.md:325`).
  - The executor compares each printed record row with the guide by hand and journals it.

**Predictions for (S), printed, not asserted.**
- Single-incidence everywhere except:
  - **`th:4`'s two profile calls (`:380`):** the 25 points sit at z = L/2, a lattice plane of the `create_box` mesh (n = 24 in z, `:184–188`), and x = a/2 at k = 12. They are multi-incidence. They sample the **CG1** export, which is continuous, so the owner value equals the in-cell value and `PROFILE_RMS_MAX`'s printed 0.200 % does not move at its printed digits.
  - **`mag:1`'s point 0** at (r_wire, 0, 0): an OCC cylinder's seam lies at φ = 0, so a multi-incidence wire/air edge is likely.
- **Printed per multi-incidence probe:** `n_inc`, the owner and legacy indices, and the delta between the owner and legacy samples ((R), ruling 537).

**Negative controls.**
- **(N0), asserted (rule (e)), on `th:3`'s first `e_real` call and `mag:5`'s penalty `B` call** (both DG; same form as t1 / t2):
  - V is the lexicographically smallest vertex of probe 0's owner cell, with coordinates shared by `Allreduce` SUM;
  - `n_inc(V) ≥ N0_MIN_INCIDENT`, and max/min `|f|` over V's incident cells (each rank's `_owned_incident_cells`, then `allgather`) > 1 + `N0_SPREAD_FLOOR`.
  - *Predicted:* > 1.007 (t1's sphere) and > 1.01 (DG B; t2 measured 4.2).
  - It is printed only for the other DG calls and skipped for CG1 calls, since a continuous field has no jump.
- **(N1), printed:** the legacy cell at V against the owner, at each width.

**R1 consequence (`mag:1`), declared before the run.**
- The switch from CG1 to DG moves the guide-table rows "Relative L2 error 51.9781 %", "Max relative error 76.7330 %" and "Decay ratio 20.31" (`01_straight_wire.md:17–23`), plus the per-radius ratio prose (`:126–128`), by construction. No assertion reads them, so this is **not** a red, and no `PREDICTED RED` line is written.
- **Re-registered only under R1:** the queueing review's ruling, written beneath this brief. Re-record from the `-n 2` log, keeping the superseded CG1 digits in an HTML comment as `EX-30` did (`01_straight_wire.md:25–36`). Also re-sync `:107–114`, which says the digits "reproduce … at fixed rank count".
- The energy row 2.630243e-08 J and the mesh row 21 830 / 4 662 cannot depend on sampling. A move in either is a finding: stop.
- *Predicted, printed:* relL2 within ±10 points of 51.98 %, direction not predicted.
- Without R1, nothing is re-recorded.

**Windows.** Each window runs in the foreground with Bash `timeout: 660000`, pytest `-s`, and container `timeout -k 30`. Each harness call is the whole bare Bash command. Before every compute window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
0. **Pre-checks:** `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up; `pgrep` reads 0.
1. **Pre-censuses:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-pre-figures "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py'"`
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
   Exit 2 means stale, exit 1 is a failure. **Predicted delta:** none (no example and no figure added).
2. **Dry run (prints lines only; host-side, never inside `docker compose exec`, never without `--dry-run`):**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-dryrun "./scripts/run_examples.sh -e 1,5,th:3,th:4,th:6,th:7,th:9,mat:2,mat:3,mri:2 -n 2 -t 180 --dry-run"`
   It confirms each example's build (mag real; th / mat / mri complex) and environment. If a printed line carries an environment variable the commands below lack, add it to that window verbatim and journal it.
3. Code: the edits and the module.
4. **Collect, both builds:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-W0-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py --collect-only -q -s'"`
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-W0-complex "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py --collect-only -q -s'"`
5. **Runs of record, one example per window, `-n 2`.** Each writes its npz.
   - Real (`mag1`, then `mag5`):
     `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-mag1-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_POST8_EXAMPLES=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py::test_mag1 -v -s --tb=short'"`
   - Complex (`th3`, `th4`, `th6`, `th7`, `th9`, in that order):
     `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-th3-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_th3 -v -s --tb=short'"`
   - The other five are these commands with the tag in the label and the node ID changed.
6. **`-n 4` groups:**
   - `POST-8-step3b-ii-e1-mag-n4`: the real form with `-n 4` and node IDs `::test_mag1 ::test_mag5`, container 180.
   - `POST-8-step3b-ii-e1-th-n4`: the complex form with `-n 4` and the five th node IDs in order, container 180.
7. **`-n 1` groups:**
   - `POST-8-step3b-ii-e1-mag-n1`: container 180.
   - `POST-8-step3b-ii-e1-th-n1`: container **400**, heavy by ceiling.
8. **Re-gates, bare examples, `-n 2`, complex:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e1-mat2-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 2 python3 examples/materials/02_lossy_sphere_sar.py'"`
   - `…-mat3-regate`: `examples/materials/03_gelled_saline_presets.py`, container 300.
   - `…-mri2-regate`: `examples/mri/02_mass_averaged_sar.py`, container 180.
9. **Guides.** Re-sync the seven guides: one sentence at the step that reads the sampled values, naming the owner rule, the single-incidence count and the `-n 2` log by full filename. `mag:1` also gets the R1 table work. Then run the post-censuses (window 1 relabelled `…-post-figures` / `…-post-docrefs`): the same counts, `broken=0`, `dead=0`.
10. `git status --porcelain examples/ tests/ src/ scripts/ docker/` lists only the seven scripts, seven guides and the new module. Run `scripts/testing/check_private_leak.py` on the staged logs.

**Traps already paid for.**
- **Batch composition (item 537's deviation, `attempts.md` ≈ 1698).** `.eval` bytes depend on the batch, at ULP level. Every byte comparison ((E), (S)) is made **call by call on that call's own point array**. Never concatenate or split captured batches, and never compare against a re-sample built differently. The (O) spread is per call.
- **Rule (i):** an edit to any edited example or to the module after window 5 starts means windows 5–7 again (all spied windows). An edit to a guide does not.
- **`./scripts/run_examples.sh --dry-run` only prints lines.** Every example runs in its own bare harness call above, never through the runner.
- **Never set `FEM_EM_SETUP_FIGURES`.** No figure changes, and a re-render dirties tracked PNGs. The examples' `paraview_output/` directories are gitignored, including `mag:1`'s cwd-relative `/workspace/paraview_output/`.
- **runpy hygiene:**
  - patch both import paths before `run_path`;
  - set `sys.argv` (`mag:4`-style argparse elsewhere would otherwise read pytest's argv);
  - a nonzero `SystemExit` is a red;
  - `undo()` in `finally`, so the grouped `-n 4` / `-n 1` windows start each example clean.
- **Rank safety:** every spy collective runs on all ranks inside the example's own collective call. Only rank 0 prints. `n_inc` and `.eval` are rank-local until reduced. Ownership comes only from `point_owner_cells` (buffer `Allreduce` MIN).
- **Complex build:** `tests/environment/test_complex_mode.py` goes first. The DG0 index is complex, so compare `.real` and assert `.imag == 0`. Use no `ufl` ordering on complex operands.
- **The bash hook** denies non-harness pytest and `pgrep -f "python3 -m pytest"`. `pgrep -c python3` is the check that runs.
- **Use new npz names only.** Never touch `step3a-n2.npz`, `step3b-n2.npz` or `step3b-ii-t1-*` / `t2-*`.
- **A kill is not a result.** Run `pgrep -c python3` first, then at 0 the rubric cache sweep as `POST-8-step3b-ii-e1-cache-sweep`. Never raise a timeout or `-n`.
- **Leak guard:** spy prints at ≤ 6 significant digits, ratios and spreads `.3e`. On a flag, move the log unchanged to `logs/leak-held/`, mark BLOCKED and report (item 526 precedent).
- Never load `ans:2`; never open `docs/private/` or `aed_results/`; no AED number. In guides, `≈`, never a tilde; full filenames only.

**Scope.** Seven example call sites, one new test module, seven guides and three re-gates. **Not claimed:**
- item 548's five sites;
- `mri:1` `:486–487` (outside the 17, for the weekly);
- 3b-iii;
- the `src/` legacy callers outside the census;
- flipping `evaluate_vector_field_parallel`'s default;
- removing the CG1 export's own partition dependence at `th:4` `:375` / `:380` (step 1's mechanism; the rule cannot reach inside an interpolant);
- any band or record;
- `POST-8` ✅.

**Negative result.** Report and stop. Park on `attempt/POST-8-3b-ii-e1-*`, mark item 547 BLOCKED in the same commit (rule (d)), and never loosen a band or re-record in the slot.
- **(X) red at `-n 2`:** stash that example's edit and run it **bare** at `-n 2` as `…-<tag>-baseline` (the window-8 form).
  - Baseline red: the failure predates the item; open a known-issues entry.
  - Baseline green: the rule moved an asserted record. Journal per probe `n_inc`, the owner and legacy indices and the delta. The review re-poses the site (a `(1*)` re-record by `record-reconciler`, or a source switch); never re-record here.
- **(X) red only at `-n 4` / `-n 1`:** the example's band is not width-free. Open a known-issues entry; the review re-poses.
- **(E) red:** the edit is not the one ruled. Park.
- **(S) red:**
  - if every red probe has owner index = legacy index and differs by a few ULP, it is the batch effect inside one call. Report it; the review re-poses (S)'s byte form;
  - otherwise it is a rule-path defect: new known-issues entry.
- **(O) red:** the indices are not partition-free. New known-issues entry.
- **(N0) red:** park for a re-pose.
- **A re-gate red, or a printed digit moved against its guide:** annotate known-issues 2026-09-27 (3b-i's flip reached an example). The review rules.
- **Census `broken` or `dead` > 0 from this item's own text:** fix the reference. Otherwise report.
- **Status 124:** not a result; park with the measured time.

There is no `PREDICTED RED:` line; a red anywhere stops the slot.

**Status it can move.** `POST-8` stays 🟡: the Done-when (`PROJECT_PLAN.md:2670`) still needs 548 and 3b-iii.
- The §7 state line's "12 examples open" becomes "3b-ii examples e1 done; 5 examples open (`mag:2`, `mag:4`, `mag:6`, `ports:1`, `ports:12`)".
- It gains: "**Step 3b-ii-e1 done** <date> (§9 item 547): seven example sites (`mag:1`, `mag:5`, `th:3/4/6/7/9`) sample under the owner rule [`mag:1` on the DG source, R1]; each example's own assertions green at `-n 1 / 2 / 4`; <k>/<N> probes single-incidence, rule = legacy bytewise there; owner indices exact across widths, spread <max> printed; `mat:2` / `mat:3` / `mri:2` re-run green after 3b-i (<logs>)."
- Known-issues 2026-09-27's "Verified at" gains e1, and the entry stays OPEN until 3b-iii.

**Census windows.** Setup figures: pre and post, counts unchanged, `broken=0`. Docrefs: pre and post, `dead=0` both sides. Journal both readings.

**Execution order.** Pre-checks → pre-censuses → dry run → code → W0 real / complex → `-n 2` × 7 → `-n 4` × 2 → `-n 1` × 2 → re-gates × 3 → guides → post-censuses → `git status` / leak check → commit → journal.

**Commit list:**
- the seven example scripts and seven guides;
- `tests/post/test_post8_step3b_ii_example_sites.py`;
- the harness logs (the `test-results.md` rows are the harness's own);
- `PROJECT_PLAN.md`: the §7 `POST-8` sentence and §9 item 547 DONE;
- `docs/testing/known-issues.md`: the 2026-09-27 "Verified at";
- this file's outcome line.

The journal goes in `attempts.md` in its own commit.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 547 (step 3b-ii-e1) per `docs/planning/chunks/POST-8.md` (heading `§9 item 547`) and the review ruling beneath it, R1 included.
- Read the brief and the ruling first; the ruling holds.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30`, in the brief's order. `pgrep -c python3` must read 0 before each compute window.
- `./scripts/run_examples.sh` is used only with `--dry-run`.
- Edit only the seven named call sites (and `mag:1`'s comment under R1), their guides, and the new module. Never set `FEM_EM_SETUP_FIGURES`. Never re-record anything except `mag:1`'s table under R1.
- Compare bytes only call by call on the same point array.
- If the docker socket is denied, stop and report the exact command string.
- If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e1-*`, mark item 547 BLOCKED in the same commit, and report.
- Return the log names with elapsed times, per example the call count / single-incidence count / (O) result / spread, the (N0) readings, the re-gate digits against their guides, and the census before and after."

**Questions for the review.**
- **R1 (`mag:1`).** Switch the sample from the CG1 `B_lag` to the DG source `B` under the rule, and license the guide-table re-record above? The CG1 interpolant carries step 1's partition dependence, which the rule cannot remove, so the current 51.9781 % is a `-n 2`-specific figure. If declined: keyword only, table held.
- **R2 (`mag:6`, item 548).** The sample is the CG1 export on purpose. The drafter recommends the keyword only, with no switch. This is the analogue of t2's `:149` / `:239` ruling, except that here the CG1 value **is** the asserted object.
- **R3.** The spied `-n 2` window is each edited example's run of record. Its own assertions run as `__main__`, and the spy only forwards and adds side calls. There is therefore no separate bare run, except for the three unedited re-gates.
- **R4.** `mri:1` `:486–487` is named for the weekly.
- **R5.** The re-classification table above supersedes item 526's α / β / γ labels for the examples.

**Outcome 2026-10-19 (09:00 CDT slot, attempt 1): BLOCKED, parked on `attempt/POST-8-3b-ii-e1-20261019T140700Z`.** Edits to the seven sites and the new module as briefed (R1 declined: `mag:1` keyword on `B_lag`, nothing re-recorded). `-n 2` runs of record all green, (X) / (E) / (S): call counts `mag:1` 1, `mag:5` 2, `th:3` 4, `th:4` 6, `th:6` 4, `th:7` 4, `th:9` 4 (all as predicted); single-incidence probes 29/30, 16/16, 54/54, 160/214, 100/100, 100/100, 192/192 (the multi-incidence probes: `mag:1` point 0, n_inc 24; `th:4` probe 20 of the DG calls and all 25 profile probes, n_inc 2, owner = legacy at the profile calls). (N0) `mag:5` 6.329747 (n_inc 22), `th:3` 1.000082 (n_inc 26; asserted floor 1 + 1e-6 holds, the 1.007 prediction missed). `-n 4`: mag group 2 passed, (O) exact, spreads 1.681e+00 (`mag:1`, CG1 interpolant) and 3.827e-06 (`mag:5`); th group red on `th:4` only (`20261019T140524Z_POST-8-step3b-ii-e1-th-n4.log`, Status 1, 36 s): `PROFILE_RMS_MAX` reads 2.25 % against 0.200 % at `-n 2`. The unedited example is red the same way at `-n 4` (baseline `20261019T140617Z_POST-8-step3b-ii-e1-th4-baseline-n4.log`, 5 s) — a failure that predates the item; known-issues entry opened. The `-n 1` groups and the three re-gates were not run (red stops the slot); the review re-poses.

### §9 item 548 — 2026-10-19 (lane 1, POST-8 3b-ii-e2)

**Outcome 2026-10-19 (13:30 CDT slot, attempt 1): BLOCKED, parked on `attempt/POST-8-3b-ii-e2-20261019T183900Z`.** Five keyword edits and five tests added to 547's module (existing code unchanged; e2 block uses its own `step3b-ii-e2-<tag>-n<w>.npz` and `_drive_e2`, since `_drive` would write and read the e1 stem). Collect 12 items on both builds. `mag:2` green at `-n 2` (`20261019T183230Z_POST-8-step3b-ii-e2-mag2-n2.log`, 156 s): 1 call, 25/25 probes single-incidence, (N0) 1.300640 (n_inc 18, owner = legacy cell 21206). `mag:4` red (`20261019T183519Z_POST-8-step3b-ii-e2-mag4-n2.log`, 100 s, Status 1): the example ran to completion, then (S) failed at call 0, probe 20, a single-incidence probe where owner and legacy sample bytes differ in the last mantissa bits of components 0 and 1 (1 ULP; component 2 equal). The assertion stops at the first red probe, so the owner and legacy indices at probe 20 and the other red probes were not printed. Per the 547 disposition this is either the batch effect inside one call (owner index = legacy index, review re-poses (S)'s byte form) or a rule-path defect; undetermined. `mag:6`, `ports:1`, `ports:12` not run (a red stops the slot). Censuses before: figures status 0, docrefs `dead=0 guide=0 stale=123` exit 2. Dry run emitted the brief's lines (mag real, ports complex).

**Executor: `example-runner`.** *Drafted read-only with item 547 for the 2026-10-19 03:00 daily review. It sets and moves no band, tolerance, Done-when or target.*

**Licence:** as item 547, including its site table and R2 / R3 / R5.

**`POST-8` step 3b-ii-e2: the five expensive example sites (`mag:2`, `mag:4`, `mag:6`, `ports:1`, `ports:12`) sample under the owner rule. Each edited example runs once at its record width under item 547's spy: its own assertions hold at their imported bands, and the rule is bitwise equal to legacy at single-incidence probes (call by call, same batch).**
- **Executor:** `example-runner`.
- **Tier:** heavy (each window ≈ 115–175 s; container 400 by ceiling). Ranks `-n 2`, except `ports:12` at `-n 4`, its record width. Measured basis:

  | Tag | Seconds | Source |
  |---|---|---|
  | `mag:2` | 136 | `20260916T094228Z_OPS-48-mag2.log` |
  | `mag:4` | ≈ 134 | 270 for `mag:2` + `mag:4` (`20260826T170305Z_EX-30-root2-run-mag2to4.log`) less `mag:2`; the guide's rungs sum to 96.4 s of mesh + solve, `04_helmholtz_analytic_comparison.md:49–51` |
  | `mag:6` | 163 | `20260826T170746Z_EX-30-root2-run-mag6.log` |
  | `ports:1` | 133 | `20260927T124920Z_EX-57-control.log` |
  | `ports:12` | 113 at `-n 4` | `20260928T112745Z_EX-57.log` |

  Predicted windows ≈ 735 s for the five spied runs, plus ≈ 60 s for collects, dry run and censuses: **≈ 13.5 min**.
- **Class:** lane 1, example. **Cost:** ≈ 8.4 measured / ≈ 45 predicted slot-minutes: fixed ≈ 5, code ≈ 8, windows ≈ 14, plus per-call overhead ≈ 3, guides ≈ 7, commit and journal ≈ 6.
- **Slot-first only.** **Park (rule (d)):**
  - the `mag:2` window not started by minute 20;
  - the `ports:1` window not started by minute 38;
  - the `ports:12` window not started by minute 44.
- **Depends on:** **547 (landed).** This item adds tests to 547's module and reuses its spy.
- **Family:** `POST-8` 3b-ii-e2; attempt 1.

**Code.**
- Add `owner="lowest_original_cell"` at `mag:2` `:269`, `mag:4` `:200`, `mag:6` `:212–214` (CG1, keyword only, R2), `ports:1` `:303`, and `ports:12` `:503–508` (two calls).
- Add `test_mag2`, `test_mag4`, `test_mag6`, `test_ports1`, `test_ports12` to `tests/post/test_post8_step3b_ii_example_sites.py`. **Add only.** If a 547 helper must change, stop and park: it would void 547's windows under rule (i).
- No other edit.

**Anchors (asserted):** item 547's (X), (E) and (S), at record width only.
- **(E)** *predicted* counts: `mag:2` 1, `mag:4` 3 (one per rung), `mag:6` 1, `ports:1` 8 (2 columns × orders 2049 / 4097 × 2 ports), `ports:12` 2.
- **(X)** covers `ports:1`'s four 0.11-image records at the 1e-6 band (`:739`) and its controls (`:750–771`), `ports:12`'s C16 identity and 3-step control (`:466–499`), and `mag:6`'s `:272` / `:285`.
- **Printed beside the guides' un-asserted tables, never re-recorded:** `mag:2`'s, `mag:4`'s centre `B_z` / errors / CV (`04_…md:47–51`), and `ports:12`'s ring RMS / worst.
- Each run writes `step3b-ii-e2-<tag>-n<w>.npz` for a future width check. There is **no (O)** in this item.

**Predictions, printed.** Single-incidence at every probe:
- `mag:2` / `mag:4` sample the z-axis of a gmsh sphere or box with no axis entity; t2's coil-phantom centerline read 9/9;
- `mag:6`'s ten +x points;
- `ports:1`'s 24 584 arc Gauss points, off any meshed curve;
- `ports:12`'s ring, unless z = 0 is a fragment plane of the birdcage mesh. In that case the rows go to (R) and the ring print may move. It is ungated, so it is re-registered in the guide only under the review's ruling.

**Negative controls.** (N0), asserted as in 547, on `mag:2`'s `B` (DG) and on `ports:12`'s first `e_complex` call (raw N1curl, at a vertex every component can jump); predicted > 1.01. Printed for `mag:4` / `ports:1`; skipped for `mag:6` (CG1). (N1) printed.

**Records at risk, declared.**
- **`ports:1`'s 1e-6 reproduction band.** A multi-incidence arc point whose owner cell differs from the legacy cell would shift V by roughly a quadrature weight times the normal-component jump, which could exceed 1e-6. The census does not predict one, so there is **no `PREDICTED RED`**.
  - If it happens it is the (X)-red disposition: baseline bare run; on a green baseline, journal the probes. The record is re-registered **only under the review's ruling**, as a `(1*)` re-record by `record-reconciler`, never in the slot.
- `mag:2` / `mag:4` / `ports:12`: un-asserted tables move only at a multi-incidence probe, and then likewise only under ruling.

**Windows.** The forms, `pgrep` checks and timeouts are item 547's. Each runs as one bare harness call per example, in this order:
0. **Pre-checks.** As in 547.
1. **Censuses.** `POST-8-step3b-ii-e2-pre-figures` / `…-pre-docrefs`, item 547's commands relabelled.
2. **Dry run (prints lines only):**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2-dryrun "./scripts/run_examples.sh -e 2,4,6,ports:1,ports:12 -n 2 -t 400 --dry-run"`
3. **Code.**
4. **Collect.** `POST-8-step3b-ii-e2-W0-real` / `…-W0-complex`, item 547's collect commands relabelled.
5. **The three magnetostatic runs, real build, container 400:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2-mag2-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py::test_mag2 -v -s --tb=short'"`
   Then `…-mag4-n2` (`::test_mag4`) and `…-mag6-n2` (`::test_mag6`).
6. **`ports:1`, complex, container 400:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2-ports1-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_ports1 -v -s --tb=short'"`
7. **`ports:12`:** the same form with `-n 4`, `::test_ports12`, label `POST-8-step3b-ii-e2-ports12-n4`.
8. **Guides.** One sentence each in the five guides. Then the post-censuses (counts unchanged, `broken=0`, `dead=0`), then `git status --porcelain examples/ tests/ src/ scripts/ docker/` (only the five scripts, five guides and the module) and the leak check.

**Traps already paid for.** All of item 547's: batch composition (byte checks call by call on the same array), rule (i) (an edit after window 5 starts means windows 5–7 again), dry-run only, never `FEM_EM_SETUP_FIGURES`, runpy hygiene, rank safety, complex-build order, the bash hook, new npz names, a kill is not a result plus the cache sweep, the leak guard, no AED number.
- **The ≈ 10-minute line:** every window here is ≤ 400 s, so no durable capture is needed. After any kill, check for orphaned ranks before anything else.
- **`ports:1` prints full-precision `(1*)` record candidates** (`:694–699`), which passed the leak guard on 09-27. The spy adds nothing at more than 6 significant digits.
- **`mag:2` passes no `comm`:** the spy's `comm` fallback (item 547) is load-bearing here.

**Scope.** Five example call sites, five tests added to 547's module, five guides. **Not claimed:**
- width invariance (no (O)) on these five fixtures, on cost (t2's straight-wire precedent);
- `mri:1` (R4);
- 3b-iii;
- `src/` legacy callers;
- removing the CG1 export's partition dependence at `mag:6`;
- any band or record;
- `POST-8` ✅.

**Negative result.** As item 547's, with `attempt/POST-8-3b-ii-e2-*` and item 548 BLOCKED. The baseline for an (X) red is the bare `-n <record width>` run of the stashed example at container 400. There is no `PREDICTED RED:` line.

**Status it can move.** `POST-8` stays 🟡 (3b-iii).
- The §7 state line's "5 examples open" becomes "examples done".
- It gains: "**Step 3b-ii-e2 done** <date> (§9 item 548): `mag:2`, `mag:4`, `mag:6` (CG1 export, keyword only, R2), `ports:1`, `ports:12` sample under the owner rule; each example's own assertions green at record width (`ports:1`'s four records inside 1e-6); <k>/<N> probes single-incidence, rule = legacy bytewise there; width invariance not run on cost. 3b-ii complete; 3b-iii open."
- Known-issues 2026-09-27's "Verified at" gains e2, and the entry stays OPEN until 3b-iii.
- With 548 green, every census candidate except `ans:2` is on the rule, and the next review may say so in §2 / §7.

**Census windows.** As 547: pre and post setup-figure counts unchanged with `broken=0`; docrefs `dead=0` both sides.

**Execution order.** Pre-checks → pre-censuses → dry run → code → W0 × 2 → `mag:2` → `mag:4` → `mag:6` → `ports:1` → `ports:12` → guides → post-censuses → `git status` / leak check → commit → journal.

**Commit list:**
- the five example scripts and five guides;
- the module (tests added);
- the harness logs;
- `PROJECT_PLAN.md`: the §7 sentence and §9 item 548 DONE;
- known-issues "Verified at";
- this file's outcome line.

The journal goes in `attempts.md` in its own commit.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 548 (step 3b-ii-e2) per `docs/planning/chunks/POST-8.md` (heading `§9 item 548`, and item 547's sections it cites) and the review ruling beneath it.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 400`, in the brief's order. `pgrep -c python3` must read 0 before each compute window.
- `./scripts/run_examples.sh` is used only with `--dry-run`.
- Edit only the five named call sites, their guides, and **add** five tests to `tests/post/test_post8_step3b_ii_example_sites.py`; change none of its existing code. Never set `FEM_EM_SETUP_FIGURES`. Never re-record anything.
- Compare bytes only call by call on the same point array.
- If the docker socket is denied, stop and report the exact command string.
- If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e2-*`, mark item 548 BLOCKED in the same commit, and report.
- Return the log names with elapsed times, per example the call count / single-incidence count, the (N0) readings, `ports:1`'s four record misses, and the census before and after."

### Review ruling 2026-10-19 03:00 (items 547, 548 — 3b-ii example half)

- **R1 declined: `mag:1` gets the keyword only, and its table is held.** Switching the sample from the CG1 `B_lag` to the DG `B` changes what the example measures and re-records three guide rows. That is a content change, not an owner-rule move. It is named for the weekly, beside `mri:1`. Item 547 therefore re-records nothing. The brief's "R1 consequence" block and the `mag:1` re-record clause in its spawn line do not apply.
- **R2 accepted:** `mag:6` gets the keyword only, with no switch (the t2 `:149` / `:239` analogue).
- **R3 accepted:** the spied `-n 2` window, running the unmodified example as `__main__` with its own assertions, is that example's run of record. The new module runs only under `FEM_EM_POST8_EXAMPLES=1`. Item 537's trap binds: every byte comparison is call by call, on that call's own point array.
- **R4 accepted:** `mri:1` `:486–487` is the weekly's.
- **R5 accepted:** the brief's re-classification table supersedes item 526's α / β / γ labels for the examples.
- **Ordering:** 548 depends on 547, because it adds tests to 547's module. If 547 is not DONE on `main`, skip 548 at no attempt cost.
- **`ports:1`'s 1e-6 records:** a red there is a reported negative (known-issues entry, item BLOCKED). It is never a re-record and never a band.

### Review ruling 2026-10-19 10:30 (item 547 blocked — diagnosis and re-pose)

- **Diagnosis (log-pathologist ruling, journal `docs/planning/reviews/2026-10-19-daily-1030.md` §4).** `th:4`'s `PROFILE_RMS_MAX` red at `-n 4` predates the item. The unedited example is red with identical digits, 2.249 % (`20261019T140617Z_POST-8-step3b-ii-e1-th4-baseline-n4.log:50`; the traceback lines 487 / 410 match the unedited source, `:52–67`). The edited `-n 2` window is green at **0.124 %** (branch log `…-th4-n2.log:78`), not the 0.200 % the slot quoted: 0.200 % is the August `-n 2` record (`20260810T020325Z_EX-7-run1.log.gz:53`). At the 25 profile probes the owner and legacy values agree to ≤ 1.016e-15 (`…-th4-n2.log:91`; bitwise at 18, not all 25), so the rule does not move the profile. The profile samples `e_cg.interpolate(fields.e_complex)` (`:250`), a direct N1curl → CG1 vertex interpolation (there is no DG split in this example). Its vertex values depend on the partition: plausible from the code, but not established by a vertex dump. The N1curl readings are width-free (γ 37.650399, L2 4.406648e-02 at both widths). `CG1_VS_NEDELEC_MAX` also moves with width, 0.027 % → 0.490 % against its 0.5 % bound, so it nearly flips at `-n 4`. `th:4`'s CG1 bands are therefore records at the default width, `-n 2` (`run_examples.sh:25`), not width-free gates. No tracked text states that width, and today's `-n 2` figure has moved from August's for an unidentified reason (the N1curl readings have not). Known-issues 2026-10-19 (`th:4`) stays OPEN; its fix (sample the profile from the DG split under the owner rule, or rule the band width-dependent) is a content change and goes to the weekly beside `mag:1`'s R1.
- **The two cited logs are landed on `main`** by this review (path checkout of two new log files, their two `test-results.md` rows appended by hand, insertions only).
- **Re-pose: item 549 (3b-ii-e1, attempt 2 of the family since t2's gate).** Same target, same seven sites, same module, same anchors, same bands. One change: **`th:4` runs at its record width `-n 2` only**; the `-n 4` and `-n 1` th groups carry `test_th3 test_th6 test_th7 test_th9`. `th:4`'s (O) across widths is not claimed. Brief: `§9 item 549` below.
- **The parked branch stays** until 549 lands; 549 takes its code by path and deletes it.
- **Item 548 now depends on 549** (it adds tests to the same module). If 549 is not DONE on `main`, skip 548 at no attempt cost. 548's brief is otherwise unchanged.

### §9 item 549 — 2026-10-19 (lane 1, POST-8 3b-ii-e1, re-pose of item 547)

**Outcome 2026-10-19 (12:00 CDT slot, attempt 2): done.** Code by path from the parked branch (eight files, keyword-only script edits, `mag:1` on `B_lag`); every window on that tree. `-n 2` spied runs of record (calls / single-incidence): `mag:1` 1 / 29 of 30, `mag:5` 2 / 16 of 16, `th:3` 4 / 54 of 54, `th:4` 6 / 160 of 214, `th:6` 4 / 100 of 100, `th:7` 4 / 100 of 100, `th:9` 4 / 192 of 192; (X) each example's assertions green (`th:4` profile 0.124 %). (O) owner indices and `n_inc` exact at `-n 4` / `-n 1`; field spread (printed) `mag:1` 1.681 / 0.6965, `mag:5` 3.284e-06 / 2.199e-05, `th:3` 2.503e-09 / 4.327e-09, `th:6` 5.206e-12 / 2.817e-12, `th:7` 2.279e-11 / 2.526e-11, `th:9` 3.957e-09 / 3.560e-09; `th:4` not run at `-n 4` / `-n 1`. (N0) `th:3` 1.000082 at all three widths, `mag:5` 6.329709 / 6.329661 / 6.329869 at `-n 2 / 4 / 1` (predicted 6.329747 ±1e-6: missed in the sixth digit, solve spread, no band). Re-gates Status 0: `mat:2` and `mat:3` digit-identical to their guides; `mri:2` differs at two printed rows (DG0 sum 2.81e-15 against the guide's 1.32e-15, both under 1e-12; 10 g kernel mass 0.0041 % against 0.0044 %, the same 0.0041 % as the 10-18 re-gate log). Censuses: figures `ok=131 broken=0` both times; docrefs `dead=0 guide=0`, exit 2 both times, stale 127 → 123 (info). Logs `20261019T170143Z_POST-8-step3b-ii-e1b-mag1-n2.log` (8 s), `…170159Z_…mag5-n2.log` (6 s), `…170214Z_…th3-n2.log` (12 s), `…170233Z_…th4-n2.log` (9 s), `…170247Z_…th6-n2.log` (27 s), `…170321Z_…th7-n2.log` (13 s), `…170341Z_…th9-n2.log` (5 s), `…170352Z_…mag-n4.log` (10 s), `…170410Z_…th-n4.log` (31 s), `…170448Z_…mag-n1.log` (36 s), `…170530Z_…th-n1.log` (50 s), `…170626Z_…mat2-regate.log` (51 s), `…170723Z_…mat3-regate.log` (61 s), `…170830Z_…mri2-regate.log` (10 s).

**Executor: `example-runner`.** *Written by the 2026-10-19 10:30 daily review as a re-pose of item 547 under rule (6) (same target, next attempt of the family). It sets and moves no band, tolerance, Done-when or target.*

**This brief is item 547's brief (heading `§9 item 547` above) with the deltas below. Read 547's brief, then the 03:00 ruling (R1 declined, R2–R5), then the 10:30 ruling, then this.** Where they differ, this brief binds.

**Deltas.**
1. **Code from the parked branch, by path.** `git checkout attempt/POST-8-3b-ii-e1-20261019T140700Z -- examples/magnetostatics/01_straight_wire.py examples/magnetostatics/05_gauge_cross_check.py examples/time_harmonic/03_dielectric_sphere_in_uniform_field.py examples/time_harmonic/04_evanescent_waveguide_decay.py examples/time_harmonic/06_larmor_lossy_sphere.py examples/time_harmonic/07_element_order_lossy_sphere.py examples/time_harmonic/09_pec_sphere_as_hole.py tests/post/test_post8_step3b_ii_example_sites.py`. **Never** path-checkout `docs/testing/test-results.md`, `PROJECT_PLAN.md`, `known-issues.md` or `POST-8.md` from the branch (rubric trap: an append-only record checked out by path deletes `main`'s later rows). Read the diff (`git diff --cached --stat` shows the eight files only) and confirm each script edit is the `owner="lowest_original_cell"` keyword alone (`mag:1` on `B_lag`, R1 declined). No other edit.
2. **All windows re-run on the committed-to-be code** (rule (i)): the branch's logs ran on an uncommitted tree; this item's own windows are the runs of record. Labels `POST-8-step3b-ii-e1b-…`. The module's npz stem stays `step3b-ii-e1-<tag>-n2.npz` (`:285`): those files are this family's own, written by the same code, and this item's `-n 2` windows overwrite them before any `-n 4` / `-n 1` window reads them. No module edit.
3. **Windows 5–7 of 547's list, changed:**
   - window 5: the seven `-n 2` spied runs of record, unchanged (mag1, mag5, th3, th4, th6, th7, th9);
   - window 6: `…-mag-n4` unchanged; `…-th-n4` with node IDs `::test_th3 ::test_th6 ::test_th7 ::test_th9` (**no `test_th4`**), container 180;
   - window 7: `…-mag-n1` unchanged; `…-th-n1` with the same four th node IDs, container 400 (heavy by ceiling).
   Every window must read N passed, 0 skipped, with N the node count it names.
4. Windows 8–10 (re-gates `mat:2` / `mat:3` / `mri:2`, guides, post-censuses, `git status`, leak check) unchanged. The `th:4` guide sentence names the `-n 2` log and says the profile band is a `-n 2` record (known-issues 2026-10-19).
5. **Finish:** `git branch -D attempt/POST-8-3b-ii-e1-20261019T140700Z` after the outcome commit (the sandbox warns on `.git/config`; the delete still happens).

**Anchor (unchanged):** (X) each example's own assertions at imported bands under the spied `__main__` run; (E) every recorded call carries the rule, returned bytes = `own`; (S) bitwise across ranks, `own` = `leg` bytes at single-incidence probes call by call on the same batch, in-cell equality at multi-incidence ≤ `CLOSED_FORM_TOL`; (O) owner original indices and `n_inc` exactly equal at `-n 4` / `-n 1` against the `-n 2` npz for `mag:1`, `mag:5`, `th:3`, `th:6`, `th:7`, `th:9`. The (O) field spread is printed only.

**Negative control (asserted, rule (e)):** (N0) on `th:3`'s first `e_real` call and `mag:5`'s penalty `B` call, max/min `|f|` over V's incident cells > 1 + `N0_SPREAD_FLOOR`. Asserted from the same comparison on the same fixture: item 547's `-n 2` windows read 1.000082 (`th:3`) and 6.329747 (`mag:5`) (attempts.md 2026-10-19T14:12Z; branch logs `…-th3-n2.log`, `…-mag5-n2.log`). *Predicted, printed:* the same digits ±1e-6.

**Tier and cost.** Standard, two windows heavy by ceiling; `-n 1 / 2 / 4`. Measured (the branch's `test-results.md` rows): item 547's seven `-n 2` windows 4–22 s each (`th:6` the longest), the `-n 4` mag group 10 s, the `-n 4` th group 36 s (with `th:4`), the baseline 5 s; censuses, dry run and collects ≤ 3 s each. The `-n 1` groups and re-gates are unmeasured (547's predictions ≈ 40 s / ≈ 160 s / ≈ 120 s stand). Predicted windows ≈ 9 min (547's ≈ 10 less `th:4`'s two widths). **Class:** lane 1, example. **≈ 10.0 measured (example median) / ≈ 40 predicted slot-minutes** (no code to write: checkout ≈ 3, windows ≈ 9, overhead ≈ 5, guides ≈ 8, commit and journal ≈ 6, fixed ≈ 5, margin 4).

**Slot-first only** (never take-next). **Park (rule (d)):** the first `-n 2` window not started by minute 15; the `-n 1` th group not started by minute 35; the re-gates not started by minute 42. **Depends on:** nothing. **Family:** `POST-8` step 3b-ii-e1; attempt 2 since t2's gate (cap 4).

**Traps.** Item 547's list binds in full (batch composition; rule (i); dry-run prints only; never `FEM_EM_SETUP_FIGURES`; runpy hygiene; rank safety; complex build ordering; the bash hook; new npz names only; a kill is not a result; leak guard; no `ans:2`, no private paths). Added: **the path checkout of an append-only record** (delta 1); **`th:4` at `-n 4` is a known red** (known-issues 2026-10-19) — a slot that runs it anyway has spent a window on a known answer.

**Scope.** As 547's, plus: **not claimed — `th:4`'s width independence**, its profile band at any width other than `-n 2`, and the known-issues 2026-10-19 entry's resolution.

**Negative result.** As 547's list. A red at `-n 4` / `-n 1` on any of the other six examples is a new finding: known-issues entry, item BLOCKED, park on `attempt/POST-8-3b-ii-e1b-*`. There is no `PREDICTED RED:` line.

**Status it can move.** As 547's (`POST-8` stays 🟡; the §7 state line gains "**Step 3b-ii-e1 done**", with "`th:4` at `-n 2` only (known-issues 2026-10-19)" added to its width clause). Known-issues 2026-09-27's "Verified at" gains e1. Item 548 becomes runnable.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 549 (step 3b-ii-e1, re-pose of item 547) per `docs/planning/chunks/POST-8.md`: read the brief `§9 item 547`, the 2026-10-19 03:00 ruling, the 2026-10-19 10:30 ruling and the brief `§9 item 549`, in that order; 549's deltas bind.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Take the code from the parked branch by path (the eight files named in delta 1 only); never path-checkout a record file.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30`. `pgrep -c python3` must read 0 before each compute window.
- `th:4` runs at `-n 2` only; the `-n 4` / `-n 1` th groups name `test_th3 test_th6 test_th7 test_th9`.
- `./scripts/run_examples.sh` only with `--dry-run`. Never set `FEM_EM_SETUP_FIGURES`. Re-record nothing.
- Compare bytes only call by call on the same point array.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e1b-*`, mark item 549 BLOCKED in the same commit, and report.
- Return the log names with elapsed times, per example the call count / single-incidence count / (O) result / spread, the (N0) readings, the re-gate digits against their guides, and the census before and after."

### Review ruling 2026-10-19 18:00 (item 548 blocked — diagnosis and re-pose)

- **Diagnosis (`log-pathologist` ruling; journal `docs/planning/reviews/2026-10-19-daily-1800.md` §4).** The log is admissible: footer Status 1, 100 s (`20261019T183519Z_POST-8-step3b-ii-e2-mag4-n2.log:755–758`), with the same tuple on both ranks (`:700–701`, `:732–733`).
  - **(1) CONFIRMED.** The red is (S) at call 0, probe 20 (`_drive_e2:410 → _check_e_s:204`, `:695–700`), after all three rungs of the example ran to completion (`:266`, `:445`, `:637`, `:678`). Components 0 and 1 differ by **exactly 1 ULP**; component 2 is equal.
  - **(2) Hypothesis A CONFIRMED from the code:** within one call, owner and legacy evaluate different inputs on this mesh. `two_torus_domain` is partitioned with `GhostMode.shared_facet` (`src/fem_em_solver/io/mesh.py:2306–2326`, docstring `:1371`); `birdcage_port_domain` is too (`:4911–4920`). The legacy batch on a rank is every point that hits a local cell, **owned or ghost**, and the gather keeps the last rank's value (`post/evaluation.py:379–388`, `:413–415`). The owner batch is only the points whose owner cell is owned on that rank (`:389–396`). So on a ghosted mesh a single-incidence probe can be evaluated in a different batch (A1), or on a different rank in a ghost copy (A2). Which of these, or a ghost-copy coefficient difference (A3), produced probe 20's ULP is UNCOUNTABLE.
  - **(3) Hypothesis B REFUTED** (from code plus bytes). `_n_incident` counts owned copies, so n_inc = 1 means every ghost copy is a copy of the one cell. A different cell on the coarse rung would move components 0 and 1 (≈ 2 % of |B|, discretisation-dominated) by far more than 1 ULP.
  - **(4)** The (S) loop stops at its first red, but probe 20 is call 0's last index (21 points). So **probe 20 is call 0's only single-incidence red**. Calls 1–2, the multi-incidence checks and (N0) are unread. *The journal's "other probes unknown" is corrected for call 0.*
  - **(5)** `mag:2` and every e1 site are ghost-free meshes (`circular_loop_domain` and the default partitioner, `mesh.py:778–780`). Their (S) greens are real for those meshes, but they cannot test A.
- **Ruling.** The owner rule is not shown to be wrong. **(S)'s premise, "call by call, same batch" (`:870`, `:963`), does not hold on a ghosted mesh**, where owner and legacy are not evaluated in the same batch even on one call. How (S) is asserted on ghosted fixtures (bytes in a same-batch re-evaluation, an in-cell tolerance, or per-copy singletons) moves an assertion's form, so it is **the weekly's**. This review does not move it. `ports:1` (`two_torus_domain`, `examples/ports/01_two_torus_port_pair.py:416`) and `ports:12` (if it is built by `birdcage_port_domain`) are predicted to go red the same way.
- **Re-pose (rule (6), target unchanged).** Split by mesh class:
  - **Item 551 (e2a, attempt 2 of the e2 family since 549's gate):** the ghost-free sites `mag:2` and `mag:6` under the unchanged (S). It also lands the e2 test code, with an e2-only checker that reports every red and the rank / ghost evidence.
  - **Item 554 (e2b, attempt 3, depends on 551):** the ghosted sites `mag:4`, `ports:1` and `ports:12`. Each runs its (X) at record width and the separating measurement, with (S) pre-registered as a predicted red. It hands the (S) form to the weekly with a full red census.
- **The parked branch stays** until item 554 lands. Item 551 takes the module and two script edits from it by path, and item 554 takes the other three.
- **Family count.** 548 was attempt 1, 551 is attempt 2 and 554 is attempt 3, against a cap of 4. A green 551 lands a gate (the §7 state line moves), which resets the count.

### §9 item 551 — 2026-10-19 (lane 1, POST-8 3b-ii-e2a, re-pose of item 548: ghost-free sites)

**Outcome 2026-10-20 (00:30 UTC slot, attempt 2): done.** Code by path from `attempt/POST-8-3b-ii-e2-20261019T183900Z` (module, `02_circular_loop.py`, `06_h_convergence_rate.py`; both script diffs the `owner="lowest_original_cell"` keyword alone). Module diff on top of the branch: added `_legacy_hits`, `_ulp_distance` and `_check_e2_s` (collects every single-incidence mismatch per call with call / probe / `n_inc` / owner and legacy indices / owning rank / per-rank legacy hits with owned-ghost flag / ULP distance, asserts none after the loop; prints `cells_ghost` per site and the (S-sb) same-batch control per call, never asserted); the one-line swap in `_drive_e2` (`_check_e_s` to `_check_e2_s`); `_check_e_s` and the e1 tests byte-identical. `cells_ghost` `[0, 0]` for both sites. `mag:2` (`20261020T003156Z_POST-8-step3b-ii-e2a-mag2-n2.log`, 139 s): 1 call, 25/25 single-incidence, 0 byte mismatches, (S-sb) 0, (N0) 1.300640 (n_inc 18, owner = legacy cell 21206), 1 passed. `mag:6` (`20261020T003421Z_POST-8-step3b-ii-e2a-mag6-n2.log`, 146 s): 1 call, 10/10 single-incidence, 0 mismatches, (S-sb) 0, (N0) skipped (CG1), 1 passed. Censuses: figures `ok=132 broken=0` before and after; docrefs `dead=0 guide=0`, exit 2 both times, stale 119 to 118 (info). Item 554 is runnable; the parked branch stays.

*Written by the 2026-10-19 18:00 daily review as a re-pose of item 548 under rule (6): same target, the next attempt of the family. It sets and moves no band, tolerance, Done-when or target.*

**Read, in this order:** item 547's brief, the 03:00 ruling, item 548's brief, the 18:00 ruling above, then this. This brief's deltas bind.

**Executor:** `example-runner`. **Tier:** heavy (two windows of ≈ 160 s each, container 400), `-n 2`. **Class:** lane 1. **Cost:** ≈ 8.8 measured / 40 predicted min. **Slot-first only.** **Depends on:** nothing (549 is DONE). **Family:** e2, attempt 2.

**Deltas from item 548.**
1. **Code from the parked branch, by path, three files only:**
   `git checkout attempt/POST-8-3b-ii-e2-20261019T183900Z -- tests/post/test_post8_step3b_ii_example_sites.py examples/magnetostatics/02_circular_loop.py examples/magnetostatics/06_h_convergence_rate.py`.
   - **Never** take `docs/testing/test-results.md`, a log, `PROJECT_PLAN.md` or a chunk file from the branch (rubric trap: an append-only record checked out by path deletes `main`'s later rows).
   - `git diff --cached --stat` must show those three files only, and each script diff must be the `owner="lowest_original_cell"` keyword alone.
2. **One module edit, e2 code only.** Item 547's code (`_check_e_s`, `_drive`, every e1 test) stays byte-identical, because changing it would void 549's windows under rule (i).
   - Add `_check_e2_s(tag, calls, comm)`. It is `_check_e_s`'s logic with two changes:
     - the single-incidence byte loop **collects every mismatch** in each call, then asserts that there are none, with the same assertion meaning;
     - for each mismatch it prints call, probe, `n_inc`, `owner_orig`, `legacy_orig`, the rank owning the cell, every rank whose legacy batch hit the point with an owned / ghost flag, and the per-component ULP distance (`np.frombuffer(…, np.int64)` difference).
   - Add one print per call: the **same-batch control** `(S-sb)`. On each rank, re-evaluate the legacy cells on exactly the owner batch's points and order, and print the byte-mismatch count against `own`. It is printed only, never asserted (rule (h): a comparand true by construction is a record).
   - Add one print per site: `cells_ghost` per rank for the sampled function's mesh (`f.function_space.mesh.topology.index_map(3).num_ghosts`).
   - `_drive_e2` calls `_check_e2_s` in place of `_check_e_s`. Disclose the diff in the journal.
3. **Mesh-class pre-condition (read by the executor from the module's print; the module prints it and never asserts it, because item 554 runs the same module on ghosted sites):** `cells_ghost == 0` on every rank. A site that reads > 0 is out of this item. Stop that site, journal the reading, and leave it for item 554. That is not a red and not a park.
4. **Windows:** item 548's pre-checks, censuses, dry run (`-e 2,6`) and collects (`…-e2a-W0-real` / `…-W0-complex`). Then:
   - `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2a-mag2-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py::test_mag2 -v -s --tb=short'"`
   - the same with `::test_mag6`, label `POST-8-step3b-ii-e2a-mag6-n2`.

   Each must read 1 passed, 0 skipped. Then guide sentences for `mag:2` and `mag:6`, post-censuses, `git status`, the leak check.
5. **Measured basis:** item 548's `mag:2` window, 156 s (`20261019T183230Z_POST-8-step3b-ii-e2-mag2-n2.log`); `mag:6` 163 s (`20260826T170746Z_EX-30-root2-run-mag6.log`).

**Anchor (unchanged from 548, two sites):**
- (X) each example's own assertions at imported bands under the spied `__main__` run;
- (E) `mag:2` 1 call, `mag:6` 1 call (predicted);
- (S) bitwise across ranks, `own` = `leg` bytes at single-incidence probes, and in-cell equality at multi-incidence ≤ `CLOSED_FORM_TOL`.

**Negative control (asserted, rule (e)):** (N0) on `mag:2`'s `B` > 1 + `N0_SPREAD_FLOOR`. Item 548 measured 1.300640 over 18 incident cells on the same fixture (`…183230Z…mag2-n2.log`, attempts.md 2026-10-19T18:41Z). `mag:6` (CG1) is skipped, as in 548.

**Predicted, printed:** `mag:2` 25/25 single-incidence, (S-sb) 0, `cells_ghost` 0 on both ranks; `mag:6`'s ten +x points single-incidence.

**Traps:** item 548's list, plus:
- the path checkout of an append-only record (delta 1);
- **a ghosted mesh breaks the "same batch" premise** (the 18:00 ruling): never read an (S) red on a ghosted site as an owner-rule defect;
- `_check_e_s` and every e1 line must not change (`git diff` on the module shows additions plus the one-line `_drive_e2` call swap only).

**Scope.** Two sites. Not claimed: `mag:4`, `ports:1`, `ports:12` (item 554); the (S) form on ghosted meshes (the weekly's); width invariance; `POST-8` ✅.

**Negative result.**
- An (S) red on a ghost-free site is a new finding. Run the bare-example baseline as 548 says, journal the red census that `_check_e2_s` prints, park on `attempt/POST-8-3b-ii-e2a-*`, and mark 551 BLOCKED in the same commit.
- An (X) red follows 548's disposition.
- There is no `PREDICTED RED:` line.

**Status it can move.** The `POST-8` §7 state line gains "**3b-ii-e2a done** <date> (§9 item 551): `mag:2`, `mag:6` on the owner rule, ghost-free, (S) bitwise; `mag:4` / `ports:1` / `ports:12` on ghosted meshes in item 554". Known-issues 2026-09-27's "Verified at" gains e2a. Item 554 becomes runnable.

**Commit:**
- the module and the two scripts;
- the two guides;
- the logs and their `test-results.md` rows (appended by hand);
- `PROJECT_PLAN.md`: the §7 state line, and §9 item 551 DONE;
- this file's outcome line.

The journal goes in its own commit. **Keep the parked branch** (item 554 takes three scripts from it).

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 551 (step 3b-ii-e2a) per `docs/planning/chunks/POST-8.md`: read item 547's brief, the 2026-10-19 03:00 ruling, item 548's brief, the 2026-10-19 18:00 ruling and item 551's brief, in that order; 551's deltas bind.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Take three files from the parked branch by path (delta 1); never path-checkout a record file.
- Edit only the e2 code in the module (delta 2); `_check_e_s` and the e1 tests stay byte-identical.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 400`. `pgrep -c python3` must read 0 before each compute window.
- `./scripts/run_examples.sh` only with `--dry-run`. Never set `FEM_EM_SETUP_FIGURES`. Re-record nothing.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e2a-*`, mark item 551 BLOCKED in the same commit, and report.
- Return the log names with elapsed times, per site `cells_ghost`, the call count / single-incidence count, (S-sb), the (N0) reading, and the census before and after."

### §9 item 554 — 2026-10-19 (lane 1, POST-8 3b-ii-e2b: ghosted sites, the separating measurement)

**Outcome 2026-10-20 (03:30 UTC slot, attempt 3): BLOCKED, unpredicted red; scripts parked on `attempt/POST-8-3b-ii-e2b-20261020T033900Z` (`ea0dc499`).** Three keyword-only script edits taken by path; module untouched. Censuses: figures `ok=132 broken=0` before and after; docrefs `dead=0 guide=0 stale=118`, exit 2 both times.
- **`mag:4` (`20261020T033119Z_POST-8-step3b-ii-e2b-mag4-n2.log`, 86 s, Status 1, `-n 2`):** `cells_ghost` [1277, 1303]; the example ran to completion; 3 calls (predicted 3); 61 of 63 probes single-incidence. (S) red as predicted: 11 byte mismatches (call 0 probe 20; call 2 probes 0, 6, 7, 8, 9, 10, 11, 13, 15, 16; call 1 none), max per-component distance **2 ULP** (range [1, 4]), `owner_orig == legacy_orig` at all 11, (S-sb) 0 in every call. Reading: probes whose legacy batch hit the point on the owning rank (owned) and on the other rank as a ghost copy (call 0 probe 20; call 2 probes 0, 6, 7, 8, 9) are A2, the gather keeping the last rank's copy; the other five (call 2 probes 10, 11, 13, 15, 16, a single rank or the last rank owning) are A1, batch composition alone; A3 not seen. Two multi-incidence probes printed (call 1 probe 15 delta 0; call 2 probe 4 owner 71716, legacy 120060, delta 1.235e-02, a record). (N0) was not printed (the assertion stops the test first).
- **`ports:1` (`20261020T033305Z_POST-8-step3b-ii-e2b-ports1-n2.log`, 58 s, Status 142; `test_complex_mode.py` 4 passed):** the mesh built (177 998 cells) and the example reached `_path_voltage` (`:561`, `:303`); the spy raised in `_in_cell` (`tests/post/test_post8_step3b_ii_example_sites.py:115`, `ValueError: could not broadcast input array from shape (3,) into shape (1,)`), computing the in-cell evaluation at the first multi-incidence point. `_ncomp` returns `index_map_bs` (1 for an N1curl function) where `f.eval` returns 3 components. A module defect, not an (S) reading; no (S), (S-sb), `cells_ghost`, (N0) or record reading exists for `ports:1`. Rank 1 then aborted in the `Allreduce` (buffer sizes 1 vs 3). Orphan check: `pgrep -c python3` 0.
- **`ports:12` not run** (an unpredicted red outside its range ends the slot's windows). Rule (i): fixing `_in_cell` (use `f.function_space.value_size` / the eval width) edits the module 551 committed, so 551's two windows and 554's run again; that is the next review's call, with the ghosted (S) form the weekly ruled.

*Written by the 2026-10-19 18:00 daily review as a re-pose of item 548 under rule (6): a measurement that separates two readings, with the target unchanged. It sets and moves no band, tolerance, Done-when or target.*

**Read, in this order:** item 548's brief, the 18:00 ruling, item 551's brief, then this.

**Executor:** `example-runner`. **Tier:** heavy, `-n 2` (`ports:12` `-n 4`), container 400. **Class:** lane 1. **Cost:** ≈ 8.8 measured / 35 predicted min. **Slot-first only.** **Depends on item 551 being DONE on `main`** (it runs 551's committed module). If 551 is not DONE, skip at no attempt cost. **Family:** e2, attempt 3 (or attempt 1 after a reset, if 551 landed its gate).

**Code:** the three script edits from the parked branch by path, and nothing else:
`git checkout attempt/POST-8-3b-ii-e2-20261019T183900Z -- examples/magnetostatics/04_helmholtz_analytic_comparison.py examples/ports/01_two_torus_port_pair.py examples/ports/12_birdcage_ring_quarter_turn.py`.
- `git diff --cached --stat` must show three files, each the keyword alone.
- **No module edit.** The module runs as 551 committed it (rule (i)).
- A site that 551's ghost print classed as ghosted (for example, `mag:6` moved here) is added to this item's windows with the same form; a site classed ghost-free is not in this item.

**Windows** (item 548's forms, labels `POST-8-step3b-ii-e2b-…`). Pre-checks and censuses as in 548. Then:
- `…-mag4-n2` (`::test_mag4`, real);
- `…-ports1-n2` (`::test_ports1`, complex, `tests/environment/test_complex_mode.py` first);
- `…-ports12-n4` (`::test_ports12`, complex, `-n 4`).

Then guide sentences for any site whose (X) is green, post-censuses, `git status`, the leak check. Measured basis: `mag:4` 100 s (`…183519Z…mag4-n2.log`), `ports:1` 133 s, `ports:12` 113 s at `-n 4` (548's table).

**Anchor:**
- (X) each example's own assertions at imported bands, including `ports:1`'s four records inside 1e-6 and `ports:12`'s C16 identity;
- (E) the predicted counts (`mag:4` 3, `ports:1` 8, `ports:12` 2);
- the ghost reading printed per rank (this item is the ghosted class; a site that reads 0 on every rank runs under 551's unchanged (S) and its predicted-red line does not apply);
- (S) as committed, expected red on bytes.

**PREDICTED RED: (S) mag4 range [1, 4]**
**PREDICTED RED: (S) ports1 range [1, 4]**
**PREDICTED RED: (S) ports12 range [1, 4]**

The range is the **maximum per-component ULP distance** over that site's single-incidence byte mismatches, as `_check_e2_s` prints it. A red is in range only when **every** mismatch also reads `owner_orig == legacy_orig` (hypothesis B refuted probe by probe), and every other assertion of the site is green. A site whose (S) comes back green is fine: a prediction not met in that direction is reported, never a red.

**The separating measurement (printed, per red probe):** the owning rank, the ranks the legacy batch hit (owned / ghost), and (S-sb), the same-batch re-evaluation mismatch count. Reading:
- (S-sb) 0 with the legacy value from the owning rank ⇒ **A1**, batch composition alone;
- the legacy value from a ghost copy on another rank ⇒ **A2**;
- (S-sb) > 0 ⇒ **A3** or something new; report it.

**Negative control:** (N0) on `ports:12`'s first `e_complex` call. It is *predicted* > 1.01 (548's brief; no measurement on this fixture), so it is printed beside the bar, not asserted (rule (e)). `mag:4` / `ports:1` print (N0) only.

**Traps:** item 548's and 551's lists. `ports:1` prints full-precision `(1*)` record candidates that passed the leak guard on 09-27; the spy adds nothing above 6 significant digits.

**Scope.** Not claimed: (S) on ghosted meshes (the item hands it to the weekly), any band or record, width invariance, `POST-8` ✅.

**Negative result.**
- **Any (S) red outside its range, or with `owner_orig ≠ legacy_orig`:** an unpredicted red. Stop, park, BLOCKED, and report.
- **An (X) red:** 548's disposition. `ports:1`'s 1e-6 records are never re-recorded in the slot.
- **Predicted reds in range, all else green:** the item **parks as predicted**. Logs and their `test-results.md` rows land on `main`, the three script edits go to `attempt/POST-8-3b-ii-e2b-*`, and item 554 is marked `BLOCKED (predicted red in range)` with the per-site A1 / A2 / A3 reading. The slot continues. The weekly rules on (S)'s form for ghosted meshes, and the scripts land with that ruling.

**Status it can move.**
- The `POST-8` §7 state line gains the three sites' (X) under the rule (or an (X) red) and the A1 / A2 / A3 reading.
- Known-issues 2026-09-27's "Verified at" gains e2b for (X).
- The weekly gains the measurement it needs to rule (S)'s form, the last open question before 3b-ii closes.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 554 (step 3b-ii-e2b) per `docs/planning/chunks/POST-8.md`: read item 548's brief, the 2026-10-19 18:00 ruling, item 551's brief and item 554's brief, in that order; 554 binds.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Take the three script edits from the parked branch by path, and edit no other file before the windows.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 400`. `pgrep -c python3` must read 0 before each.
- The three `PREDICTED RED:` lines are binding: an (S) red is in range only if every mismatch has `owner_orig == legacy_orig` and a max ULP distance in [1, 4], with every other assertion green. Run all three sites unless a red falls outside its range.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- Return the log names with elapsed times, per site `cells_ghost`, the call count / single-incidence count, every red probe with its diagnostics, (S-sb), (X) and `ports:1`'s four record misses, (N0), and the census before and after."

## §7 row moved verbatim 2026-10-19 (rotate_plan_archive.py chunks)

 Pointwise samples of a degree-1 N1curl field depend on the MPI partition: `ans:2`'s Row 1 moves −19 % … +112 % in SAR between `-n 2` and `-n 4`, while every … | 🟡 Steps 1–2 done (2026-09-27, 10-02): the partition dependence is the N1curl → CG1 interpolation, and the 09-30 weekly's owner rule (lowest `original_cell_index`) is implemented opt-in and gated at `-n 1 / 2 / 4` (16 raw SARs within 3.03e-10 of step 1's; cross-width 2.98e-10). **Step 3a done 2026-10-03** (§9 item 188): `ans:2`'s Row 1 as raw, owner-cell and 1 g averages is width-invariant at `-n 1 / 2 / 4` to ≤ 2.21e-10 (1e-9 asserted); cell/raw spans 0.64–1.60 and 1 g/raw 1.00–2.99 at the off-leg points. Step 3b classified 10-17; 3b-ii t1 and t2 done; 12 examples open, and 3b-iii (`ans:2` skip-private knob, 1 g Row 1) open. **Step 3b-i done** 2026-10-18 (§9 item 526): `point_sar` defaults to the owner rule (`owner=None` keeps the legacy path); `phantom_fields.py:163` owner-compliant by construction, asserted on 384 centroids at `-n 1 / 2 / 4`; width-invariant to 4.087e-16 (1e-9 asserted); five gate modules re-run green (`20261018T125346Z_POST-8-step3b-W5a-regate.log`, `20261018T125455Z_POST-8-step3b-W5b-regate.log`). **Step 3b-ii-t1 done** 2026-10-18 (§9 item 536): the four time-harmonic test-site calls (`test_lossy_sphere_degree2.py`, `test_lossy_sphere_fullwave.py`, `test_time_harmonic_smoke.py`, the last now sampling the DG split rather than a CG1 interpolant) use the owner rule; 25/25 sphere probes single-incidence at degree 1 and 2, where the rule is a bitwise no-op; owner indices exact and max width-spread 1.108e-11 at `-n 4 / 1` (1e-9 asserted); five gate modules green (`20261018T213254Z_POST-8-step3b-ii-t1-W1.log`). **Step 3b-ii-t2 done** 2026-10-19 (§9 item 537): the three magnetostatic DG sites (`test_coil_phantom_bfield_metrics.py` centerline / symmetry, `test_straight_wire.py::_sample_radial`) use the owner rule; `:149` / `:239` kept as the CG1 3.03× record (ruling); 91/91 probes single-incidence (coil-phantom 9 + 18 at `-n 1 / 2 / 4`, straight wire 46 + 18 at `-n 2`), the rule bytewise = legacy at every record-width probe; coil-phantom owner indices exact at `-n 1 / 2 / 4`, field spread 9.647e-04 (`-n 1`) / 5.614e-05 (`-n 4`) printed beside the 7.28e-8 LU floor (predicted ≤ 1e-6 — missed, no band); (N0) 26 cells, max/min 4.225523; both modules green (`20261019T021325Z_POST-8-step3b-ii-t2-W1.log`, `20261019T021728Z_POST-8-step3b-ii-t2-W4-wire.log`). **Step 3b-ii-e1 done** 2026-10-19 (§9 item 549, re-pose of 547): seven example sites (`mag:1` on the CG1 `B_lag`, R1 declined; `mag:5`; `th:3`, `th:4`, `th:6`, `th:7`, `th:9`) sample under the owner rule, keyword-only edits, nothing re-recorded; each runs as `__main__` under the spy with its own assertions green at `-n 2`; single-incidence 29/30, 16/16, 54/54, 160/214, 100/100, 100/100, 192/192 (rule bytewise = legacy at every single-incidence probe; `th:4`'s 54 multi-incidence are the profile points and four point-call probes, in-cell equal); owner indices and `n_inc` exact at `-n 4` / `-n 1` for six examples, **`th:4` at `-n 2` only** (known-issues 2026-10-19); (N0) `th:3` 1.000082, `mag:5` 6.329709 (`-n 2`); `mat:2` / `mat:3` / `mri:2` re-gates Status 0 (`mri:2`'s 10 g kernel-mass row prints 0.0041 % against the guide's 0.0044 %, already so on 10-18, unrelated to the rule). Logs `20261019T170143Z_POST-8-step3b-ii-e1b-mag1-n2.log` … `20261019T170830Z_POST-8-step3b-ii-e1b-mri2-regate.log`. **3b-ii-e2a done** 2026-10-20 (§9 item 551): `mag:2`, `mag:6` (CG1, keyword only) on the owner rule, ghost-free (`cells_ghost` 0 on both ranks), keyword-only edits, nothing re-recorded; each runs as `__main__` under the spy with its own assertions green at `-n 2`; 1 call each, single-incidence 25/25 and 10/10, (S) bitwise (0 byte mismatches), (S-sb) 0 (printed); (N0) `mag:2` 1.300640 (`-n 2`); `mag:4` / `ports:1` / `ports:12` on ghosted meshes in item 554. Logs `20261020T003156Z_POST-8-step3b-ii-e2a-mag2-n2.log`, `20261020T003421Z_POST-8-step3b-ii-e2a-mag6-n2.log`. **Done-when:** *(amended 2026-10-02 18:00 review, rule (j), as the 10:30 review directed; the weekly ratifies)* a sampling rule whose output is width-invariant to 1e-9 relative, asserted at `-n 1`, `-n 2` and `-n 4` on `ans:2`'s points; `ans:2`'s Row 1 quoted as a cell or mass average (3a gates both, the weekly picks the published one), width-invariant to 1e-9 relative at the same three widths; the rule adopted by every point-sampling site of the census as the weekly classifies it; and the tracked Row 1 re-issued as the averaged row ⇒ ✅. *History: `docs/planning/chunks/POST-8.md`.* 

### §9 item 555 — 2026-10-20 (lane 1, POST-8 3b-ii-e2c, re-pose of item 554: the ruled (S) form and the `_ncomp` fix)

*Drafted read-only for the 2026-10-20 03:00 daily review as a re-pose of item 554 under rule (6): same target, next attempt of the family, with the (S) form the 2026-10-19 weekly ruled (`PROJECT_PLAN.md:4694–4705`, epitaph `:4782–4785`). It sets and moves no band, tolerance, Done-when or target.*

**Read, in this order:** item 548's brief, the 2026-10-19 18:00 ruling, item 551's brief, item 554's brief and outcome, the weekly's "(S) on ghosted meshes — ruled" bullet (`PROJECT_PLAN.md:4694–4705`), then this. This brief's deltas bind.

**Executor:** `example-runner`. **Tier:** heavy (container `timeout -k 30 400` per window), `-n 2`; `ports:12` at `-n 4`, its record width. **Class:** lane 1. **Cost:** ≈ 8.8 measured (lane-1 median) / ≈ 50 predicted slot-minutes (fixed ≈ 5, checkout and code ≈ 12, windows ≈ 15, guides ≈ 7, commit, journal and branch deletes ≈ 8, margin 3). **Slot-first only.** **Depends on:** nothing (551 DONE on `main`; 554 BLOCKED, this is its re-pose). **Family:** e2, attempt 2 of 4 since 551's gate (554 was attempt 1).

**Why.** 554 ran `mag:4` with (S) red inside its predicted range ([1, 4] ULP: 11 mismatches, max 2 ULP, `owner_orig == legacy_orig` at all 11, (S-sb) 0 on all three calls; `20261020T033119Z_POST-8-step3b-ii-e2b-mag4-n2.log`). `ports:1` then crashed in the spy, not in the rule: `_in_cell` (`tests/post/test_post8_step3b_ii_example_sites.py:104–118`) sizes its buffer with `_ncomp` (`:100–101`), which returns `dofmap.index_map_bs` (1 for an N1curl space), while `f.eval` returns 3 components (`20261020T033305Z_POST-8-step3b-ii-e2b-ports1-n2.log:686`). The weekly then ruled (S): **(S-sb), the same-batch re-evaluation, is asserted; the cross-batch count (S) is printed** with each mismatch's ULP distance and rank diagnostics.

**Code.**
1. **Scripts from the parked branch, by path, three files only:**
   `git checkout attempt/POST-8-3b-ii-e2b-20261020T033900Z -- examples/magnetostatics/04_helmholtz_analytic_comparison.py examples/ports/01_two_torus_port_pair.py examples/ports/12_birdcage_ring_quarter_turn.py`
   - `git diff --cached --stat` shows those three files only. Each diff is the `owner="lowest_original_cell"` keyword alone (`mag:4` `:200`, `ports:1` `:303`, `ports:12` the two `e_complex` calls). The same three files on `attempt/POST-8-3b-ii-e2-20261019T183900Z` are byte-identical (drafter's `git diff` between the branches: empty).
   - **Never** take `docs/testing/test-results.md`, a log, `PROJECT_PLAN.md`, `known-issues.md` or a chunk file from either branch (rubric trap: an append-only record checked out by path deletes `main`'s later rows).
2. **One module edit, `tests/post/test_post8_step3b_ii_example_sites.py`, in three parts:**
   - **(a) `_ncomp` by value size.** Replace the body of `_ncomp` (`:100–101`) with the evaluator's own expression, `src/fem_em_solver/post/evaluation.py:369–370`: `shape = f.function_space.element.value_shape; return int(np.prod(shape)) if len(shape) > 0 else 1`. That is the code path the spy wraps, and it has run on this image in every e1 / e2a / e2b window (including `ports:1`'s own calls up to the crash). So the API is proven on the container's DOLFINx; do not use `index_map_bs`, `reference_value_shape` or `ufl_element()`. On every blocked Lagrange / DG vector space the e1 and e2a sites sample, `value_size == index_map_bs == 3`, so the change is a no-op there. It matters only for N1curl sources (`ports:1`, `ports:12`, and later `th:4`'s cross-check in item 556).
   - **(b) The ruled (S) form, in `_check_e2_s` only** (`:433–518`). `_check_e_s`, `_drive` and every e1 test stay byte-identical.
     - **(S-sb), asserted.** Keep the same-batch re-evaluation as written (`:474–482`: on each rank, the rank-local legacy cells `hits[i][0]` evaluated on exactly the owner batch's points and order, `idx = flatnonzero(owner_local >= 0)`, which is the batch `evaluation.py:389–403` evaluates). Count byte mismatches **at the single-incidence rows of that batch** (`c["n_inc"][idx] == 1`) and allreduce SUM. Collect `(call, count)` for every call with a non-zero count, and after the loop assert that none were collected. Also print the multi-incidence rows' same-batch count. It is printed only: at a multi-incidence point the rank-local `links[0]` can be a cell other than the owner, so a difference there is a different cell, not a different batch. Those rows stay gated by the existing in-cell check, `rel <= CLOSED_FORM_TOL` (`:509–510`). The restriction to single-incidence rows is the ruling's own reason ("owner = legacy where one cell is hit"); **the review confirms this reading when it queues the item.**
     - **(S), cross-batch, printed.** Keep the per-mismatch `S-RED` print (`:466–473`), relabelled `(S, cross-batch, printed)`, for **every** mismatch, as the ruling says. Add a per-call and a per-site summary: the mismatch count, the max |ULP| over components, and the split by hit pattern. A2-pattern means some rank other than the owning rank has a legacy hit (owned or ghost). A1-pattern means only the owning rank hit. Both are read from the `per_rank` / `owning` values already gathered. **Delete `assert not failures` (`:518`).** That is the ruled retirement of the cross-batch byte assertion (epitaph `PROJECT_PLAN.md:4782–4785`), not a loosening.
     - Update the module docstring's (S) sentence (`:17–19`) and `_check_e2_s`'s docstring to the ruled form, citing the weekly.
   - **(c) `ports:12`'s (N0) label (rule (e)).** Move `"ports12"` from `N0_CALL` (`:384`) to `N0_PRINTED_E2` (`:385`). There is no prior (N0) measurement on that fixture, and item 554's brief already labelled it predicted and printed. `_check_n0` reads `N0_CALL.get(tag, 0)`, so the call it samples is unchanged.
   - Nothing else. `git diff` on the module shows `_ncomp`'s body, `_check_e2_s`, the two docstrings and the `N0_*` lines only.

**Rule (i).** The module changes, so every window this item reports runs on the committed-to-be module:
- 551's ghost-free sites (`mag:2`, `mag:6`) re-run beside `mag:4`, `ports:1` and `ports:12`.
- **The e1 tests re-run at `-n 2`.** `_ncomp` is shared with e1's spy (`_in_cell`), so a cheap re-gate of the seven e1 tests at their record width is scheduled last (measured 4–27 s each in item 549). Their `-n 4` / `-n 1` (O) windows are not re-run: (a) is a no-op on their spaces. *If the review reads rule (i) as binding only the step being closed, it may strip these two windows; the drafter keeps them because they cost ≈ 2 min.*

**Windows.** Bare, repo-relative, foreground, Bash `timeout: 660000`, pytest `-s`. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before each compute window.
0. **Pre-checks:** `git status --porcelain` is empty; `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up; `pgrep` reads 0.
1. **Pre-censuses:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-pre-figures "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py'"`
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
2. **Dry run (prints lines only):**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-dryrun "./scripts/run_examples.sh -e 2,4,6,ports:1,ports:12 -n 2 -t 400 --dry-run"`
3. **Code:** the checkout and the module edit above.
4. **Collects:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-W0-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py --collect-only -q -s'"`
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-W0-complex "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py --collect-only -q -s'"`
5. **Ghost-free sites, real build:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-mag2-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py::test_mag2 -v -s --tb=short'"`
   then the same command with `::test_mag6`, label `POST-8-step3b-ii-e2c-mag6-n2`.
6. **`mag:4`, real build:** the same command with `::test_mag4`, label `POST-8-step3b-ii-e2c-mag4-n2`.
7. **`ports:1`, complex build:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-ports1-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_ports1 -v -s --tb=short'"`
8. **`ports:12`, complex build, `-n 4`:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-ports12-n4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 4 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_ports12 -v -s --tb=short'"`
9. **e1 re-gate, `-n 2` (rule (i), last):**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-e1-real-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_POST8_EXAMPLES=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py::test_mag1 tests/post/test_post8_step3b_ii_example_sites.py::test_mag5 -v -s --tb=short'"`
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2c-e1-complex-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_th3 tests/post/test_post8_step3b_ii_example_sites.py::test_th4 tests/post/test_post8_step3b_ii_example_sites.py::test_th6 tests/post/test_post8_step3b_ii_example_sites.py::test_th7 tests/post/test_post8_step3b_ii_example_sites.py::test_th9 -v -s --tb=short'"`
10. **Guides, then post-censuses, then `git status`.**
    - One sentence each in `mag:4`'s, `ports:1`'s and `ports:12`'s guides. Each names the owner rule, the call / single-incidence counts, (S-sb) 0, the cross-batch (S) count with its max ULP, and the log by full filename.
    - In `02_circular_loop.md:170–176` and `06_h_convergence_rate.md:172–178`, re-point the e2a sentence's log to the e2c log.
    - Post-censuses: window 1 relabelled `…-post-figures` / `…-post-docrefs`. Expect the counts unchanged, `broken=0`, `dead=0`.
    - `git status --porcelain examples/ tests/ src/ scripts/ docker/` lists only the three scripts, five guides and the module. The leak check runs on the staged logs.

Every pytest window reads N passed, 0 skipped, where N is the node count it names, plus 4 for `test_complex_mode.py` in complex windows.

**Measured basis.**

| Window | Seconds | Source |
|---|---|---|
| `mag:2` | 139 | `20261020T003156Z_POST-8-step3b-ii-e2a-mag2-n2.log` |
| `mag:6` | 146 | `20261020T003421Z_POST-8-step3b-ii-e2a-mag6-n2.log` |
| `mag:4` | 86 | `20261020T033119Z_POST-8-step3b-ii-e2b-mag4-n2.log` (ran to completion, then red) |
| `ports:1` | 133 bare | `20260927T124920Z_EX-57-control.log` |
| `ports:12` | 113 bare, `-n 4` | `20260928T112745Z_EX-57.log` |
| e1 | 4–27 each | item 549's `-n 2` logs |

`ports:1`'s spy overhead is unmeasured: one `_in_cell` Allreduce per multi-incidence probe over 8 calls, and 554 reached the first such probe at 58 s. Predicted ≤ 300 s, inside 400. Predicted total ≈ 15 min of windows.

**Park (rule (d)):**
- the `mag:2` window not started by minute 20;
- the `ports:1` window not started by minute 32;
- the `ports:12` window not started by minute 38.

If the e1 re-gates have not started by minute 44, skip them. The outcome line then says "e1 not re-run on the edited module" (rule (i)'s disclosure branch). That is not a park.

**Anchors (asserted):**
- **(X)** each example's own assertions at its imported bands under the spied `__main__` run:
  - `mag:4`: its own asserts;
  - `ports:1`: its four 0.11-image records inside `REPRODUCTION_BAND_RELATIVE` / `REPRODUCTION_BAND_SYMMETRY_ABSOLUTE` (`examples/ports/01_two_torus_port_pair.py:191–192`, asserted `:739`), and its controls `:750–771`;
  - `ports:12`: the C16 identity `:466` and its control `:494`;
  - `mag:6`: `:272` / `:286`;
  - `mag:2`: its own asserts.
- **(E)** every recorded call carries the rule, and the forwarded bytes = `own`. Call counts are printed against the prediction: `mag:2` 1, `mag:4` 3 (measured 3 in 554), `mag:6` 1, `ports:1` 8, `ports:12` 2.
- **(S-sb) = 0** at the single-incidence rows of each rank's owner batch, every call, every site (the weekly's ruled form).
- `own` is bitwise equal across ranks (`_bitwise_equal_across_ranks`).
- In-cell equality at multi-incidence probes ≤ `CLOSED_FORM_TOL` (imported, `tests/post/test_point_owner_rule.py`).

**Printed, with predictions:**
- `cells_ghost`: `mag:2` / `mag:6` [0, 0]; `mag:4` [1277, 1303] (554); `ports:1` / `ports:12` > 0 (`two_torus_domain` / `birdcage_port_domain` use `GhostMode.shared_facet`, `src/fem_em_solver/io/mesh.py:2306–2326`, `:4911–4920`).
- **Cross-batch (S):** `mag:2` / `mag:6` 0; `mag:4` 11 mismatches, max 2 ULP, A2-pattern 6 / A1-pattern 5, as 554 read them; `ports:1` / `ports:12` unknown, max ULP predicted ≤ 4 by `mag:4`'s analogy.
- (S-sb) at the multi-incidence rows.
- `ports:1`'s four record misses.

**Negative control (rule (e)):**
- *Asserted:* (N0) on `mag:2`'s `B` > 1 + `N0_SPREAD_FLOOR`. Measured 1.300640 on the same fixture and comparison twice (`20261019T183230Z_POST-8-step3b-ii-e2-mag2-n2.log`; `20261020T003156Z_…e2a-mag2-n2.log`).
- *Printed only:* (N0) for `mag:4`, `ports:1` and `ports:12` (predicted > 1.01, no prior measurement on those fixtures); `mag:6` skipped (CG1).
- *Comparator sensitivity, printed:* the cross-batch (S) on `mag:4` should re-find 554's 11 mismatches (the same comparison on the same fixture and partition). That shows the byte comparator sees 1–2 ULP. A 0 there is reported, not a red.

**No `PREDICTED RED:` line.** (S-sb) ≠ 0 anywhere is an unpredicted red: stop.

**Traps already paid for:**
- item 548's, 551's and 554's lists;
- **batch composition** (item 537, rubric element 4): bytes are compared only within one batch, which is exactly what (S-sb) does;
- **the path checkout of an append-only record** (code step 1);
- **`-k a or b` splits into stray argv** inside the quoted container command: name pytest node IDs, as window 9 does;
- **a kill is not a result.** After a kill, run `pgrep -c python3` first, then at 0 the rubric's cache sweep as `POST-8-step3b-ii-e2c-cache-sweep`. Never raise a timeout or `-n`;
- **the bash hook** denies `pgrep -f "python3 -m pytest"`: `pgrep -c python3` is the check;
- never set `FEM_EM_SETUP_FIGURES`; `./scripts/run_examples.sh` only with `--dry-run`;
- **`ports:1`'s 1e-6 records are now a live risk.** 554's crash proves at least one multi-incidence arc point, where owner ≠ legacy would shift V by about a quadrature weight times the jump. The (R) rows print owner vs legacy per multi-incidence probe, so a miss is attributable;
- **leak guard:** `ports:1` prints full-precision `(1*)` record candidates that passed on 09-27, and the spy prints at ≤ 6 significant digits (ULP distances are integers). On a flag, move the log unchanged to `logs/leak-held/`, mark BLOCKED and report (item 526 precedent);
- no AED number anywhere; never open `docs/private/` or `aed_results/`.

**Scope.** Five e2 sites at record width, the module's three-part edit, five guide sentences, and the e1 `-n 2` re-gate. **Not claimed:**
- width invariance on the e2 sites (no (O), on cost, as 548);
- `th:4` / `mag:1` / `mri:1` (item 556);
- 3b-iii;
- any band or record;
- `POST-8` ✅.

**Negative result.**
- **(S-sb) ≠ 0 at a single-incidence row:** a new finding against the rule's same-batch claim. Journal the per-row census, open a known-issues entry, park on `attempt/POST-8-3b-ii-e2c-*`, mark 555 BLOCKED in the same commit, and stop.
- **An (X) red:** 548's disposition. Run the bare `-n <record width>` example from `git stash` at container 400 as the baseline (label `…-<tag>-baseline`).
  - A green baseline means the rule moved the record: journal the multi-incidence (R) rows.
  - `ports:1`'s 1e-6 records are **never** re-recorded in the slot; that is the review's (`(1*)`, `record-reconciler`) under ruling.
- **An e1 re-gate red:** the `_ncomp` edit is not a no-op on that space. Journal the space's `index_map_bs` / value size, park, and mark BLOCKED.
- **A kill:** not a result; park with the measured time.

**Status it can move.**
- The `POST-8` §7 state line replaces "the ghosted-mesh sites (item 554, BLOCKED …)" with: "**3b-ii-e2 done** <date> (§9 item 555, e2c): `mag:2`, `mag:4`, `mag:6`, `ports:1`, `ports:12` on the owner rule at record width, each example's assertions green (`ports:1`'s four records inside 1e-6); (S-sb) = 0 asserted (same batch, weekly 2026-10-19); cross-batch (S) printed — `mag:4` <n> mismatches ≤ <u> ULP (A1 <a> / A2 <b>), `ports:1` <…>, `ports:12` <…>; e3 and 3b-iii open." The row stays 🟡.
- Known-issues 2026-09-27's "Verified at" gains e2c. The entry stays OPEN until 3b-iii.
- Item 556 becomes runnable.

**Commit:**
- the module, the three scripts and the five guides;
- the logs, with their `test-results.md` rows appended by hand (`git diff --stat` shows only insertions on record files);
- `PROJECT_PLAN.md`: the §7 state line, and §9 item 555 DONE;
- known-issues "Verified at";
- this file's outcome line.

The journal goes in its own commit. **Then delete both parked branches,** since their content is on `main`: e2's module was superseded by 551's and this item's, and the e2 / e2b scripts are the three files landed here.
`git branch -D attempt/POST-8-3b-ii-e2b-20261020T033900Z`
`git branch -D attempt/POST-8-3b-ii-e2-20261019T183900Z`
The sandbox warns on `.git/config`; the delete still happens. On a park, both branches stay.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 555 (step 3b-ii-e2c, re-pose of 554) per `docs/planning/chunks/POST-8.md`: read item 548's brief, the 2026-10-19 18:00 ruling, item 551's brief, item 554's brief and outcome, `PROJECT_PLAN.md:4694–4705` (the weekly's (S) ruling) and item 555's brief, in that order; 555 binds.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Take the three scripts from `attempt/POST-8-3b-ii-e2b-20261020T033900Z` by path, and never path-checkout a record file.
- Make the module edit exactly as the brief's parts (a)–(c) say. `_check_e_s`, `_drive` and the e1 tests stay byte-identical.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 400` (180 for the e1 re-gates), in the brief's order. `pgrep -c python3` must read 0 before each compute window.
- `./scripts/run_examples.sh` only with `--dry-run`. Never set `FEM_EM_SETUP_FIGURES`. Re-record nothing.
- There is no `PREDICTED RED:` line. A non-zero single-incidence (S-sb) anywhere stops the item.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e2c-*`, mark item 555 BLOCKED in the same commit, and report.
- Return:
  - the log names with elapsed times;
  - per site: `cells_ghost`, the call count / single-incidence count, (S-sb) single and multi, the cross-batch (S) count with max ULP and the A1 / A2 split, (N0), (X) and `ports:1`'s four record misses;
  - the e1 re-gate result;
  - the census before and after."

**Outcome (2026-10-20 04:30 slot): BLOCKED, `ports:1` (X) red; parked on `attempt/POST-8-3b-ii-e2c-20261020T095000Z`.** The `_ncomp` fix worked: `ports:1` ran past the first multi-incidence point and its example asserts ran. Raw record 0.8945163788 missed by 1.333e-04 against the 1e-06 band (`20261020T094001Z_POST-8-step3b-ii-e2c-ports1-n2.log:704`); the bare `-n 2` example from `git stash` (no keyword) reproduces it at 0.00e+00 (`20261020T094255Z_POST-8-step3b-ii-e2c-ports1-baseline.log:650`), so the owner rule moved the record. The spy's (R) rows were not printed because the example asserts before returning; nothing re-recorded. Green: `mag:2` (1 call, 0 ULP), `mag:6` (1 call), `mag:4` (3 calls, cells_ghost [1277, 1303], cross-batch (S) 8 mismatches ≤ 2 ULP, A1 3 / A2 5, (S-sb) 0), `ports:12` (-n 4, 2 calls, cells_ghost [1149, 1054, 1376, 1269], cross-batch (S) 6 mismatches ≤ 32 ULP, all A2, (S-sb) 0), e1 re-gates 2 + 9 passed. Item 556 stays blocked on 555. Re-pose needs the review's ruling on `ports:1`'s records (`(1*)`, `record-reconciler`) before the example edit lands.

### §9 item 556 — 2026-10-20 (lane 1, POST-8 3b-ii-e3: `th:4` and `mag:1` re-sampled from the source field, `mri:1` keyword; the weekly's ruled step)

*Dependency re-pointed by the 2026-10-20 10:30 review (ruling (5), below): 556 depends on item 568 Part A (the `_ncomp` fix, landed with the module), not on 555 — if 568 is not DONE on `main`, skip at no attempt cost; where this brief says "555's module" or "555's `_ncomp` fix", read 568's.*

*Drafted read-only for the 2026-10-20 03:00 daily review from the 2026-10-19 weekly's ruling (`PROJECT_PLAN.md:4706–4717`, "`th:4`, `mag:1` R1, `mri:1` `:486–487` — ruled, one step"). The content change and the `(1*)` re-record are the weekly's licence. This brief supplies execution only and sets and moves no band, tolerance, Done-when or target.*

**Read, in this order:** item 547's brief (R1 block, `:757–761`, `:830–835`), the 03:00 and 10:30 rulings of 2026-10-19, item 549's brief, known-issues 2026-10-19 (`th:4`), item 555's brief, then this.

**Executor:** `example-runner`. **Tier:** standard (container `timeout -k 30 180`), `-n 1 / 2 / 4`. **Class:** lane 1. **Cost:** ≈ 8.8 measured (lane-1 median) / ≈ 45 predicted slot-minutes:
- fixed ≈ 5;
- code ≈ 12;
- windows ≈ 6;
- guides with two table re-records ≈ 12;
- known-issues retirement ≈ 3;
- commit and journal ≈ 6;
- margin ≈ 1.

**Slot-first only.**

**Depends on item 555 DONE on `main`.** If 555 is not DONE, skip at no attempt cost. There are two reasons:
- **Shared code.** `th:4`'s re-sampled cross-check evaluates the raw N1curl `fields.e_complex`. Its axial line has a multi-incidence probe (probe 20, n_inc 2, item 547's outcome), and the spy's `_in_cell` on an N1curl space crashes exactly as `ports:1` did until 555's `_ncomp` fix lands.
- **The same module.** This item adds entries to the module 555 edits. Running it on 555's committed module keeps both items' windows on what is committed (rule (i)).

**Family:** `POST-8` 3b-ii-e3; attempt 1.

**Code.** Six edits, nothing else.
1. **`th:4` (`examples/time_harmonic/04_evanescent_waveguide_decay.py`).**
   - **Profile:** `:384` samples the DG split, `_probe_ey(fields.e_real, fields.e_imag, profile_points, comm)`, instead of the CG1 export `e_re` / `e_im`. `_probe_ey` already carries the owner rule (`:227–232`).
   - **Cross-check:** `:379–381` refits γ from the raw N1curl source, `fields.e_complex`, at `points` under the owner rule. Add a helper `_probe_ey_complex(field, points, comm)`: one `evaluate_vector_field_parallel(field, points, comm, owner="lowest_original_cell")`, raise on any invalid point as `_probe_ey` does, and return `np.real(v[:, 1]), np.imag(v[:, 1])`. The assertion stays `cg_vs_nedelec < CG1_VS_NEDELEC_MAX` (`:145`, `:408`), and `PROFILE_RMS_MAX` (`:146`, `:414`) is unmoved.
   - **CG1 export:** `_paraview_fields` (`:245–269`) and the XDMF stay as they are. The CG1 refit and CG1 profile are still computed from `e_re` / `e_im`, but **printed only**, labelled "(CG1 export — partition-dependent, printed, not asserted; `-n 2` record in the guide)".
   - Update the comment block `:135–144` and the print at `:394–406` to say which field each figure reads.
2. **`mag:1` (`examples/magnetostatics/01_straight_wire.py`).** At `:273–275`, sample the DG source `B` (from `solver.compute_b_field()`, `:245`) under the owner rule instead of `B_lag`. Update the comment at `:247–249`. `B_lag` stays the VTX export, and its round-trip check (`VTX_ROUNDTRIP_RTOL`, `:31`, `:129`) is untouched.
3. **`mri:1` (`examples/mri/01_coil_phantom_fields.py:486–487`).** Add `owner="lowest_original_cell"` to both calls. That is the keyword only (the ruling).
4. **Module `tests/post/test_post8_step3b_ii_example_sites.py`, additive only.**
   - Set `PREDICTED_CALLS["th4"]` to 9 (printed, not asserted; per run: decay 2, N1curl cross-check 1, DG profile 2, CG1 refit 2 and CG1 profile 2, the last four printed only).
   - Add `"mri1": (REPO / "examples/mri/01_coil_phantom_fields.py", 2)` to `EXAMPLES` and `"mri1": 2` to `PREDICTED_CALLS`.
   - Add `@complex_only def test_mri1(): _drive_e2("mri1")`. It goes through the e2 checker, so the ruled (S-sb) form applies whatever `mri:1`'s mesh class is.
   - No existing function body changes.

**Windows.** Bare, repo-relative, foreground, Bash `timeout: 660000`, pytest `-s`, `pgrep -c python3` 0 before each compute window. **The `-n 2` windows run before the `-n 4` / `-n 1` ones,** because they rewrite `logs/post8/step3b-ii-e1-{th4,mag1}-n2.npz`, which (O) reads.
0. **Pre-checks** as item 555's.
1. **Pre-censuses:** item 555's two commands relabelled `POST-8-step3b-ii-e3-pre-figures` / `…-pre-docrefs`.
2. **Dry run:**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e3-dryrun "./scripts/run_examples.sh -e 1,th:4,mri:1 -n 2 -t 180 --dry-run"`
3. **Code.**
4. **Collects:** item 555's W0 commands relabelled `POST-8-step3b-ii-e3-W0-real` / `…-W0-complex`.
5. **`-n 2`, the runs of record.**
   - Real:
     `scripts/testing/run_and_log.sh POST-8-step3b-ii-e3-mag1-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_POST8_EXAMPLES=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py::test_mag1 -v -s --tb=short'"`
   - Complex:
     `scripts/testing/run_and_log.sh POST-8-step3b-ii-e3-th4-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_th4 -v -s --tb=short'"`
   - Then the same complex command with `::test_mri1`, label `POST-8-step3b-ii-e3-mri1-n2`.
6. **`-n 4`:** the `mag1` command with `mpiexec -n 4` (label `…-mag1-n4`), and the `th4` command with `mpiexec -n 4` (label `…-th4-n4`).
7. **`-n 1`:** the same two commands with `mpiexec -n 1` (labels `…-mag1-n1`, `…-th4-n1`).
8. **Guides, known-issues, post-censuses, `git status`.**
   - **`th:4` guide (`04_evanescent_waveguide_decay.md`).**
     - Re-record the table rows `:47–48` and the text at `:129–140`, `:157` and `:164` from the `-n 2` log under the `(1*)` licence. Each row names its field (N1curl cross-check, DG-split profile) and the `-n 1 / 2 / 4` logs.
     - Keep the CG1 figures as a **labelled `-n 2` record**: 0.124 % / 0.027 % on 2026-10-19 (item 549's `-n 2` log, `…170233Z_…th4-n2.log`), and August's 0.200 % / 0.117 % (`20260810T020325Z_EX-7-run1.log.gz:53`).
     - Rewrite 549's sentence at `:170–176`, which says the band is a `-n 2` record: the re-sampled profile is gated at `-n 1 / 2 / 4`.
   - **`mag:1` guide (`01_straight_wire.md`).**
     - Re-record the three sampled rows (relative L2, max relative, decay ratio) at `:17–23` from the `-n 2` log. Energy and mesh are not sampled and do not move.
     - Keep the CG1 `B_lag` digits (51.9781 % / 76.7330 % / 20.31) as a labelled `-n 2` record beside them. `EX-30`'s HTML-comment pattern (`:25–36`) is acceptable.
     - Re-sync `:107–114` ("reproduce … at fixed rank count") to say the DG sample is width-invariant up to the printed (O) spread.
   - **`mri:1` guide:** one sentence naming the owner rule, the call / single-incidence counts, (S-sb) and the `-n 2` log. Its centerline record row (`01_coil_phantom_fields.md:107`) is **not** re-recorded here. If a printed digit moves, report it with the (R) rows; the review rules.
   - **Known-issues:** retire 2026-10-19 (`th:4`) in the house format (`- ✅ RETIRED <date> (POST-8 step 3b-ii-e3, §9 item 556, <slot>) — was 🟡 OPEN 2026-10-19 …`, as at `docs/testing/known-issues.md:878`). Its resolve clause, "`th:4` green at `-n 1 / 2 / 4` on the re-sampled profile", is met by windows 5–7. Add e3 to 2026-09-27's "Verified at", which stays OPEN until 3b-iii.
   - Post-censuses relabelled. `git status --porcelain examples/ tests/ src/ scripts/ docker/` shows only the three scripts, three guides and the module; check that `mri:1`'s runtime outputs are gitignored. Run the leak check.

**Measured basis:**
- item 549's `-n 2` windows: `mag:1` 8 s, `th:4` 9 s;
- its `-n 4` / `-n 1` groups: mag 10 / 36 s, th 31 / 50 s for four examples;
- `mri:1` 6 s bare at `-n 2` (`01_coil_phantom_fields.md:162`).

Predicted ≈ 6 min of windows, including container start.

**Park (rule (d)):** the first `-n 2` window not started by minute 22, or the `-n 1` windows not started by minute 35.

**Anchors (asserted):**
- **(X)** each example's own assertions under the spied `__main__` run:
  - `th:4`'s `PROFILE_RMS_MAX` and `CG1_VS_NEDELEC_MAX` (bands unmoved) on the re-sampled fields **at `-n 1`, `-n 2` and `-n 4`**, together with `GAMMA_RTOL`, `FIELD_L2_MAX` and `IMAG_RATIO_MAX`;
  - `mag:1`'s VTX round-trip;
  - `mri:1`'s own asserts.
- **(E) / (S)** as item 549 for `th:4` / `mag:1`, through `_drive` → `_check_e_s`. Both meshes are ghost-free (the 18:00 ruling, point (5)), so the same-batch premise holds there.
- **(S-sb) = 0** for `mri:1`, through `_drive_e2`.
- In-cell equality at every multi-incidence probe ≤ `CLOSED_FORM_TOL`. This covers `th:4`'s 25 profile points: they are all n_inc 2 on the z = L/2 lattice plane. `E_y` is tangential there, so N1curl continuity makes both incident cells agree up to round-off.
- **(O)** at `-n 4` / `-n 1` for `th4` and `mag1`: the call count, line numbers, `n_inc` and owner original indices exactly equal the `-n 2` npz.

**Predicted, printed (never asserted):**
- `th:4`'s re-sampled profile RMS somewhere in [0.01 %, 0.5 %], width spread ≤ 1e-6 relative;
- the N1curl-vs-DG cross-check ≤ 1e-12 (the DG1 split of a degree-1 N1curl field is exact, `src/fem_em_solver/core/time_harmonic.py:745–754`, and `th:4` solves at degree 1, `:198`);
- the (O) field spread: `th:4` ≤ 1e-6, `mag:1` ≤ 1e-4 (the class basis is `mag:5`'s DG spread, 3.3e-06 / 2.2e-05 in item 549, against the CG1 `B_lag`'s 1.681 / 0.6965);
- `mag:1`'s three re-recorded rows, direction unknown, since probe 0 is n_inc 24 at the wire edge;
- `mri:1`: 2 calls, single-incidence count unknown.

**Flag for the queueing review (rule (h)).** Re-sampled from the source as ruled, `th:4`'s cross-check compares the N1curl fit with its own exact DG1 split. It is therefore an identity at round-off, a comparand true by construction. The assertion stays as the weekly ruled, with `CG1_VS_NEDELEC_MAX` unmoved. Its informative content (the export path) now lives in the CG1 record line, so the journal should call it a record-grade check.

**Negative control (rule (e)):**
- *Printed, measured basis:* the CG1 profile at `-n 4`, printed beside the re-sampled one in the same window, is predicted at 2.249 %. The same comparison on the same fixture read 2.249 % in `20261019T140617Z_POST-8-step3b-ii-e1-th4-baseline-n4.log:50`, against 0.124 % at `-n 2`. That is the partition-dependent route this step retires, so the separation is printed from both routes in one window.
- *Not asserted:* it is not wired as an assertion in the example. It would add an asserted quantity the weekly did not rule.
- The module's asserted (N0) is not exercised here (`th:3` / `mag:5` are not re-run).

**Traps already paid for:**
- item 549's and item 555's lists, in particular batch composition (bytes only within a batch), node IDs instead of `-k`, never `FEM_EM_SETUP_FIGURES`, and dry-run only;
- **window order**: `-n 2` before `-n 4` / `-n 1`, or (O) reads a stale npz with the old 6-call structure and goes red on the call count;
- **`th:4` at `-n 4` is no longer a known red** once the profile reads the DG split. If a CG1 printed line reads ≈ 2.2 % there, that is the expected control, not a red;
- **`mag:1`'s export path is untouched**: switching `B_lag` in the VTX writer would void the round-trip control;
- **a kill is not a result** (pgrep, then the cache sweep as `POST-8-step3b-ii-e3-cache-sweep`);
- leak guard as 555; no AED number.

**Scope.** `th:4`'s profile and cross-check, `mag:1`'s R1 table, `mri:1`'s keyword, three guides, and one known-issues retirement. **Not claimed:**
- the CG1 export's own width dependence (step 1's mechanism, known-issues 2026-09-27, not reopened — the ruling);
- `mri:1`'s centerline record re-recorded;
- 3b-iii;
- any band;
- `POST-8` ✅.

**Negative result.**
- **`th:4` red at any width on the re-sampled profile or cross-check:** an unpredicted red. Stop. Known-issues 2026-10-19 stays OPEN with the reading appended. Park on `attempt/POST-8-3b-ii-e3-*` and mark 556 BLOCKED in the same commit.
- **(O) red** (owner indices or n_inc differ across widths): a new finding against the rule, with the same disposition.
- **(S) / (S-sb) red:** stop, as in items 549 / 555.
- **`mag:1` re-recorded rows** are records and cannot go red.
- There is no `PREDICTED RED:` line.

**Status it can move.**
- The `POST-8` §7 state line gains: "**3b-ii-e3 done** <date> (§9 item 556): `th:4`'s profile (DG split) and cross-check (N1curl) and `mag:1`'s R1 table (DG `B`) re-sampled under the owner rule and re-recorded `(1*)`, bands unmoved; `th:4` green at `-n 1 / 2 / 4` (profile <x> %); (O) exact for both; `mri:1` keyword, (S-sb) 0; CG1 figures kept as labelled `-n 2` records."
- With 555 DONE: "**3b-ii complete**; 3b-iii open."
- **Known-issues 2026-10-19 retires.** 2026-09-27's "Verified at" gains e3.

**Commit:**
- the three scripts, three guides and the module;
- the logs, with their `test-results.md` rows appended by hand;
- known-issues (retirement plus "Verified at");
- `PROJECT_PLAN.md`: the §7 state line, and §9 item 556 DONE;
- this file's outcome line.

The journal goes in its own commit.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 556 (step 3b-ii-e3) per `docs/planning/chunks/POST-8.md`: read item 547's brief (R1 block), the 2026-10-19 03:00 and 10:30 rulings, item 549's brief, known-issues 2026-10-19, `PROJECT_PLAN.md:4706–4717` (the weekly's ruling), item 555's brief and item 556's brief, in that order; 556 binds. If item 555 is not DONE on `main`, stop at no attempt cost.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Edit only the six things the brief names. The module change is additive.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 180`. Run the `-n 2` windows before the `-n 4` / `-n 1` ones. `pgrep -c python3` must read 0 before each compute window.
- `./scripts/run_examples.sh` only with `--dry-run`. Never set `FEM_EM_SETUP_FIGURES`. Re-record only the `th:4` and `mag:1` guide rows the brief names, keeping the CG1 figures as labelled `-n 2` records.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e3-*`, mark item 556 BLOCKED in the same commit, and report.
- Return:
  - the log names with elapsed times;
  - per example and width: the call count / single-incidence count, (O) and its spread;
  - `th:4`'s profile RMS and cross-check at `-n 1 / 2 / 4`, with the CG1 printed lines beside them;
  - `mag:1`'s new rows against the CG1 record;
  - `mri:1`'s (S-sb) and any moved digit;
  - the census before and after."

### Review ruling 2026-10-20 03:00 (items 555, 556; 3b-iii)

Binds over the briefs above. (1) **Item 554 diagnosed from code and log, no pathologist needed:** `_ncomp` (`tests/post/test_post8_step3b_ii_example_sites.py:100–101`) returns `dofmap.index_map_bs`, 1 on N1curl, against `eval`'s three components (`20261020T033305Z_POST-8-step3b-ii-e2b-ports1-n2.log:686`). A module defect, not an (S) reading. `mag:4` landed its pre-registered red in range with (S-sb) 0 on all three calls. (2) **(S-sb) is asserted at the single-incidence rows of each rank's owner batch** (item 555's reading, confirmed): the weekly's reason — "owner = legacy where one cell is hit" — is a single-incidence claim; multi-incidence rows print their same-batch count and stay gated by the in-cell `CLOSED_FORM_TOL` check. (3) The e1 re-gate at `-n 2` stays in 555 (rule (i): `_ncomp` is shared). (4) `ports:12`'s (N0) moves to printed (rule (e)). (5) **556 depends on 555** (shared `_in_cell`, same module); both are slot-first. (6) **3b-iii is not queued — "write the anchor", returned to the weekly:** `ans:2` SPEC Row 1 is pointwise and Row 4 is already the 1 g average (`examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/SPEC.md:166–168`); re-issuing Row 1 as the 1 g average is a SPEC content change, `metrics.json["pointwise"]` has a consumer (`tests/validation/test_ans9_far_side_sar.py:220–228`), and the Done-when amendment ("the 1 g mass average", `POST-8.md:346`) is owed with it. Family count e2: 555 is attempt 2 of 4 since 551's gate.

### Review ruling 2026-10-20 10:30 (item 555 blocked — diagnosis and re-pose)

(1) Diagnosis confirmed from the logs: a module-clean red — the rule, not the spy, moved `ports:1`'s record (`20261020T094001Z_POST-8-step3b-ii-e2c-ports1-n2.log:704` against `20261020T094255Z_POST-8-step3b-ii-e2c-ports1-baseline.log:650`). No (S-sb) finding.
(2) Re-pose under rule (6), family e2 attempt 3 of 4 (554 = 1, 555 = 2): item 568 (`POST-8` 3b-ii-e2d) SPLITS the step. Part A lands what 555 measured green: the module edit (a)–(c) and the `mag:4` and `ports:12` keyword edits, taken by path from the e2c branch (never a record file), with every window re-run on the module as committed (rule (i)): `mag:2`, `mag:6`, `mag:4`, `ports:12` (-n 4), the e1 -n 2 re-gates. Part B MEASURES `ports:1` under the owner rule without changing what `examples/ports/01_two_torus_port_pair.py` commits on main (it stays on the legacy path, no keyword committed): the four record values under the owner rule and their misses against the four 0.11-image records, at -n 1, -n 2 and -n 4 (the width-invariance POST-8's Done-when asks for, 1e-9 relative), plus the (R) rows per multi-incidence probe (owner vs legacy cell, value delta). Find the least-invasive way to get those prints: e.g. a measurement-only test in the e2 module that runs the parked keyword version and catches the example's AssertionError at its record check, or that calls the example's compute function before its assert — read the example to decide; no edit to the example on main beyond nothing. Nothing is re-recorded in the slot.
(3) Whether `ports:1`'s four records are re-recorded on the owner rule under the (1*) licence (record-reconciler, bands unmoved), or `ports:1` stays on the legacy path with a disclosed exemption, is RETURNED TO THE 2026-10-21 WEEKLY with Part B's measurement — the weekly ruled the analogous th:4 / mag:1 content change itself on 10-19.
(4) If Part A lands, it is a landed gate for the e2 family ((S-sb) asserted on four sites, the §7 state line moves) and the count resets; e2 itself closes only when `ports:1` lands.
(5) Item 556 (e3) now depends on 568 Part A (the `_ncomp` fix) instead of 555: "if 568 is not DONE on main, skip".
(6) The drafter's execution choices in item 568's brief are confirmed:
- Part B injects the keyword through `_make_spy(inject_owner=...)`, so the example on `main` is untouched.
- The width spread of the four `ports:1` values is **printed** beside the Done-when's 1e-9, not asserted. Asserting it would be a band decision on a record the weekly has not ruled. That the owner cells are exact across widths is asserted.
- `test_ports1`'s skip names ruling (3) and goes when the weekly rules.
- Part A commits first, and 568 is DONE on Part A. A Part B red after that commit is journaled for the weekly; it does not park.

(7) **Correction to the 04:30 journal** (attempts.md 2026-10-20T09:52Z). The journal says the other three record misses were not printed, but they were (`20261020T094001Z_POST-8-step3b-ii-e2c-ports1-n2.log:678`, re-read):

| Record | Miss | Band |
|---|---|---|
| corrected | 1.31e-04 | 1e-06 |
| symmetry | 3.01e-06 (abs) | 1e-06 |
| ‖S‖₂ | 3.77e-05 | 1e-06 |

All four records miss. What the run lost was the spy's output.

### §9 item 568 — 2026-10-20 (lane 1, POST-8 3b-ii-e2d, re-pose of item 555: land the four green sites, measure ports:1)

**Outcome 2026-10-20 (16:30 CDT slot, attempt 3): done.** Part A landed at `8ddbc9de`; Part B measured at `-n 2 / 4 / 1`. Part A: `mag:2` (1 call, 25/25 single-incidence, 0 mismatches, (N0) 1.300640), `mag:6` (1 call, 10/10), `mag:4` (3 calls, 61/63, cells_ghost [1277, 1303], cross-batch (S) 10 mismatches ≤ 2 ULP, A1 3 / A2 7, (N0) 1.577795), `ports:12` (`-n 4`, 2 calls, 126/128, cells_ghost [1149, 1054, 1376, 1269], cross-batch (S) 6 mismatches ≤ 32 ULP, all A2, (N0) 2.453210); (S-sb) 0 everywhere; e1 re-gates 2 + 9 passed. Part B (`…214423Z…ports1-measure-n2.log`, `…214712Z…-n4.log`, `…214952Z…-n1.log`): 5 passed at each width; 8 calls, 23826/24584 single-incidence, (S-sb) 0, cross-batch (S) 0, (O) exact at `-n 4` / `-n 1`. Values at `-n 2`: raw 0.8946356119, corrected 0.9399444943, symmetry 3.412669537e-05, ‖S‖₂ 0.8613244573 — misses 1.333e-04 / 1.308e-04 / 3.005e-06 abs / 3.766e-05 against the 1e-06 bands, identical at all three widths and equal to item 555's `…094001Z…ports1-n2.log:677–678`, so the injection equals the parked keyword. Pre-empted gates green (superseded 3.943e-04 / 3.869e-04, corrected −6.0056 %, blind −98.2573 %, symmetry 3.4127e-05, ‖S‖₂ 0.861324). Cross-width spread against `-n 2` (printed): `-n 1` raw 2.460e-10, corrected 2.415e-10, ‖S‖₂ 2.208e-10, symmetry 1.019e-06 relative (3.479e-11 absolute); `-n 4` raw **1.134e-09**, corrected **1.113e-09**, ‖S‖₂ 5.346e-10, symmetry 1.526e-07 relative (5.207e-12 absolute) — the `-n 4` raw / corrected are past the 1e-9 figure (predicted ≤ 1e-9, missed; not asserted). (R) rows: 758 multi-incidence probes per width; owner ≠ legacy cell at 370 (`-n 2`), 320 (`-n 4`), 342 (`-n 1`), deltas up to 0.37 (equal-cell rows delta 0); the owner indices are identical across widths, the legacy indices are not. Censuses: figures `ok=133 broken=0` before and after; docrefs `dead=0 guide=0`, exit 2, stale 117 → 118 (info). Nothing re-recorded; `ports:1`'s example and guide untouched; branches `attempt/POST-8-3b-ii-e2b-20261020T033900Z` and `attempt/POST-8-3b-ii-e2-20261019T183900Z` deleted, `attempt/POST-8-3b-ii-e2c-20261020T095000Z` kept for the weekly.

*Drafted read-only for the 2026-10-20 10:30 daily review as the re-pose of item 555 under rule (6) that the 10:30 ruling above orders. It sets and moves no band, tolerance, record, Done-when or target.*

**Read, in this order:** item 548's brief, the 2026-10-19 18:00 ruling, item 551's brief, item 555's brief and outcome, the 2026-10-20 03:00 and 10:30 rulings, then this. This brief's deltas bind. Where it says "as 555", 555's text applies with every label `…-e2c-…` read as `…-e2d-…`.

**Executor:** `example-runner`. **Tier:** heavy (container `timeout -k 30 400` per window; 180 for the e1 re-gates), `-n 2`; `ports:12` at `-n 4`, its record width; Part B's `ports:1` at `-n 2`, `-n 4` and `-n 1`. **Class:** lane 1. **Cost:** ≈ 8.8 measured (lane-1 median) / ≈ 50 predicted slot-minutes (fixed ≈ 5, checkout and code ≈ 10, windows ≈ 20, guides ≈ 4, Part A commit and branch deletes ≈ 5, Part B commit and journal ≈ 4, margin 2). The basis for the overheads is 555 itself: ten windows, then park and journal, inside ≈ 22 min wall (first log `20261020T093112Z`, journal 09:52Z, attempts.md:1572–1585). **Slot-first only.** **Depends on:** nothing (551 DONE on `main`; 555 BLOCKED, this is its re-pose). **Family:** e2, attempt 3 of 4 since 551's gate (554 = 1, 555 = 2). A landed Part A resets the count (ruling (4)).

**Why.**
- 555's module edit (a)–(c) and the `mag:4` / `ports:12` keyword edits measured green with (S-sb) 0 everywhere (attempts.md:1575–1579). Part A lands exactly that.
- `ports:1` under the owner rule missed **all four** 0.11-image records. The example prints them before it asserts (`…094001Z_…-ports1-n2.log:678`): raw 1.33e-04, corrected 1.31e-04 (relative), symmetry 3.01e-06 (absolute), ‖S‖₂ 3.77e-05 (relative), each against its 1e-06 band. The bare legacy example reproduces all four (`…094255Z_…-ports1-baseline.log:650`). *Correction to attempts.md:1582:* the three other misses were printed (`:678`); what the run lost was the spy's output — (E), (S-sb), the (R) rows and (N0) — because `_run_spied` re-raises the example's AssertionError before `_check_e2_s` runs.
- The weekly needs the (R) rows and the width behaviour to rule (ruling (3)). Part B obtains them without committing the keyword to the example.

**Mechanism chosen for Part B (ruling (2) "read the example to decide").** The example has no compute function to call before its assert: `main()` (`examples/ports/01_two_torus_port_pair.py:398–806`) solves, computes, prints the four misses (`:700–708`) and asserts (`:728–743`) in one body. Its only point-sampling call is `_path_voltage`'s (`:303`); the e2c branch's whole `ports:1` diff is `owner="lowest_original_cell"` on that call. So the measurement test **injects that keyword through the spy** on the example's own calls (main's file runs unedited, and the forwarded call is byte-for-byte the parked keyword call), catches the example's record AssertionError by its message, and reads the four values from two pass-through recorders on the package functions the example imports by name (`:90–95`). Running the parked script from a copied path was rejected: the spy keys on the caller's file path, and the copy would need its own import shims.

**Code.**
1. **Three files by path from the e2c branch:**
   `git checkout attempt/POST-8-3b-ii-e2c-20261020T095000Z -- tests/post/test_post8_step3b_ii_example_sites.py examples/magnetostatics/04_helmholtz_analytic_comparison.py examples/ports/12_birdcage_ring_quarter_turn.py`
   - `git diff --cached --stat` shows those three files only. The module diff is 555's (a)–(c): `_ncomp`'s body, `_check_e2_s`, the two docstrings and the `N0_*` lines (drafter's `git diff main...attempt/POST-8-3b-ii-e2c-20261020T095000Z`). Each script diff is the keyword alone (`mag:4` `:200`; `ports:12` the two `e_complex` calls, `:503–507`).
   - **Never** `examples/ports/01_two_torus_port_pair.py`: it stays on the legacy path on `main`.
   - **Never** a record file. The e2c branch also carries `docs/testing/test-results.md` and 13 logs that `main` already has (`386f45e9`). A path checkout of the record deletes `main`'s later rows.
2. **Part B additions to the module, (d1)–(d4).** All are made before any compute window, so every window runs on the module as committed (rule (i)).
   - **(d1) `_make_spy(path, calls, inject_owner=None)`.** Move the caller test ahead of the forward: compute `is_target` first. If `inject_owner is not None and is_target and "owner" not in kwargs`, set `kwargs = {**kwargs, "owner": inject_owner}`. Then `fwd = real(function, points, *args, **kwargs)`, then `if not is_target: return fwd`, and the rest unchanged. The record's `"owner"` then reads the injected value, so `_check_e2_s`'s (E) holds as written. With `inject_owner=None` the behaviour is unchanged. Part A's windows, the e1 re-gates included, run on it.
   - **(d2) `_run_spied_measure(tag)`, new.** It is `_run_spied`'s body with three changes:
     - it builds the spy as `_make_spy(path, calls, inject_owner=RULE)`;
     - it adds two pass-through recorders on the same `MonkeyPatch`, on `fem_em_solver.ports.systematics.mutual_systematics_ladder` and `fem_em_solver.ports.sparameters.sparameters_from_impedance`. Import both modules inside the function. Each recorder appends `(args, kwargs, return)` and returns the real result unchanged. The example binds both names at `run_path` time, so it sees the patch;
     - around `runpy.run_path`: `except AssertionError as exc:` keeps `str(exc)` only if it contains `"does not reproduce the 0.11-image record"` (the record check, `:739–743`), and otherwise re-raises. `SystemExit` is handled as in `_run_spied`.

     It returns `(calls, ladder_calls, s_calls, record_msg)`.
   - **(d3) `test_ports1_owner_measure`, new, `@complex_only`, no width assertion.**
     - `_run_spied_measure("ports1")`; then `_check_e2_s("ports1", calls, comm)`, giving (E), (S-sb), own bitwise across ranks, the in-cell check, the (R) rows and the cross-batch (S); then `_check_n0("ports1", calls, comm, asserted=False)`; then `_check_o("ports1m", calls, comm)`. At `-n 2` that writes `logs/post8/step3b-ii-e1-ports1m-n2.npz`; at `-n 1` / `-n 4` it compares against it. The `e1` prefix is `_check_o`'s, reused unchanged.
     - **Constants are imported, never restated:** `K = runpy.run_path(str(EXAMPLES["ports1"][0]), run_name="ports1_constants")`. That runs the top level only, because `main()` sits under the `__main__` guard (`:809`).
     - **The four values.** `ladder` is the first ladder call's return (the measured `im_z12`; the second call is the blind fixture, `:613–614`). `S` is the first S call's return, and `Z` its first argument. Then raw = `ladder["raw"]` and corrected = `ladder["corrected"]`; symmetry and spectral use the example's own expressions (`:617–620`).
     - **Printed:** the four values (≤ 10 significant digits), their misses against `K["RECORDED_*"]` in the example's metric (symmetry absolute), the bands, and `record_msg` (or "records reproduced" if nothing was raised).
     - **Asserted on the captured values with `K`'s constants:** the example's gates that its record assert pre-empted (`:749–774`). These are the superseded-v0.7.2 control (both misses > `REPRODUCTION_BAND_RELATIVE`), `abs(ladder["deviation"]) < MUTUAL_TOLERANCE`, the blind fixture's `>= MUTUAL_TOLERANCE`, symmetry < `S_SYMMETRY_BAND`, and spectral ≤ `S_SPECTRAL_NORM_CEILING`.
     - **Width.** At `-n 2`, rank 0 writes `logs/post8/step3b-ii-e2d-ports1-records-n2.npz` holding the four values and `Z`. At `-n 1` / `-n 4` the test reads it (`pytest.fail` if it is missing; never recompute) and prints each value's relative cross-width spread, plus symmetry's absolute spread, beside the 1e-9 Done-when figure.
   - **(d4) `test_ports1` gains a skip** with the reason "ports:1 stays on the legacy path pending the 2026-10-21 weekly (POST-8.md, review ruling 2026-10-20 10:30 (3)); measured by test_ports1_owner_measure". Without the keyword on `main`, its (E) assertion (`c["owner"] == RULE`) fails by construction. The body is unchanged. *The review confirms this when it queues the item.*
   - Nothing else. Against `main`, `git diff` on the module shows (a)–(c), `_make_spy`'s head, two new functions, one new test, one decorator, and one docstring sentence naming Part B's test.

**Windows.** As 555: bare, repo-relative, foreground, Bash `timeout: 660000`, pytest `-s`, and `pgrep -c python3` reading 0 before each compute window.
0. **Pre-checks** as 555.
1. **Pre-censuses:** 555's two commands, labelled `POST-8-step3b-ii-e2d-pre-figures` / `…-pre-docrefs`.
2. **Dry run:** `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2d-dryrun "./scripts/run_examples.sh -e 2,4,6,ports:1,ports:12 -n 2 -t 400 --dry-run"`
3. **Code:** steps 1 and 2 above.
4. **Collects:** 555's W0 commands, labelled `POST-8-step3b-ii-e2d-W0-real` / `…-W0-complex`. The complex collect lists `test_ports1_owner_measure`.

*Part A.*

5. **Ghost-free sites, then `mag:4`, real build:** 555's windows 5–6, labelled `POST-8-step3b-ii-e2d-mag2-n2`, `…-mag6-n2`, `…-mag4-n2`.
6. **`ports:12`, complex build, `-n 4`:** 555's window 8, labelled `POST-8-step3b-ii-e2d-ports12-n4`.
7. **e1 re-gate, `-n 2`:** 555's window 9, labelled `POST-8-step3b-ii-e2d-e1-real-n2` / `…-e1-complex-n2`.
8. **Guides, then post-censuses, then `git status`.**
   - One sentence each in `mag:4`'s and `ports:12`'s guides. Each names the owner rule, the call / single-incidence counts, (S-sb) 0, the cross-batch (S) count with its max ULP, and the log by full filename.
   - In `02_circular_loop.md:170–176` and `06_h_convergence_rate.md:172–178`, re-point the e2a log to the e2d log.
   - `ports:1`'s guide is not touched.
   - Post-censuses: window 1 relabelled `…-post-figures` / `…-post-docrefs`. Expect the counts unchanged, `broken=0`, `dead=0`.
   - `git status --porcelain examples/ tests/ src/ scripts/ docker/` lists only the two scripts, four guides and the module. Run the leak check on the staged logs.
9. **Part A commit (below), then the two branch deletes.**

*Part B.* Complex build, the measurement test. **`-n 2` first:** it writes the two npz files the other widths read.

10. `scripts/testing/run_and_log.sh POST-8-step3b-ii-e2d-ports1-measure-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_ports1_owner_measure -v -s --tb=short'"`
11. The same command with `mpiexec -n 4`, labelled `POST-8-step3b-ii-e2d-ports1-measure-n4`.
12. The same command with `mpiexec -n 1`, labelled `POST-8-step3b-ii-e2d-ports1-measure-n1`.
13. **Part B commit, then the journal.**

Every pytest window reads N passed, 0 skipped, where N is the node count it names, plus 4 for `test_complex_mode.py` in complex windows. Part B reads 5 passed.

**Measured basis.**

| Window | Seconds | Source |
|---|---|---|
| `mag:2` `-n 2` | 148 | `20261020T093309Z_POST-8-step3b-ii-e2c-mag2-n2.log` |
| `mag:6` `-n 2` | 150 | `20261020T093547Z_POST-8-step3b-ii-e2c-mag6-n2.log` |
| `mag:4` `-n 2` | 91 | `20261020T093824Z_POST-8-step3b-ii-e2c-mag4-n2.log` |
| `ports:12` `-n 4` | 128 | `20261020T094529Z_POST-8-step3b-ii-e2c-ports12-n4.log` |
| e1 real / complex `-n 2` | 10 / 40 | `20261020T094743Z_…-e1-real-n2.log` / `20261020T094800Z_…-e1-complex-n2.log` |
| `ports:1` `-n 2`, spied, to its record assert | 154 | `20261020T094001Z_POST-8-step3b-ii-e2c-ports1-n2.log` |
| `ports:1` `-n 2`, bare | 140 (mesh 28.7 s, solves 26.1 + 25.9 s) | `20261020T094255Z_…-ports1-baseline.log:656` |
| censuses, dry run, collects | 0–3 each | 555's logs (attempts.md:1584) |

Part A ≈ 9.7 min of windows.

Part B is predicted at ≈ 8.5 min. The `-n 2` run is ≈ 160 s (154 s, plus < 10 s for the recorders, the top-level `run_path` and the npz). The `-n 4` run is ≈ 130 s. The `-n 1` run is ≈ 230 s: the solves roughly double, while gmsh is serial either way. All three are inside 400.

**Park (rule (d)):**
- the `mag:2` window not started by minute 15;
- the `ports:12` window not started by minute 27;
- the Part A commit not made by minute 38.

Any of these parks the whole item, with nothing landed.

**Part B cut-offs (disclosed, not parks):**
- `…-measure-n2` not started by minute 42: Part B is not run. 568 is marked DONE on Part A with "Part B not run", and the review re-poses it.
- `-n 4` not started by minute 46, or `-n 1` by minute 48: that width is skipped. The outcome names the widths measured.

**Anchors (asserted).**
- **Part A,** as 555 for the four sites:
  - **(X)** each example's own assertions at its imported bands under the spied `__main__` run;
  - **(E)** every recorded call carries the rule, and forwarded bytes = `own`. Counts are printed against the prediction: `mag:2` 1, `mag:4` 3, `mag:6` 1, `ports:12` 2;
  - **(S-sb) = 0** at the single-incidence rows of each rank's owner batch, every call, every site;
  - `own` bitwise equal across ranks;
  - in-cell ≤ `CLOSED_FORM_TOL` at multi-incidence probes;
  - the e1 re-gates green (2 passed real, 9 passed complex, as 555).
- **Part B,** at each width, on the owner rule:
  - **(E)** on all calls, with the injected rule and forwarded bytes = `own`. The count is printed against the prediction of 8;
  - **(S-sb) = 0**, `own` bitwise across ranks, and in-cell ≤ `CLOSED_FORM_TOL`;
  - **the example's pre-empted gates**: superseded control, mutual band, blind control, symmetry band, passivity. All are predicted green from e2c's own `-n 2` prints: corrected −6.01 % inside 10 % (`…094001Z_…-ports1-n2.log:671`), symmetry 3.4127e-05 < 1e-3 and ‖S‖₂ 0.861324 ≤ 1 (`:676`), superseded 3.94e-04 / 3.87e-04 > 1e-06 (`:679`);
  - the only failure tolerated inside the run is the record check, matched by its message. Any other AssertionError is red;
  - **(O)** at `-n 4` / `-n 1`: the call count, call lines, `n_inc` and owner original indices exactly equal the `-n 2` npz. This is the rule's own width claim, exact at all seven e1 sites (item 549).

**Printed, with predictions.**
- **Part A:**
  - `cells_ghost`: `mag:2` / `mag:6` [0, 0]; `mag:4` [1277, 1303]; `ports:12` [1149, 1054, 1376, 1269] (555).
  - Cross-batch (S): `mag:2` / `mag:6` 0. `mag:4` non-zero and ≤ 2 ULP; the count is not reproducible to the unit (11 in 554, 8 in 555). `ports:12` 6, ≤ 32 ULP, all A2 (555 predicted ≤ 4 and read 32).
  - (N0): `mag:4` 1.577795, `ports:12` 2.453210 (555).
- **Part B:**
  - **The four values and misses at `-n 2`:** predicted to reproduce e2c's `…094001Z_…-ports1-n2.log:677–678` exactly (same image, same partition, same call). This is the check that injection equals the parked keyword. A difference is reported, and the measurement still stands.
  - **The misses at `-n 4` / `-n 1`:** predicted ≈ the `-n 2` misses.
  - **Cross-width spread against the 1e-9 figure.** Raw, corrected and spectral are predicted ≤ 1e-9 relative: the owner cells are identical by (O), and the residual is the solve's partition round-off (≤ 3.03e-10 at `ans:2`'s 3a gate). Symmetry is a 3.4e-05 residual of a near-cancellation, so its relative spread is predicted up to ≈ 3e-05; its absolute spread is predicted ≤ 1e-9.
  - **(R) rows per multi-incidence probe:** owner vs legacy original cell, and the value delta. Predicted: at least one row with owner ≠ legacy and delta ≫ 1e-12 on the undriven-port calls, which are what `ladder` reads. That row attributes the 1.33e-04 shift.
  - `_check_o`'s per-call sample spread; cross-batch (S); (N0).

**Negative control (rule (e)).**
- *Asserted:*
  - (N0) on `mag:2`'s `B` > 1 + `N0_SPREAD_FLOOR`. Measured 1.300640 three times: `20261019T183230Z_…e2-mag2-n2.log`, `20261020T003156Z_…e2a-mag2-n2.log`, `20261020T093309Z_…e2c-mag2-n2.log`.
  - Part B: the example's own superseded-v0.7.2 control on the owner-rule values, measured 3.94e-04 / 3.87e-04 under the rule at `-n 2` (`…094001Z_…:679`).
- *Printed, measured basis, not re-run:* the legacy side of Part B's separation is `…094255Z_…-ports1-baseline.log:650` (raw and corrected 0.00e+00, symmetry 1.84e-11 absolute, ‖S‖₂ 1.23e-10). The owner-rule misses are printed beside it, citing it.
- *Comparator sensitivity, printed:* `mag:4`'s cross-batch (S) non-zero, as 555.

**No `PREDICTED RED:` line.** Each of the following is an unpredicted red, and stops the item:
- (S-sb) ≠ 0 anywhere;
- a Part B (O) red;
- an AssertionError the catch does not match.

**Traps already paid for:**
- the lists in items 548, 551, 554 and 555;
- **the path checkout of an append-only record**: the e2c branch carries `test-results.md` and 13 logs, so check out the three named files only;
- **`ports:1`'s script is never taken from any branch**: Part B measures by injection;
- **the catch is narrow**: match the message substring, never a bare `except AssertionError`;
- **Part B window order**: `-n 2` first, or `_check_o` and the records comparison read a missing or stale npz;
- **`logs/post8/` is gitignored** (`.gitignore:82`): the npz files are never committed;
- a kill is not a result: `pgrep -c python3` first, then at 0 the cache sweep as `POST-8-step3b-ii-e2d-cache-sweep`. Never raise a timeout or `-n`;
- the bash hook: `pgrep -c python3` is the check;
- never set `FEM_EM_SETUP_FIGURES`; `./scripts/run_examples.sh` only with `--dry-run`;
- **leak guard**: `ports:1` prints full-precision `[EX-34]` record candidates (e2c's log passed the guard and is on `main` at `386f45e9`), and the test's own prints are ≤ 10 significant digits. On a flag, move the log unchanged to `logs/leak-held/`, mark BLOCKED and report;
- **nothing is re-recorded**: `ports:1`'s `RECORDED_*`, its guide and its bands stay untouched (ruling (3));
- no AED number anywhere; never open `docs/private/` or `aed_results/`.

**Scope.**
- Part A: `mag:2`, `mag:4`, `mag:6` and `ports:12` on the owner rule at record width; the module; two scripts; four guide sentences; the e1 `-n 2` re-gate.
- Part B: `ports:1` measured under the owner rule at `-n 1 / 2 / 4`, by injection.

**Not claimed:**
- `ports:1` on the owner rule on `main`;
- any record or band;
- width invariance on the four Part A sites (no (O), on cost, as 548);
- `th:4` / `mag:1` / `mri:1` (item 556);
- 3b-iii;
- e2 closed (ruling (4));
- `POST-8` ✅.

**Negative result.**
- **Before the Part A commit:** as 555.
  - (S-sb) ≠ 0 at a single-incidence row: a new finding, so journal the per-row census and open a known-issues entry.
  - An (X) red: 548's disposition, with the bare baseline from `git stash` labelled `…-<tag>-baseline`.
  - An e1 re-gate red: (a) or (d1) is not a no-op on that space.

  Each parks on `attempt/POST-8-3b-ii-e2d-*` with 568 BLOCKED in the same commit, and stops. All three older branches stay.
- **After the Part A commit,** Part A stays landed and 568 is marked DONE on it. Part B's result is stated in the §9 line and the outcome:
  - a Part B (S-sb) or (O) red is a new finding against the rule on this mesh: a known-issues entry and a journal entry, and it stops;
  - a pre-empted gate red under the owner rule is a measurement for the weekly: journal the value, not re-tried in the slot;
  - a harness error in the recorders, the catch or the npz is not a measurement: journal it, and the review re-poses Part B.
- **A kill:** not a result. In Part A, park with the measured time; in Part B, disclose.

**Status it can move.**
- **The `POST-8` §7 state line.** Replace "the ghosted-mesh sites (item 554 BLOCKED on a module `_ncomp` defect, `mag:4` (S-sb) 0; re-posed as item 555 under the weekly's ruled (S) form)" with: "**3b-ii-e2 four sites done** <date> (§9 item 568 Part A, e2d): `mag:2`, `mag:4`, `mag:6`, `ports:12` on the owner rule at record width, each example's assertions green; (S-sb) = 0 asserted (same batch, weekly 2026-10-19); cross-batch (S) printed — `mag:4` <n> ≤ <u> ULP, `ports:12` <n> ≤ <u> ULP; `ports:1` measured on the owner rule (Part B: four record misses <…> at `-n` <widths>, (O) <exact>, cross-width spread <…>), its records returned to the 2026-10-21 weekly; e2 closes when `ports:1` lands; e3 (item 556) and 3b-iii open." The row stays 🟡.
- Known-issues 2026-09-27's "Verified at" gains e2d (Part A). The entry stays OPEN until 3b-iii.
- The e2 family count resets (ruling (4)).
- Item 556 becomes runnable.
- The 2026-10-21 weekly receives Part B's measurement for ruling (3).

**Commits.**
- **Part A** (after window 8):
  - the module, `mag:4`'s and `ports:12`'s scripts, and the four guides;
  - the Part A logs, with their `test-results.md` rows appended by hand (`git diff --stat` shows only insertions on record files);
  - `PROJECT_PLAN.md`'s §7 state line, Part A clause;
  - known-issues "Verified at".

  Then delete the two older branches. Their `mag:4` / `ports:12` scripts are now on `main`, and their `ports:1` keyword is also on the e2c branch, which stays:
  `git branch -D attempt/POST-8-3b-ii-e2b-20261020T033900Z`
  `git branch -D attempt/POST-8-3b-ii-e2-20261019T183900Z`
  The sandbox warns on `.git/config`; the delete still happens.
- **Part B:**
  - the Part B logs with their `test-results.md` rows;
  - the §7 clause's Part B figures;
  - §9 item 568 DONE;
  - this file's outcome line.
- **The journal** goes in its own commit.
- **Keep `attempt/POST-8-3b-ii-e2c-20261020T095000Z`** until the weekly rules on `ports:1`. On a park, all three branches stay.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 568 (step 3b-ii-e2d, re-pose of 555) per `docs/planning/chunks/POST-8.md`: read item 548's brief, the 2026-10-19 18:00 ruling, item 551's brief, item 555's brief and outcome, the 2026-10-20 03:00 and 10:30 rulings and item 568's brief, in that order; 568 binds.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Take exactly three files from `attempt/POST-8-3b-ii-e2c-20261020T095000Z` by path: the module, `mag:4`'s script and `ports:12`'s script. Never take `ports:1`'s script or a record file.
- Make the module additions (d1)–(d4) exactly as the brief says, before any compute window. `_check_e_s`, `_drive`, `_check_o`, `_check_n0` and the e1 tests stay byte-identical.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 400` (180 for the e1 re-gates), in the brief's order. Commit Part A before Part B's windows. Run Part B at `-n 2` first. `pgrep -c python3` must read 0 before each compute window.
- `./scripts/run_examples.sh` only with `--dry-run`. Never set `FEM_EM_SETUP_FIGURES`. Re-record nothing; `ports:1`'s example and guide stay untouched.
- There is no `PREDICTED RED:` line. A non-zero single-incidence (S-sb), a Part B (O) red, or an unmatched AssertionError stops the item.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- Before the Part A commit, a red parks on `attempt/POST-8-3b-ii-e2d-*` with item 568 BLOCKED in the same commit. After it, Part A stays landed and Part B's result is reported.
- Return:
  - the log names with elapsed times;
  - per Part A site: `cells_ghost`, the call count / single-incidence count, (S-sb) single and multi, the cross-batch (S) count with max ULP and the A1 / A2 split, (N0) and (X);
  - the e1 re-gate result;
  - per Part B width: the four values and their misses, `record_msg`, the pre-empted gates, (S-sb), (O), the cross-width spreads, and every (R) row;
  - the census before and after;
  - the branch deletes."

### Review ruling 2026-10-20 18:00 (items 568 landed, 556 runnable)

Item 568 is DONE (`8ddbc9de` Part A, `f66f3b9e` Part B), so item 556 (3b-ii-e3) is runnable from the 19:30 slot as a slot's first item. Line drift since the brief was written: the 10-19 weekly's `th:4` / `mag:1` / `mri:1` ruling now sits at `PROJECT_PLAN.md:4744–4755` (the spawn line's `4706–4717`). Read it there; nothing else in 556's brief changes. Part A's landing resets the e-family's attempt count (a landed gate); e3 is attempt 1 of its own sub-step. `ports:1`'s four records stay with the 2026-10-21 weekly (Part B's measurement: identical misses at `-n 1 / 2 / 4`, `20261020T214423Z_…-ports1-measure-n2.log`, `…214712Z_…-n4.log`, `…214952Z_…-n1.log`); no slot re-records them. `attempt/POST-8-3b-ii-e2c-20261020T095000Z` stays until that ruling.

### Review ruling 2026-10-21 03:00 (item 556 blocked — diagnosis and re-pose)

Binds over item 556's brief. (1) **Diagnosis, from the parked log and the module, no pathologist needed (no reading is disputed):** at `-n 4`, `test_th4` failed in `_check_e_s` (`tests/post/test_post8_step3b_ii_example_sites.py:261` on `attempt/POST-8-3b-ii-e3-20261021T010145Z`), the e1 checker, which asserts bytes between the owner-batch row and the full-array legacy row — a **cross-batch** comparison. Call 0 (the decay fit), probe 24, single-incidence, differs in the last byte of the real part of two components (`20261021T010051Z_POST-8-step3b-ii-e3-th4-n4.log` on the branch, the FAILURES block). Every example assertion held at `-n 4` (profile 1.591 %, N1curl-vs-DG 0.000e+00, CG1 control 2.249 % printed as predicted); `th:4` and `mag:1` were green at `-n 1 / 2`, `mri:1` at `-n 2` with (S-sb) 0. (2) **The brief routed `th:4` through `_check_e_s` on the premise "both meshes are ghost-free" (556 brief, Anchors).** That premise was established at `-n 2` only (item 549 ran `th:4` at `-n 2`); it was never measured at `-n 4`. The 10-19 weekly retired the cross-batch byte test as an assertion ("a bitwise claim names the batch it is bitwise in", `PROJECT_PLAN.md:4822–4825`; the ruled form `:4734–4742`), and the rubric trap (item 537) says the same: whether the `-n 4` mesh is ghosted or the batch composition alone moved the ULP, the cross-batch byte equality is not a claim the rule makes. (3) **Re-pose under rule (6), family e3 attempt 2 of 4 (556 = 1): item 580.** Same target, same three sites, same code by path from the branch, with one module change: `th:4` is checked in the ruled form — `_check_e2_s` ((S-sb) asserted at single-incidence owner-batch rows, cross-batch (S) printed with ULP and the A1 / A2 split, `cells_ghost` per rank printed) — **at all three widths**, followed by `_check_o` for (O) as before. `mag:1` stays on `_drive` (green at all three widths under the stricter form; no reason to relax it). No band, tolerance, Done-when or target moves. (4) The branch's red `-n 4` log is the attempt's evidence: 580 lands it on `main` by path (a new log file, never `test-results.md`), appends its row by hand, and deletes the branch on green.

### §9 item 580 — 2026-10-21 (lane 1, POST-8 3b-ii-e3b: re-pose of item 556 — `th:4` checked in the ruled (S) form at every width)

*Written by the 2026-10-21 03:00 daily review as a re-pose of item 556 under rule (6): same target, the next attempt of the family (e3, attempt 2 of 4). It sets and moves no band, tolerance, Done-when or target.*

**Read, in this order:** item 556's brief (above), the 2026-10-21 03:00 ruling (above), item 568's brief for `_check_e2_s`'s output, then this. 556's brief binds except where this one changes it.

**Executor:** `example-runner`. **Tier:** standard (container `timeout -k 30 180`), `-n 1 / 2 / 4`. **Class:** lane 1. **Cost:** ≈ 12.1 measured (lane-1 median) / ≈ 45 predicted slot-minutes (fixed ≈ 5, code by path and one edit ≈ 6, windows ≈ 7, guides ≈ 12, known-issues ≈ 3, commit / journal / branch delete ≈ 8, margin ≈ 4). Basis: item 556 ran its nine windows and parked in ≈ 6 minutes (19:56:49 → 20:02:38). **Slot-first only.** **Depends on:** nothing (556 BLOCKED; this is its re-pose). **Family:** e3, attempt 2 of 4.

**Code.**
1. **By path from `attempt/POST-8-3b-ii-e3-20261021T010145Z` (code only, never a record file):** `examples/time_harmonic/04_evanescent_waveguide_decay.py`, `examples/magnetostatics/01_straight_wire.py`, `examples/mri/01_coil_phantom_fields.py`, `tests/post/test_post8_step3b_ii_example_sites.py`. `git diff --stat main` afterwards shows exactly these four files.
2. **One module change (disclosed):** add
   ```python
   def _drive_e3(tag):
       comm = MPI.COMM_WORLD
       calls = _run_spied(tag)
       _check_e2_s(tag, calls, comm)
       if tag in N0_CALL:
           _check_n0(tag, calls, comm, asserted=True)
       _check_o(tag, calls, comm)
   ```
   and change `test_th4`'s body from `_drive("th4")` to `_drive_e3("th4")`. Nothing else: `_drive`, `_check_e_s`, `_check_e2_s`, `_check_o` and every other test unchanged. (`_drive_e2` is not reused because it asserts the record width.)
3. **The red log by path:** `docs/testing/logs/20261021T010051Z_POST-8-step3b-ii-e3-th4-n4.log` from the branch; append its `test-results.md` row by hand (the branch's row text, verbatim). No other branch log lands — this item's own windows are the runs of record.

**Windows.** Exactly item 556's windows 0–8 (pre-checks, pre-censuses, dry run, collects, `-n 2` runs of record `mag1` / `th4` / `mri1`, then `-n 4` and `-n 1` for `mag1` and `th4`, guides, known-issues, post-censuses), with the labels `POST-8-step3b-ii-e3b-…` in place of `…-e3-…`. The **`-n 2` windows run before the `-n 4` / `-n 1` ones** ((O) reads the `-n 2` npz). Bare, repo-relative, foreground, Bash `timeout: 660000`, pytest `-s`, `pgrep -c python3` 0 before each compute window.

**Anchors (asserted):** item 556's, with one change — **`th:4`'s (S) is (S-sb) = 0 at single-incidence owner-batch rows at `-n 1`, `-n 2` and `-n 4`** (`_check_e2_s`'s final `assert not sb_failures`, `:608` in the branch module); the cross-batch count is printed. Everything else stands: (X) each example's own assertions (`PROFILE_RMS_MAX`, `CG1_VS_NEDELEC_MAX`, `GAMMA_RTOL`, `FIELD_L2_MAX`, `IMAG_RATIO_MAX` unmoved) at all three widths; (E); in-cell equality at multi-incidence probes ≤ `CLOSED_FORM_TOL`; `mag:1` through `_drive` / `_check_e_s` as 556 measured it; `mri:1` (S-sb) = 0; (O) exact for `th4` and `mag1` at `-n 4` / `-n 1`.

**Predicted, printed (never asserted):** `th:4` at `-n 4`: cross-batch (S) ≥ 1 mismatch at call 0 probe 24, 1 ULP (556's reading); `cells_ghost` per rank — **the reading 556 lacked**: nonzero on at least one rank at `-n 4` predicted, all zero at `-n 1`; profile RMS 1.591 % at `-n 4` (556's), the CG1 control ≈ 2.249 %; `th:4` (S-sb) 0 at every width.

**Negative control (rule (e)):** as 556 — the CG1 profile at `-n 4` printed beside the re-sampled one (2.249 % measured, `20261019T140617Z_POST-8-step3b-ii-e1-th4-baseline-n4.log:50`, and 556's `-n 4` window), not asserted.

**Traps already paid for:** item 556's list; **taking the branch's `test-results.md` or any record file by path** (the `8d4cf58` trap — rows by hand, `git diff --stat` on record files shows only insertions); the `-n 2`-first window order; a kill is not a result (pgrep, then the cache sweep as `POST-8-step3b-ii-e3b-cache-sweep`).

**Scope.** 556's scope, unchanged. **Not claimed:** whether `th:4`'s `-n 4` mesh is ghosted (printed); the CG1 export's own width dependence; 3b-iii; any band; `POST-8` ✅.

**Negative result.** (S-sb) red at any width on `th:4`, or an example assertion red, or (O) red: stop, park on `attempt/POST-8-3b-ii-e3b-*`, mark 580 BLOCKED in the same commit, known-issues 2026-10-19 stays OPEN with the reading appended; the e3 family is then at 2 of 4. There is no `PREDICTED RED:` line.

**Status it can move.** As 556: the `POST-8` §7 state line gains "**3b-ii-e3 done** <date> (§9 item 580, re-pose of 556): `th:4`'s profile (DG split) and cross-check (N1curl) and `mag:1`'s R1 table (DG `B`) re-sampled under the owner rule and re-recorded `(1*)`, bands unmoved; `th:4` green at `-n 1 / 2 / 4` (profile <x> %), (S-sb) 0 at every width, cross-batch (S) printed (`cells_ghost` <…> at `-n 4`); (O) exact for both; `mri:1` keyword, (S-sb) 0"; **known-issues 2026-10-19 retires**; 2026-09-27's "Verified at" gains e3. 3b-ii stays open on `ports:1` (the weekly's).

**Commit:** the three scripts, three guides, the module, the landed red log and this item's logs with their `test-results.md` rows by hand, known-issues, `PROJECT_PLAN.md` (§7 state line, §9 item 580 DONE), this file's outcome line; then `git branch -D attempt/POST-8-3b-ii-e3-20261021T010145Z`. Journal in its own commit.

**Spawn line:** "Use the example-runner agent to execute `POST-8` §9 item 580 (step 3b-ii-e3b, re-pose of item 556) per `docs/planning/chunks/POST-8.md`: read item 556's brief, the 2026-10-21 03:00 review ruling, item 568's brief (the `_check_e2_s` output) and item 580's brief, in that order; 580 binds over 556.
- You are the executor: do not spawn agents.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Take the four code files from `attempt/POST-8-3b-ii-e3-20261021T010145Z` by path, plus the one red log; never a record file. Make only the `_drive_e3` / `test_th4` edit the brief names.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 180`. Run the `-n 2` windows before the `-n 4` / `-n 1` ones. `pgrep -c python3` must read 0 before each compute window.
- `./scripts/run_examples.sh` only with `--dry-run`. Never set `FEM_EM_SETUP_FIGURES`.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e3b-*`, mark item 580 BLOCKED in the same commit, and report.
- Return: the log names with elapsed times; per example and width the call count / single-incidence count, (S-sb), the cross-batch (S) count with max ULP, `cells_ghost` per rank, (O) and its spread; `th:4`'s profile RMS and cross-check at `-n 1 / 2 / 4` with the CG1 printed lines beside them; `mag:1`'s new rows against the CG1 record; `mri:1`'s (S-sb); the census before and after; the branch delete."

**Outcome 2026-10-21 (06:00 CDT slot, item 580, e3 attempt 2): done, green first try.** Code by path from `attempt/POST-8-3b-ii-e3-20261021T010145Z` plus the module's `_drive_e3` (`test_th4` now checks (S-sb) via `_check_e2_s` at every width). `th:4` green at `-n 1 / 2 / 4` (9 calls, 200/305 single-incidence; profile 1.591 % at all widths, cross-check 1.887e-16 / 0.000e+00 / 1.887e-16, CG1 printed lines 0.124 % / 2.249 % / 0.124 %); (S-sb) 0 at every width; cross-batch (S) 0 / 0 / 4 mismatches at `-n 1 / 2 / 4` (max 1024 ULP at `-n 4`, all A2); `cells_ghost` `[0]` / `[578, 576]` / `[559, 652, 520, 670]` (the nonzero-at-`-n 4` prediction held, and also at `-n 2`); (O) exact on call count, lines, `n_inc`, owner indices. Printed misses: (O) field spread for `th:4` 1.691e-01 (`-n 4`) / 5.083e-02 (`-n 1`) against the predicted ≤ 1e-6, taken over all nine calls including the CG1 calls; profile RMS 1.591 % against the brief's predicted [0.01, 0.5] %; `mag:1` spread 1.527e-05 / 1.434e-04 against ≤ 1e-4 at `-n 1`. `mag:1` new rows 43.8071 % / 84.1711 % / 44.36 (CG1 record 51.9781 % / 76.7330 % / 20.31); `mri:1` 2 calls 10/10 single-incidence, (S-sb) 0, centerline `z = 0` `|B|` 6.261265e-07 against the recorded 6.261281e-07 (not re-recorded; the review rules). Census: figures broken 0; docrefs 134 / 134 / 0 / 0, stale 112, unchanged. Known-issues 2026-10-19 retired. Logs `20261021T110135Z_POST-8-step3b-ii-e3b-mag1-n2.log` (8 s), `…110150Z_…th4-n2.log` (9 s), `…110204Z_…mri1-n2.log` (5 s), `…110214Z_…mag1-n4.log` (8 s), `…110229Z_…th4-n4.log` (7 s), `…110242Z_…mag1-n1.log` (15 s), `…110303Z_…th4-n1.log` (11 s); the first `W0-real` call (`…110115Z`) was a transient docker-socket denial and was re-run at once (`…110123Z`). Landed red log `20261021T010051Z_POST-8-step3b-ii-e3-th4-n4.log`.

### Weekly ruling 2026-10-21 — `ports:1` (3b-ii-e4) and 3b-iii (full text PROJECT_PLAN §10)

- **`ports:1`:** a `(1*)` re-record (3b-ii-e4, lane 1). Part B measured all four 0.11-image records missing identically at `-n 1 / 2 / 4` with (S-sb) 0 and (O) exact: the move is the two-torus path-voltage samples, not the solve, and it is width-free. The four records are re-recorded from the owner-rule values at `-n 2` with the 1e-6 bands unmoved, the legacy values kept in the guide as a labelled record, the assertions green at all three widths; the cross-width spread of the S-derived records (1.1e-9 at `-n 4`, the LU floor) is printed, never asserted at the 1e-9 point band. e2 closes with it; the e2c attempt branch is deleted on landing. A legacy-path exemption is declined.
- **3b-iii, the anchor:** the published Row 1 is the **owner-cell average** (Row 4 is already the 1 g average). The step: `FEM_EM_ANS_SKIP_PRIVATE` on `ans:2` (the `ans:4` pattern); Row 1 re-issued as the owner-cell values, width-invariant to 1e-9 at `-n 1 / 2 / 4` (3a's band), the CG1 values kept as a labelled record; `SPEC.md` Row 1 re-defined as "SAR in the owner cell at the probe point" (a v-bump, the pointwise band unchanged, disclosed); `metrics.json["pointwise"]` keeps its key with a `"pointwise_rule": "owner_cell"` field beside it and its consumer `test_ans9_far_side_sar.py:220–228` is re-run (green, or `(1*)` re-recorded if it compares a stored number); known-issues 2026-09-27 retires on the width-invariance. The §7 Done-when now reads "owner-cell average". `POST-8` → ✅ needs e4 and 3b-iii both.

## §7 row moved verbatim 2026-10-21 (rotate_plan_archive.py chunks)

 Pointwise samples of a degree-1 N1curl field depend on the MPI partition: `ans:2`'s Row 1 moves −19 % … +112 % in SAR between `-n 2` and `-n 4`, while every … | 🟡 Steps 1–3a done; step 3b (the owner-cell rule as the default for every point-sampling site) is landing site by site — 3b-i and 3b-ii t1 / t2 / e1 / e2a done, **3b-ii-e2 four sites done** 2026-10-20 (§9 item 568 Part A, e2d): `mag:2`, `mag:4`, `mag:6`, `ports:12` on the owner rule at record width, each example's assertions green; (S-sb) = 0 asserted (same batch, weekly 2026-10-19); cross-batch (S) printed — `mag:4` 10 mismatches ≤ 2 ULP (A1 3 / A2 7), `ports:12` 6 ≤ 32 ULP (all A2); `ports:1` measured on the owner rule by keyword injection (Part B: all four 0.11-image records miss — raw 1.333e-04, corrected 1.308e-04, symmetry 3.005e-06 abs, ‖S‖₂ 3.766e-05 against 1e-06 — identically at `-n 2 / 4 / 1`; (S-sb) 0 and (O) exact on 8 calls at all three widths; pre-empted gates green; cross-width spread of raw / corrected / ‖S‖₂ 2.5e-10 at `-n 1`, 1.1e-09 at `-n 4`, i.e. past 1e-9 at `-n 4`, printed not asserted; 370 / 320 / 342 of 758 multi-incidence probes have owner ≠ legacy cell, deltas up to 0.37), its records returned to the 2026-10-21 weekly, `ports:1` stays on the legacy path on `main`; e2 closes when `ports:1` lands; **3b-ii-e3 done** 2026-10-21 (§9 item 580, re-pose of 556): `th:4`'s profile (DG split) and cross-check (N1curl) and `mag:1`'s R1 table (DG `B`) re-sampled under the owner rule and re-recorded `(1*)`, bands unmoved; `th:4` green at `-n 1 / 2 / 4` (profile 1.591 % at every width), (S-sb) 0 at every width, cross-batch (S) printed (0 / 0 / 4 mismatches at `-n 1 / 2 / 4`, max 1024 ULP, all A2; `cells_ghost` [559, 652, 520, 670] at `-n 4`); (O) exact for both; `mri:1` keyword, (S-sb) 0, one centerline `|B|` digit moved (reported, not re-recorded); CG1 figures kept as labelled records; and 3b-iii (`ans:2`'s 1 g Row 1 re-issued; "write the anchor", returned to the weekly 2026-10-20) open. Known-issues 2026-10-19 (`th:4` at `-n 4`) and the (S) form on ghosted meshes are ruled by the 2026-10-19 weekly (§10). **Done-when:** *(amended 2026-10-02 18:00 review, rule (j), as the 10:30 review directed; the weekly ratifies)* a sampling rule whose output is width-invariant to 1e-9 relative, asserted at `-n 1`, `-n 2` and `-n 4` on `ans:2`'s points; `ans:2`'s Row 1 quoted as a cell or mass average (3a gates both, the weekly picks the published one), width-invariant to 1e-9 relative at the same three widths; the rule adopted by every point-sampling site of the census as the weekly classifies it; and the tracked Row 1 re-issued as the averaged row ⇒ ✅. *History: `docs/planning/chunks/POST-8.md`.* 

### §9 item 586 — 2026-10-22 (`POST-8` 3b-ii-e4, the `ports:1` `(1*)` re-record; lane 1)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** The 2026-10-21 weekly, §10 rulings block: "**`ports:1` under the owner rule — ruled: a `(1*)` re-record (lane 1, `POST-8` 3b-ii-e4).** … `ports:1`'s four records are re-recorded from the owner-rule values at the record width (`-n 2`) under the `(1*)` licence with the example's 1e-6 bands **unmoved**, the legacy values kept in the guide as a labelled record, and the example's own assertions green at `-n 1 / 2 / 4`; the cross-width spread of the S-derived records (1.1e-9 at `-n 4`, the solve's LU floor — 3b-ii-t2's precedent) is printed, never asserted at the 1e-9 point band. e2 closes with it; `attempt/POST-8-3b-ii-e2c-20261020T095000Z` is deleted when it lands. A legacy-path exemption is declined" (`PROJECT_PLAN.md:4736–4754`; lane-1 target `:4997–4998`; this file's weekly ruling, `POST-8.md:1924`). It supersedes the §9 standing line "`ports:1`'s four records are not re-recorded by any slot" (`PROJECT_PLAN.md:3962`). **No band, tolerance, Done-when or target is set or moved:** `REPRODUCTION_BAND_RELATIVE` / `REPRODUCTION_BAND_SYMMETRY_ABSOLUTE` (1e-6, `examples/ports/01_two_torus_port_pair.py:191–192`), `MUTUAL_TOLERANCE`, `S_SYMMETRY_BAND` and `S_SPECTRAL_NORM_CEILING` (`:195–199`) stay as they are; only the four `RECORDED_*` records move, as the weekly ordered.

**`POST-8` 3b-ii-e4: `ports:1` commits the owner-rule keyword (by path from the e2c branch), its four records are re-based `(1*)` onto item 568 Part B's `-n 2` digits, and `test_ports1` runs the example under the spy at `-n 2 / 4 / 1`. e2 and 3b-ii close.**
- **Executor:** `example-runner`.
- **Tier:** heavy (container `timeout -k 30 400`, 568's figure), `-n 2`, then `-n 4`, then `-n 1`.
- **Class:** lane-1 revalidation. **Cost:** ≈ 9.8 measured (lane-1 median) / ≈ 35 predicted slot-minutes: fixed ≈ 5, checkout and edits ≈ 6, windows ≈ 9 (three runs at 162 / 155 / 179 s, `20261020T214423Z_POST-8-step3b-ii-e2d-ports1-measure-n2.log` footer, `…214712Z_…-n4.log`, `…214952Z_…-n1.log`; collect and censuses ≈ 1 min), guide ≈ 4, commit, branch delete and journal ≈ 8, margin ≈ 3.
- **Slot-first only** (heavy; 568 / 580 precedent). **Daytime:** no.
- **Depends on:** nothing. Item 587 (3b-iii) shares no code file with it; see "Order" below.
- **Family:** e2, attempt 1 of the e4 sub-step (568's Part A landing reset the count).

**What exists (read 2026-10-22).**
- **The example on `main`** (`examples/ports/01_two_torus_port_pair.py`):
  - `_path_voltage` samples on the legacy path, with no `owner` keyword (`:303`).
  - The records are at `:162–169`, under `OPS-34`'s comment (`:148–161`). `SUPERSEDED_V072_RECORD` is at `:177–182`.
  - Prints: `[EX-34]` full-precision candidates (`:694–699`) and the `[EX-18]` reproduction misses (`:700–708`).
  - Record assert at `:728–743`. Its message contains "does not reproduce the 0.11-image record" (`:740`).
  - Superseded control `:749–756`. Mutual, blind, symmetry and passivity gates `:758–774`.
- **The e2c branch.** `git diff main attempt/POST-8-3b-ii-e2c-20261020T095000Z -- examples/ports/01_two_torus_port_pair.py` is the keyword hunk only: `owner="lowest_original_cell"` on `_path_voltage`'s call, +3 / −1. Otherwise the branch is far behind `main`: `git diff main <branch> --stat -- examples/ tests/ src/` lists 58 files and 10 169 deletions. That includes **the module, where the branch lacks 568's Part B and 580's `_drive_e3`**, and `docs/testing/test-results.md`, which carries 13 rows.
- **The values to record (568 Part B; the injection forwards the branch's keyword byte for byte).** These are the `[EX-34]` lines:
  - `-n 2` (`…214423Z_…-measure-n2.log:677`): `raw=0.8946356119123071 corrected=0.9399444943082806 symmetry=3.4126695369442954e-05 spectral=0.8613244573231138`;
  - `-n 4` (`…214712Z_…-n4.log:701`): `0.8946356108980973 / 0.939944493262462 / 3.412669016280772e-05 / 0.8613244568626531`;
  - `-n 1` (`…214952Z_…-n1.log:665`): `0.89463561213243 / 0.9399444945352639 / 3.41267301570446e-05 / 0.8613244571329535`.
  - The misses against the legacy (0.11-image) records were identical at all three widths (`n2.log:695`, `n4.log:719`, `n1.log:683`). The pre-empted gates were green (`:697` / `:721` / `:685`).
  - Cross-width spreads: `n4.log:722`, `n1.log:686`.
- **The module** (`tests/post/test_post8_step3b_ii_example_sites.py`):
  - `test_ports1` is skipped (`:647–654`). Its body is `_drive_e2`, which asserts the record width (`:613`).
  - `_drive_e3` (`:623–629`) is `_run_spied` → `_check_e2_s` → `_check_o`, and runs at any width. Item 580 checked `th:4` with it at `-n 1 / 2 / 4`.
  - `test_ports1_owner_measure` (`:667–775`) is item 568's measurement and stays untouched.
- **The guide** (`examples/ports/01_two_torus_port_pair.md:46–47`): "`EX-18` records (177998 cells, raw 0.894516 / corrected 0.939822 / symmetry 3.1121e-05 / spectral 0.861357)".

**Code.** No `src/` edit.
1. **One file by path:** `git checkout attempt/POST-8-3b-ii-e2c-20261020T095000Z -- examples/ports/01_two_torus_port_pair.py`.
   - `git diff --cached --stat` must show that one file, 3 insertions and 1 deletion.
   - **Never the module, `test-results.md` or any log from that branch.**
2. **The records (`(1*)`).**
   - Replace the four `RECORDED_*` values with the `-n 2` `[EX-34]` digits above, verbatim.
   - Re-write the comment at `:148–169`. It should say: re-based by `POST-8` 3b-ii-e4 onto the owner rule under `(1*)` (weekly 2026-10-21, `PROJECT_PLAN.md:4736–4754`); the digits come from `20261020T214423Z_POST-8-step3b-ii-e2d-ports1-measure-n2.log:677`; the bands are unmoved.
   - Keep the `OPS-34` digits as data: `LEGACY_PATH_RECORD_0P11 = {"raw": 0.8945163788281, "corrected": 0.9398215452105435, "symmetry": 3.11213171925739e-05, "spectral": 0.8613568943949574}`. These are the four current constants, verbatim.
   - **Leave the assert message `:739–743` and the print labels unchanged.** `_run_spied_measure` catches on the substring "does not reproduce the 0.11-image record" (module `:243`).
3. **New negative control (asserted, backed; for the review, Q1).** After the superseded control (`:749–756`):
   - assert `_relative_miss(ladder[n], LEGACY_PATH_RECORD_0P11[n]) > REPRODUCTION_BAND_RELATIVE` for `raw` and `corrected`;
   - print all four legacy misses (symmetry absolute) on one `[EX-18]` line.
4. **The module.**
   - Delete the `@pytest.mark.skip(...)` at `:647–651`.
   - Change `test_ports1`'s body from `_drive_e2("ports1")` to `_drive_e3("ports1")`.
   - Add one docstring sentence at `:28–30`: "Item 586 (e4): `ports:1` carries the keyword on `main`; `test_ports1` checks it at `-n 1 / 2 / 4`; `test_ports1_owner_measure` is item 568's record (its injection is now a no-op)."
   - Nothing else changes. `_drive_e2`, `_drive_e3`, every checker and every other test stay byte-identical.
5. **The guide.**
   - Label the `:46–47` figures "legacy-path 0.11-image records (`OPS-34`), superseded by `POST-8` 3b-ii-e4; kept as a labelled record".
   - Add one sentence: `_path_voltage` samples under the owner rule; the four records were re-based `(1*)` at `-n 2` from the Part B log (full filename) and are reproduced by this item's `-n 2` log; the bands are unmoved; the cross-width spread at `-n 4` / `-n 1` is printed (from this item's logs).
   - No other guide line changes (Q3).

All edits are made before any window (rule (i)).

**Anchors (asserted; every band imported from the script itself).**
- **(X)** The example's own assertions under the spied `__main__` run, at each width:
  - the four new records at 1e-6 (relative; absolute on symmetry);
  - the superseded-v0.7.2 control;
  - the new legacy-path control;
  - `MUTUAL_TOLERANCE`, the blind fixture, `S_SYMMETRY_BAND` and `S_SPECTRAL_NORM_CEILING`.
- **(E)** 8 calls carry `owner == RULE`, and the forwarded bytes = `own`.
- **(S-sb) = 0** at the single-incidence rows of each rank's owner batch. `own` is bitwise equal across ranks. The in-cell value at multi-incidence probes is within `CLOSED_FORM_TOL`.
- **(O)** At `-n 4` / `-n 1`, the call count, call lines, `n_inc` and owner original indices equal `step3b-ii-e1-ports1-n2.npz` exactly. This is a new npz name; 568's `…-ports1m-n2.npz` is not touched. Indices only: the field spread is printed.

**Negative controls (rule (e)).**
- **Asserted, backed:** the legacy-path record against the owner-rule run. Part B measured raw 1.333e-04 and corrected 1.308e-04 at all three widths (`n2.log:695`, `n4.log:719`, `n1.log:683`). That is the same comparison on the same fixture. The floor is 1e-6, so the separation is ≈ 131×. The miss is a fixed property of the record, so there is no ceiling.
- **Asserted, existing, backed:** the superseded-v0.7.2 record. Part B measured 3.943e-04 / 3.869e-04 under the rule (`n2.log:697`).
- **Printed only:** the legacy symmetry miss 3.005e-06 absolute (3.0× the band) and the spectral miss 3.766e-05.

**Predictions (printed, never asserted).**
- **Record misses (the `[EX-18]` line).** These are the cross-width spread of the S-derived records against the `-n 2` record, which is what the weekly ordered printed:
  - `-n 2`: 0.00e+00 on all four. Same image, partition and call as Part B, which reproduced e2c's `-n 2` digits exactly.
  - `-n 4`: raw 1.13e-09, corrected 1.11e-09, symmetry 5.21e-12 abs, ‖S‖₂ 5.35e-10.
  - `-n 1`: 2.46e-10 / 2.42e-10 / 3.48e-11 abs / 2.21e-10.
  - Headroom against 1e-6 is ≥ 880×.
- **Legacy-path misses** at every width: 1.333e-04 / 1.308e-04 / 3.005e-06 / 3.766e-05.
- **Gates:** corrected −6.0056 %, blind −98.2573 %, symmetry 3.4127e-05, ‖S‖₂ 0.861324.
- **Spy:** 8 calls, 23 826 / 24 584 single-incidence. (R) owner ≠ legacy at 370 / 320 / 342 of 758.
- **Elapsed:** ≈ 162 / 155 / 179 s.

**Windows.** Bare, repo-relative and foreground, with Bash `timeout: 660000` and pytest `-s`. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before each compute window.
1. **Pre-censuses:**
   - `scripts/testing/run_and_log.sh POST-8-step3b-ii-e4-pre-figures "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py'"`
   - `scripts/testing/run_and_log.sh POST-8-step3b-ii-e4-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
2. **Code:** steps 1–5 above.
3. **Collect:** `scripts/testing/run_and_log.sh POST-8-step3b-ii-e4-W0-complex "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/post/test_post8_step3b_ii_example_sites.py --collect-only -q -s'"`
4. **`-n 2` (first: it writes the npz the other widths read):**
   `scripts/testing/run_and_log.sh POST-8-step3b-ii-e4-ports1-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_POST8_EXAMPLES=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/post/test_post8_step3b_ii_example_sites.py::test_ports1 -v -s --tb=short'"`
5. **`-n 4`:** the same command with `mpiexec -n 4`, labelled `POST-8-step3b-ii-e4-ports1-n4`.
6. **`-n 1`:** the same command with `mpiexec -n 1`, labelled `POST-8-step3b-ii-e4-ports1-n1`.
7. **Guide, then post-censuses** (window 1 relabelled `…-post-figures` / `…-post-docrefs`: counts unchanged, `broken=0`, `dead=0`).
8. **Checks:** `git status --porcelain examples/ tests/ src/ scripts/ docker/` lists only the script, its guide and the module. Then run the leak check on the staged logs, as 568 did.

Each compute window reads **5 passed, 0 skipped**: 4 from `test_complex_mode.py`, plus `test_ports1`.

**Traps already paid for.**
- **Path checkout from a stale branch** (`8d4cf58`; this branch's own diff above). Take one file only.
  - The branch's module would silently delete 568's Part B and 580's `_drive_e3`.
  - Its `test-results.md` would delete `main`'s rows. Append rows by hand; `git diff --stat` on record files shows insertions only.
- **`.eval` bytes depend on the batch** (item 537; §9 standing ruling at `PROJECT_PLAN.md:3966`). No cross-batch byte equality is asserted:
  - (S-sb) is same-batch;
  - (O) asserts indices, not field bytes;
  - the S-derived cross-width spread is printed. The 1.1e-9 at `-n 4` is past the 1e-9 point figure; never assert it.
- **Record digits come from Part B's log and go in before any window.** Never edit the script after a window to record this item's own digits (rule (i)).
  - If this item's `-n 2` `[EX-34]` line differs from Part B's at all, report the difference.
  - A difference inside 1e-6 is green and is a reading. The record stays Part B's.
- The measure test's catch keys on the assert message, so leave it unedited.
- `-n 2` runs first.
- `logs/post8/` is gitignored, so the npz is never committed.
- Never set `FEM_EM_SETUP_FIGURES`.
- **Leak guard:** `[EX-34]` full-precision prints passed it before (`386f45e9`). On a flag, move the log unchanged to `logs/leak-held/`, mark the item BLOCKED and stop.
- **A kill is not a result.** Run `pgrep -c python3` first. At 0, run the rubric's cache sweep labelled `POST-8-step3b-ii-e4-cache-sweep`. Never raise a timeout or `-n`.
- The bash hook denies `pgrep -f "python3 -m pytest"`; `pgrep -c python3` is the check.
- No AED number anywhere. Never open `docs/private/` or `aed_results/`.

**Scope.** `ports:1` on the owner rule on `main`; its four records re-based `(1*)`; the module's `test_ports1`; one guide sentence and one label. **Not claimed:**
- 3b-iii (item 587);
- any band;
- width invariance of the S-derived records at 1e-9 (printed only);
- the guide's other stale lines (Q3);
- `POST-8` ✅, unless item 587 has already landed (below).

**Negative result.** Report and stop. Park on `attempt/POST-8-3b-ii-e4-*` and mark item 586 BLOCKED in the same commit (rule (d)). The e2c branch stays.
- **(X) red on a record at `-n 2`:** the committed keyword does not reproduce the injected run. Print the four values beside Part B's.
- **(X) red on a record at `-n 4` / `-n 1`:** a spread over 1e-6, which contradicts Part B's 1.1e-9. Report both logs.
- **Legacy control under its floor:** the keyword is not active. Check `git diff` on the script and stop.
- **(S-sb) ≠ 0 or (O) red:** a finding against the rule on this mesh. Open a known-issues entry.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.**
- **The `POST-8` §7 state line** (`PROJECT_PLAN.md:2685`). Replace the `ports:1` clause with: "**3b-ii-e4 done** <date> (§9 item 586): `ports:1` on the owner rule, four records re-based `(1*)` at `-n 2`, 1e-6 bands unmoved, green at `-n 1 / 2 / 4`; cross-width spread of the S-derived records printed (<-n 4> / <-n 1>); e2 and 3b-ii closed".
- Known-issues 2026-09-27's "Verified at" gains e4, if the entry is still open.
- `git branch -D attempt/POST-8-3b-ii-e2c-20261020T095000Z`.
- **`POST-8` → ✅ only if item 587 has landed first.** Under rule (j), the row's four Done-when clauses are then all met: 3a; 3b-iii's Row 1; every census site, with `ports:1` the last; and the re-issued Row 1. Otherwise the row stays 🟡 with "3b-iii open".

**For the review.**
1. **The new legacy-path control** (code step 3). The weekly wrote "legacy values kept in the guide as a labelled record". The drafter also keeps them as data in the script and asserts raw / corrected against them, so that a silent revert to the legacy path fails the example. Confirm, or strike it (the values then print only).
2. **`test_ports1_owner_measure`** stays unchanged. Its injection is a no-op once the keyword is on `main`. Retiring it is a later sweep's call.
3. **Stale guide lines not touched:**
   - `:104` says "reproduce the record to 2e-3", but the band has been 1e-6 since `OPS-34`;
   - the §3 table digits (0.894543 …) predate `OPS-34`;
   - `:80–81` says symmetry / spectral are "printed, never asserted", but they have been asserted since `OPS-34`.

   These are for an example-sync item, not this one.
4. **Order with item 587.** The two items share no code file. They do share `PROJECT_PLAN.md`'s §7 row, known-issues 2026-09-27 and `test-results.md`. Each depends on nothing. The drafter's order is 586 first, so that 587 closes `POST-8`. Both are heavy and slot-first, so never in one slot.

**Execution order.** pgrep → pre-censuses → code (1–5) → collect → `-n 2` → `-n 4` → `-n 1` → guide → post-censuses → `git status` / leak check → commit (script, guide, module, logs with `test-results.md` rows by hand, §7 state line, §9 item 586 DONE, known-issues line, this file's outcome line) → branch delete → journal (own commit).

**Executor spawn prompt (verbatim).** "Use the example-runner agent to execute `POST-8` §9 item 586 (step 3b-ii-e4, the `ports:1` `(1*)` re-record) per `docs/planning/chunks/POST-8.md`: read item 568's brief and outcome, the 2026-10-21 weekly ruling and item 586's brief, in that order; 586 binds, and any review ruling beneath it binds over it.
- You are the executor: do not spawn agents, and never background a harness window.
- The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.
- Take exactly one file from `attempt/POST-8-3b-ii-e2c-20261020T095000Z` by path: `examples/ports/01_two_torus_port_pair.py`. Never the module, `test-results.md` or a log.
- Record the four `-n 2` digits from `20261020T214423Z_POST-8-step3b-ii-e2d-ports1-measure-n2.log:677` verbatim, before any window. Keep the `OPS-34` digits as `LEGACY_PATH_RECORD_0P11`. Leave the record assert message unchanged. Move no band.
- Run every window in the foreground with `timeout: 660000`, pytest `-s`, container `timeout -k 30 400`, in the order `-n 2`, `-n 4`, `-n 1`. `pgrep -c python3` must read 0 before each compute window.
- Never set `FEM_EM_SETUP_FIGURES`. Never open `docs/private/` or `aed_results/`.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-ii-e4-*`, mark item 586 BLOCKED in the same commit, and report. On green, delete the e2c branch after the commit.
- Return:
  - the log names with elapsed times;
  - per width: the four values and their misses against the new record, the legacy misses, the superseded misses, (E), (S-sb), the cross-batch (S) count with max ULP, `cells_ghost`, (O) and the field spread;
  - the census before and after;
  - the branch delete."

### §9 item 587 — 2026-10-22 (`POST-8` 3b-iii, `ans:2`'s Row 1 re-issued as the owner-cell average; lane 1)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** The 2026-10-21 weekly, §10: "**`POST-8` 3b-iii — the anchor, written.** … **the published Row 1 is the owner-cell average** (the DG-split SAR of the owner cell at each probe point) — Row 4 is already the 1 g mass average and a second 1 g row would say nothing. 3b-iii therefore:
- (i) adds `FEM_EM_ANS_SKIP_PRIVATE` to `ans:2`'s `main()` on `ans:4`'s pattern;
- (ii) re-issues the tracked Row 1 as the owner-cell values, width-invariant at `-n 1 / 2 / 4` to 1e-9 (3a's band, imported), the legacy CG1 values kept as a labelled record;
- (iii) amends `SPEC.md` Row 1's definition to 'SAR in the owner cell at the probe point' as a v-bump with the ±10 % pointwise band unchanged …;
- (iv) keeps `metrics.json["pointwise"]` as the key with its values re-issued and a `"pointwise_rule": "owner_cell"` field beside it, and re-runs its consumer `test_ans9_far_side_sar.py:220–228` (green, or re-recorded `(1*)` if it compares a stored number);
- (v) retires known-issues 2026-09-27 on (ii).

… `POST-8` → ✅ needs e4 and 3b-iii both" (`PROJECT_PLAN.md:4755–4774`).

The §7 row's Done-when already reads "quoted as the owner-cell average" (`PROJECT_PLAN.md:2685`). Earlier blockers, now lifted by (i): the 10-20 03:00 ruling (6) (`POST-8.md:1615`) and item 526's 3b-iii note (`:234–237`). **No band, tolerance, Done-when or target is set or moved:** `WIDTH_RTOL`, `ROUTE_IMAG_RTOL` and `PLANT` are imported from 3a's module, and the SPEC's ±10 % stays.

**The reading of "owner-cell average" (Q1).** The Done-when says "3a gates both" and "picked the owner-cell one — Row 4 is already the 1 g average". 3a's two averaged families are `cell` and `one_gram` (`tests/validation/test_post8_step3a_averaged_row1.py:15–24`). So the published Row 1 is **3a's `cell` family**: the mean of σ|E|²/(2ρ) over the owner cell, by route (ii), the DG0 indicator (`:204–236`). The weekly's parenthesis, "the DG-split SAR of the owner cell at each probe point", could instead be read as 3a's `raw` family, the owner-rule point value. The two differ: cell/raw was 0.639 … 1.596 (`20261003T033225Z_POST-8-step3a-W1.log.gz:1953`). This brief is written on the Done-when's reading, and the alternative swaps one function and one record (Q1).

**`POST-8` 3b-iii: `ans:2` honours `FEM_EM_ANS_SKIP_PRIVATE`; its tracked Row 1 becomes the owner-cell average, asserted against 3a's `-n 2` record at 3a's 1e-9 at `-n 2 / 1 / 4`; the CG1 row is kept as a labelled record and the negative control; the SPEC is v-bumped; the consumer is re-run; known-issues 2026-09-27 retires.**
- **Executor:** `example-runner`.
- **Tier:** heavy, complex build.
  - Three `ans:2` `main()` windows at `-n 2 / 1 / 4`, each with container `timeout -k 30 600`.
  - One consumer window at `-n 8` with `timeout -k 30 480`.
- **Class:** lane-1 revalidation. **Cost:** ≈ 9.8 measured (lane-1 median) / ≈ 55 predicted slot-minutes: fixed ≈ 5, code ≈ 15, SPEC and guide ≈ 6, windows ≈ 16, commit and journal ≈ 8, margin ≈ 5. The windows break down as:
  - compile ≈ 0.5 min;
  - `-n 2` ≈ 245 s (`20260926T140900Z_EX-57.log.gz`, 245 s, the same script at `-n 2`);
  - `-n 1` ≈ 305 s (scaled by 3a's `-n 1` / `-n 2` ratio, 338 / 272 s);
  - `-n 4` ≈ 195 s (`20260925T110348Z_OPS-63.log.gz`, 195 s);
  - consumer ≈ 140 s (`20261011T184301Z_ANS-9-step1b-W1-d1.log`, 136 s at `-n 8`);
  - censuses and leak check ≈ 1 min.

  The POST-8 row prices "steps 2–3a heavy ≈ 190–290 s per window" (`PROJECT_PLAN.md:2685`); the `-n 1` window is above that range.
- **Park (rule (d)):** if the `-n 2` window has not started by minute 25, the `-n 4` window by minute 40, or the commit by minute 57.
- **Slot-first only.** **Daytime:** no.
- **Depends on:** nothing. Item 586 shares no code file (see "Order").
- **Family:** 3b-iii, attempt 1.

**What exists (read 2026-10-22).**
- **`ans:2`'s script** (`examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/02_birdcage_coil_driven_sar_10MHz.py`):
  - `_pointwise_sar` (`:310–349`) uses the CG1 interpolant and `evaluate_vector_field_parallel` with no owner.
  - In `main()` (`:675`): Row 1 is computed at `:757`, written to the payload's `"pointwise"` at `:805`, and printed at `:824–833`.
  - `_write_comparison` (`:442`) writes Row 1's heading at `:503` and fills `point_rows` from both `sar_w_per_kg` and `e_magnitude_v_per_m` (`:460–471`, `:509`, `:515`).
  - **The private write is unconditional when `aed_results/` exists** (`:1056–1061`), and it exists on this box.
  - `RESOLUTION_OVERRIDDEN` (`:200`) skips only the mesh-record asserts (`:1025–1040`), with a printed line.
  - `_print_step4_readout` (`:610–660`) reads `ref["pointwise"]`.
- **`ans:4`'s pattern:** the constant at `04_birdcage_four_port_10_64_128MHz.py:192–196` and the `if skip_private: … elif AED_RESULTS_PATH.exists():` branch at `:1026–1032`.
- **3a's module:**
  - `WIDTH_RTOL = 1e-9` (`:88`), `ROUTE_IMAG_RTOL = 1e-14` (`:92`), `PLANT = 1e-8` (`:94`);
  - `point_owner_cells` (`:167`);
  - route (ii) (`:204–236`): a DG0 χ on the owning rank's `owner_local` cell, `scatter_forward`, `quadrature_degree` 2, both integrals `allreduce`d.
- **3a's logs** (tracked, gz):
  - W1 `-n 2`, `20261003T033225Z_POST-8-step3a-W1.log.gz`:
    - owner original indices `[150135, 162464, 142550, 169932]` at `:1926`;
    - `cell(ii)` rows at `.17e` at `:1936–1939`;
    - CG1 rows at `:1944–1947`;
    - 272 s.
  - W2 `-n 4` (`…033703Z_…W2.log.gz`): (W) cell 2.150e-10, 231 s.
  - W3 `-n 1` (`…034100Z_…W3.log.gz`): (W) cell 1.838e-10, 338 s.
  - The `-n 2` npz `logs/post8/step3a-n2.npz` (mtime 10-02 22:36) holds the same `cell` array. It is not read by this item.
- **Other readers of `_pointwise_sar`, which must stay CG1:**
  - the step-1 probe (`scripts/probes/post8_pointwise_partition.py`);
  - step 2's C1a / C1b (`tests/validation/test_post8_owner_rule_ans2.py`);
  - 3a's C1b (`test_post8_step3a_averaged_row1.py:264`).

  `ans:9` loads the script only for `_excitation` (`09_…py:577–579`).
- **`metrics.json`:** generated 2026-09-19 at `-n 4`. `"pointwise"` is at `:40–97`.
- **The consumer** (`tests/validation/test_ans9_far_side_sar.py:220–228`, column `d1`) **prints only; it compares no stored number.** Its last `d1` run (`20261011T184301Z_ANS-9-step1b-W1-d1.log`, prints `:2020–2023`, 7 passed, 1 skipped, 136 s) was on `25686eb7`. `50a5a43a` edited the module after that, so the `d1` column has not run on the module as committed.
- **`SPEC.md`:**
  - Row 1 is `:165`.
  - **Row 3 is `:167`, which `scripts/testing/check_private_leak.py:120–121` cites as `published_in`.**
  - The 09-18 adjudication status is at `:39–46`. It was read on the CG1 row.
- **The guide:** `02_birdcage_coil_driven_sar_10MHz.md:23` (the Row 1 table line).

**Code.** No `src/` edit.
1. **`ans:2`'s script.**
   - **(a) Skip-private (step (i)).** Add `SKIP_PRIVATE_ENV = "FEM_EM_ANS_SKIP_PRIVATE"` beside `AED_RESULTS_PATH`, with `ans:4`'s comment. Change `main()`'s private write to `ans:4`'s `if skip_private … elif AED_RESULTS_PATH.exists()` structure, and print the note. When the variable is unset, behaviour is unchanged.
   - **(b) Records, transcribed verbatim from 3a's W1 log:**
     - `OWNER_CELL_ROW1_RECORD_W_PER_KG`, a 4 × 4 array indexed [drive P1…P4][point +x, +y, −x, −y], from `:1936–1939`, all 17 digits;
     - `OWNER_CELL_ORIGINAL_INDEX_RECORD = (150135, 162464, 142550, 169932)`, from `:1926`.

     The comment cites both log lines and names 3a's npz. A transcription slip at any digit that matters fails (W) at every width, so it cannot pass silently.
   - **(c) Imports:**
     - `point_owner_cells` (with `evaluate_vector_field_parallel`);
     - `ufl`;
     - `WIDTH_RTOL`, `ROUTE_IMAG_RTOL` and `PLANT` from `tests.validation.test_post8_step3a_averaged_row1`.
   - **(d) New `_owner_cell_sar(built, points, comm)`:** 3a's route (ii), copied from `:204–236`, with the owner from `point_owner_cells(msh, points, comm)`. It returns `(rows, owner_original, imag_ratio_max)`.
     - `rows[k]` has `"sar_w_per_kg"` (the owner-cell mean) and `"e_magnitude_v_per_m"` = √(2ρ·SAR/σ), the cell-RMS |E| (Q2).
     - The docstring names `POST-8` 3b-iii and the weekly's anchor.
   - **(e) Leave `_pointwise_sar` byte-identical.**
   - **(f) `main()`:**
     - Row 1 = `_owner_cell_sar`.
     - The legacy row = `_pointwise_sar`, printed as "LEGACY CG1 route — width-dependent, a labelled record, never published".
     - Payload: `"pointwise"` holds the owner-cell rows, plus `"pointwise_rule": "owner_cell"` and `"pointwise_owner_original_cell_index"`.
     - The ROW 1 print heading names the owner-cell mean and prints `.17e` rows. It also prints the max relative miss against the record and the negative control's minimum separation.
     - The docstring's `metrics.json` bullet (`:19–22`) is updated to match.
   - **(g) Asserts.**
     - (Im) and (N) run at every rung.
     - (Own) and (W) run only when `not RESOLUTION_OVERRIDDEN`, beside the mesh-record block, with a printed skip line otherwise (the record is the 0.0025 m rung's).
   - **(h) `_write_comparison`:** owner-cell wording on Row 1's heading (`:503`) and the `|E|` caption (`:511`), plus one sentence saying Row 1 was re-issued by `POST-8` 3b-iii and that the CG1 row it replaces is a labelled record in the guide. The AED columns stay verbatim blank.
2. **`SPEC.md` (step (iii)).**
   - Re-write Row 1 at `:165` **in place, as one line**: "**SAR in the owner cell at the probe point**", the mean of `σ|E|²/(2ρ)` over the mesh cell that owns each named point (lowest original index among the cells containing it). The points stay the same.
   - Put the v-bump note **below the table** (after `:175`, before the named-points paragraph). It says: SPEC v2, 2026-10-2x, `POST-8` 3b-iii; Row 1 is re-defined; the ±10 % band is unchanged; the 2026-09-18 AGREE on rows 1–3 was read on the CG1 row; the owner-cell row's re-read is private and pending; no number.
   - **No line above `:167` may be added or removed** (the leak-guard pointer).
3. **The guide.**
   - Re-label `:23`.
   - Add one short "Row 1 re-issued" paragraph: the owner-cell rule, the three logs, and the max miss per width. Keep the **legacy CG1 row as a labelled record**: the 16 values of `"pointwise"` in `metrics.json` at HEAD before this commit (`-n 4`, 2026-09-19), at `.9e`, labelled "CG1-interpolated, width-dependent (known-issues 2026-09-27), superseded".
4. **The consumer (one print label only).** At `test_ans9_far_side_sar.py:227–228`, change "ANS-2 CG1-interpolated metrics.json" to `f"ANS-2 metrics.json ({m.get('pointwise_rule', 'cg1')})"`. Nothing else; it compares no stored number, so there is no `(1*)` re-record (Q5).
5. **`metrics.json` / `COMPARISON.md`:** every width rewrites them, and **the committed pair is the `-n 4` window's** (the tracked file's width since 09-19; Q3).

All edits are made before W0 (rule (i)).

**Anchors (asserted; every band imported).**
- **(Own)** `owner_original` equals `OWNER_CELL_ORIGINAL_INDEX_RECORD` exactly at every width. 3a's (W) asserted this at `-n 1 / 2 / 4`.
- **(W)** Each of the 16 owner-cell SARs is within `WIDTH_RTOL` (`test_post8_step3a_averaged_row1.py:88`) relative of the `-n 2` record, at `-n 2`, `-n 1` and `-n 4`. Every width is held to one `-n 2` comparand, which `-n 2` reproduces to ≤ 1e-12 (predicted). That is 3a's (W) form (each width against the `-n 2` value) and step 2's A2 logic (`test_post8_owner_rule_ans2.py:22–25`).
- **(Im)** Route (ii)'s |Im| / Re ≤ `ROUTE_IMAG_RTOL` (`:92`). 3a measured round-off, about 2.1e-15.
- **Every existing `ans:2` anchor stays asserted and unchanged:** the C4 pairs, coverage, mesh records, PORT-9, containment and the 1 g control.
- **The consumer's own assertions are green** at `-n 8`.

**Negative control (rule (e)).**
- **Asserted, backed:** min over the 16 entries of |CG1 / owner-cell − 1| > `PLANT` (`:94`, 1e-8). Measured on the same fixture from 3a's own prints of both rows in one run:
  - `-n 2`: 1.190e-02 (W1 `:1936–1947`);
  - `-n 4`: 2.022e-03 (W2);
  - `-n 1`: 1.413e-02 (W3).

  These are the drafter's stdlib arithmetic on the tracked logs. The headroom is ≥ 2.0e5×. The ceiling is not binding: the CG1 route's own error is what is being separated.
- **Printed:** the CG1 row beside 3a's CG1 row at the same width (≤ 1e-9 predicted, step 2's C1a 5.89e-10). At `-n 4`, also the CG1 row beside the tracked 09-19 `"pointwise"` values.

**Predictions (printed, never asserted).**
- (W) miss: `-n 2` ≤ 1e-12 (point 0, if `_build_mass_averaged` reproduces 3a's `case` bit for bit); `-n 4` ≈ 2.15e-10; `-n 1` ≈ 1.84e-10.
- (Own) exact.
- Im ratio ≤ 1e-14.
- CG1 min separation ≈ 1.19e-02 / 2.02e-03 / 1.41e-02.
- The published (`-n 4`) Row 1 equals 3a's `cell(ii)` to ≈ 2e-10.
- The driven-point `_spread` changes from the CG1 row's. That is a new reading, printed.
- The consumer is green, and its four `ANS-2 named point` lines carry the `owner_cell` label.
- Elapsed ≈ 245 / 305 / 195 / 140 s.

**Windows.** Bare, repo-relative and foreground, with Bash `timeout: 660000`. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before each compute window. No window exceeds ≈ 10 min, so no durable capture is needed.
1. **Pre-census:**
   `scripts/testing/run_and_log.sh POST-8-step3b-iii-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
2. **Code:** steps 1–5 above.
3. **W0 (compile):**
   `scripts/testing/run_and_log.sh POST-8-step3b-iii-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 -m py_compile examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/02_birdcage_coil_driven_sar_10MHz.py tests/validation/test_ans9_far_side_sar.py && echo compiled'"`
4. **W1 (`-n 2`):**
   `scripts/testing/run_and_log.sh POST-8-step3b-iii-ans2-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 600 mpiexec -n 2 python3 examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/02_birdcage_coil_driven_sar_10MHz.py'"`
5. **W2 (`-n 1`):** the same command with `mpiexec -n 1`, labelled `POST-8-step3b-iii-ans2-n1`.
6. **W3 (`-n 4`; its `metrics.json` / `COMPARISON.md` are the committed pair):** the same command with `mpiexec -n 4`, labelled `POST-8-step3b-iii-ans2-n4`.
7. **After each of W1–W3, as its own Bash call:** `git status --porcelain examples/`. It should list only the edited script, `SPEC.md`, the guide, `metrics.json` and `COMPARISON.md`. Each log must carry the `skipped (FEM_EM_ANS_SKIP_PRIVATE=1)` note.
8. **W4 (the consumer, after W3):**
   `scripts/testing/run_and_log.sh POST-8-step3b-iii-ans9-consumer "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/ans9_1b && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS9_1B_COLUMN=d1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 480 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_far_side_sar.py -v -s --tb=short'"`
9. **Post-census** (window 1 relabelled `…-post-docrefs`: `dead=0`), then the leak check on the staged files and logs.

**Traps already paid for.**
- **`_pointwise_sar` stays CG1.** Step 1's probe, step 2's C1a / C1b and 3a's C1b assert on it. Flipping it would turn three closed gates red.
- **The private write.** `aed_results/` exists on this box, so an unset knob rewrites `COMPARISON_private.md` (`:1056–1061`), which `EX-57.md:40` forbids. Every window sets `FEM_EM_ANS_SKIP_PRIVATE=1`.
  - The private AED re-read of the new row is **not this slot's**. No slot opens `aed_results/`, `COMPARISON_private.md` or `docs/private/`.
- **Never set `FEM_EM_ANS2_PHANTOM_RESOLUTION`**: the record is the 0.0025 m rung's. Never set `FEM_EM_SETUP_FIGURES` either.
- **`SPEC.md:167` is the leak guard's `published_in` pointer.** Edit Row 1 in place, and put the note below the table.
- **The stale-branch / record-file trap (`8d4cf58`).** Nothing is taken from any branch. `metrics.json` / `COMPARISON.md` are rewritten by every width, so commit W3's pair. `test-results.md` rows go in by hand, and `git diff --stat` on record files shows insertions only.
- **`.eval` bytes depend on the batch** (item 537). Route (ii) is an assembly, not an `.eval`. No byte equality is asserted across widths; (W) is relative at 1e-9, and the CG1 control compares within one run.
- **Pin `metadata={"quadrature_degree": 2}`** as 3a did. Reduce each `assemble_scalar` with `comm.allreduce` (they are rank-local). χ is complex in this build, so set `1.0` and read `.real`.
- **The consumer's `d1` column rewrites `/workspace/logs/ans9_1b/d1.json`** (gitignored). It has not run on the module as committed since `50a5a43a`. A red there that does not involve the print is a pre-existing finding: re-run with the label edit stashed, labelled `…-baseline`, and report.
- **Leak guard.** The new `.17e` prints are ours, but a digit-prefix collision held item 526's W5b log. On a flag, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- **A kill is not a result.** Run `pgrep -c python3`, then at 0 the rubric's cache sweep labelled `POST-8-step3b-iii-cache-sweep`. Never raise a timeout or `-n`.
- No AED number anywhere, including the SPEC note and the journal.

**Scope.** `ans:2`'s knob; its tracked Row 1 re-issued as the owner-cell average; the SPEC v-bump; the guide's labelled CG1 record; one consumer label. **Not claimed:**
- the private AED re-read of Row 1, or any verdict change: the 09-18 AGREE on rows 1–3 stands as read on the CG1 row until that re-read, which is the operator's or the weekly's, in the gitignored files only;
- rows 4–6;
- the step-4 h-rung readout;
- e4 (item 586);
- any band;
- `POST-8` ✅, unless item 586 has already landed.

**Negative result.** Report and stop. Park on `attempt/POST-8-3b-iii-*` and mark item 587 BLOCKED in the same commit (rule (d)). Restore `metrics.json` / `COMPARISON.md` on `main` from `HEAD`; the branch carries them.
- **(W) red at `-n 2`:** the transcription, or `ans:2`'s builder diverging from 3a's `case`. Print all 16 misses. Never touch `WIDTH_RTOL`.
- **(W) red at `-n 1` / `-n 4` only:** this contradicts 3a's 2.15e-10. Print the per-entry readings and annotate known-issues 2026-09-27; the entry stays open.
- **(Own) red:** the ownership has moved since 3a. That is a finding.
- **(Im) red:** a route defect.
- **Consumer red:** see the trap above.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.**
- **The `POST-8` §7 state line:** "**3b-iii done** <date> (§9 item 587): `ans:2` honours `FEM_EM_ANS_SKIP_PRIVATE`; tracked Row 1 re-issued as the owner-cell average (`pointwise_rule` `owner_cell`), within 1e-9 of 3a's `-n 2` record at `-n 2 / 1 / 4` (<three misses>); CG1 row a labelled record; SPEC v2; consumer green".
- **Known-issues 2026-09-27 retires** on (ii), with the three logs.
- **`POST-8` → ✅ if item 586 has landed.** Under rule (j), the row's Done-when (`PROJECT_PLAN.md:2685`) is then met clause by clause:
  - 3a's rule, width-invariant at three widths;
  - Row 1 quoted as the owner-cell average at 1e-9 at three widths (this item);
  - every census site on the rule (e4 last);
  - the tracked Row 1 re-issued (this item).

  The review sends the closure to the `auditor`. If 586 has not landed, the row stays 🟡 with "e4 open".
- §2 cites neither the known-issues entry nor `POST-8` (checked), so no §2 sentence moves.

**For the review.**
1. **Owner-cell average or owner-rule point value?** The drafter follows the Done-when (`:2685`): 3a's `cell` family. If the weekly's parenthesis means 3a's `raw` family, Row 1 becomes `point_sar(e_real, e_imag, points, sigma=…, rho=…)` (`post/sar.py:332`, owner rule by default), and the record moves to W1 `:1928–1931`. Everything else stands.
2. **`e_magnitude_v_per_m` under the cell rule** is the cell-RMS |E|, √(2ρ·SAR/σ). It is a derived column. Alternatively, `COMPARISON.md`'s `|E|` sub-table could be dropped.
3. **Record width of the committed `metrics.json`:** `-n 4` (the tracked file's width since 09-19, and the script's documented run) or `-n 2`.
4. **SPEC wording.** The drafter adds the owner-rule definition to the weekly's "SAR in the owner cell at the probe point". The v-bump note sits below the table so that `:167` holds.
5. **The consumer's label edit** keeps a print from calling owner-cell values "CG1-interpolated". Without it, the module runs unedited.
6. **Who reads the new Row 1 against AED privately,** and when the SPEC status line may restate the Row 1 verdict.
7. **Order with item 586:** independent code, with shared record files (the §7 row, known-issues 2026-09-27, `test-results.md`). Each depends on nothing. The drafter's order is 586 then 587, so that **587 closes `POST-8`**. Run the other way, 586 closes it and skips its known-issues line. Never in one slot.

**Execution order.** pgrep → pre-census → code (1–5) → W0 → W1 (`-n 2`) → `git status` → W2 (`-n 1`) → `git status` → W3 (`-n 4`) → `git status` → W4 (consumer) → post-census → leak check → commit (script, `SPEC.md`, guide, W3's `metrics.json` / `COMPARISON.md`, the consumer label, logs with `test-results.md` rows by hand, §7 state line, §9 item 587 DONE, known-issues retirement, this file's outcome line) → journal (own commit).

**Executor spawn prompt (verbatim).** "Use the example-runner agent to execute `POST-8` §9 item 587 (step 3b-iii, `ans:2`'s Row 1 re-issued as the owner-cell average) per `docs/planning/chunks/POST-8.md`: read item 188's brief and outcome (3a), the 2026-10-21 weekly ruling and item 587's brief, in that order; 587 binds, and any review ruling beneath it binds over it.
- You are the executor: do not spawn agents, and never background a harness window.
- The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.
- Edit only the files the brief names: the `ans:2` script, `SPEC.md` (Row 1 in place, the note below the table, no line added or removed above line 167), the `ans:2` guide, and one print label in `tests/validation/test_ans9_far_side_sar.py`. No `src/` file.
- Leave `_pointwise_sar` byte-identical. Transcribe the two records verbatim from `20261003T033225Z_POST-8-step3a-W1.log.gz:1926` and `:1936–1939`. Import `WIDTH_RTOL`, `ROUTE_IMAG_RTOL` and `PLANT` from `tests.validation.test_post8_step3a_averaged_row1`. Move no band.
- Every `ans:2` window sets `FEM_EM_ANS_SKIP_PRIVATE=1`. Never set `FEM_EM_ANS2_PHANTOM_RESOLUTION` or `FEM_EM_SETUP_FIGURES`. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. No AED number anywhere.
- Run W0, then W1 `-n 2`, W2 `-n 1`, W3 `-n 4`, then W4 (the consumer, `-n 8`), in the foreground with `timeout: 660000` and container `timeout -k 30`, with pytest `-s` in W4. `pgrep -c python3` must read 0 before each compute window. Commit W3's `metrics.json` / `COMPARISON.md`.
- If the docker socket is denied, stop and report the exact command string. If the leak guard fires, move the log unchanged to `logs/leak-held/`, mark BLOCKED and stop.
- On a red, park on `attempt/POST-8-3b-iii-*`, mark item 587 BLOCKED in the same commit, restore `metrics.json` / `COMPARISON.md` on `main` from `HEAD`, and report.
- Return:
  - the log names with elapsed times;
  - per width: the owner indices, the max (W) miss, the Im ratio, the CG1 min separation and the `_spread`;
  - the consumer's pass count and its four printed lines;
  - the census before and after;
  - the known-issues retirement text."

### Review ruling 2026-10-22 03:00 (items 586 e4, 587 3b-iii)

Binds over the briefs. Order: 586 then 587 (both heavy and slot-first, never in one slot); whichever lands second closes `POST-8` → ✅ if the other is green.
- **586:** the new negative control (the legacy record misses the owner-rule run) is **asserted** — backed by Part B's 1.333e-04 on the same comparison and fixture (`…214423Z…-measure-n2.log`). `test_ports1_owner_measure` is kept as is. The three stale guide lines are named in the outcome, not edited. Only the `ports:1` script is taken from `attempt/POST-8-3b-ii-e2c-20261020T095000Z` (its test module is stale); the branch is deleted on green.
- **587:** (1) Row 1 is **3a's owner-cell family** — the DG-split SAR of the owner cell at each probe point, as the weekly wrote and the Done-when now reads. (2) The cell value is 3a's function, the cell average of the SAR density σ\|E\|²/(2ρ), never \|⟨E⟩\|²; no new function. (3) The committed `metrics.json` is the `-n 2` one (the record width). (4) The consumer's label edit is allowed, disclosed. (5) The private AED re-read of the new row is the operator's or the weekly's, never a slot's.

### Outcome of item 586 (2026-10-22, 11:00 UTC slot)

DONE. `ports:1` carries `owner="lowest_original_cell"` (one file by path from the e2c branch; module, `test-results.md` and logs not taken); four records re-based `(1*)` to Part B's `n2.log:677` digits, bands unmoved, `LEGACY_PATH_RECORD_0P11` kept as data with the legacy-path control asserted (raw 1.33e-04 / corrected 1.31e-04, symmetry 3.01e-06 abs, ‖S‖₂ 3.77e-05 printed at every width). `test_ports1` un-skipped, body `_drive_e3`. Windows (5 passed each): `20261022T110209Z_POST-8-step3b-ii-e4-ports1-n2.log` 156 s, `…110501Z_…n4.log` 148 s, `…110738Z_…n1.log` 176 s. Misses against the new record (raw / corrected / symmetry abs / ‖S‖₂): `-n 2` 4.57e-13 / 4.48e-13 / 4.82e-11 / 1.26e-10; `-n 4` 5.47e-10 / 5.37e-10 / 2.04e-11 / 5.09e-10; `-n 1` 2.46e-10 / 2.41e-10 / 3.48e-11 / 2.21e-10. Superseded v0.7.2 misses 3.94e-04 / 3.87e-04 at all widths. (E) 8 calls, (S-sb) 0, cross-batch (S) 0 mismatches max 0 ULP at all widths; `cells_ghost` [1535, 1565] / [1649, 1594, 1099, 1071] / [0]; (O) exact at `-n 4` / `-n 1`, field spread 8.679e-07 / 7.147e-07. Reading: this item's `-n 2` `[EX-34]` differs from Part B's at ≤ 1.3e-10 relative (symmetry 3.4126743589e-05 vs 3.4126695369e-05, 4.8e-11 abs), inside the band; the record stays Part B's. Census: figures ok=136 missing=0 broken=0 before and after; docrefs dead=0 guide=0 stale=125 exit=2 before and after. Stale guide lines named, not edited: `:104` ("reproduce the record to 2e-3", band is 1e-6), the §3 table digits (0.894543 …) predate `OPS-34`, `:80–81` ("printed, never asserted" for symmetry / spectral, asserted since `OPS-34`). The e2c branch is deleted after this commit.

### Outcome of item 587 (2026-10-22, 12:30 UTC slot)

DONE. `ans:2` honours `FEM_EM_ANS_SKIP_PRIVATE`; Row 1 is 3a's owner-cell (`cell`, route (ii)) family, the cell average of σ|E|²/(2ρ), imports `WIDTH_RTOL` / `ROUTE_IMAG_RTOL` / `PLANT` from 3a's module, records transcribed from `20261003T033225Z_POST-8-step3a-W1.log.gz:1926,1936–1939`. `_pointwise_sar` byte-identical (legacy CG1 row printed as a labelled record). Windows: `20261022T123323Z_POST-8-step3b-iii-ans2-n2.log` 251 s, `…123751Z_…ans2-n1.log` 386 s, `…124425Z_…ans2-n4.log` 195 s, consumer `…124753Z_…ans9-consumer.log` 130 s (7 passed, 2 skipped). (W) max miss 8.329e-11 / 1.838e-10 / 1.890e-10 at `-n 2 / 1 / 4` (prediction ≤ 1e-12 at `-n 2` was missed — `_build_mass_averaged` did not reproduce 3a's `case` bit for bit — printed, not asserted); (Own) exact; (Im) 0; CG1 min separation 1.190e-02 / 1.413e-02 / 2.022e-03; driven-point `_spread` owner-cell 26.148 % (all widths) vs CG1 13.294 / 3.124 / 15.730 %. Committed `metrics.json` / `COMPARISON.md` are the `-n 2` pair (ruling (3)), preserved after W1 and restored after W3. Census docrefs dead=0 guide=0 stale 125 → 124 (predicted unchanged; the edited `ans:2` files are now fresh), exit 2 both.

### Audit of the closure — 2026-10-22 10:30 review

**PASS** (`auditor`, closing commit `a7a39266`, with e4 at `5b634bc5`).
- **Every cited window is Status 0 with its elapsed time:** `20261022T123323Z_POST-8-step3b-iii-ans2-n2.log:1943–1944` (251 s), `…123751Z_…ans2-n1.log:1943–1944` (386 s), `…124425Z_…ans2-n4.log:1943–1944` (195 s), `…124753Z_…ans9-consumer.log:2076–2077` (130 s).
- **The digits trace:** (W) 8.329e-11 / 1.838e-10 / 1.890e-10 and (Own) exact at `:1902–1907` of each `ans2` log. This review re-read `n2.log:1907` and `:1940` ("ALL ANCHORS GREEN").
- **The anchor is executed:** asserts at `02_birdcage_coil_driven_sar_10MHz.py:1198–1210` against the bands imported from `test_post8_step3a_averaged_row1.py:88, 92, 94`. No tolerance moved; e4's four records were re-based `(1*)` from `…20261020T214423Z_…-measure-n2.log:677`.

Caveats, recorded and not demoting:
- The row's tier cell said "≈ 190–290 s" against 3b-iii's 386 s `-n 1` window; it is amended in this review.
- The script re-implements 3a's route (ii) as `_owner_cell_sar` rather than importing it. The outcome disclosed this: the ≤ 1e-12 prediction was missed, and the misses sit ≈ 5× inside 1e-9.
- Width invariance is read against the one `-n 2` record, not pairwise; the pairwise spread is ≤ 3.8e-10.
- The consumer `test_ans9_far_side_sar` only prints Row 1.

**The consumer's skip count, item 587's question.** It rose from 1 to 2 because `test_371b_thermal_half_on_degree2_exact_source` (`@needs_column`) was added by `50a5a43a` at 13:55 CDT on 10-11. That is after the 10-11 `d1` window (`20261011T184301Z_ANS-9-step1b-W1-d1.log`, 13:43 CDT), so the earlier count could not include it. This is a reading, not a defect.

### §9 item 598 — 2026-10-22 (ops / docs; follow-on to item 586)

Queued by the 2026-10-22 10:30 review. Executor: `implementer`. Smoke, no `mpiexec`, `-n 1`. ≈ 3 measured (ops / docs median) / 20 predicted min. Depends on nothing.

**What.** Item 586's outcome named three stale lines in `examples/ports/01_two_torus_port_pair.md` and left them unedited, per the 03:00 ruling. This item corrects them from the `-n 2` e4 log and changes nothing else. **Guide only: no edit to the script, to any record or to any band.**
- (G1) The §3 ladder table (rows raw / + PEC box / + gap physics, `:68–70`) is re-written from `20261022T110209Z_POST-8-step3b-ii-e4-ports1-n2.log`. The ratio and deviation of each row are copied from the log's printed ladder: `:670` for the PEC-box rung, and `:684` for raw 0.894636 and corrected 0.939944 (−6.01 %). The raw row's deviation and the PEC-box exponent text are also copied from the log, never recomputed by hand. The current digits (0.894543 …) predate `OPS-34`.
- (G2) `:80–81`: the symmetry / spectral bullet now says "asserted" (since `OPS-34`, against `REPRODUCTION_BAND_SYMMETRY_ABSOLUTE` / `REPRODUCTION_BAND_RELATIVE` = 1e-6). The digits are re-pointed to `:684` (3.4127e-05, 0.861324), and the source log cited is the e4 `-n 2` log. The old v0.7.2 / EX-30 sentence stays as history, labelled superseded.
- (G3) `:104–105`: "raw and corrected mutual reproduce the record to 2e-3" becomes "… to `REPRODUCTION_BAND_RELATIVE` = 1e-6 (the record re-based `(1*)` onto Part B's `-n 2` digits, item 586)".

**Anchor.** Every edited figure must equal its source string. A gitignored stdlib re-read helper (`logs/post8_598_reread.py`, the H31 Part A pattern of `logs/wf10_h31_reread.py`, with its source cat'ed into the log) does the check. For each site it prints a triple (guide string, source, hit). The sources are the cited log line and the script's band constants, parsed with `ast` from `examples/ports/01_two_torus_port_pair.py:203–210`. The helper exits 1 on any miss. Window `scripts/testing/run_and_log.sh POST-8-598-reread "timeout -k 30 30 python3 logs/post8_598_reread.py"`. If a bare `python3` is denied on the host, run the same call inside the container under `docker compose exec -T fem-em-solver bash -lc 'cd /workspace && …'`.

**Negative control (asserted as built).** The same helper, run with `--guide-rev 46aa9862` (the pre-edit guide read through `git show`), must report a miss at each of G1–G3 and exit 1. Every site that hits on the edited guide must miss on the old one.

**Windows.** Docrefs pre and post, with the census command item 586 used (`…110107Z` / `…111047Z`). Expect `dead=0 guide=0` both times; the stale count is printed.

**Traps.** `≈`, never `~`. The guide is text: no edit to `01_two_torus_port_pair.py`, so no example re-run. If the script itself turns out to need a change, park the item: it would then need rule (i)'s re-run. No AED figure is involved. The harness call is the whole Bash command.

**Scope.** Text only. `POST-8` stays ✅, the record stays Part B's, and no band moves.

**Negative result.** If a site's source digit is not printed in the e4 log, leave that cell unchanged and name it in the outcome. Never compute a digit by hand.

**Status it can move.** The three stale sites named in item 586's outcome close. The guide stops contradicting the asserted 1e-6 bands. The docrefs stale count for the pair is printed.

### Outcome of item 598 (2026-10-22, 21:30 UTC slot)

DONE. Guide only (`examples/ports/01_two_torus_port_pair.md`); no script, record or band touched, so no example re-run. (G1) §3 table rows now raw 0.894636 (−10.54 %), + PEC box 0.911536 (−8.85 %, `D∞ = +0.0169`, `p = 1.657`), + gap physics 0.939944 (−6.01 %, `÷(1 − 0.030224)`), from `20261022T110209Z_POST-8-step3b-ii-e4-ports1-n2.log:669–671` (the raw deviation and the PEC-box exponent sit on `:669` / `:670`, not `:684`; `:684` carries raw / corrected and is re-read too), with a source sentence under the table naming the pre-`OPS-34` digits. (G2) the symmetry / spectral bullet reads 3.4127e-05 / 0.861324 from `:684`, "asserted" against `REPRODUCTION_BAND_SYMMETRY_ABSOLUTE` / `REPRODUCTION_BAND_RELATIVE` = 1e-6; the `EX-30` / v0.7.2 sentence kept, labelled superseded. (G3) "to 2e-3" → "to `REPRODUCTION_BAND_RELATIVE` = 1e-6 (… item 586)". Anchor: `logs/post8_598_reread.py` (gitignored, source cat'ed into the log) 10 / 10 hit — G1 5, G2 4, G3 1; bands parsed with `ast` at script `:203` / `:204` — `20261022T213157Z_POST-8-598-reread.log` (exit 0, 0 s). Negative control asserted as built: `--guide-rev 46aa9862` 10 / 10 miss, exit 1 — `20261022T213204Z_POST-8-598-reread-negctl.log`. Docrefs `dead=0 guide=0 stale=125 exit=2` pre (`20261022T213057Z_POST-8-598-pre-docrefs.log`) and post (`20261022T213205Z_POST-8-598-post-docrefs.log`); the `ports:1` pair is not among the stale. **Further stale sites, named and not edited (outside G1–G3):** `:79` reciprocity 5.8343e-04 (the e4 log prints 7.8052e-04 at `:667`); `:100` quadrature 3.9111e-04 / 1.4044e-04 (e4 `:657` / `:659`: 2.6056e-04 / 2.3522e-04); `:61–62` cost / cell count (178 055 cells, 134 s; e4 `:685`: 177 998 cells, 146.8 s).

### §9 item 600 — 2026-10-22 (ops / docs; follow-on to items 586 and 598)

*Written by the 2026-10-22 18:00 daily review from item 598's outcome note ("three more stale `ports:1` guide sites outside G1–G3").*

**Licence.** Ops / docs: a text-only correction of a guide that contradicts its own example's newest log. Not a chain step. `POST-8` stays ✅. **No band, tolerance, record, Done-when or target is set or moved; no script is edited, so the example is not re-run** (rule (i) is not engaged).

**What it is.** In `examples/ports/01_two_torus_port_pair.md` (line numbers as of `8f89f85a`, after item 598's edit), re-write the three sites item 598 named from `docs/testing/logs/20261022T110209Z_POST-8-step3b-ii-e4-ports1-n2.log`:
- **(G4) cost and cells, `:60–62`:** "**134 s** — mesh 36.1 s (178 055 cells), two solves 22.0 s + 22.5 s (`…EX-18-example-n2-v3.log`, exit 0)" → the e4 `:685` figures, "146.8 s on the script's own timer — mesh 28.4 s (177 998 cells), two solves 24.4 s + 24.3 s", citing `…e4-ports1-n2.log:685`. Keep one clause saying the earlier 134 s / 178 055 figures predate `OPS-34` (the pattern item 598 used under the §3 table).
- **(G5) reciprocity, `:84`:** 5.8343e-04 → **7.8052e-04** (`:667`), still "printed", with the same parenthesis.
- **(G6) quadrature, `:109–113`:** the undriven port's measured agreement "3.9111e-04 / 1.4044e-04" → **2.6056e-04 / 2.3522e-04** (`:657` / `:659`), and the driven diagonal's "2.3e-2 / 3.5e-2" → **3.4e-3 / 6.8e-3** (`:656` 3.4087e-03 / `:660` 6.8082e-03, written at two significant figures, as the guide does today). The disposition sentence ("printed and not gated — step 3b-x's standing disposition …") stays word for word. **Note for the executor:** at e4 the driven diagonal's order-2049 → 4097 change is 3.4e-3 / 6.8e-3, still above the 1e-3 precondition. "Does not converge in the quadrature" therefore still describes what is printed; do not re-word it.

**Anchor (asserted).** Extend item 598's gitignored helper pattern: a new `logs/post8_600_reread.py` (stdlib only; its source is `cat`'ed into the tracked log by the harness command, as 598's was), with triples (site, guide string, log line, log string). It exits 1 on any miss, where a hit needs the guide string present, the log string present on the cited line, and the numbers equal token for token. Predicted **≥ 7 / 7 hit**: G4 4 (146.8, 28.4, 177 998 written as `177 998` ↔ `177998`, 24.4 + 24.3), G5 1, G6 2 (one per log line pair, at the guide's precision). The helper's own normaliser drops thin / regular spaces inside integers. State the count printed.

**Negative control (asserted, by construction).** `--guide-rev 8f89f85a` reads the pre-edit guide and must miss every triple, exit 1. Backing: item 598's control on the same helper pattern, `20261022T213204Z_POST-8-598-reread-negctl.log`.

**Tier / cost.** Smoke; no container solve; three harness calls (`python3` through `run_and_log.sh`, as 598). Class ops / docs: **≈ 3 measured / ≈ 8 predicted** slot-minutes. Executor: the slot itself (text only, no spawn). Depends on nothing. Daytime: no. Family: `POST-8` follow-on; attempt 1 of 4.

**Commands** (item 598's, relabelled; the harness call is the whole Bash command, bare, repo-relative, nothing chained or piped; Bash `timeout: 660000`):
```
scripts/testing/run_and_log.sh POST-8-600-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"
scripts/testing/run_and_log.sh POST-8-600-reread "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && cat logs/post8_600_reread.py && timeout -k 30 30 python3 logs/post8_600_reread.py'"
scripts/testing/run_and_log.sh POST-8-600-reread-negctl "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 30 python3 logs/post8_600_reread.py --guide-rev 8f89f85a'"
scripts/testing/run_and_log.sh POST-8-600-post-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"
```
These are item 598's landed commands (`20261022T213057Z_…pre-docrefs.log:12`, `…213157Z_…reread.log:12`), relabelled. If `--guide-rev` needs `git` inside the container and it is absent there, read the pinned guide the way 598's negative control did (`…213204Z_…negctl.log:12`) and say so. Docrefs must read `dead=0` pre and post (exit 0 or 2; staleness is information).

**Traps already paid for.** `≈`, never a tilde. No AED string. The guide's `(1*)` sentence and the G1–G3 lines item 598 wrote are not touched. The helper lives under gitignored `logs/`; never commit it, and `git status --short` must be clean apart from the guide, the logs and the plan files. Do not re-run `ports:1` "to refresh" a digit: the e4 log is the source.

**Scope.** Covered: the three remaining stale sites. Not claimed: any new reading, any record move, the `EX-18` example-run timing (e4 is a pytest window; the guide says "the script's own timer").

**Negative result.** A miss after the edit means the guide string or the cited line differs. Fix once in-letter; if it still misses, revert the guide, journal the triple and stop. No BLOCKED needed unless the e4 log is absent.

**Status it can move.** Item 598's "further stale sites" list closes; the `ports:1` guide carries no pre-e4 digit outside labelled history. `POST-8` stays ✅.

### Outcome of item 600 (2026-10-23, 00:30 UTC slot)

DONE. Guide only (`examples/ports/01_two_torus_port_pair.md`); no script, record or band touched, so no example re-run. (G4) cost → "146.8 s on the script's own timer — mesh 28.4 s (177 998 cells), two solves 24.4 s + 24.3 s" from e4 `:685`, with a clause naming the earlier 134 s / 178 055 cells (`…EX-18-example-n2-v3.log`) as pre-`OPS-34`. (G5) reciprocity 7.8052e-04 from `:667`, still "printed". (G6) undriven quadrature 2.6056e-04 / 2.3522e-04 (`:657` / `:659`); driven diagonal 3.4e-3 / 6.8e-3 (`:656` / `:660`, two significant figures); the disposition sentence is unchanged. Anchor: `logs/post8_600_reread.py` (gitignored; source cat'ed into the log) 7 / 7 hit (G4 4, G5 1, G6 2) as predicted, exit 0, in `20261023T003117Z_POST-8-600-reread.log`. Negative control `--guide-rev 8f89f85a`: 7 / 7 miss, exit 1, in `20261023T003124Z_POST-8-600-reread-negctl.log` (git ran inside the container). Docrefs `dead=0 guide=0 stale=125 exit=2` pre (`20261023T003035Z_POST-8-600-pre-docrefs.log`) and post (`20261023T003132Z_POST-8-600-post-docrefs.log`). Item 598's "further stale sites" list is closed.
