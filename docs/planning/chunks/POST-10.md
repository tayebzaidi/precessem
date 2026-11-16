# POST-10 — The point-owner rule on a `create_submesh` mesh: key submesh cells by their parent's original cell index

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 254 — 2026-10-05 (lane 1, known-issues "Resolves with")

*Drafted by a general-purpose agent for the 2026-10-05 10:30 daily review; the review reviews it.*

**`POST-10` step 1: make `owner="lowest_original_cell"` work on a `create_submesh` mesh by keying each submesh cell with its parent cell's original index.** This is the "Resolves with" of the known-issues entry *2026-10-04 — 🟡 OPEN (`ANS-9` step 1a, §9 item 218) — `evaluate_vector_field_parallel(owner="lowest_original_cell")` raises `RuntimeError: Mixed topology unsupported.` on a `create_submesh` mesh* (`docs/testing/known-issues.md:1372–1379`): "A `post/evaluation.py` fix mapping submesh cells to parent original indices (a `POST` step, not yet opened), or a documented restriction." It is a lane-1 item (§9 rule (5): "open known-issues entries whose own 'Resolves with' a slot can execute", PROJECT_PLAN.md:3742–3743).
- **Executor:** implementer. One `src/` change and one new test module. No solve, no meshing beyond `create_unit_cube`.
- **Tier:** standard, three windows: W1 real build `-n 2`, W2 complex `-n 4`, W3 complex `-n 1`. Each runs `timeout -k 30 180`.
- **Class:** lane 1, ≈ 7 measured min (class median). **Predicted ≈ 40 min:** 15 fixed, ≈ 15 code (≈ 40 lines in `src/`, ≈ 220 in the test), three windows of ≈ 1 min each, and the remainder for the journal and the known-issues retirement.
- **Depends on nothing.**
- **Not a ⚠️ subsystem.** No `POST` row in §7 carries ⚠️. `POST-8`, whose rule this extends, is 🟡 and its rule has gates (PROJECT_PLAN.md:2567–2575). So §9's no-new-features rule (PROJECT_PLAN.md:3651) does not bind.

**What raises, read off the code and the log.**
- `point_owner_cells` keys the rule with `original = _original_cell_index(mesh)` (`src/fem_em_solver/post/evaluation.py:92`). That reads `mesh.topology.original_cell_index` (`:19–28`).
- On a `create_submesh` mesh, DolfinX 0.11.0.post0 (the image's version, `20261003T020658Z_POST-1-step7a-W2.log:97`) raises `RuntimeError: Mixed topology unsupported.` from `dolfinx/mesh.py:268` (`20261004T142026Z_ANS-9-step1a-W1.log:1949–1977`). A submesh carries no original cell index.
- **The raise is symmetric today.** It comes after `evaluate_vector_field_parallel`'s key allgather (`:216`) and before `point_owner_cells`' first collective (`:114`), and every rank's submesh topology refuses alike. So it does not hang, but the fix must keep it that way.
- The submesh does not reference its parent. **No zero-argument fix exists**, so the caller must pass the map. Every submesh site in the repo already holds it: `parent_ids = emap.sub_topology_to_topology(np.arange(n_sub, dtype=np.int32), False)` over owned + ghost sub cells (`tests/validation/test_mat10_sar_to_heat_sphere.py:74–88`, re-used as `_sphere_submesh` by `MAT-12` / `MAT-14` / `MAT-19` / `MAT-5`, and as `_tag_submesh` in `tests/validation/test_mat15_coil_driven_heating_phantom.py:105–120`). This idiom is green on 0.11.
- **Why the rule's meaning survives.** On the parent, the original cell index is the mesher's order, not the partition's (`POST-8` step 2, `docs/planning/chunks/POST-8.md:28`). Its restriction to the kept cells is injective. Over all ranks, the owned sub cells are exactly the kept cells, each counted once. So "lowest parent original index among the submesh cells containing the point" is one deterministic cell at every width. It is the same cell the parent's one-sided rule picks with `candidate_cells = keep` (`POST-1` step 7a, `:70–74`). **The argument does not need submesh and parent ownership to coincide**, so that is printed, not asserted.

**Code (`src/fem_em_solver/post/evaluation.py` only).**
- **New keyword `parent=None` on `point_owner_cells` and `evaluate_vector_field_parallel`.** It takes `None` or `(parent_mesh, parent_cells)`. `parent_cells` gives the parent-local cell index of every local (owned + ghost) cell of `mesh`, i.e. exactly the `parent_ids` above. Keys become `_original_cell_index(parent_mesh)[parent_cells]`, and nothing else in the rule changes: same owned-only candidates, same `candidate_cells` mask, same Allreduce MIN, same claim check.
- **One helper replaces the bare read at `:92`, collectively checked.** Compute the keys or a local error. The errors are: `original_cell_index` raised and `parent is None`; `parent_cells` has length ≠ `size_local + num_ghosts`; or `parent_cells` holds an index outside `[0, parent n_local)`. Then **one `allreduce` of a status code** runs, and the helper raises on **every** rank if any rank failed.
  - A missing index raises `RuntimeError`. Its message names `create_submesh`, `parent=(parent_mesh, parent_cells)` and `POST-10`, and chains the DolfinX error.
  - Malformed `parent_cells` raises `ValueError`.
  - This adds one integer collective to the existing owner path. It moves no value.
- **In `evaluate_vector_field_parallel`:**
  - `parent` is valid only with `owner="lowest_original_cell"`. Otherwise raise `ValueError` **before any collective**, mirroring `candidate_cells` (`:188–192`).
  - `parent is None` joins the OPS-40/55 key, which becomes `(owner, candidate_cells is None, parent is None)` (`:216–218`). The mismatch message at `:252–257` mentions it, and keeps the substrings "owner" and "collective" that K5 asserts (`tests/post/test_point_owner_rule.py:212`).
  - Docstrings for both functions gain the keyword. `OWNER_RULES` does not change. `post/__init__.py` exports nothing new.
- **With `parent=None` on a mesh that has an original index, the output must be bitwise unchanged.** The regression modules below show it.
- **No helper that builds `parent_cells`** (scope). The three test copies of the idiom stay as they are.
- **Comment-only edit, same commit:** `src/fem_em_solver/thermal/probe.py:9–13` cites the entry as live. Re-point it: "raises on `create_submesh` meshes unless `parent=` is given (`POST-10`)". The ball average stays integral-based, and `MAT-19`'s code does not change.

**Test: new `tests/post/test_submesh_owner_rule.py`, runs in both builds (as `test_point_owner_rule.py` does).**
- **Fixture.** The parent is `create_unit_cube(comm, 4, 4, 4)`: 384 cells, the `POST-8` K fixture (`tests/post/test_point_owner_rule.py:52–53`).
- **Kept cells.** `keep` is every local cell (owned + ghost) whose midpoint has x < 0.75, from `dmesh.compute_midpoints(msh, 3, np.arange(n_local, dtype=np.int32))`. A tet midpoint never lies on the grid plane x = 0.75, so the split is geometric and the same at every width.
- **The submesh** is `sub, emap, _, _ = dmesh.create_submesh(msh, 3, keep)`, with `parent_cells` taken from the idiom verbatim (`test_mat10_sar_to_heat_sphere.py:82–88`).
- **Points** (`POINTS`, shape (6, 3), identical on every rank):
  - P0 interior `(0.13, 0.27, 0.61)`.
  - P1 on the facet plane x = 0.5: `(0.5, 0.3, 0.6)`.
  - P2 on the interior vertex: `(0.5, 0.5, 0.5)`. P0–P2 are `POST-8`'s three points.
  - P3 on the submesh-boundary facet: `(0.75, 0.3, 0.6)`.
  - P4 on the submesh-boundary vertex: `(0.75, 0.5, 0.5)`.
  - P5 outside the submesh: `(0.9, 0.3, 0.6)`.
- **Planted fields.**
  - Parent: a scalar DG0 holding `_original(msh)` on owned cells, then `x.scatter_forward()` (K3's idiom, `:135–142`).
  - Submesh: a scalar DG0 holding `_original(msh)[parent_cells[c]]` on owned sub cells, then `scatter_forward()`.
  - Submesh: an N1curl degree-1 field interpolating `A_VEC + B_VEC × x` (K2's field and constants, `:46–47`, `:124–128`). N1curl represents it exactly.
- **Oracle.** The geometric oracle independent of DolfinX's collision search (`test_point_owner_rule.py:62–82`, `BARY_TOL = 1e-12`) computes each point's incident set of parent original indices, restricted to kept cells. It allgathers each owned parent cell's original index, its four vertices and its `keep` flag. Rank 0 decides and broadcasts.

**Anchors (asserted at every width).**
- **(S0) Fixture identity (arithmetic).**
  - The submesh `size_global` is 288 (= 6 · 4³ · 3/4: three of four x-slabs).
  - `len(parent_cells) == sub size_local + num_ghosts`.
  - The kept cells' parent original indices are unique across ranks and their count is 288.
  - Printed: the fraction of owned sub cells whose parent is owned on the same rank (predicted 1.0, not asserted; the rule does not need it).
- **(K) Closed form against the oracle.**
  - `point_owner_cells(sub, POINTS, comm, parent=(msh, parent_cells))[0]` equals min(oracle keep-set) exactly (integers) at P0–P4, and −1 at P5.
  - Through `evaluate_vector_field_parallel(..., owner="lowest_original_cell", parent=...)`, the submesh DG0 value equals `float(min)` with imaginary part 0, P0–P4 are valid, and P5 is invalid.
  - Each point has exactly one claiming rank (K3's claims check, `:178–180`).
- **(I) The identity with the parent's one-sided rule.**
  - The owner arrays are equal as integers: submesh rule with `parent=` vs `point_owner_cells(msh, POINTS, comm, candidate_cells=keep)`.
  - `values` and `valid` are **bitwise** equal (`np.array_equal`, band 0): submesh DG0 under the rule with `parent=` vs parent DG0 under `owner="lowest_original_cell", candidate_cells=keep`.
  - DG0's basis is 1, so the value is the dof itself, copied from the same parent array. No round-off band applies.
- **(N) N1curl closed form on the submesh under the rule.** max |eval − (A + B × p)| ≤ `CLOSED_FORM_TOL` = 1e-12 at P0–P4.
  - Define the constant with `POST-8`'s citation comment, as `test_point_owner_rule.py:44–45` does ("the guard module's bound, `test_evaluation_collective_guard.py:83`").
  - (N) shows the chosen sub cell is evaluated correctly. It does not depend on which cell is chosen.
- **(W) Partition independence.**
  - (K) and (I) are asserted at `-n 1`, `-n 2` and `-n 4`, so the owner integers equal the same oracle minimum at every width.
  - Within each run, `owner_original`, the DG0 and N1curl values are bitwise equal on every rank (`_bitwise_equal_across_ranks`, `:84–86`).
  - Print the owner integers at every width so the three logs can be diffed by eye.

**Negative controls (rule (e): labels as stated; ceiling first).**
- **(C1) Discrimination, asserted (arithmetic).** At P1, P2 and P4 the oracle keep-set has ≥ 2 members, with max − min ≥ 1.
  - **Ceiling:** the most any wrong-cell choice can be separated from the rule is exactly that set's max − min. It is an integer ≥ 1 against (K) / (I)'s band of 0, so the fixture cannot be blind to a different owner.
  - **Predicted, printed only:** set sizes 2 / 24 / ≈ 12 at P1 / P2 / P4, and 1 at P3. P1 and P2 lie inside the kept region, so their keep-sets are `POST-8`'s measured incident sets on the same cube (`20261003T020808Z_POST-1-step7a-W3-real.log:104–106`): P0 {217}; P1 {223, 232}, separation 9; P2 24 cells, min 126 / max 257, separation 131. **Predicted rule owners:** P0 217, P1 223, P2 126. They are a prediction, not a band: (K) asserts the oracle's own minimum.
- **(C2) The restriction is explicit, asserted.** The submesh call with `owner="lowest_original_cell"` and **no** `parent` raises `RuntimeError` on **every** rank, at every width.
  - The flags are allgathered as K5 does (`:192–212`). The message contains `parent=` and `create_submesh`.
  - A silent fallback to submesh-local numbering fails here.
- **(C3) Asserted.** `parent=` with `owner=None` raises `ValueError` before any collective (the `candidate_cells` analogue, `tests/post/test_one_sided_owner_rule.py:135–141`).
- **(C4) Asserted.** `parent_cells[:-1]`, one entry short, raises `ValueError` on every rank (allgathered flags). It must not hang.
- **(C5) The naive fix, printed and predicted, not asserted.** Key each submesh cell by `sub.topology.index_map(3).local_to_global` instead of the parent map: Allreduce MIN, then allgather (key, DG0 value) pairs. Print the DG0 value that key selects at P1, P2 and P4 beside the rule's value.
  - **Predicted:** differs from the rule at ≥ 1 point at ≥ 1 width. The index-map global number is partition-dependent (`POST-8.md:28`).
  - Not asserted, because nothing forces a disagreement on this fixture.
- **Printed only:** at P3 / P4, the *unrestricted* parent owner (`point_owner_cells(msh, POINTS)`), beside the submesh owner. It shows whether the parent's rule would pick a cell on the x > 0.75 side.

**Regression (rule (c)), same windows, must be green and unedited.**
- `tests/post/test_evaluation_collective_guard.py` covers the legacy `owner=None` path and the OPS-40/55 key.
- `tests/post/test_point_owner_rule.py` covers K0–K6, the unrestricted rule and the mismatch message.
- `tests/post/test_one_sided_owner_rule.py` covers `candidate_cells`.
- `tests/post/test_post1_one_sided_interface_sphere.py` covers the complex-build gmsh one-sided gate.
- **Not re-run, with reason:** `tests/validation/test_post8_owner_rule_ans2.py`, `test_post8_step3a_averaged_row1.py`, `test_post1_one_sided_interface_birdcage.py` and `test_imp3_palace_fields.py`. They call the same rule with `parent=None`. Their value path is the one the four modules above exercise. The ANS-2 module alone is a 590 s-ceiling heavy window (`20261002T200654Z_POST-8-step2-W2.log`, 191 s with it).
- Name them in the journal as untouched.

**Windows** (foreground, Bash `timeout: 660000`; each harness call is the whole Bash command, bare and repo-relative; `-s`; no pipes; no `-k` selection; `tests/environment` first).

`scripts/testing/run_and_log.sh POST-10-step1-W1-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/post/test_submesh_owner_rule.py tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py tests/post/test_one_sided_owner_rule.py -v -s --tb=short'"`

`scripts/testing/run_and_log.sh POST-10-step1-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 4 python3 -m pytest tests/environment tests/post/test_submesh_owner_rule.py tests/post/test_one_sided_owner_rule.py tests/post/test_post1_one_sided_interface_sphere.py -v -s --tb=short'"`

`scripts/testing/run_and_log.sh POST-10-step1-W3 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment tests/post/test_submesh_owner_rule.py tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py -v -s --tb=short'"`

**Cost, from logs.**
- W1's regression set (environment + guard + owner + one-sided, real, `-n 2`) ran in 53 s (`20261003T020808Z_POST-1-step7a-W3-real.log`).
- W2's (environment + one-sided + sphere, complex, `-n 4`) ran in 67 s (`20261003T020658Z_POST-1-step7a-W2.log`).
- Environment plus light modules at `-n 1` ran in 51 s (`20261002T201726Z_PORT-3-W2.log`).
- The new module is a 384-cell cube with a few collectives: ≈ 5–10 s.
- **Predicted:** W1 ≈ 60 s, W2 ≈ 75 s, W3 ≈ 60 s, all under the 180 s ceiling.
- **Box-speed gate:** a window over 120 s is a slow box. That is a price, not a defect. A window killed at 180 s stops the item BLOCKED with the reading; never raise the timeout.

**Traps already paid for.**
- **Collectives.** Every rank runs every collective, and asserts come after them. A `comm.` call inside `if comm.rank == 0:` deadlocks (item 227).
  - The must-raise tests (C2, C4) allgather their flags **after** the call (K5 idiom).
  - The helper's status allreduce runs before any raise on **both** the error and the success path. A rank-local raise ahead of the Allreduce at `:114` would hang the others.
- **Use `owner_original` to compare.** `owner_local_cell` is −1 off the owning rank, so compare `owner_original`, never local cells.
- **Never key on `local_to_global`.** The index-map global number is partition-dependent (`POST-8.md:28`). It appears only in the printed (C5).
- **Rank-local data.** The DG0 fields are set on owned cells, then `scatter_forward()`. `keep` is computed over owned + ghost cells, as the idiom requires. The oracle gathers owned cells only, so nothing is counted twice.
- **`create_submesh` on 0.11** returns a 4-tuple, and its `EntityMap` maps with `sub_topology_to_topology(cells, False)`. Copy the idiom; do not re-derive it.
- **Complex build.** DG0 values are complex128 with imaginary part 0. `np.array_equal` compares complex arrays exactly. Take `.real` only for printing and for (K)'s integer comparison, after asserting imaginary part 0.
- **Harness.** `timeout -k 30`; no `run_in_background`; `-k a or b` splits inside a quoted container command; pytest captures prints without `-s`. After a kill, run `pgrep -c python3` in the container first. Sweep orphan `.c` stems only with the rubric's one-line sweep, and only when that count reads 0.
- **Records.** Do not edit `examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/09_coil_driven_gel_heating_10MHz.py`. θ is continuous P1, so its `owner=None` call (`:286–291`) gives the same value. Switching it would re-run a run of record for nothing.

**Scope.**
- Point evaluation on a single-level `create_submesh` of an affine tetrahedral mesh, with the caller-supplied parent map.
- **Not covered:**
  - a submesh of a submesh: the parent itself lacks an index, so it raises with the new message;
  - `candidate_cells` combined with `parent=`: allowed by the code, not gated here;
  - mixed-topology meshes;
  - a library helper that builds `parent_cells`;
  - any call-site migration. `MAT-12` / `MAT-15` / `MAT-19` / `ANS-9` keep their current routes: integration, or `owner=None` on continuous fields.
- No band, record, example or existing test changes. `POST-8` stays 🟡; its open steps are untouched.

**Negative result.**
- **(K) or (I) red with (S0) and (C1)–(C4) green:** the parent-keyed rule is not the parent's one-sided rule.
  - Park on `attempt/POST-10-step1-<UTC>` and mark the item BLOCKED in the same commit.
  - Append the reading to the known-issues entry, which stays open. Do **not** land the clearer error alone as the "documented restriction" resolution. That is the review's call.
- **(S0) red:** the fixture or the `EntityMap` API drifted. Stop and report.
- **(C2) not raising:** the image now gives submeshes an original cell index. Record what it returns, beside the parent's, and stop: the review rules whether the rule may use it.
- **Any regression module red:** the change moved existing behaviour. Revert the `src/` change, park, and stop BLOCKED. Never edit the old test.
- **(C3) / (C4) or a hang:** the instrument or the collective order is wrong. Fix it once in-slot, disclose it, and re-run the window. If it is still red, BLOCKED.
- There is no PREDICTED RED line; the drafter predicts green.

**Status it can move.**
- `POST-10` ⬜ → ✅ on all three windows green: §4 needs (K) and (I) as exact identities, (W) at three widths, and elapsed times recorded.
- In the fixing commit, the known-issues entry `2026-10-04 — … Mixed topology unsupported …` (`docs/testing/known-issues.md:1372–1379`) is **retired** (its own instruction, "Remove this entry in the commit that fixes it"). Its full text goes to `known-issues-archive.md` and a `RETIRED 2026-10-05` bullet goes under "Retired entries", per the file's convention.
- The `thermal/probe.py` docstring pointer is re-pointed in the same commit.

**§9 line 2:** → docs/planning/chunks/POST-10.md#§9-item-254

**Outcome (2026-10-05 15:00 CDT slot): DONE, all three windows green on first run; `POST-10` ⬜ → ✅.** W1 real `-n 2` `20261005T200223Z_POST-10-step1-W1-real.log` (64 s), W2 complex `-n 4` `20261005T200334Z_POST-10-step1-W2.log` (80 s), W3 complex `-n 1` `20261005T200500Z_POST-10-step1-W3.log` (67 s). Owner integers P0–P5 [217, 223, 126, 229, 132, −1] equal the oracle minimum at every width (predictions P0 217 / P1 223 / P2 126 met; P3 set size 1, P4 12); (I) bitwise vs `candidate_cells=keep`; (N) ≤ 2.220e-15; (C1) separations 9 / 131 / 121; (C2) / (C3) / (C4) raise on every rank; (C5) the naive key differs at ≥ 1 point at every width (printed). Same-rank ownership 1.0000 (printed). Regression modules green, unedited. Known-issues entry retired in the same commit.

## §7 row moved verbatim 2026-10-05 (rotate_plan_archive.py chunks)

 The point-owner rule on a `create_submesh` mesh: key submesh cells by their parent's original cell index *(history and §9 briefs: `docs/planning/chunks/POST-10.md`)* | ✅ **Done 2026-10-05 15:00 CDT slot (§9 item 254).** `point_owner_cells` / `evaluate_vector_field_parallel` take `parent=(parent_mesh, parent_cells)`; keys are the parent's original index, one collectively checked status allreduce. On the 288-cell cut, owner integers [217, 223, 126, 229, 132, −1] (P0–P5) equal the oracle minimum at `-n 1 / 2 / 4`; (I) owner integers and DG0 values / validity bitwise equal to the parent's `candidate_cells=keep` rule; (N) N1curl ≤ 2.220e-15; (C1) separations 9 / 131 / 121; (C2) / (C3) / (C4) raise on every rank; regression modules green unedited. W1 real `-n 2` 64 s (`20261005T200223Z_POST-10-step1-W1-real.log`), W2 complex `-n 4` 80 s (`…200334Z_POST-10-step1-W2.log`), W3 complex `-n 1` 67 s (`…200500Z_POST-10-step1-W3.log`). Known-issues entry retired in the same commit. Opened 2026-10-05 10:30 review (lane 1) from known-issues 2026-10-04 (`ANS-9` step 1a): `owner="lowest_original_cell"` raises `Mixed topology unsupported.` on any submesh (`20261004T142026Z_ANS-9-step1a-W1.log:1949–1977`). Queued §9 item 254: a `parent=(parent_mesh, parent_cells)` keyword on `post/evaluation.py`, with a collectively checked raise. **Done-when (§4):** on a 4³ cube cut at x < 0.75 (288 cells), at P0–P5 (interior, facet, vertex, boundary facet, boundary vertex, outside): (K) the submesh owner equals the geometric oracle's minimum over the kept cells, exactly; (I) values and validity are bitwise equal to the parent's `candidate_cells=keep` rule on a planted DG0; (N) N1curl `a + b×x` meets its closed form ≤ 1e-12; (W) all of these at `-n 1/2/4`. Negative controls: (C1) the planted field separates a wrong owner by ≥ 1 at facet and vertex points; (C2) a submesh call without `parent` raises on every rank; (C3) / (C4) misuse raises. The `POST-8` / `POST-1` modules re-run green and unedited (rule (c)). The fixing commit retires the known-issues entry ⇒ ✅. 
