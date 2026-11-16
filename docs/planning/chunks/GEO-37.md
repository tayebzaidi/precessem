# GEO-37 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 222 — 2026-10-04 (chain U6)

*Drafted by a general-purpose agent for the 2026-10-04 03:00 daily review; the review reviews it.*

**`GEO-37` step 1: the 2D-facet probe.** The fixture is a unit box with one rectangle embedded in its interior, tagged as a gmsh dim-2 physical group. It is read through the project's own gmsh interop at `-n 1 / 2 / 4` in one logged window. The question is whether the interior rectangle arrives in DolfinX's `facet_tags` on the 0.11 image complete, unduplicated and without a hang.
- **Executor: implementer, not `mesh-probe`.** The §7 row and §10 U6 both say `mesh-probe`. But `mesh-probe`'s contract is "a table + a probe script, never an assertion". The row's Done-when needs the tagged area "against the rectangle's closed form at 1e-9" and the count "equal across widths after reduction" to move ⬜ → ✅, and §4 needs at least one quantitative assertion. A `mesh-probe` run therefore cannot move the row by construction (rubric element 7).
  - The deliverable keeps the probe's shape: a probe script plus a table. An implementer adds the asserting test that reads the probe's records.
  - **The review amends the row's executor word in the same commit (rule (j)).**
- **Tier:** standard ceiling (the step timeouts sum to 180 s), real build, `-n 1 / 2 / 4`. The row says "smoke by prediction". The window is predicted at 25–50 s (three launches of ≈ 6–12 s plus pytest ≈ 5 s), so the tier word is re-labelled on measurement.
- **Class:** chain step, no solve. The class median is ≈ 11 min.
- **Cost:** ≈ 11 measured / ≈ 35 predicted min. That is 15 fixed + ≈ 15 code + one window of ≈ 1 min + spare.
- **Depends on nothing.** It is measurement only: **no `src/` change** (the row's letter).

**Why this is open (history paid for).**
- On the 0.7 image, `PORT-1` step 3b-iv's dim-2 groups `201` / `202` on interior facets timed out at `-n 2` (exit 124 at 181 s) and passed at `-n 1` (22.5 s). They were filed as known-issues 9 "hangs `model_to_mesh`".
- The same day's marker probe **half-refuted that**. `model_to_mesh` and facet creation completed at `-n 2` in 14 s (`20260805T034007Z_PORT-1-step3biv-hang-localise-fine.log`). What hung was the **`dS` facet-area assembly** on a mesh with no ghost cells: the two ranks' stacks differed (a mismatched collective).
- The project then rebuilt interfaces dolfinx-side from cell tags. That is `_interface_facet_tags` (`io/mesh.py:282`), whose docstring still carries the old diagnosis. Every current sheet, including the exact-rectangle birdcage ports (`io/mesh.py:≈5520–5553`, shared faces of two piece volumes), goes through that sliver-pair route.
- No fixture has tagged an interior face natively since. The 0.11 `dolfinx.io.gmsh.model_to_mesh`, reached through the shim `_model_to_mesh` (`io/mesh.py:13–28`), has never been asked.

**Code (tests only).**
- **`tests/mesh/probe_geo37_interior_facet.py`** (probe script; `python_files = test_*.py`, so it is never collected), taking `--out <json>`. The model is built on rank 0 inside `try` / `except BaseException`, then `_raise_geometry_failure_on_every_rank` before the collective (`GEO-23`, so a gmsh throw cannot hang the other ranks).
  - **CAD:**
    - `occ.addBox(0, 0, 0, 1, 1, 1)` and `occ.addRectangle(0.1, 0.2, 0.47, 0.8, 0.6)`. The rectangle is w × h = 0.8 × 0.6 (area **0.48**) in the plane z = 0.47; it spans most of the box so a four-way partition is likely to cut near it.
    - `occ.fragment([(3, box)], [(2, rect)])`, then `synchronize`.
    - Raise unless exactly one volume remains, and unless exactly one surface has bbox z-span < 1e-12 at z = 0.47 with `occ.getMass` within 1e-12 of 0.48.
  - **Physical groups:** `addPhysicalGroup(3, [vol], 1)`; `addPhysicalGroup(2, [rect], 201)`; `addPhysicalGroup(2, the six box faces, 101)`.
  - **Mesh:** `Mesh.MeshSizeMax = Mesh.MeshSizeMin = 0.08` gives a predicted 10–20 k tets and ≈ 150 rectangle triangles. Record `N_gmsh = len(gmsh.model.mesh.getElements(2, rect)[1][0])` and broadcast it.
  - **Import:** `mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, 0, gdim=3)` (the project shim, default partitioner), then `gmsh.finalize()`.
  - **Measure, rank-safely:**
    - `fdim = 2`; `n_owned = mesh.topology.index_map(2).size_local`.
    - For tag t ∈ {201, 101}, the local facets are `facet_tags.find(t)` and the owned ones are those `< n_owned`.
    - Triangle areas in numpy: `dolfinx.mesh.entities_to_geometry(mesh, 2, facets)` → `mesh.geometry.x` → ½‖(b − a) × (c − a)‖.
    - `count_t` and `area_t` are owned sums, `allreduce`d.
    - Also `count_local_t` / `area_local_t` (owned plus non-owned local copies, reduced; printed) and the per-rank owned and local counts (gathered and printed: does the partition touch the rectangle?).
    - **The planted control quantity:** on the lowest rank owning any 201 facet, drop its first owned 201 facet, giving `count_drop`, `area_drop` (reduced).
  - **Rank 0 writes JSON:** `comm_size`, `N_gmsh`, cells (reduced), all the above, `seconds_mesh`, `seconds_import`, `seconds_measure`, and a `completed: true` stamp written last.
  - **Never use a UFL `dS(201)` / `ds` integral here.** That is the assembly that hung in 2026-08 (no ghost cells; the default partitioner gives `cells_ghost = 0`). A numpy owned-facet sum measures the same area with no collective beyond `allreduce`.
- **`tests/mesh/test_geo37_interior_facet.py`** reads `GEO37_RECORD_DIR`. **The module is skipped when that variable is unset**, so ordinary suite runs never fail on a missing record. One test per assertion; it prints the table (width | cells | N_gmsh | count | area − 0.48 | local-copy excess | per-rank owned | seconds) and `[GEO-37] VERDICT native interior facet tags on 0.11: GREEN|RED`.

**Anchors (asserted, per width w ∈ {1, 2, 4}).**
- **(A1) No hang:** `n<w>.json` exists, carries `completed: true`, and `comm_size == w`.
- **(A2) Count:** `count_201` is equal across the three widths **and** equal to `N_gmsh`. A route that drops or duplicates partition-boundary facets reads ≠ `N_gmsh` at `-n 2` / `-n 4` only.
- **(A3) Area:** |`area_201`/0.48 − 1| ≤ **1e-9**. The closed form is w · h; a planar rectangle tiled by linear triangles is exact, and the predicted reading is ≈ 1e-15.
- **(A4) Companion, the instrument check:** |`area_101`/6 − 1| ≤ 1e-9, and `count_101` equal across widths. This shows the owned-facet arithmetic on the exterior group, which has always reached the tags since `GEO-10`. The case where (A4) is green and (A3) red is exactly the interior-specific failure.

**Negative controls** (rule (e)).
- **(c1) A planted dropped facet, asserted by construction.** At every width, `count_drop == N_gmsh − 1`, and |`area_drop`/0.48 − 1| > 1e-9.
  - **Ceiling:** one facet's share ≈ 1/N_gmsh ≈ 7e-3 (predicted from the 0.08 size; printed), ≈ 10⁶ × the band. Every facet has positive area, so the miss is arithmetically certain.
  - It shows the 1e-9 band detects the failure mode `PORT-1` feared: one lost facet.
- **(c2) Local-copy double count, printed, predicted, never asserted.** `area_local_201 − 0.48` is predicted **0 at `-n 1`** and ≥ 0 at `-n 2` / `-n 4`. It is > 0 only where a partition boundary lies *on* the rectangle (both adjacent tets on different ranks), which the partitioner does not guarantee. It is printed to show what an unfiltered (rank-local) sum would have reported.

**Window** (foreground, Bash `timeout: 660000`; the harness call is the whole Bash command, bare and repo-relative; `-s`; no pipes; `;` between widths so a hang at one width does not skip the next; the pytest step is last, so the footer is pytest's status).

`scripts/testing/run_and_log.sh GEO-37-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && rm -rf /tmp/geo37 && mkdir -p /tmp/geo37 && PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 1 python3 tests/mesh/probe_geo37_interior_facet.py --out /tmp/geo37/n1.json ; PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 2 python3 tests/mesh/probe_geo37_interior_facet.py --out /tmp/geo37/n2.json ; PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 4 python3 tests/mesh/probe_geo37_interior_facet.py --out /tmp/geo37/n4.json ; PYTHONPATH=/workspace/src:/workspace GEO37_RECORD_DIR=/tmp/geo37 timeout -k 10 30 python3 -m pytest tests/mesh/test_geo37_interior_facet.py -v -s --tb=short'"`

- The `rm -rf /tmp/geo37` first is load-bearing: a stale `n2.json` from an earlier window would fake (A1).
- The command contains no `$`, so the host shell expands nothing. Do not add `echo $?`; (A1) reads completion from the record, not from the exit code.

**Traps.**
- **`dS` on an unghosted mesh** (above, `PORT-1` 3b-iv, two windows).
- **A rank-0 gmsh throw hangs the other ranks** inside the collective import unless it is broadcast first (`GEO-23`, Status 124 at 120 s × 3).
- **A non-cutting fragment.** If `occ.fragment` leaves the face unembedded, its triangles are not tet facets: (A2) reads 0 and the guard raises. Fall back to `gmsh.model.mesh.embed(2, [rect], 3, vol)` before meshing, disclosed in the journal. This is a fixture fix, not a band move.
- **Owned versus local.** Facets are rank-local, and `facet_tags` carries non-owned copies of partition-boundary facets. Filter by `index_map(2).size_local` and reduce. Never sum `facet_tags.values` rank-locally.
- **After a timeout kill, check orphans and `/dev/shm`.** Run `pgrep -c python3` in the container first. A killed `mpiexec` leaks `/dev/shm/mpich_shm_*` (known-issues), and the next `-n 4` launch can then fail for an unrelated reason.
- **Run hygiene:** `PYTHONPATH` carries `:/workspace`; print on rank 0; rule (i).

**Scope.**
- **Not claimed:** a port solve, a CAD/STEP import, curved or multi-face sheets, a face touching a material interface or the boundary, other partitioners or ghost modes, an `xl` width, or any change to the importer (the row: no `src/` change).
- This is one planar rectangle, one mesh size and one partitioner.

**Status it can move.** The `GEO-37` §7 row's Done-when: one logged window at `-n 1 / 2 / 4`; the tagged count equal across widths after reduction; the tagged area at 1e-9 against the closed form; no hang ⇒ **✅ as a probe, green or red**. Both verdicts move ⬜ → ✅:
- **Green** (A1–A3 at every width) registers native facet tags as the (9b) CAD-port route.
- **Red** at any width registers the sliver-pair importer.

The §10 U6 line and the (9b) design sentence are re-pointed by the review accordingly.

**Negative result.**
- **This is a probe whose red is a result.** A red (A1) (a width killed at 45 s, meaning a hang), (A2) or (A3) is recorded as the verdict RED with the table. Commit on main, because the row closes green or red.
  - The slot **stops** after committing, since there is no `PREDICTED RED` line and a red stops the slot.
  - It journals which width and which assertion failed, plus the stacks, if a width hung.
  - The known-issues 9 history gains a 0.11 line. `_interface_facet_tags`' docstring is **not** edited (no `src/` change).
- **(A4) or (c1) red:** the instrument is wrong, and the verdict is void. Fix the probe and re-run the window in-slot, and disclose it. If it is still red, park on `attempt/GEO-37-<UTC>` and mark BLOCKED (rule (d)).
- **The fixture guard raises** (not one volume, or no embedded face): this is a CAD fault, not a verdict. Apply the `embed` fallback once, then BLOCKED.

There is no PREDICTED RED line. The drafter predicts **green**, because the 2026-08 hang was localised to the `dS` assembly, not to `model_to_mesh`'s distribution of interior facet tags. This is predicted and not on record for 0.11.

**Outcome — 2026-10-04 22:30 slot (implementer): GREEN, ✅ as a probe.**
- **W1** `20261005T033147Z_GEO-37-W1.log` (Status 1, 6 s): the fixture guard raised at every width, `expected 1 interior face, got []`. Cause: OCC bounding boxes are padded by the shape tolerance (≈ 1e-7), so the brief's literal 1e-12 z-span classifier can never match. The `GEO-23` broadcast worked — all ranks raised, no hang. Fixture fix (disclosed, not a band): classify at z-span < 1e-6 and mid-plane within 1e-6; the exact `getMass` 1e-12 area guard is unchanged. No `embed` fallback was needed — `occ.fragment` embedded the face.
- **W2** `20261005T033214Z_GEO-37-W2.log` (Status 0, 13 passed, 5 s; smoke by measurement): 10 150 tets, `cells_ghost` 0, N_gmsh = 192 at h = 0.08.

| width | count_201 | area_201 − 0.48 | per-rank owned 201 | local-copy excess (c2) | count_101 / area_101 − 6 | (c1) area_drop/0.48 − 1 | s mesh / import |
|---|---|---|---|---|---|---|---|
| 1 | 192 | +1.1e-16 | [192] | 0 facets / +1.1e-16 | 2422 / 0 | −4.48e-3 | 0.43 / 0.15 |
| 2 | 192 | +5.6e-17 | [91, 101] | 0 / +5.6e-17 | 2422 / 0 | −4.48e-3 | 0.37 / 0.14 |
| 4 | 192 | 0 | [41, 57, 15, 79] | +9 facets / +1.94e-2 | 2422 / 0 | −5.07e-3 | 0.30 / 0.15 |

- (A1)–(A4) and (c1) green at every width. (c2) as predicted: 0 at `-n 1`, ≥ 0 elsewhere — at `-n 4` the partition runs along the rectangle and an unfiltered rank-local sum would over-report by 9 facets. The owned filter is load-bearing.
- **Verdict:** native facet tagging for CAD-drawn ports on 0.11 (9b). The `dS` assembly over such a tag on an unghosted mesh, the 2026-08 hang, was not exercised and remains unverified.

### §9 item 267 — 2026-10-06 (chain U11)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §10 U11 (grep "U11. `GEO-37` step 2"). **No band, tolerance, Done-when or validation target is set or moved**; every band is imported by name.

**`GEO-37` step 2: `dS` over the native interior tag.** Item 222's fixture (unit box, rectangle 0.8 × 0.6 at z = 0.47, tag 201, h = 0.08, 10 150 tets, N_gmsh 192), read through `_model_to_mesh` (`src/fem_em_solver/io/mesh.py:13–28`) at `-n 1 / 2 / 4`, now **assembled** with UFL `dS(201)` — the operation that hung in 2026-08 and that item 222 deliberately avoided (its Outcome: "not exercised and remains unverified").
- **Executor: implementer, not `mesh-probe`** — the same reason as item 222 (rule (j) there): the step asserts an area against a closed form, and `mesh-probe` never asserts.
- **Tier:** smoke by measurement on item 222's fixture (W2 5 s for three widths + pytest, `20261005T033214Z_GEO-37-W2.log`); ceiling standard. Real build (all integrands real). `-n 1 / 2 / 4`.
- **Class / cost:** chain step, no solve, ≈ 17 measured / **≈ 35 predicted** min: 15 fixed + ≈ 12 code + one window ≈ 1 min + ≈ 7 spare.
- **Depends on nothing.** No `src/` change.

**Code (tests only).**
- **`tests/mesh/probe_geo37_interior_dS.py`** (not collected), `--out <json>` and `--ghost {none,shared_facet}`. Reuse item 222's CAD/guard/import verbatim by import from `tests/mesh/probe_geo37_interior_facet.py` (factor its build into a function if needed — a test-file edit, disclosed). `--ghost shared_facet` passes `partitioner=dolfinx.mesh.create_cell_partitioner(dolfinx.mesh.GhostMode.shared_facet, 2)` through the shim (precedent `io/mesh.py:2299–2314`, `io/step_import.py:412`). Rank 0 writes JSON with a `completed: true` stamp last; `cells_ghost` (reduced) recorded.
- **Integrals** (each `fem.assemble_scalar(fem.form(...))` then `comm.allreduce(..., op=MPI.SUM)`; `dS = ufl.Measure("dS", domain=mesh, subdomain_data=facet_tags, metadata={"quadrature_degree": 2})`):
  - (I1) `1 * dS(201)` → area;
  - (I2) `g('+') * dS(201)` with `g = 2 + 3x + 5y` from `ufl.SpatialCoordinate` — continuous, so one-sided equals either side; closed form 0.48 · (2 + 3·0.5 + 5·0.5) = **2.88** (linear integrand on planar triangles: exact to round-off at degree 2);
  - (I3) DG0 `u` = 1 below z = 0.47, 3 above (by cell midpoint): `(u('+') + u('-')) * dS(201)` → 4 · 0.48 = **1.92** (orientation-free); `u('+') * dS(201)` **printed only** (the '+' side is per-facet cell order, not geometry — it is in [0.48, 1.44] and is not a closed form);
  - (I4) the exterior instrument `1 * ds(101)` → 6.
- **`tests/mesh/test_geo37_interior_dS.py`**, skipped unless `GEO37_DS_RECORD_DIR` is set (item 222's pattern); one test per anchor; prints the table and `[GEO-37 step2] VERDICT dS over native tag: GREEN|RED` per ghost mode.

**Anchors (asserted, every width, ghost mode `shared_facet`).**
- **(D1)** JSON `completed: true`, `comm_size == w` (no hang).
- **(D2)** |I1/0.48 − 1| ≤ `STRIP_BAND` (1e-9, `tests/validation/test_ops61_exact_strip_ports.py:58`), imported.
- **(D3)** |I2/2.88 − 1| ≤ `STRIP_BAND`; |I3/1.92 − 1| ≤ `STRIP_BAND`.
- **(D4)** instrument: |I4/6 − 1| ≤ `STRIP_BAND`.
- **(D5) Missing-tag control raises:** `validate_required_port_tags_exist` (`src/fem_em_solver/ports/definitions.py:82`) called on the reduced tag set with tag 999 raises; the same call with 201 passes.

**Negative controls.**
- **(c1) Missing tag, asserted (by construction):** `1 * dS(999)` assembles to **0.0** on every width (printed, then asserted == 0.0) — DolfinX does not raise on an absent tag, which is why (D5) needs the validator; this is the trap the control exists to show.
- **(c2) Planted dropped facet, asserted:** re-tag 201 with its first owned facet on the lowest owning rank removed (a new `meshtags`); |I1′/0.48 − 1| > `STRIP_BAND`; ceiling ≈ 1/192 ≈ 5e-3 (item 222 read −4.5e-3 / −5.1e-3), ≈ 10⁶ × the band.
- **(c3) Ghost mode `none` (the 2026-08 configuration), predicted, printed only.** Item 222's `-n 4` partition runs *along* the rectangle (9 facets / 1.94e-2 local-copy excess): on an unghosted mesh those facets have one local cell. Predicted: `-n 1` identical to `shared_facet`; at `-n 4` either a short area (≈ −2e-2 relative to −4e-2) or a raise/hang inside the 45 s per-width timeout. Never asserted; it is run **after** the asserted legs so a hang cannot cost them.

**Window** (foreground, Bash `timeout: 660000`; whole bare command; `;` between widths; `rm -rf` first; pytest last so the footer is pytest's; no `$`):
`scripts/testing/run_and_log.sh GEO-37-s2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && rm -rf /tmp/geo37ds && mkdir -p /tmp/geo37ds && PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 1 python3 tests/mesh/probe_geo37_interior_dS.py --ghost shared_facet --out /tmp/geo37ds/sf_n1.json ; PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 2 python3 tests/mesh/probe_geo37_interior_dS.py --ghost shared_facet --out /tmp/geo37ds/sf_n2.json ; PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 4 python3 tests/mesh/probe_geo37_interior_dS.py --ghost shared_facet --out /tmp/geo37ds/sf_n4.json ; PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 1 python3 tests/mesh/probe_geo37_interior_dS.py --ghost none --out /tmp/geo37ds/none_n1.json ; PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 2 python3 tests/mesh/probe_geo37_interior_dS.py --ghost none --out /tmp/geo37ds/none_n2.json ; PYTHONPATH=/workspace/src:/workspace timeout -k 10 45 mpiexec -n 4 python3 tests/mesh/probe_geo37_interior_dS.py --ghost none --out /tmp/geo37ds/none_n4.json ; PYTHONPATH=/workspace/src:/workspace GEO37_DS_RECORD_DIR=/tmp/geo37ds timeout -k 10 30 python3 -m pytest tests/environment tests/mesh/test_geo37_interior_dS.py -v -s --tb=short'"`
Predicted 30–60 s (six launches ≈ 3–8 s each incl. a cold FFCx compile of three small facet forms, + pytest ≈ 5 s); worst case (c3) hangs at -n 4 ⇒ +45 s.

**Traps already paid for.**
- **`SpatialCoordinate` facet integral on a gmsh mesh without `metadata={"quadrature_degree": …}`** sent FFCx into a > 9 min compile and poisoned the cache (`POST-5` step 1, 2026-08-18, two windows): pin degree 2 on every measure.
- **Interior facets are `dS` with `('+')` / `('-')` restrictions**; an unrestricted coefficient in `dS` fails to compile. '+' is not a geometric side (I3 printed leg).
- **`assemble_scalar` is rank-local**: allreduce every integral; never sum `facet_tags.values` rank-locally; owned filter by `index_map(2).size_local` for (c2).
- **`dS` on an unghosted mesh hung in 2026-08** (`PORT-1` 3b-iv, `20260805T034007Z_PORT-1-step3biv-hang-localise-fine.log`): asserted legs are ghosted; the unghosted leg is printed and runs last.
- **Rank-0 gmsh throw** must be broadcast before the collective (`GEO-23`); **OCC bbox padding** ≈ 1e-7 (item 222 W1: classify at 1e-6).
- **After any kill:** `pgrep -c python3` in the container first, then `/dev/shm/mpich_shm_*` (known-issues); orphan-`.c` sweep only at `pgrep` 0 (rubric (4) command).
- `-s`; no pipes; `PYTHONPATH` carries `:/workspace`; rule (i). Write ≈, never `~`.

**Scope.**
- **Covered:** one planar rectangle, h = 0.08, default partitioner with `shared_facet` ghosts (asserted) and none (printed), real build, a facet area and two facet-constant/linear integrands.
- **Not claimed:** a port solve or port operator (that is 9b's next step); complex-build forms; curved or multi-face sheets; a sheet on a material interface; STEP import (U12); other partitioners; `xl` widths.

**Status it can move.** `GEO-37` is ✅ (its Done-when is step 1's and is not re-read). On green: the §7 row gains a **step-2 record** sentence (I1–I3 residuals per width, ghost-mode reading, elapsed); §10 U11 marked DONE and U12 unblocked ("9b's port form is the native tag"). On red: U11's red branch — the (9b) port route re-points to the sliver-pair importer `_interface_facet_tags` (§10 U11 text), and known-issues 9 gains a 0.11 `dS` line. Either way the result moves a status, so the item is queueable.

**Negative result.** Report and stop; commit the record on main (a red is U11's result, as item 222's would have been), journal in attempts.md.
- **(D2)/(D3) red at `-n 2`/`-n 4` with (D4) green, ghosted:** native-tag `dS` is wrong on 0.11 → U11 red branch. **(D1) red (a 45 s kill) ghosted:** the 2026-08 hang reproduces on 0.11 → red branch, stacks journaled.
- **(D4), (c1) or (c2) red:** the instrument is wrong, verdict void; fix in the new test files once, re-run in-slot, disclose; still red ⇒ park `attempt/GEO-37-step2-<UTC>`, mark item 267 BLOCKED (rule (d)).
- **(c3) any reading** is a record, never a verdict.

No `PREDICTED RED:` line. The drafter predicts **green** on the ghosted legs.

**Outcome — 2026-10-06 13:30 slot (implementer, second item): GREEN.**
- **W1** `20261006T185004Z_GEO-37-s2-W1.log` (Status 124, 87 s): all six probe legs ran as in W2 (same numbers), but the pytest leg's 30 s ceiling was spent inside `tests/environment/test_dolfinx_api_migration.py`'s repository scan before the reader started. Instrument fix (disclosed, no band touched): **W2 narrows the environment path to `tests/environment/test_complex_mode.py`** (the scalar-type check, PASSED: `ScalarType=float64`, real build). Probe and reader unchanged between W1 and W2.
- **W2** `20261006T185157Z_GEO-37-s2-W2.log` (Status 0, 56 s, 23 passed, 3 skipped; smoke by measurement). The `STRIP_BAND` import is by name from `test_ops61_exact_strip_ports.py`; item 222's `_build_model` is imported verbatim (no edit to item 222's files).

| ghost / width | cells_ghost | I1/0.48 − 1 | I2/2.88 − 1 | I3/1.92 − 1 | u⁺ only | I4/6 − 1 | dS(999) | (c2) drop | owned / one-sided 201 per rank |
|---|---|---|---|---|---|---|---|---|---|
| sf / 1 | 0 | −5.6e-16 | +1.6e-15 | −5.6e-16 | 1.3843 | −4.6e-15 | 0.0 | −4.48e-3 | [192] / [0] |
| sf / 2 | 592 | +8.9e-16 | +2.2e-16 | +8.9e-16 | 0.9747 | +2.0e-15 | 0.0 | −4.48e-3 | [91, 101] / [6, 5] |
| sf / 4 | 1162 | +2.2e-16 | +2.2e-16 | +2.2e-16 | 1.0695 | +1.1e-15 | 0.0 | −5.07e-3 | [41, 57, 15, 79] / [4, 7, 11, 8] |
| none / 1 | 0 | as sf / 1 | | | | | | | |
| none / 2 | 0 | as sf / 2 (bit-identical) | | | | | | | [91, 101] / [0, 0] |
| none / 4 | — | **killed at 45 s, no integral printed (both windows)** | | | | | | | |

- (D1)–(D5), (c1), (c2) green at every width. (D5): the reduced tag set is {101, 201}; the validator passes on (201, 101) and raises `missing required port tags: [999]`. (c1): `dS(999)` assembles to exactly 0.0 — DolfinX is silent on an absent tag, so the validator is load-bearing.
- (c3) record: ghost `none` at `-n 2` never puts a 201 facet on a partition boundary (one-sided count 0), so it does not exercise the hazard; at `-n 4` (item 222: partition runs along the sheet) the run is killed at the 45 s ceiling before its first print — the 2026-08 hang reproduces on 0.11 (not localised: import, `compute_midpoints` or first assembly; no stacks). Assembly at the ghosted widths takes 0.13–0.34 s. `pgrep -c python3` 0 after each kill; the leaked `mpich_shm_*` segments were cleared.
- **Verdict:** §10 U11 green — 9b's port form is the native tag, **on a `shared_facet`-ghosted mesh**; U12 unblocked. Not claimed: anything in the brief's Scope list; why the unghosted leg hangs.

## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 **Does an interior 2D face reach DolfinX's facet tags on the 0.11 image? The CAD-port probe — Phase 9 (9b)** (opened 2026-10-03 by the weekly; §10 chain U6; implementer — amended from `mesh-probe` by the 2026-10-04 03:00 review, rule (j): the Done-when asserts an area, which `mesh-probe` never does) *(history and §9 briefs: `docs/planning/chunks/GEO-37.md`)* | ✅ **GREEN as a probe 2026-10-04; audit PASS 2026-10-05 03:00 review** (§9 item 222, `20261005T033214Z_GEO-37-W2.log`, Status 0, 13 passed, 5 s, smoke by measurement): on the 0.11 image, through the project shim `_model_to_mesh` with the default partitioner, the interior dim-2 group 201 arrives in `facet_tags` complete and unduplicated — owned count 192 = gmsh's 192 triangles at `-n 1 / 2 / 4`, owned area − 0.48 = +1.1e-16 / +5.6e-17 / 0 (band 1e-9), exterior instrument 2422 facets, area 6 exact; no hang (import 0.14–0.15 s). Planted dropped facet detected at every width (−4.5e-3 / −5.1e-3). Unfiltered rank-local sum over-reports at `-n 4` by 9 facets / +1.94e-2 (printed: the owned filter is load-bearing). ⇒ native facet tagging is the (9b) CAD-port route; one planar rectangle, one mesh size, one partitioner, no `dS` assembly (the 2026-08 hang) exercised. **Step-2 record (§9 item 267, 2026-10-06, `20261006T185157Z_GEO-37-s2-W2.log`, Status 0, 56 s, 23 passed, smoke by measurement):** UFL `dS(201)` assembled on the same fixture with `GhostMode.shared_facet` (cells_ghost 0 / 592 / 1162) at `-n 1 / 2 / 4`: ∫1 − 0.48 relative −5.6e-16 / +8.9e-16 / +2.2e-16, ∫g('+') (g = 2 + 3x + 5y, 2.88) +1.6e-15 / +2.2e-16 / +2.2e-16, ∫(u⁺ + u⁻) (1.92) −5.6e-16 / +8.9e-16 / +2.2e-16, exterior ∫1 ds(101) (6) ≤ 4.6e-15 — all inside `STRIP_BAND` 1e-9; absent tag `dS(999)` assembles to 0.0 silently while `validate_required_port_tags_exist` raises on it; one dropped facet reads −4.5e-3 / −4.5e-3 / −5.1e-3. Ghost mode `none` (printed, never asserted): `-n 1` / `-n 2` identical to the ghosted values (the `-n 2` partition leaves no one-sided 201 facet), **`-n 4` killed at the 45 s per-width ceiling in both windows before its first integral printed** — the 2026-08 unghosted hang reproduces on 0.11 at the width whose partition runs along the sheet (not localised). ⇒ the (9b) port form is the native tag on a ghosted mesh (§10 U11 DONE). **Done-when (§4):** a box with one embedded interior rectangle tagged as a dim-2 physical group, read at `-n 1 / 2 / 4`: the tagged facet count equal across widths after reduction and the tagged area against the rectangle's closed form at 1e-9, no hang, in one logged window ⇒ ✅ as a probe, green or red. Green ⇒ native facet tagging for CAD-drawn ports; red ⇒ the importer inserts a sliver pair and reuses `_interface_facet_tags`. Measurement only, no `src/` change. 
