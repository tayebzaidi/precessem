# GEO-7 — history (also carries GEO-2 and GEO-6)

Item briefs and execution history (operator directive 2026-10-01). The §7 rows are the status of record. One utility serves `GEO-2`, `GEO-6` and `GEO-7` (the 2026-09-30 weekly's lane-1 target); their briefs live here.

### §9 item 184 — 2026-10-02 (lane 1)

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.*

**`GEO-2` / `GEO-6` / `GEO-7`: one CAD-identity QA utility.** This is the 09-30 weekly's target.
- **Executor:** implementer. **Tier:** standard, real build, `-n 2` then `-n 1`.
- **Cost:** ≈ 4 measured / ≈ 35 predicted min.
- **Depends on:** nothing. **Attempt:** 1 of 4.

**Why these rows were rejected before** (09-29 18:00 review §6): `GEO-7` was string formatting, `GEO-6` had no identity on graded cells (its report compares cell-count ratios), and `GEO-2`'s title and code disagree. This item gives all three one volume-identity anchor.

**Code (additive, rule (c) disclosure): `src/fem_em_solver/io/mesh_qa.py`.**
- New function:

  `cad_identity_report(mesh, cell_tags, *, expected_tags, box_volume=None, exact_volumes=None, record_volumes=None, rtol=1e-9, comm=None)`
- It returns a frozen `CadIdentityReport` with these fields:
  - `tag_volumes`, `total_volume`;
  - `missing_tags`, `unexpected_tags`;
  - `partition_rel` = |Σ_expected V_t / V_mesh − 1|;
  - `box_rel`, `exact_rel{tag}`, `record_rel{tag}`;
  - `flags`, a subset of `tag_set` / `partition` / `box` / `exact_volume` / `record`;
  - `passed`.
- Tag sets come from **owned** cells only (`indices < size_local`), allgathered.
- Volumes use `tests/mesh/helpers.py:74–88`'s expressions **verbatim**: `subdomain_id=(int(tag),)`, `Constant(1.0)`, `assemble_scalar`, `allreduce` SUM, `np.real`.
- Export it from `io/__init__.py`. The default `rtol` is `VOLUME_PARTITION_BAND` (`helpers.py:71`).

**Test: new `tests/mesh/test_mesh_qa_cad_identity.py`.** It imports the existing modules' constants and calls and restates none of them.
- **(A0)** The `rtol` default `==` `VOLUME_PARTITION_BAND`.
- **(A1)** Every `tag_volumes[t]` and `total_volume` `==` `helpers.tag_volume` / `helpers.total_volume`, bitwise, on all three meshes.
- **(A2) Records reproduced, with `passed` True on each:**
  - **R1, `coil_phantom_domain`** (`test_mesh_tag_integrity.py:35–43`, called as at `:101`):
    - tag set {1, 2, 3, 4};
    - `record_volumes = UNIFORM_VOLUMES_RECORD` (`:89–94`) at 1e-9, the gate's own band (`:212–218`), reproduced on 0.11 at `20260930T123217Z_GEO-5-W1.log:605–609`;
    - `box_volume = test_coil_phantom_conforming._analytic_box_volume()` (`:61–80`), the identity at `:129`;
    - the partition identity (`:136`).
  - **R2, `two_torus_domain`** (`test_two_torus_conforming.py:36–49`, call copied from `:75–84`):
    - tag set {1, 2, 3} (`:86`);
    - box from `_analytic_box_volume(…)` (`:52–55`, identity at `:97`);
    - partition (`:104`).
  - **R3, `birdcage_port_domain`** in box-port mode at `BASELINE_CONTROL_RESOLUTION` 4.8e-3 (`test_birdcage_conductor_sizing.py:101`, call `:108–123`; constants `test_birdcage_port_tags.py:22–32`):
    - tag set {1, 2, 3, 101…104};
    - box from `test_birdcage_port_tags._analytic_box_volume()` (`:35–52`);
    - partition;
    - `exact_volumes` {101…104: 0.010·0.008·0.010 = 8.0e-7}, the identities at `:168–184`;
    - cells printed against 33 185 (`20260826T093202Z_GEO-21-step2-gate.log:11075`), not asserted.

**Negative controls (planted on the built meshes, no remesh).** Each verdict pattern is **asserted**, since it holds by construction of the plant. Each numeric miss is **predicted** and printed beside the reading (rule (e)).
- **(N1) Tag merged.** Relabel every R1 tag-3 value to 4 (owned and ghost).
  - Expect flags == {`tag_set`} (missing {3}), with partition and box still ≤ 1e-9.
  - Pass no records.
- **(N2) Ungrouped cell.** Drop one owned tag-3 cell from R1's MeshTags, on the lowest rank that owns one.
  - Expect flags == {`partition`}.
  - The miss is that cell's volume / V_mesh, predicted ≥ 1e-7 (≥ 100× the band) and printed.
- **(N3) Port box re-tagged** (`GEO-2`). Re-tag one owned tag-101 cell to 2 on R3.
  - Expect flags == {`exact_volume`}, for 101 only.
  - The miss is that cell's volume / 8e-7, predicted ≥ 1e-3 and printed.

**Windows.** Foreground, Bash `timeout: 660000`, each harness call bare and whole.
- **W1** (≈ 90 s):

  `scripts/testing/run_and_log.sh GEO-7-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_mesh_qa_cad_identity.py tests/io/test_mesh_qa_diagnostics.py tests/mesh/test_mesh_tag_integrity.py tests/mesh/test_geometry_sanity_report.py tests/mesh/test_birdcage_port_tags.py::test_birdcage_port_layout_diagnostics_match_the_closed_forms tests/mesh/test_birdcage_port_tags.py::test_birdcage_port_layout_rejects_too_small_or_overlapping_port_regions -v -s --tb=short'"`
- **W2** (≈ 40 s): the same with `-n 1`, labelled `GEO-7-W2`, and only `tests/environment tests/mesh/test_mesh_qa_cad_identity.py`.

**Traps:**
- **Planted tags:** build them with `dolfinx.mesh.meshtags` on a **copy** of `values`, with sorted, unique, int32 indices.
- **Rank-local tags:** `cell_tags.values` is rank-local, and the birdcage mesh is ghosted (`shared_facet`, `mesh.py:4052`). Read owned cells only.
- **Test selection:** never select with `-k`; node IDs are fine. Import constants and helpers only, never `test_*` names.
- **Do not run the ungraded birdcage partition test.** It is red; that is item 186's business. Items 184 and 186 both read `test_birdcage_port_tags.py`, but only 186 edits it. They are independent.
- **Run hygiene:** use `-s`. Make no edit after the last window (rule (i)).

**Scope:**
- `coil_phantom_geometry_sanity_report` (cell ratios) stays as it is: a print, not a gate.
- `helpers.py` is not rewired onto the utility.
- Port-face facet areas and sheet ports (`GEO-18`, `TH-16` 3a) are out of scope.
- No band moves.

**Status it can move:** `GEO-2`, `GEO-6` and `GEO-7` each go 🧪 → ✅ on the Done-whens this review wrote into their rows (rule (j)).

**Negative result:**
- **(A1) not bitwise:** the expression differs. Fix it once, then park.
- **An (A2) record misses 1e-9:** the generator drifted. Print every reading, open a known-issues entry, and stop. Never re-record in-slot.
- **A plant raises a second flag, or misses its own:** stop.
- **On any stop:** park on `attempt/GEO-7-<UTC>` and mark BLOCKED in the same commit (rule (d)).
- There is no PREDICTED RED line.

**Outcome — DONE 2026-10-03 (19:30 CDT slot, implementer).**

Windows:
- W1 (`-n 2`): `20261003T005340Z_GEO-7-W1.log`, 91 s, 25 passed, 4 skipped. The skips are the complex-only environment tests.
- W2 (`-n 1`): `20261003T005520Z_GEO-7-W2.log`, 73 s, 15 passed, 4 skipped.

Anchors:
- **(A0)** Green.
- **(A1)** Bitwise on R1, R2 and R3.
- **R1, `coil_phantom_domain`:** `record_rel` is 1.55e-10 / 4.1e-11 / 8.9e-11 / 1.68e-10 for tags 1–4. Partition and box are ≤ 4.2e-15.
- **R2, `two_torus_domain`:** partition and box are ≤ 9.1e-15.
- **R3, the birdcage:** the four port boxes read 8.000000000e-07, with `exact_rel` ≤ 4.4e-16. Partition and box are ≤ 1.6e-14. Cells are 33 185, equal to the record (printed only).

Controls:
- **(N1)** {`tag_set`}, missing {3}.
- **(N2)** {`partition`}. The miss is 1.88e-5 / 2.33e-5 against a predicted ≥ 1e-7. The cell differs per width.
- **(N3)** {`exact_volume`}, tag 101 only. The miss is 4.17e-2 against a predicted ≥ 1e-3.

Disclosed (rule (c)): the change to `io/mesh_qa.py` and `io/__init__.py` is additive. `src/` cannot import `tests/`, so the default band is a module constant, `CAD_IDENTITY_RTOL = 1e-9`. (A0) asserts it equals `helpers.VOLUME_PARTITION_BAND`.

`GEO-2`, `GEO-6` and `GEO-7` all go 🧪 → ✅.

**Audit, 2026-10-03 03:00 daily review (§4; one `auditor` per row).** Item 184 shares two windows across three rows.
- **`GEO-2`: PASS.** The anchors are `report.flags == frozenset()` (test `:211–212`) and the N3 flag pattern (`:289–291`). Port boxes read 8.000000000e-07 (`20261003T005340Z_GEO-7-W1.log:1998`), `exact_rel` 4.441e-16 (W1:2001), and the window took 91 s (W1:3518–3520). Caveats: the Done-when cites `test_birdcage_port_tags.py:86` / `:261`, but the tests sit at `:96` / `:304`, a line-number slip in an unchanged file. The layout tests ran at `-n 2` only, which is what the Done-when asks for.
- **`GEO-6`: PASS on substance. The auditor's DEMOTE was on the tier word alone, and this review overrules it.** The row read "smoke" against windows of 91 s (W1) and 73 s (`…005520Z_GEO-7-W2.log:2071–2073`), so the tier is relabelled **standard** with no re-run. That follows the 09-30 precedent: a tier-word DEMOTE is overruled and the cell is corrected, as for `MAT-15`. Caveat: the row's "≤ 1.6e-14 on both generators" is R3's figure. `GEO-6`'s two generators read 4.2e-15 / 9.1e-15 (W2:323, W2:481).
- **`GEO-7`: PASS.** All three generators pass with `flags=[]`: R1 at W1:343, R2 at W1:502 and R3 at W1:1997. The merged-tag control flags exactly `tag_set` (W1:2007–2009). Caveat: the row does not state its elapsed time. The 91 s is in the footer and the `test-results.md` row.
- **All three rows:** the log headers name the parent commit `7e33984d`, because the code was uncommitted at run time. The file mtimes (19:52:53 / 19:53:30 CDT) precede W1, which started at 19:53:41 CDT.
- **Review spot-check:** W1:1998 re-read, and it matches.
