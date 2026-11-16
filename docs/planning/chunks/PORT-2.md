# PORT-2 — Port data model and tagging contract

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 247 — 2026-10-05 (lane 1)

*Drafted by a general-purpose agent for the 2026-10-05 03:00 daily review; the review reviews it.*

**`PORT-2` step 1: the port tagging contract as a round trip on the exact-strip fixture.** This is the lane-1 target the 10-03 weekly named (PROJECT_PLAN.md:4645–4646): "the port tagging contract as a round-trip test on the exact-strip fixture ⇒ ✅".
- **Executor:** implementer.
- **Tier:** standard by measurement of the same build. The row's "smoke" is the pure-Python unit module's word and is re-labelled (see the row amendment). No solve. Complex env, because the fixture's records are tagged to the complex 0.11 image. Two windows: W1 at `-n 2`, W2 at `-n 1`.
- **Class:** lane 1. The median is ≈ 8 min. **Cost:** ≈ 8 measured / ≈ 30 predicted min (15 fixed, ≈ 10 code, two windows of ≈ 1.5 min each).
- **Depends on nothing.** No `src/` change: one new test module.

**The contract, read off the code.**
- `src/fem_em_solver/ports/definitions.py` contains:
  - the tagging contract comment, `:14–22` (facets primary, edges fallback, name template `port_<id>_<suffix>`);
  - `PortDefinition`, `:28–61`;
  - `required_port_tags`, `:72–79`;
  - `validate_required_port_tags_exist`, `:82–111`. Its warning (`:88–97`) requires a **globally reduced** tag set: `OPS-14` measured a rank-local read raising on 4/4 ranks at `-n 4`.
  - `validate_port_ordering`, `:114–135`, and `run_port_calibration_checks`, `:224–246`.
- On the birdcage the definitions are built in `_four_port_rung` (`tests/validation/test_port_birdcage_leg_offset_sweep.py:394–405`):
  - `positive_tag = SHEET_IFACE + i`. This is a **facet** tag, rebuilt on the dolfinx side by `_interface_facet_tags` (`src/fem_em_solver/io/mesh.py:290–373`) and narrowed by `_narrowed_radial` (`:166–196`).
  - `negative_tag = CONDUCTOR_CELL_TAG`. This is a **cell** tag.
  - `orientation = "leg_gap_axial_plus_z"`, and `z0_ohm = REFERENCE_IMPEDANCE_OHM`.
  - The matching `LumpedSheetPortSpec` objects are built at `:380–391`.
- **Nothing on the lumped-sheet route checks that a `PortDefinition`'s `positive_tag` equals its spec's `facet_tag`.** `run_n_port_sparameter_sweep` checks only validity and uniqueness of IDs (`src/fem_em_solver/ports/sparameters.py:394–402`).
- `tests/ports/test_port_definition.py` exercises the model on literal integers only. No test has read a mesh's tags back through it. The round trip is that test.

**Code (tests only): new `tests/validation/test_port2_exact_strip_tag_roundtrip.py`.**
- **Fixture (module scope).**
  - Call `_build(np.zeros(4), leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` (`tests/mesh/test_birdcage_leg_offset.py:107`).
  - Then run `_sheet_areas` and, per port, `_sheet_azimuth_deg` / `_port_frame` / `_projected_extents` / `_narrowed_radial`. This is the sequence of `_rung` (`tests/validation/test_ops61_exact_strip_ports.py:79–117`) and of `_four_port_rung` (`:297–333`), imported, never retyped.
  - Build the `PortDefinition` objects and the specs exactly as `:380–405` does. The specs are constructed only so that their `facet_tag` can be read; there is no solve.
- **Reduced tag sets.**
  - `facet_set` = the sorted union of `comm.allgather(np.unique(tags.values).tolist())`, the idiom at `src/fem_em_solver/ports/excitation.py:271–275`.
  - `cell_set` is built the same way from `cell_tags.values`.
- **Helper `_roundtrip(port_defs, specs, msh, tags, comm)`.** It runs every collective on every rank and returns per-port records. It then raises `ValueError("port tag/geometry mismatch: …")` naming every port for which either:
  - the spec's `facet_tag ≠ positive_tag`, or
  - the tagged facet set's quarter-turn index `round(az/90) mod 4 ≠ i − 1`.

**Anchors (asserted, both widths).**
- **(T0) Premise.**
  - `0.5·dx·dz` equals `SPEC_AREA_M2` (`tests/validation/test_ans4_exact_strip_records.py:131–133`, the 0.007 × 0.008 CAD rectangle) within `STRIP_BAND` (`test_ops61_exact_strip_ports.py:58`, 1e-9).
  - Cells vs `ANS4_EXACT_STRIP_CELL_RECORD` (`:64`) within `CELL_COUNT_BAND`.
- **(T1) Existence, and the right dimension.**
  - `validate_required_port_tags_exist(port_defs, facet_set ∪ cell_set)` returns.
  - Stricter, and asserted separately because the validator does not know about dimensions: every `positive_tag ∈ facet_set` and `negative_tag ∈ cell_set`.
- **(T2) Non-empty, with an exact owned-count reduction.** Per port:
  - `n_i` = the allreduced count of facets `< index_map(2).size_local` carrying the tag (`_owned_facet_count`, `test_ops61_exact_strip_ports.py:72–76`).
  - Assert `n_i > 0`.
  - Assert `n_i == |⋃_ranks local_to_global(all local facets carrying the tag)|`. This is an integer identity: a reduction that counts ghost copies, or loses one, fails it.
  - Print `n_i` beside the 40 triangles recorded at `-n 1` and `-n 2` (`20261001T215942Z_OPS-61-F10a-W1.log:3576–3579`; `20261001T220051Z_OPS-61-F10a-W2.log:3591–3594`). This is printed only: a mesher record, never registered as a constant.
- **(T3) Area.** `|A_i / SPEC_AREA_M2 − 1| ≤ STRIP_BAND` per port, with `A_i` from `_interface_area_or_zero`.
  - *Asserted*: +1.110e-15 … +1.554e-15 on this fixture at `-n 2` (`…220051Z…:3591–3594`).
- **(T4) Bijection.**
  - `port_id → positive_tag` is injective, and its image equals `{SHEET_IFACE + i}` ⊂ `facet_set`.
  - Each spec's `(port_id, facet_tag)` equals its definition's `(port_id, positive_tag)`.
  - The quarter-turn index of P_i's tagged set is `i − 1`. This is an integer test with no band. P1 reads ≈ 360.0000° (`20261001T220154Z_OPS-61-F10a-W0-regate.log:5313`), so the `mod 4` is load-bearing.
  - `validate_port_ordering(..., expected_port_order=["P1", "P2", "P3", "P4"])` returns.
  - `run_port_calibration_checks` passes on `port_face_areas = {P_i: A_i}` at its own default ratio of 1.25. Print the ratio.

**Negative controls (asserted, by construction; every rank).**
- **(C1) Missing tag.**
  - Drop every facet of `SHEET_IFACE + 3` from the tags on every rank, and rebuild `meshtags` with sorted indices (`_narrowed_radial`'s argsort idiom).
  - Re-reduce the tag set. `validate_required_port_tags_exist` must raise `ValueError` matching `r"missing required port tags: \[<SHEET_IFACE+3>\]"`.
  - The raise is deterministic and identical on every rank, because the set is reduced.
- **(C2) Swapped tag.** Relabel the facets of `SHEET_IFACE + 2` ↔ `SHEET_IFACE + 4`.
  - Assert that `validate_required_port_tags_exist` **passes**. The existence contract is blind to a swap; this is a set identity.
  - Assert that (T3) also passes. The areas are equal, as the backing above shows.
  - Assert that `_roundtrip` raises `"port tag/geometry mismatch"` with both `P2` and `P4` in the message. Their quarter-turn indices are off by 2, so the failure is arithmetically certain.
  - This is the control that shows why (T4) exists.
- **(c3) Printed, predicted, never asserted.** The rank-local read `OPS-14` paid for: whether each rank's `np.unique(tags.values)` alone holds all four sheet tags, gathered and printed.
  - At `-n 1` it does trivially. At `-n 2` it depends on the partition.
  - Compute it with no validator call inside a rank-divergent branch.

**Windows** (foreground, Bash `timeout: 660000`; each harness call is the whole Bash command, bare and repo-relative; `-s`; no pipes; no `-k` selection). The existing unit module runs in the same window.

`scripts/testing/run_and_log.sh PORT-2-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/ports/test_port_definition.py tests/validation/test_port2_exact_strip_tag_roundtrip.py -s -v --tb=short'"`

`scripts/testing/run_and_log.sh PORT-2-step1-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment tests/ports/test_port_definition.py tests/validation/test_port2_exact_strip_tag_roundtrip.py -s -v --tb=short'"`

**Cost, from logs.**
- One exact ×1 rung builds in 26.62 s at `-n 2` (`20261001T220051Z_OPS-61-F10a-W2.log:3590`) and 27.05 s at `-n 1` (`…215942Z…W1.log:3575`). The two-rung `OPS-61` module took 57 s at either width.
- `tests/environment` plus two pure ports modules took 51 s at `-n 1` (`20261002T201726Z_PORT-3-W2.log`).
- Predicted ≈ 70–100 s per window: standard, inside the 180 s ceiling.
- **Box-speed gate:** a rung print over 60 s means a slow box, which is a price. A kill at 180 s stops the item BLOCKED with the reading; never lengthen the timeout.

**Traps.**
- **Rank-local data:** `cell_tags.values` and `facet_tags.values` are rank-local; this is the very `OPS-14` defect the validator's docstring names. Filter facets by `index_map(2).size_local` before counting.
- **Reductions:** `_projected_extents` / `_sheet_azimuth_deg` reduce through `allgather`, elementwise (`PORT-9` 3c: `allreduce` on arrays falls through to Python's `min`).
- **Collectives:** every rank runs every collective, and asserts come after them. `pytest.raises` runs on every rank. A `comm.` call under `if comm.rank == 0:` deadlocks (item 227).
- **`meshtags`** needs sorted, unique indices.
- **Azimuth wrap:** 360° ≡ 0°.
- **No UFL form beyond `_interface_area_or_zero`'s `dS`.** It ran green at `-n 2` on this fixture (`OPS-61` W2).
- **Compiles:** sweep orphan `.c` stems only with `pgrep -c python3` reading 0.
- **Run hygiene:** rule (i).

**Scope.**
- One fixture (exact ×1, zero offsets), one mesh, no solve.
- **Not exercised:**
  - physical-group names (`PORT_TAG_NAME_TEMPLATE`, `definitions.py:22`). The birdcage carries integer tags only; named CAD ports are (9b)'s and `GEO-37`'s;
  - the edge fallback;
  - ragged fixtures;
  - a displaced-leg layout.
- The contract is silent on dimension. On this fixture the negative terminal is a cell region; that is recorded, not ruled on.
- A geometry-aware swap check in `src/` would be a feature, not this lane's. The test-local helper is the deliverable.

**Negative result.**
- **(T1) or (T4) red on the unplanted fixture:** the port-ID ↔ sheet contract that every four-port record relies on is broken.
  - Stop. Open a known-issues entry naming the port(s), and park on `attempt/PORT-2-step1-<UTC>`, BLOCKED in the same commit (rule (d)).
- **(T2) red:** this is a rank-safety defect in the tag reduction; known-issues entry, same disposition.
- **(T3) or (T0) red:** the exact-strip build regressed. Re-run `tests/validation/test_ops61_exact_strip_ports.py` in-slot to tell which; stop.
- **(C1) / (C2) red (a plant that does not fire):** the instrument is wrong. Fix it once, disclose, and re-run. If it is still red, BLOCKED.
- There is no PREDICTED RED line; the drafter predicts green.

**Status it can move.** `PORT-2` 🧪 → ✅.
- **Row amendment (rule (j)), in the review's own commit.** The §7 row (PROJECT_PLAN.md:2667) has no Done-when. The review writes:
  > **Done-when (2026-10-05 03:00 review, rule (j); the 10-03 weekly's lane-1 target, §9 item 247):** on the exact-strip ×1 birdcage (`leg_sheet_strip_fraction = GATED_WIDTH_FRACTION`), the `PortDefinition`s built as `_four_port_rung` builds them round-trip at `-n 2` and `-n 1`:
  > - `validate_required_port_tags_exist` passes on the globally reduced tag set, with positive tags among the facet tags and the negative tag among the cell tags;
  > - every port's facet set is non-empty, with an exact owned-count reduction;
  > - every area equals the 0.007 × 0.008 CAD rectangle within `STRIP_BAND`;
  > - port IDs ↔ tags are bijective and place P_i on leg i;
  > - a planted missing tag raises "missing required port tags", and a planted swap passes the existence check but fails the round trip, naming both ports;
  > - `tests/ports/test_port_definition.py` is green in the same window ⇒ 🧪 → ✅.
- The Tier cell becomes "smoke (unit module); standard (round trip, measured)".

**§9 line 2:** → docs/planning/chunks/PORT-2.md#§9-item-247

**Outcome — ✅ DONE 2026-10-05 09:00 CDT slot (implementer; `PORT-2` 🧪 → ✅).** New `tests/validation/test_port2_exact_strip_tag_roundtrip.py` (complex_only), no `src/` change; every helper and band imported (`_build`, `_sheet_areas`, `_sheet_azimuth_deg`, `_port_frame`, `_projected_extents`, `_narrowed_radial`, `_owned_facet_count`, `STRIP_BAND`, `CELL_COUNT_BAND`, `ANS4_EXACT_STRIP_CELL_RECORD`, `SPEC_AREA_M2`). W1 `-n 2` 27 passed, Status 0, 86.9 s (`20261005T140943Z_PORT-2-step1-W1.log`); W2 `-n 1` 27 passed, Status 0, 90.8 s (`20261005T141117Z_PORT-2-step1-W2.log`). First run of each, no edits between them. Rung 26.07 s / 27.30 s (box-speed gate 60 s: normal box).

**Step 1 narrative (2026-10-05).**
- (T0) cells 117 206 = record at both widths; 0.5·dx·dz = `SPEC_AREA_M2` within `STRIP_BAND`.
- (T1) reduced facet set [211, 212, 213, 214], cell set [1, 2, 3, 101…104, 201…204]; `validate_required_port_tags_exist` returns; every positive tag is a facet tag and the negative tag 1 (`CONDUCTOR_CELL_TAG`) a cell tag.
- (T2) n_i = 40 owned facets on every port and 40 = |⋃ local_to_global| at both widths (the mesher's 40-triangle record reproduced, printed only).
- (T3) areas 5.600000000000e-05, relative +1.110e-15 / +1.554e-15 / +1.110e-15 / +1.332e-15 at `-n 2`, +1.110e-15 / +1.110e-15 / +8.882e-16 / +1.110e-15 at `-n 1` (band 1e-9).
- (T4) positives injective with image {SHEET_IFACE + i}; spec `facet_tag` = `positive_tag` for every port; azimuths 360.0000 / 90 / 180 / 270°, quarter-turn indices [0, 1, 2, 3] (the `mod 4` carried P1); ordering and `run_port_calibration_checks` pass, max face-area ratio 1.000000000000000.
- (C1) dropping 213 raises `missing required port tags: [213]` on every rank.
- (C2) swapping 212 ↔ 214 passes the existence check and the area band; `_roundtrip` raises naming P2 (quarter 3 vs 1, az 270°) and P4 (quarter 1 vs 3, az 90°) and neither P1 nor P3.
- (c3), printed: at `-n 2` the per-rank local facet tag sets each lack at least one sheet tag ([False, False]); at `-n 1` [True]. A rank-local `tags.values` handed to the validator would have raised on both ranks of this partition — the `OPS-14` defect, here on facets.
- Not exercised (scope as briefed): named physical groups, the edge fallback, ragged or displaced fixtures. The contract stays silent on dimension; on this fixture the negative terminal is a cell region (recorded, not ruled on).
