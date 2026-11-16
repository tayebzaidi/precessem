# GEO-35 — Drop the floating leg-port boxes from unbroken-leg birdcages (history)

### §9 item 275 — 2026-10-06 (chain U17)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: the 2026-10-05 weekly's enumerated step **U17** (PROJECT_PLAN.md §10, grep "U17. `GEO-35` step 1"; target "CAD identities"; serial on nothing). **No band, tolerance, record or Done-when is set or moved.** Every band below is imported by name from the file:line cited.

**`GEO-35` step 1: the inert leg-port boxes dropped behind a knob.** `birdcage_port_domain` gains one additive keyword-only argument, `emit_leg_ports: bool = True`. With `False` on an unbroken-leg layout (`leg_gap_length=None`) the 16 floating boxes at the mid-leg azimuths (`src/fem_em_solver/io/mesh.py:4817–4846`, the `addBox` at `:4837`) are not built; everything else is unchanged. On the F-human fixture (`WF-7`'s build, `tests/validation/test_wf7_f_human_ring_matrix.py:451–456`) the step asserts the CAD identities and the ring-port census, prints the cell count, and prints how far the degree-1 32×32 moves from `WF-7`'s record. The identities are the gate. The move is expected to be at discretisation level and is printed only.
- **Executor:** implementer. Changes: additive `src/` in `io/mesh.py`, one new mesh test module, and one additive environment knob in the `WF-7` step-1 module. **No existing default path, record, band or other caller is edited, and no fixture is switched.** Switching the fixtures and re-recording under (1*) is step 2 (`record-reconciler`).
- **Tier:** W0 standard; W1 and W2 **heavy**. Ranks: `-n 2` for W0 and W1 (the rank-safety width); `-n 8` for W2 (`WF-7`'s own width, under the 12 ceiling).
- **Class:** chain step (solve). **Cost:** ≈ 19 measured (class median) / **≈ 48 predicted** min. That is 15 fixed + ≈ 15 code + W0 ≈ 2 + W1 ≈ 5.5 + W2 ≈ 6.5 + ≈ 4 spare. **If W1 is not green by minute 35, park** (rule (d)).
- **Depends on:** nothing open. It shares no file with any queued item that this draft knows of. `io/mesh.py` is also touched by the `GEO-36` / `GEO-34` lines, so add by edit and rebase, never check out by path.

**Settled from the code (read before designing).**
- **Where the boxes come from.** `port_tags` is filled once per leg in `for idx, angle in enumerate(theta_leg)` (`mesh.py:4830`). With `leg_gap_length is None` each box sits at the mid-azimuth between adjacent legs (`:4831–4836`) and is an axis-aligned `addBox` (`:4837–4846`). The ring-gap boxes are appended after it. The constructor already refuses to drive these boxes: `emit_port_sheets` without a gap raises an error (`:4346–4353`), and `leg_azimuth_offsets_rad` without a leg gap is refused because the boxes "float" (`:4379–4390`).
- **Ordinals are positional.** `port_ordinal = {tag: idx for idx, tag in enumerate(port_tags, start=1)}` (`:5218`). Every ring-port caller keys on ordinals `LEG_COUNT + 1 … 3·LEG_COUNT` = 17–48 (`_ring_ports`, `tests/validation/test_birdcage_f_human_rung.py:134–136`), that is, cell tags `PORT_LOWER + i` / `PORT_UPPER + i` = 117–148 / 217–248. **With the knob off the leg boxes, the ring ordinals must stay 17–48.** Build the map explicitly and keep the leg ordinals 1–16 as reserved gaps. A plain `enumerate` over a shorter list would renumber every ring port.
- **The group census.** `expected_groups` (`:5527–5530`) demands `100 + i` for every `port_ordinal` value, and `group_names[100 + idx]` (`:5479`) names them. Both must skip the dropped ordinals, or the build raises a missing-group error.
- **The F-human fixture.** It is built from `_params(F_HUMAN_RING_RADIUS, scale_sizing=False)` (`test_birdcage_f_human_rung.py:101–131`: `LEG_COUNT` 16 at `:76`, radius 0.15 at `:79`, `leg_gap_length=None`, ring gaps, `emit_port_sheets=True`), plus `ring_sheet_orientation="longitudinal"`, which is what `WF-7` adds (`test_wf7_f_human_ring_matrix.py:451–452`). `_build_f_human` (`test_wf7_f_human_b1_quadrature.py:64`) drives only the ring ordinals.
- **The CAD identities and their bands.** These are the four from `GEO-25` step 2 (`test_birdcage_f_human_rung.py:232–288`): cells within `CELL_COUNT_BAND` (1 %, `tests/mesh/test_birdcage_port_sheet_prerequisite.py`) of `F_HUMAN_BRANCH_B_CELL_RECORD` 504 642 (`:87`); the `GEO-18` volume partition at `EXACT` (`tests/mesh/test_birdcage_ring_gaps.py`); `GEO-19` terminal-area ratios on all 32 ring ports in `TERMINAL_RATIO_BAND` (0.95, 1.0) (`:98`); and meshed/CAD conductor mass ≥ `CAD_MASS_GATE` (`tests/mesh/test_birdcage_conductor_sizing.py`). `_build_rung` (`:139–193`) computes the readings. Its `port_cell_tags` includes the 16 leg tags `PORT_LOWER + i` (`:163`), so the new module must build its own tag list from the ring ports only. Do not reuse the leg entries: an absent tag has zero volume, so the partition would still pass and hide the drop.
- **`WF-7`'s degree-1 record.** It is the 0.11 image at `-n 8` (`20261005T003402Z_WF-7-H15a-r3-W2.log`): cells **507 266** (`:10583`, relative 5.200e-03 to 504 642; build 129.65 s at `:10521`); reciprocity 1.041217e-14 (`:10584`); σ_max **0.999914648** (`:10585`); 18 class means with worst spread 3.6997e-02 against band 0.05 (`:10587–10605`); `S_driven(P17)` **0.407423+0.344417j** (`:10606`); P17 `P_src,exact` 7.037117061e-03 W and pooled κ 4.590991e-03 (`:10524`). The window took **357.9 s** with 11.149 GiB summed `ru_maxrss` (`:10610`) and finished "3 passed, 1 skipped" (`:10613`). The module constants are `SIGMA_MAX_D1_RECORD`, `CLASS_MEAN_D1_RECORD`, `P17_P_SRC_EXACT_RECORD_W`, `P17_POOLED_KAPPA_RECORD` and `S_DRIVEN_STEP0_RECORD` (`test_wf7_f_human_ring_matrix.py:154–194`).

**Code** (the names are proposals; the executor may rename and disclose the mapping).
- (1) **`MeshGenerator.birdcage_port_domain`, additive.** Add the keyword-only `emit_leg_ports: bool = True` after `implant_rod`. The request contradictions are checked before any gmsh state exists, in the existing pattern. `False` raises `ValueError` with `leg_gap_length` set, because the leg boxes are the drive in the gapped layout. For step 1 it also raises with `as_hole` or `symmetry_quadrant`; those are later steps (the `WF-12` / `TH-20` layouts). With `False`, the leg iteration skips `addBox`, the ordinal map is explicit (ring ordinals unchanged), and `expected_groups` / `group_names` skip ordinals 1…`leg_count`. The diagnostics gain `"leg_ports_emitted": bool`. **`git diff` on `mesh.py` shows insertions plus `if emit_leg_ports` guards around the existing lines only.** Disclose anything else.
- (2) **New `tests/mesh/test_geo35_leg_port_boxes.py`.** It imports `_params`, `_ring_ports`, `F_HUMAN_RING_RADIUS`, `LEG_COUNT` and the bands above, plus `_tag_volume` / `_total_volume` / `_port_boundary_partition` exactly as `_build_rung` does. Nothing is restated. Nodes:
  - `test_geo35_knob_contract` (no mesh): the default is `True` (`inspect.signature`); the three refusals raise before gmsh.
  - `test_geo35_f_human_drops_leg_ports_keeps_identities` (W1): one default build and one `emit_leg_ports=False` build in one process, at `-n 2`, with the `WF-7` kwargs.
- (3) **`tests/validation/test_wf7_f_human_ring_matrix.py`, additive knob** `FEM_EM_GEO35_DROP_LEG_PORTS=1`, read inside the test like the `POST-7` knobs (`:393–412`). When it is set, `kwargs["emit_leg_ports"] = False` and four anchors become **printed, not asserted**, with a relative move against each record: (i) the cell record; the P17 reproduction; the save-file (L2) reproduction (refused if combined with `FEM_EM_WF7_SAVE_NPZ`); and the `S_driven(P17)` / σ_max / 18 class-mean figures. (ii) reciprocity, (iii) passivity, (iv) the C16 × mirror classes and (v) the exact identity **stay asserted at their imported bands**. When the knob is unset, no key is added to `kwargs`, so the landed run's path is unchanged. This is the `phantom_resolution` precedent (`test_wf7_f_human_b1_quadrature.py:64–80`, "adds no key … byte-identical").

**Anchors (asserted).**
- **(K) Contract, no mesh.** The default is `True`; `False` with `leg_gap_length`, `as_hole` or `symmetry_quadrant` raises `ValueError` before gmsh. Exact.
- **(D) The knob-off build is digit-identical (W1, default build).** Its global cell count `== 507 266` exactly, the `WF-7` record at `…r3-W2.log:10583`, same image and same kwargs. gmsh meshes serially, so the count does not depend on the rank count. The leg-port tags 101–116 are present and the cell count is also inside `CELL_COUNT_BAND`. **How "digit-identical" is asserted:** (a) this exact count; (b) the `mesh.py` diff touches the default path only through `if emit_leg_ports` guards; (c) the `WF-7` knob adds no key when unset; (d) W0 re-runs `tests/mesh/test_birdcage_port_tags.py`, which builds the floating-box layout with the default, green at its unchanged bands. A 32×32 re-run with the knob off is **not** bought, because (c) makes it the landed run.
- **(N) No leg-port tags on the knob-on build.** The global cell census of tags `PORT_LOWER + i` and `PORT_UPPER + i` for i = 1…16 is exactly 0 (allreduced), and `diag["leg_ports_emitted"] is False`.
- **(R) Ring-port census unchanged.** On the knob-on build, every ring ordinal 17–48 has both halves `PORT_LOWER + i` / `PORT_UPPER + i` with global cell count > 0. The 32 sheet interfaces built by `_interface_facet_tags` (as `test_wf7_f_human_b1_quadrature.py:118–121` does) each have area > 0. The sorted list of ring cell tags equals the default build's list.
- **(C) CAD identities on the knob-on build, imported bands.** Volume partition (conductor, phantom, air and ring halves only) within `EXACT` of 1. All 32 terminal ratios in `TERMINAL_RATIO_BAND`. Meshed/CAD conductor ≥ `CAD_MASS_GATE`.
- **(I) The identities gate the knob-on 32×32 (W2, degree 1).** These are `WF-7`'s (ii)–(v) at their imported bands: reciprocity ≤ `RECIPROCITY_BAND`; σ_max ≤ `MATRIX_PASSIVITY_CEILING`; exactly `N_SYMMETRY_CLASSES` (18) classes, worst spread ≤ `OPPOSITE_SPREAD_BAND`; the exact identity ≤ `DISCRETE_IDENTITY_RTOL` on the driven columns. The column-scale matrix control stays asserted as landed.

**Negative controls (rule (e)).**
- **(c1) The leg-tag detector, on the default build.** *Asserted:* the same census that reads 0 in (N) must read > 0 on all 16 leg tags of the knob-off build in the same process. This shows that (N)'s zero is not an empty comparator. It is a separation of "16 present" against "0 present"; no factor is involved.
- **(c2) The partition comparator.** *Asserted:* on the knob-off build, dropping the 16 leg tags from the tag list must break the partition by more than `EXACT`, by 16 box volumes / V_total. Predicted ≈ 16 × (0.0214 × 0.0171 × 0.0214) m³ over V_total. The boxes are `PORT_BOX_SIZE` scaled by 0.15 / 0.07. Print the measured reading. This shows that (C)'s partition sees missing volume.
- **(c3) `WF-7`'s identity control** (omit the smallest undriven sheet, 4.8e+02× the band at `…r3-W2.log:10524`) stays **printed** on the landed route, as landed. Predicted ≈ 480× on the knob-on build.

**Printed (not asserted).**
- The knob-on cell count against 507 266 and against 504 642. Predicted to fall by a few thousand cells, roughly 0.3–3 %. **A reading outside `CELL_COUNT_BAND` is not a red**: it is the record that step 2 re-records.
- Build seconds for both builds.
- The degree-1 32×32 move against `WF-7`'s record: relative |Δ| of `S_driven(P17)`, σ_max, P17 `P_src,exact`, pooled κ, and each of the 18 class means beside `CLASS_MEAN_D1_RECORD`. Predicted ≤ 1 % on every figure (the boxes are inert vacuum at the mid-leg azimuths). That is inside the 5.50 % degree-1 order caveat (`WF-7` step 0b). A move larger than 5.50 % is a finding for the review, not a red.

**Windows.** Run every window in the foreground; the harness call is the **whole** Bash command, bare and repo-relative. Bash `timeout: 660000`; `-s`; node IDs, never `-k`; `tests/environment` first; no pipe or redirect (every window is under 10 min).
- **W0, regression on the default path (rule (c): `io/mesh.py` is edited).** Standard, `-n 2`. Unpriced in this draft; predicted ≤ 120 s. **Mesh-only.**
  `scripts/testing/run_and_log.sh GEO-35-s1-W0-porttags "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_port_tags.py tests/mesh/test_geo35_leg_port_boxes.py::test_geo35_knob_contract -v -s --tb=short'"`
- **W1, mesh-only on the F-human fixture (the gate for (D), (N), (R), (C), (c1) and (c2)).** Heavy, `-n 2`. Predicted ≈ 300 s, bracket 260–420 s: two builds at ≈ 112–130 s each (`GEO-25` step 1: 112.27 s at `-n 1`, `20260904T214031Z_GEO-25.log:9733`; `WF-7`: 129.65 s at `-n 8`), plus environment ≈ 45 s. **Gate: W1 > 480 s ⇒ stop BLOCKED before W2.**
  `scripts/testing/run_and_log.sh GEO-35-s1-W1-mesh "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_geo35_leg_port_boxes.py::test_geo35_f_human_drops_leg_ports_keeps_identities -v -s --tb=short'"`
- **W2, the knob-on degree-1 32×32 (gate (I) and the printed move).** Heavy, `-n 8`. Predicted ≈ 360 s, bracket 330–450 s (`…r3-W2.log:10610`: 357.9 s, 11.149 GiB). `split` keeps the ≈ 135 s orientation-control build out of the window, as r3 did.
  `scripts/testing/run_and_log.sh GEO-35-s1-W2-wf7 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_ORIENTATION_CONTROL=split FEM_EM_GEO35_DROP_LEG_PORTS=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short'"`

**Traps already paid for.**
- **Positional ordinals** (`mesh.py:5218`). If the ring ports are renumbered to 1–32, every downstream sheet, spec and record silently re-keys. (R) is the check, and it is meant to fail on exactly this.
- **An absent tag has zero volume.** A partition computed over a tag list that still names the 16 leg tags passes vacuously. (C) uses the ring-only list, and (c2) shows that the comparator sees missing volume.
- **Rank-safety.** `cell_tags.values` is rank-local. Every census is allreduced, using `_tag_volume` / `_port_boundary_partition` with `comm` (`GEO-9` step 2b paid for this). At `-n 2` rank 0 does not own every port. `np.max` / `np.min` on an owned slice raise on an empty rank, so use `initial=` and reduce.
- **A gmsh build that fails deadlocks at rank > 1** (`GEO-23`). The (K) refusals must raise **before** any gmsh state exists, on every rank.
- **The cell record is not a red here.** The knob-on count is printed. The `WF-7` knob must demote (i) and the degree-1 reproduction anchors to printed-only, or W2 fails on records that step 2 owns. Never widen `CELL_COUNT_BAND` or `REUSE_REPRODUCTION_RTOL`.
- **`FEM_EM_WF7_SAVE_NPZ` + the knob** would overwrite or compare against the degree-1 record file. The knob refuses the combination.
- **Do not edit any exact rectangles or sheets.** The ring sheets are untouched, and κ stays a diagnostic (operator ruling 2026-10-02).
- **Orphan `.c`.** W2 compiles the degree-1 forms, which should be warm from r3. If a window is killed: `pgrep -c python3` first, then the rubric's sweep. Never use `pgrep -f "python3 -m pytest"`.
- **Status 124 is a finding, not a retry.** Never lengthen a timeout or change `-n`.
- **Rule (i):** no edit to the new module, the knob or `mesh.py` after their last green window.
- **Write ≈, never `~`, in tracked Markdown.** No AED file is touched.

**Scope.**
- **Covered:** the F-human fixture only (branch B, 0.15 m, 16 legs, ring gaps, longitudinal sheets) with the knob on; mesh identities at `-n 2`; one degree-1 32×32 at 64 MHz at `-n 8`.
- **Not claimed:** switching any fixture's default; `mesh:9`–`12`, `WF-10`, `EX-62`, `WF-12` (as_hole) or `TH-20` (quadrant); any re-recorded constant; degree 2; B1+, SAR or tuning; any absolute S.

**Status it can move.** No glyph. Step 1 meets the row's first Done-when clause: a mesh test asserts no leg-port tags and an unchanged ring-port census on the switched fixture. It does **not** meet the second clause: records re-recorded under (1*) by `record-reconciler`, with C16 / mirror / exact identities re-run green on the switched fixtures. **The Done-when does not match step 1.** `GEO-35` goes ⬜ → 🟡 at most, with a step-1 record sentence: (N) / (R) / (C) green, the knob-on cell count, and the 32×32 move. §10 U17 is marked DONE.

**Negative result.** In every case: report the measurement, park on `attempt/GEO-35-step1-*`, mark item 275 BLOCKED in the same commit (rule (d)), and journal it in attempts.md.
- **(K) or W0 red:** the default path changed. Revert (1); BLOCKED.
- **(D) count ≠ 507 266 with W0 green:** an image or environment drift since 2026-10-05 (the 10-05 `e0672b91` docker change is the first suspect). This is a finding for `record-reconciler` under (1*). Do not re-register it in-slot; stop before W2.
- **(N) or (R) red:** a renumbering or census bug in the new code. Fix and re-run in-slot only if the fix stays inside the `if not emit_leg_ports` branches; otherwise stop.
- **(C) red on the knob-on build only:** removing the boxes changed the conductor mesh. That is a real finding, so stop; no band moves.
- **(c1) or (c2) not separating:** the comparator is broken. Stop.
- **(I) red in W2:** an inert-box removal broke an identity. That is a real finding; the mesh legs stand as a record.
- **W1 over its gate, or W2 killed at 590 s:** BLOCKED without a result; the review re-poses the step.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Outcome (2026-10-06 19:30 CDT slot, second item): GREEN.** W0 `20261007T004148Z_GEO-35-s1-W0-porttags.log` (Status 0, 76 s, 17 passed); W1 `20261007T004311Z_GEO-35-s1-W1-mesh.log` (Status 0, 366 s < 480 s gate, 13 passed); W2 `20261007T005016Z_GEO-35-s1-W2-wf7.log` (Status 0, 461 s, body 385.9 s, 10.756 GiB, 16 passed / 1 skipped). `mesh.py` +40 / −0: the knob, refusals before gmsh, `leg_ports_emitted`, and three `if not emit_leg_ports` guards — **disclosed:** besides `port_ordinal`, the ring-sheet `ordinal = len(port_tags)` (`ring_sheet_of_ordinal`) was also positional and needed its own guard. The `WF-7` module's existing label strings gained a knob branch; the default-path asserts are unchanged under `if not geo35`. (K) green, no gmsh state. (D) 507 266 exactly (5.1997e-03 to 504 642). (N) all 32 leg tags 0 cells. (R) ring tags 117–148 / 217–248 equal to the default build's, min 331 cells, min sheet area 1.716154e-04 m². (C) knob-on partition 1.000000000000, terminal ratios 0.973927–0.974323, mass 0.965320, identical to the default build. (c1) leg tags > 0 (max 180). (c2) deficit 1.303288e-03 = prediction (1.303e+06 × `EXACT`). (c3) 480×, printed. (I) reciprocity 1.165976e-14, σ_max 0.999914557, 18 classes worst spread 3.6503e-02 (band 0.05), exact identity 4.190e-15. Knob-on cells 502 704 (−8.9933e-03 vs 507 266, −3.8403e-03 vs 504 642). Printed moves: `S_driven(P17)` 2.635e-04, σ_max 9.078e-08, P17 `P_src,exact` 7.706e-05, class means 3.3e-06 … 9.100e-04; **pooled κ 4.805416e-03 vs 4.590991e-03, 4.671e-02 — the ≤ 1 % prediction is missed** (inside the 5.50 % caveat, a finding for the review, not a red; κ spreads 4.46e-03 … 5.13e-03 across the 32 drives, so it reads local mesh, not field). `GEO-35` ⬜ → 🟡; U17 DONE.

## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 **Drop the floating leg-port boxes from unbroken-leg birdcages** — with `leg_gap_length=None` `birdcage_port_domain` still emits one port box per leg, "floating in the air at the *midpoint* azimuth between adjacent legs, touching no metal" (its docstring; the `PORT-9` step 3 placeholder ports), so the F-human (high-pass, ring-gap) fixture carries 16 inert vacuum boxes with their own `PORT_LOWER/UPPER + i` cell tags and interior sheets in the isocenter plane — no port or material is ever defined on them (`_build_f_human` drives only the ring ordinals 17–48), but the mesh conforms to their faces and they clutter every ParaView view (operator, 2026-09-30, interactive session: "they shouldn't be there"). **Shape:** an additive keyword on `birdcage_port_domain` (e.g. `emit_leg_ports`, default unchanged so every existing caller builds the identical mesh) that omits the leg-port boxes when the legs are unbroken; the high-pass fixtures (F-human, `mesh:9`–`12`, `WF-7`/`WF-10`, `EX-62`) switch to it. **Done-when:** a mesh test asserts no leg-port tags and an unchanged ring-port census on the switched fixture; every record the switch moves (F-human cell counts such as `F_HUMAN_BRANCH_B_CELL_RECORD` 504 642, printed `WF-7` / `EX-62` figures) is re-recorded under the (1*) licence by `record-reconciler`, with the gated anchors (C16, mirror, exact identities) re-run green at their unchanged bands — no band moves. | 🟡 ~~Deferred to Phase 9 clean-up (operator, 2026-09-30): not queueable before Phase 6 closes~~ — placed by the 2026-10-05 weekly as §10 chain U17 (step 1, §9 item 275); step 1 covers the mesh-test clause only and can move ⬜ → 🟡, the fixture switch and re-record are step 2. Step 1 (2026-10-06 19:30 slot): `emit_leg_ports=False` on F-human leaves zero leg-port cells, ring ordinals 17–48 intact and the `GEO-25` CAD identities green (partition 1.000000000000, terminal ratios 0.9739–0.9743, mass 0.9653), 502 704 cells (−0.90 % vs 507 266); the knob-on degree-1 32×32 holds reciprocity / passivity (σ_max 0.999914557) / 18 classes / exact identity, every printed figure moving ≤ 9.1e-4 except pooled κ (4.67 %, inside the 5.50 % degree-1 caveat; `20261007T004311Z_GEO-35-s1-W1-mesh.log`, `20261007T005016Z_GEO-35-s1-W2-wf7.log`). *History: `docs/planning/chunks/GEO-35.md`.* 

### §9 item 351 — 2026-10-09 (chain step U21)

**Outcome (2026-10-10, 19:30 slot): GREEN, `GEO-35` holds 🟡.** W0 17 passed 73 s (`20261010T005604Z_GEO-35-s2a-W0-contracts.log`); W1 hole 369 387 rel 0 at `:6220`, 192 s; W2 17 passed + 1 skipped, 456 s, cells 502 704, P17 rel 4.326e-11 / κ 1.136e-09, (L2) σ_max 2.250e-10 / class 4.810e-10, (c2) 77.1× / 910.0×, (N′) 0 (`20261010T010041Z_GEO-35-s2a-W2-wf7-step1.log:9710–9820`); W3 17 passed, 315 s, C16 0.1831 %, mirror 0.1264 %, identity 3.476e-16, S[P17,P17] equal to the re-record to 6 digits (`20261010T010824Z_GEO-35-s2a-W3-wf7-step2.log:10707–10715`). Mapping: `GEO35_WF7_CELL_RECORD` → `WF7_CELL_RECORD` / `WF7_BOXED_CELL_RECORD`; knob helper `_boxed_knob`; boxed constants `*_BOXED`. Re-record done in-slot under Q1 (no `record-reconciler` spawn). The (i) cell band is now asserted on both knob routes.

*Drafted by a general-purpose agent for the 2026-10-09 18:00 daily review; the review reviews it.* Licence: the 2026-10-07 weekly's enumerated step **U21** (PROJECT_PLAN.md §10, grep "U21. `GEO-35` step 2"; "serial on `xl` 10-09 (H15b's matrix is saved on the boxed mesh first)"). That prerequisite landed: `20261009T070009Z_WF-7-step0e-H15b.log` (boxed build 507 266 cells at `:10924`, `/workspace/logs/wf7/wf7-s32-degree2-n16.npz` saved and round-tripped at `:12318–12319`). **No band, tolerance or Done-when is set or moved.** Every band below is imported by name. Every new record value is copied from a footered log on the current image and cited by log:line; none is computed.

**This item is step 2a of two (U21 split; question Q4 below).** U21 as written does not fit one slot. It touches four families with moved records: the `WF-7` σ = 800 solid matrix and quadrature (steps 1 and 2); the 10 g SAR family (`WF-7` step 3, `EX-62` / `mri:4`); `GEO-25` / `mesh:12`; and the κ-corrected `WF-10` / `ANS-8` lineage. Re-running all of them is ≈ 50 min of windows alone (prices under *Remains*). **2a lands the switch and re-records the `WF-7` step 1 / step 2 family**: the records `U27` (`WF-16` step 8, "serial on U21") reproduces. Every other `_build_f_human` caller is **pinned to the boxed build in the same commit**, and the copper hole / quadrant paths are **exempt by construction**. That makes 2a independently landable: no pinned or exempt caller's mesh changes, which the W0 contract asserts without meshing and W1 asserts on a real build.

**`GEO-35` step 2a: the solid F-human fixture switches to `emit_leg_ports=False`, and `WF-7` step 1's degree-1 records are re-recorded under (1\*).**
- **Executor:** an implementer slot, then `record-reconciler`. The **implementer slot** writes the code (1)–(4) below and runs W0 and W1. It then **spawns `record-reconciler` in the foreground** (the spawn prompt states the foreground rule: no backgrounded harness window, Bash timeout 660000) with the site table below. The reconciler edits the constants and runs the anchors W2 and W3. One commit at the end holds code, logs and plan text. If Q1 is not ruled, do not spawn.
- **Tier:** W0 standard. W1 heavy by ceiling. W2 and W3 heavy. Ranks: `-n 2` for W0 and W1 (rank-safety width; W1 is xl entry 32's own build width); `-n 8` for W2 and W3 (`WF-7`'s record width, under the 12 ceiling).
- **Cost:** ≈ 13 measured (class median) / **≈ 53 predicted** min. That is 15 fixed + ≈ 12 code + W0 ≈ 1.5 + W1 ≈ 3.1 + W2 ≈ 7.7 + W3 ≈ 4.5 + reconciler spawn and edits ≈ 6 + ≈ 3 spare. **If W2 is not green by minute 42, park** (rule (d)).
- **Depends on:** H15b (landed, above); `GEO-35` step 1 (landed, item 275); **the review's ruling Q1**. It shares `test_wf10_f_human_tuned_circuit.py` with two parked branches (see *Settled* (c)). Pins are two one-line edits there; add them by edit, never by checking out a path.

**Settled from the code (read before designing).**
- **Who calls `_build_f_human`** (`tests/validation/test_wf7_f_human_b1_quadrature.py:64–208`; kwargs at `:114–119`, built from `_params(F_HUMAN_RING_RADIUS, scale_sizing=False)` + `ring_sheet_orientation="longitudinal"`):
  - **Solid route, switched in 2a:** `WF-7` step 2's fixture (`b1_quadrature.py:226`).
  - **Solid route, pinned boxed in 2a; their records move in 2b:**
    - `test_wf7_f_human_sar10g_hotspot.py:154` (`WF-7` step 3; `examples/mri/04_f_human_sar10g_hotspot.py` imports it, so `EX-62` follows the pin);
    - `test_wf10_f_human_tuned_circuit.py:159` and `:454`;
    - `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py:193`.
  - **`as_hole=True`, exempt (boxed):**
    - `test_wf12_f_human_hole.py:425`;
    - `test_wf12_f_human_copper_matrix.py:148`;
    - `test_wf12_f_human_copper_degree2.py:114` (**xl entry 32, queued for 10-13**);
    - `test_wf12_f_human_copper_c_sweep.py:279`;
    - `test_wf10_f_human_copper_esr_root.py:310`;
    - `examples/mri/06_f_human_copper_loss_partition.py:115–167` (a pass-through wrapper, `*args, **kwargs`).
  - **Its own build, switched in 2a:** `WF-7` step 1 builds its fixture itself, not through `_build_f_human` (`test_wf7_f_human_ring_matrix.py:602–606`).
- **Why the hole and quadrant paths must stay boxed:** `birdcage_port_domain` refuses `emit_leg_ports=False` with `as_hole` or `symmetry_quadrant` (`src/fem_em_solver/io/mesh.py:4272–4285`; step 1's (K)). A `_build_f_human` that forwarded the new default on the hole route would make every copper caller raise, including xl entry 32 on 10-13. **(a) Settled: the `WF-12` hole mesh, `TH-20`'s quadrant and every copper path are out of scope for U21.** Lifting the refusal is a `mesh.py` change with its own CAD gate (a later `GEO-35` step, if the weekly wants one). `WF-12`'s own solid control (`test_wf12_f_human_hole.py:286`) calls `birdcage_port_domain` directly and is untouched.
- **What the switch moves on the solid route.** These are the step 1 knob-on readings on the exact kwargs 2a switches to (`20261007T005016Z_GEO-35-s1-W2-wf7.log`, Status 0 at `:10305`):
  - cells 507 266 → **502 704** (`:9694`), −0.384 % to `F_HUMAN_BRANCH_B_CELL_RECORD` 504 642. That is inside `CELL_COUNT_BAND` 1 %, so no 1 % cell assert anywhere goes red;
  - P17 `P_src,exact` −7.706e-05 and pooled κ **+4.671e-02** (`:9733`);
  - σ_max −9.078e-08 (`:9760`);
  - class means up to 9.100e-04 (`:9762–9780`);
  - `S_driven(P17)` 2.635e-04 (`:9757`, `:9760`).
  
  `REUSE_REPRODUCTION_RTOL` is `STEP2_COLUMN_POWER_RTOL` = 1e-6 (`tests/validation/test_port_birdcage_ring_matrix.py:164`, `:171`). So the reproduction asserts at `test_wf7_f_human_ring_matrix.py:1021–1028` (P17) and `:1045–1057` (through the file) **go red on the switched default unless re-recorded**. That is the set 2a re-records.
- **Mesher stable since step 1.** H15b rebuilt the boxed fixture at exactly 507 266 on 10-09 (`…H15b.log:10924`), the same count step 1's (D) asserted on 10-07 (`20261007T004311Z_GEO-35-s1-W1-mesh.log:19872`). The 10-07 knob-on readings are therefore current-image measurements, as (1\*) requires.
- **(b) Saved-npz consumers.** Two files are on the boxed mesh, each carrying `cells`, `order`, `slot_k`, `slot_sign`, `s32` and κ keys (`…H15b.log:12318`):
  - `/workspace/logs/wf7/wf7-s32-degree1-n8.npz` (κ-off degree 1, H15a);
  - `/workspace/logs/wf7/wf7-s32-degree2-n16.npz` (H15b).

  Their consumers:
  - H15b's readout: the `ANS-8` degree-2 rows through `ring_mode_rows_v2`, for the 10-10 weekly (`xl-pending.md:1595–1597`);
  - H13b (`TH-20` step 2, whose §7 row names H15b's saved matrix);
  - the optional `ANS-8` κ-off degree-1 comparand on the degree-1 file (`docs/planning/chunks/WF-7.md:98`).

  All of them read stored arrays and rebuild nothing, so the switch cannot break them **provided** (i) neither file is ever overwritten, (ii) any fresh solve paired with them is built boxed, and (iii) a consumer that pairs a build with a file checks the file's `cells` against that build (507 266 boxed). 2a enforces (i) with a save-path guard. It provides (ii) with `_build_f_human(..., emit_leg_ports=True)` and a module knob `FEM_EM_WF7_LEG_PORTS=boxed`, and documents (iii). *Observation for the review, outside 2a:* H13b's quarter is the copper **hole** layout (`TH-20.md:22`), while H15b's file is the σ = 800 **solid** route. A 1e-9 reassembly of one against the other is not like for like, whatever GEO-35 does. `TH-20.md:568` reads the H15b number as H13b's *price* comparand.
- **(c) Open `attempt/*` branches.** `attempt/WF-10-step1b-20260929T020339Z` (+337 lines) and `attempt/WF-10-step1d-20261001T140319Z` (+712 lines) both append to `test_wf10_f_human_tuned_circuit.py`. 2a's two pins there (`:159`, `:454`) are a small textual conflict if either branch is ever revived. Their disposition stays with the review; nothing in 2a depends on them. `attempt/WF-12-step3*` add only `test_wf12_f_human_copper_c_mode.py` (an `as_hole` module, exempt) and a script. `attempt/PORT-24-step1-*` (`ports/ring_modes.py`), `attempt/OPS-50-*` (a probe) and `attempt/TH-15-step2proper-*` (two-torus) do not overlap.

**Code** (the names are proposals; the executor may rename and disclose the mapping). **No `src/` edit; `mesh.py` is not touched.**
- (1) **`_build_f_human`** gains keyword `emit_leg_ports: bool | None = None`. The kwargs construction is lifted into a pure helper `_f_human_kwargs(phantom_resolution=None, as_hole=False, emit_leg_ports=None) -> dict` (rule (a) lift, no mesh) and called from `_build_f_human`. Resolution:
  - `None` resolves to boxed under `as_hole` and to `False` (dropped) otherwise;
  - the key `emit_leg_ports=False` is added **only** when the resolved value is `False`. The boxed result **adds no key**, so the kwargs are byte-identical to the pre-switch dict (the `phantom_resolution` / `as_hole` precedent, `:116–119`);
  - `as_hole=True` with `emit_leg_ports=False` raises `ValueError` in the helper, before any gmsh state.
- (2) **Pins.** Pass `emit_leg_ports=True` at `test_wf7_f_human_sar10g_hotspot.py:154`, `test_wf10_f_human_tuned_circuit.py:159` and `:454`, and `08_f_human_quadrature_mode_content.py:193`, each with the comment `# GEO-35 step 2b: unpin with the re-record`. The `as_hole` callers are not edited.
- (3) **`test_wf7_f_human_ring_matrix.py`.**
  - **Retire** step 1's `FEM_EM_GEO35_DROP_LEG_PORTS` knob (`:231–252`). The default build (`:602–606`) adds `kwargs["emit_leg_ports"] = False`.
  - Add the inverse knob `FEM_EM_WF7_LEG_PORTS`: unset ⇒ dropped; `boxed` ⇒ no key (the pre-switch build), with the moved reproductions **printed** against the `*_BOXED` constants (the mirror of step 1's knob semantics); any other value raises.
  - The default path asserts again the (i) cell band (`:993–996`), the P17 reproduction (`:1021–1028`) and, with a save path, (L2) (`:1045–1057`). They are **unchanged lines**; only the knob guard around them changes.
  - `_save_knob` (`:255–278`) **refuses the two boxed basenames** `wf7-s32-degree1-n8.npz` and `wf7-s32-degree2-n16.npz` under every setting. Its contract test's "good" path (`:376`, `:420`) moves to `/workspace/logs/wf7/wf7-s32-degree1-n8-noleg.npz` (a literal in a contract test, not a band; disclose it).
  - New assert **(N′)**: on the solved mesh, the allreduced cell census of tags `PORT_LOWER + i` / `PORT_UPPER + i`, i = 1…16, is 0 (import step 1's census helper from `tests/mesh/test_geo35_leg_port_boxes.py`; do not restate it).
  - New no-mesh node `test_wf7_leg_port_knob_contract`.
  - `_default_orientation_extent` (the split orientation control's transverse build) is **left boxed and untouched**. Its assert concerns ring-sheet extent, which does not depend on the leg boxes.
- (4) **`tests/mesh/test_geo35_leg_port_boxes.py`** gains the no-mesh node `test_geo35_f_human_builder_contract`, which checks:
  - `_f_human_kwargs()` equals the step 1 W1 knob-on kwargs (`_wf7_kwargs()` + `emit_leg_ports=False`) as a dict;
  - `_f_human_kwargs(emit_leg_ports=True)` and `_f_human_kwargs(phantom_resolution=0.0075, emit_leg_ports=True)` equal the pre-switch dicts, with **no `emit_leg_ports` key**;
  - **`_f_human_kwargs(as_hole=True)` has no `emit_leg_ports` key and equals `_wf7_kwargs()` + `as_hole=True`** (the xl entry 32 path);
  - `_f_human_kwargs(as_hole=True, emit_leg_ports=False)` raises `ValueError`, with `gmsh.isInitialized()` False before and after.

**Record sites for `record-reconciler`.** All are in `tests/validation/test_wf7_f_human_ring_matrix.py` unless named. Under (1\*), each old value stays in a comment with its provenance log and date and the tag "boxed fixture, pre-`GEO-35`". The old P17 and class-mean values also stay as executable `*_BOXED` constants, because control (c2) and the `boxed` knob read them.

| constant | old (boxed) | new (leg-box-free) | move | new value's provenance |
|---|---|---|---|---|
| `P17_P_SRC_EXACT_RECORD_W` `:162` | 7.037117061e-03 | 7.036574789e-03 | −7.706e-05 | `20261007T005016Z_GEO-35-s1-W2-wf7.log:9733` |
| `P17_POOLED_KAPPA_RECORD` `:163` | 4.590991e-03 | 4.805416e-03 | **+4.671e-02 (Q1)** | `:9733` |
| `SIGMA_MAX_D1_RECORD` `:176` | 0.999914648 | 0.999914557 | −9.078e-08 | `:9736` |
| `CLASS_MEAN_D1_RECORD` ×18 `:177–196` | as listed | the 18 `.9e` values | ≤ 9.100e-04 | `:9762–9779` |
| `S_DRIVEN_STEP0_RECORD` `:156` (printed) | 0.407423+0.344417j | 0.407315+0.344327j | 2.635e-04 | `:9757` |
| `S_DRIVEN_STEP1_PRINTED` (`b1_quadrature.py:57`, printed) | 0.407423+0.344417j | 0.407315+0.344327j | 2.635e-04 | `:9757` |
| `GEO35_WF7_CELL_RECORD` `:228` (printed) → `WF7_CELL_RECORD` 502 704 + `WF7_BOXED_CELL_RECORD` 507 266 | 507 266 | 502 704 | **−0.899 % (Q1)** | `:9694`; `…W1-mesh.log:19883` |

**Not edited:**
- `F_HUMAN_BRANCH_B_CELL_RECORD` 504 642 (GEO-25's *transverse* build; its knob-on count is unmeasured; 2b);
- `scripts/probes/wf7_step0_f_human_cost.py:157` (a historical probe);
- `src/fem_em_solver/study.py:1170` (a measured price row for `Study.estimate()`, history);
- the `TERMINAL_RESIDUAL_RECORD_*` constants (`:208–219`; they parse committed logs);
- `ANS-8`'s `KAPPA_RECORD` (pinned, 2b).

Doc copies (record-reconciler rule 4): one sentence each on the §7 `WF-7` and `GEO-35` rows. Historical briefs and logs stay as written; they are records of the boxed fixture with their citations.

**Anchors (asserted, bands imported).**
- **(K) Contracts, no mesh (W0).**
  - Step 1's `test_geo35_knob_contract` (`test_geo35_leg_port_boxes.py:142`), unchanged.
  - (4)'s builder contract, including the **as_hole exemption**.
  - (3)'s leg-port knob contract.
  - The save-knob contract with the two boxed basenames refused.
- **(H) The exempt hole path builds as before (W1, xl entry 32's module, build-only).** `test_wf12_f_human_copper_degree2.py`'s (A0): hole cells == `HOLE_CELL_RECORD` (369 387; relative 0 at `20261008T093243Z_WF-12-step4a-W1.log:6220`) within `CELL_COUNT_BAND`, plus the module's (Deg) check, green as landed. The build goes through `_build_f_human(comm, as_hole=True)` after the edit. It must report **369 387 exactly**, as on 10-08: same image, and the kwargs are identical by (K).
- **(R1) P17 reproduction on the switched default (W2, degree 1, `-n 8`).** `P_src,exact` and pooled κ within `REUSE_REPRODUCTION_RTOL` of the re-recorded values. This is run-to-run on the same image (10-07 measured, W2 re-measures). It is not self-referential.
- **(R2) Through the file (W2).** Saved to the new `…-noleg.npz`:
  - (L1) round trip;
  - σ_max and all 18 class means within `REUSE_REPRODUCTION_RTOL` of the re-recorded values;
  - the save control ≥ `SAVE_CONTROL_FACTOR`× the band, as landed.
- **(I) `WF-7` step 1's identities at their unchanged bands (W2).**
  - (i) cells within `CELL_COUNT_BAND` of 504 642 (predicted −3.840e-03);
  - (ii) reciprocity ≤ `RECIPROCITY_BAND`;
  - (iii) σ_max ≤ `MATRIX_PASSIVITY_CEILING`;
  - (iv) exactly `N_SYMMETRY_CLASSES` (18), residual < `AZIMUTH_MATCH_DEG`, worst spread ≤ `OPPOSITE_SPREAD_BAND`;
  - (v) the exact identity ≤ `DISCRETE_IDENTITY_RTOL` on all 32;
  - the column-scale control > `CONTROL_RATIO_LIFT`×;
  - **(N′)** leg census 0.
  
  Step 1's W2 already had (ii)–(v) green on identical kwargs (`:9735–9756`).
- **(Q) `WF-7` step 2 on the switched `_build_f_human` (W3).** These are the Done-when's "C16, mirror, exact identities":
  - worst C16 image ≤ `C4_COVARIANCE_BAND` (5 %);
  - mirror ≤ `C4_COVARIANCE_BAND`;
  - exact identity ≤ `DISCRETE_IDENTITY_RTOL`;
  - unit weight ≤ `UNIT_WEIGHT_CONTROL_RTOL`;
  - (i) cells within `CELL_COUNT_BAND` of 504 642.
  
  Boxed readings: C16 0.1495 %, mirror 0.1634 %, identity 1.390e-15 (`20261002T021044Z_WF-12-step1-ruleA-WF-7-step2.log:11567–11569`).

**Negative controls (rule (e)).**
- **(c1) The pin and the exemption are not silently ignored.** *Asserted* (W0, no mesh): (4)'s dict checks separate "key present, `False`" (default) from "no key" (pinned or `as_hole`). This is a binary separation; no factor is involved.
- **(c2) The 1e-6 comparator tells the meshes apart.** *Asserted* (W2): on the switched run, the **boxed** P17 `P_src,exact` (`*_BOXED`) misses by more than `REUSE_REPRODUCTION_RTOL`, and so does the worst class mean against `CLASS_MEAN_D1_RECORD_BOXED`. Predicted ≈ 77× (7.706e-05, `:9733`) and ≈ 910× (9.100e-04, `:9780`). The ceiling is set by the mesh move itself, so the bar is "> 1× the band". This shows (R1) / (R2) are not passing on an unchanged (still boxed) fixture.
- **(c3)** `WF-7`'s omit-one-sheet identity control, *predicted, printed*: ≈ 480× the band, as landed (`…W2-wf7.log` prints it on the landed route).
- **(c4)** Step 2's landed controls, *predicted, printed*: the cw mis-paired reading ≈ 99.3 % and the single-drive P17 C16 spread ≈ 33 % (`…ruleA-WF-7-step2.log:11570`, `:11573`).

**Printed (not asserted).**
- W2: cells 502 704 beside `WF7_BOXED_CELL_RECORD`.
- W3: `S[P17,P17]` beside the re-recorded `S_DRIVEN_STEP1_PRINTED`. Predicted equal to all 6 printed digits: same drive, same mesh, `-n 8`.
- W3: the C16 / mirror moves against the boxed readings above.
- Build seconds and summed `ru_maxrss` in every window.

**Windows.** Run every window in the foreground. The harness call is the **whole** Bash command, bare and repo-relative. Bash `timeout: 660000`; `-s`; node IDs, never `-k`; `tests/environment` first. Every window is under 10 min, so no redirect is needed. Collected-count identity: W2 = step 1 W2's "16 passed, 1 skipped" (`…W2-wf7.log:9847`) + 1 (the new knob-contract node); W3 = `tests/environment` + the 5 step-2 nodes (`…ruleA-WF-7-step2.log:11610`).
- **W0, the contracts.** Standard, `-n 2`, no mesh. Predicted ≤ 60 s (step 1 W0 was 76 s including a port-tags build, `20261007T004148Z_GEO-35-s1-W0-porttags.log:2039`).
  `scripts/testing/run_and_log.sh GEO-35-s2a-W0-contracts "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_geo35_leg_port_boxes.py::test_geo35_knob_contract tests/mesh/test_geo35_leg_port_boxes.py::test_geo35_f_human_builder_contract tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_save_knob_contract tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_orientation_knob_contract tests/validation/test_wf7_f_human_ring_matrix.py::test_wf7_leg_port_knob_contract -v -s --tb=short'"`
- **W1, the exempt hole path (xl entry 32's module, build-only).** Heavy by ceiling, `-n 2`. Measured 184 s (`20261008T093243Z_WF-12-step4a-W1.log:6359`; build 114.6 s at `:6220`). This is the regression for the `_build_f_human` edit, rule (c).
  `scripts/testing/run_and_log.sh GEO-35-s2a-W1-wf12-buildonly "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF12_STEP4=1 FEM_EM_WF12_STEP4_BUILD_ONLY=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf12_f_human_copper_degree2.py::test_wf12_step4_copper_degree2_column -v -s --tb=short'"`
- **W2, `WF-7` step 1 on the switched default, with the re-recorded constants** (`record-reconciler`'s anchor). Heavy, `-n 8`. Measured 461 s on these exact kwargs (`20261007T005016Z_GEO-35-s1-W2-wf7.log:10306`; body 385.9 s, 10.756 GiB at `:9782`); bracket 430–520 s.
  `scripts/testing/run_and_log.sh GEO-35-s2a-W2-wf7-step1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf7 && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_ORIENTATION_CONTROL=split FEM_EM_WF7_SAVE_NPZ=/workspace/logs/wf7/wf7-s32-degree1-n8-noleg.npz timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_ring_matrix.py -v -s --tb=short'"`
- **W3, `WF-7` step 2 on the switched `_build_f_human`.** Heavy, `-n 8`. Measured 221 s without `tests/environment` (`20261002T021044Z_WF-12-step1-ruleA-WF-7-step2.log:11636`); predicted 240–280 s.
  `scripts/testing/run_and_log.sh GEO-35-s2a-W3-wf7-step2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_STEP2=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_b1_quadrature.py -v -s --tb=short'"`

**Traps already paid for.**
- **The hole refusal** (`mesh.py:4272–4285`). Forwarding the new default under `as_hole` raises on every copper caller, and on xl entry 32 at 02:00 on 10-13 with no session to see it. (K) and W1 are the checks.
- **"Adds no key" is the byte-identity proof** (the `phantom_resolution` precedent). A pinned `emit_leg_ports=True` key is harmless to gmsh but defeats the dict check. Omit the key.
- **Overwriting the boxed npz files.** Every old `WF-7` command in `WF-7.md` (`:66`, `:168`, `:251`) names `wf7-s32-degree1-n8.npz`. The guard refuses it. W2 saves to `…-noleg.npz` only.
- **A self-referential re-record.** New values come from the 10-07 footered log, never from W2. If W2 is red, do not copy W2's reading into the constant.
- **Rank-safety.** `cell_tags.values` is rank-local, so the (N′) census is allreduced (step 1's helper). A failing gmsh build deadlocks at rank > 1 (`GEO-23`), so the helper's refusal raises before gmsh.
- **Orphan `.c`.** The degree-1 forms are warm from step 1 and H15a. If a window is killed: `pgrep -c python3` first, then the rubric's sweep. Never use `pgrep -f "python3 -m pytest"`.
- **Status 124 is a finding, not a retry.** Never lengthen a timeout or change `-n`.
- **Rule (i):** no edit to the helper, the pins, the knob or a constant after its last green window.
- **Exact rectangles and κ (operator ruling 2026-10-02):** no sheet is touched, κ stays a diagnostic (`kappa_applied=False`), and the κ record is a run-to-run reading, not a correction.
- **Write ≈, never `~`.** No AED file is read or touched. `ANS-8` is only pinned.

**Scope.**
- **Covered:** the σ = 800 solid F-human through `WF-7` step 1's own build and `_build_f_human`'s non-hole default; `WF-7` steps 1 and 2 at degree 1, 64 MHz, `-n 8`; the exempt hole build (xl entry 32's module, build only); the save guard and the boxed rebuild knob.
- **Not claimed, and what remains (2b and the review):**
  - **2b** (one slot, predicted): unpin `WF-7` step 3 / `EX-62` and re-record the printed 10 g figures. These are `mri:4`'s guide (537 546 cells, `SAR10g(c*)`, C16 0.4904 %) and §2's citation of `20260925T093315Z_WF-7-step3.log:11615–11617`. Measured 352 s (`20260925T093315Z_WF-7-step3.log.gz:11686`) + `EX-62` 287–323 s (`20260926T095603Z_EX-62.log.gz:2214`, `20260930T213057Z_EX-62-fields.log:11613`) + census. 2b also covers `GEO-25` / `mesh:12`: `F_HUMAN_BRANCH_B_CELL_RECORD`'s knob-on *transverse* count is unmeasured (predicted ≈ −0.9 %, at the edge of `CELL_COUNT_BAND`), and `_build_rung`'s partition must take a ring-only tag list (step 1's absent-tag trap). Measured 196 s (`20260906T034048Z_EX-51-gate-rerun.log.gz:61`).
  - **`WF-10` / `ANS-8`** (Q3): κ-corrected H8-lineage records (`C_TUNED_F`, `C_MODE_F`, `A1/A2_MEAS`, `ANTI_MODE_C_PF`; `ANS-8`'s item-133 circuit records at 1e-6). Measured 334 + 423 + 405 + 292 s. That exceeds one slot, and `C_MODE_F` also feeds `test_wf12_f_human_copper_mechanism.py:215`.
  - **`mesh:9`–`11`** (Q2): r = 0.07 m GEO-20 / GEO-26 fixtures, not F-human. They are shared with `PORT-13`'s ring matrix and column, `POST-6`, `ports:10`–`12` and `mesh:7`.
  - **The hole and quadrant layouts:** out of U21 (Settled (a)).
  - No degree 2 and no absolute S.

**Status it can move.** No glyph. **The §7 Done-when cannot be met by 2a**: it needs every moved record re-recorded, which 2b and the rulings Q2 / Q3 define. That mismatch is flagged here. On green:
- `GEO-35` holds 🟡 with a step 2a record sentence: switched fixture, the six re-recorded constants with their moves, and (I) / (Q) / (N′) green;
- `WF-7` ✅ gains one sentence (records now on the leg-box-free fixture; boxed values kept and rebuildable);
- §10 U21 is annotated "2a DONE, 2b open";
- U27 is unblocked.

**Negative result.** In every case: report the measurement, park on `attempt/GEO-35-step2a-*`, mark item 351 BLOCKED in the same commit (rule (d)), and journal it in attempts.md.
- **W0 red:** fix in-slot only inside the new helper, knob or pins; otherwise stop.
- **W1 red, or a cell count ≠ 369 387:** the exemption leaked. Revert (1); BLOCKED. Entry 32 must not meet this commit.
- **(R1) / (R2) red with (I) green:** the 10-07 readings do not reproduce on the same image at 1e-6. That is a run-to-run finding (known-issues); never re-record from W2.
- **(c2) not separating, or (N′) > 0:** the switch did not take effect. Stop.
- **(I) or (Q) red:** an identity broke on the switched mesh. It is a real finding: revert the switch, and `WF-7` stands on the boxed fixture.
- **W2 over minute 42, or a window killed at 590 s:** BLOCKED without a result; the review re-poses.

No `PREDICTED RED:` line.

**For the review (rule before queueing).**
- **Q1 (binding).** `record-reconciler`'s drift sanity gate stops at ≈ 0.5 % (`.claude/agents/record-reconciler.md`, "Drift sanity gate") and was written for mesher-version drift. Two 2a sites move by more, because of geometry: pooled κ +4.671 % and the printed cell constant −0.899 %. Does `GEO-35`'s Done-when (weekly-placed) license re-recording these, with step 1's measured moves as the pre-registered reference and W2's 1e-6 reproduction as the check? If not, do not queue 351. The κ assert would be red on the switched default, and demoting it would be loosening.
- **Q2.** Should `mesh:9`–`11` be dropped from U21 (they keep the boxes and their guides say so), or opened as their own step?
- **Q3.** Should the `WF-10` / `ANS-8` κ-era records be re-recorded (≈ 2 slots), or frozen as records of the boxed fixture with a permanent pin (the analogue of binding ruling (iii))?
- **Q4.** Ratify the 2a / 2b split of U21.

**Review rulings (2026-10-09 18:00 review) on item 351 — these outrank the brief above.** (Q1) **Licensed.** `GEO-35`'s §7 Done-when says that "every record the switch moves … is re-recorded under the (1*) licence by `record-reconciler`, with the gated anchors (C16, mirror, exact identities) re-run green at their unchanged bands". That is an explicit licence for this geometry change. It outranks `record-reconciler`'s ≈ 0.5 % mesher-drift stop for these records only. **The licence is not open-ended.** The new record values are pre-registered by step 1's knob-on measurement (`20261007T005016Z_GEO-35-s1-W2-wf7.log`, the `:9694–9779` lines the brief cites; cells 502 704, pooled κ 4.805416e-03 and the rest). Each re-recorded value must reproduce that reading at `REUSE_REPRODUCTION_RTOL`. A value that does not reproduce is a red: stop and record it, with no re-record. No band moves. κ stays a diagnostic (operator ruling 2026-10-02); its record is a reproduction record, not a gate. (Q2) **`mesh:9`–`11` are not touched by 2a.** U21's text names `mesh:9`–`12`. Whether the 0.07 m fixtures belong in U21 is the weekly's call, and it is flagged on the dashboard. (Q3) **2a pins the κ-era `WF-10` / `ANS-8` records to the boxed mesh** (`emit_leg_ports=True` explicit). Whether they are re-recorded or frozen as boxed-mesh records is the weekly's call. (Q4) **The 2a / 2b split is provisional (rule (6)); the weekly ratifies or reverts it.** `GEO-35` holds 🟡 on 2a, because 2a cannot meet the Done-when. The H13b remark is already ruled in `TH-20.md` (item 350). (5) **Ordering:** item 350 (`TH-20` 2a) is independent of this item. Item 351 runs after 350 in the queue only to keep the heavy windows serial. xl entry 32 (10-13) depends on W1's build-only reading staying 369 387 cells. **If W1 is red, stop before W2 and mark the item BLOCKED.** Entry 32 must not run on a broken hole path, so the slot also writes "HOLD — item 351 W1 red" into entry 32's status line.

### §9 item 356 — 2026-10-10 (chain step U21, 2b-i, provisional)

*Drafted for the 2026-10-10 03:00 review.* Licence: U21 under rule (6). The 10-09 Q1 ruling allows re-recording only values that a prior knob-on run pre-registered. This item makes those runs and **re-records and unpins nothing**. The weekly ratifies the 2b-i / 2b-ii split tonight. Two windows would double `WF-7` step 3's cost, so the work is **split**: 356 covers `WF-7` step 3 and `EX-62`; 357 (below) covers `GEO-25` and `mesh:12`. No band, record or Done-when moves.

- **Executor:** implementer. **Tier:** W0 standard `-n 2`; W1 and W2 heavy `-n 8`, the record width. **Cost:** ≈ 12 measured (2 × 352 s, `20260925T093315Z_WF-7-step3.log.gz:11686`) / ≈ 40 predicted = 15 fixed + 8 code + W0 1.5 + W1 6.5 + W2 6.5 + 2.5 spare. **Park if W1 is not green by minute 26.**
- **Code** (test-only and additive): changes to `tests/validation/test_wf7_f_human_sar10g_hotspot.py`.
  - (1) `build_sar10g_case(comm, *, emit_leg_ports=True)` is passed through at the pinned call (`:154–156`). The comment stays, and `mri:4` (`04_f_human_sar10g_hotspot.py:236`) is unchanged.
  - (2) New knob `FEM_EM_GEO35_2B_MEASURE` ∈ {unset, `boxed`, `dropped`}, read in the fixture (`:409–416`). `dropped` passes `False`; any other value raises. Unset adds nothing.
  - (3) New node `test_geo35_2b_measure`, which skips when the knob is unset, and a no-mesh `…_knob_contract`.
  - (4) Two provenance constants, used only under the knob: 537 546 cells (`…step3.log.gz:10468`) and SAR10g(c\*) 2.085815849e-04 (`:11614`).
  - (5) The leg census imports `_owned_census` (`tests/mesh/test_geo35_leg_port_boxes.py:71`) and allreduces it.
- **Anchors (asserted).** These hold in W1 and W2, with every band imported:
  - (i) coverage power and mass ≤ `EXACT_IDENTITY_RTOL` (`test_birdcage_sar_mass_averaged.py:163`);
  - (ii) worst C16 and (iii) mirror ≤ `C4_COVARIANCE_BAND` (`test_birdcage_b1_plus_map.py:120`);
  - (iv) every ball is contained.

  **W1, `boxed`:** the pinned route is the control. Cells must equal 537 546, and SAR10g(c\*) must match 2.085815849e-04 within `REUSE_REPRODUCTION_RTOL`. Import that constant from `test_port_birdcage_ring_matrix.py:171` (1e-6), **not** from `test_port19_factor_reuse.py:82` (1e-12). This window also answers the 2a auditor's note for this caller. **W2, `dropped`:** (N′) the leg census is 0.
- **Negative controls** (rule (e); the analogue is 2a's (c2) at 77.1× / 910.0×).
  - (c1) *Asserted* in W1: all 16 leg tags have more than 0 cells.
  - (c2) *Asserted* in W2: cells ≠ 537 546, and the relative change of SAR10g(c\*) from the boxed value is > `REUSE_REPRODUCTION_RTOL`. Predicted 10–3 000×, from step 1's moves of 3.3e-06 … 9.1e-04 (`…GEO-35-s1-W2-wf7.log:9733–9780`). The ceiling is the mesh move itself, so the bar is > 1×.
  - (c3) *Printed:* the P17 control, predicted ≈ 96 %, ≈ 19× the band (`:11619`).
- **Printed values are the pre-registration** (W2), each beside its boxed value:
  - cells, predicted ≈ −0.85 % (range [−1.0 %, −0.7 %]);
  - phantom cells, predicted 22 594 (the phantom boundary is unchanged);
  - c\*, SAR10g(c\*) and the nine-ball table;
  - the 15 C16 values (boxed worst 0.4904 % at R_7);
  - mirror (0.0502 %), mean SAR, the P17 control, and the |B1+| line.

  Round-off residuals (coverage ≈ 1e-13) are band checks, never records, because they jitter between runs (`mri:4` guide §3).
- **`EX-62` gets no window.** `mri:4` prints the values `build_sar10g_case` returns. It matched the gate to every printed digit (`20260930T213057Z_EX-62-fields.log:11571–11584`), and 2b-ii re-runs it. **`WF-10` ×2 and `08` are not worth a window** (334–423 s each, and Q3 freezes them). Their pins resolve to the no-key dict that `test_geo35_f_human_builder_contract` asserts (`test_geo35_leg_port_boxes.py:180–181`), and W1 runs the same pinned line end to end.
- **Windows.** Run each in the foreground as the whole Bash command, with Bash `timeout: 660000`. Use `-s` and node IDs, never `-k`, and put `tests/environment` first. Each window is under 10 min, so no redirect is needed. Heavy windows never run side by side.
  `scripts/testing/run_and_log.sh GEO-35-s2bi-W<k> "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 <ENV> timeout -k 30 <T> mpiexec -n <N> python3 -m pytest tests/environment <NODES> -v -s --tb=short'"`
  - **W0:** T 180, N 2. Nodes: the knob contract and `test_geo35_leg_port_boxes.py::test_geo35_f_human_builder_contract`.
  - **W1:** T 590, N 8, `FEM_EM_WF7_STEP3=1 FEM_EM_GEO35_2B_MEASURE=boxed`, the whole module (6 nodes). Predicted 330–420 s.
  - **W2:** the same as W1 with `dropped`.
- **Traps.** False is reachable only through the knob. There are two `REUSE_REPRODUCTION_RTOL` constants (see Anchors). `cell_tags.values` is rank-local, so allreduce. If a window is killed, run `pgrep -c python3` first, then the orphan `.c` sweep. Status 124 is a finding: never lengthen the timeout. Rule (i): no edit after the last window. No AED file, sheet or κ is touched.
- **Scope.** Measurement only.
  - Nothing is re-recorded, unpinned or moved.
  - `GEO-35` holds 🟡.
  - `mesh:9`–`11` are untouched (Q2), and so are `WF-10` / `ANS-8` (Q3).
  - No degree 2 and no absolute SAR.
- **Negative result.** In every case: report, park on `attempt/GEO-35-step2bi-*`, mark the item BLOCKED and journal it.
  - **W0 red:** fix inside the knob or node only.
  - **W1 red:** the image has drifted since 09-30. Add a known-issues entry, raise a `record-reconciler` finding, and stop before W2.
  - **W2 (i)–(iv) red:** the switch breaks a 10 g identity. This finding holds 2b; step 3 stays boxed.
  - **(N′) > 0, or (c2) does not separate:** the switch did not take effect, or the SAR figures do not move at 1e-6 (then 2b-ii re-records cells only). Stop.

  No `PREDICTED RED:` line.
- **Status it can move.** No glyph changes.
  - GEO-35.md records W2's table, with log:line citations, as 2b-ii's pre-registered values.
  - The `WF-7` row notes the 10 g C16, mirror and coverage checks green on the leg-box-free mesh.
  - §10 U21 reads "2b-i (356) DONE".

**Outcome (2026-10-10, 06:00 CDT slot, 11:00–11:25 UTC): DONE.** Code: `build_sar10g_case(comm, *, emit_leg_ports=True)`, knob `FEM_EM_GEO35_2B_MEASURE` (`_geo35_2b_mode`), nodes `test_geo35_2b_measure` and `test_geo35_2b_measure_knob_contract`, provenance constants `GEO35_BOXED_CELLS` / `GEO35_BOXED_SAR10G_STAR`, all in `tests/validation/test_wf7_f_human_sar10g_hotspot.py`. Windows, all green: **W0** `-n 2`, 14 passed, 78 s (`20261010T110132Z_GEO-35-s2bi-W0.log`); **W1** `boxed`, `-n 8`, 18 passed, 419 s (`20261010T110256Z_GEO-35-s2bi-W1.log`); **W2** `dropped`, `-n 8`, 18 passed, 419 s (`20261010T111003Z_GEO-35-s2bi-W2.log`).

- **W1 control (asserted):** cells 537 546 = the 09-25 record; SAR10g(c\*) 2.085815849e-04, rel 1.094e-10 against `REUSE_REPRODUCTION_RTOL` 1e-6; (c1) the 16 lower leg tags hold 120–187 cells each (2 259 in total; the upper tags read 0 on the boxed build too, as in the 2a module, where (c1) is asserted on the lower tags only); worst C16 0.4904 % at R_7, mirror 0.0502 %, coverage power 9.448e-14 / mass 1.297e-13, P17 control 96.4015 % — every 09-25 figure reproduced (`…W1.log:11746–11755`). This also gives the pinned `_build_f_human` call a runtime witness (2a audit caveat, 10-10 ruling (4)).
- **W2, the 2b-ii pre-registration** (`…W2.log:10800–10834`), each beside its boxed value:

  | quantity | dropped (W2) | boxed (W1) | brief's prediction |
  |---|---|---|---|
  | cells | 532 935 (−0.8578 %) | 537 546 | −0.85 % in [−1.0, −0.7] % — **in range** |
  | phantom cells | 22 822 | 22 594 | 22 594 — **missed** (+1.01 %); the phantom boundary is unchanged but its interior mesh is not |
  | leg census (N′, asserted) | 0 over 32 tags | 2 259 | 0 |
  | c\* | r 45 mm, z = z_c (44.135, 8.779, 0.000) mm | same | — |
  | SAR10g(c\*) | 2.079992000e-04 W/kg | 2.085815849e-04 | — |
  | (c2) SAR move / `REUSE_REPRODUCTION_RTOL` | 2.792e-03 = 2792× (printed, ruling (2)) | — | 10–3 000× — **in range** |
  | worst C16 | 0.3373 % at R_3 | 0.4904 % at R_7 | ≤ 5 % (asserted) |
  | mirror | 0.0234 % | 0.0502 % | ≤ 5 % (asserted) |
  | coverage power / mass | 1.121e-14 / 1.645e-13 | 9.448e-14 / 1.297e-13 | ≤ 1e-10 (asserted); jitter, never a record |
  | mean SAR (ccw) | 1.735495022e-04 W/kg | 1.735606450e-04 | — |
  | P17 control spread | 96.4033 % (285.8× the ccw worst) | 96.4015 % | ≈ 96 % |
  | \|B1+\| C16 / mirror / valid | 0.0480 % at R_6 / 0.0392 % / 213 of 213 | 0.0513 % / 0.0456 % / 215 of 215 | — |

  Nine-ball table (W2, ccw, mm, W/kg): z_c − 50: r 15 7.463600451e-05, r 30 9.595433358e-05, r 45 1.441256449e-04; z_c: 2.889775125e-05, 9.425441341e-05, **2.079992000e-04 (c\*)**; z_c + 50: 7.451704145e-05, 9.627011835e-05, 1.438777841e-04 (`:10800–10808`). The 15 C16 values are at `:10809–10823`.
- **Reading:** the switch moves SAR10g(c\*) by 2.8e-3. That is above the 1e-6 bar, so 2b-ii re-records the SAR figures and not only the cells (ruling (2)). Two unpredicted moves: the \|B1+\| sample set drops from 215 to 213 valid points, and the phantom cell count moves +1.01 % against a predicted 0. Both are printed and neither is asserted. Nothing was re-recorded or unpinned; `GEO-35` holds 🟡.

### §9 item 357 — 2026-10-10 (chain step U21, 2b-i part 2, provisional)

The `GEO-25` / `mesh:12` transverse build. It shares no file with 356.

- **Executor:** implementer. **Tier:** heavy by ceiling, mesh only, `-n 2`, the width the logs used (`20260906T034048Z_EX-51-gate-rerun.log.gz:12`). **Cost:** ≈ 7 measured (2 × 196 s, `:61`) / ≈ 31 predicted = 15 fixed + 6 code + W1 4 + W2 4 + 2 spare.
- **Code.** `_build_rung` (`test_birdcage_f_human_rung.py:139–191`) gains the keyword `emit_leg_ports=True`.
  - The default adds no key and leaves the tag list unchanged. `mesh:12` (`:179`, `:228`) is untouched.
  - `False` adds the key and uses a **ring-only** tag list. The leg entries at `:163` would otherwise let the partition pass vacuously.
  - New node `test_geo35_transverse_rung_measure` in `test_geo35_leg_port_boxes.py` runs branches B and A knob-on.
- **Anchors (asserted, bands imported).** These apply to both knob-on builds:
  - the partition is within `EXACT` (`test_birdcage_ring_gaps.py:87`);
  - all 32 terminal ratios are in `TERMINAL_RATIO_BAND`;
  - the `CAD_MASS_GATE` separation holds (B above the gate, A below it);
  - (N′) is 0;
  - cells are within `CELL_COUNT_BAND` of each landed record.

  **W1, the control and the rule (c) regression:** the unchanged committed module, expected to give exactly 504 642 / 204 977 (`20260906T033712Z_EX-51.log.gz:20003`, `:20009`).
- **Negative controls.** *Asserted:* each knob-on count ≠ its boxed record. *Printed:* a predicted move of ≈ −0.9 % for B (step 1 moved −0.899 %); A is unpredicted.
- **Printed:** the distance of every consumer of `F_HUMAN_BRANCH_B_CELL_RECORD` from the knob-on B count. Boxed 507 266 is asserted at `test_wf7_f_human_ring_matrix.py:1092`, `test_wf12_f_human_hole.py:380` and `test_geo35_leg_port_boxes.py:243`; it is predicted at ≈ +1.4 %, **outside the band**.
- **Windows.** W1 runs `tests/validation/test_birdcage_f_human_rung.py`; W2 runs the new node. Both use T 590 and N 2. The traps, scope and dispositions are 356's. An out-of-band knob-on count is reported, never widened.
- **Status it can move.** The transverse knob-on counts are pre-registered for 2b-ii.
- **Q5 (for the review).** `F_HUMAN_BRANCH_B_CELL_RECORD` is the shared 1 % anchor of the longitudinal fixtures. Re-recording it in place would turn three boxed gates red. 2b-ii needs either a separate leg-box-free record or re-pointed consumers.

**Outcome (2026-10-10, 06:00 CDT slot, 11:19–11:31 UTC): DONE.** Code: `_build_rung(..., emit_leg_ports=True)` in `tests/validation/test_birdcage_f_human_rung.py` (default adds no key and keeps the tag list; `False` adds the key, uses a ring-only tag list, and `keep_mesh=True` now also returns `diag`); new node `test_geo35_transverse_rung_measure` in `tests/mesh/test_geo35_leg_port_boxes.py`. Both windows `-n 2`, T 590, green.

- **W1 control** (`20261010T112008Z_GEO-35-s2bi-p2-W1.log`, 10 passed, 302 s): the committed module reproduces 504 642 / 204 977 **exactly** (`:9767`, `:20054`, rel 0.000e+00); partition 1.000000000000, CAD 0.965414 / 0.893028, separation 0.072387 (`:20060–20063`). The module asserts the records at `CELL_COUNT_BAND`; exact equality is read off the printed lines.
- **W2 knob-on** (`20261010T112516Z_GEO-35-s2bi-p2-W2.log`, 9 passed, 294 s), the 2b-ii pre-registration:

  | quantity | branch B knob-on | branch A knob-on |
  |---|---|---|
  | cells (`:8771`, `:17802`) | 500 425 vs 504 642, −0.8356 % (predicted ≈ −0.9 %; in band, printed) | 203 945 vs 204 977, −0.5035 % (unpredicted; in band, printed) |
  | (N′) leg census, 32 tags (asserted) | 0 | 0 |
  | partition ring-only (asserted, `EXACT`) | 1.000000000000 | 1.000000000000 |
  | terminal ratio, 32 ports (asserted) | 0.974454791 – 0.974455230 | 0.974445583 – 0.974446474 |
  | meshed/CAD conductor (separation asserted) | 0.965414 ≥ 0.95 | 0.893028 < 0.95 |
  | build_s | 130.43 | 73.37 |

  Negative control asserted: each knob-on count ≠ its boxed record. The terminal ratios and CAD ratios are digit-identical to the boxed build (the leg boxes sit away from the ring terminals and the conductor).
- **Consumers of `F_HUMAN_BRANCH_B_CELL_RECORD`** (printed, `:17811–17813`): boxed 507 266 sits **+1.3670 %** from knob-on B 500 425, outside `CELL_COUNT_BAND` at all three sites (predicted ≈ +1.4 %, outside). This is Q5's measurement: re-recording the shared record in place to 500 425 would turn `test_wf7_f_human_ring_matrix.py:1092`, `test_wf12_f_human_hole.py:380` and `test_geo35_leg_port_boxes.py:243` red. Nothing was re-recorded; `GEO-35` holds 🟡.

**Review rulings (2026-10-10 03:00 review) on items 356 and 357 — these outrank the briefs above.** (1) **The split is accepted, provisionally** (rule (6)); 356 and 357 are independent and share no file. The weekly ratifies 2b-i / 2b-ii with the rest of U21 tonight. (2) **Item 356 (c2), rule (e):** no prior measurement backs the SAR10g(c\*) move on this fixture, so its "relative change > `REUSE_REPRODUCTION_RTOL`" clause is **predicted and printed** (10–3 000×) beside the measured value, not asserted. The cell inequality ≠ 537 546 and (N′) = 0 stay asserted: they are structural and witness that the switch took effect. A SAR move under 1e-6 is a reading for 2b-ii (then cells only are re-recorded), not a red. (3) **Item 357, the knob-on cell clause:** "cells within `CELL_COUNT_BAND` of each landed record" on the knob-on builds is **printed, not asserted**. Branch A's move is unpredicted, and an out-of-band count is the measurement Q5 needs, not a failure. The partition, terminal ratios, `CAD_MASS_GATE` separation and (N′) stay asserted, and so does W1's exact 504 642 / 204 977 control. (4) **`EX-62`, the `WF-10` pins and `08` get no window: accepted.** W1 of 356 runs the pinned `_build_f_human` line end to end; that answers the 2a audit's "no runtime witness" caveat for this caller. The boxed `FEM_EM_WF7_LEG_PORTS=boxed` rebuild route for steps 1–2 stays unexercised (an audit caveat for the weekly). (5) **Q5 goes to the weekly** with Q2 / Q3 (dashboard).

**Weekly ruling 2026-10-10 (U21, Q2 / Q3 / Q5), recorded here for the 2b-ii briefs below** (PROJECT_PLAN.md §10 rulings, "U21 (`GEO-35`), Q2 / Q3 / Q5 and the split"). The 2a / 2b split and 2b-i are ratified. **Q2:** `mesh:9`–`11` leave U21. They are the 0.07 m `GEO-20` / `GEO-26` fixtures and keep their boxes; U21 is F-human only. **Q3:** the κ-era `WF-10` / `ANS-8` records are frozen as boxed-mesh records with the explicit `emit_leg_ports=True` pin. **Q5:** the shared anchor is split, not moved. `F_HUMAN_BRANCH_B_CELL_RECORD` keeps the boxed value for its three boxed consumers (`B_RECORD_CONSUMERS`, `tests/mesh/test_geo35_leg_port_boxes.py:276–280`). A new `F_HUMAN_BRANCH_B_CELL_RECORD_UNBOXED`, item 357's pre-registered knob-on reading, anchors the switched consumers. **2b-ii** re-records `WF-7` step 3, `EX-62` and `GEO-25` transverse to it under the (1\*) licence (`record-reconciler`), and `GEO-35` → ✅ on its green anchors.

### §9 item 385a — 2026-10-11 (U21 2b-ii part 1: `WF-7` step 3 + `EX-62` switched and re-recorded; chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-11 03:00 review; the review reviews it.* Licence: the weekly ruling above, plus the 10-09 Q1 (1\*) licence (item 351 rulings). Every re-recorded value must reproduce its **pre-registered** 2b-i reading at `REUSE_REPRODUCTION_RTOL` (`tests/validation/test_port_birdcage_ring_matrix.py:171`, 1e-6; **not** `test_port19_factor_reuse.py:82`). A value that does not reproduce is a red: stop, re-record nothing. **No band, tolerance, Done-when or gate moves.**

**Why two items.** The weekly posed 2b-ii as one item. Its honest windows are W0 78 s, the `WF-7` step 3 module 419 s, `mri:4` 323 s, the `GEO-25` rung module 302 s and `mesh:12` 196 s: ≈ 1 320 s of serial compute, plus code and record text, which is over a slot's minute-35 window budget. The split follows the fixture families and shares no file. **385a: the refined longitudinal family (`WF-7` step 3, `mri:4`). 385b: the transverse family (`GEO-25`, `mesh:12`), the `…_UNBOXED` constants and the ✅.** 385b depends on 385a.

- **Executor:** `record-reconciler` (the weekly's named executor; its site list is below). **Class:** chain step (solve). **Tier:** W0 standard `-n 2`; W1 / W2 heavy `-n 8` (the record width). **Cost:** ≈ M measured / ≈ 48 predicted min (15 fixed + ≈ 8 code + W0 1.5 + W1 ≈ 7 + W2 ≈ 5.5 + ≈ 8 record text + censuses 1 + 2 spare). **Park if W1 is not green by minute 30.**
- **Depends on:** nothing open (items 356 / 357 DONE). Never beside another heavy window.
- **Price (cited):** 2b-i W2 (the same module under `dropped`) 419 s harness at `-n 8` (`20261010T111003Z_GEO-35-s2bi-W2.log`); W0 78 s (`20261010T110132Z_GEO-35-s2bi-W0.log`); `mri:4` 323 s (`20260930T213057Z_EX-62-fields.log.gz`, window 318.63 s). W1 predicted [380, 560] s at up to 1.32×, under `timeout -k 30 590`; W2 [290, 430] s.

**Site list (the reconciler's input; every other consumer is out of scope).**
1. `tests/validation/test_wf7_f_human_sar10g_hotspot.py`:
   - `build_sar10g_case(comm, *, emit_leg_ports=False)`: flip the default at `:134` and the comment at `:186`. This is the switch.
   - The knob stays: `boxed` passes `True` (the control), and `dropped` now equals the default.
   - `test_geo35_2b_measure_knob_contract` (`:454–478`): its "defaults to the pinned `True`" clause becomes `False`. That is the contract following the switch, not a gate.
   - **New provenance constants** beside `GEO35_BOXED_*` (`:82–83`), each commented with its log line: `GEO35_DROPPED_CELLS = 532935` and `GEO35_DROPPED_PHANTOM_CELLS = 22822` (`20261010T111003Z_GEO-35-s2bi-W2.log:10832`), and `GEO35_DROPPED_SAR10G_STAR = 2.079992000e-04` (the W2 nine-ball table, `:10800–10808`).
   - **New node `test_geo35_2b_default_reproduction`**, active when `FEM_EM_WF7_STEP3=1` and the knob is unset. It asserts: cells == `GEO35_DROPPED_CELLS`; phantom cells == `GEO35_DROPPED_PHANTOM_CELLS`; (N′) leg census 0 (`_owned_census`, allreduced, 2b-i's code); SAR10g(c\*) within `REUSE_REPRODUCTION_RTOL` of `GEO35_DROPPED_SAR10G_STAR`; c\* at the same lattice key.
   - The (i)–(iv) identity nodes are unchanged.
2. `examples/mri/04_f_human_sar10g_hotspot.py`: **no code edit** (it calls `build_sar10g_case(comm)`, `:236`, and so follows the default). `examples/mri/04_f_human_sar10g_hotspot.md`: the quoted records (`:120`, `:136`, `:145–148`: 537 546 / 22 594, SAR10g(c\*) 2.085815849e-04, C16 0.4904 % at R_7, mirror, mean SAR, the nine-ball table, the P17 control) are re-recorded **from W2's printed lines**, each of which must equal its 2b-i W2 value at the printed digits. Add one dated sentence: "leg port boxes dropped (`GEO-35` 2b-ii, 2026-10-11); the boxed figures are in `GEO-35.md`". Use the whole filename in `## Related` (rule (b)).
3. §7 rows `WF-7` and `EX-62`: one dated re-record sentence each, in the same commit, with log citations. If a row would pass 2 KB, the sentence goes in the chunk file with a pointer.
4. **Not touched (say so in the commit):** `F_HUMAN_BRANCH_B_CELL_RECORD` (385b adds the `…_UNBOXED` split); the three boxed 507 266 sites; `mesh:9`–`11` (Q2); the Q3 pins at `test_wf10_f_human_tuned_circuit.py:159, :454` and `08_f_human_quadrature_mode_content.py:193`. Those pins' "unpin with the re-record" comments are **superseded by Q3 and left unedited**: an edit, even to a comment, would re-open `WF-10` step 1's audited module under rule (i). This brief records the supersession.

**(1) Anchors (asserted, bands imported, unchanged).** W1, on the now-default dropped fixture: (i) coverage power / mass ≤ `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`); (ii) worst C16 and (iii) mirror ≤ `C4_COVARIANCE_BAND` (`tests/validation/test_birdcage_b1_plus_map.py:120`); (iv) containment; and the reproduction node. W2: `mri:4`'s gate-import assertions, as `EX-62` closed them.

**(2) Negative controls (rule (e)).**
- *Asserted, backed:* the reproduction node's cells ≠ `GEO35_BOXED_CELLS` (2b-i W2 read 532 935 vs 537 546); (N′) = 0 against 2b-i W1's boxed 2 259 lower-tag cells (`20261010T110256Z_GEO-35-s2bi-W1.log:11753`).
- *Printed:* the SAR10g(c\*) move from the boxed value, 2.792e-3 (2b-i W2); the P17 control ≈ 96.40 % (`EX-62`'s).

**(3) Windows.** Each is foreground, and each is the whole bare repo-relative harness call, with Bash `timeout: 660000`, `-s`, node IDs (never `-k`), and `tests/environment` first. Each is under 10 min, so the form is 2b-i's. Before W1 and W2, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
- **W0** (`-n 2`, T 180): `scripts/testing/run_and_log.sh GEO-35-s2bii-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_sar10g_hotspot.py::test_geo35_2b_measure_knob_contract tests/mesh/test_geo35_leg_port_boxes.py::test_geo35_f_human_builder_contract -v -s --tb=short'"`
- **W1** (`-n 8`, T 590, no knob): `scripts/testing/run_and_log.sh GEO-35-s2bii-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF7_STEP3=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_wf7_f_human_sar10g_hotspot.py -v -s --tb=short'"`
- **W2** (`mri:4`, flagged, so the setup PNG redraws without leg boxes; ≤ 600 KiB): first `scripts/testing/run_and_log.sh GEO-35-s2bii-mri4-dryrun "./scripts/run_examples.sh -e mri:4 -n 8 -t 590 --dry-run"` (it must print `(complex build)`), then the emitted docker line with `-T -e FEM_EM_SETUP_FIGURES=1` as its own bare harness call, labelled `GEO-35-s2bii-W2-mri4`. Then `GEO-35-s2bii-post-figures` and `GEO-35-s2bii-post-docrefs` (the `EX-109` census commands): `missing=0 broken=0`, no new dead reference.

**(4) Traps already paid for.** There are two `REUSE_REPRODUCTION_RTOL` constants: import from `test_port_birdcage_ring_matrix`. `cell_tags.values` and the leg census are rank-local, so allreduce. Never path-checkout a record file; append rows by hand and check that `git diff --stat` shows only insertions on record files. Rule (i): no edit after the window it closes on, and the guide is re-recorded before the post-censuses. Status 124: run `pgrep -c python3` first; only at 0, the orphan sweep (`GEO-35-s2bii-cache-sweep`). Never lengthen a timeout. No AED file, sheet or κ is touched. Write `≈`, never a tilde.

**(5) Scope.** The refined longitudinal family only. **Not claimed:** the transverse records and the `…_UNBOXED` constants (385b); `GEO-35` ✅ (385b); degree 2; absolute SAR; compliance.

**(6) Negative result.** Report; park on `attempt/GEO-35-step2bii-*`; mark 385a and 385b BLOCKED (rule (d)); journal.
- A reproduction red (cells, phantom cells or SAR10g(c\*) off its 2b-i value): an image or mesher drift since 10-10. File a known-issues entry and a `record-reconciler` finding; re-record nothing.
- (i)–(iv) red: the switch breaks a 10 g identity. Revert the default flip (step 3 stays boxed) and leave the finding for the weekly.
- A `mri:4` printed digit off W1 with every gate green: report it; not a stop.

No `PREDICTED RED:` line.

**(7) Status it can move.** No glyph by itself (`GEO-35` holds 🟡 until 385b). The `WF-7` and `EX-62` rows gain their re-record sentences; §10 U21 reads "2b-ii part 1 DONE".

**Executor spawn prompt (verbatim).** "Use the record-reconciler agent to execute §9 item 385a — `GEO-35` step 2b-ii part 1 (§10 U21) — per `docs/planning/chunks/GEO-35.md` §9 item 385a and the weekly ruling recorded above it. Its site list is the whole scope. Every re-recorded value must reproduce its pre-registered 2b-i value at `REUSE_REPRODUCTION_RTOL` from `test_port_birdcage_ring_matrix`; a miss is a red, not a re-record. Never touch a band, tolerance or gate. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W0, W1, W2 in order, foreground; before each heavy window run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. Never delete a log. Do not edit `tests/validation/test_wf10_f_human_tuned_circuit.py` or the ANS-8 example. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-35-step2bii-*`, mark items 385a and 385b BLOCKED in the same commit, and report."

**Outcome (2026-10-11, 21:00 CDT slot): DONE.** `build_sar10g_case` defaults to `emit_leg_ports=False` (knob `boxed` passes `True`; `dropped` equals the default); new constants `GEO35_DROPPED_CELLS` / `_PHANTOM_CELLS` / `_SAR10G_STAR` / `_KEY_STAR` with the 2b-i log lines; new node `test_geo35_2b_default_reproduction`. W0 14 passed, 81 s (`20261012T020117Z_GEO-35-s2bii-W0.log`). W1 `-n 8`, 18 passed 1 skipped (the knob-gated measure node), 424 s (`20261012T020246Z_GEO-35-s2bii-W1.log:10840`): cells 532 935, phantom 22 822, leg census 0, c\* r 45 mm dz 0, SAR10g(c\*) 2.079992000e-04 vs 2.079992000e-04 rel 9.864e-11 (0.00 × `REUSE_REPRODUCTION_RTOL`), move from boxed 2.792e-03; (i)–(iv) green. W2 `mri:4` flagged, 345 s (`20261012T021003Z_GEO-35-s2bii-W2-mri4.log`): every print equals 2b-i W2 (worst C16 0.3373 % at R_3, mirror 0.0234 %, P17 control 96.4033 %, mean SAR 1.735495022e-04, nine-ball table) except the coverage power round-off 1.965e-14 (2b-i read 1.121e-14; budget 1e-10; not a stop). Collected-count identity: W1 collected 19 (18 + 1 skipped); W0 14. Guide `examples/mri/04_f_human_sar10g_hotspot.md` re-recorded from W2's prints, boxed block kept as a record. Post-censuses `examples=106 ok=106 missing=0 broken=0`; docrefs `dead=0 guide=0 stale=85 exit=2`. **Not touched:** `F_HUMAN_BRANCH_B_CELL_RECORD`, the three boxed 507 266 sites, `mesh:9`–`11`, the Q3 pins in `test_wf10_f_human_tuned_circuit.py` and `08_f_human_quadrature_mode_content.py` (their "unpin with the re-record" comments superseded by Q3, left unedited), the ANS-8 example. **Flag for the review:** the `examples/mri/04_f_human_sar10g_hotspot.py` module docstring (line 44) still quotes the superseded boxed "196.6x the ccw worst" (now 285.8x); the brief said no code edit, so it was left. Boxed figures for reference: 537 546 / 22 594 cells, SAR10g(c\*) 2.085815849e-04, worst C16 0.4904 % at R_7, mirror 0.0502 %, P17 control 96.4015 % (`20260926T095026Z_EX-62.log`, `20260925T093315Z_WF-7-step3.log.gz`).

### §9 item 385b — 2026-10-11 (U21 2b-ii part 2: `GEO-25` transverse + `mesh:12` to `…_UNBOXED`; `GEO-35` ✅; chain step (solve))

*Drafted with 385a; same licence; no band, tolerance or gate moves.*

- **Executor:** `record-reconciler`. **Class:** chain step (solve). **Tier:** heavy by ceiling, mesh only, `-n 2` (the width of the records, `20260906T034048Z_EX-51-gate-rerun.log.gz`; `20261010T112008Z_GEO-35-s2bi-p2-W1.log`). **Cost:** ≈ M measured / ≈ 42 predicted min (15 fixed + ≈ 8 code + W3 ≈ 5.5 + W4 ≈ 4 + censuses 1 + ≈ 6 record text + 2 spare).
- **Depends on item 385a landing green; if it did not, skip.** It touches no file 385a touches.
- **Price (cited):** the rung module 302 s at `-n 2` (`20261010T112008Z_GEO-35-s2bi-p2-W1.log`); the knob-on branches B / A built in 130.43 / 73.37 s (`20261010T112516Z_GEO-35-s2bi-p2-W2.log`, 294 s); `mesh:12` 196 s at `-n 2` (`20260906T033712Z_EX-51.log.gz`).

**Site list.**
1. `tests/validation/test_birdcage_f_human_rung.py`:
   - beside `:87–88`, add **`F_HUMAN_BRANCH_B_CELL_RECORD_UNBOXED = 500425`** (`20261010T112516Z_GEO-35-s2bi-p2-W2.log:8771`) and **`F_HUMAN_BRANCH_A_CELL_RECORD_UNBOXED = 203945`** (`:17802`), each with a comment citing the weekly's Q5 and item 357;
   - `F_HUMAN_BRANCH_B_CELL_RECORD` 504 642 and `F_HUMAN_BRANCH_A_CELL_RECORD` 204 977 **keep their values** and gain a comment: "boxed; consumers `B_RECORD_CONSUMERS` and the 2b-i comparator";
   - `_build_rung(..., emit_leg_ports=False)`: flip the default at `:144`. This is the transverse switch;
   - the module's B / A record assertions (`:270`, `:345`) and their prints re-point to the `…_UNBOXED` constants.
   - *Open question for the review:* the weekly names only the B split. A's count also moves (−0.50 %) and the Done-when covers "every record the switch moves", so this brief splits A the same way. If the review declines it, A's assertion keeps 204 977, which still passes at `CELL_COUNT_BAND` (−0.50 %), and the row says A was not re-recorded.
2. `examples/meshing/12_birdcage_f_human_rung.py` (`:89–90`, `:235–281`): import and assert against the `…_UNBOXED` constants. The guide `examples/meshing/12_birdcage_f_human_rung.md` (`:29`, `:100`) quotes 500 425 / 203 945 from W4's printed lines and adds one dated sentence on the dropped leg boxes. The setup PNG is redrawn (flagged run, ≤ 600 KiB).
3. **Not touched:** `tests/mesh/test_geo35_leg_port_boxes.py`. Its `test_geo35_transverse_rung_measure` passes `emit_leg_ports=False` explicitly and compares against the **boxed** records by design (its ≠ control), so the flip does not change it. The three `B_RECORD_CONSUMERS`; `mesh:9`–`11` (Q2); the Q3 pins (385a item 4).
4. §7 rows `GEO-25` and `EX-51`: a dated re-record sentence each. The `GEO-35` row: the closure text below.

**(1) Anchors (asserted, bands imported).** W3, the rung module on the switched default: cells == `…_UNBOXED` exactly (an integer reproduction of 357's knob-on reading; asserted within `CELL_COUNT_BAND` as the module does, with exact equality read off the printed line, as item 357 did); partition within `EXACT` (`tests/validation/test_birdcage_ring_gaps.py:87`) over the ring-only tag list; 32 terminal ratios in `TERMINAL_RATIO_BAND`; the `CAD_MASS_GATE` separation (B 0.965414 above, A 0.893028 below). W4, `mesh:12`'s own assertions on the same quantities. The Done-when's mesh test (no leg-port tags, unchanged ring census on the switched fixture) is already met by step 1's (N) / (R) (`20261007T004311Z_GEO-35-s1-W1-mesh.log`) and item 357's transverse (N′) 0 (`…p2-W2.log`). Cite them; W3 is the switched-default regression.

**(2) Negative controls.** *Asserted, backed:* `test_geo35_transverse_rung_measure`'s knob-on ≠ boxed (357). It is not re-run, because its file is unchanged; cite `…p2-W2.log`. *Printed:* each `…_UNBOXED` beside its boxed record (−0.8356 % / −0.5035 %).

**(3) Windows.** The same rules as 385a.
- **W3** (`-n 2`, T 590): `scripts/testing/run_and_log.sh GEO-35-s2bii-W3-rung "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_birdcage_f_human_rung.py -v -s --tb=short'"`. This is the real build, as `EX-51`'s gate re-run was; no complex source.
- **W4:** `scripts/testing/run_and_log.sh GEO-35-s2bii-mesh12-dryrun "./scripts/run_examples.sh -e mesh:12 -n 2 -t 590 --dry-run"`, then the emitted docker line with `-T -e FEM_EM_SETUP_FIGURES=1` as one bare call labelled `GEO-35-s2bii-W4-mesh12`; then the post-figures and post-docrefs censuses (`…-post2-*`).

**(4) Traps.** 385a's. Plus: the ring-only tag list on the False path is required (item 357: leg entries would let the partition pass vacuously); never re-record a boxed constant in place (Q5: three boxed gates would go red at +1.367 %, `…p2-W2.log:17811–17813`).

**(5) Scope.** F-human transverse only. **Not claimed:** `mesh:9`–`11` (Q2); a re-record of the frozen κ-era records (Q3); degree 2.

**(6) Negative result.** Report; park on `attempt/GEO-35-step2bii-p2-*`; mark 385b BLOCKED; journal. A reproduction miss (a count ≠ 500 425 / 203 945): mesher drift, so file a known-issues entry and a reconciler finding, and re-record nothing. A partition or terminal-ratio red: the switch breaks a CAD identity on the default; revert the flip. No `PREDICTED RED:` line.

**(7) Status it can move.** **`GEO-35` 🟡 → ✅**, on 385a's and this item's green anchors plus step 1 / 357's mesh tests, if the Done-when matches. It does not match today: it names `F_HUMAN_BRANCH_B_CELL_RECORD` 504 642 among the records "re-recorded", and Q5 keeps that constant. **The row amendment (rule (j)) is in the review's return and must land before this item runs.** §10 U21 then reads DONE.

**Executor spawn prompt (verbatim).** "Use the record-reconciler agent to execute §9 item 385b — `GEO-35` step 2b-ii part 2 (§10 U21) — per `docs/planning/chunks/GEO-35.md` §9 item 385b; skip it unless item 385a is DONE. Its site list is the whole scope; boxed constants keep their values. Never touch a band, tolerance or gate. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W3, then W4, foreground. Never delete a log. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-35-step2bii-p2-*`, mark item 385b BLOCKED in the same commit, and report."

**Outcome (2026-10-11, 21:00 CDT slot): DONE.** `_build_rung` defaults to `emit_leg_ports=False`; new `F_HUMAN_BRANCH_B_CELL_RECORD_UNBOXED = 500425` and `F_HUMAN_BRANCH_A_CELL_RECORD_UNBOXED = 203945` (boxed 504 642 / 204 977 kept with a comment); the module's B / A record assertions, `mesh:12` and its guide re-pointed. W3 (`-n 2`, real build) 10 passed 4 skipped, 300 s (`20261012T021903Z_GEO-35-s2bii-W3-rung.log:8772,17804`): cells 500 425 / 203 945 exactly (rel 0.000e+00; −0.8356 % / −0.5035 % from boxed), partition 1.000000000000, terminal ratios min 0.974454791 / max 0.974455230 (B), CAD 0.965414 (+0.015414 above the 0.95 gate) / 0.893028 (0.056972 below), separation 0.072387. Collected-count identity: 14 against the census `20261010T112008Z_GEO-35-s2bi-p2-W1.log` (14, 10 passed 4 skipped). W4 `mesh:12` flagged, 230 s (`20261012T022429Z_GEO-35-s2bii-W4-mesh12.log:17761,17767`), every print equal to W3's, setup PNG redrawn (580 KiB, no leg boxes). Post-censuses `examples=106 ok=106 missing=0 broken=0` and docrefs `dead=0 guide=0 stale=84 exit=2` (age-based stale only). `GEO-35` row 🟡 → ✅ against the amended Done-when. **Not touched:** `tests/mesh/test_geo35_leg_port_boxes.py`, the three `B_RECORD_CONSUMERS`, `mesh:9`–`11`, the Q3 pins. Branch A's `…_UNBOXED` still awaits the 10-12 weekly's ratification. The `ans8` example's two imports of `F_HUMAN_BRANCH_B_CELL_RECORD` (boxed, Q3-frozen) are unchanged by design.

### Row state moved 2026-10-11 (03:00 daily review, §7 row ≤ 2 KB after the rule-(j) Done-when amendment)

🟡 **Step 2a (2026-10-10, §9 item 351, U21 split):** the solid F-human (`WF-7` step 1's own build and `_build_f_human`'s non-hole default) now drops the leg port boxes; the hole path stays boxed (369 387 cells exactly, `20261010T005723Z_GEO-35-s2a-W1-wf12-buildonly.log:6220`); `WF-10` / `WF-7` step 3 / `ANS-8` callers pinned `emit_leg_ports=True`. Six `WF-7` constants re-recorded under (1\*) from step 1's knob-on reading and reproduced at `REUSE_REPRODUCTION_RTOL` (P17 rel 4.3e-11 / κ 1.1e-9; σ_max 2.3e-10; worst class mean 4.8e-10; `20261010T010041Z_GEO-35-s2a-W2-wf7-step1.log:9750–9754`), (c2) 77× / 910× against the boxed records, (N′) 0, (ii)–(v) green; `WF-7` step 2 C16 0.1831 %, mirror 0.1264 %, identity 3.5e-16 (`20261010T010824Z_GEO-35-s2a-W3-wf7-step2.log:10708–10710`). Holds 🟡: 2b (`WF-7` step 3 / `EX-62`, `GEO-25`) and the Q2 / Q3 rulings remain. Step 1 (2026-10-06, §9 item 275, §10 U17): `emit_leg_ports=False` on F-human leaves zero leg-port cells with ring ordinals 17–48 and the `GEO-25` CAD identities intact (502 704 cells, −0.90 %), and the knob-on degree-1 32×32 holds its gates with every printed figure moving ≤ 9.1e-4 except pooled κ (4.67 %, a diagnostic). Step 2 — the fixture switch and the (1*) re-record — is open (§10 chain U). 

**Review ruling 2026-10-11 (03:00 daily) on items 385a / 385b.** The split is accepted (≈ 1 320 s of serial windows does not fit one slot). The provisional `F_HUMAN_BRANCH_A_CELL_RECORD_UNBOXED` is accepted under the Done-when's "every record the switch moves" and named for the 10-12 weekly to ratify (the weekly named branch B only). The Q3 pins' "unpin with the re-record" comments are not edited (an edit would re-open `WF-10` step 1's audited module under rule (i)); the supersession is recorded here.
