# GEO-34 — STEP import through Gmsh/OpenCASCADE with named-group tagging and a case config file

 **STEP import through Gmsh/OpenCASCADE with named-group tagging and a case config file** — the CAD-input path: complex designs come in as STEP, materials / boundaries / ports are declared on named CAD groups in one config beside the file, and the generator's own tag scheme comes out the other end (operator directive 2026-09-20: no native modeller, ever — CAD in, preflight view, solve) | 🟡 **Step 1a landed 2026-09-23** (`20260923T110405Z_GEO-34.log`): `io/step_import.import_step` + JSON group map re-import the Gmsh-OCC-written F-small STEP with the generator's tag census, all four sheet areas at `SHEET_AREA_M2` to ≤ 2.2e-16, conductor meshed/CAD 0.970165 ≥ 0.95, and both unmapped/empty-group controls raising; 116 662 cells vs the primitive 116 085 (printed, not gated). Steps 1b (`PORT-9` gates on the import) and 1c (independent kernel) open. **Step 1b attempted 2026-09-25 — negative result, parked** (`attempt/GEO-34-step1b-20260925T094541Z`, `20260925T094541Z_GEO-34-step1b.log`, 114 s standard `-n 2`): on the imported mesh the cyclic pre-gate (iv), reciprocity (i) 1.897e-14 and passivity (ii) σ_max 0.999992218 pass and the P2/P3-swapped control fails (iv) as asserted, but C4 gate (iii′) misses — self 3.3151 %, adjacent 1.1102 %, opposite 2.2713 % against 0.5 %, separation 3.0514× against 10×. The miss is two-fold (x-axis legs ≠ y-axis legs; narrowed sheets 26/26/29/26 facets), so the STEP round trip changes port-level physics beyond the bands. Leading hypothesis: the import loses the port-box triangulation congruence of `GEO-19` step B (known-issues 2026-09-25). **Rescoped 2026-09-26 03:00 review (family "step 1b": one attempt since 1a's gate; two more queued, so at most one letter remains after them).** Two independent items, bands unmoved: **1b-probe** (§9 item 46, `mesh-probe`, no solve) reads `GEO-31`'s sheet-congruence metrics — facet counts, per-facet area vectors, the centroid / boundary-vertex Hausdorff distance to sheet 1 rotated by 90°·(i − 1) — on the primitive ×1 rung and on the import side by side, so the known-issues entry's *Cause* cell gets a measurement; **1b′** (§9 item 47) is the decisive experiment and the candidate fix: an additive config key through which `import_step` binds the four imported gap sheets to sheet 1 with `gmsh.model.mesh.setPeriodic` under `_z_rotation_affine` — `GEO-32`'s congruent cut, which took the primitive's 10 MHz degree-2 self spread from 0.4451 to 0.0303 % (`ANS-4` 3c → 3f) — and then re-runs the parked module on it. A second mechanism is named for a reading of 1b′'s negative, not queued: the import's background field is the conductor Threshold alone (`io/step_import.py:224–245`), whereas the primitive's is a `Min` of the conductor Threshold and a phantom `Box` field (`io/mesh.py:4541–4604`), so the *volume* mesh differs too (116 662 vs 116 085 cells). **1b-probe landed 2026-09-26 (item 46, 🧪, `20260926T093153Z_GEO-34-step1b-probe.log`, 65 s standard `-n 2`, self-checks (A)/(B) green):** the import keeps the sheets' support and boundary (≤ 3.0e-14 m) but cuts them 66 / 58 / 66 / 58 facets against the primitive's 58 × 4 — a C2 count split matching the Z split, though P1 and P3 are not mutually congruent (d₃ 0.71 ℓ); readout (b) by the numbers, (a) excluded (known-issues 2026-09-25 *Cause* cell). **Step 1b ✅ 2026-09-26 via 1b′ (§9 item 47, `20260926T093715Z_GEO-34-step1bprime.log`, 186 s standard `-n 2`, complex, 19 passed `:7643`):** `import_step` gained the additive `congruent_sheets` config key (each lower/upper pair's post-heal shared face bound to P1's by `setPeriodic` under the snapped `_z_rotation_affine`; key absent ⇒ unchanged, rule (c) `20260926T094036Z_GEO-34-step1bprime-ruleC.log` 4 passed). On the congruent import (116 685 cells, printed): binding held — facet counts 66 × 4, centroid Hausdorff to rotated P1 ≤ 1.001e-14 m against 1e-12 (`:5698`); (iv) azimuths exact; (i) 4.286e-14; (ii) σ_max 0.999991299; **(iii′) self 0.0349 % / adjacent 0.0223 % / opposite 0.0131 %, separation 259.86×** (`:5687`), bands unmoved. Asserted negative control: 1a's config on the same STEP reproduces 1b's miss digit for digit (66/58/66/58, self 3.3151 %, `:7577–7578`); P2/P3-swap control fails (iv) at P2, P3 (`:5678`). **The sheet cut was the mechanism; the volume-sizing difference is not needed to pass.** Licence: the congruent cut is the import's port convention for `PORT-9`-class gates — the unconstrained import is not licensed. Rule (a) `20260926T094204Z_GEO-34-step1bprime-ruleA.log`: 116 085 cells, 0.0553 / 0.0353 / 0.0214 %, 166.6766× identical; ‖S − Sᵀ‖/‖S‖ read 1.065649875e-14 against the 1.044255156e-14 record — a round-off digit this fixture has printed in six earlier logs on unchanged code (8 distinct values across 28 logged readings), not moved by this item (the four-port module's code is 1b's). Steps 1c / 2 / 3 open. **Step 2 ✅ 2026-09-30 (§9 item 118, standard `-n 2`, real build, no solve):** `scripts/run_step_case.py STEP CONFIG --preflight --out DIR` → `io/step_preflight.write_step_preflight` writes `preflight.xdmf` / `.h5` (cell tags + one `facet_tags` grid: wall 1 and sheets 211–214 from the new `port_sheets` key) and `preflight_setup.png` on the congruent config. Census read back with xml / h5py / numpy equals C_solve (a second `import_step` + the solve path's sheet rebuild): counts exact on the 11 cell tags, the wall (3310 facets) and 211–214 (66 × 4), cells 116 685 (printed; primitive 116 085); worst volume / area rel diff 3.553e-15 against 1e-9 (`20260930T201605Z_GEO-34-step2-W1.log:5716–5722`); file sheet areas / `SHEET_AREA_M2` − 1 ≤ 1.11e-16 (`:5726–5729`); (n2) written without `facet_tags=`, 1 grid, 0 facet grids (`:5734`); 4 passed, 117 s. (n1, asserted) the preflight from 1a's config is unequal: 116 662 vs 116 685 cells, sheets 66 / 58 / 66 / 58 vs 66 × 4, 13 count mismatches, areas equal to 3.8e-15 (`20260930T201805Z_GEO-34-step2-W2.log:5597–5598`; 1 passed, 110 s). Stage seconds: export ≈ 28, import + tag ≈ 37, write 0.02, render ≈ 2.7. **Finding (known-issues 2026-09-30):** `import_step` returns `facet_tags = None` — its `outer_boundary_tag` group is empty on this STEP (first W1 `20260930T201313Z_GEO-34-step2-W1.log`, 2 errors). The preflight therefore tags the wall dolfinx-side as the exterior facets, and prints that it did. The test counts the wall independently (3310 = 3310); the wall area of 0.3072 m² is printed against the model bbox, rel −1.77e-6. Not ratified: the review's † readings. The fallback is **removed; the import's wall group written** (§9 item 119: exterior skin of all volumes, no tolerance; 6 surfaces, 3310 = 3310 owned exterior facets on both configs, ∫ds(1) = 0.3072 m² rel −3.3e-16; `20261001T003201Z_GEO-34-wall-W1.log`, W2 `20261001T003338Z_GEO-34-wall-W2.log:5718–5722`), so the wall re-enters step 2's identity. *(**Ruled 2026-09-30 18:00 review:** step 2 ✅ stands on the solve path's own groups — the 11 cell tags and sheets 211–214 — with (ii), (n1) and (n2). **The wall leg is a record, not part of the identity** (rule (h)): the solve path's wall group is empty, so "3310 = 3310" compares the preflight's exterior-facet fallback with the test's own exterior-facet count — a writer round trip, not C_file = C_solve. The fallback is kept as a disclosed workaround until the `import_step` wall predicate is fixed (§9 item 119, which restores the wall to the identity and removes the fallback); the † readings go to the weekly as before.)* **Step 2 Done-when (written 2026-09-30 10:30 review, rule (j); §9 item 118; the readings marked † are the review's and are named in the dashboard for ratification):** `scripts/run_step_case.py STEP CONFIG --preflight --out DIR` († "the runner" is this new entry point, preflight-only — no STEP runner exists today and a solve mode is not step 2's), in the real build and with no solve, writes the combined `.xdmf` / `.h5` (cell tags; the outer wall and the four gap sheets as facet tags, † declared by one new config key `port_sheets`) and the setup PNG for the F-small STEP under the congruent config; the census read back from that file — owned cells and volume per cell tag, facets per facet tag, area per sheet — equals that of an independent `import_step` on the same STEP and config with the solve path's sheet rebuild (`tests/validation/test_port_birdcage_four_port.py:349–352`), counts exactly and volumes and areas to 1e-9 relative; each sheet's area in the file is `SHEET_AREA_M2` within `SHEET_AREA_TOL`; and the preflight run from 1a's config fails that equality (asserted control). † The file carries the full gap sheet (66 facets), not the f = 0.5 narrowed drive sheet, which is test-level. **Why:** every fixture so far is built from Gmsh primitives inside `io/mesh.py`; the research coil, implants and anything an HFSS user would bring exist as CAD. The space of open solvers without a CAD path is full; the differentiator is a path that keeps the gates. **Step 1 — import + tag map (standard):** `io/step_import.py` (additive): `import_step(path, config)` runs `gmsh.model.occ.importShapes`, synchronises, and maps the config's named groups — CAD entity names where the STEP carries them, else geometric selectors (bounding box / centroid / axis, the same predicates the generator already uses) — to the project's physical-group tags (volumes: conductor / phantom / air / implant; surfaces: PEC, impedance, port sheets with drive direction and reference impedance); sizing fields per group; returns the same `(mesh, cell_tags, facet_tags, diagnostics)` tuple the primitive generators return. Config is one YAML/JSON file next to the STEP — the same shape as a Palace config and what an HFSS user thinks in. **Gate (asserted):** the F-small four-leg birdcage rebuilt from a STEP that a new CadQuery/build123d emitter writes from the *same* parameters (`geo/step_emit.py`, additive) reproduces, through `import_step`, the `GEO-18` CAD identities (meshed/CAD conductor mass ≥ 0.95, terminal-area ratios in the imported band) and `PORT-9`'s four-port gates at 10 MHz (reciprocity / passivity / C4 spread, imported bands) — **the same identities, not the same cell count**: a STEP round-trip changes the CAD kernel's tessellation, so cell-count records are *printed beside the primitive rung's and not asserted* (the (1*) licence question is named in the row: mesher-dependent records become routine under STEP, and a re-record needs the reconciler, never this step). **Negative control (asserted):** the config with the port-sheet group deliberately unmapped fails the import with a message naming the missing group — a silent untagged port is the failure a config invites. **Step 2 — preflight export (standard):** `--preflight` on the runner: mesh, tag, and write the combined XDMF plus the setup figure with materials, boundaries and ports as coloured regions and *no solve*; gate: the preflight's tag census equals the solve path's on the same STEP. This is the "see what I am doing before I solve" the operator asked for; an interactive viewer (pyvista/trame) is a later `EX` item, not a product surface. **Step 3 — an external STEP (heavy by prediction, cost probe first):** a STEP not written by our emitter — a 16-leg high-pass birdcage drawn in FreeCAD or CadQuery by the operator with named groups — imported, meshed, port-gated at 10 MHz; the first design that did not come from `io/mesh.py`. **Not in scope:** editing, a GUI, IGES/STL (STL has no CAD entities to name; a later step if needed). **Sequencing:** independent of the physics chains; step 1 is a slot item once the review costs it (the emitter and importer are pure Gmsh/OCC — measured in seconds on F-small); `GEO-33`'s strip coil is the first consumer. **Costed and split by the 2026-09-23 03:00 review (§9 item 37):** step **1a** = the STEP round trip, the config's group map and the *mesh-level* gates only (tag census, `GEO-18`'s analytic sheet area, the CAD-mass gate; the unmapped-group control) — no solve; step **1b** = `PORT-9`'s four-port gates at 10 MHz on the imported mesh; step **1c** = the independent-kernel emitter. **Emitter substitution in 1a:** CadQuery / build123d are not in the image (`docker/Dockerfile`), and adding them is an operator rebuild; 1a emits the STEP with Gmsh's own OCC writer (`gmsh.write("….step")` on the generator's model before meshing), which exercises the STEP reader, entity renumbering and the group map but not an independent CAD kernel — that independence is 1c's, when the image carries one (the weekly / operator decides). 

### §9 item 289 — 2026-10-06 (chain U12)

**Outcome 2026-10-07 03:30 slot — 🚫 BLOCKED at W0, branches (b) and (c) together** (`20261007T034333Z_GEO-34-step1d-W0-names.log`, Status 0, 1 s, one rank, gmsh 4.15.2, probe `tests/mesh/probe_geo34_step_names.py`). The writer reported 396 STEP entities, 18 597 bytes; the strings `PortSheet` and `AirBox` occur **0** times in the file text, although both were set by `setEntityName` and as named physical groups. On re-import the one volume is named with the translator's label (`Shapes/Open CASCADE STEP translator 7.6 1/…`), every dim-2 name is empty, and there are no physical groups (`:47–60`). **The embedded 0.48 m² rectangle is gone too:** after import there are only 6 faces, each of area 1.0 (the box skin), `getEmbedded` is empty, and nothing changes after `removeAllDuplicates` (`:61–74`). So Gmsh's own writer cannot author this fixture on this image, whether by name or by geometry. No `src/` edit and no test were written. There is no geometric-selector substitute (brief rule). Unblock: a STEP from a kernel that writes named, embedded faces, for example `cq-python` (dashboard item 10, step 1c). After that the weekly re-poses U12. Known-issues entry 2026-10-07.

*Drafted by a general-purpose agent for the 2026-10-06 18:00 daily review; the review reviews it.* Licence: §10 U12 (grep "U12. `GEO-34` step 1d"), quoted, not changed: "A STEP from Gmsh's own OCC kernel (the operator's independent-kernel STEP is step 1c, dashboard item 10, and stays open): one box, one named interior rectangle; the config maps the name to `port`; the importer tags it natively (U11) and the meshed area equals the CAD rectangle at `STRIP_BAND` at `-n 1 / 2 / 4`. Target: CAD identity. Serial on U11." U11 landed green as item 267 (`f0144df2`): the native interior tag works **on a `shared_facet`-ghosted mesh**, and the unghosted `-n 4` leg hangs. **No band, tolerance, Done-when or validation target is set or moved.**

**`GEO-34` step 1d: a named interior STEP face becomes a port tag through `import_step`.**
- **Executor:** implementer. It makes an additive `src/` change (`io/step_import.py`) and adds tests. Spawn prompt: "Use the implementer agent to execute `GEO-34` step 1d (§9 item 289) per `docs/planning/chunks/GEO-34.md` § '§9 item 289'; W0 first, before any code."
- **Tier:** standard by ceiling (180 s), smoke-sized by prediction. The fixture is item 222's: 10 150 tets, import 0.14 s, item 267's W2 56 s for seven launches. `-n 1 / 2 / 4`, real build (every integrand is real).
- **Class:** chain step (solve). **Time:** ≈ 9 measured / ≈ 45 predicted min (15 fixed + ≈ 20 code + W0 ≈ 1 + W1 ≈ 1.5 + W2 ≈ 2 + ≈ 5 spare).
- **Depends on:** U11 (item 267), landed. Nothing open.

**The importer today** (`src/fem_em_solver/io/step_import.py`, read 2026-10-06).
- `import_step` claims **volumes only**, by geometric selectors (`claim_volumes`, `:108`).
- The only surface groups it writes are the optional `outer_boundary_tag` (`:196–240`) and the `congruent_sheets` binding (`:328`). It reads **no STEP entity names**.
- It already distributes through `_model_to_mesh` with `GhostMode.shared_facet` (`:411–417`), so U11's ghost requirement holds by construction. The step must not add a path that meshes unghosted.
- The preflight's `port_sheets` key (`io/step_preflight.py:12`) rebuilds sheets dolfinx-side from volume pairs (`between`). It is not a native tag and is not touched here.

**W0 — name survival (decides the item; run before any code).** It is unknown whether Gmsh 4.15.2's OCC STEP writer carries a face's name, and whether `importShapes` surfaces it through `gmsh.model.getEntityName`. Nothing in the repo has measured this; `scripts/maintenance/install_cadquery_image.sh:111` prints names only for a CadQuery file.
- Probe `tests/mesh/probe_geo34_step_names.py` (not collected; committed as the record). On one rank:
  1. Item 222's box plus rectangle: import `RECT` / `RECT_AREA` / `H` from `tests/mesh/probe_geo37_interior_facet.py`; `occ.fragment([(3, box)], [(2, rect)])`.
  2. Name the rectangle `PortSheet` **both** ways: `gmsh.model.setEntityName(2, s, "PortSheet")` and a dim-2 physical group `setPhysicalName(2, 201, "PortSheet")`. Name the volume `AirBox` the same two ways.
  3. `gmsh.write(<tempfile.mkdtemp(prefix="geo34_1d_")>/box_port.step)`, then finalize.
  4. Re-open with a fresh session and `occ.importShapes`. Print `getEntityName` for every dim-2 / dim-3 entity, the face count, and each face's bbox and `getMass` area, both **before and after** `occ.removeAllDuplicates()`. Also print whether the 0.48 m² face is in the volume's `getBoundary` (embedded) or free.
- Command:

  `scripts/testing/run_and_log.sh GEO-34-step1d-W0-names "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 10 30 python3 tests/mesh/probe_geo34_step_names.py'"`
- **Branches:**
  - **(a)** A face named `…PortSheet` (exact, or as the last `/`-separated component, since OCC labels can carry an assembly path) exists after import, with area 0.48 to 1e-12. Proceed. Record which naming route survived and whether the name also survives `removeAllDuplicates`.
  - **(b)** No name survives. **Stop:** commit the probe and log as a record, journal the attempt, mark BLOCKED. The finding is that Gmsh's writer cannot author a named STEP on this image. A geometric selector is not "a named face", so the review re-poses U12 (for example onto `cq-python` once dashboard item 10 lands). Do not substitute a bbox selector.
  - **(c)** The name survives but the face is not re-embedded after `removeAllDuplicates` (free face, or two faces). Stop as in (b), with that finding.

**Code (on W0 branch (a) only).**
1. **`io/step_import.py`, additive config key `surface_groups`:** `[{"name": "port_1", "step_name": "PortSheet", "role": "port", "tag": 201}, …]`. If the key is absent, the import is byte-identical (rule (c), W2).
   - **Resolve names before healing.** Right after `importShapes` + `synchronize` and **before** `removeAllDuplicates`, map each `step_name` to its dim-2 entities by `getEntityName` (exact, or last path component), and record each face's bbox and `occ.getMass(2, ·)` area. Entity tags renumber under `removeAllDuplicates`, so a pre-heal tag is never used afterwards.
   - **After the heal,** re-find each face among the post-heal dim-2 entities: bbox equal per coordinate to 1e-6 (the OCC padding trap, item 222 W1) and area equal to 1e-12 relative. Require exactly one match.
   - **Errors** are `ValueError`, raised on the building rank and re-raised on every rank by the existing broadcast. They cover:
     - a `step_name` absent from the file (the message names it and lists the names present);
     - zero or several post-heal matches;
     - a `role` other than `"port"` (the only role this step licenses);
     - a tag colliding with a volume tag or `outer_boundary_tag`.
   - Then `addPhysicalGroup(2, [s], tag)` and `setPhysicalName(2, tag, name)`.
   - **Diagnostics** gain `surface_groups: {name: {"tag", "role", "step_name", "surface", "cad_area", "n_triangles"}}`. `n_triangles` is read after `generate(3)` from `gmsh.model.mesh.getElements(2, s)`.
   - The module docstring's config block gains the key.
2. **Config** `tests/mesh/data/geo34_step1d_port_sheet.json`:
   - `resolution` = item 222's `H` (0.08; the test asserts equality, as `test_config_sizes_are_the_gates` does);
   - one volume group `air` (tag 1, `bbox_contains` the unit box's interior point box);
   - `outer_boundary_tag` 101;
   - `surface_groups` `[{"name": "port_1", "step_name": "PortSheet", "role": "port", "tag": 201}]`.
3. **Test** `tests/mesh/test_step_import_port_sheet.py`, real build, any width.
   - A module fixture writes the STEP on rank 0 into a `tempfile.mkdtemp(prefix="geo34_1d_")` directory, broadcast as `test_step_import.py:113–117` does. The writer is W0's naming route, factored into a function and imported from the probe. The STEP is never written into the tree.
   - Then `import_step(step, CONFIG, comm)`.

**Anchors (asserted at every width; `STRIP_BAND` = 1e-9 imported from `tests/validation/test_ops61_exact_strip_ports.py:58`).**
- **(P0) ghosts:** `comm.allreduce(index_map(3).num_ghosts)` is > 0 at `-n 2 / 4` and 0 at `-n 1`. The width's `cells_ghost` is printed. This is U11's precondition, checked on the imported mesh.
- **(P1) CAD identity:** the owned tag-201 facet area (numpy triangle areas over facets with index < `index_map(2).size_local`, `allreduce`d) satisfies \|A/`RECT_AREA` − 1\| ≤ `STRIP_BAND` and \|A/diagnostics `cad_area` − 1\| ≤ `STRIP_BAND`.
- **(P2) completeness:** the reduced owned tag-201 facet count equals diagnostics `n_triangles` exactly. Every tag-201 facet has two incident cells on the ghosted mesh: the one-sided count, reduced, is 0. That proves the face is embedded, not a skin.
- **(P3) the port form (U11):** `assemble_scalar(form(1 * dS(201)))`, allreduced, with `metadata={"quadrature_degree": 2}`, satisfies \|I/`RECT_AREA` − 1\| ≤ `STRIP_BAND`.
- **(P4) the role reaches the port layer:** the globally reduced facet-tag set is {101, 201}, and `validate_required_port_tags_exist` (`ports/definitions.py:82`) passes on 201 with that set.

**Negative controls (rule (e)).**
- **(n1) asserted, unknown name:** the same STEP with `step_name` `"PortSheetX"` raises `ValueError` on every rank, and the message contains `PortSheetX`. This is the "silent untagged port" failure the config invites (1a's unmapped-group control, for faces).
- **(n2) asserted, the mapping makes the tag:** the same STEP with `surface_groups` removed imports, and its reduced tag set lacks 201. So the tag comes from the config, not from any physical group or label the STEP carries (W0 records which ones survive).
- **(n3) asserted, the comparator is not blind:** re-tag 201 without its first owned facet on the lowest owning rank (item 267's (c2)). \|A′/`RECT_AREA` − 1\| > `STRIP_BAND`. The ceiling is ≈ 1/`n_triangles` (item 267: −4.5e-3 / −5.1e-3), ≈ 10⁶ × the band.

**Predictions** (printed): (P1) / (P3) residuals ≤ 1e-14 (item 267: ≤ 1.6e-15); `n_triangles` near 192 (item 222's count on the primitive). The count is printed, not asserted against 192, because a STEP round trip may move the tessellation (1a's licence). There is no `PREDICTED RED:` line.

**Windows** (each bare, whole and foreground; Bash `timeout: 660000`; `-s`; pytest is the last command of each launch).
- **W0:** above (≈ 5–15 s).
- **W1** (closing). Widths are chained with `&&` so the footer carries the first failure. **`tests/environment` is narrowed to `test_complex_mode.py`**: item 267's W1 lost its whole 30 s pytest ceiling to `test_dolfinx_api_migration.py`'s repository scan.

  `scripts/testing/run_and_log.sh GEO-34-step1d-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 10 60 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_step_import_port_sheet.py -v -s --tb=short && PYTHONPATH=/workspace/src:/workspace timeout -k 10 60 mpiexec -n 2 python3 -m pytest tests/mesh/test_step_import_port_sheet.py -v -s --tb=short && PYTHONPATH=/workspace/src:/workspace timeout -k 10 60 mpiexec -n 4 python3 -m pytest tests/mesh/test_step_import_port_sheet.py -v -s --tb=short'"`

  Predicted 25–60 s.
- **W2, rule (c): the key-absent path is unchanged.** It runs 1a's module at `-n 2`. The precedent is the 1b′ rule-(c) window, 85 s, 4 passed (`20260926T094036Z_GEO-34-step1bprime-ruleC.log`).

  `scripts/testing/run_and_log.sh GEO-34-step1d-W2-ruleC "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_step_import.py -v -s --tb=short'"`
- No edit after W1 starts (rule (i)). If W2 must follow an edit, W1 is re-run.

**Traps already paid for.**
- **Unghosted `dS` hangs at `-n 4`** (item 267 (c3); 2026-08 `PORT-1` 3b-iv). Mesh only through `import_step`, which ghosts. Never call `_model_to_mesh` with the default partitioner in this test.
- **OCC bbox padding ≈ 1e-7:** match at 1e-6 and guard with the exact area (item 222 W1).
- **Pre-heal tags die in `removeAllDuplicates`:** resolve names before the heal, then re-find faces by geometry.
- **Rank-local reads:** `allreduce` every count, area and `assemble_scalar`. Gather the tag set before the validator (`OPS-14`). Use the owned filter `index_map(2).size_local`.
- **A rank-0 gmsh throw must reach every rank before the collective:** reuse `import_step`'s broadcast, never a rank-0-only raise (`GEO-23`).
- **A `SpatialCoordinate`-free measure still needs a pinned quadrature degree** on gmsh meshes (`POST-5`).
- **One gmsh session at a time in a process:** finalize the writer's session before `import_step` initializes its own.
- **After a kill:** run `pgrep -c python3` first, then clear `/dev/shm/mpich_shm_*` (known-issues). Run the orphan-`.c` sweep only at `pgrep` 0.
- No pipes. `PYTHONPATH` carries `:/workspace`. Write ≈, never `~`.

**Scope.** Covered: one planar rectangle in one box from Gmsh's own OCC writer; a name-to-`port` mapping; the native tag on the ghosted import at `-n 1 / 2 / 4`; area, completeness and `dS` against the CAD rectangle. **Not claimed:**
- an independent-kernel STEP (1c);
- any port solve, S or F14 integration line;
- roles other than `port` (PEC / copper / lumped sheets);
- curved, multi-face or interface sheets;
- the F-small STEP (U13);
- the preflight's facet grid for `surface_groups`;
- `xl` widths.

**Status it can move.**
- **Green:** the `GEO-34` §7 row gains a "step 1d" record clause (the name route that survived, P1–P3 residuals per width, elapsed). §10 U12 is marked DONE and U13 (the `Study` on a STEP-imported F-small) is unblocked.
- **No glyph:** steps 1c and 3 stay open (row read 2026-10-06), so `GEO-34` stays 🟡.
- **Rule (j):** the row has no step-1d Done-when. The review adds this one at queueing: "a Gmsh-OCC STEP of one box with one named interior rectangle, mapped by `surface_groups` to `port`, imports with the rectangle as a native facet tag whose owned area equals the CAD rectangle at `STRIP_BAND`, count equal to gmsh's triangles, `dS` at `STRIP_BAND`, at `-n 1 / 2 / 4` in one logged window; unknown-name control raises".
- **W0 (b) / (c):** a known-issues entry ("Gmsh 4.15.2's STEP writer drops entity names", or "does not keep an embedded face"), and U12 goes back to the weekly for re-posing.

**Negative result.** Report and stop. Never substitute a geometric selector, and never widen a band.
- **(P1) / (P3) red with (P2) and (n3) green:** the imported tag is wrong on 0.11. Open a known-issues entry and park on `attempt/GEO-34-step1d-<UTC>`, BLOCKED.
- **(P0) zero ghosts at `-n 2 / 4`, or a 60 s kill at `-n 4`:** this contradicts U11 on an imported mesh. Record it and mark BLOCKED.
- **(n1) / (n2) / (n3) red:** the instrument is wrong and the verdict void. Fix it once in the new test, re-run in-slot and disclose. If it is still red, park and mark BLOCKED (rule (d)).
- **W2 red:** the key-absent path moved. Revert the `src/` edit and mark BLOCKED.

### §9 item 333 — 2026-10-08 (chain U12a — GEO-34 step 1d probe; mesh-probe)

*Drafted by a general-purpose agent for the 2026-10-08 10:30 daily review; the review reviews it.* Licence: §10 U12a (grep "U12a. `GEO-34` step 1d probe"), quoted, not changed: "`mesh-probe`: the box-and-rectangle model exported as a compound (solid + the face as its own shell) and as a shell-only STEP; after `importShapes`, the number of dim-2 bodies and their areas against 0.48 m², names present or not, at `-n 1`. Measurement only; green ⇒ U12′ gains a same-kernel rehearsal branch, red ⇒ 9b waits on item 10. Serial on nothing." **No band, tolerance, Done-when or validation target is set or moved.** Measurement only (§3): nothing is asserted, and `GEO-34` stays 🟡.

**`GEO-34` step 1d probe: can Gmsh's own OCC kernel carry an interior face through a STEP file in any export form?**
- **Executor:** `mesh-probe`. It produces a table and one probe script. There is no `src/` change, no test, and no assertion. Spawn prompt: "Use the mesh-probe agent. Generator: Gmsh 4.15.2's OCC STEP writer and `occ.importShapes` on `GEO-37` item 222's box-and-rectangle fixture. Swept parameter: the export form (F-frag / F-solid / F-comp / F-shell / F-two) × the `importShapes` read (default / `highestDimOnly=False`). Question: does the 0.48 m² rectangle come back as exactly one face embedded in the one volume after `removeAllDuplicates`? Brief: `docs/planning/chunks/GEO-34.md` § '§9 item 333'. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Pass `timeout: 660000` and run in the foreground."
- **Tier:** smoke (30 s), by prior measurement. W0 wrote and re-read this exact fixture in one process in **1 s** elapsed (`20261007T034333Z_GEO-34-step1d-W0-names.log:79`). This item runs five child processes, so predict ≈ 5–10 s per window. Width `-n 1`, plain `python3` with no `mpiexec`: a gmsh throw on rank 0 at `-n 2` deadlocks the command (mesh-probe rule). Real build; no DolfinX object is built.
- **Class:** chain step (no solve). **Time:** ≈ 7 measured / ≈ 30 predicted min (15 fixed + ≈ 10 script + 2 windows ≈ 1 + ≈ 4 journal / table).
- **Depends on:** nothing.
- **Family:** the 03:00 review counted the `GEO-34` 1d family at 1 of 4 (item 289's W0 BLOCKED). The review rules whether a measurement-only probe counts. If it does, this item makes it 2 of 4 and leaves U12′ at 3.

**What exists today.**
- **The importer:** `src/fem_em_solver/io/step_import.py:154–166` (`_build`). It runs `gmsh.model.occ.importShapes(str(path))` **with the default arguments** (`:164`), then `occ.removeAllDuplicates()` (`:165`), then `synchronize`. It sets no `Geometry.OCC*` option except the optional `OCCTargetUnit` (`:159–161`).
- **The fixture:** `RECT = (0.1, 0.2, 0.47, 0.8, 0.6)`, `RECT_AREA = 0.48`, in `tests/mesh/probe_geo37_interior_facet.py:30–31`.
- **W0's writer:** `write_named_box_port_step` in `tests/mesh/probe_geo34_step_names.py:28–64`. It fragments the box with the rectangle, names both entities two ways, writes, and finalizes.
- **What W0 already paid for:** neither name string occurs in the file (`:45–46`). On re-import there is 1 volume and 6 faces, each of area 1.0. No 0.48 face exists, nothing is embedded, and `removeAllDuplicates` changes nothing (`:47–74`). The OCC bbox carries ±1e-7 padding (`:49–56`).

**Probe:** `tests/mesh/probe_geo34_step1d_export_forms.py` (new; not collected; committed as the record). It imports `RECT` / `RECT_AREA` from `probe_geo37_interior_facet` and `write_named_box_port_step`, `PORT_NAME` and `VOL_NAME` from `probe_geo34_step_names`. It restates none of them.
- **Fresh process per form.** The driver runs `subprocess.run([sys.executable, __file__, "--form", F])` once per form. Each child writes its STEP in one gmsh session, finalizes it, and reads the file back in a second session that it also finalizes.
- **Exit code.** The driver prints each child's stdout. A crashed child becomes a `FAIL` cell holding the last stderr lines. The driver exits with the largest child return code, so the footer is honest.
- **Scratch files.** STEP files go to `/workspace/logs/geo34_step1d_probe/<form>/`, which is the gitignored `/logs/`. They are never written into the tracked tree.
- **The forms.** Wherever a form has them, the face is named `PortSheet` and the volume `AirBox`, both by `setEntityName` and by a named physical group, exactly as in W0.
  - **F-frag (control, W0 reproduction):** `write_named_box_port_step`, unchanged.
  - **F-solid (control):** the unfragmented unit box only.
  - **F-comp (the compound):** the unfragmented box plus `addRectangle(*RECT)` as a free face, in one model with no `fragment`, and one `gmsh.write`.
  - **F-shell (shell-only):** the rectangle alone.
  - **F-two (two-file merge):** F-solid's file and F-shell's file are both passed to `importShapes` in one reader session.
- **The two reads.** Each file is read (r1) with the default `importShapes(path)`, which is how the importer reads at `:164`, and (r2) with `importShapes(path, highestDimOnly=False)`.
- **What each child prints:**
  - first, `inspect.signature(gmsh.model.occ.importShapes)`;
  - the value of every `Geometry.OCC*` option the reader can read: `OCCImportLabels`, `OCCSewFaces`, `OCCMakeSolids`, `OCCFixDegenerated`, `OCCFixSmallEdges`, `OCCFixSmallFaces`. These are read, never set;
  - the file's size and the number of STEP entities the writer reports;
  - how often `PortSheet` and `AirBox` occur in the file text;
  - the rest **before and after** `removeAllDuplicates()`: the dim-3 and dim-2 entity counts, and for every face its `getEntityName`, `occ.getMass`, bbox and upward adjacencies;
  - the counts of faces with \|A/0.48 − 1\| ≤ 1e-12 and bbox equal to RECT's per coordinate to 1e-6, of **free** faces (no upward adjacency), and of **embedded** faces (upward adjacency equal to the one volume);
  - `mesh.getEmbedded(3, v)` and the physical groups.
- **The table.** The last line block is one row per form × read × stage, in the mesh-probe format.

**Windows** (each bare, whole and foreground; Bash `timeout: 660000`).
- **W1:** `scripts/testing/run_and_log.sh GEO-34-step1d-probe "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 30 python3 tests/mesh/probe_geo34_step1d_export_forms.py'"`
- **W2 (repeat, the mesh-probe determinism rule):** the same command under the ID `GEO-34-step1d-probe-repeat`. Only a bit-identical table licenses the word "deterministic".

**The anchor: fixed comparands, read as brackets in the readout, never asserted.**
- A_rect = 0.8 × 0.6 = **0.48 m²** (`RECT_AREA`), matched at \|A/0.48 − 1\| ≤ 1e-12, with the bbox matched at 1e-6 (the padding trap). Box faces are 1.0 m² and the volume is 1.0 m³.
- **Expected counts per form** (dim-3 / dim-2 / 0.48-faces / free / embedded, after the heal). Rule (e) labels:

| Form | Read | Expected | Label |
|---|---|---|---|
| F-frag | r1 | 1 / 6 / 0 / 0 / 0 | **backed**: W0 `:47–74`, same writer, same fixture |
| F-frag | r2 | 1 / 6 / 0 / 0 / 0 | predicted (the face is not in the file, by W0) |
| F-solid | r1, r2 | 1 / 6 / 0 / 0 / 0 | predicted |
| F-comp | r1 | 1 / 6 / 0 / 0 / 0, the free face dropped by `highestDimOnly=True` | predicted, from the Gmsh API's documented default; unverified on this image |
| F-comp | r2 | before the heal 1 / 7 / 1 / 1 / 0; after it 1 / 7 / 1 / 0 / 1 | predicted |
| F-shell | r1, r2 | 0 / 1 / 1 / 1 / 0 | predicted |
| F-two | r2 | after the heal 1 / 7 / 1 / 0 / 1 | predicted |
| all | — | `PortSheet` / `AirBox` occurrences 0 / 0 | backed for F-frag (W0 `:45–46`), predicted for the rest |

**Decision rule (§10, unchanged).**
- **Green** means that F-comp (under either read) or F-two leaves **exactly one** 0.48 face, **embedded** in the one volume, after `removeAllDuplicates`. Then U12′ gains a same-kernel rehearsal branch. The readout names the form and the read. If the face survives only under r2, the readout says so: the rehearsal (and any CadQuery-assembly STEP) then needs an additive `import_step` read option, which is U12′'s change, not this item's.
- **Red** means no form returns an embedded 0.48 face. Then 9b waits on dashboard item 10.
- **Amber** means the face survives free but is not embedded after the heal (W0's branch-(c) shape). It is read as red for the decision; the review rules.
- **Names do not enter green.** They are recorded. U12′'s name-to-`port` map still needs a kernel that writes labels, and a geometric-selector substitute stays forbidden.

**Negative control.** F-frag reproduces W0's 6-face, no-0.48 reading. F-solid shows zero free dim-2 bodies and zero 0.48 faces. If F-frag departs from W0 `:47–74`, the instrument or the image has moved: the window is void, and the executor stops and reports without taking a decision. If F-solid shows a 0.48 face or a free face, the face-matching code is wrong; the same applies.

**Traps already paid for.**
- **`highestDimOnly` (new, unmeasured).** The importer's default read may silently drop a free face that sits beside a solid. That is why every file is read both ways and the signature is printed. Do not "fix" `step_import.py:164`.
- **One gmsh session per process.** Finalize the writer before the reader initializes (item 289 trap). Use one child process per form.
- **Healing:** `removeAllDuplicates` renumbers tags. Identify the face by area plus bbox, never by a pre-heal tag (item 289). Print the `OCCSewFaces` / `OCCMakeSolids` / fix options and never toggle them, because a sew could merge the free face into a shell.
- **Names:** W0 found that the *writer* emits no label string at all (`:45–46`), so no reader option (`OCCImportLabels`) can recover a name. No writer-side label option is known to the repo. Toggling one is a follow-up for the review, not this item.
- **OCC bbox padding ≈ 1e-7:** match at 1e-6 and use the exact area (item 222 W1; W0 `:49–56`).
- **Harness and shell:**
  - The harness call is the whole Bash command (`bash_guard` rule 3). Container-side `timeout -k 30`, ceiling 30 s.
  - No pipe, no `tee`, no redirect. `PYTHONPATH` carries `:/workspace`.
  - STEP scratch goes under `/workspace/logs/`, never `/tmp` in the tree, and never as a committed artifact. Never delete a log.
- **After a kill:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` comes first. A kill at 30 s is not a result. Shrink to one `--form` per window, never a longer timeout.
- Write ≈, never `~`.

**Scope.** Covered: Gmsh 4.15.2's OCC writer and reader on one unit box and one planar rectangle; five export forms × two reads; counts, areas, adjacency and name strings at `-n 1`. **Not claimed:**
- any gate, Done-when or glyph (`GEO-34` stays 🟡; steps 1c, 3 and U12′ stay open);
- a name route;
- a mesh, a facet tag, `dS`, or any multi-rank reading;
- a CadQuery STEP;
- the F-small STEP (U13).

U12′ stays serial on the operator's `cq-python` install. This item cannot make U12′ runnable on its own: a green opens only a *rehearsal* branch beside it, which the review or the weekly would write.

**Negative-result disposition.** Report the table and stop. The probe and both logs are committed as the record. Nothing in `src/` is touched, and no option is changed to make a face appear. A crashed child is a `FAIL` cell holding its message. It is not retried in-slot.

**Status it can move.**
- **§10 U12′:** a green adds the same-kernel rehearsal branch, naming the form and read. A red adds "9b waits on dashboard item 10, measured: no Gmsh-OCC export form carries the face".
- **Dashboard item 10:** its wording gains the measured clause ("Gmsh's own writer keeps no face names **and no interior face in any export form**", or "…keeps no names; a compound export carries the face, so a same-kernel rehearsal runs meanwhile").
- **Known-issues 2026-10-07** (the "Gmsh 4.15.2's OCC STEP writer drops entity names…" entry): a *Scope* addendum on which forms keep the face. It does **not** retire, because names are still dropped.
- **`highestDimOnly`:** if r1 drops a free face that r2 keeps, the readout flags it to the review as a finding against U12′'s importer read and against `scripts/maintenance/install_cadquery_image.sh:111–120` (step 4b's `len(faces) >= 7` read with the default). The review decides whether that is a known-issues entry.
- **The `GEO-34` §7 row** gains a one-clause U12a record (🧪 measurement, glyph unchanged). This chunk file gains the outcome paragraph.

**Outcome (2026-10-08 16:30 slot): AMBER, read as red. A measurement record; `GEO-34` stays 🟡; 1d family 2 of 4.** Executor: mesh-probe. The probe is `tests/mesh/probe_geo34_step1d_export_forms.py`. W1 `20261008T214725Z_GEO-34-step1d-probe.log` and W2 `20261008T214738Z_GEO-34-step1d-probe-repeat.log` each ran Status 0 in 5 s at `-n 1`, with tables bit-identical (deterministic). Each cell below reads dim-3 / dim-2 / 0.48 / free / embedded. The heal changed nothing, so pre and post are equal everywhere.
- **Controls valid.** F-frag reproduces W0 (1 / 6 / 0 / 0 / 0 under both reads), and F-solid is clean (1 / 6 / 0 / 0 / 0).
- **Compound (F-comp).** r1 reads 1 / 6 / 0 / 0 / 0: the free face is dropped and `importShapes` returns `[(3, 1)]`. r2 reads 1 / 7 / 1 / **1 / 0**.
- **Shell only and two-file merge.** F-shell reads 0 / 1 / 1 / 1 / 0 under both reads. F-two, which goes beyond §10's two forms, reads 1 / 7 / 1 / 1 / 0 under **both** reads, because `highestDimOnly` applies per file.
- **Departure from the prediction.** The brief predicted that `removeAllDuplicates` would embed the free face (1 / 7 / 1 / 0 / 1). It does not on this image: `getEmbedded(3, v)` is `[]`, upward adjacency is `[]`, and `OCCSewFaces` / `OCCMakeSolids` are 0, read and never set.
- **Names and options.** `PortSheet` / `AirBox` occur 0 times in every file. `OCCImportLabels` is 1; the other OCC* options are 0.
- **Decision.** The result is W0's branch-(c) shape in every form, so 9b waits on dashboard item 10.
- **`highestDimOnly` finding.** r1 drops F-comp's free face. Known-issues entry 2026-10-08 is opened per ruling (Q4), naming `step_import.py:164` and `install_cadquery_image.sh:111–120`. The dashboard line is the next review's.

**Review rulings (2026-10-08 10:30 daily review, item 333).** (Q1) ≈ 7 measured: the chain no-solve median recalibrated this review to 7 (316, 327, 308, 313, 298). (Q2) **The probe counts as an attempt in the `GEO-34` 1d family** (2 of 4), because it is a lettered sub-step of the family that lands no gate; U12′ would be 3 of 4. The cap stays honest either way. (Q3) Yes: the rehearsal branch, if green, rehearses embed → native tag → `dS` only, and names stay U12′'s gate. The weekly decides whether that branch is worth a slot. (Q4) No warning before measurement. If read r1 drops the free face as predicted, the slot opens a known-issues entry naming `step_import.py:164` and `install_cadquery_image.sh:111–120`, and the next review adds a line to dashboard item 10. (Q5) Keep F-two (≈ 1 s), printed and labelled "beyond §10's two forms". (Q6) No writer-option sweep in this item.

### §9 item 985 — 2026-11-14 (§10 U12′, chain U — rule (4), the weekly's enumerated step; GEO-34 step 1d re-posed; chain step (no solve))

*Written by the 2026-11-14 10:30 daily review under §9 rule (4); drafted by a `general-purpose` agent and reviewed here. Corrections made at review:*
- *the cost uses the current class median;*
- *no Done-when is added to the §7 row (a daily never sets one);*
- *the module hold is checked.*

*It does not count against rule (6). The target is the weekly's U12′ (`docs/planning/plan-archive.md:39275–39281`), quoted, not changed: "a CadQuery-authored STEP with a named interior face … the config maps the name to `port`, the importer tags it natively (U11, `shared_facet` ghosts) and the meshed area equals the CAD rectangle at `STRIP_BAND` at `-n 1 / 2 / 4`. Target: CAD identity." **The serial prerequisite has landed:** the operator's CadQuery install, 2026-11-14 (`fa4a6e7f`). **Family `GEO-34` 1d: attempt 3 of 4.** Item 289 was attempt 1 (BLOCKED at W0), and U12a / item 333 was attempt 2 (ruling Q2 above). A red here leaves one attempt.*

**Facts verified (by reading).**
- **Install log** (operator directive (b); gitignored `logs/maintenance/cadquery_image_20261114T091828Z.log`):
  - step 4a is green: `cadquery 2.8.0`, then `wrote …cq_sample_…step` (`:381–382`);
  - step 4b **FAILED** on `-noenv` (`:434`, `:478`), and the log stops there;
  - the corrected read (1 volume, 7 faces, both names) is in `fa4a6e7f`'s message only. **No log records it, so W0b re-measures it in a logged window.**
- **The 4d log does not test interior facets.** `20261114T093413Z_GEO-34-cadquery-image-4d.log` reads 5 passed, 13 skipped, 104 s, Status 0. The passes are all `test_step_import.py`'s. Every `test_geo37_interior_facet.py` node skipped because `GEO37_RECORD_DIR` is unset. Cite this log for the image, not for interior facets.
- **The sample's sheet is a free face, not an embedded one.** It is its own product, `PRODUCT('PortSheet')` → `OPEN_SHELL`, beside the `AirBox` solid. U12a measured that `removeAllDuplicates` does not embed a free face. So the importer must embed it explicitly, as item 222 did with `occ.fragment([(3, box)], [(2, rect)])` (`probe_geo37_interior_facet.py:40`).
- **Units.** The STEP declares millimetres but holds the numbers as drawn. Leave `occ_target_unit` unset.
- **The importer** (`src/fem_em_solver/io/step_import.py`):
  - `_build` calls `importShapes(str(path))` with defaults (`:164`), then `removeAllDuplicates` (`:165`);
  - it claims volumes only, and the wall is the combined boundary of every volume;
  - it ghosts with `shared_facet`;
  - unknown config keys are ignored (`step_preflight.py:14`), so additive keys have a seam.
- **Fixture.** `RECT = (0.1, 0.2, 0.47, 0.8, 0.6)`, `RECT_AREA = 0.48`, `H = 0.08` (`probe_geo37_interior_facet.py:30–32`). `STRIP_BAND = 1.0e-9` (`tests/validation/test_ops61_exact_strip_ports.py:58`).
- **U11.** The native-tag `dS` passed at `STRIP_BAND` only on a ghosted mesh. Unghosted `-n 4` was killed at 45 s (`GEO-37.md:165–176`).
- **Ignore rules.** `*.step` is gitignored (`.gitignore:57–58`).
- **Hold (checked at review).** No queued window (64, 65) imports `io/step_import.py`; `io/mesh.py` names it only in a docstring (`:4385`). The importer is free.

**Code** (scripts, tests, one committed STEP, two additive importer keys).
1. **`tests/mesh/cq_author_geo34_step1d.py`** (new; not collected). It runs under `cq-python`, which clears `PYTHONPATH` (`Dockerfile:41`), so it restates `RECT` as literals.
   - Build the box with `cq.Workplane("XY").box(1, 1, 1, centered=False)`.
   - Build the sheet as an explicit polygon with corners (0.1, 0.2, 0.47) … (0.1, 0.8, 0.47). Do not use `makePlane`: its orientation is ambiguous.
   - Assemble the products `AirBox` and `PortSheet`, as the sample did.
   - Assert, CadQuery-side, face area 0.48 to 1e-12 and the bbox.
   - Export with `asm.export(out)`; 2.8 warns that `save` is going. If `export` is absent, fall back to `save` and disclose it.
   - Print the `sha256` and the byte count.
2. **Commit the STEP** as `tests/mesh/data/geo34_step1d_port_sheet.step`, with a narrow `!tests/mesh/data/geo34_step1d_*.step` negation after `.gitignore:58`. Never regenerate it per run:
   - pytest runs in the DolfinX interpreter, which never imports CadQuery;
   - only `fem-em-solver` carries the venv, so other services would skip, and a skip hides a failure;
   - the STEP header carries a timestamp;
   - an externally authored file is what 9b's user brings.

   The test pins W0a's logged `sha256` for custody. This departs from item 289's "never in the tree", which was written for an in-process Gmsh STEP.
3. **`io/step_import.py`, two additive keys.** With both absent, behaviour is byte-identical (rule (c)).
   - **`"highest_dim_only": false`** makes the call `importShapes(str(path), highestDimOnly=False)`. When the key is absent, `:164` stays verbatim.
   - **`"surface_groups": [{"name": "port_1", "step_name": "PortSheet", "role": "port", "tag": 201}]`.**
     - Resolve each `step_name` by `getEntityName`: an exact match, or the last `/` component. Exactly one dim-2 entity must match. Resolve before or after the heal, whichever W0b shows the names survive.
     - Then call `out, out_map = occ.fragment(<all dim-3>, <named faces>)` and trace each named face's children through `out_map`. **This is identity by name plus topological lineage. A bbox or centroid selector stays forbidden.**
     - Require exactly one child, with area equal to the pre-fragment `getMass` to 1e-12, and upward adjacency equal to the volume. If W0b shows the child in the wall set, remove the named faces from the wall group by tag.
     - Raise `ValueError` through the existing broadcast for: a name absent (list the names present), a non-unique match, a non-embedded child, `role != "port"`, or a tag collision.
     - Add diagnostics `surface_groups[name] = {tag, role, step_name, surface, cad_area, n_triangles, route}`, and document both keys in the docstring.
4. **Config** `tests/mesh/data/geo34_step1d_port_sheet.json`:
   - `resolution` 0.08, asserted equal to `H`;
   - `air` tag 1 (`bbox_contains`);
   - `outer_boundary_tag` 101;
   - `highest_dim_only` false;
   - the `surface_groups` above.
5. **Probe and test.** `tests/mesh/probe_geo34_step1d_cq_read.py` (W0b; not collected) and `tests/mesh/test_step_import_cq_port_sheet.py` (a real build, at any width).

**W0b prints** (one rank, both reads):
- per-entity names and dims, before and after `removeAllDuplicates`;
- the `fragment` map, with each child's area, bbox and upward adjacency;
- whether the child is in `getBoundary(getEntities(3), combined=True, oriented=False, recursive=False)`;
- the bbox and the `sha256`.

**Decision:**
- Name found and the fragment child embedded ⇒ write the importer code.
- No name, or no embedded child ⇒ **stop, BLOCKED** (family 3 of 4), with no in-slot workaround.

**Anchors** (asserted at every width; every count reduced).
- **(P0)** `allreduce(index_map(3).num_ghosts)` is > 0 at `-n 2 / 4` and 0 at `-n 1`.
- **(P1) CAD identity.** `allreduce` the owned tag-201 triangle areas (filter by `index_map(2).size_local`). Assert |A / `RECT_AREA` − 1| ≤ `STRIP_BAND` and |A / `cad_area` − 1| ≤ `STRIP_BAND`, with `cad_area` = 0.48 to 1e-12.
- **(P2)** The reduced owned 201 count equals `n_triangles` exactly, and the reduced count of one-sided 201 facets is 0 (embedded, not a skin).
- **(P3)** `allreduce(assemble_scalar(form(1*dS(201))))`, with `quadrature_degree` pinned at 2, is within `STRIP_BAND` of 0.48.
- **(W)** The instrument check: the `ds(101)` area is within `STRIP_BAND` of 6.0.
- **(P4)** The gathered tag set is {101, 201}, and `validate_required_port_tags_exist` (`ports/definitions.py:82`) passes on 201.
- **Custody.** The file's `sha256` equals W0a's logged value.
- **Printed, not asserted:** the residuals (predicted ≤ 1e-14; item 267 read ≤ 1.6e-15) and `n_triangles`.

**Negative controls.**
- **(n-hDO)** The same config without `highest_dim_only` raises `ValueError` naming `PortSheet`. This is *predicted* from U12a's r1 reading (`20261008T214725Z_GEO-34-step1d-probe.log`), and *asserted* once W0b prints it on this fixture, with that line cited in the docstring (rule (e)).
- **(n1)** `step_name` `"PortSheetX"` raises on every rank, and the message contains the name.
- **(n2)** Without `surface_groups`, the reduced tag set lacks 201.
- **(n3)** Drop the first owned 201 facet on the lowest owning rank. Then |A′ / 0.48 − 1| > `STRIP_BAND`. The ceiling is ≈ 1 / `n_triangles` ≈ 5e-3 (item 267: −4.48e-3 / −5.07e-3), ≈ 5×10⁶ × the band. The factor is printed and labelled predicted.

**Windows.** Bare, foreground, Bash `timeout: 660000`. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 first. No complex build is needed, because every integrand is real.
- **W0a (author):** `scripts/testing/run_and_log.sh GEO-34-U12p-W0a-author "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 cq-python tests/mesh/cq_author_geo34_step1d.py tests/mesh/data/geo34_step1d_port_sheet.step'"` (≈ 6 s; the sample took 6 s).
- **W0b (read):** `scripts/testing/run_and_log.sh GEO-34-U12p-W0b-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 30 python3 tests/mesh/probe_geo34_step1d_cq_read.py tests/mesh/data/geo34_step1d_port_sheet.step'"` (≈ 1–5 s).
- **W1 (closing):** `scripts/testing/run_and_log.sh GEO-34-U12p-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_step_import_cq_port_sheet.py -v -s --tb=short && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 mpiexec -n 2 python3 -m pytest tests/mesh/test_step_import_cq_port_sheet.py -v -s --tb=short && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 mpiexec -n 4 python3 -m pytest tests/mesh/test_step_import_cq_port_sheet.py -v -s --tb=short'"`. Predicted 25–60 s (item 267: 56 s for seven launches). The `&&` stops at the first red width.
- **W2a (rule (c)):** the 4d command (that log's `:12`), labelled `GEO-34-U12p-W2a-ruleC`, plus `-v -s` and `tests/environment/test_complex_mode.py` on the path (104 s on record).
- **W2b (rule (c)):** the same command, labelled `GEO-34-U12p-W2b-preflight`, with path `tests/mesh/test_step_preflight.py` (117 s on record).

**Tier.** Standard by ceiling (180 s), smoke-sized by prediction; `-n 1 / 2 / 4`. **Class:** chain step (no solve). **Cost:** ≈ 10.8 measured (the class median) / ≈ 50 predicted min. About 20 min goes on code; the windows take ≈ 6 min.

**Traps already paid for.**
- **Harness discipline:**
  - `timeout -k 30` on every window;
  - every harness call is the whole Bash command (`bash_guard`): no pipes, `tee` or redirects;
  - `PYTHONPATH` carries `:/workspace`, and pytest runs with `-s`.
- **Parallel correctness:**
  - rank-local values (`assemble_scalar`, facet tag values) are `allreduce`d, and the tag set is gathered before the validator;
  - pytest verdicts are counted per rank, by token;
  - never delete a log.
- **Unghosted `dS` hangs at `-n 4`:** mesh only through `import_step`.
- **`fa4a6e7f`'s two check bugs:** gmsh 4.15.2 has no `-noenv`, and the default `highestDimOnly=True` drops a free sheet (known-issues 2026-10-08). Change `:164` only through the key.
- **Gmsh tags:**
  - `removeAllDuplicates` and `fragment` renumber tags: never carry a pre-heal tag past them, and trace through `out_map`;
  - OCC bbox padding is ≈ 1e-7: print bboxes, never select by them;
  - one gmsh session per process.
- **Forms and ranks:** pin `quadrature_degree` (`POST-5`), and broadcast a rank-0 throw to every rank (`GEO-23`).
- **Ownership:** files the container writes are owned by `nobody` on the host, so author only from inside the container.
- **Commit hygiene:**
  - after a kill, `pgrep -c python3` first;
  - no edit after W1 (rule (i));
  - never `--no-verify`. If the hook refuses the STEP, stop and report.

**Scope.** One planar rectangle in one box, authored by CadQuery 2.8.0: a name-to-`port` map, a native tag, and area, completeness and `dS` at three widths. **Not claimed:**
- U13 (the `Study` on a STEP-imported F-small);
- step 1c's F-small emitter;
- 9b closure;
- any port solve, S-parameter or band;
- roles other than `port`;
- multi-face or interface sheets;
- the preflight's facet grid.

**Negative result.** Report and stop. Park on `attempt/GEO-34-U12p-<UTC>` and mark the item BLOCKED with the unblock condition in the same commit (rule (d)). Open a known-issues entry.
- W0b finds no name or no embedded child ⇒ BLOCKED. The review may re-pose to a CadQuery-side fragment with a named sub-shape (the family's last attempt); this slot does not try it.
- (P1) or (P3) red with (P2), (n3) and (W) green ⇒ an import-tag defect: BLOCKED.
- Any (n·) or (W) red ⇒ the instrument is void. Fix it once and re-run in-slot; if it is still red, BLOCKED.
- W2 red ⇒ revert `src/`, then BLOCKED.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 U12′ → DONE, and U13 opens for the weekly to queue.
- The `GEO-34` §7 row gains a step-1d record clause. **The glyph stays 🟡** (1c's F-small emitter and step 3 stay open). No Done-when is added by this review.
- Known-issues 2026-10-08 (`highestDimOnly`) retires through its "Resolves with" clause if that clause names this read option; otherwise the review decides.
- Either outcome: dashboard item 10 retires, because the operator's install is done.

**Depends on:** nothing open. Independent of 983 / 984.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 985 — `GEO-34` step 1d / U12′ — per `docs/planning/chunks/GEO-34.md` § '§9 item 985'. Read `logs/maintenance/cadquery_image_20261114T091828Z.log` first and cite its 4a lines. Run W0a, then W0b, before any `src/` edit; if W0b finds no named face or no embedded fragment child, stop, BLOCKED. Identify the face by STEP name and the `fragment` map only — never a geometric selector. Both new config keys are additive, and key-absent behaviour is byte-identical (W2a / W2b). Mesh only through `import_step`. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Outcome (2026-11-14 12:00 CST slot): DONE, green (family 1d attempt 3 of 4).** Install log 4a: `cadquery 2.8.0` and `wrote …cq_sample_…step` (`cadquery_image_20261114T091828Z.log:381–382`). W0a (`20261114T181641Z_GEO-34-U12p-W0a-author.log:35–39`): `Assembly.export` exists, polygon face area 0.4800000000000001, 21468 B, sha256 `f729d87d88fee04de4f44524fe0edc204f273bf00123e99a6c37fce3c3fde392`. W0b (`20261114T181729Z_GEO-34-U12p-W0b-read.log`): the default read drops `PortSheet` (`:60–67`); `highestDimOnly=False` names face 7 `Shapes/geo34_step1d/PortSheet/PortSheet` before and after `removeAllDuplicates` (`:136`, `:174`); `fragment` maps (2,7) → [(2,7)] with the area unchanged, outside the wall set (`:182–183`). Upward adjacency reads `[]`, but `getEmbedded(3,1)=[(2,7)]` (`:184`) and all 115 sheet nodes lie in tets (`:237`); the importer checks embedding with `getEmbedded`. **Disclosed `src/` change:** the first W1 (`20261114T181903Z_GEO-34-U12p-W1.log:110–116`) died with SEGV in `optimize("Netgen")` ("456 illegal tets") on the embedded face, so with `surface_groups` set the Netgen pass is skipped; the key-absent path is unchanged. W1 (`20261114T181927Z_GEO-34-U12p-W1.log`, 11 s): at -n 1 / 2 / 4, (P1) +2.2e-16, (P3) ≤ 8.9e-16, (W) ≤ 5.2e-15, (P2) 192 = 192 and 0 one-sided, (P0) ghosts 0 / 581 / 1152, (n3) −4.06e-3 / −5.74e-3 / −5.48e-3 (predicted ≈ 5e-3), (n-hDO) / (n1) / (n2) green. Rule (c): W2a `20261114T181947Z_GEO-34-U12p-W2a-ruleC.log` (6 passed, 85 s) and W2b `20261114T182116Z_GEO-34-U12p-W2b-preflight.log` (5 passed, 141 s). Known-issues 2026-10-08 retired. U13 opens for the weekly.

**Not verified at review** (W0a / W0b measure these):
- the CadQuery 2.8.0 calls (`Assembly.export`, the polygon-to-face route);
- which dimension carries the product name after the gmsh import;
- whether names survive `removeAllDuplicates`;
- whether `occ.fragment` embeds an imported free face, and whether its child leaks into the wall group;
- the pre-commit hook's handling of a STEP file.

### §9 item 996 — 2026-11-15 (GEO-34 1e, lane 1 — the 11-14 weekly's target; mesh-probe)

*Posed by the 11-14 weekly at `PROJECT_PLAN.md:4995–4998` (the Netgen workaround "accepted as disclosed under rule (c), with its diagnosis posed as `GEO-34` 1e") and `:5300–5305` (lane 1: "the Netgen SEGV on the 1d fixture across a mesh-size ladder and with / without the embedded face — the illegal-tet count and whether the crash reproduces, a table and a probe script, no fix; it resolves known-issues 2026-11-14 or names its scope"). Brief drafted 2026-11-15 as execution detail. **Measurement only: no gate, no assertion, no band; a mesh-probe is 🧪 by §3 until a gate lands.** Family `GEO-34` 1e, attempt 1 of 4. **Executor:** `mesh-probe` (spawned by the slot).*

**What the step is.** U12′'s first W1 (`20261114T181903Z_GEO-34-U12p-W1.log`) died inside `gmsh.model.mesh.optimize("Netgen")` on the 1d fixture: "Remove Illegal Elements" / "456 illegal tets" (`:111–112`), "Caught signal number 11 SEGV" (`:116`), Status 15, 2 s. The workaround skips the pass whenever `surface_groups` is set (`src/fem_em_solver/io/step_import.py:320–330`). Nobody has measured whether the crash depends on h, on the embedded face, or on the CadQuery route. 1e runs the same mesh-then-optimise sequence on the 1d fixture over an h ladder, with and without the embedded face, and tabulates per rung: whether Netgen crashes, the illegal-tet count it reports, and the pre-optimise quality.

**Facts verified against the tree (this review, by `grep` / reading; unlogged, not a record).**
- **Fixture:** `tests/mesh/data/geo34_step1d_port_sheet.step` (sha256 `f729d87d…de392`, `tests/mesh/test_step_import_cq_port_sheet.py:38`), unit box `AirBox` + free face `PortSheet` = rectangle x 0.1, y 0.2, z 0.47, 0.8 × 0.6, area 0.48 (`tests/mesh/cq_author_geo34_step1d.py:19–20`). Config `tests/mesh/data/geo34_step1d_port_sheet.json`: `resolution` 0.08, `outer_boundary_tag` 101, `highest_dim_only: false`, one `surface_groups` entry `PortSheet` → 201.
- **The crash path:** `import_step` → `_build` (`importShapes`, `removeAllDuplicates`, synchronise; `step_import.py:154–166`) → fragment embeds the face (`getEmbedded(3,1)=[(2,7)]`, W0b `20261114T181729Z_GEO-34-U12p-W0b-read.log:184`) → `generate(3)` → `optimize("Netgen")` (now skipped under `surface_groups`, `:320–330`; key-absent path unchanged).
- **Without the face:** the same STEP with `surface_groups` removed takes the key-absent path and runs Netgen (U12′'s (n2), `test_step_import_cq_port_sheet.py:176–184`, green in W1 `20261114T181927Z_GEO-34-U12p-W1.log` at `-n 1/2/4`). So the "without" arm already ran Netgen once at h = 0.08 without a crash — with `highest_dim_only: false` the face is still read; whether it is still embedded on that path is a thing the probe prints, not assumes.
- **Known-issues 2026-11-14** (`docs/testing/known-issues.md:895–902`): cause "not diagnosed"; the slot's comment ("Netgen's optimiser does not respect an embedded face") is a hypothesis. **Resolves with:** "a measured quality reading on the un-optimised named-sheet mesh with the step's gates green [U13a, item 993], or an optimiser route that keeps the embedded face. Never by re-enabling the pass on the key-absent path."
- **A quality helper exists:** `_tet_quality(x)` (`tests/mesh/probe_wf6_sheet_facet_census.py:83–110`) — radius ratio and min dihedral per tet. In a gmsh-only probe use `gmsh.model.mesh.getElementQualities(tags, "minSICN")` / `"gamma"` instead (no DolfinX object needed); print which.
- **The 1d probe precedent:** `tests/mesh/probe_geo34_step1d_cq_read.py` (W0b), plain `python3`, `-n 1`, no `mpiexec` — "a gmsh throw on rank 0 at `-n 2` deadlocks the command" (`docs/planning/chunks/GEO-34.md:122`).

**The probe.** `tests/mesh/probe_geo34_step1e_netgen_ladder.py` (not collected), gmsh only, no DolfinX. One invocation = **one rung** `(h, arm)` given on argv, so a SEGV kills exactly one rung's process and the harness records its status. It replicates `import_step`'s gmsh sequence directly (`importShapes(highestDimOnly=False)`, `removeAllDuplicates`, synchronise; arm E then `occ.fragment` of the volume with the named face and synchronise), sets `Mesh.MeshSizeMin/Max = h`, `General.Terminal = 1` (so gmsh's own "illegal tets" line reaches stdout before a crash), runs `generate(3)`, prints the pre-optimise census (tets, min / p1 of minSICN, count with minSICN ≤ 0, `getEmbedded(3, vol)`), then `flush`, then `optimize("Netgen")`, then the post census. It **calls gmsh, never `import_step`**, so no `src/` file is touched and the workaround cannot mask the result. It prints `gmsh.__version__` and the OCC version.

**Arms × ladder (10 rungs + 2 repeats).**
- **Arms:** **E** — the face fragmented into the box (1d's embedded face, the crashing configuration); **N** — the same STEP read, face discarded before meshing (box only). Optional third arm, only if the first ten windows sum ≤ 120 s: **S** — the face extended to span the box (x 0..1, y 0..1 at z 0.47), splitting it into two volumes sharing the face, which is how a port sheet sits in an F-small port box (scope for U13a/U13b).
- **h ladder:** 0.16, 0.12, 0.08 (1d's), 0.06, 0.04. Stop going finer in an arm once a rung's `generate(3)` exceeds 60 s.
- **Repeats:** (E, 0.08) twice more — does the crash reproduce bit-for-bit (same illegal-tet count)?

**Table (the deliverable).** One row per rung: arm, h, tets pre, min minSICN pre, # tets minSICN ≤ 0 pre, embedded (Y/N from `getEmbedded`), Netgen "illegal tets" count (from the log line, `—` if none printed), crashed (Y/N; harness Status, 139/15/–11 named as printed), tets post and min minSICN post (if survived), elapsed s, `log:line`. Written into the probe's own stdout as a final block **and** into `docs/planning/chunks/GEO-34.md`'s Outcome paragraph.

**Windows** (bare and foreground, Bash `timeout: 660000`; `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0 first; **and again after every rung that crashed** — a SEGV can leave a half-dead rank; container-side ≤ 600 s; `timeout -k 30` mandatory). **Each rung its own window**, `-n 1`, plain `python3`, never `mpiexec`:

`scripts/testing/run_and_log.sh GEO-34-1e-<arm>-h<h> "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 90 python3 tests/mesh/probe_geo34_step1e_netgen_ladder.py --arm <E|N|S> --h <h>'"`

Label examples: `GEO-34-1e-E-h0p08`, `GEO-34-1e-E-h0p08-r2`. A crashed rung's non-zero Status is the measurement, not a failure of the window; do not re-run it to "get it green". **Rule (i):** no probe edit after its first rung's window (an edit restarts the ladder from rung 1, every earlier log kept).

**Price / Tier / Class / Cost.**
- **Smoke-to-standard by ceiling (90 s per rung)**, `-n 1`, real build, no DolfinX, no solve.
- Measured anchors: 1d's crashing window 2 s (`20261114T181903Z_…:126`), its green W1 11 s for six launches (`20261114T181927Z_…:729`); h = 0.04 is (0.08/0.04)³ ≈ 8× the tets, so ≈ 5–15 s per rung at worst. 12 windows ≈ 2–3 min of compute.
- **Class:** mesh-probe. **Cost:** **≈ 35 predicted** min (probe + 12 windows + table + journal).

**Traps already paid for.**
- **Bare harness only** (CLAUDE.md, `bash_guard.py`): the whole Bash command is the `run_and_log.sh` call; no `cd … &&`, no `| tail`, no `; date`, no redirect outside the container string.
- **No `mpiexec`:** a gmsh throw / SEGV on rank 0 at `-n 2` deadlocks (1d brief); `-n 1` plain `python3`.
- **A SEGV loses buffered stdout:** `flush=True` on every print and `General.Terminal = 1` before `optimize`; otherwise the illegal-tet line is lost with the process (the 1d log kept it, `:111–112`).
- **Never call `import_step`** for arm E: its named-surface branch skips Netgen by design, so it would read "no crash" by construction.
- **Never edit `src/fem_em_solver/io/step_import.py`**; never re-enable the pass on the key-absent path (known-issues "Never by …").
- Orphans: after a crashed rung, `pgrep -c python3` must read 0 before the next; if not, `docker compose -f docker/docker-compose.yml up -d --force-recreate` per known-issues, journal it.
- In tracked Markdown write `≈`, never a tilde.

**Scope.** Whether Netgen's optimiser crashes on the 1d fixture as a function of h and of the embedded face, with its illegal-tet count and the pre-optimise quality. **Not claimed:** a cause inside Netgen; a fix or an optimiser route (that is a second "Resolves with" half — posed by the review if the table points to one); anything about the F-small mesh (U13a, item 993); any assertion.

**Negative result / readouts (pre-registered).**
- **(a) E crashes at every h, N never:** the embedded face is necessary on this fixture; the known-issues *Cause* cell gets "embedded face, h-independent over 0.04–0.16"; the workaround's scope (every named-sheet import) stands as written.
- **(b) E crashes only at some h:** the scope is h-dependent; the cell names the band; U13a's F-small resolution is checked against it by the review.
- **(c) N also crashes at some h:** the face is not necessary — the crash is Netgen on this CAD; the key-absent path is exposed too; the review files that as a new known-issues line (the probe does not).
- **(d) No rung crashes, including (E, 0.08) repeats:** the 11-14 crash does not reproduce; report the gmsh version and the 1d log's header side by side and stop — that is the finding.
- Any of these is a complete outcome. No `PREDICTED RED:` line.

**Status it can move.** `GEO-34` 1e → 🧪 (probe landed). Known-issues 2026-11-14: its *Cause* and *Scope* cells gain the measured table (the review writes them). It **does not retire** the entry: the entry's own "Resolves with" names U13a's quality reading with gates green, or an optimiser route that keeps the face — 1e is neither. (The weekly's "resolves … or names its scope" therefore lands on "names its scope" unless the review rules (d) a resolution.)

**Depends on:** nothing open. Independent of item 993 (no shared file edited); may share a slot with it.

**Outcome commit.** The 12 (or up to 17) logs, the probe, the table in `docs/planning/chunks/GEO-34.md`, and the §9 line marked DONE with the local (CST) slot, in one commit through the hook. Known-issues is the review's edit.

**Executor spawn prompt (verbatim).** "mesh-probe: `GEO-34` 1e, §9 item 996, per `docs/planning/chunks/GEO-34.md` "§9 item 996". Generator: gmsh on `tests/mesh/data/geo34_step1d_port_sheet.step` (the 1d fixture), replicating `import_step`'s read (`highestDimOnly=False`, `removeAllDuplicates`) — never calling `import_step`. Swept: h ∈ {0.16, 0.12, 0.08, 0.06, 0.04} × arm {E: face fragmented into the box, N: box only}, plus (E, 0.08) twice more; arm S only if the first ten windows sum ≤ 120 s. Question: does `optimize(\"Netgen\")` after `generate(3)` crash, and what illegal-tet count and pre-optimise minSICN does each rung show? Write `tests/mesh/probe_geo34_step1e_netgen_ladder.py` (one rung per invocation from argv, `flush=True`, `General.Terminal=1`, prints gmsh version). Each rung its own window: the harness call is the whole Bash command, bare and foreground, Bash `timeout: 660000`, container `timeout -k 30 90`, `-n 1` plain `python3`, never `mpiexec`; `pgrep -c python3` reads 0 before every window and after every crash. A crashed rung's Status is the measurement — never re-run to get it green. Edit no `src/` file and no test; no assertion, no band, no fix. Return the table with `log:line` per row and which readout (a)–(d) it supports. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Outcome (2026-11-15 16:30 CST slot): DONE, probe landed, readout (a) (family 1e attempt 1 of 4; executor `mesh-probe`; 🧪 by §3, measurement only).** The probe is `tests/mesh/probe_geo34_step1e_netgen_ladder.py`, on gmsh `4.15.2-git-657c8e9` (OpenCASCADE, `20261115T223140Z_GEO-34-1e-E-h0p16.log:36–38`). It ran 17 windows: 5 E, 5 N, 2 E repeats and 5 S. Arm S ran because the first ten summed ≈ 20 s, under the 120 s limit. Each window was `-n 1`, plain `python3`, `timeout -k 30 90`. `pgrep -c python3` read 0 before every window and after every crash, so no recreate was needed. The probe's one edit came before rung 1's window, so rule (i) holds. Statuses were read by the slot from every log's footer; the E h = 0.08 and N h = 0.08 lines were re-read by the slot.

| arm | h | tets pre | minSICN min pre | n ≤ 0 pre | embedded | Netgen illegal tets (first pass) | first badmax | Status | tets post | minSICN min post | log:line |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E | 0.16 | 1714 | 0.3012 | 0 | `[(2,7)]` | 241 | 1e+24 | 139 | — | — | `20261115T223140Z_GEO-34-1e-E-h0p16.log:98,104,106,110` |
| E | 0.12 | 3637 | 0.2942 | 0 | `[(2,7)]` | 314 | 1e+24 | 139 | — | — | `20261115T223154Z_GEO-34-1e-E-h0p12.log:99,105,107,111` |
| E | 0.08 | 10316 | 0.2627 | 0 | `[(2,7)]` | 456 | 1e+24 | 139 | — | — | `20261115T223203Z_GEO-34-1e-E-h0p08.log:100,106,108,112` |
| E | 0.06 | 23675 | 0.2993 | 0 | `[(2,7)]` | 604 | 1e+24 | 139 | — | — | `20261115T223211Z_GEO-34-1e-E-h0p06.log:105,111,113,117` |
| E | 0.04 | 73087 | 0.2976 | 0 | `[(2,7)]` | 900 | 2.68605e+32 | 139 | — | — | `20261115T223220Z_GEO-34-1e-E-h0p04.log:118,124,126,130` |
| E r2 | 0.08 | 10316 | 0.2627 | 0 | `[(2,7)]` | 456 | 1e+24 | 139 | — | — | `20261115T223323Z_GEO-34-1e-E-h0p08-r2.log:100,106,108,112` |
| E r3 | 0.08 | 10316 | 0.2627 | 0 | `[(2,7)]` | 456 | 1e+24 | 139 | — | — | `20261115T223331Z_GEO-34-1e-E-h0p08-r3.log:100,106,108,112` |
| N | 0.16 | 1564 | 0.3104 | 0 | `[]` | 236 → 0 | 13.59 | 0 | 1801 | 0.4836 | `20261115T223231Z_GEO-34-1e-N-h0p16.log:93,99,136,143` |
| N | 0.12 | 3411 | 0.3025 | 0 | `[]` | 308 → 0 | 16.31 | 0 | 3615 | 0.4529 | `20261115T223238Z_GEO-34-1e-N-h0p12.log:94,100,137,144` |
| N | 0.08 | 10194 | 0.2981 | 0 | `[]` | 452 → 0 | 16.44 | 0 | 9906 | 0.4833 | `20261115T223246Z_GEO-34-1e-N-h0p08.log:96,102,139,146` |
| N | 0.06 | 22668 | 0.2966 | 0 | `[]` | 599 → 0 | 17.14 | 0 | 21057 | 0.4802 | `20261115T223255Z_GEO-34-1e-N-h0p06.log:99,105,142,149` |
| N | 0.04 | 72230 | 0.3002 | 0 | `[]` | 884 → 0 | 16.20 | 0 | 65197 | 0.4494 | `20261115T223306Z_GEO-34-1e-N-h0p04.log:113,119,156,163` |
| S | 0.16 | 1853 | 0.3472 | 0 | `[]` (2 vols) | 188 / 200 → 0 | finite | 0 | 2279 | 0.4635 | `20261115T223359Z_GEO-34-1e-S-h0p16.log:111,118,188,196` |
| S | 0.12 | 3690 | 0.3153 | 0 | `[]` (2 vols) | 248 / 260 → 0 | finite | 0 | 4235 | 0.5069 | `20261115T223407Z_GEO-34-1e-S-h0p12.log:111,118,192,200` |
| S | 0.08 | 10501 | 0.2895 | 0 | `[]` (2 vols) | 368 / 380 → 0 | finite | 0 | 10620 | 0.4328 | `20261115T223415Z_GEO-34-1e-S-h0p08.log:114,121,195,203` |
| S | 0.06 | 23519 | 0.2978 | 0 | `[]` (2 vols) | 489 / 501 → 0 | finite | 0 | 22656 | 0.4721 | `20261115T223425Z_GEO-34-1e-S-h0p06.log:117,124,194,202` |
| S | 0.04 | 73861 | 0.2699 | 0 | `[]` (2 vols) | 730 / 759 → 0 | finite | 0 | 67925 | 0.4801 | `20261115T223434Z_GEO-34-1e-S-h0p04.log:130,137,211,219` |

**Readout (a).**
- **E crashes at every h over 0.04–0.16, with SIGSEGV** (Status 139, "the monitored command dumped core"). N and S never crash.
- **The crash reproduces bit for bit.** (E, 0.08) gave identical tets, minSICN, 456 illegal tets and `badmax` three times, and 456 equals the 11-14 crash's count (`20261114T181903Z_GEO-34-U12p-W1.log:111–112`). The 11-14 crash read Status 15 under `mpiexec`; plain `python3` reads 139.

**Readings, not gated (for the review's *Cause* cell):**
- **The pre-optimise mesh is healthy in every arm.** No tet has minSICN ≤ 0, and the minimum is 0.26–0.35.
- **"Illegal tets" is printed by every arm.** E's first-pass count is within ≈ 2 % of N's at the same h. N and S drive it to 0 within six passes. So the count is Netgen's own criterion, not the crash signature.
- **What only E shows is a non-finite-scale first-pass `badmax`** (1e+24, and 2.68605e+32 at h = 0.04), immediately followed by the core dump. N's first badmax is 13.6–17.1.
- **The embedded face's triangles may not reach Netgen.** At h = 0.16 its `CalcLocalH` reports the same surface-element count in E as in N (708; `E-h0p16.log:102`, `N-h0p16.log:97`). This suggests Netgen does not see the embedded face's triangles as boundary. It is an inference, not a measurement.
- **Arm S is the F-small port-box shape** (a face that splits the box into two volumes sharing it). It survives Netgen at every h, which bears on U13a's option (a) / (b) (item 993).

**Caveats.**
- The N arm deletes the face before meshing. The importer's key-absent path keeps it as an unembedded free face, so the two are related but not identical.
- The crashed rungs print no final ROW block; their rows come from the census and Netgen lines cited.

**Not claimed:** a cause inside Netgen, a fix or an optimiser route. Known-issues 2026-11-14 is not retired; its *Cause* and *Scope* cells are the review's edit. No `src/` or test edit.
