# GEO-36 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 150 — 2026-10-01 (chain I5, first half — the `mesh-probe`; queued by the 10:30 review)

`GEO-36` step 1 — rod-in-phantom meshability and price over the rod radius (`mesh-probe`; standard, `-n 1`, real; ≈ 6 measured / ≈ 35 predicted min; depends on nothing). **Split by the 10:30 review:** `mesh-probe` may not touch `src/` and no rod generator exists, so this probe runs on a reduced geometry (the F-small phantom alone, rod cut out). The CAD-identity half — an additive `implant_rod=` keyword on `birdcage_port_domain` with the identities asserted at `-n 2` — is the next restock's item, written from this table.

- **Question.** Does a capped rod cut as a hole from the F-small phantom mesh with graded tips, and what does it cost, across the rod radius?
- **Sweep.** Swept: r ∈ {0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00} mm (uniform grid, never bisect). Fixed: phantom R = 0.03 m, H = 0.08 m (`tests/mesh/test_birdcage_port_tags.py:28–29`); rod 30 mm cap-to-cap with hemispherical caps, axis ẑ, centred at (0.015, 0, 0). Sizing: phantom at `RESOLUTION` 0.015; Distance→Threshold from the rod wall (r/3 → 0.015 over r → 10r); `Ball` fields at the cap centres (radius 2r, VIn r/4).
- **Script.** `tests/mesh/probe_geo36_rod_in_phantom.py`, gmsh only, one fresh subprocess per rung (template `probe_birdcage_conductor_resolution.py:127–166`).
- **Printed per rung.** OK/FAIL; tets; gmsh seconds; minimum and 1st-percentile SICN; the occ mass against πR²H − (πr²ℓ + 4πr³/3) and the rod wall area against 2πrℓ + 4πr²; the rod-wall facet centroid's y and z.
- **Predicted (never asserted).** Tets ≈ 33 k·(1 mm/r) above the phantom-only baseline; all rungs OK; each rung under 10 s. No prior rod measurement exists, so the costs are unmeasured.
- **Windows.** W1: `scripts/testing/run_and_log.sh GEO-36-probe-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 python3 -u tests/mesh/probe_geo36_rod_in_phantom.py'"`. W2: the same command relabelled `GEO-36-probe-W2`, to show the sweep is bit-identical (the determinism rule). Predicted ≈ 2–3 min each.
- **Traps.** `-n 1` only; a fresh process per rung (a failed gmsh build contaminates the process, `GEO-23`). The cavity group from the meshed boundary. Never retry a FAIL — record it verbatim with its surface pairs. Foreground, `timeout: 660000`.
- **Scope.** Reduced geometry; the full-fixture price belongs to the second half. No assertion, no `src/` change.
- **Status it can move.** `GEO-36` ⬜ → 🧪 with the table in the row; the CAD-identity half becomes queueable with a chosen r.
- **Negative result.** A FAIL is recorded verbatim with its surface pairs; the table stands as written.

**Outcome 2026-10-02 (19:30 CDT slot): DONE (measurement).** All 7 rungs are OK. Tets above a 902-tet baseline are 92 321 / 59 516 / 42 115 at r = 0.5 / 1 / 2 mm, which misses the predicted 33k·(1 mm/r) by 1.4–2.6×; the measured slope is ≈ r^-0.57. Gmsh takes ≤ 2.9 s per rung. Min SICN runs 0.0695 (0.5 mm) → 0.2914 (2 mm), with p1 at 0.356–0.370. The occ mass matches πR²H − (πr²ℓ + 4πr³/3) to ≤ 2.4e-16, and the meshed wall area is 0.36–0.41 % short of 2πrℓ + 4πr². The cavity surfaces are the same three on every rung. W1 and W2 are bit-identical bar wall-clock. Reading taken: ℓ = 30 mm is the cylindrical section between the cap centres (tip to tip 30 + 2r mm); the CAD-identity half should restate this. `20261002T003452Z_GEO-36-probe-W1.log`, `20261002T003510Z_GEO-36-probe-W2.log`.

### §9 item 171 — 2026-10-02 (chain I5, second half — CAD identities; queued by the 03:00 review)

`GEO-36` step 2. Add a keyword-only `implant_rod=None` to `MeshGenerator.birdcage_port_domain` and to `_build_birdcage_port_model`. It cuts a capped rod out of the F-small phantom as a hole with its own wall facet tag. Assert the CAD identities in a new `tests/mesh/test_birdcage_implant_rod.py`.
- Executor: implementer. Tier: standard. Build: real. Ranks: `-n 1` and `-n 2`. No solve.
- Cost: ≈ 8 measured / ≈ 45 predicted min. Depends on nothing; the probe landed `9cc8a238`, item 150.
- Drafted by a read-only `Plan` agent and reviewed by the 03:00 review. The row's Done-when matches as written.

**Restated geometry.** ℓ = 30 mm is the cylinder between the cap centres, so tip to tip is 30 + 2r mm. The rod has r = 1 mm, axis ẑ, centre (0.015, 0, 0), on leg 1's azimuth. r = 1 mm is the smallest rung with min SICN ≥ 0.22 (`20261002T003452Z_GEO-36-probe-W1.log:49`). The fixture is the four-port F-small call that `_build(True)` makes (`tests/mesh/test_birdcage_port_sheets.py:136–160`), called directly with the same imported constants. `_build` is not edited, so there is no rule (a) cascade.

- **Build.** Follow the `as_hole` cut-first sequence (`src/fem_em_solver/io/mesh.py:4467–4494`).
  - Fuse the shaft and the two spheres (probe `:77–83`) and record the rod's OCC mass.
  - Cut the rod from **both** the phantom and the air box (`removeTool=False`) before the fragment. If it is cut from the phantom alone, the air box fills the hole.
  - Read the wall off the meshed volumes' boundary: single-use faces inside the rod's padded bbox (probe `:93–101`).
  - Emit it as an exterior dim-2 group `BIRDCAGE_IMPLANT_SURFACE_TAG = 402`, beside `:91`. Confirm 402 is free first.
  - Append the probe's sizing to `size_fields` (`:4956`, Min at `:5041`): Distance → Threshold r/3 → `resolution` over r → 10r, plus two `Ball`s of radius 2r with VIn r/4.
  - Return `implant_rod_cad`: the surfaces, the OCC wall area, the OCC rod mass, the area-weighted OCC wall centroid, and the tag (`None` when off).
  - Refuse with `ValueError` on `as_hole`, on `symmetry_quadrant`, on r or ℓ ≤ 0, and on a rod not strictly inside the phantom.
- **Anchors (all asserted; bands imported, never restated).**
  - (A1) `cad_mass_by_group["phantom"]` = πR²H − (πr²ℓ + 4πr³/3) within `CAD_ANALYTIC_BAND` 1e-7 (`tests/mesh/test_birdcage_conductor_hole.py:123`). Predicted ≤ 2.4e-16 (probe W1:48).
  - (A2) Σ CAD(rod) + V_rod,occ = Σ CAD(flag off) within `EXACT` 1e-9 (`:117`).
  - (A3) OCC wall area = 2πrℓ + 4πr² within `CAD_ANALYTIC_BAND`. Predicted 0 (W1:49).
  - (A4) The meshed tag-402 area ÷ (2πrℓ + 4πr²) lies inside `TERMINAL_AREA_BAND` (0.95, 1.0) (`tests/mesh/test_birdcage_leg_gaps.py:73`), the imported inscribed-triangulation band. The probe measured 0.99625 (W1:49), so the band's lower edge sits ≈ 13× below the measured deficit.
  - (A5) The tag census:
    - the cell-tag set equals the flag-off set {1, 2, 3, 101–104, 201–204};
    - the owned census sums to `size_global` (`_owned_census`, `test_birdcage_conductor_hole.py:130`);
    - tag 402 has > 0 facets, each exterior, with its single cell in tag 3;
    - tag 402 is absent with the flag off;
    - the flag-off cells equal `SOLID_CELL_RECORD` 116 085 (`:95`) within `CELL_COUNT_BAND` 0.01.
  - (A6) Mirror symmetry: the OCC wall centroid has \|ȳ\|, \|z̄\| ≤ `EXACT`·ℓ, and x̄ = 0.015 m to the same bound. All three are printed. The meshed facet centroid is printed only.
  - (A7) Each refusal raises without meshing.
- **Printed.** Cells, vertices, tag-402 facets, `mesh_wall_time_s`, elapsed. The cell count is a first record, predicted ≈ 176 k and not asserted.
- **Negative controls (rule (e)).**
  - (C1) The flag-off build misses A1's comparand by V_rod/(πR²H − V_rod) = 4.354e-4. **Asserted** as a miss > `CAD_ANALYTIC_BAND`; the probe measured the same OCC phantom solid both ways (W1:39 against W1:48). The ceiling is that same 4.354e-4, so the largest possible separation is 4.35e3×; the factor is **predicted** and printed.
  - (C2) **Predicted:** a fragment that refilled the hole with air would miss A2 by V_rod/V_box ≈ 9.5e-6, about 9.5e3× `EXACT`. This is the failure A2 exists to catch, not a separate build.
  - (C3) **Predicted, printed:** the x-mirror is not a symmetry. \|x̄\| = 0.015 m against `EXACT`·ℓ.
- **Windows** (Bash `timeout: 660000`; foreground; bare and whole). Basis: the probe ran 15 s / 14 s, with gmsh ≤ 2.9 s per rung, and `TH-15` step 3a built this fixture twice in 54 s at `-n 2` (`20260906T213913Z_TH-15.log:2655`).
  - W1, predicted ≈ 90 s: `scripts/testing/run_and_log.sh GEO-36-step2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 1 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod.py -s -v --tb=short'"`
  - W2: the same at `-n 2`, labelled `GEO-36-step2-W2`. This is the Done-when window.
  - W3, the rule (c) re-gate, predicted ≈ 75 s: `scripts/testing/run_and_log.sh GEO-36-step2-W3-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_port_tags.py tests/mesh/test_birdcage_conductor_hole.py -s -v --tb=short'"`. The 116 085 record must print unchanged.
- **Traps already paid for.**
  - Read the wall from the meshed boundary, never from the tool's faces (`TH-15` step 0: 716 orphan nodes).
  - A failed gmsh build contaminates its process (`GEO-23`): build flag-off first, then the rod.
  - `cell_tags.values` is rank-local: use the owned census (`OPS-39`).
  - Use `-s`; no `-k` inside the quoted command; no pipes.
  - The closing window runs the module exactly as committed (rule (i)).
  - The generator also runs `optimize("Netgen")` (`:5102`), which the probe did not. A4's margin covers that.
- **Scope.** None of the following is in this item:
  - a solve, a PEC condition, or SAR;
  - `as_hole` combined with the rod (the coil-cavity rule at `:4655` would swallow the rod wall into 401; the weekly's call if `IMP-3` needs it);
  - F-human.

  `IMP-3` is the consumer.
- **Status it can move.** §7 `GEO-36` 🧪 → ✅, on the row's Done-when: the identities in one logged `-n 2` window (W2), with cells and time recorded. `IMP-3` (I6) becomes queueable.
- **Negative result.**
  - A gmsh FAIL is recorded verbatim with its surface pairs. Write a known-issues entry, park on `attempt/*`, and mark the item BLOCKED in the same commit (rule (d)).
  - An A4 miss is reported and never widened. Rule (f) does not apply, because the band is ≥ 10× the measured deficit, so a miss means the sizing was not applied.
  - If W3 moves the default path, revert `src/` and stop.
  - No `PREDICTED RED` line.

**Outcome 2026-10-02 (06:00 CDT slot): DONE, 🧪 → ✅.**
- **Code.** `BIRDCAGE_IMPLANT_SURFACE_TAG = 402` (free before; 401 was the nearest). `_normalise_implant_rod` validates the request before any gmsh state exists. `implant_rod=` is keyword-only on `birdcage_port_domain` and `_build_birdcage_port_model`, and the diagnostics carry `implant_rod_cad` (`None` with the flag off). The `src/` diff is additive only (+240 lines in `io/mesh.py`). New `tests/mesh/test_birdcage_implant_rod.py`.
- **Windows.** W1 `-n 1`: 22 passed, 4 skipped, 124 s (`20261002T110348Z_GEO-36-step2-W1.log`). W2 `-n 2`, the Done-when window: 22 passed, 4 skipped per rank, 121 s (`20261002T110559Z_GEO-36-step2-W2.log`). Both files predate W1, so the module ran as committed.
- **Anchors (W1 = W2).**
  - (A1) 1.110e-16.
  - (A2) 0.
  - (A3) 0.
  - (A4) 2.003087e-4 / 2.010619e-4 = 0.996254, which matches the probe's 0.99625; `optimize("Netgen")` did not move it.
  - (A5) Cell-tag sets are equal both ways. The owned censuses sum to 167 371 / 116 085. Tag 402 has 4 958 facets, all exterior and all on a tag-3 cell, and 0 with the flag off. Flag-off / 116 085 = 1.000000.
  - (A6) \|x̄ − 0.015\| 1.7e-18, ȳ −3.8e-20, z̄ 6.6e-20. The meshed centroid, printed only, is (0.01500005, −2.1e-8, −3.9e-7).
  - (A7) 8/8 refusals raise without meshing.
- **Controls.**
  - (C1) 4.353747e-4, 4.35e3× the band (asserted).
  - (C2) printed 8.54e-6 (8.54e3× `EXACT`). The brief's ≈ 9.5e-6 used the wrong box volume; the true box sum is 1.152e-2 m³.
  - (C3) printed 5.0e8×.
- **Records.** The rod build has 167 371 cells (≈ 176 k predicted, 5 % under), 31 308 vertices, and a 42 s mesh (45 s elapsed). The wall is 3 OCC surfaces.
- **W3 re-gate** (`20261002T110806Z_GEO-36-step2-W3-regate.log`, 98 s): 16 passed and 1 failed. The 116 085 record prints unchanged, and all five `test_birdcage_conductor_hole.py` tests pass. The one red is `test_birdcage_volumes_partition_the_box`, failing with "Invalid boundary mesh (overlapping facets) on surface 59 surface 79". That is verbatim the 🔴 OPEN known-issues entry of 2026-08-27 (`OPS-26` / `GEO-23`) on the ungapped default build, which the new code never reaches, so the default path has not moved. A "FAILED" next to the conductor-hole test at log line 392 is rank interleaving; the summary lists only the one failure.

### Step 3 — 2026-10-04 (§9 item 219, `WF-9` step 1 Part A; brief in `docs/planning/chunks/WF-9.md#§9-item-219`)

**Outcome — DONE 2026-10-04 (15:00 CDT slot, implementer).** Additive `implant_rod={..., "keep_cells": True}`: the rod is cut from the phantom and the air box with `removeTool=False`, enters the fragment once as its own tool, and its pieces carry the new cell group `BIRDCAGE_IMPLANT_CELL_TAG = 403` (grep of `src/` for 403 empty; free of 1/2/3, 100+i, 200+i and of facet groups 1/401/402/501/502). The wall is the boundary of the implant pieces (all two-use, checked; a non-two-use face raises). `implant_rod_cad` gains `keep_cells`, `cell_tag` and `kept_mass_m3`. A non-bool `keep_cells` raises in `_normalise_implant_rod` before gmsh. Absent / `False` is the step-2 hole path.

**Disclosed deviation from the brief.** The brief asked to emit 402 as a dim-2 gmsh group on interior facets. `io/mesh.py`'s own `_interface_facet_tags` docstring records that interior dim-2 gmsh groups hang `model_to_mesh` at `-n 2` (known-issues 9, 2026-08-05), and every interior facet group in the module follows that rule. So in the kept branch the gmsh group is not emitted; 402 is rebuilt dolfinx-side after `model_to_mesh` as the (3, 403) interface, merged into `facet_tags`. Downstream sees the same tag on the same facets. The hole path still emits its exterior gmsh group unchanged.

Windows (all Status 0, heavy tier):
- W1 `-n 1`: `20261004T200412Z_GEO-36-step3-W1.log`, 29 passed / 4 skipped, 188 s.
- W2 `-n 2` (of record): `20261004T200730Z_GEO-36-step3-W2.log`, 29 passed / 4 skipped, 181 s.
- W3 default-path re-gate `-n 2`: `20261004T201038Z_GEO-36-step3-W3-regate.log`, 17 passed / 4 skipped, 114 s; the 116 085 record prints unchanged.

Measured, identical at W1 and W2:
- (K1) 1.1e-16; (K2) 0 (CAD implant = OCC rod mass = V_rod); (K3) 0; (K4) OCC wall 0, meshed dS(402) / closed form 0.996254 (predicted 0.99625, the hole build's); (K6) \|x̄ − 0.015\| 1.7e-18, ȳ −3.8e-20, z̄ 6.6e-20; meshed dS centroid (printed) (0.01500005, −2.1e-8, −3.9e-7) — the hole build's to every printed digit.
- (K5) Tag set = flag-off ∪ {403}; census sums to 182 371; 15 000 implant cells; 402 has 4 958 owned facets, 0 exterior, 4 958 interior, all between one tag-3 and one implant cell (the hole build's 4 958 exterior facets, now interior); 403 absent flag-off and in the hole build.
- (K7) `keep_cells="yes"` and `keep_cells=1` refused in 0.000 s.
- A1–A7 re-run unchanged on the hole build; 167 371 prints again.
- (C1k, asserted) hole build misses K3 by 8.544841e-6 = V_rod/ΣCAD (8.5e3× `EXACT`). (C2k, predicted, printed) dropping 403 misses the census by 15 000 = 1.0000 of the implant count.

Records:
- `KEPT_ROD_CELL_RECORD = 182 371` = 167 371 + 15 000; every non-implant tag count is identical to the hole build's, cell for cell. Registered from W1, printed 1.000000 by W2. Guard 320 000.
- Meshed implant volume / V_rod = 0.986667, **outside** the printed prediction [0.990, 0.999] (printed only; moves nothing). The inscribed-faceting deficit of a 1 mm rod at h ≈ r/3 is larger than the brief guessed; the wall area deficit is 0.37 %, the volume deficit 1.33 %.
- Mesh wall time 46.6 s (W1) / 46.8 s (W2); build elapsed ≈ 50 s.

### Audit — 2026-10-02 10:30 daily review (`auditor`, closure `73f65d05`)

**PASS.**
- **Digits:** 41 / 41 traced. The review re-read W2:3627 ((A4) 0.996254, band (0.95, 1.0)).
- **W3 red attribution holds.** The red is the 🔴 OPEN known-issues entry of 2026-08-27 (`known-issues.md:145`, surfaces 59 / 79):
  - the same test and string appear in `20261002T110806Z_GEO-36-step2-W3-regate.log:3085`;
  - the same failure is in `20260828T140313Z_GEO-23-step1b-birdcagepart-n1.log.gz:357–361` on commit `200916a5`;
  - every new `src/` hunk sits behind `implant_rod is not None`.

Caveats:
- The red baseline is five weeks old, and no pre-change re-run was made.
- The Done-when's "mirror symmetry of the placement" is met at CAD level only: (A6) asserts the OCC centroid, and the meshed centroid is printed.
- (A4)'s band is imported and carries no control.

## §7 row moved verbatim 2026-10-03 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/GEO-36.md`)* **A rod in the phantom: the implant-in-phantom generator (Phase 7 subgoal (7b))** (opened 2026-09-30 by the weekly; §10 chain I5). A capped cylindrical rod inside the F-small phantom, emitted as a hole (the `TH-15` route) with its own wall tag, tips graded; parameters generic (length, radius, position, axis) — nothing from a specific product. `mesh-probe` first, on the reduced geometry (phantom alone, rod cut out; the probe writes no `src/`): meshability and price over r ∈ [0.5, 2] mm (§9 item 150, 2026-10-01 10:30 review); then the implementer adds an additive `implant_rod=` keyword on `birdcage_port_domain` and asserts the CAD identities: phantom volume minus rod volume and the rod's wall area against their closed forms at the imported bands, the tag census, the mirror symmetry of the placement. **Done-when (§4):** the identities asserted in one logged window at `-n 2`, cells and time recorded ⇒ ✅. No solve. | ✅ **Audit PASS 2026-10-02 10:30 review (chunk file).** **Step 2 2026-10-02 (§9 item 171):** `implant_rod=` on `birdcage_port_domain` (r = 1 mm, ℓ = 30 mm between cap centres, wall tag 402). At `-n 2` (W2, 22 passed, 121 s), phantom CAD vs πR²H − V_rod is 1.1e-16, the box sum is exact, the OCC wall area is exact, the meshed wall ratio is 0.996254, the census is clean (4 958 exterior tag-402 facets), and the centroid sits on the mirror planes to ≤ 1.7e-18 m; 167 371 cells, 42 s mesh (`20261002T110559Z_GEO-36-step2-W2.log`). Step 1 history below. 🧪 **Step 1 probe 2026-10-02 (§9 item 150, measurement only):** on the reduced geometry (phantom alone, 30 mm shaft + hemispherical caps at (0.015, 0, 0)), all seven rungs r = 0.50 … 2.00 mm are OK. Tets 93 223 / 60 418 / 43 017 at r = 0.5 / 1 / 2 mm, over a 902-tet baseline. That is 1.4–2.6× the predicted 33k·(1 mm/r), and it falls as ≈ r^-0.57. Gmsh takes ≤ 2.9 s per rung. Min SICN runs 0.07 → 0.29 while p1 holds 0.36–0.37. The occ mass matches its closed form to ≤ 2.4e-16, and the meshed wall area is 0.36–0.41 % short. W1 and W2 are bit-identical (`20261002T003452Z_GEO-36-probe-W1.log:40–61`, `…003510Z_GEO-36-probe-W2.log`, 15 s / 14 s; `tests/mesh/probe_geo36_rod_in_phantom.py`). The CAD-identity half is queueable with a chosen r. 

### §9 item 298 — 2026-10-07 (chain I15, part 1 — the tip-refinement knob; no solve)

*Drafted for the 2026-10-07 daily review; the review reviews it.* Licence: §10 I15 (`PROJECT_PLAN.md:4682–4687`), weekly-enumerated, serial on I14, which landed (item 284, `IMP-3.md:947–951`); §9 rule (4). Only I15's mesh half runs here: the knob, the CAD identities at two refined rungs, and cells / dofs so the weekly can price the solves (heavy vs `xl` is the weekly's, `IMP-3.md:899`). **No solve; no band, tolerance or target moves.**

- **Executor:** implementer. **Tier:** heavy, `-n 2`, real build. **Cost:** ≈ 10 measured / ≈ 50 predicted min (15 fixed + 15 code + 15 compute + 5 spare). **Depends on:** nothing.
- **Fixture.** I15's rod is the *hole* build `_build_rod(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` (`test_imp3_tip_ratio_1g.py:121–128`; 0.0025 at `test_birdcage_sar_1g_rung.py:114`, 0.5 at `test_port_lumped_narrowed_sheet.py:132`). So its identities are step 2's A1–A6 (`test_birdcage_implant_rod.py:374–474`). K1–K6 belong to the kept-cells build, which is `WF-9`'s. Record: 252 118 rod / 201 072 off cells, 4 958 tag-402 facets, rod mesh wall 58.77 s (`IMP-3.md:913`, `20261006T111102Z_IMP-3-item263-W1-build.log`). It read 252 118 again at `-n 8` (`IMP-3.md:948`).

**The knob: `implant_rod["tip_size"]` (metres), the cap-Ball `VIn` itself.** It is in metres, not a factor, because every other sizing input of `birdcage_port_domain` is an absolute length, and a length makes each step ratio a plain ratio of requested h.
- `_normalise_implant_rod` (`io/mesh.py:4589–4654`) returns `"tip_size"`: `None` when the key is absent, else a float. Before gmsh starts it refuses, with a `ValueError` naming `implant_rod`, any bool, non-number, non-finite value, value ≤ 0 or value > `radius / 4.0`. It runs twice (`:4257`, `:4690`), so `None` must pass the second call unchanged.
- Sizing (`:5900–5933`): both Balls take `VIn` = the value (`rod_r / 4.0` when `None`, today's expression). The radius 2r (`:5924`), the Threshold and the print stay as they are, and the print also names the value used.
- `implant_rod_cad` (`:5490–5498`) gains `tip_size_m` and `tip_ball_radius_m` (additive keys), so the test reads them back and never restates them.
- **Rungs:** r/(4√2) and r/8 (1.7678e-4 m and 1.25e-4 m). r/4 → r/(4√2) → r/8 has a constant step ratio √2, which an observed-order reading needs.
- **r/16 is argued out.** The cubic model's ceiling puts the whole rod increment in the balls: 51 046 cells (252 118 − 201 072). That makes r/16 ≤ 3.47 M cells, or ≈ 1.8 M if the balls hold half. Serial gmsh runs at ≈ 4.3 k cells/s (252 118 / 58.77 s), so meshing alone takes ≈ 7–13 min, and the ≈ 12 M degree-2 dofs fit no tier. The measured r/8 growth prices r/16 for the weekly; it is printed.

**New module `tests/mesh/test_birdcage_implant_rod_tip.py`.** It imports by name and never a `test_*` symbol:
- `ROD`, `ROD_RADIUS`, `ROD_LENGTH`, `ROD_CENTRE`, `_build_rod`, `_rod_volume`, `_rod_area`, `_phantom_minus_rod` and `_facet_check` (`test_birdcage_implant_rod.py`);
- `_build` (`test_birdcage_port_sheets`);
- `CAD_ANALYTIC_BAND` (`test_birdcage_conductor_hole.py:123`), `EXACT` (`:117`) and `_owned_census` (`:126`);
- `TERMINAL_AREA_BAND` (`test_birdcage_leg_gaps.py:73`), `CELL_COUNT_BAND` (`test_birdcage_port_sheet_prerequisite.py:57`) and `_exterior_facet_area`.

T1 comes from `scripts/probes/imp3_step4_palace_fields.py:41–49`, loaded by importlib as `test_imp3_palace_fields.py:55–62` does it. The module registers `ROD_1G_RUNG_CELL_RECORD = 252_118`, version-tagged and citing both logs. The env knob `FEM_EM_GEO36_STEP4` ∈ {`noop`, `rungs`} selects the builds; the refusals run ungated.

**Anchors (asserted, bands imported).**
- (R) `tip_size` ∈ {0, −1e-4, nan, inf, r/3, `"fine"`, `True`}: each raises `ValueError` matching `implant_rod`, with `gmsh.isInitialized()` False (the K7 pattern, `:575–599`).
- (N0) In `noop`, the unset build and `tip_size=ROD_RADIUS / 4.0` give **integer-equal** cells, owned censuses and tag-402 counts. The unset build against `ROD_1G_RUNG_CELL_RECORD` is asserted within `CELL_COUNT_BAND`, with the ratio printed to 6 d.p. (predicted 1.000000) and the integer printed.
- (A0) Every rod build in a window has `cad_mass_by_group`, `wall_area_m2` and `wall_centroid_m` bitwise equal to the unset build's, because sizing cannot move CAD. `tip_size_m` equals the request bitwise, or `ROD_RADIUS / 4.0` when unset.
- (A1–A6) Per rung, against the window's own off build. A5's flag-off record line becomes the off count, printed beside 201 072.
- (D) Per build, N1curl `functionspace(msh, ("N1curl", k))`, `size_global * bs`: degree 2 = 2E + 2F and degree 1 = E, exactly. This is item 284's identity (`IMP-3.md:948`), and it is scalar-agnostic.
- (G) N_ball is the number of tag-3 cells whose centroid lies within `tip_ball_radius_m` of either cap centre. It must grow strictly: N_ball(r/(4√2)) > N_ball(unset), and N_ball(r/8) > N_ball(r/(4√2)).
- Each rung stays below 1 000 000 cells, checked before (D) as a runaway guard.

**Printed, never gated (rule (e): no prior measurement).**
- Predicted N_ball growth (r/4 ÷ h)³: 2.83 at r/(4√2) and 8 at r/8.
- Mean ball-cell edge against the requested `VIn`.
- A4's meshed wall ratio, predicted to rise toward 1.
- Mean edge of the cells within 0.5 mm of T1, predicted **unchanged**: T1 (z = 0.017531) sits 2.531 mm from the cap centre, outside the 2 mm ball.
- Cells against `ROD_CELL_GUARD_W2` 320 000 (`test_imp3_rod_birdcage.py:115`).
- Predicted cells: ≈ 298 k (≤ 345 k) at r/(4√2) and ≈ 427 k (≤ 609 k) at r/8.
- **I15 price per rung.** Sweep = 284.88 s × (N₂ / 1 620 748)^e at e = 1.23 and 1.37 (`IMP-3.md:924,948`). Window = 450 s + (sweep − 284.88) + (build elapsed − 64 s). RSS = 50.97 GiB × the same ratio^e. A window over 590 s is flagged "over the heavy ceiling at `-n 8`". Predicted: ≈ 1.9 M dofs and a ≈ 530 s window at r/(4√2); ≈ 2.75 M dofs and a ≈ 760 s window at r/8. Only the rod arm changes per rung, because the off mesh carries no rod.

**Controls.** (R) and (G) are the asserted controls. A knob that never reaches gmsh leaves N_ball flat (floor 1×, no ceiling). A knob that leaks past `VIn` breaks (N0) by a nonzero integer.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; no `-k`; no pipes). gmsh builds serially on rank 0 (`io/mesh.py:3998,4480`) and broadcasts a build error to all ranks (`:4513–4525`). Each window builds the off mesh first and the rungs from coarse to fine (`GEO-23`).
- **W1 `noop`** builds off, unset and r/4. Predicted ≈ 240 s, range 200–330 (basis: 181 s for three builds, `GEO-36.md:107`; off 43.59 s and rod ≈ 64 s, `IMP-3.md:915,924`).
  `scripts/testing/run_and_log.sh GEO-36-step4-W1-noop "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_GEO36_STEP4=noop timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod_tip.py -v -s --tb=short'"`
- **W2 `rungs`** builds off, unset, r/(4√2) and r/8. Predicted ≈ 365 s, range 300–480. It is the closing window for the new module. Same command with `FEM_EM_GEO36_STEP4=rungs`, `timeout -k 30 590` and label `GEO-36-step4-W2-rungs`.
- **W3, the rule (c) re-gate.** Predicted ≈ 280 s, range 240–340 (`GEO-36.md:107–108`). It must print 167 371 / 182 371 / 116 085 unchanged. A `test_birdcage_volumes_partition_the_box` red carrying the surface 59/79 string is the 🔴 OPEN known-issue (`GEO-36.md:97`); step 3's W3 read none.
  `scripts/testing/run_and_log.sh GEO-36-step4-W3-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod.py tests/mesh/test_birdcage_port_tags.py tests/mesh/test_birdcage_conductor_hole.py -v -s --tb=short'"`
- No window passes ≈ 10 min, so no durable capture is needed. Real build: no `dolfinx-complex-mode` and no `FEM_EM_REQUIRE_COMPLEX`. Meshing and dof counts are scalar-agnostic, and step 2 ran real (`GEO-36.md:25`).

**Traps.**
- `cell_tags.values`, `assemble_scalar` and the ball / edge sums are rank-local. Reduce them over owned cells only, with no collective inside a `rank == 0` branch.
- Rule (i): no edit to the new module after W2 starts. W3 runs on the committed `src/`.
- Status 124: run `pgrep -c python3` in the container first, then the rubric's cache sweep only if it reads 0.
- A gmsh FAIL is recorded verbatim with its surface pairs and never retried.
- Never import `tests/validation/test_imp3_*` in the real build.

**Scope.** Not in this item:
- an EM solve, η_rod or SAR;
- a change to the off mesh;
- re-asserting the strips (the legs sit at r = 0.07 m, `io/mesh.py:5870`; I15's solve item re-asserts G2);
- a ball-radius knob.

`GEO-36` stays ✅ and `IMP-3` stays ✅. I15's solve windows (heavy or `xl`, `-n`, the rung set) stay the weekly's to commission.

**Negative result.** A red on (R), (N0), (A0)–(A6), (D) or (G), a gmsh FAIL, or the guard tripping: report every number, park on `attempt/GEO-36-step4-<ts>`, open a known-issues entry and mark BLOCKED in the same commit. Never loosen, never retry. If W3 moves a record, revert `src/` and stop. No `PREDICTED RED` line.

**Status it can move.** The `GEO-36` row gains a step-4 record: the knob, and per rung the cells, N_ball, mean ball edge, E / F, degree-1 / degree-2 dofs and build time. The outcome carries the I15 price table, which the review relays to the weekly. The slot makes no §10 edit.

**Review rulings (2026-10-07 10:30 daily review, item 298).** (1) **Identities A1–A6, not K1–K6:** accepted. I15's fixture is the hole build at the 1 g rung (`test_imp3_tip_ratio_1g.py:121–128`). (2) **Record equality:** the unset build equals an explicit `tip_size = r/4` build as integers in the same window. The 252 118 record is asserted at `CELL_COUNT_BAND` and the integer is printed. Accepted; this is project practice for mesher-dependent records. (3) **Rungs r/(4√2) and r/8**, not r/8 and r/16: accepted, because r/16 is beyond any tier. (4) **T1 sits outside the cap ball** (2.531 mm from the cap centre against a 2 mm Ball radius), so a `VIn`-only knob does not refine T1's own cell. The brief prints the edge length near T1 per rung. **Whether I15 also needs a ball-radius knob is the weekly's question**, and this item does not add one. (5) **Price:** the drafter's scaling puts an r/8 rod-arm degree-2 window at ≈ 760 s at `-n 8`, over the heavy ceiling. W2's measured cells and dofs set the figure the weekly uses for `xl`. (6) Class: chain step (solve), 10 measured minutes. Three windows (≈ 15 min of compute) is above the median; the item still fits one slot.

**Outcome 2026-10-07 (12:00 CDT slot, implementer): DONE — W1, W2 and W3 green; no red, no gmsh failure.**
- **Code.** `_normalise_implant_rod` returns `"tip_size"` (`None` when absent; refuses bool, non-number, non-finite, ≤ 0 and > r/4 with a `ValueError` naming `implant_rod`, before gmsh). Both cap Balls take `VIn` = the value (`r/4` when `None`); radius 2r, Threshold unchanged; the sizing print names the value. `implant_rod_cad` gains `tip_size_m` and `tip_ball_radius_m`. Additive only (`io/mesh.py`). New `tests/mesh/test_birdcage_implant_rod_tip.py` (`ROD_1G_RUNG_CELL_RECORD = 252 118`, dolfinx 0.11.0.post0 / gmsh 4.15.2).
- **Windows** (all Status 0, `-n 2`, real):
  - W1 `noop` `20261007T170256Z_GEO-36-step4-W1-noop.log`: 19 passed / 5 skipped, 238 s (predicted ≈ 240).
  - W2 `rungs` `20261007T170700Z_GEO-36-step4-W2-rungs.log` (closing window): 19 passed / 5 skipped, 330 s (predicted ≈ 365, range 300–480). Module untouched after W1.
  - W3 re-gate `20261007T171236Z_GEO-36-step4-W3-regate.log`: 39 passed / 4 skipped, 257 s; 167 371 / 182 371 / 116 085 print unchanged; no surface-59/79 red.
- **Asserted.** (R) 7/7 refused in 0.000 s, gmsh not initialised. (N0) unset = explicit r/4: 252 118 cells, identical census, 4 958 tag-402 facets; unset / 252 118 = 1.000000. (A0) CAD masses, wall area 2.0106192982974677e-4 m² and centroid bitwise equal at every rung; `tip_size_m` = request bitwise. (A1) 1.110e-16, (A2) 1.110e-16, (A3) 0, (A6) \|x̄ − 0.015\| 1.7e-18 at every rung. (A4) 0.996254 / 0.996674 / 0.996934 (unset / r/(4√2) / r/8; rises toward 1 as predicted). (A5) tag sets equal, census sums exact, 402 all exterior on tag 3, off 0. (D) exact on off and every rung. (G) N_ball 11 938 → 32 826 → 93 132, strictly increasing.
- **Per rung (W2; unset / r/(4√2) / r/8).**

  | rung | cells | N_ball (÷ unset; pred.) | mean ball edge [m] (÷ VIn) | 402 facets | E | F | d1 dofs | d2 dofs | build s (mesh s) |
  |---|---|---|---|---|---|---|---|---|---|
  | unset (r/4) | 252 118 | 11 938 | 3.398e-4 (1.359) | 4 958 | 302 002 | 508 372 | 302 002 | 1 620 748 | 61.65 (57.33) |
  | r/(4√2) | 276 088 | 32 826 (2.750; 2.83) | 2.433e-4 (1.377) | 6 256 | 331 014 | 556 963 | 331 014 | 1 775 954 | 71.40 (66.55) |
  | r/8 | 342 072 | 93 132 (7.801; 8) | 1.722e-4 (1.378) | 8 790 | 410 054 | 690 196 | 410 054 | 2 200 500 | 86.71 (80.86) |

  Off build: 201 072 cells (= the printed 201 072), E 238 627, F 403 801, d2 1 284 856, 41.1 s.
- **Printed readings.** Cells came in under the predictions (≈ 298 k / ≈ 427 k): +9.5 % and +35.7 % over unset. **r/8's 342 072 cells exceed `ROD_CELL_GUARD_W2` 320 000** — an `IMP-3`-style (G1) guard would trip at r/8. Near-T1 (0.5 mm) mean edge 6.98e-4 / 6.42e-4 / 5.08e-4 m on 11 / 19 / 26 cells — predicted unchanged, it moved −8 % / −27 %: the Ball's `VOut`/Min interplay grades past its 2 mm radius, so T1's cells refine too (on a handful of cells; printed only).
- **I15 price (rod arm, degree 2, `-n 8` basis item 284 W3: sweep 284.88 s at 1 620 748 dofs, RSS 50.97 GiB; window = 450 + Δsweep + (build − 64)).**

  | rung | d2 dofs | e | sweep s | window s | RSS GiB | flag |
  |---|---|---|---|---|---|---|
  | unset | 1 620 748 | 1.23 / 1.37 | 284.9 | 447.7 | 50.97 | — |
  | r/(4√2) | 1 775 954 | 1.23 / 1.37 | 318.8 / 322.9 | 491.3 / 495.4 | 57.04 / 57.77 | — |
  | r/8 | 2 200 500 | 1.23 / 1.37 | 415.0 / 433.1 | 602.8 / 621.0 | 74.24 / 77.49 | over the heavy ceiling at `-n 8` |

  r/16 (slot arithmetic from W2, not a logged print): from r/8's growth the ball holds ≈ 81 k of r/8's +90 k cells; a further 8× on the ball puts r/16 near ≈ 1 M cells / ≈ 6.5 M d2 dofs — still no tier (the brief's argument stands).

### §9 item 316 — 2026-10-08 (chain I18a — the coated rod in the phantom; no solve)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 I18a (`PROJECT_PLAN.md:4694–4698`), weekly-enumerated, serial on nothing; §9 rule (4). **No solve; no band, tolerance, record or Done-when moves.**

- **Executor:** implementer. **Tier:** heavy (each window ≈ 190–260 s), real build, `-n 1` then `-n 2`. **Cost:** ≈ 11 measured (chain no-solve median, `WF-9.md:155`; step 4 took ≈ 10) / ≈ 50 predicted min. **Depends on:** nothing.
- **Geometry (the weekly's volume read literally).** `IMP-5`'s sleeve a = `ROD_RADIUS` 1 mm, b = 1.5 mm, εᵣ 3 / σ 0 (`IMP-5.md:26–29`) on step 3's kept rod (`_build_kept_rod`, `tests/mesh/test_birdcage_implant_rod.py:154`; `ROD` `:93–97`). The weekly's coat volume π(b² − a²)L is a **shaft-only sleeve**: the cylinder of radius b between the two cap-centre planes, minus the rod. The hemispherical caps stay **bare** (an insulated lead with exposed tips). The sleeve's end annuli meet the cap spheres at their equators, where the sphere leaves the plane at a right angle, so there is no tangency sliver. Volumes: V_rod = πa²L + 4πa³/3 = 9.8437e-8 m³; V_coat = π(b² − a²)L = 1.17810e-7 m³; gel = πR²H − V_rod − V_coat = 2.259784e-4 m³.
- **`src/` change (one additive key; rule (c), disclosed).** `implant_rod["coat_radius"]` (metres), mirroring `IMP-5`'s `coat_radius` on `cylindrical_domain` (`io/mesh.py:71,100–134`).
  - `_normalise_implant_rod` (`io/mesh.py:4589–4672`) returns `"coat_radius"` (`None` when absent). Before gmsh it raises a `ValueError` naming `implant_rod` on a bool, non-number, non-finite value, a value ≤ `radius`, `keep_cells` not `True` (a coat around a hole is out of scope), or a sleeve not strictly inside the phantom (the existing inside-check with b in place of r). It runs twice (`:4257`, `:4690`); `None` must pass the second call unchanged.
  - Build (`:5102–5152`): add `cyl_b` = `addCylinder(end0, L·axis, b)`. Cut the phantom and the air box by `cyl_b` as well as the rod (`removeTool=False`). Make the sleeve with `occ.cut([cyl_b], rod_dimtags, removeObject=True, removeTool=False)`, raise unless exactly one solid remains, and append its tag to the fragment tools beside `rod_tag` (`:5197`). Its pieces take the new cell group **`BIRDCAGE_IMPLANT_COAT_CELL_TAG = 404`**, named `"implant_coat"`, beside `:280`, `:5289–5292`, `:5543` and `:5599`. Grep `src/` for 404 first; it read empty on 2026-10-08.
  - Interfaces: 402 stays the (3, 403) interface (`:4547–4555`), so it now covers only the bare caps. No new facet tag is emitted; the interface areas are measured by tag-pair `dS` in the test. The dolfinx-side rebuild must still find every (3, 403) facet, and the two-use check (`:5456–5476`) must count the sleeve as a neighbour, not as a hole.
  - Sizing (`:5931–5960`) is **unchanged**. The Threshold holds r/3 out to r, so the 0.5 mm sleeve meshes at ≈ 1.5 cells across. No sizing knob.
  - `implant_rod_cad` (`:5517–5528`) gains `coat_radius_m`, `coat_cell_tag` and `coat_kept_mass_m3` (additive; `None` when unset).
- **New module `tests/mesh/test_birdcage_implant_rod_coat.py`.** It imports by name and never restates a band:
  - `ROD`, `ROD_RADIUS`, `ROD_LENGTH`, `ROD_CENTRE`, `_build_kept_rod`, `_rod_volume`, `_rod_area`, `KEPT_ROD_CELL_RECORD`, `KEPT_ROD_CELL_GUARD` and `_interior_facet_check` (`test_birdcage_implant_rod.py`);
  - `_build` (`test_birdcage_port_sheets`);
  - `EXACT` and `_owned_census` (`test_birdcage_conductor_hole.py:117,126`);
  - `TERMINAL_AREA_BAND` (`test_birdcage_leg_gaps.py:73`);
  - **`STRIP_BAND`** (`tests/validation/test_ops61_exact_strip_ports.py:58`).

  It builds the flag-off, kept and coated meshes in that order (`GEO-23`).
- **Anchors (asserted).** Every relative error is taken against the comparand.
  - (P1) CAD(403) = V_rod at `STRIP_BAND`.
  - (P2) CAD(404) = π(b² − a²)L at `STRIP_BAND`.
  - (P3) CAD(3) = πR²H − V_rod − V_coat at `STRIP_BAND`.
  - (P4) Σ CAD(coated) = Σ CAD(flag off) at `STRIP_BAND`, so nothing is removed. Every group outside {3, 403, 404} is bitwise equal to the kept build's.
  - (P5) OCC boundary area of 403 = 2πaL + 4πa², and of 404 = 2π(a + b)L + 2π(b² − a²), both at `STRIP_BAND`.
  - (P6) Census: tag set = the kept set ∪ {404}. The owned census sums to `size_global`, 404 has > 0 cells, and 404 is absent in the kept and off builds. Every 402 facet lies between one tag-3 cell and one tag-403 cell (`_interior_facet_check`). Cells stay under `KEPT_ROD_CELL_GUARD` 320 000.
  - (P7) The OCC coat centroid sits at (0.015, 0, 0) to `EXACT`·ℓ, the K6 pattern.
  - (P8) The refusals above raise with `gmsh.isInitialized()` False (the K7 pattern, `:575–599`).
  - (P9) The kept build within the coated window reads `KEPT_ROD_CELL_RECORD` 182 371 at ratio 1.000000, printed with the integer.
- **Printed, never gated (rule (e): no prior measurement).**
  - Cells per tag (predicted: total 185 k–215 k; 403 within ±10 % of 15 000).
  - Meshed volume ÷ CAD for 403 and 404, beside step 3's 0.986667 for the bare rod.
  - Meshed `dS` areas ÷ closed form: (403 | 404) against 2πaL, (404 | 3) against 2πbL + 2π(b² − a²), (403 | 3) against 4πa². Each is predicted inside `TERMINAL_AREA_BAND`, and none is asserted.
  - Min and mean edge in 404, E / F, and degree-1 / degree-2 N1curl dofs (= E, = 2E + 2F; the (D) identity of step 4, asserted exactly). Mesh wall time. These figures price item 317.
- **Negative controls (rule (e)).**
  - (C1c, **asserted**) The kept (uncoated) build misses P3 by V_coat / gel = 1.17810e-7 / 2.259784e-4 = 5.213e-4, i.e. > `STRIP_BAND` by ≈ 5.2e5×. That is also the ceiling: the coat's whole volume fraction is the most this fixture can show.
  - (C2c, **predicted**, printed) A capped sleeve would add 4π(b³ − a³)/3 = 9.948e-9 m³ and miss P2 by 8.44e-2. This control tells the two sleeve readings apart.
  - (C3c, **predicted**, printed) Swapping 403 / 404 misses P1 by \|V_coat − V_rod\| / V_rod = 0.1968.
- **Windows** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; no `-k`; no pipes; real build). Basis: step 3 ran three builds in W1 188 s (`-n 1`) and W2 181 s (`-n 2`), and step 4's W3 took 257 s.
  - **W0 (smoke):** the `STRIP_BAND` import must collect in the real build, because its module pulls in `tests/validation` modules.
    `scripts/testing/run_and_log.sh GEO-36-step5-W0-collect "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/mesh/test_birdcage_implant_rod_coat.py --collect-only -q'"`
    If collection fails, **stop and ask the review**. Never retype 1e-9.
  - **W1** `-n 1`, predicted ≈ 200 s, range 170–300:
    `scripts/testing/run_and_log.sh GEO-36-step5-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 420 mpiexec -n 1 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod_coat.py -s -v --tb=short'"`
  - **W2** `-n 2`, the window of record: the same command with label `GEO-36-step5-W2` and `-n 2`.
  - **W3, the rule (c) re-gate**, predicted ≈ 260 s. It must print 167 371 / 182 371 / 116 085 unchanged.
    `scripts/testing/run_and_log.sh GEO-36-step5-W3-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod.py tests/mesh/test_birdcage_port_tags.py tests/mesh/test_birdcage_conductor_hole.py -v -s --tb=short'"`
    A `test_birdcage_volumes_partition_the_box` red carrying the surface 59 / 79 string is the 🔴 OPEN known-issue (`GEO-36.md:97`).
- **Traps already paid for.**
  - Interior dim-2 gmsh groups hang `model_to_mesh` at `-n 2` (known-issues 9). Never emit one; any interface is rebuilt dolfinx-side (step 3's deviation, `:103`).
  - Read walls from the meshed boundary, never from the tool's faces (`TH-15` step 0).
  - A failed gmsh build contaminates its process (`GEO-23`): build off → kept → coated, and never retry a FAIL.
  - `cell_tags.values` is rank-local: use `_owned_census` and the ghost-aware `_cell_values_with_ghosts`.
  - Rule (i): no edit to the new module after W1 starts. W3 runs on the committed `src/`.
  - Status 124: run `pgrep -c python3` in the container before anything else.
  - Never import `tests/validation/test_imp3_*` or `test_wf9_*` in the real build.
- **Scope.** Not in this item:
  - an EM solve, a material assignment, a TF or SAR (all of that is item 317);
  - a capped coat;
  - a coat on the hole build;
  - a coat-sizing knob;
  - re-asserting the strips;
  - the `IMP-5` closed form on this fixture (the rod is finite and the series is 2-D).

  `GEO-36` stays ✅.
- **Negative result.** A red on P1–P9 or (C1c), a gmsh FAIL on the sleeve (recorded verbatim with its surface pairs), or the guard tripping: report every number, park on `attempt/GEO-36-step5-<ts>`, open a known-issues entry and mark item 316 BLOCKED in the same commit, which makes item 317 skip. Never loosen, never retry. If W3 moves a record, revert `src/` and stop. No `PREDICTED RED` line.
- **Status it can move.** The §7 `GEO-36` row (Done-when: "the identities asserted in one logged window at `-n 2`, cells and time recorded", `PROJECT_PLAN.md:1860`) gains a step-5 record: the coat key, tag 404, P1–P7 values, cells per tag and the build time. It stays ✅. I18b (item 317) becomes runnable. No §10 edit.

**Review rulings (2026-10-08 03:00 daily review, item 316; with item 317 in `WF-9.md`).** (1) The sleeve is read literally — π(b² − a²)L over the shaft, caps bare (an insulated lead with exposed tips); (C2c) stays printed. (4) `STRIP_BAND` is imported by name; W0 `--collect-only` stops the slot on an import failure — never retype 1e-9. (5) The three interface areas stay printed against `TERMINAL_AREA_BAND`, as drafted; the weekly's target is the volume partition, not the areas. The `src/` change (`coat_radius`, tag 404) is additive and re-gated per rule (c) on the step-4 W3 set, as drafted.

**Outcome 2026-10-08 (09:00 CDT slot, take-next, implementer): DONE — W0–W3 green first run; no red, no gmsh failure.**
- **Code.** `_normalise_implant_rod` returns `"coat_radius"` (`None` when absent; refuses bool, non-number, non-finite, ≤ radius, `keep_cells` not `True`, sleeve not strictly inside the phantom, all with a `ValueError` naming `implant_rod`). Build: `cyl_b` cut from phantom and air box beside the rod (tools kept), sleeve = `cut(cyl_b, rod)` checked to be one solid, entered once into the fragment; its pieces take `BIRDCAGE_IMPLANT_COAT_CELL_TAG = 404` `"implant_coat"` (piece policy after 403, before the phantom). `implant_rod_cad` gains `coat_radius_m`, `coat_cell_tag`, `coat_kept_mass_m3`, `coat_occ_mass_m3`, `coat_wall_area_m2`, `coat_centroid_m` (all `None` when unset). Sizing unchanged; 402 stays the (3, 403) interface. New `tests/mesh/test_birdcage_implant_rod_coat.py` (19 tests).
- **Windows** (all Status 0, real build):
  - W0 collect `20261008T141318Z_GEO-36-step5-W0-collect.log`: 19 collected, 3 s.
  - W1 `-n 1` `20261008T141327Z_GEO-36-step5-W1.log`: 27 passed / 4 skipped, 202 s (predicted ≈ 200).
  - W2 `-n 2` `20261008T141656Z_GEO-36-step5-W2.log` (window of record): 27 passed / 4 skipped, 200 s; every printed number identical to W1. Module untouched after W1.
  - W3 re-gate `20261008T142024Z_GEO-36-step5-W3-regate.log`: 39 passed / 4 skipped, 254 s; 167 371 / 182 371 / 116 085 unchanged; the surface 59 / 79 string appears only as the caught `GEO-1` control.
- **Asserted.** (P1) CAD(403) rel 0; (P2) CAD(404) 1.178097245096e-07 rel 2.2e-16; (P3) CAD(3) 2.259784247641e-04 rel 2.2e-16; (P4) Σ rel 0, 10/10 groups outside {3, 403, 404} bitwise equal to the kept build's; (P5) OCC area 403 rel 0, 404 4.790928796724e-04 rel 1.1e-16; (P6) tag set = kept ∪ {404}, census 190 747/190 747, 404 absent in kept / off, 534 owned 402 facets all between one tag-3 and one 403 cell, 190 747 < 320 000; (P7) \|x̄ − 0.015\| 1.7e-18, y, z ≤ 8.2e-19; (P8) 9/9 refused in 0.000 s, gmsh not initialised; (P9) kept 182 371 = 1.000000; (C1c) 5.213317e-04 (predicted 5.213317e-04) = 5.2e5× `STRIP_BAND`; (D) d1 226 533 = E, d2 1 219 372 = 2E + 2F.
- **Printed.** Off 116 085 cells; kept 182 371; coated 190 747 (predicted 185 k–215 k) — census 3: 41 602, 403: 15 035 (predicted 15 000 ± 10 %), 404: 18 495; mesh 48.3 s (kept 47.3 s). Meshed / CAD: 403 0.986667 (= step 3's bare rod), 404 0.999821. Meshed dS / closed form: (403 \| 404) 0.996779, (404 \| 3) 0.998565, (403 \| 3) 0.988379 — all inside `TERMINAL_AREA_BAND`, as predicted. 404 edges min 1.189e-4 m, mean 4.087e-4 m. E 226 533, F 383 153. (C2c) 8.4444e-2 (predicted 8.44e-2); (C3c) 0.1968 (predicted 0.1968).
- **Item 317 price input.** Coated mesh d1 226 533 / d2 1 219 372 dofs (kept-rod 182 371-cell build for comparison: step 3).

## §7 row moved verbatim 2026-10-10 (rotate_plan_archive.py chunks)

 A rod in the phantom: the implant-in-phantom generator (Phase 7 (7b)) | ✅ Audit PASS 2026-10-02; the `implant_rod=` generator's CAD identities are asserted at `-n 2` (step 2, §9 item 171). Step 3 (2026-10-04, §9 item 219): additive `keep_cells` keeps the rod's cells (tag 403, wall 402 interior), K1–K7 green at `-n 1` / `-n 2`, 182 371 cells, `20261004T200730Z_GEO-36-step3-W2.log`, heavy, 181 s. **Step 4 (2026-10-07, §9 item 298, I15 part 1, no solve):** additive `implant_rod["tip_size"]` [m] = the cap-Ball `VIn` (unset → r/4, refused before gmsh outside (0, r/4]); on the 1 g-rung hole build at `-n 2`, unset = explicit r/4 integer-for-integer and = the 252 118 record (1.000000), CAD bitwise invariant, A1–A6 and N1curl d2 = 2E + 2F / d1 = E exact at every rung. Per rung unset / r/(4√2) / r/8: cells 252 118 / 276 088 / 342 072, N_ball 11 938 / 32 826 / 93 132 (2.75× / 7.80×, predicted 2.83 / 8), mean ball edge 3.398e-4 / 2.433e-4 / 1.722e-4 m, E 302 002 / 331 014 / 410 054, F 508 372 / 556 963 / 690 196, d2 dofs 1 620 748 / 1 775 954 / 2 200 500, build 61.7 / 71.4 / 86.7 s; r/8 is over 320 000 cells and its priced degree-2 window (603–621 s at `-n 8`) is over the heavy ceiling (`20261007T170700Z_GEO-36-step4-W2-rungs.log`, heavy, 330 s; W1 noop 238 s, W3 re-gate 257 s, records unchanged). **Step 5 (2026-10-08, §9 item 316, I18a, no solve):** additive `implant_rod["coat_radius"]` [m] (needs `keep_cells`; refused before gmsh otherwise) adds `IMP-5`'s shaft-only sleeve (a = 1 mm, b = 1.5 mm, caps bare) as cell tag 404 `"implant_coat"`; P1–P7 at `STRIP_BAND` read 0 / 2.2e-16 / 2.2e-16 / 0 / ≤ 1.1e-16 / centroid 1.7e-18 m, (C1c) the uncoated build misses P3 by 5.213e-4; 190 747 cells (3: 41 602, 403: 15 035, 404: 18 495), E 226 533, F 383 153, d2 1 219 372, mesh 48.3 s (`20261008T141656Z_GEO-36-step5-W2.log`, heavy, 200 s; W3 re-gate 254 s, records unchanged). **Done-when:** the identities asserted in one logged window at `-n 2`, cells and time recorded ⇒ ✅. No solve. *History: `docs/planning/chunks/GEO-36.md`.* 

### §9 item 372 — 2026-10-11 (chain I21a — the kept rod in the F-human phantom; chain step (no solve))

**Outcome (2026-10-11 13:30 CDT slot): ✅ record, no glyph move.** The narrow `src/` change landed (rule (c), disclosed): `_normalise_implant_rod` lifts the `as_hole` refusal only for `keep_cells: True` without `coat_radius`; the `as_hole` fragment enters the kept rod once. W0 collect 14 tests (`20261011T185749Z_GEO-36-step6-W0-collect.log`); W1 `-n 1` 14 passed, 184 s (`20261011T185757Z_GEO-36-step6-W1.log`); W2 `-n 2`, window of record, 14 passed, 180 s, every printed number equal to W1's (`20261011T190108Z_GEO-36-step6-W2.log`). Readings: 436 541 cells = 369 387 + 67 154 (ratio 1.181798, inside predicted [420 k, 455 k]), tag 403 15 025, tag 3 55 090; (F1)/(F2)/(F3)/(F4 OCC) rel 0.000; meshed 402 / A_rod 0.996276; meshed 403 / V_rod 0.986763; (F5) 240 surfaces, area rel 1.501e-10 (C2h swallow miss 9.598e-4); (F6) 4 984 owned 402 facets, all interior tag-3/403; 401 facets 78 052 (rod-free 78 066, printed); (C1h) 4.4230e-5 = 4.42e4× band; (C3h) 2.26e4; E 579 178, F 918 745, degree-2 dofs 2 995 846; mesh wall 166.17 s, build 176.09 s at `-n 2`. `F_HUMAN_KEPT_ROD_CELL_RECORD = 436541` registered from W1. Re-gates: W3a 39 passed / 4 skipped, 116 085 / 167 371 / 182 371 unchanged, (A7) as_hole still refused (`20261011T190414Z_GEO-36-step6-W3a-regate.log`, 257 s); W3b 13 passed, (A0) 369 387 and the 401 line (240, 2.094760146e-01) unchanged (`20261011T190838Z_GEO-36-step6-W3b-regate.log`, 191 s). Item 381 is runnable.

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 I21a (`PROJECT_PLAN.md:4747–4751`, the 2026-10-10 weekly), lane 3 / Phase 7 (7b), serial on nothing; §9 rule (4). **Target: CAD identities, the weekly's text. No solve; no band, tolerance, record or Done-when moves; every band is imported.** A record on a ✅ row.

- **Executor:** implementer. **Tier:** standard by expectation, heavy by ceiling: one F-human hole build per window takes ≈ 115 s, so a window is ≈ 150–200 s, which can pass the 180 s standard ceiling; `timeout -k 30 420`. Complex build (the imported F-human helpers live in `tests/validation`; the mesh is scalar-agnostic). **Ranks:** `-n 1`, then `-n 2` (the window of record). **Cost:** ≈ M measured (chain no-solve class) / ≈ 50 predicted min (15 fixed + 15 code + ≈ 14 compute across five windows + spare). **Depends on:** nothing.
- **Cost basis (measured).** The F-human hole build at `-n 2`: mesh wall 105.01 s, build 115.7 s, window 192 s (`20261010T005723Z_GEO-35-s2a-W1-wf12-buildonly.log:6219–6224, :6358`). Keeping the rod on F-small added ≈ 5 s of mesh time (hole 42 s → kept 47.3 s; step 2 / step 5 outcomes above).

**Why a `src/` change is needed (read before coding).** `_normalise_implant_rod` refuses **any** rod with `as_hole` (`src/fem_em_solver/io/mesh.py:4614–4619`: "the coil-cavity rule would swallow the rod wall"). The `as_hole` fragment (`:5258–5264`) also does not enter the kept rod as a tool; only the solid route does (`:5267–5270`). So `implant_rod` does not compose with `FHumanBirdcageGeometry(as_hole=True)` / `_build_f_human(comm, as_hole=True)` today. The weekly's text — "`keep_cells` keeping the rod" — names the composition that is safe. With the rod's cells kept, every rod-wall face bounds two meshed volumes (gel and implant). The cavity rule (`:5452–5462`, single-use non-wall faces) therefore cannot take the wall into 401.

**`src/` change (rule (c), additive, disclosed).**
- `_normalise_implant_rod`: the `as_hole` refusal fires only when `keep_cells` is not `True` **or** `coat_radius` is set. The hole rod and the coat stay refused under `as_hole`; neither is gated. `symmetry_quadrant` stays refused always. The message keeps the word `implant_rod`, because the step-2 (A7) "as_hole" case (`tests/mesh/test_birdcage_implant_rod.py:578`) builds the **hole** rod and must still raise before gmsh.
- `_build_birdcage_port_model`, `as_hole` branch: when `keep_rod_cells`, append `(3, rod_tag)` to `tool_dimtags` (`:5258`). Enter it once, as the solid route does; entering it twice is the `GEO-9` overlapping-facet failure. Nothing else changes, because each of these already keys on `keep_cells`:
  - the pre-cut (`:5187–5198`) removes the rod from the phantom and the air box with the tool kept;
  - the piece policy (`:5365`) maps rod descendants to 403;
  - the dolfinx-side 402 rebuild (`:4555–4563`) and `expected_groups` (`:5717–5721`).
- Docstrings at `:4029–4047` and the 403 comment (`:274–280`) gain one sentence each. Grep `src/` for any other `as_hole`-and-implant guard first.

**Fixture.** The H9′ mesh's generator kwargs plus the rod, called directly, with no helper edit: `kwargs = _f_human_kwargs(as_hole=True); kwargs["implant_rod"] = ROD_FH; MeshGenerator.birdcage_port_domain(comm=comm, return_diagnostics=True, **kwargs)`, where `ROD_FH = dict(ROD, keep_cells=True)`.
- `ROD` (`tests/mesh/test_birdcage_implant_rod.py:93–97`): r = 1 mm, ℓ = 30 mm between cap centres, axis ẑ, centre (0.015, 0, 0).
- `_f_human_kwargs` is pure (`tests/validation/test_wf7_f_human_b1_quadrature.py:67–98`). R = 64.2857 mm, H = 171.4286 mm and `resolution` 0.015 m come from its dict (`_params(F_HUMAN_RING_RADIUS, scale_sizing=False)`) and are never retyped.
- The rod is strictly inside: cap tips at \|z\| = 16 mm ≪ 85.7 mm, radial ≤ 16 mm ≪ 64.3 mm.
- **Placement — a review point, not a ruling.** The weekly gave no position. This brief drafts the reuse of `ROD`'s centre. That keeps `TIPS`, `_ball_mean_closed_form`'s `ROD_CENTRE` pin (`tests/validation/test_imp3_rod_birdcage.py:140, :229–272`) and `WF-9`'s masks valid for I21b. It also makes the F-human tip ratio comparable to F-small's records. The alternative is the rod at 359's hotspot c\* (r = 45 mm, az₁₇ = 11.25°, z_c). That maximises "the rod's effect on SAR10g(c\*)" but needs a re-centred copy of the closed-form ball mean in I21b. The review rules.

**New module `tests/mesh/test_birdcage_implant_rod_f_human.py`.** It imports by name and never a `test_*` symbol:
- from `tests/mesh/test_birdcage_implant_rod.py`: `ROD`, `ROD_RADIUS`, `ROD_LENGTH`, `ROD_CENTRE`, `_rod_volume`, `_rod_area`, `_interior_facet_check`, `_cell_values_with_ghosts`. **Not** `_phantom_minus_rod`, which uses F-small's R and H;
- `_f_human_kwargs` (`test_wf7_f_human_b1_quadrature.py:67`); `_ring_ports` (`test_birdcage_f_human_rung.py:134`);
- `HOLE_CELL_RECORD` 369 387 (`tests/validation/test_wf12_f_human_copper_matrix.py:55`); `_census_checks` (`tests/validation/test_wf12_f_human_hole.py:218`); `_owned_census(msh, tags, dim, comm)` and `EXACT` (`tests/mesh/test_birdcage_conductor_hole.py:126, :117`);
- `STRIP_BAND` (`tests/validation/test_ops61_exact_strip_ports.py:58`); `TERMINAL_AREA_BAND` (`tests/mesh/test_birdcage_leg_gaps.py:73`); `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`); `_outer_box_tags` (`tests/validation/test_th14_birdcage_copper.py:133`).

**Constants the module defines.**
- **Records registered here.** These are first readings, version-tagged, printed and not asserted in this item:
  - `F_HUMAN_KEPT_ROD_CELL_RECORD`: the W1 count, printed again by W2. I21b asserts its census against it at `CELL_COUNT_BAND`.
  - `F_HUMAN_KEPT_ROD_CELL_GUARD = 600_000`: checked before any further reading.
- **Record comparands copied from a log, at the digits printed** (`20261010T005723Z_GEO-35-s2a-W1-wf12-buildonly.log:88, :6219`). They carry ten significant figures, so their rounding is ≤ 2.4e-10 relative, under `STRIP_BAND`:
  - `HOLE_CAVITY_SURFACES_RECORD = 240`;
  - `HOLE_CAVITY_AREA_RECORD_M2 = 2.094760146e-01`;
  - `HOLE_PHANTOM_CAD_RECORD_M3 = 2.225676428e-03`.
- One module-scoped fixture builds once. The refusal tests and the pure tests run ungated.

**Anchors (asserted; every band imported).** Write V_rod = `_rod_volume()` = 9.8437e-8 m³, A_rod = `_rod_area()` = 2.01062e-4 m², and V_cyl = πR²H from the kwargs dict.
- (F1) CAD(403) = V_rod, and `implant_rod_cad["kept_mass_m3"]` = `rod_mass_m3` = V_rod, at `STRIP_BAND`.
- (F2) CAD(3) = V_cyl − V_rod at `STRIP_BAND` — the rod / gel partition.
- (F3) CAD(3) + CAD(403) = V_cyl at `STRIP_BAND`: nothing of the phantom is lost to air or to the cavity.
- (F4) OCC wall area `implant_rod_cad["wall_area_m2"]` = A_rod at `STRIP_BAND`. Meshed dS(402) / A_rod lies inside `TERMINAL_AREA_BAND`; step 3's K4 read 0.996254 on F-small with the same rod sizing, and the prediction is the same to 4 digits.
- (F5) **The coil cavity did not swallow the rod.** In `diag["conductor_cavity_cad"]`, surfaces == `HOLE_CAVITY_SURFACES_RECORD` (240) and \|area / `HOLE_CAVITY_AREA_RECORD_M2` − 1\| ≤ `STRIP_BAND`. A swallowed wall would add A_rod, a miss of 9.6e-4 (predicted, printed). That is the failure this anchor exists to catch, not a separate build.
- (F6) Census:
  - `_census_checks(census, cells, _ring_ports())` returns no failure: conductor tag 1 absent, air, phantom and every ring port half present, census summing to `size_global`;
  - tag 403 has > 0 cells;
  - every owned 402 facet is interior and lies between one tag-3 cell and one tag-403 cell (`_interior_facet_check`, ghost-aware), and no 402 facet is exterior;
  - the 401 facets, from `_outer_box_tags`' cavity census, are printed beside 78 066 (`…buildonly.log:6220`; mesher-dependent, printed only).
- (F7) The OCC wall centroid sits at (0.015, 0, 0) to `EXACT`·ℓ (step 2's A6 pattern).
- (F8) Refusals, ungated, with `gmsh.isInitialized()` False after each:
  - `_normalise_implant_rod` with `as_hole=True` raises a `ValueError` matching `implant_rod` for `ROD` (hole) and for `dict(ROD, keep_cells=True, coat_radius=1.5e-3)`;
  - with `symmetry_quadrant=True` it raises for `ROD_FH`;
  - positive case: `_normalise_implant_rod(ROD_FH, phantom_radius=R, phantom_height=H, as_hole=True)` returns `keep_cells` True (a pure call, no gmsh).
- (F9) Cells < `F_HUMAN_KEPT_ROD_CELL_GUARD` before any facet reading.

**Printed, never gated (rule (e): no prior measurement on this fixture).**
- Cells against `HOLE_CELL_RECORD` 369 387 — the weekly's "cells printed against 369 387" — as the ratio and the difference.
  - **Predicted:** 369 387 plus ≈ 66 k, range [420 k, 455 k]. The basis is F-small's kept-rod increment at the same 0.015 m `resolution`: 182 371 − 116 085 = 66 286.
  - Tag-403 cells within ±10 % of 15 000 (F-small read 15 000 and 15 035). Tag-3 cells printed.
- Meshed 403 volume / V_rod (F-small kept: 0.986667).
- E, F, degree-1 N1curl dofs (= E) and degree-2 dofs (= 2E + 2F): I21b's price input.
- Mesh wall time and build elapsed.
- CAD(3) of the rod build against `HOLE_PHANTOM_CAD_RECORD_M3`: the difference is V_rod, to the record's rounding.

**Negative controls (rule (e)).**
- (C1h) **The miss is asserted > `STRIP_BAND`; the factor is printed.** The rod-free phantom CAD record `HOLE_PHANTOM_CAD_RECORD_M3` misses (F2)'s comparand by V_rod / (V_cyl − V_rod) = 4.42e-5, ≈ 4.4e4× the band. That ratio is also the ceiling. The control is arithmetic on a logged record of this fixture's rod-free build, the same quantity (F2) measures.
- (C2h, **predicted**, printed) The swallow miss of (F5): A_rod / 2.094760146e-01 = 9.60e-4.
- (C3h, **predicted**, printed) Swapping tags 403 and 3 misses (F1) by ≈ V_cyl / V_rod ≈ 2.3e4.

**Windows.** Bare, whole and foreground; Bash `timeout: 660000`; `-s`; no `-k`; no pipes. Each pytest window builds once. None passes ≈ 10 min, so no durable capture is needed.
- **W0 (smoke, collect).** The imports must collect in the complex build. If collection fails, stop and report; never retype a band.
  `scripts/testing/run_and_log.sh GEO-36-step6-W0-collect "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_rod_f_human.py --collect-only -q'"`
- **W1** `-n 1`, predicted 150–200 s (one F-human build ≈ 115 s, plus ≈ 6 s for the rod, plus checks).
  `scripts/testing/run_and_log.sh GEO-36-step6-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_rod_f_human.py -s -v --tb=short'"`
  Register `F_HUMAN_KEPT_ROD_CELL_RECORD` from W1's printed count **before W2**. That is the module's last edit.
- **W2** `-n 2`, the window of record: the same command with label `GEO-36-step6-W2` and `-n 2`. Every printed number must equal W1's.
- **W3a, the rule (c) re-gate of the F-small rod paths.** Real build, step 5's W3 set, predicted ≈ 255 s (`20261008T142024Z_GEO-36-step5-W3-regate.log`). It must print 167 371 / 182 371 / 116 085 unchanged, and (A7) "as_hole" must still refuse. A `test_birdcage_volumes_partition_the_box` red carrying the surface 59 / 79 string is the 🔴 OPEN known issue (line 97 above).
  `scripts/testing/run_and_log.sh GEO-36-step6-W3a-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod.py tests/mesh/test_birdcage_port_tags.py tests/mesh/test_birdcage_conductor_hole.py -v -s --tb=short'"`
- **W3b, the rule (c) re-gate of the `as_hole` default path** (the edited branch). Predicted ≈ 192 s (`…buildonly.log:6358`). (A0) 369 387 and the 401 line (240 surfaces, area 2.094760146e-01) must print unchanged.
  `scripts/testing/run_and_log.sh GEO-36-step6-W3b-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF12_STEP4=1 FEM_EM_WF12_STEP4_BUILD_ONLY=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf12_f_human_copper_degree2.py::test_wf12_step4_copper_degree2_column -v -s --tb=short'"`
- Order: W0 → W1 → register the record → W2 → W3a → W3b. Two heavy windows never run at once.

**Traps already paid for.**
- Interior dim-2 gmsh groups hang `model_to_mesh` at `-n 2` (known-issues 9). 402 stays the dolfinx-side (3, 403) rebuild; never emit it as a gmsh group in the kept branch (step 3's disclosed deviation, line 103 above).
- Enter the rod into the fragment **once** (`GEO-9`).
- A failed gmsh build contaminates its process (`GEO-23`): one build per window, never a retry. A FAIL is recorded verbatim with its surface pairs.
- Read walls from the meshed volumes, never from the tools' faces (`TH-15` step 0, 716 orphan nodes). The cavity read runs before the conductor drop and must see the implant pieces; that is (F5)'s point.
- `cell_tags.values` is rank-local: use `_owned_census` and the ghost-aware `_cell_values_with_ghosts`, reduce every count, and keep collectives out of `rank == 0` branches.
- Rule (i): no edit to the module after W2 starts. The record is registered between W1 and W2, then frozen. W3a and W3b run the committed `src/`.
- Status 124: run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first, and the rubric's cache sweep only if it reads 0.
- Never import `tests/validation/test_imp3_*`, `test_wf9_*` or `test_wf10_*` here. They are heavy module imports, and nothing is needed from them.

**Scope.** Not in this item:
- an EM solve, a material map, SAR or a tip ratio (all of that is I21b);
- a hole rod or a coat under `as_hole`;
- `symmetry_quadrant`;
- the solid F-human route;
- a `_build_f_human` / `_f_human_kwargs` edit (I21b's);
- a placement other than `ROD`'s;
- any refinement knob.

`GEO-36` stays ✅.

**Negative result.** Report and stop. Park on `attempt/GEO-36-step6-<ts>`, open a known-issues entry and mark item 372 BLOCKED in the same commit (rule (d)); that makes item 381 skip. The reds:
- a gmsh FAIL (record the surface pairs verbatim);
- (F5) red: the rod wall reached 401, so the composition is wrong, not the band;
- (F1)–(F4), (F6) or (F7) red;
- the guard tripping.

If W3a or W3b moves a record (167 371 / 182 371 / 116 085 / 369 387 / the 401 line), revert `src/` and stop. A printed reading outside its prediction is a reading, never a red. There is no `PREDICTED RED` line.

**Status it can move.** The §7 `GEO-36` row stays ✅. Its Done-when is "the identities asserted in one logged window at `-n 2`, cells and time recorded ⇒ ✅" (`PROJECT_PLAN.md:1873`).
- The row gains a **step-6 record** in the chunk file, plus one clause in the row if it stays ≤ 2 KB. The record covers: `implant_rod` with `keep_cells` composes with `as_hole`; the F1–F5 values; cells against 369 387; tag-403 cells; build time.
- **No Done-when amendment.** No glyph moves, so rule (j) does not bind.
- I21b (item 381) becomes runnable. The slot makes no §10 edit.

**Execution order.** grep `src/` → code → W0 → W1 → register the record → W2 → W3a → W3b → commit (`src/`, module, harness logs, the row clause) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 372 — `GEO-36` step 6 (§10 I21a) — per `docs/planning/chunks/GEO-36.md` §9 item 372. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W0, W1, W2, W3a, W3b in order, foreground. Register `F_HUMAN_KEPT_ROD_CELL_RECORD` from W1 before W2 and make no module edit after W2 starts. The `src/` change is limited to lifting the `as_hole` refusal for a kept, uncoated rod and entering the rod once into the `as_hole` fragment; disclose it. Never retry a gmsh failure. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-step6-*`, mark item 372 BLOCKED in the same commit, and report."

**Review ruling 2026-10-11 (03:00 daily) on item 372.** (1) The narrow `src/` change (lift the `implant_rod` × `as_hole` refusal for a kept, uncoated rod only; add the rod to the hole fragment once) is accepted under rule (c), re-gated by W3a (the F-small rod tests) and W3b (the F-human hole build-only window reprinting 369 387 and the 401 cavity line) in the same slot. (2) Rod placement: the F-small rod's centre (0.015, 0, 0) is reused, so the tip coordinates and the closed-form helper carry to I21b; the hotspot-c\* placement is not chosen (it is not the weekly's and would change I21b's comparand). (3) Windows run at `timeout -k 30 420` (heavy ceiling) because one hole build is ≈ 115 s.

### §9 item 427 — 2026-10-13 (§10 chain I24; chain step, no solve)

`GEO-36` step 7: `implant_spheroid`, an oblate PEC body (a, a, c) in the phantom, in `implant_rod`'s pattern — hole and `keep_cells` on F-small, kept on the F-human `as_hole` phantom; CAD identities asserted at `-n 1` / `-n 2` (implementer; standard by the weekly, heavy ceiling per window as step 6's ruling (3); real build for the F-small module, complex for the F-human module; ≈ 7 measured / 70 predicted min: ≈ 22 min of windows + 15 fixed + ≈ 35 code; depends on nothing). Licence: §10 I24 (the 10-12 weekly, `PROJECT_PLAN.md:4758–4762`), Phase 7 (7b), §9 rule (4). **Target: CAD identities, the weekly's text. No solve; no band, tolerance or Done-when moves; every band imported.** A record on a ✅ row. Item 428 (`IMP-3` step 9, I25) consumes it.

**Geometry (one body for every build).** a = 6 mm (equatorial), c = 1.2 mm (along the symmetry axis), aspect 5 — `IMP-6`'s FEM-gated aspect-5 body (step 1, item 404; step 3, item 409), so I25's closed form is the one `IMP-6` validated. Centre (0.015, 0, 0), `GEO-36`'s rod centre (step 6 ruling (2): tips, masks and helpers carry over). **Symmetry axis x̂**: the disc lies in the yz-plane. Reason: the coil-driven field there is ≈ 99 % E_z under P1 (`IMP-3` step 1 measured η_z = 0.990426 on the rod axis at x = 0.015, `20261002T170843Z_IMP-3-step1-W2.log.gz:3864`), so with the normal along x̂ the dominant field lies in the disc's plane and its rim points (0.015, 0, ±a) take the in-plane enhancement 1/L_x = 8.0155 — `IMP-6`'s in-plane case, the disc-edge analogue of the rod tips (0.015, 0, ±0.016). Closed forms (computed here): e = √(1 − c²/a²) = 0.979796, V = 4πa²c/3 = 1.809557e-7 m³, A = 2πa²[1 + ((1 − e²)/e) artanh e] = 2.473638e-4 m², ρ_r = c²/a = 0.24 mm. Strictly inside both phantoms: radial extent ≤ √(16.2² + 6²) = 17.3 mm < 30 mm (F-small R), \|z\| ≤ 6 mm ≪ 40 mm.

**Why not `addSphere` + `dilate` (measured).** `IMP-6`'s builder makes the oblate body that way, and OCC's dilated sphere is a B-spline approximation: its volume misses 4πa²c/3 by **4.37e-4** at exactly this a, c (air OCC mass 5.111123404e-06 vs 5.111044263e-06 m³, `20261012T184107Z_IMP-6-step1-W0.log:503`; 7.914e-11 m³ over V), 4.4e3× `CAD_ANALYTIC_BAND`. The prolate path shows the same size of error (`:259`, `:324`). So a volume identity against the closed form is unmeetable on that construction. **And** item 408 measured that the dilated body's rim can be truncated by the mesher (31 / 93 air midpoints inside the true spheroid at finer rim sizing, `IMP-6.md` § item 408 outcome). Both are answered by one construction: revolve an exact half-ellipse profile made of two quarter-arcs, so the body is the exact ellipsoid of revolution and the rim circle is a model edge.

**`src/` change (rule (c), additive, disclosed).** In `src/fem_em_solver/io/mesh.py`:
- `birdcage_port_domain` (`:4158`) and `_build_birdcage_port_model` (`:4907`) gain keyword-only `implant_spheroid: Optional[Dict[str, object]] = None`. Keys: `semi_major` a [m], `semi_minor` c [m] (along `axis`), `centre`, `axis` (normalised), `keep_cells` (bool, default False). Docstrings gain one paragraph each.
- `_normalise_implant_spheroid` (pure, beside `_normalise_implant_rod` at `:4773`), called before any gmsh state exists. Refuse with `ValueError` (message containing `implant_spheroid`): `implant_rod` also set (one implant per build); `symmetry_quadrant`; `as_hole` unless `keep_cells is True` (step 6's rule for the rod, `:4796–4797`); missing keys, a ≤ c or c ≤ 0, non-finite values; non-3-vector or zero `axis` / `centre`; non-bool `keep_cells`; a body not strictly inside the phantom (conservative bounding-sphere test: \|centre_xy\| + a < R and \|centre_z\| + a < H/2).
- **Construction** (in the implant block at `:5340`, a sibling branch of the rod's): in the local frame, two quarter-ellipse arcs (pole (0, 0, c) → rim (a, 0, 0) → pole (0, 0, −c), centre the origin, major axis x) closed by the axis segment; `addPlaneSurface`; `occ.revolve(…, 0, 0, 0, 0, 0, 1, 2π)`; rotate π/4 about local z (the revolve seam half-way between rim points, `IMP-6`'s 45° rule); rotate local ẑ onto `axis` (about ẑ × axis; identity for +ẑ, π about x̂ for −ẑ); translate to `centre`. For axis x̂ the local rim direction x maps to world ∓z, so the seam lies at 45° from the rim points (0.015, 0, ±a) — print the seam edge's bounding box. Assert one solid (the rod's fuse check pattern, `:5359–5362`); record its OCC mass. Verify the arc API in the image first (`gmsh.model.occ.addEllipseArc(startTag, centerTag, majorTag, endTag)`, or `addEllipse` with explicit `zAxis` / `xAxis`).
- **Shared downstream, keyed exactly as the rod's** (the cut-first sequence `:5373–5384` with `removeTool=not keep`; the `as_hole` fragment entry of the kept body once — step 6's change, `GEO-9`; the piece policy → `BIRDCAGE_IMPLANT_CELL_TAG` 403 (`:5557`); the wall read `:5704–5760` — hole: single-use faces inside the padded bbox, kept: the implant pieces' boundary with the two-use check; 402 emitted as an exterior gmsh group on the hole build and rebuilt dolfinx-side on the kept build, `:4731–4739`). Recommended: generalise the rod's local variables into one "implant body" record (`kind`, `dimtags`, `bbox`, `keep`) so each shared line keys on it; every rod hunk's behaviour must be unchanged (the re-gates W3a–W3c below are the proof). Reuse tags 402 / 403 (one implant per build), so `IMP-3`'s wall-dof and census helpers work unchanged.
- **Sizing** (beside the rod's at `:6241–6281`; no Balls — item 408): a `Distance` from the body's wall surfaces → `Threshold` SizeMin c/3 = 0.4 mm, SizeMax `resolution`, DistMin c, DistMax 10c (the rod's r/3 over r → 10r with r ↔ c, the thin dimension); plus `IMP-6` step 1's rim tube, a `MathEval` distance to the rim circle written in the world frame from `centre` and `axis` (d = x − centre, z_l = d·axis, ρ_l = √(\|d\|² − z_l²), dist = √((ρ_l − a)² + z_l²); numbers baked with `repr`, string printed) → `Threshold` SizeMin ρ_r, SizeMax `resolution`, DistMin ρ_r, DistMax 4ρ_r. Print one sizing line as the rod does (`:6275–6280`). No refinement knob.
- **Diagnostics:** `implant_spheroid_cad` (`None` when off): surfaces, `wall_area_m2`, `body_mass_m3`, `wall_centroid_m`, `facet_tag`, the normalised request, `keep_cells`, `wall_size_m` (c/3), `rim_size_m` (ρ_r); with `keep_cells` also `cell_tag` and `kept_mass_m3`. `implant_rod_cad` stays `None` on a spheroid build.

**New module `tests/mesh/test_birdcage_implant_spheroid.py`** (F-small; real build; imports by name, never a `test_*` symbol):
- From `tests/mesh/test_birdcage_conductor_hole.py`: `CAD_ANALYTIC_BAND` (`:123`, 1e-7), `EXACT` (`:117`, 1e-9), `SOLID_CELL_RECORD` (`:95`, 116 085), `_owned_census` (`:126`). From `tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`: `CELL_COUNT_BAND`. From `tests/mesh/test_birdcage_leg_gaps.py:73`: `TERMINAL_AREA_BAND` (0.95, 1.0). From `tests/validation/test_ops61_exact_strip_ports.py:58`: `STRIP_BAND` (1e-9). From `tests/mesh/test_birdcage_implant_rod.py`: `_build_rod` (`:125`, called with `implant_rod=None, implant_spheroid=SPHEROID` — its kwargs are the four-port F-small call), `_interior_facet_check` (`:184`), `_cell_values_with_ghosts` (`:171`), `KEPT_ROD_CELL_GUARD` (`:110`, 320 000), `PHANTOM_CELL_TAG`. Phantom R, H from `tests/mesh/test_birdcage_port_tags.py` (`PHANTOM_RADIUS`, `PHANTOM_HEIGHT`).
- Constants: `SPHEROID_SEMI_MAJOR = 0.006`, `SPHEROID_SEMI_MINOR = 0.0012`, `SPHEROID_CENTRE = (0.015, 0.0, 0.0)`, `SPHEROID_AXIS = (1.0, 0.0, 0.0)`, `SPHEROID` (the dict), `RIM_POINTS = {"top": (0.015, 0.0, +0.006), "bottom": (0.015, 0.0, -0.006)}`; pure `_spheroid_volume()`, `_spheroid_area()`, `_inside_true_spheroid(points)` (s = ρ_l²/a² + z_l²/c² in the body frame); `_build_spheroid(**overrides)` and `_build_kept_spheroid(**overrides)` (step 3's `_build_kept_rod` pattern, guard checked before return). These are I25's imports.
- **Records registered here** (first readings, version-tagged, printed and not asserted in this item; registered from W1a, printed again by W2a): `SPHEROID_HOLE_CELL_RECORD`, `KEPT_SPHEROID_CELL_RECORD`. Build order in the module: flag-off, then hole, then kept (a failed gmsh build contaminates its process, `GEO-23`).

**New module `tests/mesh/test_birdcage_implant_spheroid_f_human.py`** (complex build; step 6's module `tests/mesh/test_birdcage_implant_rod_f_human.py` as the template): `_f_human_kwargs(as_hole=True)` (`tests/validation/test_wf7_f_human_b1_quadrature.py:67`) plus `implant_spheroid = dict(SPHEROID, keep_cells=True)`; from step 6's module the constants `HOLE_CAVITY_SURFACES_RECORD` (240), `HOLE_CAVITY_AREA_RECORD_M2` (2.094760146e-01), `HOLE_PHANTOM_CAD_RECORD_M3`, `F_HUMAN_KEPT_ROD_CELL_GUARD` (600 000) — **never its `kept` fixture or a `test_*` name**; `_census_checks` (`tests/validation/test_wf12_f_human_hole.py:218`), `_ring_ports` (`tests/validation/test_birdcage_f_human_rung.py:134`), `HOLE_CELL_RECORD` (`tests/validation/test_wf12_f_human_copper_matrix.py:55`), `_outer_box_tags` (`tests/validation/test_th14_birdcage_copper.py:133`). One module-scoped build. Record: `F_HUMAN_KEPT_SPHEROID_CELL_RECORD`, registered from W1b, printed again by W2b.

**Anchors (asserted; every band imported).** F-small module:
- (S1) OCC body mass = V at `CAD_ANALYTIC_BAND`. (S2) OCC wall area `wall_area_m2` = A at `CAD_ANALYTIC_BAND`. Predicted ≤ 1e-9 each (an exact surface of revolution; no prior on this construction).
- (S3) Hole build: CAD(3) = πR²H − V at `CAD_ANALYTIC_BAND`; ΣCAD(hole) + `body_mass_m3` = ΣCAD(flag off) at `EXACT` (step 2's A1 / A2).
- (S4) Kept build — **the body / gel partition**: CAD(403) = `kept_mass_m3` = `body_mass_m3` at `STRIP_BAND`; CAD(3) + CAD(403) = πR²H at `STRIP_BAND`; ΣCAD(kept) = ΣCAD(flag off) at `EXACT` (step 3's K2 / K3, step 6's F1–F3).
- (S5) Meshed tag-402 area / A in `TERMINAL_AREA_BAND`: `ds` on the hole build, `dS("+")` on the kept build (step 3's trap: `ds` reads 0 on an interior group).
- (S6) Census: hole — cell-tag set = flag-off set, owned census sums to `size_global`, 402 > 0 facets each exterior with its one cell in tag 3, 403 absent; kept — set = flag-off ∪ {403}, 402 all interior, each between one tag-3 and one tag-403 cell (`_interior_facet_check`, ghost-aware); flag-off cells = `SOLID_CELL_RECORD` within `CELL_COUNT_BAND`.
- (S7) The OCC wall centroid at `SPHEROID_CENTRE` to `EXACT`·a (all three components printed; the meshed centroid printed only).
- (S8) **Rim fidelity (item 408's census, made a gate).** On both builds, owned tag-3 cells whose midpoint lies within 2ρ_r of either `RIM_POINTS` entry: count > 0, and the number of those midpoints inside the true spheroid **= 0** (reduced). On the kept build, every tag-403 midpoint is inside the true spheroid (s < 1 + 1e-9; a convex body meshed with surface nodes on the CAD surface cannot do otherwise).
- (S9) Refusals, each before any gmsh state (`gmsh.isInitialized()` False after), through `_normalise_implant_spheroid`: rod + spheroid; `symmetry_quadrant`; a hole spheroid with `as_hole`; a = c; c ≤ 0; zero axis; a body poking out of the phantom; `keep_cells="yes"`. Positive pure case: the kept request with `as_hole=True` returns `keep_cells` True.
- (S10) Cells < `KEPT_ROD_CELL_GUARD` (320 000) on both builds before any facet reading.

F-human module (step 6's F1–F9 with the body for the rod): (H1) CAD(403) = `kept_mass_m3` = `body_mass_m3` at `STRIP_BAND`, and `body_mass_m3` = V at `CAD_ANALYTIC_BAND`; (H2) CAD(3) = V_cyl − `body_mass_m3` at `STRIP_BAND`; (H3) CAD(3) + CAD(403) = V_cyl at `STRIP_BAND` (V_cyl = πR²H from the kwargs dict, never retyped); (H4) OCC wall area = A at `CAD_ANALYTIC_BAND`, meshed dS(402) / A in `TERMINAL_AREA_BAND`; (H5) the coil cavity did not swallow the body: `conductor_cavity_cad` surfaces = 240 and \|area / `HOLE_CAVITY_AREA_RECORD_M2` − 1\| ≤ `STRIP_BAND`; (H6) `_census_checks(census, cells, _ring_ports())` empty, tag 403 > 0 cells, every owned 402 facet interior between one tag-3 and one tag-403 cell; (H7) the centroid as (S7); (H8) rim fidelity as (S8); (H9) cells < 600 000 before any facet reading.

**Printed, never gated (rule (e): no prior measurement on this fixture).** Cells per build against 116 085 / 369 387 — **predicted**: F-small hole 116 085 + [35 k, 85 k] (the rod added 51 286 at r/3 = 0.33 mm over 201 mm²; the disc is c/3 = 0.4 mm over 247 mm² plus the rim tube), kept + [10 k, 30 k] tag-403 cells (the rod read 15 000), F-human kept 369 387 + [45 k, 115 k]; meshed 403 volume / V predicted [0.97, 0.995] (rod: 0.986667); meshed 402 / A predicted [0.985, 0.999]; rim-mask cell counts; E, F, degree-1 and degree-2 N1curl dofs (I25's price input); mesh wall and build time; the seam bounding box.

**Negative controls (rule (e)).**
- (C1s, **asserted**) The dilated-sphere construction misses (S1)'s comparand. A pure-OCC test (scratch model, no mesh, finalised after): `addSphere(0, 0, 0, 1)`, rotate π/4 about z, `dilate(a, a, c)`, `getMass` against V. Backed by the same comparison on the same solid: 4.374e-4 (`20261012T184107Z_IMP-6-step1-W0.log:503`, the air-mass difference over V; OCC is deterministic). Assert the miss > `CAD_ANALYTIC_BAND`; print it beside 4.374e-4 and the separation (≈ 4.4e3×). It is a measured error, not a saturating quantity, so no ceiling bounds it.
- (C2s, **asserted**) The flag-off phantom CAD misses (S3)'s comparand by V/(πR²H − V) = 8.006e-4 (≈ 8.0e3× the band; that factor is also the ceiling). Backed by the flag-off phantom CAD = πR²H to 1.1e-16 and step 2's (C1), the same comparison with V_rod, reading its prediction exactly (`20261002T110559Z_GEO-36-step2-W2.log.gz:3616–3617`). Assert the miss > `CAD_ANALYTIC_BAND`; the factor printed.
- (C3s, **predicted**, printed) The hole build misses (S4)'s `EXACT` sum by `body_mass_m3`/ΣCAD (≈ 1.6e-5, step 3's C1k read 8.54e-6 for the rod).
- (C4s, **predicted**, printed) (H5)'s swallow miss: A / 2.094760146e-01 = 1.181e-3.

**Windows.** Bare, whole, foreground; Bash `timeout: 660000`; `-s`; no pytest `-k`; no pipes. Cost basis (measured): three F-small builds in one module, 188 s at `-n 1` / 181 s at `-n 2` (`20261004T200412Z_GEO-36-step3-W1.log`, `…200730Z_…-W2.log`); one F-human kept build, 184 s / 180 s (`20261011T185757Z_GEO-36-step6-W1.log`, `…190108Z_…-W2.log`); the re-gates 257 s / 191 s (`…190414Z_…-W3a-regate.log`, `…190838Z_…-W3b-regate.log`).
- **W0 (smoke, collect, complex):** `scripts/testing/run_and_log.sh GEO-36-step7-W0-collect "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_spheroid.py tests/mesh/test_birdcage_implant_spheroid_f_human.py --collect-only -q'"`
- **W1a** `-n 1`, F-small, predicted 170–260 s: `scripts/testing/run_and_log.sh GEO-36-step7-W1a "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 420 mpiexec -n 1 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_spheroid.py -s -v --tb=short'"`
- **W1b** `-n 1`, F-human, predicted 170–240 s: `scripts/testing/run_and_log.sh GEO-36-step7-W1b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_spheroid_f_human.py -s -v --tb=short'"`
- Register the three records from W1a / W1b — the modules' last edit.
- **W2a / W2b** `-n 2`, the windows of record: the W1a / W1b commands with labels `GEO-36-step7-W2a` / `GEO-36-step7-W2b` and `-n 2`. Every printed number must equal the `-n 1` window's.
- **W3a** (rule (c) re-gate, F-small rod paths; step 6's W3a verbatim, relabelled): `scripts/testing/run_and_log.sh GEO-36-step7-W3a-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod.py tests/mesh/test_birdcage_port_tags.py tests/mesh/test_birdcage_conductor_hole.py -v -s --tb=short'"` — 167 371 / 182 371 / 116 085 unchanged, (A7) still refuses. A `test_birdcage_volumes_partition_the_box` red carrying the surface 59 / 79 string is the 🔴 OPEN known issue (step 2's outcome above), not this change.
- **W3b** (the `as_hole` default path; step 6's W3b verbatim, relabelled): `scripts/testing/run_and_log.sh GEO-36-step7-W3b-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF12_STEP4=1 FEM_EM_WF12_STEP4_BUILD_ONLY=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf12_f_human_copper_degree2.py::test_wf12_step4_copper_degree2_column -v -s --tb=short'"` — 369 387 and the 401 line (240, 2.094760146e-01) unchanged.
- **W3c** (the kept rod under `as_hole`, the path this change shares): `scripts/testing/run_and_log.sh GEO-36-step7-W3c-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_rod_f_human.py -s -v --tb=short'"` — `F_HUMAN_KEPT_ROD_CELL_RECORD` 436 541 and F1–F9 unchanged.
- Order: grep `src/` → code → W0 → W1a → W1b → register → W2a → W2b → W3a → W3b → W3c → commit (`src/`, both modules, logs, the row clause) → journal. Two heavy windows never at once.

**Traps already paid for.**
- Interior dim-2 gmsh groups hang `model_to_mesh` at `-n 2` (known-issues 9): on the kept build 402 is the dolfinx-side (3, 403) rebuild, never a gmsh group (step 3's disclosed deviation).
- Enter the kept body into the fragment **once** (`GEO-9`); cut it from the phantom and the air box before the fragment, or the air box fills the hole (step 2's build note).
- A failed gmsh build contaminates its process (`GEO-23`): one build order, never a retry; a FAIL is recorded verbatim with its surface pairs.
- Walls from the meshed volumes, never from the tool's faces (`TH-15` step 0, 716 orphan nodes).
- `cell_tags.values`, `assemble_scalar` and local counts are rank-local: `_owned_census`, ghost-aware `_cell_values_with_ghosts`, every count reduced, collectives out of `rank == 0` branches. (S8)'s midpoints come from owned cells only.
- `MathEval` takes a literal expression (bake with `repr`, print it); a `Threshold` returns SizeMax beyond DistMax — SizeMax is `resolution`, never fine.
- The generator's `optimize("Netgen")` stays (step 2's (A4) margin covered it); do not add or drop a pass.
- Rule (i): no module edit after W2a starts; the records are registered between the `-n 1` and `-n 2` windows and then frozen. Status 124: `pgrep -c python3` first, the rubric's sweep only if it reads 0. Never import `tests/validation/test_imp3_*`, `test_wf9_*` or `test_wf10_*` here.

**Scope.** Not in this item: a solve, a PEC condition, SAR or a rim ratio (all I25); a hole body under `as_hole`; a coat; `symmetry_quadrant`; the solid F-human route; any refinement knob; a placement or axis other than the one above; a square-edged plate (a real plate edge carries a wedge singularity the smooth spheroid does not). `GEO-36` stays ✅.

**Status it can move.** The §7 `GEO-36` row stays ✅; its Done-when is "the identities asserted in one logged window at `-n 2`, cells and time recorded ⇒ ✅. No solve." (`PROJECT_PLAN.md:1908`), which this item meets as a record. The row gains one clause — "Step 7: `implant_spheroid` (aspect-5 oblate disc, exact revolved CAD, rim as an edge), hole and kept on F-small and kept under `as_hole` on F-human, identities green at `-n 2`, cells" — if it stays ≤ 2 KB; the record goes in this file. No Done-when amendment, so rule (j) does not bind. I25 (item 428) becomes runnable.

**Negative result.** Report and stop: park on `attempt/GEO-36-step7-<ts>`, open a known-issues entry and mark item 427 BLOCKED in the same commit (rule (d)) — that makes item 428 skip. The reds:
- a gmsh FAIL or a revolve that does not give one solid (record verbatim; no fallback to `dilate`, whose CAD cannot meet (S1) by (C1s));
- (S1) / (S2) / (H1) / (H4) red against the closed form with (S4) / (H3) green: OCC's integration on the revolved surface, not the generator — print the miss; it is a rule (f) question for the review, never a widened band;
- (S8) / (H8) red: the rim is truncated even as a model edge — the item-408 defect on this generator; print the inside count and the deepest midpoint;
- (H5) red: the body wall reached 401 (the composition, not the band);
- any other identity red, or a guard trip.
If W3a, W3b or W3c moves a record (116 085 / 167 371 / 182 371 / 369 387 / 436 541 / the 401 line), revert `src/` and stop. A printed reading outside its prediction is a reading, never a red. No `PREDICTED RED` line.

**Every-item rules** (PROJECT_PLAN §9 "Every item"): bare whole harness calls; a lost docker socket stops the slot; no edit after the last window; `timeout: 660000` and `-s`; never delete a log; BLOCKED in the parking commit; a killed window is not a result; two heavy windows never at once; no slot opens `docs/private/` or `aed_results/`.

**Executor spawn prompt (verbatim).** "Execute §9 item 427 — `GEO-36` step 7 (§10 I24) — per `docs/planning/chunks/GEO-36.md` §9 item 427. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W0, W1a, W1b, W2a, W2b, W3a, W3b, W3c in order, foreground. Register the three cell records from W1a / W1b before W2a and make no module edit after W2a starts. The `src/` change is the additive `implant_spheroid=` keyword (revolved exact body, rim as an edge, tags 402 / 403 reused); every rod path must re-gate unchanged; disclose it. Never retry a gmsh failure and never fall back to `dilate`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-step7-*`, mark item 427 BLOCKED in the same commit, and report."

**Outcome (2026-10-13, 15:00 slot): DONE.** All identities green at `-n 1` and `-n 2` with every printed number equal across widths. (S1) body mass rel 2.2e-16, (S2) wall area rel 2.45e-11 against the closed forms; (S3) / (S4) / (H1)–(H3) at ≤ 1.1e-16; meshed 402 / A 0.997686 (hole, kept, F-human); (H5) cavity 240 surfaces, rel 1.5e-10; (S8) / (H8) 120 / 128 rim-mask midpoints, 0 inside the true spheroid, every 403 midpoint inside (max s 0.98308 / 0.98206); (C1s) dilated sphere misses by 4.3735e-4 (4.37e3× the band), (C2s) 8.006405e-4 as predicted, (C3s) 1.570796e-5. Records: F-small hole 150 140 (116 085 + 34 055 — just under the predicted + [35 k, 85 k], a reading), kept 164 266 (403: 14 126), F-human kept 417 653 (369 387 + 48 266; 403: 14 047); meshed 403 / V 0.995984 (just over the predicted [0.97, 0.995]). Degree-2 N1curl dofs: F-small hole 972 888, kept 1 051 444, F-human 2 875 948. Re-gates unchanged: 116 085 / 167 371 / 182 371 (W3a), 369 387 and the 401 line (W3b), 436 541 (W3c). Disclosed: the revolve leaves its profile face as a free surface at the origin; left in, the first W1a build failed in gmsh boundary recovery ("PLC Error: A segment and a facet intersect at point", `20261013T201959Z_GEO-36-step7-W1a.log:2171`). That face is now removed, as `IMP-6`'s `oblate_rim="edge"` does. The build before that one aborted on a numpy-scalar `repr` in the `MathEval` string (`20261013T201740Z_…-W1a.log:1821`). Both were defects in the new code, fixed with their cause named, and not retries of a CAD construction. The seam-bbox diagnostic reads `None` (each seam spans c, not 2c). The seam curves are printed per curve: bbox y, z ±4.2427e-3, i.e. 45° from the rim points. Logs: `20261013T201731Z_GEO-36-step7-W0-collect.log`, `…202258Z_…-W1a.log` (192 s), `…202620Z_…-W1b.log` (177 s), `…202931Z_…-W2a.log` (193 s), `…203250Z_…-W2b.log` (172 s), `…203550Z_…-W3a-regate.log` (271 s), `…204028Z_…-W3b-regate.log` (188 s), `…204342Z_…-W3c-regate.log` (185 s).

### Review ruling 2026-10-13 03:00 (daily) on the item 427 draft

Queued as drafted: `implant_spheroid` builds the body by revolving an exact half-ellipse (two quarter-arcs), not by stretching a sphere — the stretched construction misses 4πa²c/3 by 4.37e-4 at a = 6 mm, c = 1.2 mm (`20261012T184107Z_IMP-6-step1-W0.log:503`), ≈ 4.4e3× `CAD_ANALYTIC_BAND`, so it becomes the asserted negative control. The exact revolve also makes the rim circle a model edge; whether `IMP-6`'s private builder later adopts it is item 426's licensed variant (1), not this item's.

### §9 item 530 — 2026-10-18 (§10 I33, `GEO-36` step 7 — implant_rod in C2's tagged box; chain step, build-only)

*Drafted read-only for the 2026-10-18 03:00 daily review; the review rules on it.*

**Numbering.** The chunk already has a step 7: item 427's spheroid (above). The weekly's "step 7" is the eighth step of this chunk. The header keeps the weekly's words; window labels say `rodbox` so that no log collides with step 7's. The review should renumber this as step 8 (ruling 1).

**Licence.** §10 I33 (`PROJECT_PLAN.md:4835–4840`, the 10-17 weekly); Phase 7 (7b); §9 rule (4). It is the kept rod of step 6 (item 372) in C2's ASTM box (`GEO-38` step 2, item 470), composed through the path item 522 opens. **Target: CAD identities only, in the weekly's words. No solve.** No band, tolerance, record or Done-when moves; every band is imported. It is a record on a ✅ row.

**Executor and cost.**
- **Executor:** implementer.
- **Tier:** heavy, build-only. W0 / W1 / W2 run real; R1 runs complex. Container `timeout -k 30` is the backstop, and item 522's gmsh stop rule (`FEM_EM_GMSH_STOP_RULE=1`, 450 s) is armed in W1 and W2.
- **Ranks:** `-n 4`, the window of record and the width at which C2 and I31 were measured, then `-n 2`. That is the weekly's "-n 2 / 4". gmsh runs on rank 0 either way (`io/mesh.py:4856–4890`).
- **Measured basis:**
  - C2's build: gmsh 310.01 s, window 339 s at `-n 4` (`20261015T185056Z_GEO-38-s2-W3-n4.log:11894, :11925`); 483.82 s at `-n 2` in the same slot (`…184154Z_GEO-38-s2-W2-n2.log:11878`).
  - I31 attempt 1, killed at 591 s (`20261016T020315Z_IMP-6-step4-W1-n4.log:5691`).
  - **I31 attempt 2's W1 price T**, from item 522. That price is this item's gate.
  - On F-human the kept rod cost +7.5 s of gmsh over the kept disc at `-n 2`: 166.17 s (`20261011T190108Z_GEO-36-step6-W2.log:6261`) against 158.68 s (`20261013T203250Z_GEO-36-step7-W2b.log:6264`).
- **Predicted cost per build window:** T_rod = T + [−40, +40] s, printed. The window is ≈ T_rod + 35 s, and a stopped window is ≈ 465 s.
- **Class:** chain no-solve. **Cost:** ≈ 18.6 measured (class median) / ≈ 55 predicted slot-minutes. That is 15 fixed + ≈ 15 code (≈ 25 lines of `src/`, a ≈ 300-line module lifted from item 522's disc module) + ≈ 20 of windows (W0 ≈ 15 s, R1 180–185 s, W1 and W2 ≤ ≈ 495 s each) + 5 to commit.
- **No take-next**, in either direction.
- **Slot-time gate:** item 522's, verbatim. Only in a 12:00 / 13:30 / 15:00 slot (start line ≤ 5 min late), and only as the slot's first item; any other slot skips it silently.
- **Family:** `GEO-36` step 8 (the weekly's step 7), attempt **1 of 4**.

**Depends on item 522, through its price.** "The price" means one line of item 522's W1 log (`…_IMP-6-step4-a2-W1-n4.log`) at `-n 4`. The gate is read off `main`, never off a branch:
1. **GO.** Item 522 is DONE on `main`; its W1 log prints `[gmsh-stop-rule] under the rule` and then `[birdcage-mesh] mesh wall time T s` with **T ≤ 450.00**; and (D0)–(D10) are green. Item 530 may then run in any later 12:00 / 13:30 / 15:00 slot, with its own builds under the same rule.
   - A T in (400, 450] predicts T_rod near the rule. That is printed, never a gate.
2. **BLOCKED, to `xl`.** Item 522's W1 log carries `[gmsh-stop-rule] STOPPED`, or W1 was killed at the container ceiling. Item 522's own commit marks item 530 BLOCKED: "to `xl` behind I31's 10-25 window; the weekly commissions". No daily re-poses it at heavy.
3. **WAIT.** Item 522 has not run, was skipped by the slot gate, or is OPEN for "W1 only". Item 530 is not runnable, and a slot skips it to the next item.
4. **BLOCKED with item 522.** Item 522 is BLOCKED for a defect: an R-gate not bitwise, a (D) red, or a gmsh FAIL. The composition is shared, so item 530 waits for the review's ruling on both.

**Code state** (`main`, 2026-10-18, before item 522's lift; the lift shifts these lines by ≤ +30).
- **The generator.** `MeshGenerator.birdcage_port_domain`, `src/fem_em_solver/io/mesh.py:4245`.
  - The `phantom_box` refusal is at `:4552–4555`. Item 522 narrows it to `"phantom_box cannot be combined with implant_rod"`.
  - The box is normalised at `:4563–4574`, where `phantom_radius` becomes the box's xy circumradius (`:4573`).
  - The rod's public normalise call is at `:4578–4584` and the builder's at `:5261–5267`. The builder already takes `phantom_box` (`:5257`).
- **`_normalise_implant_rod`** (`:5054–5222`) has three cylinder containment loops: rod `:5132–5141`, coat `:5166–5176`, tip ball `:5201–5212`. Its `as_hole` rule (kept and uncoated only) is at `:5074–5086`.
- **One implant per build:** `_normalise_implant_spheroid` refuses `implant_rod` (`:4995–4999`).
- **Shared rod and spheroid lines in the builder:** the kept flag `keep_rod_cells` at `:5688`; the rod mass `rod_occ_mass` at `:5728`. Item 522's change (c), parked `if implant_rod is not None or implant_spheroid is not None` under `:5910–5916`, already subtracts the rod's mass. The diagnostics `implant_rod_cad` (`wall_area_m2`, `rod_mass_m3`, `kept_mass_m3`) are at `:6283–6306`.
- **C2's tagged-box module:** `_box_kwargs()` at `tests/validation/test_geo38_f_clinical_box.py:54–66`. Its frame comment is at `:60`: Lz = 2A = 0.65 m on the coil axis, Lx = 2B = 0.42 m horizontal, Ly = 2C = 0.09 m vertical gel depth. `TABLE_OFFSET_M = 0.0` is at `:51`.
- **`implant_rod`:** `ROD` at `tests/mesh/test_birdcage_implant_rod.py:93–97` (r = 1 mm, ℓ = 30 mm between cap centres, axis ẑ); `_rod_volume` / `_rod_area` at `:113` / `:117`; `_cell_values_with_ghosts` at `:171`; `_interior_facet_check` at `:184`; pure `test_a7_refusals_raise_without_meshing` at `:593`.
- **Step 6's F-human kept-rod module:** `tests/mesh/test_birdcage_implant_rod_f_human.py`, with `F_HUMAN_KEPT_ROD_CELL_RECORD = 436541` at `:59`.

**Fixture.**
- `_box_kwargs()` plus `implant_rod = dict(ROD, keep_cells=True, centre=(rod_offset_m, off, 0.0))`.
  - The axis is ẑ, which is **the 650 mm axis**: `size_m[2]` = 0.65.
  - `off` is the box's `table_offset_m`, at mid-depth.
- **`rod_offset_m`** is a module parameter, `ROD_OFFSET_M: Optional[float] = None`, set to a value only with a cited source (the `TABLE_OFFSET_M` precedent). `None` gives the default **0.25 × `size_m[0]` = 0.105 m, toward the +x wall**, computed from `_box_kwargs()` and never retyped.
- **Inside the box:** 0.105 + 0.001 = 0.106 < 0.21; 0.001 < 0.045; cap tips at \|z\| = 0.016 < 0.325. The rod is 0.104 m from the +x wall.
- **The placement is a modelling choice and claims nothing.** The operator may set the standard's placement with a cited source.

**Edits allowed.**
1. **`src/fem_em_solver/io/mesh.py`, ruling (c).** Keyword-only, additive, disclosed, ≈ 25 lines.
   - **(a)** Delete the remaining `phantom_box` × `implant_rod` refusal. The rod's own refusals still fire: hole or coat under `as_hole`, and one implant per build.
   - **(b)** `_normalise_implant_rod` gains keyword-only `phantom_box: Optional[Dict] = None`.
     - When it is set, each of the three loops tests both cap centres e with its own radius ρ (r, b, or the tip ball): \|e_x\| + ρ < Lx/2, \|e_y − off\| + ρ < Ly/2 and \|e_z\| + ρ < Lz/2, all strict.
     - The test is exact: the capsule is the convex hull of its end balls, and the box is convex.
     - Refusal message: `"implant_rod is not strictly inside the phantom_box"`.
     - With `None`, the cylinder lines execute unchanged.
   - **(c)** Both call sites pass `phantom_box`.
   - One docstring sentence.
2. **New `tests/validation/test_geo36_rod_f_clinical.py`.**
   - The build test is selected by `FEM_EM_GEO36_ROD_BOX=1` only, never by `-k`. One build per process.
   - It imports by name, never a `test_*` symbol:
     - `_box_kwargs`;
     - `_ring_ports`, `_owned`, `_facet_coords`, `_area`, `LEG_COUNT`, `_refuse_before_gmsh` (`test_geo38_f_clinical_build.py:103, :119, :107, :125, :38, :384`);
     - `ROD`, `_rod_volume`, `_rod_area`, `_interior_facet_check`, `_cell_values_with_ghosts`;
     - `_census_checks` / `_owned_census` from `tests/validation/test_wf12_f_human_hole.py:218, :209`, as C2 imports them (`test_geo38_f_clinical_box.py:78`);
     - `gmsh_stop_armed`, `arm_gmsh_stop_rule` (`tests/validation/gmsh_stop_rule.py`, item 522);
     - `_interface_facet_tags`.
   - Constants: `F_CLINICAL_ROD_CELL_GUARD = 1_400_000`, and C2's conductor mass as a printed comparand (item 522's constant).
3. **The disc module's `test_d10_rod_with_box_still_refused`** (parked `:129`, on `main` once item 522 lands) **turns red by construction under (a).**
   - The licensed, disclosed edit: rename it `test_d10_rod_and_disc_with_box_refused`, and make its body the box with both bodies, refused with `"implant_spheroid cannot be combined with implant_rod"`.
   - Nothing else in that module changes, and W0 re-runs it (ruling 2).

**Anchors (asserted at `-n 4` and `-n 2`; owned, reduced, collectives before any rank-0 branch).** L³ = Lx·Ly·Lz from `_box_kwargs()`; V_rod = `_rod_volume()` = 9.843657e-8 m³; A_rod = `_rod_area()` = 2.010619e-4 m².
- **(G0)** Cells < `F_CLINICAL_ROD_CELL_GUARD` before any facet reading.
- **(G1)** `implant_rod_cad["rod_mass_m3"]` = `["kept_mass_m3"]` = V_rod, and `["wall_area_m2"]` = A_rod, at `STRIP_BAND` (step 6's (F1) / (F4) bands).
- **(G2)** CAD(403) = V_rod, CAD(3) = L³ − V_rod, and CAD(3) + CAD(403) = L³, at `STRIP_BAND`.
- **(G3)** The owned tag-3 + tag-403 tet volumes = L³ at `STRIP_BAND` (planar box, shared interface). The tag-3 extents = C2's (B2) extents.
- **(G4)** Σ `cad_mass_by_group` + `conductor_cad_mass_m3` = `air_domain_cad_mass_m3` at `STRIP_BAND`. This needs item 522's change (c).
- **(G5)** 402: owned count > 0, exterior 0, every facet between one tag-3 cell and one tag-403 cell (`_interior_facet_check`, after an unconditional `create_entity_permutations()`), and meshed dS(402) / A_rod in `TERMINAL_AREA_BAND`.
- **(G6)** The OCC wall centroid lies at the rod centre to `EXACT`·ℓ.
- **(G7)** `_census_checks` is empty; tags 3 and 403 are present; tag 1 is absent; `implant_spheroid_cad` is `None`.
- **(G8)** Per ring sheet: owned 210+i facets = `ring_sheet_gmsh_triangles[i]` on 64 / 64 sheets, and area = chord·w at `STRIP_BAND`. These are C2's (T1) / (T2), owned-facet numpy only.
- **(G9)** Width invariance: every integer of (G0)–(G8) and every CAD float equal at `-n 4` and `-n 2`.
- **(G10)** Refusals (pure, W0):
  - (N3r) below;
  - a coated rod with the box under `as_hole` is refused;
  - a hole rod with the box under `as_hole` is refused;
  - the kept rod is accepted by `_normalise_implant_rod(…, phantom_box=box)`.

**Imported bands:** `STRIP_BAND` 1e-9 (`tests/validation/test_ops61_exact_strip_ports.py:58`); `EXACT` 1e-9 (`tests/mesh/test_birdcage_conductor_hole.py:117`); `TERMINAL_AREA_BAND` (0.95, 1.0) (`tests/mesh/test_birdcage_leg_gaps.py:73`).

**Negative controls (rule (e)).**
- **(N1r) Asserted binary; factor predicted.** C2's comparand L³ misses CAD(3) by V_rod / (L³ − V_rod) = 4.006e-6, 4.0e3× `STRIP_BAND`. That fixed ratio is also the ceiling.
- **(N2r) Asserted by construction; factor predicted.** Item 470's pre-repair conductor reading misses (G4) by V_rod / can = 9.843657e-8 / 4.159880446e-01 = 2.366e-7, 237× `STRIP_BAND` (can volume from `…185056Z_…:11903`).
- **(N3r) The frame misreading, as a planted control. Both outcomes asserted (W0).**
  - The request: "toward a ±y wall by a quarter of the width", read in C2's frame, gives centre (0, 0.105, 0).
  - The cylinder test accepts it: hypot(0, 0.105) + 0.001 = 0.106 < R = 0.214767 m.
  - The box test refuses it: 0.106 > Ly/2 = 0.045, a 61 mm breach.
  - Separation: accept vs refuse, binary.
- **(N4r) Printed only:**
  - tag-403 cells against step 6's 15 025;
  - meshed 403 / V_rod against 0.986763, and meshed 402 / A_rod against 0.996276 (step 6 W2, `20261011T190108Z_GEO-36-step6-W2.log`);
  - the conductor mass against C2's 1.994040484843362e-03, predicted ≤ 1e-9 relative.

**Printed predictions, written first.**
- Cells: 1 149 150 + [50 k, 100 k], point 1 216 304 (+67 154, step 6's F-human kept-rod increment).
- T_rod = T + [−40, +40] s.
- Summed RSS ≈ 2.7 GiB.
- E, F and U₂ = 2E + 2F: the price input of a rod solve in F-clinical.

**Windows.**
- **General:** foreground; Bash `timeout: 660000`; each harness call is the whole Bash command; `-s`; durable capture verbatim.
- **Pre-checks before R1, W1 and W2:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0, then `docker compose -f docker/docker-compose.yml exec -T fem-em-solver cat /proc/loadavg` is journaled.
- **Order:** code → W0 → R1 → W1 → W2 → commit → journal.
- **W1 starts by minute 40 of the slot.** W2 starts only by minute 44. Otherwise commit with W1 as the record and leave item 530 OPEN for "W2 only" in the next eligible slot. That is not an attempt.
- **W0** (pure, `-n 1`, ≈ 15 s):
  `scripts/testing/run_and_log.sh GEO-36-rodbox-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/validation/test_gmsh_stop_rule.py tests/validation/test_geo36_rod_f_clinical.py tests/validation/test_imp6_disc_f_clinical.py tests/mesh/test_birdcage_implant_rod.py::test_a7_refusals_raise_without_meshing tests/mesh/test_birdcage_implant_spheroid.py::test_s9_refusals_raise_without_meshing tests/validation/test_geo38_f_clinical_build.py::test_geo38_shield_refuses_attempt1_table -v -s --tb=short'"`
- **R1** (rule-(c) re-gate of the kept rod under `as_hole` with `phantom_box=None`, the path (b) edits; also the calibrant; complex, `-n 2`, 180–185 s; step 7's W3c relabelled):
  `scripts/testing/run_and_log.sh GEO-36-rodbox-R1-ruleC-step6 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_rod_f_human.py -s -v --tb=short'"`
  - It must reprint, bitwise, the step-6 W2 numbers: 436 541 cells and (F1)–(F9) (`20261011T190108Z_GEO-36-step6-W2.log:6261–6300`).
  - Print h = gmsh / 166.17 s.
  - The F-small rod build is not re-run: no line it executes changes, and its pure (A7) refusals run in W0. This is item 482's precedent for the F-small spheroid.
- **W1** (`-n 4`, armed, the window of record):
  `scripts/testing/run_and_log.sh GEO-36-rodbox-W1-n4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_GEO36_ROD_BOX=1 FEM_EM_GMSH_STOP_RULE=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/validation/test_geo36_rod_f_clinical.py -s -v --tb=short > /workspace/logs/geo36-rodbox-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo36-rodbox-w1-raw.log; cat /workspace/logs/geo36-rodbox-w1-raw.log; exit \$rc'"`
- **W2** (`-n 2`, armed): W1's command with label `GEO-36-rodbox-W2-n2`, raw log `geo36-rodbox-w2-raw.log`, and `-n 2`. Every (G9) integer must equal W1's.

**Traps already paid for.**
- **The frame.** The weekly's "±y wall … a quarter of the box width" uses `MAT-21`'s frame, where x is 650 mm (A), y is 420 mm (B) and z is 90 mm (C) (`tests/validation/test_mat21_box_closed_form.py:5, :34`). In C2's coil frame the same placement is **+x by Lx/4**. (N3r) proves that the box refuses the literal reading.
- **(D10)** goes red by construction, so edit 3 is part of this item, not a fix in-slot.
- **Two `_owned_census` functions.** `test_wf12_f_human_hole.py:209` takes (msh, tags, comm); `test_birdcage_conductor_hole.py:126` takes (msh, tags, dim, comm). Import the first, as C2 does.
- **The default tip ball (2r) is a sizing ball, not CAD.** It is containment-tested only when `tip_ball_radius` is set.
- **Kept rod handling, unchanged:** enter the kept rod once (`GEO-9`); 402 stays the dolfinx-side rebuild (known-issues 9).
- **Do not reuse step 6's centre (0.015, 0, 0).** `IMP-3`'s `TIPS` and ball helpers are pinned to it (`tests/validation/test_imp3_rod_birdcage.py:71, :140`), so they do not carry to this placement. Never import `test_imp3_*`, `test_wf9_*` or `test_wf10_*`.
- **Item 522's traps all hold:**
  - the rule only in W1 / W2, never in R1 or an `xl` command;
  - a stop is the price, not a red;
  - Status 124 without the marker: `pgrep` first, then the sweep;
  - never retry, never fall back;
  - rule (i).

**Scope.** In this item: CAD and mesh identities of the kept 1 mm × 30 mm rod in F-clinical's ASTM box at x = Lx/4, plus the refusals. **Not claimed:** a solve, a PEC condition, a tip ratio, SAR or heating; an ASTM placement claim; a coat or a hole rod; both bodies at once; any refinement knob. `GEO-36` stays ✅ and `GEO-38` stays 🟡.

**Negative result.** Report and stop. Park on `attempt/GEO-36-rodbox-<UTC>`, open a known-issues entry for a gmsh FAIL or a composition red, and mark item 530 BLOCKED in the same commit (rule (d)).
- **A stop at the rule in W1:** the price. Record T > 450 s and the phases. Mark item 530 BLOCKED, to `xl` (the weekly's). Land the code only if R1 and W0 are green (item 494's precedent).
- **A stop in W2 after a green W1:** the `-n 4` record stands. The `-n 2` width goes to the weekly.
- **(G1)–(G3) red with W0 green:** the box cut around the rod is wrong.
- **(G4) at ≈ 2.37e-7:** change (c) does not cover the rod.
- **(G5) exterior facets:** the cavity reached the rod wall.
- **(G9) red:** a rank-local count. Print both widths.
- **R1 not bitwise:** revert `src/` and stop.
- **A printed reading outside its prediction is a reading, never a red.**

No `PREDICTED RED` line.

**Status it can move.** The §7 `GEO-36` row stays ✅. Its Done-when ("the identities asserted in one logged window at `-n 2`, cells and time recorded ⇒ ✅", `PROJECT_PLAN.md:1925`) is met as a record. The row gains, if it stays ≤ 2 KB: "Step 8 (item 530, I33): the kept rod in F-clinical's ASTM box at x = Lx/4 along the 650 mm axis — identities green at `-n 4` / `-n 2`, ⟨cells⟩ cells, gmsh ⟨T_rod⟩ s." There is no Done-when amendment, so rule (j) does not bind. Its E / F / U₂ make a rod-in-F-clinical solve poseable by the weekly (`xl`).

**For the review to rule on.**
1. Renumber as step 8.
2. Edit 3, to item 522's `test_d10_rod_with_box_still_refused`.
3. The +x reading of "±y wall".
4. Whether the `-n 2` width may defer to a second slot.

**Every-item rules** (PROJECT_PLAN §9 "Every item"):
- bare whole harness calls;
- a lost docker socket stops the slot;
- no edit after the last window;
- `timeout: 660000` and `-s`;
- never delete a log;
- BLOCKED in the parking commit;
- a killed window is not a result;
- two heavy windows never at once;
- no slot opens `docs/private/` or `aed_results/`.

**Executor spawn prompt (verbatim).** "Execute §9 item 530 — `GEO-36` rod in C2's box (§10 I33) — per `docs/planning/chunks/GEO-36.md` §9 item 530. Read the brief and any review ruling beneath it first; the rulings hold. Run it only if this session's start line reads 12:00, 13:30 or 15:00 (≤ 5 min late), it is the slot's first item, and item 522's W1 log on `main` prints `[birdcage-mesh] mesh wall time T s` with T ≤ 450 and no `[gmsh-stop-rule] STOPPED` line; otherwise skip it. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s; container `timeout -k 30`. The `src/` change is keyword-only and disclosed: the box × rod refusal deleted, box containment in `_normalise_implant_rod` for all three radii when `phantom_box` is set; every `phantom_box=None` line executes as today. New module `tests/validation/test_geo36_rod_f_clinical.py`, plus the one licensed rename-and-body edit of the disc module's rod-refusal test; nothing else. Run W0, R1, W1 (`-n 4`), W2 (`-n 2`) in order, foreground, on the final files; `FEM_EM_GMSH_STOP_RULE=1` in W1 / W2 only. R1 must read step 6's numbers bitwise or you revert `src/` and stop. Never retry a stop, a kill or a gmsh failure. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-rodbox-*`, mark item 530 BLOCKED in the same commit, and report."

### Review ruling 2026-10-18 03:00 (item 530, I33)

Accepted, with these rulings. **Run it only in the 12:00 / 13:30 / 15:00 slots, and only after item 522's W1 built inside the 450 s rule.** If item 522 stopped at the rule, or has not run, skip this item: it is BLOCKED and goes to the weekly as an `xl` candidate.
- **5. Renumbered: `GEO-36` step 8.** Step 7 is item 427's spheroid. §10 says "step 7"; this is a numbering fix only.
- **7. Placement: +x by Lx/4 = 0.105 m in C2's coil frame.** That frame is the box's own width direction. In the coil frame, y is the 90 mm depth, where a quarter-width offset cannot fit. This is the same placement as the weekly's "toward a side wall by a quarter of the box width" and is named for the weekly to ratify. The planted control, which checks that the literal y reading is refused, is accepted.
- **6. Licensed, disclosed:** renaming the parked module's `test_d10_rod_with_box_still_refused` and editing its body. This step removes that refusal by design, so the test is re-pointed to assert the new admission. No other assertion changes.
- **8. Accepted.** If the `-n 4` window lands green after minute 40, the `-n 2` width-invariance window may go to a second daytime slot. The item stays 🟡 and partial (rule (d)), and the step closes only when both widths have run.
- **4.** The shared stop-rule helper from item 522 is reused.

**Outcome 2026-10-19 (15:00 CDT slot): DONE.** `src/` change as licensed (box × rod refusal deleted; `_normalise_implant_rod(…, phantom_box=)` tests rod, coat and tip ball against the box half-extents; `phantom_box=None` unchanged). W0 (`-n 1`) 31 passed / 2 skipped, 6.4 s, including (N3r) (cylinder margin 0.108767 m, box breach 6.100e-02 m) and the re-pointed `test_d10_rod_and_disc_with_box_refused` (`20261019T200332Z_GEO-36-rodbox-W0.log`). R1 (complex, `-n 2`) 14 passed, 169.9 s; 436 541 cells and (F1)–(F9) bitwise against step 6's W2 (only difference: the `tip_ball_radius` key item 464 added to the echoed rod dict); gmsh 157.98 s, h = 157.98 / 166.17 = 0.951 (`20261019T200350Z_GEO-36-rodbox-R1-ruleC-step6.log`). W1 (`-n 4`, armed) 5 passed, 447.1 s: 1 217 025 cells (+67 875 over C2; prediction 1 216 304), gmsh 418.82 s under the 450 s rule (T_rod − T = −22.3 s, inside [−40, +40]); (G1) 4.4e-16 / 0.0; (G2) ≤ 4.4e-16; (G3) 2.2e-16, extents exact; (G4) 6.457e-11; (G5) 4 954 interior 402 facets, 0 exterior, dS/A 0.996243; (G6) 6.6e-20 m; (G7) clean; (G8) 64 / 64, worst area 2.4e-15; (N1r) 4.006372e-06 (4.0e3×); (N2r) 2.365686e-07 (237×); summed RSS 2.688 GiB (`20261019T200702Z_GEO-36-rodbox-W1-n4.log`). W2 (`-n 2`, armed) 5 passed, 444.9 s, gmsh 412.45 s, RSS 2.079 GiB; (G9) integers (1 217 025, 4 954, 0, 4 954, 4 954, 85 122, 14 834, 2 882, E 1 606 744, F 2 556 953) and every CAD float identical to W1 (`20261019T201444Z_GEO-36-rodbox-W2-n2.log`). U₂ = 2E + 2F = 8 327 394 — the price input for a rod-in-F-clinical solve (`xl`, the weekly's). Printed readings outside prediction (readings, not reds): conductor vs C2 3.498e-09 relative (predicted ≤ 1e-9); tag-403 cells 14 834 vs step 6's 15 025; meshed 403 / V 0.986659.

### §9 item 560 — 2026-10-20 (§10 I35, lane 1 — the I34 module's landing step; build-only)

*Drafted read-only for the 2026-10-20 03:00 daily review; the review rules on it.*

**Licence.** §10 I35 (`PROJECT_PLAN.md:4855–4860`, the 2026-10-19 weekly): "a column knob (`FEM_EM_GEO36_S8_COLUMN`) on `GEO-36` step 8's module in `test_geo38_f_clinical_column.py`'s shape — one degree-1 copper column at the feed sheet with the rod in, `PRICE` printed before the post-column asserts — build-only green at `-n 2` asserting I33's 1 217 025 cells; entry 45's prerequisite (i). Serial on nothing." It is `docs/testing/xl-pending.md` entry 45's prerequisite (i) (`:2057`). **The implementer marks (i) MET in the same commit and does not write the command.** The 10-22 03:00 review writes it from entry 40's measured column. No band, tolerance, record or Done-when is set or moved; every band is imported.

**Executor and cost.**
- **Executor:** implementer.
- **Tier: heavy by ceiling, not the weekly's "standard".** Tier is execution detail. The build-only proof rebuilds I33's mesh in memory: item 530 measured 447.1 s at `-n 4` and 444.9 s at `-n 2` (gmsh 418.82 / 412.45 s; `20261019T200702Z_GEO-36-rodbox-W1-n4.log`, `20261019T201444Z_GEO-36-rodbox-W2-n2.log`). That cannot fit standard's 180 s.
  - The container ceiling is `timeout -k 30 590`, not 1200. A foreground window must return its footer inside the Bash tool's 660 s, so 590 + 30 is the most a foreground window can carry (rubric element 4; item 530's W1 / W2 used the same).
- **Ranks:** W0 `-n 1`; W1 `-n 2`. Both run on the complex image, because the context builds a `TimeHarmonicSolver` and the sheet widths compile a `dS` area form (C3a's W1 ran the same way).
- **Daytime only.** This is a composed clinical-scale build (the 10-17 epitaph; the 10-19 weekly's daytime rule, `PROJECT_PLAN.md:4758–4762`).
  - It runs only in the 12:00 / 13:30 / 15:00 slots, as the slot's first item, one daytime item per slot. There is no take-next into it or out of it. Any other slot skips it at no attempt cost.
  - It never shares a slot with items 542, 553, 561 or 566, or any other clinical build. Two heavy windows never run at once.
  - **Deadline:** the 10-22 03:00 review reads entry 45's prerequisites. If both 10-20 and 10-21 daytime slots pass without this item, entry 45 stays PENDING (its own text, `:2057`). Run it ahead of the daytime examples.
- **Class:** chain no-solve (the weekly priced I35 at 6.6, `PROJECT_PLAN.md:4933`; item 519's ruling (d) re-classed a build-only clinical proof the same way). **Cost:** ≈ 6.6 measured / **≈ 45 predicted** slot-minutes. That is code ≈ 18 (≈ 150 lines added to one module, no `src/`), W0 ≈ 1, pre-checks ≈ 2, W1 ≈ 8, and commit and journal ≈ 8.
  - **If W1 has not started by minute 42:** commit the module with W0 green and mark item 560 **OPEN for "W1 only"**. The next daytime slot runs W1 on the committed file, unchanged (rule (i)). That is not an attempt.
- **Gmsh stop rule: not armed.** This build is priced (I33 built 412–419 s of gmsh under the rule, twice). The 10-19 weekly keeps the epitaph for the *next unpriced* build (`PROJECT_PLAN.md:4758–4762`). The new column path never calls `arm_gmsh_stop_rule`. `FEM_EM_GMSH_STOP_RULE` is never set in W1 or in entry 45's command (item 522's trap).
- **Family:** `GEO-36` step 8 column knob (I35), attempt 1 of 4.

**What exists (`main`, 2026-10-20).**
- **The step-8 module**, `tests/validation/test_geo36_rod_f_clinical.py` (391 lines, item 530, `c07031df`):
  - `ROD_BOX_ENV` (`:65`); `_rod_centre()` / `_kwargs()` (`:88–97`), which is `_box_kwargs()` plus the kept rod at (Lx/4, off, 0);
  - four pure refusal tests (`:103–159`);
  - the build test `test_geo36_step8_rod_in_astm_box` (`:165–391`), whose (G1) / (G2) / (G4) / (G7) are diag- or census-only reads (`:215–226`, `:265–270`, `:283`, `:383–386`).
- **Entry 40's column module**, `tests/validation/test_geo38_f_clinical_column.py`. **Never edit it: entry 40 runs it on 10-21.**
  - `DRIVEN = "P33"`, `HELD = "P34"` (`:53–54`); `_selector()` reads `FEM_EM_GEO38_STEP1B` at call time (`:72–73`).
  - `_context_from_mesh` (`:95–138`): its `material_map` maps tag 3 only (`:129–130`).
  - `_build_f_clinical(comm)` (`:141–154`) returns `{"ctx", "cells", "facets", "shield_radius_m", "outer_census"}`.
  - The column body `test_geo38_step1b_f_clinical_copper_column` (`:327–519`) calls `_build_f_clinical` as a module global (`:369`) and reads `F_CLINICAL_R1_CELL_RECORD` / `F_CLINICAL_R1_UNKNOWNS_RECORD` by a call-time `from … import` (`:348–354`).
  - Its order is (A0) / (iv) / (S) asserted (`:396–404`), then the factorising column and one held drive (`:477–478`), `PRICE` and `PROJECTIONS` (`:494–504`), then (P) / P_surf / NEG / (R1) asserted (`:515–519`). **`PRICE` already prints before the post-column asserts.**
- **C3a's helpers** (`tests/validation/test_geo38_f_clinical_matrix.py`; never edit it, since entry 44 runs it): `_census_failures` (`:148–161`) and `_sheet_geometry` (`:164–185`).
- **The rod's material precedent:** `IMP-3` step 8, the kept rod under `as_hole` on copper F-human. Tag 403 is mapped to `SIGMA_WIRE_S_PER_M` (`tests/validation/test_port_gap_voltage_impedance.py:235`, 800 S/m) with the gel's ε_r and μ_r, through `dataclasses.replace` on the problem (`tests/validation/test_imp3_f_human_tuned_rod.py:388`, `:439–446`). Its (M0) reads σ on owned tag-403 cells through `build_material_fields` (`:458–476`).
- **I33's measured numbers** (`20261019T201444Z_GEO-36-rodbox-W2-n2.log`; item 530 outcome): 1 217 025 cells; E 1 606 744; F 2 556 953; U₂ = 2E + 2F = 8 327 394; summed RSS 2.079 GiB.

**Code: one module edited, additive only.** `tests/validation/test_geo36_rod_f_clinical.py` gains constants, helpers and tests. **No existing line changes**: `git diff` on the module shows insertions only. Edit nothing else, so no `src/`, no column or matrix module, no example. The one exception is entry 45's status line, after green.
- **Constants (records imported from I33's log, each commented with it):**
  - `S8_COLUMN_ENV = "FEM_EM_GEO36_S8_COLUMN"`;
  - `F_CLINICAL_ROD_CELL_RECORD = 1_217_025`;
  - `F_CLINICAL_ROD_UNKNOWNS_RECORD = (1_606_744, 8_327_394)`.
  - These are I33's measured values, registered as the column's (A0) / PRICE comparands under the (1\*) licence, not new targets.
- **`_s8_column_mode()`:** unset or empty ⇒ `None`; `build` or `1` ⇒ itself; anything else ⇒ `ValueError` (C3a's `_mode`, `:75–81`).
- **The rod's material:** `IMP-3` step 8's convention. Tag 403 is mapped to `HomogeneousMaterial(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r)`, where `gel` is the context problem's tag-3 entry. **This is a modelling choice for the review to confirm (ruling 1).** Without it, tag 403 takes the air default: an empty rod-shaped bubble in the gel, which the (N-M) control detects.
- **`_build_rod_f_clinical(comm)`** returns `_build_f_clinical`'s dict shape (plus an extra `"problem_gel_only"` key the column body ignores). It runs, in order:
  1. `MeshGenerator.birdcage_port_domain(comm=comm, return_diagnostics=True, **_kwargs())`, timed and called through the class attribute (imported inside the function).
  2. `chord` from `diag["ring_port_layout"]`, then `ctx, outer_census, _w = col._context_from_mesh(msh, cell_tags, facets, chord, comm)`.
  3. The problem is rebuilt with the rod material via `dataclasses.replace`. `ctx["solver"]` is replaced by `TimeHarmonicSolver(problem_rod, degree=DEGREE)` after `del ctx["solver"]`; never hold two solvers (`WF-12` step 4's trap).
  4. The rod and fixture asserts (below), then the return.
- **`test_geo36_i35_rod_column(monkeypatch)`**, `@complex_only`:
  - mode `None` ⇒ skip;
  - `build` ⇒ `_build_rod_f_clinical`, then the column body's pre-column asserts (A0) / (iv) / outer census / (S), copied from `:386–404` with `F_CLINICAL_ROD_CELL_RECORD`; it prints `BUILD-ONLY return` (t_build, RSS per rank and summed, window) and returns before the Leontovich term or any assembly or factorisation;
  - `1` (entry 45 only) ⇒ `_route_column(monkeypatch)`, then `col.test_geo38_step1b_f_clinical_copper_column()` is called directly (`col` is the module object).
- **`_route_column(monkeypatch)` — the "column-shape" route. It runs entry 40's body unchanged on the rod fixture: C8's no-edit route, which item 518's ruling (g) accepted.** It applies four patches:
  - `monkeypatch.setenv("FEM_EM_GEO38_STEP1B", "1")`;
  - `monkeypatch.setattr(col, "_build_f_clinical", _build_rod_f_clinical)`;
  - `monkeypatch.setattr(s1a, "F_CLINICAL_R1_CELL_RECORD", F_CLINICAL_ROD_CELL_RECORD)`;
  - `monkeypatch.setattr(s1a, "F_CLINICAL_R1_UNKNOWNS_RECORD", F_CLINICAL_ROD_UNKNOWNS_RECORD)`.
  - `s1a` is `tests.validation.test_geo38_f_clinical_build`. The driven sheet is entry 40's P33 and the held drive P34, so entry 45's price is entry 40's shape on one more body. PRICE prints before the post-column asserts by the body's own order.
  - **Alternative, for the review:** a copied column body (≈ 110 lines) in this module. It needs no patching, but it can drift from entry 40's shape. The drafter prefers the patch route (ruling 2).
- **Coding tests** skip whenever `S8_COLUMN_ENV` is set (C3a's `_coding` marker, `:481–482`): (K0), (K1), (N-A0) below.
- **Imports:** module objects only (`import tests.validation.test_geo38_f_clinical_column as col`, `… _build as s1a`), inside function bodies. Never `from … import test_*`: pytest would collect entry 40's column test into this module.

**Anchors (asserted; owned entities, reduced; collectives before any rank-0 branch).**
- **(A0)** Cells = `F_CLINICAL_ROD_CELL_RECORD` within `CELL_COUNT_BAND` 0.01 (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`), the column body's form (`test_geo38_f_clinical_column.py:400`). Predicted rel **0 exactly**: gmsh is serial on rank 0, and the same kwargs gave 1 217 025 at `-n 4` and `-n 2`.
- **(iv), outer census, (S):** as the column body (`:401–404`). `_census_checks` (`tests/validation/test_wf12_f_human_hole.py:218`) was clean on this mesh at I33's (G7).
- **(C)** `c3._census_failures(ids, owned_counts)` is empty. There are 64 specs P33…P96, each with more than 0 owned facets, and 64 = 2·`LEG_COUNT`. Owned counts come from `c3._sheet_geometry(ctx, comm)`.
- **(D-rod)** I33's diag and census identities, re-read on this build at their imported bands. The weekly's "(D-rod)" label (entry 45, `:2063`) is mapped here to existing identities, not new ones:
  - (G1) `rod_mass_m3` = `kept_mass_m3` = V_rod, and `wall_area_m2` = A_rod, ≤ `STRIP_BAND` 1e-9 (`tests/validation/test_ops61_exact_strip_ports.py:58`);
  - (G2) CAD(403) = V_rod, CAD(3) = L³ − V_rod, CAD(3) + CAD(403) = L³;
  - (G4) Σ `cad_mass_by_group` + `conductor_cad_mass_m3` = `air_domain_cad_mass_m3`;
  - (G7) tags 3 and 403 present, tag 1 absent, `implant_spheroid_cad` is `None`.
  - V_rod and A_rod come from `_rod_volume()` / `_rod_area()` (`tests/mesh/test_birdcage_implant_rod.py:113`, `:117`). (G5), (G8) and (G3) are not repeated: they need facet or tet work, and they are I33's record.
- **(M) The rod's material in the solver's problem.** `build_material_fields` (`src/fem_em_solver/core/time_harmonic.py:219`) on the rebuilt problem. Owned tag-403 cells (via `_cell_values_with_ghosts`, `tests/mesh/test_birdcage_implant_rod.py:171`) read σ min = max = `SIGMA_WIRE_S_PER_M` and ε_r = the gel's. Owned tag-3 cells read σ = `SALINE_SIGMA` (`tests/validation/test_lossy_sphere_fullwave.py:93`). All exact, as `IMP-3` step 8's (M0) asserted.
- **(U1)** `ctx["solver"].function_space()` global size × bs = `F_CLINICAL_ROD_UNKNOWNS_RECORD[0]` = E = 1 606 744, exactly. This is C1 (a)'s (U) identity at degree 1.
- **(K0)** Knob parsing: unset, `build` and `1` accepted; `2` and `x` raise `ValueError`.
- **(K1) Routing tripwire (no gmsh).** Inside `monkeypatch.context()`:
  - Apply `_route_column`, and replace `MeshGenerator.birdcage_port_domain` with a recorder that raises a sentinel `_Routed(kwargs)`.
  - Call `col.test_geo38_step1b_f_clinical_copper_column()`. It must raise `_Routed`, and the recorded kwargs must carry `implant_rod` with `keep_cells` True and centre `_rod_centre()`, plus `phantom_box`.
  - On exit, `col._build_f_clinical` is the original, `s1a.F_CLINICAL_R1_CELL_RECORD == 1_147_366`, and `FEM_EM_GEO38_STEP1B` is unset. This proves entry 40's module is untouched in-process.

**Negative controls (rule (e)).**
- **(N-route) Asserted (K1's twin).** The same call with only `FEM_EM_GEO38_STEP1B=1` set, without `_route_column`, reaches the recorder with kwargs that have **no** `implant_rod`: it is the 1.5 T fixture. Binary, so the patch is load-bearing.
- **(N-A0) Asserted by construction.** The 1.5 T record misses the rod mesh: |1 217 025 / 1 147 366 − 1| = **6.07e-2**, which is 6.1× `CELL_COUNT_BAND`. That is the ceiling: a fixed ratio of two records. If the record patch silently failed in entry 45's window, (A0) would go red there before the column. This is a numpy-only test.
- **(N-M) Asserted, W1.** `build_material_fields` on `"problem_gel_only"`, the context's unpatched problem, reads σ = **0.0** on every owned tag-403 cell: the air default, an empty rod. The separation is 800 S/m against 0, binary.

**Predictions, written first (printed, never asserted).**
- (A0) rel 0;
- build ≈ 445 s [410, 560] at a 1-min loadavg ≤ 4 (I33: 412–419 s of gmsh at loadavg 1.94–3.81);
- context, (M) and (C) ≈ 10–20 s;
- window ≈ 465 s;
- summed RSS ≈ 2.2 GiB;
- (D-rod) at I33's levels: (G1) ≤ 4.4e-16, (G2) ≤ 4.4e-16, (G4) 6.457e-11.
- **For entry 45's review (printed in the journal, not this item's):** ρ = 1 217 025 / 369 387 = 3.2948. The body's "conservative" column 39.2·(1 606 744 / 500 568)² ≈ 404 s, against entry 40's 364 s.

**Windows.** Each runs in the foreground with Bash `timeout: 660000` and pytest `-s`, never `-k`. **Every harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.** Before W1:
1. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0;
2. then `docker compose -f docker/docker-compose.yml exec -T fem-em-solver cat /proc/loadavg`, journaled beside W1's build time.

- **W0** (complex, `-n 1`, ≈ 30 s). Expected: `test_complex_mode` green; the rod module's four step-8 refusal tests and (K0), (K1), (N-route), (N-A0) pass; the step-8 build test and the gated I35 test skip. The column module reads **0 passed, 2 skipped**, as at item 494 W1 and C3a W0, which proves entry 40's module is collectable and untouched.
  `scripts/testing/run_and_log.sh GEO-36-I35-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && unset FEM_EM_GEO36_S8_COLUMN FEM_EM_GEO36_ROD_BOX FEM_EM_GEO38_STEP1B FEM_EM_GMSH_STOP_RULE && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical.py tests/validation/test_geo38_f_clinical_column.py -v -s --tb=short'"`
- **W1** (build-only, complex, `-n 2`, the window of record). Expected: the gated I35 test and the four refusal tests pass; the step-8 build and the coding tests skip.
  `scripts/testing/run_and_log.sh GEO-36-I35-W1-build-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_GEO36_S8_COLUMN=build timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical.py -s -v --tb=short > /workspace/logs/geo36-i35-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo36-i35-w1-raw.log; cat /workspace/logs/geo36-i35-w1-raw.log; exit \$rc'"`
- **Rule (i):** no module edit after W0 starts. Any edit means W0 and W1 again.

**Traps already paid for.**
- **Never edit `test_geo38_f_clinical_column.py` or `test_geo38_f_clinical_matrix.py`.** They are entry 40's module (10-21) and entry 44's (10-24).
- **Never put `FEM_EM_GEO38_STEP1B`, `FEM_EM_GEO36_ROD_BOX` or `FEM_EM_GMSH_STOP_RULE` on a command line.** With the column module on the path, `STEP1B=1` would build the 1.5 T r1 *and* the rod mesh in one 590 s window, a certain kill. `ROD_BOX=1` would add a third build. `STEP1B` is set only by `monkeypatch` inside the I35 route, and teardown restores it.
- **Import module objects, never `test_*`.**
- **Tag 403 unmapped is a silent air bubble.** `_context_from_mesh` maps tag 3 only. (M) and (N-M) catch it.
- **One solver.** `del ctx["solver"]` before building the rod solver (`WF-12` step 4).
- **Rank safety:** owned cells and facets only. Use `_cell_values_with_ghosts` for tag reads on the ghosted mesh. `cell_tags.values` and `assemble_scalar` are rank-local. Never use UFL `dS` over the sheets for counts.
- **The sheet-width `dS` form compiles in W1.** If W1 is killed mid-compile, run `pgrep -c python3` first, then the rubric-4 sweep as `GEO-36-I35-cache-sweep`. Delete orphan stems only.
- **`pgrep -f "python3 -m pytest"` is denied by the bash guard;** use `pgrep -c python3`.
- **Durable capture** ends in `; exit \$rc`. No pipe; never delete a log.
- **A killed or backgrounded window is not a result.** Never background a harness call.

**Scope.** The column knob, its routing, and a build-only proof at `-n 2` on I33's mesh with the rod's material wired.
- **Not claimed:** any column, S entry, (P), (R1) or price; entry 45's command, tier or brackets; a PEC or insulated rod; a tip ratio, SAR or heating; an ASTM placement. `GEO-36` stays ✅. `GEO-38` and `IMP-3` are unchanged.

**Negative result.** Report and stop. Park on `attempt/GEO-36-I35-<UTC>` and mark item 560 BLOCKED in the same commit (rule (d)). Entry 45's (i) stays open.
- **(A0) red, or cells ≠ 1 217 025:** I33's mesh moved on this image. Open a known-issues entry, because entry 45 reads (A0) against the same record.
- **(D-rod) red with (A0) green:** the rod's CAD moved. BLOCKED with the reading.
- **(M) red:** the material map did not reach tag 403. BLOCKED; never re-map in-slot.
- **(C) or (U1) red:** the port layout or the dof count is not I33's or C1 (a)'s identity.
- **W1 killed at 590 s:** run `pgrep`, journal the loadavg and the last `[birdcage-mesh]` line, and park. The review may rule that entry 45's own window serves as the build proof, since the full route asserts (A0), (D-rod), (M), (C) and (U1) before its column (C3a's fallback, `GEO-38.md:874`). That is the review's call, never a retry in-slot.
- **A printed reading outside its prediction is a reading, never a red.**

No `PREDICTED RED:` line.

**Status it can move.**
- On green, in `docs/testing/xl-pending.md` entry 45, mark (i) **MET 2026-10-2x (§9 item 560, `<commit>`; `<W0 log>`, `<W1 log>`; module `tests/validation/test_geo36_rod_f_clinical.py`, knob `FEM_EM_GEO36_S8_COLUMN` = `build` / `1`)**. Edit the Status line only, never the command (step 6b rule 5).
- Entry 45 then stays PENDING on (ii), entry 40's 10-21 price.
- §10 I35 gains a DONE line.
- The `GEO-36` §7 row stays ✅, with no Done-when change. If it stays ≤ 2 KB, it may gain "I35 (item 560): column knob landed, build-only green at `-n 2`".

**For the review to rule on.**
1. The rod's material in the column: `IMP-3` step 8's σ = 800 S/m kept rod with the gel's ε_r (the drafter's reading of "with the rod in"). A PEC rod is not available on the kept-under-`as_hole` route.
2. The patch route (entry 40's body run unchanged on the rod fixture) against a copied body.
3. The tier: heavy by ceiling with container 590 s, against the weekly's "standard".
4. The ordering: I35 ahead of the daytime examples (542 / 553), because of entry 45's 10-22 03:00 deadline.

**Every-item rules** (PROJECT_PLAN §9 "Every item"):
- bare whole harness calls;
- a lost docker socket stops the slot;
- no edit after the last window;
- `timeout: 660000` and `-s`;
- never delete a log;
- BLOCKED in the parking commit;
- a killed window is not a result;
- two heavy windows never at once;
- no slot opens `docs/private/` or `aed_results/`.

**Execution order.** Read the brief and any ruling beneath it → code → W0 → pre-checks → W1 → commit (the module, harness logs, entry 45's status edit, the §9 DONE line, the §10 DONE line) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 560 — `GEO-36` I35, the column knob on step 8's module — per `docs/planning/chunks/GEO-36.md` §9 item 560. Read the brief and any review ruling beneath it first; the rulings hold. Run it only if this session's start line reads 12:00, 13:30 or 15:00 (≤ 5 min late) and it is the slot's first item; otherwise skip it. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s; container `timeout -k 30`. Edit only `tests/validation/test_geo36_rod_f_clinical.py`, additively (no existing line changes), and, after a green W1, entry 45's Status line in `docs/testing/xl-pending.md` — never its command. Never edit `tests/validation/test_geo38_f_clinical_column.py` or `tests/validation/test_geo38_f_clinical_matrix.py`. Never set `FEM_EM_GEO38_STEP1B`, `FEM_EM_GEO36_ROD_BOX` or `FEM_EM_GMSH_STOP_RULE` on a command line. Before W1, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Run W0 then W1 in the foreground. If W1 has not started by minute 42, commit the module with W0 green and mark the item OPEN for 'W1 only'. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-I35-*`, mark item 560 BLOCKED in the same commit, and report."

### Review ruling 2026-10-20 03:00 (item 560)

Binds over the brief. Tier is heavy by measurement (item 530's builds ≈ 445 s at `-n 2`), container `timeout -k 30 590`; daytime, slot-first. The rod's material on the column route is σ = 800 S/m with the gel's ε_r (`IMP-3` step 8's precedent; a PEC rod is not available on the kept route) — provisional, the weekly ratifies with entry 45. Entry 45's undefined "(D-rod)" is read as I33's diag and census identities (G1) / (G2) / (G4) / (G7) at their imported bands — provisional, the 10-22 03:00 review writes entry 45's command on that reading. No gmsh stop rule armed: a priced rebuild.

**Outcome (2026-10-20 12:00 slot): BLOCKED, W1 killed at 590 s; parked on `attempt/GEO-36-I35-20261020T171500Z`.** The module landed as written (302 insertions, 0 deletions on `test_geo36_rod_f_clinical.py`). W0 green at `-n 1` complex: 12 passed, 4 skipped; (N-A0) 6.071e-02 = 6.07× `CELL_COUNT_BAND`; (K1) routed kwargs carry the kept rod at (0.105, 0, 0) with the box; (N-route) reaches the 1.5 T kwargs with no rod; the column module reads 0 passed, 2 skipped (`20261020T170308Z_GEO-36-I35-W0.log`). Pre-checks: `pgrep -c python3` 0; loadavg **11.96 12.05 10.71**, against I33's 1.94–3.81. W1 (`-n 2`, `build`) was killed by the container timeout at 590 s inside gmsh's second optimisation pass. The last `[birdcage-mesh]` line is the implant-rod sizing line (`20261020T170326Z_GEO-36-I35-W1-build-n2.log:81`). No `[GEO-36 I35]` line was printed, so this is not a result. Gmsh phases against I33's W2 (`20261019T201444Z_GEO-36-rodbox-W2-n2.log`): 1D 163.7 / 149.3 s, 2D 30.9 / 22.0 s, 3D 110.2 / 56.3 s, first optimise 13.5 / 7.4 s. The second pass (I33: 171.2 s) had reached ≈ 1/6 of I33's line count at the kill. The build is load-sensitive: under loadavg ≈ 12, I33's ≈ 445 s exceeds 590. `pgrep -c python3` read 0 after the kill. **Anomaly, not diagnosed:** at log `:71`, one rank's interleaved stream shows `test_geo36_i35_rod_column PASSED` while the other rank was still in gmsh (`:76`). No I35 output line exists, so no reading rests on it. It is for the log-pathologist if W1 is re-run. **For the review:** re-run W1 on the parked module, unchanged, in a daytime slot with pre-check loadavg ≤ 4, or let entry 45's window serve as the build proof (the brief's killed-at-590 clause). Entry 45's (i) stays open, and its 10-22 03:00 deadline stands.

### Review ruling 2026-10-20 18:00 (item 560 BLOCKED — diagnosis and re-pose as item 575)

Binds over item 560's brief and the 03:00 ruling where they differ.

1. **Diagnosis: load, not geometry.** Every gmsh phase of the killed W1 ran 1.1–2× slower than I33's W2 (1D 163.7 / 149.3 s, 3D 110.2 / 56.3 s; `20261020T170326Z_GEO-36-I35-W1-build-n2.log` against `20261019T201444Z_GEO-36-rodbox-W2-n2.log`) at a pre-check loadavg of 11.96 against I33's 1.94–3.81. Two same-class composed builds landed the same afternoon at low load: C4a 335.5 s at loadavg 0.18 (`20261020T183431Z_GEO-38-C4a-W1-build.log`) and C11 316.2 s at 0.51 (`20261020T200157Z_GEO-38-C11-W1-build.log`). Nothing in the parked module's build path is new geometry: the build calls `birdcage_port_domain` with I33's kwargs (`attempt/GEO-36-I35-20261020T171500Z:tests/validation/test_geo36_rod_f_clinical.py:563–611`).
2. **The `:71` anomaly, read (not banked as a finding).** The `PASSED` token on `:71` is most simply rank 0's (N3r) result appended to rank 1's unterminated `test_geo36_i35_rod_column ` header by mpiexec's chunked forwarding. Both ranks printed the (N3r) header (`:62–63`) but only one `PASSED` follows it (`:66`); the second appears at `:71`. The same split already shows in the environment tests (a lone `PASSED` at `:54`, doubled `PASSEDPASSED` at `:55`, `:58`, `:61`). Rank 1's I35 line carries no output of its own, and the test body has no rank-dependent early return (`:563–611` above). No status rests on this reading, so no pathologist ruling is taken; item 575 journals whether the pattern recurs, and if rank 1 ever reports the I35 test finished before rank 0's `[GEO-36 I35]` lines, that is a red to report, not to diagnose in-slot.
3. **Re-pose: item 575 (I35 attempt 2 of 4 in the `GEO-36` step-8 column-knob family).** Same target, same module, same commands; the only new clause is a load gate on W1.

### §9 item 575 — 2026-10-20 (§10 I35, re-pose of item 560 — the column knob's W1 at low load; lane 1, build-only)

*Written by the 2026-10-20 18:00 daily review (rule (6) re-pose; target unchanged).*

**Everything in item 560's brief binds** (anchors (A0) / (iv) / outer census / (S) / (C) / (D-rod) / (M) / (U1), controls (N-route) / (N-A0) / (N-M), predictions, traps, scope, negative result, status it can move), with the 03:00 ruling and the 18:00 ruling above. This brief changes only the following.

- **Code: take the parked module by path.** `git checkout attempt/GEO-36-I35-20261020T171500Z -- tests/validation/test_geo36_rod_f_clinical.py`, then `git diff --stat HEAD -- tests/validation/test_geo36_rod_f_clinical.py` must read **302 insertions, 0 deletions**. No other file is taken from the branch (it carries no record file you need; never take `test-results.md`, `attempts.md` or `PROJECT_PLAN.md` from it). No edit to the module (rule (i)).
- **Windows** (log IDs carry `b`; commands otherwise item 560's, verbatim):
  - **W0b** — `scripts/testing/run_and_log.sh GEO-36-I35-W0b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && unset FEM_EM_GEO36_S8_COLUMN FEM_EM_GEO36_ROD_BOX FEM_EM_GEO38_STEP1B FEM_EM_GMSH_STOP_RULE && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical.py tests/validation/test_geo38_f_clinical_column.py -v -s --tb=short'"` — expected as item 560's W0 (12 passed, 4 skipped; column module 0 passed, 2 skipped).
  - **Pre-checks, in order, each its own bare call:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0; then `docker compose -f docker/docker-compose.yml exec -T fem-em-solver cat /proc/loadavg`. **Load gate: the 1-minute figure must be ≤ 4.0.**
  - **W1b** — `scripts/testing/run_and_log.sh GEO-36-I35-W1b-build-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_GEO36_S8_COLUMN=build timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical.py -s -v --tb=short > /workspace/logs/geo36-i35-w1b-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo36-i35-w1b-raw.log; cat /workspace/logs/geo36-i35-w1b-raw.log; exit \$rc'"` — Bash `timeout: 660000`, foreground.
- **If the load gate fails (1-min loadavg > 4.0):** do not run W1b. Commit the module with W0b green, mark item 575 **OPEN for "W1b only"** with the loadavg reading in the line, and journal it; the next daytime slot runs the pre-checks and W1b on the committed module, unchanged. That is not an attempt (item 560's minute-42 clause). The slot then takes no other daytime item.
- **If W1b is killed at 590 s with the load gate met:** the load reading is refuted. Run `pgrep -c python3`, journal the loadavg, each gmsh phase's time beside I33's and the last `[birdcage-mesh]` line, park on `attempt/GEO-36-I35-<UTC>` and mark 575 BLOCKED (rule (d)). Do not retry. Whether entry 45's own window serves as the build proof is the 10-22 03:00 review's call.
- **On green:** commit the module (if not already on `main`), W0b / W1b logs, entry 45's Status line (i) **MET** (item 560's wording with item 575 and the `b` logs; the Status line only, never the command), §9 item 575 DONE (one line), §10 I35 DONE line; then `git branch -D attempt/GEO-36-I35-20261020T171500Z` (its content is fully on `main`). Journal the loadavg, t_build and whether the `:71` pattern recurred.
- **Cost:** chain no-solve, ≈ 7.8 measured / ≈ 30 predicted slot-minutes (no code to write; W0b ≈ 1, pre-checks ≈ 2, W1b ≈ 8 at I33's ≈ 445 s, commit and journal ≈ 8).
- **Daytime, slot-first,** ahead of 542 / 553 (entry 45's 10-22 03:00 deadline): 12:00 / 13:30 / 15:00 only; any other slot skips it at no attempt cost. Never beside another heavy window.
- **No `PREDICTED RED:` line.**

**Outcome (2026-10-21 12:00 CDT slot): ✅ green; entry 45's (i) MET.** Module taken by path from `attempt/GEO-36-I35-20261020T171500Z` (302 insertions, 0 deletions; no edit). W0b `-n 1` complex: 12 passed, 4 skipped; (N-A0) 6.071e-02; (K1) / (N-route) as item 560's W0; column module 0 passed, 2 skipped (`20261021T170043Z_GEO-36-I35-W0b.log`). Pre-checks: `pgrep -c python3` 0; loadavg **0.36 0.85 2.13** (gate met). W1b `-n 2`, `build`: 9 passed, 5 skipped, rc 0, window 465.6 s, harness 468 s (`20261021T170059Z_GEO-36-I35-W1b-build-n2.log`). Readings: t_build 457.5 s (PREDICTED ≈ 445 [410, 560]); gmsh 426.92 s — 1D 152.8 / 2D 22.3 / 3D 59.2 / optimise 7.3 + 178.8 s against I33's W2 149.3 / 22.0 / 56.3 / 7.4 + 171.2 s; (A0) 1 217 025, rel 0; (D-rod) (G1) 4.441e-16, (G2) 4.441e-16, (G4) 6.457e-11, (G7) census clean, tag 403 14 834 cells; (M) tag 403 σ 800 / ε_r 78 exact, tag 3 σ 0.5 exact; (N-M) σ 0.0 on tag 403; (C) 64 specs, owned facets 44–46; (U1) 1 606 744 = E exact; (S) cavity max hypot 0.361 m < 0.3725 m; summed RSS 2.089 GiB (1.60 / 0.48). Printed for entry 45: ρ 3.2947; conservative column 404 s. **`:71` pattern:** the chunked-forwarding split recurred in its benign form (both ranks' I35 headers at `:75–76`; one `PASSED` at `:11994` after the `[GEO-36 I35]` lines at `:11986–11993`, the other appended as `SKIPPEDPASSED` at `:11998`). No `PASSED` preceded the I35 output, so the 18:00 ruling's red condition did not occur. Branch `attempt/GEO-36-I35-20261020T171500Z` deleted.

**Executor spawn prompt (verbatim).** "Execute §9 item 575 — `GEO-36` I35's W1 re-run at low load — per `docs/planning/chunks/GEO-36.md` §9 item 575, which defers to item 560's brief, the 2026-10-20 03:00 ruling and the 2026-10-20 18:00 ruling; read all four first. Run it only if this session's start line reads 12:00, 13:30 or 15:00 (≤ 5 min late) and it is the slot's first item; otherwise skip it. Take `tests/validation/test_geo36_rod_f_clinical.py` from `attempt/GEO-36-I35-20261020T171500Z` by path and check 302 insertions / 0 deletions; edit nothing. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s; container `timeout -k 30`. Run W0b; then `pgrep -c python3` (must read 0) and `cat /proc/loadavg` as separate bare `docker compose -f docker/docker-compose.yml exec -T fem-em-solver …` calls; run W1b only if the 1-minute loadavg is ≤ 4.0, otherwise commit the module with W0b green and mark the item OPEN for 'W1b only'. Never set `FEM_EM_GEO38_STEP1B`, `FEM_EM_GEO36_ROD_BOX` or `FEM_EM_GMSH_STOP_RULE` on a command line; never edit `tests/validation/test_geo38_f_clinical_column.py` or `tests/validation/test_geo38_f_clinical_matrix.py`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a kill at 590 s with the gate met, or any red, park on `attempt/GEO-36-I35-*`, mark item 575 BLOCKED in the same commit, and report."

### Review note 2026-10-21 03:00 (item 575; entry 45 (ii))

`xl` entry 40 ran green at 02:00 (`20261021T070008Z_GEO-38-step1b-C1b.log`: t_b 317.0 s, t_col 199.0 s at U₁ 1 525 193, 40.37 GiB at `-n 8`), so **entry 45's (ii) is MET**; where 575's brief says entry 45 "stays PENDING on (ii)", read "on (i) only". Nothing else in item 575 changes: the daytime rule and the loadavg ≤ 4.0 gate stand. The implementer landing 575 green marks entry 45's (i) MET in the same commit and does not write the command; the 10-22 03:00 review writes it from entry 40's measured column.

### §9 item 593 — 2026-10-22 (§10 I36a, chain step (no solve); lane 1)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Outcome 2026-10-22 (17:00 UTC slot, slot-first): DONE.** Module `tests/validation/test_geo36_rod_f_clinical_drive.py` (new; no other file edited). W0 (`-n 1`, complex): 25 passed / 4 skipped (`20261022T170241Z_GEO-36-I36a-W0.log`); (N-U) 5.347e-02, (N-route) implant_rod None. Pre-checks: `pgrep -c python3` 0, loadavg 0.61 0.28 0.27. W1 (`-n 2`, `build`): 5 passed / 5 skipped, rc 0, pytest 479.8 s, harness 482 s (`20261022T170259Z_GEO-36-I36a-W1-build-n2.log`). t_build 477.7 s (I35's rod build 470.9 s, gmsh 438.49 s; PREDICTED ≈ 464 [410, 560]); (A0) 1 217 025 rel 0; U₁ 1 606 744, (Deg-2) 8 327 394 = 2E + 2F exact (F 2 556 953); (D-rod) (G1) / (G2) 4.441e-16, (G4) 6.457e-11; (M) / (N-M) as I35; (C) 64, (W) 2.220e-15, (Az) 5.684e-14°; (F) P33 / P41 at 5.625° / 95.625°, bisector 50.625°; (Z) Z_RC 0.2 − 24.868j Ω, Z_th 10.012 − 19.809j Ω; **(M2)** degrees [1, 2], tag 403 σ 800 / ε_r 78 / μ_r 1 = gel's ε_r / μ_r, tag 3 σ 0.5 / ε_r 78, exact; rod azimuth 0.000° (−50.625° from the bisector); summed RSS 2.168 GiB; window 478.2 s. The `:71` split recurred benignly (the gated `PASSED` at `:11996`, after every I36a line). I36's prerequisite (i) is MET.

**Licence.** §10 chain I, the 2026-10-21 weekly (`PROJECT_PLAN.md:4922–4927`): "I36a (lane 1 → standard, build-only, `-n 2`; serial on nothing). The rod fixture's two-port drive module in C4a's pattern (the two-port feed on the bisector at P33 / P41, the rod in, a `FEM_EM_GEO36_S8_DRIVE` knob = `build` / `1`, the degree knob), build-only green at `-n 2` asserting I35's 1 217 025 cells and U₁ 1 606 744; I36's prerequisite (i). No solve, no tip ratio."
- I36 itself is "the two-port drive with / without the rod, the tip 1 g ratio at clinical scale". It is **unenumerated until entry 45 prices the family** (the 10-28 weekly reads it; `PROJECT_PLAN.md:4941–4942`; `docs/testing/xl-pending.md:2071`, `:2077`). I36a lands the module and its build proof only.
- No band, tolerance, record, Done-when or target is set or moved. Every band is imported. The two records asserted are I33's measured integers, already registered in I35's module under the (1\*) licence.

**Executor and cost.**
- **Executor:** implementer.
- **Tier: heavy by measurement, not the weekly's "standard".** The build-only proof rebuilds I33's rod-in-box mesh in memory, and that build has never fitted 180 s:
  - I35 W1b: t_build 457.5 s (gmsh 426.92 s), window 465.6 s, harness 468 s, at a pre-check loadavg of 0.36 (`20261021T170059Z_GEO-36-I35-W1b-build-n2.log:11986`, `:11992`, `:12004`);
  - I33 W2: 444.9 s (gmsh 412.45 s) (`20261019T201444Z_GEO-36-rodbox-W2-n2.log`).
  - Container `timeout -k 30 590`, the most a foreground window carries inside the Bash tool's 660 s. Predicted headroom ≈ 120 s.
- **Ranks:** W0 `-n 1`, W1 `-n 2`, both on the complex image (the fixture builds a `TimeHarmonicSolver` and compiles the sheet-width `dS` form).
- **Daytime and slot-first** (`PROJECT_PLAN.md:3964`): 12:00 / 13:30 / 15:00 only, as the slot's first item, one daytime item per slot. Any other slot skips it at no attempt cost. **Never in the same slot as item 595 (C12) or any other clinical build: two heavy windows never run at once** (`:3969`, `:4085`).
- **Load gate (item 560's lesson).** Item 560's W1 on this same mesh was killed at 590 s at a pre-check loadavg of 11.96: every gmsh phase ran 1.1–2× slower than I33's (GEO-36.md, 2026-10-20 18:00 ruling 1). Item 575 then ran it green at 0.36 under a ≤ 4.0 gate. I36a carries that gate unchanged.
- **Class:** chain step (no solve). **Cost:** ≈ 8.0 measured (class median) / **≈ 42 predicted** slot-minutes: 15 fixed, code ≈ 12 (one new module, ≈ 200 lines, composition only), W0 ≈ 1, pre-checks ≈ 1, W1 ≈ 8 (I35 W1b's 468 s harness), commit and journal ≈ 5.
  - **If W1 has not started by minute 42:** commit the module with W0 green and mark item 593 **OPEN for "W1 only"**. The next daytime slot runs the pre-checks and W1 on the committed file, unchanged (rule (i)). That is not an attempt.
- **Gmsh stop rule: not armed.** This is a priced rebuild (I33, I35). `FEM_EM_GMSH_STOP_RULE` is never set.
- **Depends on:** nothing (I35 landed, `7d065e0f`; C4a landed, `11424a82`). **Family:** `GEO-36` step 8 two-port drive module (I36a), attempt 1 of 4.

**What exists (`main`, 2026-10-22).**
- **I35's module**, `tests/validation/test_geo36_rod_f_clinical.py` (693 lines, `7d065e0f`). Never edit it:
  - `_rod_centre()` / `_kwargs()` (`:88–97`): `_box_kwargs()` plus the kept rod at (Lx/4 = 0.105, `TABLE_OFFSET_M` = 0, 0), so the rod sits at azimuth 0° in the coil frame;
  - `F_CLINICAL_ROD_CELL_RECORD = 1_217_025` (`:410`) and `F_CLINICAL_ROD_UNKNOWNS_RECORD = (1_606_744, 8_327_394)` (`:412`). These are I33's E 1 606 744 and F 2 556 953;
  - `_build_rod_f_clinical(comm)` (`:449–548`). It builds the mesh and runs `col._context_from_mesh`, then remaps tag 403 to σ = `SIGMA_WIRE_S_PER_M` with the gel's ε_r / μ_r (one solver, `del ctx["solver"]` first). It asserts (D-rod) (G1) / (G2) / (G4) / (G7), (M), (N-M), (C) and (U1). It returns `col._build_f_clinical`'s dict shape plus `problem_gel_only` and `t_build`;
  - `_route_column` (`:551–559`) sets `FEM_EM_GEO38_STEP1B=1`. **It is not reused here** (see traps).
- **C4a's module**, `tests/validation/test_geo38_f_clinical_esr_root.py` (875 lines, `11424a82`). It is `xl-pending.md` entry 46's module (10-26), so never edit it. Each of its fixture entry points reads module globals at call time:
  - `_build_c4_fixture(comm)` (`:270–364`). It calls `col._build_f_clinical(comm)` as a module attribute (`:289`). It imports `F_CLINICAL_R1_CELL_RECORD`, `F_CLINICAL_R1_UNKNOWNS_RECORD` and `_f_clinical_kwargs` from `test_geo38_f_clinical_build` at call time (`:276–283`). It asserts (f) against `FREQUENCY_HZ_C4` = 64 MHz (`:296`, `:313`), then (A0), (iv), outer census, (S) (`:314–318`), then (C), (W), (Az), (U1) (`:342–346`) and (F), which picks P33 / P41 with Δaz 90° (`:349–358`). It then calls `_apply_degree_c4` (`:360`).
  - `_run_build(comm, monkeypatch)` (`:370–428`) adds (Z) at the 100 pF wiring value (`:388–400`). It then runs (Deg-2) in the same build: degree-2 unknowns = 2E + 2F = the record's `[1]`, and the NEG identity (`:403–428`). It returns before any wall term, assembly or factorisation.
  - `_apply_degree_c4` (`:236–267`) builds the degree-2 solver on `ctx["solver"].problem`.
- **C4a's measured build** (1.5 T stand-in fixture): t_build 335.5 s, t_d2 0.3 s, window 336.6 s, RSS 2.12 GiB (`20261020T183431Z_GEO-38-C4a-W1-build.log:11859`, `:11870`). (F) P33 / P41 at 5.625° / 95.625°, bisector 50.625° (`:11864`).

**Code: one new file, `tests/validation/test_geo36_rod_f_clinical_drive.py`; edit nothing else.** No `src/`, no I35 module, no C4a / C3a / column / build module, no `xl-pending.md` (I36 has no entry yet). **C4a's build route is run unchanged on the rod fixture by composition**, which is I35's and C11's no-edit route.
- **Knobs:**
  - `S8_DRIVE_ENV = "FEM_EM_GEO36_S8_DRIVE"`: unset ⇒ the gated test skips; `build` ⇒ build-only; `1` ⇒ reserved for I36 (below); anything else ⇒ `ValueError`.
  - `S8_DRIVE_DEGREE_ENV = "FEM_EM_GEO36_S8_DRIVE_DEGREE"`: unset ⇒ 1; `2` ⇒ 2; else `ValueError` (C4a's `_degree_knob_c4` semantics under its own name). Under `build` a set degree knob is refused, because the build proof covers both degrees by construction ((Deg-2)).
  - **Foreign-variable guard:** the gated test raises `ValueError` on entry if any of these is set: `FEM_EM_GEO36_S8_COLUMN`, `FEM_EM_GEO36_ROD_BOX`, `FEM_EM_GEO38_STEP1B`, `FEM_EM_GEO38_C4`, `FEM_EM_GEO38_C4_NODES_PF`, `FEM_EM_GEO38_C4_DEGREE`, `FEM_EM_GMSH_STOP_RULE`. This way an I36a window never reads an I35 or C4 variable.
- **Constants:** none retyped. Records are read as `i35.F_CLINICAL_ROD_CELL_RECORD` / `i35.F_CLINICAL_ROD_UNKNOWNS_RECORD`, and 1.5 T records as `s1a.…` (for (N-U) only).
- **`_route_drive(monkeypatch)`**, four patches, each restored at teardown:
  1. `col._build_f_clinical` → `i35._build_rod_f_clinical`;
  2. `s1a.F_CLINICAL_R1_CELL_RECORD` → `i35.F_CLINICAL_ROD_CELL_RECORD`;
  3. `s1a.F_CLINICAL_R1_UNKNOWNS_RECORD` → `i35.F_CLINICAL_ROD_UNKNOWNS_RECORD`;
  4. `c4._apply_degree_c4` → a recording wrapper. It calls the original and then stores `(degree, solver.problem.material_map)` for (M2).

  It sets **no** environment variable.
- **`test_geo36_i36_rod_drive(monkeypatch)`**, `@complex_only`:
  - guard, then mode `None` ⇒ skip;
  - `build`: `_route_drive`, then `c4._run_build(MPI.COMM_WORLD, monkeypatch)`. After it returns, assert (M2). Print the rod's azimuth atan2(y, x) of `i35._rod_centre()` (predicted 0°) beside the bisector C4a prints in its (F) line (predicted 50.625°). Then print `BUILD-ONLY return` (window, RSS via `_rss_gib`);
  - `1`: `raise NotImplementedError("I36's drive body lands with I36; the 10-28 weekly enumerates it after entry 45 reads")`. This runs **before** any build. See for-the-review item 1.
- **Coding tests**, skipped whenever `S8_DRIVE_ENV` is set (C4a's `_coding` marker pattern): (K0), (K1), (K2), (N-route), (N-U) below.
- **Imports:** module objects only (`i35`, `c4`, `col`, `s1a`), inside function bodies. Never `from … import test_*`: pytest would collect I35's two gated builds or C4a's gated test into this module.

**Anchors (asserted in W1 at `-n 2`; owned entities, reduced; collectives before any rank-0 branch).** All but (M2) are asserted by the composed code itself, against the rod records through the patch.
- **From `_build_rod_f_clinical` (I35's asserts, unchanged):**
  - (D-rod): (G1) / (G2) / (G4) ≤ `STRIP_BAND` 1e-9 (`tests/validation/test_ops61_exact_strip_ports.py:58`), and (G7) census clean;
  - (M): tag 403 σ = 800 S/m and ε_r = the gel's, tag 3 σ = `SALINE_SIGMA` exact;
  - (C) 64 specs; (U1) = 1 606 744 exact.
- **From `c4._build_c4_fixture`:**
  - (f) ω/2π against 64 MHz ≤ 1e-12;
  - **(A0) cells = 1 217 025 within `CELL_COUNT_BAND` 0.01** (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`);
  - (iv), outer census, (S) cavity < R_s 0.3725 m;
  - (C) census 64, (W) ≤ `STRIP_BAND`, (Az) < `AZIMUTH_MATCH_DEG` 1e-6 (`tests/validation/test_port_birdcage_ring_column.py:158`);
  - **(U1) = 1 606 744 exact**;
  - (F) feed2 unique on P33's ring at Δaz 90° within `AZIMUTH_MATCH_DEG`.
- **From `c4._run_build`:**
  - (Z) at the 100 pF wiring value: Re Z_RC = `ESR_OHM` 0.2 Ω exact, Im at 1e-12, Z_th at 1e-15, κ `None` on 64;
  - **(Deg-2) = 2E + 2F = 8 327 394 exact**, with the NEG identity n(2)·E = n(1)·(2E + 2F).
- **(M2) new, I36a's own, map-level, exact.** On both recorded `_apply_degree_c4` calls (degrees [1, 2], asserted), `material_map[403]` has σ = `SIGMA_WIRE_S_PER_M` (`tests/validation/test_port_gap_voltage_impedance.py:235`) and ε_r / μ_r equal to `material_map[3]`'s, and `material_map[3]` has σ = `SALINE_SIGMA` and ε_r = `SALINE_EPSILON_R` (`tests/validation/test_lossy_sphere_fullwave.py:92–93`). This proves the rod reached the **degree-2** problem, which I35 never built. Equality of imported constants, not a band.
- **(K0)** Knob parsing: S8_DRIVE unset / `build` / `1` accepted, `2` / `x` refused; degree unset / `2` accepted, `1` / `x` refused; degree set under `build` refused; each foreign variable refused by the guard, through `monkeypatch`.
- **(K1) Routing tripwire (no gmsh).** Inside `monkeypatch.context()`, apply `_route_drive`, and replace `MeshGenerator.birdcage_port_domain` with a recorder that raises a sentinel. Then `c4._build_c4_fixture(MPI.COMM_WORLD)` must raise the sentinel with kwargs that carry `implant_rod` with `keep_cells` True at `i35._rod_centre()`, plus `phantom_box`. On exit: `col._build_f_clinical` and `c4._apply_degree_c4` are the originals, `s1a.F_CLINICAL_R1_CELL_RECORD == 1_147_366`, `s1a.F_CLINICAL_R1_UNKNOWNS_RECORD == (1_525_193, 7_885_674)`, and `FEM_EM_GEO38_STEP1B` is unset.
- **(K2) The `1` stub stops before gmsh.** With S8_DRIVE = `1` and the raising recorder in place, the gated test raises `NotImplementedError`, not the sentinel.

**Negative controls (rule (e)).**
- **(N-route) Asserted, binary (I35's (N-route) on C4a's entry point).** The same tripwire call without `_route_drive` reaches kwargs with **no** `implant_rod`: the 1.5 T stand-in fixture. Backed by I35's (N-route) on the column entry point (`20261021T170043Z_GEO-36-I35-W0b.log`).
- **(N-U) Asserted by construction.** The 1.5 T unknown records miss the rod mesh: |1 606 744 / 1 525 193 − 1| = **5.347e-2**, and 8 327 394 ≠ 7 885 674. (U1) and (Deg-2) are exact equalities, so a record patch that silently failed turns W1 red, not green. The cell-count twin is I35's (N-A0), 6.071e-02 = 6.07× `CELL_COUNT_BAND` (same W0b log). It re-runs in this W0 because I35's module is on the path.
- **(N-M) Asserted, inherited** (`_build_rod_f_clinical:543`). The unpatched problem reads σ = 0.0 on tag 403, against 800 (I35 W1b: σ 0.0).
- C4a's own (N1) / (N2) / planted (Z) controls run in W0 with C4a's module on the path.

**Predictions (written first; printed, never asserted).**
- (A0) rel 0;
- t_build (C4a's timer, which includes I35's context and rod checks) ≈ 464 s [410, 560] at 1-min loadavg ≤ 4. That is I35's 457.5 s + 6.5 s (`…I35-W1b…log:11990`);
- C4a's fixture checks ≈ 1 s; t_d2 ≈ 0.3 s (C4a);
- window ≈ 467 s, pytest ≈ 470 s, harness ≈ 473 s;
- summed RSS ≈ 2.1 GiB after the degree-2 space (I35 2.089; C4a 2.12);
- F = 2 556 953 printed beside (U₂ − 2U₁)/2;
- (F) P33 / P41, 5.625° / 95.625°, bisector 50.625°;
- rod azimuth 0.000°, i.e. −50.625° from the bisector and −5.625° from P33;
- (Z) Z_RC = 0.2 − 24.868j Ω, Z_th = 10.012 − 19.809j Ω, as C4a (`…C4a-W1-build.log:11866`);
- the `:71` chunked-forwarding pattern (2026-10-20 18:00 ruling 2) may recur. It is benign unless rank 1 reports `PASSED` before rank 0's `[GEO-38 C4] BUILD-ONLY return` line.

**Windows.** Each runs in the foreground with Bash `timeout: 660000`, pytest `-s`, never `-k`. **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**
- **W0** (complex, `-n 1`, ≈ 30 s). Expected: `test_complex_mode` green; this module's (K0), (K1), (K2), (N-route) and (N-U) pass and its gated test skips; I35's module reads as at W0b (refusal and coding tests pass, two gated tests skip); C4a's 8 coding tests pass and its gated test skips. That proves both composed modules are collectable and untouched.
  `scripts/testing/run_and_log.sh GEO-36-I36a-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && unset FEM_EM_GEO36_S8_DRIVE FEM_EM_GEO36_S8_DRIVE_DEGREE FEM_EM_GEO36_S8_COLUMN FEM_EM_GEO36_ROD_BOX FEM_EM_GEO38_STEP1B FEM_EM_GEO38_C4 FEM_EM_GEO38_C4_NODES_PF FEM_EM_GEO38_C4_DEGREE FEM_EM_GMSH_STOP_RULE && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_drive.py tests/validation/test_geo36_rod_f_clinical.py tests/validation/test_geo38_f_clinical_esr_root.py -v -s --tb=short'"`
- **Pre-checks, in order, each its own bare call:**
  1. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read **0**;
  2. `docker compose -f docker/docker-compose.yml exec -T fem-em-solver cat /proc/loadavg`. **Load gate: 1-min ≤ 4.0.** If it fails, do not run W1: commit the module with W0 green and mark item 593 **OPEN for "W1 only"** with the reading. That is not an attempt, and the slot takes no other daytime item.
- **W1** (build-only, complex, `-n 2`, the window of record; §5.1 durable capture, verbatim):
  `scripts/testing/run_and_log.sh GEO-36-I36a-W1-build-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_GEO36_S8_DRIVE=build timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_drive.py -s -v --tb=short > /workspace/logs/geo36-i36a-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo36-i36a-w1-raw.log; cat /workspace/logs/geo36-i36a-w1-raw.log; exit \$rc'"`
  Expected: the gated test passes and the coding tests skip.
- **Rule (i):** no module edit after W0 starts. Any edit means W0, the pre-checks and W1 again.

**Traps already paid for.**
- **Never reuse I35's `_route_column`:** it sets `FEM_EM_GEO38_STEP1B=1`. **Never put I35's module, the column module or C4a's module on W1's path.** Their gated tests would add a second ≈ 330–460 s build to a 590 s window, a certain kill (item 518's trap; item 560's W1 note).
- **Never edit C4a's module** (entry 46, 10-26), I35's module or the column / matrix modules. The composition depends on C4a's call-time names (`col._build_f_clinical`, the call-time `s1a` imports, `_apply_degree_c4` as a module global). (K1) is the tripwire that catches a renamed binding in W0, before the 8-minute window.
- **Load:** this build is load-sensitive (item 560: killed at 590 s at loadavg 11.96; item 575: green at 0.36). Run the gate, never skip it. Never run beside another heavy window.
- **I35-labelled lines in an I36a log.** `_build_rod_f_clinical` prints `[GEO-36 I35] …` lines. In this window they are the rod fixture's own asserts, not an I35 re-run. C4a's `_run_build` prints C3's npz path in its "never opened" line; that is a label, and no npz is read.
- **Tag 403 unmapped is a silent air bubble** (`_context_from_mesh` maps tag 3 only). (M) and (M2) catch it. **The "without the rod" twin I36 will need is tag 403 → the gel's material**, never `problem_gel_only` (air, σ 0, ε_r 1: the (N-M) control). It is named here for I36, not built.
- **The feed-bisector mirror is not available on this fixture.** The ASTM box's mirror planes are the coil-frame xz and yz planes, and the bisector at 50.625° is neither. The rod at azimuth 0° breaks the yz mirror as well. WF-10's mirror gate cannot be imported into I36 on the box. Named for the weekly.
- **One solver:** `_apply_degree_c4` deletes the degree-1 solver before building the degree-2 one (C4a). Never hold a second reference across it, including in the (M2) recorder: store the material map, not the solver.
- **Rank safety:** owned entities and reduced extrema (the composed code already does this). Never use UFL `dS` over the sheets for counts.
- **The sheet-width `dS` form compiles in W1.** If W1 is killed mid-compile, run `pgrep -c python3` first, then the rubric-4 sweep as `GEO-36-I36a-cache-sweep`. Delete orphan stems only.
- `pgrep -f "python3 -m pytest"` is denied by the bash guard; use `pgrep -c python3`. Durable capture ends in `; exit \$rc`; no pipe; never delete a log; never background a harness window. A killed window is not a result.

**Scope.** The module, its coding tests and a build-only proof of C4a's two-port build route on I33's rod-in-box mesh at degree 1 and the degree-2 count.
- **Not claimed:** any solve, drive, S entry, Y_in, root or tuning statement; a tip 1 g ratio, SAR or heating; the no-rod twin; I36's gates, command, tier or price; entry 45's reading; a PEC or insulated rod; an ASTM placement.
- `GEO-36` stays ✅, and `GEO-38` and `IMP-3` are unchanged.

**Negative result.** Report and stop. Park on `attempt/GEO-36-I36a-<UTC>` and mark item 593 BLOCKED in the same commit (rule (d)).
- **(A0), (U1) or (Deg-2) red:** the rod mesh moved on this image, or a record patch did not bind. Print both values. If (K1) was green, open a known-issues entry, because entry 45 reads the same records.
- **(D-rod), (M) or (M2) red:** the rod's CAD or material did not reach the problem. BLOCKED; never re-map in-slot.
- **(f), (C), (W), (Az) or (F) red:** the port layout or frequency is not C4a's identity on this mesh. BLOCKED with the reading. Never pick a feed by ordinal in-slot.
- **W1 killed at 590 s with the load gate met:** the load reading is refuted. Run `pgrep -c python3`, journal the loadavg, each gmsh phase's time beside I35 W1b's (152.8 / 22.3 / 59.2 / 7.3 + 178.8 s) and the last `[birdcage-mesh]` line, then park. Never retry in-slot.
- **A printed reading outside its prediction is a reading, never a red.**

No `PREDICTED RED:` line.

**Status it can move.**
- No glyph. The `GEO-36` row is ✅ and its Done-when (step 2's identities) is unaffected, so rule (j) writes no "→ ✅".
- §10 I36a gains a DONE line. **I36's prerequisite (i) is MET** for the 10-28 weekly to read when it poses I36 at entry 45's price. There is no `xl-pending.md` entry to edit, since I36 has none.
- If the `GEO-36` row stays ≤ 2 KB it may gain "I36a (item 593): the rod fixture's two-port drive module, build-only green at `-n 2` (1 217 025 cells, U₁ 1 606 744, U₂ 8 327 394)".

**For the review.**
1. **The `1` branch.** The weekly gives the knob a `1` value, but I36 is unenumerated and its target is the 10-28 weekly's. The options:
   - (a) The drafter's default: `1` is parsed, and the gated test stops with `NotImplementedError` before any build ((K2)). I36 adds its body, insertions only, and its own window re-asserts every build anchor before solving.
   - (b) Write `1` now as C4a's single-node pair (rod in, then the gel twin, at a node from a required `…_C_PF` knob) with (P_RC) asserted and the tip ratio printed. This writes I36's gate set ahead of the weekly, and C4a's node route cannot be composed as-is, because it asserts C3's 1.5 T npz cells (`:485`), which the rod mesh cannot match.
2. **Tier:** heavy by measurement (457.5 s build), container 590 s, with the ≤ 4.0 load gate, against the weekly's "standard".
3. **(M2)** as an asserted exact identity (the rod on the degree-2 problem), against printed only.
4. **The mirror and the no-rod twin** (traps above): the weekly should know both before it writes I36's gates.
5. **Daytime order with item 595 (C12).** Neither has a deadline before the 10-28 weekly; each takes its own daytime slot. The drafter suggests I36a first, since its build is the longer one and the more load-sensitive.

**Every-item rules** (PROJECT_PLAN §9 "Every item"):
- bare whole harness calls;
- a lost docker socket stops the slot;
- no edit after the last window;
- `timeout: 660000` and `-s`;
- never delete a log;
- BLOCKED in the parking commit;
- a killed window is not a result;
- two heavy windows never at once;
- no slot opens `docs/private/` or `aed_results/`.

**Execution order.** Read the brief and any ruling beneath it → code → W0 → pre-checks (pgrep, loadavg gate) → W1 → commit (the module, both harness logs with their `test-results.md` rows, the §9 DONE line, the §10 I36a DONE line, the optional row clause) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 593 — `GEO-36` I36a, the rod fixture's two-port drive module, build-only — per `docs/planning/chunks/GEO-36.md` §9 item 593. Read the brief and any review ruling beneath it first; the rulings hold. Run it only if this session's start line reads 12:00, 13:30 or 15:00 (≤ 5 min late) and it is the slot's first item; otherwise skip it. Never run it in the same slot as item 595 or beside any other heavy window. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s; container `timeout -k 30`. Create only `tests/validation/test_geo36_rod_f_clinical_drive.py`; edit no other test module and no `src/` file; never edit `tests/validation/test_geo38_f_clinical_esr_root.py`, `tests/validation/test_geo36_rod_f_clinical.py`, `tests/validation/test_geo38_f_clinical_column.py` or `tests/validation/test_geo38_f_clinical_matrix.py`; import module objects only, never `test_*` names; never reuse `_route_column`. Never set `FEM_EM_GEO36_S8_COLUMN`, `FEM_EM_GEO36_ROD_BOX`, `FEM_EM_GEO38_STEP1B`, any `FEM_EM_GEO38_C4*` or `FEM_EM_GMSH_STOP_RULE` on a command line. Run W0; then `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` (must read 0) and `docker compose -f docker/docker-compose.yml exec -T fem-em-solver cat /proc/loadavg` as separate bare calls; run W1 only if the 1-minute loadavg is ≤ 4.0, otherwise commit the module with W0 green and mark the item OPEN for 'W1 only'. If W1 has not started by minute 42, do the same. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a kill at 590 s with the gate met, or any red, park on `attempt/GEO-36-I36a-*`, mark item 593 BLOCKED in the same commit, and report."

### Review ruling 2026-10-22 03:00 (item 593, I36a)

Binds over the brief. **Heavy, daytime, slot-first** (I35's measured 457.5 s build), not the weekly's "standard"; the `pgrep` pre-check and the drafter's loadavg ≤ 4.0 gate are pre-checks, not bands. (1) The `=1` branch stops before any build (I36 is the 10-28 weekly's). (2) The rod-material-on-degree-2 check is asserted as an exact map identity. (3) The two points for I36's gates (no feed-bisector mirror in the box; the without-rod twin maps tag 403 to gel) are recorded for the weekly. (4) Runs before C12 (595), never in the same slot.

### Review note 2026-10-25 03:00 — `xl` entry 45 (I34) RUN, row amended

`20261025T070008Z_GEO-36-step8-I34.log`: 9 passed and 5 skipped on all 8 ranks, Status 0, 773 s, `memory.peak` 41.76 GiB. The skips are by design: the build test skips with `FEM_EM_GEO36_ROD_BOX` unset, and the I35 coding tests skip when the column knob is set.

Asserted green on the column path:
- (A0) 1 217 025 cells, rel 0;
- (D-rod) G1 / G2 4.441e-16, G4 6.457e-11, G7 clean;
- (M), (N-M), census 64, U₁ 1 606 744 exact (`:12093–12098`);
- (P) 4.021e-13 / 1.367e-13, P_surf 1.115e-06 W > 0, NEG 1.6237e-04, (R1) 1.858e-14 (`:12318`).

Price: t_b 458.2 s, t_col 277.8 s, t_held 2.06 s, exponent 1.643 (`:12316`). Every axis is inside entry 45's brackets.

**Pathologist caveat:** rod-in-box containment ((G3)) is **not asserted** on this path; it lives in the skipped build test. Entry 45's Readout lists it as asserted, so that wording is returned to the 10-26 weekly. Read literally, the decision rule's first branch fires: I36 (entry 48) is priced from this window, and the 10-26 03:00 review fills entry 48's price as written there. The ledger row carries the full reading. Glyph unchanged (✅).

### §9 item 656 — 2026-10-27 (§10 I36b, chain I — `GEO-36` step 9's body; chain step (no solve))

*Written by the 2026-10-27 03:00 daily review from the 10-26 weekly's I36b (`PROJECT_PLAN.md` §10 chain I). The weekly's design binds: a required node knob, arm M with the rod, arm R with tag 403 mapped to gel on the same mesh, the gates at their imported bands, no mirror gate. Three facts the weekly's text did not have (fact sheet, 2026-10-27) are written into the route below.*

**Daytime, slot-first** (as item 593): only in the 12:00, 13:30 or 15:00 slot, as its first item. Before W1, `uptime` loadavg (1 min) must read < 8; otherwise skip the item this slot (no cost).

**What exists.** `tests/validation/test_geo36_rod_f_clinical_drive.py` (352 lines): knobs `FEM_EM_GEO36_S8_DRIVE` (`:48`) and `…_DEGREE` (`:49`); `_route_drive` (`:119–140`) patches C4's builder to `i35._build_rod_f_clinical` and the records to the rod's (`F_CLINICAL_ROD_CELL_RECORD` 1 217 025, `…_UNKNOWNS_RECORD` (1 606 744, 8 327 394), `test_geo36_rod_f_clinical.py:410–412`); the gated test `:176` raises `NotImplementedError` under `1` at `:181–182`; the `build` branch (`:184–213`) runs `c4._run_build` through the route and asserts (M2).

**Route — the body, insertions only except (K2).**
1. **Node knob** `FEM_EM_GEO36_S8_DRIVE_C_PF`: **required** under `1` (missing ⇒ `ValueError` before any gmsh), refused under `build`, parsed as a float with ≤ 4 decimals.
2. **(K2) changes meaning, not strength.** `test_i36a_k2_stub_stops_before_gmsh` (`:311–320`) asserts the stub's `NotImplementedError`. The body removes the stub, so (K2) becomes: "`=1` **without** the node knob stops before gmsh" (the same gmsh recorder, now expecting `ValueError`). It is still a stop-before-gmsh test. Say so in the commit; this is the only edit to an existing line.
3. **Arm M (the rod in):** the mesh from `i35._build_rod_f_clinical(comm)` with tag 403 at the wire σ (the (M2) map), r = 0.2 Ω ESR on all 64 gaps at the node C (`h17.ESR_OHM`, `_esr_specs`, `test_geo38_f_clinical_esr_root.py:192`), feeds P33 / P41 (`_feeds_32`, `:155`), the (1, +1j) drive as one factorisation with two held columns. **`c4._run_nodes` cannot be reused:** it requires C3's npz with 1 147 366 cells (`:467`, `:485`), reads the C4 node knobs (which are in this module's `FOREIGN_ENVS`), writes C4's JSONs and computes no (R1). Its `_drive` / `_pair` / `_partition` are closures (`:497`, `:551`, `:568`), so **copy their logic into module-level helpers here, citing the source lines**, and import the module-level pieces that exist: `_feeds_32`, `_slots_32`, `_esr_specs`, `_z_failures`; `h17._thevenin_feed_impedance`, `_source_voltage`, `_z_rc`, `_active_y_in`, `SENSES`.
4. **Arm R (gel twin, same mesh, build-free second factorisation):** `replace(problem, material_map={**mmap, 403: gel})`, I41's route as written in `tests/validation/test_imp3_rod_heating_matched.py:104–118` (`_arm_with_403`). **Not** `_build_kept_rod` (that builds F-small) and **not** the generator's `"problem_gel_only"` (tag 403 unmapped, σ = 0, `:478`, `:543`). (A0m): arm R's tag-3 + tag-403 census equals arm M's, exactly.
5. **Gates (imported, unmoved):** (A0) within `CELL_COUNT_BAND` (`test_birdcage_port_sheet_prerequisite.py:57`), census 64, U₁ exact; (D-rod) as the generator asserts it (`:487–538`); (P) and NEG (`TH-14`'s identity, C4's `_partition` logic) on both arms; (R1) on the 2×2 at `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`; the precedent is `test_geo38_f_clinical_column.py:481–519`, helper `_reciprocity_ratio`); (P_RC) at the node at `DISCRETE_IDENTITY_RTOL` (`test_birdcage_power_identity.py:227`).
6. **Printed only:** the tip 1 g ball ratio Q and η_rod(T1) at degree 1, with the tips from `i35._rod_centre()` (`:88`) ± `ROD_LENGTH` / 2 (**not** `TIPS`, which are F-small coordinates, `test_imp3_rod_birdcage.py:140`), computed as `_arm_d` (`test_imp3_tip_ratio_1g.py:150`) and `_w4` (`:401–404`) do, beside `ENTRY43_READOUT_INTERVAL` (3.48, 3.56) (`test_imp3_tip_ladder_1g.py:176`); the rod's share of the dissipated power; Im Y_in at the node (the coil was tuned on C4's cylinder, so the box's detuning is printed, never a red).

**This item runs `build` only.** The `1` branch first runs on `xl` as entry 48.
- **W0** `GEO-36-I36b-W0` (`-n 1`, complex, smoke): every coding test, including the new (K2) and the node-knob parsing.
- **W1** `GEO-36-I36b-W1-build` (`-n 2`, complex, `timeout -k 30 600`, `FEM_EM_GEO36_S8_DRIVE=build`): the existing build branch, unchanged, still green with the body present (it proves the body imports and the route still patches). Predicted ≈ 470–520 s (I35 457 s; item 593's W1). Over 600 s is a kill: park, never raise the timeout.
- **W-env** (`-n 1`, `timeout -k 30 300`).

**Executor:** implementer. **Tier:** heavy (build), `-n 2`, complex. **Class:** chain step (no solve). **Cost:** ≈ 9.2 measured / 50 predicted slot-minutes. **Depends on:** nothing. **Family:** I36b (`GEO-36` step 9), attempt 1 of 4.

**Traps already paid for.** gmsh at this size is load-sensitive (the 10-20 kill was load). A body path no window executes is unverified code: keep the body a composition of cited, already-gated pieces, and make the coding tests exercise every pure-python helper on small arrays (e.g. the (R1) and (P_RC) arithmetic on a synthetic 2×2). `assemble_scalar` is rank-local. UFL comparisons containing a `sqrt` fail in the complex build. Count verdicts per rank by token.

**Scope.** The body and its build. **Not claimed:** any solve, tip ratio, power share or admittance (entry 48's), any heating.

**Negative result.** Build red ⇒ the body broke the route: park, BLOCKED. Killed ⇒ BLOCKED with the last line. A coding-test red ⇒ fix the helper, never the assertion. No `PREDICTED RED:` line.

**Status it can move.** §10 I36b → DONE; **`xl` entry 48 → READY**: the implementer marks it READY in the same commit (step 6b rule 5), and the 10-29 03:00 review writes its command with the node (entry 51's refinement, or 193.5178 pF) and queues it for 10-30; item 667 (C6a) is unblocked.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 656 — `GEO-36` I36b — per `docs/planning/chunks/GEO-36.md` "§9 item 656". Daytime, slot-first. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Build-only: never run the `1` branch. The only edit to an existing line is (K2), as the brief states. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome 2026-10-27 (12:00 slot): DONE.** Body landed (`_run_drive`, helpers `_node_knob`, `_partition`, `_p_rc_split`, `_gel_twin_map`, `_a0m_failures`, `_probe_points`, `_drive`, `_arm`; `_drive` / `_pair` / `_partition` lifted from C4a `:497–572` with citations). Loadavg 0.53 and `pgrep -c python3` 0 before W1. W0 `20261027T170539Z_GEO-36-I36b-W0.log` (`-n 1`): 30 passed / 4 skipped, 5 s — (K2) under its new meaning, K3 node knob (3 accepted, 7 refused, build + knob stops before gmsh), K4 partition / (P_RC) arithmetic, K5 (R1) 3.223e-18 on a reciprocal synthetic circuit vs planted 1.029e-02 > 1e-3, K6 probe points, K7 gel twin map / (A0m). W1 `20261027T170559Z_GEO-36-I36b-W1-build.log` (`-n 2`): 5 passed / 10 skipped on both ranks, Status 0, 502 s; t_build 498.3 s; (A0) 1 217 025 rel 0; census 64; U₁ 1 606 744 exact; (D-rod) G1 / G2 4.441e-16, G4 6.457e-11; (M2) none; U₂ 8 327 394. W-env `20261027T171433Z_GEO-36-I36b-Wenv.log`: 12 passed, 52 s. **(K2) changed meaning** (stub `NotImplementedError` → missing node knob `ValueError`, still stop-before-gmsh). Not exactly insertions-only besides (K2): the three stub lines in the gated test were replaced by the knob read and the `_run_drive` dispatch (the brief's "the body removes the stub"); the module docstring's "reserved" line is left and superseded by an inserted paragraph. Tips: the brief's ± `ROD_LENGTH`/2 (the cap centres) are printed beside ±(L/2 + r), which is where F-small's `TIPS` sit. `xl` entry 48 → READY. The `1` branch has not been executed.

### §9 item 695 — 2026-10-30 (§10 I36c, chain I — daily, 2026-10-30, provisional; the `=1` glue on a small mesh before entry 48's window; chain step (no solve))

*Written by the 2026-10-30 18:00 daily review. Depends on nothing. `xl` entry 48 runs Monday 11-02. The 10-26 weekly and the 10-28 03:00 review both flagged its risk: "the `=1` glue (`e.copy()`, superposition, the 1 g SAR, the second solver) has never run on any mesh; no cheap smoke exists (the module has no sizing knob, and its cell records are exact)". That sentence is true of the **module's fixture**. `_arm`'s post-processing takes its fields from `_drive`, so it can run on a small mesh with `_drive` stubbed. This item removes the part of the risk that is glue; the second factorisation stays entry 48's.*

**Route — insertions only, `tests/validation/test_geo36_rod_f_clinical_drive.py`:** one new `_coding` test, `test_i36c_k8_arm_glue_on_small_mesh(monkeypatch)`. No line of `_arm`, `_drive`, `_vol_powers` or `_run_drive` changes.
- **Mesh:** `dolfinx.mesh.create_box` on [−0.05, 0.05]³ m (12 × 12 × 12 cells, tetrahedra), cell tags 3 everywhere and 403 on cells whose centroid has max(|x|, |y|, |z|) < 0.012, passed through `ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)`.
- **Fields:** N1curl degree 1, E_u interpolated from the uniform (E₀, 0, 0), E₀ = 3.0 V/m. σ DG0 0.5 S/m. `rho_field` from `build_density_field` with ρ = `PHANTOM_RHO_KG_PER_M3` on every tag.
- **`_drive` stub** (`monkeypatch.setattr(i36, "_drive", stub)`): returns the dict `_arm` reads, `{"currents", "kind", "t_solve", "t_post", "p_*", "e", "sigma"}`, with `kind` "solve" then "held", `e` = `E_u.copy()` on both columns, and currents from the reciprocal synthetic circuit K5 already uses (`test_i36b_k5_r1_and_current_route_on_synthetic_2x2`).
- **`env`:** the dict `_arm` and `_vol_powers` read (`comm`, `dx`); the `solver` argument is `None`, which the stub ignores.
- **`tips`** two points in tag 3 at least one 1 g ball radius from tag 403 and from the walls (print the radius); **`t1`** = (0.03, 0.03, 0.0).

**Anchors (closed form, asserted).** `_arm`'s `w` = (1, +1j) on e₁ = e₂ = E_u gives e_sup = (1 + j) E_u, so |e_sup|² = 2 E₀²:
- (K8-P) `p_sup["total"]` = ½ σ · 2 E₀² · V_mesh ≤ `EXACT_IDENTITY_RTOL`; `p_sup["rod"]` / `p_sup["total"]` = V₄₀₃ / V_mesh ≤ `EXACT_IDENTITY_RTOL`, with V₄₀₃ reduced from the cell tags.
- (K8-S) each tip's 1 g SAR = σ E₀² / ρ ≤ `EXACT_IDENTITY_RTOL`, and `mass` = ρ × the ball's quadrature volume > 0 (uniform field, uniform ρ in the ball).
- (K8-E) `e_t1` = ((1 + j) E₀, 0, 0) ≤ `DISCRETE_IDENTITY_RTOL`, `t1_valid` True.
- (K8-k) `kinds` == ["solve", "held"]; `s2` reciprocal ≤ `RECIPROCITY_BAND` on K5's synthetic circuit.
- **Negative control (asserted):** the same stub with e₂ = 0 gives |e_sup|² = E₀², so (K8-S)'s closed form misses by 0.5 > `EXACT_IDENTITY_RTOL` (ceiling 1.0; the miss is exact, not statistical).

**Windows.** W0 `GEO-36-I36c-W0` (`-n 1`, complex, `timeout -k 30 120`) and W0b `GEO-36-I36c-W0b` (`-n 2`, same): the module's coding tests with `FEM_EM_GEO36_S8_DRIVE` unset. W-env (`-n 1`, `timeout -k 30 300`). No gmsh. **Class:** chain step (no solve), ≈ 8.8 measured / 20 predicted minutes. Tier smoke, measured against item 656's W0 (30 passed, 5 s).

**Traps already paid for.** `_arm` calls `gc.collect()` and `del` on its own fields, so the stub hands it copies. `mass_averaged_sar` in `_arm` has no `comm=` argument: if `-n 2` shows a rank-local value, that is a **finding about entry 48's window**, not a test to adjust. Every point read goes through `evaluate_vector_field_parallel`, already so in `_arm`. Complex build: no comparison on a `sqrt`.

**Scope.** `_arm`'s superposition, volume powers, 1 g balls and T1 read at `-n 1` and `-n 2`. **Not claimed:** `_drive`, the held solve, the second factorisation, any clinical reading. The risk note in entry 48 is reduced, not retired.

**Negative result.** A red here is a defect that entry 48 would have met after ≈ 750 s of window: **BLOCKED in the same commit**, parked on `attempt/GEO-36-I36c-*`, with the failing line. The slot does not touch the queue file; the 10-31 reviews hold `xl-queue.d/20261102-GEO-36-step9-I36.env` (daily-review step 6b) until a fix lands. A fix to `_arm` is **not** this item's. No `PREDICTED RED:` line.

**Status it can move.** §10 I36c → DONE; `xl` entry 48's "never run on any mesh" note narrows to the second solver; or entry 48 is held by a found defect.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 695 — `GEO-36` I36c — per `docs/planning/chunks/GEO-36.md` "§9 item 695". Any slot. Insertions only: no line of `_arm`, `_drive`, `_vol_powers` or `_run_drive` changes. Never run the `1` branch. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome 2026-10-30 19:30 slot — DONE, green; entry 48 stays queued.** W0 (`20261031T003838Z_GEO-36-I36c-W0.log`, `-n 1`, 17 s) and W0b (`20261031T003859Z_GEO-36-I36c-W0b.log`, `-n 2`, 6 s): 11 passed, 1 skipped, Status 0, which is the whole module (12 tests; the brief's "30 passed" for item 656's W0 is not this module's count alone). W-env (`20261031T003914Z_GEO-36-I36c-Wenv.log`): 12 passed, 48 s. Insertions only: `git diff` on the module is 132 added lines at its end, so no line of `_arm`, `_drive`, `_vol_powers` or `_run_drive` changed. No gmsh and no solve.

Readings, `-n 1` / `-n 2` (W0 `:75–78`, W0b `:114–117`):
- **(K8-P):** total 6.0e-13 / 4.8e-14, rod share 6.0e-13 / 5.1e-14, against 1e-10; V₄₀₃ = 1.041666667e-05 m³, V_mesh = 1.000000000e-03 m³.
- **(K8-S):** each tip's 1 g SAR ≤ 3.3e-16 / 4.7e-15; mass against ρ × the ball's quadrature volume ≤ 6.7e-15 / 7.1e-15; a_1g = 6.203505e-03 m.
- **(K8-E):** miss 1.15e-15 against 1e-6 (W0b `:116`), `t1_valid` True.
- **(K8-k):** kinds ["solve", "held"], (R1) 3.2e-18 against 1e-3.
- **Control, e₂ = 0:** misses 0.5 to 2.3e-15, exact as predicted.

`mass_averaged_sar` without `comm=` is not a finding: it defaults to the mesh's communicator (`post/sar.py:267–268`), and the `-n 2` readings agree with `-n 1`.

**Disclosed.** The test takes the mass comparand from a second `mass_averaged_sar` call at the same centre and radius, because `_arm` does not return the ball's volume. The tips are (0.03, −0.03, 0) and (−0.03, 0.03, 0) m. The stub is installed through `sys.modules[_arm.__module__]`, so it patches the module `_arm` actually reads under pytest's import.

**Status moved.** §10 I36c DONE; entry 48's note narrows to `_drive` under `=1` and the second solver. Reduced, not retired.

### §9 item 702 — 2026-10-31 (§10 I36d, chain I — daily, 2026-10-31, provisional, chain I's sixth; `_run_drive` end to end on a small mesh before entry 48's window; chain step (no solve))

*Written by the 2026-10-31 03:00 daily review. Depends on nothing. `xl` entry 48 runs Monday 11-02 02:00 and is QUEUED. Not daytime; no gmsh and no solve, so the `xxl` window does not block it. I36c (item 695) ran `_arm`; the rest of `_run_drive` (`test_geo36_rod_f_clinical_drive.py:472–630`: the arm order, the census, the printing, the JSON) has never run on any mesh, and it follows ≈ 750 s of window.*

**Route — insertions only, the same module:** one new `_coding` test, `test_i36d_k9_run_drive_end_to_end_on_stubs(monkeypatch, tmp_path)`, on `_k8_arm`'s mesh, field and synthetic circuit (factor out nothing: call or copy, no line of `_k8_arm` changes). It calls `_run_drive(MPI.COMM_WORLD, monkeypatch, None, 193.5996)`, and `_arm` runs real.
- **`_drive` stub:** arm M's columns are `E_u` and `E_u`; arm R's are `0.5 E_u` and `0.5 E_u`, chosen by which solver token it receives. Powers close `TH-14`'s identity as in `GEO-38.md` item 701.
- **Other stubs,** through the module the body reads (`sys.modules[_arm.__module__]`, item 695's disclosure):
  - `_route_drive` → `[(1, mmap)]`, `mmap` with `SimpleNamespace` entries carrying `SIGMA_WIRE_S_PER_M` on 403 and `SALINE_SIGMA` / `SALINE_EPSILON_R` on 3, so `_m2_failures` reads clean;
  - `c4._build_c4_fixture`, `_drive_env`, `c4._esr_specs`, `c4._z_failures`, `fem_em_solver.core.TimeHarmonicSolver`, as `GEO-38.md` item 701 specifies;
  - `_probe_points` → `_k8_arm`'s `tips, t1`; `I36_JSON_DIR` → `comm.bcast(str(tmp_path), root=0)`; `_say` → a recorder that still prints.
- `build_density_field(...)`, `global_cell_tag_set`, `col._owned_tag_counts`, `_gel_twin_map`, `_a0m_failures` and `_rss_gib` run real.

**Anchors (closed form, asserted, at `-n 1` and `-n 2`).** Arm R's field is half arm M's on both columns, so:
- (K9-Q) every `q[n]` in the JSON = **4** ≤ `EXACT_IDENTITY_RTOL`.
- (K9-E) `eta_t1` = **2** ≤ `DISCRETE_IDENTITY_RTOL`.
- (K9-A) (A0m) reads no failures, and `census_m` is the mesh's own count.
- (K9-J) the JSON parses on rank 0 with `M` and `R` blocks, `s2` equal between them at `EXACT_IDENTITY_RTOL` (same circuit). **`census_m`'s values must serialise:** if `_owned_tag_counts` returns numpy integers, `json.dump` raises here rather than in the window.
- (K9-L) the recorder holds `(P_RC)`, `(A0m)`, `tip 1 g ratio` and `price:`.
- **Negative control (asserted):** arm R's column 2 with `p_src` perturbed by 1e-3 makes `_run_drive` raise `AssertionError` naming `R (gel twin) (P) P2`, and the JSON is still written. Separation six orders; the ceiling is 1.

**Windows** (harness, bare, foreground, Bash `timeout: 660000`; `FEM_EM_GEO36_S8_DRIVE` unset):
- W0 `GEO-36-I36d-W0`: `docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/validation/test_geo36_rod_f_clinical_drive.py -v -s --tb=short'`.
- W0b `GEO-36-I36d-W0b`: the same at `-n 2`.
- W-env `GEO-36-I36d-Wenv`: `tests/environment`, `-n 1`, `timeout -k 30 300`.

**Class:** chain step (no solve), ≈ 8.8 measured / 30 predicted minutes. **Tier:** smoke, from item 695's W0 (17 s) and W0b (6 s).

**Traps already paid for.** `_arm` calls `gc.collect()` and `del` on its own fields, so the stub hands it copies. `_run_drive` does `del ctx["solver"]` before arm R: keep `ctx` a real dict. `monkeypatch` is passed in and used on `C4_DEGREE_ENV`: harmless with the fixture stubbed. Count verdicts per rank, by token. `.4f` of the node only (entry 49's note).

**Scope.** `_run_drive`'s control flow, print lines and JSON. **Not claimed:** `_drive` under `=1`, the second factorisation, any clinical reading. Entry 48's Decision rule is unchanged.

**Negative result.** As `GEO-38.md` item 701: a one-line print or JSON defect may be fixed here, first log kept, line disclosed, windows re-run once. Anything touching a formula, a gate or the drive route is **BLOCKED in the same commit**, parked on `attempt/GEO-36-I36d-*`; the slot does not touch the queue file, and the next review holds `xl-queue.d/20261102-GEO-36-step9-I36.env`. **This item must land before Monday 02:00** to protect that window; a slot that cannot land it says so in the journal. No `PREDICTED RED:` line.

**Status it can move.** §10 I36d DONE; entry 48's note narrows to `_drive` and the second factorisation, edited in the landing commit. Chain I is then at 6 of 6.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 702 — `GEO-36` I36d — per `docs/planning/chunks/GEO-36.md` "§9 item 702". Insertions only, except the one-line fix the Negative result allows. No gmsh, no solve. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome 2026-10-31 (04:30 slot): DONE, `test_i36d_k9_run_drive_end_to_end_on_stubs` green at `-n 1` and `-n 2`; `_run_drive` needed no fix, and `census_m` serialises (`_owned_tag_counts` returns Python ints).** Windows: W0-r2 `20261031T093641Z_GEO-36-I36d-W0-r2.log` (10 s), W0b `20261031T093652Z_GEO-36-I36d-W0b.log` (10 s; 12 passed, 1 skipped on each rank), W-env `20261031T093703Z_GEO-36-I36d-Wenv.log` (49 s, 12 passed). The first W0 (`20261031T093619Z_GEO-36-I36d-W0.log`) was green too; the test's own fixture lost an unused `None` slot after it, so W0 was re-run on the module as committed (rule (i)).

Readings, the same at `-n 1` and `-n 2` (W0-r2 `:113–115`, W0b `:153–155`):
- **(K9-Q):** both tips miss 4 by 0.0, against 1e-10.
- **(K9-E):** miss from 2 = 0.0, against 1e-6.
- **(K9-A):** `census_m` {3: 10260, 403: 108} equals an independent count of owned cells, summing to the mesh's global cells; (A0m) reads no failures.
- **(K9-J):** `M` and `R` blocks parse; `s2` M against R 0.0, against `_circuit_current_route_s` 0.0, at 1e-10.
- **(K9-L):** no line missing.
- **Control:** raised `['R (gel twin) (P) P2 9.990e-04']`; the JSON was still written, R `parts[1]` `p_src` 1.001, M's 1.0.

**Disclosed.** `_k8_arm` ends in a call to `_arm`, so its mesh and field are copied into `_k9_fixture`, and no line of it changed. Specs are `SimpleNamespace`s, because `_run_drive` never calls `dataclasses.replace` on one; the problem is a dataclass. (K9-J)'s check against the circuit route goes beyond the brief. Exact 0.0 readings follow from one stub field scaled by 0.5 and one circuit, so they test the body's arithmetic and not any physics.

**Status moved.** §10 I36d DONE; entry 48's first-run risk narrows to `_drive` and the second factorisation, its Status edited in this commit. Chain I is at 6 of 6.

### §9 item 707 — 2026-11-01 (§10 I32a, chain I — the clinical implant-heating body, composed and stub-run before entry 55; chain step (no solve))

*Enumerated by the 10-31 weekly. Depends on nothing. Next: `xl` entry 55. Family I32a, attempt 1 of 4. The weekly's 10-31 epitaph binds: a window body's geometry, sense and tags are named before it is queued, and its glue runs on stubs as part of the step.*

**Three execution details the weekly's text leaves open, fixed here. Each is returned to the 11-02 weekly; no target moves.**
1. **The rod is not a PEC hole on this mesh.** §10's sentence says so; the fixture keeps the rod as tag 403 cells at σ_wire (`_kwargs`, `keep_cells=True`; `_run_drive`'s "arm M (rod, tag 403 sigma_wire)", `test_geo36_rod_f_clinical_drive.py:501`). The thermal statement is unchanged: tag 403 carries no source and is outside the thermal domain, and its wall is adiabatic, I39's route.
2. **Arm letters are reversed between the two precedents.** In I36 arm M is the rod and arm R the gel twin. In I41 (`test_imp3_rod_heating_matched.py:17`) arm M is the merged gel twin and arm R the rod. Entry 55's "arm M / arm R beside I41's 1.10" means **rod / gel twin**, so print it under those words and let I36's letters stand.
3. **The thermal domain is matched as I41 matched it:** the rod arm heats tag 3 only; the gel-twin arm heats tags 3 + 403 as one gel (`_merge_rod_into_phantom`, `test_imp3_rod_heating_matched.py:87`). One G = 1 W / P_gel(gel twin) on both arms. Sense: (1, +1j), `h17.SENSES`, so the θ ratio reads beside entry 48's Q on the same drive. Entry 52's lattice reads ccw, so a θ figure here is not comparable with C6's; print the sense on every θ line.

**Route.** New module `tests/validation/test_geo36_rod_f_clinical_heating.py`, imports by name from `test_geo36_rod_f_clinical_drive.py` (`_route_drive`, `_drive_env`, `_drive`, `_gel_twin_map`, `_a0m_failures`, `_probe_points`, `_node_knob`). **That module stays byte-identical** (`_arm` frees its fields, so a new `_heat_columns(env, solver, specs, feeds)` holds the two `_drive` outputs and forms `e_sup` as `_arm` does at `:444–447`). Knob `FEM_EM_GEO36_S8_HEAT` = `1`; node knob `FEM_EM_GEO36_S8_DRIVE_C_PF` read through `_node_knob`; degree 1 only.

`_run_heat(comm, monkeypatch, c_pf)`: the build and specs as `_run_drive` `:503–530`; rod arm, then gel twin (`:550–553`), each `_heat_columns`; free both solvers (`gc.collect()`); then `_thermal_arm` on each, I39's `_heat_arm` (`test_imp3_rod_heating_phantom.py:150`) with these changes:
- `_outer_robin_tags` re-centred: a local `_dist_to_rod_axis` on `i35._rod_centre()`, half-segment `ROD_LENGTH / 2`, `WALL_SELECT_M` imported. The wall is at 1 mm from the axis and the nearest box face at 45 mm. On the merged arm there is no hole, so every exterior facet is outer: assert the wall count is 0 there.
- **MAT-21's time base:** `T_STAR` = 0.2 T_C, Δt = t\*/`N_MAIN`, `DT_1` from `test_mat21_box_fem.py:72–73`. MAT-15's `_cn`, `_be_one` and `_e_prime` close over MAT-15's own t\* (`test_mat15_coil_driven_heating_phantom.py:154–170`), so call `_transient` directly with these values, and write a local `_e_prime_box`. `H_GEL` and `K_GEL` are the box's H and K: assert equality as a precondition (`test_mat21_box_closed_form.H`, `.K`).
- **Probes:** the tip 1 g ball centres from `_probe_points`, plus P37's six points mapped into the fixture frame. `PROBE_PTS` is (x, y, z) in `MAT-21`'s frame (A, B, C on x, y, z); the fixture has Lz = 2A, Lx = 2B, Ly = 2C centred at (0, off, 0). So a point (a, b, c) maps to (b, off + c, a); name the map `_p37_fixture` and assert every mapped point lies in the box. They lie on ∂Ω except P0, and a sub-mesh read needs `parent=` (item 700's fix, `test_imp7_duty_dose_fem.py:103`).
- Print per-phase time; print θ on the sense label.

**Gates (imported, unmoved).** On each arm: (B) ∫q = G·P ≤ `B_BAND`; the scheme's balance ≤ `E_BAND`; (E′) stored + surface = G·P ≤ `E_BAND`; (S2) the calorimetric slope on the gel ≤ `S2_BAND`; (S2-tip) the same on the tip 1 g ball, ρc·ball-mean θ¹ / (Δt₁·ball-mean q) on the rod arm, ≤ `S2_BAND`; (S0) census as `_a0m_failures`, plus the merged arm's tag-3 count = rod arm's tag 3 + tag 403. I36b's EM gates: (R1), (P), NEG on both arms, exactly as `_run_drive` asserts them (the two-column pair is the same); collect every gate, then assert (item 668). **Printed:** θ at t\* at both tip balls and the six P37 points on both arms; tip θ ratio rod / gel twin beside I41's 1.10; rod share; per-phase time.

**Anchor for this slot (the stub proof, knob unset, `@_coding`).** On a `create_box` with `_drive` stubbed, the `_k8_arm` pattern (`test_geo36_rod_f_clinical_drive.py:961`: 12³ cells, inner cube tag 403, the rest tag 3):
- (K10-q) a uniform stub field E₀ x̂ on both columns gives e_sup = (1 + j) E₀ x̂, so ½σ|e_sup|² = σE₀² on every cell, ≤ `DISCRETE_IDENTITY_RTOL`, and q/ρ is the SAR `mass_averaged_sar` reads;
- (K10-E) `_transient` on that source with Robin on the outer faces: balance ≤ `E_BAND`, (E′) ≤ `E_BAND`;
- (K10-S) the slope identity on the uniform source ≤ `S2_BAND` (the source is uniform, so it holds to the scheme, far inside the band);
- (K10-C) a constant-curl stub field (E₀(−y, x, 0)/L) with (B) ≤ `B_BAND` against the closed-form ∫½σ|E|²;
- (K10-p) the fixture map: `_p37_fixture` on P0 … P5 and on three planted points;
- (K10-M) the merged arm's census identity on the small mesh;
- `_run_heat` end to end with fixture and drive stubbed, items 701 / 702's pattern, so the window's glue has run before 02:00 on 11-06.

**Negative control (asserted).** A perturbed source q(1 + 1e-3) misses (B) by 1e-3 and `B_BAND` by the stated ratio; compute that ratio from the imported band and name it. A single-column drive (e₂ = 0) misses (K10-q)'s closed form by 0.5, `_k8_arm`'s ceiling. Named controls only, no band set.

**Windows.** W0 `GEO-36-I32a-W0` (`-n 1`, complex, `timeout -k 30 120`); W0b `GEO-36-I32a-W0b` (`-n 2`, same); W-env. Knob unset throughout, no gmsh, no solve on a clinical mesh. Item 702's W0b is the size guide. **Class:** chain step (no solve), ≈ 4.5 measured / 40 predicted minutes; a composed body with a stub proof took 8–16 minutes before (items 656, 694, 702), so expect ≈ 10–15. **Tier:** smoke.

**Traps already paid for.** Item 700: a `create_submesh` mesh has no original cell index, so `evaluate_vector_field_parallel` needs `parent=`. Item 668: collect every gate before asserting. `_transient`'s time base is yours, not MAT-15's: a silent MAT-15 t\* would pass the stub and heat for 1 254 s rather than 2 821.5 s; assert `dt == T_STAR_BOX / N_MAIN` printed. UFL: no ordering comparison on a `sqrt`, the ball predicate is `ufl.lt(ufl.real(dot), ufl.real(r) ** 2)` as I39 has it. Pin `metadata={"quadrature_degree": …}` on any `SpatialCoordinate` form. Pytest's double-import: stub `_drive` through `sys.modules[_heat_columns.__module__]`. `.eval` bytes depend on the batch, so print, never compare across batches. Release fields before the thermal phase: entry 48's 42–50 GiB is the EM half alone.

**Scope.** The body and its stub proof. **Not claimed:** any clinical window, a price, a θ figure, a tissue or published comparand (7c, W3). The module's `=1` branch does not run here.

**Negative result.** A stub identity red is a defect in the composed body: fix it in the slot, keep the first log, re-run once under rule (i). A red that needs a change to `_arm` or `_drive` ⇒ stop, park on `attempt/GEO-36-I32a-*`, BLOCKED, because the drive module is entry 48's and is frozen until Monday's window has read it. No `PREDICTED RED:` line.

**Status it can move.** §10 I32a → DONE; **entry 55's prerequisite (ii) MET, and the implementer marks that in `xl-pending.md` in the landing commit** (step 6b rule 5). Entry 55 goes READY only once entry 48 is green too, and the 11-03 03:00 review writes its command and price.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 707 — `GEO-36` I32a — per `docs/planning/chunks/GEO-36.md` "§9 item 707". Stub proof only: the `=1` branch never runs here. `test_geo36_rod_f_clinical_drive.py` stays byte-identical. Mark entry 55's prerequisite (ii) MET in the landing commit. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome 2026-11-01 (04:30 slot, item 707): DONE, run inline, no executor spawned, ≈ 7 min of slot.** New module `tests/validation/test_geo36_rod_f_clinical_heating.py`; `test_geo36_rod_f_clinical_drive.py` byte-identical; no `src` edit; the `=1` branch never ran, and no gmsh ran. First run green at both widths, no fix needed. W0 (`-n 1`): 7 passed, 1 skipped, 17 s (`20261101T103903Z_GEO-36-I32a-W0.log:107`). W0b (`-n 2`): 7 passed, 1 skipped per rank, 14 `PASSED` tokens, 14 s (`20261101T103929Z_GEO-36-I32a-W0b.log:126`, `:129`). W-env: 12 passed, 48 s (`20261101T103943Z_GEO-36-I32a-Wenv.log:118`). All Status 0. Readings at `-n 1` (W0 `:48`, `:63–64`, `:67`, `:104`): (K10-q) miss 4.4e-15, SAR·ρ miss 2.2e-16; (K10-E) (B) 7.9e-13, (E) 1.5e-10, (E′) 6.5e-13, Δt 70.5375 s; (K10-S) slope 2.02e-05 and tip 8.0e-06 against `S2_BAND` 1e-4; (K10-C) rel 1.6e-15 against `B_BAND`; `_run_heat` end to end: one G, deposited power 4 V₃/V at 6.9e-15, walls 144 / 0, probes valid. Controls: q(1 + 1e-3) misses (B) by 1.000000001e-03, 1.0000e+07 × `B_BAND`, the named ratio `PERTURBED_B_RATIO`; e₂ = 0 misses (K10-q) by 0.500000000000.

**Four details for the 11-02 weekly; no target moved.** (1) `_heat_columns` takes `order, z_th, v_s` beyond the brief's signature: (R1) needs S, and S needs them. (2) The drive names are rebound in the stub (`_stubs` runs `_k9_stubs`, then copies the four stubbed names here), because this module imports them by name, as the brief says. (3) The stub sets `WALL_SELECT_M` to 30 mm: `_k8_arm`'s jagged block reaches ≈ 24 mm from the axis, against the real rod's 1 mm wall at 5 mm. So the rod wall's selection radius is exercised, its real value is not. (4) **Unmeasured until entry 55:** (S2-tip) ≤ `S2_BAND` is asserted on a source that is uniform here and steep at the real tip, and nothing has priced conduction into a 1 g ball over Δt₁ there. A red on that gate alone would be a finding about the gate. The body writes no JSON: entry 55's record is its log. **Entry 55's prerequisite (ii) marked MET** in `xl-pending.md`; (i), entry 48, stays open.

### §9 item 779 — 2026-11-03 (§10 I57, chain I — the F-small rod's tip Q by entry 48's definition; chain step (solve))

**Executes** §10 I57, with the fixture premise corrected. **The weekly's fixture cannot run as written.** Item 235's 168 507-cell mesh is the *hole* build. Tag 402 is a *facet* tag (`BIRDCAGE_IMPLANT_SURFACE_TAG`, `src/fem_em_solver/io/mesh.py:272`), so it cannot be mapped to gel.

Entry 48's definition needs three things on one mesh: arm M with the rod as σ_wire cells, arm R with the same cells mapped to gel, and (A0m), the matched census. On F-small that exists only on **I41's kept build**: `_build_kept_rod`, with cell tag 403 = `BIRDCAGE_IMPLANT_CELL_TAG` (`mesh.py:280`) and `KEPT_ROD_CELL_RECORD` 182 371 (`tests/mesh/test_birdcage_implant_rod.py:108`). That is the like-for-like fixture. The target is kept: both definitions printed side by side on one F-small mesh.

**Files.** New `tests/validation/test_imp3_tip_q_entry48_def.py`, knob `FEM_EM_IMP3_I57=1`.
- **The two arms,** both through I41's `_arm_with_403(built, label, field_hook, map_gel)` (`test_imp3_rod_heating_matched.py:104`):
  - the rod: `map_gel=False`, so 403 keeps its built σ (print that σ);
  - the gel twin: `map_gel=True`.
  - **I41's arm letters are the reverse** (I41's M is the gel twin, `:5–18`). Label arms by content ("rod", "gel twin"), never by letter.
- **The drive:** I41's four-port drive at 64 MHz, degree 1. Read `_run` (`:131`) and use its excitation unchanged.
- **The sheets:** exact rectangles, as I41 built them. Confirm `_build_kept_rod`'s sheet arguments first, and print them.

**Readouts (a record, no band).**
- **Q = SAR_1g(rod)/SAR_1g(gel twin)** on the four balls of `_probe_points(ROD_CENTRE)` (`test_geo36_rod_f_clinical_drive.py:273–292`): top, bottom, top apex and bottom apex. Use `averaging_ball_radius` + `mass_averaged_sar` (`:449–460`).
- **η_rod(T1) = ‖E_rod(T1)‖/‖E_gel(T1)‖** (`:599`), with T1 = `PTS["T1"]` (`test_imp3_palace_fields.py:66`).
- **Printed beside them:**
  - entry 48's Q 4.460568 / 4.428040 / 4.413534 / 4.341688 and η 3.642324 (`20261102T080008Z_GEO-36-step9-I36.log:12501`);
  - `IMP-3`'s η interval (3.48, 3.56) (`ENTRY43_READOUT_INTERVAL`, `test_imp3_tip_ladder_1g.py:176`), on the hole-vs-gel definition;
  - I41's hole-vs-gel Q 4.090 / 4.114 (`20261025T201613Z_IMP-3-I41-W1b.log:3776–3777`).

**Gates (imported, unmoved).**
- (A0): the cell count equals `KEPT_ROD_CELL_RECORD` exactly.
- (A0m): the census matches between arms.
- On each arm's 4×4: (i) `RECIPROCITY_BAND` 1e-3, and (ii) σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`.
- (E): `PEC_WALL_MAX_E` 1e-12, wherever I41 asserts it.
- (P): the power partition on both arms, at `DISCRETE_IDENTITY_RTOL` 1e-6.

**Controls.**
- The rod arm run twice (the same `map_gel` on both sides) reads Q = 1 exactly on all four balls. It fires if state leaks between arms.
- The weekly's off-route control: `_build_kept_rod()` alone reproduces 182 371 cells.

**Windows.** W1: `-n 2`, complex, `timeout -k 30 600`. Window form: `20261025T201613Z_IMP-3-I41-W1b.log:12`. I41's three arms took 288.5 s at `-n 2`; two arms here should take ≈ 200 s. Pre-check `pgrep -c python3` reads 0.

**Tier.** Heavy, `-n 2`. **Class:** chain step (solve). **Cost:** ≈ 6.1 measured / 35 predicted min.

**Traps already paid for.**
- The arm-letter reversal (above).
- Count pytest verdicts per rank, by token.
- The leak guard on printed S digits: on a hit, move the log unchanged to `logs/leak-held/` and report, as in `EX-177.md`.

**Scope.** One F-small fixture, degree 1, one drive. **Not claimed:** convergence, a tip-Q band, or anything clinical.

**Negative result.** A gate red: stop, BLOCKED. A Q far from both comparands is a reading, not a red.

**Status it can move.** §10 I57, the like-for-like reading that entry 48's ruling asked for, and the `GEO-36` row's I57 clause. **Depends on:** nothing.

**Outcome 2026-11-03 (12:00 UTC slot, item 779): BLOCKED, a gate red; parked on `attempt/GEO-36-I57-20261103T1203Z` (`7ffbeee1`).** Module `tests/validation/test_imp3_tip_q_entry48_def.py`, knob `FEM_EM_IMP3_I57=1`. **Premise corrected by measurement:** `_arm_with_403(map_gel=False)` does not keep a σ for 403 — the build's map has no tag 403 (I41 W1b `:3775`, re-read here `20261103T120306Z_GEO-36-I57-W1.log:1966`, "had 403 False"), so that arm is vacuum. The rod arm therefore maps 403 to `HomogeneousMaterial(SIGMA_WIRE_S_PER_M = 800, gel ε_r 78, μ_r 1)` as I35 does (`test_geo36_rod_f_clinical.py:481–483`), through a wrapper of `_arm_with_403`'s shape; the gel twin is `_arm_with_403(map_gel=True)` unchanged. Sheets printed (`:71`): `_build_rod`'s `emit_port_sheets=True`, `port_box_size=(0.01, 0.008, 0.01)`, `leg_gap_length=0.008`, 26 facets each — I41's sheets as built. W0 (`20261103T120253Z_GEO-36-I57-W0.log`): 3 passed, 6 skipped, Status 0. **W1 (`20261103T120306Z_GEO-36-I57-W1.log`, -n 2, 287 s, Status 1): 11 passed, 2 failed per rank.** Green: (A0) 182 371 cells (`:1871`); (A0m) the per-tag owned census identical on all three arms (`:1986`); (i) 2.190e-15 / 2.290e-15; (ii) σ_max 0.999729735 / 0.999723588; the repeat arm's Q = 1.0 exactly on all four balls (`:1993`). **Red: (P) on the rod arm, P1 rel 2.882e-06 against `DISCRETE_IDENTITY_RTOL` 1e-6** (`:1968`); P2–P4 3.699e-07 / 3.973e-07 / 3.702e-07. The gel twin is green but its P1 also stands out (2.853e-07 against 9.0e-09 / 2.7e-08 / 8.9e-09, `:1974–1977`), so σ_wire raises an existing P1 excess ≈ 10×. The repeat arm reproduces it bit for bit. (E) is n/a: I41 asserts `PEC_WALL_MAX_E` on no arm, and the kept build has no PEC wall. **Readings, unlicensed under the red** (`:1987–1992`): ccw, balls of 0.98 g; Q = 3.951922 top / 3.973311 bottom / 3.899322 top apex / 3.921665 bottom apex, beside entry 48's 4.460568 / 4.428040 / 4.413534 / 4.341688 and I41's hole-vs-gel 4.090 / 4.114; η_rod(T1) = 2.157532 beside entry 48's 3.642324 and `IMP-3`'s (3.48, 3.56). **Hypothesis, not diagnosed:** `_exact_shares`' volume-loss integral under-resolves σ_wire's loss in the rod's cells near P1 (az 0°, the rod at x = +15 mm), or the identity's discrete form does not carry a 1600× σ contrast at 1e-6; entry 48's (P) on F-clinical uses a different helper (`_partition`, with P_surf) and was green. **Unblock:** the review rules — diagnose the P1 excess (print P_vol by tag at two quadrature degrees on the rod arm), or re-pose the rod arm's (P) with that measurement. Never a re-band in-slot.

### Review ruling 2026-11-03 10:30 (item 779 BLOCKED — the material confirmed, the (P) red re-posed as item 801)

1. **The rod arm's material is entry 48's definition, so the slot's mapping stands.** It maps tag 403 to `HomogeneousMaterial(SIGMA_WIRE_S_PER_M, gel ε_r, gel μ_r)`. The brief's own premise paragraph already names it ("arm M with the rod as σ_wire cells", line 1279 above). I35 maps the rod the same way (`test_geo36_rod_f_clinical.py:481–483`), and I36a's (M2) asserts it on entry 48's own maps (this file, the I36a brief). Only the brief's helper premise was wrong: `_arm_with_403(map_gel=False)` keeps no σ, and the map has no 403 (`20261103T120306Z_GEO-36-I57-W1.log:1966`). The executor replaced that helper call and kept the target. The slot note's "in-slot material" flag is answered: no improvised target.
2. **The (P) red stands as measured: P1 2.882e-06 against `DISCRETE_IDENTITY_RTOL` 1e-6 (`:1968`).** No band moves.
   - The identity is form-level (`_exact_shares`, `tests/validation/test_birdcage_power_identity.py:313`). On degree-1 Nédélec with a DG0 σ, the loss integrand σ|E|² is quadratic per cell, so a degree-2 rule is exact. **A quadrature cause is therefore unlikely, though not ruled out:** the form's estimated degree is not printed.
   - The alternative is the algebraic solve's floor at a σ contrast of ≈ 1.6e3 against the gel, with ωε₀ small beside 800 S/m. The gel twin's P1 is also ≈ 30× its other columns (`:1974–1977`), which fits a P1-local conditioning effect that σ_wire amplifies.
   - Item 801 measures both readings and decides between them.
3. **The family count:** I57 is at attempt 2 of 4. The unlicensed readings (Q 3.95–3.97, η 2.158) stay unlicensed until (P) is disposed of.

### §9 item 801 — 2026-11-03 (§10 I57, chain I — the rod arm's (P) red separated into quadrature or solve floor, attempt 2; chain step (solve), a measurement)

**Executes** the ruling above. It does not close I57. It decides which reading the (P) red is.

**Code.** Take `tests/validation/test_imp3_tip_q_entry48_def.py` from `attempt/GEO-36-I57-20261103T1203Z` (`7ffbeee1`) by path checkout. Add one knob-selected test, `FEM_EM_IMP3_I57_DIAG=1`, that runs **the rod arm only** (no gel twin and no repeat arm) on I41's kept build. Use the same material wrapper, the same drive and the same four single-port drives. Print, per drive:
- **(q) quadrature.** `_loss_power_w`'s two terms, recomputed by tag (3 = gel, 403 = rod, and the rest), with `metadata={"quadrature_degree": d}` for d = 2 and d = 4. Print beside them the unpinned form's value and its estimated degree (`ufl.algorithms.estimate_total_polynomial_degree`). Print the identity residual (P_src − Σ sheet − P_vol)/P_src at each d.
- **(r) the algebraic residual.** ‖A e − b‖ / ‖b‖ on the assembled complex system of each drive's solve, through the module's `LinearProblem`'s `A` / `b` (PETSc `mult`, all-reduced norms). If that operator cannot be reached within 15 minutes of reading the drive's code path, print `(r) n/a` and say why. Do not rebuild the solver.

**Gates.** All imported and unmoved, re-asserted on the rod arm: (A0) 182 371 cells, (i) `RECIPROCITY_BAND`, (ii) σ_max, and (P) at `DISCRETE_IDENTITY_RTOL` with the unpinned form. That last one is pre-registered red:

PREDICTED RED: (P) rod arm P1 range [2.5e-06, 3.3e-06]

**Control.** The unpinned form's P1 residual reproduces W1's 2.882e-06 to three digits. If it does not, the run is not the parked run: stop.

**Decision rule.** Write it into the outcome. The review rules; the slot does not.
- **Reading Q:** the d = 4 residual is ≤ 1e-6 on P1, and d = 2 / unpinned reproduce the red. The next review queues I57 with the loss form's quadrature pinned at the resolving degree on both arms. The identity is then evaluated exactly, and the band is unmoved. The 11-04 weekly ratifies.
- **Reading S:** the residual is flat in d (P1 equal to two digits), and (r) on P1 is ≥ 1e-7. That is a solve floor at this contrast. The review files a known-issues entry. The rod arm's (P) is the 11-04 weekly's to rule, because it is a band question.
- **Neither:** report the numbers and stop. The question goes to the weekly.

**Windows.** One W1 at `-n 2`, complex, `timeout -k 30 600`, bare and foreground. Window form `20261103T120306Z_GEO-36-I57-W1.log:12`, with the knob changed. Measured: the build 50.2 s and one arm 75–82 s (item 779), so ≈ 150–250 s with the extra assemblies. Pre-check `pgrep -c python3` must read 0.

**Tier.** Heavy, `-n 2`. **Class:** chain step (solve). **Cost:** ≈ 6.1 measured / 25 predicted min.

**Traps already paid for.**
- Arms are labelled by content.
- Count pytest verdicts per rank, by token.
- Apply the leak guard on S digits.
- A facet or volume form with a `SpatialCoordinate` needs its quadrature degree pinned. The (q) forms are pinned by construction.
- Norms and integrals are rank-local until all-reduced.

**Scope.** One arm, one mesh, degree 1. **Not claimed:** Q, η, or I57's closure.

**Status it can move.** Which reading the red is, and from that the `GEO-36` row's I57 clause. On reading Q, I57's re-run is queued by the next review. On reading S, a known-issues entry and the 11-04 weekly's ruling. **Depends on:** nothing.

**Outcome 2026-11-03 (13:30 CST slot, item 801, taken next at minute 20 after item 791): BLOCKED, predicted red inside its range, all else green.** Parked on `attempt/GEO-36-I57-diag-20261103T1952Z` (`c1efe280`). The branch holds the module, taken from `7ffbeee1` by path with the knob-selected `FEM_EM_IMP3_I57_DIAG=1` test added, plus the log and its test-results row.
- **The tests:** the rod arm only, with the same material wrapper and the same drive. (q) is assembled by tag at `quadrature_degree` 2 and 4 and unpinned. (r) wraps `TimeHarmonicSolver.solve` / `solve_with_held_factorization` and reads `‖A e − b‖/‖b‖` on the held `LinearProblem`'s `A` / `b`; PETSc norms are global. No solver was rebuilt.
- **W1:** `20261103T195044Z_GEO-36-I57-diag-W1.log`, `-n 2`, 124 s, Status 1. 9 passed, 1 failed, 6 skipped per rank. The failure is `test_diag_p_rod_arm`, the pre-registered red.
- **Gates:**
  - (A0) 182 371 cells (`:1879`).
  - (i) 2.190e-15 and (ii) σ_max 0.999729735; column powers ≤ 1 (`:1911`).
  - (P) P1 **2.882e-06, inside the PREDICTED RED range [2.5e-06, 3.3e-06]**; P2–P4 3.699e-07 / 3.973e-07 / 3.702e-07 (`:1912–1915`).
  - **Control:** P1 \|resid\| 2.8818e-06 reproduces W1's 2.882e-06 to three digits, so this is the parked run (`:1942`).
- **(q):** flat in d to every printed digit, on every drive (`:1919–1921`). The estimated degree is 2, and the unpinned, d = 2 and d = 4 values agree in all ten digits on every tag.
- **(r):** 6.444e-13 (factorise), then 6.343e-13 / 6.211e-13 / 6.393e-13 (held) (`:1922`, `:1927`, `:1932`, `:1937`).
- **The red is the identity's own accounting.**
  - `_loss_power_w` (`test_port_drive_superposition.py:205–221`) sums P_vol over `PHANTOM_CELL_TAG` (3) and `CONDUCTOR_CELL_TAG` (1) only, so **tag 403's loss is not in it**.
  - On P1, tag 403 dissipates 1.505772888e-08 W against P_src 5.225076965e-03 W. That ratio is 2.882e-06, exactly the residual.
  - With P_vol over **all cells**, the residual is −6.640e-16 / +9.960e-16 on P1 and ≤ 9.6e-15 on every drive and degree (`:1919–1936`).
  - "Rest" (all − 3 − 403) equals the conductor tag to ten digits, so no other lossy tag exists.
  - P1's excess is P1-local because tag 403's loss is ≈ 8× larger under P1 than under P2–P4 (1.506e-08 vs 1.93–2.08e-09 W). The rod sits at x = +15 mm, nearest P1 at az 0°.
  - The gel twin's P1 excess in item 779 (2.853e-07) is the same omission at σ_gel: 403's loss there is in neither term.
- **Decision rule, read:** **neither.** It is not reading Q, because the residual is flat in d and d = 4 does not resolve it. It is not reading S, because (r) is 6.4e-13, below the 1e-7 the rule needs. Per the rule, the numbers are reported and the question goes to the weekly.
- **The slot's reading for the review, not a ruling:** the (P) identity is green to ≤ 1e-14 once P_vol includes tag 403. Entry 48's (P) passed on F-clinical because its volume term is the all-cells integral (`p_total = _vol(dx)`, `test_geo36_rod_f_clinical_drive.py:381`; `p_rest` contains `p_rod`). The fix would put tag 403 into the rod arm's P_vol, a change to *what is summed*, with the band unmoved. Whether that is an in-letter correction of I57's (P), or a band question, is the review's call.
- **Not claimed:** Q, η, or I57's closure. Item 779's readings (Q 3.95–3.97, η 2.158) stay unlicensed until the review disposes of (P).

### Review ruling 2026-11-03 18:00 (item 801 — the (P) red is the identity's omitted tag; I57 re-posed as item 804)

**Reading, from `20261103T195044Z_GEO-36-I57-diag-W1.log` (on `attempt/GEO-36-I57-diag-20261103T1952Z`, `c1efe280`).** Neither pre-written reading holds, and the measurement names a third: (q) is flat in d in all ten printed digits and the estimated degree is 2 (`:1919–1921`); (r) is 6.2–6.4e-13 on all four drives (`:1922–1937`). The residual is **exactly tag 403's loss**: `_loss_power_w` (`tests/validation/test_port_drive_superposition.py:205–221`) sums P_vol over tags 3 and 1 only; on P1, 403 dissipates 1.505772888e-08 W of P_src 5.225076965e-03 W = 2.882e-06, the residual; over all cells the identity closes to ≤ 9.6e-15 on every drive and degree (`:1919–1936`); "rest" = the conductor tag to ten digits, so no other lossy tag exists.

**Ruling.**
1. **The (P) gate is the power partition, and a partition sums every lossy cell.** Item 779's brief imports (P) "on both arms at `DISCRETE_IDENTITY_RTOL`" as entry 48's identity, and entry 48's form is the all-cells volume integral (`p_total = _vol(dx)`, `test_geo36_rod_f_clinical_drive.py:381`, with `p_rod` inside `p_rest`). The helper's two-tag sum was written for fixtures whose only lossy volumes are tags 3 and 1; on I41's kept build it is a helper premise that does not hold, as `_arm_with_403`'s material premise did not (10:30 ruling 1). Correcting *what is summed* to entry 48's form keeps the target and the band: `DISCRETE_IDENTITY_RTOL` 1e-6, imported, unmoved. **This is the review's reading of the brief, not a band ruling;** it is named for the 11-04 weekly to ratify or revert (dashboard).
2. **The decision rule's "neither" branch is honoured:** the numbers are reported above and in the journal, and the question as posed (Q vs S) goes to the weekly as answered — "neither; an accounting omission".
3. **I57 is re-posed as item 804, attempt 3 of 4**: item 779's full window (rod arm, gel twin, repeat arm) with (P) over all cells on every arm, per-tag P_vol printed.
4. **A disclosure, filed in known-issues (2026-11-03):** `_exact_shares` / `_loss_power_w` under-sum on any mesh with lossy cells outside tags 3 and 1. No test on `main` is red; its callers on kept-rod or kept-disc builds have not been audited here.

### §9 item 804 — 2026-11-03 (§10 I57, chain I — the F-small rod's tip Q by entry 48's definition, attempt 3: (P) over all cells; chain step (solve))

**Executes** §10 I57 under the ruling above. Item 779's brief binds wherever this one is silent (fixture, arms by content, drive, readouts, gates, controls, traps); the 10:30 ruling's material mapping stands (403 → `HomogeneousMaterial(SIGMA_WIRE_S_PER_M, gel ε_r, gel μ_r)` on the rod arm; `_arm_with_403(map_gel=True)` on the gel twin).

**Code.** Take `tests/validation/test_imp3_tip_q_entry48_def.py` from `attempt/GEO-36-I57-diag-20261103T1952Z` (`c1efe280`) by path checkout — code only. Copy the attempt branches' logs — `20261103T120253Z_GEO-36-I57-W0.log` and `20261103T120306Z_GEO-36-I57-W1.log` from `attempt/GEO-36-I57-20261103T1203Z`, `20261103T195044Z_GEO-36-I57-diag-W1.log` from the diag branch — by path as the records of attempts 1–2 and **append their test-results rows by hand**; `git diff --stat` on `test-results.md` shows insertions only. Edits, and only these:
- **(P) on every arm:** P_vol = the all-cells loss integral (the diag test's all-cells form, which closed to ≤ 9.6e-15), asserted at `DISCRETE_IDENTITY_RTOL`. Do not edit `_loss_power_w` or `_exact_shares` in `test_port_drive_superposition.py` — other modules import them. Compute the all-cells term in the I57 module.
- **Print per drive and arm:** P_vol by tag (3, 403, the conductor) and the two-tag sum's residual beside the all-cells one, so the omitted term stays visible.
- The `FEM_EM_IMP3_I57_DIAG` test stays in the module, skipped unless its knob is set.

**Windows.** W1 only: `-n 2`, complex, `FEM_EM_IMP3_I57=1`, bare, foreground, `timeout -k 30 600`, Bash timeout 660000 ms; window form `20261103T120306Z_GEO-36-I57-W1.log:12`. Measured 287 s (item 779). Pre-check `pgrep -c python3` reads 0.

**Gates (imported, unmoved).** (A0) 182 371 = `KEPT_ROD_CELL_RECORD`; (A0m) the census matched between arms; (i) `RECIPROCITY_BAND` 1e-3; (ii) σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`; (P) all cells at `DISCRETE_IDENTITY_RTOL` on every arm. (E) n/a on the kept build, as item 779 found (no PEC wall) — say so. **No `PREDICTED RED:` line.**

**Controls.** Item 779's: the repeat arm reads Q = 1 exactly on all four balls; `_build_kept_rod()` alone reproduces 182 371 cells. And the ruling's: the two-tag (P) residual on the rod arm's P1 reproduces 2.882e-06 to three digits (asserted — `20261103T195044Z_GEO-36-I57-diag-W1.log:1942` is the prior measurement on this comparison and fixture), so the omitted term is the one this edit adds.

**Readouts (a record, no band).** As item 779: Q on the four balls and η_rod(T1), printed beside entry 48's 4.460568 / 4.428040 / 4.413534 / 4.341688 and 3.642324, `IMP-3`'s (3.48, 3.56) and I41's hole-vs-gel 4.090 / 4.114. Item 779 read Q 3.95–3.97 and η 2.158 under the red; a reading that reproduces them is expected, not a gate.

**Tier.** Heavy, `-n 2`, complex. **Class:** chain step (solve). **Cost:** ≈ 6.1 measured / 20 predicted min.

**Traps already paid for.** Item 779's and item 801's: arms labelled by content; `_arm_with_403(map_gel=False)` keeps no σ; per-rank verdicts counted by token; the leak guard on S digits; path checkouts take code only.

**Scope.** One F-small fixture, degree 1, one drive. **Not claimed:** convergence, a tip-Q band, anything clinical; the like-for-like statement is a record.

**Negative result.** A red: stop, park, BLOCKED. A Q or η far from both comparands is a reading, not a red.

**Status it can move.** §10 I57 DONE (the like-for-like tip-Q record entry 48's ruling asked for); the `GEO-36` row's I57 clause; the known-issues 2026-11-03 entry gains its first corrected caller. **Depends on:** nothing.

**Outcome 2026-11-04 (19:30 CST slot, item 804, started 01:30 UTC): BLOCKED, a control red — the repeat arm's Q = 1 *exactly*; every gate green.** Parked on `attempt/GEO-36-I57-20261104T0132Z` (`af5c0458`: the module from `c1efe280` with (P) over all cells on every arm and the new two-tag control, `20261104T013204Z_GEO-36-I57-W1c.log`, the diag log carried by path, both test-results rows). Pre-check `pgrep -c python3` 0. **W1c** (`-n 2`, complex, 291 s, Status 1): **13 passed, 1 failed, 3 skipped per rank** (two identical summary lines, `:2036`, `:2048`). Green: (A0) 182 371 (`:1871`); (A0m) the owned census identical on all three arms (`:1997`); (i) 2.116e-15 / 2.290e-15 / 1.901e-15 (rod / gel twin / repeat); (ii) σ_max 0.999729735 / 0.999723588 / 0.999729735; **(P) all cells ≤ 1.5e-15 / 8.3e-15 / 7.6e-15 on every drive** (`:1967–1995`); the **two-tag control: rod P1 2.8818e-06 reproduces 2.882e-06** (`:2010`), tag 403 dissipating 1.505772888e-08 W of P_src 5.225076965e-03 W on P1 (`:1968`); the gel twin's two-tag P1 2.853e-07 is the same omitted 403 term (1.490e-09 W). (E) n/a: no PEC wall on the kept build. **Red:** `test_control_rod_twice_q_is_one` — repeat / rod Q = 0.9999999999997838 top, 1.0000000000002625 bottom, 0.999999999999793 top apex, 1.0000000000002354 bottom apex (`:2004`, `:2031`), against `== 1.0`; item 779 read 1.0 exactly. The repeat arm's per-tag losses equal the rod arm's in all ten printed digits, but its (i) differs (1.901e-15 vs 2.116e-15), and the rod arm's own (i) differs from item 779's (2.190e-15) — so the solve is not bit-reproducible arm to arm at roundoff in this window. **Readings, unlicensed under the red** (unchanged from item 779 to six digits, `:1999–2003`): Q = 3.951922 / 3.973311 / 3.899322 / 3.921665 (balls 0.98 g), η_rod(T1) = 2.157532, beside entry 48's 4.460568 / 4.428040 / 4.413534 / 4.341688 and 3.642324, `IMP-3`'s (3.48, 3.56), I41's 4.090 / 4.114. **Hypothesis, not diagnosed:** the control's "exactly" is not a property of this solve path at `-n 2` (a roundoff-level arm-to-arm difference, ≈ 2.6e-13 on Q, whose source — the direct solver's parallel reduction order, or state left by the added per-drive assemblies — is not separated). **Unblock:** the review rules whether the repeat control may read |Q − 1| at a roundoff bound with this measurement recorded (a control-bound question, not a slot call), or queues a diagnosis. Never a re-band in-slot. I57 is at attempt 3 of 4.

### Review ruling 2026-11-04 03:00 (item 804 BLOCKED — the repeat control's bound; I57 re-posed as item 810)

**Reading, from `20261104T013204Z_GEO-36-I57-W1c.log` (on `attempt/GEO-36-I57-20261104T0132Z`, `af5c0458`).** Every gate is green: (A0) 182 371 (`:1871`); (A0m) matched on all three arms (`:1997`); (i) ≤ 2.290e-15; (ii) σ_max 0.999729735 / 0.999723588; (P) over all cells ≤ 8.3e-15 on every arm and drive (`:1967–1995`); the ruling's two-tag control reproduces 2.882e-06 (`:2010`). The one red is item 779's repeat control, `test_control_rod_twice_q_is_one`, written as `== 1.0`: it read 1 − 2.16e-13 / 1 + 2.63e-13 / 1 − 2.07e-13 / 1 + 2.35e-13 on the four balls (`:2004`, `:2031`). Item 779 had read exactly 1.0, and W1c's rod arm (i) differs from item 779's (2.116e-15 against 2.190e-15) on the same mesh and drive. So the solve path is not bit-reproducible arm to arm at `-n 2`, at the 1e-13 level. That is a round-off statement about the solve, not a state leak: the per-tag losses agree in all ten printed digits.

**Ruling.**
1. **The control's "exactly" was a bound no rule backs.** The 18:00 standing rule says a reproduction control names its numeric bound, or the print whose last digit sets it. Item 779 wrote "Q = 1 exactly" from one run's digits. **The control is re-registered at the module's own identity band, factor 1** (the rule (e) reminder, 2026-11-03 18:00): max over the four balls of |Q_repeat / Q_rod − 1| ≤ `DISCRETE_IDENTITY_RTOL` (1e-6, imported in the module from `test_birdcage_power_identity`, the (P) gate's band). The measured 2.63e-13 is printed beside it, with item 779's 0.
2. **The separation is computed, not assumed.** The leak this control guards against is the gel twin's material carried into the rod arm, or the reverse. Its signature is the rod-vs-gel contrast itself: Q ≈ 3.9, so |Q − 1| ≈ 2.9. The bound sits ≈ 3.8e6 × above the measured round-off and ≈ 2.9e6 × below that ceiling. No gate band moves; (A0), (A0m), (i), (ii) and (P) are unchanged.
3. **Not diagnosed, and not needed for the gate:** the source of the arm-to-arm difference, whether MUMPS's parallel reduction order or state from the added per-drive assemblies. It is printed, not asserted.
4. **Family count: item 810 is attempt 4 of 4 in I57's family.** Any red other than an in-range predicted one freezes the family. The weekly then re-scopes I57 as a new numbered step or banks the negative (§9 rule (2)).
5. **Named for the 11-04 weekly to ratify or revert,** beside the 18:00 all-cells reading of (P).

### §9 item 810 — 2026-11-04 (§10 I57, chain I — the F-small rod's tip Q by entry 48's definition, attempt 4 of 4: the repeat control at the identity band; chain step (solve))

**Executes** §10 I57 under the ruling above. Where this brief is silent, item 804's brief binds, then item 779's beneath it: fixture, arms by content, the material mapping, drive, readouts, gates, (P) over all cells, the two-tag control, traps.

**Code.** Take `tests/validation/test_imp3_tip_q_entry48_def.py` from `attempt/GEO-36-I57-20261104T0132Z` (`af5c0458`) by path checkout, code only. Copy the attempts' logs that are not on `main` by path: `20261103T195044Z_GEO-36-I57-diag-W1.log` and `20261104T013204Z_GEO-36-I57-W1c.log`, both from that branch. **Append their test-results rows by hand**; `git diff --stat` on `test-results.md` shows insertions only. Make these edits and no others:
- `test_control_rod_twice_q_is_one`: assert `max_n |q_ctl[n] − 1| ≤ DISCRETE_IDENTITY_RTOL`, with the per-ball deviation in the failure message.
- Print one line: the measured max deviation, beside "item 779 0.0, item 804 2.63e-13 (`20261104T013204Z_GEO-36-I57-W1c.log:2004`)".
- In the docstring, label the control *asserted at the (P) band, factor 1* and cite this ruling.
- Rename nothing.

**Windows.** W1 only: `-n 2`, complex, `FEM_EM_IMP3_I57=1`. Bare and foreground; `timeout -k 30 600` container-side; Bash timeout 660000 ms. Window form: `20261104T013204Z_GEO-36-I57-W1c.log:12`. Measured: 291 s. Pre-check: `pgrep -c python3` reads 0.

**Gates (imported, unmoved).** As item 804: (A0) 182 371, (A0m), (i) `RECIPROCITY_BAND`, (ii) σ_max, and (P) over all cells at `DISCRETE_IDENTITY_RTOL` on every arm. (E) is n/a (no PEC wall); say so. **No `PREDICTED RED:` line.**

**Controls (asserted).**
- The repeat control at `DISCRETE_IDENTITY_RTOL` (ruling 1). Measured 2.63e-13; the ceiling for a leak is ≈ 2.9.
- The two-tag (P) residual on the rod arm's P1 reproduces 2.882e-06 to three digits.
- `_build_kept_rod()` alone reproduces 182 371 cells.

**Readouts (a record, no band).** As item 779: Q on the four balls and η_rod(T1). Print them beside entry 48's 4.460568 / 4.428040 / 4.413534 / 4.341688 and 3.642324, beside `IMP-3`'s (3.48, 3.56), and beside I41's 4.090 / 4.114. Item 804 read Q 3.951922 / 3.973311 / 3.899322 / 3.921665 and η 2.157532. Reproducing those values is expected; it is not a gate.

**Tier.** Heavy, `-n 2`, complex. **Class:** chain step (solve). **Cost:** ≈ 11.2 measured / 20 predicted min.

**Traps already paid for.** Items 779, 801 and 804:
- Label arms by content.
- `_arm_with_403(map_gel=False)` keeps no σ.
- `_loss_power_w` sums tags 3 and 1 only. The module's all-cells term is the gate, and `test_port_drive_superposition.py` is never edited.
- Count verdicts per rank, by token.
- Apply the leak guard on S digits.
- Path checkouts take code only.
- **An exact float equality on a solve output is not reproducible at `-n 2`** (this ruling). No new equality assertion is added.

**Scope.** One F-small fixture, degree 1, one drive. **Not claimed:** convergence, a tip-Q band, anything clinical, or the arm-to-arm round-off's source.

**Negative result.** Any red: stop, park on `attempt/GEO-36-I57-*`, BLOCKED. That is the family's fourth attempt, so the family freezes and goes to the weekly (ruling 4). A Q or η far from both comparands is a reading, not a red.

**Status it can move.** §10 I57 DONE: the like-for-like tip-Q record that entry 48's ruling asked for. The `GEO-36` row's I57 clause. The known-issues 2026-11-03 entry gains its first corrected caller. **Depends on:** nothing.

**Outcome 2026-11-04 (04:30 CST slot, item 810, window started 10:48 UTC): DONE — every gate and every control green; §10 I57 DONE as a record.** Module path-checked out from `af5c0458` (code only); edits only the ruling's: the repeat control asserts max_n |Q_repeat/Q_rod − 1| ≤ `DISCRETE_IDENTITY_RTOL`, prints its deviation beside item 779's 0.0 and item 804's 2.63e-13, docstring labels it asserted at the (P) band, factor 1. The carried logs `20261103T195044Z_GEO-36-I57-diag-W1.log` and `20261104T013204Z_GEO-36-I57-W1c.log` are on `main` by path, their test-results rows appended by hand (insertions only). Pre-check `pgrep -c python3` 0, and 0 after. **W1d** `20261104T104820Z_GEO-36-I57-W1d.log` (`-n 2`, complex, 353 s harness / 350.61 s pytest, Status 0): **14 passed, 3 skipped per rank** (two identical summary lines, `:2028`, `:2031`). (A0) 182 371 (`:1871`); (A0m) the owned census identical on all three arms (`:1997`); (i) 2.190e-15 / 2.290e-15 / 2.190e-15 (rod / gel twin / repeat; `:1966`, `:1976`, `:1986`); (ii) σ_max 0.999729735 / 0.999723588 / 0.999729735; (P) over all cells ≤ 8.300e-15 on every arm and drive (`:1968–1994`); (E) n/a, no PEC wall on the kept build (`:1996`). Controls: two-tag rod P1 2.8818e-06 reproduces 2.882e-06 (`:2010`); repeat control max |Q − 1| = 0.000e+00 against 1e-6 (`:2012`) — this window the repeat was bit-identical (rod and repeat (i) both 2.190e-15, item 779's value), so the arm-to-arm round-off item 804 saw did not recur; still not diagnosed. **Disclosed:** the fixture's own summary print at `:2004` still reads "(ASSERTED exact)" — a label left from item 779 in non-gate code the brief did not list for editing; the assertion is the `:2012` one at the identity band. **Record (no band):** Q = 3.951922 / 3.973311 / 3.899322 / 3.921665 (top / bottom / top apex / bottom apex, balls 0.98 g, `:1999–2002`), η_rod(T1) = 2.157532 (`:2003`), equal to item 804's to six digits; beside entry 48's 4.460568 / 4.428040 / 4.413534 / 4.341688 and 3.642324, `IMP-3`'s (3.48, 3.56), I41's 4.090 / 4.114. Not claimed: convergence, a tip-Q band, the round-off's source, anything clinical.

### Review ruling 2026-11-05 03:00 — I67's fixture is I41's matched pair; C32's body

**I67, "entry 48's definition on `IMP-3`'s PEC-hole rod build".** A hole build (`_build_rod`, 167 371 cells) has no rod cells, so entry 48's gel twin ("the same cells mapped to gel") cannot exist on it. The like-for-like fixture is **I41's matched pair**:
- **rod arm:** `_rb._arm(_build_rod(), rod=True)`, the PEC wall on tag 402;
- **gel twin:** I57's `_gel_seen` on `_build_kept_rod()`, 182 371 cells, with 403 → gel.

The two share their tag-3 cells exactly (I41's (A0m)). I41 already printed the **apex** Q on that pair, 4.090342 / 4.113695 (`20261025T201613Z_IMP-3-I41-W1b.log:3776–3777`), so that becomes a reproduction control. The new content is the two tip-centre balls and η_rod(T1) on entry 48's readout. The target is kept.

**Disclosed for the 11-07 weekly.** `IMP-3`'s (3.48, 3.56) was read on entry 43's tip-refined ladder, not on this default hole build. So the 11-04 weekly's "≈ 40 % rod-model lever" mixes the rod model with the mesh, and I67 separates the two. On one definition the apex Q lever is already in the logs: 4.090342 / 3.899322 = 1.049.

### §9 item 822 — 2026-11-05 (§10 C32, chain C — `GEO-36` step 11, I36-3T's body: I36b's arm M / arm R route on C20's 3 T rod mesh; chain step (no solve))

**Executes** §10 C32. The target is kept: tag 403 kept as σ_wire cells, 128 MHz, the node knob, and the two-port drive on C16's feeds P33 / P41 (`c4._feeds_32`, `esr_root.py:155`, the same at both field strengths), stub-run as I36c / I36d were. `GEO-38.md`'s "Review ruling 2026-11-05 03:00" covers the chain's other steps. Family C32, attempt 1 of 4.

**Code.** A new module, `tests/validation/test_geo36_rod_f_clinical_3t_drive.py`. I36b, C21, C20 and `src` stay byte-identical.
- **Knobs:**
  - `FEM_EM_GEO36_S11_DRIVE` = `1`: unset skips; anything else raises `ValueError`.
  - The node, `FEM_EM_GEO36_S8_DRIVE_C_PF`, read through `i36._node_knob`. It is required under `1`.
  - Degree 1 only: `FEM_EM_GEO36_S8_DRIVE_DEGREE` is in the foreign list.
- **Foreign guard:** I36b's `FOREIGN_ENVS` (`drive.py:68–76`), plus `FEM_EM_GEO36_S8_DRIVE`, `FEM_EM_GEO36_S8_DRIVE_DEGREE`, `FEM_EM_GEO38_C6`, `FEM_EM_GEO38_C21`, `FEM_EM_GEO38_C20` and the C11 / C12 knobs.
- **`_run_drive_3t(comm, monkeypatch, c_pf)`:**
  - Apply `c21._c21_patches(monkeypatch)`: C20's mesh route, C20's records, and 128 MHz on `wf7` / `c4`.
  - Patch 4 wraps `box._box_readouts`, which I36 never calls. It is inert; disclose that.
  - Then call `i36._run_drive(comm, monkeypatch, None, c_pf)`.
- **Print-only constants:** entry 48's Q 4.341688–4.460568 (`20261102T080008Z_GEO-36-step9-I36.log:12501`) and η_rod(T1) 3.642, printed beside the 3 T values. I36b's "tuned on C4's cylinder" print is a stale label at 3 T; disclose it and do not edit it.

**Windows.** Bare and foreground, Bash timeout 660000, knobs unset. Window form: `20261031T093652Z_GEO-36-I36d-W0b.log:12`.
- **W0** `GEO-36-C32-W0`:
  `scripts/testing/run_and_log.sh GEO-36-C32-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/validation/test_geo36_rod_f_clinical_3t_drive.py tests/validation/test_geo36_rod_f_clinical_drive.py -v -s --tb=short'"`
- **W0-n1:** the same command at `-n 1`.

**Anchors / K-gates.**
- **K1:** the knobs and the foreign list.
- **K2:** `1` without the node knob stops before gmsh (C21's K2 pattern, `3t_box_loaded.py:214`).
- **K7-3T:** `1` plus node 44.0615 runs `_run_drive_3t` into C21's `_builder_sentinel`, and the sentinel sees:
  - rung Ø 0.704;
  - C20's rod `implant_rod` with `keep_cells=True`;
  - f 128e6 on both globals;
  - records 1 128 430 / (1 508 256, 7 783 028).
  - Also (Bore): `c20._bore_gap(i35._kwargs())` > 0 under the patch.
- **K9-3T:** I36d's K9 (`drive.py:1218`) on `_k9_stubs`, with the 3 T overrides applied after the call: ω = 2π·128e6, z_rc at 44.0615 pF and 128 MHz, and `_drive` on the 128 MHz loads.
  - (K9-Q): q = 4 ≤ `EXACT_IDENTITY_RTOL`.
  - (K9-E): η = 2 ≤ `DISCRETE_IDENTITY_RTOL`.
  - (K9-A): the census.
  - (K9-J): the JSON at `i36-node-44.0615.json`. M's and R's s2 are equal, and equal to the 128 MHz circuit route, ≤ `EXACT_IDENTITY_RTOL`.
  - (K9-L): the heads.

**Control (asserted).**
- **Arm R's power:** with arm R column 2's P_src × (1 + 1e-3), the run raises naming `R (gel twin) (P) P2`, and the JSON is still written. I36d read 9.990e-04 against the 1e-6 band.
- **K7-3T negative:** I36b's `_run_drive` without C21's patches reaches rung Ø 0.710 / 64 MHz / 1 217 025 cells (`20261103T103420Z_GEO-38-C21-W0.log:83`). The separation is categorical.
- **The frequency discriminator:** s2 against the 64 MHz route misses by 2.515e-01 on this seed-656 circuit at node 44.0615 (`20261103T103420Z_GEO-38-C21-W0.log:149`). It is asserted > `EXACT_IDENTITY_RTOL`; the log backs the size.

**Tier.** Smoke, `-n 2` (plus `-n 1`). Measured: I36d W0b 10 s, C21 W0 14 s. W0 is predicted at ≈ 25 s. **Class:** chain step (no solve). **Cost:** ≈ 6.3 measured / 25 predicted min.

**Traps already paid for.**
- `_f_clinical_kwargs` is bound by name (item 775).
- `_k8_arm` and `_k9_stubs` hard-code 64 MHz and 193.5xxx pF (`drive.py:998`, `:1136`, `:1144`). Override after the call; never edit them.
- `_arm` frees its fields, so hand it copies.
- `_run_drive` deletes `ctx["solver"]`, so keep `ctx` a real dict.
- Stub through `sys.modules[_arm.__module__]`.
- Count verdicts per rank, by token.
- Print the node at `.4f` only.
- Box containment is not asserted on this route.

**Scope.** The body and its stub proof. **Not claimed:** any 128 MHz tip Q or rod reading (that is entry 61), the second factorisation, or a published comparand.

**Negative result.** Any K red: stop, park on `attempt/GEO-36-C32-*`, BLOCKED. No `PREDICTED RED:` line.

**Status it can move.**
- §10 C32 becomes DONE.
- **Entry 61's prerequisite (i) is MET:** the implementer marks it in the landing commit, on the Status line only, with the module and knob names. READY still needs entry 57 green.
- The `GEO-36` row gains the body.

**Depends on:** nothing.

**Outcome 2026-11-05 (13:30 slot, item 822): DONE — every K-gate and every control green; §10 C32 DONE; entry 61's prerequisite (i) MET.** New module `tests/validation/test_geo36_rod_f_clinical_3t_drive.py` (I36b, C21, C20 and `src` byte-identical); knob `FEM_EM_GEO36_S11_DRIVE` (unset skips, `1` runs, else `ValueError`), node through `i36._node_knob`; foreign guard = I36b's 7 + 10 (17). **W0** `20261105T133319Z_GEO-36-C32-W0.log` (`-n 2`, 20 s harness / 17.65 s pytest, Status 0): **17 passed, 2 skipped per rank** (34 PASSED / 4 SKIPPED tokens; two identical summary lines, `:250`, `:253`). **W0-n1** `20261105T133349Z_GEO-36-C32-W0-n1.log` (18 s / 16.09 s, Status 0): 17 passed, 2 skipped (`:197`). K1: the degree knob (`2` and empty) is asserted refused **by the foreign guard** (match on its message plus the knob's name), which is what raises for it on this route — item 821's defect shape avoided, no gate weakened. K2 green. K7-3T (`:67–68`): rung Ø 0.704000000000 m, C20's rod (centre (0.105, 0, 0), `keep_cells` True), f 1.28e8 on `wf7` / `c4`, records 1128430 / (1508256, 7783028), (Bore) corner 0.214767316 m, gap 0.137232684 m; negative (I36b's `_run_drive` unpatched) 0.710 m / 6.4e7 / 1217025. K9-3T (`:114–116`): (K9-Q) 0.0 on both balls, (K9-E) 0.000e+00, census {3: 10260, 403: 108} = the mesh's, (A0m) clean, s2 M vs R 0.000e+00 and vs the 128 MHz circuit route 0.000e+00, vs the 64 MHz route 2.515e-01 (C21 W0's figure, asserted > 1e-10), heads present; control raised `R (gel twin) (P) P2 9.990e-04` with the JSON still written. **Disclosed:** K9-3T calls `_run_drive_3t` (the body, C21's patches included), not bare `_run_drive`; the 3 T arm R stub uses its own solver token because `_k9_stubs`' token is local; `_c21_patches`' patch 4 is inert on this route and I36b's "tuned on C4's cylinder" print is a stale label at 3 T (both in the docstring, neither edited); entry 48's Q 4.341688–4.460568 / η 3.642 are printed beside the stub's 4 / 2 (`:95`), never asserted. Not claimed: any 128 MHz tip Q or rod reading (entry 61).

### §9 item 825 — 2026-11-05 (§10 I67, chain I — `GEO-36` I67, the rod model's lever on one definition; chain step (solve))

**Executes** §10 I67 on I41's matched pair, under the ruling above. Family I67, attempt 1 of 4.

**Code.** A new module, `tests/validation/test_imp3_tip_q_hole_def.py`, with knob `FEM_EM_IMP3_I67=1`. No landed module is edited.
- **Imports:**
  - From `test_imp3_tip_q_entry48_def`: `_make_hook`, `_gel_seen`, `_gate_failures`, `_print_gates`, `_two_tag_resid`, `BALLS`, `ENTRY48_Q`, `ENTRY48_ETA_T1`.
  - From `test_imp3_rod_birdcage`: `_arm`, `_control_c1`, `PEC_WALL_MAX_E`, `PORT_IDS`, `ROD_CELL_RECORD_W1`.
  - I41's `_census`; `_build_rod`, `_build_kept_rod`, `KEPT_ROD_CELL_RECORD`, `ROD_CENTRE`; and `SOURCE` (ccw).
- **Order:**
  1. Build the hole mesh and run the hole arm with `_make_hook(store, "hole", …)`.
  2. `gc.collect()`.
  3. Build the kept mesh and run the gel twin.
  4. Set `arm["losses"]` from the store, as I57's `_run` does (`:265–266`). On the hole mesh, `dx(403)` integrates to 0.

**Gates (imported, unmoved).**
- **(A0):** the hole mesh has exactly 167 371 cells and the kept mesh exactly 182 371.
- **(A0m):** the owned tag-3 counts are equal, and the kept mesh has 15 000 tag-403 cells.
- **On each arm:**
  - (i) `RECIPROCITY_BAND`;
  - (ii) σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every column power ≤ 1;
  - (P) over all cells at `DISCRETE_IDENTITY_RTOL`. On the hole arm the two-tag residual is printed beside it, and is expected to be equal: a hole build has no lossy tag outside 3 and 1 (`OPS-80`).
- **(E):** on the hole arm, `wall_max` ≤ `PEC_WALL_MAX_E` for every drive, with wall dofs > 0.

**Controls (asserted).**
- **(C1):** `_control_c1({"rod": hole, "off": gel})` reads max|E| on the 402 dofs > `PEC_WALL_MAX_E`, with its (P) ≤ `DISCRETE_IDENTITY_RTOL`. `IMP-3` W1 measured 1.29e-4 on this build and drive (`20261002T170428Z_IMP-3-step1-W1.log`), a separation of ≈ 1.3e8. C1's "natural / without" prints carry predicted ranges written for the no-rod mesh; they are mislabelled here, so they are printed only.
- **(Rep):** the apex Q equals 4.090342 / 4.113695 within `DISCRETE_IDENTITY_RTOL` relative. The bound is named: the `.6f` print's half-ulp is ≈ 1.2e-7.

**Printed (record, no band).**
- Q = SAR_1g(hole) / SAR_1g(gel twin) on the four balls, beside I57's kept-cell 3.951922 / 3.973311 / 3.899322 / 3.921665 (`20261104T104820Z_GEO-36-I57-W1d.log:1999–2002`), with the hole / kept ratio per ball.
- η_rod(T1) = |E_hole| / |E_gel| at T1, beside 2.157532 (`:2003`), `IMP-3`'s (3.48, 3.56) and entry 48's figures.
- The ball masses on both arms, which should be equal.

**Windows.**
- **W0:** the coding tests, knob unset, `-n 2`, `timeout -k 30 120`.
- **W1:** `-n 2`, complex. Pre-check `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0. Bare and foreground, `timeout -k 30 600`, Bash timeout 660000 ms. The command is `20261104T104820Z_GEO-36-I57-W1d.log:12`'s, with the knob and the module swapped.
- **Price:** I41's two builds and three arms took 286.7 s at `-n 2` (`20261025T201613Z_IMP-3-I41-W1b.log:3780`). This run has two arms and a one-drive C1, so ≈ 230–300 s.

**Tier.** Heavy, complex, `-n 2`. **Class:** chain step (solve). **Cost:** ≈ 11.2 measured / 30 predicted min.

**Traps already paid for.** Items 779, 801, 804 and 810:
- Label the arms by content.
- `_arm_with_403(map_gel=False)` is vacuum.
- `_loss_power_w` is two-tag, so the all-cells form is the gate. If item 828 (`OPS-81`) has landed, the hole arm needs no acknowledgment (no lossy tag outside 1, 3), and `_gel_seen` already carries one.
- An exact float equality on a solve output does not reproduce at `-n 2`.
- Count verdicts per rank, by token.
- Run the leak guard on S digits.

**Scope.** One F-small pair, degree 1, 64 MHz, the ccw drive. **Not claimed:** convergence, which rod model is right, or anything clinical.

**Negative result.**
- (Rep) red: the pair is not I41's. Stop; it is a fixture finding.
- A gate red: BLOCKED.
- η is a record either way:
  - ≈ 3.5 means the lever is the rod model;
  - ≈ 2.2 means the lever was the mesh (entry 43's tip refinement).
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 I67 becomes DONE.
- The `GEO-36` row gains its I67 clause.
- The 11-04 weekly's "≈ 40 % lever" sentence is split into its model and mesh parts.

**Depends on:** nothing.

**Outcome 2026-11-05 (09:00 slot, item 825): DONE — every gate and both controls green; η is a record and reads ≈ 2.2, so the lever was the mesh.** New module `tests/validation/test_imp3_tip_q_hole_def.py` (no landed module edited; gates and bands imported unmoved), knob `FEM_EM_IMP3_I67` (unset skips, `1` runs, else `ValueError`). Item 828 (`OPS-81`) had not landed, so no acknowledgment was needed on either arm. **W0** `20261105T150251Z_GEO-36-I67-W0.log` (`-n 2`, 4 s, Status 0): 3 passed / 7 skipped per rank. **W1** `20261105T150305Z_GEO-36-I67-W1.log` (`-n 2`, complex, 268 s harness / 266.8 s pytest, Status 0; pgrep read 0 before it): **14 passed per rank** (two identical summary lines, `:3766`, `:3769`). (A0) `:1847`, `:3678`: hole 167 371, kept 182 371; (A0m) owned tag-3 51 693 = 51 693, kept tag-403 15 000. Hole arm (`:3711–3719`): (i) 1.201e-15, σ_max 0.999730355, (P) all-cells ≤ 4.150e-15, with the two-tag residual equal to it to ≤ 4.3e-15 and all − 3 − conductor ≤ 1.7e-18 W (`:3731–3734`, `OPS-80`'s premise confirmed). Gel twin (`:3721–3729`): (i) 2.290e-15, σ_max 0.999723588, (P) all-cells ≤ 8.300e-15 (two-tag 2.853e-07 on P1, printed). (E) `:3735`: 7 437 wall dofs, max|E| 0.000e+00 on every drive. (C1) `:3739`: natural wall max|E| 1.287e-04 (`IMP-3` W1 1.29e-4), (P) 2.158e-15. (Rep) `:3755`: apex Q 4.090342 / 4.113695, rel 3.210e-08 / 1.180e-07 against I41, and `_arm`'s own TIPS ratio prints the same. **Record** (`:3750–3756`): ball masses equal on both arms (0.980284 / 0.976690 / 0.981846 / 0.982903 g); Q = 4.119254 / 4.141242 / 4.090342 / 4.113695 (top / bottom / top_apex / bottom_apex), hole / kept-cell 1.042342 / 1.042265 / 1.048988 / 1.048966; η_rod(T1) 2.228368 against the kept rod's 2.157532 (+3.3 %), `IMP-3`'s (3.48, 3.56) and entry 48's 3.642324. **Reading:** on one definition and one mesh the rod-model lever is ≈ 3–5 % (Q) and ≈ 3 % (η). The 11-04 weekly's ≈ 40 % η lever is therefore mostly entry 43's tip-refined mesh, not the rod model. Not claimed: convergence, which rod model is right, anything clinical.

### Review note 2026-11-06 03:00 — `xl` entry 55 (I32) RUN, red on (S2-tip) alone

- **Window:** `20261106T080008Z_GEO-36-step10-I32.log`, `-n 8`, 1015 s, `memory.peak` 44.12 GiB, Status 1, orphans 0. The `log-pathologist` ruled it admissible and countable: 1 failed, 4 passed and 7 skipped on all eight ranks by token. The failing test is the body test only; the K-tests skip under the knob.
- **The red:** `rod (S2-tip)` 6.420e-03 against `S2_BAND` 1e-04, worst of the two tip balls (top 2.4589e-03, bottom 6.4201e-03; `:12512–12513`, `:12523`).
- **Green:** every EM gate ((A0), U₁, census 64, (D-rod), (M2) / (Z), (R1), (P) / NEG on both arms), (S0) with (A0m), (B), (E), (E′), and the global (S2) on both arms (1.2537e-05 / 1.2541e-05). The full table is in the ledger row.
- **What the review reads (a finding, for the weekly):** the gel twin, with no rod and no wall, prints tip misses of 3.1051e-03 / 1.3389e-03. They are not asserted there, since `_thermal_failures` gates the tip on the rod arm only. So the red is not about the rod; it is the same gate form that U127 failed. Filed as known-issues 2026-11-06. `OPS-83` (§9 item 859) measures whether the comparand ⟨q⟩_B should be the projected ⟨Πq⟩_B.
- **Disposition under the pre-registered rule:** "an energy or slope identity red with the EM gates green ⇒ a known-issues entry on the thermal route at clinical scale, the EM record stands". The θ figures are a record (tip ball-mean θ rod / gel twin ratio 1.0890 / 1.0898 against I41's 1.10). I68 waits on the weekly's reading.
- **Readout corrections for the weekly:** the tip θ is printed as a ball mean, not a centre value, and this module prints no Im Y_in (entry 48 carries it).

## §7 row moved verbatim 2026-11-07 (rotate_plan_archive.py chunks)

 A rod in the phantom: the implant-in-phantom generator (Phase 7 (7b)) | ✅ Audit PASS 2026-10-02 (CAD identities at `-n 2`); steps 3–8, I34 (`xl` entry 45), I36a–d and I32a (rod as σ_wire cells, ratified 11-02) are records in the chunk file — the kept rod on F-small, F-human and in F-clinical's ASTM box (1 217 025 cells, U₁ 1 606 744). **Entry 48 (I36, `xl` 11-02) is a record, read first branch by the 11-02 weekly with the pathologist's correction:** 5 passed / 12 skipped, tip 1 g Q 4.341688–4.460568, η_rod(T1) 3.642324, rod power share 3.740924e-05, arms M / R agreeing in S to ≈ 1e-7 (`20261102T080008Z_GEO-36-step9-I36.log:12050–12503`); box containment **not asserted** there, (A0m) / (P_RC) records not gates; entry 55 (I32, `xl` 11-06) **RUN red on (S2-tip) alone, EM gates green** (known-issues 2026-11-06; the 11-07 weekly reads it); I57 prints this tip-Q definition on `IMP-3`'s F-small rod mesh — **landed 2026-11-04 (§9 item 810) as a record:** on I41's kept-rod F-small build (182 371 cells) tip 1 g Q 3.90–3.97, η_rod(T1) 2.157532, every identity green with (P) over all cells (`20261104T104820Z_GEO-36-I57-W1d.log:1966–2012`); no band — read by the 11-04 weekly (≈ 10–12 % under the clinical box's Q). **I67 (§9 item 825) landed 2026-11-05 as a record:** entry 48's definition with the PEC-hole rod on I41's matched pair reads tip 1 g Q 4.09–4.14 (hole / kept-cell 1.042–1.049) and η_rod(T1) 2.228368 against the kept rod's 2.157532, every gate and control green, the apex Q reproducing I41's 4.090342 / 4.113695 to ≤ 1.2e-7 (`20261105T150305Z_GEO-36-I67-W1.log:3711–3757`). So on one definition and one mesh the rod-model lever is ≈ 3–5 %, and the 11-04 weekly's ≈ 40 % η lever was mostly entry 43's tip-refined mesh. **C32 (item 822) is I36-3T's body — landed 2026-11-05 as a stub proof:** `test_geo36_rod_f_clinical_3t_drive.py` runs I36b's `_run_drive` under C21's 3 T patches (K7-3T reaches 0.704 m / 128 MHz / 1 128 430 cells; K9-3T Q = 4, η = 2, s2 = the 128 MHz route, `20261105T133319Z_GEO-36-C32-W0.log`), entry 61's (i) MET. **Done-when:** the identities asserted in one logged window at `-n 2`, cells and time recorded ⇒ ✅. No solve. *History: `docs/planning/chunks/GEO-36.md`.* 

### §9 item 870 — 2026-11-08 (§10 I32b, chain I — the (S2-tip) re-pose on the consistent-mass projection; chain step (no solve))

*Enumerated by the 11-07 weekly (PROJECT_PLAN §10 I32b, `PROJECT_PLAN.md:4983–4999`; the ruling at `:4801–4817`: option (b), the comparand, never the band). Depends on nothing. Next: `xl` entry 63 (Friday 11-13), then U127′ (`WF-15.md` §9 item 884). Family I32b, attempt 1 of 4. The 11-07 18:00 rule "no slot edits `_thermal_failures`" is lifted for this item only, by the weekly's re-pose.*

**Executes** §10 I32b. Two premises made exact; no target moves.
1. **"The thermal helper it shares with `WF-15` step 14" does not exist on `main`.** Step 14's module (`tests/package/test_wf15_study_rod_heating.py`, on `attempt/WF-15-s14-20261105T212843Z` only) imports nothing from this module. It computes its tip comparand test-side (`_ball_qbar`). So I32b **writes the shared helper here**, as two module-level functions, and U127′ imports them by name. The comparand's code path is `OPS-83`'s `_projections` (`tests/validation/test_ops83_s2tip_form_on_stub.py:92–116`), which `OPS-84` already reuses (`test_ops84_s2tip_form_rod_pair.py:77`). It is imported, never copied.
2. **Other callers share `_thermal_arm` / `_thermal_failures`.** P61 calls them (`test_geo38_f_clinical_box_heating.py:158`, `:251`), and so does `xl` entry 59, QUEUED for 11-10 on that module. C31 calls `_thermal_failures` (`test_geo38_f_clinical_3t_box_heating.py:265`). P61's `_thermal_half` monkeypatches `i32._transient` and keeps the **last θ = 0.5 call** as its CN history (`:142–148`). So every new behaviour sits behind keywords that default to off, and P61, C31 and K10es reach byte-identical paths.

**Code.** Edit `tests/validation/test_geo36_rod_f_clinical_heating.py` only. No `src/` edit. `test_geo36_rod_f_clinical_drive.py` and `test_ops83_s2tip_form_on_stub.py` stay byte-identical.
- **The shared helper (new):**
  - `_projected_ball_comparands(sub, q, V1, centres, radius, tag)` → `{name: {"vol", "q", "pq"}}`. It calls `pq, _ = _projections(q, V1, tag)`, then reads each centre with `_ball_reader(sub, {"q": q, "pq": pq})`. **Import `_projections` inside the function body.** `test_ops83_s2tip_form_on_stub` imports this module at its top (`:32–55`), so a top-level import is a cycle.
  - `_tip_slope_rows(theta1_means, dt, comparands, rho_c)` → per ball: `R_Pi` = \|ρc⟨θ¹⟩/(Δt⟨Πq⟩) − 1\|, `R_DG0` = \|ρc⟨θ¹⟩/(Δt⟨q⟩) − 1\|, and `pq_q` = ⟨Πq⟩/⟨q⟩ − 1 (signed). It is plain arithmetic with no mesh, so U127′ can pass the `Study`'s own ball means.
- **`_thermal_arm(..., *, merged, tip_projection=False, f2182=None)`.** With `tip_projection=True`:
  - after `res1`, run `res10 = _transient(sub, ft, q, dt=DT_1_BOX / 10.0, n_steps=1, theta=1.0)`. It is BE, so P61's θ = 0.5 intercept never sees it;
  - add `"th10"` to the existing `read`;
  - call the helper on `V1` (`:298`) with `tag = "i32b_" + label.replace(" ", "_")`;
  - store `tip_rows` (`{"dt1": …, "dt1/10": …}`) and `tip_s2_pi` (R_Π at Δt₁/10, **the gate**). Keep `tip_s2` (R_DG0 at Δt₁) unchanged, printed as **the retired form**.
- **The `f2182` keyword** (a dict name → point). It adds one CN run, `_transient(sub, ft, q, dt=F2182_T_S / N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2)`, **run before the t\* CN**, so the t\* CN stays the last θ = 0.5 call. Read it with `_probe_theta(res900.temperature, pts, (msh, parent_ids))`.
- **`_thermal_failures`.** On the non-merged arm, gate `tip_s2_pi` when the key is present, otherwise the old `tip_s2`. Only `_run_heat` builds a non-merged arm. The merged-arm branch is unchanged.
- **`_run_heat`** passes `tip_projection=True` on both arms (the twin's tip stays ungated and printed, as in entry 55). It also passes `f2182`: per tip `n`, `apex + F2182_OFFSET_M · (apex − tip)/|apex − tip|`, with apex = `tips_all[f"{n}_apex"]` (`_probe_points`, `test_geo36_rod_f_clinical_drive.py:285–290`: ±(L/2 + `ROD_RADIUS`) on ẑ, so the clinical probes sit at the rod centre ± 17 mm). A tip with no apex key prints "no apex on this fixture" and is skipped.
- **New constants:** `F2182_T_S = 900.0` and `F2182_OFFSET_M = 1.0e-3`, each commented "quoted from memory, unconfirmed (W3)". Both are printed, never asserted.
- **`_stubs`:** wrap the stubbed `_probe_points` so it adds the planted keys `a_apex` / `b_apex` = tip + (0, 0, 1 mm). K10-R′ then runs the F2182 branch on the stub. `_k9_stubs` is not edited.
- **No knob is added.** Entry 63's command (`xl-pending.md` entry 63) stands byte for byte.

**Gates (asserted; bands imported, none set).**
- **K11** (new `test_i32b_k11_projected_tip_gate_on_checkerboard`, `@complex_only @_coding`). Use K10es's set-up verbatim (`:656–666`: the merged arm, `STUB_CENTRE`), with tips `{"a": fx[4]["a"], "interior": STUB_PROBES["interior"]}` (`OPS-83`'s two balls). Run two arms, labels `"K11 u"` and `"K11 r"`. Each label gives its own tag, because a reused PETSc prefix collides. On arm `r`, `this._heat_source` is patched inside `monkeypatch.context()` (as K10es's NEG is) to return `OPS-83`'s q_r = `_checkerboard(real(*a))[0]`. Import that inside the test body. Each arm calls `_thermal_arm(..., merged=True, tip_projection=True)`.
  - **(K11-Π)** On arm r, `tip_s2_pi` ≤ `S2_BAND` at both balls.
  - **(K11-rep)** On arm r, R_Π at Δt₁/10 equals `OPS-83`'s 6.8552e-07 (tip) and 8.5219e-07 (interior), within **5e-12** (`20261106T120232Z_OPS-83-W1.log:75–76`). That bound is half a unit of the `.4e` print's last digit, the 18:00 rule on a reproduction control. \|`pq_q`\| equals 1.573e-03 / 2.115e-03 within 5e-07 (the `.3e` half unit, same lines).
  - **(K11-u)** On arm u, R_Π = R_DG0 at Δt₁/10 within `PROJECTION_ID_TOL` (1e-10, imported from `test_ops83_s2tip_form_on_stub`), and R_Π equals 8.0980e-07 (tip, `:69`) within 5e-12.
  - **(E)** on both arms ≤ `E_BAND`; **global (S2)** on both arms ≤ `S2_BAND`. The checkerboard keeps ∫q: `OPS-83` read 2.0247e-05 on both sources (`:68`, `:74`). **(B)** is asserted on arm u only. On arm r it is printed: `OPS-83` never measured ∫q_r against G·P.
- **K11-f** (new, no mesh). `_thermal_failures` on two synthetic rod-arm dicts:
  - `tip_s2` 1.6e-03 with `tip_s2_pi` 7e-07 → no (S2-tip) failure;
  - `tip_s2_pi` 2e-04 → exactly one failure, named "(S2-tip)".
- **K10-R′.** `test_i32a_k10r_run_heat_end_to_end_on_stubs` stays green through the new path, with its asserts unedited.
- **K0, K10-q, K10es, K10-C, K10-p, K10-M** are unchanged and re-run. K10-S's tip 8.0400e-06 is re-asserted by `OPS-83`'s own reproduction control in Wreg.

**Negative control (asserted).** **(K11-c)** On arm r, the retired DG0 comparand on the same θ¹ misses: `tip_s2["a"]` > `S2_BAND`. `OPS-83` read 1.5793e-03 at Δt₁ (`:72`), ≈ 15.8× the band. The ceiling is the checkerboard's own form gap (\|⟨Πq⟩/⟨q⟩ − 1\| 1.573e-03, `:72`), so assert only `> S2_BAND` (rule (e)) and print the factor.

**Printed (≤ 7 digits, the sense label on every θ line).** On every arm the module runs, per tip ball:
- ⟨q⟩_B, ⟨Πq⟩_B, ⟨Πq⟩_B/⟨q⟩_B − 1;
- R_Π at Δt₁ and at Δt₁/10;
- R_DG0 at Δt₁, labelled "retired form";
- the F2182-shaped θ at each probe, labelled "(placement and time quoted from memory, unconfirmed; W3)".

**The new printed forms must pass the leak guard.** The daily checks them before queueing entry 63.

**Windows.** Before each window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Every window is a bare foreground harness call, Bash timeout 660000 ms. Use the durable capture of `20261104T030134Z_IMP-8-I55-W2-item803.log:12`, including the trailing `; exit $rc`. Put the complex prefix of `20261101T103929Z_GEO-36-I32a-W0b.log:12` in front (`source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1`). Knobs stay unset throughout; no gmsh runs.
- **W0** `GEO-36-I32b-W0`: `-n 1`, `timeout -k 30 120`, `tests/validation/test_geo36_rod_f_clinical_heating.py -v -s --tb=short`. I32a's W0 took 17 s (`20261101T103903Z_GEO-36-I32a-W0.log:111`).
- **W0b** `GEO-36-I32b-W0b`: the same at `-n 2`. I32a's W0b took 14 s (`20261101T103929Z_GEO-36-I32a-W0b.log:133`). K11 adds four BE solves, a CN pair and two projections on the 12³ stub; `OPS-83`'s whole module took 13 s (`20261106T120232Z_OPS-83-W1.log:94`). Expect ≈ 25–45 s.
- **Wreg** `GEO-36-I32b-Wreg`: `-n 2`, `timeout -k 30 180`, `tests/validation/test_ops83_s2tip_form_on_stub.py tests/validation/test_ops84_s2tip_form_rod_pair.py tests/validation/test_geo38_f_clinical_3t_box_heating.py tests/validation/test_geo38_f_clinical_box_heating.py`, knobs unset. Expected:
  - `OPS-83`: 5 passed per rank, reproducing K10-S and its table (`OPS-83-W1.log:87`);
  - `OPS-84` W0: 3 passed, 7 skipped (`20261106T180947Z_OPS-84-W0.log:71`);
  - P61 + C31: 14 passed, 2 skipped per rank (`20261105T181050Z_GEO-38-C31-W0r.log:333`, 33 s).
  
  Total ≈ 60–90 s. I58 (`_p37_fixture`) and I55 (`PERTURBED_SOURCE`) import unedited names only, so they are not re-run.

**Tier.** Smoke, complex, `-n 1` / `-n 2`. **Class:** chain step (no solve). **Cost:** ≈ 4.8 measured / 35 predicted min.

**Traps already paid for.**
- The ops83 import cycle (above): a function-level import only.
- P61's `_transient` intercept: the new keywords default off, and the F2182 CN is ordered before the t\* CN in any case.
- `_projections` builds its options prefix from the tag (`fem_em_ops83_proj_{tag}_`). Reuse collides (`OPS-84`'s trap), and a space in an arm label breaks it.
- `ArityMismatch` is a `BaseException` raised on rank 0 alone (item 849). `_projections` already conjugates (`ufl.inner(u, v)`). Never put a refusal inside `pytest.raises` around `fem.form`.
- `ufl.lt` on complex operands goes through `ufl.real`, and every ball form pins `quadrature_degree` (`_ball_reader` does both).
- The checkerboard is computed in numpy on cell midpoints (`_checkerboard`). Do the ball reads inside `monkeypatch.context()` (`OPS-83`'s trap).
- Reduce before asserting; count verdicts per rank by token; use `-s`.
- `.eval` bytes depend on the batch: the F2182 read is one batch per arm, printed, never compared.
- A leak-guard refusal moves the log unchanged to `logs/leak-held/` and the slot reports it. Never `--no-verify`.

**Scope.** The module's tip comparand is re-posed and stub-proved, and the shared helper is written. **Not claimed:**
- any clinical number (entry 63's);
- U127′;
- a band;
- the retirement of known-issues 2026-11-06 (it resolves only when entry 63 and U127′ both land green in this form);
- the F2182 standard's placement or time.

**Negative result.**
- (K11-Π) or (K11-rep) red: the projected read on this route is not `OPS-83`'s. Print every row, park on `attempt/GEO-36-I32b-*`, BLOCKED. Entry 63 stays PENDING.
- (K11-c) red (DG0 under the band): the patched source never reached `_thermal_arm`. BLOCKED.
- Wreg red on P61, C31 or `OPS-83`: the default-off path is not byte-identical. Revert and stop, BLOCKED, because entry 59 is queued on this module.
- A crash in non-gate code before any gate value exists: fix it and re-run once at the same `-t`, keeping the first log (the 10-31 narrow ruling).
- No `PREDICTED RED:` line.

**Status it can move.** §10 I32b DONE; `GEO-36` gains a record line (the re-posed tip gate, stub-proved). **The implementer who lands it marks `docs/testing/xl-pending.md` entry 63's status line `READY <date> (I32b landed, <commit>)` in the same commit, and runs nothing else** (step 6b rule 5): no queue file, no `xl` command. A commit cannot carry its own sha, so `<commit>` is written as "§9 item 870 landing commit". The next review substitutes the sha, checks the new printed forms against the leak guard and writes `xl-queue.d/20261113-GEO-36-step10-I32-prime.env`.

**How U127′ picks it up.** U127′ (item 884) runs on `main` after this lands and imports `_projected_ball_comparands` / `_tip_slope_rows` from this module by name. The park branch `attempt/WF-15-s14-20261105T212843Z` is never rebased or merged. Item 884 takes its code onto `main`, so the helper reaches U127′ by import, not through the branch.

**Depends on:** nothing.

**Outcome (04:30, 2026-11-08): BLOCKED, parked on `attempt/GEO-36-I32b-20261108T103328Z` (commit `01e58e60`).** W0 (`-n 1`, `20261108T103328Z_GEO-36-I32b-W0.log`, Status 1) went red on (K11-rep) only: R_Π against `OPS-83`'s literals off by 8.598e-10 (tip) / 5.407e-11 (interior) on the checkerboard and 2.641e-10 on the uniform tip, against 5e-12. W0b (`-n 2`, `20261108T103407Z_GEO-36-I32b-W0b.log`, Status 0) passed the same code at 4.402e-12 / 1.436e-12 / 3.954e-12. The executor then scoped the literal checks to `-n 2` (`OPS83_RECORD_RANKS = 2`, printed at other widths, tolerance unmoved). Its diagnosis is that the solve stores absolute T, so a residual near 1e-9 sits at a round-off floor that shifts with the partition (K10-S prints 8.0402e-06 at `-n 1` and 8.0400e-06 at `-n 2`). On that code `20261108T103526Z_GEO-36-I32b-W0r.log` (`-n 1`), `20261108T103552Z_GEO-36-I32b-W0b2.log` (`-n 2`: K11-Π 6.855244e-07 / 8.521914e-07, K11-c 15.79× band, K11-f, K10-R′ green) and `20261108T103615Z_GEO-36-I32b-Wreg.log` (18 passed / 9 skipped per rank, P61 + C31 + `OPS-83` + `OPS-84` green) were all Status 0. The slot did not land this: rule (e) and this brief's negative result both say a red on (K11-rep) parks. The review rules on the width scoping. The logs are on the branch.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 870 — `GEO-36` I32b — per `docs/planning/chunks/GEO-36.md` "§9 item 870". Stub proof only; no knob, no gmsh, no clinical window. `test_geo36_rod_f_clinical_drive.py` and `test_ops83_s2tip_form_on_stub.py` stay byte-identical. Mark `xl-pending.md` entry 63 READY in the landing commit and queue nothing. Every harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

### §9 item 872 — 2026-11-08 (§10 I68′, chain I — the lumped reading against the clinical FEM on one normalisation; lane 1 → chain step (no solve))

*Enumerated by the 11-07 weekly (`PROJECT_PLAN.md:5008–5018`; re-posed from I68 at `:4819–4823`). Depends on nothing. Family I68′, attempt 1 of 4.*

**Executes** §10 I68′. Two premises made exact; the target is kept and returned to the weekly.
1. **There is no saved JSON for entry 55.** `_run_heat` writes none (item 707's Outcome: "The body writes no JSON: entry 55's record is its log", `GEO-36.md:1273`). The only JSON under `/workspace/logs/geo36/` is entry 48's `i36-node-193.5996.json`, written by `_run_drive` (`test_geo36_rod_f_clinical_drive.py:607–620`). It is gitignored (`.gitignore:82`, `/logs/`) and present on this host's clone only, so a test that read it would be red on any fresh clone. **The reader therefore reads the tracked log `docs/testing/logs/20261106T080008Z_GEO-36-step10-I32.log`.** Its printed strings are the record that the ledger row and the 11-07 weekly read. (K) becomes: the log's values re-read equal, by string equality.
2. **No accepted power is in a tracked record for the superposed (1, +1j) drive.** Entry 55 prints P_vol and its gel / rod split (`:12497`) and the per-column (P) residuals (`:12492–12496`). Entry 48 prints the same split beside Y_in (`20261102T080008Z_GEO-36-step9-I36.log:12494`). Neither prints the superposition's source power. So "per accepted watt" is printed as **"not in a tracked record"**. The nearest tracked normalisation, per watt of P_vol (gel + rod), is printed beside it and labelled as such.

**Code.** New `tests/validation/test_geo36_i68_lumped_vs_clinical.py`: numpy / scipy, no mesh, both builds. No landed module is edited.
- **Imports, by name:**
  - from `test_imp8_entry48_lumped`: `_ball_average_axis`, `_line`, `_tips`, `_a_1g`, `N_BQ`, `PER_WATT`, `P_ROD_LOG`;
  - from `test_imp8_rod_line_source_closed_form`: `L`, `T_PRINT`;
  - from `test_mat21_box_fem`: `T_STAR` as `T_STAR_BOX`.
  
  **Never import a `test_*` name.** Pytest re-collects it in this module; I58 does that deliberately (`test_imp8_entry48_lumped.py:52–55`), and this module must not.
- **Locating the log.** `ENTRY55_LOG` = the basename above, resolved from `Path(__file__).resolve().parents[2] / "docs/testing/logs"`. Open the `.log`, or the `.log.gz` with `gzip.open` when the `.log` is absent. `scripts/maintenance/housekeeping.py` gzips cited logs at `COMPRESS_DAYS` = 7 (`:65`), so the `.log` is gone on ≈ 11-13. Never gunzip the record.
- **The literals.** Module constants, each with its line, transcribed from the log:
  - `:12497`: rod share `3.740924e-05`, P_vol `2.538288262e-03`, gel `2.538193307e-03`, rod `9.495544579e-08` W;
  - `:12499`: twin P_gel `2.538800824240e-03` W, G `3.938867478110e+02` 1/W, G·P `1.000000000000e+00` W;
  - `:12501–12502`: twin tip θ(t\*) `3.082473767e-02` / `3.090736898e-02` K, ball radius `6.2035` mm;
  - `:12509`: t\* `2821.500` s, dt `70.537500` s;
  - `:12510`: rod G·P `9.997607068690e-01` W;
  - `:12512–12513`: rod tip θ(t\*) `3.356701564e-02` / `3.368366738e-02` K;
  - `:12521`: ratio rod / gel twin `1.088964` / `1.089826`.
- **The parser.** `_entry55_record()` parses the `[GEO-36 I32]` lines by fixed regexes, and each pattern must match exactly once.

**Gates (asserted; no band set).**
- **(K)** Every parsed token string equals its literal (string equality), and every pattern matches exactly once.
- **(K-d)** Re-derived at the print's last digit:
  - θ_rod / θ_twin per tip within 1e-6 of the `.6f` ratios (`:12521`);
  - rod / P_vol within 5e-12 of 3.740924e-05 (half a `.6e` unit, `:12497`);
  - `f"{T_STAR_BOX:.3f}"` == `2821.500` (`:12509`; the FEM t\* is MAT-21's, not I58's 900 s).
- **(N)** I58's records reproduced from its own module. `_ball_average_axis(_line, ±L/2, _a_1g(), 2 * N_BQ) * PER_WATT` formatted `.7g` == `7.307753`; with `_tips`, == `12.93694` (`20261103T211525Z_IMP-8-item786-W1.log:63`, `:65`). `f"{_a_1g() * 1e3:.4f}"` == `6.2035`, entry 55's ball radius (`:12501`): one ball radius on both sides.

**Negative controls (asserted).**
- **(N-c)** I58's line kernel read at entry 55's t\* (`T_STAR_BOX`) instead of `T_PRINT` misses the 7.307753 record by more than (N)'s own bound, one unit of the 7th digit (rule (e): asserted at the gate's own bound). Its value is printed as the **time-matched** lumped figure.
- **(K-c)** A copy of line `:12512` with one digit changed, made in memory, fails (K).

**Printed (a record, ≤ 7 digits, ball radius 6.2035 mm stated on every line).** Per tip and arm:
- θ_FEM(t\*) per watt deposited in the gel. That is the print itself (G·P = 1 W on the twin; 0.9997607 W on the rod arm).
- θ_FEM / (G·P_rod), per watt deposited in the rod, with the conversion G·P_rod printed: watts in the rod per watt in the gel, G from `:12499`, P_rod from `:12497`.
- θ_FEM / (G·P_vol), per watt of P_vol.
- "per accepted watt: not in a tracked record".
- Beside them, I58 per rod watt at 900 s (ball average and surface point, `…item786-W1.log:63–66`) and at t\* (from (N-c)). The FEM ball mean sits beside the lumped ball average and the lumped surface point.
- The ratios FEM / lumped per rod watt, labelled with what differs:
  - FEM: SAR in the gel, the rod unheated (tag 403 outside the thermal domain);
  - lumped: P_rod in the rod's line and tips;
  - times 2821.5 s against 900 s;
  - a walled box against an infinite medium.

**Window.** Before the window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. **W1** `GEO-36-I68p-W1` is a bare foreground harness call, Bash timeout 660000 ms. It runs the real build, `-n 1`, `timeout -k 30 120`, `tests/validation/test_geo36_i68_lumped_vs_clinical.py tests/validation/test_imp8_entry48_lumped.py -v -s --tb=short`, in the durable capture of `20261104T030134Z_IMP-8-I55-W2-item803.log:12` (a real build) verbatim, including the trailing `; exit $rc`. I58's module rides along as its own regression. **Price:** I58's module ran in 8.33 s (`…item786-W1.log:87`). The parse adds well under 1 s and (N-c) adds two ball averages, so expect ≈ 15–25 s.

**Tier.** Standard (the weekly's label), real, `-n 1`. **Class:** lane 1 → chain step (no solve). **Cost:** ≈ 4.8 measured / 20 predicted min.

**Traps already paid for.**
- **The leak guard.** Print log-derived values in the log's own strings, never reformatted (I58's trap). The converted figures are new digit strings; a hook refusal moves the log unchanged to `logs/leak-held/`. Never `--no-verify`.
- The `.gz` rotation (above).
- Never read `logs/geo36/*.json`: it is gitignored and host-only.
- `quad` raises an `IntegrationWarning` inside `_phi` near the axis (I58's disclosure). It is not a red.
- Never open `docs/private/`.
- `-s`.

**Scope.** A reader with no band. **Not claimed:**
- a like-for-like comparison: the FEM has no rod source, and the lumped model has no box and no gel SAR;
- a FEM read at 900 s (that is entry 63's F2182 print and I76's);
- P_acc;
- any tissue or published comparand.

**Negative result.**
- (K) red: the tracked log differs from the record the weekly quoted. Report it, stop, BLOCKED: a finding. Never edit a literal to match the log.
- (N) red: I58's module no longer reproduces its record. Stop, BLOCKED.
- No `PREDICTED RED:` line.

**Status it can move.** §10 I68′ DONE. `GEO-36` gains a record line: what a lumped rod reading misses at clinical scale, on stated normalisations. The "per accepted watt" gap goes back to the weekly: it needs an entry-63 print, or none.

**Depends on:** nothing (entry 55's log and I58's module are on `main`).

**Outcome (06:00, 2026-11-08): DONE.** W1 `20261108T121322Z_GEO-36-I68p-W1.log` (real, `-n 1`, Status 0, 17 s elapsed, pytest 15.41 s): 10 passed, 0 failed, 0 skipped. Seven of those are I58's module riding along as its own regression. Gates: (K) 9 patterns, each matched once, 17 tokens, failures none; (K-c) a one-digit change to `:12512` is caught; (K-d) the rod / twin ratios re-derive to 4.539e-07 / 4.230e-07 (band 1e-6), rod / P_vol to 4.337e-12 (band 5e-12), and MAT-21's T_STAR prints `2821.500`; (N) the ball radius is `6.2035` mm on both sides and the line / tips values are `7.307753` / `12.93694` at both tips; (N-c) the line kernel read at t\* gives 9.367587 K/W, missing the record by 2.060 against a bound of 1e-6. Record (ball a 6.2035 mm), with G·P_rod = 3.740169e-05 W in the rod per W in the gel and G·P_vol = 0.9997981 W. θ_FEM per rod watt is 824.1536 / 826.3629 K/W (twin, top / bottom) and 897.4732 / 900.5921 K/W (rod). FEM / lumped per rod watt is 87.97928 / 88.21513 (twin) and 95.80623 / 96.13918 (rod) against the t\* ball line, and 112.778 to 123.2379 against the 900 s ball line. Per accepted watt: not in a tracked record, so that gap goes back to the weekly. The ratio is not like-for-like: the FEM heats by gel SAR with the rod unheated, while the lumped model puts P_rod in the rod's line and tips in an infinite medium.

**Review ruling (2026-11-08 10:30 daily review) on item 870's red — provisional, the 11-09 weekly ratifies or reverts.** (K11-rep)'s 5e-12 literal reproduction of `OPS-83`'s `.4e` prints is a **same-width** reproduction: `OPS-83`'s record was taken at `-n 2`, and §10 I32b's own step text is "smoke, `-n 2`", with its validation target quoted at three figures (R_Π 6.86e-07 / 8.52e-07). At `-n 1`, W0 reads 6.863798e-07 / 8.521359e-07 / 8.100641e-07 (`20261108T103328Z_GEO-36-I32b-W0.log:139–141`, `:135`). These equal the weekly's three-figure target, and they miss the four-figure literals by 8.598e-10 / 5.407e-11 / 2.641e-10. That spread between widths is five decades under `S2_BAND`. It is the partition-dependent round-off of an absolute-T solve, the same class as K10-S's tip at 8.0402e-06 (`-n 1`) against 8.0400e-06 (`-n 2`). At `-n 2` the same code reproduces the literals to ≤ 4.4e-12 (`20261108T103407Z_GEO-36-I32b-W0b.log`). **So the parked code's scoping stands:** the literals are asserted at `OPS83_RECORD_RANKS = 2` and printed at other widths, with the tolerance unmoved and every physics gate (K11-Π, (E), global (S2), (B) on arm u, K11-c) asserted at both widths. No band moves. The `-n 1` spread is a record. A `log-pathologist` read the five branch logs before this ruling (journal `2026-11-08-daily-1030.md` §4) and returned CONFIRMED, with three things not established:
- **Which code ran.** All five logs carry `Commit: 5054e5fd`, and the code landed in `01e58e60` after them, so the code is identified by its printed strings only. Item 890's re-run on committed code closes this.
- **The mechanism.** It is not shown: relative to R_Π, the r-a spread (1.25e-3) is far larger than K10-S's (2.3e-5).
- **Margins.** (B) on arm r is printed, not asserted, and the |pq_q| checks pass at 0.93× / 0.97× of their tolerance at both widths.

Item 890 lands it, and adds the cheap decider the pathologist named: a printed `-n 3` window.

### §9 item 890 — 2026-11-08 (§10 I32b, chain I — item 870 landed under the 10:30 ruling; chain step (no solve))

*Written by the 2026-11-08 10:30 daily review. It is attempt 2 of family I32b, the landing of attempt 1's parked commit under the ruling above, with no new code. Depends on nothing. Next: U127′ (item 884), then `xl` entry 63.*

**Executes** §10 I32b by landing `01e58e60` (on `attempt/GEO-36-I32b-20261108T103328Z`) on `main`, then re-running its windows on `main`.

**Steps.**
1. **Precondition (asserted by `git`, before any window).** `git diff 5054e5fd HEAD -- tests/validation/test_geo36_rod_f_clinical_heating.py tests/validation/test_geo36_rod_f_clinical_drive.py tests/validation/test_ops83_s2tip_form_on_stub.py tests/validation/test_ops84_s2tip_form_rod_pair.py tests/validation/test_geo38_f_clinical_box_heating.py tests/validation/test_geo38_f_clinical_3t_box_heating.py` must be empty. If it is not, stop BLOCKED: `main` has moved under the branch.
2. **Take the code by path:** `git checkout attempt/GEO-36-I32b-20261108T103328Z -- tests/validation/test_geo36_rod_f_clinical_heating.py`. Then take the five attempt-1 logs by path. They are new files, the records of attempt 1. Their headers read `Commit: 5054e5fd` because the code was uncommitted when they ran:
   - `…103328Z_GEO-36-I32b-W0.log`
   - `…103407Z_…-W0b.log`
   - `…103526Z_…-W0r.log`
   - `…103552Z_…-W0b2.log`
   - `…103615Z_…-Wreg.log`

   Append their five test-results rows **by hand** from `git show 01e58e60 -- docs/testing/test-results.md`. Never check out a record file by path. `git diff --stat` on `docs/testing/test-results.md` must show insertions only.
3. **Add one docstring line** next to `OPS83_RECORD_RANKS` (`:940` on the branch): "Width-scoped by the 2026-11-08 10:30 review ruling (`GEO-36.md`, provisional; the 11-09 weekly ratifies)". Make no other edit.
4. **Re-run on `main`.** These are item 870's W0b / W0 / Wreg commands, byte for byte except the labels: `GEO-36-I32b-L-W0b` (`-n 2`), `GEO-36-I32b-L-W0` (`-n 1`) and `GEO-36-I32b-L-Wreg`. A fourth window, **`GEO-36-I32b-L-W0n3`**, is L-W0b's command at `-n 3`, with the module path narrowed to `::test_i32b_k11_projected_tip_gate_on_checkerboard`. It is the pathologist's decider, and nothing in it is asserted beyond the module's own gates. Expected: K11-Π green, and the literals printed with a miss in [1e-10, 1e-9]. A miss outside that range is recorded for the weekly, not a stop. Each is a bare foreground harness call with Bash timeout 660000 ms and `pgrep -c python3` reading 0 first.

**Gates (asserted, imported, none set).** The module's own gates are unchanged:
- (K11-Π) R_Π ≤ `S2_BAND` at both balls;
- (K11-rep) at `-n 2` within 5e-12 of 6.8552e-07 / 8.5219e-07, and (K11-u) of 8.0980e-07;
- (E), global (S2), and (B) on arm u;
- K11-f and K10-R′.

L-W0b must print K11-Π equal to the branch's W0b2 line by `.6e` string (6.855244e-07 / 8.521914e-07). L-Wreg must read 18 passed / 9 skipped per rank, by token.

**Negative control (asserted, backed).** (K11-c): the retired DG0 tip on q_r is greater than `S2_BAND`. `OPS-83` and W0b2 read 1.5793e-03, 15.79× the band.

**Tier.** Smoke, complex, `-n 2` / `-n 1`. The windows took 20–38 s each on the branch. **Class:** chain step (no solve). **Cost:** ≈ 5.0 measured / 20 predicted min.

**Traps already paid for.**
- A path checkout of a record file drops rows that `main` appended later (`8d4cf58`).
- The complex prefix and durable capture are item 870's.
- Count verdicts per rank by token, and use `-s`.
- A leak-guard refusal moves the log unchanged to `logs/leak-held/`. Never `--no-verify`.

**Scope.** The landing and a re-run only. **Not claimed:**
- any clinical number (entry 63's);
- U127′;
- the retirement of known-issues 2026-11-06;
- ratification of the width ruling (the weekly's).

**Negative result.** Stop BLOCKED and revert to `main` in any of these cases:
- the precondition fails;
- L-W0b is red;
- K11-Π's `.6e` differs from W0b2's;
- L-Wreg is red.

Then print the rows and leave entry 63 PENDING. No `PREDICTED RED:` line.

**Status it can move.**
- §10 I32b DONE, and `GEO-36` gains its record line: the branch's `PROJECT_PLAN.md` row edit, applied by hand.
- §9 item 870 → DONE via 890.
- **`docs/testing/xl-pending.md` entry 63's status line becomes `READY 2026-11-08 (I32b landed, §9 item 890 landing commit)`** in the same commit (step 6b rule 5). Queue nothing.
- Item 884 is unblocked.

**Depends on:** nothing.

**Outcome (2026-11-08 12:00 slot, take-next second item, implementer, not delegated; attempt 2 of 4): DONE.**
- **Precondition.** `git diff 5054e5fd HEAD --` on the six named modules was empty. The module was taken from `attempt/GEO-36-I32b-20261108T103328Z` by path, and so were the five attempt-1 logs, which are new files. Their five test-results rows were appended by hand from `git show 01e58e60`, insertions only. The only edit is the docstring line beside `OPS83_RECORD_RANKS`. `pgrep -c python3` read 0 before the first window.
- **L-W0b** `20261108T181211Z_GEO-36-I32b-L-W0b.log` (`-n 2`, Commit `505d8ef1`, which closes the pathologist's "which code ran"): 9 passed / 1 skipped on each rank, 22 s, `[capture] rc=0`.
  - K11-Π r: tip 6.855244e-07, interior 8.521914e-07. These equal W0b2's `.6e` strings (`:159–160` against `…103552Z_…-W0b2.log:159–160`).
  - (K11-rep) |diff| 4.402e-12 / 1.436e-12, u 3.954e-12, asserted (`:161`).
  - (K11-c) 1.5793e-03 = 15.79 × `S2_BAND` (`:162`).
- **L-W0** `20261108T181245Z_GEO-36-I32b-L-W0.log` (`-n 1`): 9 passed / 1 skipped, 19 s, rc 0. R_Π r 6.863798e-07 / 8.521359e-07, the literals printed with |diff| 8.598e-10 / 5.407e-11 / u 2.641e-10, equal to attempt 1's W0. K11-c is green (`:139–142`).
- **L-W0n3** `20261108T181312Z_GEO-36-I32b-L-W0n3.log` (`-n 3`, the K11 node only, the decider): 1 passed on each of 3 ranks, 11 s, rc 0. R_Π r 6.859740e-07 / 8.534698e-07. Literal misses: a 4.540e-10, interior **1.280e-09**, u 7.460e-10 (`:83–85`). K11-Π and K11-c are green.
  - Two of the three lie inside the brief's [1e-10, 1e-9]. The interior miss lies just above it, which the brief makes a record for the weekly, not a stop.
  - Across the three widths, the interior R_Π moves by 1.3e-09 against `S2_BAND` 1e-04. That is five decades under the band and consistent with the partition-dependent round-off reading. The mechanism is still not shown.
- **L-Wreg** `20261108T181335Z_GEO-36-I32b-L-Wreg.log` (`-n 2`): 18 passed / 9 skipped on each rank, by token, 37 s, rc 0 (`:377–381`).
- **Status moved:** §10 I32b DONE; §9 870 DONE via 890; `xl-pending.md` entry 63 READY in this commit, with nothing queued. Item 884 is unblocked. The width ruling's ratification stays the 11-09 weekly's. The attempt branch is kept for the review to dispose of.

### §9 item 893 — 2026-11-08 (§10 I79, chain I — daily, 2026-11-08, provisional, rule (6); chain step (no solve))

*Written by the 2026-11-08 18:00 daily review under §9 rule (6): toward I76's stated readout ("entry 63's probe prints beside I75's bracket arms and I78's closed form, per watt deposited in the rod; JSON written. A record", `PROJECT_PLAN.md:5068–5071`), the weekly's target, not moved. It prepares I76's table from the landed records, with entry 63's column blank, as C36 did for `ANS-11`'s entry-62 rows (item 871). Family I79, attempt 1 of 4. **Executor:** the implementer. Depends on nothing. I76 stays serial on entry 63.*

**Why now.** I76 is a reader over a log that does not exist until Friday 11-13. Every piece it needs except entry 63's own lines is on `main`: the print form entry 63 will write (I32b's stub proof), the bracket arms (I75) and the closed form (I78). If the reader is built and proved on those now, I76 shrinks to setting one basename and transcribing literals, and a print-form mismatch is found before the window, not after it.

**One geometric fact the table must carry (unlogged review arithmetic, asserted below as (P)).** The two probes are each "1 mm beyond the tip" of their own fixture, and they are not the same point relative to the source:
- I75's bracket rod is a flat-ended source (r ≤ `ROD_RADIUS`, |z| ≤ L/2), so its probe sits at z = ±16 mm (`test_imp8_rod_probe_f2182_fem.py:74`).
- The clinical rod is a capsule. `_probe_points` puts the apex at ±(`ROD_LENGTH`/2 + `ROD_RADIUS`) (`test_geo36_rod_f_clinical_drive.py:288–289`; `ROD_LENGTH` 0.030, `ROD_RADIUS` 1.0e-3, `tests/mesh/test_birdcage_implant_rod.py:93–94`), and `_run_heat` puts the probe `F2182_OFFSET_M` = 1 mm beyond it (`test_geo36_rod_f_clinical_heating.py:293`, `:582–584`). So it sits at ±17 mm from the rod centre: **2 mm beyond the line's end, which is I78's `Z_C`, not its `Z_P`**.
- On I78's line source that 1 mm moves θ(900 s) by 1 − 6.810428/9.648662 = 0.2941583 (`20261108T133157Z_IMP-8-item874-W1.log:54`). That is a 29 % confound if the rows are set side by side unlabelled. The table therefore carries the clinical column beside I78 at **both** placements, and labels I75's arms "16 mm, 1 mm beyond a flat end".

**Code.** New `tests/validation/test_geo36_i76_f2182_table.py`: numpy / scipy, no mesh, both builds. No landed module is edited.
- **Imports, by name (never a `test_*` name):**
  - from I68′ (`test_geo36_i68_lumped_vs_clinical`): `PATTERNS`, `_H`, `_S`, `_p`. **Not** `_parse` / `_log_text`: `_parse` checks tokens against *its own* module globals (`:118–133`) and `_log_text` reads one fixed basename (`:110–115`), so neither generalises. Use only the pattern strings `PATTERNS["share"][0]`, `PATTERNS["twin_g"][0]` and `PATTERNS["rod_g"][0]`;
  - from I78 (`test_imp8_probe_f2182_closed_form`): `SHAPES`, `_reads`, `Z_P`, `Z_C`;
  - from I75 (`test_imp8_rod_probe_f2182_fem`): `Z_P` as `Z_P_BRACKET_MM`, `ARMS`;
  - `ROD_LENGTH`, `ROD_RADIUS` from `tests.mesh.test_birdcage_implant_rod`; `MM` from `test_imp8_rod_line_source_closed_form`.
- **A generic log reader** `_read_log(basename)`: `docs/testing/logs/<basename>`, or its `.log.gz` through `gzip.open` when the `.log` is absent. `housekeeping.py` gzips cited logs at 7 days, so the stub and I75 logs rotate on ≈ 11-15. Never gunzip a record.
- **The F2182 pattern (new, one string):** `_H + r"\[{arm}\] F2182-shaped probe {tip} \(([^)]*)\) \((\S+) mm beyond the apex\), " + _S + r": theta\((\S+) s\) (\S+) K \(valid (\S+);"`. It is the print at `test_geo36_rod_f_clinical_heating.py:336–339`.
- **`_clinical_rows(text, tips)`** → per arm (`gel twin`, `rod`) and tip, (point string, offset, t, θ, valid), with each pattern required to match **exactly once**. It also returns the twin G (`twin_g`), the rod G·P (`rod_g`) and the `share` tokens. Per watt deposited in the rod is θ / (G·P_rod), with G the twin's (as I68′ converts, `test_geo36_i68_lumped_vs_clinical.py:195–205`).
- **`ENTRY63_LOG = None` and `ENTRY63_TIPS = ("top", "bottom")`** are module constants. While the basename is `None`, the clinical column prints "blank: entry 63 not run (`xl` Friday 11-13; I76 sets `ENTRY63_LOG`)". **No glob:** a reader whose output changes when a file appears is not a record. I76's whole edit is to set the basename, transcribe its literals and write the JSON.
- **I75's lines (new patterns):**
  - `^\[imp8-i75\] probe \(0, ([+-]16)\) mm arm {arm}: theta\(900 s\) (\S+) K/W; .*D=0\.5 superposed theta\(900 s\) (\S+) K/W$`, exactly twice per log;
  - `^\[imp8-i75\] probe arm {arm}: mean (\S+) K/W, .*D=0\.5 mean (\S+)$`, exactly once.
  - The logs are W1 `20261108T194323Z_IMP-8-item881-W1.log` (ti), **W2b** `20261108T195853Z_IMP-8-item881-W2b.log` (gel) and W3 `20261108T200636Z_IMP-8-item881-W3.log` (hole), at `:120–122` of each.

**The literals (module constants, each with its line).**
- **Stub** `20261108T181211Z_GEO-36-I32b-L-W0b.log`:
  - `:103` share `1.041667e-02`, P_vol `4.500000000e-03`, gel `4.453125000e-03`, rod `4.687500000e-05`;
  - `:106` twin P_gel `1.125000000000e-03`, G `8.888888888889e+02`, G·P `1.000000000000e+00`;
  - `:117` rod G·P `3.958333333333e+00`;
  - `:108–109` gel twin, a `(0.03, -0.03, 0.002)` and b `(-0.03, 0.03, 0.002)`, offset `1.0`, t `900`, θ `+2.121079410e-01` both, valid `True`;
  - `:119–120` rod, the same points, θ `+8.484317694e-01` both, valid `True`.
- **I75** (`:120–122` of each log):
  - ti: `1.064117e+01` / `1.049960e+01`, D `2.109553e+00` / `2.099602e+00`, means `1.057039e+01` / `2.104578e+00`;
  - gel: `9.220827e+00` / `9.102961e+00`, `1.979610e+00` / `1.970586e+00`, means `9.161894e+00` / `1.975098e+00`;
  - hole: `8.815610e+00` / `8.750725e+00`, `1.895292e+00` / `1.888582e+00`, means `8.783168e+00` / `1.891937e+00`.
- **I78** (`20261108T133157Z_IMP-8-item874-W1.log:58–59`, `:61–62`), as `.7g` strings: line at Z_P `9.648662` / `1.978384` / `15.18156`, at Z_C `6.810428` / `1.822237` / `12.25753`; twin at Z_P `63.14` / `2.498341` / `68.45374`, at Z_C `29.97039` / `2.451249` / `35.22961`.

**Gates (asserted; no band set).**
- **(K-s)** On the stub log, the F2182 pattern matches exactly once for each of (`gel twin`, `rod`) × (`a`, `b`), and `share` / `twin_g` / `rod_g` once each. Every token string equals its literal (string equality). This is the proof that I76's parser reads the form entry 63 will print.
- **(K-b)** On I75's three logs, every token string equals its literal and every pattern has its stated count.
- **(K-d)** Re-derived at the print's last digit:
  - each I75 mean against (a + b)/2 of its own two prints, within **one** unit of the mean's last printed digit (two half-unit roundings meet). Unlogged review arithmetic: the worst is ti θ at 5.0e-06 against 1e-05; gel θ and both gel / hole duty means are exact.
  - the stub's P_rod / P_vol against `1.041667e-02`, within half a `.6e` unit (5e-09); predicted 3.33e-09.
- **(N)** I78 reproduced from its own module: `_reads(kern, steady, z)` for both `SHAPES` at `Z_P` and `Z_C`, formatted `.7g`, equals the twelve literals.
- **(P)** The placement identity: |`ROD_LENGTH`/2 + `ROD_RADIUS` + 1.0·`MM` − `Z_C`| ≤ 1e-15 m, with the 1.0 parsed from the stub's offset token. Also |`Z_P_BRACKET_MM`·`MM` − `Z_P`| ≤ 1e-15 m. The clinical column is then mapped to I78's Z_C row by these constants, never by a literal 17.

**Negative controls (asserted).**
- **(K-c)** A copy of the stub's line `:119` with one digit of θ changed, made in memory, fails (K-s), as I68′'s (K-c) does.
- **(N-c)** I78's line read at `Z_P`, put in the clinical-placement row, misses the `Z_C` record `6.810428` by more than (N)'s bound (one unit of the 7th digit, 1e-06). Its separation is 2.838 K/W, unbounded above (a difference, not a ratio). A reader that keyed the clinical column on "1 mm beyond the tip" would show it.

**Printed (≤ 7 digits; a record, no band; the source and normalisation stated on each row).** The I76 table:
- rows ti / gel / hole (I75: 16 mm, a flat end, 1 W in the rod, self-heated), θ(900 s) mean and D = 0.5 mean;
- line / twin at 16 mm and at 17 mm (I78: infinite gel, 1 W per rod);
- the clinical rod / gel twin at 17 mm, per watt deposited in the rod, **blank**;
- beside them, the stub dry run of that same column path: G·P_rod 4.166667e-02 W, so θ per rod watt is rod 20.36236 / twin 5.090591 K/W (unlogged review arithmetic, predicted), labelled "**stub, form proof, not a record**".

Each clinical cell carries what differs from the bracket: gel SAR with the rod unheated (tag 403 outside the thermal domain), against 1 W in the rod body; a capsule against a flat end; the ASTM box against I55's gel cube. The FEM-over-closed-form ratios are printed for the bracket rows only.

**Window.** Before it, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. **W1** `GEO-36-I79-W1`, smoke, real, `-n 2`. Make one bare foreground harness call, Bash timeout 660000 ms; I78's module rides along as its own regression:

`scripts/testing/run_and_log.sh GEO-36-I79-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/validation/test_geo36_i76_f2182_table.py tests/validation/test_imp8_probe_f2182_closed_form.py -v -s --tb=short'"`

**Price.** I68′'s reader with I58's module took 17 s at `-n 1` (`20261108T121322Z_GEO-36-I68p-W1.log`, pytest 15.41 s). I78's module took 5 s (`20261108T133157Z_IMP-8-item874-W1.log:74`). Four log parses and twelve `_reads` calls add a few seconds, so expect ≈ 15–30 s.

**Tier.** Smoke, real, `-n 2`. **Class:** chain step (no solve). **Cost:** ≈ 5.0 measured / 15 predicted min.

**Traps already paid for.**
- **Never read `logs/geo36/*.json`.** It is gitignored and host-only (item 872's premise 1).
- **Never read the crashed W2** (`20261108T195128Z_IMP-8-item881-W2.log`): the gel arm is W2b.
- The `.gz` rotation: every cited log, not only entry 55's.
- **Never import a `test_*` name**; pytest re-collects it in this module.
- `_say` lines appear once per log (rank 0), but pytest's own lines appear once per rank. Anchor every pattern on `^\[GEO-36 I32\] ` or `^\[imp8-i75\] `, with `re.MULTILINE`.
- **The leak guard.** Print parsed values as the log's own strings. The converted figures are new digit strings, so a hook refusal moves the log unchanged to `logs/leak-held/` and the slot reports it. Never `--no-verify`.
- `-s`, or the prints are lost.
- Never open `docs/private/` or any `aed_results/`.

**Scope.** A reader and its stub proof. **Not claimed:**
- any clinical number (entry 63's) and I76 itself (its basename, literals and JSON);
- a like-for-like comparison: the bracket heats the rod's body with 1 W, while the clinical FEM heats the gel by SAR with the rod unheated, so "per watt deposited in the rod" is a normalisation, not a matched source;
- a band; the F2182 placement or time (W3);
- the placement ladder (I80 / I81, items 894 / 895). The reader imports neither, and I76 may add them if they land.

**Negative result.**
- **(K-s) red:** the reader does not read the form entry 63 will print. Report the mismatching line, stop, BLOCKED. It is a finding for the next review **before** Friday's window. The module's print and the queued command stay untouched, because entry 63 is QUEUED on them.
- **(K-b) or (N) red:** a landed record no longer reads as `IMP-8.md` quotes it. BLOCKED. Never edit a literal to match.
- **(P) red:** the probe geometry is not what this brief derived. Report both placements and stop.
- No `PREDICTED RED:` line.

**Outcome (21:00 CST slot 2026-11-08, attempt 1 of 4) — DONE, green.** Executor: `implementer`, foreground; `pgrep -c python3` read 0. One new module, `tests/validation/test_geo36_i76_f2182_table.py`; no landed module or `src/` edited; every literal checked against its cited line before it was written. **W1** `20261109T032235Z_GEO-36-I79-W1.log` (smoke, real, `-n 2`, 3 s): 13 passed per rank (26 `PASSED` tokens), I78's module riding along. (K-s) the stub's four F2182 rows and share / twin_g / rod_g once each, string-equal (`:64`); (K-b) I75's ti / gel (W2b) / hole lines string-equal at their counts (`:68–70`); (K-d) ti θ 5.00e-06 ≤ 1e-05, D 5.00e-07 ≤ 1e-06, gel 0 / 0, hole 5.00e-07 / 0, stub P_rod / P_vol 3.33e-09 ≤ 5e-09 (`:72–78`); (N) the twelve I78 `.7g` strings equal (`:80–83`); (P) 0.0 m on both identities (`:87`). Controls: (K-c) a changed θ digit on the stub's `:119` fails (K-s) (`:66`); (N-c) Z_P read in the Z_C row misses by 2.838 K/W (`:85`). The I76 table is printed at `:89–102`, the clinical column blank; the stub's per-rod-watt dry run reads rod 20.36236 / twin 5.090591 K/W, as predicted, labelled "stub, form proof, not a record". **Disclosed:** `STUB_TIPS` and `STUB_ROD_LINE` constants added; (K-c) asserts the line's prefix before mutating it.

**Status it can move.** §10 I79 DONE. `GEO-36` gains a record line: the I76 table reader, stub-proved on I32b's print form, with entry 63's column blank. I76 reduces to setting `ENTRY63_LOG` and its literals, then writing the JSON.

**Depends on:** nothing (the stub log, I75's three logs, I78's module and I68′'s module are on `main`).

### §9 item 932 — 2026-11-10 (§10 I85, chain I — rule (4), the 11-09 weekly's step; chain step (no solve))

*Enumerated by the 11-09 weekly (`PROJECT_PLAN.md:4947–4958`; the chain I header `:4932–4934`; entry 67 `:4959–4963` and `docs/testing/xl-pending.md` entry 67). The 2026-11-10 03:00 daily review wrote this brief as execution detail. The target, bands and done-when are the weekly's; none is set or moved. Family I85, attempt 1 of 4. **Executor:** the implementer. **Depends on nothing:** this is a **new module that only imports** the queued windows' modules. Entry 61 (11-12) runs C32's module, entry 63 (11-13) runs the heating module, and entry 60 (11-11) imports the heating module and C21's. The 11-09 ruling "a queued window's module is not touched" (`PROJECT_PLAN.md:4824–4827`) is met by construction, and checked below.*

**What it is.** Entry 67 (I32-3T) is entry 63's twin at 128 MHz. Its body is I32b's `_run_heat` run under C21's 3 T patch set, exactly as C32 runs I36b's `_run_drive` under it (`test_geo36_rod_f_clinical_3t_drive.py:120–144`).

**Two premises made exact (no target moves).**
1. **The thermal half needs no 3 T override.**
   - `_run_heat` takes its frequency from `_drive_env`'s ω (`test_geo36_rod_f_clinical_heating.py:522–523`). C4's fixture builds that ω from the `wf7` / `c4` globals, which C21's patch 3 rebinds (`test_geo38_f_clinical_3t_box_loaded.py:130–131`).
   - It takes the rod centre from `i35._rod_centre()` (`heating.py:538`). That reads `s2._f_clinical_kwargs` at call time (`test_geo36_rod_f_clinical.py:88–91`, through `s2._box_kwargs`), and patch 1 rebinds it (`3t_box_loaded.py:125`).
   - C32 printed C20's rod centre (0.105, 0.0, 0.0) under the patch (`20261105T133319Z_GEO-36-C32-W0.log:67`). It is the same ASTM box at both field strengths, so the F2182 probes sit at the same ±17 mm (I79's (P)).
2. **On the stub, the thermal half sees the same two fields at both frequencies.** The rod arm gets E_u and the twin 0.5 E_u (`_k9_stubs_3t`'s `_drive`, `3t_drive.py:315–325`). So I32b's stub prints must reproduce by string under the 3 T route. Only the (M2) / (Z) line and the s2 change, and C32's stub block anchors those.

**Code.** New module `tests/validation/test_geo36_rod_f_clinical_3t_heating.py`. **Byte-identical:** `test_geo36_rod_f_clinical_heating.py`, `test_geo36_rod_f_clinical_3t_drive.py`, `test_geo36_rod_f_clinical_drive.py`, `test_geo38_f_clinical_3t_box_loaded.py` (C21), `test_geo38_f_clinical_3t_rod_box.py` (C20), and `src/`.
- **Knob.** The implementer names it. Recommended: `FEM_EM_GEO36_S12_HEAT`, because entry 67 is `GEO-36` step 12 and C32's knob is `S11_DRIVE`.
  - Unset ⇒ the gated test skips and the coding tests run. `1` ⇒ `_run_heat_3t`. Anything else ⇒ `ValueError`.
  - The node is `FEM_EM_GEO36_S8_DRIVE_C_PF`, read through `i36._node_knob(mode)`. It is required under `1`; if missing, `ValueError` before gmsh.
  - Degree 1 only.
- **Foreign guard: 19 names.**
  - C32's `_foreign_envs()`: I36b's 7 plus `FOREIGN_ENVS_C32_EXTRA`'s 10, which already hold `FEM_EM_GEO36_S8_DRIVE` and `FEM_EM_GEO36_S8_DRIVE_DEGREE` (`3t_drive.py:52–63`, `:83–84`).
  - Plus `heat.HEAT_ENV` (`FEM_EM_GEO36_S8_HEAT`) and `c32.S11_DRIVE_ENV` (`FEM_EM_GEO36_S11_DRIVE`).
  - The node knob is not in it. The guard's message head is this module's own, and K1-12 matches it.
- **`_run_heat_3t(comm, monkeypatch, c_pf)`**, three lines:
  1. `c21._c21_patches(monkeypatch)`;
  2. one `_say` line in C32's form (`3t_drive.py:129–132`: rung Ø, f, records, node, "patch 4 inert on this route");
  3. `return heat._run_heat(comm, monkeypatch, c_pf)`.
- **Gated test** `test_geo36_i32_3t_rod_heating` (`@complex_only`): the knobs, a skip when unset, then `_run_heat_3t(MPI.COMM_WORLD, monkeypatch, c_pf)`.
- **The 3 T stub, by import and never copied.** Inside `monkeypatch.context()`:
  - set `this = sys.modules[heat._heat_columns.__module__]`;
  - `m.setattr(this, "_k9_stubs", wrapper)`. The wrapper calls `c32._k9_stubs_3t(m_, fx_, jd_, p_src_r2=p)`, keeps `(lines, z_sym, z_rc)` in a closure dict, and returns the first two;
  - then call `heat._stubs(m, fx, root)`. It rebinds `_route_drive`, `_drive_env`, `_drive` and `_probe_points` from the drive module, which are now C32's 3 T versions (`heating.py:643–647`);
  - a recorder wraps `this._heat_columns` and keeps each call's `s2`.
- **Print-only.** Entry 55's 1.5 T tip θ ratio rod / gel twin, `1.088964` / `1.089826` (`20261106T080008Z_GEO-36-step10-I32.log:12521`), printed as the log's strings. Label: "1.5 T, entry 55, before I32b's gate".

**Gates (asserted; bands imported, none set).** Literal string comparisons are asserted at `-n 2` only and printed at other widths (`OPS83_RECORD_RANKS`, `heating.py:941`; the width ruling, `PROJECT_PLAN.md:3966`, ratified at `:4828`).
- **(K-q), the queued-module check, before the commit:** `git diff --stat HEAD` lists none of `tests/validation/test_geo36_rod_f_clinical_heating.py`, `tests/validation/test_geo36_rod_f_clinical_3t_drive.py`, `tests/validation/test_geo36_rod_f_clinical_drive.py`, `tests/validation/test_geo38_f_clinical_3t_box_loaded.py`, `tests/validation/test_geo38_f_clinical_3t_rod_box.py`, and `docs/testing/xl-queue.d/` is unchanged. Quote the stat in the Outcome.
- **K1-12:** unset and `1` accepted; `2` / `build` / `x` refused. Each of the 19 foreign names is refused, matched by the message head and the name (C32's K1, `3t_drive.py:163–201`). The node knob is required under `1`.
- **K2-12:** `1` without the node raises a `ValueError` naming `FEM_EM_GEO36_S8_DRIVE_C_PF`, with `MeshGenerator.birdcage_port_domain` replaced by a raising sentinel (C32's K2, `:205–221`).
- **K7-12:** `1` plus node `c32.K_NODE_3T_PF` (44.0615) runs the gated test into `c21._builder_sentinel`. The sentinel must see:
  - rung Ø 0.704;
  - `implant_rod` == `c20._kwargs_3t()["implant_rod"]`, with `keep_cells` True;
  - f 128e6 on both `wf7` and `c4`;
  - records 1 128 430 / (1 508 256, 7 783 028).
  
  These are C32's values (`…C32-W0.log:67`). Inside the patched context, also assert `i35._rod_centre()` == `c20._rod_centre_3t()`. The unpatched centre is printed beside it; it is equal, since it is the same box.
- **K12-R, the anchor:** end to end on the 3 T stubs at node `c32.K_NODE_3T_PF`. `_run_heat_3t` returns green on its own collected gates: I36b's EM gates on both arms, (B), (E), (E′), the global (S2), (S0), and the re-posed (S2-tip) on the rod arm. Then:
  - **(Z-3T):** the (M2) / (Z) line's tokens equal C32's stub strings: Z_RC `2.000000000e-01-2.821960199e+01j`, Z_th `1.215760095e+01-2.127285736e+01j`, V_s `1.012556726e+00+1.771729051e+00j` (`…C32-W0.log:81`).
  - **(s2-3T):** both of these are ≤ `EXACT_IDENTITY_RTOL`; this is C32's K9-J comparison (`3t_drive.py:381–389`; C32 read 0.000e+00 on both, `…C32-W0.log:115`):
    - rod s2 against twin s2;
    - rod s2 against `h17._circuit_current_route_s(z_sym, [0, 1], z_rc, z_rc, 50.0)`.
    
    The frequencies seen during the run are (128e6, 128e6).
  - **(T-rep), I32b's stub literals by string** (`20261108T181211Z_GEO-36-I32b-L-W0b.log`):
    - `:103`: the share `1.041667e-02`, with P_vol / gel / rod `4.500000000e-03` / `4.453125000e-03` / `4.687500000e-05`;
    - `:104`: the census `{3: 10260, 403: 108}` and merged tag 3 `10368`;
    - `:106`, the twin: P_gel `1.125000000000e-03`, G `8.888888888889e+02`, G*P `1.000000000000e+00`, V_sub `9.999999999999e-04`, `exterior 1728 = outer 1728 + rod wall 0`;
    - `:108–109`, the twin's F2182 θ: `+2.121079410e-01` twice;
    - `:110–111`, the twin's ball θ: `6.062014386e-01` / `6.062135760e-01`;
    - `:112–113`, the twin's R_Pi dt1/10: `8.0980e-07` / `8.0969e-07`;
    - `:114–115`, the twin's probes: `+6.070276390e-01` / `+5.217085743e-01`;
    - `:117`, the rod: G*P `3.958333333333e+00`, V_sub `9.895833333333e-04`, `rod wall 144`;
    - `:119–120`, the rod's F2182 θ: `+8.484317694e-01` twice;
    - `:121–122`, the rod's ball θ: `2.424561997e+00` / `2.424610520e+00`;
    - `:123–124`, the rod's R_Pi dt1/10: `8.0963e-07` / `8.1058e-07`;
    - `:125–126`, the rod's probes: `+2.427815119e+00` / `+2.086793306e+00`;
    - `:128`: the ratio `a 3.999598, b 3.999598`;
    - `:130`: `failures none`.
  - **Not compared by string: the round-off residual tokens and every timing.** These are the (B) rel, (E), (E′) rel and (S2) tokens. They move with the call history at a fixed width: the same BE solve prints `(E) 1.521e-10` in K10-E and `BE 1.520e-10` in K10-R of one log (`L-W0b.log:76` against `:116`). They stay asserted at their bands inside `_run_heat`.
  - K10-R's six heads are present (`heating.py:900`), and every θ line carries `(1,+1j)`.
- **K12-11, I32b's stub literals on the composed route (`-n 2`):**
  - **Set-up:** K11's (`heating.py:946–1018`): `_k9_fixture`, `_columns_e_sup`, the merged arm, and two `_thermal_arm(..., tip_projection=True)` arms, u and r, with `_checkerboard` and `PROJECTION_ID_TOL` imported inside the test body. It runs in a context where both the 3 T stubs and `c21._c21_patches(m)` are live.
  - **Labels:** `K11_3T u` / `K11_3T r`. Each label gets its own PETSc prefix: `i32b_K11_3T_u`, with no hyphen.
  - **Asserted:**
    - (K11-Π): R_Π ≤ `S2_BAND` at both balls on arm r;
    - arm r's R_Π dt1/10 equals `6.855244e-07` / `8.521914e-07` by `.6e` string (`L-W0b.log:159–160`);
    - (K11-rep): |R_Π − `OPS83_R_PI_R`| ≤ `REP_TOL_4E`, and arm u's tip |R_Π − `OPS83_R_PI_U_TIP`| ≤ `REP_TOL_4E`. I32b's 4.402e-12 / 1.436e-12 / 3.954e-12 are printed beside (`:161`);
    - |pq_q| within `REP_TOL_3E`;
    - R_Π − R_DG0 on arm u ≤ `PROJECTION_ID_TOL`;
    - (E) and the global (S2) on both arms; (B) on arm u.

**Negative controls (asserted; each backed by a log).**
- **(N-64):** s2 against the 64 MHz route misses by more than `EXACT_IDENTITY_RTOL`. C32 read 2.515e-01 at this node on this seed-656 circuit (`…C32-W0.log:115`); printed beside.
- **(N-K7):** `heat._run_heat` without C21's patches reaches the builder at 0.710 m / 64 MHz / 1 217 025 cells (`…C32-W0.log:68`). The separation is categorical.
- **(N-P):** with `p_src_r2 = 1 + 1e-3`, `_run_heat_3t` raises `AssertionError`, and its failure list contains `gel twin (P) P2 9.990e-04`.
  - That is `_em_failures`' label (`heating.py:476`). C32 / I36d read 9.990e-04 (`…C32-W0.log:116`).
  - Ceiling: the perturbation itself, 1e-3/1.001. The precedent is C32's control (`3t_drive.py:372–375`).
- **(N-c):** K11-c on arm r. The retired DG0 `tip_s2["a"]` reads `1.5793e-03` = `15.79` × `S2_BAND` by string (`L-W0b.log:162`), asserted > `S2_BAND`.

**Windows.** Before each window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Every window is a bare foreground harness call, with Bash timeout 660000 ms, knobs unset, and no gmsh. The complex prefix and durable capture are `L-W0b.log:12`'s, verbatim, including the trailing `; exit $rc`.
- **W0b** `GEO-36-I85-W0b` (`-n 2`):
  `scripts/testing/run_and_log.sh GEO-36-I85-W0b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/validation/test_geo36_rod_f_clinical_3t_heating.py tests/validation/test_geo36_rod_f_clinical_heating.py tests/validation/test_geo36_rod_f_clinical_3t_drive.py -v -s --tb=short > /workspace/logs/geo36-i85-w0b-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo36-i85-w0b-raw.log; cat /workspace/logs/geo36-i85-w0b-raw.log; exit \$rc'"`
  - **Expected per rank, by token:** the new module 5 passed / 1 skipped; the heating module 9 passed / 1 skipped (`L-W0b.log:171`); C32's 4 passed / 1 skipped (its K1 / K2 / K7 / K9 plus the gated skip).
  - The two ride-alongs run unedited code: they are regressions, not edits.
- **W0n1** `GEO-36-I85-W0n1`: the new module alone at `-n 1`, `timeout -k 30 120`, with the label and raw-log name changed.
  - The literal strings print as "printed at -n 1 (record is -n 2)".
  - Every physics gate and control is asserted.

**Price.**
- **Measured:**
  - I32b's L-W0b took 22 s (`L-W0b.log:179`), and its `_run_heat` stub window 7.9 s (`:129`);
  - C32's W0 took 20 s with I36b's module beside it (this file's item 822 Outcome).
- **Predicted:**
  - the new module runs `_run_heat` twice (the good run and N-P) plus K11's two arms, ≈ 25–35 s;
  - W0b ≈ 60–90 s against 180; W0n1 ≈ 30–40 s.

**Tier.** Standard (the weekly's), complex, `-n 2` (plus `-n 1`). **Class:** chain step (no solve). **Cost:** ≈ 5.0 measured / 30 predicted min.

**Traps already paid for.**
- **Pytest's double import.** Patch through `sys.modules[heat._heat_columns.__module__]` (I32a's trap, this file's item 707).
- **Never import a `test_*` name.** Pytest would re-collect it here (I68′'s trap).
- **Patch the heating module's `_k9_stubs`, never `i36._k9_stubs`.** `heat._stubs` resolves `_k9_stubs` at call time from the heating module's namespace, which imported the name.
- **`_columns_e_sup` hard-codes the 64 MHz node for its own s2** (`heating.py:689–691`). K12-11 reads only its fields. Disclose this; never edit it.
- **`_k8_arm` / `_k9_stubs` hard-code 64 MHz and 193.5xxx pF.** `_k9_stubs_3t` overrides them after the call (item 822's trap).
- **`_projections` builds its PETSc prefix from the arm label** (`OPS-84`'s trap). A reused label collides, and a space or hyphen breaks it.
- **`ArityMismatch` is a `BaseException` raised on rank 0 alone** (item 849). Never put a `fem.form` refusal inside `pytest.raises`.
- Count verdicts per rank by token, and run with `-s`.
- **The leak guard.**
  - The new forms are the 3 T patch line and node 44.0615's (Z) tokens; C32's logs already carry the latter.
  - If the hook refuses, move the log unchanged to `logs/leak-held/`. Never `--no-verify`.

**Scope.** The 3 T body and its stub proof. **Not claimed:**
- any 128 MHz θ or tip figure (entry 67's);
- the 44.0525 pF node's printed forms (the READY review's leak check);
- entry 61's EM record;
- a band;
- a published comparand (W3);
- known-issues 2026-11-06 (it resolves on entry 63).

**Negative result.**
- **(K-q) red**, a queued module touched: revert, stop, BLOCKED.
- **Any other K-gate red:** park on `attempt/GEO-36-I85-*`, mark item 932 BLOCKED in the same commit, and leave entry 67 untouched.
- **A (T-rep) red** names its token. A thermal token that moves under the 3 T patch set is a coupling defect, and it is reported before entry 67 is priced.
- **An (s2-3T) red** means `_heat_columns`' s2 is not C32's route. That is a finding.
- **A crash in non-gate code before any gate value exists:** fix it and re-run once at the same `-t`, keeping the first log (the 10-31 narrow ruling).
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 I85 DONE.
- `GEO-36` gains a record line: the I32-3T body, stub-proved on C32's 3 T stub with I32b's literals.
- **In the landing commit, `docs/testing/xl-pending.md` entry 67's Status line becomes:**
  > `**PENDING on entries 61 / 63 only** — (i) **§10 I85 MET <date> (§9 item 932 landing commit; module `tests/validation/test_geo36_rod_f_clinical_3t_heating.py`, gated test `test_geo36_i32_3t_rod_heating`, knob `<name>=1` plus the node `FEM_EM_GEO36_S8_DRIVE_C_PF` (required); degree 1 only — the S8 / S11 knobs refused by its foreign guard; `<W0b log>`, `<W0n1 log>`)**`
  - Items (ii) and (iii) and the "first review after all three writes the command" sentence stay unchanged.
  - The old status line is kept as `*Was:*`. Entry 61's (i) is the form (`xl-pending.md:2449`); this follows step 6b rule 5.
- Queue nothing and write no `xl` command. The READY review writes it in entry 63's shape, with this module, the knob and node 44.0525.

**Depends on:** nothing (a new import-only module; entry 63's, 61's and 60's modules stay byte-identical).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 932 — `GEO-36` I85 — per `docs/planning/chunks/GEO-36.md` "§9 item 932". Create only `tests/validation/test_geo36_rod_f_clinical_3t_heating.py`; edit no landed module and no `src/` file. The heating module (entry 63's), C32's (entry 61's) and C21's (imported by entry 60) are queued windows' modules and stay byte-identical; quote `git diff --stat` in the Outcome. No gmsh, no clinical solve. Mark `xl-pending.md` entry 67 PENDING on 61 / 63 only in the landing commit; queue nothing. Before each window `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; pytest with -s; container `timeout -k 30`. Never move a band. You are the executor — do not spawn agents. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-I85-*`, mark item 932 BLOCKED in the same commit, and report."

**Outcome — 2026-11-10 (04:30 slot): GREEN, item 932 DONE.** Module `tests/validation/test_geo36_rod_f_clinical_3t_heating.py`; knob `FEM_EM_GEO36_S12_HEAT` (the recommended name); gated test `test_geo36_i32_3t_rod_heating`; 19-name foreign guard (C32's 17 + `FEM_EM_GEO36_S8_HEAT` + `FEM_EM_GEO36_S11_DRIVE`).
- **W0b** `20261110T105744Z_GEO-36-I85-W0b.log` (`-n 2`, 50 s against 180): 18 passed / 3 skipped per rank — new module 5 / 1, the heating module 9 / 1, C32's 4 / 1, as predicted. **W0n1** `20261110T105845Z_GEO-36-I85-W0n1.log` (`-n 1`, 21 s against 120): 5 passed / 1 skipped.
- **(K-q):** `git diff --stat HEAD` before the commit read `docs/testing/test-results.md | 3 +`, `.../test_geo36_rod_f_clinical_3t_heating.py | 481 +++` — "2 files changed, 484 insertions(+)" (plus the untracked logs); none of the five queued / landed modules, no `src/`, `docs/testing/xl-queue.d/` unchanged.
- **K1-12 / K2-12 / K7-12** green (`W0b.log:51–63`): rung Ø 0.704, C20's rod with `keep_cells` True, f 128e6 on `wf7` / `c4`, records 1128430 / (1508256, 7783028), `i35._rod_centre()` = `c20._rod_centre_3t()` = (0.105, 0.0, 0.0) in the patched context (unpatched the same, printed). (N-K7): 0.710 m / 64 MHz / 1 217 025.
- **K12-R** (`W0b.log:155–158`): `_run_heat_3t` green on its own gates; (s2-3T) rod vs twin 0.000e+00, rod vs the 128 MHz circuit route 0.000e+00; frequencies (128e6, 128e6) on both `_heat_columns` calls; (Z-3T) and all 22 (T-rep) lines match by string; K10-R's six heads present, every θ line carries `(1,+1j)`. (N-64) 2.515e-01 > band; (N-P) raised `['gel twin (P) P2 9.990e-04']`. Entry 55's 1.5 T ratio 1.088964 / 1.089826 printed beside the stub's 3.999598.
- **K12-11** (`W0b.log:180–188`): arm r R_Π dt1/10 `6.855244e-07` / `8.521914e-07` by string; (K11-rep) |diff| 4.402e-12 / 1.436e-12, u tip 3.954e-12 (= I32b's); (E), (S2), (B), (K11-Π), |pq_q|, identity green; (N-c) `1.5793e-03` = 15.79 × `S2_BAND`.
- **At `-n 1`** (`W0n1.log:137`, `:167`) the literal strings are printed, not asserted, and some miss in their last digits, as the width ruling expects: the twin's P_gel / G, the rod's V_sub, the four R_Π dt1/10 `.4e` tokens on the K12-R route, and K12-11's R_Π (|diff| 8.598e-10 / 5.407e-11 / 2.641e-10 — I32b's own `-n 1` figures). Every physics gate and control is asserted at `-n 1` and green.
- **Deviation (disclosed):** the first W0b (`20261110T105646Z_GEO-36-I85-W0b.log`, kept) failed on two test-code defects with no gate value at stake: K1-12 matched the guard message as a regex without `re.escape` (its parentheses), and K12-R's module resolution went through `_heat_columns.__module__` after the recorder had wrapped `_heat_columns`, so it resolved to this module (`AttributeError: ... no attribute '_run_heat'`). Both were fixed (`re.escape`; resolve through `_s0_failures`, which no stub patches) and re-run once at the same `-t` (the 10-31 narrow ruling). No band moved.
- Entry 67's Status set to PENDING on 61 / 63 only (old line kept as *Was:*). Nothing queued.

### §9 item 933 — 2026-11-10 (§10 I87, chain I — rule (4), the 11-09 weekly's step; chain step (no solve); write the anchor)

*Enumerated by the 11-09 weekly (`PROJECT_PLAN.md:4975–4985`). The 2026-11-10 03:00 daily review wrote this brief as execution detail.*

*It **writes the anchor and returns the step's three gates to the weekly**. As posed, (K-gel) is false by construction, and (K-hole)'s sign is predicted reversed on I32a's stub (premises below). They are printed, not asserted: asserting them would leave a known-red node on `main` (`PROJECT_PLAN.md:3985`). No band, target or done-when is set or moved. Family I87, attempt 1 of 4. **Executor:** the implementer.*

*It **depends on `xl` entry 63 having run**: take it no earlier than the 11-13 04:30 slot, with the 02:00 window's commit on `main`. This item edits entry 63's own module, which entry 60 (11-11) also imports (`test_geo38_f_clinical_3t_box_heating.py:36`). The 11-09 ruling "a queued window's module is not touched" (`PROJECT_PLAN.md:4824–4827`) binds, and a default-off knob does not lift it. The 11-11 weekly can rule the premises first.*

**Why the three gates cannot be asserted as posed (premises, for the weekly).**
1. **The module's default rod arm is not "gel".**
   - It excludes tag 403 from the thermal domain behind an adiabatic wall (I32a's detail 1, this file's item 707). The stub prints `V_sub 9.895833333333e-04` and `rod wall 144 (adiabatic)` (`20261108T181211Z_GEO-36-I32b-L-W0b.log:117`).
   - That is I69's hole in the ε → 0 limit, not a gel rod. An arm that puts the rod cells in the domain at gel properties adds their heat capacity, so "(K-gel) the gel arm reproduces the default's prints to the digit" fails by construction.
   - **So the default-off knob that keeps entry 63's path is unset = the present exclusion.** `gel` is a fourth configuration.
2. **On the stub, the in-domain arms differ by heat capacity, below `INTERIOR_BAND`.**
   - **The stub geometry.** The stub's "rod" is a 108-cell block, the cells with midpoints at |x|∞ < 12 mm (`test_geo36_rod_f_clinical_drive.py:1194`). It is 1.0417e-05 m³, 1.04 % of the box (`L-W0b.log:103`, `:106`, `:117`). Its tip balls sit 42 mm off, at `_k9_fixture`'s (±0.03, ∓0.03, 0) (`drive.py:1201`). The source is uniform in the gel, and the block is unheated.
   - **The prediction (unlogged drafter arithmetic, a point-sink estimate, ±50 %).** An in-domain block is a heat sink of strength ≈ ρc_rod·V_b·dθ/dt. At the balls, the gel arm's deficit against the exclusion is ≈ q·V_b/(4πk r)·erfc(r/(2√(κt\*))) ≈ 1.8e-2 K on 2.42 K, i.e. 0.4–0.7 %. Inputs: q 4.0e3 W/m³ (`L-W0b.log:123`), k 0.6, r 42.4 mm, √(κt\*) 20.1 mm.
   - **Ti's sink is ≈ 0.57 of gel's.** Ti's ρc is 0.557 of gel's (`ROD_RHO_C` 2.33018e6 against `RHO_C` 4.18e6), and Ti equilibrates in ≈ 50 s.
   - **Predicted ordering at the balls: hole ≈ default > Ti > gel.** Ti/gel − 1 ≈ +1.5e-3 … +3.5e-3, and hole/gel − 1 ≈ +4e-3 … +7e-3.
   - **So (K-Ti) keeps its sign but not the band, and (K-hole) is reversed.** I69's source sits at the tips, a different configuration. It reads Ti > gel > hole at z = 0, but hole > gel > Ti at the tip point (2, 15) mm (`IMP-8.md:919–921`). The weekly's "(I69's sign)" matches the z = 0 row.
3. **(S2-tip) on balls that contain rod cells.** This is a finding for the Ti clinical variant; the stub cannot reach it.
   - With the rod in the domain at ρc ≠ `RHO_C`, the BE identity becomes ⟨θ¹⟩/Δt ≈ ⟨M_ρc⁻¹F⟩, not the scalar-ρc ρc⟨θ¹⟩/Δt ≈ ⟨Πq⟩.
   - The clinical tip ball (a 6.2035 mm, centred on the rod's end) holds ≈ 22 mm³ of rod, ≈ 2.2 % of its 1000 mm³ (drafter arithmetic: π·1²·6.2 + 2π/3).
   - So `_tip_slope_rows`' scalar form is predicted red there at ≈ 1e-3 … 1e-2 on the Ti and hole arms. It stays exact on the gel arm.
   - The stub's balls contain no rod cells. **The Ti clinical window is not priceable until the weekly rules this comparand.**

**Code.** Edit `tests/validation/test_geo36_rod_f_clinical_heating.py` **by addition only**.
- **Byte-identical:** `_thermal_arm`, `_thermal_failures` and every existing test. P61 and C31 call the first two, and P61 intercepts `i32._transient` (`test_geo38_f_clinical_box_heating.py:140–158`). No `src/` edit; `test_geo36_rod_f_clinical_drive.py` is untouched.
- **The knob.**
  - `ROD_THERMAL_ENV = "FEM_EM_GEO36_S8_ROD_THERMAL"` and `ROD_THERMAL_ARMS = ("gel", "ti", "hole")`.
  - `_rod_thermal()`: unset or empty ⇒ `None` (the present path); a listed arm ⇒ itself; anything else ⇒ `ValueError`.
  - Read it once, in `_run_heat`.
- **Values: I66's and I69's, imported inside the function body** (see the cycle trap).
  - `ROD_K` 6.7 and `ROD_RHO_C` 2.33018e6, both nominal, from `test_imp8_composite_cylinder_closed_form` (`:52–53`).
  - `EPS_HOLE` 1e-6 from `test_imp8_rod_conducting_body_fem` (`:84`).
  - The coefficient rule is I66's `_coefficients` (`:149–166`): k = `K_GEL` and ρc = `RHO_C` everywhere. On rod cells, `ti` sets `ROD_K` / `ROD_RHO_C`; `hole` sets `EPS_HOLE`·`K_GEL` / `EPS_HOLE`·`RHO_C`; `gel` overwrites nothing.
- **`_thermal_arm_rod_cells(label, msh, cell_tags, e_sup, sigma_f, G, centre, tips, probes, *, rod_thermal, f2182=None)`, new.** It mirrors `_thermal_arm`'s order line for line:
  - P_gel and G as there;
  - the F2182 CN before the t\* CN;
  - the t\* CN with its two snapshots;
  - BE Δt₁, then BE Δt₁/10;
  - one ball batch, then one probe batch.
  
  **The changes:**
  - **Domain:** `_tag_submesh(msh, _merge_rod_into_phantom(msh, cell_tags), PHANTOM_CELL_TAG)`, i.e. tags 3 + 403.
  - **Source:** q is `_heat_source(msh, cell_tags, …, sub, parent_ids)` on the **original** tags, so it is zero on the rod cells. `_cell_average_dg0_tag` zeroes cells off the tag (`test_mat15_coil_driven_heating_phantom.py:139–140`).
  - **Rod cells and coefficients:** rod cells are those whose parent tag is 403, using `_tag_submesh`'s `parent_values` pattern (`:111–113`). DG0 `k` and `rho_c` live on the submesh.
  - **Solves:** every solve is the call `_transient` wraps (`:145–151`; import the solver from where that module imports it), with `k=` / `rho_c=` set to the DG0 Functions and every other argument unchanged.
  - **(E′)ρc:** ∫ρc_DG0(θᴺ − θᴺ⁻¹)/Δt + ∮hθ^(n+½), against G·P.
  - **(S2)ρc:** |∫ρc_DG0 θ¹ dV / (Δt₁·G·P) − 1|.
  - **Ball reader:** gains a rod-indicator field, so V(B ∩ rod) is printed per ball.
  - **(S2-tip):** `_projected_ball_comparands` plus `_tip_slope_rows` with `RHO_C`, exactly as `_thermal_arm` computes it (premise 3: the form is the weekly's).
  - **Label:** `f"rod {rod_thermal}"`, so the tag is `i32b_rod_ti`.
- **`_rod_cells_failures(arm)`, new:**
  - (B), (E), (E′)ρc, (S2)ρc;
  - wall facets == 0 (the rod is interior now);
  - the projected (S2-tip) ≤ `S2_BAND`.
- **`_run_heat`.**
  - `rt = _rod_thermal()` at the top.
  - When `rt is None`, the rod-arm call and `_thermal_failures(heat_rod)` are the present lines, untouched.
  - Otherwise the rod arm is `_thermal_arm_rod_cells(…)`, with failures from `_rod_cells_failures`.
  - Print one line: `rod thermal: {rt or 'excluded (adiabatic wall; entry 55 / 63 route)'}`.
  - **No knob value changes the twin arm, G, the EM half, or the F2182 / probe placement.**

**Gates (asserted; bands imported, none set).**
- **(K-def), the weekly's "entry 63's queued command is unchanged":**
  - **(a)** Before the commit, `git diff --stat HEAD` lists only `tests/validation/test_geo36_rod_f_clinical_heating.py` among the tests. `docs/testing/xl-queue.d/` is unchanged. The Grep tool finds no `ROD_THERMAL` in `docs/testing/xl-queue.d/` or `xl-pending.md` entry 63; that command sets only `FEM_EM_GEO36_S8_HEAT=1` and the node.
  - **(b)** A new coding test runs `_run_heat` on `_stubs` with the knob **unset**. Every value token item 932's (T-rep) lists (`L-W0b.log:103–130`) equals the log's string at `-n 2`. Residual and timing tokens are excluded, as there.
  - **(c)** Wreg and Wimp are green.
- **(K-arm), per arm `gel` / `ti` / `hole`.** One new coding test sets `monkeypatch.setenv(ROD_THERMAL_ENV, arm)` and runs `_run_heat` on `_stubs` at `K9_NODE_PF`:
  - `_run_heat`'s collected gates are green, including the rod arm's (B), (E), (E′)ρc, (S2)ρc and projected (S2-tip), and the twin's as now;
  - the rod arm's G·P equals the default's `3.958333333333e+00` by string (`L-W0b.log:117`), so the source is unchanged by the knob;
  - wall facets are 0, and V_sub equals the twin's `9.999999999999e-04` by string (`:106`);
  - V(B ∩ rod) = 0 at both stub balls. This is the precondition of the scalar (S2-tip) form on the stub (premise 3).
- **(K-gel=):** on the gel arm the DG0 forms equal the scalar ones, since ρc is uniform:
  - |(E′)ρc − (E′)scalar| ≤ `E_BAND`;
  - |(S2)ρc − (S2)scalar| ≤ `S2_BAND`.

**Negative controls (asserted; ceilings computed).**
- **(N-src):** q built on the merged tags (rod cells heated) misses (B)'s G·P by share/(1 − share) = 1.052632e-02.
  - The share `1.041667e-02` is logged (`L-W0b.log:103`); the share itself is the ceiling.
  - Asserted > `B_BAND`.
- **(N-ρc):** on `ti` and `hole`, the scalar-`RHO_C` (E′) misses `E_BAND`. This shows the coefficient reached the storage term.
  - Predicted ≈ (1 − ρc_rod/ρc)·(V_b/V)·(stored / G·P): ti ≈ 0.4425 × 0.010417 × 0.711 ≈ 3.3e-3, hole ≈ 7.4e-3.
  - Stored / G·P is 2.814/3.958 (`L-W0b.log:127`).
  - This is unlogged arithmetic, so it is asserted only as > `E_BAND` (rule (e)), with the factor printed.
- **(K-c):** a copy of `L-W0b.log:121` with one θ digit changed, made in memory, fails (K-def)(b), as I68′'s (K-c) does.

**Printed: the weekly's three gates as posed (record; predictions beside; not asserted).** Per arm and ball, the rod-arm tip θ(t\*) and the F2182 θ(900 s), with:
- (K-gel) |gel/default − 1|, predicted 4e-3 … 7e-3;
- (K-Ti) Ti/gel − 1, predicted +1.5e-3 … +3.5e-3;
- (K-hole) hole/gel − 1, predicted +4e-3 … +7e-3, with the sign reversed;
- hole/default − 1, predicted ≤ 1e-4 (I69's ε-hole against the exact exclusion);
- the measured ordering line;
- V(B ∩ rod) and the (S2-tip) rows per arm.

**Windows.** Before each window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Every window is a bare foreground harness call with Bash timeout 660000 ms, and no knob is set on any command line.
- **W0b** `GEO-36-I87-W0b`: L-W0b's command (`L-W0b.log:12`) with the label and raw-log name changed (`-n 2`, `timeout -k 30 180`, the heating module alone). Expected: 9 plus the new tests passed / 1 skipped per rank.
- **W0n1** `GEO-36-I87-W0n1`: the same at `-n 1`, `timeout -k 30 120`. The literals are printed only.
- **Wreg** `GEO-36-I87-Wreg`: L-Wreg's command byte for byte except the label and raw-log name (`20261108T181335Z_GEO-36-I32b-L-Wreg.log:12`: `OPS-83`, `OPS-84`, C31, P61). Expected: **18 passed / 9 skipped per rank** (`:377`, `:380`), ≈ 37 s (`:385`).
- **Wimp** `GEO-36-I87-Wimp`: the real build, without the complex prefix, `-n 1`, `timeout -k 30 120`:
  `python3 -m pytest --collect-only -q tests/validation/test_imp8_rod_line_source_fem.py tests/validation/test_imp8_entry48_lumped.py tests/validation/test_geo36_i76_f2182_table.py tests/validation/test_mat21_clinical_probe_reader.py tests/validation/test_geo36_rod_f_clinical_heating.py`
- **Wimp-c** `GEO-36-I87-Wimp-c`: the same in the complex build, plus `tests/package/test_wf15_study_rod_heating.py`.

Both Wimp windows expect 0 collection errors.

**Price.**
- **Measured:** L-W0b took 22 s (`L-W0b.log:179`). One `_run_heat` stub takes 7.9 s (`:129`).
- **Predicted:**
  - the new tests run `_run_heat` four times (default plus three arms) plus the controls, so W0b ≈ 60–80 s;
  - Wreg ≈ 40 s;
  - each Wimp ≈ 10–20 s.

**Tier.** Standard, complex, `-n 2` (plus `-n 1` and one real collect-only). **Class:** chain step (no solve). **Cost:** ≈ 5.0 measured / 40 predicted min.

**Traps already paid for.**
- **Import cycle.** `test_imp8_rod_line_source_fem` imports this module, and `test_imp8_rod_conducting_body_fem` imports that one, so a top-level `EPS_HOLE` import is a cycle into a half-initialised module. Import inside the function body, as `_projected_ball_comparands` imports `OPS-83` (`heating.py:267–268`).
- **The module must still import in the real build.** I58, I55 and I79 import it there; Wimp is the check.
- **P61's `_transient` intercept keeps the last θ = 0.5 call** (this file's item 870). The knob path is reached only when the knob is set, and P61 never sets it.
- **The PETSc prefix comes from the label.** Never put parentheses, a colon or a hyphen in one.
- **`cell_tags.values` and the rod-cell selection are rank-local.** Reduce before printing or asserting.
- **A new ball field** uses `ufl.lt(ufl.real(…), …)` and a pinned quadrature degree (`_ball_reader` does both).
- **Collect every gate, then assert** (item 668). Count verdicts per rank by token. Run with `-s`.
- **Literal reproductions are asserted at `-n 2` only** (`PROJECT_PLAN.md:3966`).
- **The leak guard.** If the hook refuses, move the log unchanged to `logs/leak-held/`. Never `--no-verify`.

**Scope.** The knob, its in-domain thermal arm and its stub identities. **Not claimed:**
- any clinical θ;
- the weekly's (K-gel), (K-Ti) and (K-hole) as gates;
- a ruling on the (S2-tip) comparand for balls containing rod cells (premise 3);
- the Ti clinical window's price;
- any Ti property (the values are nominal, from memory).

**Negative result.**
- **(K-def), Wreg or Wimp red:** the default path or an importer moved. Revert, stop, BLOCKED.
- **(K-arm) or a control red:** park on `attempt/GEO-36-I87-*`, mark item 933 BLOCKED in the same commit, and print every value.
- **A printed (K-…) value outside its predicted range** is a record for the weekly, not a stop.
- No `PREDICTED RED:` line: the three gates are printed, not asserted.

**Status it can move.**
- `GEO-36` gains a record line: the `ROD_THERMAL` knob (default off), stub-proved, with the three arms' tip θ beside the default's.
- §10 I87 stays OPEN, with the anchor written. The three gates and premise 3 go back to the weekly.
- `xl-pending.md` is untouched, and nothing is queued.

**Depends on:** `xl` entry 63 having run (the 11-13 02:00 window committed). Take it no earlier than the 11-13 04:30 slot.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 933 — `GEO-36` I87 — per `docs/planning/chunks/GEO-36.md` "§9 item 933". First confirm `xl` entry 63's window (11-13) is committed on `main`; if not, stop. Edit only `tests/validation/test_geo36_rod_f_clinical_heating.py`, by addition: `_thermal_arm`, `_thermal_failures` and every existing test stay byte-identical; no `src/` edit. The three (K-…) gates are printed, never asserted. Before each window `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; pytest with -s; container `timeout -k 30`. Never move a band. You are the executor — do not spawn agents. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-I87-*`, mark item 933 BLOCKED in the same commit, and report."

**Review ruling — 2026-11-10 03:00 daily, on item 933: not queued; returned to the 11-11 weekly as "write the anchor".** The brief above would land the knob while asserting gates other than the step's own. It asserts the default path, the per-arm identities and two controls, and prints the weekly's (K-gel), (K-Ti) and (K-hole). That changes the step's validation target, which a daily may not do (§9 rule (4)). Hand estimates, with no log behind them, predict the step's gates fail as posed on the stub:
- **(K-gel)** cannot reproduce the default "to the digit". The default arm excludes tag 403 behind an adiabatic wall (`20261108T181211Z_GEO-36-I32b-L-W0b.log:117`, `rod wall 144`), so a gel arm adds the rod's heat capacity.
- **(K-Ti)** keeps its sign but is predicted below `INTERIOR_BAND`.
- **(K-hole)** is predicted with the opposite sign: on the stub the rod is unheated and the source is uniform in the gel.

I69's ordering is row-dependent: Ti > gel > hole at z = 0, but hole > gel > Ti at the tip point (`IMP-8.md:919–921`). The Ti clinical variant's (S2-tip) is also predicted red on a tip ball that is ≈ 2.2 % rod. **Nothing is lost by waiting:** the item could not run before entry 63's window (11-13 02:00, the queued-module rule), and the 11-11 weekly comes first. The weekly re-poses the gates, or rules the brief's printed form, and the next review queues it. The brief stays as the proposal.

### §9 item 946 — 2026-11-10 (daily, 2026-11-10 — provisional; chain I step I88)

*Written by the 2026-11-10 10:30 daily review under §9 rule (6), provisional; the 11-11 weekly ratifies or reverts it. The target is entry 67's pre-registered readout, unmoved: "the 3 T / 1.5 T ratio of the tip θ per accepted watt (C40's twin, a record)" (`docs/testing/xl-pending.md:2569`; `PROJECT_PLAN.md` §10 chain I). Family I88, attempt 1 of 4. Chain I's rule-(6) count is 1 of 6. **Executor:** implementer. **Depends on nothing:** it reads logs already on `main` and is one new import-only module. The 03:00 drafter's "I88" (a reader for entry 61) duplicated C40 / C45 and is withdrawn. Item 906's (Stub61) is already proved on C32's stub (`GEO-38.md:3425`).*

**Executes.** I88 builds the reader that turns entry 67's readout into one knob, as C45 did for C40 and I79 for I76.
- It computes the rod-arm and gel-twin tip 1 g θ(t\*) per accepted watt at each field strength, and their 3 T / 1.5 T ratio.
- **The normaliser is P61's** (`test_geo38_f_clinical_box_heating.py:21–26`): θ_acc = θ · P_gel / P_acc(cw), with P_acc the lumped line cut (ruling 2026-10-12, `_p_acc_spec`).
- **It reaches the rod route by string custody between logs.** The I32 heating log prints θ at G = 1 / P_gel(twin), with no S-matrix and no P_acc. Entry 48 (same mesh, drive and node) prints S2 on both arms. Entry 59 prints the twin's cw factor.
- **1.5 T** is computed now from entries 55, 48 and 59. Entry 55's figure is a pre-figure; entry 63 supersedes it.
- **3 T** is stub-proved now. Its real figure needs entries 60 (the factor), 61 (S2) and 67 (θ).

**Facts verified.** These come from the drafter; the review re-traced the two digits marked ✓.
- **Sense.** The I32 route runs (1,+1j) (`GEO-36.md:1236`), which C6a labels **cw** (`20261103T080009Z_GEO-38-C6.log:12287`; `GEO-38.md:1981`). Entries 52 / 58's P_acc are **ccw** and are never used here.
- **Entry 59** (`20261110T080008Z_GEO-38-P33.log`):
  - cw `:12343`: P_gel `2.538800824240e-03` W, P_acc `6.095252867444e-02` W, factor `4.165210008e-02` ✓;
  - ccw `:12325`: factor `5.285093370e-02` ✓;
  - (Z) `:12065`.
- **Entry 55** (`20261106T080008Z_GEO-36-step10-I32.log`):
  - (M2)/(Z) `:12062`: Z_RC `2.000000000e-01-1.284504713e+01j`, Z_th `3.259458821e+00-1.195984969e+01j`, V_s `1.060593121e+00+3.891607472e+00j`;
  - rod share `:12497`;
  - twin P_gel `:12499` `2.538800824240e-03`, which equals entry 59's cw P_gel by string ✓;
  - twin θ `:12501–12502`; rod θ `:12512–12513`;
  - gates `:12523` `['rod (S2-tip) 6.420e-03']`.
- **Entry 48** (`20261102T080008Z_GEO-36-step9-I36.log`):
  - (Z) `:12062`: the same three tokens;
  - arm M S2 `:12491`; arm R S2 `:12495`;
  - arm M P_vol / gel / rod `:12494` = entry 55 `:12497`;
  - arm R P_vol `:12498` `2.538800824e-03`.
- **3 T stubs:**
  - C32 `20261105T133319Z_GEO-36-C32-W0.log`: (Z) `:81`; S2 arm M `:82` = arm R `:86`; P_vol `:85` / `:89`;
  - C31 `20261105T181050Z_GEO-38-C31-W0r.log`: (Z) `:87` = C32 `:81`; S2 `:91` = C32 `:86`; cw P_acc `:117` `1.485167398058e-02`;
  - I85 `20261110T105744Z_GEO-36-I85-W0b.log`: green K12-R block `:80`–`:115` ((Z) `:81`, rod share `:88`, twin P_gel `:91`, `failures none` `:115`); `:118`–`:153` is the N-P block, red by design;
  - entry 63's form, I32b: `20261108T181211Z_GEO-36-I32b-L-W0b.log:95`–`:130`.
- **Imports (module objects only):**
  - `_p_acc_spec(s2, w, v_s, z0)` (`test_geo38_f_clinical_box_loaded.py:510–515`);
  - C45's `parse_i36_log`, `C32_STUB_LOG` / `C32_STUB_LINE = 81`, `ENTRY61_GLOB` (`test_geo36_c40_rod_frequency_reader.py:44`, `:53–54`, `:130–210`);
  - C46's `parse_heat_log`, `ENTRY59_GLOB` / `ENTRY60_GLOB` (`test_geo38_c39_thermal_reader.py:52–54`, `:144–145`, `:158`);
  - C33's (D-Pacc) form, rel ≤ 1e-8 (`test_geo38_c33_field_strength_reader.py:161–164`).

**Code.** One new module, `tests/validation/test_geo36_i88_rod_heat_per_acc_reader.py`. No existing file and no `src/` file is edited.
- **Knob `FEM_EM_GEO36_I88`:**
  - unset ⇒ the gated read skips and the coding tests run;
  - `1` ⇒ entries 48 / 55 / 59 are required; 60 / 61 / 63 / 67 are optional, and each absent one prints "pending" with nothing written;
  - `both` ⇒ all of 48 / 59 / 60 / 61 / 63 / 67 are required, 63 replaces 55 as the 1.5 T θ source (55 printed beside), and rank 0 writes `/workspace/logs/geo36/i88-rod-heat-acc-1p5t-3t.json`;
  - anything else ⇒ `ValueError`.
- **Globs:** `ENTRY63_GLOB = "*_GEO-36-step10-I32-prime.log"`. `ENTRY67_GLOB = "*_GEO-36-step12-I32-3T.log"` is provisional: entry 67's READY review confirms it against its `XL_CHUNK` or edits this one constant.
- **`parse_i32_block(path, at_line=None)`:**
  - a block runs from a `start: -n N, I32 body` header to the next `[GEO-36 I32] gates: failures` line;
  - tokens are kept as strings;
  - one θ regex covers the I32a and I32b forms, through `re.search`;
  - a missing or repeated line raises `ValueError`.
- **Figures:**
  - twin θ_acc = θ × factor(cw);
  - rod θ_acc = θ × P_gel_twin / P_acc_rod, with P_acc_rod = `c6._p_acc_spec(S2_M, (1, 1j), V_s, 50.0)`. The twin-P_acc version and P_acc_rod / P_acc_twin − 1 are printed beside it;
  - the 3 T / 1.5 T ratio per tip and arm, ≤ 7 digits, with the label `(1,+1j) = cw` on every line.

**Gates.** No band is set. Each check is string equality, or agreement to half a unit of the printed digit.
- **(X48) / (X55) / (X59):** one block each, node `193.5996`; entry 55's balls `top` / `bottom` and every θ line `(1,+1j)`; entry 55's failures string equals `['rod (S2-tip) 6.420e-03']` (custody of a red window, so its figure is labelled a pre-figure).
- **(Cust-Z):** Z_RC / Z_th / V_s are string-equal across entry 55 `:12062`, entry 48 `:12062` and entry 59 `:12065`.
- **(Cust-P):**
  - entry 55 `:12499` == entry 59 cw P_gel `:12343`;
  - entry 55 `:12497` == entry 48 arm M `:12494`;
  - entry 48 arm R P_vol against entry 55's twin P_gel, within half a unit of its 10th digit.
- **(D-Pacc):** `c6._p_acc_spec(S2_R(48 :12495), (1, 1j), V_s(48), 50.0)` against entry 59's cw P_acc at C33's 1e-8.
- **(D-f):** entry 59's cw factor equals P_gel / P_acc within C46's half-unit bound.
- **(Stub-3T):**
  - (Z) C32 `:81` == I85 `:81` == C31 `:87`;
  - C32 arm R S2 `:86` == C31 `:91`;
  - `_p_acc_spec` against C31 `:117` at 1e-8;
  - P_acc_rod / P_acc_twin − 1 == 0.0 exactly (the S2 strings are identical);
  - I85 `:88` == C32 `:85`;
  - I85 `:91` against C32 `:89` within half a unit;
  - I85 `failures none`.

  The stub table is printed as "stub, form proof, not a record". **(Cust-P) at 3 T cannot be stub-proved,** because C31's P61 stub feeds a different box field. Its first real check is entry 60 against entry 67 at `both`; the Outcome says so.
- **(Form-63):** `parse_i32_block` on I32b L-W0b at `at_line=95` parses, and its θ tokens equal `:121–122` / `:110–111`.
- **(Self):** the 1.5 T side read against itself gives every ratio `1.000000e+00` exactly.
- **(Pend):** with entries 60 / 61 / 63 / 67 absent, it prints "pending" and writes nothing. The check monkeypatches the log roots to empty `tmp_path` subdirectories (item 918) and compares rank 0's JSON `(mtime_ns, bytes)`.
- **(K-q):** `git diff --stat HEAD` before the commit lists only the new module and the record files. It names none of: the three queued modules, C45's, C46's, C6a's, `src/` or `docs/testing/xl-queue.d/`. The Outcome quotes it.

**Negative controls.**
- **(N-sense):** entry 55's twin θ normalised with entry 59's **ccw** factor. The separation is 5.285093370e-02 / 4.165210008e-02 − 1 = **0.268866**, fixed by two logged factors, so it is the ceiling. The check asserts the separation exceeds the two tokens' combined half-unit bound (rule (e)) and prints 0.268866.
- **(N-drive):** (Cust-Z) on C32 `:81`'s 3 T tokens against entry 48's tokens fails.
- **(N-digit):** a one-digit in-memory change to entry 55's `:12499` token fails (Cust-P).
- **(N-blocks):** `parse_i32_block` on I85's W0b without `at_line` raises `ValueError`.
- **Printed, predicted, unlogged** (the drafter's arithmetic):
  - P_acc_rod / P_acc_twin − 1 ≈ −5.8e-6;
  - twin ≈ 1.2839e-03 / 1.2874e-03 K/W;
  - rod ≈ 1.3981e-03 / 1.4030e-03 K/W (top / bottom).

**Windows.** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0 before each window. Each is one bare foreground harness call, Bash timeout 660000 ms. C45 and C46 ride along at W0 as regressions (priors 15 / 1 and 13 / 1 per rank, `GEO-38.md:3650`).
- **W0** (knob unset, `-n 2`):
  `scripts/testing/run_and_log.sh GEO-36-I88-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_i88_rod_heat_per_acc_reader.py tests/validation/test_geo36_c40_rod_frequency_reader.py tests/validation/test_geo38_c39_thermal_reader.py -v -s --tb=short'"`
- **W1** (the real 1.5 T read, `-n 2`):
  `scripts/testing/run_and_log.sh GEO-36-I88-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_GEO36_I88=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_i88_rod_heat_per_acc_reader.py -v -s --tb=short'"`
- **W1-n1:** W1 at `mpiexec -n 1`, `timeout -k 30 120`, label `GEO-36-I88-W1-n1`.

**Tier.** Standard, complex, `-n 2` (plus `-n 1`). Predicted 5–30 s per window, from C45's 2 / 3 s plus parsing three ≈ 12 500-line logs. **Class:** chain step (no solve). **Cost:** ≈ 6.5 measured / 30 predicted min.

**Traps already paid for.**
- (1,+1j) is C6's **cw**; never borrow the ccw P_acc.
- Never import a `test_*` name (I68′).
- Interleaved pytest tokens sit on log lines (entry 59 `:12343` begins `PASSED…`): use `re.search`, never `^` or `startswith`.
- I85's W0b holds five I32 headers. The green block is `at_line=80`.
- Logs gzip at 7 days: read through `_read_log`, and never gunzip a record.
- The host JSON is never a source.
- Premise tests read a tmp log root (item 918). `find_log` raises on two matches.
- Parse on rank 0 and broadcast. Count verdicts per rank by token, and use `-s`.
- The leak guard: these are new digit strings. If the hook refuses, move the log unchanged to `logs/leak-held/`; never `--no-verify`.
- Queued windows' modules stay byte-identical (11-09 ruling (c)).
- Never delete a log. A non-gate crash before any gate value may be fixed and re-run once at the same `-t`, keeping the first log.

**Scope.** A reader and its stub proof. **Not claimed:**
- any 128 MHz θ or 3 T P_acc;
- entry 55's figure as a record (that window was red on (S2-tip); entry 63 decides);
- a band, a tissue model, compliance, or the published comparand (W3).

**Negative result.** Report and stop, BLOCKED in the same commit, parked on `attempt/GEO-36-I88-*`.
- (Cust-Z) / (Cust-P) red ⇒ the rod route's twin is not entry 59's fixture or drive. That is a finding for the weekly before entry 67 is priced.
- (D-Pacc) red ⇒ S2 → P_acc does not reproduce P61's print.
- `c46.parse_heat_log` raising on entry 59 ⇒ report it before entry 60's read.
- (Form-63) red ⇒ the reader cannot read entry 63's form.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 I88 DONE.
- `GEO-36` gains a record line: the clinical rod's tip 1 g θ per accepted watt at 1.5 T, as a pre-figure. It supersedes I68′'s "per accepted watt: not in a tracked record" (`GEO-36.md:1735`, `:1802`), for the weekly to ratify as a route.
- Entry 67's readout becomes one `FEM_EM_GEO36_I88=both` window once entries 60 / 61 / 63 / 67 have landed (11-17 at the earliest).
- `xl-pending.md` is untouched.

**Review notes (route readings for the weekly; no target moved).**
1. Carrying P61's normaliser to the rod by custody is a route choice.
2. The rod arm is normalised by its own P_acc, with the twin-P_acc version printed beside it (they differ by ≈ 6e-6, predicted).
3. (Cust-P) at 3 T is deferred to entries 60 / 67.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 946 — `GEO-36` I88 — per `docs/planning/chunks/GEO-36.md` "§9 item 946". Create only `tests/validation/test_geo36_i88_rod_heat_per_acc_reader.py`; edit no existing module and no `src/` file. `test_geo36_rod_f_clinical_heating.py`, `test_geo38_f_clinical_3t_box_loaded.py` and `test_geo36_rod_f_clinical_3t_drive.py` are queued windows' modules and stay byte-identical; quote `git diff --stat` in the Outcome. No solve, no mesh. Import module objects, never a `test_*` name. Before each window `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; pytest with -s; container `timeout -k 30`. Never move a band. Never delete a log. You are the executor — do not spawn agents. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-I88-*`, mark item 946 BLOCKED in the same commit, and report."

**Outcome (12:00 slot 2026-11-10) — DONE, green; a reader and its stub proof, 1.5 T figure a pre-figure.** New module `tests/validation/test_geo36_i88_rod_heat_per_acc_reader.py` (imports `c45` / `c46` / `c6` as module objects; no existing file edited). Windows, each preceded by `pgrep -c python3` = 0: **W0** `20261110T180513Z_GEO-36-I88-W0.log` (5 s) 1 failed / 34 passed / 3 skipped per rank — the red was the reader's own line-pin, not a gate value: (Stub-3T)'s P_acc was read through `c46.parse_heat_log(C31, at_line=130)`, whose cw sources line is `:174`, while the brief's `:117` (with its `:87` (Z) and `:91` S2) sits in C31's **cw-only** block (header `:83`), which `parse_heat_log` refuses by design; the value already matched (rel 2.682e-10, same token string at `:117` and `:174`). Fixed by reading `:117` directly with C46's `_RX_P61_SRC` (`parse_p61_cw_sources`) and re-run once at the same `-t`, first log kept: **W0r** `20261110T180551Z_GEO-36-I88-W0r.log` (4 s) 35 passed / 3 skipped per rank (complex-mode 4, I88 11 + gated skip, C45 11 + skip, C46 9 + skip — C45 / C46 at their 15 / 1 and 13 / 1 priors); **W1** `20261110T180613Z_GEO-36-I88-W1.log` (3 s, `FEM_EM_GEO36_I88=1`, `-n 2`) 16 passed per rank; **W1-n1** `20261110T180623Z_GEO-36-I88-W1-n1.log` (4 s) 16 passed, every figure identical to `-n 2`. Measured: (X48) / (X55) / (X59) green (node `193.5996`, balls `top` / `bottom`, every θ line `(1,+1j)`, entry 55 failures `['rod (S2-tip) 6.420e-03']`, rod G == twin G, the brief's eleven cited line numbers all hold); (Cust-Z) string-equal 55 `:12062` / 48 `:12062` / 59 `:12065`; (Cust-P) 55 `:12499` == 59 `:12343`, 55 `:12497` == 48 `:12494`, 48 arm R P_vol vs twin P_gel \|d\| 2.400e-13 (≤ 5.000e-13); **(D-Pacc) rel 1.082e-10** (≤ 1e-8); (D-f) rel 9.329e-11 (≤ 1.203e-10). **1.5 T θ per accepted watt, (1,+1j) = cw, pre-figure (entry 55 red on (S2-tip); entry 63 decides):** twin 1.283915e-03 / 1.287357e-03 K/W, rod 1.398145e-03 / 1.403004e-03 K/W (top / bottom; twin-P_acc version 1.398137e-03 / 1.402995e-03), P_acc_rod / P_acc_twin − 1 = −5.850512e-06 — the drafter's predictions held. This supersedes I68′'s "per accepted watt: not in a tracked record" as a route reading for the weekly to ratify. (Stub-3T) green on C32 `:81` == I85 `:81` == C31 `:87`, C32 arm R S2 == C31 `:91`, `_p_acc_spec` vs C31 `:117` rel 2.682e-10, rod / twin P_acc − 1 = 0.0 exactly, I85 `:88` == C32 `:85`, I85 `:91` vs C32 `:89` \|d\| 0, I85 `failures none` (stub, form proof, not a record). **(Cust-P) at 3 T is not stub-provable** (C31's P61 stub feeds a different box field); its first real check is entry 60 against entry 67 at `both`. (Form-63) green at `at_line=95` (rod `:121–122`, twin `:110–111`); (Self) all four ratios exactly 1; (Pend) entries 60 / 61 / 63 / 67 pending on empty roots, JSON stamp unchanged; at `both` on empty roots it fails naming them. Controls: (N-sense) separation 0.268866 > bound 2.146e-10; (N-drive) (Cust-Z) C32 vs 48 fails; (N-digit) `2.538800824240e-03 → …241e-03` fails exactly the (Cust-P) P_gel string check; (N-blocks) I85 W0b without `at_line` raises. Entry 67's readout is now one `FEM_EM_GEO36_I88=both` window once entries 60 / 61 / 63 / 67 land; `ENTRY67_GLOB` stays provisional. (K-q) `git diff --stat HEAD` before the commit: `PROJECT_PLAN.md 4 +-`, `docs/planning/chunks/GEO-36.md 2 +`, the four `GEO-36-I88-*` logs (W0 807, W0r 875, W1 603, W1-n1 323 +), `docs/testing/test-results.md 4 +`, `tests/validation/test_geo36_i88_rod_heat_per_acc_reader.py 689 +` — 8 files, 3305 insertions, 2 deletions; none of the three queued modules, C45's, C46's, C6a's, `src/` or `docs/testing/xl-queue.d/`.

### §9 item 968 — 2026-11-12 (§10 I87′, chain I — rule (4), the 11-11 weekly's step; chain step (no solve))

*Enumerated by the 11-11 weekly: the ruling `PROJECT_PLAN.md:4919–4945` ("I87 → I87′") and the step `:5022–5030` (line numbers at `28372b5b`; the 11-12 03:00 review's uncommitted §9 edit shifts both by +8). The 2026-11-12 03:00 daily review wrote this brief as execution detail. Gates, bands, controls and the validation target are the weekly's; none is set or moved. It supersedes item 933's brief (above) where they differ; item 933's premises 1–3 are the weekly's accepted premises. **Family I87, attempt 1 of 4:** item 933 was written but never queued or run (its review returned it to the weekly before any slot took it), so it consumed no attempt; the weekly's "re-pose" restates the gates. The review may rule otherwise. **Executor:** the implementer.*

*It **depends on `xl` entry 63 having run**, because of the module hold (`PROJECT_PLAN.md:3970–3974`, the 11-09 ruling (c)). The slot confirms both of these, or it skips the item at no cost:*
- *`git log --oneline` shows `chore(xl): scheduled xl window — GEO-36-step10-I32-prime 2026-11-13` on `main`;*
- *`docs/testing/xl-queue.d/20261113-GEO-36-step10-I32-prime.env` is gone.*

*Entry 63 need not be green: this item reads none of its numbers. Earliest slot: **11-13 04:30**. Entry 55's window took 1015 s (`xl-ledger.md` 2026-11-06 row), so entry 63 is predicted to commit by ≈ 02:20.*

**Facts verified against the tree (11-12).**
- **The validation target is all there.** Item 932's (T-rep) tokens are at `20261108T181211Z_GEO-36-I32b-L-W0b.log:103–130`, every one as item 932 lists it. They are already a module constant: `T_REP` (22 heads), with its matcher `_find(lines, head, tokens)`, in `tests/validation/test_geo36_rod_f_clinical_3t_heating.py:48–77` and `:324–326`. Import both inside the test body; never copy them.
- **The default rod arm is the exclusion.** It reads `rod wall 144 (adiabatic)` and V_sub `9.895833333333e-04` (`L-W0b.log:117`); the twin reads `rod wall 0` and `9.999999999999e-04` (`:106`). Thermal domain: `_tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)` (`test_geo36_rod_f_clinical_heating.py:312`).
- **The stub.** The rod block is the 108 cells with midpoints at \|x\|∞ < 12 mm (`test_geo36_rod_f_clinical_drive.py:1194`), with census `{3: 10260, 403: 108}` (`L-W0b.log:104`). The tip balls sit at (±0.03, ∓0.03, 0) (`drive.py:1201`), radius 6.2035 mm (`L-W0b.log:121`). The ball nearest the block's centre is ≥ 36 mm away, and the block reaches ≈ 24 mm (`heating.py:632`). So V(B ∩ rod) = 0 holds on the stub.
- **The solver takes DG0 coefficients.** `solve_transient_conduction(mesh, facet_tags, *, k, rho_c, q, h, t_inf, t_init, robin_tag, dt, n_steps, theta, startup_half_steps, degree, snapshot_steps, …)` accepts a `fem.Function` for `k` and `rho_c`. The same `rho_c` weights the matrix, the right-hand side and the balance's storage term (`src/fem_em_solver/thermal/transient.py:94–123`). `_transient` hard-codes `k=K_GEL, rho_c=RHO_C` (`test_mat15_coil_driven_heating_phantom.py:145–151`).
- **The source is zero on rod cells.** `_heat_source` averages on the tag it is given (`heating.py:212–216`), and `_cell_average_dg0_tag` leaves untagged cells at zero (`mat15:139–141`). `_tag_submesh`'s `parent_values` pattern is at `mat15:111–113`.
- **The values.**
  - `ROD_K` 6.7 and `ROD_RHO_C` 2.33018e6, both nominal, from memory (`test_imp8_composite_cylinder_closed_form.py:52–53`);
  - `EPS_HOLE` 1e-6 (`test_imp8_rod_conducting_body_fem.py:84`);
  - I66's coefficient rule `_coefficients(msh, rod_idx, arm)` (`:149–167`): gel everywhere; `ti` sets the rod cells to `ROD_K` / `ROD_RHO_C`; `hole` sets them to `EPS_HOLE`·`K_GEL` / `EPS_HOLE`·`RHO_C`.

  The I66 and heating modules' `RHO_C` are the same object, `test_mat12_transient_sar_heating_sphere.RHO_C` (`mat15:73–77`; `imp8_composite:41`).
- **The import cycle is real.** `test_imp8_rod_line_source_fem.py:64` imports this module, and `test_imp8_rod_conducting_body_fem.py:51` imports that one.
- **Who imports the heating module** (grep, 11-12):
  - P61 `test_geo38_f_clinical_box_heating.py:58`;
  - `test_geo38_f_clinical_3t_box_heating.py:36`;
  - I85 `test_geo36_rod_f_clinical_3t_heating.py:92`, inside a function;
  - `test_mat21_clinical_probe_reader.py`, through P61;
  - `test_imp8_entry48_lumped.py:43`;
  - `test_imp8_rod_line_source_fem.py:64`;
  - `test_ops83_s2tip_form_on_stub.py:32`;
  - `test_mat21_p86_three_pictures_3t.py:41–42`;
  - `tests/package/test_wf15_study_rod_heating.py:46`.

  C46's reader `test_geo38_c39_thermal_reader.py` imports only P61's `RUNNING_WINDOW_S`, `SENSE_TUPLE`, `T_STAR_BOX` and `_probes` (`:42–47`). I79's reader `test_geo36_i76_f2182_table.py` does not import it at all.
- **What a default-unset knob can change in the importers: nothing, by construction.**
  - P61, P61-3T and C46 never call `_run_heat`. They use `_thermal_arm`, `_thermal_failures`, `_probe_theta`, `_box_preconditions`, `_refuse_degree`, `_p37_fixture`, `_merged_q`, `_heat_source`, the `STUB_*` constants and the `i32._transient` intercept (`box_heating.py:140–158`). None of these is edited.
  - The only importer that reaches `_run_heat` is I85's body: its K12-R stub and entry 67's gated test call `heat._run_heat` and `heat._stubs`. The knob is read only inside `_run_heat`, and its print appears only when set, so I85's lines are byte-identical while the knob is unset.
- **Entry 63's command sets no `ROD_THERMAL`.** The queue file sets `FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_SOLVER_PROGRESS=2 FEM_EM_GEO36_S8_HEAT=1 FEM_EM_GEO36_S8_DRIVE_C_PF=193.5996` and runs the heating module at `mpiexec -n 8` (`docs/testing/xl-queue.d/20261113-GEO-36-step10-I32-prime.env`). Item 933's "sets only `S8_HEAT` and the node" omitted the first two. No file under `tests/` or in `xl-pending.md` names `ROD_THERMAL`.
- **No foreign guard refuses the knob.** `_guard_foreign` refuses seven names (`drive.py:68–76`, `:97–100`), and the new name is not among them. That leaves entry 70′'s gated run free to set it. I85's 19-name guard does not list it either (Scope).

**Code. One file, `tests/validation/test_geo36_rod_f_clinical_heating.py`, by addition only: every hunk of `git diff` is `+` lines.**
- **Byte-identical:** every existing function, including `_run_heat`'s existing lines, `_thermal_arm`, `_thermal_failures`, `_e_prime_box`, `_ball_reader`, `_tip_slope_rows`, `_stubs`, and every existing test. No `src/` edit. No other test module is touched.
- **Knob.**
  - `ROD_THERMAL_ENV = "FEM_EM_GEO36_S8_ROD_THERMAL"` and `ROD_THERMAL_ARMS = ("gel", "ti", "hole")`.
  - `_rod_thermal()`: unset or empty ⇒ `None`; a listed arm ⇒ itself; anything else (including `TI`, `1`, `x`) ⇒ `ValueError` naming the variable.
- **`_rod_coefficients(sub, parent_ids, cell_tags, arm)`, new.** It returns `(k, rho_c, rod_ind, n_rod)`, three DG0 Functions on the thermal submesh plus the rod-cell count.
  - Rod cells are the submesh cells whose parent tag is 403. Build `parent_values` as `mat15:111–113` does, then take `parent_values[parent_ids] == BIRDCAGE_IMPLANT_CELL_TAG` over local + ghost sub cells.
  - The values follow I66's rule. `rod_ind` is 1 on rod cells and 0 elsewhere.
  - `n_rod` is allreduced over owned cells only.
  - The imports of `ROD_K`, `ROD_RHO_C` and `EPS_HOLE` sit **inside the function body**, because of the cycle.
- **`_thermal_arm_rod_cells(label, msh, cell_tags, e_sup, sigma_f, G, centre, tips, probes, *, rod_thermal, f2182=None)`, new.** It follows `_thermal_arm`'s order line for line (`heating.py:296–426`). Only these differ:
  - **P** is `mean_sar(…, cell_tags=cell_tags, subdomain_ids=PHANTOM_CELL_TAG)` on the **original** tags, as `_thermal_arm` computes it. So G·P is the default rod arm's.
  - **Domain:** `merged = _merge_rod_into_phantom(msh, cell_tags)`, then `sub, parent_ids = _tag_submesh(msh, merged, PHANTOM_CELL_TAG)`. `_outer_robin_tags(sub, centre)` follows.
  - **Source:** `q = _heat_source(msh, cell_tags, e_sup, sigma_f, G, sub, parent_ids)`, again on the original tags.
  - **Solves:** every `_transient(...)` call becomes `solve_transient_conduction(sub, ft, k=k, rho_c=rho_c, q=q, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN_TAG, dt=…, n_steps=…, theta=…, startup_half_steps=…, degree=1, snapshot_steps=…)`. The other arguments are exactly `_transient`'s, and the call order is `_thermal_arm`'s: F2182 CN, t\* CN, BE Δt₁, BE Δt₁/10. Import the solver inside the body from `fem_em_solver.thermal`. Never call `_transient` itself (P61 intercepts it).
  - **(E′)ρc:** ∫ρc_DG0(θᴺ − θᴺ⁻¹)/Δt + ∮hθ^(n+½), against G·P. The scalar `_e_prime_box` on the same `res` is kept beside it.
  - **(S2)ρc:** \|∫ρc_DG0 θ¹ dV / (Δt₁·G·P) − 1\|. The scalar form `RHO_C`·∫θ¹/(Δt₁·G·P) is kept beside it.
  - **The ball reader's fields** are `_thermal_arm`'s plus `rc_th1 = rho_c*th1`, `rc_th10 = rho_c*th10` and `rod = rod_ind`. V(B ∩ rod) = vol·⟨rod⟩_B is printed per ball.
  - **(S2-tip), the gate, the ρc-weighted form as ruled:** `_tip_slope_rows({n: ⟨ρc θ¹⁰⟩_B}, DT_1_BOX/10, comps, 1.0)`.
    - `comps = _projected_ball_comparands(sub, q, V1, tips, a_1g, "i32b_" + label.replace(" ", "_"))`.
    - **Premise corrected:** `_tip_slope_rows` multiplies a scalar `rho_c` by a ball mean (`heating.py:279–288`), so it cannot take a DG0 Function. The weekly's "with the DG0 ρc" enters as the ball mean of ρc·θ¹⁰ with `rho_c = 1.0`, which is algebraically ⟨ρc θ¹⟩_B/Δt against ⟨Πq⟩_B.
    - The scalar form `_tip_slope_rows({n: ⟨θ¹⁰⟩_B}, …, RHO_C)` is printed beside it.
  - **(N-src), the control:** `q_m = _heat_source(msh, merged, …, sub, parent_ids)`, the source built on the merged tags. Its rel miss \|∫q_m − G·P\|/(G·P) is returned.
  - **The return dict** has `_thermal_arm`'s keys plus `Ep_rc`, `Ep_sc`, `S2_rc`, `S2_sc`, `tip_s2_pi` (the ρc form), `tip_s2_pi_sc`, `v_rod_ball`, `n_rod` and `n_src`.
  - Every print uses `_thermal_arm`'s form with `[{label}]`. The label is `f"rod {rod_thermal}"`, giving the PETSc tag `i32b_rod_ti`.
- **`_rod_cells_failures(arm)`, new**, collected as names. It checks (B), (E), (E′)ρc and (S2)ρc at their bands, `counts["wall"] == 0`, and max `tip_s2_pi` ≤ `S2_BAND`.
- **`_run_heat`: four inserted statements, no line changed.**
  1. After `_box_preconditions()`: `rt = _rod_thermal()`. A bad value raises before the fixture or gmsh. If set, the slot also prints `rod thermal: {rt} (rod cells in the thermal domain, I66's coefficients; printed after the exclusion arm)`. **Nothing is printed when unset.**
  2. Before the `heat_rod = _thermal_arm("rod", …)` line: `rod_fields = (arm_rod["e_sup"], arm_rod["sigma"]) if rt is not None else None`. The existing call `pop`s them, so this keeps them alive.
  3. After `fails += _thermal_failures(heat_twin) + _thermal_failures(heat_rod)`, an `if rt is not None:` block:
     - `heat_rod_cells = _thermal_arm_rod_cells(f"rod {rt}", msh, cell_tags, *rod_fields, heat_twin["G"], centre, tips, probes, rod_thermal=rt, f2182=f2182)`;
     - `fails += _rod_cells_failures(heat_rod_cells)`;
     - one print line: the arm's tip θ(t\*) and F2182 θ(900 s) over the exclusion arm's, minus 1, per tip.
  4. At the return, `out["rod_cells"] = heat_rod_cells`, only when set. The return keys are otherwise unchanged.

  **So the default path runs exactly as now**, and a set knob adds one rod arm after the exclusion arm on the same EM fields. On the clinical mesh that costs about one more rod thermal arm (entry 55's took 7.0 s, `20261106T080008Z_GEO-36-step10-I32.log:12522`). Every knob window therefore carries the exclusion and the in-domain arm side by side.
- **New coding tests (`@complex_only @_coding`):**
  - `test_i87p_k0_knob_parsing`: unset, `""`, the three arms, and refusals.
  - `test_i87p_rod_thermal_arms_on_stubs`: (K-def)(b), (K-arm), (K-gel=), (N-src), (N-ρc), (K-c), and the printed former gates.

  Every run happens inside its own `with monkeypatch.context() as m:` with `lines = _stubs(m, fx, root)` and `_run_heat(comm, m, K9_NODE_PF)`, as in K10-R (`heating.py:880–892`). The order is: knob unset, then `m.setenv(ROD_THERMAL_ENV, arm)` for `gel`, `ti` and `hole`. Collect every verdict, print, then assert.

**Gates (asserted; bands imported, none set). These are the weekly's.**
- **(K-def): entry 63's queued command is unchanged, and the default path is byte-identical.**
  - **(a)** Before the commit, `git diff --stat HEAD` lists no test module except the heating module, and no `src/` file. `docs/testing/xl-queue.d/` is unchanged by this slot. Grep shows no `ROD_THERMAL` in `docs/testing/xl-queue.d/`, `xl-pending.md` entry 63, or entry 63's log `- Command:` line. Quote the stat in the Outcome.
  - **(b)** In the unset run, every `T_REP` token is found by `_find` in the recorded `lines`. It is asserted at `comm.size == OPS83_RECORD_RANKS` (`-n 2`) and printed otherwise (the width ruling, `heating.py:928–941`, ratified by the 11-09 weekly). The residual and timing tokens are excluded, as in item 932. Also: no recorded line starts with `rod thermal`, and `"rod_cells" not in out`.
  - **(c)** Wreg, Wimp and Wimp-c are green (Windows).
- **(K-arm), for each of `gel` / `ti` / `hole`:**
  - `_run_heat` returns, so its collected gates are green. These include the exclusion arm's and the twin's gates as now, and the rod-cells arm's (B), (E), (E′)ρc, (S2)ρc, wall 0 and ρc-form (S2-tip);
  - the rod-cells arm's G·P equals `3.958333333333e+00` by `.12e` string (`L-W0b.log:117`) at `-n 2`. At every width, `float` equality with the same run's exclusion arm;
  - wall facets are 0, and V_sub equals the same run's twin V_sub exactly, and `9.999999999999e-04` by string at `-n 2` (`:106`);
  - `n_rod == 108` (`:104`);
  - V(B ∩ rod) == 0.0 at both stub balls.
- **(K-gel=), gel arm only:** \|Ep_rc − Ep_sc\| ≤ `E_BAND` and \|S2_rc − S2_sc\| ≤ `S2_BAND`. On the stub, the ρc form of (S2-tip) equals the scalar one within `S2_BAND` on every arm, because V(B ∩ rod) = 0; it is predicted at the round-off floor.

**Negative controls (asserted, with the ceiling computed).**
- **(N-src):** on every arm, `n_src` > `B_BAND`.
  - Ceiling = the exact separation. With uniform E and σ on the stub, ∫q_m/(G·P) − 1 = V₄₀₃/V₃ = 108/10260 = **1.052632e-02** = share/(1 − share), from `L-W0b.log:103–104`.
  - Printed beside the measured value, as *predicted*: the form is new.
- **(N-ρc):** on `ti` and `hole`, \|Ep_sc\| > `E_BAND`. This shows the coefficient reached the storage term.
  - Predicted factors (item 933's unlogged arithmetic, printed beside): ti ≈ 0.4425 × 0.010417 × 0.711 ≈ 3.3e-3 and hole ≈ 7.4e-3, using stored / G·P = 2.814/3.958 (`L-W0b.log:127`).
  - Asserted only at the gate's own band (rule (e)).
- **(K-c):** in memory, change one θ digit of `T_REP`'s `[rod] 1 g ball … at a,` token (`L-W0b.log:121`, `2.424561997e+00 → …998e+00`). `_find` on the unset run's `lines` must then report it missing, as I68′'s (K-c) does.

**Printed, not asserted: the three former gates with the weekly's predicted ranges, copied verbatim** (`PROJECT_PLAN.md:4935–4938`). Per tip ball a / b, on θ(t\*) and on the F2182 θ(900 s):
- "(K-gel) \|gel / default − 1\| 4e-3 … 7e-3";
- "(K-Ti) Ti / gel − 1 +1.5e-3 … +3.5e-3";
- "(K-hole) hole / gel − 1 +4e-3 … +7e-3";
- "hole / default − 1 ≤ 1e-4";
- "and the measured ordering line" (`ordering at a: …`, sorted by θ).

"default" is the unset run's `[rod]` arm, which is also each knob run's own exclusion arm. Print both: they must agree by string at `-n 2`, and that agreement is also asserted as part of (K-def)(b). V(B ∩ rod) and the two (S2-tip) forms are printed per arm. A value outside its predicted range is a record for the weekly, not a stop.

**Windows.** Before each window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Each window is one bare foreground harness call with Bash timeout 660000 ms. No knob is set on any command line; the knob is set only by `monkeypatch`.
- **Wimp0 / Wimp0-c (before the edit, baseline):** Wimp's and Wimp-c's commands below on the unedited tree, labels `GEO-36-I87p-Wimp0` / `GEO-36-I87p-Wimp0-c`. A collection error that is already there is a pre-existing condition, not this item's.
- **W0b** (`-n 2`), L-W0b's command (`L-W0b.log:12`) with the label and raw log renamed and `timeout -k 30 180`:
  `scripts/testing/run_and_log.sh GEO-36-I87p-W0b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/validation/test_geo36_rod_f_clinical_heating.py -v -s --tb=short > /workspace/logs/geo36-i87p-w0b-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo36-i87p-w0b-raw.log; cat /workspace/logs/geo36-i87p-w0b-raw.log; exit \$rc'"`
  Expected per rank, by token: 11 passed / 1 skipped (the 9 / 1 prior, `L-W0b.log:171`, plus the two new tests).
- **Wreg** (`-n 2`): L-Wreg's command (`20261108T181335Z_GEO-36-I32b-L-Wreg.log:12`) plus I85's module, the one importer that calls `_run_heat`:
  `scripts/testing/run_and_log.sh GEO-36-I87p-Wreg "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/validation/test_ops83_s2tip_form_on_stub.py tests/validation/test_ops84_s2tip_form_rod_pair.py tests/validation/test_geo38_f_clinical_3t_box_heating.py tests/validation/test_geo38_f_clinical_box_heating.py tests/validation/test_geo36_rod_f_clinical_3t_heating.py -v -s --tb=short > /workspace/logs/geo36-i87p-wreg-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo36-i87p-wreg-raw.log; cat /workspace/logs/geo36-i87p-wreg-raw.log; exit \$rc'"`
  Expected per rank: **23 passed / 10 skipped**. That is L-Wreg's 18 / 9 (`:377`, `:380`) plus I85's 5 / 1 (`20261110T105744Z_GEO-36-I85-W0b.log`; item 932 Outcome). I85's K12-R must still find all 22 `T_REP` lines.
- **Wimp** (real build, `-n 1`, collect-only):
  `scripts/testing/run_and_log.sh GEO-36-I87p-Wimp "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 python3 -m pytest --collect-only -q tests/validation/test_geo36_rod_f_clinical_heating.py tests/validation/test_imp8_rod_line_source_fem.py tests/validation/test_imp8_rod_conducting_body_fem.py tests/validation/test_imp8_entry48_lumped.py tests/validation/test_imp8_probe_cylinder_source_closed_form.py tests/validation/test_geo36_i76_f2182_table.py tests/validation/test_mat21_clinical_probe_reader.py'"`
- **Wimp-c** (complex build, `-n 1`, collect-only): the same, with `source /usr/local/bin/dolfinx-complex-mode && ` after `cd /workspace && `, and `FEM_EM_REQUIRE_COMPLEX=1` beside `PYTHONPATH`. Its list adds `tests/validation/test_ops83_s2tip_form_on_stub.py tests/validation/test_geo38_f_clinical_box_heating.py tests/validation/test_geo38_f_clinical_3t_box_heating.py tests/validation/test_geo38_c39_thermal_reader.py tests/validation/test_geo36_rod_f_clinical_3t_heating.py tests/validation/test_mat21_p86_three_pictures_3t.py tests/package/test_wf15_study_rod_heating.py`. Labels `GEO-36-I87p-Wimp-c`.

  Each Wimp window expects the same collection-error count as its Wimp0 baseline, predicted 0.

**Price.**
- **Measured:**
  - L-W0b took 22 s (`L-W0b.log:179`), and one stub `_run_heat` 7.9 s (`:129`), each thermal arm 3.8 s;
  - I85's W0b ran three `_run_heat`s plus K11 in 50 s;
  - L-Wreg took 37 s (`L-Wreg.log:385`).
- **Predicted:**
  - W0b ≈ 21 s (the existing tests) + 4 × 8 s + 3 × ≈ 5 s (the extra rod arm) ≈ 70–90 s, against 180;
  - Wreg ≈ 65–80 s;
  - each collect-only window ≈ 10–20 s.

**Tier.** Standard, complex, `-n 2`, plus four `-n 1` collect-only windows. **Class:** chain step (no solve). **Cost:** ≈ 5.1 measured (class median) / 45 predicted min. This is about 200 lines of code, mirroring `_thermal_arm`; I85's 481-line module is the comparable, built in one slot.

**Traps already paid for.**
- **The import cycle.** Import `ROD_K`, `ROD_RHO_C` and `EPS_HOLE` inside a function body, never at module top (`heating.py:267–268` is the pattern). `T_REP` and `_find` come from I85's module, also inside the test body: I85's module imports this one.
- **The module must still import in the real build** (I58, I55, I79, I83 via `entry48_lumped`). Wimp is the check, against its Wimp0 baseline.
- **Never call `_transient` in the new arm** (P61's intercept keeps the last θ = ½ call, item 870). The F2182 CN stays before the t\* CN.
- **PETSc prefixes come from the label.** `rod ti` → `i32b_rod_ti`. No parenthesis, colon or hyphen. **Never label the new arm `rod`:** I79's `_once` (`test_geo36_i76_f2182_table.py:108–113`) and I88's `parse_i32_block` require one `[rod]` block per log. Under `rod <arm>` a knob window's log stays readable by both. Sequential `_run_heat` calls reusing the `gel twin` tag in one process are proved fine (I85 W0b, three calls).
- **Rank-locality.** `cell_tags.values`, the rod-cell selection and `n_rod` are rank-local: reduce over owned cells. The ball reader allreduces. The DG0 rod indicator is set on local + ghost cells.
- **A new ball field** goes through `_ball_reader` (`ufl.lt(ufl.real(…))`, pinned `QUADRATURE_DEGREE`, `heating.py:228–246`). Never write a bare `<=` on complex operands.
- **Collect every gate, then assert** (item 668). Count verdicts per rank by token. Run with `-s`.
- **Literal (`.12e` / `.9e`) reproductions are asserted at `-n 2` only.** At `-n 1`, item 932's W0n1 missed last digits on V_sub and on the twin's P_gel / G.
- **The leak guard.** If the hook refuses, move the log unchanged to `logs/leak-held/`. Never `--no-verify`.
- **Never delete a log.** A crash in non-gate code before any gate value exists may be fixed and re-run once at the same `-t`, keeping the first log (the 10-31 ruling).

**Scope.** The knob, its in-domain thermal arm beside the exclusion, and the stub identities. **Not claimed:**
- any clinical θ;
- (K-gel), (K-Ti) and (K-hole) as gates;
- any Ti property (nominal, from memory);
- the Ti clinical window's command or price (the daily's variant licence);
- I85's 3 T body under a set knob. Its foreign guard does not list the knob, so whoever makes entry 67 or 70′ READY names it explicitly.

**Negative result.**
- **(K-def) red, or Wreg / Wimp red against its baseline:** the default path or an importer moved. Revert, stop, mark BLOCKED.
- **(K-arm), (K-gel=) or a control red:** park on `attempt/GEO-36-I87p-*`, mark item 968 BLOCKED in the same commit (rule (d)), and print every value.
- **A ρc-form (S2-tip) red on a stub ball with V(B ∩ rod) = 0** is a code defect, not a comparand question: report it.
- **No `PREDICTED RED:` line.** Every asserted gate holds by construction or by identity on the stub, and the three former gates are printed only.

**Status it can move.**
- §10 I87′ DONE.
- `GEO-36` gains a record line: the `ROD_THERMAL` knob (default unset = the adiabatic exclusion), stub-proved, with the gel / Ti / hole tip θ beside the exclusion's and the ρc-weighted tip form gated on the stub.
- **Entry 70′ becomes priceable** under the daily's variant licence once entry 63 reads green (the weekly's "Landing ⇒").
- This slot leaves `xl-pending.md` untouched and queues nothing.

**Depends on:** `xl` entry 63 having run (its `chore(xl)` commit on `main`, its queue file gone). **Earliest slot: 11-13 04:30.** **For the review (not verified as a hold):** the §9 text says "every module importing it is held" (`PROJECT_PLAN.md:3970–3971`). Read literally, that also holds this module through two other entries, because the module itself imports their window modules:
- entry 62's `test_geo36_rod_f_clinical_drive.py` at its top (`heating.py:60`; `xxl` 11-14);
- entry 65's `test_geo38_f_clinical_matrix` inside `_em_failures` (`:458`; `xl` 11-16).

The precedents read the other direction. Entry 60's module imported this one and held it; C21's module was freed once entry 60 ran (`:3972`). In that direction, no queued window's module (62 / 64 / 65) imports this one, directly or through its importers (grep, 11-12). Under it, the module is free once entry 63 has run. The review rules which direction binds.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 968 — `GEO-36` I87′ — per `docs/planning/chunks/GEO-36.md` "§9 item 968". First confirm `xl` entry 63's `chore(xl)` commit (GEO-36-step10-I32-prime, 2026-11-13) is on `main` and `docs/testing/xl-queue.d/20261113-GEO-36-step10-I32-prime.env` is gone; if not, stop. Run Wimp0 / Wimp0-c before editing. Edit only `tests/validation/test_geo36_rod_f_clinical_heating.py`, by addition only: every diff hunk is `+` lines, every existing function and test byte-identical, nothing printed when the knob is unset; no `src/` edit. The three former (K-…) gates are printed, never asserted. Before each window `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; pytest with -s; container `timeout -k 30`. Never move a band. Never delete a log. You are the executor — do not spawn agents. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-I87p-*`, mark item 968 BLOCKED in the same commit, and report."

**Outcome (04:30 slot, 2026-11-13): DONE — I87′ landed; the `ROD_THERMAL` knob is stub-proved.** The implementer added 419 lines to `tests/validation/test_geo36_rod_f_clinical_heating.py`, all additions (`git diff -U0` shows 0 removed lines). There is no `src/` edit. One deviation from the letter: the fourth insertion is an `if rt is not None: return {...}` placed just before the existing return. The existing return is a literal dict, so it cannot be extended without changing that line. The unset path reaches the unchanged return.
- **Windows** (`-n 2` complex unless noted; `pgrep -c python3` read 0 before each):
  - Wimp0 (real, `-n 1`): 53 collected, 0 errors (`20261113T103040Z_GEO-36-I87p-Wimp0.log`).
  - Wimp0-c (complex, `-n 1`): 95 collected, 0 errors (`20261113T103053Z_GEO-36-I87p-Wimp0-c.log`).
  - **W0b:** 11 passed / 1 skipped per rank, Status 0, 74 s (`20261113T103423Z_GEO-36-I87p-W0b.log:388–396`).
  - **Wreg:** 1 failed / 22 passed / 10 skipped per rank, Status 1, 59 s (`20261113T103546Z_GEO-36-I87p-Wreg.log:526`, `:542`). The red is I85's `test_i85_k1_12_knob_parsing`, `assert 20 == 19`.
  - **Wbase:** that node alone on clean `main` `cebdb8fc` (edit parked) reads the same red, Status 1, 3 s (`20261113T103903Z_GEO-36-I87p-Wbase.log:58`, `:74`). So the red is pre-existing (C52's `7ffc523d` added a 20th foreign name) and not against the baseline. It is filed in known-issues 2026-11-13 and not fixed. In the same Wreg, I85's K12-R still finds all 22 `T_REP` lines (`:478`).
  - Wimp (real, `-n 1`): 55 collected, 0 errors, 4 s (`20261113T103929Z_GEO-36-I87p-Wimp.log`).
  - Wimp-c (complex, `-n 1`): 97 collected, 0 errors, 4 s (`20261113T103937Z_GEO-36-I87p-Wimp-c.log`).

  The two Wimp counts are each the baseline plus the two new tests.
- **(K-def)(a):** `git diff --stat HEAD` showed `test_geo36_rod_f_clinical_heating.py | 419 +` and `test-results.md | 4 +`, with no other test module and no `src/` file. `ROD_THERMAL` occurs 0 times in `xl-queue.d/`, `xl-pending.md` and entry 63's log. **(K-def)(b):** green (`W0b.log:385`):
  - `T_REP` missing none at `-n 2`;
  - no `rod thermal` line, and `rod_cells` not in `out`;
  - every knob run's exclusion arm agrees by string with the unset `[rod]`.
- **(K-arm), on gel / ti / hole:**
  - G·P `3.958333333333e+00`, float-equal to the exclusion arm;
  - V_sub `9.999999999999e-04`, equal to the twin's;
  - wall 0, `n_rod` 108;
  - V(B ∩ rod) 0.0 at a and at b.
- **Rod-cells arm:**
  - (B) 8.942e-14 on every arm;
  - (E) 1.058e-10 / 1.154e-10 / 1.165e-10;
  - (S2)ρc 2.0461e-05 on every arm;
  - ρc-form (S2-tip) a / b: gel 8.1631e-07 / 8.1706e-07, ti 4.4918e-06 / 4.4493e-06, hole 2.1220e-06 / 2.0706e-06, each equal to the scalar form.
- **(K-gel=):** \|Ep_rc − Ep_sc\| is 0.000e+00 and \|S2_rc − S2_sc\| is 2.776e-15.
- **Controls:**
  - (N-src) measured 1.052632e-02 on all three arms, against 1.052632e-02 predicted;
  - (N-ρc) \|Ep_sc\| ti 4.159412e-03 (predicted ≈ 3.3e-3) and hole 9.527800e-03 (predicted ≈ 7.4e-3), both > `E_BAND`;
  - (K-c) the altered token is reported missing.
- **Former gates, printed only (records for the weekly).** At tip a, θ(t\*) reads default 2.424561997e+00, gel 2.404589138e+00, ti 2.412392258e+00, hole 2.424561974e+00. The ratios, against the weekly's predicted ranges:
  - (K-gel) 8.2377e-03 against 4e-3 … 7e-3 — **outside**;
  - (K-Ti) +3.2451e-03 against +1.5e-3 … +3.5e-3 — inside;
  - (K-hole) +8.3061e-03 against +4e-3 … +7e-3 — **outside**;
  - hole/default − 1 is −9.2e-09, against ≤ 1e-4.

  On the F2182 θ(900 s) reads, (K-gel) is 2.0901e-03, (K-Ti) +3.9080e-04 and (K-hole) +2.0945e-03, so all three fall below their ranges. Tip b reads within ≈ 1e-6 of tip a. The ordering is gel < ti < hole < default on both tips and both reads.

  These are stub readings only; no clinical θ is claimed. Entry 70′ stays on entry 63 green, which is the weekly's to decide.

### §9 item 970 — 2026-11-12 (§10 I76, chain I — rule (4), the 11-11 weekly's step; chain step (no solve))

*Enumerated by the 11-11 weekly (`PROJECT_PLAN.md:5017–5021` at `28372b5b`, +8 in the review's working tree). The 2026-11-12 03:00 daily review wrote this brief as execution detail. It is a **record**: no band, tolerance or target is set or moved. Family I76, attempt 1 of 4. **Executor:** the implementer.*

*It is serial on **entry 63 green**, by the 03:00 11-13 review's mechanical read: the after-window licence (`PROJECT_PLAN.md:4856–4868`; C39 / P76 are the precedents). The slot checks both of these, or it skips the item at no cost:*
- *`docs/testing/xl-ledger.md`'s `| 2026-11-13 | GEO-36-step10-I32-prime |` row has been filled by a review: memory, verdict and gate columns are non-empty, unlike the cron skeleton of the 11-12 row, `:82`;*
- *the row reads `Status 0`, `[capture] rc=0`, orphans 0 / 0, and every asserted gate green, including `gates: failures none`.*

*If (S2-tip) is red again, entry 63's rule changes the tip read to printed, and the 11-14 weekly decides the branch, so this item waits. Earliest slot: **11-13 04:30**.*

**What the step is.** The clinical F2182-shaped readout: entry 63's probe θ(900 s), 1 mm beyond each capsule apex, per watt deposited in the rod. It is set beside three rows:
- **I75's bracket arms** at 16 mm (1 W in the rod body, flat end);
- **I81's placement ladder** at 16 / 17 / 18 / 20 mm per arm;
- **I83's finite-radius form Φ_a** at the same rungs.

It runs through **I79's reader**, whose one knob is the log name.

**Facts verified against the tree (11-12), with premises corrected.**
- **I79's reader** is `tests/validation/test_geo36_i76_f2182_table.py` (item 893). It does not import the heating module, so no queued window holds it.
  - **Its knob is the module constant `ENTRY63_LOG = None`** (`:45`; "never globbed", `:7`). It is not an environment variable.
  - `ENTRY63_TIPS = ("top", "bottom")` (`:46`). `test_printed_table` (`:252–269`) prints the column as `BLANK` until the knob is set, and then calls `_clinical_cells(_read_log(ENTRY63_LOG), ENTRY63_TIPS, …)` (`:236–249`).
  - `_clinical_cells` asserts that the four F2182 rows and `share` / `twin_g` / `rod_g` each match exactly once (`_once`, `:108–113`). It prints θ / (G_twin·P_rod) "per rod watt".
- **Premise corrected: I79's reader carries I75 and I78 only, not I81 or I83.** Its docstring says "I76 may add them if they land" (`GEO-36.md:1969`). Both have landed:
  - **I81** is `tests/validation/test_imp8_rod_probe_ladder_fem.py`. Its reads are already held as string-custody literals in I83's module: `FEM[arm]` with `mean` / `dmean` / `plus` / `minus` at 16 / 17 / 18 / 20 mm, logs `20261109T06*_IMP-8-item895-W{1,2,3}.log` (`test_imp8_probe_cylinder_source_closed_form.py:57–82`). They are asserted against the logs by `test_k_reads_from_logs` (`:190–205`).
  - **I83** is the same module: `_phi_cyl` (`:106`), and `_closed(z)` → (SCALE·Φ_a at 900 s, D = 0.5 `_duhamel`) (`:208–209`).
  - `RUNGS` / `LADDER_MM` = 16, 17, 18, 20 mm (`test_imp8_probe_ladder_closed_form.py:52–53`), and `Z_C` = L/2 + 2 mm = 17 mm (`test_imp8_probe_f2182_closed_form.py:61`). So the clinical probe sits at `RUNGS[1]`.
  - So "one knob" is the entry-63 column. The I81 / I83 rows are an addition, which I79's own brief licensed.
- **Premise checked and found false: entry 55's log is not this reader's shape.**
  - `20261106T080008Z_GEO-36-step10-I32.log` has **no** F2182 line (grep count 0): entry 55 ran I32a's code, and the F2182 print arrived with I32b (item 890).
  - Its `share` (`:12497`), `twin_g` (`:12499`) and `rod_g` (`:12510`) lines do match I68′'s `^`-anchored patterns.
  - So `_clinical_cells` would fail on entry 55 with four "0 matches". It reads entry 63 only because entry 63 runs the I32b module, whose print form I79's (K-s) proved on `L-W0b.log:108–109`, `:119–120`.
- **Clean lines in entry 55's family.** All of entry 55's `[GEO-36 I32]` lines start at column 0 (`_say` is rank 0 only, `heating.py:116–118`; `start:` begins with `\n`, `:512`). So the anchored patterns are expected to hold on entry 63's `-n 8` log too.
- **The clinical fixture has apexes.** `_probe_points` carries `top_apex` / `bottom_apex` at ±(L/2 + `ROD_RADIUS`) (`drive.py:288–289`), so entry 63 prints four F2182 lines.
- **The log's name is not yet known.** The queue file's `XL_CHUNK` is `GEO-36-step10-I32-prime`, so entry 63's log is predicted at `docs/testing/logs/20261113T0800??Z_GEO-36-step10-I32-prime.log` (I88's `ENTRY63_GLOB`). Take the real basename from the ledger row, never by glob in code.
- **No JSON.** The 11-11 text drops the older "JSON written". The tracked log is the source (§9 ruling, `PROJECT_PLAN.md:3977`), and none is written.

**Code. One file, `tests/validation/test_geo36_i76_f2182_table.py`; every existing function and test stays byte-identical.**
- **The knob:** `ENTRY63_LOG = "<the ledger row's basename>"`. This is the only changed line.
- **Literals, transcribed from entry 63's log** with their line numbers in comments, in `STUB_*`'s form:
  - `ENTRY63_SHARE` (4 tokens), `ENTRY63_TWIN_G` (3), `ENTRY63_ROD_GP` (1);
  - `ENTRY63_F2182[(arm, tip)]` = (point, `"1.0"`, `"900"`, θ, `"True"`) for `gel twin` / `rod` × `top` / `bottom`.
- **`_entry63_failures(text)`, new:** `_stub_failures`' body over the `ENTRY63_*` literals and `ENTRY63_TIPS`. `_stub_failures` itself is untouched.
- **New tests:**
  - **`test_i76_k63_entry63_reads_as_recorded`: (K-63).**
  - **`test_i76_kc63_one_digit_changed_fails`: (K-c63).**
  - **`test_i76_p63_placement`: (P63) and (L17).**
  - **`test_i76_readout_beside_ladder_and_phi_a`: the printed record, plus (N-row).** Import I83's module **inside the test body** (`import tests.validation.test_imp8_probe_cylinder_source_closed_form as i83`). I83's module imports `test_imp8_entry48_lumped`, which imports the heating module at its top (`entry48_lumped.py:43`). A module-level import would put I79's reader into the heating module's import closure, and hold it whenever a heating window is queued (entry 70′).

**Gates (asserted; no band; string custody and identities).**
- **(K-63):** on entry 63's log, each of the four F2182 rows and `share` / `twin_g` / `rod_g` matches exactly once, and every token string-equals its `ENTRY63_*` literal. This is I79's (K-s) on the real log.
- **(P63), the placement identity on the real log.**
  - The parsed top and bottom point strings share their x and y tokens.
  - \|(z_top − z_bottom)/2 − `Z_C`\| ≤ 1e-6 m. The bound is two half-units of the print's 6-decimal rounding (`heating.py:336`).
  - The offset token is `"1.0"`.
- **(L17):** \|`i83`'s `RUNGS[1]` − `Z_C`\| ≤ 1e-15 m (I79's `P_BOUND`). I81's 17 mm row is the clinical row by constants, never by a literal 17.
- **Ride-alongs:** I83's own `test_k_reads_from_logs` (I81 custody) and I78's module, in the same window.

**Negative controls (asserted).**
- **(K-c63):** in memory, change one θ digit of entry 63's `[rod] F2182-shaped probe top` line. `_entry63_failures` must then be non-empty (I79's (K-c) form).
- **(N-row):** I81's gel 16 mm mean placed in the 17 mm cell misses I81's 17 mm mean by more than one unit of its 7th digit.
  - The separation is 1 − 6.662160/9.161894 = 0.272840, from two logged strings (`i83.FEM["gel"]["mean"]`, the W2 log `:124` / `:127`).
  - It is unbounded above (a difference). This is the 1 mm confound that I79's brief named.

**Printed (≤ 7 digits; a record, no band; each row states its source and normalisation).**
- **The clinical row** (`_clinical_cells`' arithmetic): rod and gel twin at top and bottom, θ(900 s) in K, and per rod watt θ / (G_twin·P_rod). Label: "gel SAR, rod unheated, ASTM box; not a matched source".
- **Beside it:**
  - I75's means at 16 mm (`I75`);
  - I81's ± mean and D = 0.5 mean per arm at 16 / 17 / 18 / 20 mm (`i83.FEM`), labelled "1 W in the rod, I55's cube";
  - I83's SCALE·Φ_a and D = 0.5 at the same rungs (`i83._closed(z)`), labelled "uniform cylinder, infinite gel, 1 W";
  - I78's Z_C row (`I78`).
- **Ratios at 17 mm:** clinical rod per rod watt / I81's ti / gel / hole mean, and / SCALE·Φ_a. The rod / twin tip ratio is also printed. The existing `test_printed_table` now prints the entry-63 column in place of `BLANK`.

**Window.** Before it, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. It is one bare foreground harness call, Bash timeout 660000 ms, real build, `-n 1` (the weekly's):

`scripts/testing/run_and_log.sh GEO-36-I76-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/validation/test_geo36_i76_f2182_table.py tests/validation/test_imp8_probe_f2182_closed_form.py tests/validation/test_imp8_probe_cylinder_source_closed_form.py -v -s --tb=short'"`

Expected: **24 passed**: I79's 8 existing + 4 new, I78's 5, I83's 7. The priors are I79 W1 `20261109T032235Z_GEO-36-I79-W1.log` (13 per rank with I78) and I83 W1 `20261109T110006Z_IMP-8-item911-W1.log`.

**Price.** I79's W1 took 3 s and I83's W1 4 s (`…item911-W1.log:218`). Parsing a ≈ 12 700-line `xl` log several times adds a few seconds. Predicted ≈ 10–20 s. **Tier.** Standard (the weekly's), real, `-n 1`. **Class:** chain step (no solve). **Cost:** ≈ 5.1 measured (class median) / 25 predicted min.

**Traps already paid for.**
- **Never glob for the log in code; never read `logs/geo36/*.json`** (host-only). Read through `_read_log`: `xl` logs gzip at 7 days, and a record is never gunzipped.
- **Never import a `test_*` name.** I83's module objects (`FEM`, `RUNGS`, `LADDER_MM`, `_closed`) are imported inside the test body (see Code).
- **Anchor every pattern on `^\[GEO-36 I32\] ` with `re.MULTILINE`.** pytest's own tokens appear once per rank. Entry 59's `:12343` began `PASSED…` (I88's trap), so a 0-match on entry 63 is a finding, not a regex to loosen.
- **Never write a literal from memory.** Copy each token from the log line it cites, and check the line before writing.
- **The leak guard.** The per-rod-watt and ratio figures are new digit strings. If the hook refuses, move the log unchanged to `logs/leak-held/`. Never `--no-verify`. Never open `docs/private/` or any `aed_results/`.
- Run with `-s`.

**Scope.** Entry 63's F2182-shaped read beside the bracket, the ladder and the closed form. **Not claimed:**
- a like-for-like comparison: the clinical FEM heats the gel by SAR with the rod unheated, while I75 / I81 / I83 put 1 W in the rod, so "per rod watt" is a normalisation;
- a band; compliance; a tissue model;
- the F2182 placement or time (quoted from memory, W3);
- the published comparand (7c, W3);
- I88's per-accepted-watt figure (its own `FEM_EM_GEO36_I88` read).

**Negative result.**
- **(K-63) red** (a pattern matching 0 or 2 times, or a token mismatch): the reader does not read entry 63. Name the line, stop, mark BLOCKED (rule (d)). Never edit a pattern or a literal to match.
- **(P63) or (L17) red:** the clinical probe is not at I78's / I81's 17 mm. Print both placements, stop, mark BLOCKED.
- **A ride-along red:** a landed record moved. Stop, mark BLOCKED.
- **No `PREDICTED RED:` line.**

**Status it can move.**
- §10 I76 DONE.
- `GEO-36` gains a record line: the clinical rod's F2182-shaped θ(900 s) per rod watt at 17 mm, beside I75's bracket, I81's ladder and I83's Φ_a. It is a record, not a matched comparison.
- No `xl` file is touched.

**Depends on:** `xl` entry 63 green by the 03:00 11-13 review's mechanical read, with the ledger row filled. **Earliest slot: 11-13 04:30.** It is independent of item 968: this file is not the heating module, so the two can share a slot in either order.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 970 — `GEO-36` I76 — per `docs/planning/chunks/GEO-36.md` "§9 item 970". First confirm `docs/testing/xl-ledger.md`'s 2026-11-13 `GEO-36-step10-I32-prime` row is filled by a review and reads Status 0, `[capture] rc=0`, orphans 0 / 0 and every asserted gate green; if not, stop. Edit only `tests/validation/test_geo36_i76_f2182_table.py`: set `ENTRY63_LOG` to the ledger row's basename, add the `ENTRY63_*` literals copied from that log with their line numbers, and the four new tests; every existing function and test byte-identical; import I83's module inside the test body. No `src/` edit, no glob, no JSON. Before the window `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; pytest with -s; container `timeout -k 30`. Never edit a pattern or literal to rescue a red. You are the executor — do not spawn agents. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-I76-*`, mark item 970 BLOCKED in the same commit, and report."

**Review note (2026-11-13 03:00 daily) — item 970 BLOCKED before any slot took it; `xl` entry 63 read RED.** The mechanical read (ledger row 2026-11-13, from a `log-pathologist` reading, COUNTABLE) finds every EM gate, (B), (E′), the global (S2) and (S0) green, and three thermal reds in `test_geo36_i32_rod_heating` (`20261113T080008Z_GEO-36-step10-I32-prime.log:12532`): (E) 9.862e-09 / 9.220e-09 on both arms — the residual of I32b's added BE step at Δt₁/10, folded into `out["E"]` (module `:354`, `:418–420`), the absolute-T floor by class — and the rod's (S2-tip) R_Π at Δt₁/10, 1.0920e-04 / 1.1039e-04, falling ≈ 10× from Δt₁ (the O(Δt) step term; a zero Δt-independent part not shown on two points). θ(t\*), the probes and the rod share equal entry 55's at every printed digit. The decision rule (`xl-pending.md` entry 63) has no branch for an (E) red beside an (S2-tip) red, so **the 11-14 weekly reads the branch**; this item's own premise ("every asserted gate green") is unmet, and it stays BLOCKED — attempt count unchanged, it never ran. The same holds for U153, 70′ and `xl` entry 67. Known-issues 2026-11-06 gains "Entry 63's half". Journal: `docs/planning/reviews/2026-11-13-daily.md` §4–§5.

**Review note (2026-11-13 10:30 daily) — item 968 (I87′) accepted as DONE; its branch deleted; the I85 count literal ruled.**
- **968 stands DONE.** W0b green at `-n 2` (`20261113T103423Z_GEO-36-I87p-W0b.log`); the one Wreg red (`20261113T103546Z_GEO-36-I87p-Wreg.log`) is reproduced on clean `main` `cebdb8fc` with the edit absent (`20261113T103903Z_GEO-36-I87p-Wbase.log:58`, `:74`, `assert 20 == 19`), so it is not 968's. The former gates printed only (K-gel / K-hole above the weekly's predicted ranges at tip a, K-Ti inside) are records for the 11-14 weekly, not reds: the brief demoted them to prints.
- **`attempt/GEO-36-I87p-20261113T1038Z` deleted.** Its one commit `c834ad75` (a pre-baseline WIP) carries the test module byte-identical to `ca90290f` (`git diff attempt/… ca90290f` touches only plan, records and the three post-baseline logs that `ca90290f` adds), so nothing on it is uncaptured.
- **Ruling on known-issues 2026-11-13 (I85's foreign-guard count).** The literal `19` at `test_geo36_rod_f_clinical_3t_heating.py:206` is a **record** of how many names the guard refuses, and C52 (`7ffc523d`) changed that number by construction when it added `FEM_EM_GEO38_C21_DEGREE` to C32's `FOREIGN_ENVS_C32_EXTRA` (`test_geo36_rod_f_clinical_3t_drive.py:62`). It is **re-recorded at 20 in its own step (§9 item 978)**, with the new member named in an assertion — **not** derived from C32's list, because a count derived from the list it counts asserts nothing. This is a record re-record, not a loosened assertion: the uniqueness check and the four named members stay.

### §9 item 978 — 2026-11-13 (lane 1 — a known-issues entry's own "Resolves with"; `GEO-36` I85; revalidation)

*Written by the 2026-11-13 10:30 daily review under §9 rule (5), lane 1 (open known-issues entries whose "Resolves with" a slot can execute). Not a chain step. Attempt 1 of 4. Depends on `xxl` entry 62 having run (the module hold, below).*

**Anchor.** A count identity on the foreign guard of I85's module: `len(foreign) == len(set(foreign)) == 20`, where `foreign = _foreign_envs()` (`:106–108`) is C32's `_foreign_envs()` plus `HEAT_ENV` and `S11_DRIVE_ENV`, and `"FEM_EM_GEO38_C21_DEGREE" in foreign` (the member C52 added). Every name in `foreign` is still refused by `_knobs()` with `FOREIGN_MSG` (the node's existing loop, unedited).

**Negative control (already logged, pinned to the pre-change commit).** The same node at `cebdb8fc` with the literal 19 reads `assert 20 == 19` on both ranks (`20261113T103903Z_GEO-36-I87p-Wbase.log:58`, `:74`). No re-run: that log is the control, cited by its full basename. **Reading (by string, not asserted):** the K1-12 print in W1 lists the guard's names; set against the 19 names printed in `20261110T105744Z_GEO-36-I85-W0b.log`'s K1-12 line, the difference is exactly `{FEM_EM_GEO38_C21_DEGREE}` — the executor reads the two printed tuples and records the difference in the Outcome.

**Code.** `tests/validation/test_geo36_rod_f_clinical_3t_heating.py` only: at `:206` the literal `19` → `20`, and one new line after it, `assert "FEM_EM_GEO38_C21_DEGREE" in foreign, foreign`. No other line of the module, no other module, no `src/`.

**Module hold (the 11-12 03:00 ruling, provisional for the 11-14 weekly).** I85's module imports `tests.validation.test_geo36_rod_f_clinical_drive` (`:91`, as `i36`), which is `xxl` entry 62's own module (`docs/testing/xxl-queue.d/20261114-GEO-38-C19.env`), so I85's module is an importer of a queued window's module and is held until 62 runs. It imports neither entry 64's (`test_geo38_f_clinical_3t_esr_root.py`) nor entry 65's (`test_geo38_f_clinical_matrix.py`) module, and neither window imports it. **The slot confirms entry 62's `chore(xxl)` commit is on `main` and `20261114-GEO-38-C19.env` is gone, or skips the item at no cost.**

**Windows.** Bare and foreground, Bash timeout 660000 ms.
- **W1** (the node; prior 3 s, Wbase):
  `scripts/testing/run_and_log.sh GEO-36-I85c-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/validation/test_geo36_rod_f_clinical_3t_heating.py::test_i85_k1_12_knob_parsing -v -s --tb=short'"`
- **W2** (the whole module; Wreg ran it among five modules in 59 s):
  `scripts/testing/run_and_log.sh GEO-36-I85c-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo36_rod_f_clinical_3t_heating.py -v -s --tb=short'"`
  Expect 0 failed per rank by token; the knob-gated body skips by design.

**Tier.** Standard, complex, `-n 2`, no mesh. **Class:** lane-1 revalidation. **Cost:** ≈ 3.5 measured / 10 predicted min.

**Traps already paid for.** Count verdicts per rank by token, never by adjacency; `-s` or the K1-12 print is captured; `pytest -k a or b` splits inside the quoted container command — name the node by `::`; do not "fix" the count by deriving it from C32's list (the ruling above); never `--no-verify`; record rows appended by hand.

**Scope.** One count literal re-recorded and one membership assertion added. **Not claimed:** anything about entry 67's knob-gated body (it skips), I87′, or any band.

**Negative result.** W1 red with a count other than 20 ⇒ a guard name was added or lost since `cebdb8fc` — report the printed tuple against the 11-10 log's, stop, leave known-issues 2026-11-13 open with the new count. W2 red on any other node ⇒ a new finding: file it in known-issues, stop; the literal edit still lands only if W1 is green and the other red reproduces with the edit reverted. No `PREDICTED RED:` line.

**Status it can move.** Known-issues 2026-11-13 (I85's foreign-guard count) **retired** in the same commit; §9's "`main` has no known red test node in an ordinary window" is true again for I85's module; `GEO-36` (✅) records step I85c.

**Depends on:** `xxl` entry 62 having run (Saturday 11-14).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 978 — `GEO-36` I85c, the I85 foreign-guard count re-record — per `docs/planning/chunks/GEO-36.md` "§9 item 978". First confirm `xxl` entry 62's `chore(xxl)` commit is on `main` and `docs/testing/xxl-queue.d/20261114-GEO-38-C19.env` is gone; otherwise stop at no cost. Edit only `tests/validation/test_geo36_rod_f_clinical_3t_heating.py` `:206` (19 → 20) and add the one membership assertion after it. Run W1 then W2. Retire known-issues 2026-11-13 in the same commit if both are green. Move no band. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`. Do not spawn agents."

**Outcome (07:30 CST slot, 2026-11-14) — DONE, in-slot (not delegated).** Pre-check met: `xxl` entry 62's `chore(xxl)` commit `3185acfd` is on `main` and `docs/testing/xxl-queue.d/20261114-GEO-38-C19.env` is gone; `pgrep -c python3` read 0. Code: `tests/validation/test_geo36_rod_f_clinical_3t_heating.py` `:206` `19` → `20`, plus `:207` `assert "FEM_EM_GEO38_C21_DEGREE" in foreign, foreign`; nothing else. **W1** `20261114T133050Z_GEO-36-I85c-W1.log`: 1 passed per rank (`:57`, `:60`), Status 0, 5 s. **W2** `20261114T133105Z_GEO-36-I85c-W2.log`: 9 passed / 1 skipped per rank (`:203`, `:206`), 0 failed, Status 0, 28 s. **Reading (by string):** W1's K1-12 tuple (`:52`) has 20 names, and the 11-10 tuple (`20261110T105744Z_GEO-36-I85-W0b.log:51`) has 19. The difference is exactly `{FEM_EM_GEO38_C21_DEGREE}`, which sits between `FEM_EM_GEO38_C21` and `FEM_EM_GEO38_C20`. Control: the logged `assert 20 == 19` (`20261113T103903Z_GEO-36-I87p-Wbase.log:58`, `:74`), not re-run. Known-issues 2026-11-13 retired in the same commit. No band moved.

## §7 row moved verbatim 2026-11-14 (rotate_plan_archive.py chunks)

 A rod in the phantom: the implant-in-phantom generator (Phase 7 (7b)) | ✅ Audit PASS 2026-10-02 (CAD identities at `-n 2`); the kept rod on F-small, F-human and in F-clinical's ASTM box is a family of records in the chunk file — entry 48 (I36, read first branch by the 11-02 weekly), I57 and I67 (the rod-model lever ≈ 3–5 % on one definition), C32 (I36-3T's body, entry 61's (i) MET) and entry 55 (I32, red on (S2-tip) alone with the EM gates green; read by the 11-07 weekly: the EM record stands and the tip gate is re-posed on the consistent-mass projection — I32b, §9 item 870, then `xl` entry 63 (I32′, 11-13); I68′ is item 872). **I32b landed 2026-11-08 (item 870 via 890) as a stub proof:** the tip gate reads R_Π against ⟨Πq⟩_B at Δt₁/10 at the unmoved `S2_BAND`. It reproduces `OPS-83`'s checkerboard R_Π 6.8552e-07 / 8.5219e-07 to ≤ 4.4e-12 at `-n 2` (literals width-scoped by the 10:30 ruling, provisional), and the retired DG0 form misses at 15.79 × band (`20261108T181211Z_GEO-36-I32b-L-W0b.log:159–162`). The shared helper `_projected_ball_comparands` / `_tip_slope_rows` is on `main` for U127′. **I85 landed 2026-11-10 (item 932) as a stub proof:** the I32-3T body (`test_geo36_rod_f_clinical_3t_heating.py`, knob `FEM_EM_GEO36_S12_HEAT`) runs I32b's `_run_heat` under C21's 3 T patch set on C32's 3 T stub; s2 equals C32's 128 MHz route to 0.000e+00 (64 MHz route misses by 2.515e-01) and I32b's 22 stub-literal lines reproduce by string at `-n 2` (`20261110T105744Z_GEO-36-I85-W0b.log:155–158`); no 128 MHz θ claimed (entry 67's). **I87′ landed 2026-11-13 (item 968) as a stub proof:** the knob `FEM_EM_GEO36_S8_ROD_THERMAL` (default unset = the adiabatic exclusion, byte-identical) adds an in-domain rod arm (gel / Ti / hole coefficients) beside the exclusion, with the ρc-weighted tip form gated on the stub (`20261113T103423Z_GEO-36-I87p-W0b.log`). **I85c landed 2026-11-14 (item 978):** I85's foreign-guard count re-recorded 19 → 20 with `FEM_EM_GEO38_C21_DEGREE` asserted as a member (`20261114T133105Z_GEO-36-I85c-W2.log`, module 0 failed). *History: `docs/planning/chunks/GEO-36.md`.* 

### §9 item 987 — 2026-11-15 (§10 I76′, chain I — rule (4), the 11-14 weekly's re-pose of item 970; chain step (no solve))

*Enumerated by the 11-14 weekly (`PROJECT_PLAN.md:5153–5164`, landed in `b23fb062`). It re-poses §9 item 970, which was BLOCKED by the 03:00 11-13 review and never ran. The 2026-11-15 review wrote this brief as execution detail. It reuses item 970's brief (`GEO-36.md` "§9 item 970", `:2654–2762`) and changes only what the re-pose changes: the pre-check, the table header, and the ID. It is a **record**: no band, tolerance or target is set or moved. Family I76, **attempt 1 of 4**, because 970 never ran. **Executor:** the implementer.*

*It is serial on nothing (the weekly). The slot checks both of these, or it skips the item at no cost:*
- *`docs/testing/xl-ledger.md`'s `| 2026-11-13 | GEO-36-step10-I32-prime |` row (`:83`) has been filled by a review: memory, verdict and gate columns are non-empty. Verified on 11-15: it is filled, reads RED, and was filled by the 03:00 11-13 review from a `log-pathologist` reading.*
- *The 11-14 weekly's branch ruling is on `main`: `PROJECT_PLAN.md:4857` "Rulings, 2026-11-14", summarised at `:3975` ("I76 is re-posed onto entry 63's log as I76′ (987)"). It landed in `b23fb062`.*

*This replaces 970's "Status 0 … every asserted gate green" pre-check. Entry 63 is RED by record. That red is the window's status and goes in the table header. It is not a stop condition.*

**What the step is.** It is item 970's step, unchanged in substance. The clinical F2182-shaped readout is entry 63's probe θ(900 s), 1 mm beyond each capsule apex, per watt deposited in the rod, on both arms. It is set beside three rows:
- **I75's bracket arms** at 16 mm (1 W in the rod body, flat end);
- **I81's placement ladder** at 16 / 17 / 18 / 20 mm per arm;
- **I83's finite-radius form Φ_a** at the same rungs.

It runs through **I79's reader**, with `ENTRY63_LOG = "20261113T080008Z_GEO-36-step10-I32-prime.log"`.

**Anchor (the validation target).** Two parts:
- the log's printed strings, reproduced by string custody: (K-63) on the real log;
- I83's Φ_a at its imported band: I83's own module rides along green in the same window.

No new number is gated.

**Facts verified against the tree (2026-11-15, `main` at `fb9c0ea6`).**
- **The log exists:** `docs/testing/logs/20261113T080008Z_GEO-36-step10-I32-prime.log`, 816 955 B, not yet gzipped. Read it through `_read_log` anyway, because it rotates at 7 days (≈ 11-20).
- **The reader exists:** `tests/validation/test_geo36_i76_f2182_table.py`, 269 lines, last touched by `5be8e5e6` (I79, item 893).
  - `ENTRY63_LOG = None` at `:45`; `ENTRY63_TIPS = ("top", "bottom")` at `:46`; `BLANK` at `:47`.
  - `_read_log` `:99`, `_once` `:108`, `_clinical_rows` `:116`, `_stub_failures` `:130`, `_clinical_cells` `:236`.
  - `test_printed_table` `:252`; its `ENTRY63_LOG is None` branch is at `:265–268`.
  - It has 8 `test_` functions. `P_BOUND = 1.0e-15` at `:90`.
- **Nothing imports the reader.** A grep for `i76_f2182_table` over `tests/`, `src/` and `scripts/` finds no importer.
  - So the one edited file is outside every queued window's import closure. The queued windows are `xl-queue.d/20261116-GEO-38-C35w.env` (entry 65, `test_geo38_f_clinical_matrix.py`) and `xxl-queue.d/20261128-GEO-38-P33-d2.env` (entry 69, `test_geo38_f_clinical_box_heating.py`).
  - The reader's top-level imports (`:35–40`) touch neither module.
- **The entry 63 lines the reader parses all match their `^\[GEO-36 I32\] ` anchors (`_H`, `test_geo36_i68_lumped_vs_clinical.py:85`).**
  - `share` (`PATTERNS["share"]`, `i68:89`): line `:12497`.
  - `twin_g` (`i68:92`): line `:12500`.
  - `rod_g` (`i68:97`): line `:12515`.
  - F2182 rows: gel twin top `:12502`, bottom `:12503`; rod top `:12517`, bottom `:12518`. There are four lines, plus the header `:12499`, which the `F2182` pattern (`:92–93`) does not match because it has no `[arm]`.
  - The point tokens are `0.105, 0.0, 0.017` (top) and `0.105, 0.0, -0.017` (bottom); offset `1.0`, t `900`, valid `True`.
  - Copy the θ tokens from those lines. Never take them from this brief.
- **(P63) is satisfiable on the printed strings.** The two points share x and y. (z_top − z_bottom)/2 = 0.017 m, and `Z_C = 0.5 * L + 2.0 * MM` (`test_imp8_probe_f2182_closed_form.py:61`) is 17 mm.
- **The window's status, for the header (ledger `:83`, log `:12532`).**
  - The single gates line reads `failures ['gel twin (E) 9.862e-09', 'rod (E) 9.220e-09', 'rod (S2-tip) 1.104e-04']`.
  - **The Δt₁/10 (E) residual is never printed alone.** It is inferred from the max: CN / Δt₁ are printed green at `:12514` and `:12529`. So the header cites `:12532` and the ledger, and makes no claim of a separate print.
  - (S2-tip) at Δt₁/10 is `:12521–12522`.
  - EM gates green: ledger `:83`; (M2) / (Z) "failures none" `:12062`.
- **I83's module** is `tests/validation/test_imp8_probe_cylinder_source_closed_form.py`. Its objects: `FEM` `:60`, `_phi_cyl` `:106`, `test_k_reads_from_logs` `:190`, `_closed` `:208`. `RUNGS` / `LADDER_MM` are imported from the ladder module at `:47`.
  - It imports `test_imp8_entry48_lumped` at its top (`:45`), which pulls in the heating module. So I83's import stays **inside the test body**, as in 970.
  - It has 7 tests; I78's module has 5. Expected: 8 + 4 + 5 + 7 = **24 passed**.
- **Item 970's other facts stand unchanged:**
  - I79 carries I75 / I78 only; I81 / I83 are added rows;
  - entry 55's log has no F2182 line;
  - no JSON;
  - the apexes are at ±(L/2 + `ROD_RADIUS`).

**Code. One file, `tests/validation/test_geo36_i76_f2182_table.py`. Item 970's Code section applies verbatim, with two changes:**
- **The knob:** `ENTRY63_LOG = "20261113T080008Z_GEO-36-step10-I32-prime.log"`. Confirm the basename against the ledger row before writing it.
- **The table header.** The entry-63 column's printed header states the window's status:
  - "entry 63 RED as a window (ledger 2026-11-13; `:12532`): EM gates green; (E) at Δt₁/10 and rod (S2-tip) above band as readings (I32c, item 986; the asserted set is re-bought by `xl` entry 71)";
  - and "θ / probe records equal to entry 55's at every printed digit".
  - The header is one `_p` line, printed only, beside the existing `test_printed_table` call. The header string is new; the function itself stays byte-identical except for the `BLANK` branch, which is no longer taken.

Everything else is 970's:
- the `ENTRY63_SHARE` / `_TWIN_G` / `_ROD_GP` / `_F2182` literals, with line comments `:12497`, `:12500`, `:12515`, `:12502–12503`, `:12517–12518`;
- `_entry63_failures`;
- the four tests (`test_i76_k63_entry63_reads_as_recorded`, `test_i76_kc63_one_digit_changed_fails`, `test_i76_p63_placement`, `test_i76_readout_beside_ladder_and_phi_a`);
- every existing function and test byte-identical.

**Gates (asserted; no band; string custody and identities).** These are 970's, unchanged:
- **(K-63):** each pattern matches exactly once and every token string-equals its literal.
- **(P63):** shared x / y; \|(z_top − z_bottom)/2 − `Z_C`\| ≤ 1e-6 m; offset `"1.0"`.
- **(L17):** \|`RUNGS[1]` − `Z_C`\| ≤ 1e-15 m.
- **Ride-alongs:** I83's `test_k_reads_from_logs` and I78's module.

**Negative controls (asserted).**
- **(K-c63):** change one θ digit of the `:12517` rod top line in memory; `_entry63_failures` must then be non-empty.
- **(N-row):** I81's gel 16 mm mean placed in the 17 mm cell misses the 17 mm mean by more than one unit of the 7th digit. The separation is 0.272840, from `i83.FEM["gel"]["mean"]`.

**Printed (≤ 7 digits; a record).** As in 970:
- the clinical row per rod watt;
- I75 at 16 mm;
- I81 ± / D = 0.5 at 16 / 17 / 18 / 20 mm;
- I83's SCALE·Φ_a / D = 0.5;
- I78's Z_C;
- the 17 mm ratios;
- the rod / twin tip ratio.

Add the header line above.

**Tier / ranks / cost.**
- Standard (the weekly's), real build, `-n 1`. **Class:** chain step (no solve).
- Window: predicted ≈ 10–20 s. The priors are I79 W1 3 s and I83 W1 4 s, plus several parses of a ≈ 12 700-line log.
- **Cost:** ≈ 10.8 measured (class median) / **25 predicted min**.

**Window.** Before it, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Run one bare foreground harness call, Bash `timeout: 660000`:

`scripts/testing/run_and_log.sh GEO-36-I76p-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/validation/test_geo36_i76_f2182_table.py tests/validation/test_imp8_probe_f2182_closed_form.py tests/validation/test_imp8_probe_cylinder_source_closed_form.py -v -s --tb=short'"`

Expected: **24 passed**.

**Traps already paid for.** These are 970's:
- No glob and no `logs/geo36/*.json`. Read through `_read_log`, and never gunzip.
- Never import a `test_*` name. Import I83 in the test body.
- Use the `^`-anchored patterns with `re.MULTILINE`. Entry 63's line `:12533` is `FAILEDFAILED…`, so a 0-match is a finding.
- Never write a literal from memory.
- Use the leak guard. On a hook refusal, move the log unchanged to `logs/leak-held/`. Never `--no-verify`, and never open `docs/private/` or any `aed_results/`.
- Run with `-s`.

There are three new traps:
- **The ledger's verdict is RED, and that is not this item's stop.** The pre-check is "filled + ruling on main". A slot that stops on "Status 1" re-creates 970's block.
- **The log's "(E) balance … ASSERTED" text (`:12514`, `:12529`) prints green while the asserted (E) is red.** It is a label trap (ledger `:83`). The header quotes `:12532`, never those lines.
- **If I32c (item 986) lands earlier in the same slot,** the heating module changes, and I83's ride-along imports it through `entry48_lumped`. That is harmless at import. A ride-along red after 986 is named against 986, not fixed here. This item never edits `test_geo36_rod_f_clinical_heating.py` (I32c's, and entry 71's once queued).

**Scope.** Entry 63's F2182-shaped read beside the bracket, the ladder and the closed form, on a window that is red as a whole. **Not claimed:**
- that entry 63 is green;
- a like-for-like comparison ("per rod watt" is a normalisation: SAR-heated gel, rod unheated, against 1 W in the rod);
- a band; compliance; a tissue model;
- the F2182 placement or time (quoted from memory, W3);
- the published comparand (7c);
- I88's per-accepted-watt figure.

**Negative result.** The same as 970:
- (K-63) red: name the line, stop, BLOCKED; never edit a pattern or a literal.
- (P63) or (L17) red: print both placements, then BLOCKED.
- A ride-along red: BLOCKED.

No `PREDICTED RED:` line.

**Status it can move.**
- §10 I76′ DONE. §9 item 970 is superseded (DONE-by-987).
- `GEO-36` gains a record line: the clinical rod's F2182-shaped θ(900 s) per rod watt at 17 mm, from a window red on two one-step thermal readings, beside I75 / I81 / I83. It is a record, not a matched comparison.
- No `xl` file is touched. Entry 71's re-run does not re-open this item: its log is a new basename and I76″ is the weekly's to pose.

**Files.** `tests/validation/test_geo36_i76_f2182_table.py` (the only edit), `docs/testing/logs/<W1 log>`, `PROJECT_PLAN.md` (§9 item 987, §10 I76′ annotation), `docs/planning/chunks/GEO-36.md` (Outcome).

**Depends on:**
- the ledger row filled (met, `:83`);
- the 11-14 weekly's ruling on `main` (met, `b23fb062`).

It is independent of item 986 (a different file) and can share a slot with it in either order. Earliest slot: the next one.

**Outcome commit.** `test(GEO-36 I76′): entry 63's F2182 column through I79's reader, beside I75 / I81 / I83 (§9 item 987)`. It contains the module, the W1 log, the §9 / §10 / chunk-file annotations, and `git diff --stat` quoted in the Outcome (one test file).

**Outcome (12:00 slot, 2026-11-15): DONE, 24 passed.** W1 (`20261115T180232Z_GEO-36-I76p-W1.log`, 24 passed, 6 s) was green but labelled the I81 / I83 rows with `LADDER_MM` (1 / 2 / 3 / 5 mm beyond the flat end) instead of the |z| rungs; the labels were changed to `RUNGS` (print only, no assertion touched) and the window re-run as W1b (`20261115T180304Z_GEO-36-I76p-W1b.log`, **24 passed in 2.34 s, elapsed 4 s**, `:248`, `:252`). (K-63) failures none (`:89`); (K-c63) non-empty on `:12517` 7→8; (P63) |(z_top − z_bottom)/2 − Z_C| 0.0 m and offsets `['1.0']` on both arms; (L17) 0.0 m; (N-row) miss 2.499734 K/W (> 1e-06), separation 0.272840. Status header printed (`ENTRY63_STATUS`, citing `:12532`). Record, no band, not a matched source: per rod watt (G·P_rod 3.740169e-05 W) rod top / bottom 397.8928 / 401.3429 K/W, gel twin 250.369 / 250.8551 K/W; rod / twin F2182 θ ratio 1.589225 / 1.599899; rod tip mean 399.6178 K/W over I81 gel / ti / hole at 17 mm = 59.98322 / 53.61394 / 61.72729, over SCALE·Φ_a(17 mm) 6.682151 = 59.80377. Ride-alongs (I78's 5, I83's 7) green. `git diff --stat`: `tests/validation/test_geo36_i76_f2182_table.py | 128 ++++++++++++++++++++++++-` / `1 file changed, 127 insertions(+), 1 deletion(-)`. §9 item 970 superseded (DONE-by-987).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 987 — `GEO-36` I76′ — per `docs/planning/chunks/GEO-36.md` "§9 item 987". First confirm that `docs/testing/xl-ledger.md`'s 2026-11-13 `GEO-36-step10-I32-prime` row is filled by a review (its RED verdict is expected, not a stop) and that the 11-14 weekly's "Rulings, 2026-11-14" is on `main`; if not, stop. Edit only `tests/validation/test_geo36_i76_f2182_table.py`: set `ENTRY63_LOG` to the ledger row's basename, add the `ENTRY63_*` literals copied from that log with their line numbers, the one printed status header (EM green; (E) at Δt₁/10 and rod (S2-tip) above band as readings, citing `:12532`), and the four new tests; every existing function and test byte-identical; import I83's module inside the test body. Do not edit `test_geo36_rod_f_clinical_heating.py`. No `src/` edit, no glob, no JSON. Before the window `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped; foreground, Bash `timeout: 660000`; pytest with -s; container `timeout -k 30`. Never edit a pattern or literal to rescue a red. You are the executor — do not spawn agents. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/GEO-36-I76p-*`, mark item 987 BLOCKED in the same commit, and report."

### §9 item 986 — 2026-11-15 (§10 I32c, chain I — rule (4), the 11-14 weekly's step; chain step (no solve))

*Enumerated by the 11-14 weekly at `PROJECT_PLAN.md:5127–5152`, under the entry 63 ruling at `:4909–4952` (b23fb062). This brief adds execution detail only. No band is set or moved: `E_BAND` stays 1e-09 and `S2_BAND` stays 1e-04. Family I32c, attempt 1 of 4. **Executor:** the implementer, in-slot. The weekly's target binds: "I32b's stub literals reproduced by string except the lines this step changes (listed in the brief), `20261108T181211Z_GEO-36-I32b-L-W0b.log:103–130`; W0b / W0 / W0n3 green at `-n 2 / 1 / 3`; Wreg green; `_tip_slope_rows`'s existing rows unchanged (U127′ imports it). No `src/` edit; no band moves. Landing ⇒ entry 71 READY."*

> **HOLD: read this before taking the item.** The weekly says the step is "serial on nothing — entry 63 has run, the module is free" (`:5127–5128`). **That premise is contradicted by the queue the same weekly commit wrote.**
> - `xxl` entry 69 (P33-d2, Saturday 11-28) is queued at `docs/testing/xxl-queue.d/20261128-GEO-38-P33-d2.env`, written in b23fb062. Its module is `tests/validation/test_geo38_f_clinical_box_heating.py`.
> - That module imports this one directly at `:61` (`from tests.validation import test_geo36_rod_f_clinical_heating as i32`).
> - It calls the two functions I32c edits: `i32._thermal_arm(…, merged=True)` at `:173–174` and `i32._thermal_failures(arm)` at `:273`. It also uses `i32._transient`, `_probe_theta`, `_refuse_degree`, `_box_preconditions`, `_p37_fixture` and the STUB constants.
> - An AST forward-closure walk at this review: entry 69's closure has 116 modules and **contains** `test_geo36_rod_f_clinical_heating.py`. Entry 65's `test_geo38_f_clinical_matrix.py` (55 modules) and `tests/environment/test_complex_mode.py` do not.
>
> By the module-hold rule (`PROJECT_PLAN.md:3982`, "a queued window's own module, every module that window imports … are held until it runs"; `:3986`, "a default-off knob does not lift it"), **the module is held until entry 69 runs on 11-28.** I32c is **not takeable as written.**
>
> The slot does not rule on this. It skips the item, or marks it BLOCKED (hold) with this paragraph, and the review rules. The options the review has, none chosen here:
> - **(α)** A ruled exception. Every I32c change sits behind `tip_projection=True` or the stub-only test nodes, and entry 69 calls `_thermal_arm` with the default `tip_projection=False`. So (i)–(iii) leave entry 69's executed path byte-equivalent. That would need proof: run Wreg, which includes entry 69's module, and show `_thermal_failures` unchanged on arms without `tip_projection`. See the note under (i) below.
> - **(β)** Put the repair in a new module that composes `_run_heat` and patches `_thermal_arm`, as I85 does. Entry 71's command would then name that module, not "the same … module" as entry 71 is written (`xl-pending.md:2627`).
> - **(γ)** Wait until after 11-28. Entry 71 (11-17) and entry 67 (11-19) would slip by about two weeks.
>
> Everything below assumes the hold is lifted by (α) or by date.

**What the step is.** It repairs two of I32b's thermal gates, so that entry 71 can re-run entry 63 with the asserted set entry 55 passed:
- **(i)** The asserted (E) returns to CN plus the Δt₁ BE step. The Δt₁/10 residual is printed as `(E10)`, beside its floor estimate.
- **(ii)** A third one-step BE at Δt₁/100 joins the tip rows, with the two-point and three-point Richardson intercepts printed. (S2-tip) is printed on the clinical body and stays asserted on the stub.
- **(iii)** The separating measurement on the stub: the same Δt₁/10 step with T carried from 0.

I85's composed module inherits the change. Its stub log is re-run and the changed print lines are re-recorded under the (1\*) licence.

**Facts verified against the tree (this review, by `grep` / reading; unlogged, not a record).**
- **The module** is `tests/validation/test_geo36_rod_f_clinical_heating.py` (1452 lines). It is entry 63's module (`xl` 11-13, run) and entry 71's.
  - `_thermal_arm` `:296–426`.
  - The third BE step: `res10 = _transient(sub, ft, q, dt=DT_1_BOX / 10.0, n_steps=1, theta=1.0)` at **`:354`**, under `if tip_projection:` (`:352–355`).
  - `out` is built at `:411–416`. `"E": max(res.max_balance_residual, res1.max_balance_residual)` at `:413` is entry 55's definition.
  - The fold: **`:418`** `if tip_projection:` and **`:419`** `out["E"] = max(out["E"], res10.max_balance_residual)`. `:420` already stores `out["E10"] = res10.max_balance_residual`, but nothing prints it.
  - The ruling's line citations (`:354`, `:418–420`) are **exact on the tree.**
- **The printed (E) line** is at `:403–409`: `(E) balance CN … BE … (band 1e-09, ASSERTED)`. It prints only CN and the Δt₁ BE step. That is why entry 63's printed (E) was green while the asserted maximum went red.
- **`_thermal_failures`** is at `:429–451`.
  - (E) at `:435–436`, against `arm["E"]`.
  - (S2-tip) at `:446–449`: `worst = max(arm.get("tip_s2_pi", arm["tip_s2"]).values())`, asserted at `S2_BAND`, rod arm only (`not arm["merged"]`).
  - `tip_s2_pi` is set from `tip_rows["dt1/10"]` at `:392`.
- **The tip rows** are at `:386–398`. `tip_rows = {"dt1": …, "dt1/10": …}` at `:389–391`, and the tip-ball print is at `:393–398`. The print is the literal at W0b `:116–117` / `:123–124`.
- **`_tip_slope_rows`** is at `:279–288`: plain arithmetic, rows `R_Pi`, `R_DG0`, `pq_q`. **It is not edited.** Its importers:
  - `tests/package/test_wf17_implant_rod_heating_yaml_tutorial.py:61` (U127′);
  - `tests/package/test_wf15_study_rod_heating.py:46–49`, together with `_projected_ball_comparands`;
  - the rod-cells arm at `:652–653`.
- **`T_INF`** is imported from `test_mat15_coil_driven_heating_phantom` (`:86–93`). The θ subtractions are at `:349`, `:359–360` and `:366`. The rod-cells arm passes `t_inf=T_INF, t_init=T_INF` explicitly at `:595`.
  - `_thermal_arm` calls `_transient(sub, ft, q, …)` with no `t_inf` / `t_init` kwargs. Read `_transient`'s signature before choosing the route for (iii): the helper's kwargs, or a module-local monkeypatch of `T_INF`. Monkeypatching `T_INF` also moves the subtractions above, so prefer explicit kwargs on the one extra call.
- **Importers of the module** (by grep):
  - P61's `test_geo38_f_clinical_box_heating.py:61` (**entry 69, queued**);
  - C31's `test_geo38_f_clinical_3t_box_heating.py:39` (entry 70, not queued);
  - I85's `test_geo36_rod_f_clinical_3t_heating.py:92` (entry 67, pending);
  - `test_ops83_s2tip_form_on_stub.py:32` (imports `E_BAND`, …);
  - `test_imp8_entry48_lumped.py:43`, `test_imp8_rod_line_source_fem.py:64`, `test_mat21_p86_three_pictures_3t.py:41–42`;
  - WF-15 / WF-17 (above);
  - readers that cite line numbers in prose only: `test_geo36_i76_f2182_table.py:5` (cites `:336–339`, the F2182 print form) and `test_mat21_clinical_probe_reader.py:13` (`:145–152`).

  Edits will shift line numbers that those docstrings cite. Leave the docstrings alone. Their readers parse logs, not the module's line numbers.
- **I85's module** composes `heat._run_heat` (`test_geo36_rod_f_clinical_3t_heating.py:143–153`, `:267`). Its stub log is `20261110T105744Z_GEO-36-I85-W0b.log`, with `-n 1` at `20261110T105845Z_GEO-36-I85-W0n1.log`. The W0b command runs three paths at `-n 2`, `timeout -k 30 180`: I85's module, this module and `test_geo36_rod_f_clinical_3t_drive.py`.
- **The validation target**, `20261108T181211Z_GEO-36-I32b-L-W0b.log:103–130`, is the `[GEO-36 I32]` stub block of `test_i32a_k10r_run_heat_end_to_end_on_stubs`. These literal lines change:
  - **:116, :117, :123, :124:** the tip-ball lines gain an `R_Pi dt1/100 …` token and the intercepts. On the rod arm (`:123–124`) the `(S2-tip) gate` token stays only if the stub branch keeps it asserted (it does, by the weekly).
  - **:120 and :127:** the CN (E) lines gain `(E10) … floor …`, or (E10) gets a new adjacent line. Prefer a **new line**, so that `:120` / `:127` stay byte-identical.
  - **New lines:** the T-from-0 measurement (rod and twin, or the stub's arm only).

  Every other line in `:103–130` must be byte-identical: the EM, (P), (S0), F2182, (B), 1 g ball, probe, ratio and K10-R summary lines. The `price:` line `:129` carries timings, so it is a reading by form only.
- **The W0b footer:** Status 0, **22 s elapsed**, `9 passed, 1 skipped in 21.44s` per rank, `[capture] rc=0`. The command is `-n 2`, `timeout -k 30 120`, the whole module, `-v -s --tb=short`, raw log `/workspace/logs/geo36-i32b-l-w0b-raw.log`.
  - W0 is the same at `-n 1` (`20261108T181245Z_GEO-36-I32b-L-W0.log`).
  - W0n3 is `-n 3` on `::test_i32b_k11_projected_tip_gate_on_checkerboard` (`20261108T181312Z_GEO-36-I32b-L-W0n3.log`).
  - Wreg is `-n 2`, `timeout -k 30 180`, over `test_ops83_s2tip_form_on_stub.py`, `test_ops84_s2tip_form_rod_pair.py`, `test_geo38_f_clinical_3t_box_heating.py` and `test_geo38_f_clinical_box_heating.py` (`20261108T181335Z_GEO-36-I32b-L-Wreg.log`).
- **The stub's K11 tests:** `test_i32b_k11_projected_tip_gate_on_checkerboard` is at `:1178`, and `test_i32b_k11f_thermal_failures_gate_projected_tip` is at `:1254`. K11-f feeds `_thermal_failures` synthetic arms and expects `red ['rod (S2-tip) 2.000e-04']` and the old-form `1.600e-03` (W0b tail). Its expectation follows from the stub branch (S2-tip still asserted). If the clinical branch makes (S2-tip) print-only through a flag on the arm, K11-f's synthetic arm must still carry the asserted form.
- **Entry 71** (`docs/testing/xl-pending.md:2627–`) is PENDING PREREQUISITE on I32c. Its command is "entry 63's command byte for byte except `XL_CHUNK` and the raw-log name; I32c changes the module, not the knobs". **I32c adds no knob**, so entry 71's command stays as written.
- **Queued windows today:** `xl-queue.d/20261116-GEO-38-C35w.env` (entry 65, `test_geo38_f_clinical_matrix.py`, closure free of this module) and `xxl-queue.d/20261128-GEO-38-P33-d2.env` (entry 69, closure **contains** this module; see HOLD).

**Code.** One file, plus I85's module only if one of its literal pins names a changed print. Read its test bodies for pinned `R_Pi` / `(E)` strings first. Every other file stays byte-identical: `src/`, `_tip_slope_rows`, `_projected_ball_comparands`, `_transient`'s home module, P61's / C31's modules, OPS-83 / OPS-84, WF-15 / WF-17.
1. **(i) The asserted (E).**
   - Delete `:419` (the fold), so `out["E"]` stays `:413`'s maximum over CN and the Δt₁ BE step.
   - Keep `out["E10"]` (`:420`).
   - Add `E10_FLOOR(…)`, an arithmetic helper returning ε·T∞·ρc·V_sub / (Δt·∫q), with `ε = np.finfo(float).eps`, `T∞ = T_INF + 273.15` if the route carries T in °C as stated (**check:** the ruling's ≈ 1.6e-08 at ≈ 40.7 W/m³ per watt must reproduce to one significant figure from the clinical numbers; otherwise use T_INF as carried), Δt = Δt₁/10 and `∫q = G·P` (the arm's `dep`).
   - Print a **new** line after `:403–409`: `[{label}] (E10) BE dt_1/10 {E10:.3e} beside floor ε·T∞·ρc·V_sub/(Δt·∫q) {floor:.3e} (ratio {E10/floor:.3g}; reading, no band)`.
   - **Note for option (α):** with `tip_projection=False`, `out["E"]` is unchanged by this edit, because `:418–420` are inside `if tip_projection:`. Entry 69's arms (merged, default) therefore execute identically.
2. **(ii) The Δt₁/100 row.**
   - Inside `if tip_projection:`, add `res100 = _transient(…, dt=DT_1_BOX / 100.0, n_steps=1, theta=1.0)`, with `imag` updated, a `th100` field and `"dt1/100"` added to `tip_rows` through `_tip_slope_rows`, called **unedited**.
   - Intercepts are readings. The two-point intercept is the linear Richardson extrapolation on the (Δt₁/10, Δt₁/100) pair: R₀ = (10·R₁₀₀ − R₁₀)/9. Print the (Δt₁, Δt₁/10) pair's intercept as well, for comparison with the ruling's −1.0e-06 / +2.2e-06. The three-point intercept is a quadratic fit through (Δt₁, Δt₁/10, Δt₁/100). Compute both on the **signed** ρc⟨θ¹⟩/(Δt⟨Πq⟩) − 1, because `R_Pi` is an absolute value and an intercept of |·| cannot show a sign. Derive the signed value as `lhs/pq − 1` locally, without editing `_tip_slope_rows`.
   - Print them on a new line per ball, at ≤ 7 significant digits.
   - **Clinical versus stub:** add a keyword `tip_gate=True` (default; stub and K11 keep asserting) to `_thermal_arm` and `_run_heat`. When the clinical `_run_heat` body path (not the stub) is active, set `out["tip_s2_pi_asserted"] = False`, and have `_thermal_failures` skip (S2-tip) only when that key is present and `False`. **Do not** delete `tip_s2_pi`. The key's absence means asserted, so K11-f and P61's merged arms are unchanged.
   - How the clinical body is told apart from the stub: read `_run_heat` `:711–846` and the stub test `:1112–1176`. Both take the same path except for `_stubs`. The cleanest discriminator is an explicit keyword from `test_geo36_i32_rod_heating` (`:847–854`, the gated clinical node), `tip_gate=False`. I85's `_run_heat_3t` passes it too (`False`) for entry 67. That would be **an I85 edit**, so the brief allows it, with I85's W0b re-run.
   - Without that, entry 67 inherits the asserted (S2-tip), which is against the ruling ("entry 67 rides I32c's repaired gates"). **Name the choice in the Outcome.**
3. **(iii) The separating measurement (stub only).**
   - In `test_i32a_k10r_run_heat_end_to_end_on_stubs`, or a new node on the same stub, rerun the Δt₁/10 one-step with T carried from 0, through `_transient`'s `t_inf=0.0, t_init=0.0` kwargs if they exist. Otherwise use `monkeypatch.setattr(<this module>, "T_INF", 0.0)` scoped to that call only.
   - Print `(E10) T from 20 °C {a:.3e}, T from 0 {b:.3e}, ratio {a/b:.3g} (floor class predicts ≳ 1e3; printed, no band)`.
   - On the stub (≈ 25× the clinical density) both residuals may sit near ε. Print, never assert. A ratio near 1 is the negative result below, not a red.

**Gates (asserted; no band set or moved).**
- Every W0b gate stays ASSERTED at its imported band: (M2) / (Z), (R1), (P) and its NEG, (S0), (B), (E) on CN + Δt₁, (E′), (S2), the stub (S2-tip) on `R_Pi dt1/10` at `S2_BAND`, and K10-R's six heads.
- **(Lit) the byte rule:** every literal line of W0b `:103–130` other than those listed above, by string, against the new W0b. This is the weekly's validation target.
- **(Rows) `_tip_slope_rows`'s existing rows unchanged:** the `R_Pi dt1`, `R_Pi dt1/10` and `R_DG0 dt1` tokens are byte-identical to W0b `:116–117` / `:123–124`.

**Negative control (asserted).**
- **(N-fold)**, a new stub node or an extension of K11-f. On a synthetic arm with `E` = 5e-10 and `E10` = 9.9e-09 (entry 63's shape), `_thermal_failures` returns no (E) red. With `E` = 2e-09, it returns `… (E) 2.000e-09`.
- **Ceiling:** the control must read red at 2× `E_BAND` and green at 0.5×, so it straddles the band within a factor of 2 on both sides.
- K11-f's existing `red ['rod (S2-tip) 2.000e-04']` stays (the stub branch). **(N-clin)** is the same synthetic arm with `tip_s2_pi_asserted=False` and `tip_s2_pi` 2.0e-04: no (S2-tip) red, with the value printed.

**Windows** (bare and foreground, CLAUDE.md form, Bash `timeout: 660000`, container `timeout -k 30`, each ≤ 600 s container-side; read `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first, expect 0). Commands are the I32b windows' exact form with labels and raw-log names changed:
1. Code. **Rule (i):** no module edit after its last window.
2. **W0b** `GEO-36-I32c-W0b`: `-n 2`, `timeout -k 30 120`, the whole module, `-v -s --tb=short`, raw `/workspace/logs/geo36-i32c-w0b-raw.log`. Expect 9 + new nodes passed and 1 skipped per rank. The prior run was 22 s; expect ≈ 25–30 s with the Δt₁/100 steps.
3. **W0** `GEO-36-I32c-W0`: the same at `-n 1`.
4. **W0n3** `GEO-36-I32c-W0n3`: `-n 3` on `::test_i32b_k11_projected_tip_gate_on_checkerboard`.
5. **Wreg** `GEO-36-I32c-Wreg`: Wreg's four paths at `-n 2`, `timeout -k 30 180`, plus `tests/package/test_wf15_study_rod_heating.py` and `tests/package/test_wf17_implant_rod_heating_yaml_tutorial.py`, the `_tip_slope_rows` importers (check that they are complex-mode stub tests that run under this command; if WF-17 needs its own harness form, give it its own window). **Entry 69's module is in this window.** Its stub nodes green is the (α) evidence.
6. **I85** `GEO-36-I32c-I85-W0b`: I85 W0b's three paths, `-n 2`, `timeout -k 30 180`. Re-record I85's changed print lines under (1\*): every other line is byte-identical to `20261110T105744Z_GEO-36-I85-W0b.log`, and the changed lines are listed by number in the Outcome.

**Price / Tier / Class / Cost.**
- Standard tier, `-n 1 / 2 / 3`, complex build, no solve beyond the stub (12³ box).
- **Class:** chain step (no solve).
- **Cost:** windows ≈ 22 + 22 + 10 + ≈ 60 + ≈ 40 s ≈ 2.5–3 min measured-scale; predicted ≈ 45 min of slot time (code plus six windows plus literal diffing).

**Traps already paid for.**
- Count verdicts per rank, by token. `re.search`, never `startswith`. Derived prints at ≤ 7 significant digits.
- The leak hook: never `--no-verify`. On a refusal, move the log unchanged to `logs/leak-held/`.
- Never delete a log, red or not (item 907).
- **The 11-10 18:00 ruling:** no 1e-12 (or 1e-9) identity on a one-step Δt₁/10 or Δt₁/100 increment. (E10), the Δt₁/100 residual and the T-from-0 line are **readings**.
- **"Every record at its prior"** (10:30 11-10): R_Pi and the intercepts are readings, and θ prints are records at ≤ 11 digits.
- **Pytest's double import** (I85 `:99`): I85 resolves names through `m["heat"]`. A new keyword on `_run_heat` must have a default, or I85's call breaks.
- `T_INF` monkeypatching leaks into `:349` / `:359–366` if it is not scoped. Prefer kwargs.
- Do not edit `_tip_slope_rows`, `_projected_ball_comparands`, `_transient` or any `src/` file. Do not move `E_BAND` / `S2_BAND`.
- Do not edit the F2182 print (`:336–339`, which I76′ reads on entry 63's log), and do not edit the docstrings that cite line numbers.

**Scope.** The thermal gate repair in one module, its stub proof, and I85's inheritance. **Not claimed:**
- any clinical θ, (E10) or tip-row value (entry 71's);
- an explanation of the 23× / 28× growth (iii only prints a separating measurement);
- a zero Δt-independent tip part (the intercepts are readings);
- the comparand question;
- entry 67's readings;
- any `src/` floor fix (declined twice).

**Negative result.**
- **(Lit) mismatch on a line not listed above:** something shared moved. Stop, BLOCKED, park on `attempt/GEO-36-I32c-*`.
- **An existing gate red on the stub, or Wreg / I85 red:** stop, BLOCKED.
- **T-from-0 ratio ≲ 10 on the stub:** the floor class is not separated on the stub. Record it as printed and continue. This is a finding for the review (the clinical body may still separate it in entry 71), not a red.
- **Floor estimate not reproducing the ruling's ≈ 1.6e-08 to one significant figure** (clinical inputs): print what it gives and flag it. Do not change the formula to match.
- No `PREDICTED RED:` line.

**Status it can move.**
- §10 I32c → DONE. `GEO-36` stays as it is.
- `docs/testing/xl-pending.md` entry 71 goes **READY in the same commit**: Status line updated to name the landing commit and logs, the command unchanged (no knob), and the new print forms listed for the leak check (`(E10)`, the Δt₁/100 row, the intercepts, the T-from-0 line).
- Entry 67 rides the repaired gates only if `tip_gate=False` is threaded through I85 (Code 2). Otherwise flag it.
- Known-issues 2026-11-06 stays open (until entry 71 reads green).

**Depends on:** the HOLD (entry 69's closure) being ruled, (α), or the date passing 11-28. Otherwise nothing open. Never in the same slot as an item that edits this module, I85's module, or P61's / C31's modules.

**Outcome commit.** The module (plus I85's if it was edited), the six logs, entry 71's READY Status in `xl-pending.md`, this file's Outcome paragraph (listing changed literal lines by number, for W0b and I85), and the §9 line marked DONE with the local slot: one commit, through the hook.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 986 — `GEO-36` I32c — per the brief "§9 item 986". FIRST confirm the review has ruled the HOLD (entry 69's queued `xxl` module `test_geo38_f_clinical_box_heating.py:61` imports this module); if it has not, stop and mark the item BLOCKED (hold), with no edit. Edit only `tests/validation/test_geo36_rod_f_clinical_heating.py` and, if the brief's Code 2 is taken, `tests/validation/test_geo36_rod_f_clinical_3t_heating.py`. Never edit `src/`, `_tip_slope_rows`, `_projected_ball_comparands`, `_transient`, `E_BAND` or `S2_BAND`, the F2182 print, P61 / C31 / OPS-83 / OPS-84 / WF-15 / WF-17. Remove the fold at `:419`; print `(E10)` beside its floor on a new line; add the Δt₁/100 row and the intercepts (readings); keep (S2-tip) asserted on the stub; print the T-from-0 measurement on the stub. Every W0b `:103–130` line not listed in the brief must stay byte-identical. Windows: W0b / W0 / W0n3 / Wreg / I85-W0b, in the I32b form; harness call bare, repo-relative, foreground, Bash `timeout: 660000`; container `timeout -k 30`; `pgrep -c python3` reads 0 first. Never loosen or move a band to rescue a red. On a red, park on `attempt/GEO-36-I32c-*` and mark §9 item 986 BLOCKED in the same commit. Landing marks `xl-pending.md` entry 71 READY in the same commit. Never delete a log. Never `--no-verify`. Do not spawn agents."

**Review ruling (2026-11-15 10:30 daily) — item 986 is HELD, for the 11-16 weekly.** The brief's hold box is confirmed at this review: `xxl` entry 69's module `tests/validation/test_geo38_f_clinical_box_heating.py` imports this module at `:61` and calls `i32._thermal_arm` (`:173`) and `i32._thermal_failures` (`:273`, `:425`, `:540`, `:607`, `:706`, `:930`), and entry 69 is queued (`xxl-queue.d/20261128-GEO-38-P33-d2.env`). The 11-09 ruling (c) says a default-off knob does not lift a hold, so option (α) is not the daily's to grant; option (β) changes entry 71's pre-registered command ("entry 63's command byte for byte except the chunk and raw-log names"), which is the weekly's; option (γ) is a date. **No slot takes 986 until the weekly rules.** Two facts for the weekly: the 11-14 ruling already lets the launcher take entry 69 on 11-21 if entry 66 is not READY by the 11-20 18:00 review, which would free the module a week early; and entry 71 cannot run before the review after 986 lands, so each day of hold moves entry 71 (and entry 67 behind it) by a day.
