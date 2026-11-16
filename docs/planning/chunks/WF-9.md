# WF-9 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 183 — 2026-10-02 (chain I, I7 — daily provisional, rule (6))

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.*

**Licence.** I6 (`IMP-3`) landed ✅, and subgoal (7d) is `WF-9` "after (7b)". The target is the §7 row's own: heating from ∫TF·E_tan gated against the direct lead-in-phantom solve. This brief writes no band, Done-when or target. **(7d) has no validation target of its own in §10, so the review took the row's; the weekly ratifies that reading or reverts the step.**

- **Executor:** implementer. **Tier:** heavy, complex; W1 at `-n 2`, W2 at `-n 4`.
- **Cost:** ≈ 14 measured / ≈ 50 predicted min.
- **Depends on:** nothing.

**Why this fixture.** The `IMP-2` / `GEO-36` rods are PEC holes: the without-lead field does not exist on the same unknowns. A with/without comparison across two meshes carries interpolation error, and no imported band covers it.

A lead of finite conductivity on **one** mesh avoids that. The only difference between the two systems is ΔA = jωμ₀Δσ·M_lead, supported on the lead cells (the lead is ≥ 7 mm from every wall). So q_with − q_without is the response, with the lead present, to an impressed current J_eq = Δσ·E_inc on the lead. That is the row's "local impressed source", split into steps along the path.

A is complex-symmetric (real basis, scalar ε_c; `TH-1` asserts A = Aᵀ at 1e-10). So ℓᵀA⁻¹b_j = b_jᵀA⁻¹ℓ: 24 stepped-source solves equal one tip-source adjoint solve read along the lead (Jin §6.3 eq (6.67), §6.4 eqs (6.78)–(6.85)). The Hermitian route (6.86)–(6.88) does not apply to lossy problems; it is control C1.

**Exact versus not:**
- **Exact (LU round-off):** superposition, reciprocity and the volume fold.
- **Not exact:** the ISO path fold Σ_k TF_k·Ē_k with unit-ẑ sources and segment means. It misses by the transverse term ∫w_⊥·E_inc,⊥, which does not shrink with N, and by the within-segment covariance. The along-axis part should fall roughly as N⁻²; the cross-section part leaves a thin-wire floor.
- **Bands:** the exact parts use `DISCRETE_IDENTITY_RTOL` 1e-6, imported. The path fold has no imported band; this step records it, and its band is the weekly's.

**Fixture.** Tests only: new `tests/validation/test_wf9_transfer_function_identity.py`, with `test_wf9_census` (no solve) and `test_wf9_identities`.
- **Mesh.** `dolfinx.mesh.create_box` on [−8, 8]² × [−20, 20] mm, (16, 16, 40) hexes split into tetrahedra: 61 440 cells. N1curl degree 1, predicted 76 360 dofs.
- **Host.** `IMP-1`'s: `SIGMA_HIGH` 0.57 and `FREQUENCY_HZ` 64e6 (`test_lossy_sphere_sar.py:74,77`), `EPSILON_R_HOST` 78 (`test_imp1_pec_sphere_lossy_host.py:65`).
- **Lead.**
  - Cells with midpoint |x|, |y| < 1 mm and |z| < 12 mm (`dolfinx.mesh.compute_midpoints`), cell tag 2.
  - σ = `SIGMA_WIRE_S_PER_M` 800 (`test_port_gap_voltage_impedance.py:235`), εᵣ 78. Since Δεᵣ = 0, J_eq = Δσ·E.
  - "Without the lead" is the same mesh with no `material_map`.
- **Layers.** j = 0…23, one per mm of z. N ∈ {3, 6, 12} is formed by summing layers.
- **Tip functional.** q(e) = (1/V_b)∫χ_b e_z dx, where χ_b is a DG0 mask of host cells whose midpoint lies within 3 mm of (0, 0, 12 mm).
- **Drives.**
  - **D1:** the imposed trace E₀ẑ on all walls (`IMP-2`'s axial drive).
  - **D2:** an impressed ẑ current block (cells with midpoint 4 < x < 6, |y| < 1, 2 < z < 10 mm), with load −jωμ₀∫J·v̄, unprojected.
- **Solves.** One factorisation each; everything else is a held back-substitution. Free each solver after use.
  - S_oD and S_wD: Dirichlet, lead out / in, for D1.
  - S_o0: zero trace, lead out, for D2's E_inc.
  - S_w0: zero trace, lead in. A full solve on the adjoint load ℓ = ẑχ_b/V_b gives w. Then held solves for D2 with the lead, the 24 unit sources Δσẑχ_j (TF_j^step), and the 24 + 24 partition sources Δσ·E_inc·χ_j (D1, D2).

**Anchors (asserted; bands imported):**
- **(A0) Census.**
  - 61 440 cells, 576 lead cells and 192 block cells, exact.
  - Σ_j χ_j = χ_lead cell by cell.
  - Lead volume 9.6e-8 m³ at `EXACT_IDENTITY_RTOL` 1e-10 (`tests/mesh/test_birdcage_symmetry_quadrant.py:47`).
  - χ_b, χ_lead and χ_src pairwise disjoint; ball count > 0 (predicted [400, 700]).
- **(A1) Symmetry.** ‖A − Aᵀ‖_F ≤ 1e-10‖A‖_F on the with-lead form (`test_time_harmonic_mms.py:255`). This is (T2)'s premise.
- **(T1) Superposition, D1 and D2.** |Σ_j q(A⁻¹b[Δσ·E_inc·χ_j]) − (q_with − q_without)| ≤ `DISCRETE_IDENTITY_RTOL` 1e-6 (`test_birdcage_power_identity.py:227`) × |q_with − q_without|.
- **(T2) Reciprocity.** max_j |TF_j^step − TF_j^adj| ≤ 1e-6·max_j |TF_j^adj|, with TF_j^adj = −jωμ₀Δσ∫χ_j w_z dx.
- **(T3) The row's identity, exact form, D1 and D2.** −jωμ₀Δσ∫χ_lead w·E_inc dx = q_with − q_without within 1e-6 relative.
- **(H) Held vs fresh.** TF_23 from a fresh `solve(current_density=…, project_source=False)` equals the held value within 1e-6.
- W2 re-asserts every anchor at `-n 4`. The W1/W2 agreement of TF_j is printed (predicted ≤ 1e-10).

**Record (printed only; the path-fold band is the weekly's):**
- Q_N = Σ_k TF_k·Ē_k for N = 3, 6, 12 and 24, beside Q_exact, with the error split into the transverse term and the within-segment covariance.
  - D1: predicted ≤ 5 % at every N, flat within a factor of 3.
  - D2: predicted to fall from N = 3 to 12; N = 24 in [1e-4, 1e-1].
- |Q_N|² against |Q_exact|², the ISO form.
- q_with/q_without per drive (D1 predicted [1.2, 10]).
- The |TF_j| and arg TF_j profile; dofs, times and summed RSS.

**Negative controls (rule (e)):**
- **(C1) Hermitian read-back:** TF_j computed with conj(w_z), which is what `ufl.inner` would give.
  - *Asserted:* wired, deviation > 1e-6.
  - *Predicted:* [50 %, 200 %], with a ceiling of 200 % since |a − ā| ≤ 2|a|.
- **(C2) Total-field fold:** −jωμ₀Δσ∫χ_lead w·E_with.
  - *Asserted:* wired, > 1e-6.
  - *Predicted:* a miss of [80 %, 105 %], because E inside the lead is ≈ 4 % of E_inc.

**Windows.** Bash `timeout: 660000`, foreground, each harness call bare and whole.

Price basis: `IMP-1` a/6 (93 146 degree-2 dofs) took 13.4 s for solve + post at 2.01 GiB, `-n 2` (`20261001T033100Z_IMP-1-land-W1.log:371–372`). A held drive is ≤ 3.57 s at 369 387 cells (`WF-12` step 2). This step does five factorisations and ≈ 74 back-substitutions.

- **W1,** predicted [90, 300] s:

  `scripts/testing/run_and_log.sh WF-9-step0-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf9_transfer_function_identity.py::test_wf9_census tests/validation/test_wf9_transfer_function_identity.py::test_wf9_identities -v -s --tb=short'"`
- **W2,** predicted [70, 280] s: the same command with label `WF-9-step0-W2` and `-n 4`.

**Traps:**
- **`boundary_condition` defaults to NATURAL.** Set `"pec_zero_tangential_a"` on every problem, or the trace is ignored.
- **`solve(current_density=…)` projects by default.** Projection removes ∇ψ and breaks (T1) and (H), so pass `project_source=False`. The held loads carry −jωμ₀ explicitly (`time_harmonic.py:508`).
- **The held path re-applies its own BC values.** Zero-trace sources run on S_w0, never on S_wD (`PORT-19`).
- **`ufl.inner` conjugates its second argument.** Read q, TF and the fold through components or `ufl.dot`.
- **`_assemble_real` drops the imaginary part.** Assemble and reduce complex values.
- **No UFL comparisons or conditionals.** Build the masks as DG0 from numpy midpoints (`OPS-22`).
- **Owned vs ghost cells.** Midpoints and masks cover owned + ghost cells, but counts sum owned cells only. Keep collectives outside `rank == 0`.
- **Harness discipline:** always `-s`; node IDs only, no `-k`; no pipes. On Status 124, run `pgrep -c python3`, then the cache sweep; never lengthen the timeout. No module edit after W1 (rule (i)).

**Scope.** Degree 1, one square lead, two drives, a linear tip functional. **Not covered:**
- a PEC or insulated lead;
- F-small or a coil drive;
- the path fold's band;
- ISO 10974's measured TF or its calibration;
- published heating data;
- a `src/` post-processor;
- h-convergence;
- thermal.

**Status it can move.** §7 `WF-9` ⬜ → 🟡 when W1 and W2 are green. There is no ✅ move: the row has no Done-when, so rule (j) is not engaged. The weekly writes the Done-when and step 1.

**Negative result:**
- **(A0) or (A1) misses:** stop and mark BLOCKED.
- **(T1), (T2), (T3) or (H) red:** park on `attempt/WF-9-step0-<ts>`, write a known-issues entry and mark BLOCKED (rule (d)). Never widen a band.
- **(T2) red with (T1) and (T3) green:** the read-back and the source functional disagree. Print both vectors.
- **A control not wired:** a test bug. Fix it in-slot and re-run W1.
- **A record outside its prediction:** report it; it moves no status.
- There is no PREDICTED RED line.

**Outcome — DONE 2026-10-03 (19:30 CDT slot, implementer).** Every anchor is green at both widths:
- W1 (`-n 2`): `20261003T003327Z_WF-9-step0-W1.log`, 129 s.
- W2 (`-n 4`): `20261003T003546Z_WF-9-step0-W2.log`, 104 s.

Measured values, W1 / W2:
- **(A0)** 61 440 cells, 576 lead cells, 192 block cells and 624 ball cells. The lead volume is 9.6e-8 m³. Partition and disjointness hold.
- **(A1)** 1.102e-16.
- **(T1)** D1 1.7e-13 / 4.3e-13; D2 4.8e-14 / 2.1e-12.
- **(T2)** 1.48e-12 / 1.43e-12.
- **(T3)** D1 1.4e-12 / 1.7e-12; D2 1.1e-12 / 3.4e-12.
- **(H)** 1.2e-14 / 6.0e-13.
- **W1/W2 TF agreement:** 1.25e-12.

Controls:
- **C1** is wired at 199.87 %. That is just under its 200 % ceiling, because TF_adj is nearly real (phase −2° to −5° along most of the lead).
- **C2** is wired at 94.6 % (D1) and 97.3 % (D2).

Record:
- **D1 path fold:** 2.4e-6 at every N. It is flat, because E_inc is near-uniform and axial.
- **D2 path fold:** the error is 0.49 / 0.21 / 0.035 / 9.7e-3 at N = 3 / 6 / 12 / 24. At N = 24 the transverse term (8.1e-3) dominates, which is the thin-wire floor. \|Q_N\|²/\|Q\|² is 0.26 → 0.98.
- **\|q_with/q_without\|:** D1 2.03 and D2 0.21.
- **TF_j:** \|TF_j\| peaks at j = 19. The last two layers, j = 22 and 23, flip phase (≈ −170°) next to the tip ball. This is a record only.
- **Size and cost:** 76 360 dofs. Summed RSS is 1.74 GiB at `-n 2` and 2.47 GiB at `-n 4`.

A1's 1e-10 is restated as `SYMMETRY_RTOL`, citing `test_time_harmonic_mms.py:255`, which holds only a literal. Every record is inside its prediction. Row ⬜ → 🟡. The Done-when and the path-fold band are the weekly's to write.

### §9 item 219 — 2026-10-04 (chain I8)

*Drafted by a general-purpose agent for the 2026-10-04 03:00 daily review; the review reviews it.*

**The target is the 2026-10-03 weekly's `WF-9` ruling (PROJECT_PLAN §10, "`WF-9` — Done-when written"), verbatim.** Step 1 puts the lead in the coil-driven F-small phantom at 64 MHz, on one mesh with the lead's cells kept and switched by material (σ = 800 S/m against gel), so the adjoint fold is an exact discrete identity. Asserted, for each of the four port drives:
- the adjoint TF folded with the lead-free field equals the direct with-lead tip response, at `DISCRETE_IDENTITY_RTOL`;
- the 1-D path fold ∫TF(z)·E_tan(z) dz reproduces it within **5 %**;
- the lead-free Born fold misses.

The result is ⇒ ✅ as a fixture. Not claimed: a measured TF, heating, or a real lead. Chain I8 adds that this needs an additive keep-the-rod-cells option on `GEO-36`'s generator, with the CAD identities re-asserted. This brief supplies execution detail and changes no target, band or Done-when.

**Split: one slot cannot hold this.** On paper it is ≈ 40 min of code (a `src/` mesher branch, a CAD module and a solve module) plus ≈ 17 min of windows: about 170 + 170 + 100 s for the mesher half and about 200 + 180 s for the solve half. That is ≈ 60–90 min predicted against a 1 h box. The measured halves are smaller: `GEO-36` step 2 took 8 min and `WF-9` step 0 took 14. Even so, a mesher FAIL would strand the solve half mid-slot.
- **Part A** (`GEO-36` step 3, the generator option and its CAD identities) is this item.
- **Part B** (the solve) is a second §9 item, serial on Part A.
- If Part A commits before minute 30, take-next lets the same slot run Part B.

#### Part A — keep the rod's cells (`GEO-36` step 3; no solve)

- **Executor:** implementer. **Tier:** heavy (each window ≈ 170 s), real build, `-n 1` then `-n 2`. **Cost:** ≈ 11 measured (chain no-solve median) / ≈ 40 predicted min. **Depends on:** nothing.
- **Mesher change (additive; rule (c), disclosed).** Accept an optional key `"keep_cells": bool` in the `implant_rod` mapping. Validate it in `_normalise_implant_rod` (`src/fem_em_solver/io/mesh.py:4538`): a non-bool raises before any gmsh call. Absent or `False` takes today's hole branch untouched.
  - When it is `True`:
    - Cut the rod from the phantom with `removeTool=False`, as now. Cut it from the air box also with **`removeTool=False`** (`:5041–5043` drops it today).
    - Pass the rod solid into the fragment as its own group `"implant"`, physical cell tag **`BIRDCAGE_IMPLANT_CELL_TAG`**. This is a new module constant beside `:272`. Choose a value free of the birdcage's cell groups (1/2/3, 100+i, 200+i) and facet groups (1, 401, 402, 501, 502), and confirm it is free in the code; 403 is suggested.
    - The wall is the boundary of the `"implant"` pieces, all of which are **two-use** faces. The single-use filter at `:5326–5341` must not be applied in this branch, or it finds nothing and raises.
    - Emit it as dim-2 group 402, which is now an **interior** facet group. Sizing fields are unchanged; inside the rod the Distance field already gives r/3.
  - `implant_rod_cad` gains `"keep_cells": True`, `"cell_tag"`, and the OCC mass of the kept solid.
  - The `as_hole` / `symmetry_quadrant` refusals stand.
- **Tests:** extend `tests/mesh/test_birdcage_implant_rod.py`.
  - Add a third build, `_build_rod(implant_rod=dict(ROD, keep_cells=True))`, **after** flag-off and the hole build in the `builds` fixture (`GEO-23`: a failed gmsh build contaminates its process).
  - Add a helper `_build_kept_rod()` for Part B.
  - Register the cell count as `KEPT_ROD_CELL_RECORD` after W2. It is a first record, printed, with a guard of 320 000 checked before returning.
- **Anchors (asserted; bands imported from the module's own imports).**
  - (K1) `cad_mass_by_group["phantom"]` = πR²H − V_rod at `CAD_ANALYTIC_BAND`, as A1 does.
  - (K2) `cad_mass_by_group["implant"]` = πr²ℓ + 4πr³/3 (`_rod_volume()`) at `CAD_ANALYTIC_BAND`.
  - (K3) Σ CAD(kept) = Σ CAD(flag off) at `EXACT`. Nothing is removed, which is the kept/hole difference.
  - (K4) The OCC wall area = `_rod_area()` at `CAD_ANALYTIC_BAND`. The meshed tag-402 area ÷ `_rod_area()` lies in `TERMINAL_AREA_BAND`, with the area taken by a **dS** integral (see Traps). Predicted 0.99625, the hole build's.
  - (K5) The census:
    - the cell-tag set equals flag-off's ∪ {`BIRDCAGE_IMPLANT_CELL_TAG`};
    - the owned census sums to `size_global`;
    - tag 402 has > 0 facets, **0 exterior**, and every interior one is bounded by exactly one tag-3 cell and one implant cell;
    - the implant tag is absent with the flag off and in the hole build.
  - (K6) The OCC wall centroid at (0.015, 0, 0) to `EXACT`·ℓ, as A6 does.
  - (K7) A non-bool `keep_cells` raises without meshing (added to `test_a7_…`'s parametrisation).
  - A1–A7 on the hole build re-run unchanged. That is the hole path's rule (c) re-gate, so 167 371 must print again.
- **Printed:** the meshed implant volume ÷ V_rod (predicted [0.990, 0.999], an inscribed-faceting deficit), kept cells (predicted [180 k, 230 k]: 167 371 plus [15 k, 60 k] rod cells at h ≈ r/3), mesh time, and the meshed wall centroid.
- **Negative controls (rule (e)).**
  - (C1k) **Asserted.** The hole build misses (K3) by V_rod/ΣCAD. `GEO-36` W2 measured the same comparison on the same fixture: 8.54e-6, its (C2) print (`20261002T110559Z_GEO-36-step2-W2.log`). It is asserted as a miss > `EXACT`.
  - (C2k) **Predicted, printed.** The implant tag is dropped from the census and the owned census misses by the implant cell count, against a ceiling of 100 % of that count.
- **Windows** (bare, whole, foreground, Bash `timeout: 660000`; basis `GEO-36` W1 124 s / W2 121 s for two builds, plus a ≈ 45–50 s third build):
  - W1 (≈ 175 s): `scripts/testing/run_and_log.sh GEO-36-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod.py -s -v --tb=short'"`
  - W2, the same at `-n 2`, labelled `GEO-36-step3-W2`. This is the window of record.
  - W3, the default-path re-gate (≈ 100 s). It excludes the 🔴 OPEN pre-existing red `test_birdcage_volumes_partition_the_box` (known-issues 2026-08-27, `OPS-26` / `GEO-23`; `GEO-36` W3 attributed it):
    `scripts/testing/run_and_log.sh GEO-36-step3-W3-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_conductor_hole.py tests/mesh/test_birdcage_port_tags.py::test_birdcage_port_layout_diagnostics_match_the_closed_forms tests/mesh/test_birdcage_port_tags.py::test_birdcage_port_layout_rejects_too_small_or_overlapping_port_regions tests/mesh/test_birdcage_port_tags.py::test_ungraded_build_still_hits_the_documented_floor -s -v --tb=short'"`
    The 116 085 record must print unchanged.
- **Traps.**
  - **`_exterior_facet_area` and `_facet_check`'s moments integrate `ds`, so they read 0 on an interior group.** Measure 402's area and centroid with `dS(402)` on a `("+")` restriction; each interior facet is counted once.
  - Count interior facets per owned facet, and look up both cells' tags including ghosts (`shared_facet` ghosting holds the neighbours).
  - A FAIL is recorded verbatim with its surface pairs, never retried; the `GEO-9` overlapping-facet failure is the one to expect if the rod enters the fragment twice.
  - Use `-s`, node IDs only, no pipes. Make no edit after W2 (rule (i)).
- **Scope.** No solve or material; F-small only; no `as_hole`.
- **Status it can move.** The `GEO-36` row gains a step-3 annotation and keeps its glyph (✅). `KEPT_ROD_CELL_RECORD` is registered, and Part B becomes queueable.
- **Negative result.**
  - A gmsh FAIL: write a known-issues entry, park on `attempt/GEO-36-step3-<ts>`, and mark the item BLOCKED (rule (d)). Part B stays blocked.
  - (K4) out of band: report it; never widen.
  - W3 moves the default path: revert `src/` and stop.
  - There is no `PREDICTED RED` line.

**Outcome — DONE 2026-10-04 (15:00 CDT slot, implementer).** K1–K7 green at `-n 1` (`20261004T200412Z_GEO-36-step3-W1.log`, 188 s) and `-n 2` (`20261004T200730Z_GEO-36-step3-W2.log`, 181 s, of record); W3 re-gate green (`20261004T201038Z_GEO-36-step3-W3-regate.log`, 114 s, 116 085 unchanged). `KEPT_ROD_CELL_RECORD = 182 371` (167 371 + 15 000 implant cells) and `_build_kept_rod()` exist; Part B is queueable. One disclosed deviation: 402 in the kept build is rebuilt dolfinx-side from the (3, 403) cell tags, not emitted as an interior gmsh group (known-issues 9 rule, `_interface_facet_tags`). Meshed implant volume / V_rod 0.986667 is below the printed [0.990, 0.999] prediction (record only). Detail: `docs/planning/chunks/GEO-36.md` § Step 3.

### §9 item 224 — 2026-10-04 (chain I8, Part B)

*Renumbered by the 2026-10-04 03:00 review: the drafter called this item 220, which `WF-13` holds. Part A is §9 item 219.*

**Outcome — DONE 2026-10-04 (15:00 CDT slot, implementer, second item).** `tests/validation/test_wf9_lead_in_coil_phantom.py`; W1 `-n 2` 14 passed, 190 s (`20261004T201655Z_WF-9-step1-W1.log`); W2 `-n 4` 14 passed, 165 s (`20261004T202017Z_WF-9-step1-W2.log`). (A0) 182 371 cells = record; 15 000 implant, 7 605 ball cells, all 32 slabs populated. (M0) σ 800 / 0.5, εᵣ 78 / 78 exact. (A1) 9.757e-17. (i′) with 1.5e-15, without 1.8e-15 (W1) / 1.4e-15 (W2). (H) 0.0 (W1) / 2.7e-14 (W2). (T3) per drive P1–P4: W1 2.9e-13 / 3.7e-12 / 6.7e-13 / 1.0e-12; W2 1.6e-13 / 2.9e-13 / 4.3e-13 / 4.9e-14. (PF) N = 32: 1.44e-4 / 1.47e-3 / 1.51e-4 / 1.02e-3 (band 0.05). (BF) 1 399 / 1 435 / 1 389 / 1 344 % (must exceed 5 %; predicted [100 %, 3000 %]). C1 ≈ 199.5 %, C2 ≈ 93.5 % on every drive (both IN). |q_with/q_without| 2.33 / 2.37 / 2.39 / 2.29 (all IN). W1/W2 TF agreement 2.86e-12. **Records outside prediction (no status effect):** (PF) on P2 / P4 is 1.47e-3 / 1.02e-3, *below* the predicted [3e-3, 4e-2] — the transverse term is 1.5e-3 / 1.1e-3 of Δq, smaller than the on-axis η_z argument implied; summed RSS 3.58 GiB at `-n 4` (the ≤ 3.5 GiB prediction was for `-n 2`, where it read 2.91). The path-fold error is not monotone in N on P2 / P4 (N = 8 ≈ 1.7e-4 on P4) because the transverse floor and the within-slab covariance partly cancel.

#### Part B — the solve (`WF-9` step 1; serial on Part A)

- **Executor:** implementer. **Tier:** heavy, complex, W1 at `-n 2` and W2 at `-n 4` (step 0's widths). **Cost:** ≈ 18 measured (chain-solve median) / ≈ 50 predicted min. **Depends on Part A:** `KEPT_ROD_CELL_RECORD` and `_build_kept_rod` exist, and `GEO-36-step3-W2` is green in `test-results.md`.
- **Fixture.** Create a new `tests/validation/test_wf9_lead_in_coil_phantom.py`. Its module-scoped fixture `_case()` builds once, and two nodes read it: `test_wf9_step1_census` and `test_wf9_step1_identities`.
  - **Build.** `built = _build_kept_rod()`. Check the cell guard before any solve. Then call `case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)`, which is `IMP-3`'s entry (`test_imp3_rod_birdcage.py:299–300`).
  - **Two problems on one mesh.** `gel = case["problem"].material_map[PHANTOM_CELL_TAG]`. Then:
    - `p_with = replace(problem, material_map={**map, IMPLANT: HomogeneousMaterial(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r)})`;
    - `p_without = replace(…, {**map, IMPLANT: gel})`.
    - Δσ = `SIGMA_WIRE_S_PER_M − gel.sigma`, computed and never typed. With Δε = Δμ = 0, J_eq = Δσ·E.
  - **The tip functional** is q(e) = (1/V_b)∫χ_b e_z dx. χ_b is a DG0 mask of **tag-3** cells whose numpy midpoint lies within `BALL_RADIUS` of `TIPS["top"]` (both imported; (0.015, 0, 0.016) is the top cap tip).
  - **Layers.** χ_k marks implant cells by midpoint z in 32 equal 1 mm slabs over [−16, 16] mm. That is step 0's 1 mm pitch; its 0.97 % was measured at that pitch. N ∈ {4, 8, 16} are formed by summing slabs and are recorded.
- **Solves.** Follow `_solve_driven`'s pattern (`test_birdcage_b1_plus_map.py:142`): a `TimeHarmonicSolver` per arm, with **every** sheet's `lumped_port_bilinear_term` in `extra_bilinear_terms`.
  - **With arm:**
    - `solve(extra_linear_terms=[ℓ])` gives w, where ℓ = ẑχ_b/V_b and is read with `ufl.dot`;
    - four held solves, each carrying one driven sheet's `lumped_port_linear_term`, give E_with,d for d = P1…P4.
  - **Without arm:** four drive solves (one factor plus three held) give E_without,d, and one held solve on ℓ gives w₀ for the Born fold.
  - **(H) fresh:** `_solve_driven(dict(case, problem=p_with), "P1")` on its own solver.
  - Free each solver after use (step 0's `_free`).
  - Total: three factorisations and about 10 back-substitutions.
- **Anchors (asserted; bands imported).**
  - **(A0) Census.**
    - The cells equal `KEPT_ROD_CELL_RECORD` within `CELL_COUNT_BAND`.
    - The implant, ball and slab counts are each > 0; Σ_k χ_k = χ_lead cell by cell; χ_b and χ_lead are disjoint.
    - **(M0)** DG0 σ and εᵣ read back on the implant cells are rank-reduced min = max = 800 / gel εᵣ (with) and gel σ / εᵣ (without).
  - **(A1)** ‖A − Aᵀ‖_F ≤ `SYMMETRY_RTOL`·‖A‖_F on the with-arm volume form plus all four sheet terms (step 0's `:235–245` pattern). Predicted ≈ 1e-16.
  - **(i′) Loaded reciprocity, both arms.** Y_ij = I_i(drive j)/V_src (`sheet_terminal_current`, reduced) gives ‖Y − Yᵀ‖/‖Y‖ ≤ `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`). Predicted ≤ 2e-15.
  - **(T3) The Done-when identity, per drive d.** |−jωμ₀Δσ∫χ_lead w·E_without,d dx − (q_with,d − q_without,d)| ≤ `DISCRETE_IDENTITY_RTOL` (`test_birdcage_power_identity.py:227`) × |q_with,d − q_without,d|. Predicted ≤ 5e-12 (step 0: ≤ 3.4e-12). The value q_without + fold against q_with is printed.
  - **(H)** Fresh q_with,P1 equals held within `DISCRETE_IDENTITY_RTOL`, predicted ≤ 1e-12. This guards the `PORT-19` premise: the four drives share one operator only for matched terminations.
  - **(PF) Path fold, per drive, N = 32.**
    - Q₃₂ = Σ_k TF_k·Ē_k, with TF_k = −jωμ₀Δσ∫χ_k w_z dx and Ē_k = ∫χ_k E_without,z dx / V_k.
    - Asserted: |Q₃₂ − Δq_d|/|Δq_d| ≤ 0.05. This is the weekly's band, a module constant `PATH_FOLD_BAND` citing PROJECT_PLAN §10's `WF-9` ruling.
    - Predicted: P1 / P3 in [1e-4, 2e-2], since their on-axis background is axial (`IMP-3` W2: η_z 0.990 / 0.902). P2 / P4 in [3e-3, 4e-2], since it is mostly transverse (η_z 0.405 / 0.443).
    - The error split is printed, as step 0 did: the transverse term ∫w_⊥·E_⊥, and the within-slab covariance.
  - **(BF) Born, per drive.** B_d = −jωμ₀Δσ∫χ_lead w₀·E_without,d.
    - Asserted: |B_d − Δq_d|/|Δq_d| > `PATH_FOLD_BAND`. This is the Done-when's "misses", and the band separates it from (PF).
    - Predicted: a miss in [100 %, 3000 %], which has no ceiling. The quasi-static internal axial factor is 1/(1 + L_z(ε_c,rod/ε_c,gel − 1)) ≈ 0.07 at aspect 16 and contrast ≈ 1.4e3, so Born overstates the in-lead field ≈ 14×.
- **Negative controls (rule (e)), carried from step 0; both are asserted only as wired (> `DISCRETE_IDENTITY_RTOL`), and their factors are predicted.**
  - (C1) Hermitian read-back with conj(w). Predicted [10 %, 200 %], with a ceiling of 200 %. Step 0 measured 199.87 % on another fixture.
  - (C2) Total-field fold with E_with. Predicted [80 %, 105 %]. Step 0 measured 94.6 / 97.3 %.
  - (BF) is the third control, and it is asserted by the Done-when.
- **Printed:** |q_with/q_without| per drive (predicted P1 [1.5, 6], P3 [1.3, 5], P2 / P4 [0.8, 3], scaled from `IMP-3`'s PEC 1 g Q); |TF_k| and arg TF_k along z; Q_N for N = 4…32; the W1/W2 TF agreement (predicted ≤ 1e-10); dofs, phase times and summed RSS (predicted ≤ 3.5 GiB at `-n 2`, against `IMP-3` W1's 2.61 GiB).
- **Windows** (bare, whole, foreground, Bash `timeout: 660000`). Basis:
  - `IMP-3` W1: 241 s at `-n 2` for two builds, two sweeps and a fresh control solve;
  - the hole build's mesh took 42 s;
  - step 0's five factorisations plus 74 held solves took 129 s.

  The commands:
  - W1, predicted [120, 300] s: `scripts/testing/run_and_log.sh WF-9-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf9_lead_in_coil_phantom.py::test_wf9_step1_census tests/validation/test_wf9_lead_in_coil_phantom.py::test_wf9_step1_identities -v -s --tb=short'"`
  - W2, predicted [100, 280] s: the same command with label `WF-9-step1-W2` and `-n 4`.
- **Traps already paid for.**
  - **Unmapped kept cells silently become vacuum.** `build_four_port_sweep`'s default material is σ = 0, εᵣ = 1 (`test_port_birdcage_four_port.py:254`). (T3) would still hold on a vacuum rod, so (M0) is the guard.
  - **Every sheet's bilinear term goes into both factors.** A held solve with the wrong termination solves a different problem silently (`PORT-19`, `WF-12`), which is why (H) exists.
  - **`ufl.inner` conjugates its second argument.** Read q, TF and the folds with `ufl.dot` or components.
  - **`_assemble_real` drops imaginary parts.** Reduce complex values.
  - **Masks are DG0 from numpy midpoints, never UFL conditionals or `sqrt` comparisons** (`OPS-22`, `WF-6` 4c). Masks cover owned + ghost cells; counts sum owned cells only.
  - **Tag lookups must include ghost cells.** The facet tags for the sheets come from the build tuple and the sweep, as `IMP-3`'s do.
  - Keep collectives outside `rank == 0` branches.
  - **Harness discipline:** `-s`; node IDs only, no `-k`; no pipes. On Status 124, run `pgrep -c python3`, then the any-orphan-`.c` cache sweep, and never lengthen the timeout. Make no module edit after W1 (rule (i)).
- **Scope.** Degree 1, 64 MHz, the untuned solid-coil F-small 4×4 with 50 Ω sheets, one σ = 800 S/m capped rod on leg 1's azimuth, and a linear tip functional at one tip. **Not covered:**
  - a measured or ISO 10974-calibrated TF;
  - heating or SAR;
  - a PEC or insulated lead;
  - a real lead;
  - h-convergence;
  - quadrature drives;
  - F-human;
  - Palace (`IMP-3` step 3 is the twin).
- **Status it can move.** §7 `WF-9` 🟡 → ✅ "as a fixture" when W1 and W2 are green, on the Done-when the queueing review copies into the row (rule (j)).
- **Negative result.**
  - **(A0), (M0) or (A1) red:** the fixture is wrong. Stop and mark the item BLOCKED.
  - **(T3), (H) or (i′) red:** park on `attempt/WF-9-step1-<ts>`, write a known-issues entry, and mark the item BLOCKED (rule (d)).
  - **(PF) red on P2 / P4 only, with (T3) green on all four drives:** the thin-wire transverse term on a transverse on-axis background. It is a finding about the 1-D fold, not the code. Report the error split, park, and mark the item BLOCKED. The weekly rules, and the band is never widened.
  - **(BF) not missing:** the lead is not a strong scatterer at this contrast. Report it, because it reverses the Born premise, and stop.
  - **A record outside its prediction:** report it; it has no status effect.
  - There is no `PREDICTED RED` line.

### §7 row text moved verbatim — 2026-10-04 (18:00 daily review, row over 2 KB)

Step 0 ✅ 2026-10-03: the exact identities (superposition, reciprocity, adjoint fold, held = fresh) all hold to ≤ 3.4e-12 at `-n 2` and `-n 4`. The path fold is recorded: D2's error is 9.7e-3 at N = 24, the thin-wire transverse floor (`20261003T003327Z_WF-9-step0-W1.log`, `…003546Z_WF-9-step0-W2.log`). *Step 0 (I7, daily 2026-10-02 — provisional, rule (6); §9 item 183): the row's identity in exact discrete form — a σ = 800 S/m lead on one structured mesh in `IMP-1`'s host, stepped-source TF = adjoint TF (reciprocity, Jin §6.4) and the adjoint fold of the lead-free E = the direct with-lead tip field, two drives, at `DISCRETE_IDENTITY_RTOL`; the path fold against N recorded, no band. W1 + W2 green ⇒ 🟡 only.*

### §9 item 276 — 2026-10-06 (chain I16)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: PROJECT_PLAN §10 chain I16 (grep "I16. `WF-9` step 2"), whose target is the identities. **No band, tolerance, Done-when or validation target is set or moved.** Every band below is imported by name from the file:line cited.

**`WF-9` step 2: step 1's fixture at degree 2.** This is the same kept-rod F-small mesh (182 371 cells, `KEPT_ROD_CELL_RECORD`), the same 64 MHz, four 50 Ω sheets, σ = 800 S/m rod, tip functional and 32 slabs (§9 item 224, `tests/validation/test_wf9_lead_in_coil_phantom.py`). The only change is that every solve is in degree-2 Nédélec space. The step re-asserts step 1's identities and prints how much the TF moves between orders on each drive.

- **Executor:** implementer. **Tier:** W0 standard; W1 and W2 **heavy**. **Ranks:** `-n 2` for W0, `-n 8` for W1, `-n 12` for W2. **Cost:** ≈ 19 measured (chain-solve median) / **≈ 50 predicted** min. That is 15 fixed + ≈ 10 code + W0 3.2 + sweep 0.5 + W1 ≈ 8 + W2 ≈ 7 + ≈ 6 spare. **Park if W1 is not green by minute 40** (rule (d)).
- **Depends on:** nothing open.
- **Why not `-n 2` at degree 2.** Degree-1 dofs are 216 731 (`20261004T201655Z_WF-9-step1-W1.log:1895`). Degree-2 first-kind Nédélec carries 2 dofs per edge + 2 per face, so ≈ 5.3× that, ≈ 1.15 M. The closest measured degree-2 solve on the same birdcage family is `ANS-4` item 225 W3 (`20261004T183338Z_ANS-4-item225-W3-p2-cf0p5-128-koff.log`):
  - 259 524 cells, ≈ 1.6 M dofs by the same ratio;
  - one factorisation + three held solves in 138.0 s at `-n 8` (`:2083`);
  - 36.28 GiB summed RSS (`:2102`);
  - a 261 s window (`:2712`).

  Scale by cells (0.703) and take factor work ∝ n² and memory ∝ n^(4/3). This step needs three factorisations: the without arm, the with arm, and the fresh (H) solve. At `-n 8` that is ≈ 50–70 s each, ≈ 22–30 GiB summed, and one solver alive at a time (`_free`). At `-n 2` the three factorisations would cost ≈ 3–4× more and risk the 590 s ceiling. Rank-safety of the module's reductions is already proved at `-n 2` / `-n 4` (step 1), and W0 re-proves it at `-n 2`.

**Code (additive; names are proposals).** All changes go in `tests/validation/test_wf9_lead_in_coil_phantom.py` only. No `src/` change.
- (1) Add `DEGREE_ENV = "FEM_EM_WF9_DEGREE"` with values in {1, 2}. Unset means 1, so step 1's nodes behave byte for byte as before. Add a `_degree()` reader on the pattern of `test_ans4_exact_strip_records.py:113, 336`, which raises on any other value.
- (2) Thread the degree through the module:
  - `build_four_port_sweep(..., degree=_degree())` at `:124`; the signature is at `test_port_birdcage_four_port.py:219–230`. `sweep["degree"]` then carries it, so the fresh (H) solve `_solve_driven(...)` at `:391` takes it without an argument (`test_birdcage_b1_plus_map.py`, the `degree=None` branch).
  - `TimeHarmonicSolver(problem, degree=...)` at `:288` (`_arm`) and `:327` (`_symmetry`).
- (3) **TF exchange.**
  - At degree 1, keep the file name `wf9_step1_tf_n{size}.npy` and the partner loop `(2, 4)` (`:443–447`).
  - At degree 2, write `wf9_step2_tf_p2_n{size}.npy` with partners `(8, 12)`. Print the cross-width TF agreement, predicted ≤ 1e-10.
- (4) **Order move** (printed only). Add `STEP1_Q_RECORD`, holding q_without / q_with per drive at `.10e`, copied from `20261004T201655Z_WF-9-step1-W1.log:1900, 1913, 1926, 1939` with that citation. Print per drive:
  - |Δq₂/Δq₁ − 1|;
  - arg(Δq₂/Δq₁) in degrees;
  - |q_with,2/q_with,1 − 1|.

  If `wf9_step1_tf_n2.npy` exists (W0 writes it fresh this slot), also print max_k |TF₂,k − TF₁,k| / max |TF₁|. Otherwise print "absent", which is not an error.
- (5) **New node `test_wf9_step2_degree2`.** It is `@complex_only` and skips unless `_degree() == 2`. It calls the same identities body as `test_wf9_step1_identities`, factored into a shared helper with the degree passed in, so no logic is duplicated. It adds (Deg) and (4).

**Anchors (asserted; bands imported, none restated).**
- **(Deg)** The solver's space reads back as degree 2 on every factor, and global dofs > 4 × 216 731. Predicted ≈ 1.1–1.2 M. This guards against the silent degree-1 path (trap 1).
- **(A0) / (M0)** Step 1's census, unchanged: the cells equal `KEPT_ROD_CELL_RECORD` within `CELL_COUNT_BAND` (`test_birdcage_port_sheet_prerequisite.py`), and DG0 σ / εᵣ on the implant read exactly 800 / gel.
- **(A1)** ‖A − Aᵀ‖_F ≤ `SYMMETRY_RTOL` (1e-10, `test_wf9_transfer_function_identity.py:59`) · ‖A‖_F on the degree-2 with-arm form plus the four sheet terms. Predicted ≈ 1e-16.
- **(i′)** ‖Y − Yᵀ‖/‖Y‖ ≤ `RECIPROCITY_BAND` (1e-3, `test_port_lumped_sheet_sweep.py:86`), both arms. Predicted ≤ 1e-14.
- **(T3) The adjoint fold identity, per drive.** |−jωμ₀Δσ∫χ_lead w·E_without,d − Δq_d| ≤ `DISCRETE_IDENTITY_RTOL` (1e-6, `test_birdcage_power_identity.py:227`) · |Δq_d|. Predicted ≤ 5e-11: step 1 measured ≤ 3.7e-12, and degree 2 has a worse condition number.
- **(H)** The fresh q_with,P1 equals the held value within `DISCRETE_IDENTITY_RTOL`. Predicted ≤ 1e-11.
- **(PF)** |Q₃₂ − Δq_d|/|Δq_d| ≤ `PATH_FOLD_BAND` (0.05, `test_wf9_lead_in_coil_phantom.py:75`) on every drive. Predicted ≈ 1e-4–3e-3 (step 1: 1.44e-4 … 1.47e-3). The fold's error is set by geometry, through the transverse term and the within-slab covariance, not by order.
- **(BF)** The Born fold misses: |B_d − Δq_d|/|Δq_d| > `PATH_FOLD_BAND` on every drive. Predicted [100 %, 3000 %] (step 1: 1 344–1 435 %).

**Negative controls (rule (e)).**
- (BF) is the Done-when's own control. It is **asserted**, its factor is **predicted**, and it has no ceiling.
- (C1), the Hermitian read-back with conj(w), and (C2), the total-field fold with E_with, are **asserted only as wired** (> `DISCRETE_IDENTITY_RTOL`). Their factors are **predicted**: C1 in [10 %, 200 %] (ceiling 200 %; step 1 read 199.5 %) and C2 in [80 %, 105 %] (step 1 read 93.5 %).
- (Deg) is **asserted**. It is a planted-path guard: at degree 1 its dof count would read 216 731, below the 4× bar.

**Printed (records, not asserted).**
- The order move (4), per drive. Predicted |Δq₂/Δq₁ − 1| in [0.5 %, 30 %]. The 1 mm rod is resolved by a few degree-1 cells across, so the tip ball's response is where the order should show most.
- |q_with/q_without| per drive.
- The TF profile (|TF_k|, arg TF_k).
- Q_N for N = 4…32 with the error split.
- The W1/W2 TF agreement.
- Dofs, phase times and summed RSS. Predicted 20–32 GiB at `-n 8`.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`, `-s`, node IDs only, no pipe; each window is < 10 min, so no durable capture is needed).
- **W0, step 1 re-gate at degree 1** (rule (c): the module is edited). Standard. Predicted ≈ 190 s (step 1's W1 was 190 s). It must reproduce step 1's printed q values at `.10e` within `L0_BAND`, and it writes `wf9_step1_tf_n2.npy`:
  `scripts/testing/run_and_log.sh WF-9-step2-W0-deg1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf9_lead_in_coil_phantom.py::test_wf9_step1_census tests/validation/test_wf9_lead_in_coil_phantom.py::test_wf9_step1_identities -v -s --tb=short'"`
- **Orphan sweep before W1.** The degree-2 forms compile cold. `docker compose exec -T fem-em-solver pgrep -c python3` must read 0, then run the rubric's one-line any-orphan-`.c` sweep under tag `WF-9-step2-cache-sweep`.
- **W1, degree 2, `-n 8`.** Heavy. Predicted 320–480 s:
  - environment ≈ 40 s;
  - mesh ≈ 50 s;
  - cold degree-2 JIT ≈ 30–90 s;
  - A1 assembly ≈ 25 s;
  - three factorisations ≈ 150–210 s;
  - nine held solves ≈ 20 s.

  `scripts/testing/run_and_log.sh WF-9-step2-W1-p2-n8 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF9_DEGREE=2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_wf9_lead_in_coil_phantom.py::test_wf9_step1_census tests/validation/test_wf9_lead_in_coil_phantom.py::test_wf9_step2_degree2 -v -s --tb=short'"`
- **Gate.** If W1 is > 540 s or Status 124, record the phase times and stop BLOCKED; do not run W2. **Never lengthen the 590 s timeout or raise `-n` above 12.**
- **W2, degree 2, `-n 12`.** Heavy. Predicted 280–450 s. It runs the same command with label `WF-9-step2-W2-p2-n12` and `-n 12`, and reads the W1/W2 TF agreement.

**Traps already paid for.**
1. **Silent degree 1.** `TimeHarmonicSolver(problem)` and `build_four_port_sweep()` both default to degree 1 (`WF-15` step 1; `PORT-21` step 3). There are three sites to change (`:124`, `:288`, `:327`), and `:391` follows `sweep["degree"]`. (Deg) is the guard.
2. **Stale exchange files.** `/workspace/logs/wf9_step1_tf_n{2,4}.npy` exist from 2026-10-04. A degree-2 run reading them would print a meaningless "agreement", so the degree-2 file name and partner widths are distinct, per (3).
3. **No module edit after W0** (rule (i)). W0, W1 and W2 run on one module state.
4. **Unmapped kept cells become vacuum** (`test_port_birdcage_four_port.py:254`). (M0) stays first.
5. **Use `ufl.dot` or components, never `ufl.inner`**, for q / TF / folds. Reduce complex assembles.
6. **Masks are DG0 from numpy midpoints, owned + ghost.** Counts are over owned cells only.
7. **Orphan `.c`.** Sweep only with `pgrep -c python3` = 0. After any kill, run `pgrep -c python3` first and check for orphaned ranks before anything else.
8. **Status 124 is a finding, not a retry.**
9. Write ≈, never `~`. No AED file is touched.

**Scope.**
- **Covered:** degree 2 on step 1's single fixture: one mesh, 64 MHz, untuned solid-coil F-small, 50 Ω sheets, one σ = 800 S/m capped rod, a linear tip functional at one tip, at `-n 8` / `-n 12`.
- **Not claimed:**
  - h-convergence or an order-converged TF (two orders on one mesh is not a convergence rate);
  - a measured or ISO 10974-calibrated TF;
  - heating or SAR;
  - a PEC, insulated or real lead;
  - quadrature drives;
  - F-human;
  - Palace.

**Status it can move.** No glyph: `WF-9` stays ✅ "as a fixture". Its Done-when is step 1's (I8) and it does not name degree 2. On green, the §7 row gains a **step-2 record** sentence giving:
- (T3) / (H) / (PF) / (BF) maxima at degree 2;
- the order move per drive;
- dofs and elapsed.

§10 I16 is then marked DONE.

**Negative result.** In every case, report the measurement, park on `attempt/WF-9-step2-*`, mark item 276 BLOCKED in the same commit (rule (d)), and journal it.
- **W0 red:** the edit broke step 1. Revert and stop.
- **(Deg), (A0), (M0) or (A1) red:** the fixture or the threading is wrong. Stop.
- **(T3), (H) or (i′) red at degree 2 with W0 green:** a degree-2 discrete-identity failure (sheet terms or held termination at p = 2). Write a known-issues entry and stop.
- **(PF) red with (T3) green:** a finding about the 1-D fold at degree 2. Report the error split; the weekly rules, and the band is never widened.
- **(BF) not missing:** report it and stop.
- **Gate or 590 s kill:** BLOCKED without a result; the review re-poses (e.g. `-n 12` only, or a split census / identities window).
- **A record outside its prediction** has no status effect.

No `PREDICTED RED:` line.

**Outcome (2026-10-06 21:00 CDT slot, 2026-10-07 UTC): ✅ GREEN.** W0 (degree-1 re-gate, `-n 2`, `20261007T020147Z_WF-9-step2-W0-deg1.log`, Status 0, 196 s) reproduced all eight step-1 q values at `.10e`, 216 731 dofs. Cache sweep `20261007T020512Z_WF-9-step2-cache-sweep.log` (no orphans). W1 (`-n 8`, `20261007T020646Z_WF-9-step2-W1-p2-n8.log`, Status 0, 512 s, under the 540 s gate) and W2 (`-n 12`, `20261007T021530Z_WF-9-step2-W2-p2-n12.log`, Status 0, 425 s), 14 passed each:
- (Deg) four factors at degree 2, 1 166 260 dofs (5.381 × degree 1; bar 866 924); (A0) 182 371 cells, rel 0; (M0) exact.
- (A1) 2.7e-16; (i′) ≤ 2.19e-14 both arms; (H) ≤ 6.85e-12; (T3) max 1.922e-11 (P2, W1).
- (PF) P1 1.473e-4, P2 1.550e-3, P3 1.419e-4, P4 9.078e-4 (band 0.05). (BF) 1 319.6–1 464.5 % (misses, as asserted). C1 199.56 %, C2 94.03–94.04 % (wired).
- Printed: order move |Δq₂/Δq₁ − 1| P1 13.15 %, P2 28.86 %, P3 2.63 %, P4 29.10 % (all inside [0.5 %, 30 %]); arg −0.63° / +0.56° / −1.46° / +0.44°; |q_with,2/q_with,1 − 1| 10.8 / 27.3 / 3.8 / 27.3 %; TF profile move max_k |TF₂ − TF₁| / max |TF₁| = 9.06e-2; |q_with/q_without| 2.22–2.31; W1/W2 TF agreement 2.118e-11; summed RSS 30.95 GiB (`-n 8`), 36.76 GiB (`-n 12`). The module's older step-1 per-drive PF prediction text prints "OUT" for P2 / P4 at both orders (record only). W2's RSS line reuses the "-n 8" prediction label (cosmetic; module frozen after W0 by rule (i)).

## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 **Implant transfer-function excitation (ISO 10974 Tier 3)** — a local impressed source stepped along a lead path + the tangential-`E` path-integral post-processor; heating from `∫ TF(z)·E_tan(z) dz` gated against the direct lead-in-phantom solve — **feature ladder C1** (operator directive 2026-09-04; Phase 7; `WF-7` is to be scoped with this in mind) | ✅ **as a fixture, 2026-10-04 (step 1, §9 item 224, chain I8; audit PASS, 18:00 daily review):** the σ = 800 S/m capped rod kept in the coil-driven F-small phantom at 64 MHz (`GEO-36` step 3 mesh, 182 371 cells), switched by material only. On all four port drives the adjoint fold equals the direct with-lead tip response to ≤ 3.7e-12 (`DISCRETE_IDENTITY_RTOL` 1e-6), the N = 32 path fold to ≤ 1.47e-3 (band 5 %), and the Born fold misses by 1 344–1 435 %; A = Aᵀ 9.8e-17, loaded reciprocity ≤ 1.8e-15 both arms, held = fresh ≤ 2.7e-14, C1 / C2 wired (199.5 % / 93.5 %); W1/W2 TF agreement 2.9e-12 (`20261004T201655Z_WF-9-step1-W1.log` `-n 2` 190 s, `…202017Z_WF-9-step1-W2.log` `-n 4` 165 s). Not claimed: a measured TF, heating, a real lead, h-convergence. **Step-2 record (degree 2, same mesh, §9 item 276, chain I16, 2026-10-06 21:00 slot):** 1 166 260 dofs (5.38 × degree 1); (T3) ≤ 1.92e-11, (H) ≤ 6.8e-12, (PF) ≤ 1.55e-3, (BF) 1 320–1 465 %, A1 2.7e-16, (i′) ≤ 2.2e-14; order move |Δq₂/Δq₁ − 1| = 13.2 / 28.9 / 2.6 / 29.1 % (P1–P4), TF profile move 9.1e-2; W1/W2 TF agreement 2.1e-11 (`20261007T020646Z_WF-9-step2-W1-p2-n8.log` 512 s, `…021530Z_WF-9-step2-W2-p2-n12.log` 425 s; W0 degree-1 re-gate `…020147Z_WF-9-step2-W0-deg1.log` 196 s). Two orders on one mesh — not a convergence claim. Step 0 ✅ 2026-10-03 (history in the chunk file). **Done-when (written by the 2026-10-03 weekly, §10 rulings; copied by the 10-04 03:00 daily):** *step 1 (I8)* — the lead in the coil-driven F-small phantom at 64 MHz, on one mesh with the lead's cells kept and switched by material (σ = 800 S/m against gel), so the adjoint fold is an exact discrete identity; asserted: for each of the four port drives the adjoint TF folded with the lead-free field equals the direct with-lead tip response at `DISCRETE_IDENTITY_RTOL`, the 1-D path fold ∫TF(z)·E_tan(z) dz reproduces it within **5 %** on every drive, and the lead-free Born fold misses ⇒ ✅ as a fixture. Not claimed: a measured transfer function, heating, a real lead. *History: `docs/planning/chunks/WF-9.md`.* 

### §9 item 317 — 2026-10-08 (chain I18b — the insulated lead's transfer function)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 I18b (`PROJECT_PLAN.md:4699–4704`), weekly-enumerated, serial on I18a. **Depends on item 316; if 316 did not land, skip.** Check: `grep -n 'BIRDCAGE_IMPLANT_COAT_CELL_TAG' src/fem_em_solver/io/mesh.py` is non-empty and 316's W2 is green in `test-results.md`. No band, tolerance, record or Done-when moves.

- **Executor:** implementer. **Tier:** heavy (the weekly's 190 s class), complex, degree 1, `-n 2` (step 1's width, where a missing reduction shows), then `-n 4`. **Cost:** ≈ 14 measured (`WF-9` step 0) / ≈ 55 predicted min. **Depends on:** item 316.
- **Fixture.** Item 316's coated build: the rod is tag 403, the sleeve tag 404, the gel tag 3. The F-small four-port sweep runs at `STEP2D_FREQUENCY_HZ` 64 MHz through `build_four_port_sweep(..., build_only=True, degree=1)`, as `_build_case` does (`tests/validation/test_wf9_lead_in_coil_phantom.py:151–219`). There are three material states on **one** mesh, every tag set explicitly in each `material_map`. A tag missing from the map must raise, never default.
  - **without:** 403 and 404 are gel. This is the lead-free field, the ISO 10974 background.
  - **bare:** 403 is σ = `SIGMA_WIRE_S_PER_M` 800 with the gel's εᵣ (step 1's `lead_mat`, `:168–170`); 404 is gel.
  - **coated:** 403 is as in bare; 404 is εᵣ 3, σ 0, μᵣ 1 (`IMP-5.md:29`).
- **The anchor: the adjoint fold with a complex, per-cell contrast.** In step 1 the contrast is `LOAD·Δσ` on χ_lead (`:450,459`; `LOAD` = −jωμ₀ at `:78`). With the sleeve, the lead is rod ∪ coat, and the operator difference A_coated − A_without has a σ term on 403 and an εᵣ **and** σ term on 404. So the fold weight becomes the DG0 field c = `LOAD`·Δσ + ω²μ₀ε₀·Δεᵣ, where Δσ and Δεᵣ are read cell by cell as the difference of two `build_material_fields` calls (`:229–250`, `core/time_harmonic.py`).
  - **Derive the sign from the module's own operator, never by hand.** This is the identity check (T3c) below. With εᵣ ≡ gel on 403 the ε term vanishes there, and c reduces to step 1's `LOAD·Δσ` cell for cell. This is asserted, in (R1).
  - The bare arm keeps step 1's real Δσ on 403 only.
- **New module `tests/validation/test_wf9_insulated_lead.py`** (step 1's module is **not** edited, so there is no W0 re-gate). It imports by name: `_arm`, `_symmetry`, `_scalar`, `_material_minmax`, `PATH_FOLD_BAND` (`:75`), `N_SLABS`, `SLAB_Z0`, `SLAB_DZ`, `FOLD_NS`, `LOAD`, `OMEGA` and `PORT_IDS` from `test_wf9_lead_in_coil_phantom`; `DISCRETE_IDENTITY_RTOL` from `test_birdcage_power_identity.py:227`; and the coated builder from item 316's module. Its own case dict carries every key that `_arm` (`:300`) and `_symmetry` (`:360`) read. If `_arm` cannot take a DG0 contrast, write a local fold and keep `_arm` for the solves. The ball, q and the slabs follow `:181–185` and `:409–414`; the slabs now cover rod ∪ coat cells.
- **Anchors (asserted, per drive P1–P4, for the coated arm and again for the bare arm).**
  - (A0) The cell count equals item 316's printed record at `CELL_COUNT_BAND`. 404 has > 0 cells, and the ball is disjoint from 403 ∪ 404.
  - (M0) The read-back is exact on every state: σ / εᵣ minmax on 403 and on 404.
  - (R1) c on 403 equals `LOAD·Δσ` exactly.
  - (A1) A = Aᵀ for each state at `DISCRETE_IDENTITY_RTOL`.
  - (T3c) \|Σ c·w·E_without − Δq\| ≤ `DISCRETE_IDENTITY_RTOL`·\|Δq\|. This is step 1's Done-when identity on the coated lead; predicted ≤ 1e-11.
  - (H) Held equals fresh on the coated P1.
  - (PF) The N = 32 path fold is within `PATH_FOLD_BAND`.
  - (BF) The Born fold misses by more than `PATH_FOLD_BAND`.
- **Printed, never gated (rule (e)).**
  - \|TF_coated(z)\| ÷ \|TF_bare(z)\| per slab, and the shaft-to-tip TF ratio for each arm. The prediction (not asserted) is that the sleeve suppresses the shaft TF and concentrates it at the bare tips.
  - \|q_with/q_without\| for coated and bare, beside step 1's 2.29–2.39 (PRED_RATIO, `:130`). This is read as the weekly's "tip ratio" (see Q2).
  - The bare arm's Δq beside step 1's `STEP1_Q_RECORD` (`:100`). The mesh differs, so this is printed only.
  - (PF) on the coated lead, and the fraction of the fold carried by 404.
  - Per-phase times and summed RSS.
- **Negative controls (rule (e)).**
  - (C3c, **asserted** only as wired, > `DISCRETE_IDENTITY_RTOL`) The rod-only contrast (drop 404's term) misses T3c. The factor is **predicted** and printed. It has no prior measurement, and its ceiling is unknown because the coat's share of Δq is the quantity being read.
  - (C1) Hermitian conj(w) and (C2) total-field E_with, carried from step 1 (`:461–462`): wired, factors predicted in [10 %, 200 %] and [80 %, 105 %].
- **Windows** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; no `-k`; no pipes; complex).
  - Price basis: step 1 ran two arms plus a fresh solve in W1 190 s at `-n 2` (`20261004T201655Z_WF-9-step1-W1.log`) and W2 165 s at `-n 4`. Item 316 adds a third arm and a predicted 2–18 % in cells. Predicted W1 ≈ 260 s, range 200–380; W2 ≈ 220 s.
  - Run `tests/environment/test_complex_mode.py` first, as part of `tests/environment`.
  - **W1:**
    `scripts/testing/run_and_log.sh WF-9-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 480 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_wf9_insulated_lead.py -v -s --tb=short'"`
  - **W2:** the same with label `WF-9-step3-W2` and `-n 4`; it reads the W1 / W2 TF agreement.
  - **Gate:** if W1 runs > 450 s or ends with Status 124, record the phase times and stop BLOCKED. Never lengthen the timeout.
- **Traps already paid for.**
  - Read-backs use components or `ufl.dot`, never `ufl.inner` (step 1's C1 is exactly that conjugation bug).
  - No `ufl.max_value` and no ordering on complex operands (`OPS-22`).
  - `assemble_scalar` is rank-local, so reduce through `_scalar`.
  - If a window is killed, run `pgrep -c python3` first, then the rubric's cache sweep only if it reads 0.
  - Only one factorised solver is alive at a time (`_free`, step 2).
  - Rule (i): no edit after W1 starts.
- **Scope.** Not in this item:
  - a measured or ISO-calibrated TF, heating, SAR, or 1 g SAR;
  - a real lead (the coat is order-of-magnitude silicone, `IMP-5.md:29`);
  - degree 2 or h-convergence;
  - the `IMP-5` series.

  The coat-versus-bare difference is a record on one mesh, not a claim. `IMP-5` and `GEO-36` do not move.
- **Negative result.** A red on (T3c), (A1) or (H) means the contrast has the wrong sign or misses a term. Print c's per-tag means and C3c, park on `attempt/WF-9-step3-<ts>`, and mark BLOCKED. A (PF) red on the coated lead is a finding: the 1-D fold does not survive a dielectric sleeve at N = 32. Report every drive and keep it as a record; never widen `PATH_FOLD_BAND`. A (BF) green-as-miss failure (the Born fold lands inside the band) is reported verbatim. No `PREDICTED RED` line.
- **Status it can move.** The §7 `WF-9` row stays ✅ as a fixture. Its Done-when (`PROJECT_PLAN.md:2903`) covers the bare σ = 800 lead and is not extended by this item. The row gains a step-3 record: the identities on the insulated lead, plus the coated / bare TF profile and tip ratios. This unblocks I19 (`imp:7` coated-lead rows).

**Review rulings (2026-10-08 03:00 daily review, item 317).** (2) The tip readout is `WF-9`'s ball-averaged E_z ratio |q_with / q_without| for the coated and bare arms, printed; no 1 g SAR is computed — the weekly's "tip 1 g ratio" is read in `WF-9`'s own quantity and the reading is named for the 10-10 weekly. (3) **Accepted:** the without-lead arm is all gel (the lead-free background, as in step 1), and the whole coated lead is the perturbation with the per-cell complex contrast; step 1's module is not edited. (6) The cost and the 450 s W1 stop are accepted; a kill at 590 s parks the item. Skip if item 316 did not land.

**Outcome (2026-10-08 09:00 slot, take-next): ✅ green first run, both windows.** New `tests/validation/test_wf9_insulated_lead.py`; c = k₀²(ε_c,state − ε_c,without) read from `bilinear_form`'s own ε_c = εᵣ − jσ/(ωε₀) (mean c on 404 = −134.94 + 252.66j m⁻², on 403 −4.040063e5j = LOAD·Δσ). (A0) 190 747 = record; (M0) exact; (R1) 1.44e-16 (rtol 1e-13, the ulp floor of k₀²/(ωε₀) vs ωμ₀), c ≡ 0 off the lead; (A1) ≤ 9.86e-17; (H) 0; (T3c) coated 1.3e-13 … 1.9e-12, bare 1.8e-14 … 6.4e-13; (PF) coated 1.22e-3 / 2.54e-3 / 2.87e-3 / 1.88e-3 (P1–P4), bare 3.09e-4 / 1.54e-3 / 8.68e-4 / 1.72e-3; (BF) coated 2 713 / 2 870 / 2 729 / 2 578 %, bare 1 575 / 1 667 / 1 578 / 1 493 %; (C3c) 1.216 / 1.438 / 1.369 / 0.864 % (= the coat's share of the fold). Printed: C1 coated 200.5–200.8 % (just OUT of [10 %, 200 %]), C2 coated 363–371 % (OUT of [80 %, 105 %] — the coated lead's total field differs far more from the lead-free field than the bare one's), bare 199.5 % / 93.6–93.7 % (IN); \|q_with/q_without\| coated 1.901 / 1.865 / 1.955 / 1.887, bare 2.516 / 2.456 / 2.611 / 2.494; bare Δq / step-1 record 1.11–1.15 (other mesh); \|TF_coated/TF_bare\| per slab 3.76 (z = −15.5 mm) → 1.0 (−8.5) → 0.36 (+12.5) → 0.92 (+15.5); shaft (\|z\| < 12 mm) / end-slab \|TF\| 3.04 coated vs 8.67 bare — the sleeve flattens the profile toward the bare tips, as predicted; W1/W2 TF agreement 4.16e-11. Times (`-n 2`): three factorisations 12.3–12.7 s each, fresh 12.7 s, A1×3 6.3 s, body 104 s, summed RSS 3.19 GiB. `20261008T142909Z_WF-9-step3-W1.log` (233 s), `20261008T143309Z_WF-9-step3-W2.log` (205 s, `-n 4`, RSS 3.89 GiB). A record on one mesh, degree 1; unblocks I19.

**Review rulings (2026-10-08 10:30 daily review, item 317 audit PASS).** (R1-tol) The module asserts (R1) at `R1_RTOL` = 1e-13 (`tests/validation/test_wf9_insulated_lead.py:98–100, :472`) where the brief said "exactly"; measured 1.441e-16 (`20261008T142909Z_WF-9-step3-W1.log:1960–1961`). **Accepted as the floating-point rendering of "exactly"**, not a loosened band: no prior bound existed, the constant predates W1 (mtime), the run was green first time, and (R1) is a sign-and-term check whose failure modes miss by order one, not the §4 anchor (that is (T3c) at `DISCRETE_IDENTITY_RTOL`). It sets no precedent; the weekly may tighten it toward the few-ulp level the code comment names. (C1 / C2) The coated arm's C1 / C2 controls fell outside step 1's predictions (W1 `:2024–2060`): rule (e) records, named for the 10-10 weekly.
