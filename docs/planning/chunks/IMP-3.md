# IMP-3 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 174 — 2026-10-02 (chain I6; queued by the 10:30 review)

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 10:30 daily review. Chain I6 is a weekly-enumerated step (§9 rule (4)), and its serial prerequisite `GEO-36` landed ✅ in `73f65d05` (item 171). The review amends the §7 Done-when under rule (j) in this commit: W1 ⇒ 🟡, W2 ⇒ ✅. That **tightens** the row, because the record mesh cannot carry the 1 g ratio. The weekly ratifies.*

`IMP-3` step 1 compares the loaded F-small birdcage at 64 MHz with and without a capped PEC rod in the phantom.
- **Executor:** implementer. **Tier:** heavy, complex: W1 at `-n 2`, W2 at `-n 4`.
- **Cost:** ≈ 11 measured / ≈ 55 predicted min.
- **Depends on:** nothing.

**Two findings drive the design:**
1. **The flag-off record mesh has 537 phantom cells** (`20261002T110559Z_GEO-36-step2-W2.log:3611`). A 6.2 mm 1 g ball therefore holds less than one cell, so the ratio of record needs a second window on the 2.5 mm phantom rung.
2. **The phantom absorbs only 3.0e-4 of the source power at 64 MHz.** The rod moves the 4×4 below the mesh-noise floor, so the asserted mirror check is a consistency identity. The readings that do separate the rod's mirror are the tip-SAR pairs, and those are printed.

**Fixture and geometry:**
- **Rod.** `GEO-36`'s rod: r = 1 mm, ℓ = 30 mm between cap centres (32 mm tip to tip), axis ẑ, centre (0.015, 0, 0), on leg 1's azimuth (`tests/mesh/test_birdcage_implant_rod.py:75–79`). It is built by `_build_rod(**overrides)` (`:97`).
- **Without the rod:** `tests/mesh/test_birdcage_port_sheets._build(True, …)` (`:80`).
- **Entry point.** Both meshes enter `build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=<tuple>, build_only=True)` (`tests/validation/test_port_birdcage_four_port.py:219`, keyword doc `:234–241`, build-only return `:516–529`). `_build_rod()` returns exactly the `(msh, cell_tags, facet_tags, diag, t)` shape it needs. `STEP2D_FREQUENCY_HZ` = 64e6 (`test_birdcage_power_identity.py:182`).
- **Coil.** PORT-11 / PORT-16's coil: solid σ = 800 S/m (`SIGMA_WIRE_S_PER_M`), 50 Ω sheets, untuned at 64 MHz. `ANS-7`'s tuned copper coil is out of reach, because `implant_rod` refuses `as_hole`.

**Code: tests only, no `src/`, no shared-helper edit.**
- New file `tests/validation/test_imp3_rod_birdcage.py`, with `test_imp3_default_phantom` (W1) and `test_imp3_one_gram_rung` (W2).
- Both call a module-level `_case(phantom_resolution)`, which runs these steps in order:
  1. **Sweep the without-rod arm first.** Call `run_n_port_sparameter_sweep(problem, port_defs, lumped_sheet_ports=specs, lumped_sheet_facet_tags=facet_tags, keep_fields=True)` (pattern at `test_port_drive_superposition.py:317–323`): one factorisation plus three held solves.
  2. **Matrix quantities.** Z, S, σ(S), column powers and ‖S − Sᵀ‖/‖S‖, by `test_port_birdcage_four_port.py:543–557`'s arithmetic. Classes and spreads come from `_circulant_classes` (`:193`) and `_class_spread` (`:176`).
  3. **Power partition per drive.** `_exact_shares(case, {"driven": pid, "omega": ω, "fields": result.fields[pid]})` (`test_birdcage_power_identity.py:313`).
  4. **Quadrature drives.** `superpose_drives(result, w)` (`src/fem_em_solver/ports/superposition.py:118`), with weights from `quadrature_phase_weights(ks, "ccw"/"cw")` (`test_birdcage_b1_quadrature.py:175`) and `ks` from `_port_index` (`:158`).
  5. **1 g SAR.** `mass_averaged_sar(e, sigma=fields.sigma_field, rho=rho_field, center=tip, radius=a_1g, quadrature_degree=QUADRATURE_DEGREE)` (`src/fem_em_solver/post/sar.py:231`).
     - ρ comes from `build_density_field` (`:147`), zero outside tag 3, as in `test_birdcage_sar_mass_averaged.py:251–260`.
     - a_1g = `averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)` = 6.204 mm.
     - The tips are at (0.015, 0, ±0.016) m. Assert containment: 15 + 6.2 < 30 mm and 16 + 6.2 < 40 mm.
  6. **Background field, without-rod arm only.** Sample E at 9 points on the rod axis, z ∈ [−16, 16] mm, through `evaluate_vector_field_parallel` on each drive's `e_real` / `e_imag`, and combine the superposed drives in numpy.
  7. **Free the without-rod arm.** Store its scalars, then `del` the result, fields and solver.
  8. **Rod arm.** Build it and repeat steps 1–5.
  9. **PEC dofs on the rod wall** come from the **build tuple's** facet tags, tag `BIRDCAGE_IMPLANT_SURFACE_TAG` = 402 (`src/fem_em_solver/io/mesh.py:97`), via a local three-line copy of `_cavity_dofs` (`test_imp1_pec_sphere_lossy_host.py:191`).
  10. **Closed form.** `pec_prolate_spheroid_field(p − ROD_CENTRE, semi_major=0.016, semi_minor=0.001, e0=1.0)` (`src/fem_em_solver/utils/analytical.py:1576`), for an equivalent spheroid of aspect 16.
      - Evaluate it only where the FEM has cells.
      - Integrate it by a numpy basix degree-8 rule over owned cells, with a ball indicator at the quadrature points, in chunks of 5 000 cells. This is a local `_ball_mean_closed_form` modelled on `test_imp2_spheroid_fem.py:181–201`, with χ_ball added.
  11. **Print before asserting.** Every reading is printed before the first assert. W1 also runs the control.

**The mirror symmetry the rod keeps:**
- **σ_y (y → −y), kept.** This is the plane through legs 1 and 3 and the axis. On the 4×4 it swaps P2 (90°) ↔ P4 (270°) and fixes P1 (360°) and P3 (180°) (sheet azimuths: `20260912T123236Z_PORT-14-step2d-64mhz.log.gz:1861–1864`).
- **σ_z, kept.** The rod is centred at z = 0 and the gaps are centred (`io/mesh.py:3337`). σ_z maps each port onto itself with its drive reversed, so it does not constrain S. It does force the top and bottom tips to agree.
- **C4 and σ_x (P1 ↔ P3), broken.**

**Anchors (asserted; every band imported):**
- **(A0) Census.** All within `CELL_COUNT_BAND` 0.01 (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`):
  - W1, without the rod: 116 085 = `SOLID_CELL_RECORD` (`tests/mesh/test_birdcage_conductor_hole.py:95`);
  - W1, with the rod: 167 371 (`20261002T110559Z_GEO-36-step2-W2.log:3612`);
  - W2, without the rod: 199 920 = `ONE_GRAM_RUNG_CELL_RECORD` (`test_birdcage_sar_1g_rung.py:125`);
  - W2's rod-arm count is a first record, predicted [225 k, 275 k]. **Guard: 320 000, checked before any solve.**
- **(R0) Reproduction, W1 without-rod arm.** Z reproduces the 64 MHz record entrywise within `REUSE_REPRODUCTION_RTOL` 1e-6 (`test_port_birdcage_ring_matrix.py:171`).
  - The record is a new module constant, `IMP3_FLAG_OFF_Z_64MHZ`, taken from `…123236Z_PORT-14-step2d-64mhz.log.gz:1866–1869`.
  - It is digit-identical in `20261001T171613Z_PORT-21-step5-W2.log:3670–3673`.
- **(i) Reciprocity, both arms:** ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND` 1e-3 (`test_port_lumped_sheet_sweep.py:86`). Without the rod it measured 9.95e-16 (`…123236Z…:1875`).
- **(ii) Passivity, both arms:** σ_max(S) ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` 1e-9 (`test_port_birdcage_four_port.py:129`), and each column's Σ|S|² ≤ 1. Without the rod: 0.999721388 (`:1876`).
- **(P) Power partition, every drive of both arms:** P_src,exact = P_vol + P_sheet,exact within `DISCRETE_IDENTITY_RTOL` 1e-6 (`test_birdcage_power_identity.py:227`). Measured 3.3e-15 to 7.0e-15 (`…123236Z…:1882–1891`).
- **(M) σ_y mirror, rod arm.** `_class_spread` of each of the five pairs (Z₂₂, Z₄₄), (Z₁₂, Z₁₄), (Z₂₁, Z₄₁), (Z₃₂, Z₃₄), (Z₂₃, Z₄₃) is ≤ `ADJACENT_SPREAD_BAND` 0.005 (`test_port_birdcage_lumped_column.py:121`).
  - Without the rod: 0.034 / 0.027 / 0.047 / 0.025 / 0.018 %.
  - Predicted with the rod: [0.01 %, 0.15 %].
- **(E) n × E = 0 on the rod wall, rod arm.** `pec_facet_tags=None` pins every exterior facet (`src/fem_em_solver/core/time_harmonic.py:365–367`).
  - The tag-402 dof count, reduced, is > 0.
  - Max |E| on those dofs is ≤ 1e-12 on every kept drive (`IMP-2` (iii)'s bound, `test_imp2_spheroid_fem.py:446`).
- **(S) Separation, without-rod arm:** the pooled off-diagonal spread is ≥ `POOLED_SEPARATION_FLOOR` 10× the worst class spread (`test_port_birdcage_four_port.py:173`).
  - Asserted, because the same comparison on the same fixture measured 671× (`20261001T171613Z_PORT-21-step5-W2.log:3682`).
  - Printed on the rod arm, predicted [300, 1000]×.

**Negative controls (rule (e)):**
- **(C1) The rod wall left natural.** W1, P1 only.
  - Build exterior facets minus 402, tagged `OUTER_BOX_TAG` 499 (`test_th14_birdcage_copper.py:114`), with a local copy of `_outer_box_tags` (`:133–157`), which subtracts 401 as written.
  - Assert first that 499 is absent from the mesh, and that the census reads exterior = outer + 402.
  - Run `_solve_driven(dict(case, problem=replace(problem, facet_tags=outer, pec_facet_tags=(OUTER_BOX_TAG,))), "P1")` (`test_birdcage_b1_plus_map.py:142`) on a fresh solver.
  - *Asserted:* (P) still closes at 1e-6 (it is a theorem), and the knob is wired: max |E| on the 402 dofs is > 1e-12.
  - *Predicted, printed:* the natural-wall tip 1 g over the without-rod tip 1 g in [0.7, 1.3], and PEC over natural in [2, 10]. A ratio of positive SARs has no ceiling, so no factor is claimed.
- **(C2) What the 4×4 sees of the rod. Predicted, printed.**
  - ΔP_phantom/P_phantom in [1 %, 15 %].
  - ΔΣ_i|S_i1|² in [5e-6, 1e-4], and max |ΔS| in [1e-6, 1e-3].
  - The σ_x pairs sit at the same mesh floor as the mirror pairs (without the rod: 0.051 %). This 4×4 cannot separate σ_y from C4; that is a scope fact.
- **(C3) Field-level mirror, both arms and windows. Printed.**
  - The tip 1 g pairs P2 / P4, ccw / cw and top / bottom: each predicted ≤ 1 %. `C4_COVARIANCE_BAND` 5 % is shown as a reference, never asserted.
  - P1 / P3, related only by the broken σ_x: predicted [50 %, 95 %], ceiling 100 %.

**Printed beside predictions written first.** None of these is asserted; their accuracy is `ANS-10`'s.
- **The ratio of record (W2; W1's coarse value beside it):** Q = SAR_1g,tip(rod) / SAR_1g,tip(without) at both tips, for P1–P4, ccw and cw. Predicted [2, 10].
- **Ball masses:** 1.000 g without the rod; predicted 0.9817 g with it.
- **The closed form:**
  - M_cf = ⟨|E_cf|²⟩ / E₀² on the same quadrature-masked region, predicted [2.5, 6];
  - L_z(aspect 16) = 0.009693 (`prolate_spheroid_depolarisation`, `analytical.py:1493`);
  - the axial fraction η_z, predicted [0.5, 1.0];
  - the estimate η_z·M_cf + (1 − η_z), with Q / estimate predicted in [0.6, 2.5].
- **The quasi-static scale:** (|k|·c_eq)² ≈ 7.4 % (|k| = 17.0 1/m, from `SALINE_SIGMA` 0.5 and `SALINE_EPSILON_R` 78).
- **Also:** P_phantom per drive, the narrowed facet counts per sheet, ΔS, cells, times and the summed peak RSS.

**Windows.** Bash `timeout: 660000`, foreground, each harness call bare and whole.
- **Basis:**
  - without the rod, 64 MHz, `-n 2`: mesh 21.5 s, four solves 7.0 s, module 93 s (`…123236Z…:1859,2084`);
  - the rod build: mesh 42.1 s, module 121 s (`GEO-36` W2:3612,3792);
  - the 0.0025 rung: mesh 34.2 s, sweep 75.4 s, module 251 s at `-n 4` (`20260908T140457Z_MAT-4-step5.log.gz:1903,1934,2290`).
- **W1** (the 🟡 window), predicted [150, 300] s:

  `scripts/testing/run_and_log.sh IMP-3-step1-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp3_rod_birdcage.py::test_imp3_default_phantom -v -s --tb=short'"`
- **W2** (the ✅ window), predicted [220, 480] s, with §5.1 durable capture and the trailing `; exit $rc`:

  `scripts/testing/run_and_log.sh IMP-3-step1-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_imp3_rod_birdcage.py::test_imp3_one_gram_rung -v -s --tb=short > /workspace/logs/imp3-step1-w2-raw.log 2>&1; rc=\$?; echo [capture] rc=\$rc >> /workspace/logs/imp3-step1-w2-raw.log; cat /workspace/logs/imp3-step1-w2-raw.log; exit \$rc'"`

**Traps already paid for:**
- **Build order.** Build without the rod first, then with it, in both windows: a failed gmsh build contaminates the process (`GEO-23`). The rod + 0.0025 combination has never been built; a FAIL is recorded verbatim and never retried.
- **`built=` ignores mesh-sizing keywords.** Pass `phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG` to `_build` / `_build_rod` themselves.
- **Tag 402 lives in the build tuple, not in the sweep's facet tags.** Keep that handle.
- **A held factorisation silently solves a different problem** (`PORT-19`, `WF-12`). The control uses `_solve_driven`'s fresh solver.
- **Point evaluation inside the hole raises** (`post/sar.py:357–361`). The axis samples therefore come from the without-rod arm only.
- **Rank-local quantities.** `assemble_scalar`, `cell_tags.values` and `locate_dofs_topological` are rank-local: reduce them, and count owned dofs only. Keep collectives outside any `rank == 0` branch.
- **No UFL comparisons, and no `sqrt` comparisons.** Do the ball mask and the closed form in numpy.
- **Free the without-rod arm before the rod build in W2.**
- **Harness discipline:**
  - always `-s`; node IDs only; no pytest `-k`; no pipes;
  - Status 124 is the finding: `pgrep -c python3`, then the cache sweep, and never lengthen the timeout;
  - no module edit after W1 starts (rule (i)).

**Scope:** degree 1, 64 MHz, the untuned solid-coil 4×4 with 50 Ω sheets, and one rod on leg 1's azimuth, axial. **Not covered:**
- the accuracy of the tip ratio (`ANS-10`);
- absolute or pointwise SAR, heating, compliance or convergence;
- the tuned or copper coil;
- F-human;
- other rod placements;
- ISO 10974 (`WF-9`) and ASTM F2182 (7c).

Aspect 16 lies beyond `IMP-2`'s FEM-checked aspects (2 and 5).

**Status it can move:**
- §7 `IMP-3` ⬜ → 🟡 on W1, then 🟡 → ✅ on W2, per the amended Done-when.
- §10 chain I6 is then done, and the weekly that reads it commissions `ANS-10`.

**Negative result:**
- **(A0) or (R0) misses:** the fixture moved. Stop and mark BLOCKED.
- **(i), (ii), (P), (M) or (E) red:** park on `attempt/IMP-3-step1-<ts>`, write a known-issues entry and mark BLOCKED (rule (d)). Never widen a band.
- **(M) red with the without-rod pairs green:** print the per-sheet facet counts, and file a known-issues entry on the rod grading perturbing the P2 / P4 sheet triangulations.
- **The C1 knob not wired:** a test bug. Fix it in-slot and re-run W1.
- **W2 over the cell guard, a gmsh FAIL, or Status 124:** W1 stands and `IMP-3` holds 🟡. Mark BLOCKED on the mesher.
- **A printed ratio outside its prediction:** report it; it has no status effect.
- There is no PREDICTED RED line.

**Outcome (2026-10-02, implementer slot 12:00 CDT): W1 and W2 green; `IMP-3` ⬜ → ✅.** New module `tests/validation/test_imp3_rod_birdcage.py` (mtime 17:04:04Z, before W1 at 17:04:28Z; unedited after).
- **W1** `20261002T170428Z_IMP-3-step1-W1.log`, `-n 2`, 241 s (13 passed; predicted [150, 300] s). (A0) 116 085 / 167 371 exact; (R0) max entrywise 2.98e-10 (band 1e-6); (i) 9.95e-16 / 1.20e-15; (ii) σ_max 0.999721388 / 0.999730355; (P) ≤ 6.97e-15 on all eight drives; (M) rod arm 0.0054 / 0.0373 / 0.0341 / 0.0542 / 0.0567 % (band 0.5 %); (E) 7 437 wall dofs, max |E| 0 on every drive; (S) 671.1× without the rod (486.6× with). C1: census exterior 8 272 = outer 3 314 + 402 4 958, tag 499 absent; (P) 2.16e-15; max |E| on the 402 dofs 1.29e-4 (knob wired); natural/without 1.085 / 0.913, PEC/natural 4.65 / 4.67 (both IN). Summed peak RSS 2.61 GiB.
- **W2** `20261002T170843Z_IMP-3-step1-W2.log`, `-n 4`, 260 s (13 passed; predicted [220, 480] s). (A0) 199 920 exact; the rod arm's first record is **250 866** cells (predicted [225 k, 275 k], guard 320 000). (i) 1.46e-15 / 1.51e-15; (ii) σ_max 0.999752421 / 0.999757728; (P) ≤ 4.98e-15; (M) rod arm 0.0772 / 0.0085 / 0.0428 / 0.0357 / 0.0819 %; (E) 7 437 dofs, max |E| 0; (S) 402.7× (456.6× with). Summed peak RSS 5.58 GiB.
- **Printed vs predicted (no status effect):**
  - Q on the 1 g rung, top / bottom: P1 4.496 / 4.498, P3 3.806 / 3.796, ccw 3.965 / 3.962, cw 3.959 / 3.958 (IN [2, 10]); **P2 1.580 / 1.576 and P4 1.586 / 1.582 OUT** (below 2). The P2/P4 drives put a mostly transverse field on the axis: η_z 0.405 / 0.443 (OUT of [0.5, 1.0]); P1 0.990, P3 0.902, ccw/cw 0.939. W1's coarse Q: P1 5.04 / 4.26, ccw 4.03 / 3.64.
  - M_cf 2.915 / 2.921 (W1 2.938 / 2.936), IN [2.5, 6]; L_z(16) 0.009692; (|k|c_eq)² 7.37 %. Q/estimate 0.855–1.552 on W2, all IN [0.6, 2.5].
  - Ball masses W2: 0.99967 / 1.00007 g without, 0.98208 / 0.98270 g with (predicted 1.000 / 0.9817).
  - C2: ΔP_phantom/P_phantom +4.80 / −0.41 / +2.07 / −0.53 % (P1–P4, W2) — P2 and P4 OUT of [1 %, 15 %]; ΔΣ|S_i1|² −3.55e-5, max |ΔS| 1.32e-4 (IN). σ_x pairs sit at the mirror pairs' floor (≤ 0.096 %), as predicted.
  - C3 on W2: kept-mirror pairs (P2/P4, ccw/cw, top/bottom) ≤ 0.69 % without and ≤ 0.31 % with the rod (IN ≤ 1 %). On W1's record mesh they reach 7.7 % (P2/P4) and 35 % (top/bottom, P4) without the rod — the 537-cell phantom cannot hold a 1 g ball, as finding 1 said. **Print-definition defect:** the module's P1/P3 spread divides by the pair mean (ceiling 200 %), not the maximum the prediction assumed (ceiling 100 %); it printed 153–160 %. Re-normalised by the maximum from the printed SARs, W2 reads 86.8 % without and 88.8 % with the rod, inside [50 %, 95 %].
- **Next:** chain I6 is done; the weekly that reads this commissions `ANS-10`.

**Audit, 2026-10-02 18:00 daily review — PASS** (`auditor`; the review re-read W2:3822, (M) 0.0819 %). Evidence: `20261002T170428Z_IMP-3-step1-W1.log:3686–3720,3765–3770,3840,3912–3914`; `20261002T170843Z_IMP-3-step1-W2.log:3756,3798–3839,3860–3877,3947,4160–4162`; asserts `tests/validation/test_imp3_rod_birdcage.py:619–634,700–705,720–722`. Both log headers carry the parent `8359727c`; module mtime 17:04:04Z precedes W1 (17:04:28Z). Caveats recorded: reciprocity's worst case is 1.507e-15 (the row's "≤ 1.5e-15" is rounding); (M) does not discriminate the rod (the without-rod arm shows the same spread scale; the phantom absorbs ≈ 3e-4 of the source power); the P2/P4 Q, η_z and the mean-normalised C3 print are printed outside their predictions with no status effect.

### §9 item 235 — 2026-10-04 (chain I9)

*Drafted by a general-purpose agent for the 2026-10-04 18:00 daily review; the review reviews it.* Licence: §10 chain I9 (PROJECT_PLAN.md ≈ line 4598), a weekly-enumerated step (§9 rule (4)); `ANS-10` was re-homed to this Palace twin by the 10-03 weekly (≈ line 4483). No band, tolerance, Done-when or validation target moves.

**`IMP-3` step 3 — the Palace twin.** Palace's 4×4 at 64 MHz, order 1, on the rod fixture with and without the rod, its own identities asserted; ours at degree 1 on the same exported meshes; ΔS(rod) per class, ours against Palace, printed.
- **Executor:** implementer. New probe + new test module; no `src/`, no edit to any existing file.
- **Tier:** heavy. Export `-n 1`; Palace `-np 8`; ours `-n 2`. Complex container for every non-Palace window.
- **Class:** chain step (solve), measured median ≈ 18–23 min.
- **Cost:** ≈ 20 measured / ≈ 45 predicted min (15 fixed + ≈ 15 code + ≈ 7 compute + ≈ 8 spare).
- **Depends on:** nothing. `GEO-36`'s hole build and `leg_sheet_strip_fraction` are both on `main`.

**Settled from the code (read before designing).**
- **The rod is a PEC hole, not a meshed region.** `_build_rod()` (`tests/mesh/test_birdcage_implant_rod.py:125`) cuts the rod out of the phantom; its wall is exterior facet tag `BIRDCAGE_IMPLANT_SURFACE_TAG` = 402, every 402 facet exterior with its cell in tag 3 (`GEO-36` (A5)). Ours pins it by the default exterior PEC (`pec_facet_tags=None`, `IMP-3` (E), 7 437 dofs, max |E| 0). Palace therefore needs the rod wall as its **own boundary attribute listed under `PEC`**. The kept-cells build (`keep_cells=True`, item 219) is `WF-9`'s, not this item's.
- **`IMP-3`'s fixture uses ragged port sheets.** `test_birdcage_port_sheets._build(True)` / `_build_rod()` pass no `leg_sheet_strip_fraction`, so `build_four_port_sweep` narrows each sheet by the midpoint filter: 26 triangles per port (`20261002T170428Z_IMP-3-step1-W1.log:1837,3646`). **Palace refuses that set** (`MFEM abort: Specified direction does not align …`, `20260929T110504Z_OPS-61-step1a-W3.log:86`), and the 10-02 binding ruling makes a new fixture's ports exact CAD rectangles. So this twin runs on the **exact-rectangle sibling** of `IMP-3`'s fixture: the same two builders with `leg_sheet_strip_fraction=GATED_WIDTH_FRACTION`. Nothing in the generator refuses the rod + strip combination (`io/mesh.py:4244–4271`), but **it has never been built**. `IMP-3`'s ragged records stay valid as records.
- **Ours' ΔS(rod) per class at degree 1 is not a tracked record.** `IMP-3` W1 printed both arms' S to 10 digits (`…170428Z_IMP-3-step1-W1.log:3682–3685, 3704–3707`) and max |ΔS| = 2.663e-4 (`:3727`; W2 1.32e-4). No per-class ΔS and no module constant exist, and those are ragged meshes anyway. Computed from those lines for scale (this brief's arithmetic, not a record): |Δ mean_k| = 1.25e-4 / 3.48e-5 / 3.80e-5, i.e. 2.1e-4 / 8.3e-5 / 1.1e-4 of |mean_k| (self / adjacent / opposite). The same arms' class spreads are 0.057 / 0.060 / 0.037 % without the rod (`:3687`). **ΔS(rod) between two separately meshed arms is therefore at or below the remesh floor.** Ours against Palace on the same mesh pair is still a sharp test, because both codes see the same two meshes. But the remesh-free rod effect is ΔS_wall = S(rod mesh, PEC wall) − S(rod mesh, natural wall), on one mesh. That is printed beside ΔS(rod).

**Code (new files only; additive; disclosed).**
- **New probe `scripts/probes/imp3_step3_palace_rod.py`.** It loads `scripts/probes/ops61_step1a_palace_ans4.py` by `importlib` (precedent `tests/validation/test_wf5_q_two_ways.py:71–78`) for `_tri_areas`, `_read_port_s`, `_identities`, `_p_check`, `_fmt`, `NARROW_BASE`, `EXTERIOR_TAG`. **Do not edit the OPS-61 probe.** Its `export` cannot be re-run as a regression without overwriting step 1b's `f05x.json`. Output directory `logs/imp3/` (gitignored).
  - **`export`** (serial; `assert comm.size == 1`). Build without the rod first: `_build(True, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`. Then the rod: `_build_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`. For each arm:
    - Take the drive set from **`build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=<tuple>, build_only=True)`** (`tests/validation/test_port_birdcage_four_port.py:219,516–529`): N_i = `case["facet_tags"].find(SHEET_IFACE + i)`. These are exactly the facets ours drives, so no second narrowing code exists.
    - Split `exterior_facet_indices` into the rod wall (the **build tuple's** `facet_tags.find(402)`) and the rest. Write triangles with attribute 1 (outer), 20 + i (N_i) and **402** (rod wall).
    - Write the tets with their cell tags, positively oriented. Copy the Gmsh 2.2 writer of `export` (`:248–283`); do not lift it.
    - Meshes are `imp3x_off.msh` and `imp3x_rod.msh`.
    - Configs use the materials, `LumpedPort` (R 50, `+Z`, Attributes `[20+i]`), `Order 1`, `Freq [0.064]` and `Linear` block of `export`'s `config()` (`:300–342`):
      - `imp3x_off.json`: `PEC [1]`;
      - `imp3x_rod.json`: `PEC [1, 402]`;
      - `imp3x_rodnat.json`: `PEC [1]`, `PMC [402]`. Check the `PMC` key spelling in the image's `config-schema.json`, as step 1a's W0 did (`:296–297`). If it is absent, leave 402 unassigned (Palace's natural default) and journal which form ran.
    - `Output "postpro_<name>"`. Assert no target file exists before writing (never overwrite).
    - Write `imp3x_meta.json`: cells, per-port N_i facets / area / h, and the wall facet count and area, per arm.
    - Anchors, asserted:
      - (X1) without-rod cells within `CELL_COUNT_BAND` of `ANS4_EXACT_STRIP_CELL_RECORD` 117 206 (`test_ops61_exact_strip_ports.py:64`); exact equality printed;
      - (X2) rod-arm cells < `ROD_CELL_GUARD_W2` 320 000 (`test_imp3_rod_birdcage.py:115`), checked before writing;
      - (X3) every N_i area = 5.6e-5 within `STRIP_BAND` (`:58`), both arms;
      - (X4) outer exterior area 0.3072 within the probe's `EXTERIOR_AREA_BAND` 1e-9, both arms;
      - (X5) rod arm: every 402 facet exterior, count > 0, meshed wall area / (2πrℓ + 4πr²) inside `TERMINAL_AREA_BAND` (`tests/mesh/test_birdcage_leg_gaps.py:73`); without-rod arm: 402 count 0.
  - **`read`.**
    - Run `_p_check(name, path, "[IMP-3 3 W5]")` on `imp3x_off` and `imp3x_rod`; their fails are asserted. On `imp3x_rodnat` it is printed only, because that run is the control.
    - Print per class: Palace ΔS(rod) = mean_k(S_rod) − mean_k(S_off), and Palace ΔS_wall = mean_k(S_rod) − mean_k(S_rodnat), as |Δ|/|mean_k(S_off)|.
    - Print Palace's σ_y pairs on S (`MIRROR_PAIRS`, `test_imp3_rod_birdcage.py:148`) beside `ADJACENT_SPREAD_BAND`.
    - Print Palace `imp3x_off` against `logs/ops61/postpro_f05x/port-S.csv`, read-only, per class.
    - Never read `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
  - Dispatch is `export` | `read`, nothing else.
- **New module `tests/validation/test_imp3_palace_twin.py`** (`@complex_only`, one test):
  - Same two builds, same order, through `build_four_port_sweep(…, built=…, build_only=True)`, then `run_n_port_sparameter_sweep(problem, port_defs, lumped_sheet_ports=specs, lumped_sheet_facet_tags=facet_tags, degree=1)` (pattern `test_imp3_rod_birdcage.py:299–321`). `keep_fields=True` on the two rod-arm sweeps only. Free the without-rod arm before the rod build.
  - **The natural-wall sweep** is a fresh sweep on `replace(problem, facet_tags=outer, pec_facet_tags=(OUTER_BOX_TAG,))`, with `outer` from `IMP-3`'s `_outer_box_tags(msh, build_facet_tags, comm)` (`:205`). This is `IMP-3` C1's construction (`:655`), applied to all four drives.
  - It reads Palace's three CSVs through the new probe (`importlib`), prints the comparisons below, and saves ours' three S matrices to `logs/imp3/ours_imp3x_p1.json` on rank 0.
  - κ stays off (the default `width_correction_kappa=None`; binding ruling 10-02).

**Anchors (asserted; bands imported, never restated).**
- **(P) Palace, W5,** on `imp3x_off` and `imp3x_rod`. This is item 159's set as `_p_check` applies it:
  - ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`);
  - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every column power ≤ 1 + 1e-9;
  - every S-class spread ≤ `ADJACENT_SPREAD_BAND` (`test_port_birdcage_lumped_column.py:121`);
  - pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst (`test_port_birdcage_four_port.py:129,173`).
  - On the rod arm the C4 spreads are backed by `IMP-3`'s rod arm (≤ 0.083 % on the ragged sibling, `:3709`).
- **Palace log check (executor, journalled), W2–W4:**
  - `elements` total = the export's count;
  - `Rs = 4.375e+01 Ω/sq` on attributes 21–24 (= 50 · 0.007 / 0.008);
  - ω/2π = 6.400e-02 GHz;
  - four `GMRES solver converged` lines.
- **Export (X1)–(X5), W1** (above).
- **Ours, W6:**
  - (F) every N_i = 5.6e-5 within `STRIP_BAND`, with w = A/h = 7.0e-3, on both arms. Cells equal `imp3x_meta.json`'s **exactly**, per arm: same builder, deterministic gmsh, serial mesh.
  - (i) and (ii) as `IMP-3` asserts them (`:621–623`) on all three sweeps.
  - (E) rod PEC arm: max |E| on the reduced tag-402 dofs ≤ `PEC_WALL_MAX_E` (`:135`) for every drive.
  - (Kw) natural arm, P1: max |E| on those dofs > `PEC_WALL_MAX_E`. This is backed on the ragged sibling by `IMP-3` C1, which measured 1.29e-4 against 1e-12.
  - (M) σ_y pairs on the rod PEC Z ≤ `ADJACENT_SPREAD_BAND`.
  - (S) without-rod separation ≥ `POOLED_SEPARATION_FLOOR`.

**Printed, never gated** (ours against Palace, per class, |Δ mean|/|mean|, via `_circulant_classes`):
- **(T) The matched-order tie, per arm** (off, rod PEC, rod natural): ours S against Palace S. Predicted ≤ 1e-6. `f05x` measured 2e-10 / 1e-10 / 6e-11 (item 159).
- **ΔS(rod), ours against Palace:** |Δ_k,ours − Δ_k,Palace| / |Δ_k,Palace|, plus the entrywise max |ΔS_ours − ΔS_Palace| / max |ΔS_Palace|. Predicted ≤ 1e-4. The derivation: a matched-order S gap of ≈ 1e-10 · |S| ≈ 6e-11 absolute, against |Δ_k| ≥ 1e-5, gives ≈ 1e-5 relative; the margin ×10 is for the smallest class.
- **ΔS_wall, ours against Palace:** the same form. Predicted ≤ 1e-3, because |Δ_wall| may be as small as ≈ 1e-6 · |S|.
- **Sizes:**
  - |Δ_k(rod)|/|mean_k| predicted in [1e-5, 1e-3]; the ragged pair read 2.1e-4 / 8.3e-5 / 1.1e-4.
  - |Δ_k,wall|/|mean_k| predicted in [1e-6, 1e-3].
- **Free ties:**
  - Ours without-rod against `ANS4_EXACT_STRIP_S_P1_RECORDS["64"]` (`test_ans4_exact_strip_records.py:141`). Predicted ≤ 1e-6 if `test_birdcage_port_sheets._build` and `test_birdcage_leg_offset._build` coincide; both read 116 085 with the flag off.
  - Palace `imp3x_off` against Palace `f05x`, predicted ≤ 1e-9 on the same condition.
  - A miss on either is a reading about the builders, not a stop.
- **Palace prices:** total time, iterations, peak memory and ND unknowns per run.

**Negative controls (rule (e); unbacked on this fixture, so predicted and printed, with separations computed).**
- **(C1) The wrong wall:** ours rod (PEC) against Palace `imp3x_rodnat`, per class.
  - This gap equals |Δ_wall|/|mean|, predicted in [1e-6, 1e-3].
  - Against the matched (T) gap (≤ 1e-9, measured ≈ 1e-10), the separation is predicted ≥ 1e3×. The least case is 1e-6 / 1e-9.
  - Ceiling ≈ 1e-3 / 6e-11 ≈ 1.7e7×. The round-off floor of (T) bounds it.
  - It shows that the 402 attribute reached Palace's solve. A PEC list that silently dropped 402 would make `imp3x_rod` equal `imp3x_rodnat`, so (T) on the rod arm would read ≈ |Δ_wall|, not ≈ 1e-10.
- **(C2) The arm swap:** ours without-rod against Palace `imp3x_rod`. This gap equals |Δ(rod)|/|mean| plus the remesh change, predicted [1e-5, 1e-3]. Separation against (T) is ≥ 1e4×.
- **(Kw)** is the asserted knob (above). It is backed on the same quantity on the ragged sibling with ≈ 1.3e8× separation.
- **(P)'s power to discriminate** was shown on `f05x` by item 159's flip (200.02 % spreads against the 0.5 % band, ceiling 400×). It is not re-run.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`).
- **Pre-check** (not a harness call): `ls logs/imp3/`. If any `imp3x_*` file exists, stop and report: never overwrite. `logs/ops61/` is read-only to this item.
- **W1, export** (≈ 90 s, range 60–150). Basis: strip export 31 s at `-n 1` (F10a W3); rod mesh 42.1–47.4 s (`GEO-36` W2; `IMP-3` W1 `:3646`).

  `scripts/testing/run_and_log.sh IMP-3-item235-W1-export "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/imp3_step3_palace_rod.py export'"`
- **W2, Palace without the rod** (≈ 25 s; `f05x` o1: 23.35 s Palace total, 1.8 G, item 159). This is **the box-speed gate**.

  `scripts/testing/run_and_log.sh IMP-3-item235-W2-palace-off "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/imp3 && timeout -k 30 590 palace -np 8 imp3x_off.json'"`
- **Box-speed gate (binding, after W2).** r = W2's Palace `Total` / 23.35 s (`20261002T021838Z_OPS-61-F10b-W2-palace-f05x.log`).
  - **r ≤ 2.5:** run everything.
  - **r > 2.5:** run W3–W5 (the Palace record still fits 590 s), skip W6, and mark BLOCKED, slow box. W6 projects to ≥ 375 s against its 420 s ceiling.
- **W3, Palace with the rod** (predicted 35–60 s Palace total, peak 2.3–3.5 G).
  - Cells ≈ 168 500 (167 371 + the strip's +1 121), ×1.44 `f05x`; ND (p = 1) predicted 190–215 k.
  - The command is W2's with `imp3x_rod.json`, labelled `IMP-3-item235-W3-palace-rod`.
- **W4, Palace with the natural wall** (as W3): `imp3x_rodnat.json`, labelled `IMP-3-item235-W4-palace-rodnat`. After any kill, check `docker compose -f docker/docker-compose.yml ps` for a lingering `fem-em-palace` container.
- **W5, read** (≈ 3 s). This is the probe's last window; no probe edit after it starts.

  `scripts/testing/run_and_log.sh IMP-3-item235-W5-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/imp3_step3_palace_rod.py read'"`
- **W6, ours** (≈ 150 s, range 110–250). Basis, `IMP-3` W1 at `-n 2`: builds 25.9 + 47.4 s, four solves ≈ 7 s per arm, module 241 s with SAR and C1. The closing window; no module edit after it starts.

  `scripts/testing/run_and_log.sh IMP-3-item235-W6-ours-p1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp3_palace_twin.py -v -s --tb=short'"`

**Traps already paid for.**
- **Never hand Palace a ragged sheet set** (the step 1a abort) and never build the rod fixture without `leg_sheet_strip_fraction`. An abort on any port set is a finding: never change `Direction` in-slot.
- **The exterior includes the rod wall.** Exporting `exterior_facet_indices` as attribute 1 would make the wall PEC in all three configs and break (X4) and (C1). Split it first.
- **Tag 402 lives in the build tuple's facet tags**, not in the sweep's (`IMP-3` trap).
- **Build order and width.** Without the rod first, then the rod (`GEO-23`). The new keyword combination builds at `-n 1` first (W1), because a gmsh failure deadlocks wider runs (F10a). A gmsh FAIL is recorded verbatim and never retried.
- **Import only helpers and constants from `test_imp3_rod_birdcage`.** Never import a `test_*` name: pytest would collect it and run `IMP-3`'s W1 inside W6.
- **A held factorisation silently solves a different problem** (`PORT-19`, `WF-12`). The natural-wall arm is its own sweep on the replaced problem, and (Kw) proves the knob.
- Areas, counts, `assemble_scalar` and dof counts are rank-local: reduce. The export is serial.
- Compare S classes only. Never read Palace's `z_matrix_ohm`.
- A `GMRES … did not converge` line is a finding: that run's S is not measured. Never change `Order`, `MaxIts` or `Tol` in-slot.
- Harness discipline: `-s`, `tests/environment` first, `FEM_EM_REQUIRE_COMPLEX=1`, foreground, no pipes. On Status 124, run `pgrep -c python3` first and never lengthen the timeout.
- Palace's numbers go in the record at full detail (version, config, mesh cited). No AED number, and no Palace figure set against an AED one, in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Shared files.** Item 234 (`OPS-61` step 2, the `ANS-6` hole mesh in Palace) is expected to add modes to `scripts/probes/ops61_step1a_palace_ans4.py`. This item only **imports** from that file, so the two do not overlap. If item 234 lands first, check that `_p_check(name, path, tag)` and `_read_port_s(path)` keep their signatures. If this item ever needs to touch that file, the change is additive in separate functions with separate dispatch branches: rebase and keep both sides.

**Scope.** Covered: Palace order 1 and ours at degree 1, 64 MHz, on the exact-rectangle sibling of `IMP-3`'s rod fixture (σ = 800 S/m untuned coil, 50 Ω, one axial rod on leg 1's azimuth): S per arm, ΔS(rod) and ΔS_wall per class. Agreement at matched order is a check of the PEC-hole boundary and port bookkeeping between two independent codes. It is **not** a check of discretisation error, and **not** of the rod's physical effect on S, which sits at the remesh floor. Not covered:
- **The tip SAR, the tip ratio Q and its accuracy.** That was `ANS-10`'s original question (§2 ≈ line 437). Palace's fields are not read here, so the question stays open.
- Tissue, heating, compliance, ISO 10974 (`WF-9`), ASTM F2182 (7c).
- Absolute S / `Z_in`, any h- or p-statement (one mesh, one order; Palace order 2 would cost ≈ 270 s × 1.44 ≈ 6–7 min per run, which the review can price).
- The tuned or copper coil (`implant_rod` refuses `as_hole`), F-human, other rod placements, the kept-cells rod.
- AED, `src/` and any band.

**Status it can move.**
- §7 `IMP-3` **stays ✅** (no glyph). Its Done-when names step 1 only, so under rule (j) the review adds a step-3 record sentence. The candidate text, which holds only if W5 and W6 are green and every printed (T) gap is ≤ 1e-6: "Palace order 1 matches ours at degree 1 on the exact-rectangle rod fixture to ≤ 1e-6 per class on both arms; ΔS(rod) and ΔS_wall agree at matched order".
- §10 I9 is done.
- `ANS-10` has no §7 row. The weekly reads this as its Palace reading: agree / disagree / inconclusive on the rod's 4×4 effect. It does not answer the tip ratio's accuracy.
- §2 ≈ line 437 may gain one clause ("the rod's 4×4 effect agrees with Palace at matched order"). Its "ratio's accuracy" clause stays.

**Negative result.**
- **W1 gmsh FAIL on the rod + strip build:** record it verbatim, open a known-issues entry, BLOCKED on the mesher. Never retry.
- **(X1)–(X5) red:** the export is not the fixture. Stop before any Palace window and mark BLOCKED.
- **A Palace abort, or a log lacking the export's element count or `Rs = 4.375e+01`:** record it verbatim, open a known-issues entry, BLOCKED (rule (d)).
- **The `PMC` key refused:** re-write `imp3x_rodnat.json` with 402 unassigned. If Palace then rejects an unassigned exterior attribute, (C1) is "not measured". Not a stop.
- **(P) red on `imp3x_off` or `imp3x_rod`:** report every identity's number, park, BLOCKED. Never widen.
- **Ours (F), (i), (ii), (E), (Kw), (M) or (S) red:** park on `attempt/IMP-3-step3-<ts>`, open a known-issues entry, BLOCKED.
- **(T) without the rod > 1e-6:** ours did not solve the exported mesh. Stop, BLOCKED; do not print ΔS.
- **(T) without the rod green, but the rod arm > 1e-6:** a reading, not a failure. It would be the first disagreement between the two codes, located at the PEC hole. Print it per class and report; the review reads it.
- **Box-speed gate trips:** a partial record (Palace only), BLOCKED (slow box). The review re-queues W6.
- **A printed size outside its prediction:** report it; no status effect.
- No `PREDICTED RED` line.

**Outcome, 2026-10-04 22:30 slot — ✅ DONE, all six windows green.** New files `scripts/probes/imp3_step3_palace_rod.py`, `tests/validation/test_imp3_palace_twin.py`; the OPS-61 probe was only imported.
- W1 `20261005T034848Z_IMP-3-item235-W1-export.log` (79 s, `-n 1`): (X1) 117 206 = record exactly; (X2) rod arm 168 507 < 320 000; (X3) every N_i 5.600000000000e-05 (|Δ| ≤ 5.4e-20), h 8.0e-3, w 7.0e-3, 40 facets per port on both arms; (X4) outer 0.3072 to 5.6e-17 (3 318 / 3 316 tris); (X5) 4 958 wall facets, all exterior, area ratio 0.996254. The rod + strip combination built first time.
- `imp3x_rodnat` ran with 402 left **unassigned** (Palace's natural wall; item 234's `h05x_pmc` form), not a `PMC` key.
- Palace v0.18.1 `-np 8`, all checks met (elements = export count, `Rs = 4.375e+01 Ω/sq` on 21–24, 6.400e-02 GHz, four `GMRES solver converged`): W2 off `20261005T035015Z` Total 24.364 s, 1.8 G, ND 140 460, 31–32 its — box-speed r = 1.043; W3 rod `20261005T035048Z` 34.19 s, 2.5 G, ND 204 156; W4 rodnat `20261005T035124Z` 35.46 s, 2.5 G.
- W5 `20261005T035220Z_IMP-3-item235-W5-read.log` (2 s): (P) PASS on off (reciprocity 1.556e-10, σ_max 0.999760093669, spreads 0.0852 / 0.0783 / 0.0921 %, separation 439.6×) and rod (1.390e-10, 0.999767876823, 0.0937 / 0.0885 / 0.0930 %, 432.2×); rodnat printed (all PASS). Palace imp3x_off vs f05x ≤ 2.9e-14 (free tie). σ_y pairs on S ≤ 7.4e-4.
- W6 `20261005T035322Z_IMP-3-item235-W6-ours-p1.log` (172 s, `-n 2`, 13 passed): (F) cells = meta exactly, areas 5.6e-5; (i)/(ii) reciprocity ≤ 1.8e-15, σ_max ≤ 0.99976788; (E) 0.0 on 7 437 wall dofs, every drive; (Kw) natural P1 1.297e-4 > 1e-12; (M) σ_y Z ≤ 6.5e-4; (S) 451.7×.
- **(T)** ours vs Palace per class (self / adjacent / opposite): off 2.168e-10 / 1.412e-10 / 6.412e-11; rod 1.952e-10 / 9.638e-11 / 1.392e-10; rodnat 2.103e-10 / 1.446e-10 / 4.585e-11.
- |Δ_k(rod)|/|mean_k| 1.477e-4 / 7.477e-5 / 6.268e-5 (both codes; in [1e-5, 1e-3]); |Δ_k,wall|/|mean_k| 2.223e-5 / 5.317e-6 / 3.380e-5 (in [1e-6, 1e-3]).
- ΔS(rod) ours vs Palace 2.481e-7 / 5.994e-7 / 1.304e-6, entrywise 6.080e-7 (≤ 1e-4); ΔS_wall 2.239e-6 / 9.342e-6 / 2.872e-6, entrywise 3.331e-6 (≤ 1e-3).
- (C1) wrong wall 2.223e-5 / 5.317e-6 / 3.380e-5 — separation over the rod (T) 1.1e5 / 5.5e4 / 2.4e5×; (C2) arm swap 1.477e-4 / 7.477e-5 / 6.269e-5 — 6.8e5 / 5.3e5 / 9.8e5× over the off (T). Free tie ours off vs `ANS4_EXACT_STRIP_S_P1_RECORDS["64"]` ≤ 1.2e-10.
- Reading: the PEC-hole boundary and port bookkeeping agree between the two codes at round-off; the rod's 4×4 effect agrees at matched order. Not a tip-SAR, h- or p-statement.

### §9 item 242 — 2026-10-05 (chain I12, daily — provisional)

*Drafted by a general-purpose agent for the 2026-10-05 03:00 daily review; the review reviews it.* Licence: §10 chain I12 **(daily, 2026-10-05 — provisional)** under rule (6) (PROJECT_PLAN.md ≈ lines 3694–3704), the next step after I9 toward chain I's stated target (Palace). `ANS-10` was re-homed to Palace by the 10-03 weekly (≈ line 4493). Item 235's journal named this step (`docs/testing/attempts.md:2346`): the tip question "needs Palace field output, not more S". No band, tolerance, Done-when or validation target moves.

**`IMP-3` step 4 — the near-tip field, ours against Palace.** This runs on item 235's two exact-rectangle meshes, at order 1 / degree 1, at 64 MHz, with drive P1 as the reading of record. It compares E at six fixed phantom points and the phantom's electric energy, with the rod (PEC wall), with the natural wall, and without the rod. Palace's identities are asserted. The rod's tip-field effect is printed beside the two codes' disagreement. **Target: Palace — a reading, not a gate**, as in I9.
- **Executor:** implementer. One new probe and one new test module. No `src/`, and no edit to any existing file.
- **Tier:** heavy. Configs `-n 1`; Palace `-np 8`; ours `-n 2`. Use the complex container for every non-Palace window.
- **Class:** chain step (solve). Measured median ≈ 18–23 min (item 235's brief).
- **Cost:** ≈ 18 measured / ≈ 40 predicted min (15 fixed + ≈ 12 code + ≈ 7 compute + ≈ 6 spare). Basis: item 235 took ≈ 10 min of compute in a ≈ 30 min slot (`attempts.md:2341–2345`; §9 line "≈ 16 measured"). This item rebuilds no mesh.
- **Depends on:** item 235's gitignored artifacts in `logs/imp3/`, which are present on this box (pre-check below).

**Feasibility, settled from the repo (read before designing).**
- **Palace v0.18.1 writes both quantities as CSV, so no field-format reader is needed.** The image's own driven example `/opt/palace/examples/cpw/cpw_lumped_uniform.json` carries `Domains.Postprocessing.Energy: [{Index, Attributes}]` and `Domains.Postprocessing.Probe: [{Index, Center}]` (`docs/testing/logs/20260929T110106Z_OPS-61-step1a-W0.log:65,97–117`). Palace then writes `probe-E.csv` (complex E at each `Center`, per excitation) and per-domain columns in `domain-E.csv`.
- Item 235's runs already write `domain-E.csv`, with total energies only and one column set per excitation: `E_elec[1] (J)` … `E_ind[4] (J)` (`logs/imp3/postpro_imp3x_rod/domain-E.csv:1`). `OutputFormats` defaults to `Paraview: true`, but nothing is saved without `Save`. ParaView / GridFunction output is **not** used.
- **Ours evaluates on a full mesh, not a submesh.** `_build` / `_build_rod` cut the rod as a gmsh hole; no `create_submesh` appears in `io/mesh.py`, `test_birdcage_implant_rod.py`, `test_birdcage_port_sheets.py`, `test_port_birdcage_four_port.py` or `test_imp3_palace_twin.py`. The owner rule `evaluate_vector_field_parallel(f, pts, owner="lowest_original_cell")` (`src/fem_em_solver/post/evaluation.py:131,279–286`) is therefore usable. Known-issues entry 2026-10-04 (owner rule vs `create_submesh`) does not bite.
- **Normalisation.** Palace drives each port at V_inc = 7.071 V on 50 Ω (`logs/imp3/postpro_imp3x_rod/port-V.csv:2`). Ours has its own incident scale, the same in every sweep of one code. So every comparison below is either a **ratio inside one code** (arm to arm, same incident) or **one fitted complex constant g** shared by all arms and points. Absolute V/m are never compared.
- **The time convention is shared.** Item 235's `_read_port_s` builds S as |S|·e^{i·arg} with no conjugation (`scripts/probes/ops61_step1a_palace_ans4.py:391–414`) and matched ours to 2e-10, so the phase convention agrees. The conjugate form is printed as a trap check.

**Probe points** (metres, L0 = 1; all in phantom tag 3, outside the rod hole; generic offsets keep them off symmetry planes and mesh faces). The rod has r = 1 mm, centre (0.015, 0, 0), and tips at z = ±0.016 (`tests/mesh/test_birdcage_implant_rod.py:93–96`). The phantom is R 0.03 × H 0.08 (`test_birdcage_port_tags.py:28–29`).
- T1 (0.015023, 0.000017, 0.017531): 1.53 mm beyond the top tip.
- T2 (0.015023, 0.000017, 0.019031): 3.03 mm beyond.
- T3 (0.015023, 0.000017, 0.022031): 6.03 mm beyond.
- S1 (0.017523, 0.000017, 0.000031): beside the shaft, 1.52 mm off the wall.
- F1 (−0.014977, 0.000017, 0.000031): 30 mm from the rod axis. This is the reference point of g.
- F2 (0.000023, −0.017983, 0.010031): ≈ 23 mm from the axis.

**Code (new files only; additive; disclosed).**
- **New probe `scripts/probes/imp3_step4_palace_fields.py`.** It loads `scripts/probes/imp3_step3_palace_rod.py` by `importlib` for `OUT_DIR`, `RUNS`, `csv_of`, `class_means`, `line` and `_ops61()` (and through it `_p_check`, `_read_port_s`). **Do not edit either earlier probe.** Dispatch is `configs` | `read`, nothing else.
  - **`configs`** (serial; `assert comm.size == 1`):
    - **(X0)** Parse the `$Elements` blocks of `logs/imp3/imp3x_off.msh` and `imp3x_rod.msh` read-only. The type-4 (tet) counts must equal `imp3x_meta.json`'s `cells` **exactly** (117 206 / 168 507), and the type-2 triangle attribute set must be {1, 21–24}, plus {402} on the rod arm.
    - Load `imp3x_off.json`, `imp3x_rod.json` and `imp3x_rodnat.json` read-only. Assert their PEC lists are [1] / [1, 402] / [1], as item 235 wrote them (`imp3_step3_palace_rod.py:257–262`).
    - Deep-copy each config, add `Domains.Postprocessing = {"Energy": [{"Index": 1, "Attributes": [PHANTOM_CELL_TAG]}], "Probe": [{"Index": k, "Center": p_k}, k = 1…6]}` and set `Problem.Output = "postpro_imp3p_<arm>"`. Write `imp3p_off.json`, `imp3p_rod.json` and `imp3p_rodnat.json`.
    - Also write the energy-only fallbacks `imp3e_<arm>.json` (no `Probe`, Output `postpro_imp3e_<arm>`). They run **only** if Palace refuses `Probe`.
    - Open every target with mode `"x"`: never overwrite. `PHANTOM_CELL_TAG` = 3 is imported (`test_port_birdcage_lumped_column.py`, as item 235 imports it).
  - **`read`:**
    - **(H)** Parse `probe-E.csv` and `domain-E.csv` in each `postpro_imp3p_*` **by header**, mapping each column to (probe, component, Re/Im, excitation) and (domain, excitation). The map must be a bijection onto 6 × 3 complex × 4 excitations and onto one Energy domain × 4 excitations. If it is ambiguous, print the header verbatim and stop.
    - Run `_p_check` on `imp3p_off` / `imp3p_rod` (fails asserted). On `imp3p_rodnat`, print only.
    - **(RP)** Each `imp3p_*` S must match item 235's `postpro_imp3x_*` S entrywise.
    - Print Palace's own η (below). Save the parsed probe and energy arrays to `logs/imp3/palace_imp3p_fields.json` (mode `"x"`).
    - Never read `z_matrix_ohm`, `docs/private/` or `aed_results/`.
- **New module `tests/validation/test_imp3_palace_fields.py`** (`@complex_only`, one test):
  - Use the same builds, the same order and the same three sweeps as `test_imp3_palace_twin.py:110–159`, through `probe.build_arm(arm)` → `build_four_port_sweep(…, built=…, build_only=True)` → `run_n_port_sparameter_sweep(…, degree=1)`. **`keep_fields=True` on all three sweeps.** The natural-wall sweep is the same `replace(problem, facet_tags=outer, pec_facet_tags=(OUTER_BOX_TAG,))` on `_outer_box_tags(msh, built[2], comm)`. κ stays off.
  - Import helpers only, never a `test_*` name: `_sweep`, `_ids` and `_sheet_areas` from `test_imp3_palace_twin`, and `_matrix_quantities`, `_wall_dofs`, `_wall_max_e` and `PEC_WALL_MAX_E` from `test_imp3_rod_birdcage`.
  - Per arm, **before** `del res`, do three things for every drive:
    1. E at the six points: `evaluate_vector_field_parallel(res.fields[p].e_complex, PTS, owner="lowest_original_cell")`, with `PTS` a module constant identical on every rank.
    2. W_ph ∝ ∫_tag3 e·ē dx: `assemble_scalar` over `ufl.Measure("dx", domain=msh, subdomain_data=built[1], metadata={"quadrature_degree": 2})`, then `allreduce(SUM)` and take `.real`.
    3. **(B)** For each point, the number of **owned** colliding cells (`geometry.bb_tree` / `compute_collisions_points` / `compute_colliding_cells`, as in `evaluation.py:270–278`, links < `size_local`), reduced with SUM.
  - Write `logs/imp3/ours_imp3p_fields.json` on rank 0 with mode `"w"`, a new name. **Never write `ours_imp3x_p1.json`.**
  - Read `palace_imp3p_fields.json` and print the comparisons below.

**Anchors (asserted; bands imported, never restated).**
- **(X0)**, configs window: exact tet counts and attribute sets (above). This is item 235's (X1)/(X2) count, re-read from the file that Palace reads.
- **(P) Palace**, W5, on `imp3p_off` / `imp3p_rod`, via `_p_check` (`ops61_step1a_palace_ans4.py:742–782`):
  - `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`);
  - `PASSIVITY_SIGMA_TOLERANCE`, with column power ≤ 1 + 1e-9 (`test_port_birdcage_four_port.py:129`);
  - `ADJACENT_SPREAD_BAND` (`test_port_birdcage_lumped_column.py:121`);
  - `POOLED_SEPARATION_FLOOR` (`test_port_birdcage_four_port.py:173`).
  - Backed on these meshes by item 235's W5 (reciprocity 1.556e-10 / 1.390e-10, spreads ≤ 0.0937 %).
- **(RP) Palace reproduction:** every `imp3p_*` S against item 235's `imp3x_*` S, entrywise, within `REUSE_REPRODUCTION_RTOL` 1e-6 (`tests/validation/test_port_birdcage_ring_matrix.py:171`, = `STEP2_COLUMN_POWER_RTOL` `:164`). Postprocessing does not enter the solve. The same Palace on one mesh tied at 2.9e-14 (item 235 W5, `imp3x_off` vs `f05x`).
- **(H)** The header bijection (above). It is a structure check, with no band.
- **Palace log check** (executor, journalled), W2–W4:
  - `elements` total equals the (X0) count;
  - `Rs = 4.375e+01 Ω/sq` on 21–24;
  - ω/2π = 6.400e-02 GHz;
  - four `GMRES solver converged` lines.
- **Ours, W6:**
  - (F) cells equal `imp3x_meta.json` exactly, and N_i equals 5.6e-5 within `STRIP_BAND` (`test_ops61_exact_strip_ports.py:58`).
  - **(R)** S per arm reproduces item 235's `logs/imp3/ours_imp3x_p1.json` entrywise within `REUSE_REPRODUCTION_RTOL`. It is backed by `IMP-3` (R0), 2.98e-10 against that band (`20261002T170428Z_IMP-3-step1-W1.log`). Repeated complex MUMPS solves differ at ≈ 2.7e-12 (item 237), far inside it.
  - (i)/(ii) on all three sweeps (`_ids`).
  - (E) rod PEC wall max |E| ≤ `PEC_WALL_MAX_E` (`test_imp3_rod_birdcage.py:135`), every drive.
  - (Kw) natural wall, P1: > `PEC_WALL_MAX_E` (item 235 W6: 1.297e-4).
  - **(V)** `valid_mask` is all True for the six points on all three arms. A False means a coordinate fell in the hole or outside the mesh, so the fixture is not the one designed.

**Printed, never gated** (ours against Palace, P1 the reading of record, P2–P4 printed beside):
- **(B)** The owned-cell count per point, predicted 1 everywhere. A point reading ≠ 1 sits on a face, where the two codes may evaluate different cells. It is flagged *face-ambiguous* and left out of (T_E). This is a filter, not a stop.
- **(T_E) The field tie.** g = E^ours_c*(F1) / E^Pal_c*(F1) on the without-rod arm, with c* = argmax_c |E^Pal_c(F1)|.
  - Per arm, the reading is max over (p, c) of |E^ours_c(p) − g·E^Pal_c(p)| / ‖g·E^Pal(p)‖ over the five non-reference points plus F1's other components. **Predicted ≤ 1e-6 on all three arms.**
  - Derivation: same mesh, same N1curl order 1 and the same discrete operator (S tied at ≤ 2.2e-10, item 235). Evaluation is in the same cell, so the residual is solver tolerance (Palace `Tol` 1e-10) times a local amplification. The 1e4 margin covers that amplification.
  - Print g itself (the ratio of incident scales; a constant, no prediction), and the same reading against conj(E^Pal): predicted O(1), which confirms the convention.
- **(T_W) The energy tie.** q_a = W^ours_ph,a / W^Pal_ph,a for P1. The reading is max_a |q_a / q_off − 1|, **predicted ≤ 1e-6**. Print q_off / |g|² (the codes' energy-prefactor convention; a constant, no prediction).
- **The rod's tip-field effect, inside each code, then ours against Palace:**
  - **η_wall(p)** = |E_rodPEC(p)| / |E_rodnat(p)|: one mesh, no remesh. Predicted for P1:
    - T1 [1.5, 20], T2 [1.2, 8], T3 [1.05, 4];
    - S1 [0.2, 1.5] (n × E = 0 suppresses the axial field beside the shaft);
    - F1, F2 |η − 1| ∈ [3e-3, 1e-1]. This is a PEC-spheroid dipole estimate: a·b²/(3 L_z) = 5.5e-7 m³ over r³ = 2.7e-5 / 1.7e-5 m³ ≈ 2 % / 3 %, with L_z(16) = 0.009692 (item 174 W2).
  - **η_rod(p)** = |E_rod(p)| / |E_off(p)|: two meshes, the remesh included. The tip predictions are as for η_wall. At F1 and F2, |η − 1| ∈ [1e-3, 0.3]; the without-rod phantom has 537 cells, so a point value moves with the remesh.
  - **Ours against Palace:** d_η = |η^ours − η^Pal|, predicted ≤ 2e-6·η, i.e. two ratio errors at (T_E)'s level.
- **ΔW:** W_ph(rod) / W_ph(rodnat) − 1 and W_ph(rod) / W_ph(off) − 1, per code and per drive, ours against Palace. Predicted |Δ| ∈ [1e-3, 0.1] and code gap ≤ 1e-6. `IMP-3` C2 read +4.80 % (P1) on the ragged 1 g rung; that is a record on another mesh, not a comparand.
- **Closed form beside it:** |E_cf(p)| / E₀ from `pec_prolate_spheroid_field(p − ROD_CENTRE, semi_major=0.016, semi_minor=0.001, e0=1.0)` (`src/fem_em_solver/utils/analytical.py:1576`) at T1–T3 and S1, beside η_rod(p) / η_rod(F1). The ratio is predicted in [0.3, 3]. It carries no accuracy claim: the record mesh cannot resolve a 1 mm rod's apex.
- **Palace prices:** Total time, peak memory and ND unknowns per run, as in item 235.

**Negative controls (rule (e); there is no prior field-level comparison on this fixture, so every factor is predicted and printed, with its ceiling computed).**
- **(C1) The tip effect against the two codes' disagreement (the without-rod arm sets the floor).** Separation = |η^Pal_wall(T1) − 1| / d_η(T1).
  - **Predicted ≥ 1e5.** Least case η = 1.5: 0.5 / (2e-6 · 1.5) ≈ 1.7e5.
  - **Ceiling ≈ 2.4e9:** (20 − 1) / (2 · 2e-10 · 20), with the tie at item 235's S round-off.
  - Tie_off (T_E on the rod-free arm) is printed beside it as the floor the rod has to clear.
- **(C2) Tip against far: the effect is local.** Separation = |η_wall(T1) − 1| / max_F |η_wall(F) − 1|, per code.
  - **Predicted ≥ 5** (least: 0.5 / 0.1).
  - **Ceiling ≈ 6e3** ((20 − 1) / 3e-3).
  - Below 5, the far points are not "far" on this mesh. That is a reading about placement, not a stop.
- **(C3) The arm swap: the probe reading sees which arm it is on.** Take ours without the rod at T1 against Palace with the rod at T1, with g applied: |E^ours_off(T1) − g·E^Pal_rod(T1)| / ‖g·E^Pal_rod(T1)‖.
  - **Predicted ≥ 0.3**, i.e. ≥ 3e5× over (T_E).
  - **Ceiling ≈ 1** for a relative difference against the rod field when η > 1. In ×, it is bounded by (T_E)'s round-off (≈ 1e9).
- **(Kw)** is the asserted knob on the wall (above). (RP) and (R) prove that both codes solved item 235's problems, not new ones.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`).
- **Pre-check** (not a harness call): `ls logs/imp3/ logs/imp3/postpro_imp3x_off logs/imp3/postpro_imp3x_rod logs/imp3/postpro_imp3x_rodnat`.
  - Required: `imp3x_off.msh`, `imp3x_rod.msh`, the three `imp3x_*.json`, `imp3x_meta.json`, `ours_imp3x_p1.json`, and three `port-S.csv`.
  - Absent: every `imp3p_*`, `imp3e_*`, `postpro_imp3p_*` and `*_imp3p_fields.json`.
  - Any miss ⇒ stop and report. Never re-run item 235's export.
- **W1, configs** (≈ 15 s; two ASCII `.msh` scans):

  `scripts/testing/run_and_log.sh IMP-3-item242-W1-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/imp3_step4_palace_fields.py configs'"`
- **W2, Palace without the rod** (≈ 25 s; item 235 W2 Total 24.364 s, 1.8 G, `20261005T035015Z_IMP-3-item235-W2-palace-off.log:326,342`). This is **the box-speed gate**.

  `scripts/testing/run_and_log.sh IMP-3-item242-W2-palace-off "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/imp3 && timeout -k 30 590 palace -np 8 imp3p_off.json'"`
- **Box-speed gate (binding, after W2):** r = W2 Total / 24.364 s.
  - r ≤ 2.5 ⇒ run everything.
  - r > 2.5 ⇒ run W3–W5, skip W6, BLOCKED (slow box). W6 projects to ≥ 430 s against its 420 s ceiling.
- **W3 / W4, Palace with the rod and with the natural wall.** These are W2's command with `imp3p_rod.json` / `imp3p_rodnat.json`, labelled `IMP-3-item242-W3-palace-rod` / `IMP-3-item242-W4-palace-rodnat`. Predicted ≈ 35 s each; item 235 measured 34.19 / 35.46 s and 2.5 G. Probe and Energy postprocessing add < 1 s (`Postprocessing` 0.28 s, `…035048Z_IMP-3-item235-W3-palace-rod.log:324`). After any kill, check `docker compose -f docker/docker-compose.yml ps` for a lingering `fem-em-palace` container.
- **W5, read** (≈ 3 s). This is the probe's last window: no probe edit after it starts.

  `scripts/testing/run_and_log.sh IMP-3-item242-W5-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/imp3_step4_palace_fields.py read'"`
- **W6, ours** (≈ 190 s, range 150–300). Basis: item 235 W6 took 172 s at `-n 2` for the same three sweeps (`20261005T035322Z_IMP-3-item235-W6-ours-p1.log`). Keeping the without-rod fields, 72 point evaluations and 12 energy integrals add ≈ 10–20 s. This is the closing window: no module edit after it starts.

  `scripts/testing/run_and_log.sh IMP-3-item242-W6-ours-fields "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp3_palace_fields.py -v -s --tb=short'"`
- **Fallback W2e–W4e**, only if W2 aborts on `Probe`: the same commands with `imp3e_<arm>.json`, labelled `IMP-3-item242-W2e-…` and so on. W5 and W6 then run with the probe readings marked "not measured".

**Traps already paid for.**
- **Pointwise N1curl depends on the MPI width without the owner rule** (known-issues 2026-09-27, `ans:2` Row 1: −19 % … +112 % between `-n 2` and `-n 4`).
  - Always pass `owner="lowest_original_cell"` and evaluate `e_complex` once.
  - **Do not copy `IMP-3` step 1's pattern** (`test_imp3_rod_birdcage.py:385–386`: `e_real` / `e_imag`, default owner).
  - Never `f.eval(points, np.arange(n))`. The point array must be bitwise identical on every rank (OPS-40/55 digest).
- **The owner rule raises on `create_submesh` meshes** (known-issues 2026-10-04). It does not apply here (full mesh, above). If it raises anyway, that is a finding: stop and report; never fall back to the default owner.
- **Never overwrite item 235's artifacts.** Write with mode `"x"`, use new `imp3p_` / `imp3e_` names, and leave `ours_imp3x_p1.json` alone (item 235's module writes it with `"w"`; this module must not call that code path). `logs/ops61/` is read-only.
- **Tag 402 lives in the build tuple's facet tags.** The cell tags for the energy integral are the build tuple's `built[1]`, too.
- **Build order.** Build without the rod first, then with it (`GEO-23`). Free the without-rod arm's fields after its readings and before the rod build. Evaluate before `del res`.
- **A held factorisation silently solves a different problem** (`PORT-19`, `WF-12`). The natural-wall arm is its own sweep, and (Kw) proves the knob.
- **Rank-local quantities.** `assemble_scalar`, collision counts and wall dofs are rank-local: reduce them, and count owned cells only. Keep collectives outside any `rank == 0` branch.
- **Pin `quadrature_degree`** on the energy form (POST-5). Do no UFL or `sqrt` comparisons; take ratios in numpy.
- **Palace CSV headers are parsed, never assumed.** The excitation index sits in brackets: `V_inc[1][1]` (`port-V.csv:1`), `E_elec[1]` (`domain-E.csv:1`). An unparseable header is a stop, not a guess.
- **A `GMRES … did not converge` line is a finding.** Never change `Order`, `Tol` or `MaxIts` in-slot.
- **Harness discipline:**
  - `-s`, `tests/environment` first, `FEM_EM_REQUIRE_COMPLEX=1`, foreground, no pipes, no `-k`;
  - Status 124: run `pgrep -c python3` first and never lengthen the timeout;
  - rule (i): no edit after a file's last window.
- **Palace numbers go in at full detail** (version, config, mesh). No AED number, and no Palace figure set against one, in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace order 1 and ours at degree 1, 64 MHz, on item 235's exact-rectangle record meshes (117 206 / 168 507 cells, σ = 800 S/m untuned coil, 50 Ω, κ off). The comparison is E at six fixed points and phantom electric energy, on three arms, with drive P1 of record. Agreement is a **matched-order check of the near-tip field** between two independent codes. It is **not** a discretisation statement. Not covered:
- **the 1 g tip SAR, the ratio Q and its accuracy** (`ANS-10`'s question). Palace has no mass-averaged SAR, and this mesh cannot hold a 1 g ball (finding 1, item 174). That needs the rod + strip build on `PHANTOM_RESOLUTION_1G_RUNG` (never built) and an order-2 or h-rung Palace run: a later I-step for the weekly to price;
- any h- or p-statement, or absolute |E| / SAR;
- tissue, heating, compliance, ISO 10974 (`WF-9`), ASTM F2182 (7c);
- the tuned or copper coil, F-human, other rod placements, the kept-cells rod;
- AED, `src/`, any band.

**Status it can move.**
- **§7 `IMP-3` stays ✅** (no glyph). Its Done-when names step 1 only (rule (j) read: the row's Done-when is "identities green … on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture"). The review appends a **step-4 record sentence** beside item 235's. Candidate text, which holds only if W5 and W6 are green, (T_E) ≤ 1e-6 on all three arms and (C1) ≥ 1e5×: "Palace order 1 matches ours at degree 1 on E at six fixed phantom points and on phantom electric energy to ≤ ⟨x⟩ per arm (P1); the rod's near-tip field effect (PEC against natural wall, one mesh) exceeds the two codes' disagreement by ≥ ⟨y⟩×. Record mesh, no 1 g or accuracy claim."
- **§10 I12 is done**, as a provisional step for the weekly to ratify or revert.
- **`ANS-10`** has no §7 row. The weekly reads this as its first **field-level** Palace reading (agree / disagree / inconclusive on the near-tip effect at matched order). The Q-accuracy clause at §2 ≈ line 437 stays.
- §2 ≈ line 437 may gain one clause ("the near-tip field agrees with Palace at matched order").

**Negative result.**
- **Pre-check miss, or (X0) red:** item 235's artifacts are not the fixture. BLOCKED, and the review decides whether to re-export.
- **Palace refuses `Probe`** (no GSLIB in the image is the likely cause; the Dockerfile does not set `PALACE_WITH_GSLIB`, `docker/palace.Dockerfile:40–50`): record it verbatim, open a known-issues entry, and run the `imp3e_*` fallback. The item lands with the energy readings only, and the row sentence drops its point clause. If `Energy` is refused too: BLOCKED (rule (d)).
- **(H) unparseable:** print the header, BLOCKED pending a header ruling.
- **(P) or (RP) red:** report every number, park, BLOCKED. Never widen.
- **Ours (F), (R), (i), (ii), (E), (Kw) or (V) red:** park on `attempt/IMP-3-step4-<ts>`, open a known-issues entry, BLOCKED.
- **(T_E) or (T_W) without the rod > 1e-6, with (RP) and (R) green:** this is not a stop. It is the first field-level disagreement at matched order, on the arm without the rod. Print it per point and component beside (B) and the conj reading, and report. The candidates are cell choice at the point, Palace's probe interpolation, or the convention. No record sentence; the review reads it.
- **(T_E) green without the rod but > 1e-6 on a rod arm:** a reading located at the PEC hole. Report it.
- **Box-speed gate trips:** a partial record (Palace only), BLOCKED (slow box). The review re-queues W6.
- **A printed reading outside its prediction:** report it; no status effect.
- There is no `PREDICTED RED` line.

**Outcome, 2026-10-05 06:00 CDT slot — ✅ DONE, all six windows green.** New files are `scripts/probes/imp3_step4_palace_fields.py` and `tests/validation/test_imp3_palace_fields.py`. Neither earlier probe was edited. Item 235's artifacts were read only, and nothing was overwritten. The fallback configs `imp3e_*` were written but not run, because Palace accepted `Probe` (it also writes `probe-B.csv`).
- **W1** `20261005T110133Z_IMP-3-item242-W1-configs.log` (4 s, `-n 1`): (X0) 117 206 / 168 507 tets equal the meta counts exactly; triangle attributes {1, 21–24} and {1, 21–24, 402}; the PEC lists are [1] / [1, 402] / [1].
- **Palace v0.18.1, `-np 8`, all log checks met** (elements equal the (X0) count, `Rs = 4.375e+01 Ω/sq` on 21–24, 6.400e-02 GHz, four `GMRES solver converged` lines at 31–32 its):
  - W2 off `20261005T110142Z`: Total 25.782 s, peak 2.2 G, ND 140 460. Box-speed r = 1.058.
  - W3 rod `20261005T110223Z`: 38.333 s, 3.0 G, ND 204 156.
  - W4 rodnat `20261005T110306Z`: 37.319 s, 3.0 G.
  - The probe header form is `Re{E_x[<probe>][<excitation>]} (V/m)`. The energy form is `E_elec[<domain>][<excitation>] (J)`. The `read` regex was fitted to the W2 header before W5.
- **W5** `20261005T110353Z_IMP-3-item242-W5-read.log` (2 s): (H) bijections PASS; (P) PASS on `imp3p_off` / `imp3p_rod` (rodnat printed, PASS). (RP) against item 235's S read 8.2e-14 / 1.3e-13 / 1.7e-13.
- **W6** `20261005T110512Z_IMP-3-item242-W6-ours-fields.log` (179 s, `-n 2`, 13 passed):
  - (F) cells equal the meta counts exactly, N_i = 5.6e-5;
  - (R) 1.5e-13 / 3.5e-15 / 2.5e-15;
  - (V) all True;
  - (E) 0.0 on every drive;
  - (Kw) P1 1.297e-4;
  - (B) 1 owned hit at every point on both meshes, so no point is face-ambiguous.
- **(T_E)** (P1; off / rod / rodnat): 1.32e-8 / 1.46e-8 / 8.3e-9, ≤ 9.2e-8 over P1–P4. The conj form reads 0.49 / 0.25 / 0.23, so the convention is confirmed. g = −7.071067853e-2 + 6e-11 j, i.e. 1/(10 √2) and real.
- **(T_W)**: 1.2e-9. q_off/|g|² = 2.896e9.
- **η_wall P1** (identical in both codes, d_η/η ≤ 5e-9): T1 4.053, T2 2.370, T3 1.573, S1 0.439, F1 1.255, F2 0.838. **η_rod P1**: T1 3.782, T2 2.359, T3 1.461, S1 0.454, F1 1.051, F2 1.461.
- **ΔW P1**: rod/rodnat +4.799 %, rod/off +3.809 %. The code gap is ≤ 1.1e-9 on every drive.
- **Closed form** |E_cf|/E₀ against η_rod(p)/η_rod(F1): ratios 0.71 / 0.74 / 0.90 / 0.85 at T1 / T2 / T3 / S1, inside [0.3, 3].
- **Negative controls:**
  - (C1) tip effect / code gap at T1 = 3.053 / 5.85e-9 = **5.2e8×** (≥ 1e5; ceiling ≈ 2.4e9).
  - (C2) 11.99× in both codes (≥ 5).
  - (C3) arm swap 0.750, i.e. 5.7e7× over the off tie (≥ 0.3).
- **Outside prediction (reported, no status effect):**
  - η_wall P1 at F1 / F2: |η − 1| = 0.255 / 0.162, against [3e-3, 0.1].
  - η_rod at F2: 0.461, against ≤ 0.3.
  - So the "far" points are not dipole-far on this record mesh: a 1 mm rod in 537 phantom cells perturbs the field well beyond the PEC-spheroid estimate. That is a reading about placement and resolution, not a stop.
- **Reading:** the near-tip field and the phantom electric energy agree between Palace order 1 and ours at degree 1 on all three arms, to ≤ 1.5e-8 (field, P1) and ≤ 1.2e-9 (energy). The rod's tip effect exceeds the codes' disagreement by 5.2e8×. This is not a 1 g SAR, h- or p-statement.

### §9 item 252 — 2026-10-05 (chain I13, daily provisional)

*Drafted by a general-purpose agent for the 2026-10-05 10:30 daily review; the review reviews it.* Licence: §10 chain I13 **(daily, 2026-10-05 — provisional)** under rule (6) (PROJECT_PLAN.md ≈ lines 3705–3720). This is the second provisional step on chain I since the 10-03 weekly (I12 was the first), so it is inside the three-per-chain cap. It is the next step toward chain I's stated target, Palace. Item 242's journal named it (`docs/testing/attempts.md:2402`): "`ANS-10`'s tip-ratio accuracy is now purely a discretisation question … Palace order 2 on these meshes". Degree 2 is the production order (operator order rule, 2026-10-01). No band, tolerance, Done-when or validation target moves.

**`IMP-3` step 5: I12 at matched order 2.** Palace order 2 against ours at degree 2, on item 235's exact-rectangle meshes, at 64 MHz, with P1 as the drive of record. The comparison is E at item 242's six points plus the phantom electric energy. Each code's S on the arm without the rod must reproduce that code's own order-2 record on the same mesh (`f05x`). Palace's identities are asserted. Each code's order 1 → 2 shift of the near-tip field is printed beside the gap between the two codes. **Target: Palace — a reading, not a gate**, as in I9 and I12.

**Two arms, not three (slot price, below).** The item runs the arm without the rod (`off`) and the rod PEC arm (`rod`). It **drops the natural-wall arm (`rodnat`)**. The three-arm plan prices at ≈ 60 min predicted, which is over the slot. The off arm stays because:
- it is the only arm with order-2 records in both codes;
- it is the cheapest Palace run;
- it is the reference arm for g.

So the order-2 tip effect is read as η_rod (two meshes, remesh included), not η_wall. Adding `rodnat` at order 2 costs ≈ 8 min of Palace plus ≈ 2 min of ours. That is a later step for the weekly to price.

- **Executor:** implementer. One new probe and one new test module. No `src/`, and no edit to any existing file.
- **Tier:** heavy. Configs `-n 1`; Palace `-np 8`; ours **`-n 8`**. Use the complex container for every non-Palace window.
  - Why `-n 8` and not `-n 2`: a single degree-2 solve of 482 158 dofs took 291.7 s at `-n 2` (`WF-15` cost probe, `attempts.md:2419`). This fixture has 753 062 / 1 090 614 dofs (below), so `-n 2` cannot fit 590 s.
  - The rank-local reductions are item 242's code, already exercised at `-n 2` (`…110512Z_IMP-3-item242-W6-ours-fields.log`).
  - The (R2) record was measured at `-n 8` (`ours_f05x_p2.json` `"ranks": 8`).
- **Class:** chain step (solve).
- **Cost:** ≈ 9 measured / ≈ 50 predicted min (15 fixed + ≈ 12 code + ≈ 19 compute + ≈ 4 spare). Price arithmetic below.
- **Depends on:** nothing. Item 235's and item 242's gitignored artifacts in `logs/imp3/`, and items 166 / 194's in `logs/ops61/`, are present on this box (pre-check below).

**Price, from measurement** (each per-run figure is divided by the solves actually in its interval; item 237's lesson).
- **Palace off, order 2 (W2): ≈ 285 s, range 265–300 s.**
  - Basis: `f05x` order 2 at 64 MHz is the same 117 206-cell mesh. One run (one frequency, four excitations) at `-np 8` took Total 265.630 s, peak 6.1 G, ND (p = 2) 753 062, with 37 GMRES iterations per excitation (`20261003T110249Z_OPS-61-item194-W4-palace-o2.log:110,164,351,368`).
  - Box factor: item 242's order-1 off run with `Probe` read 25.782 s, against `f05x` order 1's 23.35 s (×1.104), and against item 235's 24.364 s (×1.058).
- **Palace rod, order 2 (W3): ≈ 460 s, range 400–520 s, peak ≈ 9–11 G.**
  - Split `f05x` order 2 into Preconditioner 149.377 s and the rest 116.25 s (`…item194-W4-palace-o2.log:343,351`).
  - Rod against off at order 1, on these meshes (item 242): Preconditioner 7.097 / 4.179 = ×1.698 and Total 38.333 / 25.782 = ×1.487 (`20261005T110223Z_IMP-3-item242-W3-palace-rod.log:320,327`; `20261005T110142Z_IMP-3-item242-W2-palace-off.log:319,326`).
  - The order-2 preconditioner grew faster on the hole mesh: Preconditioner 88.158 s at ND 570 332 (`20261005T021936Z_OPS-61-item234-W4-palace-hole-pec-o2.log:110,316`), against 149.377 s at 753 062, i.e. exponent ≈ 1.9 in ND.
  - Low: 149.377 × 1.698 + 116.25 × 1.487 ≈ 427 s. High: 149.377 × 1.444^1.9 (= 2.01) + 173 ≈ 473 s.
  - Multiply by the box factor ≈ 1.06.
  - Item 241's copper hole run at order 2 (175.95 s, 4.5 G, `20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log:327,344`) is the same 80 798-cell hole mesh as item 234. It confirms the order-2 price is stable run to run (+9.6 % with the conductivity term).
- **Ours, degree 2, both arms in one window (W5): ≈ 290 s, range 240–450 s.**
  - Off sweep: one four-drive sweep on the 117 206-cell mesh at `-n 8` printed "four driven solves in 64.14 s wall" at 64 MHz (`20261003T094036Z_ANS-4-item192-W3-p2.log:2080`). That window holds **four** such sweeps (10 / 64 / 128 / 128 κ, `:2031,2080,2129,2178`; 429 s elapsed), so 429 / 4 is **not** the per-sweep price. 64–73 s is.
  - Item 166's sweep read 72.33 s (`20261002T095914Z_OPS-61-F10b-ii-W1-ours-p2.log:2031`).
  - Rod sweep: 1 090 614 / 753 062 = ×1.448 dofs. The scaling exponent from item 192 (70.70 s at 117 206 cells, 10 MHz) to item 212 W1 (186.93 s at 259 524 cells, 10 MHz, 36.12 GiB summed; `attempts.md:2145`) is ln 2.644 / ln 2.214 = 1.23. That gives 64.14 × 1.438^1.23 ≈ 100 s (range 95–135 s), and ≈ 20 GiB summed.
  - Builds: off mesh 23.90 s (`…item192-W3-p2.log:2031`); rod mesh ≈ 47 s (`IMP-3` W1 `:3646`, item 235 brief).
  - Total: environment tests ≈ 10 + 24 + 64–72 + 47 + 100–135 + readings ≈ 15 ⇒ ≈ 260–305 s.
- **Configs (W1)** ≈ 5 s; **read (W4)** ≈ 3 s (item 242: 4 s / 2 s).
- **Compute** ≈ 5 + 285 + 460 + 3 + 290 ≈ 1 043 s, plus ≈ 1 min of container overhead, ≈ **18.5 min**. The tops of every range give 1 278 s ≈ 21.3 min. The slot is then ≈ 15 + 12 + 21 ≈ 48 min, with spare left.
- **The third arm does not fit.** `rodnat` would add one more Palace rod run (≈ 460 s) and one more rod sweep (≈ 110 s), ≈ +9.7 min, for ≈ 60 min predicted.

**Feasibility, settled from the repo.**
- **The order knob is `Solver.Order`.** `ops61_step1a_palace_ans4.py:734–741` writes `f05x_o2` by asserting `Order == 1` on the source and setting it to 2. Nothing else changes except `Output`. That run solved at ND (p = 2) 753 062 (`…item194-W4-palace-o2.log:110`).
- **Palace's `Probe` / `Energy` output and its header form are known** (item 242 outcome; `imp3_step4_palace_fields.py:156–206` parses them). The parsers are imported, never re-written.
- **The two meshes are `f05x`'s and item 235's.** Palace `imp3x_off` against `f05x` tied at 2.9e-14 at order 1 (item 235 W5). Ours without the rod against `ANS4_EXACT_STRIP_S_P1_RECORDS["64"]` tied at ≤ 1.2e-10 (item 235 W6). So the order-2 records on `f05x` are this fixture's off-arm records:
  - Palace: `logs/ops61/postpro_f05x_o2/port-S.csv`;
  - ours: `logs/ops61/ours_f05x_p2.json`, `"degree": 2, "ranks": 8, "cells": 117206`, key `s_kappa_off`. This is the file `ANS-4`'s (X) reproduces at 8.9e-13 (`attempts.md:1951`).
- **The ND (p = 2) count is predictable exactly.** It is 2E + 2F, with F = 2·cells + F_bdry/2.
  - Off: E = 140 460, F = 2 · 117 206 + 3 318/2 = 236 071, giving **753 062**. This equals item 194's count, so the formula checks.
  - Rod: E = 204 156, F = 2 · 168 507 + (3 316 + 4 958)/2 = 341 151, giving **1 090 614**.
  - DOLFINx N1curl degree 2 spans the same space, so our global dof count is predicted equal.

**Code (new files only; additive; disclosed).**
- **New probe `scripts/probes/imp3_step5_palace_o2.py`.** It loads `imp3_step4_palace_fields.py` by `importlib` for `_step3`, `_msh_census`, `parse_probe`, `parse_energy`, `PROBE_NAMES` and `PROBE_POINTS`. Through `_step3()._ops61()` it gets `_p_check`, `_read_port_s` and `_fmt`. **Edit no earlier probe.** Dispatch is `configs` | `read`, nothing else.
  - **`configs`** (serial; `assert comm.size == 1`):
    - **(X0)** Re-run item 242's census: tets 117 206 / 168 507 equal `imp3x_meta.json` exactly, and the triangle attributes are {1, 21–24} / {1, 21–24, 402}.
    - Load `imp3p_off.json` / `imp3p_rod.json` read-only. Assert `Solver.Order == 1`, PEC [1] / [1, 402], `Domains.Postprocessing` keys {Energy, Probe}, and Output `postpro_imp3p_<arm>`.
    - Deep-copy each, set `Solver.Order = 2` and `Problem.Output = "postpro_imp3q_<arm>"`, and write `imp3q_off.json` / `imp3q_rod.json` with mode `"x"`.
    - **(X0q, printed)** List every key path where `imp3q_off.json` differs from `logs/ops61/f05x_o2.json`, excluding `Model.Mesh`, `Problem.Output` and `Domains.Postprocessing`. Predicted: empty. A non-empty list is printed verbatim and is not a stop, because (RP2) below is the asserted test of the same fact.
  - **`read`:**
    - (H) Run `parse_probe` / `parse_energy` on each `postpro_imp3q_*` (bijection or stop).
    - Run `_p_check` on `imp3q_off` / `imp3q_rod` (fails asserted).
    - **(RP2)** Below.
    - Print Palace's own order 1 → 2 shifts against `logs/imp3/palace_imp3p_fields.json` (read-only).
    - Save the parsed arrays to `logs/imp3/palace_imp3q_fields.json` (mode `"x"`).
    - Never read `z_matrix_ohm`, `docs/private/` or `aed_results/`.
- **New module `tests/validation/test_imp3_palace_fields_o2.py`** (`@complex_only`; one test, `test_imp3_palace_fields_degree2`):
  - Builds: `probe.build_arm("off")`, then `build_arm("rod")` (`GEO-23` order), each through `build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=…, build_only=True)`.
  - Sweeps: its **own** wrapper `_sweep2(case)` = `run_n_port_sparameter_sweep(case["problem"], case["port_defs"], lumped_sheet_ports=case["specs"], lumped_sheet_facet_tags=case["facet_tags"], keep_fields=True, degree=2)`. `test_imp3_palace_twin._sweep` hard-codes `degree=1` (`:81–89`) and must not be used. κ stays off.
  - Import helpers only, never a `test_*` name:
    - `_ids` and `_sheet_areas` from `test_imp3_palace_twin`;
    - `_owned_hits`, `_pal_fields`, `PTS`, `PNAMES` and `FP` from `test_imp3_palace_fields`;
    - `_matrix_quantities`, `_wall_dofs`, `_wall_max_e`, `PEC_WALL_MAX_E` and `PORT_IDS` from `test_imp3_rod_birdcage`.
  - **A local `_readings2`**: item 242's `_readings` with `quadrature_degree` **4**. |e|² of a degree-2 field is degree 4, so quadrature 4 is exact on affine cells. Item 242's pinned 2 was exact only at degree 1.
  - Per arm, before `del res`, for every drive:
    - E at `PTS` via `evaluate_vector_field_parallel(e_complex, PTS, owner="lowest_original_cell")`;
    - W_ph, allreduced;
    - (B);
    - the global dof count `V.dofmap.index_map.size_global * V.dofmap.index_map_bs`.
  - Write `logs/imp3/ours_imp3q_fields.json` on rank 0 (mode `"w"`, a new name).
  - Read `palace_imp3q_fields.json` (order 2), and item 242's `palace_imp3p_fields.json` and `ours_imp3p_fields.json` (order 1), all read-only, and print the readings below.
  - **Never write** `ours_imp3p_fields.json`, `ours_imp3x_p1.json` or anything under `logs/ops61/`.

**Anchors (asserted; bands imported, never restated).**
- **(X0)**, W1: exact tet counts and attribute sets, as item 242 (X0) measured them.
- **(K2) The order reached Palace's solve**, read by the executor from the W2 / W3 logs and journalled:
  - W2's `ND (p = 2)` equals **753 062** exactly. This is backed by the same mesh at the same order (`…item194-W4-palace-o2.log:110`).
  - W3's `ND (p = 2)` is printed, predicted 1 090 614.
- **Palace log check, W2 / W3:**
  - `elements` equals the (X0) count;
  - `Rs = 4.375e+01 Ω/sq` on 21–24;
  - 6.400e-02 GHz;
  - four `GMRES solver converged` lines.
- **(P) Palace, W4**, on `imp3q_off` / `imp3q_rod` via `_p_check` (`ops61_step1a_palace_ans4.py:742–782`):
  - `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`);
  - `PASSIVITY_SIGMA_TOLERANCE`, with column power ≤ 1 + 1e-9 (`test_port_birdcage_four_port.py:129`);
  - `ADJACENT_SPREAD_BAND` (`test_port_birdcage_lumped_column.py:121`);
  - `POOLED_SEPARATION_FLOOR` (`test_port_birdcage_four_port.py:173`).
  - Backed at order 2 on the off mesh by item 194 (P) green (`20261003T110728Z_OPS-61-item194-W5-read.log`). On the rod arm it is backed only at order 1 (item 235 W5).
- **(RP2) Palace reproduces its own order-2 record:** `imp3q_off` S against `logs/ops61/postpro_f05x_o2/port-S.csv`, entrywise max|ΔS| / max|S|, within `REUSE_REPRODUCTION_RTOL` 1e-6 (`tests/validation/test_port_birdcage_ring_matrix.py:171`). This is backed by the same comparison at order 1, 2.9e-14 (item 235 W5). Postprocessing does not enter the solve. There is no rod-arm record at order 2, so the rod arm has no (RP2).
- **Ours, W5:**
  - (F): cells equal `imp3x_meta.json` exactly, and N_i equals 5.6e-5 within `STRIP_BAND` (`test_ops61_exact_strip_ports.py:58`).
  - **(R2)** Ours without the rod at degree 2 reproduces `logs/ops61/ours_f05x_p2.json` `s_kappa_off`, entrywise, within `REUSE_REPRODUCTION_RTOL`.
    - Assert the file's `"ranks"` equals `comm.size` (8) and `"degree"` equals 2. A mismatch is a stop, not a skip, because (R2) is this item's anchor.
    - Backing: the same builder coincidence at degree 1 (≤ 1.2e-10, item 235 W6), and `ANS-4` (X) at degree 2 against this file (8.9e-13, item 192 W3).
  - (i)/(ii) on both sweeps (`_ids`).
  - (E): rod PEC wall max |E| ≤ `PEC_WALL_MAX_E` (`test_imp3_rod_birdcage.py:135`) on every drive, on the degree-2 space's wall dofs (`_wall_dofs(msh, built[2], V)`).
  - (V): `valid_mask` all True for the six points on both arms.

**Printed, never gated.** P1 is the reading of record; P2–P4 are printed beside it. I12's (T_E) / (T_W) ≤ 1e-6 were **predictions**, printed and never asserted (item 242 brief, "Printed, never gated"). They are not code-agreement bands, so here they stay predictions.
- **(D)** Our global dof count per arm, against Palace's ND (p = 2). Predicted equal: 753 062 / 1 090 614.
- **(B)** Owned-cell hits per point. Predicted 1 everywhere, as item 242 measured on these meshes.
- **(T_E2) The field tie at order 2.** g₂ = E^ours_c*(F1) / E^Pal_c*(F1) on the off arm, with c* = argmax over c of |E^Pal_c(F1)|. Per arm, the reading is the max over (p, c) of |E^ours − g₂·E^Pal| / ‖g₂·E^Pal(p)‖, with item 242's exclusions.
  - **Predicted ≤ 1e-6.** The derivation is I12's: same mesh, same order, same operator, Palace `Tol` 1e-10. Order 1 measured ≤ 1.5e-8.
  - Print g₂ against item 242's g = −7.071067853e-2. It is the ratio of incident scales, which does not depend on order: predicted |g₂/g − 1| ≤ 1e-6.
  - Print the conj form, predicted O(1).
- **(T_W2) The energy tie.** q_rod / q_off − 1, predicted ≤ 1e-6. Also q_off / |g₂|² against item 242's 2.896e9 (a convention constant), predicted equal to ≤ 1e-6 relative.
- **(T_S2) The S tie per class, per arm** (`_circulant_classes`). Predicted ≤ 1e-6. On the off arm it is implied by (RP2) + (R2) + item 194's 0.0000 %.
- **The order 1 → 2 shift, per code (the item's reading for `ANS-10`).** Same code and same incident scale at both orders, so the ratios are meaningful:
  - δ_p = |E^o2(p)| / |E^o1(p)| − 1, per arm, at each point. Predicted |δ| ∈ [1e-3, 1] at T1–T3 and S1, and [1e-3, 0.3] at F1 / F2. For scale, the S class shift on the off arm is 6.667 / 2.544 / 4.522 % (item 194).
  - η_rod^o2(p) / η_rod^o1(p) − 1 at every point. Predicted [1e-3, 0.5] at T1.
  - ΔW(rod/off) at order 2, beside order 1's +3.809 %. Predicted in [1e-3, 0.1].
  - Per class, the S shift for both arms.
  - Ours against Palace on each shift: |δ^ours − δ^Pal|, predicted ≤ 2e-6 (two ties).
- **Closed form beside it:** `pec_prolate_spheroid_field` (`utils/analytical.py:1576`) against η_rod^o2(p) / η_rod^o2(F1) at T1–T3 and S1. Predicted [0.3, 3], with no accuracy claim, as in I12.
- **Palace prices:** Total, Preconditioner, peak memory and ND per run.

**Negative controls (rule (e)).** No order-2 field comparison exists on this fixture, so every factor below is **predicted** and printed, with its ceiling computed. (K2), (RP2) and (R2) are the asserted knobs.
- **(C1′) The order shift of the tip effect against the matched-order code gap.** Separation = |η^o2_rod(T1) / η^o1_rod(T1) − 1| / max(d_η^o1(T1), d_η^o2(T1)) / η.
  - **Predicted ≥ 1e3.** Least case: 1e-3 / 1e-6.
  - **Ceiling ≈ 1e8:** a shift ≤ 1 against I12's measured d_η/η ≤ 5e-9, doubled.
  - Below 1e3, the two codes do not agree on the order shift at the tip: report it.
- **(C2′) The arm swap at order 2.** |E^ours,o2_off(T1) − g₂·E^Pal,o2_rod(T1)| / ‖g₂·E^Pal,o2_rod(T1)‖.
  - **Predicted ≥ 0.3.** Order 1 measured 0.750 (item 242), but that is a different order, so the factor is predicted, not asserted.
  - **Ceiling ≈ 1** relative, i.e. ≈ 1e8× over (T_E2) at the round-off floor.
- **(C3′) The order swap: the probe reading sees the order.** On the rod arm at T1: |E^ours,o2 − g·E^Pal,o1| / ‖g·E^Pal,o1‖, using item 242's Palace array.
  - **Predicted ≥ 1e-3**, so ≥ 1e3× over (T_E2).
  - **Ceiling** ≈ |δ_p(T1)| / 1e-8 ≤ 1e8×.
  - It catches a Palace run that ignored `Order` in the solve but not the ND count, and an ours that silently ran degree 1.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`).
- **Pre-check** (not a harness call): `ls logs/imp3/ logs/ops61/postpro_f05x_o2`.
  - Required:
    - `imp3x_off.msh`, `imp3x_rod.msh`, `imp3x_meta.json`;
    - `imp3p_off.json`, `imp3p_rod.json`;
    - `palace_imp3p_fields.json`, `ours_imp3p_fields.json`;
    - `logs/ops61/f05x_o2.json`, `logs/ops61/ours_f05x_p2.json` and `logs/ops61/postpro_f05x_o2/port-S.csv`.
  - Absent: every `imp3q_*`, `postpro_imp3q_*` and `*_imp3q_fields.json`.
  - Any miss: stop and report.
- **W1, configs** (≈ 5 s):

  `scripts/testing/run_and_log.sh IMP-3-item252-W1-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/imp3_step5_palace_o2.py configs'"`
- **W2, Palace without the rod, order 2** (≈ 285 s). This is **the box-speed gate**:

  `scripts/testing/run_and_log.sh IMP-3-item252-W2-palace-off-o2 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/imp3 && timeout -k 30 590 palace -np 8 imp3q_off.json'"`
- **Box-speed gate (binding, after W2):** r = W2 Total / 265.630 s. W3 projects to 1.78 × W2 Total (the high price, 473 / 265.63).
  - r ≤ 1.18 (W3 ≤ 560 s): run everything.
  - r > 1.18: stop after W2, journal W2's Total / ND / iterations, and mark BLOCKED (slow box). The review re-queues. Never move W3 to a longer timeout or to `-np 12` in-slot: neither is priced.
- **W3, Palace with the rod, order 2** (≈ 460 s, range 400–520; peak ≈ 9–11 G). The command is W2's with `imp3q_rod.json`, labelled `IMP-3-item252-W3-palace-rod-o2`. After any kill, check `docker compose -f docker/docker-compose.yml ps` for a lingering `fem-em-palace` container before anything else.
- **W4, read** (≈ 3 s). This is the probe's last window: no probe edit after it starts.

  `scripts/testing/run_and_log.sh IMP-3-item252-W4-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/imp3_step5_palace_o2.py read'"`
- **W5, ours at degree 2** (≈ 290 s, range 240–450). This is the closing window: no module edit after it starts.

  `scripts/testing/run_and_log.sh IMP-3-item252-W5-ours-p2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_imp3_palace_fields_o2.py -v -s --tb=short'"`

**Traps already paid for.**
- **Degree is hard-coded in item 235's `_sweep`** (`degree=1`). Importing it would silently run order 1 against Palace order 2, and (T_E2) would read the order shift, not a tie. Use `_sweep2`. (D) and (C3′) are the tell.
- **Quadrature pinned for the degree.** Use `quadrature_degree=4` on the degree-2 energy form (POST-5's rule; item 242's 2 under-integrates here and would put a quadrature error into (T_W2)). Do no UFL or `sqrt` comparisons; take ratios in numpy.
- **Pointwise N1curl depends on the MPI width without the owner rule** (known-issues 2026-09-27). Always pass `owner="lowest_original_cell"` and evaluate `e_complex` once. Keep `PTS` bitwise identical on every rank (it is imported). Never call `f.eval(points, np.arange(n))`.
- **The owner rule raises on `create_submesh` meshes** (known-issues 2026-10-04). It does not apply here (full mesh, item 242). If it raises anyway, stop: never fall back.
- **Never overwrite earlier artifacts.**
  - Write with mode `"x"` for the configs and Palace arrays.
  - Use new `imp3q_` names.
  - `logs/ops61/` and every `imp3x_*` / `imp3p_*` file are read-only.
  - Do not call any code path that writes `ours_imp3x_p1.json` or `ours_imp3p_fields.json`.
- **Tag 402 and the cell tags live in the build tuple** (`built[2]` / `built[1]`), not in the sweep's.
- **Build order and memory.** Build without the rod first (`GEO-23`). Evaluate before `del res`, and free the off arm's fields, case and mesh before the rod build.
- **Rank-local quantities.** `assemble_scalar`, collision counts and wall dofs are rank-local: reduce them. Keep collectives outside any `rank == 0` branch. At `-n 8`, a rank-0-only collective deadlocks the window (item 227's W1).
- **Palace CSV headers are parsed, never assumed.** Reuse item 242's parsers. An unparseable header is a stop.
- **A `GMRES … did not converge` line is a finding.** That run's fields are not measured. Never change `Order`, `Tol`, `MaxIts` or `-np` in-slot.
- **Divide window prices by the solves in them** (item 237). The 64.14 s is per sweep. Item 192 W3's 429 s holds four sweeps.
- **Harness discipline:**
  - `-s`, `tests/environment` first, `FEM_EM_REQUIRE_COMPLEX=1`, foreground, no pipes, no `-k`;
  - Status 124: run `pgrep -c python3` in the container first and never lengthen the timeout. A killed compile leaves an orphan `.c`: run the rubric's cache sweep only when nothing is running;
  - rule (i): no edit after a file's last window.
- **Palace numbers go in at full detail** (version, config, mesh). No AED number, and no Palace figure set against one, in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace order 2 and ours at degree 2, 64 MHz, on item 235's exact-rectangle record meshes (117 206 / 168 507 cells, σ = 800 S/m untuned coil, 50 Ω, κ off), on two arms (without the rod, and rod PEC). The comparison is E at six fixed points and phantom electric energy, P1 of record. The order 1 → 2 shift of the near-tip field is printed per code. Agreement is a **matched-order check at the production order** between two independent codes. The order shift is a **one-mesh p-reading**, not an h-statement and not a convergence rate (two orders, one mesh). Not covered:
- **the 1 g tip SAR, the ratio Q and its accuracy.** These need the rod + strip build on `PHANTOM_RESOLUTION_1G_RUNG` (never built) or an h-ladder: a later I-step for the weekly to price;
- the natural-wall arm at order 2 (η_wall), dropped for price;
- absolute |E| / SAR;
- tissue, heating, compliance, ISO 10974 (`WF-9`), ASTM F2182 (7c);
- the tuned or copper coil, F-human, other rod placements, the kept-cells rod;
- AED, `src/`, any band.

**Status it can move.**
- **§7 `IMP-3` stays ✅** (no glyph; rule (j) read: the row's Done-when, "identities green … on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture", is met by step 1 and says nothing about order). The review appends a **step-5 record sentence** beside steps 3 and 4. Candidate text, valid only if W4 and W5 are green, (T_E2) ≤ 1e-6 on both arms and (C1′) ≥ 1e3×: "Palace order 2 matches ours at degree 2 on E at six fixed phantom points and on phantom electric energy to ≤ ⟨x⟩ per arm (P1; off, rod PEC), each code reproducing its own order-2 `f05x` record; the order 1 → 2 shift of the rod's near-tip effect, η_rod(T1), is ⟨y⟩ % in both codes, agreeing to ⟨z⟩. Record mesh, one-mesh p-reading, no 1 g or accuracy claim."
- **§10 I13 is done**, as a provisional step for the weekly to ratify or revert.
- **`ANS-10`** (no §7 row; §10 ≈ line 4460): the weekly reads this as the first **order-2** field-level Palace reading, and as the first measured size of the tip ratio's order sensitivity on the record mesh. If |η_rod^o2(T1) / η_rod^o1(T1) − 1| is small (say < 1 %), the order-1 tip ratio is p-stable on this mesh. If it is large, `ANS-10`'s accuracy question needs the 1 g rung or an h-ladder before anything else. The weekly decides; this item sets no threshold.
- §2 ≈ line 437 keeps "the ratio's accuracy is `ANS-10`'s". It may gain one clause: "the near-tip field agrees with Palace at matched order 1 and 2; its order shift is ⟨y⟩ %".

**Negative result.**
- **Pre-check miss, or (X0) red:** not the fixture. BLOCKED; the review decides.
- **(K2) red** (W2's ND ≠ 753 062): Palace did not solve at order 2 on `f05x`'s mesh. Stop before W3, BLOCKED.
- **A Palace abort or a non-converged GMRES:** record it verbatim, open a known-issues entry, BLOCKED (rule (d)). Never re-run with changed solver settings.
- **(P) or (RP2) red:** report every number, park, BLOCKED. Never widen.
- **Ours (F), (R2), (i), (ii), (E) or (V) red:** park on `attempt/IMP-3-step5-<ts>`, open a known-issues entry, BLOCKED. (R2) red at matched width means ours at degree 2 did not reproduce item 166, which is a finding about degree-2 reproducibility, not about the rod.
- **(T_E2) or (T_W2) on the off arm > 1e-6, with (RP2) and (R2) green:** not a stop. It is the first field-level disagreement at order 2. Print it per point and component beside (B), (D) and the conj form, and report. The candidates are Palace's point evaluation of a p = 2 field, our degree-2 evaluation, or the energy quadrature. No record sentence; the review reads it.
- **(T_E2) green on the off arm but > 1e-6 on the rod arm:** a reading located at the PEC hole. Report it.
- **(D) unequal:** the two codes discretise different spaces. Report the counts; with (RP2) / (R2) green it is a counting-convention reading, not a stop.
- **Box-speed gate trips:** a partial record (W2 only), BLOCKED (slow box). The review re-queues.
- **A printed reading outside its prediction:** report it; no status effect.
- There is no `PREDICTED RED` line.

**Outcome (2026-10-05 12:00 CDT slot, second item, implementer): GREEN, W1–W5 all Status 0, first run.**
- **W1** `20261005T172148Z_IMP-3-item252-W1-configs.log` (1 s): (X0) 117 206 / 168 507 tets, attributes {1, 21–24} / {1, 21–24, 402}; `imp3q_off.json` / `imp3q_rod.json` written (Order 2); (X0q) 0 differing key paths against `f05x_o2.json`.
- **W2** `20261005T172153Z_IMP-3-item252-W2-palace-off-o2.log` (263 s): Total 261.083 s, Preconditioner 145.112 s, peak 6.7 G, **(K2) ND (p = 2) 753 062** exact, 4 × GMRES converged in 37 its, 117 206 elements, Rs 4.375e+01 on 21–24, 6.400e-02 GHz. Box-speed gate r = 261.083 / 265.630 = 0.983 ≤ 1.18: run everything.
- **W3** `20261005T172625Z_IMP-3-item252-W3-palace-rod-o2.log` (432 s): Total 429.928 s (priced 400–520), Preconditioner 256.740 s, peak 9.8 G, ND (p = 2) 1 090 614 (predicted exactly), GMRES 49 / 49 / 47 / 48 its, 168 507 elements.
- **W4** `20261005T173342Z_IMP-3-item252-W4-read.log` (2 s): (H) bijection; (P) PASS on `imp3q_off` and `imp3q_rod`; **(RP2) 7.852e-14** ≤ 1e-6.
- **W5** `20261005T173453Z_IMP-3-item252-W5-ours-p2.log` (307 s, 13 passed at `-n 8`, 305.3 s pytest): (F) cells exact, N_i 5.6e-05; **(R2) 8.670e-13** against `ours_f05x_p2.json` (ranks 8, degree 2); (i)/(ii) both arms; (E) wall max |E| 0.0 on every drive; (V) all True; (B) 1 everywhere; **(D) 753 062 / 1 090 614, equal to Palace's ND**.
- **Ties (printed):** g₂ = −7.071067778e-02 (c* = z), |g₂/g − 1| 1.06e-08. (T_E2) P1 2.625e-08 off (F2.x) / 1.179e-08 rod (S1.z); P1–P4 off ≤ 5.26e-08, rod 1.18e-08 / **1.377e-06** / 4.93e-08 / **1.344e-06** — P2 and P4 on the rod arm sit just above the 1e-6 prediction (a printed reading outside its prediction, no status effect; P1 of record is inside). Conj form 0.306 / 0.317. (T_W2) q_rod/q_off − 1 +9.33e-11; q_off/|g₂|² 2.895920711e9.
- **Order 1 → 2 shift (the `ANS-10` reading):** η_rod^o2(T1) P1 = 3.594534 in both codes (d_η/η 1.7e-09); **η_rod^o2/η_rod^o1 − 1 at T1 = −4.961 % in both codes** (T2 −5.95 %, T3 −8.95 %, S1 −12.17 %); δ_p rod P1 T1 −6.458 %, off T1 −1.575 %, F2 off +54.7 % (outside the [1e-3, 0.3] F prediction: printed, no effect); ours-vs-Palace shift gaps ≤ 1.1e-08. ΔW(rod/off) P1 +4.927 % at order 2 vs +3.809 % at order 1, gap 9.8e-11.
- **Closed form** ratio T1 0.830, T2 0.874, T3 1.097, S1 1.083 (within [0.3, 3]).
- **Controls:** (C1′) 4.961e-02 / 2.598e-09 = 1.91e7× (≥ 1e3); (C2′) 0.726 (≥ 0.3); (C3′) 0.123, 1.04e7× over (T_E2 rod).
- Artifacts (gitignored): `logs/imp3/imp3q_{off,rod}.json`, `postpro_imp3q_{off,rod}/`, `palace_imp3q_fields.json`, `ours_imp3q_fields.json`. §7 `IMP-3` unchanged (✅); record sentence left to the review.

## §7 row moved verbatim 2026-10-05 (rotate_plan_archive.py chunks)

 **A rod in the coil-driven phantom: with / without on the loaded F-small birdcage at 64 MHz (Phase 7 subgoal (7b))** (opened 2026-09-30 by the weekly; §10 chain I6; serial on `GEO-36`). The rod as a PEC hole in the phantom. Asserted, identities only: 4×4 reciprocity and passivity at the imported bands, the exact discrete power partition (`PORT-16`), the mirror symmetry the rod's placement keeps. Printed: the 1 g SAR at the tip with and without the rod, beside `IMP-2`'s closed-form estimate from the local background field. *(Made precise 2026-10-02 by the 10:30 review under rule (j), §9 item 174; ratified by the 2026-10-03 weekly. Fixture: `build_four_port_sweep(frequency_hz=64e6, built=…)` on `GEO-36`'s rod build against the flag-off build. The rod wall's n × E = 0 is the default exterior PEC, asserted. The mirror is σ_y (P2 ↔ P4), asserted on Z's five mirror pairs at `ADJACENT_SPREAD_BAND`. The control leaves the rod wall natural. The 1 g ratio of record is on `PHANTOM_RESOLUTION_1G_RUNG`, because the record mesh's 537 phantom cells cannot hold a 1 g ball.)* **Done-when (§4):** the identities green with the rod in on `GEO-36`'s 167 371-cell record mesh ⇒ 🟡; the same identities green on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture. The ratio's accuracy is `ANS-10`'s (commissioned by the weekly that reads this). No heating or compliance claim. *History: `docs/planning/chunks/IMP-3.md`.* | ✅ **Audit PASS 2026-10-02 18:00 review.** 2026-10-02 (§9 item 174, W1 `20261002T170428Z` 241 s `-n 2` ⇒ 🟡, W2 `20261002T170843Z` 260 s `-n 4` ⇒ ✅): (A0)/(R0 2.98e-10)/(i) ≤ 1.5e-15/(ii)/(P) ≤ 7e-15/(M) ≤ 0.082 %/(E) 0 on 7 437 dofs/(S) 403–671× green on both meshes, C1 natural wall wired (1.3e-4) with (P) closed. Printed on the 1 g rung: tip Q = 4.50 (P1), 3.96 (ccw/cw), 1.58 (P2/P4, below its [2, 10] prediction), Q/estimate 0.86–1.55 — a fixture, no accuracy claim (`ANS-10`). **Step 3 record (§9 item 235, 2026-10-04 22:30 slot; W6 `20261005T035322Z` 172 s `-n 2`):** Palace v0.18.1 order 1 matches ours at degree 1 on the exact-rectangle rod fixture (117 206 / 168 507 cells, κ off) to ≤ 2.2e-10 per class on all three arms (off, rod PEC, rod natural); ΔS(rod) agrees at matched order to ≤ 1.3e-6 and ΔS_wall to ≤ 9.3e-6 per class; the rod's 4×4 effect is ≤ 1.5e-4 of \|S\| (remesh floor). No tip-SAR or discretisation claim. **Step 4 record (§9 item 242, 2026-10-05 06:00 CDT slot; W6 `20261005T110512Z` 179 s `-n 2`):** Palace order 1 matches ours at degree 1 on E at six fixed phantom points and on phantom electric energy to ≤ 1.5e-8 / ≤ 1.2e-9 per arm (P1; off, rod PEC, rod natural); the rod's near-tip field effect (PEC against natural wall, one mesh; η_wall(T1) = 4.05) exceeds the two codes' disagreement by ≥ 5.2e8×. Record mesh, no 1 g or accuracy claim. 

### §9 item 263 — 2026-10-06 (chain I14)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §10 chain I14 (PROJECT_PLAN.md ≈ line 4655), a weekly-enumerated step (§9 rule (4)), posed by the 10-05 weekly beside its ruling bullet "I13 — the tip ratio is not p-stable at 1 %" (≈ line 4490). No band, tolerance, Done-when or validation target moves.

**`IMP-3` step 6 — the tip ratio on the 1 g rung at degree 2.** This runs item 235's exact-rectangle rod fixture at `phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG` (0.0025 m, `tests/validation/test_birdcage_sar_1g_rung.py:114`), at 64 MHz with κ off. The fixture is `_build(True, …)` / `_build_rod(…)` with `leg_sheet_strip_fraction=GATED_WIDTH_FRACTION` (`tests/validation/test_port_lumped_narrowed_sheet.py`), built as `scripts/probes/imp3_step3_palace_rod.py:54–58` builds it. Two arms (off and rod PEC) are solved at degree 1 and at degree 2. The identities are asserted at the imported bands. Printed: η_rod(T1) at both orders, its p-move against I13's −4.961 % on the record mesh, and the 1 g tip SAR ratio Q at both orders. **Target: the identities. The p-stability question is a reading.**

- **Executor:** implementer. One new test module (`tests/validation/test_imp3_tip_ratio_1g.py`). No `src/`, and no edit to any existing file.
- **Tier:** heavy. W1 `-n 1`; W2 / W3 `-n 8`; W4 `-n 4`. Every window runs in the complex container.
- **Class:** chain step (solve).
- **Cost:** ≈ 17 measured / ≈ 55 predicted min (15 fixed + ≈ 15 code + ≈ 18 compute + ≈ 7 spare).
- **Depends on:** nothing. Every module named below is on `main`. The rod + strip + 1 g rung combination has **never been built**: it is new geometry, so W1 builds it serially first.

**Settled from the code.**
- **The modules items 235 / 252 used:**
  - `tests/validation/test_imp3_rod_birdcage.py` (step 1): `_arm` `:295`, `_case` `:503`, `_assert_identities` `:617`, `_control_c1` `:639`, `_wall_dofs` `:188`, `_wall_max_e` `:199`, `_outer_box_tags` `:205`, `_matrix_quantities` `:275`, `_ball_mean_closed_form` `:229`, `_peak_rss_gib` `:498` and `TIPS` `:140`. The 1 g ball goes through `mass_averaged_sar` at `QUADRATURE_DEGREE` 16 (`test_mass_averaged_sar_standard_masses.py:121`).
  - `tests/validation/test_imp3_palace_twin.py` (item 235): `_ids` and `_sheet_areas`.
  - `tests/validation/test_imp3_palace_fields.py` (item 242): `PTS`, `PNAMES`, `FP` and `_owned_hits`.
  - `tests/validation/test_imp3_palace_fields_o2.py` (item 252): `_readings2` `:68` (quadrature 4), and the η_rod / p-move print at `:207–225`.
- **`_arm` cannot run degree 2.** It calls `build_four_port_sweep(…, build_only=True)` without `degree`, then sweeps at `case["degree"]`, whose default is 1 (`test_port_birdcage_four_port.py:219,230`).
  - The new module therefore carries a local `_arm_d(built, label, rod, degree, field_hook)`: `_arm` copied, with `degree=degree` passed to `build_four_port_sweep`.
  - It also carries a local `_case_d(degree, arms)`, which builds with **both** `phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG` and `leg_sheet_strip_fraction=GATED_WIDTH_FRACTION`. Step 1's `_case` passes neither the strip nor a degree.
  - Import helpers and constants only, never a `test_*` name.
- **η_rod is item 252's definition:** η_rod(p) = ‖E_rod(p)‖ / ‖E_off(p)‖ per drive, where ‖·‖ is the Euclidean norm over components.
  - E comes from `evaluate_vector_field_parallel(e_complex, PTS, owner="lowest_original_cell")`, evaluated inside `field_hook` before the fields are freed.
  - The record-mesh values are η_rod^o2(T1) P1 = 3.594534 and η^o2/η^o1 − 1 = −4.961114e-02 (`20261005T173453Z_IMP-3-item252-W5-ours-p2.log:3910–3911`).
- **Q** = SAR_1g,tip(rod) / SAR_1g,tip(off) per tip and drive, which is step 1's reading. Step 1 read P1 top 4.496087 on the **ragged** 1 g rung at degree 1 (`20261002T170843Z_IMP-3-step1-W2.log:3864`).
- **Sizes, from the ragged sibling:**
  - 199 920 / 250 866 cells (`…170843Z_IMP-3-step1-W2.log:1918,3758`).
  - On the record mesh the exact strip added 1 121 cells (116 085 → 117 206, item 235 W1), so the prediction is ≈ 201 000 / ≈ 252 000 cells.
  - Degree-2 dofs are predicted at ≈ 1.3 M / ≈ 1.6 M (2E + 2F, item 252's counting identity).

**One knob, four windows.** The env knob `FEM_EM_IMP3_STEP6` takes one of {`build`, `d2off`, `d2rod`, `d1`}, and one test function dispatches on it (never `-k`).
- Each degree-2 arm gets its own window. One arm at degree 2 is priced at ≈ 240–320 s, so the pair would cross 590 s.
- Each window saves its arrays to `logs/imp3/ours_imp3g_<knob>.json` on rank 0 (mode `"x"`).
- The closing window (W4) reads the two degree-2 files read-only.

**Anchors (asserted; bands imported, never restated).**
- **W1, census (`build`, `-n 1`):**
  - (G0) Both builds complete serially. Off-arm cells are printed against `ONE_GRAM_RUNG_CELL_RECORD` 199 920 (`test_birdcage_sar_1g_rung.py:125`), beside `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py`). This is **printed, not asserted**: it is the first record of the exact-rectangle 1 g rung, and no record exists yet.
  - (G1) Rod cells < `ROD_CELL_GUARD_W2` 320 000 (`test_imp3_rod_birdcage.py:115`), asserted before any solve.
  - (G2) Every N_i area = 5.6e-5 within `STRIP_BAND` 1e-9 (`test_ops61_exact_strip_ports.py:58`), on both arms.
  - (G3) Rod arm: tag-402 facets > 0, all of them exterior. Off arm: 402 count 0.
  - (G4) The 1 g ball fits: the two geometric asserts at the end of `_assert_identities` (`:636–637`).
- **W2 / W3 / W4, per arm and degree.** This is `_assert_identities`'s set, applied to each arm the window solves:
  - (i) ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`).
  - (ii) σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` (`test_port_birdcage_four_port.py:129`), and every column power ≤ 1.
  - (PP) The `PORT-16` exact power partition, |p_src − p_vol − sheet| / |p_src| ≤ `DISCRETE_IDENTITY_RTOL` 1e-6 (`test_birdcage_power_identity.py:227`), per drive, via `_exact_shares`. It is backed at degree 2 by `test_birdcage_power_identity_degree2.py:168`.
  - (M) Rod arm: σ_y mirror pairs ≤ `ADJACENT_SPREAD_BAND` (`test_port_birdcage_lumped_column.py:121`).
  - (S) Off arm: separation ≥ `POOLED_SEPARATION_FLOOR` (`test_port_birdcage_four_port.py:173`).
  - (E) Rod arm: max |E| on the reduced tag-402 dofs ≤ `PEC_WALL_MAX_E` 1e-12 (`test_imp3_rod_birdcage.py:135`), on every drive, on the space of the degree being solved.
  - (V) `valid_mask` all True at the six `PTS`.
  - (D) The global dof count equals 2E + 2F at degree 2, and E at degree 1, where E and F are that mesh's global edge and face counts. This is an exact counting identity and it is the order knob: item 252 measured it exactly (753 062 / 1 090 614).

**Negative controls.**
- **(Kw) asserted.** The rod wall is left natural, P1, at degree 1, in W4 only. It runs on a fresh solver on `replace(problem, facet_tags=outer, pec_facet_tags=(OUTER_BOX_TAG,))`, using `_outer_box_tags` and `_control_c1`'s construction (`:639`).
  - Assertion: max |E| on the 402 dofs > `PEC_WALL_MAX_E`.
  - Backing: the same exact-rectangle fixture at the record resolution, item 235 W6, read 1.297e-4 against 1e-12 (`20261005T035322Z_IMP-3-item235-W6-ours-p1.log:3797`).
  - Ceiling ≈ 1e8×.
- **(D) asserted** is the order control. A run that silently solved degree 1 fails it by a factor of ≈ 5.4.
- **(C-p) predicted (rule (e)), printed:** |η_rod^o2(T1)/η_rod^o1(T1) − 1| on the 1 g rung ≥ 1e-3.
  - Ceiling ≈ 1, since a shift cannot exceed the field itself. The record mesh read 4.96e-2.
  - Below 1e-3 is reported, not a stop. With (D) green, that would be the p-stability answer, not a defect.
- **(C2) predicted, printed:** |η_rod(T1) − 1| / max over F1, F2 of |η_rod(F) − 1|, at both orders. Predicted ≥ 3. The record mesh at degree 2 read 2.594 / 0.171 ≈ 15 (log `:3910`).

**Printed, never gated.** P1 is the reading of record; P2–P4 and ccw / cw are printed beside it.
- η_rod(p) at both orders, at all six points.
- **η_rod^o2(T1)/η_rod^o1(T1) − 1 on the 1 g rung, beside −4.961 % on the record mesh.**
- The move from the record mesh to the 1 g rung at fixed order, η_rod^{1g}(T1)/η_rod^{rec}(T1) − 1, at degree 1 and at degree 2.
  - The record values come from item 242's `ours_imp3p_fields.json` and item 252's `ours_imp3q_fields.json`, both read-only.
  - The two meshes differ in resolution and are remeshed, so this is a reading, not a rate.
- Q (top and bottom tip) at degree 1 and 2, and Q^o2/Q^o1 − 1. Q^o1 is printed beside step 1's ragged 4.496 (P1 top).
- `pec_prolate_spheroid_field` (`utils/analytical.py:1576`) beside η_rod(p)/η_rod(F1) at T1–T3 and S1, as items 242 / 252 print it. No accuracy claim.
- Per window: cells, global dofs, sweep seconds and summed peak RSS. `_peak_rss_gib` is called by every rank, outside any `rank == 0` branch.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; container-side `timeout -k 30` ≤ 590 s, the carried §9 rule; ranks ≤ 12).
- **Pre-check** (not a harness call): `ls logs/imp3/`.
  - Any `ours_imp3g_*` present: stop and report. Never overwrite.
  - `ours_imp3q_fields.json` and `ours_imp3p_fields.json` must be present.
- **W1, census** (≈ 130 s, range 100–200).
  - Basis: the ragged builds' mesh wall was 35.4 + 58.3 s (`…170843Z_IMP-3-step1-W2.log:1918,3758`).
  - It runs serially first because a gmsh failure deadlocks a wide run (F10a), and this combination is new.

  `scripts/testing/run_and_log.sh IMP-3-item263-W1-build "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP6=build timeout -k 30 300 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ratio_1g.py -v -s --tb=short'"`
- **W2, degree 2 without the rod: the cost probe** (≈ 240 s, range 190–330; ≈ 25–35 GiB summed).
  - Basis: one four-drive degree-2 sweep on 259 524 cells took 138.01 s at `-n 8` (`20261004T183338Z_ANS-4-item225-W3-p2-cf0p5-128-koff.log:2083`; the window took 261 s).
  - At 64 MHz on 117 206 cells it took 64.14 s (item 192).
  - Scaling exponent 1.23 gives ≈ 125–170 s for ≈ 201 000 cells. Add the build (≈ 40 s), 12 SAR balls at quadrature 16, and the readings.

  `scripts/testing/run_and_log.sh IMP-3-item263-W2-d2off "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP6=d2off timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ratio_1g.py -v -s --tb=short'"`
- **Cost gate (binding, after W2).** W3 projects to ≈ 1.45 × W2's elapsed time.
  - Basis: the rod build took 61.8 s against 38.7 s on the ragged rung. The sweep scales by (252/201)^1.23 ≈ 1.32, plus the hole's preconditioner growth.
  - **W2 elapsed ≤ 380 s:** run W3.
  - **W2 elapsed > 380 s:** stop. Journal W2's sweep seconds, dofs and RSS as the measured price, park, and mark BLOCKED (priced over the ceiling). The review re-queues; `xl` is the weekly's to commission (I15 already names it).
  - **Never lengthen a timeout, and never raise `-n` in-slot.**
- **W3, degree 2 with the rod** (≈ 320 s, range 250–480; ≈ 35–45 GiB summed; item 212 W1 read 36.12 GiB at 259 524 cells).

  `scripts/testing/run_and_log.sh IMP-3-item263-W3-d2rod "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP6=d2rod timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ratio_1g.py -v -s --tb=short'"`
- **W4, degree 1 on both arms, plus (Kw) and the readout** (≈ 270 s, range 230–360).
  - Basis: step 1 W2 ran the same two builds and two degree-1 sweeps with SAR in 260 s at `-n 4` (`…170843Z_IMP-3-step1-W2.log:4162`). (Kw) adds one P1 solve.
  - This is the closing window: no module edit after it starts.

  `scripts/testing/run_and_log.sh IMP-3-item263-W4-d1-readout "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP6=d1 timeout -k 30 420 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ratio_1g.py -v -s --tb=short'"`
- **Compute** ≈ 130 + 240 + 320 + 270 ≈ 960 s ≈ 16 min, plus ≈ 2 min of container overhead. The tops of the ranges give 1 370 s ≈ 23 min.

**Traps already paid for.**
- **Degree defaults to 1 in two places:** `_arm` (via `build_four_port_sweep`) and item 235's `_sweep`. Pass `degree` explicitly; (D) is the tell.
- **Fixture build.**
  - Never build the rod fixture without `leg_sheet_strip_fraction`: that gives ragged sheets (binding ruling 2026-10-02).
  - Build without the rod first, then with the rod (`GEO-23`).
  - A gmsh FAIL at W1 is recorded verbatim and never retried.
- **Quadrature pinned for the degree.**
  - The energy and reading forms take `quadrature_degree=4` at degree 2 (`_readings2`).
  - A `SpatialCoordinate`-bearing integral without pinned quadrature can stall FFCx (POST-5).
  - No UFL ordering comparison, and no `sqrt` comparison, on complex operands. Take ratios in numpy.
- **Point evaluation:** use `owner="lowest_original_cell"` with `PTS` imported bitwise. Never call `f.eval(points, np.arange(n))`.
- **Rank-local quantities** (`assemble_scalar`, wall dofs, cell tags, RSS) must be reduced. Keep every collective outside any `rank == 0` branch: a rank-0-only allreduce deadlocked item 227's W1 at `-n 8`.
- **Tag 402 and the cell tags live in the build tuple** (`built[2]` / `built[1]`), not in the sweep's.
- **Memory:** evaluate before `del res`. Free the off arm's fields, case and mesh before the rod build (`_arm`'s step 7).
- **A held factorisation solves a different problem** (`PORT-19`): (Kw) gets its own solver on the replaced problem.
- **Divide window prices by the solves in them** (item 237).
- **Harness discipline.**
  - `-s`, `tests/environment` first, `FEM_EM_REQUIRE_COMPLEX=1`, foreground, no pipes, no `-k`.
  - Status 124: run `pgrep -c python3` in the container first.
  - A killed compile leaves an orphan `.c`: run the rubric's cache sweep, and only when nothing is running.
  - Rule (i): no edit after a file's last window.
- **Never write** `ours_imp3x_p1.json`, `ours_imp3p_fields.json`, `ours_imp3q_fields.json` or anything under `logs/ops61/`. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: ours at degree 1 and 2, 64 MHz, on the exact-rectangle rod fixture at the 1 g rung (σ = 800 S/m untuned coil, 50 Ω, κ off), on the off and rod-PEC arms. It measures the identities, η_rod at six points and Q at both tips. The order move is a **one-mesh p-reading**, and the move from the record mesh to the 1 g rung is a remesh reading. Neither is an h-statement or a rate. Not covered:
- **any tip-SAR accuracy claim.** The accuracy of Q / η_rod is `ANS-10`'s, which the 10-03 weekly re-homed to Palace. This item runs no Palace;
- the near-tip h-ladder (I15, priced by this item);
- the natural-wall arm at degree 2 (η_wall);
- absolute |E| / SAR, tissue, heating, compliance, ISO 10974 (`WF-9`), ASTM F2182 (7c);
- the tuned or copper coil, F-human, other rod placements; AED, `src/`, any band.

**Status it can move.**
- **§7 `IMP-3` stays ✅** (no glyph). Rule (j) read: the Done-when, "the same identities green on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture", is already met by step 1.
  - The review appends a **step-6 record sentence**. Candidate text, valid only if W1–W4 are green: "On the exact-rectangle 1 g rung (⟨n_off⟩ / ⟨n_rod⟩ cells) the identities hold at degree 1 and 2; η_rod(T1) moves ⟨y⟩ % from order 1 to 2 there (record mesh −4.961 %) and Q (P1, top) reads ⟨q1⟩ / ⟨q2⟩. One-mesh p-reading, no accuracy claim."
- **§10 I14 is done.** Its measured W2 / W3 prices are I15's price; the weekly decides between ordinary heavy and `xl`.
- **`ANS-10`** (no §7 row):
  - if |p-move| < 1 % on the 1 g rung, the weekly may read the tip ratio as p-stable at the rung;
  - if not, I15 is needed before any accuracy reading.
  - This item sets no threshold.

**Negative result.**
- **W1 gmsh FAIL, or (G1) / (G2) / (G3) red:** this is not the fixture. Record it verbatim, open a known-issues entry, and mark BLOCKED on the mesher. Never retry.
- **Cost gate trips, or W2 / W3 returns Status 124:** run `pgrep -c python3` first. The result is a partial record (W1 plus W2's price). Park on `attempt/IMP-3-step6-<ts>` and mark BLOCKED (over the ceiling). Never lengthen the timeout or widen `-n`.
- **(i), (ii), (PP), (M), (S), (E), (V), (D) or (Kw) red:** report every number, park on `attempt/IMP-3-step6-<ts>`, open a known-issues entry, BLOCKED. Never widen. (PP) red at degree 2 only would be a finding about `_exact_shares` on a hole mesh at degree 2.
- **A printed reading outside its prediction** (including (C-p) < 1e-3 or (C2) < 3): report it. No status effect.
- There is no `PREDICTED RED` line.

**Outcome (2026-10-06 06:00 CDT slot, second item, implementer): BLOCKED — the binding cost gate tripped after W2 (397 s > 380 s). No red anchor.** Module parked on `attempt/IMP-3-step6-20261006T1114Z` (`c7de9eae`).
- **W1** `20261006T111102Z_IMP-3-item263-W1-build.log` (178 s, `-n 1`, 13 passed): first build of the exact-rectangle 1 g rung with the rod, no gmsh failure. (G0) off 201 072 cells (+0.576 % on the ragged 199 920; printed), mesh wall 39.05 s; rod 252 118 cells, mesh wall 58.77 s. (G1) 252 118 < 320 000; (G2) all eight N_i areas 5.600000000000e-05; (G3) rod 4 958 tag-402 facets, all exterior; off 0; (G4) a_1g 6.2035 mm, 0.021204 < 0.03 and 0.022204 < 0.04.
- **W2** `20261006T111406Z_IMP-3-item263-W2-d2off.log` (397 s, `-n 8`, 13 passed): every anchor green on the degree-2 off arm — (i) 2.561e-15; (ii) σ_max 0.999823786, column powers ≤ 0.834906; (PP) ≤ 8.389e-15 per drive; (S) 1 270.3× (floor 10); (V) True at all six points; **(D) 1 284 856 = 2E + 2F exact** (E 238 627, F 403 801). Ball masses 1.000139 / 0.999488 g. SAR_1g top P1 1.183938e-05 W/kg (off arm, degree 2).
- **Price (the measured W2 price, I15's input):** sweep 246.15 s (priced 125–170 s: 1.45–1.97× over), 12 SAR balls 26.26 s, arm 281.6 s, build 43.59 s, test body 325.3 s, pytest 395.3 s; summed peak RSS 44.17 GiB at `-n 8` (priced 25–35); 1 284 856 dofs (predicted ≈ 1.3 M). The gate projects W3 to ≈ 1.45 × 397 ≈ 575 s, at the 590 s container ceiling.
- **Not run:** W3 (degree-2 rod), W4 (degree 1, (Kw), readout). No η_rod, p-move, Q or (C-p) / (C2) reading exists yet.
- Artifacts (gitignored): `logs/imp3/ours_imp3g_build.json`, `logs/imp3/ours_imp3g_d2off.json` — a re-run must move or reuse them (the module opens with mode `"x"`).
- **Hypothesis / next:** the degree-2 sweep on a 0.0025 m phantom is costlier than the 259 524-cell 128 MHz scaling implied (the fine phantom adds dofs where the preconditioner works hardest). Re-pose options for the review: (a) the d2off result stands — re-run only W3 and W4, with W3's sweep and SAR in separate windows (save fields, or drop the ccw / cw SAR balls); (b) `xl` (I15) runs W3.

**Review disposition of item 263 (2026-10-06 10:30 daily review): re-posed as §9 item 284 (attempt 2 of 4, rule (6) "split windows"). W2 stands, and only W3 and W4 run, on the module exactly as parked.**
- **The gate was a crude projection.** It priced W3 at 1.45 × W2's *whole window* (≈ 575 s). Re-pricing from W2's own component clocks, on the same fixture and at the same width, gives ≈ 510 s:
  - fixed overhead outside the test body (environment tests, collection, MPI start) ≈ 70 s (pytest 395.3 s − test body 325.3 s);
  - rod build ≈ 64 s (W2's off build 43.59 s, plus W1's rod-minus-off mesh wall, 58.77 − 39.05 s);
  - rod sweep ≈ 325–335 s: 246.15 s × (252 118 / 201 072)^e. The brief's e = 1.23 gives ×1.32. Palace's own hole growth on the record mesh, 429.9 / 261.1 s for 168 507 / 117 206 cells (item 252 W2 / W3), gives e ≈ 1.37 and ×1.36;
  - SAR balls ≈ 29 s and readings ≈ 10 s.
- **Range [470, 560] s.** That is inside 590 with a 30–120 s margin.
- **Neither of the slot's sub-options is needed.** The SAR balls are 26 s of W2, so splitting them out buys little. Saving degree-2 fields is new code, and new code would void W2 under rule (i).
- **`-n 12` is not a measured lever.** On one 259 524-cell degree-2 mesh, four driven solves took 183.12 s at `-n 12` (`20261004T094455Z_ANS-4-item212-W2-p2-cf0p5-64.log:2170`) against 138.01 s at `-n 8` (item 225 W3).
- **If W3 is killed at 590 s, that is I15's price record** ("rod arm at degree 2 on the 1 g rung > 590 s at `-n 8`"), and the weekly commissions `xl`.

### §9 item 284 — 2026-10-06 (chain I14, re-pose, attempt 2 of 4)

*Written by the 2026-10-06 10:30 daily review.* Licence: §9 rule (6), re-posing the blocked step I14 by splitting windows. Item 263's W1 / W2 are kept as run, and its W3 / W4 run in this item, unchanged. The target and every anchor, band, control, prediction, scope line and negative-result line of item 263's brief carry verbatim. Only the cost gate is replaced by the re-price above. Family `IMP-3` step 6, attempt 2 of 4. **No band, tolerance, Done-when or target is set or moved.**

- **Executor:** implementer. **Tier:** heavy; W3 `-n 8`, W4 `-n 4`. **Class:** chain step (solve), ≈ 19 measured / ≈ 30 predicted min (≈ 8 fixed + W3 ≈ 8.5 + W4 ≈ 4.5 + ≈ 9 spare).
- **Depends on:** nothing. Never concurrent with another heavy window.
- **Land** `tests/validation/test_imp3_tip_ratio_1g.py` by path from `attempt/IMP-3-step6-20261006T1114Z` (`c7de9eae`), **unedited**: it is the module W1 / W2 ran, so W2 stands under rule (i). Check that `git diff c7de9eae -- tests/validation/test_imp3_tip_ratio_1g.py` is empty before W3. W1 / W2's logs and rows are already on `main`.
- **Pre-check** (not a harness call): `ls logs/imp3/`. `ours_imp3g_build.json` and `ours_imp3g_d2off.json` must be present, and `ours_imp3g_d2rod.json` and `ours_imp3g_d1.json` absent (the module writes with mode `"x"`). `ours_imp3q_fields.json` and `ours_imp3p_fields.json` must be present. Any mismatch: stop and report, never delete.
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0, then run the rubric's cache sweep (label `IMP-3-item284-cache-sweep`).
- **W3** is item 263's W3 command verbatim, label `IMP-3-item284-W3-d2rod`, predicted ≈ 510 s, range [470, 560], ≈ 55–65 GiB summed. Bash `timeout: 660000`.
- **W4** is item 263's W4 command verbatim, label `IMP-3-item284-W4-d1-readout`, predicted ≈ 270 s, range 230–360. This is the closing window, and it reads both degree-2 files.
- **Fallback:** if W3 returns Status 124, run `pgrep -c python3` first, then the cache sweep. The killed window **is** I15's price record. Do not run W4. Park, mark BLOCKED (over the ceiling), and journal the last printed phase. **Never lengthen the timeout, never raise `-n`, never split in-slot.**
- **Traps:** item 263's list. Also, the module's `"x"` writes mean a second W3 attempt in the same slot would fail on its own d2rod file. There is no second attempt.

**Status it can move:** item 263's, unchanged. `IMP-3` stays ✅, and the review appends the step-6 record sentence from item 263's candidate text. §10 I14 is done, and its W2 / W3 prices are I15's price. `ANS-10`'s p-stability reading goes to the weekly.

**Outcome (2026-10-06 12:00 CDT slot, third item, implementer): GREEN — every anchor green on W3 and W4; no red, no control outside its prediction.** Module landed unedited from `c7de9eae` (empty diff checked before W3). Pre-check matched; `pgrep -c python3` 0; cache sweep `20261006T172410Z_IMP-3-item284-cache-sweep.log` (85 s, no orphan).
- **W3** `20261006T172542Z_IMP-3-item284-W3-d2rod.log` (450 s, `-n 8`, 13 passed; **below** the [470, 560] range): rod arm at degree 2, 252 118 cells — (i) 3.937e-15; (ii) σ_max 0.999828021, column powers ≤ 0.834970; (PP) ≤ 5.0e-15 per drive; (M) ≤ 0.0275 % (band 0.50 %); (E) 24 790 tag-402 dofs, max |E| 0 on all four drives; (V) True; **(D) 1 620 748 = 2E + 2F exact** (E 302 002, F 508 372). Price: sweep 284.88 s, 12 SAR balls 29.91 s, arm 319.2 s, body 381.0 s; summed peak RSS 50.97 GiB at `-n 8` (priced 55–65).
- **W4** `20261006T173320Z_IMP-3-item284-W4-d1-readout.log` (322 s, `-n 4`, 13 passed, in range): degree 1, off / rod — (i) 6.235e-16 / 7.147e-16; (ii) σ_max 0.999789653 / 0.999796114; (PP) ≤ 7.7e-15 / ≤ 3.7e-15; (S) off 604.9× (floor 10); (M) rod ≤ 0.0312 %; (E) rod 7 437 dofs, max |E| 0; (V) True on both; (D) 238 627 = E / 302 002 = E exact. **(Kw)** natural wall P1: max |E| on 402 dofs 1.263e-04 > 1e-12 (1.263e8×), (PP) 4.0e-15. Summed peak RSS 5.75 GiB; body 251.4 s.
- **Readings (printed, never gated; P1 of record):** η_rod(T1) o1 3.704420 / o2 3.482328; **p-move −5.995 %** on the 1 g rung beside −4.961 % on the record mesh. (C-p) 5.995e-2 ≥ 1e-3: IN. (C2) o1 10.243, o2 16.321 (≥ 3: IN). Remesh move record → 1 g rung, η_rod(T1): o1 −2.056 %, o2 −3.122 %. Q P1 top o1 4.476027 / o2 4.463326 (−0.284 %), bottom 4.469228 / 4.462910 (−0.141 %); Q^o1 top beside step 1's ragged 4.496087. Q P2 / P4 move ≈ −10 % between orders; ccw / cw ≈ −0.3 to −0.5 %. Closed form / η ratio at T1–T3, S1 (o2): 0.843 / 0.881 / 1.050 / 1.073 — no accuracy claim. Rod-arm 1 g ball masses 0.982 g (off 1.000 g).
- **For the weekly (`ANS-10`, no threshold set here):** |p-move| of η_rod(T1) is 6.0 % (> 1 %), so by item 263's reading η_rod is not p-stable at the rung; the 1 g Q (P1) moves 0.28 % between orders. One-mesh p-reading, no accuracy claim. I15's price: degree-2 sweeps 246.15 s (off, 1 284 856 dofs, 44.17 GiB) and 284.88 s (rod, 1 620 748 dofs, 50.97 GiB) at `-n 8`.

### §9 item 305 — 2026-10-07 (chain I15 part 2 — the r/(4√2) rung at degree 2; conditional on the 10-07 weekly)

*Drafted by a general-purpose agent for the 2026-10-07 18:00 daily review; the review reviews it.* Licence: §10 I15 (`PROJECT_PLAN.md:4672–4679`), a weekly-enumerated step (§9 rule (4)), serial on I14 (item 284, landed, `IMP-3.md:947–951`) and on I15 part 1 (`GEO-36` step 4, item 298, `a49f997c`, `GEO-36.md:211–236`). I15 leaves the solve tier to the weekly ("`xl` if I14 prices it above the ordinary ceiling (the 10-07 or 10-10 weekly commissions)"; `IMP-3.md:899`). **No band, tolerance, Done-when or target is set or moved.**

**CONDITIONAL — read this first.** The item runs only if **all three** hold on `main` when the slot starts:
1. a commit `docs(plan): weekly review 2026-10-07` exists, dated after 21:00 CDT on 2026-10-07 (`git log --oneline main --grep "docs(plan): weekly review 2026-10-07"`);
2. in that commit, §10 I15's bullet rules **the r/(4√2) rod-arm degree-2 window heavy** (the ordinary tier, `-n 8`, container `timeout -k 30 590`), **not `xl`**;
3. the same bullet **keeps I15's rung set** (r/4 → r/(4√2) → r/8 on `implant_rod["tip_size"]`) **and its question** (the η_rod(T1) sequence at degree 2 on the 1 g rung). A weekly that adds a ball-radius knob, swaps a rung, or re-words the question has changed it.

If the commit is absent, rules `xl`, changes the rung set or question, or reads ambiguously on any of the three, **the slot marks this item BLOCKED ("waits on the weekly's I15 ruling") with nothing run** — no build, no pre-check window, no cache sweep — and takes the next item. The executor never interprets the weekly's text in its own favour.

**`IMP-3` step 7, part 2 — the r/(4√2) rung's rod arm at degree 2.** I14's exact-rectangle hole-rod fixture on the 1 g rung (`phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG`, `leg_sheet_strip_fraction=GATED_WIDTH_FRACTION`, 64 MHz, κ off, P1 of record), built with `implant_rod=dict(ROD, tip_size=TIP_R4S2)` (ROD_RADIUS / (4√2) = 1.7678e-4 m). Only the **rod arm** is solved. The off arm carries no rod, so `tip_size` never reaches it (GEO-36 W2 printed the off build at 201 072 cells with the knob in the code, `GEO-36.md:226`). Its degree-2 record is I14 W2 (`20261006T111406Z_IMP-3-item263-W2-d2off.log:2084`, every anchor green, `IMP-3.md:914`), read back from `logs/imp3/ours_imp3g_d2off.json`, never re-solved. **Target: the identities on the refined mesh. The first step of the η_rod(T1) ladder is a printed reading.** The step ratio and the convergence statement need r/8, which is `xl` (`GEO-36.md:234`), and neither is claimed here.

- **Executor:** implementer. One new test module, `tests/validation/test_imp3_tip_ladder_1g.py`. No `src/`, and **no edit to any existing file**. It imports I14's and GEO-36 step 4's helpers by name, so rules (a) and (c) are not triggered and no I14 re-gate (a heavy window) is owed.
- **Tier:** heavy. W0 `-n 1`, W1 `-n 8` (I14's W3 width). Both windows run in the complex container.
- **Class:** chain step (solve).
- **Cost:** ≈ 10 measured / ≈ 45 predicted min (15 fixed + 12 code + ≈ 11 compute + ≈ 7 spare).
- **Depends on:** the 10-07 weekly ruling above, and nothing else. Every module named below is on `main`. Never run concurrently with another heavy window.

**Settled from the code.**
- **Imported from `tests/validation/test_imp3_tip_ratio_1g.py`** (I14's module, last committed `0f6bbbc7`; helpers only, never `test_imp3_tip_ratio_1g`):
  - `_arm_d` `:150`, which builds the four-port case at the given degree, asserts (G1) before the solve, sweeps, computes the 12 SAR balls and evaluates E at `PTS`;
  - `_identities` `:241`, which runs (i), (ii), (PP), (M), (E), (V), (D) and (G2) for one arm;
  - `_ser` `:283`, `_cplx` `:302`, `_eval_points` `:138` and `_edge_face_counts` `:131`.
  - **Not** `_built` `:121`: it takes no `tip_size` and silently builds the r/4 rung (control (Kc) is the tell). **Not** `_save` `:295`: its prefix is fixed to `ours_imp3g_`.
- **Imported from `tests/mesh/test_birdcage_implant_rod_tip.py`** (GEO-36 step 4, `a49f997c`): `TIP_R4S2` `:113`, `FIXTURE_KW` `:131`, `ROD_1G_RUNG_CELL_RECORD` `:90`, `T1` `:84`, `T1_RADIUS` `:100`, `_ball_and_t1` `:158`, `_cell_geometry` `:142` and `_mean_edge` `:150`. Never the `builds` fixture and never a `test_*` name.
- **Also imported:** `ROD` and `_build_rod` (`tests/mesh/test_birdcage_implant_rod.py:125`); `_outer_box_tags`, `_peak_rss_gib`, `TIPS`, `DRIVES`, `PORT_IDS` (`tests/validation/test_imp3_rod_birdcage.py`); `PNAMES`, `PTS` (`test_imp3_palace_fields.py`); `point_owner_cells` (`fem_em_solver.post.evaluation:133`).
- **The build** is `built = _build_rod(implant_rod=dict(ROD, tip_size=TIP_R4S2), **FIXTURE_KW)`, the call GEO-36 W2 made (`test_birdcage_implant_rod_tip.py:257–258`). Then `a = _arm_d(built, "degree 2 rod r/(4√2)", True, 2)`.
- **The writer** is local: rank 0 writes `logs/imp3/ours_imp3h_r4s2.json` (mode `"x"`) with `_ser(a)` plus `rss_gib`, `ranks`, `tip_size_m`, the near-T1 and T1-cell edges, and a barrier outside the rank branch.
- **Env knob:** `FEM_EM_IMP3_STEP7` ∈ {`r4s2`}. One test function dispatches on it and skips otherwise. Never `-k`.
- **Record constants** the module registers, each citing its log line:
  - `I14_ETA_O2_T1_P1 = "3.482328"` (`20261006T173320Z_IMP-3-item284-W4-d1-readout.log:3894`);
  - `I14_Q_O2_P1_TOP = 4.463326` (`…:3905`);
  - `GEO36_R4S2_CELLS = 276_088`, `GEO36_R4S2_D2_DOFS = 1_775_954`, `GEO36_R4S2_N402 = 6_256`, `GEO36_R4S2_NEAR_T1 = (19, 6.422237e-04)` and `GEO36_UNSET_NEAR_T1_EDGE = 6.983224e-04` (`20261007T170700Z_GEO-36-step4-W2-rungs.log:5612–5613, 7516`). These are printed comparands, not bands.
- **η_rod is item 252's definition**, as I14 used it (`IMP-3.md:784`): η_rod^h(p) = ‖E_rod^h(p)‖ / ‖E_off(p)‖ per drive, with the Euclidean norm over components.
  - E_off(p) comes from `ours_imp3g_d2off.json` (I14 W2) via `_cplx`, read-only.
  - η_rod^{r/4} comes from `ours_imp3g_d2rod.json` (I14 W3), read-only.
  - E_rod^{r/(4√2)} is this window's `a["E"]`.
- **Geometry already known** (`GEO-36.md:175, 209 (4), 227`). T1 sits 2.531 mm from the cap centre, outside the 2 mm cap ball, so the knob refines T1's neighbourhood only by grading spill-over. GEO-36 W2 measured the near-T1 (0.5 mm) mean edge at 6.983224e-04 → 6.422237e-04 m (−8.0 %) at r/(4√2), and −27 % at r/8 (log `:7516`). **The local step at T1 is therefore ≈ 1.087, not the requested √2.** This item prints the edge; it names no step ratio.

**Anchors (asserted; bands imported by name, never restated).** On the rod arm at r/(4√2), degree 2, through the imported helpers:
- (G1) cells < `ROD_CELL_GUARD_W2` (`test_imp3_rod_birdcage.py:115`), inside `_arm_d` before the solve. The predicted 276 088 passes. r/8's 342 072 would trip it, but r/8 is not run here.
- (G2) every N_i area = 5.6e-5 within `STRIP_BAND` (`test_ops61_exact_strip_ports.py:58`), inside `_identities`.
- (G3) tag-402 facets > 0 and all exterior, via `_outer_box_tags(msh, built[2], comm)`, as I14's `_census` does it (`test_imp3_tip_ratio_1g.py:317, 342–343`).
- (i) ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`).
- (ii) σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` (`test_port_birdcage_four_port.py:129`), and every column power ≤ 1.
- (PP) the `PORT-16` exact power partition ≤ `DISCRETE_IDENTITY_RTOL` (`test_birdcage_power_identity.py:227`) per drive, via `_exact_shares`. It is backed at degree 2 by `test_birdcage_power_identity_degree2.py:168`.
- (M) σ_y mirror pairs ≤ `ADJACENT_SPREAD_BAND` (`test_port_birdcage_lumped_column.py:121`).
- (E) max |E| on the reduced tag-402 dofs ≤ `PEC_WALL_MAX_E` (`test_imp3_rod_birdcage.py:135`), on every drive.
- (V) `valid_mask` is all True at the six `PTS`.
- (D) global dofs = 2E + 2F exactly (the counting identity I14 W3 and GEO-36 (D) measured).
- (S) is off-arm only and is not re-run. Its record is I14 W2's 1 270.3× (`IMP-3.md:914`).

**Negative controls (rule (e) labels).**
- **(Kc) asserted — the knob reached gmsh.** Two checks:
  - `built[3]["implant_rod_cad"]["tip_size_m"] == TIP_R4S2` bitwise;
  - rod cells > `ROD_1G_RUNG_CELL_RECORD`.
  - Backing, same comparison on the same build: GEO-36 W2 read `tip_size_m` equal to the request bitwise (log `:7503`) and 276 088 / 252 118 = 1.0951 (`:7516`).
  - A build that silently took I14's `_built("rod")`, or a key that never reached gmsh, reads exactly 252 118 and fails by an integer. There is no ceiling: the separation is the integer itself.
- **(D) asserted** is the order control. A run that silently solved degree 1 fails it by a factor of ≈ 5.4 (331 014 against 1 775 954, GEO-36 W2 `:5612`).
- **(Rb) asserted — the off-arm reuse is wired.** η_rod^{r/4}(T1), P1, is recomputed in this window from the two I14 JSONs. It must be **string-equal at 6 d.p.** to `I14_ETA_O2_T1_P1` "3.482328", the same comparison on the same files (W4 log `:3894`). The JSONs' `cells` / `ndofs` must equal 201 072 / 1 284 856 and 252 118 / 1 620 748 as integers. A red means the reading path is wrong, and no η at r/(4√2) may be printed as a reading.
- **(C-h) predicted, printed:** |η_rod^{r/(4√2)}(T1)/η_rod^{r/4}(T1) − 1| (P1) ≥ 1e-3.
  - Ceiling ≈ 1, since the move cannot exceed the field.
  - Basis: the near-T1 edge moved −8 %, and I14's remesh move record → 1 g rung was −3.122 % at degree 2 (`IMP-3.md:950`).
  - Below 1e-3 is reported, not a stop. With (Kc) and (D) green, that would be the first rung's answer, not a defect.
- **(C2) predicted, printed:** |η_rod(T1) − 1| / max over F1, F2 of |η_rod(F) − 1| at r/(4√2) ≥ 3 (r/4 read 16.321, `IMP-3.md:950`).
- **Mesh identity predicted, printed:**
  - cells, tag-402 facets and degree-2 dofs equal `GEO36_R4S2_*` as integers;
  - near-T1 count and mean edge equal `GEO36_R4S2_NEAR_T1` to the printed digits.
  - Basis: gmsh builds serially on rank 0 (`io/mesh.py:3998,4480`), and 252 118 read the same at `-n 1` and `-n 8` (`GEO-36.md:145`).
  - A mismatch is reported, not a stop. (Kc) is the asserted knob check.

**Printed, never gated.** P1 is the reading of record; P2–P4 and ccw / cw are printed beside it.
- **η_rod^{r/(4√2)}(p) at all six points, beside η_rod^{r/4}(p)** (recomputed under (Rb)). Per point, the relative move η^{r/(4√2)}/η^{r/4} − 1.
- **The first ladder step at T1, P1:** η_rod(T1) at r/4 (3.482328) → r/(4√2), and its relative move. It is printed beside I14's p-move −5.995 % and the remesh move −3.122 % for scale. **The step ratio is not computed**: it needs r/8 (`xl`).
- The F1 / F2 moves, printed as this remesh's floor beside T1's move.
- **T1's local mesh:**
  - the near-T1 (`T1_RADIUS`) cell count and mean edge via `_ball_and_t1(msh, built[1], cad["tip_ball_radius_m"], comm)`, beside GEO-36's r/4 6.983224e-04 and r/(4√2) 6.422237e-04 m, with the local edge ratio h(r/4)/h(r/(4√2)) beside the requested √2;
  - **the T1 cell's own mean edge**: `point_owner_cells(msh, T1[None, :], comm)`, then `_mean_edge(_cell_geometry(msh)[[cell]])` on the owning rank, made global by an allreduce (MAX, 0 elsewhere) outside any rank branch, printed with its original cell index. This is a first record at this rung: no r/4 value exists.
- N_ball and the mean ball edge from the same call, beside GEO-36 W2's 32 826 / 2.433e-4 m.
- Q = SAR_1g,tip(rod) / SAR_1g,tip(off) per tip and drive, with SAR_off from `ours_imp3g_d2off.json`. It is printed beside the r/4 value recomputed from `ours_imp3g_d2rod.json` (P1 top 4.463326, log `:3905`), together with the rod-arm 1 g ball masses (r/4: 0.982 g).
- `pec_prolate_spheroid_field` (`utils/analytical.py:1576`) beside η_rod(p)/η_rod(F1) at T1–T3 and S1, as I14 prints it. No accuracy claim.
- Cells, E, F, global dofs, build seconds, sweep seconds, SAR seconds and summed peak RSS (`_peak_rss_gib`, called on every rank) against GEO-36's price row (`GEO-36.md:233`: sweep 318.8–322.9 s, window 491.3–495.4 s, 57.04–57.77 GiB).

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; container-side `timeout -k 30` ≤ 590 s; `-s`; `tests/environment` first; `FEM_EM_REQUIRE_COMPLEX=1`; no `-k`; no pipes). **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**
- **Gate (not a harness call):** the three conditions at the top of this item. If any fails, the item is BLOCKED and nothing below runs.
- **Pre-check (not a harness call):**
  - `ls logs/imp3/`: `ours_imp3g_d2off.json` and `ours_imp3g_d2rod.json` must be present, and `ours_imp3h_r4s2.json` absent (mode `"x"`).
  - `git diff 0f6bbbc7 -- tests/validation/test_imp3_tip_ratio_1g.py` and `git diff a49f997c -- tests/mesh/test_birdcage_implant_rod_tip.py` must both be empty: the imported helpers are the ones that ran.
  - Any mismatch: stop and report. Never delete or move a file.
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Then run the rubric's cache sweep, label `IMP-3-item305-cache-sweep` (≈ 85 s; item 284's took 85 s).
- **W0, import smoke** (`-n 1`, knob unset, so the test skips; ≈ 40 s, range 20–90). It catches an import or collection error before an eight-minute window.

  `scripts/testing/run_and_log.sh IMP-3-item305-W0-import "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short'"`
- **W1, the rod arm at r/(4√2), degree 2, plus the readout** (≈ 495 s, range [450, 570]; ≈ 57–58 GiB summed at `-n 8`).
  - Basis: GEO-36's price row, 491.3 / 495.4 s at e = 1.23 / 1.37 (`GEO-36.md:233`), scaled from I14 W3's 450 s window (sweep 284.88 s at 1 620 748 dofs, 50.97 GiB). The readout adds ≈ 5 s of JSON reads and norms. **This is a prediction, not a measurement on the refined mesh.**
  - This is the closing window (rule (i)): no module edit after it starts.
  - Durable capture is §5.1's idiom verbatim, trailing `; exit $rc` included:

  `scripts/testing/run_and_log.sh IMP-3-item305-W1-d2rod-r4s2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP7=r4s2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short > /workspace/logs/imp3-item305-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp3-item305-w1-raw.log; cat /workspace/logs/imp3-item305-w1-raw.log; exit \$rc'"`
- **Compute** ≈ 85 + 40 + 495 ≈ 620 s ≈ 10.3 min, plus ≈ 1 min of container overhead. The tops of the ranges give 745 s ≈ 12.4 min.
- **Fallback.** If W1 returns Status 124:
  - run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` first, then the cache sweep (the daily-review.md rubric's trap list), and only when it reads 0;
  - if the wrapper itself died, recover with `scripts/testing/run_and_log.sh --capture-orphan IMP-3-item305-W1-d2rod-r4s2 logs/imp3-item305-w1-raw.log` once the container command has finished;
  - **the killed window is I15's price record**: "r/(4√2) rod arm at degree 2 on the 1 g rung > 590 s at `-n 8`". Journal the last printed phase (build / sweep / SAR / readout) with its seconds;
  - park on `attempt/IMP-3-step7-<ts>` and mark BLOCKED (over the ceiling; the weekly commissions `xl` for both rungs);
  - **never lengthen the timeout, never raise `-n`, never split in-slot.**

**Traps already paid for.**
- **Item 263's list carries in full** (`IMP-3.md:866–887`):
  - degree defaults to 1 in two places;
  - never build without `leg_sheet_strip_fraction` (binding ruling 2026-10-02);
  - quadrature 4 at degree 2;
  - no UFL ordering or `sqrt` comparison on complex operands;
  - `owner="lowest_original_cell"` with `PTS` imported bitwise;
  - rank-local quantities reduced, with no collective inside a `rank == 0` branch (item 227's `-n 8` deadlock);
  - tag 402 and the cell tags live in `built[2]` / `built[1]`;
  - evaluate before freeing fields;
  - divide window prices by the solves in them;
  - Status 124 means `pgrep` first;
  - an orphan `.c` means the cache sweep, only with nothing running.
- **New here:**
  - `_built` and `_save` are not reusable (above). Calling `_built("rod")` gives the r/4 mesh and a run that looks green; (Kc) is the tell.
  - `_ball_and_t1` takes the **cell tags** `built[1]` and reduces internally. Call it on every rank.
  - The T1 owner lookup is collective (`point_owner_cells` allreduces). Call it on every rank, outside any rank branch.
  - Importing `test_birdcage_implant_rod_tip` runs its module-level probe load (`:83`), which is harmless. Never import its `builds` fixture: pytest would collect it into the new module.
  - **A requested `VIn` is not the step at T1**: print the measured local edge and never write "step ratio √2" for T1.
  - The module writes with mode `"x"`, so a second W1 in the same slot fails on its own file. There is no second attempt.
  - If the permission layer denies the `$rc` capture form (it did once, `PORT-21` step 5 W0b, `docs/planning/chunks/PORT-21.md:38`), stop and report. Do not improvise a different capture.
- **Never write** `ours_imp3g_*`, `ours_imp3p_fields.json`, `ours_imp3q_fields.json` or anything under `logs/ops61/`. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: ours at degree 2, 64 MHz, on the exact-rectangle hole-rod fixture at the 1 g rung (σ = 800 S/m untuned coil, 50 Ω, κ off), with the cap-ball `VIn` = r/(4√2), **rod arm only**. The off arm is I14's record, reused. It measures the identities on the refined mesh and prints the first step of the near-tip ladder. Not covered:
- **any convergence claim.** One step gives no rate, and the step ratio and observed order need r/8 (`xl`, the weekly's). Because T1 sits outside the cap ball, even the step's size in local h is a measured ≈ 1.087, not √2;
- **any accuracy claim** on η_rod or Q (`ANS-10`'s);
- the r/8 rung; the degree-1 arm at r/(4√2); the natural-wall control at this rung (I14's (Kw) at r/4 stands, `IMP-3.md:949`, and is not re-run);
- a ball-radius knob (`GEO-36.md:209` (4), the weekly's question);
- absolute |E| / SAR, tissue, heating, compliance, ISO 10974 (`WF-9`), ASTM F2182 (7c); F-human, the tuned or copper coil; AED, `src/`, any band.

**Status it can move.**
- **§7 `IMP-3` stays ✅** (no glyph); rule (j) read as in item 263 (`IMP-3.md:897`). The review appends a **step-7 record sentence**. Candidate text, valid only if W1 is green: "Step 7 part 2 (§9 item 305): with the cap-ball tip size at r/(4√2) on the exact-rectangle 1 g rung (⟨n⟩ cells, ⟨d⟩ degree-2 dofs) the rod-arm identities hold at degree 2; η_rod(T1) (P1) reads ⟨η⟩ against 3.482328 at r/4 (⟨Δ⟩ %), with the near-T1 edge ⟨e⟩ % shorter. One ladder step, no rate, no accuracy claim."
- **§10 I15 part 2 is done; I15 stays open** until r/8 lands. The review annotates I15's bullet. The slot makes no §10 edit. W1's measured window, sweep and RSS re-price r/8 for the weekly's `xl` commissioning, replacing the scaled prediction.
- **`ANS-10`:** no reading moves. The convergence statement waits for r/8.

**Negative result.**
- **Gate fails:** BLOCKED ("waits on the weekly's I15 ruling"), nothing run. This is not a negative result about the fixture.
- **Pre-check mismatch, or `pgrep` non-zero:** stop and report. Mark BLOCKED with the unblock condition (rule (d)). Never delete.
- **W0 red** (import or collection): this is a module defect. Fix it in-slot only if W1 has not started, then re-run W0. Otherwise park.
- **gmsh FAIL at W1:** record it verbatim with its surface pairs, open a known-issues entry, and mark BLOCKED on the mesher. Never retry. GEO-36 W2 built this rung cleanly, so a FAIL here would be a rank-count finding.
- **(Kc) red:** the knob did not reach the build. No η reading. Park, mark BLOCKED, and open a `GEO-36` known-issues entry.
- **(Rb) red:** the reading path is wrong. Report the identities as measured but print no η / Q reading. Park and mark BLOCKED.
- **(G1)–(G3), (i), (ii), (PP), (M), (E), (V) or (D) red:** report every number, park on `attempt/IMP-3-step7-<ts>`, open a known-issues entry, and mark BLOCKED. Never widen. (PP) red only here would be a finding about `_exact_shares` on the refined hole mesh.
- **Status 124:** the fallback above. The killed window is the price record.
- **A printed reading outside its prediction** ((C-h) < 1e-3, (C2) < 3, a mesh-identity mismatch, W1 outside [450, 570] s): report it. No status effect.
- There is no `PREDICTED RED` line.

**Review rulings (2026-10-07 18:00 daily review, item 305).** (1) **"Weekly not yet committed" is a dependency that has not landed, not a block:** a slot that finds no `docs(plan): weekly review 2026-10-07` commit on `main` skips item 305 and leaves it ⬜ (implementer-run.md step 2: never an item whose stated dependency has not landed). Only a weekly commit that rules `xl`, changes the rung set or the question, or reads ambiguously makes the slot mark it BLOCKED. This replaces the brief's "if the commit is absent … BLOCKED" for the absent case only. (2) **The unequal local step at T1** (edge ratio ≈ 1.087 then ≈ 1.264, because T1 sits outside the 2 mm cap ball) goes to the weekly with I15's ball-radius question; the item claims no step ratio and no observed order, as drafted. (3) **Capture:** if the durable-capture form is denied, stop and report, as drafted; do not improvise a bare window. (4) **W1's predicted range [450, 570] s** is a scaled prediction with 20 s of margin; a kill at 590 s is the price record, never a reason to lengthen, raise `-n` or split. (5) The new module `tests/validation/test_imp3_tip_ladder_1g.py` imports I14's and part 1's helpers and edits neither, so no re-gate window is owed; the brief's reasons for not reusing `_built` (no `tip_size`) and `_save` (fixed prefix) are accepted.

**Outcome (2026-10-07 22:30 CDT slot, implementer): GREEN — every anchor and both asserted controls green on W1; no reading outside its prediction.** Gate checked by the slot (weekly `732f960b`). Pre-check matched (both I14 JSONs present, `ours_imp3h_r4s2.json` absent, both helper diffs empty); `pgrep -c python3` 0; cache sweep `20261008T033202Z_IMP-3-item305-cache-sweep.log` (93 s, no orphan); W0 `20261008T033340Z_IMP-3-item305-W0-import.log` (71 s, 12 passed, 1 skipped).
- **W1** `20261008T033458Z_IMP-3-item305-W1-d2rod-r4s2.log` (**459 s**, `-n 8`, 13 passed, `[capture] rc=0`, in [450, 570]): rod arm at r/(4√2), degree 2, **276 088 cells** — (G1) < 320 000; (G2) all four N_i 5.600000000000e-05; (G3) 6 256 tag-402 facets, all exterior; (i) 2.829e-15; (ii) σ_max 0.999828051, column powers ≤ 0.834917; (PP) ≤ 3.4e-15 per drive; (M) ≤ 0.0221 % (band 0.50 %); (E) 31 280 tag-402 dofs, max |E| 0 on all four drives; (V) True; **(D) 1 775 954 = 2E + 2F exact** (E 331 014, F 556 963).
- **Controls:** (Kc) `tip_size_m` 0.00017677669529663688 == `TIP_R4S2` bitwise, cells 276 088 > 252 118 (1.0951); (Rb) η_rod^{r/4}(T1) P1 recomputed from the I14 JSONs reads "3.482328", string-equal; sizes 201 072 / 1 284 856 and 252 118 / 1 620 748 equal. Mesh identity against GEO-36 W2: cells, 402 facets, d2 dofs and near-T1 (19, 6.422237e-04 m) all equal. N_ball 32 826, mean ball edge 2.433349e-04 m.
- **Readings (printed, never gated; P1 of record):** **η_rod(T1) r/4 3.482328 → r/(4√2) 3.556262, move +2.123 %** (beside I14's p-move −5.995 % and remesh move −3.122 %; no step ratio — it needs r/8). F1 / F2 moves +0.027 % / −0.116 % (this remesh's floor). (C-h) 2.123e-2 ≥ 1e-3: IN. (C2) r/4 16.321, r/(4√2) 16.773 (≥ 3: IN). Near-T1 edge −8.033 %, local ratio h(r/4)/h(r/(4√2)) 1.0874 (the requested VIn ratio 1.4142 is not the step at T1). T1 cell (original index 187 478) mean edge 7.732801e-04 m, a first record. Q P1 top 4.463326 → 4.468599 (+0.118 %), bottom +0.112 %; every drive moves ≤ +0.12 %. Rod-arm 1 g masses 0.9817 / 0.9820 g. Closed form / η ratio at T1–T3, S1: 0.825 / 0.911 / 1.020 / 1.080 — no accuracy claim.
- **Price (re-prices r/8 for the weekly):** build 74.92 s, sweep 267.03 s, 12 SAR balls 35.52 s, arm 308.0 s, body 383.2 s, summed peak RSS 53.85 GiB at `-n 8` — under GEO-36's scaled row (sweep 318.8–322.9 s, 57.04–57.77 GiB). One ladder step, no rate, no accuracy claim.

### §9 item 327 — 2026-10-08 (chain I15 part 3a — the r/8 knob, xl entry 33's prerequisite)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.*

**What it is.** Add `r8` to `FEM_EM_IMP3_STEP7` in `tests/validation/test_imp3_tip_ladder_1g.py`, so `xl` entry 33 (`docs/testing/xl-pending.md`, entry 33) has a module to run. No solve. Item 305's knob set is `KNOBS = ("r4s2",)` (`:68`) and the dispatch at `:240–244` calls `_r4s2()` only.

**How (G1) is applied: a module global, not a parameter.**
- `_arm_d(built, label, rod, degree)` (`test_imp3_tip_ratio_1g.py:150`) does `assert ncells < ROD_CELL_GUARD_W2` at `:165`, when `rod` is true, after `build_four_port_sweep(build_only=True)` and before the solve.
- The name resolves at call time in `test_imp3_tip_ratio_1g`'s namespace. It is imported at `:74` from `test_imp3_rod_birdcage.py:115` (`320_000`, "W2 rod-arm guard", §9 item 174 (A0)). It is a resource pre-check, not a physics band.
- `rod=False` is not a way round it: it also drops `wall_dofs`, `wall_max`, `case` and `built` (`:227–234`), which `_identities(a, True)` needs for (M) and (E).
- **Proposal:** the r8 path only does `monkeypatch.setattr(tests.validation.test_imp3_tip_ratio_1g, "ROD_CELL_GUARD_W2", ROD_CELL_GUARD_XL)`. This changes neither file nor the constant's definition, and pytest restores it at teardown. The weekly ruled exactly this (`xl-pending.md` entry 33 Price: "the window's module raises the guard for the knobbed rung only, disclosed, never for the default path"). So this item is **not** "write the anchor".
- **The constant:** `ROD_CELL_GUARD_XL = 400_000`, defined in the new module with a comment calling it an `xl` resource guard, not a band. It is derived from item 305 W1's measured 53.85 GiB summed at `-n 8` for 276 088 cells (`20261008T033458Z_…W1-d2rod-r4s2.log`).
  - At the guard, ×1.449 cells gives 78–89.5 GiB at e = 1.0 / 1.37.
  - r/8 (342 072, ×1.239) gives 66.7–72.2 GiB. That is far inside 512 GiB, and r/8 sits at 0.855 of the guard.
  - r/16 (≈ 1 M cells, `GEO-36.md:236`) still trips it, so the guard keeps its purpose of catching a wrong rung or a runaway mesh before the solve.
- The r8 path prints one disclosure line: patched value, original 320 000, measured cells.

**Code (additive).**
- `KNOBS = ("r4s2", "r8", "r8build")`. Add `TIP_R8` to the `test_birdcage_implant_rod_tip` import. The test takes `monkeypatch`.
- A new `_r8(build_only)` beside `_r4s2()`. **`_r4s2()`'s body is not touched.** Duplicate rather than refactor. `git diff` on the module may delete only the `KNOBS` line, the docstring's knob sentence, the test signature and the dispatch lines.
- `_r8` reads `ours_imp3g_d2off.json`, `ours_imp3g_d2rod.json` **and** `ours_imp3h_r4s2.json`, and writes `logs/imp3/ours_imp3h_r8.json` with mode `"x"`.
- GEO-36 W2 r/8 comparands become constants with log-line cites from `20261007T170700Z_GEO-36-step4-W2-rungs.log`: cells 342 072, d2 dofs 2 200 500, 402 facets 8 790, near-T1 (26, 5.08e-4 m, exact digits from the log), N_ball 93 132 / 1.722e-4 m. Add item 305's `ITEM305_ETA_T1_P1 = "3.556262"`, sizes 276 088 / 1 775 954, and the price row (sweep 267.03 s, body 383.2 s, 53.85 GiB at `-n 8`). All of these are printed comparands except where asserted below.
- **Readout:** η_rod(T1) P1 at r/4 → r/(4√2) → r/8. Print Δ₁ (+0.073934), Δ₂, their signs and \|Δ₂\| / \|Δ₁\| beside entry 33's 0.7 threshold (0.7 \|Δ₁\| = 0.051754). Printed, never asserted; the 10-14 weekly reads the rule. Also print the near-T1 edge ratio h(r/(4√2)) / h(r/8) beside 1.264, and Q at three rungs.
- `r8build` builds the r/8 fixture and runs (Kc8), (G3), the census, `_ball_and_t1` and the T1 owner on every rank. It does not call `_arm_d` and writes no JSON.

**Anchor (asserted).**
- (Kc8) `tip_size_m == TIP_R8` bitwise, and rod cells > 276 088.
- (G1x) inside the patch, `_ratio.ROD_CELL_GUARD_W2 == ROD_CELL_GUARD_XL` and cells < it.
- (G3) as item 305.
- (Rb8) η_rod^{r/(4√2)}(T1) P1 recomputed from `ours_imp3h_r4s2.json` and the off JSON is string-equal to "3.556262", with sizes equal as integers.
- At `r8build`: mesh identity cells == 342 072, 402 facets == 8 790 (printed; predicted equal).

**Negative controls (rule (e)).**
- (N1, *asserted*) r/8 cells > the unpatched `ROD_CELL_GUARD_W2`, so the unpatched (G1) would trip. Predicted 342 072 / 320 000 = 1.069.
- (N2, *asserted*) the r4s2 dispatch never enters the patch: `_ratio.ROD_CELL_GUARD_W2 == 320_000` is checked at the top of every test call before dispatch.
- (N3, *predicted*) an unset knob skips, as W0 shows.

**Tier / cost.** Smoke W0, then standard W1b. Measured basis: item 305's cache sweep 93 s and W0 71 s; GEO-36 W2's r/8 build 86.71 s at `-n 2`.
- Cache sweep (rubric label `IMP-3-item327-cache-sweep`), after `pgrep -c python3` reads 0.
- **W0** (knob unset, `-n 1`, ≈ 70 s, range 30–120):
  `scripts/testing/run_and_log.sh IMP-3-item327-W0-import "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short'"`
- **W1b** (`r8build`, `-n 2`, ≈ 160 s, range 110–300):
  `scripts/testing/run_and_log.sh IMP-3-item327-W1b-r8build "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP7=r8build timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short'"`
- Total ≈ 93 + 70 + 160 ≈ 325 s ≈ 5.4 min of compute; predicted ≈ 25 slot minutes.
- **Never run `r8` in-slot.** At `-n 8` its scaled body is ≈ 475–515 s plus the 87 s build, and its window is over 590 s.

**Traps.**
- Item 305's list carries in full (`IMP-3.md:1050–1068`).
- `monkeypatch` must target the **module object `_arm_d` lives in**. Patching `test_imp3_rod_birdcage.ROD_CELL_GUARD_W2`, or the new module's own imported name, has no effect on `:165`.
- Never import `test_birdcage_implant_rod_tip`'s `builds` fixture.
- `ours_imp3h_r8.json` uses mode `"x"`: a pre-existing file fails the `xl` window, so the slot confirms it is absent.
- The `xl` container shares the `..:/workspace` mount (`docker/docker-compose.yml:54`), so the r4s2 JSON is visible there.
- A **cold FFCx cache in the `xl` container** pays degree-2 JIT inside the 4 h window. That is a price line, not a defect.

**Scope.** Covers the knob, the guard path and the r/8 build identity. It does not cover the solve, any η / Q reading, convergence, a bracket, the ball-radius knob (I17), `src/`, any band, I14's module or `ROD_CELL_GUARD_W2`. `git diff 0f6bbbc7 -- tests/validation/test_imp3_tip_ratio_1g.py` and `git diff -- tests/validation/test_imp3_rod_birdcage.py` must both stay empty.

**Negative result.**
- W0 red: a module defect. Fix it in-slot and re-run W0.
- W1b (Kc8) / (N1) / (G3) red, or a gmsh FAIL: report every number, park on `attempt/IMP-3-step7-<ts>`, mark BLOCKED, and open a `GEO-36` known-issues entry. Entry 33 stays PENDING.
- A mesh-identity mismatch with GEO-36 W2: printed and reported. It does not block. The window's (Kc8) and (G1x) still hold.
- No `PREDICTED RED` line.

**Status it can move.** `xl` entry 33 goes from PENDING to **READY** in the landing commit. The implementer does not write the command; the review issues it in entry 30's shape and queues it for Wednesday 10-14 02:00. §7 `IMP-3` stays ✅ with no record sentence. The slot makes no §10 edit.

**Review rulings (2026-10-08 03:00 daily review, item 327).** (1) **The monkeypatch is the disclosed raise entry 33's price paragraph licenses** ("raises the guard for the knobbed rung only, disclosed, never for the default path"): on the `r8` / `r8build` paths only, printed, undone at teardown; I14's module and `ROD_CELL_GUARD_W2` are not edited. (2) Keep `r8build`. (3) `ROD_CELL_GUARD_XL = 400_000` is accepted as a **resource guard** (the 128 G-container sizing logic re-applied to the `xl` window), not a physics band. (4) `FEM_EM_SOLVER_PROGRESS=2` stays in entry 33's command only if the n-port sweep honours it — the executor reports which; the review issuing the command decides. (5) The READY note states that r/8 needs `xl` for the 590 s ceiling, not for memory (≈ 67–72 GiB predicted at `-n 8`). (6) The executor copies the r/8 near-T1 digits with line numbers from `20261007T170700Z_GEO-36-step4-W2-rungs.log`, never from `GEO-36.md`. The implementer landing 327 marks entry 33 READY and does not write the command.

**Outcome (2026-10-08 07:30 CDT slot, implementer): GREEN — every anchor and both asserted negative controls green on W1b; mesh identity with GEO-36 W2 exact.** `pgrep -c python3` 0; cache sweep `20261008T123307Z_IMP-3-item327-cache-sweep.log` (85 s, no orphan); W0 `20261008T123440Z_IMP-3-item327-W0-import.log` (5 s, 4 passed, 1 skipped — (N3) as predicted; (N2) ran before the skip; below the 30–120 s range, a warm cache). `ours_imp3h_r8.json` absent; both helper diffs empty; the module diff deletes only the docstring knob sentence, the `KNOBS` line, the test signature and the dispatch line.
- **W1b** `20261008T123452Z_IMP-3-item327-W1b-r8build.log` (**88 s**, `-n 2`, 5 passed; under the 110–300 s range): (Rb8) "3.556262" string-equal, sizes 276 088 / 1 775 954 and 201 072 / 1 284 856 equal (`:64`); (Kc8) `tip_size_m` 0.000125 == `TIP_R8` bitwise, 342 072 cells > 276 088 (1.2390) (`:1952`); DISCLOSURE line (`:1953`); (G1x) patched 400 000, 0.855 of the guard; (N1) 1.069 × 320 000 (predicted 1.069); (G3) 8 790 tag-402 facets, all exterior (`:1954–1956`).
- **Mesh identity vs GEO-36 W2** (r/8 comparands copied from `20261007T170700Z_GEO-36-step4-W2-rungs.log:7500–7501, 7514`): cells 342 072, 402 facets 8 790, d2 2 200 500 = 2E + 2F (E 410 054, F 690 196), near-T1 26 cells / 5.081837e-04 m, N_ball 93 132 / 1.722398e-04 m — all equal (`:1957, 1960`). Near-T1 edge ratio h(r/(4√2))/h(r/8) 1.2638 beside 1.264; T1 cell original index 163 123, mean edge 7.850364e-04 m (a first record at r/8). Build 83.57 s, summed peak RSS 0.87 GiB at `-n 2`.
- **Ruling (4):** the n-port sweep honours `FEM_EM_SOLVER_PROGRESS` (code reading, not run): the factorising first drive's `TimeHarmonicSolver.solve` sets `mat_mumps_icntl_4` (`time_harmonic.py:580–583`); drives 2–4 use `solve_with_held_factorization`, which prints `[solve] held-factor back-substitution done in … s` (`:715–723`).

### §9 item 381 — 2026-10-11 (chain I21b — the rod under the tuned copper F-human's two-port drive; chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 I21b (`PROJECT_PLAN.md:4753–4761`, the 2026-10-10 weekly), lane 3 / Phase 7 (7b), serial on I21a; §9 rule (4).
- **Target: identities.** The tip ratio is a record with no published comparand until (7c).
- **No band, tolerance, Done-when or validation target is set or moved; every band is imported.** It is a record on a ✅ row.
- The 2026-10-02 sheets-and-κ ruling, the capacitor-ESR ruling and the 10-10 10:30 rulings on items 358 / 359 bind (`WF-10.md`, "Review rulings (2026-10-10 10:30 review)"). So do the 2026-10-10 18:00 rulings on item 361 (`MAT-15.md`, 361-1 … 361-5), which govern how 358's builder is reused.

`IMP-3` step 8: `WF-9`'s lead-in-phantom construction on H9′'s tuned copper F-human. The fixture is:
- I21a's kept-rod hole mesh, at degree 1, 64 MHz, κ off;
- r = 0.2 Ω on all 32 gaps, with the two Thevenin feeds P17 / P21 at C_in,F-human = 118.2437 pF;
- the rod switched by material: σ = 800 S/m (`SIGMA_WIRE_S_PER_M`) against the gel.

**Executor and cost.**
- **Executor:** implementer. **Tier:** heavy, `-n 8`, complex build.
- **Windows:** cache sweep, W0 (smoke, `-n 1`), W1a (rod-free reproduction) and W1b (rod arms), all foreground.
- **Class and cost:** chain step (solve); ≈ M measured / ≈ 55 predicted min. That is 15 fixed + ≈ 25 code + ≈ 10 compute (W1a ≈ 3 min, W1b ≈ 5–7 min) + spare. **If W1b has not started by minute 45, park** (rule (d)).
- **Depends on item 372 (I21a) landing.** This means `tests/mesh/test_birdcage_implant_rod_f_human.py` exists with `ROD_FH` and `F_HUMAN_KEPT_ROD_CELL_RECORD`, and `GEO-36-step6-W2` is green in `test-results.md`. **If it did not land, skip this item; it stays ⬜.**
- **Gitignored inputs, both on disk when this brief was drafted:**
  - `logs/wf10/h9p-node-118.2437.json` (5 468 bytes, rewritten 2026-10-10 19:50 CDT by item 361's builder call);
  - `logs/wf12/h12c-cu-s32-a2.npz` (20 812 bytes; `build_tuned_quadrature` asserts it, `test_wf10_f_human_tuned_b1_sar.py:208–209`).

  **If either is absent, stop and mark the item BLOCKED; never re-run 358 or re-solve the 32×32 in-slot.**
- **Family:** `IMP-3` step 8, attempt 1 of 4.

**The reproduction question, answered: the rod-free mesh is not item 358's mesh.** The weekly's "the rod-free drive reproduces item 358's node record at `REUSE_REPRODUCTION_RTOL`" is posed as follows.
- **On I21a's mesh, "rod-free" means the rod's cells switched to gel.** That is `WF-9`'s *without* arm, which the fold needs. It is a different triangulation from item 358's 369 387-cell hole mesh: ≈ +66 k cells, the r/3 grading out to 10 mm and the cap balls. Its node S therefore differs from 358's record by a discretisation move. On this family that move is O(1e-4 – 1e-2), so it cannot be held to 1e-6. It is **printed** as the remesh move (`_residual(s_without, s_rec)`), predicted in [1e-5, 1e-2].
- **The reproduction is posed on 358's own mesh, a separate build (W1a).** `build_tuned_quadrature(comm, with_z50=False)` builds `_build_f_human(comm, as_hole=True)` at 369 387 cells, exactly as items 359 and 361 did. (S_rep) is `_residual(s_node, s_rec)` ≤ `REUSE_REPRODUCTION_RTOL`, with s_rec read from the JSON **before** the builder rewrites it. The same comparison on the same fixture read 6.027e-13 (`20261010T171705Z_WF-10-step2s-W1.log:6312`) and 2.705e-13 (item 361).
- It is a separate window from the rod arms. The rod arms do not need anything W1a computes; they re-read the JSON themselves. Splitting keeps each window well inside the 590 s container ceiling (a single window would be ≈ 480–650 s).

**Helper edit (rule (a) / (c) style, additive, disclosed).** `_f_human_kwargs` and `_build_f_human` (`tests/validation/test_wf7_f_human_b1_quadrature.py:67, :101`) gain a keyword `implant_rod=None`.
- **`None`** adds no key and leaves the material map unchanged, so the call is byte-identical to today.
- **A rod** requires `as_hole=True`; anything else raises `ValueError` before gmsh, because the solid route is not gated. With a rod:
  - `kwargs["implant_rod"] = implant_rod`;
  - the material map gains `BIRDCAGE_IMPLANT_CELL_TAG: <the gel's HomogeneousMaterial>`. This is the *without* arm. The tag is never left unmapped: unmapped cells are vacuum, `WF-9`'s trap (`WF-9.md:262`);
  - the returned dict gains an additive key `problem`.
- **Re-gate.**
  - W0 asserts `_f_human_kwargs(as_hole=True) == _f_human_kwargs(as_hole=True, implant_rod=None)`, the same for `as_hole=False`, and that `"implant_rod"` is in neither.
  - W1a runs the default path through the edited `_build_f_human` on the H9′ consumer and asserts (A0) 369 387 and (S_rep) at 1e-6.
  - `_build_f_human`'s own gate module, `WF-7` step 2 (`test_wf7_f_human_b1_quadrature.py`, heavy `-n 8`), is **not** re-run. That is disclosed for the review to rule on: W1a exercises the edited default path at the same comparand grain.
- 358's, 359's and H17's modules are **not edited**, and neither is any `MAT-15` or `WF-9` module.

**New module `tests/validation/test_imp3_f_human_tuned_rod.py`.** It is `@complex_only`.
- **Knob `FEM_EM_IMP3_STEP8` ∈ {`repro`, `rod`}.** Unset ⇒ the gated test skips; any other value ⇒ `ValueError`.
- **Ungated test `test_imp3_s8_coding`** (numpy only):
  - the knob parser;
  - W0's kwargs identity;
  - `ROD_FH == dict(ROD, keep_cells=True)`;
  - `TIPS["top"]` = `ROD_CENTRE` + (0, 0, ℓ/2 + r) bit for bit;
  - LOAD = −jωμ₀ computed from `FREQUENCY_HZ` (not `WF-9`'s `OMEGA` constant), equal to `test_wf9_lead_in_coil_phantom.LOAD` to 1e-15 relative.

**Imports.** Names only, never a `test_*` symbol.
- **358's module** (`test_wf10_f_human_tuned_b1_sar.py`): `build_tuned_quadrature` (`:147`), `C_IN_F_HUMAN_PF` (`:58`), `NODE_JSON` (`:82`), `_sense_weights` (`:108`).
- **H17's module** (`test_wf10_f_human_copper_esr_root.py`): `ESR_OHM`, `_z_rc`, `_thevenin_feed_impedance`, `_source_voltage`, `_current_route_s`, `_active_y_in`, `Y0_CODING_BAND`.
- **`test_wf10_f_human_tuned_circuit.py`:** `QUADRATURE_SLOT_OFFSET` (`:52`), `_respec` (`:64`).
- **`test_wf10_f_human_tuned_sar10g.py`:** `AZ17_RECORD_DEG` (`:45`), `_lattice` (`:58`).
- **`test_wf7_f_human_sar10g_hotspot.py`:** `_phantom_geometry` (`:95`), `COVERAGE_BALL_PAD_M`, `CONDUCTOR_CLEARANCE_M` (`:65–66`).
- **`test_wf7_f_human_b1_quadrature.py`:** `FREQUENCY_HZ`, `DEGREE`, `_build_f_human`.
- **`test_wf12_f_human_copper_matrix.py`:** `COPPER_SIGMA`, `BLIND_FACTOR`.
- **`test_wf12_f_human_hole.py`:** `_census_checks`, `_owned_census`, `_rss_gib`.
- **`test_th14_birdcage_copper.py`:** `_leontovich_term` (`:160`), `_surface_loss_w` (`:174`).
- **`test_imp3_rod_birdcage.py`:** `TIPS` (`:140`), `SPHEROID_SEMI_MAJOR`, `SPHEROID_SEMI_MINOR`, `QDEG_CF`, `CF_CHUNK`, `AXIS_SAMPLES`, `PRED_Q_OVER_ESTIMATE` (`:138–159`).
- **`test_wf9_transfer_function_identity.py`:** `BALL_RADIUS` (`:72`, 3 mm), `_free` (`:165`).
- **`tests/mesh/test_birdcage_implant_rod.py`:** `ROD`, `ROD_CENTRE`, `_cell_values_with_ghosts` (`:171`).
- **I21a's module:** `ROD_FH`, `F_HUMAN_KEPT_ROD_CELL_RECORD`, `F_HUMAN_KEPT_ROD_CELL_GUARD`.
- **Bands and helpers:**
  - `DISCRETE_IDENTITY_RTOL`, `_sheet_field_dissipation_w`, `_source_power_w` (`test_birdcage_power_identity.py:227`);
  - `EXACT_IDENTITY_RTOL` (`test_birdcage_sar_mass_averaged.py:163`);
  - `REUSE_REPRODUCTION_RTOL` (`test_port_birdcage_ring_matrix.py:171`) and `_residual` (`test_port_circuit_layer_field.py:142`);
  - `UNIT_WEIGHT_CONTROL_RTOL` (`test_port_drive_superposition.py:303`) and `_ring_quadrature_slots` (`:1208`);
  - `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`);
  - `PHANTOM_RHO_KG_PER_M3` (`test_birdcage_b1_plus_map.py:139`);
  - `ONE_GRAM_KG`, `TEN_GRAM_KG`, `QUADRATURE_DEGREE` (`test_mass_averaged_sar_standard_masses.py:121–124`);
  - `PHANTOM_CELL_TAG` (`test_port_birdcage_lumped_column.py`);
  - `SIGMA_WIRE_S_PER_M` (`test_port_gap_voltage_impedance.py:235`);
  - `AZIMUTH_MATCH_DEG` (`test_port_birdcage_ring_column.py`);
  - `REFERENCE_IMPEDANCE_OHM`, `TERMINATED_PORT_IMPEDANCE_OHM`;
  - `averaging_ball_radius`, `build_density_field`, `mass_averaged_sar`, `mean_sar` (`src/fem_em_solver/post/sar.py`);
  - `evaluate_vector_field_parallel`, `pec_prolate_spheroid_field`, `BIRDCAGE_IMPLANT_CELL_TAG`, `MU_0`.

**Lifted, never imported (provenance comment: base commit and line span).**
- **(L1)** 358's `_drive`, `_pair` and `_partition` closures (`test_wf10_f_human_tuned_b1_sar.py:245–349`) and its order block (`:305–326`), verbatim except one change: `_drive` takes the solver as an argument instead of `ctx["solver"]`. They are closures inside `build_tuned_quadrature`, so they cannot be imported. 358's module stays unedited (rule (a) cascade into 359 / 361 / `EX-75`).
- **(L2)** `_ball_mean_closed_form` (`test_imp3_rod_birdcage.py:229–272`), with one added filter: only owned cells whose `_cell_values_with_ghosts` value is `PHANTOM_CELL_TAG` contribute. On F-small the hole excluded the rod by construction; here the rod's 403 cells exist and must be excluded the same way.

**W1a — `repro` (rod-free reproduction on 358's own mesh).**
1. Read 358's JSON (rank 0, `bcast`). Assert `c_in_pf == C_IN_F_HUMAN_PF` exactly and keep `s_rec`.
2. `out = build_tuned_quadrature(comm, with_z50=False)`.
3. **(S_rep)** ≤ `REUSE_REPRODUCTION_RTOL`.
4. Re-assert from `out["anchors"]`: (A0) 369 387 within `CELL_COUNT_BAND` with the census clean; (P_RC) ≤ `DISCRETE_IDENTITY_RTOL`; node NEG ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL`; (P_w) both senses ≤ `DISCRETE_IDENTITY_RTOL`. The builder already asserts them; this re-states them in this module's log.
5. **Printed:** SAR10g(c\*) under ccw on this mesh (σ = the builder's field, as 359 did) beside 359's 4.852895959e-03 W/kg (`…step2s-W1.log:6372`), with c\* the 18-point lattice point "r=45 dz=+0 az17+0" of `_lattice(az₁₇, z_c)`. Predicted equal to ≤ 1e-9 relative (same mesh, same drive). It is a reading, not an anchor.

**W1b — `rod` (the two arms on I21a's mesh).**
1. Read the JSON first (as W1a).
2. **Build.** `built = _build_f_human(comm, as_hole=True, implant_rod=ROD_FH)`.
   - Check the guard (`F_HUMAN_KEPT_ROD_CELL_GUARD`) first.
   - **(A0r)** cells vs `F_HUMAN_KEPT_ROD_CELL_RECORD` within `CELL_COUNT_BAND`; `_census_checks` clean; tag 403 > 0 owned cells.
3. **Problems.** `gel = built["problem"].material_map[PHANTOM_CELL_TAG]`.
   - `p_without = built["problem"]`.
   - `p_with = dataclasses.replace(p_without, material_map={**map, BIRDCAGE_IMPLANT_CELL_TAG: HomogeneousMaterial(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r)})`.
   - Δσ = `SIGMA_WIRE_S_PER_M − gel.sigma`, computed and never typed.
   - **(M0)** DG0 σ and εᵣ on owned 403 cells, rank-reduced min = max: 800 / gel εᵣ (with) and gel σ / εᵣ (without). Use `WF-9`'s `_material_minmax` pattern (`test_wf9_lead_in_coil_phantom.py:229–248`).
4. **Masks** (DG0 from numpy midpoints over owned + ghost cells, `WF-9`'s `_build_case` pattern `:174–195`):
   - χ_lead = tag-403 cells;
   - χ_b = tag-3 cells whose midpoint lies within `BALL_RADIUS` of `TIPS["top"]`;
   - V_b = ∫χ_b dx;
   - χ₃ = tag-3 cells.
   - σ₃ = a DG0 field equal to `gel.sigma` on tag 3 and 0 elsewhere, used for every SAR ball and the coverage (see Traps).
   - Asserted: χ_b and χ_lead are disjoint, and both counts are > 0.
5. **The without arm first.**
   - `TimeHarmonicSolver(p_without, degree=DEGREE)`.
   - The node pair at C_in via (L1): P17 factorises, P21 is held, fields kept.
   - The P17 drive's partition gives **(P_RC)** and its NEG.
   - The ccw / cw superpositions via 358's `:414–461` body give **(P_w)** on both senses; (Lin) is printed.
   - Then `_free(solver)`.
6. **The with arm.**
   - The same solver and node pair on `p_with`, then one more held solve on ℓ = `ufl.inner((χ_b / V_b)·ẑ, v)·dx` (`WF-9` `:414`), giving w on the same factorisation.
   - (P_RC), NEG, (P_w) and (Lin) as above.
   - Then `_free(solver)`.
   - Only one factorisation is alive at a time.
7. **(T3), per node drive d ∈ {P17, P21}.**
   - q(e) = (1/V_b)∫χ_b e_z dx and Δq_d = q(E_with,d) − q(E_without,d).
   - fold_d = LOAD·Δσ·∫χ_lead w·E_without,d dx (`ufl.dot`, never `inner`), as `WF-9` `:456–459`.
   - Asserted: \|fold_d − Δq_d\| ≤ `DISCRETE_IDENTITY_RTOL`·\|Δq_d\|. The ccw value follows by linearity and is printed.
8. **(Σ_P), both arms, ccw:** `mean_sar(e_ccw, sigma=<solver σ field>, rho=PHANTOM_RHO_KG_PER_M3, cell_tags=…, subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]` against the partition's `p_phantom`, ≤ `EXACT_IDENTITY_RTOL`.
9. **Coverage, both arms, ccw** — 359's anchor (i) (`test_wf10_f_human_tuned_sar10g.py:258–269`), with **σ₃** in place of the solver σ.
   - The whole-phantom ball (radius hypot(R, H/2) + `COVERAGE_BALL_PAD_M` < `F_HUMAN_RING_RADIUS − CONDUCTOR_CLEARANCE_M`) has dissipated power = `mean_sar(…, sigma=σ₃, subdomain 3)` and mass = ρ·V(tag 3), both ≤ `EXACT_IDENTITY_RTOL`.
   - ρ = `build_density_field(…, density_map={t: 0 for t ≠ 3})` (359's construction; it zeroes 403 too).
10. **Readings** (printed, never asserted):
    - **Tip 1 g ratio.** Q_tip,d = SAR1g(tip; with) / SAR1g(tip; without) for tip ∈ {top, bottom} and d ∈ {ccw, P17, P21}, with σ₃ and ρ in both arms and a₁g = `averaging_ball_radius(ONE_GRAM_KG, …)`. Ball masses are printed and are identical in both arms.
    - **The local background field's `IMP-2` estimate,** `IMP-3` step 1's construction (`test_imp3_rod_birdcage.py:379–394, :604`):
      - η_z,d = Σ\|E_z\|² / Σ\|E\|² over `AXIS_SAMPLES` points on the rod axis (x = 0.015 m, z ∈ [−0.016, 0.016]) of the **without** arm's field, read through `evaluate_vector_field_parallel` with the rank-0 `bcast` of `test_imp3_tip_ratio_1g.py:138–147`;
      - M_cf(tip) from (L2) at a₁g;
      - est = η_z·M_cf + (1 − η_z);
      - Q / est is printed beside `PRED_Q_OVER_ESTIMATE` (0.6, 2.5).
    - **The rod's effect on SAR10g(c\*).** SAR10g(c\*; with) / SAR10g(c\*; without) − 1 under ccw, with σ₃, with c\* as in W1a. The without-arm value is printed beside 359's 4.852895959e-03 as the remesh move.
    - \|q_with / q_without\| per drive.
    - `_residual(s_without, s_rec)` (remesh move) and `_residual(s_with, s_without)` (the rod's effect on the coil's 2×2).
    - P_rod = ½∫σ\|E\|² over 403 in the with arm, and P_rod / P_phantom.
    - Cells, E, degree-1 dofs, phase times and `_rss_gib`.

**Anchors (asserted; every band imported by name, none restated).**
- **W1a:** (S_rep); (A0); (P_RC); node NEG; (P_w) ccw / cw.
- **W1b:**
  - (A0r); (M0); mask disjointness;
  - (P_RC) and node NEG per arm;
  - (P_w) ccw / cw per arm;
  - (T3) on P17 and P21;
  - (Σ_P) per arm;
  - coverage power and mass per arm.

**Bands, read against rule (f).** None was registered by measurement on this fixture.
- `DISCRETE_IDENTITY_RTOL` (1e-6), `EXACT_IDENTITY_RTOL` (1e-10) and `REUSE_REPRODUCTION_RTOL` (1e-6) are round-off bands that 358 / 359 / 361 and `WF-9` assert.
- `CELL_COUNT_BAND` is the record band.
- Nothing here is compared with an analytic comparand that the fixture's mesher cannot meet.

**Negative controls (rule (e)).**
- **Asserted (backed, same fixture class):**
  - The node NEG, the P_surf-blind miss of (P_RC), at ≥ `BLIND_FACTOR` × the band. It read 7.0125e-3 on 358's mesh (`…step2s-W1.log:6304`). The rod-mesh arms are predicted ≈ 7.0e-3, since the rod barely moves the copper loss. The ceiling is P_surf/P_src ≈ 7e-3.
  - **W1b's own NEG on the rod mesh has no prior measurement.** It is asserted at the same ≥ `BLIND_FACTOR` × band bar, which sits ≈ 700× under the prediction. That is the bar 358 / 361 assert on a rod-free mesh. The review may relabel it printed.
- **(C1) Hermitian read-back:** the fold with conj(w). **(C2) Total-field fold:** the fold with E_with in place of E_without. Both are asserted only as wired (> `DISCRETE_IDENTITY_RTOL`), as `WF-9` step 1 did (`WF-9.md:248–251`). The factors are **predicted** and printed: C1 in [10 %, 200 %] (ceiling 200 %); C2 in [80 %, 105 %] (`WF-9` step 1 read 94.6 / 97.3 % on F-small, a different fixture).
- **(C3, predicted, printed) Solver-σ coverage:** the with-arm whole ball computed with the solver's σ field instead of σ₃ misses by ≈ P_rod / P_phantom. That is the trap the σ₃ construction avoids.

**Predictions** (written first in the module, printed):
- (S_rep) ≤ 1e-12.
- (P_w), (P_RC) ≤ 1e-12. (T3) ≤ 5e-12 (`WF-9` step 1: ≤ 3.7e-12). (Σ_P) ≤ 1e-14. Coverage ≤ 1e-13.
- Cells = I21a's record. Tag-403 cells ≈ 15 000. Degree-1 dofs ≈ E of I21a's printout.
- Q_top,ccw in [1.5, 8]: `IMP-3`'s PEC Q on F-small's 1 g rung read 4.46–4.50 (P1 top), and `WF-9`'s σ = 800 S/m 3 mm-ball field ratio read 2.46–2.61.
- η_z,ccw in [0.2, 1.0].
- Q / est in (0.6, 2.5).
- \|SAR10g(c\*) with / without − 1\| in [1e-4, 3e-2]: c\* sits ≈ 30.4 mm from the rod axis, and its 10 g ball's edge ≈ 17 mm.
- The remesh move `_residual(s_without, s_rec)` in [1e-5, 1e-2]; `_residual(s_with, s_without)` ≤ 1e-3.
- 1 g ball masses ≈ 0.98 g at both tips (`IMP-3` I14's rod arm: 0.982 g).

**Windows.** Foreground. The harness call is the whole Bash command, bare. Bash `timeout: 660000`; `-s`; `tests/environment/test_complex_mode.py` first.
- **Price, measured on this route.**
  - Item 361's W1: the builder with `with_z50=False` took 152.7 s (build 112.6 s, node pair 34.4 s, superposition 6.5 s); harness 162 s; summed RSS 9.59 GiB at `-n 8` (`20261011T003520Z_MAT-15-step5-W1.log`).
  - Item 359's 10 g balls: mean 3.44 s, max 4.69 s (`…step2s-W1.log:6384`).
- **W1a**, predicted 150–220 s.
- **W1b**, predicted 270–450 s, ≤ 560 s on the 1.32× box. Its components:
  - the rod build ≈ 125–140 s;
  - two node pairs at ≈ 1.24× 361's (≈ 43 s each, scaled by (435/369)^1.3);
  - two superpositions ≈ 13 s;
  - ℓ held ≈ 3 s;
  - ≈ 16 balls ≈ 55 s;
  - two coverage balls ≈ 20 s;
  - (L2) ≈ 5–10 s;
  - pytest start ≈ 15 s.

  Summed RSS ≈ 10–13 GiB, with one solver alive at a time.
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Then run the any-orphan sweep as one bare call:

  `scripts/testing/run_and_log.sh IMP-3-step8-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **W0** (smoke, `-n 1`, numpy only): `test_imp3_s8_coding` green and the gated test skipped.

  `scripts/testing/run_and_log.sh IMP-3-step8-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_f_human_tuned_rod.py -v -s --tb=short'"`
- **W1a** (heavy, `-n 8`). §5.1's durable capture idiom is copied verbatim, trailing `; exit $rc` included, as in 358 / 359 / 361:

  `scripts/testing/run_and_log.sh IMP-3-step8-W1a-repro "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf10 && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP8=repro PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_f_human_tuned_rod.py -s -v --tb=short > /workspace/logs/imp3-step8-w1a-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp3-step8-w1a-raw.log; cat /workspace/logs/imp3-step8-w1a-raw.log; exit \$rc'"`
- **Before W1b:** `pgrep -c python3` reads 0.
- **W1b** (heavy, `-n 8`): the same command with label `IMP-3-step8-W1b-rod`, `FEM_EM_IMP3_STEP8=rod` and raw file `/workspace/logs/imp3-step8-w1b-raw.log` (all three occurrences).
- **No edit after W0 starts** (rule (i)). Any edit re-runs W0, W1a and W1b. W1a and W1b plus one re-run of either fit the hour. No `src/` file is edited. The helper edit's re-gate is W0 plus W1a, as stated above.

**Traps already paid for.**
- **Orphan `.c` of any size** (item 160, 44 KB). Sweep only when `pgrep -c python3` reads 0, and delete orphan stems only. After a kill or status 124: run `pgrep -c python3` first, then the sweep, then park. Never raise the timeout or `-n`.
- **Complex build.**
  - No `ufl` ordering on complex operands and no `sqrt` in a comparison (`OPS-22`; `WF-6` 4c).
  - Masks are DG0 from numpy midpoints, never UFL conditionals.
  - Ball centres go in as plain floats.
  - Call `np.real` before every `float`.
- **`ufl.inner` conjugates its second argument.** q, fold, C2 and TF read with `ufl.dot` or components. Only ℓ's linear form uses `inner` against the test function, as `WF-9` `:414` does.
- **Unmapped kept cells silently become vacuum** (`WF-9.md:262`). (M0) is the guard, and the helper maps 403 to gel.
- **Every sheet's bilinear term and the Leontovich wall go into both factorisations.** A held solve under a different termination solves a different problem silently (`PORT-19`, `WF-12`). ℓ is held on the with-arm factorisation, which the P17 node solve made.
- **The kept σ = 800 S/m rod is inside every whole-phantom ball and inside the tip balls' reach.** Every SAR ball and the coverage use **σ₃**, never the solver's σ field. The F-small `IMP-3` ratio had a hole there and needed nothing.
- **`_ball_mean_closed_form` integrates every owned cell in the ball.** (L2) restricts it to tag 3. The `ROD_CENTRE` it pins is the rod's centre only because I21a placed the rod at `ROD`. If the review moves the rod, (L2) must take the centre as an argument.
- **Rank-local values.** `cell_tags.values` and `assemble_scalar` are rank-local; reduce everything. Every rank enters every collective (item 227's deadlock). Point evaluation goes through `evaluate_vector_field_parallel` with `bcast`.
- **Sense and plane.** ccw = (1, −1j) comes from the convention code (`_sense_weights`, 10:30 ruling (1)). Never re-choose it, and never recombine fields with hand-typed weights.
- **The builder rewrites `logs/wf10/h9p-node-118.2437.json`.** Read it before W1a's builder call. W1b reads it before anything. Write nothing else outside the harness logs and the raw capture files.
- **Memory.** Free each solver (`_free`) before the next factorisation. Two heavy windows never run at once.
- **Harness hygiene.** Calls are bare and whole; no `-k` inside the quoted command; no pipe; `-s` always; foreground with `timeout: 660000`. Never delete a log; never path-checkout a record file.

**Scope.** Covered:
- I21a's kept-rod copper F-human hole mesh at degree 1, κ off, r = 0.2 Ω on all 32 gaps, the two-port 90° drive at C_in,F-human = 118.2437 pF, 64 MHz;
- one σ = 800 S/m capped rod (r = 1 mm, ℓ = 30 mm, axis ẑ, centre (0.015, 0, 0)) in the 0.015 m cylindrical gel phantom;
- `WF-9`'s exact fold identity and H9′'s power identities on this fixture.

**Not claimed:**
- an accurate tip ratio: no published comparand until (7c), and `ANS-10` is the accuracy route;
- a PEC or insulated lead, or a real lead;
- ISO 10974 or ASTM F2182;
- heating;
- absolute SAR or compliance;
- degree 2 or h-convergence (the 5.50 % human-scale order caveat applies, `WF-7` step 0b);
- f₀, Q, match or a re-tune with the rod in (C_in is the rod-free root);
- C16 or the mirror;
- AED or Palace.

**Negative result.** Report and stop. Park on `attempt/IMP-3-step8-*` and mark item 381 BLOCKED in the same commit (rule (d)).
- **(S_rep) or a W1a builder anchor red:** the copper route or 358's fixture moved since 361. Read nothing else. File a known-issues entry; never re-band.
- **(A0r), (M0) or a mask check red:** the fixture is wrong (I21a's mesh, or the helper's material map). Stop.
- **(P_RC), (P_w) or a NEG red on a rod arm, green on W1a:** the rod's partition accounting (P_rod in P_rest), or a termination mismatch between arms. File a known-issues entry.
- **(T3) red with the power identities green:** read the fold's premise first. A = Aᵀ with the wall and the 32 sheets is assumed, not asserted. Print the (A1) Frobenius asymmetry diagnostic of `WF-9` `:360–389` in a follow-up, never in-slot. File a known-issues entry.
- **(Σ_P) or coverage red:** σ₃ / ρ construction. File a known-issues entry.
- **(C1) or (C2) not wired:** w or the fields are swapped. Stop.
- **Status 124 or > 590 s:** BLOCKED without a result; the killed window is the price record.
- Every Q, η_z, M_cf, estimate, ΔSAR10g(c\*), remesh move, P_rod and mass is a reading, never a red.

No `PREDICTED RED:` line. Every asserted anchor is predicted green; a red anywhere stops the slot.

**Status it can move.**
- **`IMP-3` stays ✅.** Its Done-when is "identities green with the rod in on `GEO-36`'s 167 371-cell record mesh ⇒ 🟡; the same identities green on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture" (`PROJECT_PLAN.md:2569`). This step reads nothing that touches it, so rule (j) does not bind and there is no Done-when amendment.
- **On green**, the row (or this chunk file, if the row is at 2 KB) gains a **"Step 8 record"**:
  - H9′'s tuned copper F-human with a σ = 800 S/m rod;
  - the fold identity and the power identities asserted on both arms;
  - (S_rep) on 358's own mesh;
  - printed: the tip 1 g Q against the `IMP-2` estimate, the rod's effect on SAR10g(c\*), and the remesh move.
- §10 I21b gains a DONE line, which the weekly marks.
- `WF-9` stays ✅. Its Done-when names F-small; F-human is a new record of its identity and moves nothing.

**Execution order.** pgrep → cache sweep → W0 → W1a → pgrep → W1b → commit (module, the helper edit, harness logs, the Step 8 record, the §10 I21b DONE line together) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 381 — `IMP-3` step 8 (§10 I21b) — per `docs/planning/chunks/IMP-3.md` §9 item 381. It depends on item 372: if `tests/mesh/test_birdcage_implant_rod_f_human.py` with `F_HUMAN_KEPT_ROD_CELL_RECORD` is not on main, skip it. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run the cache sweep, W0, W1a, then W1b, in order, foreground. Before W1b, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. Never delete a log. No edit to the module after W0 starts. Do not edit item 358's, 359's or H17's modules, any `MAT-15` or `WF-9` module, or any `src/` file; the only helper edit is the additive `implant_rod=None` keyword on `_f_human_kwargs` / `_build_f_human`, disclosed. If `logs/wf10/h9p-node-118.2437.json` or `logs/wf12/h12c-cu-s32-a2.npz` is absent, stop (BLOCKED); never re-run 358. Every rank enters every collective. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-3-step8-*`, mark item 381 BLOCKED in the same commit, and report."

**Review ruling 2026-10-11 (03:00 daily) on item 381.** (1) The reproduction is posed on item 358's own mesh (W1a, a separate `build_tuned_quadrature(comm, with_z50=False)` build); the rod-mesh rod-as-gel residual against the record is printed as the remesh move, never asserted. (2) **`_build_f_human` / `_f_human_kwargs` are not edited**: the rod-as-gel mapping lives in the item's own module (a copy-and-extend wrapper), so `WF-7` step 2's gate module is not touched and owes no re-run. (3) W1b's node NEG has no same-fixture prior: it is **printed, not asserted** (rule (e)). (4) SAR balls and coverage use the gel σ masked to tag 3, as briefed.

**Outcome (2026-10-11 15:00 CDT slot, implementer): ✅ — every asserted anchor green, `IMP-3` stays ✅.** New module `tests/validation/test_imp3_f_human_tuned_rod.py`; no `src/` file and no imported module edited. Per ruling (2) the rod-as-gel mapping is a wrapper in the item's own module: `_rod_kwargs` (pure, W0-gated) and `_build_f_human_rod`, which swaps `_f_human_kwargs` in the `WF-7` module namespace for the duration of one `_build_f_human(comm, as_hole=True)` call (restored in `finally`), then rebuilds the problem with tag 403 → gel and a fresh solver on it — disclosed here; `WF-7` step 2's module is byte-unchanged. Windows: cache sweep (no orphans) `20261011T200451Z_IMP-3-step8-cache-sweep.log`; W0 smoke `-n 1` 5 passed / 1 skipped `20261011T200637Z_IMP-3-step8-W0.log`; W1a heavy `-n 8` 156 s `20261011T200648Z_IMP-3-step8-W1a-repro.log`; W1b heavy `-n 8` 323 s, summed RSS 11.49 GiB `20261011T200932Z_IMP-3-step8-W1b-rod.log`. Measured ≈ 15 min.

**Step 8 record** (H9′'s tuned copper F-human, degree 1, 64 MHz, κ off, r = 0.2 Ω, C_in 118.2437 pF, one σ = 800 S/m rod kept as tag 403):
- **W1a (358's own mesh, 369 387 cells):** (S_rep) 2.658e-13 (`…W1a-repro.log:6297`); (P_RC) 7.807e-13, node NEG 7.0125e-3; (P_w) 4.860e-13 / 6.055e-13 (ccw / cw); SAR10g(c\*) 4.852895959e-3 W/kg, 4.3e-11 from 359.
- **W1b (I21a's mesh):** (A0r) 436 541 cells rel 0, tag 403 15 025; (M0) 800 / 78 with, 0.5 / 78 without; masks lead 15 025, ball 7 974, disjoint. Per arm (without / with): (P_RC) 1.067e-12 / 1.064e-12; (P_w) ccw 8.045e-13 / 7.858e-13, cw 7.493e-13 / 7.639e-13; (Σ_P) 0 / 0; coverage power 2.2e-16 / 4.9e-15, mass 1.2e-14 / 1.2e-14. (T3) 8.428e-13 (P17), 1.954e-13 (P21), ccw by linearity 2.357e-13 (`…W1b-rod.log:6361–6363`). Controls wired: C1 misses 1.9933 / 1.9932 (inside its printed (0.10, 2.00), at the ceiling), C2 0.9306 / 0.9305 (inside (0.80, 1.05)); C3 solver-σ whole ball misses 8.3070e-5 = P_rod / P_phantom, as predicted. Node NEG printed (ruling (3)): 7.0154e-3 / 7.0153e-3.
- **Printed readings:** tip 1 g Q (top / bottom): ccw 3.804 / 3.642, P17 1.777 / 1.648, P21 4.085 / 4.014; η_z (ccw / P17 / P21) 0.918 / 0.459 / 0.954; M_cf 2.905 / 2.903; Q / est ccw 1.384 / 1.325 (all six in (0.6, 2.5)); 1 g masses 0.9806 / 0.9821 g, identical in both arms. SAR10g(c\*) with / without − 1 = −1.65 % (in [1e-4, 3e-2]); the without-arm SAR10g(c\*) sits 0.506 % above 359's. |q_with / q_without| 2.233 (P17), 2.308 (P21), 2.305 (ccw). Rod move `_residual(s_with, s_without)` 3.56e-4 (≤ 1e-3). P_rod / P_phantom 8.31e-5.
- **One printed prediction outside its range (a reading, not a red):** the remesh move `_residual(s_without, s_rec)` read **1.24e-2** against the predicted [1e-5, 1e-2] — I21a's grading around the rod moves the coil's node 2×2 slightly more than this family's prior O(1e-4–1e-2) spread. No band is involved.
- Not claimed: an accurate tip ratio (no published comparand until (7c); `ANS-10` is the accuracy route), any lead other than this σ rod, heating, absolute SAR, degree 2, h-convergence, f₀ / Q with the rod in, AED or Palace.

### §9 item 428 — 2026-10-13 (§10 chain I25; chain step, solve; depends on item 427)

`IMP-3` step 9: the PEC oblate disc in the coil-driven F-small phantom — the 64 MHz four-port fixture, with / without, exactly as step 1 (item 174) did for the rod (implementer; heavy, complex: W1 at `-n 2`, W2 at `-n 4`; ≈ 12 measured / 60 predicted min: ≈ 10 min of windows + 15 fixed + ≈ 35 code; **depends on item 427 — skip if it did not land**). Licence: §10 I25 (the 10-12 weekly, `PROJECT_PLAN.md:4763–4767`), §9 rule (4). A record on a ✅ row; no band, record or Done-when moves.

**Skip check (before anything else).** `tests/mesh/test_birdcage_implant_spheroid.py` on `main` with `SPHEROID_HOLE_CELL_RECORD` registered and `_build_spheroid` defined (item 427's outcome). Absent: skip item 428, leave it ⬜, take the next item (implementer-run.md step 2: never an item whose stated dependency has not landed). Item 427 BLOCKED: the same.

**Fixture.** Step 1's, with `GEO-36` step 7's disc in place of the rod: the four-port F-small birdcage, solid σ = 800 S/m coil (`SIGMA_WIRE_S_PER_M`), 50 Ω lumped sheets, untuned, `STEP2D_FREQUENCY_HZ` = 64e6 (`tests/validation/test_birdcage_power_identity.py:182`), degree 1, phantom σ / εᵣ as step 1. The disc: a = 6 mm, c = 1.2 mm (aspect 5, `IMP-6`'s gated aspect), centre (0.015, 0, 0), symmetry axis x̂, cut as a **hole** (PEC by its exterior wall, tag 402) — `_build_spheroid(**kw)` from item 427's module, the twin of step 1's `_build_rod`. Without the disc: `tests/mesh/test_birdcage_port_sheets._build(True, **kw)`. Both enter `build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=<tuple>, build_only=True)` (`tests/validation/test_port_birdcage_four_port.py:219`). The rim points of interest are `RIM_POINTS` = (0.015, 0, ±0.006): the disc lies in the yz-plane, the P1 field there is ≈ 99 % E_z (η_z = 0.990426 on the rod axis, `20261002T170843Z_IMP-3-step1-W2.log.gz:3864`), i.e. in the disc's plane and normal to the rim at those points — `IMP-6`'s in-plane rim case (1/L_x = 8.0155).

**Mirror symmetry the disc keeps** (as the rod): σ_y (y → −y) kept — it swaps P2 ↔ P4 and fixes P1, P3; σ_z kept — it maps the top rim point onto the bottom one; C4 and σ_x broken. The revolve seam (item 427) breaks σ_y only in where the mesher puts nodes, as any mesh does.

**Code: tests only, no `src/` and no imported module edited.** New module `tests/validation/test_imp3_disc_birdcage.py` with `test_imp3_disc_default_phantom` (W1) and `test_imp3_disc_one_gram_rung` (W2). It imports by name and never a `test_*` symbol:
- from `tests/validation/test_imp3_rod_birdcage.py` (private, rod-agnostic names only): `IMP3_FLAG_OFF_Z_64MHZ` (`:120`), `PEC_WALL_MAX_E` (`:135`), `PORT_IDS`, `DRIVES`, `MIRROR_PAIRS`, `SIGMA_X_PAIRS`, `QDEG_CF`, `CF_CHUNK`, `AXIS_SAMPLES`, `ROD_CELL_GUARD_W2` (320 000), `_say`, `_in`, `_wall_dofs` (`:188`), `_wall_max_e` (`:199`), `_outer_box_tags` (`:205`), `_matrix_quantities` (`:275`), `_pair_spread` (`:467`), `_peak_rss_gib` (`:498`). **Not** `_arm`, `_case`, `_report`, `_assert_identities`, `_control_c1` or `_ball_mean_closed_form`: each hard-codes `TIPS`, `ROD_CENTRE` or the rod's 0.016 m.
- from item 427's module: `SPHEROID_SEMI_MAJOR`, `SPHEROID_SEMI_MINOR`, `SPHEROID_CENTRE`, `RIM_POINTS`, `SPHEROID_HOLE_CELL_RECORD`, `_build_spheroid`, `_inside_true_spheroid`;
- every band, record and helper step 1 imports, from the same files (`test_imp3_rod_birdcage.py:53–109`): `CELL_COUNT_BAND`, `SOLID_CELL_RECORD`, `ONE_GRAM_RUNG_CELL_RECORD` / `PHANTOM_RESOLUTION_1G_RUNG` (`tests/validation/test_birdcage_sar_1g_rung.py`), `REUSE_REPRODUCTION_RTOL`, `RECIPROCITY_BAND`, `PASSIVITY_SIGMA_TOLERANCE`, `DISCRETE_IDENTITY_RTOL`, `ADJACENT_SPREAD_BAND`, `POOLED_SEPARATION_FLOOR`, `C4_COVARIANCE_BAND`, `_solve_driven`, `_exact_shares`, `OUTER_BOX_TAG`, `PHANTOM_RHO_KG_PER_M3`, `ONE_GRAM_KG`, `QUADRATURE_DEGREE`, the quadrature-weight helpers, `PHANTOM_RADIUS` / `PHANTOM_HEIGHT`.
- `_arm_disc(built, label, disc)`: step 1's `_arm` (`:295–411`) copied with `TIPS` → `RIM_POINTS`, and the without-arm background sampled on the disc's in-plane diameter — `AXIS_SAMPLES` = 9 points at (0.015, 0, z), z ∈ [−a, a] (gel in that arm) — returning per drive η_∥ = Σ\|E_z\|²/Σ\|E\|² (world z, toward the rim points), η_t (world y, in-plane tangential there) and η_n (world x, the disc normal). On the disc arm: the wall dofs from the **build tuple's** facet tags, tag 402, max \|E\| per drive.
- `_disc_ball_means(msh, centre, radius)`: step 1's `_ball_mean_closed_form` (`:229–272`) with the local map d = p − `SPHEROID_CENTRE`, p_l = (d_z, d_y, d_x) (local x = world z toward the rim, local z = the normal) and three closed forms from `pec_oblate_spheroid_field` (`src/fem_em_solver/utils/analytical.py:1831`, `semi_major=a`, `semi_minor=c`, `e0=1`, `interior="continue"`): M_∥ (drive `"inplane"` at p_l), M_t (drive `"inplane"` at (y_l, −x_l, z_l), the 90°-rotated in-plane drive) and M_n (drive `"axial"` at p_l), each ⟨\|E_cf\|²⟩/E₀² over the ball's quadrature points on existing (gel) cells. Also return the number of those quadrature points inside the true spheroid (printed).
- **The estimate**, derived in the docstring: for a uniform local background E_bg, Q_cf = η_∥M_∥ + η_t M_t + η_n M_n exactly in the quasi-static limit, because every cross term ⟨e_i·e_j⟩ vanishes by the ball's two mirrors about a rim point (y → −y and the disc's mid-plane). This is step 1's η_z·M_cf + (1 − η_z) with the disc's three-way split.
- Step 1's order otherwise: the without-disc arm first, freed (`del`, `gc.collect()`), then the disc arm; print everything before the first assert; the C1 control in W1.

**Anchors (asserted; every band imported — step 1's set, `test_imp3_rod_birdcage.py:617–636, :702–707, :725–727`):**
- (A0) Census within `CELL_COUNT_BAND`: W1 without = `SOLID_CELL_RECORD` 116 085; W1 with = `SPHEROID_HOLE_CELL_RECORD` (item 427); W2 without = `ONE_GRAM_RUNG_CELL_RECORD` 199 920. W2's disc-arm count is a first record, predicted [235 k, 300 k] (199 920 + item 427's hole increment of [35 k, 85 k]; step 1's rod arm read 250 866). **Guard `ROD_CELL_GUARD_W2` 320 000, checked before any solve.**
- (R0) W1 without-arm Z reproduces `IMP3_FLAG_OFF_Z_64MHZ` entrywise within `REUSE_REPRODUCTION_RTOL` (step 1 measured 2.98e-10, `docs/planning/chunks/IMP-3.md` § item 174 outcome).
- (i) ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND`, both arms. (ii) σ_max(S) ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` and every column Σ\|S\|² ≤ 1, both arms. (P) P_src,exact = P_vol + P_sheet,exact within `DISCRETE_IDENTITY_RTOL` on every drive of both arms.
- (M) σ_y mirror on the disc arm: each of the five `MIRROR_PAIRS` spreads ≤ `ADJACENT_SPREAD_BAND` (rod arm read ≤ 0.0819 %; predicted here [0.005 %, 0.15 %]).
- (E) n × E = 0 on the disc wall: tag-402 dof count (reduced) > 0, max \|E\| ≤ `PEC_WALL_MAX_E` on every drive.
- (S) Separation, without-disc arm: pooled off-diagonal spread ≥ `POOLED_SEPARATION_FLOOR` × the worst class spread — asserted because the same comparison on the same fixture measured 671× (`20261001T171613Z_PORT-21-step5-W2.log:3682`) and step 1's W1 / W2 read 671.1× / 402.7×. Printed on the disc arm.
- Containment: 0.015 + a_1g < R and a + a_1g < H/2 (12.2 mm < 40 mm).

**Negative controls (rule (e)):**
- **(C1) The disc wall left natural** (W1, P1 only, step 1's `_control_c1` pattern, `:639–685`, re-written for `RIM_POINTS`): exterior facets minus 402 tagged `OUTER_BOX_TAG` 499 via `_outer_box_tags`; first assert 499 absent and exterior = outer + 402; a fresh `_solve_driven`. *Asserted:* (P) still closes at `DISCRETE_IDENTITY_RTOL` (a theorem) and the knob is wired — max \|E\| on the 402 dofs > `PEC_WALL_MAX_E` (step 1 read 1.29e-4 on the rod wall). *Predicted, printed:* natural / without rim 1 g in [0.7, 1.3]; PEC / natural in [1.1, 3.0] (a ratio of positive SARs has no ceiling, so no factor is claimed).
- **(C2) What the 4×4 sees of the disc — predicted, printed:** ΔP_phantom/P_phantom (P1) in [0.5 %, 10 %] (the disc's volume is 1.84× the rod's; step 1 read +4.80 % for the rod at P1); ΔΣ_i\|S_i1\|² and max \|ΔS\| beside the rod's −3.55e-5 / 1.32e-4.
- **(C3) Field-level mirror — printed:** rim 1 g pairs P2 / P4, ccw / cw and top / bottom, each predicted ≤ 1 % on W2 (the rod read ≤ 0.69 %); P1 / P3, related only by the broken σ_x, **normalised by the maximum** (step 1's print-definition defect, `IMP-3.md` § item 174 outcome), predicted [50 %, 95 %], ceiling 100 %. `C4_COVARIANCE_BAND` shown as a reference only.

**Printed beside predictions written first — never asserted (the ratio's accuracy is `ANS-10`'s / (7c)'s):**
- **The rim 1 g ratio** Q = SAR_1g,rim(disc) / SAR_1g,rim(without) at both rim points, for P1–P4, ccw, cw (W2 of record; W1's coarse value beside).
- **The closed form from the local background field** (computed here by a stdlib Monte Carlo of `IMP-6`'s closed form over the continuous 6.204 mm ball centred on a rim point, 200 000 gel points; the discrete mask will differ by a few %): M_∥ ≈ 1.650, M_t ≈ 0.880, M_n ≈ 1.068; gel fraction 0.923. Predicted η_∥ ≈ 0.99 for P1 (the rod-axis η_z at the same x), [0.3, 0.6] for P2 / P4; estimate ≈ 1.64 for P1 and ≈ 1.2–1.3 for P2 / P4. **Predicted Q:** P1 in [1.3, 2.6], P2 / P4 in [0.9, 2.0]; **Q / estimate in [0.7, 1.6]** (step 1's rod read 0.855–1.552, `…170843Z_IMP-3-step1-W2.log.gz:3864–3875`; the disc's quasi-static scale (\|k\|a)² = 1.04 % against the rod's 7.37 %). The 1 g ratio is a weak observable of the rim: the pointwise rim value (1/L_x)² = 64.25 is spread over a ball ≈ a in radius.
- **Ball masses:** 1.000 g without; predicted 0.923 g with (the gel fraction).
- **Rim fidelity on each disc mesh:** gel cells within 2ρ_r = 0.48 mm of each rim point and the number whose midpoints lie inside the true spheroid (item 427's (S8) asserts 0 on the default-resolution build; the 1 g rung is new — printed), and the count of closed-form quadrature points inside the body. A non-zero count on W2 means the estimate there is read against a truncated rim; the outcome says so and the estimate is not interpreted.
- **Also:** P_phantom per drive, sheet facet counts, ΔS, cells, times, summed peak RSS, L_x = 0.124758 (`oblate_spheroid_depolarisation`, `analytical.py:1756`).

**Windows.** Bash `timeout: 660000`, foreground, each harness call bare and whole. Basis: step 1's W1 241 s at `-n 2` and W2 260 s at `-n 4` (`20261002T170428Z_IMP-3-step1-W1.log.gz`, `…170843Z_IMP-3-step1-W2.log.gz:4162`); the rod build 42 s mesh; item 427 prices the disc build (its W2a).
- **W0** (smoke, collect, ≈ 10 s): `scripts/testing/run_and_log.sh IMP-3-step9-W0-collect "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 30 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_disc_birdcage.py --collect-only -q'"`
- **W1** (the default phantom, predicted [180, 360] s): `scripts/testing/run_and_log.sh IMP-3-step9-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_imp3_disc_birdcage.py::test_imp3_disc_default_phantom -v -s --tb=short'"`
- **W2** (the 1 g rung, the ratio of record, predicted [220, 500] s), with §5.1 durable capture and the trailing `; exit $rc`: `scripts/testing/run_and_log.sh IMP-3-step9-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_imp3_disc_birdcage.py::test_imp3_disc_one_gram_rung -v -s --tb=short > /workspace/logs/imp3-step9-w2-raw.log 2>&1; rc=\$?; echo [capture] rc=\$rc >> /workspace/logs/imp3-step9-w2-raw.log; cat /workspace/logs/imp3-step9-w2-raw.log; exit \$rc'"`
- Order: skip check → code → W0 → W1 → W2 → commit (module, logs, the record) → journal. No module edit after W1 starts (rule (i)).

**Traps already paid for.**
- **Build order:** without the disc first, then with it, in both windows (`GEO-23`). The disc + 0.0025 m phantom combination has never been built: a FAIL is recorded verbatim, never retried.
- **`built=` ignores mesh-sizing keywords:** pass `phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG` to `_build` / `_build_spheroid` themselves.
- **Tag 402 lives in the build tuple, not the sweep's facet tags** — keep that handle.
- **A held factorisation silently solves a different problem** (`PORT-19`, `WF-12`): the C1 control uses `_solve_driven`'s fresh solver.
- **Point evaluation inside the hole raises** (`post/sar.py:357–361`): the background diameter samples come from the without-disc arm only.
- **The basis follows the comparison:** M_∥ and M_t are the in-plane closed form at p_l and at the 90°-rotated p_l; M_n is the axial one. Copying step 1's prolate call, its `ROD_CENTRE` shift or its (1 − η_z) term is the trap.
- **`interior="continue"` near the focal ring** (item 408): the continued field is singular at the focal ring (a − 0.121 mm); the inside-point count above is how a truncated rim shows.
- Rank-local quantities reduced; owned dofs only; collectives outside `rank == 0`; masks and the closed form in numpy, no UFL `<` or `sqrt` comparison.
- Free the without-disc arm before the disc build in W2.
- `-s`; node IDs only; no pytest `-k`; no pipes; Status 124 is the finding — `pgrep -c python3` first, then the rubric's orphan-`.c` sweep only if it reads 0; never lengthen a timeout.

**Scope:** degree 1, 64 MHz, the untuned solid-coil 4×4 with 50 Ω sheets, one aspect-5 disc, normal x̂, at (0.015, 0, 0). **Not covered:** the accuracy of the rim ratio (`ANS-10` / (7c)); absolute or pointwise SAR, heating, compliance or convergence; a square-edged plate (a wedge singularity the spheroid lacks); the tuned or copper coil; F-human (item 427's kept F-human build is for a later step); other placements or orientations; a thin plate (aspect ≥ 20).

**Status it can move.** The §7 `IMP-3` row stays ✅: its Done-when ("the identities green with the rod in on `GEO-36`'s 167 371-cell record mesh ⇒ 🟡; the same identities green on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture", `PROJECT_PLAN.md:2606`) is the rod's and is not amended (rule (j) does not bind). This file gains a **step 9 record**: the identities on both disc meshes, the rim Q per drive, the three η's, M_∥ / M_t / M_n and Q / estimate; one row clause if the row stays ≤ 2 KB. `IMP-6`'s Not-claimed "a plate in a coil" is **not** lifted by the slot — this is a fixture with no accuracy claim; whether to amend is the review's. §10 I25 is done for the weekly to read.

**Negative result:**
- **(A0) or (R0) misses:** the fixture moved. Stop, BLOCKED.
- **(i), (ii), (P), (M) or (E) red:** park on `attempt/IMP-3-step9-<ts>`, write a known-issues entry and mark item 428 BLOCKED in the same commit (rule (d)). Never widen a band. (M) red with the without-disc pairs green: print the per-sheet facet counts and file the entry on the disc grading perturbing the P2 / P4 sheet triangulations.
- **The C1 knob not wired:** a test bug; fix in-slot and re-run W1.
- **W2 over the guard, a gmsh FAIL, or Status 124:** W1 stands as the record of this step; W2's half is marked BLOCKED on the mesher.
- **A printed reading outside its prediction** (Q, η, Q / estimate, masses): report it; no status effect.
- No `PREDICTED RED` line.

**Every-item rules** (PROJECT_PLAN §9 "Every item"): bare whole harness calls; a lost docker socket stops the slot; no edit after the last window; `timeout: 660000` and `-s`; never delete a log; BLOCKED in the parking commit; a killed window is not a result; two heavy windows never at once; complex mode, `FEM_EM_REQUIRE_COMPLEX=1`, `tests/environment` first; no slot opens `docs/private/` or `aed_results/`.

**Executor spawn prompt (verbatim).** "Execute §9 item 428 — `IMP-3` step 9 (§10 I25) — per `docs/planning/chunks/IMP-3.md` §9 item 428. It depends on item 427: if `tests/mesh/test_birdcage_implant_spheroid.py` with `SPHEROID_HOLE_CELL_RECORD` and `_build_spheroid` is not on main, skip it. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W0, W1, then W2, in order, foreground. No edit to the module after W1 starts. Do not edit `tests/validation/test_imp3_rod_birdcage.py`, item 427's module or any `src/` file; never import a `test_*` name or step 1's `_arm` / `_case` / `_assert_identities` / `_control_c1` / `_ball_mean_closed_form`. Every rank enters every collective. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-3-step9-*`, mark item 428 BLOCKED in the same commit, and report."

**Outcome (2026-10-13 16:30 CDT slot, implementer): GREEN — every asserted anchor and both asserted C1 clauses green on W1 and W2; `IMP-3` stays ✅.** New module `tests/validation/test_imp3_disc_birdcage.py`; no `src/` file, step 1's module or item 427's module edited; no `test_*` name and no step 1 helper on the forbidden list imported. A container-side Monte Carlo of the three closed forms before W0 (200 000 points, unlogged) read M_∥ 1.6573 / M_t 0.8823 / M_n 1.0671, gel fraction 0.9231 — the frame map agrees with the brief's figures. Windows: W0 smoke 4 s, 6 collected `20261013T213336Z_IMP-3-step9-W0-collect.log`; W1 heavy `-n 2` 264 s (13 passed), summed RSS 2.40 GiB `20261013T213344Z_IMP-3-step9-W1.log`; W2 heavy `-n 4` 283 s (13 passed, capture rc 0), summed RSS 5.31 GiB `20261013T213816Z_IMP-3-step9-W2.log`. Measured ≈ 13 min.

**Step 9 record** (the aspect-5 PEC disc, a = 6 mm, c = 1.2 mm, normal x̂ at (0.015, 0, 0), in the untuned solid-coil F-small 4×4 at 64 MHz, degree 1; W1 `…213344Z_…-W1.log:3673–3770`, W2 `…213816Z_…-W2.log:3784–3876`):
- (A0) W1 116 085 / 150 140 = both records exactly; W2 without 199 920 = record; **W2 disc arm 233 782 — first record**, just under the predicted [235 k, 300 k] (a reading outside its prediction, no status effect), guard 320 000. (R0) 2.983e-10 (band 1e-6).
- (i) ≤ 1.965e-15; (ii) σ_max ≤ 0.999753069, every column power ≤ 0.8049; (P) ≤ 7.0e-15 on every drive of every arm; (M) disc-arm σ_y pairs max 0.0425 % (W1) / 0.0441 % (W2), predicted [0.005 %, 0.15 %] IN; (E) 6 669 tag-402 dofs, max \|E\| 0 on every drive; (S) without-arm separation 671.1× (W1) / 402.7× (W2) — step 1's digits; disc arm 770.9× / 458.7× printed.
- (C1, W1) natural wall: (P) 3.154e-15, max \|E\| on the 402 dofs 2.627e-4 (knob wired); natural / without 1.0659 / 1.0610 (pred [0.7, 1.3] IN), PEC / natural 1.5363 / 1.5517 (pred [1.1, 3.0] IN).
- (C2) ΔP_phantom/P at P1: W1 −0.544 % (outside the predicted [0.5 %, 10 %] by sign), **W2 +0.558 %** (IN); ΔΣ\|S_i1\|² +1.41e-5 / +3.55e-5; max \|ΔS\| 6.07e-5 / 2.29e-4 (rod −3.55e-5 / 1.32e-4). The 4×4 barely sees the disc.
- (C3, W2) kept-mirror pairs ≤ 0.232 % on the disc arm, ≤ 0.875 % without (pred ≤ 1 %, IN); P1/P3 normalised by the max 88.7 % (pred [50 %, 95 %] IN). On W1's coarse mesh P2/P4 and top/bottom pairs read up to 9.3 % (disc) and 71.7 % (without, P4): the weak-drive rim SARs (≈ 30× below P1) are mesh-noisy at the default phantom; P1, ccw and cw pairs stay ≤ 0.7 %.
- **Rim 1 g Q (W2 of record, top / bottom):** P1 1.5921 / 1.5934 (pred [1.3, 2.6] IN); P2 1.1369 / 1.1298, P4 1.1479 / 1.1364 (pred [0.9, 2.0] IN); P3 1.5617 / 1.5637; ccw 1.5524 / 1.5531; cw 1.5545 / 1.5552. W1 beside: P1 1.6376 / 1.6463, ccw 1.6080 / 1.5334; W1's P2 / P4 bottom 0.849 / 0.715 are the mesh noise above.
- **η's (W2, in-plane diameter, without arm; ∥ / t / n):** P1 0.998816 / 0.000464 / 0.000720; P2 0.878558 / 0.117120 / 0.004322; P4 0.847404 / 0.145189 / 0.007407; P3 η_∥ 0.989060; ccw 0.992580; cw 0.992122. **η_∥ for P2 / P4 is above the predicted [0.3, 0.6]**: the side-port background at the disc is mostly z-directed too. Printed, no status effect.
- **M's on the W2 disc mesh (top / bottom):** M_∥ 1.656547 / 1.652899, M_t 0.883409 / 0.883000, M_n 1.068570 / 1.068132 (pred 1.650 / 0.880 / 1.068); gel ball volume / (4/3)πa_1g³ 0.9187 / 0.9255; ball masses 0.999884 / 0.999874 g without, 0.923539 / 0.923097 g with (pred 1.000 / 0.923). L_x 0.124758, (\|k\|a)² 1.0405 %.
- **Estimate and Q / estimate (W2):** P1 estimate 1.6558 / 1.6521, **Q / estimate 0.9615 / 0.9644**; P3 0.9465 / 0.9498; ccw 0.9402 / 0.9427; cw 0.9417 / 0.9441; P2 0.7272 / 0.7242; P4 0.7454 / 0.7395 — all inside [0.7, 1.6]. The P2 / P4 estimates (1.56 / 1.54) sit above the predicted [1.2, 1.3] because their η_∥ does. W1's P1 Q / estimate reads 0.9984 / 1.0029.
- **Rim fidelity:** gel cells within 2ρ_r of a rim point 54 / 66 (W1), 62 / 69 (W2); **0 midpoints inside the true spheroid on both disc meshes** (the 1 g rung is new, and it reads 0). The closed-form ball means include 2 388 / 2 367 (W1) and 2 390 / 2 423 (W2) quadrature points inside the true body (gel tets whose flat faces chord the curved wall); none is non-finite. The estimate's comparand there is the continued field, as briefed.
- Not claimed: the accuracy of the rim ratio (`ANS-10` / (7c)); absolute or pointwise SAR; convergence; a plate in a coil (whether to amend `IMP-6`'s Not-claimed is the review's call).

### §9 item 464 — 2026-10-15 (§10 chain I17; chain step, no solve)

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.*

**What it is.** A cap-ball radius knob on `implant_rod`: the key `tip_ball_radius` [m], absent / `None` = today's `2·radius` (2 mm for `ROD`). The I17 builds then show that with the knob at 4 mm the cells the pointwise tip sample is read from sit inside the ball, and print the I30 price. No solve.
- **Licence:** §10 I17 (`PROJECT_PLAN.md:4800–4804`, the 2026-10-14 weekly), a weekly-enumerated step under §9 rule (4). The weekly's ruling on entry 33 (`PROJECT_PLAN.md:4595–4602`) opens it: "I17 opens (the ball-radius knob on `implant_rod`, no solve); no finer cap-ball rung is bought". Known-issues 2026-10-14 "Next" (`docs/testing/known-issues.md:955–961`) agrees. **No band, tolerance, Done-when or target is set or moved.**
- **Executor:** implementer.
- **Class:** chain step, no solve.
- **Tier:** **heavy, not standard as the weekly wrote.** The measured basis does not fit 180 s: one 1 g-rung off build takes 41 s, the 2 mm rod build 62 s and the r/8 build 87 s at `-n 2` (`20261007T170700Z_GEO-36-step4-W2-rungs.log:1917, 3752, 7500`). Four builds with their dof counts took 330 s (`…:7674`). Each I17 window builds off + 2 mm + 4 mm, so ≈ 200 s of builds plus ≈ 60 s of collection and dof spaces. The review may re-tier it.
- **Ranks:** `-n 1` and `-n 2`, as the weekly names.
- **Cost:** ≈ 20 measured / ≈ 50 predicted min (15 fixed + 15 code + ≈ 20 of windows). That is above the chain no-solve median (10 min) because rule (c) owes three re-gate windows (≈ 11 min). See "Problems for the review".
- **Depends on:** nothing. Serial on nothing. Never run concurrently with another heavy window.

**Settled from the code.**
- **The hard-coded 2 mm.** `MeshGenerator._build_birdcage_port_model` (`src/fem_em_solver/io/mesh.py:5086`) sets each cap `Ball` field's `Radius` to `2.0 * rod_r` (`:6598`). Its print says "radius 2r" (`:6602–6608`).
  - `implant_rod_cad["tip_ball_radius_m"] = 2.0 * float(implant_rod["radius"])` (`:6120`) reports it.
  - The ball is centred on each **cap centre** `rod_c ± ½·rod_l·rod_a` (`:6593`), z = ±0.015 m for `ROD`. It is not centred on the tip apex at z = ±0.016 m.
- **The normaliser** `_normalise_implant_rod` (`:4952–5083`) validates every key and returns a **fixed dict** (`:5075–5083`). A new key that is not added to that return dict never reaches the sizing code. `tip_size` (`:5011–5029`, `GEO-36` step 4) is the precedent to copy.
- **Code (`src/`, additive, rule (c)):**
  - In `_normalise_implant_rod`: read `tip_ball_radius = implant_rod.get("tip_ball_radius", None)`.
  - If it is not `None`, refuse it before gmsh with a `ValueError` matching "implant_rod" when it is a bool, not a number, not finite, `≤ radius`, or `> 4·radius`. A ball no bigger than the cap is pointless; the upper bound is a runaway guard. Also refuse it unless both cap centres satisfy `hypot(x, y) + tip_ball_radius < phantom_radius` and `|z| + tip_ball_radius < ½·phantom_height`, following the coat check at `:5064–5074`.
  - Add `"tip_ball_radius": tip_ball_radius` to the return dict.
  - At `:6598`: `ball_r = 2.0 * rod_r if implant_rod.get("tip_ball_radius") is None else float(implant_rod["tip_ball_radius"])`, then `setNumber(ball, "Radius", ball_r)`.
  - At `:6120`, report the same `ball_r`.
  - The print at `:6602–6608` may append `ball radius={ball_r:.6e} m`. It is a print, not a gate; disclose it.
  - Nothing else in `src/` changes.
  - **Default bitwise:** `2.0 * 1e-3 == 0.002` in IEEE double (scaling by 2 is exact), so the unset path and an explicit `0.002` feed gmsh the same double.
- **New test module** `tests/mesh/test_birdcage_implant_rod_ball.py`, real build, no solve.
  - Env knob `FEM_EM_IMP3_I17` ∈ {`run`}. A module fixture builds off → 2 mm → 4 mm, coarse to fine (`GEO-23`). The refusals run ungated.
  - Imported by name, never a `test_*` name and never the `builds` fixture:
    - from `tests/mesh/test_birdcage_implant_rod_tip.py`: `FIXTURE_KW` (`:131`), `ROD_1G_RUNG_CELL_RECORD` (`:90`), `RUNAWAY_CELL_GUARD` (`:98`), `T1` (`:84`), `T1_RADIUS` (`:100`), `_ball_and_t1` (`:158–183`), `_cell_geometry` (`:142`), `_mean_edge` (`:150`), `_dofs` (`:186`), `_summarise` (`:199`), `AIR_CELL_TAG`, `PHANTOM_CELL_TAG`;
    - the A1–A6 helpers and bands, exactly as that module imports them (`:49–71`);
    - `ROD`, `ROD_RADIUS`, `ROD_LENGTH`, `ROD_CENTRE`, `_build_rod` (`tests/mesh/test_birdcage_implant_rod.py:93–97, :125`);
    - `PHANTOM_RADIUS` / `PHANTOM_HEIGHT` (`tests/mesh/test_birdcage_port_tags.py:28–29`);
    - `point_owner_cells` (`fem_em_solver.post.evaluation`).
  - `_summarise` already passes `cad["tip_ball_radius_m"]` to `_ball_and_t1` (`:221–223`). N_ball is therefore counted inside each build's **own** ball radius.
  - The build calls are `_build_rod(implant_rod=dict(ROD, tip_ball_radius=TIP_BALL_2MM), **FIXTURE_KW)` and the same with `TIP_BALL_4MM`. Here `TIP_BALL_2MM = 2.0 * ROD_RADIUS` and `TIP_BALL_4MM = 4.0 * ROD_RADIUS`. `tip_size` stays unset, so VIn = r/4, "the r/4 grain" I30 needs.
- **"Near-T1 edge", defined exactly as the r/8 window printed it.** The definition is `_ball_and_t1` (`test_birdcage_implant_rod_tip.py:158–183`): every owned cell, of any tag, whose **centroid** lies within `T1_RADIUS` = 0.5 mm of `T1` = (0.015023, 0.000017, 0.017531) m (`scripts/probes/imp3_step4_palace_fields.py:43`), reduced over ranks.
  - The **near-T1 edge** h is the mean over those cells of each cell's mean of its six edge lengths. `20261014T070008Z_IMP-3-step7-I15-r8.log:2639` prints "near-T1 (5.0e-04 m) 26 cells, mean edge 5.081837e-04 m".
  - The **near-T1 edge ratio** is h(build A) / h(build B). `…:2640` prints "local edge ratio h(r/(4√2))/h(r/8) 1.2638" = 6.422237e-04 / 5.081837e-04, from `tests/validation/test_imp3_tip_ladder_1g.py:384–388`.
  - The **T1 cell** is `point_owner_cells(msh, T1[None, :], comm)`, with its mean edge allreduced (MAX, 0 elsewhere) (`…:2641`, code `:365–369`).
- **Why the 2 mm ball misses T1: geometry, not a measurement.** \|T1 − top cap centre\| = √(0.023² + 0.017² + 2.531²) mm = **2.5312 mm**. So every near-T1 centroid lies ≥ 2.031 mm from the cap centre, outside a 2 mm ball. With the ball at 4 mm, the farthest near-T1 centroid is ≤ 3.031 mm out, and a vertex sits within about one local edge of its centroid. The `phantom_resolution` = 2.5 mm (`tests/validation/test_birdcage_sar_1g_rung.py:114`) does not reach T1's cells.

**Anchors (asserted; every band imported by name, never restated).**
- **(R-b)** refusals before gmsh: `tip_ball_radius` ∈ {0.0, −1e-3, nan, inf, `ROD_RADIUS` (the boundary), 4.0001e-3, "big", True} each raises `ValueError` matching "implant_rod", and `gmsh.isInitialized()` is False afterwards. This is the (R) pattern at `test_birdcage_implant_rod_tip.py:263–273`.
- **(B0) default unchanged.** On the 2 mm build:
  - `tip_ball_radius_m == 2.0 * ROD_RADIUS` bitwise;
  - **cells == `ROD_1G_RUNG_CELL_RECORD` 252 118 as integers**. The basis is the same build measured at `-n 1` and `-n 8` (`test_birdcage_implant_rod_tip.py:87–89`) and at `-n 2` (`20261007T170700Z_GEO-36-step4-W2-rungs.log:3752`); gmsh builds serially on rank 0.
  - Printed with True / False beside, not asserted: the tag census, tag-402 facets, N_ball and near-T1 against GEO-36 W2's unset line `…:3752–3753` (census {1: 36258, 2: 99975, 3: 109724, …}, 402 facets 4 958, N_ball 11 938, near-T1 (11, 6.983224e-04 m)).
- **(A0)** CAD masses, OCC wall area and wall centroid are bitwise equal between the 2 mm and 4 mm builds (the knob is sizing-only), as `:295–311`.
- **(A1–A6)** per rod build, against the window's own off build, exactly as `:315–364`, at `CAD_ANALYTIC_BAND` / `EXACT` / `TERMINAL_AREA_BAND`.
- **(D)** N1curl d2 = 2E + 2F and d1 = E, exactly, on every build (`:368–378`).
- **Guard:** every build has fewer than `RUNAWAY_CELL_GUARD` 1 000 000 cells, checked before (D).
- **(Kb) the knob reached gmsh**, on the 4 mm build:
  - `tip_ball_radius_m == TIP_BALL_4MM` bitwise;
  - cells > 252 118;
  - N_ball(4 mm) > N_ball(2 mm).
- **(T1in) the near-T1 cells sit inside the 4 mm ball**, on the 4 mm build:
  - `n_t1 > 0`;
  - the largest distance from the **top cap centre** (0.015, 0, 0.015) m of any vertex of any near-T1 cell is ≤ `tip_ball_radius_m`;
  - the same holds for the T1 owner cell's four vertices.
  - Every rank computes this (MAX allreduce) outside any rank branch. This is the weekly's "the near-T1 edge sits inside the 4 mm ball", made strict: whole cells, not centroids.

**Negative controls (rule (e) labels).**
- **(T1out), asserted — the same predicate on the 2 mm build is false for every near-T1 cell.** The count of near-T1 cells with all four vertices inside the 2 mm ball is 0, with `n_t1 > 0`.
  - The backing is a convexity identity, not a measured factor: a cell inside a ball has its centroid inside it, and every near-T1 centroid is ≥ 2.031 mm from the cap centre.
  - `n_t1 > 0` at this build is backed by GEO-36 W2's 11 near-T1 cells at the same build (`…:3753`).
  - The review may relabel it *predicted* if it reads rule (e) as measurement-only.
- **(Kb)** is the knob control. A key that never reached gmsh (dropped by the normaliser's fixed return dict) reads exactly 252 118 and fails by an integer.
- **Predicted, printed beside the measured values, never asserted** (no prior 4 mm measurement exists):
  - 4 mm cells ≈ 350 k, range [320 k, 390 k]. The basis is N_ball 11 938 in the 2 × ≈ 25.8 mm³ of gel inside the 2 mm balls (≈ 232 cells/mm³ at VIn r/4). The 4 mm balls hold 2 × ≈ 254 mm³ of gel, so N_ball ≈ 118 k, range [95 k, 135 k], replacing ≈ 17 k cells.
  - Near-T1 count range [60, 250]; near-T1 mean edge ≈ 3.4e-4 m, range [2.9e-4, 4.0e-4]. The basis is GEO-36's mean ball edge at VIn r/4, 3.397517e-04 m (`…:3753`).
  - Near-T1 edge ratio h(2 mm)/h(4 mm) ≈ 2.05, range [1.7, 2.4].
  - Largest near-T1 vertex distance from the cap centre range [3.1, 3.6] mm.
  - **Predicted: the 4 mm build is over `ROD_CELL_GUARD_W2` = 320 000** (`tests/validation/test_imp3_rod_birdcage.py:115`). See the I30 price below.

**Printed, never gated.**
- Per build: cells, census, E / F / d1 / d2, build seconds, N_ball and mean ball edge, near-T1 count and mean edge, the T1 cell's original index and mean edge, the largest near-T1 / T1-cell vertex distance from the cap centre, and cells against `ROD_CELL_GUARD_W2`.
- The near-T1 edge ratio h(2 mm)/h(4 mm), beside the r/8 window's 1.2638 (`…r8.log:2640`).
- **The I30 price row** (rod arm, degree 2, `-n 8`). The measured basis is item 305 W1: 276 088 cells, d2 1 775 954, build 74.92 s, sweep 267.03 s, window 459 s, 53.85 GiB summed (`20261008T033458Z_IMP-3-item305-W1-d2rod-r4s2.log:2078, 2167, 2726`).
  - The formula at e ∈ {1.23, 1.37}: sweep = 267.03·(d2/1 775 954)^e; window = 459 + (sweep − 267.03) + (build − 74.92); RSS = 53.85·(d2/1 775 954)^e.
  - The decision line: **"I30 slot-feasible at `-n 8`: cells < 320 000 and window ≤ 545 s at e = 1.37"**, yes / no. The 545 s cap leaves 45 s of margin under 590 s, where item 305's ruling (4) kept 20 s on a prediction that landed 36 s inside its range.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; container `timeout -k 30`; `-s`; `tests/environment` first; no `-k`; no pipes; real build, no complex sourcing, as GEO-36 step 4). **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**
- **Pre-check:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0. No cache sweep is owed: there are no forms and no JIT.
- **W0** collect (smoke, ≈ 10 s):
  `scripts/testing/run_and_log.sh IMP-3-I17-W0-collect "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 30 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_rod_ball.py --collect-only -q'"`
- **W1** (`-n 1`, predicted ≈ 260 s, range [200, 380]):
  `scripts/testing/run_and_log.sh IMP-3-I17-W1-n1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_IMP3_I17=run timeout -k 30 420 mpiexec -n 1 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod_ball.py -v -s --tb=short'"`
- **W2** (`-n 2`, the window of record, predicted ≈ 260 s, range [200, 380]):
  `scripts/testing/run_and_log.sh IMP-3-I17-W2-n2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_IMP3_I17=run timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod_ball.py -v -s --tb=short'"`
- **Rule (c) re-gates, in the same slot, verbatim from the step 4 / step 7 precedents:**
  - **R1**, step 4's tip-knob gate (noop), 238 s at `-n 2` (`GEO-36.md:184, 214`). The unset build must equal explicit r/4 and the 252 118 record.
    `scripts/testing/run_and_log.sh IMP-3-I17-R1-noop-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace FEM_EM_GEO36_STEP4=noop timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod_tip.py -v -s --tb=short'"`
  - **R2**, the F-small rod paths, 257 / 254 / 271 s (`GEO-36.md:216, 315, 531`). 167 371 / 182 371 / 116 085 must print unchanged and (A7) must still refuse. A `test_birdcage_volumes_partition_the_box` red carrying the surface 59 / 79 string is the 🔴 OPEN known issue (`GEO-36.md:97`), not this change.
    `scripts/testing/run_and_log.sh IMP-3-I17-R2-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment tests/mesh/test_birdcage_implant_rod.py tests/mesh/test_birdcage_port_tags.py tests/mesh/test_birdcage_conductor_hole.py -v -s --tb=short'"`
  - **R3**, the kept rod under `as_hole`, which shares `_normalise_implant_rod` and the Ball code. 185 s (`GEO-36.md:502, 531`). `F_HUMAN_KEPT_ROD_CELL_RECORD` 436 541 and F1–F9 must be unchanged.
    `scripts/testing/run_and_log.sh IMP-3-I17-R3-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 420 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/mesh/test_birdcage_implant_rod_f_human.py -s -v --tb=short'"`
- **Order:** pre-check → code → W0 → W1 → W2 → R1 → R2 → R3 → commit. **No `src/` or module edit after W1 starts** (rule (i)). An edit forced later re-runs W1 / W2 and every re-gate.
- **Compute:** ≈ 10 + 260 + 260 + 238 + 257 + 185 ≈ 1 210 s ≈ 20 min. Every window is under 10 min, so durable capture is not owed.

**Traps already paid for.**
- **The normaliser's fixed return dict.** A key validated but not returned never reaches `:6598`; (Kb) is the tell.
- **Never import a `test_*` name or the `builds` fixture** from `test_birdcage_implant_rod_tip.py`: pytest would collect it into the new module (item 305's trap, `IMP-3.md:1078`). Importing the module runs its probe load (`:83`), which is harmless.
- **`_ball_and_t1` takes the cell tags `built[1]`** and reduces internally, so call it on every rank. `point_owner_cells` is collective. **No collective inside a `rank == 0` branch** (item 227's `-n 8` deadlock; item 440's W1q, `IMP-6.md:879`).
- **Owned cells only** (`size_local`). The vertex-distance max is reduced with MAX.
- **Build order is off first, then rods coarse to fine** (`GEO-23`). A gmsh FAIL is recorded verbatim and never retried. The 4 mm ball has never been built.
- **A requested ball radius is not the local step at T1.** Print the measured near-T1 edge, and never write "step ratio 2" for it (item 305's trap, `IMP-3.md:1079`).
- `tip_size` is never set in this item. Combining knobs is I30's question, not this one's.
- Never `-k a or b`; `-s` on every window; Status 124 means `pgrep -c python3` first; never lengthen a timeout.
- **Rule (i):** the re-gates run on the `src/` as committed. Never path-checkout a record file from a branch.

**Scope.** Covered: the knob, its refusals, the default-path identity, the CAD identities at both widths, and the near-T1 location at 2 mm and 4 mm on the 1 g-rung hole build at VIn r/4. Not covered:
- any solve or η / Q reading (that is I30);
- convergence; accuracy (`ANS-10`);
- a ball at any radius other than 2 mm and 4 mm;
- the kept or coated rod at 4 mm; F-human at 4 mm;
- any band; `ROD_CELL_GUARD_W2`.

**Status it can move** (rule (j)). §7 `IMP-3` stays ✅. Its Done-when (`PROJECT_PLAN.md:2611`) is the step-1 / 1 g-rung fixture's, met, and not amended.
- **Record clause:** the row gains "I17 (item 464): `tip_ball_radius` knob, default bitwise (252 118); at 4 mm the near-T1 cells lie inside the ball, ⟨cells⟩ cells, I30 ⟨slot-feasible / over the guard → `xl`⟩".
- **Known-issues 2026-10-14** (`known-issues.md:932–961`) gains one "Item 464" line with the 4 mm geometry and the I30 decision. It is **not** retired; I30 retires it or keeps it.
- **§10 I17 is done.** I30's tier is decided by W2's decision line (see item 465).
- `GEO-36`'s row is not edited. The knob is recorded here, as the `tip_size` knob's I15 consequences were.

**Negative result.**
- **W0 red:** a module defect. Fix it in-slot, then re-run W0.
- **(R-b) red:** a refusal path is missing. Fix it in-slot before W1.
- **(B0) red** (the 2 mm build is not 252 118): the default path moved. Stop. Revert the `src/` change if it is the cause; otherwise it is an image move: report and park, BLOCKED (rule (d)).
- **(Kb) red:** the key did not reach gmsh. Fix the return dict in-slot and re-run W1 / W2.
- **(T1in) red at 4 mm**, with (Kb) green: a finding. The 4 mm ball does not contain the near-T1 cells (cells larger than predicted near the ball edge). Report the largest vertex distance. I30 is **not** queueable as posed; the weekly re-poses the radius. No band moves.
- **(T1out) red:** a test bug (the predicate or the reduction). Fix it and re-run.
- **A1–A6 / (D) / any re-gate red:** report every number, park on `attempt/IMP-3-I17-<ts>`, mark item 464 BLOCKED (rule (d)), and open a `GEO-36` known-issues entry. **Never widen a band.**
- **A gmsh FAIL on the 4 mm build:** record it verbatim, BLOCKED on the mesher; never retry.
- **4 mm cells over `ROD_CELL_GUARD_W2`:** *not* a red. It is the predicted reading, and it rules I30 out of the slot tier (item 465).
- No `PREDICTED RED` line.

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

**Executor spawn prompt (verbatim).** "Execute §9 item 464 — `IMP-3` §10 I17 (the `tip_ball_radius` knob) — per `docs/planning/chunks/IMP-3.md` §9 item 464. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. The `src/` change is additive and confined to `_normalise_implant_rod` (validate + return the key) and the cap-Ball radius / CAD report in `_build_birdcage_port_model`; the default must stay bitwise; disclose it. Run W0, W1, W2, R1, R2, R3 in order, foreground; no edit after W1 starts. Never import a `test_*` name or the `builds` fixture. Every collective is called on every rank, outside any rank-0 branch. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-3-I17-*`, mark item 464 BLOCKED in the same commit, and report."

**Problems for the review.**
- (1) The weekly wrote "Standard". The measured builds put each window at ≈ 260 s, so it is drafted heavy at `timeout -k 30 420`.
- (2) Rule (c) owes R1–R3, ≈ 11 min. Dropping R3 is defensible only if the review rules that the F-human kept-rod path never reaches `:6598`. It does reach it: the Ball loop is in the shared `_build_birdcage_port_model`.
- (3) **The 4 mm ball is predicted at ≈ 350 k cells, over `ROD_CELL_GUARD_W2`.** I30 as the weekly posed it ("heavy, `-n 8`") is therefore predicted **not** slot-feasible; see item 465. A 3.5 mm ball would contain the near-T1 cells only narrowly (largest vertex distance predicted ≈ 3.3–3.4 mm) and is predicted at ≈ 305–325 k cells, straddling the guard. That is the weekly's to choose; this item builds only 2 mm and 4 mm.

**Outcome (2026-10-15 07:30 slot, implementer): DONE, first try, no edit after W1.** `src/` change disclosed: `_normalise_implant_rod` validates and returns `tip_ball_radius` (refusals as the brief), and `_build_birdcage_port_model` uses it for the cap-Ball `Radius` and the CAD report (`ball_r`, default `2.0 * rod_r`); the sizing print gains `ball radius=…`. W1 (`-n 1`, 273 s) and W2 (`-n 2`, 268 s) print identical numbers: (R-b) 8/8 refused before gmsh; (B0) 2 mm = **252 118** exact, census / tag 402 4 958 / N_ball 11 938 / near-T1 11 at 6.983224e-04 m all True against GEO-36 W2; (A0) bitwise; (A1–A6) green (A4 0.996254 / 0.996444); (D) exact (d2 1 620 748 / 2 251 660); (Kb) 4 mm **351 478** cells (predicted [320 k, 390 k]), N_ball 111 415 (predicted [95 k, 135 k]); (T1in) **130** near-T1 cells (predicted [60, 250]), mean edge 3.383011e-04 m, largest vertex distance **3.151970e-03 m** ≤ 4 mm (predicted [3.1, 3.6] mm), T1 cell (original 227010) 2.768817e-03 m / edge 3.282035e-04 m; (T1out) 0 of 11 inside 2 mm (nearest centroid 2.065680e-03 m); edge ratio h(2 mm)/h(4 mm) **2.0642** (predicted [1.7, 2.4]). **I30 decision line: no** — 351 478 ≥ 320 000 and window 564.7 s at e = 1.37 (552.6 s at 1.23; RSS 72.10 / 74.54 GiB) → I30 goes to the weekly as an `xl` candidate (item 465 stays unqueued). Re-gates green: R1 252 118 unset = explicit r/4 (247 s); R2 39 passed, 167 371 / 182 371 / 116 085, (A7) refusing (263 s); R3 14 passed, 436 541 (173 s). Logs: `20261015T123211Z_IMP-3-I17-W0-collect.log`, `20261015T123222Z_IMP-3-I17-W1-n1.log`, `20261015T123700Z_IMP-3-I17-W2-n2.log`, `20261015T124134Z_IMP-3-I17-R1-noop-regate.log`, `20261015T124547Z_IMP-3-I17-R2-regate.log`, `20261015T125015Z_IMP-3-I17-R3-regate.log`.

### §9 item 465 — 2026-10-15 (§10 chain I30, `IMP-3` step 10; chain step, solve) — CONDITIONAL, predicted NOT slot-feasible

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.*

**Queue verdict (read first).** **Do not queue this item now.** Item 464 must land first. I30 is slot-feasible only if item 464's W2 decision line reads **yes**: the 4 mm build has fewer than `ROD_CELL_GUARD_W2` = 320 000 cells **and** its scaled window is ≤ 545 s at `-n 8`. The prediction is ≈ 350 k cells (item 464, "Negative controls"). That is the size of the r/8 rung (342 072), which GEO-36 W2 priced over the heavy ceiling at `-n 8`: 602.8–621.0 s (`20261007T170700Z_GEO-36-step4-W2-rungs.log:7524–7525`).
- Item 327's ruling (5) read that rung as needing `xl` "for the 590 s ceiling, not for memory" (`IMP-3.md:1180`). It also trips `_arm_d`'s (G1) assert before any solve (`test_imp3_tip_ratio_1g.py:165`, `ROD_CELL_GUARD_W2` 320 000). No slot is licensed to raise that guard; only `xl` entry 33's price paragraph did, for its own window (`IMP-3.md:1126, 1180`).
- **If item 464 reads "no",** I30 goes to the weekly as an `xl` candidate. Entry 33's shape is the precedent: 342 072 cells, `-n 16`, 511 s, 80.55 GiB, under a disclosed `ROD_CELL_GUARD_XL` patch (`known-issues.md:936–937`; `…r8.log:2663`). The daily does not commission it (§5.1).
- **If item 464 reads "yes",** the review queues this brief as written.

**What it is.** I14's rod arm at degree 2 on the exact-rectangle 1 g-rung fixture, built with `implant_rod = dict(ROD, tip_ball_radius=TIP_BALL_4MM)`: VIn r/4 ("the r/4 grain") inside a 4 mm cap ball. η_rod(T1) is printed beside the three-rung record 3.482328 / 3.556262 / 3.530461.
- **Licence:** §10 I30 (`PROJECT_PLAN.md:4805–4811`), weekly-enumerated, §9 rule (4), serial on I17.
- **Executor:** implementer.
- **Class:** chain step, solve.
- **Tier:** heavy, `-n 8`, complex.
- **Cost:** ≈ 12 measured / ≈ 45 predicted min.
- **Depends on item 464 landing with its decision line "yes". If it did not land, or read "no", skip and leave this item ⬜** (implementer-run.md step 2).
- **No band, tolerance, Done-when or target moves.**

**Code (test module only; no `src/`).**
- Add `b4` to `KNOBS` in `tests/validation/test_imp3_tip_ladder_1g.py` (`:84`) and a new `_b4()` beside `_r4s2()` / `_r8()`. **Those two bodies are not touched**; duplicate rather than refactor, as item 327 did (`IMP-3.md:1135`). The `git diff` may delete only the docstring knob sentence, the `KNOBS` line and the dispatch lines.
- `_b4()`:
  - reads `ours_imp3g_d2off.json`, `ours_imp3g_d2rod.json`, `ours_imp3h_r4s2.json` and `ours_imp3h_r8.json` (all present in `logs/imp3/`);
  - builds `_build_rod(implant_rod=dict(ROD, tip_ball_radius=TIP_BALL_4MM), **FIXTURE_KW)`, with `TIP_BALL_4MM` imported from item 464's module;
  - calls `a = _arm_d(built, "degree 2 rod b4", True, 2)` **unpatched**;
  - writes `logs/imp3/ours_imp3h_b4.json` with mode `"x"`.
- Record constant: `I17_B4_CELLS` = item 464's W2 4 mm cell count, cited with its log line.

**Anchors (asserted; bands imported through I14's helpers, never restated).**
- (G1) cells < `ROD_CELL_GUARD_W2`, unpatched, inside `_arm_d`.
- (G2), (i), (ii), (PP), (M), (E), (V), (D) via `_identities(a, True)`.
- (G3) tag-402 facets > 0 and all exterior.
- (Kb4) the knob and the mesh:
  - `tip_ball_radius_m == TIP_BALL_4MM` bitwise and `tip_size_m == ROD_RADIUS / 4` bitwise;
  - cells == `I17_B4_CELLS` as integers. The backing is item 464's `-n 1` / `-n 2` builds of the same call; gmsh builds serially on rank 0, and 252 118 read the same at `-n 1/2/8`.
- (T1in) item 464's predicate re-asserted on this mesh.
- (Rb) "3.482328" recomputed from I14's JSONs, string-equal, sizes as integers (item 305's (Rb)).
- (Rb8) "3.556262" from `ours_imp3h_r4s2.json`, as item 327.
- (Rb16) "3.530461" from `ours_imp3h_r8.json` + off, string-equal (`…r8.log:2643`); sizes 342 072 / 2 200 500.

**Negative controls (rule (e)).**
- (Kb4) and (T1in) are *asserted*, backed by item 464's measurement on the same build.
- (D) is the order control: degree 1 would fail it by ≈ 5.4×.
- Predicted, printed: |η_rod^{b4}(T1)/η_rod^{r/4}(T1) − 1| ≥ 1e-3. The basis is the T1 cell's mean edge falling from ≈ 0.77–0.79 mm (`…r8.log:2641`; item 305 `…:2140`) to ≈ 0.3–0.4 mm.

**The readout (printed, never asserted; the weekly's rule, `PROJECT_PLAN.md:4806–4811`).**
- η_rod(T1), P1, at b4, beside r/4 3.482328, r/(4√2) 3.556262 and r/8 3.530461.
- **Inside [3.48, 3.56]:** the ladder's spread reads as the ball-edge effect, and known-issues 2026-10-14 retires (`known-issues.md:958–959`).
- **Outside:** the pointwise tip ratio is recorded as unresolved at this grain, `IMP-3`'s row says so, and the entry stays as the record.
- Predicted: inside, in [3.48, 3.60].
- Also printed:
  - η at all six points × four drives;
  - the F1 / F2 remesh floors against r/4;
  - Q per tip and drive against r/4 (predicted \|ΔQ/Q\| ≤ 0.5 %);
  - near-T1 and T1-cell edges against item 464's W2;
  - the price against item 464's prediction.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; `-k 30`; `-s`; `tests/environment` first; `FEM_EM_REQUIRE_COMPLEX=1`; no `-k`; no pipes). **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**
- **Pre-check:**
  - `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0;
  - `ls logs/imp3/` shows all four input JSONs and no `ours_imp3h_b4.json`;
  - `git diff 0f6bbbc7 -- tests/validation/test_imp3_tip_ratio_1g.py` is empty.
  - Then the rubric's cache sweep, label `IMP-3-I30-cache-sweep` (93 s, item 305).
- **W0** import (`-n 1`, knob unset, ≈ 70 s, `IMP-3.md:1110`):
  `scripts/testing/run_and_log.sh IMP-3-I30-W0-import "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short'"`
- **W1** (`b4`, `-n 8`, the closing window). The range comes from item 464's decision line, ≤ 545 s by construction. It uses §5.1 durable capture verbatim, trailing `; exit $rc` included:
  `scripts/testing/run_and_log.sh IMP-3-I30-W1-d2rod-b4 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP7=b4 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short > /workspace/logs/imp3-I30-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/imp3-I30-w1-raw.log; cat /workspace/logs/imp3-I30-w1-raw.log; exit \$rc'"`
- **Status 124:** `pgrep -c python3` first. The killed window is the price record. Park and mark BLOCKED (`xl`, the weekly's). **Never lengthen the timeout, raise `-n` or patch the guard.**

**Traps.**
- Item 305's list carries in full (`IMP-3.md:1061–1082`).
- (G1) is **unpatched** on this path. Copying `_r8`'s monkeypatch is the trap: that licence was `xl` entry 33's alone.
- The mode-`"x"` writer means there is no second W1.

**Scope.** One radius (4 mm) at one VIn (r/4), degree 2, 64 MHz, P1 of record. Not covered: convergence or an observed order; accuracy (`ANS-10`); a bracket; other radii; F-human.

**Status it can move.**
- §7 `IMP-3` stays ✅; the Done-when is met and not amended.
- The row's step-7 sentence gains the I30 reading, inside or outside the band.
- Known-issues 2026-10-14 retires if the reading is inside [3.48, 3.56]; otherwise it stays as the record.
- No band moves.

**Negative result.**
- An anchor red: report every number, park on `attempt/IMP-3-I30-<ts>`, mark BLOCKED (rule (d)) and open a known-issues entry; never widen.
- (Rb*) red: no η printed as a reading.
- A reading outside its prediction: report it, no status effect beyond the row's record.
- No `PREDICTED RED` line.

**Executor spawn prompt (verbatim).** "Execute §9 item 465 — `IMP-3` §10 I30 — per `docs/planning/chunks/IMP-3.md` §9 item 465, only if item 464 landed and its W2 decision line reads 'I30 slot-feasible at -n 8: yes'; otherwise skip and leave it unchecked. Read the brief and any review rulings beneath it first. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Do not monkeypatch ROD_CELL_GUARD_W2. No edit after W1 starts. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-3-I30-*`, mark item 465 BLOCKED in the same commit, and report."

### Review ruling 2026-10-15 03:00 (items 464, 465)

- **Item 464 accepted, heavy** (each window builds off / 2 mm / 4 mm, ≈ 260 s; the ruling-(c) re-gates of the three builder gates, the F-human kept-rod path included, ≈ 11 min). The 2 mm control is **asserted as a CAD identity** (T1 lies 2.531 mm from the cap centre, so no near-T1 cell can sit inside a 2 mm ball) — an exact geometric fact, not a separation factor, so rule (e) does not apply. The default build must reproduce 252 118 cells exactly.
- **Item 465 (I30) is not queued.** Its 4 mm ball is predicted at ≈ 350 k cells ([320 k, 390 k]) against `ROD_CELL_GUARD_W2` = 320 000 (`tests/validation/test_imp3_rod_birdcage.py:115`, checked), which `_arm_d` asserts before any solve, and the r/8 rung at 342 072 cells priced 603–621 s at `-n 8`, over the 590 s ceiling. Item 464 prints the decision line; **if it reads "slot-feasible: yes", the 10:30 review queues item 465 as drafted; otherwise I30 returns to the 10-17 weekly as an `xl` candidate in entry 33's shape** (a smaller ball is the weekly's choice, not a daily's). The number stays reserved.

## §7 row moved verbatim 2026-10-17 (rotate_plan_archive.py chunks)

 A rod in the coil-driven phantom: with / without on the loaded F-small birdcage at 64 MHz (Phase 7 subgoal (7b)) (opened 2026-09-30 by the weekly; §10 chain … | ✅ Audit PASS 2026-10-02 18:00 review: the rod's identities are green on the record mesh and the 1 g rung, tip Q printed (4.50 at P1), a fixture with no accuracy claim. Steps 3 and 4 (items 235 / 242) are Palace records: ours = Palace to ≤ 2.2e-10 per class on all three arms and to ≤ 1.5e-8 on near-tip E at matched order; the tip ratio's accuracy is `ANS-10`'s, re-homed to Palace. Step 5 (item 252, 2026-10-05): at matched order 2 the codes still agree on near-tip E, and η_rod(T1) moves −4.961 % from order 1 to 2 in both — not p-stable at 1 %; steps 6 / 7 (§10 I14 / I15, 2026-10-05 weekly) carry the question on item 235's exact-rectangle fixture, now the fixture of record. **Step 7 / I15 record (parts 1–3; `xl` entry 33, 2026-10-14):** η_rod(T1) at degree 2 on the cap-ball rungs r/4 / r/(4√2) / r/8 reads 3.482328 / 3.556262 / 3.530461 (Δ₁ +0.0739, Δ₂ −0.0258, opposite signs) with I14's identities and the 1 g ball ratio green on every rung (`20261014T070008Z_IMP-3-step7-I15-r8.log`); entry 33's pre-registered rule read its "otherwise" branch (2026-10-14 weekly): **the pointwise tip ratio is not resolved by the cap ball alone** — the near-T1 edge sits outside the 2 mm ball (known-issues 2026-10-14); no bracket, no finer cap-ball rung; §10 I17 (the ball-radius knob) and I30 (the ratio with the edge inside the ball) follow. I17 (item 464, 2026-10-15): `tip_ball_radius` knob, default bitwise (252 118); at 4 mm the near-T1 cells lie inside the ball (130 cells, every vertex ≤ 3.151970e-03 m from the cap centre; none wholly inside 2 mm), 351 478 cells, I30 over the guard → `xl` (`20261015T123700Z_IMP-3-I17-W2-n2.log`). **Done-when:** the identities green with the rod in on `GEO-36`'s 167 371-cell record mesh ⇒ 🟡; the same identities green on the 1 g rung with the tip ratio printed there ⇒ ✅ as a fixture. The ratio's accuracy is `ANS-10`'s (commissioned by the weekly that reads this). No heating or compliance claim. *History: `docs/planning/chunks/IMP-3.md`.* 

### §9 item 532 — 2026-10-18 (§10 I30a, `IMP-3` step 10a — the `b4` knob for `xl` entry 43; chain no-solve)

*Drafted by a read-only general-purpose agent for the 2026-10-18 03:00 daily review; the review reviews it.*

**What it is.** Add the selectors `b4` and `b4build` to `FEM_EM_IMP3_STEP7` in `tests/validation/test_imp3_tip_ladder_1g.py`, so that `xl` entry 43 (`docs/testing/xl-pending.md`, "### 43.") has a module path to run. Today `KNOBS = ("r4s2", "r8", "r8build")` (`:84`) and the test skips on anything else (`:458–459`), so entry 43's `FEM_EM_IMP3_STEP7=b4` would read Status 0 with one skip. `b4` is the full path: build, solve, readout. It is written here and **run only by entry 43**. `b4build` builds the 4 mm mesh and stops. It is this item's window of record. **No solve.**

- **Licence:** §9 rule (6) (`PROJECT_PLAN.md:3868–3877`). This re-poses a blocked step of live chain I toward the weekly's pre-registered target: §10 I30 = `xl` entry 43 (`PROJECT_PLAN.md:4676–4678`, `:4834`). Entry 43's own text licenses the guard patch: "the guard is patched to the `xl` value on this path only, as entry 33 did, and disclosed in the log" (`xl-pending.md`, entry 43 Status). Item 327's review ruling (1) (`IMP-3.md:1180`) is the precedent; it covered the build-only twin (`r8build`) in-slot. **No band, tolerance, Done-when or target is set or moved.** `ROD_CELL_GUARD_W2` (320 000) and `ROD_CELL_GUARD_XL` (400 000) keep their values and their definitions.
- **Supersedes** item 465's code plan (`IMP-3.md:1723–1730`) on exactly one point. Item 465 called `_arm_d` **unpatched**, because it was written for an `-n 8` slot. On the `xl` path the 351 478-cell mesh is over 320 000, so the unpatched (G1) would trip before the solve. Everything else in item 465's anchors and readout carries over, but item 465 itself stays unqueued.
- **Executor:** implementer.
- **Class:** chain no-solve, lane 1.
- **Family:** I30, attempt 1 of 4.
- **Depends on:** nothing in §9. Inputs on disk, all present 2026-10-18: `logs/imp3/ours_imp3g_d2off.json`, `ours_imp3g_d2rod.json`, `ours_imp3h_r4s2.json`, `ours_imp3h_r8.json` (342 072 / 2 200 500). `ours_imp3h_b4.json` is absent.
- **Tier per window** (all measured on the same family):
  - **W0** import, `-n 1`, knob unset, `timeout -k 30 120`. Measured 5 s on a warm cache (`20261008T123440Z_IMP-3-item327-W0-import.log`, footer) and 71 s cold (item 305, `IMP-3.md:1110`).
  - **W1** `b4build`, `-n 2`, heavy by ceiling (`timeout -k 30 300`). The basis is the same 4 mm build at `-n 2`: 77.98 s (mesh 72.79 s), `20261015T123700Z_IMP-3-I17-W2-n2.log:5648`. The build-only twin `r8build` ran 88 s harness with an 83.57 s build, body 84.3 s and 0.87 GiB summed (`20261008T123452Z_IMP-3-item327-W1b-r8build.log:1961`, footer `:1985`). **Predicted ≈ 90 s, range [70, 200].**
  - **W2** `r8build` re-gate, `-n 2`, heavy by ceiling (`timeout -k 30 300`). It is item 327's W1b verbatim except for the label; 88 s measured (same log). **Predicted ≈ 90 s, range [70, 200].**
  - Compute is ≈ 3–5 min. Every window is under 10 min, so durable capture is not owed.
- **Cost:** ≈ 8 measured / ≈ 20 predicted slot-minutes.
  - The measured figure is the twin, item 327: slot start 07:30 CDT to commit `6c65dee5` at 07:37 CDT, with the same shape (one build-only window plus W0).
  - The prediction adds the re-gate window, a larger body (three JSON reads and a dry-run readout) and the entry-43 edit. Class medians: chain no-solve 18.6, lane 1 ≈ 3.5.
  - **If W2 has not started by minute 40, park** (rule (d)).

**What the module is today** (`git diff 6c65dee5 -- tests/validation/test_imp3_tip_ladder_1g.py` is empty, checked 2026-10-18).
- **Dispatch, `:452–464`:**
  - (N2) asserts `_ratio.ROD_CELL_GUARD_W2 == ROD_CELL_GUARD_W2 == 320_000` before anything else (`:456–457`);
  - an unknown knob skips (`:458–459`);
  - `r4s2` returns before the patch (`:460–462`);
  - `monkeypatch.setattr(_ratio, "ROD_CELL_GUARD_W2", ROD_CELL_GUARD_XL)` (`:463`) then covers every later path;
  - `_r8(build_only=KNOB == "r8build")` (`:464`).
- **The guard.** `_arm_d` asserts `ncells < ROD_CELL_GUARD_W2` when `rod` is true, after `build_four_port_sweep(build_only=True)` and before the solve (`test_imp3_tip_ratio_1g.py:165`). It resolves the name at call time in `_ratio`'s namespace; the value is imported from `test_imp3_rod_birdcage.py:115`.
- **`ROD_CELL_GUARD_XL = 400_000`** (`:117`, comment `:109–116`). At 351 478 cells the 4 mm build sits at 0.879 of it.
- **`_r8(build_only)` (`:282–449`)** is the shape to duplicate:
  - mode-`"x"` pre-check only when not `build_only` (`:287–289`);
  - JSON reads and (Rb8) before the build (`:290–309`);
  - (Kc8), (G1x) and (N1) with the DISCLOSURE print (`:318–340`);
  - mesh anchors stop the solve path (`:341–346`);
  - `build_only` never calls `_arm_d` (`:348–352`);
  - collective mesh measurements on every rank (`:361–369`);
  - `_ratio._edge_face_counts` for 2E + 2F at build-only (`:370–374`).
- **Item 464's helpers** (`tests/mesh/test_birdcage_implant_rod_ball.py`; `git diff fee90cf8 --` on it and on `test_birdcage_implant_rod_tip.py` / `test_birdcage_implant_rod.py` is empty, checked):
  - `TIP_BALL_4MM = 4.0 * ROD_RADIUS` (`:73`);
  - `TOP_CAP` (`:75–77`);
  - `_t1_geometry(msh, comm)` (`:108–145`). This is the (T1in) predicate as item 464 asserted it (`:332–335`): `n_t1 > 0`, `t1_vmax <= r`, `owner >= 0`, `owner_vmax <= r`. It is collective.
- **`tip_size_m` default** is `float(implant_rod["radius"]) / 4.0` (`src/fem_em_solver/io/mesh.py:6292–6296`). This equals `TIP_R4 = ROD_RADIUS / 4.0` (`tests/mesh/test_birdcage_implant_rod_tip.py:112`) bitwise.

**Code. Edit only `tests/validation/test_imp3_tip_ladder_1g.py`, and only additively.** Do not edit `src/`, `test_imp3_tip_ratio_1g.py`, `test_imp3_rod_birdcage.py`, `test_birdcage_implant_rod_ball.py`, `test_birdcage_implant_rod_tip.py` or `test_birdcage_implant_rod.py`.
- **What `git diff` may delete:** only the docstring knob sentence (`:7–8`) and the `KNOBS` line (`:84`). Everything else is insertions. **`_r4s2()` and `_r8()` are untouched: duplicate, never refactor** (item 327's rule, `IMP-3.md:1135`).
- **`KNOBS = ("r4s2", "r8", "r8build", "b4", "b4build")`.** The docstring gains one paragraph naming item 532, the two paths and the patch.
- **Imports, insertions only:**
  - `TIP_R4` joins the `test_birdcage_implant_rod_tip` import block;
  - a new line `from tests.mesh.test_birdcage_implant_rod_ball import TIP_BALL_4MM, TOP_CAP, _t1_geometry`.
  - **Never import `builds` or any `test_*` name from it.** Pytest would collect it (`IMP-3.md:1078`).
  - A module-level `assert TIP_BALL_4MM == 4.0 * ROD_RADIUS` is allowed as an import-time identity.
- **Dispatch.** Insert two lines between `:463` and `:464`, so `b4` / `b4build` run **inside the existing patch** and `r8` / `r8build` are byte-for-byte unchanged:
  `if KNOB in ("b4", "b4build"):` / `_b4(build_only=KNOB == "b4build"); return` (two lines, `return` on its own line).
- **New record constants**, each with its log-line citation in a comment; printed comparands except where asserted:
  - `I17_B4_CELLS = 351_478` (asserted). From `20261015T123700Z_IMP-3-I17-W2-n2.log:5648`, and identical at `-n 1` (`…I17-W1-n1.log`, item 464 Outcome, `IMP-3.md:1703`).
  - `I17_B4_D2_DOFS = 2_251_660`, E 418 532, F 707 298 (`:5648`, `:5661`).
  - `I17_B4_N402 = 5_366` (`:5649`).
  - `I17_B4_N_BALL, I17_B4_BALL_EDGE = 111_415, 3.489008e-04` (`:5649`).
  - `I17_B4_NEAR_T1 = (130, 3.383011e-04)` (`:5649`).
  - `I17_B4_T1_VMAX = 3.151970e-03` and the T1 cell `(227_010, 3.282035e-04, 2.768817e-03)` (original index, mean edge, largest vertex distance; `:5650`).
  - `I17_B4_BUILD_S = 77.98` (`:5648`).
  - `I17_2MM_T1_VMAX = 3.068100e-03` (`:3755`), for the control below.
  - `R8_ETA_T1_P1 = "3.530461"` (`20261014T070008Z_IMP-3-step7-I15-r8.log:2643`). The r/8 JSON sizes are asserted against the existing `GEO36_R8_CELLS` / `GEO36_R8_D2_DOFS` (`:121–122`), the same integers.
  - `ENTRY33_PRICE` = "build 93.47 s, sweep 278.25 s, SARs 24.57 s, body 402.2 s, 82.47 GiB summed at -n 16 (`…r8.log:2663`); window 511 s, memory.peak 86 486 282 240 B = 80.55 GiB (`:3773`, `:3778`)".
  - `ENTRY43_PRICE_PRED = "500–700 s, 75–100 GiB"`.
  - `ENTRY43_READOUT_INTERVAL = (3.48, 3.56)`. This is entry 43's pre-registered decision interval, **printed IN / OUT, never asserted, never called a band**.
- **New `_b4(build_only)`** beside `_r8`, in this order:
  1. If not `build_only`, assert `ours_imp3h_b4.json` does not exist (mode `"x"`: fail before the build).
  2. Read the four JSONs (off, r/4 rod, r/(4√2), r/8). Assert **(Rb)** "3.482328" with I14's sizes, **(Rb8)** "3.556262" with item 305's sizes, and **(Rb-r8)** "3.530461" with sizes 342 072 / 2 200 500, each string-equal at 6 d.p. and with integer sizes. These run on **both** paths, so W1 proves the reading wiring before the `xl` window. Item 465 called the r/8 check "(Rb16)"; that was a misnomer, so it is renamed here.
  3. `built = _build_rod(implant_rod=dict(ROD, tip_ball_radius=TIP_BALL_4MM), **FIXTURE_KW)`. `tip_size` stays unset, which gives the r/4 grain. Then `cad = built[3]["implant_rod_cad"]` and `cells = int(msh.topology.index_map(3).size_global)`.
  4. **Mesh anchors** (below), the DISCLOSURE line and the (G1x) / (N1) prints copied from `_r8:329–336` with the labels changed. On the `b4` path, any mesh red asserts before `_arm_d`, as `_r8:341–346` does.
  5. If not `build_only`: `a = _arm_d(built, "degree 2 rod b4", True, 2)`, `fails += _identities(a, True)`, `msh = a["case"]["mesh"]`.
  6. Every rank, outside any rank branch:
     - (G3) via `_outer_box_tags(msh, built[2], comm)`;
     - `_ball_and_t1(msh, built[1], cad["tip_ball_radius_m"], comm)`;
     - `_t1_geometry(msh, comm)`;
     - at build-only, `_ratio._edge_face_counts(msh)`, with d2 = 2E + 2F printed beside 2 251 660;
     - `_peak_rss_gib()`.
  7. Print the mesh identity against item 464's W2 constants, each with True / False.
  8. **On `b4` only: write `logs/imp3/ours_imp3h_b4.json` (mode `"x"`, rank 0, then `Barrier`) before the readout prints.** This is a change of order from `_r8`, so that a readout typo after a ≈ 500 s `xl` solve cannot cost the data. Payload: `_ser(a)` plus `rss_gib`, `ranks`, `tip_size_m`, `tip_ball_radius_m`, `near_t1`, `t1_cell`, `n_ball`, `t1_vmax`, `guard_patched`.
  9. **Readout, through a new helper `_b4_readout(eta_new, sar_new, mass_new, label)`**, printed and never asserted:
     - η_rod at six points × four drives as r/4 / r/(4√2) / r/8 / b4;
     - η_rod(T1), P1, beside 3.482328 / 3.556262 / 3.530461, with the move against r/4 and `ENTRY43_READOUT_INTERVAL` IN / OUT;
     - the F1 / F2 remesh floors against r/4;
     - Q per tip and drive at the four meshes;
     - the 1 g masses;
     - the near-T1 edge ratio h(r/4 unset, 6.983224e-04 = `GEO36_UNSET_NEAR_T1_EDGE`, `:102`) / h(b4).

     **On `b4build`, call it once with the r/8 JSON as a stand-in**, every line prefixed `DRY-RUN (r/8 JSON stand-in — not a reading)`. That executes every readout line in-slot, so the `xl` window cannot be the first run of that code. On `b4`, call it with `a`.
  10. Price. On `b4`: cells, E, F, dofs, build, sweep, SARs, arm, body, summed RSS and `-n`, beside `ENTRY33_PRICE` and `ENTRY43_PRICE_PRED`. On `b4build`: build, body, RSS and `-n`, beside `I17_B4_BUILD_S`.
  11. `_say(anchors failed (all) …)`, `del`, `gc.collect()`, `assert not fails`.

**Anchors (asserted on both paths; on `b4` they run before the solve).**
- **(Kb4) the knob and the grain reached gmsh:**
  - `cad["tip_ball_radius_m"] == TIP_BALL_4MM` bitwise;
  - `cad["tip_size_m"] == TIP_R4` bitwise (= `ROD_RADIUS / 4`, the r/4 grain);
  - **cells == `I17_B4_CELLS` (351 478) as integers.**
  - The cell-count backing: item 464's `-n 1` and `-n 2` builds of the same call (`…I17-W2-n2.log:5648`; W1 identical). gmsh builds serially on rank 0. Mode independence (item 464 ran real, this runs complex) is measured twice in this family: 276 088 real (GEO-36 W2) = complex (item 305 W1, `IMP-3.md:1112`), and 342 072 real = complex (item 327 W1b, `IMP-3.md:1183–1184`).
- **(T1in)** item 464's predicate on this mesh via `_t1_geometry`: `n_t1 > 0`, `t1_vmax <= tip_ball_radius_m`, `owner >= 0`, `owner_vmax <= tip_ball_radius_m`. Printed beside 130 / 3.151970e-03 / 227 010 / 2.768817e-03 (`:5650`, `:5665`).
- **(G1x)** `_ratio.ROD_CELL_GUARD_W2 == ROD_CELL_GUARD_XL` and cells < it (0.879 predicted).
- **(G3)** tag-402 facets > 0, all exterior. Printed beside 5 366.
- **(Rb), (Rb8), (Rb-r8)** as in step 2 above.
- **On `b4` only (entry 43's window, not run here):** (G1) inside `_arm_d` at the patched value; (G2), (i), (ii), (PP), (M), (E), (V), (D) via `_identities(a, True)` at the bands I14's helpers import. None of these is restated.

**Negative controls (rule (e)).**
- **(N1), asserted:** cells > the unpatched `ROD_CELL_GUARD_W2` 320 000, so the unpatched (G1) would trip. Ratio 1.098. Backed by item 464's measurement of the same build: "cells 351478 vs ROD_CELL_GUARD_W2 320000 (over)" (`…I17-W2-n2.log:5650`).
- **(Kb4-blind), asserted through (Kb4):** a selector that dropped the key (`dict(ROD)`) builds the 2 mm mesh at **252 118** cells with `tip_ball_radius_m` 0.002 (`…I17-W2-n2.log:3753–3754`). It fails the integer equality by 99 360 cells and the bitwise radius. Measured on the same fixture.
- **(T1in-blind), asserted through (T1in):** on that 2 mm mesh the largest near-T1 vertex distance is 3.068100e-03 m > 2 mm (`:3755`), so the predicate fails against its own radius. Same fixture, same predicate.
- **(N2), asserted (existing, `:456–457`):** the unpatched guard reads 320 000 before every dispatch, so `r4s2` and the unset path never run inside the patch.
- **(N3), predicted:** unset knob ⇒ skip (W0).
- **Predicted, printed only:**
  - W1 summed RSS ≤ 1.5 GiB at `-n 2` (r8build 0.87 GiB);
  - W1 near-T1 / N_ball / T1-cell / 402 values equal to item 464 W2's (True on every line);
  - the dry-run's η(T1) reads 3.530461, the stand-in's own value.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; container `timeout -k 30`; `-s`; complex build; `FEM_EM_REQUIRE_COMPLEX=1`; only `tests/environment/test_complex_mode.py` first, as item 327; no `-k`; no pipes). **The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**
- **Pre-check** (read-only):
  - `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0;
  - `ls logs/imp3/` shows the four input JSONs and **no** `ours_imp3h_b4.json`;
  - `git diff 0f6bbbc7 -- tests/validation/test_imp3_tip_ratio_1g.py tests/validation/test_imp3_rod_birdcage.py` is empty;
  - `git diff fee90cf8 -- tests/mesh/test_birdcage_implant_rod_ball.py tests/mesh/test_birdcage_implant_rod_tip.py tests/mesh/test_birdcage_implant_rod.py` is empty;
  - `git diff 6c65dee5 -- tests/validation/test_imp3_tip_ladder_1g.py` is empty before the edit.
  - On any mismatch, stop and report.
  - **No cache sweep is owed:** no window compiles a form (W0 imports, W1 / W2 never call `_arm_d`), as item 464 ruled.
- **W0** (import, `-n 1`, knob unset; (N2) runs, then the test skips):
  `scripts/testing/run_and_log.sh IMP-3-I30a-W0-import "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short'"`
- **W1** (`b4build`, `-n 2`, the window of record; predicted ≈ 90 s, range [70, 200]):
  `scripts/testing/run_and_log.sh IMP-3-I30a-W1-b4build "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP7=b4build timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short'"`
- **Before W2:** `pgrep -c python3` reads 0.
- **W2** (`r8build` re-gate, `-n 2`; item 327 W1b's command with only the label changed; must reproduce (Kc8) 342 072, (G1x) 0.855, (N1) 1.069, (G3) 8 790, (Rb8) "3.556262" and the GEO-36 mesh identity, as at `…item327-W1b-r8build.log:1952–1960`; predicted ≈ 90 s, range [70, 200]):
  `scripts/testing/run_and_log.sh IMP-3-I30a-W2-r8build-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_IMP3_STEP7=r8build timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_imp3_tip_ladder_1g.py -v -s --tb=short'"`
- **Not re-gated, and why:**
  - `r4s2` and `r8` solve and write mode-`"x"` JSONs that already exist, so both would fail on their own files by design.
  - Their bodies are diff-unchanged, and the dispatch lines they reach (`:460–464`) are byte-unchanged. W0 proves the module imports and the unset path skips.
- **Order:** pre-check → code → W0 → W1 → `pgrep` → W2 → commit. **No module edit after W0 starts** (rule (i)); any edit means W0–W2 again.
- **Never run `b4` in-slot.** Item 464's scaled `-n 8` window is 552.6–564.7 s (`…I17-W2-n2.log:5671–5672`), over the 545 s cap, and the window is entry 43's.

**Traps already paid for.**
- **Item 305's list carries in full** (`IMP-3.md:1061–1082`): degree defaults to 1 in two places; `leg_sheet_strip_fraction` always; quadrature 4 at degree 2; no complex ordering or `sqrt` comparison; `owner="lowest_original_cell"`; rank-local values reduced and **no collective inside a `rank == 0` branch**; tags in `built[1]` / `built[2]`; evaluate before freeing; Status 124 means `pgrep` first; mode `"x"` means no second window; a denied `$rc` capture means stop and report.
- **Item 327's list** (`IMP-3.md:1162–1168`): `monkeypatch` must target the module object `_arm_d` lives in (`_ratio`). Patching the new module's imported name, or `test_imp3_rod_birdcage`'s, has no effect on `:165`. The existing `:463` line already does it right, so do not add a second `setattr`.
- **The patch covers only the `b4` / `b4build` / `r8` / `r8build` paths.** Never move the dispatch above `:463`, never patch on `r4s2` or on the unset path, and keep (N2) first. **Item 465's "unpatched, do not monkeypatch" instruction is superseded on this path by entry 43's licence. It is not a trap to obey.** Copying it would make entry 43 red at (G1) before the solve.
- **`_t1_geometry`, `_ball_and_t1`, `point_owner_cells` and `_outer_box_tags` are collective.** Call them on every rank.
- **`b4build` never calls `_arm_d`.** A build-only path that solves is an unpriced heavy window.
- **The 4 mm build has not been built since `src/fem_em_solver/io/mesh.py` moved.** `git diff fee90cf8 -- src/fem_em_solver/io/mesh.py` is +183/−10 from `GEO-38` items 478 / 470 (`090178c2`, `c1b8debc`, 2026-10-15; `phantom_box`-gated by design). A (Kb4) cell mismatch is therefore a finding about the shared builder, not a constant to re-record in-slot (see Negative result).
- **`tip_size` is never set.** The r/4 grain is the default, and (Kb4) asserts it.
- **`b4build` never appears in entry 43's command.** The implementer does not edit that command.
- Never `-k`; `-s` on every window; never lengthen a timeout or raise `-n`; never delete a log; never write `ours_imp3g_*`, `ours_imp3h_r4s2.json` or `ours_imp3h_r8.json`; no slot opens `docs/private/` or `aed_results/`.

**Scope.** Covered: the selectors, the disclosed patch path, the 4 mm / r/4-grain build identity in the complex build at `-n 2`, the near-T1-inside-the-ball predicate on that mesh, the three stored-reading paths (Rb / Rb8 / Rb-r8), a dry run of every readout line, and the `r8build` path re-gated. **Not covered:**
- any solve, η / Q reading or identity at b4 (all of that is entry 43's);
- convergence or an observed order; accuracy (`ANS-10`);
- entry 43's command, tier, date or price;
- other radii; F-human;
- `src/`; any band; the two guard values.

**Negative result.** Report and stop. Park on `attempt/IMP-3-I30a-<ts>`, mark item 532 BLOCKED in the same commit (rule (d)), and leave entry 43 PENDING PREREQUISITE.
- **W0 red:** a module defect. Fix it before W1 and re-run W0.
- **(Kb4) cells ≠ 351 478**, with the radius and grain bitwise: the 4 mm build moved since `fee90cf8`. Report both counts and the census.
  - Open a `GEO-36` / `GEO-38` known-issues entry naming `090178c2` / `c1b8debc`, and stop.
  - Re-recording is the review's (via `record-reconciler`), never in-slot. Entry 43 cannot run against an unverified mesh.
  - W2 still runs: it tells whether r/8 moved too.
- **(Kb4) radius or grain red:** the key or the default did not reach gmsh. Fix in-slot only if it is the new module's call; otherwise park.
- **(T1in) red** with (Kb4) green: a finding (the 4 mm ball no longer contains the near-T1 cells). Report the vertex distances; the weekly re-poses I30.
- **(G1x) / (N1) / (G3) red:** report every number and park.
- **(Rb*) red:** the reading path or a stored JSON is wrong. Print no dry-run reading and park.
- **W2 red:** the edit is not additive on the `r8build` path. Revert nothing on `main`; the review rules.
- **Status 124:** not a result. Run `pgrep` first, then park with the measured time.
- There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.**
- **§7 `IMP-3` stays ✅**, with no record sentence (item 327's precedent, `IMP-3.md:1178`). Its Done-when is met and not amended.
- **`docs/testing/xl-pending.md` entry 43, in the same commit:** mark **(i′) MET 2026-10-18 (§9 item 532, `<commit>`; <W0 log>, <W1 log>, <W2 log>)** in its Status paragraph, and report the selector `FEM_EM_IMP3_STEP7=b4` as landed.
  - **Never edit the entry's command.**
  - **Never write `docs/testing/xl-queue.d/20261023-IMP-3-step10-I30.env`.** The next review confirms the knob and queues.
- The slot makes no §10 edit; the review annotates I30a.
- No known-issues entry retires.

**Executor spawn prompt (verbatim).** "Execute §9 item 532 — `IMP-3` §10 I30a (the `b4` / `b4build` selectors for `xl` entry 43) — per `docs/planning/chunks/IMP-3.md` §9 item 532. Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Edit only `tests/validation/test_imp3_tip_ladder_1g.py`, additively: `_r4s2()` and `_r8()` untouched, the `:463` monkeypatch line untouched, the two `b4` dispatch lines inserted after it. Never import `builds` or a `test_*` name from `test_birdcage_implant_rod_ball.py`. Every collective is called on every rank, outside any rank-0 branch. Never run `FEM_EM_IMP3_STEP7=b4`. Run the pre-check, then W0, W1, `pgrep -c python3`, W2 in order, foreground; no edit after W0 starts. In `docs/testing/xl-pending.md` make only the entry-43 (i′) MET edit; never edit its command and never write a queue file. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/IMP-3-I30a-*`, mark item 532 BLOCKED in the same commit, and report."

**Problems for the review.**
- (1) **Patch on `b4build` in-slot.** Entry 43 licenses the patch for its own window. The build-only path inherits it through the existing `:463` line, as `r8build` did under item 327's ruling (1). The review should say so explicitly.
- (2) **`mesh.py` moved after item 464** (`GEO-38`, `090178c2` / `c1b8debc`). W1 is the first post-move 4 mm build, and W2 the first post-move r/8 build. A (Kb4) cell red blocks entry 43 for 10-23.
- (3) **Entry 43's Status asks the clerk to confirm the knob "against item 464's landing commit".** It should read item 532's landing commit; that is a text edit for the review. The command itself needs no change (see the report).
- (4) The dry-run readout, and the write-before-readout order on `b4`, depart from `_r8`'s shape. Both are proposed to protect the `xl` window; the review may strike either.

### Review ruling 2026-10-18 03:00 (item 532, I30a)

Accepted. Spot-check: `20261015T123700Z_IMP-3-I17-W2-n2.log:5671–5672`, the 4 mm build, 351 478 cells, 77.98 s, re-read.
- **The guard patch on `b4build`** is licensed on the build-only path, as item 327's ruling (1) licensed it for `r8build`. That path dispatches inside the existing `:463` patch, and the (G1x) / (N1) disclosure lines print. Item 465's "do not monkeypatch" is superseded on these two paths only, because entry 43 licenses the patch for its window.
- **Both departures from `_r8`'s shape are accepted.** `b4build` dry-runs the readout lines with the r/8 JSON as a stand-in. `b4` writes its JSON before printing.
- **A cell mismatch at W1 is a finding, not a red to work around.** `io/mesh.py` has changed since `fee90cf8`. Report the number, mark BLOCKED, and entry 43 waits.
- **Entry 43's command is unchanged.** Its Status (edited by this review) now names this item as prerequisite (i′).

**Outcome (2026-10-18 11:06Z, implementer slot 11:00 UTC) — DONE, green, attempt 1 of 4.** `b4` / `b4build` added to `tests/validation/test_imp3_tip_ladder_1g.py` additively (2 lines deleted: the docstring knob sentence and `KNOBS`; `_r4s2` / `_r8` and the `:463` patch line untouched; dispatch inside the patch). Pre-check clean (pgrep 0, four input JSONs, no `ours_imp3h_b4.json`, the three `git diff` checks empty). **W0** (`20261018T110240Z_IMP-3-I30a-W0-import.log`): 4 passed, 1 skipped (N3), 4 s. **W1** (`20261018T110251Z_IMP-3-I30a-W1-b4build.log`, 5 passed, 82 s, inside [70, 200]): (Rb) 3.482328, (Rb8) 3.556262, (Rb-r8) 3.530461, sizes integer-equal (`:63–65`); (Kb4) tip_ball_radius_m 0.004 and tip_size_m 0.00025 bitwise, **351 478 cells exact** (Δ 0, `:1959`) — so the 4 mm build did not move under `GEO-38`'s `090178c2` / `c1b8debc`; (G1x) 0.879 of 400 000, (N1) 1.098 over 320 000 with the DISCLOSURE line (`:1960–1962`); (T1in) 130 near-T1 cells, vertex max 3.151970e-03 m ≤ 4 mm, T1 owner 227 010 with vertex max 2.768817e-03 m (`:1963`); (G3) 5 366 tag-402 facets, all exterior; mesh identity vs item 464 W2 True on every line (d2 dofs 2 251 660 = 2E + 2F with E 418 532 / F 707 298, N_ball 111 415 / 3.489008e-04 m, near-T1 mean edge 3.383011e-04 m, T1 cell 227 010 / 3.282035e-04 / 2.768817e-03, `:1965–1967`); dry-run readout executed on every line with the r/8 JSON stand-in, η(T1) 3.530461 as predicted (`:1968–1988`); build 78.08 s (item 464 77.98 s), body 78.9 s, 0.90 GiB summed at `-n 2` (predicted ≤ 1.5). One cosmetic note for entry 43's reader: the dry-run's near-T1 edge-ratio line divides by the **b4 mesh's** measured near-T1 edge (2.0642), not the stand-in's; its label says "stand-in mesh". On `b4` the label is "b4" and the line is right. **pgrep** 0. **W2** (`20261018T110425Z_IMP-3-I30a-W2-r8build-regate.log`, 5 passed, 81 s): (Rb8) 3.556262, (Kc8) 342 072, (G1x) 0.855, (N1) 1.069, (G3) 8 790, GEO-36 mesh identity True — item 327 W1b reproduced, so r/8 did not move either. No `b4` run. `xl-pending.md` entry 43: (i′) marked MET; command untouched; no queue file written. §7 `IMP-3` stays ✅ with no record sentence.

### Review note 2026-10-25 03:00 — row amended from the 10-24 weekly's reading of `xl` entry 43 (I30)

The §7 row's state line now carries the four-mesh record. Source: the 10-24 weekly's §10 ruling on entry 43 (`PROJECT_PLAN.md` § 10, "Rulings, 2026-10-24"). η_rod(T1), P1, degree 2: 3.482328 / 3.556262 / 3.530461 (the three earlier rungs) and **3.539570** at the 4 mm cap ball, r/4 grain (`20261023T070008Z_IMP-3-step10-I30.log:2654`), inside the pre-registered [3.48, 3.56] ⇒ the decision rule's first branch; every I14 identity green (`:2630–2672`). Known-issues 2026-10-14 retires by its own clause (the 10-24 weekly). Two caveats travel with the record, both from the `log-pathologist`:
- the interval was derived from the three rungs it brackets, so "inside" separates little (10-24 epitaph);
- at T1 the b4 mesh is finer than r/8's (edge 3.383e-4 against 5.082e-4 m), so "r/4 grain" names the requested tip size, not the local cell.

The 1 g ball mass and the degree-2 dof match were printed, not asserted, on this path. **Reading of record:** the pointwise tip ratio at degree 2 is 3.48–3.56 across four meshes (±1.1 %), a record and not a converged value. No finer rung is bought. The 1 g ball ratio Q (4.4617 / 4.4610 top / bottom) is what (7c) compares when the F2182 reference arrives (W3). Glyph unchanged (✅ as a fixture).

### §9 item 626 — 2026-10-25 (§10 I39, `IMP-3` step 11; lane 3)

*Written by the 2026-10-25 03:00 daily review from the 10-24 weekly's enumerated step I39 (§9 rules (4) / (5)). Two execution readings below go to the 10-26 weekly; neither changes the target.*

**Question.** What does the rod's own SAR do to the F-small gel phantom's temperature, compared with the same coil drive without the rod? This is the first implant-driven heating on a coil-driven phantom.

**Execution readings (this review's):**
1. **"The two-port drive" on this mesh.** `IMP-3`'s F-small module drives four single ports P1–P4 and their **ccw / cw superpositions** (`tests/validation/test_imp3_rod_birdcage.py:145–146`). No two-port drive exists on F-small; that drive belongs to the F-human tuned rod (step 8). The source here is the **ccw superposition**, the coil's circular drive and the one whose 1 g tip ratio `IMP-3` records.
2. **Ranks.** The weekly wrote `-n 8`. The measured with-rod case is 171 s at `-n 2`: build 47.4 s, arm 53.2 s, 2.61 GiB summed peak RSS (`20261002T170428Z_IMP-3-step1-W1.log.gz:3646`, `:3719`, `:3764`). The thermal CN on a ≈ 60 k-cell submesh is 2.4–2.8 s per solve at `-n 2` (`20261001T123309Z_MAT-15-step1-W1.log.gz:1929–1944`). The smallest width that fits is `-n 2`, which is also the width that exposes a missing reduction.

**Route.** One new module, `tests/validation/test_imp3_rod_heating_phantom.py`, gated by `FEM_EM_IMP3_I39=1`. It runs `IMP-3`'s with-rod and without-rod arms on the 167 371-cell record mesh (`ROD_CELL_RECORD_W1`, `test_imp3_rod_birdcage.py:113`) at 64 MHz through `_case` / `_arm`, with `run_n_port_sparameter_sweep(…, keep_fields=True)`. Fields are not saved to disk, so the SAR is taken through the module's `field_hook(label, rod, msh, cell_tags, e_of, sigma_field, rho_field)` (`:405–406`) before the fields are freed.
- The ccw source is q = SAR · ρ on the phantom (tag 3) as a DG0 / DG1 Function interpolated onto the P1 thermal submesh.
- Then `MAT-15` step 1's scheme (`tests/validation/test_mat15_coil_driven_heating_phantom.py`): `_cn` θ = ½, N = 40, two BE half-steps (`:154–157`); `_be_one` at `DT_1` (`:89`, `:160`); t_R, t\* as `:86–87`.
- **Import these by name, never copy them:** `_cn`, `_be_one`, the energy helper, `E_BAND`, `S2_BAND`, `C4_COVARIANCE_BAND`, `B_BAND`. If a helper is bound to `MAT-15`'s own 10 MHz sweep (`_run_heated` is), call the pieces below it and say which.

**Anchors (imported bands, both arms).**
- **(E′)** Stored + surface = G · P_k at `E_BAND` (1e-9; `MAT-15` module `:164–171`, `:378–379`).
- **(S2)** The calorimetric first-step slope against q at `S2_BAND` (1e-4; `:384–385`).
- Both on the with-rod and on the without-rod source.

**Control (printed, predicted to miss).** **(C4)** The C4 covariance at `C4_COVARIANCE_BAND` (5e-2; `test_birdcage_b1_plus_map.py:120`). The rod breaks the coil's C4 symmetry, so the with-rod θ field should miss. It is printed with its miss, never asserted. The **without-rod arm's (C4) is asserted** (`MAT-15` asserted it on the uniform phantom, `:390–392`), so the comparison is backed.

**Printed (the record, no band):** the tip-region 1 g ball-mean θ(t\*) over the without-rod value at the same ball, beside `IMP-3`'s 1 g SAR ratio Q on the same mesh; the whole-phantom ΔT_max ratio; P_phantom per arm.

**Windows.**
1. W0: knob unset, coding tests and the knob refusal.
2. W1 `IMP-3-I39-W1`: both arms, `-n 2`, complex build, heavy, `timeout -k 30 1200`.

**Predicted** ≈ 2 × 171 s + ≈ 2 × 120 s of thermal (40 CN + the BE step + the energy forms) ≈ 550–700 s, ≤ 6 GiB summed. Run the `pgrep -c python3` pre-check (must read 0) first.

**Tier:** heavy, `-n 2`. **Class:** lane-3 closed form (the class median applies). **Cost:** ≈ 7.5 measured / 50 predicted slot-minutes. **Executor:** implementer. **Depends on:** nothing. **Family:** `IMP-3` step 11, attempt 1 of 4. **Daytime:** no.

**Traps already paid for.**
- **The rod wall must be adiabatic.** `_robin_tags` puts Robin on *every* exterior facet of the submesh (`test_mat10_…:92–98`), so the rod's hole wall (facet tag 402) would become a Robin sink. Apply Robin only on the phantom's outer surface (`IMP-7`'s trap 2, `IMP-7.md:382–395`), and say which tags.
- `MAT-9`'s `float(np.sqrt…)` is not complex-safe (`MAT-15.md:19`). Take `np.real` before any `float`.
- `ufl.max_value` and complex comparisons do not compile.
- `cell_tags.values` and `assemble_scalar` are rank-local.
- A killed compile leaves an orphan `.c`. Sweep only with `pgrep -c python3` reading 0 (`MAT-15.md:44`, `:73`).
- Never background the window.

**Scope boundary.** Not closed here: any published comparand (none exists until (7c), W3), compliance, the F-human or F-clinical rod, perfusion. `IMP-3` stays ✅ as a fixture. An `imp:` example follows only if this lands (the weekly's examples list).

**A negative result means:**
- (E′) or (S2) red on the without-rod arm: the port from `MAT-15`'s 10 MHz route is wrong. Report it and park; nothing is read.
- Red on the with-rod arm only: the rod wall's boundary is the suspect. Record it in `IMP-3.md` and on an `attempt/*` branch.
- A without-rod (C4) red: the source transfer is wrong.

**Status it can move.** §10 I39 → DONE. `IMP-3`'s row gains the first coil-driven implant-heating record (step 11). Phase 7's (7b) subgoal gains its heating half.

**Outcome (2026-10-25 11:06Z, implementer slot 11:00 UTC) — DONE, green, attempt 1 of 4.** New module `tests/validation/test_imp3_rod_heating_phantom.py` (knob `FEM_EM_IMP3_I39=1`), `MAT-15`'s pieces imported by name (`_tag_submesh`, `_cell_average_dg0_tag`, `_cn`, `_be_one`, `_e_prime`, `_max_imag`, the gel constants, `E_BAND`, `S2_BAND`; `B_BAND`, `ROBIN_TAG`, `_allsum`, `_to_sub` from `MAT-10`; `_run_heated` not called — it is bound to the 10 MHz sweep). Source: ccw superposition, one G = 1 W / P_phantom(without rod) for both arms. **Robin tags:** exterior facets of the tag-3 submesh farther than 5 mm from the rod's capsule axis segment; the rod wall is adiabatic, and the excluded count is asserted equal to the parent's owned tag-402 count. **W0** (`20261025T110246Z_IMP-3-I39-W0.log`): 4 passed, 14 skipped (knob unset), 2.3 s. **pgrep** 0. **W1** (`20261025T110300Z_IMP-3-I39-W1.log`, 18 passed, 177 s elapsed, heavy tier, `-n 2`, complex, `timeout -k 30 590` — one window, since the measured arms (`20261002T170428Z_IMP-3-step1-W1.log.gz`) priced it at ≈ 350 s, inside the 660 s foreground ceiling; arms 43.5 / 60.8 s, thermal 3.9 / 7.0 s): census 116 085 / 167 371 exact (`:3690`); rod wall 4 958 excluded = 4 958 tag-402, without-rod 0 (`:1841`, `:3662`); **(B)** 8.9e-16 / 1.7e-15 (`:1842`, `:3663`); **(E)** CN 5.3e-14 / 3.8e-14, BE 1.2e-11 / 3.6e-11 (`:1843`, `:3664`); **(E′)** 3.1e-15 / 5.0e-15 vs 1e-9 (`:1844`, `:3665`); **(S2)** 5.21e-6 / 5.00e-6 vs 1e-4 (`:1845`, `:3666`); **(C4) without rod** pairs 2.07 / 0.91 / 0.35 / 1.55 % vs 5 % (`:1848–1851`). **Record (printed):** P_phantom 1.034026e-5 W without, 1.054706e-5 W with (+2.00 %); 1 g ball θ(t\*) rod / without **1.1109 top / 1.0884 bottom** beside `IMP-3`'s 1 g SAR ratio Q 4.0263 / 3.6434 on the same meshes (`:3687–3688`); whole-phantom max θ ratio 1.0412 (1.6203 / 1.5562 K), mean θ ratio 1.0216 (`:3689`). Two findings for the review:
- **The with-rod (C4) control did not miss.** Its pairs read 2.96 / 3.66 / 2.61 / 3.79 %, all inside the 5 % band (`:3669–3672`), against a prediction of a miss. It was printed, never asserted, and it carried a "predicted" label. It is reported here as a wrong prediction, not a red. The consequence: at this band and at 10 g, (C4) does not separate the rod from mesh asymmetry (the worst without-rod pair is 2.07 %). The without-rod (C4) assertion therefore checks the source transfer, not rod sensitivity.
- **The without-rod thermal submesh has only 537 cells** (`:1840`), against 51 693 with the rod (`:3661`). The W1 solid mesh leaves the phantom coarse, while the rod mesh refines it around the rod. The θ ratios therefore compare a coarse P1 solve with a fine one; the 1 g ball (a = 6.2 mm) is not resolved on the 537-cell side. Read 1.11 / 1.09 as a record on this mesh pair, not a mesh-independent heating ratio. (The SAR ratio Q has always carried the same pairing.)

### Review ruling 2026-10-25 10:30 (item 626, I39)

Audit PASS (`20261025T110300Z_IMP-3-I39-W1.log:1840–1852`, `:3661–3673`, `:3687–3690`, `:3706`, `:3711–3713`). Two text corrections, no status change: the row's Done-when tail read "No heating or compliance claim", which step 11 contradicts — it now reads "No compliance claim; step 11's heating is a record"; and the item-626 line's "≤" figures are the maxima rounded to two digits (1.742e-15, 5.007e-15, 5.2118e-6, 2.0732 %), not bounds. The 537-cell finding stands as the executor wrote it. The like-for-like reading is posed below as I41 (rule (6), provisional).

### §9 item 639 — 2026-10-25 (§10 chain I, I41 — daily, 2026-10-25, provisional; chain step (solve))

*Written by the 2026-10-25 10:30 daily review under §9 rule (6), toward chain I's subgoal (7b) — the rod in the coil-driven phantom — and I39's validation target, unchanged: `MAT-15`'s identities at their imported bands on both arms, the tip 1 g θ ratio printed. The 10-26 weekly ratifies or reverts it.*

**Question.** I39's 1 g tip θ ratio (1.111 / 1.088) compares a 537-cell phantom (the solid record mesh, `phantom_resolution=None`) with a 51 693-cell one (the rod mesh, refined around the rod). What is the ratio when both arms share one mesh?

**The matched mesh exists.** The rod mesh of I39 is a **hole** build (no rod cells; tag 402 an exterior PEC facet group, `src/fem_em_solver/io/mesh.py:268–272`), so there is nothing to map to gel on it. Its kept-cells twin `_build_kept_rod()` (`tests/mesh/test_birdcage_implant_rod.py:154–168`, `KEPT_ROD_CELL_RECORD` 182 371 at `:102–108`) is "the hole build's 167 371 + 15 000 implant cells, the rest of the census identical cell for cell"; on it tag 402 is an interior interface (`mesh.py:274–280`) and the rod's cells carry tag 403 (`BIRDCAGE_IMPLANT_CELL_TAG`). With tag 403 mapped to the gel, the rod is gone and the mesh is the with-rod arm's, cell for cell outside the rod. `WF-9` step 1 already solved exactly this rod-as-gel problem on F-small (`tests/validation/test_wf9_lead_in_coil_phantom.py:151–172`: `p_without = replace(problem, material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: gel})`); item 381 did the same on F-human (`tests/validation/test_imp3_f_human_tuned_rod.py:110–146`).

**Route.** New module `tests/validation/test_imp3_rod_heating_matched.py`, knob `FEM_EM_IMP3_I41=1` (unset ⇒ skip with a reason; any other value fails a knob test, as I39's `:290–298`). Import I39's pieces by name from `tests/validation/test_imp3_rod_heating_phantom.py` (`SOURCE`, `WALL_SELECT_M`, `_heat_arm`, `_outer_robin_tags`, `_parent_wall_count`) and `MAT-15`'s thermal pieces as I39 does. No landed module edited, unless composing the rod-as-gel EM arm needs an override on `_arm` (`tests/validation/test_imp3_rod_birdcage.py:295`): then add one keyword, default `None`, under rule (c), and re-run that module's own gate window green in the same slot (`IMP-3` step 1, ≈ 241 s at `-n 2`) — disclose it. Prefer the `WF-9` route (build, `replace`, solve through the same ccw superposition and `G` normalisation I39 uses).
- **Arm M (matched, rod as gel):** `_build_kept_rod()`, tag 403 → the phantom's gel material, EM solve under `SOURCE`; **G = 1 W / P_phantom(arm M)**, where the phantom is tags {3, 403} — this arm replaces I39's without-rod arm as the normaliser.
- **Arm R (with rod):** I39's with-rod arm unchanged (`_build_rod()`, the PEC hole), normalised by arm M's G.
- **Thermal:** arm R exactly as I39 (submesh tag 3, rod-wall facets excluded from Robin). Arm M on a submesh over tags **{3, 403}** (a local two-tag variant of `_tag_submesh`; the DG0 source average over both tags), Robin on its whole outer boundary — the rod wall is interior there.

**Anchors (asserted; every band imported as I39 imports it).**
- **(A0m) the matched census:** arm M's global cells = `KEPT_ROD_CELL_RECORD` exactly; tag-403 cells = 182 371 − 167 371 = 15 000 exactly; arm M's tag-3 cell count equals arm R's tag-3 count (51 693, I39 `:3661`) exactly — the "cell for cell" claim, now asserted. Reduced over ranks before asserting.
- **(B), (E), (E′), (S2)** on both arms at `B_BAND` / `E_BAND` / `E_BAND` / `S2_BAND` (I39 `:345–368`).
- **(C4) on arm M** at `C4_COVARIANCE_BAND` (uniform gel, mesh refined off-axis — the test that the refinement does not break the drive's fourfold covariance). Arm R's (C4) printed, as I39.
- **(W0) Robin census:** arm M excludes 0 facets (no exterior rod wall); arm R excludes 4 958 = its parent tag-402 count (I39 `:3662`).

**Control (rule (e)).** *Asserted, backed:* arm R's (B)–(S2) reproduce I39's prints to the printed digits' order (`:3663–3666`) — the same arm, same mesh. *Predicted, printed:* arm M with tag 403 left unmapped (vacuum, the default at `test_port_birdcage_four_port.py:254`) changes P_phantom — print the change; do not assert it.

**Printed (the record this step exists for).** The tip 1 g ball-mean θ(t\*) ratio R / M at both tips, beside I39's mesh-pair 1.1109 / 1.0884 and beside the SAR ratio Q on the matched pair; the max-θ ratio; P_phantom(R) / P_phantom(M) − 1 beside I39's +2.00 %; both arms' phantom cell counts (51 693 and 66 693 predicted).

**Window.** W0: `--collect-only` plus the knob test, `-n 2`, smoke. W1 `IMP-3-I41-W1`: `-n 2`, complex, heavy, container `timeout -k 30 590` (the 660 s foreground cap), bare harness call, foreground, Bash `timeout: 660000`, after `pgrep -c python3` reads 0. **Predicted** ≈ 45 s kept-rod build + 60 s arm-M EM + 44 s rod build + 61 s arm-R EM + 2 × ≈ 8 s thermal ≈ 230–300 s.

**Tier:** heavy, `-n 2`, complex. **Class:** chain step (solve). **Cost:** ≈ 7.0 measured / 45 predicted slot-minutes. **Executor:** implementer. **Depends on:** nothing. **Family:** I41, attempt 1 of 4. **Daytime:** no. Never beside another heavy window.

**Traps already paid for.** Unmapped tag-403 cells silently become vacuum (`WF-9.md:262`, `:367`) — map them explicitly and assert (A0m) before any solve. `cell_tags.values` and `assemble_scalar` are rank-local — reduce. Point samples through `evaluate_vector_field_parallel`; call by call on one batch (item 537). `-k 30` on every `timeout`; never background a window; `pgrep -c python3` after any kill. A `SpatialCoordinate` facet integral needs its quadrature degree pinned. Never pipe pytest.

**Scope boundary.** A like-for-like record of I39's ratio. Not claimed: a converged θ ratio (one mesh), the rod's own conduction (the wall stays adiabatic on arm R), a heating band, `ANS-10`, compliance.

**A negative result means:** (A0m) red ⇒ the kept build is not the hole build plus the rod — stop and report (the comparison is void, a `GEO-36` finding). (C4) red on arm M ⇒ the off-axis refinement breaks the drive's covariance at this band — a finding, recorded with its pairs; stop. Any identity red ⇒ stop and report. Park on `attempt/*`, BLOCKED in the same commit (rule (d)). No `PREDICTED RED:` line.

**Status it can move.** §10 I41 → DONE (provisional until the weekly ratifies); the `IMP-3` row's step-11 sentence gains the matched ratio ("on a matched mesh pair, I41"), replacing the 537-cell caveat. `IMP-3` stays ✅. An `imp:` example of the rod heating its phantom is opened after this lands (the weekly's "an `imp:` example after I39").

**Outcome (2026-10-25 20:21Z, implementer slot 15:00 CDT) — DONE, green, attempt 1 of 4.** New module `tests/validation/test_imp3_rod_heating_matched.py` (knob `FEM_EM_IMP3_I41=1`); no landed module edited. **Route:** `WF-9`'s `replace(problem, material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: gel})`, applied by wrapping `build_four_port_sweep` as `test_imp3_rod_birdcage._arm` sees it, for arm M's call only (restored in `finally`); the problem's map had no 403 entry before the swap (printed `had False -> mapped True`). Arm M's heating calls I39's `_heat_arm` on a copy of the cell tags with 403 relabelled to 3, so P_phantom, the DG0 source and the thermal submesh span {3, 403}. Both arms' `mass_averaged_sar` balls carry `_arm`'s density map (ρ = 0 off tag 3), so the 1 g balls exclude the rod volume on both arms alike (masses 0.981846 / 0.982903 g on both). **W0** (`20261025T201033Z_IMP-3-I41-W0.log`): 3 passed, 9 skipped (knob unset), 1.4 s. **W1** (`20261025T201043Z_IMP-3-I41-W1.log`, 300 s): 15 passed, 1 failed — the module's own control `abs(S2 − I39) ≤ 5e-12` read 2.97e-11 (`:3795–3796`). That bound was set by this executor, not the brief, which asks for reproduction "to the printed digits' order". Arm R differs from I39's only in G (+0.26 %), and the BE step solves for T = 20 K + θ, so θ̄¹ carries round-off from cancellation. Re-set to 1e-10 (the last printed digit of I39's 4.9999e-06), with the measurement in a code comment. P_phantom(R) reproduced I39's 13 printed digits (1.054706186864e-05 W). **pgrep** 0 before each window. **W1b** (`20261025T201613Z_IMP-3-I41-W1b.log`, 16 passed, 290 s harness, heavy, `-n 2`, complex, `timeout -k 30 590`; arms M 71.0 / R 61.8 / control 61.4 s): **(A0m)** 182 371 cells = `KEPT_ROD_CELL_RECORD`, tag-403 15 000, tag-3 51 693 = arm R's 51 693, asserted before any solve (`:3645`); **(W0)** arm M rod-wall 0 (exterior 258 = outer 258), arm R 4 958 = parent tag-402 4 958 (`:3679`, `:3724`); **(B)** 5.9e-15 / 2.4e-15; **(E)** CN 3.2e-14 / 5.1e-14, BE 4.0e-11 / 5.4e-12; **(E′)** 2.4e-14 / 1.2e-14; **(S2)** 5.13e-6 / 5.00e-6 (M / R, `:3680–3683`, `:3725–3728`); **(C4) arm M** pairs 0.93 / 0.05 / 0.29 / 0.61 % vs 5 % (`:3686–3689`); arm R's (C4) printed 2.96 / 3.66 / 2.61 / 3.79 %, the same as I39's. **Control (asserted):** arm R's P_phantom to 1e-9 and (S2) to 1e-10 of I39's prints, PASS (`:3789`). **Control (predicted, printed):** with 403 left unmapped (vacuum), P_phantom moves −0.0263 % from arm M (`:3779`). **Record (printed, `:3776–3778`):** the 1 g tip θ(t\*) ratio R / M is **1.104515 top / 1.108684 bottom**, beside I39's mesh-pair 1.1109 / 1.0884. The 1 g SAR ratio Q on the matched pair is 4.090342 / 4.113695 (I39's pair: 4.0263 / 3.6434). The max-θ ratio is 1.013722 (I39: 1.0412), the mean-θ ratio 1.024921, and P(R)/P(M) − 1 = +2.2683 % (I39: +2.00 %). Phantom cells 66 693 (M) / 51 693 (R), as predicted. Reading: on the matched pair the top and bottom tips agree to 0.4 %, where I39's pair split them by 2.1 %. Arm M's (C4) is ≈ 2× tighter than I39's 537-cell without-rod arm (worst 0.93 % vs 2.07 %). One mesh, a record, not a converged ratio.

### Review ruling 2026-10-25 18:00 (item 639, I41) — audited PASS; the re-set reproduction bound accepted, with its reading recorded

The `auditor` returned **PASS**. Evidence re-cited:
- **(A0m)** exact before any solve: `test_imp3_rod_heating_matched.py:147–149`, log `20261025T201613Z_IMP-3-I41-W1b.log:3645`.
- **(B) / (E) / (E′) / (S2)** on both arms at the imported `B_BAND` (`test_mat10_sar_to_heat_sphere.py:61`) and `E_BAND` / `S2_BAND` (`test_mat12_transient_sar_heating_sphere.py:75–76`): source `:253–256`, log `:3680–3683` and `:3725–3728`.
- **(C4)** on arm M at `C4_COVARIANCE_BAND` (`test_birdcage_b1_plus_map.py:120`): source `:263`, log `:3686–3689`.
- 16 passed (`:3791`); Status 0, 290 s (`:3797–3798`).

**The W1 red and the re-set.** W1 failed on `|S2 − I39_ROD_S2| ≤ 5e-12`, reading 2.969e-11 (`20261025T201043Z_IMP-3-I41-W1.log:3794–3797`). That bound was the executor's own, in a new and unlanded module. The brief asked only for reproduction "to the printed digits' order". I39 prints `|ratio − 1|` as 4.9999e-06, five significant digits, so 1e-10 is that order. The relaxation is **not a loosening in the CLAUDE.md sense**: no imported band, no §4 anchor and no brief-named gate moved, and the measurement is recorded in the code (`:270–274`). Arm M's own (S2) ratio also differs in the 12th digit between W1 and W1b on identical inputs, so a floor near 1e-12 exists anyway. Two caveats are recorded, and neither changes the result:
- Arm R is not literally I39's arm. Its normalisation G is arm M's (+0.263 %, `W1b.log:3723` against `I39-W1.log:3661`), so only P and the mesh are I39's.
- The brief's "(B)–(S2) reproduce I39" control holds only for P_phantom and (S2). (B), (E) and (E′) are held to their bands, not to I39's values.

**Standing reading for future briefs:** a reproduction control names its numeric bound in the brief, or names the print whose last digit sets it. An executor that tightens beyond that and then relaxes after a red writes both values and the measurement, as here. I41 stays provisional until the 10-26 weekly ratifies chain I's step.
