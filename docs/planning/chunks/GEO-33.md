# GEO-33 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 149 — 2026-10-01 (chain H6a′ 0a; queued by the 10:30 review)

`GEO-33` step 0a — the TE₁₀₁ box with the y = b/2 facet plane duplicated, both lips tagged in DolfinX (implementer; smoke, `-n 1` / `-n 2`, real; ≈ 6 measured / ≈ 41 predicted min; depends on nothing).

- **Why.** A single-valued two-sided term on an interior facet is a resistive sheet; at the septum of `TH-14` step 1's box E_t = 0, so that sheet reads Q ratio 1.000 (2026-09-27 daily §6). Thick tape needs decoupled lips. Step 0a builds them; 0b solves on them with Q_sep/Q₀ = 0.771914 against Pozar.
- **Construction (numpy, no gmsh, no crack plugin).** The lower half [0, 1.0] × [0, 0.4] × [0, 0.6] from `dmesh.create_box(MPI.COMM_SELF, …, (9, 4, 6), tetrahedron)` (`EDGES` from `test_cavity_leontovich_q.py:42`). Mirror it across y = 0.4 on `TH-16`'s pattern (`test_th16_symmetry_plane_driven.py:116–154`, cells reordered `[1, 0, 2, 3]` for the reflection), **without merging the plane vertices** — skip the `np.unique` lattice step; the upper half's vertices are offset by n_v. Distribute with a local copy of that module's `_distribute` (do not import the complex-only module into this real-build test). Tag on the DolfinX side: boundary facets with \|y − 0.4\| < 1e-12, classified by the single adjacent cell's midpoint y (below `LIP_LOWER`, above `LIP_UPPER`); tag the six box faces too. `merge=True` builds the control (plane vertices merged, single-valued).
- **Code (tests only).** New `tests/mesh/test_geo33_cracked_cavity_lips.py` holding `_cracked_box(merge: bool)` (0b moves it into `src/`).
- **Anchors (asserted; exact construction identities).** (L1) `assemble_scalar(1·ds(LIP_LOWER))` and `…(LIP_UPPER)`, allreduced, each 0.6 within 1e-9 relative. (L2) 108 facets per lip (9·6·2); sorted facet midpoints gathered on rank 0 pair up between the lips within 1e-12 m (conforming twins). (L3) every lip facet has exactly one adjacent cell; no cell touches both lips. (L4) total boundary ∫ds = 2(ab + bd + ad) + 2ad = 4.96 within 1e-9; volume 0.48 within 1e-9. (L5) global counts 2 592 cells, 700 vertices. (L6) the same readings at `-n 1` and `-n 2`.
- **Negative control (predicted, binary; asserted as the binary).** `merge=True`: 0 lip facets, boundary 3.76 (24.2 % short of 4.96), 630 vertices — the single-valued sheet 0b uses as its 1.000 control.
- **Windows** (real build, ≈ 15 s each).
  - W1: `scripts/testing/run_and_log.sh GEO-33-step0a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/mesh/test_geo33_cracked_cavity_lips.py -s -v --tb=short'"`
  - W2: the same, labelled `GEO-33-step0a-W2`, at `-n 2`.
- **Traps.** `create_box` on `COMM_SELF`, then pass the arrays from rank 0 only (else duplicated). `ds` integrals and `exterior_facet_indices` are rank-local, and interprocess facets are not exterior: allreduce. The plane coordinate 0.4 must come out exact (4 × 0.1); if not, print the y-spread and stop — never widen 1e-12. `meshtags` indices sorted. y = 0.8 needs 8 divisions, not `TH-14`'s 7; 0b inherits (9, 8, 6). Harness bare, foreground, `timeout: 660000`.
- **Scope.** The mesh and its tags. Not: the Leontovich lip term, any eigen-solve, Pozar, a `src/` change (0b), the strip-sheet birdcage.
- **Status it can move.** The `GEO-33` row records step 0a (Done-when written by the 10:30 review); no glyph. H6a′ 0b becomes queueable on lips proven decoupled.
- **Negative result.** (L1)–(L3) red: DolfinX merges or misplaces the duplicated plane — known-issues entry, park, BLOCKED; the gmsh-crack route then becomes a `mesh-probe` item. (L6) width variance: a distribution finding, report.

**Outcome 2026-10-02 (19:30 CDT slot): DONE.** (L1)–(L6) green at `-n 1` and `-n 2`: each lip 0.6 to 1.7e-15, 108 / 108 twins at 0.0 m, 0 shared lip vertices, boundary 4.96 / volume 0.48 to ≤ 2.9e-14, 2 592 cells / 700 vertices; control `merge=True` 0 lip facets / 3.76 / 630 as predicted. `20261002T003221Z_GEO-33-step0a-W1.log`, `20261002T003236Z_GEO-33-step0a-W2.log`. 0b is queueable.

## §7 row moved verbatim 2026-10-01 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/GEO-33.md`)* **F-human-32: the 32-leg birdcage rung** — `birdcage_port_domain` at `n_legs = 32` (64 ring-gap ports on the high-pass topology), human-scale radius and length in the F-human class | ⬜ **Step 0a done 2026-10-02 (§9 item 149):** cracked box lips decoupled — each 0.6 m² to 1.7e-15, 108/108 twins at 0.0 m, no shared vertex, merged control 3.76 / 630 (`tests/mesh/test_geo33_cracked_cavity_lips.py`; `20261002T003221Z_GEO-33-step0a-W1.log`, `…003236Z_…-W2.log`, 8 s / 3 s); 0b next. Future work, operator directive 2026-09-15: the 3T research coil the project is ultimately compared against has 32 legs, so the production fixture eventually needs a 32-leg rung beside the 16-leg F-human. **Parameters are generic by rule:** leg count, and a coil length and ring radius in the ordinary published clinical-birdcage range, nothing else — no rung/strip widths, capacitor values, shield or end-ring geometry taken from any specific product (the operator holds NDA-covered sizing and it never enters this repo, tracked or not; treat it like the AED numbers). Sequence: after the F-human degree-order question (`WF-7` step 0b, 2026-09-19 XXL) is read — step 1 is a `GEO-25`-style cost probe (mesh only, CAD identities, the C32 terminal-ratio gate, cell count and price), then a single degree-1 drive priced for the XL tier. Not dated; the weekly places it in §10 once the F-human order price is known. Its downstream deliverable — the tuned 3 T body coil at 128 MHz — is `WF-11` (operator directive 2026-09-28). **Conductor model note (operator, 2026-09-20):** the research coil is built from copper *strips*, and the operator's own AED workflow models them as 2-D sheets with a Finite Conductivity boundary — so this rung should carry the coil as **zero-thickness sheets**, not round tubes: the same Leontovich `Z_s = (1 + j)/(σδ)` as `TH-14`, but applied on **both sides of an interior facet** (the tape is many skin depths thick at 64/128 MHz, so the two sides decouple and each sees plain `Z_s`; HFSS's thickness/layered option is not needed). That two-sided interior-facet form is a small extension of `TH-14`'s exterior-hole surface term and is this chunk's first step once scoped. `ANS-6` deliberately stays on round tubes to match the gated fixture; the strip geometry is a different fixture (edge current crowding changes `L` and the loss), not a different treatment. **Step 0 (H6a′, scoped 2026-10-01 10:30 review, rule (j)):** (0a, §9 item 149) the 1.0 × 0.8 × 0.6 box with the y = b/2 plane duplicated — both lips tagged in DolfinX, each lip a·d to 1e-9, twins to 1e-12 m, merged-plane control asserted; (0b) the two-sided Leontovich term on the lips, Q_sep/Q₀ = 0.771914 (Pozar, two half-cavities) at `TH-14`'s 5 % band, the merged single-valued sheet reading 1.000 as the control. **Done-when (§4):** 0a + 0b green ⇒ 🟡 (the lip term exists and is gated); ✅ additionally needs step 1 (the 32-leg cost probe) and the strip-sheet rung's single priced drive. 

### §9 item 170 — 2026-10-02 (chain H6a′ 0b; queued by the 03:00 review)

`GEO-33` step 0b: the two-sided Leontovich term on 0a's duplicated lips. The anchor is Q_sep/Q₀ = 0.771914 against Pozar. The control is the merged single-valued sheet, which should read 1.000. Executor: implementer. Tier: standard, `-n 2`, complex build. Cost: ≈ 16 measured / ≈ 45 predicted min. Depends on nothing; 0a landed `85739849`. Drafted by a read-only `Plan` agent and reviewed by the 03:00 review, which re-did the Pozar arithmetic and cut W2's ceiling to the standard tier.

- **Why.** Item 149 showed the lips are separate: 108/108 twin facets sit at 0.0 m, and no vertex is shared. The cut spans the whole y = 0.4 plane, so the cracked box is two separate 1.0 × 0.4 × 0.6 half-cavities. Each lip is an exterior facet with its own outward normal, so `TH-14`'s existing `ds` term (Jin §5.8.3; impedance condition Jin §1.5.3, eqs 1.54–1.57) already acts on both sides of the sheet. The new pieces are two:
  - an entry point that takes a supplied mesh;
  - an interior `dS` sheet, used only by the control.
- **Closed form.** `pozar_te10l_q_conductor` (`core/cavity.py:549`). For TE₁₀₁, k does not depend on b, so R_s is common to both boxes. Q_c ∝ b / [b(2a³ + 2d³) + ad(a² + d²)], with 2a³ + 2d³ = 2.432 and ad(a² + d²) = 0.816. That gives 0.4/1.7888 ÷ 0.8/2.7616 = **0.7719141** (checked by the review). At σ = 1e4, Q₀ = 801.6849, `TH-14`'s logged Q_c. The test computes the ratio from the function and asserts the literal against it to 1e-6.
- **Code (additive, opt-in; §9 rule (c)).**
  - Move `_cracked_box`, `_distribute` and `_tag` from `tests/mesh/test_geo33_cracked_cavity_lips.py` into a new `src/fem_em_solver/core/cracked_box.py` as `cracked_box_mesh(merge, comm, septum_tag=None)`.
    - The module must stay importable in the real build.
    - Keep `GhostMode.shared_facet`.
    - Fix the `geometry.dofmap` DeprecationWarning with the accessor the warning names.
    - `septum_tag` tags the interior y = 0.4 facets of the merged mesh (`locate_entities`, indices sorted and unique).
  - Add `solve_impedance_mesh_mode(msh, *, sigma_s_per_m, k0_sq, sheet=(tags, tag) | None, degree=2, nev=6)` to `cavity.py`.
    - It assembles curl-curl + γ⟨n×u, n×v⟩ ds over all exterior facets.
    - With `sheet` set it adds γ Σ± ⟨n±×u±, n±×v±⟩ dS(tag): a single-valued resistive sheet, R = Z_s/2 (Jin §1.5.4, eqs 1.60–1.63).
    - It returns the eigenvalues ordered by \|λ − k₀²\|.
    - Do not touch `_cavity_forms` or `solve_impedance_wall_cavity_mode`.
  - Point 0a's test at the new `src` module and delete its local copies.
  - New test: `tests/validation/test_geo33_two_sided_leontovich.py`. It skips on the real build, as `TH-14`'s does.
- **Anchors (asserted; `TH-14`'s 5 % band and `TH-9`'s 1 % band, both imported, never restated).**
  - (Q0) the merged mesh with walls only: Q₀ within 5 % of 801.6849. Predicted ≤ 0.5 %; the mesh is (9, 8, 6) mirrored, not `TH-14`'s (9, 7, 6).
  - (Qs) the cracked mesh: Q_sep within 5 % of `pozar_te10l_q_conductor((1.0, 0.4, 0.6), 1, 1e4)`.
  - (R) Q_sep/Q₀ within 5 % of 0.771914.
  - (F) Re f within 1 % of 291.3459 MHz.
  - (D) Im ω > 0 on all three solves.
- **Printed.**
  - Both members of the degenerate pair (one per half).
  - 1/Q_sep − 1/Q₀ (predicted 3.686e-4).
  - The ratio at σ = 1e6 (predicted 0.771914 again).
  - Lip areas and dof counts.
- **Negative control (rule (e)).** The merged mesh with the `septum_tag` sheet. TE₁₀₁ has E_t = 0 on that plane, so Q_sheet/Q₀ = 1.000.
  - Asserted: within 5 % of 1.000, and outside (R)'s band.
  - Predicted: \|Δ\| ≤ 1e-3, and a 29.5 % miss from 0.771914.
  - Ceiling: a sheet only adds loss, so the ratio is ≤ 1. The largest possible miss is therefore 1/0.771914 − 1 = 29.5 %, which is 5.9× the band. Assert no larger factor.
  - What it does not prove: it shows the sheet is single-valued, but a sheet that does nothing also reads 1.000. (Q0) against (Qs) is the half that bites.
- **Cost.** `TH-14` solves take 1.3–1.8 s each at `-n 2` (2 268 cells, 15 998 dofs), and its module windows took 58–77 s. 0a's windows took 3–8 s.
- **Windows** (Bash `timeout: 660000`; foreground; bare and whole).
  - W1, real build, 0a re-gate, must reproduce item 149's readings: `scripts/testing/run_and_log.sh GEO-33-step0b-W1-0a "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 2 python3 -m pytest tests/mesh/test_geo33_cracked_cavity_lips.py -s -v --tb=short'"`
  - W2, closing, includes rule (c)'s `TH-14` re-gate, which must reproduce Q = 8.017679e+02: `scripts/testing/run_and_log.sh GEO-33-step0b-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_geo33_two_sided_leontovich.py tests/validation/test_cavity_leontovich_q.py -s -v --tb=short'"`
- **Traps already paid for.**
  - No ordering comparison or `sqrt` on complex operands inside a form; tag classification stays in numpy.
  - `assemble_scalar`, `exterior_facet_indices` and `meshtags` are rank-local, so allreduce. At `-n 2` each rank may hold one half.
  - Run pytest with `-s`, and never `-k`.
  - Harness calls are bare, foreground, Bash `timeout: 660000`, never piped.
  - Sweep orphan `.c` files only after a killed compile, and only when `pgrep -c python3` reads 0.
  - Pin `quadrature_degree` if a `SpatialCoordinate` enters any facet integral.
  - Make no edit after the last window (rule (i)).
- **Scope.** The lip term and its gate on separate half-cavities only. Not covered:
  - a partial crack (tape edges inside the domain), which is the later strip step's;
  - a strip mesh, gmsh, birdcage sheets, copper gating, `GEO-33` step 1.
- **Status it can move.** §7 `GEO-33` ⬜ → 🟡, on the row's Done-when "0a + 0b green ⇒ 🟡", which this item matches.
- **Negative result.** Report and stop; no band moves.
  - (Q0) red: the mirrored mesh is at fault, not the lip term. Write a known-issues entry.
  - (R) red with (Q0) green: the lip path is at fault. Park on `attempt/GEO-33-step0b-…` and mark the item BLOCKED.
  - Control red: the sheet restriction or its normal is wrong. Write a known-issues entry.
  - One pair member unconverged: print it; that is not a failure.
  - No red is expected, so the item carries no `PREDICTED RED` line.

**Outcome 2026-10-02 (04:30 CDT slot): DONE.** Code: `src/fem_em_solver/core/cracked_box.py` (`cracked_box_mesh(merge, comm, septum_tag=None)`, moved from 0a's test, `dofmaps[0]` accessor, `shared_facet` kept, the septum tagged by `locate_entities`); `solve_impedance_mesh_mode` + `MeshImpedanceModes` appended to `core/cavity.py` (`_cavity_forms` and `solve_impedance_wall_cavity_mode` untouched); 0a's test repointed; new `tests/validation/test_geo33_two_sided_leontovich.py`. W1 (0a re-gate, real, `-n 2`, 1.5 s): every item-149 reading reproduced (lips 0.6 to 1.7e-15, 108/108 twins at 0.0 m, 0 shared vertices, 2 592 / 700, control 3.76 / 630). W2 (complex, `-n 2`, 18 passed in 67 s; the 0b module's solves 12.4 s): Pozar ratio 0.77191413 against the literal to < 1e-6; (Q0) merged walls Q₀ = 801.7516 vs 801.6849 (0.0083 %, 18 202 dofs); (Qs) cracked Q_sep = 619.2100 on both pair members vs 618.8319 (0.0611 %, 18 772 dofs, 10 converged); (R) 0.772322 vs 0.771914 (0.053 %); (F) Re f 291.1745 / 291.1208 MHz (≤ 0.077 %); (D) Im ω > 0 on all three; 1/Q_sep − 1/Q₀ = 3.677e-4 (predicted 3.686e-4); σ = 1e6 ratio 0.771954 (predicted 0.771914). Control (septum sheet, 108 owned facets): Q_sheet/Q₀ = 1.00000000 (\|Δ\| 6.7e-12, predicted ≤ 1e-3), miss from 0.771914 29.548 % = the ceiling, 5.91× the band. `TH-14` re-gate Q = 8.017679e+02 reproduced. Every prediction landed. `GEO-33` ⬜ → 🟡. `20261002T095630Z_GEO-33-step0b-W1-0a.log`, `20261002T095638Z_GEO-33-step0b-W2.log`.
