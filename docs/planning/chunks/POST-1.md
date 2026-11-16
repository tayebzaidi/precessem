# POST-1 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record. Steps 1–6 are archived in `docs/planning/plan-archive.md` § "§7 POST-1 adjudication + steps 1–6".

### §9 item 187 — 2026-10-02 (lane 1)

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 18:00 daily review.*

This is the sphere half (7a) of the 09-30 weekly's last-step target for `POST-1`. 7b — the same residuals on the loaded birdcage as a C4-gated record — reuses 187's extraction, is serial on it, and is not written.
- **Executor:** implementer. **Tier:** heavy (W1 carries a ≈ 195 s `POST-8` re-gate; the item's own module is standard). `-n 2` / `-n 4`, complex, plus one real-build window.
- **Cost:** ≈ 4 measured / ≈ 35 predicted min. **Depends on:** nothing.

**Code (rule (c), additive, disclosed).** In `src/fem_em_solver/post/evaluation.py`:
- Add a keyword-only `candidate_cells=None` (local cell indices) to `point_owner_cells` and to `evaluate_vector_field_parallel`.
- It is valid only with `owner="lowest_original_cell"`; anything else raises `ValueError` before any collective.
- Candidates are the owned cells ∩ the set; the tie-break is unchanged.
- `candidate_cells is None` joins the OPS-40/55 allgather key.
- With `None`, both functions are byte-identical.
- A point with no candidate on its side returns `valid=False`.

**Modules.**
- **(K) `tests/post/test_one_sided_owner_rule.py`** (both builds). On `create_unit_cube(4,4,4)` at the face point (0.5, 0.3, 0.6), restricted to cells with centroid x < 0.5 and then x > 0.5, each call returns the min of a DG0 `original_cell_index` field over the incident set on that side, with imaginary part 0.
- **`tests/post/test_post1_one_sided_interface_sphere.py`** (`complex_only`).
  - Imports `_series`, `_mesh_and_solve`, `SPHERE_TAG`, `RESOLUTIONS_64` and `RESOLUTIONS_128` from `tests/validation/test_lossy_sphere_fullwave.py` (unedited).
  - Rungs: `RESOLUTIONS_64` + `RESOLUTIONS_128[-1]` (5 866 / 17 667 / 55 241 cells), 64 MHz, degree 1, `fields.e_complex`. Air is tag 2 (`io/mesh.py:2251`).
  - **Interface points:** on each rank, every face of an owned tag-1 cell whose three vertices satisfy ||v| − a| ≤ 1e-9·a.
    - The vertices are sorted lexicographically; the centroid p is their mean.
    - n_f is from the cross product of the sorted vertices, oriented so n_f·p > 0.
    - The points are allgathered and sorted by (inside cell's `original_cell_index`, sorted coordinates), so the array is bitwise identical on all ranks.
  - **E⁻, E⁺:** the owner rule restricted to `cell_tags.find(1)` / `find(2)`.
  - **Series values:** S⁻ = `internal_field(p)`; S⁺ = `incident_field(p) + scattered_field(p)`.
  - **Residuals:**
    - R_n = ‖ε_c E⁻·n_f − E⁺·n_f‖₂ / ‖E⁺·n_f‖₂;
    - e∓ = ‖E∓ − S∓‖₂ / ‖S∓‖₂ over all three components.

**Asserted.**
- **(S0) Premise, closed form only.** At p̂ = a·p/|p|, both of these are ≤ 1e-9:
  - ‖(ε_c S⁻ − S⁺)·r̂‖/‖S⁺·r̂‖;
  - ‖r̂×(S⁻ − S⁺)‖/‖S⁺‖.

  Predicted ≤ 1e-12: each series term matches termwise.
- **(X0a)** F = 2V − 4 exactly. F is the number of interface faces and V the number of distinct surface vertices (a closed genus-0 surface).
- **(X0b)** Every point has exactly one owned containing cell of tag 1 and one of tag 2, counted globally with `Allreduce` SUM.
- **(X1)** max |n_f×(E⁻ − E⁺)| ≤ 1e-10·max|E⁺|. The N1curl tangential trace is continuous, so this proves the two calls read the two cells of the same face.
- **(C)** R_n, e⁻ and e⁺ each strictly decrease over the three rungs (`TH-10`'s "decreasing rung to rung" precedent).
  - The levels are a record, predicted in [0.01, 0.30] at the finest rung.
  - The pairwise rates are printed, predicted ≈ 1.
- **(B)** The legacy `owner=None` route for both sides gives R_n = |ε_c − 1| (160.155 at 64 MHz) to 1e-12 relative.

**Controls (rule (e)).**
- (B) is algebraic, and asserted above.
- The separation |ε_c − 1|/R_n is **predicted** ≥ 500× and printed, not asserted.
- Printed beside R_n:
  - the σ-blind residual (εᵣ in place of ε_c), **predicted** 0.874;
  - the conjugated convention (ε_c*), **predicted** 1.748 — `TH-10`'s 173.8 % signature;
  - the fraction of points where the legacy route lands on the tag-1 side, **predicted** ≈ 0.5.

**Windows.** Bare, foreground, Bash `timeout: 660000`.

Basis: the `TH-10` module ran five solves in 21 s (`20260822T123746Z_OPS-18-step3-th10-rerun.log`), and `POST-8` step 2's W1 took 195 s.
- **W1 (`-n 2`, ≈ 270 s):**

  `scripts/testing/run_and_log.sh POST-1-step7a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py tests/post/test_one_sided_owner_rule.py tests/validation/test_post8_owner_rule_ans2.py tests/post/test_post1_one_sided_interface_sphere.py -v -s --tb=short'"`
- **W2 (`-n 4`, ≈ 60 s):**

  `scripts/testing/run_and_log.sh POST-1-step7a-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 4 python3 -m pytest tests/environment tests/post/test_one_sided_owner_rule.py tests/post/test_post1_one_sided_interface_sphere.py -v -s --tb=short'"`
- **W3 (real, ≈ 50 s):**

  `scripts/testing/run_and_log.sh POST-1-step7a-W3-real "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/post/test_evaluation_collective_guard.py tests/post/test_point_owner_rule.py tests/post/test_one_sided_owner_rule.py -v -s --tb=short'"`

After W2, print the R_n difference between `-n 2` and `-n 4` as a record (predicted ≤ 1e-9).

**Traps:**
- **Never `total_field` for S±.** It switches on r < a, and each centroid sits inside by the chord sagitta.
- **Point lists must be bitwise equal on every rank,** or OPS-55 raises. Build from sorted vertices, and never round.
- **`cell_tags.values` is rank-local;** candidates are owned cells only.
- **A face point found on one side only** is a collision-tolerance finding. Print it, and never tune the 1e-9.
- **The `POST-8` re-gate loads `ans:2` without `main()`** (item 177's precedent). No slot opens `docs/private/` or `aed_results/`.
- **Rule (i):** no edit after W1 without re-running it.
- **Harness hygiene:** `-s`, no `-k`, no pipes. After a kill, run `pgrep -c python3` in the container first.
- **File overlap:** item 188 imports from `post/evaluation.py` but does not edit it. The two items are independent.

**Scope.** 7a only: degree 1, 64 MHz, pointwise interface residuals. No peak, SAR, birdcage or ✅ claim, and no band on the levels.

**Negative result.** Report and stop; never loosen a bound.
- **(S0) red:** a known-issues entry against `LossySphereSeries` (the `TH-10` records are implicated).
- **(X0a) / (X0b) red:** a mesh or extraction finding. Park and mark BLOCKED (rule (d)).
- **(X1), (K) or a re-gate red:** revert `src/` and park.
- **(C) red:** record the levels on the row, mark BLOCKED, and the review rules.
- There is no `PREDICTED RED` line.

**Status it can move.** The §7 `POST-1` state line gains "7a done" and stays 🟡 (or a known-issues entry if red). 7b becomes writable.

**Outcome — DONE 2026-10-03 (21:00 CDT slot, executor `implementer`).** Every asserted check green; nothing edited after W1 (rule (i)).
- **Windows:** W1 `20261003T020316Z_POST-1-step7a-W1.log` (`-n 2`, 215 s, 27 passed, 1 skipped — the pre-existing `test_prechange_control_was_silent`; `POST-8` `ans:2` re-gate green); W2 `20261003T020658Z_POST-1-step7a-W2.log` (`-n 4`, 67 s, 16 passed); W3 `20261003T020808Z_POST-1-step7a-W3-real.log` (real, `-n 2`, 53 s, 20 passed, 5 skipped). A smoke dev run precedes them: `20261003T020248Z_POST-1-step7a-dev.log`.
- **(K):** the face point resolves to cell 223 (x < 0.5) and 232 (x > 0.5), each the oracle minimum on its side; the unrestricted rule gives 223. Same at `-n 2`, `-n 4` and real.
- **(S0):** normal 5.5e-16 … 6.1e-16, tangential 2.8e-16 … 2.9e-16 (bound 1e-9).
- **(X0a):** F / V = 494 / 249, 1104 / 554, 2496 / 1250 — F = 2V − 4 exact. **(X0b):** one tag-1 and one tag-2 cell per point on every rung.
- **(X1):** 2.3e-15, 3.8e-15, 5.8e-15 (bound 1e-10).
- **(C), a record** (5 866 / 17 667 / 55 241 cells): R_n 0.5007 / 0.2821 / 0.1778; e⁻ 0.0718 / 0.0366 / 0.0202; e⁺ 0.2810 / 0.1900 / 0.1328. All strictly decreasing, finest-rung levels inside the predicted [0.01, 0.30]. Pairwise rates: R_n 1.414, 1.142; e⁻ 1.658, 1.473; e⁺ 0.963, 0.886.
- **(B):** legacy R_n = |ε_c − 1| = 160.155313 to ≤ 3.5e-16.
- **`-n 2` vs `-n 4` R_n** (record): 2.2e-13, 4.8e-13, 2.1e-12.
- **Controls (printed):** separation 319.9× / 567.8× / 900.8× (finest meets the predicted ≥ 500×); σ-blind 0.893 / 0.876 / 0.873 (predicted 0.874); legacy tag-1 fraction 0.478 / 0.488 / 0.506. **The ε_c\* control misses its prediction:** 2.153 / 1.974 / 1.889 against the predicted 1.748 (+8 % at the finest rung). It still separates from R_n by ≈ 10×, but `TH-10`'s 173.8 % signature value is not reproduced. Under rule (e) this is reported, not decided; the review rules.
- **Code (rule (c), disclosed):** `candidate_cells=None` keyword on `point_owner_cells` / `evaluate_vector_field_parallel` in `post/evaluation.py`, per the brief's letter. New: `tests/post/test_one_sided_owner_rule.py`, `tests/post/test_post1_one_sided_interface_sphere.py`.

### §9 item 203 — 2026-10-03 (lane 1)

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.*

This is 7b, the birdcage half of the 09-30 weekly's last-step target ("the same residuals on the loaded birdcage as a C4-gated record"). It reuses 7a's one-sided extraction (`candidate_cells`, item 187) unchanged.
- **Executor:** implementer. **Tier:** heavy (W1 is predicted at 170–260 s, above the standard ceiling). `-n 2` (W1) and `-n 4` (W2), complex only. No real-build window: no `src/` change.
- **Class:** lane-1 revalidation (measured median 8 min). **Cost:** ≈ 8 measured / ≈ 40 predicted min. **Depends on:** nothing (7a is done).

**Fixture.** The loaded F-small four-leg gapped birdcage on the 1 g rung at 64 MHz, degree 1:
- `build_four_port_sweep(frequency_hz=FREQUENCY_64_HZ, phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG, build_only=True)`, then four `_solve_driven(sweep, pid)` (`test_birdcage_b1_plus_map.py:142`). `build_only` + `_solve_driven` is `test_port_birdcage_four_port.py:321`'s idiom. If it raises on a missing key, use POST-8 step 2's full-sweep `case` construction (`test_post8_owner_rule_ans2.py:114–148`) with the frequency keyword added. Edit no helper.
- **Why this rung.** C4 gates already run on it: `MAT-15`'s 10 g pairs (`test_mat15_coil_driven_heating_phantom.py:50–51`, `:390`). `IMP-3` W2 solved it at 64 MHz (`20261002T170843Z_IMP-3-step1-W2.log:1918`, 199 920 cells; build 38.7 s, sweep 21.3 s at `-n 4`, `:3809`). Its 58 866 tag-3 cells (h ≈ 2.5 mm) give ≈ 7 700 interface faces. The default mesh has 537 tag-3 cells (`test_birdcage_phantom_resolution.py:64`, ≈ 1.5 cm), which is too few for a quadrant-level reading.
- **The band, imported and never restated:** `C4_COVARIANCE_BAND` = 5 % (`test_birdcage_b1_plus_map.py:120`). It is the band the 64 MHz field gate (`test_birdcage_b1_larmor.py:69–83`), the coil-driven SAR integral gate (`test_birdcage_sar_integral.py:117–121`) and `MAT-15` import.

**Module `tests/validation/test_post1_one_sided_interface_birdcage.py`** (`complex_only`). Nothing under `src/`.
- **Imported, never restated:**
  - `C4_COVARIANCE_BAND` and `_solve_driven` (`test_birdcage_b1_plus_map.py`);
  - `PARTITION_EPS_M` and `PARTITION_RTOL` (`test_birdcage_sar_integral.py:142`, `:149`);
  - `QUADRATURE_STEP_DEG` and `_port_index` (`test_birdcage_b1_quadrature.py`);
  - `ONE_GRAM_RUNG_CELL_RECORD`, `ONE_GRAM_RUNG_PHANTOM_CELL_RECORD` and `PHANTOM_RESOLUTION_1G_RUNG` (`test_birdcage_sar_1g_rung.py:114–126`);
  - `PHANTOM_CELL_TAG` (`test_port_birdcage_lumped_column.py:101`; air is tag 2 per its comment at `:96–99`, so define `AIR_CELL_TAG = 2` citing it);
  - `PHANTOM_RADIUS` and `PHANTOM_HEIGHT` (`test_birdcage_port_tags.py:28–29`; the phantom is centred on the origin with axis z);
  - `_phantom_sigma_range` and `PHANTOM_SIGMA_RTOL` (`test_birdcage_sar_map.py:239`, `:201`);
  - `_series`, `FREQUENCY_64_HZ` and `SALINE_SIGMA` (`test_lossy_sphere_fullwave.py`). Take ε_c = `_series(FREQUENCY_64_HZ).epsilon_c`, exactly 7a's convention.
  - `CAD_AREA_M2` (`test_mat15_coil_driven_heating_phantom.py:99`, printed only).
- **Drive order:** `sorted` sheet ids, with k from `_port_index` (P1 at 360° is k = 0). The rotation sense is imported from `test_birdcage_sar_integral.py:337–345`: drive k+1 at Rx equals drive k at x, so P_{j+1}^{(k+1)} = P_j^{(k)}.
- **Interface faces (partition-free, rank-safe):**
  - On each rank, take every face of an **owned** tag-3 cell whose three vertices are all on the lateral surface (‖(x, y)‖ − R ≤ 1e-9·R in absolute value, |z| ≤ (H/2)(1 + 1e-9)) **or** all on one cap (||z| − H/2| ≤ 1e-9·H/2, ‖(x, y)‖ ≤ R(1 + 1e-9)). Apply the test per face, never per vertex.
  - Store each face's lexicographically sorted vertices, and orient n_f outward: n_f·(p − c) > 0, where c is the tag-3 cell's centroid.
  - Allgather the faces and lexsort them on the nine sorted coordinates, as at `test_post1_one_sided_interface_sphere.py:107–111`. **Drop every face that appears twice:** those are interior chords between two tag-3 cells with all three vertices on the surface. An interface face has exactly one tag-3 cell, so it appears exactly once. Print the dropped count (predicted 0–30).
  - **Points:** three per face at barycentric (2/3, 1/6, 1/6) and its permutations, all interior to the face. This rule is exact for the quadratic |·|² of a degree-1 N1curl trace. The array is bitwise identical on every rank (OPS-55).
- **E⁻, E⁺:** `evaluate_vector_field_parallel(e, pts, comm, owner="lowest_original_cell", candidate_cells=cell_tags.find(3))` for E⁻, and the same with `find(2)` for E⁺. The **legacy route** is the default `owner=None`, read once.
- **Residual integrals** per drive k, with quadrant weights w_j as a numpy twin of `_quadrant_weight` (`test_birdcage_sar_integral.py:207–220`, same `PARTITION_EPS_M`):
  - r = ε_c E⁻·n_f − E⁺·n_f;
  - N_j^{(k)} = Σ_f (A_f/3) Σ_q w_j(x_q)|r(x_q)|²;
  - D_j^{(k)} = the same sum with |E⁺·n_f|²;
  - N^{(k)} = Σ_f (A_f/3) Σ_q |r|², and R_n^{(k)} = √(N^{(k)}/D^{(k)}).

**Asserted.**
- **(A0)** Global cells == `ONE_GRAM_RUNG_CELL_RECORD` (199 920) and tag-3 cells == `ONE_GRAM_RUNG_PHANTOM_CELL_RECORD` (58 866), exactly, both reduced. The four sheets sit on the 90° grid (`_port_index`), and f = `FREQUENCY_64_HZ`.
- **(A1)** The phantom σ is uniform at `SALINE_SIGMA` to `PHANTOM_SIGMA_RTOL` (`_phantom_sigma_range`, reduced). |ε_c − 1| is printed beside 7a's 160.155313.
- **(X0a)** F = 2V − 4 exactly over the kept faces: the phantom is a closed genus-0 cylinder. V is counted as distinct `input_global_indices`, as in 7a.
- **(X0b)** Every point has exactly one owned containing cell of tag 3, one of tag 2, and two in total, counted globally with `Allreduce` SUM (7a's `_side_counts`, plus the total). Both one-sided calls return `valid` everywhere.
- **(X0c)** Σ A_f equals `MAT-15`'s ∮ds on the same mesh, 2.072479761e-02 m² (`20261001T123309Z_MAT-15-step1-W1.log:1926`), to 1e-8 relative. It is a second route to the same facet set; store it as a module constant with that provenance. Print the deviation from `CAD_AREA_M2`.
- **(X1)** Per drive, max |n_f×(E⁻ − E⁺)| ≤ 1e-10·max|E⁺| (7a measured ≤ 5.8e-15).
- **(P)** Partition: |Σ_j N_j^{(k)} − N^{(k)}| ≤ `PARTITION_RTOL`·N^{(k)}, and likewise for D, for all four drives.
- **(B)** The legacy route gives R_n = |ε_c − 1| to 1e-12 relative for every drive. This is algebraic: both sides read the same cell. 7a measured ≤ 3.5e-16.
- **(G1) The C4 gate on whole-interface residuals:** the four cyclic pairs |N^{(k+1)} − N^{(k)}|/N^{(k)} are each ≤ `C4_COVARIANCE_BAND`. The whole surface is its own C4 image. **Predicted** ≤ 2 %.
- **(G2) The C4 gate on the residual map at quadrant resolution:** the sixteen cyclic pairs |N_{j+1}^{(k+1)} − N_j^{(k)}|/N_j^{(k)} (k, j = 0…3, indices mod 4) are each ≤ `C4_COVARIANCE_BAND`, one parametrised test per pair so that a miss names its (k, j). **Predicted** in [0.3 %, 8 %]. The upper end is above the band; see Negative result.

**Records (printed, no band):**
- R_n^{(k)} per drive, **predicted** in [0.05, 0.6] (7a's sphere: 0.50 → 0.18 over three rungs);
- the 4 × 4 tables N_j^{(k)}, D_j^{(k)} and quadrant R_n;
- F, V, the number of points;
- elapsed time for the build, each solve, the extraction and the window;
- the D-pair table: the same sixteen pairings on D, **predicted** ≤ 1 % (`test_birdcage_sar_integral.py:172–185`, which measured 0.2–1.5 % on the coarse mesh);
- the single-drive mirror |N_{k+1}^{(k)} − N_{k−1}^{(k)}|/N_{k−1}^{(k)}, **predicted** ≤ the band;
- P1–P4 phantom power ½∫σ|E|² dx(3) beside `IMP-3` W2's 1.530006187e-06 / 1.527879591e-06 / 1.527774714e-06 / 1.528062711e-06 W (`…IMP-3-step1-W2.log:3802–3805`; a different solve path, `built=` + sweep `keep_fields`), **predicted** ≤ 1e-6 relative.

**Why the gate is at quadrant resolution, not pointwise (a deviation from the 03:00 design letter, for the review).** A pointwise residual map read against its rotated image is ill-posed on this fixture. The mesh is not C4-congruent: `c4_congruent_sheets` binds only the port sheets (`io/mesh.py:3748`). So a rotated face point is neither a face point nor on the faceted surface, and one of the two one-sided reads would be invalid. The fixture's own history also prices the pointwise construction: the pointwise SAR C4 identities missed this same 5 % band by 25–41 % (`test_birdcage_sar_map.py:24–32`). The same solves met it as cell integrals at ≤ 1.52 % (`test_birdcage_sar_integral.py:172–185`). R_n is itself a cell-level discretisation error. The smooth partition makes w_{j+1}(x) = w_j(R⁻¹x) exact, so (G2) is the rotated map at the finest resolution where the identity is exact on a non-congruent mesh.
- **Printed only:** a nearest-face pointwise reading, for the record. For each face f of drive k, take the face g whose centroid is nearest R^{−k}p_f, by brute force on rank 0 in row chunks (no scipy dependency). Print the median and 90th percentile of |ρ_g^{(0)} − ρ_f^{(k)}|/ρ_f^{(k)} (ρ = the face's mean |r|²) and the median match distance divided by h. **Predicted** median ≥ 10 % (the pointwise floor).

**Negative controls (rule (e)).**
- **(B) is asserted** (algebraic, above). Its separation |ε_c − 1|/R_n^{(k)} is **predicted** ≥ 250× and printed. **Ceiling:** the separation *is* 160.155/R_n by construction, so it is whatever R_n allows: 267× at R_n = 0.6 and 900× at 7a's finest sphere level.
- **The mis-rotated pairing, predicted and printed:** |N_{j+2}^{(k+1)} − N_j^{(k)}|/N_j^{(k)}, i.e. the 180° image read where the 90° one belongs (`test_birdcage_sar_integral.py:345`'s control, on N). There is no prior measurement on N, so it is not asserted. Print per k its mean, its ratio to the (G2) mean, and the **ceiling**: under C4 the control reads drive k's own adjacent-quadrant contrast |N_{j+1}^{(k)} − N_j^{(k)}|/N_j^{(k)}, so it cannot exceed max_j N_j^{(k)}/min_j N_j^{(k)} − 1. Print that bound per drive. **Predicted:** control in [10 %, 150 %] and ratio ≥ 3×, never above the printed bound.
  - **Trap:** never use the wrong-*sense* pairing (j − 1) as a control. Under C4 it maps to N_{j−2}^{(k)} vs N_j^{(k)}, and for j = k+1 that is the mirror pair N_{k−1}^{(k)} = N_{k+1}^{(k)}, which the symmetry makes equal. The control would sit inside the band by construction on a quarter of its pairs.
- **Printed beside R_n, 7a's two sign sentinels:**
  - σ-blind (εᵣ for ε_c): **predicted** [0.80, 0.95]. The leading order is |σ/ωε₀|/|ε_c| = 140.43/160.64 = 0.874; 7a measured 0.873–0.893.
  - the conjugated ε_c\*: **predicted** [1.6, 2.4]. The leading order is 1.748; 7a measured 1.889–2.153, which **missed** its point prediction, so this is a range.
  - The legacy route's tag-3 fraction is **predicted** ≈ 0.5.

**Rank invariance (W2).**
- At `-n 2`, rank 0 writes `/workspace/logs/post1/step7b-n2.npz` (gitignored) with the sorted face-vertex array, F, V, Σ A_f, N, D and R_n, then `Barrier`. At `-n 4` the module loads it; if it is missing, `pytest.fail` and never recompute.
- **(W0), asserted:** the face-vertex array is bitwise equal, and F and V are equal. Mesh coordinates do not depend on the partition.
- **Printed:** max relative N, D and R_n difference, **predicted** ≤ 1e-8. N is a difference of two O(1) normal fields, so the ≈ 2e-10 solve drift `POST-8` 3a measured across widths is amplified ≈ 1/R_n.

**Windows.** Bare, foreground, Bash `timeout: 660000`, in order W1 → W2.

Basis: `MAT-15` W1 on this rung at `-n 2` took mesh + rung 61.3 s and EM 142.6 s (`20261001T123309Z_MAT-15-step1-W1.log:1924`, `:1954`). `POST-8` 3a's W1 / W2 took 272 / 231 s on the same mesh. The extraction is 12 evaluator calls on ≈ 23 000 points.
- **W1 (`-n 2`, predicted 170–260 s):**

  `scripts/testing/run_and_log.sh POST-1-step7b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/post/test_one_sided_owner_rule.py tests/validation/test_post1_one_sided_interface_birdcage.py -v -s --tb=short'"`
- **W2 (`-n 4`, predicted 140–240 s):**

  `scripts/testing/run_and_log.sh POST-1-step7b-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_post1_one_sided_interface_birdcage.py -v -s --tb=short'"`
- A smoke dev run before W1 is allowed, as 7a did (`…POST-1-step7a-dev.log`). Rule (i) still binds W1.

**Traps:**
- **The face test is per face, not per vertex.** A face with one rim, one lateral and one cap-interior vertex is an interior face. Faces seen twice are chords: drop them. Never tune the 1e-9 surface tolerance. A point found on one side only is a collision-tolerance finding: print it.
- **Point lists must be bitwise equal on every rank** (OPS-55): build from sorted vertex coordinates and never round.
- **Rank-local data:** `cell_tags.values` / `.find` are rank-local, and the candidates are owned ∩ the set (the evaluator does this). Every face sum is built from the allgathered global arrays, so do not reduce it again. The phantom-power `assemble_scalar` *is* rank-local: `allreduce` it.
- **Complex build:** w_j is numpy, so there is no UFL comparison. If any UFL form carries `SpatialCoordinate`, pin `quadrature_degree` (POST-5), and never use `ufl.max_value` or `<=` on complex operands (OPS-22).
- **ε_c** comes from `_series(...).epsilon_c`, never re-typed. The conjugate sentinel is the sign check.
- **A killed compile leaves an orphan `.c`:** run the rubric's cache sweep (`daily-review.md:426`) only when `pgrep -c python3` in the container reads 0.
- **Harness hygiene:** `-s`, no `-k`, no pipes. The harness call is the whole Bash command. After a kill, run `pgrep -c python3` first.
- **Rule (i):** an edit after W1 means W1 and W2 re-run.
- **Private data:** nothing here touches `examples/ansys_benchmarks/`. No slot opens `docs/private/` or `aed_results/`.
- **File overlap:** none. The item only imports from `post/evaluation.py` and the gate modules, and edits no helper.

**Scope.** 7b only: degree 1, 64 MHz, one fixture (the 1 g rung), interface-residual C4 identities plus records.
- No band on the R_n levels, and no convergence, SAR, B₁⁺, tuning or absolute claim.
- No 128 MHz, no degree 2, and nothing on chains F or T, κ or the sheets.

**Negative result.** Report and stop; never loosen a bound or move the band.
- **(A0) / (A1) red:** fixture drift. Add a known-issues entry, park, and mark BLOCKED (rule (d)).
- **(X0a) / (X0b) / (X0c) red:** a mesh or extraction finding. Park and mark BLOCKED.
- **(X1) or (B) red:** a one-sided-path defect. Park; `src/` is untouched, so 7a's W1 is the re-gate the review orders.
- **(P) red:** a defect in the sums. Fix the code, never `PARTITION_RTOL`, and re-run W1.
- **(G1) or (G2) red:** record both tables on the row, mark BLOCKED, and state which reading holds:
  - (G1) green, (G2) red, D pairs ≤ band ⇒ the residual estimator's quadrant floor (the `WF-6` step 3 analogue). A finer phantom rung or a review re-pose decides.
  - D pairs above the band too ⇒ a field / mesh C4 finding on this rung at 64 MHz.

  Either way the review rules.
- **(W0) red:** a partition-dependent face set, which is a finding. W1 stands; mark BLOCKED.
- **A window killed at 590 s** is not a result: park and mark BLOCKED. Never swap the fixture in-slot.
- There is no `PREDICTED RED` line: (G2)'s prediction straddles the band, so a red stops the slot.

**Status it can move.** The §7 `POST-1` row's Done-when (rule (j), read 2026-10-03) is "7b recorded on the loaded birdcage with its C4 images gated at the band its brief imports, elapsed recorded ⇒ ✅". This brief matches it: **(G1) and (G2) green at `C4_COVARIANCE_BAND`**, with every anchor green and the elapsed time recorded on the row ⇒ `POST-1` 🟡 → ✅. A (G1)/(G2) red holds 🟡 with the tables recorded as 7b's measurement, for the review to rule on.

**Outcome — BLOCKED 2026-10-03 (09:00 CDT slot, second item, executor `implementer`).** The C4 gates are red on the `-n 2` dev run, so the item stopped there. W1 and W2 were not run: there is no `PREDICTED RED` line, and any red stops the slot. The test is parked on `attempt/POST-1-step7b-20261003T141500Z` with `20261003T140719Z_POST-1-step7b-dev.log` (166 s; 28 passed, 9 failed, 1 skipped — the skip is (W0) at the writer width). Nothing is under `src/`.
- **Anchors, all green.**
  - (A0): 199 920 cells, 58 866 tag-3, sheets at k = 0…3, f = 64 MHz.
  - (A1): σ ∈ [0.5, 0.5] over 58 866 cells; |ε_c − 1| = 160.15531336764565 (7a 160.155313).
  - (X0a): F 7 776, V 3 890, 2V − 4 = 7 776, 0 chords dropped, 23 328 points.
  - (X0b): tag 3 [1, 1], tag 2 [1, 1], total [2, 2], 0 one-sided points; every one-sided read is valid.
  - (X0c): Σ A_f 2.0724797614e-02 against MAT-15's 2.072479761e-02, rel 2.0e-10. Against CAD it is −4.92e-4.
  - (X1): ratio 5.0e-15 / 5.5e-15 / 6.4e-15 / 6.0e-15.
  - (P): N ≤ 1.85e-15, D ≤ 2.15e-15.
  - (B): legacy R_n = |ε_c − 1| to ≤ 1.8e-16.
- **(G1), red at k = 0:** |N^{(k+1)} − N^{(k)}|/N^{(k)} = 5.6410 / 4.0536 / 3.4576 / 1.5543 % (band 5 %, predicted ≤ 2 %).
- **(G2), red on 8 of 16.** The table gives N pairs and D pairs, rows k, columns j = 0…3:

  | k | N j=0 | N j=1 | N j=2 | N j=3 | D j=0 | D j=1 | D j=2 | D j=3 |
  |---|---|---|---|---|---|---|---|---|
  | 0 | **11.73 %** | **9.05 %** | 0.03 % | 4.42 % | 0.92 % | 0.67 % | 0.07 % | 0.76 % |
  | 1 | **13.67 %** | **13.21 %** | 0.09 % | 0.80 % | 0.71 % | 0.32 % | 0.18 % | 1.20 % |
  | 2 | **10.38 %** | **18.75 %** | 2.90 % | 3.74 % | 0.86 % | 0.25 % | 0.51 % | 1.34 % |
  | 3 | **11.68 %** | **12.51 %** | 2.06 % | 2.75 % | 0.87 % | 0.41 % | 0.19 % | 1.73 % |

- **N_j^{(k)}** (rows k): 7.6446e-2 2.1799e-2 2.0426e-2 2.5051e-2 / 2.3944e-2 6.7480e-2 2.3771e-2 2.0420e-2 / 2.0256e-2 2.0670e-2 7.6394e-2 2.3791e-2 / 2.4681e-2 1.8155e-2 2.4545e-2 7.8611e-2.
- **D_j^{(k)}**: 6.1460e-2 3.0125e-2 1.9447e-2 3.0178e-2 / 2.9950e-2 6.2026e-2 3.0327e-2 1.9434e-2 / 1.9202e-2 3.0162e-2 6.2222e-2 3.0272e-2 / 2.9867e-2 1.9368e-2 3.0236e-2 6.2540e-2.
- **Quadrant R_n:** 1.115 0.851 1.025 0.911 / 0.894 1.043 0.885 1.025 / 1.027 0.828 1.108 0.887 / 0.909 0.968 0.901 1.121.
- **Records.**
  - R_n per drive: 1.0089 / 0.9782 / 0.9974 / 1.0139. All four miss the predicted [0.05, 0.6], so the separation is 158–164× against the predicted ≥ 250×.
  - Phantom power is P1–P4 = IMP-3 W2 to ≤ 3.2e-10.
  - The single-drive mirror reads 12.98 / 0.72 / 15.10 / 0.55 %; k = 0 and 2 miss the predicted ≤ band.
  - Elapsed: build 39.1 s, solves 20.9 / 20.4 / 19.9 / 20.5 s, extraction 40.6 s, fixture 163.0 s.
- **Controls (printed).**
  - The mis-rotated mean is 65.5 / 75.0 / 98.2 / 84.7 %, inside the predicted [10, 150] %. Its ratio to the (G2) mean is 10.4–11.7× (predicted ≥ 3), under the ceilings 274 / 230 / 277 / 333 %.
  - σ-blind reads 1.008 / 1.001 / 1.005 / 1.009, which misses the predicted [0.80, 0.95].
  - ε_c\* reads 2.315 / 2.297 / 2.316 / 2.323, inside [1.6, 2.4].
  - The legacy tag-3 fraction is 0.516.
  - Nearest-face pointwise: median 69.5 / 70.0 / 68.8 %, p90 ≈ 400–437 %, match distance ≈ 0.21 h (predicted median ≥ 10 %).
- **Reading, for the review.** The brief's two readings both assume a case that did not occur: (G1) is also red, at one pair.
  - The D pairs are ≤ 1.73 %, so the field itself is C4 at this rung. This is closer to the residual estimator's quadrant floor than to a field or mesh C4 finding.
  - The red pairs sit in the fixed absolute columns j = 0, 1 for every k, not at a drive-relative position. That pattern points at a mesh feature in fixed space, quadrants 0–2, rather than at the drive. The mirror split (k = 0, 2 red; k = 1, 3 green) is consistent with this.
  - R_n ≈ 1, roughly 5× the predicted level, says the normal ε̃E jump is not resolved on this rung.
  - The review rules on what happens next (a finer phantom rung or a re-pose). (W0) and the `-n 4` record were not measured.

### Ruling — 2026-10-03 10:30 daily review (item 203's red)

**Reading: the residual estimator's mesh floor in fixed space, not a field C4 finding.** The review ruled from the tables above and from the dev log. That log is now on `main` (`20261003T140719Z_POST-1-step7b-dev.log`), landed by path from the attempt branch; its test-results row was appended by hand.
- **The field is C4 on this rung.** The sixteen D pairs (the \|E⁺·n\|² quadrant integrals) are all ≤ 1.73 %, and phantom power reproduces `IMP-3` W2 to ≤ 3.2e-10.
- **The residual is not resolved.** R_n ≈ 1.0 per drive, about 5× the brief's predicted [0.05, 0.6], and the σ-blind sentinel saturates (1.00–1.01 against a predicted [0.80, 0.95]). At R_n ≈ 1, N measures the local discretisation error of the normal trace, not the interface physics.
- **The red sits in fixed space.** Every red (G2) pair is in absolute columns j = 0 / 1, for every k. Each drive's own-quadrant N, the dominant entry, reads 7.64 / 6.75 / 7.64 / 7.86e-2: quadrant 1 is ≈ 12 % lower whatever the drive.
  - On a mesh that is C4-congruent only on the port sheets (`io/mesh.py:3748`), a discretisation error inherits the mesh's quadrant differences. A C4 identity on it holds only as far as the mesh is congruent.
  - (G1) fails at k = 0 for the same reason: the whole-interface sum is dominated by the driven quadrant.
- **What would not obviously fix it.** A finer phantom rung lowers R_n, but it does not make the mesh congruent, so the fixed-space contrast need not fall below 5 %. This is a prediction, not a measurement.

**Disposition.**
- Item 203 stays BLOCKED. The parked module stays on `attempt/POST-1-step7b-20261003T141500Z`; every anchor in it is green and reusable. `POST-1` holds 🟡.
- This is attempt 1 of the step-7b family: 7a registered its record, so the count restarted there.
- No band moves, and (G1) / (G2) are not converted to prints.

**The re-pose is the 10-03 weekly's.** These are the candidates it weighs:
- **(i) A C4-congruent phantom mesh,** so that the identity is exact on the mesh. This is a new geometry capability: a quadrant-built mesh like `TH-16`'s cutter is D2, not C4, so it would need an eighth-cut or a rotated-copy build.
- **(ii) A finer rung in the `xl` tier,** as a record of whether the contrast falls.
- **(iii) Banking the measured negative** and closing 7b as a record. That needs a Done-when change, which is the weekly's.

## §7 row moved verbatim 2026-10-03 (rotate_plan_archive.py chunks)

 Interface-aware field extraction reliability | 🟡 Steps 1–6 ✅ (2026-08-04 … 08-06: ghost-cell and guardrail rank-safety fixed, `prefer_interior=False` the production default, CSV/stats parity gated). The 09-30 weekly's last step is split into 7a, the one-sided interface extraction on `TH-10`'s sphere against the series (§9 item 187), and 7b, the same residuals on the loaded birdcage as a C4-gated record (not written; serial on 7a). **7a done 2026-10-03** (item 187): one-sided owner rule exact (F = 2V − 4, trace ≤ 5.8e-15), R_n / e⁻ / e⁺ strictly decreasing to 0.178 / 0.020 / 0.133 (a record); 7b is queued as §9 item 203, with C4 gated on whole-interface and quadrant integrals, since the mesh is not C4-congruent pointwise. **7b measured 2026-10-03, red, held 🟡** (item 203 BLOCKED; `-n 2` dev run, 166 s; elapsed: build 39 s, four solves ≈ 20 s each, extraction 41 s). Every anchor is green, but (G1) k = 0 reads 5.64 % against the 5 % band, (G2) is red on 8 of 16 quadrant pairs (max 18.75 %), and the D pairs are ≤ 1.73 %. R_n ≈ 1.0 per drive. Both 4 × 4 tables are in the chunk file. **Ruled 2026-10-03 10:30 review:** the estimator's mesh floor in fixed space (the red pairs sit in absolute quadrants j = 0 / 1 for every drive), not a field C4 finding. The re-pose is the weekly's (attempt 1 of the 7b family). **Done-when:** *(written 2026-10-02 18:00 review, rule (j), from the 09-30 weekly's lane-1 target; the weekly ratifies)* 7a asserted — a partition-free interface set (one tag-1 and one tag-2 cell per point, F = 2V − 4), the N1curl tangential trace continuous to 1e-10, and the normal ε̃E-jump residual and both one-sided fields against `LossySphereSeries` decreasing over three rungs, with levels recorded ⇒ stays 🟡; 7b recorded on the loaded birdcage with its C4 images gated at the band its brief imports, elapsed recorded ⇒ ✅. *History: `docs/planning/chunks/POST-1.md`; `docs/planning/plan-archive.md` § "§7 POST-1 adjudication + steps 1–6".* 

### §9 item 245 — 2026-10-05 (lane 1)

*Drafted by a general-purpose agent for the 2026-10-05 03:00 daily review; the review reviewed it.*

This is 7b′, the 10-03 weekly's re-pose of 7b (`PROJECT_PLAN.md:4451–4456`): "the same one-sided residuals on `TH-16`'s **mirror-built** F-small mesh, where the D2 images exist cell for cell, gating the D2 image pairs at the brief's imported band; the C4 statement is dropped, not loosened". It is **attempt 2 of 4** in the 7b family (attempt 1 is item 203, ruled above).
- **Executor:** implementer. **Tier:** standard (W1 predicted 75–130 s). `-n 2` (W1) and `-n 4` (W2), complex only. There is no real-build window, because nothing under `src/` changes.
- **Class:** lane-1 revalidation (measured median 8 min). **Cost:** ≈ 8 measured / ≈ 30 predicted min. **Depends on:** nothing (7a is done; `TH-16` is ✅).

**Fixture: `TH-16` step 3c's mirror-built full mesh** (`tests/validation/test_th16_quarter_birdcage.py`). The module imports that file's helpers and edits nothing in it.
- **Mesh.** Build the quarter with `_build(True, symmetry_quadrant=True)` (`tests/mesh/test_birdcage_port_sheets.py:80`): 36 034 cells, mesh 6.90 s, rung 8.09 s at `-n 2` (`20260930T093451Z_TH-16-step3c-W1.log:934`).
- **Mirror.** Gather, mirror and distribute exactly as `full` does (`test_th16_quarter_birdcage.py:514–598`): `_gather_owned` (`:234`); `_mirror_by_index` across axis 1, then axis 0 (`:255`, `:518–519`); `_distribute` (`test_th16_symmetry_plane_driven.py:145`); cell tags through `distribute_entity_data` + `meshtags_from_entities` (`:528–539`). This gives 144 136 cells and 26 253 vertices, built in 1.03 s; per-tag cells 1:50 624, 2:87 016, 3:1 256, 102:2 704, 201:2 536 (`…W1.log:990–991`). The phantom is tag 3, with **1 256 cells** (≈ 1.1 cm).
- **Port sheets: the whole cut-plane box face, not the narrowed sheet (a deviation from `TH-16`'s letter, for the review).**
  - `TH-16` narrows each quarter sheet to `GATED_WIDTH_FRACTION` 0.5 (`:318–327`). That is a sub-sheet triangulated out of a box, which the 2026-10-02 binding ruling forbids in new fixtures.
  - The exact-strip knob refuses `symmetry_quadrant` (`src/fem_em_solver/io/mesh.py:4266–4271`).
  - So the quarter sheets are the **un-narrowed** faces `_owned_exterior_on_plane(msh, cells, SHEET_PLANE_AXIS[p], KEPT_HALF[p])` (`:154`, `:306`). Each is the whole CAD face of the port box on the cut plane, 50 / 56 facets (`SHEET_FACET_RECORD`, `:109`). The box is 10 mm wide against the 12 mm leg (`test_birdcage_port_tags.py:24,30`). This whole-mid-plane sheet has been solved before: it is `PORT` step 2's own sheet, the f = 1.0 control (`test_port_lumped_narrowed_sheet.py:127–131`).
  - Width is w = A/h with κ off. Copy `:304–317` and `:328–341` and omit `:318–327`.
  - The mesh is unchanged, and the D2 symmetry holds with either sheet. This item reads no port quantity.
- **Full-mesh sheets** 211–214 come from the midpoint match at `:565–598`: P1 = P1q; P2 = P2q; P3 = P1q flipped on axis 0; P4 = P2q flipped on axis 1 (`FULL_SHEETS`, `:134`; `SHEET_PLANE_AXIS` `:140`).
- **Solve.** The problem is as at `:603–617` (all exterior PEC; conductor `SIGMA_WIRE_S_PER_M`; tag 3 at `SALINE_SIGMA` / `SALINE_EPSILON_R`). Two changes: `frequency_hz=FREQUENCY_64_HZ`, and the sweep is `run_n_port_sparameter_sweep(problem, defs, lumped_sheet_ports=specs, lumped_sheet_facet_tags=tags_s, degree=1, keep_fields=True)` (`ports/sparameters.py:310,325`). The specs are as at `:618–641`, with `interior=True` and `FULL_Z0_OHM`. Drive k's field is `result.fields[f"P{k+1}"].e_complex`, and its σ is `.sigma_field`.
- **Code by path.** Land the parked module with `git checkout attempt/POST-1-step7b-20261003T141500Z -- tests/validation/test_post1_one_sided_interface_birdcage.py`. It is code that `main` has never had. Never check out its log or `test-results.md` by path.
  - **Re-target in place.** Keep `_FACES`, `_BARY`, `_w`, `_interface_faces`, `_side_counts`, the per-drive extraction loop and its sentinels, and the (X0a) / (X0b) / (X1) / (P) / (B) / (W0) tests.
  - **Delete:** `build_four_port_sweep`, the 1 g records in (A0), (X0c) with `MAT15_SURFACE_AREA_M2`, `IMP3_W2_PHANTOM_POWER_W`, `_rotz` and the nearest-face reading, and the (G1) / (G2) C4 tests.
  - **Docstring:** 7b′. Say that the 7b version is at `cbd32287`.
  - **New npz path:** `/workspace/logs/post1/step7b2-n2.npz`.

**The D2 image map (derived; mesh and drives).**
- **Port positions.** The quarter keeps x, y ≥ 0. P1q lies in y = 0, so it is at θ = 0; P2q lies in x = 0, so it is at θ = 90°. P3 is the axis-0 flip of P1 (θ = 180°) and P4 the axis-1 flip of P2 (θ = 270°). So drive k = (port number − 1) sits at θ_k = kπ/2, the same centring as `_w`'s θ_j = jπ/2.
- **The three maps** (indices mod 4):
  - σ_x, the mirror x → −x: P1 ↔ P3, P2 and P4 fixed; k ↦ 2 − k, j ↦ 2 − j, since c_j(M_x p) = c_{2−j}(p) in `_w`.
  - σ_y, the mirror y → −y: P2 ↔ P4, P1 and P3 fixed; k ↦ −k, j ↦ −j.
  - σ_π = σ_xσ_y: k ↦ k + 2, j ↦ j + 2.
- **Why the residual is invariant.** Every drive is ẑ-directed (`:375`, `:637`), and ẑ is invariant under both mirrors, so the image of drive k is drive σk with the same sign. E is a polar vector, so the image field is M·E(Mx). The image face's outward normal is M n_f, so E·n is invariant. Hence |r|² at Mx under drive σk equals |r|² at x under drive k.
- **The cell-for-cell premise is `TH-16`'s measurement.** On this mesh the S-level D2 pairs read S₁₁/S₃₃ 1.491e-13, S₂₂/S₄₄ 9.769e-13 and adjacent 1.663e-13 (`…W1.log:946–947,1021–1022`).
- **No face is its own image.** A self-image would need three collinear vertices, so every orbit has 4 faces and F ≡ 0 (mod 4).

**Imported, never restated.**
- **The band:** `C4_COVARIANCE_BAND` = 5 % (`test_birdcage_b1_plus_map.py:120`), the band item 203 imported (its brief, above). It is used here on D2 pairs only.
- From `TH-16`: `MIDPOINT_MATCH_TOL_M` (`:132`), `PLANE_TOL`, `QUARTER_CELL_RECORD`, `FULL_CELL_RECORD`, `SHEET_FACET_RECORD`, `SHEET_AREA_RTOL` and `FULL_SHEETS` / `FULL_SHEET_TAG` / `SHEET_TAG` / `KEPT_HALF` / `SHEET_PLANE_AXIS` / `SHEET_TRANSVERSE_AXIS`.
- `EXACT_IDENTITY_RTOL` (`test_birdcage_symmetry_quadrant.py:47`).
- `LEG_D0_REPRODUCTION_BAND` and `PASSIVITY_SIGMA_TOLERANCE` (`test_port_birdcage_four_port.py:165`).
- `_class_means` is not needed. Copy `TH-16`'s `d2` dict lines (`:659–664`) verbatim.
- Everything else as in the parked module: `_series`, `FREQUENCY_64_HZ`, `SALINE_SIGMA`, `PHANTOM_CELL_TAG`, `PHANTOM_RADIUS` / `PHANTOM_HEIGHT` (the same 0.03 / 0.08 m that `_build` passes), `PARTITION_EPS_M` / `PARTITION_RTOL`, `_phantom_sigma_range` / `PHANTOM_SIGMA_RTOL`, `CAD_AREA_M2` (printed). For `_phantom_sigma_range`, pass `{"mesh": msh, "cell_tags": cells}` as the sweep argument.

**Asserted.**
- **(M0) The mesh is `TH-16`'s,** as `test_mirror_build_conforms` (`:758–768`) with the sheet record swapped:
  - quarter cells == `QUARTER_CELL_RECORD`;
  - full cells == `FULL_CELL_RECORD`; vertices == builder count;
  - per-tag cells ×4 exactly; |V_f/(4V_q) − 1| ≤ `EXACT_IDENTITY_RTOL`; zero owned exterior facets on either plane;
  - un-narrowed quarter sheets == `SHEET_FACET_RECORD`; full sheets 50 / 56 / 50 / 56 owned facets, with `dS` areas equal to the quarter's at `SHEET_AREA_RTOL`;
  - sheet k's mean midpoint on the predicted side: P1 x > 0, P2 y > 0, P3 x < 0, P4 y < 0.
- **(A1)** The tag-3 σ is uniform at `SALINE_SIGMA` to `PHANTOM_SIGMA_RTOL` over 1 256 cells (reduced). f = `FREQUENCY_64_HZ`.
- **(I3) The drives are D2 images.** `S_f`'s pairs S₁₁/S₃₃, S₂₂/S₄₄ and the adjacent spread are each ≤ `LEG_D0_REPRODUCTION_BAND`. **Predicted** ≤ 1e-11 (1.5e-13 … 9.8e-13 at 10 MHz with the narrowed sheet, above). (I1) reciprocity and (I2) σ_max are printed.
- **(X0a)** F = 2V − 4 exactly. **(X0b)** One owned tag-3 cell, one tag-2 cell and two in total per point (`Allreduce`), and every one-sided read is valid. **(X1)** max |n_f×(E⁻ − E⁺)| ≤ 1e-10·max|E⁺| per drive. **(P)** The partition identity at `PARTITION_RTOL`. **(B)** The legacy route gives |ε_c − 1| to 1e-12. All five are item 203's, unchanged.
- **(M1) The face set is D2-closed.** For σ ∈ {σ_x, σ_y}, every face's image centroid has exactly one face centroid within `MIDPOINT_MATCH_TOL_M`.
  - The map is a fixed-point-free involution.
  - Matched areas agree to `EXACT_IDENTITY_RTOL`, and so do matched normals (n_g = M n_f, absolute on unit vectors).
  - The matching is brute force on every rank: the arrays are allgathered and identical, and F ≈ 600.
- **(G1′) Whole-interface D2 pairs.** |N^{(σk)} − N^{(k)}|/N^{(k)} ≤ `C4_COVARIANCE_BAND` for (σ_x, k = 0, 2) and (σ_y, k = 1, 3), 4 readings (σ_π gives the same pairs).
- **(G2′) Quadrant D2 pairs.** For each σ ∈ {σ_x, σ_y, σ_π} and each (k, j) with (σk, σj) ≠ (k, j): |N_{σj}^{(σk)} − N_j^{(k)}|/N_j^{(k)} ≤ `C4_COVARIANCE_BAND`. That is 12 + 12 + 16 = 40 parametrised readings, each naming (σ, k, j). They include the single-drive mirrors (σ_y on k = 0 is N₁⁽⁰⁾ against N₃⁽⁰⁾), which read 12.98 / 0.72 / 15.10 / 0.55 % on the non-mirror 1 g mesh (`20261003T140719Z_POST-1-step7b-dev.log:1897–1900`; a different mesh, so context only).
- **(G3′) Face-level D2 pairs (pointwise).** For σ ∈ {σ_x, σ_y} and each k: max_f |ρ_{σf}^{(σk)} − ρ_f^{(k)}|/ρ_f^{(k)} ≤ `C4_COVARIANCE_BAND`, with ρ the face mean of |r|². That is 8 readings. Also print min ρ_f / mean ρ.
- **Predictions** for (G1′) – (G3′): ≤ 1e-9 for the integrals and ≤ 1e-7 face-wise. Basis: the D2 identity is exact on this mesh up to solve round-off. `TH-16` (I3) measured ≤ 9.8e-13, and `-n 2` in-process solves differ by ≈ 2e-12 (`MAT-5` W1, ruled 2026-10-05). |ε_c| ≈ 160 amplifies E⁻'s noise in r.
  - **The band is ≈ 1e7 above the prediction. Reporting rule (not a band):** any reading above 1e-8 is journaled as unexplained, even when green.

**Records (printed, no band).**
- R_n per drive, **predicted** [0.8, 1.5]: 7b read ≈ 1.0 at h ≈ 2.5 mm, and this phantom is coarser.
- The N_j^{(k)}, D_j^{(k)} and quadrant R_n tables, and the same D2 pairs on D.
- F (**predicted** [400, 900], ≡ 0 mod 4), V and the number of points.
- Σ A_f against `CAD_AREA_M2`.
- Phantom power per drive through `mean_sar` (reduced), with its σ pairs.
- `S_f`.
- Elapsed: quarter build, mirror build, sweep, extraction, window.
- **The dropped C4 statement, as a record only:** item 203's (G1) / (G2) C4 pairings on this mesh. They are never gated and never read as a field finding: this is the mesh's own C4 contrast. **Predicted** 0–20 %.

**Negative controls (rule (e)).**
- **(B) is asserted** (algebraic: both sides read one cell). 7b measured ≤ 1.775e-16 (`…POST-1-step7b-dev.log:1861–1873`). Its separation |ε_c − 1|/R_n is **predicted** 100–200× and printed. **Ceiling:** 160.155/R_n by construction, so ≈ 160× at R_n ≈ 1.
- **(N1) The half-mapped pairing, predicted and printed:** map the quadrant but not the drive, only on drives σ moves.
  - Under σ_x: k ∈ {0, 2}, j ∈ {0, 2}. Under σ_y: k ∈ {1, 3}, j ∈ {1, 3}.
  - The reading is |N_{σj}^{(k)} − N_j^{(k)}|/N_j^{(k)}, the driven quadrant against the opposite one.
  - **Predicted** [50 %, 400 %], from 7b's own-to-opposite ratios on the 1 g mesh: 0.733 / 2.743 at k = 0 through 0.769 / 3.330 at k = 3 (N table above).
  - **Ceiling** per drive: max_j N_j^{(k)}/min_j N_j^{(k)} − 1, printed. The reading cannot exceed it.
  - Separation over (G2′)'s max: **predicted** ≥ 1e6×, printed.
  - **Trap:** never take the half-mapped pairing on a drive σ fixes (σ_y on k = 0). That *is* the identity, and it sits inside the band by construction. This is the D2 analogue of item 203's wrong-sense trap.
- **(N2) The face-level half-mapped pairing** ρ_{σf}^{(k)} against ρ_f^{(k)}, on drives σ moves: median and p90 printed, **predicted** median ≥ 10 % (7b's nearest-face C4 median read 68.8–70.0 %).
- Item 203's sentinels, printed:
  - σ-blind, **predicted** [0.95, 1.10] (7b 1.001–1.009, `:1862–1874`);
  - ε_c\*, **predicted** [1.6, 2.6] (7b 2.297–2.323);
  - the legacy tag-3 fraction, ≈ 0.5.

**Rank invariance (W2, `-n 4`), as item 203.**
- **(W0), asserted:** the sorted face-vertex array is bitwise equal to the `-n 2` npz, F and V are equal, and (M0) holds at `-n 4`.
- **Printed:** max relative N / D / R_n difference, **predicted** ≤ 1e-8. If the npz is missing, `pytest.fail`; never recompute.
- **Note:** every `TH-16` window ran at `-n 2` (`…TH-16-step3a-W1/W2/W3`, `…3c-W1/W2`), so the quarter mesh at `-n 4` is unmeasured. It is **predicted** identical, because gmsh meshes on rank 0.

**Windows.** Bare, whole, repo-relative harness calls; foreground; Bash `timeout: 660000`; `-s`; no `-k`; no pipes.

Pricing: `TH-16` 3c W1 was 85 s at `-n 2` with three extra quarter sweeps (`…W1.log:1174`). Here: quarter ≈ 8 s + mirror ≈ 1 s + the four-drive sweep ≈ 20 s (20.30 s at 10 MHz, `:995`) + extraction ≈ 10 s (item 203 took 40.6 s for 23 328 points; here ≈ 1 800) + JIT.
- **First,** `pgrep -c python3` in the container must read 0. Run the rubric's orphan sweep (`daily-review.md:426`) only if a previous window was killed.
- **A smoke dev run** (`POST-1-step7b2-dev`) is allowed before W1. Rule (i) binds W1.
- **W1 (`-n 2`, predicted 75–130 s):**

  `scripts/testing/run_and_log.sh POST-1-step7b2-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/post/test_one_sided_owner_rule.py tests/validation/test_post1_one_sided_interface_birdcage.py -v -s --tb=short'"`
- **W2 (`-n 4`, predicted 60–120 s):**

  `scripts/testing/run_and_log.sh POST-1-step7b2-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 180 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_post1_one_sided_interface_birdcage.py -v -s --tb=short'"`

**Traps.**
- **Rank-local data.** `cell_tags.values` / `.find` and `assemble_scalar` are rank-local. Face sums are built from allgathered arrays, so never reduce them twice. `mean_sar` power and the σ range are reduced.
- **Point evaluation** goes only through `evaluate_vector_field_parallel` (`owner="lowest_original_cell"`, `candidate_cells=…`), never `f.eval`. The owner rule raises on a `create_submesh` mesh (known-issues 2026-10-04, `ANS-9` 1a). The mirrored mesh comes from `create_mesh`, so it is fine. Never build the interface or evaluate on a submesh.
  - The mirrored mesh's `original_cell_index` follows the rank-order gather, so it is **width-dependent**. That is harmless only because (X0b) asserts one candidate per side. Never read cell indices across widths.
- **Pair faces, not points.** The three barycentric points follow the sorted vertices, and a mirror permutes them. Gate on ρ_f (the face mean) and on quadrant sums, never per-point by index.
  - On-plane vertices are reused, not negated, so their coordinates may be ≈ 1e-17 rather than 0. Match with `MIDPOINT_MATCH_TOL_M`, never bitwise.
  - `_w`'s `np.cos(π/2)` is 6e-17, so the w-identity holds to ≈ 1e-16, not bitwise.
- **Complex:** every face and quadrant sum is numpy, with no UFL comparison. If any UFL form carries `SpatialCoordinate`, pin `quadrature_degree` (POST-5). Never apply `ufl.max_value` or `<=` to complex operands (OPS-22). ε_c comes from `_series(FREQUENCY_64_HZ).epsilon_c`, never re-typed.
- **Known-issues 9:** cell tags only through `distribute_entity_data`; hoist `create_connectivity(fdim, tdim)` and `create_entity_permutations()` before any facet form (`:541–543`). Copy, don't re-derive.
- **Parked code by path, record files never.** After landing, `git diff --stat` shows the module only.
- **Orphan `.c`:** sweep only with `pgrep -c python3` reading 0. After a kill, run `pgrep -c python3` first.
- **Rule (i):** an edit after W1 re-runs W1 and W2.
- **Private data:** nothing here touches `examples/ansys_benchmarks/`. No slot opens `docs/private/` or `aed_results/`.
- **File overlap:** the item imports from `test_th16_quarter_birdcage.py`, `test_birdcage_b1_plus_map.py` and `post/evaluation.py`, and edits none of them.

**Scope.** 7b′ only: degree 1, 64 MHz, one fixture (`TH-16`'s mirror-built F-small mesh, 1 256-cell phantom, whole-face sheets), interface-residual D2 identities plus records.
- No band on R_n, and no convergence.
- No C4 statement: dropped by the weekly and printed only.
- No SAR, B₁⁺, S, tuning or absolute claim, and nothing about the sheet or κ.
- No 128 MHz, no degree 2, no 1 g quadrant rung (`phantom_resolution` with `symmetry_quadrant` is unmeasured), and nothing on chains F / T.

**Negative result.** Report and stop. Never loosen a bound or move the band. Park on `attempt/POST-1-step7b2-*`, mark the item BLOCKED (rule (d)) and journal it in attempts.md. This is attempt 2 of 4: two attempts remain before the family goes to the weekly to bank (item 203 ruling, candidate (iii)).
- **(M0) red:** `TH-16`'s fixture has drifted. Add a known-issues entry; `TH-16`'s W1 is the re-gate the review orders.
- **(M1) or (I3) red:** the mesh or the drives are not D2 as built. This is a fixture finding, not a POST-1 one.
- **(X0a) / (X0b) / (X1) / (B) red:** a mesh or one-sided-path finding, as item 203.
- **(P) red:** fix the sums, never `PARTITION_RTOL`, and re-run W1.
- **(G1′) / (G2′) / (G3′) red with (M0), (M1) and (I3) green:** the field and the mesh are D2, yet the one-sided extraction is not mirror-consistent. That is a **POST-1 extraction finding** (the owner rule, candidate restriction or face selection). Record every pair table; the review rules.
- **(W0) red:** a width-dependent face set or quarter mesh, which is a finding. W1 stands; mark BLOCKED.
- **A window killed at 180 s** is not a result. Park BLOCKED with the printed clocks; never lengthen the timeout and never swap the fixture in-slot.
- There is no `PREDICTED RED` line: a red anywhere stops the slot.

**Status it can move.** The §7 `POST-1` Done-when (rule (j), read 2026-10-05) is "7b recorded on the loaded birdcage with its C4 images gated at the band its brief imports, elapsed recorded ⇒ ✅".
- **This brief matches it except for one word.** The fixture is a loaded birdcage (the saline phantom, tag 3), the band is the one item 203 imported, and the elapsed time is recorded. But the gated images are **D2 image pairs**, per the weekly's re-pose.
- **The review must amend the row in its queueing commit:** "with its C4 images gated" → "with its D2 image pairs gated (on `TH-16`'s mirror-built F-small mesh; the C4 statement dropped by the 10-03 weekly)". If the row is not amended, the executor closes on the row: it holds 🟡 and flags the mismatch.
- **With the row amended:** (M0), (M1), (I3), (G1′), (G2′) and (G3′) green in W1, (W0) green in W2, and every anchor green ⇒ `POST-1` 🟡 → ✅, with the elapsed time on the row.
- A red holds 🟡, with the tables recorded as attempt 2's measurement.


**Outcome — DONE 2026-10-05 (07:30 CDT slot, executor `implementer`). Every gate green; `POST-1` 🟡 → ✅.** The parked module was landed by path from `attempt/POST-1-step7b-20261003T141500Z` and re-targeted in place (`tests/validation/test_post1_one_sided_interface_birdcage.py`; nothing under `src/`, no imported file edited).
- **Windows** (all Status 0): dev `20261005T123358Z_POST-1-step7b2-dev.log` (36 s, 71 passed, 1 skipped — (W0) at the writer width); W1 `20261005T123446Z_POST-1-step7b2-W1.log` (`-n 2`, 88 s, 86 passed, 1 skipped); W2 `20261005T123627Z_POST-1-step7b2-W2.log` (`-n 4`, 83 s, 84 passed). No edit between the dev run and W1.
- **(M0):** quarter 36 034 cells; full 144 136 cells, 26 253 vertices (= builder); per-tag ×4 exactly (1:50 624, 2:87 016, 3:1 256, 102:2 704, 201:2 536); |V_f/(4V_q) − 1| 6.9e-15 (`-n 2`) / 1.1e-15 (`-n 4`); no exterior facet on either plane. Un-narrowed sheets 50 / 56 / 50 / 56; dS areas against the quarter ≤ 4.4e-16; mean midpoints P1 x = +6.99e-2, P2 y = +7.03e-2, P3 x = −6.99e-2, P4 y = −7.03e-2 m. Each whole-face sheet reads A = 1.12e-4 m², h = 8 mm, w = A/h = 14.0 mm = w_bbox (the brief's "10 mm" box width is not what the mesh carries; no gate reads it).
- **(A1):** σ ∈ [0.5, 0.5] over 1 256 cells. **(I3):** S₁₁/S₃₃ 1.94e-13, S₂₂/S₄₄ 3.03e-14, adjacent 5.23e-14 (band 1e-9, predicted ≤ 1e-11). (I1) reciprocity 1.59e-14; (I2) σ_max 0.999745.
- **Anchors:** (X0a) F 296, V 150, 2V − 4 = 296, 0 chords dropped, 888 points; (X0b) tag 3 [1, 1], tag 2 [1, 1], total [2, 2]; (X1) ratio ≤ 1.34e-15; (P) N ≤ 1.69e-15, D ≤ 1.74e-15; (B) legacy R_n = |ε_c − 1| to ≤ 1.8e-16.
- **(M1):** σ_x and σ_y each one hit per face (worst distance 9.8e-18 m), fixed-point-free involutions, matched areas ≤ 2.2e-15, normals ≤ 2.8e-16.
- **(G1′)** 2.17e-13 (σ_x, k = 0, 2), 1.95e-13 (σ_y, k = 1, 3). **(G2′)** all 40 readings ≤ 1.63e-12 on N (D pairs ≤ 1.92e-12). **(G3′)** face-level max 5.04e-11 (`-n 2`); the `-n 4` worst D2 reading is 6.39e-11. All under the 1e-8 reporting level, so nothing is journaled as unexplained. Band 5 %.
- **(W0):** sorted face-vertex array bitwise equal at `-n 4`; F, V equal; max relative N 7.8e-13, D 3.8e-12, R_n 2.3e-12 (predicted ≤ 1e-8).
- **Records.**
  - R_n per drive 3.947 / 3.312 / 3.947 / 3.312. N_j^{(k)} rows: 0.6510 0.1449 0.1791 0.1449 / 0.1246 0.3839 0.1246 0.1547 / 0.1791 0.1449 0.6510 0.1449 / 0.1246 0.1547 0.1246 0.3839. D_j^{(k)} rows: 3.245e-2 1.489e-2 9.648e-3 1.489e-2 / 1.505e-2 3.229e-2 1.505e-2 9.415e-3 / (k = 2, 3 the D2 images). Quadrant R_n 2.88–4.48.
  - Σ A_f 2.0504e-2 against CAD 2.0735e-2 (−1.11 %). Phantom power 1.5789e-6 / 1.5621e-6 W (P1 / P2; P3 / P4 equal to ≤ 3.9e-14).
  - Dropped C4 record (never gated): G1 29.66 / 42.17 %; G2 N pairs 13.6–69.6 %, D pairs 0.47–2.47 %.
  - Elapsed (W1): quarter 8.1 s (gmsh 6.8 s), mirror 1.0 s, sweep 11.2 s, extraction 2.7 s, fixture 23.1 s; W2 fixture 18.5 s.
- **Controls.** (N1) half-mapped 72.5 / 263.5 % (σ_x, ceiling 349.4 %) and 59.7 / 148.2 % (σ_y, ceiling 208.2 %), inside the predicted [50, 400] %, separation ≥ 3.7e11× (predicted ≥ 1e6). (N2) face median 84.2 / 91.7 %, p90 771 / 1 210 % (predicted median ≥ 10 %). Legacy tag-3 fraction 0.497.
- **Printed predictions that missed (records, no status effect, no assertion touched):** F 296 against [400, 900] (the phantom is coarser than priced); R_n 3.31–3.95 against [0.8, 1.5]; (B) separation 40.6 / 48.4× against 100–200; σ-blind 1.89–2.18 against [0.95, 1.10]; ε_c\* 4.04–4.53 against [1.6, 2.6]; the dropped C4 contrast up to 69.6 % against 0–20 %. All follow from R_n ≈ 3–4 on a 1 256-cell phantom (h ≈ 1.1 cm): the normal-jump residual is not resolved here either, so these levels describe the discretisation, not the interface physics. The D2 identity holds at round-off regardless, which is what 7b′ gates.
