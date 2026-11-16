# GEO-27 — Phantom-resolution cost probe for the 1 g SAR rung

 **Phantom-resolution cost probe for the 1 g SAR rung** — opened 2026-09-07 10:30 review from `MAT-4` step 4's printed 1 g column (worst C4 pair 6.8383% on the 7.5 mm rung, `20260907T003548Z_MAT-4-step4.log:1934–1947`: a 6.2 mm ball spanning ≈1.7 cells is a mesh statement, not an operator one, and the route to a 1 g gate is `h`). Executor **`mesh-probe`**, measurement-only, asserts nothing — 🧪 by the §3 rule until a gate consumes it. A ladder on `MeshGenerator.birdcage_port_domain(phantom_resolution=…)` exactly as `tests/mesh/test_birdcage_phantom_resolution.py` calls it (F-small defaults, real build, no solve): rungs **0.0075** (the record rung — the control that must reproduce 120 499 cells / 2 746 phantom cells, `20260907T140658Z_EX-53.log:1875`), **0.005**, **0.00375**, **0.0025**; per rung print total cells, phantom cells (tag 3), mesh wall time, the phantom share of the cells, and the cells across the two C95.3 balls, `2a_1g/h` and `2a_10g/h` (`a_1g` = 6.2 mm, `a_10g` = 13.4 mm at ρ = 1000 kg/m³, the gate's own `averaging_ball_radius`); the table lands in this row and the script as `scripts/probes/geo27_phantom_resolution_probe.py`. **Stop rule:** a rung above 300 k cells or 120 s to mesh is recorded and the ladder stops there — no finer rung is attempted. **Tier / cost:** the 120 499-cell mesh builds in ≈23–29 s (`20260906T213913Z_TH-15.log:2628`, `…093741Z_TH-15.log` solid 28.70 s); four rungs ≈ 200 s class ⇒ `-n 2`, `timeout -k 30 600`, heavy by ceiling, standard by expectation. **Traps:** the meshing happens on rank 0 and the counts are `size_global` after distribution — print the global count, never `cell_tags.values` lengths; the `PHANTOM_GROWTH_BAND` in the resolution module is a gate on the 0.0075 rung only and is not re-asserted here. **Question the probe answers:** the first rung with ≥ 4 cells across the 1 g ball and its cell count — the input to a `MAT-4` step 5 (the 1 g column gated at that rung, or ruled ungateable on this box) that a review scopes from the table, never the slot that runs it. **Negative result:** the 0.0075 control not reproducing 120 499 is a mesher-drift finding (`OPS-18` class) — record, stop. **Executed 2026-09-08 (22:30 slot, `mesh-probe`, one window `-n 2`, `timeout -k 30 600`, elapsed 132 s, `20260908T033218Z_GEO-27.log`).** All four rungs meshed; the ladder ran to completion. **The control REPRODUCED to the integer** — 120 499 cells / 2 746 tag-3 cells against the record (`:1803–1810`), so no `OPS-18`-class drift and the ladder proceeded. **The stop rule was never approached** — the worst rung is 199 920 cells (67% of the 300 k ceiling) and 36.27 s to mesh (30% of the 120 s ceiling); no rung failed to mesh, and the log carries no "overlapping facets" and no `MeshAdapt` fallback, so there is no `GEO-23`-class record either. The table (`:7152–7155`): | `h_p` (m) | cells (`size_global`) | phantom cells (tag 3) | phantom share | mesh s | `2a_1g/h` | `2a_10g/h` | — | **0.0075** | 120 499 | 2 746 | 0.0228 | 25.26 | **1.65** | 3.57 | — | **0.005** | 129 505 | 8 497 | 0.0656 | 25.94 | **2.48** | 5.36 | — | **0.00375** | 144 212 | 18 570 | 0.1288 | 28.81 | **3.31** | 7.15 | — | **0.0025** | 199 920 | 58 866 | 0.2944 | 36.27 | **4.96** | 10.72 |. **Answer:** the first rung with ≥ 4 cells across the 1 g ball is **`h_p` = 0.0025 m** at **4.96** across (10.72 across the 10 g ball), costing **199 920 cells** (1.66× the record rung), 58 866 phantom cells (21.4× it), a 29.4% phantom share, and 36.27 s to mesh (1.44× it) — the three coarser rungs reach only 1.65 / 2.48 / 3.31, the 0.0075 figure confirming `MAT-4` step 4's "≈1.7 cells" reading. **Two ungated observations for the review that scopes `MAT-4` step 5:** (i) the cost is cheap in *total* cells because the refinement stays confined to the phantom — total cells grow 1.66× and mesh time 1.44× while phantom cells grow 21.4× — so a 1 g-resolved F-small mesh is a ≈200 k-cell / ≈36 s build, not the runaway the phantom's `1/h³` alone suggests; (ii) the phantom count tracks `(h_ref/h)³` from slightly below throughout (3.09× vs 3.38, 6.76× vs 8, 21.4× vs 27) and cell growth is monotone in `1/h` on every rung (no `GEO-22`-class non-monotonicity here). **Caveat the consuming step must carry:** the three finer rungs are **single readings** — only the 0.0075 rung has a cross-process repeat (today against 2026-09-02) — so none of them may be pinned as a record without its own repeat. Nothing was asserted, no band moved, no `src/` change, `PHANTOM_GROWTH_BAND` not re-asserted; the row stays 🧪 until a gate consumes the table. | 🧪 *(opened 2026-09-07; **table delivered 2026-09-08**; measurement-only, no band, no `src/`)* 

### §9 item 377 — 2026-10-11 (`GEO-27` step 2 — the 1 g rung's price at degree 2 on exact rectangles; lane-1 revalidation)

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.*

**Weekly re-pose being executed.** The 2026-10-10 weekly (`PROJECT_PLAN.md:4810–4814`, listed under "Lane 2 (revalidation / Tier B)") reads: "`GEO-27` 🧪 re-posed as step 2 — the 1 g rung's price on the exact-rectangle degree-2 F-small fixture (`h_p` = 0.005 / 0.0025 m, `-n 8`): cells, unknowns, time, memory in a table, `MAT-4` step 5's 1 g C4 identity asserted at the imported band on the rung that fits ⇒ ✅ (a measurement with one imported gate)".
- It is 🧪 → ✅ work, so it falls in the revalidation class (§9 rule (5)'s lane 1 by kind).
- **Target: the weekly's text.** No band, tolerance or validation target is set or moved. `C4_COVARIANCE_BAND` is `MAT-4`'s, imported.
- The row's Done-when must be amended to match (rule (j); exact text below).

**Executor and cost.**
- **Executor:** implementer. It is not `mesh-probe`: the item solves and asserts.
- **Tier:** heavy, `-n 8`, complex build. **Class:** lane-1 revalidation (median 17).
- **Cost:** ≈ M measured / ≈ 50 predicted min: 15 fixed + 15 code + W0 ≈ 1 + W1 ≈ 4 + W2 ≈ 7 + spare.
- **Depends on:** nothing.
- **Memory guard: coarse first, with the stop rule below.**

**The fixture is already priced at 0.0025 m (read this first).** The exact-rectangle degree-2 F-small at `h_p` = 0.0025 m is `IMP-3` I14's *off* arm.
- **Construction:** `_build(True, phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` → `build_four_port_sweep(built=…, build_only=True, degree=2)` → `run_n_port_sparameter_sweep(…, keep_fields=True, degree=2)` (`tests/validation/test_imp3_tip_ratio_1g.py:121–178`).
- **Measured** at `-n 8` on 2026-10-06 (`20261006T111406Z_IMP-3-item263-W2-d2off.log:2029, :2080–2084`): 201 072 cells, 1 284 856 degree-2 dofs = 2E + 2F; build 43.59 s; one-factorisation 4-port sweep 246.15 s; 12 ball SARs 26.26 s; window 397 s; **summed peak RSS 44.17 GiB**.
- The integers are registered as `I14_OFF_CELLS, I14_OFF_DOFS` (`tests/validation/test_imp3_tip_ladder_1g.py:94`).
- **So 0.0025 m is predicted to fit the heavy tier at `-n 8`.** 44 GiB is well under the container's 128 GiB limit (`docker/docker-compose.yml:34`), and the window is under the 590 s container ceiling even on the 1.32× box (≈ 525 s).
- That run was at 64 MHz. This item runs at `MAT-4`'s 10 MHz (`build_four_port_sweep`'s default `FREQUENCY_HZ` = 10.0e6, `tests/validation/test_port_reaction_impedance.py:108`). It is the same matrix pattern, so the price carries over.
- The 0.005 m rung has **no** degree-2 measurement on this fixture; it is predicted below.

**Why not `_build_mass_averaged`.** `MAT-4` steps 4–5 call `_build_mass_averaged(phantom_resolution=…)` (`tests/validation/test_birdcage_sar_mass_averaged.py:209–331`). That routine:
- builds the ragged default fixture (no `built=`, no strip fraction);
- runs at degree 1 (`build_four_port_sweep`'s default);
- solves with **four fresh factorisations** (`_solve_driven` per port, `:249`). At degree 2 on 1.28 M dofs that is ≈ 4 × 200 s, which does not fit.

Editing it (an additive `built=` / `degree=` / held-factor path) would owe re-runs of `MAT-4` steps 4 and 5b under rule (a). Instead this item copies its centres-and-pairing block (`:239–331`) with a provenance comment (base commit, line span) and feeds it I14's one-factorisation sweep fields. `run_n_port_sparameter_sweep` holds one factorisation; with matched 50 Ω terminations its held drives equal fresh single-drive solves (the `PORT-19` premise). The pairing, centres, radii, density map, operator and `QUADRATURE_DEGREE` stay `MAT-4`'s.

**New module `tests/validation/test_geo27_one_gram_rung_degree2.py`.**
- `@complex_only`. Knob `FEM_EM_GEO27_STEP2` ∈ {`h005`, `h0025`}: unset ⇒ skip; any other value ⇒ `ValueError`. One test dispatches on it (never `-k`).
- **Imports (names only, never `test_*`):**
  - `_build`, `SHEET_IFACE` (`tests/mesh/test_birdcage_port_sheets.py`);
  - `GATED_WIDTH_FRACTION` (`tests/validation/test_port_lumped_narrowed_sheet.py:132`);
  - `build_four_port_sweep` (`test_port_birdcage_four_port.py:219`) and `run_n_port_sparameter_sweep`;
  - `PHANTOM_CELL_TAG` (`test_port_birdcage_lumped_column.py`);
  - `_sheet_areas` (`test_imp3_palace_twin.py:69`); `NOMINAL_STRIP_AREA` (`test_imp3_tip_ratio_1g.py:105`); `STRIP_BAND` (`test_ops61_exact_strip_ports.py:58`);
  - `I14_OFF_CELLS`, `I14_OFF_DOFS` (`test_imp3_tip_ladder_1g.py:94`); `PHANTOM_RESOLUTION_1G_RUNG` (`test_birdcage_sar_1g_rung.py:114`);
  - `C4_COVARIANCE_BAND` (`test_birdcage_b1_plus_map.py:120`);
  - `CENTRE_RADIUS_M`, `WHOLE_PHANTOM_BALL_RADIUS_M`, `EXACT_IDENTITY_RTOL`, `PHANTOM_RHO_KG_PER_M3`, `_global_tag_cell_count` (`test_birdcage_sar_mass_averaged.py:149, :158, :163`);
  - `ONE_GRAM_KG`, `TEN_GRAM_KG`, `QUADRATURE_DEGREE` (`test_mass_averaged_sar_standard_masses.py:121–124`);
  - `_port_index`, `QUADRATURE_STEP_DEG` (`test_birdcage_b1_quadrature.py`);
  - `DISCRETE_IDENTITY_RTOL`, `_exact_shares` (`test_birdcage_power_identity.py:227`);
  - `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py`);
  - `_matrix_quantities`, `_peak_rss_gib` (`test_imp3_rod_birdcage.py:275, :498`);
  - `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`).
- **Constants defined here:**
  - `PHANTOM_RESOLUTION_005 = 0.005`, step 1's rung (`20260908T033218Z_GEO-27.log:7152–7155`);
  - `MAT4_STEP5B_ONE_GRAM_PAIRS = {0: 0.000957, 1: 0.001199, 2: 0.001305, 3: 0.001065}` and `MAT4_STEP5B_ONE_GRAM_CONTROL_K0 = 0.870143` (the k = 0 mis-paired control), both from `20260908T200629Z_MAT-4-step5b.log.gz:1990–1993`. They are printed comparands (degree 1, ragged sheets, same `h_p`), never bands.
- **Body (both modes):**
  1. Build at the mode's `h_p`. Then `build_four_port_sweep(built=…, build_only=True, degree=2)` with the default frequency.
  2. **Before the sweep:** cells and (D) are checked, and the mode `h0025` checks the cells-and-dofs record.
  3. `run_n_port_sparameter_sweep(…, keep_fields=True, degree=2)`, timed.
  4. The copied `MAT-4` block:
     - azimuths → `_port_index` → `by_k`; centres at `CENTRE_RADIUS_M` on each drive's azimuth, z = 0;
     - ρ = `build_density_field` (0 off tag 3); σ = `result.fields[pid].sigma_field`;
     - balls: 1 g diagonal (k, k) ×4; 1 g mis-paired (k, k+2) ×4; 10 g diagonal ×4; the whole-phantom ball and `mean_sar` on tag 3 for P1.
  5. Per-drive exact power shares (`_exact_shares`, I14's `:181–184`).
  6. Print the price line: `h_p`, cells, tag-3 cells, E, F, d2 dofs, build s, sweep s, balls s, test-body s, `_peak_rss_gib()` and ranks.
  7. The `h005` window prints a **stop-rule line** (below).

**Anchors (asserted; every band imported).**
- **Both windows:**
  - (D) d2 dofs = 2E + 2F exactly.
  - (G2) the four N_i strip areas equal `NOMINAL_STRIP_AREA` within `STRIP_BAND`. This is the exact-rectangle fixture tie; I14 read 5.600000000000e-05 ×4 (`…d2off.log:2029`).
  - (i) ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND`.
  - (P) P_src = P_vol + P_sheet per drive ≤ `DISCRETE_IDENTITY_RTOL`. I14 read ≤ 8.4e-15 at 64 MHz on the 0.0025 mesh.
  - (ii) the whole-phantom ball: power = the tag-3 `mean_sar` power and mass = ρ·V(tag 3), both ≤ `EXACT_IDENTITY_RTOL`. This is `MAT-4` step 5's coverage anchor; it read 5.26e-14 at degree 1 (`…step5b.log.gz:1999`).
- **W2 (`h0025`) only — the gate:**
  - (R) cells within `CELL_COUNT_BAND` of `I14_OFF_CELLS`, and d2 dofs == `I14_OFF_DOFS`. Both integers are printed; mesh determinism on this fixture has been shown five times (GEO-36 W2, IMP-3 W1–W4, item 327 W1b).
  - **(C4-1g) `MAT-4` step 5's four cyclic 1 g C4 pairs**, \|SAR₁g(c_{k+1}; P_{k+1}) − SAR₁g(c_k; P_k)\| / SAR₁g(c_k; P_k), each ≤ `C4_COVARIANCE_BAND`.
  - The 10 g pairs are printed beside them.
- **W1 (`h005`):** the 1 g and 10 g pairs are **printed, not asserted.** There is no prior measurement at this rung and order, and 2.48 cells across the 1 g ball is the regime where degree-1 ragged read 2.9–6.8 % at 1.65 cells (`MAT-4` step 4).

**Negative controls (rule (e)).**
- **(N1g) the mis-paired 1 g control:** the far-side ball under the same drive, \|SAR₁g(c_{k+2}; P_k) − SAR₁g(c_k; P_k)\| / SAR₁g(c_k; P_k).
  - **The miss is asserted > `C4_COVARIANCE_BAND`** in W2: a C4 "identity" met by a field that is constant across the phantom would fail it. The ceiling is 100 % (a ratio of positive SARs).
  - **Its size is predicted, not asserted:** ≈ 87 %. Degree-1 ragged read 87.01–87.06 % at this `h_p` and frequency (`…step5b.log.gz:1990–1993`), on a different fixture (ragged sheets, degree 1). So the 5 % floor sits ≈ 17× under the prediction.
- **(N-sep, predicted, printed)** worst 1 g C4 pair ÷ smallest N1g control: the separation the gate shows, predicted ≥ 100×.

**Predictions** (written first in the module, printed).
- **0.005 m.** Cells in [128 k, 136 k]: ragged 129 505 (`GEO-27.log:7153`), plus the strip delta seen at 0.0025 (201 072 − 199 920 = 1 152).
- Tag-3 cells ≈ 8 500 (ragged 8 497).
- d2 dofs in [780 k, 900 k] (I14's 6.39 dofs per cell).
- Sweep ≈ 246.15 s × (N/1 284 856)^{1.23…1.37}, about 130–150 s. Summed RSS ≈ 44.17 GiB × the same ratio, about 22–30 GiB. Window 190–320 s.
- 1 g C4 pairs at 0.005 in [0.05 %, 3 %].
- **0.0025 m.** Cells 201 072 and d2 dofs 1 284 856 (equal).
- Sweep 230–330 s; summed RSS 42–50 GiB; window 330–520 s.
- 1 g C4 pairs in [0.01 %, 0.5 %] (degree-1 ragged read 0.0957–0.1305 %, exact rectangles at degree 2 should not be worse); N1g ≈ 87 %; (ii) ≤ 1e-13.

**Stop rule (coarse first; printed by W1, applied by the executor before W2).**
- Write N₁, T₁ and R₁ for W1's measured d2 dofs, test-body seconds and summed peak RSS.
- W1 prints T̂₂ = T₁·(1 284 856/N₁)^1.37 and R̂₂ = R₁·(1 284 856/N₁)^1.37, beside I14's measured 397 s / 44.17 GiB.
- **Do not start W2 if R̂₂ > 100 GiB** (≈ 78 % of the 128 GiB container limit, headroom for MUMPS peaks) **or if T̂₂ > 540 s**, which leaves room for the 590 s container ceiling.
- In that case record the row "0.0025 m: does not fit heavy at `-n 8` (predicted T̂₂ s, R̂₂ GiB)" in the table and commit W1 as a record. The item **holds 🧪** and reports.
- Whether W1's printed 0.005 pairs may then be the gate is the review's ruling (a rule (j) row matter), never in-slot.
- **Predicted: the rule does not fire** (I14 measured this rung at 44 GiB and 397 s).

**Windows** (foreground; the harness call is the whole Bash command, bare; Bash `timeout: 660000`; `-s`; `tests/environment/test_complex_mode.py` first; §5.1 durable capture copied verbatim with the trailing `; exit $rc`).
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W0** (smoke, `-n 1`): the module collects, the knob parser raises on a bad value, and the gated test skips when unset.

  `scripts/testing/run_and_log.sh GEO-27-step2-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo27_one_gram_rung_degree2.py -v -s --tb=short'"`
- **W1** (`h005`, heavy, `-n 8`, predicted 190–320 s):

  `scripts/testing/run_and_log.sh GEO-27-step2-W1-h005 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_GEO27_STEP2=h005 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_geo27_one_gram_rung_degree2.py -s -v --tb=short > /workspace/logs/geo27-step2-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/geo27-step2-w1-raw.log; cat /workspace/logs/geo27-step2-w1-raw.log; exit \$rc'"`
- **Then apply the stop rule; then run pgrep again (must read 0).**
- **W2** (`h0025`, heavy, `-n 8`, predicted 330–520 s; the window this step closes on): the same command with label `GEO-27-step2-W2-h0025`, `FEM_EM_GEO27_STEP2=h0025` and raw file `/workspace/logs/geo27-step2-w2-raw.log` (all three occurrences).
- **No module edit after W0 starts** (rule (i)). No `src/` change and no edit to any imported module.

**Traps already paid for.**
- **Four factorisations at degree 2 do not fit.** Use the one-factorisation sweep; never `_solve_driven` per port. I14's sweep is the measured price.
- **`build_four_port_sweep(built=…)` ignores the sizing keywords** (`:241`). `h_p` and the strip fraction go to `_build`, never to the sweep builder. (G2) is the guard that the exact strips arrived.
- **Orphan `.c` of any size** (item 160). After a kill or status 124, run `pgrep -c python3` first, then the rubric's cache sweep (label `GEO-27-step2-cache-sweep`), then park. Never raise the timeout or `-n`; a killed W2 is the price record ("does not fit heavy").
- **Rank-local values.** Use the reduced `_global_tag_cell_count`, `mass_averaged_sar` and `mean_sar`; never `cell_tags.values` lengths. Keep collectives out of `rank == 0` branches.
- **Complex build:** ball centres as plain floats; `np.real` before every `float`.
- **`-s` always; no `-k` inside the quoted command; no pipe.** Copy the durable-capture idiom verbatim; a dropped `; exit $rc` reads 0 over a red (`WF-6` 4h).
- **Two heavy windows never run at once.**

**Scope.** Covered: F-small, exact-rectangle sheets (`GATED_WIDTH_FRACTION`, κ off), the solid σ = 800 S/m coil with 50 Ω terminations, 10 MHz, degree 2, phantom rungs 0.005 and 0.0025 m, `-n 8`. **Not claimed:**
- 1 g convergence: two rungs are not a rate;
- absolute SAR, C95.3 compliance or Larmor;
- the copper / hole fixture or F-human;
- a 1 g gate at 0.005 m (printed);
- any `MAT-4` row change (`MAT-4` stays ✅ on its own degree-1 gate).

**Negative result.** Report and stop. Park on `attempt/GEO-27-step2-*` and mark item 377 BLOCKED in the same commit (rule (d)).
- **(D), (G2), (R), (i), (P) or (ii) red:** the fixture or the sweep is not I14's. File a known-issues entry, never re-band.
- **(C4-1g) red:** the 1 g C4 pairs at 0.0025 m, degree 2, exact rectangles miss the band that the degree-1 ragged fixture met by 38×. That is a finding about the degree-2 fixture's ball integration (`QUADRATURE_DEGREE` on P2 fields at 4.96 cells across). Record the per-pair table and the 10 g pairs, file a known-issues entry and hold 🧪. **No band moves.**
- **(N1g) not missing:** the far-side ball is not far, so the centres or the drive map are wrong. Stop.
- **The stop rule firing:** this is not a red. Record the row and hold 🧪 as above.
- **Status 124 / > 590 s on W2:** the "does not fit heavy" record; hold 🧪.
- No `PREDICTED RED:` line.

**Status it can move.** §7 `GEO-27` 🧪 → ✅ on W2 green, but only with the Done-when amended in the queueing commit (rule (j)). The row today says "stays 🧪 until a gate consumes the table" (`PROJECT_PLAN.md:1864`). **Recommended exact replacement for the row's state cell** (≤ 2 KB):

> 🧪 Step 1 (2026-09-08, `mesh-probe`): the phantom-resolution cost table on the ragged F-small — `h_p` 0.0075 / 0.005 / 0.00375 / 0.0025 m → 120 499 / 129 505 / 144 212 / 199 920 cells, 4.96 cells across the 1 g ball at 0.0025; only the 0.0075 rung repeated (`20260908T033218Z_GEO-27.log`, 132 s). **Done-when (re-posed by the 2026-10-10 weekly; amended 2026-10-11 under rule (j)):** step 2 — the 1 g rung priced on the exact-rectangle degree-2 F-small fixture at `h_p` = 0.005 and 0.0025 m, `-n 8` (cells, degree-2 unknowns, wall time, summed peak memory in one table; a rung that does not fit the heavy tier is recorded as such), **and** `MAT-4` step 5's four cyclic 1 g C4 pairs asserted at the imported `C4_COVARIANCE_BAND` on the finest rung that fits, in one logged window ⇒ ✅. A measurement with one imported gate; no band moves, no convergence claim. *History: `docs/planning/chunks/GEO-27.md`.*

On green, the executor appends the two-row price table and the W2 pairs to this chunk file. The row gets one clause ("✅ 2026-10-1x: priced at 0.005 / 0.0025 m, degree 2, `-n 8`; 1 g C4 pairs ≤ x % at 0.0025 m (`<W2 log>`)"). Nothing else moves; §2 and `MAT-4` are untouched.

**Execution order.** pgrep → W0 → pgrep → W1 → stop rule → pgrep → W2 → commit (module, harness logs, the price table, the row amendment if the review has not already applied it) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 377 — `GEO-27` step 2 (lane revalidation) — per `docs/planning/chunks/GEO-27.md` §9 item 377. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run W0, W1 (h005), then apply the stop rule, then W2 (h0025), foreground, with `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reading 0 before each heavy window. Never use per-port fresh factorisations; use the one-factorisation sweep. Never delete a log. No edit to the module after W0 starts; no edit to any imported module or `src/`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. If the stop rule fires, record the row and hold 🧪 (not a red). On a red, park on `attempt/GEO-27-step2-*`, mark item 377 BLOCKED in the same commit, and report."

**Outcome (2026-10-12, 22:30 CDT slot, implementer) — GREEN, `GEO-27` ✅.** Module `tests/validation/test_geo27_one_gram_rung_degree2.py`; W0 `20261012T034721Z_GEO-27-step2-W0.log` (3 s), W1 `20261012T034731Z_GEO-27-step2-W1-h005.log` (136 s), W2 `20261012T034956Z_GEO-27-step2-W2-h0025.log` (323 s), all `-n 8` complex except W0. Stop rule (W1): T̂₂ = 239.0 s, R̂₂ = 37.27 GiB — did not fire.

| `h_p` (m) | cells | tag-3 | E | F | d2 dofs | build s | sweep s | balls s | body s | summed peak RSS (GiB) | ranks |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.005 | 130 427 | 8 498 | 155 934 | 262 511 | 836 890 | 30.49 | 77.02 | 22.17 | 132.86 | 20.72 | 8 |
| 0.0025 | 201 072 | 58 960 | 238 627 | 403 801 | 1 284 856 | 41.67 | 246.74 | 27.04 | 319.07 | 44.05 | 8 |

W2 pairs (`:1990–1995`), k→k+1: 1 g C4 0.0237 / 0.0446 / 0.0349 / 0.0140 % (**asserted** ≤ `C4_COVARIANCE_BAND` 5 %); 10 g 0.0114 / 0.0264 / 0.0361 / 0.0211 % (printed); N1g 84.7411 / 84.7312 / 84.7243 / 84.7287 % (asserted > 5 %; predicted ≈ 87 %, read 2.3 points under); N-sep 1900.9× (predicted ≥ 100×). Anchors: (R) 201 072 / 1 284 856 equal I14; (G2) 5.600000000000e-05 ×4; (i) 1.265e-14; (P) ≤ 1.146e-14; (ii) power 1.432e-14, mass 3.331e-14. W1 (printed only): 1 g C4 0.0939 / 0.0323 / 0.1327 / 0.0066 %, 10 g 0.1269 / 0.0231 / 0.0233 / 0.1266 %, N1g 84.77–84.81 %. **Prediction misses (recorded, none asserted):** W1 sweep 77.02 s vs predicted 130–150 s, RSS 20.72 vs 22–30 GiB, window 136 vs 190–320 s (the 1.23–1.37 dof-power extrapolation from I14 over-predicted the coarse rung); W2 window 323 s under the 330–520 s band; W2 tag-3 58 960 (no prediction). No band moved, no `src/` or imported module edited; no convergence, absolute-SAR or `MAT-4` claim.
