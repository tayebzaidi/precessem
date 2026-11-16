# OPS-61 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 148 — 2026-10-01 (chain F10a; queued by the 10:30 review)

`OPS-61` step 1b — the exact f = 0.5 port strip on the `ANS-4` ×1 rung (implementer; standard, `-n 1` / `-n 2`, complex container; ≈ 6 measured / ≈ 54 predicted min; depends on nothing).

- **What "the exact-rectangle port set" means in code.** Today `_narrowed_radial` (`test_port_birdcage_leg_offset_sweep.py:166`) keeps the facets of the full 0.014 × 0.008 sheet F_i whose midpoints lie within f·w/2 of the sheet centre (`GATED_WIDTH_FRACTION` 0.5). The mesh does not conform to the strip edges, so N_i is ragged: 26 triangles, **5.835298880e-05 m² (fraction 0.521009)** on every port (`20260929T110345Z_OPS-61-step1a-W1.log`, the P1–P4 lines after `:1784`). Palace refuses that set (`Specified direction does not align … bounding box axes`, `…110504Z_OPS-61-step1a-W3.log:86`). The fix makes the strip edges CAD edges; the same midpoint filter then selects exactly the centre rectangle.
- **Code (`src/` keyword, disclosed, rule (c)).** `leg_sheet_strip_fraction: Optional[float] = None` on `MeshGenerator.birdcage_port_domain`. When set, in the leg-sheet block (`io/mesh.py:4112–4131`), before placement: `addRectangle(port_radius − f·port_dx/2, −port_dz/2, 0, f·port_dx, port_dz)` and `occ.fragment([(2, sheet)], [(2, centre)])` at azimuth 0 (so the cut is exact); place every resulting piece with `_place_sheet_in_leg_frame` and append all to `sheet_tags`. `sheet_of_ordinal` unchanged (same plane); the port box still splits in two halves, so F_i is unchanged. Raise `ValueError` unless `emit_port_sheets=True`, `leg_gap_length` set and 0 < f < 1; refuse `c4_congruent_sheets` (`:4876` needs one surface per sheet), `as_hole` and `symmetry_quadrant`. With `None`, no new line executes. The same keyword on `tests/mesh/test_birdcage_leg_offset._build` (`:107`), forwarded only when not `None` (the `c4_congruent_sheets` precedent).
- **Test.** New `tests/validation/test_ops61_exact_strip_ports.py`: builds through `_build(np.zeros(4), …)`, then `_sheet_areas`, `_sheet_azimuth_deg`, `_port_frame`, `_projected_extents`, `_narrowed_radial`, all imported, as the probe does at `:133–168`.
- **Anchors (asserted, bands imported).** (A1) flag on, every F_i = dx·dz = 1.12e-4 m² within `SHEET_AREA_BAND` (`test_birdcage_leg_offset.py:104`). (A2) every N_i = 0.5·dx·dz = 5.6e-5 m² within 1e-9; width re-measured as A/h = 0.007 within 1e-9. (A3) N_i ⊂ F_i and F_i∖N_i = 5.6e-5 within 1e-9. (A4) flag off, cells = `STEP2_CELL_COUNT` 116 085 within `CELL_COUNT_BAND` (the default path re-gated). (A5) `GEO-18`'s volume partition within `EXACT`, flag on. (A6) flag-on cell count printed in W1, registered after W1 as `ANS4_EXACT_STRIP_CELL_RECORD` (version-tagged, §10's licence for a new tagged count), asserted at `CELL_COUNT_BAND` from W2 on.
- **Negative control (asserted; backed on the same fixture and quantity).** Flag-off N_i misses 5.6e-5 by +4.20e-2 relative (step 1a W1's reading) — asserted > 1e-3, printed.
- **Predictions (printed).** Flag-on cells 116 085 × (1.001–1.01); F_i gains ≈ 10–30 % triangles (58 now); N_i facet count equal on all four ports.
- **Windows** (complex container, as the probe used; ×1 mesh 26.4 s, `…110345Z…W1.log:1784`).
  - W1 (≈ 90 s): `scripts/testing/run_and_log.sh OPS-61-F10a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 -m pytest tests/validation/test_ops61_exact_strip_ports.py -s -v --tb=short'"`
  - Register the record, then W2 (≈ 90 s): the same command labelled `OPS-61-F10a-W2` at `-n 2` — the closing window (rule (i)).
  - W0 (rule (c) re-gate, ≤ 600 s): the same shape running `tests/mesh/test_birdcage_leg_offset.py`, labelled `OPS-61-F10a-W0-regate`.
  - W3 (≈ 40 s): an `export --exact-strip` mode in `scripts/probes/ops61_step1a_palace_ans4.py` writing `logs/ops61/ops61_ans4_x1_strip.msh` and `f05x.json` (`Model.Mesh` that file; port i = {20+i}, `+Z`), asserting the new record instead of 116 085 plus the exterior area 0.3072; command as the probe docstring's W1, labelled `OPS-61-F10a-W3-export`.
- **Traps.** Never overwrite step 1a's `ops61_ans4_x1.msh`, `f1.json`, `f05.json` (F10b reads them). Fragment before rotating — fragmenting rotated rectangles moves edges by a rounding error that can break 1e-9. New-keyword builds at `-n 1` first (a gmsh failure deadlocks wider runs). Areas are rank-local. No edit after the last window. No Palace run. Harness bare, foreground, `timeout: 660000`.
- **Scope.** The port set, its identities and the export. Not: Palace (F10b), any S-parameter, κ, the full F_i drive.
- **Status it can move.** The `OPS-61` row records step 1b (Done-when written by the 10:30 review); no glyph. F10b becomes queueable; Palace's width objection to f = 0.5 removed by construction.
- **Negative result.** (A2) red with (A1) green: the strip edges did not reach the mesh — report N_i's span and facet count, park, BLOCKED. (A4) red: the default path moved — revert, known-issues entry. The build fails to mesh: a `GEO-21`-class finding, park.

**Outcome (2026-10-01 16:30 CDT slot, third item): DONE.** `leg_sheet_strip_fraction` on `MeshGenerator.birdcage_port_domain` (fragment at azimuth 0, then place every piece; `None` executes the original lines, re-gated), forwarded by `tests/mesh/test_birdcage_leg_offset._build`; new `tests/validation/test_ops61_exact_strip_ports.py`; `export --exact-strip` in `scripts/probes/ops61_step1a_palace_ans4.py`. W1 `20261001T215942Z_OPS-61-F10a-W1.log` (`-n 1`, 57 s, 7 passed, A6 skipped pre-record, `:3602`); W2 `20261001T220051Z_OPS-61-F10a-W2.log` (`-n 2`, 57 s, 8 passed, `:3617,3620`); W0 `20261001T220154Z_OPS-61-F10a-W0-regate.log` (`-n 2`, 82 s, 6 passed, `:5333`); W3 `20261001T220353Z_OPS-61-F10a-W3-export.log` (`-n 1`, 31 s, ANCHORS PASS `:1941`). (A1) F_i 1.120000000e-04, ratio 1.000000000000; (A2) N_i 5.600000000e-05 (+0.9 … +1.6e-15 relative), A/h 7.000000000e-03, radial span 7.000e-03; (A3) N ⊂ F, F∖N 5.6e-5 (≤ 2.2e-15); (A4) flag off 116 085 (`W1:3581`); (A5) box halves sum to 1.000000000000; (A6) flag on 117 206 cells (ratio 1.009657, `W1:3575`) registered as `ANS4_EXACT_STRIP_CELL_RECORD`, asserted in W2 (`:3608`). Control: flag-off N_i +4.202e-2 on every port (26 triangles, ragged span 9.167e-3; `W1:3595–3598`). W3: 117 206 tets, exterior 0.3072 to 5.6e-17, `logs/ops61/ops61_ans4_x1_strip.msh` + `f05x.json` written (gitignored); step 1a's files untouched. Prediction miss (printed): F_i triangles +52–55 % (88–90 vs 58, P2 88), not 10–30 %; N_i 40 on all four ports as predicted. Deviations: `EXACT` imported from `test_birdcage_conductor_hole.py:117`, `CELL_COUNT_BAND` from `test_birdcage_port_sheet_prerequisite.py:57`; one extra refusal test (f = 1.0 raises). No edit after each file's last window (mtimes checked).

### §9 item 159 — 2026-10-01 (chain F10b, first half; queued by the 18:00 review)

`OPS-61` step 1c, first half (F10b-i): Palace `f05x` at 64 MHz, its identities asserted at our bands, and ours on the same strip mesh at degree 1, κ off and in-run (implementer; heavy, Palace `-np 8` / ours `-n 2`, complex container; ≈ 15 measured / ≈ 45 predicted min; depends on nothing, since F10a landed as item 148). **Split by the 18:00 review on price, not re-targeted.** Our degree-2 half is F10b-ii, queued as item 166 and serial on this item. Together the two halves make the row's step 1c.

- **Framing.** Step 1a's cross-solver class gaps (Palace `f1` against our ragged f = 0.5: 1.315 / 0.709 / 0.593 %, `20260929T110544Z_OPS-61-step1a-W5.log:48–50`) mix up port width with solver. Palace refused the ragged set (`…110504Z_OPS-61-step1a-W3.log:86`). F10a made N_i an exact CAD rectangle: 40 triangles and 5.600000000e-05 m² on every port, 117 206 tets, written to `logs/ops61/ops61_ans4_x1_strip.msh` and `f05x.json` (`20261001T220353Z_OPS-61-F10a-W3-export.log:1817,1819–1822,1826,1828`). Both solvers now see the same mesh, the same port set and the same 50 Ω, so the remaining gap is between the solvers.
- **Code (tests and probe only; no `src/`; additive; disclosed).**
  - Probe `scripts/probes/ops61_step1a_palace_ans4.py`:
    - New `flip-config`. It loads `logs/ops61/f05x.json` and asserts that `Model.Mesh` is `STRIP_MESH_NAME` and that there are four ports with Attributes `[21]…[24]`, all `+Z`. It writes `flipx.json` only: port 1 `"-Z"`, `Problem.Output` `"postpro_flipx"`.
    - New `read --exact-strip`, dispatched like `export --exact-strip` at `:524`. It loops over `("f05x", "flipx")` only, reusing `_read_port_s` (`:381`) and `_identities` (`:409`). Today `read()` loops over `f1` / `f05` / `flip` only (`:444`) and asserts on `f1` only (`:504–520`).
    - Default `read` stays byte-identical.
  - `_four_port_rung` (`tests/validation/test_port_birdcage_leg_offset_sweep.py:199`; today it cannot build the strip mesh, `:274–279`): new keyword `leg_sheet_strip_fraction=None`.
    - It goes to `_build` only when not `None`, following the `c4_congruent_sheets` precedent. `_build` already accepts it (`tests/mesh/test_birdcage_leg_offset.py:112,143–144`).
    - It is ignored under `reuse`. Docstring paragraph as at `:221–223`. Re-gated by W5.
  - New `tests/validation/test_ops61_palace_strip_comparison.py`:
    - `@complex_only`, controlled by `FEM_EM_OPS61_DEGREE ∈ {1, 2}`. When it is unset every test skips (the `PORT-22` precedent).
    - Rung κ-off: `_four_port_rung("strip 64 MHz", np.zeros(4), frequency_hz=FREQUENCY_64_HZ, degree=d, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`. `FREQUENCY_64_HZ` is from `test_lossy_sphere_fullwave.py:94`.
    - Rung κ-on: the same call with `reuse=_reuse_of(off)` (`test_port_birdcage_larmor_gate.py:184`) and `width_correction_kappa="in_run"`.
    - Palace's S is read through the probe, loaded by `importlib` (precedent `tests/validation/test_wf5_q_two_ways.py:71–78`; `scripts/` has no `__init__`).
    - Rank 0 saves the S matrices, κ and the readings to `logs/ops61/ours_f05x_p{d}.json` (gitignored).
- **Anchors (asserted; bands imported, never restated).** `RECIPROCITY_BAND` is at `tests/validation/test_port_lumped_sheet_sweep.py:86`. `PASSIVITY_SIGMA_TOLERANCE` and `POOLED_SEPARATION_FLOOR` are at `tests/validation/test_port_birdcage_four_port.py:129,173`. `ADJACENT_SPREAD_BAND` is at `tests/validation/test_port_birdcage_lumped_column.py:121`. The gates are as `test_port_birdcage_larmor_gate.py:33–39,445–513` defines them.
  - **(P) Palace `f05x`, checked in W4 by `read --exact-strip`** (step 1a's `f1` assertions at `:503–520`, extended to `PORT-11`'s full (ii) and (iii′)):
    - ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND`.
    - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every column power sum ≤ 1 + 1e-9.
    - Every S-class spread ≤ `ADJACENT_SPREAD_BAND`.
    - Pooled off-diagonal spread ≥ `POOLED_SEPARATION_FLOOR` × the worst class spread.
    - The classes come from `_circulant_classes` / `_class_spread` (`test_port_birdcage_four_port.py:193,176`). Palace is checked on S classes, because `z_matrix_ohm` is the terminated transimpedance.
  - **(F) Same fixture, W6:**
    - Our cell count equals `ANS4_EXACT_STRIP_CELL_RECORD` 117 206 (`test_ops61_exact_strip_ports.py:64`) within `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`).
    - Every narrowed sheet area equals 5.6e-5 within `STRIP_BAND` 1e-9 (`:58`), and w = A/h = 7.0e-3 within the same band.
  - **(G) Our gates on the strip mesh, W6, κ-off and κ-on rungs:** (i) `RECIPROCITY_BAND`; (ii) σ_max and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`; (iii′) every Z-class spread ≤ `ADJACENT_SPREAD_BAND`, with pooled ≥ `POOLED_SEPARATION_FLOOR` × worst.
  - **(K) κ reaches the solve:** the κ-on vs κ-off self-class S move is > `KNOB_REACHES_THE_SOLVE` 1e-6 (`test_port21_feed_sensitivity.py:156`).
- **Printed, never gated (W4 and W6).**
  - Per class, \|Δ mean\|/\|mean\| for ours-strip p1 κ-off vs Palace `f05x`, and ours-strip p1 κ-on vs Palace `f05x`. Uses `_class_means` / `_relative_move`, `:287,304`.
  - Our S-class spreads beside Palace's.
  - Ours-strip vs ours-ragged `PORT11_S_64MHZ_RECORD` (`test_ans4_resolution_ladder.py:209`). This is the width confound on our side.
  - Palace `f05x` vs Palace `f1`. This is the width confound on Palace's side, also confounded by the mesh (116 085 vs 117 206 tets).
  - Palace `f05x` vs `ans4` `metrics.json` 64 MHz `s_matrix`, which is ours, ragged. The `aed` key is never read.
  - The in-run κ, pooled and per sheet.
- **Negative controls (rule (e)).**
  - **(C1) Palace flip on this fixture — predicted, printed, never asserted.** Step 1a's flip was measured on `f1`, a different port set and mesh, so it does not back an assertion here.
    - Predicted ‖S_flipx − D·S_f05x·D‖/‖S_f05x‖ ≤ 1e-6, with D = diag(−1, 1, 1, 1). Step 1a measured 4.279e-11 on `f1` (`…W5.log:69`).
    - Predicted flipx adjacent and opposite spreads ≈ 200 %, which is the ceiling (a pure sign flip of one row and column gives exactly 2.0). Step 1a measured 200.0159 / 200.0000 % (`:59`). The most the control can separate from the band is 2 / 0.005 = 400×.
  - **(C2) κ move — the floor is asserted (K); the size is predicted** ≈ 1.55e-2, from the ragged degree-1 move 1.552110e-02 (`20261001T171613Z_PORT-21-step5-W2.log:5654`).
- **Predictions (written first, printed).**
  - **Palace `f05x`:**
    - Log shows 117 206 elements and Rs = 4.375e+01 Ω/sq on attributes 21–24 (= 50·0.007/0.008; `f1` printed 8.750e+01, `…110423Z…W2.log:79–86`). This is a free geometry check: it confirms Palace sees the 7 mm strip.
    - About 25 s and about 1.8 G peak.
    - Class spreads ≲ 0.1 % (`f1`: 0.0712 / 0.0562 / 0.0370 %).
  - **Ours at degree 1:**
    - κ-off spreads ≈ 0.04–0.08 % and σ_max ≈ 0.9997 (ragged: 0.0573 / 0.0599 / 0.0370 %, 0.999721388, `…171613Z…W2.log:5601`).
    - The in-run κ is of the same order as the ragged value 1.0607621546e-02 (`:5592`).
  - **Our strip-vs-ragged self-class move:** 0.5–6 %, with the coupling classes moving less.
  - No prediction is made for the ours-vs-Palace size. That reading is the purpose of the item.
- **No `PREDICTED RED:` line.** Any red stops the slot.
- **Windows** (in order; compose service up; every Bash call `timeout: 660000`; each harness call is the whole command, bare and repo-relative).
  - **Pre-check (not a harness call):** `ls logs/ops61/`. If `ops61_ans4_x1_strip.msh` or `f05x.json` is missing, run W0 first. If both are present, do not run W0.
  - **W0 (conditional, ≈ 31 s)** — F10a W3 verbatim, relabelled:
    `scripts/testing/run_and_log.sh OPS-61-F10b-W0-export "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py export --exact-strip'"`
  - **W1 (≈ 5 s):**
    `scripts/testing/run_and_log.sh OPS-61-F10b-W1-flipcfg "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 60 python3 scripts/probes/ops61_step1a_palace_ans4.py flip-config'"`
  - **W2 (≈ 25 s)** — step 1a's Palace command with `f05x.json`:
    `scripts/testing/run_and_log.sh OPS-61-F10b-W2-palace-f05x "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x.json'"`
  - **W3 (≈ 25 s):** the same with `flipx.json`, labelled `OPS-61-F10b-W3-palace-flipx`.
  - **W4 (≈ 3 s):**
    `scripts/testing/run_and_log.sh OPS-61-F10b-W4-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip'"`
    This is the probe's last window; no probe edit after it.
  - **W5 (the re-gate for the helper edit, rules (a)/(c); ≈ 220 s, `20261001T171226Z_PORT-21-step5-W1.log:5675`)** — that log's `:12` verbatim, relabelled:
    `scripts/testing/run_and_log.sh OPS-61-F10b-W5-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_port_birdcage_larmor_gate.py -s -v --tb=short'"`
  - **W6 (the closing window, rule (i); ≈ 80–110 s)** — strip mesh ≈ 28 s; four drives 10.34 s and κ-on 11.02 s at `-n 2` (`…171613Z…W2.log:3664,5523`):
    `scripts/testing/run_and_log.sh OPS-61-F10b-W6-ours-p1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_OPS61_DEGREE=1 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ops61_palace_strip_comparison.py -s -v --tb=short'"`
- **Traps.**
  - Never overwrite step 1a's `ops61_ans4_x1.msh`, `f1` / `f05` / `flip` or their `postpro_*`. Never re-export when the strip files are present.
  - Palace aborts on any port set that is not aligned with the axes. An abort on `f05x` is a finding: never change `Direction` in-slot.
  - `-np` ≤ 12, and the compose command is quoted whole inside `run_and_log.sh`.
  - Under `reuse` the build keywords are not read, so κ-on must reuse the κ-off rung.
  - Compare S classes only; never read `z_matrix_ohm` or the `aed` key.
  - Areas and facet counts are per-rank: reduce before asserting.
  - `-s` on every pytest window, with `tests/environment` first and `FEM_EM_REQUIRE_COMPLEX=1`.
  - Run foreground; on a kill, `pgrep -c python3` first and never lengthen the timeout.
  - No edit to the probe after W4, or to the new module or the helper after W6.
  - Record row and journal: relative gaps and identity readings only. No absolute `S₁₁` / `Z_in`. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.
- **Scope.** Covered: 64 MHz, Palace order 1 (`f05x` as written), ours at degree 1. Not covered: our degree 2 (item 166), Palace order 2 (not in the row; the weekly's to add), 10 / 128 MHz, the PEC column, step 2, any band move, `src/`, and the (a)/(b)/(c) verdict.
- **Status it can move.** The `OPS-61` row records step 1c's first half: Palace-strip identities and ours-vs-Palace at degree 1, κ on and off. No glyph moves; 🟡 holds. Item 166 becomes runnable. The (a)/(b)/(c) reading needs the private AED column, so it belongs to the weekly.
- **Negative result.**
  - Palace aborts on `f05x`, or its log does not show 117 206 elements and Rs 43.75 Ω/sq: record verbatim, BLOCKED (rule (d)), file a known-issues entry, and do not run W5–W6.
  - (P) red: report the numbers per identity and park as BLOCKED. Never widen.
  - W5 red: the helper edit moved the default. Revert and stop before W6.
  - (F) red: the solved mesh is not Palace's mesh. Stop.
  - (G) red with the ragged `PORT-11` green: the strip changed the fixture. File a known-issues entry and park on `attempt/*`.
  - (K) red: fix once; if still red, park.
  - C1 not as predicted: print it and report; this does not stop the slot.

**Outcome (2026-10-01 21:00 CDT slot, second item): DONE — step 1c first half recorded, no glyph; item 166 runnable.** Pre-check: both strip files were present, so W0 was not run. Windows, all Status 0: W1 `20261002T021829Z_OPS-61-F10b-W1-flipcfg.log` (1 s), W2 `…021838Z…W2-palace-f05x.log` (26 s), W3 `…021913Z…W3-palace-flipx.log` (24 s), W4 `…021945Z…W4-read.log` (2 s), W5 `…022055Z…W5-regate.log` (148 s; 18 passed, 1 skipped — the same skip as the reference `20261001T171226Z_PORT-21-step5-W1.log:5517`), W6 `…022338Z…W6-ours-p1.log` (103 s; 16 passed). File times precede their last windows (probe 21:18:23 < W4 21:19:45; helper 21:17:59 < W5/W6; new module 21:20:50 < W6 21:23:38).
- **Palace log check:** 117 206 elements, Rs 4.375e+01 Ω/sq on attributes 21–24, 23.35 s, peak 1.8 G (as predicted).
- **(P) Palace `f05x`:** reciprocity 1.556e-10; σ_max 0.999760093668; column powers ≤ 1 + 1e-9; S-class spreads 0.0852 / 0.0783 / 0.0921 % (predicted ≲ 0.1 %); pooled 40.49 % ≫ 10 × worst.
- **(F):** 117 206 cells (ratio 1.000000); each sheet has 40 facets, area 5.6e-5 and w 7.0e-3, both to ≤ 1.6e-15.
- **(G):** κ off — reciprocity 4.5e-16, σ_max 0.999760093664, Z spreads 0.0747 / 0.0691 / 0.0897 %, separation 451.7×. In-run κ — reciprocity 6.8e-16, σ_max 0.999760150829, the same spreads.
- **(K):** self-class move 2.253922e-05 > 1e-6.
- **Ours vs Palace (printed):** κ off **0.0000 %** in all three classes (`W6.log:1975`; the executor read ≈ 2e-10 / 1e-10 / 6e-11 from the gitignored `ours_f05x_p1.json` — not a logged figure). σ_max agrees to 4e-12. In-run κ: 0.0023 / 0.0013 / 0.0010 %. The step 1a class gaps (1.3 / 0.7 / 0.6 %) were port width. On the same mesh at order 1, the two degree-1 solves agree to solver precision. **For the review:** agreement this exact suggests both codes assemble the same discrete system. That makes it a strong check of our port and impedance-sheet bookkeeping against an independently written code, but not of discretisation error; the order-mismatched comparison (item 166) is the one that probes that.
- **C1 (predicted):** 1.023e-10 ≤ 1e-6; flipx spreads 200.0204 / 200.0000 % — as predicted.
- **C2 (predicted, missed):** the κ move is 2.25e-5, not ≈ 1.55e-2. The in-run κ is pooled 1.5515e-05 (per sheet 1.59 / 1.30 / 1.74 / 1.38e-05), ≈ 680× below the ragged 1.06e-2: on an exact rectangle, the terminal-form deficit κ corrects nearly vanishes. Interpretation, not tested: κ on the ragged set was compensating the ragged width.
- **Width confounds (printed):** ours strip vs ragged `PORT11_S_64MHZ_RECORD` 1.539 / 0.858 / 0.713 % (inside the predicted 0.5–6 %, couplings moving less); Palace `f05x` vs `ans4` ragged-ours the same figures; Palace `f05x` vs Palace `f1` 2.085 / 1.112 / 0.919 % (width and mesh both change).

### §9 item 166 — 2026-10-01 (chain F10b, second half; queued by the 18:00 review)

`OPS-61` step 1c, second half (F10b-ii): ours on the strip mesh at degree 2, κ off and in-run (implementer; heavy, `-n 8`, complex; ≈ 15 measured / ≈ 25 predicted min; **depends on item 159 landing — if it did not, skip this item**). It is one window with no code change: item 159's W6 with `FEM_EM_OPS61_DEGREE=2`, `timeout -k 30 600`, `mpiexec -n 8`, labelled `OPS-61-F10b-ii-W1-ours-p2`.

- **Anchors.** Item 159's (F), (G) and (K) at degree 2 (the same module, the same imported bands).
- **Cost.** ≈ 280 s and ≈ 19.4 GiB peak summed, from `20261001T094205Z_PORT-21-step4-W1.log`: 746 144 dofs (`:2083`), mesh 23.55 s and four drives 66.48 s (`:2403`), κ route 143.39 s (`:2943`), 19.384 GiB (`:2941`), window 284 s (`:3576`).
- **Printed.**
  - Ours p2 vs Palace p1, with an order-mismatch note.
  - Ours p2 vs ours p1, read from item 159's JSON. Predicted ≈ 8.2 / 4.1 / 4.4 % from the C4-congruent mesh (`…094205Z…:3049`).
  - κ at degree 2, predicted ≈ 5.2e-2.
- **Trap.** This is the first `-n 8` build of the strip mesh. A gmsh failure deadlocks, so a Status 124 is the finding, not a retry.
- **Status it can move.** The `OPS-61` row's step 1c is complete as a record (no glyph).
- **Negative result.** Any (F) / (G) / (K) red at degree 2: park, BLOCKED; the review reads it.

**Outcome (2026-10-02 04:30 CDT slot): DONE — step 1c complete as a record, no glyph.** No code change. W1 `20261002T095914Z_OPS-61-F10b-ii-W1-ours-p2.log` (Status 0, 280 s, 16 passed at `-n 8`; the first `-n 8` strip build ran cleanly: mesh 23.47 s, four drives 72.33 s κ off and 70.61 s κ in-run, `:2031,2080`). A reading window, W2 `20261002T100407Z_OPS-61-F10b-ii-W2-p2-vs-p1-read.log` (Status 0, 2 s), computed ours p2 vs ours p1 from the two gitignored JSONs with the module's own `_class_means` / `_relative_move`.
- **(F):** 117 206 cells vs record (ratio 1.000000); each sheet has 40 facets, area 5.6e-5 and w 7.0e-3 to ≤ 8.9e-16 (`W1.log:2137–2142`).
- **(G):** κ off — reciprocity 2.569e-15, σ_max 0.999819887372, max column power 0.834996, Z spreads 0.0333 / 0.0338 / 0.0441 %, separation 957.45×. In-run κ — reciprocity 2.571e-15, σ_max 0.999819943311, the same spreads (`:2148,2151`).
- **(K):** self-class move 2.425986e-05 > 1e-6 (`:2154`).
- **Ours p2 vs Palace p1 (printed, order mismatch):** κ off 6.6674 / 2.5442 / 4.5224 %; in-run κ 6.6678 / 2.5442 / 4.5228 % (self / adjacent / opposite, `:2126–2127`). This is a degree-2 vs order-1 comparison, so it measures discretisation, not port-model disagreement.
- **Ours p2 vs ours p1 (printed, W2):** 6.6674 / 2.5442 / 4.5224 % κ off, 6.6675 / 2.5442 / 4.5225 % in-run. Predicted ≈ 8.2 / 4.1 / 4.4 % from the C4-congruent mesh: self and adjacent landed lower, opposite on. These figures equal the p2-vs-Palace ones to the printed digit, which item 159's 0.0000 % p1 agreement implies.
- **κ at degree 2 (predicted ≈ 5.2e-2, missed):** pooled 1.7919e-05 (per sheet 1.81 / 1.60 / 2.07 / 1.69e-05), the same order as degree 1's 1.55e-5. As at degree 1, the exact rectangle leaves almost no terminal-form deficit for κ to correct.
- **Ours strip p2 vs ours ragged `PORT11` record (printed):** 7.1265 / 2.6949 / 4.8011 %. This mixes degree and width, because the record is degree 1.
- **For the review:** the 6.7 % self-class shift from p1 to p2 on this mesh is the order-1 discretisation error, and Palace at order 1 shares it. Palace order 2 (not in this step) is the like-for-like check.

### §9 item 194 — 2026-10-03 (chain F, F17 — daily provisional, rule (6))

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.* Licence: rule (6), toward chain F's F10b target (Palace's own reciprocity / passivity asserted; ours vs Palace per class printed at both orders) and the operator directive 2026-10-03 (a) / (c): Palace carries the routine external checks. No band, tolerance, Done-when or validation target moves.

`OPS-61` step 1d: the `ANS-4` exact-strip Palace twin at 10 / 128 MHz (order 1) and at 64 MHz (order 2, a cost probe).
- **Executor:** implementer. Probe and test module only; no `src/`.
- **Tier:** heavy. Palace at `-np 8`; ours at `-n 2`. Complex build.
- **Class:** chain step (solve).
- **Time:** ≈ 14 measured / ≈ 35 predicted min.
- **Depends on:** nothing. It uses F10a's exported mesh, not item 192's records. It runs beside 192 and 193.

**Why.**
- Step 1c compared ours with Palace on the exact f = 0.5 strip mesh at 64 MHz only.
  - At matched order 1 they agree to 0.0000 % in every class (item 159, `20261002T022338Z_OPS-61-F10b-W6-ours-p1.log:1975`).
  - Ours at degree 2 sits 6.67 / 2.54 / 4.52 % from Palace order 1. That is an order mismatch (item 166).
- Two pieces of F10b's "both orders" reading are missing:
  - Palace at the fixture's other two frequencies;
  - Palace order 2, the like-for-like partner of our degree 2 (both use 20 unknowns per tet).
- Both run on the mesh already exported: `logs/ops61/ops61_ans4_x1_strip.msh`, 117 206 tets, port i = {20+i}, where N_i is the SPEC rectangle (F10a W3, `20261001T220353Z_OPS-61-F10a-W3-export.log`). No new export is needed.

**Code** (additive; disclosed):
- **Probe `scripts/probes/ops61_step1a_palace_ans4.py`.**
  - **New `freq-configs`.** It loads `f05x.json` and checks it as `flip_config` does (`:525–534`), and also checks `Solver.Order == 1` and `Freq == [0.064]`. It then writes three files, and only these three:
    - `f05x_f10.json`: `Freq [0.010]`, `Output "postpro_f05x_f10"`;
    - `f05x_f128.json`: `[0.128]`, `"postpro_f05x_f128"`;
    - `f05x_o2.json`: `Order 2`, `[0.064]`, `"postpro_f05x_o2"`.
  - **New `read --exact-strip-set`,** a new function with no edit to `read_exact_strip`. For each of the three runs present it applies item 159's (P) set at the same imported bands:
    - reciprocity ≤ `RECIPROCITY_BAND`;
    - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`;
    - every column power ≤ 1 + 1e-9;
    - every S-class spread ≤ `ADJACENT_SPREAD_BAND`;
    - pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst spread.

    It prints:
    - Palace o2 vs Palace o1 (`postpro_f05x`) per class;
    - Palace o2 vs our degree 2 from `logs/ops61/ours_f05x_p2.json` `s_kappa_off` (item 166, gitignored), or "absent";
    - the order-mismatch separation, ratio = \|ours p1 − Palace o2\| / \|ours p2 − Palace o2\| per class.

    Never read `z_matrix_ohm` or the `aed` key.
  - `export`, `read`, `read --exact-strip` and `flip-config` stay byte-identical.
- **`tests/validation/test_ops61_palace_strip_comparison.py`.** New knob `FEM_EM_OPS61_FREQUENCY_MHZ` ∈ {10, 64, 128}; any other value raises.
  - **Unset:** today's behaviour, byte-identical — 64 MHz, `postpro_f05x`, `ours_f05x_p{d}.json`.
  - **Set**, frequency from the constants:
    - 10 MHz: `FREQUENCY_HZ`. 64 MHz: `FREQUENCY_64_HZ`. 128 MHz: `FREQUENCY_128_HZ`.
    - Palace CSV by (MHz, degree): (10, 1) → `postpro_f05x_f10`; (64, 1) → `postpro_f05x`; (64, 2) → `postpro_f05x_o2`; (128, 1) → `postpro_f05x_f128`. Any other pair prints "no order-matched Palace run" and skips the Palace prints.
    - Output JSON `ours_f05x_p{d}_f{MHz}.json`.
    - Ragged comparand: `LEG_D_S_MATRIX_10MHZ` at 10 MHz, `PORT11_S_64MHZ_RECORD` at 64 MHz, "n/a" at 128 MHz.
    - (K) is printed, not asserted, at 10 and 128 MHz. Its size there is unmeasured (rule (e)); at 64 MHz it stays asserted as today.
    - (F) and (G) are unchanged.

**Anchors** (asserted; bands imported, never restated):
- **(P)** on Palace `f05x_f10`, `f05x_f128` and `f05x_o2`, in W5.
- **Palace log check, by the executor, recorded in the journal.** Each W2–W4 log shows 117 206 elements and Rs 4.375e+01 Ω/sq on attributes 21–24 (item 159's free geometry check). W4 also prints Palace's ND order 2 / unknown count.
- **(F)** in W6 / W7. Our cells equal `ANS4_EXACT_STRIP_CELL_RECORD` within `CELL_COUNT_BAND`. Every N_i is 5.6e-5 and w is 7.0e-3, within `STRIP_BAND`.
- **(G)** `PORT-11` (i), (ii) and (iii′) on our κ-off and κ-in-run rungs at 10 and 128 MHz.
- **W8 (default path):** item 159's W6 assertions as they stand.

**Negative control (rule (e)).**
- **(C2′) Order-mismatch separation — predicted, printed.** \|ours p1 − Palace o2\| per class ≈ 6.67 / 2.54 / 4.52 %, the size of our own p1-to-p2 move on this mesh (item 166 W2, `20261002T100407Z_OPS-61-F10b-ii-W2-p2-vs-p1-read.log`). The matched-order gap is predicted ≤ 1e-3 relative, so the separation should be ≥ 25×.
- **(P)'s power to discriminate** was shown on this mesh at 64 MHz: flipped port-1 spreads read 200.0204 % against the 0.5 % band, the sign-flip ceiling (item 159 C1). The flip is not re-run at 10 / 128 MHz.

**Predictions** (written first, printed):
- **Palace `f05x_f10` / `f05x_f128`:** 117 206 elements, Rs 43.75 Ω/sq; ≈ 25 s each and ≈ 1.8 G peak (`f05x` 23.35 s, 1.8 G); S-class spreads ≲ 0.15 %.
- **Ours p1 vs Palace o1 at 10 and 128 MHz:** ≤ 1e-6 relative in every class, i.e. "0.0000 %" printed; σ_max agreeing to ≤ 1e-10. Both codes assemble the same discrete system (item 159's reading).
- **Palace o2 (cost probe — first Palace order-2 run here):** 60–400 s, peak 10–40 G at `-np 8`. Unknowns ≈ 5.4× order 1; our degree-2 count on the sibling mesh is 746 144 (`20261001T094205Z_PORT-21-step4-W1.log:2083`). The service limit is 128 G.
- **Palace o2 vs ours p2:** ≤ 1e-3 relative per class, because the two spaces are the same 20-dof-per-tet first-kind Nédélec. No tighter prediction is made: the size of this gap is the purpose.
- **Palace o2 vs Palace o1:** ≈ 6.7 / 2.5 / 4.5 %.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12).
- **Pre-check** (not a harness call): `ls logs/ops61/`. If `ops61_ans4_x1_strip.msh` or `f05x.json` is missing, first run item 159's W0 export verbatim, relabelled `OPS-61-item194-W0-export`. If `postpro_f05x/` is missing, also re-run item 159's W2.
- **W1** (≈ 5 s):

  `scripts/testing/run_and_log.sh OPS-61-item194-W1-freqcfg "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 60 python3 scripts/probes/ops61_step1a_palace_ans4.py freq-configs'"`
- **W2** (≈ 25 s):

  `scripts/testing/run_and_log.sh OPS-61-item194-W2-palace-f10 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_f10.json'"`
- **W3** (≈ 25 s): the same with `f05x_f128.json`, labelled `OPS-61-item194-W3-palace-f128`.
- **W4** (the order-2 cost probe, ≤ 590 s): the same with `f05x_o2.json`, labelled `OPS-61-item194-W4-palace-o2`.
- **W5** (≈ 3 s), the probe's last window:

  `scripts/testing/run_and_log.sh OPS-61-item194-W5-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-set'"`
- **W6** (≈ 100 s; item 159's W6 103 s):

  `scripts/testing/run_and_log.sh OPS-61-item194-W6-ours-p1-f10 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_OPS61_DEGREE=1 FEM_EM_OPS61_FREQUENCY_MHZ=10 timeout -k 30 300 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ops61_palace_strip_comparison.py -s -v --tb=short'"`
- **W7:** the same with `=128`, labelled `OPS-61-item194-W7-ours-p1-f128`.
- **W8** (the closing window and default-path regression, ≈ 103 s): item 159's W6 command verbatim, with no frequency knob, relabelled `OPS-61-item194-W8-regress-p1`.
  - It must reproduce 117 206 cells, the (G) Z spreads 0.0747 / 0.0691 / 0.0897 %, σ_max 0.999760093664 to its printed digits, and ours vs Palace 0.0000 %.
  - It rewrites `ours_f05x_p1.json` with reproduced content. Item 192's (X) reads that file, and the drift is ≤ 5e-11.

**Traps.**
- Never overwrite `f05x.json`, `flipx.json`, `f1` / `f05` / `flip`, any existing `postpro_*`, or `ours_f05x_p2.json`.
- Palace aborts on a port set that is not aligned with the axes. An abort is a finding: never change `Direction` in-slot.
- **W4 is the first Palace order-2 run in this repo.** A Status 124 there is the finding (Palace order 2 does not fit a slot window at `-np 8`):
  - record it, and do not retry or raise `-np` in-slot;
  - check `docker compose -f docker/docker-compose.yml ps` for a lingering `fem-em-palace` container before W5;
  - W5–W8 still run, with o2 "not measured".
- Compare S classes only; never read Palace's `z_matrix_ohm`.
- The test module's default must stay byte-identical, which W8 proves. No edit to the probe after W5, or to the module after W8.
- Areas are rank-local.
- Record row: relative gaps and identity readings only. No absolute `S₁₁` / `Z_in`. Never set a Palace number against an AED one in any tracked file (CLAUDE.md, the Palace directive).

**Scope.** Covered: Palace order 1 at 10 / 128 MHz, Palace order 2 at 64 MHz, and our degree 1 at 10 / 128 MHz. Not covered:
- Palace order 2 at 10 / 128 MHz: the review prices it from W4.
- The PEC / copper hole columns (Palace's PEC-on-tag-401 and conductivity boundaries, step 2, which needs item 193's hole-mesh export).
- `ANS-7`'s lumped RLC.
- The (a) / (b) / (c) private reading, which is the weekly's.
- `src/` and any band.

**Status it can move.**
- The `OPS-61` row records step 1d. No glyph: 🟡 stays.
- **Rule (j):** the row's Done-when covers step 1c only, and step 1c's text says "Palace order 2 is not part of step 1c". The review adds a step-1d Done-when sentence in its commit.

**Negative result.**
- **Palace aborts, or a log does not show 117 206 elements and Rs 43.75 Ω/sq:** record verbatim, open a known-issues entry, BLOCKED (rule (d)).
- **(P) red on any run:** report the numbers per identity, park, BLOCKED. Never widen.
- **W8 not reproducing:** the module edit moved the default. Revert, BLOCKED.
- **(G) red at 10 or 128 MHz:** open a known-issues entry and park on `attempt/*`.
- **Matched-order gap ≫ 1e-3 at order 2:** this is a reading, not a failure. Print it and report; it is the item's purpose and the review reads it.
- No `PREDICTED RED:` line.

**Outcome (2026-10-03 06:00 CDT slot): DONE — step 1d recorded, no glyph.** Probe gained `freq-configs` and `read --exact-strip-set` (existing modes untouched); the test module gained `FEM_EM_OPS61_FREQUENCY_MHZ` (unset = item 159's path, proven by W8). Pre-check: strip mesh, `f05x.json` and `postpro_f05x/` present, so W0 was not run. Windows, all Status 0: W1 `20261003T110135Z_OPS-61-item194-W1-freqcfg.log` (1 s), W2 `…110141Z…W2-palace-f10.log` (21 s), W3 `…110208Z…W3-palace-f128.log` (35 s), W4 `…110249Z…W4-palace-o2.log` (268 s), W5 `…110728Z…W5-read.log` (3 s), W6 `…110833Z…W6-ours-p1-f10.log` (102 s, 16 passed), W7 `…111026Z…W7-ours-p1-f128.log` (101 s, 16 passed), W8 `…111214Z…W8-regress-p1.log` (108 s, 16 passed). File times precede their last windows (probe 06:01:32 < W5 06:07:28; module 06:08:25 < W6 06:08:33).
- **Palace log check:** each of W2–W4 shows 117 206 elements and Rs 4.375e+01 Ω/sq on attributes 21–24. o1 runs: ND (p = 1) 140 460, 19.3 s / 32.9 s Palace total, 1.8 G / 1.9 G peak (predicted ≈ 25 s, 1.8 G). **o2 cost probe:** ND (p = 2) 753 062 (5.36× o1; ours p2 on the sibling mesh 746 144), 265.6 s Palace total (four solves ≈ 64 s each), 6.1 G peak — inside the predicted 60–400 s, below the predicted 10–40 G.
- **(P), W5:** `f05x_f10` reciprocity 2.125e-10, σ_max 0.999993740607, spreads 0.0877 / 0.0326 / 0.0499 %, separation 106.2×; `f05x_f128` 9.949e-11, 0.999136167642, 0.0674 / 0.1369 / 0.1420 %, 413.1×; `f05x_o2` 6.663e-11, 0.999819887373, 0.0371 / 0.0248 / 0.0500 %, 845.0×. Column powers ≤ 0.8703. The predicted spread ceiling ≲ 0.15 % held.
- **Ours vs Palace at matched order (printed):** degree 1 vs o1, κ off, at 10 MHz (`W6.log:1976`) and 128 MHz (`W7.log:1977`): **0.0000 %** in every class; σ_max 0.999993740606 vs 0.999993740607, and 0.999136167640 vs …642. **Degree 2 vs o2 at 64 MHz (`W5.log:56`): 0.0000 %** in every class (σ_max 0.999819887372 vs …373). In-run κ vs Palace: 0.0012 / 0.0004 / 0.0003 % (10 MHz), 0.0026 / 0.0026 / 0.0016 % (128 MHz).
- **C2′ (predicted ≥ 25×):** ours p1 vs Palace o2 6.4132 / 2.5730 / 4.6383 % (predicted ≈ 6.67 / 2.54 / 4.52; self lower because the denominator is now o2's mean); separation 2.838e8 / 2.408e8 / 4.696e8×, so the matched-order gap is ≈ 2e-10 relative. Palace o2 vs o1 6.6674 / 2.5442 / 4.5224 %, the same as ours p2 vs p1 to the printed digit.
- **(F)/(G) at 10 / 128 MHz:** 117 206 cells; every N_i 40 facets, 5.6e-5 and w 7.0e-3 to ≤ 1.6e-15. 10 MHz κ off: reciprocity 2.2e-15, σ_max 0.999993740606, Z spreads 0.0818 / 0.0397 / 0.0433 %, separation 113.8×. 128 MHz κ off: 5.0e-16, 0.999136167640, 0.0816 / 0.1048 / 0.1669 %, 351.4×. κ in-run green at both.
- **(K) printed, not asserted (rule (e)):** 1.191884e-05 at 10 MHz and 2.570109e-05 at 128 MHz; the in-run κ is pooled 3.79e-6 and 4.34e-5.
- **Width confound at 10 MHz (printed):** ours strip vs ragged `LEG_D_S_MATRIX_10MHZ` 3.3926 / 0.9969 / 0.9409 %; not available at 128 MHz.
- **W8:** reproduced item 159 digit for digit: 117 206 cells, Z spreads 0.0747 / 0.0691 / 0.0897 %, σ_max 0.999760093664, ours vs Palace 0.0000 %, κ 1.5514965581e-05. Only reciprocity moved, at round-off (3.7e-16 vs 4.5e-16). `ours_f05x_p1.json` rewritten with the reproduced content.
- **For the review:** at both orders and all three frequencies, ours and Palace assemble the same discrete system. On this fixture the lumped-port and impedance-sheet bookkeeping are therefore independently confirmed. The order-1 → order-2 self-class shift is discretisation error that the two codes share. Palace order 2 at 10 / 128 MHz is priced by W4: ≈ 270 s and ≈ 6 G per frequency at `-np 8`.

### §9 item 229 — 2026-10-04 (chain F21)

*Drafted by a general-purpose agent for the 2026-10-04 10:30 daily review; the review reviews it.*

**`OPS-61` step 1e (§10 F21): Palace order 2 on `f05x` at 10 and 128 MHz.** Palace's own (P) identities are asserted on both new runs. Ours at degree 2 is set against Palace order 2 per class, printed, from our **tracked** degree-2 record on the same mesh. Ours needs no new solve.
- **Executor:** implementer. Probe only; no test module, no `src/`.
- **Tier:** heavy. Palace `-np 8`; the read at `-n 1`. Complex container for the read.
- **Class:** chain step (solve), measured median ≈ 18 min.
- **Cost:** ≈ 18 measured / ≈ 40 predicted min (15 fixed + ≈ 8 code + W1 < 1 + W2 ≈ 3–4 + W3 ≈ 6–8 + W4/W5 < 1 + ≈ 5 spare).
- **Depends on:** nothing. Every input is on disk or tracked:
  - `logs/ops61/` holds `ops61_ans4_x1_strip.msh`, `f05x.json`, `postpro_f05x_f10/`, `postpro_f05x_f128/` and `postpro_f05x_o2/` (checked 2026-10-04).
  - Ours is `ANS4_EXACT_STRIP_S_P2_X1` / `ANS4_EXACT_STRIP_S_P1_RECORDS` in `tests/validation/test_ans4_exact_strip_records.py:134,226`.

**Why the comparand is a record, not a solve.** `ANS4_EXACT_STRIP_S_P2_X1` is item 192's degree-2 κ-off 4×4 at 10 / 64 / 128 MHz on the exact-rectangle ×1 rung: 117 206 cells, `-n 8`, `20261003T094036Z_ANS-4-item192-W3-p2.log:2199–2216`. Palace's mesh is that rung exported (`ANS4_EXACT_STRIP_CELL_RECORD`, F10a W3).
- Item 192's (X) tied its 64 MHz S to `ours_f05x_p2.json` within `STEP3A_RECORD_RTOL` 1e-6.
- Item 194 W5 printed `ours_f05x_p2.json` vs Palace o2 at **0.0000 %** (≈ 2e-10).
- So the record is the discrete system Palace solves. W4 re-proves that at 64 MHz as a printed tie.
- Re-solving ours at degree 2 would cost ≈ 280 s per frequency (item 166 W1), ≈ 500 s on today's 1.4–2× box, and add nothing the record does not carry.

**Code (probe `scripts/probes/ops61_step1a_palace_ans4.py` only; additive; disclosed).**
- **`FREQ_SET_O2`** = `{"f05x_o2_f10": (0.010, 2), "f05x_o2_f128": (0.128, 2)}`, beside `FREQ_SET` (`:671`). `FREQ_SET` itself is not edited.
- **New `o2-configs` mode.** A new function mirroring `freq_configs` (`:679–704`) with the same checks on `f05x.json`:
  - checks: `Model.Mesh`, four ports with Attributes `[21]…[24]`, all `+Z`, `Order == 1`, `Freq == [0.064]`;
  - writes `f05x_o2_f10.json` and `f05x_o2_f128.json` only (`Order 2`, `Freq` set, `Output "postpro_<name>"`);
  - asserts neither file exists (never overwrite).
- **New `read --exact-strip-o2` mode.**
  - Lift the per-run (P) block of `read_exact_strip_set` (`:735–786`) into `_p_check(name, path, tag)`, returning `(s, ident, fails)` and printing with the caller's tag. `read_exact_strip_set` calls it with its existing `[OPS-61 1d W5]` tag, so its output stays byte-identical. W5 proves that.
  - The new function runs `_p_check` on each present `FREQ_SET_O2` run, and on `f05x_o2` as the 64 MHz tie.
  - It imports the two record constants from `tests.validation.test_ans4_exact_strip_records`. The probe already puts `REPO` on `sys.path` (`:59–61`) and already imports `tests.validation.*` in this mode family.
  - Class means use the probe's `_circulant_classes`, as `ours_means` does (`:724–732`).
  - It never reads `z_matrix_ohm`, the `aed` key or any file under `docs/private/` / `aed_results/`.
- **Dispatch:** two new `sys.argv` branches beside `:837–846`. `export`, `read`, `read --exact-strip`, `flip-config`, `freq-configs` and `read --exact-strip-set` keep their behaviour.

**Anchors (asserted; bands imported, never restated).**
- **(P)** on Palace `f05x_o2_f10` and `f05x_o2_f128`, item 159's set, as `read_exact_strip_set` applies it:
  - reciprocity ≤ `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`);
  - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every column power ≤ 1 + 1e-9;
  - every S-class spread ≤ `ADJACENT_SPREAD_BAND` (`test_port_birdcage_lumped_column.py:121`);
  - pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst (`test_port_birdcage_four_port.py:129,173`).
- **Palace log check (executor, journalled).** Each of W2 / W3 shows:
  - 117 206 elements;
  - `Rs = 4.375e+01 Ω/sq` on attributes 21–24;
  - `ND (p = 2): 753062`;
  - `ω/2π = 1.000e-02` / `1.280e-01` GHz;
  - every `GMRES solver converged` line (four per run; `Linear.MaxIts` is 100).
- **W5 (regression of the lift):** `read --exact-strip-set` reproduces `20261003T110728Z_OPS-61-item194-W5-read.log`'s (P) lines and order-2 readings digit for digit.

**Negative controls (rule (e)).**
- **(C2″) Order-mismatch separation, predicted ≥ 25×, printed.** Per class: \|ours p1 (`ANS4_EXACT_STRIP_S_P1_RECORDS[k]`) − Palace o2\| / \|ours p2 (`…_S_P2_X1[k]`) − Palace o2\|.
  - The numerator is predicted to be our own p2-vs-p1 move: 1.98 / 0.43 / 0.72 % at 10 MHz, 5.85 / 4.50 / 7.22 % at 128 MHz (`…item192-W3-p2.log:2222,2224`).
  - The denominator is predicted ≤ 1e-6. 64 MHz measured 2.4e8–4.7e8× (item 194 W5).
  - The ceiling is set by the denominator's round-off, so ≥ 25× is ≥ 10⁶× under what the 64 MHz case showed.
  - The 10 MHz adjacent class (0.43 %) is the least separated numerator. It still clears 25× unless the matched gap exceeds 1.7e-4.
- **(T) 64 MHz tie, predicted ≤ 1e-6 relative per class, printed.** `ANS4_EXACT_STRIP_S_P2_X1["64"]` vs `postpro_f05x_o2`. A miss means the record is not Palace's discrete system: stop and report, and do not read 10 / 128.
- **(P)'s power to discriminate** was shown on this mesh by item 159's flip: 200.0204 % spreads against the 0.5 % band, the sign-flip ceiling, 400×. It is not re-run.

**Predictions (written first, printed).**
- **Price**, from item 194 W4 (`20261003T110249Z_OPS-61-item194-W4-palace-o2.log:110,164–336,341–351`). o2 at 64 MHz: Palace total 265.6 s.
  - Fixed ≈ 85 s: setup 48.1, estimation ≈ 33, the rest ≈ 4.
  - Iteration-bound ≈ 180 s: linear solve, preconditioner and coarse solve over 4 × 37 GMRES iterations, i.e. ≈ 4.87 s per iteration-of-four.
  - Peak 6.1 G.
- **Iterations** scale with frequency as at order 1: 16–17 at 10 MHz, 31–32 at 64, 53–55 at 128 (item 194 W2 / F10b W2 / item 194 W3). The o2/o1 ratio at 64 MHz is 1.17.
  - ⇒ o2 **10 MHz** ≈ 19 iterations (15–25), Palace total ≈ 180 s (160–210).
  - ⇒ o2 **128 MHz** ≈ 63 iterations (55–75), Palace total ≈ 390 s (350–450).
  - Peak 6.1–6.5 G on both, against the 128 G service limit.
- **(P):** spreads ≲ 0.15 % (o2 64 MHz: 0.0371 / 0.0248 / 0.0500 %); reciprocity ≈ 1e-10; σ_max < 1.
- **Ours p2 vs Palace o2** at 10 and 128 MHz: **0.0000 %** per class (≤ 1e-6), the 64 MHz result.
- **Palace o2 vs Palace o1** (`postpro_f05x_f10` / `_f128`) equals ours p2 vs p1 to the printed digit: 1.9821 / 0.4323 / 0.7216 % (10 MHz) and 5.8471 / 4.5018 / 7.2217 % (128 MHz).

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`).
- **Pre-check** (not a harness call): `ls logs/ops61/`.
  - If `f05x_o2_f10.json` or `f05x_o2_f128.json` already exists, stop: never overwrite. Report it.
  - If `postpro_f05x_o2/` is missing, (T) is "not measured"; continue.
- **W1** (≈ 2 s):
  `scripts/testing/run_and_log.sh OPS-61-item229-W1-o2cfg "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 60 python3 scripts/probes/ops61_step1a_palace_ans4.py o2-configs'"`
- **W2**, 10 MHz o2: **the price probe.** Predicted 160–210 s:
  `scripts/testing/run_and_log.sh OPS-61-item229-W2-palace-o2-f10 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_o2_f10.json'"`
- **Price gate (binding, after W2).** Read W2's `Total` (Elapsed Time Report) T₁₀ and its iteration count n₁₀.
  - The box factor is r = T₁₀ / (85 + 4.87 × n₁₀) (≈ 178 s at n₁₀ = 19). The projected 128 MHz total is r × 450 s, the upper prediction.
  - **Run W3 iff r × 450 ≤ 540 s**, i.e. r ≤ 1.2 (T₁₀ ≤ ≈ 215 s at 19 iterations). Print r, n₁₀ and the projection in the journal.
  - **T₁₀ > 215 s** ⇒ do **not** run W3. 128 MHz o2 is "not measured — price". W4 / W5 still run, and the item lands the 10 MHz half as a record.
- **W3**, 128 MHz o2. Predicted 350–450 s on a normal box. The same command with `f05x_o2_f128.json`, labelled `OPS-61-item229-W3-palace-o2-f128`. After any kill, check `docker compose -f docker/docker-compose.yml ps` for a lingering `fem-em-palace` container before W4.
- **W4** (≈ 3 s), the new read; the probe's last window together with W5:
  `scripts/testing/run_and_log.sh OPS-61-item229-W4-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-o2'"`
- **W5** (≈ 3 s), the regression of the lift: item 194's W5 command verbatim, labelled `OPS-61-item229-W5-regress-read`. **No probe edit after W4 starts.**

**Traps.**
- Never overwrite `f05x*.json`, `flipx.json`, `f1` / `f05` / `flip`, any existing `postpro_*`, or `ours_f05x_p*.json`.
- `freq-configs` asserts its outputs do not exist, so **do not re-run it**. Item 194's files are present and it would fail.
- An abort on a port set not aligned with the axes is a finding: never change `Direction`, `Order`, `MaxIts` or `Tol` in-slot.
- A `GMRES … did not converge` line is a finding, not a retry: report the iteration count, and treat that run's S as not measured.
- Compare S classes only; never read Palace's `z_matrix_ohm` or the `aed` key.
- Palace's `postpro_*` and the configs are gitignored under `logs/`, so nothing new is tracked. The harness logs under `docs/testing/logs/` are.
- Record row and journal carry relative gaps and identity readings only, plus Palace's own numbers in full (version, config, mesh cited). No absolute `S₁₁` / `Z_in` claim, and no Palace figure set against an AED one in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.
- Status 124 on W2 or W3: the price is the finding. Never lengthen the timeout or raise `-np` in-slot (no measured `-np 12` Palace point exists).

**Scope.** Covered: Palace order 2 at 10 / 128 MHz on the exact strip mesh, its (P) set, and ours-vs-Palace from the tracked records. Not covered:
- our degree-2 (F) / (G) at 10 / 128 MHz, re-solved (item 192 asserted them);
- the PEC / copper hole columns (F22, step 2);
- `ANS-7`'s lumped RLC;
- the private (a) / (b) / (c) reading, which is the weekly's;
- any h-statement: one mesh, two orders;
- `src/` and any band.

**Status it can move.**
- The `OPS-61` row records step 1e. No glyph: 🟡 stays.
- The row's §2-facing sentence can become "ours matches Palace to 0.0000 % at matched order in every class at 10 / 64 / 128 MHz, both orders". That holds only if W2 and W3 both land and every printed matched gap is ≤ 1e-6.
- **Rule (j):** the Done-when names step 1c only. The review adds a step-1e sentence in its commit.

**Negative result.**
- **Palace aborts, or a log lacks 117 206 elements / Rs 43.75 Ω/sq / ND 753 062:** record verbatim, open a known-issues entry, BLOCKED (rule (d)).
- **(P) red on either run:** report every identity's number, park, BLOCKED. Never widen.
- **(T) misses ≤ 1e-6:** the record and Palace's mesh have parted. Report the per-class gap and BLOCKED. Do not print 10 / 128 against the record.
- **Matched gap at 10 or 128 MHz ≫ 1e-6 with (T) green:** a reading, not a failure. Print it; the review reads it. It would be the first frequency-dependent disagreement between the two codes.
- **W5 not reproducing item 194's digits:** the lift moved the old mode. Revert the probe, BLOCKED.
- **Price gate trips:** a partial record (10 MHz only). Not BLOCKED; the review re-queues 128 MHz at a quiet hour.
- No `PREDICTED RED` line.

**Outcome (2026-10-04 16:30 CDT slot): DONE — step 1e recorded, no glyph.** The probe gained `FREQ_SET_O2`, `o2-configs` and `read --exact-strip-o2`. The (P) block of `read_exact_strip_set` was lifted into `_p_check(name, path, tag)`. Pre-check: neither o2 config existed and `postpro_f05x_o2/` was present. `freq-configs` was not re-run. Windows, all Status 0: W1 `20261004T213348Z_OPS-61-item229-W1-o2cfg.log` (0 s), W2 `20261004T213356Z_OPS-61-item229-W2-palace-o2-f10.log` (213 s), W3 `20261004T213741Z_OPS-61-item229-W3-palace-o2-f128.log` (351 s), W4 `20261004T214342Z_OPS-61-item229-W4-read.log` (3 s), W5 `20261004T214346Z_OPS-61-item229-W5-regress-read.log` (2 s). The probe was last edited before W1. Palace is v0.18.1, configs `logs/ops61/f05x_o2_f{10,128}.json`, mesh `ops61_ans4_x1_strip.msh`.
- **Palace log check:** W2 and W3 each show 117 206 elements (`:68`), Rs = 4.375e+01 Ω/sq on attributes 21–24 (`:79–82`), `ND (p = 2): 753062` (`:110`) and ω/2π = 1.000e-02 / 1.280e-01 GHz. Each has four `GMRES solver converged` lines: 26 iterations at 10 MHz (W2 `:153,194,235,276`) and 55 at 128 MHz (W3 `:182,252,322,392`). No line reports a failure to converge.
- **Price gate:** T₁₀ = 211.423 s Palace total (W2 `:307`), n₁₀ = 26. Both are just above prediction: 15–25 iterations and 160–210 s. r = 211.423 / (85 + 4.87 × 26) = 0.999, and the projection is r × 450 = 449.6 s ≤ 540 s, so W3 ran. W3 measured 349.3 s Palace total (`:423`) at 55 iterations, against a prediction of 350–450 s and 55–75 iterations. Peak was 6.0 G and 6.3 G, against a prediction of 6.1–6.5 G.
- **(P), W4:** `f05x_o2_f10` reciprocity 1.390e-10, σ_max 0.999995411583, spreads 0.0963 / 0.0270 / 0.0619 %, separation 99.51×. `f05x_o2_f128` reciprocity 3.261e-11, σ_max 0.999307666282, spreads 0.0167 / 0.0319 / 0.0329 %, separation 1834.92×. Column powers are ≤ 0.8835. The predicted ceiling of spreads ≲ 0.15 % held. `f05x_o2` was re-checked green, with numbers identical to item 194.
- **(T) 64 MHz tie:** `ANS4_EXACT_STRIP_S_P2_X1["64"]` vs `postpro_f05x_o2` reads 2.464e-10 / 9.668e-11 / 9.883e-11, against a prediction of ≤ 1e-6. GREEN.
- **Ours p2 record vs Palace o2, matched order:** **0.0000 %** in every class. At 10 MHz the gaps are 2.542e-10 / 9.997e-11 / 6.236e-11; at 128 MHz they are 9.930e-11 / 8.685e-11 / 1.495e-10. That is the record's 10-digit rounding, and the frequency-dependent disagreement is nil.
- **C2″ (predicted ≥ 25×):** ours p1 record vs Palace o2 reads 1.9923 / 0.4335 / 0.7239 % at 10 MHz and 5.6659 / 4.6653 / 7.7237 % at 128 MHz. Separations are 7.836e7 / 4.336e7 / 1.161e8× and 5.706e8 / 5.371e8 / 5.167e8×.
- **Palace o2 vs Palace o1** reads 1.9821 / 0.4323 / 0.7216 % (10 MHz) and 5.8471 / 4.5018 / 7.2217 % (128 MHz). That equals ours p2 vs p1 to the printed digit, as predicted.
- **W5:** the output section is byte-identical to `20261003T110728Z_OPS-61-item194-W5-read.log` (a `diff` of the two output sections is empty), so the lift did not move the old mode.
- **For the review:** at both orders and all three frequencies, ours and Palace solve the same discrete system on this mesh. The §2-facing sentence the brief licences now holds. Still one mesh: this is not an h-statement.

### §9 item 234 — 2026-10-04 (chain F22)

*Drafted by a general-purpose agent for the 2026-10-04 18:00 daily review; the review reviews it.* Licence: §10 chain F22 as the 2026-10-03 weekly wrote it ("a Gmsh export of `ANS-6`'s exact hole mesh and Palace's PEC column on it at 64 MHz, order 1 (order 2 after F19). Target: Palace. The first Palace run on a hole mesh — a capability probe, red is a record"). F19's 64 MHz degree-2 record has landed (item 226), so the order-2 leg runs here behind a price gate. No band, tolerance, Done-when or validation target moves. The sheets-and-κ ruling (§9, 2026-10-02) binds: exact CAD rectangles, κ off.

**`OPS-61` step 2a (§10 F22): the exact-rectangle hole mesh exported to Gmsh 2.2, and Palace's PEC column on it at 64 MHz.** Palace's own (P) identities are asserted. Ours is set against Palace per class from our **tracked** records on the same mesh at matched order, printed. Ours needs no new solve.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Export and read at `-n 1` in the complex container; Palace at `-np 8`.
- **Class:** chain step (solve), measured median ≈ 23 min.
- **Cost:** ≈ 23 measured / ≈ 40 predicted min (15 fixed + ≈ 10 code + W1 ≈ 1 + W2 / W3 ≈ 1 + W4 ≈ 3–5 + W5 / W6 < 1 + ≈ 5 spare).
- **Depends on:** nothing. Every input is on `main` or on disk:
  - the hole route with exact rectangles: `_sheets_build(True, as_hole=True, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`, the call `th14._build_ladder` makes (`tests/validation/test_th14_birdcage_copper.py:313`; item 193);
  - `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` = 80 798 (`tests/validation/test_ans6_exact_strip.py:100`);
  - ours, degree 1: `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]` (`:119–121`; `-n 2`, `20261003T095650Z_ANS-6-item193-W2-solve.log:3096`, reproduced to 10 digits by item 193 W3 and items 216 / 226 W0);
  - ours, degree 2: `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["pec"]` (`:154–156`; `-n 8`, `20261004T184757Z_ANS-6-item226-W2a-p2-64-sweeps.log:1256–1257`, reproduced by W2b);
  - `logs/ops61/f05x.json` as the config template, and `logs/ops61/postpro_f05x/` for one printed cross-fixture line.

**What "the PEC column" is on our side, and how it maps to Palace.**
- Our PEC column is `TH-15` step 3's fixture verbatim: `_problem(..., pec_facet_tags=None)` pins **every exterior facet** (`test_th14_birdcage_copper.py:339`). On the hole mesh those are the outer box faces *and* the cavity wall, facet tag 401 (`BIRDCAGE_CONDUCTOR_SURFACE_TAG`, `src/fem_em_solver/io/mesh.py:266`). The coil is a hole: there are no conductor cells (tag 1 absent, `test_birdcage_conductor_hole.py` (iv)).
- Palace's equivalent is `Boundaries.PEC.Attributes = [1, 401]`, with the box written as attribute 1 (`EXTERIOR_TAG`, probe `:67`) and the wall as attribute 401.
- **Why step 1a's writer cannot be reused as is.** `export()` writes every `exterior_facet_indices` triangle as attribute 1 (`:205–208`). On a hole mesh that set includes the wall. So the wall would be PEC by accident, and the 0.3072 m² exterior anchor (`:241–245, :374–376`) would fail. The split is `th14._outer_box_tags` (`:133`): box = exterior ∖ 401, with a reduced census (`exterior`, `cavity`, `cavity_exterior`, `outer`, `tag_collision`).
- **Every exterior triangle carries an attribute.** Palace applies its natural condition (PMC, n × H = 0) to any exterior boundary with no listed BC. An unlisted wall is therefore a silent magnetic wall, not a conductor. That is exactly control (K′) below, run on purpose.
- **Ports.** `th14._build_ports` (`:185`) is the construction our PEC column was solved with. Its returned `tags_f` holds exactly N_i, the 0.007 × 0.008 SPEC rectangle (item 193 (M2)). The export writes N_i as 20 + i. The template's four `LumpedPort` entries are unchanged: Attributes `[20+i]`, `R` 50 Ω (`REFERENCE_IMPEDANCE_OHM` = `TERMINATED_PORT_IMPEDANCE_OHM`, `test_port_birdcage_four_port.py:125`), `+Z`. So Palace reads Rs = 50 · 0.007 / 0.008 = 43.75 Ω/sq, as on `f05x`. F_i ∖ N_i is not written, because no BC reads it (`f05x`'s 30 + i triangles were inert).
- **Materials** are built from the mesh's own cell tags: phantom (σ = 0.5, εᵣ = 78) plus every other tag as air. Never list `CONDUCTOR_CELL_TAG`, which is absent from a hole mesh. Whether Palace aborts on a material attribute it cannot find is unmeasured, so do not test it.

**Code (probe only; additive; disclosed).**
- **Constants** beside `FREQ_SET_O2` (`:679`):
  - `HOLE_MESH_NAME = "ops61_ans6_x1_hole_strip.msh"`;
  - `HOLE_WALL_ATTR = BIRDCAGE_CONDUCTOR_SURFACE_TAG` (imported, never restated);
  - `HOLE_SET = {"h05x_pec": ([EXTERIOR_TAG, HOLE_WALL_ATTR], 1), "h05x_pmc": ([EXTERIOR_TAG], 1), "h05x_pec_o2": ([EXTERIOR_TAG, HOLE_WALL_ATTR], 2)}`, mapping name to (PEC attributes, Order).
- **New `export_hole_pec()`, dispatched by `export --hole-pec`.** This is a new function. **`export()` is not edited and not lifted**: the only way to regression-run it is `export --exact-strip`, which rewrites F10a's strip mesh and `f05x.json` without an existence check.
  1. `assert comm.size == 1`. Then check `logs/ops61/` and assert that `HOLE_MESH_NAME` and the three `HOLE_SET` JSONs do not exist (never overwrite). Load `f05x.json` and check it as `o2_configs` does (`:866–876`): `Model.Mesh`, four ports with Attributes `[21]…[24]`, all `+Z`, `Order == 1`, `Freq == [0.064]`.
  2. Build `_sheets_build(True, as_hole=True, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`. Then `create_connectivity(fdim, tdim)` and `create_entity_permutations()`, as `_build_ladder` does (`:313–317`).
  3. `outer_tags, census = th14._outer_box_tags(msh, facet_tags, comm)` and `tags_f, _defs, specs = th14._build_ports(msh, cell_tags, comm)`.
  4. Write the triangle blocks with step 1a's `entities_to_geometry` pattern (`:196–202`): box → 1, `facet_tags.find(401)` → 401, N_i → 20 + i. Write the tets with step 1a's reorientation and Gmsh 2.2 writer, copied (`:250–283`). Assert no tet is untagged.
  5. Write three configs from the `f05x.json` template, changing only these keys: `Model.Mesh` → `HOLE_MESH_NAME`; `Domains.Materials` → as above; `Boundaries.PEC.Attributes` → the `HOLE_SET` list; `Solver.Order`; `Problem.Output` → `postpro_<name>`. Print each config's `Boundaries` and `Domains` blocks.
  6. Anchors (E1)–(E5) below are asserted after writing, as step 1a does, so a failure still leaves the artefacts to inspect.
- **New `read_hole_pec()`, dispatched by `read --hole-pec`.**
  - It runs `_p_check(name, path, "[OPS-61 2a W5]")` (`:714`) on `h05x_pec` (required) and on `h05x_pec_o2` if present. It does the same on `h05x_pmc` but **prints only**: the PMC wall is a control, not a fixture.
  - It imports the two ours records and `ANS4_EXACT_STRIP_S_P1_RECORDS` from tests. Class means come from `_circulant_classes`, as `read_exact_strip_o2`'s `means` does (`:914–915`).
  - If the gitignored `logs/ans6_exact/ans6_exact_p1.json` is present, it also prints the 64 MHz PEC 4×4's entry-level max\|ΔS\| and σ_max beside Palace's, labelled "gitignored source, printed only".
  - It never reads `z_matrix_ohm`, the `aed` key or any file under `docs/private/` / `aed_results/`.
- **Dispatch:** two new `sys.argv` branches beside `:992–1005`, plus a docstring paragraph. Every existing mode keeps its behaviour, which W6 proves.

**Anchors (asserted; bands imported, never restated).**
- **Export, W1:**
  - **(E1)** tets = `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` exactly (as step 1b's export asserted its record).
  - **(E2)** `CONDUCTOR_CELL_TAG` is absent from the cell tags, and the phantom tag is present.
  - **(E3)** census: `cavity_exterior == cavity > 0` (the wall is all exterior), `outer == exterior − cavity`, `tag_collision == 0`.
  - **(E4)** box (attribute 1) area = `EXTERIOR_AREA_M2` 0.3072 within `EXTERIOR_AREA_BAND` (`:70–71`).
  - **(E5)** for each port: the numpy area of N_i, and the spec's `sheet_width_m` and `gap_height_m`, equal `SPEC_AREA_M2` / `SPEC_WIDTH_M` / `SPEC_HEIGHT_M` (`test_ans4_exact_strip_records.py:131–133`) within `STRIP_BAND`, relative, as (M2) applies it.
  - Printed only: the wall's triangle count and area, with the area's ratio to `diag["conductor_cavity_cad"]["area_m2"]` (an inscribed triangulation of tori and cylinders, so a little under 1), and the N_i facet count per port.
- **(P), W5:** item 159's set via `_p_check`, on `h05x_pec` and on `h05x_pec_o2` if run:
  - reciprocity ≤ `RECIPROCITY_BAND`;
  - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every column power ≤ 1 + 1e-9;
  - every S-class spread ≤ `ADJACENT_SPREAD_BAND`;
  - pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst.
- **Palace log check (executor, journalled).** W2 / W4 each show:
  - 80 798 elements, and the bounding box ±0.12 / ±0.12 / ±0.10 m;
  - `Rs = 4.375e+01 Ω/sq` on attributes 21–24;
  - `Configuring Dirichlet PEC BC at attributes:` listing **1 and 401** (W3: **1 only**);
  - `ND (p = 1): 111946` (W2 / W3) or `ND (p = 2): 570332` (W4). These are our degree-1 / degree-2 unknown counts on this mesh (`…item193-W2-solve.log:3089`; item 226). Palace counts every edge, as `f05x`'s 140 460 = its edge count showed (`20261002T021838Z_OPS-61-F10b-W2-palace-f05x.log:66,110`). A free geometry check;
  - `ω/2π = 6.400e-02 GHz`;
  - four `GMRES solver converged` lines (`Linear.MaxIts` 100).

**Ours against Palace (printed, never gated).** A same-mesh, same-order record of ours exists. Its band is not imported, though: `STEP3A_RECORD_RTOL` is a re-recording band for our own code, not a cross-solver one. So the readings are printed, as item 229 printed its own.
- **(T1)** `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]` vs Palace `h05x_pec`, \|Δ mean\|/\|mean\| per class. **Predicted ≤ 1e-6** ("0.0000 %"), because the same discrete system is solved. The record's 10-digit rounding sets the floor near 1e-10.
- **(T2)** if W4 runs: `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["pec"]` vs `h05x_pec_o2`, predicted ≤ 1e-6.
- **(C3)** Palace `h05x_pec` vs Palace `postpro_f05x` (σ = 800 solid, sibling mesh) per class, beside the same move between ours' two records (`ANS6_…_P1["64"]["pec"]` vs `means(ANS4_EXACT_STRIP_S_P1_RECORDS["64"])`). Predicted equal to the printed digit. This is the conductor-model move as each code sees it. Prints "absent" without `postpro_f05x/`.

**Negative controls (rule (e)).**
- **(K′) The wall attribute reaches Palace's operator — floor asserted, size predicted** (item 159 (K)'s pattern). `h05x_pmc` drops 401 from `PEC`, which makes the coil wall Palace's natural PMC wall: the broken export, run on purpose.
  - Asserted: every class mean moves from `h05x_pec` by > `KNOB_REACHES_THE_SOLVE` (`test_port21_feed_sensitivity.py:156`, imported).
  - Predicted, printed: ≥ 10 % on at least one class, because a magnetic wall where the current path was conducting rewrites the coil. No ours-PMC record exists, so the size is not banded.
  - The read also prints `h05x_pmc`'s own (P) readings, never asserted.
- **(C2‴) Order-mismatch separation — predicted ≥ 25×, printed** (item 229's C2″ on this fixture). Per class: \|ours p2 − Palace o1\| / \|ours p1 − Palace o1\|. If W4 ran, the reverse is also printed: \|ours p1 − Palace o2\| / \|ours p2 − Palace o2\|.
  - The numerator is predicted to be our own PEC p1-to-p2 move at 64 MHz, computed from the two tracked records: ≈ 4.35 / 1.67 / 3.11 % (item 226's (N3) Cu move, 4.35 / 1.67 / 3.11 %).
  - The denominator is (T1), predicted ≤ 1e-6. So the predicted separation is ≥ 1.6e4×, and ≥ 25× is the item-229 bar. The ceiling is set by the record's rounding.
- (P)'s power to discriminate was shown by item 159's flip (200.0204 % spreads against the band). It is not re-run.

**Predictions (written first, printed).**
- **W1:** hole mesh ≈ 20–25 s at `-n 1` (item 226 P: 23.30 / 22.49 s); window 45–90 s. Box 0.3072 to ≈ 1e-16. N_i facet count equal on all four ports.
- **W2 (o1, PEC):** Palace total 15–30 s (`f05x` 23.35 s on 140 460 ND; this mesh 111 946); peak ≈ 1.3–1.9 G. GMRES iterations printed, not predicted (a PEC hole changes the operator; `f05x` 31–32).
- **W3 (o1, PMC):** the same price.
- **(P):** spreads ≲ 0.15 %, reciprocity ≈ 1e-10, σ_max < 1.
- **W4 (o2):** `f05x`'s o2 / o1 ratio was 265.6 / 23.35 = 11.4 (item 194 W4, item 159 W2). That projects 11.4 × T₁, ≈ 170–340 s; peak 4–6 G.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`).
- **Pre-check** (not a harness call): `ls logs/ops61/`.
  - If `ops61_ans6_x1_hole_strip.msh` or any `h05x_*.json` exists, stop and report: never overwrite.
  - If `f05x.json` is missing, run item 159's W0 export verbatim first, relabelled `OPS-61-item234-W0-export`.
- **W1** (≈ 45–90 s), the export:

  `scripts/testing/run_and_log.sh OPS-61-item234-W1-export-hole "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py export --hole-pec'"`
- **W2** (≈ 15–30 s), PEC column, order 1:

  `scripts/testing/run_and_log.sh OPS-61-item234-W2-palace-hole-pec "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 h05x_pec.json'"`
- **W3** (≈ 15–30 s), the (K′) control: W2's command with `h05x_pmc.json`, labelled `OPS-61-item234-W3-palace-hole-pmc`.
- **Price gate (binding, after W2).** Read W2's `Total` (Elapsed Time Report) T₁ and its iteration count.
  - **Run W4 iff 11.4 × T₁ ≤ 450 s**, i.e. T₁ ≤ 39.5 s. Print T₁, the iterations and the projection in the journal.
  - Otherwise order 2 is "not measured — price". W5 / W6 still run, and the item lands the order-1 half.
- **W4**, order 2: W2's command with `h05x_pec_o2.json`, labelled `OPS-61-item234-W4-palace-hole-pec-o2`. After any kill, run `docker compose -f docker/docker-compose.yml ps` to look for a lingering `fem-em-palace` container before W5.
- **W5** (≈ 3 s), the new read; together with W6, the probe's last window:

  `scripts/testing/run_and_log.sh OPS-61-item234-W5-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --hole-pec'"`
- **W6** (≈ 3 s), the dispatch regression: item 229's W4 command verbatim (`20261004T214342Z_OPS-61-item229-W4-read.log:12`), labelled `OPS-61-item234-W6-regress-read`. Its output section must `diff` empty against that log's. **No probe edit after W1 starts.**

**Traps.**
- **Never run `export` or `export --exact-strip`.** They rewrite step 1a's / F10a's mesh and configs without an existence check. Do not re-run `flip-config`, `freq-configs` or `o2-configs`, which assert that their outputs do not exist and would fail.
- **The first `as_hole` build in this probe runs at `-n 1`** (the probe asserts it). A gmsh failure deadlocks wider runs.
- **Palace on a hole mesh is untried.** An abort is the finding: record the MFEM / Palace message verbatim. Never change `Direction`, attributes, `Order`, `MaxIts` or `Tol` in-slot. A `GMRES … did not converge` line is a finding, not a retry: that run's S is "not measured".
- Write the wall as its own attribute. Folding it into attribute 1 passes the PEC tie but fails (E4) and voids (K′).
- `cell_tags.values` and facet sets are rank-local. The export is serial by assertion; `_outer_box_tags` and `_build_ports` reduce internally.
- Compare S classes only. Never read Palace's `z_matrix_ohm` or the `aed` key.
- The record widths differ: P1 is `-n 2`, P2 is `-n 8`. Both are κ off. Never re-register a record from Palace.
- Palace's `postpro_*`, the configs and the mesh are gitignored under `logs/`. The harness logs under `docs/testing/logs/` are tracked.
- Record row and journal: identity readings and relative gaps, plus Palace's own numbers in full (version, config, mesh cited). No absolute `S₁₁` / `Z_in` claim, and no Palace figure set against an AED one in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.
- Status 124: the price is the finding. Never lengthen a timeout or raise `-np` in-slot.

**Scope.** Covered: the exact-rectangle hole mesh as a Gmsh 2.2 export with box / wall / port attributes; Palace's PEC column on it at 64 MHz, order 1, plus order 2 if the gate allows; its (P) set; ours vs Palace from the tracked records; the PMC-wall control. Not covered:
- the copper column: Palace's finite-conductivity / impedance boundary on 401 against our Leontovich term. That is step 2's other half, for the weekly to enumerate;
- 10 / 128 MHz;
- the F-human rung;
- any h-statement (one mesh);
- the private (a) / (b) / (c) reading;
- `ans:6`'s tracked artefacts;
- `src/` and any band.

**Status it can move.**
- The `OPS-61` row gains a "step 2a record (F22, §9 item 234)" sentence: Palace's first hole-mesh run, (P) readings, and the matched-order (T1) / (T2) per class. **No glyph: 🟡 stays.** The row's Done-when names step 1c only, and step 2 is "at closure".
- **Rule (j)**, proposed Done-when sentence for the review: "Step 2a is done when Palace `-np 8` solves the exported exact hole mesh's PEC column at 64 MHz, order 1 (80 798 elements, PEC on attributes 1 and 401, Rs 43.75 Ω/sq), its 4×4 meets item 159's (P) set at the imported bands, the PMC-wall control moves every class mean by more than `KNOB_REACHES_THE_SOLVE`, and ours-vs-Palace per class is printed against `ANS6_EXACT_STRIP_CLASS_MEANS_P1['64']['pec']`."
- If (T1) prints ≤ 1e-6, the §2 / CLAUDE.md Palace clause ("matches ours to 0.0000 % at matched order") can extend to the PEC hole route. That is the weekly's call.
- §10 F22 completes, and the copper half becomes enumerable.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(E1)–(E5) red:** the exported mesh is not the solved mesh. Stop before any Palace window, park on `attempt/*`, BLOCKED (rule (d)). If (E1) alone is red, the mesher moved: open a known-issues entry. Re-recording under (1*) is `record-reconciler`'s job, not this item's.
- **Palace aborts on `h05x_pec`, or the W2 log lacks 80 798 elements / Rs 43.75 / PEC on 1 and 401:** this is the capability finding — "red is a record". Record the message verbatim in the journal, open a known-issues entry ("Palace v0.18.1 on the exported hole mesh: <message>"), park the probe on `attempt/*` and mark the item BLOCKED (implementer-run.md step 4). The review writes the finding into the row as a record.
- **(P) red:** report every identity's number, park, BLOCKED. Never widen. Our (G) is green on this mesh at degree 1 (item 193), so this would be the first identity-level disagreement between the codes.
- **(K′) red** (the PMC wall moves S by ≤ 1e-6): attribute 401 is not reaching Palace's operator, so the PEC tie would be unexplained. Park, BLOCKED.
- **W6 not byte-identical:** the dispatch edit moved an old mode. Revert the probe, BLOCKED.
- **(T1) or (T2) ≫ 1e-6 with (P) and (K′) green:** a reading, not a failure. The item lands DONE as a record, and the review reads it. It would be the first disagreement between the codes. The ports, phantom and box were already shown to agree on `f05x`, so it would point at the PEC hole boundary.
- **W3 aborts or fails to converge:** (K′) "not measured". (C2‴) still carries the control; journal it and continue to W5.
- **Price gate trips:** a partial record (order 1 only). Not BLOCKED; the review re-queues order 2.


## §7 row moved verbatim 2026-10-01 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/OPS-61.md`)* **Palace as a second independent reference** — AWS's open-source curl-curl FEM (Nédélec, lumped ports, S-parameters, AMR, hypre AMS) run on the `ANS-4` fixture beside HFSS (operator directive 2026-09-20: keep our solver as the thing being validated, with two judges that share no code with it or with each other; **for the 2026-09-26 weekly to commission**, §5.4) | 🟡 **Step 1a 2026-09-29 (06:00 slot, §9 item 75):** the `ANS-4` ×1 rung exported from dolfinx to Gmsh 2.2 (`scripts/probes/ops61_step1a_palace_ans4.py`; W1 116 085 tets, Fᵢ areas nominal to 1e-9, exterior 0.3072 m² to 5.6e-17). Palace `-np 8` `f1` at 64 MHz, σ = 800, in 25 s wall with 1.8 G peak (139 140 ND dofs; its log confirms 116 085 elements, the ±0.12/±0.10 m bbox, 0.064 GHz and Rs = 87.5 Ω/sq = the full 0.014 × 0.008 sheet). Identities on Palace's 4×4 are green: reciprocity 8.480e-11, σ_max 0.999733, class spreads 0.071 / 0.056 / 0.037 % (`20260929T110544Z_OPS-61-step1a-W5.log:42,73`). The `flip` control matched D·S·D to 4.3e-11, with spreads 200.0 % as predicted. **Port-geometry answer for the weekly:** Palace **refuses** the ragged f = 0.5 set (`f05`) with `MFEM abort: Specified direction does not align sufficiently with bounding box axes (10.614020, …)` (`…110504Z_OPS-61-step1a-W3.log:86`), so a like-for-like step 1 needs an exact f = 0.5 rectangle in our mesh. The class-mean gaps printed against `ans4` (1.3 / 0.7 / 0.6 %) are confounded by port width: no verdict. **Step 0 image + smoke landed 2026-09-28 (interactive session, operator-authorised):** `docker/palace.Dockerfile` (`aa0168e`) builds Palace **v0.18.1** by its CMake superbuild (CPU; MUMPS + SuperLU_DIST; SLEPc) — not Spack, whose builtin package stops at 0.16.0 (spack/spack-packages#4447 open). Image `fem-em-palace:0.18.1` built in 948 s compile at `-j6`; smoke: the bundled `spheres` example, `palace -np 2`, exit 0, 27 s wall, 815 MB peak, binary reports `Git changeset ID: v0.18.1` (`20260929T041759Z_OPS-61-step0-smoke.log`). **Step 0 ✅ 2026-09-28:** the `fem-em-palace` compose service (operator-authorized mount, `palace` profile, 128 G) runs the same smoke through compose with the repo visible at `/workspace`, exit 0, 22 s (`20260929T042009Z_OPS-61-step0-service.log`). **Step 1 is queueable** — run it as `docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace …` under `run_and_log.sh`, `-np` ≤ 12. **Why:** `ANS-6` isolated a ≈ 3 % self-class offset on the four-port birdcage to the lumped-port model / our degree-1 mesh (it is identical on a PEC coil and a copper coil). HFSS is one reference; a second general-purpose code with an *independently written* lumped-port model says whether the offset is ours or a port-model convention HFSS shares with nobody. Unlike AED numbers, **Palace results are publishable** — open source, no licence clause — so this reference can sit in tracked `COMPARISON.md` columns and in gates, which the AED columns never can. **Step 0 — the service (ops, cost probe first):** a `fem-em-palace` compose service from Palace's official container image if one is published, else a build stage (Spack or CMake against MFEM/hypre) pinned to a release; smoke: Palace's own bundled example runs to completion in the container; `-n` ≤ 12 ranks; the price of the smoke printed. **Step 1 — the `ANS-4` fixture:** export the gated four-leg birdcage mesh (`birdcage_port_domain`, the ×1 rung) to a format Palace reads (Gmsh `.msh` via the generator's own gmsh session; element tags carried), write the Palace JSON config — PEC box faces, phantom σ = 0.5 / εᵣ = 78, coil σ = 800 solved inside (the `ANS-4` fixture) *and* a PEC-coil variant (the `ANS-6` PEC column), four lumped ports on the gap sheets at 50 Ω, driven at 10 / 64 / 128 MHz — and run it under `run_and_log.sh` at the tier the probe sets. **Readout:** the 4×4 S per frequency; the three C4 classes against ours (tracked, public) and against the private AED column (`docs/private/`); Palace's own reciprocity / passivity / C4-spread identities beside our gate bands. **Decision rule, pre-registered:** (a) Palace agrees with HFSS on the self class and both differ from us by ≈ 3 % ⇒ the offset is ours (`PORT-21`'s target, now with two references); (b) Palace agrees with *us* and differs from HFSS ⇒ the offset is HFSS's port convention, and the `ANS-4` self-class verdict is re-read; (c) Palace agrees with neither ⇒ the lumped-port model is convention-sensitive at this fixture and the case moves to a wave-port or a de-embedded comparison. **Step 2 (at closure):** Palace on the `ANS-6` copper column (its surface-impedance boundary) and, later, the F-human rung at degree 2 as a price and S cross-check. **Not a migration:** nothing in `src/` changes; Palace is a reference, like AED, never the solver. **Step 1b (F10a, written 2026-10-01 10:30 review, rule (j); §9 item 148):** `birdcage_port_domain(leg_sheet_strip_fraction=0.5)` makes the narrowed port a CAD rectangle; asserted N_i = 0.5·dx·dz and F_i = dx·dz to 1e-9, flag-off 116 085 reproduced, a version-tagged cell record; Gmsh 2.2 export `ops61_ans4_x1_strip.msh` + `f05x.json` written. *(**Step 1b DONE 2026-10-01 (16:30 CDT slot, §9 item 148):** F_i 1.12e-4 and N_i 5.6e-5 m² to ≤ 2.2e-15 relative on all four ports, width 7.000e-3, F∖N exact; flag-off 116 085 re-gated; flag-on 117 206 cells registered as `ANS4_EXACT_STRIP_CELL_RECORD`; control (flag off) N_i +4.202e-2; `20261001T220051Z_OPS-61-F10a-W2.log` (8 passed, `-n 2`); history `docs/planning/chunks/OPS-61.md`. F10b queueable.)* **Step 1c (F10b, serial on 1b):** Palace `f05x` at 64 MHz — Palace's reciprocity / passivity / class spreads asserted at our imported bands; ours-vs-Palace per class printed at degrees 1 and 2, κ-corrected and not. **Step 1c Done-when (F10b, written 2026-10-01 18:00 review, rule (j); §9 items 159 + 166):** Palace `-np 8` solves `f05x` at 64 MHz: 117 206 elements, Rs 43.75 Ω/sq, port i = N_i, `+Z`, 50 Ω. Its 4×4 meets `RECIPROCITY_BAND`, σ_max and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every S-class spread ≤ `ADJACENT_SPREAD_BAND` with pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst — each asserted with the band imported. Our four-port on the same strip mesh meets `PORT-11`'s (i) / (ii) / (iii′) at degree 1 (item 159) and degree 2 (item 166), κ off and `"in_run"`, with the fixture tie (cells, N_i to 1e-9) asserted. Ours-vs-Palace per-class \|Δ mean\|/\|mean\| is printed for all four of ours. Step 1c is a record and moves no glyph. The (a)/(b)/(c) reading needs the private column and is the weekly's. Palace order 2 is not part of step 1c. *(**Step 1c first half DONE 2026-10-01 (21:00 CDT slot, §9 item 159):** Palace `f05x` 117 206 elements / Rs 43.75 Ω/sq, (P) reciprocity 1.556e-10, spreads 0.085 / 0.078 / 0.092 %; ours at degree 1 on the same mesh (F) and (G) green κ off and in-run, (K) 2.25e-5; **ours κ-off vs Palace prints 0.0000 % in every class (σ_max 0.999760093664 vs 0.999760093668)** — the step 1a gaps were port width; in-run κ 1.55e-5, ≈ 680× below the ragged value; `20261002T021945Z_OPS-61-F10b-W4-read.log`, `20261002T022338Z_OPS-61-F10b-W6-ours-p1.log`. Item 166 runnable.)* 

**Outcome (2026-10-04 21:00 slot): ✅ DONE.** W1–W6 Status 0 (`20261005T021821Z` … `20261005T022234Z_OPS-61-item234-*`). (E1)–(E5) green; Palace order 1 (13.4 s) and order 2 (160.6 s, price gate passed) both meet (P); (K′) moves every class ≥ 100 %; (T1) / (T2) ≤ 2.2e-10 per class; W6 diff-empty. F22 complete; `OPS-61` stays 🟡.


### §9 item 241 — 2026-10-05 (chain F23, daily provisional)

*Drafted by a general-purpose agent for the 2026-10-05 03:00 daily review; the review reviews it.* Licence: §9 rule (6), a provisional next step on live chain F. Item 234's outcome said "the copper (impedance-boundary) half is now enumerable" (`docs/testing/attempts.md`, entry "2026-10-05 02:26 UTC", "Hypothesis for next"). Chain F's Palace steps carry "Target: Palace" (PROJECT_PLAN.md §10, F21 / F22). This item sets no target, and it moves no band, tolerance or Done-when. The sheets-and-κ ruling (§9, 2026-10-02) binds: exact CAD rectangles, κ off.

**`OPS-61` step 2b (§10 F23): Palace's copper column on the exact hole mesh at 64 MHz.** Item 234's export is reused unchanged. Attribute 401 (the conductor wall) carries Palace's finite-conductivity boundary for copper, σ = 5.8e7 S/m. The run is order 1, plus order 2 behind a price gate. Palace's own (P) identities are asserted. The comparison with ours on the same mesh at matched order comes from ours' **tracked** copper records and is printed. Ours needs no new solve.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Configs and the read run at `-n 1` in the complex container; Palace runs at `-np 8`.
- **Class:** chain step (solve), measured median ≈ 16 min (the 2026-10-04 18:00 restock, PROJECT_PLAN.md:3875).
- **Cost:** ≈ 16 measured / ≈ 35 predicted min. That is 15 fixed + ≈ 8 code + W1 < 1 + W2 < 1 + W3 / W4 ≈ 1 + W5 ≈ 3–6 + W6 / W7 < 1 + ≈ 5 spare.
- **Depends on:** nothing. Every input is on `main` or on disk:
  - item 234's artefacts in the gitignored `logs/ops61/`: `ops61_ans6_x1_hole_strip.msh` (80 798 tets, box = attribute 1, wall = 401, N_i = 20 + i), `h05x_pec.json` as the template, and `postpro_h05x_pec/port-S.csv` (plus `postpro_h05x_pec_o2/port-S.csv`) as the PEC column;
  - ours, degree 1: `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]` and `["pec"]` (`tests/validation/test_ans6_exact_strip.py:116–121`; `-n 2`, `20261003T095650Z_ANS-6-item193-W2-solve.log:3091–3101`);
  - ours, degree 2: `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["cu"]` and `["pec"]` (`:152–157`; `-n 8`, `20261004T184757Z_ANS-6-item226-W2a-p2-64-sweeps.log:1256–1257`);
  - `th14.COPPER_SIGMA` = 5.8e7 and `th14.SIGMA_LADDER` (`tests/validation/test_th14_birdcage_copper.py:109, :111`), `th14.FREQUENCIES_HZ["64"]` (`:107`, = `FREQUENCY_64_HZ`), and `surface_resistance_ohm` (`src/fem_em_solver/core/cavity.py:543`).

**Our Leontovich term, and the Palace equivalent.**
- **Ours.** The copper column is `th14._build_ladder`'s σ loop (`test_th14_birdcage_copper.py:345–353`). It runs on `imp_problem`, which pins only the outer box (`OUTER_BOX_TAG`, the exterior minus 401). The wall is left unpinned and given `_leontovich_term` (`:160–171`): γ ∮ (n × u)·(n × v̄) ds(401), with γ = jωμ₀/Z_s. This is Jin §5.8.3's third-kind condition, the same as `core/boundary_terms.py:13–30`.
- **Our Z_s** = (1 + j)·R_s, where R_s = `surface_resistance_ohm(ω, σ)` = √(ωμ₀/(2σ)) (`cavity.py:543–546`; ω in rad/s; μ₀ = 4π·10⁻⁷, `utils/constants.py:6`).
- **The equivalence.** δ = √(2/(ωμ₀μ_rσ)), so 1/(σδ) = √(ωμ₀/(2σ)) = R_s identically, and Z_s = (1 + j)/(σδ).
- **At 64 MHz, copper** (ω = 2π·64.0e6 = 4.0212e8 rad/s): δ = 8.2607e-6 m, R_s = 2.087161e-3 Ω, Z_s = 2.087161e-3 (1 + j) Ω. Ours logged exactly this value: `sigma = 5.8e+07 S/m: Z_s = 2.087161e-03(1+j) ohm` (`…item193-W2-solve.log:2911`).
- **The (K2) control at σ = 5.8e11:** R_s = 2.087161e-5 Ω, which is 100× smaller (`:2926`).
- **Palace's keyword.** It is **not confirmed from any file in the repo.** No probe-written config uses it; the resolved JSON Palace writes (`postpro_h05x_pec/h05x_pec_resolved.json`) echoes only the blocks a config used; and the image ships no schema. `docker/palace.Dockerfile` copies only `/opt/palace` plus `examples/`. The 2026-10-03 daily notes that "Palace's lumped RLC and conductivity boundaries have never run here" (`docs/planning/reviews/2026-10-03-daily.md:93`).
  - Palace's documentation (v0.18 config reference) describes a `Boundaries.Conductivity` array with entries `{"Attributes", "Conductivity" [S/m], "Permeability" (μ_r, default 1), optional "Thickness"}`. It implements the thick-conductor Z_s = (1 + j)/(σδ) above. That is the keyword this item uses.
  - `Boundaries.Impedance` (`Rs` / `Ls` / `Cs`) is **not** used. Whether Palace combines those in series or in parallel is unverified, so an equivalent `Rs` / `Ls` pair cannot be written down with confidence.
  - W2 confirms the keyword by file, or the item stops BLOCKED with the reading.
- **The Palace configs.** Copy `h05x_pec.json` and change only these keys:
  - `Boundaries.PEC.Attributes` → `[1]`;
  - add `Boundaries.Conductivity` = `[{"Attributes": [401], "Conductivity": σ, "Permeability": 1.0}]`. **Never write `Thickness` or `External`**: the thick limit is Leontovich's;
  - `Solver.Order`;
  - `Problem.Output` → `postpro_<name>`.
  - Attribute 401 goes in exactly one boundary block. Every exterior triangle still carries a BC, so Palace's PMC fallback warning must not appear (item 234 W3 shows what that warning looks like: `20261005T021915Z_OPS-61-item234-W3-palace-hole-pmc.log:52–56`).

**Code (probe only; additive; disclosed).**
- **New `cu_configs()`, dispatched by `cu-configs`.**
  1. `assert comm.size == 1` (or no MPI import at all).
  2. Assert that none of `h05x_cu.json`, `h05x_cu_o2.json`, `h05x_cu11.json` exists (never overwrite).
  3. Load `h05x_pec.json` and check it the way `o2_configs` checks its template (`:892–904`): `Model.Mesh == HOLE_MESH_NAME`, four ports with Attributes `[21]…[24]`, all `+Z`, `Order == 1`, `Freq == [0.064]`, and `PEC.Attributes == [EXTERIOR_TAG, _hole_wall_attr()]`.
  4. Write three configs with `CU_SET = {"h05x_cu": (COPPER_SIGMA, 1), "h05x_cu_o2": (COPPER_SIGMA, 2), "h05x_cu11": (SIGMA_LADDER[-1], 1)}`, mapping name to (σ, Order). Both σ values are imported, never restated.
  5. Print each config's `Boundaries` block.
  6. Then assert (Z) below, reading σ, μ_r and Freq back from each written JSON.
- **New `read_hole_cu()`, dispatched by `read --hole-cu`.**
  - It reuses `_p_check` (`:742`), `_identities` (`:418`) and the `rel` / `means` idioms of `read_hole_pec` (`:1260–1273`).
  - It prints each run's `postpro_<name>/<name>_resolved.json` `Boundaries.Conductivity` block, which is file evidence of the keyword and of the defaults Palace filled in.
  - If the gitignored `logs/ans6_exact/ans6_exact_p1.json` is present, it prints the entry-level max\|ΔS\| and σ_max of ours' 64 MHz `cu` 4×4 against `h05x_cu`, labelled "gitignored source, printed only".
  - It never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** two new `sys.argv` branches beside `:1389–1392`, plus a docstring paragraph. `export --hole-pec`, `read --hole-pec` and every older mode keep their behaviour, which W7 proves.

**Anchors (asserted; bands imported, never restated).**
- **(Z), W1: the surface-impedance tie, pure arithmetic.**
  - For each `CU_SET` config: \|Z_ours − Z_doc\| / \|Z_ours\| ≤ `STRIP_BAND` (`tests/validation/test_ops61_exact_strip_ports.py:58`, imported), applied as a relative tie the way item 234's (E5) applies it.
  - Z_ours = (1 + j)·`surface_resistance_ohm(2π·th14.FREQUENCIES_HZ["64"], σ_json)`.
  - Z_doc = (1 + j)/(σ_json·δ), with δ = √(2/(ω_json·μ₀·μ_r,json·σ_json)) and ω_json = 2π·Freq_json·1e9.
  - Also asserted exactly: σ_json == `COPPER_SIGMA` (or `SIGMA_LADDER[-1]`), μ_r,json == 1.0, and 2π·Freq_json·1e9 equal to 2π·`FREQUENCY_64_HZ` within the same band.
  - Printed: δ, R_s, and \|Z_s\| beside ours' logged `2.087161e-03(1+j)` (`…item193-W2-solve.log:2911`). Also printed: the size of the CODATA-μ₀ variant (2.7e-10 relative), which Palace might use internally. It sits under the band, so a μ₀ convention cannot trip (Z).
- **Palace log check (executor, journalled), W3 / W5:**
  - 80 798 elements and the bounding box ±0.12 / ±0.12 / ±0.10 m;
  - `Rs = 4.375e+01 Ω/sq` on attributes 21–24;
  - `Configuring Dirichlet PEC BC at attributes:` listing **1 only**;
  - a Robin / conductivity configuration line listing **401**, quoted verbatim, with any σ or δ Palace prints compared against W1's printed values to the printed digits;
  - **no** "no associated boundary condition … PMC assumed" warning;
  - `ND (p = 1): 111946` (W3 / W4) or `ND (p = 2): 570332` (W5), as in item 234;
  - `ω/2π = 6.400e-02 GHz`;
  - four `GMRES solver converged` lines.
- **(P), W6:** item 159's set via `_p_check`, on `h05x_cu` (required) and on `h05x_cu_o2` if it ran: reciprocity ≤ `RECIPROCITY_BAND`; σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` and every column power ≤ 1 + 1e-9; every S-class spread ≤ `ADJACENT_SPREAD_BAND`; pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst. These are the bands item 234 imported.

**Ours against Palace (printed, never gated — Target: Palace, as F22).** `STEP3A_RECORD_RTOL` is a re-recording band for our own code, not a cross-solver one. So no band is imported for these readings.
- **(T3)** `ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]` vs Palace `h05x_cu`, \|Δ mean\|/\|mean\| per class. **Predicted ≤ 1e-6**, if Palace's conductivity boundary is the same discrete term. Note the ceiling: if the two boundary models differed entirely, (T3) could not exceed ours' own copper-vs-PEC move (below), which is ≤ 6.4e-4. So (T3) alone is weak evidence, and (T3m) is the discriminating reading.
- **(T3m) The copper move, as each code sees it.** Per class, Δ_Pal = mean(`h05x_cu`) − mean(`h05x_pec`) and Δ_ours = P1 `cu` − P1 `pec` (complex). Print \|Δ_Pal − Δ_ours\| / \|Δ_ours\| and arg(Δ_Pal/Δ_ours). **Predicted ≤ 1e-3 and \|arg\| < 0.1°.**
  - Ours' moves are 6.420e-4 / 2.811e-4 / 4.861e-4 relative (computed from `:116–121`).
  - The floor is set by the 10-digit record rounding and item 234's (T1) PEC tie (≤ 1.3e-10, `20261005T022223Z_OPS-61-item234-W5-read.log:48`). Both are about 1e-6 of a move this size.
  - A conjugated or 1/σ-scaled Palace impedance would show as a phase or magnitude miss here.
- **(T4)** If W5 runs: (T3) and (T3m) at order 2, against `P2["64"]["cu"]` / `["pec"]` and `h05x_cu_o2` / `h05x_pec_o2`. Ours' order-2 moves are 6.330e-4 / 2.866e-4 / 4.947e-4.
- **(C2⁗) Order-mismatch separation — predicted ≥ 25×, printed** (item 229's bar). Per class: \|ours P2 cu − Palace `h05x_cu`\| / \|ours P1 cu − Palace `h05x_cu`\|. The numerator is predicted to be ours' copper p1-to-p2 move, 4.35 / 1.67 / 3.11 % (computed from the records).

**Negative controls (rule (e)).**
- **(K1) The σ boundary reaches Palace's operator — floor asserted, size predicted.**
  - Asserted: every class mean of `h05x_cu` moves from `h05x_pec` by > `KNOB_REACHES_THE_SOLVE` (`tests/validation/test_port21_feed_sensitivity.py:156`, imported).
  - Backing: the same comparison on the same mesh, measured by ours. The class-mean moves are 6.42e-4 / 2.81e-4 / 4.86e-4 (records `:116–121`), and the entry-level max\|S_σ − S_PEC\| is 4.243e-4 / 1.285e-4 / 1.869e-4 (`…item193-W2-solve.log:2918`). The floor sits ≥ 280× below the smallest of them.
  - The size is (T3m), predicted only.
  - A Palace that silently ignored the block, or treated 401 as PEC, reads ≤ 1e-10 (item 234's (T1) floor).
- **(K3) Loss direction — asserted:** σ_max(`h05x_cu`) < σ_max(`h05x_pec`).
  - Backing: ours on the same mesh at 64 MHz, degree 1. Copper is 0.999851545100 (`…item193-W2-solve.log:2917`) and PEC is 0.999904317474 (`…item234-W5-read.log:65`), a gap of 5.3e-5. Ours' σ ladder rises monotonically toward PEC: 0.999868053 at 5.8e9 and 0.999900690 at 5.8e11 (`:2924, :2931`).
  - A wall with a sign error (active rather than lossy) would put σ_max above PEC's.
  - Predicted, printed: Palace's copper σ_max equal to ours' 0.999851545 to ≈ 1e-9.
- **(K2) σ scaling — predicted, printed.** `h05x_cu11` (σ = 5.8e11) vs `h05x_pec`, per class move, and the ratio (K1) / (K2).
  - Predicted ≈ 100 (Z_s ∝ σ^(−1/2); the Leontovich perturbation is linear at this small Z_s).
  - Comparand: ours' measured entry-level ratio on this fixture, 4.243e-4 / 4.250e-6 = 99.83, with 99.83 and 99.79 for the other two classes (`…item193-W2-solve.log:2918, :2932`).
  - A Palace σ that never reached the operator reads ≈ 1, and a Z ∝ 1/σ reads ≈ 1e4.
  - It is not asserted: Palace has no prior measurement on this comparison. The (K2) run's own (P) is printed only.
- (P)'s power to discriminate was shown by item 159's flip (200.0204 % spreads). It is not re-run.

**Predictions (written first, printed).**
- **W1:** < 5 s. (Z) holds to ≈ 1e-16 at μ₀ = 4π·10⁻⁷.
- **W3 (o1, copper):** Palace total 13–30 s (`h05x_pec` 13.408 s, `h05x_pmc` ≈ 14 s; `…item234-W2-palace-hole-pec.log:271`); peak 1.4–1.9 G. GMRES iterations are printed, not predicted (PEC took 17–19).
- **W4 (o1, σ = 5.8e11):** the same price.
- **(P):** reciprocity ≈ 1e-11. The spreads are predicted to equal ours' copper spreads to the printed digit: 0.0657 / 0.0612 / 0.0524 % (`…item193-W2-solve.log:2917`).
- **W5 (o2):** item 234's measured o2 / o1 ratio was 160.575 / 13.408 = 12.0 (`…item234-W4-palace-hole-pec-o2.log:324`), which projects 12.0 × T₁, ≈ 160–360 s; peak 4.5–6 G.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`).
- **Pre-check** (not a harness call): `ls logs/ops61/ logs/ops61/postpro_h05x_pec/`.
  - Stop and report if `h05x_pec.json`, the hole mesh or `postpro_h05x_pec/port-S.csv` is missing. Never re-run `export --hole-pec`: it asserts its outputs are absent.
  - Stop and report if any `h05x_cu*.json` or `postpro_h05x_cu*` exists.
- **W1** (< 10 s), the configs and (Z):

  `scripts/testing/run_and_log.sh OPS-61-item241-W1-cu-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py cu-configs'"`
- **W2** (smoke, ≤ 30 s), the keyword check. Palace's launcher options, the bundled examples' `Conductivity` blocks, and a configuration-only parse of `h05x_cu.json`:

  `scripts/testing/run_and_log.sh OPS-61-item241-W2-palace-keyword "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && palace --help; grep -rn -B2 -A6 -F Conductivity /opt/palace/examples; timeout -k 10 25 palace -dry-run h05x_cu.json'"`

  Read it before W3:
  - **Confirmed:** an example shows a `Conductivity` array under `Boundaries`, or the dry run exits 0. Journal the log line.
  - **Rejected:** the dry run aborts naming the `Conductivity` boundary keyword (for example, an unsupported-keyword message). Stop BLOCKED with the message verbatim (negative result below).
  - **Neither:** the launcher has no `-dry-run`, and no example shows the block. W3 is the validation: a configuration-time abort in W3 is the same BLOCKED reading.
- **W3** (≈ 13–30 s), copper, order 1:

  `scripts/testing/run_and_log.sh OPS-61-item241-W3-palace-hole-cu "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 h05x_cu.json'"`
- **W4** (≈ 13–30 s), the (K2) control: W3's command with `h05x_cu11.json`, labelled `OPS-61-item241-W4-palace-hole-cu11`.
- **Price gate (binding, after W3).** Read W3's `Total` T₁ and the iteration count.
  - **Run W5 iff 12.0 × T₁ ≤ 450 s**, i.e. T₁ ≤ 37.5 s. Print T₁, the iterations and the projection in the journal.
  - Otherwise order 2 is "not measured — price". W6 / W7 still run.
- **W5**, order 2: W3's command with `h05x_cu_o2.json`, labelled `OPS-61-item241-W5-palace-hole-cu-o2`. After any kill, run `docker compose -f docker/docker-compose.yml ps` and look for a lingering `fem-em-palace` container before W6.
- **W6** (≈ 3 s), the new read:

  `scripts/testing/run_and_log.sh OPS-61-item241-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --hole-cu'"`
- **W7** (≈ 3 s), the dispatch regression: item 234's W5 command verbatim (`20261005T022223Z_OPS-61-item234-W5-read.log:12`), labelled `OPS-61-item241-W7-regress-read`. Its output section must `diff` empty against that log's. **No probe edit after W1 starts** (rule (i)).

**Traps.**
- **Never run `export`, `export --exact-strip` or `export --hole-pec`.** The first two rewrite meshes and configs without an existence check; the third asserts that item 234's outputs are absent and would fail. Do not re-run `flip-config`, `freq-configs` or `o2-configs` either.
- **Attribute 401 goes in exactly one block.** If it stays in `PEC`, the run reproduces `h05x_pec`, and (K1) fails at ≈ 1e-10. If it is dropped from both, it is item 234's PMC control, and the W3 log shows the fallback warning.
- **Units.** `Model.L0` = 1.0, so the boundary σ is in S/m like the phantom's 0.5. `Freq` is in GHz, and `surface_resistance_ohm` takes ω in rad/s, not Hz (the trap `core/layered_impedance.py:23` names). Never write `Thickness` / `External`.
- **The conductivity boundary is untried here** (2026-10-03 daily, `:93`). An abort is the finding: record the MFEM / Palace message verbatim. Never change `Direction`, attributes, `Order`, `MaxIts`, `Tol` or the keyword in-slot, and never fall back to `Impedance` with a guessed Rs / Ls. A `GMRES … did not converge` line is a finding, not a retry: that run's S is "not measured".
- Compare S classes only. Never read Palace's `z_matrix_ohm` or the `aed` key. The record widths differ (P1 is `-n 2`, P2 is `-n 8`), and both are κ off. Never re-register a record from Palace.
- Palace's `postpro_*`, the configs and the mesh are gitignored under `logs/`. The harness logs are tracked.
- Record row and journal: identity readings, relative gaps, and Palace's own numbers in full (version, config, mesh cited). No absolute `S₁₁` / `Z_in` claim, and no Palace figure set against an AED one in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.
- Status 124: the price is the finding. Never lengthen a timeout or raise `-np` in-slot.

**Scope.** Covered:
- Palace's `Boundaries.Conductivity` on the exported exact hole mesh's wall at 64 MHz, order 1, plus order 2 if the gate allows;
- the (Z) arithmetic tie between our Leontovich Z_s and the thick-conductor (1 + j)/(σδ);
- (P) on Palace's copper 4×4;
- ours vs Palace per class from the tracked copper records;
- the (K1) / (K3) / (K2) controls.

Not covered:
- 10 / 128 MHz;
- Palace's `Impedance` (Rs / Ls / Cs) keyword and any `Thickness` model;
- the loss partition (P_coil / P_in; Palace's surface-loss postprocessing);
- the tuned `ANS-7` twin (needs Palace's lumped RLC);
- the F-human rung;
- any h-statement (one mesh);
- the private (a) / (b) / (c) reading;
- `src/` and any band.

**Status it can move.**
- The `OPS-61` row (`PROJECT_PLAN.md:1278`) gains a "Step 2b record (F23, §9 item 241)" sentence: Palace's first conductivity-boundary run, the (P) readings, and the matched-order (T3) / (T3m) / (T4) per class.
- **No glyph: 🟡 stays.** The row's Done-when names step 1c only (Palace `f05x`). Step 2 is "at closure", and the (a) / (b) / (c) reading is the weekly's (rule (j): read, no mismatch).
- **Rule (j)**, a proposed Done-when sentence for the review: "Step 2b is done when Palace `-np 8` solves the exported exact hole mesh at 64 MHz, order 1, with PEC on attribute 1 and `Boundaries.Conductivity` σ = `COPPER_SIGMA`, μ_r = 1 on 401 (80 798 elements, Rs 43.75 Ω/sq); the config's thick-conductor Z_s equals `(1 + j)·surface_resistance_ohm` within `STRIP_BAND`; its 4×4 meets item 159's (P) set at the imported bands; every class mean moves from the PEC column by more than `KNOB_REACHES_THE_SOLVE` with σ_max below the PEC column's; and ours-vs-Palace per class and per copper move is printed against `ANS6_EXACT_STRIP_CLASS_MEANS_P1['64']['cu']`."
- If (T3m) prints ≤ 1e-3 with (T3) ≤ 1e-6, the §2 / CLAUDE.md Palace clause ("matches ours to 0.0000 % at matched order") could extend to the copper hole route. That would be the first public external check of our Leontovich term at the S level. Both calls are the weekly's.
- §10 F23 completes.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(Z) red:** a config wrote the wrong σ, μ_r or Freq. Stop before any Palace window, revert the probe, BLOCKED (rule (d)).
- **The keyword is rejected (W2 or W3 abort at configuration):** this is the capability finding. Record the message verbatim, and open a known-issues entry ("Palace v0.18.1 `Boundaries.Conductivity` on the exported hole mesh: <message>"). Park the probe on `attempt/*` and mark the item BLOCKED. Unblock: the review re-poses with a keyword confirmed from Palace's own schema or an example file, or with an `Impedance` route whose series / parallel convention is confirmed first.
- **The W3 log lacks 401 in a non-PEC BC line, or shows the PMC fallback warning:** the wall is not carrying the conductor. Park, BLOCKED.
- **(P) red:** report every identity's number, park, BLOCKED. Never widen. Ours' copper column is green on this mesh (`…item193-W2-solve.log:2917`), so this would be the first identity-level disagreement on the copper route.
- **(K1) red** (copper moves S by ≤ 1e-6 from PEC): the σ block is not reaching the operator. Park, BLOCKED.
- **(K3) red** (σ_max(cu) ≥ σ_max(pec)): the wall adds energy, which is a sign-convention finding. Report it, park, BLOCKED; the review rules on it.
- **(T3) / (T3m) far from prediction, with (P), (K1) and (K3) green:** a reading, not a failure. The item lands DONE as a record, and the review reads it. PEC agrees to 1e-10 on this mesh, so a miss isolates the conductor-boundary model. (K2)'s ratio then says whether the σ-dependence differs, or only a constant factor or phase does.
- **W4 aborts or fails to converge:** (K2) is "not measured". (K1) / (K3) still carry the controls; journal it and continue.
- **Price gate trips:** a partial record (order 1 only). Not BLOCKED; the review re-queues order 2.
- **W7 not diff-empty:** the dispatch edit moved an old mode. Revert the probe, BLOCKED.

**Outcome (2026-10-05 04:30 CDT slot): ✅ DONE — step 2b recorded, no glyph; F23 complete.** The probe gained `_cu_set()`, `cu-configs` and `read --hole-cu` (additive; last edited 09:31:54 UTC, before W1 at 09:32:01). Pre-check: `h05x_pec.json`, the hole mesh and `postpro_h05x_pec/port-S.csv` present; no `h05x_cu*`. Palace v0.18.1, configs `logs/ops61/h05x_cu{,_o2,11}.json`, mesh `ops61_ans6_x1_hole_strip.msh`. Windows, all Status 0:
- W1 `20261005T093201Z_OPS-61-item241-W1-cu-configs.log` (3 s): (Z) \|Z_ours − Z_doc\|/\|Z_ours\| = 0.000e+00 on all three configs; δ 8.260687e-06 m, R_s 2.087161e-03 Ω (2.087161e-05 Ω at 5.8e11), equal to ours' logged value; CODATA-μ₀ variant 2.722e-10 (printed).
- W2 `20261005T093211Z_OPS-61-item241-W2-palace-keyword.log` (2 s): **keyword confirmed** — `palace --dry-run` exits 0, "Dry-run: No errors detected in configuration file "h05x_cu.json"" (`:303`); resolved block `Conductivity` / `External: false` / `Permeability: 1.0` / `Thickness: 0.0` (`:86–95`). The bundled examples carry `Conductivity` only as a material key (coaxial), not as a boundary block.
- W3 `20261005T093223Z_OPS-61-item241-W3-palace-hole-cu.log` (17 s; Palace Total 15.532 s, peak 1.4 G; GMRES 22 / 20 / 20 / 20): 80 798 elements, box ±0.12 / ±0.12 / ±0.10 m, `Configuring Robin finite conductivity BC at attributes:` / ` 401: σ = 5.800e+07 S/m, n = (+0.9,-0.5,+0.2)` (`:78–79`), Rs 4.375e+01 Ω/sq on 21–24, Dirichlet PEC on `1` only, no PMC-fallback warning, `ND (p = 1): 111946`, ω/2π 6.400e-02 GHz, four converged lines. Palace prints no δ.
- W4 `20261005T093251Z_OPS-61-item241-W4-palace-hole-cu11.log` (17 s; Total 15.478 s; GMRES 20 / 19 / 19 / 19): same checks, `401: σ = 5.800e+11 S/m`.
- Price gate: T₁ 15.532 s × 12.0 = 186 s ≤ 450 s → W5 ran.
- W5 `20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log` (177 s; Total 175.950 s, peak 4.6 G; GMRES 31 / 30 / 30 / 30): same checks, `ND (p = 2): 570332`. o2 / o1 = 11.3.
- W6 `20261005T093623Z_OPS-61-item241-W6-read.log` (2 s): (P) PASS on `h05x_cu` (reciprocity 1.327e-11, σ_max 0.999851545101, spreads 0.0714 / 0.0730 / 0.0349 %, separation 559.57×) and `h05x_cu_o2` (2.062e-11, 0.999901278001, 0.0175 / 0.0214 / 0.0044 %, 1956.25×). (T3) 1.276e-10 / 6.423e-11 / 2.524e-10; (T3m) 7.289e-09 / 4.126e-07 / 4.390e-07 with arg ±0.0000°; (T4) 1.612e-10 / 1.586e-10 / 2.033e-10; (T4m) 1.031e-07 / 2.263e-07 / 1.860e-07, arg ±0.0000°; (C2⁗) numerator 4.3517 / 1.6744 / 3.1058 % = ours p1→p2 move, separation ≥ 1.2e8×. (K1) moves 6.420e-04 / 2.811e-04 / 4.861e-04 > 1e-6: PASS. (K3) σ_max cu 0.999851545101 < pec 0.999904317452 (gap 5.277e-05): PASS. (K2) cu11 moves 6.431e-06 / 2.816e-06 / 4.871e-06, ratio 99.83 / 99.85 / 99.79, entry-level 99.82 (ours 99.83); cu11 σ_max 0.999900690069 (ours 0.999900690). Gitignored source, printed only: ours p1 cu 4×4 vs Palace entry max\|ΔS\| 8.933e-11.
- W7 `20261005T093640Z_OPS-61-item241-W7-regress-read.log` (2 s): output section diff-empty against `20261005T022223Z_OPS-61-item234-W5-read.log` (35 lines each).
- Prediction misses (printed, no effect): the (P) spreads read 0.0714 / 0.0730 / 0.0349 %, not the predicted 0.0657 / 0.0612 / 0.0524 % from ours' log `:2917`. Ours' and Palace's 4×4s agree entrywise to 8.9e-11, so the miss is the spread definition (th14's harness vs `_class_spread`), not a solver difference — not diagnosed further. The o2 price was 176 s against the ≈ 186 s projection.

## §7 row moved verbatim 2026-10-05 (rotate_plan_archive.py chunks)

 Palace as a second independent reference | 🟡 Steps 0, 1a and 1b are done, and the first half of step 1c (item 159, 2026-10-01) shows ours at degree 1 agreeing with Palace on the exact f = 0.5 strip mesh to 0.0000 % in every class, so the step-1a gaps were port width. The second half (item 166, 2026-10-02) completes step 1c as a record: ours at degree 2 on the same mesh is (F)/(G)/(K) green κ off and in-run, and sits 6.667 / 2.544 / 4.522 % from Palace order 1 (self / adjacent / opposite) — an order mismatch, numerically equal to ours p2 vs ours p1 (`20261002T095914Z_OPS-61-F10b-ii-W1-ours-p2.log`). Open: the (a)/(b)/(c) reading against the private column is the weekly's. Step 1d record (F17, §9 item 194, 2026-10-03): Palace on the same strip mesh at 10 / 128 MHz order 1 and at 64 MHz order 2 (the first Palace order-2 run, 268 s / 6.1 G at `-np 8`) passes its own reciprocity / passivity / class-spread identities at our imported bands; ours matches Palace to 0.0000 % per class at matched order in all four cases (degree 1 at 10 / 64 / 128 MHz, degree 2 at 64 MHz), so the order-1 → order-2 shift of 6.667 / 2.544 / 4.522 % is shared discretisation error, not a port-model difference (`20261003T110728Z_OPS-61-item194-W5-read.log`). A record, no glyph. Step 1e record (F21, §9 item 229, 2026-10-04): Palace order 2 at 10 / 128 MHz (211 s / 349 s, 6.0 / 6.3 G at `-np 8`) passes the same identities, and our tracked degree-2 record matches it to 0.0000 % per class (≤ 2.5e-10), so ours matches Palace to 0.0000 % at matched order in every class at 10 / 64 / 128 MHz, both orders (`20261004T214342Z_OPS-61-item229-W4-read.log`). A record, no glyph. Step 2a record (F22, §9 item 234, 2026-10-04): Palace's first hole-mesh run — `ANS-6`'s exact-rectangle hole mesh (80 798 tets, PEC on box 1 + wall 401, Rs 43.75 Ω/sq) at 64 MHz, order 1 (13.4 s, 1.4 G) and order 2 (160.6 s, 4.5 G) at `-np 8` — passes its own reciprocity / passivity / class-spread identities at our imported bands; our tracked PEC records match it to 0.0000 % per class at both matched orders (≤ 2.2e-10), and the PMC-wall control moves every class by ≥ 100 % (`20261005T022223Z_OPS-61-item234-W5-read.log`). A record, no glyph; the copper (impedance-boundary) half is the weekly's to enumerate. Step 2b record (F23, §9 item 241, 2026-10-05): Palace's first conductivity-boundary run — the same hole mesh with PEC on box 1 only and `Boundaries.Conductivity` σ = 5.8e7 S/m, μ_r 1, thick (resolved Thickness 0.0) on wall 401, whose thick-conductor Z_s equals our `(1 + j)·surface_resistance_ohm` to 0.0 (W1) — at 64 MHz, order 1 (15.5 s, 1.4 G) and order 2 (176.0 s, 4.5 G) at `-np 8` passes the same (P) identities at our imported bands; our tracked copper records match it per class to ≤ 2.6e-10 at both matched orders (T3 / T4), the copper-vs-PEC move as each code sees it agrees to ≤ 4.4e-7 relative with 0.0000° phase (T3m), the move clears `KNOB_REACHES_THE_SOLVE` with σ_max below the PEC column's (K1 / K3), and σ = 5.8e11 scales the move by 99.79–99.85 (ours 99.83) (`20261005T093623Z_OPS-61-item241-W6-read.log`). A record, no glyph — the first public external check of our Leontovich term at the S level, on one mesh at one frequency. **Done-when:** Palace `-np 8` solves `f05x` at 64 MHz: 117 206 elements, Rs 43.75 Ω/sq, port i = N_i, `+Z`, 50 Ω. Its 4×4 meets `RECIPROCITY_BAND`, σ_max and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every S-class spread ≤ `ADJACENT_SPREAD_BAND` with pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst — each asserted with the band imported. *History: `docs/planning/chunks/OPS-61.md`.* 

### §9 item 418 — 2026-10-13 (§10 chain F34; lane 1 → heavy)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* Licence: §10 chain F, F34 (PROJECT_PLAN.md, grep "**F34 (lane 1 → heavy)**"), quoted, not changed: "Palace at matched order on the h = 0.0075 exact mesh (entry 27's), 10 MHz first: `OPS-61`'s configuration with the resolution knob; target: matched-order agreement with our 0.0075 rung at `OPS-61`'s band (the 0.0000 % precedent at ×1); the measured price (`OPS-61` W4 268 s at ×1; predicted ≈ 800 s at 2.4× cells) decides whether 64 / 128 MHz run heavy or as an `xl` cost probe the 10-14 weekly writes. Serial on nothing." The Palace directive (§5.4 (i)–(v), PROJECT_PLAN.md grep "Palace is the public external reference, long term") and the sheets-and-κ ruling (exact rectangles, κ off) bind. **No band, tolerance, Done-when or validation target is set or moved.**

**`OPS-61` step 1f: Palace order 2 on `ANS-4`'s exact-strip h = 0.0075 rung at 10 MHz, against our degree-2 record of that rung.** Our side is **not re-solved**. It is read from `xl` entry 27's tracked log.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Export and read run at `-n 1` in the complex container; Palace runs at `-np 8` (the only width Palace has been measured at here). **Class:** chain step (solve), ≈ 12 measured min (18:00 restock table, PROJECT_PLAN.md grep "| chain solve |").
- **Time:** ≈ 55 predicted min: 15 fixed + ≈ 12 code + W1 ≈ 2.5 + W2–W5 ≈ 4 × 4 (range 3.5–7 each) + W6 / W7 < 1 + ≈ 5 spare.
- **Depends on:** nothing. Inputs on disk (checked 2026-10-13) or tracked:
  - `logs/ops61/f05x_o2_f10.json` (the ×1 order-2 10 MHz config, item 229) and `logs/ops61/postpro_f05x_o2_f10/port-S.csv`;
  - ours at h = 0.0075, degree 2, κ off, 10 MHz: the 10-digit S lines `docs/testing/logs/20261005T070009Z_ANS-4-step6b.log:3702–3706` (header `:3702`, `S_1k` … `S_4k` at `:3703–3706`), σ_max 0.999995526884 at `:3707`, cells 282 198 and 1 822 094 degree-2 unknowns at `:2167` / `:3700`;
  - ours at ×1, degree 2: `ANS4_EXACT_STRIP_S_P2_X1["10"]` (`tests/validation/test_ans4_exact_strip_records.py:233`).

**Why the window is split four ways (the price, and the binding ceiling).**
- The weekly's ≈ 800 s is for one Palace command with all four excitations. That is inside the heavy tier's 1200 s with 1.5× margin. **The binding ceiling is not 1200 s, though.** A slot's harness call is a foreground Bash call capped at 660 s. PROJECT_PLAN.md:952 says so, and §5.1's standing caution at PROJECT_PLAN.md:2231 is that ≈ 480 s is the safe container ceiling for a foreground slot. So no single slot command can hold an 800 s Palace run.
- **Item 229's breakdown** of the ×1 order-2 10 MHz run (`20261004T213356Z_OPS-61-item229-W2-palace-o2-f10.log:296–308`, 753 062 unknowns, four excitations at 26 GMRES iterations each, Total 211.4 s):
  - Setup 46.7 s, paid once;
  - Linear Solve + Preconditioner + Coarse Solve 126.8 s, i.e. ≈ 31.7 s per excitation;
  - Estimation 33.2 s, i.e. ≈ 8.3 s per excitation.
- **Scaled to the rung** (1 822 094 / 753 062 = 2.42× unknowns):
  - Setup ≈ 113–192 s (2.42^1.0 to 2.42^1.6; it includes the coarse factorisation);
  - ≈ 77–97 s per excitation (linear in unknowns, with 0–25 % more iterations on the finer mesh);
  - ≈ 20 s estimation per excitation.
  - **All four excitations in one command: ≈ 500–700 s.** That is over the slot ceiling.
  - **One excitation per command: ≈ 215–400 s, central ≈ 245 s.** Peak memory ≈ 14–22 G against the service's 128 G.
- **So the window is split.** W1 is the mesh export (our container, ≈ 90–150 s). W2–W5 are four Palace commands, each exciting one port with the other three as passive 50 Ω loads. W6 is the read. No command exceeds 590 s, and the item's whole compute is ≈ 17–29 min. Palace's factorisation is repeated four times. That costs ≈ 3 × Setup, the price of fitting the slot. It is disclosed, not hidden.

**Code (probe only; additive; disclosed).**
- **Constants**, beside `HOLE_MESH_NAME` (`:705`):
  - `H0075_RES = 0.0075`;
  - `H0075_MESH_NAME = "ops61_ans4_h0075_strip.msh"`;
  - `H0075_CELL_RECORD = 282198`, a literal with provenance `20261005T070009Z_ANS-4-step6b.log:2167, :3700`. No tracked constant holds this count; entry 27's (F′) asserted only "> 1.02 × 117 206".
  - `H0075_SET = {f"f05x_h0075_o2_f10_e{k}": k for k in (1, 2, 3, 4)}`.
- **New `export_exact_strip_h0075()`, dispatched by `export --exact-strip-h0075`.** A new function. **`export()` is neither edited nor lifted.** Its only regression run is `export --exact-strip`, which rewrites F10a's ×1 strip mesh and `f05x.json` without an existence check.
  1. Assert `comm.size == 1`. Assert that `H0075_MESH_NAME` and the four `H0075_SET` JSONs do not exist (never overwrite).
  2. Load `f05x_o2_f10.json` and check it the way `o2_configs` checks its template: `Model.Mesh == STRIP_MESH_NAME`, four ports with Attributes `[21]…[24]` all `+Z` and `R` 50, `Solver.Order == 2`, `Freq == [0.010]`.
  3. Build `_build(np.zeros(LEG_COUNT), resolution=H0075_RES, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` (`tests/mesh/test_birdcage_leg_offset.py:107–113`). This is the build `test_ans4_exact_strip_records.exact_rungs` passes under `FEM_EM_ANS4_EXACT_RESOLUTION=0.0075`. **Assert ncells == `H0075_CELL_RECORD`** before anything is written.
  4. The sheet steps, the N_i strip selection, the triangle and tet blocks and the Gmsh 2.2 writer are **copied** from `export(exact_strip=True)` (`:143–383`), with a provenance comment (commit and line span). Never call the function itself.
  5. Write four configs from the template, changing only:
     - `Model.Mesh` → `H0075_MESH_NAME`;
     - `Excitation` → `true` on port k and `false` on the other three (R 50 Ω stays on all four);
     - `Problem.Output` → `postpro_<name>`.
  6. Then assert (E1)–(E4) below. As in step 1b, a failure still leaves the artefacts to inspect.
- **New `read_exact_strip_h0075()`, dispatched by `read --exact-strip-h0075`.**
  - It reads each present `postpro_f05x_h0075_o2_f10_e{k}/port-S.csv` with `_read_port_s` (`:400`). It **asserts the column it keeps is labelled `S[i][k]`** for the excited port k. If Palace labels a single excitation `S[i][1]`, the read keeps that column by position and prints the labelling as a reading.
  - It assembles the 4×4. A new `_p_check_s(name, s, tag)` copies `_p_check`'s (`:752`) assertions onto a matrix. **`_p_check` is not edited.**
  - It parses ours from the tracked log. Find the header `[ANS-4 step5] 10 MHz kappa-off S (10 digits):` and take the next four `S_ik =` lines. Strip any leading `PASSED` (pytest output is interleaved, `:3703`), and read the four complex tokens per line with the regex `[+-]\d\.\d+e[+-]\d+[+-]\d\.\d+e[+-]\d+j`. Open the log with `gzip` if only `….log.gz` exists, since rotation may compress it.
  - It never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** two new `sys.argv` branches beside `:1693–1714`, plus a docstring paragraph. Every existing mode keeps its behaviour, which W7 proves.

**Anchors (asserted; bands imported, never restated).**
- **Export, W1:**
  - **(E1)** tets == `H0075_CELL_RECORD` (282 198) exactly.
  - **(E2)** exterior (attribute 1) area == `EXTERIOR_AREA_M2` 0.3072 within `EXTERIOR_AREA_BAND` (probe `:89–90`).
  - **(E3)** for each port, N_i area == `SPEC_AREA_M2` and A/h == `SPEC_WIDTH_M`, within `STRIP_BAND` relative (`tests/validation/test_ans4_exact_strip_records.py:131–133`; `tests/validation/test_ops61_exact_strip_ports.py:58`).
  - **(E4)** N_i ⊂ F_i, as step 1b's (A3).
- **(Tie) the parse, W6:** the parsed 4×4's σ_max equals the log's printed 0.999995526884 (`:3707`) within 1e-9 absolute. This is a parse check, not a physics band. The 10-digit rounding moves each entry by ≤ 5e-11, so it moves σ_max by ≤ ‖ΔS‖_F ≈ 3e-10. A wrong block misses by ≥ 1.7e-4 (the 64 MHz block's σ_max is 0.999823496773, `:3713`). A miss means the wrong block was parsed: stop.
- **(P), W6, on Palace's assembled 4×4** (item 159's set, as `_p_check` applies it):
  - ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py:86`);
  - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every column power ≤ 1 + 1e-9 (`tests/validation/test_port_birdcage_four_port.py:129`);
  - every S-class spread ≤ `ADJACENT_SPREAD_BAND` (`tests/validation/test_port_birdcage_lumped_column.py:121`);
  - pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst (`test_port_birdcage_four_port.py:173`).
  - Reciprocity here is a real cross-window check: S_ij and S_ji come from two different Palace commands.
- **Palace log check (executor, journalled), W2–W5:**
  - 282 198 elements, and the bounding box ±0.12 / ±0.12 / ±0.10 m;
  - `Rs = 4.375e+01 Ω/sq` on attributes 21–24;
  - PEC on 1 only;
  - ω/2π = 1.000e-02 GHz;
  - exactly one excited port, k;
  - one `GMRES solver converged` line;
  - `ND (p = 2)` printed. It is predicted to equal our 1 822 094: on the hole mesh the two counts matched, 570 332. A mismatch is a reading, not a stop.

**Ours against Palace (printed, never gated — `OPS-61`'s precedent).**
- **No cross-solver band exists.** Items 159 / 194 / 229 / 234 / 241 printed the matched-order gap and predicted ≤ 1e-6. `STEP3A_RECORD_RTOL` (`tests/validation/test_ans4_resolution_ladder.py:198`) is a re-recording band for our own code, not a cross-solver one (item 234). "0.0000 %" in the §2 wording is the `.4f %` print, i.e. < 5e-7.
- **(T5)** Per class and per entry, \|Palace o2 − ours p2\|/\|ours p2\| on the rung. **Predicted ≤ 1e-6** in every class, with the 10-digit log rounding as the floor (≈ 1e-10). Printed as "x.xxxx %" and in e-notation.
- **(C-h)** Palace o2 at h = 0.0075 against Palace o2 at ×1 (`postpro_f05x_o2_f10`), per class, printed beside ours' h-move 0.2527 / 0.0447 / 0.1375 % (`20261005T070009Z_ANS-4-step6b.log:3726`). Predicted equal to the printed digit: both codes see the same h-move.

**Negative controls (rule (e)).**
- **(C2-h) Mesh-mismatch separation — predicted, printed.** Per class: \|ours ×1 p2 (`ANS4_EXACT_STRIP_S_P2_X1["10"]`) − Palace o2 h\| / \|ours h p2 − Palace o2 h\|. The numerator is predicted to be the h-move above, 2.5e-3 / 4.5e-4 / 1.4e-3. The denominator is (T5), ≤ 1e-6. So the predicted separation is ≥ 447× on the weakest (adjacent) class, against item 229's ≥ 25× bar. The ceiling is set by the denominator's rounding floor: ≈ 4.5e-4 / 2.5e-10 ≈ 1.8e6×.
- **(P)'s power to discriminate** was shown by item 159's flip: 200.0204 % spreads against the 0.5 % band, the sign-flip ceiling, 400×. It is not re-run.
- **A wrong rung reads as a mesh-move-sized gap.** If the export built ×1 by mistake, (E1) fails first. If it built some other rung, (T5) prints the h-move (≈ 2.5e-3), not ≤ 1e-6.

**Predictions (written first, printed).**
- **W1:** mesh ≈ 60–110 s at `-n 1` (×1 strip 26.4 s; GEO-29 measured 281 728 cells at 0.0075 without the strip, so the cost scales ≈ 2.4×). Window 90–150 s. Box area 0.3072 to ≈ 1e-16; N_i 5.6e-5 to ≤ 2e-15.
- **W2–W5, each:**
  - Palace Total ≈ 245 s (215–400);
  - GMRES 26–33 iterations (×1 26);
  - peak 14–22 G.
- **(P):** reciprocity ≈ 1e-10. Spreads ≲ 0.1 % (ours on the rung: Z-class spreads 0.0346 / 0.0166 / 0.0168 %, `:3707`). σ_max 0.999995527 ± 1e-9.
- **(T5)** ≤ 1e-6; **(C-h)** equal to ours' h-move to the printed digit.

**Price readout and the 64 / 128 MHz decision (pre-registered; printed by W6 from the four W2–W5 logs).**
- Let T_e be the Palace `Total` and S the `Setup` line of the excitation-k run (Elapsed Time Report, Avg column). Let n₁₀ be that run's iteration count.
- The projected single-excitation total at frequency f is S + (T_e − S) × r_f, with r₆₄ = 37 / 26 = 1.423 and r₁₂₈ = 55 / 26 = 2.115. These are the ×1 order-2 iteration ratios: item 194 W4 (37) and item 229 W2 / W3 (26 / 55).
- The projected four-excitation single-command total is S + 4 (T_e − S) r_f. It is printed for the 10-14 weekly's `xl` entry.
- **Rule:** if max_k projected T_e(128) ≤ 450 s (590 / 450 = 1.31× margin), then 64 / 128 MHz are **heavy**, as four single-excitation windows per frequency, and the next review writes that item. Otherwise they are an **`xl` cost probe** written by the 10-14 weekly. Print the rule's inputs and outcome.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`, one-shot `run --rm` containers as every `OPS-61` Palace window has used — no `up` / `down`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call): `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up, and `ls logs/ops61/`:
  - stop and report if `ops61_ans4_h0075_strip.msh`, any `f05x_h0075_*.json` or any `postpro_f05x_h0075_*` exists;
  - `f05x_o2_f10.json` and `postpro_f05x_o2_f10/port-S.csv` must be present.
- **W1** (≈ 90–150 s), the export:

  `scripts/testing/run_and_log.sh OPS-61-item418-W1-export-h0075 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 400 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py export --exact-strip-h0075'"`
- **W2** (≈ 245 s, ≤ 590), excitation 1, **the price probe**:

  `scripts/testing/run_and_log.sh OPS-61-item418-W2-palace-h0075-e1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_h0075_o2_f10_e1.json'"`
- **Price gate (binding, after W2).** Read W2's Palace `Total` T₁.
  - **Run W3–W5 iff T₁ ≤ 450 s.** Print T₁, Setup, the iterations and the projection in the journal.
  - Otherwise columns 2–4 are "not measured — price". W6 still runs on column 1: it prints (T5) per entry of column 1 only, and (P) is "not measured". The item lands as a partial record, not BLOCKED, and the review re-queues or prices it.
- **W3 / W4 / W5:** W2's command with `…_e2.json` / `…_e3.json` / `…_e4.json`, labelled `OPS-61-item418-W3-palace-h0075-e2` / `-W4-…-e3` / `-W5-…-e4`.
  - **Minute gate:** if W5 cannot start by minute 50, stop launching. The `postpro_*` directories are durable and gitignored, so the next slot runs the missing excitations with no edit and then W6. Commit what landed, marked "partial; columns k… pending" (rule (d)).
  - After any kill, run `docker compose -f docker/docker-compose.yml ps` and look for a lingering `fem-em-palace` run container before the next window.
- **W6** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item418-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-h0075'"`
- **W7** (≈ 3 s), the dispatch regression: item 241's W6 command verbatim (`20261005T093623Z_OPS-61-item241-W6-read.log:12`), labelled `OPS-61-item418-W7-regress-read`. Its output section must `diff` empty against that log's.

**Traps already paid for.**
- **Never run `export`, `export --exact-strip` or `export --hole-pec`.** The first two rewrite meshes and configs without an existence check; the third asserts that item 234's outputs are absent. Never re-run `flip-config`, `freq-configs`, `o2-configs` or `cu-configs`; each asserts its outputs are absent.
- **The first 0.0075 export runs at `-n 1`** (the function asserts it). A gmsh failure deadlocks wider runs.
- **A slot command is capped at 660 s.** Never put two excitations in one Palace command, and never raise `-np` (no `-np 12` Palace point exists) or lengthen a timeout in-slot. A Status 124 is the price finding: record it and run no further e-window.
- **Docker socket.** Every harness call is the whole Bash command, bare and repo-relative (`bash_guard.py`). An executor that loses the docker socket mid-item stops and reports.
- **Orphans.** After any killed window, check before anything else: for the solver container, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0; for Palace, `docker compose -f docker/docker-compose.yml ps` must show no lingering `fem-em-palace` container.
- An abort on a port set that is not aligned with the axes is a finding: never change `Direction`, `Order`, `MaxIts`, `Tol` or attributes in-slot. A `GMRES … did not converge` line is a finding, not a retry: that column is "not measured".
- **Parse trap:** line `:3703` starts with `PASSED`, and other `S_1k` blocks exist earlier in the log (`:2498` …). Anchor the parse on the 10 MHz κ-off header only. (Tie) catches a wrong block.
- Compare S only. Never read Palace's `z_matrix_ohm` or the `aed` key. Never re-register a record from Palace.
- Record row and journal carry Palace's numbers in full (version 0.18.1, the four configs, the mesh) beside ours and their gaps (§5.4 (i)–(iii)). **No Palace figure is ever set against an AED one in a tracked file.** Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace order 2 at 10 MHz on the exported h = 0.0075 exact-strip rung; its (P) set on a 4×4 assembled from single excitations; the matched-order gaps against entry 27's record; the measured price and the 64 / 128 MHz decision. Not covered:
- 64 / 128 MHz;
- Palace order 1 on the rung;
- the h = 0.005 rung;
- any h-statement for Palace;
- the copper / hole columns;
- the private reading;
- `src/` and any band.

**Status it can move.**
- The `OPS-61` row (✅ 2026-10-05) gains a "step 1f record (F34, §9 item 418)" sentence: Palace on a second mesh, (P), and (T5) per class. **No glyph change.**
- **Rule (j):** the row's Done-when is step 1c's and is met. A record sentence needs no amendment.
- If (T5) ≤ 1e-6 in every class, the §2 / CLAUDE.md Palace clause ("matches ours to 0.0000 % at matched order") may extend to "on two rungs (×1, h = 0.0075)". That is the weekly's call.
- §10 F34's price question is answered by the pre-registered rule above.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(E1) red:** the export is not entry 27's mesh. Stop before any Palace window, park on `attempt/*`, BLOCKED (rule (d)). Open a known-issues entry only if the count moved on an unchanged image; a re-record is `record-reconciler`'s, under (1*).
- **(E2)–(E4) red:** the strip did not reach the mesh at this resolution. Report spans and facet counts, park, BLOCKED.
- **Palace aborts, or a W2–W5 log lacks 282 198 elements / Rs 43.75 Ω/sq:** record verbatim, open a known-issues entry, BLOCKED.
- **(Tie) red:** a parse bug. Fix it once in-slot, before W6 is re-run; still red, BLOCKED.
- **(P) red:** report every identity's number, park, BLOCKED. Never widen.
- **(T5) ≫ 1e-6 with (P) and (Tie) green:** a reading, not a failure. The item lands DONE as a record, and the review reads it. It would be the first matched-order disagreement between the codes on a refined mesh.
- **W7 not diff-empty:** the dispatch edit moved an old mode. Revert the probe, BLOCKED.
- **Price gate trips:** a partial record. Not BLOCKED.

**Outcome (2026-10-13, 07:30 slot, DONE — record):** every anchor green first time, nothing killed. (E1) 282 198 tets exact; (E2) exterior 0.3072 m² to 5.6e-17; (E3) each N_i 44 triangles, area / A/h to ≤ 8.9e-16 relative; (E4) N_i ⊂ F_i. W2–W5 Palace logs: 282 198 elements, Rs 4.375e+01 Ω/sq on 21–24, PEC on 1, 1.000e-02 GHz, one excited port each, one GMRES converged line each (25 iterations), `ND (p = 2): 1822094` = ours; Total / Setup (Avg) 146.2 / 34.2, 143.5 / 33.6, 143.4 / 33.5, 144.5 / 33.9 s, peak ≈ 14.4 G — under the predicted 215–400 s, price gate passed (T₁ 146.2 ≤ 450). Reading, not a stop: Palace prints attribute 22's normal as (−1, 0, 0) on this mesh (+1 on ×1); S unaffected. W6: (Tie) σ_max 0.999995526892 vs printed 0.999995526884, |Δ| 8.4e-12; each excitation wrote its own `S[i][k]` column; (P) PASS — reciprocity 1.525e-11, σ_max 0.999995526884, spreads self 0.0374 / adjacent 0.0173 / opposite 0.0161 %, column powers 0.813217154401 / 0.813253013484 / 0.812988621828 / 0.813196344670, pooled 9.5104 % = 253.98× the worst; (T5) per class self 1.917e-10 / adjacent 1.861e-11 / opposite 7.373e-11 ("0.0000 %"), per entry 1.8e-11 … 2.7e-10; (C-h) Palace's own h-move 0.2527 / 0.0447 / 0.1375 % = ours to the printed digit; (C2-h) separation 1.315e7 / 2.403e7 / 1.862e7× (the adjacent numerator is 2.112e-4 absolute = the brief's 4.5e-4 relative × |S| ≈ 0.47). Price readout: projected single-excitation 189.8–193.7 s at 64 MHz, 265.9–271.3 s at 128 MHz; four-excitation single command 658.7–672.1 / 962.8–982.5 s; rule: max_k T_e(128) 271.3 ≤ 450 s ⇒ **64 / 128 MHz heavy, four single-excitation windows per frequency**. W7 output section diff-empty against item 241's W6. Disclosed gap: the module-docstring dispatch paragraph and `usage:` string were not added (rule (i) forbade the edit after W1 started); the new functions carry their own docstrings. Logs `20261013T123326Z_OPS-61-item418-W1-export-h0075.log`, `…123438Z_…-W2-palace-h0075-e1.log`, `…123719Z_…-W3-palace-h0075-e2.log`, `…123955Z_…-W4-palace-h0075-e3.log`, `…124228Z_…-W5-palace-h0075-e4.log`, `…124505Z_OPS-61-item418-W6-read.log`, `…124517Z_OPS-61-item418-W7-regress-read.log`.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 418 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; one Palace excitation per window, never two heavy windows at once. Probe only — no `src/`, no test module; never run `export`, `export --exact-strip` or `export --hole-pec`; never open `docs/private/` or `aed_results/`."

### §9 item 420 — 2026-10-13 (§10 chain T39; lane 1 → heavy)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* Licence: §10 chain T, T39 (PROJECT_PLAN.md, grep "**T39 (lane 1 → heavy)**"), quoted, not changed: "Palace twin of the ESR-tuned F-small at degree 2, one frequency: `OPS-61`'s hole-route configuration with the four gap sheets as unexcited Palace lumped elements (R = 0.2 Ω in series with C_in) and the feed as the lumped port, at 64.000 MHz; target: S at the tuned node against T26's at `OPS-61`'s matched-order band; the price is the readout. Serial on nothing."

The rulings that apply:
- the capacitor-ESR ruling (§9, PROJECT_PLAN.md grep "Binding ruling — capacitor ESR (operator, 2026-10-06"), (i)–(iii);
- the sheets-and-κ ruling (exact rectangles, κ off);
- the Palace directive (§5.4 (i)–(v)).

**No band, tolerance, Done-when or validation target is set or moved.**

**What "the four gap sheets as unexcited lumped elements" means against our T26.** Our T26 Route A (item 311; `tests/validation/test_ans7_exact_strip.py:2401–2417`) works like this:
- Port 1 (the feed) is the **bare 50 Ω port**.
- P2–P4 carry series-RC sheets `series_rlc_impedance(f, r_ohm=0.2, c_f=C_in(r))` (`src/fem_em_solver/ports/lumped.py:138`).
- The feed's own (r + C) branch is added **in parallel in post-processing**: Z_in = 1/(1/Z_bare + 1/(r + 1/(jωC))), as in ruling 311-Q5 and `_z_branch` (`:2223`).

So the twin is:
- Palace port 1 = `LumpedPort` R 50 Ω, excited;
- Palace ports 2–4 = `LumpedPort`, `Excitation: false`, carrying the gap impedance;
- **the comparand is the bare S11 at port 1**, our S11_A.

The feed branch is circuit post-processing on both sides. It enters T40, not this item. The mesh is item 234's `ops61_ans6_x1_hole_strip.msh` (80 798 tets). That is the mesh `ANS-7`'s `_exact_geo` builds: the same `_sheets_build(True, as_hole=True, leg_sheet_strip_fraction=…)` call (`test_ans7_exact_strip.py:409–425`), and `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` = 80 798 (`tests/validation/test_ans6_exact_strip.py:112`) equals T26's printed cells (`20261012T094112Z_ANS-7-item395-W0-esr64.log:2137`). The copper wall is item 241's `Boundaries.Conductivity`, which matched our Leontovich term to ≤ 2.6e-10.

**Does Palace's lumped element do series R–C? Not confirmed by any file in the repo.**
- **What is confirmed:** Palace v0.18.1's resolved `LumpedPort` block carries `R`, `L`, `C`, `Rs`, `Ls`, `Cs` (defaults 0.0) and `Excitation` (`logs/ops61/postpro_h05x_cu_o2/h05x_cu_o2_resolved.json`, gitignored; the tracked echo is `docs/testing/logs/20261005T093211Z_OPS-61-item241-W2-palace-keyword.log:103–151`). The keys exist.
- **What is not confirmed:** how Palace combines them. The image ships no schema (`docker/palace.Dockerfile` copies `/opt/palace` and `examples/` only). The drafter's reading of Palace's public documentation is that a lumped port's R, L, C form a **parallel** RLC (Y = 1/R + 1/(jωL) + jωC), as HFSS's Lumped RLC does (ESR ruling (iii)).
- **The brief does not rely on that reading.** It discriminates in-run (control (K-short) below). The exact single-frequency series → parallel map is ESR ruling (iii)'s:
  - Q = 1/(ωCr), R_p = r(1 + Q²), C_p = C·Q²/(1 + Q²), with C = C_in(r), r = 0.2 Ω, ω = 2π·64.000 MHz;
  - then 1/(1/R_p + jωC_p) = r + 1/(jωC) **exactly** at that ω.
- **Values** (drafter's arithmetic from `T26_C_IN_R_X1` = 6.359301781809668e-11 F, `test_ans7_exact_strip.py:2184`): Q = 195.52429426, R_p = 7 646.14993 Ω, C_p = 6.35913544e-11 F (C_p/C − 1 = −2.6157e-5). Palace's sheet values are R_p·w/l ≈ 6.690e+03 Ω/sq and C_p·l/w ≈ 7.268e-11 F/sq (w = 7 mm, l = 8 mm, as Rs = 50 · 7/8 = 43.75 Ω/sq on `f05x`). The identity's residual in double precision is 3.6e-16 relative.
- **Entering 0.2 Ω as Palace's R would short the gap** if the convention is parallel (ruling (iii)): \|Z\| 0.2 Ω against the series branch's 39.1 Ω. That misread is run on purpose as the control.

`OPS-61` step 3a: Palace order 2 on the exported hole mesh at 64.000 MHz, port 1 excited at 50 Ω, ports 2–4 carrying the ESR-tuned capacitors, against our T26 Route A S11.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Configs and the read run at `-n 1` in the complex container; Palace runs at `-np 8`. **Class:** chain step (solve), ≈ 12 measured.
- **Time:** ≈ 40 predicted min: 15 fixed + ≈ 12 code + W1 / W2 < 1 + W3–W5 ≈ 3 × 1.5 (1–2.2 each) + W6 / W7 < 1 + ≈ 5 spare.
- **Depends on:** nothing. Inputs on disk (checked 2026-10-13) or tracked:
  - `logs/ops61/ops61_ans6_x1_hole_strip.msh` and `h05x_cu_o2.json` (item 241);
  - `logs/ans7_exact/tune_p2.json` (`c_in` = C_in(0)) and `tune_p2_esr0.2.json` (`c_in_r`);
  - ours' S11_A at (r, C_in(r)), 64.0 MHz: −2.392625126044e-01+9.637828488809e-01j (`20261012T094112Z_ANS-7-item395-W0-esr64.log:2128`; 12 digits; the latest window that rewrote item 311's ESR files);
  - ours' lossless S11_A at C_in(0), 64.0 MHz: −2.404062978868e-01+9.685876774136e-01j (`20261011T093544Z_ANS-7-item311-W1-cin64-lossless.log:2055`; C_in(0) 6.359009101426113e-11 F at `:2058`).

**Code (probe only; additive; disclosed).**
- **New `ans7_configs()`, dispatched by `ans7-configs`.**
  1. Assert `comm.size == 1`.
  2. Assert none of `t39_esr.json`, `t39_lossless.json`, `t39_short.json` exists (never overwrite).
  3. Load `h05x_cu_o2.json` and check it as `cu_configs` checks its template:
     - `Model.Mesh == HOLE_MESH_NAME`;
     - four ports, Attributes `[21]…[24]`, all `+Z`, `R` 50;
     - `Order == 2`, `Freq == [0.064]`;
     - `PEC.Attributes == [1]`;
     - `Conductivity` σ == `COPPER_SIGMA` on `[401]`.
  4. Import `T26_C_IN_R_X1` from `tests.validation.test_ans7_exact_strip` (a module import, never a `test_*` name). Assert it equals `tune_p2_esr0.2.json["c_in_r"]` bit for bit (pure quantity). Assert `tune_p2.json["c_in"]` equals the literal 6.359009101426113e-11 (`…item311-W1-cin64-lossless.log:2058`) bit for bit.
  5. Write three configs from the template, changing only the keys named. In every config, port 1 keeps `R` 50 with `Excitation: true`, and ports 2–4 get `Excitation: false` and `R` removed unless stated. `Problem.Output` → `postpro_<name>`.
     - **`t39_esr`:** ports 2–4 `R` = R_p and `C` = C_p from the map at 64.000 MHz.
     - **`t39_lossless`:** ports 2–4 `C` = C_in(0), no `R`. This tie is convention-free: a pure C is the same in series or in parallel.
     - **`t39_short`:** ports 2–4 `R` = 0.2 and `C` = C_in(r). This is ruling (iii)'s misread, the control.
  6. Then assert (Z) below, reading `R`, `C` and `Freq` back from each written JSON. Print each config's `Boundaries.LumpedPort` block.
- **New `read_ans7_esr()`, dispatched by `read --ans7-esr`.**
  - It reads each present `postpro_t39_*/port-S.csv` with `_read_port_s` (`:400`) and keeps **S[1][1] only**. Palace's S[i][1] for i = 2–4 are referenced to impedances that are not ours' 50 Ω and are never compared.
  - It prints each run's resolved `LumpedPort` block, which is file evidence of the keys and defaults Palace filled in.
  - It prints the readings below.
  - It never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** two new `sys.argv` branches beside `:1693–1714`. Every existing mode keeps its behaviour, which W7 proves.

**Anchors (asserted; bands imported, never restated).**
- **(Z), W1: the map's own identity, pure numpy.** For `t39_esr`, \|Z_par − Z_ser\| / \|Z_ser\| ≤ `STRIP_BAND` (1e-9, `tests/validation/test_ops61_exact_strip_ports.py:58`), applied as item 241's (Z) applied it.
  - Z_par = 1/(1/R_json + jω_json C_json) and Z_ser = `series_rlc_impedance(64.0e6, r_ohm=0.2, c_f=T26_C_IN_R_X1)`, with ω_json = 2π·Freq_json·1e9.
  - Also asserted: 2π·Freq_json·1e9 equals 2π·`TUNE_TARGET_HZ` (`test_ans7_exact_strip.py:928`) within the same band.
  - Printed (not asserted): R_p, C_p, Q, the residual (predicted ≈ 4e-16), and the short config's \|Z\| 0.2000 Ω beside the series 39.105 Ω.
- **(F) the mesh:** each W3–W5 Palace log shows:
  - 80 798 elements, and the bounding box ±0.12 / ±0.12 / ±0.10 m;
  - `Configuring Robin finite conductivity BC` on 401 with σ = 5.800e+07, and Dirichlet PEC on 1 only;
  - `ND (p = 2): 570332`;
  - ω/2π = 6.400e-02 GHz;
  - one excitation, one `GMRES solver converged` line;
  - port 1 `Rs = 4.375e+01 Ω/sq`;
  - for ports 2–4, whatever surface R / C line Palace prints, quoted verbatim and compared with the predicted 6.690e+03 Ω/sq and 7.268e-11 F/sq (`t39_esr`), 1.750e-01 Ω/sq (`t39_short`) to the printed digits. **Executor-journalled**, as items 234 / 241's log checks.

**Ours against Palace (printed, never gated — `OPS-61`'s precedent).** No cross-solver band exists: items 159 / 194 / 229 / 234 / 241 all printed and predicted ≤ 1e-6. `STEP3A_RECORD_RTOL` is our own re-recording band (item 234), and the §2 "0.0000 %" is the `.4f %` print, < 5e-7. If the review wants T39's "at `OPS-61`'s matched-order band" asserted, it must name the band (Q1 below); this brief does not choose one.
- **(T6) the tuned node:** \|S11_Pal(`t39_esr`) − S11_A\| / \|S11_A\|, against `…item395-W0-esr64.log:2128`. **Predicted ≤ 1e-6** if Palace's lumped R ‖ C is parallel: same discrete system, map exact at 64.000 MHz, floor ≈ 1e-12 from the 12-digit print, GMRES `Tol` 1e-10.
- **(T6-0) the C path alone:** \|S11_Pal(`t39_lossless`) − S11_A,lossless\| / \|S11_A,lossless\| against `…item311-W1-cin64-lossless.log:2055`. **Predicted ≤ 1e-6**, whatever the R convention.
- Printed beside them: \|S11\| of each run; the Z_bare each implies through `_ANS7._z_from_s11`; Palace Total and iterations per run.

**Negative controls (rule (e); every factor predicted, since Palace has never run a lumped R / C here).**
- **(K-short), the convention discriminator and the misread control.** \|S11_Pal(`t39_short`) − S11_A\| / \|S11_A\|.
  - If the convention is parallel, it is **predicted ≥ 5e-2, in [5e-2, 2.0]**: three gaps at 0.2 Ω instead of 39 Ω rewrite the coil.
  - The ceiling is (1 + \|S11_A\|)/\|S11_A\| = 2.007 (\|S11_A\| = 0.99303).
  - **Exactly one of (T6) / (K-short) is predicted to read ≤ 1e-6**: (T6) if Palace is parallel, (K-short) if Palace is series. Printed as the "convention reading": parallel / series / neither.
- **(K-esr) r reaches the operator:** \|S11_Pal(`t39_esr`) − S11_Pal(`t39_lossless`)\| / \|S11_Pal(`t39_lossless`)\|, printed beside ours' same move. Ours' move comes from the two 12-digit S11_A lines above: ≈ 4.9e-3, the drafter's arithmetic. The executor recomputes it and prints it.
  - Predicted equal to ours to ≤ 1e-6 relative.
  - A Palace that dropped `R` (or a typo key) reads ≈ 1e-10 here.

**Predictions (written first, printed).**
- **W1:** < 5 s; (Z) ≈ 4e-16.
- **W3–W5, each:**
  - Palace Total ≈ 70 s (60–130). Item 241 W5's breakdown, hole mesh, order 2, four excitations, 175.95 s (`20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log:312–327`): Setup 30.6 s once, 29.7 s iterative and 5.8 s estimation per excitation, so one excitation ≈ 69 s.
  - GMRES 30–80 iterations. The tuned node may slow convergence; `MaxIts` is 100.
  - Peak 3–4.6 G.
- **(T6), (T6-0):** ≤ 1e-6. **(K-short):** in [5e-2, 2.0]. **Convention reading:** parallel.
- **The price readout for T40:** `t39_esr`'s Palace Total, predicted ≤ 120 s.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`, one-shot `run --rm` as every `OPS-61` Palace window). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call): `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up, and `ls logs/ops61/ logs/ans7_exact/`:
  - stop and report if any `t39_*.json` or `postpro_t39_*` exists;
  - `h05x_cu_o2.json`, the hole mesh, `tune_p2.json` and `tune_p2_esr0.2.json` must be present. Never re-create one.
- **W1** (< 10 s), configs and (Z):

  `scripts/testing/run_and_log.sh OPS-61-item420-W1-ans7-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py ans7-configs'"`
- **W2** (smoke, ≤ 30 s), the keyword check, item 241 W2's shape:

  `scripts/testing/run_and_log.sh OPS-61-item420-W2-palace-dryrun "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 10 25 palace -dry-run t39_esr.json'"`

  It must print "No errors detected". A configuration-time rejection of `R` / `C` on an unexcited `LumpedPort` is the capability finding (negative result below).
- **W3** (≈ 70 s), `t39_esr`:

  `scripts/testing/run_and_log.sh OPS-61-item420-W3-palace-t39-esr "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 t39_esr.json'"`
- **W4 / W5:** W3's command with `t39_lossless.json` / `t39_short.json`, labelled `OPS-61-item420-W4-palace-t39-lossless` / `OPS-61-item420-W5-palace-t39-short`. After any kill, run `docker compose -f docker/docker-compose.yml ps` and look for a lingering `fem-em-palace` container before the next window.
- **W6** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item420-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --ans7-esr'"`
- **W7** (≈ 3 s), the dispatch regression: item 241's W6 command verbatim (`20261005T093623Z_OPS-61-item241-W6-read.log:12`), labelled `OPS-61-item420-W7-regress-read`. Its output section must `diff` empty against that log's.

**Traps already paid for.**
- **Never enter 0.2 Ω as a parallel R** (ESR ruling (iii)), except in the labelled `t39_short` control.
- **The map is single-frequency.** R_p and C_p are exact at 64.000 MHz only. R_p follows 1/f² away from it (ruling (iii)), so never reuse `t39_esr`'s values at another frequency. That is T40's trap.
- **The feed's (r + C) branch is post-processing, not a Palace element.** Putting a capacitor on port 1 would make Palace's S11 a different quantity from ours' S11_A (ruling 311-Q5 applies in T40's read).
- **Compare S[1][1] only.** S[i][1] for i ≥ 2 are referenced to Palace's port impedances, not ours' 50 Ω. Never read Palace's `z_matrix_ohm` or the `aed` key.
- **Never run `export`, `export --exact-strip` or `export --hole-pec`**, and never re-run `flip-config`, `freq-configs`, `o2-configs` or `cu-configs`.
- **Lumped R / C in Palace is untried here.** An abort is the finding: record the MFEM / Palace message verbatim. Never change `Direction`, attributes, `Order`, `MaxIts`, `Tol` or a key name in-slot. A `GMRES … did not converge` line is a finding, not a retry: that run's S is "not measured".
- **Docker socket and orphans:** every harness call is the whole Bash command, bare and repo-relative, and an executor that loses the socket stops and reports. After any killed window: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0, and `docker compose -f docker/docker-compose.yml ps` must show no lingering `fem-em-palace`. Never lengthen a timeout or raise `-np`.
- Record row and journal: Palace's numbers in full (v0.18.1, the three configs, the mesh) beside ours, and the gaps. **No Palace figure set against an AED one in a tracked file**, and no AED ESR case exists to tempt it. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered:
- the bare S11 at the ESR-tuned node at 64.000 MHz, order 2, on the ×1 hole mesh, Cu;
- the lossless C-path tie;
- the convention reading and the misread control;
- the price.

Not covered:
- f₀, R_in, Q (T40);
- 127.74 MHz;
- the h = 0.01 rung;
- order 1;
- any AED reading;
- a measured ESR;
- `src/` and any band.

**Status it can move.**
- The `OPS-61` row (✅) gains a "step 3a record (T39, §9 item 420)" sentence: Palace's first lumped R / C run, the convention reading, and (T6) / (T6-0).
- The `ANS-7` row gains "T39 (item 420, record): Palace twin of the ESR-tuned node at matched order". It holds 🟡: its Done-when is Q's h-statement (T30 / F31 / T41), which this does not write.
- **Rule (j):** both rows read; no mismatch.
- T40 (item 421) becomes runnable or is returned to the weekly by its own first step.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(Z) red:** a config wrote the wrong R, C or Freq. Stop before any Palace window, revert the probe, BLOCKED (rule (d)).
- **W2 rejects the config, or W3 aborts at configuration:** the capability finding. Record verbatim, and open a known-issues entry ("Palace v0.18.1 `LumpedPort` R / C on an unexcited port: <message>"). Park the probe on `attempt/*` and mark BLOCKED. Unblock: the review re-poses on Palace's `Boundaries.Impedance` (Rs / Ls / Cs), with its series / parallel convention confirmed first.
- **(F) red** (wrong element count, PEC on 401, PMC-fallback warning): the export is not the solved mesh. Park, BLOCKED.
- **Convention reading "series"** ((K-short) ≤ 1e-6, (T6) ≫): a finding about Palace's convention, not a failure. The item lands DONE as a record and journals it. T40's config route switches to (r, C_in(r)) directly at each frequency (item 421's step 0 reads it).
- **Convention reading "neither"** with (T6-0) ≤ 1e-6: the C path agrees and the R path matches neither convention. That is a reading for the review. Land DONE as a record, and item 421 stops at its step 0.
- **(T6-0) ≫ 1e-6:** the C path or the port geometry disagrees on a convention-free quantity, the first lumped-element disagreement between the codes. Land DONE as a record with a known-issues entry; item 421 stops at its step 0.
- **W7 not diff-empty:** revert the probe, BLOCKED.

**Question for the review.** (Q1) T39's target says "at `OPS-61`'s matched-order band". No such constant exists, and every `OPS-61` precedent printed the gap. The brief prints (T6) / (T6-0) against the ≤ 1e-6 prediction and asserts only (Z). Registering a cross-solver band (for example 5e-7, the "0.0000 %" print threshold) is the review's or the weekly's call, not the slot's.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 420 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; never two heavy windows at once. Probe only — no `src/`, no test module; never enter 0.2 Ω as a parallel R outside `t39_short`; never open `docs/private/` or `aed_results/`."

**Outcome (2026-10-13, 09:00 slot, DONE — record; no glyph):** every anchor and control green first run, nothing killed. The probe gained `T39_*` constants, `_t39_map`, `ans7_configs()` (`ans7-configs`) and `read_ans7_esr()` (`read --ans7-esr`), additive; last edited before W1 (14:02:18 UTC). Pre-check: `fem-em-solver` Up; inputs present; no `t39_*` config or `postpro_t39_*`. Palace v0.18.1, configs `logs/ops61/t39_{esr,lossless,short}.json`, mesh `ops61_ans6_x1_hole_strip.msh`.
- **W1 (Z):** \|Z_par − Z_ser\|/\|Z_ser\| 3.634e-16 ≤ `STRIP_BAND`; \|dω\|/ω 1.482e-16; Q 195.52429426, R_p 7 646.14992955 Ω, C_p 6.3591354418e-11 F (C_p/C − 1 −2.6157e-5); `T26_C_IN_R_X1` = `tune_p2_esr0.2.json` `c_in_r` and `tune_p2.json` `c_in` = 6.359009101426113e-11 bit for bit; the short control's \|Z\| 0.2000 Ω vs series 39.105 Ω (printed).
- **W2:** `Dry-run: No errors detected in configuration file "t39_esr.json"`; the resolved block carries R_p / C_p on ports 2–4 with `Excitation: false`.
- **(F), W3–W5 (executor-journalled):** 80 798 elements; box (±1.200e-01, ±1.200e-01, ±1.000e-01) m; `401: σ = 5.800e+07 S/m`; Dirichlet PEC on `1` only; `ND (p = 2): 570332`; ω/2π = 6.400e-02 GHz; one excitation (index 1 on 21); port 1 `Rs = 4.375e+01 Ω/sq`. Ports 2–4 verbatim: `t39_esr` `Rs = 6.690e+03 Ω/sq, Cs = 7.268e-11 F/sq` (predicted 6.690e+03 / 7.268e-11, equal); `t39_lossless` `Cs = 7.267e-11 F/sq` (C_in(0)·8/7 = 7.2674e-11); `t39_short` `Rs = 1.750e-01 Ω/sq, Cs = 7.268e-11 F/sq` (predicted 1.750e-01, equal). Circuit lines `Index = 2: R = 7.646e+03 Ω, C = 6.359e-11 F` etc. GMRES 34 / 35 / 33 iterations; Palace Total 51.765 / 53.164 / 51.132 s (log Avg; `palace.json` 51.778 / 53.176 / 51.143 s); peak 4.2–4.3 G; harness 53 / 55 / 53 s.
- **W6 (printed, none asserted):** S11_Pal(`t39_esr`) −2.392625129e-01+9.637828488e-01j, \|S11\| 0.993037728250; (T6) **2.569e-10**; (T6-0) **2.565e-10**; (K-short) **1.437e-01** (S11 −9.9427e-02+9.9207e-01j, Z_bare +0.1348+45.24j Ω, in [5e-2, 2.0]); (K-esr) Palace 4.949105e-03 = ours 4.949105e-03, \|ratio − 1\| 1.994e-9; **convention reading: parallel** — Palace v0.18.1's `LumpedPort` R, C form R ‖ C, as HFSS's Lumped RLC (ESR ruling (iii)). Z_bare(`t39_esr`) +0.2815019608+39.10426459j Ω vs ours +0.2815019607+39.10426460j Ω.
- **Price readout for T40:** `t39_esr` Palace Total 51.8 s ≤ 120 s per frequency — item 421 runs heavy.
- **W7:** item 241's W6 command, output section diff-empty (41 = 41 lines, 0 diff lines).
- Logs `20261013T140218Z_OPS-61-item420-W1-ans7-configs.log`, `…140228Z_…-W2-palace-dryrun.log`, `…140240Z_…-W3-palace-t39-esr.log`, `…140343Z_…-W4-palace-t39-lossless.log`, `…140447Z_…-W5-palace-t39-short.log`, `…140547Z_…-W6-read.log`, `…140556Z_…-W7-regress-read.log`.
- Disclosed gap (as item 418): the module-docstring dispatch paragraph and `usage:` string were not extended (rule (i)); the new functions carry their own docstrings.

### §9 item 421 — 2026-10-13 (§10 chain T40; lane 1 → heavy; depends on item 420)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* Licence: §10 chain T, T40 (PROJECT_PLAN.md, grep "**T40. Palace's driven sweep**"), quoted, not changed: "Palace's driven sweep on T39's configuration — T26's five nodes; `f₀` and `R_in` / Q read as ours are; target: `f₀` within T26's zero resolution, Q within the 5 % the copper contribution is quoted at (the Q clause's public external check while F31 settles its h). Heavy if T39 reads ≤ 120 s per frequency, else the 10-14 weekly tiers it. Serial on T39." The ESR ruling (iii), the sheets-and-κ ruling and the Palace directive apply. **No band, tolerance, Done-when or validation target is set or moved.**

**Step 0 (first, before any code or window): the tier gate and the convention gate.**
1. Read item 420's W3 log (`…OPS-61-item420-W3-palace-t39-esr.log`), the Elapsed Time Report `Total`, Avg column: T₄₂₀.
2. **If T₄₂₀ > 120 s:** stop. No code, no window. Mark item 421 BLOCKED with "tier returned to the 10-14 weekly (T₄₂₀ = … s > 120 s per frequency; §10 T40)" (rule (d)), journal it, and take the next item.
3. **If item 420 did not land**, or its convention reading is "neither", or (T6-0) read ≫ 1e-6: stop the same way, with that reason.
4. **Otherwise the route follows item 420's convention reading.**
   - **Parallel:** each frequency gets its own map (R_p(f), C_p(f)) at C_in(r).
   - **Series:** each frequency gets (R = 0.2, C = C_in(r)) directly. Under that convention one Palace multi-point sweep would be exact, but the brief still runs one point per window, for the price.

**Why three frequencies, not five (a deviation from the §10 letter; the review rules — Q1 below).**
- T26 has no five-node set. Its Route A, the in-model quantity this twin reproduces, is read at the three tune points 63.9 / 64.0 / 64.1 MHz (`TUNE_MHZ`, `test_ans7_exact_strip.py:929`) through `parallel_rows` (`ans:7`, `07_tuned_copper_birdcage_driven_sweep.py:238–262`). Its Route B uses item 196's seven ×1 nodes, not five. Item 214's five nodes are the h = 0.01 rung's.
- `parallel_rows` reads f₀, R_in and Q from the bracketing pair alone (the − → + crossing interval and its linear interpolants). Points outside (63.9, 64.1) cannot change a row.
- Five nodes 1.3 MHz apart (item 214's spacing) would interpolate Im Y_in across ≈ 3 half-power bandwidths (f₀/Q ≈ 0.39 MHz), against a band of 1.287e-5 (≈ 824 Hz). That is unreadable by construction.
- So "read as ours are" means T26's triple. Item 420 supplies 64.000 MHz; this item adds 63.9 and 64.1.

`OPS-61` step 3b: Palace order 2 at 63.9 and 64.1 MHz on T39's configuration, and with item 420's 64.000 MHz run, Palace's f₀, R_in and Q at the ESR-tuned node read through `parallel_rows`, against T26's Route A rows.
- **Executor:** implementer. Probe only; no test module, no `src/`.
- **Tier:** heavy. Configs and the read run at `-n 1`; Palace runs at `-np 8`. **Class:** chain step (solve), ≈ 12 measured.
- **Time:** ≈ 30 predicted min: 15 fixed + ≈ 6 code + W1 < 1 + W2 / W3 ≈ 2 × 1.2–2 (from T₄₂₀) + W4 / W5 < 1 + ≈ 5 spare.
- **Depends on:** item 420, through step 0.

**Code (probe only; additive; disclosed).**
- **New `ans7_sweep_configs()`, dispatched by `ans7-sweep-configs`.**
  - It writes `t40_esr_f63.9.json` and `t40_esr_f64.1.json` from `t39_esr.json`, changing only:
    - `Freq` → [0.0639] / [0.0641];
    - ports 2–4 `R` / `C`, from the map at that frequency (parallel route) or (0.2, C_in(r)) (series route);
    - `Problem.Output`.
  - It asserts neither file exists.
  - It asserts (Z′) per config, reading the JSON back.
- **New `read_ans7_sweep()`, dispatched by `read --ans7-sweep`.**
  1. Take S11 at 63.9 / 64.0 / 64.1 from `postpro_t40_esr_f63.9`, `postpro_t39_esr`, `postpro_t40_esr_f64.1` (`_read_port_s`, S[1][1] only).
  2. Form Z_bare = `_ANS7._z_from_s11(S11)`, then Z_in = 1/(1/Z_bare + 1/`_z_branch`(ω, C_in(r), 0.2)). The feed carries r too (ruling 311-Q5). Import `_z_branch` and `parallel_rows` from `tests.validation.test_ans7_exact_strip`; never a `test_*` name.
  3. Read rows = `parallel_rows([63.9e6, 64.0e6, 64.1e6], Z_in)`.
  4. Print the per-point S11 ties against ours' S11_A at 63.9 / 64.1 (`20261012T094112Z_ANS-7-item395-W0-esr64.log:2127`, `:2129`). Each is predicted ≤ 1e-6.

**Anchors (asserted; bands imported or quoted from the target, never chosen).**
- **(Z′), W1:** per config, \|Z_json − Z_ser(f)\| / \|Z_ser(f)\| ≤ `STRIP_BAND` (`tests/validation/test_ops61_exact_strip_ports.py:58`), with Z_ser(f) = `series_rlc_impedance(f, r_ohm=0.2, c_f=T26_C_IN_R_X1)`, as item 420's (Z).
- **(C1) one crossing:** Palace's Im Y_in changes sign − → + exactly once on (63.9, 64.0, 64.1), and never + → −. f₀,Pal is inside (63.9, 64.1) MHz.
- **(f₀)** \|f₀,Pal / f₀,A − 1\| ≤ `T21P_ZERO_RESOLUTION` (1.287e-5, `tests/validation/test_ans7_exact_strip.py:1018`, imported). f₀,A = 63 999 834.680890724 Hz (`…item395-W0-esr64.log:2130`). This is T40's "within T26's zero resolution", the band T26 (T2-r) asserted.
- **(Q)** \|Q_Pal / Q_A − 1\| ≤ 0.05. Q_A = 164.5434331272106 (`:2130`).
  - **The 5 % is T40's target, quoted from §10, not a code constant.** It is the copper contribution's AGREE level as the 2026-09-30 weekly quoted it (PROJECT_PLAN.md §6 Phase 6 row, "AGREE on the mode and the copper contribution (5 %)").
  - The probe holds it as `T40_Q_TARGET = 0.05` with that provenance comment. It is not a band this item chose. The review confirms (Q2).

**Printed, never gated.**
- R_in,Pal against R_in,A = 3 175.9794661876795 Ω (`:2130`), the relative gap;
- f₀ and Q gaps in e-notation;
- r_AB-style f₀,Pal / f₀,A − 1;
- Palace's Im Y_in triple beside ours (−1.616350144e-04 / +2.676560526e-07 / +1.624710261e-04 S, `:2130`).

All are predicted ≤ 1e-6 relative.

**Negative controls (rule (e); predicted, printed — no prior measurement of these comparisons on Palace data).**
- **(N-feed) the feed's r omitted**, pure post-processing on the same three Palace S11. Z_in′ = 1/(1/Z_bare + jωC_in(r)) gives Q′.
  - Predicted Q′/Q_A − 1 in [+20 %, +80 %].
  - Reasoning, closed form: 1/Q′ = 1/Q_0 + (1 − s_feed)/Q_C, where Q_0 = 971.07 and Q_C = 195.52 (`:2133`). s_feed, the feed capacitor's share of capacitor energy, is ≈ 0.25–0.5 in the k = 1 mode.
  - This misses the 5 % target by ≥ 4×. The ceiling is Q′ ≤ Q_0, i.e. ≤ +490 %.
- **(N-map) the single-frequency map reused:** S11 at 63.9 / 64.1 recomputed from `t39_esr`'s 64.000 MHz (R_p, C_p) is **not** run (it would cost two windows). Its size is printed as the gap-impedance error alone, \|Z_par(64-map, f) − Z_ser(f)\| / \|Z_ser(f)\| at 63.9 / 64.1 MHz: pure arithmetic, 1.60e-5 at both (the drafter's arithmetic). The parallel route's (Z′) sits at ≈ 1e-16, so the separation is ≈ 1e11×. That is ≥ 1e4× under `STRIP_BAND`, which the misused map (1.6e-5) would fail.
- **(N-nofeed) the feed branch dropped:** Z_in″ = Z_bare. Predicted to have **no** − → + crossing in (63.9, 64.1): the feed capacitor is part of the tuned loop. Printed "no crossing / f₀ undefined".

**Predictions (written first, printed).**
- W2 / W3 ≈ T₄₂₀ each (60–130 s).
- (Z′) ≈ 1e-16.
- Per-point S11 ties ≤ 1e-6.
- f₀ / Q / R_in gaps ≤ 1e-6.
- (N-feed) in [+20 %, +80 %]. (N-nofeed) no crossing.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; Palace as item 420). No probe edit after W1 starts (rule (i)).
- **Pre-check** (not a harness call): step 0 done. `ls logs/ops61/` must show `postpro_t39_esr/port-S.csv` and `t39_esr.json`, and no `t40_*`.
- **W1** (< 10 s):

  `scripts/testing/run_and_log.sh OPS-61-item421-W1-ans7-sweep-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py ans7-sweep-configs'"`
- **W2** (≈ T₄₂₀):

  `scripts/testing/run_and_log.sh OPS-61-item421-W2-palace-t40-f63.9 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 t40_esr_f63.9.json'"`
- **W3:** the same with `t40_esr_f64.1.json`, labelled `OPS-61-item421-W3-palace-t40-f64.1`.
- **W4** (≈ 3 s):

  `scripts/testing/run_and_log.sh OPS-61-item421-W4-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --ans7-sweep'"`
- **W5** (≈ 3 s), the dispatch regression: item 420's W6 command verbatim, labelled `OPS-61-item421-W5-regress-read`. Its output section must `diff` empty against item 420's W6 log.

**Traps already paid for.** All of item 420's apply, plus:
- **Never one Palace sweep with fixed R_p / C_p** under the parallel route (ruling (iii): R_p ∝ 1/f² off the conversion frequency). Every frequency gets its own map, in its own config.
- **The feed's (r + C) branch is in the read**, not in Palace (ruling 311-Q5). Dropping its r is (N-feed), the control, never the reading.
- Never overwrite `t39_*` or its `postpro_*`. Item 420's 64.000 MHz run is the middle point; it is not re-run.
- Docker socket and orphans: as item 420. Never lengthen a timeout or raise `-np`.

**Scope.** Covered: Palace's f₀ / R_in / Q at the ESR-tuned F-small node on the ×1 hole mesh, degree 2, read as T26 reads them. Not covered:
- an h-statement for Q (F31 / T41's);
- 127.74 MHz;
- the lossless rows;
- AED;
- an absolute tuned C, f₀, R_in or Q claim beyond "Palace = ours at matched order";
- `src/` and any band.

**Status it can move.**
- The `ANS-7` row's Q clause gains "T40 (item 421, record): Palace at matched order reads f₀ within `T21P_ZERO_RESOLUTION` and Q within T40's 5 % of ours at the ESR-tuned node, ×1, degree 2". That is the public external check §10 names. **No glyph:** the Done-when is Q's h-statement, which T41 / the 10-17 weekly writes. Rule (j): read; no mismatch.
- The `OPS-61` row (✅) gains a "step 3b record (T40)" sentence.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(Z′) red:** revert, BLOCKED.
- **(C1) red** (no crossing, or two): report Palace's Im Y_in triple. Park, BLOCKED; the review reads it against (T6).
- **(f₀) or (Q) red with (C1) green:** a disagreement on the external check. Report the gaps and the per-point S11 ties, open a known-issues entry, park, BLOCKED. Never widen, never re-map, never change r.
- **A GMRES non-convergence or abort:** that point is "not measured". With a point missing, (C1) cannot be read: BLOCKED with the log line.
- **W5 not diff-empty:** revert the probe, BLOCKED.

**Questions for the review.** (Q1) "T26's five nodes" is read as T26's Route A triple; the reasons are above. If the weekly meant a wider sweep, it would be a printed curve only, since `parallel_rows` reads the bracketing pair. (Q2) The 5 % Q target is quoted from §10 T40; confirm it as the item's asserted bound, or demote it to printed.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 421 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). Do step 0 first: if item 420's `t39_esr` Palace Total exceeds 120 s, or its convention reading is not 'parallel' or 'series', stop and mark the item BLOCKED for the 10-14 weekly. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; one frequency per Palace window, its own map per frequency. Probe only — no `src/`; never open `docs/private/` or `aed_results/`."

**Outcome (2026-10-13, 09:00 slot take-next, DONE — record; no glyph):** every anchor and control green first run, nothing killed. **Step 0:** item 420's W3 Elapsed Time Report `Total` Avg 51.765 s ≤ 120 s (`20261013T140240Z_OPS-61-item420-W3-palace-t39-esr.log:184`); W6 convention reading **parallel**, (T6-0) 2.565e-10 (`…140547Z_…-W6-read.log:46,49`) — the parallel route, each frequency its own exact map. The probe gained `T40_*` constants (with `T40_Q_TARGET = 0.05`'s provenance comment), `ans7_sweep_configs()` (`ans7-sweep-configs`) and `read_ans7_sweep()` (`read --ans7-sweep`), additive, dispatch branches first in `__main__`; last edited before W1 (14:10:42 UTC). Pre-check: `fem-em-solver` Up; `t39_esr.json` and `postpro_t39_esr/port-S.csv` present; no `t40_*`. (f₀) and (Q) are asserted with the bands imported / quoted (ruling (3), provisional for the weekly); the S11 ties and R_in are printed (ruling (1)).
- **W1 (Z′):** 63.9 MHz Q 195.83027908, R_p 7 670.09964059 Ω, C_p 6.3591359612e-11 F, \|Z_par − Z_ser\|/\|Z_ser\| 7.087e-19; 64.1 MHz Q 195.21926416, R_p 7 622.31222019 Ω, C_p 6.3591349216e-11 F, 1.422e-18 — both ≤ `STRIP_BAND`; each config differs from `t39_esr.json` only in `Freq`, ports 2–4 `R` / `C` and `Problem.Output` (asserted). (N-map, printed) the 64-MHz map at 63.9 / 64.1 MHz: 1.599e-5 / 1.597e-5 (predicted 1.60e-5), ≥ 1e13× the (Z′) residual.
- **W2 / W3 (executor-journalled):** 80 798 elements; `401: σ = 5.800e+07 S/m`; ND (p = 2) 570 332; ω/2π = 6.390e-02 / 6.410e-02 GHz; port 1 `Rs = 4.375e+01 Ω/sq`; ports 2–4 `Rs = 6.711e+03 Ω/sq, Cs = 7.268e-11 F/sq` / `Rs = 6.670e+03 Ω/sq, Cs = 7.268e-11 F/sq` (R_p·7/8 = 6 711.3 / 6 669.5); circuit lines `R = 7.670e+03 Ω` / `7.622e+03 Ω, C = 6.359e-11 F`; GMRES 33 / 33 its; Palace Total Avg 52.270 / 51.575 s (`palace.json` 52.286 / 51.593 s); peak 4.2 G; harness 54 / 53 s (predicted 60–130 s: under).
- **W4:** per-point \|S11_Pal − S11_A\|/\|S11_A\| 2.301e-10 / 2.569e-10 / 2.580e-10 (printed; predicted ≤ 1e-6); Palace Im Y_in −1.616350205e-04 / +2.676492848e-07 / +1.624710193e-04 S beside ours −1.616350144e-04 / +2.676560526e-07 / +1.624710261e-04 S (relative 3.7e-8 / 2.5e-5 / 4.2e-8 — the middle one is the near-zero susceptance, absolute 6.8e-12 S). **(C1)** one − → + crossing (bracket 63.9–64.0 MHz), none + → −: green. **(f₀)** f₀,Pal 63 999 834.68507022 Hz, f₀,Pal/f₀,A − 1 = **+6.530e-11** ≤ `T21P_ZERO_RESOLUTION` 1.287e-5; f₀,Pal/64 MHz − 1 = −2.583046e-6 (ours −2.583111e-6). **(Q)** Q_Pal 164.54343232572575, Q_Pal/Q_A − 1 = **−4.871e-9** ≤ T40's 5 %. R_in,Pal 3 175.979463819766 Ω, −7.456e-10 (printed). Ours' rows recomputed from the 12-digit S11_A prints: f₀ 63 999 834.68089134 Hz, Q 164.54343312645585 (the print-truncation floor, ≈ 1e-11 on f₀).
- **Controls:** (N-feed) feed's r omitted: Q′ 281.4602, Q′/Q_A − 1 = **+71.06 %** (predicted [+20 %, +80 %]: met; ≈ 1.5e8× the (Q) gap); f₀′ 63 999 421.55 Hz. (N-nofeed) feed branch dropped: Im Y_bare −2.5693e-2 / −2.5571e-2 / −2.5449e-2 S, no crossing, f₀ undefined (predicted: met).
- **W5:** item 420's W6 command, output section diff-empty (20 = 20 lines, 13 tagged lines identical).
- Logs `20261013T141042Z_OPS-61-item421-W1-ans7-sweep-configs.log`, `…141050Z_…-W2-palace-t40-f63.9.log`, `…141149Z_…-W3-palace-t40-f64.1.log`, `…141248Z_…-W4-read.log`, `…141257Z_…-W5-regress-read.log`.
- Disclosed gap (as items 418 / 420): the module-docstring dispatch paragraph and `usage:` string were not extended; the new functions carry their own docstrings.

### Review rulings 2026-10-13 03:00 (daily) on the item 418 / 420 / 421 drafts

1. **No matched-order band is registered by this review.** No constant for "OPS-61's matched-order band" exists, and a daily does not set bands (rule (6)). The cross-solver ties in items 418, 420 and 421 are **printed**, each beside the ≤ 1e-6 precedent; none is asserted. Item 420's question Q1 (register e.g. 5e-7) goes to the 10-14 weekly. What the items do assert is our imported identities on Palace's matrix and the ESR map's own identity.
2. **Item 418 is split** into a mesh export plus single-excitation Palace windows, each ≤ 480 s, with its price gate and minute-50 gate. Every command must stay within the 660 s foreground cap. The pre-registered rule for 64 / 128 MHz (heavy, or an `xl` cost probe) is a recommendation for the 10-14 weekly; the weekly commissions any `xl`.
3. **Item 421 uses T26's three tune points** (63.9 / 64.0 / 64.1 MHz; no five-node set exists), with the 5 % Q bound taken from §10 T40's text. This is accepted provisionally for the weekly to ratify (Q1 / Q2). Its step 0 returns the tier to the weekly if item 420 measured > 120 s per frequency.
4. **Palace availability** could not be checked from the review session. Each item's first command is a bare `docker compose -f docker/docker-compose.yml ps` / image check, and a missing image means the item stops and reports. No compose edit is made (the §9 compose constraint).

### §9 item 434 — 2026-10-13 (§10 chain F34, 64 MHz; chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-13 10:30 daily review; the review reviews it.* Licence: §10 chain F, F34 (PROJECT_PLAN.md, grep "**F34 (lane 1 → heavy)**"), quoted, not changed: "Palace at matched order on the h = 0.0075 exact mesh (entry 27's), 10 MHz first: `OPS-61`'s configuration with the resolution knob; target: matched-order agreement with our 0.0075 rung at `OPS-61`'s band (the 0.0000 % precedent at ×1); the measured price (`OPS-61` W4 268 s at ×1; predicted ≈ 800 s at 2.4× cells) decides whether 64 / 128 MHz run heavy or as an `xl` cost probe the 10-14 weekly writes. Serial on nothing." Item 418 ran the 10 MHz half and its pre-registered price rule read **heavy, four single-excitation windows per frequency** (max_k projected T_e(128) 271.3 s ≤ 450 s; Outcome above). This item is F34 as written at 64 MHz — a chain step, not provisional. The Palace directive (§5.4 (i)–(v)), the sheets-and-κ ruling (exact rectangles, κ off) and the 03:00 review's rulings (1) / (2) above bind. **No band, tolerance, Done-when or validation target is set or moved.**

**`OPS-61` step 1g: Palace order 2 on `ANS-4`'s exact-strip h = 0.0075 rung at 64 MHz, against our degree-2 record of that rung.** Our side is **not re-solved**: it is read from `xl` entry 27's tracked log. Item 418's mesh is **not re-exported**: it is on disk and reused.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Configs and the read run at `-n 1` in the complex container; Palace runs at `-np 8`. **Class:** chain step (solve), ≈ 9 measured min (class median).
- **Time:** ≈ 45 predicted min: 15 fixed + ≈ 10 code + W1 < 1 + W2–W5 4 × ≈ 3.3 (192 s central, 160–330 s each) + W6 / W7 < 1 + ≈ 5 spare.
- **Depends on:** nothing (item 418 landed; item 435 is independent — see "Shared code"). Inputs on disk (checked 2026-10-13) or tracked:
  - `logs/ops61/ops61_ans4_h0075_strip.msh` and `f05x_h0075_o2_f10_e{1..4}.json` (item 418 W1), `postpro_f05x_o2/port-S.csv` (Palace ×1 order 2, 64 MHz, item 194);
  - ours at h = 0.0075, degree 2, κ off, 64 MHz: the 10-digit S lines `docs/testing/logs/20261005T070009Z_ANS-4-step6b.log:3708–3712` (header `:3708`, `S_1k` … `S_4k` at `:3709–3712`; `:3709` carries eight leading `PASSED`), σ_max 0.999823496773 at `:3713` (and `(G)` `:3851`), ours' ×1 → 0.0075 h-move 0.2345 / 0.0950 / 0.2098 % at `:3727`;
  - ours at ×1, degree 2: `ANS4_EXACT_STRIP_S_P2_X1["64"]` (`tests/validation/test_ans4_exact_strip_records.py:262`).

**Price (from item 418's measurement, not a scaling guess).** Item 418 W2–W5 measured Total / Setup (Avg) 146.2 / 34.2 … 143.4 / 33.5 s at 10 MHz, 25 GMRES iterations, peak ≈ 14.4 G. Its read projected a 64 MHz single excitation at 189.8–193.7 s (r₆₄ = 37 / 26, the ×1 order-2 iteration ratio). Each window is one excitation, ≤ 590 s container, well inside the 660 s foreground cap and §5.1's ≈ 480 s safe ceiling. Four windows ≈ 13 min compute.

**Shared code (written once by whichever of items 434 / 435 runs first; probe only; additive; disclosed).**
- **Pre-check for the code:** `grep -n "h0075-configs" scripts/probes/ops61_step1a_palace_ans4.py`. If the dispatch exists (item 435 landed first), **make no probe edit** and go to the windows. Otherwise write the following, covering both frequencies:
- **Constants**, beside `H0075_OURS_H_MOVE` (`:722`):
  - `H0075_FREQS = {"64": 0.064, "128": 0.128}` (GHz);
  - `H0075_FREQ_SET = {f: {f"f05x_h0075_o2_f{f}_e{k}": k for k in (1, 2, 3, 4)} for f in H0075_FREQS}`;
  - `H0075_OURS_BLOCKS`, per frequency: the header string (`"[ANS-4 step5] 64 MHz kappa-off S (10 digits):"` / `"… 128 MHz …"`), the printed σ_max (0.999823496773, `:3713`; 0.999320852426, `:3720` / `(G)` `:3859`), ours' printed h-move (`:3727` / `:3728`) and the Palace ×1 order-2 output (`postpro_f05x_o2` / `postpro_f05x_o2_f128`), each with its log-line provenance comment.
- **New `h0075_freq_configs(f)`, dispatched by `h0075-configs --f64` / `--f128`.**
  1. Assert `comm.size == 1` (or no MPI). Assert the four target JSONs do not exist (never overwrite).
  2. **(E1′)** parse `H0075_MESH_NAME`'s Gmsh 2.2 `$Elements` block and assert the type-4 (tet) count == `H0075_CELL_RECORD` (282 198); print the file's sha256 and size.
  3. For each k, load `f05x_h0075_o2_f10_e{k}.json` and check it: `Model.Mesh == H0075_MESH_NAME`, four ports `[21]…[24]` `+Z` R 50, `Excitation` true on k only, `Order == 2`, `Freq == [0.010]`. Write `f05x_h0075_o2_f{f}_e{k}.json` changing **only** `Freq` → `[H0075_FREQS[f]]` and `Problem.Output` → `postpro_<name>`.
  4. Read each written JSON back and assert it differs from its template in those two keys only (item 421's template-diff pattern).
- **New `_ours_rung_s(log, header)`**, a copy of `_ours_h0075_s` (`:1999`) with the log and header as arguments. **`_ours_h0075_s` is not edited.**
- **New `read_exact_strip_h0075_freq(f)`, dispatched by `read --exact-strip-h0075 --f64` / `--f128`.** A copy of `read_exact_strip_h0075` (`:2022`, provenance comment) with the frequency's set, ours' block, σ_max, h-move and ×1 Palace output substituted, and the price section reading `*_OPS-61-item{434|435}-W[2-5]-palace-h0075-f{f}-e*.log` and printing measured Total / Setup / iterations beside item 418's projection for that frequency. **`read_exact_strip_h0075` is not edited** (W7 proves it).
- It never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** four new `sys.argv` branches beside `:2598–2630`. Every existing mode keeps its behaviour.

**Anchors (asserted; bands imported, never restated).**
- **(E1′), W1:** the reused mesh's tet count == 282 198, parsed from the file.
- **(Tie), W6:** the parsed 4×4's σ_max equals the printed 0.999823496773 within `H0075_TIE_BAND` 1e-9 absolute (a parse check, not a physics band). A wrong block misses by ≥ 1.7e-4 (10 MHz block 0.999995526884; 128 MHz 0.999320852426). Red ⇒ wrong block: stop.
- **(P), W6, on Palace's assembled 4×4** (`_p_check_s`, item 418's): ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND`; σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` and every column power ≤ 1 + 1e-9; every S-class spread ≤ `ADJACENT_SPREAD_BAND`; pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst. Reciprocity is cross-window (S_ij and S_ji from different Palace commands).
- **Palace log check (executor, journalled), W2–W5:** 282 198 elements, box ±0.12 / ±0.12 / ±0.10 m; `Rs = 4.375e+01 Ω/sq` on 21–24; PEC on 1 only; ω/2π = 6.400e-02 GHz; exactly one excited port, k; one `GMRES solver converged` line; `ND (p = 2): 1822094` (= ours, as item 418).

**Ours against Palace (printed, never gated — ruling (1) above).**
- **(T5)** per class and per entry, \|Palace o2 − ours p2\|/\|ours p2\| on the rung. **Predicted ≤ 1e-6** (item 418 measured ≤ 1.9e-10 at 10 MHz; floor the 10-digit print, ≈ 1e-10).
- **(C-h)** Palace o2 at h = 0.0075 against Palace o2 at ×1 (`postpro_f05x_o2`), per class, beside ours' 0.2345 / 0.0950 / 0.2098 %. Predicted equal to the printed digit.

**Negative controls (rule (e)).**
- **(C2-h) mesh-mismatch separation — predicted, printed.** Per class \|ours ×1 p2 (`ANS4_EXACT_STRIP_S_P2_X1["64"]`) − Palace o2 h\| / \|ours h p2 − Palace o2 h\|. Numerator predicted ≈ the h-move (weakest, adjacent, 9.5e-4 relative); denominator (T5) ≤ 1e-6 ⇒ **predicted ≥ 950×** against item 229's ≥ 25× bar. Ceiling set by the denominator's rounding floor: ≈ 9.5e-4 / 2.5e-10 ≈ 4e6×. Item 418 measured 1.3–2.4e7× at 10 MHz.
- **(P)'s power to discriminate** was shown by item 159's flip (200.0204 % spreads against the 0.5 % band, 400×); not re-run.
- **A wrong-frequency config** reads as (Tie) green with (T5) ≈ the 10 → 64 MHz S change (order 1), and the Palace log check's ω/2π line fails first.

**Predictions (written first, printed).**
- **W1:** < 10 s; (E1′) 282 198.
- **W2–W5, each:** Palace Total ≈ 192 s (160–330); GMRES 33–45 iterations (×1 order 2 at 64 MHz: 37; h = 0.0075 at 10 MHz: 25 × 1.423 ≈ 36); peak 14–20 G.
- **(P):** reciprocity ≈ 1e-11; spreads ≲ 0.05 % (ours on the rung: 0.0154 / 0.0157 / 0.0158 %, `:3713`); σ_max 0.999823497 ± 1e-9.
- **(T5)** ≤ 1e-6; **(C-h)** equal to ours' h-move to the printed digit; **(C2-h)** ≥ 950×.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call): `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up; `ls logs/ops61/`: stop and report if any `f05x_h0075_o2_f64_*.json` or `postpro_f05x_h0075_o2_f64_*` exists; `ops61_ans4_h0075_strip.msh`, `f05x_h0075_o2_f10_e{1..4}.json` and `postpro_f05x_o2/port-S.csv` must be present.
- **W1** (< 10 s), configs and (E1′):

  `scripts/testing/run_and_log.sh OPS-61-item434-W1-h0075-configs-f64 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py h0075-configs --f64'"`
- **W2** (≈ 192 s, ≤ 590), excitation 1, the price check:

  `scripts/testing/run_and_log.sh OPS-61-item434-W2-palace-h0075-f64-e1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_h0075_o2_f64_e1.json'"`
- **Price gate (binding, after W2).** Run W3–W5 iff W2's Palace `Total` (Avg) T₁ ≤ 450 s. Otherwise columns 2–4 are "not measured — price"; W6 still runs on column 1 ((T5) per entry, (P) "not measured"); the item lands as a partial record, not BLOCKED.
- **W3 / W4 / W5:** W2's command with `…_f64_e2.json` / `…_e3.json` / `…_e4.json`, labelled `OPS-61-item434-W3-palace-h0075-f64-e2` / `-W4-…-e3` / `-W5-…-e4`. **Minute gate:** if W5 cannot start by minute 50, stop launching; the `postpro_*` directories are durable, so the next slot runs the missing excitations with no edit and then W6 (rule (d), "partial; columns k… pending"). After any kill, `docker compose -f docker/docker-compose.yml ps` before the next window.
- **W6** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item434-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-h0075 --f64'"`
- **W7** (≈ 3 s), the dispatch regression: item 418's W6 command verbatim (`20261013T124505Z_OPS-61-item418-W6-read.log:12`), labelled `OPS-61-item434-W7-regress-read`. Its output section must `diff` empty against that log's. (If item 435 landed first and this item made no probe edit, W7 still runs.)

**Traps already paid for.** All of item 418's apply, plus:
- **Never re-run `export --exact-strip-h0075`** — the mesh exists and the function asserts its outputs are absent; a fresh mesh would also break (E1′)'s identity with item 418's columns. Never run `export`, `export --exact-strip`, `export --hole-pec`, or any `*-configs` mode other than `h0075-configs --f64`.
- **Parse trap:** `:3709` starts with eight `PASSED`; the 128 MHz block's spread line (`:3720`) is pytest-prefixed. Anchor on the 64 MHz header only; (Tie) catches a wrong block.
- One excitation per Palace command; never raise `-np`, never lengthen a timeout in-slot. A Status 124 is a price finding: record it, run no further e-window. A `GMRES … did not converge` line is a finding, not a retry: that column is "not measured".
- Orphans after any killed window: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0; `docker compose -f docker/docker-compose.yml ps` shows no lingering `fem-em-palace`.
- Compare S only; never read `z_matrix_ohm` or the `aed` key; never re-register a record from Palace. Palace's numbers go in full (v0.18.1, the four configs, the mesh); **no Palace figure set against an AED one in a tracked file**. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace order 2 at 64 MHz on the h = 0.0075 exact-strip rung; (P) on the assembled 4×4; the matched-order gaps against entry 27's record; the measured price beside item 418's projection. Not covered: 128 MHz (item 435); 10 MHz (item 418); Palace order 1 on the rung; the h = 0.005 rung (item 437); any h-statement for Palace; the hole / copper columns; the private reading; `src/` and any band.

**Status it can move.**
- The `OPS-61` row (✅) gains a "step 1g record (F34 64 MHz, §9 item 434)" clause: (P) and (T5) per class. **No glyph change**; rule (j): the row's Done-when is step 1c's and is met.
- With item 435, §10 F34 is complete at 10 / 64 / 128 MHz; the §2 / CLAUDE.md Palace clause may extend to "on two rungs (×1, h = 0.0075) at 10 / 64 / 128 MHz" — the weekly's call.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(E1′) red:** the on-disk mesh is not item 418's. Stop before any Palace window, park on `attempt/*`, BLOCKED (rule (d)); report the count and sha256.
- **Palace aborts, or a W2–W5 log lacks 282 198 elements / Rs 43.75 Ω/sq / 6.400e-02 GHz:** record verbatim, known-issues entry, BLOCKED.
- **(Tie) red:** a parse bug; fix once in-slot before W6 is re-run; still red, BLOCKED.
- **(P) red:** report every identity's number, park, BLOCKED. Never widen.
- **(T5) ≫ 1e-6 with (P) and (Tie) green:** a reading, not a failure — DONE as a record; the first matched-order disagreement at 64 MHz on a refined mesh, for the review.
- **W7 not diff-empty:** revert the probe, BLOCKED.
- **Price gate trips:** a partial record, not BLOCKED; the review re-prices.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 434 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). Check whether `h0075-configs` is already dispatched (item 435 may have landed first); if so, make no probe edit. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; one Palace excitation per window, never two heavy windows at once. Probe only — no `src/`, no test module; never re-export the mesh, never run `export`, `export --exact-strip` or `export --hole-pec`; never open `docs/private/` or `aed_results/`."

**Outcome (2026-10-13, 12:00 CDT slot, DONE — record; no glyph):** every anchor green on the first run, nothing killed. This item landed first, so it wrote the shared code (additive, covering both frequencies): `H0075_FREQS`, `H0075_FREQ_SET`, `H0075_OURS_BLOCKS` (with provenance comments, plus item 418's projected price per frequency), `H0075_ITEM`, `h0075_freq_configs(f)`, `_ours_rung_s`, `read_exact_strip_h0075_freq(f)`, and four dispatch branches (`h0075-configs --f64/--f128`, `read --exact-strip-h0075 --f64/--f128`). Item 435 makes no probe edit. The last probe edit was before W1 (17:02:01 UTC). Pre-check: `fem-em-solver` Up; mesh, the 10 MHz configs and `postpro_f05x_o2/port-S.csv` were present; no `f64` config or output existed. Palace v0.18.1, configs `logs/ops61/f05x_h0075_o2_f64_e{1..4}.json`, mesh `ops61_ans4_h0075_strip.msh`.
- **W1** (2 s, Status 0): (E1′) 282 198 type-4 tets parsed from `$Elements` (295 462 elements in total); mesh sha256 `bd93fb9c4dc25cb4294a38ad231aded3a4b5fcc040e50957550a49da9d4c22ae`. Each config differs from its 10 MHz template in `Freq` → [0.064] and `Problem.Output` only (asserted). **The W1 log is quarantined, not tracked:** the private-leak hook refused it because the mesh's printed byte size (an 8-digit integer) prefix-collides with a reference figure. This was confirmed by checking that token alone. Following item 417's precedent, the log moved unedited to the gitignored `logs/ops61_item434_quarantine/` (never deleted). Its test-results row stays. The (E1′) count is corroborated in every tracked W2–W5 Palace log (282 198 elements). Note for item 435: its W1 prints the same size and will hit the same refusal, so quarantine it the same way.
- **W2–W5 Palace log check:** each log shows 282 198 elements, box ±0.12 / ±0.12 / ±0.10 m, `Rs = 4.375e+01 Ω/sq` on 21–24, PEC on 1, ω/2π = 6.400e-02 GHz, exactly the one excited port k, one `GMRES solver converged` line (43 iterations each, inside the predicted 33–45), and `ND (p = 2): 1822094`. Total / Setup (Avg): 185.8 / 32.9, 181.0 / 32.1, 184.7 / 32.7, 184.2 / 33.2 s, against item 418's projection of 189.8–193.7 s. Peak per node ≈ 14.8 G. The price gate passed (T₁ 185.8 ≤ 450).
- **W6:** (Tie) σ_max 0.999823496793 against the printed 0.999823496773, \|Δ\| 1.96e-11. Each excitation wrote its own `S[i][k]` column. (P) PASS: reciprocity 2.018e-11, σ_max 0.999823496779, spreads self 0.0172 / adjacent 0.0161 / opposite 0.0139 %, column powers 0.836242679288 / 0.836259384510 / 0.836084337326 / 0.836216920023, pooled 42.1693 % = 2458.5× the worst. (T5) per class: self 4.757e-11 / adjacent 8.022e-11 / opposite 3.630e-11 ("0.0000 %"); per entry 2.6e-11 … 2.0e-10. (C-h): Palace's own ×1 → 0.0075 move is 0.2345 / 0.0950 / 0.2098 %, equal to ours to the printed digit. (C2-h) separation 4.926e7 / 1.183e7 / 5.768e7× (predicted ≥ 950×; this is above the brief's ≈ 4e6× ceiling estimate because the denominator landed at ≈ 3e-11 absolute, below the assumed 2.5e-10).
- **W7:** item 418's W6 command run verbatim; its output section is diff-empty against `20261013T124505Z_OPS-61-item418-W6-read.log`.
- **Disclosed gap** (same as item 418): the module docstring and the `usage:` string do not list the new modes.
- Logs: `20261013T170201Z_OPS-61-item434-W1-h0075-configs-f64.log` (quarantined, untracked), `…170207Z_…-W2-palace-h0075-f64-e1.log`, `…170523Z_…-W3-palace-h0075-f64-e2.log`, `…170828Z_…-W4-palace-h0075-f64-e3.log`, `…171138Z_…-W5-palace-h0075-f64-e4.log`, `…171447Z_OPS-61-item434-W6-read.log`, `…171458Z_OPS-61-item434-W7-regress-read.log`.

### §9 item 435 — 2026-10-13 (§10 chain F34, 128 MHz; chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-13 10:30 daily review; the review reviews it.* Licence: §10 chain F, F34, quoted under item 434 and not changed; item 418's price rule read **heavy** for 128 MHz (projected single excitation 265.9–271.3 s ≤ 450 s). F34 as written at 128 MHz — a chain step, not provisional. The same rulings bind. **No band, tolerance, Done-when or validation target is set or moved.**

**`OPS-61` step 1h: Palace order 2 on the h = 0.0075 exact-strip rung at 128 MHz, against our degree-2 record of that rung.** Item 434's brief **in full**, with only the following changed:
- **Independent of item 434.** If `h0075-configs` is not dispatched on main, write item 434's "Shared code" section verbatim (both frequencies); if it is, make no probe edit.
- **Class / time:** chain step (solve), ≈ 9 measured; ≈ 50 predicted min: 15 fixed + ≈ 10 code (0 if 434 landed) + W1 < 1 + W2–W5 4 × ≈ 4.5 (270 s central, 230–420 s each) + W6 / W7 < 1 + ≈ 5 spare. The minute-50 gate is more likely to bind here than in item 434; a partial landing is the planned fallback.
- **Inputs:** ours at 128 MHz, the 10-digit S lines `20261005T070009Z_ANS-4-step6b.log:3714–3718` (header `:3714`, `S_1k` … `S_4k` at `:3715–3718`; none `PASSED`-prefixed); σ_max 0.999320852426 printed on the pytest-prefixed line `:3720` (and `(G)` `:3859`); ours' h-move 0.2204 / 0.1590 / 0.2863 % at `:3728`; `ANS4_EXACT_STRIP_S_P2_X1["128"]` (`test_ans4_exact_strip_records.py:290`); Palace ×1 order 2 at 128 MHz `postpro_f05x_o2_f128/port-S.csv` (item 229).
- **(Tie):** σ_max vs 0.999320852426 within 1e-9; a wrong block misses by ≥ 5.0e-4 (64 MHz 0.999823496773).
- **Palace log check:** ω/2π = 1.280e-01 GHz; everything else as item 434.
- **(C-h):** against `postpro_f05x_o2_f128`, beside ours' 0.2204 / 0.1590 / 0.2863 %.
- **(C2-h):** numerator weakest (adjacent) ≈ 1.59e-3 relative ⇒ **predicted ≥ 1 590×** against the ≥ 25× bar; ceiling ≈ 1.6e-3 / 2.5e-10 ≈ 6e6×.
- **Predictions:** W2–W5 Palace Total ≈ 270 s (230–420); GMRES 48–65 iterations (×1 order 2 at 128 MHz: 55; 25 × 2.115 ≈ 53); peak 14–22 G. (P): spreads ≲ 0.05 % (ours 0.0140 / 0.0220 / 0.0150 %, `:3720`); σ_max 0.999320852 ± 1e-9. (T5) ≤ 1e-6.
- **Price gate:** as item 434 (T₁ ≤ 450 s). A W2 above 450 s would contradict item 418's projection (≤ 271.3 s) by ≥ 1.66×: record it as a price finding.
- **Windows:** item 434's with `f64` → `f128` and `item434` → `item435` throughout:
  - Pre-check: stop if any `f05x_h0075_o2_f128_*` config or output exists; `postpro_f05x_o2_f128/port-S.csv` must be present.
  - **W1:** `scripts/testing/run_and_log.sh OPS-61-item435-W1-h0075-configs-f128 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py h0075-configs --f128'"`
  - **W2:** `scripts/testing/run_and_log.sh OPS-61-item435-W2-palace-h0075-f128-e1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_h0075_o2_f128_e1.json'"`
  - **W3 / W4 / W5:** W2's command with `…_f128_e2.json` / `…_e3.json` / `…_e4.json`, labelled `OPS-61-item435-W3-palace-h0075-f128-e2` / `-W4-…-e3` / `-W5-…-e4`.
  - **W6:** `scripts/testing/run_and_log.sh OPS-61-item435-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-h0075 --f128'"`
  - **W7:** item 418's W6 command verbatim, labelled `OPS-61-item435-W7-regress-read`, output section diff-empty against `20261013T124505Z_OPS-61-item418-W6-read.log`. If item 434 has landed, also diff a re-run of item 434's W6 command (labelled `OPS-61-item435-W7b-regress-read-f64`) against item 434's W6 log.
- **Traps:** item 434's, with the parse trap restated: the 128 MHz block has **no** spread line directly after `:3718` — σ_max sits on `:3720` behind a pytest test id; the (Tie) literal is the constant, never re-parsed from that line.
- **Scope / Status / Negative result:** item 434's with 128 MHz for 64 MHz. With item 434, §10 F34 is complete at 10 / 64 / 128 MHz.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 435 per `docs/planning/chunks/OPS-61.md` (read item 434's brief and the latest Review ruling paragraph first). Check whether `h0075-configs` is already dispatched (item 434 may have landed first); if so, make no probe edit. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; one Palace excitation per window, never two heavy windows at once. Probe only — no `src/`, no test module; never re-export the mesh; never open `docs/private/` or `aed_results/`."

**Outcome (2026-10-13, 12:00 CDT slot second item, DONE — record; no glyph):** every anchor green on the first run, nothing killed. Item 434 had landed the shared code (`h0075-configs` dispatched), so **no probe edit** was made. Pre-check: `fem-em-solver` Up; the mesh, the 10 MHz configs and `postpro_f05x_o2_f128/port-S.csv` were present; no `f128` h0075 config or output existed. Palace v0.18.1, configs `logs/ops61/f05x_h0075_o2_f128_e{1..4}.json`, mesh `ops61_ans4_h0075_strip.msh`.
- **W1** (1 s, Status 0): (E1′) 282 198 type-4 tets parsed. The mesh sha256 is `bd93fb9c…4c22ae`, identical to item 434's. Each config differs from its 10 MHz template in `Freq` → [0.128] and `Problem.Output` only (asserted). **The W1 log is quarantined, not tracked**, as item 434 predicted: it prints the mesh byte size. The log moved unedited to the gitignored `logs/ops61_item435_quarantine/` (never deleted). Its test-results row stays, and the 282 198 count is corroborated in every tracked W2–W5 Palace log.
- **W2–W5 Palace log check:** each log shows 282 198 elements, box ±0.12 / ±0.12 / ±0.10 m, `Rs = 4.375e+01 Ω/sq` on 21–24, PEC on 1, ω/2π = 1.280e-01 GHz, exactly the one excited port k, one `GMRES solver converged` line, and `ND (p = 2): 1822094`. Total / Setup (Avg) were 266.5 / 35.2, 271.3 / 32.4, 274.6 / 33.8 and 284.3 / 33.4 s, against item 418's projection of 265.9–271.3 s (W4 and W5 ran 1–5 % above it). Peak per node ≈ 15.7 G. The price gate passed (T₁ 266.5 ≤ 450). **Prediction miss (printed, not an anchor):** GMRES took 72 iterations per excitation, outside the predicted 48–65. Even so, the wall time stayed in the predicted 230–420 s, so the iteration-ratio scaling (r₁₂₈) under-predicts iterations but not time on this rung.
- **W6:** (Tie) σ_max 0.999320852444 against the printed 0.999320852426, \|Δ\| 1.85e-11. Each excitation wrote its own `S[i][k]` column. (P) PASS: reciprocity 3.296e-11, σ_max 0.999320852432, spreads self 0.0110 / adjacent 0.0220 / opposite 0.0106 %, column powers 0.885491327339 / 0.885478209277 / 0.885363495014 / 0.885437070852, pooled 60.3004 % = 2745.6× the worst. (T5) per class: self 1.275e-10 / adjacent 1.189e-10 / opposite 7.497e-11 ("0.0000 %"); per entry 2.8e-11 … 3.0e-10. (C-h): Palace's own ×1 → 0.0075 move is 0.2204 / 0.1590 / 0.2863 %, equal to ours to the printed digit. (C2-h) separation 1.726e7 / 1.336e7 / 3.807e7× (predicted ≥ 1 590×; this is above the ≈ 6e6× ceiling estimate for the same reason as item 434: the absolute denominator is ≈ 1.6–9.9e-11).
- **W7:** item 418's W6 command run verbatim; its output section is diff-empty against `20261013T124505Z_OPS-61-item418-W6-read.log`. **W7b:** item 434's W6 command re-run; its output section is diff-empty against `20261013T171447Z_OPS-61-item434-W6-read.log`.
- With item 434, **§10 F34 is complete at 10 / 64 / 128 MHz** on the h = 0.0075 rung.
- Logs: `20261013T171923Z_OPS-61-item435-W1-h0075-configs-f128.log` (quarantined, untracked), `…171931Z_…-W2-palace-h0075-f128-e1.log`, `…172406Z_…-W3-palace-h0075-f128-e2.log`, `…172842Z_…-W4-palace-h0075-f128-e3.log`, `…173323Z_…-W5-palace-h0075-f128-e4.log`, `…173812Z_OPS-61-item435-W6-read.log`, `…173817Z_OPS-61-item435-W7-regress-read.log`, `…173823Z_OPS-61-item435-W7b-regress-read-f64.log`.

### §9 item 436 — 2026-10-13 (§10 chain T42, daily, 2026-10-13 — provisional, §9 rule (6); chain step, solve)

**Outcome 2026-10-13 (16:30 slot, item 2; DONE, record, provisional).** Every asserted anchor green; every prediction met. W1 (`20261013T214657Z_OPS-61-item436-W1-ans7-128-configs.log`): C_in(r) = `T35_C_IN_R_X1` 1.5275048999871143e-11 F bit-equal to `t127.74/tune_p2_esr0.2.json`, C_in(0) bit-equal to `T42_C_IN_0_LITERAL`; (Z′) 1.741e-16 / 1.742e-16 / 6.811e-19 against `STRIP_BAND`; Q_C 408.15087 / 407.83135 / 407.51234, R_p 33 317.627 / 33 265.483 / 33 213.461 Ω, C_p/C − 1 −6.00e-6; (Z0) bit-equal, no R; (N-map) 3.841e-6 / 3.838e-6 (printed); no `6.359` string in the log. W2–W5 (`…T214704Z_…W2-palace-t42-f127.74.log`, `…T214809Z_…W3-palace-t42-f127.64.log`, `…T214909Z_…W4-palace-t42-f127.84.log`, `…T215010Z_…W5-palace-t42-lossless.log`): 80 798 elements, σ = 5.800e+07 on 401, PEC on 1 only, ND (p = 2) 570 332, ω/2π 0.1276 / 0.1277 / 0.1278 / 0.1277 GHz, sheets Rs 2.915e+04 / 2.911e+04 / 2.906e+04 Ω/sq with Cs 1.746e-11 F/sq (lossless Cs only), GMRES 42 / 42 / 42 / 43 iterations, Palace Total 55.983 / 55.994 / 56.628 / 56.886 s, peak 4.4 G. W6 (`20261013T215115Z_OPS-61-item436-W6-read.log`): per-point \|S11_Pal − S11_A\|/\|S11_A\| 2.373e-10 / 2.306e-10 / 2.259e-10; (T6-0′) 2.308e-10; (K-esr′) move 2.181790e-03 both, \|Pal/ours − 1\| 5.648e-9; Palace Im Y_in −4.057389419e-05 / +1.196490895e-07 / +4.085113005e-05 S; **(C1)** one − → + crossing at 127 739 705.97524849 Hz; **(f₀)** f₀,Pal/f₀,A − 1 = +6.073e-11 ≤ `T21P_ZERO_RESOLUTION` (also within T22P), f₀,Pal/127.74 MHz − 1 = −2.301744e-06 (ours −2.301805e-06); **(Q)** Q_Pal 254.46926048070927, Q_Pal/Q_A − 1 = +5.444e-9 ≤ `T40_Q_TARGET` (vs `T35_Q_R_X1` +5.355e-9); R_in,Pal 9 790.703634300547 Ω, +2.906e-9. Controls: (N-feed) Q′ 360.604, +41.71 % (predicted [+20 %, +70 %]: met); (N-nofeed) no crossing (met). W7 (`20261013T215124Z_OPS-61-item436-W7-regress-read.log`): output section diff-empty against item 421 W4's. (Q1) stays with the weekly: the gap is ≈ 2e5× under both bands, so T21P vs T22P does not change this reading.

*Drafted by a general-purpose agent for the 2026-10-13 10:30 daily review; the review reviews it.* Licence: §9 rule (6) (PROJECT_PLAN.md, grep "**Rule (6), the chain-extension licence**"): the next numbered step of live §10 chain T toward its stated subgoal — Q's public external check (T40's "the Q clause's public external check while F31 settles its h") — tagged *(daily, 2026-10-13 — provisional)*; the 10-14 weekly ratifies or reverts. **T42 (provisional):** items 420 + 421 combined at the second Larmor frequency — Palace's twin and driven read of the ESR-tuned F-small at 127.74 MHz, against our ×1 r = 0.2 Ω tune there (T35, item 387; reproduced to the printed digits by item 419 W0). Target: **T40's, unchanged** — f₀ within the zero resolution, Q within T40's 5 % (`T40_Q_TARGET`). The capacitor-ESR ruling (i)–(iii), the sheets-and-κ ruling and the Palace directive apply. **No band, tolerance, Done-when or validation target is set or moved.**

**What exists on our side (the comparand; not re-solved).** The ×1 127.74 MHz ESR tune is T35 (`ANS-7` step 11b, item 387), registered literals `T35_C_IN_R_X1` / `T35_Q_R_X1` / `T35_R_IN_R_X1` (`tests/validation/test_ans7_exact_strip.py:2201–2205`, provenance `20261012T004412Z_ANS-7-item387-W2-esr128.log`). Item 419 W0 rewrote the ESR JSONs and reproduced them; its log is the latest window that wrote the files this item reads, so it is the comparand:
- `20261013T124951Z_ANS-7-item419-W0-esr128.log:393` — C_in(r) 1.5275048999871143e-11 F (= `T35_C_IN_R_X1`), C_in(0) 1.5274843607471887e-11 F from `t127.74/tune_p2.json`; `:394` (N-r0) bit-equal;
- `:2127–2129` — S11_A at (r, C_in(r)), 12 digits: 127.64 MHz 4.493753890384e-01+8.875155598136e-01j; 127.74 MHz 4.513684481906e-01+8.864778722343e-01j; 127.84 MHz 4.533642518309e-01+8.854328851717e-01j;
- `:2130` — Im Y_in,A −4.057389094e-05 / +1.196522461e-07 / +4.085113315e-05 S; rows A: f₀ 127 739 705.96749091 Hz, R_in 9 790.70360584452 Ω, Q 254.46925909529418;
- `:2131` — (T2-r′) −2.301805e-06 ≤ `T22P_ZERO_RESOLUTION`; `:2133` Q_0 618.338591513781, Q_C 407.8314; `:1279` cells 80 798.
- The lossless C path at 127.74 MHz, C_in(0): S11_A 4.523660701161e-01+8.884107356215e-01j (`20261006T171745Z_ANS-7-item283-W2-cin128.log:2055`; C_in 1.5274843607471887e-11 F at `:2058`).
- Route A's tune triple is `TUNE128_MHZ` = (127.64, 127.74, 127.84) (`test_ans7_exact_strip.py:943`); `parallel_rows` reads the bracketing pair, as item 421's reasoning (its Q1, ruled (3) above).

`OPS-61` step 3c: Palace order 2 on the ×1 hole mesh at 127.64 / 127.74 / 127.84 MHz with ports 2–4 carrying each frequency's exact parallel map of (r, C_in(r)), plus the lossless C path at 127.74 MHz; Palace's f₀, R_in and Q through `parallel_rows`, against T35's Route A rows.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Configs and the read at `-n 1`; Palace at `-np 8`. **Class:** chain step (solve), ≈ 9 measured.
- **Time:** ≈ 40 predicted min: 15 fixed + ≈ 12 code + W1 < 1 + W2–W5 4 × ≈ 1.3 (73 s central, 60–130 s each) + W6 / W7 < 1 + ≈ 5 spare.
- **Price (measured):** items 420 / 421 measured 51.1–53.2 s Palace Total per point at 64 MHz on this mesh (Setup 8.7 s, 33–35 GMRES iterations, peak 4.2–4.3 G). At 128 MHz the strip mesh's order-2 iterations rose 37 → 55 (×1.49), so ≈ 8.7 + 43 × 1.49 ≈ 73 s per point.
- **Depends on:** nothing (items 420 / 421 landed; their configs and outputs are templates, never overwritten). Inputs on disk (checked 2026-10-13): `logs/ops61/t39_esr.json`, `t39_lossless.json`, `ops61_ans6_x1_hole_strip.msh`; `logs/ans7_exact/t127.74/tune_p2.json` (`c_in`) and `tune_p2_esr0.2.json` (`c_in_r`).

**ESR map values (drafter's arithmetic, for the predictions; the probe recomputes).** At C = `T35_C_IN_R_X1`, r = 0.2 Ω: Q_C(127.64 / 127.74 / 127.84 MHz) = 408.15087 / 407.83135 / 407.51234; R_p = 33 317.63 / 33 265.48 / 33 213.46 Ω; C_p = 1.5274957e-11 F at each (C_p/C − 1 ≈ −6.01e-6). Palace's sheet values (w 7 mm, l 8 mm): Rs ≈ 2.915e+04 / 2.911e+04 / 2.906e+04 Ω/sq, Cs ≈ 1.746e-11 F/sq; lossless Cs = C_in(0)·8/7 ≈ 1.746e-11 F/sq. Series \|Z\| at 127.74 MHz 81.57 Ω.

**Code (probe only; additive; disclosed).**
- **Constants** beside `T40_*`: `T42_SWEEP = (("t42_esr_f127.64", 0.12764, 127.64e6), ("t42_esr_f127.74", 0.12774, 127.74e6), ("t42_esr_f127.84", 0.12784, 127.84e6))`; `T42_LOSSLESS = ("t42_lossless_f127.74", 0.12774, 127.74e6)`; `T42_S11_A` (the three `:2127–2129` values), `T42_F0_A`, `T42_R_IN_A`, `T42_Q_A`, `T42_IM_Y_A` (`:2130`), `T42_S11_A_LOSSLESS` (`…item283-W2-cin128.log:2055`), `T42_C_IN_0_LITERAL = 1.5274843607471887e-11` (`…item283-W2…:2058`; `…item419-W0…:393`) — each with its log-line provenance comment. `T40_Q_TARGET` is reused, not copied.
- **New `ans7_128_configs()`, dispatched by `ans7-128-configs`.**
  1. Assert `comm.size == 1` (or no MPI). Assert none of the four `t42_*.json` exists.
  2. Import `T35_C_IN_R_X1` from `tests.validation.test_ans7_exact_strip` (a module import, never a `test_*` name). Assert it equals `t127.74/tune_p2_esr0.2.json["c_in_r"]` bit for bit; assert `t127.74/tune_p2.json["c_in"]` equals `T42_C_IN_0_LITERAL` bit for bit.
  3. From `t39_esr.json` (checked as item 421 checks it: `Freq == [0.064]`, port 1 R 50 excited, ports 2–4 unexcited with R / C), write the three `t42_esr_*` configs changing **only** `Freq`, ports 2–4 `R` / `C` (`_t39_map(ω_f, T35_C_IN_R_X1, 0.2)` at each frequency) and `Problem.Output`.
  4. From `t39_lossless.json` (checked: ports 2–4 `C` only), write `t42_lossless_f127.74` changing only `Freq`, ports 2–4 `C` → C_in(0) at 127.74 MHz, and `Problem.Output`.
  5. Assert (Z′) / (Z0) below and the template diff (item 421's pattern), reading each JSON back.
- **New `read_ans7_128()`, dispatched by `read --ans7-128`.** Item 421's `read_ans7_sweep` (`:2475`) copied with the 127.74 sets and literals: S[1][1] only; Z_bare = `_ANS7._z_from_s11`; Z_in = 1/(1/Z_bare + 1/`_z_branch`(ω, `T35_C_IN_R_X1`, 0.2)) (the feed's own (r + C) branch in post-processing, ruling 311-Q5); rows = `parallel_rows([127.64e6, 127.74e6, 127.84e6], Z_in)`. Imports `T21P_ZERO_RESOLUTION`, `T22P_ZERO_RESOLUTION`, `T35_C_IN_R_X1`, `T35_Q_R_X1`, `_ANS7`, `_z_branch`, `parallel_rows`. **`read_ans7_sweep` is not edited.** Never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** two new `sys.argv` branches beside `:2598`.

**Anchors (asserted; bands imported or quoted from the target, never chosen).**
- **(Z′), W1:** per ESR config, \|Z_par − Z_ser(f)\| / \|Z_ser(f)\| ≤ `STRIP_BAND` (1e-9, `tests/validation/test_ops61_exact_strip_ports.py:58`), Z_ser(f) = `series_rlc_impedance(f, r_ohm=0.2, c_f=T35_C_IN_R_X1)`; and \|dω\|/ω ≤ `STRIP_BAND`.
- **(Z0), W1:** the lossless config's `C` == C_in(0) bit for bit, no `R` key.
- **(C1) one crossing:** Palace's Im Y_in changes sign − → + exactly once on the triple and never + → −; f₀,Pal inside (127.64, 127.84) MHz.
- **(f₀)** \|f₀,Pal / f₀,A − 1\| ≤ `T21P_ZERO_RESOLUTION` (1.287e-5, imported) — T40's band as the commission states it; printed beside `T22P_ZERO_RESOLUTION` (1.616e-5), the band T35's (T2-r′) asserted at 127.74 MHz. T21P is the tighter, so asserting it loosens nothing (question (Q1) below). f₀,A = 127 739 705.96749091 Hz.
- **(Q)** \|Q_Pal / Q_A − 1\| ≤ `T40_Q_TARGET` (0.05, T40's quoted target; ruling (3) above, provisional for the weekly). Q_A = 254.46925909529418.

**Printed, never gated (ruling (1)).** Per-point \|S11_Pal − S11_A\|/\|S11_A\| at the three points (predicted ≤ 1e-6; items 420 / 421 measured ≤ 2.6e-10); **(T6-0′)** the lossless C path at 127.74 MHz against `T42_S11_A_LOSSLESS` (predicted ≤ 1e-6); R_in,Pal against 9 790.70360584452 Ω; Palace's Im Y_in triple beside ours; f₀,Pal / 127.74 MHz − 1 beside ours' −2.301805e-06; Q_Pal beside `T35_Q_R_X1` (item 387's digits); Palace Total and iterations per run.

**Negative controls (rule (e); predicted, printed — no prior Palace measurement at 127.74 MHz).**
- **(N-feed) the feed's r omitted** (post-processing on the same three S11): Q′ from Z_in′ = 1/(1/Z_bare + jωC_in(r)). Closed form 1/Q′ = 1/Q_0 + (1 − s_feed)·(1/Q_A − 1/Q_0); with Q_0 618.34, Q_A 254.47 and s_feed 0.35–0.65 (64 MHz measured ≈ 0.5), **predicted Q′/Q_A − 1 in [+20 %, +70 %]**, ≥ 4× the 5 % target. Ceiling Q′ ≤ Q_0, i.e. ≤ +143 %. (64 MHz measured +71.06 %, item 421; at 127.74 MHz the ESR share of 1/Q is smaller, 0.59 against 0.83, so the move is predicted smaller.)
- **(N-map) the 127.74 MHz map reused at 127.64 / 127.84 MHz:** pure arithmetic, \|Z_par(127.74-map, f) − Z_ser(f)\| / \|Z_ser(f)\| = 3.84e-6 at both (drafter's arithmetic) — fails `STRIP_BAND` by ≈ 3 800×, against (Z′) ≈ 1e-16 (separation ≈ 1e10×). Not run as Palace windows.
- **(N-nofeed) the feed branch dropped** (Z_in″ = Z_bare): predicted **no** − → + crossing on the triple, "f₀ undefined" (64 MHz: met).
- **(K-esr′) r reaches the operator:** \|S11_Pal(esr, 127.74) − S11_Pal(lossless, 127.74)\| / \|S11_Pal(lossless)\| beside ours' same move from the two 12-digit lines, ≈ 2.18e-3 (drafter's arithmetic; the executor recomputes). Predicted equal to ≤ 1e-6 relative; a Palace that dropped `R` reads ≈ 1e-10 here.

**Predictions (written first, printed).** W1 < 5 s, (Z′) ≈ 1e-16. W2–W5 each Palace Total ≈ 73 s (60–130), GMRES 40–70 iterations (`MaxIts` 100), peak 4–5 G; Palace's sheet lines `Rs ≈ 2.915e+04 / 2.911e+04 / 2.906e+04 Ω/sq, Cs = 1.746e-11 F/sq` (lossless `Cs = 1.746e-11 F/sq`, no Rs). Per-point S11 ties ≤ 1e-6; (T6-0′) ≤ 1e-6; f₀ / Q / R_in gaps ≤ 1e-6. (N-feed) in [+20 %, +70 %]; (N-nofeed) no crossing; (K-esr′) ratio 1 ± 1e-6.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; Palace v0.18.1, image `fem-em-palace:0.18.1`, profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call): `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up; `ls logs/ops61/ logs/ans7_exact/t127.74/`: stop if any `t42_*.json` or `postpro_t42_*` exists; `t39_esr.json`, `t39_lossless.json`, the hole mesh, `tune_p2.json` and `tune_p2_esr0.2.json` must be present. Never re-create one.
- **W1** (< 10 s), configs, (Z′), (Z0):

  `scripts/testing/run_and_log.sh OPS-61-item436-W1-ans7-128-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py ans7-128-configs'"`
- **W2** (≈ 73 s), the tuned point first, the price check:

  `scripts/testing/run_and_log.sh OPS-61-item436-W2-palace-t42-f127.74 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 t42_esr_f127.74.json'"`
- **W3 / W4 / W5:** W2's command with `t42_esr_f127.64.json` / `t42_esr_f127.84.json` / `t42_lossless_f127.74.json`, labelled `OPS-61-item436-W3-palace-t42-f127.64` / `OPS-61-item436-W4-palace-t42-f127.84` / `OPS-61-item436-W5-palace-t42-lossless`. Price gate: run W3–W5 iff W2's Palace Total ≤ 300 s (≤ 4 × 300 s fits the slot with the read); otherwise the item lands partial and the review re-prices. After any kill, `docker compose -f docker/docker-compose.yml ps` before the next window.
- **W6** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item436-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --ans7-128'"`
- **W7** (≈ 3 s), the dispatch regression: item 421's W4 command verbatim (`20261013T141248Z_OPS-61-item421-W4-read.log:12`), labelled `OPS-61-item436-W7-regress-read`; its output section must `diff` empty against that log's.

**Traps already paid for.** All of items 420 / 421's apply, plus:
- **Never reuse a 64 MHz value at 127.74 MHz:** C_in(r) is `T35_C_IN_R_X1`, not `T26_C_IN_R_X1` (item 419's "silent misread": a 64 MHz literal under a 127.74 MHz print fails nothing). Grep the W1 log for `6.3593` / `6.3590`: there must be none.
- **Every frequency gets its own exact map** (ruling (iii), R_p ∝ 1/f² off the conversion frequency); never one Palace sweep with fixed R_p / C_p. Never enter 0.2 Ω as a parallel R (no short control is re-run: item 420 read the convention, parallel).
- **The feed's (r + C) branch is in the read, not in Palace** (ruling 311-Q5); dropping its r is (N-feed), the control, never the reading.
- **Compare S[1][1] only**; never read `z_matrix_ohm` or the `aed` key. Never overwrite `t39_*` / `t40_*` or their outputs; never re-run `ans7-configs` or `ans7-sweep-configs` (they assert their outputs absent).
- **Comparand digits:** item 387 W2 and item 419 W0 differ in the 11th–12th digit (Q 254.469259118 vs 254.469259095); the read uses item 419 W0's (the window that wrote the files) and prints `T35_Q_R_X1` beside it.
- Docker socket and orphans: as item 420. Never lengthen a timeout or raise `-np`. A `GMRES … did not converge` line is a finding: that point is "not measured", and with a tune point missing (C1) cannot be read.
- Record row and journal: Palace's numbers in full beside ours; **no Palace figure set against an AED one in a tracked file**. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace's bare S11, f₀, R_in and Q at the 127.74 MHz ESR-tuned F-small node, ×1 hole mesh, degree 2, Cu, read as T35 reads them; the lossless C-path tie at 127.74 MHz. Not covered: an h-statement for Q (F31 / T41); the h = 0.01 rung (T38's); 64 MHz (items 420 / 421); order 1; AED; a measured ESR; an absolute tuned C / f₀ / R_in / Q claim beyond "Palace = ours at matched order"; `src/` and any band.

**Status it can move.**
- The `ANS-7` row's Q clause gains "T42 (item 436, record, provisional): Palace at matched order reads f₀ within `T21P_ZERO_RESOLUTION` and Q within T40's 5 % of ours at the 127.74 MHz ESR-tuned node, ×1, degree 2" — the public external check at the second Larmor frequency. **No glyph:** the Done-when is Q's h-statement (T41). Rule (j): read; no mismatch.
- The `OPS-61` row (✅) gains a "step 3c record (T42)" clause.
- §10 T42 recorded DONE (provisional; the weekly ratifies or reverts).

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(Z′) / (Z0) red, or a bit-equality assert red in W1:** the wrong C or map. Stop before any Palace window, revert the probe, BLOCKED (rule (d)).
- **A W2–W5 log lacking 80 798 elements, σ = 5.800e+07 on 401, PEC on 1 only, ND (p = 2) 570 332, the right ω/2π:** park, BLOCKED.
- **(C1) red:** report Palace's Im Y_in triple; park, BLOCKED; the review reads it against the per-point ties.
- **(f₀) or (Q) red with (C1) green:** a disagreement on the external check. Report the gaps and the per-point ties, open a known-issues entry, park, BLOCKED. Never widen, re-map or change r.
- **(f₀) green against `T22P_ZERO_RESOLUTION` but red against `T21P_ZERO_RESOLUTION`:** not predicted (gaps ≈ 1e-10); if it happens, BLOCKED as above, and (Q1) goes to the weekly with the number.
- **(T6-0′) ≫ 1e-6:** a convention-free disagreement at 127.74 MHz: DONE as a record with a known-issues entry (as item 420's disposition), (f₀) / (Q) still read.
- **W7 not diff-empty:** revert the probe, BLOCKED.

**Questions for the review / weekly.** (Q1) The commission names T40's `T21P_ZERO_RESOLUTION`; the matched 127.74 MHz band is `T22P_ZERO_RESOLUTION` (T35's (T2-r′)). The brief asserts the tighter T21P and prints T22P; the weekly confirms which one T42 states. (Q2) T42 is chain T's first provisional step since the 10-12 weekly; it counts toward rule (6)'s three.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 436 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; one frequency per Palace window, its own map per frequency, C_in(r) = `T35_C_IN_R_X1` (never a 64 MHz literal). Probe only — no `src/`, no test module; never overwrite `t39_*` / `t40_*`; never open `docs/private/` or `aed_results/`."

### §9 item 437 — 2026-10-13 (§10 chain F36, daily, 2026-10-13 — provisional, §9 rule (6); chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-13 10:30 daily review; the review reviews it.* Licence: §9 rule (6) (PROJECT_PLAN.md, grep "**Rule (6), the chain-extension licence**"): the next numbered step of live §10 chain F, tagged *(daily, 2026-10-13 — provisional)*; the 10-14 weekly ratifies or reverts. **F36 (provisional):** F34 on the next rung — Palace at matched order (2) on `ANS-4`'s h = 0.005 exact-strip rung (`xl` entry 30, F27's finest rung), 10 MHz. Target: **F34's, unchanged** — matched-order agreement with our record of that rung, printed beside the ≤ 1e-6 precedent (ruling (1) above), our identities asserted on Palace's 4×4. It extends the Palace match from two rungs to the three F27's h-statement rests on. The Palace directive and the sheets-and-κ ruling bind. **No band, tolerance, Done-when or validation target is set or moved.**

**`OPS-61` step 1i: Palace order 2 on the h = 0.005 exact-strip rung at 10 MHz, against our degree-2 record of that rung.** Our side is **not re-solved**: it is read from `xl` entry 30's tracked log.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Export and read at `-n 1` in the complex container; Palace at `-np 8`. **Class:** chain step (solve), ≈ 9 measured.
- **Time:** ≈ 60 predicted min: 15 fixed + ≈ 10 code + W1 ≈ 3 (140–250 s) + W2–W5 4 × ≈ 4.5 (270 s central, 225–370 s each) + W6 / W7 < 1 + ≈ 5 spare. The minute-50 gate is likely to bind; a partial landing (columns pending, next slot finishes with no edit) is the planned fallback.
- **Depends on:** nothing. Inputs on disk or tracked:
  - `logs/ops61/f05x_o2_f10.json` (template, item 229), `postpro_f05x_h0075_o2_f10_e{1..4}/port-S.csv` (item 418, for (C-h));
  - ours at h = 0.005, degree 2, κ off, 10 MHz: `docs/testing/logs/20261011T070008Z_ANS-4-step6c.log:3815–3819` (header `:3815`, `S_1k` … `S_4k` at `:3816–3819`; `:3816` carries two leading `PASSED`), σ_max 0.999995535227 at `:3820` (and `(G)` `:3893`), cells 593 478 and 3 833 208 degree-2 unknowns at `:2280` / `:3814`, ours' ×1 → 0.005 move 0.2519 / 0.0377 / 0.1528 % at `:3840`;
  - ours at h = 0.0075 (item 418's comparand, `20261005T070009Z_ANS-4-step6b.log:3702–3706`) and `ANS4_EXACT_STRIP_S_P2_X1["10"]`.

**Can the mesh be exported in a slot? Yes.** The h = 0.005 exact-strip mesh is not on disk; it is built by the same `_build(np.zeros(LEG_COUNT), resolution=…, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` call item 418 used, at resolution 0.005. Entry 30 measured that build at 112.46 s (`step6c.log:2600`; gmsh runs serially on rank 0 at any width), against 52.15 s at 0.0075 (`step6b.log:2487`), whose `-n 1` export measured 61 s harness (item 418 W1). Predicted W1 ≈ 140–250 s at `-n 1`, under a 400 s container timeout; serial memory a few GB.

**Price (scaled from item 418's measurement; no step 0 probe needed).**
- Unknowns 3 833 208 / 1 822 094 = 2.104×. Item 418 measured 146 s per excitation (Setup 34 s) at 1 822 094; the ×1 → 0.0075 step (753 062 → 1 822 094, 2.42×) cost 1.68× per excitation (86.7 s → 146 s, sublinear).
- **Projected per excitation: ≈ 225 s (the sublinear ×1.55) to 307 s (linear) to 371 s (unknowns^1.25), central ≈ 270 s**; Setup ≈ 70–110 s of it. Every bracket is under 480 s (§5.1's safe ceiling) and far under the 600 s threshold for a separate probe, so W2 doubles as the price probe through the gate below.
- **Memory:** item 418 peaked at ≈ 14.4 G (8 ranks) ⇒ ≈ 30–40 G here, against the 128 GiB container: safe. A peak above 64 G is printed as a reading.

**Code (probe only; additive; disclosed).**
- **Constants** beside `H0075_*`: `H0050_RES = 0.005`; `H0050_MESH_NAME = "ops61_ans4_h0050_strip.msh"`; `H0050_CELL_RECORD = 593478` (literal, provenance `20261011T070008Z_ANS-4-step6c.log:2280, :3814`; no tracked constant holds it); `H0050_SET = {f"f05x_h0050_o2_f10_e{k}": k for k in (1, 2, 3, 4)}`; ours' log, header, σ_max 0.999995535227 and ×1 move, each with provenance.
- **New `export_exact_strip_h0050()`, dispatched by `export --exact-strip-h0050`.** A **copy** of `export_exact_strip_h0075` (`:1711`, provenance comment with commit and line span) with the H0050 constants substituted. **`export_exact_strip_h0075` is neither edited nor called.** Same steps: `comm.size == 1`; never overwrite; template check of `f05x_o2_f10.json`; build; **assert ncells == `H0050_CELL_RECORD` before anything is written**; write the mesh and four single-excitation configs (only `Model.Mesh`, `Excitation`, `Problem.Output` changed); (E1)–(E4).
- **New `read_exact_strip_h0050()`, dispatched by `read --exact-strip-h0050`.** A copy of `read_exact_strip_h0075` (`:2022`) using `_ours_rung_s(log, header)` (item 434's helper if it has landed; otherwise add that helper here, as item 434 specifies it) for the 0.005 block, with (C-h) against Palace's h = 0.0075 columns (`postpro_f05x_h0075_o2_f10_e{k}`, assembled) and ours' 0.0075 → 0.005 move computed from the two parsed logs (not a literal), and the price section reading `*_OPS-61-item437-W[2-5]-palace-h0050-e*.log`. **`read_exact_strip_h0075` is not edited.**
- Never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`. **Dispatch:** two new `sys.argv` branches.

**Anchors (asserted; bands imported, never restated).**
- **Export, W1:** **(E1)** tets == 593 478 exactly; **(E2)** exterior (attribute 1) area == `EXTERIOR_AREA_M2` within `EXTERIOR_AREA_BAND`; **(E3)** per port, N_i area == `SPEC_AREA_M2` and A/h == `SPEC_WIDTH_M` within `STRIP_BAND`; **(E4)** N_i ⊂ F_i.
- **(Tie), W6:** parsed σ_max vs printed 0.999995535227 within `H0075_TIE_BAND` 1e-9. A wrong block in this log misses by ≥ 1.7e-4; reading the 0.0075 log by mistake misses by 8.3e-9 (8.3× the band) — caught.
- **(P), W6,** on Palace's assembled 4×4: item 418's set (`RECIPROCITY_BAND`, `PASSIVITY_SIGMA_TOLERANCE` + column powers ≤ 1 + 1e-9, `ADJACENT_SPREAD_BAND`, `POOLED_SEPARATION_FLOOR`).
- **Palace log check (executor, journalled), W2–W5:** 593 478 elements; box ±0.12 / ±0.12 / ±0.10 m; `Rs = 4.375e+01 Ω/sq` on 21–24; PEC on 1 only; ω/2π = 1.000e-02 GHz; one excited port; one `GMRES solver converged`; `ND (p = 2): 3833208` predicted (= ours; a mismatch is a reading).

**Ours against Palace (printed, never gated).**
- **(T5)** per class and per entry on the rung; **predicted ≤ 1e-6** (item 418 ≤ 1.9e-10).
- **(C-h)** Palace o2 h = 0.005 against Palace o2 h = 0.0075, per class, beside ours' same move — 0.0144 / 0.0088 / 0.0161 % (self / adjacent / opposite; drafter's arithmetic from the two 10-digit blocks; the read recomputes). Predicted equal to the printed digit: the third rung of F27's sequence seen by both codes.

**Negative controls (rule (e)).**
- **(C2-h) mesh-mismatch separation — predicted, printed.** Per class \|ours h=0.0075 p2 − Palace o2 h=0.005\| / \|ours h=0.005 p2 − Palace o2 h=0.005\|. Numerator ≈ the 0.0075 → 0.005 move (weakest, adjacent, 8.8e-5 relative, 4.2e-5 absolute); denominator (T5) ≤ 1e-6 ⇒ **predicted ≥ 88×** against item 229's ≥ 25× bar. Ceiling ≈ 8.8e-5 / 2.5e-10 ≈ 3.5e5×. Printed also against ours ×1 (numerator ≈ the ×1 → 0.005 move, ≥ 3.8e-4 ⇒ ≥ 380×).
- **(P)'s power to discriminate:** item 159's flip; not re-run.
- **A wrong rung:** (E1) fails first; any other rung reads (T5) ≈ a mesh move (≥ 8.8e-5), not ≤ 1e-6.

**Predictions (written first, printed).** W1 140–250 s, (E2) to ≈ 1e-16, N_i to ≤ 2e-15. W2–W5 each Palace Total ≈ 270 s (225–370), GMRES 25–32 iterations (0.0075: 25), peak 30–40 G. (P): spreads ≲ 0.05 % (ours 0.0259 / 0.0110 / 0.0145 %, `:3820`), σ_max 0.999995535 ± 1e-9. (T5) ≤ 1e-6; (C-h) equal to ours' move to the printed digit; (C2-h) ≥ 88×.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call): `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up; `ls logs/ops61/`: stop if `ops61_ans4_h0050_strip.msh`, any `f05x_h0050_*.json` or any `postpro_f05x_h0050_*` exists; `f05x_o2_f10.json` and the four `postpro_f05x_h0075_o2_f10_e*/port-S.csv` must be present.
- **W1** (≈ 140–250 s), the export:

  `scripts/testing/run_and_log.sh OPS-61-item437-W1-export-h0050 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 400 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py export --exact-strip-h0050'"`
- **W2** (≈ 270 s, ≤ 590), excitation 1, **the price probe**:

  `scripts/testing/run_and_log.sh OPS-61-item437-W2-palace-h0050-e1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_h0050_o2_f10_e1.json'"`
- **Price gate (binding, after W2).** Run W3–W5 iff T₁ ≤ 450 s. Otherwise columns 2–4 are "not measured — price"; W6 runs on column 1 ((T5) per entry; (P) "not measured"); the item lands as a partial record and **64 / 128 MHz and the remaining columns return to the weekly** (an `xl` question). A Status 124 is the price finding: record it, no further e-window.
- **W3 / W4 / W5:** W2's command with `…_e2.json` / `…_e3.json` / `…_e4.json`, labelled `OPS-61-item437-W3-palace-h0050-e2` / `-W4-…-e3` / `-W5-…-e4`. **Minute gate:** if W5 cannot start by minute 50, stop launching; the next slot runs the missing excitations with no edit, then W6 (rule (d)).
- **W6** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item437-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-h0050'"`
- **W7** (≈ 3 s), the dispatch regression: item 418's W6 command verbatim (`20261013T124505Z_OPS-61-item418-W6-read.log:12`), labelled `OPS-61-item437-W7-regress-read`, output section diff-empty against that log's.

**Traps already paid for.** All of item 418's apply, plus:
- **The first 0.005 export runs at `-n 1`** (the function asserts it); a gmsh failure deadlocks wider runs. Never re-run `export --exact-strip-h0075` or any older export / configs mode.
- **Parse trap:** `:3816` starts with `PASSED PASSED`; other 10 MHz `S_1k` blocks exist earlier in the log (the PORT-9 prints at `:2605–2614`). Anchor on the header string only; (Tie) catches a wrong block or the wrong log.
- One excitation per Palace command; never raise `-np` (no `-np 12` Palace point exists) or lengthen a timeout in-slot. GMRES non-convergence is a finding: that column is "not measured".
- Orphans after any killed window: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0; `docker compose -f docker/docker-compose.yml ps` shows no lingering `fem-em-palace`.
- Compare S only; never re-register a record from Palace; **no Palace figure set against an AED one in a tracked file**; never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace order 2 at 10 MHz on the exported h = 0.005 exact-strip rung; (P) on the assembled 4×4; the matched-order gaps against entry 30's record; Palace's own 0.0075 → 0.005 move; the measured price. Not covered: 64 / 128 MHz on this rung (the weekly prices them from this item's readout); Palace order 1; a Palace h-statement or extrapolation; the hole / copper columns (F31's rung); the private reading; `src/` and any band.

**Status it can move.**
- The `OPS-61` row (✅) gains a "step 1i record (F36, §9 item 437, provisional)" clause: Palace on the third rung, (P), (T5) per class, (C-h). **No glyph change**; rule (j): Done-when met.
- If (T5) ≤ 1e-6 in every class, the §2 / CLAUDE.md Palace clause may read "on all three rungs of F27's h-statement" at 10 MHz — the weekly's call. §10 F36 recorded DONE (provisional).

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(E1) red:** the export is not entry 30's mesh. Stop before any Palace window, park on `attempt/*`, BLOCKED; a known-issues entry only if the count moved on an unchanged image (a re-record is `record-reconciler`'s, under (1*)).
- **(E2)–(E4) red:** the strip did not reach the mesh at this resolution; report spans and facet counts, park, BLOCKED.
- **W1 Status 124 (export > 400 s):** a price finding — record it, park, return the export's tier to the weekly. Never lengthen the timeout.
- **Palace aborts, or a W2–W5 log lacks 593 478 elements / Rs 43.75 Ω/sq:** record verbatim, known-issues entry, BLOCKED.
- **(Tie) red:** fix the parse once in-slot before W6 is re-run; still red, BLOCKED.
- **(P) red:** report every identity's number, park, BLOCKED. Never widen.
- **(T5) ≫ 1e-6 with (P) and (Tie) green:** a reading, DONE as a record, for the review.
- **W7 not diff-empty:** revert the probe, BLOCKED.
- **Price gate trips:** a partial record, not BLOCKED.

**Outcome 2026-10-13 16:30 slot (DONE, record, provisional).** All windows ran, none edited after W1. W1 export 118 s harness (mesh 112.54 s): (E1) 593 478 tets, (E2) exterior |diff| 5.551e-17, (E3) N_i area rel ≤ 8.882e-16 and A/h rel ≤ 8.882e-16, (E4) subset — PASS (`20261013T215455Z_OPS-61-item437-W1-export-h0050.log`). W2–W5 Palace Total 334.2 / 311.5 / 332.8 / 304.2 s (Setup 121.9 / 94.0 / 104.7 / 96.4 s), 25 GMRES iterations each, peak 32.8 / 32.3 / 32.2 / 32.2 G; every log: 593 478 elements, box ±0.12 / ±0.12 / ±0.10 m, Rs 4.375e+01 Ω/sq on 21–24, PEC on 1, ω/2π = 1.000e-02 GHz, one excited port, `ND (p = 2): 3833208` (= ours). Price gate: T₁ 334 s ≤ 450 s. W6 (`20261013T221851Z_OPS-61-item437-W6-read.log`): (Tie) |diff| 4.837e-11 PASS; (P) PASS — reciprocity 1.654e-11, σ_max 0.999995535243, spreads self 0.0280 / adjacent 0.0081 / opposite 0.0174 %; (T5) per class 6.308e-11 / 5.760e-11 / 1.679e-11, per entry ≤ 1.674e-10 (predicted ≤ 1e-6); (C-h) Palace 0.0075 → 0.005 move 0.0144 / 0.0088 / 0.0161 %, ours' computed move identical to the printed digit (predicted); (C2-h) 2.282e6 / 1.53e6 / 9.611e6× against ours h = 0.0075 (predicted ≥ 88×), 3.984e7 / 6.536e6 / 9.088e7× against ours ×1. Printed 64 / 128 MHz single-excitation projections (item 418's iteration ratios) 392–429 s / 536–587 s — the weekly tiers them (ruling (5)). W7 (`20261013T221857Z_OPS-61-item437-W7-regress-read.log`) output section diff-empty against item 418's W6.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 437 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; the export at `-n 1`, one Palace excitation per window, never two heavy windows at once; the price gate after W2 and the minute-50 gate bind. Probe only — no `src/`, no test module; never run an older `export` mode; never open `docs/private/` or `aed_results/`."

### Review ruling 2026-10-13 10:30 (daily) on the item 434–437 drafts

1. **434 / 435 execute F34 as the weekly wrote it** — item 418's measured price put 64 / 128 MHz on the heavy tier, which is F34's own decision rule; not provisional. Independent of each other; whichever runs first lands the shared reader modes, the second only checks for them.
2. **436 (T42) and 437 (F36) are rule (6) provisional steps** (the first for each chain since the 10-12 weekly; one of three each), written into §10 in this review's commit. F36 serves chain F's live public-table subgoal (the h-table's 0.005 rung read by the public external reference); the weekly ratifies or reverts both.
3. **436's f₀ band** — accepted as drafted: asserted at `T21P_ZERO_RESOLUTION` (T40's band, as chain T states the target), `T22P_ZERO_RESOLUTION` (the band T35 asserted at 127.74 MHz) printed beside; a red only against T21P parks the item BLOCKED for the weekly, never re-registered in-slot. The weekly confirms which band T42 states.
4. **Cross-solver band** — still none: every Palace-vs-ours gap is printed, as in items 418 / 420 / 421; the weekly's (carried from item 418).
5. **If 437 lands green**, F36 at 64 / 128 MHz on the 0.005 rung is the weekly's to tier (heavy at ≈ 4.5–7 min per excitation, or an `xl` cost probe); not queued here.

### §9 item 444 — 2026-10-14 (§10 chain T43, daily — provisional)

**Outcome 2026-10-14 (04:30 CDT slot, DONE — record, provisional; no glyph).** Every asserted anchor green on the first run; nothing killed; no probe edit after W1 (last edit before 09:32:58 UTC). The probe gained `H001_*` / `T43_*` constants (provenance comments), `export_hole_h001()` (`export --hole-h001`, a mesh-only copy of `export_hole_pec` with `resolution=H001_RES`), `ans7_h001_configs()` (`ans7-h001-configs`) and `read_ans7_h001()` (`read --ans7-h001`), additive, dispatch branches first in `__main__`. W1 (`20261014T093258Z_OPS-61-item444-W1-export-hole-h001.log`): mesh built in 31.39 s, export 32.70 s; (E1) 140 270 tets before the write and 140 270 type-4 elements parsed back; (E3) census exterior 31 702 / cavity 24 464 / outer 7 238 / collision 0; (E4) box area \|d\| 1.110e-16; (E5) N_i area rel 8.882e-16, width rel ≤ 8.882e-16, height 0; (E6) tag set [2, 3, 101–104, 201–204] = `t39_esr.json`'s Materials. W2 (`…T093341Z_…-W2-ans7-h001-configs.log`): C_in(r, h) / C_in(0, h) bit-equal to both tune JSONs; (Z′) 1.816e-16 / 7.104e-19 / 2.134e-18; Q_C 195.66692417 / 195.36119460 / 195.05641894, R_p 7 657.30904 / 7 633.39927 / 7 609.60131 Ω; (Z0) bit-equal, no R; (N-map) 1.601e-5 / 1.598e-5; no `6.3593` / `6.3590` string. W3–W6 (`…T093349Z_…-W3-palace-t43-f64.0.log`, `…T093534Z_…-W4-palace-t43-f63.9.log`, `…T093713Z_…-W5-palace-t43-f64.1.log`, `…T093852Z_…-W6-palace-t43-lossless.log`; Palace v0.18.1, `-np 8`, mesh `ops61_ans6_h0010_hole_strip.msh`): 140 270 elements, box ±0.12 / ±0.12 / ±0.10 m, σ = 5.800e+07 on 401, PEC on 1 only, port 1 Rs 4.375e+01 Ω/sq, ports 2–4 Rs 6.679e+03 / 6.700e+03 / 6.658e+03 Ω/sq with Cs 7.274e-11 F/sq (lossless Cs only), ω/2π 6.400e-02 / 6.390e-02 / 6.410e-02 / 6.400e-02 GHz, `ND (p = 2): 968016` (= ours), GMRES 40 / 39 / 40 / 43 iterations, Palace Total 94.502 / 95.097 / 93.784 / 96.910 s (Setup ≈ 16.2–16.8 s), peak 7.5–7.6 G. Price gate: W3 94.5 s ≤ 300 s. W7 (`20261014T094041Z_OPS-61-item444-W7-read.log`): (T6′) 2.591e-10 / 2.533e-10 / 2.507e-10; (T6-0″) 2.499e-10; (K-esr″) 4.954501e-03 both, \|Pal/ours − 1\| 1.717e-9; Palace Im Y_in −1.617828228e-04 / +2.684615079e-07 / +1.626210717e-04 S; **(C1)** one − → + crossing at 63 999 834.335464224 Hz; **(f₀)** +6.437e-11 ≤ `T21P_ZERO_RESOLUTION`, f₀,Pal/64 MHz − 1 −2.588508e-06 (ours −2.588573e-06); **(Q)** Q_Pal 164.21576013522352, +1.829e-9 ≤ `T40_Q_TARGET`; R_in,Pal 3 166.747989302149 Ω, +6.507e-10. **(C-h)** Palace's ×1 → 0.01 move f₀ −5.462608e-09 / Q −1.991402e-03 / R_in −2.906654e-03 beside ours' −5.461672e-09 / −1.991409e-03 / −2.906656e-03: equal to five significant digits (−0.199 %), not to the seventh printed digit as predicted. The last-digit gap (7e-9 absolute on Q) is the size of the two per-rung matched-order Q ties (×1 −4.9e-9, item 421; h +1.8e-9, here), so it is read as tie arithmetic, not a matched-order disagreement — (T6′) is 2.5e-10, ≪ 1e-6 — and no known-issues entry is opened; the review may rule otherwise. (C2-h) 3.198e6× on S11 (8.102e-4 / 2.533e-10; predicted ≥ 810×, central ≈ 3e6×) and 1.091e6× on Q (predicted ≥ 2 000×). Controls: (N-feed) Q′ 280.652, +70.90 % (predicted [+20 %, +80 %]: met); (N-nofeed) no crossing (met). W8 (`20261014T094100Z_OPS-61-item444-W8-regress-read.log`): output section diff-empty (16 lines) against item 421 W4's.

*Drafted by a general-purpose agent for the 2026-10-14 03:00 daily review; the review reviews it.* Licence: §9 rule (6) (PROJECT_PLAN.md, grep "**Rule (6), the chain-extension licence**"): the next numbered step of live §10 chain T toward its stated subgoal, Q's h-statement (T41, the 10-17 weekly), tagged *(daily, 2026-10-14 — provisional)*; the weekly ratifies or reverts. It is chain T's second provisional step since the 10-12 weekly (T42 the first). **T43 (provisional):** items 420 + 421 on T37's refined rung — Palace's twin and driven read of the r = 0.2 Ω tune at 64.000 MHz on the h = 0.01 exact-rectangle hole mesh, against T37's Route A rows (item 395 W1). Target: **T40's, unchanged** — f₀ within `T21P_ZERO_RESOLUTION`, Q within `T40_Q_TARGET` (5 %) — with the matched-order ties printed beside the ≤ 1e-6 precedent (the 2026-10-13 rulings (1) / (4): no cross-solver band exists). The capacitor-ESR ruling (i)–(iii), the sheets-and-κ ruling and the Palace directive (§5.4 (i)–(v)) apply. **No band, tolerance, Done-when or validation target is set or moved.**

**Why h = 0.01, not h = 0.005 (the rung whose Q T30 calls unconverged).**
- **No tuned node exists at h = 0.005.** `logs/ans7_exact/` holds only `node_{62.0,63.2,64.5,65.8,67.0}_p2_h0.005_cu.json` (`xl` `ANS-7-step4p`, 523 880 cells, 3 467 542 degree-2 unknowns, xl-ledger row "ANS-7-step4p"): lossless, untuned, at `C_F`. There is no `tune_p2_h0.005*.json` and no C_in at that rung (T29 inferred it; item 395's scope line). The only refined tuned rung is h = 0.01: T28 (lossless, item 310) and T37 (r = 0.2 Ω, item 395) at 64 MHz, T34 / T38 at 127.74 MHz.
- **T30's h = 0.005 quantity is a Route B row** (`sparse_rows`: `bary_interp` of the full 4×4 Z over the five nodes, then `route_b_z` at `C_F`; `tests/validation/test_ans7_sparse_route_b.py:321–338`). Twinning it needs Palace's 4×4 at all five nodes: 20 single-excitation windows on a ≈ 524 k-cell mesh. Scaled from item 437 (304–334 s per excitation at 10 MHz on 3 833 208 unknowns; its printed 64 MHz projection 392–429 s) by 3 467 542 / 3 833 208 = 0.905: **≈ 330–480 s per excitation, ≈ 1.8–2.7 h in all.** Each command fits heavy; the item does not fit a slot. step4p's three-node control missed the five-node rows (`:4140`), so dropping nodes changes the quantity. The xl runner has no Palace route (no `palace` in `scripts/automation/`). **That twin is the weekly's to tier** (five per-node heavy items, or a Palace route on the xl runner, which does not exist yet); this item does not take it.
- **What h = 0.01 buys T41.** T30's first successive change (×1 → 0.01) and T37's Q(r) move (−0.199 %) both sit on this mesh. If Palace at matched order on the same rung reproduces T37's rows to ≈ 1e-9 and its own ×1 → 0.01 move equals ours to the printed digit, then T41 can write "unconverged at h = 0.005" as **the discretisation's, seen by two codes on the first refined rung** — not a defect of ours. One rung does not write an h-statement.

**What exists on our side (the comparand; not re-solved).** T37, `20261012T094711Z_ANS-7-item395-W1-esr64-h001.log` (`-n 12`, 140 270 cells; the latest window that wrote `tune_p2_h0.01_esr0.2.json`, 2026-10-12 04:53 local):
- `:486` — C_in(r, h) 6.364610922203437e-11 F (= `tune_p2_h0.01_esr0.2.json` `c_in_r`; `c_in_0` 6.364316664701229e-11 F = `tune_p2_h0.01.json` `c_in`; `resolution` 0.01);
- `:2203` — rung mesh cells 140 270 (node cells and `tune_p2_h0.01.json` `cells` equal);
- `:2204–2206` — S11_A at (r, C_in(r, h)), 12 digits: 63.9 MHz −2.445046675335e-01+9.625434769724e-01j; 64.0 MHz −2.400384573551e-01+9.635701508391e-01j; 64.1 MHz −2.355297117830e-01+9.645820886247e-01j;
- `:2207` — Im Y_in,A −1.617828160e-04 / +2.684681835e-07 / +1.626210783e-04 S; rows A: f₀ 63 999 834.331344604 Hz, R_in 3 166.7479872416206 Ω, Q 164.21575983485775;
- `:2210–2211` — Q_0(h) 964.432729560342, Q_C(h) 195.3612; Q(r, h)/Q(r, ×1) − 1 = −1.991409e-03.
- Lossless C path on the rung at C_in(0, h), 64.0 MHz: S11_A −2.411870945081e-01+9.683792879359e-01j (`20261008T111904Z_ANS-7-item310-W3-p64.0.log:2288`; C_in(h) 6.364316664701229e-11 F at `:490`).
- Degree-2 unknowns on the rung: 968 016 (`20261004T112945Z_ANS-7-item214-W3-rung.log:2236`).
- Palace's ×1 rows for (C-h): f₀,Pal 63 999 834.68507022 Hz, Q_Pal 164.54343232572575, R_in,Pal 3 175.979463819766 Ω (`20261013T141248Z_OPS-61-item421-W4-read.log:48–50`); ours' ×1: f₀,A 63 999 834.680890724 Hz, Q_A 164.5434331272106, R_in,A 3 175.9794661876795 Ω (`20261012T094112Z_ANS-7-item395-W0-esr64.log:2130`).

`OPS-61` step 3d: export the h = 0.01 exact-rectangle hole mesh, then Palace order 2 on it at 63.9 / 64.0 / 64.1 MHz with ports 2–4 carrying each frequency's exact parallel map of (r, C_in(r, h)), plus the lossless C path at 64.0 MHz; Palace's f₀, R_in and Q through `parallel_rows`, against T37's Route A rows.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Export, configs and read at `-n 1` in the complex container; Palace at `-np 8`. **Class:** chain step (solve), ≈ 18 measured.
- **Time:** ≈ 45 predicted min: 15 fixed + ≈ 15 code + W1 ≈ 1.5 (45–120 s) + W2 < 1 + W3–W6 4 × ≈ 1.5 (80 s central, 65–150 s each) + W7 / W8 < 1 + ≈ 5 spare.
- **Price (scaled from measurement).** Items 420 / 421 measured 51.1–53.2 s Palace Total per point on the ×1 hole mesh (570 332 unknowns, Setup ≈ 8.7 s, 33–35 GMRES iterations, peak 4.2–4.3 G). The rung has 968 016 / 570 332 = 1.697× the unknowns. Item 418's ×1 → 0.0075 step (2.42× unknowns) cost 1.68× per excitation. So per point ≈ 71 s (that sublinear rate) to 88 s (linear) to 102 s (unknowns^1.25), central ≈ 80 s; peak ≈ 6–9 G. Every bracket is ≪ 590 s container and ≪ 1 200 s heavy. The rung mesh built in 32.8 s inside T37's W1 (gmsh is serial on rank 0 at any width); the ×1 hole export (item 234 W1) built in 22.78 s.
- **Depends on:** nothing. Inputs on disk (checked 2026-10-14): `logs/ops61/t39_esr.json`, `t39_lossless.json` (templates, never overwritten), `postpro_t39_esr/`, `postpro_t40_esr_f63.9/`, `postpro_t40_esr_f64.1/`; `logs/ans7_exact/tune_p2_h0.01.json`, `tune_p2_h0.01_esr0.2.json`. Absent, as required: any `ops61_ans6_h0010_hole_strip.msh`, `t43_*.json`, `postpro_t43_*`.

**ESR map values (drafter's arithmetic, for the predictions; the probe recomputes).** At C = C_in(r, h) = 6.364610922203437e-11 F, r = 0.2 Ω: Q_C(63.9 / 64.0 / 64.1 MHz) = 195.66692417 / 195.36119460 / 195.05641894; R_p = 7 657.309043 / 7 633.399271 / 7 609.601314 Ω; C_p = 6.36444469e-11 / 6.36444417e-11 / 6.36444364e-11 F (C_p/C − 1 ≈ −2.62e-5). Palace's sheet values (w 7 mm, l 8 mm): Rs ≈ 6.700e+03 / 6.679e+03 / 6.658e+03 Ω/sq, Cs ≈ 7.274e-11 F/sq; lossless Cs = C_in(0, h)·8/7 ≈ 7.274e-11 F/sq (7.2735e-11). Series \|Z\| at 64.0 MHz 39.073 Ω. Map residual in double precision ≤ 2e-16.

**Code (probe only; additive; disclosed).**
- **Constants** beside `T42_*`, each with its log-line provenance comment: `H001_RES = 0.01` (= `H_RES`, `test_ans7_exact_strip.py:1059`, asserted equal on import); `H001_HOLE_MESH_NAME = "ops61_ans6_h0010_hole_strip.msh"`; `H001_HOLE_CELL_RECORD = 140270` (literal; `…item395-W1…:2203`, `…item214-W3-rung.log:2236`; no tracked constant holds it); `T43_SWEEP = (("t43_esr_f63.9", 0.0639, 63.9e6), ("t43_esr_f64.0", 0.064, 64.0e6), ("t43_esr_f64.1", 0.0641, 64.1e6))`; `T43_LOSSLESS = ("t43_lossless_f64.0", 0.064, 64.0e6)`; `T43_C_IN_R_H_LITERAL = 6.364610922203437e-11`, `T43_C_IN_0_H_LITERAL = 6.364316664701229e-11`; `T43_S11_A` (the three `:2204–2206` values), `T43_F0_A`, `T43_R_IN_A`, `T43_Q_A`, `T43_IM_Y_A` (`:2207`), `T43_Q_H_MOVE_OURS = -1.991409e-03` (`:2211`), `T43_S11_A_LOSSLESS` (`…item310-W3-p64.0.log:2288`), `T43_PAL_X1_ROWS` (item 421 W4 `:48–50`). `T40_Q_TARGET` is reused, not copied.
- **New `export_hole_h001()`, dispatched by `export --hole-h001`.** A **copy** of `export_hole_pec` (`:1095–1319`; provenance comment with commit and line span) that writes **the mesh only, no configs**. `export_hole_pec` is neither edited nor called. Steps:
  1. Assert `comm.size == 1`; assert `H001_HOLE_MESH_NAME` does not exist (never overwrite).
  2. Build with `_sheets_build(True, as_hole=True, leg_sheet_strip_fraction=float(GATED_WIDTH_FRACTION), resolution=H001_RES)`: the call `_exact_geo(th14, comm, GATED_WIDTH_FRACTION, res=H_RES)` makes inside `tune_h` (`test_ans7_exact_strip.py:409–425`, `:2091`). Assert `th14._sheets_build is _sheets_build`.
  3. **Assert ncells == `H001_HOLE_CELL_RECORD` before anything is written.**
  4. Write the Gmsh 2.2 file as `export_hole_pec` does (box → 1, wall 401 → its attribute, N_i → 20 + i).
  5. Anchors after writing: (E1)–(E5) as `export_hole_pec`'s, with `H001_HOLE_CELL_RECORD`, plus (E6) below.
- **New `ans7_h001_configs()`, dispatched by `ans7-h001-configs`.**
  1. Assert `comm.size == 1`; assert none of the four `t43_*.json` exists.
  2. Assert `tune_p2_h0.01_esr0.2.json` `c_in_r` == `T43_C_IN_R_H_LITERAL` and `c_in_0` == `T43_C_IN_0_H_LITERAL` bit for bit, `resolution` == 0.01, `esr_ohm` == 0.2; assert `tune_p2_h0.01.json` `c_in` == `T43_C_IN_0_H_LITERAL` and `cells` == `H001_HOLE_CELL_RECORD`.
  3. From `t39_esr.json` (checked as item 421 checks it), write the three `t43_esr_*` configs changing **only** `Model.Mesh` → `H001_HOLE_MESH_NAME`, `Freq`, ports 2–4 `R` / `C` (`_t39_map(ω_f, T43_C_IN_R_H_LITERAL, 0.2)` at each frequency) and `Problem.Output`. From `t39_lossless.json`, write `t43_lossless_f64.0` changing only `Model.Mesh`, ports 2–4 `C` → C_in(0, h) and `Problem.Output`.
  4. Assert (Z′) / (Z0) and the template diff (item 421's pattern), reading each JSON back.
- **New `read_ans7_h001()`, dispatched by `read --ans7-h001`.** Item 421's `read_ans7_sweep` copied with the T43 sets and literals: S[1][1] only; Z_bare = `_ANS7._z_from_s11`; Z_in = 1/(1/Z_bare + 1/`_z_branch`(ω, C_in(r, h), 0.2)) (the feed's own (r + C) branch in post-processing, ruling 311-Q5); rows = `parallel_rows([63.9e6, 64.0e6, 64.1e6], Z_in)`. Plus (C-h) and (C2-h) below. **`read_ans7_sweep` and `read_ans7_128` are not edited.** It never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** three new `sys.argv` branches, first in `__main__` as items 421 / 436 placed theirs. Every existing mode keeps its behaviour (W8).

**Anchors (asserted; bands imported or quoted from the target, never chosen).**
- **(E1), W1:** tets == 140 270, asserted before the write and again parsed from the written file's `$Elements` type-4 count.
- **(E2)–(E5), W1:** as `export_hole_pec`'s: no conductor cell tag and the phantom tag present; the census identities (`cavity_exterior == cavity > 0`, `outer == exterior − cavity`, `tag_collision == 0`); box area == `EXTERIOR_AREA_M2` within `EXTERIOR_AREA_BAND`; per port, N_i area / width / height == `SPEC_AREA_M2` / `SPEC_WIDTH_M` / `SPEC_HEIGHT_M` within `STRIP_BAND`. Exact rectangles at the rung: the sheets-and-κ ruling's precondition.
- **(E6), W1:** the global cell-tag set equals the ×1 hole mesh's set as `t39_esr.json`'s `Domains.Materials` lists it (phantom `[3]` plus the air list). The templates' materials then carry over unchanged. A rung that renumbered a region fails here, not in Palace.
- **(Z′), W2:** per ESR config, \|Z_par − Z_ser(f)\| / \|Z_ser(f)\| ≤ `STRIP_BAND` (1e-9, `tests/validation/test_ops61_exact_strip_ports.py:58`), Z_ser(f) = `series_rlc_impedance(f, r_ohm=0.2, c_f=T43_C_IN_R_H_LITERAL)`; and \|dω\|/ω ≤ `STRIP_BAND`.
- **(Z0), W2:** the lossless config's `C` == C_in(0, h) bit for bit, no `R` key.
- **(C1), W7:** Palace's Im Y_in changes sign − → + exactly once on the triple and never + → −; f₀,Pal inside (63.9, 64.1) MHz.
- **(f₀), W7:** \|f₀,Pal / f₀,A − 1\| ≤ `T21P_ZERO_RESOLUTION` (1.287e-5, imported), f₀,A = 63 999 834.331344604 Hz. This is T40's band; it is the band T37 printed its offset against at 64 MHz.
- **(Q), W7:** \|Q_Pal / Q_A − 1\| ≤ `T40_Q_TARGET` (0.05, T40's quoted target, imported from the probe), Q_A = 164.21575983485775.
- **Palace log check (executor, journalled), W3–W6:** 140 270 elements; box ±0.12 / ±0.12 / ±0.10 m; `401: σ = 5.800e+07 S/m`, Dirichlet PEC on `1` only; ω/2π = 6.390e-02 / 6.400e-02 / 6.410e-02 / 6.400e-02 GHz; one excitation (port 1 on 21), port 1 `Rs = 4.375e+01 Ω/sq`; ports 2–4 `Rs` / `Cs` quoted verbatim against the values above to the printed digits; one `GMRES solver converged` line; `ND (p = 2): 968016` predicted (= ours; a mismatch is a reading, not a red).

**Printed, never gated (rulings (1) / (4), 2026-10-13).**
- **(T6′)** per-point \|S11_Pal − S11_A\|/\|S11_A\| at 63.9 / 64.0 / 64.1 MHz. Predicted ≤ 1e-6; items 420 / 421 / 436 measured 2.26–2.58e-10 at ×1.
- **(T6-0″)** the lossless C path at 64.0 MHz against `T43_S11_A_LOSSLESS`. Predicted ≤ 1e-6, whatever the R convention.
- R_in,Pal against 3 166.7479872416206 Ω; f₀ and Q gaps in e-notation; Palace's Im Y_in triple beside ours; f₀,Pal/64 MHz − 1 beside ours' −2.588573e-06.
- **(C-h) Palace's own ×1 → 0.01 move**, from Palace's rows here against `T43_PAL_X1_ROWS`, beside ours'. Ours' moves are Q −1.991409e-03, R_in −2.906656e-03, f₀ −5.4617e-09 (drafter's arithmetic from the cited rows; the read recomputes from literals). **Predicted equal to ours to the printed digit**, as items 434 / 435 / 437's (C-h) measured on the exact-strip rungs. This is the line T41 needs.
- Palace Total, Setup, iterations and peak per run.

**Negative controls (rule (e); predicted, printed).**
- **(C2-h) mesh-mismatch separation.** At 64.0 MHz, \|S11_A(ours ×1, `…item395-W0-esr64.log:2128`) − S11_Pal(h)\| / \|S11_A(ours h) − S11_Pal(h)\|. The numerator ≈ the ×1 → 0.01 S11 move, 8.10e-4 relative (drafter's arithmetic from the two 12-digit lines). The denominator is (T6′) ≈ 2.6e-10. **Predicted ≥ 810×** against item 229's ≥ 25× bar; ceiling ≈ 8.1e-4 / 1e-12 (the 12-digit print floor) ≈ 8e8×, central ≈ 3e6×. On Q the same separation is \|Q_A(×1) − Q_Pal(h)\| / \|Q_A(h) − Q_Pal(h)\|, numerator 1.99e-3 relative; predicted ≥ 2 000×. A Palace run that silently read the ×1 mesh reads ≈ 1 here and fails (E1) / the log check first.
- **(K-esr″) r reaches the operator:** \|S11_Pal(esr, 64.0) − S11_Pal(lossless, 64.0)\| / \|S11_Pal(lossless)\| beside ours' same move, 4.954501e-03 (drafter's arithmetic from `:2205` and `…item310-W3…:2288`; ×1 4.949105e-03). Predicted equal to ≤ 1e-6 relative. A Palace that dropped `R` reads ≈ 1e-10 here.
- **(N-feed) the feed's r omitted** (post-processing on the same three S11): Q′ from Z_in′ = 1/(1/Z_bare + jωC_in(r, h)). With 1/Q′ = 1/Q_0 + (1 − s_feed)/Q_C, Q_0(h) 964.43 and Q_C(h) 195.36, **predicted Q′/Q_A − 1 in [+20 %, +80 %]** (×1 measured +71.06 %, item 421), ≥ 4× the 5 % target; ceiling Q′ ≤ Q_0(h), i.e. ≤ +487 %.
- **(N-map) the 64.0 MHz map reused at 63.9 / 64.1 MHz**, pure arithmetic: 1.6008e-5 / 1.5983e-5 (drafter's arithmetic). It fails `STRIP_BAND` by ≈ 1.6e4×, against (Z′) ≤ 2e-16. Not run as Palace windows.
- **(N-nofeed) the feed branch dropped** (Z_in″ = Z_bare): predicted **no** − → + crossing on the triple, "f₀ undefined" (×1 and 127.74 MHz: met).

**Predictions (written first, printed).** W1 45–120 s, (E1) 140 270, (E4) to ≈ 1e-16, (E5) ≤ 2e-15. W2 < 5 s, (Z′) ≤ 2e-16. W3–W6 each Palace Total ≈ 80 s (65–150), GMRES 33–45 iterations (`MaxIts` 100), peak 6–9 G. (T6′) / (T6-0″) ≤ 1e-6; f₀ / Q / R_in gaps ≤ 1e-6; (C-h) equal to ours' moves to the printed digit; (C2-h) ≥ 810× on S11 and ≥ 2 000× on Q; (K-esr″) ratio 1 ± 1e-6; (N-feed) in [+20 %, +80 %]; (N-nofeed) no crossing.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call): `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up; `ls logs/ops61/ logs/ans7_exact/`. Stop and report if `ops61_ans6_h0010_hole_strip.msh`, any `t43_*.json` or any `postpro_t43_*` exists. The templates, the three ×1 `postpro_t39/t40_*` outputs and both rung tune JSONs must be present. Never re-create one.
- **W1** (45–120 s), the export, (E1)–(E6):

  `scripts/testing/run_and_log.sh OPS-61-item444-W1-export-hole-h001 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py export --hole-h001'"`

  If the private-leak hook refuses the W1 log on a printed byte size or node count (items 434 / 435), quarantine it unedited under the gitignored `logs/ops61_item444_quarantine/` as they did. The test-results row stays, and (E1) is corroborated by every W3–W6 Palace log.
- **W2** (< 10 s), configs, (Z′), (Z0):

  `scripts/testing/run_and_log.sh OPS-61-item444-W2-ans7-h001-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py ans7-h001-configs'"`
- **W3** (≈ 80 s), the tuned point first, the price check:

  `scripts/testing/run_and_log.sh OPS-61-item444-W3-palace-t43-f64.0 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 t43_esr_f64.0.json'"`
- **Price gate (binding, after W3).** Run W4–W6 iff W3's Palace Total (Avg) ≤ 300 s (≤ 4 × 300 s fits the slot with the read). Otherwise the item lands as a partial record, and the review re-prices.
- **W4 / W5 / W6:** W3's command with `t43_esr_f63.9.json` / `t43_esr_f64.1.json` / `t43_lossless_f64.0.json`, labelled `OPS-61-item444-W4-palace-t43-f63.9` / `OPS-61-item444-W5-palace-t43-f64.1` / `OPS-61-item444-W6-palace-t43-lossless`. After any kill, run `docker compose -f docker/docker-compose.yml ps` before the next window.
- **W7** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item444-W7-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --ans7-h001'"`
- **W8** (≈ 3 s), the dispatch regression: item 421's W4 command verbatim (`20261013T141248Z_OPS-61-item421-W4-read.log:12`), labelled `OPS-61-item444-W8-regress-read`. Its output section must `diff` empty against that log's.

**Traps already paid for.** All of items 420 / 421 / 436's apply, plus:
- **Never reuse a ×1 value on the rung:** C_in(r, h) is 6.364610922203437e-11 F, not `T26_C_IN_R_X1` 6.359301781809668e-11 F (item 419's "silent misread": a ×1 literal under a rung print fails nothing). Grep the W2 log for `6.3593` / `6.3590`: there must be none.
- **The mesh name is the only mesh link.** A config that keeps `Model.Mesh == ops61_ans6_x1_hole_strip.msh` solves the ×1 problem at the rung's C and reads (T6′) ≈ 8e-4. The template diff asserts `Model.Mesh` changed, and the Palace log check's 140 270 catches it first.
- **Every frequency gets its own exact map** (ruling (iii), R_p ∝ 1/f² off the conversion frequency); never one Palace sweep with fixed R_p / C_p. Never enter 0.2 Ω as a parallel R (no short control is re-run: item 420 read the convention, parallel).
- **The feed's (r + C) branch is in the read, not in Palace** (ruling 311-Q5). Dropping its r is (N-feed), the control, never the reading.
- **The export runs at `-n 1`** (the function asserts it). Never run `export`, `export --exact-strip`, `export --hole-pec`, `export --exact-strip-h0075` / `-h0050` or any older `*-configs` mode. Never overwrite `t39_*` / `t40_*` / `t42_*` or their outputs.
- **Compare S[1][1] only**; never read `z_matrix_ohm` or the `aed` key.
- Docker socket and orphans: as item 420. After any killed window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0, and `docker compose -f docker/docker-compose.yml ps` must show no lingering `fem-em-palace`. Never lengthen a timeout or raise `-np`. A `GMRES … did not converge` line is a finding: that point is "not measured", and with a tune point missing (C1) cannot be read.
- Record row and journal: Palace's numbers in full (v0.18.1, the four configs, the mesh) beside ours; **no Palace figure set against an AED one in a tracked file**. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace's bare S11, f₀, R_in and Q at the 64 MHz ESR-tuned F-small node on the h = 0.01 exact-rectangle hole mesh, degree 2, Cu, read as T37 reads them; the lossless C-path tie on the rung; Palace's own ×1 → 0.01 move beside ours'. Not covered:
- an h-statement for Q (T41's, from F31 + T30);
- the h = 0.005 rung, whose Route B twin at `C_F` is the weekly's to tier (above);
- T30's `C_F` rows at h = 0.01 (lossless, untuned; a Route B twin of item 214's five nodes, 20 excitations);
- 127.74 MHz on the rung (T38's twin);
- order 1; AED; a measured ESR;
- an absolute tuned C, f₀, R_in or Q claim beyond "Palace = ours at matched order on one refined rung";
- `src/` and any band.

**Status it can move.**
- The `ANS-7` row's Q clause gains "T43 (item 444, record, provisional): Palace at matched order on the h = 0.01 rung reads f₀ within `T21P_ZERO_RESOLUTION` and Q within T40's 5 % of T37's rows, and its own ×1 → 0.01 Q move equals ours (−0.199 %)". That is a record for T41. **No glyph:** the Done-when is Q's h-statement (T41). Rule (j): read; no mismatch.
- The `OPS-61` row (✅) gains a "step 3d record (T43)" clause.
- §10 T43 recorded DONE (provisional; the weekly ratifies or reverts).

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(E1) red:** the export is not T37's rung mesh. Stop before any Palace window, park on `attempt/*`, BLOCKED (rule (d)). Open a known-issues entry only if the count moved on an unchanged image; a re-record is `record-reconciler`'s, under (1*).
- **(E2)–(E6) red:** the strip or the region tags did not reach the rung. Report the spans, facet counts and tag sets; park, BLOCKED.
- **W1 Status 124:** a price finding; record it, park, and return the export's tier to the review. Never lengthen the timeout.
- **(Z′) / (Z0) red, or a W2 bit-equality assert red:** the wrong C or map. Stop before any Palace window, revert the probe, BLOCKED.
- **A W3–W6 log lacking 140 270 elements, σ = 5.800e+07 on 401, PEC on 1 only, or the right ω/2π:** park, BLOCKED.
- **(C1) red:** report Palace's Im Y_in triple; park, BLOCKED; the review reads it against (T6′).
- **(f₀) or (Q) red with (C1) green:** a disagreement on the external check on the refined rung. That would be the first one, and it is the reading T41 most needs. Report the gaps, (T6′) and (C-h), open a known-issues entry, park, BLOCKED. Never widen, re-map or change r.
- **(T6′) ≫ 1e-6, or (C-h) ≠ ours' move beyond the print, with (f₀) / (Q) green:** a reading, not a failure. DONE as a record with a known-issues entry: the codes agree on T40's targets but not at matched order on the rung, so T41 cannot cite this item as "the discretisation's".
- **(T6-0″) ≫ 1e-6:** a convention-free disagreement on the rung. DONE as a record with a known-issues entry (item 420's disposition); (f₀) / (Q) are still read.
- **W8 not diff-empty:** revert the probe, BLOCKED.
- **Price gate trips:** a partial record, not BLOCKED.

**Questions for the review / weekly.**
- (Q1) The h = 0.005 twin of T30's `C_F` row (20 excitations, ≈ 1.8–2.7 h of Palace) is the weekly's to tier: five per-node heavy items, or a Palace route on the xl runner, which does not exist yet.
- (Q2) T43 is chain T's second rule (6) step since the 10-12 weekly.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 444 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; the export at `-n 1`; one frequency per Palace window, its own map per frequency, C_in(r, h) = 6.364610922203437e-11 F (never a ×1 literal); the price gate after W3 binds. Probe only — no `src/`, no test module; never run an older `export` or `*-configs` mode; never overwrite `t39_*` / `t40_*` / `t42_*`; never open `docs/private/` or `aed_results/`."

### Review ruling 2026-10-14 03:00 (item 444)

The brief is accepted on **T37's h = 0.01 rung**, not T30's h = 0.005. h = 0.005 has no tuned node (`logs/ans7_exact/` holds only the five untuned, lossless `C_F` nodes). Its Route B twin is 20 single-excitation Palace runs (≈ 1.8–2.7 h), which does not fit a slot and is returned to the 10-14 weekly to tier.

The `export --hole-h001` exporter copy is accepted: the same call `tune_h` makes, with `resolution=0.01`. The configs come from `t39_esr` / `t39_lossless` with only the mesh, frequency, port R/C and output folder changed. The lumped convention is item 420's parallel R ‖ C map per frequency.

### §9 item 448 — 2026-10-14 (chain step T44, daily — provisional)

**Outcome 2026-10-14 (12:00 CDT slot, DONE — record, provisional; no glyph).** Every asserted anchor green on the first run; nothing killed; no probe edit after W1 (last edit before 17:02:23 UTC). The probe gained `T44_*` constants (provenance comments), `ans7_h001_128_configs()` (`ans7-h001-configs --f128`) and `read_ans7_h001_128()` (`read --ans7-h001 --f128`), additive copies of item 444's two functions with every 64 MHz comparand replaced (lossless `Freq` → [0.12774], `T44_PAL_X1_ROWS` from item 436 W6, `T42_*` / `T42_S11_A[1]` for the ×1 side); dispatch branches first in `__main__`; no export run, the mesh reused. W1 (`20261014T170223Z_OPS-61-item448-W1-ans7-h001-128-configs.log`, 4 s): (E1′) `$Elements` type-4 count 140 270 before any write; C_in(r, h) 1.5288749690061608e-11 F and C_in(0, h) 1.5288541964455146e-11 F bit-equal to `t127.74/` JSONs (`f_hz` 127740000.0); (Z′) 1.742e-16 / 0 / 1.363e-18; (Z0) bit-equal, `Freq` [0.12774]; (N-map) 3.844e-6 / 3.841e-6; no `1.52750` / `6.3646` / `6.3593` / `6.3590` in the log. W2–W5 Palace v0.18.1 order 2 at `-np 8` (`…170232Z…W2-palace-t44-f127.74.log`, `…170434Z…W3-palace-t44-f127.64.log`, `…170629Z…W4-palace-t44-f127.84.log`, `…170820Z…W5-palace-t44-lossless.log`): 140 270 elements, box ±0.12 / ±0.12 / ±0.10 m, `401: σ = 5.800e+07`, PEC on `1` only, port 1 `Rs = 4.375e+01`, ports 2–4 `Rs` 2.906e+04 / 2.910e+04 / 2.901e+04 Ω/sq, `Cs` 1.747e-11 F/sq (lossless `Cs` only), ω/2π 0.1277 / 0.1276 / 0.1278 / 0.1277 GHz, ND 968 016 (= ours), GMRES 50 / 50 / 50 / 52 iterations, Palace Total 111.206 / 110.272 / 106.551 / 112.285 s, peak 7.6–7.7 G; the price gate after W2 (111.2 s ≤ 300 s) passed. W6 (`20261014T171034Z_OPS-61-item448-W6-read.log`): (C1) one − → + crossing at 127 739 705.27084099 Hz; **(f₀) +5.912e-11** ≤ `T21P_ZERO_RESOLUTION` (within `T22P_ZERO_RESOLUTION` too), f₀,Pal/127.74 MHz − 1 −2.307258e-06 vs ours −2.307317e-06; **(Q) Q_Pal 252.30120807645983, +1.027e-08** ≤ `T40_Q_TARGET`; R_in,Pal 9 698.113367618138 Ω, +4.872e-09; (T6′) 2.408e-10 / 2.267e-10 / 2.128e-10; (T6-0‴) 2.093e-10; (K-esr‴) 2.185645e-03 both, ratio gap 4.356e-09. **(C-h)** Palace's ×1 → 0.01 move f₀ −5.514397e-09 / Q −8.519899e-03 / R_in −9.456957e-03 beside ours' −5.512787e-09 / −8.519904e-03 / −9.456959e-03, \|ΔQ-move\| 4.786e-09 (predicted ≲ 1e-8: met). (C2-h) 3.538e6× on S11, 8.366e5× on Q; (N-feed) +41.26 % (met; predicted central +41.3 %); (N-nofeed) no crossing (met). W7 / W8 (`…171041Z…W7-regress-read-h001.log`, `…171044Z…W8-regress-read-128.log`) diff-empty against item 444 W7 and item 436 W6. Every prediction met, including the iteration bracket (45–60) and the Total bracket (90–150 s).

*Drafted by a general-purpose agent for the 2026-10-14 10:30 daily review; the review reviews it.* Licence: §9 rule (6) (PROJECT_PLAN.md, grep "**Rule (6), the chain-extension licence**", `:3835–3845`): the next numbered step of live §10 chain T toward its stated subgoal, Q's h-statement (T41, the 10-17 weekly), tagged *(daily, 2026-10-14 — provisional)*; the weekly ratifies or reverts. It is chain T's **third** provisional step since the 10-12 weekly (T42 = item 436, T43 = item 444), so it is the **last** rule (6) allows on chain T before the next weekly. **T44 (provisional):** T43's recipe at the second Larmor frequency, as T42 did items 420 + 421's ×1 recipe at 127.74 MHz — Palace's twin and driven read of T38's r = 0.2 Ω tune at **127.74 MHz** on the **same** h = 0.01 exact-rectangle hole mesh item 444 exported, against T38's Route A rows. Target: **T40's, unchanged** — f₀ within `T21P_ZERO_RESOLUTION`, Q within `T40_Q_TARGET` (5 %) — with Palace's ×1 → 0.01 move printed beside ours. The capacitor-ESR ruling (i)–(iii), the sheets-and-κ ruling, the Palace directive (§5.4 (i)–(v)) and the 2026-10-13 rulings (1) / (3) / (4) (no cross-solver band; at 127.74 MHz T21P asserted, T22P printed beside) apply. **No band, tolerance, Done-when or validation target is set or moved.**

**Why this step.** T43 showed both codes move Q by the same −0.199 % from ×1 to h = 0.01 at 64 MHz (`20261014T094041Z_OPS-61-item444-W7-read.log:54`). At 127.74 MHz ours moves Q four times as far, −8.519904e-03 (`20261013T125603Z_ANS-7-item419-W1-esr128-h001.log:2211`), and the Q(r)-move / Q₀-move ratio differs between the frequencies (0.4525 at 127.74 MHz against 0.2915 at 64 MHz, `:2211`). T41 cannot call the 127.74 MHz move "the discretisation's, seen by two codes" without the same twin at that frequency. One rung still writes no h-statement.

**The mesh is reused as-is (checked 2026-10-14, read-only).** `logs/ops61/ops61_ans6_h0010_hole_strip.msh` exists (written 2026-10-14 04:33 local by item 444 W1, 140 270 tets, (E1)–(E6) green, `20261014T093258Z_OPS-61-item444-W1-export-hole-h001.log`). The rung mesh does not depend on frequency. `_exact_geo(th14, comm, strip, res)` takes no frequency argument (`tests/validation/test_ans7_exact_strip.py:409–425`). T38's rung mesh is the same call, `_exact_geo(th14, comm, GATED_WIDTH_FRACTION, res=H_RES)` (`:2434`), and it printed 140 270 cells (`…item419-W1-esr128-h001.log:2203`), the count item 444's export asserted. So **no re-export**: item 448 never runs `export --hole-h001`, and that function asserts the file absent anyway. The configs mode re-reads the file's `$Elements` type-4 count as its own (E1′) before writing anything (below).

**What exists on our side (the comparand; not re-solved).** T38 = `ANS-7` step 14, §9 item 419. Its W1, `20261013T125603Z_ANS-7-item419-W1-esr128-h001.log` (`-n 12`, Status 0, 378 s at `:2278–2279`), is the window that wrote `logs/ans7_exact/t127.74/tune_p2_h0.01_esr0.2.json` (mtime 2026-10-13 08:02 local = 13:02Z, the W1 window). The log holds every row this item needs:
- `:476` — C_in(r = 0.2 Ω, h 0.01) **1.5288749690061608e-11 F**, secant from C_in(0, h) **1.5288541964455146e-11 F**; both bit-equal to the JSON's `c_in_r` / `c_in_0` (`resolution` 0.01, `esr_ohm` 0.2, `f_hz` 127740000.0). `t127.74/tune_p2_h0.01.json`: `c_in` 1.5288541964455146e-11, `cells` 140270, `resolution` 0.01.
- `:2203` — rung mesh cells 140 270 (node cells and `tune_p2_h0.01.json` 140 270).
- `:2204–2206` — S11_A at (r, C_in(r, h)), 12 digits: 127.64 MHz 4.486361641616e-01+8.878179861686e-01j; 127.74 MHz 4.506307628802e-01+8.867812811584e-01j; 127.84 MHz 4.526280544402e-01+8.857371566245e-01j.
- `:2207` — Im Y_in,A −4.061198730e-05 / +1.200522737e-07 / +4.089005599e-05 S. Rows A: f₀ **127 739 705.26328915 Hz**, R_in **9 698.113320373346 Ω**, Q **252.30120548497052**.
- `:2208` — one − → + crossing; f₀,A/127.74 MHz − 1 = −2.307317e-06, printed beside `T22P_ZERO_RESOLUTION` 1.616e-05. No zero band is registered for the rung (rulings 310-Q1 / 365-Q1 / 395-Q3), so T38 asserted only the crossing.
- `:2210` — Q_0(h) 606.6969050131794 (`t127.74/tune_inmodel_cin_p2_h0.01.json`, item 365b A); Q_C(h) 407.4659.
- `:2211` — Q(r, h)/Q(r, ×1) − 1 = **−8.519904e-03**, against `T35_Q_R_X1` 254.46925911792314. `:2212` — R_in(r, h)/R_in(r, ×1) − 1 = **−9.456959e-03**, against `T35_R_IN_R_X1` 9 790.703605847839.
- Lossless C path on the rung at C_in(0, h), 127.74 MHz: S11_A 4.516287177353e-01+8.887181496060e-01j (`20261011T111940Z_ANS-7-item365b-D2-p127.74.log:2288`; C_in(h) 1.5288541964455146e-11 F read at `:493`; cells 140 270 on `:2288`).
- Degree-2 unknowns on the rung: 968 016. Ours: `20261004T112945Z_ANS-7-item214-W3-rung.log:2236`. Palace: `20261014T093349Z_OPS-61-item444-W3-palace-t43-f64.0.log:104`.
- **The ×1 rows for (C-h).** Palace ×1 at 127.74 MHz: f₀ 127 739 705.97524849 Hz, R_in 9 790.703634300547 Ω, Q 254.46926048070927 (`20261013T215115Z_OPS-61-item436-W6-read.log:48`). Ours ×1: `T42_F0_A` 127 739 705.96749091 Hz, `T42_R_IN_A` 9 790.70360584452 Ω, `T42_Q_A` 254.46925909529418 (probe `:2863–2865`, from `20261013T124951Z_ANS-7-item419-W0-esr128.log:2130`). `T42_S11_A[1]` (×1, 127.74 MHz) 4.513684481906e-01+8.864778722343e-01j (`…item419-W0…:2128`, probe `:2859–2861`).
- **Ours' ×1 → 0.01 move at 127.74 MHz (drafter's arithmetic from the cited rows; the read recomputes from literals):** f₀ −5.512787e-09, Q −8.519904e-03, R_in −9.456959e-03. Against `T42_*` (item 419 W0) these agree with T38's printed `:2211–2212` (against item 387's `T35_*`) to all seven printed digits.

`OPS-61` step 3e: Palace order 2 on item 444's h = 0.01 exact-rectangle hole mesh at 127.64 / 127.74 / 127.84 MHz. Ports 2–4 carry each frequency's exact parallel map of (r, C_in(r, h)), plus the lossless C path at 127.74 MHz. Palace's f₀, R_in and Q come through `parallel_rows`, against T38's Route A rows.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Configs and the read at `-n 1` in the complex container; Palace at `-np 8`. **Class:** chain step (solve), ≈ 10 measured.
- **Time:** ≈ 40 predicted min. That is 15 fixed + ≈ 12 code (two copies with constants substituted) + W1 < 1 + W2–W5 4 × ≈ 1.8 (≈ 100 s central, 90–150 s each) + W6–W8 < 1 + ≈ 5 spare. No export window.
- **Price (scaled from measurement).** On the rung at 64 MHz item 444 measured Palace Total 94.502 s with Setup 16.631 s, 40 GMRES iterations and peak 7.5 G (`…item444-W3-palace-t43-f64.0.log:190`, `:181`, `:159`, `:173`); W4–W6 measured 93.784–96.910 s. On the ×1 mesh, going from 64 to 127.74 MHz moved the per-point Total from 51.1–53.2 s (items 420 / 421) to 55.982 s (Setup 8.106 s, 42 iterations, peak 4.4 G; `20261013T214704Z_OPS-61-item436-W2-palace-t42-f127.74.log:192`, `:183`, `:161`, `:175`). The non-Setup part rose ≈ ×1.10 and the iterations rose 33–35 → 42 (×1.24). So per point on the rung ≈ 16.6 + 77.9 × 1.10 ≈ 102 s (time-scaled) to 16.6 + 77.9 × 1.24 ≈ 113 s (iteration-scaled). Predicted ≈ 100 s central (90–150), GMRES 45–60 iterations, peak 7.5–8.5 G. Every bracket is ≪ the 590 s container timeout and ≪ the 1 200 s heavy ceiling.
- **Depends on:** nothing. Inputs on disk (checked 2026-10-14): `logs/ops61/ops61_ans6_h0010_hole_strip.msh`, `t39_esr.json`, `t39_lossless.json` (templates, never overwritten); `logs/ans7_exact/t127.74/tune_p2_h0.01.json`, `tune_p2_h0.01_esr0.2.json`. Absent, as required: any `t44_*.json`, any `postpro_t44_*`.

**ESR map values (drafter's arithmetic, for the predictions; the probe recomputes).** At C = C_in(r, h) = 1.5288749690061608e-11 F and r = 0.2 Ω:
- Q_C(127.64 / 127.74 / 127.84 MHz) = 407.78511594 / 407.46588538 / 407.14715424. The 127.74 value matches T38's 407.4659 at `:2210`.
- R_p = 33 257.940156 / 33 205.889549 / 33 153.961041 Ω.
- C_p = 1.52886577e-11 / 1.52886576e-11 / 1.52886575e-11 F (C_p/C − 1 ≈ −6.01e-6 to −6.03e-6).

Palace's sheet values (w 7 mm, l 8 mm): Rs ≈ 2.9101e+04 / 2.9055e+04 / 2.9010e+04 Ω/sq, Cs ≈ 1.7473e-11 F/sq. Lossless Cs = C_in(0, h)·8/7 ≈ 1.7473e-11 F/sq. Series \|Z\| at 127.74 MHz 81.4934 Ω. ×1 for scale (item 436): R_p 33 265.483 Ω, Rs 2.911e+04 Ω/sq. The rung's values differ from the ×1 values only in the fourth significant digit. **That is the silent-misread hazard**, so the (Z′) anchor and the bit-equality asserts, not a reader's eye, are what catch a ×1 literal.

**Code (probe only; additive; disclosed).**
- **Constants** beside `T43_*`, each with its log-line provenance comment:
  - `T44_SWEEP = (("t44_esr_f127.64", 0.12764, 127.64e6), ("t44_esr_f127.74", 0.12774, 127.74e6), ("t44_esr_f127.84", 0.12784, 127.84e6))`;
  - `T44_LOSSLESS = ("t44_lossless_f127.74", 0.12774, 127.74e6)`;
  - `T44_C_IN_R_H_LITERAL = 1.5288749690061608e-11`, `T44_C_IN_0_H_LITERAL = 1.5288541964455146e-11` (`…item419-W1…:476`);
  - `T44_S11_A` (the three `:2204–2206` values); `T44_F0_A`, `T44_R_IN_A`, `T44_Q_A`, `T44_IM_Y_A` (`:2207`); `T44_Q_H_MOVE_OURS = -8.519904e-03` (`:2211`); `T44_R_IN_H_MOVE_OURS = -9.456959e-03` (`:2212`); `T44_Q_0_H = 606.6969050131794` (`:2210`, printed ceiling only);
  - `T44_S11_A_LOSSLESS` (`…item365b-D2-p127.74.log:2288`);
  - `T44_PAL_X1_ROWS = (127739705.97524849, 254.46926048070927, 9790.703634300547)` (f₀, Q, R_in; `…item436-W6-read.log:48`).

  `H001_HOLE_MESH_NAME`, `H001_HOLE_CELL_RECORD`, `T42_F0_A` / `T42_Q_A` / `T42_R_IN_A` / `T42_S11_A`, `T40_Q_TARGET` and `T39_ESR_OHM` are **reused, not copied**.
- **New `ans7_h001_128_configs()`, dispatched by `ans7-h001-configs --f128`.** A copy of `ans7_h001_configs` (`:3795–3949`; provenance comment with commit and span) with the T44 sets and literals. `ans7_h001_configs` is neither edited nor called. Steps:
  1. Assert `comm.size == 1`; assert none of the four `t44_*.json` exists.
  2. **(E1′)** Assert `H001_HOLE_MESH_NAME` exists, parse its `$Elements` type-4 count, and assert == `H001_HOLE_CELL_RECORD` (140 270) **before any config is written**. Print the count only, never the file's byte size (the private-leak hook trap, items 434 / 435).
  3. Assert `t127.74/tune_p2_h0.01_esr0.2.json` `c_in_r` == `T44_C_IN_R_H_LITERAL` and `c_in_0` == `T44_C_IN_0_H_LITERAL` bit for bit, `resolution` == 0.01, `esr_ohm` == 0.2, `f_hz` == 127740000.0. Assert `t127.74/tune_p2_h0.01.json` `c_in` == `T44_C_IN_0_H_LITERAL`, `cells` == `H001_HOLE_CELL_RECORD` and `f_hz` == 127740000.0. **The paths are under `t127.74/`**: the top-level `tune_p2_h0.01*.json` are T37's 64 MHz files.
  4. From `t39_esr.json` (checked as item 444 checks it), write the three `t44_esr_*` configs. Change **only** `Model.Mesh` → `H001_HOLE_MESH_NAME`, `Freq`, ports 2–4 `R` / `C` (`_t39_map(ω_f, T44_C_IN_R_H_LITERAL, 0.2)` at each frequency) and `Problem.Output`. From `t39_lossless.json`, write `t44_lossless_f127.74`, changing only `Model.Mesh`, `Freq` → [0.12774], ports 2–4 `C` → C_in(0, h) and `Problem.Output`. (Item 444's lossless kept the template's 64 MHz `Freq`; this one must change it, as item 436's did.)
  5. Assert (Z′) / (Z0) and the template diff (item 444's `strip` pattern, with `Freq` stripped on the lossless config too), reading each JSON back. Print (N-map) at 127.64 / 127.84 MHz from the 127.74 MHz map.
- **New `read_ans7_h001_128()`, dispatched by `read --ans7-h001 --f128`.** A copy of `read_ans7_h001` (`:3952–4109`) with the T44 sets and literals:
  - S[1][1] only; Z_bare = `_ANS7._z_from_s11`; Z_in = 1/(1/Z_bare + 1/`_z_branch`(ω, C_in(r, h), 0.2)), the feed's own (r + C) branch in post-processing (ruling 311-Q5); rows = `parallel_rows([127.64e6, 127.74e6, 127.84e6], Z_in)`.
  - (C1) window (127.64, 127.84) MHz. (f₀) also prints `T22P_ZERO_RESOLUTION` beside T21P, as `read_ans7_128` does (`:3120`), and prints f₀,Pal/127.74 MHz − 1 beside ours' −2.307317e-06.
  - (C-h) against `T44_PAL_X1_ROWS` (Palace) and `T42_*` (ours ×1 → `T44_*`). (C2-h) S11 numerator uses `T42_S11_A[1]`, **not** `T39_S11_A_ESR` (that is the 64 MHz ×1 value `read_ans7_h001` uses); the Q numerator uses `T42_Q_A`.
  - (N-feed) band [+20 %, +70 %] (item 436's), ceiling `T44_Q_0_H`.
  - Imports `T21P_ZERO_RESOLUTION`, `T22P_ZERO_RESOLUTION`, `_ANS7`, `_z_branch`, `parallel_rows`.
  - **`read_ans7_h001`, `read_ans7_128` and `read_ans7_sweep` are not edited.** It never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** two new exact-list `sys.argv` branches, first in `__main__` as items 436 / 444 placed theirs. Every existing mode keeps its behaviour (W7 / W8).

**Anchors (asserted; bands imported or quoted from the target, never chosen).**
- **(E1′), W1:** the reused mesh's `$Elements` type-4 count == 140 270 (`H001_HOLE_CELL_RECORD`; item 444 W1 measured it before and after its write).
- **(Z′), W1:** per ESR config, \|Z_par − Z_ser(f)\| / \|Z_ser(f)\| ≤ `STRIP_BAND` (1e-9, `tests/validation/test_ops61_exact_strip_ports.py:58`), with Z_ser(f) = `series_rlc_impedance(f, r_ohm=0.2, c_f=T44_C_IN_R_H_LITERAL)`; and \|dω\|/ω ≤ `STRIP_BAND`.
- **(Z0), W1:** the lossless config's `C` == C_in(0, h) bit for bit, no `R` key, `Freq` == [0.12774].
- **(C1), W6:** Palace's Im Y_in changes sign − → + exactly once on the triple and never + → −; f₀,Pal inside (127.64, 127.84) MHz.
- **(f₀), W6:** \|f₀,Pal / f₀,A − 1\| ≤ `T21P_ZERO_RESOLUTION` (1.287e-5, imported), f₀,A = 127 739 705.26328915 Hz. This is T40's band as chain T states the target, asserted as item 436 asserted it at 127.74 MHz (ruling 2026-10-13 10:30 (3)). `T22P_ZERO_RESOLUTION` (1.616e-5) is printed beside.
- **(Q), W6:** \|Q_Pal / Q_A − 1\| ≤ `T40_Q_TARGET` (0.05, imported from the probe), Q_A = 252.30120548497052.
- **Palace log check (executor, journalled), W2–W5:** 140 270 elements; box ±0.12 / ±0.12 / ±0.10 m; `401: σ = 5.800e+07 S/m`, Dirichlet PEC on `1` only; one excitation (port 1 on 21), port 1 `Rs = 4.375e+01 Ω/sq`; ω/2π = 1.2774e-01 / 1.2764e-01 / 1.2784e-01 / 1.2774e-01 GHz (Palace prints 0.1277 / 0.1276 / 0.1278 / 0.1277, item 436 W2–W5); ports 2–4 `Rs` / `Cs` quoted verbatim against the values above to the printed digits; one `GMRES solver converged` line; `ND (p = 2): 968016` predicted (= ours and item 444's; a mismatch is a reading, not a red).

**Printed, never gated (rulings (1) / (4), 2026-10-13).**
- **(T6′)** per point, \|S11_Pal − S11_A\|/\|S11_A\| at 127.64 / 127.74 / 127.84 MHz. Predicted ≤ 1e-6. Measured precedents: item 436 2.259–2.373e-10 (×1, 127.74 MHz); item 444 2.507–2.591e-10 (rung, 64 MHz).
- **(T6-0‴)** the lossless C path at 127.74 MHz against `T44_S11_A_LOSSLESS`. Predicted ≤ 1e-6, whatever the R convention.
- R_in,Pal against 9 698.113320373346 Ω; f₀ and Q gaps in e-notation; Palace's Im Y_in triple beside ours.
- **(C-h) Palace's own ×1 → 0.01 move at 127.74 MHz**, from Palace's rows here against `T44_PAL_X1_ROWS`, beside ours' (f₀ −5.512787e-09, Q −8.519904e-03, R_in −9.456959e-03). **Predicted equal to ours to ≈ five significant digits** (item 444 measured −1.991402e-03 against −1.991409e-03 at 64 MHz, `…item444-W7-read.log:54`). The last-digit gap is predicted to be the size of the two per-rung matched-order Q ties: ×1 +5.444e-9 (`…item436-W6-read.log:52`) plus this item's (T6′)-sized tie, so \|ΔQ-move\| ≲ 1e-8. This is the line T41 needs.
- Palace Total, Setup, iterations and peak per run.

**Negative controls (rule (e); all *predicted*, printed — no prior Palace measurement exists at 127.74 MHz on the rung).**
- **(C2-h) mesh-mismatch separation.** At 127.74 MHz, \|S11_A(ours ×1) − S11_Pal(h)\| / \|S11_A(ours h) − S11_Pal(h)\|. The numerator ≈ the ×1 → 0.01 S11 move, 8.0189e-04 relative (drafter's arithmetic, `T42_S11_A[1]` against `:2205`). The denominator is (T6′), ≈ 2.3–2.6e-10. **Predicted ≥ 800×** (if (T6′) ≤ 1e-6), central ≈ 3e6× (item 444 measured 3.198e6× at 64 MHz, `:56`). On Q, \|Q_A(×1) − Q_Pal(h)\| / \|Q_A(h) − Q_Pal(h)\|, numerator 8.5931e-03 relative: **predicted ≥ 8 000×** (item 444 measured 1.091e6× on a 1.995e-3 numerator, `:55`). A Palace run that silently read the ×1 mesh reads ≈ 1 here and fails the log check's 140 270 first.
- **(K-esr‴) r reaches the operator:** \|S11_Pal(esr, 127.74) − S11_Pal(lossless, 127.74)\| / \|S11_Pal(lossless)\| beside ours' same move on the rung, 2.185645e-03 (drafter's arithmetic from `:2205` and `…item365b-D2…:2288`; ×1 2.181790e-03, item 436). Predicted equal to ≤ 1e-6 relative (item 444 measured 1.717e-9 at 64 MHz, `:43`). A Palace that dropped `R` reads ≈ 1e-10 here.
- **(N-feed) the feed's r omitted** (post-processing on the same three S11): Q′ from Z_in′ = 1/(1/Z_bare + jωC_in(r, h)). Closed form: 1/Q′ = 1/Q_0 + (1 − s_feed)(1/Q_A − 1/Q_0). The s_feed implied by item 436's ×1 measurement (Q′ 360.604, Q_0 618.339, Q_A 254.469; `…item436-W6-read.log:54`) is 0.5002. With Q_0(h) 606.697 and Q_A 252.301 it predicts Q′ ≈ 356.44, **+41.3 %**; the band is item 436's **[+20 %, +70 %]** (s_feed 0.35 / 0.65 → +25.7 % / +61.2 %), ≥ 4× the 5 % target. Ceiling Q′ ≤ Q_0(h), i.e. ≤ +140.5 %.
- **(N-map) the 127.74 MHz map reused at 127.64 / 127.84 MHz,** pure arithmetic: 3.8440e-6 / 3.8410e-6 (drafter's arithmetic). It fails `STRIP_BAND` by ≈ 3 800×, against (Z′) ≈ 1e-16. Not run as Palace windows.
- **(N-nofeed) the feed branch dropped** (Z_in″ = Z_bare): predicted **no** − → + crossing on the triple, "f₀ undefined" (met at ×1 / 64 MHz, item 421; ×1 / 127.74 MHz, item 436 `:55`; rung / 64 MHz, item 444).

**Predictions (written first, printed).**
- W1 < 5 s; (E1′) 140 270; (Z′) ≤ 2e-16.
- W2–W5 each Palace Total ≈ 100 s (90–150), GMRES 45–60 iterations (`MaxIts` 100), peak 7.5–8.5 G.
- (T6′) / (T6-0‴) ≤ 1e-6 (central ≈ 2.5e-10); f₀ / Q / R_in gaps ≤ 1e-6; f₀,Pal/127.74 MHz − 1 ≈ −2.307e-06.
- (C-h) equal to ours' moves to ≈ five significant digits.
- (C2-h) ≥ 800× on S11 and ≥ 8 000× on Q; (K-esr‴) ratio 1 ± 1e-6; (N-feed) in [+20 %, +70 %]; (N-nofeed) no crossing.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call):
  - `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up.
  - `ls logs/ops61/ logs/ans7_exact/t127.74/`. Stop and report if any `t44_*.json` or `postpro_t44_*` exists.
  - The rung mesh, both templates and both `t127.74/tune_p2_h0.01*.json` must be present. Never re-create one; **never re-export the mesh**.
- **W1** (< 10 s), configs, (E1′), (Z′), (Z0):

  `scripts/testing/run_and_log.sh OPS-61-item448-W1-ans7-h001-128-configs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py ans7-h001-configs --f128'"`
- **W2** (≈ 100 s), the tuned point first, the price check:

  `scripts/testing/run_and_log.sh OPS-61-item448-W2-palace-t44-f127.74 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 t44_esr_f127.74.json'"`
- **Price gate (binding, after W2).** Run W3–W5 iff W2's Palace Total (Avg) ≤ 300 s (≤ 4 × 300 s fits the slot with the read). Otherwise the item lands as a partial record, and the review re-prices.
- **W3 / W4 / W5:** W2's command with `t44_esr_f127.64.json` / `t44_esr_f127.84.json` / `t44_lossless_f127.74.json`, labelled `OPS-61-item448-W3-palace-t44-f127.64` / `OPS-61-item448-W4-palace-t44-f127.84` / `OPS-61-item448-W5-palace-t44-lossless`. After any kill, run `docker compose -f docker/docker-compose.yml ps` before the next window.
- **W6** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item448-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --ans7-h001 --f128'"`
- **W7** (≈ 3 s), dispatch regression 1: item 444's W7 command verbatim (`20261014T094041Z_OPS-61-item444-W7-read.log:12`), labelled `OPS-61-item448-W7-regress-read-h001`. Its output section must `diff` empty against that log's.
- **W8** (≈ 3 s), dispatch regression 2: item 436's W6 command verbatim (`20261013T215115Z_OPS-61-item436-W6-read.log`, its `- Command:` line), labelled `OPS-61-item448-W8-regress-read-128`. Its output section must `diff` empty against that log's.

**Traps already paid for.** All of items 420 / 421 / 436 / 444's apply, plus:
- **Never reuse a ×1 or a 64 MHz value on the 127.74 MHz rung.** C_in(r, h) is 1.5288749690061608e-11 F. It is not `T35_C_IN_R_X1` 1.5275048999871143e-11 F (×1, 127.74 MHz), not `T43_C_IN_R_H_LITERAL` 6.364610922203437e-11 F (rung, 64 MHz) and not `T26_C_IN_R_X1` (item 419's "silent misread": a wrong literal under a right print fails nothing). Grep the W1 log for `1.52750` / `6.3646` / `6.3593` / `6.3590`: there must be none.
- **Two `tune_p2_h0.01*.json` pairs exist.** The top-level pair in `logs/ans7_exact/` is T37's (64 MHz); this item reads `logs/ans7_exact/t127.74/`. The `f_hz` == 127740000.0 assert catches the wrong directory.
- **The mesh name is the only mesh link.** A config that keeps `Model.Mesh == ops61_ans6_x1_hole_strip.msh` solves the ×1 problem at the rung's C and reads (T6′) ≈ 8e-4. The template diff asserts `Model.Mesh` changed, and the Palace log check's 140 270 catches it first.
- **The lossless template's `Freq` is 64 MHz** (`t39_lossless.json` `[0.064]`). Item 444 kept it; this item must change it to [0.12774], and (Z0) asserts so. A lossless run at 64 MHz reads (T6-0‴) ≈ O(1) against a 127.74 MHz comparand.
- **Every frequency gets its own exact map** (ruling (iii), R_p ∝ 1/f² off the conversion frequency); never one Palace sweep with fixed R_p / C_p. Never enter 0.2 Ω as a parallel R (item 420 read the convention: parallel).
- **The feed's (r + C) branch is in the read, not in Palace** (ruling 311-Q5). Dropping its r is (N-feed), the control, never the reading.
- **(C-h) / (C2-h) at 127.74 MHz use the 127.74 MHz ×1 rows.** Palace's are `…item436-W6-read.log:48`, ours' are `T42_*`, and the S11 numerator is `T42_S11_A[1]`. Never use `T43_PAL_X1_ROWS`, `T40_*` or `T39_S11_A_ESR`: those are 64 MHz values, and copying `read_ans7_h001` unedited would carry them in.
- **Comparand digits:** T38's `:2211–2212` moves are printed against item 387's `T35_Q_R_X1` / `T35_R_IN_R_X1`. The read recomputes ours' move against `T42_*` (item 419 W0); the two agree to the seven printed digits.
- Never run `export --hole-h001` or any `export*` or older `*-configs` mode. Never overwrite `t39_*` / `t40_*` / `t42_*` / `t43_*` or their outputs.
- **Compare S[1][1] only**; never read `z_matrix_ohm` or the `aed` key.
- **Docker socket and orphans:** as item 420. After any killed window, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0, and `docker compose -f docker/docker-compose.yml ps` must show no lingering `fem-em-palace`. Never lengthen a timeout or raise `-np`. A `GMRES … did not converge` line is a finding: that point is "not measured", and with a tune point missing (C1) cannot be read.
- **Record row and journal:** Palace's numbers in full (v0.18.1, the four configs, the mesh) beside ours; **no Palace figure set against an AED one in a tracked file**. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace's bare S11, f₀, R_in and Q at the 127.74 MHz ESR-tuned F-small node on the h = 0.01 exact-rectangle hole mesh, degree 2, Cu, read as T38 reads them; the lossless C-path tie on the rung at 127.74 MHz; Palace's own ×1 → 0.01 move at 127.74 MHz beside ours'. Not covered:
- an h-statement for Q (T41's, from F31 + T30);
- the two-frequency ratio reading (Q(r)-move / Q₀-move 0.4525 against 0.2915), which is the 10-17 weekly's;
- the h = 0.005 rung, whose `C_F` Route B twin is the weekly's to tier (item 444 (Q1));
- Route B rows (T38's `:2209`);
- order 1; AED; a measured ESR;
- an absolute tuned C, f₀, R_in or Q claim beyond "Palace = ours at matched order on one refined rung at 127.74 MHz";
- `src/` and any band.

**Status it can move.**
- The `ANS-7` row's Q clause gains "T44 (item 448, record, provisional): Palace at matched order on the h = 0.01 rung reads f₀ within `T21P_ZERO_RESOLUTION` and Q within T40's 5 % of T38's rows at 127.74 MHz, and its own ×1 → 0.01 Q move equals ours (−0.852 %)". With T43 that gives both Larmor frequencies. A record for T41. **No glyph:** the Done-when is Q's h-statement (T41). Rule (j): read; no mismatch.
- The `OPS-61` row (✅) gains a "step 3e record (T44)" clause.
- §10 T44 recorded DONE (provisional; the weekly ratifies or reverts).

**Negative result** (no `PREDICTED RED:` line; nothing here is expected red; any red stops the slot).
- **(E1′) red:** the file on disk is not item 444's rung mesh. Stop before any Palace window, revert the probe, park on `attempt/*`, BLOCKED (rule (d)). Never re-export in-slot; the review decides.
- **(Z′) / (Z0) red, or a W1 bit-equality / `f_hz` assert red:** the wrong C, map or directory. Stop before any Palace window, revert the probe, BLOCKED.
- **A W2–W5 log lacking 140 270 elements, σ = 5.800e+07 on 401, PEC on 1 only, or the right ω/2π:** park, BLOCKED.
- **(C1) red:** report Palace's Im Y_in triple; park, BLOCKED; the review reads it against (T6′).
- **(f₀) or (Q) red with (C1) green:** a disagreement on the external check on the refined rung at 127.74 MHz. Report the gaps, (T6′) and (C-h), open a known-issues entry, park, BLOCKED. Never widen, re-map or change r.
- **(f₀) green against `T22P_ZERO_RESOLUTION` but red against `T21P_ZERO_RESOLUTION`:** not predicted (gaps ≈ 1e-10). If it happens, BLOCKED as above, and the band question goes to the weekly with the number (ruling 2026-10-13 10:30 (3)).
- **(T6′) ≫ 1e-6, or (C-h) ≠ ours' move beyond ≈ 1e-8 on Q, with (f₀) / (Q) green:** a reading, not a failure. DONE as a record with a known-issues entry: the codes agree on T40's targets but not at matched order on the rung, so T41 cannot cite this item as "the discretisation's" at 127.74 MHz.
- **(T6-0‴) ≫ 1e-6:** a convention-free disagreement on the rung. DONE as a record with a known-issues entry (item 420's disposition); (f₀) / (Q) are still read.
- **W7 or W8 not diff-empty:** revert the probe, BLOCKED.
- **Price gate trips:** a partial record, not BLOCKED.

**Questions for the review / weekly.**
- (Q1) T44 is chain T's third rule (6) step since the 10-12 weekly, the last one allowed; no further provisional chain-T step can be written before the weekly.
- (Q2) If T44 lands green, both Larmor frequencies carry a matched-order Palace twin on the h = 0.01 rung. Whether T41 then needs Palace at h = 0.005 (item 444 (Q1)) or can write "unconverged, the discretisation's" from these two records is the weekly's ruling.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 448 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; configs and the read at `-n 1`, Palace at `-np 8`. Reuse item 444's `ops61_ans6_h0010_hole_strip.msh` (never re-export). One frequency per Palace window, each with its own map; C_in(r, h) = 1.5288749690061608e-11 F from `logs/ans7_exact/t127.74/` (never a ×1 or 64 MHz literal); the lossless config's `Freq` set to 127.74 MHz. The price gate after W2 binds. Probe only — no `src/`, no test module; never run an `export` or older `*-configs` mode; never overwrite `t39_*` / `t40_*` / `t42_*` / `t43_*`; never open `docs/private/` or `aed_results/`."

### Review ruling 2026-10-14 10:30 (item 448)

Accepted as written (T44, chain T's third and last provisional step before the 10-14 weekly).
- **Mesh guard.** (E1′) binds: re-read the exported h = 0.01 mesh's tet count and require 140 270 before writing any config. Do not re-export.
- **Copy traps.** The 127.74 MHz configs and read are additive copies (`--f128`). Every 64 MHz comparand item 444 hard-coded must be replaced by its 127.74 MHz counterpart (item 436's Palace ×1 rows, T38's rows): the lossless template's frequency, `T43_PAL_X1_ROWS`, `T40_*` and `T39_S11_A_ESR`. A comparand left at 64 MHz is a red to stop on, not to fix after a window.
- **Controls.** Every control is *predicted* (nothing has been measured at 127.74 MHz on this rung) and printed beside its measured value. (C-h) is printed, never gated.
- **Price gate.** If W2's first Palace point exceeds 300 s, park after it and return the price to the weekly.

### §9 item 454 — 2026-10-14 (§10 chain F37, daily — provisional; chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-14 18:00 daily review; the review reviews it.* Licence: §9 rule (6) (PROJECT_PLAN.md, grep "**Rule (6), the chain-extension licence**"): the next numbered step of live §10 chain F toward its remaining subgoal (the copper-loss rung and the public tables, PROJECT_PLAN.md grep "Chain F — closed on the h-statement"), tagged *(daily, 2026-10-14 — provisional)*; the weekly ratifies or reverts. It is chain F's **second** rule (6) step since the 10-12 weekly (F36 = item 437 the first). **F37 (provisional):** `OPS-61` step 2a / 2b's ×1-hole Palace twin (items 234 / 241) moved onto F26's copper-loss h-rung — `ANS-6` step 5b's h = 0.0075 exact hole mesh (`xl` entry 31, `docs/testing/logs/20261012T070009Z_ANS-6-step5b.log`) — exactly as F34 (items 418 / 434 / 435) moved step 1's strip twin onto the 0.0075 strip rung. Target: **F34's, imported, not new** — our identities (P) asserted on Palace's 4×4 at the imported bands; the per-class matched-order gap to entry 31's degree-2 class means printed beside the ≤ 1e-6 precedent (no cross-solver band exists: 2026-10-13 rulings (1) / (4)); Palace's own ×1 → 0.0075 move printed beside ours (C-h); the rung's Cu − PEC move (step 2b's (T3m) analogue) printed. The Palace directive (§5.4 (i)–(v)), the sheets-and-κ ruling (exact rectangles, κ off) and the 2026-10-13 rulings bind. **No band, tolerance, Done-when or validation target is set or moved.**

**Split (the price forces it; derivation under "Price").** Both columns need 4 excitations × 2 conductor models = eight single-excitation Palace windows plus a fresh export. That is ≈ 57 predicted min, past the minute-50 gate. **This item is the copper half**, because the copper-loss rung is the chain's subject. It writes the shared code for both columns (as item 434 wrote item 435's) and runs Cu only. **The PEC half (four windows, no probe edit) is the 2026-10-15 03:00 review's to queue**; its shape is in (Q1) below. No second item number is drawn here.

**`OPS-61` step 2c (§10 F37, Cu half): export the h = 0.0075 exact-rectangle hole mesh, then Palace order 2 with `Boundaries.Conductivity` (σ = 5.8e7 S/m) on wall 401 at 64 MHz, one excitation per window, against entry 31's degree-2 copper 4×4.** Ours is **not re-solved**: it is read from entry 31's tracked log.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`); no test module, no `src/`.
- **Tier:** heavy. Export, configs and read run at `-n 1` in the complex container; Palace runs at `-np 8` (the only width Palace has been measured at here). **Class:** chain step (solve), ≈ 8 measured min (the 10:30 restock's class median).
- **Time:** ≈ 48 predicted min. That is 15 fixed + ≈ 15 code + W1 ≈ 1.5 (60–90 s) + W2 < 1 + W3–W6 4 × ≈ 2.4 (140 s central, 115–180 s each) + W7–W9 < 1 + ≈ 5 spare.
- **Depends on:** nothing. Inputs on disk (checked read-only 2026-10-14) or tracked:
  - `logs/ops61/h05x_cu_o2.json` and `h05x_pec_o2.json` (the ×1 hole order-2 configs, items 241 / 234; templates, never overwritten);
  - `postpro_h05x_cu_o2/port-S.csv` and `postpro_h05x_pec_o2/port-S.csv` (Palace ×1 order 2, four excitations each).
  - Absent, as required: `ops61_ans6_h0075_hole_strip.msh`, any `h05x_h0075_*.json`, any `postpro_h05x_h0075_*`.

**What exists on our side (the comparand; not re-solved).** Entry 31, `20261012T070009Z_ANS-6-step5b.log`: `-n 16` `xl`, degree 2, κ off, 64 MHz only, Status 0. Ledger row "ANS-6-step5b" in `docs/testing/xl-ledger.md`. Line numbers below are into this log.
- **Mesh.** `:2093` reads `hole flag on: cells=226904 build 51.53 s`. `:3101` reads `hole mesh 226904 cells, mesh 52.0 s, -n 16; census {'exterior': 40212, 'cavity': 27352, 'cavity_exterior': 27352, 'outer': 12860, 'tag_collision': 0}`. `:3102` and `:4587` give **1 538 030** degree-2 unknowns. The mesh is `th14._build_ladder`'s call `_sheets_build(True, as_hole=True, resolution=0.0075, leg_sheet_strip_fraction=float(GATED_WIDTH_FRACTION))` (`tests/validation/test_th14_birdcage_copper.py:310–314`; knob `FEM_EM_EXACT_ANS6_RESOLUTION`, `tests/validation/test_ans6_exact_strip.py:189`, forwarded at `:491–500`).
- **PEC 4×4, 10 digits:** `:4556–4559` (rows `S_PEC_1k =` … `S_PEC_4k =`), under the header `:4555` `[TH-14 step2] f = 6.400000e+07 Hz: PEC sweep 250.2 s, …`.
- **Cu 4×4, 10 digits:** `:4561–4564` (rows `S_1k =` … `S_4k =`), under the header `:4560` `  sigma = 5.8e+07 S/m: Z_s = 2.087161e-03(1+j) ohm, sweep 260.4 s`. That is the same Z_s as item 241's (Z) tie, which read 0.0. σ_max **0.999904939998** at `:4565`.
- **Class means, 10 digits:** Cu at `:4589`, PEC at `:4590`.
- **Cu − PEC relative gap per class:** 0.0615 / 0.0279 / 0.0479 % (`:4592`).
- **Ours' ×1 → 0.0075 move:** Cu 0.1206 / 0.0506 / 0.1091 % (`:4593`) and PEC 0.1193 / 0.0501 / 0.1081 % (`:4594`). Both are against `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]` (`tests/validation/test_ans6_exact_strip.py:163–169`, `-n 8`).
- **σ ladder on the rung (rising, the (K3′) backing):** 0.999904939998 (Cu, `:4565`) < 0.999905161481 (5.8e9, `:4572`) < 0.999908196413 (5.8e11, `:4579`). The PEC-limit max\|dS\| is 4.139e-06 (`:4580`).
- **No tracked constant holds the rung's class means or its cell count.** `grep -rn 'h0075\|0p0075\|STEP5B\|226904' tests/` finds none for the hole route. `test_ans6_exact_strip.py` asserts only cells > `KNOB_CELL_FACTOR` × 80 798 (`:448`, `:769`). So the comparand is the tracked log, as item 418 read entry 27's.

**Price (from measurement; F34's rule).** The non-Setup Palace cost per (ND unknown × GMRES iteration) is measured on three fixtures:
- ×1 hole Cu o2: (175.950 − 30.604) / 4 = 36.3 s per excitation at 30–31 iterations on 570 332 ND, i.e. 2.10e-6 s (`20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log:318, :327`, iterations `:161–296`).
- 0.0075 strip, 64 MHz, e1: 185.767 − 32.911 = 152.9 s at 43 iterations on 1 822 094 ND, i.e. 1.95e-6 s (`20261013T170207Z_OPS-61-item434-W2-palace-h0075-f64-e1.log:101, :159, :181, :190`).
- h = 0.01 hole, single excitation: 94.502 − 16.631 = 77.9 s at 40 iterations on 968 016 ND, i.e. 2.01e-6 s (item 444 W3).

Scaled to the 0.0075 hole (1 538 030 ND, predicted equal to ours, as on every prior mesh):
- **Iterations:** 30–31 at ×1 hole. On the strip, ×1 → 0.0075 raised 37 → 43 (×1.16), so ≈ 35 here (bracket 30–45).
- **Per excitation, non-Setup:** 1 538 030 × 35 × 2.05e-6 ≈ 110 s (90–145).
- **Setup:** 16.6 s at 968 k ND and 32.9 s at 1.82 M ND, so ≈ 25–35 s.
- **One excitation per command: ≈ 140 s (115–180).** Peak ≈ 12.5 G (11–16 G; 14.8 G × 0.844 from item 434), against the service's 128 G.
- **All four excitations in one command: ≈ 470 s (390–615).** That is above F34's 450 s rule and §5.1's ≈ 480 s safe foreground ceiling, so this item runs one excitation per window. The cost is ≈ 3 × Setup (≈ 90 s), disclosed.
- **Both columns:** the PEC column is ≈ 0.91× Cu (×1 hole: 160.575 vs 175.950 s, `20261005T021936Z_OPS-61-item234-W4-palace-hole-pec-o2.log:324`), so ≈ 9 more minutes. With the export and the code, that is ≈ 57 min, which fails the minute-50 gate. Hence the split above.
- **Export:** the rung mesh built in 51.53 s at `-n 16` (gmsh is serial on rank 0; `:2093`). Item 418's 0.0075 strip export took 56.77 s at `-n 1` (window 61 s, `20261013T123326Z_OPS-61-item418-W1-export-h0075.log:1862, :1867`), and item 444's h = 0.01 hole export took 32.70 s (window 37 s). So W1 ≈ 55–75 s, window 60–90 s. Every window is ≪ the 590 s container timeout and the 1 200 s heavy ceiling.

**Code (probe only; additive; disclosed; written once, covering both columns).**
- **Pre-check for the code:** `grep -n "hole-h0075" scripts/probes/ops61_step1a_palace_ans4.py`. If it is dispatched, make no probe edit.
- **Constants** beside `T44_*`, each with its log-line provenance comment:
  - `H0075_HOLE_MESH_NAME = "ops61_ans6_h0075_hole_strip.msh"`;
  - `H0075_HOLE_CELL_RECORD = 226904`, a literal (`:2093`, `:3101`, `:4587`; no tracked constant holds it);
  - `H0075_HOLE_ND_P2 = 1538030` (`:3102`; a predicted Palace reading, not an anchor);
  - `H0075_HOLE_SETS = {col: {f"h05x_h0075_{col}_o2_e{k}": k for k in (1, 2, 3, 4)} for col in ("cu", "pec")}`;
  - `H0075_HOLE_TEMPLATES = {"cu": "h05x_cu_o2.json", "pec": "h05x_pec_o2.json"}`;
  - `H0075_HOLE_OURS_LOG` (the entry 31 log path);
  - `H0075_HOLE_OURS`, per column:
    - the header prefix (`"sigma = 5.8e+07 S/m: Z_s ="` for Cu, `"[TH-14 step2] f = 6.400000e+07 Hz: PEC sweep"` for PEC);
    - the row label (`"S_"` / `"S_PEC_"`);
    - the 10-digit class means (`:4589` / `:4590`);
    - the printed σ_max (Cu 0.999904939998; PEC `None`, since it is not printed at 10 digits);
    - the h-move (`:4593` / `:4594`);
    - the Palace ×1 output (`postpro_h05x_cu_o2` / `postpro_h05x_pec_o2`);
  - `H0075_HOLE_CU_PEC_GAP_OURS = {"self": 6.15e-4, "adjacent": 2.79e-4, "opposite": 4.79e-4}` (`:4592`).

  `H0075_RES`, `H0075_TIE_BAND`, `HOLE_MESH_NAME`, `T43_GLOBAL_TAGS` and `_hole_wall_attr()` are **reused, not copied**.
- **New `export_hole_h0075()`, dispatched by `export --hole-h0075`.** A **copy** of `export_hole_h001` (`:3595–3792`; provenance comment with commit and span) that writes **the mesh only**. `export_hole_h001` is neither edited nor called. Changes from that function:
  - `resolution=float(H0075_RES)`, `H0075_HOLE_MESH_NAME`, `H0075_HOLE_CELL_RECORD`, and the tag `[OPS-61 2c W1]`;
  - drop `assert H001_RES == H_RES`;
  - (E6) reads the cell-tag set from `h05x_cu_o2.json`'s `Domains.Materials`, the template actually used, instead of `t39_esr.json`, and also asserts it equals `T43_GLOBAL_TAGS`;
  - print the tet and element counts, never the file's byte size (the private-leak hook trap, items 434 / 435).

  Like its source, it asserts `comm.size == 1`, asserts the mesh is absent, and **asserts ncells == 226 904 before anything is written**.
- **New `hole_h0075_configs(col)`, dispatched by `hole-h0075-configs --cu` / `--pec`.**
  1. Assert `-n 1` (or no MPI). Assert the four `H0075_HOLE_SETS[col]` JSONs do not exist.
  2. **(E1′)** Parse `H0075_HOLE_MESH_NAME`'s `$Elements` type-4 count and assert it == 226 904 before any write (item 434's parse; print the count only).
  3. Load `H0075_HOLE_TEMPLATES[col]` and check it:
     - `Model.Mesh == HOLE_MESH_NAME`;
     - four ports, Attributes `[21]…[24]`, `+Z`, `R` 50, `Excitation` 1…4;
     - `Order == 2`, `Freq == [0.064]`;
     - **Cu:** `PEC.Attributes == [1]` and `Conductivity == [{"Attributes": [401], "Conductivity": th14.COPPER_SIGMA, "Permeability": 1.0}]`, with no `Thickness` / `External` key;
     - **PEC:** `PEC.Attributes == [1, 401]` and no `Conductivity` key.
  4. Write four configs, changing **only** `Model.Mesh` → `H0075_HOLE_MESH_NAME`, each port's `Excitation` → `true` on k and `false` on the other three (item 418's single-excitation form), and `Problem.Output` → `postpro_<name>`.
  5. Read each config back and assert the template diff is exactly those keys. This carries item 241's (Z) tie over unchanged (σ, μ_r and Freq are untouched; (Z) read 0.0, `20261005T093201Z_OPS-61-item241-W1-cu-configs.log`).
- **New `_ours_hole_rung_s(log, header_prefix, label)`.** A copy of `_ours_rung_s` (`:2268–2288`) with two changes: the header is matched by `startswith` on the stripped line (exactly one match asserted), and the row prefix is `f"{label}{k}k ="`. `_ours_rung_s` asserts `S_{k}k =` and **cannot parse the `S_PEC_` block**; it is not edited. Leading `PASSED` runs are stripped as before, and a `.log.gz` is opened if rotation compressed the log.
- **New `read_hole_h0075(col)`, dispatched by `read --hole-h0075 --cu` / `--pec`.** A copy of `read_exact_strip_h0075_freq` (`:2291–2402`, provenance comment) with the hole set, ours' block, the ×1 Palace output and `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"][col]` substituted. Its additions:
  - (Tie) on the class means (and on σ_max for Cu);
  - ours' and Palace's `_class_spread` printed side by side;
  - (K1′) / (K3′) in `--cu` mode, against ours' parsed PEC 4×4;
  - in `--pec` mode, if all four Cu columns are present, Palace-internal (K1) / (K3) and (T3m-h) (below);
  - the price section reading `*_OPS-61-item*-W*-palace-hole-h0075-{col}-e*.log` and printing measured Total / Setup / iterations / peak beside this brief's 115–180 s.

  `read_exact_strip_h0075_freq`, `read_hole_cu` and `read_hole_pec` are not edited. It never reads `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** five new exact-list `sys.argv` branches, first in `__main__` as items 436 / 444 / 448 placed theirs. Every existing mode keeps its behaviour, which W8 / W9 prove.

**Anchors (asserted; bands imported, never restated).**
- **Export, W1:**
  - **(E1)** tets == 226 904, asserted before the write and again parsed from the written file.
  - **(E2)–(E5)** as `export_hole_pec`'s:
    - no `CONDUCTOR_CELL_TAG`, and the phantom tag is present;
    - the census holds `cavity_exterior == cavity > 0`, `outer == exterior − cavity` and `tag_collision == 0` (predicted 40 212 / 27 352 / 12 860, as `:3101`);
    - box area == `EXTERIOR_AREA_M2` within `EXTERIOR_AREA_BAND`;
    - per port, the N_i area / width / height == `SPEC_AREA_M2` / `SPEC_WIDTH_M` / `SPEC_HEIGHT_M` within `STRIP_BAND` (`tests/validation/test_ops61_exact_strip_ports.py:58`). These are the exact rectangles at the rung, which is the precondition of the sheets-and-κ ruling.
  - **(E6)** the cell-tag set equals `h05x_cu_o2.json`'s Materials set, which equals `T43_GLOBAL_TAGS`.
- **(E1′), W2:** the type-4 count parsed from the file == 226 904.
- **(Tie), W7: a parse check, not a physics band.**
  - The parsed Cu 4×4's σ_max equals 0.999904939998 (`:4565`), and its class means equal `:4589`, each within `H0075_TIE_BAND` 1e-9 absolute.
  - The parsed PEC 4×4's class means equal `:4590` within the same band.
  - A wrong block misses by much more. Cu against the 5.8e9 block differs by ≈ 3.7e-4 on self (`:4566` / `:4573`), and Cu against PEC by ≥ 1.2e-4 (adjacent).
- **(P), W7, on Palace's assembled Cu 4×4** (`_p_check_s`, item 418's):
  - ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py:86`);
  - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every column power ≤ 1 + 1e-9 (`tests/validation/test_port_birdcage_four_port.py:129`);
  - every S-class spread ≤ `ADJACENT_SPREAD_BAND` (`tests/validation/test_port_birdcage_lumped_column.py:121`);
  - pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst (`test_port_birdcage_four_port.py:173`).
  - Reciprocity is cross-window: S_ij and S_ji come from two Palace commands.
- **Palace log check (executor, journalled), W3–W6:**
  - 226 904 elements, and the bounding box ±0.12 / ±0.12 / ±0.10 m;
  - `Rs = 4.375e+01 Ω/sq` on 21–24;
  - `Configuring Dirichlet PEC BC at attributes:` listing **1 only**;
  - `Configuring Robin finite conductivity BC at attributes:` / `401: σ = 5.800e+07 S/m`, quoted verbatim (item 241 W3's form, `20261005T093223Z_OPS-61-item241-W3-palace-hole-cu.log:78–79`);
  - **no** PMC-fallback warning;
  - ω/2π = 6.400e-02 GHz;
  - exactly one excited port, k;
  - one `GMRES solver converged` line;
  - `ND (p = 2): 1538030` predicted (= ours; a mismatch is a reading, not a red).

**Ours against Palace (printed, never gated — 2026-10-13 rulings (1) / (4)).**
- **(T5-Cu)** per class and per entry, \|Palace o2 Cu h − ours p2 Cu h\|/\|ours\|. **Predicted ≤ 1e-6.** The measured precedents are item 241 (T4) 1.6–2.0e-10 on the ×1 hole Cu and items 434 / 435 ≤ 1.3e-10 on the strip rung. The floor is the 10-digit print, ≈ 1e-10. Ours is `-n 16`; Palace is `-np 8`.
- **(C-h-Cu)** Palace o2 Cu at h = 0.0075 against Palace o2 Cu at ×1 (`postpro_h05x_cu_o2`), per class, printed beside ours' 0.1206 / 0.0506 / 0.1091 % (`:4593`). **Predicted equal to the printed digit** (items 434 / 435 / 437 measured exactly that on the strip rungs). Palace ×1 Cu o2 equals our ×1 record to ≤ 2.0e-10 (item 241 (T4)), so both codes start from the same point.
- Ours' and Palace's `_class_spread` per class, side by side. **Predicted equal to the printed digit** (the two 4×4s tie entrywise). The xl log's `:4565` spreads (0.0009 / 0.0013 / 0.0006 %) use th14's definition, not `_class_spread`; item 241's prediction miss was exactly this. Never compare across the two definitions.
- Palace Total, Setup, iterations and peak per window, beside 115–180 s.

**Negative controls (rule (e)).**
- **(K1′) The conductivity boundary reaches Palace's operator — floor *asserted*, size predicted.**
  - Asserted: per class, \|Palace Cu h − ours PEC h\| / \|ours PEC h\| > `KNOB_REACHES_THE_SOLVE` (1e-6, `tests/validation/test_port21_feed_sensitivity.py:156`, imported).
  - **Same-fixture prior:** ours' Cu − PEC gap on this mesh is 6.15e-4 / 2.79e-4 / 4.79e-4 (`:4592`), ≥ 279× above the floor. Palace's PEC equals ours to ≤ 2.2e-10 at ×1 (item 234 (T1)).
  - Printed: the separation \|Palace Cu − ours PEC\| / \|Palace Cu − ours Cu\|. **Predicted ≥ 279×** if (T5-Cu) ≤ 1e-6, central ≈ 3e6×.
  - A Palace that treated 401 as PEC reads ≈ 1e-10 here. A wall left out of every block reads ≥ 100 % (item 234's (K′)), and the PMC-fallback warning shows first.
- **(K3′) Loss direction — *asserted*:** σ_max(Palace Cu h) < σ_max(ours PEC h). The read computes the latter from `:4556–4559` and prints it.
  - **Same-fixture prior:** ours' σ ladder on the rung rises monotonically toward PEC (`:4565` < `:4572` < `:4579`).
  - Predicted margin ≥ 3.3e-6 (0.999908196 − 0.999904940), against a ≈ 1e-9 tie.
  - A wall with a sign error (active, not lossy) puts σ_max above PEC's.
- **(C2-h) Mesh-mismatch separation — *predicted*, printed** (no prior Palace measurement on this mesh).
  - Per class: \|ours ×1 p2 Cu (`ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["cu"]`) − Palace Cu h\| / \|ours h p2 Cu − Palace Cu h\|.
  - The numerator ≈ the h-move. The weakest class is adjacent, 5.06e-4 relative (≈ 2.3e-4 absolute).
  - **Predicted ≥ 506×** against item 229's ≥ 25× bar, central ≈ 5e6×. Items 434 / 435 measured 1.2e7–5.8e7× on the strip rung.
  - A config left on the ×1 hole mesh reads ≈ 1 here, and the 226 904 element check catches it first.
- (P)'s power to discriminate was shown by item 159's flip (200.0204 % spreads against the 0.5 % band); not re-run.

**Predictions (written first, printed).**
- W1: 55–75 s of mesh and export (window 60–90 s); (E1) 226 904; (E4) to ≈ 1e-16; (E5) ≤ 2e-15.
- W2: < 5 s.
- W3–W6, each: Palace Total ≈ 140 s (115–180), GMRES 30–45 iterations (`MaxIts` 100), peak 11–16 G.
- (P): reciprocity ≈ 1e-11; σ_max 0.999904940 ± 1e-9; spreads equal to ours' `_class_spread` to the printed digit.
- (T5-Cu) ≤ 1e-6, central ≈ 1e-10; (C-h-Cu) equal to 0.1206 / 0.0506 / 0.1091 % to the printed digit; (K1′) separation ≥ 279×; (K3′) margin ≥ 3.3e-6; (C2-h) ≥ 506×.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not a harness call):
  - `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up and no lingering `fem-em-palace`. If an `xl` window is still computing, start nothing: journal and stop.
  - `ls logs/ops61/`. Stop and report if `ops61_ans6_h0075_hole_strip.msh`, any `h05x_h0075_*.json` or any `postpro_h05x_h0075_*` exists. Both `h05x_*_o2.json` templates and both `postpro_h05x_*_o2/port-S.csv` must be present.
  - `docs/testing/logs/20261012T070009Z_ANS-6-step5b.log` (or `.log.gz`) must be present.
- **W1** (60–90 s), the export, (E1)–(E6):

  `scripts/testing/run_and_log.sh OPS-61-item454-W1-export-hole-h0075 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 400 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py export --hole-h0075'"`

  If the private-leak hook refuses the W1 log on a printed integer, quarantine it unedited under the gitignored `logs/ops61_item454_quarantine/` (items 434 / 435). The test-results row stays, and every W3–W6 Palace log corroborates (E1).
- **W2** (< 10 s), the Cu configs and (E1′):

  `scripts/testing/run_and_log.sh OPS-61-item454-W2-hole-h0075-configs-cu "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py hole-h0075-configs --cu'"`
- **W3** (≈ 140 s, ≤ 590), excitation 1, the price check:

  `scripts/testing/run_and_log.sh OPS-61-item454-W3-palace-hole-h0075-cu-e1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 h05x_h0075_cu_o2_e1.json'"`
- **Price gate (binding, after W3).** Run W4–W6 iff W3's Palace `Total` (Avg) T₁ ≤ 300 s. Otherwise columns 2–4 are "not measured — price": W7 still runs on column 1 ((T5-Cu) per entry; (P) "not measured"), and the item lands as a partial record, not BLOCKED; the review re-prices. A T₁ above 300 s would contradict this brief's 115–180 s by ≥ 1.67×: record it as a price finding.
- **W4 / W5 / W6:** W3's command with `…_cu_o2_e2.json` / `…_e3.json` / `…_e4.json`, labelled `OPS-61-item454-W4-palace-hole-h0075-cu-e2` / `-W5-…-cu-e3` / `-W6-…-cu-e4`.
  - **Minute gate:** if W6 cannot start by minute 50, stop launching. The `postpro_*` directories are durable, so the next slot runs the missing excitations with no edit and then W7 (rule (d), "partial; columns k… pending").
  - After any kill, run `docker compose -f docker/docker-compose.yml ps` before the next window.
- **W7** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item454-W7-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --hole-h0075 --cu'"`
- **W8** (≈ 3 s), dispatch regression 1: item 241's W6 command verbatim (`20261005T093623Z_OPS-61-item241-W6-read.log:12`), labelled `OPS-61-item454-W8-regress-read-hole-cu`. Its output section must `diff` empty against that log's.
- **W9** (≈ 3 s), dispatch regression 2: item 448's W6 command verbatim (`20261014T171034Z_OPS-61-item448-W6-read.log`, its `- Command:` line), labelled `OPS-61-item454-W9-regress-read-ans7-h001-128`. Its output section must `diff` empty against that log's.

**Traps already paid for.** All of items 234 / 241 / 418 / 434 / 444's apply, plus:
- **The harness call is the whole Bash command**: bare, repo-relative, nothing chained or piped (`bash_guard.py`; the §9 "Every item" block). An executor that loses the docker socket mid-item stops and reports.
- **Never run** `export`, `export --exact-strip`, `--hole-pec`, `--exact-strip-h0075`, `--exact-strip-h0050` or `--hole-h001`, or any older `*-configs` mode, and **never `hole-h0075-configs --pec`** (that is the PEC half's W1). Never overwrite `h05x_*`, `f05x_*`, `t39_*`–`t44_*` or their outputs.
- **The export runs at `-n 1`** (the function asserts it); a gmsh failure deadlocks wider runs. 226 904 is the count from `-n 16`. `-n 1` exports reproduced the wider-width counts twice: item 418 (282 198, entry 27) and item 444 (140 270, T37's `-n 12`). A miss here is (E1) red, never a re-record.
- **Parse trap.** The log holds **three** `S_1k =` blocks: Cu at `:4561`, 5.8e9 at `:4568` and 5.8e11 at `:4575`. It holds one `S_PEC_1k =` block, at `:4556`. Anchor Cu on the unique `sigma = 5.8e+07 S/m: Z_s =` header (`:4560`) and PEC on the `[TH-14 step2] f = 6.400000e+07 Hz: PEC sweep` prefix (`:4555`); never on the timings in either line. (Tie) catches a wrong block.
- **Attribute 401 goes in exactly one block.** For Cu that is `Conductivity` with PEC on 1 only. If 401 stays in `PEC`, the run reproduces the PEC column and (K1′) reads ≈ 1e-10. If it is dropped from both, the wall is PMC and the log shows the fallback warning. Never write `Thickness` / `External`, and never fall back to `Impedance`.
- **The mesh name is the only mesh link.** A config left on `ops61_ans6_x1_hole_strip.msh` solves ×1 and reads (T5-Cu) ≈ the h-move. The template diff asserts the change, and the log check's 226 904 catches it first.
- **One excitation per Palace command.** Never raise `-np` or lengthen a timeout in-slot. A Status 124 is a price finding: record it and run no further e-window. A `GMRES … did not converge` line is a finding, not a retry: that column is "not measured".
- **Orphans after any killed window:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0, and `docker compose -f docker/docker-compose.yml ps` must show no lingering `fem-em-palace`.
- **Spread definitions differ** (item 241's miss): compare `_class_spread` with `_class_spread` only.
- **Compare S only.** Never read Palace's `z_matrix_ohm` or the `aed` key. Never register a record from Palace or from this read. The record row and journal give Palace's numbers in full (v0.18.1, the four configs, the mesh) beside ours; **no Palace figure is set against an AED one in any tracked file**. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: the h = 0.0075 exact-rectangle hole mesh as a Gmsh 2.2 export; Palace order 2 with `Boundaries.Conductivity` (Cu) on 401 at 64 MHz on it; (P) on the assembled Cu 4×4; (T5-Cu) / (C-h-Cu) against entry 31's degree-2 rung; the (K1′) / (K3′) / (C2-h) controls; the measured price. Not covered:
- **the PEC column on the rung and the Palace-internal Cu − PEC move (T3m-h)**, which is the PEC half for the 03:00 review (Q1);
- 10 / 128 MHz on the rung;
- the h = 0.005 hole rung (F31, `xl` entry 34, no record yet);
- order 1; the σ-scaling (`cu11`) control; Palace's surface-loss postprocessing and the loss partition;
- any h-statement for Palace (two hole meshes are not a ladder);
- the tuned twin; AED; `src/` and any band.

**Status it can move.**
- The `OPS-61` row (✅, PROJECT_PLAN.md grep '^| `OPS-61`') gains a "step 2c record (F37 Cu, §9 item 454)" clause: Palace's conductivity boundary on a second, refined hole mesh, with (P), (T5-Cu) per class and (C-h-Cu). **No glyph change.** Rule (j): the row's Done-when is step 1c's and is met; a record sentence needs no amendment.
- If (T5-Cu) ≤ 1e-6 and (C-h-Cu) equals ours to the printed digit, the §2 / CLAUDE.md wording "Palace … matches ours to 0.0000 % at matched order" may extend to "the copper (Leontovich) hole route on two meshes (×1, h = 0.0075)". That extends the public external S-level check of the Leontovich term to F26's rung. The call is the weekly's.
- §10 F37 is recorded "Cu half DONE (provisional)"; the PEC half is the 03:00 review's to queue.

**Negative result** (no `PREDICTED RED:` line; nothing here is expected red; any red stops the slot).
- **(E1) red:** the export is not entry 31's mesh. Stop before any Palace window, park on `attempt/*`, BLOCKED (rule (d)). Open a known-issues entry only if the count moved on an unchanged image; a re-record is `record-reconciler`'s, under (1*).
- **(E2)–(E6) red:** the strip or the region tags did not reach the rung. Report the spans, facet counts and tag sets; park, BLOCKED.
- **W1 Status 124:** a price finding. Park and return the export's tier to the review; never lengthen the timeout.
- **(E1′) or a W2 template-diff / conductor-block assert red:** the wrong mesh or template. Stop before any Palace window, revert the probe, BLOCKED.
- **Palace aborts, or a W3–W6 log lacks 226 904 elements / `401: σ = 5.800e+07` / PEC on 1 only / 6.400e-02 GHz, or shows the PMC-fallback warning:** record it verbatim, open a known-issues entry, park, BLOCKED.
- **(Tie) red:** a parse bug. Fix it once in-slot before W7 is re-run (the read's own code, not the windows); still red, BLOCKED.
- **(P) red:** report every identity's number, park, BLOCKED. Never widen. Ours' (G) is green on this mesh (`:4565`, the 2026-10-13 closure audit), so this would be the first identity-level disagreement on a refined hole mesh.
- **(K1′) red** (Palace Cu within 1e-6 of ours' PEC): the conductivity block is not reaching the operator. Park, BLOCKED.
- **(K3′) red:** the wall adds energy, which is a sign-convention finding. Report it, park, BLOCKED; the review rules on it.
- **(T5-Cu) ≫ 1e-6, or (C-h-Cu) ≠ ours beyond the print, with (P), (Tie), (K1′) and (K3′) green:** a reading, not a failure. DONE as a record with a known-issues entry: the first matched-order disagreement on the copper route on a refined mesh. PEC agreed to ≤ 2.2e-10 at ×1, so it would point at the Robin term's discretisation on the finer wall; the PEC half then separates wall model from mesh.
- **W8 or W9 not diff-empty:** revert the probe, BLOCKED.
- **Price gate or minute gate trips:** a partial record, not BLOCKED.

**Questions for the review / weekly.**
- **(Q1) The PEC half, for the 2026-10-15 03:00 review.** Its shape is this brief with `--cu` → `--pec`, **no probe edit** (the dispatch exists), and no export (the mesh is reused; (E1′) re-reads 226 904).
  - Windows: W1 `hole-h0075-configs --pec`; four windows `h05x_h0075_pec_o2_e{1..4}.json` at ≈ 105–165 s each (0.91× Cu); W6 `read --hole-h0075 --pec`; a regression re-run of this item's W7, diff-empty.
  - Its anchors: (P) and (Tie) on PEC. Palace-internal (K1) (Cu vs PEC > `KNOB_REACHES_THE_SOLVE`) and (K3) (σ_max cu < pec) are *asserted*, with this item's columns as the same-fixture prior.
  - Printed: (T5-PEC); (C-h-PEC) beside ours' 0.1193 / 0.0501 / 0.1081 % (`:4594`); and **(T3m-h)**, Δ_Pal = Palace Cu − Palace PEC against Δ_ours = ours Cu − ours PEC on the rung, \|Δ_Pal − Δ_ours\|/\|Δ_ours\| and arg, predicted ≤ 1e-3 and \|arg\| < 0.1° (×1: item 241 measured ≤ 4.4e-7, 0.0000°).
  - ≈ 30 predicted min.
  - The 10-15 02:00 `xl` window (entry 34, F31) may still be running at 04:30, and the pre-check's stop rule then binds.
- (Q2) F37 is chain F's second rule (6) step since the 10-12 weekly. The 10-14 weekly (21:00 tonight) ratifies or reverts it. If it reverts, the exported mesh and the Cu configs stay as gitignored artefacts, and nothing tracked depends on them beyond this item's record.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 454 per `docs/planning/chunks/OPS-61.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; the export and the configs at `-n 1`, one Palace excitation per window at `-np 8`, never two heavy windows at once; the price gate after W3 and the minute-50 gate bind. Copper only: never run `hole-h0075-configs --pec`. Probe only — no `src/`, no test module; never run an older `export` or `*-configs` mode; never overwrite `h05x_*`, `f05x_*` or `t39_*`–`t44_*`; never open `docs/private/` or `aed_results/`."

### Review ruling 2026-10-14 18:00 (item 454)

Accepted as drafted, as the **copper half** of §10 F37 (daily 2026-10-14 — provisional, §9 rule (6)). F37 is chain F's second provisional step since the 10-12 weekly (F36, F37).
- **Citations checked by this review.** The comparand rows exist as cited: the Cu / PEC class means at `20261012T070009Z_ANS-6-step5b.log:4589–4590`, the Cu − PEC gap at `:4592`, the h-moves at `:4593–4594`, and the 226 904-cell build in 51.53 s at `:2093`.
- **The split binds.** The item runs the shared code, the export and the four Cu windows only. The PEC half is the same step's second half, not a new step. The 2026-10-15 03:00 review queues it if the weekly has not reverted F37. No PEC window runs in this slot, even with time left.
- **Records only.** The rung has no tracked record constant, so (Tie) reads the xl log's 10-digit block, as item 418 read entry 27's. (E1) red, meaning the `-n 1` export missing 226 904 cells, stops the item. It never re-records.
- **xl overlap.** If a slot reaches this item while an `xl` window is computing (02:00–≈ 04:00), it stops at the brief's pre-check. Two heavy windows never run at once.
- **Withdrawn unrun** if tonight's weekly reverts F37.

**Outcome 2026-10-15 00:46 UTC (00:00 slot): DONE (Cu half).** W1 `-n 1` export 226 904 tets, (E1)–(E6) PASS (51 s window); W2 (E1′) PASS, template diff = Mesh / Excitation / Output; W3–W6 Palace v0.18.1 `-np 8`, every log check as briefed (226 904 elements, `401: σ = 5.800e+07 S/m`, PEC on 1 only, no warning, ND 1 538 030 = ours), Total 142.9–144.9 s, Setup 26.6–27.3 s, 39–40 GMRES its, peak 12.3–12.4 G. W7: (Tie) PASS (5.59e-11); (P) PASS (reciprocity 1.546e-11, σ_max 0.999904940000, `_class_spread` = ours' to the printed digit); (T5-Cu) per class 1.274e-10 / 4.521e-11 / 6.621e-11; (C-h-Cu) 0.1206 / 0.0506 / 0.1091 % = ours'; (C2-h) ≥ 9.47e6×; (K1′) PASS, separation ≥ 4.83e6×; (K3′) PASS, margin 6.828e-6 (ours' parsed PEC σ_max 0.999911767900, larger than the brief's 5.8e11-ladder proxy, so the margin is ≈ 2× the predicted floor). W8 / W9 diff-empty. Logs `20261015T003312Z_…-W1-…` through `20261015T004452Z_…-W9-…`. The PEC half (Q1) needs no probe edit.

## §7 row moved verbatim 2026-10-14 (rotate_plan_archive.py chunks)

 Palace as a second independent reference | ✅ **2026-10-05 (weekly, reading (a); audit stamp owed to the 03:00 daily on item 194's window `20261003T110728Z_OPS-61-item194-W5-read.log`).** Steps 0–1e and 2a–2b are records: Palace matches ours to 0.0000 % per class at matched order on the exact strip mesh at 10 / 64 / 128 MHz at both orders, and on `ANS-6`'s exact hole mesh at 64 MHz with the PEC and copper (`Boundaries.Conductivity`) columns to ≤ 2.6e-10. Step 1f record (F34, §9 item 418): Palace order 2 on `ANS-4`'s h = 0.0075 exact-strip rung (282 198 cells, 1 822 094 ND unknowns in both codes) at 10 MHz, assembled from four single-excitation runs, passes (P) (reciprocity 1.5e-11 across windows) and matches our degree-2 rung to ≤ 1.9e-10 per class (T5); price 146 s per excitation ⇒ 64 / 128 MHz run heavy. Step 1g record (F34 64 MHz, §9 item 434): the same rung at 64 MHz, four single-excitation runs (181–186 s Palace Total each, 43 GMRES iterations) — (P) green (reciprocity 2.0e-11 across windows, σ_max 0.999823496779, spreads ≤ 0.0172 %), (T5) per class ≤ 8.0e-11 against our degree-2 rung (printed, no cross-solver band), Palace's own ×1 → 0.0075 h-move equal to ours to the printed digit (`20261013T171447Z_OPS-61-item434-W6-read.log`). Step 1h record (F34 128 MHz, §9 item 435): the same rung at 128 MHz, four single-excitation runs (266–284 s Palace Total each, 72 GMRES iterations) — (P) green (reciprocity 3.3e-11 across windows, σ_max 0.999320852432, spreads ≤ 0.0220 %), (T5) per class ≤ 1.3e-10 against our degree-2 rung (printed, no cross-solver band), Palace's own ×1 → 0.0075 h-move 0.2204 / 0.1590 / 0.2863 % equal to ours to the printed digit (`20261013T173812Z_OPS-61-item435-W6-read.log`); F34 is complete at 10 / 64 / 128 MHz. Step 3a record (T39, §9 item 420): Palace's first lumped R / C run — the ESR-tuned F-small (r = 0.2 Ω, C_in(r)) on `ANS-6`'s ×1 hole mesh, order 2, Cu, 64.000 MHz, port 1 the bare 50 Ω feed and ports 2–4 unexcited `LumpedPort` R / C; convention reading **parallel** (Palace's lumped R, C form R ‖ C; the 0.2 Ω-as-R control reads 1.437e-1), bare S11 against our T26 Route A (T6) 2.569e-10, the convention-free C path (T6-0) 2.565e-10, 51.8 s Palace Total per frequency (printed, no cross-solver band — Q1 to the 10-14 weekly). Step 3b record (T40, §9 item 421): Palace's driven read of that configuration at T26's three tune points (63.9 / 64.0 / 64.1 MHz, one exact parallel map per frequency, ≈ 52 s each) through `parallel_rows`: one − → + crossing, f₀ +6.530e-11, Q −4.871e-9, R_in −7.456e-10 relative to T26's Route A rows, per-point S11 ties ≤ 2.580e-10. Step 3c record (T42, §9 item 436, provisional): the same at 127.74 MHz against T35's Route A rows (127.64 / 127.74 / 127.84 MHz, ≈ 56 s Palace Total each, 42 GMRES iterations): one − → + crossing, f₀ +6.073e-11, Q +5.444e-9, R_in +2.906e-9, per-point S11 ties ≤ 2.373e-10, lossless C path 2.308e-10. Reading (a): the two FEM codes share the port model and both sit within 1 % of AED First Order on the exact mesh, so the residual against AED is discretisation, not a port-model difference (§10 rulings 2026-10-05). Step 1i record (F36, §9 item 437, provisional): Palace order 2 on the h = 0.005 exact-strip rung (593 478 elements, ND 3 833 208 = ours) at 10 MHz, four single-excitation runs (304–334 s Palace Total each, 25 GMRES iterations, peak 32.8 G) — (P) green (reciprocity 1.654e-11, σ_max 0.999995535243, spreads 0.0280 / 0.0081 / 0.0174 %); (T5) per class against `ANS-4` step 6c's degree-2 rung 6.308e-11 / 5.760e-11 / 1.679e-11 (every entry ≤ 1.68e-10); (C-h) Palace's 0.0075 → 0.005 move 0.0144 / 0.0088 / 0.0161 %, equal to ours' to the printed digit; (C2-h) ≥ 1.53e6× (`20261013T221851Z_OPS-61-item437-W6-read.log`). Step 3d record (T43, item 444, provisional): the h = 0.01 exact-rectangle hole mesh exported (140 270 tets, (E1)–(E6) green), Palace v0.18.1 order 2 at 63.9 / 64.0 / 64.1 MHz (ESR maps at C_in(r, h)) and the 64.0 MHz lossless C path, 93.8–96.9 s Total each at `-np 8`, ND 968 016 = ours; (C1) / (f₀) / (Q) green against T37's rows, (C2-h) 3.198e6× on S11 and 1.091e6× on Q (`20261014T094041Z_OPS-61-item444-W7-read.log`). Step 3e record (T44, item 448, provisional): the same h = 0.01 hole mesh reused, not re-exported ((E1′) 140 270 re-read), Palace v0.18.1 order 2 at 127.64 / 127.74 / 127.84 MHz (ESR maps at C_in(r, h) 1.5288749690061608e-11 F) and the 127.74 MHz lossless C path, 106.6–112.3 s Total each at `-np 8`, 50–52 GMRES iterations, ND 968 016 = ours; (C1) / (f₀) / (Q) green against T38's rows, (C2-h) 3.538e6× on S11 and 8.366e5× on Q (`20261014T171034Z_OPS-61-item448-W6-read.log`). **Step 2c record (F37 Cu, §9 item 454, 2026-10-15):** Palace's conductivity boundary on a second, refined hole mesh — the h = 0.0075 exact-rectangle hole rung (`ANS-6` step 5b, `xl` entry 31; `-n 1` export 226 904 tets, (E1)–(E6) green), Palace v0.18.1 order 2, `Boundaries.Conductivity` σ = 5.8e7 S/m on 401, 64 MHz, configs `h05x_h0075_cu_o2_e{1..4}.json`, one excitation per window at `-np 8` (142.9–144.9 s Total, 39–40 GMRES its, ND 1 538 030 = ours): (P) green on the assembled Cu 4×4 (reciprocity 1.546e-11, σ_max 0.999904940000); (T5-Cu) per class 1.274e-10 / 4.521e-11 / 6.621e-11 against entry 31's degree-2 Cu block; (C-h-Cu) 0.1206 / 0.0506 / 0.1091 %, equal to ours' h-move to the printed digit; (K1′) / (K3′) asserted green (`20261015T004434Z_OPS-61-item454-W7-read.log`). No glyph change. **Done-when:** Palace `-np 8` solves `f05x` at 64 MHz: 117 206 elements, Rs 43.75 Ω/sq, port i = N_i, `+Z`, 50 Ω. Its 4×4 meets `RECIPROCITY_BAND`, σ_max and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every S-class spread ≤ `ADJACENT_SPREAD_BAND` with pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst — each asserted with the band imported. *History: `docs/planning/chunks/OPS-61.md`.* 

### §9 item 457 — 2026-10-15 (§10 chain F, F37 PEC; chain step, solve)

**Outcome 2026-10-15 09:52 UTC (09:30 slot): DONE (PEC half; F37 complete).** No probe edit. W1 (`20261015T094106Z_OPS-61-item457-W1-hole-h0075-configs-pec.log`, 2 s): (E1′) 226 904 tets PASS; template PEC [1, 401], Conductivity null, ports 21–24 `+Z` 50 Ω, Order 2, Freq [0.064] PASS; four configs, diff = Mesh / Excitation / Output. W2–W5 Palace v0.18.1 `-np 8` (`…T094114Z_…-W2-…-pec-e1.log`, `…T094341Z_…-W3-…-pec-e2.log`, `…T094604Z_…-W4-…-pec-e3.log`, `…T094828Z_…-W5-…-pec-e4.log`): every log 226 904 elements, box ±0.12 / ±0.12 / ±0.10 m, Rs 4.375e+01 Ω/sq on 21–24, `Configuring Dirichlet PEC BC at attributes:` `1, 401`, zero `Robin finite conductivity` / PMC lines, ω/2π 6.400e-02 GHz, one excited port k, `ND (p = 2): 1538030` (= Cu, = ours), GMRES converged 40 / 39 / 40 / 40 its, Total 138.1 / 138.1 / 139.7 / 135.5 s, Setup 25.1–26.8 s, peak 12.1–12.3 G (predicted 115–165 s: met; the price gate, T₁ 138.1 s ≤ 300 s, passed). W6 (`20261015T095049Z_OPS-61-item457-W6-read.log`): (Tie) PASS on both blocks (5.590e-11); **(P) PASS** — reciprocity 1.529e-11, σ_max 0.999911767864 (ours' computed 0.999911767900, Δ 3.6e-11), column powers 0.99929625–0.99929635, spreads 0.0009 / 0.0013 / 0.0002 %, pooled separation 32 046×; **(T5-PEC)** per class 1.378e-10 / 1.153e-10 / 1.484e-10, per entry ≤ 2.512e-10 (predicted ≤ 1e-6, central ≈ 1e-10: met); `_class_spread` Palace = ours to the printed digit; **(C-h-PEC)** 0.1193 / 0.0501 / 0.1081 % = ours' (predicted); (C2-h) 8.661e6 / 4.346e6 / 7.28e6× (predicted ≥ 501×); **(K1) PASS** 6.155e-4 / 2.793e-4 / 4.789e-4 (= ours' Cu − PEC gap to the printed digit); **(K3) PASS** 0.999904940000 < 0.999911767864 (margin 6.828e-6); **(T3m-h)** 3.581e-8 / 2.515e-7 / 1.790e-7, arg +0.0000° (predicted ≤ 1e-3, central ≈ 3e-7: met) — the Leontovich term's own contribution matches across codes on the refined mesh. W7 (`20261015T095055Z_OPS-61-item457-W7-regress-read-hole-cu.log`) output section diff-empty against item 454 W7's. Disclosed, not edited (ruling 2026-10-15 03:00): the shared read path prints Cu-worded labels — the (C2-h) "predicted >= 506x" (PEC's floor is ≥ 501×) and the price "brief's projection 115-180 s" (this brief's bracket 115–165 s). Configs `logs/ops61/h05x_h0075_pec_o2_e{1..4}.json` are gitignored, as item 454's.

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.* Licence: §9 rule (4), the enumerated-chain licence (PROJECT_PLAN.md:3805–3813). F37 PEC is a numbered step the 2026-10-14 weekly enumerated in §10 chain F (PROJECT_PLAN.md:4701–4705) and assigned to this review ("the 10-15 03:00 review queues it"). Its rulings (PROJECT_PLAN.md:4633–4634) read "four windows, no probe edit. F37 Cu ratified". **Not provisional.** Its serial prerequisite has landed: item 454 (F37 Cu, DONE 10-15 00:00 slot, PROJECT_PLAN.md:4032; outcome this file :2076) wrote the shared code for both columns and exported the mesh. Target, the weekly's: F37 Cu's — (P) green, (T5-PEC) per class against entry 31's degree-2 PEC block printed beside the ≤ 1e-6 precedent, (C-h-PEC) beside ours'. Item 454's (Q1) (this file :2057–2062) gives the shape. **No band, tolerance, Done-when or validation target is set or moved.**

**`OPS-61` step 2c, PEC half: Palace order 2 with wall attribute 401 in `Boundaries.PEC` on the h = 0.0075 exact-rectangle hole mesh at 64 MHz, one excitation per window, against entry 31's degree-2 PEC 4×4; then the Palace-internal Cu − PEC move against ours'.** Ours is not re-solved.
- **Executor:** implementer. **No probe edit, no test module, no `src/`.**
- **Tier:** heavy. Configs and read at `-n 1` in the complex container; Palace at `-np 8`. **Class:** chain step (solve), ≈ 12 measured min (the class median).
- **Time:** ≈ 32 predicted min. That is 15 fixed + W1 < 1 + W2–W5 4 × ≈ 2.4 (Palace Total ≈ 131 s plus ≈ 10 s of docker / harness) + W6 / W7 < 1 + ≈ 5 spare.
- **Price (measured basis).**
  - Item 454's four Cu windows on this mesh: Total 142.9–144.9 s, Setup 26.6–27.3 s, 39–40 GMRES iterations, peak 12.3–12.4 G (`20261015T004434Z_OPS-61-item454-W7-read.log:77–80`; e.g. `20261015T003418Z_OPS-61-item454-W3-palace-hole-h0075-cu-e1.log:159, :181, :190`). Their wall clock was ≈ 150–158 s per window (log stamps 00:34:18 / 00:36:55 / 00:39:25 / 00:41:56 / 00:44:34 UTC).
  - At ×1 the PEC column cost 0.91× the Cu column at equal iterations: 160.575 s against 175.950 s for four excitations each, 30–31 iterations in both (`20261005T021936Z_OPS-61-item234-W4-palace-hole-pec-o2.log:158–293, :324`; `20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log:161–296, :327`).
  - **Per PEC excitation ≈ 131 s (115–165)**, 35–45 iterations, peak 11–13 G. Every window ≪ 590 s.
- **Depends on:** nothing open. Inputs, checked read-only at drafting (2026-10-15 ≈ 03:00), in `logs/ops61/`:
  - `ops61_ans6_h0075_hole_strip.msh` (item 454 W1);
  - `h05x_pec_o2.json`, the template: Mesh `ops61_ans6_x1_hole_strip.msh`, `PEC.Attributes` [1, 401], no `Conductivity` key, Order 2, Freq [0.064];
  - `postpro_h05x_pec_o2/port-S.csv` (×1 Palace PEC order 2, for (C-h-PEC));
  - the four `postpro_h05x_h0075_cu_o2_e{1..4}/port-S.csv` (item 454, for (K1) / (K3) / (T3m-h));
  - absent, as required: any `h05x_h0075_pec_*.json` and any `postpro_h05x_h0075_pec_*`.
  - Tracked: `docs/testing/logs/20261012T070009Z_ANS-6-step5b.log`.

**The code it runs (item 454's, unedited; lines in `scripts/probes/ops61_step1a_palace_ans4.py` at HEAD `fdd2d03f`).**
- Constants :4510–4544. `H0075_HOLE_SETS["pec"]` = `h05x_h0075_pec_o2_e{1..4}` (:4515–4516); `H0075_HOLE_TEMPLATES["pec"] = "h05x_pec_o2.json"` (:4517); `H0075_HOLE_OURS["pec"]` (:4531–4541) — header `"[TH-14 step2] f = 6.400000e+07 Hz: PEC sweep"`, label `"S_PEC_"`, the 10-digit class means of `:4590`, `sigma_max None`, h-move 1.193e-3 / 5.01e-4 / 1.081e-3 (`:4594`), ×1 output `postpro_h05x_pec_o2`.
- `hole_h0075_configs("pec")` (:4749–4829), dispatched by `hole-h0075-configs --pec` (:5056–5057). It asserts `-n 1` and never-overwrite, parses (E1′) 226 904 type-4 tets before any write (:4767–4778), checks the template (:4780–4799; the PEC branch asserts `PEC.Attributes == [1, 401]` and no `Conductivity` key, :4797–4799), and writes four configs changing only Mesh, each port's Excitation and Output, asserting the template diff (:4804–4828).
- `read_hole_h0075("pec")` (:4883–5048), dispatched by `read --hole-h0075 --pec` (:5060–5061):
  - (Tie) on both of ours' blocks (:4911–4932). The PEC block ties class means only; its σ_max is computed and printed (:4917).
  - Palace's PEC columns assembled (:4934; `_hole_h0075_palace_s` :4856–4880), (T5-pec) per entry (:4939–4942).
  - With all four columns: (P) via `_p_check_s` (:4945–4946); (T5-pec) per class (:4948–4950); `_class_spread` side by side (:4951–4954); (C-h-pec) against `postpro_h05x_pec_o2` beside ours' `:4594` (:4955–4961); (C2-h) against `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["pec"]` (:4962–4971).
  - In PEC mode with all four Cu columns present: Palace-internal (K1) and (K3) asserted, (T3m-h) printed (:5000–5025).
  - The price section globs `*_OPS-61-item*-W*-palace-hole-h0075-pec-e*.log` (:5030–5043).
- **Two printed strings on this path were written for the Cu column.** The (C2-h) line says "predicted >= 506x"; PEC's own floor is ≥ 501× (below). The price line says "brief's projection 115-180 s"; this brief's PEC bracket of 115–165 s lies inside it. They are labels. The executor journals this disclosure and **does not edit** (no edit at all in this item).

**Anchors (asserted; bands imported, never restated).**
- **(E1′), W1:** 226 904 type-4 tets parsed from the mesh before any write (the function asserts it).
- **Template asserts, W1:** `h05x_pec_o2.json` has PEC [1, 401], no `Conductivity`, ports 21–24 `+Z` 50 Ω, Order 2, Freq [0.064]. Each written config differs from it only in Mesh / Excitation / Output.
- **(Tie), W6 — a parse check:** ours' Cu and PEC class means tie `:4589` / `:4590` within `H0075_TIE_BAND` 1e-9, and the Cu σ_max ties `:4565`. Item 454 measured 5.590e-11 on both blocks (`…item454-W7-read.log:35, :40`). A wrong block misses by ≥ 1.2e-4 (Cu against PEC, adjacent).
- **(P), W6, on Palace's assembled PEC 4×4:** item 454's set, the bands imported as that brief cites them (this file :1944–1949) — `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` with every column power ≤ 1 + 1e-9, every spread ≤ `ADJACENT_SPREAD_BAND`, pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst. Reciprocity is cross-window.
- **(K1), W6, Palace-internal — floor *asserted*:** per class, \|Palace Cu h − Palace PEC h\| / \|Palace PEC h\| > `KNOB_REACHES_THE_SOLVE` (1e-6, imported at :4891).
  - **Same-fixture prior:** item 454's (K1′) on this mesh read \|Palace Cu h − ours PEC h\| / \|ours PEC h\| = 6.155e-4 / 2.793e-4 / 4.789e-4 (`…item454-W7-read.log:75`). Ours' own Cu − PEC gap is 6.15e-4 / 2.79e-4 / 4.79e-4 (`…ANS-6-step5b.log:4592`). Both are ≥ 279× above the floor.
- **(K3), W6 — *asserted*:** σ_max(Palace Cu h) < σ_max(Palace PEC h).
  - **Same-fixture prior:** Palace Cu σ_max 0.999904940000 (`…item454-W7-read.log:67`) against ours' PEC 0.999911767900 (`:40`), margin 6.828e-6 (`:76`).
- **Palace log check (executor, journalled), W2–W5:**
  - 226 904 elements, and the bounding box ±0.12 / ±0.12 / ±0.10 m;
  - `Rs = 4.375e+01 Ω/sq` on 21–24;
  - `Configuring Dirichlet PEC BC at attributes:` followed by `1, 401` (item 234 W4's form, `20261005T021936Z_OPS-61-item234-W4-palace-hole-pec-o2.log:96–97`);
  - **no** `Configuring Robin finite conductivity BC` line and no PMC-fallback warning;
  - ω/2π = 6.400e-02 GHz; exactly one excited port, k; one `GMRES solver converged` line;
  - `ND (p = 2): 1538030` predicted, the Cu column's figure. At ×1, Palace printed the same ND for PEC (570 332, item 234 W4 `:110`) as item 454's brief used for Cu (this file :1869). A mismatch is a reading, not a red.

**Ours against Palace (printed, never gated — 2026-10-13 rulings (1) / (4); no cross-solver band exists).**
- **(T5-PEC)** per class and per entry against ours' PEC block (`:4555–4559`; class means `:4590`). **Predicted ≤ 1e-6, central ≈ 1e-10.** Precedents: item 234 (T1) ≤ 2.2e-10 on the ×1 PEC hole; item 454 ≤ 1.274e-10 on this mesh's Cu column (`…item454-W7-read.log:70`).
- **(C-h-PEC):** Palace PEC h against Palace PEC ×1 (`postpro_h05x_pec_o2`), per class, beside ours' 0.1193 / 0.0501 / 0.1081 % (`:4594`). **Predicted equal to the printed digit.** Item 454's Cu read 0.1206 / 0.0506 / 0.1091 % = ours' (`…item454-W7-read.log:73`).
- **`_class_spread`**, Palace and ours' PEC side by side. Predicted equal to the printed digit. Never compare across th14's definition (item 241's miss).
- **(T3m-h):** Δ_Pal = Palace Cu − Palace PEC against Δ_ours = ours Cu − ours PEC, per class: \|Δ_Pal − Δ_ours\| / \|Δ_ours\| and arg. **Predicted ≤ 1e-3 and \|arg\| < 0.1°, central ≈ 3e-7.** Arithmetic: \|Δ_ours\| ≈ 4.1e-4 absolute on self (6.15e-4 × 0.671) and 1.3e-4 on adjacent (2.79e-4 × 0.451). Each code's tie is ≈ 1e-10 relative, ≈ 3e-11–9e-11 absolute. That puts the ratio at ≈ 2e-7–4e-7. At ×1, item 241 measured ≤ 4.4e-7 and 0.0000° (this file :2060).
- The price per window, beside 115–165 s.

**Negative controls (rule (e)).**
- **(K1) and (K3) above** — asserted, each with its prior cited.
  - A config that kept the Cu template's `Conductivity` block reads ≈ 1e-10 against item 454's Cu columns, so (K1) goes red. The W1 template assert catches it first.
  - A wall left out of every block (PMC) reads ≥ 100 % (item 234's (K′)), and the fallback warning shows first.
- **(C2-h) mesh-mismatch separation — *predicted*, printed.** Per class: \|ours ×1 p2 PEC (`ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["pec"]`) − Palace PEC h\| / \|ours h p2 PEC − Palace PEC h\|.
  - The numerator ≈ the PEC h-move; the weakest class is adjacent, 5.01e-4 relative.
  - **Predicted ≥ 501×** if (T5-PEC) ≤ 1e-6, central ≈ 1e7×. Item 454's Cu measured 9.469e6–1.648e7× (`…item454-W7-read.log:74`). The printed "506x" is Cu's label.
  - A config left on the ×1 hole mesh reads (C-h-PEC) ≈ 0 and (T5-PEC) ≈ the h-move. The 226 904-element log check catches it first.

**Predictions (written first, printed).**
- W1: < 10 s; (E1′) 226 904.
- W2–W5, each: Total ≈ 131 s (115–165), Setup 24–28 s, GMRES 35–45 iterations (Cu: 39–40), peak 11–13 G.
- (P): reciprocity ≈ 1e-11–2e-11; σ_max 0.999911768 ± 1e-9; spreads equal to ours' `_class_spread` to the printed digit.
- (T5-PEC) ≤ 1e-6 (central ≈ 1e-10); (C-h-PEC) 0.1193 / 0.0501 / 0.1081 %; (K1) equal to ours' gap 6.15e-4 / 2.79e-4 / 4.79e-4 to the printed digits; (K3) margin ≈ 6.83e-6; (T3m-h) ≤ 1e-3, central ≈ 3e-7, \|arg\| < 0.1°; (C2-h) ≥ 501×.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; Palace v0.18.1, image `fem-em-palace:0.18.1`, compose profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit at any point.**
- **Pre-check** (not harness calls):
  - `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up and no lingering `fem-em-palace`. If an `xl` window is still computing (tonight's entry 34, F31, from 02:00), start nothing: journal and stop.
  - `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
  - `ls logs/ops61/`: the four inputs above are present. Stop and report if any `h05x_h0075_pec_*.json` or `postpro_h05x_h0075_pec_*` exists.
  - The probe dispatches `hole-h0075-configs --pec` (:5056). If not, stop: the brief's premise failed.
- **W1** (< 10 s), the PEC configs and (E1′):

  `scripts/testing/run_and_log.sh OPS-61-item457-W1-hole-h0075-configs-pec "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py hole-h0075-configs --pec'"`
- **W2** (≈ 131 s, ≤ 590), excitation 1, the price check:

  `scripts/testing/run_and_log.sh OPS-61-item457-W2-palace-hole-h0075-pec-e1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 h05x_h0075_pec_o2_e1.json'"`
- **Price gate (binding, after W2; item 454's gate).** Run W3–W5 iff W2's Palace `Total` (Avg) T₁ ≤ 300 s. Otherwise columns 2–4 are "not measured — price"; W6 still runs ((T5-pec) per entry on column 1; (P) not measured). The item lands as a partial record, not BLOCKED.
- **W3 / W4 / W5:** W2's command with `…_pec_o2_e2.json` / `…_e3.json` / `…_e4.json`, labelled `OPS-61-item457-W3-palace-hole-h0075-pec-e2` / `-W4-…-pec-e3` / `-W5-…-pec-e4`.
  - **The `-palace-hole-h0075-pec-e` label stem is load-bearing:** the read's price glob (:5030–5031) finds nothing else.
  - **Minute gate:** if W5 cannot start by minute 50, stop launching. The next slot runs the missing excitations with no edit and then W6 (rule (d), "partial; columns k… pending").
  - After any kill, the `pgrep` check first.
- **W6** (≈ 3 s), the read:

  `scripts/testing/run_and_log.sh OPS-61-item457-W6-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --hole-h0075 --pec'"`
- **W7** (≈ 3 s), the regression: item 454's W7 command verbatim (`20261015T004434Z_OPS-61-item454-W7-read.log:12`), labelled `OPS-61-item457-W7-regress-read-hole-cu`. Its output section must `diff` empty against that log's. The Cu read never opens a PEC output, so the new files must not move it.

**Traps already paid for.** All of items 234 / 241 / 418 / 434 / 444 / 454's apply, plus:
- **The harness call is the whole Bash command**: bare, repo-relative, nothing chained or piped. An executor that loses the docker socket mid-item stops and reports.
- **Never run** `export --hole-h0075` (the mesh exists; the function asserts absence), `hole-h0075-configs --cu` (its configs exist), or any other `export` / `*-configs` mode. Write nothing but the four new PEC configs; never overwrite `h05x_*`, `f05x_*`, `t39_*`–`t44_*` or their outputs.
- **Attribute 401 goes in exactly one block** — `PEC` here, with no `Conductivity`. The template assert guards it. A PEC log showing a `Robin finite conductivity` line means the wrong template: stop.
- **Parse trap (entry 31's log).** It holds one `S_PEC_1k =` block (`:4556`) under the `[TH-14 step2] … PEC sweep` header (`:4555`) and three `S_1k =` blocks (Cu `:4561`, 5.8e9 `:4568`, 5.8e11 `:4575`). The read anchors on the header prefixes; (Tie) catches a wrong block. Ours' PEC σ_max is not printed at 10 digits in that log: it is computed (0.999911767900, `…item454-W7-read.log:40`), never a literal.
- **One excitation per Palace command.** Never raise `-np` or lengthen a timeout in-slot. A Status 124 is a price finding. A `GMRES … did not converge` line is a finding: that column is "not measured".
- **Orphans after any killed window:** `pgrep -c python3` in the container reads 0, and `docker compose -f docker/docker-compose.yml ps` shows no lingering `fem-em-palace`.
- **Spread definitions differ:** `_class_spread` against `_class_spread` only.
- **Compare S only.** Never read `z_matrix_ohm` or the `aed` key. Never register a record from Palace. No Palace figure is set against an AED one in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered: Palace order 2 with PEC on 401 on the h = 0.0075 exact hole rung at 64 MHz; (P) on the PEC 4×4; (T5-PEC) / (C-h-PEC) against entry 31; Palace-internal (K1) / (K3) and (T3m-h). This completes F37. Not covered:
- 10 / 128 MHz on the rung;
- the h = 0.005 hole rung (F31, `xl` entry 34; T41 reads it);
- order 1, the σ-scaling control, Palace's surface-loss postprocessing;
- any Palace h-statement (two hole meshes are not a ladder);
- the tuned twin; AED; `src/`; any band;
- F38's public tables (lane 1, serial on this item).

**Status it can move.**
- The `OPS-61` row (✅, PROJECT_PLAN.md:1351): "the h = 0.0075 copper hole rung (F37 Cu, item 454)" becomes "the h = 0.0075 hole rung, PEC and copper (F37, items 454 / 457)". **No glyph change.** Rule (j): the row's Done-when (step 1c's `f05x` solve) is met; a record clause needs no amendment.
- §10 F37 is recorded DONE (both halves). F38 (lane 1, "serial on F37 PEC", PROJECT_PLAN.md:4706–4710) becomes queueable, and F40's sentence (PROJECT_PLAN.md:4718–4721, "×1, 0.0075 hole") gains its PEC half.
- If (T3m-h) ≤ 1e-3 with (T5-PEC) ≤ 1e-6, the Leontovich term's own contribution — not only the Cu 4×4 — is matched across codes on a refined mesh. Extending the §2 / CLAUDE.md Palace clause is the weekly's call.

**Negative result** (no `PREDICTED RED:` line; nothing here is expected red; any red stops the slot).
- **(E1′) or a W1 template assert red:** the wrong mesh or template. Stop before any Palace window, report the count or keys, BLOCKED (rule (d)). No re-export in-slot.
- **Palace aborts, or a W2–W5 log lacks 226 904 elements / PEC on `1, 401` / 6.400e-02 GHz, or shows a Robin conductivity line or the PMC-fallback warning:** record it verbatim, open a known-issues entry, park, BLOCKED.
- **(Tie) red:** a parse failure on a log item 454 tied at 5.59e-11 — the log changed (rotation to `.log.gz` is already handled, :4840–4842). This item may not edit the probe: stop, BLOCKED; the review licenses any fix.
- **(P) red:** report every identity's number, park, BLOCKED. Never widen.
- **(K1) red** (Palace Cu within 1e-6 of Palace PEC): either the PEC configs carry a conductivity block or item 454's Cu outputs were overwritten — it contradicts item 454's green (K1′). Park, BLOCKED.
- **(K3) red:** a loss-direction finding. Report it, park, BLOCKED; the review rules.
- **(T5-PEC) ≫ 1e-6, or (C-h-PEC) ≠ ours beyond the print, with (P), (Tie), (K1) and (K3) green:** a reading, not a failure. DONE as a record with a known-issues entry. PEC tied at ≤ 2.2e-10 at ×1 and Cu at ≤ 1.3e-10 on this mesh, so it would isolate a PEC-elimination difference on the refined wall.
- **(T3m-h) > 1e-3 with everything else green:** a reading (printed, never gated); DONE as a record, for the review.
- **W7 not diff-empty:** an input moved without an edit. Stop, report the diff, BLOCKED.
- **Price gate or minute gate trips:** a partial record, not BLOCKED.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 457 per `docs/planning/chunks/OPS-61.md` (read item 454's brief, its outcome and the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; the configs and the read at `-n 1`, one Palace excitation per window at `-np 8`, never two heavy windows at once; the price gate after W2 and the minute-50 gate bind. No probe edit at all — item 454 wrote the PEC dispatch; never run `export --hole-h0075`, `hole-h0075-configs --cu` or any older `export` / `*-configs` mode; never overwrite `h05x_*`, `f05x_*` or `t39_*`–`t44_*`; never open `docs/private/` or `aed_results/`."

### §9 item 458 — 2026-10-15 (§10 chain F, F39 first half: code, configs, excitations 1–2; chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.* Licence: §9 rule (4) (PROJECT_PLAN.md:3805–3813). F39 is a numbered step the 2026-10-14 weekly enumerated in §10 chain F (PROJECT_PLAN.md:4711–4717) and tiered in its rulings (:4627–4632: "64 MHz **heavy** (F39 — four single-excitation windows, predicted 400–430 s each by F34's 10 → 64 MHz ratio 1.24–1.27×); 128 MHz **declined**"). Serial on nothing. **Not provisional.** Target, the weekly's (F36's): (P) green; (T5) per class against `ANS-4` step 6c's 64 MHz degree-2 rung; (C-h) equal to ours' 0.0075 → 0.005 move at 64 MHz (0.0569 / 0.0235 / 0.0611 %) to the printed digit. **If any window exceeds 590 s the step stops and records the price.** **No band, tolerance, Done-when or validation target is set or moved.**

**Split (the minute-50 gate forces it).** One item would cost ≈ 65 min:
- 15 fixed;
- ≈ 15 code (three additive functions, constants, two dispatch branches; item 437 budgeted ≈ 10, item 454 ≈ 15);
- W1 configs ≈ 1;
- four excitation windows at ≈ 7 min wall each (Palace Total 380–435 s plus ≈ 9 s of harness / docker — item 437's W2 → W3 stamps 21:56:59 → 22:02:42 are 343 s of wall for a 334.2 s Total);
- the read and two regressions ≈ 1, and ≈ 5 to journal and commit.

Excitation 4 would start at ≈ minute 52 (15 + 15 + 1 + 3 × 7), past the minute-50 gate. With excitations 1–2 only, e2 starts at ≈ minute 38 (pessimistic 41) and the item closes at ≈ 52 min. So: **item 458 = the code, the configs, e1 (the price probe), e2, a partial read and two dispatch regressions; item 458b = e3, e4 and the full read, with no edit.** It is one chain step in two slots. The `postpro_*` outputs are durable (item 454 / 457's precedent for halves).

**`OPS-61` step 1j, first half: Palace order 2 on `ANS-4`'s h = 0.005 exact-strip rung (the mesh item 437 exported) at 64 MHz, excitations 1 and 2, against our degree-2 record of that rung.** Ours is not re-solved: it is read from `xl` entry 30's tracked log.
- **Executor:** implementer. Probe only (`scripts/probes/ops61_step1a_palace_ans4.py`), additive and disclosed; no test module, no `src/`.
- **Tier:** heavy. Configs and read at `-n 1` in the complex container; Palace at `-np 8`. **Class:** chain step (solve), ≈ 12 measured min.
- **Time:** ≈ 52 predicted min (above).
- **Price (from measurement).**
  - **This rung at 10 MHz (item 437):** Total 304.2–334.2 s, Setup 94.0–121.9 s, 25 GMRES iterations, peak 32.2–32.8 G (`20261013T221851Z_OPS-61-item437-W6-read.log:87–90`; outcome this file :1524).
  - **The 0.0075 strip rung, 10 → 64 MHz (items 418 / 434):** Total 143.4–146.2 → 181–186 s (×1.24–1.27), Setup ≈ 33–34 s at both, iterations 25 → 43 (this file :906; item 434 e1 `20261013T170207Z_OPS-61-item434-W2-palace-h0075-f64-e1.log:101, :159, :181, :190`, as item 454's brief cites it at :1870).
  - **Model.** Fit non-Setup = a + b·iterations on the 0.0075 pair: 112.0 s at 25 and 152.9 s at 43 give b ≈ 2.27 s per iteration and a ≈ 55 s. Scale both by the 0.005 / 0.0075 non-Setup ratio at 10 MHz (207.8–228.1 / 112.0 = 1.86–2.04). At 43 iterations, non-Setup ≈ 284–312 s; plus Setup 94–122 s, **Total ≈ 380–435 s, central ≈ 410 s.** That matches item 437's printed projection, 392.1–429.3 s (`…item437-W6-read.log:87–90`), and the weekly's 400–430 s.
  - **Envelopes:** at 50 iterations in the same model, ≈ 410–465 s. With non-Setup fully proportional to iterations (43 / 25), ≈ 450–515 s. Every bracket is < 590 s.
  - **Memory:** peak ≈ 33–34 G (10 MHz 32.2–32.8 G × the 0.0075 rung's 64 / 10 MHz peak ratio 14.8 / 14.4), against the service's 128 G.
- **Depends on:** nothing. Inputs, checked read-only at drafting (2026-10-15 ≈ 03:00), in `logs/ops61/`:
  - `ops61_ans4_h0050_strip.msh` (item 437 W1);
  - `f05x_h0050_o2_f10_e{1..4}.json` (templates);
  - `postpro_f05x_h0075_o2_f64_e{1..4}/port-S.csv` (item 434, for (C-h));
  - absent, as required: any `f05x_h0050_o2_f64_*.json` and any `postpro_f05x_h0050_o2_f64_*`.
  - Tracked: `docs/testing/logs/20261011T070008Z_ANS-4-step6c.log` (ours, h = 0.005) and `docs/testing/logs/20261005T070009Z_ANS-4-step6b.log` (ours, h = 0.0075).

**What exists on our side (the comparand; not re-solved).** `xl` entry 30, `20261011T070008Z_ANS-4-step6c.log`: degree 2, `-n 16`, κ off, 593 478 cells (`:3814`).
- **The 64 MHz block:** header `[ANS-4 step5] 64 MHz kappa-off S (10 digits):` at `:3821`; `S_1k` at `:3822`, `S_2k` at `:3823`; **`:3824` is blank; `S_3k` sits at `:3825` behind a pytest node-ID prefix (`tests/validation/test_ans4_exact_strip_records.py::test_g_port11_gates_on_the_exact_rectangles[10]` plus spaces); `:3826` is blank; `S_4k` sits at `:3827` behind the same prefix.** Z-class spreads, σ_max 0.999823581365 and reciprocity at `:3828`; the (G) line at `:3907` repeats σ_max.
- Ours' ×1 → 0.005 move at 64 MHz: 0.2868 / 0.1144 / 0.2699 % (`:3841`), against `ANS4_EXACT_STRIP_S_P2_X1["64"]`.
- **Ours at h = 0.0075, 64 MHz:** `20261005T070009Z_ANS-4-step6b.log`, header `:3708`, rows `:3709–3712` (contiguous; `:3709` PASSED-prefixed), σ_max 0.999823496773 at `:3713`. `_ours_rung_s` already parses it through `H0075_OURS_BLOCKS["64"]` (probe :730–741; item 434).
- **Ours' 0.0075 → 0.005 move at 64 MHz,** recomputed by this drafter from the two 10-digit blocks: self 0.05689 % (3.491e-4 absolute), adjacent 0.02353 % (9.851e-5), opposite 0.06113 % (2.020e-4). That is the weekly's 0.0569 / 0.0235 / 0.0611 %. The read recomputes it from the parsed logs, never from a literal.

**The parse trap this item must not buy.** `_ours_rung_s` (:2268–2288) takes the four lines after the header and asserts that each, after stripping leading `PASSED`s, starts with `S_{k}k =`. On step 6c's 64 MHz block, line `:3824` is blank, so it raises `AssertionError: (3, '')` before any reading. **`_ours_rung_s` is not edited** — items 418 / 434 / 435 / 437's reads depend on it. A tolerant copy is added (below) and proven against it.

**Code (probe only; additive; disclosed).**
- **Pre-check for the code:** at HEAD `fdd2d03f` only `export --exact-strip-h0050` (:5074) and `read --exact-strip-h0050` (:5076) are dispatched for this rung. If a 64 MHz `h0050` mode already exists, make no edit and report.
- **Constants** beside `H0050_*` (:759–771), each with a provenance comment:
  - `H0050_F64_SET = {f"f05x_h0050_o2_f64_e{k}": k for k in (1, 2, 3, 4)}`;
  - `H0050_F64_OURS_HEADER = "[ANS-4 step5] 64 MHz kappa-off S (10 digits):"` (step 6c `:3821`);
  - `H0050_F64_OURS_SIGMA_MAX = 0.999823581365` (`:3828`, `:3907`; the (Tie) target);
  - `H0050_F64_OURS_X1_MOVE = {"self": 2.868e-3, "adjacent": 1.144e-3, "opposite": 2.699e-3}` (`:3841`);
  - `H0050_F64_PROJ = (380.0, 435.0)` (this brief's bracket; item 437 W6 printed 392.1–429.3).

  `H0050_MESH_NAME`, `H0050_CELL_RECORD`, `H0050_OURS_LOG`, `H0050_OURS_HEADER`, `H0075_FREQ_SET`, `H0075_OURS_LOG`, `H0075_OURS_BLOCKS` and `H0075_TIE_BAND` are **reused, not copied**.
- **New `h0050_f64_configs()`, dispatched by `h0050-configs --f64`.** A copy of `h0075_freq_configs` (:2202–2265; provenance comment with commit and span) with `H0050_MESH_NAME`, `H0050_CELL_RECORD`, `H0050_F64_SET`, templates `f05x_h0050_o2_f10_e{k}.json` and Freq → [0.064].
  - Same steps as its source: assert `-n 1`; never overwrite; parse (E1′) the type-4 count == 593 478 before any write.
  - Per template, assert Mesh == `H0050_MESH_NAME`, ports 21–24 `+Z` 50 Ω, Excitation true on k only, Order 2, Freq [0.010].
  - Write changing only `Freq` and `Problem.Output`; read back and assert the template diff is exactly those two keys.
  - **Drop the source's sha256 / byte-size print (:2227–2229)**; print element and tet counts only (the private-leak hook trap, items 434 / 435; this file :1905).
- **New `_ours_rung_s_scan(log, header)`.** A copy of `_ours_rung_s` that:
  1. requires exactly one stripped-equal header;
  2. walks forward at most 8 lines, skipping blank lines;
  3. on each non-blank line finds `S_{k}k =` by `re.search`, so a pytest node-ID or `PASSED` prefix is tolerated, and takes the text from that match, asserting 4 complex tokens;
  4. asserts k runs 1, 2, 3, 4 in order, and that no line beginning `[` (a new tagged header) is met before row 4;
  5. returns (array, header line number). Like its source, it opens `.log.gz` if the log was rotated.
- **New `read_exact_strip_h0050_f64()`, dispatched by `read --exact-strip-h0050 --f64`.** A copy of `read_exact_strip_h0050` (:3407 to the end of that function; provenance comment) with:
  - **(Parse-0), asserted first:** `_ours_rung_s_scan(H0050_OURS_LOG, H0050_OURS_HEADER)` equals `_ours_rung_s(H0050_OURS_LOG, H0050_OURS_HEADER)` bitwise on the contiguous 10 MHz block (`:3815–3819`) — the new parser reproduces the old one where both apply;
  - ours at h = 0.005 from `_ours_rung_s_scan(H0050_OURS_LOG, H0050_F64_OURS_HEADER)`, with (Tie) σ_max against `H0050_F64_OURS_SIGMA_MAX` within `H0075_TIE_BAND`;
  - ours at h = 0.0075, 64 MHz, from `_ours_rung_s(H0075_OURS_LOG, H0075_OURS_BLOCKS["64"]["header"])`, with (Tie75) σ_max against `H0075_OURS_BLOCKS["64"]["sigma_max"]` within the same band;
  - Palace assembled from `H0050_F64_SET` (h = 0.005) and `H0075_FREQ_SET["64"]` (h = 0.0075, item 434's, for (C-h)) through the copied `assemble`;
  - (T5) per entry on the present columns;
  - with four columns: (P) via `_p_check_s`; (T5) per class; `_class_spread` for Palace and ours side by side (an addition: item 437's read did not print ours'); ours' 0.0075 → 0.005 move computed from the two parsed blocks, printed beside the drafter's 0.0569 / 0.0235 / 0.0611 %; (C-h) Palace h = 0.005 against Palace h = 0.0075 beside it; (C2-h) against ours at h = 0.0075 and against `ANS4_EXACT_STRIP_S_P2_X1["64"]`;
  - a price section globbing `*_OPS-61-item458*-W*-palace-h0050-f64-e*.log` (both halves), printing Total / Setup / iterations / peak beside `H0050_F64_PROJ`. No 128 MHz projection (declined).

  `read_exact_strip_h0050`, `read_exact_strip_h0075_freq`, `_ours_rung_s` and `h0075_freq_configs` are not edited. The read never touches `z_matrix_ohm`, the `aed` key, `docs/private/` or `aed_results/`.
- **Dispatch:** two new exact-list `sys.argv` branches, first in `__main__`, as items 444 / 448 / 454 placed theirs. W5 / W6 prove every existing mode unchanged.

**Anchors (asserted; bands imported, never restated).**
- **(E1′), W1:** 593 478 type-4 tets parsed before any write; the template diff is Freq / Output only.
- **(Parse-0), W4:** bitwise equality on the 10 MHz block.
- **(Tie), W4 — a parse check:** the parsed 64 MHz h = 0.005 block's σ_max against 0.999823581365 within 1e-9. Wrong blocks miss by much more: step 6b's 64 MHz block (0.999823496773) by 8.46e-8, 85× the band; the 10 MHz block by 1.7e-4; the 128 MHz block (`:3834`, 0.999320969630) by 5.0e-4. **(Tie75)** as item 434 asserted it.
- **(P)** needs four columns: item 458b.
- **Palace log check (executor, journalled), W2 / W3:**
  - 593 478 elements, and the box ±0.12 / ±0.12 / ±0.10 m;
  - `Rs = 4.375e+01 Ω/sq` on 21–24; PEC on 1 only;
  - ω/2π = **6.400e-02** GHz;
  - exactly one excited port, k; one `GMRES solver converged` line;
  - `ND (p = 2): 3833208` predicted (= ours, item 437; a mismatch is a reading).

**Ours against Palace (printed, never gated).** (T5) per entry on columns 1–2. **Predicted ≤ 1e-6, central ≈ 1e-10.** Precedents: item 434 ≤ 8.0e-11 per class at 0.0075 and 64 MHz; item 437 ≤ 1.68e-10 per entry at 0.005 and 10 MHz.

**Negative controls (rule (e)).**
- **(Tie)'s wrong-block margins — *asserted*.** The backing is the printed σ_max values in the same two logs (`step6c.log:3820, :3828, :3834`; `step6b.log:3713`).
- **(Parse-0) — *asserted*.** A parser that mis-assigned or dropped rows fails bitwise equality on the contiguous block. It is identity by construction on the same log, so it is a parser control, not physics.
- **A config left at 10 MHz:** the log shows ω/2π = 1.000e-02 and (T5) reads O(1). The log check catches it first (*predicted*).
- **A wrong rung:** (T5) per entry ≈ the 0.0075 → 0.005 move (≥ 2.35e-4 relative on class means) — *predicted*. The 593 478-element check catches it first.
- **(C2-h)** is item 458b's reading — *predicted* ≥ 235× against ours at h = 0.0075 (weakest adjacent 2.35e-4 / 1e-6), central ≈ 2e6× (9.85e-5 / ≈ 4e-11 absolute), and ≥ 1 144× against ours ×1 (weakest adjacent 1.144e-3, `:3841`).

**Predictions (written first, printed).**
- W1: < 30 s; (E1′) 593 478.
- W2 / W3, each: Total ≈ 410 s (380–435; envelope to 515), Setup 90–125 s, GMRES 40–50 iterations (0.0075 at 64 MHz: 43), peak 32–35 G.
- (Parse-0) bitwise; (Tie) \|diff\| ≲ 5e-11 (item 437 measured 4.837e-11); (T5) per entry ≤ 1e-6, central ≈ 1e-10.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`; `-np` ≤ 12; Palace v0.18.1, image `fem-em-palace:0.18.1`, profile `palace`, one-shot `run --rm`). Never two heavy windows at once. **No probe edit after W1 starts** (rule (i)).
- **Pre-check** (not harness calls):
  - `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up and no lingering `fem-em-palace`. If an `xl` window is still computing, start nothing: journal and stop.
  - `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
  - `ls logs/ops61/`: the inputs above are present. Stop if any `f05x_h0050_o2_f64_*` or `postpro_f05x_h0050_o2_f64_*` exists.
  - If item 457 ran earlier in the same slot, its last window has returned.
- **W1** (< 30 s), the configs and (E1′):

  `scripts/testing/run_and_log.sh OPS-61-item458-W1-h0050-configs-f64 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py h0050-configs --f64'"`
- **W2** (≈ 410 s, ≤ 590), excitation 1, **the price probe**:

  `scripts/testing/run_and_log.sh OPS-61-item458-W2-palace-h0050-f64-e1 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_h0050_o2_f64_e1.json'"`
- **Price gate (binding, after W2; the weekly's rule).**
  - **Status 124 (T₁ > 590 s):** the step stops and records the price — "> 590 s", the last iteration line, no Elapsed Time Report. W3 is not run. W4 still runs (columns present []), and so do W5 / W6. The item lands as a partial record, not BLOCKED. Item 458b is withdrawn, and F39 returns to the weekly as an `xl` question.
  - **Status 0 with T₁ in (435, 590]:** a reading outside the bracket. Continue and print it.
- **W3** (≈ 410 s), excitation 2: W2's command with `…_e2.json`, labelled `OPS-61-item458-W3-palace-h0050-f64-e2`. **Minute gate:** if W3 cannot start by minute 50, do not launch it; item 458b then runs e2–e4.
- **W4** (≈ 3 s), the partial read:

  `scripts/testing/run_and_log.sh OPS-61-item458-W4-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-h0050 --f64'"`
- **W5** (≈ 3 s), dispatch regression 1: item 437's W6 command verbatim (`20261013T221851Z_OPS-61-item437-W6-read.log:12`), labelled `OPS-61-item458-W5-regress-read-h0050-f10`. Its output section must `diff` empty against that log's.
- **W6** (≈ 3 s), dispatch regression 2: item 434's W6 command verbatim (`20261013T171447Z_OPS-61-item434-W6-read.log:12`), labelled `OPS-61-item458-W6-regress-read-h0075-f64`. Its output section must `diff` empty against that log's.
- **The label stems are load-bearing.** Item 437's price glob is `item437-W[2-5]-palace-h0050-e*` (:3529–3530) and item 434's is `item434-W[2-5]-palace-h0075-f64-e*` (:2386–2387). `item458…-palace-h0050-f64-e*` collides with neither, so W5 / W6 stay diff-empty, and the new glob finds both halves.

**Traps already paid for.** All of items 418 / 434 / 437 / 454's apply, plus:
- **The harness call is the whole Bash command**: bare, repo-relative, nothing chained or piped. An executor that loses the docker socket mid-item stops and reports.
- **Never run** `export --exact-strip-h0050` (the mesh exists; that function also writes the 10 MHz configs) or any older `export` / `*-configs` mode. Never overwrite `f05x_h0050_o2_f10_*`, `f05x_h0075_*`, `h05x_*`, `t39_*`–`t44_*` or their outputs.
- **The step 6c parse trap** (above): a blank line and a pytest node-ID prefix inside the 64 MHz block. Never edit `_ours_rung_s`.
- **One excitation per Palace command.** Never raise `-np` (no `-np 12` Palace point exists) or lengthen a timeout in-slot. A `GMRES … did not converge` line is a finding: that column is "not measured".
- **Orphans after any killed window:** `pgrep -c python3` in the container reads 0, and `docker compose -f docker/docker-compose.yml ps` shows no lingering `fem-em-palace`.
- **Spread definitions differ:** the log's Z-class spreads (`:3828`) are not `_class_spread`.
- **Compare S only.** Never register a record from Palace. No Palace figure is set against an AED one in any tracked file. Never open `docs/private/`, `aed_results/` or `COMPARISON_private.md`.

**Scope.** Covered, across both halves: Palace order 2 at 64 MHz on the h = 0.005 exact-strip rung; (P) on the assembled 4×4; (T5) / (C-h) / (C2-h) against entry 30's record; the measured 64 MHz price on this rung. **This half:** the code, the configs, columns 1–2 with their per-entry (T5), and the price. Not covered:
- 128 MHz on this rung (declined, PROJECT_PLAN.md:4631–4632);
- order 1; a Palace h-statement or extrapolation;
- the hole / copper rungs; the private reading; `src/`; any band.

**Status it can move.** By itself: the §10 F39 line gains "e1–e2 measured, T₁ … s". That confirms or refutes the weekly's 400–430 s prediction, the measurement its tier ruling rests on. A Status 124 instead returns F39 to the weekly as an `xl` question. The §7 clause lands with item 458b.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(E1′) or template-diff red:** stop before any Palace window, revert the probe, BLOCKED (rule (d)).
- **(Parse-0) or (Tie) red:** the parser. Fix it once in-slot, in the read's own code, then re-run W4–W6 on the final file (rule (i)). Still red: BLOCKED.
- **W2 Status 124:** the price finding, as above — a partial record, not BLOCKED.
- **Palace aborts, or a W2 / W3 log lacks 593 478 elements / Rs 43.75 Ω/sq / 6.400e-02 GHz:** record it verbatim, open a known-issues entry, BLOCKED.
- **W5 or W6 not diff-empty:** revert the probe, BLOCKED.
- **(T5) per entry ≫ 1e-6 with (Tie) green:** a reading. DONE as a partial record; item 458b still runs, because (P) and (C-h) need four columns.

**Outcome (2026-10-15 09:30 slot, implementer): DONE, first half.** Code additive at the probe (`h0050_f64_configs`, `_ours_rung_s_scan`, `read_exact_strip_h0050_f64`, five `H0050_F64_*` constants, two dispatch branches first in `__main__`); `_ours_rung_s`, `read_exact_strip_h0050`, `h0075_freq_configs` not edited. W1 `20261015T095500Z_OPS-61-item458-W1-h0050-configs-f64.log`: (E1′) 593 478 PASS, template diff Freq / Output only on all four, 1 s. **W2 (e1, the price probe) `20261015T095507Z_…-W2-palace-h0050-f64-e1.log`: Status 0, Palace Total 419.1 s (harness Elapsed 421 s), Setup 98.0 s, 44 GMRES its, peak 33.2 G — within 590 s and inside the 380–435 s bracket, so item 458b is not withdrawn.** W3 (e2) `20261015T100219Z_…-W3-palace-h0050-f64-e2.log`: Status 0, Total 433.6 s (435 s), Setup 102.3 s, 45 its, peak 33.9 G. Log checks on both: 593 478 elements, box ±0.12 / ±0.12 / ±0.10 m, Rs 4.375e+01 Ω/sq on 21–24, PEC on 1 only, ω/2π = 6.400e-02 GHz, one excited port (k), one `GMRES solver converged` line, `ND (p = 2): 3833208` (= predicted). W4 `20261015T100941Z_OPS-61-item458-W4-read.log` (2 s): (Parse-0) PASS bitwise on the 10 MHz block (`:3815`); (Tie) σ_max 0.999823581396 vs 0.999823581365, |diff| 3.091e-11 (predicted ≲ 5e-11); (Tie75) 1.961e-11; (T5) per entry, column 1: 1.428e-10 / 4.905e-11 / 1.098e-10 / 1.098e-10, column 2: 1.140e-10 / 1.254e-10 / 9.192e-11 / 1.299e-10 (predicted ≤ 1e-6, central ≈ 1e-10). W5 `20261015T100950Z_…-W5-regress-read-h0050-f10.log` and W6 `20261015T100953Z_…-W6-regress-read-h0075-f64.log` diff-empty (on `[OPS-61` / `ours` / `Palace` / Status lines) against item 437 W6 and item 434 W6. Item 458b runs e3–e4 and the full read with no edit.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 458 per `docs/planning/chunks/OPS-61.md` (read item 437's brief and outcome and the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; configs and read at `-n 1`, one Palace excitation per window at `-np 8`, never two heavy windows at once; excitations 1 and 2 only; the price gate after W2 and the minute-50 gate bind. Probe only, additive — never edit `_ours_rung_s`, `read_exact_strip_h0050` or `h0075_freq_configs`; never run `export --exact-strip-h0050` or any older `export` / `*-configs` mode; never overwrite `f05x_*`, `h05x_*` or `t39_*`–`t44_*`; never open `docs/private/` or `aed_results/`."

### §9 item 458b — 2026-10-15 (§10 chain F, F39 second half: excitations 3–4 and the full read; chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.* Licence and target: item 458's (§9 rule (4); F39, PROJECT_PLAN.md:4711–4717). It is the second slot of the same step, not a new step. The review may renumber it.

**`OPS-61` step 1j, second half: excitations 3 and 4 on the h = 0.005 exact-strip rung at 64 MHz, then the full read.** No edit of any kind: every function it runs is the one item 458's W4–W6 ran, as committed (rule (i)).
- **Executor:** implementer. No probe edit, no test module, no `src/`.
- **Tier:** heavy, Palace `-np 8`; read at `-n 1`. **Class:** chain step (solve), ≈ 12 measured min.
- **Time:** ≈ 36 predicted min: 15 fixed + 2 × ≈ 7 (item 458's measured T₁ / T₂ replace this brief's 380–435 s) + read < 1 + ≈ 5 spare. If item 458's minute gate left e2 unrun, add ≈ 7 (≈ 43 min).
- **Depends on:** item 458 — its code at HEAD with W4–W6 green, the four `f05x_h0050_o2_f64_e*.json` on disk, and `postpro_f05x_h0050_o2_f64_e1` (and `_e2`) holding `port-S.csv`. **Withdrawn** if item 458 recorded W2 Status 124.

**Anchors (asserted; bands imported, never restated).**
- (Parse-0), (Tie) and (Tie75), re-asserted by the read, as in item 458.
- **(P) on Palace's assembled 4×4:** item 454's set, the bands imported as that brief cites them (this file :1944–1949). Reciprocity is cross-window, across two slots.
- **The Palace log check** for e3 / e4, as in item 458 (593 478 elements, 6.400e-02 GHz, one excited port, `ND (p = 2): 3833208`).

**Ours against Palace (printed, never gated).**
- **(T5)** per class and per entry against step 6c's 64 MHz degree-2 block (`:3821–3827`). **Predicted ≤ 1e-6, central ≈ 1e-10.**
- **(C-h):** Palace h = 0.005 against Palace h = 0.0075 at 64 MHz (item 434's `postpro_f05x_h0075_o2_f64_e*`), per class, beside ours' move computed from the two parsed blocks. **Predicted 0.0569 / 0.0235 / 0.0611 %, equal to ours' to the printed digit.** Item 437 read exactly that at 10 MHz, and item 434 at the 0.0075 rung.
- `_class_spread` for Palace and ours, side by side; predicted equal to the printed digit.
- The price for all four excitations, beside 380–435 s.

**Negative controls (rule (e)).**
- (Tie)'s wrong-block margins and (Parse-0) — *asserted*, as in item 458.
- **(C2-h) — *predicted*, printed:** ≥ 235× against ours at h = 0.0075, central ≈ 2e6×; ≥ 1 144× against ours ×1. A config left on the 0.0075 mesh reads ≈ 1 here, and the element check catches it first.
- (P)'s power to discriminate was shown by item 159's flip; not re-run.

**Predictions (written first, printed).**
- e3 / e4, each: Total within ±10 % of item 458's T₁ (central ≈ 410 s).
- (P): reciprocity ≈ 1e-11–3e-11; σ_max 0.999823581 ± 1e-9; spreads ≲ 0.05 % (item 434's 0.0075 rung read ≤ 0.0172 %).
- (T5) ≤ 1e-6; (C-h) to the printed digit; (C2-h) ≥ 235×.

**Windows** (as item 458's conventions). **No edit at any point.**
- **Pre-check:**
  - `docker compose -f docker/docker-compose.yml ps` shows `fem-em-solver` Up, no lingering `fem-em-palace` and no `xl` window computing.
  - `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0.
  - `ls logs/ops61/`: the four `f05x_h0050_o2_f64_e*.json` are present, `postpro_…_f64_e1` holds `port-S.csv`, and no `postpro_…_f64_e3` / `_e4` exists.
  - The probe dispatches `read --exact-strip-h0050 --f64`.
- **W0** (only if item 458 left e2 unrun), labelled `OPS-61-item458b-W0-palace-h0050-f64-e2`:

  `scripts/testing/run_and_log.sh OPS-61-item458b-W0-palace-h0050-f64-e2 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_h0050_o2_f64_e2.json'"`
- **W1** (≈ 410 s, ≤ 590), excitation 3:

  `scripts/testing/run_and_log.sh OPS-61-item458b-W1-palace-h0050-f64-e3 "docker compose -f docker/docker-compose.yml --profile palace run --rm -T fem-em-palace bash -lc 'cd /workspace/logs/ops61 && timeout -k 30 590 palace -np 8 f05x_h0050_o2_f64_e3.json'"`
- **W2** (≈ 410 s), excitation 4: W1's command with `…_e4.json`, labelled `OPS-61-item458b-W2-palace-h0050-f64-e4`.
  - **Minute gate:** W2 must start by minute 50; otherwise stop launching, and the next slot runs it with no edit, then W3 (rule (d)).
  - **Any Status 124:** the step stops and records the price (the weekly's rule). No further excitation window; W3 runs on what is present; a partial record, not BLOCKED.
- **W3** (≈ 3 s), the full read:

  `scripts/testing/run_and_log.sh OPS-61-item458b-W3-read "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace timeout -k 30 300 mpiexec -n 1 python3 scripts/probes/ops61_step1a_palace_ans4.py read --exact-strip-h0050 --f64'"`

**Traps already paid for.** Item 458's list; plus, **no edit of any kind** — a parse or read defect found here stops the item (BLOCKED) for the review to license, because the read was closed on item 458's committed file.

**Scope.** As item 458's, covered across both halves. Not covered: 128 MHz (declined); order 1; any Palace h-statement; the hole rungs; `src/`; any band.

**Status it can move.**
- The `OPS-61` row (✅, PROJECT_PLAN.md:1351): "the h = 0.005 strip rung at 10 MHz (F36, item 437)" becomes "the h = 0.005 strip rung at 10 / 64 MHz (F36 / F39, items 437 / 458 / 458b)". **No glyph change.** Rule (j): the row's Done-when (step 1c's) is met; a record clause needs no amendment.
- §10 F39 is recorded DONE. F38's `ans:4` half (PROJECT_PLAN.md:4709–4710, "gains F34 / F36 the same way") can carry the 64 MHz 0.005 column, and F40's "0.005 strip" (PROJECT_PLAN.md:4718–4721) holds at 64 MHz too.
- Extending the §2 / CLAUDE.md Palace clause to "all three rungs of F27's h-statement at 10 and 64 MHz" is the weekly's call.

**Negative result** (no `PREDICTED RED:` line; any red stops the slot).
- **(P) red:** report every identity's number, park, BLOCKED. Never widen.
- **(Parse-0) / (Tie) red:** cannot happen on an unchanged log and file. If it does, stop, BLOCKED.
- **(T5) ≫ 1e-6, or (C-h) ≠ ours beyond the print, with (P) and (Tie) green:** a reading, not a failure. DONE as a record with a known-issues entry: the first matched-order disagreement on the strip route. It tied at ≤ 8.0e-11 at 0.0075 and 64 MHz, and at ≤ 6.3e-11 at 0.005 and 10 MHz, so it would point at the solve at 3.8 M unknowns and 64 MHz, not at the mesh.
- **Palace aborts or a log check misses:** record it verbatim, open a known-issues entry, BLOCKED.
- **Status 124 or a minute gate trip:** a partial record, not BLOCKED.

**Outcome (2026-10-15 06:00 slot, implementer, run directly by the slot): DONE — F39 complete.** No edit of any kind. Pre-check: `fem-em-solver` Up, no `fem-em-palace` / `xl` computing (xl Exited 137), `pgrep -c python3` 0, the four `f05x_h0050_o2_f64_e*.json` present, `postpro_…_e1` / `_e2` holding `port-S.csv`, no `_e3` / `_e4`. W1 (e3) `20261015T110032Z_OPS-61-item458b-W1-palace-h0050-f64-e3.log`: Status 0, Total 423.4 s (harness 425 s), Setup 98.2 s, 45 its, peak 33.3 G. W2 (e4) `20261015T110745Z_OPS-61-item458b-W2-palace-h0050-f64-e4.log`: Status 0, Total 403.4 s (405 s), Setup 91.0 s, 45 its, peak 33.3 G. Log checks on both: 593 478 elements, box ±0.12 / ±0.12 / ±0.10 m, Rs 4.375e+01 Ω/sq on 21–24, PEC on 1, ω/2π = 6.400e-02 GHz, one excited port (k), one `GMRES solver converged` line, `ND (p = 2): 3833208`. All four excitations inside the 380–435 s bracket (prediction ±10 % of T₁ held). W3 `20261015T111439Z_OPS-61-item458b-W3-read.log` (4 s): (Parse-0) PASS; (Tie) 3.091e-11, (Tie75) 1.961e-11 PASS; **(P) PASS** — reciprocity 4.296e-11, σ_max 0.999823581372, spreads self 0.0112 % / adjacent 0.0078 % / opposite 0.0151 %, pooled separation 2790.87× (`:64–66`); (T5) per entry ≤ 1.598e-10 on all four columns, per class 9.151e-11 / 8.137e-11 / 1.143e-10 (`:67`); `_class_spread` Palace / ours equal to the printed digit (`:68`); **(C-h) Palace 0.0075 → 0.005 at 64 MHz 0.0569 / 0.0235 / 0.0611 %, equal to ours' to the printed digit** (`:86–87`); (C2-h) 6.218e6 / 2.891e6 / 5.344e6× vs ours h = 0.0075 and ≥ 1.405e7× vs ours ×1 (predicted ≥ 235× / ≥ 1 144×). The `OPS-61` row now reads "the h = 0.005 strip rung at 10 / 64 MHz (F36 / F39, items 437 / 458 / 458b)"; no glyph change.

**Spawn line:** "Use the implementer agent to execute `OPS-61` §9 item 458b per `docs/planning/chunks/OPS-61.md` (read item 458's brief and outcome and the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, container `timeout -k 30`; one Palace excitation per window at `-np 8`, the read at `-n 1`, never two heavy windows at once; any Status 124 stops the excitations; the minute-50 gate binds. No edit of any kind; never run an `export` or `*-configs` mode; never open `docs/private/` or `aed_results/`."

### Review ruling 2026-10-15 03:00 (items 457, 458, 458b)

- **Item 457 accepted as written.** The PEC dispatch exists (`scripts/probes/ops61_step1a_palace_ans4.py:5056`, `:5060`, checked by this review), so there is no probe edit. The two copper-worded labels on the shared read path are disclosed in the journal, not edited.
- **F39 is split, as the drafter's arithmetic requires** (≈ 65 min whole; excitation 4 would start past the minute-50 gate). **Item 458** writes the additive 64 MHz h = 0.005 code — `h0050-configs --f64`, the gap-tolerant parser `_ours_rung_s_scan` (step 6c's 64 MHz block has blank lines and node-ID prefixes, `20261011T070008Z_ANS-4-step6c.log:3824–3827`, checked) and `read --exact-strip-h0050 --f64` — and runs e1–e2. `_ours_rung_s` is not edited; (Parse-0) proves the new parser bitwise on the 10 MHz block. **Item 458b** runs e3–e4 and the full read, no edit; it depends on 458 landing, and is withdrawn if 458's e1 window hit Status 124 or exceeded 590 s (the price is then the record, and F39 returns to the 10-17 weekly).
- Both are weekly-enumerated (rule (4)). No band, tolerance or target is set or moved.
