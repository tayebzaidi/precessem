# ANS-6 — Copper birdcage, phantom-loaded, four lumped ports: the `ANS-4` fixture with a σ = 5.8e7 S/m coil

 **Copper birdcage, phantom-loaded, four lumped ports: the `ANS-4` fixture with a σ = 5.8e7 S/m coil — the first AED case on a realistic conductor** (operator directive 2026-09-04, interactive session; **serial on `TH-15` ✅ and `TH-14` ✅**; AED side: HFSS *Finite Conductivity* boundary on the coil with Solve Inside off, plus a PEC-coil column) | 🟡 **Step 3 private-mode comparison writer ✅ 2026-09-28 (04:30 slot, §9 item 60, §10 F6b landed; glyph stays 🟡 — adjudication open).** `_write_comparison(m, *, private=False)` in `06_…py`; private cells from `aed[f"{column}_{order}"][key]` (class_*, S<ij> Cu, identities, loss_partition Cu, metadata rows Cu + PEC in one cell, AED `_column_difference` in two extra private-only columns labelled as the different, complex quantity — nothing computed across); `_aed_fmt` / `_aed_cell` imported by path from `04_…py`; `_load_aed_results(path=None)` resolves at call time; `main()` copies `ANS-4`'s `FEM_EM_ANS_SKIP_PRIVATE` rule via `_write_private_if_available`. New `tests/unit/test_ans6_private_comparison.py`, synthetic AED JSON only (`ans:6` never run, `aed_results/` never opened): **364 sentinel keys, 180 rendered exactly once in their (section, frequency, column, entry, order) cell, 184 unrendered keys absent**; tracked render **byte-identical** to the committed `COMPARISON.md` (9 292 bytes, 0 sentinels); scratch-repo `check_private_leak.py --audit` **exit 0 `clean:`, 622 reference figures**; knob ⇒ no private file, note names it; **negative control (asserted): the leaking variant exits 1 at `COMPARISON.md`, 330 hits = the 330 rendered sentinel tokens predicted**. 2 passed, 1.70 s, 4 s wall, smoke, one process (`20260928T100228Z_ANS-6-step3.log`). Tracked `COMPARISON.md` / `metrics.json` unchanged. The "private-mode comparison writer (none exists)" follow-up below is retired. **Step 2b resolution knob ✅ 2026-09-23 (04:30 slot, §9 item 35).** `FEM_EM_ANS6_RESOLUTION` (m) on the degree-knob route only (alone ⇒ `RuntimeError`, `20260923T094749Z_ANS-6-step2b-resolution-alone-raises.log:39`); additive `resolution=None` on `_build_ladder()` / `_hole_rung()`; `ans6_reproduction_check.py --knob`; no `src/`. **(a) control** degree 1, h = 0.015, 64 MHz, `-n 4`: 80 181 cells, 111 121 unknowns, tracked `metrics.json` S leaves reproduced to 2.647e-15 (≤ 1e-9; the move floor replaced there by this stronger assert, per the item), step-2 gates green (`20260923T094235Z_ANS-6-step2b-a-h0p015-d1-64MHz.log:940–944`, 36 s). **(b) knob reaches the mesh** h = 0.0075: 226 642 cells, 294 763 unknowns, class moves 1.15e-2 – 2.86e-2 (> 1e-6, asserted), gates green; price point PEC / Cu / P1 sweeps 38.2 / 32.8 / 28.9 s, 147.4 s, summed `ru_maxrss` 6.39 GiB (`20260923T093956Z_…-b-h0p0075-d1-64MHz.log:950–953`, 149 s). Rule (a): `TH-14` 11 passed / 83 s (`20260923T094319Z_…ruleA-TH-14.log:1378`), `TH-15`@64 6 passed / 4 skipped / 32 s (`20260923T094451Z_…ruleA-TH-15-64MHz.log:995`). Censuses: setup figures 56 / 25 / 31 / 0 before (`20260923T080435Z_EX-57-next.log:35`) and after (`20260923T094549Z_…census-after.log:112`); docrefs after dead 0 / guide 0 (`:53`), its 14 stale entries all age-driven, none `ans6`. `xl-pending.md` entry 14 → READY. History: **Runnable half ✅ 2026-09-20 (07:30 slot, §9 item 14).** `examples/ansys_benchmarks/ans6_copper_birdcage_four_port_10_64_128MHz/06_copper_birdcage_four_port_10_64_128MHz.py` + same-stem guide, dispatched through `./run_examples.sh -e ans:6`. Both columns (Cu, PEC) at all three frequencies come from one call to `TH-14`'s own `_build_ladder()` — its `ladder` fixture body lifted, unchanged, to a module-level function (rule (a); `TH-14`'s 11-test gate suite re-run green in the same slot, `20260920T123746Z_ANS-6-TH14-gate-relift.log`, 101.62 s). Every imported gate (reciprocity, passivity, C4 spreads, both columns, all three frequencies) and the copper surface-loss identity (≤ `DISCRETE_IDENTITY_RTOL`) pass. **PEC column cross-check:** `TH-15` step 3's own `_hole_rung`, independently meshed and solved at 10 MHz, matches `TH-14`'s embedded PEC reference to a worst relative deviation of 3.331e-13 (band 1e-8) — two separate builds of "the coil as a hole" agreeing to 13 digits. `P_coil/P_in` (0.9288392 / 0.4483664 / 0.2186740 at 10/64/128 MHz) is not re-derived — read directly off the live `_build_ladder()` call, reproducing the SPEC's reference table and `TH-14`'s own log to every printed digit. Negative control (predicted, rule (e)): the σ = 5.8e11 S/m rung's `max|S−S_PEC|` (2.4–4.2e-06) sits ≈25–40× under `PREDICTED_PEC_LIMIT_MAX_ABS_DS` (1e-4) at all three frequencies. `metrics.json` and `COMPARISON.md` written (AED columns blank by construction, no `aed_results/`/`docs/private/` opened); setup figure and combined XDMF (128 MHz, copper, P1-driven) land alongside. Docrefs census 68→70 guides (the new guide + the new `COMPARISON.md`), dead=0, exit=2 (pre-existing staleness, unrelated); setup-figure census 54→55 examples, 20→21 `ok`, broken=0. Run: `20260920T124620Z_ANS-6.log`, Status 0, 167 s at `-n 2`, degree 1. **Not claimed:** no adjudication, no resonance/tuning/absolute-accuracy claim, no absolute `S₁₁`/`Z_in` at 10 or 64 MHz (`PORT-21`). Operator AED replication and adjudication remain open (§5.4 Waiting-on-you).  **AED numbers landed 2026-09-20** (operator, HFSS 2026 R1, Cu and PEC × Zero and First Order; Linux half re-run `20260920T175038Z_ANS-6-aed-compare.log`, 140 s at `-n 4`, all gates green). **Adjudication owed to the next review** (weekly 2026-09-26, or a daily if asked); preliminary reading (interactive session, `docs/private/ans6-operator-notes-2026-09-20.md`): the **loss partition agrees** at every frequency and both AED orders — the first external check of the Leontovich surface-loss model; the **column difference S(Cu) − S(PEC) agrees** with HFSS First Order; and the absolute S classes carry the *same* self-class offset on the copper and the PEC columns alike, which **attributes `ANS-4`'s self-class disagreement to the lumped-port model / our degree-1 mesh and clears the resolved-conductor interior** (the SPEC's discriminator, answered). AED's own energy balance closes on every design. Follow-ups for the review: ~~a private-mode comparison writer for this example (none exists; the reading above is by hand)~~ (retired by step 3, 2026-09-28), and the self-class offset's home (`PORT-21`, or a degree-2 `xl` rung of this case). **2026-09-21 03:00 review:** adjudication is the 09-26 weekly's (§5.4 — no daily was asked); **step 2 (chain T9b, slot half) queued as §9 item 23** — `FEM_EM_ANS6_DEGREE` knob with the unset run digit-identical at `-n 4` (imported `LEG_D0_REPRODUCTION_BAND`), and the first degree-2 price on the 80 181-cell hole mesh at h = 0.015, 64 MHz, two columns, `-n 8` (predicted ≈ 515 k unknowns, ≈ 45 s, ≈ 12 GiB from `TH-19` step 3); the h = 0.005 `xl` cost probe is written by the review after that price exists. **Step 2 (T9b) knob ✅ + degree-2 price at h = 0.015 (2026-09-21 09:30 slot):** `FEM_EM_ANS6_DEGREE` (+ required `FEM_EM_ANS6_FREQ_MHZ`) on the script; additive `degree` on `_build_ladder()` / `_hole_rung()` and a `sigmas` subset on `_build_ladder()` (default the full ladder), no `src/`. Control: the unset run at `-n 4` reproduces every tracked S leaf to 1.395e-14 against the imported `LEG_D0_REPRODUCTION_BAND` 1e-9 (`scripts/testing/ans6_reproduction_check.py`, `20260921T095749Z_ANS-6-step2-control-check.log:34`); the regenerated `metrics.json` / `COMPARISON.md` were restored, not committed. Rule (a): `TH-14` module green, 11 passed (`20260921T095800Z_ANS-6-step2-gates.log:1476`; the `TH-15` half of that run errored at setup — `TH15_STEP3_FREQ_MHZ` unset, an invocation slip); `TH-15` re-run at 64 MHz green, 6 passed / 4 skipped (`20260921T095937Z_ANS-6-step2-gate-th15.log:995`). **Degree-2 price (64 MHz, `-n 8`, 80 181 cells, Cu + PEC + P1 = three factorisations):** PEC sweep 61.9 s, Cu sweep 53.9 s, P1 44.9 s, 183.9 s total, summed ru_maxrss 13.78 GiB — STOP rule not approached; reciprocity ≤ 4.866e-14, σ_max 0.999740, spreads ≤ 0.0202 %, copper identity 7.812e-13 (all asserted); negative control (asserted > 1e-6, ≤ 2): class moves vs degree 1 self 7.40e-02 / adjacent 4.06e-02 / opposite 4.13e-02 on both columns — self just above the printed-only 3–7e-2 prediction, adjacent/opposite inside it (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:924–942`). The global unknown count was **not printed** (the ≈ 515 k prediction stays unmeasured); the h = 0.005 `xl` entry can price from time and memory. 

### §9 item 193 — 2026-10-03 (chain F, F16 — daily provisional, rule (6))

**Outcome (2026-10-03, 04:30 CDT slot): complete.** `src/fem_em_solver/io/mesh.py`: `as_hole` removed from `leg_sheet_strip_fraction`'s refused set (rule (c), disclosed; error text and docstring updated, no other line). Additive keywords on `tests/mesh/test_birdcage_port_sheets._build` and `test_th14_birdcage_copper._build_ladder` (rule (a)); new `tests/validation/test_ans6_exact_strip.py`. **W0** default-path re-gate: 37 passed (`TH-14` 11, `OPS-61` 1b 8, conductor hole, environment), 262 s, `-n 2` (`20261003T095019Z_ANS-6-item193-W0-regate.log:7612`). **W1** (`-n 1`, 99 s): the strip hole mesh builds, 80 798 cells (predicted [80 600, 81 800], met; ratio 1.007695), (M1)–(M6) green, cavity CAD 40 surfaces, area equal flag on/off within `EXACT`; (N1) flag-off N_i miss +4.202e-2 on every port (predicted ≈ +4e-2) (`20261003T095457Z_ANS-6-item193-W1-mesh.log:1740–1765`). **W2** (`-n 2`, 201 s, 31 passed / 3 (R) skipped): 111 946 unknowns, ladder 103.9 s; identity residual 1.597e-12 / 8.018e-13 / 9.354e-13; P_coil/P_in 0.9287082 / 0.4478617 / 0.2182683 (−0.014 / −0.113 / −0.186 % vs ragged, predicted ±2 %, met); PEC-limit max\|dS\| 2.4–4.3e-6 (≤ 1e-4, met); (c) 10 MHz bracket against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]` green (copper self 2.405e-4 vs exact solid-800 9.976e-2); (N2) 64 MHz Cu move vs the tracked ragged `metrics.json` 1.651 / 0.930 / 0.829 % (self inside the predicted 0.5–2 %, couplings smaller, as predicted); 10 / 128 MHz moves printed 3.165 / 1.055 / 1.051 % and 0.683 / 0.686 / 0.492 % (`20261003T095650Z_ANS-6-item193-W2-solve.log:3089–3120`). **Printed miss:** κ pooled 1.61e-7 at 10 MHz sits below its predicted [5e-6, 5e-5] (64 / 128 MHz 6.58e-6 / 2.64e-5, inside) — a diagnostic, never gated. **W3** closing window with both records registered (M7, (R) asserted): 34 passed, 210 s (`20261003T100029Z_ANS-6-item193-W3-solve.log:3189`). Heavy tier; ≈ 13 min of compute. Rule (j) proposed Done-when sentence for the review: "Step 4 is done when the hole-route exact rectangle passes (M1)–(M7) and the `TH-14` gates at 10 / 64 / 128 MHz with κ off, class means registered."

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.* Licence: rule (6), toward chain F's F10a / F10b target and the operator directive 2026-10-03 (a). The sheets-and-κ ruling (§9, 2026-10-02, items (i)–(iii)) applies. No band, tolerance, Done-when or validation target moves.

`ANS-6` step 4: exact SPEC port rectangles on the hole / copper route, then the copper four-port re-recorded at 10 / 64 / 128 MHz.
- **Executor:** implementer. One `src/` validator line (rule (c), disclosed), two additive helper keywords (rule (a)) and a new module.
- **Tier:** heavy. W0, W2 and W3 at `-n 2`; W1 at `-n 1`. Complex build.
- **Class:** chain step (solve).
- **Time:** ≈ 14 measured / ≈ 40 predicted min.
- **Depends on:** nothing. Item 195 (T17) and item 196 (T18) depend on this one.

**Why: what the sheets are today.**
- **The fixture.** `ANS-6` and `ANS-7` share it. It is `TH-14`'s hole mesh: `tests/mesh/test_birdcage_port_sheets._build(True, as_hole=True)`, 80 181 cells, built inside `tests/validation/test_th14_birdcage_copper._build_ladder` (`:263–303`). Its ports come from `_build_ports` (`:185–226`), which does the same thing `ANS-4` does: one whole F_i face, narrowed by `_narrowed_radial` at `GATED_WIDTH_FRACTION`. That gives a ragged sub-sheet triangulated out of the box. The SPEC sheet is the same 0.007 × 0.008 rectangle (`examples/ansys_benchmarks/ans6_copper_birdcage_four_port_10_64_128MHz/SPEC.md:76`).
- **The blocker.** The exact option exists only on the solid builder. `birdcage_port_domain` refuses `leg_sheet_strip_fraction` together with `as_hole` (`src/fem_em_solver/io/mesh.py:3948–3953`; docstring `:3716–3723`), and the hole builder `_build` has no such keyword.
- **Why the refusal can go.** It is conservative, not structural:
  - the strip block (`:4553–4568`) only fragments the sheet before placement and appends every piece to `sheet_tags`;
  - under `as_hole`, `sheet_tags` enter the cut-air fragment as dim-2 tools exactly like the whole sheet (`:4783–4786`);
  - the half-box labelling (`:4882–4901`) and the cavity read (a face that bounds one meshed volume, `:4925–4953`) do not depend on how many faces a sheet has.

  The refusal's own message gives the one-surface-per-sheet need as its reason. That need belongs to `c4_congruent_sheets` (item 148's brief: it binds one sheet surface per port) and to `symmetry_quadrant`, and both of those stay refused.

**Code.**
1. **`src/fem_em_solver/io/mesh.py`** (rule (c), disclosed):
   - In the validator at `:3948`, remove `as_hole` from the refused set. Keep `c4_congruent_sheets` and `symmetry_quadrant`.
   - Update the error text and the docstring sentence at `:3721–3722`.
   - No other line changes. With `leg_sheet_strip_fraction=None`, no new line executes.
   - Pre-check: run `git log --oneline 341fe583..HEAD -- src/fem_em_solver/io/mesh.py`. If anything since then touches lines 3936–3953 or 4537–4573, stop and report.
2. **`tests/mesh/test_birdcage_port_sheets._build`:** add the keyword `leg_sheet_strip_fraction=None`, forwarded only when it is not `None` (the precedent at `tests/mesh/test_birdcage_leg_offset.py:143–144`). Add a docstring paragraph.
3. **`tests/validation/test_th14_birdcage_copper._build_ladder`:** add the keyword `leg_sheet_strip_fraction=None`, forwarded to `_sheets_build` only when it is not `None`. Add a docstring paragraph.
   - `_build_ports` is unchanged. Its κ default is `None`, and the midpoint filter selects exactly the centre rectangle once the strip edges are CAD edges (item 148, A2).
4. **New `tests/validation/test_ans6_exact_strip.py`**, `@complex_only`, opt-in through `FEM_EM_EXACT_ANS6` ∈ {`mesh`, `solve`}. When it is unset every test skips; any other value raises.
   - **Not `FEM_EM_ANS6_*`.** `ans:7`'s `main()` refuses every such variable (`07_…py:849–851`).
   - **Mesh fixture** (both modes). Build the strip hole mesh (`leg_sheet_strip_fraction=GATED_WIDTH_FRACTION`) and the default hole mesh. Measure them with `_sheet_areas`, `_narrowed_radial`, `_projected_extents` and `_tag_volume`, the way `test_ops61_exact_strip_ports._rung` does (`:79–143`), and also read `diag["conductor_cavity_cad"]`.
   - **Solve fixture** (`solve` only). `lad = th14._build_ladder(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION, measure_kappa=True)`. Every other argument stays at its default: degree 1, the full σ ladder, the PEC column, the P1 field solve, and 10 / 64 / 128 MHz through `TH14_STEP2_FREQS_MHZ` unset. κ is therefore off on every sweep, and `measure_kappa` only reads the terminal-form deficit from the P1 field already being solved.
   - Rank 0 writes `logs/ans6_exact/ans6_exact_p1.json` (gitignored): the S matrices per column per frequency, κ, the power records and the census.
   - **Two record constants**, registered in turn, each with a version-tag comment and the (1*) note as in item 192:
     - `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, from W1;
     - `ANS6_EXACT_STRIP_CLASS_MEANS_P1`, from W2: `{"10" | "64" | "128": {"cu" | "pec": {"self" | "adjacent" | "opposite": complex}}}` at 10 digits.

**Anchors** (asserted; bands imported, never restated).

*Mesh*, at W1 and again inside W2 and W3:
- **(M1)** Every F_i equals dx·dz within `SHEET_AREA_BAND`.
- **(M2)** For P1–P4, each of the following matches within `STRIP_BAND`. This is the SPEC rectangle.
  - N_i area against 5.6e-5;
  - w = A/h against 7.0e-3;
  - radial span against 7.0e-3;
  - h against 8.0e-3;
  - centre radius against `RING_RADIUS`.
- **(M3)** F_i ∖ N_i = 5.6e-5 within `STRIP_BAND`, and N_i ⊂ F_i.
- **(M4)** The box halves partition the box within `EXACT` (`tests/mesh/test_birdcage_conductor_hole.py:117`).
- **(M5)** The cavity CAD area with the flag on equals the area with the flag off within `EXACT`. The strip must not touch the coil wall.
- **(M6)** Flag-off hole cells equal the `n_cells` 80 181 in `ans:7`'s tracked `metrics.json` within `CELL_COUNT_BAND`. This shows the default path is unchanged.
- **(M7)** Flag-on cells equal `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` within `CELL_COUNT_BAND`. Skipped while the record is `None`.

*Solve*, per frequency, as `TH-14` asserts them:
- **(G)** For the Cu and PEC columns: reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and Z-class spreads ≤ `ADJACENT_SPREAD_BAND`.
- **(b)** The copper identity residual ≤ `DISCRETE_IDENTITY_RTOL`.
- **(σ)** The σ-ladder per-class max\|S_σ − S_PEC\| is strictly decreasing.
- **(c)** The 10 MHz bracket. **If item 192 has landed**, assert it against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]`: the σ = 800 S/m solid on the same exact rectangles, and the only like-for-like `S_800`. **Otherwise** print it against the ragged `LEG_D_S_MATRIX_10MHZ`, labelled "port-width confounded, not asserted".
- **(R)** At `-n 2` only, skipped while the record is `None`: every class mean matches `ANS6_EXACT_STRIP_CLASS_MEANS_P1` within `STEP3A_RECORD_RTOL`.

**Negative controls (rule (e)).**
- **(N1) Mesh — predicted, printed.** The flag-off hole N_i misses 5.6e-5 by more than `NEGATIVE_CONTROL_FLOOR`, predicted ≈ +4e-2. The hole route's ragged N_i area has never been logged as such; the solid route read +4.202e-2.
- **(N2) Solve — floor asserted, size predicted** (item 159 (K)'s pattern). The 64 MHz Cu class means on the strip mesh move from the tracked ragged `ans:6` `metrics.json` `rungs["64 MHz"]["cu"]["s_matrix"]` by more than `KNOB_REACHES_THE_SOLVE` in every class.
  - A no-op keyword reproduces those leaves to ≈ 1e-14 (1.395e-14, `20260921T095749Z_ANS-6-step2-control-check.log:34`).
  - The size is predicted at 0.5–2 % on the self class, with the couplings moving less. On the solid sibling the move was 1.539 / 0.858 / 0.713 %.

**Predictions** (printed):
- Strip hole cells in [80 600, 81 800]. The solid route gained 1 121 cells, ratio 1.009657.
- κ pooled in [5e-6, 5e-5] at each frequency. The ragged value was 1.06e-2.
- Copper identity residual ≲ 1e-11.
- P_coil/P_in within ±2 % of the ragged 0.9288392 / 0.4483664 / 0.2186740.
- The PEC-limit control (σ = 5.8e11) still ≤ `PREDICTED_PEC_LIMIT_MAX_ABS_DS`.
- The exact-vs-ragged class moves are printed for both columns at all three frequencies.

**Windows** (bare, whole, foreground, Bash `timeout: 660000`, `-s`, `tests/environment` first).
- **W0** (the rule (a)/(c) re-gate on the default path, after the code edits; ≈ 200–300 s).
  - Reference timings: `TH-14` 83 s (`20260923T094319Z_ANS-6-step2b-ruleA-TH-14.log:1378`), the exact-strip module 57 s (`20261001T220051Z_OPS-61-F10a-W2.log`), and the conductor-hole module (two builds, ≈ 60 s, predicted).
  - It must reproduce `TH-14`'s pass set and `OPS-61` 1b's 8 passed.

  `scripts/testing/run_and_log.sh ANS-6-item193-W0-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_th14_birdcage_copper.py tests/mesh/test_birdcage_conductor_hole.py tests/validation/test_ops61_exact_strip_ports.py -v -s --tb=short'"`
- **W1** (the first strip hole build; mesh only; `-n 1`; ≈ 60–90 s, two hole builds of ≈ 20 s each):

  `scripts/testing/run_and_log.sh ANS-6-item193-W1-mesh "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=mesh timeout -k 30 300 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short'"`

  Then register `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`.
- **W2** (`solve`, `-n 2`, ≈ 160–240 s). `ans:6`'s runnable half ran the same ladder in 167 s at `-n 2` (`20260920T124620Z_ANS-6.log`), plus one extra hole build here. Command: W1's, with `FEM_EM_EXACT_ANS6=solve`, `timeout -k 30 590` and `mpiexec -n 2`, labelled `ANS-6-item193-W2-solve`. Then register `ANS6_EXACT_STRIP_CLASS_MEANS_P1`.
- **W3:** W2 relabelled `ANS-6-item193-W3-solve`. This is the closing window (rule (i)), and (R) and (M7) are asserted.

**Traps.**
- **The `as_hole` + strip combination has never been built.** Build it at `-n 1` first: a gmsh failure deadlocks wider runs (item 148's trap). A meshing failure is a `GEO-21`-class finding: record it verbatim, revert, BLOCKED.
- Keep `pec_column` and `field_solve` at their defaults. `pec_column=False` together with `field_solve=True` raises a rank-0-only `KeyError` (`ANS-7` item 157's trap).
- Never pass a float κ, and never `"in_run"`. `_build_ports` casts `float(…)` (`:192–193`).
- Set no `FEM_EM_ANS6_*` variable. Never run `ans:6` or `ans:7`, and leave their `metrics.json` / `COMPARISON.md` untouched. No `aed_results/`, `COMPARISON_private.md` or `docs/private/`.
- Areas, facet counts and `cell_tags.values` are rank-local: reduce before asserting.
- No edit after W3. A module edit after any window re-runs W3 (rule (i)).
- On Status 124: run `pgrep -c python3`, then the cache sweep; report. Never lengthen the timeout.

**Scope.** Covered: the hole-route exact rectangle (opt-in, default digit-identical) and `ANS-6`'s degree-1 copper / PEC four-port re-recorded at three frequencies. Not covered:
- degree 2. The h = 0.015 price is 183.9 s per frequency at `-n 8` for three factorisations (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:942`), so three frequencies are a next step for the weekly;
- the adjudication and the self-class offset re-read against the private columns (the weekly's);
- Palace's copper column (`OPS-61` step 2);
- any κ correction or band;
- the examples' tracked artefacts.

**Status it can move.**
- The `ANS-6` row gains a "step 4" record sentence. No glyph: 🟡 stays until the weekly adjudicates.
- Items 195 and 196 become runnable.
- The ruling's (i) becomes buildable on the hole route.
- **Rule (j):** the row has no step-4 Done-when, so the review adds the proposed sentence.

**Negative result.**
- **W0 red:** the default path moved. Revert `src/` and both helpers, then stop, BLOCKED.
- **W1 build fails:** park on `attempt/*`, with a known-issues entry naming the gmsh message, BLOCKED. Items 195 and 196 stay blocked.
- **(M2) red with (M1) green:** the strip edges did not reach the hole mesh. Report spans and facets, park.
- **(M5) or (M6) red:** stop, BLOCKED.
- **(G), (b) or (σ) red on the strip mesh:** the ragged gates are green, so open a known-issues entry ("exact-rectangle hole mesh fails `<gate>`") and park. Never widen.
- **(c) red against item 192's record:** report the per-class numbers, then stop.
- No `PREDICTED RED:` line. Any red stops the slot.

### §9 item 216 — 2026-10-04 (chain F19)

*Drafted by a general-purpose agent for the 2026-10-04 03:00 daily review; the review reviews it.* Licence: §10 chain F19 as the 2026-10-03 weekly wrote it ("the exact hole route at degree 2, Cu and PEC, 10 / 64 / 128 MHz, κ off: `TH-14`'s gates and the copper identity at the imported bands, class means registered under (1*)"). The target, the bands and the private row (a) reading are the weekly's and are not restated or moved here. The sheets-and-κ ruling (§9, 2026-10-02) binds: exact CAD rectangles, κ off on every sweep, κ printed only.

`ANS-6` step 5: item 193's exact-rectangle hole route at degree 2, both columns, three frequencies, one window per frequency.
- **Executor:** implementer. Tests only (`tests/validation/test_ans6_exact_strip.py`). No `src/` change and no helper edit: `th14._build_ladder` already takes `degree` and reads its frequency selector.
- **Tier:** heavy. W0 at `-n 2`; W1–W4 at `-n 8` (fallback `-n 12`, below). Complex build.
- **Class:** chain step (solve). **Time:** ≈ 18 measured / ≈ 55 predicted min (≈ 12 code + ≈ 27 of windows + ≈ 5 registration + ≈ 10 fixed).
- **Depends on:** nothing. Item 193's module, records and hole route are on `main`.

**Price, and why one window per frequency.** All figures are measured on this mesh (80 798 cells, exact rectangles, degree 2, `-n 8`). Item 196 W2 measured them (`20261003T124111Z_ANS-7-item196-W2-p2.log:2025, :2030, :2955, :2960`):
- PEC sweep 53.1 / 56.0 s;
- σ = 5.8e7 sweep 49.9 / 41.6 s;
- P1 copper field solve 40.4 / 39.9 s.

Item 196 swept Cu only. `TH-14`'s σ-monotone gate needs the full ladder (5.8e7, 5.8e9, 5.8e11), so per frequency the cost is PEC ≈ 55 + three sweeps × ≈ 46 + P1 ≈ 40 ≈ **233 s**.
- **Fixed cost.** Solve mode builds three hole meshes (`meshes` on and off, plus `_build_ladder`'s own). With environment and JIT, item 193 W2 measured 201 s at `-n 2` against a 103.9 s ladder, so ≈ 100 s fixed. Add ≈ 20 s of degree-2 JIT.
- **Three frequencies in one window** ≈ 3 × 233 + 120 ≈ 820 s at `-n 8`, and ≈ 690 s at `-n 12` (≈ 0.8× on the solves). **Neither fits the ≤ 590 s container timeout** that a 660 s foreground call allows. The 1200 s heavy ceiling is not the binding limit; the slot's tool cap is. So the item runs one window per frequency: ≈ 233 + 120 ≈ **355 s** at `-n 8` (bracket 300–480 s).
- **Memory.** Summed `ru_maxrss` was 13.78 GiB at degree 2, `-n 8`, on the 80 181-cell ragged hole mesh with three factorisations (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:942`). Predicted 14–20 GiB, far inside 128 G. Not `xl`. The `xl` rows `ANS-6-step2c` / `2d` (h = 0.005 / 0.0075, 523 816 cells, 267 GiB) are a different mesh and price nothing here.
- **Unknowns.** Predicted [505 000, 530 000] (6.39–6.40 × cells on `ANS-4`'s two degree-2 meshes). Printed by `lad["unknowns"]`, never asserted.

**Code** (tests only; additive; with `FEM_EM_EXACT_ANS6_DEGREE` and `TH14_STEP2_FREQS_MHZ` unset, `ladder` makes the identical call):
1. **`FEM_EM_EXACT_ANS6_DEGREE`** ∈ {1, 2}, read in the solve fixture. Unset passes no `degree` (`_build_ladder`'s default 1). Any other value raises `ValueError` on every rank before any mesh.
   - **Not `FEM_EM_ANS6_*`.** `ans:7`'s `main()` refuses that prefix (item 193's trap), and the `ans:6` script owns `FEM_EM_ANS6_DEGREE`.
   - Parse it only inside the fixture or a helper, never at import. `tests/validation/test_ans7_exact_strip.py` imports `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` from this module, so an import-time raise would break it.
2. **The frequency selector already exists.** `_build_ladder` loops `_selected_keys()`, which reads `TH14_STEP2_FREQS_MHZ` (comma list ⊂ {10, 64, 128}; `test_th14_birdcage_copper.py:120–126`). Pass nothing new. In this module:
   - the parametrized (G) / (b) / (σ) / (R) tests skip a key absent from `ladder["freqs"]`;
   - (c) skips without `"10"`; (N2) skips without `"64"`;
   - `moves` loops only `lad["freqs"]`, as now.
3. **Degree-2 records.** Add `ANS6_EXACT_STRIP_CLASS_MEANS_P2 = None`, then register it from W1–W3's `cu` / `pec` "class means (10 digits)" prints, keyed like `_P1`.
   - Add `ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH = 8`, with a version-tag comment: "complex 0.11 image, `-n 8`, `<log>:<lines>` per frequency; (1*) licence — re-record on an image change, never re-band".
   - `test_r_class_mean_records` picks the record and width by degree, and skips while the record is `None` or at another width.
4. **(N3), degree 2 only.** Every class mean of each present frequency and column moves from `ANS6_EXACT_STRIP_CLASS_MEANS_P1` by more than `KNOB_REACHES_THE_SOLVE`. Print the moves.
5. **Output.** With both variables unset, the JSON path stays `logs/ans6_exact/ans6_exact_p1.json`. Otherwise it is `ans6_exact_p{d}_f{keys}.json`, so the per-frequency windows and item 193's file are never overwritten. Print the degree in the solve header line.
6. **Everything else as item 193:** `leg_sheet_strip_fraction=GATED_WIDTH_FRACTION`, `measure_kappa=True`, the full σ ladder, `pec_column` / `field_solve` at their defaults, and κ never passed.

**Anchors** (asserted; bands imported, never restated). These are `TH-14`'s gates, per present frequency, at degree 2:
- **(G)** For the Cu and PEC columns: reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and Z-class spreads ≤ `ADJACENT_SPREAD_BAND`.
- **(b) The copper identity.** The residual is ≤ `DISCRETE_IDENTITY_RTOL` with `P_surf` > 0. Degree 2 measured 7.812e-13 at h = 0.015 (ragged) and 6.772e-13 at h = 0.005 (`20260929T070008Z_ANS-6-step2c.log:1034`). Predicted ≲ 1e-11.
- **(σ)** The σ-ladder per-class max\|S_σ − S_PEC\| is strictly decreasing. **It has never been run at degree 2**: step 2 and `xl` 2c swept Cu only. Predicted ≈ 10× per 100× σ, as at degree 1.
- **(c) The 10 MHz bracket** against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]`, as committed. **The order is mismatched:** that S_800 is degree 1 and S_PEC is degree 2, and the print says so.
  - Predicted: copper 1.5e-4 – 5e-4 per class, against solid ≥ 6e-2. The degree-1 values were 2.405e-4 against 9.976e-2, and the PEC column's own order move at 10 MHz is a few e-2.
  - The bracket's margin (≈ 300×) survives the mismatch. It stays asserted as `TH-14`'s gate.
- **(M1)–(M7)** on the mesh, unchanged, every window.
- **(R)** At degree 1 / `-n 2` (W0): the item 193 records, all three frequencies. At degree 2 / `-n 8` (W4): the new `_P2` record on W4's frequency.

**Negative controls** (rule (e)):
- **(N3) The degree knob reaches the solve: floor asserted, size predicted** (item 193 (N2)'s pattern). The same comparison (degree 2 against degree 1, both columns) on this route with ragged sheets read 7.40e-2 / 4.06e-2 / 4.13e-2 at 64 MHz (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:924–942`), and `xl` 2c read 4.48 / 2.00 / 2.95 % at h = 0.005. A no-op knob reads ≈ 1e-14.
  - **Size predicted, printed:** 64 MHz self 5–9 %, couplings 2–6 %; 10 MHz self 1–4 %; 128 MHz self 4–9 %. The basis is `ANS-4`'s exact p2-vs-p1 moves of 1.98 / 6.67 / 5.85 % on the self class (item 192 outcome) and the ragged hole route's 7.40 %.
- **(N2)** stays as committed: floor asserted against the ragged degree-1 `metrics.json` at 64 MHz. At degree 2 it measures order plus sheet; the print is labelled so.
- **The PEC-limit control** (σ = 5.8e11, ≤ 1e-4) stays predicted and printed, as in `TH-14`.

**Predictions** (printed, never asserted):
- Z-class spreads ≤ 0.1 % (ragged degree 2 at h = 0.015: ≤ 0.0202 %).
- σ_max within 1e-3 of 1.
- κ pooled in [1e-7, 1e-4] (degree 1 on this route: 1.61e-7 / 6.58e-6 / 2.64e-5; `ANS-4` exact degree 2 at 128 MHz: 4.55e-5).
- P_coil/P_in 0–4 % above item 193's degree-1 0.9287082 / 0.4478617 / 0.2182683 (`xl` 2c degree 2 at 64 MHz: +2.0 % on the ragged degree-1 figure).
- Cu and PEC class means within 0.05 pp of each other on every class (F16: 0.02 pp at degree 1).

**Windows.** Each harness call is the whole Bash command: bare, foreground, Bash `timeout: 660000`. Every window runs `-s` with `tests/environment` first. W1–W4 use §5.1's durable capture verbatim, trailing `; exit \$rc` included. **Order (rule (i)):** code, then W1–W3 (records `None`, (R) skips), then register `_P2`, then W0 and W4 on the module as committed. Nothing is edited after W0.
- **W1** (degree 2, 10 MHz, `-n 8`, ≈ 355 s, bracket 300–480 s):

  `scripts/testing/run_and_log.sh ANS-6-item216-W1-p2-10 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve FEM_EM_EXACT_ANS6_DEGREE=2 TH14_STEP2_FREQS_MHZ=10 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short > /workspace/logs/ans6-item216-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans6-item216-w1-raw.log; cat /workspace/logs/ans6-item216-w1-raw.log; exit \$rc'"`
- **The price gate**, before W2: read W1's PEC / σ / P1 sweep times and the window's elapsed time.
  - Window ≤ 480 s ⇒ W2 / W3 at `-n 8`.
  - Window in (480, 560] s ⇒ W2 / W3 at `-n 12`.
  - Window over 560 s, or a kill ⇒ stop after W1. The other two frequencies become the review's to re-queue from W1's measured price.
- **W2:** W1's command with `TH14_STEP2_FREQS_MHZ=64`, raw file `ans6-item216-w2-raw.log`, labelled `ANS-6-item216-W2-p2-64`.
- **W3:** with `TH14_STEP2_FREQS_MHZ=128`, raw file `ans6-item216-w3-raw.log`, labelled `ANS-6-item216-W3-p2-128`.
- **Register** `ANS6_EXACT_STRIP_CLASS_MEANS_P2` from W1–W3's prints, then make no further edit.
- **W0** (the unset control on the module as committed; degree 1, all three frequencies, `-n 2`). It is item 193 W3 relabelled (34 passed, 210 s). `test_ans7_exact_strip.py` is added: with no `FEM_EM_ANS7_*` set every test skips, which proves the import still works.

  `scripts/testing/run_and_log.sh ANS-6-item216-W0-p1-unset "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve timeout -k 30 450 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`

  (R) and (M7) are asserted at degree 1, as item 193 W3.
- **W4: the degree-2 closing window** (W2's command relabelled `ANS-6-item216-W4-p2-64-close`, raw file `ans6-item216-w4-raw.log`). (R) at degree 2 is asserted on 64 MHz. The 10 / 128 MHz `_P2` entries are registered from W1 / W3 and first asserted by a later run. The row says so (rule (i)'s disclosure clause).
- **Clock rule.** If the slot is past minute 40 when W3 would start, skip W3. Register 10 / 64 only, run W0 and W4, and journal 128 MHz as owed. The review re-queues it. Never shorten a timeout to make a window fit.

**Traps.**
1. **κ at degree 2 on this route is unmeasured.** Item 196 ran `_build_ladder` at degree 2 with `measure_kappa=False`; `ANS-4`'s in-run κ at degree 2 ran on the solid route (`20261003T094036Z_ANS-4-item192-W3-p2.log:2217`).
   - If W1 raises inside `terminal_form_deficit`, κ is a diagnostic and feeds no gate, so the slot may pass `measure_kappa=(degree == 1)`, disclose it, and re-run W1. Any other exception stops the slot.
2. **Do not drop σ rungs to save time.** (σ) is `TH-14`'s gate, and item 196's Cu-only price is not this item's.
3. **Keep `pec_column` and `field_solve` at their defaults.** `pec_column=False` with `field_solve=True` raises a rank-0-only `KeyError` (item 157). Never pass a float κ or `"in_run"`.
4. **`TH14_STEP2_FREQS_MHZ` takes the keys `10` / `64` / `128`**, not Hz. Anything else raises `RuntimeError` in `_selected_keys()`, which is pure env logic, so every rank raises and nothing hangs.
5. **Never overwrite** `logs/ans6_exact/ans6_exact_p1.json`. Do not run `ans:6` or `ans:7`, and do not touch their `metrics.json` / `COMPARISON.md`. Set no `FEM_EM_ANS6_*` variable.
6. **Register degree-2 records at `-n 8` only.** Degree 2 at `-n 2` is ≈ 3× slower and is a different record width.
7. **Copying 10-digit prints.** Interleaved `PASSED` prefixes from other ranks can land mid-line (item 192 W3:2214). Copy the values only, and spot-check one per column per frequency.
8. **Rank-local quantities.** Areas, facet counts and `cell_tags.values` are rank-local and reduced before asserting (the module does this). There must be no ordering comparison or `sqrt` on complex operands.
9. **Status 124.** Run `pgrep -c python3` in the container, then the rubric's cache sweep as one bare harness call; report. Never lengthen a timeout. A W2 / W3 kill does not void W1.
10. **Private data.** No `aed_results/`, `COMPARISON_private.md` or `docs/private/`. The journal and the row carry identities, records and relative moves only, and no absolute S₁₁ or Z_in at 10 / 64 MHz (known-issues 2026-09-19).

**Scope.** Covered:
- the exact hole route at degree 2 on the ×1 mesh, Cu and PEC, at 10 / 64 / 128 MHz (or as many as the price gate and clock allow);
- `TH-14`'s gates and the copper identity;
- degree-2 class means registered under (1*);
- the degree-1 path proved unchanged.

Not covered:
- the private row (a) re-read (the weekly's);
- any h-rung of the hole route at degree 2 (the `xl` 2c / 2d family is ragged and h-refined);
- Palace (F22, which takes order 2 after this item);
- `ANS-7`;
- the examples' tracked artefacts;
- any κ correction or band.

**What a negative result means.**
- **W0 red:** the unset path moved. Revert the module, stop, BLOCKED (rule (d)).
- **W0 red on `test_ans7_exact_strip.py` collection:** the import contract broke. Revert, stop.
- **(G), (b) or (σ) red at degree 2:** the degree-1 gates are green on this mesh, so the order changed the fixture's identities. Open a known-issues entry ("exact hole route fails `<gate>` at degree 2"), park on `attempt/*`, BLOCKED. Never widen. A (σ) red in particular is the first degree-2 σ-ladder measurement on any mesh. Print the ladder per class; it is a finding about the Leontovich-vs-PEC limit at degree 2, not about the band.
- **(c) red:** report per class, with the order mismatch named; stop. The review rules whether an order-matched S_800 is owed.
- **(N3) red:** the degree knob did not reach the solve. Park.
- **W1 over its price gate:** stop as written. That is a price, not a failure.
- No `PREDICTED RED:` line. Any red stops the slot.

**Status it can move.**
- The `ANS-6` row gains a "step 5" record sentence: the gates per frequency, the identity residuals, the degree-2 class means registered, and the order moves. No glyph change; 🟡 stays until the weekly adjudicates.
- It supplies the weekly's private F16 / row (a) re-read at matched order (F16: "Row (a) is re-read on F19").
- F22's order-2 leg ("order 2 after F19") becomes runnable.
- **Rule (j):** the row carries no step-5 Done-when. Proposed: "Step 5 is done when the exact-rectangle hole route at degree 2 passes `TH-14`'s gates and the copper identity at 10 / 64 / 128 MHz with κ off for both columns, class means registered under (1*)."

**Outcome (2026-10-04 07:30 slot): BLOCKED at the price gate.** W1 (degree 2, 10 MHz, `-n 8`, `20261004T123147Z_ANS-6-item216-W1-p2-10.log`) hit Status 124 at 590 s with no traceback on any rank (both PETSc signal-15 blocks are the timeout; `pgrep -c python3` 0 afterwards). It had built the three meshes (42.9 / 41.3 / 42.1 s) and finished the PEC, 5.8e7 and 5.8e9 sweeps, and died inside the 5.8e11 sweep before the P1 field solve; the ladder prints its sweep timings only after P1, so no per-sweep price was read. The box ran ≈ 1.7× slow: W0, the same command as item 193 W3 (210 s, ladder 103.9 s), took 356 s with ladder 169.9 s (`20261004T124227Z_ANS-6-item216-W0-p1-unset.log`, 34 passed / 14 skipped; (R) reproduces item 193's records to 10 digits at all three frequencies; `test_ans7_exact_strip.py` collects and skips). Cache sweep `20261004T124212Z_ANS-6-item216-cache-sweep.log` removed one 0-byte FFCx stub. Code (the `FEM_EM_EXACT_ANS6_DEGREE` knob, frequency skips, (N3), `_P2 = None`, per-degree (R), per-key JSON path) parked on `attempt/ANS-6-step5-20261004T1249Z`; its added cu-vs-PEC gap print is labelled "pp" but computes a relative % — fix the label before re-running W0. Next: at today's speed even `-n 12` predicts ≈ 500 s per frequency, close to the cap; re-queue on a quieter box, or split each frequency's ladder (PEC + σ sweeps in one window, P1 in another) — never a longer timeout.

### §9 item 226 — 2026-10-04 (chain F19, re-pose of item 216, attempt 2 of 4)

*Drafted by a general-purpose agent for the 2026-10-04 10:30 daily review; the review reviews it.* Licence: §10 chain F19 as the 2026-10-03 weekly wrote it ("the exact hole route at degree 2, Cu and PEC, 10 / 64 / 128 MHz, κ off: `TH-14`'s gates and the copper identity at the imported bands, class means registered under (1*)"). This re-poses item 216, BLOCKED at its price gate (`attempts.md`, 2026-10-04 12:49 UTC). **The target is unchanged: every quantity item 216 gated stays gated; only the window plan changes.** Attempt 2 of the step-5 family (cap four). The sheets-and-κ ruling (§9, 2026-10-02) binds.

`ANS-6` step 5: land item 216's parked module, fix one label, split each frequency into a sweeps window and an identity window, and drop the two extra mesh builds from the split windows.
- **Executor:** implementer. Tests only (`tests/validation/test_ans6_exact_strip.py`). No `src/` change and no edit to `test_th14_birdcage_copper.py`: `_build_ladder` already takes `sigmas`, `field_solve` and `pec_column`.
- **Tier:** heavy. P at `-n 1`; W1a–W3b at `-n 8`; W0 at `-n 2`. Complex build.
- **Class:** chain step (solve). **Time:** ≈ 18 measured / ≈ 60 predicted min (≈ 15 code + ≈ 30 of windows at normal speed + ≈ 5 registration + ≈ 10 fixed). The clock rule sheds 128 MHz before the slot overruns.
- **Depends on:** nothing. Branch `attempt/ANS-6-step5-20261004T1249Z` (`e4909c28`) is based on `c2f352ea`, still the last `main` commit to touch the module: `git diff main...attempt/…` names only this file (+103 / −11), checked 2026-10-04 10:30.

**What item 216's W1 measured** (`20261004T123147Z_ANS-6-item216-W1-p2-10.log`, degree 2, 10 MHz, `-n 8`):
- **Three hole meshes**, mesh wall 37.02 / 36.79 / 36.69 s (`:1064, :1892, :2821`), against 19.8 / 19.4 s on the same mesh at `-n 8` at normal speed (`20261003T124111Z_ANS-7-item196-W2-p2.log:1963, :2889`). **The box ran ≈ 1.87× slow.**
- **Ladder progress:** the PEC sweep (`:2844–2852`, column power 1.000), 5.8e7 (`:2873–2881`) and 5.8e9 (`:2902–2910`) completed; the 5.8e11 sweep had finished drive P1 and was printing drive P2 at the kill (`:2911–2919`). No P1 field solve, no κ.
- **No sweep timings**, because the field-solve branch of `_build_ladder` prints them only after P1 (`test_th14_birdcage_copper.py:410–419`). The `field_solve=False` branch prints the PEC and every σ sweep time as soon as the sweeps finish (`:362–368`); the sweeps windows below use it, so W1a gives the first per-σ degree-2 price on this mesh.
- **Projection of item 216's whole per-frequency window at normal speed:**
  - fixed ≈ 60 s for environment, import and JIT (item 196 W2: 409 s less 327.6 s of two ladders less one 23 s mesh; `ANS-4` item 212 W1: ≈ 67 s);
  - the `meshes` fixture's two extra builds plus the M-tests ≈ 50 s;
  - ladder mesh 23 s, PEC 55 s, three σ sweeps 3 × 46 s, P1 40 s;
  - total ≈ **365 s**. At 1.87× that is ≈ 680 s. The kill at 590 s came ≈ 307 normal-seconds in, which is consistent with ≈ 1.9×.
- **W0** (unset, degree 1, `-n 2`): 356 s against item 193 W3's 210 s, ladder 169.9 against 107.2 s, mesh 32.7 against 19.6 s, so ≈ 1.6–1.7× at `-n 2` (`20261004T124227Z_ANS-6-item216-W0-p1-unset.log:932, :3082, :3292`). The same box read ≈ 1.1× slower at `-n 8` than at `-n 2`.

**The re-pose.**
- **`-n 12`: rejected by measurement.** `ANS-4` item 212 W2 at `-n 12` swept no faster than `-n 8` (183.12 against 186.93 s, a mesh with 3× the unknowns of this one) and doubled the mesh and fixed cost (101.97 against 55.74 s; ≈ 129 against ≈ 67 s).
- **Fewer σ points: rejected.** (σ) is `TH-14`'s gate and needs PEC and all three rungs in one record (`test_sigma_ladder_monotone`).
- **The whole window on a quiet box:** 365 s fits only to ≈ 1.55× (≤ 560 s). Not chosen; the split costs ≈ 10 min of code and survives ≈ 1.9×.
- **Chosen: split each frequency** into a *sweeps* part (PEC + 5.8e7 + 5.8e9 + 5.8e11, `field_solve=False`) and an *identity* part (PEC + 5.8e7 + P1, `sigmas=(COPPER_SIGMA,)`), and run neither with the `meshes` fixture.
  - PEC stays in the identity part because `pec_column=False` with `field_solve=True` raises a rank-0-only `KeyError` at `test_th14_birdcage_copper.py:412` (item 157). It costs ≈ 55 s and buys the cross-window (R) check below.
  - **Sweeps** ≈ 60 + 23 + 55 + 3 × 46 ≈ **276 s** at normal speed (bracket 240–330 s); ≈ 515 s at 1.87×.
  - **Identity** ≈ 60 + 23 + 55 + 46 + 40 + ≈ 5 (κ deficit, powers) ≈ **229 s** (bracket 200–280 s); ≈ 430 s at 1.87×, fits to ≈ 2.3×.
  - Every gated quantity keeps a window: (G), (σ), (c), (N2), (N3) in the sweeps part; (b), (G), (c), (N3) and (R) at degree 2 in the identity part.

**Code** (tests only; additive; with every new variable unset the module makes the branch's calls, which are item 193's, in the same order):
0. **Land the branch module by path:** `git checkout attempt/ANS-6-step5-20261004T1249Z -- tests/validation/test_ans6_exact_strip.py`. First confirm `git log -1 --format=%h main -- tests/validation/test_ans6_exact_strip.py` still reads `c2f352ea`. If `main` has moved, merge by hand: a path checkout would drop main's change. (This is code, not an append-only record, so the path checkout is the right tool.)
1. **Fix the mislabelled print** (branch `:405`). It computes \|S̄_cu − S̄_PEC\| / \|S̄_PEC\| in %, not pp. Relabel it `cu vs pec class-mean relative gap |cu-pec|/|pec| (printed; degree 1 read <= 0.0642 %, item 216 W0)`.
   - Degree 1 read 0.0475 / 0.0125 / 0.0235 % at 10 MHz, 0.0642 / 0.0281 / 0.0486 % at 64 MHz, 0.0397 / 0.0305 / 0.0457 % at 128 MHz (`…124227Z…:3089, :3094, :3099`).
   - Predicted at degree 2: ≤ 0.1 % on every class. Printed only.
2. **`FEM_EM_EXACT_ANS6_PART`** ∈ {unset, `sweeps`, `identity`}. Not `FEM_EM_ANS6_*`: `ans:7`'s `main()` refuses that prefix.
   - Validate it inside the `ladder` fixture: a bad value raises `ValueError` on every rank before any mesh. Never raise at import, because `test_ans7_exact_strip.py` imports this module.
   - `sweeps` ⇒ `_build_ladder(..., field_solve=False)`; `identity` ⇒ `_build_ladder(..., sigmas=(th14.COPPER_SIGMA,))`. Both keep `leg_sheet_strip_fraction=GATED_WIDTH_FRACTION`, `measure_kappa=True`, the degree forward, and `pec_column` / `field_solve` otherwise at their defaults. Unset ⇒ the branch's call.
   - **Drop `meshes` from the split windows.** Change `ladder(meshes)` to `ladder(request)`, calling `request.getfixturevalue("meshes")` first when the part is unset; the M-tests have already built it by then, so the order is unchanged. Give M1–M7 and N1 a `skipif` marker that reads the variable without raising (`bool(os.environ.get(PART_ENV, "").strip())`), reason "asserted in window P and W0". This saves ≈ 50 s at normal speed and ≈ 95 s at 1.9×.
   - **Rank-0 guards** in the fixture's print and JSON blocks. `rec["power"]` and `rec["kappa"]` exist only with the field solve, and `rec["sigma"][th14.SIGMA_LADDER[-1]]` only under `sweeps`. Print what exists; write the JSON with `.get`. An unguarded access is a rank-0-only `KeyError`, so the other ranks wait in the next collective until Status 124, which looks like a price overrun (item 157's class).
   - Tests: (b) skips under `sweeps`, and (σ) skips under `identity`. Everything else runs as on the branch.
   - **(P′) the part reaches the ladder: asserted.** One small test, parametrized by key, skips when the part is unset or the key absent. Under `sweeps`, `set(rec["sigma"]) == set(th14.SIGMA_LADDER)` and `"power" not in rec`. Under `identity`, `set(rec["sigma"]) == {th14.COPPER_SIGMA}` and `"power" in rec`.
   - JSON path: with a part set, `ans6_exact_p{d}_f{keys}_{part}.json`. The unset paths stay the branch's.
3. **Register** `ANS6_EXACT_STRIP_CLASS_MEANS_P2` from W1a–W3a's `cu` / `pec` "class means (10 digits)" prints, keyed like `_P1`. Fill in the branch's version-tag comment with `<log>:<lines>` per frequency. The width stays `ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH = 8`.

**Anchors** (asserted; bands imported, never restated). These are item 216's, per present frequency, at degree 2, κ off:
- **(G)** For the Cu and PEC columns: reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, Z-class spreads ≤ `ADJACENT_SPREAD_BAND`. Sweeps windows, and again in the identity windows on the same S.
- **(σ)** Per-class max\|S_σ − S_PEC\| strictly decreasing over 5.8e7 / 5.8e9 / 5.8e11. Sweeps windows. **Never yet run at degree 2.** Predicted ≈ 10× per 100× σ (degree 1 at 128 MHz: self 3.201e-4 > 3.205e-5 > 3.205e-6, `…124227Z…:3120`).
- **(c)** The 10 MHz bracket against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]`, with the order mismatch printed. Item 216's prediction stands: copper 1.5e-4 to 5e-4 per class, against solid ≥ 6e-2.
- **(b) The copper identity.** Residual ≤ `DISCRETE_IDENTITY_RTOL` with `P_surf` > 0. Identity windows. Predicted ≲ 1e-11 (degree 2 on the ragged route: 7.812e-13 / 6.772e-13).
- **(R) at degree 2, `-n 8`:** the identity windows reproduce the PEC and Cu class means the sweeps windows registered, at `STEP3A_RECORD_RTOL`. It is the same mesh, rank count and sweep code, so predicted ≤ 1e-10. This replaces item 216's W4: it runs at all three frequencies instead of one, and it ties the two halves to one fixture.
- **(M1)–(M7)** in window P and in W0 (both with the part unset).
- **The unset path** (W0): item 193's (R) records at `-n 2`, all three frequencies. Item 216 W0 reproduced them to 10 digits.

**Negative controls** (rule (e)):
- **(N3) The degree knob reaches the solve:** floor asserted, sizes predicted as in item 216 (64 MHz self 5–9 %, couplings 2–6 %; 10 MHz self 1–4 %; 128 MHz self 4–9 %).
- **(P′) The part knob reaches the ladder:** asserted (above). A misspelt part runs the whole ladder plus the extra meshes; on a slow box that dies near 590 s, so read the `[ANS-6 step4] solve:` line and the `pec_column=True … field_solve=False` line first.
- **(R) across windows:** a degree or mesh mismatch between the halves moves class means by ≥ 1 % (N3's size) against a 1e-6 band.
- **The PEC-limit control** (σ = 5.8e11, ≤ 1e-4) stays predicted and printed, in the sweeps windows. In the identity windows `_build_ladder`'s own "negative control" line reads the 5.8e7 rung and says "NOT met". That is cosmetic: disclose it, and never read it.
- **(N2)** stays as committed (64 MHz, floor asserted, labelled order plus sheet).

**Predictions** (printed, never asserted):
- Z-class spreads ≤ 0.1 %; σ_max within 1e-3 of 1.
- κ pooled in [1e-7, 1e-4].
- P_coil/P_in 0–4 % above item 193's degree-1 0.9287082 / 0.4478617 / 0.2182683.
- Unknowns in [505 000, 530 000]; memory 14–20 GiB summed.
- At normal speed, PEC 50–60 s and each σ sweep 40–55 s.

**Windows.** Each harness call is the whole Bash command: bare, foreground, Bash `timeout: 660000`. Every window runs `-s` with `tests/environment` first. W1a–W3b use §5.1's durable capture verbatim, trailing `; exit \$rc` included. **Order (rule (i)):** code steps 0–2 → P → W1a, W2a, W3a → register (step 3) → W1b, W2b, W3b → W0. The registration is the only edit after P, and nothing is edited after it. The sweeps windows ran before the record constant existed; the row says so (rule (i)'s disclosure clause).
- **P: the mesh window and box-speed probe** (part unset, mesh mode, `-n 1`; item 193 W1 relabelled). Normal speed is 99 s with mesh walls 19.03 / 19.81 s (`20261003T095457Z_ANS-6-item193-W1-mesh.log`). Expect 19 passed, as there.

  `scripts/testing/run_and_log.sh ANS-6-item226-P-mesh "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=mesh timeout -k 30 300 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short'"`
- **Gate 1, after P.** s_P = mean mesh wall / 19.4 s.
  - s_P ≤ 1.6 (mean ≤ 31 s; W1a ≈ 276 × 1.1 × 1.6 ≈ 486 s) ⇒ go.
  - Otherwise **stop.** The module may land with P as its record (M-tests asserted on the landed code, no solve claimed). Journal the reading. Item 216's W0 speed (1.67×) would have stopped here.
- **W1a: sweeps, 10 MHz** (`-n 8`; ≈ 276 s, bracket 240–330 s at normal speed):

  `scripts/testing/run_and_log.sh ANS-6-item226-W1a-p2-10-sweeps "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve FEM_EM_EXACT_ANS6_DEGREE=2 FEM_EM_EXACT_ANS6_PART=sweeps TH14_STEP2_FREQS_MHZ=10 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short > /workspace/logs/ans6-item226-w1a-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans6-item226-w1a-raw.log; cat /workspace/logs/ans6-item226-w1a-raw.log; exit \$rc'"`
- **Gate 2, after every window at `-n 8`.** s = PEC sweep time / 54.5 s (item 196 W2: 53.1 / 56.0 s). Sweeps windows print it in the `pec_column=True (PEC … s)` line; identity windows print it in the `PEC sweep … s, P1 copper field solve … s` line.
  - The next sweeps window runs only if s ≤ 1.85 (≲ 510 s) and the last window took ≤ 520 s.
  - The next identity window runs only if s ≤ 2.2 (≲ 505 s).
  - A kill ⇒ stop, as below.
- **W2a, W3a:** W1a's command with `TH14_STEP2_FREQS_MHZ=64` / `128`; raw files `ans6-item226-w2a-raw.log` / `ans6-item226-w3a-raw.log`; labels `ANS-6-item226-W2a-p2-64-sweeps` / `ANS-6-item226-W3a-p2-128-sweeps`.
- **Register** (code step 3).
- **W1b, W2b, W3b: identity** (≈ 229 s each, bracket 200–280 s). W1a's command with `FEM_EM_EXACT_ANS6_PART=identity` and the frequency; raw files `ans6-item226-w1b-raw.log` / `…w2b…` / `…w3b…`; labels `ANS-6-item226-W1b-p2-10-identity` / `…W2b-p2-64-identity` / `…W3b-p2-128-identity`. (R) at degree 2 is asserted here.
- **W0: the unset control** on the module as committed (degree 1, all three frequencies, `-n 2`; item 216 W0 relabelled). Normal speed 210 s; 356 s at this morning's speed. Expect item 216 W0's 34 passed, with three more skips for (P′). `test_ans7_exact_strip.py` collects and skips.

  `scripts/testing/run_and_log.sh ANS-6-item226-W0-p1-unset "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve timeout -k 30 450 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **Clock rule.** Start W3a only if the slot minute + 20 × s + 8 ≤ 60. (W3a, three identity windows and W0 come to ≈ 19.5 min at s = 1; 8 min is bookkeeping.) Otherwise 128 MHz is owed: register 10 / 64, run W1b, W2b and W0, and journal it. Identity windows outrank a third sweeps window. Never shorten a timeout to make a window fit.

**Traps.**
1. **κ at degree 2 on this route is still unmeasured**: item 216 died before P1. If an identity window raises inside `terminal_form_deficit`, κ is a diagnostic. Pass `measure_kappa=(degree == 1)`, disclose it, and re-run that window. Any other exception stops the slot.
2. **The rank-0-only `KeyError`** (code step 2's guards). A Status 124 whose log ends after a complete sweep, with no traceback and no further sweep output, is this trap until shown otherwise. Read the guards before calling it a price.
3. **Never drop σ rungs.** Never pass `pec_column=False` with the field solve (item 157). Never pass a float κ or `"in_run"`.
4. **`TH14_STEP2_FREQS_MHZ` takes `10` / `64` / `128`**, not Hz. Anything else raises in `_selected_keys()` on every rank.
5. **No `-k` deselection.** The part markers do the skipping; a `-k a or b` splits into stray argv inside the quoted container command.
6. **Never overwrite** `logs/ans6_exact/ans6_exact_p1.json`, except through W0's unset path (as item 193 W3 and item 216 W0 did). Do not run `ans:6` or `ans:7`. Set no `FEM_EM_ANS6_*` variable.
7. **Width 8.** Register at `-n 8` only, and run the identity windows at `-n 8`: (R) at degree 2 is a width-8 record.
8. **Copying 10-digit prints.** Interleaved `PASSED` prefixes from other ranks can land mid-line; copy the values only and spot-check one per column per frequency against the log.
9. **Rank-local quantities** are reduced before asserting. No ordering comparison or `sqrt` on complex operands.
10. **Status 124.** Run `pgrep -c python3` in the container, then the rubric's cache sweep as one bare harness call; report. A later kill does not void earlier green windows. Item 216's sweep removed one 0-byte FFCx stub (`20261004T124212Z_ANS-6-item216-cache-sweep.log`); if any window has been killed since, sweep before P.
11. **Private data.** No `aed_results/`, `COMPARISON_private.md` or `docs/private/`. The journal and the row carry identities, records and relative moves only, and no absolute S₁₁ or Z_in at 10 / 64 MHz (known-issues 2026-09-19).

**Scope.** Covered:
- item 216's target, unchanged: the exact hole route at degree 2 on the ×1 mesh, Cu and PEC, at 10 / 64 / 128 MHz (as many as the gates and clock allow);
- `TH-14`'s gates and the copper identity;
- degree-2 class means registered under (1*) and reproduced across windows;
- the degree-1 path proved unchanged;
- the split knob and the label fix.

Not covered:
- the private row (a) re-read (the weekly's);
- any h-rung at degree 2;
- Palace (F22);
- `ANS-7`;
- the examples' tracked artefacts;
- any κ correction or band.

**What a negative result means.**
- **P or W0 red:** the landing or the part knob moved the unset path. Revert, stop, BLOCKED (rule (d)).
- **Gate 1 stop:** a price on a loaded box. The module may land with P green. The item is BLOCKED with the reading; attempt 3 is the review's call (a quiet-hour re-queue, or a weekly-commissioned `xl` night).
- **(G), (b) or (σ) red at degree 2:** open a known-issues entry ("exact hole route fails `<gate>` at degree 2"), park on `attempt/*`, BLOCKED. Never widen. A (σ) red is the first degree-2 σ-ladder measurement on any mesh; print the ladder per class.
- **(R) red across windows:** the halves are not one fixture, or the sweep is not deterministic at `-n 8`. Report the relative gaps and stop. Never re-register from the identity windows.
- **(c) red:** report per class, with the order mismatch named; stop.
- **(N3) or (P′) red:** a knob did not reach the solve; park.
- **A kill that gate 2 allowed:** a price, not a failure; journal it.
- No `PREDICTED RED:` line. Any red stops the slot.

**Status it can move.**
- The `ANS-6` row gains its step-5 record sentence: the gates per frequency, the identity residuals, `_P2` registered and reproduced across windows, and the order moves. No glyph change; 🟡 stays until the weekly adjudicates.
- §10 F19 completes (or completes at 10 / 64 MHz with 128 owed). It supplies the weekly's private row (a) re-read at matched order, and F22's order-2 leg becomes runnable.
- §9 item 216's BLOCKED line is discharged. Once the module lands, the review may delete `attempt/ANS-6-step5-20261004T1249Z`.
- **Rule (j):** item 216's proposed Done-when stands unchanged.

**Outcome (2026-10-04 13:30 slot, second item): DONE at 10 / 64 MHz, 128 MHz sweeps owed.** Module landed by path from `attempt/ANS-6-step5-20261004T1249Z` (`main` still at `c2f352ea`), label fixed, `FEM_EM_EXACT_ANS6_PART` added with rank-0 guards and (P′); `_P2` registered for "10" / "64" from W1a / W2a `:1256–1257` (file mtime 18:53:17 UTC, before W1b). Windows, all Status 0: P `20261004T184052Z_ANS-6-item226-P-mesh.log` 105 s, 20 passed (mesh walls 23.30 / 22.49 s ⇒ s_P 1.18); W1a `…184300Z_…-W1a-p2-10-sweeps.log` 287 s; W2a `…184757Z_…-W2a-p2-64-sweeps.log` 293 s; W1b `…185322Z_…-W1b-p2-10-identity.log` 236 s; W2b `…185732Z_…-W2b-p2-64-identity.log` 227 s; W3b `…190127Z_…-W3b-p2-128-identity.log` 244 s; W0 `…190537Z_…-W0-p1-unset.log` 216 s, 34 passed. PEC sweep 55.9–60.2 s (s 1.03–1.11); σ sweeps 44.8–52.9 s. Readings at degree 2: 570 332 unknowns (above the predicted [505 000, 530 000], printed only); (G) Cu σ_max 0.999997 / 0.999901 / 0.999774, spreads ≤ 0.0064 / 0.0214 / 0.0282 %, reciprocity ≤ 1.0e-13; (σ) 10 MHz self 2.472e-4 > 2.486e-5 > 2.488e-6, 64 MHz self 4.251e-4 > 4.261e-5 > 4.263e-6 (couplings likewise ≈ 10× per rung); (c) 10 MHz copper 2.5e-4 / 6.4e-5 / 1.2e-4 vs solid ≥ 2.5e-2; (b) 2.426e-12 / 4.748e-14 / 2.032e-12; (R) W1b / W2b reproduce W1a / W2a to 10 digits (skipped at 128, no record); κ pooled 1.52e-7 / 6.23e-6 / 2.50e-5; P_coil/P_in +0.73 / +6.29 / +8.83 % above item 193's degree-1 (prediction 0–4 %, printed only); (N3) Cu move from degree 1 1.34 / 0.29 / 0.58 % (10), 4.35 / 1.67 / 3.11 % (64), 3.32 / 2.76 / 4.58 % (128), floors green; Cu-vs-PEC gap ≤ 0.0633 %. Identity windows' "negative control … NOT met" line is the cosmetic 5.8e7 read the brief warned of. **Shed:** W3a (128 MHz sweeps) at 18:52 UTC under the slot's minute-42 window-start cap, the brief's priority putting identity windows first; owed: W3a as written, then register `_P2["128"]`, then W3b re-run for (R) at 128.

### §9 item 232 — 2026-10-04 (chain F19, 128 MHz completion of item 226)

**Outcome (2026-10-04 19:30 slot, second item, implementer in-session):** DONE — §10 F19 complete at 10 / 64 / 128 MHz. Module last touched `53ac5fc6` (checked). **W3a** `20261005T004117Z_ANS-6-item232-W3a-p2-128-sweeps.log` Status 0, 285 s, 16 passed / 22 skipped (as predicted): PEC 56.5 s (`:1248`) ⇒ s = 0.97, s_PEC = 0.94, gate passed at slot minute 16; (σ) self 3.120e-4 > 3.126e-5 > 3.127e-6, adjacent 1.133e-4 > 1.135e-5 > 1.135e-6, opposite 1.092e-4 > 1.095e-5 > 1.095e-6 (`:1323`); PEC-limit 3.127e-6 (printed); class means `:1256–1257` equal item 226 W3b's `:1209–1210` to every printed digit. **Registered** `_P2["128"]` from W3a only (spot-checked Cu self real, PEC opposite imag); module mtime 00:46:20Z < W3b start 00:46:26Z. **W3b** `20261005T004626Z_ANS-6-item232-W3b-p2-128-identity.log` Status 0, 248 s, 17 passed / 21 skipped: (R)[128] passes (class means `:1215–1216` identical to the record to 10 digits); (G) reciprocity 2.258e-15, σ_max 0.999773673, spreads ≤ 0.0282 %; (b) 2.033e-12; (N3) Cu 3.3246 / 2.7557 / 4.5846 %, PEC 3.3257 / 2.7559 / 4.5850 % (reproduced); Cu-vs-PEC gap ≤ 0.0458 %; (P′) `{5.8e7}`, power present. Disclosed: the cosmetic 5.8e7 "NOT met" line (`:1202`, trap 8) and the printed ragged P_coil/P_in move +8.626 % vs "PREDICTED within ±2 %" (`:1214`, printed only, κ-off exact vs degree-1 ragged — same as item 226's 128 MHz identity window). **W0** `20261005T005047Z_ANS-6-item232-W0-p1-unset.log` Status 0, 219 s, 34 passed / 24 skipped (24 = item 226's 23 + item 228's new test).

*Drafted by a general-purpose agent for the 2026-10-04 18:00 daily review; the review reviews it.* Licence: §10 chain F19 as the 2026-10-03 weekly wrote it ("the exact hole route at degree 2, Cu and PEC, 10 / 64 / 128 MHz, κ off: `TH-14`'s gates and the copper identity at the imported bands, class means registered under (1*)"). **This is the remainder of item 226, not a new attempt.** Item 226 landed green at 10 / 64 MHz and shed W3a on the slot clock (`attempts.md`, 2026-10-04 19:12 UTC). The attempt count stays at 2 of 4. The target, the bands, the windows and the module are item 226's. The only new work is the 128 MHz sweeps window, one record entry, and the two closing windows that rule (i) owes after that edit. The sheets-and-κ ruling (§9, 2026-10-02) binds: exact CAD rectangles, κ off, κ printed only.

`ANS-6` step 5 at 128 MHz: W3a as item 226 wrote it, register `ANS6_EXACT_STRIP_CLASS_MEANS_P2["128"]`, re-run W3b so (R) is asserted at 128, then W0.
- **Executor:** implementer. Tests only: one dict entry and its version-tag comment in `tests/validation/test_ans6_exact_strip.py`. There is no `src/` change and no logic edit. `FEM_EM_EXACT_ANS6_PART`, the rank-0 guards and (P′) are already on `main` (`53ac5fc6`).
- **Tier:** heavy. W3a and W3b at `-n 8`, W0 at `-n 2`. Complex build.
- **Class:** chain step (solve). **Time:** ≈ 23 measured (class median) / ≈ 20 predicted min. That is ≈ 12.6 min of windows at normal speed (293 + 244 + 216 s), ≈ 3 min of registration and spot-check, and ≈ 5 min of bookkeeping. It fits as a take-next second item.
- **Depends on:** nothing. Before W3a, confirm that `git log -1 --format=%h main -- tests/validation/test_ans6_exact_strip.py` reads `53ac5fc6`. If `main` has moved the module since, read the diff before editing. The registration is an edit, never a path checkout.

**Price, measured on this mesh and this code path** (80 798 cells, 570 332 unknowns, degree 2, `-n 8`, box at normal speed):
- **W3a ≈ 293 s** (bracket 270–330 s). This is item 226's W2a: 293 s elapsed, ladder 229.8 s, PEC 60.2 s, and σ sweeps of 52.9 / 45.5 / 47.7 s (`20261004T184757Z_ANS-6-item226-W2a-p2-64-sweeps.log:1248–1252, :1925`). W1a, the 10 MHz twin, took 287 s. Frequency does not change the operator size, so 128 MHz is priced as 64 MHz.
- **W3b ≈ 244 s** (bracket 220–280 s). This is item 226's own W3b at 128 MHz: ladder 181.2 s, PEC 57.7 s, P1 46.1 s (`20261004T190127Z_ANS-6-item226-W3b-p2-128-identity.log:1190, :1205, :1877`).
- **W0 ≈ 216 s** (bracket 200–260 s). Item 226's W0: ladder 108.5 s, 111 946 unknowns at degree 1 (`20261004T190537Z_ANS-6-item226-W0-p1-unset.log:3082, :3089, :3310`).
- **Memory:** as item 226, far inside 128 G. Not `xl`.

**The box-speed gate (on W3a against W2a's 293 s).**
- **Before W3a.** Start W3a only if the slot minute is ≤ 32: W3a + registration + W3b + W0 + bookkeeping ≈ 24 min at s = 1. If an earlier `-n 8` window in the same slot read a PEC sweep over 1.9 × 60.2 s (≈ 114 s), do not start W3a. Journal the reading; the item stays ⬜ for the next slot.
- **After W3a**, read two speeds:
  - s = W3a elapsed / 293 s;
  - s_PEC = the PEC time on W3a's `[TH-14 step2] f = 1.280000e+08 Hz: pec_column=True (PEC … s)` line / 60.2 s.
  - Take the larger.
- **Go to registration** only if both hold:
  - s ≤ 1.9. At that speed W3b ≈ 464 s fits its 590 s container timeout and W0 ≈ 410 s fits its 450 s.
  - The slot minute + 3 + 7.7 × s + 8 ≤ 60. That is registration, then W3b + W0 at s, then bookkeeping.
- **Otherwise, do not register.** W3a ran on the module as committed (`53ac5fc6`, unedited), so its (σ), (G) and (N3) at 128 are a complete record. Commit the log and the test-results row, journal "registration + W3b + W0 owed (≈ 12 min)", and leave the module untouched. An edit with no closing window after it is the one outcome this gate exists to prevent.
- **A kill at 590 s** (≈ 2.0×) is a price on a loaded box, not a result. Run `pgrep -c python3` in the container, then the rubric's cache sweep as one bare harness call. Register nothing, mark the item BLOCKED with the reading, and never lengthen a timeout.

**Anchors** (asserted; bands imported, never restated). These are item 226's, at 128 MHz, degree 2, κ off:
- **(G)**, Cu and PEC columns: reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and Z-class spreads ≤ `ADJACENT_SPREAD_BAND`. Asserted in W3a and again in W3b. Item 226's W3b already read them green at 128 (σ_max 0.999774, spreads ≤ 0.0282 %).
- **(σ) The σ-monotone ladder at 128 MHz,** asserted in W3a and owed by item 226 (`test_sigma_ladder_monotone[128]`). Per class, max\|S_σ − S_PEC\| must be strictly decreasing over 5.8e7 / 5.8e9 / 5.8e11.
  - Predicted ≈ 3.1e-4 > 3.1e-5 > 3.1e-6 on the largest class. The basis is item 226's W3b, whose 5.8e7 rung at 128 read max 3.120e-4 (`…190127Z…:1202`), together with ≈ 10× per 100× σ at 10 / 64 MHz at degree 2. Degree 1 at 128 read self 3.201e-4 > 3.205e-5 > 3.205e-6.
  - The PEC-limit control (5.8e11 rung ≤ 1e-4) is printed only; predicted ≈ 3e-6.
- **(b) The copper identity at 128** stays as item 226's W3b asserted it (2.032e-12). The re-run W3b asserts it again; predicted ≲ 1e-11, and the same value to ≈ 3 digits.
- **(R) reproduction at 128.** The re-run W3b's PEC and Cu class means reproduce W3a's registered `_P2["128"]` within `STEP3A_RECORD_RTOL`, at width `ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH = 8`. Predicted ≤ 1e-10: W1b / W2b reproduced W1a / W2a to 10 digits.
  - **Pre-registration cross-check** (printed, read by the slot, not asserted): W3a's two "class means (10 digits)" lines should equal item 226 W3b's (`…190127Z…:1209–1210`; `:1209` carries an interleaved `SKIPPED` prefix) to the printed digit. If they differ beyond the last digit, stop before registering. That is the "(R) red across windows" disposition, below.
- **The unset path (W0):** item 193's degree-1 records at `-n 2`, all three frequencies. `test_ans7_exact_strip.py` collects and skips.

**Negative controls** (rule (e); as item 226's):
- **(N3) The degree knob reaches the solve:** floor `KNOB_REACHES_THE_SOLVE` asserted at 128 in W3a and W3b. Predicted to reproduce item 226 W3b's Cu 3.3246 / 2.7557 / 4.5846 % and PEC 3.3257 / 2.7559 / 4.5850 % (`…190127Z…:1212–1213`). A no-op knob reads ≈ 1e-14.
- **(P′) The part knob reaches the ladder:** asserted.
  - W3a: `set(rec["sigma"]) == set(th14.SIGMA_LADDER)`, with no `power`.
  - W3b: `{COPPER_SIGMA}`, with `power` present.
  - Read the `[ANS-6 step4] solve: degree 2, part sweeps, … unknowns 570332, -n 8` line first. A misspelt part or a lost degree shows there.
- **(R) across windows:** a degree or mesh mismatch between W3a and W3b moves the class means by ≥ 2.7 % ((N3)'s size at 128) against `STEP3A_RECORD_RTOL`.

**Predictions** (printed, never asserted):
- **Pass counts.** `tests/environment` contributes 12 passes in every window.
  - W3a: **16 passed / 22 skipped**, from (G), (σ), (N3) and (P′) at 128. (c) and (N2) skip because their frequency is absent, (b) skips under `sweeps`, and (R)[128] skips "not yet registered".
  - Re-run W3b: **17 passed / 21 skipped**. Item 226's W3b read 16 / 22; (R)[128] now passes.
  - W0: **34 passed**. The skips are item 226's 23 plus any tests the 18:00 review's landing of item 228 added to `test_ans7_exact_strip.py`; the working tree has one more top-level test than `53ac5fc6`. All of them skip with no `FEM_EM_ANS7_*` set, so a changed skip count is not a red. A changed pass count is.
- Unknowns 570 332, exactly as item 226. This is now a measured value, not a prediction.
- PEC sweep 55–62 s and σ sweeps 44–54 s at normal speed.
- Cu-vs-PEC relative gap ≤ 0.05 % on every class (item 226 W3b at 128: ≤ 0.0458 %).

**Windows.** Each harness call is the whole Bash command: bare, foreground, Bash `timeout: 660000`. Every window runs `-s` with `tests/environment` first. W3a and W3b use §5.1's durable capture verbatim, trailing `; exit \$rc` included. **Order (rule (i)):** W3a (module as on `main`) → gate → register (the only edit) → W3b → W0. Nothing is edited after the registration. W3a ran before the 128 entry existed; the row says so, as item 226's did for W1a / W2a.
- **W3a: sweeps, 128 MHz** (`-n 8`; ≈ 293 s, bracket 270–330 s):

  `scripts/testing/run_and_log.sh ANS-6-item232-W3a-p2-128-sweeps "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve FEM_EM_EXACT_ANS6_DEGREE=2 FEM_EM_EXACT_ANS6_PART=sweeps TH14_STEP2_FREQS_MHZ=128 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short > /workspace/logs/ans6-item232-w3a-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans6-item232-w3a-raw.log; cat /workspace/logs/ans6-item232-w3a-raw.log; exit \$rc'"`
- **Gate** (above).
- **Register** `"128"` in `ANS6_EXACT_STRIP_CLASS_MEANS_P2` from W3a's `cu` / `pec` "class means (10 digits)" lines, keyed and formatted like `"64"`.
  - In the version-tag comment, replace "128 MHz owed (W3a shed by the slot clock)" with `20261004T…Z_ANS-6-item232-W3a-p2-128-sweeps.log:<lines>`, the "registered by §9 item 226" clause gaining "and item 232 (128 MHz)".
  - Spot-check one value per column against the log.
  - Record the module's mtime; it must precede W3b's start.
- **W3b: identity, 128 MHz, re-run** (`-n 8`; ≈ 244 s, bracket 220–280 s). (R) at 128 is asserted here.

  `scripts/testing/run_and_log.sh ANS-6-item232-W3b-p2-128-identity "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve FEM_EM_EXACT_ANS6_DEGREE=2 FEM_EM_EXACT_ANS6_PART=identity TH14_STEP2_FREQS_MHZ=128 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short > /workspace/logs/ans6-item232-w3b-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans6-item232-w3b-raw.log; cat /workspace/logs/ans6-item232-w3b-raw.log; exit \$rc'"`
- **W0: the unset control** on the module as committed (degree 1, all three frequencies, `-n 2`; ≈ 216 s, bracket 200–260 s):

  `scripts/testing/run_and_log.sh ANS-6-item232-W0-p1-unset "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve timeout -k 30 450 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **If W3b is green and the clock dies before W0**, which the gate should prevent: the edit is uncommitted with one closing window. Commit nothing from the module. Journal W3b as a green reading on an uncommitted entry, park the edit and both logs on `attempt/ANS-6-step5-128-<UTC>`, and mark the item BLOCKED. The next slot re-runs W0 then lands. Never commit the edit without W0.

**Traps.**
1. **Register from W3a only, never from an identity window** (item 226's rule). The identity windows reproduce the record; they do not source it. Item 226's W3b prints are a cross-check, not a source.
2. **Interleaved rank prefixes** land mid-line: item 226 W3b `:1209` reads `SKIPPED    cu class means …`. Copy the values only, and spot-check.
3. **Width 8.** W3a and W3b run at `-n 8`. At any other width (R) skips, and that is not a pass.
4. **`TH14_STEP2_FREQS_MHZ=128`**, the key, not Hz.
5. **No `-k` deselection.** The part markers skip; `-k a or b` splits into stray argv inside the quoted container command.
6. **Durable capture verbatim** (`; exit \$rc` last). Item 216's lesson stands: a window near 590 s on a slow box is a price; never lengthen.
7. **Status 124 or a kill:** run `pgrep -c python3` in the container first, then the rubric's cache sweep as one bare harness call. A Status 124 that ends after a complete sweep with no traceback is the rank-0-only `KeyError` class until shown otherwise (item 226 trap 2); the guards are on `main`, but read them before calling it a price.
8. **The identity window's "negative control … NOT met" line** is the cosmetic 5.8e7 read (item 226). Disclose it; never read it.
9. **Never overwrite `logs/ans6_exact/ans6_exact_p1.json`** except through W0's unset path. Do not run `ans:6` or `ans:7`, and set no `FEM_EM_ANS6_*` variable. W3b's re-run overwrites the gitignored `ans6_exact_p2_f128_identity.json`; that is expected.
10. **Shared working tree.** The 18:00 review's landing of item 228 touches `test_ans7_exact_strip.py`, which W0 collects. If that landing is not yet committed when the slot starts, the tree is dirty and the slot follows its dirty-tree rule. Never stash or check out that file.
11. **Private data.** No `aed_results/`, `COMPARISON_private.md` or `docs/private/`. The journal and the row carry identities, records and relative moves only, with no absolute S₁₁ or Z_in (known-issues 2026-09-19). No AED number anywhere.

**Scope.** Covered:
- (σ) at 128 MHz at degree 2;
- `_P2["128"]` registered under (1*);
- (R) at 128 across the sweeps and identity windows;
- (G), (b) and (N3) re-asserted at 128;
- the degree-1 path re-proved on the module as committed.

With item 226, this completes F19's target at 10 / 64 / 128 MHz.

Not covered:
- **the private row (a) re-read**: the weekly's, done in the gitignored files, never here;
- any h-rung at degree 2;
- Palace (F22);
- `ANS-7`;
- the examples' tracked artefacts;
- any κ correction or band.

**`ANS-6` stays 🟡.**

**What a negative result means.**
- **(σ) red at 128:** a non-monotone degree-2 ladder at 128 after green 10 / 64. Print the ladder per class and open a known-issues entry ("exact hole route fails (σ) at degree 2, 128 MHz"). Register nothing, park W3a's log on `attempt/ANS-6-step5-128-<UTC>`, BLOCKED. Never widen and never drop a rung. It is a finding about the Leontovich-vs-PEC limit at 128, not about the band.
- **(G) or (N3) or (P′) red in W3a:** item 226's W3b was green on the same fixture, so a red here means the sweeps path differs from the identity path at 128. Report, park, BLOCKED.
- **Cross-check or (R) red:** W3a and W3b are not one fixture, or the sweep is not deterministic at `-n 8`. Report the relative gaps per class and column and stop. Never re-register from W3b, and never move `STEP3A_RECORD_RTOL`. Park, BLOCKED.
- **W0 red:** the registration moved the unset path, which a dict entry should not do. Revert the module, stop, BLOCKED (rule (d)).
- **A gate stop or a kill:** a price, not a failure (above).
- No `PREDICTED RED:` line. Any red stops the slot. **Never move a band to land the item.**

**Status it can move.**
- **§10 F19 completes** at 10 / 64 / 128 MHz. Item 226's "128 owed" is discharged, and the review strikes it from the item 226 line.
- **The `ANS-6` row's "Step 5 record" sentence covers 128 MHz.** Replace "the σ-ladder strictly monotone at 10 / 64 MHz … registered at 10 / 64 … **128 MHz sweeps (σ, registration) owed**" with three-frequency wording: (σ) at 128 with its self-class ladder, `_P2` at 10 / 64 / 128 reproduced across windows, and the item 232 logs cited. Keep the rule (i) disclosure: the sweeps windows ran before their record entries existed.
- **No glyph change.** The row carries no Done-when of its own; rule (j) has only item 216's proposed step-level line ("Step 5 is done when the exact-rectangle hole route at degree 2 passes `TH-14`'s gates and the copper identity at 10 / 64 / 128 MHz with κ off for both columns, class means registered under (1*)"). This item meets that step line. The row's own text keeps 🟡 "because the weekly's adjudication is still open", and nothing here adjudicates. The review may adopt the proposed line as the row's step-5 Done-when (rule (j)), marked met.
- **It supplies** the weekly's private row (a) re-read at matched order at all three frequencies (F16: "Row (a) is re-read on F19").

### §9 item 266 — 2026-10-06 (chain F24)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: PROJECT_PLAN §10 chain F, step F24 as the 2026-10-05 weekly wrote it (grep "F24. \`ANS-6\` step 5b′"): "a resolution knob on the exact hole module (`FEM_EM_EXACT_ANS6_RESOLUTION`, beside `_PART` / `_DEGREE` / `_MODE`), unset digit-identical to F19's records (items 226 / 232) at the imported reproduction band; the knob reaches the mesh (cell count printed at h = 0.0075). No solve beyond the reproduction. Target: identity." The contract is `docs/testing/xl-pending.md` entry 31 (F26). **No band, tolerance, Done-when or validation target is set or moved.** The sheets-and-κ ruling (§9, 2026-10-02) binds: exact CAD rectangles, κ off, κ printed only.

**`ANS-6` step 5b′: one mesh knob on `tests/validation/test_ans6_exact_strip.py`, the unset path reproduced, the knob shown to reach the mesh. No solve at h = 0.0075.**
- **Executor:** implementer. Tests only. No `src/` change. No edit to `test_th14_birdcage_copper.py`: `_build_ladder` already takes `resolution` and forwards it verbatim to `_sheets_build` (`test_th14_birdcage_copper.py:263, :313–314`; `ANS-6` step 2b). `_hole_rung` already takes `resolution` (`test_ans6_exact_strip.py:228–238`; item 214). The knob only has to read one variable and pass it to both.
- **Tier:** heavy for every window (WK ≈ 175 s and W0 ≈ 219 s sit at or over the 180 s standard ceiling; W1 / W2 are one degree-2 identity window each, as items 226 / 232 ran them). Every container timeout is ≤ 590 s, and the ceiling is 12 ranks. WK at `-n 2`, W1 / W2 at `-n 8` (the `_P2` record width), W0 at `-n 2` (the `_P1` record width). Complex build.
- **Class:** chain step (code plus reproduction). **Time:** ≈ 17 measured / **≈ 45 predicted** min: ≈ 10 fixed, ≈ 12 code, windows ≈ 14.5 at normal speed (WK ≈ 3 + W1 ≈ 3.8 + W0 ≈ 3.7 + W2 ≈ 4.1), ≈ 6 bookkeeping.
- **Depends on:** nothing. Before editing, confirm `git log -1 --format=%h main -- tests/validation/test_ans6_exact_strip.py` still names item 232's registration commit. If `main` has moved the module since, read that diff first. `test_ans7_exact_strip.py` imports this module (its `_hole_rung`, records and constants), so nothing new may raise at import.

**The contract with `xl` entry 31.** Entry 31 is PENDING PREREQUISITE on exactly one new name: `FEM_EM_EXACT_ANS6_RESOLUTION`. Spell it that way. The weekly's entry text runs it with `FEM_EM_EXACT_ANS6_DEGREE=2`, `TH14_STEP2_FREQS_MHZ=64`, `-n 16` and **no** `_PART`. That means the `meshes` fixture runs too, and it must work under the knob (code step 3).
- **The implementer that lands this item green marks entry 31 READY in the same commit.** Change only its `**Status:**` line to: "READY (2026-10-06, §9 item 266): knob present, unset digit-identical (W0 / W1 / W2 logs), resolution knob reaches the mesh (WK log, cells)". **Do not write, issue, edit or queue its command**, and do not touch its price, readout or decision rule. The review issues the command.
- If any window is red, entry 31 stays PENDING PREREQUISITE, untouched.

**Code** (tests only; every addition additive; with the knob unset, the module makes the same calls in the same order, with the same arguments):
1. **`RESOLUTION_ENV = "FEM_EM_EXACT_ANS6_RESOLUTION"`** beside `DEGREE_ENV` / `PART_ENV` (`:170–175`). Add `_resolution()` in the `_degree()` pattern (`:195–204`). Unset or blank ⇒ `None`. Otherwise `float(raw)`, which must be finite and > 0; anything else raises `ValueError` naming the variable. Also add `_resolution_set()`, a non-raising read for skip markers, in the `_part_set()` pattern.
   - Parse it **in fixtures only, never at import**. Call it first in **both** `meshes` (beside `_degree()` at `:281`) and `ladder` (beside `_degree()` / `_part()` at `:411–412`), so a bad value raises on every rank before any mesh.
   - Do not name it `FEM_EM_ANS6_*`: `ans:7`'s `main()` refuses that prefix. `FEM_EM_ANS7_EXACT_RESOLUTION` is a different module's knob. Leave it alone.
2. **Forward it.**
   - `meshes`: `_hole_rung(GATED_WIDTH_FRACTION, resolution=res)` and `_hole_rung(None, resolution=res)`. `res is None` forwards nothing, so `_hole_rung` is byte-identical to before (its own docstring).
   - `ladder`: when `res is not None`, set `extra["resolution"] = res`. Unset adds no key, so the `_build_ladder` call is item 232's.
   - Rank 0 echoes one line before the ladder call: `[ANS-6 step5b′] knobs: resolution <res|unset>, degree <d|unset>, part <p|unset>, frequencies <TH14_STEP2_FREQS_MHZ|all>`. Read it first in every window.
   - When the knob is set, call `report_peak_rss(comm, label="ANS-6 step5b′ after ladder")` on every rank after `_build_ladder` (import from `fem_em_solver.utils.instrumentation`, as `test_ans4_exact_strip_records.py:74, :449` does). Entry 31's readout needs the memory figure. Do not add it to the unset path.
3. **Tests under the knob.** "Knobbed" means `_resolution()` is not `None`.
   - **M1–M5** (CAD identities: F/dxdz, the SPEC rectangle at `STRIP_BAND`, the complement, the box-halves partition, the cavity) stay asserted unchanged on the refined mesh. They are CAD identities, not mesh ones (item 212 code step 4's rule).
   - **M6 and M7 become (F′) when knobbed** (item 212's pattern, `test_ans4_exact_strip_records.py:119–120, :640–673`). Import `KNOB_CELL_FACTOR` (1.02) from `tests.validation.test_ans4_exact_strip_records`; the module already imports from it at `:52`. Knobbed M7: assert `cells > KNOB_CELL_FACTOR × ANS6_EXACT_STRIP_HOLE_CELL_RECORD`. Knobbed M6: assert flag-off `cells > KNOB_CELL_FACTOR × tracked`. Print each count against its predicted bracket (below), marked "printed". Unset: both tests exactly as now.
   - **N1** (printed only) runs unchanged.
   - **New `test_f_prime_resolution_reaches_the_ladder(ladder)`** (`needs_solve`). It skips when unset; knobbed, it asserts `ladder["cells"] > KNOB_CELL_FACTOR × ANS6_EXACT_STRIP_HOLE_CELL_RECORD`. It is first asserted in entry 31's own window, and skips in every window of this item.
   - **(R)** (`test_r_class_mean_records`): add a skip when knobbed, reason "(R) records are ×1 (h = 0.015); not asserted under `FEM_EM_EXACT_ANS6_RESOLUTION`". The skip goes before the width check.
   - **(N3)** (`test_n3_…`, and `order_moves` in the fixture): when knobbed, suppress it. Leave `order_moves` empty and skip the test with reason "(N3) would mix h and order". This is item 212's suppression of `p2_vs_p1`.
   - **Printed h-move (entry 31's readout, no gate).** Knobbed, rank 0 prints per frequency and per column the class-mean move against the ×1 record at the same degree (`_P2` at degree 2, `_P1` at degree 1), using `_relative_move`. It also prints the Cu − PEC gap per class on the rung beside the same gap computed from the ×1 record. Label the lines "printed, never asserted".
   - **(c) and (N2)** keep their assertions. When knobbed, their print gains a tag such as `[MESH MISMATCHED: comparand is ×1]`, in the style of (c)'s existing order-mismatch tag (`:579–580`).
   - **(G), (b), (σ), (P′)** are unchanged.
4. **JSON path.** Unset: exactly as now. Knobbed: `ans6_exact_p{d}_r{res:g}_f{keys}[_{part}].json`, so F19's `ans6_exact_p2_f*` files and `ans6_exact_p1.json` are never overwritten.
5. **Docstring.** Document the knob under the module's existing knob notes, with "§9 item 266, `ANS-6` step 5b′, xl entry 31's contract".

**Anchors** (asserted; bands imported, never restated):
- **(R) the unset path reproduces F19's records** at `STEP3A_RECORD_RTOL` (`tests/validation/test_ans4_resolution_ladder.py:198`, imported at `test_ans6_exact_strip.py:58`).
  - W1: `ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]` (item 226), at width `ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH = 8`.
  - W2: `_P2["128"]` (item 232).
  - W0: `ANS6_EXACT_STRIP_CLASS_MEANS_P1` at 10 / 64 / 128, at width 2.
  - Backing, same fixture and width: item 226 W2b reproduced `_P2["64"]` to 10 digits (`20261004T185732Z_ANS-6-item226-W2b-p2-64-identity.log:1339`); item 232 W3b reproduced `_P2["128"]` to 10 digits (`20261005T004626Z_ANS-6-item232-W3b-p2-128-identity.log:1324`); item 232 W0 passed `_P1` (`20261005T005047Z_ANS-6-item232-W0-p1-unset.log:3243`). Predicted ≤ 1e-10.
- **M7 / M6 unset:** 80 798 / 80 181 cells within `CELL_COUNT_BAND` (W0).
- **(G), (b) at degree 2** in W1 / W2, as items 226 / 232 asserted them. **M1–M5 on the h = 0.0075 mesh** (WK).

**Negative controls** (rule (e)):
- **(F′) the knob reaches the mesh: asserted** in WK on both hole meshes (`cells > 1.02 ×` the ×1 count). An ignored or misspelt knob reproduces 80 798 / 80 181 exactly, as it does in every unset window (`20261004T184052Z_ANS-6-item226-P-mesh.log:1740, :1746`; `…item232-W0…:1762, :1768`). Predicted separation ≈ 2.8× against the 1.02 floor.
  - **Predicted, printed:** flag-off **226 642** exactly. It is the same generator call (`_sheets_build(True, as_hole=True, resolution=0.0075)`) that `ANS-6` step 2b measured (`20260923T093956Z_ANS-6-step2b-b-h0p0075-d1-64MHz.log:950`). Flag-on in **[226 642, 232 000]**: the ×1 strip ratio 1.007695 gives ≈ 228 400.
  - A flag-off count that differs from 226 642 is an image or mesher drift. Print it and disclose it; it is not a red.
- **The ladder-level (F′)** is armed but skips here (no solve at 0.0075). It reads the same `_resolution()` as `meshes` through one keyword, so check that by reading the code. Entry 31's window asserts it first.
- **Unset reproduction: asserted** through (R), backed above.
- **The echo line** shows a misspelt knob before anything else.

**Predictions** (printed, never asserted):
- **Pass counts** (`tests/environment` contributes 12 in every window). A changed pass count is a red; a changed skip count is not.
  - WK: 20 passed (item 226 P's count, `…184052Z…:1847`).
  - W1: 18 passed (item 226 W2b, `…185732Z…:1417`).
  - W2: 17 passed (item 232 W3b, `…004626Z…:1417`).
  - W0: 34 passed (`…005047Z…:3243`).
  - Skips: each window's prior count + 1 for the new ladder (F′) test. **Read the pass count, not the skips.**
- **WK mesh walls** 50–65 s each at normal speed. Basis: the ragged 0.0075 hole mesh 46.9 s (step 2b); the exact 0.01 hole rung 53.6–55.4 s on a ≈ 1.4× box (item 214, `ANS-7.md` item 214 outcome); the ×1 exact mesh 23.30 s against ragged 22.49 s.

**Windows.** Each harness call is the whole Bash command: bare, foreground, Bash `timeout: 660000`. Every window runs `-s` with `tests/environment` first. W1 / W2 use §5.1's durable capture verbatim, trailing `; exit \$rc` included. **Order (rule (i)):** code steps 1–5 → WK → W1 → W0 → W2. **No edit to the module after WK starts.** If anything is edited later, every window re-runs.
- **WK: the knob reaches the mesh** (mesh mode, degree unset, `-n 2`; no solve). ≈ 175 s, bracket 150–260 s: item 226 P's 105 s less its 45.8 s of mesh, plus 2 × 50–65 s.

  `scripts/testing/run_and_log.sh ANS-6-item266-WK-r0p0075-mesh "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=mesh FEM_EM_EXACT_ANS6_RESOLUTION=0.0075 timeout -k 30 400 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short'"`
- **W1: unset reproduction at degree 2, 64 MHz** (identity part, `-n 8`). This is item 226 W2b relabelled: 227 s (`…185732Z…:1877`), bracket 210–280 s.

  `scripts/testing/run_and_log.sh ANS-6-item266-W1-p2-64-identity "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && mkdir -p /workspace/logs && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve FEM_EM_EXACT_ANS6_DEGREE=2 FEM_EM_EXACT_ANS6_PART=identity TH14_STEP2_FREQS_MHZ=64 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py -v -s --tb=short > /workspace/logs/ans6-item266-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans6-item266-w1-raw.log; cat /workspace/logs/ans6-item266-w1-raw.log; exit \$rc'"`
- **W0: unset control at degree 1, all three frequencies** (`-n 2`). This is item 232 W0 relabelled: 219 s (`…005047Z…:3312`), bracket 200–260 s. `test_ans7_exact_strip.py` collects and skips. It also proves the import-time rule.

  `scripts/testing/run_and_log.sh ANS-6-item266-W0-p1-unset "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_EXACT_ANS6=solve timeout -k 30 450 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans6_exact_strip.py tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **W2: unset reproduction at degree 2, 128 MHz** (identity part, `-n 8`). This is item 232 W3b relabelled: 248 s (`…004626Z…:1877`), bracket 220–290 s. It is W1's command with `TH14_STEP2_FREQS_MHZ=128`, raw file `ans6-item266-w2-raw.log`, and label `ANS-6-item266-W2-p2-128-identity`.
- **Clock rule.** Start W2 only if the slot minute + 5 + 6 ≤ 60 (the window plus bookkeeping). Otherwise shed it. W0 and W1 already carry the unset path at all three frequencies (degree 1) and at degree 2 on entry 31's frequency. Journal "W2 shed (item 232's `_P2["128"]` not re-reproduced)". The READY mark still stands, because the forwarding code does not depend on frequency. Never shorten a timeout to make a window fit.
- **Box speed.** The Status 0 baselines above were measured at s ≈ 1.0–1.1. At ≈ 2× slow, W1 / W2 still fit 590 s (≈ 455 / 500 s) and W0 its 450 s only to ≈ 2.05×. A kill is a price (below).

**Traps.**
1. **Nothing may raise at import.** `test_ans7_exact_strip.py` imports this module. A raising `_resolution()` at module level reds every `ANS-7` window (item 226's `_part` rule). W0 collects `ans:7` and proves it.
2. **An unset path must forward no key.** Passing `resolution=None` explicitly to `_build_ladder` is equivalent there (its default is `None`). Do not "normalise" it to `0.015`: `_sheets_build`'s `RESOLUTION` default is the record's mesh, and an explicit float may not take the identical gmsh path. Forward nothing when unset.
3. **M6 / M7 under the knob** are the trap entry 31 would otherwise buy. Unhandled, they red at 0.0075 in an `xl` window (the cell-record band against ≈ 2.8× the cells). Code step 3 is that guard. WK is where it is proven.
4. **The rank-0-only `KeyError` class** (item 226 trap 2). The new prints read only what exists (`_P2` / `_P1` keys for the frequencies present). A Status 124 with no traceback after a complete ladder is this trap until shown otherwise.
5. **`TH14_STEP2_FREQS_MHZ` takes `10` / `64` / `128`**, not Hz.
6. **No `-k` deselection.** A `-k a or b` splits into stray argv inside the quoted container command.
7. **Width.** (R) at degree 2 is a width-8 record and skips at any other `-n`. That is not a pass. W1 / W2 run at `-n 8` only.
8. **Durable capture verbatim** in W1 / W2 (`; exit \$rc` last).
9. **Status 124 or a kill:** run `pgrep -c python3` in the container, then the rubric's cache sweep as one bare harness call. A later kill does not void earlier green windows.
10. **Never overwrite** `logs/ans6_exact/ans6_exact_p1.json` except through W0's unset path. Do not run `ans:6` or `ans:7`. Set no `FEM_EM_ANS6_*` variable. W1 / W2 overwrite the gitignored `ans6_exact_p2_f64_identity.json` / `…_f128_identity.json`; that is expected.
11. **Rank-local quantities** are reduced before asserting. `_hole_rung`'s counts are already `size_global` / allreduced. No ordering comparison or `sqrt` on complex operands.
12. **Private data.** No `aed_results/`, `COMPARISON_private.md` or `docs/private/`. The journal and the row carry identities, records and cell counts only, with no absolute S₁₁ or Z_in (known-issues 2026-09-19).

**Scope.** Covered:
- the knob, its validation and its forwarding to both fixtures;
- (F′) asserted on both hole meshes at h = 0.0075, with the CAD identities M1–M5 on that mesh;
- F19's records reproduced on the edited module (degree 1 at 10 / 64 / 128; degree 2 at 64 and, clock permitting, 128);
- entry 31 made READY.

Not covered:
- **any solve at h = 0.0075** (that is entry 31);
- the h-statement and the Cu − PEC h-move (F26 / F27, the 10-14 weekly);
- entry 31's command, price or decision rule;
- the private row (a);
- Palace; `ANS-7`; κ; any band.

**`ANS-6` stays 🟡.**

**What a negative result means.**
- **WK (F′) red:** the knob did not reach `_hole_rung`'s mesh. Park on `attempt/ANS-6-step5b-prime-<UTC>`, BLOCKED. Entry 31 stays PENDING.
- **WK M1–M5 red at 0.0075:** the exact rectangles are not exact on the refined hole mesh. That is a finding about the generator, not the knob. Open a known-issues entry ("exact hole route fails M<n> at h = 0.0075"), park, and leave entry 31 PENDING. It would void F26's premise; the weekly reads it. Never widen `STRIP_BAND` or `SHEET_AREA_BAND`.
- **W0, W1 or W2 red:** the edit moved the unset path. Revert the module, stop, BLOCKED (rule (d)). Never re-register a record and never move `STEP3A_RECORD_RTOL`.
- **A kill:** a price, not a failure. Journal the reading; the windows already green stand. Never lengthen a timeout.
- No `PREDICTED RED:` line. Any red stops the slot.

**Status it can move.**
- **`xl` entry 31: PENDING PREREQUISITE → READY** (Status line only, text above). That makes §10 F26 runnable; the review issues the command and the clerk queues it (10-12).
- **§10 F24 completes.** No glyph change. The `ANS-6` row gains at most one clause, "step 5b′: `FEM_EM_EXACT_ANS6_RESOLUTION`, unset reproduction and (F′) at h = 0.0075 (item 266 logs)", if the review wants it there.
- **Rule (j):** no Done-when is proposed. The step is a prerequisite, and its completion is entry 31's READY mark.

**Outcome (2026-10-06 13:30 CDT slot): ✅ GREEN, all four windows, W2 not shed.** Knob `FEM_EM_EXACT_ANS6_RESOLUTION` landed per code steps 1–5 (tests only; no `src/` change); no edit after WK started.
- **WK** `20261006T183201Z_ANS-6-item266-WK-r0p0075-mesh.log` (`-n 2`, Status 0, 161 s; 20 passed, prediction met): (F′) asserted on both hole meshes — flag-off **226 642** (predicted exactly, no mesher drift), flag-on **226 904** (inside [226 642, 232 000]; strip ratio 1.001156, below the ×1 1.007695); M1–M5 green on the 0.0075 mesh; mesh walls 45.0 / 47.0 s (under the 50–65 s prediction).
- **W1** `20261006T183450Z_ANS-6-item266-W1-p2-64-identity.log` (`-n 8`, Status 0, 242 s; 18 passed, prediction met): (R) `_P2["64"]` green, class means identical to the record at 10 digits; (G), (b) green.
- **W0** `20261006T183901Z_ANS-6-item266-W0-p1-unset.log` (`-n 2`, Status 0, 225 s; 34 passed, prediction met): (R) `_P1` green at 10 / 64 / 128; M7 / M6 80 798 / 80 181; `test_ans7_exact_strip.py` collected (import-time rule held).
- **W2** `20261006T184254Z_ANS-6-item266-W2-p2-128-identity.log` (`-n 8`, Status 0, 246 s; 17 passed, prediction met): (R) `_P2["128"]` green, class means identical at 10 digits.
- The echo line `[ANS-6 step5b′] knobs: …` read correctly first in every solve window. `xl` entry 31 marked READY in this commit (Status line only).

## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 Copper birdcage, phantom-loaded, four lumped ports: the `ANS-4` fixture with a σ = 5.8e7 S/m coil | 🟡 The runnable half, the degree and resolution knobs and the step 3 private-mode comparison writer (2026-09-28) have all landed, and the operator's AED numbers are in. The glyph stays 🟡 because the weekly's adjudication is still open. **Step 4 record (F16, §9 item 193, 2026-10-03):** exact SPEC rectangles now build on the hole route (opt-in `leg_sheet_strip_fraction` with `as_hole`, default path re-gated, 37 passed); 80 798 cells, every N_i the 0.007 × 0.008 rectangle to `STRIP_BAND`; `TH-14` gates, copper identity (≤ 1.6e-12) and σ-ladder monotone at 10 / 64 / 128 MHz, degree 1, κ off; Cu / PEC class means registered; 64 MHz Cu move vs ragged 1.651 / 0.930 / 0.829 % (`20261003T100029Z_ANS-6-item193-W3-solve.log`). A record, no glyph. **Step 5 record (F19, §9 items 226 and 232, 2026-10-04):** item 216's module landed with `FEM_EM_EXACT_ANS6_PART` ∈ {unset, `sweeps`, `identity`}, unset path re-proved (W0 34 passed, twice); at degree 2 (570 332 unknowns, `-n 8`, κ off) `TH-14`'s (G) is green at 10 / 64 / 128 MHz (σ_max ≤ 0.999997, spreads ≤ 0.028 %), the copper identity ≤ 2.4e-12, the σ-ladder strictly monotone at 10 / 64 / 128 MHz (first degree-2 measurement, ≈ 10× per 100× σ; 128 MHz self 3.120e-4 > 3.126e-5 > 3.127e-6), `ANS6_EXACT_STRIP_CLASS_MEANS_P2` registered at 10 / 64 / 128 and reproduced across windows to 10 digits (`20261004T184300Z_ANS-6-item226-W1a-p2-10-sweeps.log` … `…190537Z_…-W0-p1-unset.log`; 128 MHz `20261005T004117Z_ANS-6-item232-W3a-p2-128-sweeps.log`, `…004626Z_…-W3b-p2-128-identity.log`, `…005047Z_…-W0-p1-unset.log`). The sweeps windows ran before their record entries existed (rule (i) disclosure). Step 5b′ (F24, §9 item 266, 2026-10-06): `FEM_EM_EXACT_ANS6_RESOLUTION`, unset reproduction (`_P1` at 10 / 64 / 128, `_P2` at 64 / 128, 10 digits) and (F′) at h = 0.0075 (226 642 / 226 904 cells; `20261006T183201Z_ANS-6-item266-WK-r0p0075-mesh.log` … `…184254Z_…-W2-p2-128-identity.log`). *History: `docs/planning/chunks/ANS-6.md`.* 

## §7 state line replaced 2026-10-13 03:00 review (closure; previous text verbatim)

🟡 The runnable half, the knobs and the private comparison writer landed and the operator's AED columns are in; the weekly's adjudication (AGREE on the copper loss partition, §10) holds the glyph for the exact-rectangle re-read. Steps 4 / 5 / 5b′ (F16 / F19 / F24, 2026-10-03 … 06) rebuild the fixture on exact SPEC rectangles through the hole route — `TH-14`'s gates and the copper identity (≤ 2.4e-12) green at degree 1 and 2 at 10 / 64 / 128 MHz, the σ-ladder strictly monotone, class means registered — and `FEM_EM_EXACT_ANS6_RESOLUTION` reaches the mesh (226 904 cells at h = 0.0075), making its `xl` hole-route rung READY.

### Closure audit 2026-10-13 03:00 (daily) — PASS

`auditor` on the 2026-10-12 weekly's F26 / F19 ruling. Our side's evidence meets §4: (G) at 64 MHz, the copper identity 6.628e-13 and (F′) on the h = 0.0075 rung (`20261012T070009Z_ANS-6-step5b.log:4565–4582`, 25 passed, Status 0, 1532 s); the Cu − PEC gap on the rung 0.0615 / 0.0279 / 0.0479 % against ×1 0.0633 / 0.0287 / 0.0495 % (`:4595`; ×1 `20261004T185732Z_ANS-6-item226-W2b-p2-64-identity.log:1233`), within 3.2 % of itself; largest degree-2 identity residual 2.426e-12 (`20261004T185322Z_ANS-6-item226-W1b-p2-10-identity.log:1203`) — the old row's "≤ 2.4e-12" is corrected to "≤ 2.43e-12". Asserts at `tests/validation/test_ans6_exact_strip.py:630–633, :644–645, :664, :739, :769`. The AED-relative levels rest on the weekly's private adjudication (not audited). Caveats recorded: the ledger row's "`test_g_gates_both_columns[64]` SKIPPED" was a misread of interleaved `-s` output — corrected in the ledger; the copper h-statement rests on one rung at one frequency (F31 open); the chunk header's "serial on `TH-15` ✅" should read "on `TH-15`'s hole-route steps 3–3c" (`TH-15` is 🟡 on its own row); the ✅ covers degree-2 absolute S on this fixture only — no tuning, Q or loss-h claim.

### §9 item 462 — 2026-10-15 (§10 chain F38; lane 1, no solve; depends on §9 item 457)

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.*

**Licence.**
- §10 F38 (PROJECT_PLAN.md:4706–4710), weekly-enumerated, §9 rule (4): "`ans:6` `COMPARISON.md` gains a public 'Palace at matched order' section — the ×1 and h = 0.0075 Cu / PEC columns (F26 / F37) with (T5) per class and both codes' h-moves side by side (Palace figures public; no AED value) — and `ans:4`'s gains F34 / F36 the same way; `docrefs` 0 dead. Standard. Serial on F37 PEC."
- §5.4 Palace rules (i)–(iv) (PROJECT_PLAN.md:1173–1197) — publishable with version / config / mesh; matched order on a shared mesh; reproducible; **never a bridge to AED**.
- The pattern is items 416 / 417 (`docs/planning/chunks/ANS-4.md` § "§9 item 416", § "§9 item 417", and their 2026-10-13 03:00 ruling).
- No band, tolerance, Done-when or target is set or moved. No §7 / §2 edit (F40, item 476, does that).

**Depends on §9 item 457 (F37 PEC) landing DONE.** The PEC h = 0.0075 Palace column does not exist until it does. **If item 457 is not DONE when this item comes up, skip to the next item** (do not run a partial `ans:6` section; `ans:4`'s half alone is not this item). Everything else below is pre-traced; the item-457 cells are read from its W-read log at run time.

**What.** Two public sections, each rendered by the case's tracked writer, never hand-written; render-only (no `ans:4` / `ans:6` run; `metrics.json` untouched):
- `examples/ansys_benchmarks/ans6_copper_birdcage_four_port_10_64_128MHz/COMPARISON.md` — appended after `## Provenance` (the file's tail).
- `examples/ansys_benchmarks/ans4_birdcage_four_port_10_64_128MHz/COMPARISON.md` — appended after item 417's `## h-sequence of the adjudicated pairing …` section (the file's tail).
- **Writer edits.** In each script a module-level `PALACE_SECTION` constant (every figure a literal with a `# <log>:<line>` comment) appended **on the public render only** (`if not private`): `06_…py` after the template write in `_write_comparison` (`:486`, template tail `:643–…`); `04_…py` after the h-table append (`:670–674`). The private render's bytes do not change, so no Palace figure enters `COMPARISON_private.md` from this item (§5.4 (iv) is about tracked files, but keeping the private bytes unchanged also keeps `test_ans6_private_comparison.py`'s sentinel positions untouched).
- **Content: ties and moves only — no Palace absolute class mean.** Our class means are already in each file; Palace's equal them to 1e-10, and printing them adds nothing but 10-digit strings for the leak hook.

**`ans:6` section text (the slot transcribes).** Heading `## Palace at matched order (public; \`OPS-61\` steps 2a–2c)`.
- Scope paragraph: "Palace v0.18.1, order 2 — our N1curl degree — on the shared exact-rectangle hole meshes at 64 MHz: ×1 (80 798 elements, `ND (p = 2)` 570 332) and h = 0.0075 (226 904 elements, 1 538 030). Cu is Palace's `Conductivity` boundary at σ = 5.8e7 S/m (a surface impedance, as our Leontovich condition); PEC is Dirichlet. Mesh export, configs and reads: the tracked `scripts/probes/ops61_step1a_palace_ans4.py` (`export --hole-pec`, `export --hole-h0075`, `cu-configs`, `hole-h0075-configs --cu` / `--pec`; `read --hole-pec`, `read --hole-cu`, `read --hole-h0075 --cu` / `--pec`). Ours: the degree-2 class means of `tests/validation/test_ans6_exact_strip.py` (×1) and `xl` entry 31's rung (`20261012T070009Z_ANS-6-step5b.log`). 64 MHz only — no Palace hole run exists at 10 or 128 MHz."
- **Table A — class-mean ties |m_Palace − m_ours| / |m_ours| (self / adjacent / opposite):**

  | Rung | Cu | PEC |
  |---|---|---|
  | ×1 | 1.612e-10 / 1.586e-10 / 2.033e-10 | 1.340e-10 / 1.120e-10 / 2.197e-10 |
  | h = 0.0075 | 1.274e-10 / 4.521e-11 / 6.621e-11 | ⟨item 457 (T5-pec)⟩ |

- **Table B — each code's own ×1 → 0.0075 move |m(h) − m(×1)| / |m(×1)|, % (self / adjacent / opposite):**

  | Column | Palace | Ours |
  |---|---|---|
  | Cu | 0.1206 / 0.0506 / 0.1091 | 0.1206 / 0.0506 / 0.1091 |
  | PEC | ⟨item 457 (C-h-pec), Palace part⟩ | 0.1193 / 0.0501 / 0.1081 |

- **Table C — the copper move (Cu − PEC class-mean difference), |d_Palace − d_ours| / |d_ours| (self / adjacent / opposite):** ×1 1.031e-07 / 2.263e-07 / 1.860e-07; h = 0.0075 ⟨item 457 (T3m-h)⟩ (or "—" if item 457 printed "not measured").
- Closing line: "Ties are printed by the read windows, never gated (the ≤ 1e-6 precedent, 2026-10-13 ruling (1)). No AED value."

**`ans:4` section text.** Heading `## Palace at matched order (public; \`OPS-61\` F21 / F34 / F36)`.
- Scope paragraph: "Palace v0.18.1, order 2 = our degree 2, on the shared exact-strip meshes (σ = 800 S/m solid conductor, κ off): ×1 (117 206 elements, `ND (p = 2)` 753 062), h = 0.0075 (282 198, 1 822 094), h = 0.005 (593 478, 3 833 208). Reads: `scripts/probes/ops61_step1a_palace_ans4.py` (`read --exact-strip-o2`, `read --exact-strip-h0075 [--f64 | --f128]`, `read --exact-strip-h0050`). Ours: `ANS4_EXACT_STRIP_S_P2_X1` and the rung prints of the three logs the h-sequence above cites."
- **Table A — class-mean ties (self / adjacent / opposite):**

  | f (MHz) | ×1 | h = 0.0075 | h = 0.005 |
  |---|---|---|---|
  | 10 | 2.542e-10 / 9.997e-11 / 6.236e-11 | 1.917e-10 / 1.861e-11 / 7.373e-11 | 6.308e-11 / 5.760e-11 / 1.679e-11 |
  | 64 | 2.464e-10 / 9.668e-11 / 9.883e-11 | 4.757e-11 / 8.022e-11 / 3.630e-11 | ⟨F39, §9 item 458, if DONE; else "not bought"⟩ |
  | 128 | 9.930e-11 / 8.685e-11 / 1.495e-10 | 1.275e-10 / 1.189e-10 / 7.497e-11 | not bought (declined, 2026-10-14 weekly) |

- **Table B — each code's own moves, % (self / adjacent / opposite), Palace / ours:**

  | f (MHz) | ×1 → 0.0075 | 0.0075 → 0.005 |
  |---|---|---|
  | 10 | 0.2527 / 0.0447 / 0.1375 — ours 0.2527 / 0.0447 / 0.1375 | 0.0144 / 0.0088 / 0.0161 — ours 0.0144 / 0.0088 / 0.0161 |
  | 64 | 0.2345 / 0.0950 / 0.2098 — ours 0.2345 / 0.0950 / 0.2098 | ⟨F39 (C-h) if DONE⟩ — ours 0.0569 / 0.0235 / 0.0611 |
  | 128 | 0.2204 / 0.1590 / 0.2863 — ours 0.2204 / 0.1590 / 0.2863 | — |

- Closing line as `ans:6`'s.

**The anchor — every figure traced (drafter-read 2026-10-15; logs under `docs/testing/logs/`, `.log` or `.log.gz`).** All checks are verbatim substring on the cited line.

| Cell(s) | Log | Line |
|---|---|---|
| `ans:6` ×1 Cu ties 1.612e-10 / 1.586e-10 / 2.033e-10 | `20261005T093623Z_OPS-61-item241-W6-read.log` | 56 |
| `ans:6` ×1 copper-move 1.031e-07 / 2.263e-07 / 1.860e-07 | same | 57 |
| `ans:6` ×1 PEC ties 1.340e-10 / 1.120e-10 / 2.197e-10 | `20261005T022223Z_OPS-61-item234-W5-read.log` | 49 |
| `ans:6` h Cu ties 1.274e-10 / 4.521e-11 / 6.621e-11 | `20261015T004434Z_OPS-61-item454-W7-read.log` | 70 |
| `ans:6` Cu moves, Palace and ours 0.1206 / 0.0506 / 0.1091 | same | 73 |
| `ans:6` ours PEC move 0.1193 / 0.0501 / 0.1081 | `20261012T070009Z_ANS-6-step5b.log` | 4594 |
| `ans:6` v0.18.1 / 80798 / 570332 | `20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log` | 45 / 68 / 113 |
| `ans:6` v0.18.1 / 226904 / 1538030 | `20261015T003418Z_OPS-61-item454-W3-palace-hole-h0075-cu-e1.log` | 45 / 68 / 104 |
| `ans:4` ×1 64 ties | `20261004T214342Z_OPS-61-item229-W4-read.log` | 41 |
| `ans:4` ×1 10 ties | same | 48 |
| `ans:4` ×1 128 ties | same | 60 |
| `ans:4` h0075 10 ties; moves | `20261013T124505Z_OPS-61-item418-W6-read.log` | 65; 67 |
| `ans:4` h0075 64 ties; moves | `20261013T171447Z_OPS-61-item434-W6-read.log` | 65; 67 |
| `ans:4` h0075 128 ties; moves | `20261013T173812Z_OPS-61-item435-W6-read.log` | 65; 67 |
| `ans:4` h0050 10 ties; moves | `20261013T221851Z_OPS-61-item437-W6-read.log` | 65; 84 |
| `ans:4` ours 64 Δ₂ 0.0569 / 0.0235 / 0.0611 | `examples/ansys_benchmarks/ans4_birdcage_four_port_10_64_128MHz/COMPARISON.md` | 266–268 |
| `ans:4` v0.18.1 / 117206 / 753062 | `20261004T213356Z_OPS-61-item229-W2-palace-o2-f10.log` | 45 / 68 / 110 |
| `ans:4` v0.18.1 / 282198 / 1822094 | `20261013T123438Z_OPS-61-item418-W2-palace-h0075-e1.log` | 45 / 68 / 101 |
| `ans:4` v0.18.1 / 593478 / 3833208 | `20261013T215659Z_OPS-61-item437-W2-palace-h0050-e1.log` | 45 / 68 / 101 |

**Item-457 cells — read at run time from its read log** (the `read --hole-h0075 --pec` window; name from §9 item 457's outcome in `docs/planning/chunks/OPS-61.md`). The probe prints them under the tag `[OPS-61 2c read pec]` (`scripts/probes/ops61_step1a_palace_ans4.py:4895`, `:4949`, `:4958–4961`, `:5020–5022`):
- `(T5-pec) per class |Palace o2 h - ours p2 h|/|ours p2 h|: self <s>% (<a>)  adjacent <s>% (<b>)  opposite <s>% (<c>)` → Table A PEC h-cell = `<a> / <b> / <c>` (the parenthesised `.3e` values).
- `(C-h-pec) Palace o2 h=0.0075 vs Palace o2 x1 (postpro_h05x_pec_o2) per class: self <p>% (…)  adjacent <q>% (…)  opposite <r>% (…); ours' h-move self 0.1193%  adjacent 0.0501%  opposite 0.1081% …` → Table B PEC Palace cell = `<p> / <q> / <r>`; the "ours'" figures must read 0.1193 / 0.0501 / 0.1081 (they are `H0075_HOLE_OURS["pec"]["h_move"]`, `:4539`) — if not, stop.
- `(T3m-h) |d_pal - d_ours|/|d_ours| self <x>  adjacent <y>  opposite <z>; arg …` → Table C h-cell; if the line reads `(K1) / (K3) / (T3m-h): not measured`, write "—".
- Palace elements / ND from item 457's first Palace window (`… elements … 226904`, `ND (p = 2): 1538030`) — must equal item 454's; if not, stop.
- **F39 (§9 item 458)** cells, if it is DONE: its read log's `(T5) per class` line (ties) and `(C-h)` line (Palace 0.0075 → 0.005 move beside ours'); if not DONE, the cells read "not bought" and the item does not wait for it.

**Assertions (new `tests/unit/test_ans4_ans6_palace_sections.py`, on item 416's `_exec_as` / `_render` shape; tmp copies only).**
- **(X) Traced.** Every table row above (typed `(cell, log, line)` tuples, item 457 / 458 rows added with the lines actually read) holds by substring against the tracked log, and each cell string occurs in the committed section.
- **(R4) / (R6) Byte-for-byte.** Each live writer, public mode, on the committed `metrics.json` (and, for `ans:4`, the committed `h_table_s.json`) reproduces the committed `COMPARISON.md`; the Palace heading occurs exactly once, as the file's tail.
- **(P) Private unchanged.** `ans:6`'s `_write_comparison(m, private=True)` is **not** called (it reads `aed_results/`); instead assert from source that `PALACE_SECTION` is referenced only under the `not private` branch (an AST or string check on the function body).
- **(B) No bridge.** Neither section contains "AED" except the closing "No AED value.", nor `aed_results` / `COMPARISON_private` / `docs/private`.
- **(R0)** The repo files are unchanged by the test.

**Negative controls (rule (e); asserted, by construction).**
- (N-cite) `1.917e-10` cited at `418 W6:64` fails (X).
- (N-digit) a tmp copy of each committed file with one Palace digit changed (`ans:6` "2.197e-10" → "2.198e-10"; `ans:4` "1.917e-10" → "1.918e-10") differs from the render.
- (N-old) the pinned pre-change writers render files with Palace-heading count 0 that differ from the committed files.

**Pre-registered change to item 416's test (the review licenses it or not).** `tests/unit/test_ans4_ans6_verdict_sections.py:143–146` cuts `ans:4`'s committed text at item 417's h-section heading before (Δ4); the new `ans:4` Palace section sits after that heading, so (Δ4) needs nothing. For `ans:6`, (Δ6) (`:162–168`) would count every Palace line as stray. The licensed change: for tag "6", cut the committed text at the Palace heading (asserted present exactly once) exactly as `:143–146` does for tag "4"; nothing else exempted. Without the licence the item stops before any edit.

**Executor, tier, cost.** Implementer. Standard tier; no MPI; complex build (the scripts import dolfinx). Lane 1, class median 3 measured minutes; predicted ≈ 40 min (15 fixed + ≈ 18 two writers and a test + ≈ 5 windows). Depends on §9 item 457 (skip if not DONE); reads item 458 if DONE.

**Windows.** Foreground, bare, repo-relative, Bash `timeout: 660000`, pytest with `-s`, no `-k`, no pipe.
1. **Pins (host, plain `git show`):** the two scripts at their last commits at drafting — `04_…py` at `4989ebdd`, `06_…py` at `c7b1b682`; the slot re-reads `git log -1 --format=%h -- <script>` and types the result as a literal — into `logs/f38_pin/<repo-relative path>` (mirrored depth), loaded with `__file__` set to the live path.
2. **Edits** (two scripts, new test, item 416's test as licensed).
3. **W-regen** (gitignored `logs/f38_regen.py`, Write tool; calls each `_write_comparison(m)` with `private=False` on the committed `metrics.json`):
   `scripts/testing/run_and_log.sh ANS-6-item462-regen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 logs/f38_regen.py'"`
   Then host `git diff --stat`: both `COMPARISON.md` changed, no `metrics.json`, no `h_table_s.json`.
4. **W1** (item 416 / 417's modules re-run, rule (c)):
   `scripts/testing/run_and_log.sh ANS-6-item462-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans4_ans6_palace_sections.py tests/unit/test_ans4_ans6_verdict_sections.py tests/unit/test_ans4_h_table.py tests/unit/test_ans6_private_comparison.py -s -v --tb=short'"`
5. **Docrefs (gate `dead=0`, exit ≠ 1; exit 2 = stale report severity):**
   `scripts/testing/run_and_log.sh ANS-6-item462-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
6. **Commit** both scripts, both tests, both `COMPARISON.md` and the logs. The pre-commit leak hook must print `clean:`; `SKIPPED` is not a pass.

**Traps already paid for.**
- **Leak hook: ≥ 8 significant digits** (`check_private_leak.py:76`, item 417's refusal). Every figure above is ≤ 4 significant digits; never add a Palace class mean or σ_max.
- **Docrefs reads `*.json` / `*.msh` / `*.csv` tokens as artifacts** (`check_example_doc_references.py:116–122`): never name `h05x_cu_o2.json`, a `.msh`, `port-S.csv` or a `postpro_…` path in the sections. (`h_table_s.json` is declared in `LISTED_EXAMPLE_ARTIFACTS`; nothing new is.)
- **The 64 MHz-only hole runs.** `OPS-61`'s §7 row says "the ×1 strip and hole meshes at 10 / 64 / 128 MHz" (PROJECT_PLAN.md:1351); the hole runs exist at 64 MHz only (items 234 / 241, `ω/2π 6.400e-02 GHz`). Write 64 MHz; flag the row to the review, do not edit it.
- **Never run `ans:4` / `ans:6`** (a run re-solves and rewrites `metrics.json`); never hand-edit `COMPARISON.md`; never call `_write_comparison(private=True)`.
- **Interleaved output** on Palace W-read lines is not expected (`-n 1`, no pytest), but parse with `re.search`, never `startswith` (item 417).
- **Pinned controls:** literal sha, mirrored depth (`OPS-43`, item 400).
- **No Palace-vs-AED juxtaposition.** The sections sit outside both verdict sections and quote no AED level.
- Formatting: `≈`, never a tilde. No edit after W1 (rule (i)); never delete a log.

**Scope.** Two public sections on two cases from nine to eleven printed read windows. **Not claimed:** an h-statement beyond what F27 / the 10-12 weekly wrote; copper at 10 / 128 MHz on Palace; Palace absolute S; any AED value; the §7 / §2 sentence (F40, item 476).

**Negative result.** Report and stop; park on `attempt/ANS-6-item462-<UTC>` and mark item 462 BLOCKED in the same commit (rule (d)).
- (X) red: the trace or a log moved; print the line; never retype a figure. The review re-traces.
- An item-457 "ours'" figure ≠ 0.1193 / 0.0501 / 0.1081, or its element count ≠ 226 904: item 457 read the wrong block; stop, the review rules.
- (R4) / (R6) red: writer / regen mismatch; fix in-slot, re-run W-regen and W1.
- `test_ans6_private_comparison.py` red: the public-only branch leaked into the private path; fix in-slot if so, else park.
- Leak-hook refusal: unstage and park; report the location only.
- No `PREDICTED RED:` line.

**Status it can move.** No glyph: `ANS-4` and `ANS-6` are ✅ with Done-whens met (PROJECT_PLAN.md:3527, :3530); rule (j) read, no mismatch. §10 F38 → DONE (the review marks it), which unblocks F40 (item 476). Known-issues: none.

**Outcome (2026-10-15 06:00 slot, take-next, implementer): DONE.** Preconditions read first: item 457's read `20261015T095049Z_OPS-61-item457-W6-read.log` — (T5-pec) 1.378e-10 / 1.153e-10 / 1.484e-10 (`:70`), (C-h-pec) Palace 0.1193 / 0.0501 / 0.1081 with "ours'" 0.1193 / 0.0501 / 0.1081 as required (`:73`), (T3m-h) 3.581e-08 / 2.515e-07 / 1.790e-07 (`:89`); its first Palace window `20261015T094114Z_OPS-61-item457-W2-palace-hole-h0075-pec-e1.log` v0.18.1 / 226 904 / `ND (p = 2): 1538030` (= item 454's). F39 from `20261015T111439Z_OPS-61-item458b-W3-read.log` `:67` / `:87`. Pins `4989ebdd` (`04_…py`) and `c7b1b682` (`06_…py`), re-read by `git log -1`. Each writer gains a `PALACE_SECTION` constant appended under `if not private` only; both `COMPARISON.md` regenerated append-only (`20261015T112311Z_ANS-6-item462-regen.log`; no `metrics.json` / `h_table_s.json` change). New `tests/unit/test_ans4_ans6_palace_sections.py`; item 416's test changed only by the licensed tag-6 cut. First W1 `20261015T112320Z_ANS-6-item462-W1.log` 1 failed / 11 passed — (X) on three `ans:4` ×1 tuples whose log line prints `self X  adjacent Y  opposite Z`, not the slash form; the tuples were changed to the log's own form (no figure retyped). W1b `20261015T112340Z_ANS-6-item462-W1b.log` 12 passed: (X) 73 cells, (R4) / (R6) byte-for-byte, heading once as tail, (P) AST, (B), (R0), N-cite / N-digit / N-old; items 416 / 417's tests and `test_ans6_private_comparison.py` green. Docrefs `20261015T112349Z_ANS-6-item462-docrefs.log` `dead=0` (exit 2, stale report). Disclosed: the `ans:4` heading keeps the brief's "F21 / F34 / F36"; F39 is credited in the scope paragraph and `[--f64]` added to the h0050 read listing. **Flag for the review (not edited):** the `OPS-61` §7 row's "×1 strip and hole meshes at 10 / 64 / 128 MHz" — the hole runs exist at 64 MHz only (items 234 / 241 / 454 / 457). §10 F38 → DONE is the review's mark; it unblocks item 476.

**Spawn line:** "Use the implementer agent to execute `ANS-6` §9 item 462 per `docs/planning/chunks/ANS-6.md` (read the latest Review ruling paragraph first). First confirm §9 item 457 is DONE; if not, skip this item. The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Never run `ans:4` / `ans:6`, never hand-edit a `COMPARISON.md`, never call `_write_comparison(private=True)`; change item 416's test only by the licensed tag-6 cut; never write a figure of 8+ significant digits, a `.json` / `.msh` / `.csv` name or any AED value into a section; never open `aed_results/`, `COMPARISON_private.md` or `docs/private/`."

### Review ruling 2026-10-15 03:00 (item 462)

Accepted. **Depends on item 457 landing; if it has not, skip to the next item.** If item 458b has landed, the `ans:4` section also carries the 64 MHz h = 0.005 strip column; otherwise it says that rung is not yet bought. **Licensed:** the (Δ6) tag-6 cut at the new Palace heading in `tests/unit/test_ans4_ans6_verdict_sections.py`, disclosed, the item-416 pattern. No `.json` / `.msh` / `.csv` name in a `COMPARISON.md` (docrefs reads them as run artifacts). No Palace number is set against an AED one.
