# ANS-7 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record. Step 3 (i) and (ii) (§9 items 127, 132) are archived verbatim in `docs/planning/plan-archive.md` § "§9 digest and restock of the 2026-10-01 03:00 review".

### §9 item 157 — 2026-10-01 (chain T11b; queued by the 18:00 review)

`ANS-7` step 3 (iii): T11's degree-2 leg, run with the in-run degree-2 κ.
- **Executor:** implementer. Tests only, in this item's own module.
- **Tier:** heavy; `-n 8` for degree 2, `-n 4` for the degree-1 controls; complex.
- **Time:** ≈ 15 measured / ≈ 60 predicted min.
- **Depends on:** nothing. F13 landed as item 140.

**Why now:** the 10:30 review's `log-pathologist` ruling (§7 `ANS-7` row) found the degree-2 R_in / Q fall of ≈ −98 % is not established as an order move. With κ held at its degree-1 value, the raw degree-2 S loses ≈ 8 % of the incident power. `xl` entry 24 is HELD until this item supplies a degree-matched κ (`docs/testing/xl-pending.md:1241–1250`).

- **How the in-run κ reaches this path.**
  - `th14._build_ladder` cannot take F13's string `"in_run"`. Its κ goes through `_build_ports`, which casts `float(width_correction_kappa)` (`tests/validation/test_th14_birdcage_copper.py:185–193`), so the string raises `ValueError`. F13 added it only to the two four-port constructors and excluded `_build_ladder` by name (`docs/planning/chunks/PORT-21.md:31`).
  - **No TH-14 keyword is needed.** `_build_ladder`'s existing `measure_kappa=True` (`:264,381–385`) runs the same route: `PORT-14` step 3's pooled `terminal_form_deficit` on an uncorrected Cu P1 solve. `ans:7` used it for its degree-1 κ (`07_tuned_copper_birdcage_driven_sweep.py:16–17,284–285,314`), and F13's `"in_run"` runs it too (`…leg_offset_sweep.py:401–425`). It solves at the ladder's degree (`:364–366`).
  - **Rejected alternative,** for API parity with F13: accept `width_correction_kappa="in_run"` on `_build_ladder`. That would need a rule (a)/(c) disclosure and a TH-14 re-gate, and it costs a window while saving nothing.
- **Code.** Tests only, in `tests/validation/test_ans7_sparse_route_b.py`, which no other file imports. No `src/` change and no TH-14 edit.
  - **New knob `FEM_EM_ANS7_SPARSE_KAPPA`.**
    - **Unset:** today's `_kappa()` (`metrics.json` `"kappa"`, 1.0599570850048146e-02). Behaviour, JSON names and print text stay byte-identical.
    - **`in_run`:** before the grid, call `th14._build_ladder(deg, sigmas=(th14.COPPER_SIGMA,), resolution=res, frequencies_hz=[KAPPA_FREQUENCY_HZ], width_correction_kappa=None, measure_kappa=True)`. **Leave `pec_column` and `field_solve` at their defaults (True) whatever the module knobs say.** κ is `rec["kappa"]["pooled"]`. Rank 0 writes `logs/ans7_sparse/kappa_p<deg>_h<tag>.json` containing κ (repr), the per-sheet κ, degree, resolution, ranks, cells, the uncorrected Cu S and `rec["power"]`; then a barrier.
    - **`stored`:** read that JSON. Degree and resolution must match. If the file is missing, `pytest.fail`; never fall back to the held κ.
    - **Any other value** raises `ValueError` (F13's pattern).
  - **Imported constants.**
    - `KAPPA_FREQUENCY_HZ = _ANS7.ANCHOR1_FREQUENCY_HZ` = 64.0e6 (`07_…py:104`)
    - `KAPPA_RECOMPUTE_ATOL` = 1e-9 (`07_…py:487`)
    - `REDUCTION_BAND` = 1e-3 (`tests/validation/test_port_lumped_rlc_termination.py:88`)
    - `REFERENCE_IMPEDANCE_OHM` = 50 (`test_port_package_sparameters.py:140`)
    - `PASSIVITY_SIGMA_TOLERANCE` = 1e-9 (`test_port_birdcage_four_port.py:129`)
  - **Other additions to the module.**
    - Per-node JSON gets a `_kin` suffix whenever the knob is set.
    - The fixture keeps `rec["power"]` per node when the field solve runs.
    - Print Im Y_in at every node, both columns.
    - New `HELD_KAPPA_P2_ROWS` (item 132's rows, `20261001T131348Z_ANS-7-step3d2-W3.log:1186–1192`):
      - Cu: f₀ 62657821.491631046 Hz, R_in 355.1839277557437 Ω, Q 18.136752918228133.
      - PEC: 62672359.551981606 Hz, 357.8640796561647 Ω, 18.286075024732494.
    - New `test_kappa_in_run_discriminator`, which runs only with the knob at `in_run` and the field solve on.
- **Anchor: the loss identity on the P1 column.** L_S = 1 − Σᵢ\|S_i1\|² from the Cu sweep's S. L_F = (P_phantom + P_(Ω∖phantom) + P_surf)/P_inc from the same specs' P1 field solve. P_inc = \|V_src\|²/(8 Z₀) = 2.5e-3 W, because V_src = 1 V (`test_th14_birdcage_copper.py:224`) and a₁ = (V + Z₀I)/(2√Z₀) (`src/fem_em_solver/ports/sparameters.py:72–80`). By hand, V₁ + 50 I₁ = 1.0000 V at `…125435Z_…-W0.log:1137`. Both sides come from the window's own solve.
- **Anchors (asserted).**
  - **(K1)** At degree 1: \|κ_in_run − 1.0599570850048146e-02\| ≤ `KAPPA_RECOMPUTE_ATOL`. Backed by 4.441e-16 (`20260929T100040Z_ANS-7-item68-W3-A.log:865`).
  - **(R2)** At degree 2: the uncorrected κ pass reproduces `ANS-6` step 2's degree-2 64 MHz line within `NODE_S_RTOL` 1e-9, i.e. P_phantom 1.574450007e-06 W and P_surf 1.419459787e-06 W (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:937`). The 10-digit print limits agreement to ≈ 3.5e-10.
  - **(L)** For the corrected node: \|L_S − L_F\| ≤ `REDUCTION_BAND`, at degree 1 (W0b) and at degree 2 (W0). Degree-1 backing is 7.1e-7 (hand arithmetic on `20261001T095511Z_ANS-7-step3-W1.log:1249,1256`).
  - **Per node:**
    - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`
    - reciprocity ≤ `RECIPROCITY_BAND` 1e-3
    - identity residual ≤ `DISCRETE_IDENTITY_RTOL` 1e-6
  - **Grid:** (ii)'s assertions, unchanged, at the in-run κ:
    - one − → + zero per column
    - k = ±1 share > `K_SHARE_FLOOR` 0.5
    - degree-2 Δ_LOO ≤ (i)'s `BANDS` 1e-6
    - all seven JSONs carry the same κ to 1e-15
- **Negative control (rule (e)).** **(N)** The uncorrected pass's gap \|L_S − L_F\| > 10 × `REDUCTION_BAND`.
  - **Degree 2, asserted:** backed at ≈ 9.8e-2 (98× the band). Max column power is 9.005e-01 (`20260921T100022Z_…:921`) with L_F 1.1976e-3 (`:937`).
  - **Degree 1, asserted:** backed at 2.169e-2 (21.7×), from `20261001T125435Z_ANS-7-step3d2-W0.log:1230,1251` at `-n 4`.
  - **Ceiling:** the gap is ≤ 1 − L_F ≈ 0.9988, so at most ≈ 999×.
- **Predictions, written first (printed, never asserted).**
  - **κ_deg2 in [4.5e-2, 6.0e-2], point estimate 5.2e-2.** Two independent routes give it:
    1. The sibling fixture's degree-2 / degree-1 ratio: 5.210595564e-02 against 1.060762155e-02 (`20261001T094205Z_PORT-21-step4-W1.log:2943,3049`).
    2. A linear extrapolation of this fixture's degree-2 Cu column power, 0.9005 at κ = 0 and 0.9205 at κ = 1.06e-2, to closure: κ ≈ 5.19e-2.
  - **Column power at 64.5 MHz with κ_deg2:** Cu in [0.995, 0.9995]; PEC in [0.996, 1.000].
  - **Gap with κ_deg2:** ≤ 1e-4.
  - **Grid rows:**
    - Q_Cu in [700, 1250]; Q_PEC in [1100, 2200].
    - R_in,Cu in [13 000, 24 000] Ω; R_in,PEC in [21 000, 43 000] Ω.
    - The copper part of 1/Q survives the order step. The phantom's field loss falls ≈ 15 % from degree 1 to degree 2 at κ = 0.
  - **f₀ in [61.1, 64.2] MHz;** the sign is not pre-registered.
  - **k-share** ≥ 0.999.
  - **Δ_LOO** ≤ 1e-9.
- **Reading rule (pre-registered, for the review and T13).** Let m = (p2,in-run − p1)/p1 for Q and R_in, in both columns.
  - **"κ confound":** (L) green at W0 and every \|m\| ≤ 30 %.
  - **"Genuine order move":** (L) green and every m ≤ −90 %. Predicted improbable: the field-side loss changes only ≈ −10 % between orders at κ = 0.
  - **"Partial":** (L) green, and any \|m\| between 30 % and 90 % or m > +30 %.
  - **"Not κ, port extraction":** (L) red at degree 2. W0 stops the item.
  - **"Mixed":** rows or columns in different categories. Report and stop.
  - Printed beside each row: the degree-1 value, the held-κ value, both moves, and the 1/Q split.
- **No `PREDICTED RED` line.** A red anywhere stops the slot.
- **Windows** (complex, in order, each < 10 min, so no durable capture, as item 132 ran).
  - **W0a** (knob-unset regression, ≈ 74 s): `20261001T125801Z_ANS-7-step3d2-W0b.log:12` verbatim, relabelled `ANS-7-step3iii-W0a`. It must print 1.136e-14 again (`:997`).
  - **W0b** (degree-1 κ route, positive control for (L) and (N) at degree 1, ≈ 120 s): `scripts/testing/run_and_log.sh ANS-7-step3iii-W0b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_SPARSE_MHZ=64.5 FEM_EM_ANS7_SPARSE_DEGREE=1 FEM_EM_ANS7_SPARSE_KAPPA=in_run timeout -k 30 400 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_ans7_sparse_route_b.py::test_wiring_control tests/validation/test_ans7_sparse_route_b.py::test_kappa_in_run_discriminator tests/validation/test_ans7_sparse_route_b.py::test_node_anchors -s -v --tb=short'"`
  - **W0** (the degree-2 discriminator, ≈ 420 s): two degree-2 calls at 183.9 s each, 13.78 GiB summed (`20260921T100022Z_…:942`), plus ≈ 48 s. `scripts/testing/run_and_log.sh ANS-7-step3iii-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_SPARSE_MHZ=64.5 FEM_EM_ANS7_SPARSE_DEGREE=2 FEM_EM_ANS7_SPARSE_KAPPA=in_run timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_sparse_route_b.py -s -v --tb=short'"`. W0's 64.5 MHz node JSON is the grid's 64.5 node.
  - **W1** (≈ 421 s): W0's command with `FEM_EM_ANS7_SPARSE_MHZ=62.0,62.3,63.2`, `FEM_EM_ANS7_SPARSE_FIELD_SOLVE=0` and `FEM_EM_ANS7_SPARSE_KAPPA=stored`, labelled `ANS-7-step3iii-W1`. **Window-edge guard:** if Im Y_in at 62.0 MHz is > 0 in either column (f₀ below the window), stop before W2.
  - **W2** (≈ 428 s): the same with `=65.8,66.7,67.0`, labelled `ANS-7-step3iii-W2`. It holds 67.0, so the assembly over all seven `_kin` JSONs runs here.
  - **Total:** ≈ 1 463 s of windows.
- **Traps.**
  - **`pec_column=False` together with `field_solve=True` raises `KeyError` on rank 0 only** (`rec["pec"]` at `test_th14_birdcage_copper.py:401`; `class_ds` / `max_abs_ds` at `:417–420`). The other ranks wait and the window burns to the timeout. Found by the drafter while reading; not yet paid for. Never combine `_COLUMNS=cu` with `_FIELD_SOLVE=1` or `_KAPPA=in_run`.
  - The κ pass must use `width_correction_kappa=None`: a corrected pass returns a residual deficit, not κ. `rec["kappa"]` exists only inside the field-solve branch (`:381`).
  - Item 132's held-κ JSONs share the p2 names. The `_kin` suffix keeps them apart, and the `:447` κ check catches any mix.
  - `stored` reads only after rank 0's barrier. Write κ with `repr`.
  - Never import `_assemble`, `_write_outputs` or `_stage_*` from `ans:7`. `git diff --stat` shows no change to `metrics.json` or `COMPARISON.md`. Every `FEM_EM_ANS6_*` knob stays unset.
  - Degree 2 at `-n 8`; degree 1 at `-n 4`. No `-k`; node ids only.
  - Harness bare and whole, foreground, `timeout: 660000`, `-s`. No edit after W0a (rule (i)).
  - A window killed at 590 s is not a result: `pgrep -c python3` first, then report. Never lengthen the timeout.
  - Sweep only 0-byte `.c` stubs.
- **Scope.**
  - This is a record: no glyph (T13 decides), no band moved, no `src/` change, no TH-14 edit.
  - Not covered: h = 0.005 (T12 / entry 24), Route A at degree 2, F12, `PORT-21`'s fixture, T6, AED.
  - Item 132's (ii) record stands as written.
- **Status it can move.**
  - The `ANS-7` row's step-3 order-move sentence: the review replaces "not established as an order move" with (iii)'s reading, κ_deg2 and the in-run rows.
  - **On a green item** (any reading but "not κ"), the implementer marks `xl` entry 24 "PREREQUISITE LANDED <date> — knob `FEM_EM_ANS7_SPARSE_KAPPA=in_run` exists; command to be issued by the next review" in the same commit (daily-review step 6b.5). The implementer never writes an `xl` command. κ is mesh-dependent, so entry 24 must run `in_run`, never `stored`, and its κ pass must include the PEC sweep (the trap above).
- **Negative result.**
  - **W0a not digit-identical:** revert, BLOCKED.
  - **W0b red:** stop, BLOCKED; no degree-2 window.
  - **W0 (R2) red:** known-issues entry, stop.
  - **W0 (N) red:** report, stop.
  - **W0 (L) red (the "not κ" reading):**
    - print the gap, κ_deg2 and the column powers;
    - park on `attempt/*`, BLOCKED (rule (d));
    - known-issues entry "degree-2 lumped-sheet S not power-consistent at a degree-matched κ";
    - entry 24 stays HELD.
  - **Window-edge guard fires, or a column has ≠ 1 crossing:** BLOCKED; the review re-poses the grid.
  - **Δ_LOO above the bands:** rows printed as unlicensed; the review re-poses the grid.
  - In every case, report and stop. Never loosen a band.

**Outcome 2026-10-02 (19:30 CDT slot): DONE — reading "κ confound".**
- **Windows.** All Status 0:
  - W0a `20261002T003944Z_…-W0a.log` (75 s): 1.079e-14.
  - W0a re-run after a lazy-import edit, `20261002T004136Z_…-W0a.log` (78 s): 1.157e-14 at `:997`. This is the closing W0a; the module mtime precedes it.
  - HEAD baseline `20261002T004311Z_…-W0a-baseline-HEAD.log` (30 s): the unmodified HEAD module, run from a deleted temporary copy, prints **5.410e-14** (`:933`).
  - W0b `…004349Z` (116 s), W0 `…004553Z` (400 s), W1 `…005242Z` (384 s), W2 `…005913Z` (392 s).
- **κ.** κ_deg1 1.0599570850048146e-02, with (K1) Δ = 0. κ_deg2 = 5.2094420914264195e-02 (W0 `:1124`), inside the predicted [4.5e-2, 6.0e-2].
- **Anchors.**
  - (R2): 2.3e-10 / 6.1e-11.
  - (L): gap 1.744e-6 at p2 (W0 `:2092`) and 7.085e-7 at p1.
  - (N): 9.829e-2 (98×) at p2 and 2.169e-2 (21.7×) at p1.
  - Column power: Cu 0.998714, PEC 0.999322.
  - Per node: σ_max ≤ 0.99993, reciprocity ≤ 3.8e-14, identity ≤ 5.2e-13.
  - Window edge: Im Y_in(62.0) < 0 in both columns.
- **Grid** (W2 `:1324–1351`), p2 in-run vs p1:

  | Column | f₀ (MHz) | R_in (Ω) | Q |
  |---|---|---|---|
  | Cu | 62.704216 vs 63.983031 (−2.0 %) | 18 570 vs 16 999 (+9.2 %) | 984.2 vs 918.4 (+7.2 %) |
  | PEC | 62.717977 vs 63.997103 (−2.0 %) | 32 545 vs 28 425 (+14.5 %) | 1 724.8 vs 1 536.0 (+12.3 %) |

  k-share ≥ 0.999995; Δ_LOO ≤ 8.0e-10. The copper part of 1/Q is 4.36e-4 at p2 against 4.38e-4 at p1. The phantom field loss at κ = 0 moves −14.6 % from p1 to p2 (≈ −15 % predicted). Every printed prediction landed in its range.
- **Deviation flagged for the review.** The brief's "W0a not digit-identical: revert, BLOCKED" letter was not met: 1.157e-14 against 1.136e-14. The HEAD module reproduces neither figure (5.41e-14), consistent with the carried ruling that records tighter than ≈ 5e-11 drift run to run. The asserted ≤ 1e-9 check is green in every run. The slot landed the item and asks the review to ratify or revert. `xl` entry 24 is marked PREREQUISITE LANDED, with the same flag.

**Ruled 2026-10-02 03:00 review — ratified; the landing stands.** The brief's "digit-identical" letter named a figure (1.136e-14) that is a round-off residual, not a record. The unmodified HEAD module prints 5.41e-14 on the same configuration (`20261002T004311Z_ANS-7-step3iii-W0a-baseline-HEAD.log`), so no run of any version reproduces a 1e-14-level residual to its digits. That is the carried ruling: records tighter than ≈ 5e-11 drift from run to run. The asserted wiring check, node S vs stored ≤ 1e-9, is green in every run, and its quantity and band are unchanged. The edit between the two W0a runs (lazy imports) is behaviour-neutral, and the closing W0a ran on the committed module (rule (i)). **Lesson for briefs:** a "default digit-identical" check names a record at ≥ 1e-11 or the asserted band, never a round-off residual. Separately, `xl` entry 24 is still not issuable: the module's assembly needs all seven nodes and both columns. At h = 0.005 and degree 2 that is ≈ 14 factorisations, over the 4 h cap. The five-node Cu-only assembly path is §9 item 173 (below).

### §9 item 173 — 2026-10-02 (chain T12 prerequisite — `xl` entry 24's assembly path; queued by the 03:00 review)

`ANS-7` step 4 prerequisite: a five-node, Cu-only, in-memory assembly path, so that entry 24 fits one `xl` window. Executor: implementer, tests only. Tier: heavy; `-n 4` (W0a, W2), `-n 1` (W1), `-n 8` (W3); complex. Cost: ≈ 17 measured / ≈ 60 predicted min. Depends on nothing. Drafted by a read-only `Plan` agent and reviewed by the 03:00 review. The review adopted the drafter's memory guard (doubt 3) and confirmed W0a's identity wording, which is consistent with the item-157 ruling.

- **Why.** Entry 24 at seven nodes × two columns plus the κ pass is 17 factorisations, ≈ 22 000 s. Seven Cu nodes is ≈ 14 440 s at the top of the range. Both exceed the 14 400 s cap. The module also cannot assemble a Cu-only run today:
  - `in_run` refuses `cu` (`:330–332`).
  - Node JSON is written only for both columns (`:387`).
  - Assembly needs 67.0 MHz, both columns and all seven JSONs (`:587–594`).
  - `_subset`, `_loo` and `_kin_rows_report` hard-code `("cu","pec")`.
- **Code** (`tests/validation/test_ans7_sparse_route_b.py` only; nothing imports it; no `src/`, no `TH-14` edit).
  1. `ASSEMBLY5_NODES_MHZ = (62.0, 63.2, 64.5, 65.8, 67.0)`: five of the seven defaults, keeping both window ends and 64.5 MHz, so `CONTROL_NODES_MHZ` ⊂ it. `SUBSET5_NODES_MHZ` stays untouched: its 62.7 / 66.3 MHz nodes are solved nowhere.
  2. Knob `FEM_EM_ANS7_SPARSE_ASSEMBLY`.
     - Unset: today's behaviour.
     - `solved`: a new `test_assemble_solved` assembles this run's own nodes and columns in memory. 62.0 and 67.0 MHz are required, otherwise `ValueError`. `test_assemble_from_node_json` skips.
     - Anything else: `ValueError`.
  3. Under `solved` only:
     - The guard at `:330` becomes `"pec" not in columns and field_solve`. That combination is exactly the rank-0 `KeyError` trap, and it still raises.
     - The κ pass is unchanged: uncorrected, with its PEC sweep and field solve.
     - Rank 0 writes Cu-only node JSON with a `_kin_cu` suffix, then a barrier.
     - **After the κ pass, `del klad` and `gc.collect()`** inside `_kappa_in_run`, before the grid factorises. Otherwise a held κ-pass factor plus the grid's would peak at ≈ 530 GiB against the 512 GiB limit. No number changes.
  4. `_subset`, `_loo` and `_kin_rows_report` take `columns=("cu","pec")` as the default, so the unset path is byte-identical. The 1/Q-split lines print only when PEC is present.
  5. `test_assemble_solved` **asserts** exactly one − → + crossing and k-share > `K_SHARE_FLOOR` per column. It **prints**:
     - the rows and the order move against p1;
     - `_kin_rows_report`;
     - Im Y_in(62.0);
     - Δ_LOO over the set's interior nodes and B₅ = ceil_1sf(max(3·Δ_LOO, 1e-6)), printed, not asserted: dropping a node from five measures four-node interpolants;
     - the 3-node control rows;
     - at resolution 0.005, the h-move against item 157's ×1 in-run p2 Cu rows (a new constant from `20261002T005913Z_ANS-7-step3iii-W2.log:1325–1327`).

     Where a non-`_cu` JSON set exists at the same degree, resolution and κ (to 1e-15), it also asserts that the solved Cu rows equal the JSON-assembled rows within `NODE_S_RTOL`.
  6. `test_kappa_in_run_discriminator`.
     - With `FIELD_SOLVE=0` it runs (K1), (R2) and (N) from the κ pass and skips (L) with a print.
     - The (R2) reference is keyed by resolution:
       - `None`: today's reference;
       - `0.005`: entry 14's uncorrected line, P_phantom 1.696818288e-06 W and P_surf 1.431230505e-06 W (`20260929T070008Z_ANS-6-step2c.log:1034`; same ladder, 64 MHz, `-n 16`);
       - any other value: prints "no reference".
  7. `test_assembly_subset_anchor`: pure, no fixture; runs only with the knob set.
- **Anchor (asserted).** Take the Cu column from the seven stored `f_*_p1_hdef.json`. Its five-node rows reproduce the seven-node rows for f₀, R_in and Q: \|row₅ − row₇\|/\|row₇\| ≤ `BANDS` (1e-6, frozen from W1 by `TH-18` step 2's rule).
  - Repeat on the seven `f_*_p2_hdef_kin.json` (κ 5.2094e-2).
  - Predicted: ≈ 2e-9 at p1 and ≈ 7e-8 at p2. The p2 figure is scaled from its six-node Δ_LOO, 8.0e-10 at `…005913Z_…-W2.log:1294`.
  - All fourteen JSONs sit under the gitignored `logs/ans7_sparse/`. W0a regenerates the p1 set. If the p2 set is absent, stop.
- **Negative control (rule (e)).** The 3-node control's Cu rows miss B on R_in or Q.
  - **Asserted at p1**, backed by 1.915e-6 (1.9×) and 1.311e-5 (13×) at `20261001T095802Z_ANS-7-step3-W2.log:1548`.
  - **Predicted at p2**, printed.
  - Ceiling: this is an interpolation error, not a saturating fraction. The thin side is 1.9×, so no factor is asserted, only "misses".
- **Windows** (each < 10 min; bare; foreground; Bash `timeout: 660000`).
  - **W0a** (≈ 157 s): `…095802Z_ANS-7-step3-W2.log:12` verbatim, knob unset, labelled `ANS-7-item173-W0a`.
    - The pass set must be the same, with the new tests skipped.
    - The 3-node control and LOO prints at `:1482–1549` must reproduce at their printed digits for every figure ≥ 1e-6. Smaller figures drift, per the item-157 ruling.
  - **W1** (≈ 30 s): `ASSEMBLY=solved`, `mpiexec -n 1`, node `…::test_assembly_subset_anchor`.
  - **W2** (≈ 140 s): `scripts/testing/run_and_log.sh ANS-7-item173-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_SPARSE_MHZ=62.0,63.2,64.5,65.8,67.0 FEM_EM_ANS7_SPARSE_DEGREE=1 FEM_EM_ANS7_SPARSE_COLUMNS=cu FEM_EM_ANS7_SPARSE_FIELD_SOLVE=0 FEM_EM_ANS7_SPARSE_KAPPA=in_run FEM_EM_ANS7_SPARSE_ASSEMBLY=solved timeout -k 30 400 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_ans7_sparse_route_b.py -s -v --tb=short'"`
    - Checks: (K1); (N) at p1 (21.7×, backed); solved rows = p1 JSON rows within 1e-9.
  - **W3** (≈ 340 s; item 157's W2 ran six sweeps in 342 s): W2's command with `DEGREE=2`, `KAPPA=stored`, `-n 8`, `timeout -k 30 590`.
    - Check: solved rows = p2_kin JSON rows within 1e-9.
    - It uses `stored` because `in_run` would take ≈ 520 s against 590.
- **Traps.**
  - No ordering comparison or `sqrt` on complex operands.
  - `-s` on every window.
  - Rows are computed on every rank; rank 0 writes, then a barrier. `assemble_scalar` is rank-local.
  - `_COLUMNS=cu` with `_FIELD_SOLVE=1` gives a rank-0-only `KeyError`, so the guard keeps that combination raising.
  - The κ pass stays uncorrected (`width_correction_kappa=None`).
  - Harness calls bare and whole, foreground; no `-k`, node IDs only.
  - `FEM_EM_ANS_SKIP_PRIVATE=1` on every window, and no `FEM_EM_ANS6_*` knob.
  - Never import `_assemble`, `_write_outputs` or `_stage_*`. `git diff --stat` must show no change to `metrics.json` or `COMPARISON.md`.
  - An edit after any window re-runs W0a; make no edit after the last window.
  - A killed window is not a result: run `pgrep -c python3`, then report.
  - **No slot writes or edits an `xl` command.**
- **Scope.** A record: no band moved, no glyph. The five-node licence is measured at the default mesh (p1, p2) and carried to h = 0.005 by assumption; the `xl` window's B₅ print is its check. Not covered: the PEC column at h = 0.005, Route A, AED.
- **Status it can move.** `xl` entry 24 PENDING PREREQUISITE → READY. On green, the implementer marks it in the same commit (daily-review step 6b.5) and leaves the command alone. The `ANS-7` §7 row is unchanged.
- **Negative result.**
  - W0a not identical: revert and mark BLOCKED.
  - Anchor red at p1, or the control inside B: BLOCKED.
  - Anchor red at p2 only: report the gap. The review then picks six Cu nodes (≈ 12 100–13 000 s, still `xl`) or seven on `xxl`.
  - W2 or W3 red: BLOCKED.
  - In every case, report and stop. Never loosen `BANDS`.

**Outcome 2026-10-02 (06:00 CDT slot): DONE (a record, no glyph); `xl` entry 24 → READY.** Tests only (`test_ans7_sparse_route_b.py`, +254 / −; no `src/`, `metrics.json` / `COMPARISON.md` untouched). The module predates W0a (06:14:48 vs 06:14:59 CDT).

- **W0a** `-n 4`, knob unset (`20261002T111459Z_ANS-7-item173-W0a.log`, 160 s, Status 0):
  - 18 passed, 3 skipped. That is the reference's 15 plus three tests added by items 132 / 157. The skips are `test_kappa_in_run_discriminator` and the two new tests.
  - The 3-node control reproduces digit for digit: Cu 1.005e-07 / 1.915e-06 / 1.311e-05 and PEC 9.794e-08 / 2.718e-06 / 1.390e-05 (`:1559–1560`). Every LOO figure is ≤ 2.3e-11, so all of them sit under the 1e-6 reproduction rule.
- **W1** `-n 1`, `test_assembly_subset_anchor` (`…111756Z_ANS-7-item173-W1.log`, 3 s): 1 passed. The anchor is \|row₅ − row₇\|/\|row₇\| on f₀ / R_in / Q:
  - p1: 7.306e-13 / 1.064e-10 / 1.465e-10.
  - p2_kin (κ 5.2094e-2): 1.233e-12 / 1.267e-10 / 2.780e-10.
  - Both are ≤ `BANDS` 1e-6. The predictions were ≈ 2e-9 and ≈ 7e-8, so the measured anchor is smaller.
  - The 3-node Cu control misses. At p1, asserted: 1.005e-07 / 1.914e-06 / 1.311e-05. At p2, predicted and printed: 1.705e-07 / 3.474e-06 / 2.240e-06. Q misses at p2 by only 2.2×.
- **W2** `-n 4`, degree 1, Cu, κ in-run, `ASSEMBLY=solved` (`…111807Z_ANS-7-item173-W2.log`, 127 s): 17 passed, 4 skipped.
  - Checks: (K1) 4.441e-16. (N) 2.168874e-02, which is 21.7× the band. (L) is skipped with a print.
  - Solved rows equal the p1 JSON rows to 3.8e-15 / 1.6e-12 / 3.0e-12.
  - Rows: f₀ 63.983031 MHz, R_in 16 998.71 Ω, Q 918.441. One − → + crossing; k-share 0.99999385.
  - Five-node Δ_LOO is 3.5e-10 / 1.8e-8 / 4.7e-8, so B₅ = 1e-6.
- **W3** `-n 8`, degree 2, κ stored (`…112024Z_ANS-7-item173-W3.log`, 310 s): 15 passed, 6 skipped.
  - Solved rows equal the p2_kin JSON rows to 4.2e-14 / 2.8e-11 / 3.4e-11.
  - Rows: f₀ 62.704216 MHz, R_in 18 569.61 Ω, Q 984.223. The order move against p1 is −1.999 % / +9.24 % / +7.16 %. k-share 0.99999501.
  - Five-node Δ_LOO is 6.5e-10 / 3.6e-8 / 2.8e-8, so B₅ = 1e-6.
  - The h-move is not printed, because it fires only at resolution 0.005, which is the `xl` window's job.
- **Disclosed.**
  - The solved-vs-JSON κ match "to 1e-15" is implemented as an *absolute* tolerance. A relative one would have excluded the p1 set, because the κ difference is 4.4e-16 on κ ≈ 1.06e-2.
  - `test_assemble_solved` reuses `_kin_rows_report`'s reading-rule line ("κ confound …"), which means nothing on a Cu-only or p1 run. It is cosmetic and was left in place because no edit is allowed after the last window.
  - W2 rewrote the gitignored `logs/ans7_sparse/kappa_p1_hdef.json` (κ unchanged to 4.4e-16) and added `*_kin_cu.json` node files.

### §9 item 195 — 2026-10-03 (chain T, T17 — daily provisional, rule (6))

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.* Licence: rule (6), toward chain T's subgoal (the tuned F-small resonance read on a rung, T13) and the operator directive 2026-10-03 (a). The sheets-and-κ ruling (§9, 2026-10-02) applies: the port **and** RLC sheets are exact rectangles, and κ is off. No band, tolerance, Done-when or validation target moves.

`ANS-7` step 5: the tuned fixture on exact port and RLC rectangles at degree 1 — Route A = Route B, and the rows re-recorded.
- **Executor:** implementer. Tests only, in a new module.
- **Tier:** heavy, `-n 4`, complex build.
- **Class:** chain step (solve).
- **Time:** ≈ 14 measured / ≈ 40 predicted min.
- **Depends on item 193 landing:** the hole-route `leg_sheet_strip_fraction` and its two helper keywords. If 193 did not land or is BLOCKED, skip this item and item 196, and take the next item.

**Why: what the sheets are today.**
- The SPEC's four gap sheets are `ANS-4`'s rectangle, 0.007 × 0.008 (`examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/SPEC.md:77–78`). Gap 1 carries the 50 Ω port; gaps 2–4 carry Lumped RLC capacitors, C = 66.29 pF.
- **Our side builds both kinds of sheet with one call.**
  - Route B (`_build_ladder` → `_build_ports`) and Route A (`07_…py:690–700`, `_route_a_solve` `:515–545`) both use `th14._build_ports` on the hole mesh.
  - In Route A, the capacitor sheets are the same narrowed N_2…N_4 with `port_impedance_ohm = Z_C` (`:527–529`).
  - So the ragged sheet geometry of item 193 is also the RLC sheet geometry, and one keyword fixes both.
- **The κ evidence.**
  - Today's rows need κ = 1.0599570850e-02 on every sheet. Without it, Route A misses Route B by 16.6× `REDUCTION_BAND` at 64 MHz Cu (control (ii), 1.657e-02, `20260929T100040Z_ANS-7-item68-W3-A.log`; the 13.5× of `docs/planning/chunks/PORT-14.md:650`).
  - On exact rectangles κ falls to ≈ 1.6e-5 (items 159 / 166). The prediction: Route A = Route B holds **with κ off**.
  - That is the ruling's claim, tested on the lumped-element path for the first time.

**Code** (tests only; new `tests/validation/test_ans7_exact_strip.py`; nothing imports it; no `src/`, no `TH-14` edit, no edit to `ans:7`):
- `@complex_only`. Knob `FEM_EM_ANS7_EXACT_ROUTE` ∈ {`B`, `A`}; unset, every test skips. Also `FEM_EM_ANS7_EXACT_DEGREE` (default 1), which item 196 uses.
- **Imports, never restated:**
  - from `tests/validation/test_ans7_sparse_route_b.py`: `DEFAULT_NODES_MHZ`, `GRID_MHZ`, `BANDS`, `N_FACTOR`, `NODE_S_RTOL`, `K_SHARE_FLOOR`, `sparse_rows`, `bary_interp`, `z_from_s`, `_ANS7` (the loaded `ans:7` script, hence `route_b_s`, `route_b_z`, `_route_a_solve`, `_kappa_in_run`, `_s11_from_z`, `parallel_rows`, `add_c1`, `C_F`);
  - `REDUCTION_BAND` (`tests/validation/test_port_circuit_layer_field.py`);
  - `GATED_WIDTH_FRACTION`;
  - `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` (item 193);
  - `CELL_COUNT_BAND`, `STRIP_BAND`, `RECIPROCITY_BAND`, `PASSIVITY_SIGMA_TOLERANCE`.
- **Route `B`.**
  - The call is `lad = th14._build_ladder(d, sigmas=(th14.COPPER_SIGMA,), frequencies_hz=[m*1e6 for m in nodes], width_correction_kappa=None, measure_kappa=False, pec_column=True, field_solve=False, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`, with nodes = `DEFAULT_NODES_MHZ` (the 7-node grid licensed at both degrees by step 3 (i) / (ii)).
  - Per node, Z = `z_from_s(S)` for the Cu and PEC columns; then `sparse_rows(DEFAULT_NODES_MHZ, z_nodes)`, giving f₀ / R_in / Q / k-share per column.
  - Rank 0 writes `logs/ans7_exact/route_b_p{d}.json` (gitignored), containing the node S, Z, rows and cells; then a barrier.
- **Route `A`** (degree 1 only).
  - Build the strip hole mesh (`th14._sheets_build(True, as_hole=True, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`), then `_outer_box_tags` and `_build_ports` with κ `None`.
  - Run `_route_a_solve` at all 51 `GRID_MHZ` points × {cu, pec}.
  - The comparand is Route B's bare S11 on the same grid: read W1's JSON, interpolate the network Z entrywise (`bary_interp`, as `sparse_rows` does), take `route_b_z(…)[1]`, then `_s11_from_z`.
  - Also `_kappa_in_run` (64 MHz Cu, uncorrected) on the exact sheets, printed only.
  - **The ragged control:** on the default hole mesh, one uncorrected `_route_a_solve` at 64 MHz Cu against `route_b_s(tracked metrics["check_64mhz_uncorrected"]["cu_s"], 64e6)`. This is `_stage_a`'s control (ii) verbatim (`07_…py:734–737`).
- **Record constant.** `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1 = None` until W1, then registered from W1's prints. It holds f₀ / R_in / Q per column at full repr, with a version-tag comment and the (1*) note.

**Anchors** (asserted; bands imported, never restated):
- **(F) Fixture tie** (W1, W2):
  - cells equal `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` within `CELL_COUNT_BAND`;
  - on all four sheets, N_i area = 5.6e-5 and w = 7.0e-3 within `STRIP_BAND`.

  This covers the port sheet and the three RLC sheets.
- **(B) Route B, per node:**
  - reciprocity ≤ `RECIPROCITY_BAND`;
  - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`;
  - exactly one − → + zero of Im Y_in per column on the interpolated grid;
  - k-share > `K_SHARE_FLOOR`;
  - the 7-node Δ_LOO on f₀ / R_in / Q ≤ `BANDS`. This is the grid's licence, re-checked on the new sheets.
- **(AB) Route A = Route B,** the row's v2 identity: \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` at **all 102** points (51 per column), with κ off on both routes. Also exactly one − → + zero of Im Y_in per column on Route A.
- **(R)** At `-n 4` only (the record width): Route B's rows equal `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1` within `BANDS`. Skipped while the record is `None`.

**Negative control (rule (e)).**
- **(N) Ragged κ off — asserted.** The control (ii) residual exceeds `N_FACTOR` × `REDUCTION_BAND`.
  - Backed on this fixture, comparison and mesh by 1.657e-02 = 16.6× (item 68 W3).
  - The ceiling is \|ΔS\|/\|S\| ≲ 2, so no saturation.
  - The exact-sheet counterpart of the same comparison is the 64.0 MHz Cu point of (AB), printed beside it.

**Predictions** (written first, printed):
- **Exact cells:** item 193's record.
- **κ_exact (64 MHz Cu):** in [5e-6, 5e-5].
- **(AB) residual:** in [5e-6, 3e-4] at every point. This is the ragged κ-off 1.657e-02 scaled by κ_exact / κ_ragged ≈ 1.5e-3, gives ≈ 2.4e-5, plus the corrected-ragged floor of 3.4e-5 to 9.3e-5.
- **Rows, both routes, both columns:**
  - \|Δf₀\|/f₀ against the ragged κ-corrected rows ≤ 1 % (Cu 63.983031 / PEC 63.997103 MHz, `metrics.json`). The sign is not pre-registered.
  - Q and R_in within ±15 % of 918.44 / 1535.95 and 1.6999e4 / 2.8425e4 Ω.
  - k-share ≥ 0.9999. Δ_LOO ≤ 1e-9.
- **Printed beside each row:** the ragged degree-1 value and the relative move. The f₀ move is the reading the weekly takes privately against the known-issues 2026-09-30 gap. **No AED value appears here.**

**Windows** (bare, whole, foreground, Bash `timeout: 660000`, `-s`, `tests/environment` first, `FEM_EM_ANS_SKIP_PRIVATE=1` on every window).
- **W1** (Route B, ≈ 110–160 s). Step 3 (i)'s 7-node degree-1 window ran 157 s at `-n 4` *with* the field solve (`20261001T095511Z_ANS-7-step3-W1.log`).

  `scripts/testing/run_and_log.sh ANS-7-item195-W1-routeB "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_ROUTE=B timeout -k 30 400 mpiexec -n 4 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`

  Then register `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1`.
- **W2** (Route A, ≈ 330–420 s). `ans:7` stage A ran 102 drives, plus κ, the controls and XDMF, in 367 s at `-n 4` (`20260929T100040Z_ANS-7-item68-W3-A.log`). Add one ragged hole build (≈ 20 s) and one drive.
  - Command: W1's, with `FEM_EM_ANS7_EXACT_ROUTE=A` and `timeout -k 30 590`, labelled `ANS-7-item195-W2-routeA`.
  - It runs on the module as committed, after the registration.
- **W3:** W1 relabelled `ANS-7-item195-W3-routeB`. This is the closing window (rule (i)), and (R) is asserted. Make no edit after W2.

**Traps.**
- Never combine `pec_column=False` with `field_solve=True`: rank 0 alone raises `KeyError` and the window burns (item 157). Route B here runs with the field solve off and the PEC column on.
- Every sheet spec keeps `width_correction_kappa=None` (ruling (ii)). κ is read only through `_kappa_in_run` on uncorrected specs.
- Never call `ans:7`'s `_stage_*`, `_assemble` or `_write_outputs`. `git diff --stat` must show no change to the tracked `metrics.json` or `COMPARISON.md`. Set no `FEM_EM_ANS6_*` variable. Never run `ans:1`–`ans:3`, and do not run `ans:7`.
- Write only under `logs/ans7_exact/`. Never write into `logs/ans7/` or `logs/ans7_sparse/` (item 173's and the `xl` entry's files).
- Rows are computed on every rank; rank 0 writes, then a barrier. `assemble_scalar` is rank-local.
- No ordering comparison or `sqrt` on complex operands.
- On a kill: run `pgrep -c python3`, report. Never lengthen the timeout.

**Scope.** Covered: degree 1, h = 0.015, Cu and PEC columns, both routes, κ off. Not covered:
- degree 2 (item 196);
- the h = 0.005 rung (`xl` entry 24, the weekly's to re-issue on exact sheets);
- T14 / T15 (`WF-5`);
- the f₀ / Q re-read against the private columns (the weekly's, directive (a));
- any band, and the example's tracked artefacts.

**Status it can move.**
- The `ANS-7` row gains a "step 5" record sentence. No glyph: the row's Done-when is T13's re-read, and that stands.
- Item 196 becomes runnable.
- **Rule (j):** the review adds a step-5 Done-when sentence.

**Negative result.**
- **(F) red:** stop, BLOCKED.
- **(B) red, or a column with ≠ 1 crossing:** BLOCKED; the review re-poses the grid.
- **(AB) red with (N) green:** the ruling's prediction fails on the RLC path. Print the residual curve and κ_exact, open a known-issues entry ("exact-rectangle Route A ≠ Route B at κ = 0"), and park on `attempt/*`, BLOCKED. Never apply κ to rescue it.
- **(N) red:** the ragged control did not reproduce. Report and stop.
- No `PREDICTED RED:` line.

**Why there is no Palace twin of this fixture.** It is not queued as a chain-T step. `OPS-61.md` has exercised Palace only with 50 Ω lumped ports on the σ = 800 solid fixture. `ANS-7`'s twin needs two Palace features that have never been run here:
- the lumped R / L / C termination of non-excited ports, plus a check of its port-S convention with a capacitive port;
- the finite-conductivity (surface-impedance) boundary on tag 401 for the copper column, which is `OPS-61` step 2.

It also needs a Gmsh export of the hole mesh that does not yet exist. Each is a capability probe for the weekly to place, not a step toward T13. The nearest cheap Palace leg on this fixture is a PEC-column twin of item 193's mesh (PEC on tag 401, 50 Ω ports), which needs only that export.

**Outcome (2026-10-03, 06:00 CDT slot): DONE, all green.** New `tests/validation/test_ans7_exact_strip.py` (tests only; no `src/`, `TH-14` or `ans:7` edit). W1 Route B (`20261003T112739Z_ANS-7-item195-W1-routeB.log`, 164 s, -n 4): (F) 80 798 cells = item 193's record, all four N_i 5.6e-5 / w 7.0e-3 to ≤ 1.6e-15; (B) reciprocity ≤ 1.01e-15, σ_max ≤ 0.99991, one − → + crossing per column, k-share 0.999993 / 0.999998, Δ_LOO 2.99e-14 / 9.44e-12 / 9.99e-12 (f₀ / R_in / Q). Rows (registered as `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1`): Cu f₀ 63.971138 MHz, R_in 16 907.70 Ω, Q 912.717; PEC f₀ 63.985193 MHz, R_in 28 088.54 Ω, Q 1516.439 — moves against the ragged κ-corrected rows f₀ −0.019 % (both), R_in −0.54 % / −1.18 %, Q −0.62 % / −1.27 %, all inside the predictions. W2 Route A (`20261003T113037Z_ANS-7-item195-W2-routeA.log`, 546 s): κ_exact 6.583e-6 (predicted [5e-6, 5e-5], met); (AB) residual 9.9e-6 … 3.11e-5 (Cu) and 9.9e-6 … 3.09e-5 (PEC) at all 102 points, inside the predicted [5e-6, 3e-4]; one crossing per column on Route A (Cu f₀ 63.970953 MHz); (N) ragged κ-off control (ii) 1.656839e-02 = 16.57 × `REDUCTION_BAND` (backing 1.657e-02 reproduced), beside the exact-sheet 64.0 MHz Cu (AB) point 1.409e-5. W3 closing Route B (`20261003T113959Z_ANS-7-item195-W3-routeB.log`, 161 s): 15 passed, 3 skipped, (R) asserted. The ruling's prediction holds on the lumped-element path: with exact rectangles Route A = Route B with κ off. Item 196 is runnable.

### §9 item 196 — 2026-10-03 (chain T, T18 — daily provisional, rule (6))

*Drafted by a general-purpose agent for the 2026-10-03 03:00 daily review; the review reviewed it.* Licence: rule (6), toward chain T's T13 re-read (the order-matched rung). The ruling applies: exact rectangles, κ off. No band, tolerance, Done-when or validation target moves.

`ANS-7` step 5b: item 195's Route B at degree 2 on exact rectangles, κ off, plus the loss identity at 64.5 MHz.
- **Executor:** implementer. Tests only, in item 195's module.
- **Tier:** heavy, `-n 8`, complex build.
- **Class:** chain step (solve).
- **Time:** ≈ 14 measured / ≈ 45 predicted min.
- **Depends on item 195 landing.** If it did not, skip this item.

**Why.** The ragged degree-2 history:
- Item 132 held κ at its degree-1 value: R_in and Q fell 98 %, and the column power was 0.919.
- Item 157 found the cause, a "κ confound": with the in-run degree-2 κ = 5.2094e-02, the loss identity (L) closed to 1.744e-6. The degree-2 rows came out Cu f₀ 62.704216 MHz, R_in 18 570 Ω, Q 984.2.

On exact rectangles κ at degree 2 is ≈ 1.8e-5 (item 166). The degree-2 rows should therefore stand **without any correction**, and (L) is the identity that says so.

**Code** (item 195's module only; disclosed; no `src/`):
- `FEM_EM_ANS7_EXACT_DEGREE=2` with `ROUTE=B`, plus a node-subset knob `FEM_EM_ANS7_EXACT_MHZ` (a comma list ⊂ `DEFAULT_NODES_MHZ`). Rank 0 writes one JSON per node, `logs/ans7_exact/node_<mhz>_p2.json`. Assembly through `sparse_rows` runs once all seven node JSONs exist with the same degree and cells. This is item 132's split.
- **The 64.5 MHz node:** when it is in the subset, that node's call runs with `field_solve=True` (and `pec_column=True`), and the test asserts (L).
- **Knob `FEM_EM_ANS7_EXACT_HELDKAPPA=1`:** at 64.5 MHz only, one more Cu sweep plus P1 field solve with `width_correction_kappa=5.2094420914264195e-02`, the ragged degree-2 κ of item 157 W0 `:1124` (restated with its cite). This is the negative control.
- Degree-2 rows are printed only; no constant is registered (a follow-up for the review).

**Anchors** (asserted; bands imported):
- **(L) at 64.5 MHz, degree 2, κ off.** \|L_S − L_F\| ≤ `REDUCTION_BAND`.
  - L_S = 1 − Σ\|S_i1\|² from the Cu S.
  - L_F = (P_phantom + P_(Ω∖phantom) + P_surf)/P_inc from the same call's P1 field.
  - Use `_loss_pair` from `test_ans7_sparse_route_b`. P_inc = 2.5e-3 W.
- **Per node:** reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and the copper identity residual ≤ `DISCRETE_IDENTITY_RTOL` at 64.5 MHz.
- **Assembly:**
  - exactly one − → + zero per column;
  - k-share > `K_SHARE_FLOOR`;
  - degree-2 Δ_LOO ≤ `BANDS`;
  - (F) cells = `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, carried in every node JSON.

**Negative control (rule (e)), predicted and printed.** The same (L) with the ragged κ_deg2 applied to the exact sheets misses by a gap in [2e-2, 2e-1]. A correction sized for the ragged sheet over-corrects the exact one by about the amount it once fixed: 9.829e-2 at κ = 0 on the ragged sheets, item 157 W0 (N). It has not been measured on the exact fixture, so it is never asserted. The separation from (L) is printed.

**Predictions:**
- (L) gap in [1e-7, 1e-4]. The ragged in-run gap was 1.744e-6.
- Column power at 64.5 MHz: Cu in [0.995, 0.9995].
- Rows against item 157's ragged in-run p2 rows (Cu 62.704216 MHz / 18 570 Ω / 984.2; PEC 62.717977 MHz / 32 545 Ω / 1 724.8): \|Δf₀\|/f₀ ≤ 1 %, R_in and Q within ±20 %.
- p2-vs-p1-exact f₀ move ≈ −2 %, the same as on the ragged sheets (−1.999 %).
- k-share ≥ 0.9999; Δ_LOO ≤ 1e-9.

**Windows** (`-n 8`, `FEM_EM_ANS_SKIP_PRIVATE=1`, bare, foreground, Bash `timeout: 660000`). The price is item 132 / 157's degree-2 node windows: 3 nodes in 421 s at `-n 8`, field solve off (`20261001T125923Z_ANS-7-step3d2-W1.log`); one degree-2 call is ≈ 60 s per factorisation (183.9 s for three, `20260921T100022Z_ANS-6-step2-degree2-64MHz.log:942`).
- **W1** (≈ 420 s): nodes `62.0,62.3,63.2`.

  `scripts/testing/run_and_log.sh ANS-7-item196-W1-p2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_ROUTE=B FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_MHZ=62.0,62.3,63.2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **W2** (≈ 380–480 s): nodes `64.5` with `FEM_EM_ANS7_EXACT_HELDKAPPA=1`, labelled `ANS-7-item196-W2-p2`. That is five factorisations: PEC, Cu, P1 and the held-κ Cu + P1.
- **W3** (≈ 450–540 s): nodes `65.8,66.7,67.0` and the assembly, labelled `ANS-7-item196-W3-p2`.
- **Guards:**
  - **Window edge:** if Im Y_in at 62.0 MHz is > 0 in either column, stop before W2 (item 157's guard).
  - **Size:** if W1 runs above 560 s, split W3 into two windows of 1 + 2 nodes.

**Traps.**
- Item 195's traps all apply: the rank-0 `KeyError`, κ only as stated, no `ans:7` stage call, and `logs/ans7_exact/` only.
- The held-κ pass is the control. It never feeds the rows.
- An edit after W1 re-runs W1 (rule (i)).
- A window killed at 590 s is not a result: run `pgrep -c python3`, report, never lengthen the timeout. The node JSONs already written stand.

**Scope.** Degree 2 at h = 0.015, Route B only. Not covered:
- Route A at degree 2 (102 degree-2 drives, `xl`-class);
- h = 0.005 (`xl` entry 24, re-issued on exact sheets by the weekly);
- degree-2 record constants;
- the T13 re-read, which is the weekly's.

**Status it can move.**
- The `ANS-7` row's step-5 sentence gains the degree-2 rows and (L) at κ = 0. No glyph.
- It supplies the order move on exact sheets for the T13 / weekly private re-read.
- It tells the weekly whether `xl` entry 24 should run with κ off on exact sheets.

**Negative result.**
- **(L) red:** the degree-2 sheet is not power-consistent at κ = 0 on an exact rectangle. Print the gap, κ_exact and the column powers; open a known-issues entry; park, BLOCKED. Never apply κ to rescue it.
- **The window-edge guard fires, or a column has ≠ 1 crossing:** BLOCKED.
- **Δ_LOO above `BANDS`:** the rows are printed as unlicensed, and the review re-poses the grid.
- No `PREDICTED RED:` line.

**Outcome (2026-10-03, 07:30 CDT slot): DONE, all green.** Tests only, in item 195's module: knobs `FEM_EM_ANS7_EXACT_MHZ` (per-node split; the whole-grid `route_b` fixture skips when it is set) and `FEM_EM_ANS7_EXACT_HELDKAPPA`; rank 0 writes `logs/ans7_exact/node_<mhz>_p2.json` (cells, sheet tie, S, Z, κ = null, power). *Disclosed:* the held-κ control pass runs with `pec_column=True` (six factorisations in W2, not the brief's five), because item 157's trap forbids `pec_column=False` with `field_solve=True`; the window-edge guard is asserted (`test_p2_window_edge_guard`, Im Y_in(62.0 MHz) < 0 in both columns). No edit after W1's start (module mtime 07:33:47 CDT, W1 start 07:33:56). All windows `-n 8`, degree 2, 80 798 cells, (F) tie green on all four sheets in each window.
- **W1** {62.0, 62.3, 63.2} (`20261003T123356Z_ANS-7-item196-W1-p2.log`, 421 s, 15 passed / 8 skipped): reciprocity ≤ 3.65e-15, σ_max ≤ 0.99992; guard Im Y_in(62.0) −1.1779e-3 (Cu) / −1.2010e-3 S (PEC); column power Cu 0.99874, PEC 0.99935. Below 560 s, so W3 was not split.
- **W2** {64.5} + held κ (`20261003T124111Z_ANS-7-item196-W2-p2.log`, 409 s, 15 passed): **(L) κ off: L_S 1.295634668e-3, L_F 1.282896503e-3, \|L_S − L_F\| = 1.273817e-5 ≤ `REDUCTION_BAND` 1e-3** (predicted [1e-7, 1e-4], met; the ragged in-run gap was 1.744e-6). Cu column power 0.998704 (predicted [0.995, 0.9995], met), PEC 0.999311; copper identity residual 1.355e-12 ≤ 1e-6. **Control (printed, never asserted):** κ = 5.2094420914264195e-02 on the exact sheets gives Cu column power 1.10536 (> 1: S-derived loss −0.1054), \|L_S − L_F\| = 1.067250e-1 (predicted [2e-2, 2e-1], met), 8.4e3 × the κ-off gap.
- **W3** {65.8, 66.7, 67.0} + assembly (`20261003T124814Z_ANS-7-item196-W3-p2.log`, 405 s, 15 passed): one − → + zero per column; k-share 0.9999950 (Cu) / 0.9999973 (PEC), predicted ≥ 0.9999, met; **Δ_LOO f₀ / R_in / Q = 6.30e-13 / 1.64e-9 / 1.55e-9 ≤ `BANDS` 1e-6** (predicted ≤ 1e-9: not met on R_in / Q by ≈ 1.6×, printed only — the grid is licensed). Rows (printed, no constant): Cu f₀ 62.701081939 MHz, R_in 18 441.755 Ω, Q 977.261; PEC f₀ 62.714806597 MHz, R_in 32 123.933 Ω, Q 1 701.733. Against item 157's ragged in-run p2 rows: f₀ −0.0050 % / −0.0051 %, R_in −0.69 % / −1.29 %, Q −0.71 % / −1.34 % (predictions ≤ 1 % / ±20 %, met). Order move p2 − p1 on exact sheets: f₀ −1.985 % (both columns; predicted ≈ −2 %, the ragged −1.999 %), R_in +9.07 % / +14.37 %, Q +7.07 % / +12.22 %.
- **Reading:** on exact rectangles the degree-2 sheets are power-consistent with no width correction, and the degree-2 rows reproduce the ragged κ-in-run rows to ≤ 1.3 % — item 157's κ-confound reading holds; the correction sized for the ragged sheet over-corrects the exact one. For the weekly: `xl` entry 24 can run κ off on exact sheets.

## §7 row moved verbatim 2026-10-01 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/ANS-7.md`)* **The tuned copper F-small birdcage: a capacitor on every gap, the drive across one of them, `Z_in(f)` across the k = 1 resonance — the first AED case with lumped capacitors** (**SPEC v2, operator directive 2026-09-28**: C = 66.29 pF on all four gaps, drive across C₁, parallel resonance, Route A in-model + Route B four-port/circuit, rows on `Y_in` over 62–67 MHz; v1 — 15.57 pF on three gaps, bare feed, series resonance — withdrawn before any AED result; the v1 text and item 50's measurement below are the v1 record; the runnable half is §10 T10 as re-targeted) — *v1 as commissioned:* (commissioned 2026-09-26 by the weekly on `PORT-22` ✅; §10 chain T10; SPEC `examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/SPEC.md`; **serial on `PORT-22` step 2**, chain T7) | 🟡 **Step 3 complete — (ii), the degree-2 leg, landed 2026-10-01 (07:30 CDT slot, §9 item 132, chain T11) — a record, no glyph (T13 decides).** Rule (a) lift: `_build_ladder(pec_column=True, field_solve=True)` keywords, defaults byte-identical, re-gated by W0 (`TH-14` module 23 passed, 133 s, `20261001T125435Z_ANS-7-step3d2-W0.log:1704`). Knobs `FEM_EM_ANS7_SPARSE_FIELD_SOLVE` / `_COLUMNS` added. **Wiring control (asserted):** degree 1, Cu only, no PEC column, no P1 solve, 64.5 MHz reproduces 127's stored Cu S to 1.136e-14 ≤ 1e-9 (`…125801Z_ANS-7-step3d2-W0b.log:997`, 74 s, `-n 4`). Degree 2, κ held at the degree-1 value 1.0599570850e-02, `-n 8`, 80 181 cells, field solve off: W1 {62.0, 62.3, 63.2} 421 s, W2 {64.5, 65.8, 66.7} 426 s, W3 {67.0} + assembly 189 s, all Status 0 (`…125923Z_…-W1.log`, `…130632Z_…-W2.log`, `…131348Z_…-W3.log`); per-node reciprocity ≤ 4.44e-14 (band 1e-3). **Assembly (asserted, W3 `:1185–1195`):** one − → + zero per column; k = ±1 share 0.98664 (Cu) / 0.98687 (PEC) > 0.5; **degree-2 Δ_LOO f₀ / R_in / Q = 5.71e-13 / 1.34e-11 / 8.21e-11 ≤ (i)'s bands 1e-6 — the 7-node grid is licensed at degree 2** (predicted pass, met). **Printed rows at degree 2:** Cu f₀ 62.657821 MHz, R_in 355.18 Ω, Q 18.137; PEC f₀ 62.672360 MHz, R_in 357.86 Ω, Q 18.286. **Order move (p2 − p1)/p1:** f₀ −2.071 % (Cu) / −2.070 % (PEC); R_in −97.9 % / −98.7 %; Q −98.0 % / −98.8 %. The f₀ move is the size of the 2026-09-30 known-issues gap; the ≈ 50× fall in R_in and Q at one order step (with κ held at its degree-1 value, and the k = ±1 share down from 0.99999 to 0.987) is **not interpreted here** — a record for T11b / T13, not a ruling. **Step 3 (i) landed 2026-10-01 (04:30 CDT slot, §9 item 127, chain T11) — a record, no glyph; tests only.** New `tests/validation/test_ans7_sparse_route_b.py` (knobs `FEM_EM_ANS7_SPARSE_MHZ` / `_DEGREE` / `_RESOLUTION`): the 7-node degree-1 sparse Route B (62.0 … 67.0 MHz, barycentric degree 6 on `z_from_s`) reproduces the tracked 51-point rows — Cu f₀ / R_in / Q relative **9.7e-14 / 3.22e-10 / 2.06e-10**, PEC **4.0e-14 / 1.04e-11 / 2.08e-11** — inside bands **B = 1e-6** frozen from W1 (every 3·Δ_LOO under the floor: Δ_LOO ≤ 2.5e-14 / 1.94e-11 / 1.98e-11); node S vs stored `B1`/`B2` ≤ 2.27e-13 (band 1e-9); one − → + zero per column; k = ±1 share 0.99999385 / 0.99999793. **Negative control (3 nodes) misses as asserted but far closer than predicted:** R_in 1.9e-6 / 2.7e-6, Q 1.3e-5 / 1.4e-5 (predicted ≈ 1e-3–1e-2 on R_in) — Z(f) is smoother in the window than the prediction assumed; a record for the weekly. Printed: off-node max entrywise ΔZ 2.63e-11, 51-point max \|ΔZ_in\|/\|Z_in\| 2.86e-10; **the 5-node subset was not printed** (its 62.7 / 66.3 MHz nodes are not solved nodes) — open printed figure. Windows `20261001T095511Z_ANS-7-step3-W1.log` (157 s, Status 0) and `…095802Z_ANS-7-step3-W2.log` (157 s, Status 0, 15 passed), `-n 4`, 80 181 cells. `xl` entry 24's grid has a measured degree-1 licence; item 132 (degree 2) unblocks. **Adjudicated 2026-09-30 (weekly, §10 rulings; known-issues 2026-09-30): AGREE on the mode (k = ±1 share ≥ 0.9999 in both codes), on Route A = Route B and on the copper contribution (the conductor's 1/Q within 5 %, `Q_Cu/Q_PEC` within 1.1–2.4 %); `f₀` 1.2–1.4 % above the order-matched AED columns and 1.9 % above First Order, Q 6–10 % and `R_in` 7–11 % low — the absolute rows are INCONCLUSIVE at degree 1 on one rung. The row stays 🟡 until §10 T11 / T12 supply the order- and h-matched rung and the weekly re-reads (T13).** **Runnable half done 2026-09-29 (09:30 slot, §9 item 68: Route A, the band, the private writer) — the v2 Done-when is met; the row stays 🟡 pending adjudication (T10b, the weekly's).** Stage `A` of `07_…py`, `-n 4`, `FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1`, 367 s (`20260929T100040Z_ANS-7-item68-W3-A.log`, Status 0; 80 181 cells = tracked): κ recomputed in-run 1.059957085005e-02, \|Δκ\| 4.4e-16 ≤ 1e-9 (asserted); **band (asserted, imported `REDUCTION_BAND` 1e-3): all 51 points × 2 columns, `S11_bare` Route A vs Route B worst 9.250426e-05 (Cu, 67.0 MHz; PEC worst 9.203e-05), min/median 3.38e-05 / 4.97e-05 Cu, 3.38e-05 / 4.96e-05 PEC** — the Cu residual grows monotonically in f across the window (3.4e-05 → 9.3e-05), straddling the predicted 6–8e-5 (printed, not a gate). Exactly one − → + zero of `Im Y_in` per column on Route A (asserted): Cu 63.982289 MHz, PEC 63.996361 MHz. Printed rows (Route A): Cu R_in 1.700853e+04 Ω, Q 919.0374; PEC R_in 2.845367e+04 Ω, Q 1537.6010; A − B: Cu `df0_rel` −1.161e-05, `dq_rel` +6.494e-04, `dr_in_rel` +5.777e-04, `max_dzin_rel` +1.905e-02; PEC −1.160e-05 / +1.074e-03 / +1.003e-03 / +3.516e-02 (the Z_in near-pole amplifies the ≤ 1e-4 S residual; not gated, per the item); (d) Cu vs PEC Δf₀/f₀ −2.198887e-04. **Controls (predicted, rule (e)):** (i) circuit at 0.5·C vs in-model at C, 64.0 MHz Cu: 4.836e-01 = 483.6× the band (predicted ≥ 10×; met); (ii) uncorrected Route A Cu at 64 MHz vs the circuit on item 69's uncorrected 4×4: 1.657e-02 = 16.6× the band (predicted ≈ 1e-2; the band can fail). Route A cells + band + controls written into `metrics.json`, Route B entries unchanged; `COMPARISON.md` untouched (byte-identical); combined XDMF (P1 vector, Cu + PEC E at 64.0 MHz) under the gitignored `/workspace/logs/ans7/`. Private writer (`_write_private_if_available`, keys per `aed/ans7_hfss_pyaedt.py:44–61`, names only): `tests/unit/test_ans7_private_comparison.py` — sentinels in the AED row, tracked files byte-identical with no sentinel, knob ⇒ nothing, no JSON ⇒ nothing; negative control: the leaking variant is caught (all four sentinels in the tracked copy) — 4 passed with the pin test's 20 (`20260929T100722Z_ANS-7-item68-W4-unit.log`, 11 s). Censuses after: docrefs dead 0 / guide 0 / stale 0, setup figures 58 / 57 / 1 / 0 — identical to item 69's after-censuses (`20260929T094749Z`/`094751Z`), which were this item's before-state. *Disclosed:* the first W3 (`20260929T095248Z`) never reached the container (docker socket permission denied, 0 s; that call had been piped through `grep`, the bare re-issue ran); the second (`20260929T095257Z`, 443 s) went green on every asserted anchor across all 102 solves then died in the XDMF writer on a DG1 space (dof layout ≠ mesh), before metrics were written — fixed to P1, the third window is the record. **Route B landed 2026-09-29 (09:30 slot, §9 item 69).** Re-run as committed at `-n 4`: W0 `TH-14` module 23 passed, 121 s (`20260929T093153Z_ANS-7-item69-W0-TH-14.log:1505`); W1 `B1` 353 s — anchor 1 worst rel 5.042e-14 ≤ `LEG_D0_REPRODUCTION_BAND` 1e-9, negative control (asserted) 4.189583e-04 / 1.267300e-04 / 1.847852e-04 > 1e-6, κ 1.059957085e-02 (`20260929T093403Z_ANS-7-item69-W1-B1.log:941–943`); W2 `B2` 331 s — Cu `f₀` 63.983031 MHz, R_in 1.699871e+04 Ω, Q 918.4410, k-share 0.999994; PEC `f₀` 63.997103 MHz, R_in 2.842515e+04 Ω, Q 1535.9516, k-share 0.999998; exactly one − → + zero per column, 2·C control met; S-vs-Z route record 7.891e-13 (`20260929T094009Z_ANS-7-item69-W2-B2.log:2819–2830`). Every predicted reading reproduced to the printed digits; regenerated `metrics.json` vs the attempt-branch copy: worst relative 5.482e-11 excluding timings and sub-1e-12 residuals (`20260929T094622Z_ANS-7-item69-metrics-vs-branch.log:38`), `COMPARISON.md` byte-identical. `metrics.json` declared in `LISTED_EXAMPLE_ARTIFACTS`. **Parked 2026-09-28 21:00 (§9 item 63), ruled by the 2026-09-29 03:00 review — land it (§9 item 69):** every asserted anchor green on the fixed code (`20260929T022842Z_ANS-7-W2-B2.log:2819–2822`: Cu `f₀` 63.983031 MHz, PEC 63.997103 MHz, one − → + zero each, k = ±1 share 0.999994 / 0.999998; anchor 1 2.04e-15; κ 1.059957e-02). **Anchor (3) is a record, not a gate (rule (h)):** the first W1's 6–307× miss (`…021845Z_ANS-7-W1-B1.log:3624–3648`) fed the sweep's `z_matrix` — on the lumped-sheet route the *terminated* transimpedance, documented as a diagnostic (`ports/sparameters.py:347–356`, `:484–493`) — to the open-circuit Z formula (the item's drafting error); the re-routed check compares two algebraic routes from one S (5.334e-13, W2 `:2818`) and gates the circuit layer's algebra only. **Printed Q 918.44 (Cu) / 1535.95 (PEC) vs "≈ 68":** the prediction was imported from v1's series resonance (a different circuit) — not a failure. The v1 → v2 phantom:conductor loss ratio falls ≈ 24× (36 → 1.49) against (66.29 / 15.57)² ≈ 18 — a hypothesis (gap-voltage E-field loss in the phantom), for T10b, not a ruling. **Measured 2026-09-27 (§9 item 50, a design measurement — no gate, band or record):** fixed C = 15.57 pF on P2–P4, P1 driven, degree 1, F-small hole mesh (80 181 cells), `_build_ladder(frequencies_hz=…)` + `PORT-15`'s `tuned_input`, `scripts/probes/port22_step2_anchor_q.py`, `-n 4`. **Window A** (SPEC grid, 44 → 84 MHz): one sign change of `Im Z_in` per column, a zero (− → +) in [66, 68] MHz, no pole; Cu f₀ 66.298 MHz, R_in 0.7197 Ω, Q 61.93, width 1.070 MHz; PEC f₀ 66.313 MHz, R_in 0.7003 Ω, Q 63.66, width 1.042 MHz (`20260927T093227Z_ANS-7-item50-windowA.log:2445–2448`, 340 s). **Window B** (65.3 → 67.3 MHz, 0.1 MHz): Cu f₀ 66.2951 MHz, R_in 0.7197 Ω, Q 62.52, width 1.060 MHz (half-power 65.766–66.827 MHz); PEC f₀ 66.3095 MHz, R_in 0.7003 Ω, Q 64.18, width 1.033 MHz (`20260927T093815Z_ANS-7-item50-windowB.log:2444–2447`, 356 s). **SPEC 2 MHz grid points inside f₀ ± f₀/(2Q): 1 per column** (66 MHz); 11 of window B's 21. So SPEC box 1's "does the resonance sit inside the grid": yes, at ≈ 66.3 MHz; T7 (c)'s "≥ 5 grid points in the half-power band" is unreachable on the 2 MHz grid (it needs a step ≲ 0.2 MHz) — confirmed, though the estimate behind it was wrong: R_in is ≈ 0.72 Ω, not ≈ 0.06 Ω, because the phantom dominates the loss (Cu − PEC ΔR_in ≈ 0.019 Ω is the conductor's share), so Q ≈ 62–64, not hundreds. Cu − PEC: Δf₀ −14 kHz, ΔQ −1.7. Self-checks green: 64 MHz S vs tracked `ans:6` `metrics.json` 3.817e-15 (≤ 1e-9); reciprocity ≤ 5.78e-15; identity residual ≤ 7.18e-13; negative control \|S(Cu) − S(PEC)\| at 64 MHz 4.19e-4 / 1.27e-4 / 1.85e-4 (> 1e-6). Sheet widths uncorrected (κ ≈ 1 % moves f₀ by ≈ κ/2). The weekly re-poses T7 (c) and, if it wants (c) on the SPEC, the grid. **Runnable half, v2 (amended 2026-09-28 18:00 review under standing rule (j) to §10 T10 as re-targeted by the operator the same day; queued as §9 items 63 + 68):** `ans:7` (`07_tuned_copper_birdcage_driven_sweep.py` + same-stem guide), C = 66.29 pF on every gap, drive across C₁, degree 1, Cu and PEC columns, the SPEC's 62–67 MHz / 0.1 MHz FineWindow; κ-corrected sheet widths on both routes (uncorrected sheets miss `REDUCTION_BAND` by 13.5× at 64 MHz, `docs/planning/chunks/PORT-14.md:650`). Route B = `_build_ladder` four-port through `PORT-15`, C₁ across the drive; Route A = in-model (bare port + `PORT-14` sheets on gaps 2–4, C₁ added across the port). Rows (a)–(d) on `Y_in` and Route B's k = ±1 gap-current share at `f₀`, printed; `metrics.json` + `COMPARISON.md` with AED columns blank; private-mode writer on `ans:6`'s pattern, `FEM_EM_ANS_SKIP_PRIVATE` honoured. **Asserted (imported, never restated):** Route A = Route B on `S11_bare` within `REDUCTION_BAND` at every window frequency, both columns; exactly one − → + zero of `Im Y_in` per column in the window; reciprocity and the copper surface-loss identity at every point. **Done-when (§4, v2):** those asserts green through the harness, elapsed recorded, censuses before and after ⇒ the runnable half is done and the row goes 🟡 pending adjudication (T10b, the weekly's). **Not claimed:** any absolute `Z_in`, a match. *v1 runnable half (superseded, kept as the record):* `ans:7`, an example + same-stem guide dispatched through `./run_examples.sh -e ans:7`, importing every constant and the sweep from `PORT-22` step 2's gate module (rule (a)); Cu and PEC columns, degree 1, the SPEC's frequency grid; `metrics.json` and `COMPARISON.md` with AED columns blank; combined XDMF at the grid point nearest `f₀`. **Asserted (imported, never restated):** circuit = field within `REDUCTION_BAND` at every grid frequency; exactly one sign change of `Im Z_in`; the copper surface-loss identity. *v1 Done-when (superseded):* the example runs green through the harness with those asserts, elapsed recorded, censuses before and after. **Not claimed:** any absolute `Z_in`, a match, a mode. Adjudication is the weekly's, at a stated relative level; no AED value enters a tracked file. **Step 3 (§10 T11; written 2026-10-01 03:00 review, rule (j)):** (i) the 7-node sparse Route B (62.0 … 67.0 MHz, entrywise barycentric interpolation of the S-derived network Z, terminated on the 51-point grid) reproduces the dense f₀ / Q / `R_in` rows of both columns inside bands frozen from its own leave-one-out residual (`TH-18` step 2's rule), the 3-node control misses (§9 item 127); (ii) the same at degree 2 on the same mesh, κ held at the degree-1 value, rows and the order move printed, the grid licensed if its degree-2 Δ_LOO sits inside (i)'s bands (§9 item 132); elapsed recorded ⇒ step 3 done, no glyph (T13 decides). (ii) supplies `xl` entry 24's knob names. **(iii) (§10 T11b; written 2026-10-01 18:00 review, rule (j); §9 item 157):** the degree-2 leg of (ii), re-run with the in-run degree-2 κ (`FEM_EM_ANS7_SPARSE_KAPPA=in_run` / `stored`; `_build_ladder(measure_kappa=True)` on an uncorrected Cu P1 solve at 64.0 MHz, the route `ans:7` used for its degree-1 κ). At degree 1 the route reproduces the tracked κ within `KAPPA_RECOMPUTE_ATOL`. At degree 2 and 64.5 MHz, the P1 column's S-derived loss 1 − Σ\|S_i1\|² equals the field-balance loss within `REDUCTION_BAND`, and the uncorrected pass misses it by more than 10×. The 7-node grid at that κ asserts (ii)'s per-node reciprocity, one − → + zero per column, k-share > `K_SHARE_FLOOR`, and Δ_LOO inside (i)'s bands. Rows are printed beside (ii)'s held-κ rows and the degree-1 rows, with the pre-registered reading (κ confound / genuine order move / partial / not κ / mixed) recorded. Elapsed recorded ⇒ step 3 (iii) done, no glyph (T13 decides), and `xl` entry 24's prerequisite is met; the next review issues its command. **Step 3 (ii)'s order move read by the 2026-10-01 10:30 review (`log-pathologist`, COUNTABLE):** the degree-2 R_in / Q fall (≈ −98 %) is **not established as an order move** — the raw degree-2 4-port S already carries column power ≈ 0.919 (≈ 8 % unaccounted, `20261001T125923Z_ANS-7-step3d2-W1.log:1077,1106`) against 0.9986–0.9992 at degree 1 (`20261001T095511Z_ANS-7-step3-W1.log:1211,1240`), the size a held degree-1 κ predicts (degree-2 κ ≈ 5.2e-2 on `PORT-21`'s sibling fixture; on this mesh a 1 % κ error alone made ≈ 2.2 pp of S-derived loss absent from the field balance, `…125435Z_ANS-7-step3d2-W0.log:1134–1221,1251`); the copper part of 1/Q survives the order step (4.38e-4 → 4.50e-4). The f₀ −2.07 % is a degree + κ composite. **`xl` entry 24 is held until T11b supplies the in-run degree-2 κ** (T11b serial on F13, §9 item 140). 

## §7 row moved verbatim 2026-10-03 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/ANS-7.md`)* The tuned copper F-small birdcage: a capacitor on every gap, the drive across one of them, `Z_in(f)` … | 🟡 Adjudicated 2026-09-30: AGREE on the mode, on Route A against Route B and on the copper contribution, but the absolute rows are INCONCLUSIVE at degree 1 on one rung. Step 3 completed 2026-10-01 as records (the 7-node sparse Route B licensed at degree 1 and degree 2, plus the degree-2 rows), and the row stays open until T11/T12 supply the order- and h-matched rung for the weekly's T13 re-read. Steps 5 / 5b (T17 / T18, exact port and capacitor rectangles, κ off, records) are queued as §9 items 195 / 196 (daily provisional, 2026-10-03). **Step 5 record (T17, §9 item 195, 2026-10-03):** exact port and capacitor rectangles on the 80 798-cell hole mesh, degree 1, κ off — Route A = Route B to ≤ 3.11e-5 at all 102 points (band `REDUCTION_BAND`), the ragged κ-off control misses by 16.57 × band; Route B rows registered under (1*), f₀ −0.019 % against the ragged κ-corrected rows in both columns (`20261003T113037Z_ANS-7-item195-W2-routeA.log`, `20261003T113959Z_ANS-7-item195-W3-routeB.log`). A record, no glyph. **Step 5b record (T18, §9 item 196, 2026-10-03):** the same exact rectangles at degree 2, κ off, 7-node Route B at `-n 8` — the loss identity (L) at 64.5 MHz closes to 1.274e-5 ≤ `REDUCTION_BAND` with no correction (the ragged degree-2 κ applied to the exact sheets misses by 1.067e-1, printed control); one − → + zero per column, k-share ≥ 0.999995, Δ_LOO ≤ 1.64e-9 inside `BANDS`; degree-2 rows Cu f₀ 62.701082 MHz / R_in 18 441.75 Ω / Q 977.26, PEC 62.714807 MHz / 32 123.93 Ω / 1 701.73 (−0.005 % f₀ against item 157's ragged in-run degree-2 rows; order move p2 − p1 f₀ −1.985 %, R_in +9.1 % / +14.4 %, Q +7.1 % / +12.2 %), printed, no constant registered (`20261003T124111Z_ANS-7-item196-W2-p2.log`, `20261003T124814Z_ANS-7-item196-W3-p2.log`). A record, no glyph. **Done-when:** The row stays 🟡 until §10 T11 / T12 supply the order- and h-matched rung and the weekly re-reads (T13). *History: `docs/planning/chunks/ANS-7.md`.*

### §9 item 214 — 2026-10-04 (chain T19)

*Drafted by a general-purpose agent for the 2026-10-04 03:00 daily review; the review reviews it.* Licence: §10 chain T19 (weekly 2026-10-03), whose target is quoted, not changed: "`FEM_EM_ANS7_EXACT_RESOLUTION` on the exact-strip module, unset digit-identical to item 196's rows; one intermediate rung, sized by a one-node cost probe so the five-node Cu grid fits one heavy window, at degree 2: the h-move of `f₀`, Q, `R_in` printed. Target: T11's identity on each rung." It is also the prerequisite contract of `xl` entry 28 (`docs/testing/xl-pending.md`, § 28): the knob **and** "a Cu-only five-node path". The 2026-10-02 sheets-and-κ ruling binds: exact rectangles, κ off. No band, tolerance, Done-when or validation target moves.

`ANS-7` step 6: a resolution knob and a Cu-only five-node path on the exact-strip module; the ×1 rows reproduced through them; one intermediate rung at degree 2.
- **Executor:** implementer. Tests only, plus one additive keyword on a test helper (disclosed).
- **Tier:** heavy; `-n 8` (W0), `-n 1` (W0b), `-n 12` (W1–W3); complex build.
- **Class:** chain step (solve). **Time:** ≈ 18 measured / 60 predicted min.
- **Depends on nothing.** It shares `tests/validation/test_ans7_exact_strip.py` with item 217. Whichever lands second rebases onto the first and runs its own windows on the module as committed (rule (i)).

**What "T11's identity" means here.** T11 (step 3 (i)/(ii), items 127 / 132) licensed the sparse Route B with five checks, and those are the anchor on every rung:
- per-node reciprocity ≤ `RECIPROCITY_BAND` and σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`;
- exactly one − → + zero of Im `Y_in` per column;
- k = ±1 share > `K_SHARE_FLOOR`;
- Δ_LOO on f₀ / R_in / Q ≤ `BANDS` (1e-6, frozen from item 127's W1, never moved).

"Each rung" means two rungs. The ×1 rung (h = 0.015) is checked on item 196's stored seven degree-2 JSONs for the Cu-only five-node subset, a pure test. The intermediate rung is solved here.

**Code** (`tests/validation/test_ans7_exact_strip.py`, plus one additive keyword in `tests/validation/test_ans6_exact_strip.py`; no `src/`, no `TH-14` edit, no `ans:7` edit):
1. **Knobs.** These two names, plus the existing `_ROUTE`, `_DEGREE` and `_MHZ`, are entry 28's contract:
   - `FEM_EM_ANS7_EXACT_RESOLUTION`: a float in metres. Unset ⇒ `None` ⇒ nothing is forwarded.
   - `FEM_EM_ANS7_EXACT_COLUMNS`: unset ⇒ `("cu", "pec")`; `cu` ⇒ Cu only. Any other value raises `ValueError`.
2. **Forwarding.**
   - `route_b` and `p2_nodes` call `_build_ladder(..., resolution=res)`. With `None` this is today's call.
   - `_fixture_tie(res)` calls `_hole_rung(GATED_WIDTH_FRACTION, resolution=res)`.
   - In `_hole_rung`, add a `resolution=None` keyword and pass it to `_sheets_build` only when it is not `None`. The default stays byte-identical. Rule (a): the `ANS-6` exact module's mesh tests are re-run green in W0b.
3. **The Cu-only path** (in `p2_nodes`):
   - Set `pec_column = "pec" in cols` and `field = (LOSS_NODE_MHZ in nodes) and pec_column`. This keeps out item 157's rank-0 `KeyError` trap (`pec_column=False` with `field_solve=True`).
   - `HELDKAPPA=1` together with `COLUMNS=cu` raises `ValueError`. So does the whole-grid `route_b` fixture with `COLUMNS=cu`.
   - Every per-node loop and test iterates over `cols`, never a hard-coded `("cu", "pec")`. That covers node anchors, the window-edge guard, the Im `Y_in` print and the JSON write.
4. **JSON names.**
   - Unset resolution with both columns keeps `node_<mhz>_p<deg>.json`, item 196's names.
   - Otherwise: `node_<mhz>_p<deg>[_h<res:g>][_cu].json`.
   - New keys `"resolution"` and `"columns"` are added. When they are absent, as on item 196's files, they read as `None` and both columns.
   - Rank 0 writes, then a barrier.
5. **Assembly** (`test_p2_assembly`, generalised):
   - **Required node set:** `DEFAULT_NODES_MHZ` for both columns. For `cu` it is `ASSEMBLY5_NODES_MHZ = (62.0, 63.2, 64.5, 65.8, 67.0)`, imported from `test_ans7_sparse_route_b` (entry 24's set, licensed at the ragged mesh by item 173).
   - It reads only JSONs at this run's (degree, resolution, columns). It asserts one degree, one cell count, κ `None` and one resolution.
   - `_loo` takes `columns=("cu", "pec")` as its default, so the unset path is byte-identical. `sparse_rows` is called with `columns=cols`.
   - **Asserts:** one − → + zero per column, k-share > `K_SHARE_FLOOR`, and Δ_LOO ≤ `BANDS`. With five nodes there are three interior drops, so the interpolants have four nodes.
   - **Prints:**
     - the rows;
     - B₅ = ceil_1sf(max(3·Δ_LOO, 1e-6));
     - Im `Y_in` at every node;
     - the 3-node control (`CONTROL_NODES_MHZ`) against the set;
     - the per-node solve time;
     - summed `ru_maxrss` (the `06_…py:300` idiom);
     - when `res` is set, the h-move (rung − ×1)/×1 of f₀ / R_in / Q, against `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2["cu"]` and against the ×1 five-node Cu rows from the stored JSONs.
6. **(F) at a resolution.**
   - Unset: today's (F), unchanged.
   - Set:
     - all four N_i = `SPEC_AREA_M2` and w = `SPEC_WIDTH_M` within `STRIP_BAND` (CAD rectangles, independent of h; asserted);
     - the tie's cells equal the ladder's cells exactly (asserted: the tie measured the mesh the ladder solved);
     - cells differ from `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` by more than `CELL_COUNT_BAND` (asserted: (N-wire) below).
7. **Record constant, registered before W0.** `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2` comes from item 196 W3 (`20261003T124814Z_ANS-7-item196-W3-p2.log:2189,2193`, full repr):
   - Cu: f₀ 62701081.93936459 Hz, R_in 18441.754600739256 Ω, Q 977.2613699232847;
   - PEC: 62714806.5966287 Hz, 32123.932870177392 Ω, 1701.733130785279.

   It carries a version-tag comment (complex 0.11 image, `-n 8`) and the (1*) note. If item 217 landed first, reuse its constant of the same name.

   New `test_r2_p2_rows_record`: at degree 2, resolution unset, both columns, once the assembly has run, rows vs record ≤ `NODE_S_RTOL` (1e-9; item 173's solved-vs-JSON precedent measured 3.4e-11). This is the "unset digit-identical" assertion. The relative differences are printed.
8. **The Cu-only five-node licence at ×1.** New pure `test_p2_cu5_subset_anchor`, which needs `ROUTE=B` and runs no solve. When the seven ×1 two-column p2 JSONs exist (else skip):
   - the five-node Cu rows vs the seven-node Cu rows ≤ `BANDS` (asserted);
   - the five-node Δ_LOO ≤ `BANDS` (asserted), with B₅ printed;
   - the 3-node control, printed.

**Anchors** (asserted; every band imported, none restated):
- (F);
- per node (B);
- the assembly's identity on the rung;
- (R2) at ×1;
- (S5), the ×1 Cu-5 licence.

**Negative controls (rule (e)).**
- **(N-wire): asserted, structural, no factor.** Rung cells ≠ the ×1 record beyond `CELL_COUNT_BAND`. A knob that never reaches the mesh reproduces 80 798 cells and item 196's rows. The same generator path moved 80 181 → 226 642 cells at h = 0.0075 on the ragged sibling (`20260923T093956Z_ANS-6-step2b-b-h0p0075-d1-64MHz.log:950`).
- **(N3): predicted, printed.** The 3-node Cu control misses B on R_in or Q at ×1 (exact p2) and on the rung.
  - It is backed only on other sets: ragged p2_kin 3.474e-6 / 2.240e-6 (`20261002T111756Z_ANS-7-item173-W1.log`), so the thin side is 2.2×.
  - It has never been measured on exact sheets, so it is not asserted.
  - It is an interpolation error, not a saturating fraction, so no ceiling applies.

**Predictions** (written first, printed):
- **(R2):** relative ≤ 1e-11. W0 runs at item 196's `-n 8`; one node is re-solved and six are read.
- **(S5) at ×1:** \|row₅ − row₇\|/\|row₇\| ≈ 1e-12 / 1e-10 / 3e-10 (item 173's ragged p2_kin: 1.233e-12 / 1.267e-10 / 2.780e-10). Five-node Δ_LOO ≈ 1e-9 / 4e-8 / 3e-8 (item 173 W3: 6.5e-10 / 3.6e-8 / 2.8e-8).
- **Rung h = 0.01 price** (scaled from the measured ×1 and h = 0.005 points, below):
  - cells in [135 000, 175 000];
  - degree-2 unknowns ≈ 6.6 × cells;
  - Cu sweep per node at `-n 12` in [95, 175] s;
  - mesh 30–40 s;
  - summed `ru_maxrss` in [30, 60] GiB.
- **h-move against ×1** (degree 2, Cu): f₀ in [−0.5 %, +0.5 %]; R_in and Q in [−10 %, +10 %]; k-share ≥ 0.9999; Im `Y_in`(62.0 MHz) < 0; one crossing. Rung Δ_LOO(5) ≤ 1e-7.

**Price, measured.**
- **×1, degree 2, `-n 8`, 80 798 cells:**
  - PEC sweep 52.9–62.2 s and Cu sweep 44.3–50.0 s per node;
  - mesh 22–23 s;
  - fixed cost ≈ 78 s per window (418.7 s total − 340.6 s solve, `20261003T123356Z_ANS-7-item196-W1-p2.log:2025–2156,2255`);
  - summed `ru_maxrss` 13.78 GiB for three degree-2 factorisations (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:942`).
- **h = 0.005, degree 2, `-n 16`:** 523 816 cells, 3 466 952 unknowns, Cu sweep 1 302.3 s, memory.peak 286.5 GB (`20260929T070008Z_ANS-6-step2c.log:961,1026,1037,1045`).
- **Ragged hole mesh at h = 0.0075:** 226 642 cells, mesh 46.9 s.
- **Scaling between the two degree-2 points:** wall ∝ N^≈1.9, memory ∝ N^≈1.6.
- **Consequence for h = 0.0075:** unknowns ≈ 1.5 M; a Cu node ≈ 300–400 s at `-n 12` and ≈ 75 GiB. Five Cu nodes come to ≈ 1 700–2 100 s, **over** the 1 200 s heavy ceiling. So 0.0075 fails the sizing criterion on prediction, and the probe starts at 0.01.

**Windows** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1`; every window ≤ 590 s, so §5.1 durable capture is not needed). **The five-node Cu grid does not fit one ≤ 590 s foreground window** at the predicted price: 5 × (95–175) s plus ≈ 150 s fixed exceeds 590 s. "Fits one heavy window" is therefore read against the 1 200 s tier ceiling, and the nodes are run in split windows.

Before W0, on the host, back up the stored set: `cp -a logs/ans7_exact logs/ans7_exact.item196`. The path is gitignored. W0 overwrites `node_62.0_p2.json` with its reproduction; on a red, restore it from the copy.
- **W0, the ×1 reproduction** (`-n 8`, knob unset, ≈ 200–260 s): node 62.0 MHz, both columns. The assembly then runs on the seven JSONs, and (R2), (S5) and the window-edge guard are checked.

  `scripts/testing/run_and_log.sh ANS-7-item214-W0-x1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_ROUTE=B FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_MHZ=62.0 timeout -k 30 400 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **W0b, rule (a) for `_hole_rung`** (≈ 100 s): item 193 W1's command verbatim (`20261003T095457Z_ANS-6-item193-W1-mesh.log:12`; 19 passed, 15 skipped, 97 s), labelled `ANS-7-item214-W0b-ans6mesh`. The pass set must be the same.
- **W1, the cost probe** (`-n 12`, h = 0.01, Cu, node 62.0, ≈ 200–330 s):

  `scripts/testing/run_and_log.sh ANS-7-item214-W1-probe "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_ROUTE=B FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_MHZ=62.0 timeout -k 30 590 mpiexec -n 12 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`

  **Rung rule.** t is the Cu sweep wall time at 62.0 MHz, from `TH-14`'s "sweep … s" line. The chosen h is written into the outcome before W2.
  - t ≤ 85 s: 0.0075 is predicted to fit (≈ 2.35× per node). Run **W1b**, W1 at `RESOLUTION=0.0075`, labelled `…-W1b-probe`. Take 0.0075 if its t ≤ 200 s; otherwise keep 0.01.
  - 85 s < t ≤ 200 s: **h = 0.01 is the rung.** 5 × 200 s + 2 meshes + ≈ 80 s fixed ≤ 1 200 s, so the five-node Cu grid fits one heavy window.
  - t > 200 s: run **W1b** at 0.0125. If that is still > 200 s, stop: BLOCKED, with both prices.
  - Summed `ru_maxrss` > 100 GiB at any probe: stop, BLOCKED (128 G container).
  - The window-edge guard is red (Im `Y_in`(62.0) ≥ 0): stop. The rung's f₀ has left the window below 62 MHz, and that bound is the finding (the review re-poses the grid; entry 28 stays PENDING).
  - The probe's node JSON is the rung's 62.0 MHz node. It is never re-solved.
- **W2** (≈ 2t + 150 s): W1's command at the chosen h with `FEM_EM_ANS7_EXACT_MHZ=63.2,64.5`, labelled `ANS-7-item214-W2-rung`.
- **W3** (≈ 2t + 160 s): `MHZ=65.8,67.0`, labelled `ANS-7-item214-W3-rung`. The assembly runs here.
- **If the chosen h has t > 210 s** (possible only on the 0.0075 branch), run one node per window (W2a/W2b/W3a/W3b), still with no edit.
- **No edit to the module after W0** (rule (i)). The rung is chosen by the knob value only.
- **If the slot reaches minute 50 before W3:** park per rule (d) — the node JSONs stand, and the §9 line says which nodes remain. The next slot runs only the remaining windows, with no edit.

**Traps already paid for.**
- Item 157's rank-0 `KeyError` (`pec_column=False` with `field_solve=True`): guarded by code step 3.
- κ stays `None` on every spec; no `ans:7` `_stage_*` / `_assemble` / `_write_outputs` call.
- Write only under `logs/ans7_exact/`, never `logs/ans7/` or `logs/ans7_sparse/` (entry 24's files).
- Item 196's seven JSONs are the ×1 comparand and T21's (item 217's) input. Never delete them, and never let a rung run write the unsuffixed names.
- `cell_tags.values` and `assemble_scalar` are rank-local. Rows are computed on every rank; rank 0 writes.
- No ordering comparison or `sqrt` on complex operands.
- A window killed at its timeout is not a result. Run `pgrep -c python3` in the container first, then report; never lengthen the timeout. Run the orphan-`.c` sweep only when that count reads 0.
- Harness calls are bare and whole. No `-k` selection inside the quoted command.
- **No slot writes or edits an `xl` command.**

**Scope.** Covered: one intermediate rung, degree 2, Cu, five nodes, κ off, exact rectangles; plus the ×1 reproduction and the ×1 five-node licence. Not covered:
- the PEC column on the rung;
- (L) on the rung (it needs the field solve, hence the PEC column);
- Route A on the rung;
- h = 0.005, which is entry 28's;
- a record constant for the rung's rows (a follow-up for the review, to avoid an edit after the last window);
- the T13 / T6 re-read and any private comparison (the weekly's);
- any band;
- any claim about the absolute f₀, Q or R_in, or about convergence. Two rungs are not an h-statement.

**Status it can move.**
- The `ANS-7` row gains a "step 6" record sentence (no glyph; the row's Done-when is the weekly's T13 re-read after T12′).
- **`xl` entry 28 PENDING PREREQUISITE → READY.** The implementer landing this item green marks entry 28 READY in the same commit, naming the knobs (`FEM_EM_ANS7_EXACT_ROUTE=B`, `_DEGREE=2`, `_RESOLUTION=0.005`, `_COLUMNS=cu`, `_MHZ=62.0,63.2,64.5,65.8,67.0`) and the rung's measured per-node price. **It does not issue, write or edit entry 28's command**; the next review issues it. READY by the 2026-10-07 18:00 review ⇒ entry 28 takes entry 24's night (10-08).
- **Rule (j):** the review adds a step-6 Done-when sentence to the row.

**Negative result.**
- **(R2) red:** the knob changed the unset path. Restore the JSON backup, revert, BLOCKED.
- **(S5) red:** the Cu-only five-node path is unlicensed at ×1 on exact sheets. BLOCKED, and entry 28 stays PENDING. The review picks six or seven Cu nodes (entry 24's own fallback).
- **(F) or (N-wire) red:** BLOCKED.
- **Rung (B) red, crossings ≠ 1, or the rung's Δ_LOO above `BANDS`:** print the rows as unlicensed and stop, BLOCKED. Entry 28 stays PENDING, and the review re-poses the grid.
- **h-move outside its predicted ranges with every assertion green:** a record, not a red. Report it, and the weekly reads it.
- No `PREDICTED RED:` line.

**Outcome (2026-10-04 06:00 slot): DONE, green.** W0 ×1 reproduction 18 passed (R2 ≤ 2.15e-12 vs record; seven-node Δ_LOO 6.3e-13 / 1.6e-9 / 1.6e-9); W0b `ANS-6` exact mesh 20 passed / 14 skipped (item 193's set plus M7, now running against its since-registered record); (S5) five- vs seven-node Cu ≤ 2.2e-12 / 9.2e-10 / 6.6e-10, five-node Δ_LOO 6.5e-10 / 3.5e-8 / 2.7e-8, B₅ = 1e-6. W1 probe t = 155.4 s ⇒ **rung h = 0.01** (W1b not run): 140 270 cells, 968 016 unknowns, (N-wire) cell ratio 1.736, sheets exact to 8.9e-16; Cu sweep 143.5–159.8 s per node at `-n 12`, mesh 53.6–55.4 s (over the predicted 30–40 s), summed `ru_maxrss` 32.8–36.5 GiB. Rung Cu rows f₀ 62 726 846.660 Hz, R_in 18 313.602 Ω, Q 970.439, k-share 0.99999176, one crossing (Im `Y_in` −1.22e-3 S at 62.0 → +8.03e-4 S at 63.2), Δ_LOO 6.5e-10 / 3.8e-8 / 3.1e-8. h-move vs ×1: f₀ +0.0411 %, R_in −0.6949 %, Q −0.6981 % (inside predictions). (N3) missed as predicted (printed). The box ran ≈ 1.4× slow (W0 Cu sweep 70.9 s vs 44–50 s); W2 used 536 s of 590. Logs `20261004T110446Z_ANS-7-item214-W0-x1.log`, `…111046Z_…-W0b-ans6mesh.log`, `…111415Z_…-W1-probe.log`, `…112040Z_…-W2-rung.log`, `…112945Z_…-W3-rung.log`. `xl` entry 28 → READY (Status line only).

### §9 item 217 — 2026-10-04 (chain T21)

*Drafted by a general-purpose agent for the 2026-10-04 03:00 daily review; the review reviews it.* Licence: §10 chain T21 (weekly 2026-10-03), whose target is quoted, not changed: "the capacitor that puts the k = 1 mode at **64.000 MHz**, from the k = 1 block of the degree-2 exact 4×4 (`PORT-22`'s `C_mode` construction), then the in-model sweep with that C: the zero of Im `Y_in` within the sparse grid's own band of 64 MHz and circuit = field at `REDUCTION_BAND`. Target: circuit ↔ field identity." The ruling binds: exact rectangles, κ off. No band, tolerance, Done-when or validation target moves.

`ANS-7` step 7: tune the exact-sheet F-small birdcage to 64.000 MHz at degree 2 on the ×1 mesh, then check the tuning in-model.
- **Executor:** implementer. Tests only, plus two additive keywords (disclosed).
- **Tier:** heavy, `-n 8`, complex build.
- **Class:** chain step (solve). **Time:** ≈ 18 measured / 45 predicted min.
- **Depends on nothing.** It shares the module with item 214; whichever lands second rebases onto the first.

**Where the pieces live (verified 2026-10-04).**
- **The `C_mode` construction.** It is `TH-17` step 3's, at `tests/validation/test_th17_birdcage_eigenmodes.py:1093–1100` (repeated in `examples/ports/20_birdcage_eigenmodes.py:191–200`): `zk, off = _modal_z(Z)` (`:846`, diag of F Z Fᴴ with the 4-point DFT), z₁ = ½(zk[1] + zk[3]), assert Im z₁ > 0, then `C_mode = 1/(ω₆₄ · Im z₁)`. It is the capacitance that zeroes the k = 1 loop reactance z₁ + 1/(jωC) at 64 MHz. The SPEC's 66.29 pF is that number at degree 1, rounded (`SPEC.md:80–83`).
  - `tests/validation/test_port22_driven_sweep_resonance.py` contains no `C_mode`; its C is `PORT-15`'s `C_tuned` (`build_c_tuned`, `:189`). The chain's attribution to `PORT-22` refers to this construction, and the brief uses `TH-17`'s lines (imported `_modal_z`; the two formula lines restated with that cite).
  - `ports/ring_modes.k_blocks` is the 32-port version and does not apply to four ports.
- **`REDUCTION_BAND`** = 1e-3 is defined at `tests/validation/test_port_lumped_rlc_termination.py:88` (`PORT-14`'s pre-stated band). The exact module imports it through `test_port_circuit_layer_field`.
- **`BANDS`** = 1e-6 on f₀ / R_in / Q is defined at `tests/validation/test_ans7_sparse_route_b.py:93`. This is "the sparse grid's own band": \|f₀/64 MHz − 1\| ≤ 1e-6, i.e. 64 Hz.
- **The degree-2 exact 4×4 is persisted** at seven nodes, not at 64.0 MHz.
  - Location: `logs/ans7_exact/node_{62.0,62.3,63.2,64.5,65.8,66.7,67.0}_p2.json` (gitignored, on disk, written 2026-10-03 07:40–07:54 CDT by item 196 W1–W3).
  - Contents: cells 80 798, the tie, S and Z for both columns, κ `null`. On disk, `route_a_p1.json` / `route_b_p1.json` (item 195) are there too.
  - No 64.0 MHz node exists. TH-17's construction uses a 4×4 solved at 64 MHz, so W1 solves it: one degree-2 node, both columns, ≈ 112 s at `-n 8`.
  - The network Z does not depend on C (only the circuit does), so the seven stored nodes serve any C. If they are absent, stop: re-running item 196 costs ≈ 21 min.

**Code** (`tests/validation/test_ans7_exact_strip.py`; two additive keywords, both defaults byte-identical; no `src/`, no `TH-14` edit):
1. **The additive keywords.**
   - **`c_f=C_F` on `sparse_rows`** (`test_ans7_sparse_route_b.py:321`), forwarded to `route_b_z`.
   - **`degree=1` on `ans:7`'s `_route_a_solve`** (`07_…py:515`, which today hard-codes `degree=1` at `:540`), forwarded to `run_lumped_sheet_port_case`. `ans:7` itself is not run.
   - `metrics.json` and `COMPARISON.md` stay untouched (`git diff --stat`).
   - Each keyword is covered by a wiring control: (N1) for `c_f` and (W) for `degree`.
2. **Knob `FEM_EM_ANS7_EXACT_TUNE`** ∈ {`cmode`, `inmodel`}. Unset ⇒ the new tests skip; anything else raises `ValueError`. It uses `FEM_EM_ANS7_EXACT_DEGREE` (2 here), and `ROUTE` is left unset in both windows. Constants:
   - `TUNE_TARGET_HZ = 64.0e6`;
   - `TUNE_MHZ = (63.9, 64.0, 64.1)`, three `GRID_MHZ` points, so Route A and Route B are read on the same 0.1 MHz bracket.
3. **`ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2`:** item 196's rows at full repr, as in item 214's brief, code step 7. Reuse it if item 214 landed first.
4. **`tune_cmode` fixture (W1).**
   - Fail if the seven ×1 p2 JSONs are absent or inconsistent (degree, cells, κ).
   - Run `_fixture_tie()`, then `th14._build_ladder(2, sigmas=(th14.COPPER_SIGMA,), frequencies_hz=[64.0e6], width_correction_kappa=None, measure_kappa=False, pec_column=True, field_solve=False, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`. Per column, `Z64 = z_from_s(S)`.
   - Apply the construction to Cu (the tuned coil's capacitor) and to PEC (printed). Print `off`, \|z₁ − z₃\|/\|z₁\|, all z_k and C_mode.
   - Printed beside it: the entrywise max \|Z64_interp − Z64\|/\|Z64\|, where `Z64_interp = bary_interp` of the seven stored nodes at 64.0 MHz.
   - Printed diagnostic: **C_in**, the C that zeroes Im `Y_in`(64.0 MHz) on the full solved Cu 4×4, i.e. the root of g(C) = Im[1/`route_b_z`(Z64, 64e6, C)[0]]. It is found by secant from C_mode on a real-valued g, to \|g\| ≤ 1e-15 S in ≤ 20 iterations. Print (C_mode − C_in)/C_in: the k = 1-block construction's own offset, against the full network.
   - Route B at C_mode: `sparse_rows(DEFAULT_NODES_MHZ, z_nodes, columns=("cu", "pec"), c_f=C_mode_cu)` on the stored nodes; Δ_LOO via `_loo` with `c_f`.
   - Rank 0 writes `logs/ans7_exact/tune_p2.json`: C_mode per column, C_in, z_k, Z64 per column, the rows at C_mode, cells, degree. Then a barrier.
5. **`tune_inmodel` fixture (W2).**
   - Read `tune_p2.json` and fail if it is absent or its degree differs. C_mode is read from the JSON, never recomputed.
   - Build geo as `route_a`'s `_geo(GATED_WIDTH_FRACTION)`, plus `_fixture_tie()`.
   - **(W) first:** `_route_a_solve(th14, geo, specs, 64.0e6, "cu")` (defaults: degree 1, C_F). Compare its residual against stored Route B (`route_b_p1.json`, at 64.0 MHz) and against `route_a_p1.json["residual"]["cu"][j64]` (item 195 W2 printed 1.409e-5 there).
   - Then, for each g in `TUNE_MHZ`: `s_a, _ = _route_a_solve(th14, geo, specs, g*1e6, "cu", c_f=C_mode, degree=2)`.
     - The Route B comparand is `_s11_from_z(route_b_z(bary_interp(xn, z_nodes_cu, g − 64.5), g*1e6, C_mode)[1])`. At 64.0 MHz the comparand from the solved Z64 is printed beside it.
     - Z_in,A = `add_c1(_z_from_s11(s_a), ω, C_mode)`.
   - Read f₀,A, R_in,A and Q,A from `parallel_rows(TUNE_MHZ·1e6, Z_in,A)`.
   - (N2) is pure circuit: `_s11_from_z(route_b_z(Z64_cu, 64e6, 0.5·C_mode)[1])` against the 64.0 MHz in-model S11 at C_mode.
   - Rank 0 writes `tune_inmodel_p2.json`.

**Anchors** (asserted; bands imported):
- **(C):** Im z₁(64.0 MHz) > 0 (TH-17's assert); per-column reciprocity ≤ `RECIPROCITY_BAND` and σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` at 64.0 MHz.
- **(T1) Route B at C_mode,** both columns: exactly one − → + zero, k-share > `K_SHARE_FLOOR`, Δ_LOO ≤ `BANDS`. This is the grid's licence at the new C. Also printed: \|f₀,B/64 MHz − 1\| for Cu, plus R_in, Q and k-share at C_mode.
- **(AB) circuit = field,** degree 2, Cu: \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` at all three `TUNE_MHZ` points. This is the chain's stated target identity, and it has never been measured at degree 2. Item 195's degree-1 measurement on exact sheets was 9.9e-6 … 3.11e-5.
- **(T2) the in-model zero:** exactly one − → + crossing of Im `Y_in`,A in (63.9, 64.1) MHz, and \|f₀,A/64 MHz − 1\| ≤ `BANDS["f0_hz"]`.

**Negative controls (rule (e)).**
- **(N1): asserted,** backed by item 196 W3 `:2189` (same fixture, same comparison).
  - At C = `C_F` (66.29 pF) the stored-node rows reproduce `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2` within `NODE_S_RTOL`; this is the wiring of `c_f`.
  - \|f₀(C_F)/64 MHz − 1\| = 2.03e-2 > `N_FACTOR` × `BANDS["f0_hz"]`. Without the tuning the zero sits ≈ 2 × 10⁴ bands from 64 MHz, so (T2) is not vacuous.
- **(W): asserted,** \|Δresidual\| ≤ `NODE_S_RTOL` (absolute) against item 195's stored 64.0 MHz Cu point. It is a deterministic re-solve on the same mesh, with rank-count round-off ≈ 1e-13. A changed default degree or C would move the residual by ≥ 1e-4.
- **(N2): predicted, printed.** The circuit at 0.5·C_mode against the in-model at C_mode, 64.0 MHz Cu, is ≥ 10× `REDUCTION_BAND`. `ans:7`'s control (i) measured 483.6× at degree 1 on ragged sheets (`20260929T100040Z_ANS-7-item68-W3-A.log`), a different degree and sheet, hence predicted. The ceiling is \|ΔS\|/\|S\| ≲ 2, i.e. 2 000× the band, so no saturation.

**Predictions** (written first, printed):
- **C_mode,Cu in [63.4, 63.8] pF.** This is C_F·(62.701082/64)² = 63.62 pF, less the ≈ 6e-4 L_eff dispersion that TH-17 step 3 noted. (C_mode,PEC − C_mode,Cu)/C_mode,Cu ≈ +4e-4. Off-diagonal ≤ 1e-6; interpolated-vs-solved Z64 ≤ 1e-8.
- **(C_mode − C_in)/C_in in [1e-7, 5e-5].** The estimate G²X_bg/(ωC) + 1/(2Q²) gives ≈ 3e-7 – 1e-6 (G ≈ 1/(2 R_in) ≈ 2.7e-5 S; X_bg is the k = 0 / 2 background, of order 10–100 Ω). The degree-1 pair (TH-17's 66.2877 pF → `ans:7` PEC f₀ 63.997103 MHz at 66.29 pF) reads ≈ 2.8e-5 after rounding, but on sheets that are not matched.
- **\|f₀,B/64 MHz − 1\| in [1e-7, 5e-5],** printed only. Tuned rows at C_mode: R_in and Q within ±10 % of item 196's (18 441.75 Ω / 977.26); k-share ≥ 0.9999.
- **(AB) in [5e-6, 3e-4]** at every point.
- **f₀,A − f₀,B ≈ −3e-6 relative:** item 195 measured Route A 63.970953 MHz against Route B 63.971138 MHz on exact sheets at degree 1 (W2 / W3), the same fixture and the same identity.
- **Pre-registered negative.** The in-model zero is expected to miss the 64 Hz band by the A–B gap plus the construction's offset:

PREDICTED RED: (T2) |f₀,A/64 MHz − 1| ≤ BANDS["f0_hz"] range [1e-6, 5e-5]

  The one-crossing half of (T2) is predicted green.

**Windows** (`-n 8`; bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1`; no `ROUTE`; each window < 10 min).
- **Price:**
  - degree-2 PEC / Cu sweeps 52.9–62.2 / 44.3–50.0 s per node at `-n 8` (`20261003T123356Z_ANS-7-item196-W1-p2.log:2025–2148`);
  - one degree-2 in-model drive (`run_lumped_sheet_port_case`, the Route A call), costed as the degree-2 P1 field solve, 44.9 s (`20260921T100022Z_ANS-6-step2-degree2-64MHz.log:924`);
  - mesh 22–23 s;
  - fixed ≈ 78 s.
- **W1** (≈ 180–240 s): the 64.0 MHz node (PEC + Cu ≈ 112 s), two meshes, C_mode, C_in, (C), (T1), (N1).

  `scripts/testing/run_and_log.sh ANS-7-item217-W1-cmode "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_TUNE=cmode timeout -k 30 400 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **W2** (≈ 230–320 s): three degree-2 drives plus one degree-1 drive for (W), two meshes; then (AB), (T2), (W), (N2).

  `scripts/testing/run_and_log.sh ANS-7-item217-W2-inmodel "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_TUNE=inmodel timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **No edit after W1** (rule (i)). W2 is the closing window. With no edit between the two windows, there is no W3.

**Traps already paid for.**
- Items 195 / 196's traps all apply:
  - item 157's rank-0 `KeyError` (W1 runs `pec_column=True`, `field_solve=False`);
  - κ `None` everywhere;
  - no `ans:7` `_stage_*` / `_assemble` / `_write_outputs` call, and `ans:7` is not run;
  - `logs/ans7_exact/` only.
- `_route_a_solve`'s Cu column builds its own Leontovich term on `outer_tags`. Pass `degree` through and change nothing else in it.
- Item 69's drafting error: never feed the sweep's `z_matrix` (the terminated transimpedance) to the Z route. Use `z_from_s(S)` only.
- The secant for C_in runs on the real scalar g. No ordering comparison or `sqrt` on complex operands.
- Never write to item 196's seven node files; this item only reads them.
- Rows are computed on every rank; rank 0 writes, then a barrier.
- A killed window is not a result: `pgrep -c python3`, report, never lengthen the timeout.
- Harness calls are bare and whole. `-s` on every window.

**Scope.** Covered: the ×1 exact-sheet mesh (h = 0.015), degree 2, κ off, the Cu column tuned (PEC printed at Cu's C_mode), Route B on the stored seven-node grid, and Route A at three FineWindow points around 64 MHz. **What T21 does not claim:**
- **no absolute tuned f₀, Q, R_in or C.** T6 and T13 stay open, and the h-rung (T12′, entry 28) has not landed. The tuned C is a ×1-mesh, degree-2 number, and its h-move is unmeasured;
- not the 51-point degree-2 Route A (102 drives, `xl`-class);
- not the eigen route (TH-17's pencil at C_mode is not run);
- not 127.74 MHz (T22);
- no change to the SPEC's C, which stays 66.29 pF for the AED side;
- no AED or Palace comparison.

**Status it can move.**
- The `ANS-7` row gains a "step 7" record sentence: C_mode at degree 2 on exact sheets, the first degree-2 Route A = Route B, and the in-model zero's offset from 64 MHz. No glyph; the row's Done-when is the weekly's T13 re-read after T12′.
- T22 (127.74 MHz) becomes queueable, using either C_mode or C_in per the review's reading of (T2).
- **Rule (j):** the review adds a step-7 Done-when sentence.

**Negative result.**
- **(C) red** (Im z₁ ≤ 0), or (N1) / (W) red: BLOCKED. For (W), revert the keyword.
- **(T1) red** (≠ 1 crossing at C_mode, or Δ_LOO above `BANDS`): BLOCKED, and the review re-poses the grid around 64 MHz.
- **(AB) red with (W) green:** the degree-2 in-model and circuit routes disagree on exact sheets at κ = 0. Print the three residuals, open a known-issues entry ("degree-2 exact-rectangle Route A ≠ Route B"), park on `attempt/*`, BLOCKED. Never apply κ to rescue it.
- **(T2)'s band half red inside [1e-6, 5e-5] with everything else green:** the pre-registered negative. It is a record: the k = 1-block C_mode puts the in-model zero within x of 64 MHz, with C_in and f₀,A − f₀,B printed beside it so the offset is split into construction and route parts. The slot continues and commits DONE (a record). The review decides what T22 tunes with.
- **(T2) red outside the range, or no crossing in (63.9, 64.1) MHz:** the construction does not tune this fixture. Report the measurement and stop, BLOCKED. 

**Outcome (2026-10-04, 09:00 CDT slot): BLOCKED, predicted red, inside range (the pre-registered negative).** The brief says to commit this as DONE (a record), but implementer-run.md step 2 parks a pre-registered negative like an incomplete item. The slot followed the protocol. Code is on `attempt/ANS-7-step7-20261004T1413Z` (`84da27b7`). The logs and test-results rows are on `main`. If the review rules it a record, it lands the branch as is.
- **Logs.** W1 `20261004T140347Z_ANS-7-item217-W1-cmode.log`: 15 passed, 16 skipped, 215 s. W2 `20261004T140738Z_ANS-7-item217-W2-inmodel.log`: 14 passed, 16 skipped, 1 failed, 261 s. Both at `-n 8`, with no kill. The box ran at normal speed (Cu sweep 46.8 s). No edit after W1 (modules 14:03:21Z, W1 start 14:03:47Z).
- **(C) green.** Im z₁ is +39.107 Ω (Cu) and +39.090 Ω (PEC). Reciprocity ≤ 1.94e-15, σ_max 0.999901 / 0.999917, 80 798 cells.
- **(T1) green.** Cu has one crossing at 64 000 470.55 Hz and PEC one at 64 014 337.03 Hz. k-share 0.9999949 / 0.9999974. Δ_LOO ≤ 1.92e-9.
- **(AB) green.** 1.237e-5 / 1.258e-5 / 1.287e-5 at 63.9 / 64.0 / 64.1 MHz. This is the first degree-2 circuit = field reading.
- **(N1) green.** The rows at C_F are within 2.2e-12, and untuned \|f₀/64 − 1\| = 2.0296e-2.
- **(W) green.** Δresidual = 2.77e-14.
- **(T2).** One crossing (green) at f₀,A = 64 000 308.55 Hz. \|f₀,A/64 MHz − 1\| = 4.821e-6 is red, inside the predicted [1e-6, 5e-5].
- **Records.** C_mode,Cu = 63.589141 pF (prediction met). (PEC − Cu)/Cu = +4.406e-4. C_in = 63.590091 pF, so (C_mode − C_in)/C_in = −1.495e-5: inside the range, but well above the 3e-7 – 1e-6 estimate. \|f₀,B/64 − 1\| = 7.352e-6. Rows at C_mode: R_in 18 668.70 Ω (+1.23 %), Q 969.21 (−0.82 %). Route A: R_in 18 739.36 Ω, Q 972.90. f₀,A − f₀,B = −2.53e-6 relative. (N2) 0.4938, which is 494× the band.
- **Offset split.** The construction alone puts the circuit zero at +7.35e-6. The route gap of −2.53e-6 brings the in-model zero back to +4.82e-6. C_in is the candidate that tunes the full network.
- **Prediction missed (printed only).** The k-block off-diagonal is 4.56e-6, against the predicted ≤ 1e-6.

**Review disposition of item 217 (2026-10-04 10:30 daily).** Landed as a record. The three code files were checked out by path from `84da27b7` (`test_ans7_exact_strip.py`, `test_ans7_sparse_route_b.py`, the `ans:7` script). `main` had not touched them since the branch point. The logs and rows were already on `main` (`4d2b8135`). The branch is deleted. The red (T2) is filed in known-issues (2026-10-04), and `BANDS["f0_hz"]` is unchanged. The T21 tests run only under `FEM_EM_ANS7_EXACT_TUNE`. Which capacitor defines "tuned", `C_mode` or `C_in`, is the weekly's question. The brief's "DONE (a record)" wording conflicted with implementer-run.md step 2. The slot's park was correct, and §9's digest now fixes that template.

### §9 item 228 — 2026-10-04 (chain T22)

*Drafted by a general-purpose agent for the 2026-10-04 10:30 daily review; the review reviews it.* Licence: §10 chain T22 (weekly 2026-10-03), quoted, not changed: "T21 at **127.74 MHz** (3 T). Serial on T21." T21's construction and target therefore carry over verbatim with 64.000 → 127.74 MHz: the capacitor from the k = 1 block of the degree-2 exact 4×4 (`TH-17` step 3's `C_mode`), Route B at that C, the in-model sweep with that C, the zero of Im `Y_in` within `BANDS` of the target and circuit = field at `REDUCTION_BAND`. Target: circuit ↔ field identity. The sheets-and-κ ruling binds (exact rectangles, κ off). No band, tolerance, Done-when or validation target moves; `C_in` is a printed record and never a substitute target.

`ANS-7` step 8: T21's tuning construction at 127.74 MHz on the ×1 exact-sheet mesh at degree 2, then the in-model check.
- **Executor:** implementer. Tests only (`tests/validation/test_ans7_exact_strip.py`), plus one additive keyword on a helper of the same module (disclosed). No `src/`, no `ans:7` edit, no `TH-14` edit.
- **Tier:** heavy, `-n 8` (the record width of T21 and item 196), complex build.
- **Class:** chain step (solve). **Time:** ≈ 18 measured / 50 predicted min.
- **Depends on** item 217's T21 section being on `main` (the 2026-10-04 10:30 review lands `attempt/ANS-7-step7-20261004T1413Z`). If `TUNE_TARGET_HZ` is not in the module on `main`, stop and take the next item.

**What T21 measured, and why (T2) is pre-registered red** (item 217's outcome above; `20261004T140347Z_ANS-7-item217-W1-cmode.log:2084–2095`, `20261004T140738Z_ANS-7-item217-W2-inmodel.log:2046–2052`).
- The in-model zero sat at \|f₀,A/64 MHz − 1\| = 4.821e-6 against the 1e-6 band (64 Hz). It has two parts:
  - **construction (c):** the k = 1-block `C_mode` puts the Route B zero at +7.352e-6, i.e. −½ (C_mode − C_in)/C_in with (C_mode − C_in)/C_in = −1.495e-5;
  - **route gap (r):** f₀,A − f₀,B = −2.531e-6 (degree 1 on the same exact sheets: −2.9e-6, item 195 W2 / W3).
- **C_in does not rescue (T2) on its own.** C_in zeroes Im `Y_in` on the *circuit* at the target by construction, so f₀,B(C_in) sits on the target to interpolation error (≈ 1e-10). The in-model zero then sits at the route gap r alone: at 64 MHz that reads ≈ −2.5e-6, which is still outside the 1e-6 band. The slot's hypothesis ("C_in tunes the full network") is right for the circuit and wrong for the in-model zero. This item measures that (a printed record, code step 6) and never swaps C_in in for C_mode.
- **How c and r scale to 127.74 MHz.** Two hypotheses, printed and discriminated by the reading:
  - **c flat.** The absolute k-block coupling over Im z₁ is flat 64 → 128 MHz on `ans:6`'s tracked degree-1 4×4 (≈ 3.9e-4 at both; z_from_s + 4-point DFT on `ans6_…/metrics.json` "rungs", computed for this brief). Then c ≈ +7.4e-6.
  - **c tracks the printed off-diagonal.** That metric grows ×2.47 64 → 128 MHz on `ANS-4`'s degree-2 exact records (`ANS4_EXACT_STRIP_S_P2_X1`: 8.0e-6 → 2.0e-5) and ×4.2 on `ans:6`'s degree-1 records (1.27e-5 → 5.35e-5). Then c ≈ +1.8e-5 … +3.1e-5.
  - **r** is allowed ×1–4 growth with f: [−1e-5, −2.5e-6]. The degree-1 exact (AB) residual grew ×3.1 across the 62–67 MHz window alone (item 195 W2).
  - So \|c + r\| lies in ≈ [0, 4.9e-6] if c is flat and ≈ [0.8e-5, 2.9e-5] if c tracks. A cancellation below 1e-6 is possible under the flat hypothesis. That is green, not red, and is reported as such.

**The capacitor, and its priors** (computed for this brief from tracked, public records — our numbers only).
- **`ans:6` degree-1 ragged 4×4** (`metrics.json` "rungs", Cu): C_mode 66.258 pF at 64 MHz (TH-17's 66.2877 pF, consistent) and 15.811 pF at 128 MHz. Ratio 0.23863, so L₁ grows +4.8 %.
- **`ANS-4` exact records** (σ = 800 S/m solid fixture, not the copper hole; a prior only):
  - the 128/64 C_mode ratio is 0.24544 at degree 1 and 0.24310 at degree 2, so the degree factor on the ratio is ×0.9905;
  - the order move on C_mode grows from −6.0 % (64 MHz) to −6.9 % (128 MHz).
- **Prediction:** 63.589 pF × 0.23863 × 0.9905 × (128/127.74)² ≈ 15.09 pF from the same fixture. The `ANS-4` degree-2 ratio alone gives 15.52 pF. **C_mode,Cu(127.74) in [14.8, 15.6] pF**, so Im z₁ is in [79.8, 84.2] Ω.
- **Inside the capacitor range the fixture has used.** SPEC v1's 15.57 pF on three gaps (item 50, degree 1, ragged, bare feed: a different circuit, which put a series zero at 66.3 MHz) and v2's 66.29 pF bracket it. C enters only the circuit layer (Route B) and the sheet impedance Z_C (Route A). \|Z_C\| ≈ 82 Ω at 127.74 MHz against ≈ 39 Ω at 64 MHz, on a 50 Ω reference, so this is no new conditioning regime.
- **No mode crowding.** k = 2 with C ≈ 15.1 pF resonates near 161 MHz (Im z₂ ≈ 51.5 Ω at 128 MHz on `ans:6`), +26 % against +22 % at 64 MHz. k = 0 is capacitive (z₀ ≈ −567j Ω). k-share ≥ 0.9999 is predicted.

**Code** (`tests/validation/test_ans7_exact_strip.py` only; with the knob unset, every T21 code path and object is used as is, byte-identical):
1. **Knob `FEM_EM_ANS7_EXACT_TUNE_MHZ`** ∈ {unset, `127.74`}. Any other value raises `ValueError`; it is not a free frequency.
   - Unset ⇒ T21's objects verbatim: `TUNE_TARGET_HZ`, `TUNE_MHZ`, `DEFAULT_NODES_MHZ`, `GRID_MHZ`, both columns, and the files `logs/ans7_exact/tune_p2.json` / `tune_inmodel_p2.json` / item 196's `node_*_p2.json`. Do not recompute them from offsets: `64.0 + (62.3 − 64.0)` is not guaranteed to be `62.3` in floating point.
   - `127.74` ⇒ literal constants, registered before W1:
     - `TUNE128_TARGET_HZ = 127.74e6`;
     - `TUNE128_MHZ = (127.64, 127.74, 127.84)`;
     - `TUNE128_NODES_MHZ = (125.74, 126.04, 126.94, 128.24, 129.54, 130.44, 130.74)`. This is `DEFAULT_NODES_MHZ` + 63.74: the same offsets from the target as T21's nodes from 64.0, a 5 MHz window with the target 2.0 MHz above its low edge;
     - `TUNE128_GRID_MHZ = tuple(round(125.74 + 0.1 * k, 6) for k in range(51))`, which contains the three `TUNE128_MHZ` points (k = 19, 20, 21).
   - An active-set helper returns (target, tune points, nodes, grid, out dir) for the knob.
2. **New tune mode `FEM_EM_ANS7_EXACT_TUNE=nodes`** (so `_tune()` accepts `cmode`, `inmodel`, `nodes`). It requires `TUNE_MHZ=127.74` and raises otherwise, so item 196's 64 MHz node files are never re-solved or re-written.
   - **Knob `FEM_EM_ANS7_EXACT_TUNE_NODES`:** a comma list ⊂ `TUNE128_NODES_MHZ`, validated after `round(…, 6)`.
   - **Fixture `tune_nodes`:** `_fixture_tie()`, then `th14._build_ladder(2, sigmas=(th14.COPPER_SIGMA,), frequencies_hz=[m*1e6 …], width_correction_kappa=None, measure_kappa=False, pec_column=("pec" in _cols()), field_solve=False, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`. The field solve is always off, so `COLUMNS=cu` is safe from item 157's trap.
   - **Per node, rank 0 writes** `logs/ans7_exact/t127.74/node_<mhz:.2f>_p2[_cu].json`, using `.2f` because `.1f` would round 125.74 to "125.7". The keys are p2_nodes' keys plus `"target_mhz": 127.74`. Then a barrier.
   - **Asserted:** reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and the tie equal to the ×1 record.
   - **Printed:** Im `Y_in` at C_F (predicted > 0 at every node, above the C_F resonance), the Cu sweep time and the summed `ru_maxrss`.
3. **Additive keyword `grid_mhz=GRID_MHZ` on `_loo`,** forwarded to both `sparse_rows` calls. The default is byte-identical. Without it, `_loo` on the 128 MHz nodes silently evaluates on 62–67 MHz, an extrapolation 60 MHz off the nodes.
4. **`tune_cmode` / `tune_inmodel` read the active set.** Replace every 64-specific literal with the active value:
   - `_stored_x1_p2()` → load the active node set: item 196's at unset; at 127.74, the `t127.74/` files at `_cols()`, checked for degree 2, κ `None`, resolution `None`, one cell count = `ANS6_EXACT_STRIP_HOLE_CELL_RECORD`, and `target_mhz`;
   - `assert abs(rec["f_hz"] − TUNE_TARGET_HZ) < 1.0` → the active target;
   - `_loo(…, c_f=c_cu)` / `sparse_rows` → `grid_mhz=` the active grid;
   - Route B rows and (T1) over `_cols()`;
   - the (T2) bracket `63.9e6 < f₀ < 64.1e6` → (tune points[0], tune points[−1]);
   - the prediction prints → a dict keyed by target. The 64 MHz entries are T21's text unchanged; the 127.74 entries are this brief's.
   - Rank 0 writes `tune_p2.json` / `tune_inmodel_p2.json` under the active out dir. T21's two JSONs are never overwritten.
   - The cmode window still solves the target node with **both** columns (`pec_column=True`), as T21 did: (C) on both columns and C_mode,PEC printed.
   - `bary_interp` keeps its `− 64.5` centring. The barycentric weights depend only on node differences and the grid points hit nodes exactly after the same `round(…, 6)`, so no edit is needed.
5. **Controls.**
   - **(N1)** is unchanged and still asserted: `c_f=C_F` on item 196's stored 64 MHz nodes reproduces `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2`. It is pure and costs milliseconds, and it is the `c_f` wiring at either target.
   - **(N0)** is new, at 127.74, predicted and printed: at C_F the 127.74 MHz Cu grid has no − → + zero of Im `Y_in`. The k = 1 zero is at 62.7 MHz and the k = 2 zero near 76 MHz. Without the tuning there is no resonance in the window.
   - **Printed:** the sensitivity \|f₀,B(C_mode·(1 + 1e-5))/f₀,B(C_mode) − 1\| (≈ 5e-6), so the band is shown to be resolvable.
6. **C_in record in W4** (printed, never asserted, never a target).
   - One extra degree-2 drive at the target with c_f = C_in (read from the active `tune_p2.json`).
   - Print Im `Y_in`,A(target; C_in) and the estimated in-model offset δf/f ≈ −Im `Y_in`,A / (f · `dbdf_s_per_hz`), the slope taken from the C_mode rows' bracket. Label it "estimate (one drive + slope)".
   - Also print \|S11_A(C_in) − S11_B(C_in)\|/\|S11_B\| on the solved target 4×4.
   - The (W) control stays exactly as in T21: the degree-1 default drive at 64.0 MHz against item 195's stored residual.

**Anchors** (asserted; every band imported, none restated): T21's, at the active target.
- **(C):** Im z₁(127.74) > 0. Per column, reciprocity ≤ `RECIPROCITY_BAND` and σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` at the target. Cells equal the stored nodes' cells.
- **Per node (W1, W2):** reciprocity, passivity, and the (F) tie.
- **(T1)** Route B at C_mode,Cu on the 127.74 MHz Cu grid: exactly one − → + zero, k-share > `K_SHARE_FLOOR`, Δ_LOO ≤ `BANDS`. This is the new grid's licence.
- **(AB)** circuit = field at degree 2, Cu: \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` at all three `TUNE128_MHZ` points. This is the chain's target identity and has never been measured at 3 T.
- **(T2)** exactly one − → + crossing of Im `Y_in`,A in (127.64, 127.84) MHz, and \|f₀,A/127.74 MHz − 1\| ≤ `BANDS["f0_hz"]`.
- **(N1)** and **(W)**, as in T21.

**Negative controls (rule (e)).**
- **(N1): asserted,** backed by item 217 W1 `:2096–2100` (≤ 4.4e-13).
- **(W): asserted,** backed by item 217 W2 `:2045` (2.77e-14).
- **(N0): predicted, printed.** It is unmeasured at 3 T, so it is not asserted.
- **(N2): predicted, printed.** The circuit at 0.5·C_mode against the in-model at C_mode at 127.74 MHz is ≥ 10× `REDUCTION_BAND`. T21 measured 494×; the ceiling \|ΔS\|/\|S\| ≲ 2 is 2 000× the band, so there is no saturation.

**Predictions** (written first, printed):
- **C_mode,Cu** in [14.8, 15.6] pF. (C_mode,PEC − C_mode,Cu)/C_mode,Cu in [+1e-4, +1e-3] (T21: +4.4e-4).
- **k-block off-diagonal** in [5e-6, 3e-5]. Interpolated-vs-solved Z at the target ≤ 1e-8.
- **(C_mode − C_in)/C_in:** negative, \|·\| in [3e-6, 1e-4].
- **\|f₀,B/127.74 MHz − 1\|** in [3e-6, 5e-5], printed with its sign.
- **(T1):** Δ_LOO ≤ 1e-9 (half T21's relative span); k-share ≥ 0.9999. Q at C_mode in [400, 2 500], printed with R_in.
- **(AB)** in [5e-6, 3e-4] at every point.
- **f₀,A − f₀,B relative** in [−1.5e-5, −1e-6].
- **C_in record:** \|δf/f\| in [1e-7, 2e-5]. Predicted outside `BANDS` too, because the route gap alone was 2.5e-6 at 64 MHz.
- **Hypothesis readout,** printed: \|f₀,A/127.74 − 1\| < 1e-5 favours "c flat"; > 1.2e-5 favours "c tracks the off-diagonal"; between, undecided.
- **Pre-registered negative:**

PREDICTED RED: (T2) |f₀,A/127.74 MHz − 1| ≤ BANDS["f0_hz"] range [1e-6, 5e-5]

  The one-crossing half of (T2) is predicted green. A reading ≤ 1e-6 (cancellation of c and r) is green, a record, and not a contradiction.

**Windows** (`-n 8`; bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE` unset; `FEM_EM_ANS7_EXACT_COLUMNS=cu` on all four, so the node files and their readers agree).
- **Price** (measured at `-n 8`, ×1, 80 798 cells; the 4×4 factorisation cost does not depend on f):
  - Cu sweep 44.3–50.0 s and PEC sweep 52.9–62.2 s per node (`20261003T123356Z_ANS-7-item196-W1-p2.log:2025–2148`); ladder mesh 22–23 s;
  - window overhead ≈ 86 s (item 217 W1: 215 s total, 128.6 s solve);
  - one degree-2 drive 41.8–50.3 s (item 217 W2 `:2052`).
  - A ×1.4 slow box (item 214) is inside every window's 590 s.
- **W1** (≈ 300 s; ≈ 415 s on a slow box): nodes 125.74, 126.04, 126.94, 128.24.

  `scripts/testing/run_and_log.sh ANS-7-item228-W1-nodes "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_TUNE=nodes FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_NODES=125.74,126.04,126.94,128.24 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **W2** (≈ 250 s): W1's command with `TUNE_NODES=129.54,130.44,130.74`, labelled `ANS-7-item228-W2-nodes`.
- **W3** (≈ 215 s, item 217 W1's measured price): `TUNE=cmode`, with `TUNE_NODES` dropped, labelled `ANS-7-item228-W3-cmode`. It runs (C), (T1), (N1), (N0), C_in and the predictions.
- **W4** (≈ 310 s; ≈ 435 s on a slow box): `TUNE=inmodel`, labelled `ANS-7-item228-W4-inmodel`. It runs (W), three C_mode drives, one C_in drive, (AB), (T2) and (N2). This is the closing window.
- **Totals:** ≈ 1 075 s of windows (≈ 18 min).
- **No edit after W1** (rule (i)). If the slot reaches minute 50 before W4, park per rule (d): the node JSONs and `tune_p2.json` stand, and the next slot runs only the remaining windows, with no edit.
- **Review option:** for the PEC grid at 3 T, drop `COLUMNS=cu` and split the nodes 2 + 2 + 2 + 1. That is +≈ 420 s and one more window. It is not the default here because (T1)'s PEC half is not the target.

**Traps already paid for.**
- All of item 217's traps:
  - κ `None`;
  - no `ans:7` `_stage_*` / `_assemble` / `_write_outputs` call, and `ans:7` is not run;
  - `z_from_s(S)`, never the sweep's `z_matrix`;
  - the C_in secant on the real scalar g;
  - rows computed on every rank, rank 0 writes, then a barrier;
  - a killed window is not a result: `pgrep -c python3`, report, never lengthen the timeout;
  - bare, whole harness calls, no `-k` selection.
- **New here:**
  - **Never write** T21's `tune_p2.json` / `tune_inmodel_p2.json` or item 196's `node_*_p2.json`. All 127.74 MHz output goes under `logs/ans7_exact/t127.74/`.
  - **Don't route the 127.74 nodes through `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ`.** `_nodes_env` validates against `DEFAULT_NODES_MHZ`, and `p2_asm` / (R2) would assemble the 127.74 MHz nodes at C_F, find no crossing and go red.
  - **The `_loo` default grid**, code step 3.
  - **The hard-coded 64 MHz literals,** code step 4. A missed one shows up as a (T2) bracket or `f_hz` assert failing on a correct measurement.
- **Item 217's reconciliation:** a (T2) red inside the range parks per implementer-run.md step 2 (journaled `blocked`, "predicted red, inside range"), and the slot continues. The review lands the branch as a record, as it did for item 217. The brief does not say "commit DONE".

**Scope.** Covered: the ×1 exact-sheet mesh (h = 0.015), degree 2, κ off, the Cu column tuned at 127.74 MHz (PEC printed at the target node), Route B on a new seven-node Cu grid, and Route A at three points. **Not covered:**
- any absolute tuned f₀, Q, R_in or C at 3 T;
- any h-statement (T12′ is 64 MHz only);
- the PEC grid at 3 T;
- (L) at 127.74 MHz;
- the 51-point Route A;
- the eigen route;
- any SPEC change (the SPEC is a 64 MHz AED case, and no 3 T AED or Palace case exists, so there is no private reading);
- any choice of construction for later steps (C_mode / C_in / an in-model C), which is the weekly's.

**Status it can move.**
- The `ANS-7` row's chain-T record gains a "step 8" sentence: C_mode at 127.74 MHz (degree 2, exact), the first 3 T circuit = field reading, the in-model zero's offset and its c / r split, the C_in record, and which scaling hypothesis the reading favours.
- No glyph: `ANS-7` stays 🟡 per its Done-when (the weekly's T13 re-read after T12′).
- It hands the weekly the evidence to choose a tuning construction. That choice is not made here.
- **Rule (j):** the review adds a step-8 Done-when sentence.

**Negative result.**
- **(C) red, or (N1) / (W) red:** BLOCKED.
- **A per-node anchor red in W1 / W2:** BLOCKED; the node JSONs stand.
- **(T1) red** (≠ 1 crossing at C_mode, Δ_LOO above `BANDS`, or k-share at or below the floor): BLOCKED. The review re-poses the grid around 127.74 MHz.
- **(AB) red with (W) green:** print the residuals and open a known-issues entry ("degree-2 exact-rectangle Route A ≠ Route B at 127.74 MHz"). Park, BLOCKED. Never apply κ.
- **(T2)'s band half red inside [1e-6, 5e-5], everything else green:** the pre-registered negative. It is a record, with the c / r split, the C_in record and the hypothesis readout printed. Park per implementer-run.md step 2; the slot continues.
- **(T2) red above 5e-5, or no crossing in (127.64, 127.84) MHz:** the construction does not tune this fixture at 3 T to its predicted class. Report the measurement and stop, BLOCKED.

**Outcome (2026-10-04 15:00 slot; landed as a record by the 18:00 daily review).** The pre-registered negative fired: (T2) = 1.390e-6, inside [1e-6, 5e-5], every other anchor green (`20261004T202814Z_ANS-7-item228-W1-nodes.log`, `…203339Z…-W2-nodes.log`, `…203817Z…-W3-cmode.log`, `…204206Z…-W4-inmodel.log`). c = +3.653e-6, r = −2.263e-6 (flat in f against −2.53e-6 at 64 MHz). `C_mode` − `C_in` = −7.754e-6 of `C_in`. (AB) 1.589 / 1.601 / 1.616e-5 at 127.64 / 127.74 / 127.84 MHz, the first 3 T circuit = field reading at degree 2. The review took `test_ans7_exact_strip.py` and the four logs by path from `attempt/ANS-7-step8-20261004T2049Z` (`22de12f4`; `main` had not touched the module since the branch point `a0d675c5`), added the four test-results rows by hand, and deleted the branch. Known-issues 2026-10-04 (T21) carries the 3 T reading. **No step-8 Done-when sentence is added under rule (j)**: the step moves no glyph, and the construction choice is the weekly's.

### §9 item 230 — 2026-10-04 (chain T20; example-runner)

*Drafted by a general-purpose agent for the 2026-10-04 10:30 daily review; the review reviews it.* Licence: §10 chain T20 (weekly 2026-10-03), quoted, not changed: "`ans:7` re-issued at degree 2 on exact rectangles — the example, its `metrics.json` and guide carry the degree-2 rows; the degree-1 rows are labelled superseded. `example-runner`. Serial on nothing." The sheets-and-κ ruling binds. No band, tolerance, Done-when or validation target moves. The SPEC (`SPEC.md`, C = 66.29 pF, the operator's AED case definition) is not edited.

`ANS-7` step 9: the `ans:7` example gains degree-2 exact-rectangle stages. Its tracked `metrics.json`, `COMPARISON.md` and guide carry the degree-2 rows, and the degree-1 ragged κ-corrected rows stay in place, labelled superseded.
- **Executor:** `example-runner`. Example script, same-stem guide, `metrics.json`, `COMPARISON.md`. No test, `src/` or SPEC edit.
- **Tier:** heavy, `-n 8` (the width of the imported record), complex build.
- **Class:** example re-issue (solve). **Time:** ≈ 16 measured / 50 predicted min.
- **Depends on nothing** in the queue. It is independent of item 228: different files, C_F around 62.7 MHz against C_mode at 127.74 MHz, and no shared output path. It uses item 217's additive `_route_a_solve(degree=)` and the `_exact_geo` helper, landed on `main` by this review. If they are absent, stop.

**What the degree-2 exact rows are and where they came from.**
- **Route B, degree 2, exact rectangles, κ off, seven nodes** (`DEFAULT_NODES_MHZ`), assembled by `sparse_rows` on the 51-point `GRID_MHZ`. These are item 196's W1–W3 at `-n 8` (`20261003T123356Z_…-W1-p2.log`, `…124111Z_…-W2-p2.log`, `…124814Z_…-W3-p2.log`).
- **The record** is `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2` (`tests/validation/test_ans7_exact_strip.py`; item 196 W3 `:2189,2193`; version tag complex 0.11 image, `-n 8`; the (1*) licence). It was reproduced by item 214's (R2) at ≤ 2.15e-12 and by item 217's (N1) at ≤ 4.4e-13.
- **The loss identity (L)** at 64.5 MHz, κ off: 1.274e-5 ≤ `REDUCTION_BAND` (item 196 W2).
- **Degree-2 Route A = Route B** was measured once, at C_mode around 64 MHz (item 217 W2: 1.237e-5 … 1.287e-5). It is an identity at an imported band, not a registered record. At the SPEC's C_F it has not been measured; this item measures it at three points.

**Code** (`examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py`; additive stages, with B1 / B2 / A unchanged except code step 5):
1. **New stages** `FEM_EM_ANS7_STAGE` ∈ {`X1`, `X2`, `X3`, `X4`, `XA`}, beside B1 / B2 / A.
   - Stage files go to `/workspace/logs/ans7_p2/` (gitignored), with the `FEM_EM_ANS7_RESULTS_DIR` override honoured. They never go to `logs/ans7/` (the B stages' files), `logs/ans7_exact/` (the tests' files) or `logs/ans7_sparse/` (`xl` files).
   - An X stage raises if any `FEM_EM_ANS7_EXACT_*` or `FEM_EM_ANS7_SPARSE_*` variable is set, as the existing `FEM_EM_ANS6_*` guard does.
   - **Imports are lazy, inside the X-stage functions.** `test_ans7_sparse_route_b` re-executes this script's module body on import (`_load(SCRIPT, "ans7_for_sparse")`), so a module-level import makes a partially initialised cycle. Imported, never restated:
     - from `test_ans7_exact_strip`: `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2`, `_fixture_tie`, `_assert_tie`, `_exact_geo`, `_loo`;
     - from `test_ans7_sparse_route_b`: `DEFAULT_NODES_MHZ`, `GRID_MHZ`, `BANDS`, `NODE_S_RTOL`, `N_FACTOR`, `K_SHARE_FLOOR`, `sparse_rows`, `bary_interp`, `_loss_pair`, `_p_inc_w`;
     - `GATED_WIDTH_FRACTION`, `RECIPROCITY_BAND`, `PASSIVITY_SIGMA_TOLERANCE`, `REDUCTION_BAND`, `DISCRETE_IDENTITY_RTOL`.
2. **X1–X4: node solves.**
   - X1 = {62.0, 62.3} plus (F) via `_fixture_tie()`; X2 = {63.2, 64.5}; X3 = {65.8, 66.7}; X4 = {67.0} plus the assembly.
   - Each is item 196's call: `th14._build_ladder(2, sigmas=(th14.COPPER_SIGMA,), frequencies_hz=[…], width_correction_kappa=None, measure_kappa=False, pec_column=True, field_solve=(64.5 in nodes), leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`.
   - Per node, rank 0 writes `node_<mhz:.1f>.json` (S, `z_from_s` Z, reciprocity, σ_max, cells; at 64.5 MHz also the power and identity residual), then a barrier.
   - Per node, asserted: reciprocity, passivity, and at 64.5 MHz the copper identity residual and **(L)** \|L_S − L_F\| ≤ `REDUCTION_BAND`.
3. **X4: the assembly.**
   - It reads all seven node files (stop if any is absent, naming the missing stage) and asserts one cell count = `ANS6_EXACT_STRIP_HOLE_CELL_RECORD` within `CELL_COUNT_BAND`.
   - `_loo(DEFAULT_NODES_MHZ, z_nodes)` (C_F default) gives the rows per column and Δ_LOO.
   - Asserted: one − → + zero per column, k-share > `K_SHARE_FLOOR`, Δ_LOO ≤ `BANDS`, and **(R2x)**.
   - It then writes **`metrics.json` additively.** It loads the tracked file and snapshots `json.dumps(v, sort_keys=True)` for every pre-existing top-level key. It adds:
     - `"degree2_exact"`: degree 2, κ `null`, sheets "exact CAD rectangles (operator ruling 2026-10-02)", cells, unknowns, `mpi_ranks`, nodes, grid, rows per column (f₀, R_in, Q, k-share, crossings, bracket), Δ_LOO, the (R2x) relative differences, (L), and per-stage seconds;
     - `"superseded"`: names the degree-1 keys (`rows`, `rows_route_a`, `points`, `kappa`, `kappa_per_sheet`, `n_cells`, `degree`, `band`, `controls`, `route_comparison_a_minus_b`), says they are superseded by `degree2_exact` (§10 T20), and records that they are kept in place for their readers (`test_ans7_sparse_route_b`, `test_wf5_q_two_ways`, `test_ans6_exact_strip`, `test_port23_skrf_route_b`).
   - It asserts **(K)**: every pre-existing key re-serialises identically.
   - It regenerates `COMPARISON.md` from `metrics.json`. The degree-2 table comes first: Ours Cu (B), Ours PEC (B), Ours Cu (A) (blank until XA), AED blank. Below it sits "Superseded — degree 1, ragged sheets, κ-corrected (stage B2, 2026-09-29)", the existing table's numbers read from `metrics.json` and not restated in code. "Not claimed" is unchanged.
4. **XA: Route A at degree 2, Cu, at C_F.**
   - The geometry is `_exact_geo(th14, comm, GATED_WIDTH_FRACTION)`. The drives are `_route_a_solve(th14, geo, specs, f, "cu", degree=2)` at 62.6, 62.7, 62.8 MHz, three `GRID_MHZ` points bracketing the degree-2 Cu f₀ of 62.701 MHz.
   - The comparand is `_s11_from_z(route_b_z(bary_interp(xn, z_cu_nodes, g − 64.5), f, C_F)[1])` on the X node files.
   - Asserted: **(AB)** at every point, exactly one − → + crossing of Im `Y_in`,A in (62.6, 62.8) MHz, and **(K)**.
   - Rows A come from `parallel_rows` and are printed. Control (i) is printed.
   - It writes `degree2_exact["route_a"]` (points, S11_A, S11_B, residuals, rows A) and regenerates `COMPARISON.md`.
5. **Stage B2's `_write_outputs` carries over** the `degree2_exact` / `superseded` keys if present (additive, two lines). Re-running B2 would otherwise silently drop the degree-2 block. B2 is not run here.
6. **Docs.**
   - The script docstring gets an X-stage paragraph.
   - The guide gets a new first section, "Degree 2 on exact rectangles (current; §10 T20)": the stages, commands and anchors.
   - §§2–3's degree-1 text moves under "Superseded: degree 1, ragged sheets, κ-corrected (kept as the record)".
   - "Not claimed" gains "no h-statement (T12′); no absolute f₀ / Q / R_in".
   - Setup figure: no re-render, since no X stage calls `write_setup_figure`.

**Anchors** (asserted; every gate record and band imported, none restated):
- **(R2x):** the example's seven-node degree-2 rows reproduce `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2` (Cu and PEC; f₀, R_in, Q) within `NODE_S_RTOL`, at `-n 8`. At any other width it prints and does not assert (the `EX-56` rule).
- **(F)** in X1.
- **Per node (B)** in X1–X4.
- **(L)** in X2.
- **The assembly** in X4: one crossing per column, k-share > floor, Δ_LOO ≤ `BANDS`.
- **(AB)** plus the Route A crossing in XA.
- **(K)** in X4 and XA.

**Negative controls (rule (e)).**
- **(N-deg): asserted.** \|f₀,X(Cu)/f₀,superseded(Cu) − 1\| > `N_FACTOR` × `BANDS["f0_hz"]` (1e-5), where the superseded value is `metrics.json` `rows.cu.f0_hz`, the ragged degree-1 63.983031 MHz.
  - Backed by item 196's −2.0 % against the ragged rows (−1.985 % against the exact degree 1).
  - A stage that silently ran degree 1 or ragged sheets would sit within 0.02 % (item 195's exact-vs-ragged degree-1 move). The separation is ≈ 2 000×, and no ceiling applies.
- **(K):** structural, asserted. A writer that rewrites `metrics.json` wholesale, as `_write_outputs` does, fails it.
- **Control (i): predicted, printed.** The circuit at 0.5·C_F against the in-model at C_F at 62.7 MHz, Cu, is ≥ 10× `REDUCTION_BAND`. The measured analogue is item 217's (N2) at 494×; the ceiling is 2 000×.

**Predictions** (written first, printed):
- **(R2x):** ≤ 1e-11 (item 214's W0: 2.15e-12).
- **Δ_LOO:** ≈ 6.3e-13 / 1.6e-9 / 1.6e-9.
- **(L):** in [1e-6, 1e-4] (item 196: 1.274e-5).
- **(AB):** in [5e-6, 3e-4] (item 217: 1.24–1.29e-5 at C_mode).
- **f₀,A − f₀,B relative:** in [−1e-5, −1e-6] (item 217: −2.53e-6).
- **Route A Cu:** R_in and Q within ±2 % of the Route B rows.

**Price** (`ans:7` run history and the degree-2 node prices).
- **`ans:7` at degree 1, `-n 4`:**
  - B1 353 s, B2 331 s, A 367 s (`20260929T093403Z_ANS-7-item69-W1-B1.log`, `…094009Z_…-W2-B2.log`, `…100040Z_ANS-7-item68-W3-A.log`);
  - item 79's re-run: B2 375 s, A 427 s (`20260929T172027Z`, `…172651Z`);
  - the pin test 9 s (`…173423Z_ANS-7-item79-W3-pin.log`); the census 1 s.
- **Degree-2 ×1, `-n 8`:**
  - PEC 52.9–62.2 s plus Cu 44.3–50.0 s per node; ladder mesh 22–23 s;
  - a degree-2 P1 field solve ≈ 45 s;
  - one degree-2 drive 41.8–50.3 s (item 217 W2 `:2052`);
  - item 196's three-node windows ran 405–421 s.
- **Stages:** X1 ≈ 300 s, X2 ≈ 365 s, X3 ≈ 280 s, X4 ≈ 170 s, XA ≈ 200 s. Each is under 590 s on a ×1.4 slow box (X2 ≈ 510 s). Total ≈ 1 315 s (≈ 22 min).

**Windows** (`-n 8`; bare, whole, foreground; Bash `timeout: 660000`; `FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1`; no `FEM_EM_ANS6_*` / `FEM_EM_ANS7_EXACT_*` / `FEM_EM_ANS7_SPARSE_*`; spawn and run foreground, per the 2026-09-01 rule).
1. **Censuses before:** docrefs and setup-figure (`check_example_doc_references.py`, `check_example_setup_figures.py --next 3`). Predicted unchanged by this item.
2. **X1:**

   `scripts/testing/run_and_log.sh ANS-7-item230-X1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_STAGE=X1 timeout -k 30 590 mpiexec -n 8 python3 examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py'"`
3. **X2, X3, X4, XA:** the same command with the stage swapped, labelled `ANS-7-item230-X2` … `-XA`. No edit after X1 (rule (i)).
4. **Pin and unit tests,** run on the committed `metrics.json`. The path is already in `LISTED_EXAMPLE_ARTIFACTS`, so the `ANS-1` note is satisfied with no list edit. The rule is still that the commit that changes `metrics.json` carries the green pin run. No other new tracked file under `examples/` is allowed; if one appears, it must be listed in the same commit.
   - The pin test is `tests/unit/test_doc_reference_exit_codes.py` (predicted 20 passed).
   - `tests/unit/test_ans7_private_comparison.py` runs because the `metrics.json` shape changed (predicted 4 passed).
   - Both go through the harness at `-n 1`, `timeout -k 30 60`, labelled `ANS-7-item230-pin` / `-unit`.
5. **Censuses after:** predicted equal to the before readings.
6. **Commit the code, guide, `metrics.json`, `COMPARISON.md`, logs and test-results rows together.** Check `git diff --stat` shows no `COMPARISON_private.md`, `aed_results/` or SPEC change.

**Traps already paid for.**
- **Lazy imports,** code step 1.
- **The write path:** `_write_outputs` rewrites `metrics.json` wholesale. Code step 5 and (K) guard against that.
- **The degree-1 key set has readers.** `kappa`, `n_cells`, `points` and `rows` are read by four test modules. They are never renamed or moved; they are labelled.
- **The ragged degree-1 XDMF and figure paths stay untouched.**
- **Never run B1 / B2 / A here.** They would re-solve degree 1 for ≈ 18 min and rewrite the degree-1 keys.
- **Item 157's rank-0 `KeyError`:** X2 runs `field_solve=True` with `pec_column=True`.
- **Item 69's drafting error:** `z_from_s(S)`, never the sweep's `z_matrix`.
- **κ `None` on every spec.**
- **AED values never enter a tracked file;** `FEM_EM_ANS_SKIP_PRIVATE=1` on every window.
- **A killed window is not a result:** `pgrep -c python3`, report, never lengthen the timeout. The node files from finished stages stand; the next slot resumes at the missing stage with no edit.
- **Harness calls are bare and whole;** no pipe.

**Scope.** Covered: the example re-issued at degree 2 on exact rectangles, κ off, ×1 mesh, SPEC C = 66.29 pF, Route B seven-node on both columns, Route A at three Cu points. **Not covered:**
- the 51-point degree-2 Route A (102 drives, `xl`-class);
- the PEC Route A;
- the h-rung (T12′);
- any tuned C (T21 / T22 stay in the test module);
- the private writer's degree-2 rows (the weekly's private re-read);
- the SPEC;
- any AED or Palace reading;
- any absolute claim.

**Status it can move.**
- §10 T20 closes; the weekly removes it from the live list.
- The `ANS-7` row gains a "step 9 (T20)" record sentence: `ans:7`'s tracked outputs now carry the degree-2 exact rows, with (R2x), (L) and (AB) re-asserted through the example, and the degree-1 rows labelled superseded.
- No glyph: `ANS-7` stays 🟡 per its Done-when.
- **Rule (j):** the review adds a step-9 Done-when sentence.

**Negative result.**
- **(R2x) red at `-n 8`:** the example's path is not item 196's. Print the differences. Revert `metrics.json` / `COMPARISON.md` (`git checkout` those two paths only), park, BLOCKED. Never re-register the constant.
- **(F), per-node (B), (L) or the assembly red:** BLOCKED. This contradicts a landed record, so open a known-issues entry if it reproduces.
- **(AB) red with (N-deg) and (K) green:** degree-2 Route A ≠ Route B at C_F on exact sheets. Open a known-issues entry and park, BLOCKED. Never apply κ.
- **(K) or (N-deg) red:** revert both tracked outputs, BLOCKED.
- No `PREDICTED RED:` line.

**Outcome (2026-10-04 16:30 CDT slot): complete.** `-n 8`, ×1, 80 798 cells, 570 332 unknowns, no kill, script unedited after X1 (rule (i)). X1 `20261004T215054Z_ANS-7-item230-X1.log` Status 0, 251 s; X2 `…215514Z_…-X2.log` 323 s; X3 `…220047Z_…-X3.log` 225 s; X4 `…220443Z_…-X4.log` 133 s; XA `…220710Z_…-XA.log` 186 s; pin `…221037Z_…-pin.log` 20 passed; unit `…221048Z_…-unit.log` 4 passed. (R2x) Cu f₀ 1.3e-14, R_in 1.6e-12, Q 3.1e-12; PEC f₀ 9.0e-14, R_in 4.2e-13, Q 2.5e-11 (asserted ≤ `NODE_S_RTOL`; the printed prediction ≤ 1e-11 was missed on PEC Q only). Δ_LOO 8.3e-13 / 1.64e-9 / 1.53e-9 (predicted 6.3e-13 / 1.6e-9 / 1.6e-9). (L) 1.274e-5 at 64.5 MHz (1.239e-5 at 63.2 MHz, also solved because the ladder's field solve covers the stage). (AB) 1.221e-5 / 1.242e-5 / 1.272e-5; Route A one − → + crossing at 62.700923 MHz, f₀,A − f₀,B = −158.96 Hz (−2.535e-6, predicted [−1e-5, −1e-6]); R_in and Q A/B − 1 = +0.37 % / +0.38 %. (N-deg) 2.0036e-2 = 2 003.6× the floor. (K) green at X4 (25 keys) and XA. Control (i) 4.894e-1 = 489× `REDUCTION_BAND` (predicted ≥ 10×). Census before and after identical (docrefs dead=0 guide=0 stale=38 exit=2; setup figures examples=83 ok=82 missing=1 broken=0), as predicted. The SPEC, `COMPARISON_private.md` and `aed_results/` are untouched.

### §9 item 262 — 2026-10-06 (chain T21′)

*Drafted by a general-purpose agent for the 2026-10-06 03:00 daily review; the review reviews it.* Licence: §10 chain T21′ (grep "T21′. `ANS-7` step 7′") under the weekly ruling "T21 / T22 — \"tuned\" means the in-model zero" (PROJECT_PLAN.md ≈ :4468–4481), quoted, not changed: "C_in is the in-model root (Route A …); asserted, the in-model zero of Im `Y_in` at C_in sits at 64.000 … MHz within the sparse grid's own frequency resolution, measured from T21's … node residuals and pre-registered in the brief before the run (it must exceed r …), never 1e-6; circuit = field at `REDUCTION_BAND` at C_in; printed, c = C_mode / C_in − 1 … and r." Exact rectangles, κ off. `BANDS` is not edited; the new band is a resolution on a new quantity, per the ruling.

`ANS-7` step 7′: drive the exact-sheet F-small birdcage at degree 2 on the ×1 mesh with **C_in** (item 217's root of Im `Y_in`(64.0 MHz) on the solved Cu 4×4) at T21's three FineWindow points, and read the in-model zero.
- **Executor:** implementer. Tests only (`tests/validation/test_ans7_exact_strip.py`); no `src/`, `ans:7`, `TH-14` or `test_ans7_sparse_route_b.py` edit.
- **Tier:** heavy, `-n 8` (T21's record width), complex build. **Class:** chain step (solve). **Time:** ≈ 17 measured / ≈ 35 predicted min (15 fixed + ≈ 10 code + W1 ≈ 4.5 + ≈ 5 spare).
- **Depends on nothing.** Item 217's `logs/ans7_exact/tune_p2.json` (64 MHz; carries `c_in` = 6.359009101426113e-11 F and `c_mode`) and item 196's seven ×1 p2 node JSONs are on disk (checked 2026-10-06). If `tune_p2.json` is absent, first run item 217's W1 command verbatim (215 s, `20261004T140347Z_ANS-7-item217-W1-cmode.log`), with no edit.

**The resolution band, pre-registered** (measured, not chosen).
- **Source: T21's node residuals.** The (AB) field-vs-circuit residual \|S11_A − S11_B\|/\|S11_B\| at C_mode at the three nodes 63.9 / 64.0 / 64.1 MHz: 1.237e-5 / 1.258e-5 / 1.287e-5 (`20261004T140738Z_ANS-7-item217-W2-inmodel.log:2046–2048`). This is the relative distance at which the grid's field read and its circuit read cannot be told apart at the nodes; a zero read off those nodes is resolved no finer.
- **Band:** `T21P_ZERO_RESOLUTION = 1.287e-5`. It is the maximum, at logged precision, registered as a literal constant in the module before W1 and never recomputed from the new run. Asserted: \|f₀,A(C_in)/64 MHz − 1\| ≤ `T21P_ZERO_RESOLUTION`.
- **It exceeds r.** r = f₀,A/f₀,B − 1 = −2.531e-6 at 64 MHz (`…140738Z…:2050`), −2.535e-6 at 62.70 MHz (item 230 X-windows, above) and −2.263e-6 at 127.74 MHz (`20261004T204206Z_ANS-7-item228-W4-inmodel.log:2058–2059`). The band is ×5.1 / ×5.7 \|r\|, ≈ 1e-5 as the ruling expected, and 12.9× the retired 1e-6.
- **It is resolvable and not vacuous.** The untuned zero (C_F) sits 2.0296e-2 off (≈ 1 580 bands, `…140347Z…-W1-cmode.log:2102`). A 1e-4 error in C moves f₀ by ≈ 5e-5, ≈ 3.9 bands.
- **What it does not separate** (stated so no one reads more into it): C_mode's in-model zero (4.821e-6, item 217) also lies inside it. This assertion licenses "C_in is the in-model zero to the grid's resolution". It does not rank C_in over C_mode; the ruling already did that.
- **T22′ needs its own band:** item 228 W4 `:2054–2056` gives a max of 1.616e-5 (×7.1 \|r\|). It is not used here.

**Code** (`tests/validation/test_ans7_exact_strip.py`; with the new knob unset, every existing path is byte-identical).
1. **Knob `FEM_EM_ANS7_EXACT_TUNE_C`** ∈ {unset, `cin`}; anything else raises `ValueError`. `cin` requires `FEM_EM_ANS7_EXACT_TUNE=inmodel` (`_tune()`, `:977`) with `FEM_EM_ANS7_EXACT_TUNE_MHZ` unset (64 MHz only; T22′ extends it), and raises `ValueError` otherwise.
2. **Constant** `T21P_ZERO_RESOLUTION = 1.287e-5`, with a comment citing the three log lines above.
3. **`tune_inmodel` (`:1369`) under `cin`:**
   - read `c_in` from `tune_p2.json` (never recomputed);
   - (W) first, unchanged (`:1398–1413`);
   - the existing three-point loop (`:1416–1437`) with the capacitor swapped to `c_in`: `_route_a_solve(…, c_f=c_in, degree=deg)`, Route B comparand `route_b_z(bary_interp(xn, z_nodes_cu, g − 64.5), f, c_in)`, Z_in,A = `add_c1(_z_from_s11(s_a), ω, c_in)`; rows from `parallel_rows`;
   - print c = C_mode/C_in − 1, f₀,B(C_in) from `sparse_rows(…, c_f=c_in)` on the stored nodes, and r = f₀,A/f₀,B(C_in) − 1;
   - (N2) at 0.5·C_in; (N3′) below;
   - rank 0 writes `tune_inmodel_cin_p2.json`. Never overwrite item 217's `tune_inmodel_p2.json`.
4. **Tests** (skip unless `TUNE_C=cin`): `test_t21p_ab_circuit_equals_field_cin`, `test_t21p_inmodel_zero_cin`, `test_t21p_n1_untuned`. The existing `test_t21_ab_circuit_equals_field_p2` / `test_t21_t2_inmodel_zero` skip under `cin`, because they assert at C_mode. Their bodies and `BANDS` are untouched.

**Anchors** (asserted; bands imported by name).
- **(AB′) circuit = field at C_in,** degree 2, Cu: \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` (1e-3, `tests/validation/test_port_lumped_rlc_termination.py:88`, imported through `test_port_circuit_layer_field`) at all three `TUNE_MHZ` points (`test_ans7_exact_strip.py:925`).
- **(T2′) the in-model zero at C_in:** exactly one − → + crossing of Im `Y_in`,A in (63.9, 64.1) MHz, and \|f₀,A/`TUNE_TARGET_HZ` − 1\| ≤ `T21P_ZERO_RESOLUTION` (target at `:924`).
- **(W)** unchanged: \|Δresidual\| ≤ `NODE_S_RTOL` against item 195's stored 64.0 MHz point.

**Negative controls (rule (e)).**
- **(W): asserted,** backed by item 217 W2 `:2045` (2.772e-14; same fixture, same comparison).
- **(N1′) untuned zero: asserted,** pure, on the stored nodes: \|f₀,B(C_F)/64 MHz − 1\| ≥ 10 × `T21P_ZERO_RESOLUTION`. Backed by item 217 W1 `:2102` (2.029559e-2, ≈ 1 580×). No ceiling binds: an O(1) relative offset is reachable.
- **(N2′): predicted, printed.** The circuit at 0.5·C_in against the in-model at C_in, 64.0 MHz, is ≥ 10× `REDUCTION_BAND`. T21 measured 493.78× at C_mode (`…140738Z…:2051`); the ceiling is ≲ 2 000×.
- **(N3′): predicted, printed.** Route B at C_in·(1 + 1e-4) moves f₀ by ≈ −5e-5, ≈ 3.9 × `T21P_ZERO_RESOLUTION`. This shows the band resolves a 1e-4 capacitor error.

**Predictions** (written first, printed).
- f₀,B(C_in)/64 MHz − 1 within ±1e-9. C_in zeroes the circuit at 64.0 MHz by construction; interpolation error is ≈ 1e-10.
- r in [−4e-6, −1.5e-6] (flat in f: −2.531e-6 / −2.535e-6 / −2.263e-6). Hence \|f₀,A(C_in)/64 MHz − 1\| ≈ 2.5e-6, **predicted green** against 1.287e-5.
- (AB′) in [1.0e-5, 1.5e-5] at each point. C moved by 1.5e-5 relative; T21 read 1.237–1.287e-5.
- c = C_mode/C_in − 1 = −1.494697e-5, read back from the JSON (item 217 W1 `:2089`).
- R_in and Q at C_in within 0.1 % of T21's Route A rows (18 739.36 Ω / 972.90, `…140738Z…:2049`).

**Window** (`-n 8`; bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE` unset).
- **Price:** item 217 W2 is the same window with C swapped: 261 s, drives 41.8–50.3 s, (W) 5.0 s (`…140738Z…:2045–2048`; test-results.md row `ANS-7-item217-W2-inmodel`). A ×1.4 slow box gives ≈ 365 s, inside 590.
- **W1** (closing window):

  `scripts/testing/run_and_log.sh ANS-7-item262-W1-cin "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- There is no separate default-suite window. With the knob unset, the paths are byte-identical and the new tests skip; `git diff` shows additions plus the skip guards only.
- **No edit after W1** (rule (i)), so there is no W2.

**Traps already paid for.**
- Items 217 / 228's list applies whole:
  - κ `None`;
  - no `ans:7` `_stage_*` / `_assemble` / `_write_outputs` call, and `ans:7` is not run;
  - `z_from_s(S)`, never the sweep's `z_matrix`;
  - real-scalar secant, and no ordering comparison on complex operands;
  - rows computed on every rank, rank 0 writes, then a barrier;
  - a killed window is not a result: `pgrep -c python3`, report, never lengthen the timeout;
  - bare, whole harness calls, no `-k` selection.
- Never write `tune_p2.json`, `tune_inmodel_p2.json`, `t127.74/*` or item 196's node files. This item only reads them.
- Don't route through `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ` (item 228's trap: `_nodes_env` / `p2_asm` assemble at C_F).
- The band is a literal. Computing it from this run's (AB′) would make the assertion self-referential.

**Scope.** Covered: the ×1 exact-sheet mesh (h = 0.015), degree 2, κ off, the Cu column, 64.000 MHz, Route A at three points with C_in, and Route B on item 196's stored grid. **Not claimed:**
- no absolute tuned f₀, Q, R_in or C. T6 / T13 are open and the h-rung T12′ has not landed; C_in is a ×1-mesh, degree-2 number whose h-move is unmeasured;
- no tuning-accuracy claim against any external tool;
- not 127.74 MHz (T22′);
- not the 51-point Route A, and not the eigen route;
- no SPEC change (66.29 pF stays);
- no AED or Palace comparison.

**Status it can move.**
- On a green (T2′), retires known-issues **2026-10-04** (`ANS-7` step 7, §9 item 217, chain T21; `docs/testing/known-issues.md:1353`), as the ruling states. The C_mode readings (4.821e-6 / 1.390e-6) stay records, superseded by the ruling, and `BANDS` is unchanged.
- The `ANS-7` row gains a "step 7′" record sentence: C_in, the in-model zero, c, r and (AB′).
- **No glyph.** The §7 row's Done-when is the weekly's T13 re-read after T12′ (§10 T13 / T24). This step does not reach it, so `ANS-7` stays 🟡 per §2.
- T22′ (127.74 MHz) and T23 (`ans:7` tuned rows) become queueable.
- **Rule (j):** the review adds a step-7′ Done-when sentence.

**Negative result.**
- **(W) or (N1′) red:** BLOCKED (wiring); revert the knob.
- **(AB′) red with (W) green:** circuit and field disagree at C_in, where they agreed at C_mode 1.5e-5 away. Print the three residuals, open a known-issues entry, park, BLOCKED. Never apply κ.
- **(T2′) with no crossing in (63.9, 64.1) MHz:** BLOCKED; report the Im `Y_in` triple.
- **(T2′) band half red:** \|r\| at C_in has grown ≥ ×5 against its flat 62.7 / 64 / 127.74 MHz records. That contradicts the ruling's premise, so it is reported, never re-banded. Park BLOCKED with r and c printed; the weekly decides.
- No `PREDICTED RED:` line. The prediction is green with a ×5 margin, so a red here is a finding.

**Outcome (2026-10-06 06:00 CDT slot): GREEN.** W1 `20261006T110217Z_ANS-7-item262-W1-cin.log`, Status 0, 257 s (pytest 255 s; W2-part 161.1 s: mesh 23.5 s, (W) 2.9 s, drives 48.2 / 42.5 / 43.9 s), `-n 8`, 80 798 cells; 16 passed, 19 skipped. Module mtime 11:02:02Z precedes the log start 11:02:17Z (rule (i)). C_in = 6.359009101426113e-11 F read from `tune_p2.json`. (W) \|Δ\| 4.034e-14 (≤ `NODE_S_RTOL`). (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5 at 63.9 / 64.0 / 64.1 MHz (≤ `REDUCTION_BAND`; predicted [1.0e-5, 1.5e-5], met). (T2′) Im `Y_in`,A = −1.616747e-4 / +2.628149e-7 / +1.625093e-4 S, one − → + crossing at 63 999 837.71 Hz: f₀,A/64 MHz − 1 = −2.535844e-6 = 0.197 × `T21P_ZERO_RESOLUTION` (predicted ≈ 2.5e-6, met). f₀,B(C_in)/64 MHz − 1 = +1.88e-12 (predicted ±1e-9, met). r = −2.535846e-6 (predicted [−4e-6, −1.5e-6], met). c = −1.494697e-5 (met). (N1′) 2.029559e-2 = 1 577 × band (asserted ≥ 10×). (N2′) 4.938e-1 = 493.78 × `REDUCTION_BAND` (predicted ≥ 10×, met). (N3′) −4.928e-5 = 3.83 bands (predicted ≈ −5e-5 / 3.9×, met). R_in 18 739.27 Ω (−0.0005 % vs T21 Route A, met); **Q 971.07, −0.1885 % vs T21's 972.90 — the printed 0.1 % prediction is missed** (printed only; a 1.5e-5 move in C moving Q by 1.9e-3 is not explained here). Known-issues 2026-10-04 retired. `ANS-7` stays 🟡.

**Review note (2026-10-06 10:30 daily review).** No §4 audit applies: there is no glyph move. The brief's "rule (j): the review adds a step-7′ Done-when sentence" is **declined**, as item 228's was: the step moves no glyph, and the row's Done-when (the weekly's T13 re-read) is unchanged. The slot's record sentence in the row stands. T22′ (item 283) is now queued. T23 waits on T22′, because its rows need both fields.

### §9 item 283 — 2026-10-06 (chain T22′)

*Written by the 2026-10-06 10:30 daily review.* Licence: §10 chain T22′ (grep "T22′. T21′ at 127.74 MHz"), quoted, not changed: "T21′ at 127.74 MHz. Serial on T21′." T21′ landed green (item 262, `9bb7625c`). Item 262's construction, target and ruling carry verbatim with 64.000 → 127.74 MHz: C_in, the in-model zero within the sparse grid's own measured resolution, pre-registered and > r, and circuit = field at `REDUCTION_BAND`. The band was already named by item 262's brief ("T22′ needs its own band: item 228 W4 `:2054–2056` gives a max of 1.616e-5"). This brief registers it; it does not choose it. Exact rectangles, κ off. `BANDS` is not edited. **No band, tolerance, Done-when or target is set or moved.**

`ANS-7` step 8′: drive the exact-sheet F-small birdcage at degree 2 on the ×1 mesh with item 228's **C_in at 127.74 MHz** (1.5274843607471887e-11 F, `logs/ans7_exact/t127.74/tune_p2.json`; `20261004T203817Z_ANS-7-item228-W3-cmode.log:2097`), at T22's three points 127.64 / 127.74 / 127.84 MHz, and read the in-model zero.
- **Executor:** implementer. Tests only (`tests/validation/test_ans7_exact_strip.py`). No `src/`, `ans:7`, `TH-14` or `test_ans7_sparse_route_b.py` edit.
- **Tier:** heavy, `-n 8`, complex. **Class:** chain step (solve), ≈ 19 measured / ≈ 40 predicted min (15 fixed + ≈ 10 code + W1 ≈ 4.5 + W2 ≈ 5 + ≈ 5 spare).
- **Depends on:** nothing open. The seven `t127.74/node_*_p2_cu.json`, `t127.74/tune_p2.json` and `t127.74/tune_inmodel_p2.json` are on disk (checked 2026-10-06 10:30). If `t127.74/tune_p2.json` is absent, stop and take the next item. Never re-create it.

**The resolution band, pre-registered** (measured, not chosen). `T22P_ZERO_RESOLUTION = 1.616e-5` is the maximum, at logged precision, of T22's (AB) node residuals at C_mode: 1.589e-5 / 1.601e-5 / 1.616e-5 at 127.64 / 127.74 / 127.84 MHz (`20261004T204206Z_ANS-7-item228-W4-inmodel.log:2054–2056`). It is a literal with that citation, written before W1 and never recomputed. It exceeds r = −2.263e-6 (`:2058`) by ×7.1. The untuned grid has no zero in the window at all ((N0′) below), and a 1e-4 C error moves f₀ by ≈ 4.7e-5 ≈ 2.9 bands (sensitivity 4.711572e-6 per 1e-5, `…203817Z…:2100`), so the band resolves a capacitor error. As at 64 MHz, it does not rank C_in over C_mode (T22's C_mode zero, 1.390e-6, also lies inside it); the ruling did that.

**Code** (`tests/validation/test_ans7_exact_strip.py`; every knob-unset path is byte-identical, and item 262's `cin` path at 64 MHz is behaviourally identical).
1. **`cin` accepts the 127.74 set.** `FEM_EM_ANS7_EXACT_TUNE_C=cin` currently raises unless `FEM_EM_ANS7_EXACT_TUNE_MHZ` is unset (`:1011`). Allow `TUNE_MHZ` ∈ {unset, `127.74`}. Anything else still raises `ValueError`.
2. **Constant** `T22P_ZERO_RESOLUTION = 1.616e-5`, with the citation above. The active band is `T21P_ZERO_RESOLUTION` when the knob is unset and `T22P_ZERO_RESOLUTION` at 127.74, chosen by the knob and never by a reading.
3. **The `cin` path reads the active set** (item 228's active-set helper: target, tune points, nodes, grid, out dir):
   - C_in comes from the active `tune_p2.json`;
   - the three Route A drives run at the active tune points;
   - Route B comes from the active stored Cu nodes (`COLUMNS=cu`) on the active grid (`_loo` / `sparse_rows` with `grid_mhz=`);
   - the (T2″) bracket is (tune points[0], tune points[−1]);
   - rank 0 writes `t127.74/tune_inmodel_cin_p2.json`. It never overwrites item 228's `t127.74/tune_inmodel_p2.json` or item 262's `tune_inmodel_cin_p2.json`.
4. **Controls by target.** At 64 MHz, (N1′) is unchanged. At 127.74, (N0′) replaces it: on the stored 127.74 MHz Cu grid at C_F, the up-crossings of Im `Y_in` are `[]`. It is pure and takes milliseconds.
5. Prints carry the target in their tag (`T22′`), so the two targets' logs never read alike.

**Anchors** (asserted; bands imported by name), at 127.74 MHz:
- **(AB″) circuit = field at C_in:** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` (`test_port_lumped_rlc_termination.py:88`) at all three points.
- **(T2″) the in-model zero:** exactly one − → + crossing of Im `Y_in`,A in (127.64, 127.84) MHz, and \|f₀,A/127.74 MHz − 1\| ≤ `T22P_ZERO_RESOLUTION`.
- **(W):** unchanged, the degree-1 default drive at 64.0 MHz against item 195's stored residual, ≤ `NODE_S_RTOL`.

**Negative controls (rule (e)).**
- **(W): asserted,** backed by item 262 W1 (4.034e-14, same fixture, same comparison).
- **(N0′): asserted,** backed by the same comparison on the same grid: up-crossings `[]` at C_F (`20261004T203817Z_ANS-7-item228-W3-cmode.log:2104`). This is binary, so there is no ceiling.
- **(N2″): predicted, printed.** The circuit at 0.5·C_in against the in-model at C_in is ≥ 10× `REDUCTION_BAND`. T22 read 493.82× at C_mode (`…204206Z…:2060`).
- **(N3″): predicted, printed.** Route B at C_in·(1 + 1e-4) moves f₀ by ≈ −4.7e-5, ≈ 2.9 bands.

**Predictions** (written first, printed).
- f₀,B(C_in)/127.74 MHz − 1 within ±1e-9.
- r in [−4e-6, −1.5e-6], hence \|f₀,A(C_in)/127.74 MHz − 1\| ≈ 2.26e-6 ≈ 0.14 bands: **predicted green**. Item 228's one-drive estimate at C_in read δf/f = −2.262515e-6 (`…204206Z…:2061`).
- (AB″) in [1.4e-5, 1.8e-5] at each point.
- c = C_mode/C_in − 1 = −7.754048e-6, read back from the JSON.
- R_in and Q at C_in within 0.3 % of T22's Route B rows at C_mode (23 688.41 Ω / 616.17, `…203817Z…:2099`). Item 262 missed its 0.1 % Q prediction by −0.19 % for a 1.5e-5 C move, so this one is wider.

**Windows** (`-n 8`; bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE` unset). No edit after W1 starts (rule (i)), so W1 and W2 run the same bytes.
- **W1, the 64 MHz re-gate of item 262's path** (the edit touches its gate): item 262's W1 command verbatim, label `ANS-7-item283-W1-cin64`. It must reproduce (T2′) −2.535844e-6 and (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5 to the printed digits. Price 257 s (item 262). **Red ⇒ revert the edit, BLOCKED.**
- **W2, 127.74 MHz** (closing window), label `ANS-7-item283-W2-cin128`. Price ≈ 260–300 s: item 228 W4 ran (W), three C_mode drives and one C_in drive in 297 s.

  `scripts/testing/run_and_log.sh ANS-7-item283-W2-cin128 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_C=cin timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`

**Traps already paid for.**
- Items 217 / 228 / 262's lists, whole:
  - κ `None`;
  - no `ans:7` call;
  - `z_from_s(S)`;
  - a real-scalar secant;
  - rank 0 writes, then a barrier;
  - bare, whole harness calls;
  - a killed window is not a result.
- **`COLUMNS=cu` is required at 127.74.** The stored nodes are `_cu`-suffixed, and without it the loader looks for files that do not exist.
- **Never route through `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ`** (item 228's trap: `_nodes_env` validates against the 64 MHz nodes).
- **A missed 64 MHz literal in the `cin` path** shows up as a W2 bracket or `f_hz` assert failing on a correct measurement (item 228's code step 4).
- **The band is a literal.** Never compute it from this run.

**Scope.** Covered: the ×1 exact-sheet mesh, degree 2, κ off, the Cu column, 127.74 MHz, Route A at three points with C_in, and Route B on item 228's stored Cu grid. **Not claimed:** any absolute tuned f₀, Q, R_in or C at 3 T, or any h-statement (T12′ is 64 MHz only). Also not: the PEC grid at 3 T, a SPEC change, or an AED or Palace comparison (none exists at 3 T).

**Status it can move.**
- The `ANS-7` row gains a "step 8′" record (narrative here; one clause in the row): C_in at 127.74 MHz, the in-model zero, c, r and (AB″).
- **No glyph:** the Done-when is the weekly's T13.
- T23 (`ans:7`'s tuned-to-Larmor rows at both fields) becomes queueable.

**Negative result.**
- **W1 red:** the edit changed item 262's path. Revert it and mark BLOCKED.
- **(W) or (N0′) red:** BLOCKED (wiring).
- **(AB″) red with (W) green:** open a known-issues entry, park and mark BLOCKED. Never apply κ.
- **No crossing in (127.64, 127.84):** BLOCKED. Report the Im `Y_in` triple.
- **(T2″) band half red:** \|r\| has grown ≥ ×7 against its flat record, which contradicts the ruling's premise. Report r and c, park BLOCKED, and the weekly decides.
- There is no `PREDICTED RED:` line: the prediction is green with a ×7 margin.

**Outcome (2026-10-06 12:00 CDT slot, second item): GREEN.** `-n 8`, 80 798 cells, both windows 16 passed, 19 skipped, Status 0. Module mtime 17:13:04Z precedes W1's log start 17:13:12Z (rule (i)); W1 and W2 ran the same bytes. **W1** `20261006T171312Z_ANS-7-item283-W1-cin64.log`, 258 s (W2-part 164.1 s): reproduces item 262 to the printed digits — (T2′) −2.535844e-6 (0.197 bands; Im `Y_in`,A −1.616747e-4 / +2.628149e-7 / +1.625093e-4 S, crossing 63 999 837.71 Hz), (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5, r −2.535846e-6, c −1.494697e-5, (N1′) 1 577×, (N2′) 493.78×, (N3′) −4.928048e-5, R_in 18 739.27 Ω, Q 971.07; (W) 1.127e-14. **W2** `20261006T171745Z_ANS-7-item283-W2-cin128.log`, 261 s (W2-part 165.9 s: mesh 23.5 s, (W) 2.8 s, drives 45.2 / 46.6 / 47.9 s). C_in = 1.5274843607471887e-11 F read from `t127.74/tune_p2.json`. (W) 1.470e-13 (≤ `NODE_S_RTOL`). (AB″) 1.590e-5 / 1.601e-5 / 1.616e-5 at 127.64 / 127.74 / 127.84 MHz (≤ `REDUCTION_BAND`; predicted [1.4e-5, 1.8e-5], met). (T2″) Im `Y_in`,A = −4.057100e-5 / +1.177054e-7 / +4.084487e-5 S, one − → + crossing at 127 739 710.72 Hz: f₀,A/127.74 MHz − 1 = −2.264621e-6 = 0.140 × `T22P_ZERO_RESOLUTION` (predicted ≈ 2.26e-6, met). f₀,B(C_in)/127.74 MHz − 1 = −9.16e-13 (predicted ±1e-9, met). r = −2.264620e-6 (predicted [−4e-6, −1.5e-6], met). c = −7.754048e-6 (met). (N0′) up-crossings `[]` at C_F (asserted). (N2″) 4.938232e-1 = 493.82 × `REDUCTION_BAND` (predicted ≥ 10×, met). (N3″) −4.714727e-5 = 2.92 bands (predicted ≈ −4.7e-5 / 2.9×, met). R_in 23 793.40 Ω, Q 618.34: **the printed 0.3 % prediction against T22's Route B rows at C_mode is missed** (+0.443 % / +0.352 %; printed only). The comparand mixes route and capacitor: item 230 measured the Route A / Route B R_in and Q gap at +0.37 % / +0.38 % on this fixture, which is the size of the miss, so the C move itself is not resolved from it. Under `cin` at 127.74 item 228's one-drive C_in record is not repeated (the three drives already run at C_in). `ANS-7` stays 🟡; T23 is queueable.

### §9 item 288 — 2026-10-06 (chain T23)

*Drafted by a general-purpose agent for the 2026-10-06 18:00 daily review; the review reviews it.* Licence: §10 chain T23 (grep "T23. `ans:7` gains the tuned-to-Larmor rows"), quoted, not changed: "`ans:7` gains the tuned-to-Larmor rows (T21′ / T22′: C_in, the in-model zero, c, r at both fields); the degree-1 rows stay labelled superseded. `example-runner`. Serial on T21′." T21′ landed as item 262 (`9bb7625c`), T22′ as item 283 (`43f836bc`). **No band, tolerance, Done-when, SPEC or validation target is set or moved**; every band is imported by name.

**Housing.** No EX chunk owns `ans:7`: the §7 EX table has no `ans:7` row, and `examples/ansys_benchmarks/README.md` lists the case under **ANS-7**. The precedent is T20 (item 230), which re-issued `ans:7` as `ANS-7` step 9 with `example-runner`. T23 is therefore **`ANS-7` step 10**, briefed here.

`ANS-7` step 10: `ans:7` gains two additive stages, `XT64` and `XT128`. Each re-drives Route A at degree 2 on the ×1 exact-sheet mesh at **C_in** at its field's three tune points. It re-asserts T21′ / T22′'s gates through the example and writes the tuned-to-Larmor rows into `metrics.json` (`degree2_exact["tuned_larmor"]`) and `COMPARISON.md`. The degree-1 keys stay labelled superseded, and (K) proves they are untouched.
- **Executor:** `example-runner`. Spawn prompt: "Use the example-runner agent to execute `ANS-7` step 10 (§9 item 288) per `docs/planning/chunks/ANS-7.md` § '§9 item 288'. Run every harness window foreground with Bash `timeout: 660000`. Never background a window or end your turn while one is running. Run `XT64` and `XT128` serially, never concurrently. Set `FEM_EM_ANS_SKIP_PRIVATE=1` on every window. Never run `ans:1`–`ans:3`, never run B1 / B2 / A / X1–X4 / XA, and never use `./scripts/run_examples.sh -e ans`."
- **Edits:** the example script, its same-stem guide, `metrics.json` and `COMPARISON.md` only. No `src/`, test-module, SPEC or `aed/` edit.
- **Tier:** heavy, `-n 8` (T21′ / T22′'s record width), complex build.
- **Class:** example. **Time:** ≈ 12 measured / ≈ 50 predicted min (15 fixed + ≈ 15 code and guide + `XT64` ≈ 4.5 + `XT128` ≈ 4.5 + ≈ 3 census / pin + ≈ 8 spare).
- **Depends on:** items 262 and 283, both landed. Stored inputs, all on disk (checked 2026-10-06):
  - `logs/ans7_exact/tune_p2.json`, `tune_inmodel_p2.json`, `tune_inmodel_cin_p2.json` and item 196's seven `node_<mhz>_p2.json`;
  - `logs/ans7_exact/t127.74/tune_p2.json`, `tune_inmodel_p2.json`, `tune_inmodel_cin_p2.json` and the seven `node_<mhz>_p2_cu.json`.
  - If any is absent, stop and take the next item. Never re-create one; that is the test module's job.

**Code** (`examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py`; additive; B1 / B2 / A / X1–X4 / XA are byte-unchanged except for the `main()` stage list).
1. **Stages `XT64`, `XT128`** are added to `main()`'s accepted set. They start with `X`, so `_x_guard()` applies: every `FEM_EM_ANS7_EXACT_*` / `FEM_EM_ANS7_SPARSE_*` knob must be unset. That is why the stage **never calls the knob-driven helpers** `_active()`, `_stored_active()`, `_tune_mhz()` or `_zero_resolution()`. It selects the field set by stage name alone.
2. **Lazy imports inside the stage function** (item 230's cycle trap; never module-level). Imported, never restated:
   - from `tests.validation.test_ans7_exact_strip`: `TUNE_TARGET_HZ`, `TUNE_MHZ`, `TUNE128_TARGET_HZ`, `TUNE128_MHZ`, `TUNE128_NODES_MHZ`, `TUNE128_GRID_MHZ`, `OUT_DIR`, `TUNE128_OUT_DIR`, `T21P_ZERO_RESOLUTION`, `T22P_ZERO_RESOLUTION`, `_node128_json`, `_stored_x1_p2`, `_tie_from_json`, `_assert_tie`, `_fixture_tie`, `_exact_geo`, `_uj`;
   - from `tests.validation.test_ans7_sparse_route_b`: `DEFAULT_NODES_MHZ`, `GRID_MHZ`, `NODE_S_RTOL`, `sparse_rows`, `bary_interp`;
   - `REDUCTION_BAND` (`test_port_circuit_layer_field`), `GATED_WIDTH_FRACTION`.
3. **Inputs per stage.**
   - `XT64`: `OUT_DIR / "tune_p2.json"`. Route B nodes come from `_stored_x1_p2()` (Cu Z), which are T21′'s exact inputs. Before calling it, check that each file exists and raise `RuntimeError` naming any that is missing.
   - `XT128`: `TUNE128_OUT_DIR / "tune_p2.json"`. Route B nodes come from `_node128_json(m, ("cu",))` over `TUNE128_NODES_MHZ`. Check them inline as `_stored_active` does: degree 2, κ `None`, `target_mhz == 127.74`, one cell count, and `_assert_tie` per file.
   - From `tune_p2.json`, read `c_in`, `c_mode["cu"]`, `degree == 2` and `abs(f_hz − target) ≤ 1 Hz`. Never recompute C_in.
4. **Drives** follow item 262's loop verbatim. Geometry is `_exact_geo(th14, comm, GATED_WIDTH_FRACTION)`. At each tune point g:
   - `s_a = _route_a_solve(th14, geo, specs, g·1e6, "cu", c_f=c_in, degree=2)`;
   - `s_b = _s11_from_z(route_b_z(bary_interp(xn, z_cu, g − 64.5), f, c_in)[1])`, with `xn = nodes − 64.5` exactly as the test does at both fields;
   - Z_in,A = `add_c1(_z_from_s11(s_a), ω, c_in)`;
   - rows from `parallel_rows`;
   - f₀,B(C_in) from `sparse_rows(nodes, {"cu": z_cu}, grid_mhz=grid, columns=("cu",), c_f=c_in)`;
   - c = C_mode/C_in − 1 and r = f₀,A/f₀,B(C_in) − 1.
   - There is no (W) drive: it is the test module's wiring check, and the reproduction anchor (R-T) below replaces it here.
5. **Write** (rank 0, then a barrier).
   - Snapshot every top-level key and every `degree2_exact` key except `tuned_larmor`, plus the other field's entry if it is present.
   - Set `degree2_exact["tuned_larmor"]["64.0" | "127.74"]` = {`target_hz`, `c_in_f`, `c_mode_f`, `c_rel`, `tune_mhz`, `s11_a`, `s11_b`, `residual_ab`, `reduction_band`, `rows_a`, `f0a_rel_target`, `zero_resolution` (value and constant name), `bands_fraction`, `f0b_cin_hz`, `r`, `n1_off` or `n0_crossings_up_hz`, `n2`, `rt` (the reproduction residuals), `cells`, `mpi_ranks`, `t_mesh_s`, `t_drive_s`, `source` (the item 262 / 283 log names)}.
   - Assert (K), and revert `metrics.json` on a miss, as XA does.
   - Regenerate `COMPARISON.md`. `_write_comparison` gains a third section, "Tuned to the Larmor frequency: C_in, the in-model zero (×1 mesh, degree 2; §10 T23)". It has one column per field (64.000 / 127.74 MHz) and rows C_in [pF], f₀,A [MHz], f₀,A/target − 1, bands used, c, r, (AB) max, R_in, Q. A field not yet run shows "—". **There is no AED column:** no AED tuned case exists, and the SPEC's C is 66.29 pF. The existing two tables are byte-unchanged.
6. **Docs.**
   - The docstring gets an `XT64` / `XT128` paragraph.
   - The guide gets a section after "Degree 2 on exact rectangles (current; §10 T20)", titled "Tuned to the Larmor frequency (current; §10 T23)": the two stages, the commands and the anchors, in a table like the X-stage one.
   - "Not claimed" gains: "C_in, f₀, R_in, Q at C_in are ×1-mesh degree-2 records — no absolute tuned C or tuning-accuracy claim (T12′ / T13)".
   - **Setup figure and combined XDMF (`EX-93` conventions):** they do not apply to this modified existing example. `XT*` calls neither `write_setup_figure` nor `_write_xdmf` (the drives run `return_fields=False`). The tracked PNG and stage A's degree-1 Cu/PEC XDMF stand, so the census is predicted unchanged. Do not set `FEM_EM_SETUP_FIGURES`.

**Anchors** (asserted in-stage; every band imported):
- **(AB′) / (AB″):** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` at all three tune points.
- **(T2′) / (T2″):** exactly one − → + crossing of Im `Y_in`,A inside (tune[0], tune[−1]), and \|f₀,A/target − 1\| ≤ `T21P_ZERO_RESOLUTION` (XT64) or `T22P_ZERO_RESOLUTION` (XT128). The stage chooses the band, never a reading.
- **(R-T) reproduction of the registered run:** against the stored `tune_inmodel_cin_p2.json` of its field:
  - \|f₀,A,ex/f₀,A,stored − 1\| ≤ `NODE_S_RTOL`;
  - max \|res_ex − res_stored\| ≤ `NODE_S_RTOL` (absolute, as (W));
  - max \|S11_A,ex − S11_A,stored\|/\|S11_A,stored\| ≤ `NODE_S_RTOL`;
  - also asserted: the stored file's `ranks == 8`, `degree == 2` and `c_in` equals `tune_p2.json`'s bit for bit.
  - At any width other than 8 it prints and does not assert (the `EX-56` rule).
- **(K)** as in code step 5.

**Negative controls (rule (e)).**
- **(N1′) asserted, XT64:** \|f₀,B(C_F)/64 MHz − 1\| ≥ 10 × `T21P_ZERO_RESOLUTION` on the stored nodes. Backing 2.029559e-2 = 1 577× (item 262). This is pure: milliseconds, no ceiling binds.
- **(N0′) asserted, XT128:** at C_F the stored 127.74 MHz Cu grid has up-crossings `[]` (item 283 W2; item 228 W3 `:2104`). A missing key reads `None` and fails, never a vacuous `[]`.
- **(N-C) asserted, both stages:** the (R-T) comparator separates C_in from C_mode. \|f₀,A,ex(C_in)/f₀,A,stored(C_mode) − 1\| > 100 × `NODE_S_RTOL`, where the C_mode reading is `rows_a.f0_hz` of the field's stored `tune_inmodel_p2.json` (items 217 / 228). Backing: the C_mode zeros sit at 4.821e-6 / 1.390e-6 against C_in's −2.536e-6 / −2.265e-6. Whatever the sign, the separation is ≥ 8.75e-7, i.e. ≥ 875× `NODE_S_RTOL`.
- **(N2′) / (N2″) predicted, printed:** the circuit at 0.5·C_in against the in-model at C_in at the target is ≥ 10 × `REDUCTION_BAND` (measured 493.78× / 493.82×).

**Predictions** (written first, printed):
- (R-T) ≤ 1e-12. Item 283 W1 reproduced item 262 to every printed digit at `-n 8`.
- (T2′) −2.535844e-6 = 0.197 bands; (T2″) −2.264621e-6 = 0.140 bands.
- (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5; (AB″) 1.590e-5 / 1.601e-5 / 1.616e-5.
- c −1.494697e-5 / −7.754048e-6.
- r −2.535846e-6 / −2.264620e-6.
- R_in / Q at C_in: 18 739.27 Ω / 971.07 and 23 793.40 Ω / 618.34. These are records, printed.

**Price.**
- Item 262 W1: 257 s, W2-part 161 s (mesh 23.5 s, (W) 2.9 s, drives 48.2 / 42.5 / 43.9 s).
- Item 283 W2: 261 s, W2-part 166 s.
- With no pytest, no `tests/environment` and no (W): `XT64` ≈ 180–260 s and `XT128` ≈ 180–265 s. Each is under 590 s on a ×1.4 slow box (≈ 370 s).

**Windows** (`-n 8`; each bare, whole and foreground; Bash `timeout: 660000`; the two never concurrent; no edit after `XT64` starts, rule (i)).
0. Orphan check: `docker compose exec -T fem-em-solver pgrep -c python3` must read 0.
1. **Censuses before**, labels `ANS-7-item288-pre-figures` / `-pre-docrefs`:
   - `scripts/testing/run_and_log.sh ANS-7-item288-pre-figures "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py'"`
   - `scripts/testing/run_and_log.sh ANS-7-item288-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
   - Last reading (item 273 post): figures `examples` 88 with `broken=0`; docrefs runnable 88, `dead=0`, `guide=0`, exit 2.
2. **XT64:**

   `scripts/testing/run_and_log.sh ANS-7-item288-XT64 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_STAGE=XT64 timeout -k 30 590 mpiexec -n 8 python3 examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py'"`
3. **XT128:** the same command with `FEM_EM_ANS7_STAGE=XT128`, label `ANS-7-item288-XT128`. Run it only after XT64 has returned its footer.
4. **Pin and unit tests** on the committed `metrics.json`, at `-n 1`, labels `ANS-7-item288-pin` / `-unit` (item 230's pattern; the path is already in `LISTED_EXAMPLE_ARTIFACTS`):

   `scripts/testing/run_and_log.sh ANS-7-item288-pin "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/unit/test_doc_reference_exit_codes.py tests/unit/test_ans7_private_comparison.py -v -s --tb=short'"`

   Predicted 20 + 4 passed.
5. **Censuses after**, labels `-post-figures` / `-post-docrefs`. Predicted equal to the before readings: figures +0 with `broken=0`; docrefs runnable +0, `dead=0`, `guide=0`, `exit != 1`.
6. Commit the script, guide, `metrics.json`, `COMPARISON.md`, logs and test-results rows together. `git diff --stat` must show no `COMPARISON_private.md`, `aed_results/`, `SPEC.md` or test-module change.

**Traps already paid for.**
- Item 230's list, whole:
  - lazy imports;
  - `_write_outputs` rewrites wholesale, so never run B2 here;
  - the degree-1 keys have four readers, so never rename or move them;
  - `z_from_s(S)`, never `z_matrix`;
  - κ `None` on every spec;
  - AED values never enter a tracked file;
  - a killed window is not a result: run `pgrep -c python3`, report, never lengthen the timeout;
  - bare, whole harness calls with no pipe.
- **Knob-free selection:** setting `FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74` to reach the 127.74 set trips `_x_guard`. Choose by stage name only.
- **`ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ` are never routed** (item 228's trap: the 64 MHz node validator).
- **The 127.74 nodes are `_cu`-suffixed:** pass `("cu",)` to `_node128_json`.
- **A missed 64 MHz literal in the 127.74 path** shows as a bracket or `f_hz` failure on a correct measurement (item 228 code step 4).
- **Docrefs exit 1 on a bare file name in the guide** (item 273's first post-docrefs). Name the stored inputs as "the gitignored `logs/ans7_exact/` stage files", as the existing section does for `logs/ans7_p2/`.
- Write ≈, never `~`.

**Scope.** Covered: the example re-asserting T21′ / T22′ at C_in on the ×1 exact-sheet mesh, degree 2, κ off, Cu, three points per field, and recording the rows in its tracked outputs. **Not claimed:**
- any absolute tuned C, f₀, R_in or Q;
- any h-statement (T12′);
- the PEC column at C_in;
- the 51-point Route A;
- a SPEC change;
- any AED or Palace comparison;
- the private writer's rows.

**Status it can move.**
- §10 T23 closes, and the weekly removes it from the live list.
- The `ANS-7` row gains a "step 10 (T23)" clause: `ans:7`'s tracked outputs carry the tuned-to-Larmor rows at both fields, with (T2′) / (T2″), (AB′) / (AB″) and (R-T) re-asserted through the example.
- **No glyph:** the row's Done-when is the weekly's T13 re-read after T12′ (read 2026-10-06), so `ANS-7` stays 🟡.
- **Rule (j):** declined, as for items 228 / 262 / 283. No glyph moves, so no step Done-when is added.

**Negative result.**
- **(R-T) red at `-n 8` with the bands green:** the example's path is not the registered run's. Print the differences, revert `metrics.json` / `COMPARISON.md` (`git checkout` those two paths only), park on `attempt/ANS-7-step10-<UTC>`, BLOCKED. Never re-register a record.
- **(T2′) / (T2″) or (AB′) / (AB″) red:** this contradicts a landed green record. Report the Im `Y_in` triple and residuals, open a known-issues entry if it reproduces, then park and mark BLOCKED. Never re-band, never apply κ.
- **(N1′), (N0′) or (N-C) red:** wiring. BLOCKED.
- **(K) red:** revert both tracked outputs, BLOCKED.
- **`XT128` killed or red, `XT64` green:** commit `XT64`'s `metrics.json` / `COMPARISON.md` entry with its log as a partial record (the 127.74 column reads "—"). T23 stays open, and the next slot runs `XT128` alone with no edit.
- No `PREDICTED RED:` line.

**Outcome (2026-10-06 19:30 CDT slot): GREEN.** `-n 8`, 80 798 cells, both stages Status 0 (`20261007T003240Z_ANS-7-item288-XT64.log` 154 s, mesh 24.8 s, drives 39.0 / 43.6 / 42.8 s; `20261007T003522Z_ANS-7-item288-XT128.log` 165 s, mesh 24.1 s, drives 42.9 / 46.2 / 49.8 s). Script mtime 00:32:22Z precedes XT64's start (rule (i)). (AB′) / (AB″), (T2′) −2.535844e-6 (0.197 bands) / (T2″) −2.264621e-6 (0.140 bands), c, r, R_in, Q all equal the predictions to the printed digits. (R-T) ≤ 9.97e-13 (predicted ≤ 1e-12, met). (N1′) 1 577×, (N0′) `[]`, (N-C) 7 357× / 3 655× `NODE_S_RTOL`, (N2′) / (N2″) 493.78× / 493.82×. (K) green. Censuses unchanged (figures 88 `broken=0`; docrefs `dead=0`, `guide=0`, exit 2). Pin + unit ran as one window (`20261007T003811Z_ANS-7-item288-pin.log`, 24 passed). `ANS-7` stays 🟡.

## §7 row moved verbatim 2026-10-07 (rotate_plan_archive.py chunks)

 The tuned copper F-small birdcage: `Z_in(f)`, f₀ and Q against AED | 🟡 Adjudicated 2026-09-30: AGREE on the mode, on Route A against Route B and on the copper contribution, with the absolute rows INCONCLUSIVE at degree 1. Steps 5 / 5b (2026-10-03) re-recorded the fixture on exact rectangles with κ off at degree 1 and degree 2, as records. Step 6 (2026-10-04, item 214) added the resolution / Cu-only knobs (unset path = item 196 to ≤ 2.2e-12) and held T11's identity on an h = 0.01 rung at degree 2 (h-move f₀ +0.041 %, R_in −0.69 %, Q −0.70 %; a record, not an h-statement). Step 7 (T21, item 217, landed as a record by the 10-04 10:30 review): `C_mode` tunes the degree-2 in-model zero to 4.82e-6 of 64.000 MHz, a pre-registered red against the 1e-6 band (known-issues 2026-10-04); circuit = field green. Step 8 (T22, item 228, landed as a record by the 10-04 18:00 review): the same construction at 127.74 MHz reads 1.390e-6, red inside its pre-registered range; the Route A − Route B gap alone (−2.26e-6, flat in f) exceeds the band. Step 9 (T20, item 230, 2026-10-04 16:30 slot): `ans:7`'s tracked outputs now carry the degree-2 exact-rectangle rows (`metrics.json` `degree2_exact`, `COMPARISON.md`), with (R2x) re-asserted through the example at 1.6e-12 / 3.1e-12 (Cu R_in / Q), (L) at 1.274e-5 and (AB) at 1.22–1.27e-5, and the degree-1 ragged rows labelled superseded. Steps 7 / 8 re-posed 2026-10-05 (weekly) as 7′ / 8′ (§10 T21′ / T22′): tune on C_in, the in-model root, with the frequency band set from the sparse grid's measured resolution — the 1e-6 band was an identity band on a frequency target and the Route A − Route B gap (≈ −2.3e-6) is a floor on this mesh. Step 7′ (T21′, §9 item 262, 2026-10-06 06:00 slot): driven at C_in = 63.590091 pF (item 217's root, read from `tune_p2.json`), the degree-2 in-model zero sits at f₀,A/64 MHz − 1 = −2.536e-6, inside the pre-registered grid resolution `T21P_ZERO_RESOLUTION` 1.287e-5 (0.197 bands); c = C_mode/C_in − 1 = −1.494697e-5, r = f₀,A/f₀,B(C_in) − 1 = −2.536e-6 (flat against −2.531e-6 / −2.535e-6 / −2.263e-6); (AB′) circuit = field 1.237e-5 / 1.258e-5 / 1.287e-5 at `REDUCTION_BAND`; a ×1-mesh, degree-2 record of C_in, not an absolute tuning claim (`20261006T110217Z_ANS-7-item262-W1-cin.log`). Step 8′ (T22′, §9 item 283, 2026-10-06 12:00 slot): at 127.74 MHz on item 228's C_in = 15.274844 pF the in-model zero sits at f₀,A/127.74 MHz − 1 = −2.265e-6, inside `T22P_ZERO_RESOLUTION` 1.616e-5 (0.140 bands); c = −7.754048e-6, r = −2.265e-6; (AB″) 1.590e-5 / 1.601e-5 / 1.616e-5 at `REDUCTION_BAND`; same record status as 7′ (`20261006T171745Z_ANS-7-item283-W2-cin128.log`). Step 10 (T23, §9 item 288, 2026-10-06 19:30 slot): `ans:7` stages `XT64` / `XT128` re-assert (T2′) / (T2″) and (AB′) / (AB″) through the example and reproduce the stored C_in runs to ≤ 9.97e-13 (R-T); the tuned-to-Larmor rows are in `metrics.json` `degree2_exact.tuned_larmor` and `COMPARISON.md`, records only (`20261007T003240Z_ANS-7-item288-XT64.log`, `20261007T003522Z_ANS-7-item288-XT128.log`). **Done-when:** The row stays 🟡 until §10 T11 / T12 supply the order- and h-matched rung and the weekly re-reads (T13). *History: `docs/planning/chunks/ANS-7.md`.* 

### §9 item 310 — 2026-10-08 (chain T28)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 chain T28 (PROJECT_PLAN.md :4595–4601), quoted, not changed: "T21′'s construction (`FEM_EM_ANS7_EXACT_TUNE_C=cin`) on step 6's h = 0.01 degree-2 mesh: the in-model zero within that grid's own measured resolution (pre-registered from the rung's node residuals), circuit = field at `REDUCTION_BAND`; printed: C_in(h = 0.01)/C_in(×1) − 1 beside step 6's `f₀` +0.041 %. Target: identity; the first h-move of C_in. Serial on nothing. Heavy." Exact rectangles, κ off (ruling 2026-10-02). Lossless C: the ESR ruling (§9 :3749) scopes r to new steps, and T28 is enumerated lossless, so it is labelled lossless. **No band, tolerance, Done-when or target is set or moved.**

`ANS-7` step 12: on item 214's h = 0.01 rung (140 270 cells, 968 016 unknowns, degree 2), compute C_in(h). C_in(h) is the root of Im `Y_in`(64.0 MHz) on a solved h = 0.01 Cu 4×4, found by `_c_in` (`tests/validation/test_ans7_exact_strip.py:1199`, unchanged). Then drive Route A at C_in(h) at T21′'s three points 63.9 / 64.0 / 64.1 MHz and read the in-model zero.
- **Executor:** implementer. Tests only (`tests/validation/test_ans7_exact_strip.py`). No `src/`, `ans:7`, `TH-14` or `test_ans7_sparse_route_b.py` edit.
- **Tier:** heavy, `-n 12` (item 214's rung width), complex build. **Class:** chain step (solve).
- **Time:** ≈ 0 measured (no Route A drive has run on this rung) / ≈ 80 predicted min (15 fixed + ≈ 15 code + six windows ≈ 25–30, W0 included + ≈ 15 spare). On a slow box this runs past the 60-min timebox. The per-window split below makes it parkable partway through.
- **Depends on nothing.** Inputs on disk (checked 2026-10-08):
  - item 214's five `logs/ans7_exact/node_<mhz>_p2_h0.01_cu.json` (62.0 / 63.2 / 64.5 / 65.8 / 67.0 MHz; `resolution` 0.01, `ranks` 12, `field_solve` false, `cells` 140 270);
  - the ×1 `tune_p2.json` (C_in(×1) = 6.359009101426113e-11 F).
  - Never write the ×1 `tune_*.json` files or the h = 0.01 node files.

**Code** (every knob-unset path stays byte-identical; item 262 / 283's `cin` paths stay behaviourally identical).
1. **The h = 0.01 combination becomes legal.** `FEM_EM_ANS7_EXACT_RESOLUTION=0.01` (`_resolution()`, `:159`) is accepted with `FEM_EM_ANS7_EXACT_TUNE` ∈ {`cmode`, `inmodel`}, at 64 MHz only (`TUNE_MHZ` unset) and Cu only (`COLUMNS=cu`). Any other combination raises `ValueError`.
2. **`_active()` / `_stored_active()`** (`:954` / `:1066`) gain the h = 0.01 set: the five `_h0.01_cu` nodes, Cu only, with the tie checked by `_assert_tie(..., res=0.01)`. Output files are suffixed `tune_p2_h0.01.json` and `tune_inmodel_cin_p2_h0.01.json`.
3. **`tune_cmode`** (`:1218`) under the h set runs with `pec_column=False`, `field_solve=False` and `resolution=0.01` passed to `_build_ladder`, as item 214 did. It writes C_mode(h), C_in(h) and z64. Item 157's trap forbids `pec_column=False` only together with `field_solve=True`; Q3 asks the review to confirm.
4. **`tune_inmodel`** (`:1404`):
   - a per-point knob `FEM_EM_ANS7_EXACT_TUNE_POINT` ∈ {unset, 63.9, 64.0, 64.1, `assemble`}. **Unset gives the existing behaviour** (ruling 310-Q2), and W0 shows it digit-identical;
   - a numeric value runs one Route A drive per window, and rank 0 writes `tune_inmodel_cin_p2_h0.01_<mhz>.json`;
   - `assemble` is a pure path: it reads the three point files and asserts. Any other value raises `ValueError`;
   - the (W) degree-1 ×1 wiring drive runs in the 64.0 window only.

**Anchors** (asserted; bands imported by name).
- **(AB‴) circuit = field at C_in(h):** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` at all three points. The band is 1e-3, `tests/validation/test_port_lumped_rlc_termination.py:88`, imported through `test_port_circuit_layer_field`. Route B is the five-node h = 0.01 grid, interpolated with `bary_interp`.
- **(T2‴) the in-model zero exists inside the bracket** (review ruling 310-Q1): Im `Y_in`,A changes sign, − → +, exactly once between the bracketing points 63.9 and 64.1 MHz. Asserted.
  - **The zero's offset is printed, never asserted.** \|f₀,A/64 MHz − 1\| is printed with `T21P_ZERO_RESOLUTION` (1.287e-5, `test_ans7_exact_strip.py:1001`) beside it, for scale only. No band is transplanted.
  - **Why T28's "within that grid's own measured resolution" could not be pre-registered.** T21′'s band came from T21's (AB) node residuals at C_mode (`20261004T140738Z_ANS-7-item217-W2-inmodel.log:2046–2048`). The rung's node files carry no identity residuals (`field_solve` false, `identity_residual` null), so there is nothing to register it from. This gap is named for the 10-10 weekly.
- **(W)** unchanged: ≤ `NODE_S_RTOL` against item 195's stored point.

**Negative controls (rule (e)).**
- **(W): asserted.** Backed by item 262 W1 (4.034e-14) and item 283 (1.127e-14 / 1.470e-13): same fixture, same comparison.
- **(N1‴) untuned zero: predicted, printed.** \|f₀,B(C_F)/64 MHz − 1\| on the h = 0.01 Cu grid is ≥ 10 × `T21P_ZERO_RESOLUTION`.
  - It is not asserted: the comparison has run on the ×1 grid (2.029559e-2) but never on this grid.
  - Step 6's rung f₀ at C_F is 62.726847 MHz (≈ −1.99e-2), so a ratio of ≈ 1 500× is expected.
- **(N2‴): predicted, printed.** The circuit at 0.5·C_in(h) against the in-model at C_in(h) is ≥ 10 × `REDUCTION_BAND` (×1 measured 493.78×).
- **(N-h): predicted, printed.** C_in(h)/C_in(×1) − 1 is nonzero, by at least 10 × the secant tolerance. The h-move is the step's reading, so it is not gated.

**Predictions** (written first, printed).
- C_in(h)/C_in(×1) − 1 in [−1.5e-3, −0.3e-3]. Step 6 moved f₀ by +0.0411 % at fixed C. Since f ∝ C^(−1/2) at first order, the expected move is ≈ −8.2e-4.
- (AB‴) in [5e-6, 3e-5]. ×1 read 1.237–1.287e-5; the h-trend is unmeasured.
- r(h) = f₀,A/f₀,B(C_in(h)) − 1 in [−5e-6, −1e-6]. ×1 is flat at −2.5e-6 / −2.3e-6.
- R_in and Q at C_in(h) within 1 % of step 6's rung rows (18 313.602 Ω / 970.439). Printed only; item 262 missed a tighter Q prediction.

**Windows** (`-n 12`; each bare, whole and foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE` unset; no edit after W1 starts, rule (i)).
- **Orphan check first:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0.
- **The trap check comes before any window** (ruling 310-Q3). Quote item 157's trap text and confirm that `pec_column=False` with `field_solve=False` lies outside it. If the executor cannot confirm this, it stops and reports.
- Each window is under 10 min, so the §5.1 capture idiom is not required.
- **W0 (the knob-unset path, digit-identical; ruling 310-Q2):** item 262's W1 command verbatim, label `ANS-7-item310-W0-cin64`, `-n 8`, 257 s (measured, item 262; reproduced by item 283 W1 at 258 s). It must reproduce, to the printed digits:
  - (T2′) −2.535844e-6;
  - (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5;
  - r −2.535846e-6.
  - **Red ⇒ revert, BLOCKED.**
- **Price (nearest measurement).** Item 214 W1–W3 at h = 0.01, `-n 12`, on a box running ≈ 1.4× slow: mesh 53.6–55.4 s; one Cu 4×4 at 143.5–159.8 s; `ru_maxrss` 32.8–36.5 GiB; W2 used 536 s of 590 for its node set.
  - A Route A drive is one factorisation of the same operator, so ≈ 150 s is predicted.
  - W1 ≈ 220–280 s; each point window ≈ 230–300 s; W5 takes seconds.
- **W1 (C_in(h)):**

  `scripts/testing/run_and_log.sh ANS-7-item310-W1-cin-h001 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=cmode timeout -k 30 590 mpiexec -n 12 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
- **W2 / W3 / W4 (one drive each):** W1's command with `FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_TUNE_POINT=63.9` (then 64.0, then 64.1). Labels `ANS-7-item310-W2-p63.9` / `-W3-p64.0` / `-W4-p64.1`. Run them serially, never concurrently (≈ 36 GiB each).
- **W5 (assembly, closing):** `FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_TUNE_POINT=assemble`, `-n 2`, `timeout -k 30 180`, label `ANS-7-item310-W5-assemble`.

**Traps already paid for.**
- Items 217 / 228 / 262 / 283's lists, whole:
  - κ `None`;
  - no `ans:7` call;
  - `z_from_s(S)`, never `z_matrix`;
  - a real-scalar secant;
  - rank 0 writes, then a barrier;
  - a killed window is not a result: run `pgrep -c python3`, report, never lengthen the timeout;
  - never route through `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ` (`_nodes_env` validates against the ×1 nodes).
- **The h = 0.01 nodes are `_cu`-suffixed,** so `COLUMNS=cu` is required.
- **A missed ×1 literal** (cells 80 798, the seven-node grid, `tune_p2.json`) shows up as a tie or grid failure on a correct measurement.
- **Price on the measured mesh time.** Item 214's mesh ran over its 30–40 s prediction (54 s).

**Scope.** Covered: h = 0.01, degree 2, κ off, Cu, 64.000 MHz, lossless C. **Not claimed:**
- any absolute tuned C, f₀, R_in or Q;
- an h-statement: this is one rung, and T29 / T13 read it beside T12′;
- 127.74 MHz;
- ESR;
- any AED or Palace comparison.

**Status it can move.**
- **No glyph.** The §7 row's Done-when is the weekly's T13 re-read (PROJECT_PLAN.md :3444), so `ANS-7` stays 🟡.
- The row gains a "step 12 (T28)" record clause: C_in(h)/C_in(×1) − 1, (AB‴) and the zero.
- The 10-10 weekly can read T29.
- **Rule (j):** declined, as for items 228 / 262 / 283 / 288. No glyph moves.

**Negative result.**
- **(W) red:** BLOCKED (wiring).
- **(AB‴) red with (W) green:** circuit and field disagree on the rung. Print the three residuals, open a known-issues entry, park and mark BLOCKED. Never apply κ.
- **No crossing in (63.9, 64.1):** BLOCKED. Report the Im `Y_in` triple. With C_in(h) taken from the rung's own grid, this would contradict r's flat record.
- **A point window killed:** commit nothing for that point. The next slot reruns that point alone, with no edit.
- No `PREDICTED RED:` line.

**Review rulings (2026-10-08 03:00 daily review).** The rulings below are applied to the text above.
- **(310-Q1)** Neither (b), three C_mode(h) drives first, nor (c), an imported `T21P_ZERO_RESOLUTION`. T28's asserted anchors are the ones that can be asserted without a new band:
  - Im `Y_in` changes sign inside the bracket, so the in-model zero exists there;
  - circuit = field at `REDUCTION_BAND` at C_in(h).
  - The zero's offset is printed with `T21P_ZERO_RESOLUTION` beside it for scale, and never asserted. "Within that grid's own measured resolution" could not be pre-registered, because the rung's node files carry no identity residuals. It is named for the 10-10 weekly. No band is transplanted.
- **(310-Q2)** The per-point split and the new `FEM_EM_ANS7_EXACT_TUNE_POINT` knob are accepted. Unset gives the existing behaviour, shown digit-identical by an existing window (W0). Never one long window.
- **(310-Q3)** The executor confirms the trap from item 157's trap text before the first window. If it cannot, it stops and reports.

**Outcome (2026-10-08 06:00 CDT slot): GREEN, all six windows.** Trap check (310-Q3), from item 157's text (`:86`): "`pec_column=False` together with `field_solve=True` raises `KeyError` on rank 0 only (`rec["pec"]` …)". In `_build_ladder` the `rec["pec"]` access (`test_th14_birdcage_copper.py:412`) comes after the `if not field_solve: … continue` at `:362–370`, and the `class_ds` path is gated by `if pec_column:` (`:357`). So `pec_column=False` with `field_solve=False` is outside the trap, which is the combination item 214 ran on this rung. Orphan check 0 before W0 and after W5. Module mtime 11:05:15Z precedes W1's start at 11:10:45Z (rule (i)).
- **W0** (`20261008T110604Z_ANS-7-item310-W0-cin64.log`, `-n 8`, 262 s): reproduces to the printed digits: (T2′) −2.535844e-6 (0.197 bands), (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5, r −2.535846e-6. The ×1 `tune_inmodel_cin_p2.json` it rewrites differs from the stored copy only at the 4e-8 relative level in the residuals.
- **W1** (`…111045Z…-W1-cin-h001.log`, `-n 12`, 252 s; ladder 144.8 s; 140 270 cells; 15 passed): C_mode(h) 63.643437 pF; C_in(h) = 6.364316664701229e-11 F (\|g\| 1.0e-16 S after 2 iterations). **C_in(h)/C_in(×1) − 1 = +8.346526e-4.** The predicted range [−1.5e-3, −0.3e-3] is **NOT met on sign**; the magnitude matches the brief's 8.2e-4. The brief's sign was wrong: f₀ rose +0.0411 % at fixed C, and f ∝ C^(−1/2), so tuning back down to 64 MHz needs more C. The drafter wrote the f₀ move, not the C move needed to undo it. (N-h) 8.35e-4 is far above 10 × the secant tolerance (3.9e-13). (T1) on the five h nodes at C_mode(h): Δ_LOO 3.3e-10 / 1.9e-8 / 4.4e-8; interpolated-vs-solved Z64 4.3e-9. (N1) ×1 wiring ≤ 2.2e-12.
- **W2 / W3 / W4** (`…111513Z…-W2-p63.9.log` 216 s, `…111904Z…-W3-p64.0.log` 236 s, `…112308Z…-W4-p64.1.log` 214 s; mesh 32.7–34.0 s, drive 98.2–105.7 s): (AB‴) 1.239e-5 / 1.262e-5 / 1.289e-5, inside the predicted [5e-6, 3e-5]. Im `Y_in`,A −1.618360e-4 / +2.640932e-7 / +1.626669e-4 S. (W) 1.720e-13 (×1 degree-1 drive, 28.0 s, 64.0 window).
- **W5** (`…112649Z…-W5-assemble.log`, `-n 2`, 70 s, 15 passed): (T2‴) one − → + change and no + → − change in (63.9, 64.1). f₀,A 63 999 837.080 Hz, so \|f₀,A/64 − 1\| = 2.545621e-6 = 0.198 × `T21P_ZERO_RESOLUTION`; this is printed, never asserted. f₀,B(C_in(h)) / 64 − 1 = −5.7e-13. r(h) −2.545621e-6, inside the predicted [−5e-6, −1e-6]. R_in 18 592.58 Ω is +1.52 % against the step 6 rung row, which **misses** the printed 1 % prediction; Q 964.43 is −0.62 %, which meets it. The comparand is step 6's rows at C_F, a different C, so the miss is on the comparand, not a gate. (N1‴) 1.989302e-2 = 1 545.7 × (predicted ≈ 1 500×). (N2‴) 493.62 × `REDUCTION_BAND`.
- **No glyph.** `ANS-7` stays 🟡, and the §7 row gains the step 12 clause. The gap is named for the 10-10 weekly: the rung has no registered "own resolution" (310-Q1).

### §9 item 311 — 2026-10-08 (chain T26)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* The 03:00 review's rulings 311-Q4 to Q6 are applied to the text below.

**Licence.** §10 T26 (PROJECT_PLAN.md :4586–4592), quoted, not changed: "the in-model zero of Im `Y_in` at C_in(r) at 64.000 MHz within `T21P_ZERO_RESOLUTION`, circuit = field at `REDUCTION_BAND` at C_in(r); printed: Q(r)/Q(0) against the closed-form 1/Q = 1/Q_0 + 1/Q_C with Q_C = 1/(ωC_in r), C_in(r)/C_in(0) − 1. Target: identities. Serial on T25 and `PORT-24` step 2. Heavy (257 s class)." The ESR ruling (§9 :3735–3754) sets r = 0.2 Ω on both sides at once; lossless records stay, labelled lossless.

**Depends on items 308 and 309.** If item 308 (`PORT-24` step 1, `esr_ohm` in the circuit layer) or item 309 (`PORT-24` step 2, the series-RC sheet) did not land, skip this item and leave it ⬜. Their Done-when is imported from the `PORT-24` row (PROJECT_PLAN.md :2716), not restated. Lossless records stay valid and are labelled lossless. r never re-reads a lossless reading (ruling, §9 :3749–3750).

`ANS-7` step 11: tune the exact-sheet F-small birdcage on C_in with r = 0.2 Ω on both sides. The work is on the ×1 mesh (h = 0.015, 80 798 cells), degree 2, κ off, Cu, at 64.000 MHz.
- **Circuit prediction (ruling 311-Q4):** C_in(r) comes from the stored ×1 `z64` in `tune_p2.json` through the circuit layer with `esr_ohm = 0.2`. It is pure and takes seconds.
- **The check:** the in-model solve with series-RC sheets at (r, C_in(r)).
- **Executor:** implementer. Tests only (`tests/validation/test_ans7_exact_strip.py`). No `src/` edit beyond what items 308 / 309 landed; no `ans:7` edit (T27 does that).
- **Tier:** heavy, `-n 8` (T21′'s record width), complex build. **Class:** chain step (solve).
- **Time:** ≈ 9 measured (W1 + W2 ≈ 257 s each; item 262, item 283 W1) / ≈ 45 predicted min (15 fixed + ≈ 15 code + ≈ 9 windows + ≈ 6 spare).

**Code** (every knob-unset path byte-identical; item 262 / 283 / 310's paths behaviourally identical).
1. **Knob `FEM_EM_ANS7_EXACT_ESR_OHM`** ∈ {unset, `0.2`}. Anything else raises `ValueError`. The knob requires:
   - `TUNE=inmodel`, `TUNE_C=cin`, `TUNE_MHZ` unset;
   - `RESOLUTION` unset;
   - `TUNE_POINT` unset.
2. **C_in(r), pure.** `_c_in` (`:1199`) finds the root of Im[1/Z_in] at 64.0 MHz on the stored `z64` (Cu) from `tune_p2.json`. Every gap termination is item 308's `esr_ohm`-carrying termination, Z = r + 1/(jωC), and the feed is driven across its gap's (r + C) branch (ruling 311-Q5; ESR ruling (i)).
   - Use item 308's landed API by name, with no local re-implementation of the termination. If `route_b_z` does not take `esr_ohm`, thread it through as a keyword whose default is 0.
   - The secant start is `c_in` (lossless).
   - Rank 0 writes `tune_p2_esr0.2.json`, holding `c_in_r`, `c_in_0`, `esr_ohm` and `f_hz`. It never overwrites `tune_p2.json`.
3. **Route A at (r, C_in(r)).** `_route_a_solve(…, c_f=c_in_r, degree=2)` places item 309's series-RC sheet, `series_rlc_impedance(R=r, L=0, C)`, on the three non-feed gaps.
   - Z_in,A adds the feed gap's (r + 1/(jωC_in(r))) branch in parallel, in place of `add_c1`'s pure C.
   - Route B is `route_b_z(bary_interp(xn, z_nodes_cu, g − 64.5), f, c_in_r, esr_ohm=r)` at the three T21′ points 63.9 / 64.0 / 64.1 MHz.
   - Rank 0 writes `tune_inmodel_cin_esr0.2_p2.json`.
4. **Tests** (skip unless the ESR knob is set): `test_t26_ab_circuit_equals_field_esr`, `test_t26_inmodel_zero_esr`, `test_t26_n_lossless_separates`.

**Anchors** (asserted; bands imported by name).
- **(AB-r) circuit = field at C_in(r):** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` (1e-3, `tests/validation/test_port_lumped_rlc_termination.py:88`) at all three points.
- **(T2-r) the in-model zero:** exactly one − → + crossing of Im `Y_in`,A in (63.9, 64.1) MHz, and \|f₀,A/64 MHz − 1\| ≤ `T21P_ZERO_RESOLUTION` (1.287e-5, `test_ans7_exact_strip.py:1001`). This is the band T26 names, imported, not chosen.
- **(W)** unchanged: ≤ `NODE_S_RTOL` against item 195's stored point.

**Negative controls (rule (e)).**
- **(W): asserted.** Backed by item 262 W1 (4.034e-14) and item 283 W1 (1.127e-14): same fixture, same comparison.
- **(N-L) the band separates ESR from lossless: predicted, printed.** The lossless circuit (r = 0 at C_in(0)) against the in-model at (r, C_in(r)) at 64.0 MHz is ≥ 10 × `REDUCTION_BAND`.
  - Not asserted: no ESR comparison has run on this fixture.
  - Expected large: with R_in falling by roughly Q(r)/Q(0), S11 at the zero moves by O(1).
- **(N-r0) r = 0 reproduces C_in: predicted, printed.** `_c_in` with `esr_ohm = 0` on the same `z64` equals `tune_p2.json`'s `c_in` bit for bit. Not asserted: the comparison has never run through item 308's path.
- **(N2-r): predicted, printed.** The circuit at 0.5·C_in(r) against the in-model at C_in(r) is ≥ 10 × `REDUCTION_BAND` (lossless 493.78×).

**Printed, never asserted (ruling 311-Q6).**
- C_in(r)/C_in(0) − 1. **Predicted** in [−1e-3, +1e-3], magnitude O(1/Q_C²) ≈ 3e-5. ESR ruling (i) says uniform r leaves `C_mode` unchanged, but the C_in root of Im `Y_in` may move at second order.
- Q(r)/Q(0) beside the closed form 1/Q = 1/Q_0 + 1/Q_C, with Q_C = 1/(ωC_in r) per capacitor, as the weekly wrote it (Q_0 is the lossless-capacitor coil's loaded Q, Q(0); renamed from "Q_Cu" by the 2026-10-10 weekly).
  - At C_in = 63.59 pF and 64 MHz, Q_C ≈ 195.5.
  - With Q(0) = 971.07 (item 262), the closed form gives Q(r) ≈ 163 and Q(r)/Q(0) ≈ 0.168. **Predicted**, labelled predicted in the print.
- R_in(r), f₀,B(C_in(r))/64 MHz − 1 (predicted within ±1e-9, by construction), and r_AB = f₀,A/f₀,B − 1 (predicted in [−5e-6, −1e-6]; lossless −2.536e-6).
- (AB-r) predicted in [5e-6, 3e-5] (lossless 1.237–1.287e-5). (T2-r) predicted ≈ 2.5e-6, about 0.2 bands: **predicted green**.

**Windows** (`-n 8`; each bare, whole and foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE` unset; no edit after W1 starts, rule (i)).
- **Orphan check first:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0. Each window is under 10 min, so the §5.1 capture idiom is not required.
- **W1 (knob unset, digit-identical):** item 262's W1 command verbatim, label `ANS-7-item311-W1-cin64-lossless`. It must reproduce (T2′) −2.535844e-6 and (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5 to the printed digits. Price 257 s. **Red ⇒ revert, BLOCKED.**
- **W2 (closing):**

  `scripts/testing/run_and_log.sh ANS-7-item311-W2-cin64-esr "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_ESR_OHM=0.2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`

  Price ≈ 257 s (the same drive count as item 262; the C_in(r) root is pure). On a 1.4× slow box ≈ 365 s, inside 590.

**Traps already paid for.**
- Items 217 / 228 / 262 / 283's lists, whole:
  - κ `None`;
  - no `ans:7` call;
  - `z_from_s(S)`;
  - a real-scalar secant;
  - rank 0 writes, then a barrier;
  - a killed window is not a result;
  - never `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ`.
- **Never enter 0.2 Ω as a parallel R.** ESR ruling (iii): a parallel R ‖ C with R = 0.2 Ω shorts the gap. The sheet is the series RC.
- **The feed's branch carries r too** (ruling 311-Q5). Leaving the feed lossless gives a non-uniform r, which item 308's planted control shows breaks the uniform shift on every mode.
- **Never write `tune_p2.json` / `tune_inmodel*_p2.json`.** Read them only.
- **r never re-reads a lossless record.** The lossless rows stay as they are.

**Scope.** Covered: ×1 mesh, degree 2, κ off, Cu, 64.000 MHz, r = 0.2 Ω uniform on all four capacitors on both sides. **Not claimed:**
- a measured ESR or its frequency dependence;
- any absolute tuned C, f₀, R_in or Q;
- an h-statement;
- 127.74 MHz;
- any AED re-run or Palace comparison.

**Status it can move.**
- **No glyph.** The `ANS-7` Done-when is the weekly's T13 re-read (PROJECT_PLAN.md :3444), so the row stays 🟡.
- The row gains a "step 11 (T26)" record clause: C_in(r), (AB-r), (T2-r), Q(r)/Q(0) beside the closed form.
- T27 (`ans:7` ESR rows) becomes queueable.
- **Rule (j):** declined, as for items 262 / 283 / 288 / 310.

**Negative result.**
- **W1 red:** the edit changed item 262's path. Revert it and mark BLOCKED.
- **(W) red:** BLOCKED (wiring).
- **(AB-r) red with (W) and W1 green:** the series-RC sheet and the circuit's (r + C) termination disagree. Print the three residuals, open a known-issues entry, park and mark BLOCKED. Never apply κ, and never adjust r.
- **No crossing in (63.9, 64.1):** BLOCKED. Report the Im `Y_in` triple.
- **(T2-r) band half red:** report r_AB and C_in(r)/C_in(0) − 1, park BLOCKED, and the weekly decides. Never re-band.
- No `PREDICTED RED:` line.

**Review rulings, second pass (2026-10-08 03:00 daily review, items 310 / 311).** (310) The new W0 (item 262's W1 verbatim, digits reproduced) is accepted; the per-drive ≈ 150 s is inferred, so a 590 s kill parks the item as the price record. (311) If item 308's landed `route_b_z` (or its test-side equivalent) lacks `esr_ohm`, the executor threads it through as a test-side keyword defaulting to 0, with the knob-unset path shown unchanged by W1; no `src/` signature beyond item 308's is added in this item.

**Review amendment 2026-10-11 (item 311; 2026-10-10 weekly, §10 T25′ / T26 and the "Q_Cu" naming ruling).**
- **Dependency.** Item 311 now depends on **item 363** (T25′, `PORT-24` step 1 landed under the re-stated Im clause) **instead of 308**. Item 308 is superseded and never lands as written. If item 363 has not landed, skip and leave 311 ⬜.
- **Item 309 has landed.** Commit `91254c5a`, "test(port24): PORT-24 step 2 -- in-model series-RC sheet, three arms green (§9 item 309)", journaled by `7f388752`. Its one-line DONE record is in `docs/planning/plan-archive.md` ("309 DONE 10-08 12:00 slot"), and its readings are in `PORT-24.md` § "§9 item 309". The 309 half of the dependency is met.
- **Naming.** "Q_Cu" is now **Q_0**, the lossless-capacitor coil's loaded Q. That is Q(0), item 262's 971.07 at the ×1 C_in. The two places this brief printed the formula (the licence quote and "Printed, never asserted") are edited in place to 1/Q = 1/Q_0 + 1/Q_C. Nothing else in the text changes.
- **One fact the executor needs, read from the code on 2026-10-11.**
  - Item 363 lands `esr_ohm` only in `src/fem_em_solver/ports/ring_modes.py`, plus the test-local `capacitor_impedance_esr` / `tuned_input_esr` in `tests/unit/test_port24_esr_circuit.py`.
  - `ans:7`'s `_route_a_solve` (`07_tuned_copper_birdcage_driven_sweep.py:541`) sets every terminated sheet to `port_impedance_ohm = z_cap(ω, c_f)`. Its `route_b_z` (`:181`) and `add_c1` (`:171`) take no r either.
  - So, under the second-pass ruling (311) above, item 311 builds its ESR Route A and Route B **test-side** in `test_ans7_exact_strip.py`. The terminated sheets take `dataclasses.replace(spec, port_impedance_ohm=series_rlc_impedance(f, r_ohm=r, l_h=0.0, c_f=C))`, item 309's sheet. `_route_a_solve`'s body is copied with a provenance comment (commit and line span). The feed's (r + C) branch is formed locally.
  - `ans:7` is not edited; T27 (item 378) lifts the ESR terminations into the example.
- **Entry:** `311. ⬜ — ANS-7 step 11 / T26 … (chain T; implementer; heavy, -n 8; ≈ M measured / 45 predicted min; depends on 363)`.

**Outcome (2026-10-11 04:30 CDT slot): GREEN.** Code test-side only (`_esr_ohm`, `_route_b_z_esr`, `_c_in_esr`, `_route_a_solve_esr` copied from `ans:7` `:541–574` at `2f0f83d1`, fixture `tune_esr`, three `test_t26_*`; `tune_inmodel` skips under the knob). W1 `20261011T093544Z_ANS-7-item311-W1-cin64-lossless.log`, Status 0, 264 s, 16 passed: (T2′) −2.535844e-6, (AB′) 1.237e-5 / 1.258e-5 / 1.287e-5 — digit-identical to item 262. W2 `20261011T094028Z_ANS-7-item311-W2-cin64-esr.log`, Status 0, 260 s (W2-part 157.4 s; drives 46.1 / 38.9 / 44.1 s), `-n 8`, 80 798 cells, 15 passed. C_in(r) = 6.359301781809668e-11 F (3 secant iterations): C_in(r)/C_in(0) − 1 = +4.602610e-5 (predicted in ±1e-3, met; larger than the O(1/Q_C²) ≈ 3e-5 estimate). (N-r0) bit-equal, met. (W) 3.927e-14 (asserted). (AB-r) 1.233e-5 / 1.258e-5 / 1.283e-5 (asserted ≤ `REDUCTION_BAND`; predicted [5e-6, 3e-5], met). (T2-r) Im `Y_in`,A −1.616350e-4 / +2.676561e-7 / +1.624710e-4 S, one − → + crossing, f₀,A/64 MHz − 1 = −2.583111e-6 = 0.201 × `T21P_ZERO_RESOLUTION` (asserted; predicted ≈ 2.5e-6, met). f₀,B/64 MHz − 1 = +1.89e-12 (met); r_AB = −2.583113e-6 (predicted [−5e-6, −1e-6], met). R_in(r) 3 175.98 Ω (ratio 0.169483 to R_in(0)); Q(r) 164.543, Q(r)/Q(0) = 0.169446 beside the closed form 1/Q = 1/Q_0 + 1/Q_C (Q_C 195.524, Q 162.754, ratio 0.167603): Q(r) sits +1.10 % above it (printed). (N2-r) 494.39 × (met). **(N-L) 4.97 × `REDUCTION_BAND` on the bare S11, predicted ≥ 10 × — NOT met** (printed, unlabelled control: a negative result reported, not ruled on; on Z_in the same comparison moves 4.88, O(1) as the brief expected — the bare S11 excludes the feed branch, so r enters it only through the three terminated gaps). Note: W1 (item 262's verbatim command) rewrites `tune_inmodel_cin_p2.json` with the same digits, as items 283 / 310's W0 did. `ANS-7` stays 🟡.

### §9 item 364 — 2026-10-11 (T30 — `ANS-7`: the three-rung R_in / Q sequence from logged data; lane-1 revalidation)

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain T, step T30 (grep "T30 (lane 1, no solve). The three-rung"), quoted, not changed: "Pull Cu `R_in` and Q at h = 0.01 from step 6's rows (`ANS-7.md`, its log) beside ×1 (18 441.75 Ω / 977.26) and h = 0.005 (18 223.23 Ω / 965.55, T12′); print the two successive changes and their ratio; if the second change is at most half the first, write Q's bracket as a record on the row, else 'unconverged at h = 0.005' and F31 is the next rung. Target: a record from logged data. Serial on nothing." The T13 ruling (§10 rulings, 2026-10-10) binds: `R_in` and Q are not an agreement, and `ANS-7` holds 🟡 on the Q clause until T30 / F31 write Q's h-statement. **No band, tolerance, Done-when or target is set or moved.** The 0.5 test is the weekly's rule; the only band used is `NODE_S_RTOL`, imported.

**Why a window at all.** The record is arithmetic on three logged rows, but host `python3` is denied, and a record "from logged data" should assert that the data are what the logs say. So the item is one smoke window running a new pure module. It re-derives each rung's Cu rows from the stored node JSONs through the production `sparse_rows` path and asserts that they reproduce the logged digits. The successive changes and their ratio are then printed from the reproduced values.

- **Executor:** implementer. **Class:** lane-1 revalidation (no solve). **Tier:** smoke, `-n 1`.
- **Cost:** ≈ 10 s compute, priced from item 214's import smoke of `test_ans7_exact_strip.py` (7 s, `20261004T110430Z_ANS-7-item214-import-smoke.log`, "13 skipped in 4.02s") plus three `sparse_rows` calls on five nodes over the 51-point grid. Predicted ≈ 30 min: 10 fixed, 15 code, 1 window, 4 records. **Depends on nothing.**
- **Inputs (gitignored, checked on the host 2026-10-11):**
  - `logs/ans7_exact/node_{62.0,62.3,63.2,64.5,65.8,66.7,67.0}_p2.json` (×1, item 196);
  - `node_{62.0,63.2,64.5,65.8,67.0}_p2_h0.01_cu.json` (item 214);
  - `node_{62.0,63.2,64.5,65.8,67.0}_p2_h0.005_cu.json` (T12′, `xl` entry, `20261008T070008Z_ANS-7-step4p.log`).
  - If any is absent, mark BLOCKED. Never re-solve one; the h = 0.005 rung is an `xl` window.

**The logged rows (Cu, degree 2, κ off, 5-node Route B at `C_F`, full repr).** These are the comparands, entered as module constants with their citations:

| Rung | Cells | f₀ [Hz] | R_in [Ω] | Q | Log |
|---|---:|---:|---:|---:|---|
| ×1 (h = 0.015) | 80 798 | 62701081.939501435 | 18441.75458382276 | 977.2613692833095 | `20261004T112945Z_ANS-7-item214-W3-rung.log:2259–2261` ((S5) 5-node; also `20261008T070008Z_ANS-7-step4p.log:4171–4173`) |
| h = 0.01 | 140 270 | 62726846.66023456 | 18313.602153804146 | 970.4389973340836 | `…item214-W3-rung.log:2218` |
| h = 0.005 | 523 880 | 62751784.313079365 | 18223.232609192608 | 965.5456118046693 | `…step4p.log:4128` |

- The ×1 comparand is the **5-node** row, which matches the rungs' 5-node assembly. The 7-node record `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2` (`test_ans7_exact_strip.py:151–156`) differs by 9.2e-10 / 6.6e-10 ((S5), `:2260–2261`) and is printed beside it, not used.
- §10's "18 441.75 / 977.26" and "18 223.23 / 965.55" are these values rounded. Step 6's printed comparand `STEP6_RUNG_ROWS` (`test_ans7_exact_strip.py:1057`, 18 313.602 / 970.439) is the h = 0.01 row rounded.
- The ×1 resolution 0.015 is the module's own label (`test_ans7_exact_strip.py:147–150`, "×1 (h = 0.015)"). It is used only for a printed order estimate.

**Code.** One new pure module, `tests/validation/test_ans7_t30_rin_q_sequence.py`, `@complex_only`. No edit to any existing file.
- **Imports, by name, none restated:**
  - from `tests.validation.test_ans7_exact_strip`: `_x1_cu5_rows` (`:768`, the ×1 five-node Cu rows from item 196's seven JSONs; it returns `(rows, recs)`, and the rows are `rows["cu"]["rows"]`), `_load_nodes` (`:719`), `_uj`, `_assert_tie`, `_tie_from_json`, `ROWS`, `ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2`, `STEP6_RUNG_ROWS`;
  - from `tests.validation.test_ans7_sparse_route_b`: `ASSEMBLY5_NODES_MHZ`, `NODE_S_RTOL` (1e-9, `:87`, "backing 5.042e-14"), `sparse_rows`.
- **Rung rows:** `recs, missing = _load_nodes(ASSEMBLY5_NODES_MHZ, 2, res, ("cu",))` for res ∈ {0.01, 0.005}, then `sparse_rows(ASSEMBLY5_NODES_MHZ, {"cu": [_uj(r["z"]["cu"]) for r in recs]}, columns=("cu",))["cu"]["rows"]`. This is exactly the call `_loo` makes for `full` (`:268–272`, defaults `GRID_MHZ`, `C_F`), so the logged rows came from it. The ×1 row comes from `_x1_cu5_rows()`.
- **Record literals** (provenance comments citing the table's log lines; records, not bands): `T30_LOGGED_ROWS` (the nine values above), `T30_CELLS` (80 798 / 140 270 / 523 880; the h = 0.005 count from `step4p.log` and `docs/testing/xl-ledger.md` row "ANS-7-step4p"), and `T30_H_M = (0.015, 0.01, 0.005)`.
- **Tests**:
  - `test_t30_inputs`: (F);
  - `test_t30_reproduces_logged_rows`: (R-rep);
  - `test_t30_rung_swap_control`: (N-swap);
  - `test_t30_sequence_record`: prints only, plus the arithmetic identity in (A).

**Anchors (asserted).**
- **(F) inputs are the logged rungs.** For each rung's JSONs:
  - `degree == 2`, `kappa is None`, `resolution` equal to the rung's (None at ×1), and `columns == ["cu"]` on the h rungs;
  - one cell count, equal to `T30_CELLS` exactly;
  - `_assert_tie(cells, _tie_from_json(r["tie"]), res)` per file, as `test_p2_assembly` does (`:795`);
  - `f_hz` within 1 Hz of its node.
- **(R-rep) reproduction of the logged rows.** For each of the nine (rung, quantity) pairs, \|recomputed/logged − 1\| ≤ `NODE_S_RTOL`, the project's "unset digit-identical" reproduction band. It is applied as (R2) applies it (`:856–875`). **Predicted:** exactly 0, the same function on the same bytes. Printed with `!r`.
- **(A) the sequence arithmetic.** For R_in and Q, Δ₁ = v(0.01) − v(×1) and Δ₂ = v(0.005) − v(0.01) are computed from the reproduced values. Asserted: Δ₁ + Δ₂ = v(0.005) − v(×1) to ≤ `NODE_S_RTOL` relative. This is a bookkeeping guard on the pulled numbers, not a physics claim.

**Negative controls (rule (e)).**
- **(N-swap): asserted.** Each rung's recomputed R_in and Q against the *other two* rungs' logged values gives \|rel\| > 100 × `NODE_S_RTOL`. That is the factor item 288's (N-C) used, an existing precedent. It is backed by the same comparison on the same data: the logged h-moves (h = 0.01 − ×1)/×1 of −0.6949 % / −0.6981 % (`…item214-W3-rung.log:2234–2235`) and (h = 0.005 − ×1)/×1 of −1.1849 % / −1.1988 % (`…step4p.log:4146–4147`). Their difference, the smallest separation (h = 0.01 vs h = 0.005), is ≈ 4.9e-3 ≈ 4.9e6 × `NODE_S_RTOL`. A loader that read one rung's files under another's name fails here.
- **Printed only:** the 7-node ×1 record's R_in / Q against the 5-node ×1 row (predicted 9.2e-10 / 6.6e-10, as (S5)).

**Printed record (the step's reading; never asserted):**
- Δ₁, Δ₂ and ρ = Δ₂/Δ₁ for R_in and Q, and for f₀ for comparison with the T13 ruling's "equal to 3 %".
- The §10 branch, evaluated per quantity: `|Δ₂| ≤ 0.5·|Δ₁|`. If it holds for Q (and R_in), print the bracket. Otherwise print "**unconverged at h = 0.005 — F31 is the next rung**".
- The order implied by unequal spacing: p solving (h₁ᵖ − h₂ᵖ)/(h₀ᵖ − h₁ᵖ) = ρ on `T30_H_M`, by bisection on [0.1, 6]. It is labelled "printed; not the §10 rule". At this spacing, ρ = 0.5 corresponds to p ≈ 2.3, and first order gives ρ = 1.
- **Predictions, by the drafter, from the logged digits:**
  - R_in: Δ₁ = −128.152430 Ω, Δ₂ = −90.369545 Ω, ρ = 0.70517.
  - Q: Δ₁ = −6.822372, Δ₂ = −4.893386, ρ = 0.71726.
  - f₀: Δ₁ = +25 764.72 Hz, Δ₂ = +24 937.65 Hz, ρ = 0.96790.
  - Branch: **unconverged at h = 0.005** (both ρ > 0.5).
  - p: R_in ≈ 1.68, Q ≈ 1.64, f₀ ≈ 1.06, the last consistent with the T13 ruling's "first order in h".

**Harness (one window).** The Bash call's `timeout` is 660000. The harness call is the whole Bash command: bare and repo-relative, with nothing chained or piped.

`scripts/testing/run_and_log.sh ANS-7-item364-T30 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 30 mpiexec -n 1 python3 -m pytest tests/environment tests/validation/test_ans7_t30_rin_q_sequence.py -v -s --tb=short'"`

Every `FEM_EM_ANS7_EXACT_*` knob stays unset: importing `test_ans7_exact_strip` with one set can raise at fixture time.

**Traps already paid for.**
- **A missed literal shows as a tie failure on a correct measurement** (item 228). Pass `res` and `("cu",)` to `_load_nodes` for the h rungs, and `None` / both columns at ×1 (via `_x1_cu5_rows`).
- **Never `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ`.** This module calls no fixture of `test_ans7_exact_strip`; it only imports helpers.
- **Rows are at `C_F`, untuned and lossless.** They are not the tuned R_in / Q at C_in (item 310's 18 592.58 / 964.43) and not AED's. Never mix the two in the sequence.
- pytest `-s`; harness foreground; no pipe; rule (i): no edit after the window.
- No AED value enters the module or the record. The −1.6 → −2.7 % / −1.9 → −3.1 % gaps in the T13 ruling are relative levels already public there, and this item does not compute them.

**Scope boundary.** Covered: the h-sequence of the Cu R_in / Q / f₀ rows at `C_F` over three logged rungs (degree 2, κ off, exact rectangles, F-small, 64 MHz band). **Not claimed:** any h-converged R_in or Q; any AED agreement (T13's verdict stands); the tuned rows; the copper-loss partition (F31's); 128 MHz.

**What a negative result means.**
- **(F) or (R-rep) red:** the stored JSONs are not what the logs printed. A file was rewritten since, or the path is wrong. Report the per-value residuals and stop, BLOCKED. Never re-solve, and never re-register a logged value.
- **(N-swap) red:** wiring (the wrong files were read). BLOCKED.
- **ρ ≤ 0.5 against the prediction:** the bracket branch fires. Record the bracket on the row, **hold 🟡**, and flag it for the review and auditor. Whether a bracket from this rule is "Q's h-statement" (T13's ✅ condition) is the review's call, not the slot's.
- No `PREDICTED RED:` line.

**Status it can move.** The `ANS-7` row gains a "T30 (item 364)" record clause: "Cu R_in / Q at `C_F`, ×1 → 0.01 → 0.005: successive-change ratios `<ρ_R>` / `<ρ_Q>`; `<branch>`." On the predicted branch, F31 (`xl` entry 34, 10-15) is named as the next rung. **No glyph** on the predicted branch. **Rule (j):** the `ANS-7` Done-when must first carry the T13 ruling (weekly flag 2); see the row amendment recommended with this brief. Against that text, the unconverged branch moves nothing.

**Questions for the review.**
- (364-Q1) The §10 text does not define "Q's bracket". On the ρ ≤ 0.5 branch (predicted not to fire) the brief would print the geometric tail, Q∞ ∈ [v₂ + Δ₂·ρ/(1 − ρ), v₂], ordered. The review confirms that formula or strikes the branch's print.

**Outcome 2026-10-11 (04:30 CDT slot): ✅ done, record.** `tests/validation/test_ans7_t30_rin_q_sequence.py`, `20261011T094910Z_ANS-7-item364-T30.log` (smoke, `-n 1`, 8 passed, 3 s). (F) green on all three rungs (7 / 5 / 5 JSONs, cells 80 798 / 140 270 / 523 880); (R-rep) rel 0.000e+00 on all nine (rung, quantity) pairs; (N-swap) smallest separation 4.935e-3 (h = 0.005 R_in vs h = 0.01 logged), > 1e-7; (A) 0.000e+00. Record: R_in Δ₁ −128.152430 Ω, Δ₂ −90.369545 Ω, ρ 0.705172, p 1.679; Q Δ₁ −6.822372, Δ₂ −4.893386, ρ 0.717256, p 1.645; f₀ Δ₁ +25 764.72 Hz, Δ₂ +24 937.65 Hz, ρ 0.967899, p 1.063 — every prediction met. Branch: **unconverged at h = 0.005 — F31 is the next rung** (R_in and Q). Printed only: 7-node ×1 record vs 5-node 9.173e-10 / 6.549e-10 (as (S5)). Deviation, disclosed: the brief's command ran all of `tests/environment` first; the first window (`20261011T094813Z_ANS-7-item364-T30.log`) overran the 30 s ceiling inside `test_dolfinx_api_migration.py::test_negative_control_detects_reverted_migrations` (a src + tests sweep, before any T30 test ran), so the case was shrunk to `tests/environment/test_complex_mode.py` (the complex guard) and re-run at the same ceiling. No glyph moves.

### §9 item 365 — 2026-10-11 (T34 part 1 — `ANS-7` step 12b: code, regressions, h = 0.01 nodes at 127.74 MHz (1 of 2); chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain T, step T34 (grep "T34. `ANS-7` step 12b"), quoted, not changed: "T28's construction with the 128 MHz knob (T22′'s node set): the in-model zero within the grid's own resolution, circuit = field at `REDUCTION_BAND`; printed C_in(h = 0.01)/C_in(×1) − 1 at 128 MHz beside 64 MHz's +8.347e-4. Target: identity; the second frequency's first h-move. Heavy. Serial on nothing." T29's ruling binds: no third C_in rung at 64 MHz is bought, and this item is the second *frequency*. Exact rectangles, κ off (ruling 2026-10-02). The ESR ruling scopes r to new steps, but T34 is enumerated as T28's construction, so it is **lossless** and labelled lossless. **No band, tolerance, Done-when or target is set or moved.**

**Why two items.** T34 needs seven h = 0.01 node solves at 127.74 MHz before the T28 sequence (C_in, three drives, assembly) can run. None exist: `logs/ans7_exact/t127.74/` holds only ×1 files, checked 2026-10-11. Priced from measured windows, the whole step is ≈ 49 min of compute, more than one 60-minute slot. Item 365 lands the code, proves every existing path unchanged, and solves the first four nodes. Item 365b (below) solves the last three nodes and runs the sequence with **no code edit** (rule (i)). The node JSONs are durable gitignored products, so a park between the two loses nothing.

- **Executor:** implementer. **Edits:** tests only, `tests/validation/test_ans7_exact_strip.py`. No `src/`, `ans:7` or other test-module edit.
- **Tier:** heavy, `-n 12` (item 214's / 310's rung width; ≈ 33–37 GiB summed `ru_maxrss` per window, `…item214-W2-rung.log:2169`), complex build. **Class:** chain step (solve).
- **Cost:** ≈ 52 predicted min: 10 fixed, ≈ 20 code, then W0a ≈ 4.4, W0b ≈ 1.2 and N1 + N2 ≈ 15, plus spare.
  - W0a is item 283 W2, 261 s at `-n 8`.
  - W0b is item 310 W5, 70 s at `-n 2`.
  - N1 and N2 are two h = 0.01 Cu nodes each. Item 214 W2 measured 536 s for nodes (63.2, 64.5) at `-n 12` on a box ≈ 1.4× slow (`20261004T112040Z_ANS-7-item214-W2-rung.log`: Cu sweeps 156.0 / 159.8 s, "solved in 370.8 s"), and item 214 W3 512 s (`…112945Z…`). On 10-08 one h = 0.01 Cu 4×4 took 144.8 s including the mesh (`20261008T111045Z_ANS-7-item310-W1-cin-h001.log:2248`).
  - **If W0a has not started by minute 25, park** (rule (d)).
- **Depends on nothing.** Inputs, all on disk (checked 2026-10-11):
  - `logs/ans7_exact/t127.74/tune_p2.json` (×1 C_in at 127.74 = 1.5274843607471887e-11 F) and `t127.74/tune_inmodel_cin_p2.json`;
  - the seven ×1 `t127.74/node_*_p2_cu.json`;
  - item 214's five `node_*_p2_h0.01_cu.json`, `tune_p2_h0.01.json` and the four `tune_inmodel_cin_p2_h0.01*.json`.
  - If any is absent, stop, BLOCKED. Never re-create one outside this item's own windows.

**Code** (every knob-unset path byte-identical; item 262 / 283 / 310's paths behaviourally identical, shown by W0a / W0b).
1. **`_h_tune()`** (`:1060`) accepts `FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74` with `RESOLUTION=0.01`, `COLUMNS=cu` and `TUNE` ∈ {`nodes`, `cmode`, `inmodel`}. The 64 MHz set stays {`cmode`, `inmodel`}. Every other combination with a resolution raises `ValueError`. (`_tune()` already refuses `nodes` without 127.74.)
2. **`_active()`** (`:958`): in the 127.74 branch, `if _h_tune()` returns `{"target_hz": TUNE128_TARGET_HZ, "tune": TUNE128_MHZ, "nodes": TUNE128_NODES_MHZ, "grid": TUNE128_GRID_MHZ, "out": TUNE128_OUT_DIR, "lab": "T34", "tmhz": "127.74", "res": H_RES, "sfx": H_SFX}`. This is **T22′'s seven-node set**, as §10 says, not a five-node subset. No five-node assembly has been validated at 127.74.
3. **`_node128_json(mhz, cols, res=None)`** (`:972`): with `res` None the name is unchanged. Otherwise it is `t127.74/node_<mhz:.2f>_p2_h0.01_cu.json`.
4. **`tune_nodes`** (`:1176`):
   - `res = H_RES if _h_tune() else None`, then `_fixture_tie(res)` and `_build_ladder(..., resolution=res, ...)`. `resolution=None` is the default (`test_th14_birdcage_copper.py:263`), so the unset path is the same call.
   - The JSON gains `"resolution": res`, written via `_node128_json(mhz, cols, res)`.
   - The prints say "h 0.01" in place of "x1" on the h path only.
   - `test_t22_node_anchors` passes `res` to `_assert_tie`.
5. **`_stored_active()`** (`:1136`) gains a 127.74-h branch: the seven `_node128_json(m, ("cu",), H_RES)` files, each checked as the 64-h branch checks its own (`:1140–1154`): degree 2, κ `None`, `resolution == H_RES`, `field_solve is False`, `columns == ["cu"]`, `target_mhz == 127.74`, one cell count, `f_hz` within 1 Hz, and `_assert_tie(res=H_RES)`. Absent files ⇒ `pytest.fail` naming them.
6. **`tune_cmode`** (`:1386–1401`): the ×1 comparand is read from `act["out"] / f"tune_p{deg}.json"`, which on the 64 MHz path is `OUT_DIR`, unchanged. Prints use `lab`.
   - At 127.74 the C_in h-move is printed beside a new record literal `T28_C_IN_H_REL_64 = 8.346526e-4`, provenance `…item310-W1-cin-h001.log:2234`. Its prediction (below) replaces the 64 MHz text; the 64 MHz print bytes stay.
7. **`_tune_point()`** (`:1076`): `v in _active()["tune"]`.
8. **`_h_point_json(mhz)`** (`:1852`): on the 64 path it stays `OUT_DIR / f"tune_inmodel_cin_p2{H_SFX}_{mhz:.1f}.json"`, byte for byte. At 127.74 it is `TUNE128_OUT_DIR / f"tune_inmodel_cin_p2{H_SFX}_{mhz:.2f}.json"`.
9. **`_h_assemble`, `tune_h`, `test_t28_*`** (`:1856–2056`): every 64 MHz literal is read from `_active()`:
   - `TUNE_MHZ` → `act["tune"]`;
   - `TUNE_TARGET_HZ` → `act["target_hz"]`, including the `tj["f_hz"]` check at `:1953`;
   - `H_NODES_MHZ` → `act["nodes"]`;
   - `GRID_MHZ` → `act["grid"]`;
   - `index(64.0)` and `abs(pt − 64.0)` → the middle tune point `act["tune"][1]`;
   - `T21P_ZERO_RESOLUTION` → `_zero_resolution()`, in the print only;
   - the assembled output → `act["out"] / f"tune_inmodel_cin_p2{H_SFX}.json"`;
   - the "T28" label → `act["lab"]`.
   - **At 127.74 only:**
     - the `STEP6_RUNG_ROWS` comparand print is replaced by the ×1 127.74 rows at C_in, read from `t127.74/tune_inmodel_cin_p2.json` `rows_a` (23 793.40 Ω / 618.34; item 283);
     - (N1‴) becomes **(N0⁗)**: the up-crossings at `C_F` on the h grid, printed and never asserted;
     - (W), item 262's degree-1 64 MHz ×1 drive, runs unchanged in the 127.74 window;
     - `xn = nodes − 64.5` stays, as at both fields today (item 288's code step 4).

**Windows (item 365).** Each runs at `-n` as stated, bare, whole and foreground, with Bash `timeout: 660000`. Every pytest window runs with `-s`, `tests/environment` first, `FEM_EM_REQUIRE_COMPLEX=1` and `FEM_EM_ANS_SKIP_PRIVATE=1`, with `ROUTE` and `FEM_EM_ANS7_EXACT_MHZ` unset. Every window is under 10 min, so §5.1's capture idiom is not required. No edit after W0a starts (rule (i)). Never two heavy windows at once.
0. **Orphan check:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
1. **W0a, ×1 127.74 C_in path unchanged.** Item 283's W2 command verbatim, label `ANS-7-item365-W0a-cin128`, `-n 8`, 261 s (`20261006T171745Z_ANS-7-item283-W2-cin128.log`).
   - It must reproduce to the printed digits (T2″) −2.264621e-6 (0.140 bands, `:2060`), (AB″) 1.590e-5 / 1.601e-5 / 1.616e-5, and r −2.264620e-6.
   - It rewrites `t127.74/tune_inmodel_cin_p2.json`, as item 310's W0 rewrote the ×1 64 MHz file (that file differed at 4e-8 only in residuals). Item 288's (R-T) compares residuals absolutely at 1e-9, so this is harmless.
   - **Red ⇒ revert, BLOCKED.**

   `scripts/testing/run_and_log.sh ANS-7-item365-W0a-cin128 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_C=cin timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
2. **W0b, 64 MHz h path unchanged (pure).** Item 310's W5 command verbatim, label `ANS-7-item365-W0b-assemble64h`, `-n 2`, `timeout -k 30 180`, 70 s (`20261008T112649Z_ANS-7-item310-W5-assemble.log`).
   - It must reproduce (T2‴) \|f₀,A/64 − 1\| 2.545621e-6 (`:138`), r(h) −2.545621e-6 (`:139`), (AB‴) 1.239e-5 / 1.262e-5 / 1.289e-5, and C_in(h)/C_in(×1) − 1 = +8.346526e-4.
   - **Red ⇒ revert, BLOCKED.**

   `scripts/testing/run_and_log.sh ANS-7-item365-W0b-assemble64h "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_TUNE_POINT=assemble timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
3. **N1 (nodes 125.74, 126.04):**

   `scripts/testing/run_and_log.sh ANS-7-item365-N1-nodes-h001 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=nodes FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_NODES=125.74,126.04 timeout -k 30 590 mpiexec -n 12 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
4. **N2 (nodes 126.94, 128.24):** the same command with `FEM_EM_ANS7_EXACT_TUNE_NODES=126.94,128.24`, label `ANS-7-item365-N2-nodes-h001`.
   - **Pre-registered contingency, not a timeout change:** if N1's footer exceeds 480 s, N2 and item 365b's nodes run **one node per window** (labels `…-N2a` / `…-N2b`, …).
5. **Commit** the module, the W0a / W0b / N1 / N2 logs, their test-results rows (appended by hand) and the §9 marking together. Item 365 is DONE as "code + regressions + 4 of 7 nodes". Item 365b becomes runnable.

**Anchors (asserted; every band imported by name).**
- **W0a / W0b:** the existing tests, green, with the readings above reproduced to the printed digits.
- **N1 / N2 (`test_t22_node_anchors`, per node):**
  - reciprocity ≤ `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py`);
  - σ_max ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE` (`tests/validation/test_port_birdcage_four_port.py`);
  - (F) `_assert_tie(cells, tie, res=0.01)`, with cells equal to item 214's 140 270 within `CELL_COUNT_BAND`.

**Negative controls (rule (e)).** W0a's (N0′) is asserted, backed by item 283 W2 (`n0_crossings_up_hz: []`). W0b's (N1‴) and (N2‴) are printed, as item 310 W5 printed them.

**Predictions (printed).**
- N1 / N2: Im `Y_in` at `C_F` at each node > 0, as on the ×1 grid (item 228's "PREDICTED > 0").
- Cells 140 270; summed `ru_maxrss` 33–37 GiB; each window 380–540 s.

**Outcome (2026-10-11 04:30 slot): ✅ code + regressions + 4 of 7 nodes.** Code steps 1–9 landed in `tests/validation/test_ans7_exact_strip.py` before W0a (rule (i); no edit after). W0a (`20261011T095408Z_ANS-7-item365-W0a-cin128.log`, `-n 8`, 261 s, 16 passed): (T2″) −2.264621e-6 (0.140 bands), (AB″) 1.590e-5 / … / 1.616e-5, r −2.264620e-6, (N0′) `[]` — reproduced. W0b (`20261011T095838Z_ANS-7-item365-W0b-assemble64h.log`, `-n 2`, 79 s, 15 passed): (T2‴) 2.545621e-6 (0.198 × T21P), r(h) −2.545621e-6, (AB‴) 1.239e-5 / 1.262e-5 / 1.289e-5, C_in(h)/C_in(×1) − 1 = +8.346526e-4 — reproduced byte for byte. N1 (`20261011T100007Z_ANS-7-item365-N1-nodes-h001.log`, 378 s, 13 passed) and N2 (`20261011T100632Z_ANS-7-item365-N2-nodes-h001.log`, 377 s, 13 passed; footers under 480 s, so the one-node contingency is not triggered for 365b): cells 140 270, Im `Y_in` at `C_F` +4.0553e-2 / +4.0717e-2 / +4.1207e-2 / +4.1910e-2 S at 125.74 / 126.04 / 126.94 / 128.24 MHz (PREDICTED > 0: met ×4), reciprocity ≤ 3.3e-15, σ_max ≤ 0.999783, summed `ru_maxrss` 36.74 / 36.62 GiB. Four `t127.74/node_*_p2_h0.01_cu.json` written. Item 365b is runnable.

### §9 item 365b — 2026-10-11 (T34 part 2 — `ANS-7` step 12b: nodes 5–7, C_in(h), three drives, assembly; chain step (solve))

*Same licence, executor, file and tier as item 365.* **Depends on item 365 landing; if it did not, skip.** **No code edit:** the module runs as item 365 committed it (rule (i)). If an edit proves necessary, re-run W0a / W0b before any further window and say so in the commit.
- **Cost:** ≈ 45 predicted min: 10 fixed, then N3 ≈ 7.5, N4 ≈ 4.2, C1 ≈ 4.2, D1–D3 ≈ 11.3 and A ≈ 1.2, plus spare. C1 is priced from item 310 W1 (252 s); D1–D3 from item 310 W2 / W3 / W4 (216 / 236 / 214 s); A from item 310 W5 (70 s).
- **Inputs:** item 365's four `t127.74/node_{125.74,126.04,126.94,128.24}_p2_h0.01_cu.json`. If they are absent, BLOCKED.

**Windows** (same rules as item 365; orphan check first; one heavy window at a time).
1. **N3:** nodes `129.54,130.44`, label `ANS-7-item365b-N3-nodes-h001`, item 365's N1 command otherwise.
2. **N4:** node `130.74`, label `ANS-7-item365b-N4-nodes-h001`.
3. **C1, C_in(h) at 127.74:**

   `scripts/testing/run_and_log.sh ANS-7-item365b-C1-cin128-h001 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=cmode FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 timeout -k 30 590 mpiexec -n 12 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`

   It writes `t127.74/tune_p2_h0.01.json` and never touches the ×1 `t127.74/tune_p2.json`.
4. **D1 / D2 / D3, one Route A drive each:** C1's command with `FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_TUNE_POINT=127.64` (then `127.74`, then `127.84`). Labels `ANS-7-item365b-D1-p127.64` / `-D2-p127.74` / `-D3-p127.84`. Run them serially; (W) runs in D2.
5. **A, assembly (closing):** `FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_TUNE_POINT=assemble`, `-n 2`, `timeout -k 30 180`, label `ANS-7-item365b-A-assemble`.

**Anchors (asserted; bands imported).**
- **N3 / N4:** as item 365's N windows.
- **C1:**
  - (C) Im z₁ > 0, reciprocity and passivity, with cells equal to the stored node cells;
  - (T1) at C_mode,Cu(h): one − → + crossing on `TUNE128_GRID_MHZ`, k-share > `K_SHARE_FLOOR`, and Δ_LOO over the **seven** h nodes ≤ `BANDS`;
  - (N1) the `c_f` wiring on item 196's stored 64 MHz ×1 nodes, ≤ `NODE_S_RTOL` and > `N_FACTOR` × `BANDS`. These are existing assertions, backed by items 217 / 228 / 310 on the same comparison.
- **D1–D3:** (AB⁗) \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` (1e-3, `tests/validation/test_port_lumped_rlc_termination.py:88`). (F) the drive mesh equals the node files' cells and is off the ×1 record (N-wire). **(W)** ≤ `NODE_S_RTOL` in D2.
- **A:** (AB⁗) at all three points. **(T2⁗)** Im `Y_in`,A changes sign − → + exactly once in (127.64, 127.84) MHz and never + → −, with f₀,A inside the bracket. This is 310-Q1's form.

**"Within the grid's own resolution": the same gap as 310-Q1 (365-Q1 below).** T22P_ZERO_RESOLUTION (1.616e-5) was registered from the ×1 grid's (AB) residuals at C_mode (`test_ans7_exact_strip.py:1011–1017`). The h = 0.01 node files carry no identity residuals (`field_solve` false), so nothing exists to register the h-rung's own band from before the drives. The brief therefore asserts the crossing inside the bracket, and prints \|f₀,A/127.74 MHz − 1\| beside `T22P_ZERO_RESOLUTION` **for scale only, never asserted**. No band is transplanted.

**Negative controls (rule (e)).**
- **(W): asserted.** Backed by item 262 W1 (4.034e-14), item 283 (1.470e-13) and item 310 W3 (1.720e-13): same ×1 degree-1 drive, same comparison.
- **(N0⁗): predicted, printed.** At `C_F` the h 127.74 grid has no − → + up-crossing. It is asserted at ×1 as (N0′) but never measured on this grid.
- **(N2⁗): predicted, printed.** The circuit at 0.5·C_in(h) against the in-model at C_in(h), at 127.74, is ≥ 10 × `REDUCTION_BAND` (×1 127.74: 493.82×).
- **(N-h): predicted, printed.** \|C_in(h)/C_in(×1) − 1\| ≥ 10 × the secant tolerance.

**Predictions (written first, printed).**
- **C_in(h)/C_in(×1) − 1 at 127.74 in [+2e-4, +2e-3], sign positive.** It is printed beside 64 MHz's +8.346526e-4. At 64 MHz the f₀ h-move at fixed C was +0.0411 %, and tuning back needed more C. If the 127.74 MHz f₀ moves the same way, the sign is the same. Item 310's brief got this sign wrong; the reasoning is stated so the miss, if any, is legible.
- C_mode,Cu(h) in [14.8, 15.6] pF (`_PRED["127.74"]`).
- (T1) Δ_LOO ≤ 1e-7.
- (AB⁗) in [8e-6, 4e-5] (×1 127.74: 1.590–1.616e-5; h at 64 MHz moved +0.2 %).
- r(h) in [−5e-6, −1e-6] (×1 127.74: −2.265e-6).
- \|f₀,A/127.74 − 1\| ≈ 2.3e-6, about 0.14 × T22P (printed).
- R_in / Q at C_in(h) within ±2 % of ×1 127.74's 23 793.40 Ω / 618.34. At 64 MHz, h = 0.01 moved them −0.78 % / −0.68 % against ×1 at C_in.

**Traps already paid for** (both items).
- Items 217 / 228 / 262 / 283 / 310's lists, whole:
  - κ `None`;
  - no `ans:7` call;
  - `z_from_s(S)`, never `z_matrix`;
  - a real-scalar secant;
  - rank 0 writes, then a barrier;
  - a killed window is not a result: run `pgrep -c python3`, report, never lengthen the timeout;
  - never `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ` (`_nodes_env` validates against the 64 MHz nodes).
- **A missed 64 MHz literal in the 127.74 path** shows as a bracket, `f_hz` or tie failure on a correct measurement (item 228 code step 4). Code step 9 lists them.
- **`.1f` file names collide at 127.74.** 127.64 / 127.74 / 127.84 would be written as "127.6" / "127.7" / "127.8". Use `.2f` on the 127.74 path, and keep `.1f` on the 64 MHz path (W0b reads those names).
- **Item 157's trap:** `pec_column=False` with `field_solve=True` raises on rank 0 only. The h path runs `pec_column=False`, `field_solve=False`, which is outside the trap (item 310's 310-Q3 reading of `test_th14_birdcage_copper.py:357–412`).
- **Never write** the ×1 `t127.74/tune_p2.json` / `node_*_p2_cu.json`, or item 214's / 310's 64 MHz h files. W0a / W0b's rewrites of their own outputs are the accepted exception, as item 310 W0's was.
- Price on the measured node time (item 214 ran ≈ 1.4× slow), not on 10-08's faster single-node figure.
- Exact rectangles only. Never a sub-sheet out of a box.

**Scope** (both items). Covered: h = 0.01, degree 2, κ off, Cu, 127.74 MHz, lossless C, T22′'s seven nodes. **Not claimed:**
- any absolute tuned C, f₀, R_in or Q;
- an h-statement at 127.74 (one rung; T29's record is the weekly's to extend);
- ESR;
- any AED or Palace comparison.

**Status it can move.** **No glyph.** The `ANS-7` row gains a "step 12b (T34, items 365 / 365b)" record clause: C_in(h)/C_in(×1) − 1 at 127.74 beside +8.347e-4 at 64 MHz, (AB⁗), and the zero's printed offset. **Rule (j):** declined, as for items 228 / 262 / 283 / 288 / 310. The Q clause is T30 / F31's, and this step moves no glyph.

**What a negative result means.**
- **W0a / W0b red:** the edit changed an existing path. Revert it, BLOCKED.
- **A node window killed:** commit nothing for that node. The next slot reruns it alone, one node per window, with no edit.
- **(T1) Δ_LOO red with the nodes green:** the seven-node h grid does not resolve Route B at 127.74 within `BANDS`. Print Δ_LOO, open a known-issues entry, park and mark BLOCKED. Never drop to five nodes, and never re-band.
- **(AB⁗) red with (W) green:** circuit and field disagree on the rung at 127.74. Print the residuals, open a known-issues entry, park and mark BLOCKED. Never apply κ.
- **No crossing in (127.64, 127.84):** BLOCKED. Report the Im `Y_in` triple and C_in(h).
- No `PREDICTED RED:` line.

**Questions for the review.**
- (365-Q1) T34's "within the grid's own resolution" cannot be pre-registered on the h rung (the gap 310-Q1 named for the 10-10 weekly, which did not rule on it). The brief applies 310-Q1's treatment: crossing asserted, offset printed for scale. The review confirms, or names the band's source.
- (365-Q2) The split into 365 / 365b is the drafter's. Item 310 ran its six windows in one slot, but needed no node solves. The review may merge them if it prefers one item that parks between N2 and N3.

**Outcome (2026-10-11 06:00 slot): ✅ (record, lossless; no glyph).** No code edit (rule (i)). N3 (`20261011T110054Z_ANS-7-item365b-N3-nodes-h001.log`, `-n 12`, 358 s, 13 passed): Im `Y_in` at `C_F` +4.2610e-2 / +4.3092e-2 S at 129.54 / 130.44 MHz, recip ≤ 2.4e-15, σ_max ≤ 0.999777, 35.69 GiB. N4 (`20261011T110701Z_ANS-7-item365b-N4-nodes-h001.log`, 259 s, 13 passed): +4.3252e-2 S at 130.74 MHz, 33.76 GiB; cells 140 270 on all three (PREDICTED > 0: met ×3). C1 (`20261011T111130Z_ANS-7-item365b-C1-cin128-h001.log`, 260 s, 15 passed): C_mode,Cu(h) 15.288627 pF (PREDICTED [14.8, 15.6]: met); C_in(h) = 1.5288541964455146e-11 F; **C_in(h)/C_in(×1) − 1 = +8.967920e-4** (PREDICTED [+2e-4, +2e-3], positive: met) beside 64 MHz's +8.346526e-4; (T1) one up-crossing, k-share 0.999987, Δ_LOO 1.8e-13 / 2.6e-10 / 2.8e-10 (PREDICTED ≤ 1e-7: met); (N1) ≤ 2.2e-12; (N0) `[]`; (N-h) 8.968e-4 vs 8.1e-13 (met). Two inherited ×1-text prints in C1 read "NOT met" — (C_mode − C_in)/C_in = +5.548e-6 (printed prediction negative) and f₀,B(C_mode)/127.74 − 1 = −2.616e-6 (printed \|.\| ∈ [3e-6, 5e-5]); neither is this brief's prediction nor asserted, both recorded. D1 / D2 / D3 (`20261011T111601Z…D1-p127.64`, `…111940Z…D2-p127.74`, `…112349Z…D3-p127.84`; 209 / 239 / 215 s, 14 passed each): (AB⁗) 1.593e-5 / 1.605e-5 / 1.619e-5 (PREDICTED [8e-6, 4e-5]: met); (W) 1.454e-13 ≤ `NODE_S_RTOL`. A (`20261011T112734Z_ANS-7-item365b-A-assemble.log`, `-n 2`, 80 s, 15 passed): (T2⁗) Im `Y_in`,A −4.0615e-5 / +1.1818e-7 / +4.0889e-5 S, one − → +, none + → −; f₀,A 127.739710 MHz, \|f₀,A/127.74 − 1\| = 2.271224e-6 (0.141 × `T22P_ZERO_RESOLUTION`, printed only); r(h) −2.271224e-6 (PREDICTED [−5e-6, −1e-6]: met); R_in / Q at C_in(h) 23 319.94 Ω / 606.70, −1.99 % / −1.88 % against ×1 (PREDICTED within ±2 %: met, narrowly); (N0⁗) `[]`; (N2⁗) 494.00 × `REDUCTION_BAND`. One rung at 127.74, not an h-statement.

### §9 item 378 — 2026-10-11 (T27 — `ANS-7` step 13: `ans:7` gains the ESR rows at 64 MHz; example)

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain T, step T27 (grep "T27.\*\* `ans:7` gains the ESR rows"), quoted, not changed: "`ans:7` gains the ESR rows (T26 at 64 MHz; lossless rows stay, labelled). `example-runner`. Serial on T26." The ESR ruling (§9, 2026-10-06) sets r = 0.2 Ω on every capacitor, both sides. Lossless records stay valid and labelled, and r never re-reads a lossless record. **No band, tolerance, Done-when, SPEC or validation target is set or moved; every band is imported by name.**

**Housing.** No EX chunk owns `ans:7`: the §7 EX table has no `ans:7` row, and `examples/ansys_benchmarks/README.md` lists the case under **ANS-7**. Items 230 (T20, step 9) and 288 (T23, step 10) set the precedent: `example-runner` re-issues `ans:7` as an `ANS-7` step, briefed here. T27 is **`ANS-7` step 13**; T26 is step 11 and T28 / T34 are steps 12 / 12b.

**Depends on item 311 (T26) landing; if it did not, skip and leave 378 ⬜.** Item 311 writes the stored inputs this stage reproduces: `logs/ans7_exact/tune_p2_esr0.2.json` (`c_in_r`, `c_in_0`, `esr_ohm`, `f_hz`) and `tune_inmodel_cin_esr0.2_p2.json`. Neither exists today (checked 2026-10-11). If either is absent when the slot opens, stop and take the next item. Never re-create one; that is the test module's job.

`ANS-7` step 13: `ans:7` gains one additive stage, **`XT64R`**. It re-drives Route A at degree 2 on the ×1 exact-sheet mesh at (r = 0.2 Ω, C_in(r)) at 63.9 / 64.0 / 64.1 MHz, with the series-RC sheet on P2–P4 and the feed gap's (r + C) branch in parallel at P1. It re-asserts T26's gates through the example and writes the ESR rows into `metrics.json` (`degree2_exact["tuned_larmor_esr"]["64.0"]`) and `COMPARISON.md`. The lossless `tuned_larmor` rows stay byte-identical and are labelled lossless.
- **Executor:** `example-runner`. Spawn prompt: "Use the example-runner agent to execute `ANS-7` step 13 (§9 item 378) per `docs/planning/chunks/ANS-7.md` § '§9 item 378'. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run every harness window foreground with Bash `timeout: 660000`. Never background a window or end your turn while one is running. Set `FEM_EM_ANS_SKIP_PRIVATE=1` on every window that runs the `ans:7` script. Never run `ans:1`–`ans:3`, never run B1 / B2 / A / X1–X4 / XA / XT64 / XT128, and never use `./scripts/run_examples.sh -e ans`."
- **Edits:** the example script, its same-stem guide, `metrics.json` and `COMPARISON.md` only. No `src/`, test-module, SPEC or `aed/` edit.
- **Standing rule (a) applies.** `z_cap`, `add_c1`, `route_b_z` and `_route_a_solve` are imported by gate modules: `test_ans7_exact_strip.py:83–90` and `test_ans7_sparse_route_b.py` via `_ANS7`. So the keyword is additive with a bitwise default, disclosed in the commit, and the gate module is re-run green in the same slot by window **W-gate** below.
- **Tier:** heavy, `-n 8` (T21′ / T26's record width), complex build. **Class:** example.
- **Cost:** ≈ 8 measured (class median) / ≈ 45 predicted min: 15 fixed, ≈ 15 code and guide, `XT64R` ≈ 3–4.5, ≈ 3 census / pin, ≈ 8 spare. Price: item 288's `XT64`, 154 s (mesh 24.8 s, drives 39.0 / 43.6 / 42.8 s; `20261007T003240Z_ANS-7-item288-XT64.log`), which has the same drive count. The series-RC sheet changes only a port impedance, and item 309 measured no cost difference between its r = 0.2 and r = 0 arms beyond noise (94 / 140 s, `PORT-24.md`). Predicted 150–260 s, under 590 s on a ×1.4-slow box (≈ 370 s).

**Code** (`examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py`; additive; every existing stage byte-unchanged in behaviour).
1. **The SPEC circuit layer gains `esr_ohm=0.0`**, additive, on `z_cap` (`:167`), `add_c1` (`:171`), `route_b_z` (`:181`) and `_route_a_solve` (`:541`).
   - With r = 0 each keeps its **old expression by a branch, not arithmetic**: `1/(0 + 1/(jωC))` is not bitwise `jωC` (item 308's trap).
   - With r ≠ 0, the gap impedance is `series_rlc_impedance(f, r_ohm=r, l_h=0.0, c_f=C)` (`src/fem_em_solver/ports/lumped.py:138`), imported, never `r + 1/(jωC)` restated.
   - `add_c1`'s feed branch is `1/(1/z_bare + 1/Z_RC)`.
   - `route_b_s` is untouched; it is not on the XT path.
   - **Never enter 0.2 Ω as a parallel R** (ESR ruling (iii)).
2. **Stage `XT64R`** is added to `main()`'s accepted set (`:1503`) and the error message. It starts with `X`, so `_x_guard()` applies: `FEM_EM_ANS7_EXACT_ESR_OHM` and every other `FEM_EM_ANS7_EXACT_*` knob must be **unset**. The stage selects by name alone (item 288's knob-free trap) and dispatches to `_xt_stage` or a sibling `_xt_esr_stage`, whichever keeps XT64 / XT128 byte-unchanged.
3. **Lazy imports inside the stage** (item 230's cycle trap), imported, never restated:
   - item 288's list (`:1293–1298`);
   - whatever constant item 311 lands for r (its landing commit names it). If item 311 lands none, read `esr_ohm` from `tune_p2_esr0.2.json` and assert it equals 0.2, the operator ruling's value, never a band.
4. **Inputs.**
   - `OUT_DIR / "tune_p2_esr0.2.json"`: `c_in_r`, `c_in_0`, `esr_ohm`, `f_hz` within 1 Hz of 64 MHz. **C_in(r) is never recomputed.**
   - Route B nodes come from `_stored_x1_p2()` (Cu), T26's inputs.
   - `OUT_DIR / "tune_inmodel_cin_esr0.2_p2.json"` is the (R-T-r) comparand.
5. **Drives:** item 288's loop at `c_f = c_in_r` and `esr_ohm = r`.
   - `s_a = _route_a_solve(th14, geo, specs, g·1e6, "cu", c_f=c_in_r, degree=2, esr_ohm=r)`;
   - `s_b = _s11_from_z(route_b_z(bary_interp(xn, z_cu, g − 64.5), f, c_in_r, esr_ohm=r)[1])`;
   - Z_in,A = `add_c1(_z_from_s11(s_a), ω, c_in_r, esr_ohm=r)`;
   - rows from `parallel_rows`;
   - f₀,B(C_in(r)): `sparse_rows` takes no `esr_ohm` (`test_ans7_sparse_route_b.py:321`, a test module this item does not edit). So f₀,B comes from the Route B Im `Y_in` over the 51-point `GRID_MHZ` through the keyworded `route_b_z` and the example's `parallel_rows`, exactly as item 311 forms it. Copy 311's construction and cite it; never invent one.
   - **Printed, never asserted:** Q(r)/Q_0 against the closed form 1/Q = 1/Q_0 + 1/Q_C, Q_C = 1/(ωC_in(r)·r). Q_0 is the lossless `tuned_larmor["64.0"]` `rows_a.q_loaded`, the coil's loaded Q with lossless capacitors. Also printed: C_in(r)/C_in(0) − 1.
6. **Write** (rank 0, then a barrier).
   - Snapshot every top-level key and every `degree2_exact` key **including `tuned_larmor`** (it must stay byte-identical).
   - Set `degree2_exact["tuned_larmor_esr"]["64.0"]` = {`esr_ohm`, `target_hz`, `c_in_r_f`, `c_in_0_f`, `c_in_rel`, `tune_mhz`, `s11_a`, `s11_b`, `residual_ab`, `reduction_band`, `rows_a`, `f0a_rel_target`, `zero_resolution` (value and name), `bands_fraction`, `r_ab`, `q_ratio`, `q_closed_form`, `lossless_comparand: "tuned_larmor['64.0']"`, `rt`, `cells`, `mpi_ranks`, `t_mesh_s`, `t_drive_s`, `source` (item 311's log names)}.
   - Assert (K), and revert `metrics.json` on a miss.
7. **`COMPARISON.md`:**
   - The T23 section's heading gains "— lossless capacitors (r = 0)". That is the only change to existing text, and it is disclosed in the commit.
   - A fourth section, "Tuned with capacitor ESR (r = 0.2 Ω on every capacitor; §10 T26 / T27)", has one 64.000 MHz column with rows C_in(r) [pF], C_in(r)/C_in(0) − 1, f₀,A [MHz], f₀,A/64 − 1, bands used, (AB) max, R_in, Q(r), Q(r)/Q_0 and the closed-form Q.
   - **There is no AED column.** No AED ESR case exists.
8. **Docs.**
   - The docstring gets an `XT64R` paragraph.
   - The guide gets a section after "Tuned to the Larmor frequency (current; §10 T23)": "Capacitor ESR (current; §10 T27)", with the stage, its command and its anchors.
   - "Not claimed" gains "a measured ESR or its frequency dependence; ESR rows are ×1-mesh degree-2 records".
   - No setup figure or XDMF (item 288's reading of the `EX-93` conventions). Do not set `FEM_EM_SETUP_FIGURES`.

**Anchors** (asserted in-stage; every band imported):
- **(AB-r)** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` at all three points.
- **(T2-r)** Exactly one − → + crossing of Im `Y_in`,A in (63.9, 64.1), and \|f₀,A/64 MHz − 1\| ≤ `T21P_ZERO_RESOLUTION`. These are T26's anchors, as item 311 asserts them.
- **(R-T-r) reproduction of item 311's registered run**, against `tune_inmodel_cin_esr0.2_p2.json`, applied as item 288's (R-T):
  - \|f₀,A,ex/f₀,A,stored − 1\| ≤ `NODE_S_RTOL`;
  - max \|res_ex − res_stored\| ≤ `NODE_S_RTOL`;
  - max \|S11_A,ex − S11_A,stored\|/\|S11_A,stored\| ≤ `NODE_S_RTOL`;
  - the stored `ranks == 8` and `degree == 2`.
  - At any other width it prints and does not assert (the `EX-56` rule). This anchor shows that the example's lifted ESR terminations are item 311's test-side ones.
- **(K)** Every pre-existing `metrics.json` key, `tuned_larmor` included, is byte-identical.
- **(R0) the r = 0 branch is bitwise:** `z_cap`, `add_c1` and `route_b_z` with `esr_ohm=0.0` equal the unkeyworded call, bit for bit, on the stored ×1 `z64` at 64.0 MHz. This is pure and takes milliseconds.

**Negative controls (rule (e)).**
- **(N-L): predicted, printed.** The lossless Route B at C_in(0) against the in-model ESR S11 at 64.0 MHz is ≥ 10 × `REDUCTION_BAND`. **If item 311 measured (N-L) on the same fixture, it is asserted here, citing item 311's log line; otherwise it is printed.** The executor reads 311's outcome before the window and labels it.
- **(N-C-r): asserted.** The (R-T-r) comparator separates the ESR run from the lossless run: \|f₀,A,ex(r)/f₀,A,stored(lossless C_in) − 1\| > 100 × `NODE_S_RTOL`, against `tune_inmodel_cin_p2.json` `rows_a.f0_hz`. **Asserted only if item 311's outcome records f₀,A(r) and f₀,A(0) differing by > 1e-7. Otherwise it is printed.** It is predicted O(1/Q_C²) ≈ 3e-5 (311's brief).

**Predictions (printed).** (R-T-r) ≤ 1e-12 (item 288 reproduced item 262 to 9.97e-13 at `-n 8`). Everything else is item 311's readings, to the printed digits.

**Windows** (`-n 8`; each bare, whole and foreground; Bash `timeout: 660000`; no edit after `XT64R` starts, rule (i)). The window is under 10 min, so §5.1 capture is not required.
0. Orphan check: `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
1. **Censuses before**, labels `ANS-7-item378-pre-figures` / `-pre-docrefs`, item 288's commands (`check_example_setup_figures.py`, `check_example_doc_references.py`, `timeout -k 30 60`).
1b. **W-gate (rule (a)), pure:** item 310's W5 command verbatim (`ANS-7.md` § "§9 item 310"; `20261008T112649Z_ANS-7-item310-W5-assemble.log`, `-n 2`, `timeout -k 30 180`, 70 s), label `ANS-7-item378-W-gate`.
   - It drives `route_b_z`, `z_cap`, `_s11_from_z`, `parallel_rows` and `sparse_rows` through the gate module at the default r.
   - It must reproduce (T2‴) 2.545621e-6, r(h) −2.545621e-6 and (AB‴) 1.239e-5 / 1.262e-5 / 1.289e-5 to the printed digits.
   - **Red ⇒ revert the keyword, BLOCKED.**
   - `_route_a_solve`'s r = 0 branch is covered by (R0): factor the gap impedance into one helper used by both `_route_a_solve` and `route_b_z`, and (R0) checks that helper bitwise.
2. **XT64R:**

   `scripts/testing/run_and_log.sh ANS-7-item378-XT64R "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_STAGE=XT64R timeout -k 30 590 mpiexec -n 8 python3 examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py'"`
3. **Pin and unit tests** on the committed `metrics.json`, `-n 1`, label `ANS-7-item378-pin`, item 288's command (`tests/unit/test_doc_reference_exit_codes.py tests/unit/test_ans7_private_comparison.py`, `-s`). Predicted 24 passed. `test_ans7_private_comparison.py` is the tracked-output guard: it must stay green with the new key.
4. **Censuses after**, labels `-post-figures` / `-post-docrefs`. Predicted equal to the before readings (`broken=0`; docrefs `dead=0`, `guide=0`, exit ≠ 1).
5. Commit the script, guide, `metrics.json`, `COMPARISON.md`, logs and test-results rows (appended by hand) together. `git diff --stat` must show no `COMPARISON_private.md`, `aed_results/`, `SPEC.md` or test-module change.

**Traps already paid for.**
- Items 230 / 288's lists, whole:
  - lazy imports;
  - never run B2 (`_write_outputs` rewrites wholesale);
  - never rename or move the degree-1 keys;
  - `z_from_s(S)`;
  - κ `None`;
  - AED values never enter a tracked file;
  - a killed window is not a result;
  - bare harness calls.
- **Knob-free selection:** setting `FEM_EM_ANS7_EXACT_ESR_OHM` to reach the ESR set trips `_x_guard`. Select by stage name.
- **`FEM_EM_ANS_SKIP_PRIVATE=1` on every window that runs the script.** `_write_private_if_available` honours it (`:511`). Without it the private writer reads `aed_results/`. No slot opens `docs/private/` or `aed_results/`.
- **The r = 0 path is a branch** (item 308). **Never enter 0.2 Ω as a parallel R** (ESR ruling (iii)). **The feed's branch carries r too** (ruling 311-Q5).
- Docrefs exit 1 on a bare file name in the guide (item 273): name stored inputs as "the gitignored `logs/ans7_exact/` stage files". Write ≈, never `~`.

**Scope.** Covered: the example re-asserting T26 at (r, C_in(r)) on the ×1 exact-sheet mesh, degree 2, κ off, Cu, 64.000 MHz, and recording the ESR rows beside the labelled lossless rows. **Not claimed:**
- a measured ESR or its frequency dependence;
- 127.74 MHz with ESR (not enumerated);
- any absolute tuned C, f₀, R_in or Q;
- any h-statement;
- any AED or Palace comparison;
- the private writer's rows.

**Status it can move.** §10 T27 closes, and the weekly removes it. The `ANS-7` row gains a "step 13 (T27)" clause: `ans:7`'s tracked outputs carry the 64 MHz ESR rows, with (AB-r), (T2-r) and (R-T-r) re-asserted through the example and the lossless rows labelled. **No glyph.** **Rule (j):** declined, as for item 288.

**What a negative result means.**
- **(R-T-r) red at `-n 8` with the bands green:** the example's ESR path is not item 311's. Print the differences, revert `metrics.json` / `COMPARISON.md` (those two paths only), park on `attempt/ANS-7-step13-<UTC>` and mark BLOCKED. Never re-register.
- **(R0) red:** the r = 0 branch is not bitwise. Revert, BLOCKED.
- **(AB-r) / (T2-r) red:** this contradicts item 311's green. Report, open a known-issues entry if it reproduces, BLOCKED. Never re-band or adjust r.
- **(K) red:** revert both tracked outputs, BLOCKED.
- No `PREDICTED RED:` line.

**Review ruling 2026-10-11 18:00 (daily) on item 378's (N-L).** Item 311 measured (N-L) on this fixture and it **missed** its ≥ 10× bar (4.97×, the 10:30 review's "printed predictions missed"). The brief's condition "if item 311 measured (N-L) … it is asserted here" therefore reads **"measured and met"**: (N-L) stays *printed* in this item, beside 311's 4.97×, and is never asserted (rule (e)). No band moves.

**Outcome 2026-10-12 (00:00 CDT slot, example-runner): DONE.** `XT64R` at `-n 8`, 80 798 cells, 187 s window (mesh 28.3 s, drives 47.5 / 53.9 / 54.3 s; `20261012T050520Z_ANS-7-item378-XT64R.log`): (R0) `z_cap` / `_gap_z` / `add_c1` / `route_b_z` bitwise at r = 0; (AB-r) 1.233e-5 / 1.258e-5 / 1.283e-5 ≤ `REDUCTION_BAND`; (T2-r) one − → + crossing, f₀,A/64 MHz − 1 = −2.583111e-6 (0.201 × `T21P_ZERO_RESOLUTION`); (R-T-r) \|f₀ ratio − 1\| 4.330e-15, max \|res diff\| 4.150e-13, max S11_A rel 5.182e-13 ≤ `NODE_S_RTOL`; (K) 26 + 19 + 2 keys identical. Printed: (N-L) 4.966e-3 = 4.97 × `REDUCTION_BAND` (≥ 10× not met, as 311); (N-C-r) 4.727e-8 = 47.3 × `NODE_S_RTOL`, not asserted since 311's stored separation is 4.73e-8 < 1e-7 (100× bar not met); Q(r)/Q_0 0.169446 vs closed form 0.167603 (+1.10 %); C_in(r)/C_in(0) − 1 = +4.602610e-5; f₀,B/64 MHz − 1 = +1.89e-12. W-gate 15 passed / 26 skipped 97 s reproduces (AB‴) 1.239e-5 / 1.262e-5 / 1.289e-5 and (T2‴) 2.545621e-6 (`20261012T050333Z_ANS-7-item378-W-gate.log`); pin 24 passed (`…050836Z_ANS-7-item378-pin.log`); censuses before = after (`missing=0 broken=0`; docrefs `dead=0 guide=0 stale=84 exit=2`). Deviation: f₀,B uses item 311's three-point `parallel_rows` construction (`tune_esr` in `test_ans7_exact_strip.py`), not the 51-point `GRID_MHZ` the brief named — the brief's "copy 311's construction" governs. r is read from the stored JSON and asserted equal to 0.2 and to `ESR_VALUES[0]` (311 landed no example-side constant).

### Row text moved 2026-10-11 (03:00 daily review, §7 row ≤ 2 KB; replaced on the row by the T13 / T29 state line and the weekly's re-stated Done-when)

Step 12 (T28, §9 item 310, 2026-10-08): on step 6's h = 0.01 degree-2 rung (140 270 cells) C_in(h)/C_in(×1) − 1 = +8.347e-4 (lossless; −2× step 6's +0.0411 % f₀ move is 8.22e-4 in size); (AB‴) circuit = field 1.239e-5 / 1.262e-5 / 1.289e-5 at `REDUCTION_BAND`; (T2‴) Im `Y_in`,A changes sign − → + once in (63.9, 64.1) MHz, the zero printed at \|f₀,A/64 MHz − 1\| = 2.546e-6 (0.198 × `T21P_ZERO_RESOLUTION`, for scale only, never asserted); one rung, not an h-statement (`20261008T111045Z_ANS-7-item310-W1-cin-h001.log`, `20261008T112649Z_ANS-7-item310-W5-assemble.log`). **Done-when:** The row stays 🟡 until §10 T11 / T12 supply the order- and h-matched rung and the weekly re-reads (T13). 

**Review ruling 2026-10-11 (03:00 daily) on the weekly's T13 / T29 rulings (flag 2).** The row now carries: T13 decided — `f₀` AGREE with its h-statement (0.11 % at ×1, 0.03 % at h = 0.005, first order in h), `R_in` / Q moving away (Q −3.1 %, unconverged); T29's C_in record (+0.083 % at h = 0.01 measured, ≈ +0.16 % at h = 0.005 inferred); Done-when "✅ when T30 / F31 write Q's h-statement; 🟡 on the Q clause until then". Item 364 (T30) — 364-Q1 ruled: on the ≤ ½ branch the bracket is written in entry 34's form, [x(0.005) − \|Δ₂\|, x(0.005) + \|Δ₂\|], no rate; the drafter's geometric tail is struck. Item 365 / 365b (T34) split accepted (no h = 0.01 node exists at 127.74 MHz); 365-Q1: the crossing is asserted and the offset printed beside `T22P_ZERO_RESOLUTION` for scale, never asserted — as 310-Q1.

### §9 item 387 — 2026-10-11 (chain T35, daily provisional under §9 rule (6); chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-11 18:00 daily review; the review reviews it.*

**Licence.** §9 rule (6) (PROJECT_PLAN.md :3812–3822): the 18:00 review writes T35 as the next numbered step of live §10 chain T, tagged *(daily, 2026-10-11 — provisional)*; the weekly ratifies or reverts. T35 is T26 (§10 chain T, PROJECT_PLAN.md :4651–4655, quoted, not changed: "the in-model zero at C_in(r) at 64.000 MHz within `T21P_ZERO_RESOLUTION`, circuit = field at `REDUCTION_BAND`; printed Q(r)/Q_0 against 1/Q = 1/Q_0 + 1/Q_C, Q_C = 1/(ωC_in r). Target: identities.") repeated at **127.74 MHz** with the 128 MHz knob that items 283 (T22′), 288 (T23 `XT128`) and 365 / 365b (T34) used: `FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74`, T22′'s seven ×1 Cu nodes. The band is the one T22′ registered for this grid: "within `T21P_ZERO_RESOLUTION`" becomes **`T22P_ZERO_RESOLUTION`** (1.616e-5). Item 283 registered it before its W1 (`test_ans7_exact_strip.py:1019–1025`), and the knob selects it through `_zero_resolution()` (`:1028–1030`), never a reading. The ESR ruling (§9, PROJECT_PLAN.md :3791–3810) applies: r = 0.2 Ω, frequency-independent, on every capacitor, on both sides at once. Lossless records stay, labelled lossless, and r never re-reads them. The sheets ruling (:3789) also applies: exact rectangles, κ `None`. **No band, tolerance, Done-when or target is set or moved; every band is imported by name.**

**Fixture check (done 2026-10-11 18:00).** Item 311 ran on the ×1 exact-rectangle F-small hole mesh: h = 0.015, 80 798 cells, degree 2, Cu, `-n 8` (`20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:2137`). Item 283's 127.74 MHz `cin` window ran on the same mesh, also 80 798 cells (`20261006T171745Z_ANS-7-item283-W2-cin128.log:2066`). The stored 127.74 inputs record 80 798 cells (`t127.74/tune_inmodel_cin_p2.json`). The fixture is item 311's.

`ANS-7` step 11b: tune the exact-sheet F-small birdcage on C_in with r = 0.2 Ω on both sides, at 127.74 MHz.
- **Circuit prediction:** C_in(r) is the root of Im[1/Z_in(C; r)] at 127.74 MHz on the stored ×1 `t127.74/tune_p2.json` `z64["cu"]`, found by item 311's `_c_in_esr`. It is pure and takes seconds.
- **The check:** the in-model solve with item 309's series-RC sheets at (r, C_in(r)), at 127.64 / 127.74 / 127.84 MHz.
- **Executor:** implementer. **Edits:** tests only, `tests/validation/test_ans7_exact_strip.py`. No `src/` edit and no `ans:7` edit (T36, item 391, does that).
- **Tier:** heavy, `-n 8` (item 311's / 283's width), complex build. **Class:** chain step (solve).
- **Time:** ≈ 12 measured (chain-solve class median, 10:30 restock table) / ≈ 42 predicted min: 15 fixed, ≈ 12 code, W1 ≈ 4.4, W2 ≈ 4.5, ≈ 6 spare. W1 is priced from item 311 W2 (260 s). W2 has the same drive count as item 283 W2 (261 s; drives 45.2 / 46.6 / 47.9 s, `…171745Z…:2066`), plus pure secants.
- **Depends on nothing open.** Item 311 has landed (`8e78c6b6`). Inputs, all on disk (checked 2026-10-11 18:00):
  - `logs/ans7_exact/t127.74/tune_p2.json`: degree 2, `f_hz` 127 740 000.0, `c_in` 1.5274843607471887e-11 F, `c_mode`, `z64`.
  - `t127.74/tune_inmodel_cin_p2.json`: lossless `rows_a`, Q 618.338591513781 and R_in 23 793.40 Ω; ranks 8; cells 80 798.
  - The seven ×1 `t127.74/node_*_p2_cu.json`.
  - Item 311's `logs/ans7_exact/tune_p2_esr0.2.json` and `tune_inmodel_cin_esr0.2_p2.json`, which W1 rewrites.
  - If any is absent, stop and mark BLOCKED. Never re-create one.

**Code** (`tests/validation/test_ans7_exact_strip.py`). The knob-unset path is untouched by construction: `_esr_ohm()` returns `None` before the edited condition, so item 262's / 283's lossless paths do not run any changed line, and there is no item 262 / 283 re-gate. Item 311's 64 MHz ESR path must stay behaviourally identical, and W1 shows it (rules (c) / (i)).
1. **`_esr_ohm()`** (`:2174–2185`): `_tune_mhz() is not None` becomes `_tune_mhz() not in (None, 127.74)`. At 127.74, also require `_cols() == ("cu",)` (item 283's trap: the stored nodes are `_cu`-suffixed). `RESOLUTION` and `TUNE_POINT` must still be unset; the h = 0.01 rung is not this item. Anything else raises `ValueError`, and the message names both legal sets.
2. **`test_t26_inmodel_zero_esr`** (`:2417–2424`): `T21P_ZERO_RESOLUTION` becomes `_zero_resolution()`. At 64 MHz that returns the same object, so the assertion is bitwise the same. The docstring loses its "(63.9, 64.1) MHz" literal. **This is the band-transplant trap:** without this edit, the 127.74 window would assert the 64 MHz grid's band (1.287e-5), which was never registered for this grid.
3. **`tune_esr` prints** (`:2280–2386`). Add `t35 = _tune_mhz() == 127.74`. On the 64 path every print string stays **byte-identical**. At 127.74:
   - the tag is `T35` in place of `T26`;
   - `/64 MHz` and `64.0 MHz` become `/127.74 MHz` and `127.74 MHz`;
   - the (AB-r) point format is `.2f`, never `.1f` (item 365's trap: 127.64 / 127.74 / 127.84 would all print as "127.6" / "127.7" / "127.8");
   - the band is printed by name and value from `_zero_resolution()`;
   - "item 262" becomes "item 283", "lossless -2.536e-6" becomes "lossless -2.265e-6", and "lossless 493.78 x" becomes "lossless 493.82 x";
   - every PREDICTED clause is the 127.74 text under **Predictions** below. No 64 MHz prediction is inherited: item 365b's C1 printed two inherited "NOT met" clauses, and this item must not.
4. **One record literal**, printed at 127.74 only and never asserted: `T26_C_IN_R_REL_64 = 4.602610e-5`, with provenance `20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:359`. The C_in(r) print shows it beside the 127.74 reading.
5. **Nothing else changes.**
   - The target, tune points, nodes, out dir and JSON names already flow from `_active()`, `_tune_json()` and `_stored_active()`. At 127.74, rank 0 writes **`t127.74/tune_p2_esr0.2.json`** and **`t127.74/tune_inmodel_cin_esr0.2_p2.json`**, and never the 64 MHz ESR files or any lossless file.
   - The `xn = nodes − 64.5` offset stays, as at both fields today (item 288's code step 4; item 283 interpolated the same way).
   - The tune-point index `j` resolves to 127.74.
   - (W), item 262's degree-1 64 MHz ×1 drive, runs unchanged in the 127.74 window, as in item 283 W2.

**Anchors** (asserted; every band imported by name), at 127.74 MHz:
- **(AB-r′) circuit = field at C_in(r):** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` (1e-3, `tests/validation/test_port_lumped_rlc_termination.py:88`) at 127.64 / 127.74 / 127.84 MHz (`test_t26_ab_circuit_equals_field_esr`).
- **(T2-r′) the in-model zero:** exactly one − → + crossing of Im `Y_in`,A in (127.64, 127.84) MHz, and \|f₀,A/127.74 MHz − 1\| ≤ `T22P_ZERO_RESOLUTION` via `_zero_resolution()` (`test_t26_inmodel_zero_esr`).
- **(W):** unchanged; ≤ `NODE_S_RTOL` against item 195's stored point.
- Cells equal the stored node cells (80 798), and `_assert_tie` passes.

**Negative controls (rule (e)).**
- **(W): asserted.** Backed by item 283 W2 (1.470e-13, `…171745Z…:2053`) and item 365 W0a: the same degree-1 64 MHz ×1 drive, the same comparison, in a 127.74 window.
- **(N-r0′): predicted, printed.** `_c_in_esr(z64, 127.74 MHz, C_mode, 0)` equals `t127.74/tune_p2.json`'s `c_in` bit for bit (64 MHz: bit-equal, `…094028Z…:360`). It is not asserted: that comparison has never run on the 127.74 `z64`.
- **(N2-r′): predicted, printed.** The circuit at 0.5·C_in(r) against the in-model at C_in(r) is ≥ 10 × `REDUCTION_BAND`. Lossless 127.74 MHz read 493.82×; ESR at 64 MHz read 494.39×.
- **(N-L′): predicted, printed — predicted NOT to separate at 10×.** The lossless circuit at C_in(0) against the in-model at (r, C_in(r)) on the bare S11 is predicted in [0.5, 5] × `REDUCTION_BAND`.
  - At 64 MHz it read 4.97× against a ≥ 10× prediction (item 311, a negative reported).
  - r's weight falls as Q_C(64)/Q_C(127.74) ≈ 0.48.
  - The separating reading is the Z_in form the code already prints (`:2372–2377`): \|dZ\|/\|Z\| predicted ≈ Q_0/Q_closed − 1 = 1.52, in [1.3, 1.8]. The same formula gives 4.90 at 64 MHz, where the reading was 4.878.

**Predictions** (written first, printed, never asserted).
- **C_in(r)/C_in(0) − 1 positive, in [+2e-6, +5e-5], central ≈ +1.06e-5**, printed beside 64 MHz's +4.602610e-5.
  - Reasoning, so that a miss is legible: a series-RC branch keeps its susceptance with C·(1 + 1/Q_C²), and 1/Q_C² is 6.01e-6 at 127.74 against 2.62e-5 at 64.
  - The 64 MHz reading was 1.76 × 1/Q_C², and carrying that factor gives 1.06e-5.
- **Q_C and the closed form.** Q_C = 1/(ωC_in r) = 407.84 at 127.74 MHz, C_in 15.2748 pF and r = 0.2 Ω. Q_0 = 618.34 is the lossless Route A Q at C_in (item 283; `t127.74/tune_inmodel_cin_p2.json` `rows_a`; `…171745Z…:2057`). The closed form 1/Q = 1/Q_0 + 1/Q_C gives Q ≈ 245.75, ratio 0.3974.
- **Q(r) in [240, 255], Q(r)/Q_0 in [0.38, 0.42], and Q(r)/closed form − 1 in [0, +3 %].** At 64 MHz the last read +1.10 %.
- R_in(r)/R_in(0) in [0.37, 0.43] (64 MHz: 0.169483, against Q ratio 0.169446).
- f₀,B(C_in(r))/127.74 MHz − 1 within ±1e-9, by construction.
- r_AB = f₀,A/f₀,B − 1 in [−5e-6, −1e-6]. Lossless 127.74 read −2.264620e-6; ESR moved it −1.9 % at 64 MHz.
- (AB-r′) in [1.0e-5, 3e-5] at each point (lossless 127.74: 1.590 / 1.601 / 1.616e-5).
- (T2-r′) \|f₀,A/127.74 − 1\| ≈ 2.3e-6, about 0.14 × `T22P_ZERO_RESOLUTION`: **predicted green with a ×7 margin.**
- W2 260–300 s; cells 80 798.

**Windows** (`-n 8`; each bare, whole and foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first, which includes `test_complex_mode.py`, as in items 283 / 311's heavy windows; `FEM_EM_REQUIRE_COMPLEX=1`; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE`, `FEM_EM_ANS7_EXACT_MHZ`, `RESOLUTION` and `TUNE_POINT` unset). No edit after W1 starts (rule (i)), so W1 and W2 run the same bytes. Each window is under 10 min, so §5.1's capture idiom is not required. Never two heavy windows at once.
0. **Orphan check:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
1. **W1, the 64 MHz ESR re-gate** (the edit touches item 311's path): item 311's W2 command verbatim, label `ANS-7-item387-W1-esr64`, 260 s.
   - It must reproduce, to the printed digits: C_in(r) 6.359301781809668e-11 F; C_in(r)/C_in(0) − 1 +4.602610e-05; (AB-r) 1.233e-5 / 1.258e-5 / 1.283e-5; (T2-r) −2.583111e-6 (0.201 bands); r_AB −2.583113e-6; Q(r) 164.5434331710199; (N2-r) 494.39×.
   - The (W) residual is not digit-stable across windows (4.034e-14 / 1.127e-14 / 3.927e-14 in items 262 / 283 / 311); it must only stay ≤ `NODE_S_RTOL`.
   - It rewrites item 311's two 64 MHz ESR JSONs with the same digits at `-n 8`. That is the accepted exception, as items 283 / 310 / 365's W0 were. Item 378's (R-T-r) reads those files and checks `ranks == 8`, so the rewrite is harmless.
   - **Red ⇒ revert the edit, BLOCKED.**

   `scripts/testing/run_and_log.sh ANS-7-item387-W1-esr64 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_ESR_OHM=0.2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
2. **W2, 127.74 MHz (closing):** label `ANS-7-item387-W2-esr128`, ≈ 260–300 s; on a ×1.4-slow box ≈ 400 s, inside 590.

   `scripts/testing/run_and_log.sh ANS-7-item387-W2-esr128 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_ESR_OHM=0.2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
3. **Commit** the module, both logs, their test-results rows (appended by hand) and the §9 marking together. `git diff --stat` on record files must show insertions only.

**Traps already paid for.**
- Items 217 / 228 / 262 / 283 / 310 / 311 / 365's lists, whole:
  - κ `None`;
  - no `ans:7` call;
  - `z_from_s(S)`, never `z_matrix`;
  - a real-scalar secant;
  - rank 0 writes, then a barrier;
  - a killed window is not a result: run `pgrep -c python3`, report it, and never lengthen the timeout;
  - never `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ` (`_nodes_env` validates against the 64 MHz nodes).
- **`COLUMNS=cu` is required at 127.74.** Without it the loader looks for unsuffixed node files that do not exist. Code step 1 makes the omission raise.
- **A missed 64 MHz literal on the 127.74 path** shows up as a bracket, `f_hz` or band failure on a correct measurement (item 228 code step 4). The known one is `T21P_ZERO_RESOLUTION` at `:2424` (code step 2).
- **`.1f` collides at 127.74** (item 365). Prints use `.2f` on this path.
- **Never enter 0.2 Ω as a parallel R** (ESR ruling (iii)): the sheet is the series RC. **The feed's branch carries r too** (ruling 311-Q5).
- **Never write** `t127.74/tune_p2.json`, any `tune_inmodel*_p2.json`, or the h = 0.01 files. W1's rewrite of item 311's ESR files is the one accepted exception.
- **r never re-reads a lossless record.** Q_0 is read and never rewritten.
- Exact rectangles only. Never a sub-sheet out of a box.

**Scope.** Covered: the ×1 exact-sheet mesh, degree 2, κ off, Cu, 127.74 MHz, r = 0.2 Ω uniform on all four capacitors on both sides, Route B on item 228's stored ×1 Cu grid. **Not claimed:**
- a measured ESR or its frequency dependence (r is the ruling's input at both frequencies);
- any absolute tuned C, f₀, R_in or Q at 3 T;
- an h-statement, or ESR on the h = 0.01 rung;
- a SPEC change;
- any AED re-run or Palace comparison (none exists at 3 T with ESR).

**Status it can move.**
- **No glyph.** The `ANS-7` Done-when is "✅ when T30 / F31 write Q's h-statement; 🟡 on the Q clause until then" (row, 2026-10-11 03:00 ruling), so the row stays 🟡.
- The row gains a "step 11b (T35)" record clause: C_in(r)/C_in(0) − 1 at 127.74 beside +4.60e-5 at 64 MHz, (AB-r′), (T2-r′) in `T22P_ZERO_RESOLUTION` bands, and Q(r)/Q_0 beside the closed form.
- §10 T35 is recorded DONE (provisional; the weekly ratifies).
- T36 (item 391, `ans:7`'s 127.74 MHz ESR rows) becomes runnable once item 378 has also landed.
- **Rule (j):** declined, as for items 283 / 311.

**What a negative result means.**
- **W1 red:** the edit changed item 311's path. Revert it and mark BLOCKED.
- **(W) red:** BLOCKED (wiring).
- **(AB-r′) red with (W) and W1 green:** the series-RC sheet and the circuit's (r + C) termination disagree at 127.74, but agree at 64 MHz. Print the three residuals, open a known-issues entry, park on `attempt/ANS-7-step11b-<UTC>` and mark BLOCKED. Never apply κ, never adjust r, never re-band.
- **No crossing in (127.64, 127.84):** BLOCKED. Report the Im `Y_in` triple and C_in(r).
- **(T2-r′) band half red:** report r_AB and C_in(r)/C_in(0) − 1, park BLOCKED, and the weekly decides. Never re-band.
- **A printed prediction missed:** record it as a reading, with no disposition (the missed-prediction rule of item 311's (N-L)).
- No `PREDICTED RED:` line: every assertion is predicted green, with a ×7 margin on the band.

**Spawn prompt (if the slot delegates; the slot may also execute directly):** "Use the implementer agent to execute `ANS-7` step 11b (§9 item 387) per `docs/planning/chunks/ANS-7.md` § '§9 item 387'. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run every harness window foreground with Bash `timeout: 660000`; never background a window or end your turn while one is running. Run the orphan check first. Edit only `tests/validation/test_ans7_exact_strip.py`, and make no edit after W1 starts. Every window runs with `-s`, `FEM_EM_REQUIRE_COMPLEX=1` and `FEM_EM_ANS_SKIP_PRIVATE=1`; `FEM_EM_ANS7_EXACT_ROUTE`, `FEM_EM_ANS7_EXACT_MHZ`, `FEM_EM_ANS7_EXACT_RESOLUTION` and `FEM_EM_ANS7_EXACT_TUNE_POINT` stay unset. Never write `logs/ans7_exact/t127.74/tune_p2.json` or any `tune_inmodel*_p2.json`. Never open `docs/private/` or any `aed_results/`. Report every printed reading, met or NOT met."

**Outcome 2026-10-12 (00:30Z implementer slot): DONE (record; row stays 🟡; T35 DONE, provisional).**
- W1 (64 MHz ESR re-gate, `20261012T003929Z_ANS-7-item387-W1-esr64.log`, 268 s, 15 passed): C_in(r) 6.359301781809668e-11 F, +4.602610e-05, (AB-r) 1.233e-5 / 1.258e-5 / 1.283e-5, (T2-r) −2.583111e-6 (0.201 bands), r_AB −2.583113e-6, (N2-r) 494.39×, (N-L) 4.97× — all to the printed digits. Q(r) printed 164.5434331273751 against item 311's 164.5434331710199 (2.7e-10 relative): S11_A moved at 1e-12 relative between windows while S11_B is bitwise equal, so this is the solve's run-to-run round-off (like (W)), not the edit. (W) 7.091e-14.
- W2 (127.74 MHz, `20261012T004412Z_ANS-7-item387-W2-esr128.log`, 279 s, 15 passed; 80 798 cells, `-n 8`; drives 51.1 / 45.9 / 47.6 s):
  - C_in(r) 1.5275048999871143e-11 F; C_in(r)/C_in(0) − 1 = +1.344645e-05 (prediction [+2e-6, +5e-5] met; 0.2921 × the 64 MHz +4.602610e-5).
  - (N-r0′) bit-equal (met). (W) 1.062e-13 ≤ `NODE_S_RTOL`.
  - (AB-r′) 1.587e-5 / 1.600e-5 / 1.614e-5 ≤ `REDUCTION_BAND` (prediction [1.0e-5, 3e-5] met).
  - (T2-r′) Im `Y_in`,A −4.057389094e-05 / +1.196522503e-07 / +4.085113314e-05 S, one crossing at 127 739 705.97 Hz; f₀,A/127.74 MHz − 1 = −2.301805e-06 = 0.142 × `T22P_ZERO_RESOLUTION` (met; prediction ≈ 2.3e-6 met). f₀,B − 1 = −9.16e-13 (met); r_AB −2.301804e-6 (met).
  - Q(r) 254.469 (prediction [240, 255] met), Q(r)/Q_0 0.411537 (prediction [0.38, 0.42] met); closed form Q 245.7467, ratio 0.397431; Q(r)/closed form − 1 = **+3.5494 %, prediction [0, +3 %] NOT met** — a reading, no disposition. R_in(r)/R_in(0) 0.411488 (prediction met). Route B: R_in 9 772.95 Ω, Q_B 254.003.
  - (N-L′) 2.17 × `REDUCTION_BAND` (prediction [0.5, 5]× met; does not separate at 10×); Z_in form 1.419 (prediction [1.3, 1.8] met). (N2-r′) 494.07× (met).
- Wrote `t127.74/tune_p2_esr0.2.json` and `t127.74/tune_inmodel_cin_esr0.2_p2.json` only; no lossless file touched.

### §9 item 391 — 2026-10-11 (chain T36, daily provisional under §9 rule (6); example)

*Drafted by a general-purpose agent for the 2026-10-11 18:00 daily review; the review reviews it.*

**Licence.** §9 rule (6) (PROJECT_PLAN.md :3812–3822): T36 is a provisional next step of live chain T, tagged *(daily, 2026-10-11 — provisional)*, and the weekly ratifies or reverts. T36 is T27's sibling at 127.74 MHz. T27 (§10 chain T, PROJECT_PLAN.md :4656–4657, quoted, not changed) reads: "`ans:7` gains the ESR rows (T26 at 64 MHz; lossless rows stay, labelled). `example-runner`. Serial on T26." So `ans:7` gains the **127.74 MHz** ESR rows from item 387's (T35's) stored run. The ESR ruling (§9, :3791–3810) and the sheets / κ ruling (:3789) apply. **No band, tolerance, Done-when, SPEC or validation target is set or moved; every band is imported by name.**

**Housing.** Housing follows item 378: no EX chunk owns `ans:7`, so the case re-issues as an `ANS-7` step. T36 is **`ANS-7` step 13b**, beside step 13 (T27, item 378).

**Depends on items 387 and 378 — skip if either did not land, and leave 391 ⬜.**
- Item 378 lands `esr_ohm` on `z_cap` / `add_c1` / `route_b_z` / `_route_a_solve`, the 64 MHz stage `XT64R`, and `degree2_exact["tuned_larmor_esr"]["64.0"]`. Item 391 extends that stage by name and never re-implements it.
- Item 387 writes this stage's stored inputs: `logs/ans7_exact/t127.74/tune_p2_esr0.2.json` (`c_in_r`, `c_in_0`, `esr_ohm`, `f_hz`) and `t127.74/tune_inmodel_cin_esr0.2_p2.json`. Neither exists today (checked 2026-10-11 18:00). If either is absent when the slot opens, stop and take the next item. Never re-create one; that is the test module's job.

`ANS-7` step 13b: `ans:7` gains one additive stage, **`XT128R`**. It re-drives Route A at degree 2 on the ×1 exact-sheet mesh at (r = 0.2 Ω, C_in(r)) at 127.64 / 127.74 / 127.84 MHz, with item 378's series-RC terminations. It re-asserts T35's gates through the example and writes `degree2_exact["tuned_larmor_esr"]["127.74"]` and a 127.74 MHz column in `COMPARISON.md`'s ESR section. The lossless `tuned_larmor` rows and the 64 MHz ESR rows stay byte-identical.
- **Executor:** `example-runner`. Spawn prompt: "Use the example-runner agent to execute `ANS-7` step 13b (§9 item 391) per `docs/planning/chunks/ANS-7.md` § '§9 item 391'. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run every harness window foreground with Bash `timeout: 660000`. Never background a window or end your turn while one is running. Set `FEM_EM_ANS_SKIP_PRIVATE=1` on every window that runs the `ans:7` script. Never run `ans:1`–`ans:3`, never run B1 / B2 / A / X1–X4 / XA / XT64 / XT128, run XT64R only as the brief's W-reg window when its condition holds, and never use `./scripts/run_examples.sh -e ans`. Never open `docs/private/`, any `aed_results/` or `COMPARISON_private.md`."
- **Edits:** the example script, its same-stem guide, `metrics.json` and `COMPARISON.md` only. No `src/`, test-module, SPEC or `aed/` edit.
- **Standing rule (a):** this item does **not** edit the gate-imported helpers (`z_cap`, `add_c1`, `route_b_z`, `_route_a_solve`), because item 378 already keyworded them. If the executor finds it must, item 378's **W-gate** (item 310's W5 verbatim) runs in the same slot and the diff is disclosed.
- **Tier:** heavy, `-n 8` (T35's record width), complex build. **Class:** example.
- **Cost:** ≈ 8 measured (example class median) / ≈ 40 predicted min: 15 fixed, ≈ 12 code and guide, `XT128R` ≈ 3, ≈ 3 census and pin, ≈ 7 spare.
  - Price: item 288's `XT128`, 165 s (`20261007T003522Z_ANS-7-item288-XT128.log`), which has the same drive count. Predicted 150–280 s, ≈ 390 s on a ×1.4-slow box, inside 590.
  - If W-reg is needed, add ≈ 3 min (item 378's `XT64R`, same class).

**Code** (`examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py`; additive).
1. **Stage `XT128R`** is added to `main()`'s accepted set (`:1503`, as item 378 left it) and the error message. It starts with `X`, so `_x_guard()` applies: every `FEM_EM_ANS7_EXACT_*` knob, `FEM_EM_ANS7_EXACT_ESR_OHM` included, must be **unset**. The field set is chosen by stage name alone (item 288's knob-free trap).
2. **A 127.74 branch by stage name** inside item 378's ESR stage, mirroring `_xt_stage`'s `is64` split (`:1301–1316`):
   - key `"127.74"`, target `TUNE128_TARGET_HZ`, tune `TUNE128_MHZ`;
   - nodes `TUNE128_NODES_MHZ`, read as the seven ×1 `_node128_json(m, ("cu",))` files with `_xt_stage`'s consistency checks (`:1331–1340`);
   - grid `TUNE128_GRID_MHZ`, out `TUNE128_OUT_DIR`;
   - band `T22P_ZERO_RESOLUTION`, named `"T22P_ZERO_RESOLUTION"`;
   - the tag strings `T35` / `(AB-r′)` / `(T2-r′)`;
   - `.2f` point formats.
   - Everything is imported lazily from `tests.validation.test_ans7_exact_strip` (item 230's cycle trap), never restated. **Preferred shape:** every line `XT64R` executes stays byte-unchanged, with the new branch reached only by the new stage name. Then no `XT64R` re-run is needed, and (K) proves its entry unchanged.
3. **Inputs.**
   - `TUNE128_OUT_DIR / "tune_p2_esr0.2.json"`: `c_in_r`, `c_in_0`, `esr_ohm` (assert it equals 0.2, the ruling's value, never a band), and `f_hz` within 1 Hz of 127.74 MHz. **C_in(r) is never recomputed.**
   - `TUNE128_OUT_DIR / "tune_inmodel_cin_esr0.2_p2.json"` is the (R-T-r′) comparand.
   - Q_0 is the lossless `degree2_exact["tuned_larmor"]["127.74"]` `rows_a.q_loaded` (item 288's `XT128` row; the coil's loaded Q with lossless capacitors).
4. **Drives and rows:** item 378's loop at the 127.74 tune points. f₀,B comes from the Route B Im `Y_in` over `TUNE128_GRID_MHZ` through the keyworded `route_b_z`, as item 378 forms it at 64 MHz. **Printed, never asserted:** Q(r)/Q_0 against 1/Q = 1/Q_0 + 1/Q_C with Q_C = 1/(ωC_in(r)·r), and C_in(r)/C_in(0) − 1 beside the 64 MHz `tuned_larmor_esr["64.0"]` value.
5. **Write** (rank 0, then a barrier). Snapshot every top-level key and every `degree2_exact` key, **including `tuned_larmor` and `tuned_larmor_esr["64.0"]`**. Set `degree2_exact["tuned_larmor_esr"]["127.74"]` with item 378's field list, with `lossless_comparand: "tuned_larmor['127.74']"` and `source` naming item 387's W2 log. Assert (K), and revert `metrics.json` on a miss.
6. **`COMPARISON.md`.** Item 378's section "Tuned with capacitor ESR (r = 0.2 Ω on every capacitor; §10 T26 / T27)" gains a **127.74 MHz column**, the same rows, and "/ T35 / T36" in its heading. That is the only change to existing text, and it is disclosed in the commit. **There is no AED column.** No AED ESR case exists.
7. **Docs.**
   - The docstring's `XT64R` paragraph names `XT128R`.
   - The guide's "Capacitor ESR" section gains the `XT128R` command and anchors.
   - Any 64-MHz-only clause item 378 wrote in "Not claimed" is narrowed to the two frequencies (disclosed).
   - Stored inputs are named as "the gitignored `logs/ans7_exact/` stage files" (item 273's docrefs trap). Write ≈, never a tilde. No setup figure or XDMF.

**Anchors** (asserted in-stage; every band imported):
- **(AB-r′)** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` at 127.64 / 127.74 / 127.84 MHz.
- **(T2-r′)** Exactly one − → + crossing of Im `Y_in`,A in (127.64, 127.84), and \|f₀,A/127.74 MHz − 1\| ≤ `T22P_ZERO_RESOLUTION`. These are T35's anchors as item 387 asserts them.
- **(R-T-r′) reproduction of item 387's registered run** against `t127.74/tune_inmodel_cin_esr0.2_p2.json`, applied as item 378's (R-T-r):
  - \|f₀,A,ex/f₀,A,stored − 1\| ≤ `NODE_S_RTOL`;
  - max \|res_ex − res_stored\| ≤ `NODE_S_RTOL`;
  - max \|S11_A,ex − S11_A,stored\|/\|S11_A,stored\| ≤ `NODE_S_RTOL`;
  - the stored `ranks == 8` and `degree == 2`.
  - At any other width it prints and does not assert (the `EX-56` rule).
- **(K)** Every pre-existing `metrics.json` key, `tuned_larmor` and `tuned_larmor_esr["64.0"]` included, is byte-identical.
- **(R0′)** The r = 0 branch is bitwise at 127.74: `z_cap`, `add_c1` and `route_b_z` with `esr_ohm=0.0` equal the unkeyworded call, bit for bit, on the stored ×1 `t127.74/tune_p2.json` `z64["cu"]` at 127.74 MHz. This is item 378's (R0) at the second frequency; it is pure and takes milliseconds.

**Negative controls (rule (e)).**
- **(N-C-S′) the comparator separates ESR from lossless: asserted only if backed.** max \|S11_A,ex(r) − S11_A,stored,lossless\|/\|S11_A,stored,lossless\| > 100 × `NODE_S_RTOL`, against `t127.74/tune_inmodel_cin_p2.json` `s11_a`.
  - It is **asserted only if** the two stored JSONs (item 387's `s11_a` and the lossless `s11_a`) differ by > 1e-7 relative at every point.
  - Before the window the executor computes that with a pure read, cites it in the commit, and labels the control. Otherwise it is printed.
  - Predicted O(1e-3) at every point.
- **(N-C-r′) the f₀ form of item 378's (N-C-r): predicted, printed.** \|f₀,A,ex(r)/f₀,A,stored(lossless) − 1\| > 100 × `NODE_S_RTOL`. It is asserted only under item 378's condition: item 387's outcome records f₀,A(r) and f₀,A(0) differing by > 1e-7. That difference is predicted ≈ 1e-8 to 5e-8 (64 MHz: 4.7e-8), so it is expected to be printed.
- **(N-L′): predicted, printed.** It is asserted only if item 387's reading **met** the ≥ 10× bar, and item 387 predicts it will not.

**Predictions (printed).** (R-T-r′) ≤ 1e-12 (item 288 reproduced item 262 to 9.97e-13 at `-n 8`). Every other reading is item 387's, to the printed digits.

**Windows** (`-n 8`; each bare, whole and foreground; Bash `timeout: 660000`; no edit after `XT128R` starts, rule (i)). Each window is under 10 min, so §5.1 capture is not required.
0. **Orphan check:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
1. **Censuses before**, labels `ANS-7-item391-pre-figures` / `-pre-docrefs`: item 378's / 288's commands (`check_example_setup_figures.py`, `check_example_doc_references.py`, `timeout -k 30 60`).
1b. **W-reg, only if the edit changed any line `XT64R` executes.** Run item 378's `XT64R` command verbatim, label `ANS-7-item391-W-reg-XT64R`.
   - The rewritten `tuned_larmor_esr["64.0"]` must equal the committed entry in every field except `t_mesh_s` / `t_drive_s` (pure JSON compare, printed).
   - Any other difference ⇒ revert, BLOCKED.
   - If the preferred shape held, skip W-reg and say so in the commit.
2. **XT128R:**

   `scripts/testing/run_and_log.sh ANS-7-item391-XT128R "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_STAGE=XT128R timeout -k 30 590 mpiexec -n 8 python3 examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py'"`
3. **Pin and unit tests** on the committed `metrics.json`, `-n 1`, label `ANS-7-item391-pin`: item 378's command (`tests/unit/test_doc_reference_exit_codes.py tests/unit/test_ans7_private_comparison.py`, `-s`). The tracked-output guard must stay green with the new key.
4. **Censuses after**, labels `-post-figures` / `-post-docrefs`. Predicted equal to the before readings (`broken=0`; docrefs `dead=0`, `guide=0`, exit ≠ 1).
5. **Commit** the script, guide, `metrics.json`, `COMPARISON.md`, logs and test-results rows (appended by hand) together. `git diff --stat` must show no `COMPARISON_private.md`, `aed_results/`, `SPEC.md` or test-module change.

**Traps already paid for.**
- Items 230 / 288 / 378's lists, whole:
  - lazy imports;
  - never run B2 (`_write_outputs` rewrites wholesale);
  - never rename or move the degree-1 keys;
  - `z_from_s(S)`;
  - κ `None`;
  - AED values never enter a tracked file;
  - a killed window is not a result;
  - bare harness calls;
  - knob-free selection: setting `FEM_EM_ANS7_EXACT_ESR_OHM` or `FEM_EM_ANS7_EXACT_TUNE_MHZ` trips `_x_guard`;
  - `FEM_EM_ANS_SKIP_PRIVATE=1` on every window that runs the script;
  - the r = 0 path is a branch, not arithmetic (item 308);
  - never enter 0.2 Ω as a parallel R;
  - the feed's branch carries r too.
- **A 64 MHz literal on the 127.74 branch** (`64.5` is the interpolation origin at both fields and stays, as in `XT128`; anything else is a defect). The known candidates are `T21P_ZERO_RESOLUTION`, the `GRID_MHZ` f₀,B grid, `.1f` formats and the Q_0 key `"64.0"`.
- **(N-L′) / (N-C-r′) conditions are "measured and met", not "measured".** Asserting a control that item 387 recorded as missed turns a recorded negative into a red window.

**Scope.** Covered: the example re-asserting T35 at (r, C_in(r)) on the ×1 exact-sheet mesh, degree 2, κ off, Cu, 127.74 MHz, and recording the ESR rows beside the labelled lossless rows. **Not claimed:**
- a measured ESR or its frequency dependence;
- any absolute tuned C, f₀, R_in or Q;
- any h-statement;
- any AED or Palace comparison;
- the private writer's rows.

**Status it can move.**
- §10 T36 is recorded DONE (provisional; the weekly ratifies).
- The `ANS-7` row gains a "step 13b (T36)" clause: `ans:7`'s tracked outputs carry the ESR rows at both fields, with (AB-r′), (T2-r′) and (R-T-r′) re-asserted through the example.
- **No glyph.** **Rule (j):** declined, as for items 288 / 378.

**What a negative result means.**
- **(R-T-r′) red at `-n 8` with the bands green:** the example's 127.74 ESR path is not item 387's. Print the differences, revert `metrics.json` / `COMPARISON.md` (those two paths only), park on `attempt/ANS-7-step13b-<UTC>` and mark BLOCKED. Never re-register.
- **(R0′) red:** revert, BLOCKED.
- **(AB-r′) / (T2-r′) red:** this contradicts item 387's green. Report it, open a known-issues entry if it reproduces, and mark BLOCKED. Never re-band or adjust r.
- **(K) red, or W-reg red:** revert both tracked outputs, BLOCKED.
- No `PREDICTED RED:` line.

**Outcome 2026-10-12 (00:00 CDT slot, example-runner): BLOCKED (W-reg).** The 127.74 MHz column was rendered by editing `_write_comparison`, which `XT64R` executes, so the preferred shape did not hold and W-reg ran (`20261012T051318Z_ANS-7-item391-W-reg-XT64R.log`, `-n 8`, 177 s; mesh 27.5 s, drives 50.2 / 51.8 / 45.4 s). Every in-stage assertion was green: (R0) bitwise, (AB-r) 1.233e-5 / 1.258e-5 / 1.283e-5, (T2-r) −2.583111e-6 (0.201 bands), (R-T-r) 7.96e-14 / 2.57e-13 / 3.14e-13, (K) yes. The pure JSON compare of the rewritten `tuned_larmor_esr["64.0"]` against the committed entry differs in 12 fields beyond `t_mesh_s` / `t_drive_s`: `s11_a` 3.07e-12 rel (`s11_b` bitwise), f₀,A 7.5e-14 rel, f₀,A-derived fields (`rows_a`, `f0a_rel_target`, `bands_fraction`, `r_ab`, `q_ratio`) ≤ 1.75e-10, `rt` re-measured. By the brief's letter ("any other difference ⇒ revert, BLOCKED") the item stops; `metrics.json` / `COMPARISON.md` reverted, `XT128R` not run. Pure reads before the window: (N-C-S′) item 387's `s11_a` vs lossless differ 2.169e-3 / 2.182e-3 / 2.197e-3 > 1e-7 at every point ⇒ assertable; (N-C-r′) separation 3.718e-8 < 1e-7 ⇒ printed; (N-L′) 387 read 2.17× ⇒ printed. Code (`_xt_esr128_stage`, two-column ESR writer, guide section; f₀,B by item 378's three-point construction; row label "f₀,A/64 MHz − 1" → "f₀,A/target − 1") parked on `attempt/ANS-7-step13b-20261012T0517Z`. **For the review:** W-reg's equality criterion cannot hold across two solves (run-to-run round-off, as item 387 noted); either re-read it as equality at `NODE_S_RTOL` for solve-derived fields, or re-pose with a separate `XT128R` comparison writer so no `XT64R` line changes and W-reg is skipped. No band moved.

### §9 item 395 — 2026-10-12 (chain T37, daily provisional under §9 rule (6); chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-12 03:00 daily review; the review reviews it.*

**Licence.** §9 rule (6) (PROJECT_PLAN.md :3813–3823): the 03:00 review writes T37 as the next numbered step of live §10 chain T, tagged *(daily, 2026-10-12 — provisional)*. This is chain T's third and last provisional step between weeklies, after T35 and T36, and the weekly ratifies or reverts it. T37 is T26 (§10 chain T, PROJECT_PLAN.md :4649–4653, quoted, not changed: "the in-model zero at C_in(r) at 64.000 MHz within `T21P_ZERO_RESOLUTION`, circuit = field at `REDUCTION_BAND`; printed Q(r)/Q_0 against 1/Q = 1/Q_0 + 1/Q_C, Q_C = 1/(ωC_in r). Target: identities.") run **on T28's h = 0.01 exact-rectangle rung** (item 310, `ANS-7` step 12): item 214's five Cu nodes, T28's C_in(h) and z64(h), and item 311's construction. It serves chain T's open question, Q's h-statement (T30 / F31; the weekly's Phase 6 paragraph, PROJECT_PLAN.md :4905–4907). With ESR, the tuned coil's Q is dominated by capacitor loss, so this step measures whether Q(r) moves with h less than Q_0 does.

The rulings that apply:
- the ESR ruling (§9 :3792–3811): r = 0.2 Ω on every capacitor, both sides at once; lossless records stay, labelled lossless;
- the sheets ruling (:3790): exact rectangles, κ `None`.

**The zero band on this rung.** T26's band is `T21P_ZERO_RESOLUTION`, registered from the ×1 grid's node residuals. Item 310 found that **no band exists for the h = 0.01 rung**: its node files carry no identity residuals (`field_solve` false). Ruling 310-Q1 therefore asserted only the crossing inside the bracket. It printed the offset beside `T21P_ZERO_RESOLUTION` for scale only (`20261008T112649Z_ANS-7-item310-W5-assemble.log:138`), and 365-Q1 repeated this at 127.74 MHz. No weekly has registered a rung band since. Rule (6) never sets a band, so T37 takes **310-Q1's form**: the zero half of T26's target is the crossing, and the offset is printed. **No band, tolerance, Done-when or target is set or moved; every band is imported by name.**

**Feasibility (checked 2026-10-12 against the code and the disk).** The two paths combine with a tests-only edit and no new node solve:
- **The h-path inputs are all on disk** (`logs/ans7_exact/`):
  - item 214's five `node_{62.0,63.2,64.5,65.8,67.0}_p2_h0.01_cu.json`;
  - item 310's `tune_p2_h0.01.json`: `c_in` 6.364316664701229e-11 F, `c_mode["cu"]`, `z64["cu"]`, `tie`, cells 140 270, `resolution` 0.01, `c_in_h_rel` 8.346525677918049e-4;
  - item 310's assembled lossless `tune_inmodel_cin_p2_h0.01.json`: `rows_a` Q 964.432729560342, R_in 18 592.58100053142 Ω, `c_in` and `cells`.
- **ESR's C_in(r, h) is pure.** It is `_c_in_esr` (`test_ans7_exact_strip.py:2215`) applied to that `z64`, so no 4×4 is re-solved.
- **What costs compute** is three degree-2 Route A drives on the rung, plus the (W) wiring drive.
- **Five literals and guards in the current code would block the h path or misread it** (code steps 1–6 below):
  - `_esr_ohm()` (`:2178–2194`) refuses `RESOLUTION`;
  - `_tune_point()` (`:1090–1109`) raises on the h `inmodel` path when `TUNE_POINT` is unset;
  - `tune_esr` builds the ×1 tie and mesh (`:2312`, `:2314`), reads the ×1 lossless Q_0 (`:2280`) and writes the ×1 ESR file names (`:2457–2466`);
  - `test_t26_inmodel_zero_esr` (`:2494`) asserts `_zero_resolution()`, which is `T21P_ZERO_RESOLUTION` when `TUNE_MHZ` is unset.
- **Price:** one ESR window on the rung is three ≈ 100 s drives (item 310's per-drive measurement), so it fits one window under 590 s at `-n 12`, and the whole item fits one ≤ 1 h slot. No a / b split is needed (see the price and the kill contingency below).

`ANS-7` step 11c: tune the exact-sheet F-small birdcage on C_in with r = 0.2 Ω on both sides, at 64.000 MHz, on the h = 0.01 rung (140 270 cells, degree 2, κ off, Cu).
- **Executor:** implementer. Spawn sentence: "The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped."
- **Edits:** tests only, `tests/validation/test_ans7_exact_strip.py`. No `src/`, `ans:7`, `test_ans7_sparse_route_b.py` or `test_ans7_t30_rin_q_sequence.py` edit.
- **Tier:** heavy. W1 runs at **`-n 12`**, the rung's measured width (items 214 / 310 / 365): summed `ru_maxrss` 33.8–36.7 GiB, and h = 0.01 is unpriced at `-n 8`. The re-gates run at their records' widths (`-n 8` / `-n 2`). Complex build. **Class:** chain step (solve).
- **Time:** ≈ 13 measured (chain-solve class median, 18:00 restock table, PROJECT_PLAN.md :4002) / **≈ 50 predicted min**: 15 fixed + ≈ 15 code + W0 ≈ 4.5 + W0b ≈ 1.5 + W1 ≈ 7 + ≈ 7 spare.
- **Depends on nothing open.** Items 310 and 311 have landed (outcomes above). If any input file named above is absent, stop and mark BLOCKED. Never re-create one.

**Code** (`tests/validation/test_ans7_exact_strip.py`). Every knob-unset path stays behaviourally identical. Item 311's ×1 64 MHz ESR path is shown unchanged by W0, item 387's 127.74 ESR path is untouched by construction (the edits branch on `h37` only), and item 310's h lossless path is shown unchanged by W0b.
1. **`_esr_ohm()`** (`:2186–2193`) gains a third legal set: `_tune_mhz() is None`, `_resolution() == H_RES` and `_cols() == ("cu",)`. `TUNE_POINT` must still be unset in every set. Any other `RESOLUTION` raises `ValueError`, and the message names all three legal sets.
2. **`_tune_point()`** (`:1093–1096`): the "unset on the h rung" `ValueError` fires only when `_esr_ohm() is None`. Under the ESR knob, unset returns `None`. There is no recursion: `_esr_ohm()` reads `TUNE_POINT_ENV` directly.
3. **`tune_h`** (`:2025–2031`): `if _esr_ohm() is not None: pytest.skip(...)` before `pt = _tune_point()`. Otherwise `pt * 1e6` fails on `None`.
4. **`tune_esr`** (`:2264–2469`): `h37 = _h_tune()`. `_active()`, `_tune_json()` and `_stored_active()` already return the T28 set under these knobs, giving `tune_p2_h0.01.json` and the five checked `_h0.01_cu` nodes. Then, on the h path only:
   - **Q_0 comparand:** `jl = act["out"] / f"tune_inmodel_cin_p{deg}{act.get('sfx', '')}.json"`. That is `tune_inmodel_cin_p2_h0.01.json` on the h path, and the name is unchanged on the other two paths (their `act` has no `sfx`).
   - **(F-Q0), asserted on the h path:** `lossless["c_in"] == tj["c_in"]` and `lossless["cells"] == recs[0]["cells"]`. A Q_0 taken from the wrong rung fails here, not silently.
   - **Tie without a mesh:** `tie = _tie_from_json(tj["tie"])` in place of `_fixture_tie()`. The ×1 call builds an extra ×1 mesh and asserts the ×1 record.
   - **(W) first, on the ×1 mesh.** Run `_exact_geo(th14, comm, GATED_WIDTH_FRACTION)`, then `_route_a_solve(…, TUNE_TARGET_HZ, "cu")`, then `del` the ×1 geo, as `tune_h` does (`:2066–2085`).
   - **Then build the rung:** `geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION, res=H_RES)`. The three drives are `_route_a_solve_esr` unchanged.
   - **Outputs:** rank 0 writes **`tune_p2_h0.01_esr0.2.json`** and **`tune_inmodel_cin_esr0.2_p2_h0.01.json`**, the ×1 payloads plus `"resolution": H_RES`. It never writes item 311's `tune_p2_esr0.2.json` / `tune_inmodel_cin_esr0.2_p2.json`, any lossless file, or any node file.
   - **The return dict** gains `"h": h37` and `"hres": H_RES if h37 else None`. The key `"res"` already holds the residuals.
   - **Prints:** tag `T37`, with the 64 MHz text. Every PREDICTED clause is from **Predictions** below, and no ×1 or 127.74 prediction is inherited (item 365b C1's two inherited "NOT met" lines). The ×1 and 127.74 print strings stay byte-identical.
5. **Record literals**, printed on the h path only and never asserted, each with its provenance comment:
   - `T26_C_IN_R_X1 = 6.359301781809668e-11` (`20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:359`);
   - `T26_Q_R_X1 = 164.5434331710199` (`…094028Z…:2130`);
   - `T26_Q_CLOSED_EXCESS_X1 = 1.0996e-2` (`…094028Z…:2133`);
   - `T21P_Q0_X1 = 971.0667834101893` (Q_0 at the ×1 C_in, `…094028Z…:2133`);
   - `T30_Q_H_REL_CF = -6.9811e-3` (Q at `C_F`, h = 0.01 against ×1, `20261011T094910Z_ANS-7-item364-T30.log:76`).
   - The C_in(h) move uses the existing `T28_C_IN_H_REL_64` (`:1068`).
   - These are literals because W0 rewrites item 311's ×1 ESR JSONs in this slot. A literal from the registered log cannot drift with the run-to-run round-off.
6. **Tests:**
   - `test_t26_ab_circuit_equals_field_esr` (`:2473–2481`): `_assert_tie(cells, tie, res=tune_esr["hres"])`.
   - `test_t26_inmodel_zero_esr` (`:2485–2494`): on the h path it asserts exactly one up-crossing, no down-crossing (`rows_a["crossings_down_hz"] == []`) and `lo < f₀,A < hi`. It **skips the band line**, and its docstring says why (310-Q1). The ×1 and 127.74 paths are unchanged.

**Anchors** (asserted; every band imported by name), at 64.000 MHz on the rung:
- **(AB-rh) circuit = field at C_in(r, h):** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` (1e-3, `tests/validation/test_port_lumped_rlc_termination.py:88`, imported through `test_port_circuit_layer_field` at `test_ans7_exact_strip.py:81`) at 63.9 / 64.0 / 64.1 MHz. Route B is the five-node h grid through `bary_interp` and `_route_b_z_esr`.
- **(T2-rh) the in-model zero (310-Q1's form):** Im `Y_in`,A changes sign − → + exactly once, and never + → −, in (63.9, 64.1) MHz, with f₀,A inside the bracket.
  - \|f₀,A/64 MHz − 1\| is **printed** beside `T21P_ZERO_RESOLUTION` (1.287e-5, `test_ans7_exact_strip.py:1018`) for scale only, never asserted.
- **(F) the rung:**
  - the drive-mesh cells equal the stored node cells (140 270) and `tune_p2_h0.01.json`'s;
  - `_assert_tie(…, res=H_RES)` passes;
  - its (N-wire) clause holds: the cells are off the ×1 record by > `CELL_COUNT_BAND` (0.01, `tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`).
- **(F-Q0):** the lossless comparand's `c_in` and `cells` are the rung's.
- **(W):** unchanged; ≤ `NODE_S_RTOL` (1e-9, `tests/validation/test_ans7_sparse_route_b.py:87`) against item 195's stored point.

**Negative controls (rule (e)).**
- **(W): asserted.** Same degree-1 ×1 64 MHz drive, same comparison. Backed by:
  - item 310 W3, on an h-rung window at `-n 12`: 1.720e-13 (`20261008T111904Z_ANS-7-item310-W3-p64.0.log:1430`);
  - item 311 W2: 3.927e-14 (`…094028Z…:2126`);
  - item 387 W1: 7.091e-14 (`20261012T003929Z_ANS-7-item387-W1-esr64.log:2126`).
  - The residual is not digit-stable across windows. It must only stay ≤ `NODE_S_RTOL`.
- **(N-wire): asserted** inside `_assert_tie(res=H_RES)`. The same comparison on the same rung was green in item 310 W2–W4 (`test_t28_rung_and_wiring`, 140 270 cells) and item 365's nodes.
- **(N2-rh): predicted, printed.** The circuit at 0.5·C_in(r, h) against the in-model at C_in(r, h) is ≥ 10 × `REDUCTION_BAND`. On this rung the lossless reading was 493.62× (`…112649Z…:143`), and ×1 ESR read 494.39× (`…094028Z…:2136`). It has never been measured with r on this rung.
- **(N-Lh): predicted, printed — predicted NOT to separate at 10×.** The lossless circuit at C_in(0, h) against the in-model at (r, C_in(r, h)), on the bare S11, is predicted in [3, 7] × `REDUCTION_BAND`.
  - At ×1 it read 4.97× against a ≥ 10× prediction (item 311, a recorded negative).
  - The Z_in form \|dZ\|/\|Z\| is predicted ≈ Q_0/Q_closed − 1 = 4.94, in [4.6, 5.2]. The ×1 reading was 4.878.
- **(N-r0h): predicted, printed.** `_c_in_esr(z64(h), 64 MHz, C_mode(h), 0)` equals `tune_p2_h0.01.json`'s `c_in` bit for bit. That was bit-equal at ×1 and 127.74, but never run on the rung's `z64`.
- **(N-h): predicted, printed.** \|C_in(r, h)/C_in(r, ×1) − 1\| ≥ 10 × the secant tolerance (\|g\| ≤ 1e-15 S → ≈ 1e-12 relative). This is the step's reading, so it is not gated.

**Predictions** (written first, printed, never asserted; drafter's arithmetic from the stored JSONs).
- **C_in(r, h)/C_in(0, h) − 1 in [+3e-5, +6e-5], central +4.610e-5.** The ×1 value +4.602610e-5 scales as C² (1/Q_C²) by (1 + 8.35e-4)².
- **C_in(r, h)/C_in(r, ×1) − 1 in [+8.0e-4, +8.7e-4], central +8.3473e-4**, printed beside `T28_C_IN_H_REL_64` +8.346526e-4. r moves the lossless h-move only at ≈ 8e-8.
- **Q_C(h)** = 1/(ωC_in(r, h)·r) ≈ 195.36 (×1 195.524). With Q_0(h) = 964.4327 (item 310 W5), the closed form 1/Q = 1/Q_0 + 1/Q_C gives Q ≈ 162.454, ratio 0.16844.
- **Q(r, h) in [161, 167], central 164.24**, which is the closed form carrying ×1's +1.10 % excess.
  - Q(r, h)/Q_0(h) in [0.166, 0.174].
  - Q(r, h)/closed form − 1 in [+0.7 %, +1.5 %] (×1 +1.0996 %; 127.74 MHz +3.5494 %).
- **The step's reading:** Q(r, h)/Q(r, ×1) − 1 is **predicted in [−4e-3, 0], central −1.84e-3**, printed beside two references:
  - Q_0's h-move at C_in: Q_0(h)/Q_0(×1) − 1 = −6.832e-3, from the two stored lossless rows;
  - T30's move at `C_F`: −6.981e-3.
  - The ratio of the Q(r) move to Q_0's is **predicted ≈ 0.27, in [0, 0.6]**. In the closed form, Q(r)'s relative h-move is 0.168 × Q_0's plus 0.832 × Q_C's (−8.34e-4, the C_in move).
- R_in(r, h) ≈ 3 166 Ω. R_in(r, h)/R_in(r, ×1) − 1 is predicted in [−6e-3, −1e-3] (×1 R_in(r) 3 175.98 Ω; lossless R_in h-move −0.78 %).
- (AB-rh) in [5e-6, 3e-5] at each point (lossless h 1.239 / 1.262 / 1.289e-5; ESR at ×1 moved them −0.3 %).
- \|f₀,A/64 − 1\| ≈ 2.59e-6, about 0.20 × `T21P_ZERO_RESOLUTION`, printed (lossless h 2.545621e-6; ESR at ×1 moved +1.9 %).
- f₀,B(C_in(r, h))/64 MHz − 1 within ±1e-9, by construction. r_AB = f₀,A/f₀,B − 1 in [−5e-6, −1e-6].
- W1 in 330–450 s at `-n 12` (≈ 560 s on a ×1.4-slow box); cells 140 270.

**Price** (measured):
- **h-rung drives at `-n 12`:** 98.2 / 105.5 / 105.7 s, mesh 32.7–34.0 s (`20261008T111513Z_ANS-7-item310-W2-p63.9.log`, `…111904Z…-W3-p64.0.log:2288`, `…112308Z…-W4-p64.1.log:1441`).
- **(W) with its ×1 mesh:** 28.0 s at `-n 12` (`…111904Z…:1430`).
- **Pytest overhead:** item 310's point windows spent ≈ 74 s outside the fixture with the full `tests/environment` (216 s footer against a 139.6 s fixture, `…111513Z…:1442`, `:2250`). W1 puts only `test_complex_mode.py` on the path, so ≈ 15–25 s is predicted.
- **W1 ≈** 20 + 28 + 34 + 3 × 103 ≈ 390 s.
- **Re-gates:** W0 is item 311 W2 / item 387 W1 (260 / 268 s, `…094028Z…:2677`, `…003929Z…:2677`). W0b is item 310 W5 (70 s, `…112649Z…:282`) and item 378's W-gate (97 s).

**Windows** (each bare, whole and foreground; Bash `timeout: 660000`; `-s`; `FEM_EM_REQUIRE_COMPLEX=1`; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE`, `FEM_EM_ANS7_EXACT_MHZ` and `TUNE_MHZ` unset). No edit after W0 starts (rule (i)), so all three windows run the same bytes. Each window is under 10 min, so §5.1's capture idiom is not required. Never two heavy windows at once.
0. **Orphan check:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
1. **W0, the ×1 64 MHz ESR re-gate** (the edit touches `tune_esr` / `_esr_ohm`): item 311's W2 command verbatim, label `ANS-7-item395-W0-esr64`, `-n 8`, ≈ 265 s.
   - **Must reproduce to the printed digits:**
     - C_in(r) 6.359301781809668e-11 F (pure, so full repr; item 387 reproduced it);
     - +4.602610e-05;
     - (AB-r) 1.233e-5 / 1.258e-5 / 1.283e-5;
     - (T2-r) −2.583111e-6 (0.201 bands);
     - r_AB −2.583113e-6;
     - Q(r) **164.543** (six significant figures, **not** full repr: it moved 2.7e-10 relative between item 311 and item 387 on solve round-off);
     - (N2-r) 494.39×;
     - (N-L) 4.97×.
   - (W) must only stay ≤ `NODE_S_RTOL`.
   - It rewrites item 311's two ×1 ESR JSONs at `-n 8` with the same printed digits. That is the accepted exception, as for item 387 W1. Item 378's (R-T-r) reads them at `NODE_S_RTOL`.
   - **Red ⇒ revert, BLOCKED.**

   `scripts/testing/run_and_log.sh ANS-7-item395-W0-esr64 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_ESR_OHM=0.2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
2. **W0b, the h lossless path unchanged (pure)** (the edit touches `_tune_point` / `tune_h`): item 310's W5 command verbatim, label `ANS-7-item395-W0b-assemble64h`, `-n 2`, `timeout -k 30 180`, ≈ 70–97 s.
   - **Must reproduce** (pure, so full digits as printed):
     - (T2‴) 2.545621e-6 (0.198 ×);
     - r(h) −2.545621e-6;
     - (AB‴) 1.239e-5 / 1.262e-5 / 1.289e-5;
     - C_in(h)/C_in(×1) − 1 = +8.346526e-4;
     - Q 964.432729560342.
   - It rewrites `tune_inmodel_cin_p2_h0.01.json` from the same point files, byte for byte as items 365 W0b and 378 W-gate did. That is the file W1 reads as Q_0(h), and (F-Q0) checks it.
   - **Red ⇒ revert, BLOCKED.**

   `scripts/testing/run_and_log.sh ANS-7-item395-W0b-assemble64h "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_TUNE_POINT=assemble timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
3. **W1, the ESR tune on the rung (closing):** label `ANS-7-item395-W1-esr64-h001`, `-n 12`, ≈ 390 s.

   `scripts/testing/run_and_log.sh ANS-7-item395-W1-esr64-h001 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_ESR_OHM=0.2 timeout -k 30 590 mpiexec -n 12 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
4. **Commit** the module, the three logs, their test-results rows (appended by hand) and the §9 marking together. `git diff --stat` on record files must show insertions only.

**Traps already paid for.**
- Items 217 / 228 / 262 / 283 / 310 / 311 / 365 / 387's lists, whole:
  - κ `None`;
  - no `ans:7` call;
  - `z_from_s(S)`, never `z_matrix`;
  - a real-scalar secant;
  - rank 0 writes, then a barrier;
  - a killed window is not a result: run `pgrep -c python3`, report it, and never lengthen the timeout;
  - never `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ` (`_nodes_env` validates against the ×1 nodes).
- **A missed ×1 literal on the h path** (item 228 code step 4; item 310's "a missed ×1 literal shows up as a tie or grid failure"). There are four here, and one of them is **silent**:
  - the lossless comparand `jl` would read the ×1 Q_0 971.07 in place of the rung's 964.43, and every Q print would be wrong with no failure. (F-Q0) exists for this;
  - `_fixture_tie()` would mesh ×1 and assert the ×1 record;
  - `_exact_geo` without `res` would drive the ×1 mesh, which the cells check catches;
  - the output names would overwrite item 311's ESR files, which item 378's (R-T-r) and item 391's re-pose read.
- **The band-transplant trap** (item 387 code step 2, in reverse): `_zero_resolution()` returns `T21P_ZERO_RESOLUTION` whenever `TUNE_MHZ` is unset. Asserting it on the rung would transplant the ×1 grid's band onto a grid it was never registered for (310-Q1). The h path asserts the crossing only.
- **Run-to-run solve round-off** (items 387 / 391): S11_A moves at ≈ 1e-12 relative between two solves of the same operator, and solve-derived rows (Q, f₀,A) at ≈ 1e-10.
  - Never demand bitwise equality across two solves.
  - W0 "reproduces to the printed digits", meaning the precision listed above, not the full repr. Item 391's W-reg was blocked by exactly this letter.
  - Only pure quantities (C_in(r) from a stored `z64`, and W0b's assembly) are compared at full digits.
- **`_tune_point()` and `tune_h` under the ESR knob** (code steps 2–3): without them the window raises at collection of the first `tune_inmodel` / `tune_h` test, before any drive.
- **`COLUMNS=cu` is required on the rung** (item 214's nodes are `_cu`-suffixed). Code step 1 makes the omission raise.
- **Free the ×1 geo before meshing the rung** (`tune_h`'s `del geo1, specs1`); ≈ 36 GiB summed at `-n 12`.
- **Price on the h drive** (≈ 100 s at `-n 12`), never on the ×1 ESR drive (≈ 45 s at `-n 8`).
- **Never enter 0.2 Ω as a parallel R** (ESR ruling (iii)). **The feed's branch carries r too** (ruling 311-Q5). `_route_a_solve_esr` / `_route_b_z_esr` are reused unchanged.
- **Never write** `tune_p2_h0.01.json`, the h point files, item 214's node files, any ×1 tune file or item 311's ESR files. The accepted exceptions are W0 rewriting item 311's ESR files and W0b rewriting the h assembly.
- **r never re-reads a lossless record.** Q_0(h) and C_in(0, h) are read and never rewritten by W1.
- Exact rectangles only; never a sub-sheet out of a box.

**Scope.** Covered: the h = 0.01 exact-rectangle rung (140 270 cells), degree 2, κ off, Cu, 64.000 MHz, r = 0.2 Ω uniform on all four capacitors on both sides, Route B on item 214's five stored h nodes. **Not claimed:**
- an h-statement for Q(r), Q_0 or C_in(r): this is one rung beside ×1. Q's h-statement stays T30 / F31's, and h = 0.005 has no C_in (T29 inferred it);
- any registered zero band for the rung (310-Q1's gap stands);
- any absolute tuned C, f₀, R_in or Q;
- 127.74 MHz on the rung with ESR;
- a measured ESR or its frequency dependence;
- any AED or Palace comparison.

**Status it can move.**
- **No glyph.** The `ANS-7` Done-when (PROJECT_PLAN.md :3501, re-stated by the 2026-10-10 weekly) is "✅ when T30 / F31 write Q's h-statement; the row holds 🟡 on the Q clause". One ESR rung does not write that statement, so the row stays 🟡.
- **Rule (j):** declined, as for items 310 / 311 / 387.
- The row gains a "step 11c (T37, item 395, record)" clause:
  - C_in(r, h)/C_in(r, ×1) − 1 beside +8.347e-4;
  - Q(r, h)/Q(r, ×1) − 1 beside Q_0's −0.683 % (at C_in) and T30's −0.698 % (at `C_F`);
  - Q(r, h)/Q_0(h) beside the closed form;
  - (AB-rh), and the (T2-rh) crossing with its offset in `T21P_ZERO_RESOLUTION` units (scale only).
- §10 T37 is recorded DONE (provisional; the weekly ratifies).
- The weekly reads the Q(r)-vs-Q_0 h-sensitivity ratio beside F31 (`xl` 10-15) when it decides which Q carries the h-statement.

**What a negative result means.**
- **W0 or W0b red:** the edit changed an existing path. Revert it and mark BLOCKED.
- **(W) red:** BLOCKED (wiring).
- **(F) or (F-Q0) red:** the inputs are not item 310's rung. Report the mismatching field and mark BLOCKED. Never re-solve a node or re-run item 310's W1.
- **(AB-rh) red with W0, W0b and (W) green:** the series-RC sheet and the circuit's (r + C) termination disagree on the rung, although both ×1 ESR and h lossless agree. Print the three residuals, open a known-issues entry, park on `attempt/ANS-7-step11c-<UTC>` and mark BLOCKED. Never apply κ, never adjust r, never re-band.
- **No crossing in (63.9, 64.1), or a + → − change:** BLOCKED. Report the Im `Y_in` triple and C_in(r, h).
- **W1 killed at 590 s:** not a result. Run the orphan check first, commit nothing for W1, park the code on `attempt/ANS-7-step11c-<UTC>` and mark BLOCKED with the unblock condition "the review re-poses W1 as item 310's per-point split under the ESR knob". Never raise the timeout.
- **The step's question answered "no"** (\|Q(r) move\| ≥ \|Q_0 move\|, or a positive Q(r) move): a reading, recorded as measured, with no disposition. The weekly reads it.
- **A printed prediction missed:** recorded as a reading, with no disposition (item 311's (N-L) rule).
- No `PREDICTED RED:` line: every assertion is predicted green, and the zero half carries no band.

**Spawn prompt (if the slot delegates; the slot may also execute directly):** "Use the implementer agent to execute `ANS-7` step 11c (§9 item 395) per `docs/planning/chunks/ANS-7.md` § '§9 item 395'. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run every harness window foreground with Bash `timeout: 660000`; never background a window or end your turn while one is running. Run the orphan check first. Edit only `tests/validation/test_ans7_exact_strip.py`, and make no edit after W0 starts. Every window runs with `-s`, `FEM_EM_REQUIRE_COMPLEX=1` and `FEM_EM_ANS_SKIP_PRIVATE=1`; `FEM_EM_ANS7_EXACT_ROUTE`, `FEM_EM_ANS7_EXACT_MHZ` and `FEM_EM_ANS7_EXACT_TUNE_MHZ` stay unset. Never write `tune_p2_h0.01.json`, any node file, any lossless tune file or item 311's ESR files except through W0 / W0b's verbatim commands. Never open `docs/private/` or any `aed_results/`. Report every printed reading, met or NOT met."

**Outcome 2026-10-12 (09:30 UTC slot, implementer): DONE (record; no glyph).** Tests-only edit per code steps 1–6 (`_esr_ohm` third legal set; `_tune_point` / `tune_h` under the ESR knob; `tune_esr` h path with `jl` → `tune_inmodel_cin_p2_h0.01.json`, (F-Q0), tie from `tj["tie"]`, (W) on the ×1 mesh freed before the rung, outputs `tune_p2_h0.01_esr0.2.json` / `tune_inmodel_cin_esr0.2_p2_h0.01.json`; six record literals incl. `T26_R_IN_R_X1` from `…094028Z…:2130`). (T2-rh)'s "no down-crossing" is asserted as `_h_assemble`'s sign-change counts (`n_up == 1`, `n_dn == 0`), since `parallel_rows` returns no `crossings_down_hz`. **W0** (`20261012T094112Z_ANS-7-item395-W0-esr64.log`, `-n 8`, 257 s, 15 passed): every listed digit reproduced — C_in(r) 6.359301781809668e-11, +4.602610e-05, (AB-r) 1.233 / 1.258 / 1.283e-5, (T2-r) −2.583111e-6 (0.201), r_AB −2.583113e-6, Q(r) 164.543 (full repr 164.5434331272106, round-off), (N2-r) 494.39×, (N-L) 4.97×; (W) 2.121e-14. **W0b** (`20261012T094542Z_ANS-7-item395-W0b-assemble64h.log`, `-n 2`, 79 s, 15 passed): (T2‴) 2.545621e-6 (0.198 ×), r(h) −2.545621e-6, (AB‴) 1.239 / 1.262 / 1.289e-5, +8.346526e-4, Q 964.432729560342 — all reproduced. **W1** (`20261012T094711Z_ANS-7-item395-W1-esr64-h001.log`, `-n 12`, 360 s pytest / 358.1 s fixture: rung mesh 32.8 s, (W) 2.5 s, drives 98.3 / 99.1 / 102.3 s; 7 passed): asserted — (F) cells 140 270 = node cells = `tune_p2_h0.01.json`, `_assert_tie(res=H_RES)` incl. (N-wire); (F-Q0) c_in and cells equal; (W) 1.373e-13 ≤ `NODE_S_RTOL`; (AB-rh) 1.235e-5 / 1.260e-5 / 1.285e-5 ≤ `REDUCTION_BAND`; (T2-rh) one − → + (f₀,A 63 999 834.33 Hz), none + → −. Printed (all predictions met unless stated): C_in(r, h) 6.364610922203437e-11 F, C_in(r, h)/C_in(0, h) − 1 = +4.623552e-5; C_in(r, h)/C_in(r, ×1) − 1 = +8.348622e-4 (beside `T28_C_IN_H_REL_64` +8.346526e-4; (N-h) met); **(N-r0h) NOT met** — ESR-path root at r = 0 from C_mode 6.364316664701223e-11 vs stored 6.364316664701229e-11, rel −9.992e-16 (not bit-equal on the rung's `z64`; a reading); offset f₀,A/64 MHz − 1 = −2.588573e-6 = 0.201 × `T21P_ZERO_RESOLUTION` (scale only); f₀,B/64 MHz − 1 = −5.657e-13; r_AB −2.588572e-6; Q(r, h) 164.21575983485775, Q(r, h)/Q_0(h) 0.170272, Q_C(h) 195.3612, closed form 162.4536 (ratio 0.168445), excess +1.0847 % (×1 +1.0996 %); **step reading** Q(r, h)/Q(r, ×1) − 1 = −1.991e-3 beside Q_0's −6.832e-3 (at C_in) and T30's −6.981e-3 (at `C_F`), ratio 0.2915 — Q(r) moves with h ≈ 0.29 × as much as Q_0; R_in(r, h) 3 166.748 Ω, R_in(r, h)/R_in(r, ×1) − 1 = −2.907e-3; (N-Lh) 4.97 × `REDUCTION_BAND` (in [3, 7], not separating at 10 ×, as predicted), Z_in form 4.849; (N2-rh) 494.23 ×. One rung beside ×1: not an h-statement; no band registered for the rung (310-Q1 stands).

**Review ruling 2026-10-12 03:00 (daily) on items 387 and 391.**
(387) W1's 64 MHz re-gate is green as the executor read it: "reproduce to the printed digits" means the printed precision; Q(r)'s full repr moving 2.7e-10 relative with S11_A at 1e-12 is solve round-off. The +3.55 % excess of Q(r)/Q_0 over the closed form at 127.74 MHz (64 MHz: +1.10 %; prediction [0, +3 %] missed) is a reading for the weekly, not a red.
(391) Attempt 1 of the step-13b family is BLOCKED on the brief's W-reg letter, not on the code: exact field equality across two solves cannot hold (run-to-run round-off ≈ 1e-12 on S11_A, ≤ 1.75e-10 in f₀-derived fields), and W-reg was owed only because the 127.74 MHz column was rendered by editing `_write_comparison`, which `XT64R` executes. **Re-posed as item 392 (attempt 2 of 4), same target, the brief's own preferred shape.** No band, tolerance or target moves.

### §9 item 392 — 2026-10-12 (chain T36 re-pose, daily provisional under §9 rule (6); example)

*Written by the 2026-10-12 03:00 daily review. Item 391's brief above stands in full except where this section amends it; read both.*

**Same target as item 391** — `ANS-7` step 13b: stage `XT128R` re-asserts T35's (AB-r′) / (T2-r′) and (R-T-r′) against item 387's stored run, (K), (R0′), and writes `degree2_exact["tuned_larmor_esr"]["127.74"]` beside the untouched 64 MHz entry. Executor `example-runner`; heavy, `-n 8`; ≈ 10 measured (example class) / 35 predicted min (code is mostly parked); depends on nothing (387 and 378 landed; the stored inputs exist — check, as item 391's brief says).

**Amendments to item 391's brief.**
1. **Take from `attempt/ANS-7-step13b-20261012T0517Z` by path** only what `XT128R` alone executes: the `_xt_esr128_stage` stage, its `main()` accepted-set entry and the guide's `XT128R` command / anchors. **Do not take the parked `_write_comparison` edit** (the two-column ESR table and the "f₀,A/64 MHz − 1" → "f₀,A/target − 1" relabel): every line `XT64R` executes stays byte-identical to `main` (`git diff main -- <script>` shows no hunk inside any function `XT64R` calls — print the list of changed functions in the commit).
2. **Step 6 (`COMPARISON.md`) is re-stated:** `XT128R` writes its own sibling table, "Tuned with capacitor ESR at 127.74 MHz (§10 T35 / T36)", directly after item 378's 64 MHz ESR section, through a new writer function that only `XT128R` calls; item 378's section stays byte-identical. Same rows as the 64 MHz table, no AED column.
3. **W-reg (window 1b) is not owed and not run** — say so in the commit. If, against amendment 1, any `XT64R`-executed line must change, stop and park: do not run W-reg, mark item 392 BLOCKED (attempt 2), and report which line and why.
4. **(K)** is unchanged and now binding: every pre-existing `metrics.json` key, `tuned_larmor_esr["64.0"]` included, byte-identical after `XT128R`; and `COMPARISON.md`'s pre-existing text byte-identical (`git diff` shows insertions only in that file).
5. **(N-C-S′) is asserted** — item 391's pure read backs it (item 387's `s11_a` vs the lossless `s11_a` differ 2.169e-3 / 2.182e-3 / 2.197e-3 > 1e-7 at every point; `docs/planning/chunks/ANS-7.md`, item 391 outcome). **(N-C-r′)** stays printed (separation 3.718e-8 < 1e-7) and **(N-L′)** stays printed (387 read 2.17×), both as item 391's outcome recorded.

**Windows** as item 391's 0, 1, 2, 3, 4, 5, with 1b skipped; labels `ANS-7-item392-…` (the `XT128R` command verbatim from item 391's window 2 with the label changed).

**Status it can move.** As item 391's: §10 T36 DONE (provisional), the `ANS-7` row's step-13b clause, no glyph. **What a negative result means:** as item 391's; a red here is attempt 2 of 4 in the family.

**Spawn prompt:** item 391's verbatim, with "§9 item 391" → "§9 item 392" and "run XT64R only as the brief's W-reg window when its condition holds" → "never run XT64R".

**Review ruling 2026-10-12 03:00 (daily) on the item 395 draft.**
(Q1) The closing window runs at `-n 12`, the rung's measured width (item 310), inside the project ceiling; the two re-gates keep their recorded widths. (Q2) Accepted: a window killed at 590 s parks the item BLOCKED (attempt 1 of the step-11c family) for the review to re-pose as item 310's per-point split; no pre-split now. (Q3) The rung's zero band stays unruled — only the sign change inside (63.9, 64.1) MHz is asserted, the offset printed beside `T21P_ZERO_RESOLUTION` for scale (310-Q1); registering a rung band is the weekly's. (Q4) Q_0 at C_in(h) (964.43) is the closed form's comparand; the ×1 comparands enter as record literals with log-line provenance, as drafted. (Q5) The printed-prediction range stays wide. (Q6) The re-gates' rewrites of item 311's 64 MHz ESR files and of the h lossless assembly are accepted (items 387 / 365 / 378 precedent); item 392 reads only the 127.74 MHz files, so the two items do not interact. **Order:** item 392 runs before 395 in §9; they touch different stage files.

**Outcome 2026-10-12 (04:30 CDT slot, example-runner): DONE.** `XT128R` at `-n 8`, 80 798 cells, 163 s window (mesh 23.6 s, drives 44.1 / 48.3 / 43.7 s; `20261012T093211Z_ANS-7-item392-XT128R.log`), all asserted readings met: (R0′) bitwise; (AB-r′) 1.587e-5 / 1.600e-5 / 1.614e-5 ≤ `REDUCTION_BAND`; (T2-r′) one − → + crossing, f₀,A/127.74 MHz − 1 = −2.301805e-6 (0.142 × `T22P_ZERO_RESOLUTION`); (R-T-r′) |f₀ ratio − 1| 1.275e-13, max |res diff| 3.353e-13, max S11_A rel 6.145e-13 ≤ `NODE_S_RTOL`; (N-C-S′) asserted, 2.197e6 × `NODE_S_RTOL`; (K) 26 + 19 + 2 + 1 keys identical. Printed: (N-C-r′) 3.718e-8 = 37.2 × (not asserted, 387 separation < 1e-7); (N-L′) 2.17× (not met, printed); Q(r)/Q_0 0.411537 vs closed form 0.397431 (+3.55 %); C_in(r)/C_in(0) − 1 = +1.344645e-5. Code from the parked branch by path (`_xt_esr128_stage`, `main()` accepted set, guide) plus a new `_write_comparison_esr128` that only `XT128R` calls; `_write_comparison` is byte-identical to `main` (changed functions: `_write_comparison_esr128` and `_xt_esr128_stage` new, `main` accepted-set line, module docstring); no `XT64R`-executed line changed, W-reg not owed and not run. `COMPARISON.md` diff is 17 insertions, 0 deletions; `metrics.json` 109 insertions, 0 deletions. Pin 24 passed; censuses before = after (`missing=0 broken=0`; docrefs `dead=0 guide=0 stale=84 exit=2`). Note: a later `XT64R` run regenerates `COMPARISON.md` via `_write_comparison` and would drop the sibling section until `XT128R` re-runs. Deviation: the pre-census ran after the script/guide edits (not before any file was written); the census reads artifacts, and before = after.

### §9 item 400 — 2026-10-12 (`ans:7` comparison writer keeps the 127.74 MHz ESR section; lane 1, ops)

*Written by the 2026-10-12 10:30 daily review from item 392's caveat (its attempts entry, "Next-attempt hypothesis").*

**What.** `_write_comparison` (`examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py:1155–1237`) rewrites `COMPARISON.md` wholesale from `metrics.json`. Four stages call it (`:1149`, `:1370`, `:1582`, `:1807`). It renders the 64 MHz ESR section. It does not render the 127.74 MHz sibling, which only `XT128R`'s `_write_comparison_esr128` (`:1240–1272`) inserts. So if any of those stages runs after item 392, the section disappears from the file while its data stays in `metrics.json` (known-issues "🟡 OPEN 2026-10-12 (found by `ANS-7` step 13b …)").

**Edit (one function, no solve).** At the end of `_write_comparison`, after its `write_text`, call `_write_comparison_esr128(m)` when `m["degree2_exact"].get("tuned_larmor_esr", {})` holds `"127.74"`. The inserter already removes an earlier copy and places the section before the closing "Not claimed" line, so the result is the bytes `XT64R` + `XT128R` produced. Edit nothing else in the script:
- not the stage-B2 degree-1 writer at `:517–540`, which is superseded and out of scope (say so in the outcome);
- not the inserter;
- no stage logic.

**Anchors.** These go in a new unit test, `tests/unit/test_ans7_comparison_writer.py`. It imports the script by path with `importlib.util` (grep `spec_from_file_location` under `tests/` for the pattern).
- (R) Render-only reproduction.
  - Copy the committed `metrics.json` and `COMPARISON.md` into `tmp_path`.
  - Set the loaded module's `CASE_DIR = tmp_path`.
  - Call `_write_comparison(json.loads(metrics))`.
  - Assert the rendered `COMPARISON.md` is **byte-identical** to the committed one, including the `## Tuned with capacitor ESR at 127.74 MHz (§10 T35 / T36)` section.
  - Print both lengths and the count of that section heading (1).
- (R0) The committed files are untouched: the repo's `COMPARISON.md` bytes are equal before and after the test.

**Negative controls (both asserted).**
- (N-old) The pre-change writer drops the section.
  - Load the script as of the **pinned pre-change commit `b01dc790`**, never `HEAD:` (`OPS-43`'s trap).
  - Get that copy before the window with a plain host git command (not a harness call): `git show b01dc790:examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py > logs/ans7_writer_b01dc790.py`. `logs/` is gitignored.
  - Load it from `/workspace/logs/` and run its `_write_comparison` on the same copied metrics.
  - Assert the 127.74 MHz heading count is 0 and that the bytes differ from the committed file.
- (N-key) Delete `tuned_larmor_esr["127.74"]` from the copied metrics.
  - The new writer's output then has no 127.74 MHz heading.
  - That output equals (N-old)'s output byte for byte: the call is conditional and adds nothing else.

**Tier, cost.** Smoke, no MPI. The test itself takes ≈ 5 s. The module import pulls in dolfinx, so it runs in complex mode. ≈ 3 measured (ops / docs median) / 15 predicted min. Depends on nothing.

**Window (bare, whole, foreground, Bash `timeout: 660000`):**

`scripts/testing/run_and_log.sh ANS-7-item400-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans7_comparison_writer.py -s -v --tb=short'"`

Then run the examples census and docrefs as item 392 ran them (`check_example_setup_figures.py`, `check_example_doc_references.py`). Docrefs gates on `exit != 1`. The script changed, so `ans:7`'s own pin module re-runs green in the same slot (standing rule (c)). Take the command from item 392's pin window, `20261012T093503Z_ANS-7-item392-pin.log` (24 passed, 14 s).

**Rule (i).** No stage of `ans:7` is re-run. The edit touches only the rendering path, and the closing window is (R), a pure function of the committed `metrics.json`. The outcome says "no `ans:7` stage re-run; render-only". W-reg is not owed, because nothing solves and nothing writes `metrics.json`.

**Traps already paid for.**
- Never write the repo's `COMPARISON.md` or `metrics.json` from the test; use tmp copies only.
- The pinned control commit is `b01dc790`, never `HEAD:`.
- `-k a or b` splits inside the quoted container command.
- Pytest prints need `-s`.
- No slot opens `aed_results/` or `docs/private/`. The private writer at `:743` is not called.

**Scope.** Rendering only. Not claimed: anything about the solve, the records or the 64 MHz entry. The stage-B2 writer (`:540`) still writes its own superseded file, as it always has.

**Status it can move.** The known-issues entry "🟡 OPEN 2026-10-12 (found by `ANS-7` step 13b …)" is retired in the same commit (✅ RETIRED with the W1 log). `ANS-7` stays 🟡. Its row gains nothing: it is already over 2 KB, so the outcome goes in this file.

**What a negative result means.**
- (R) red means the committed `COMPARISON.md` is not what the two writers produce from the committed metrics, i.e. a record already out of step. Stop and do not regenerate the file. Mark item 400 BLOCKED with the first differing line printed, and leave the entry open.
- If (N-old) finds the section present, the issue is mis-stated. Report it and stop.

There is no `PREDICTED RED` line.

**Executor.** The slot runs this item directly (ops). If it spawns an executor, the spawn prompt carries the §9 "Every item" bare-command sentence verbatim.

**Outcome 2026-10-12 (12:00 CDT slot, run directly): DONE.** `_write_comparison` ends by calling `_write_comparison_esr128(m)` when `tuned_larmor_esr` holds `"127.74"`. That is two lines, and nothing else in the script changed: not the stage-B2 degree-1 writer at `:517–540` (superseded, out of scope), not the inserter, no stage logic. New `tests/unit/test_ans7_comparison_writer.py`, all asserted: (R) the committed metrics render the committed `COMPARISON.md` byte for byte, 3543 B, 127.74 MHz heading count 1; (N-old) the pinned `b01dc790` writer renders 2886 B, count 0, and differs; (N-key) with the key deleted the new writer renders 2886 B, count 0, equal to (N-old) byte for byte; (R0) the repo's `COMPARISON.md` / `metrics.json` are unchanged. W1 `20261012T170114Z_ANS-7-item400-W1.log`, 5 passed, 3 s. Pin `20261012T170127Z_ANS-7-item400-pin.log`, 24 passed, 14 s. Censuses `missing=0 broken=0` (`20261012T170142Z_ANS-7-item400-figures.log`); docrefs `dead=0 guide=0 stale=90 exit=2` (`20261012T170143Z_ANS-7-item400-docrefs.log`). Known-issues guard `test_rotate_known_issues_entries.py` 8 passed after the retirement (`20261012T170223Z_ANS-7-item400-ki.log`). No `ans:7` stage re-run; render-only; W-reg not owed. **Deviation:** the pinned copy lives at `logs/ans7_writer_b01dc790/<the script's repo-relative path>`, not `logs/ans7_writer_b01dc790.py`. The module body evaluates `Path(__file__).resolve().parents[3]`, which raised `IndexError: 3` at the flat path (first W1, `20261012T170051Z_ANS-7-item400-W1.log`, a setup error before any assertion). The mirrored depth changes only the `sys.path` entry the module inserts, and `/workspace` is already on `PYTHONPATH`. Known-issues entry retired.

## §7 row moved verbatim 2026-10-12 (rotate_plan_archive.py chunks)

 The tuned copper F-small birdcage: `Z_in(f)`, f₀ and Q against AED | 🟡 Adjudicated 2026-09-30: AGREE on the mode and the copper contribution, absolute rows INCONCLUSIVE at degree 1; on exact rectangles at degree 2 (steps 5–6) the tuned `f₀` sits within 0.11 % of AED First Order on one mesh (2026-10-03 weekly). Steps 7′ / 8′ / 10 (T21′ / T22′ / T23, 2026-10-06) tune on C_in — the in-model zero sits at 64.000 / 127.74 MHz within the grid's measured resolution (−2.536e-6 / −2.265e-6 against 1.287e-5 / 1.616e-5), circuit = field at `REDUCTION_BAND`, and `ans:7` carries the tuned-to-Larmor rows — records on a ×1 mesh, not a tuning claim. **T13 decided by the 2026-10-10 weekly:** `f₀` AGREE with its h-statement (0.11 % from AED First Order at ×1, 0.03 % at h = 0.005, first order in h); `R_in` / Q move away with h (Q −3.1 %, unconverged). T29 (record): C_in moves +0.083 % at h = 0.01 (T28, item 310), ≈ +0.16 % inferred at h = 0.005. Step 11 (T26, item 311, record): with r = 0.2 Ω on all four capacitors, C_in(r)/C_in(0) − 1 = +4.60e-5, (AB-r) 1.233–1.283e-5, (T2-r) −2.583e-6 (0.201 bands), Q(r)/Q(0) 0.1694 beside the closed form's 0.1676 (+1.10 %). T30 (item 364): Cu R_in / Q at `C_F`, ×1 → 0.01 → 0.005: successive-change ratios 0.705 / 0.717; unconverged at h = 0.005 — F31 (`xl` entry 34) is the next rung. Step 12b (T34, items 365 / 365b, record, lossless): at h = 0.01 and 127.74 MHz C_in(h)/C_in(×1) − 1 = +8.968e-4 (64 MHz: +8.347e-4), (AB⁗) 1.593–1.619e-5, zero printed at 2.271e-6 (0.141 × `T22P_ZERO_RESOLUTION`, scale only); one rung, not an h-statement. Step 11b (T35, item 387, record): at 127.74 MHz, r = 0.2 Ω, C_in(r)/C_in(0) − 1 = +1.34e-5 (64 MHz +4.60e-5), (AB-r′) 1.587–1.614e-5, (T2-r′) −2.302e-6 (0.142 × `T22P_ZERO_RESOLUTION`), Q(r)/Q_0 0.4115 beside the closed form's 0.3974 (+3.55 %). Step 13 (T27, item 378, record): `ans:7`'s tracked outputs carry the 64 MHz ESR rows (`XT64R`), (AB-r), (T2-r) and (R-T-r) (≤ 5.2e-13 at `-n 8`) re-asserted through the example, the lossless rows labelled. Step 13b (T36, item 392, record): the 127.74 MHz ESR rows (`XT128R`, own comparison section), (AB-r′) ≤ 1.614e-5, (T2-r′) −2.302e-6 (0.142 bands) and (R-T-r′) ≤ 6.145e-13 re-asserted through the example at `-n 8`. Step 11c (T37, item 395, record): r = 0.2 Ω on the h = 0.01 rung at 64 MHz, C_in(r, h)/C_in(r, ×1) − 1 = +8.349e-4 (lossless +8.347e-4), Q(r, h)/Q(r, ×1) − 1 = −0.199 % beside Q_0's −0.683 % (C_in) and T30's −0.698 % (`C_F`); Q(r, h)/Q_0(h) 0.1703 vs closed form 0.1684 (+1.08 %); (AB-rh) 1.235–1.285e-5, one − → + crossing, offset 0.201 × `T21P_ZERO_RESOLUTION` (scale only). **Done-when (re-stated 2026-10-10 weekly):** ✅ when T30 / F31 write Q's h-statement; the row holds 🟡 on the Q clause. *History: `docs/planning/chunks/ANS-7.md`.* 

### §9 item 419 — 2026-10-13 (§10 chain T38; chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.*

**Licence.** §10 chain T, T38 (PROJECT_PLAN.md, grep "**T38. `ANS-7` step 14**"), quoted, not changed: "`ANS-7` step 14 — T37 at 127.74 MHz: r = 0.2 Ω on T34's h = 0.01 node set; T37's identities and prints (C_in(r, h)/C_in(r, ×1) − 1 beside +8.968e-4, Q(r, h)/Q(r, ×1) − 1, Q(r)/Q₀ against 1/Q = 1/Q₀ + 1/Q_C). A record. Heavy. Serial on nothing."

The rulings that apply:
- the ESR ruling (§9, grep "Binding ruling — capacitor ESR (operator, 2026-10-06"): r = 0.2 Ω on every capacitor, both sides at once; the feed's branch carries r too (311-Q5); lossless records stay, labelled lossless;
- the sheets ruling: exact rectangles, κ `None`;
- **365-Q1 / 310-Q1 / 395-Q3:** no zero band is registered for the h = 0.01 rung, so the zero half asserts the crossing only, and the offset is printed beside `T22P_ZERO_RESOLUTION` for scale.

**No band, tolerance, Done-when or target is set or moved; every band is imported by name.**

**Feasibility (checked 2026-10-13 against the code and the disk).** T38 is item 395's construction, with T34's 127.74 MHz node set in place of item 214's 64 MHz set. **No node is re-solved.**
- **Inputs on disk** (`logs/ans7_exact/t127.74/`):
  - T34's seven `node_{125.74,…,130.74}_p2_h0.01_cu.json`;
  - `tune_p2_h0.01.json` (C_in(0, h) = 1.5288541964455146e-11 F, `z64["cu"]` at 127.74 MHz, `tie`, cells 140 270);
  - the lossless assembly `tune_inmodel_cin_p2_h0.01.json`: Q_0(h) 606.6969050131794, R_in(0, h) 23 319.940890546182 Ω (`20261011T112734Z_ANS-7-item365b-A-assemble.log:106`);
  - item 387's ×1 ESR files `tune_p2_esr0.2.json` and `tune_inmodel_cin_esr0.2_p2.json`, which W0 rewrites.
- **`_active()` and `_tune_json()` already return T34's set** under `TUNE_MHZ=127.74` with `RESOLUTION=0.01` (`test_ans7_exact_strip.py:968–972`, `:1121–1123`). `tune_esr`'s generic h path (`h37 = _h_tune()`, `:2309`) already does what T38 needs:
  - the tie from `tj["tie"]`;
  - (W) on the ×1 mesh, then `del`;
  - the rung mesh;
  - the outputs under `act["out"]` = `t127.74/` with the `_h0.01` suffix;
  - `jl` → `t127.74/tune_inmodel_cin_p2_h0.01.json`;
  - (F-Q0).
- **Three things block or misread the 127.74-h path today:**
  1. `_esr_ohm()` (`:2198–2220`) **raises**: `h37` there requires `_tune_mhz() is None`, so `RESOLUTION` set with `TUNE_MHZ=127.74` fails `(res_raw and not h37)`.
  2. In `tune_esr` both `t35` (`_tune_mhz() == 127.74`) and `h37` are True. So the tag reads "T35". The C_in prints take the h37 branch, with 64 MHz text and the ×1 64 MHz literal `T26_C_IN_R_X1`. The (AB) / (T2) / Q / R_in / (N-L) / (N2) prints take the **t35** branch, with ×1 127.74 predictions, and the Q print compares against the rung's Q_0 under ×1 text. **This one is silent:** every print would read plausibly with the wrong comparand.
  3. The h37 Q and R_in prints divide by `T26_Q_R_X1` / `T21P_Q0_X1` / `T26_R_IN_R_X1`, which are 64 MHz literals.
- **Price:** one ESR window on the rung is (W) plus three ≈ 100 s degree-2 drives at `-n 12`. It fits one window under 590 s.

`ANS-7` step 14: tune the exact-sheet F-small birdcage on C_in with r = 0.2 Ω on both sides at 127.74 MHz, on the h = 0.01 rung (140 270 cells, degree 2, κ off, Cu, T34's seven h nodes for Route B).
- **Executor:** implementer.
- **Edits:** tests only, `tests/validation/test_ans7_exact_strip.py`. No `src/`, `ans:7`, or other module edit.
- **Tier:** heavy. W1 runs at **`-n 12`**, the rung's measured width (items 365 / 365b / 395; summed `ru_maxrss` 33.8–36.7 GiB). The re-gates run at their records' widths: W0 `-n 8`, W0b `-n 2`. Complex build.
- **Class:** chain step (solve), ≈ 12 measured (18:00 restock table).
- **Time:** ≈ 45 predicted min: 15 fixed + ≈ 12 code + W0 ≈ 4.7 + W0b ≈ 1.4 + W1 ≈ 7 + ≈ 5 spare.
- **Depends on nothing open.** Items 365 / 365b / 387 / 395 have landed. If any input named above is absent, stop and mark BLOCKED. Never re-create one.

**Code** (`tests/validation/test_ans7_exact_strip.py`). Every existing path is behaviourally identical **by construction**: every new line sits under a new `h38` flag that is False on all of them, and W0 / W0b prove it.
1. **`_esr_ohm()`:** a fourth legal set, `h38 = bool(res_raw) and _tune_mhz() == 127.74 and _resolution() == H_RES and _cols() == ("cu",)`. The guard becomes `(res_raw and not (h37 or h38))`. `TUNE_POINT` stays unset in every set, and the message names all four sets ("… or `RESOLUTION`=0.01, `TUNE_MHZ`=127.74, `COLUMNS`=cu (127.74 MHz on the h rung, T38)").
2. **`tune_esr` flags:**
   - `h37 = _h_tune()` is unchanged: True on both rungs; its non-print logic is shared.
   - New `h38 = h37 and _tune_mhz() == 127.74`.
   - `t35 = _tune_mhz() == 127.74 and not h37`. This is unchanged on the ×1 127.74 path, where `h37` is False.
   - `tag = "T38" if h38 else ("T35" if t35 else ("T37" if h37 else "T26"))`.
3. **Every print chain** in `tune_esr` that has an `h37` or `t35` branch gets **`if h38:` first**, and the existing first branch becomes `elif`. The chains are: C_in, (AB), (T2), Q, R_in, (N-L), (N2). The h38 texts carry this brief's predictions and the 127.74 literals below, and no 64 MHz or ×1-127.74 prediction is inherited (item 365b C1's two inherited "NOT met" lines). On every other path the branch taken and its bytes are unchanged.
   - **The single-branch `if h37:` prints also need the tag.** The (F-Q0) print (`:2361`) and the rung-mesh print (`:2399`) hard-code "T37". Their literal becomes `{tag}`, which prints "T37" on the 64 MHz rung byte for byte and "T38" here.
   - The (T2) h38 branch prints the offset beside `T22P_ZERO_RESOLUTION`, not the h37 branch's `T21P_ZERO_RESOLUTION`.
4. **Record literals**, printed on the h38 path only and never asserted, each with its provenance comment:
   - `T35_C_IN_R_X1 = 1.5275048999871143e-11` (`20261012T004412Z_ANS-7-item387-W2-esr128.log:291`);
   - `T35_Q_R_X1 = 254.46925911792314` (`…004412Z…:2130`);
   - `T35_R_IN_R_X1 = 9790.703605847839` (`…004412Z…:2130`);
   - `T35_Q_CLOSED_EXCESS_X1 = 3.5494e-2` (`…004412Z…:2133`);
   - `T22P_Q0_X1 = 618.338591513781` (Q_0 at the ×1 127.74 C_in, `…004412Z…:2133`);
   - `T34_C_IN_H_REL_128 = 8.967920e-4` (`20261011T112734Z_ANS-7-item365b-A-assemble.log:113`), the "+8.968e-4" §10 names.
   - Literals, not reads, because W0 rewrites item 387's ×1 ESR JSONs in this slot (item 395's reason).
5. **The tests are unchanged.** `test_t26_ab_circuit_equals_field_esr` already asserts the rung's tie, (F) and (F-Q0) whenever `tune_esr["h"]`. `test_t26_inmodel_zero_esr`'s h branch already asserts the crossing only: one − → + via `n_up` / `n_dn` and lo < f₀ < hi. Both are generic in `tune` / `tgt`. Their docstrings gain "and T38 at 127.74 MHz (§9 item 419)".
6. **No output-name edit.** On the h38 path `act["out"]` is `t127.74/`, so the writes go to `t127.74/tune_p2_h0.01_esr0.2.json` and `t127.74/tune_inmodel_cin_esr0.2_p2_h0.01.json`. These are new names, and nothing existing is touched.

**Anchors** (asserted; every band imported by name), at 127.64 / 127.74 / 127.84 MHz on the rung:
- **(AB-rh′) circuit = field at C_in(r, h):** \|S11_A − S11_B\|/\|S11_B\| ≤ `REDUCTION_BAND` (1e-3, `tests/validation/test_port_lumped_rlc_termination.py:88`, imported through `test_port_circuit_layer_field` at `test_ans7_exact_strip.py:81`). Route B is T34's seven h nodes through `bary_interp` and `_route_b_z_esr` (the series r + C on P2–P4, and the feed's (r + C) in parallel).
- **(T2-rh′) the in-model zero (310-Q1's form):** Im `Y_in`,A changes sign − → + exactly once and never + → − in (127.64, 127.84) MHz, with f₀,A inside the bracket.
  - \|f₀,A/127.74 MHz − 1\| is **printed** beside `T22P_ZERO_RESOLUTION` (1.616e-5, `test_ans7_exact_strip.py:1025`) for scale only, never asserted (365-Q1).
- **(F) the rung:** drive-mesh cells = the seven node files' cells = `tune_p2_h0.01.json`'s (140 270). `_assert_tie(…, res=H_RES)` passes, including its (N-wire) clause: off the ×1 record by > `CELL_COUNT_BAND` (0.01, `tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`).
- **(F-Q0):** the lossless comparand `t127.74/tune_inmodel_cin_p2_h0.01.json` has `c_in` and `cells` equal to the rung's.
- **(W):** ≤ `NODE_S_RTOL` (1e-9, `tests/validation/test_ans7_sparse_route_b.py:87`), unchanged.

**Negative controls (rule (e)).**
- **(W): asserted.** It is the same degree-1 ×1 64 MHz drive and comparison as items 310 W3 (1.720e-13), 387 W2 (1.062e-13) and 395 W1 (1.373e-13). Run-to-run round-off means it is not digit-stable; it must only stay ≤ `NODE_S_RTOL`.
- **(N-wire): asserted** inside `_assert_tie(res=H_RES)`. The same comparison on the same rung was green in items 365 / 365b / 395.
- **(N2-rh′): predicted, printed.** The circuit at 0.5·C_in(r, h) against the in-model at C_in(r, h), at 127.74 MHz, is ≥ 10 × `REDUCTION_BAND`. Earlier readings: lossless h 494.00× (`…365b-A-assemble.log:112`); ×1 ESR 494.07× (item 387).
- **(N-Lh′): predicted, printed — predicted NOT to separate at 10×.** The lossless circuit at C_in(0, h) against the in-model at (r, C_in(r, h)) on the bare S11 is predicted in [1, 4] × `REDUCTION_BAND` (×1 127.74: 2.17×). The Z_in form \|dZ\|/\|Z\| is predicted ≈ Q_0(h)/Q_closed − 1 = 1.489, in [1.3, 1.7] (×1 1.419).
- **(N-r0h′): predicted, printed.** `_c_in_esr(z64(h), 127.74 MHz, C_mode(h), 0)` equals `tune_p2_h0.01.json`'s `c_in` bit for bit, or within 1e-15 relative. At 64 MHz on the rung it read −9.99e-16, not bit-equal (item 395).
- **(N-h′): predicted, printed.** \|C_in(r, h)/C_in(r, ×1) − 1\| ≥ 10 × the secant tolerance (1e-11).

**Predictions** (written first, printed, never asserted; the drafter's arithmetic from the logged values above).
- **C_in(r, h)/C_in(0, h) − 1 in [+0.8e-5, +2.0e-5], central +1.351e-5.** This is ×1's +1.344645e-5 scaled by the 64 MHz h / ×1 ratio, 4.623552e-5 / 4.602610e-5.
- **C_in(r, h)/C_in(r, ×1) − 1 in [+8.6e-4, +9.3e-4], central +8.9685e-4**, printed beside `T34_C_IN_H_REL_128` +8.967920e-4. r moves the lossless h-move only at ≈ 6e-8.
- **The closed form.** Q_C(h) = 1/(ωC_in(r, h)·r) ≈ 407.47 (×1 407.83). With Q_0(h) = 606.697, 1/Q = 1/Q_0 + 1/Q_C gives Q ≈ 243.756, ratio 0.40178.
- **Q(r, h) in [246, 258], central 252.41**, which is the closed form carrying ×1's +3.55 % excess.
  - Q(r, h)/Q_0(h) in [0.405, 0.425].
  - Q(r, h)/closed form − 1 in [+2.5 %, +4.5 %] (×1 +3.5494 %; at 64 MHz the excess moved −0.015 pp from ×1 to h).
- **The step's reading:** Q(r, h)/Q(r, ×1) − 1 is **predicted in [−1.2e-2, −4e-3], central −8.10e-3**, printed beside Q_0's h-move at C_in, Q_0(h)/Q_0(×1) − 1 = −1.8827e-2 (`…365b-A-assemble.log:110`).
  - The ratio of the Q(r) move to Q_0's is **predicted ≈ 0.43, in [0.2, 0.7]**. In the closed form, Q(r)'s relative h-move is 0.402 × Q_0's plus 0.598 × Q_C's (−8.97e-4).
  - At 64 MHz the ratio read 0.2915 (item 395). At 127.74 MHz Q_0 carries a larger share of 1/Q, so the ratio is predicted larger.
- **R_in(r, h) ≈ 9 703 Ω.** R_in(r, h)/R_in(r, ×1) − 1 is predicted in [−1.4e-2, −4e-3], central −9.0e-3: the Q move less the C_in move, as at 64 MHz, where −1.99e-3 − 8.35e-4 predicted the measured −2.907e-3.
- (AB-rh′) in [8e-6, 4e-5] at each point (lossless h 1.593 / 1.605 / 1.619e-5; ×1 ESR 1.587 / 1.600 / 1.614e-5).
- \|f₀,A/127.74 − 1\| ≈ 2.29e-6, ≈ 0.142 × `T22P_ZERO_RESOLUTION`, in [1e-6, 5e-6] (lossless h 2.271224e-6; ×1 ESR 2.301805e-6). Printed.
- f₀,B(C_in(r, h))/127.74 MHz − 1 within ±1e-9, by construction. r_AB in [−5e-6, −1e-6].
- W1 in 340–450 s at `-n 12`; cells 140 270.

**Price** (measured):
- **h-rung windows at 127.74, `-n 12`:**
  - C1: 260 s (`20261011T111130Z_ANS-7-item365b-C1-cin128-h001.log`);
  - D1 / D2 / D3: 209 / 239 / 215 s, one drive each plus the mesh and the full `tests/environment`.
- **The 64 MHz ESR rung window (item 395 W1, `20261012T094711Z_ANS-7-item395-W1-esr64-h001.log`):** 360 s, with rung mesh 32.8 s, (W) 2.5 s and drives 98.3 / 99.1 / 102.3 s.
- **W1 ≈** 20 (pytest with only `test_complex_mode.py`) + 28 (the ×1 (W) mesh) + 34 + 3 × 105 ≈ 400 s.
- **Re-gates:**
  - W0 is item 387 W2, 279 s at `-n 8` (`20261012T004412Z_ANS-7-item387-W2-esr128.log:2677`);
  - W0b is item 365b A, 80 s at `-n 2` (`20261011T112734Z_ANS-7-item365b-A-assemble.log`).

**Windows** (each bare, whole and foreground; Bash `timeout: 660000`; `-s`; `FEM_EM_REQUIRE_COMPLEX=1`; `FEM_EM_ANS_SKIP_PRIVATE=1`; `ROUTE` and `FEM_EM_ANS7_EXACT_MHZ` unset). No edit after W0 starts (rule (i)), so all three windows run the same bytes. Each window is under 10 min, so §5.1's capture idiom is not required. Never two heavy windows at once.
0. **Orphan check:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
1. **W0, the ×1 127.74 ESR re-gate** (the edit touches `_esr_ohm` and the `t35` definition): item 387's W2 command verbatim (`…004412Z…:12`), label `ANS-7-item419-W0-esr128`, `-n 8`, ≈ 280 s.
   - **Must reproduce to the printed digits:**
     - C_in(r) 1.5275048999871143e-11 (pure, so full repr);
     - +1.344645e-05;
     - (AB-r′) 1.587e-5 / 1.600e-5 / 1.614e-5;
     - (T2-r′) −2.301805e-06 (0.142 bands);
     - r_AB −2.301804e-6;
     - Q(r) 254.469 (six significant figures, not full repr: run-to-run round-off);
     - (N2-r′) 494.07×;
     - (N-L′) 2.17×;
     - (N-r0′) bit-equal.
   - (W) must only stay ≤ `NODE_S_RTOL`.
   - It rewrites item 387's two ×1 ESR JSONs with the same printed digits. That is the accepted exception, as item 395 W0 was for item 311's. Item 392's `XT128R` (R-T-r′) reads them at `NODE_S_RTOL`.
   - **Red ⇒ revert, BLOCKED.**

   `scripts/testing/run_and_log.sh ANS-7-item419-W0-esr128 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_ESR_OHM=0.2 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
2. **W0b, the 127.74 h lossless path unchanged (pure)**: item 365b A's command verbatim (`…112734Z…:12`), label `ANS-7-item419-W0b-assemble128h`, `-n 2`, `timeout -k 30 180`, ≈ 80 s.
   - **Must reproduce** (pure, so full digits as printed):
     - (T2⁗) Im `Y_in`,A −4.061501062e-05 / +1.181776383e-07 / +4.088914063e-05 S;
     - \|f₀,A/127.74 − 1\| 2.271224e-6 (0.141 ×);
     - r(h) −2.271224e-6;
     - (AB⁗) 1.593e-5 / 1.605e-5 / 1.619e-5;
     - Q 606.6969050131794;
     - C_in(h)/C_in(×1) − 1 = +8.967920e-4.
   - It rewrites `t127.74/tune_inmodel_cin_p2_h0.01.json` from the same point files, byte for byte (the items 365 W0b / 395 W0b precedent). That is the file W1 reads as Q_0(h), and (F-Q0) checks it.
   - **Red ⇒ revert, BLOCKED.**

   `scripts/testing/run_and_log.sh ANS-7-item419-W0b-assemble128h "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_POINT=assemble timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
3. **W1, the ESR tune on the rung at 127.74 MHz (closing):** label `ANS-7-item419-W1-esr128-h001`, `-n 12`, ≈ 400 s.

   `scripts/testing/run_and_log.sh ANS-7-item419-W1-esr128-h001 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS7_EXACT_DEGREE=2 FEM_EM_ANS7_EXACT_COLUMNS=cu FEM_EM_ANS7_EXACT_RESOLUTION=0.01 FEM_EM_ANS7_EXACT_TUNE=inmodel FEM_EM_ANS7_EXACT_TUNE_MHZ=127.74 FEM_EM_ANS7_EXACT_TUNE_C=cin FEM_EM_ANS7_EXACT_ESR_OHM=0.2 timeout -k 30 590 mpiexec -n 12 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans7_exact_strip.py -v -s --tb=short'"`
4. **Commit** the module, the three logs, their test-results rows (appended by hand) and the §9 marking together. `git diff --stat` on record files must show insertions only.

**Traps already paid for.**
- Items 217 / 228 / 262 / 283 / 310 / 311 / 365 / 387 / 395's lists, whole:
  - κ `None`;
  - no `ans:7` call;
  - `z_from_s(S)`, never `z_matrix`;
  - a real-scalar secant;
  - rank 0 writes, then a barrier;
  - a killed window is not a result: run `pgrep -c python3`, report, and never lengthen the timeout;
  - never `ROUTE=B` / `FEM_EM_ANS7_EXACT_MHZ`.
- **The silent misread (code steps 2–4).** Without `h38`, the 127.74-h window would print "T35" and divide the rung's Q by ×1 literals. It would also print 64 MHz text (`T26_C_IN_R_X1`) under the C_in line, and nothing would fail. Grep the W1 log for `T35` and `T37` tags: there must be none. A stray tag is a missed branch, so fix it before W0 starts.
- **The band-transplant trap** (310-Q1 / 365-Q1 / 395-Q3): `_zero_resolution()` returns `T22P_ZERO_RESOLUTION` under `TUNE_MHZ=127.74`. Asserting it on the rung would transplant the ×1 grid's band. The h branch of `test_t26_inmodel_zero_esr` returns before the band line; keep it so.
- **`.2f` names at 127.74** (item 365's trap): this item writes no per-point file. The rung's outputs take only `H_SFX` and the ESR suffix (code step 6).
- **Run-to-run solve round-off** (items 387 / 391 / 395): S11_A moves ≈ 1e-12 relative between two solves of the same operator, and the solve-derived rows (Q, f₀,A) ≈ 1e-10. W0 "reproduces to the printed digits", meaning the precision listed above. Only pure quantities are compared at full digits: C_in(r) from a stored `z64`, and W0b's assembly.
- **`COLUMNS=cu` is required** (T34's nodes are `_cu`-suffixed); code step 1 makes its omission raise.
- **Free the ×1 geo before meshing the rung** (`tune_esr`'s `del geo1, specs1`); ≈ 36 GiB summed at `-n 12`.
- **Price on the h drive** (≈ 100 s at `-n 12`), never on the ×1 ESR drive (≈ 48 s at `-n 8`).
- **Never enter 0.2 Ω as a parallel R** (ESR ruling (iii)). **The feed's branch carries r too** (ruling 311-Q5). `_route_a_solve_esr` / `_route_b_z_esr` are reused unchanged.
- **Never write** `tune_p2_h0.01.json`, T34's node or point files, any lossless tune file, item 311's or item 395's ESR files. The accepted exceptions are W0 rewriting item 387's ESR files and W0b rewriting the 127.74 h assembly.
- **r never re-reads a lossless record.** Q_0(h) and C_in(0, h) are read and never rewritten by W1.
- **Docker socket:** every harness call is the whole Bash command, bare and repo-relative, and an executor that loses the socket stops and reports.

**Scope.** Covered: the h = 0.01 exact-rectangle rung (140 270 cells), degree 2, κ off, Cu, 127.74 MHz, r = 0.2 Ω uniform on all four capacitors on both sides, Route B on T34's seven stored h nodes. **Not claimed:**
- an h-statement for Q(r), Q_0 or C_in(r) at 127.74 MHz: this is one rung beside ×1;
- any registered zero band for the rung (365-Q1's gap stands);
- any absolute tuned C, f₀, R_in or Q;
- a measured ESR or its frequency dependence;
- any AED or Palace comparison (there is no AED ESR case, and the Palace twin is T39 / T40 at 64 MHz).

**Status it can move.**
- **No glyph.** The `ANS-7` Done-when (§7 row, re-stated 2026-10-10 weekly) is "✅ when T30 / F31 write Q's h-statement; the row holds 🟡 on the Q clause". One ESR rung at a second frequency does not write it. **Rule (j):** read; no mismatch; declined as for items 310 / 311 / 365 / 387 / 395.
- The row gains a "step 14 (T38, item 419, record)" clause:
  - C_in(r, h)/C_in(r, ×1) − 1 beside +8.968e-4;
  - Q(r, h)/Q(r, ×1) − 1 beside Q_0's −1.883 % at C_in;
  - the Q(r)-vs-Q_0 h-sensitivity ratio beside 64 MHz's 0.2915;
  - Q(r, h)/Q_0(h) beside the closed form;
  - (AB-rh′), and the crossing with its offset in `T22P_ZERO_RESOLUTION` units (scale only).
- §10 T38 is recorded DONE. The 10-17 weekly reads the two frequencies' ratios beside F31 when T41 decides which Q carries the h-statement.

**What a negative result means.**
- **W0 or W0b red:** the edit changed an existing path. Revert it, BLOCKED.
- **(W) red:** BLOCKED (wiring).
- **(F) or (F-Q0) red:** the inputs are not T34's rung. Report the mismatching field, BLOCKED. Never re-solve a node or re-run item 365b's C1.
- **(AB-rh′) red with W0, W0b and (W) green:** the series-RC sheet and the circuit's (r + C) termination disagree on the rung at 127.74, although both ×1 ESR at 127.74 and h lossless at 127.74 agree. Print the three residuals, open a known-issues entry, park on `attempt/ANS-7-step14-<UTC>` and mark BLOCKED. Never apply κ, never adjust r, never re-band.
- **No crossing in (127.64, 127.84), or a + → − change:** BLOCKED. Report the Im `Y_in` triple and C_in(r, h).
- **W1 killed at 590 s:** not a result. Run the orphan check first, commit nothing for W1, park the code on `attempt/ANS-7-step14-<UTC>`, and mark BLOCKED with the unblock condition "the review re-poses W1 as item 310's per-point split under the ESR knob". Never raise the timeout.
- **The step's reading lands outside its predicted range** (for example \|Q(r) move\| ≥ \|Q_0 move\|, or a positive move): a reading, recorded as measured, with no disposition. The weekly reads it.
- **A printed prediction missed:** recorded as a reading (item 311's (N-L) rule).
- No `PREDICTED RED:` line: every assertion is predicted green, and the zero half carries no band.

**Outcome (2026-10-13, 07:30 slot take-next, DONE — record, no glyph):** the `h38` flag, six cited literals and `{tag}` prints landed in `tests/validation/test_ans7_exact_strip.py` (+115 / −16); every window green first time, orphan checks 0 before and after. W0 (`-n 8`, 266 s, 15 passed) reproduced item 387 W2 to the printed digits (C_in(r) 1.5275048999871143e-11, +1.344645e-05, (AB-r′) 1.587 / 1.600 / 1.614e-5, (T2-r′) −2.301805e-06, Q(r) 254.469, (N2-r′) 494.07×, (N-L′) 2.17×, (N-r0′) bit-equal; (W) 1.107e-13). W0b (`-n 2`, 83 s, 15 passed) reproduced item 365b A in full digits (Im Y_in,A −4.061501062e-05 / +1.181776383e-07 / +4.088914063e-05 S, 2.271224e-6, (AB⁗) 1.593 / 1.605 / 1.619e-5, Q 606.6969050131794, +8.967920e-4). W1 (`-n 12`, 378 s, 7 passed, 140 270 cells; zero T35 / T37-tagged lines): C_in(r, h) 1.5288749690061608e-11 F, /C_in(0, h) − 1 = +1.358701e-05, /C_in(r, ×1) − 1 = +8.969327e-04 (lossless +8.967920e-04, difference +1.407e-07); (W) 1.068e-13; (F) / (F-Q0) green; (AB-rh′) 1.591 / 1.604 / 1.617e-05; Im Y_in,A −4.061198730e-05 / +1.200522737e-07 / +4.089005599e-05 S, one − → +, f₀,A/127.74 − 1 = −2.307317e-06 (0.143 × `T22P_ZERO_RESOLUTION`, printed); r_AB −2.307318e-06; Q(r, h) 252.30120548497052, /Q_0(h) 0.415860, closed form 243.7560 (Q_C 407.4659), excess +3.5056 %; **step reading Q(r, h)/Q(r, ×1) − 1 = −8.519904e-03** beside Q_0's −1.882737e-02, ratio 0.4525 (64 MHz 0.2915); R_in(r, h) 9698.113 Ω, /R_in(r, ×1) − 1 = −9.456959e-03; (N-Lh′) 2.18× bare / 1.394 Z_in form; (N2-rh′) 494.25×. Every prediction met except (N-r0h′) bit-equality (−2.220e-16 relative, inside the brief's 1e-15 fallback; as at 64 MHz). Logs `20261013T124951Z_ANS-7-item419-W0-esr128.log`, `…125429Z_ANS-7-item419-W0b-assemble128h.log`, `…125603Z_ANS-7-item419-W1-esr128-h001.log`.

**Spawn line:** "Use the implementer agent to execute `ANS-7` §9 item 419 per `docs/planning/chunks/ANS-7.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; run the orphan check first; never two heavy windows at once. Edit only `tests/validation/test_ans7_exact_strip.py`, every new line under the `h38` flag, and make no edit after W0 starts. `FEM_EM_ANS7_EXACT_ROUTE` and `FEM_EM_ANS7_EXACT_MHZ` stay unset. Never write `t127.74/tune_p2_h0.01.json`, any node or point file, any lossless tune file, or item 311's / 395's ESR files. Never open `docs/private/` or any `aed_results/`. Report every printed reading, met or NOT met."

### Review ruling 2026-10-13 03:00 (daily) on the item 419 draft

Queued as drafted. T38 cannot run on today's code: `_esr_ohm()` raises at 127.74 MHz on the h = 0.01 rung (`test_ans7_exact_strip.py:2209–2213`), and `tune_esr` would silently mislabel ("T35" tag, 64 MHz text, ×1 literals at `:2325–2326, :2361, :2399`). The item lands that fix under a new `h38` flag first, with new literals cited to their log lines, then re-gates item 387's W2 (`-n 8`) and item 365b's A (`-n 2`) before the T38 window. `-n 12` is the ceiling; nothing above it.

## §7 row moved verbatim 2026-10-14 (rotate_plan_archive.py chunks)

 The tuned copper F-small birdcage: `Z_in(f)`, f₀ and Q against AED | 🟡 `f₀` AGREE with its h-statement (T13, 2026-10-10 weekly); `R_in` / Q move away with h and are unconverged at h = 0.005 (T30, successive-change ratios 0.705 / 0.717 > ½, 2026-10-12 weekly), so the row holds 🟡 on the Q clause until F31 (`xl` entry 34, 10-15) and T41 (the 10-17 weekly); the ESR steps 11–13b (T26, T27, T34–T37) are records on the ×1 and h = 0.01 rungs; step 14 (T38, item 419, record) at 127.74 MHz on the h = 0.01 rung: C_in(r, h)/C_in(r, ×1) − 1 = +8.969e-4 beside the lossless +8.968e-4, Q(r, h)/Q(r, ×1) − 1 = −0.852 % beside Q_0's −1.883 % (ratio 0.4525; 64 MHz 0.2915), Q(r, h)/Q_0(h) 0.4159 (closed form 0.4018, excess +3.51 %), (AB-rh′) ≤ 1.62e-5, one − → + crossing at 0.143 × `T22P_ZERO_RESOLUTION` (scale only). T39 (item 420, record): Palace twin of the ESR-tuned node at matched order — Palace v0.18.1 order 2 on the same 80 798-cell hole mesh, gaps as parallel (R_p, C_p) from ESR ruling (iii)'s map, bare S11 within 2.569e-10 of T26's Route A at 64.000 MHz (lossless C path 2.565e-10; Palace's lumped R, C read as parallel). T40 (item 421, record): Palace at matched order reads f₀ within `T21P_ZERO_RESOLUTION` (+6.530e-11) and Q within T40's 5 % (−4.871e-9) of ours at the ESR-tuned node, ×1, degree 2 — Palace's own driven runs at 63.9 / 64.0 / 64.1 MHz (one exact (R_p, C_p) map per frequency) read through `parallel_rows` with the feed's (r + C) branch in post-processing; R_in −7.456e-10 (printed). The public external check §10 names for the Q clause — not Q's h-statement. T42 (item 436, record, provisional): Palace at matched order reads f₀ within `T21P_ZERO_RESOLUTION` (+6.073e-11) and Q within T40's 5 % (+5.444e-9) of ours at the 127.74 MHz ESR-tuned node, ×1, degree 2 — Palace's driven runs at 127.64 / 127.74 / 127.84 MHz (one exact (R_p, C_p) map per frequency at `T35_C_IN_R_X1`), R_in +2.906e-9, per-point S11 ties ≤ 2.373e-10, lossless C path 2.308e-10 (printed); the public external check at the second Larmor frequency. T43 (item 444, record, provisional): Palace at matched order on the h = 0.01 rung reads f₀ within `T21P_ZERO_RESOLUTION` (+6.437e-11) and Q within T40's 5 % (+1.829e-9) of T37's rows, R_in +6.507e-10, per-point S11 ties ≤ 2.591e-10, lossless C path 2.499e-10, and its own ×1 → 0.01 Q move −1.991402e-03 against ours' −1.991409e-03 (−0.199 %; equal to five significant digits, the last-digit gap the size of the two per-rung Q ties) (`20261014T094041Z_OPS-61-item444-W7-read.log`); a record for T41, no glyph. T44 (item 448, record, provisional): Palace at matched order on the h = 0.01 rung reads f₀ within `T21P_ZERO_RESOLUTION` (+5.912e-11) and Q within T40's 5 % (+1.027e-8) of T38's rows at 127.74 MHz, R_in +4.872e-09, per-point S11 ties ≤ 2.408e-10, lossless C path 2.093e-10, and its own ×1 → 0.01 Q move −8.519899e-03 against ours' −8.519904e-03 (−0.852 %; equal to six significant digits, |ΔQ-move| 4.8e-9) (`20261014T171034Z_OPS-61-item448-W6-read.log`); with T43 both Larmor frequencies; a record for T41, no glyph. **Done-when:** ✅ when T30 / F31 write Q's h-statement; the row holds 🟡 on the Q clause. *History: `docs/planning/chunks/ANS-7.md`.* 

### §9 item 460 — 2026-10-15 (§10 chain T45; lane 1, no solve)

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.*

**Licence.**
- §10 T45 (PROJECT_PLAN.md:4729–4732), weekly-enumerated, §9 rule (4): "`ans:7` `COMPARISON.md` gains a public 'Palace at matched order' section: T39–T44's `f₀` / Q / `R_in` ties and both codes' ×1 → 0.01 moves at 64 and 127.74 MHz (Palace figures public; the AED cells stay blank); `docrefs` 0 dead. Standard. Serial on nothing."
- The 2026-10-14 weekly ruling "T42 / T43 / T44 — ratified; what T41 reads" (PROJECT_PLAN.md:4620–4632) and the ratification list (:4661–4664).
- §5.4 Palace rules (i)–(iv) (PROJECT_PLAN.md:1173–1197): full Palace numbers publishable with version, config and mesh cited; matched order on a shared mesh; reproducible (the generator is tracked); **never a bridge to AED**.
- No band, tolerance, Done-when or target is set or moved. No §7 / §2 edit (the review's).

**What.** One new public section in `examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/COMPARISON.md`, rendered by the tracked writer (`07_tuned_copper_birdcage_driven_sweep.py`), never hand-written. Render-only, as item 400: no `ans:7` stage runs, `metrics.json` is not touched.
- **"AED cells stay blank" is met by having no AED column.** These are ESR-tuned rows (r = 0.2 Ω) and no AED ESR case exists (`COMPARISON.md:47`, `:64`); an empty AED column would imply a comparison that is not coming. The review may instead rule a blank column; the brief assumes none.
- **Writer design (closes a trap, below).** A module-level constant `PALACE_SECTION` (the section's lines, every figure a literal with a `# <log>:<line>` comment) and a new idempotent inserter `_write_comparison_palace()`: it removes an earlier copy of its own section (from `"\n" + its heading` to the closing "Not claimed" tail) and inserts the section directly before that tail. It is called (1) at the end of `_write_comparison`, after the existing `_write_comparison_esr128(m)` call (`07_…py:1238–1239`), and (2) as the last line of `_write_comparison_esr128` (`:1242–1274`). Reason for (2): `_write_comparison_esr128` cuts from its own heading to the tail (`:1268–1271`), so a stand-alone `XT128R` rewrite would otherwise delete a Palace section placed after it. Nothing else in the script changes.

**The section text (the slot transcribes; every figure is traced below).** Heading `## Palace at matched order (public; §10 T39–T44, \`OPS-61\`)`. Body:
- Scope paragraph: "Palace v0.18.1, order 2 — the same N1curl degree as ours, `ND (p = 2)` 570 332 on the ×1 hole mesh (80 798 elements) and 968 016 on the h = 0.01 hole mesh (140 270 elements), equal to our unknown counts — on the meshes our Route A solves, exported and configured by the tracked `scripts/probes/ops61_step1a_palace_ans4.py` (`export --hole-pec`, `export --hole-h001`; configs `ans7-configs`, `ans7-sweep-configs`, `ans7-128-configs`, `ans7-h001-configs [--f128]`; reads `read --ans7-esr`, `--ans7-sweep`, `--ans7-128`, `--ans7-h001 [--f128]`). The capacitor ESR r = 0.2 Ω enters Palace as its lumped R ‖ C — one exact (R_p, C_p) map per frequency (ESR ruling (iii); Palace's lumped R, C read as parallel, T39) — and the feed gap's (r + C) branch is added in post-processing on both sides. f₀, Q and R_in come from `parallel_rows` on the same three tune points (63.9 / 64.0 / 64.1 and 127.64 / 127.74 / 127.84 MHz) on both sides. No AED column: no AED ESR case exists."
- **Table A — ties, Palace / ours − 1 (signed; S₁₁ columns are |S_Pal − S_ours| / |S_ours|):**

  | Rung, frequency | f₀ | Q | R_in | S₁₁ per tune point (max of 3) | lossless-C S₁₁ | ESR-vs-lossless move |
  |---|---|---|---|---|---|---|
  | ×1, 64.000 MHz | +6.530e-11 | -4.871e-09 | -7.456e-10 | 2.580e-10 | 2.565e-10 | 1.994e-09 |
  | ×1, 127.74 MHz | +6.073e-11 | +5.444e-09 | +2.906e-09 | 2.373e-10 | 2.308e-10 | 5.648e-09 |
  | h = 0.01, 64.000 MHz | +6.437e-11 | +1.829e-09 | +6.507e-10 | 2.591e-10 | 2.499e-10 | 1.717e-09 |
  | h = 0.01, 127.74 MHz | +5.912e-11 | +1.027e-08 | +4.872e-09 | 2.408e-10 | 2.093e-10 | 4.356e-09 |

- **Table B — the values side by side (Palace / ours):**

  | Rung, frequency | Q | f₀/target − 1 | R_in [Ω] |
  |---|---|---|---|
  | ×1, 64.000 MHz | 164.5434 / 164.5434 | -2.583046e-06 / -2.583111e-06 | 3175.979 / 3175.979 |
  | ×1, 127.74 MHz | 254.4693 / 254.4693 | -2.301744e-06 / -2.301805e-06 | 9790.704 / 9790.704 |
  | h = 0.01, 64.000 MHz | 164.2158 / 164.2158 | -2.588508e-06 / -2.588573e-06 | 3166.748 / 3166.748 |
  | h = 0.01, 127.74 MHz | 252.3012 / 252.3012 | -2.307258e-06 / -2.307317e-06 | 9698.113 / 9698.113 |

- **Table C — each code's own ×1 → 0.01 move, X(h = 0.01) / X(×1) − 1 (Palace / ours):**

  | Frequency | f₀ | Q | R_in |
  |---|---|---|---|
  | 64.000 MHz | -5.462608e-09 / -5.461672e-09 | -1.991402e-03 / -1.991409e-03 | -2.906654e-03 / -2.906656e-03 |
  | 127.74 MHz | -5.514397e-09 / -5.512787e-09 | -8.519899e-03 / -8.519904e-03 | -9.456957e-03 / -9.456959e-03 |

- Reading paragraph (the weekly's, PROJECT_PLAN.md:4620–4626, quoted not extended): "The Q move is equal to five significant digits at 64 MHz and six at 127.74 MHz (−0.199 % / −0.852 %). The h-sequence of Q is therefore a property of the discretisation both codes share; no cross-solver reading shortens Q's h-statement, which is T41's, on ours alone."
- Not-claimed line: "Not claimed: an h-statement (two rungs, one move each); Q, f₀ or R_in accuracy against any external value; the h = 0.005 Palace twin (declined, 2026-10-14 weekly); any AED value."

**The anchor — every figure traced (the slot copies this table into the test; drafter-read 2026-10-15).** Logs are under `docs/testing/logs/` (rotation may gzip them: open `.log` or `.log.gz`, cite the base name). "Verbatim" = the cell string occurs on that line; "Rounded" = the float after the label, formatted as shown, equals the cell.

| Cell(s) | Log | Line | Check |
|---|---|---|---|
| ×1 64: f₀ +6.530e-11 | `20261013T141248Z_OPS-61-item421-W4-read.log` | 48 | verbatim |
| ×1 64: Q -4.871e-09 | same | 49 | verbatim |
| ×1 64: R_in -7.456e-10 | same | 50 | verbatim |
| ×1 64: S₁₁ 2.301e-10 / 2.569e-10 / 2.580e-10 → max 2.580e-10 | same | 38–40 | verbatim per line, max computed |
| ×1 64: lossless 2.565e-10; ESR move 1.994e-09; convention "parallel" | `20261013T140547Z_OPS-61-item420-W6-read.log` | 46; 48; 49 | verbatim |
| ×1 64: Q 164.5434 / 164.5434 (the `Q_Pal … vs Q_A …` pair) | 421 W4 | 49 | rounded `.4f` |
| ×1 64: f₀/target −2.583046e-06 (Palace) | 421 W4 | 48 | verbatim |
| ×1 64: f₀/target −2.583111e-06 (ours) | `examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/COMPARISON.md` | 54 | verbatim (metrics-rendered) |
| ×1 64: R_in 3175.979 / 3175.979 | 421 W4 | 50 | rounded `.3f` |
| ×1 127.74: f₀ +6.073e-11, Palace / ours f₀/target −2.301744e-06 / −2.301805e-06 | `20261013T215115Z_OPS-61-item436-W6-read.log` | 51 | verbatim |
| ×1 127.74: Q +5.444e-09; Q 254.4693 / 254.4693 | 436 W6 | 52 | verbatim; rounded `.4f` |
| ×1 127.74: R_in +2.906e-09; 9790.704 / 9790.704 | 436 W6 | 53 | verbatim; rounded `.3f` |
| ×1 127.74: S₁₁ 2.373e-10 / 2.306e-10 / 2.259e-10 → 2.373e-10 | 436 W6 | 39–41 | verbatim, max |
| ×1 127.74: lossless 2.308e-10; ESR move 5.648e-09 | 436 W6 | 42; 43 | verbatim |
| h 64: f₀ +6.437e-11, f₀/target −2.588508e-06 / −2.588573e-06 | `20261014T094041Z_OPS-61-item444-W7-read.log` | 51 | verbatim |
| h 64: Q +1.829e-09; 164.2158 / 164.2158 | 444 W7 | 52 | verbatim; rounded |
| h 64: R_in +6.507e-10; 3166.748 / 3166.748 | 444 W7 | 53 | verbatim; rounded |
| h 64: S₁₁ 2.591e-10 / 2.533e-10 / 2.507e-10 → 2.591e-10 | 444 W7 | 39–41 | verbatim, max |
| h 64: lossless 2.499e-10; ESR move 1.717e-09 | 444 W7 | 42; 43 | verbatim |
| 64 moves: f₀ −5.462608e-09 / −5.461672e-09; Q −1.991402e-03 / −1.991409e-03; R_in −2.906654e-03 / −2.906656e-03 | 444 W7 | 54 | verbatim |
| h 127.74: f₀ +5.912e-11, f₀/target −2.307258e-06 / −2.307317e-06 | `20261014T171034Z_OPS-61-item448-W6-read.log` | 51 | verbatim |
| h 127.74: Q +1.027e-08; 252.3012 / 252.3012 | 448 W6 | 52 | verbatim; rounded |
| h 127.74: R_in +4.872e-09; 9698.113 / 9698.113 | 448 W6 | 53 | verbatim; rounded |
| h 127.74: S₁₁ 2.408e-10 / 2.267e-10 / 2.128e-10 → 2.408e-10 | 448 W6 | 39–41 | verbatim, max |
| h 127.74: lossless 2.093e-10; ESR move 4.356e-09 | 448 W6 | 42; 43 | verbatim |
| 127.74 moves: f₀ −5.514397e-09 / −5.512787e-09; Q −8.519899e-03 / −8.519904e-03; R_in −9.456957e-03 / −9.456959e-03 | 448 W6 | 54 | verbatim |
| Palace v0.18.1 | `20261013T140240Z_OPS-61-item420-W3-palace-t39-esr.log`, `20261014T093349Z_OPS-61-item444-W3-palace-t43-f64.0.log` | 45 | verbatim `Git changeset ID: v0.18.1` |
| 80 798 / 140 270 elements | same two logs | 68 | verbatim `80798` / `140270` |
| ND (p = 2) 570 332 / 968 016 | same two logs | 104 | verbatim `570332` / `968016` |

Signs: the logs print `-7.456e-10` with an ASCII hyphen; the table may typeset the same digits. (X) compares the digits and the sign, not the glyph.

**Assertions (new `tests/unit/test_ans7_palace_section.py`, built on item 400's `tests/unit/test_ans7_comparison_writer.py`; tmp copies only).**
- **(X) Traced.** Every row of the table above, as a typed `(cell, log, line, mode)` tuple, holds against the tracked log (string containment or rounded-float equality). Each cell string also occurs inside the committed section exactly where `PALACE_SECTION` puts it.
- **(R) Byte-for-byte.** `_write_comparison(m)` on the committed `metrics.json` (`CASE_DIR` = tmp) reproduces the committed `COMPARISON.md`; the Palace heading occurs once, after the 127.74 MHz ESR section and before the "Not claimed" tail.
- **(I) Idempotence / trap closed.** Calling `_write_comparison_esr128(m)` on a tmp copy of the committed file leaves it byte-identical (the stand-alone `XT128R` path keeps the Palace section).
- **(B) No bridge.** The section contains no "AED" except inside "No AED column: no AED ESR case exists" and the not-claimed line's "any AED value", and no `docs/private` / `aed_results` / `COMPARISON_private` string.
- **(R0)** The repo's `COMPARISON.md` / `metrics.json` are unchanged by the test.

**Negative controls (rule (e); asserted, by construction).**
- (N-cite) One citation moved by one line (Q +1.027e-08 cited at `448 W6:53`) fails (X).
- (N-digit) A tmp copy of the committed file with "+6.530e-11" → "+6.540e-11" differs from the render.
- (N-old) The pinned pre-change writer (`3c359dd0`, the script's last commit at drafting; the slot re-reads `git log -1 --format=%h -- <script>`) renders a file with Palace-heading count 0 that differs from the committed file.

**Pre-registered change to item 400's test (the review licenses it or not).** `tests/unit/test_ans7_comparison_writer.py:96` asserts `nokey == old`: with the 127.74 key deleted, the new writer equals the pinned `b01dc790` writer byte for byte. With an unconditional Palace section that becomes false by construction. The licensed change: cut the Palace section (heading asserted present exactly once) from `nokey` before the comparison; nothing else in that test moves. This is item 417's ruling-(2) pattern. Without the licence the item stops before any edit.

**Executor, tier, cost.** Implementer (or the slot directly, as item 400). Standard tier; no MPI; complex build (the module imports dolfinx). Lane 1, class median 3 measured minutes; predicted ≈ 35 min (15 fixed + ≈ 15 writer and test + ≈ 5 windows). Depends on nothing.

**Windows.** Foreground, bare, repo-relative, Bash `timeout: 660000`, pytest with `-s`, no `-k`, no pipe.
1. **Pin (host, plain `git show`, not a harness call):** `git show 3c359dd0:examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py` into `logs/t45_pin/examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py` (mirrored depth: the module evaluates `parents[3]`, item 400's `IndexError: 3`). Load it with `__file__` set to the live path (item 416's fix).
2. **Edits** (script, new test, item 400's test line as licensed).
3. **W-regen** (gitignored helper `logs/t45_regen.py`, Write tool, calls `_write_comparison(json.loads(metrics))` on the committed `metrics.json`):
   `scripts/testing/run_and_log.sh ANS-7-item460-regen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 logs/t45_regen.py'"`
   Then host `git diff --stat`: `COMPARISON.md` changed, `metrics.json` not.
4. **W1:**
   `scripts/testing/run_and_log.sh ANS-7-item460-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans7_comparison_writer.py tests/unit/test_ans7_palace_section.py -s -v --tb=short'"`
5. **Pin module (rule (c); item 400's command, `20261012T170127Z_ANS-7-item400-pin.log`, 24 passed):**
   `scripts/testing/run_and_log.sh ANS-7-item460-pin "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/unit/test_doc_reference_exit_codes.py tests/unit/test_ans7_private_comparison.py -v -s --tb=short'"`
6. **Census:** `scripts/testing/run_and_log.sh ANS-7-item460-figures "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py'"` — `missing=0 broken=0`.
7. **Docrefs (gate `dead=0`, exit ≠ 1; exit 2 is the stale-artifact report severity, `20261015T005411Z_EX-114-post-docrefs.log:153`):**
   `scripts/testing/run_and_log.sh ANS-7-item460-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
8. **Commit** script, both tests, `COMPARISON.md` and the logs. The pre-commit private-leak hook (`scripts/testing/check_private_leak.py`) must print `clean:`; `SKIPPED` is not a pass — report it.

**Traps already paid for.**
- **The leak hook matches ≥ 8 significant digits** (`check_private_leak.py:76`); item 417's first commit was refused for an 8-digit figure. Every figure in the section is ≤ 7 significant digits (the table above already is); never paste a full-precision Q, f₀ in Hz or R_in.
- **Docrefs reads any `*.json` / `*.msh` / `*.csv` token as a run artifact** (`check_example_doc_references.py:116–122`) and flags it dead if absent beside the guide. Never write a Palace config name (`t39_esr.json`), mesh file name or `port-S.csv` into the section; cite the probe's subcommands instead (the probe `.py` exists).
- **`_write_comparison_esr128` cuts to the tail** (`07_…py:1268–1271`): hence design call (2) and assertion (I).
- **Never run an `ans:7` stage**: every stage rewrites `metrics.json`. Never hand-edit `COMPARISON.md`.
- **Pinned control:** a literal sha, never `HEAD:` (`OPS-43`); mirrored depth.
- **Rotation gzips logs:** open either suffix; the 421 / 420 logs are still `.log` at drafting.
- **No Palace-vs-AED sentence anywhere** (§5.4 (iv)); `aed_results/`, `COMPARISON_private.md` and `docs/private/` are never opened; the private writer (`:691–743`) is not called.
- Formatting: `≈`, never a tilde. No edit after W1 (rule (i)); never delete a log.

**Scope.** One public section on one case; arithmetic already printed by five read windows. **Not claimed:** an h-statement; Q's bracket (T41); h = 0.005 (declined); any AED comparison; the degree-1 or lossless `C_F` rows (those have no Palace twin); a §7 / §2 edit.

**Negative result.** Report and stop; park on `attempt/ANS-7-item460-<UTC>` and mark item 460 BLOCKED in the same commit (rule (d)).
- (X) red on a figure: the drafter's trace is wrong or a log moved. Print the line; never retype the figure to pass. The review re-traces.
- (R) red: a writer / regen mismatch; fix in-slot, re-run W-regen and W1.
- (I) red: the inserter pair is not idempotent; fix in-slot.
- Item 400's test red after the licensed cut: the cut is wrong, or the writer changed more than the Palace call; fix in-slot only if it is the cut.
- Leak-hook refusal: unstage and park; report the location only.
- No `PREDICTED RED:` line.

**Status it can move.** No glyph: `ANS-7`'s Done-when (§7 row, PROJECT_PLAN.md:3529) is Q's h-statement, which this does not write; rule (j) read, no mismatch. §10 T45 → DONE (the review marks it). The §7 row's T39–T44 clause may gain "(public table: `COMPARISON.md` § Palace at matched order)" — the review's edit. Known-issues: none.

**Outcome (2026-10-15 06:00 slot, take-next, implementer): DONE.** The writer gains `PALACE_HEAD`, `PALACE_SECTION` (each row's figures carry a log:line comment) and the idempotent `_write_comparison_palace()`, called at the end of `_write_comparison` and as the last line of `_write_comparison_esr128`; nothing else in the script moved. `COMPARISON.md` regenerated by the writer (`20261015T111816Z_ANS-7-item460-regen.log`; heading count 1, `metrics.json` unchanged), no AED column. New `tests/unit/test_ans7_palace_section.py`: (X) all 53 traced cells hold against their log lines and sit in the section; (R) render = committed, 6 569 B, heading once, after the 127.74 MHz ESR section and before the tail; (I) stand-alone `_write_comparison_esr128` leaves the file byte-identical; (B) no "AED" outside the two licensed phrases, no private path; (R0) repo files untouched; N-cite fails (X), N-digit differs, N-old (`3c359dd0`, 3 543 B, heading count 0) differs — all asserted. Item 400's test changed only by the licensed cut (heading asserted once in `nokey`, then cut; its pre-cut `(N-key)` print still reads `equals (N-old): False`, left as is — disclosed). Windows: first W1 `20261015T111824Z_ANS-7-item460-W1.log` 1 failed / 6 passed — the new test's expected slice omitted the inserter's trailing newline (test-only fix, no figure or writer line changed); W1 `20261015T111835Z_ANS-7-item460-W1.log` 7 passed; pin `20261015T111844Z_ANS-7-item460-pin.log` 24 passed; census `20261015T111902Z_ANS-7-item460-figures.log` `missing=0 broken=0`; docrefs `20261015T111903Z_ANS-7-item460-docrefs.log` `dead=0`, exit 2 (stale report severity, stale=119). §10 T45 → DONE is the review's mark.

**Spawn line:** "Use the implementer agent to execute `ANS-7` §9 item 460 per `docs/planning/chunks/ANS-7.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Never run an `ans:7` stage, never hand-edit `COMPARISON.md`, never write a figure with 8 or more significant digits, a `.json` / `.msh` / `.csv` name, or any AED value into the section; change item 400's test only by the licensed cut; never open `aed_results/`, `COMPARISON_private.md` or `docs/private/`."

### Review ruling 2026-10-15 03:00 (item 460)

Accepted. **Licensed:** the change to `tests/unit/test_ans7_comparison_writer.py` that cuts the new Palace section out of the rendered text before the (N-old) comparison at `:96` — the item-417 ruling-(2) pattern; the comparison still covers every other byte, and the change is disclosed. The section carries **no** AED column (no AED ESR case exists), not a blank one. Every figure ≤ 7 significant digits (the leak hook refused 8 on item 417). No band moves.

### §9 item 514 — 2026-10-18 (§10 T47, `ANS-7` verdict carried to `ans:7`; lane 1, no solve)

*Drafted by a read-only general-purpose agent for the 2026-10-18 03:00 daily review; the review reviews it. No AED figure was opened; every number below is ours or a relative level already public in `PROJECT_PLAN.md`.*

**Licence.**
- §10 chain T, T47 (`PROJECT_PLAN.md:4759–4763`), weekly-enumerated, §9 rule (4), quoted, not changed: "`ANS-7` ✅ carried to `ans:7`'s `COMPARISON.md` verdict section and `SPEC.md` — the h-statement row (the three rungs, the named steps, the printed rates, the Richardson limits, the ≈ 3 % level against AED First Order as INCONCLUSIVE pending T48) — and §2.1's wording checked; `docrefs` 0 dead. Standard. Serial on nothing."
- The 2026-10-17 weekly ruling T41 (`PROJECT_PLAN.md:4615–4643`) and the §7 row (`:3556`).
- §5.4 rules on AED (relative levels only, operator ruling 2026-09-24) and Palace (never set against AED, (iv)).
- **No band, tolerance, Done-when or target is set or moved.** No §7 row, §2 sentence, known-issues or dashboard edit — those are the review's (item 416's split). The slot *checks* §2.1 and reports; it does not edit it.

**Executor, tier, class, cost.**
- **Executor:** implementer (or the slot directly, as items 400 / 460).
- **Tier:** standard; no MPI (`python3` direct, `-n 1`); complex build (the script imports dolfinx).
- **Class:** lane 1 (median ≈ 3.5 measured min; ops/docs ≈ 3). **Cost:** ≈ 3.5 measured / ≈ 40 predicted slot-minutes — 15 fixed + ≈ 12 writer and new test + ≈ 6 for `SPEC.md` / guide text + ≈ 5 windows. Comparables: item 460 (3 measured / 35 predicted), item 416 (6 / 35).
- **Depends on:** nothing in §9. On disk: the gitignored pins `logs/ans7_writer_b01dc790/…/07_tuned_copper_birdcage_driven_sweep.py` and `logs/t45_pin/…/07_tuned_copper_birdcage_driven_sweep.py` (both present 2026-10-18 03:00; the two existing tests `pytest.fail` without them).
- **Family:** `ANS-7` chain T, attempt 1.

**Current state (verified 2026-10-18).**
- `examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/COMPARISON.md` (112 lines, 6 569 B) has **no verdict section**. Its header (`:3`) says "the weekly adjudicates at a stated relative level (T10b)", and its closing line (`:112`) says "no h-statement (T12′); no absolute f₀ / Q / R_in" — **stale since T13 (f₀) and T41 (Q / `R_in`)**. Its tables are ×1-mesh records; the h = 0.01 / 0.005 rungs appear nowhere in the example.
- The file is writer-rendered and never hand-edited: `_write_comparison` (`07_tuned_copper_birdcage_driven_sweep.py:1155–1240`) builds the header and tables, writes the closing line at `:1235–1236`, then calls `_write_comparison_esr128` (`:1243–1276`) and `_write_comparison_palace` (`:1351–1363`). **Both inserters assert the closing line as a literal** (`:1268`, `:1356`) and insert before it.
- `SPEC.md` § Status (`:23–58`) ends on the 2026-09-30 degree-1 adjudication (`:46–57`: "INCONCLUSIVE on the absolute rows … No tuning or match claim is licensed by this case yet"). It is hand-written; no test reads it.
- The guide `07_tuned_copper_birdcage_driven_sweep.md` § Not claimed (`:201–207`) carries the same stale verdict clauses: "no h-statement (T12′); no absolute f₀ / Q / R_in", "(T12′ / T13)", "the weekly adjudicates at a stated relative level". **This is a verdict line that must change.** No test pins the guide; docrefs and the figure census read it.
- Tests that pin `COMPARISON.md` byte for byte: `tests/unit/test_ans7_comparison_writer.py` ((R), and (N-key) `nokey == old` at `:96–104` after the licensed Palace cut) and `tests/unit/test_ans7_palace_section.py` ((R) / (I); `TAIL` literal at `:40`; **the (X) trace cites `COMPARISON.md` line 54 at `:59`**; ordering asserts at `:184–185`). `tests/unit/test_ans7_private_comparison.py` checks the tracked files stay byte-identical under the private writer.
- The script's last commit is `248416c9` (2026-10-15, item 460). The slot re-reads `git log -1 --format=%h -- <script>` and types the result as a literal.

**What. Edit only:**
1. `07_tuned_copper_birdcage_driven_sweep.py`:
   - a module-level `VERDICT_HEAD` / `VERDICT_SECTION` (a tuple of lines, as `PALACE_SECTION`), every figure a literal with a `# <log>:<line>` or `# PROJECT_PLAN.md:<line>` comment;
   - rendered **inside `_write_comparison`'s `lines` list, directly after the header paragraph (`:1178–1181`) and before `## Degree 2 …`**, so it sits above every table and the two before-tail inserters keep it untouched (they keep every pre-existing character, `:1246`);
   - under ruling R1 option (A), the closing line re-pointed through one module constant `CLOSING_LINE`, used at `:1235–1236`, `:1268` and `:1356`. Nothing else in the script moves: no stage logic, no table, not the header line `:3`, not `metrics.json`, not the private writer.
2. The new `tests/unit/test_ans7_verdict_section.py`.
3. The licensed edits to the two existing tests (below).
4. `SPEC.md` § Status: one new bullet after `:57`, and on the 09-30 bullet (`:57`) only the words "(re-read at degree 2 below)" appended after "yet". Nothing else.
5. The guide's § Not claimed (`:201–207`): the stale clauses only (text below).

**The verdict section (the slot transcribes; every figure is traced below).** Heading `## Verdict (public; adjudicated by the weekly review)`. Body:
- "**What was adjudicated:** our degree-2 Route B Cu rows on the exact SPEC rectangles, κ off, C = 66.29 pF on every gap — the first table below, on the ×1 mesh (h = 0.015 m, 80 798 cells), and the same rows on two refined rungs (h = 0.01 m, 140 270 cells; h = 0.005 m, 523 880 cells; `tests/validation/test_ans7_t30_rin_q_sequence.py`) — against AED First Order. The degree-1, C_in, ESR and Palace sections below are records, not the adjudicated pairing."
- "**Mode and copper contribution** (2026-09-30 weekly, degree 1): AGREE — one k = ±1 parallel resonance in the window in both codes (share ≥ 0.9999); the conductor's part of 1/Q, `Q_Cu/Q_PEC` and the Cu − PEC shift of f₀ within 5 %."
- "**f₀** (T13, 2026-10-10 weekly): AGREE with an h-statement — within 0.11 % of AED First Order at ×1 and 0.03 % at h = 0.005. Ours: 62.70108 / 62.72685 / 62.75178 MHz on h = 0.015 / 0.01 / 0.005 m; changes +25.76 / +24.94 kHz, ratio 0.968, printed rate p = 1.0625 (first order in h); Richardson limit ≈ 62.775 MHz (+0.04 % from h = 0.005)."
- "**Q and R_in** (T41, 2026-10-17 weekly): converge in h. Q 977.26 / 970.44 / 965.55 and R_in 18 441.75 / 18 313.60 / 18 223.23 Ω on h = 0.015 / 0.01 / 0.005 m; successive-change ratios 0.717 (Q) / 0.705 (R_in) on these steps, where a first-order quantity gives 1.00 and a second-order one 0.60; printed rates p = 1.6451 (Q) / 1.6786 (R_in); Richardson limits at those rates Q∞ ≈ 963.2 (−0.24 % from h = 0.005, −1.43 % from ×1), R_in∞ ≈ 18 182 Ω (−0.23 %). The T30 log's BRANCH lines (`20261011T094910Z_ANS-7-item364-T30.log:87`, `:91`) print the retired 'second change at most half the first' test; on these steps first order itself gives 1.00, so that test does not read the sequence and was retired (2026-10-17 weekly epitaph)."
- "**Against AED:** at h = 0.005 Q sits 3.1 % and R_in 2.7 % below AED First Order (T13), and the limits sit ≈ 3 % below it. The absolute Q / R_in is **INCONCLUSIVE — not an agreement and not a disagreement**: AED First Order is one discretisation whose own pass-to-pass trend is unread (T48, the operator's). `ANS-7` ✅ on its Done-when (Q's h-statement)."
- "Not claimed: any AED value; an absolute Q, R_in, `Z_in` or match; tuning accuracy at C_in or with ESR; human scale (`ANS-8`)."

**Deliberate departures from the plan's wording (each is a defect report below):**
- f₀'s rate is printed as the log's **1.0625**, not "1.063". The rate is 1.06245, so "1.063" is a double rounding.
- Q∞ is **963.2**, not "963.3". It is 963.2455 at the printed rate and 963.2458 at the solved rate.
- **No bracket** is printed: the letter does not list one, and the plan's bracket does not match the 10-11 form (D3).
- **No Palace sentence** appears in the verdict (D1).

**Closing line (ruling R1, option (A), the default).** The new `CLOSING_LINE` is: "Not claimed: any absolute `Z_in`, a match; an absolute Q or R_in (INCONCLUSIVE, T48); the tables above are one mesh each — the h-statements are the verdict section's." Option (B), if the review declines (A): keep the closing line byte-identical; then the palace test's `TAIL` and item 400's tail handling do not change.

**`SPEC.md` bullet (after `:57`):** "- [x] Re-adjudication at degree 2 on exact rectangles, κ off — **`ANS-7` ✅ (2026-10-17 weekly, T41).** f₀ AGREE with its h-statement (T13, 2026-10-10: within 0.11 % of AED First Order at ×1, 0.03 % at h = 0.005, first order in h). Q and `R_in` converge in h (p = 1.6451 / 1.6786 on h = 0.015 / 0.01 / 0.005 m) to a level ≈ 3 % below AED First Order: INCONCLUSIVE as an absolute Q, pending AED's own pass history (T48). Rungs, rates and limits: `COMPARISON.md` § Verdict."

**Guide § Not claimed (`:201–207`), replacement text:** "Any absolute `Z_in`, any match to AED; an absolute Q or R_in (INCONCLUSIVE against AED First Order, T48). C_in, f₀, R_in, Q at C_in are ×1-mesh degree-2 records — no absolute tuned C or tuning-accuracy claim. A measured ESR or its frequency dependence; ESR rows (64 and 127.74 MHz) are ×1-mesh degree-2 records. No AED value is read or written; the weekly's adjudication (f₀ AGREE with its h-statement, T13; Q / R_in converged in h, T41) is `COMPARISON.md` § Verdict."

**The anchor — every figure traced.** The slot copies this table into the new test as typed `(cell, source, line, mode)` tuples. `T30` = `docs/testing/logs/20261011T094910Z_ANS-7-item364-T30.log`, still `.log` at drafting; rotation may gzip it, so open either suffix and cite the base name.

| Cell(s) | Source | Line | Check |
|---|---|---|---|
| 80 798 / 140 270 / 523 880 cells; degree 2; kappa None; resolution None / 0.01 / 0.005 | T30 | 51 / 52 / 53 | verbatim `cells 80798` / `140270` / `523880` |
| h = 0.015 / 0.01 / 0.005 m (the named steps) | T30 | 86 | verbatim `on h (0.015, 0.01, 0.005)` |
| Q 977.26 / 970.44 / 965.55 | T30 | 89 | `.2f` of `v(x1)` / `v(0.01)` / `v(0.005)` |
| R_in 18 441.75 / 18 313.60 / 18 223.23 Ω | T30 | 85 | `.2f` (thousands space typeset) |
| f₀ 62.70108 / 62.72685 / 62.75178 MHz | T30 | 93 | `.5f` of value × 1e-6 |
| +25.76 / +24.94 kHz | T30 | 93 | `.2f` of Delta1 / Delta2 × 1e-3 |
| ratios 0.717 / 0.705 / 0.968 | T30 | 89 / 85 / 93 | `.3f` of `rho` |
| p = 1.6451 / 1.6786 / 1.0625 | T30 | 90 / 86 / 94 | verbatim |
| the BRANCH lines named retired | T30 | 87 / 91 | verbatim `BRANCH r_in_ohm` / `BRANCH q_loaded` |
| 1.00 / 0.60 (first / second order on the named steps) | computed | — | (h₂ − h₃)/(h₁ − h₂), (h₂² − h₃²)/(h₁² − h₂²) |
| Q∞ ≈ 963.2, −0.24 %, −1.43 %; R_in∞ ≈ 18 182 Ω, −0.23 %; f₀∞ ≈ 62.775 MHz, +0.04 % | computed (L) | — | see (L) |
| 0.11 % / 0.03 % (f₀, T13) | `PROJECT_PLAN.md` | 316–317 | verbatim |
| 3.1 % / 2.7 % (Q / R_in at h = 0.005, T13) | `PROJECT_PLAN.md` | 4632 | verbatim `Q −3.1 %` / `R_in` −2.7 % |
| ≈ 3 %; INCONCLUSIVE; T48 | `PROJECT_PLAN.md` | 4633–4637 | verbatim |
| share ≥ 0.9999; within 5 % | `SPEC.md` | 50 / 51–52 | verbatim |

**Assertions (new `tests/unit/test_ans7_verdict_section.py`, built on item 460's test; tmp copies only).**
- **(X) Traced.** Every row above holds against its source, and every cell string occurs inside the committed verdict section.
- **(L) The numpy smoke.** Parse `v(x1)`, `v(0.01)` and `v(0.005)` from `T30:85`, `:89` and `:93`. Solve (h₂ᵖ − h₃ᵖ)/(h₁ᵖ − h₂ᵖ) = ρ for p by bisection on [0.1, 5] and assert |p − printed| ≤ 5e-5 against `:86` / `:90` / `:94` (drafter-computed: 1.645099 / 1.678624 / 1.062450). Form v∞ = v₃ − (v₂ − v₃)·h₃ᵖ/(h₂ᵖ − h₃ᵖ) at the printed p, and assert the section's formatted strings equal the recomputation. Drafter values: Q∞ 963.2455, R_in∞ 18 182.19 Ω, f₀∞ 62.774678 MHz; −0.2382 %, −1.4342 %, −0.2252 %, +0.0365 %. Also assert the same limit to ≤ 1e-3 relative from the coarse pair (h₁, h₂) at the printed p: 963.2452 / 18 182.22 / 62.774664 MHz.
- **(R) Byte-for-byte.** `_write_comparison(m)` on the committed `metrics.json` (`CASE_DIR` = tmp) reproduces the committed `COMPARISON.md`. The verdict heading occurs once, after the header paragraph and before `## Degree 2`.
- **(I)** A stand-alone `_write_comparison_esr128` and a stand-alone `_write_comparison_palace` each leave the committed file byte-identical.
- **(B) No bridge, no leak.** The section contains:
  - no "Palace";
  - no `docs/private`, `aed_results` or `COMPARISON_private`;
  - no numeric token with ≥ 8 significant digits (digits counted across a typeset thousands space), with log filenames (`\d{8}T\d{6}Z_…\.log`) excluded before counting. The `ans:6` verdict carries one (`…/COMPARISON.md:37`, `20261012T070009Z_ANS-6-step5b.log:4595`) and passed the hook;
  - no "unconverged".
- **(R0)** The repo's `COMPARISON.md` and `metrics.json` are unchanged by the test.

**Negative controls (rule (e); asserted, by construction).**
- (N-p1) and (N-p2): Richardson at p = 1 / p = 2 gives Q∞ 960.65 / 963.91 (drafter-computed). Both formatted values ≠ "963.2" — the limit is rate-sensitive at the printed digit.
- (N-cite): p = 1.6451 cited at `T30:89` (one line early) fails (X).
- (N-digit): a tmp copy with "963.2" → "963.3" (the plan's rounding) differs from the render.
- (N-old): the pinned `248416c9` writer renders verdict-heading count 0 and differs from the committed file.

**Licensed changes to the two existing tests (the review licenses them or not; without the licence the item stops before any edit; item-417 ruling-(2) pattern, as item 460).**
- `test_ans7_palace_section.py:59`: the `CMP` line literal `54` re-pointed to the line the regenerated file puts "-2.583111e-06" on (54 + the verdict block's line count). (X) itself proves it.
- **Option (A) only:** `:40` `TAIL` becomes the new closing line.
- `test_ans7_comparison_writer.py:96–104`: cut the verdict section (asserted present once) as the Palace section is cut. Under (A), also restore the old closing line (asserted once) before `nokey == old`. Every other byte is still compared.

**Windows.** Foreground, bare, repo-relative, Bash `timeout: 660000`, pytest with `-s`, no `-k`, no pipe. **The harness call is the whole Bash command.**
1. **Pin (host, plain `git show`, not a harness call):** `248416c9:<script>` into `logs/t47_pin/examples/ansys_benchmarks/ans7_tuned_copper_birdcage_driven_sweep/07_tuned_copper_birdcage_driven_sweep.py`. Mirror the depth (`parents[3]`, item 400's `IndexError: 3`), and load with `__file__` set to the live path (item 416's fix).
2. **Edits** (script, new test, the two licensed lines, `SPEC.md`, guide).
3. **W-regen** (gitignored helper `logs/t47_regen.py`, Write tool; calls `_write_comparison(json.loads(<committed metrics.json>))`):
   `scripts/testing/run_and_log.sh ANS-7-item514-regen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 logs/t47_regen.py'"`
   Then host `git diff --stat`: `COMPARISON.md` changed, `metrics.json` not.
4. **W1:**
   `scripts/testing/run_and_log.sh ANS-7-item514-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 180 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans7_comparison_writer.py tests/unit/test_ans7_palace_section.py tests/unit/test_ans7_verdict_section.py -s -v --tb=short'"`
5. **Pin module** (rule (c); item 460's command, `20261015T111844Z_ANS-7-item460-pin.log`, 24 passed):
   `scripts/testing/run_and_log.sh ANS-7-item514-pin "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/unit/test_doc_reference_exit_codes.py tests/unit/test_ans7_private_comparison.py -v -s --tb=short'"`
6. **Census:** `scripts/testing/run_and_log.sh ANS-7-item514-figures "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py'"` — `missing=0 broken=0`.
7. **Docrefs** (gate `dead=0`, exit ≠ 1; exit 2 is the stale-artifact report severity — last read `RESULT: dead=0 guide=0 stale=133 stale_severity=report exit=2`, `20261017T172100Z_EX-128-post-docrefs.log:173`, command at `:12`):
   `scripts/testing/run_and_log.sh ANS-7-item514-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
8. **§2.1 check (no window, journal only).** Read `PROJECT_PLAN.md:309–322` and `:588–602` against the committed verdict and report match / mismatch per clause. The drafter's reading is D1 / D2 below. The review edits.
9. **Commit:** script, new test, the two licensed lines, `COMPARISON.md`, `SPEC.md`, guide, logs. The pre-commit leak hook (`scripts/testing/check_private_leak.py`, `MIN_SIGNIFICANT_DIGITS = 8` at `:76`) must print `clean:`; `SKIPPED` is not a pass — report it.

**Traps already paid for.**
- **The (X) line citation into `COMPARISON.md`** (`test_ans7_palace_section.py:59`, line 54): any insertion above it shifts it. Hence the licensed re-point. Never move the verdict below line 54 to dodge it; the inserters' cut-to-tail would then delete it (`:1269–1272`, `:1357–1360`).
- **Three closing-line literals** (`:1235–1236`, `:1268`, `:1356`) plus `TAIL` (`test_ans7_palace_section.py:40`). Under (A) they change together through `CLOSING_LINE` or not at all; one missed site turns `assert text.endswith(tail)` red on the next XT128R run.
- **T30's BRANCH lines read "unconverged"** (`T30:87`, `:91`): that is the retired "≤ ½" rule (epitaph 2026-10-17, `PROJECT_PLAN.md:4732–4737`). The verdict *names* them as retired, as the 03:00 audit's caveat asks (`ANS-7.md`, Audit 2026-10-18, first caveat), but never quotes their reading; (B) asserts the word "unconverged" is absent. The `.log` filename in that sentence is a `.log` token, not `.json` / `.msh` / `.csv`, so docrefs does not treat it as an artifact. W7 confirms this.
- **Leak hook: ≥ 8 significant digits.** Item 417's first commit was refused. f₀ goes as 62.70108 (7), never the table's 62.701082 (8), and never Hz.
- **Docrefs reads `*.json` / `*.msh` / `*.csv` tokens as artifacts** (`check_example_doc_references.py:116–122`). Cite the test module (`.py`), never a node JSON.
- **Never run an `ans:7` stage** (each rewrites `metrics.json`); never hand-edit `COMPARISON.md`; never call the private writer.
- **No Palace-vs-AED sentence** (§5.4 (iv)). See D1: the Palace reads are not of these rows.
- Pinned control: a literal sha, never `HEAD:` (`OPS-43`); mirrored depth. `≈`, never a tilde. No edit after W1 (rule (i)); never delete a log. `aed_results/`, `COMPARISON_private.md` and `docs/private/` are never opened.

**Scope.** One public verdict section, one re-pointed closing line under (A), one `SPEC.md` bullet, one guide paragraph, and a §2.1 check reported. **Not claimed:**
- any new arithmetic beyond T30's log and the weekly's levels;
- any AED value;
- a bracket;
- the `ANS-6` / F31 loss-partition reading;
- the T48 / T46 outcome;
- `examples/ansys_benchmarks/README.md` (R3);
- any §7, §2, known-issues or dashboard edit.

**Negative result.** Report and stop. Park on `attempt/ANS-7-item514-<UTC>` and mark item 514 BLOCKED in the same commit (rule (d)).
- (X) or (L) red: the drafter's trace or arithmetic is wrong, or the log moved. Print the line; never retype a figure to pass. The review re-traces.
- (R) or (I) red: a writer, inserter or regen mismatch; fix in-slot, then re-run W-regen and W1.
- An existing test red beyond the licensed lines: the writer changed more than licensed; fix in-slot only if it is the verdict or closing-line code.
- Leak-hook refusal: unstage and park; report the location only.
- No `PREDICTED RED:` line.

**Status it can move.** No glyph: `ANS-7` is ✅ with its Done-when met (`PROJECT_PLAN.md:3556`), so rule (j) does not bind. §10 T47 → DONE (the review marks it). The review may add "(public verdict: `COMPARISON.md` § Verdict)" to the §7 row. Known-issues 2026-10-17 (`docs/testing/known-issues.md:873`) stays open on T48.

**Defects found while drafting (for the review; the slot does not fix them).**
- **D1 — a Palace bridge on a configuration Palace never solved.** These sites set "Palace ties ours to six digits" beside the ≈ 3 % AED level:
  - the §7 row (`PROJECT_PLAN.md:3556`);
  - §2.2 (`:600`, "the level Palace reproduces to six digits");
  - the T41 ruling (`:4633–4635`, "so both FEM codes converge to that level");
  - known-issues 2026-10-17 (`:873`).

  T43 / T44 (and T39–T42) are the **ESR-tuned** reads at C_in(r) (Q ≈ 164 / 252; `COMPARISON.md:94–99`) on ×1 and h = 0.01 only. T41's rows are the **lossless-capacitor `C_F`** rows (Q ≈ 977), which item 460 recorded as having no Palace twin (`ANS-7.md:2825`), and Palace has no h = 0.005 rung (declined). So the shared-discretisation reading is an inference, not a measurement of this level, and the sentence pairs Palace with AED (§5.4 (iv)).
- **D2 — §2.1 `:313–314`** still reads "Degree 1, one rung; the absolute `f₀` / Q / `R_in` are §2.2's". The `f₀` clause is superseded by T13 in the next sentence (`:315–317`). Otherwise §2.1's T41 clause (`:317–322`) matches the verdict.
- **D3 — the bracket form.** The §7 row and `:4625` print `[960.7, 965.6]` / `[18 133, 18 223]` Ω "in the 10-11 ruling's form". That ruling (`ANS-7.md:1956`) set the form [x(0.005) − |Δ₂|, x(0.005) + |Δ₂|]. That gives [960.7, 970.4] / [18 133, 18 314] Ω; the printed bracket is its lower half.
- **D4 — two rounding slips** in `:4622–4624` and the §7 row: f₀'s p "1.063" (the rate is 1.06245) and Q∞ "963.3" (the limit is 963.2455).
- **D5 — stale lines outside T47's letter:**
  - `examples/ansys_benchmarks/README.md:110–115` still describes `ANS-7` v1 (15.57 pF on three gaps, 44–84 MHz, "Not ready for AED until the runnable half is ✅");
  - `examples/mri/16_study_copper_birdcage.md:87` and `examples/mri/24_study_tune_then_capacitor_sheets.md:91` say "`ANS-6` / `ANS-7` stay 🟡".

**Rulings needed.**
- **R1:** the closing line — (A) re-point it, the default, or (B) keep it byte-identical.
- **R2:** license the two existing-test edits above.
- **R3:** whether the README's `ANS-7` entry joins this item, or goes to a separate lane-1 line.

**§5.4 example check (for the review).** No example demonstrates `ANS-7`'s closed h-statement. The h = 0.01 / 0.005 rungs exist only in `tests/validation/test_ans7_t30_rin_q_sequence.py` and its log. `ans:7`'s `metrics.json` and `COMPARISON.md` carry the ×1 rows only; this item adds the rungs as traced text, not as an example run.

**Spawn line (if delegated):** "Use the implementer agent to execute `ANS-7` §9 item 514 per `docs/planning/chunks/ANS-7.md` (read the brief and any review ruling beneath it first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Never run an `ans:7` stage, never hand-edit `COMPARISON.md`, never call the private writer; write no figure with 8 or more significant digits, no `.json` / `.msh` / `.csv` name, no 'Palace' and no AED value into the verdict; change the two existing tests only as licensed; no §7 / §2 / known-issues edit; never open `aed_results/`, `COMPARISON_private.md` or `docs/private/`. You are the executor: do not spawn agents. On a red, park on `attempt/ANS-7-item514-*`, mark item 514 BLOCKED in the same commit, and report."

### Audit 2026-10-18 03:00 review — the `ANS-7` ✅ closure (the 2026-10-17 weekly, `169a81b4`): PASS

An `auditor` traced every digit of the row and the T41 ruling. The review re-read `20261011T094910Z_ANS-7-item364-T30.log:85–91`. Q 977.2613692833095 / 970.4389973340836 / 965.5456118046693, ρ 0.7172557529598782, p 1.6451 (`:89–90`); `R_in` ρ 0.7051723061233568, p 1.6786 (`:85–86`); identity (A) rel 0 (`:88`). The footer reads Status 0 (`:100–103`), and the log's commit `8e78c6b6` is an ancestor of `169a81b4`. F31's partition 0.4576 is at `20261015T070008Z_ANS-6-step5c.log:4894`. Its 0.4611 is at `20261012T070009Z_ANS-6-step5b.log:4583`, and its 0.4760 at `20261006T183450Z_ANS-6-item266-W1-p2-64-identity.log.gz:1206` (zcat line). The ruling cites only `:4894`, so the other two are recorded here. No AED value entered a tracked hunk.

**Ruled on the auditor's judgement call (check 5).** The executed asserts pin the inputs: the logged-digit reproduction at `NODE_S_RTOL` 1e-9 (module `:143`), the rung-swap control (`:163`) and the Δ₁ + Δ₂ identity (`:221`). The convergence classification is arithmetic on those inputs plus the weekly's ruling: ρ = 0.717 between the first-order 1.00 and the second-order 0.60 on these steps. §4 accepts a measured convergence rate, and the Done-when is "write the h-statement", so the review **accepts** it. No demotion.

**Caveats recorded, no status change:**
- The T30 log's `BRANCH` lines (`:87`, `:91`) still print the retired "> 0.5 → unconverged" rule. The 2026-10-17 epitaph retires it. A reader following the row to the log meets it, so T47 (§9 item 514) says so in `ans:7`'s verdict section.
- The row's Done-when trailed "the row holds 🟡 on the Q clause". It now reads "met 2026-10-17 (T41) …; until then the row held 🟡". This is a state-line edit only; the criterion is unchanged.
- `CLAUDE.md`'s onboarding paragraph still says "`ANS-7` stays 🟡 on its Q clause until F31 / T41". That file is the operator's, and the point is on the dashboard.

### Review ruling 2026-10-18 03:00 (item 514, T47)

Accepted, with these rulings. Spot-check: `tests/unit/test_ans7_palace_section.py:59` pins `COMPARISON.md` line 54, re-read. **The path is `tests/unit/`, not `tests/validation/`, wherever the brief says otherwise.**
- **R1:** the default. The stale closing line is re-pointed through one shared constant.
- **R2: licensed, disclosed.** Two existing-test edits: the line-54 citation in `tests/unit/test_ans7_palace_section.py` (plus `TAIL` under R1), and the verdict-section cut in `test_ans7_comparison_writer.py:96–104`. Neither may weaken an assertion. Each test is re-run green in the same slot (rule (c)), and the `git diff` of each is quoted in the journal.
- **R3:** the stale `ANS-7` v1 entry in `examples/ansys_benchmarks/README.md:110–115` joins this item: the same verdict, the same docrefs window. The `mri:16` / `mri:24` "stay 🟡" lines (D5) do not. They are guide pages with their own example gates, so they are named for a later lane-1 line.
- **D1 — named for the weekly, not edited here.** The drafter found that T43 / T44's six-digit Palace tie is on the ESR-tuned rows (Q ≈ 164 / 252), not the lossless-capacitor rows T41 reads (Q ≈ 977). That makes "Palace ties our Q", as used beside the ≈ 3 % AED level, an inference. It also sits close to §5.4's rule that no tracked file sets Palace against AED. The four sites (the §7 row, §2.2, the T41 ruling, known-issues 2026-10-17) are the weekly's text. The verdict this item writes **leaves Palace out**, as drafted.
- **D2–D4:** corrected in the verdict only (f₀'s rate as the log's 1.0625, Q∞ 963.25, no bracket printed). The §7 row's rounding is left for the weekly's next row pass.

**Outcome (item 514, 2026-10-18 implementer slot): DONE.** `VERDICT_SECTION` / `CLOSING_LINE` in the script; `COMPARISON.md` regenerated from the committed `metrics.json` (`20261018T093828Z_ANS-7-item514-regen.log`: heading count 1, closing line re-pointed, `metrics.json` unchanged); the verdict adds 9 lines, so the palace test's CMP line moved 54 → 63. W1 `20261018T094012Z_ANS-7-item514-W1.log` 10 passed in 3 s: (X) 18 traced + 10 computed cells hold; (L) p solved 1.678624 / 1.645099 / 1.062450 against printed 1.6786 / 1.6451 / 1.0625 (|Δp| 2.36e-5 / 1.48e-6 / 4.97e-5 — f₀'s sits just inside the 5e-5 bound because the log prints 4 decimals); limits Q∞ 963.2458, R_in∞ 18 182.18 Ω, f₀∞ 62.77469 MHz, coarse pair within 9e-8 relative; (N-p1/N-p2) 960.65 / 963.91 absent; (B) max 7 significant digits; (N-old) `248416c9` renders no verdict. Pin `20261018T094021Z_ANS-7-item514-pin.log` 24 passed; census `…-figures.log` missing=0 broken=0; docrefs `…-docrefs.log` dead=0 (stale=141, report severity, exit 2). **Executor departures:** (1) the brief's "Palace sections below are records" phrase contradicts its own (B) "no Palace" assert and the ruling's "leaves Palace out" — written "cross-solver sections"; (2) Q∞ printed 963.25 per the ruling (the brief's body said 963.2); (3) the plan-sourced (X) rows are re-located by a unique anchor string with the recorded offset, so a review's edit above them does not turn the test red. **§2.1 check (journal only):** `:313–314` "Degree 1, one rung; the absolute `f₀` / Q / `R_in` are §2.2's" still mismatches (D2); `:315–322` (T13 f₀, T41 Q / R_in, T48) match the verdict. §2.2 `:600` still carries "the level Palace reproduces to six digits" (D1, the weekly's). README `ANS-7` entry re-written (R3); `mri:16` / `mri:24` untouched.
