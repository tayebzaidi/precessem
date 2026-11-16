# ANS-8 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record; the SPEC is `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/SPEC.md`.

### §9 item 141 — 2026-10-01 (chain H, operator-commissioned `ANS-8` step 1 runnable half; queued by the 10:30 review)

`ANS-8` step 1, runnable half (`ans:8`): one shared rows function for both codes' 32×32 (implementer — every `ANS-7` runnable half ran under it; `example-runner` is EX-only — applying `example-runner`'s template rules: import, never restate; same-stem guide with the EX-15 headings; censuses predicted first; heavy, `-n 8` for stage S, `-n 1` for rows and unit tests, complex; ≈ 10 measured / ≈ 60 predicted min; depends on nothing). **If the slot reaches W1 after minute 45 of its hour,** commit the module, unit tests and W0 as a partial, mark the item BLOCKED "stage S pending", and stop — the next slot runs W1–W4 from the brief.

- **Shared function.** New pure-NumPy `src/fem_em_solver/ports/ring_modes.py`: `ring_mode_rows(s, *, port_map, freq_hz, c_f, feeds, drives, z0=50.0)` — S → Z with `circuit.s_to_z`, then SPEC rows 1–3: the k-blocks, the off-block norm, C from the antisymmetric eigenpair of k = 1 and k = 2, the Norton-terminated currents `I = Y_net V`, per-ring 16-bin k-fractions with the k = ±1 share and k = 0 leak, the feed admittances and the feed 2×2. `ring_port_map_from_names()` maps P(17+j) → (j, bottom), P(33+j) → (j, top). Our side uses it in stage R; the private writer uses it on AED's exported `"S"` (key layout `aed/ans8_hfss_pyaedt.py:619, 699`). Not in the case directory, because `run_examples.sh:129` lists every depth-2 `*.py` there as an `ans:` script.
- **Script.** `08_f_human_quadrature_mode_content.py` + guide, stages by `FEM_EM_ANS8_STAGE`; results in the gitignored `/workspace/logs/ans8/`. **Stage S** builds the matrix as `test_wf10_step1c_c_sweep` W1 does (`test_wf10_f_human_tuned_circuit.py:452–495`, from imported pieces: `_build_f_human`, `_ring_quadrature_slots`, `_solve_one_drive`, `terminal_form_deficit`, `_respec`) and saves an npz. **Stage R** loads that npz, never recomputes (missing ⇒ stop, the `:448` precedent), writes `metrics.json` (added to `LISTED_EXAMPLE_ARTIFACTS`) and `COMPARISON.md` with the AED column blank. The private writer checks the knob first, then whether the JSON exists.
- **Anchors (asserted, on our matrix, at `c_f = C_MODE_F` imported from `:310`)** — item 133's circuit-only records, cited to `20261001T140319Z_WF-10-step1d-W1.log`: (A1) C_mode, relative ≤ `LEG_D0_REPRODUCTION_BAND` (`:10652`); (A2) C₂ = 7.545440170601e-11, same band (`:10665`); (A3) off-block 4.182887e-03, relative ≤ 1e-6 (`:10652`); (A4) **cw** (1, +j) two-port 0.790993 / 0.952414, ≤ 1e-6 absolute (`:10666`), \|k\|-folds both rings ≤ 1e-4 (`:10667–68`); (A5) P17-alone 0.789857 / 0.951115 and its folds (`:10660, :10663–64`); (A6) opposite pair 0.988772 / 0.987500 and four-port 0.989921 / 0.988704, even-\|k\| 0.0000 (`:10669–74`). (G1) cells within `CELL_COUNT_BAND`, reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ `MATRIX_PASSIVITY_CEILING`; κ printed beside 4.590991e-03 (`:10644`). (G2) the fixture's measured slot map equals the SPEC naming map. (G3) if item 128's held `STEP1C_NPZ` exists: same `order`, entrywise ≤ `LEG_D0_REPRODUCTION_BAND` (item 128 A3 3.804e-13); absent ⇒ print and skip.
- **Synthetic (asserted ≤ 1e-12), `tests/unit/test_ring_mode_rows.py`.** (U1) a block-circulant Z under a shuffled port order is recovered, off-block ≈ 0. (U2) the C_mode closed form. (U3) on a D16 network with the SPEC loads: four-port k = 0 leak and every even-k bin ≈ 0; opposite-pair even-k ≈ 0; ccw and cw \|k\|-fold shares equal (mirror symmetry). (U4) every non-±1 mode open (j·1e9 Ω): four-port share ≥ 1 − 1e-9, leak ≤ 1e-12. (U5) ideal-source two-port: the dense feed-ring share equals the k-domain closed form from the per-k 2×2 solves. (U6) sense: drive (ii)'s weights equal `quadrature_phase_weights([0, 1], "ccw")`; ccw sits wholly in NumPy bin 15, cw in bin 1.
- **Private writer test** (`tests/unit/test_ans8_private_comparison.py`, `ANS-7` item 68's four sentinel anchors plus the leak control); its synthetic AED S reaches the table through `ring_mode_rows`.
- **Negative control.** Asserted: rows at `c_f` = C₂ under P17-alone give feed-ring share 0.099704 (`:10665`, measured in item 133), asserted < 0.5 against 0.789857 at C_mode. Asserted (synthetic): a shuffled slot map gives off-block > 0.1 (U1). Predicted, printed: \|ccw − cw\| on our matrix ≈ 2e-3 (SPEC).
- **Windows** (bare, repo-relative, foreground, `timeout: 660000`; every window sets `FEM_EM_ANS_SKIP_PRIVATE=1`).
  - (pre) `ANS-8-census-docrefs-before` / `ANS-8-census-setupfig-before`, the exact commands of `20260929T093115Z` / `…093120Z` line 12.
  - **W0** (≈ 30 s): `scripts/testing/run_and_log.sh ANS-8-step1-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 180 python3 -m pytest tests/environment tests/unit/test_ring_mode_rows.py tests/unit/test_ans8_private_comparison.py -s -v --tb=short'"`
  - **W1** stage S: `FEM_EM_ANS8_STAGE=S timeout -k 30 590 mpiexec -n 8 python3 examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py`, in §5.1's durable-capture idiom verbatim in the shape of `20260929T094009Z_ANS-7-item69-W2-B2.log:12`, trailing `; exit $rc` included. Predicted 260–330 s, ≈ 12.1 GiB (item 133's 32×32 done at t = 246 s, `:10645`; item 128 W1 423 s, 12.17 GiB with extra solves).
  - **W2** stage R, `-n 1`, `timeout -k 30 120` (≈ 30 s).
  - **W3** the pin test (`…094733Z…` line 12's command), after `git add metrics.json`.
  - **W4** `scripts/testing/run_and_log.sh ANS-8-list "./scripts/run_examples.sh --list"`, predicting the `ans:8 -> …/08_…py` line.
  - (post) docrefs: runnable +1, `exit != 1`; setup figures `missing` +1 (`ans:8`; the `EX-57` fallback draws it, the `ANS-7` precedent).
- **Traps.** Anchors run at `C_MODE_F` and **cw**; published rows at the SPEC's 1.1311e-10 F and **ccw** — mixing them shifts results ≈ 1e-5 and turns A4 red. Item 133 used Thevenin feeds (Z_th), the SPEC Norton: shares agree, V/I does not. NumPy `fft` bins as defined. Do not import, edit or run `aed/` (git-excluded). **Take no code from `attempt/WF-10-step1d-…`** — H8b stays BLOCKED, and landing its helpers would partly land a blocked item. Stage S calls are collective; rank 0 writes the npz, then Barrier. Heavy imports lazy, so the unit test can import the script. The private test must `monkeypatch.delenv` the skip knob. The pin test reads `git ls-files`. `-s`; no pipes; after a kill, `pgrep -c python3` first; no edit after the last window; never delete a log.
- **Scope.** Not: AED numbers or any verdict; sweep row 4; Route A; the loss partition (left blank and said so); copper; the H8b identity band. "One function for both codes" holds because our private writer recomputes AED's rows from AED's exported S — the AED script cannot import `fem_em_solver` (its `__init__` imports dolfinx).
- **Status it can move.** `ANS-8` ⬜ → 🟡 "runnable half landed", on the Done-when written by the 10:30 review; the SPEC's "Runnable half" box ticked.
- **Negative result.** A1–A6 red with G3 green: the shared algebra differs from item 133's — a divergence finding: park, BLOCKED, no re-record. G3 red: the rebuild differs from item 128's matrix — known-issues entry; stop before rows. G2 red: an operator finding (AED port naming) — stop. A U-red: park. W1 killed at 590 s: not a result — park with the measured time; the review re-prices (S/R split). No `PREDICTED RED:` line: every assertion is predicted green, so any red stops the slot.
- **Outcome (2026-10-01, 13:30 CDT slot): DONE, every assertion green.** Pre-census `20261001T183730Z_ANS-8-census-docrefs-before-clean.log` (67 runnable, exit 0) and `…183732Z_ANS-8-census-setupfig-before.log` (67 ok, 0 missing); a first `…183717Z_ANS-8-census-docrefs-before.log` ran with the new script already in the tree (68, guide violation, exit 1) and is superseded by the clean one (script moved aside for it). W0 `…183808Z_ANS-8-step1-W0.log` 22 passed, 48 s (U1 off-block 4.5e-15, shuffled-map control 0.823; U2 ≤ 2.2e-16; U3 even-k ≤ 4.3e-30, mirror 1.6e-15; U4 share 1 − 0, leak 4.7e-22; U5 identical). W1 `…183904Z_ANS-8-step1-W1-S.log` 289.7 s at `-n 8`, 11.08 GiB summed (predicted 260–330 s, ≈ 12.1 GiB): cells 507266, κ 4.590991e-03, reciprocity 8.409e-15, σ_max 0.999923038. W2 `…184403Z_ANS-8-step1-W2-R.log`: G2 0 mismatches; G3 8.976e-13; A1 4.381e-13, A2 2.021e-14, A3 2.699e-08 relative; A4–A6 shares ≤ 4.49e-07 absolute, folds ≤ 4.98e-05, even-\|k\| ≤ 9.5e-07; negative control 0.099704; \|ccw − cw\| 1.912e-03 (predicted ≈ 2e-3). Published (SPEC C, ccw): feed ring 0.789066, leak 0.173623, other ring 0.950224. W3 `…184419Z_ANS-8-step1-W3-pin-staged.log` 20 passed (pinned 75). W4 `…184433Z_ANS-8-list.log` lists `ans:8`. Post-census `…184434Z` 68 runnable, exit 0; `…184435Z` setup figures 68 / ok 67 / missing 1 (`ans:8`, exit 2 — the `ANS-7` precedent). The private writer's AED JSON layout (a block carrying `"S"` plus `port_names`, or `{"re", "im"}` / `[re, im]` entries) was written without opening `aed/` (slot rule); the operator's first private run checks it against `aed/ans8_hfss_pyaedt.py:619, 699`.

### §9 item 161 — 2026-10-01 (ANS-8 SPEC v2 rows, operator-commissioned; queued by the 18:00 review)

`ANS-8` step 1 runnable half, **SPEC v2 rows 2b / 2c / 4a** — pure algebra on the held 32×32, no new solve (SPEC.md:223–226). Executor: implementer, because the item extends `src/` (`example-runner` handles EX chunks only; item 141 is the precedent). `example-runner`'s template rules still apply: import, never restate; same-stem guide; censuses predicted first. Standard tier at `-n 1`; heavy `-n 8` only if the held npz is absent. ≈ 7 measured / ≈ 40 predicted min (+6 if stage S must regenerate). Depends on nothing (item 141 is DONE). **45-minute rule:** if the slot reaches W-S after minute 45, it commits the module, unit tests and W0 as a partial, marks the item BLOCKED "stage S pending", and stops.

**Two premises corrected by the drafter (read before coding).**
- (1) SPEC C 113.11 pF is *not* our C_mode to ≈ 4e-13. The 4.381e-13 compares our matrix with the `C_MODE_F` record. SPEC C is C_mode rounded to 5 figures, so it sits **+1.700e-5** relative above it. 2b on our matrix therefore does **not** reproduce row 2. The identity is **2b = item 141's A4 / A5 / A6-opposite**, which were computed at `C_MODE_F`. Row 2 − 2b is our own detuning record: printed, not asserted.
- (2) Under exact C16 symmetry the controls' k = 0 leak and other forbidden bins are exactly 0, but their k = ±1 share is 1 only when every other mode is open (U4's construction). The ≈ 0.99 is real k = ±3, ±5, ±7 content. The gauge sees C16 breaking only, not mirror breaking.

**Held matrix:** `/workspace/logs/ans8/s32.npz` (host `logs/ans8/s32.npz`, gitignored). Stage R re-runs without stage S. On a fresh tree it stops at `08_…py:238–240` and the slot runs stage S.

- **Shared function.** Add `ring_mode_rows_v2(s, *, port_map, freq_hz, c_spec_f, feeds, drives, controls, z0=50.0)` to `src/fem_em_solver/ports/ring_modes.py`. It composes `ring_mode_rows` (:126) and does not edit it, so the v1 output and U1–U6 stay bit-identical.
  - It returns `{"fixed_c": rows at c_spec_f, "own_c_mode": rows at c_f = fixed_c["c_mode_f"], "controls": {"fixed_c": {label: rows}, "own_c_mode": {label: rows}}, "c_own_f", "f0_one_point_hz": freq_hz·√(c_own/c_spec)}`.
  - **Own C_mode is wired from the function's own row-1 output and is never a parameter.** That way the private writer cannot hand AED our C_mode.
  - Every control is its own `ring_mode_rows` call with its own `feeds`. The Norton 1/z0 sits only on driven ports (:117–121): the opposite pair is P17 and P25 alone, with P21 a plain C gap, matching item 141's A6 calls (`08_…py:204–207`).
  - The feeds' capacitor is already C (:117), as SPEC 2b requires.
  - Add `"anti"` for k = 1 beside 2b, so a wrong eigenpair pick shows (ours 1.000000, `20261001T140319Z_WF-10-step1d-W1.log:10648`).
- **Script** (`08_f_human_quadrature_mode_content.py`):
  - `SPEC_CONTROLS = {"opposite": (("P17","P25"), (1, −1)), "four_ccw": (("P17","P21","P25","P29"), (1, −1j, −1, +1j))}`, mapped through `port_map` as `SPEC_FEED_NAMES` is (:336).
  - Stage R adds the B, C and F assertions below **after** A1–A6. It writes `metrics.json` with v1's `rows` unchanged, a new `rows_v2` key, and `"spec": "SPEC.md v2 (2026-10-01)"`.
  - `COMPARISON.md` gains rows (2b), (3b), (2c) at both C, and (4a), with the AED column blank. Row 2 is captioned "detuning record (SPEC v2)".
  - `aed_rows` (:446–454) returns `ring_mode_rows_v2` per AED design, and `_private_comparison_text` (:457) gains the v2 columns. **Control leaks go in the gap table as absolute differences**, because the relative formula (:475) divides by zero on a synthetic D16 AED.
  - **P7 (clobber guard):** `_write_private_if_available` leaves alone an existing `COMPARISON_private.md` whose first line is not the writer's own header (:466). It reports "kept" and never reads past line 1. That path runs only with the knob unset, so never in this slot.
  - Update the guide §1, §3 (new step 5) and "Not claimed", and tick the SPEC Status line 28.
- **Anchors on our matrix.** Asserted unless labelled otherwise, at own C_mode through `ring_mode_rows_v2`. R = `20261001T184403Z_ANS-8-step1-W2-R.log`, H = `20261001T140319Z_WF-10-step1d-W1.log`.
  - **B0.** `own_c_mode.c_f` is bitwise equal to `fixed_c.c_mode_f`, and within `LEG_D0_REPRODUCTION_BAND` (1e-9) of `C_MODE_F` (`test_wf10_f_human_tuned_circuit.py:310`; R:41 measured 4.381e-13).
  - **B1.** 2b P17-alone shares 0.789857 / 0.951115 to ≤ 1e-6; \|k\| folds to ≤ 1e-4 (A5 literals `08_…py:76–78`; R:47–49).
  - **B2.** 2b cw shares 0.790993 / 0.952414, folds likewise (A4, `:73–75`; R:44–46).
  - **B3.** The 2b P17-alone and cw blocks equal `anchor_rows`' `C_MODE_F` blocks to ≤ 1e-9 absolute on every share, leak and fraction. The expected difference is ≈ 4e-13: the two C values differ by 4.381e-13 relative, and the measured sensitivity is ≈ 0.88 share per unit relative C.
  - **B4.** Feed-ring \|2b ccw − 2b cw\| = 1.912e-3 ± 1e-6 (R:55), with ccw < cw.
  - **B5 (predicted, printed).** 2b ccw feed-ring share ≈ 0.789081, range [0.789075, 0.789087]. Row 2 − 2b: feed ring in [−2e-5, −1e-5], other ring in [0, +5e-6], for all three drives.
  - **C1.** At own C, the opposite pair's shares are 0.988772 / 0.987500 to ≤ 1e-6 (A6, R:50). Its k = 0 leak and every even \|k\| fold are ≤ 5e-5 (item 141's bar; R:52 measured 9.46e-7).
  - **C2.** At own C, the four-port ccw's leak and even \|k\| folds are ≤ 5e-5. The C4 selection rule does not depend on sense; cw measured 3.38e-7 at R:53. Predicted, printed: share within 5e-3 of A6-cw's 0.989921 / 0.988704, and wrong-sense bins {1, 5, 9, 13} ≈ 0.
  - **C3.** Both controls at SPEC C have leak ≤ 5e-5. Predicted, printed: shares within 1e-4 of the own-C values.
  - **C4.** At each C, the control shares exceed every two-port share. Predicted, printed: opposite / 2b-ccw ≈ 1.253, against a ceiling of 1/0.789081 = 1.267.
  - **F1.** `f0_one_point_hz` = `64e6·√(C_MODE_F/SPEC_C_F)` computed in-test to ≤ 1e-9 relative, and f₀ < 64 MHz. Printed: ≈ 63 999 455.96 Hz (−544.04 Hz).
  - **F2 (printed).** The one-point form's first-order error is (1 − p)·δ/2, with p = d ln Im λ₁/d ln ω unmeasured and δ = +8.50e-6.
  - **Unchanged and re-run:** G1–G3 and A1–A6 (`08_…py:250–321`).
  - **D1 (host `git diff`, not a window).** `metrics.json` has insertions only plus the one `spec` line, and v1 lines of `COMPARISON.md` are unchanged — on the held npz only. After a regeneration, D1 becomes "v1 rows within 1e-9 relative".
- **Synthetic tests** (asserted ≤ 1e-12; appended to `tests/unit/test_ring_mode_rows.py`; U1–U6 untouched).
  - **U7 (2b).** On `_d16_blocks(c_mode=c0)`, c0 = 1.0731e-10:
    - `own_c_mode.c_f` = c0, and is bitwise equal to `fixed_c.c_mode_f`;
    - the own-C blocks are bitwise equal to a direct `ring_mode_rows(c_f=…)` call;
    - they are bitwise invariant when `c_spec_f` is 1.1311e-10 or 0.5e-10;
    - U5's ideal-source dense result equals the k-domain closed form at c0.
  - **U8 (2c selection rules, C16 not D16).** On U1's mirror-broken C16 matrix, at both C:
    - opposite pair: every even bin ≈ 0;
    - four-port ccw: every bin outside {3, 7, 11, 15} ≈ 0;
    - four-port cw: every bin outside {1, 5, 9, 13} ≈ 0.
    - On U4's open-others network, both controls' share ≥ 1 − 1e-9.
  - **U9 (gauge negative control).** A single-gap 50 % self-impedance defect on non-feed slot 2 gives a control leak > 1e-9 (predicted ≥ 1e-4, printed).
  - **U10 (4a).** With λ₁ = R + jωL, R > 0: f₀ = 1/(2π√(LC)), and rebuilding at f₀ gives `c_mode_f` = C. Control, asserted as a closed form: λ₁ = R + j(ωL − 1/(ωC_s)), C_s = 10·C, true root 5 % below 64 MHz — f₀ equals its own one-point closed form and misses the true root by > 1e-3 relative.
- **Private writer tests** (`tests/unit/test_ans8_private_comparison.py`; the existing four keep passing):
  - **P5.** The v2 sentinels land in the AED row and none appears in tracked files. The AED 2b `c_f` is bitwise equal to the AED's own row-1 C_mode (109.87654 pF synthetic, :112), not ours.
  - **P6 (cross-wiring control).** AED 2b recomputed at our synthetic C_mode differs from its sentinel by > 1e-6.
  - **P7 (clobber guard).** A pre-existing `COMPARISON_private.md` with a foreign header is left byte-identical.
  - The fixture keeps `monkeypatch.delenv` on the knob (:101). The slot never opens the real private file or `aed_results/`.
- **Negative controls (asserted).**
  - The C₂ control through the v2 function (`c_spec_f = C2_RECORD_F`): `fixed_c` P17-alone feed-ring share 0.099704 ± 1e-6 (R:54; H:10665), with `own_c_mode` bitwise equal to the SPEC-C call. Separation 7.92×, against an arithmetic ceiling of 10.03×.
  - Gauge sensitivity on our matrix: the 2b two-port ccw leak is > 0.1 (0.1734 at C_mode, R:45), against controls ≤ 5e-5. Predicted, printed: separation ≈ 1.8e5.
  - Also U9, U10's control, and P6.
- **Windows** (bare, repo-relative, foreground, `timeout: 660000`; every window sets `FEM_EM_ANS_SKIP_PRIVATE=1`; a slot that spawns the implementer puts the §9 "Every item" harness sentence verbatim in the spawn prompt).
  - **(pre, before any edit to the script or guide)** `ANS-8-v2-census-docrefs-before` / `…-setupfig-before`, using the exact commands at `20261001T183730Z_ANS-8-census-docrefs-before-clean.log:12` / `20261001T183732Z_…setupfig-before.log:12`. Predict from item 141's post-censuses (`…184434Z…:38–40`, `…184435Z…:82,104–107`); never hard-code a count.
  - **(npz check, host)** `ls -la logs/ans8/s32.npz logs/wf10/h8a-held-s32.npz`. Both were present at 18:00.
  - **W-S (only if `s32.npz` is absent).** `20261001T183904Z_ANS-8-step1-W1-S.log:12` verbatim, chunk `ANS-8-v2-WS`, raw log `ans8-v2-ws-raw.log`, trailing `; exit $rc` kept. Predicted 260–330 s and ≈ 11.1 GiB (measured 289.7 s, 11.083 GiB). G3 must then be green against `h8a-held-s32.npz`.
  - **W0 (≈ 50 s, standard):** `scripts/testing/run_and_log.sh ANS-8-v2-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 180 python3 -m pytest tests/environment tests/unit/test_ring_mode_rows.py tests/unit/test_ans8_private_comparison.py -s -v --tb=short'"` — predicted ≈ 29 passed.
  - **W1, stage R (≈ 2 s):** R:12's command with chunk `ANS-8-v2-W1-R` (`timeout -k 30 120`, `mpiexec -n 1`).
  - **W2, pin test:** `20261001T184419Z_ANS-8-step1-W3-pin-staged.log:12`'s command (20 passed, 8 s).
  - **(post)** Both censuses again. Predicted unchanged (N → N).
- **Traps.**
  - **SPEC C ≠ our C_mode** (+1.700e-5), so 2b ≠ row 2.
  - **Sense:** A6's four-port is **cw** (`08_…py:207`); 2c's is **ccw**. Only the opposite pair reproduces a record.
  - **f₀ ratio** is √(C_own/C_spec); inverted, it gives +544 Hz.
  - **Lossless synthetics** at own C_mode are near-singular. Keep R > 0.
  - **Relative-gap division** by an exactly-0 synthetic leak.
  - **Pre-census contamination:** census before editing (item 141's first pre-census ran with the script already in the tree).
  - **The operator's private file** on disk is host-owned (likely written outside our writer). Un-knobbed writes replace it whole (:490), so the knob is set in every window and P7 guards it.
  - Keep v1's `rows` and every A-literal. Do not import, edit or run `aed/`. Leave `_aed_designs` (:422) unchanged.
  - Lazy heavy imports, `-s`, no pipes, no edit after the last window (rule (i)), `pgrep -c python3` after a kill, never delete a log.
- **Scope.** Our rows and the shared function only. Not covered: any AED number or verdict, the row-4 sweep interpolation, the degree-2 32×32 (a new solve the weekly prices), copper (`WF-12` step 3), `P_coil`, Route A, measuring p, adjudication, ✅.
- **Status it can move.** On an all-green run, the §7 `ANS-8` row's 🟡 runnable-half sentence gains "v2 rows (2b / 2c / 4a) landed 2026-10-01 (§9 item 161)", which replaces the clause "rows 2b / 2c / 4a are a pure-algebra follow-up". The SPEC Status line 28 is ticked for v2. The glyph stays 🟡.
- **Negative result.** Every outcome is recorded with its log.
  - U7–U10 or P5–P7 red: the composition is wrong. Park on `attempt/ANS-8-v2-<UTC>`, BLOCKED in the same commit (rule (d)), and regenerate no metrics.
  - A1–A6 or G1–G3 red: v1 regression. Stop before writing rows.
  - B0–B4 red with A green: 2b is wired to the wrong C. Park.
  - C1–C3 leak red on our matrix: a finding. Report it, BLOCKED; the review rules on the bar, which is never widened in-slot.
  - F1 red: wiring error. Park.
  - W-S killed: not a result; apply the 45-minute rule.
  - No `PREDICTED RED:` line: any red stops the slot.
- **Outcome (2026-10-02, 00:00 CDT slot): DONE, every assertion green.** Held `s32.npz` present, so no W-S. Pre-census `20261002T050048Z_ANS-8-v2-census-docrefs-before.log` (68 runnable, exit 0) / `…050050Z_…setupfig-before.log` (68 / ok 67 / missing 1 `ans:8`, exit 2), run before any script or guide edit; post-census `…050750Z` / `…050751Z` identical, as predicted. W0 `…050559Z_ANS-8-v2-W0.log` 29 passed, 58 s: U7 own C relative 2.2e-16, bitwise, c_spec-invariant, ideal-source dense = closed form; U8 forbidden bins ≤ 3.4e-29, open-others shares 1; U9 clean leak 8.4e-30, defect 5.43e-2; U10 f₀ 2.2e-16, rebuilt C 1.55e-15, control misses the true root by 4.84e-3; P5 AED 2b c_f bitwise the AED's own row 1, P6 \|d\| 7.23e-4, P7 foreign file kept byte-identical. W1 `…050705Z_ANS-8-v2-W1-R.log` 2 s: G1–G3, A1–A6 unchanged from item 141; B0 own C 1.131080769863e-10 F bitwise row 1, 4.381e-13 from `C_MODE_F`, SPEC C +1.7002e-5 above, anti(k = 1) 1.000000; B1 0.789857 / 0.951115 (\|d\| 4.41e-7), B2 0.790993 / 0.952414 (4.49e-7), folds ≤ 4.98e-5; B3 4.486e-13; B4 1.912287e-3, ccw < cw; B5 ccw 0.789081, row 2 − 2b feed ring −1.49e-5, other ring +2.4e-6 for all three drives (in range); C1 0.988772 / 0.987500, leak 9.455e-7; C2 leak 9.307e-7 (printed shares 0.987648 / 0.986317, wrong-sense bins 1.2e-6); C3 leak 9.456e-7, shares within 5.3e-7 of own C; C4 controls ≥ 0.986317 > two-port ≤ 0.952417 at both C, ratio 1.2531 vs ceiling 1.2673; F1 63 999 455.96 Hz, 2.189e-13 relative, −544.04 Hz. Negative controls: C₂ 0.099704, own C bitwise the SPEC-C call, separation 7.92×; gauge two-port ccw leak 0.1736 vs controls ≤ 9.5e-7, separation 1.836e5. F2 note: the log's printed δ is SPEC C / C_own − 1 = 1.700e-5, the brief's 8.50e-6 is the f₀ offset (half of it) — a label difference on a printed-only line, script not edited after its last window. W2 pin `…050738Z_ANS-8-v2-W2-pin.log` 20 passed (pinned 75). D1 (host `git diff`): `metrics.json` +2745 / −1 (the `spec` line), `COMPARISON.md` +33 / −0.


## §7 row moved verbatim 2026-10-01 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/ANS-8.md`)* **Human scale: the tuned F-human birdcage, two-port 90° quadrature — is the non-k = 1 content seen by AED too?** (operator directive 2026-10-01, interactive session, pulling the reserved human-scale case forward; was "commissioned by the weekly that reads H9"). **Step 1 — mode content, before H9:** F-human (16 legs, 32 ring gaps, the `WF-7` phantom) at 64 MHz with **C on all 32 gaps** at our H8b `C_mode` (fixed, not tuned on the AED side) and the **two feeds 90° apart** (the clinical drive). Coil **σ = 800 S/m solved inside**, exactly our fixture (as `ANS-4` did), so the comparison is like-for-like now; a copper column follows `WF-12`. **Primary route: 32-port S export + NumPy termination** (`ANS-7` Route B's pattern — ports 3…32 in C, the feeds in 50 Ω), which yields the gap currents without a field calculator; an in-model route optional. **Rows (both codes, same algebra):** the feed-ring and opposite-ring k-spectrum of the gap currents under the quadrature drive (k = ±1 share, k = 0 leak), the k = 1 block's `C_mode`, and the resonance of the terminated network near 64 MHz; Zero Order (= our degree 1) with First Order optional. **Readout:** agreement on the share / leak at a stated relative level ⇒ the ≈ 0.79 / 0.17 is the drive's physics; disagreement ⇒ our fixture or port model. **Step 2 (after H9):** tuned `\|B₁⁺\|` and the 10 g SAR hotspot, as originally reserved. AED values private (§5.4); relative levels public. Solve resources 12 cores / 1 TB (operator 2026-09-29). | 🟡 **Runnable half landed 2026-10-01 (13:30 CDT slot, §9 item 141).** `ans:8` (`ans8_…/08_f_human_quadrature_mode_content.py` + same-stem guide) and the shared pure-NumPy `fem_em_solver.ports.ring_modes.ring_mode_rows`. Stage S (`-n 8`, 289.7 s, 11.08 GiB summed RSS, `20261001T183904Z_ANS-8-step1-W1-S.log`): 507266 cells (5.200e-03 vs record), κ 4.590991e-03, reciprocity 8.409e-15, σ_max 0.999923038. Stage R (`-n 1`, `20261001T184403Z_ANS-8-step1-W2-R.log`): G2 slot map = SPEC naming (0 mismatches); G3 vs item 128's held matrix entrywise 8.976e-13; (A1) C_mode relative 4.381e-13, (A2) C₂ 2.021e-14, (A3) off-block 2.699e-08 relative; (A4) cw 0.790993 / 0.952414, (A5) P17-alone 0.789857 / 0.951115, (A6) opposite pair 0.988772 / 0.987500, four-port 0.989921 / 0.988704 — all ≤ 4.5e-07 absolute, folds ≤ 4.98e-05; control at C₂ 0.099704 < 0.5; \|ccw − cw\| 1.912e-03 (predicted ≈ 2e-3). Published rows (SPEC C, ccw): feed-ring share 0.789066, k = 0 leak 0.173623, other ring 0.950224. Unit U1–U6 + private writer 22 passed (W0 `20261001T183808Z`). `metrics.json` / `COMPARISON.md` written, AED column blank. ✅ only after the AED half lands and the weekly adjudicates. *Was:* ⬜ Commissioned (operator). **SPEC v1 written 2026-10-01** (interactive session acting for the weekly at the operator's request): `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/SPEC.md` — Route B (32 lumped ports + NumPy termination) primary, C = 113.11 pF on all 32 gaps, feeds P17 / P21, drives P17-alone / ccw / cw, rows: `C_mode` from the k = 1 block, per-ring k-spectrum (k = ±1 share, k = 0 leak), feed admittances; σ = 800 solve-inside, copper optional. **Next:** the runnable half (`ans:8`, one shared rows function for both codes' 32×32) — queued 2026-10-01 10:30 review as §9 item 141; the operator's AED replication in parallel. **AED half, v1 column (Zero Order, σ = 800) landed 2026-10-01** (operator; private, not yet adjudicated). **SPEC v2 2026-10-01** (operator directive, interactive session): row 2b — every drive at each code's *own* `C_mode` — becomes the primary readout, because v1's fixed-C row confounds mode content with the codes' `C_mode` gap; symmetric-drive controls 2c; row 4 redefined as the k = 1 block resonance (v1's `Im Y_in` zero has no root on this fixture); `P_coil` by power balance; copper and First Order AED columns commissioned (sweep Zero Order only). Row 2b is post hoc for the v1 column and blind for the three v2 columns. Item 141 landed against v1's rows; rows 2b / 2c / 4a are a pure-algebra follow-up on `ring_mode_rows`. **Step 1 runnable half, v2 rows, Done-when (written 2026-10-01 18:00 review, rule (j); §9 item 161):** `fem_em_solver.ports.ring_modes.ring_mode_rows_v2` is the one function that stage R applies to our held 32×32 and that the private writer applies to each AED design's 32×32. It adds rows 2b (rows 2–3 at the matrix's own row-1 `C_mode`, wired from its own output), 2c (the opposite pair and four-port ccw at both C) and 4a (f₀ = f·√(C_mode/C)) with no new solve. On our matrix it (a) reproduces item 141's A5, A4 and A6 opposite-pair shares at own `C_mode` to 1e-6 and \|ccw − cw\| 1.912e-3 to 1e-6, holds the controls' k = 0 leak ≤ 5e-5 at both C, and returns the C₂ control 0.099704. It (b) passes the synthetic C16 selection-rule and λ₁ ∝ ω closed forms ≤ 1e-12. It (c) shows each AED design's 2b at that design's own `C_mode` in the sentinel unit test, with `FEM_EM_ANS_SKIP_PRIVATE` honoured. `metrics.json` / `COMPARISON.md` are regenerated with the AED column blank and v1's `rows` unchanged. Then the runnable-half sentence reads "v2 rows landed", and the glyph stays 🟡. Note: SPEC C 113.11 pF sits +1.700e-5 above our `C_mode`, so 2b ≠ row 2 on our side; the difference is our own detuning record. **Step 1 runnable half, Done-when (written 2026-10-01 10:30 review, rule (j)):** `ans:8` (`ans8_…/08_f_human_quadrature_mode_content.py` + same-stem guide) rebuilds our κ-corrected F-human 32×32 at 64 MHz, degree 1, from imported construction, and runs SPEC rows 1–3 through the single `fem_em_solver.ports.ring_modes.ring_mode_rows`, which (a) reproduces item 133's records on our matrix — C_mode and C₂ to `LEG_D0_REPRODUCTION_BAND`, the off-block norm, and the P17-alone / two-port cw / opposite-pair / four-port circuit shares at C_mode to their printed digits; (b) passes the synthetic C16 / D16 closed forms ≤ 1e-12; (c) is the function the private writer applies to the AED 32×32 (sentinel unit test, `FEM_EM_ANS_SKIP_PRIVATE` honoured); `metrics.json` and `COMPARISON.md` written with the AED column blank ⇒ ⬜ → 🟡 "runnable half landed"; ✅ only after the AED half lands and the weekly adjudicates share and leak at a stated relative level. 

### §9 item 366 — 2026-10-11 (§10 H15c, `ANS-8` step 1c; lane-1 revalidation)

**Outcome 2026-10-11 07:30 slot: ✅ DONE, first try.** Module `tests/validation/test_ans8_degree2_rows.py`; W0 `20261011T123321Z_ANS-8-step1c-W0.log` (13 passed, 1 skipped, 74.7 s; synthetic D16: 18 classes, spread 0), W1 `20261011T123443Z_ANS-8-step1c-W1.log` (14 passed, 72.2 s; standard, `-n 1`). Anchors: (F) / (G2′) green on both npz (0 slot-map mismatches); (G1′) reciprocity 9.581484e-14, σ_max 0.999913975; (L-d1) σ_max 2.596e-10, class means 4.279e-10 (the az synthesis holds, worst azimuth residual 0); (L-d2) σ_max vs npz key 0.000e+00, vs record 4.297e-10, class means 2.473e-10; NEG P17 × 1.01 3.147e-04 = 314.7× on both (d2 inside the predicted [2e-4, 4e-4]). Rows (d2 κ off | d1 κ off | d1 κ on): `C_mode` 1.082287732e-10 | 1.131081186e-10 | 1.131080770e-10 F; `C_k2` 7.175798070e-11 | 7.545442055e-11 | 7.545440171e-11 F; off-block 3.909947e-03 | 4.182881e-03 | 4.182887e-03; anti 1.000000 all; f₀ one-point vs SPEC C 62 603 820.95 | 63 999 467.72 | 63 999 455.96 Hz; 2b ccw feed-ring share 0.713591 | 0.754987 | 0.789081, leak 0.236516 | 0.201230 | 0.173606, other ring 0.928063 | 0.939635 | 0.950222; 2c controls' leak at own C ≤ 4.549e-07 (d2); C₂ control (d2, P17 alone at `C_k2`) 0.111506 < 0.5. Order move at matched κ: `C_mode` −4.3139e-02, `C_k2` −4.8989e-02, f₀ −2.1807e-02 (ratio 0.5055); 2b ccw Δshare −0.041396, Δleak +0.035286. κ move at degree 1: `C_mode` −3.7e-07, but 2b ccw share +0.034094, leak −0.027624. **Predictions:** `C_mode` move in [0.5 %, 6 %] — hit; anti ≥ 0.99 — hit; 2b ccw share / leak within 0.03 of degree 1's κ-on 0.789081 / 0.1736 — **missed** (−0.0755 / +0.0629; printed predictions, readings, not reds): the move is order (−0.041) plus κ (−0.034), each about the size of the predicted window. Records only; no AED value read.

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain H, "H15c (lane 1, no solve). `ANS-8` step 1c — `C_mode` / `f₀` / k-share from the saved degree-2 32×32" (PROJECT_PLAN.md:4660–4668), and the 2026-10-10 weekly's entry-29 ruling (:4543–4545): the saved degree-2 32×32 has not been through `ring_mode_rows_v2`; H15c computes the rows without a solve; the 10-12 weekly reads rows 1 / 4 privately. **This item prints our numbers only.** No band, tolerance, Done-when or target is set or moved; every band is imported by name.

**`ANS-8` step 1c: SPEC v2 rows (1, 2b, 2c, 4a) on the saved degree-2 σ = 800 32×32, beside degree 1, with the loaded matrix tied to its window.**
- **Executor:** implementer (item 141 / 161 precedent; `example-runner` is EX-only). **Tier:** standard, `-n 1`, complex build (pure NumPy, but the item-161 W0 environment is kept). **Class:** lane-1 revalidation (median 17). **Cost:** ≈ M measured / ≈ 30 predicted min (15 fixed + ≈ 12 code + W0 ≈ 1 + W1 ≈ 1 + spare).
- **Depends on nothing open.** Needs two gitignored files, both on the host when this brief was drafted (`ls -la logs/wf7/`, 2026-10-11 ≈ 03:15 CDT): `logs/wf7/wf7-s32-degree2-n16.npz` (**20 823 bytes**, written 2026-10-09 02:13 by `xl` entry 29, `20261009T070009Z_WF-7-step0e-H15b.log:12318`) and `logs/wf7/wf7-s32-degree1-n8.npz` (20 823 bytes, H15a, `20261005T003402Z_WF-7-H15a-r3-W2.log:10560`). Container path `/workspace/logs/wf7/`. **If the degree-2 npz is absent, stop and mark item 366 BLOCKED ("npz absent: entry 29 is re-posed by the weekly"); never re-solve.** If only the degree-1 npz is absent, run degree 2 alone and say so.
- **Never open** `aed_results/`, `COMPARISON_private.md` or `docs/private/`.

**What the two files are (read from the code and the logs; the review checks).** Both were saved by `tests/validation/test_wf7_f_human_ring_matrix.py:917–931`: raw S at 50 Ω (`kappa_applied=False`, **κ off**), the **boxed** F-human σ = 800 solid mesh (507 266 cells, `WF7_BOXED_CELL_RECORD` `:270`; `BOXED_NPZ_BASENAMES` `:274`), keys `s32`, `order`, `slot_k`, `slot_sign`, `cells`, `degree`, `freq_hz`, `sigma_max`, `reciprocity`, `kappa`, `pooled_kappa`, `ranks`. **The published degree-1 rows are not this convention:** `ans:8` stage S applied κ (4.590991e-03) to the sheets (`08_…py:216–222`), so `20261002T050705Z_ANS-8-v2-W1-R.log` and the tracked `metrics.json` `rows_v2` are **degree 1, κ on**. The order move is therefore read at matched κ — degree 2 κ off against the degree-1 κ-off npz — and the κ-on rows are printed as a third column, so κ and order are never mixed in one difference.

**Code.** One new module, `tests/validation/test_ans8_degree2_rows.py`; **no `src/` change, no edit to `ans:8`, `metrics.json` or `COMPARISON.md`** (H15c is a record, not a SPEC row).
- Knob `FEM_EM_ANS8_STEP1C=1` (the `_enabled` / `_need` pattern, `tests/validation/test_mat15_tuned_copper_heating.py:112–123`): unset ⇒ the gated test skips; `"1"` ⇒ an absent npz **fails** (never skips); anything else ⇒ `ValueError`.
- Load the `ans:8` script through `importlib.util.spec_from_file_location`, as `tests/unit/test_ans8_private_comparison.py:41–50` does; use its `rows_v2`, `_slot_port_map`, `SPEC_C_F`, `SPEC_FEED_NAMES`, `SPEC_CONTROLS`, `CONTROL_LEAK_BAR`, `NEGATIVE_CONTROL_CEILING` and `fem_em_solver.ports.ring_modes.ring_port_map_from_names`. Never restate a constant.
- From `tests/validation/test_wf7_f_human_ring_matrix.py`: `_class_means` (`:331`), `SIGMA_MAX_D1_RECORD_BOXED` (`:214`), `CLASS_MEAN_D1_RECORD_BOXED` (`:215`), `WF7_BOXED_CELL_RECORD` (`:270`), `SAVE_CONTROL_FACTOR` (`:239`). Bands: `REUSE_REPRODUCTION_RTOL` (`tests/validation/test_port_birdcage_ring_matrix.py:171`, = 1e-6), `MATRIX_PASSIVITY_CEILING` (`:150`), `AZIMUTH_STEP_DEG` (`:185`), `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py:86`), `_reciprocity_ratio` (`tests/validation/test_port_birdcage_ring_column.py:776`).
- **`_class_means` needs per-port azimuth and ring.** The npz has no azimuths; synthesise `az = slot_k · AZIMUTH_STEP_DEG`, `zc = slot_sign`. Only differences enter `_azimuth_class` (`test_port_birdcage_ring_matrix.py:206–215`), and the windows' own residual was 5.7e-14° (`…H15b.log:12344`), so the classes are the slot differences. The degree-1 anchor below validates this synthesis against an independent record before degree 2 is read.
- **Record constants (comparands, not bands; 361-3 precedent):** `SIGMA_MAX_D2_RECORD_BOXED = 0.999913975` (`20261009T070009Z_WF-7-step0e-H15b.log:12320`, `:12343`) and `CLASS_MEAN_D2_RECORD_BOXED`, the 18 means of `:12321–12338` at their printed 10 digits, each with its line number, in the `CLASS_MEAN_D1_RECORD_BOXED` dict shape.

**Anchors (asserted, every band imported).**
- **(F) file identity, both npz:** `degree` 2 / 1; `kappa_applied` False; `cells == WF7_BOXED_CELL_RECORD`; `freq_hz` 64e6; 32 ports, 16 slots per ring.
- **(G2′) slot map:** `_slot_port_map(npz)` equals `ring_port_map_from_names(order)` — 0 mismatches (stage R's G2, `08_…py:346–355`).
- **(G1′) degree 2:** reciprocity ≤ `RECIPROCITY_BAND`, σ_max ≤ `MATRIX_PASSIVITY_CEILING`. The window read 9.581e-14 and 0.999913975 (`…H15b.log:12342–12343`).
- **(L-d1) degree-1 reproduction — the synthesis check:** σ_max of `wf7-s32-degree1-n8.npz` vs `SIGMA_MAX_D1_RECORD_BOXED`, and all 18 class means vs `CLASS_MEAN_D1_RECORD_BOXED`, each ≤ `REUSE_REPRODUCTION_RTOL`. The same comparison on the same file read 2.596e-10 / 4.279e-10 (`20261005T003402Z_WF-7-H15a-r3-W2.log:10562`).
- **(L-d2) the H15c anchor:** σ_max recomputed (`np.linalg.norm(s, 2)`) vs the npz's own `sigma_max` key and vs `SIGMA_MAX_D2_RECORD_BOXED`, and all 18 class means vs `CLASS_MEAN_D2_RECORD_BOXED`, each ≤ `REUSE_REPRODUCTION_RTOL`. Predicted ≤ 1e-9 (10-digit print rounding ≤ 5e-10; the window computed the same numbers from the same file, `:12320`).

**Negative controls (rule (e)).**
- **Asserted (backed):** on the degree-1 array, P17 × 1.01 moves the worst class-mean reproduction to ≥ `SAVE_CONTROL_FACTOR` × `REUSE_REPRODUCTION_RTOL`. Same comparison, same file: 3.147e-04 = 314.7× (`20261005T003402Z_WF-7-H15a-r3-W2.log:10581`; ceiling ≈ 0.01/32 ≈ 3.1e-4).
- **Degree 2 — the miss asserted against the same floor, the factor predicted:** P17 × 1.01 on the degree-2 array against `CLASS_MEAN_D2_RECORD_BOXED`. Predicted [2e-4, 4e-4], printed. The window's own control (`…H15b.log:12339`, 1.092e-01) was against the degree-1 record, so it is not the same comparison and backs nothing here.
- **Printed only:** the C₂ control through `rows_v2(…, c_spec_f=<degree-2 C_k2>)` — P17-alone feed-ring share predicted < `NEGATIVE_CONTROL_CEILING` (degree 1 κ on: 0.099704, R:79); the 2c controls' k = 0 leak at own C, predicted ≤ `CONTROL_LEAK_BAR` (degree 1: 9.46e-7, R:70). Neither is asserted on degree 2; a miss is a reading.

**Printed (the record), three columns — degree 2 κ off | degree 1 κ off | degree 1 κ on (from the tracked `metrics.json` `rows_v2`, loaded, never retyped):**
- row 1: `c_mode_f` (= `c_own_f`), `c_k2_f`, `off_block`, `anti_k1`;
- row 2b at each matrix's **own** C_mode: P17-alone / ccw / cw feed-ring (bottom) share and leak, other-ring (top) share, and the |k| folds;
- row 2c: the opposite pair and four-port ccw at own C (share and leak per ring);
- row 4a: `f0_one_point_hz` against SPEC C (= 64 MHz·√(C_mode / 113.11 pF));
- the order move at matched κ: (d2 − d1)/d1 for C_mode, C_k2 and f₀, and d2 − d1 for each 2b share and leak; the κ move at degree 1 beside it.

**Predictions (written first, printed).** C_mode order move \|Δ\|/C in [0.5 %, 6 %] (the class means moved 0.4–10.9 % d1 → d2, `…H15b.log:12321–12338`), so f₀ moves by about half that, sign printed. The 2b ccw feed-ring share stays within 0.03 of degree 1's 0.789081 (R:65), and the leak within 0.03 of 0.1736. The anti fraction stays ≥ 0.99.

**Windows.** Foreground, bare, repo-relative; Bash `timeout: 660000`; `-s`; `tests/environment` first; `FEM_EM_ANS_SKIP_PRIVATE=1` set (nothing here writes a private file, but the knob stays set in every `ANS-8` window).
- **W0** (smoke: the gated test skips; an ungated `test_ans8_step1c_coding` asserts the knob parser through `monkeypatch`, plus `_class_means` on a synthetic D16-symmetric 32×32 with synthesised az returns 18 classes with zero spread):
  `scripts/testing/run_and_log.sh ANS-8-step1c-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 120 python3 -m pytest tests/environment tests/validation/test_ans8_degree2_rows.py -s -v --tb=short'"`
- **W1** (standard, `-n 1`; predicted ≤ 60 s: stage R's whole v2 algebra took 2 s, `20261002T050705Z_ANS-8-v2-W1-R.log:81`, and here it runs on two matrices plus the C₂ call):
  `scripts/testing/run_and_log.sh ANS-8-step1c-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS8_STEP1C=1 timeout -k 30 180 python3 -m pytest tests/environment tests/validation/test_ans8_degree2_rows.py -s -v --tb=short'"`
- No durable capture: both windows are well under ≈ 10 min.

**Traps already paid for.**
- **κ convention:** the npz files are κ off; the published v2 rows are κ on (stage S `_respec(base_specs, kappa)`). Never difference d2 against the κ-on rows and call it the order move.
- **SPEC C ≠ any matrix's C_mode** (+1.700e-5 on degree 1, item 161 trap 1). The 2b rows are at each matrix's own C_mode, wired inside `ring_mode_rows_v2`; never pass one matrix's C_mode to another.
- **Sense:** ccw = (1, −1j) (`SPEC_DRIVES`); A6's four-port is cw, 2c's is ccw (item 161).
- **Boxed files are never overwritten** (`_save_knob` refuses both basenames, `:318–322`). Open them read-only (`np.load`, `allow_pickle=False`).
- **Do not import or run `aed/`** (git-excluded). Lazy heavy imports in the script are what make the `importlib` load cheap; do not add module-level dolfinx imports to the new module.
- `-s`; no `-k`; no pipe; never delete a log; no module edit after W0 starts (rule (i)).

**Scope.** Rows from two saved σ = 800 boxed-mesh matrices, degree 2 and 1, κ off, 64 MHz; no solve. **Not claimed:** any AED value or verdict (the 10-12 weekly reads rows 1 / 4 privately); an h-statement; copper; the tuned (ESR) fixture; any `ans:8` artifact change; an absolute C_mode, f₀ or tuning figure.

**Negative result.** Report and stop; park on `attempt/ANS-8-step1c-*` and mark item 366 BLOCKED in the same commit (rule (d)).
- **(F) or (G2′) red:** the file is not what entry 29 wrote. Stop; the weekly re-poses entry 29.
- **(L-d1) red with (F) green:** the az synthesis is wrong. Print the per-class table; park; do not read degree 2.
- **(L-d2) red with (L-d1) green:** the saved degree-2 matrix does not reproduce its window. File a known-issues entry; the 10-12 weekly does not read rows 1 / 4 from it.
- **(G1′) red:** a passivity or reciprocity finding on the degree-2 matrix. File a known-issues entry.
- Every row, move and printed control is a reading, never a red.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Status it can move.** No glyph (rule (j) does not apply). On green, the §7 `ANS-8` row gains one sentence (replacement text in the 2026-10-11 review's return); §10 H15c gains a DONE line; the dashboard's Waiting-on-you W2 ("H15c") can say the degree-2 rows exist for the 10-12 weekly's private read. Known-issues changes only on a red.

**Execution order.** Host `ls -la logs/wf7/` → W0 → W1 → commit (module, both harness logs, the §7 sentence, §10 H15c DONE) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 366 — `ANS-8` step 1c (§10 H15c) — per `docs/planning/chunks/ANS-8.md` §9 item 366. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. If `logs/wf7/wf7-s32-degree2-n16.npz` is absent on the host, stop and mark item 366 BLOCKED; never re-solve. Run W0, then W1, in order, foreground. Never delete a log. No edit to the module after W0 starts. Do not edit `ans:8`, its `metrics.json` / `COMPARISON.md`, or any `src/` file; never open `aed_results/`, `COMPARISON_private.md` or `docs/private/`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ANS-8-step1c-*`, mark item 366 BLOCKED in the same commit, and report."

### §9 item 369 — 2026-10-11 (§10 H21, `ANS-8` step 2 runnable half; chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain H, "H21. `ANS-8` step 2 runnable half — `ans:8` v3 rows from the SPEC section written tonight (the tuned two-port drive with ESR: mode share and leak at the root, `|B₁⁺|` at the five probes per √W accepted, the 10 g peak and location, phantom power fraction, the ESR term), every figure imported from items 358 / 359's records, nothing restated; `COMPARISON.md` with the AED columns blank. Heavy. Serial on nothing" (PROJECT_PLAN.md:4700–4705), and the weekly's commissioning paragraph (:4908–4913). The tracked SPEC v3 section is `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/SPEC.md:256–309`. The 10-10 weekly rulings on items 354 / 358 / 359 (`WF-10.md`, "Weekly rulings 2026-10-10") and the 10:30 review rulings beneath item 359 bind. **No band, tolerance, Done-when or target is set or moved.** The SPEC's readout percentages (1 %, 10 %, 15 %, 5 %) are the weekly's AED-adjudication levels; this item never asserts them.

**`ANS-8` step 2, our side: SPEC v3 rows 1–6 on H9′'s tuned copper F-human (degree 1, C_in,F-human = 118.2437 pF, r = 0.2 Ω on all 32 gaps, ccw two-port drive), per √W / W of the SPEC's accepted power, tied to items 358 / 359 by reproduction.**
- **Executor:** implementer, applying `example-runner`'s template rules (import, never restate; same-stem guide; censuses predicted first). Items 141 / 161 are the precedent. **Tier:** heavy, `-n 8`, complex, for stage V3; `-n 1` for W0, stage R and the pin test. **Class:** chain step (solve), median 14. **Cost:** ≈ M measured / ≈ 65 predicted min (15 fixed + ≈ 35 code + W0 ≈ 1 + V3 ≈ 5 + R ≈ 1 + pin and censuses ≈ 2 + spare).
- **45-minute rule (item 141's):** if stage V3 has not started by minute 45, commit the script, guide edits and unit tests with W0 as a partial, mark item 369 **BLOCKED "stage V3 pending"** in the same commit, and stop; the next slot runs V3 → R → pin → post-censuses from this brief.
- **Depends on nothing open.** Never run with `FEM_EM_WF10_DEGREE` set (item 368's knob; unset is the degree-1 path). Never concurrent with another heavy window.
- **Inputs (gitignored, on the host at drafting):** `logs/wf10/h9p-node-118.2437.json` (5 468 bytes, rewritten 2026-10-10 19:50 CDT by item 361's builder call); `logs/wf12/h12c-cu-s32-a2.npz` (20 812 bytes); `logs/ans8/s32.npz` (19 802 bytes, stage R's v1 / v2 input). **Any absent ⇒ stop, BLOCKED; never re-solve a 32×32.**
- **Family:** `ANS-8` step 2, attempt 1 of 4.

**Two premises the drafter settled from the code (the review checks them).**
1. **"Accepted power" is the SPEC's, not item 361's.** SPEC v3 says P_acc = incident − reflected at the two feeds, at the 50 Ω line. In our model each feed is a 50 Ω source in parallel with its own gap's r + C, reduced to a Thevenin pair (w_i, Z_th = 50 Ω ∥ Z_RC) (`tests/validation/test_wf10_f_human_copper_esr_root.py:76–83`). The line's accepted power therefore **includes** the feed gaps' ESR loss: P_acc,SPEC = ½ Re Σ_feeds V_port,i · conj(I_line,i), with V_port,i = w_i − Z_th·I_i and I_line,i = I_i + V_port,i / Z_RC. Item 361's P_acc excludes that term (its "caveat, printed", `MAT-15.md` item 361). Drafter's arithmetic from the node JSON (unlogged): P_acc,SPEC ≈ 2.906e-2 + 1.265e-3 ≈ 3.032e-2 W at |w| = 1 V.
2. **Row 2 on our side is the in-model field's sheet currents, not the circuit reduction of the stored 32×32.** SPEC row 2 says "step 1's row 2b algebra on the 32-port S at this C". On our copper fixture the circuit ↔ field identity with local loads is red at 1.58 (H12c, the pre-registered red; printed on every H17 node, e.g. `20261010T093200Z_WF-10-step1e2-W1a.log:6324`). Our row 2 is therefore `ring_modes.k_spectrum` applied to the builder's `cur["ccw"]`, which 358's (Lin) tied to the field to ≤ 6.6e-16. The circuit-route value on `h12c-cu-s32-a2.npz` (`_parallel_currents` with Z_RC loads and Z_th feeds, `:103–112`) is **printed beside it, labelled "circuit-only; H12c identity red at 1.58; not our row"**. **Flag for the 10-12 weekly:** SPEC v3 row 2's wording names the circuit route. Whether the AED column uses a 32-port export or the in-model sheets is the weekly's to settle in the SPEC; this item does not edit the SPEC's rows.

**Code.** `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py` and its same-stem guide. **No `src/` change, no edit to any `tests/validation/test_wf10_*` module** (import their module-level helpers only; never a `test_*` name).
- **New stage `FEM_EM_ANS8_STAGE=V3`** (`-n 8`), dispatched in `main()` beside S / R; `_stage_v3(comm, outdir)` with lazy imports:
  1. **Before the builder** (it rewrites the JSON, 359 `:173–182` / 361 step 1): rank 0 reads `NODE_JSON` (`tests/validation/test_wf10_f_human_tuned_b1_sar.py:82`), `bcast`, asserts `c_in_pf == C_IN_F_HUMAN_PF` (`:58`) exactly, keeps `s_node`, `y_in`, `partitions`.
  2. `out = build_tuned_quadrature(comm, with_z50=False)` (`:147`). It asserts (A0), (iv), the 90° feed spacing, (P_RC) and the node NEG, and returns `e`, `cur`, `w`, `parts`, `y_in`, `s_node`, `order`, `slots`, `feeds`, `az_feeds`, `az_mirror`, `ctx`, `sigma`, `anchors`.
  3. **Row 1:** `out["y_in"]["ccw"]` per feed (Re, Im in S), beside the JSON value.
  4. **Row 2:** per ring, `fem_em_solver.ports.ring_modes.k_spectrum` (`src/fem_em_solver/ports/ring_modes.py:100`) on `cur["ccw"]` mapped through `ring_port_map_from_names(order)` (feed ring = bottom = P17…P32). It reports share_k1, leak_k0 and the |k| folds. The circuit-route value is printed beside it, as above.
  5. **Power (P_acc):** P_acc,SPEC from the currents (premise 1). The partition is P_RC,30 = Σ `p_sheet_by_id` over the 30 non-feed sheets, P_ESR,feeds = Σ ½ r \|V_port,i / Z_RC\|², `p_phantom`, `p_surf`, `p_rest`. **Row 5:** P_phantom / P_acc, P_surf / P_acc (copper), (P_RC,30 + P_ESR,feeds) / P_acc (ESR).
  6. **Rows 3 and 6:** B_ccw through `magnetic_flux_density_from_e` (`src/fem_em_solver/post/faraday.py:45`), then `project_to_cg1` (`:69`), then `_read_senses` (`tests/validation/test_birdcage_b1_quadrature.py:214`). Read at the SPEC probes O (0, 0, z_c), P₁ (0.045, 0, z_c), P₂ (0, 0.045, z_c), P₃ (−0.045, 0, z_c), P₄ (0, −0.045, z_c), with z_c from `_phantom_geometry` (`tests/validation/test_wf7_f_human_sar10g_hotspot.py:95`; predicted 0.0, as 359's balls sat at z = 0, `20261010T171705Z_WF-10-step2s-W1.log:6314`). Row 3: \|B₁⁺\| / √P_acc,SPEC. Row 6: \|B₁⁺\| / \|B₁⁻\| at O. Also re-read 358's centre sample (`_sample_points`, `tests/validation/test_birdcage_b1_plus_map.py:246`) for the (Rep-B) anchor below.
  7. **Row 4:** 359's 18-point lattice (`_lattice`, `test_wf10_f_human_tuned_sar10g.py:58`) under ccw through `mass_averaged_sar` (`src/fem_em_solver/post/sar.py:231`) with `QUADRATURE_DEGREE`, `TEN_GRAM_KG` (`tests/validation/test_mass_averaged_sar_standard_masses.py:121, :124`) and ρ as 359 builds it. c\* = lattice argmax; SAR10g(c\*) / P_acc,SPEC; location in metres. The averaging volume is **a sphere** (the SPEC asks the codes to state which). Containment via `_containment` (`:85`).
  8. Rank 0 writes `/workspace/logs/ans8/v3-rows.json` (gitignored), then `Barrier`. Prints the price per phase and the summed RSS.
- **Stage R** (`-n 1`) gains one block after the v1 / v2 assertions: if `v3-rows.json` exists, load it, re-check the (P_acc) identities from its stored scalars, and add `rows_v3` plus `"spec": "SPEC.md v3 (2026-10-10)"` to `metrics.json`. If it is absent, print "v3 rows absent" and write v1 / v2 exactly as today. v1's `rows` and `rows_v2` are bit-unchanged.
- **`COMPARISON.md`** gains a section "Step 2 (SPEC v3) — tuned two-port drive with ESR, degree 1 (ours; AED column blank)": rows 1–6 in a `| Row | ours | AED |` table with the AED cells empty. It is labelled degree 1 on the 0.015 m phantom, c\* is labelled a lattice argmax, and rows are per √W / W of P_acc,SPEC. The step-1 tables stay as they are, labelled **"lossless (step 1)"** (SPEC :308–309).
- **The private writer is not extended** (no AED v3 results exist; the operator's session follows H21). `_write_private_if_available` is not called with v3 content; `FEM_EM_ANS_SKIP_PRIVATE=1` is set in every window anyway.
- **Unit tests** (append to `tests/unit/test_ring_mode_rows.py` or a new `tests/unit/test_ans8_v3_rows.py`; U1–U10 and P5–P7 untouched):
  - **V1 (P_acc algebra):** on a synthetic 2-feed network, ½ Re Σ V conj(I_line) equals the wave form (1/(8 z0))·(‖v_s‖² − ‖S v_s‖²), with v_s = w·`_source_voltage(Z_RC, z0)` and S from `_current_route_s`, to ≤ `Y0_CODING_BAND` (H17 `:55`).
  - **V2 (fractions):** the row-5 fractions of a synthetic partition sum to 1 exactly.
  - **V3 (comparison text):** `comparison_text` with a synthetic `rows_v3` puts the six rows in the table with every AED cell empty, keeps the v1 / v2 lines byte-identical, and contains no number from a planted AED sentinel.

**Anchors (asserted in stage V3, every band imported).**
- **(S_rep)** `_residual(out["s_node"], json s_node)` ≤ `REUSE_REPRODUCTION_RTOL` (`_residual`, `tests/validation/test_port_circuit_layer_field.py:142`). Same comparison, same fixture: 359 read 6.027e-13, 361 read 2.705e-13.
- **Builder anchors**, re-asserted from `out["anchors"]`: (A0) 369 387 cells, (P_RC) ≤ `DISCRETE_IDENTITY_RTOL`, node NEG ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL`, (P_w) ccw ≤ `DISCRETE_IDENTITY_RTOL`.
- **(P_acc-1) line cut = wave form:** P_acc,SPEC from the currents vs (1/(8 z0))(‖v_s‖² − ‖s_node v_s‖²) ≤ `DISCRETE_IDENTITY_RTOL`. Predicted ≤ 1e-12 ((Lin) ≤ 6.6e-16 ties the currents to `s_node`).
- **(P_acc-2) energy identity at the line cut:** P_acc,SPEC = P_RC,30 + P_ESR,feeds + P_phantom + P_surf + P_rest ≤ `DISCRETE_IDENTITY_RTOL`. This is (P_w) moved across the feed gaps. Predicted ≤ 1e-12. Row 5's fractions then sum to 1 by construction; that sum is printed.
- **(Rep-B)** mean \|B₁⁺\|_ccw over 358's masked centre sample reproduces 358's 2.884091e-07 T (`20261010T170356Z_WF-10-step2-W1.log:6339`) to ≤ `REUSE_REPRODUCTION_RTOL`. Record constant `MEAN_B1P_CCW_H9P_RECORD_T`, with the citation; 7-digit print rounding ≤ 1.8e-7.
- **(Rep-S)** SAR10g at 359's c\* (r = 45 mm, z_c, az₁₇) reproduces 4.852895959e-03 W/kg (`20261010T171705Z_WF-10-step2s-W1.log:6316`) to ≤ `REUSE_REPRODUCTION_RTOL`. Record constant `SAR10G_CSTAR_H9P_RECORD_W_PER_KG`. c\* equals 359's lattice point.
- **(Rep-Y)** row 1 `y_in["ccw"]` equals the JSON's to ≤ `REUSE_REPRODUCTION_RTOL` relative to \|Y_in\|.
- **(Sense)** \|B₁⁺\| / \|B₁⁻\| at O > 1 (a sign check, 358's).
- **(iv)** every ball clears the phantom by `CONTAINMENT_MARGIN_M`.

*Imported, nothing restated:* the three record constants above are 361-3-type comparands (logged readings, cited). Every published v3 number comes from the computation, never from a constant.

**Negative controls (rule (e)).**
- **Asserted (backed, same fixture):** the node NEG inside the builder. 359 / 361 read 7.0125e-3 against the bar 10 × 1e-6 (`…step2s-W1.log:6304`).
- **Asserted (backed by 358's same comparison):** the P_surf-blind miss of (P_w) ccw ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL`; 358 read 7.0397e-3 (`…step2-W1.log:6320`).
- **Predicted, printed (no prior of the same comparison):**
  - the ESR-blind miss of (P_acc-2) (drop P_RC,30 + P_ESR,feeds), predicted ≈ 0.71 (ceiling 1);
  - the feed-ESR-blind miss (P_ESR,feeds dropped), predicted ≈ 4.2e-2;
  - item 361's P_acc in place of P_acc,SPEC in (P_acc-1), predicted ≈ 4.2e-2.

**Predictions (written first, printed).** Drafter's arithmetic from the node JSON and the 358 / 359 logs, unlogged:
- P_acc,SPEC ≈ 3.03e-2 W; row 5 fractions: phantom 0.270, copper 0.020, ESR 0.710 (each ± 0.01);
- row 1 Im Y_in ccw: P17 ≈ +1.02e-5 S, P21 ≈ −2.16e-4 S (the JSON);
- row 2 field route: feed-ring share ≈ 0.9938 (358's unsigned ±1 on ring −, `:6330`), leak ≤ 5e-3;
- row 3 \|B₁⁺\|(O) ≈ 1.66e-6 T/√W (2.884e-7 T / √3.03e-2), P₁ … P₄ within ± 10 % of O;
- row 4 ≈ 0.160 W/kg per W at c\* = (44.135, 8.779, 0) mm;
- row 6 ≈ 21.5 (358's centre purity 21.5146, `:6334`).

**Windows.** Foreground, bare, repo-relative; Bash `timeout: 660000`; `-s`; `FEM_EM_ANS_SKIP_PRIVATE=1` in **every** window; complex windows source `/usr/local/bin/dolfinx-complex-mode` with `FEM_EM_REQUIRE_COMPLEX=1` and `tests/environment` first.
- **(pre, before any edit to the script or guide)** `ANS-8-v3-census-docrefs-before` / `ANS-8-v3-census-setupfig-before`. Use the exact commands at `20261002T050048Z_ANS-8-v2-census-docrefs-before.log:12` / `20261002T050050Z_ANS-8-v2-census-setupfig-before.log:12`, chunk labels changed only. Predict from the current census; never hard-code a count.
- **(inputs, host)** `ls -la logs/wf10/h9p-node-118.2437.json logs/wf12/h12c-cu-s32-a2.npz logs/ans8/s32.npz`.
- **W0** (≈ 60 s; item 161's W0 command, chunk `ANS-8-v3-W0`, with the new unit file added to the path list):
  `scripts/testing/run_and_log.sh ANS-8-v3-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 180 python3 -m pytest tests/environment tests/unit/test_ring_mode_rows.py tests/unit/test_ans8_private_comparison.py tests/unit/test_ans8_v3_rows.py -s -v --tb=short'"`
- **Before V3:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **V3** (heavy, `-n 8`; §5.1 durable capture verbatim, trailing `; exit $rc` included):
  `scripts/testing/run_and_log.sh ANS-8-v3-V3 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/ans8 && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS8_STAGE=V3 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py > /workspace/logs/ans8-v3-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans8-v3-raw.log; cat /workspace/logs/ans8-v3-raw.log; exit \$rc'"`
  - **Price, from measured parts:** builder with `with_z50=False` 152.7 s (361, `20261011T003520Z_MAT-15-step5-W1.log`); B projection, points and power identity ≈ 10 s (358's B phase 3.3 s, `…step2-W1.log`); 18 balls at 3.44 s mean ≈ 62 s (359, `…step2s-W1.log`); harness ≈ 10 s. Total ≈ 235 s, ≤ 310 s at the 1.32× box reading, inside 590 s. Summed RSS ≈ 10 GiB.
  - The stage must exit non-zero on any red, as stage R's `fails` list does.
- **W-R** (stage R, `-n 1`; R:12's command, i.e. `20261002T050705Z_ANS-8-v2-W1-R.log:12`, chunk `ANS-8-v3-R`, `timeout -k 30 120`). Predicted ≈ 3 s; every v1 / v2 check prints exactly as at R:35–79.
- **W-pin:** `20261002T050738Z_ANS-8-v2-W2-pin.log:12`'s command, after `git add` of `metrics.json` and `COMPARISON.md` (20 passed there).
- **(post)** Both censuses again. Predicted unchanged N → N: no new script, the guide is the same stem.
- **D1 (host `git diff --stat`, not a window):** `metrics.json` has insertions plus the one changed `spec` line; `COMPARISON.md` v1 / v2 lines unchanged (insertions only).

**Traps already paid for.**
- **The builder rewrites `h9p-node-118.2437.json`.** Read it first (item 359 / 361).
- **Sense:** ccw = (1, −1j) from `quadrature_phase_weights`; never re-chosen; the mirror plane is az₁₇ + 45°.
- **Accepted power:** P_acc,SPEC (premise 1), never item 361's P_acc and never P_src. Print all three, so the per-W figures compare with 358 / 359 / 361.
- **The circuit route on the stored 32×32 is not our row 2** (premise 2).
- **B through `project_to_cg1`** (an L² projection); points only through `evaluate_vector_field_parallel` inside `_read_senses`; never `f.eval` (known-issues 2026-09-27).
- **Every rank enters every collective** (balls, `_phantom_geometry`, `_read_senses`); never guard one with "owns phantom cells" (item 227). Rank 0 alone writes files, then `Barrier`.
- **Complex build:** no `ufl` ordering or `sqrt` comparison on complex operands; ball centres as plain floats; `np.real` before every `float`.
- **Orphan `.c` after a kill:** `pgrep -c python3` first; at 0, run the rubric-4 sweep under label `ANS-8-v3-cache-sweep`. Never raise the timeout or `-n`.
- **Pre-census contamination:** census before editing (item 141's first pre-census ran with the script already in the tree).
- **The operator's private file** stays untouched: knob set everywhere, P7 guard in place, `aed/` never imported, edited or run.
- **Import, never restate:** no reading of 358 / 359 is retyped into the published rows. The three record constants are comparands only.
- Never import a `test_*` name; never `-k`; no pipe; never delete a log; no edit to the script after V3 starts (rule (i)); never path-checkout a record file.

**Scope.** Our side of SPEC v3 rows 1–6 at degree 1, on H9′'s copper F-human: hole mesh, κ off, r = 0.2 Ω on all 32 gaps, ccw at C_in,F-human = 118.2437 pF, 64 MHz, the 0.015 m phantom (≈ 1.8 cells per 10 g ball).
- **Not claimed:**
  - any AED value, verdict or agreement level;
  - the degree-2 column (`xl` entry 35);
  - f₀, Q, match or any tuning-accuracy figure (C_in is a record at 0.10 pF resolution);
  - c\* as a global maximum (it is a lattice argmax);
  - homogeneity;
  - compliance;
  - the private v3 writer;
  - a SPEC edit (the row-2 wording is flagged, not changed).

**Negative result.** Report and stop. Park on `attempt/ANS-8-v3-*` and mark item 369 BLOCKED in the same commit (rule (d)).
- **(S_rep), a builder anchor, (Rep-Y), (Rep-B) or (Rep-S) red:** the fixture moved since 358 / 359. Write no metrics; file a known-issues entry for a builder anchor.
- **(P_acc-1) red:** the line-cut algebra or the source convention is wrong (V1 should have caught it). Park.
- **(P_acc-2) red with (P_w) green:** the feed-ESR term is mis-assembled. Park.
- **(Sense) red:** a convention finding; the review rules.
- **(iv) red:** stop.
- **Stage R v1 / v2 check red:** a regression. Stop before writing.
- **V3 killed or status 124:** not a result. Apply the 45-minute rule; the killed window is the price record.
- Every row value, fraction and circuit-route print is a reading, never a red.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Outcome 2026-10-11 09:00 slot (take-next): 🚫 BLOCKED, attempt 1 of 4 — (P_acc-2) red.** Pre-censuses `20261011T141113Z_ANS-8-v3-census-docrefs-before.log` (dead 0, guide 0, stale 76, report-severity exit 2) / `…141115Z_…setupfig-before.log` (105 / ok 105). W0 `…141506Z_ANS-8-v3-W0.log` 32 passed, 78 s (V1 line cut = wave form ≤ 1.5e-15 on four drives, item-361 form misses ≥ 1.2e-3; V2 sum exact; V3 six rows, AED cells empty, v1 / v2 prefix unchanged). `pgrep` 0. V3 `…141650Z_ANS-8-v3-V3.log` 202 s harness (199.2 s in-script: builder 160.2, B + points 3.9, 18 balls 34.7; summed RSS 9.58 GiB, `-n 8`), exit 1. **Green:** (S_rep) 3.092e-13, (A0) 369 387, (P_RC) 7.850e-13, node NEG 7.0125e-3, (P_w) ccw 5.055e-13, its P_surf-blind NEG 7.0397e-3, (Rep-Y) 2.280e-13, names map = measured slot map, (P_acc-1) 6.808e-16, (Rep-B) 1.101e-7 (41 masked points), probes valid, (Sense) 20.7092, (Rep-S) 4.362e-11 with c\* = 359's point, (iv) worst −4.8813 mm. **Red:** (P_acc-2) 8.317e-3 against `DISCRETE_IDENTITY_RTOL`. P_acc,SPEC (line cut) 3.057630e-2 W; partition sum P_RC,30 2.024950e-2 + P_ESR,feeds 1.265338e-3 + P_phantom 8.197883e-3 + P_surf 6.092762e-4 + P_rest 0 = 3.032186e-2 W *[corrected by the 2026-10-11 10:30 review on the `log-pathologist`'s ruling: the terms at `…V3.log:6195` add to **3.032200e-2 W**; the row-5 sum 0.991683 at `:6199` confirms it]*. **Cause, measured (diagnostic `…142050Z_ANS-8-v3-diag.log`, arithmetic on the node JSON's `y_in` / `w` / partition, no solve):** premise 1 takes the field's P_src − P_feed sheets (item 361's P_acc, 2.905666e-2 W) to equal the terminal power ½ Re Σ V_port·conj(I) from the sheet currents (2.931096e-2 W); they differ by 2.543e-4 W (8.7e-3). Field P_src 8.654821e-2 vs Thevenin ½ Re Σ w·conj(I) 8.016767e-2 W; field feed-sheet dissipation 3.036767e-2 / 2.712387e-2 vs ½ Re(Z_th)\|I\|² 2.655324e-2 / 2.430347e-2 and the Norton form ½ Re(1/Z_th)\|V\|² 3.023474e-2 / 2.700251e-2 W (P17 / P21). The field's feed-sheet accounting is not the circuit's terminal V·I\* at 1e-6; (P_acc-2) as posed mixes the two. The drafter's 3.03e-2 prediction was the partition sum, which is why it matched the partition and not the line cut. Negative controls printed: ESR-blind 0.71196 (predicted ≈ 0.71), feed-ESR-blind 4.970e-2 and item 361's P_acc in (P_acc-1) 4.970e-2 (both predicted ≈ 4.2e-2). Readings (unpublished, no metrics written): row 2 feed ring 0.993845 / leak 4.739e-3, other ring 0.998448 (circuit-only 0.980957 / 1.527e-2); row 3 \|B₁⁺\|(O) / √P_acc 1.668096e-6 T/√W, P₁…P₄ 1.6107–1.6528e-6; row 4 0.1587143 W/kg/W at (0.044135, 0.008779, 0) m; row 5 phantom 0.268112, copper 0.019926, ESR 0.703644 (of the line-cut P_acc; sum 0.991683); row 6 20.7092. Code parked on `attempt/ANS-8-v3-20261011T1420Z` (script stage V3, R block, `_v3_section`, `tests/unit/test_ans8_v3_rows.py`); no post-census (nothing on main changed under `examples/`). **Unblock:** the review rules which identity (P_acc-2) states — e.g. line cut = ½ Re Σ V·I\* (terminal) + P_ESR,feeds, with the field partition held to (P_w) as now, or the line cut re-based on the sheets' field power — and re-poses; the code needs no re-solve beyond one V3 window.

**Status it can move.** No glyph. The `ANS-8` Done-when ("✅ only after the AED half lands and the weekly adjudicates") is unaffected, so rule (j) does not bind. On green:
- the §7 row gains "step 2 runnable half (v3 rows) landed" (replacement text in the 2026-10-11 review's return);
- the SPEC Status list gains a ticked line "Runnable half, v3 rows (2026-10-…, §9 item 369)" (item 161 ticked line 28 the same way; no other SPEC edit);
- §10 H21 gains a DONE line;
- **the dashboard's Waiting-on-you list gains `ANS-8` step 2: the operator's one AED session (SPEC v3)**, per the weekly (:4912–4913).

**Execution order.** Pre-censuses → inputs `ls` → code → W0 → `pgrep` → V3 → W-R → `git add` artifacts → W-pin → post-censuses → D1 → commit (script, guide, unit tests, `metrics.json`, `COMPARISON.md`, SPEC status line, harness logs, the §7 sentence, §10 H21 DONE) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 369 — `ANS-8` step 2 runnable half (§10 H21) — per `docs/planning/chunks/ANS-8.md` §9 item 369. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Set FEM_EM_ANS_SKIP_PRIVATE=1 in every window; never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. Run the pre-censuses before any edit. If `logs/wf10/h9p-node-118.2437.json`, `logs/wf12/h12c-cu-s32-a2.npz` or `logs/ans8/s32.npz` is absent, stop (BLOCKED); never re-solve. Before stage V3, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. If V3 has not started by minute 45, commit the partial and mark item 369 BLOCKED 'stage V3 pending'. Never set FEM_EM_WF10_DEGREE. Never delete a log. No edit to the script after V3 starts. Do not edit any `tests/validation/test_wf10_*` module or any `src/` file. Every rank enters every collective. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ANS-8-v3-*`, mark item 369 BLOCKED in the same commit, and report."

**Review ruling 2026-10-11 10:30 on item 369 (provisional, for the 10-12 weekly).** A `log-pathologist` ruling (CONFIRMED on all four sub-questions) found the following.

- **The red is not the feed-ESR term.** P_ESR,feeds = ½ r Σ\|V/Z_RC\|² sits identically on both sides of (P_acc-2) and cancels, so no consistent mis-assembly of it can produce a residual. This ruling **overrides the brief's routing line** "(P_acc-2) red with (P_w) green ⇒ the feed-ESR term is mis-assembled".
- **The red is premise 1's tacit equality**, field net (P_src − P_feed sheets) = terminal ½ Re Σ V·I\*, and that equality is false by 2.5430e-4 W. The pathologist's arithmetic on `…diag.log:37–39` (not printed by any run) locates it exactly:
  - The lumped Norton source power, terminal + Σ ½ Re(1/Z_th)\|V_i\|², equals the field P_src (8.654820635e-2 W) to every printed digit. The terms linear in the field agree with the circuit.
  - The whole gap is the field feed-sheet loss in excess of ½ Re(1/Z_th)\|V\|²: 1.32932e-4 W on P17 and 1.21365e-4 W on P21.
  - That is the shape of the Cauchy–Schwarz excess ½ Re(Y_s)(∫\|n̂×E\|² dA − \|∫E·ĥ dA\|²/A) that a field varying across the sheet produces (`_sheet_total_field_dissipation_w`'s docstring, `tests/validation/test_birdcage_power_identity.py:272–293`). The mechanism is unmeasured.
  - The diag's Thevenin columns (½ Re(Z_th)\|I\|², ½ Re(w·I\*)) are internal quantities of a different equivalent, and their 7–14 % differences are not evidence.
- **Correction:** the partition sum is 3.032200e-2 W, not 3.032186e-2 (corrected inline above). The gap is 8.317e-3 of the line cut and 8.75e-3 of the field net, matching item 361's 8.75e-3 (`MAT-20.md`, "Which accepted watt"). The same fixture and the same term have now been measured twice.
- **Ruling:**
  - (a) **P_acc,SPEC (the line cut) stays the per-W normaliser of rows 3–5.** It is the terminal definition a lumped-port code reports, and (P_acc-1) ties it to the wave form at 6.8e-16.
  - (b) **(P_acc-2) is withdrawn as posed.** No assertion equates field net with terminal power; the field partition stays held by (P_w).
  - (c) Row 5's fractions of P_acc,SPEC are published with their sum (0.991683) and the residual labelled "feed-sheet field-vs-terminal gap, degree 1 (known-issues 2026-10-11)". Fractions of the field partition are printed beside them.
  - (d) The gap's mechanism is measured in the same window, printed and never asserted.
  - This matches the 03:00 review's ruling 2 for P20: both normalisations are printed and neither is asserted as the other.
- **Family:** `ANS-8` step 2, attempt 2 of 4 (item 386). The target is H21's, unchanged, and no band moves.

### §9 item 386 — 2026-10-11 (§10 H21 re-pose, `ANS-8` step 2 runnable half, attempt 2 of 4; chain step (solve))

*Written by the 2026-10-11 10:30 daily review under rule (6)'s re-pose clause (same target, same family count).* The item 369 brief above binds in full — inputs, code layout, anchors, traps, scope, windows, execution order and spawn prompt — **except** where this section amends it. Read the 10:30 ruling directly above first.

**Land `attempt/ANS-8-v3-20261011T1420Z` (`6e005c22`) with three amendments, then run item 369's windows.**
- **Take by path** the script `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py` and `tests/unit/test_ans8_v3_rows.py` (`git checkout attempt/ANS-8-v3-20261011T1420Z -- <path>`, code files only). **Do not take `scratch_v3_diag.py`**: amendment 3 supersedes it. Never path-checkout a record file.
- **Amendment 1 — (P_acc-2) withdrawn.** Delete its assertion and its `fails` entry. Keep its two numbers as prints: "line cut − partition = gap", relative to the line cut and to the field net.
- **Amendment 2 — reproduction anchors replace it, asserted at `REUSE_REPRODUCTION_RTOL` (1e-6, the constant (Rep-B) imports):**
  - **(Rep-P)** P_acc,SPEC reproduces 3.057629563251e-02 W (`…V3.log:6194`).
  - **(Rep-Part)** P_RC,30, P_ESR,feeds, P_phantom and P_surf reproduce `…V3.log:6195` term by term (9-digit print, rounding ≤ 1e-9 relative).
  - Same code, same mesh and same fixture as item 369's V3. (S_rep) read 3.092e-13 there, so a reproduction at 1e-6 is a fixture check, not a tautology.
  - Record constants: `P_ACC_SPEC_H21_RECORD_W` and `PARTITION_H21_RECORD_W` (a dict), each with its log:line citation.
- **Amendment 3 — the gap's mechanism, printed, pre-registered and never asserted.** On each feed sheet s ∈ {P17, P21}, with E = `out["e"]["ccw"]` and every rank entering every collective:
  - A_s = ∫ 1 dA over the sheet.
  - J_s = ∫ E·ĥ dA, with `ufl.dot` and ĥ = `sheet.unit_drive()`, real and never conjugated.
  - L_s = ∫ \|n̂×E\|² dA, the integrand `_sheet_field_dissipation_w` already carries.
  - T_s = ∫ \|E·ĥ\|² dA.
  - Print CS_s = ½ Re(Y_s)(L_s − \|J_s\|²/A_s) with Y_s = 1 / `sheet.sheet_resistivity`, beside the excess X_s = P_sheet,field(s) − ½ Re(1/Z_th)\|V_s\|².
  - Print the split ½ Re(Y_s)(L_s − T_s) (the in-plane component across the drive) and ½ Re(Y_s)(T_s − \|J_s\|²/A_s) (the drive-component variance).
  - Also print P_src against the Norton lumped power, terminal + Σ ½ Re(1/Z_th)\|V\|².
  - **Predictions (the pathologist's arithmetic):** X = 1.32932e-4 / 1.21365e-4 W; CS_s / X_s in [0.99, 1.01] on each sheet if the field-variation reading holds; P_src = Norton lumped to ≤ 1e-12 relative.
  - Assembled forms go through `_reduce` (`test_birdcage_power_identity.py:248`), never a rank-local `assemble_scalar`. Measures and the restricted side follow `_sheet_total_field_dissipation_w` (`:272–293`): `("+")` restriction, `ufl.dot` for the real normal.
- **Rows:**
  - Row 5 is published per (c) of the ruling: fractions of P_acc,SPEC, their sum, and the labelled residual. The partition-normalised fractions are printed beside them, and `COMPARISON.md`'s v3 section carries both columns.
  - Rows 3 and 4 are per √W / W of P_acc,SPEC, with the partition-normalised value printed beside each.
  - `metrics.json`'s `rows_v3` stores both normalisations, keyed `per_p_acc_spec` and `per_partition`.
- **Unit tests:** V1–V3 unchanged. Add **V4**: on a synthetic sheet field with a known transverse component and a known drive-component variance, the split's two terms recover both to ≤ 1e-12. This is pure numpy on sampled values, the discrete analogue of amendment 3's algebra.

**Anchors (asserted):** every item 369 anchor except (P_acc-2): (S_rep), the builder anchors ((A0), (P_RC), node NEG, (P_w)), (P_acc-1), (Rep-B), (Rep-S), (Rep-Y), (Sense), (iv). Add (Rep-P) and (Rep-Part). The bands are imported, and none moves.

**Negative controls (rule (e)):**
- **Asserted, now backed on this fixture:**
  - The node NEG and the P_surf-blind miss of (P_w), as in item 369.
  - **Item 361's P_acc in place of P_acc,SPEC misses (P_acc-1)** by ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL`. Item 369 read 4.9700e-02 on this comparison (`…V3.log:6198`).
- **Printed:** the ESR-blind miss (0.71196 at `:6198`) is no longer attached to an asserted identity, so it is printed only.

**Tier, ranks, cost:**
- W0 (`-n 1`, ≈ 80 s; item 369's W0 measured 78 s).
- V3 (`-n 8`, heavy): item 369's window measured 202 s. Amendment 3 adds four facet integrals per sheet, ≈ 2 s. Predicted ≤ 230 s, inside `timeout -k 30 590`.
- W-R and W-pin are as item 369's.
- **Class:** chain step (solve), median 12. **Cost:** ≈ 12 measured / ≈ 30 predicted min. The code exists, so the cost is three amendments, W0, V3, R, the pin and the censuses.
- **Depends on nothing.** Never beside another heavy window, and never with `FEM_EM_WF10_DEGREE` set.

**Traps already paid for:** item 369's list, plus the following.
- `log-pathologist` correction: the partition sum is 3.032200e-2 W. Never type 3.032186e-2 into a constant.
- The builder rewrites `h9p-node-118.2437.json`, and item 368's W2 rewrote it at 09:08. Read it before the builder.
- The smoke tier's environment suite (ruling 2 of this review): W0 is standard tier, so the full `tests/environment` stays in its path list.

**Scope:** item 369's scope; the "Not claimed" list holds. The following are not claimed either:
- the mechanism of the gap as a gate (printed only);
- any statement about the gap at degree 2 or on another fixture;
- which normalisation the AED column should read (the 10-12 weekly's, with SPEC v3's row-2 wording).

**Negative result:** report and stop.
- **(Rep-P) or (Rep-Part) red:** the fixture or the solver moved since item 369's window. Park and file nothing new; the review compares the two logs.
- **Any other anchor red:** item 369's dispositions.
- **CS_s / X_s outside [0.99, 1.01]:** a reading, never a red. It re-heads the known-issues entry 2026-10-11 ("mechanism not the sheet variance").

No `PREDICTED RED:` line.

**Status it can move:** as item 369's — the §7 row's step 2 sentence, the SPEC Status tick, §10 H21 DONE, and the dashboard Waiting-on-you item `ANS-8` step 2 (the operator's AED session, SPEC v3). In addition, the known-issues entry 2026-10-11 (feed-sheet field vs terminal power) gains its mechanism reading. It **retires** only if the weekly rules the gap a defined, documented sheet property rather than an issue.

**Execution order:** item 369's, with "take by path + amendments 1–3 + V4" in place of "code". Post-censuses are predicted unchanged: the same stem, no new script. The commit carries the script, the unit tests, `metrics.json`, `COMPARISON.md`, the SPEC status line, the logs, the §7 sentence, §10 H21 DONE, and the known-issues mechanism line. Leave the attempt branch for the review to delete.

**Executor spawn prompt (verbatim).** "Execute §9 item 386 — `ANS-8` step 2 runnable half (§10 H21, re-pose, attempt 2 of 4) — per `docs/planning/chunks/ANS-8.md` §9 item 386, which amends §9 item 369's brief above it; read the 2026-10-11 10:30 review ruling, item 369's brief and item 386's section first. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Set FEM_EM_ANS_SKIP_PRIVATE=1 in every window; never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. Run the pre-censuses before any edit. Take the script and the unit-test file by path from `attempt/ANS-8-v3-20261011T1420Z`, never `scratch_v3_diag.py`, never a record file. If `logs/wf10/h9p-node-118.2437.json`, `logs/wf12/h12c-cu-s32-a2.npz` or `logs/ans8/s32.npz` is absent, stop (BLOCKED); never re-solve. Before stage V3, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. If V3 has not started by minute 45, commit the partial and mark item 386 BLOCKED 'stage V3 pending'. Never set FEM_EM_WF10_DEGREE. Never delete a log. No edit to the script after V3 starts. Do not edit any `tests/validation/test_wf10_*` module or any `src/` file. Every rank enters every collective. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ANS-8-v3b-*`, mark item 386 BLOCKED in the same commit, and report."

**Outcome 2026-10-11 12:00 slot: ✅ landed, attempt 2 of 4.** Script and `tests/unit/test_ans8_v3_rows.py` taken by path from `attempt/ANS-8-v3-20261011T1420Z` (no `scratch_v3_diag.py`), amendments 1–3 + V4 applied, no `src/` or `test_wf10_*` edit. Pre-censuses `20261011T170124Z_ANS-8-v3b-census-docrefs-before.log` (dead 0, guide 0, stale 85, report-severity exit 2) / `…170129Z_…setupfig-before.log` (105 / ok 105). Inputs present (node JSON 5 463 B, rewritten 09:19 CDT; both npz). W0 `…170500Z_ANS-8-v3b-W0.log` 33 passed, 75 s (V4: transverse rel 8.5e-14, variance rel 0, CS rel 8.4e-15). `pgrep` 0. V3 `…170627Z_ANS-8-v3b-V3.log` 197 s harness, exit 0 (194.3 s in-script: builder 151.7, mechanism 5.3, B + points 3.4, 18 balls 33.6; summed RSS 9.88 GiB, `-n 8`). **Asserted, green:** (S_rep) 5.013e-13, (A0) 369 387, (P_RC) 7.720e-13, node NEG 7.0125e-3, (P_w) ccw 4.839e-13 and its P_surf-blind NEG 7.0397e-3, (Rep-Y) 3.750e-13, names map = slots, (P_acc-1) 3.404e-16, **(Rep-P) 1.205e-13**, **(Rep-Part) P_RC,30 8.497e-11 / P_ESR,feeds 3.660e-10 / P_phantom 2.699e-12 / P_surf 6.924e-11**, item 361's P_acc in (P_acc-1) NEG 4.9700e-2, (Rep-B) 1.101e-7, probes valid, (Sense) 20.7092, (Rep-S) 4.337e-11 with c\* = 359's point, (iv) worst −4.8813 mm. **Printed:** gap 2.542976202e-4 W = 8.3168e-3 of the line cut, 8.7518e-3 of the field net; ESR-blind 0.71196. **Mechanism (pre-registered, printed, never asserted):** P17 CS_s 1.329324349e-4 W vs X_s 1.329324349e-4 W, **CS_s / X_s 1.000000000** (split: transverse 9.401572e-5 + drive variance 3.891671e-5 W); P21 CS_s 1.213651853e-4 vs X_s 1.213651853e-4 W, **CS_s / X_s 1.000000000** (8.562073e-5 + 3.574446e-5 W); ΣX_s = the gap; P_src = Norton lumped 8.654820635490e-2 W to 3.3e-16; A_s = w·h on both sheets. The pathologist's predictions (X 1.32932e-4 / 1.21365e-4 W, ratio in [0.99, 1.01], Norton ≤ 1e-12) all hold. **Rows (per W of P_acc,SPEC = 3.057630e-2 W; per W of the field partition beside):** row 1 Y_in ccw P17 +1.010006e-2 +1.017991e-5j S, P21 +1.184993e-2 −2.164135e-4j S; row 2 feed-ring share 0.993845 / leak 4.739e-3, other ring 0.998448 (circuit-only 0.980957 / 1.527e-2, not our row); row 3 \|B₁⁺\|(O) 1.668096e-6 T/√W (1.675076e-6), P₁…P₄ 1.6107–1.6528e-6; row 4 0.1587143 W/kg/W (0.1600454) at (0.044135, 0.008779, 0) m, sphere a = 13.365 mm; row 5 phantom / copper / ESR 0.268112 / 0.019926 / 0.703644, sum 0.991683, residual 0.008317 = the gap (of the partition: 0.270361 / 0.020094 / 0.709546, sum 1); row 6 20.7092. W-R `…171000Z_ANS-8-v3b-R.log` (2 s, every v1 / v2 check PASS, (v3) (P_acc-1) 3.404e-16 from the stored scalars). W-pin `…171035Z_ANS-8-v3b-pin.log` 20 passed. Post-censuses `…171052Z_…docrefs-after.log` / `…171054Z_…setupfig-after.log` unchanged (stale 85 → 85; 105 → 105). **D1 reading:** `COMPARISON.md` is insertions only (v1 / v2 lines unchanged), but `metrics.json`'s v1 / v2 floats are **not bit-unchanged**: stage R's recomputation moved them at ≤ 1.4e-10 relative (7.5e-8 absolute, last digits) against the 2026-10-02 file, with every v1 / v2 check green — a stage R reproducibility reading on today's container, not caused by the v3 block (which only adds `rows_v3` and the `spec` line); recorded, not diagnosed. Attempt branch left for the review to delete.

### §9 item 411 — 2026-10-13 (§10 chain H22; lane 1, no solve)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* **Licence:** §10 H22 (PROJECT_PLAN.md:4696–4701), the 10-12 weekly's accepted-power ruling (:4552–4559) and its flag 4 (:4929). **No band, tolerance, Done-when or target is set or moved. No SPEC row, readout, fixture value or definition changes.** No stage of `ans:8` runs.

**`ANS-8` step 2 to the operator.** Two text deliverables. Both are render-only or insertions-only.
1. **`COMPARISON.md`'s step-2 section, re-stated with SPEC v3's accepted-power definition** (`SPEC.md:282–295`, added by the 10-12 weekly).
   - The change goes through the tracked writer `_v3_section` (`08_f_human_quadrature_mode_content.py:1082–1132`).
   - The file is rendered from the **committed** `metrics.json`, never hand-edited. Item 400 is the precedent.
2. **The operator's one-session brief**, as a new subsection of SPEC.md's step-2 section, inserted before `## Revision history` (`SPEC.md:327`). It comes with one Status checkbox and one revision-history line.
   - The slot **returns** a W2 text of 4 lines or fewer for the review to paste.
   - It does **not** edit `docs/status/dashboard.md`: moving W2 to the top of Waiting-on-you is the review's job (flag 4).

**Executor and cost.**
- **Executor:** implementer, or the slot directly (item 400's precedent for a writer-only change).
- **Tier:** standard (container `timeout -k 30 180`), no MPI (`-n 1`), complex build, because the script's import pulls in dolfinx.
- **Class:** lane 1 (median 6).
- **Cost:** ≈ 6 measured / ≈ 30 predicted min: 15 fixed + ≈ 10 for the writer, SPEC and test + ≈ 3 for windows + spare. The window times come from item 386's W0, 33 passed in 75 s (`20261011T170500Z_ANS-8-v3b-W0.log`), and item 400's render-only W1 at 3 s.
- **Depends on nothing.**
- **Inputs, all tracked:** `metrics.json` `rows_v3`. It carries `power_w.p_acc_spec` / `p_acc_wave` and `gap.rel_line_cut` (0.008317), `gap.rel_field_net` and `gap.mechanism.{P17,P21}.cs_over_x`. Read on the host while drafting.

**Edit 1: the writer (`_v3_section` only).**
- **Replace the intro clause** "Per-W / per-√W figures are per W of P_acc = incident − reflected at the two 50 Ω feed lines (… W at unit Thevenin drive)" with a definition paragraph in the SPEC's terms. It states five things:
  - **What P_acc is.** P_acc is the sum over the two feed ports of the accepted power at the 50 Ω lumped port: incident − reflected at the port plane, HFSS's default per-port accepted power.
  - **How our side computes it.** Ours is the line cut ½ Re Σ_feeds V_port·conj(I_line), rendered from `power_w["p_acc_spec"]`. Beside it goes its relative difference from `power_w["p_acc_wave"]`, computed from the two stored floats and never typed.
  - **The feed-gap load.** The feed gaps' capacitor and 0.2 Ω ESR are load inside P_acc; they appear as row 5's ESR fraction.
  - **The field partition.** The field's feed-sheet partition exceeds the line cut by `gap["rel_line_cut"]`, rendered as a percentage to two decimals. That excess is the discrete feed sheet's Cauchy–Schwarz excess. Render CS/X from `cs_over_x` on both sheets. It is a documented sheet property (weekly ruling 2026-10-12; known-issues 2026-10-11), not a defect, and not a quantity the AED column is asked for.
  - **The per-partition values** stay printed beside each row, as now.
- **Re-label row 5's residual line** as "residual = feed-sheet Cauchy–Schwarz excess (documented sheet property, weekly 2026-10-12; known-issues 2026-10-11)". The 0.991683 / 0.008317 values do not change.
- **Nothing else changes:**
  - not a v1 / v2 line;
  - not `comparison_text`'s other branches;
  - no stage logic;
  - not `metrics.json` (its stored `gap.label` and `spec` strings stay as written).

**Edit 2: SPEC.md, insertions only.**
- **Status list:** `- [ ] Operator replication in AED, **v3 column** (step 2, one session — see "Operator session brief")`.
- **Revision history:** `- **v3.1 (2026-10-13, §9 item 411):** operator session brief added; no row, readout, fixture value or definition changed.`
- **New subsection `### Operator session brief (step 2 — one AED session)`.** Each bullet points to the SPEC line it executes. Only C, R_p, the port pair and the order / frequency are written out.
  - **Project.** Open the operator's existing F-human project and copy the step-1 design to a new design. Never overwrite a step-1 design.
  - **Conductor.** Copper: Finite Conductivity σ = 5.8e7 S/m, Solve Inside off (`SPEC.md:98–99`, `:266–268`).
  - **Gaps.** All 32 gaps carry C = 118.2437 pF in parallel with R_p = 2 212 Ω (`SPEC.md:268–271`). Enter C with all its digits and its unit, `118.2437pF` (step 1's Route A lesson, `:128`).
    - The 30 non-feed gaps are Lumped RLC sheets.
    - The two feed gaps (P17 at 0°, P21 at 90°) are 50 Ω lumped ports. **The same R_p ‖ C must be present across each feed gap in the solved field.** Rows 3–5 need it, because the feed capacitor's dissipation is inside P_acc (`:285–287`). A capacitor added "in post", as step 1's Route A did (`:127–130`), does not change the field.
    - The realisation is the operator's choice; state it in the results file.
  - **Drive.** Equal magnitudes, P21 lagging P17 by 90° (`ccw` = (1, −1j), `:272–273`). Row 6 is the sense check (`:319–320`).
  - **Solve.** 64 MHz only, First Order. Adaptive as in "Frequency and mesh" (`:134–137`). **No re-tune**; say so in the results file (`:275–277`).
  - **Export.** Rows 1–6 with all printed digits (`:279–304`), per the accepted-power note (`:282–295`), plus the solve metadata (`:165`). Results go only to the gitignored `aed_results/`. No AED number ever enters a tracked file.

**Anchors (asserted).** They live in a new `tests/unit/test_ans8_v3_definition.py`, which loads the script by `importlib`, as `tests/unit/test_ans8_v3_rows.py:40–45` does. Every comparison runs on tmp copies.
- **(R) Byte-for-byte regeneration.** `comparison_text(json.loads(committed metrics.json))` equals the committed `COMPARISON.md`, bytes for bytes.
- **(P) Unchanged prefix.** Everything up to the `## Step 2 (SPEC v3)` heading equals the pinned pre-change writer's render up to that heading, byte for byte. The edit is confined to the v3 section.
- **(D) The definition landed once.** The v3 section contains the definition paragraph exactly once: the substrings "accepted power at the 50 Ω lumped port", "Cauchy–Schwarz" and the rendered `rel_line_cut` ("0.83 %"). It no longer contains "incident − reflected at the two 50 Ω feed lines". Every v3 table row still ends `| |`, so the AED cells are empty.
- **(Rp) The operator subsection's arithmetic.** Take ω = 2π·64e6, C = `C_IN_F_HUMAN_PF`·1e-12 (`tests/validation/test_wf10_f_human_tuned_b1_sar.py:68`) and r = `ESR_OHM` (`tests/validation/test_wf10_f_human_copper_esr_root.py:53`), all imported. Then:
  - round(1/(ω²C²r)) == 2212 and round(1/(ωCr)) == 105, the values `SPEC.md:270–271` prints;
  - the new subsection contains "118.2437" and "2 212".
  - **Printed only:** the series-vs-parallel admittance factor 1/(1 + Q_C⁻²) − 1. The drafter's unlogged arithmetic gives −9.0e-5, ≈ 0.11 of the 0.10 pF root grid, so it is inside row 1's resolution.
- **(R0) The test touches no repo file.** The repo's `COMPARISON.md`, `metrics.json` and `SPEC.md` bytes are unchanged by the test.

**Negative controls (rule (e); all asserted, all by construction).**
- **(N-digit).** Make a tmp copy of the committed `COMPARISON.md` and change one digit inside the v3 section (row 5's `0.268112` → `0.268113`). The copy must differ from the render. Assert the inequality and print the first differing line.
- **(N-old).** Copy the pre-change script to `logs/ans8_writer_8f50751e/<the script's repo-relative path>` with a host `git show 8f50751e:<script>` before W1. `8f50751e` is the script's last commit at drafting; the slot re-reads `git log -1 --format=%h -- <script>` before its edit, uses that value and types it into the test as a literal, never `HEAD:`. The mirrored depth is item 400's deviation (`parents[3]`). That writer's render must differ from the committed file, and its v3 section must lack "Cauchy–Schwarz".
- **(N-Rp).** Setting r = 0.02 Ω (one digit changed) gives round(R_p) ≠ 2212.

**Windows.** Foreground, bare, repo-relative, Bash `timeout: 660000`, `-s`, `FEM_EM_ANS_SKIP_PRIVATE=1` in every `ANS-8` window. The steps run in this order:
1. Pin copy (host).
2. Edits.
3. **W-regen.** Use the Write tool to create the gitignored helper `logs/ans8_regen_comparison.py`. It loads the script by `importlib`, reads the committed `metrics.json`, and writes **only** `COMPARISON.md`, from `comparison_text`. Then run:
   `scripts/testing/run_and_log.sh ANS-8-item411-regen "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 logs/ans8_regen_comparison.py'"`
   Afterwards, host `git diff --stat` on the case directory must show `COMPARISON.md` changed and `metrics.json` untouched.
4. **W1.** The script changed, so its unit modules re-run in the same window (rule (c)):
   `scripts/testing/run_and_log.sh ANS-8-item411-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 180 python3 -m pytest tests/environment tests/unit/test_ring_mode_rows.py tests/unit/test_ans8_private_comparison.py tests/unit/test_ans8_v3_rows.py tests/unit/test_ans8_v3_definition.py -s -v --tb=short'"`
   This is item 386's W0 shape (33 passed there) plus the new module.
5. **Wdoc.** Gate on dead = 0 and exit ≠ 1; exit 2 on the stale report is the standing state (`20261012T170143Z_ANS-7-item400-docrefs.log:130`):
   `scripts/testing/run_and_log.sh ANS-8-item411-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
6. Commit: script, new test, `COMPARISON.md`, SPEC.md and the logs.

**Traps already paid for.**
- **The synthetic `rows_v3` in `test_v3_comparison_text`** (`tests/unit/test_ans8_v3_rows.py:91–103`) has no `gap` key and no `power_w["p_acc_wave"]`. The new paragraph must render with `.get` guards, or V3 raises `KeyError` in W1.
- **Never run stage R, S or V3.** Stage R rewrites `metrics.json`; on 10-11 it moved the v1 / v2 floats by up to 1.4e-10 (item 386, D1 reading).
- **`COMPARISON.md` and `metrics.json` are owned by `nobody`** (written by the container). Only a container window rewrites them, and never a host hand-edit.
- **The private side is never touched.** The private writer is never called. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`.
- **No AED number goes into SPEC.md** and no ours-vs-AED level beyond the readout levels the SPEC already carries. No Palace figure is ever set against an AED one.
- **SPEC.md:319's row-6 readout says "ours 21.5".** That is item 358's centre-sample purity (item 369's brief, above). Our row 6 at O reads 20.7092 (`COMPARISON.md:94`). **Do not edit the readout.** Return it as a flag for the review / weekly; the SPEC is the weekly's.
- **SPEC.md:270's R_p ‖ C is the high-Q equivalent.** The exact parallel C is lower by 1/Q_C² ≈ 9e-5. Printed by (Rp) and never changed.
- **Pinned copies:** the mirrored path, never `HEAD:` (`OPS-43`).
- **The pre-commit leak hook:** a refusal is a stop. Never use `--no-verify`.
- **Formatting:** `≈`, never a tilde. No `-k`, no pipe, never delete a log, no edit after W1 (rule (i)).

**Scope.** Text only: one writer paragraph, one relabel, one SPEC subsection.

**Not claimed:**
- any AED value or step-2 verdict;
- the degree-2 v3 column (`xl` entry 35);
- any change to a SPEC row, readout, definition or fixture value;
- the dashboard edit (the review's);
- the retirement of known-issues 2026-10-11. It retires when entry 35 prints the excess at degree 2.

**Negative result.** Report and stop. Park on `attempt/ANS-8-item411-*` and mark item 411 BLOCKED in the same commit (rule (d)).
- **(P) red:** the edit touched a v1 / v2 line. This is a coding error: fix it in-slot, then re-run W-regen and W1.
- **(R) red after W-regen:** the helper did not use the committed writer, or `metrics.json` moved. Restore nothing by hand, and park.
- **A pre-existing `ANS-8` unit module red in W1:** check known-issues first. If it is new, park.
- **A leak-hook refusal:** unstage and park. Report the location only.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Outcome (2026-10-13, 09:30 slot): DONE.** Writer edit as briefed (`_v3_p_acc_definition` + row-5 relabel, nothing else); SPEC.md insertions as briefed. W-regen `20261013T094922Z_ANS-8-item411-regen.log` (metrics.json sha256 unchanged); W1 `20261013T094942Z_ANS-8-item411-W1.log` 35 passed in 82.6 s — (R) 8 250 B equal, (P) 4 794-char prefix equal to the pinned `8f50751e` render, (D) definition once / "0.83 %" once / old clause absent / 18 rows with empty AED cells, (Rp) R_p 2 211.54 Ω → 2 212, Q_C 105.16 → 105, parallel-C factor −9.043e-5 (−0.0107 pF, printed); (N-digit), (N-old), (N-Rp) caught; Wdoc `20261013T095110Z_ANS-8-item411-docrefs.log` dead 0, exit 2 (standing stale report). **Reading:** the regeneration also re-orders the `p17_alone` rows of the two step-1 per-ring k-spectra tables (values unchanged): `_fmt_share_table` iterates the `drives` dict and `metrics.json` is stored with `sort_keys`, so the 10-11 file (rendered from the in-memory dict) was not a pure function of the JSON; the pinned writer gives the same prefix, so (P) holds and the edit is confined to the v3 section. Flags left for the 10-14 weekly unchanged: SPEC.md row-6 "ours 21.5" vs 20.7092; R_p ‖ C is the high-Q form.

**Status it can move.** No glyph. `ANS-8` is ✅ on step 1; its Done-when is met, and step 2 stays open on the row (`PROJECT_PLAN.md:3522`), so rule (j) does not bind.
- On DONE, §10 H22 gains a DONE line (the review marks it).
- The dashboard's W2 moves to the top of Waiting-on-you, with the slot's returned text (the review's edit, flag 4).
- The SPEC Status list gains the unticked operator line.
- Known-issues 2026-10-11 is unchanged.

**Spawn line:** "Use the implementer agent to execute `ANS-8` §9 item 411 per `docs/planning/chunks/ANS-8.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`, `FEM_EM_ANS_SKIP_PRIVATE=1` in every window. No stage of `ans:8` runs; never hand-edit `COMPARISON.md`; never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`; do not edit `docs/status/dashboard.md` — return the W2 text."


## §7 state line replaced 2026-10-13 03:00 review (closure on step 1; previous text verbatim)

🟡 The runnable half landed 2026-10-01 against the v1 rows and the operator's v1 AED column is in but not adjudicated. SPEC v2 makes the own-`C_mode` row 2b primary; v2 rows (2b / 2c / 4a) landed 2026-10-02 (§9 item 161, `ring_mode_rows_v2`, all green, `20261002T050705Z_ANS-8-v2-W1-R.log`), and the copper and First Order AED columns are open. Step 1c (§10 H15c, §9 item 366, 2026-10-11, no solve) ran the v2 rows on the saved degree-2 σ = 800 32×32 (κ off): the loaded matrix reproduces its window (σ_max 4.297e-10, 18 class means ≤ 2.473e-10, `REUSE_REPRODUCTION_RTOL`), and at matched κ degree 1 → 2 moves `C_mode` −4.314 % (one-point f₀ −2.181 %) and the 2b ccw feed-ring share 0.754987 → 0.713591 (k = 0 leak 0.201230 → 0.236516) — records for the 10-12 weekly's private read, no tuning or h claim (`20261011T123443Z_ANS-8-step1c-W1.log`). **Step 2 runnable half (v3 rows) landed 2026-10-11** (§10 H21, §9 item 386, attempt 2 of 4): SPEC v3 rows 1–6 on H9′'s tuned copper F-human at degree 1, per W of P_acc,SPEC (the 50 Ω line cut) with the field-partition figure beside — feed-ring k = ±1 share 0.993845, \|B₁⁺\|(O) 1.668096e-6 T/√W, SAR10g(c\*) 0.1587143 W/kg/W, phantom / copper / ESR 0.268112 / 0.019926 / 0.703644 (sum 0.991683, residual = the feed-sheet gap, measured as the sheets' Cauchy–Schwarz excess, CS_s / X_s = 1.000000000), \|B₁⁺\| / \|B₁⁻\| 20.7092; AED column blank (`20261011T170627Z_ANS-8-v3b-V3.log`, `20261011T171000Z_ANS-8-v3b-R.log`). **Done-when:** ✅ only after the AED half lands and the weekly adjudicates share and leak at a stated relative level. *History: `docs/planning/chunks/ANS-8.md`.*

### Closure audit and review ruling 2026-10-13 03:00 (daily)

`auditor` PASS on the 2026-10-12 weekly's entry 29 / H15c ruling (our side): v2 rows reproduce on today's code (`20261011T171000Z_ANS-8-v3b-R.log:41–80`), the saved degree-2 σ = 800 32×32 reproduces its window to 4.297e-10 (`20261011T123443Z_ANS-8-step1c-W1.log:94–95`), asserts at `tests/validation/test_ans8_degree2_rows.py:256–296, :370`. The row text was checked against the private file's verdict lines (verdict words only). Caveat for the operator: the Done-when's "share and leak at a stated relative level" rests on the 10-03 adjudication (mode content within 0.6 %), whose comparand order is not stated in tracked text; the headline degree-2 `C_mode` / `f₀` are printed readings, gated by the reproduction identities.

**On the item 411 draft:** queued as drafted. Two SPEC readings the drafter found go to the 10-14 weekly, unchanged by the slot: `SPEC.md:319` quotes "ours 21.5" (item 358's centre-sample purity) where row 6 at O reads 20.7092; and `SPEC.md:270`'s R_p ‖ C is the high-Q form (the exact series→parallel map lowers C by 1/Q_C² ≈ 9e-5, inside row 1's 0.10 pF resolution). The operator brief states that the feed gaps' R_p ‖ C must be present in the solved field (rows 3–5) and leaves the method to the operator.

### §9 item 513 — 2026-10-18 (§10 H29, `ANS-8` SPEC v4; lane 1, no solve)

**Outcome 2026-10-18 09:36Z: DONE.** SPEC v4 landed (PRE_V4_SHA `e5d287d4`). (Rp4) R_p 2 317.016850 Ω high-Q / 2 317.216850 Ω exact (difference r to 1e-12), Q_C 107.634 → 108, C_p/C − 1 = −8.6310e-5 (printed), R_p ratio 1.047695056 = (C_d1/C_d2)²; (B4) / (F4) / (Rec) / (N-pin) / (N-Rp4) / (R0′) and item 411's (R) / (P) / (D) / (Rp) / (R0) green — 19 passed in 92.8 s, 94 s harness (`20261018T093350Z_ANS-8-H29-W1.log`); docrefs `dead=0 … exit=2` (`20261018T093527Z_ANS-8-H29-docrefs.log`). Two deviations from the brief's text, both disclosed: (i) the brief's new "Project" bullet would have put "118.2437" inside the span (B4) asserts free of it, so it reads "built at the v3 capacitor (\"Degree-1 record\" below)"; (ii) the brief's "Our side" filename `h9p-node-115.5209-d2.json` in a code span tripped docrefs as dead (`20261018T093327Z_ANS-8-H29-docrefs.log:37`, the record is not yet produced), so it is named in words ("entry 41's degree-2 node record at 115.5209 pF"). Superseded windows kept: `20261018T093146Z_ANS-8-H29-W1.log` (green, pre-wording-fix), `20261018T093340Z_ANS-8-H29-docrefs.log` (dead=0, before a parenthetical correction). The exact C_p 115.5109 pF is printed in the record subsection (ruling Q4). The dashboard W2 text is returned in the slot report (the review's edit).

*Written by the 2026-10-18 03:00 daily review from a read-only drafter's report. The review checked every command against its source before queueing.*

**Licence.** §10 chain H, "H29 (lane 1, no solve). `ans:8` SPEC v4 — step 2's fixture capacitor moves from 118.2437 pF to the degree-2 root C_in(d2) = 115.5209 pF … R_p = 1/(ω²C²r) recomputed at 64 MHz and r = 0.2 Ω (≈ 2 317 Ω), the operator session brief updated, the degree-1 column kept as a labelled record; `docrefs` 0 dead" (`PROJECT_PLAN.md:4783–4790`). The 10-17 weekly's benchmarks paragraph re-poses step 2 at the root before the operator's session (`:5004–5006`) and puts H29 first in the queue (`:5024`). The root is H26's: P17 (1, +1j) secant 115.520946 pF in the bracket [115.5000, 116.0000] pF, resolution 0.5000 pF (`20261017T070008Z_WF-10-step1e-H26.log:8664`, `:8674`), landed to 4 decimals as `C_IN_F_HUMAN_D2_PF = 115.5209` (`tests/validation/test_wf10_f_human_tuned_b1_sar.py:78`, §9 item 511). The 10-12 epitaph binds: a benchmark column is adjudicated at the production convention only (`docs/planning/plan-archive.md:41797–41801`). **This is a fixture-value edit the weekly ordered. No row, readout, band or definition changes.** The operator's session (dashboard W2) waits on it.

**`ANS-8` SPEC v4: step 2's capacitor at C_in(d2) = 115.5209 pF, R_p ‖ C and Q_C recomputed, the operator brief re-pointed, the degree-1 values kept as a labelled record. One additive unit test proves the arithmetic and the text.**
- **Executor:** implementer, or the slot directly (item 411's precedent for a docs-and-test change).
- **Tier:**
  - W1 standard (`-n 1`, complex build, because the test loads the `ans:8` script, which imports dolfinx). Item 411's W1 ran this module with four others: 35 passed in 82.6 s, 84 s harness (`20261013T094942Z_ANS-8-item411-W1.log:185–189`).
  - Wdoc smoke (1 s; `20261013T095110Z_ANS-8-item411-docrefs.log:134`).
- **Class:** lane-1 revalidation (median 3.5). **Cost:** ≈ 3.5 measured / ≈ 25 predicted slot-minutes: 15 fixed, ≈ 6 for the SPEC and the test, ≈ 2 for windows, the rest commit and journal.
- **Depends on:**
  - nothing in §9;
  - on disk: the pinned writer copy `logs/ans8_writer_8f50751e/examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py` (present, 69 770 B, 2026-10-13). `test_comparison_definition`'s (N-old) fails without it (`tests/unit/test_ans8_v3_definition.py:61–65`).
  - If it is absent, restore it with a host `git show 8f50751e:examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py`, redirected to that path (item 411's step 1). Nothing else is read from `logs/`.
- **Family:** `ANS-8` step 2 at the root; attempt 1 of 4.

**The arithmetic (drafter's, unlogged; the test re-computes it).** With ω = 2π·64e6, r = 0.2 Ω and C = 115.5209 pF:
- Q_C = 1/(ωCr) = 107.634, so it rounds to **108** (v3: 105.156 → 105).
- The step's form gives R_p = 1/(ω²C²r) = 2 317.017 Ω. The ESR ruling's exact form gives R_p = r(1 + Q_C²) = 2 317.217 Ω (`PROJECT_PLAN.md:3856–3858`). The two forms differ by exactly r = 0.2 Ω (8.6e-5 relative), and both round to **2 317 Ω** (v3: 2 211.54 / 2 211.74 → 2 212). The step's "≈ 2 317 Ω" is correct.
- Exact parallel C: C_p = C·Q²/(1 + Q²) = 115.51093 pF, a factor of −8.631e-5 (−0.0100 pF). That is 0.020 of the 0.50 pF degree-2 bracket. The high-Q print "‖ 115.52 pF" is v3's convention; v3 printed "‖ 118.24 pF" where the exact value was 118.2330 pF.
- Scaling check: R_p(d2) / R_p(d1) = (118.2437 / 115.5209)² = 1.047695 in the high-Q form, exact as a closed form.

**What changes in `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/SPEC.md` (line numbers at `e5d287d4`, the file's last commit; re-read them before editing).**
- **Status (`:22–36`).** Insert after `:30`: `- [x] **SPEC v4 (2026-10-18, §9 item 513, §10 H29)** — step 2's capacitor at the degree-2 root, C = 115.5209 pF; the degree-1 values are a record.` Re-word `:35`'s bold label from "**v3 column**" to "**v4 column** (C = 115.5209 pF)". The line stays unticked.
- **Step 2 precondition (`:260–265`).** Unchanged: it is v3's history (118.2437 pF, resolution 0.10 pF). Insert one paragraph after it, `**v4 (2026-10-18, §10 H29): the degree-2 root.**`, stating:
  - at degree 2 (the production order) the degree-1 node is detuned;
  - `WF-10`'s H26 found the in-model root at **C_in(d2) = 115.5209 pF**, P17 (1, +1j) secant, bracket [115.50, 116.00] pF, **resolution 0.50 pF**, four-reading spread 0.036 pF (`20261017T070008Z_WF-10-step1e-H26.log:8674`);
  - the step-2 fixture capacitor is that value, because a benchmark column is adjudicated at the production convention only (the 10-12 ruling).
- **Fixture (`:267–278`).** In `:270`, "C = 118.2437 pF" becomes "C = 115.5209 pF". In `:271–272`, "**2 212 Ω ‖ 118.24 pF**, Q_C = 1/(ωCr) = 105" becomes "**2 317 Ω ‖ 115.52 pF**, Q_C = 1/(ωCr) = 108". Nothing else in the paragraph changes.
- **Accepted power (`:293–294`).** "at C = 118.2437 pF" becomes "at C = 115.5209 pF". The definition text and the "≈ 0.8 % … at degree 1" clause stay as written: that clause describes the degree-1 record, and H30 prints the degree-2 excess.
- **Row 1 (`:300`).** "with C = 118.2437 pF" becomes "with C = 115.5209 pF".
- **Readout (`:307–321`).** Untouched (see ruling questions Q2 and Q3 below).
- **Our side (`:323–326`).** Replace the paragraph with three statements:
  - our adjudicated column is `ans:8`'s v3 rows at **degree 2 at C_in(d2)**, from `xl` entry 41's node record `h9p-node-115.5209-d2.json` (§10 H30);
  - the degree-1 rows at 118.2437 pF (§10 H21, `COMPARISON.md`'s step-2 section, which already labels itself "degree 1 … C_in = 118.2437 pF", `COMPARISON.md:73`) are **a labelled record, not the adjudicated column**;
  - "Lossless step-1 columns stay as they are, labelled lossless."
- **Operator session brief (`:328–355`).**
  - In `:336`, the C and R_p values become **C = 115.5209 pF in parallel with R_p = 2 317 Ω**.
  - In `:338`, the entry string becomes `` `115.5209pF` ``.
  - Add one bullet after "Project": "If a design was already built at 118.2437 pF, copy the step-1 design again; never re-use it for v4."
  - Nothing else in the brief changes.
- **New subsection, between the brief and `## Revision history` (`:357`):** `### Degree-1 record (v3 / v3.1 values — do not enter)`. It holds one paragraph:
  - v3's capacitor C = 118.2437 pF (resolution 0.10 pF), R_p = 2 212 Ω, Q_C = 105;
  - our degree-1 column at that node (`COMPARISON.md` step 2) is kept as a record;
  - no AED column is commissioned at it.
- **Revision history.** Insert after `:362`: `- **v4 (2026-10-18, §9 item 513, §10 H29):** step 2's capacitor moves to the degree-2 root C_in(d2) = 115.5209 pF (R_p ‖ C = 2 317 Ω ‖ 115.52 pF, Q_C = 108, high-Q form, 64 MHz); the operator brief re-pointed; the degree-1 values kept as a labelled record; no row, readout or definition changed.`
- **Not edited:** `COMPARISON.md` (test (R) requires it byte-identical), `metrics.json`, the script, the guide, `aed/`, and the dashboard (the review's job: the slot returns the W2 text).

**Code: append only to `tests/unit/test_ans8_v3_definition.py`.** The two existing test functions and the imports at `:34–35` stay byte-identical. **The trap this avoids:** `test_operator_brief_arithmetic` asserts that `"118.2437"` and `"2 212"` appear in `spec[BRIEF_HEAD : "## Revision history"]` (`:152`, `:156`). Placing the degree-1 record subsection inside that span keeps the test green without editing it. That is disclosed in the new test's docstring, and the review rules on it (Q1).
- **New import:** `C_IN_F_HUMAN_D2_PF` from `tests.validation.test_wf10_f_human_tuned_b1_sar` (`:78`), beside the existing two.
- **New constants:**
  - `PRE_V4_SHA = "<the sha>"`: re-read with `git log -1 --format=%h -- examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/SPEC.md` **before** the edit (`e5d287d4` at drafting) and typed in as a literal, never `HEAD:`;
  - `SPEC_REL`, the SPEC's repo-relative path;
  - `RECORD_HEAD = "### Degree-1 record (v3 / v3.1 values — do not enter)"`.
- **New `_brief(text)`:** returns `text[text.index(BRIEF_HEAD):text.index(RECORD_HEAD)]`. That is the brief alone, without the record.
- **New test `test_operator_brief_arithmetic_v4`:**
  - **(Rp4)** `_rp(C_IN_F_HUMAN_D2_PF·1e-12, ESR_OHM)` (`:136–138`) gives round(R_p) == 2317 and round(Q_C) == 108. The exact form r(1 + Q_C²) also rounds to 2317, and it differs from the high-Q form by exactly r to ≤ 1e-9 relative. C_p / C − 1 is printed (predicted −8.631e-5). The ratio R_p(d2) / R_p(d1) equals (C_d1 / C_d2)² to ≤ 1e-12.
  - **(B4)** `_brief(spec)` contains `"115.5209"`, `` "`115.5209pF`" `` and `"2 317"`. It contains neither `"118.2437"` nor `"2 212"`.
  - **(F4)** The fixture paragraph (from `"**Fixture (both codes):**"` to `"**Quantities to export"`) contains `"C = 115.5209 pF"`, `"2 317 Ω ‖ 115.52 pF"` and `"= 108"`, and does not contain `"C = 118.2437 pF"`. Row 1's table line contains `"C = 115.5209 pF"`.
  - **(Rec)** `spec.count(RECORD_HEAD) == 1`. The record span (`RECORD_HEAD` to `"## Revision history"`) contains `"118.2437"` and `"2 212"`. Exactly one revision-history line starts with `"- **v4 "`.
  - **(R0′)** The repo's `SPEC.md` bytes are unchanged by the test.

**Anchors (asserted).**
- (Rp4): the closed forms, with C and r imported, never retyped.
- (B4), (F4) and (Rec): the text, on the committed SPEC.
- Unchanged and re-run green: item 411's (R), (P), (D), (Rp) and (R0) (`test_ans8_v3_definition.py:73–161`). (R) proves `COMPARISON.md` was not touched.
- **(Doc)** The docrefs checker's `RESULT:` line reads `dead=0` with `exit=2` (the standing stale report; `RESULT: dead=0 guide=0 stale=133 stale_severity=report exit=2` at `20261017T172100Z_EX-128-post-docrefs.log:173`). An exit of 1 is red.

**Negative controls (rule (e)).** All are asserted, all by construction.
- **(N-pin)** The pre-v4 SPEC, read by `subprocess.run(["git", "-c", f"safe.directory={REPO}", "-C", str(REPO), "show", f"{PRE_V4_SHA}:{SPEC_REL}"])` (the `PORT-24` pattern, `tests/unit/test_port24_esr_circuit.py:70–73`), fails (B4): its brief span contains "118.2437" and lacks "115.5209". It also has no `RECORD_HEAD`, so `_brief` is taken to `"## Revision history"`. Assert both, and print the first offending substring.
- **(N-Rp4)** r = 0.02 Ω gives round(R_p) ≠ 2317. C = `C_IN_F_HUMAN_PF` (the v3-blind edit) gives round(R_p) == 2212 ≠ 2317. The separation is 2317 / 2212 = 1.047, the closed form (C_d1 / C_d2)².
- None is a predicted factor.

**Windows.** Each runs in the foreground with Bash `timeout: 660000`, with `-s`, and with `FEM_EM_ANS_SKIP_PRIVATE=1` in the `ANS-8` window. **The harness call is the whole Bash command.** The order:
1. Host `git log -1 --format=%h -- <SPEC>`, which sets `PRE_V4_SHA`. Then `ls -la` the pinned writer copy.
2. Edit the SPEC, then the test.
3. **W1** (standard, `-n 1`; predicted 80–90 s; item 411's W1 path list minus `test_ring_mode_rows.py` and `test_ans8_private_comparison.py`, which this item does not touch):
   `scripts/testing/run_and_log.sh ANS-8-H29-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 180 python3 -m pytest tests/environment tests/unit/test_ans8_v3_rows.py tests/unit/test_ans8_v3_definition.py -s -v --tb=short'"`
4. **Wdoc** (item 411's command, label changed):
   `scripts/testing/run_and_log.sh ANS-8-H29-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
- No edit to the SPEC or the test after W1 starts (rule (i)). Any edit means W1 and Wdoc again.

**Traps already paid for.**
- **Item 411's (Rp) substring assert** (`:156`) reads the whole span to `## Revision history`. Moving the v3 values out of that span turns a green test red, and editing the test to rescue it is the review's to rule (Q1).
- **(R) is byte-for-byte on `COMPARISON.md`** (`:85`). Never re-render it, and never run stage R, S or V3: stage R moved v1 / v2 floats by up to 1.4e-10 on 10-11 (item 386, D1 reading).
- **The pinned copy:** a pinned sha, never `HEAD:`. A control copied from `HEAD:` goes stale the moment the change lands (`OPS-43`).
- **Step 2 is a fixture value, not a definition.** The 10-12 weekly's accepted-power paragraph (`:283–296`) stays verbatim. Its "≈ 0.8 % … at degree 1" is a degree-1 reading, and H30 prints the degree-2 one.
- **The leak hook:** no AED number and no ours-vs-AED level enters the SPEC. A refusal is a stop; never `--no-verify`.
- **Formatting:** `≈`, never a tilde (CLAUDE.md, 2026-09-28 sweep). Thin-space thousands ("2 317") match v3's "2 212", which the tests look for.
- **`COMPARISON.md` and `metrics.json` are `nobody`-owned**, and neither is edited here.
- No `-k`, no pipe, never delete a log.

**Scope.** The SPEC's step-2 capacitor and its derived values, the operator brief, a labelled degree-1 record, and one additive unit test. **Not claimed:**
- any AED value, verdict or level;
- the degree-2 v3 rows (H30, serial on entry 41);
- a readout change (Q2, Q3);
- a change to the accepted-power definition;
- the dashboard edit (the review's);
- the exact-vs-high-Q C_p choice (Q4).

**Negative result.** Report and stop. Park on `attempt/ANS-8-H29-*` and mark item 513 BLOCKED in the same commit (rule (d)).
- **(Rp4) red:** the imported constant is not 115.5209, or the formula is wrong. Never retype the value.
- **(B4), (F4) or (Rec) red:** an edit was missed. Fix it in-slot before W1 re-runs (the test exists to catch exactly this).
- **(N-pin) red:** the pinned sha already carries v4. Re-read the sha.
- **Item 411's (R) / (P) red:** `COMPARISON.md` or the writer moved. Check `git status`; restore nothing by hand; park.
- **Wdoc exit 1:** a dead reference. Report the location.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** No glyph. `ANS-8` is ✅ on step 1, and step 2 stays open on the row (`PROJECT_PLAN.md:3558`), so rule (j) does not bind.
- §10 H29 gains a DONE line.
- The `ANS-8` row's step-2 sentence gains "SPEC v4 (2026-10-18, §9 item 513): C = 115.5209 pF, R_p ‖ C = 2 317 Ω ‖ 115.52 pF; degree-1 values a labelled record".
- **Dashboard W2 flips from "hold" to "ready".** That is the review's edit, made from the slot's returned W2 text (≤ 4 lines: C = `115.5209pF`, R_p = 2 317 Ω on all 32 gaps and across both feed gaps in the solved field, First Order, 64 MHz, no re-tune).

**Execution order.** sha re-read → `ls` → SPEC edit → test append → W1 → Wdoc → commit (SPEC.md, the test, both harness logs, §10 H29 DONE, the §7 sentence) → journal → return the W2 text.

**For the review to rule on before queueing.**
- **(Q1)** The geography fix keeps item 411's test unedited by leaving "118.2437 / 2 212" in its span, inside a "do not enter" record subsection. The alternative edits `:156` to read the record span instead. That edit points the assert at the text the SPEC now holds, and adds the (B4) absence check, so it is not a loosening. The drafter recommends the geography fix (additive), disclosed.
- **(Q2)** Row 1's readout "within one grid resolution of zero" was written against degree 1's 0.10 pF grid. At degree 2 the bracket is 0.50 pF (4.3e-3 relative). The v4 paragraph states the 0.50 pF, which widens what the readout tolerates fivefold. Either the review accepts that, or the weekly re-poses the readout (for example, AED's own Im Y_in slope converted to pF, as the readout's second clause already does).
- **(Q3)** `SPEC.md:320` "(ours 21.5)" is still item 358's centre-sample purity at the degree-1 node (flagged 10-13, never ruled). At v4 it is a stale fixture value inside a readout line. The drafter leaves it untouched.
- **(Q4)** The high-Q print (C entered as 115.5209 pF) against the exact C_p = 115.5109 pF. The 0.0100 pF difference is 0.02 of the bracket. Flagged 10-13 and never ruled; v4 keeps v3's convention unless the review says otherwise.

**Spawn line (if delegated):** "Use the implementer agent to execute `ANS-8` §9 item 513 (H29) per `docs/planning/chunks/ANS-8.md` (read the brief and any review ruling beneath it first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W1 then Wdoc, in the foreground, with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`, and `FEM_EM_ANS_SKIP_PRIVATE=1`. Edit only `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/SPEC.md` (the lines the brief names) and append to `tests/unit/test_ans8_v3_definition.py`, leaving its existing functions byte-identical. No stage of `ans:8` runs. Never edit `COMPARISON.md`, `metrics.json`, the script, the guide or `docs/status/dashboard.md`; return the W2 text. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. You are the executor: do not spawn agents, and never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ANS-8-H29-*`, mark item 513 BLOCKED in the same commit, and report."

### §9 item 516 — 2026-10-18 (§10 H25, `ans:8` v3 rows 1 / 2 at degree 2; chain no-solve)

**Outcome 2026-10-18 10:00Z: DONE, attempt 1 of 4.** New module `tests/validation/test_ans8_copper_degree2_rows.py` only. Host `ls`: all five npz present. W0 `20261018T095551Z_ANS-8-H25-W0.log` 5 passed / 1 skipped (shuffled-D16 names = slot synthesis bitwise, 18 classes; (Σk) 0; (S) synthetic 1.563e-13). W1 `20261018T095559Z_ANS-8-H25-W1.log` 14 passed, 96 s harness (test body 0.19 s). **Asserted, green:** (F) degree 2, κ off, 369 387 cells, 2 570 024 unknowns, 16 ranks, P17 first, names bijection, set(order) = a2; (F1) both `h12c` files, 0 slot-map mismatches; (Syn) bitwise; (L-d1) class means worst 2.176e-10, C_mode 118.118580 pF (rel 1.776e-15 vs npz key, 3.515e-13 vs record), off-block 4.195925e-3, 2b shares/leaks worst 4.809e-7; (G1′) reciprocity 9.964981e-14, σ_max 0.999983227022; **(L-d2) σ_max rel 1.010e-13, 18 class means worst 3.416e-10**; S17,17 0.402723+0.391780j = record; (L-σ) 3.053e-10 / 3.423e-10; NC P17 × 1.01: σ800 d1 314.7×, Cu d2 314.7×, Cu d1 314.8× (bar 10×; Cu factors inside `PRED_D2_COLUMN_CONTROL`); (B0) bitwise; (Σk) 1.110e-16, bitwise; **(S) ESR shift on Cu d2: dRe 4.477e-13, dIm 2.522e-15, dC_mode 4.441e-16** (≤ 1e-12); (P) planted 0.4 Ω on P17 misses by 3.125e-2 (predicted ≈ 3e-2). **Record (circuit route, Norton feeds, not the in-model row):** Cu d2 | Cu d1 (a2) | σ800 d2 | σ800 d1 — C_mode 114.5221 | 118.1186 | 108.2288 | 113.1081 pF; C_k2 76.383 | 78.817 | 71.758 | 75.454 pF; off-block 3.909696e-3 | 4.195925e-3 | 3.909947e-3 | 4.182881e-3; anti 1.000000 all; 2b ccw feed share 0.968088 | 0.992768 | 0.713591 | 0.754987, leak 0.025825 | 0.005671 | 0.236516 | 0.201230; with r = 0.2 Ω ccw 0.948363 / leak 0.041953 (C_mode unchanged to all 13 digits). **H25:** C_mode(d2, Cu) = 114.5221 pF beside 118.1186 / 115.5209 pF, (C_mode − root)/root = −8.646e-3 — inside `PRED_CMODE_CU_D2_PF_WINDOW` (HIT); point predictions 115.399 (−0.877 pF) / 113.023 (+1.499 pF); copper order move −3.045 % (σ800 −4.314 %). Predictions missed (readings, printed): ccw feed share ≥ 0.99 (0.968) and leak in [0.002, 0.02] (0.0258). C₂ control share 0.032015 (< 0.5); 2c leak max 4.127e-7 (≤ 5e-5).

*Written by the 2026-10-18 03:00 daily review from a read-only drafter's report. The review checked every command against its source before queueing.*

**Licence.** §10 chain H, "H25 (lane 1 → heavy, no solve). `ans:8` v3 rows 1 / 2 at degree 2 on copper from entry 36's saved matrix (`ring_mode_rows_v2`) beside the degree-1 and σ = 800 rows; target: the k = ±1 share and k = 0 leak printed at their imported identities, `C_mode`(d2, copper) printed beside H12c's degree-1 118.1186 pF and the root 115.5209 pF (predicted within 1 % of the root — printed, never asserted). Serial on entry 36" (`PROJECT_PLAN.md:4777–4782`). It also discharges entry 36's own readout, "`ring_mode_rows_v2` on the saved matrix (`C_mode`, k-share) beside the degree-1 and σ = 800 rows" (`docs/testing/xl-pending.md:1804–1809`), which the window did not print. **No band, tolerance, Done-when or target is set or moved.** Template: item 366 (H15c, `tests/validation/test_ans8_degree2_rows.py`).

**Entry 36 ran (2026-10-18 02:00).** `20261018T070008Z_WF-12-step5-H16c.log`: 15 passed, `[capture] rc=0` (`:9106`), Status 0 (`:9109`). Readings (a `log-pathologist` reads the log separately; the slot uses only the lines cited):
- (L1) save/reload bitwise (`:8032`);
- (ii) reciprocity 9.964981e-14 (`:8035`);
- (iii) σ_max 0.999983227022 (`:8036`);
- 18 classes with d2 / d1 means (`:8038–8055`), worst spread 3.2946e-2 (`:8056`);
- S17,17 0.402723+0.391780j (`:8057`);
- price: build 124.8 s, factorising drive 299.4 s, window 701.2 s (`:8058`).

**The saved matrix is on disk:** `logs/wf12/h16b-cu-s32-d2.npz`, 23 077 B, 2026-10-18 02:13 (container path `/workspace/logs/wf12/h16b-cu-s32-d2.npz`, `ENTRY36_SAVE`, `tests/validation/test_wf12_f_human_copper_s32_degree2.py:56`).

**`ANS-8` H25: SPEC v2 rows 1 / 2b / 2c / 4a on the saved degree-2 copper 32×32, beside degree-1 copper and both σ = 800 matrices, tied to their windows by reproduction, with the ESR uniform-shift identity asserted on the new matrix.**
- **Executor:** implementer (item 366's precedent).
- **Tier:** standard, `-n 1`, complex build (pure NumPy; item 366's environment is kept). W0 is smoke-shaped. Item 366's W1 ran the same algebra on two matrices in 72.2 s, 74 s harness (`20261011T123443Z_ANS-8-step1c-W1.log:221–225`), and its rows phase took 0.09 s (`:156`). **The step text's "→ heavy" is not needed:** no window here exceeds the standard ceiling.
- **Class:** chain no-solve (median 18.6). **Cost:** ≈ 18.6 measured / ≈ 38 predicted slot-minutes: 15 fixed, ≈ 15 code, W0 ≈ 0.5, W1 ≈ 1.3, commit and journal ≈ 6.
- **Depends on:**
  - §10 entry 36 (ran, green); nothing open in §9.
  - On disk, all gitignored and all opened read-only (`np.load(…, allow_pickle=False)`):
    - `logs/wf12/h16b-cu-s32-d2.npz`;
    - `logs/wf12/h12c-cu-s32-a2.npz` (20 812 B; H12c attempt 2, saved at `20261006T170222Z_WF-12-step3-a2-W1.log.gz:6399`);
    - `logs/wf12/h12c-cu-s32.npz` (19 636 B; entry 36's d1 comparand, `test_wf12_f_human_copper_s32_degree2.py:57`, `:380–389`);
    - `logs/wf7/wf7-s32-degree2-n16.npz` and `logs/wf7/wf7-s32-degree1-n8.npz` (item 366's).
  - All five were present at drafting (`ls -la logs/wf12/ logs/wf7/`, 2026-10-18).
- **BLOCKED-if-absent:**
  - `h16b-cu-s32-d2.npz` absent ⇒ stop, and mark item 516 BLOCKED "entry 36's matrix absent; the review re-poses". Never re-solve.
  - Any degree-1 or σ = 800 file absent ⇒ run without that column and say so. The anchors on the absent file are not run, and the item stays open for the review.
- **Family:** `ANS-8` step 2's degree-2 copper record; attempt 1 of 4.

**What the files are (read from the code and the zip members; the slot re-checks).**
- `h16b-cu-s32-d2.npz` was written by `_save_s32` (`test_wf12_f_human_copper_s32_degree2.py:112–132`).
  - Keys: `s32`, `order` (port names, P17 first), the per-column `t_solve` / `t_post` / `held` / `p_*`, and the scalars `degree`, `cells`, `n_dofs`, `n_edges`, `n_faces`, `freq_hz`, `ranks`, `kappa_applied` (`:69–72`).
  - **It has no `slot_k` / `slot_sign`,** so `ans:8`'s `_slot_port_map` (`08_…py:725–728`) cannot run on it. Its port map comes from names only, through `ring_port_map_from_names(order)` (`src/fem_em_solver/ports/ring_modes.py:42`).
  - Load it with `_load_s32` (`:135–139`).
- `h12c-cu-s32-a2.npz` carries `c_mode`, `slot_k`, `slot_sign` and `order`. Its window printed:
  - C_mode_Cu = 1.181185799070e-10 F = 118.1186 pF and off-block 4.195925e-03 (`20261006T170222Z_WF-12-step3-a2-W1.log.gz:6366`);
  - `ring_mode_rows_v2`'s own c_mode_f equal to that to 1.776e-15 (`:6367`);
  - row 2b p17_alone at own C feed ring 0.993012, k = 0 leak 0.005612, top 0.998482 (`:6382`);
  - cw at own C: share 0.993217 / 0.998703, leak 0.005618 / 0.000013 (`:6387`).
  - The constant `C_MODE_CU_PF = 118.1186` is that value to 4 decimals (`tests/validation/test_wf10_f_human_copper_esr_root.py:59`).
- **Both copper matrices are κ off on the 369 387-cell hole mesh** (`HOLE_CELL_RECORD`, `tests/validation/test_wf12_f_human_copper_matrix.py:55`).

**Code. One new module, `tests/validation/test_ans8_copper_degree2_rows.py`.** No `src/` change. No edit to `ans:8`, `metrics.json`, `COMPARISON.md`, `test_ans8_degree2_rows.py`, any `test_wf12_*` or `test_wf10_*` module, or `tests/unit/test_port24_esr_circuit.py`.
- **Knob** `FEM_EM_ANS8_H25` (item 366's `_enabled` / `_need` pattern, `test_ans8_degree2_rows.py:73–84`): unset ⇒ the gated test skips; `"1"` ⇒ runs, and an absent d2 npz **fails**; anything else ⇒ `ValueError`.
- **Loading.** Load the `ans:8` script by `importlib` under a distinct module name (`"ans8_h25"`; the `:87–92` pattern) and use its `rows_v2`, `spec_controls`, `SPEC_DRIVES`, `SPEC_FEED_NAMES`, `SPEC_C_F`, `CONTROL_LEAK_BAR`, `NEGATIVE_CONTROL_CEILING`, `SHARE_ATOL`, `OFF_BLOCK_RTOL` and `_slot_port_map`. For ESR rows, call `ring_mode_rows_v2(…, esr_ohm=ESR_OHM)` directly (`ring_modes.py:233–267`; `rows_v2` has no ESR argument, `08_…py:757–763`).
- **Imported, never restated:**
  - `REUSE_REPRODUCTION_RTOL`, `MATRIX_PASSIVITY_CEILING`, `AZIMUTH_STEP_DEG`, `_azimuth_class` (`tests/validation/test_port_birdcage_ring_matrix.py:171`, `:150`, `:185`, `:206`);
  - `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py:86`);
  - `_reciprocity_ratio`, `CONTROL_COLUMN_SCALE` (`tests/validation/test_port_birdcage_ring_column.py:776`, `:180`);
  - `_class_means`, `SAVE_CONTROL_FACTOR` (`tests/validation/test_wf7_f_human_ring_matrix.py:331`, `:239`);
  - `ENTRY36_SAVE`, `H12C_NPZ`, `UNKNOWNS_RECORD`, `_load_s32` (`test_wf12_f_human_copper_s32_degree2.py:56`, `:57`, `:63`, `:135`);
  - `HOLE_CELL_RECORD`;
  - `ESR_OHM`, `C_MODE_CU_PF` (`test_wf10_f_human_copper_esr_root.py:53`, `:59`);
  - `C_IN_F_HUMAN_D2_PF` (`test_wf10_f_human_tuned_b1_sar.py:78`);
  - `SHIFT_RTOL` (`tests/unit/test_port24_esr_circuit.py:62`);
  - `_summary`, `CLASS_MEAN_D2_RECORD_BOXED`, `SIGMA_MAX_D2_RECORD_BOXED` from `test_ans8_degree2_rows.py` (`:150`, `:44`, `:43`; module-level helpers, never a `test_*` name);
  - item 366's σ = 800 bands as that module imports them.
- **The names-based azimuth synthesis.** For the d2 file: `(j, ring) = ring_port_map_from_names(order)[i]`, `az = j · AZIMUTH_STEP_DEG`, `zc = −1` (bottom) / `+1` (top). Only differences enter `_azimuth_class`, as item 366 showed. It is validated on `h12c-cu-s32.npz`, where both `slot_k` / `slot_sign` and names exist (anchor (Syn)).
- **Record constants** (comparands, not bands; the 361-3 precedent), each with its log:line:
  - `CLASS_MEAN_CU_D2_RECORD` and `CLASS_MEAN_CU_D1_RECORD`: the 18 d2 and d1 means at `20261018T070008Z_WF-12-step5-H16c.log:8038–8055`, 10 digits, in item 366's `(same_ring, steps)` dict shape;
  - `SIGMA_MAX_CU_D2_RECORD = 0.999983227022` (`:8036`);
  - `S1717_CU_D2_RECORD = 0.402723 + 0.391780j` (`:8057`, printed beside);
  - `C_MODE_CU_D1_RECORD_F = 1.181185799070e-10` and `OFF_BLOCK_CU_D1_RECORD = 4.195925e-03` (`…step3-a2-W1.log.gz:6366`);
  - `P17_2B_CU_D1_RECORD = (0.993012, 0.005612, 0.998482)` (`:6382`) and `CW_2B_CU_D1_RECORD = (0.993217, 0.998703, 0.005618, 0.000013)` (`:6387`);
  - `C_MODE_S800_RECORD_F = {2: 1.082287732e-10, 1: 1.131081186e-10}` (`20261011T123443Z_ANS-8-step1c-W1.log:116`).
- **Predictions written first (printed, never asserted):**
  - `PRED_CMODE_CU_D2_PF_WINDOW = (114.37, 116.68)`, within 1 % of the root (the step's);
  - point predictions printed beside it:
    - 115.399 pF (the degree-1 C_mode / root offset 118.1186 / 118.2437 carried to the d2 root);
    - 113.023 pF (σ = 800's order move −4.314 %, item 366, applied to 118.1186; this one falls outside the 1 % window);
  - anti ≥ 0.99;
  - 2b ccw feed-ring share ≥ 0.99, with leak in [0.002, 0.02];
  - `PRED_D2_COLUMN_CONTROL = (2.0e-4, 4.0e-4)` (the ≈ 0.01/32 ceiling);
  - the planted-ESR miss ≈ 3e-2 (`PORT-24`'s (P) order, on σ = 800).

**Anchors (asserted, every band imported).**
- **(F) file identity, d2 copper:**
  - `degree` 2, `kappa_applied` False, `cells == HOLE_CELL_RECORD`, `n_dofs == UNKNOWNS_RECORD` (2 570 024), `freq_hz` 64e6, `ranks` 16;
  - shape (32, 32), `order[0] == "P17"`;
  - `ring_port_map_from_names(order)` is a bijection onto 16 slots × 2 rings;
  - `set(order)` equals the a2 file's `set(order)`.
- **(F1) d1 copper:**
  - both `h12c` files have shape (32, 32);
  - `_slot_port_map(saved) == ring_port_map_from_names(order)`, with 0 mismatches (item 366's (G2′));
  - cells within `HOLE_CELL_RECORD`.
- **(Syn) the names synthesis.** On `h12c-cu-s32.npz`, `_class_means` with the names-synthesised (az, zc) equals `_class_means` with `slot_k · AZIMUTH_STEP_DEG` / `slot_sign` bitwise (both 18 classes). This validates the only synthesis the d2 file allows.
- **(L-d1) degree-1 reproduction:**
  - the 18 class means of `h12c-cu-s32.npz` (names synthesis, read by name) reproduce `CLASS_MEAN_CU_D1_RECORD` to ≤ `REUSE_REPRODUCTION_RTOL`, predicted ≤ 5e-10 (10-digit print rounding);
  - `rows_v2` on `h12c-cu-s32-a2.npz` gives `c_own_f` vs the npz `c_mode` key and vs `C_MODE_CU_D1_RECORD_F`, each ≤ `REUSE_REPRODUCTION_RTOL`. The same comparison read 1.776e-15 at `:6367`.
  - `|c_own_f·1e12 − C_MODE_CU_PF| ≤ 5e-5` pF (the constant's 4-decimal resolution);
  - off-block vs `OFF_BLOCK_CU_D1_RECORD` ≤ `OFF_BLOCK_RTOL`;
  - the p17_alone and cw 2b shares and leaks vs their records ≤ `SHARE_ATOL` (6-digit prints).
- **(G1′) d2:** reciprocity ≤ `RECIPROCITY_BAND`, printed beside 9.964981e-14 (`:8035`), and σ_max ≤ `MATRIX_PASSIVITY_CEILING`.
- **(L-d2) the H25 anchor:**
  - σ_max recomputed vs `SIGMA_MAX_CU_D2_RECORD` ≤ `REUSE_REPRODUCTION_RTOL`;
  - all 18 class means vs `CLASS_MEAN_CU_D2_RECORD`, each ≤ `REUSE_REPRODUCTION_RTOL` (predicted ≤ 1e-9);
  - S17,17 printed beside its record.
- **(L-σ) σ = 800 reproduction.** `rows_v2` on the two `logs/wf7/` files reproduces `C_MODE_S800_RECORD_F` ≤ `REUSE_REPRODUCTION_RTOL` (10-digit prints; item 366 computed them from the same files with the same function).
- **The identities the share and leak are printed at:**
  - **(B0) own-C wiring:** on d2 copper, `own_c_mode`'s `c_f` is bitwise `fixed_c["c_mode_f"]`, and `c_own_f` is bitwise the same (`ring_modes.py:252–264`). This is item 161's B0, the wiring that makes 2b each matrix's own resonance.
  - **(Σk) the k-fractions:** for every drive and control at own C, on each ring, Σ `fractions` = 1 to ≤ 1e-12. `share_k1 == fractions[1] + fractions[15]` and `leak_k0 == fractions[0]` bitwise (`ring_modes.py:100–112`).
  - **(S) the ESR uniform shift, on the d2 copper matrix:**
    - `ring_mode_rows(…, esr_ohm=ESR_OHM)` against `esr_ohm=0` gives |ΔRe λ_k − r| / r ≤ `SHIFT_RTOL` for k = 0…8;
    - |ΔIm λ_k| / |Im λ_k| ≤ `SHIFT_RTOL`, and |ΔC_mode| / C_mode ≤ `SHIFT_RTOL`.
    - This is the ESR ruling's (i) ("uniform r adds r to every mode's eigen-impedance, so `C_mode` is unchanged", `PROJECT_PLAN.md:3853–3854`), which `PORT-24` asserted on σ = 800 at 2.889e-15 / 4.571e-16 (`tests/unit/test_port24_esr_circuit.py:276`).
    - It is what licenses quoting the lossless `C_mode` beside the ESR root.

**Negative controls (rule (e)).**
- **Asserted, backed (same comparison, same file):** P17 × `CONTROL_COLUMN_SCALE` on `wf7-s32-degree1-n8.npz` moves the worst class mean to ≥ `SAVE_CONTROL_FACTOR` × `REUSE_REPRODUCTION_RTOL`. Item 366 read 3.147e-04 = 314.7× (`20261011T123443Z_ANS-8-step1c-W1.log`, its (NC-d1) line); the ceiling is ≈ 0.01/32 ≈ 3.1e-4.
- **Asserted against the same floor, factor predicted** (no prior on these files; item 366's (NC-d2) precedent): P17 × 1.01 on `h16b-cu-s32-d2.npz` against `CLASS_MEAN_CU_D2_RECORD`, and on `h12c-cu-s32.npz` against `CLASS_MEAN_CU_D1_RECORD`. Each ≥ `SAVE_CONTROL_FACTOR` × `REUSE_REPRODUCTION_RTOL`, with the factor printed beside `PRED_D2_COLUMN_CONTROL`.
- **Asserted against the floor, factor predicted:** (P) a planted 0.4 Ω on P17 (the rest 0.2 Ω) misses (S)'s Re shift by ≥ `SAVE_CONTROL_FACTOR` × `SHIFT_RTOL` on at least one k. Predicted ≈ 3e-2, printed. `PORT-24` asserted this on σ = 800, a different fixture.
- **Printed only:**
  - the C₂ control, `rows_v2(…, c_spec_f=<d2 Cu c_k2_f>)` p17_alone feed-ring share, predicted < `NEGATIVE_CONTROL_CEILING`. H12c's degree-1 circuit read 0.105707 (`…a2-W1.log.gz:6397`), but through the Thevenin route, so it is not the same comparison.
  - the 2c controls' k = 0 leak at own C, predicted ≤ `CONTROL_LEAK_BAR`.

**Printed (the record).** Four columns, **Cu d2 | Cu d1 (a2) | σ800 d2 κ off | σ800 d1 κ off**, through item 366's `_summary` and line format:
- row 1: `c_mode_f`, `c_k2_f`, off-block, anti;
- row 4a: f₀ one-point against SPEC C, and also 64 MHz·√(C_mode / C_in(d2));
- row 2b at own C, for p17_alone, ccw and cw: feed / other share and leak, plus the |k| folds;
- row 2c, both controls;
- the order moves d2 − d1 on copper and on σ = 800, side by side.
- The H25 lines:
  - `C_mode`(d2, Cu) in pF, beside **118.1186** (`C_MODE_CU_PF`) and **115.5209** (`C_IN_F_HUMAN_D2_PF`);
  - (C_mode − root) / root;
  - the window `PRED_CMODE_CU_D2_PF_WINDOW` hit or missed, and both point predictions;
  - beside them, the degree-1 offset (118.1186 − 118.2437) / 118.2437 = −1.058e-3.
- The 2b rows again at r = `ESR_OHM` (`ring_mode_rows_v2(…, esr_ohm=ESR_OHM)`), labelled "ESR (production convention)", beside the lossless ones.
- **Every 2b share and leak is labelled "circuit route (Norton feeds); on copper the circuit ↔ field identity at C_mode is the pre-registered H12c record (1.581167); not the in-model row".** Item 369's premise 2 applies: our SPEC row 2 is the field route, which H30 prints at degree 2. H12c's d1 in-model share was 0.998860 against the circuit's 0.993012 (`:6382`).

**Windows.** Each runs in the foreground with Bash `timeout: 660000`, with `-s`, and with `FEM_EM_ANS_SKIP_PRIVATE=1` set (nothing here writes a private file, but the knob is set in every `ANS-8` window). **The harness call is the whole Bash command.**
- **Host first:** `ls -la logs/wf12/ logs/wf7/`, then the BLOCKED-if-absent rule above.
- **W0** (smoke, `-n 1`; the gated test skips). The ungated `test_ans8_h25_coding` checks:
  - the knob parser, through `monkeypatch`;
  - names synthesis = slot synthesis on a synthetic D16 32×32 under a **shuffled** port order (18 classes, spread 0, bitwise);
  - (Σk) on a synthetic current vector;
  - (S) on a synthetic passive 32×32 at ≤ `SHIFT_RTOL`.

  `scripts/testing/run_and_log.sh ANS-8-H25-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans8_copper_degree2_rows.py -s -v --tb=short'"`
- **W1** (standard, `-n 1`; predicted ≤ 90 s, item 366's W1 shape):
  `scripts/testing/run_and_log.sh ANS-8-H25-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS8_H25=1 timeout -k 30 180 python3 -m pytest tests/environment tests/validation/test_ans8_copper_degree2_rows.py -s -v --tb=short'"`
- No durable capture: both windows are well under ≈ 10 min. No module edit after W0 starts (rule (i)).

**Traps already paid for.**
- **Two orders.** The d2 file is ordered with P17 first; the d1 files carry their own orders. Read every comparison **by name** (`test_wf12_f_human_copper_s32_degree2.py:380–389`) and never by index.
- **No slot keys on the d2 file.** `_slot_port_map` raises `KeyError` on it. Use the names map, validated by (Syn).
- **κ.** Both copper files and both `wf7` files are κ off. The tracked `metrics.json` `rows_v2` are degree 1 **κ on** σ = 800 (item 366's trap). They are not a column here, and never differenced against κ-off rows.
- **SPEC C ≠ any matrix's C_mode** (item 161 trap 1). 2b is at each matrix's own C_mode, wired inside `ring_mode_rows_v2`, and one matrix's C_mode is never passed to another.
- **The circuit route on copper is not our row 2** (item 369 premise 2; H12c 1.581167).
- **ESR.** `rows_v2` is lossless. The root 115.5209 pF carries r = 0.2 Ω; `C_mode` does not, which (S) licenses. Never compare a lossless share with the ESR in-model share as if both were one quantity.
- **Sense:** ccw = (1, −1j) (`SPEC_DRIVES`). H26's root is the P17 (1, +1j) reading (ruling (2)).
- **The `h12c` files are refused by the step-5 writer** (`REFUSED`, `:55`). Opening them read-only is safe. Never write under `logs/wf12/`.
- **Gzipped logs:** cite `…a2-W1.log.gz` line numbers of the decompressed text. Never decompress in place.
- Never import a `test_*` name; never `-k`; no pipe; never delete a log; `-s`; do not import or run `aed/`.

**Scope.** Rows from four saved 32×32 matrices (copper d2 / d1 on the hole mesh, σ = 800 d2 / d1 boxed), κ off, 64 MHz, circuit route, no solve. **Not claimed:**
- any AED value, verdict or level;
- the in-model (field-route) row 2 at degree 2 (H30);
- f₀, Q, match or tuning accuracy;
- that `C_mode`(d2, Cu) is the root (it is printed beside the root, never asserted);
- the `WF-12` glyph;
- an h-statement.

**Negative result.** Report and stop. Park on `attempt/ANS-8-H25-*` and mark item 516 BLOCKED in the same commit (rule (d)).
- **(F) red:** the d2 file is not what entry 36 wrote. Stop; the review compares it with `:8032`.
- **(Syn) red:** the names synthesis is wrong. Read no d2 class mean.
- **(L-d1) red with (Syn) green:** the d1 records do not reproduce. Print the table and park.
- **(L-d2) red with (L-d1) green:** the saved d2 matrix does not reproduce its window. File a known-issues entry, and H30 / the weekly do not read it.
- **(G1′) red:** a passivity or reciprocity finding. File a known-issues entry.
- **(B0) / (Σk) red:** a `ring_modes` regression. Stop.
- **(S) red:** the ESR identity fails on copper. That is a `PORT-24` finding: known-issues, and the weekly rules.
- **A floor-asserted negative control below its floor:** report the factor; the review rules.
- Every row, move, prediction and printed control is a reading, never a red. A `C_mode`(d2, Cu) outside the 1 % window is a reading for the weekly.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** No glyph; rule (j) does not bind.
- `ANS-8` (✅ on step 1; step 2 open) gains one record sentence: `C_mode`(d2, Cu) beside 118.1186 / 115.5209 pF, and the circuit-route 2b ccw share / leak, labelled.
- §10 H25 gains a DONE line.
- `docs/testing/xl-pending.md` entry 36's readout clause "`ring_mode_rows_v2` on the saved matrix" is marked met by this item (a status line, the review's edit; never the command).
- The record feeds the 10-19 weekly's reading of whether the degree-2 copper `C_mode` sits with the root (H30's comparand context). It moves no known-issues entry unless a red files one.

**Execution order.** `ls` → code → W0 → W1 → commit (module, both harness logs, the §7 sentence, §10 H25 DONE) → journal.

**Spawn line (if delegated):** "Use the implementer agent to execute `ANS-8` §9 item 516 (H25) per `docs/planning/chunks/ANS-8.md` (read the brief and any review ruling beneath it first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W0 then W1, in the foreground, with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`, and `FEM_EM_ANS_SKIP_PRIVATE=1`. If `logs/wf12/h16b-cu-s32-d2.npz` is absent on the host, stop and mark item 516 BLOCKED; never re-solve. Create only `tests/validation/test_ans8_copper_degree2_rows.py`; no `src/` edit, and no edit to `ans:8`, `metrics.json`, `COMPARISON.md` or any existing test module. Open every npz read-only. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. You are the executor: do not spawn agents, and never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ANS-8-H25-*`, mark item 516 BLOCKED in the same commit, and report."

### Review ruling 2026-10-18 03:00 (items 513 H29, 516 H25)

Both accepted. Spot-check: `tests/unit/test_ans8_v3_definition.py:156` (`assert "118.2437" in sub and "2 212" in sub`) re-read. Entry 36's matrix `logs/wf12/h16b-cu-s32-d2.npz` is on disk (23 077 B, 02:13); the window is green (`20261018T070008Z_WF-12-step5-H16c.log:8032–8036`, re-read).
- **Q1:** the drafter's placement. The degree-1 values go in a labelled "do not enter" record subsection inside the span, so the item-411 test stays unedited.
- **Q4:** enter **C = 115.5209 pF**, the weekly's text and the high-Q convention SPEC v3 used. The exact C_p (115.5109 pF, 8.7e-5 below) is printed beside it in the record subsection.
- **Q2 / Q3:** named for the weekly, not edited.
  - Q2: row 1's readout says "within one grid resolution", written for a 0.10 pF grid. The degree-2 root's bracket is 0.50 pF.
  - Q3: `SPEC.md:320` "(ours 21.5)".
- **H25 tier: standard, `-n 1`.** It is NumPy on saved matrices; the matching step-1c window took 74 s.
- **H25's 1 % prediction stays printed, never asserted.** The drafter's two extrapolations (115.399 / 113.02 pF) are printed beside it.
- **Entry 36's decision-rule wording** ("step 2's degree-2 column is read from this matrix") is named for the weekly. Step 2's row 2 at the root is H30, from entry 41. Item 516 meets entry 36's `ring_mode_rows_v2` readout clause.

### §9 item 576 — 2026-10-20 (§10 H30a, chain step (rule (6), provisional); lane 1)

*Drafted read-only for the 2026-10-20 18:00 daily review; the review rules on it.*

**Licence.** §10 chain H, "H30 (lane 1 → heavy, no solve; serial on entry 41). `ans:8` v3 rows 3–5 at degree 2 at the root from entry 41's record (`_stage_v3` reading `h9p-node-115.5209-d2.json` and the lattice table), per accepted watt as SPEC v4 defines it; target: (Lin) and (P_acc-1) at their imported bands, the rows printed beside the degree-1 column, plus the v3b diagnostic at degree 2" (`PROJECT_PLAN.md:4849–4854`). The 10-19 weekly re-posed known-issues 2026-10-11's trigger onto H30: retire if ΣX_s / P_acc < 8.3e-3 at degree 2, otherwise re-title the entry a finding on the lumped-sheet model (`PROJECT_PLAN.md:4789–4795`; `docs/testing/known-issues.md:814`, `:824`). **H30a is the pre-wiring step toward H30's own target, under rule (6). It is the first provisional step on chain H since the 10-19 weekly.** It lands the code H30 needs and proves it on records that already exist, so H30 on 10-22 is one read. **No band, tolerance, Done-when, validation target or retirement condition is set or moved.** H30a reads no degree-2 root figure as H30's result.

**Why H30 needs code (read from the code, 2026-10-20).** `_stage_v3` cannot do what H30 asks, for four reasons:
1. It asserts `FEM_EM_WF10_DEGREE` unset (`examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py:411`). It reads only the degree-1 `NODE_JSON` (`:415–424`; `tests/validation/test_wf10_f_human_tuned_b1_sar.py:102`).
2. **It is a solve.** It calls `build_tuned_quadrature` (`08_…py:427`). At degree 2 that is entry 35's builder: build 145.2 s plus a node pair of 367.1 s at `-n 16`, peak 86.81 GiB (`docs/testing/xl-pending.md:1959`). That is not "no solve".
3. Its reproduction anchors are degree-1 literals:
   - (Rep-P) / (Rep-Part) (`:125–127`, `:514–522`);
   - (Rep-B) (`:117`, `:608–612`);
   - (Rep-S) (`:120–121`, `:649–653`).
4. It hard-codes `"degree": 1` (`:674`) and writes `v3-rows.json` (`:112`, `:708`).

The degree-2 node record is the JSON the builder writes. Its keys (`test_wf10_f_human_tuned_b1_sar.py:616–630`) are `c_f`, `c_in_pf`, `r_ohm`, `freq_hz`, `order`, `feeds`, `s_node`, `w`, `y_in`, `partitions[sense]` (`p_src`, `p_sheet_by_id`, `p_sheets`, `p_phantom`, `p_rest`, `p_surf`), `degree` and `node_knob`. Together with the 10 g module's c\* line (`tests/validation/test_wf10_f_human_tuned_sar10g.py:448–455`), these carry what H30 needs for:
- P_acc,SPEC (line cut and wave form);
- row 4 per accepted watt;
- row 5 (both normalisations);
- the gap;
- X_s per feed sheet, ΣX_s and ΣX_s / P_acc.

The feed currents are not stored, but the record determines them. `s_node` is built by the current route, S_ij = 2V_i/V_s − δ_ij (`tests/validation/test_wf10_f_human_copper_esr_root.py:154–158`). So under the drive w, **V_port = ½ V_s (𝟙 + S) w** and **I = (w − V_port) / Z_th**, with Z_th and V_s from `:144–151`.

**Not in any record — row 3, row 6, and CS_s.**
- Row 3 needs the five SPEC probes at z_c (`08_…py:133–134`, `:594–621`). The b1 module prints only the masked-sample mean \|B₁⁺\| and the purity at the sample point nearest the centre, which is not O (`test_wf10_f_human_tuned_b1_sar.py:757–781`).
- Row 6 is the purity at O.
- CS_s needs the field on the feed sheet: A_s, J_s, L_s, T_s (`08_…py:550–561`).

H30a cannot wire these without a field. The choice is the review's or the weekly's (for-the-review, item 1).

**Drafter's check (unlogged stdlib Python on the three node records on this box, a reading, not a window).** The record route reproduces item 386's degree-1 window to every printed digit:
- P_acc,SPEC 3.057629563251e-02 W (`20261011T170627Z_ANS-8-v3b-V3.log:6194`);
- the partition terms (`:6195`);
- X_s 1.329324349e-04 / 1.213651853e-04 W (`:6205`, `:6207`);
- ΣX_s = the gap, 2.542976202e-04 W;
- row 5 0.268112 / 0.019926 / 0.703644, sum 0.991683;
- the item-361 negative control 4.9700e-02 and the ESR-blind miss 0.71196.

So the code is small, and the proof is a reproduction.

**`ANS-8` H30a: the v3 record route — rows 4–5 per accepted watt, the gap, and X_s / ΣX_s / P_acc from a node record and its harness log(s), no builder. It is proven on the degree-1 record against item 386's `v3-rows.json` and rehearsed on the two degree-2 records already on this box. H30's own gated test lands skipped.**
- **Executor:** implementer.
- **Tier:** standard, `-n 1`, no `mpiexec`. The work is NumPy, JSON and regex, with no collective, so no rank-local quantity exists (item 516's precedent).
- **Class:** chain step (no solve). **Cost:** ≈ 6.5 measured / ≈ 35 predicted slot-minutes: 15 fixed, ≈ 12 code, W0 ≈ 1.5, W1 ≈ 1.5, commit and journal ≈ 5. The window prices come from item 516's W0 / W1: 96 s harness, test body 0.19 s (`20261018T095551Z_ANS-8-H25-W0.log`, `20261018T095559Z_ANS-8-H25-W1.log`).
- **Depends on:** nothing. It runs before or after entry 41 (`xl`, Thursday 10-22 02:00) and reads nothing entry 41 writes.
- **Family:** `ANS-8` step 2's degree-2 column; attempt 1 of 4.

**What exists (host, 2026-10-20; every input gitignored or tracked-log, opened read-only).**
- **Degree-1 record:** `logs/wf10/h9p-node-118.2437.json` (5 495 B, written 2026-10-18 04:50). Its writer is `MAT-15` step 6b W3:
  - `wrote` at `20261018T094815Z_MAT-15-step6b-W3-regate.log:6299`;
  - P_w line `:6295`, (Lin) `:6296`.
  - It is **not** item 386's window: the builder rewrites the file on every call.
- **Degree-1 lattice log:** `20261013T111006Z_WF-10-step4b-W1-regate.log`:
  - table header `:6371`, c\* `:6390` (4.852895959e-03 W/kg, label `r=45 dz=+0 az17+0`);
  - P_w `:6325`, (Sigma_P) `:6396`, phantom bbox line present.
- **Degree-1 comparand:** `logs/ans8/v3-rows.json` (9 696 B, 2026-10-11 12:09). This is item 386's V3 output (`DEFAULT_RESULTS_DIR` + `V3_JSON_NAME`, `08_…py:61`, `:112`). Its keys are at `:673–705`.
- **Entry 35's degree-2 record at the degree-1 node (detuned):** `logs/wf10/h9p-node-118.2437-d2.json` (5 471 B, 2026-10-16 02:19). It has **no `node_knob` key**, because it was written before item 511. Its log `20261016T070008Z_WF-10-step4-H20.log` has **entry 41's exact shape** (both modules, two builder calls):
  - `wrote` at `:6834` and `:13392`; P_w `:6830` / `:13388`; (Lin) `:6831` / `:13389`;
  - table `:13434`, c\* `:13453` (2.873208988e-03 W/kg).
- **Entry 42's degree-2 record at the root:** `logs/wf10/h9p-node-115.5209-d2.json` (5 475 B, 2026-10-19 02:09; `degree` 2, `node_knob` `"d2root"`). Its log `20261019T070008Z_MAT-15-step6-P22-d2root.log`:
  - `wrote` `:6790`, P_w `:6786`, (Lin) `:6787`;
  - **no 10 g table** (the `MAT-15` module).
  - **Entry 41 will overwrite this file at 10-22 ≈ 02:10–02:30.** H30a reads it tonight; H30 reads entry 41's rewrite.
- **Reusable code:**
  - `p_acc_line`, `p_acc_wave`, `row5_fractions` (`08_…py:286–313`);
  - `_z_rc`, `_thevenin_feed_impedance`, `_source_voltage`, `_active_y_in`, `ESR_OHM`, `JSON_DIR`, `Y0_CODING_BAND` (`test_wf10_f_human_copper_esr_root.py:188`, `:144`, `:149`, `:198`, `:53`, `:51`, `:55`);
  - `_node_json_path`, `NODE_JSON`, `C_IN_F_HUMAN_PF`, `C_IN_F_HUMAN_D2_PF` (`test_wf10_f_human_tuned_b1_sar.py:155–171`, `:102`, `:68`, `:78`);
  - `parse_table`, `CSTAR_RE`, `HEAD_RE`, `PPH_CCW_RE`, `EXPECTED_POINTS` (`tests/validation/test_wf10_lattice_d2.py:102–170`, `:73`, `:72`, `:88`, `:59`);
  - `REFERENCE_IMPEDANCE_OHM` (`tests/validation/test_port_package_sparameters.py:140`).

**Code.** No `src/` edit. No edit to `_stage_v3`, `_stage_r`, `main()`, any `V3_*` constant, `metrics.json`, `COMPARISON.md`, `SPEC.md`, any `test_wf10_*` module, `tests/unit/test_ans8_v3_rows.py` or the H27 module.
- **`08_…py`, additive only (pure functions beside `row5_fractions`):**
  - `record_port_state(s_node, w, z_th, v_s)` returns `(V_port, I_feed)` by the two formulas above.
  - `v3_record_rows(rec, sar10g_cstar_w_per_kg=None)` works on the `ccw` sense only. It asserts `rec["r_ohm"] == ESR_OHM` and `rec["feeds"] == rec["order"][:2]`.
    - Z_RC = `_z_rc(rec["c_f"], rec["freq_hz"], rec["r_ohm"])`; z0 = `REFERENCE_IMPEDANCE_OHM`.
    - Returns:
      - P_acc,SPEC = `p_acc_line(w, I_feed, z_th, z_rc)`, and P_acc,wave = `p_acc_wave(s_node, w, v_s, z0)`;
      - r_acc1;
      - p_rc30, p_esr_feeds, p_phantom, p_surf, p_rest, p_src, p_361;
      - gap and r_acc2;
      - row 5 of P_acc and of the partition (`row5_fractions`);
      - per feed: `x_w` = P_sheet,field − ½ Re(1/Z_th)\|V\|², `lumped_norton_w`, `p_sheet_field_w`;
      - `sum_x_w`, `sum_x_over_p_acc`, `terminal_w`, `norton_lumped_w`;
      - the negative-control misses (item 361's P_acc against the wave form; ESR-blind);
      - when a SAR is given, row 4 per W of each normaliser;
      - `degree`, `c_in_pf`, and `not_in_record: ["row3", "row6", "cs_s"]`.
    - These are `_stage_v3`'s formulas at `:491–537` and `:562–582`, with the record route's V / I in place of the field currents.
  - `node_log_readings(log_path, node_basename)` reads the harness log:
    - the **last** `[WF-10 step2] wrote …/<basename>` line (asserted present) and its full path;
    - the **last** `[ccw] w = … (P_w)` line before it (P_src, sumP_sheet, P_phantom, P_surf, at `.9e`; format `test_wf10_f_human_tuned_b1_sar.py:604–608`);
    - **every** `[ccw]` / `[cw] (Lin) currents rel … Y_in current route vs _active_y_in …` pair (`:609–611`).
- **New module `tests/validation/test_ans8_v3_record_rows.py`.** It loads the script by `importlib` under `"ans8_h30a"` (the item-516 pattern). Knobs follow item 366's `_enabled` / `_need` pattern (`tests/validation/test_ans8_degree2_rows.py:73–84`): unset ⇒ skip, `"1"` ⇒ run with absent inputs **failing**, anything else ⇒ `ValueError`.
  - `test_ans8_h30a_coding`, ungated (W0):
    - (C1) on item 369's synthetic two-feed network, rebuilt as `tests/unit/test_ans8_v3_rows.py:48–58` builds it (seed 369, copied, not imported), `record_port_state` from `_current_route_s`'s S equals the current-route V = w − Z_th·I and I for four drives, ≤ `Y0_CODING_BAND`.
    - (C2) `node_log_readings` on a synthetic three-line text built from the format strings above. It must pick the last `wrote`, the last P_w before it, and both (Lin) lines.
    - (C3) the knob parser, through `monkeypatch`.
  - `test_ans8_h30a_degree1_reproduction` (`FEM_EM_ANS8_H30A=1`): the d1 record, both d1 logs and `v3-rows.json`.
  - `test_ans8_h30a_degree2_rehearsal` (`FEM_EM_ANS8_H30A=1`) runs two records:
    - (R35) `_node_json_path(2)` with entry 35's log, used as both node log and lattice log;
    - (R42) `_node_json_path(2, c_pf=C_IN_F_HUMAN_D2_PF)` with entry 42's log, which has no lattice, so row 4 prints "not in this log".
  - `test_ans8_h30_degree2_root` (`FEM_EM_ANS8_H30=1`; **lands skipped, H30 runs it**):
    - record `_node_json_path(2, c_pf=C_IN_F_HUMAN_D2_PF)`;
    - log: exactly one match of the glob `20261022T*_WF-10-step4-H28.log` under `docs/testing/logs/` (H27's condition-2 pattern, `test_wf10_lattice_d2.py:96–99`), used as both node log and lattice log;
    - it asserts the record's `node_knob == "d2root"` and `degree == 2`;
    - it prints every row beside the degree-1 column read from `v3-rows.json`, and writes `/workspace/logs/ans8/v3-rows-d2-115.5209.json` (asserted ≠ `V3_JSON_NAME`).
  - H30a's two tests write nothing.

**Anchors (asserted, every band imported).** For each record:
- **(Prov) the JSON is the one the log wrote.**
  - The log's last `wrote` path equals the expected container path (`_node_json_path(…)` under `JSON_DIR`).
  - `rec["degree"]`, `rec["c_in_pf"]` and `rec.get("node_knob", "")` are as expected: d1 1 / 118.2437 / ""; R35 2 / 118.2437 / absent; R42 2 / 115.5209 / "d2root".
  - The P_w line's P_src, sumP_sheet, P_phantom and P_surf match the JSON's `partitions["ccw"]`, each ≤ `REUSE_REPRODUCTION_RTOL` (`tests/validation/test_port_birdcage_ring_matrix.py:171`, = 1e-6 via `:164`, the definition `_stage_v3` imports at `:371`). Predicted ≤ 6e-10 (`.9e` prints).
- **(Rec-Y)** the JSON's `y_in["ccw"]` against `_active_y_in(s_node, w_ccw, z0)` ≤ `Y0_CODING_BAND`. Predicted 0: same function, same floats, exact JSON round-trip.
- **(Lin)** every (Lin) currents reading in the node log ≤ `UNIT_WEIGHT_CONTROL_RTOL` (`tests/validation/test_port_drive_superposition.py:303`), and its Y_in route reading ≤ `Y0_CODING_BAND`. This re-asserts the window's own field-side gate from its log. Logged values: d1 6.287e-16 / R35 1.212e-15, 9.610e-16 / R42 9.692e-16.
- **(P_acc-1)** record-route line cut against the wave form ≤ `DISCRETE_IDENTITY_RTOL` (`tests/validation/test_birdcage_power_identity.py:227`), as `_stage_v3:506–507`. Drafter: 1.1e-16 / 1.1e-15 / 6.0e-16; predicted ≤ 2e-15.
- **(Prov-L) where a lattice log is read:**
  - `parse_table` succeeds (one header, `EXPECTED_POINTS` rows);
  - its `p_phantom_ccw` matches the JSON's ≤ `REUSE_REPRODUCTION_RTOL`.
- **(Rep-d1) degree 1 only, against `v3-rows.json`, each ≤ `REUSE_REPRODUCTION_RTOL`:**
  - `power_w`: `p_acc_spec`, `p_acc_wave`, `p_acc_item361`, `p_src`, `p_rc30`, `p_esr_feeds`, `p_phantom`, `p_surf`;
  - `row5` and `row5_partition`: phantom, copper and ESR;
  - `gap.line_cut_minus_partition_w`, `gap.terminal_w`, `gap.norton_lumped_w`;
  - per feed, `gap.mechanism[pid]`: `x_w`, `lumped_norton_w`, `p_sheet_field_w`;
  - `row4.sar10g_per_w` and `per_partition.row4_sar10g_per_w`.

  Predicted ≤ 1e-9. The drafter's stdlib values agree to every printed digit. X_s is a difference of two ≈ 3e-2 W terms, a ×≈ 230 cancellation, which is still far inside 1e-6.
- **(Rep-S) degree 1:** the c\* label equals `V3_CSTAR_LABEL_RECORD`, and SAR10g(c\*) matches `SAR10G_CSTAR_H9P_RECORD_W_PER_KG` ≤ `REUSE_REPRODUCTION_RTOL` (`08_…py:120–121`).
- **Printed, as `_stage_v3` printed them (`:577–582`):** P_src against the Norton lumped total, and ΣX_s against the gap and against p_361 − terminal. Predicted ≤ 1e-12. On the record route these follow from P_src = Norton.

**Negative controls (rule (e)).**
- **Asserted, backed (item 369 read 4.9700e-02, `_stage_v3:526–528`):** item 361's P_acc in (P_acc-1) misses by ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL` = 1e-5 (`tests/validation/test_wf12_f_human_copper_matrix.py:61`) on all three records. Drafter: 4.9700e-02 / 9.7512e-02 / 6.2316e-02. The ceiling is the feed-ESR share of P_acc, ≈ 0.04–0.1, so the bar has three orders of headroom.
- **Asserted, floor imported, factor predicted (new comparison):** (Prov) wrong-record. The d1 JSON against entry 35's last P_w line, and the R42 JSON against entry 35's, each miss by ≥ `BLIND_FACTOR` × `REUSE_REPRODUCTION_RTOL` = 1e-5. Drafter, on P_src: 3.18e-02 and 3.36e-03. This shows (Prov) can tell entry 41's record from entry 35's.
- **Asserted, floor imported, factor predicted:** feed-label swap on d1. With X_s(P17) ↔ X_s(P21) exchanged, (Rep-d1)'s `x_w` misses by ≥ 1e-5. Drafter: \|X₂₁/X₁₇ − 1\| = 8.70e-02.
- **Printed only:** the ESR-blind miss (`_stage_v3:529–530`). Drafter: 0.71196 / 0.73012 / 0.71896.

**Predictions (written first; printed, never asserted).** Values are d1 | R35 (degree 2 at 118.2437 pF, detuned, a labelled record) | R42 (degree 2 at the root, entry 42's record):
- P_acc,SPEC: 3.057629563251e-02 | 1.907593125981e-02 | 2.910591023476e-02 W.
- X_s P17 / P21:
  - d1 1.329324349e-04 / 1.213651853e-04 W;
  - R35 2.050749253e-04 / 1.899495613e-04 W;
  - R42 3.262132656e-04 / 2.994809404e-04 W.
- **ΣX_s / P_acc: 8.316822e-03 | 2.070801e-02 | 2.149715e-02.**
- Row 5 (phantom / copper / ESR, sum):
  - d1 0.268112 / 0.019926 / 0.703644, 0.991683;
  - R35 0.250469 / 0.019414 / 0.709409, 0.979292;
  - R42 0.260954 / 0.020083 / 0.697465, 0.978503.
- Row 4 per W of P_acc: d1 0.1587143 | R35 ≈ 0.150620 | R42 not in its log.
- **For H30 (printed by H30a beside R42, never asserted, never H30's reading):**
  - Entry 41 runs the same builder, knob, mesh, degree and rank count as entry 42, so its node record should reproduce R42's at `REUSE_REPRODUCTION_RTOL`. **ΣX_s / P_acc at the degree-2 root is then predicted ≈ 2.15e-02, ≈ 2.6× the degree-1 value: the excess grows with order.**
  - Row 4 at the root is predicted in (0.145, 0.170) W/kg/W, point ≈ 0.159: d1 0.1587 moved by entry 42's per-accepted-watt SAR10g d1 → d2 move, +0.258 %.
  - Known-issues 2026-10-11 does not move on these numbers. Its trigger is H30's.

**Windows.** Each runs in the foreground with Bash `timeout: 660000`, pytest `-s`, container `timeout -k 30`, and `FEM_EM_ANS_SKIP_PRIVATE=1`. **The harness call is the whole Bash command.**
- **Host first:** Glob the four gitignored inputs above (the three JSONs and `logs/ans8/v3-rows.json`; the four logs are tracked). Any absent ⇒ BLOCKED, never re-solve.
- **W0** (smoke, `-n 1`; the gated tests skip):
  `scripts/testing/run_and_log.sh ANS-8-H30a-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans8_v3_record_rows.py -s -v --tb=short'"`
- **W1** (standard, `-n 1`; predicted ≤ 100 s harness):
  `scripts/testing/run_and_log.sh ANS-8-H30a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS8_H30A=1 timeout -k 30 180 python3 -m pytest tests/environment tests/validation/test_ans8_v3_record_rows.py -s -v --tb=short'"`
- **H30's window** is recorded here for the review and **not run by H30a**. It is W1 with `FEM_EM_ANS8_H30=1` in place of `FEM_EM_ANS8_H30A=1` and the ID `ANS-8-H30-W1`.
- No durable capture is needed (both windows are well under ≈ 10 min). No module edit after W0 starts (rule (i)).

**Traps already paid for.**
- **`_stage_v3` is a solve.** Never run stage V3 with `FEM_EM_WF10_DEGREE=2`, and never edit it. The 10-19 text's "→ heavy" came from reading it as the route. The record route is standard, `-n 1`.
- **The builder rewrites the node JSON on every call.**
  - The d1 file is `MAT-15` 6b W3's, not item 386's, so (Rep-d1) is at `REUSE_REPRODUCTION_RTOL`, never bitwise.
  - Entry 35's (and entry 41's) log writes its file twice. Provenance is to the **last** `wrote` line and the last P_w line before it. Assert **every** (Lin).
- **Entry 35's JSON has no `node_knob` key** (written before item 511). Read it with `.get`.
- **H27 trap 2:** in a two-module log, the lattice log's P_phantom is the sar10g module's, the last P_w line before the table. `parse_table` already does this; never re-implement it.
- **(P_acc-1) on the record route is the S ↔ V / I algebra,** not `_stage_v3`'s field-current comparison. The field side is the window's (Lin), re-asserted from its log. Label both.
- **Rows 3 / 6 and CS_s are not in any record.** Never fill them from the mean \|B₁⁺\| or the "centre purity" line. The first is a masked-sample mean; the second is read at the nearest sample point, not at O.
- **`REUSE_REPRODUCTION_RTOL` has two definitions:** 1e-6 at `test_port_birdcage_ring_matrix.py:171` and 1e-12 at `tests/validation/test_port19_factor_reuse.py:82`. Import the first, as `_stage_v3` does.
- **`w["ccw"][1]` is 6.123e-17 − 1j** (the convention code). Read it from the record; never substitute −1j.
- **Do not read ΣX_s / P_acc against 8.3e-3 here.** The degree-1 value 8.3168e-03 sits just above the rounded threshold, which is the weekly's wording, and the comparison is H30's.
- Never import a `test_*` name; never `-k`; no pipe; never delete a log; `-s`; never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`.

**Scope.** This item is a record route on existing records, with no solve. **Not claimed:**
- H30's reading, or any degree-2 root row as a result;
- row 3, row 6 or CS_s at degree 2;
- the known-issues 2026-10-11 retirement or re-title;
- any AED value, verdict or level;
- entry 41's gates;
- an h- or order-convergence statement.

**Negative result.** Report and stop. Park on `attempt/ANS-8-H30a-*` and mark item 576 BLOCKED in the same commit (rule (d)).
- **(Prov) red on d1:** the record is not what the cited log wrote. Name both mtimes; the review re-points the log.
- **(Rep-d1) red with (Prov) green:** the record route does not reproduce item 386, and H30 is not wired. Print the per-key table.
- **(P_acc-1) or (Rec-Y) red:** an algebra or coding error in `record_port_state`. Stop.
- **(Lin) red in a logged window:** a contradiction with that window's own green footer. Report it to the review (`log-pathologist`).
- **R35 / R42 red with d1 green:** the degree-2 path or the naming. H30 is not wired; park.
- **A negative control below its floor:** report the factor.

Every prediction is a reading. There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** No glyph; rule (j) does not bind.
- §10 gains an H30a DONE line.
- **H30's risk is retired:** H30 becomes one standard `-n 1` window of ≈ 2 min, not "heavy", and its command is fixed above.
- The 10-21 weekly gets two inputs before entry 41 runs:
  - the R42 prediction that the excess grows with order (known-issues 2026-10-11 is expected to be re-titled, not retired, at H30);
  - the fact that row 3 / row 6 / CS_s cannot come from entry 41's record as queued.
- No known-issues entry moves on this item.

**For the review.**
1. **Rows 3 / 6 and CS_s at degree 2** are in H30's text ("rows 3–5", "the v3b diagnostic") but in no record entry 41 writes. Options:
   - (a) H30 prints rows 4–5, X_s and ΣX_s / P_acc and names the rest "not in entry 41's record". This narrows H30's target, so it is the 10-21 weekly's call.
   - (b) A printed-only addition to `test_wf10_f_human_tuned_b1_sar.py` before 10-22 02:00: the five SPEC probes \|B₁±\| at z_c, and the feed sheets' A_s / J_s / L_s / T_s under ccw, written into the node JSON. This edits a queued window's module and needs a degree-1 re-gate; the 10-21 weekly sits before the window.
   - (c) H30 as a degree-2 stage V3 solve at heavy `-n 12` in the 128 GiB service. Entry 35's 86.81 GiB peak was at `-n 16`; this is unpriced at `-n 12` and contradicts "no solve".

   The drafter's default is (a), named for the weekly. CS_s / X_s is predicted 1 by the sheet model's own definitions (½Re(Y_s)L_s = P_sheet,field and the lumped term read equal at d1, `20261011T170627Z_ANS-8-v3b-V3.log:6205`, `:6207`). X_s carries the content.
2. **The R42 prediction** (ΣX_s / P_acc ≈ 2.15e-02 at the degree-2 root) says known-issues 2026-10-11's "shrinks with order" branch will not fire at H30. The review may name this for the 10-21 weekly. It moves nothing tonight.
3. **H30's tier:** standard `-n 1`, not "→ heavy".
4. Whether stage R should fold H30's `v3-rows-d2-115.5209.json` into `metrics.json` / `COMPARISON.md`. That is H31's or a later step's; it is not in H30a.
5. The 8.3e-3 threshold against the degree-1 8.3168e-03 is rounding in the weekly's wording. It is named, not edited.

**Execution order.** Host Glob → code → W0 → W1 → commit (script, new module, both harness logs, §10 H30a DONE) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 576 — `ANS-8` H30a (§10 chain H, provisional) — per `docs/planning/chunks/ANS-8.md` §9 item 576. Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Run W0 then W1, in the foreground, with `timeout: 660000` on every Bash call that runs a window, pytest with -s, container `timeout -k 30`, and `FEM_EM_ANS_SKIP_PRIVATE=1`. Glob `logs/wf10/h9p-node-118.2437.json`, `logs/wf10/h9p-node-118.2437-d2.json`, `logs/wf10/h9p-node-115.5209-d2.json` and `logs/ans8/v3-rows.json` first; if any is absent, stop and mark item 576 BLOCKED; never re-solve. Never run stage V3 and never set `FEM_EM_WF10_DEGREE`; do not edit `_stage_v3`, `_stage_r`, `main()`, any `V3_*` constant, `metrics.json`, `COMPARISON.md`, `SPEC.md`, any `tests/validation/test_wf10_*` module, `tests/unit/test_ans8_v3_rows.py` or any `src/` file. Add only the pure functions named in the brief to the `ans:8` script and create only `tests/validation/test_ans8_v3_record_rows.py`; `test_ans8_h30_degree2_root` lands skipped and is never run with `FEM_EM_ANS8_H30=1` by you. Open every JSON read-only and never write under `logs/wf10/`. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ANS-8-H30a-*`, mark item 576 BLOCKED in the same commit, and report."

### Review ruling 2026-10-20 18:00 (item 576, H30a; rule (6), provisional)

Binds over the brief. (1) **Queued as written:** the record route (no builder) for rows 4–5 per accepted watt, X_s per feed sheet and ΣX_s / P_acc, proven on the degree-1 record against item 386 at `REUSE_REPRODUCTION_RTOL` with (P_acc-1) at `DISCRETE_IDENTITY_RTOL`, rehearsed on entries 35 / 42; H30's own gated test lands skipped. Standard, `-n 1`; ≈ 7.8 measured (chain no-solve median) / 35 predicted. **It must run before entry 41 overwrites `logs/wf10/h9p-node-115.5209-d2.json` (Thursday 10-22 02:00);** if it has not landed by the 10-22 03:00 review, the entry-42 rehearsal is dropped, not re-solved. (2) **Rows 3 / 6 and CS_s are in no record entry 41 writes.** The drafter's option (a) is this item's scope: H30a prints "not in record" for them. Whether H30's target narrows to rows 4–5 + X_s + ΣX_s / P_acc, or entry 41's module gains a printed-only probe / sheet-field write before its window (option (b), which edits a queued window's module and is not licensed to a slot), is **returned to the 2026-10-21 weekly**. Option (c), a degree-2 stage-V3 solve, is not queued. (3) **The entry-42 reading** (ΣX_s / P_acc ≈ 2.15e-02 at the degree-2 root against 8.32e-03 at degree 1, the drafter's unlogged arithmetic) is a prediction only until H30a's rehearsal prints it through the harness; it is named to the weekly as a likely outcome of H21's trigger, not banked. (4) H30's tier (standard `-n 1`) and the "< 8.3e-3" threshold wording are the weekly's; folding H30's output into `metrics.json` / `COMPARISON.md` is H31's or later.

### Outcome 2026-10-20 (19:30 slot, take-next; item 576 DONE)

W0 `20261021T004356Z_ANS-8-H30a-W0.log` (5 passed / 3 skipped, 3 s); W1 `20261021T004407Z_ANS-8-H30a-W1.log` (15 passed / 1 skipped, 100 s harness, `-n 1`, standard). Every anchor green; every prediction landed to its printed digit.
- **d1:** (Prov) P_w vs JSON ≤ 6.9e-11; (Rec-Y) 0; (Lin) 6.287e-16 / 4.627e-16; (P_acc-1) 4.539e-16; (Prov-L) 2.715e-12; (Rep-S) 0 on label r=45 dz=+0 az17+0; **(Rep-d1) worst 4.339e-11 over 25 keys** (row 4 per W; the rest ≤ 5.6e-13). P_acc,SPEC 3.057629563251e-02 W, X_s 1.329324349e-04 / 1.213651853e-04 W, ΣX_s / P_acc 8.316822e-03, row 5 0.268112 / 0.019926 / 0.703644 (sum 0.991683), row 4 0.1587143 W/kg/W.
- **R35** (degree 2 at 118.2437 pF, detuned; last of two `wrote` at `:13392`, all four (Lin) ≤ 1.212e-15): P_acc 1.907593125981e-02 W, ΣX_s / P_acc 2.070801e-02, row 5 sum 0.979292, row 4 0.1506196.
- **R42** (degree 2 at the root, entry 42's record): P_acc 2.910591023476e-02 W, X_s 3.262132656e-04 / 2.994809404e-04 W, **ΣX_s / P_acc 2.149715e-02** (a rehearsal on entry 42's record, not H30's reading), row 5 sum 0.978503, row 4 not in its log.
- **Negative controls (asserted):** item 361's P_acc in (P_acc-1) 4.9700e-02 / 9.7512e-02 / 6.2316e-02; wrong record (P_src vs entry 35's P_w) 3.0775e-02 (d1), 3.3445e-03 (R42); feed-label swap 8.7016e-02. ESR-blind (printed) 0.71196 / 0.73012 / 0.71896.
- **Printed lines that differ from the brief's prediction (printed only, never asserted):** ΣX_s against the gap reads ≈ 1e-10, not ≤ 1e-12. That is the field P_w partition residual (≈ 5e-13 of P_src) amplified by the ×≈ 340 cancellation in the gap. "ΣX_s vs p_361 − terminal" reads rel 2.000: p_361 − terminal = −ΣX_s by construction. Item 386's own line prints the same sign (`20261011T170627Z_ANS-8-v3b-V3.log:6208`: +2.542976202e-04 vs −2.542976202e-04). The magnitudes agree.
- Rows 3 / 6 and CS_s printed "not in record". `test_ans8_h30_degree2_root` landed skipped; H30's window is W1 with `FEM_EM_ANS8_H30=1`.

### Weekly ruling 2026-10-21 — H30 narrowed; H32 posed (full text PROJECT_PLAN §10)

- **H30 (target narrowed):** rows 4-5 per accepted watt, X_s per feed sheet and ΣX_s / P_acc from the record entry 41 writes (`h9p-node-115.5209-d2.json`), through the H30a route (one standard `-n 1` window, `FEM_EM_ANS8_H30=1`); (Lin) and (P_acc-1) at their imported bands, the rows printed beside the degree-1 column. Rows 3 / 6 and CS_s are "not in record" at degree 2 — recorded, not posed: CS_s / X_s = 1 by the sheet model definitions at degree 1, so X_s carries the content, and rows 3 / 6 would need a probe write the queued module does not have. No queued window module is edited before it runs (option (b) declined; option (c) not queued).
- **H32 (lane 1, no solve; serial on H30):** known-issues 2026-10-11 retired or re-titled on the printed ΣX_s / P_acc against 8.3e-3 — the entry resolve clause executed; on the re-title branch §2.2 says what the excess is (the sheet-field variance, X_s = CS_s) and that it does not shrink with order. The entry-42 rehearsal (2.149715e-02) is the likely branch, not a reading.
- **H31** gains: the `ans:8` `COMPARISON.md` / `metrics.json` v3 block takes the degree-2 rows as insertions only.

### §9 item 582 — 2026-10-22 (§10 H30, chain step (no solve); lane 1)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** §10 chain H, "H30 (lane 1, standard `-n 1`, no solve; serial on entry 41). As narrowed tonight (rulings): W1 with `FEM_EM_ANS8_H30=1` reading entry 41's `h9p-node-115.5209-d2.json` — rows 4–5 per accepted watt as SPEC v4 defines it, X_s per feed sheet, ΣX_s / P_acc; target: (Lin) and (P_acc-1) at their imported bands, the rows printed beside the degree-1 column; rows 3 / 6 / CS_s printed "not in record"; H21's branch stated in the outcome" (`PROJECT_PLAN.md:4879–4885`), and the 10-21 weekly's ruling "H30 — narrowed" (`PROJECT_PLAN.md:4775–4786`; this file, "Weekly ruling 2026-10-21"). The window is the one item 576 recorded: "W1 with `FEM_EM_ANS8_H30=1` in place of `FEM_EM_ANS8_H30A=1` and the ID `ANS-8-H30-W1`" (this file, item 576 "Windows"). Entry 41 has run (`20261022T070008Z_WF-10-step4-H28.log`, Status 0, 20 passed / 2 skipped on all 16 ranks; the review's reading). **No band, tolerance, Done-when, validation target or retirement condition is set or moved.** The 8.3e-3 comparison is not asserted anywhere; the outcome states which branch the printed number falls on, and §10 H32 (item 584) executes it.

- **Executor:** implementer.
- **Tier:** standard, `-n 1`, no `mpiexec` (NumPy / JSON / regex; no collective, no rank-local quantity — item 516's and item 576's precedent).
- **Class:** chain step (no solve). **Cost:** ≈ 8.0 measured (chain no-solve median) / ≈ 22 predicted slot-minutes: 15 fixed, W1 ≈ 2 (H30a's W1 took 100 s harness, pytest 99.02 s, `20261021T004407Z_ANS-8-H30a-W1.log:221`; this window parses one log twice instead of four logs), outcome / commit / journal ≈ 5. No code.
- **Depends on:** nothing (entry 41 has run). **Daytime:** no. **Slot-first:** no.
- **Family:** `ANS-8` step 2's degree-2 column; attempt 1 of 4.

**What exists (host, 2026-10-22; every input opened read-only).**
- **The module, as landed by item 576 (`87750550`), unedited since** (`git log` on the module and on the `ans:8` script shows `87750550` as the last commit for both): `tests/validation/test_ans8_v3_record_rows.py` (401 lines).
  - `test_ans8_h30_degree2_root` `:371–401`, gated on `H30_KNOB = "FEM_EM_ANS8_H30"` (`:62`) through `_need` (`:107–109`; unset ⇒ skip, `"1"` ⇒ run, anything else ⇒ `ValueError`).
  - It globs `H30_LOG_GLOB = "20261022T*_WF-10-step4-H28.log"` (`:71`) under `LOG_DIR = "docs/testing/logs"` (`tests/validation/test_wf10_lattice_d2.py:43`) and asserts exactly one match (`:378–379`). **Exactly one match exists:** `docs/testing/logs/20261022T070008Z_WF-10-step4-H28.log` (a listing of `docs/testing/logs/` for `20261022T*_WF-10-step4-H28.log`, 2026-10-22).
  - Record `_node_json_path(2, c_pf=C_IN_F_HUMAN_D2_PF)` (`:381`; `C_IN_F_HUMAN_D2_PF = 115.5209`, `tests/validation/test_wf10_f_human_tuned_b1_sar.py:78`), i.e. `/workspace/logs/wf10/h9p-node-115.5209-d2.json`.
  - It runs `_check_record` (`:126–183`) on that record with the entry-41 log as both node log and lattice log, prints the rows through `_print_rows("H30", …)` (`:186–218`), prints the degree-1 column from `V3_ROWS_JSON` (`:63`, `:388–392`), and writes `H30_JSON = /workspace/logs/ans8/v3-rows-d2-115.5209.json` (`:64`, `:393–401`).
- **Script functions** (`examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py`): `record_port_state` `:325–336`, `v3_record_rows` `:339–409` (asserts `r_ohm == ESR_OHM` and `feeds == order[:2]`, `:354–356`), `node_log_readings` `:415–441` (last `wrote …/<basename>` line, last ccw `(P_w)` line before it, every `(Lin)` line).
- **Entry 41's record:** `logs/wf10/h9p-node-115.5209-d2.json`, 5 486 B, mtime 2026-10-22 02:16:37 CDT (07:16:37Z, inside entry 41's 07:00:08Z + 1075 s window). Keys as the review's reading lists; `degree` 2, `node_knob` `"d2root"`, `c_in_pf` 115.5209, `r_ohm` 0.2, `feeds` `["P17", "P21"]` = `order[:2]`, `w["ccw"][1]` = `[6.123233995736766e-17, -1.0]`. **Its `partitions.ccw.rel` is 7.555134057306291e-13, which is the 10 g module's P_w line (`:13521`, "rel 7.555e-13")**, not the b1 module's (`:6831`, 7.441e-13) and not entry 42's (`20261019T070008Z_MAT-15-step6-P22-d2root.log:6786`, 7.324e-13): the file on disk is entry 41's second write (drafter's reading, not asserted by the test).
- **Entry 41's log, the lines the test reads:**
  - `wrote /workspace/logs/wf10/h9p-node-115.5209-d2.json` at `:6835` and `:13525` (two writes; the test takes the last);
  - ccw `(P_w)` `:6831` and `:13521` (the last before `:13525` is `:13521`); cw `:6833`, `:13523`;
  - `(Lin)` `:6832`, `:6834`, `:13522`, `:13524`: currents 1.309e-15 / 1.241e-15 / 1.118e-15 / 9.667e-16, Y_in route 1.406e-15 / 1.321e-15 / 1.714e-15 / 1.290e-15;
  - lattice: header `:13567` ("ccw 10 g table (18 points)", one header, not REDUCED), rows `:13568–13585` (18 `[WF-10 step2s]   r=` rows), c\* `:13586` (`r=45 dz=+0 az17+0`, SAR10g(c\*) 4.547175546e-03 W/kg), (iii-a) `:13587`, (iii-b) `:13588`, spreads `:13589`, (Sigma_P) `:13592`, phantom bbox with mirror plane `:13527` — every line `parse_table` (`test_wf10_lattice_d2.py:102–170`) needs.
- **The degree-1 column:** `logs/ans8/v3-rows.json` (9 696 B, 2026-10-11 12:09; item 386's V3 output); it carries `gap.rel_line_cut` (= 8.316822392098522e-03, file line `:71`), `row4.sar10g_per_w`, `row5`, `power_w.p_acc_spec`. The degree-1 X_s per feed are not in the H30 test's d1 line; they are printed at `20261021T004407Z_ANS-8-H30a-W1.log:83–84` (1.329324349e-04 / 1.213651853e-04 W), ΣX_s / P_acc at `:85` (8.316822e-03), row 4 at `:87` (0.1587143 W/kg/W).
- `logs/ans8/v3-rows-d2-115.5209.json` does **not** exist yet (listing of `logs/ans8/`: `s32.npz`, `v3-rows.json`).

**Code.** None. No edit to the module, the `ans:8` script, any `test_wf10_*` module or any `src/` file — **rule (i): the window runs the module as committed at `87750550`**, and no module edit is licensed before it. No W0: the module is unchanged since H30a's green W0 / W1, and its ungated `test_ans8_h30a_coding` runs in this window anyway.

**Anchors (asserted by the test as landed — re-read from the code; every band imported).**
- **Provenance of the log:** exactly one `20261022T*_WF-10-step4-H28.log` (`:379`); `MPI.COMM_WORLD.size == 1` (`:376`).
- **Record identity:** `node_knob == "d2root"` and `degree == 2` (`:383`); then in `_check_record` the tuple (degree, c_in_pf, node_knob) == (2, 115.5209, "d2root") (`:136`).
- **(Prov):** the log's last `wrote` path equals `_node_json_path(2, c_pf=115.5209)` (`:133`); the last ccw P_w line's P_src / sumP_sheet / P_phantom / P_surf against the JSON's `partitions["ccw"]`, each ≤ `REUSE_REPRODUCTION_RTOL` (`:143`; `tests/validation/test_port_birdcage_ring_matrix.py:171`, = 1e-6).
- **(Rec-Y):** the JSON's `y_in["ccw"]` against `_active_y_in(s_node, w_ccw, z0)` ≤ `Y0_CODING_BAND` (`:152`; `tests/validation/test_wf10_f_human_copper_esr_root.py:55`, = 1e-12).
- **(Lin):** every `(Lin)` reading in the log — four here — currents ≤ `UNIT_WEIGHT_CONTROL_RTOL` (`tests/validation/test_port_drive_superposition.py:303`, = 1e-12) and Y_in route ≤ `Y0_CODING_BAND` (`:154–158`). The field side of the target, re-asserted from entry 41's own log.
- **(P_acc-1):** the record-route line cut against the wave form ≤ `DISCRETE_IDENTITY_RTOL` (`:164`; `tests/validation/test_birdcage_power_identity.py:227`, = 1e-6). The S ↔ V / I algebra; label it so.
- **(Prov-L):** `parse_table` on entry 41's log succeeds (one header, 18 rows, `not lat["reduced"]`, `:180`), and its ccw P_phantom (the last ccw `w =` line before the table, `:13521`) matches the JSON ≤ `REUSE_REPRODUCTION_RTOL` (`:181`).
- **Output name:** `v3-rows-d2-115.5209.json` ≠ `V3_JSON_NAME` (`:397`; `08_…py:112`).
- **Printed, never asserted (the readings H30 exists for):** P_acc,SPEC and the partition terms; row 5 of P_acc and of the partition; X_s per feed sheet (P17, P21) with P_sheet,field and the lumped Norton term; the gap, ΣX_s, **ΣX_s / P_acc**; row 4 per W of P_acc and of the partition (from c\* `:13586`); "not in record: ['row3', 'row6', 'cs_s']"; the degree-1 column line (P_acc, row 5, row 4, ΣX_s / P_acc).

**Negative controls (rule (e)).**
- **Asserted, backed:** item 361's P_acc in (P_acc-1) misses by ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL` = 1e-5 (`:171`; `tests/validation/test_wf12_f_human_copper_matrix.py:61`). The same comparison on the same fixture (degree-2 root, entry 42's record) read 6.2316e-02 (`20261021T004407Z_ANS-8-H30a-W1.log:145`); predicted 6.2316e-02. Ceiling: the feed-ESR share of P_acc, ≈ 0.04–0.1, so the bar has three orders of headroom.
- **Printed only:** the ESR-blind miss (predicted 0.71896, `…H30a-W1.log:145`) and the feed-ESR-blind miss.
- **Not in this test:** the wrong-record and feed-label-swap controls live in H30a's two tests, which skip in this window. (Prov)'s power to tell a degree-2 record from another fixture's was shown by H30a (3.3445e-03, `…H30a-W1.log:154`). **(Prov) cannot tell entry 41's record from entry 42's at the printed digits**: the two windows' ccw P_w lines agree to every printed digit (`:13521` vs `20261019T070008Z_MAT-15-step6-P22-d2root.log:6786`). Provenance to entry 41 rests on the last `wrote` line, the mtime and the JSON's `rel` (above), which the outcome quotes.

**Predictions (written first; printed, never asserted).** The record should reproduce entry 42's (same builder, knob, mesh, degree, ranks; the P_w lines and the P17 / P21 sheet powers agree to every printed digit, `…H30a-W1.log:148–149` against the JSON's `p_sheet_by_id`):
- pytest 14 passed / 2 skipped (12 in `tests/environment` as in H30a's W1, plus `test_ans8_h30a_coding` and this test; H30a's two gated tests skip); harness ≈ 100 s.
- (Prov) ≈ 3.6e-11 / 3.8e-12 / 1.8e-11 / 3.7e-11 (drafter's arithmetic: JSON against the `.9e` prints at `:13521`); (Rec-Y) 0; (Lin) the four logged values above; (P_acc-1) ≤ 2e-15 (R42: 4.768e-16); (Prov-L) ≈ 1.8e-11, c\* `r=45 dz=+0 az17+0`.
- P_acc,SPEC 2.910591023476e-02 W; X_s 3.262132656e-04 / 2.994809404e-04 W; gap = ΣX_s 6.256942060e-04 W.
- **ΣX_s / P_acc 2.149715e-02 — the R42 rehearsal (`…H30a-W1.log:150`), a prediction, not H30's reading.** It is ≈ 2.58× the degree-1 8.316822e-03.
- Row 5 0.260954 / 0.020083 / 0.697465, sum 0.978503.
- Row 4 per W of P_acc ≈ 0.156229 W/kg/W, per W of the partition ≈ 0.159661 (drafter's arithmetic: entry 41's printed SAR10g(c\*) 4.547175546e-03 over R42's P_acc and partition). Inside H30a's printed window (0.145, 0.170); ≈ −1.57 % from the degree-1 0.1587143.
- NC 361 6.2316e-02; ESR-blind 0.71896.
- **H21's branch, predicted:** ΣX_s / P_acc ≥ 8.3e-3 ⇒ known-issues 2026-10-11 is re-titled (the excess does not shrink with order), not retired.

There is no `PREDICTED RED:` line: every assertion is predicted green, so a red anywhere stops the slot.

**Windows.** Foreground, Bash `timeout: 660000`, pytest `-s`, container `timeout -k 30`, `FEM_EM_ANS_SKIP_PRIVATE=1`. **The harness call is the whole Bash command.**
- **Host first:** Glob `logs/wf10/h9p-node-115.5209-d2.json` and `logs/ans8/v3-rows.json` (absent ⇒ BLOCKED, never re-solve), and Glob `docs/testing/logs/20261022T*_WF-10-step4-H28.log` (exactly one, else BLOCKED — the test would assert it anyway).
- **W1** (standard, `-n 1`; predicted ≈ 100 s harness):
  `scripts/testing/run_and_log.sh ANS-8-H30-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_ANS8_H30=1 timeout -k 30 180 python3 -m pytest tests/environment tests/validation/test_ans8_v3_record_rows.py -s -v --tb=short'"`
- No durable capture (well under ≈ 10 min). One window; it is the closing window.

**Traps already paid for.**
- **The provenance lines are labelled `[H30a R42]`.** `test_ans8_h30_degree2_root` calls `_check_record(mod, "R42", …)` (`:385`), so its (Prov) / (Rec-Y) / (Lin) / (P_acc-1) / NC / (Prov-L) lines print with the `R42` tag and R42's predictions beside the NC line. In this window they are **entry 41's record**. The rows print as `[H30a H30]` and the degree-1 column as `[H30]`. Do not edit the tag (rule (i)); the outcome says so.
- **The degree-1 column's "ΣX_s / P_acc" is `gap.rel_line_cut`** (`:392`), i.e. gap / P_acc from `v3-rows.json`. At degree 1 that equals ΣX_s / P_acc (ΣX_s = gap to ≈ 1e-10, item 576 outcome). Quote the degree-1 X_s per feed from `…H30a-W1.log:83–84`, not from this window.
- **Never set `FEM_EM_ANS8_H30A=1` in this window.** Entry 41 overwrote `h9p-node-115.5209-d2.json` at 02:16; H30a's R42 rehearsal would now read entry 41's JSON against entry 42's log. That is not item 576's rehearsal, and this item does not run it.
- **Never create another log matching `20261022T*_WF-10-step4-H28.log`.** The test asserts exactly one.
- **The builder rewrites the node JSON on every call**: provenance is to the **last** `wrote` line (`:13525`) and the last P_w before it (`:13521`); every (Lin) is asserted.
- **"ΣX_s vs p_361 − terminal" prints rel 2.000** — a sign by construction (item 576 outcome; item 386's own line, `20261011T170627Z_ANS-8-v3b-V3.log:6208`). "ΣX_s vs gap" prints ≈ 1e-10 (the partition residual amplified by the cancellation). Both are printed only, not reds.
- **Rows 3 / 6 and CS_s are not in the record.** Never fill them from the mean \|B₁⁺\| or the "centre purity" line.
- **Do not assert ΣX_s / P_acc against 8.3e-3** in code. State the branch in the outcome by reading the printed value against the literal 8.3e-3, with the degree-1 8.316822e-03 beside it. If the value lands in [8.3e-3, 8.316822e-03] (the rounding gap item 576 named, Q5), state that and do not choose a branch — the review rules.
- The test writes `/workspace/logs/ans8/v3-rows-d2-115.5209.json` (gitignored); never commit it, never write under `logs/wf10/`.
- Never run stage V3, never set `FEM_EM_WF10_DEGREE`; never import a `test_*` name; never `-k`; no pipe; `-s`; never delete a log; never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`.

**Scope.** A record-route reading of entry 41's record, no solve. **Not claimed:**
- rows 3 / 6 or CS_s at degree 2;
- the known-issues 2026-10-11 retirement or re-title (H32 executes it);
- `COMPARISON.md` / `metrics.json` rows (H31);
- any AED value, verdict or level;
- an h- or order-convergence statement, or entry 41's 10-24 weekly reading.

**Negative result.** Report and stop. Park on `attempt/ANS-8-H30-*` and mark item 582 BLOCKED in the same commit (rule (d)).
- **Glob ≠ 1 match:** a second H28 log exists. Name it; the review re-points.
- **(Prov) red:** the JSON is not the last write of entry 41's log. Name the mtime and the JSON's `rel`; never re-solve.
- **(Lin) red:** a contradiction with entry 41's green footer. Report it to the review (`log-pathologist`).
- **(P_acc-1) or (Rec-Y) red:** an algebra or record error. Stop.
- **(Prov-L) red / `parse_table` assertion:** the lattice block of entry 41's log is not the shape H27's parser reads. Print the header and row count.
- **NC 361 below 1e-5:** report the factor.
- **Any exception before the rows print** (a key or format mismatch): stop the brief there; no module edit in-slot.

**Status it can move.** No glyph; rule (j) does not bind (`ANS-8`'s §7 row, `PROJECT_PLAN.md:3583`, is ✅ on step 1; step 2 stays open on the row, its Done-when is the AED half).
- §10 gains an H30 DONE line with ΣX_s / P_acc and the branch.
- **H32 (item 584) becomes executable** on the printed value: known-issues 2026-10-11 retired or re-titled.
- H31 gains its input: the degree-2 rows for the `ans:8` v3 block, insertions only (`v3-rows-d2-115.5209.json`).

**For the review.**
1. **The `[H30a R42]` tag on entry 41's provenance lines** (`:385`) is a labelling quirk in the landed test, harmless to every assertion. Rule (i) forbids fixing it before the window. Whether a later step relabels it is named, not posed.
2. **The rounding gap [8.3e-3, 8.316822e-03]** in H21's threshold (item 576, Q5) stays the weekly's wording. The prediction (2.15e-02) is far from it.
3. **"Printed beside the degree-1 column":** the test's d1 line carries P_acc, row 5, row 4 and ΣX_s / P_acc, but not X_s per feed. The outcome quotes `…H30a-W1.log:83–84` for those. The drafter reads this as meeting the weekly's text without an edit.
4. **(Prov) and entries 41 / 42:** at the printed digits the two records are indistinguishable by P_w. The JSON's `rel` 7.555134057306291e-13 = `:13521` is the cheapest discriminator on record. The test does not assert it.

**Execution order.** Host Globs → W1 → outcome under this brief (every printed row, the d1 column, ΣX_s / P_acc and **H21's branch stated**, the `[H30a R42]` note) → §9 item 582 DONE → §10 H30 DONE line → commit (harness log, chunk outcome, PROJECT_PLAN lines) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 582 — `ANS-8` H30 (§10 chain H) — per `docs/planning/chunks/ANS-8.md` §9 item 582. Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Glob `logs/wf10/h9p-node-115.5209-d2.json`, `logs/ans8/v3-rows.json` and `docs/testing/logs/20261022T*_WF-10-step4-H28.log` (exactly one) first; if any check fails, stop and mark item 582 BLOCKED; never re-solve. Run W1 only, in the foreground, with `timeout: 660000` on the Bash call, pytest with -s, container `timeout -k 30`, `FEM_EM_ANS_SKIP_PRIVATE=1` and `FEM_EM_ANS8_H30=1`; never set `FEM_EM_ANS8_H30A`. Edit no code: not `tests/validation/test_ans8_v3_record_rows.py`, not the `ans:8` script, no `tests/validation/test_wf10_*` module, no `src/` file. Never run stage V3 and never set `FEM_EM_WF10_DEGREE`. Open every JSON read-only, never write under `logs/wf10/`, never commit `logs/ans8/v3-rows-d2-115.5209.json`. Do not assert ΣX_s / P_acc against 8.3e-3; state H21's branch in the outcome from the printed value. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/ANS-8-H30-*`, mark item 582 BLOCKED in the same commit, and report."

### Outcome 2026-10-22 (09:30 UTC slot; item 582 DONE)

W1 `20261022T093054Z_ANS-8-H30-W1.log` (standard, `-n 1`, no `mpiexec`; 14 passed / 2 skipped, pytest 98.16 s, harness 100 s, Status 0). Module run as committed at `87750550` (no edit, rule (i)). Host Globs: the node JSON (5 486 B, mtime 02:16:37 CDT = 07:16:37Z), `v3-rows.json` and exactly one `20261022T*_WF-10-step4-H28.log` present. Every assertion green; every prediction landed to its printed digit.
- **Label trap (review ruling (1)):** the provenance lines `:74–84` print with the tag `[H30a R42]` (`test_ans8_v3_record_rows.py:385` hard-codes it). In this window they are **entry 41's record**, not entry 42's. Not fixed here.
- **Provenance (review ruling (4)):** last `wrote` at `:13525` of entry 41's log (2 writes), JSON mtime 07:16:37Z (after `:13525`), JSON `partitions.ccw.rel` 7.555134057306291e-13 = the `:13521` P_w line. At the printed digits (Prov) cannot tell entry 41's record from entry 42's; provenance rests on those three facts.
- **Asserted:** (Prov) p_src 3.603e-11 / p_sheets 3.783e-12 / p_phantom 1.818e-11 / p_surf 3.730e-11 (band 1e-6, `:76`); (Rec-Y) 0 (`:77`); (Lin) currents 1.309e-15 / 1.241e-15 / 1.118e-15 / 9.667e-16, Y_in route 1.406e-15 / 1.321e-15 / 1.714e-15 / 1.290e-15 (band 1e-12, `:78–81`); (P_acc-1) 3.576e-16 (band 1e-6, `:82`; the S ↔ V / I algebra); NC 361 miss 6.2316e-02 ≥ 1e-5 (`:83`); (Prov-L) 18 rows, c\* `r=45 dz=+0 az17+0`, SAR10g 4.547175546e-03 W/kg, P_phantom rel 1.818e-11 (`:84`).
- **Printed rows, degree 2 at the root (`:85–92`)** beside the degree-1 column (`:93`; degree-1 X_s from `20261021T004407Z_ANS-8-H30a-W1.log:83–84`):

  | quantity | degree 2, 115.5209 pF | degree 1, 118.2437 pF |
  |---|---|---|
  | P_acc,SPEC (W) | 2.910591023475e-02 | 3.057629563251e-02 |
  | row 5 phantom / copper / ESR (of P_acc) | 0.260954 / 0.020083 / 0.697465, sum 0.978503 | 0.268112 / 0.019926 / 0.703644 |
  | X_s P17 / P21 (W) | 3.262132656e-04 / 2.994809404e-04 | 1.329324349e-04 / 1.213651853e-04 |
  | ΣX_s / P_acc | **2.149715e-02** | 8.316822e-03 |
  | row 4 per W of P_acc (W/kg/W) | 0.1562286 (of the partition 0.1596609) | 0.1587143 |

  Rows 3 / 6 / CS_s: "not in record" (`:92`). Printed-only lines as item 576 named them: ΣX_s vs gap rel 1.016e-10, ΣX_s vs p_361 − terminal rel 2.000 (sign by construction), ESR-blind 0.71896.
- **H21's branch:** the printed ΣX_s / P_acc 2.149715e-02 (`:89`) is above the literal 8.3e-3 and above the degree-1 8.316822e-03 (≈ 2.58×), outside the rounding gap ⇒ **re-title branch** (the excess does not shrink with order, two points). H32 (item 584) executes it; nothing asserted against 8.3e-3.
- `logs/ans8/v3-rows-d2-115.5209.json` written (gitignored, not committed); H31 Part B (item 597) is now executable.

### §9 item 584 — 2026-10-22 (§10 H32, chain step (no solve); lane 1)

*Drafted read-only for the 2026-10-22 03:00 daily review; the review rules on it.*

**Licence.** §10 chain H, "H32 (lane 1, no solve; serial on H30). Known-issues 2026-10-11 retired or re-titled on H30's printed ΣX_s / P_acc against 8.3e-3 — the entry's own resolve clause, executed; on the re-title branch §2.2's lumped-sheet line says what the excess is (the sheet-field variance, X_s = CS_s) and that it does not shrink with order" (`PROJECT_PLAN.md:4893–4897`), and the 10-21 weekly's "H21's trigger stands" (`PROJECT_PLAN.md:4783–4786`; this file, "Weekly ruling 2026-10-21"). **The resolve clause** is known-issues 2026-10-11's ruling (3) as re-posed: "The entry **retires if the excess shrinks with order (ΣX_s / P_acc < 8.3e-3 at degree 2)**; if it does not shrink, the entry is re-titled a finding on the lumped-sheet model. Nothing else in ruling (3) moves" (`docs/testing/known-issues.md:814`), with ruling (3)'s "chain F owns it" (`:812`). **No band, tolerance, Done-when, validation target or retirement condition is set or moved;** the re-title branch sets no new retirement condition (rule (6): the daily may not).

- **Executor:** implementer.
- **Tier:** none — documentation only, no window (see "Windows").
- **Class:** chain step (no solve). **Cost:** ≈ 8.0 measured (chain no-solve median) / ≈ 18 predicted slot-minutes: ≈ 10 fixed reading, ≈ 5 edits, ≈ 3 commit and journal.
- **Depends on:** item 582 (H30's printed ΣX_s / P_acc). Serial. **Daytime:** no. **Slot-first:** no.

**What exists.**
- **The entry:** `docs/testing/known-issues.md:804` heading ("## 2026-10-11 — 🟡 OPEN (`ANS-8` step 2 / H21, …) — on the tuned copper F-human at degree 1, the **field's feed-sheet power exceeds the lumped port's terminal accounting by 2.543e-4 W (8.3e-3 of the accepted power)** …"); body `:806–:814`: test id, what is measured, the reading, consequence, mechanism (CS_s / X_s = 1.000000000 on P17 / P21 at degree 1, `20261011T170627Z_ANS-8-v3b-V3.log:6204–6208`), next, ruling 2026-10-12 (1)–(3) (`:812`), the 10-18 check (`:813`), and the 10-19 re-pose (`:814`).
- **The file's retire convention** (`docs/testing/known-issues.md:727–735`): "Remove the entry in the same commit that fixes the test. If the retirement itself needs recording …, rewrite the entry's **heading** to begin `RETIRED YYYY-MM-DD` and leave the body: the weekly review moves every such entry, verbatim, to `known-issues-archive.md` and leaves one index line under 'Retired entries'". No test is fixed here, so the heading form applies. The rotation tool admits `## [✅ ]RETIRED YYYY-MM-DD` (`KI_ENTRY_RE`, `scripts/maintenance/rotate_plan_archive.py:116–117`). It treats as closed **any heading containing the uppercase words RETIRED, RESOLVED, FIXED or CLOSED** (`_ki_closed`, `:186`).
- **§2.2** (`PROJECT_PLAN.md:537–858`) has **no sentence on the feed-sheet excess today**: a grep of the span for "Cauchy", "excess", "X_s" or "known-issues 2026-10-11" returns nothing in it. Its lumped-sheet lines are the first bullet's S licence, "degree 2 on exact CAD rectangles, κ off, is the convention for every quoted lumped-sheet S figure" (`:575–580`), inside the bullet that lists the port-model findings to read "before quoting any S, `f₀` or Q entry" (`:539–615`, items (1)–(3)), and the port-model bullet `:617–750`. The F-small analogue is in §2.1 (`:216–222`, "its own Cauchy–Schwarz deficit … a fixed 1.06% sheet-field non-uniformity"; `:283–290`).
- **Readings to cite:** degree 1, ΣX_s / P_acc 8.316822e-03 and X_s 1.329324349e-04 / 1.213651853e-04 W (`20261021T004407Z_ANS-8-H30a-W1.log:83–85`); degree 2 at the root, **H30's harness log** (`<ts>_ANS-8-H30-W1.log`, item 582): the `[H30a H30] X_s P17` / `X_s P21` lines and the `[H30a H30] gap … ΣX_s / P_acc <v>` line.

**Edits (exactly these; nothing else).** First read `<v>` from H30's log line (as printed, `.6e`), never from the outcome prose. Compare it with the literal 8.3e-3.

*Re-title branch (`<v>` > 8.316822e-03; predicted).*
1. `docs/testing/known-issues.md:804` — replace the heading line with:
   `## 2026-10-11 — 🟡 OPEN (`ANS-8` step 2 / H21, §9 item 369; re-titled 2026-10-22 by §10 H32, §9 item 584, on §10 H30's reading; chain F owns it, weekly ruling 2026-10-12 (3)) — a finding on the lumped-sheet model: **the field's dissipation on a lumped feed sheet exceeds the port's lumped Norton term by the sheet-field variance X_s = CS_s, and the excess does not shrink with order** — ΣX_s / P_acc 8.316822e-03 at degree 1 (118.2437 pF), <v> at degree 2 (115.5209 pF), on the tuned copper F-human`
2. Append one bullet after `:814`:
   `- **Re-titled 2026-10-22 (§10 H32, §9 item 584; the re-pose's "does not shrink" branch, executed):** §10 H30 (§9 item 582) printed, at the degree-2 root from `xl` entry 41's record (`logs/wf10/h9p-node-115.5209-d2.json`, written at `20261022T070008Z_WF-10-step4-H28.log:13525`), X_s <x17> / <x21> W on P17 / P21 and **ΣX_s / P_acc <v>** (`<H30 log>:<lines>`, printed, never asserted) against 8.3e-3 — <v / 8.316822e-03, 3 s.f.>× the degree-1 8.316822e-03 (`20261021T004407Z_ANS-8-H30a-W1.log:85`). CS_s is not in entry 41's record; X_s = CS_s is the sheet model's own identity (CS_s / X_s = 1.000000000 at degree 1, `20261011T170627Z_ANS-8-v3b-V3.log:6204–6208`). Ruling (1)'s reading (a documented property of the discrete feed sheet, not of the solve, the ESR term or the port model) and ruling (2)'s normaliser stand; no assertion is added or changed. Chain F owns the finding (ruling (3)); this re-title sets no new retirement condition — that is the weekly's.`
3. `PROJECT_PLAN.md` §2.2 — insert, as continuation lines of the first bullet directly after `:615` (the end of item (3), "…no Q, `Z_in` or match figure is licensed.**"), wrapped to the bullet's two-space indent:
   `(4) **The lumped feed sheet's terminal accounting misses the field's sheet dissipation by the sheet-field variance, and the miss does not shrink with order** (known-issues 2026-10-11, re-titled 2026-10-22 by §10 H32 on H30's reading): on the tuned copper F-human each driven sheet dissipates X_s more in the field than its lumped Norton term ½ Re(1/Z_th)\|V\|², where X_s = CS_s, the Cauchy–Schwarz excess ½ Re(Y_s)(∫\|n̂×E\|² dA − \|∫E·ĥ dA\|²/A_s) of a field that varies across the sheet (CS_s / X_s = 1.000000000 at degree 1, `20261011T170627Z_ANS-8-v3b-V3.log:6204–6208`). ΣX_s / P_acc reads 8.316822e-03 at degree 1 and <v> at degree 2 at the root (`<H30 log>:<line>`). It is a property of the lumped-sheet model, not of the solve or the ESR term; every per-watt figure states its normaliser — the lumped line cut P_acc,SPEC (weekly ruling 2026-10-12), with the field partition printed beside it — and no assertion equates field and terminal feed power.`

*Retire branch (`<v>` < 8.3e-3).*
1. `docs/testing/known-issues.md:804` — replace the heading line with:
   `## ✅ RETIRED 2026-10-22 (§10 H32, §9 item 584, on §10 H30's reading) — was 🟡 OPEN 2026-10-11 (`ANS-8` step 2 / H21, §9 item 369) — on the tuned copper F-human at degree 1 the field's feed-sheet power exceeded the lumped port's terminal accounting by 2.543e-4 W (8.3e-3 of the accepted power); at degree 2 the excess shrinks, ΣX_s / P_acc <v> < 8.3e-3`
2. Append one bullet after `:814`:
   `- **Retired 2026-10-22 (§10 H32, §9 item 584; the re-pose's retire branch, executed):** §10 H30 (§9 item 582) printed ΣX_s / P_acc <v> at the degree-2 root (`<H30 log>:<line>`; X_s <x17> / <x21> W; entry 41's record), below 8.3e-3 and below the degree-1 8.316822e-03 (`20261021T004407Z_ANS-8-H30a-W1.log:85`): the excess shrinks with order, as ruling (3) said a sheet-discretisation property should. Ruling (2)'s normaliser stands for every per-watt row. The weekly moves this entry to `known-issues-archive.md`.`
3. No §2.2 edit, and do not run `rotate_plan_archive.py` (the weekly's).

*Both branches:* the outcome under this brief, §9 item 584 DONE, a §10 H32 DONE line naming the branch. Nothing else: not the dashboard, not §2.1, not `COMPARISON.md` / `metrics.json` / `SPEC.md` (H31), not ruling (3)'s other clauses.

**Anchor.** The entry's own resolve clause, executed on a printed number: `<v>` from H30's log against the literal 8.3e-3 (`known-issues.md:814`), with the degree-1 8.316822e-03 beside it. A docs step closes nothing under §4 and moves no glyph.

**Negative control.** None asserted (no window). The mechanical check that matters is the heading form: on the re-title branch the new heading still begins `## 2026-10-11` (so `KI_ENTRY_RE` sees an entry) and contains none of the uppercase words RETIRED / RESOLVED / FIXED / CLOSED (so `_ki_closed` does not rotate an open entry out). On the retire branch it begins `## ✅ RETIRED 2026-10-22`.

**Prediction.** The re-title branch: H30 is predicted to print ΣX_s / P_acc ≈ 2.149715e-02 (the R42 rehearsal, `20261021T004407Z_ANS-8-H30a-W1.log:150` — a prediction, not a reading), ≈ 2.58× the degree-1 value.

**Windows.** None. The `docrefs` checker (`scripts/testing/check_example_doc_references.py --docs-root docs/guide`) does not scan `known-issues.md` or `PROJECT_PLAN.md`, and no `docs/guide` page cites this entry (grep, 2026-10-22). Two tracked example files name it: `ans:8`'s `SPEC.md:291` ("known-issues 2026-10-11") and `COMPARISON.md:95` ("feed-sheet Cauchy–Schwarz excess (documented sheet property …"). Both stay true on either branch, and neither is edited here (`COMPARISON.md` is H31's; `SPEC.md` is the weekly's). `tests/unit/test_rotate_known_issues_entries.py` reads the pinned `KI_PINNED_SHA` copy once its four headings are gone from the live file (`:102–107`; they are gone), so a heading edit cannot move it.

**Traps already paid for.**
- **`_ki_closed` matches the uppercase words anywhere in the heading** (`rotate_plan_archive.py:186`). A re-titled heading that says "not RETIRED" would be rotated out of the live file as closed by the next weekly. Use none of the four words on the re-title branch.
- **Entries leave only with the commit that fixes them** (`known-issues.md:727`). A re-title keeps the entry 🟡 OPEN in place; only the retire branch rewrites the heading to `RETIRED`, and the weekly, not this item, moves it.
- **The 8.3e-3 threshold sits below the degree-1 8.316822e-03** (item 576, Q5). If `<v>` lands in [8.3e-3, 8.316822e-03], neither branch is unambiguous ("not < 8.3e-3" yet "shrinks"): stop and mark item 584 BLOCKED for the review.
- **Copy `<v>` and the X_s values from H30's log as printed**; never round them again, never take them from the outcome prose.
- In tracked Markdown write `≈`, never `~`; escape `|` as `\|` inside running text, as the entry already does.
- No AED number anywhere; never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`.

**Scope.** Not claimed:
- any new retirement condition for the entry (the weekly's);
- chain F's next step on the finding;
- rows 3 / 6 or CS_s at degree 2 (X_s = CS_s is quoted as the sheet model's identity, measured at degree 1);
- §2.1's human-scale line, `COMPARISON.md` / `metrics.json` (H31);
- any AED value or verdict.

**Negative result.**
- **Item 582 not DONE, or its log has no `[H30a H30] … ΣX_s / P_acc` line:** do nothing; mark item 584 BLOCKED "on item 582".
- **`<v>` in the rounding gap:** as above, BLOCKED for the review.
- No `attempt/*` branch is needed for a docs item that has not edited anything.

**Status it can move.**
- Known-issues 2026-10-11: re-titled and kept 🟡 OPEN as a finding on the lumped-sheet model (predicted), or `RETIRED` for the weekly's rotation.
- §2.2: gains clause (4) on the re-title branch.
- §10: an H32 DONE line. No §7 glyph; rule (j) does not bind.

**For the review.**
1. **"§2.2's lumped-sheet line" does not exist yet.** The drafter places a new clause (4) at the end of the first §2.2 bullet (after `:615`), the bullet of port-model findings to read before quoting a figure. The alternatives are a sentence beside the ANS-2 SAR bullet's "normalised to accepted power … the residual is the lumped-port feed systematic" (`:811–816`), or the port-model bullet (`:617–750`). The review picks one.
2. **"Chain F owns it"** is ruling (3)'s word (`known-issues.md:812`). Whether chain F is live to own it, and what retires the re-titled entry, are the weekly's.
3. **"Does not shrink with order" is two points** at two tuned roots (degree 1 at 118.2437 pF, degree 2 at 115.5209 pF), on one mesh. The clause says only that.

**Execution order.** Read item 582's outcome and H30's log line → choose the branch → edits → outcome under this brief → §9 item 584 DONE → §10 H32 DONE line → commit (`known-issues.md`, `PROJECT_PLAN.md`, chunk file) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 584 — `ANS-8` H32 (§10 chain H) — per `docs/planning/chunks/ANS-8.md` §9 item 584. Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. This item runs no window. Read ΣX_s / P_acc and X_s per feed from H30's harness log (`*_ANS-8-H30-W1.log`, item 582) as printed; if item 582 is not DONE or the line is absent, mark item 584 BLOCKED and stop. Compare it with the literal 8.3e-3. If it lies in [8.3e-3, 8.316822e-03], mark item 584 BLOCKED for the review and stop. Otherwise make exactly the edits of the branch it selects, in `docs/testing/known-issues.md` (the 2026-10-11 entry's heading and one appended bullet) and, on the re-title branch only, `PROJECT_PLAN.md` §2.2 clause (4). Add the §9 DONE and §10 H32 lines and the outcome, and edit nothing else. Never put RETIRED, RESOLVED, FIXED or CLOSED in a re-titled heading; never run `rotate_plan_archive.py`. Never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`; no AED number anywhere. You are the executor — do not spawn agents."

### Outcome 2026-10-22 (09:30 UTC slot, third item; item 584 DONE)

No window (documentation only, as the brief prices it). Item 582 DONE; H30's log `20261022T093054Z_ANS-8-H30-W1.log` (Status 0) prints `X_s P17: 3.262132656e-04 W` (`:87`), `X_s P21: 2.994809404e-04 W` (`:88`) and `ΣX_s / P_acc 2.149715e-02` (`:89`). Against the literal 8.3e-3 and the degree-1 8.316822e-03 (`20261021T004407Z_ANS-8-H30a-W1.log:85`) it is 2.58× higher, outside the rounding gap ⇒ **re-title branch**.
- **Known-issues 2026-10-11:** heading re-titled (still begins `## 2026-10-11 — 🟡 OPEN`, contains none of RETIRED / RESOLVED / FIXED / CLOSED — checked by grep, 0 matches) and one bullet appended after the 10-19 re-pose. **Review ruling 2026-10-22 03:00 applied over the brief's text:** the brief's "chain F owns it" became "ruling (3) names chain F; ownership of the re-titled finding is the weekly's", and the entry stays under `ANS-8` / known-issues; "does not shrink with order" is written as two points (degree 1, degree 2; one mesh), not a trend.
- **§2.2:** clause (4) inserted as continuation lines of the first bullet directly after item (3)'s "…`Z_in` or match figure is licensed.**" (the brief's `:615` is `:624` today; same place). It carries the same two-point wording.
- **§10 H32 DONE line** and **§9 item 584 DONE**. No band, tolerance, Done-when, validation target or retirement condition set or moved; no new retirement condition for the re-titled entry (the weekly's). `rotate_plan_archive.py` not run.

### §9 item 597 — 2026-10-22 (§10 H31 Part B, chain step (no solve); lane 1)

*Drafted read-only for the 2026-10-22 03:00 daily review as Part B of item 583 (`WF-10.md`); moved here verbatim by the review, which numbers it 597. Serial on item 582.*

**Outcome 2026-10-22 (09:30 UTC slot, fourth item): DONE.** Dependencies globbed first: the H30 JSON present, one `*_ANS-8-H30-W1.log` with `Status: 0`. Code as briefed: `_v3_d2_section` in `08_…py` and the call at the end of `comparison_text`; new `tests/unit/test_ans8_v3_d2_block.py`; gitignored helper `logs/ans8_h31_insert_d2.py`. W-ins `20261022T094349Z_ANS-8-H31-insert.log` (2 s): key absent before, (Deg) 2 / 115.5209 / d2root, (P_acc-1) r1 3.576e-16 ≤ 1e-6, both files re-rendered; host `git diff --numstat` reads `20 0` (`COMPARISON.md`) and `94 0` (`metrics.json`), so insertions only. W1 `20261022T094610Z_ANS-8-H31-W1.log` (standard, `-n 1`, 97 s, 39 passed / 3 skipped — the three knob-gated H30 validation tests): (R2) byte-equal at 9 743 B; (Ins) prefix = `e5d287d4`'s 8 250 B, remainder begins `"\n## Step 2 at the degree-2 root"`, metrics less `rows_v3_d2` deep-equals `8f50751e`; (Dig) 9 / 9 degree-2 numbers found in the 8 `[H30a H30]` lines, and the degree-1 cells equal `rows_v3`'s floats in the same format (degree-1 ΣX_s / P_acc is computed as Σ x_w / p_acc_spec from `rows_v3.gap`, 8.316822e-03); item 411's (R) / (P) / (D) green. Negative controls: N-absent equals `e5d287d4`'s bytes; N-digit breaks (R2); N-d1 misses exactly `0.1587143`. Wdoc `20261022T094756Z_ANS-8-H31-Wdoc.log`: `dead=0 guide=0 stale=125 exit=2` (the standing state). **Disclosed deviations:** (a) the renderer takes `(d2, v3)`, not `(d2)` alone, so it can put the degree-1 value beside each row, and the block renders only when `rows_v3` is present. The first W1 (`20261022T094401Z_ANS-8-H31-W1.log`, exit 1) failed `test_ans8_v3_rows.py::test_v3_comparison_text` with `KeyError: 'gap'`: that test's synthetic `rows_v3` has no `gap`, and its prefix check pops `rows_v3`. The fix was in the renderer (`.get` guards; the block goes only after a v3 section); the test was not edited. W1 re-ran green after the edit, and (R2) shows the committed bytes unchanged. (b) The brief's (Dig) is read as covering the degree-2 cells. The degree-1 beside-values are checked against `rows_v3` instead, because they are not `[H30a H30]` digits (N-d1 depends on that). Not claimed: any AED value, verdict or level; rows 1–3 / 6 / CS_s at degree 2; an order-convergence statement.


*Outline by the same drafter. The review assigns the number and moves the text.*

**`ANS-8` H31b: `COMPARISON.md` / `metrics.json` gain H30's degree-2 rows at the root as a new block after the v3 section, insertions only, rendered by the tracked writer from the committed JSON.**
- **Executor:** implementer.
- **Tier:** standard, `-n 1`, complex build (the script imports dolfinx).
- **Class:** chain step (no solve). **Cost:** ≈ 8.0 measured / ≈ 35 predicted. The price comes from item 411 (W-regen seconds; W1 82.6 s, `20261013T094942Z_ANS-8-item411-W1.log`) and item 513 (W1 94 s harness).
- **Depends on: §9 item 582 (H30) DONE.** Host Glob first:
  - `logs/ans8/v3-rows-d2-115.5209.json`;
  - exactly one `docs/testing/logs/*_ANS-8-H30-W1.log` with `Status: 0`.

  If either is absent, mark the item BLOCKED and never re-run H30.

**What exists.**
- **H30's JSON** is `test_ans8_h30_degree2_root`'s `out` (`tests/validation/test_ans8_v3_record_rows.py:393–400`). Its keys are those of `v3_record_rows` (`08_…py:386–408`) less `v_port` / `i_feed` as complex, which are re-added as pairs, plus `node_log`:
  - `power_w`, `row5`, `row5_partition`, `gap` (with `line_cut_minus_partition_w`, `mechanism[pid]`, `terminal_w`, `norton_lumped_w`);
  - `sum_x_w`, `sum_x_over_p_acc`, `row4` (when the lattice log is read; it is for entry 41);
  - `not_in_record`, `degree`, `c_in_pf`, `node_knob`.
  - **It has no `rel_line_cut`, no `cs_over_x`, no `row1_y_in_ccw_s`, no `row2_field`, no `row3`, no `row6_purity_o` and no `cells`.**
- **H30's printed digits** are the `[H30a H30] …` lines and the `[H30] degree-1 column …` line. Their formats are at `test_ans8_v3_record_rows.py:189–218` and `:389–392`: P_acc `.12e`, powers `.9e`, row 5 `.6f`, X_s `.9e`, ΣX_s / P_acc `.6e`, row 4 `.7f`, SAR10g `.9e`.
- **The committed writer** is `comparison_text` / `_v3_section` / `_v3_p_acc_definition` (`08_…py:1163–1290`) with `write_outputs` (`:1315–1318`).
  - `COMPARISON.md` was last changed at `e5d287d4`.
  - `metrics.json` was last changed at `8f50751e`. Its top-level keys run `… rows_v2, rows_v3, s32_reciprocity, …`, so a `rows_v3_d2` key sorts between `rows_v3` and `s32_reciprocity` and serialises as a pure insertion.

**Code.**
- **In `08_…py`, additive only:**
  - a pure `_v3_d2_section(d2)` renderer;
  - three lines at the end of `comparison_text`: `d2 = metrics.get("rows_v3_d2")` and, if present, `lines += [""] + _v3_d2_section(d2)`.
  - With the key absent, the render is byte-identical to today's.
- **The heading must not begin `## Step 2 (SPEC v3)`.** Proposed: `## Step 2 at the degree-2 root (SPEC v4: rows 4–5 from the node record, degree 2; ours; AED column blank)`.
- **The table**, with every row ending `| |`:
  - P_acc,SPEC;
  - row 4 per W of P_acc;
  - row 5 phantom / copper / ESR and the sum;
  - X_s per feed sheet;
  - ΣX_s / P_acc;
  - each **in the H30 log's own format spec**, with the degree-1 value from `rows_v3` beside it in the same cell, labelled "(degree 1 record: …)".
  - One line states that rows 1–3, row 6 and CS_s are "not in the node record at degree 2 (§10 H30 as narrowed by the 2026-10-21 weekly)".
- **No edit to** `_v3_section`, `_v3_p_acc_definition`, any stage, `main()`, `SPEC.md`, the private writer, or any existing test.
- **New `tests/unit/test_ans8_v3_d2_block.py`**, loaded by `importlib` as `test_ans8_v3_definition.py` does.

**Windows.**
1. Code.
2. **W-ins.** A gitignored helper `logs/ans8_h31_insert_d2.py` loads the committed `metrics.json`. It:
   - asserts `"rows_v3_d2" not in metrics`;
   - loads H30's JSON and asserts degree 2, `c_in_pf == C_IN_F_HUMAN_D2_PF` (imported) and `node_knob == "d2root"`;
   - asserts `check_v3_power(d2, DISCRETE_IDENTITY_RTOL)` (stage R's own re-check, imported);
   - sets the key and calls `write_outputs`;
   - prints the sha256 of both files before and after.

   Command: `scripts/testing/run_and_log.sh ANS-8-H31-insert "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 logs/ans8_h31_insert_d2.py'"`
3. **Host check:** `git diff --numstat` on the case directory reads `N 0` (no deletions) for both files.
4. **W1:** `scripts/testing/run_and_log.sh ANS-8-H31-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 timeout -k 30 180 python3 -m pytest tests/environment tests/unit/test_ring_mode_rows.py tests/unit/test_ans8_private_comparison.py tests/unit/test_ans8_v3_rows.py tests/unit/test_ans8_v3_definition.py tests/unit/test_ans8_v3_d2_block.py tests/validation/test_ans8_v3_record_rows.py -s -v --tb=short'"`
5. **Wdoc:** docrefs as in Part A (dead = 0, exit ≠ 1).

**Anchors (asserted).**
- **(R2)** `comparison_text(committed metrics.json)` equals the committed `COMPARISON.md`, byte for byte.
- **(Ins)** The committed `COMPARISON.md` starts with `git show e5d287d4:<case>/COMPARISON.md` (the `git -c safe.directory` pattern, `test_ans8_v3_definition.py:252–256`), and the remainder begins `"\n## Step 2 at the degree-2 root"`. `json.loads(committed metrics)` with `rows_v3_d2` popped deep-equals `json.loads(git show 8f50751e:<case>/metrics.json)`.
- **(P_acc-1)** `check_v3_power(metrics["rows_v3_d2"], DISCRETE_IDENTITY_RTOL)` (`tests/validation/test_birdcage_power_identity.py:227`).
- **(Deg)** The block has degree 2, `c_in_pf` = `C_IN_F_HUMAN_D2_PF` and `node_knob` "d2root".
- **(Dig)** Every number the block renders occurs verbatim in a `[H30a H30]` line of the tracked H30 log.
- **Item 411's existing (R) / (P) / (D) stay green in W1.**
  - (D) reads from the v3 heading to the end, so the new block must contain no "0.83 %" and no "accepted power at the 50 Ω lumped port", and every `| (` row in it must end `| |`.
  - `test_ans8_v3_definition.py:101–116`. Never edit that test to pass.

**Negative controls (asserted, by construction).**
- **(N-absent)** The render without `rows_v3_d2` equals `e5d287d4`'s bytes.
- **(N-digit)** A planted digit in the block breaks (R2).
- **(N-d1)** The block rendered with `rows_v3`'s degree-1 row 4 (0.1587143) in place of H30's fails (Dig).

**Prediction (printed only).** If entry 41's record reproduces entry 42's at `REUSE_REPRODUCTION_RTOL`, as H30a predicted:
- P_acc ≈ 2.910591023476e-02 W;
- ΣX_s / P_acc ≈ 2.149715e-02;
- row 5 sum ≈ 0.978503;
- row 4 ≈ 4.547175546e-03 / 2.910591e-02 ≈ 0.15623 W/kg/W.

These are H30's numbers to print, not this item's.

**Traps.**
- Never run stage R, S or V3: stage R rewrites `metrics.json` and moved the v1 / v2 floats by up to 1.4e-10 (item 386).
- Never hand-edit `COMPARISON.md` or `metrics.json`; both are owned by `nobody`.
- Never route H30's JSON through `_v3_section` / `_v3_p_acc_definition`: there is no `cs_over_x` / `row3`, so it raises `KeyError`.
- Never name `v3-rows-d2-115.5209.json` or the node JSON in `COMPARISON.md` (docrefs; H29).
- The v3 section's "Not claimed (v3): … the degree-2 column" line stays (insertions only). The new block's own first sentence says it is the degree-2 column at the root, labelled per SPEC v4 (`SPEC.md:331–336`). This is flagged, not edited.
- `FEM_EM_ANS_SKIP_PRIVATE=1` in every window; never open `aed_results/`, `COMPARISON_private.md`, `aed/` or `docs/private/`. The leak hook's refusal is a stop.

**Scope.** Not claimed: any AED value, verdict or level; rows 1–3 / 6 / CS_s at degree 2; known-issues 2026-10-11's disposition (H32's); an order-convergence statement.

**Status it can move.** No glyph. `ANS-8`'s step-2 Done-when waits on the AED half. The v3 degree-2 column exists in `COMPARISON.md` for the operator's W2 session, which SPEC v4 names as "our adjudicated column" (`SPEC.md:331–333`), and §10 H31 Part B gains a DONE line.

**For the review (Part B).**
- **Rows 1 and 2 at degree 2.** SPEC v4 row 1 (Im Y_in at both feeds) is printed in entry 41's log (`:6839–6840`, the ccw = (1, −1j) readings). It is not in H30's JSON, and the weekly's narrowed H30 did not pose it.
- **Row 2 at degree 2** would come from H25's saved degree-2 matrix at the root, which is not posed.
- Whether the block carries row 1 from the log is the weekly's question, not a slot's.

### Review ruling 2026-10-22 03:00 (items 582 H30, 584 H32, 597 H31 Part B)

Binds over the briefs. **Entry 41 is read** (ledger row 2026-10-22, `log-pathologist` claims (a)–(g) CONFIRMED; `:6837`, `:6852`, `:13586`, `:13597`, `:14681–14686` re-read): every gate green at the degree-2 root, the rule read literally fires its **second branch** (SAR10g(c\*) −6.30 %, \|B₁⁺\| −3.17 % root to root) — the 10-24 weekly decides.
- **582 (H30): queued as drafted.** No module edit before W1 (rule (i)). (1) The hard-coded `R42` tag on the root test's provenance lines (`test_ans8_v3_record_rows.py:385`) is a **label trap**, named in the outcome, not fixed in this item. (2) **Rounding gap:** a printed ΣX_s / P_acc inside [8.3e-3, 8.316822e-03] decides neither branch; H32 then stops and returns the reading to the 10-24 weekly. (3) The degree-1 X_s per feed is quoted from `20261021T004407Z_ANS-8-H30a-W1.log`. (4) (Prov) cannot tell entry 41's record from entry 42's at printed digits; provenance rests on the last `wrote` line and the file's mtime (07:16Z, after entry 41's `:13525`), stated in the outcome.
- **584 (H32): queued as drafted, serial on 582.** The §2.2 clause goes after `:615` (the drafter's default). Chain ownership of a re-titled finding is the weekly's; the entry stays under `ANS-8` / known-issues. "Does not shrink with order" is written as a two-point statement (degree 1 and degree 2), never as a trend.
- **597 (H31 Part B): the split is accepted.** Serial on 582; if 582 has not landed, the slot skips it at no attempt cost. Insertions only on `COMPARISON.md` / `metrics.json`; degree-2 rows 1 / 2 stay the weekly's question.

### §9 item 623 — 2026-10-25 (§10 H35, lane 1, no solve)

*Written by the 2026-10-25 03:00 daily review from the 10-24 weekly's enumerated step H35 (§9 rule (4)). Two execution readings below are this review's and go to the 10-26 weekly; neither changes the target.*

**Question.** Item 597 left the `ans:8` degree-2 block without row 2, the feed-ring share and the k = 0 leak. Can it carry them from the circuit route on entry 36's saved degree-2 copper 32×32 (H25), labelled as that route?

**Execution readings (this review's, for the weekly):**
1. **Row numbering.** In `COMPARISON.md`, row (1) is Y_in and row (2) is the feed-ring share / k = 0 leak (`COMPARISON.md:79–82`). The weekly's "rows 1 / 2 (feed-ring share, k = 0 leak)" is therefore **row 2 only**. Row 1 (Im Y_in at both feeds, printed in entry 41's log `:6839–6840`) stays item 597's open question for the weekly. It is not inserted here.
2. **The condition of H25's values.** 0.968088 / 0.025825 are the **lossless** circuit-route values at H25's own `C_mode`(d2, Cu) = 114.5221 pF. H25 also prints the **ESR production convention** at the same capacitance, 0.948363 / 0.041953 (`20261018T095559Z_ANS-8-H25-W1.log:142–143`, `:165–166`). Neither value is at the 115.5209 pF in-model root, which the block's header names. The insertion carries **both**, each with its condition stated. That keeps the binding ESR ruling (production = 0.2 Ω; lossless = labelled record). H25 wrote no committed JSON, so the log lines are the source.

**Anchor.** A transcription identity: every inserted digit string is byte-equal to its H25 log line (`:142–143`, `:165–166`). `git diff --numstat` on the two files reads insertions only, 0 deletions (H31 Part B's precedent, `ANS-8.md:1319`). Docrefs reads 0 dead.

**What is inserted.**
- In `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/COMPARISON.md`'s degree-2 block (`:102–120`), after the existing rows, one row-2 line:
  - "row 2 (degree 2, **circuit route — Norton feeds on entry 36's saved 32×32 at C_mode(d2, Cu) = 114.5221 pF; in-model row not in record**): ccw feed-ring share 0.968088 / k = 0 leak 0.025825 (lossless); 0.948363 / 0.041953 (r = 0.2 Ω, the production convention)", with the log path.
  - Edit nothing on lines `:118` / `:120`. **Append**, beneath them, one sentence saying row 2 is now carried from the circuit route as labelled, and that rows 1, 3, 6 and CS_s stay not in record.
- In `metrics.json` under `rows_v3_d2` (`:3711`), one new key `row2_circuit_route_h25`. It holds `{"route": "circuit (Norton feeds)", "c_mode_pf": "114.5221", "lossless": {"ccw_feed_share": "0.968088", "ccw_feed_leak_k0": "0.025825"}, "esr_0p2_ohm": {"ccw_feed_share": "0.948363", "ccw_feed_leak_k0": "0.041953"}, "source": "docs/testing/logs/20261018T095559Z_ANS-8-H25-W1.log:142-143,165-166"}`. Store the values as strings copied from the log, and insert the key before `not_in_record`. The JSON must re-parse (`python3 -m json.tool` in the container, as its own harness window).

**Negative control.** A planted transcription error must be caught: the check script, run on a copy with one digit changed, reports a mismatch. Write it as `logs/ans8_h35_check.py`, gitignored scratch, run at `-n 1` in the harness. It greps each inserted string against the log line and exits 1 on any miss. W-c runs it on a copy with 0.968088 → 0.968089 and must exit 1. W1 runs it on the real files and must exit 0.

**Windows.**
1. W1 `ANS-8-H35-W1` (the check script).
2. W-c `ANS-8-H35-Wc` (the planted copy).
3. Wjson (`json.tool`).
4. Wdoc: `docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'` (the H31 Wdoc form, `20261022T094756Z_ANS-8-H31-Wdoc.log`).

Each is a bare harness call.

**Tier:** smoke / standard, `-n 1`, no MPI, no solve. **Class:** lane-1 revalidation. **Cost:** ≈ 11.4 measured / 20 predicted slot-minutes. **Executor:** implementer. **Depends on:** nothing. **Family:** H35, attempt 1 of 4.

**Traps.**
- No AED figure enters these files. The block is ours only, and the AED column is the operator's (W2). Do not open `COMPARISON_private.md`.
- Insertions only: check `git diff --numstat` before committing.
- The leak guard may flag a digit string. If it does, move the log to `logs/leak-held/` unchanged and report it (§9 digest); never move it to `docs/private/`.

**Scope boundary.** No new number is computed. Not closed here: row 1 at degree 2, the in-model row 2, `ANS-8` step 2 (W2). `ANS-8` stays ✅ on step 1.

**A negative result means:** a mismatch is a transcription defect. Fix the transcription and re-run W1; that is not a re-pose. If the leak guard blocks the commit, stop and report.

**Status it can move.** Item 597's open question closes for row 2 (§10 H35 → DONE). The `ans:8` degree-2 block gains row 2 as a labelled circuit-route record.

**Outcome 2026-10-25 (04:30 slot, item 623): DONE.** Inserted in `COMPARISON.md` (one table row after ΣX_s / P_acc, one sentence after the "Not claimed (degree 2)" line) and `metrics.json` (`rows_v3_d2.row2_circuit_route_h25`, before `not_in_record`): ccw feed-ring share 0.968088 / k = 0 leak 0.025825 (lossless) and 0.948363 / 0.041953 (r = 0.2 Ω, the production convention), circuit route at C_mode(d2, Cu) = 114.5221 pF, labelled; row 2 only, per the review's reading 1. Windows, all `-n 1`, smoke: W1 `20261025T094023Z_ANS-8-H35-W1.log` — the check script parses the four values from H25 `:142–143` / `:165–166` (and cross-checks the two lossless prints against each other, and `c_mode_pf` against `:171`'s 1.145221310828e-10 F) and finds every string byte-equal in both files, 0 misses, exit 0, 0 s. W-c `20261025T094028Z_ANS-8-H35-Wc.log` — on a copy with 0.968088 → 0.968089 the script reports 3 misses (COMPARISON value, COMPARISON phrase, metrics value) and exits 1: the negative control bites. Wjson `20261025T094033Z_ANS-8-H35-Wjson.log` — `json.tool` re-parses, exit 0. Wdoc `20261025T094034Z_ANS-8-H35-Wdoc.log` — dead=0 guide=0; exit 2 from the stale-artifact report only (stale=145, `--stale-severity report`), the same form as H31's Wdoc (`20261022T094756Z_ANS-8-H31-Wdoc.log:165`, stale=125, exit 2). `git diff --numstat`: COMPARISON.md 3 / 0, metrics.json 13 / 0. No AED figure entered; the leak guard flagged nothing. Not closed: row 1 at degree 2, the in-model row 2, step 2 (W2).

### §9 item 740 — 2026-11-01 (lane 1, known-issues 2026-11-01: `main`'s red node)

**Ruling (18:00 review): the renderer moves, and (R2) stays.** H35 hand-inserted three lines into `COMPARISON.md`, and `metrics.json` already holds their values under `rows_v3_d2.row2_circuit_route_h25`. (R2), render = committed bytes, is the identity that keeps the page generated. Re-posing it would remove a gate to save a renderer edit. No band moves.

**Anchor.** (R2) `comparison_text(metrics)` equals the committed `COMPARISON.md` byte for byte, in `test_ans8_v3_d2_block`, with every assertion below it running **for the first time since 2026-10-25**: (Ins), (P_acc-1), (Deg), (Dig), the degree-1 cells, the (D) compatibility lines. New (Row2): the four values and `c_mode_pf` on the rendered row equal `row2_circuit_route_h25`'s strings.

**Negative control.** (N-row) metrics with `row2_circuit_route_h25` popped render with neither the row nor the closing sentence, and differ from the committed page. (N-digit′) a planted digit in that key's `ccw_feed_share` renders a page ≠ the committed one. (N-absent) on the older metrics stays green, so the render must be conditional on the key.

**Execution.** In `examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/08_f_human_quadrature_mode_content.py`, `_v3_d2_section`:
- After the `ΣX_s / P_acc` row, when the key is present, the "row 2 (degree 2, **circuit route — …**)" row, values from the key.
- After the "Not claimed (degree 2)" line, a blank line and the "Row 2 is now carried at degree 2 …" sentence.
- The fixed wording is copied from the page; `git show 704025e6 -- <COMPARISON.md>` is the source.

**Licensed in `tests/unit/test_ans8_v3_d2_block.py`, and only this:** `_block_cells` skips a line starting `| row 2`, so `len(cells) == 9` and `n_tok == 9` stay as they are and that row's digits, which are not `[H30a H30]` numbers, are gated by (Row2) in their place. The two new controls are additive. No other landed line moves.

**Tier, ranks, wall clock.** Standard, `-n 2`, complex, `timeout -k 30 170`, the width of `20261101T182231Z_OPS-73-Wreg.log`; ≈ 5 s. Then that window's module list again, Wdoc, Wenv. **Class:** lane-1 revalidation. **Cost:** ≈ 6.0 measured / 20 predicted min. **Depends on:** nothing.

**Traps already paid for.** **The page's source reads `log:142–143,165–166` with en dashes, and the key's `source` string has hyphens:** render the page's form, and do not edit `metrics.json` or `COMPARISON.md`; both are byte-pinned records. The page's row ends `| |`, and the (D) line only checks rows starting `| (`, so it passes untouched. (Ins) pins the page's prefix at `e5d287d4`: the additions sit after it. `H30_LOG` at `:129` has not executed since `OPS-73` routed it (audit caveat), so a `.gz` or line failure there is a finding. AED numbers never enter this file: the block's columns stay empty.

**Scope.** The renderer and one parser line. **Not claimed:** anything about the row's values, which are H25's records; `ANS-8`'s status.

**Negative result.** Any assertion below (R2) red, other than the licensed line: stop, keep the log, BLOCKED, a new known-issues entry, because that red has been hiding behind (R2) since 10-25. No `PREDICTED RED:` line.

**Status it can move.** The known-issues entry retires in this commit; `main` has no known red node; the §9 digest's line goes.

**Outcome 2026-11-02T02:02Z (19:30 CST slot on 11-01, sixth item): DONE.** Run inline, ≈ 3 min, take-next at minute 29. `_v3_d2_section` renders the row and the closing sentence when `row2_circuit_route_h25` is present, from the key; `metrics.json` and `COMPARISON.md` are untouched. W1 `20261102T015938Z_ANS-8-item740-W1.log` (`-n 2`, complex, Status 0, 3 s): 6 passed per rank. **(R2) true, 10302 bytes,** and every assertion below it ran for the first time since 10-25, all green: (Row2) rendered and key strings equal on five values; (Ins) 8250-byte prefix and metrics equal; (P_acc-1) 3.576e-16; (Deg) degree 2, 115.5209 pF, `d2root`; (Dig) 9 rows, 9 numbers, none missing, on the `.gz` `H30_LOG` read through `read_log_text`; the degree-1 cells; the (D) lines. Controls: (N-absent) true, (N-digit), (N-d1), (N-row) neither row nor sentence, (N-digit′) unequal. Wreg `…015951Z_ANS-8-item740-Wreg.log`, the `OPS-73` window's module list, 10 passed, 4 skipped per rank, Status 0, its one red gone. Wdoc `…015955Z_ANS-8-item740-Wdoc.log` `dead=0 guide=0 stale=0`, Status 0. Wenv `…015957Z_ANS-8-item740-Wenv.log` 12 passed, Status 0. Beyond the brief, Wcal `…020053Z_ANS-8-item740-Wcal.log`, the four other modules that load this script (`tests/unit/test_ans8_v3_rows.py`, `tests/unit/test_ans8_v3_definition.py`, `tests/validation/test_ans8_degree2_rows.py`, `tests/validation/test_ans8_copper_degree2_rows.py`): 9 passed, 2 skipped per rank, Status 0.

**Disclosures.** The source's line range takes its en dashes in the renderer, split at the last colon so the log's name keeps its hyphens. The closing sentence reads `c_mode_pf` from the key and the root from `d2['c_in_pf']`. Wreg's skips are nodes needing gitignored `logs/wf10` records; they did not run. The known-issues entry is retired in this commit. The §9 digest line on `main`'s red node is left for the review, which owns the digest; `main` has no known red node from here.
