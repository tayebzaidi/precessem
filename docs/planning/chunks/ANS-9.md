# ANS-9 — history

Item briefs and execution history (operator directive 2026-10-01). The §7 row is the status of record; the SPEC is `examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md`.

### §9 item 218 — 2026-10-04 (chain P4)

*Drafted by a general-purpose agent for the 2026-10-04 03:00 daily review; the review reviews it.*

**`ANS-9` step 1a: the runnable half of the SPEC (`ans:9`).** Rows 1–7 come from `MAT-15` step 1's fixture, with every band imported. Row 6, the energy balance, is asserted at `MAT-15`'s `E_BAND`. The target is the weekly's (§10 P4); this brief adds execution detail only.
- **Executor:** implementer. The SPEC asks for a script dispatched by `./scripts/run_examples.sh -e ans:9`, not a test module. `example-runner` handles EX chunks only, and every `ANS-*` runnable half has run under the implementer (`ANS-8` item 141 is the precedent). `example-runner`'s template rules still apply: import, never restate; a same-stem guide; censuses predicted first.
- **Tier:** heavy (windows of 160–220 s), `-n 2`, complex. Use `-n 2` only, because the `P_k` digits and `G` are `-n 2`'s (`20261002T114023Z_EX-79.log:1878`).
- **Cost:** ≈ 18 measured (chain solve class median) / ≈ 50 predicted min. That is 20 min of code plus about 9 min of windows: the 206 s re-gate, the ≈ 190 s run, the pin test and the censuses.
- **Depends on:** nothing. `MAT-15` is ✅ and `EX-79` already exercises `_run_heated(return_fields=True)`.

**Why one fixture answers every row.** The SPEC's thermal half is `MAT-15` step 1's thermal half: k 0.6, ρc 4.18e6, h 10, T∞ = T(0) = 20, duration 1254 s = `T_STAR`, Δt 31.35 s = `T_STAR / N_MAIN`, and the source normalised to 1 W by `G = 1 W / P₁`. Its EM half is `ANS-2`'s 10 MHz fixture on the 0.0025 m rung, which is `MAT-15`'s. Drive k = 0 is P1 (`EX-79.log:1877`). So every row is drive k = 0 of `_run_heated`, plus three things the fixture does not already keep: snapshots at every step, the first half-step state, and the P1 excitation record.

**Code.**
- **Gate helper (additive return keys; rule (a)/(c), disclosed).** Under `return_fields=True`, `_run_heated` (`tests/validation/test_mat15_coil_driven_heating_phantom.py:180`) also returns `sweep` and `solves`, which are already local variables at `:190` and `:199`. With `False`, the prints and keys stay byte-identical. Re-gate the module in the same slot (W0).
- **Script:** `examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/09_coil_driven_gel_heating_10MHz.py`, plus the same-stem guide `09_coil_driven_gel_heating_10MHz.md`.
  - The `09_` prefix is what `ans:9` selects (`scripts/run_examples.sh:247–248`).
  - Put no other `*.py` in the case directory: the runner lists every depth-2 script there (`run_examples.sh:129`).
  - Copy the `_REPO_ROOT` insertion and the complex guard from `examples/materials/09_coil_driven_heating_phantom.py`.
- **Imports, never restated.** From the `MAT-15` module: `_run_heated`, `_transient`, `_be_one`, `K_GEL`, `H_GEL`, `T_INF`, `A`, `T_STAR`, `N_MAIN`, `DT_1`, `RHO_C`, `V_SUB_RECORD_M3`, `V_BAND`, `B_BAND`, `E_BAND`, `S2_BAND`, `D_CONTROL_FLOOR`, `PHANTOM_CELL_TAG`, `ROBIN_TAG`. Also `PHANTOM_HEIGHT` from `tests/mesh/test_birdcage_port_tags.py` and `evaluate_vector_field_parallel`.
  - Load row 7's convention from `ANS-2`'s own `_excitation` / `_incident_power_w` (`02_birdcage_coil_driven_sar_10MHz.py:256–307`) by `importlib`, the way that script loads its own siblings (`:93`). The file is `__main__`-guarded (`:1073`). Row 7 is `ANS-2`'s question, so it must be `ANS-2`'s arithmetic.
- **SPEC table, used only for a check.** A dict `SPEC_THERMAL` holds the SPEC's typed numbers: k, ρc, h, T∞, T(0), duration, Δt bound, a, the half-height, the four times, and the five probes. It is used **only** in (A0); every computation uses the imports.
- **Computation (drive k = 0).**
  1. Call `out = _run_heated(return_fields=True)`.
  2. Re-run `_transient(sub, ft, q_h[0], dt=T_STAR/N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, snapshot_steps=range(N_MAIN + 1))`. That is `_cn`'s schedule (`:154–157`) with every full step kept.
  3. Get θ^{½} from `_transient(…, dt=T_STAR/N_MAIN/2, n_steps=1, theta=1.0)`, the same operator as the first Rannacher half-step (`src/fem_em_solver/thermal/transient.py:323–326`).
  4. Get row 3 from `_be_one(sub, ft, q_h[0])`.
  5. The doubled-source and dropped-surface controls (below).
  6. Write `metrics.json`, `COMPARISON.md` and the XDMF.
- **Rows.**
  - **Row 1:** ∫_sub q_h dV.
  - **Row 2:** ⟨θ⟩ = ∫θ/V_sub at steps 10 / 20 / 30 / 40.
  - **Row 3:** the slope θ̄¹/Δt₁ in K/s from `_be_one` (one BE step of `DT_1`, `MAT-15`'s calorimetric step), printed beside 1 W/(ρc·V_sub).
  - **Row 4:** θ at A–E through `evaluate_vector_field_parallel(θ_40, pts, owner="lowest_original_cell")`, with the points identical on every rank.
  - **Row 5:** the maximum θ over **owned** dofs, with its coordinates, gathered by `allgather` and picked on every rank.
  - **Row 6:** the surface heat as an independent sum:
    `(Δt/2)·∮hθ^{½} + (Δt/2)·∮hθ¹ + Σ_{n=2}^{40} Δt·∮h·½(θⁿ + θⁿ⁻¹)`,
    which is the scheme's own θ-weighted `conv_form` (`transient.py:249–251`) per (half-)step.
  - **Row 7:** `G · _excitation(solves["P1"])["incident_power_w"]`, with supplied and accepted power beside it.
  - **Metadata:** solver and dolfinx version; EM N1curl degree 1 and thermal P1; 199 920 / 58 866 cells; the scheme ("CN θ = 0.5, Δt = 31.35 s, 40 steps, Rannacher start of 2 BE half-steps"); wall time.
- **Artifacts.**
  - `metrics.json` is tracked. Add it to `LISTED_EXAMPLE_ARTIFACTS` (`tests/unit/test_doc_reference_exit_codes.py:469`) and update the pin count.
  - `COMPARISON.md` carries the SPEC's row table: our column filled, the AED columns blank (`—`), and the SPEC's readout copied verbatim.
  - The combined XDMF `paraview_output/ans9_coil_driven_gel_heating_10mhz_combined.xdmf` (gitignored) is written by `write_xdmf_with_tags(…, sub, None, …)` + `adopt_host_ownership`. It holds P1 `theta_t0313`, `theta_t0627`, `theta_t0940`, `theta_t1254` (K) and DG0 `q_h` (W/m³).
  - **No private writer.** The adjudication files are the weekly's.

**Anchors (asserted; bands imported).**
- **(A0) Structure.**
  - Every `SPEC_THERMAL` entry equals its `MAT-15` import to 1e-12 relative: the four times are `T_STAR · {¼, ½, ¾, 1}` at steps 10 / 20 / 30 / 40; Δt = `T_STAR/N_MAIN` ≤ 31.35; a = `A`; the half-height = `PHANTOM_HEIGHT/2`.
  - `by_k[0] == "P1"`.
  - These `is` checks hold: `m15.E_BAND is m12.E_BAND`, `m15.S2_BAND is m12.S2_BAND`, `m15.B_BAND is m10.B_BAND` (EX-79's).
  - All five probes are valid.
- **(R1) Row 1:** |∫q_h − 1 W| ≤ `B_BAND`. Predicted 2.2e-14 (`EX-79.log:1881`). V_sub against `V_SUB_RECORD_M3` ≤ `V_BAND`, predicted 3.4e-13.
- **(R0) Re-run = gated run:** ‖θ_40 − `T_star[0]`‖/‖θ_star‖ ≤ `E_BAND`. Predicted ≤ 1e-14, because the operator and schedule are identical. Each new run's `max_balance_residual` is ≤ `E_BAND`.
- **(R3) Row 3:** |ρc·θ̄¹/(Δt₁·∫q/V_sub) − 1| ≤ `S2_BAND`. Predicted 5.7543e-6 (`EX-79.log:1884`). The slope is predicted at 1.058563e-3 K/s against 1/(ρcV) = 1.058569e-3 K/s.
- **(R6) Row 6, the energy balance:** |T·∫q − ρc∫θ_40 − row 6| / (T·∫q) ≤ `E_BAND`, with T = `T_STAR`.
  - Predicted ≤ 1e-12, which is about 42 steps × the ≤ 2e-14 per-step residual (`EX-79.log:1882`).
  - Row 6 itself is predicted at **164.50 J**. That is 1254 − ρcV_sub·⟨θ⟩(t*) = 1254 − 944.671 × 1.153314 (`:1882`), a surface fraction of 0.1312.

**Printed (predicted, never asserted).**
- **Row 2:** monotone, and below the no-loss bound t/(ρcV) = 0.3319 / 0.6637 / 0.9956 / 1.3274 K. ⟨θ⟩(1254) is 1.153314 K.
- **Row 4:** order B > C > D. P1 sits at azimuth 360° (+x), and the far-side ball is 75 % cooler (`:1902`). E is a record.
- **Row 5:** at r ∈ [0.020, 0.030] m, within ±30° of +x.
- **Row 7:** 45 428.08 W = 2.5e-3 W × G (`ANS-2` `metrics.json` `excitation/0/incident_power_w`; `EX-79.log:1878`). It is printed beside the tracked `metrics.json` ratio `incident_power_w / phantom_power/0/dissipated_power_w`, predicted to agree to ≈ 5e-12 (`-n 4` against `-n 2`).
- **The main run's own first half-step slope** θ̄^{½}/(Δt/2), the comparand closest to "the first time step": deficit predicted [3e-3, 7e-3] (κ_s·hAΔt/(2ρcV)). It is printed so the SPEC's 5 % readout sees our scheme's first step. It is not row 3.

**Negative controls (rule (e)).**
- **(d6) Doubled source**, the factor-of-two time-average error the SPEC exists to catch. The snapshot schedule re-runs on 2·q_h.
  - (R6) against T·∫q_h misses by **predicted 1.000**. The ceiling is exactly 1 by linearity.
  - Asserted only > `D_CONTROL_FLOOR`, `MAT-15`'s sign-only form (`:101–103`). The (E′) analogue measured 1.000000 at `EX-79.log:1901`, but that is a different comparison, so the factor is labelled predicted.
- **(s6) Surface term dropped.** (R6) with row 6 set to 0 misses by row 6/(T·∫q). **Predicted 0.1312.** Printed, and asserted only > `E_BAND`.

**Windows.** Each harness call is bare and whole. Run foreground with Bash `timeout: 660000`.
- **Pre-censuses:**
  - `scripts/testing/run_and_log.sh ANS-9-pre-figures "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_setup_figures.py'"`
  - `scripts/testing/run_and_log.sh ANS-9-pre-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
  - Run both with the new script **not yet** in the tree (`ANS-8`'s superseded-census trap).
- **W0, the re-gate after the additive keys** (predicted 206 s, `20261002T112729Z_EX-79-mat15-regate.log`):
  `scripts/testing/run_and_log.sh ANS-9-step1a-mat15-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 480 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat15_coil_driven_heating_phantom.py -s -v --tb=short'"`
- **Dry run:** `scripts/testing/run_and_log.sh ANS-9-dryrun "./scripts/run_examples.sh -e ans:9 -n 2 -t 480 --dry-run"`. It must print `(complex build)` and the `09_…py` path.
- **W1, the run** (the line the dry run emits, in `20261002T114312Z_EX-79-control.log:12`'s shape; predicted 170–220 s, from the 162 s of `EX-79` plus ≈ 10 s of extra transients):
  `scripts/testing/run_and_log.sh ANS-9-step1a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 480 mpiexec -n 2 python3 examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/09_coil_driven_gel_heating_10MHz.py'"`
  The script exits non-zero on any red anchor.
- **W2, the pin test**, after `git add …/metrics.json`:
  `scripts/testing/run_and_log.sh ANS-9-pin-staged "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src timeout -k 30 30 python3 -m pytest tests/unit/test_doc_reference_exit_codes.py -s -v --tb=short'"`
- **W3:** `scripts/testing/run_and_log.sh ANS-9-list "./scripts/run_examples.sh --list"`. Write its prediction first: the `ans:9 -> …/09_…py` line.
- **Post-censuses** `ANS-9-post-figures` / `ANS-9-post-docrefs`. Write the prediction first: runnable +1, docrefs `exit != 1`, setup figures `missing` +1 for `ans:9` (the `ANS-7` / `ANS-8` precedent, drawn later by the `EX-57` fallback), `broken = 0`.

**Traps already paid for.**
- **Name every `Function` before the XDMF write.** Unnamed θ / q_h Functions collided in HDF5 and hung `EX-79`'s first run to Status 124 (`20261002T113202Z_EX-79.log`).
- **No dots in the XDMF basename.** `write_xdmf_with_tags` treats the text after the last dot as a suffix (`02_…py:210–211`).
- **Snapshot indices are full steps only**, and index 1 lands after **two** BE half-steps (`transient.py:325`). θ^{½} needs its own run. Never reconstruct row 6 from the balance identity itself: that makes (R6) true by construction (rule (h)).
- **Every reduction takes `np.real`.** `assemble_scalar` and owned dofs are rank-local: reduce before printing, and keep collectives out of `rank == 0` branches.
- **`evaluate_vector_field_parallel` is collective**, so the point array must be bitwise identical on every rank (`OPS-40/55`).
- **`@complex_only` does not apply on import**, so the script carries its own guard. Import no `test_*` names.
- **Never open, import or create `aed/` or `aed_results/`.** No AED value goes in any file. `COMPARISON.md`'s AED columns stay blank.
- **The pin test reads `git ls-files`**, so stage `metrics.json` first.
- **After a kill:** run `pgrep -c python3` in the container, then the daily-review cache sweep (any `.c` without its `.so`) as one bare call. Status 124 is a finding, not a retry.
- **Rule (i):** no edit to the gate module after W0, and none to the script after W1. Never delete a log. `-s` everywhere; no `-k` and no pipes inside the quoted command.

**Scope.** Drive P1 only, degree-1 EM (`MAT-15`'s fixture), P1 thermal, one rung, constant σ, and the 1 W normalisation. **Not covered:**
- any AED value, verdict or private file;
- a degree-2 EM column (the order rule's default; the 1 W normalisation leaves the EM order only the source *shape*, so a degree-2 column is a later step's);
- (8c), meaning σ(T) or perfusion;
- an absolute temperature;
- tissue;
- `ANS-9` steps 2 / 3.

No `PREDICTED RED` line: every assertion is predicted green.

**Status it can move.**
- §7 `ANS-9` ⬜ → 🟡 "runnable half landed", and the SPEC's first Status box ticked.
- The row has no Done-when yet. Under rule (j) the queueing review writes the 🟡 / ✅ ladder into the row: runnable half ⇒ 🟡; the weekly's adjudication of the operator's run ⇒ ✅.

**Negative result.**
- **W0 red, or its digits off `EX-79`'s beyond noise:** the additive keys changed behaviour. Revert them, report, and mark the item BLOCKED.
- **(A0) red:** the SPEC and `MAT-15` disagree on a constant. That is a SPEC defect for the weekly, so stop before any row is written.
- **(R1), (R3) or (R0) red with W0 green:** a divergence between the script and the gate. Park on `attempt/ANS-9-step1a-<ts>`, write a known-issues entry, and mark the item BLOCKED (rule (d)).
- **(R6) red with (R0) green:** the surface sum mis-weights a (half-)step. Print the per-step residuals, park, and mark the item BLOCKED. Never widen `E_BAND`.
- **A control not wired:** a script bug. Fix it once, then re-run W1.
- **A printed record outside its prediction:** report it; it has no status effect.

**Outcome (2026-10-04, 09:00 CDT slot, take-next): DONE, every anchor green; `ANS-9` ⬜ → 🟡.**
- **W0 re-gate** `20261004T141433Z_ANS-9-step1a-mat15-regate.log`: 35 passed, 220 s. P1 matches `EX-79`'s record to 2.6e-14 relative, ⟨θ⟩ = 1.153314.
- **Final W1** `20261004T142653Z_ANS-9-step1a-W1.log`: exit 0, 175 s, `-n 2`. Every anchor passed:
  - (A0) all entries;
  - (R1) 7.3e-15, V_sub 3.447e-13;
  - (R0) 2.59e-15, balances ≤ 4.1e-12;
  - (R3) 5.7543e-6;
  - (R6) 7.09e-15;
  - (d6) 1.000000 and (s6) 0.1312, both above their floors.
- **Rows (all as predicted except row 4's ordering).**
  - Row 1: 1.000000000000007 W.
  - Row 2: ⟨θ⟩ 0.318096 / 0.614315 / 0.892104 / 1.153314 K.
  - Row 3: slope 1.058563e-3 K/s, against 1/(ρcV) = 1.058569e-3.
  - Row 4: A 0.905798, B 3.564854, C 0.685305, D 0.697444, E 0.958639 K.
  - Row 5: max 4.457206 K at r = 0.0275 m, azimuth +3.3°.
  - Row 6: 164.497813 J.
  - Row 7: 45 428.075518 W, 4.8e-12 from `ANS-2`'s `metrics.json` ratio.
  - The main run's first half-step deficit is 5.06e-3.
- **Prediction miss (printed only).** Row 4 orders B > D > C (D above C by 1.8 %), against the predicted B > C > D. This bears on the SPEC's ordering readout, which is the weekly's to read.
- **Two script fixes after the first W1, both disclosed.** The brief allows one. Neither moved a computed number, and the final W1 ran after the last edit (script 14:26:51Z, W1 start 14:26:53Z).
  - `20261004T142026Z_…W1.log`: `owner="lowest_original_cell"` raises `Mixed topology unsupported` on a `create_submesh` mesh. θ is continuous P1, so the slot uses the default collective owner, with points identical on every rank.
  - `…142345Z_…W1.log`: every anchor passed, but rank 1 crashed printing the `None` that `write_xdmf_with_tags` returns on ranks other than 0.
- **The other windows.**
  - Pin test `…143000Z_ANS-9-pin-staged.log`: 20 passed, pinned = git = 91.
  - `--list` `…143012Z_ANS-9-list.log` shows `ans:9`.
  - Censuses pre `…141422Z` / `…141424Z` and post `…143016Z` / `…143018Z`. Runnable went 82 → 83. Setup figures `missing` went 0 → 1 (ans:9, left for the `EX-57` fallback), `broken` 0. Docrefs exit 2 with dead 0 and guide 0; stale went 36 → 38 from older outputs ageing, not this item.

### §7 row text moved here verbatim — 2026-10-04 03:00 daily review (row near 2 KB)

> **Progression (steps written now, placed by the weekly):** **1 — one-way, gel:** HFSS drives an existing four-leg fixture (`ANS-2`'s 10 MHz loaded birdcage or `ANS-7`'s tuned coil at 64 MHz — the weekly picks), its SAR maps into the thermal solve as the volumetric source; transient conduction in the gel with a stated outer boundary (Robin, as `MAT-9`/`MAT-12`) for a stated duration. Compared: temperature rise at pre-registered probe points and times, the initial slope `dT/dt|₀ = SAR/c` (the calorimetry quantity, `MAT-12`'s anchor), and the volume-averaged rise; our side is `MAT-10` ✅ (FEM SAR as the heat source) + `MAT-12` ✅ (transient). **2 — two-way, gel:** temperature fed back into σ(T) and the EM re-solved to the coupled fixed point (AEDT's HFSS ↔ thermal iteration; ours `MAT-5` step 2's fixed point). **3 — perfusion (Pennes), tissue:** only if the Ansys thermal solver carries a perfusion term natively or by a stated temperature-dependent sink — **confirm on the Ansys side before the SPEC promises it**; ours `MAT-11` ✅. AED/Ansys numbers private as every ANS case (§5.4); relative agreement levels public (operator ruling 2026-09-24). **Independent of `ANS-8`** — it needs only the four-leg fixture and landed Phase 8 chunks, so it may be commissioned before the human-scale case. Solve resources on the Ansys side: 12 cores, 1 TB (operator 2026-09-29).

### Adjudication record and known-issues entry (2026-10-10 weekly; noted 2026-10-11)

The 2026-10-10 weekly adjudicated step 1 (§10 "Benchmarks, 2026-10-10 (step 5)"). It found AGREE on the EM → thermal hand-off, with rows 1 / 2 / 3 / 6 within 0.8 % / 0.5 % / 1.4 % (bands 5 %) and probes A–C within 2 % with B > D > C held. It found DISAGREE at probe D: ours is 11 % below, against a 10 % band. The numbers are in `docs/private/t13-t29-ans9-adjudication-2026-10-10.md`. `ANS-9` holds 🟡.

The known-issues entry the weekly filed is **"2026-10-10 — 🟡 OPEN (`ANS-9` step 1 adjudication, weekly review; §10 chain P21) — the coil-driven gel phantom's far-side probe D is 11 % colder in our thermal solve than in AED's (band 10 %), while the integral rows and the near-feed probes agree within 2 %"** (`docs/testing/known-issues.md`, the last open entry before "Retired entries"). Its reading names three candidates, in order:
1. the EM order (AED First Order against our degree 1);
2. the ragged pre-10-03 port sheets;
3. the termination model.

Diagnosis is step 1b (§10 P21), below.

### §9 item 371a — 2026-10-11 (chain P21 (a) — `ANS-9` step 1b, local SAR at the monitor points on two EM columns; chain step (solve))

**Outcome (2026-10-11, 13:30 CDT implementer slot): DONE, every anchor green, ρ = 2.119 ⇒ "no branch".** Module `tests/validation/test_ans9_far_side_sar.py` (imports only, no `src/` or existing-module edit). Cache sweep `20261011T184108Z_ANS-9-step1b-cache-sweep.log` (no orphans). W0 `20261011T184252Z_ANS-9-step1b-W0.log` (5 passed, 3 skipped, smoke). W1 `20261011T184301Z_ANS-9-step1b-W1-d1.log` (heavy, `-n 8`, 136 s, 7 passed). W2 `20261011T184526Z_ANS-9-step1b-W2-d2x.log` (heavy, `-n 8`, 180 s, 7 passed, `[capture] rc=0`).
- **(P1W):** d1 8.438e-15, d2x 1.110e-16, d1x 2.220e-16 (band `EXACT_IDENTITY_RTOL` 1e-10). Every point evaluation valid.
- **(F0):** d2x against `ANS4_EXACT_STRIP_S_P2_X1["10"]` 7.310e-11; d1x against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]` 6.597e-11 (band 1e-6). **(F0-ctl):** d2x against the degree-1 record misses by 1.020e-2 (predicted ≈ 2e-2; asserted > band). 117 206 cells.
- **(ORD), d1:** SAR_B 17.18 > SAR_D 2.651 > SAR_C 0.0129 W/kg. Held.
- **SAR at A / B / C / D / E, W/kg at 1 W in the phantom:**
  - d1: 0.4642 / 17.18 / 0.01289 / 2.651 / 2.552; D/C 205.56, D/A 5.709; G = 1.8171e7 1/W.
  - d1x: 0.5859 / 17.71 / 0.3028 / 3.094 / 2.283; D/C 10.217, D/A 5.281.
  - d2x: 0.3697 / 18.49 / 0.007034 / 3.064 / 2.585; D/C 435.63, D/A 8.288.
- **Decision (pre-registered, unchanged):** ρ = 435.63 / 205.56 = **2.119**, which lies in (1.5, 3), so **no branch** fires. The 10-14 weekly reads it.
- **Printed prediction missed:** (D/C)_d1 was predicted in [3, 12] and reads 205.6. C at (0, 0.020, 0) sits near a field null on the adjacent side: its SAR is 23× higher at d1x than at d1 and 43× lower at d2x than at d1x. D/C therefore mostly measures C's resolution, not the far-side source. **Recorded for the weekly, not acted on:** SAR_D moves +15.6 % d1 → d2x. The sheets and phantom mesh step (d1 → d1x) carries +16.7 % of that and the order step (d1x → d2x) −1.0 %. D/A moves 5.71 → 8.29.
- **Other prints:**
  - Blind-σ: the ratio is exactly σ_wire/σ_saline = 1600.
  - θ₄₀ at A–E: equal to `ROW4_POINT_K` at all six printed digits.
  - `ANS-2` r = 0.015 m named points: our N1curl owner-rule value against `ANS-2`'s CG1-interpolated `metrics.json` (both at 1 W): +x 10.02 vs 10.68, +y 0.146 vs 0.190, −x 1.117 vs 1.192, −y 0.182 vs 0.169 W/kg. These are route differences, printed only.
- **Not done:** 371b (the thermal half). F–I are still untracked.

*Drafted by a general-purpose agent for the 2026-10-11 03:00 daily review; the review reviews it.* Licence: §10 chain P21 ("Chain P — Phase 8, gel first", the 2026-10-10 weekly), the weekly's text executed as written. **No band, tolerance, Done-when, decision rule or target is set or moved.** Bands are imported by name. **Our values only: no slot opens `docs/private/`, `aed_results/`, `aed/` or `COMPARISON_private.md`.** P21 is split in two because one slot cannot hold both halves' code plus three heavy windows (≈ 63 predicted min). **371a is (a)**, which includes the pre-registered decision ratio. **371b is (b)**, the thermal re-read, and it depends on 371a.

**Precondition finding: the monitor points.** The weekly names "AED's nine monitor points A–I".
- **The tracked SPEC carries five:** row 4, A (0, 0, 0), B (0.020, 0, 0), C (0, 0.020, 0), D (−0.020, 0, 0), E (0, 0, 0.030) m (`examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md:98`).
- **F–I appear in no tracked file.** They are defined only in the operator's untracked PyAEDT script under the case's `aed/` directory, which is excluded by `.git/info/exclude:9`. The drafter confirmed only that the names F–I occur there and read no coordinate.
- By the letter, therefore, **the F–I half is "write the anchor"**: the weekly or the operator copies the four coordinates into SPEC.md (they are model inputs, not AED results, but whether they may go public is the operator's call).
- **Recommended:** queue 371a / 371b on **A–E**, the SPEC's own points. The pre-registered decision rule reads only D, C and A, so it is fully executable. F–I join the printed column as a one-window follow-up once they are tracked.
- **If the review prefers the letter, mark 371a BLOCKED** with the unblock condition "F–I coordinates in the tracked SPEC".

**`ANS-9` step 1b (a): local SAR ½σ|E|²/ρ at A–E, 1 W in the phantom, on two columns.**
- **Executor:** implementer.
- **Tier:** heavy, `-n 8` (the weekly's two windows), complex. **Class:** chain step (solve).
- **Cost:** M measured / ≈ 50 predicted min (15 fixed + 20 code + W0 ≈ 1 + W1 ≈ 3 + W2 ≈ 5 + spare). **If W2 has not started by minute 45, park** (rule (d)).
- **Depends on:** nothing open.
- **Gitignored inputs:** none. W1 writes `logs/ans9_1b/d1.json`, which W2 reads for the ratio. If that file is absent, W2 prints "W1 JSON absent" and the journal computes the ratio from the two logs.

**Columns.**
- **d1, `ans:9`'s own EM solve:** `_run_heated(return_fields=True)` (`tests/validation/test_mat15_coil_driven_heating_phantom.py:180`). This is `ANS-2`'s 10 MHz fixture with the pre-10-03 ragged sheets, 199 920 cells, degree 1. Drive k = 0 is P1 at +x (`20261004T142653Z_ANS-9-step1a-W1.log:1895`), with G = 1/P₁ as `ans:9`.
- **d2x, T6's convention:** item 192's exact-rectangle fixture at degree 2 and 10 MHz. Build it as `rung = _four_port_rung("ANS-9 1b exact 10 MHz", np.zeros(4), frequency_hz=FREQUENCY_HZ, degree=2, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`, κ off, 117 206 cells (`tests/validation/test_port_birdcage_leg_offset_sweep.py:199`; `test_ans4_exact_strip_records.py:435–448`). Then take P1's field from `_solve_driven(rung, "P1")` (`tests/validation/test_birdcage_b1_plus_map.py:142`). The rung dict carries `mesh`, `problem`, `facet_tags`, `specs` and `degree` (`…leg_offset_sweep.py:513–551`). G = 1 / `mean_sar(...)["dissipated_power_w"]` over `PHANTOM_CELL_TAG`.
- **d1x, a decomposition column on d2x's mesh, printed only:** `_four_port_rung(..., degree=1, reuse=rung)` followed by `_solve_driven(·, "P1")`. It costs ≈ 15 s at `-n 8`. It splits the d1 → d2x move into sheets + phantom mesh (d1 → d1x) and order (d1x → d2x). This column is not in the weekly's text and nothing is asserted on it.

**Code.** One new module, `tests/validation/test_ans9_far_side_sar.py`: `@complex_only`, opt-in knob `FEM_EM_ANS9_1B_COLUMN` ∈ {unset, `d1`, `d2x`}. Unset ⇒ the gated tests skip. Any other value raises `ValueError`; this is the `FEM_EM_ANS4_EXACT_DEGREE` pattern (`test_ans4_exact_strip_records.py:336–395`). There is no `src/` change and no edit to any existing module.
- **Imports, never restated.**
  - **Probes:** `_load_spec_probes` (`tests/validation/test_mat19_probe_coil_phantom.py:96`; it reads `ans:9`'s `SPEC_THERMAL["probes"]`, A–E), and `ROW4_POINT_K` / `ROW2_MEAN_K` (`:91–92`) for the d1 reproduction print.
  - **Fixtures:** `_run_heated`, `PHANTOM_CELL_TAG` (`tests/validation/test_port_birdcage_lumped_column.py:101`), `_four_port_rung`, `_solve_driven`, `GATED_WIDTH_FRACTION` (`tests/validation/test_port_lumped_narrowed_sheet.py:132`), `FREQUENCY_HZ` (as `test_ans4_exact_strip_records.py:106` imports it), and `SALINE_SIGMA` (`tests/validation/test_lossy_sphere_fullwave.py`).
  - **Other constants:** `PHANTOM_RHO_KG_PER_M3` (`test_birdcage_b1_plus_map.py`); `ANS4_EXACT_STRIP_S_P2_X1`, `ANS4_EXACT_STRIP_S_P1_RECORDS` (`test_ans4_exact_strip_records.py:141, :233`); `STEP3A_RECORD_RTOL` (`tests/validation/test_ans4_resolution_ladder.py:198`); `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`).
  - **Helpers:** `mean_sar` (`src/fem_em_solver/post/sar.py:57`); `evaluate_vector_field_parallel` (`src/fem_em_solver/post/evaluation.py:224`); `_allsum`, `_say` (`tests/validation/test_mat10_sar_to_heat_sphere.py:65–70`).
- **Point SAR, per column.** Evaluate the **N1curl field itself** at A–E on the parent mesh: `evaluate_vector_field_parallel(e, pts, owner="lowest_original_cell")`, which gives a width-invariant owning cell (`POST-8` step 2).
  - The point array must be identical on every rank. Assert every `valid`.
  - |E|² = Σ|Eᵢ|² and SAR = G·`SALINE_SIGMA`·|E|²/(2ρ), with σ uniform in the phantom.
  - Print SAR per point, D/C and D/A.
- **d1 also prints:**
  - `ANS-2`'s r = 0.015 m named points under P1, as a fixture sanity record: our N1curl value beside `ANS-2`'s CG1-interpolated `metrics.json` `pointwise["0"]` × G, labelled by route;
  - θ₄₀ at A–E from `out["T_star"][0]` through the default owner on the submesh (item 218's paid-for route), beside `ROW4_POINT_K`.
- **W1 writes** `logs/ans9_1b/d1.json` (rank 0) with the d1 SAR at A–E, D/C, D/A and G.
- **W2 reads it** and prints the decision ratio ρ = (D/C)_d2x / (D/C)_d1 and its branch.
- **The ungated `test_ans9_1b_coding`** (W0) asserts three things:
  - the knob parser;
  - that `_load_spec_probes()` holds exactly A–E, with D = (−0.020, 0, 0) and C = (0, 0.020, 0);
  - the branch function on synthetic ratios: 3.0 ⇒ "degree-1 ragged source", 1.0 and 1/1.5 ⇒ "far-leg current model", 2.0 and 0.5 ⇒ "no branch".

**Anchors (asserted; bands imported).**
- **(P1W) Phantom power = 1 W per column.** ∫_phantom G·½σ|E|² dx is assembled on its own form, `_allsum` over `dx(PHANTOM_CELL_TAG)` (item 361's "_vol" route), and must equal 1 W within `EXACT_IDENTITY_RTOL`. This is an independent assembly against G from `mean_sar`.
  - Predicted ≤ 1e-14. Item 218's (R1) read 7.3e-15 on d1, and item 361's (Σ_P) read 0.0.
- **(ORD), d1:** SAR_B > SAR_D > SAR_C. This is the SPEC's B / C / D readout ordering as `ans:9`'s rises have it (`ROW4_POINT_K`). Read "reproduces `ans:9`'s probe rises' ordering" as this: the far-side D sits between the driven-side B and the adjacent-side C.
  - The prediction comes from a nearby comparison, so it is *predicted*, not same-points. `ANS-2`'s P1 SAR at r = 0.015 m reads +x 5.878e-7 > −x 6.558e-8 > +y 1.045e-8 W/kg per drive (`ans2…/metrics.json` `pointwise["0"]`), so D/C ≈ 6.3 and B/D ≈ 9.0.
  - **Question for the review:** if it intended "the d1 θ₄₀ reproduces `ROW4_POINT_K`'s ordering" instead, that is the printed θ line, predicted equal to ≈ 1e-9.
- **(F0), d2x and d1x: the fixture is item 192's.**
  - d2x: the rung's 10 MHz S against `ANS4_EXACT_STRIP_S_P2_X1["10"]` within `STEP3A_RECORD_RTOL`, ‖ΔS‖/‖S‖.
  - d1x: against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]`, same band.
  - Backing for the same comparison on this family: item 192's (X) read 8.9e-13 at `-n 8` (`ANS-4.md`, item 192 outcome), and `ANS-4` step 4 W1 reproduced the 2a″ digits to 7.1e-11 across widths at `-n 8`.

**Printed (the weekly's readout; nothing asserted).**
- **Per column:** SAR at A–E (W/kg at 1 W in the phantom), D/C, D/A.
- **Predicted (D/C)_d1 in [3, 12]** (`ANS-2`'s 6.3 at r = 0.015 m, moved to r = 0.020 m).
- The d2x and d1x values are the measurement, so they have no prediction.

**Pre-registered decision rule (the weekly's, unchanged).** ρ = (D/C)_d2x / (D/C)_d1.
- **ρ ≥ 3** ⇒ the finding is the degree-1 ragged source. Our `ANS-9` column is re-recorded at degree 2 exact, and the 10-14 weekly re-adjudicates D.
- **"Within 1.5×", read as max(ρ, 1/ρ) ≤ 1.5** ⇒ the far-leg current model (the 50 Ω terminations, the `PORT` lane), and known-issues stays open.
- **Anything else (1.5 < ρ < 3, or ρ < 1/1.5):** neither branch fires. The slot prints "no branch" and the 10-14 weekly reads it.
- **The slot only reports the branch. It re-records nothing:** under the first branch the re-record at degree 2 exact is (b) plus a later step (see 371b's scope).
- **D/A is printed beside, not used.** The rule's "within 1.5×" reading is the drafter's, and the review confirms it.

**Negative controls (rule (e)).**
- **(F0-ctl), asserted against a floor; the factor is predicted.** The d2x S against the **degree-1** record `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]` must miss by > `STEP3A_RECORD_RTOL`.
  - Predicted ≈ 2 %: p2-vs-p1 exact at 10 MHz read 1.98 / 0.43 / 0.72 % by class (`ANS-4.md`, item 192, "Printed").
  - This shows (F0) can tell the orders apart.
- **(Blind-σ), printed only:** the d1 point SAR with σ taken as the conductor's `SIGMA_WIRE_S_PER_M` instead of `SALINE_SIGMA`. It shows the σ lookup is live, by the ratio of the two σ, and is printed.

**Windows.** Foreground; the harness call is the whole Bash command, bare, repo-relative; Bash `timeout: 660000`; `-s`; `tests/environment/test_complex_mode.py` first. **Two heavy windows never run at once.**
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Then run the any-orphan cache sweep as one bare call, ID `ANS-9-step1b-cache-sweep` (daily-review.md rubric item 4).
- **W0** (smoke, `-n 1`):

  `scripts/testing/run_and_log.sh ANS-9-step1b-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_far_side_sar.py -v -s --tb=short'"`
- **W1, d1** (heavy, `-n 8`). Predicted harness 110–190 s, against `20261004T142653Z_ANS-9-step1a-W1.log`: 175 s at `-n 2`, mesh 61.6 s serial, EM 137 s, total 161 s (`:1895, :1925`).
  - EM and thermal scale with ranks; the gmsh build does not.
  - This is the first `-n 8` run of `_run_heated`. `MAT-19` step 2 ran it at `-n 4`.

  `scripts/testing/run_and_log.sh ANS-9-step1b-W1-d1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/ans9_1b && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS9_1B_COLUMN=d1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 480 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_far_side_sar.py -v -s --tb=short'"`
- **W2, d2x + d1x** (heavy, `-n 8`). Predicted harness 200–330 s, against `20261003T094036Z_ANS-4-item192-W3-p2.log`: 429 s for four degree-2 rungs at `-n 8`, with mesh 23.9 s and rung 26.9 s, and the 10 MHz four-drive sweep 70.7 s (`:2031`).
  - This window adds one rung, one degree-2 `_solve_driven` (≈ 20–60 s, one fresh factorisation), the degree-1 reuse rung (≈ 15 s), and degree-2 form compiles if they are not cached (0–120 s).
  - RSS: ≈ 15–20 GiB summed. For comparison, the conductor-halved rung was 36 GiB.
  - The container ceiling is 590 s, so §5.1's durable-capture idiom is copied verbatim, trailing `; exit $rc` included:

  `scripts/testing/run_and_log.sh ANS-9-step1b-W2-d2x "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/ans9_1b && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS9_1B_COLUMN=d2x PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_far_side_sar.py -s -v --tb=short > /workspace/logs/ans9-step1b-w2-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans9-step1b-w2-raw.log; cat /workspace/logs/ans9-step1b-w2-raw.log; exit \$rc'"`

**Traps already paid for.**
- **Never interpolate the degree-2 field to CG1 for the point read** (`ANS-2`'s `_pointwise_sar` route, `02_…py:310–349`). On d2x's mesh the phantom cells are ≈ 15 mm across (the global 0.015 m size, with no phantom size field). A CG node on the phantom surface takes the air-side value, and E sits 10 mm from the top face. Evaluate the N1curl field directly on the parent mesh with `owner="lowest_original_cell"`. That rule raises "Mixed topology unsupported" only on a `create_submesh` mesh (item 218, `…142026Z_…W1.log`), and this is the parent.
- `evaluate_vector_field_parallel` is collective: give it identical points on every rank (`OPS-40/55`).
- **Every scalar is `np.real(_allsum(...))` on every rank.** `mean_sar` is already reduced. No collective inside `rank == 0`; only the JSON is rank-0.
- `_four_port_rung` always runs its four-drive sweep, so (F0) costs nothing extra. Do not pass `width_correction_kappa`: κ stays off (the 2026-10-02 sheets-and-κ ruling).
- **Degree 2 is set on the rung** (`degree=2`), and `_solve_driven` takes `rung["degree"]` (`test_birdcage_b1_plus_map.py:161–165`). Never re-solve at degree 1 by default.
- After a kill: `pgrep -c python3`, then the any-orphan sweep. Status 124 is a finding, not a retry.
- **Rule (i):** no edit to the module after W0 starts; any edit re-runs W0 → W2.
- Never import `test_*` names, and never open `aed/`, `aed_results/` or `docs/private/`.

**Scope.** Our local SAR at the SPEC's probe points on two EM columns, plus a printed decomposition, and the weekly's ratio. **Not covered:**
- the thermal re-read (371b);
- any AED value, verdict or comparison;
- F–I;
- a re-record of `ANS-9`'s column;
- 1 g or 10 g averaging;
- the termination model itself;
- 64 / 128 MHz.

**Negative result.**
- **(P1W) red:** a normalisation or assembly defect. Stop, park on `attempt/ANS-9-step1b-<ts>`, and mark BLOCKED.
- **(F0) red:** the d2x rung is not item 192's. Stop before any SAR is read and mark BLOCKED. A red here invalidates T6's column.
- **(ORD) red:** the d1 source does not order as the rises do. That is a finding, the far-side rise being diffusion-led. **Report it and stop:** no PREDICTED RED range is registered, so a red stops the slot. The 10-14 weekly reads it, and W2 waits for the review.
- Whatever ρ comes out is the measurement; no branch is forced.

**Status it can move.**
- The known-issues entry above gets a "Measured (step 1b (a))" bullet with ρ and the branch.
- The `ANS-9` row's state line gets "step 1b (a) landed: ρ = …, branch …". The row stays 🟡, because its Done-when (✅ = the weekly's adjudication of the operator's run) is not the target here. Rule (j) is not engaged, and no row amendment is needed.
- The 10-14 weekly re-adjudicates D on it.

**Spawn prompt (verbatim sentences).** "Execute §9 item 371a per docs/planning/chunks/ANS-9.md#§9-item-371a. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Every Bash call that runs a window passes `timeout: 660000`; never background a window; run pytest with `-s`; two heavy windows never run at once. Never open docs/private/, aed_results/ or aed/. Report a red and stop; never widen a band or force a branch."

**Same-commit duties.** Mark the §9 item DONE or BLOCKED. Add the known-issues bullet, the row state clause, an outcome paragraph here, and the journal entry.

### §9 item 371b — 2026-10-11 (chain P21 (b) — `ANS-9` step 1b, the thermal half on the degree-2 exact source; chain step (solve))

**Outcome (2026-10-11, 13:30 CDT implementer slot): DONE, every anchor green.** Same module, extended under `d2x` (imports only; `_e_prime` and `D_CONTROL_FLOOR` also imported from step 1's module). W0b `20261011T185136Z_ANS-9-step1b-W0b.log` (5 passed, 4 skipped, smoke). W2′ `20261011T185145Z_ANS-9-step1b-W2b-d2x-thermal.log` (heavy, `-n 8`, 172 s, 8 passed, `[capture] rc=0`; 371a's anchors re-run green with the same digits — (P1W) 2.2e-16 / 2.2e-16, (F0) 7.346e-11 / 6.597e-11, (F0-ctl) 1.020e-2, ρ = 2.119196; thermal block 4.0 s).
- **Tag-3 submesh:** 537 cells (predicted 500–1 500), V_sub 2.195295e-4 m³, ∮ds 2.040655e-2 m²; κ_s 1.338 (d2x) / 1.433 (d1x).
- **(R1′):** 4.441e-16 / 0.0 (band `B_BAND` 1e-10). **(E):** CN 4.95e-14 / 3.61e-14, BE 2.85e-12 / 3.92e-12, control 1.69e-14 (band `E_BAND` 1e-9). (E′) against 1 W, printed: 3.0e-14 both. **(R3′):** 4.8876e-6 / 5.1784e-6 (band `S2_BAND` 1e-4; predicted 5e-6 to 2e-5 — d2x sits just below the printed range). **(d6′):** 1.000000 > `D_CONTROL_FLOOR`.
- **Row 2 ⟨θ⟩ at n = 10 / 20 / 30 / 40, K:** d2x 0.327896 / 0.633639 / 0.920643 / 1.190592; d1x 0.327423 / 0.632311 / 0.918342 / 1.187285 (`ans:9` d1 on its 58 866-cell phantom: 0.318096 / 0.614315 / 0.892104 / 1.153314 — a different thermal mesh, printed only).
- **Row 4 θ₄₀ at A / B / C / D / E, K:** d2x 1.070257 / 3.166010 / 0.929158 / 0.723341 / 1.195038 (D/C 0.7785, D/A 0.6759); d1x 1.050658 / 3.189674 / 0.914637 / 0.620227 / 1.203214 (D/C 0.6781, D/A 0.5903).
- **Order effect d2x / d1x on the same thermal mesh:** A 1.0187, B 0.9926, C 1.0159, **D 1.1663**, E 0.9932; ⟨θ⟩₄₀ 1.0028. At degree 2 the far-side probe D warms by 16.6 % while the source integral and every other probe move ≤ 1.9 % — although 371a's *point* SAR_D moved −1.0 % with the order. The thermal D is read on a 537-cell phantom (≈ 15 mm cells, the P10 trap), so this is the cell-averaged source near D, not a resolved statement; the re-record on `ans:9`'s 0.0025 m phantom needs the fixture knob named in Scope.
- **Not done:** the `ANS-9` column re-record (Scope), F–I.

*Drafted for the 2026-10-11 03:00 daily review.* **Depends on item 371a landing, with (F0) and (P1W) green. If it did not land, skip.** It is the same module, extended. **No band, rule or target is set or moved.**

**(b) Re-run the thermal half on d2x's source and re-read rows 2–4.**
- **Executor:** implementer.
- **Tier:** heavy, `-n 8`, complex. **Class:** chain step (solve).
- **Cost:** M measured / ≈ 35 predicted min (15 fixed + 10 code + W0 ≈ 1 + W2′ ≈ 5–6 + spare). W2′ repeats 371a's degree-2 solve: rule (i) requires the module, as edited, to be re-run.

**Code.** In `tests/validation/test_ans9_far_side_sar.py`, under `FEM_EM_ANS9_1B_COLUMN=d2x`, add `ans:9`'s thermal half for **both d2x and d1x on d2x's mesh**.
- **Imports** from step 1's module: `_tag_submesh`, `_cell_average_dg0_tag`, `_transient`, `_be_one`, `T_INF`, `T_STAR`, `N_MAIN`, `DT_1`. From `MAT-10`: `_robin_tags`, `_to_sub`, `B_BAND`. From `MAT-12`: `RHO_C`, `E_BAND`, `S2_BAND`.
- **Per column:**
  1. Form `q_h = _to_sub(*_cell_average_dg0_tag(msh, cell_tags, PHANTOM_CELL_TAG, 0.5*sigma*inner(e, e)) × G)`.
  2. Run `_transient(sub, ft, q_h, dt=T_STAR/N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, snapshot_steps=range(N_MAIN+1))`, `ans:9` step 1a's schedule (item 218, computation step 2).
  3. Read the rows:
     - **row 2:** ⟨θ⟩ at steps 10 / 20 / 30 / 40;
     - **row 3:** `_be_one` slope against 1/(ρc·V_sub);
     - **row 4:** θ₄₀ at A–E through the default owner on the submesh.
- **Print:**
  - the tag-3 submesh cell count (predicted 500–1 500, against `MAT-15` step 3's 537 cells at the same default sizing on `ANS-7`'s ×1 exact-strip fixture);
  - rows 2–4 per column, beside `ans:9`'s degree-1 values (`metrics.json`, our own);
  - row 4's D/C and D/A, and the order effect (d2x / d1x) on the same thermal mesh.

**Anchors (asserted; bands imported):**
- **(R1′):** |∫_sub q_h − 1 W| ≤ `B_BAND` per column;
- **(E′):** every run's `max_balance_residual` ≤ `E_BAND`;
- **(R3′):** |ρc·θ̄¹/(Δt₁·∫q/V_sub) − 1| ≤ `S2_BAND` per column. Predicted 5e-6 to 2e-5 (κ_s-dependent; item 218 read 5.75e-6 on d1's fine phantom).

**Negative control:** (d6′), doubled source, with (E′) against 1 W missing by the predicted 1.000. It is asserted only > `D_CONTROL_FLOOR` (`ANS-9` 1a's form).

**Window.** W0 as in 371a, ID `ANS-9-step1b-W0b`. Then **W2′**: 371a's W2 command with ID `ANS-9-step1b-W2b-d2x-thermal` and raw file `ans9-step1b-w2b-raw.log`, otherwise verbatim (durable capture, `; exit $rc`). Predicted 210–345 s (W2 + ≈ 10 s of thermal).

**Traps.**
- All of 371a's traps apply.
- **A 537-cell-class phantom reads the mesh at a point.** P10's red (2026-10-07 weekly) was a 5 mm probe on the 537-cell tuned phantom. Row 4 on this mesh is therefore printed, and the order effect is read only as d2x / d1x on the *same* thermal mesh, never d2x against `ans:9`'s 58 866-cell column.
- `_transient` snapshot 1 is after two half-steps.
- Never reconstruct row 6 from the balance (rule (h)). Row 6 is not asked for.

**Scope.** These are readings on item 192's coarse-phantom mesh. **This step does not re-record `ANS-9`'s column**, even if 371a's ρ ≥ 3 fires the first branch. A re-record fit for the 10-14 weekly's re-adjudication needs exact sheets *and* `ans:9`'s 0.0025 m phantom. `_four_port_rung` / `_build` take `leg_sheet_strip_fraction` but no `phantom_resolution`, and `build_four_port_sweep` takes the reverse. That is an additive fixture knob, a later step the weekly enumerates. It is flagged here and not built. **Not covered:** AED, F–I, rows 5–7, σ(T), perfusion.

**Negative result.**
- **(R1′) or (E′) red:** a source-transfer defect on the new mesh. Park and mark BLOCKED.
- **(R3′) red with (E′) green:** a time-truncation reading on a coarse phantom (κ_s large). Report κ_s and stop; `S2_BAND` is never widened.
- No `PREDICTED RED` line.

**Status it can move.** The same as 371a: a known-issues bullet and the row state line, with no glyph. Rows 2–4 at degree 2 exact are printed for the 10-14 weekly's private re-read.

**Spawn prompt (verbatim sentences).** "Execute §9 item 371b per docs/planning/chunks/ANS-9.md#§9-item-371b; skip it unless item 371a is DONE with (F0) and (P1W) green. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Every Bash call that runs a window passes `timeout: 660000`; never background a window; run pytest with `-s`. Never open docs/private/, aed_results/ or aed/."

**Same-commit duties.** Mark the §9 item DONE or BLOCKED. Add the known-issues bullet, the row state clause, an outcome paragraph here, and the journal entry.

**Review ruling 2026-10-11 (03:00 daily) on items 371a / 371b.** (1) Only monitor points A–E are in the tracked SPEC (`SPEC.md` row 4); F–I are not read from any untracked file. The decision rule reads D, C and A only, so 371a / 371b run in full on A–E; F–I are a one-window follow-up once their coordinates are in the SPEC (named for the 10-12 weekly). (2) "Within 1.5×" is read as a D/C ratio in [1/1.5, 1.5]; a ratio in (1.5, 3) or below 1/1.5 fires no branch and goes to the 10-14 weekly. (3) The asserted ordering is the degree-1 local SAR ordering B > D > C. (4) 371b does not re-record our `ANS-9` column (item 192's coarse phantom); that re-record needs a fixture knob the weekly enumerates. Readings (2) and (3) are provisional, for the weekly to ratify.

### §9 item 396 — 2026-10-12 (chain P23, daily provisional under §9 rule (6); chain step (solve))

*Drafted by a general-purpose agent for the 2026-10-12 03:00 daily review; the review reviews it.*

**Licence.** §10 chain P, step P23, written by the 2026-10-12 03:00 daily under §9 rule (6) as a provisional next step on a live chain. The weekly ratifies or strikes it. It re-runs P21 (a) + (b)'s d2x window (371a / 371b, above) with the phantom refined. **No band, tolerance, Done-when, decision rule or target is set or moved.** Every band is imported by name. **Our values only: no slot opens `docs/private/`, `aed_results/`, `aed/` or `COMPARISON_private.md`.**

**Why this step.** P21 left the weekly two readings of the far-side finding, and P21's own data cannot separate them:
1. 371b's θ₄₀ at D moved **+16.6 %** with the order (d1x → d2x) on the same thermal mesh, while 371a's *point* SAR_D moved **−1.0 %** (`20261011T185145Z_ANS-9-step1b-W2b-d2x-thermal.log:2053`, `:1981`, `:2038`). The 2026-10-11 18:00 review read this as the coarse phantom's **cell-averaged** source near D. Item 192's mesh has 537 tag-3 cells, ≈ 15 mm across (`:2041`). On that reading, the order sensitivity is an artefact of the DG0 source on 15 mm cells (the hypothesis line of 371b's journal entry, `docs/testing/attempts.md:1554`).
2. The other reading: the degree-2 source really is different near D, and a refined phantom keeps the θ_D order effect.

A refined phantom tells them apart. Under (1), θ_D's order effect falls toward the point-SAR order effect as h_p falls. Under (2), it stays. **That is a reading for the weekly, not a decision rule this item applies.**

**Precondition finding: the knob.** GEO-27 step 2's knob does **not** reach the ANS-9 fixture builder as it stands.
- 371a / 371b build d2x through `_four_port_rung` (`tests/validation/test_port_birdcage_leg_offset_sweep.py:199`). That calls `tests/mesh/test_birdcage_leg_offset._build` (`:107–171`), which takes `leg_sheet_strip_fraction` but **no `phantom_resolution`**. 371b's Scope already flagged this.
- GEO-27 step 2 (item 377) builds through `tests/mesh/test_birdcage_port_sheets._build(True, phantom_resolution=h_p, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` (`:80–177`). It then calls `build_four_port_sweep(built=…, build_only=True, degree=2)` (`tests/validation/test_port_birdcage_four_port.py:219`) and `run_n_port_sparameter_sweep(…, keep_fields=True, degree=2)` (`tests/validation/test_geo27_one_gram_rung_degree2.py:150–192`). Call this **route B**.

**Smallest change: none to any existing module.** P23 uses route B as GEO-27 runs it and ties route B to item 192's fixture in its own window, `x1`.
- The two builders pass the identical `MeshGenerator.birdcage_port_domain` arguments. The geometry constants are the same imports from `tests/mesh/test_birdcage_port_tags`, and `CONDUCTOR_RESOLUTION` comes from the same module.
- There are two differences:
  - the leg-offset builder passes `leg_azimuth_offsets_rad=zeros`, which the generator documents as bit-for-bit the default `None` (`src/fem_em_solver/io/mesh.py:4174–4177`);
  - the narrowing helper (`_narrowed_radial` vs `_narrowed_transverse`). Both select the exact centre rectangle on exact strips, and for legs at 0 / 90 / 180 / 270° the radial and transverse axes coincide.
- The materials, the 50 Ω terminations, the PEC boundary and κ off are the same lines in both builders (`test_port_birdcage_four_port.py:469–513` against `test_port_birdcage_leg_offset_sweep.py:357–392`).
- **Rejected alternative:** an additive `phantom_resolution=None` keyword on the leg-offset `_build` and `_four_port_rung`. That is two lines each, but it edits a fixture that a large set of `ANS-4` / `PORT-9` gate modules import, and rule (a)/(c)'s re-run duty would cost more than the item. The review may prefer it; if so, this brief's anchors stand unchanged.

**`ANS-9` step 1c: point SAR at A–E, D/C, D/A and θ₄₀ (rows 2–4) on the exact-rectangle degree-2 fixture with the phantom refined.**
- **Executor:** implementer. **Tier:** heavy, `-n 8`, complex. **Class:** chain step (solve).
- **Cost:** ≈ 13 measured (class median) / ≈ 50 predicted min. That is 15 fixed + 15 code + W0 ≈ 1 + W1 ≈ 3 + W2 ≈ 3.5 + W3 ≈ 7 + spare.
- **If W3 has not started by minute 45:** commit the module and the W1 / W2 logs. Mark item 396 "⬜ W3 remaining (module at `<sha>`)": the item stays runnable, so rule (d) is met without BLOCKED. The next slot runs W3 alone, as committed (rule (i)).
- **Depends on:** nothing. 371a / 371b and 377 have landed.
- **Gitignored inputs:** none. The coarse and d1 comparands are module constants with log provenance (the GEO-27 `MAT4_STEP5B_*` precedent).

**Code.** One new module, `tests/validation/test_ans9_far_side_sar_refined.py`.
- `@complex_only`. Knob `FEM_EM_ANS9_1C_RUNG` ∈ {unset, `x1`, `h005`, `h0025`}. Unset ⇒ the gated tests skip; any other value ⇒ `ValueError` (the `FEM_EM_GEO27_STEP2` pattern).
- **Rung map:** `{"x1": None, **RUNGS}`, with `RUNGS` imported from `tests/validation/test_geo27_one_gram_rung_degree2.py:86`. That gives h005 → `PHANTOM_RESOLUTION_005` (`:85`) and h0025 → `PHANTOM_RESOLUTION_1G_RUNG` (`tests/validation/test_birdcage_sar_1g_rung.py:114`).
- **No `src/` change and no edit to any existing module.**

**Imports (names only, never `test_*`).**
- **From 371a / 371b's module** (`tests/validation/test_ans9_far_side_sar.py`, unedited): `_column_readout` (`:125`), `_thermal_readout` (`:148`), `_rel` (`:321`), `PROBE_NAMES` (`:90`).
  - Its import-time `needs_column` reads `FEM_EM_ANS9_1B_COLUMN`, which is unset in these windows, so it is inert.
- **Route B:** `_build`, `SHEET_IFACE` (`tests/mesh/test_birdcage_port_sheets.py`); `build_four_port_sweep`; `run_n_port_sparameter_sweep`; `GATED_WIDTH_FRACTION` (`test_port_lumped_narrowed_sheet.py:132`); `PHANTOM_CELL_TAG`.
- **Thermal:** `_tag_submesh` (`test_mat15_coil_driven_heating_phantom.py:106`), `_robin_tags`; `_global_tag_cell_count` (`test_birdcage_sar_mass_averaged.py`).
- **Bands and records:**
  - `EXACT_IDENTITY_RTOL` (`test_birdcage_sar_mass_averaged.py:163`, 1e-10);
  - `STEP3A_RECORD_RTOL` (`test_ans4_resolution_ladder.py:198`, 1e-6);
  - `ANS4_EXACT_STRIP_S_P2_X1`, `ANS4_EXACT_STRIP_S_P1_RECORDS` (`test_ans4_exact_strip_records.py:233, :141`);
  - `ANS4_EXACT_STRIP_CELL_RECORD`, `STRIP_BAND` (`test_ops61_exact_strip_ports.py:64, :58`);
  - `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`);
  - `I14_OFF_CELLS`, `I14_OFF_DOFS` (`test_imp3_tip_ladder_1g.py:94`);
  - `NOMINAL_STRIP_AREA` (`test_imp3_tip_ratio_1g.py:105`); `_sheet_areas` (`test_imp3_palace_twin.py:69`);
  - `RECIPROCITY_BAND` (`test_port_lumped_sheet_sweep.py:86`);
  - `B_BAND` (`test_mat10_sar_to_heat_sphere.py:61`); `E_BAND`, `S2_BAND` (`test_mat12_transient_sar_heating_sphere.py:76, :75`); `D_CONTROL_FLOOR` (`test_mat15_coil_driven_heating_phantom.py:102`);
  - `ROW4_POINT_K`, `ROW2_MEAN_K` (`test_mat19_probe_coil_phantom.py:91–92`).
- **Printed comparands defined in the module, never bands, each with its log line:**
  - `P21_D2X_SAR`, `P21_D1X_SAR` at A–E: `20261011T185145Z_ANS-9-step1b-W2b-d2x-thermal.log:1978–1982, :2035–2039`.
  - `P21_D2X_THETA40`, `P21_D1X_THETA40`: `:2046, :2052`.
  - `P21_ROW2_D2X`, `P21_ROW2_D1X`: `:2043, :2049`.
  - `ANS9_D1_SAR` at A–E: `20261011T184301Z_ANS-9-step1b-W1-d1.log:2013–2017`, `ans:9`'s own degree-1 ragged field on its 0.0025 m phantom.
  - `GEO27_H005_CELLS = (130427, 8498)`: `20261012T034731Z_GEO-27-step2-W1-h005.log:1942`.

**Body (all rungs, one module-scoped fixture, one window per knob value).**
1. `built = _build(True, phantom_resolution=h_p, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`. For `x1`, h_p is `None`.
2. `case = build_four_port_sweep(built=built, build_only=True, degree=2)`. Print the cells, the tag-3 count from `_global_tag_cell_count`, E, F and the (G2) strip areas.
3. `r2 = run_n_port_sparameter_sweep(case["problem"], case["port_defs"], lumped_sheet_ports=case["specs"], lumped_sheet_facet_tags=case["facet_tags"], keep_fields=True, degree=2)`, timed.
   - The S matrix is `np.asarray(r2.s_matrix)` in `r2.port_ids` order (P1…P4).
   - P1's field is `r2.fields["P1"]`: one factorisation, the other three ports terminated by their 50 Ω sheets. This is the GEO-27 / I14 route.
4. **d2x readout:** `_column_readout("<rung>-d2x", f.e_complex, f.sigma_field, case["cell_tags"], comm)`. It gives G, (P1W), SAR at A–E, D/C and D/A.
5. **Thermal half:** `sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)`, `ft = _robin_tags(sub)`, then `_thermal_readout("<rung>-d2x", …, control=True)` (371b's call, verbatim). Keep only the returned numbers, then `del r2` and `gc.collect()`.
6. **d1x on the same mesh:** `case1 = build_four_port_sweep(reuse={k: case[k] for k in ("mesh", "cell_tags", "facet_tags", "sheets", "halves", "cells")}, build_only=True, degree=1)`. This is the in-function κ route's own reuse dict (`test_port_birdcage_four_port.py:370`).
   - Run the sweep at `degree=1`, `keep_fields=True`.
   - Then `_column_readout("<rung>-d1x", …)` and `_thermal_readout("<rung>-d1x", …)` on the **same** `sub`.
7. **Prints:**
   - The PRICE line: rung, h_p, cells, tag-3, E, F, the d2 dofs, build s, d2 sweep s, d1 sweep s, thermal s, body s, `_peak_rss_gib()` (`test_imp3_rod_birdcage.py:498`, imported) and ranks.
   - Per column: SAR_A…E, D/C, D/A, and θ₄₀ at A–E with its D/C and D/A.
   - **The order effect at this rung:** point SAR d2x / d1x and θ₄₀ d2x / d1x at A–E, and ⟨θ⟩₄₀ d2x / d1x.
   - **The h_p move:** each d2x and d1x value ÷ its `P21_*` comparand − 1.
   - **At `h0025` only, the three-way split at matched nominal h_p.** `ans:9`'s d1 (ragged, 199 920 cells) → d1x (exact, degree 1) is the sheet step; d1x → d2x is the order step. It covers point SAR at A–E against `ANS9_D1_SAR`, and θ₄₀ at A–E against `ROW4_POINT_K`.
   - The weekly's ratio re-read, **printed only**: ρ(h) = (D/C)_d2x(h) / 205.564455, beside 2.119196. `_branch` is **not** called; P21's branch is the weekly's (2026-10-11 03:00 ruling (2)).
8. **`test_ans9_1c_coding` (W0, ungated):** the knob parser (`x1` / `h005` / `h0025` accepted, `h001` / `d2x` / `X1` raise); the rung map gives `None`, 0.005 and 0.0025; `PROBE_NAMES == ("A", "B", "C", "D", "E")`.

**Anchors (asserted; every band imported).**
- **Every rung:**
  - **(P1W)** per column (d2x, d1x): ∫_phantom G·½σ\|E\|² = 1 W within `EXACT_IDENTITY_RTOL`, with every point evaluation `valid` (371a's anchor, by name). Predicted ≤ 1e-14: 371a/b read 2.2e-16 at x1 (`…W2b…log:1977, :2034`).
  - **(G2):** the four N_i strip areas equal `NOMINAL_STRIP_AREA` within `STRIP_BAND`. GEO-27 W1 / W2 read 5.600000000000e-05 ×4.
  - **(i):** ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND`. GEO-27 W2 read 1.265e-14.
  - **(R1′)** ≤ `B_BAND`, **(E)** every run's `max_balance_residual` ≤ `E_BAND`, and **(R3′)** ≤ `S2_BAND`, per column (371b's anchors, by name). (R3′) is predicted at 2e-6 to 6e-6: κ_s falls with refinement; 371b read 4.9e-6 / 5.2e-6 at κ_s 1.34 / 1.43, and `ans:9` read 5.75e-6 on its 58 866-cell phantom.
- **`x1` only — the route tie.**
  - **(F0′):** d2x S against `ANS4_EXACT_STRIP_S_P2_X1["10"]`, and d1x S against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]`, each ≤ `STEP3A_RECORD_RTOL` (371a's (F0), by name).
  - **(X1c):** cells within `CELL_COUNT_BAND` of `ANS4_EXACT_STRIP_CELL_RECORD`; predicted 117 206 exactly.
  - Predicted for (F0′): ≤ 1e-10, the record-digit floor; 371a / 371b read 7.3e-11 / 6.6e-11 through `_four_port_rung`.
  - **Printed, not asserted:** x1's SAR at A–E and θ₄₀ at A–E against the `P21_*` comparands. Predicted equal to ≤ 1e-8. The held-factor field against 371a's fresh `_solve_driven` is the PORT-19 premise, unmeasured on this fixture. On a green (F0′), these prints are what lets the h005 / h0025 windows read their h_p move against the `P21_*` constants.
- **`h0025` only — (R):** cells within `CELL_COUNT_BAND` of `I14_OFF_CELLS`, and d2 dofs == `I14_OFF_DOFS` (GEO-27 W2's anchor, by name). This is where the knob is shown to have reached the solve.
- **`h005`:** cells and tag-3 are **printed** beside `GEO27_H005_CELLS`. There is no registered record at this rung.

**Negative controls (rule (e)).**
- **(F0-ctl), x1, asserted against a floor:** d2x S against the **degree-1** record misses by > `STEP3A_RECORD_RTOL`. The same comparison on the same fixture read 1.020e-2 (`…W2b…log:2080`), ≈ 1e4× the floor. The factor is printed, and predicted 1.02e-2.
- **(d6′), every rung, asserted > `D_CONTROL_FLOOR`:** the doubled source's (E′) against 1 W.
  - The same comparison read 1.000000 on item 192's mesh (`…W2b…log:2047`) and 1.000 on `ans:9`'s 0.0025 m ragged phantom (`ANS-9` 1a's d6, row state). It is linearity, so the value is mesh-independent; predicted 1.000.
  - It shows (E′) and the thermal source are live on the refined submesh.
- **(Live-h), printed only:** ‖S(h_p) − S_x1‖/‖S_x1‖ against `ANS4_EXACT_STRIP_S_P2_X1["10"]` at h005 / h0025. It shows the phantom knob moved the solve. Predicted nonzero and above 1e-6, with no size predicted. If it reads ≤ 1e-9, report: the knob did not reach the EM solve.

**Predictions (written first in the module; printed, never asserted).**
- **Cells / tag-3:**
  - x1: 117 206 / 537;
  - h005: 130 427 / 8 498 (GEO-27 W1 `:1942`, same builder call);
  - h0025: 201 072 / 58 960 (GEO-27 W2 `:1953`).
- **κ_s** falls from 1.34 (x1) toward 1.0–1.15 at h0025.
- **⟨θ⟩₄₀ order effect** stays within ±0.5 % of 1.0028 at every rung: the integral source agrees to < 1 % across order.
- **The θ_D order effect:** no prediction is set. Under reading (1) it moves from 1.166 toward the point-SAR_D order effect (0.990 at x1) as h_p falls; under reading (2) it stays near 1.17. Neither is a band. The weekly reads where it lands.
- **SAR_C stays a near-null at every rung,** so D/C and ρ(h) swing with C's resolution. 371a's own finding is that D/C mostly measures C.

**Windows.** Foreground. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Bash `timeout: 660000`; `-s`; `tests/environment/test_complex_mode.py` first. §5.1 durable capture is copied verbatim, trailing `; exit $rc` included. **Two heavy windows never run at once.**
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Then run the any-orphan cache sweep as one bare call, label `ANS-9-step1c-cache-sweep` (daily-review.md rubric item 4).
- **W0** (smoke, `-n 1`):

  `scripts/testing/run_and_log.sh ANS-9-step1c-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_far_side_sar_refined.py -v -s --tb=short'"`
- **W1, `x1`, the route tie** (heavy, `-n 8`). Predicted harness 120–170 s, against 371b's W2′: 172 s, which included one extra fresh degree-2 factorisation that route B does not pay (`…W2b…log:1955`: mesh 23.8 s, sweep 67.6 s, thermal 4.0 s). RSS ≈ 15–20 GiB summed.

  `scripts/testing/run_and_log.sh ANS-9-step1c-W1-x1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS9_1C_RUNG=x1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_far_side_sar_refined.py -s -v --tb=short > /workspace/logs/ans9-step1c-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans9-step1c-w1-raw.log; cat /workspace/logs/ans9-step1c-w1-raw.log; exit \$rc'"`
- **If W1 is red, stop here** (Negative result). Run pgrep again; it must read 0.
- **W2, `h005`** (heavy, `-n 8`). Predicted harness 130–200 s, against GEO-27 W1's 136 s (build 30.49 s, sweep 77.02 s, balls 22.17 s not paid here, 20.72 GiB; `…034731Z…W1-h005.log:1984`). This window adds a degree-1 reuse sweep (≈ 5–15 s) and three CN/BE runs on 8 498 cells (≈ 2–6 s). The command is W1's with label `ANS-9-step1c-W2-h005`, knob `FEM_EM_ANS9_1C_RUNG=h005` and raw file `/workspace/logs/ans9-step1c-w2-raw.log` (all three occurrences).
- **Then pgrep (must read 0).**
- **W3, `h0025`** (heavy, `-n 8`). Predicted harness 310–420 s, against GEO-27 W2's 323 s (build 41.67 s, sweep 246.74 s, balls 27.04 s not paid here, 44.05 GiB; `…034956Z…W2-h0025.log:1996`). It adds a degree-1 sweep on E = 238 627 dofs (≈ 15–40 s) and CN/BE on 58 960 cells (≈ 2.5 s per run at `-n 2` in `ans:9` step 1a, `20261004T142653Z_ANS-9-step1a-W1.log:1900`). Predicted RSS 44–50 GiB summed, under the 128 GiB container limit.
  - **Do not start W3 if W2's body exceeded 300 s.** That would be ≥ 1.5× its prediction, and W3's ≈ 2.4× scaling would then pass 590 s. Record "h0025: not attempted (W2 body … s)", which is not a red.
  - The command is W1's with label `ANS-9-step1c-W3-h0025`, knob `h0025` and raw file `/workspace/logs/ans9-step1c-w3-raw.log`.

**Traps already paid for.**
- **`_four_port_rung` has no `phantom_resolution`:** passing one is a `TypeError`. This is why the item uses route B. `build_four_port_sweep(built=…)` **ignores** its sizing keywords (`:241`), so h_p goes to `_build` only (GEO-27's trap).
- **Never `_solve_driven` per port, and not even once for P1 at h0025.** One extra degree-2 factorisation there costs ≈ 200–250 s and breaks the 590 s ceiling. P1's field comes from the one-factorisation sweep with `keep_fields=True` (GEO-27 / I14).
- **Degree 2 goes on both calls:** `build_four_port_sweep(..., degree=2)` *and* `run_n_port_sparameter_sweep(..., degree=2)`. d1x passes `degree=1` to both. Never pass `width_correction_kappa`: κ stays off (2026-10-02 ruling).
- **Free the degree-2 result before the d1x sweep** (`del`, `gc.collect()`). Keep only the readout numbers; the h0025 peak is 44 GiB without it.
- **Never interpolate the N1curl field to CG1 for the point read** (371a's trap). Evaluate on the parent mesh with `owner="lowest_original_cell"`. The θ read on the `create_submesh` mesh uses the **default** owner, because `lowest_original_cell` raises "Mixed topology unsupported" there (item 218). `_column_readout` / `_thermal_readout` already do both: import them, never copy-and-edit.
- `evaluate_vector_field_parallel` is collective, so give it identical points on every rank. Every scalar is `np.real(_allsum(...))` on every rank, and no collective goes inside `rank == 0`. Cell counts come from `size_global` / `_global_tag_cell_count`, never `cell_tags.values` lengths.
- **A 537-cell phantom reads the mesh at a point** (P10's red, 2026-10-07 weekly; 371b). This item is the measurement on finer phantoms. Its x1 window is the control, not a result.
- Importing from `test_ans9_far_side_sar.py` imports no `test_*` name and no fixture: pytest would collect them.
- **Orphan `.c` of any size** (item 160). After a kill or status 124, run `pgrep -c python3`, then the cache sweep, then park. Never raise the timeout or `-n`. A killed W3 is the price record "h0025 does not fit heavy here", not a retry.
- `-s` always. No `-k` inside the quoted command, no pipe. Copy the durable-capture idiom verbatim: a dropped `; exit $rc` reads 0 over a red (`WF-6` 4h).
- **Rule (i):** no module edit after W0 starts. Any edit re-runs W0 → W3.
- Never open `aed/`, `aed_results/`, `docs/private/` or `COMPARISON_private.md`.

**Scope.** This item measures our point SAR at A–E and `ans:9`'s thermal rows 2–4 on the exact-rectangle degree-2 (and degree-1) F-small fixture, 10 MHz, κ off, at three phantom resolutions (default / 0.005 / 0.0025 m). It adds a route tie to item 192 and the h_p and order moves. **It separates P21's two readings for the weekly; it does not decide P21's branch, re-evaluate ρ against the rule, or re-record `ANS-9`'s column.** The h0025 d2x column is the candidate re-record (exact sheets on `ans:9`'s own 0.0025 m phantom). Whether it becomes `ANS-9`'s column, and the re-adjudication of D, are the weekly's. **Not covered:**
- any AED value, verdict or comparison;
- F–I;
- 1 g / 10 g averaging;
- the termination model;
- 64 / 128 MHz;
- an h-convergence claim (three rungs, no rate asserted);
- rows 5–7;
- σ(T).

**Negative result.** Report the measurement and stop. Never widen a band, and never force a reading.
- **(F0′) or (X1c) red at x1:** route B is not item 192's fixture, so the refined windows' h_p moves cannot be read against P21. Stop before W2. Park on `attempt/ANS-9-step1c-<ts>` and mark item 396 BLOCKED, with the unblock condition "route-B / `_four_port_rung` difference located". Add a known-issues bullet.
- **(P1W), (R1′) or (E) red:** a normalisation or source-transfer defect on the refined mesh. Park and mark BLOCKED.
- **(G2), (i) or (R) red:** the refined fixture is not GEO-27's. Stop and mark BLOCKED.
- **(R3′) red with (E) green:** a time-truncation reading. Report κ_s and stop; `S2_BAND` is never widened.
- **(d6′) or (F0-ctl) not missing:** a control is blind. Stop.
- **Any θ_D order effect, D/C or h_p move** is the measurement, whatever it reads. No reading is forced.
- **W3 not attempted under its start rule, or killed at 590 s:** a price record, not a red. h005 still stands as the refined reading.
- No `PREDICTED RED:` line.

**Status it can move.** The status changes are records only; no glyph changes anywhere.
- **Known-issues:** the 2026-10-10 entry ("`ANS-9` step 1 adjudication … far-side probe D …", `docs/testing/known-issues.md:824`) gains a **"Measured (step 1c, §9 item 396)"** bullet. It carries, per rung, the θ_D and point-SAR_D order effects and the h_p moves of SAR_D, θ_D, D/C and D/A. The entry stays open.
- **`ANS-9` row:** the state line gains "step 1c landed: θ_D order effect x1 / h005 / h0025 = … / … / …; point SAR_D order effect … ; `<log>`". The row stays 🟡: its Done-when (✅ = the weekly's adjudication of the operator's run) is not this item's target. Rule (j) is not engaged, and no row amendment is needed.
- **§10 P23** gets a DONE line. The 10-14 weekly reads the two readings off it.

**Spawn prompt (verbatim sentences).** "Execute §9 item 396 per docs/planning/chunks/ANS-9.md#§9-item-396. Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Every Bash call that runs a window passes `timeout: 660000`; never background a window; run pytest with `-s`; two heavy windows never run at once; `pgrep -c python3` in the container reads 0 before each heavy window. No edit to any existing module or `src/`; no edit to the new module after W0 starts. Never open docs/private/, aed_results/, aed/ or COMPARISON_private.md. You are the executor: do not spawn agents. Report a red and stop; never widen a band, call `_branch`, or decide P21's reading."

**Same-commit duties.** Mark §9 item 396 DONE, BLOCKED, or "⬜ W3 remaining". Commit the module and the harness logs. Add the known-issues bullet, the row state clause, the §10 P23 DONE line, an outcome paragraph here with the per-rung table (rung, cells, tag-3, κ_s, SAR_D d2x / d1x, θ_D d2x / d1x, both order effects, D/C, D/A, h_p moves, body s, GiB), and the journal entry.

**Review ruling 2026-10-12 03:00 (daily) on the item 396 draft.**
(Q1) Route B (`GEO-27` step 2's builder) behind the x1 tie window is accepted; no keyword is added to `_four_port_rung`. (Q2) The x1 point-SAR / θ₄₀ reproduction stays printed. (Q3) Three windows accepted. If the h0025 window cannot start before minute 45, the slot commits the x1 and h005 windows with their logs as part (a), marks the §9 entry "⬜ — part (a) DONE, h0025 window remaining" — runnable, so rule (d)'s BLOCKED mark is not owed — and the next slot takes the h0025 window at its start. (Q4) ρ(h) is printed beside 2.119 without calling `_branch` (no branch is read in a slot; P21's reading is the weekly's). (Q5) The doubled-source linearity control (d6′) is **printed, not asserted**, on both refined rungs: rule (e) asks for a prior measurement of the same comparison on the same fixture, and neither refined mesh has one. (Q6) The h0025 d2x column is printed as a candidate only; adopting it as `ANS-9`'s column is the weekly's.

**Outcome 2026-10-12 (09:30 slot, implementer): DONE — all three rungs green.** Module `tests/validation/test_ans9_far_side_sar_refined.py`. W0 (`-n 1`, 3 s) 5 passed / 4 skipped; W1 x1, W2 h005, W3 h0025 (`-n 8`, heavy; harness 113 / 121 / 333 s) each 8 / 7 / 8 passed. Asserted and met: (P1W) ≤ 7.0e-15 every column (predicted ≤ 1e-14, met); (G2) 5.600000000000e-05 ×4; (i) ≤ 3.1e-14; (R1′) ≤ 1.7e-14; (E) ≤ 1.7e-11; (R3′) 4.9e-6 … 5.8e-6; x1: (F0′) 7.308e-11 / 6.601e-11 (predicted ≤ 1e-10, met), (X1c) 117 206 exactly, (F0-ctl) 1.020e-2 (predicted, met), (d6′) 1.000000; h0025: (R) 201 072 cells / 1 284 856 dofs exactly. Printed: x1 reproduces 371a / 371b (point SAR to ≤ 2e-9, θ₄₀ to ≤ 6.5e-7 — predicted ≤ 1e-8 for SAR, met; θ at the thermal solver's digit floor); (Live-h) 4.225e-5 / 6.663e-5 (predicted > 1e-6, met); (d6′) 1.000000 on both refined rungs (ruling Q5, printed); cells / tag-3 equal the predictions at every rung; ⟨θ⟩₄₀ order effect 1.0028 / 1.0005 / 1.0001 (predicted within ±0.5 % of 1.0028, met); **κ_s rose, not fell** (predicted 1.34 → 1.0–1.15; read 1.34 → 1.54 → 1.61 at d2x — NOT met, a prediction only). ρ(h) = 2.119 / 2.527 / 2.687 beside 2.119196, `_branch` not called. Logs: `20261012T100006Z_ANS-9-step1c-W0.log`, `20261012T100019Z_ANS-9-step1c-W1-x1.log`, `20261012T100223Z_ANS-9-step1c-W2-h005.log`, `20261012T100436Z_ANS-9-step1c-W3-h0025.log`; cache sweep `20261012T095822Z_ANS-9-step1c-cache-sweep.log` (no orphan).

| rung | cells | tag-3 | κ_s d2x / d1x | SAR_D d2x / d1x [W/kg] | θ_D d2x / d1x [K] | SAR_D order | θ_D order | D/C d2x / d1x | D/A d2x / d1x | h_p move d2x / d1x: SAR_D, θ_D, D/C, D/A | body s | GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| x1 | 117 206 | 537 | 1.3378 / 1.4326 | 3.0641 / 3.0940 | 0.7233 / 0.6202 | 0.9903 | 1.1663 | 435.63 / 10.217 | 8.2883 / 5.2810 | (comparand) | 109.1 | 18.18 |
| h005 | 130 427 | 8 498 | 1.5415 / 1.6063 | 2.9838 / 2.6854 | 0.7537 / 0.6875 | 1.1111 | 1.0964 | 519.45 / 2757.4 | 8.4125 / 5.4395 | −2.62 % / −13.2 %; +4.20 % / +10.8 %; +19.2 % / ×270; +1.50 % / +3.00 % | 117.9 | 21.03 |
| h0025 | 201 072 | 58 960 | 1.6112 / 1.6543 | 2.9797 / 2.4979 | 0.7569 / 0.7012 | 1.1929 | 1.0794 | 552.31 / 91.175 | 8.4823 / 6.3566 | −2.75 % / −19.3 %; +4.64 % / +13.1 %; +26.8 % / ×8.92; +2.34 % / +20.4 % | 328.9 | 46.50 |

At h0025 the three-way split prints the sheet step (d1x / `ans:9` d1) at SAR_D 0.9424 and θ_D 1.0054 (row 2 0.99995), and the order step (d2x / d1x) at SAR_D 1.1929 and θ_D 1.0794. SAR_C stays a near-null (d2x 5.4e-3 W/kg at h0025), so D/C swings with C. No reading of P21 is made here; the weekly reads where θ_D's order effect landed.

## §7 row moved verbatim 2026-10-12 (rotate_plan_archive.py chunks)

 **EM–thermal coupling against AEDT HFSS + Ansys thermal (Icepak FEA): SAR from the coil heats the gelled-saline phantom** (operator directive 2026-09-30, interactive session; Phase 8's first external check — every Phase 8 gate so far is a closed form). **Gel first, by the operator's instruction** (the configuration most relevant to the operator's research): no perfusion. Progression (one-way gel, two-way gel, perfusion) and resources in the chunk file. *(history and §9 briefs: `docs/planning/chunks/ANS-9.md`)* | 🟡 **Step 1 adjudicated 2026-10-10 (weekly): AGREE on the EM → thermal hand-off — rows 1 / 2 / 3 / 6 within 0.8 % / 0.5 % / 1.4 % of the operator's AEDT 2026 R1 HFSS + IcepakFEA run (bands 5 %); probes A–C within 2 % with the B > D > C ordering held; DISAGREE at probe D (ours 11 % below AED, band 10 %) — known-issues 2026-10-10, diagnosis step 1b (§10 P21); rows 5 / 7 records. AED's stepping (one load step per report time, 13 fixed substeps) and the First Order vs degree-1 EM mismatch are stated with the verdict; numbers in `docs/private/`.** Commissioned 2026-10-03 by the weekly (SPEC tracked: `examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md`, one-way gel heating on `ANS-2`'s 10 MHz EM half). **Done-when (daily 2026-10-04, rule (j); ratified by the 2026-10-05 weekly, which also added the SPEC's probe-volume note after `MAT-19` step 2):** the runnable half — `ans:9`, SPEC rows 1–7, rows 1 / 3 / 6 asserted at `MAT-15`'s imported bands — ⇒ 🟡; the weekly's adjudication of the operator's run ⇒ ✅. **🟡 2026-10-04 (09:00 slot, §9 item 218): runnable half landed.** `ans:9` is green at `-n 2`, 175 s. Row 1 is within 7.3e-15 of 1 W, (R3) is 5.7543e-6, and the row 6 energy balance closes to 7.09e-15 (surface 164.50 J); controls d6 1.000 and s6 0.1312. Step narrative in the chunk file. **Step 1b (a) landed 2026-10-11 (§9 item 371a): ρ = (D/C)_d2x / (D/C)_d1 = 2.119, branch "no branch"** — the 10-14 weekly reads it; local SAR at A–E at 1 W on d1 / d1x / d2x, (P1W) ≤ 8.4e-15, item 192's S reproduced to 7.3e-11; SAR_D moves +15.6 % d1 → d2x, almost all of it in the sheets + mesh step (d1 → d1x +16.7 %, order −1.0 %); C is a near-null (0.013 / 0.30 / 0.007 W/kg) that dominates D/C (`20261011T184526Z_ANS-9-step1b-W2-d2x.log`). **Step 1b (b) landed 2026-10-11 (§9 item 371b):** `ans:9`'s thermal half on item 192's 537-cell coarse phantom, sources d2x and d1x — (R1′) ≤ 4.4e-16, (E) ≤ 3.9e-12, (R3′) 4.9e-6 / 5.2e-6, (d6′) 1.000; θ₄₀ at D moves +16.6 % with the order on the same thermal mesh (0.620 → 0.723 K), the other probes and ⟨θ⟩₄₀ ≤ 1.9 %; readings on a coarse phantom, `ANS-9`'s column not re-recorded (`20261011T185145Z_ANS-9-step1b-W2b-d2x-thermal.log`). **Step 1c landed 2026-10-12 (§9 item 396):** θ_D order effect x1 / h005 / h0025 = 1.166 / 1.096 / 1.079; point SAR_D order effect 0.990 / 1.111 / 1.193; route B tied to item 192 at x1 (7.3e-11); every identity green per rung; the h0025 d2x column is a candidate only (`20261012T100436Z_ANS-9-step1c-W3-h0025.log`). 

### §9 item 410 — 2026-10-13 (§10 chain P24; chain step, solve)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* Licence: §10 chain P, step P24 (the 2026-10-12 weekly), executed as written. **No band, tolerance, Done-when or target is set or moved; every band is imported by name. Our values only: no slot opens `docs/private/`, `aed_results/`, `aed/` or `COMPARISON_private.md`, and the AED cells stay blank.**

`ANS-9` step 1d: `ans:9`'s column re-recorded at degree 2 on exact rectangles with the h_p = 0.0025 m phantom (route B, item 396's h0025 window as the EM half, the thermal half on that phantom), every SPEC row 1–7 as the primary column, the degree-1 column kept and labelled (implementer; heavy, `-n 8`, complex; ≈ 12 measured / 60 predicted min: 15 fixed + 25 code + ≈ 14 of windows + spare; depends on nothing).

- **Why one window holds it.** Item 396's W3 measured the whole EM half at this fixture: build 40.87 s, degree-2 sweep 261.47 s, thermal 10.16 s, body 328.91 s, 46.50 GiB summed, harness 333 s (`20261012T100436Z_ANS-9-step1c-W3-h0025.log:2054`, `:2109`). This item drops W3's degree-1 sweep (13.99 s) and adds `ans:9`'s snapshot run, half-step, BE step, doubled-source run, row 7 and the XDMF (≈ 15–30 s on 58 960 cells). **EM + thermal fit one 590 s window; no split.**
- **Gitignored inputs:** none.

**Precondition findings (read before coding).**
1. **The writer deletes the weekly's verdict.** The script rewrites `COMPARISON.md` wholesale (`examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/09_coil_driven_gel_heating_10MHz.py:426`). The 2026-10-10 weekly appended a hand-written "## Verdict (weekly review 2026-10-10; AED values private)" section (`COMPARISON.md:36–67`). A bare re-run erases it. This is item 400's `ans:7` defect in another file, and the same fix applies (below).
2. **Two of the script's degree-1 anchors do not transfer.** (R0) compares the re-run against `MAT-15`'s gated `T_star` (`09_…py:246–252`), which exists only on the ragged fixture. The V_sub check is against `V_SUB_RECORD_M3`, the 58 866-cell phantom's record (`tests/validation/test_mat15_coil_driven_heating_phantom.py:96`, `V_BAND` 1e-10 at `:98`). Route B's phantom has 58 960 cells. W3 printed its V_sub as 2.259978941e-04, which is equal at 10 digits (`…W3-h0025.log:1993`) but not measured at 1e-10. On the new column both are **printed**: (R0) becomes a reproduction print against W3, and V_sub is printed to 13 digits beside the record.
3. **The script is imported by path elsewhere.** `MAT-19`'s `_load_spec_probes` (`tests/validation/test_mat19_probe_coil_phantom.py:96–105`) executes the script as a module, and so do item 371a's and item 396's modules through `_column_readout` (`tests/validation/test_ans9_far_side_sar.py:133`). `SPEC_THERMAL` stays byte-identical. The top level stays import-only, with all compute inside `main()`.

**Code (script edit, a disclosed rule (c) change to an example; no `src/` edit; no edit to any test module).**
- **The EM half is route B, imported and never restated.** Item 396's body steps 1–4 (`tests/validation/test_ans9_far_side_sar_refined.py:161–200`):
  - `_build(True, phantom_resolution=RUNGS["h0025"], leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`;
  - `build_four_port_sweep(built=…, build_only=True, degree=2)`;
  - `_sweep(case, 2)` (`…refined.py:138`), with `_recip` (`:149`);
  - `_column_readout("h0025-d2x", f.e_complex, f.sigma_field, case["cell_tags"], comm)` (`test_ans9_far_side_sar.py:125`) for G and (P1W);
  - `_tag_submesh` + `_robin_tags`, then q_h exactly as `_thermal_readout` builds it (`test_ans9_far_side_sar.py:151–153`): `_to_sub(*_cell_average_dg0_tag(msh, ct, PHANTOM_CELL_TAG, G·½σ·inner(e, e)), sub, parent_ids)`.
  - Imports by name: `RUNGS` (`tests/validation/test_geo27_one_gram_rung_degree2.py:86`), `_build`, `build_four_port_sweep`, `GATED_WIDTH_FRACTION`, `_sweep`, `_recip`, `_column_readout`, `_to_sub`, `_robin_tags`, `_tag_submesh`, `_cell_average_dg0_tag`, and the bands below.
- **Rows on the new column.** The script's own row code (`09_…py:206–338`) is refactored into a function of `(sub, ft, q_h, v_sub)` and applied to the degree-2 source unchanged. That gives rows 1–6, (R1), (R3), (R6), (d6) and (s6) with the same snapshot, half-step and BE schedule.
  - **Row 7 uses `ANS-2`'s own arithmetic.** Call `_excitation` (`02_birdcage_coil_driven_sar_10MHz.py:271`) on the adaptor `{"driven": "P1", "source_voltage_v": complex(P1 spec's drive_voltage_v), "currents": {"P1": r2.excitation_results["P1"].responses["P1"].current_a}}`.
  - That current is the package's `sheet_terminal_current` (`src/fem_em_solver/ports/lumped.py:571`), the quantity `_solve_driven` returns (`tests/validation/test_birdcage_b1_plus_map.py:189–194`). The lumped route's source voltage is the spec's `drive_voltage_v`, default 1 V (`lumped.py:400`, `:422`).
- **The degree-1 column is kept as a labelled record, not re-solved.** A module constant `D1_RAGGED_COLUMN` copies every value of the tracked `metrics.json` at HEAD `59809ba7` (rows 1–7, the anchors, cells 199 920 / 58 866, `-n 2`). Its provenance is `20261004T142653Z_ANS-9-step1a-W1.log`. Label: "N1curl degree 1, ragged pre-10-03 sheets, `MAT-15` step 1 fixture, 2026-10-04 run, adjudicated 2026-10-10".
  - The 10-04 run's (R0), (R6), (d6) and (s6) are not re-asserted; disclose that in the outcome.
  - **Alternative the review may prefer:** keep the degree-1 path behind a `FEM_EM_ANS9_COLUMN=d1` knob as its own window. That costs ≈ 140 s more (item 371a's W1 ran `_run_heated` in 136 s at `-n 8`) and re-asserts 1a's anchors.
- **Artifacts.**
  - **`metrics.json`.** The top-level keys keep their names and now carry the degree-2 column. Metadata: `"em_element": "N1curl degree 2"`, `em_cells` 201 072, `thermal_cells` 58 960, `"fixture": "route B, exact rectangles, κ off, h_p 0.0025 m"`, `mpi_ranks` 8. Add `"columns": {"d1_ragged": D1_RAGGED_COLUMN}`.
  - **`COMPARISON.md`.** Columns: Row | Quantity | Ours (degree 2, exact, h_p 0.0025) | Ours (degree 1, ragged; 2026-10-04) | AED | ours vs AED. The AED cells are `—`, and the SPEC readout is copied verbatim.
  - **The writer becomes `_write_comparison(case_dir, metrics)`.** It reads the existing file and re-appends everything from the first line starting `## Verdict` **byte-for-byte**.
  - **XDMF.** Same basename, written from the degree-2 run. Every `Function` is named.
  - **Guide** `09_coil_driven_gel_heating_10MHz.md`:
    - run line `./scripts/run_examples.sh -e ans:9 -n 8 -t 590`, with ≈ 6 min stated;
    - "degree 2 primary, degree-1 column carried";
    - Not claimed: one rung, no h-statement on the thermal side.
  - **Setup figure:** not regenerated (geometry only; opt-in; disclosed).
- **New unit test `tests/unit/test_ans9_comparison_writer.py`.** It loads the script by `spec_from_file_location` (item 400's pattern) and covers three things:
  - (K1) `D1_RAGGED_COLUMN` equals the committed `metrics.json`'s degree-1 values exactly (top-level keys before W2, `columns.d1_ragged` after);
  - (R-keep) render into `tmp_path` from copies of the committed files; the tail from `## Verdict` is byte-identical to the committed tail (lengths and heading count 1 printed);
  - (N-tail), asserted: render from a copy with the Verdict section removed, and the output's `## Verdict` count is 0. This shows the tail is read, not hard-coded.

**Anchors (asserted; every band imported).**
- **Route-tie window W1 (x1).** Item 396's module re-run as committed (unchanged since `e3edf5de`):
  - (F0′) d2x S against `ANS4_EXACT_STRIP_S_P2_X1["10"]` and d1x S against `ANS4_EXACT_STRIP_S_P1_RECORDS["10"]`, each ≤ `STEP3A_RECORD_RTOL` (`tests/validation/test_ans4_resolution_ladder.py:198`). Predicted ≤ 1e-10; W1 of item 396 read 7.308e-11 / 6.601e-11 (`20261012T100019Z_ANS-9-step1c-W1-x1.log:2039`).
  - (X1c) cells within `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`) of `ANS4_EXACT_STRIP_CELL_RECORD`. Predicted 117 206 exactly.
  - (P1W), (G2), (i), (R1′), (E) and (R3′), as that module asserts them (`…refined.py:286–307`).
- **The column window W2 (h0025, the script).**
  - (P1W) ≤ `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`), all five points valid. Predicted ≤ 1e-14; W3 read 6.994e-15 (`…W3-h0025.log:1986`).
  - (G2) the four strip areas equal `NOMINAL_STRIP_AREA` within `STRIP_BAND`.
  - (i) ‖S − Sᵀ‖/‖S‖ ≤ `RECIPROCITY_BAND` (`tests/validation/test_port_lumped_sheet_sweep.py:86`). W3 read 1.335e-14 (`:1984`).
  - (R) cells within `CELL_COUNT_BAND` of `I14_OFF_CELLS`, and d2 dofs == `I14_OFF_DOFS` (`tests/validation/test_imp3_tip_ladder_1g.py:94`). W3 read 201 072 / 1 284 856 (`:2080`).
  - (R1) = (R1′): |∫_sub q_h − 1 W| ≤ `B_BAND` (`tests/validation/test_mat10_sar_to_heat_sphere.py:61`). Predicted ≤ 2e-14; W3 read 1.688e-14 (`:1993`).
  - (E): every run's `max_balance_residual` ≤ `E_BAND` (`tests/validation/test_mat12_transient_sar_heating_sphere.py:76`). W3 read ≤ 9.9e-12 (`:1995`).
  - (R3) = (R3′): ≤ `S2_BAND` (`…mat12…py:75`). Predicted 5.6446e-6 (`:1996`).
  - (R6): the energy balance with row 6 built as the independent per-(half-)step surface sum, ≤ `E_BAND`. Predicted ≤ 1e-12 (1a read 7.09e-15).
  - **Why (R6) is asserted beyond P24's list:** the §7 Done-when's runnable-half clause reads "rows 1 / 3 / 6 asserted at `MAT-15`'s imported bands" (rule (j)). The primary column has to meet it.
  - (A0): the SPEC entries equal their imports (unchanged), `r2.port_ids[0] == "P1"`, and the three `is` checks.
  - W0 / W0b: (K1), (R-keep), (N-tail).

**Printed (predicted, never asserted).**
- **Rows (predictions from W3's h0025-d2x lines, the same fixture and route):**
  - row 1: 1 W;
  - row 2 ⟨θ⟩: 0.318123 / 0.614370 / 0.892184 / 1.153415 K (`:1994`);
  - row 3: slope 1.058563e-3 K/s, beside 1 W/(ρcV) = 1.058569e-3;
  - row 4 θ₄₀: A 0.889137, B 3.514465, C 0.683153, D 0.756873, E 0.959953 K (`:1997`);
  - row 5: max ≈ 4.3–4.5 K at r ∈ [0.020, 0.030] m, \|az\| ≤ 30° (1a: 4.457206);
  - row 6: 1254 − ρcV·⟨θ⟩₄₀ = **164.40 J** (surface fraction 0.1311);
  - row 7: G·2.5e-3 W = **45 844.81 W** at G = 1.833792433645e7 (`:1986`; +0.92 % from the degree-1 45 428.08), with supplied and accepted beside.
- **Reproduction against W3** (G, rows 2 / 3 / 4, SAR at A–E): predicted equal to ≤ 1e-6 relative. Item 396's x1 reproduced 371b's θ₄₀ to ≤ 6.5e-7, which is the thermal solver's digit floor; MUMPS at `-n 8` is not bitwise.
- **Beside each row:** the degree-1 value and the ratio degree 2 / degree 1.
- **Also:** V_sub to 13 digits against `V_SUB_RECORD_M3`, κ_s, the first half-step deficit, and the PRICE line (build / sweep / thermal / body s, `_peak_rss_gib()`, ranks).

**Negative controls (rule (e)).**
- **(F0-ctl), x1 (W1), asserted > `STEP3A_RECORD_RTOL`; the factor is asserted.** The same comparison on the same fixture read 1.020e-2 (`…W1-x1.log:2039`), ≈ 1e4× the band.
- **(d6), W2, the doubled source.** (R6) against T·∫q_h misses by **predicted 1.000** (the ceiling is exactly 1, by linearity). It is asserted only > `D_CONTROL_FLOOR` (`test_mat15_coil_driven_heating_phantom.py:102`; 1a's sign-only form).
  - The factor is *predicted*. The nearest measurement on this fixture is (d6′), the (E′) form, at 1.000000 (`…W3-h0025.log:1998`), which is a different comparison.
- **(s6), W2, surface term dropped.** (R6) misses by row 6 / (T·∫q), **predicted 0.1311** (its own ceiling; 1a read 0.1312). Asserted only > `E_BAND`.
- **(N-tail), W0:** asserted, as above.

**Windows.** Foreground. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Bash `timeout: 660000`. **Two heavy windows never run at once.**
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Then the any-orphan cache sweep as one bare call, label `ANS-9-step1d-cache-sweep` (daily-review.md rubric item 4).
- **Pre-censuses** (with item 218's commands): `ANS-9-step1d-pre-figures` and `ANS-9-step1d-pre-docrefs`.
- **W0 (smoke, `-n 1`), before W2:**

  `scripts/testing/run_and_log.sh ANS-9-step1d-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/unit/test_ans9_comparison_writer.py -s -v --tb=short'"`
- **W1 (heavy, `-n 8`, the x1 route tie).** Predicted ≈ 113 s, item 396's W1 (`…W1-x1.log:2078`). This is 396's W1 command with a new label and raw file:

  `scripts/testing/run_and_log.sh ANS-9-step1d-W1-x1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS9_1C_RUNG=x1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_far_side_sar_refined.py -s -v --tb=short > /workspace/logs/ans9-step1d-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans9-step1d-w1-raw.log; cat /workspace/logs/ans9-step1d-w1-raw.log; exit \$rc'"`
- **If W1 is red, stop** (see Negative result). Otherwise run pgrep again; it must read 0.
- **Dry run:** `scripts/testing/run_and_log.sh ANS-9-step1d-dryrun "./scripts/run_examples.sh -e ans:9 -n 8 -t 590 --dry-run"`. It must print `(complex build)` and the `09_…py` path.
- **W2 (heavy, `-n 8`, the column).** Predicted harness 330–400 s, RSS ≤ 46.5 GiB summed. It is the runner's line in durable capture:

  `scripts/testing/run_and_log.sh ANS-9-step1d-W2-h0025 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 8 python3 examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/09_coil_driven_gel_heating_10MHz.py > /workspace/logs/ans9-step1d-w2-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans9-step1d-w2-raw.log; cat /workspace/logs/ans9-step1d-w2-raw.log; exit \$rc'"`

  The script exits non-zero on any red anchor.
- **W0b:** W0's command with label `ANS-9-step1d-W0b`, after W2. (K1) and (R-keep) then read the new committed-to-be files.
- **Pin**, after `git add …/metrics.json`: item 218's W2 command, label `ANS-9-step1d-pin`. The count is unchanged, since the path is already listed (`tests/unit/test_doc_reference_exit_codes.py:477`).
- **Post-censuses** `ANS-9-step1d-post-figures` / `ANS-9-step1d-post-docrefs`. Write the prediction first: runnable unchanged, figures `broken = 0`, docrefs `exit != 1`.
- **Time rule:** if W2 has not started by minute 45, commit the edits and the W0 / W1 logs and mark the §9 entry "⬜ W2 remaining (script at `<sha>`)". The entry stays runnable, so rule (d) is met without BLOCKED (item 396's Q3 form).

**Traps already paid for.**
- **COMPARISON.md overwrite** (finding 1). Never run W2 before W0's (R-keep) is green.
- **`_four_port_rung` has no `phantom_resolution`**, and `build_four_port_sweep(built=…)` ignores sizing keywords, so h_p goes to `_build` only. Pass `degree=2` on both the build and the sweep. Never pass `width_correction_kappa`: κ is off (2026-10-02 ruling).
- **Never call `_solve_driven` at h0025, degree 2.** One extra factorisation costs ≈ 230–260 s and breaks 590 s. P1's field and current come from the one-factorisation sweep with `keep_fields=True`.
- **Free the sweep before the thermal runs.** After taking q_h and row 7's current: `del r2`, `gc.collect()`.
- **Point reads.** The N1curl point read uses `owner="lowest_original_cell"` on the parent (inside `_column_readout`). θ on the `create_submesh` mesh uses the **default** owner: `lowest_original_cell` raises "Mixed topology unsupported" there (item 218). Points are identical on every rank.
- **Snapshots and row 6.** `_transient` snapshot 1 lands after two half-steps, and θ^{½} needs its own run. Never reconstruct row 6 from the balance itself (rule (h)).
- **XDMF.** Name every `Function` before the write. No dots in the basename. `write_xdmf_with_tags` returns `None` off rank 0.
- **Reductions.** Every one takes `np.real` and is allreduced on every rank; no collective sits inside `rank == 0`. Cell counts come from `size_global` / `_global_tag_cell_count`.
- **Imports.** Import no `test_*` name and no pytest fixture. The 1c module's import-time knob read is inert with `FEM_EM_ANS9_1C_RUNG` unset (W2 leaves it unset).
- **Wrong rank count.** Do not run `ans:9` at `-n 2`: the degree-2 h0025 solve was priced at `-n 8` only (46.5 GiB summed). The guide's old `-n 2 -t 480` would not fit.
- **The pin test reads `git ls-files`.** Stage `metrics.json` first.
- **Durable capture is copied verbatim**, trailing `; exit $rc` included. Use `-k 30`. No pipe inside the quoted command.
- **After a kill:** run `pgrep -c python3`, then the cache sweep, then park. Status 124 is a price record, not a retry. Never raise the timeout or `-n`.
- **Rule (i):** no script or test edit after its last window. An edit after W2 re-runs W0 → W2.

**Scope.** `ans:9`'s our-side column at degree 2 on the exact-rectangle F-small fixture, κ off, 10 MHz, drive P1, one phantom rung (h_p = 0.0025 m), P1 thermal, CN Δt = 31.35 s. The degree-1 column is carried as a record. **Not covered:**
- any AED value, verdict or adjudication (the 10-14 weekly's, private);
- an h-statement on either half;
- F–I;
- σ(T) (item 425);
- perfusion and tissue;
- an absolute temperature;
- the setup figure.

**Negative result.** Report the measurement and stop. Never widen a band.
- **(K1) red:** a transcription error. Fix it in-slot before W2.
- **(R-keep) or (N-tail) red:** a writer defect. Fix it before W2.
- **(F0′) or (X1c) red at x1:** route B is no longer item 192's fixture on this image. Stop before W2, park on `attempt/ANS-9-step1d-<ts>` and mark BLOCKED (rule (d)) with the unblock condition "route-B / item 192 difference located". Add a known-issues bullet.
- **(R), (G2) or (i) red:** the fixture is not item 396's h0025. Stop and mark BLOCKED.
- **(P1W), (R1), (E) or (R6) red:** a source-transfer or surface-sum defect. Park and mark BLOCKED.
- **(R3) red with (E) green:** a time-truncation reading. Report κ_s and stop. `S2_BAND` is never widened.
- **(d6) or (s6) not missing:** a blind control. Stop.
- **A reproduction print against W3 off by more than 1e-6:** report it. If any row moves ≥ 1e-3, flag it for the weekly, which adjudicates on these numbers.
- **W2 killed at 590 s:** a price record. The review re-poses it (for example a JSON hand-off of q_h between an EM window and a thermal window).
- No `PREDICTED RED` line.

**Status it can move.**
- **`ANS-9` row.** The state line gains "step 1d landed: `ans:9`'s primary column is degree 2 on exact rectangles, h_p 0.0025 m (route B), rows 1 / 3 / 6 asserted, degree-1 column kept and labelled; `<log>`".
  - The runnable-half clause of the Done-when ("`ans:9`, SPEC rows 1–7, rows 1 / 3 / 6 asserted at `MAT-15`'s imported bands ⇒ 🟡") is met again on the new column. The glyph stays 🟡.
  - ✅ is the weekly's adjudication. The 10-14 weekly re-adjudicates every row against the operator's existing run. Rule (j) is matched, and no row amendment is needed.
- **Known-issues.** The 2026-10-10 entry ("`ANS-9` step 1 adjudication … far-side probe D …") gains a "Measured (step 1d, §9 item 410)" bullet with rows 2–4 at degree 2 against degree 1. It retires on the weekly's adjudication, not here.
- **§10 P24** gets a DONE line.

**Spawn line:** "Use the implementer agent to execute `ANS-9` §9 item 410 per `docs/planning/chunks/ANS-9.md` (read the brief and the newest Review ruling beneath it first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; two heavy windows never run at once; `pgrep -c python3` reads 0 before each heavy window. No `src/` edit and no edit to any existing test module; the `COMPARISON.md` Verdict section survives byte-for-byte. Never open docs/private/, aed_results/, aed/ or COMPARISON_private.md, and write no AED number anywhere. You are the executor: do not spawn agents. Report a red and stop; never widen a band."

**Outcome 2026-10-13 (04:30 slot, implementer): DONE — every anchor and control green, first run of each window.** Script `09_coil_driven_gel_heating_10MHz.py` refactored (route B EM half imported inside `main()`, `_thermal_rows`, `_write_comparison`, `D1_RAGGED_COLUMN`); new `tests/unit/test_ans9_comparison_writer.py`; no `src/` edit, no existing test module edited, `SPEC_THERMAL` byte-identical. Windows: W0 (2.7 s) 7 passed — (K1) 21 shared keys equal, (R-keep) tail 1 938 B in / out, Verdict count 1 / 1, (N-tail) count 0; import check (3 s); W1 x1 (`-n 8`, harness 118 s) 8 passed / 1 skipped — (F0′) 7.337e-11 / 6.597e-11 (predicted ≤ 1e-10, met), (F0-ctl) 1.020e-2 (asserted > band, met), (X1c) 117 206 exactly, (d6′) 1.000000; dry run prints `(complex build)` and the `09_…py` path; **W2 h0025 (`-n 8`, harness 330 s, predicted 330–400 s, met)** Status 0; W0b 7 passed — (K1) on `columns.d1_ragged`, 23 keys, (R-keep) tail byte-identical on the new files; pin 20 passed; post-censuses figures `broken=0`, docrefs exit 2 (stale 90, unchanged from the pre-census; predicted `exit != 1`, met). The 10-04 run's (R0), (R6), (d6) and (s6) are not re-asserted on the carried degree-1 column (disclosed). (R0) and V_sub become prints, per the brief.

Asserted on W2 (bands imported): (P1W) 1.377e-14 ≤ `EXACT_IDENTITY_RTOL`, all five points valid; (G2) 5.600000000000e-05 ×4 (rel ≤ 1.6e-15); (i) 1.321e-14; (R) 201 072 cells, 1 284 856 d2 dofs exactly; (R1) 1.121e-14 ≤ `B_BAND`; (E) CN 2.731e-14, half-step 2.731e-14, BE 6.864e-12 ≤ `E_BAND`; (R3) 5.6446e-6 ≤ `S2_BAND` (predicted 5.6446e-6, met); (R6) 2.244e-15 ≤ `E_BAND` (predicted ≤ 1e-12, met); (A0) all green, `port_ids[0] == "P1"`. Controls: (d6) 1.000000 > `D_CONTROL_FLOOR` (predicted 1.000, met); (s6) 0.1311 > `E_BAND` (predicted 0.1311, met). Printed: V_sub 2.2599789411112e-04 vs the 58 866-cell record, rel 3.532e-13; reproduction against item 396's W3 — G 5.2e-11, SAR A–E ≤ 5.5e-9, row 2 ≤ 7.2e-7, row 4 ≤ 2.8e-7 (predicted ≤ 1e-6, met; no row moved ≥ 1e-3); first half-step deficit 5.0262e-3. PRICE: build 41.31 s, sweep 272.12 s, thermal 3.61 s, body 327.53 s, peak RSS 46.87 GiB summed, 8 ranks.

| row | quantity | degree 2 (primary) | degree 1 (carried) | d2 / d1 |
|---|---|---|---|---|
| 1 | power, W | 1 − 1.1e-14 | 1 − 2.9e-15 | 1.000000 |
| 2 | ⟨θ⟩ 313.5 / 627 / 940.5 / 1254 s, K | 0.318123 / 0.614370 / 0.892184 / 1.153415 | 0.318096 / 0.614315 / 0.892104 / 1.153314 | 1.000085 / 1.000090 / 1.000089 / 1.000088 |
| 3 | slope, K/s | 1.058563e-3 | 1.058563e-3 | 1.000000 |
| 4 | θ₄₀ A / B / C / D / E, K | 0.889137 / 3.514465 / 0.683153 / 0.756873 / 0.959953 | 0.905798 / 3.564854 / 0.685305 / 0.697444 / 0.958639 | 0.9816 / 0.9859 / 0.9969 / 1.0852 / 1.0014 |
| 5 | max θ₄₀, K (location m) | 4.404696 (+0.0279, −0.0005, +0.0008) | 4.457206 (+0.0275, +0.0016, −0.0003) | 0.9882 |
| 6 | surface heat, J | 164.4017 | 164.4978 | 0.9994 |
| 7 | incident power for 1 W, W | 45 844.81 | 45 428.08 | 1.0092 |

Logs: `20261013T093058Z_ANS-9-step1d-cache-sweep.log` (no orphan), `…093246Z_…-pre-figures.log`, `…093248Z_…-pre-docrefs.log`, `…093547Z_ANS-9-step1d-W0.log`, `…093604Z_ANS-9-step1d-importcheck.log`, `…093614Z_ANS-9-step1d-W1-x1.log`, `…093819Z_ANS-9-step1d-dryrun.log`, `…093826Z_ANS-9-step1d-W2-h0025.log`, `…094427Z_ANS-9-step1d-W0b.log`, `…094433Z_ANS-9-step1d-pin.log`, `…094446Z_…-post-figures.log`, `…094448Z_…-post-docrefs.log`. No PREDICTED RED line in the brief; none fired.

**Same-commit duties.** Mark §9 item 410 DONE, BLOCKED or "⬜ W2 remaining". Commit the script, guide, `metrics.json`, `COMPARISON.md`, the unit test and every harness log together. Add the row state clause, the known-issues bullet, the §10 P24 DONE line, an outcome paragraph here (rows table: degree 2 / degree 1 / ratio; anchors; controls; PRICE) and the journal entry.

### §9 item 425 — 2026-10-13 (§10 chain P27; chain step, solve; depends on item 410)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* Licence: §10 chain P, step P27 (the 2026-10-12 weekly). **No band, tolerance, Done-when or target is set or moved; every band is imported by name. Skip if item 410 did not land.** Our values only: no slot opens `docs/private/`, `aed_results/`, `aed/` or `COMPARISON_private.md`.

**Precondition findings (the review rules before queueing).**
1. **"`MAT-14`'s route" names the wrong chunk.** `MAT-14` is the SAR-driven transient Pennes on the sphere and carries no σ(T) (§7 `MAT-14` row; `docs/planning/chunks/MAT-14.md:3`). The coil-driven σ(T) route is **`MAT-5` step 3's** Picard map: `tests/validation/test_mat5_coil_driven_sigma_feedback.py`, with `_Fixture.phi` at `:196–205`, σ written through `_solve_driven(…, sigma_field_hook=…)` and `_sub_to_parent`, and the contraction identity (K) at `:416` / `:474`. P27's "contraction identity at `MAT-5` step 3's imported band" agrees with this reading. **The brief reads "`MAT-14`" as "`MAT-5` step 3"; the review confirms it.**
2. **As posed (degree 2 on P24's h0025 fixture), P27 does not fit heavy.**
   - The contraction identity re-solves E at every Picard iterate and power-iteration step: N = 8 + M = 5 = 13 `_solve_driven` calls (`:99–100`, `:327–390`), plus 3 at α = 0.
   - A degree-2 factorisation at h0025 costs ≈ 230–260 s at `-n 8`: the four-port sweep took 261.47 s with one factorisation (`20261012T100436Z_ANS-9-step1c-W3-h0025.log:2054`). σ(T) changes the matrix, so no factorisation can be held.
   - That is ≈ 16 × 245 ≈ 65 min, against a 590 s window and the 1200 s heavy ceiling. Splitting across windows would need θ handed between processes at `-n 8`, which no helper does.
   - **The letter-faithful form is an `xl` entry, the weekly's to commission:** ≈ 2 h of degree-2 factorisations, inside xl's 4 h.
   - **Re-pose written below, under §9 rule (6), for the review to adopt or strike:** both one-way and two-way at **degree 1 (d1x) on P24's mesh**, built by route B exactly as item 410. W3 measured the degree-1 four-port sweep, factorisation included, at 13.99 s (`:2054`). The order moves ⟨θ⟩₄₀ by only 1.0001 at h0025 (item 396's table), so a degree-1 σ(T) *difference* prices the two-way effect. Item 410's degree-2 one-way column is printed beside.
3. **`solve_transient_conduction` cannot restart from a field.** It takes only a uniform scalar `t_init` (`src/fem_em_solver/thermal/transient.py:103`, `:123`). A staggered transient two-way (σ updated at the report times) is therefore built by **superposition** over `source_schedule` (`:112`, `:136–143`). The scheme is linear in q with θ(0) = 0 (t_init = T∞), so a run with source q_j switched on only in segment j, summed over j, is the staggered run. No `src/` edit; a field `t_init` would be a rule (c) re-pose.

**`ANS-9` step 2 (re-posed): two-way σ(T) on P24's fixture at degree 1. The steady contraction identity is asserted; the transient one-way vs two-way ⟨θ⟩(t\*) and θ_D are printed as a record** (implementer; heavy, `-n 8`, complex; ≈ 12 measured / 65 predicted min: 15 fixed + 25 code + ≈ 15 of windows + spare; depends on item 410).

**Code.** One new module, `tests/validation/test_ans9_two_way_sigma.py`. No `src/` edit and no edit to any existing module.
- `@complex_only`. Knob `FEM_EM_ANS9_STEP2` ∈ {unset, `alpha0`, `steady`, `transient`}. Unset ⇒ the gated tests skip; any other value ⇒ `ValueError` (item 396's `_knob` pattern).
- **`_RouteBFixture(m5._Fixture)`** overrides `__init__` only, with `import tests.validation.test_mat5_coil_driven_sigma_feedback as m5`.
  - It copies `m5._Fixture.__init__` (`:120–142`) with one change: the build is `built = _build(True, phantom_resolution=RUNGS["h0025"], leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)` and then `self.sweep = build_four_port_sweep(built=built, build_only=True, degree=1)`.
  - Every other method (`em`, `phi`, `q_h`, `steady`, `sigma_sub`, `theta`, `records`, `ip`, `norm`, `p1`) is inherited unchanged.
  - Print cells and tag-3 (predicted 201 072 / 58 960, `…W3-h0025.log:1953`).
- **`_coupled(fx)`** is `m5`'s `coupled` body (`:320–442`) as a plain function of a fixture, verbatim except for the fixture construction. That is a disclosed copy: a pytest fixture cannot be called directly. Add θ_D at A–E (default owner on the sub) for iterates 1 and 8.
- **Imports by name:**
  - from `m5`: `K_BAND` (`:96`), `KCONV_BAND` (`:97`), `M0_FLOOR` (`:98`), `N_PICARD`, `M_POWER`, `FD_RMS_K`, `R0C_BAND` (`:95`), `_owned_cell_dofs`, `_real`, `_rss_summed_gib`;
  - `ALPHA` (`tests/materials/test_mat5_step2_lumped_root.py:54`);
  - `L0_BAND`, `S_BAND`, `D_BAND`, `_sub_to_parent` (`tests/validation/test_mat5_em_coupled_sphere.py:104, :108, :107, :117`);
  - `B_BAND`, `C_BAND` (`test_mat10_sar_to_heat_sphere.py:61–62`); `E_BAND`;
  - route B's names as in item 410; `CELL_COUNT_BAND`, `I14_OFF_CELLS`.
  - **Never** `alpha0`, `coupled` or any `test_*` name: importing a fixture registers it.
- **Superposition stagger, `_stagger(fx, n_seg, alpha_law)`.** Segments j = 0 … n_seg − 1 cover [jT\*/n_seg, (j+1)T\*/n_seg]. Segment boundaries fall on full steps 10j (n_seg = 4) or 5j (n_seg = 8).
  - For each j: θ(t_j) = Σ_{i<j} θ_i(t_j), from each earlier run's snapshot at step 40j/n_seg. Then σ_j = σ(Π₀(T∞ + θ(t_j))) and E_j = `fx.em(σ_j)` (segment 0: σ₀, i.e. the one-way solve; G fixed from it and never renormalised, as `m5`). Then q_j = `fx.q_h(E_j)`.
  - θ_j comes from `_transient(sub, ft, q_j, dt=T_STAR/N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, snapshot_steps=(10, 20, 30, 40), source_schedule=χ_j)` with χ_j(t_mid) = 1 on segment j and 0 otherwise. The two Rannacher half-steps fall in segment 0.
  - The two-way θ(t) = Σ_j θ_j(t).
- **W0 coding test (`test_ans9_step2_coding`, ungated):**
  - the knob parser;
  - the partition: Σ_j χ_j(t_mid) = 1 at every (half-)step midpoint of the Rannacher schedule for n_seg = 4 and 8, which is pure Python.

**Anchors (asserted; every band imported).**
- **W1 `alpha0`:**
  - (R) cells within `CELL_COUNT_BAND` of `I14_OFF_CELLS`: this is P24's mesh.
  - (R0a) hook path = P1 path: P and q ≤ `L0_BAND` (`m5:448–449`).
  - (D0) the hook touches only tag 3, exact 0 (`m5:468`).
  - (R0c) the α = 0 map residual ≤ `R0C_BAND` (`m5:463`).
  - (R1′) ∫q₀ = 1 W ≤ `B_BAND`.
  - (E) every CN run's `max_balance_residual` ≤ `E_BAND`. The residual is normalised by the full-on ∫Q under a schedule (`transient.py:141–142`).
- **W2 `steady` (α = `ALPHA`): the contraction identity, `MAT-5` step 3's gate set on this fixture.**
  - **(K)** \|ρ₄/μ_pred − 1\| ≤ `K_BAND`;
  - (K-conv) ≤ `KCONV_BAND`;
  - (M0) \|μ_pred\| ≥ `M0_FLOOR`;
  - (F) ‖δ₈‖/‖θ₈‖ ≤ `S_BAND`;
  - (E\*) ≤ `C_BAND`, (B′) ≤ `B_BAND`, (C′) ≤ `C_BAND`, (D′) ≤ `D_BAND` at iterate 8.
  - Every factor is **predicted**: `MAT-5` step 3 measured these on the ragged degree-1 twin, which is not this fixture. Its readings were ρ₄ +0.099100, μ_pred +0.099081, \|ρ₄/μ_pred − 1\| 1.83e-4, (K-conv) 2.6e-4, (F) 8.5e-8 (`20261005T094412Z_MAT-5-step3-r2-W2-coupled.log:1887, :1899–1901`).
  - Predicted here: μ_pred ∈ [+0.08, +0.12], (K) ≲ 1e-3, (F) ≲ 1e-6.
- **W3 `transient` (α = `ALPHA`):** (R1′) on q₀, and (E) on every run. **Nothing else is asserted. The two-way differences are the record.**

**Printed (predicted, never asserted).**
- **W1:**
  - one-way rows 2 / 4 against W3's h0025-d1x: ⟨θ⟩ 0.318088 / 0.614292 / 0.892065 / 1.153258 K, and θ₄₀ D 0.701208 K (`…W3-h0025.log:2038`, `:2041`). Predicted ≤ 1e-6.
  - P₁ from `_solve_driven` against W3's d1x sweep, 5.600076563536e-08 W (`:2030`). Predicted ≤ 1e-9; that is the held-factorisation premise, unmeasured on this mesh.
  - (Sup): the four schedules with every q_j = q₀, summed, against the one-way θ₄₀ (relative L2). Predicted ≤ 1e-12.
  - EM s/solve.
- **W2:**
  - ρ₁…ρ₇, μ₁…μ₅, μ_E;
  - θ̄₈/θ̄₁ − 1, predicted ≈ +10 to +12 % (`m5` read +11.15 %, `:1906`);
  - P_diss,8/P₁; θ_D at iterates 1 and 8;
  - EM s/solve and RSS.
- **W3, the price of the operator's two-way run:**
  - ⟨θ⟩ at the four report times and θ_D(t\*): one-way, two-way with 4 segments and two-way with 8 segments, plus the relative differences;
  - P_diss,j/P₁ per segment;
  - the two-way energy total Σ_j G·P_j·Δt_j against ρc∫θ(t\*) + the surface integral;
  - item 410's degree-2 one-way ⟨θ⟩(t\*) and θ_D beside.
  - Predicted two-way / one-way − 1 at t\*: ⟨θ⟩ +0.5 % to +1.5 % (4 segments), with 8 segments ≥ 4 segments (smaller lag); θ_D of the same sign.
  - Basis: the steady coupling μ ≈ 0.1 at θ̄ ≈ 5.4 K gives δP/P ≈ 0.018·θ̄. The segment-start ⟨θ⟩ values 0 / 0.32 / 0.61 / 0.89 K average ≈ 0.46 K, which gives ≈ +0.8 %.

**Negative controls (rule (e); all printed, none asserted).** No fixture-matched prior measurement exists; item 396's Q5 ruling is the precedent.
- **(N1) frozen σ₀.** The miss of (K) is 1.000 by construction (`m5:423`), backed by (R0c).
- **(N2) frozen \|E\|.** \|μ_E/μ_pred − 1\|, predicted 0.25–0.45 (the twin read 0.350, `…W2-coupled.log:1903`). It shows (K) separates a re-solved E from a frozen one. The factor is predicted.
- **(T-α0), W1:** the 4-segment stagger at α = 0 (3 more solves) against the one-way θ₄₀. Predicted ≤ 1e-9.
  - Against W3's α = `ALPHA` difference of ≈ 1e-2, the predicted separation is ≈ 1e7.
  - It shows the printed two-way difference is σ(T), not a stagger or superposition artefact.

**Windows.** Foreground. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Bash `timeout: 660000`; `-s`; durable capture verbatim. **Two heavy windows never run at once; `pgrep -c python3` must read 0 before each.**
- **First:** pgrep, then the cache sweep, label `ANS-9-step2-cache-sweep`.
- **W0 (smoke, `-n 1`):**

  `scripts/testing/run_and_log.sh ANS-9-step2-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_two_way_sigma.py -s -v --tb=short'"`
- **W1 `alpha0`** (heavy, `-n 8`, 6 degree-1 solves). Predicted harness 120–190 s: build 40.87 s (`…W3-h0025.log:1953`) plus 6 × 8–20 s per solve (`m5` W1 read 19.3–19.7 s per solve at `-n 2` on 199 920 cells, `20261005T094115Z_MAT-5-step3-r2-W1-alpha0.log:1883`), plus ≈ 10 CN runs at ≈ 2 s.

  `scripts/testing/run_and_log.sh ANS-9-step2-W1-alpha0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS9_STEP2=alpha0 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_ans9_two_way_sigma.py -s -v --tb=short > /workspace/logs/ans9-step2-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/ans9-step2-w1-raw.log; cat /workspace/logs/ans9-step2-w1-raw.log; exit \$rc'"`
- **Start rule for W2 and W3:** only if W1's **max** EM s/solve is ≤ 35 s. W2 is ≈ 41 + 13·s + ≈ 50 s and passes 590 s at s ≈ 38 s. Otherwise record "W2 / W3 not attempted (EM … s/solve)", which is a price record, not a red.
- **W2 `steady`** (13 solves + 13 steady + 5 thermal-only). Predicted 200–360 s. The command is W1's with knob `steady`, label `ANS-9-step2-W2-steady` and raw file `/workspace/logs/ans9-step2-w2-raw.log` (all three occurrences).
- **W3 `transient`** (1 + 3 + 7 = 11 solves, ≈ 16 CN runs). Predicted 180–320 s. Knob `transient`, label `ANS-9-step2-W3-transient`, raw file `/workspace/logs/ans9-step2-w3-raw.log`.
- **Time rule:** if W3 has not started by minute 45, commit with "⬜ W3 remaining (module at `<sha>`)" (item 396's Q3 form).

**Traps already paid for.**
- **Never reuse a degree-2 case dict.** `_solve_driven` takes `sweep["degree"]` (`tests/validation/test_birdcage_b1_plus_map.py:164–165`), and `build_four_port_sweep` returns the degree it was given (`tests/validation/test_port_birdcage_four_port.py:528`). Build with `degree=1`. One accidental degree-2 solve costs ≈ 245 s.
- **Never call `m5._Fixture(alpha)`.** Its `__init__` builds the ragged 199 920-cell fixture (`m5:122`), which is the wrong mesh and 40 s lost.
- **G is fixed from the σ₀ solve and never renormalised** (fixed source amplitude, `m5`). The two-way deposits more than 1 W, and that excess is the effect.
- **No `src/` restart.** Superposition replaces a field `t_init`; never add one in-slot. χ is evaluated at step midpoints (`transient.py:136–140`), and boundaries sit on full steps, never on a midpoint.
- **Complex build.** Masks and schedules are Python floats, with no UFL comparison. `law.sigma` gets `np.real` arrays (`m5:166–167`).
- **Reductions.** Every one is real and allreduced on every rank, with no collective under `rank == 0`. Cell counts come from `size_global`. θ_D is read with the default owner on the sub, and the points are identical on every rank.
- **(R0c) at `-n 8` is not bitwise.** Item 237 measured 2.747e-12 against an older 1e-12; `L0_BAND` 1e-9 is the ruled band (`m5:91–95`).
- **Imports.** Import no fixture or `test_*` name. Use `-s`, `-k 30` and node-free paths; keep `; exit $rc`.
- **After a kill:** pgrep, then the sweep. Status 124 is a price record.
- **Rule (i):** no module edit after W0 starts. An edit re-runs W0 → W3.

**Scope.** The in-house half of (8d) on P24's mesh at degree 1 (re-posed): one drive, 10 MHz, linear σ(T) with α = `ALPHA` 0.02 /K as an input (not a measured gel coefficient), fixed source amplitude, a steady Picard fixed point, and an explicit stagger at the report times (lagged, with an 8-segment halving printed, not a converged transient coupling). **Not covered:**
- degree 2 (the `xl` form);
- ε(T);
- perfusion and tissue;
- any AED value or the operator's two-way Icepak run;
- an absolute temperature;
- h-convergence;
- C4.

**Negative result.** Report the measurement and stop. Never widen a band, and never raise N or M after a red.
- **(R), (R0a), (D0) or (R0c) red:** a fixture or hook wiring defect. Park on `attempt/ANS-9-step2-<ts>` and mark BLOCKED (rule (d)).
- **(K) red with (F), (K-conv) and (M0) green:** the contraction does not match its linearised prediction on this fixture. That is a finding. Park and mark BLOCKED with the unblock condition "review rules"; add a known-issues entry.
- **(K-conv), (F) or (M0) red:** report the μ_j sequence and stop.
- **(E) or (R1′) red:** a source-transfer defect. Park and mark BLOCKED.
- **W3's differences, whatever they read:** the record.
- No `PREDICTED RED` line.

**Status it can move.**
- **`ANS-9` row.** The state line gains "step 2 record (two-way σ(T), degree 1 on P24's mesh): (K) …, ⟨θ⟩(t\*) two-way / one-way − 1 = … (4 / 8 segments), θ_D …; `<log>`". No glyph: the Done-when concerns step 1's adjudication. Rule (j) is not engaged.
- **§10 P27** gets a DONE line, and (8d)'s in-house half is recorded on P24's fixture for the weekly.
- **`MAT-5` row:** unchanged (✅ on steps 1–2).
- **If the review keeps the letter** (degree 2), this item is not queueable as heavy. It goes to the weekly as an `xl` candidate priced above, and §9 carries nothing.

**Spawn line:** "Use the implementer agent to execute `ANS-9` §9 item 425 per `docs/planning/chunks/ANS-9.md` (read the brief and the newest Review ruling beneath it first; skip unless item 410 is DONE). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`; two heavy windows never run at once; `pgrep -c python3` reads 0 before each heavy window. Build route B at degree 1; never solve at degree 2 in this item. No `src/` edit; no edit to any existing module; never import a fixture or `test_*` name. Never open docs/private/, aed_results/, aed/ or COMPARISON_private.md. You are the executor: do not spawn agents. Report a red and stop; never widen a band."

**Same-commit duties.** Mark §9 item 425 DONE, BLOCKED or "⬜ W3 remaining". Commit the module and every harness log together. Add the row state clause, the §10 P27 DONE line, an outcome paragraph here (the (K) table, the W3 one-way / two-way table, EM s/solve, RSS) and the journal entry.

### Review ruling 2026-10-13 03:00 (daily) on the item 410 / 425 drafts

1. **Item 410 is queued as drafted.** Row 6's energy balance is asserted at `E_BAND` as well, because the §7 Done-when names rows 1 / 3 / 6 on the primary column. The degree-1 column is carried as a labelled constant taken from the tracked `metrics.json` at `59809ba7` and is not re-solved. The guide's run line moves to `-n 8 -t 590`. **Binding:** the slot first lands a `COMPARISON.md` writer that keeps the hand-added "## Verdict" section byte for byte, the item-400 pattern, with its smoke unit test green before the solve window. `SPEC_THERMAL` stays unchanged, and all compute stays inside `main()`, because `MAT-19` and the step 1b / 1c modules import the script by path. The two degree-1-only checks become prints on the new fixture: re-run equals gated run, and the phantom volume against the 58 866-cell record.
2. **Item 425 (P27) is not queued; it goes back to the 10-14 weekly.** As posed — two-way σ(T) on P24's **degree-2** column — it needs ≈ 13 + 3 field re-solves with no reusable factorisation. That is ≈ 230–260 s each at `-n 8` on the h0025 mesh (`20261012T100436Z_ANS-9-step1c-W3-h0025.log:2054`), ≈ 65 min in all, more than one heavy window. The drafted degree-1 re-pose would run on a different column than the step names, which changes its target, and this daily may not do that (rule (6)). The weekly chooses between an `xl` cost probe at degree 2 and the degree-1 re-pose drafted above. The draft stays in this file as the brief either choice can use. The step text also names "`MAT-14`'s route"; `MAT-14` has no σ(T). The coil-driven σ(T) route is `MAT-5` step 3 (`tests/validation/test_mat5_coil_driven_sigma_feedback.py:196–205`), which the draft reads. The weekly confirms.

### §9 item 475 — 2026-10-15 (§10 chain P31; lane 1, no solve)

*Drafted by a general-purpose agent for the 2026-10-15 03:00 daily review; the review reviews it.*

**Licence.**
- §10 P31 (PROJECT_PLAN.md:4781–4784), weekly-enumerated, §9 rule (4): "The `mat:`-family guides and the `ANS-9` case `.md` that cite 'DISAGREE at D' re-read as 'AGREE within 2.5 % on the degree-2 column' (the `COMPARISON.md` verdict, the §7 row and §2.1 were written tonight); `docrefs` 0 dead. Standard. Serial on nothing."
- The 2026-10-14 weekly ruling "`ANS-9` — every row AGREE on the production convention; `ANS-9` ✅" (PROJECT_PLAN.md:4578–4594).
- No band, tolerance, Done-when or target is set or moved. No §7 / §2 / `COMPARISON.md` edit (the weekly wrote them).

**What the sweep found (drafter, 2026-10-15; the scope is smaller than the step assumed).** Searched: `grep -rn "DISAGREE"` and `grep -rni "disagree\|probe D"` over `examples/`, `docs/guide/`, `docs/*.md`, `docs/status/START_HERE.md`, `README.md`, `.md` and `.py`, private files excluded.
- **`examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz/SPEC.md:49–54`** — the status checklist's ticked "Adjudication … done 2026-10-10 … DISAGREE at probe D (ours 11 % below AED, band 10 %) … step 1b (PROJECT_PLAN §10 P21) diagnoses D." **The only live citation; this item edits it.** (`:46–48`, the interactive first pass, is history and stays.)
- `…/ans9_coil_driven_gel_heating_10MHz/COMPARISON.md:82`, `:90–95` — inside `## Verdict (weekly review 2026-10-10; … — superseded above on the degree-1 column)` (`:64`): a labelled record. **Not edited.** `:55` is the 10-14 verdict's own "the 2026-10-10 DISAGREE at probe D was discretisation …" — correct as written.
- **The `mat:` family cites no DISAGREE and no probe-D gap.** `examples/materials/*.md` hits for "disagree" are only `01_dodd_deeds_coil_loading.md:172`, `:223` (generic, unrelated); the `ANS-9` mentions (`09_coil_driven_heating_phantom.md:33`, `10_finite_cylinder_robin_series.md:30`, `14_coil_driven_phantom_in_container_wall.md:46`, `15_duty_cycled_coil_heating_phantom.md:74`, `16_probe_readouts_in_coil_phantom.md:13`, `:65`, `:125`) are "not claimed: `ANS-9`" or SPEC-probe-centre references — no verdict wording. **No `mat:` guide edit.**
- **`…/ans9_coil_driven_gel_heating_10MHz/09_coil_driven_gel_heating_10MHz.md`** (the case guide) cites no DISAGREE; `:14` says the degree-1 column was "adjudicated 2026-10-10" (true) and `:78` "No AED value or verdict" (a deliberate scope line). **Not edited** — adding a verdict pointer would contradict `:78`; the review may rule otherwise.
- **`.py` docstrings:** `09_coil_driven_gel_heating_10MHz.py:244` (the SPEC readout's generic "A disagreement on rows 2 / 3 …") and `:497` (an (A0) error string) are unrelated. **No `.py` edit, so rule (i) is not engaged and no example re-runs.**
- `docs/guide/`, `docs/status/START_HERE.md`, `README.md`: no hit. `docs/testing/known-issues.md:816` is already ✅ RETIRED 2026-10-14.
- **Outside this item's reach, flagged to the review:** `CLAUDE.md` (project instructions) still says "probe D's gap was read as our discretisation and is being re-adjudicated (🟡)". Agents do not edit `CLAUDE.md`; the operator or an interactive session does.

**The edit (one file, two touches; the slot transcribes).** In `SPEC.md`:
1. `:54`, end of the ticked 10-10 line: "… step 1b (PROJECT_PLAN §10 P21) diagnoses D." → "… step 1b (PROJECT_PLAN §10 P21) diagnoses D. *Superseded 2026-10-14 (next item).*"
2. A new ticked line directly after `:54`:
   "- [x] Re-adjudication on the degree-2 column (weekly, **2026-10-14**; numbers in the gitignored `docs/private/ans9-readjudication-2026-10-14.md`; no new AED run) — on step 1d's primary column (N1curl degree 2, exact rectangles, κ off, h_p 0.0025 m, the production convention) every row **AGREES**: probe D within 2.5 % and probes A–C within 2.3 % (band 10 %), B > D > C held in both codes; ⟨θ⟩ within 0.8 % at every report time, the initial slope within 0.5 %, the surface heat within 1.3 % on the Simpson figure (bands 5 %); rows 5 / 7 records. The 2026-10-10 DISAGREE at D was our discretisation — 13 % phantom resolution and 8 % element order (steps 1b / 1c). `ANS-9` ✅ (step 1); known-issues 2026-10-10 retired; verdict table in `COMPARISON.md`."

**The anchor — every figure in the new text traced (relative levels only; no AED value).**

| Figure | Source (tracked) |
|---|---|
| probe D within 2.5 % (band 10 %) | PROJECT_PLAN.md:4579–4580; `COMPARISON.md:49` |
| probes A–C within 2.3 %, B > D > C | PROJECT_PLAN.md:4581–4582 |
| ⟨θ⟩ within 0.8 % at every report time | PROJECT_PLAN.md:4582; `COMPARISON.md:47` |
| initial slope within 0.5 % | PROJECT_PLAN.md:4582–4583; `COMPARISON.md:48` |
| surface heat within 1.3 % on the Simpson figure (bands 5 %) | PROJECT_PLAN.md:4583; `COMPARISON.md:52` |
| rows 5 / 7 records | PROJECT_PLAN.md:4584; `COMPARISON.md:51`, `:53` |
| 13 % phantom resolution, 8 % element order | PROJECT_PLAN.md:4586–4587; `COMPARISON.md:55–57` |
| h_p 0.0025 m, degree 2, exact, κ off | `COMPARISON.md:3`, `:38–40` |
| private file name | `COMPARISON.md:42` |

Plus docrefs `dead=0`.

**Negative control: none, and why.** No code, test or generated file changes; the text is a transcription of figures already public. The private-leak hook cannot see a 2–3-digit relative level by design (`scripts/testing/check_private_leak.py:45–47`) — the written rule (relative levels only, never an AED value) is the control here. The before / after sweep (W0 / W2) is the record that the one live citation moved and nothing else did.

**Executor, tier, cost.** Implementer, or the slot directly. Standard tier; no MPI; no dolfinx. Lane 1, class median 3 measured min; predicted ≈ 15 min. Depends on nothing.

**Windows.** Foreground, bare, repo-relative, Bash `timeout: 660000`.
1. **W0 (before the edit), the sweep — private files excluded by flag:**
   `scripts/testing/run_and_log.sh ANS-9-item475-sweep-pre "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && grep -rn --include=*.md --exclude=COMPARISON_private.md --exclude-dir=aed_results --exclude-dir=aed DISAGREE examples/materials examples/ansys_benchmarks/ans9_coil_driven_gel_heating_10MHz'"`
   Predicted hits: `SPEC.md:52`, `COMPARISON.md:55`, `COMPARISON.md:82` — nothing in `examples/materials`.
2. **Edit** `SPEC.md` (Edit tool).
3. **W2 (after):** the same command, label `ANS-9-item475-sweep-post`. Predicted **four** hits: W0's three plus the new SPEC line (it says "The 2026-10-10 DISAGREE at D was our discretisation") — `SPEC.md:52`, `SPEC.md:55`, `COMPARISON.md:55`, `COMPARISON.md:82`. Any other hit is a miss of this sweep — report it.
4. **Docrefs (gate `dead=0`, exit ≠ 1):**
   `scripts/testing/run_and_log.sh ANS-9-item475-docrefs "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && timeout -k 30 60 python3 scripts/testing/check_example_doc_references.py'"`
5. **Commit** `SPEC.md` and the three logs; the pre-commit leak hook must print `clean:` (`SKIPPED` is not a pass).

**Traps already paid for.**
- **The sweep must exclude the private files by flag.** A bare `grep -rn DISAGREE` over the case directory prints `COMPARISON_private.md` lines into a **tracked** harness log — a leak through the log. Keep `--exclude=COMPARISON_private.md --exclude-dir=aed_results --exclude-dir=aed`; read the W0 log before anything else and, if any private path appears, do not commit it: move it to the gitignored `logs/` and stop.
- **Never quote an AED value**; only the relative levels above. Never name a Palace figure in this file (none is involved).
- **`COMPARISON.md` is the weekly's and its verdict tail is preserved byte for byte by the writer** (`09_…py:200–214`); this item does not touch it.
- No `.py` edit (rule (i) would require the `ans:9` window, ≈ 6 min heavy at `-n 8`).
- Formatting: `≈`, never a tilde (a single `~text~` strikes through in GFM). Never delete a log.

**Scope.** One status line in one SPEC. **Not claimed:** any new adjudication; a `mat:` guide change (none cites the old verdict); step 2 (P27, not posed); `CLAUDE.md` (flagged).

**Negative result.** If W0 shows a live DISAGREE citation the drafter missed, edit it the same way only if it is a guide `.md`; a `.py` hit is flagged to the review, not edited (rule (i)). If docrefs reports `dead > 0` on a line this item did not touch, record it and stop (not this item's to fix). No `PREDICTED RED:` line.

**Status it can move.** No glyph (`ANS-9` ✅, PROJECT_PLAN.md:3532; rule (j) read, no mismatch). §10 P31 → DONE (the review marks it). The `CLAUDE.md` sentence is the operator's to update.

**Spawn line:** "Use the implementer agent to execute `ANS-9` §9 item 475 per `docs/planning/chunks/ANS-9.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped; foreground, `timeout: 660000`. Edit only `SPEC.md` as written; the sweep window keeps its private-file exclusions; never quote an AED value; never open `aed_results/`, `COMPARISON_private.md` or `docs/private/`."

### Review ruling 2026-10-15 03:00 (item 475)

Accepted at its measured scope: the only live citation of the 10-10 "DISAGREE at probe D" is the ans9 `SPEC.md` status line; no `mat:` guide or example script cites it, so no rule-(i) re-run arises. The sweep window keeps its private-file exclusions so no `COMPARISON_private.md` line reaches a tracked log. `CLAUDE.md`'s stale "re-adjudicated (🟡)" clause is the operator's file and stays on the dashboard.

### Outcome — §9 item 475, 2026-10-15 21:30 slot (take-next)

**DONE.** `SPEC.md` edited as written (two touches): the 10-10 line ends "*Superseded 2026-10-14 (next item).*", and the ticked re-adjudication line follows it. Every relative level traced: PROJECT_PLAN.md:4609–4619 (the brief's 4578–4594 had shifted ≈ 30 lines) and `COMPARISON.md:47–57`; "A–C within 2.3 %" is from the plan ruling only (`COMPARISON.md:49` states the 2.5 % envelope over A–D). No AED value written. W0 3 hits (`SPEC.md:52`, `COMPARISON.md:55`, `:82`); W2 4 hits, the new one at `SPEC.md:65` (the brief predicted `:55`; the line wraps over twelve lines) — no other hit, no private path in either log. Docrefs `dead=0` (exit 2: 120 stale `paraview_output/` artifacts, report severity, untouched). Logs `20261015T214915Z_ANS-9-item475-sweep-pre.log`, `20261015T214923Z_ANS-9-item475-sweep-post.log`, `20261015T214925Z_ANS-9-item475-docrefs.log`.
