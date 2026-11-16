# MAT-15 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record. Steps 1 and 2(a) (§9 items 130, 137) are archived verbatim in `docs/planning/plan-archive.md` § "§9 digest and restock of the 2026-10-01 03:00 review".

### §9 item 145 — 2026-10-01 (chain P2b; queued by the 10:30 review)

`MAT-15` step 2(b) — the FEM on the phantom's own mesh against the finite Robin cylinder series (implementer; standard by ceiling, `-n 2`, complex; ≈ 10 measured / ≈ 55 predicted min; depends on nothing — item 137 landed `0f5bb02c`).

- **Anchor.** No EM solve. The mesh is `_build(True, phantom_resolution=h_p)` (`tests/mesh/test_birdcage_port_sheets.py:80`), then `_tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)` and `_robin_tags(sub)` — step 1's thermal domain (`test_mat15_coil_driven_heating_phantom.py:106–121,211–212`). Imported from the step-1 module: `K_GEL`, `H_GEL`, `T_INF`, `A`, `T_R`, `T_STAR` (1254 s), `N_MAIN`, `_tag_submesh`; from `MAT-12`: `RHO_C`, `E_BAND`, `INF_BAND`; from `MAT-13`: `INTERIOR_BAND`; from 2(a)'s test: `Q`, `L`, `KW`. Comparand: P2a's series in numpy, interpolated into P2 in chunks of 20 000 points (`MAT-9`'s `_exact_theta` pattern), r = √(x² + y²).
- **Precondition (asserted).** 2(a)'s (A, L, K, H, RHO_C) equal step 1's (A, `PHANTOM_HEIGHT`/2, K_GEL, H_GEL, RHO_C). The submesh is centred and aligned, from rank-reduced bounds: z ∈ [−L, L] to 1e-9·L, max r ≤ a·(1 + 1e-9), centroid within 1e-9·a of the origin.
- **Anchors (asserted).**
  - (S1) Steady P1 at h_p ∈ {0.0075, 0.005, 0.00375, 0.0025} — `GEO-27`'s rungs: 2 746 / 8 497 / 18 570 / 58 866 phantom cells (`20260908T033218Z_GEO-27.log.gz:1804,3573,5346,7142`). Whole-cylinder relative L2 on θ; the LSQ slope of log e against log h_eff, h_eff = (V_sub/N)^{1/3}, ≥ 1.7 (`MAT-9` step 1's bar; measured 2.039 on the ball, `20260929T034247Z_MAT-9-step1.log`). At 0.0025, interior relative L2 (r < 0.9a, \|z\| < 0.9L) ≤ `INTERIOR_BAND`.
  - (S2) Gauss balance on every rung: \|∮hθ dS − qV_sub\|/(qV_sub) ≤ 1e-9.
  - (T) CN with the Rannacher start, N = 40 to t* at 0.0025: interior ≤ `INTERIOR_BAND`.
  - (E) every run's `max_balance_residual` ≤ `E_BAND`. (∞) backward Euler 20 × 5 t_R equals the steady solve to ≤ `INF_BAND`.
- **Printed (predicted).** Finest steady interior ≈ 1–3e-3; θ_h at the centre beside P2a's 2.983150 K; (T) ≈ 2e-3; pairwise rates. V_sub against πa²·2L (≈ −8.7e-4 from the record 2.259978941112e-04 — the wall faceting). max \|Im θ\| (predicted 0).
- **Negative controls (rule (e): predicted; asserted only as misses of `INTERIOR_BAND`).** (c1) the comparand with h → 2h: steady miss predicted 0.3–0.5. (c2) Dirichlet radial roots, `radial_roots=scipy.special.jn_zeros(0, 200)`, in the transient comparand at t*: miss predicted 0.05–0.5.
- **Window.** W1 (≈ 4 min; rung builds 28.1 + 28.7 + 31.8 + 39.9 s, `GEO-27.log.gz:1807,3576,5349,7145`; CN at 0.0025 is 2.5 s, `20261001T123309Z_MAT-15-step1-W1.log:1929`): `scripts/testing/run_and_log.sh MAT-15-step2b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 480 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat15_finite_cylinder_fem.py -s -v --tb=short'"`
- **Traps.** P2a's functions build (points × 400) and (points × 200 × 200) arrays: always chunk. Never evaluate the transient near t = 0. In the complex build every reduction takes `np.real` (`MAT-13`'s `:134–149`); `MAT-9`'s `float(np.sqrt(e/n))` is not complex-safe. `cell_tags` and `assemble_scalar` are rank-local. Rung cell counts other than 58 866 are printed, not asserted. A killed compile leaves 0-byte `.c` stubs: sweep them. FD ratios on summed series pre-register δ ≳ 1e-2 of the length scale (items 136 / 137). Foreground, `timeout: 660000`; no edit after the last window.
- **Scope.** Uniform source, Robin walls, P1, four rungs. Not: EM, perfusion, degree 2, `ANS-9`.
- **Status it can move.** §7 `MAT-15` 🟡 → ✅ (the row's 2(b) Done-when, pinned by the 10:30 review).
- **Negative result.** Rate < 1.7 with the finest interior green: print the pairwise rates; if only the 0.0075 pair is low (≈ 4 cells across a), do not drop the rung in-slot (rules (f), (h)) — record, park on `attempt/*`, BLOCKED; the review rules. Interior > 1 % at 0.0025: report beside the faceting deficit and stop. (S2) red: a code bug.
- **Outcome 2026-10-01 (15:00 slot) — 🚫 BLOCKED, parked `attempt/MAT-15-step2b-20261001T203800Z` (`640f2594`).** W1 `20261001T202818Z_MAT-15-step2b-W1.log` (`-n 2`, complex) was killed by the container `timeout -k 30 480` (exit 124, 481 s) after the (c2) line (`:7167`), inside the (∞) backward-Euler tail; no orphaned ranks; no assertion on the fixture reached. Two blockers: (1) **the pre-registered centroid bar is unmeetable on this mesh** — the submesh volume centroid is (2.09e-8, 6.33e-9, −5.91e-8) m, \|c\| = 6.3e-8 m = 2.1e-6·a against 1e-9·a (`:7163`; z ∈ [−L, L] and max r = a hold to 1e-12) — the unstructured lateral faceting is not mirror-symmetric, so the bar needs a review ruling (measurement recorded, not moved in-slot); (2) the four rung builds take 185 s and the finest-rung tail (two P2 interpolations of the 200 × 200 transient series in 1 000-point chunks, CN, BE) did not finish in the remaining ≈ 290 s — the next window needs the tail split or priced. Measured, all printed (`:1852,3611,5375,7162–7167`): cells 2 746 / 8 497 / 18 570 / 58 866; h_eff 4.340e-3 / 2.982e-3 / 2.299e-3 / 1.566e-3 m; whole relative L2 6.744e-3 / 2.943e-3 / 1.686e-3 / 7.520e-4 — LSQ slope ≈ 2.15 (computed from the printed values; the fixture's own slope line was not reached), pairwise 2.21 / 2.14 / 2.10; (S2) Gauss ≤ 1.3e-14 every rung; finest interior 7.37e-4 (predicted 1–3e-3); θ_h(0) 2.980951 K vs 2.983150 K (−7.4e-4); (T) CN N = 40 interior 3.49e-4 (predicted ≈ 2e-3), balance 7.2e-14; V_sub deficit −8.6995e-4 (predicted ≈ −8.7e-4); controls (c1) 2h 0.739 (predicted 0.3–0.5 — outside its predicted range, high), (c2) Dirichlet roots 0.779 (predicted 0.05–0.5 — outside, high), both miss 1 %. Not measured: (∞), (E) on the BE run, max \|Im θ\|.

**Ruled 2026-10-01 18:00 review (§9 rule (f); the 10-03 weekly ratifies or reverts).** (1) **Centroid bar.** The precondition's 1e-9·a centroid bar is the brief's own idealisation (a perfectly centred cylinder), and an unstructured lateral faceting cannot meet it: the measured \|c\| is 6.3e-8 m = 2.1e-6·a (`…202818Z_MAT-15-step2b-W1.log:7163`). It is re-registered by measurement at **\|c\| ≤ 1e-5·a**, 4.8× the measured value. A centroid offset δ moves the comparand by O(δ/a) ≈ 2e-6 relative on θ, four orders under `INTERIOR_BAND` 1e-2, so the bar still excludes a displaced or mis-scaled submesh (a 1 mm shift reads 3.3e-2·a). The z-extent and max-r bars stay at 1e-9 (they hold to 1e-12). (2) **Window split.** The four rung builds cost 36.3 / 39.2 / 37.7 / 71.5 s (`:1852,3611,5375,7162`). The tail then printed (S1), (c1), (T) (solve 2.6 s) and (c2), and was killed inside (∞) with ≈ 260 s spent past the finest rung and no per-phase clock. A backward-Euler run of 20 steps on 58 866 P1 cells is seconds, so the item splits the steady and transient halves into two windows and prints a wall clock per phase, so that a stall is located rather than re-bought. (3) **Controls.** (c1) 0.739 and (c2) 0.779 fell outside their predicted ranges on the high side. Both are printed predictions per rule (e), asserted only as misses of 1 %, and both missed: no ruling needed. The re-run predicts the measured values. (4) **Records.** The W1 log and its `test-results.md` row exist only on the parked branch. The landing takes the log by path, a new file, and appends the row **by hand** — never a path checkout of `test-results.md`.

### §9 item 160 — 2026-10-01 (chain P2b, item 145 re-run on the ruled bar; queued by the 18:00 review)

`MAT-15` step 2(b): the FEM on the phantom's own mesh against the finite Robin cylinder series, run as two windows (implementer; heavy, `-n 2`, complex; ≈ 15 measured / ≈ 40 predicted min; depends on nothing). Family `MAT-15` step 2: 2(a) registered its records (item 137), so this is **attempt 2 of 4** since then.

- **Code.** Take `tests/validation/test_mat15_finite_cylinder_fem.py` **by path** from `attempt/MAT-15-step2b-20261001T203800Z`, together with its W1 log `docs/testing/logs/20261001T202818Z_MAT-15-step2b-W1.log` (a new file). Append that log's `test-results.md` row by hand from the branch's diff. Then make the following edits. (1) The centroid assert becomes `<= 1e-5 * A`, with a comment citing this ruling and `:7163`. (2) Split the module-scoped `run` fixture into two module-scoped fixtures. `steady` builds the four rungs and computes (S1), (S2), the precondition and (c1). `transient` builds the 0.0025 rung alone and computes (T), (E), (∞) and (c2). A pytest node only instantiates the fixtures it uses. Wrap each phase in `time.perf_counter()` and print the timing with `flush=True`: build, steady solve, the series interpolation, CN, the BE (∞) run, and each control. (3) Keep the parked chunk sizes (`STEADY_CHUNK` 20 000, `TRANSIENT_CHUNK` 1 000, `:67–68`): the (points × 200 × 200) complex einsum at 20 000 points would be ≈ 12.8 GB per rank. The per-phase clocks will say whether the transient interpolation is the cost. Assert nothing new, and move no band other than the centroid bar.
- **Anchors (asserted), unchanged from item 145 apart from the bar.** Precondition: 2(a)'s constants equal step 1's; z ∈ [−L, L] to 1e-9·L; max r ≤ a(1 + 1e-9); centroid ≤ 1e-5·a. (S1) LSQ rate ≥ 1.7 over the four rungs (predicted 2.15 from W1's printed values) and interior ≤ `INTERIOR_BAND` at 0.0025 (predicted 7.37e-4). (S2) Gauss ≤ 1e-9 on every rung (predicted ≤ 1.3e-14). (T) CN N = 40 interior ≤ `INTERIOR_BAND` (predicted 3.49e-4). (E) every run's `max_balance_residual` ≤ `E_BAND`. (∞) BE 20 × 5 t_R = steady ≤ `INF_BAND`. Cell count 58 866 at 0.0025.
- **Negative controls (rule (e)).** (c1) comparand at 2h: miss **predicted 0.739** (W1 `:7165`). (c2) Dirichlet radial roots: miss **predicted 0.779** (`:7167`). Each is asserted only above `INTERIOR_BAND`.
- **Printed.** Pairwise rates (W1: 2.21 / 2.14 / 2.10), θ_h(0) beside 2.983150 K, the V_sub deficit, max \|Im θ\| (predicted 0), and every phase clock.
- **Windows** (complex; foreground; Bash `timeout: 660000`).
  - W1, steady (≈ 200–260 s; builds 185 s measured): `scripts/testing/run_and_log.sh MAT-15-step2b-W2-steady "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 480 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat15_finite_cylinder_fem.py::test_precondition_2a_constants_are_step1s tests/validation/test_mat15_finite_cylinder_fem.py::test_precondition_submesh_centred_and_aligned tests/validation/test_mat15_finite_cylinder_fem.py::test_finest_rung_cell_count tests/validation/test_mat15_finite_cylinder_fem.py::test_s1_steady_rate tests/validation/test_mat15_finite_cylinder_fem.py::test_s1_steady_interior_finest tests/validation/test_mat15_finite_cylinder_fem.py::test_s2_gauss_balance_every_rung tests/validation/test_mat15_finite_cylinder_fem.py::test_control_c1_doubled_h_misses -s -v --tb=short'"`
  - W2, transient (predicted 120–400 s, unmeasured past (c2); heavy ceiling): `scripts/testing/run_and_log.sh MAT-15-step2b-W3-transient "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 900 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat15_finite_cylinder_fem.py::test_t_transient_interior tests/validation/test_mat15_finite_cylinder_fem.py::test_e_every_run_balances tests/validation/test_mat15_finite_cylinder_fem.py::test_inf_long_backward_euler_is_steady tests/validation/test_mat15_finite_cylinder_fem.py::test_control_c2_dirichlet_roots_miss -s -v --tb=short'"`
  - The test names are the parked module's (`:245–305`). If the split renames them, use the new node IDs and say so in the journal. Never use `-k` selection inside the quoted command.
- **Traps.** **Run this first: `find /root/.cache/fenics -name '*.c' -size 0` in the container.** W1 was killed mid-run and may have left a 0-byte stub, which is a live lock. Delete stubs only, never the whole cache. Then `pgrep -c python3` to check for orphans. Item 145's traps also stand. P2a's functions build (points × 400) and (points × 200 × 200) arrays, so always chunk. Never evaluate the transient near t = 0. Every reduction takes `np.real`. `cell_tags` and `assemble_scalar` are rank-local. If (∞) prints no phase clock within 120 s of (c2), suspect a rank-divergent collective (an `if` on a rank-local value ahead of an allreduce) before suspecting cost. Status 124 is a finding, not a retry. No edit after a file's last window (rule (i)).
- **Scope.** Uniform source, Robin walls, P1, four rungs. Not covered: EM, perfusion, degree 2, `ANS-9`.
- **Status it can move.** §7 `MAT-15` 🟡 → ✅ on the row's 2(b) Done-when (pinned by the 10:30 review; the centroid bar is a precondition, not a Done-when term). `attempt/MAT-15-step2b-20261001T203800Z` is deleted after the landing.
- **Negative result.** W2 killed at 900 s: record which phase clock last printed, park, BLOCKED — the review prices or re-poses (∞). That would be attempt 3 of 4. A rate below 1.7 or an interior red: item 145's rule — record, park, BLOCKED, the review rules. A red centroid at the new bar means the mesh has changed under the item: stop.

**Outcome 2026-10-01 (22:30 CDT slot) — 🚫 BLOCKED, parked `attempt/MAT-15-step2b-20261002T034707Z` (`7120eeed`).** Slot deviation, disclosed: W2 ran at container `timeout -k 30 590`, not 900 (the foreground Bash window is 660 s; see the slot note below). Trap step: one 0-byte `.c` stub deleted, no orphans. **W1 steady green** — `20261002T033219Z_MAT-15-step2b-W2-steady.log` (`-n 2`, complex, 244 s, 19 passed): cells 2 746 / 8 497 / 18 570 / 58 866; LSQ slope 2.1511 (bar 1.7), pairwise 2.2106 / 2.1415 / 2.1015; finest interior 7.374e-4; Gauss ≤ 1.3e-14; centroid \|c\|/a 2.101e-6 (bar 1e-5), z and r exact; θ_h(0) 2.980951 vs 2.983150 K; V_sub −8.6995e-4; (c1) 0.7394 (predicted 0.739); max \|Im θ\| 0; phase clocks: builds 26.2 / 27.8 / 29.9 / 40.0 s, steady solves ≤ 0.6 s, steady series interpolation 0.9 / 3.4 / 5.8 / 27.8 s, (c1) 29.1 s. **W2 transient, first run** `20261002T033630Z_MAT-15-step2b-W3-transient.log` (286 s, 4 errors): FFCx JIT lock on `libffcx_forms_9f4f0568….c`. That was a **44 345-byte** `.c` with no `.o` / `.so`, written at 20:36:19 UTC, the instant item 145's W1 was killed. The brief's 0-byte sweep cannot see a killed compile that already wrote its source. Only that file was deleted and W2 re-run once (no code change). **W2 re-run** `20261002T034149Z_MAT-15-step2b-W3-transient.log` (273 s, 15 passed / 1 failed): (T) CN N = 40 interior 3.4945e-4 (predicted 3.49e-4), whole 4.16e-4; (c2) 0.7795 (predicted 0.779); (E) worst 7.234e-14 (CN), BE 1.3e-14; max \|Im θ\| 0; **(∞) 0.9999999999999808 > `INF_BAND` 1e-8 — red outside the predictions, so the slot stopped.** Phase clocks: build 39.4 s, steady 0.6, CN 2.3, transient series interpolation 95.1, (c2) 86.1, BE (∞) 1.2 s, fixture 225.9 s (item 145's overrun was the two 1 000-point-chunk interpolations, not (∞)). **Diagnosis (code reading, the arithmetic fits): a test bug, inherited verbatim from item 145's module.** `_rel(sub, a, b)` returns ‖a − b‖/‖b‖, but the (∞) call passes a = T_be − T_st and b = T_st − T∞. So it reads ‖(T_be − T_st) − θ‖/‖θ‖ = 1 − ⟨T_be − T_st, θ⟩/‖θ‖², which is 1 − 1.9e-14 when T_be = T_st. The 1.9e-14 deficit is the BE-vs-steady agreement. **Unblock:** one line, `out["inf"] = _rel(sub, res_inf.temperature - T_INF, st.temperature - T_INF)` (`MAT-12`'s `:238` form, with T∞ ≠ 0), then re-run W1 and W2 per rule (i). Both fit (244 + 273 s). Attempt 3 of 4. Pre-register the trap as "any FFCx `.c` without a sibling `.so`", not "0-byte".

**Slot note (2026-10-01 21:00 CDT slot, not started — no attempt spent).** W2's `timeout -k 30 900` exceeds the 660 000 ms foreground Bash window. implementer-run.md caps container time at ≈ 590 s, so the call would be killed host-side mid-window. A review should re-price W2 at ≤ 590 s, or split (∞) into a third window, before a slot takes this item.

**Read 2026-10-02 03:00 review.** The diagnosis holds on the parked module (`7120eeed`): `_rel(sub, th_h, th_ex)` is ‖th_h − th_ex‖/‖th_ex‖ (`:125`), and `:290` passes th_h = T_be − T_st and th_ex = T_st − T∞. So (∞) reads ‖(T_be − T_st) − θ_st‖/‖θ_st‖ = 1 − O(‖T_be − T_st‖/‖θ_st‖), and the printed 1 − 1.9e-14 says T_be = T_st to ≈ 1.9e-14. Correcting the comparand's arguments is a bug fix, not a band move: `INF_BAND` 1e-8 and the quantity (BE(∞) against the steady solve, `MAT-12`'s `:238` form) are unchanged. Every other anchor and both controls were green on their predictions. Re-queued as item 168 (attempt 3 of 4).

### §9 item 168 — 2026-10-02 (chain P2b, item 160 re-run with the (∞) argument fix; queued by the 03:00 review)

`MAT-15` step 2(b): the FEM on the phantom's own mesh against the finite Robin cylinder series, on item 160's split module with the one-line (∞) fix (implementer; heavy by the transient window, `-n 2`, complex; ≈ 16 measured / ≈ 25 predicted min; depends on nothing). Family `MAT-15` step 2: **attempt 3 of 4** since 2(a) registered its records (item 137). A red here leaves one attempt; the family cap then binds.

- **Code.** Take `tests/validation/test_mat15_finite_cylinder_fem.py` **by path** from `attempt/MAT-15-step2b-20261002T034707Z` (`7120eeed`). Also take, by path, the four log files that branch adds and `main` lacks: `20261001T202818Z_MAT-15-step2b-W1.log`, `20261002T033219Z_MAT-15-step2b-W2-steady.log`, and the two `…W3-transient.log` files. Append their four `test-results.md` rows **by hand** from `git diff main...attempt/MAT-15-step2b-20261002T034707Z -- docs/testing/test-results.md`, and never path-checkout `test-results.md`. Then make the one edit: `:290` becomes `out["inf"] = _rel(sub, res_inf.temperature - T_INF, st.temperature - T_INF)`, with a comment citing this item and the W3 log's (∞) line. Nothing else moves: no band, no node ID, no chunk size.
- **Anchors (asserted), unchanged from item 160.**
  - Precondition: 2(a)'s constants equal step 1's. z ∈ [−L, L] to 1e-9·L. max r ≤ a(1 + 1e-9). Centroid ≤ 1e-5·a.
  - (S1): LSQ rate ≥ 1.7, predicted 2.1511. Interior ≤ `INTERIOR_BAND` at 0.0025, predicted 7.374e-4.
  - (S2): Gauss ≤ 1e-9, predicted ≤ 1.3e-14.
  - (T): CN N = 40 interior ≤ `INTERIOR_BAND`, predicted 3.4945e-4.
  - (E): ≤ `E_BAND`, predicted ≤ 7.3e-14.
  - **(∞): ‖T_be − T_st‖/‖θ_st‖ ≤ `INF_BAND` 1e-8, predicted ≈ 1.9e-14** (the deficit item 160's misreading printed).
  - Cell count 58 866.
- **Negative controls (rule (e)).**
  - (c1): 2h comparand, **asserted** as a miss of 1 %, predicted 0.7394 (`…033219Z…W2-steady.log`).
  - (c2): Dirichlet roots, **asserted** as a miss, predicted 0.7795 (`…034149Z…W3-transient.log`).
  - Both are backed by item 160's measurement on this fixture.
  - Ceiling for (∞)'s own control: the pre-fix expression *is* the broken comparand, and it reads 1.0. Print it beside the fixed value as `(∞-control)`, not asserted.
- **Windows** (complex; foreground; Bash `timeout: 660000`; both container-side `timeout -k 30 590`; measured 244 s and 273 s).
  - First: `pgrep -c python3` in the container must read 0. Then run the any-orphan sweep as one bare call: `scripts/testing/run_and_log.sh MAT-15-step2b-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`.
  - W1, steady: item 160's W1 command verbatim, with label `MAT-15-step2b-W4-steady` and `timeout -k 30 590`. Use the node IDs as the parked module names them.
  - W2, transient: item 160's W2 command verbatim, with label `MAT-15-step2b-W5-transient` and `timeout -k 30 590`.
  - Both windows run after the edit, as rule (i) requires: W1's nodes do not touch `:290`, but the module changed.
- **Traps.** Item 160's list stands. The orphan sweep is "any `.c` without its `.so`", not 0-byte (`6388d5cd`). Every reduction takes `np.real`. Never `-k` inside the quoted command. Status 124 is a finding, not a retry. Make no edit after the last window.
- **Scope.** Uniform source, Robin walls, P1, four rungs. This item does not cover EM, perfusion, degree 2, or `ANS-9`.
- **Status it can move.** §7 `MAT-15` 🟡 → ✅, on the row's 2(b) Done-when, which this item matches (the centroid bar is a precondition, not a Done-when term). On the landing, delete both `attempt/MAT-15-step2b-*` branches.
- **Negative result.**
  - (∞) red after the fix: the BE run does not reach the steady state. That is physics or a scheme fault, not the comparand. Record the value and the phase clocks, park, and mark the item BLOCKED. Attempt 4 of 4 is then the review's to re-pose at a higher altitude.
  - Any other red: item 160's dispositions.

**Outcome 2026-10-02 (04:30 CDT slot) — ✅ COMPLETE.** Took the module and the four logs by path from `attempt/MAT-15-step2b-20261002T034707Z`, appended their four `test-results.md` rows by hand, and made the one edit: (∞) now reads `_rel(sub, res_inf.temperature - T_INF, st.temperature - T_INF)`, with the pre-fix expression printed beside it as `(inf-control)`. No band, node ID or chunk size moved. Pre-run: `pgrep -c python3` read 0, and the any-`.c`-without-`.so` sweep found no orphans (`20261002T093054Z_MAT-15-step2b-cache-sweep.log`). **W4 steady** `20261002T093152Z_MAT-15-step2b-W4-steady.log` (`-n 2`, complex, 252 s, 19 passed): cells 2 746 / 8 497 / 18 570 / 58 866; LSQ slope 2.1511, pairwise 2.2106 / 2.1415 / 2.1015; finest interior 7.373685e-4; Gauss ≤ 1.326e-14; centroid \|c\|/a 2.101e-6, z and r exact; θ_h(0) 2.980951 vs 2.983150 K; V_sub −8.6995e-4; (c1) 0.73935; max \|Im θ\| 0 (`:7171–7181`). **W5 transient** `20261002T093610Z_MAT-15-step2b-W5-transient.log` (`-n 2`, complex, 273 s, 16 passed): (T) CN N = 40 interior 3.494532e-4, whole 4.1617e-4; (c2) 0.77946; **(∞) 2.204e-14** (predicted ≈ 1.9e-14, bar 1e-8), (∞-control) 0.9999999999999808; (E) worst 7.234e-14 (CN), BE 1.326e-14; max \|Im θ\| 0 (`:1873–1883`). Every prediction landed. §7 `MAT-15` 🟡 → ✅; both `attempt/MAT-15-step2b-*` branches deleted.


## §7 row moved verbatim 2026-10-01 (rotate_plan_archive.py chunks)

 *(history and §9 briefs: `docs/planning/chunks/MAT-15.md`)* **Coil-driven heating of the gel phantom on F-small — Phase 8 subgoal (8b)** (opened 2026-09-30 by the weekly; §10 chain P1 / P2; the fixture is `ANS-2`'s 10 MHz loaded birdcage so that `ANS-9` can reuse it). *Step 1 (P1):* the coil-driven SAR of each of the four drives as the volumetric source of `MAT-12`'s transient solve on a phantom-only thermal domain with a Robin surface, gel properties, no perfusion. Asserted: the θ-scheme's discrete energy identity (stored + surface loss = deposited) to round-off; the first-step volume-mean slope against ⟨SAR⟩ / c inside `MAT-12`'s measured time-truncation bound; the four drives' temperature maps related by the C4 rotation at the imported 5 % band; a doubled-source control that misses the balance. *Step 2 (P2a / P2b):* a uniform source in the finite cylinder with Robin walls — the double eigen-series verified in numpy by residuals before any solve, then the FEM on the phantom's own mesh against it, band from `MAT-9`'s measured rate. **Done-when (§4):** step 1 asserted green ⇒ 🟡; step 2 ⇒ ✅. Not claimed: tissue, perfusion, an external check (`ANS-9`), human scale. *Pinned 2026-10-01 03:00 review (rule (j)):* step 1 on `ANS-2`'s 0.0025 m rung, gel constants k 0.6, ρc 4.18e6, h 10, T∞ = T(0) = 20, one normalisation G = 1 W / P₁ for all four drives; "deposited" = G·P_k from `mean_sar` on the parent mesh; the slope bound is `MAT-12`'s `S2_BAND` 1e-4 at Δt₁ = 2.5e-6 t_R with the surface-deficit estimate printed; "temperature maps" = the 10 g-ball means of θ at t\* = 0.2 a²ρc/k at the four rotated centres (§9 item 130). Step 2 splits: *2(a)* the numpy double eigen-series verified by residuals — recorded, no glyph (§9 item 137); *2(b)* the FEM against it ⇒ ✅ — *pinned 2026-10-01 10:30 review (rule (j), §9 item 145):* uniform q with step 1's gel constants; P1 on the tag-3 submesh at `phantom_resolution` 0.0075 / 0.005 / 0.00375 / 0.0025 (`GEO-27`'s rungs); comparand `robin_finite_cylinder_{steady,transient}_temperature`; steady whole-cylinder L2 LSQ rate against (V/N)^{1/3} ≥ 1.7 (`MAT-9` step 1's bar); interior ≤ `INTERIOR_BAND` (1 %) at 0.0025, steady and CN N = 40 at t*; Gauss balance ≤ 1e-9 per rung; per-step balance ≤ `E_BAND`; long backward Euler = steady ≤ `INF_BAND`; controls h → 2h and Dirichlet radial roots each miss `INTERIOR_BAND`. *Step 1 asserted green 2026-10-01 (§9 item 130):* (B) ≤ 1.7e-14, (E) ≤ 8.5e-12, (E′) ≤ 1.5e-14, (S2) 5.75e-6 (= the surface deficit; κ_s 1.65, below the predicted 2–5), (C4) θ 10 g pairs 0.020–0.049 %, controls (d) 1.000 / (f) 75.25–75.28 % (`20261001T123309Z_MAT-15-step1-W1.log:1924–1954`, `-n 2`, 218 s); step 2 open. *Step 2(a) recorded 2026-10-01 (§9 item 137, no glyph):* the numpy double eigen-series verified — Robin residuals ≤ 4.3e-13·qa, Gauss balance 7.6e-15, Duhamel 3.5e-14, limits ≤ 2.3e-9, θ_ss(0,0) 2.983 K (`20261001T143248Z_MAT-15-step2a-W1b.log`, `-n 2`, 49 s; (R0) norm and (R2) δ moved with measurement, disclosed); P2b queueable. *Step 2(b) attempt 1 parked 2026-10-01 (15:00 slot, §9 item 145, `attempt/MAT-15-step2b-20261001T203800Z`).* Every anchor it measured was green: LSQ rate ≈ 2.15, finest interior 7.37e-4, (T) 3.49e-4, Gauss ≤ 1.3e-14. Two things stopped it. The pre-registered centroid bar of 1e-9·a was unmeetable (measured 2.1e-6·a). And W1 was killed at 480 s inside (∞). **Ruled 18:00 review (rule (f); the weekly ratifies):** the centroid bar is re-registered at 1e-5·a, because a centroid offset moves θ by O(δ/a) ≈ 2e-6 against the 1e-2 band. The z and r bars stay at 1e-9. The windows split into steady and transient with per-phase clocks. Re-queued as §9 item 160 (attempt 2 of 4); history in `docs/planning/chunks/MAT-15.md`. | 🟡 

### Audit — 2026-10-02 10:30 daily review (`auditor`, closure `caa4fab6`)

**DEMOTE(8), on the tier word only. The review relabels the row instead of demoting it** (the precedent is the `OPS-27` pattern, which the 09-30 weekly overruled as "a tier word"). The row declared `standard` (180 s), but the chunk's own windows measured 252 s (`20261002T093152Z_MAT-15-step2b-W4-steady.log:7324–7325`) and 273 s (`20261002T093610Z_MAT-15-step2b-W5-transient.log:2022–2023`). Items 160 and 168 had already declared them heavy. The row's tier cell now reads heavy.

Checks 1–7 are clean:
- **Digits:** every claimed digit traces. The review re-read W5:1882, (∞) 2.204e-14 against the 1e-8 bound, with the pre-fix control at 0.9999999999999808.
- **Pre-fix red:** `20261002T034149Z_MAT-15-step2b-W3-transient.log:1882`.
- **Bands:** the imported bands are unchanged. The one loosened bar is the precondition centroid (1e-9·a → 1e-5·a), ruled 10-01 18:00 and awaiting the weekly.
- **Done-when:** met across steps 1, 2(a) and 2(b). The footers are `20261001T123309Z_MAT-15-step1-W1.log:2131–2132` and `20261001T143248Z_MAT-15-step2a-W1b.log:275–276`.

Caveat: the test file was untracked when the harness recorded HEAD `a959b466`. attempts.md records no edit after the windows.

### §9 item 236 — 2026-10-04 (chain P5)

*Drafted by a general-purpose agent for the 2026-10-04 18:00 daily review; the review reviews it.* Licence: §10 chain P5 (weekly 2026-10-03), quoted, not changed: "`MAT-15` step 3 — heating per accepted watt at 64 MHz on copper. `ANS-7`'s tuned fixture as the source, degree 2, exact rectangles: the energy identity, the calorimetric slope and the σ_y mirror image of the one-port drive at `MAT-15`'s imported bands. Target: identities. The first thermal figure on a tuned coil. Serial on nothing." The 2026-10-02 sheets-and-κ ruling binds: exact rectangles, κ off. No band, tolerance, Done-when or validation target moves.

`MAT-15` step 3: step 1's thermal identities, run once on the SAR of `ANS-7`'s tuned copper F-small birdcage, driven across C₁ at 64.0 MHz, normalised per accepted watt.
- **Executor:** implementer. One new test module; no `src/` change and no edit to any existing module.
- **Tier:** heavy, `-n 8`, complex build.
- **Class:** chain step (solve). **Time:** ≈ 23 measured / ≈ 45 predicted min.
- **Depends on nothing on `main`.** It needs two gitignored files that T21 wrote, `logs/ans7_exact/tune_p2.json` and `logs/ans7_exact/tune_inmodel_p2.json` (item 217 W1 / W2, 2026-10-04 09:07 / 09:11 CDT). Both were on disk when this brief was drafted. If either is absent, stop and take the next item. Do not re-create them: T21's W2 is red by pre-registration.

**What changes against steps 1–2.** Step 2 used no EM, so nothing from it is imported here.

| | Step 1 (§9 item 130) | This step |
| --- | --- | --- |
| Fixture | `ANS-2`'s solid-conductor F-small, `build_four_port_sweep(phantom_resolution=0.0025)`, 199 920 cells, 58 866 tag-3 cells | The ×1 exact-strip **hole** mesh (`_exact_geo(th14, comm, GATED_WIDTH_FRACTION)`, h = 0.015), 80 798 cells (`ANS6_EXACT_STRIP_HOLE_CELL_RECORD`), default phantom sizing. The tag-3 count has never been printed on this mesh |
| Conductor | Volumetric, σ = 800 S/m | Copper by Leontovich on tag 401, σ = 5.8e7 (`_route_a_solve(…, "cu")` builds the term) |
| Sheets | Ragged f = 0.5 cut, κ none | Exact CAD rectangles, κ off (`_build_ports(width_correction_kappa=None)`) |
| EM order | Degree 1 | **Degree 2** (thermal stays P1) |
| Circuit | Untuned: every port at 50 Ω, four single-leg drives at 10 MHz | Tuned: Z_C(C_mode,Cu) on P2–P4, C₁ across the P1 port (circuit layer, lossless), **one** drive (P1) at 64.0 MHz |
| Normalisation | G = 1 W / P₁, one deposited watt | G = 1 W / P_acc, **one accepted watt** at the coil terminals |
| Map identity | C4 across the four drives | **σ_y mirror** of the one drive |

**Which C, and what "tuned" means here.** The capacitor is **C_mode,Cu = 6.358914053477356e-11 F (63.589141 pF)**: the k = 1-block construction at 64 MHz on the degree-2 exact 4×4. It is read from `tune_p2.json["c_mode"]["cu"]` and never recomputed (`20261004T140347Z_ANS-7-item217-W1-cmode.log:2088`). Known-issues 2026-10-04 (T21) is open: with this C, the in-model zero of Im `Y_in` sits at 64 000 308.55 Hz, 4.82e-6 above 64 MHz against the 1e-6 `BANDS` band (`20261004T140738Z_ANS-7-item217-W2-inmodel.log:2049–2050`). Whether "tuned" means C_mode or C_in is the weekly's question. **This item does not claim the coil is tuned to 64.000 MHz.** It drives at 64.0 MHz at C_mode,Cu. That frequency is 308.55 Hz from the in-model zero, i.e. 9.4e-3 of the half-linewidth f₀/(2Q_A) = 32.9 kHz (Q_A 972.90). Normalising per accepted watt removes the mismatch, and the detuning moves the field shape at O((9.4e-3)²) ≈ 1e-4 (printed estimate). Every asserted anchor below is an identity that holds for any C.

**Code.** New `tests/validation/test_mat15_tuned_copper_heating.py`, `@complex_only`, with one module-scoped fixture so the EM solve runs once.
- **Opt-in knob `FEM_EM_MAT15_STEP3=1`.** Unset ⇒ every test in the module skips; any other value raises `ValueError`. That keeps default suites from collecting a `-n 8` heavy job that needs gitignored inputs. No existing module, helper default or record moves, so the change is unset-identical by construction.
- **Imports only, never `test_*` names.**
  - From step 1's module: `_tag_submesh`, `_cell_average_dg0_tag`, `_cn`, `_be_one`, `_e_prime`, `_max_imag`, `K_GEL`, `H_GEL`, `T_INF`, `A`, `T_R`, `T_STAR`, `N_MAIN`, `DT_1`, `D_CONTROL_FLOOR`, `F_CONTROL_FLOOR`, `CAD_AREA_M2`.
  - From `MAT-10`: `_robin_tags`, `_to_sub`, `_allsum`, `_say`, `B_BAND`, `ROBIN_TAG`.
  - From `MAT-12`: `E_BAND`, `S2_BAND`, `RHO_C`.
  - Ball means: `C4_COVARIANCE_BAND` and `PHANTOM_RHO_KG_PER_M3` (`test_birdcage_b1_plus_map`), `CENTRE_RADIUS_M`, `QUADRATURE_DEGREE`, `TEN_GRAM_KG`, `averaging_ball_radius`, `mean_sar`.
  - From `TH-14`: the module as `th14` (for `_exact_geo` / `_route_a_solve`), plus `_surface_loss_w`, `DRIVEN`, `COPPER_SIGMA` and `PHANTOM_CELL_TAG`.
  - From `test_birdcage_power_identity`: `_source_power_w`, `_sheet_field_dissipation_w`, `DISCRETE_IDENTITY_RTOL`.
  - From `test_ans7_exact_strip`: `_exact_geo`, `_fixture_tie`, `_assert_tie`, `_route_a_solve`, `OUT_DIR`, `GATED_WIDTH_FRACTION`.
  - From `test_ans7_sparse_route_b`: `_ANS7` (`z_cap`, `TUNING_TERMINATED_INDICES`), `NODE_S_RTOL`, `_p_inc_w`.
  - `_sheet_azimuth_deg` from `tests.mesh.test_birdcage_leg_offset`, and `SHEET_IFACE`.
- **Do not edit** step 1's module, `test_ans7_exact_strip.py`, `test_ans7_sparse_route_b.py` or the `ans:7` script. `ans:9` imports step 1's `_run_heated`, and item 230 has just re-issued `ans:7`. An edit to any of them owes those gates' windows under rule (a).
- **Fixture body, in this order.**
  1. Read the two T21 JSONs and fail if either is absent. Check: `degree == 2`; \|`f_hz` − 64e6\| < 1; `cells == 80798`; `tune_inmodel_p2.json["tune_mhz"] == [63.9, 64.0, 64.1]`; and its `c_mode` equal to `tune_p2.json["c_mode"]["cu"]` bit for bit. Take `S11_ref = s11_a[1]`, the 64.0 MHz point.
  2. `tie = _fixture_tie()`, then `geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)`.
  3. φ₁ = `_sheet_azimuth_deg(msh, tags_f, SHEET_IFACE + 1, comm)`, the bbox centre of an exact rectangle. This precondition runs **before** the solve.
  4. `s11, fields = _route_a_solve(th14, geo, specs, 64.0e6, "cu", c_f=c_mode, return_fields=True, degree=2)`.
  5. **Power.**
     - Rebuild the solve's specs exactly as `_route_a_solve` does: `dataclasses.replace(s, port_impedance_ohm=_ANS7.z_cap(ω, c_mode))` for i ∈ `TUNING_TERMINATED_INDICES`.
     - Build the sheets as `sp.sheet(driven=(sp.port_id == DRIVEN))` and set `acct = {"mesh": msh, "facet_tags": tags_f}`.
     - Compute P_src, ΣP_sheet, P_total = ½∫σ\|E\|², P_phantom (both `_vol(dx(PHANTOM_CELL_TAG))` and `mean_sar(...)["dissipated_power_w"]`), and P_surf = `_surface_loss_w(msh, facet_tags, z_s, e, comm)`, with z_s = (1 + j)·`surface_resistance_ohm`(ω, `COPPER_SIGMA`). This is `TH-14`'s (b) block, `test_th14_birdcage_copper.py:375–407`, verbatim.
     - P_acc = P_src − ΣP_sheet. C₁ is lossless, so P_acc is the power accepted at the coil terminals. Then G = 1 W / P_acc.
  6. **Thermal.** Step 1's per-drive loop body (`test_mat15_coil_driven_heating_phantom.py:217–301`) for this one drive, with dep = G·P_phantom (= η W):
     - submesh and Robin tags, q_h through `_cell_average_dg0_tag` → `_to_sub`;
     - (B), κ_s, `_cn`, (E′) through `_e_prime`, then `_be_one` and (S2);
     - control (d) on a doubled q_h;
     - max\|Im θ\|.
  7. **Ball means at t\*.** Eight centres at radius `CENTRE_RADIUS_M`, azimuths φ₁ + m·45° (m = 0 … 7), with step 1's compiled ball forms (`:231–261`). Print each ball volume beside (4/3)π·a₁₀³.

**Anchors** (asserted; every band imported, none restated):
- **(F)** `_assert_tie(cells, tie)`: 80 798 cells within `CELL_COUNT_BAND`, and all four N_i rectangles within `STRIP_BAND`. This is the ruling's CAD identity.
- **(Φ) precondition:** \|sin φ₁\| ≤ 1e-9, so the mirror through the driven leg is σ_y (y → −y). Step 1 printed P1 at 360.0° on its fixture. The bbox of an exact rectangle is a set of CAD corners.
- **(W) source wiring:** \|S11 − S11_ref\|/\|S11_ref\| ≤ `NODE_S_RTOL`. This is a deterministic re-solve of T21's own 64.0 MHz C_mode drive on the same mesh, at the same width and degree. It is backed by the same comparison on the same fixture, T21 (W) at 2.77e-14 (`…140738Z…-W2-inmodel.log:2045`), with S11_A printed at `:2047`.
- **(P) accepted-power identity:** \|P_src − ΣP_sheet − P_phantom − P_(Ω∖phantom) − P_surf\|/\|P_src\| ≤ `DISCRETE_IDENTITY_RTOL`. This is `TH-14`'s (b) on the tuned configuration, and it ties the accepted watt to the phantom's deposited share. Item 196 measured it at degree 2 on this mesh, in the 50 Ω-terminated configuration: 1.355e-12 (`20261003T124111Z_ANS-7-item196-W2-p2.log:2038`).
- **(B)** \|∫_sub q_h dV − G·P_phantom\|/(G·P_phantom) ≤ `B_BAND`.
- **(E)** Every thermal run's `max_balance_residual` ≤ `E_BAND`.
- **(E′) the energy identity:** ρc∫Δθ/Δt + ∮hθ^(n+½) = G·P_phantom, **from the EM side**, ≤ `E_BAND`.
- **(S2) the calorimetric slope:** \|ρc·θ̄¹/(Δt₁·G·P_phantom/V_sub) − 1\| ≤ `S2_BAND`, one BE step Δt₁ from θ = 0.
- **(M) σ_y mirror image:** 10 g ball means of θ at t\*. (M1) θ̄(φ₁ + 90°) against θ̄(φ₁ − 90°), and (M2) θ̄(φ₁ + 45°) against θ̄(φ₁ − 45°), each relative to the first ≤ `C4_COVARIANCE_BAND`.

**Imported vs printed — the bands, read against rule (f).**
- **Imported:** `B_BAND` 1e-10 (`MAT-10`) and `E_BAND` 1e-9 (`MAT-12`) are pre-stated round-off bands. Neither depends on EM order.
- **Imported:** `S2_BAND` 1e-4 is `MAT-12`'s pre-stated time-truncation bound, on the P1 thermal side. EM order enters only through κ_s. The deficit is ≈ 3.44e-6·κ_s, so the bound holds to κ_s ≈ 29.
- **Imported:** `C4_COVARIANCE_BAND` 5 % was pre-registered before any reading. Its provenance was *measured* at degree 1 (CG1 \|B₁⁺\| 2.19 % worst, `WF-6` step 1b, `test_birdcage_b1_plus_map.py:108–120`), but the value was never registered by measurement under rule (f). It can therefore be imported at degree 2. The asymmetry it bounds here comes from the phantom's P1 mesh, which does not depend on EM order.
- **Imported:** `DISCRETE_IDENTITY_RTOL` 1e-6 and `NODE_S_RTOL` 1e-9.
- **No band in `MAT-15` was registered by measurement at degree 1.** The rule (f) moves in its history (2(a)'s (R0) / (R2) parameters, 2(b)'s centroid bar 1e-5·a) belong to step 2 and are not used here.
- **Printed only:** step 1's **records** are `ANS-2`-rung numbers: `V_SUB_RECORD_M3`, `P1_RECORD_W`, κ_s 1.65, and mean θ(t\*) 1.153 K per deposited watt. On this mesh they appear only as references beside the new values, never as comparands, and this item registers no new record.

**Negative controls (rule (e)).**
- **(d) doubled source:** (E′) against G·P_phantom reads 1.000, by linearity on any fixture. Step 1 measured 1.000000 (`20261001T123309Z_MAT-15-step1-W1.log:1948`). Like step 1's `test_control_d_…`, it is asserted only as a miss, > `D_CONTROL_FLOOR`.
- **(f′) off-mirror pair:** θ̄(φ₁ + 90°) against θ̄(φ₁), relative to θ̄(φ₁).
  - Prediction: 40–80 %. Leg currents go as cos(φ − φ₁), so the gap voltages at φ₁ ± 90° are ≈ 0, and the eddy E_z ∝ x vanishes on x = 0.
  - Never measured, so the factor is printed beside the reading. Only the miss is asserted (> `F_CONTROL_FLOOR`), as step 1's (f) was.
  - Ceiling: ≤ 1 when θ̄(φ₁ + 90°) ≤ θ̄(φ₁), so the 5 % floor sits 20× under it.
- **(P-ctl):** the identity with ΣP_sheet dropped misses by ≈ P_sheet/P_src ≈ 0.998, because R_bare ≈ 0.08 Ω sits against the 50 Ω source sheet. Predicted, printed, never asserted.

**Predictions** (written first, printed):
- Tag-3 cells in [300, 3 000]. V_sub against πa²·2L = 2.26195e-4 m³: deficit in [−3e-2, 0].
- (W) ≤ 1e-12. (P) ≤ 1e-10.
- P_phantom, `mean_sar` against `_vol`: ≤ 1e-12.
- η = P_phantom/P_acc in [0.52, 0.62]. The estimate is Q_Cu/Q_PEC at C_mode,Cu, 969.21/1664.90 = 0.582 (T21 (T1), `…140347Z…-W1-cmode.log:2090–2091`). P_surf/P_acc = 1 − η to (P).
- L_S = 1 − \|S11\|² ≈ 4.05e-3 (from `:2047`'s S11_A) against L_F = P_acc/`_p_inc_w()`: \|L_S − L_F\| ≤ 1e-4 absolute. Printed only: in the 50 Ω configuration this is item 196's (L), and it is unmeasured here.
- (B) ≤ 1e-13. (E), (E′) ≤ 1e-12. κ_s in [1.5, 6]: a uniform E_z ∝ x alone gives 1.73, and the gap-voltage E adds surface weight.
- (S2) ≈ 3.44e-6·κ_s, in [5e-6, 2.1e-5].
- **The first thermal figure on a tuned coil:** mean θ(t\*) per accepted watt in [0.55, 0.75] K, which is η × step 1's 1.153 K per deposited watt (`:1929`). The peak 10 g θ̄ per accepted watt is printed with its azimuth.
- (M1), (M2) ≤ 1 %. The q_h 10 g mirror pairs are printed, predicted ≤ 3 %. The C2 pair θ̄(φ₁) against θ̄(φ₁ + 180°) is not a symmetry of a one-port drive; it is printed, predicted ≤ 2 %.
- max\|Im θ\| = 0.

**Windows** (`-n 8`; bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first).
- **Price**, measured at `-n 8` on this mesh:
  - pytest start-up + `tests/environment` + `_fixture_tie` ≈ 90 s (T21 W2: 261 s window less its 171.4 s fixture, `:2052`);
  - `_exact_geo` mesh 24.4 s;
  - one degree-2 C_mode drive 41.8–50.3 s (`:2052`);
  - degree-2 post-processing forms ≈ 10–40 s (first compile of the `mean_sar` and DG0-average forms at degree 2);
  - the thermal side ≈ 5 s on ≈ 10³ phantom cells. Step 1 ran four drives plus the control on 58 866 cells in 29 s at `-n 2`.
  - Total ≈ 200–260 s, or ≤ 365 s on item 214's ×1.4 slow box. That is inside the container's 590 s.
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0. Then run the any-orphan sweep (item 160's 44 KB partial `.c` trap) as one bare call:

  `scripts/testing/run_and_log.sh MAT-15-step3-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **W1** (≈ 200–260 s), the closing window:

  `scripts/testing/run_and_log.sh MAT-15-step3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS_SKIP_PRIVATE=1 FEM_EM_MAT15_STEP3=1 timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment tests/validation/test_mat15_tuned_copper_heating.py -v -s --tb=short'"`
- **No edit after W1** (rule (i)). Any edit re-runs W1. One window and one re-run fit the hour.

**Traps already paid for.**
- **Orphan `.c`:** a stem with no `.so`, not only 0-byte (item 160, 44 345 bytes). Sweep only when `pgrep -c python3` reads 0, and delete orphan stems only.
- **Never write under `logs/ans7_exact/`.** This item only reads T21's two JSONs and writes nothing outside its harness log.
- **The sheet terms of (P)** must use the Z_C-replaced specs. With the 50 Ω specs, P2–P4 would show a dissipation they do not have and (P) would go red on a correct field.
- **Complex build:** `ufl.real` inside `ufl.lt` for the ball indicator, with no `sqrt` in any comparison (`WF-6` 4c). `np.real` before every `float`.
- **Rank-local values:** `cell_tags.values` and `assemble_scalar`; every reduction goes through `_allsum` / `allreduce`.
- **No point evaluation on the submesh.** `evaluate_vector_field_parallel`'s owner rule raises on `create_submesh` meshes (known-issues 2026-10-04, `ANS-9`). All θ readings are assembled ball integrals.
- **Inherited from items 195 / 196 / 217:** κ `None`, `z_from_s` never applied to the sweep's `z_matrix`, and no `ans:7` `_stage_*` call.
- **Harness hygiene:** harness calls are bare and whole; no `-k` selection inside the quoted command; a killed window is not a result (run `pgrep -c python3`, report, never lengthen the timeout).

**Scope.** Covered:
- the ×1 exact-strip hole mesh, degree-2 EM, copper by Leontovich, κ off;
- C_mode,Cu on P2–P4 and C₁ across the driven port, one drive at 64.0 MHz;
- P1 thermal on the default-sized phantom submesh, with step 1's gel constants (k 0.6, ρc 4.18e6, h 10, T∞ = T(0) = 20), no perfusion, constant σ;
- normalisation per accepted watt.

**Not claimed:**
- that the coil is tuned to 64.000 MHz (T21 is open, and the weekly chooses C_mode or C_in);
- any absolute SAR, temperature or compliance figure;
- a closed form for the coil-driven θ, or any phantom or h convergence. The phantom is at default sizing, and the h-rung T12′ has not landed;
- the quadrature two-port drive, 127.74 MHz, or the PEC column;
- `ANS-9`, whose SPEC is `ANS-2`'s 10 MHz case, or any AED / Palace reading;
- tissue, ASTM F2182.

**Status it can move.**
- `MAT-15` is ✅ on steps 1–2 and stays ✅ whatever this step reads.
- On green, the §7 row gains a "step 3 record" sentence: the identities asserted on the tuned copper fixture at degree 2 on exact rectangles, and the per-accepted-watt figures printed (η, mean θ(t\*), the peak 10 g θ̄ and its azimuth). §10 P5 then closes; the weekly marks it.
- **Rule (j):** the reviewing daily adds one Done-when clause to the row: "step 3 (§10 P5): identities green ⇒ a record, no glyph".

**Negative result.**
- **(W), (F) or (Φ) red:** the fixture is not the one T21 measured, because a JSON, the mesh or `_route_a_solve` moved. Read nothing thermal, stop, BLOCKED.
- **(P) red:** the tuned degree-2 Route A configuration leaks power. Print every term, open a known-issues entry and park BLOCKED. Never rescale G to absorb it.
- **(B) red:** the DG0 transfer of a degree-2 field on the hole mesh. Known-issues entry, BLOCKED.
- **(E) or (E′) red:** the thermal path on this submesh. BLOCKED, and step 1's ✅ stands.
- **(S2) red with (E) green:** print κ_s and its estimate. If the miss equals the deficit (κ_s > 29), it is a step-size finding: report it, move no band, BLOCKED.
- **(M) red:** print the q_h pairs, ball volumes and tag-3 cells.
  - If the q_h pair is also red, the EM field is not mirror-symmetric: known-issues entry, BLOCKED.
  - If only θ is red, on fewer than ≈ 1 000 phantom cells, it is the phantom mesh: record it, BLOCKED. The review may re-pose on a phantom-refined exact hole mesh. That is a new fixture, so it needs its own C_mode.
  - No band moves.
- **(d) or (f′) does not miss:** stop.
- No `PREDICTED RED:` line. Every asserted anchor is predicted green.

**Outcome (2026-10-04 22:30 slot): GREEN, 13/13** (`20261005T040135Z_MAT-15-step3-W1.log`, 156 s window, fixture 90.2 s, `-n 8`, complex; cache sweep `20261005T040008Z_MAT-15-step3-cache-sweep.log`, no orphans). Asserted: (F) 80 798 cells; (Φ) φ₁ = 0.000000000000°; (W) 1.293e-12; (P) 7.387e-13; (B) 7.687e-16; (E) worst 5.310e-11 (BE); (E′) 3.171e-14; (S2) 4.8742e-6; (M1) 1.3599 %, (M2) 1.8431 %; controls (d) 1.000000, (f′) 35.78 %. Printed: P_acc 1.0068e-5 W, η = 0.577681 (estimate 0.582), P_surf/P_acc = 1 − η to (P), P outside phantom 0 (copper is a surface term); (P-ctl) 0.9973; L_S − L_F 1.54e-5; tag-3 cells 537, V_sub −2.947e-2 against πa²·2L; κ_s 1.3284; **mean θ(t\*) 0.687806 K per accepted watt**; peak 10 g θ̄ 0.9426 K at φ₁ = 0°; C2 pair 1.29 %. **Predictions missed (printed, none asserted):** (W) 1.29e-12 over its ≤ 1e-12; κ_s 1.33 under [1.5, 6] and so (S2) 4.87e-6 just under [5e-6, 2.1e-5] (deficit = 3.486e-6·κ_s, as derived); (M1)/(M2) over the ≤ 1 % prediction (inside the 5 % band); the q_h (M1) pair 22.7 % against ≤ 3 % — at φ₁ ± 90°, where q_h is at its minimum (665 / 817 W/m³ against 4 170 at φ₁), on 537 cells; the q_h (M2) pair 1.01 %; (f′) 35.8 % under the 40–80 % prediction. Not a tuning claim (T21 open).

### §9 item 315 — 2026-10-08 (chain P15)

*Drafted by a general-purpose agent for the 2026-10-08 03:00 daily review; the review reviews it.* Licence: §10 P15 (weekly 2026-10-07, PROJECT_PLAN.md:4669–4674), quoted, not changed: "`MAT-15` step 4 — the Pennes sink on the coil-driven phantom. `MAT-11`'s perfusion term on P1's fixture: the energy identity with the sink at `E_BAND`, w → 0 reproducing P1's ⟨θ⟩(t\*) at `B_BAND`, θ monotone decreasing in w at three w values, the C4 identity; printed, ⟨θ⟩(t\*) per w. Target: identities and a limit. Serial on nothing. Heavy (P11's class, 266 s)." No band, tolerance, Done-when or target moves. A record on a ✅ row: no glyph moves.

`MAT-15` step 4 puts the Pennes sink w(T − T∞) of `solve_transient_conduction(..., perfusion_w_per_m3_k=w)` (`src/fem_em_solver/thermal/transient.py:111`, balance with the ∫w(T^(n+θ) − T∞) term at `:125–128`; `MAT-11` / `MAT-13`) on step 1's coil-driven phantom: `ANS-2`'s F-small, four single-leg drives at 10 MHz, σ = 800 S/m, the 0.0025 m rung, 58 866 tag-3 cells, G = 1 W / P₁ (`_run_heated`, `tests/validation/test_mat15_coil_driven_heating_phantom.py:180`).
- **Executor:** implementer. One new module `tests/validation/test_mat15_pennes_coil_driven.py`, `@complex_only`, one module-scoped fixture. No `src/` change (the sink is already a solver argument). **Do not edit** step 1's module (`ans:9` imports `_run_heated`; rule (a)) nor `test_mat18_wall_coil_driven.py`.
- **Tier:** heavy, `-n 2`, complex. **Time:** ≈ 17 measured / ≈ 30 predicted min. **Depends on nothing.**
- **Template:** item 270's module (`MAT-18` step 3, `test_mat18_wall_coil_driven.py:134–271`) is the same fixture with h replaced; copy its structure (local `_transient_h` → local `_transient_w`, ball-mean block `:162–191`), never import its `test_*` names.

**The w values (fixture parameters, not bands).** w ∈ {`W`/10, `W`, `W_STRONG`} = {24, 240, 2400} W/m³/K, `W` and `W_STRONG` imported from `test_mat11_perfusion_sphere.py:55–56`; and the limit value w_ε = 1e-8 W/m³/K. w_ε is non-zero on purpose: it takes the `perfused` branch of `_solve_transient` (`transient.py:219–220`), so (L0) is an independent limit, not the default-path identity item 270's (L0) turned out to be (its outcome caveat, `MAT-18.md:481`). **Gel has no perfusion:** these w are a numerical exercise of the sink on a coil-driven source, not a gel or tissue property.

**Code.**
- Fixture: `base = _run_heated(return_fields=True)` once (`sub`, `ft`, `G`, `P`, `v_sub`, `centres`, `q_h`, `T_star` by k).
- Local `_transient_w(sub, ft, q, w, …)`: step 1's `_transient` body with `perfusion_w_per_m3_k=w`, `h=H_GEL`; `_cn_w` = step 1's `_cn` (CN, N = `N_MAIN` = 40 to `T_STAR`, two BE half-steps, snapshots N−1, N).
- Local `_e_prime_w(sub, ds, res, deposited, w)`: step 1's `_e_prime` (`:164–171`) plus perf = w∫(½(θᴺ + θᴺ⁻¹))dx. Step 1's `_e_prime` has no sink term, so using it would fail (E′) on a correct field.
- Per drive k ∈ 0…3: CN at w_ε and at the three w; ⟨θ⟩(t\*) = ∫θ/V_sub; 10 g ball means at c_k (item 270's block). Print a `time.perf_counter()` per phase, `flush=True`.

**Anchors** (asserted; bands imported by name):
- **(E′) energy identity with the sink:** ρc∫(θᴺ − θᴺ⁻¹)/Δt + ∮hθ^(n+½) + ∫wθ^(n+½) = G·P_k from the EM side, ≤ `E_BAND` (1e-9, `test_mat12_transient_sar_heating_sphere.py:76`), all four drives × three w. **(E):** every run's `max_balance_residual` ≤ `E_BAND`.
- **(L0) w → 0:** ⟨θ⟩(t\*) at w_ε against the mean of `base["T_star"][k]`, relative ≤ `B_BAND` (1e-10, `test_mat10_sar_to_heat_sphere.py:61`), per drive. Predicted ≈ w_ε·t\*/(2ρc) ≈ 1.5e-12 plus round-off.
- **(Mo) monotone in w:** per drive, ⟨θ⟩(t\*) strictly decreasing over w ∈ {0 (P1), 24, 240, 2400}, and the own-centre 10 g θ̄ likewise. Asserted as strict inequalities between consecutive values, each gap > 1e-6 relative (≥ 10³× under the smallest predicted gap below).
- **(C4) with the sink:** the four drives' own-centre 10 g θ̄ in step 1's cyclic pairs, ≤ `C4_COVARIANCE_BAND` (5 %, `test_birdcage_b1_plus_map.py:120`) at each w. A uniform w does not break C4; predicted the order of step 1's 0.020–0.049 % (`20261001T123309Z_MAT-15-step1-W1.log:1924–1954`).

**Negative controls (rule (e)).**
- **(d) doubled source, asserted** (prior same-fixture measurement): (E′) with 2·q_h at w = `W` against G·P₁ misses; step 1 and item 270 read 1.000000 on this fixture. Asserted > `D_CONTROL_FLOOR`.
- **(P-ctl) sink dropped from the identity, predicted** (no prior measurement): (E′) at w = `W` with the ∫wθ term omitted misses by the perfused share ≈ w⟨θ⟩V_sub/(G·P_k) ≈ 240 × 1.11 × 2.262e-4 ≈ 0.060. Range [0.04, 0.08]; ceiling < 1 (the sink cannot exceed the deposit while θ rises). Asserted only as a miss > `D_CONTROL_FLOOR`.

**Predictions** (printed; lumped estimate f(x) = (1 − e⁻ˣ)/x with x = (λ_R + w/ρc)t\*, λ_R t\* ≈ 0.28 fitted to P1's 1.153 K against the adiabatic 1.326 K per watt, item 270's arithmetic): ⟨θ⟩(t\*)/P1 ≈ 0.996–0.998 (w = 24), 0.955–0.975 (240), 0.68–0.77 (2400); the same ratios on the 10 g θ̄ within a few %. (E′), (E) ≤ 1e-13. max|Im θ| = 0. **Printed records:** ⟨θ⟩(t\*) per w per drive and the ratio to P1, peak 10 g θ̄ per w, the perfused power at t\* (`result.perfused_power`), the C4 and q_h pairs.

**Window** (bare, whole, foreground; Bash `timeout: 660000`; `-s`; `tests/environment` first).
- **Price:** item 270's window on this fixture, 266 s at `-n 2` with 8 CN runs + (d) (`20261006T200649Z_MAT-18-step3-W1.log`); step 1's CN at 0.0025 is 2.5 s (`20261001T123309Z_MAT-15-step1-W1.log:1929`). This item runs 16 CN + (d) + (P-ctl, no extra solve): ≈ +20–30 s. Predicted 270–330 s.
- **First:** `docker compose exec -T fem-em-solver pgrep -c python3` must read 0; then the any-orphan sweep (rubric item 4) as one bare call:

  `scripts/testing/run_and_log.sh MAT-15-step4-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **W1:**

  `scripts/testing/run_and_log.sh MAT-15-step4-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat15_pennes_coil_driven.py -v -s --tb=short'"`
- No edit after W1 (rule (i)).

**Traps already paid for.**
- The `$` in the sweep must stay escaped (`\$`): item 270's first sweep failed host-side on it (`MAT-18.md:481`).
- Orphan `.c` of any size (item 160, 44 KB): sweep only at `pgrep` 0, orphan stems only.
- Step 1's `_e_prime` and `_transient` have no w argument (above); the public solvers raise on w < 0 (`transient.py:145`, `steady.py:190`; `MAT-11`'s negative-w control uses a private path) — never pass a negative w here.
- `cell_tags.values` / `assemble_scalar` rank-local: reduce through `_allsum`. `np.real` before every `float`.
- No point evaluation on the `create_submesh` mesh (known-issues 2026-10-04): every θ reading is an assembled integral. `ufl.real` inside `ufl.lt` for the ball indicator; no `sqrt` in a comparison (`WF-6` 4c).
- No `-k` inside the quoted command; a killed window is not a result (`pgrep -c python3`, report, never lengthen the timeout).

**Scope.** Covered: the Pennes sink as a uniform numerical term at four w on P1's untuned 10 MHz σ = 800 S/m fixture; P1 thermal; gel constants; arterial = ambient. **Not claimed:** any gel or tissue perfusion value, a Pennes closed form on the cylinder, σ(T) or temperature-dependent w, a tuned-coil / copper / 64 MHz figure, an absolute temperature or compliance figure, `ANS-9`.

**Negative result.**
- **(E′) or (E) red:** print stor / conv / perf / deposited. If (E) is green and only (E′) red, the local `_e_prime_w` is wrong — fix in the new module, re-run. If (E) is red too, the solver's sink path on a submesh: known-issues entry, BLOCKED; `MAT-15` ✅ stands.
- **(L0) red:** the local helper is not step 1's run (args, half-steps, snapshots) — fix in the new module only. If it reads ≈ 1e-9–1e-8, report the value beside the w_ε·t\*/(2ρc) estimate; no band moves.
- **(Mo) red:** a sink that raises θ is a sign error — stop, known-issues entry, BLOCKED.
- **(C4) red:** print the q_h pairs; red too ⇒ the fixture moved (known-issues), else BLOCKED.
- **(d) or (P-ctl) does not miss:** the sink never reached the solve. Stop.
- No `PREDICTED RED:` line; every asserted anchor is predicted green.

**Status it can move.** None: `MAT-15` stays ✅ (Done-when: step 1 ⇒ 🟡, step 2 ⇒ ✅, PROJECT_PLAN.md:2534). On green the §7 row gains a "step 4 record" sentence (the identities with the sink, ⟨θ⟩(t\*)/P1 per w) and §10 P15 closes; the weekly marks it. **Rule (j):** the reviewing daily adds "step 4 (§10 P15): identities green ⇒ a record, no glyph" to the row's Done-when (step 3's clause, item 236, is the precedent).

**Review rulings (2026-10-08 03:00 daily review, item 315).** (6) w ∈ {0, 24, 240, 2400} accepted (`MAT-11`'s 240 / 2400 plus 24); the values are a numerical exercise, gel has no perfusion, and the brief says so. (7) No Done-when amendment: no glyph moves, so rule (j) does not apply; the outcome is recorded here and on the row's state line.

**Outcome (2026-10-08 09:00 slot): GREEN, 42/42** first run (`20261008T140322Z_MAT-15-step4-W1.log`, 303 s window, 54 passed with `tests/environment`, P1 fixture 168.8 s, sink block 62.9 s, each CN 2.7–3.1 s, `-n 2`, complex; cache sweep `20261008T140145Z_MAT-15-step4-cache-sweep.log`, no orphans). Fixture reproduced step 1: P₁ rel 2.65e-13 to record, V_sub 3.45e-13, C4 0.0403 / 0.0202 / 0.0287 / 0.0487 %. Asserted: (E′) with the sink ≤ 1.63e-14 on all 12 (k, w) runs *(corrected 2026-10-08 10:30 review from "1.67e-14", a transcription slip — the largest printed value is 1.623e-14 at W1 `:1985`; audit PASS)*; (E) worst 8.659e-14 over 16 runs (k3 at w_ε); (L0) 1.526e-12 / 1.506e-12 / 1.509e-12 / 1.522e-12 against the w_ε·t\*/(2ρc) estimate 1.50e-12; (Mo) relative gaps ⟨θ⟩ 3.43e-3 / 3.03e-2 / 2.51e-1 and 10 g θ̄ 3.14e-3 / 2.78e-2 / 2.32e-1 (identical across drives to ≈ 1e-7); (C4) θ̄ pairs ≤ 0.0489 % at w = 24 / 240 / 2400 (0.0389 / 0.0187 / 0.0287 / 0.0489 % at 2400); (d) 1.000000; (P-ctl) 0.059809–0.059810 on all four drives. Printed: ⟨θ⟩(t\*)/P1 = 0.996565 (w = 24), 0.966366 (240), 0.724055 (2400), all inside the predictions; 10 g θ̄ ratios 0.996856 / 0.96918 / 0.74481; perfused power at t\* (k = 0) 6.166e-3 / 5.981e-2 / 0.4495 W per watt deposited; peak 10 g θ̄ on k = 0: 2.5793 / 2.5077 / 1.9272 K; q_h pairs 0.0133 / 0.0084 / 0.0466 / 0.0683 %; max\|Im θ\| 0. No prediction missed. Gel has no perfusion: these w are a numerical exercise, not a property.

### §9 item 361 — 2026-10-10 (chain step P17, daily provisional; lane 3 / Phase 8 (8e))

*Drafted by a general-purpose agent for the 2026-10-10 18:00 daily review; the review reviews it.* Licence: §10 chain P17, written by this review as a **provisional** step under §9 rule (6) (daily, 2026-10-10), toward Phase 8 subgoal (8e) "an ASTM-style gel phantom at human scale, serial on H9" (PROJECT_PLAN.md:4284; §10 "(8e) stays serial on H9′", :4765–4766). H9′ is DONE (items 358 / 359, `20261010T170356Z_WF-10-step2-W1.log`, `20261010T171705Z_WF-10-step2s-W1.log`). **Target: identities — P5's set (item 236) transferred to H9′'s fixture, no new validation target.** The 2026-10-02 sheets-and-κ ruling, the capacitor-ESR ruling and the 10-10 10:30 rulings on items 358 / 359 (`WF-10.md`, "Review rulings (2026-10-10 10:30 review)") bind. **No band, tolerance, Done-when or validation target is set or moved; every band is imported.** A record on a ✅ row: no glyph moves.

`MAT-15` step 5: step 3's thermal identities (item 236), run on the coil-driven SAR of H9′'s tuned copper F-human — the two-port 90° drive at C_in,F-human — heating its own gel phantom, normalised per accepted watt, with the feed-bisector mirror of the two drive senses as the map identity.
- **Executor:** implementer. One new test module, `tests/validation/test_mat15_f_human_tuned_heating.py`; **no `src/` change** (standing rule (c) applies only if one proves necessary — disclosed, its gate re-run green in the same slot) and no edit to any existing module.
- **Tier:** heavy, `-n 8`, complex build; W0 (smoke, `-n 1`), then W1. **Class:** chain step (solve), median 21 measured min. **Cost:** ≈ 21 measured / ≈ 50 predicted min (15 fixed + ≈ 25 code + W0 ≈ 1 + W1 ≈ 4–5.5 + ≈ 5 spare). **If W1 has not started by minute 40, park** (rule (d)).
- **Depends on nothing open.** It needs two gitignored files, both on disk when this brief was drafted: `logs/wf12/h12c-cu-s32-a2.npz` (the builder asserts it, 10:30 ruling (5)) and 358's node JSON `logs/wf10/h9p-node-118.2437.json` (5 466 bytes, written 2026-10-10 12:19 CDT by item 359's builder call). Neither is tracked. **If either is absent, stop and mark the item BLOCKED; never re-run 358 or re-solve the 32×32 in-slot.**
- **Family:** `MAT-15` step 5, attempt 1 of 4.

**What changes against step 3.**

| | Step 3 (§9 item 236) | This step |
| --- | --- | --- |
| Fixture | `ANS-7`'s ×1 exact-strip hole F-small, 80 798 cells, 537 tag-3 cells | H9′'s copper F-human hole mesh, `_build_f_human(comm, as_hole=True)` through `build_tuned_quadrature`, 369 387 cells (`HOLE_CELL_RECORD`), phantom at the default 0.015 m: R 64.2857 mm, H/2 85.7143 mm, meshed volume 2.210693836095e-03 m³ (`20261010T171705Z_WF-10-step2s-W1.log:6313, :6377`) |
| EM order | Degree 2 | **Degree 1** (`DEGREE = 1`, `tests/validation/test_wf7_f_human_b1_quadrature.py:53`); thermal stays P1 |
| Circuit | C_mode,Cu on P2–P4, C₁ across P1, lossless, one drive at 64.0 MHz | C_in,F-human = 118.2437 pF with r = 0.2 Ω on all 32 gaps; P17 / P21 driven through Z_th; **two senses** ccw = (1, −1j), cw = (1, +1j), as superpositions of one node pair |
| Normalisation | G = 1 W / (P_src − ΣP_sheet) | G = 1 W / P_acc,ccw with P_acc defined below; **one G for both senses** (a per-sense G would hide a power asymmetry, step 1's rule) |
| Map identity | σ_y mirror of the one drive, 10 g θ̄ pairs at ±90° / ±45° | **Feed-bisector mirror M** (plane az₁₇ + 45° = 56.25°): θ_cw(Mx, t\*) = θ_ccw(x, t\*), 10 g θ̄ at c\*_θ and over 359's 18-point lattice |

**What "accepted watt" means here.** P_acc,s := P_src,s − P_sheet,s(P17) − P_sheet,s(P21), i.e. ½Re Σ_feeds V_port,i I_i\* with V_port,i = w_i − Z_th I_i: the power the two Thevenin feeds deliver across their sheets into the coil network. By (P_w) it equals ΣP_sheet over the 30 RC gaps + P_phantom + P_rest + P_surf. **Caveat, printed:** each feed gap's own r + C sits inside Z_th = 50 Ω ∥ Z_RC (`tests/validation/test_wf10_f_human_copper_esr_root.py:76`), so the two feed gaps' ESR loss is on the source side of this cut. P_acc is "accepted by the coil network beyond the feed gaps", not the 50 Ω-line accepted power. The omitted term Σ_feeds ½·r·\|V_port,i/Z_RC\|² is printed, together with P_acc + that term. θ is also printed per W of P_src and per W of P_phantom, so the record sits beside 359's SAR per W of P_src.

**Time grid — a fixture choice, not a band.** Step 1's helpers `_cn`, `_be_one` and `_e_prime` are imported unchanged, so the grid is step 1's **in seconds**: t\* = `T_STAR` = 0.2·a²ρc/k at F-small's a = 0.03 m, i.e. 1 254 s (≈ 21 min); N = 40; Δt₁ = `DT_1` = 0.015675 s (`tests/validation/test_mat15_coil_driven_heating_phantom.py:85–89`). It is not rescaled to F-human's t_R,H = a_H²ρc/k ≈ 28 791 s. Rescaling would need local copies of the helpers, and no identity depends on t\*. Printed: t\*/t_R,H ≈ 0.044 and Bi_H = h·a_H/k ≈ 1.07.

**Code.** New `tests/validation/test_mat15_f_human_tuned_heating.py`, `@complex_only`, with one module-scoped fixture so the builder runs once.
- **Opt-in knob `FEM_EM_MAT15_STEP5=1`** (step 3's `_enabled` / `_need` pattern, `test_mat15_tuned_copper_heating.py:112–123`). Unset ⇒ every gated test skips; any other value raises `ValueError`.
- **Imports only, never `test_*` names.**
  - From item 358's module (`tests/validation/test_wf10_f_human_tuned_b1_sar.py`): `build_tuned_quadrature` (`:147`), `C_IN_F_HUMAN_PF` (`:58`), `NODE_JSON` (`:82`).
  - From item 359's module (`tests/validation/test_wf10_f_human_tuned_sar10g.py`): `_lattice` (`:58`), `_containment` (`:85`), `_plane_distance` (`:80`), `AZ17_RECORD_DEG` (`:45`).
  - From H17's module: `ESR_OHM` (`:53`), `_z_rc` (`:120`), `_thevenin_feed_impedance` (`:76`). From `test_wf7_f_human_b1_quadrature`: `FREQUENCY_HZ`, `DEGREE`.
  - From step 1's module: `_tag_submesh` (`:106`), `_cell_average_dg0_tag` (`:124`), `_cn` (`:154`), `_be_one` (`:160`), `_e_prime` (`:164`), `_max_imag` (`:174`), `K_GEL`, `H_GEL`, `T_INF`, `T_R`, `T_STAR`, `N_MAIN`, `DT_1` (`:82–89`), `V_BAND` (`:98`), `D_CONTROL_FLOOR`, `F_CONTROL_FLOOR` (`:102–103`).
  - From `MAT-10` (`tests/validation/test_mat10_sar_to_heat_sphere.py`): `_say` (`:65`), `_allsum` (`:70`), `_robin_tags` (`:92`), `_to_sub` (`:122`), `ROBIN_TAG` (`:56`), `B_BAND` (`:61`). From `MAT-12`: `RHO_C` (`:68`), `S2_BAND` (`:75`), `E_BAND` (`:76`).
  - Bands and helpers: `C4_COVARIANCE_BAND`, `PHANTOM_RHO_KG_PER_M3`, `_relative_l2` (`tests/validation/test_birdcage_b1_plus_map.py:120, :139, :600`); `EXACT_IDENTITY_RTOL` (`tests/validation/test_birdcage_sar_mass_averaged.py:163`); `DISCRETE_IDENTITY_RTOL` (`tests/validation/test_birdcage_power_identity.py:227`); `REUSE_REPRODUCTION_RTOL` (`tests/validation/test_port_birdcage_ring_matrix.py:171`, = 1e-6 via `:164`); `UNIT_WEIGHT_CONTROL_RTOL` (`tests/validation/test_port_drive_superposition.py:303`); `Y0_CODING_BAND` (H17 `:55`); `PHANTOM_CELL_TAG` (`tests/validation/test_port_birdcage_lumped_column.py:101`); `_residual` (`tests/validation/test_port_circuit_layer_field.py:142`); `QUADRATURE_DEGREE`, `TEN_GRAM_KG` (`tests/validation/test_mass_averaged_sar_standard_masses.py:121, :124`); `_phantom_geometry` (`tests/validation/test_wf7_f_human_sar10g_hotspot.py:95`); `averaging_ball_radius`, `mean_sar` (`src/fem_em_solver/post/sar.py:214, :57`); `_rss_gib` (`tests/validation/test_wf12_f_human_hole.py:191`).
- **Do not edit** 358's module, 359's module, H17's module (`test_wf10_f_human_copper_esr_root.py`), step 1's module (`ans:9` imports `_run_heated`; rule (a)) or step 3's module. The ball-mean block of step 3 (`test_mat15_tuned_copper_heating.py:251–282`) sits inside `_run`, so it is **copied** into a local helper with a provenance comment (base commit and line span), never imported.
- **Record constant (printed comparand, not a band):** `V_PHANTOM_H9P_RECORD_M3 = 2.210693836095e-03`. It is the meshed tag-3 volume that 359's coverage line printed as mass / ρ (`…step2s-W1.log:6377`, 13 digits, ρ = 1000).
- **Fixture body, in this order.**
  1. **Read 358's JSON before the builder** (the builder rewrites it, as in 359, `test_wf10_f_human_tuned_sar10g.py:173–182`): rank 0 reads, `bcast`, assert `c_in_pf == C_IN_F_HUMAN_PF` exactly, keep `s_node`.
  2. `out = build_tuned_quadrature(comm, with_z50=False)`. It asserts (A0), (iv), az₂₁ − az₁₇ = 90°, (P_RC) and the node NEG inside (`:232–234, :325, :401–402`), and returns `e` (per sense), `e17`, `e21`, `sigma`, `parts` (per sense, with `p_sheet_by_id`), `cur`, `w`, `anchors`, `az_feeds`, `az_mirror`, `ctx`, `times`, `t_window` (`:482–489`).
  3. **Power per sense.** P_acc,s from `parts[s]` as defined above. G = 1 / P_acc,ccw, one value for all thermal runs. η_s = P_phantom,s / P_acc,s. Then (Σ_P): `mean_sar(e["ccw"], sigma=…, rho=PHANTOM_RHO_KG_PER_M3, cell_tags=…, subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]` against `parts["ccw"]["p_phantom"]`. Print the feed-gap ESR term.
  4. **Thermal.** `sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)`; `ft = _robin_tags(sub)`; V_sub, ∮ds. Print the per-rank submesh owned / ghost counts (one `gather`).
     - For each field f ∈ {ccw, cw, P17, P21}: compute `W, avg = _cell_average_dg0_tag(msh, cell_tags, PHANTOM_CELL_TAG, 0.5*sigma*ufl.inner(e_f, e_f))` **without G inside the form**, then `avg *= G` (see Traps), and `q_f = _to_sub(W, avg, sub, parent_ids)`.
     - Then (B), κ_s, `_cn`, (E′) through `_e_prime(sub, ds, res, G·P_phantom,f)`, `_be_one` and (S2), as step 3 (`:236–299`).
     - The P_phantom of P17 / P21 comes from `_vol` on `e17` / `e21` (two assemblies) or `mean_sar`, never from the ccw partition.
     - Control (d) runs on 2·q_ccw. Record max\|Im θ\|.
  5. **Ball means of θ at t\*** (and of q) on the submesh, through the copied step-3 forms (indicator `ufl.lt(ufl.real(ufl.dot(off, off)), a10**2)`, `QUADRATURE_DEGREE`, centre a `fem.Constant(sub, zeros(3))` whose value is set per ball):
     - the lattice `_lattice(az₁₇, z_c)` (18 points, z_c from `_phantom_geometry`) under ccw;
     - their images `_mirror_xy(·, az_mirror)` under cw;
     - c\*_θ = argmax of θ̄_ccw over the lattice;
     - the controls at c\*_θ and Mc\*_θ.

     Each ball volume is printed beside (4/3)πa₁₀³. Containment uses `_containment` against the measured bbox.

**Anchors** (asserted; every band imported by name, none restated):
- **(S_rep) source wiring:** `_residual(s_node, s_rec)` ≤ `REUSE_REPRODUCTION_RTOL` (1e-6). Same comparison, same fixture: 359 read 6.027e-13 (`…step2s-W1.log:6312`). Step 3's (W) analogue.
- **Builder anchors, asserted inside `build_tuned_quadrature` and re-asserted here from `out["anchors"]`:** (A0) 369 387 within `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`) and census none; (P_RC) ≤ `DISCRETE_IDENTITY_RTOL`; node NEG ≥ `BLIND_FACTOR` × `DISCRETE_IDENTITY_RTOL` (`tests/validation/test_wf12_f_human_copper_matrix.py:61`).
- **(P_w) accepted-power identity, both senses:** ≤ `DISCRETE_IDENTITY_RTOL`, as 358 / 359 (measured 5.087e-13 / 5.968e-13, `:6307, :6309`). **(Lin)** ≤ `UNIT_WEIGHT_CONTROL_RTOL` / `Y0_CODING_BAND`. With P_acc defined above, (P_w) is the step-3 (P) identity on this circuit.
- **(Σ_P)** P_phantom,ccw = `mean_sar` dissipated ≤ `EXACT_IDENTITY_RTOL` (1e-10); 359 read 0.0 (`:6378`).
- **(V) fixture tie:** \|V_sub / `V_PHANTOM_H9P_RECORD_M3` − 1\| ≤ `V_BAND` (1e-10, step 1 `:98`). Predicted ≤ 1e-13 (the submesh is the tag-3 cells).
- **(B), per field:** \|∫_sub q_f dV − G·P_phantom,f\| / (G·P_phantom,f) ≤ `B_BAND` (1e-10).
- **(E):** every thermal run's `max_balance_residual` ≤ `E_BAND` (1e-9).
- **(E′) energy identity from the EM side, ccw and cw:** ρc∫Δθ/Δt + ∮hθ^(n+½) = G·P_phantom,s, ≤ `E_BAND`.
- **(S2) calorimetric slope, ccw and cw:** \|ρc·θ̄¹ / (Δt₁·G·P_phantom,s/V_sub) − 1\| ≤ `S2_BAND` (1e-4). The deficit estimate is h·Δt₁·(∮ds/V_sub)/ρc·κ_s ≈ 1.60e-6·κ_s on the CAD ratio 42.8 m⁻¹, so the band holds to κ_s ≈ 62.
- **(M_θ-a) mirror at c\*_θ:** \|θ̄_cw(Mc\*_θ) − θ̄_ccw(c\*_θ)\| / θ̄_ccw(c\*_θ) ≤ `C4_COVARIANCE_BAND` (5 %).
- **(M_θ-b) mirror on the lattice:** `_relative_l2(θ̄_cw(M·), θ̄_ccw(·), all-true mask)` over the 18 points ≤ `C4_COVARIANCE_BAND`. Every point is off the plane (≥ 15 mm·sin 45°), which meets the 2026-09-26 weekly's off-plane condition.
- **(iv) containment:** every ball, lattice, images and controls, is ≤ 0 by `_containment` (1 mm margin, measured bbox). 359 read −4.881 mm worst over the same centres (`:6382`).

**Bands, read against rule (f).** Every band above is imported, and none was registered by measurement on this fixture. `C4_COVARIANCE_BAND` is the band 358 / 359 asserted on this mesh, and step 3 asserted it on θ̄ pairs. `B_BAND`, `E_BAND`, `V_BAND`, `EXACT_IDENTITY_RTOL` and `DISCRETE_IDENTITY_RTOL` are round-off bands; `S2_BAND` is `MAT-12`'s time-truncation bound. `V_PHANTOM_H9P_RECORD_M3` is a record comparand from a logged reading, not a band.

**Negative controls (rule (e)).**
- **Asserted (backed, same fixture):** the builder's node NEG, the P_surf-blind miss of (P_RC). 359 read 7.0125e-3 against the bar 10 × 1e-6 (`:6304`, its NEG line). The ceiling is P_surf/P_src ≈ 7.0e-3, so the bar sits ≈ 700× under it.
- **(d) doubled source — the miss asserted against a floor, the factor predicted.** (E′) on 2·q_ccw against G·P_phantom,ccw. By linearity it reads exactly \|2·dep − dep\|/dep = 1.000, which is also its ceiling. It is asserted only > `D_CONTROL_FLOOR` (1e-3), 1 000× under the ceiling. The factor 1.000 has never been measured on this fixture; it is printed. Step 1, 3 and 4 each read 1.000000 on F-small.
- **(f′) single-feed control — the miss asserted against a floor, the factor predicted.** The P17-only field with the same G, θ̄_P17(Mc\*_θ) against θ̄_P17(c\*_θ), relative to the latter. It is asserted only > `F_CONTROL_FLOOR` (= `C4_COVARIANCE_BAND`, 5 %).
  - **Predicted ≥ 100 %, printed.** The SAR analogue on this fixture read 3 675 % (359, `…step2s-W1.log`, control line). At t\* the diffusion length √(k t\*/ρc) ≈ 13.4 mm is ≈ a₁₀, well under the ≈ 63 mm between c\* and Mc\*, so θ̄ keeps most of the SAR contrast.
  - The quantity has no ceiling above. The floor sits ≥ 20× under the prediction.
- **Printed only (no same-fixture prior):**
  - the M-partner reading θ̄_P21(Mc\*_θ) against θ̄_P17(c\*_θ). It is predicted ≤ 10 %; 359's SAR analogue read 7.38 %, over its own prediction.
  - the **mis-paired** reading θ̄_ccw(Mc\*_θ) against θ̄_ccw(c\*_θ), which drops the sense swap. It is predicted at 3–25 %, from the 16.7–21.0 % azimuthal SAR spread of the z_c rows (359).
  - the P_surf-blind miss of each (P_w), predicted ≈ 7.0e-3.

**Predictions** (written first in the module, printed):
- Tag-3 submesh cells in [2 000, 8 000]: 537 cells on 2.26e-4 m³ at the same default sizing (step 3) scales to ≈ 5 300; 22 594 at 0.0075 m (`20261010T110256Z_GEO-35-s2bi-W1.log:10599`) scales to ≈ 2 800. V_sub against the CAD 2.225676e-3 m³: −6.73e-3, which is the coverage mass already logged. ∮ds against the CAD 2πa·2L + 2πa² = 9.521e-2 m²: deficit in [−2e-2, 0].
- (S_rep) ≤ 1e-12. (P_w) ≤ 1e-12. (Σ_P) ≤ 1e-14. (B) ≤ 1e-13. (E), (E′) ≤ 1e-12.
- η_ccw = P_phantom/P_acc in [0.22, 0.34]. The basis is the P17-alone node partition, 4.094e-3 / (4.094e-3 + 3.040e-4 + 1.0107e-2) = 0.282 (drafter's arithmetic on `…step2s-W1.log:6304`, unlogged). P_acc,cw / P_acc,ccw − 1 is predicted ≤ 1 %; the P_phantom pair reads 0.15 % (`:6307, :6309`).
- κ_s in [1.2, 2.5]: q ∝ r² gives 1.73 on this cylinder, and a coarse DG0 lowers it (step 3: 1.33 against 1.73). (S2) ≈ 1.6e-6·κ_s, in [1.9e-6, 4.0e-6].
- **Mean θ(t\*) per accepted watt** in [0.025, 0.048] K: η × t\*/(ρc·V) = η × 0.1348 K, less ≈ 5–10 % surface loss. Per deposited watt it is ≈ 0.12–0.13 K. The peak 10 g θ̄ per accepted watt and c\*_θ are printed. c\*_θ is predicted at r = 45 mm, z_c, az₁₇, where 359's SAR c\* sat. The θ̄ peak-to-mean is predicted in [1.1, 1.35] (SAR's 1.31, smoothed).
- (M_θ-a), (M_θ-b) in [0.2 %, 3 %]: 359's SAR10g mirror read 0.64 % / 1.28 %, and θ at t\* is that source smoothed over ≈ a₁₀. The q-ball mirror pairs are printed, predicted at 0.5–3 %.
- max\|Im θ\| = 0.

**Windows.** Foreground. The harness call is the whole Bash command, bare. Bash `timeout: 660000`; `-s`; `tests/environment/test_complex_mode.py` first.
- **Price, measured on this route:** the builder with `with_z50=False` takes ≈ 155 s: build 112.6 s, node pair 34.4 s, superposition 6.5 s (`…step2s-W1.log:6384`), plus pytest start and the harness at ≈ 10–15 s. The new work is:
  - four DG0 averages on the 369 387-cell parent, ≈ 1 s each;
  - five CN runs, two BE runs and ≈ 80 ball assemblies on ≈ 5 × 10³ submesh cells, ≈ 15–25 s (step 1's CN on 58 866 cells took 2.5 s at `-n 2`);
  - first compiles of any form not already cached, 0–60 s.

  Total ≈ 185–250 s, ≤ 330 s on the 1.32× box, inside the container's 590 s. Summed RSS ≈ 10 GiB (359: 9.95 GiB).
- **First:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Then run the any-orphan sweep as one bare call:

  `scripts/testing/run_and_log.sh MAT-15-step5-cache-sweep "docker compose exec -T fem-em-solver bash -lc 'cd /root/.cache/fenics && for c in libffcx_*.c; do s=\${c%.c}; ls \$s.*.so >/dev/null 2>&1 || { echo orphan \$s; rm -f \$s.c \$s.c.cached \$s.o; }; done; echo sweep done'"`
- **W0** (smoke, `-n 1`, numpy only). The ungated `test_mat15_s5_coding` must be green and the gated test must skip. The coding test asserts:
  - the knob parser: unset ⇒ off, `"1"` ⇒ on, anything else ⇒ `ValueError`, via `monkeypatch`;
  - on `_lattice(AZ17_RECORD_DEG, 0.0)`, `_mirror_xy(·, AZ17_RECORD_DEG + 45)` is an involution to 1e-12 m, and no point lies within 15 mm·sin 45° of the plane;
  - `T_STAR == 0.2*T_R` and `DT_1 == 2.5e-6*T_R` bit for bit.

  It prints t\*/t_R,H and Bi_H.

  `scripts/testing/run_and_log.sh MAT-15-step5-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -v -s --tb=short'"`
- **W1** (heavy, `-n 8`, predicted 190–330 s harness). Its durable capture is §5.1's idiom copied verbatim, trailing `; exit $rc` included, as in 358 / 359:

  `scripts/testing/run_and_log.sh MAT-15-step5-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf10 && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_MAT15_STEP5=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -s -v --tb=short > /workspace/logs/mat15-step5-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/mat15-step5-w1-raw.log; cat /workspace/logs/mat15-step5-w1-raw.log; exit \$rc'"`
- **No edit after W0 starts** (rule (i)); any edit re-runs W0 and W1. One W1 plus one re-run fit the hour. No rule (a) / (c) re-gate is owed: no `src/` file and no existing helper is edited.

**Traps already paid for.**
- **Orphan `.c` of any size** (item 160, 44 KB). Sweep only when `pgrep -c python3` reads 0, and delete orphan stems only. After a kill or status 124, run `pgrep -c python3` first, then the sweep under the label above; park, and never raise the timeout or `-n`.
- **Complex build.**
  - No `ufl` ordering on complex operands, and no `sqrt` inside a comparison (`OPS-22`; `WF-6` 4c). The ball indicator compares `ufl.real(dot(off, off))` with a₁₀², as step 3 does.
  - Ball centres go in as plain floats through `Constant.value`.
  - Call `np.real` before every `float`.
- **A Python float inside a UFL form is a literal in the generated kernel.** Step 3 put G inside the DG0-average integrand, so each G value is its own FFCx compile. Here G is applied to the averaged array instead, after the form, so the four fields share one compiled form.
- **Rank-local values.** `cell_tags.values` and `assemble_scalar` are rank-local; every reduction goes through `_allsum` / `allreduce`.
- **The `create_submesh` caveats (`MAT-19`).**
  - Every θ reading is an assembled ball integral; never `f.eval` or `evaluate_vector_field_parallel` on the submesh.
  - The phantom submesh may be owned by only some of the 8 ranks (`MAT-19` D4 at `-n 4`: owned 23 401 / 35 465 / 0 / 0, `20261006T005546Z_MAT-19-step2-a2-D4.log:1871–1875`). **Every rank enters every collective** — the CN / BE solves, `_allsum`, `mean_sar`, `_phantom_geometry` and each ball assembly. Never guard a collective with "owns phantom cells" (item 227's deadlock).
  - Ghosts propagate through `_tag_submesh`, so any owned-cell loop iterates `size_local` only.
- **Sense and plane.** ccw = (1, −1j) comes from the convention code (10:30 ruling (1)); never re-choose it. The mirror plane is `out["az_mirror"]` (az₁₇ + 45°), never `WF-7`'s az_ref plane. Hold E_ccw and E_cw as the builder returns them; never recombine e17 / e21 with hand-typed weights.
- **The builder rewrites `logs/wf10/h9p-node-118.2437.json`.** Read it first (step 1 of the fixture body). Write nothing else outside the harness log and the raw capture file.
- **Memory.** The builder holds E₁₇, E₂₁, E_ccw and E_cw (≈ 10 GiB summed). Free the DG0 parent arrays after `_to_sub`. Two heavy windows never run at once.
- **Harness hygiene.** Calls are bare and whole; no `-k` inside the quoted command; no pipe; `-s` always; foreground with `timeout: 660000`; never delete a log; never path-checkout a record file.

**Scope.** Covered:
- H9′'s fixture, unchanged: copper F-human hole mesh, degree-1 EM, κ off, r = 0.2 Ω on all 32 gaps, the two-port 90° drive at C_in,F-human = 118.2437 pF, 64 MHz;
- the existing **cylindrical** F-human phantom at the default 0.015 m (≈ 1.8 cells per 10 g ball);
- P1 thermal on its tag-3 submesh, with step 1's gel constants (k 0.6, ρc 4.18e6, h 10, T∞ = T(0) = 20), Robin on every exterior facet, no perfusion, constant σ, and step 1's time grid in seconds;
- normalisation per accepted watt as defined above.

**Not claimed:**
- **the ASTM F2182 box phantom**: that geometry is the weekly's to commission, and this step builds no new fixture;
- **(8e) closed**: (8e) asks for an ASTM-style phantom, and this is a step toward it;
- tissue, perfusion, or σ(T);
- an absolute temperature, ΔT or compliance figure;
- a closed form for the coil-driven θ, or h- or phantom-convergence;
- degree 2 (the 5.50 % human-scale order caveat, `WF-7` step 0b);
- f₀, Q, match or any tuning-accuracy figure (C_in is a record at 0.10 pF resolution);
- C16 invariance (the drive breaks it);
- `ANS-9`, or any AED or Palace reading.

**Negative result.** Report and stop. Park on `attempt/MAT-15-step5-*` and mark item 361 BLOCKED in the same commit (rule (d)).
- **(S_rep) red, or any builder anchor red** ((A0), (P_RC), node NEG, (P_w), (Lin)): the fixture or the copper route moved since 359. Read nothing thermal and stop. For a builder anchor, file a known-issues entry. Never re-band.
- **(Σ_P) or (V) red:** a defect in the SAR chain, or a phantom that is not 359's. File a known-issues entry; BLOCKED.
- **(B) red:** the DG0 transfer on the F-human hole mesh. File a known-issues entry; BLOCKED.
- **(E) or (E′) red:** the thermal path on this submesh partition. BLOCKED. Print the per-rank owned / ghost table. Step 1's ✅ stands.
- **(S2) red with (E) green:** print κ_s and its estimate. If the miss equals the deficit (κ_s > 62), it is a step-size finding: report it, move no band, BLOCKED.
- **(M_θ-a) or (M_θ-b) red:** print the per-point θ̄ and q̄ tables, the ball volumes and the submesh counts.
  - If the q̄ mirror is also red, the source is not mirror-symmetric on this phantom. That is P10's mechanism at 10 g; file a known-issues entry.
  - If only θ̄ is red, it is the P1 thermal mesh at ≈ 1.8 cells per ball: record it.
  - Either way BLOCKED, and no band moves. The re-pose (a refined phantom plus a re-root at C_in, an H17-class step) is the weekly's; never refine in-slot.
- **(iv) red:** a lattice or bbox error that 359 did not meet. Stop.
- **(d) or (f′) does not miss:** the source never reached the thermal solve, or the fields are swapped. Stop.
- **Status 124 or > 590 s:** BLOCKED without a result; the killed window is the price record.
- Every η, θ, κ_s, peak, mis-paired, partner or q̄ value is a reading, never a red.

No `PREDICTED RED:` line. Every asserted anchor is predicted green; a red anywhere stops the slot.

**Status it can move.**
- `MAT-15` stays ✅: its Done-when is step 1 ⇒ 🟡, step 2 ⇒ ✅ (§7 row, PROJECT_PLAN.md:2561), and this step reads nothing that touches it.
- **On green:** the §7 row gains a **"Step 5 record"** sentence: the identities on H9′'s tuned copper F-human at degree 1; the feed-bisector thermal mirror at c\*_θ and over the lattice; and, printed, η, mean θ(t\*) and the peak 10 g θ̄ with c\*_θ per accepted watt, on the 0.015 m cylindrical phantom. §10 P17 gains a DONE line; the weekly marks it.
- **No Done-when amendment.** No glyph moves, so rule (j) does not apply. Item 315's ruling (7) is the precedent.
- Nothing else moves: not the `WF-10` row (🟡, the weekly's rulings), not (8e), not §2. Known-issues changes only on a red.

**Execution order.** pgrep → cache sweep → W0 → pgrep → W1 → commit (module, harness logs, the row's Step 5 record and the §10 P17 DONE line together) → journal.

**Executor spawn prompt (verbatim).** "Execute §9 item 361 — `MAT-15` step 5 (§10 P17, provisional) — per `docs/planning/chunks/MAT-15.md` §9 item 361. Read the brief and any review rulings beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Pass timeout: 660000 on every Bash call that runs a window; run pytest with -s. Run the cache sweep, W0, then W1, in order, foreground. Before W1, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. Never delete a log. No edit to the module after W0 starts. Do not edit item 358's or 359's modules, `tests/validation/test_wf10_f_human_copper_esr_root.py`, or either existing `MAT-15` module. If `logs/wf10/h9p-node-118.2437.json` or `logs/wf12/h12c-cu-s32-a2.npz` is absent, stop (BLOCKED); never re-run 358. Every rank enters every collective. You are the executor — do not spawn agents, never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/MAT-15-step5-*`, mark item 361 BLOCKED in the same commit, and report."

**Review rulings (2026-10-10 18:00 daily review, item 361).** (361-1) **Normalisation:** per accepted watt, as step 3 (P5) did, with P_acc as the brief defines it. The feed gaps' ESR loss and the per-W-of-P_src figures are printed beside it so they compare with items 358 / 359. (361-2) **The time grid is step 1's** (t\* = 1 254 s). No identity depends on it, and t\*/t_R is printed. (361-3) **(V) stands.** It is a reproduction of the fixture's own logged volume (`20261010T171705Z_WF-10-step2s-W1.log:6377`, 13 significant digits) at step 1's `V_BAND`, so it is a comparand and not a new band. A red means a phantom other than 359's, and the item stops. (361-4) **Applying G after the DG0 average is an execution choice.** If the form compiles differently than the brief predicts, step 3's in-form construction is the fallback, and the journal says which one ran. (361-5) **Scope:** this is the first enumerated step toward (8e). It uses the cylindrical F-human phantom; the ASTM box geometry is a fixture the weekly commissions. `MAT-15` holds ✅ unchanged, and the result is a record in the row. The step is provisional under rule (6): the weekly ratifies or reverts it.

**Outcome 2026-10-10 (19:30 CDT slot) — ✅ DONE, first-run green.** Both inputs present; `pgrep -c python3` read 0 before the sweep and before W1. Sweep `20261011T003325Z_MAT-15-step5-cache-sweep.log` (no orphans, 97 s). **W0** `20261011T003506Z_MAT-15-step5-W0.log` (`-n 1`, 5 s): 5 passed, 20 skipped; involution 6.9e-18 m, min plane distance 10.6066 mm, T_STAR / DT_1 identities bit for bit, t\*/t_R,H 0.0436, Bi_H 1.0714. **W1** `20261011T003520Z_MAT-15-step5-W1.log` (`-n 8`, complex, 25 passed, pytest 159.9 s, harness 162 s; builder 152.7 s, thermal 4.3 s, 40 balls 0.8 s, summed RSS 9.59 GiB). Ruling 361-4: G applied **after** the DG0 average (one form for the four fields), not step 3's in-form fallback. Anchors (`:6296–6383`): (S_rep) 2.705e-13; (A0) 369 387; (P_RC) 7.739e-13; node NEG 7.0125e-3 (7012×); (P_w) 4.965e-13 / 5.842e-13 (blind 7.040e-3 / 7.035e-3); (Lin) 4.67e-16, Y 1.0e-15 / 1.2e-15; (Σ_P) 0.0; (V) 7.749e-14; (B) 0 / 7.9e-16 / 5.9e-16 / 5.9e-16 (ccw / cw / P17 / P21); (E) worst 3.269e-10 (be_ccw) over seven runs; (E′) 2.168e-13 / 3.513e-13 (P17 / P21 5.5e-13 / 1.2e-13, printed); (S2) 2.4658e-6 both senses (estimate 1.609e-6·κ_s = 2.41e-6); (M_θ-a) 1.1599 %; (M_θ-b) 0.7623 %; (iv) −4.8813 mm over 40 balls; (d) 1.000000; (f′) 697.40 %; max\|Im θ\| 0. Readings: P_acc,ccw 2.905666e-2 W (cw/ccw − 1 = +4.02e-4), feed-gap ESR term 1.265e-3 W, ½Re ΣV_port I\* differs from P_acc by 8.75e-3 rel (printed, not an anchor); η 0.2821 / 0.2824; sub cells 3 306, ∮ds −3.66e-3 vs CAD; κ_s 1.49; mean θ(t\*) 3.553e-2 K per accepted W (1.193e-2 per W of P_src, 0.1259 per deposited W); peak 10 g θ̄ 4.926e-2 K at c\*_θ = r 45 mm, z_c, az₁₇ (peak-to-mean 1.3865, over its predicted [1.1, 1.35]); q-ball mirror 1.47 % at c\*, 2.35 % rel L2; M-partner 3.53 %; mis-paired 13.21 %. MAT-15 holds ✅; a record.

### §9 item 413 — 2026-10-13 (§10 chain P22a; lane 1, no solve)

*Drafted by a general-purpose agent for the 2026-10-13 03:00 daily review; the review reviews it.* Licence: §10 chain P, "P22a (lane 1, no solve). The degree knob on `MAT-15` step 5's module (`FEM_EM_MAT15_STEP5_DEGREE`, build-only at `-n 2`, unknowns = 2E + 2F asserted, unset bitwise) — entry 38's prerequisite (i). Serial on nothing" (PROJECT_PLAN.md:4725–4727), and `xl` entry 38 (`docs/testing/xl-pending.md:1819–1829`). The pattern is H20a, `WF-10` step 4a, §9 item 368 (`docs/planning/chunks/WF-10.md` § "§9 item 368"; `tests/validation/test_wf10_f_human_tuned_b1_sar.py:94–169`, `:771–845`). **No band, tolerance, Done-when or target is set or moved.**

**What H20a already did, and what is left for `MAT-15`.** Step 5's module is `tests/validation/test_mat15_f_human_tuned_heating.py` (item 361; `MAT-20` step 3, item 380, added `FEM_EM_MAT15_STEP5_SERIES`). Its EM fixture is H20a's `build_tuned_quadrature(comm, with_z50=False)` (`:256`), which **already reads `FEM_EM_WF10_DEGREE` at call time** (`test_wf10_f_human_tuned_b1_sar.py:307–309`, `:332–343`): setting that env alone would run this module's EM at degree 2 today. What would then go wrong, and is this item's work:
- (S_rep) (`:244–264`, asserted at `:607`) compares the degree-2 `s_node` with the **degree-1** record `NODE_JSON` — a guaranteed red at degree 2. H20a ruled (S_rep) a degree-1 record comparand, printed only at degree 2 (`DEGREE2_RECORD_GATES`, `test_wf10_f_human_tuned_b1_sar.py:111–115`).
- The series knob writes `logs/mat15_s5/theta-series-ccw.json` (`:134`, `:589–592`), which `MAT-20` step 3 and `mat:18` (`EX-110` Part B) read — a degree-2 run would overwrite the degree-1 record.
- Nothing ties the label to the build: `FEM_EM_WF10_DEGREE=2` leaking into a run of this module (entry 38 may fold into entry 35's command, which sets it) would produce degree-2 numbers under degree-1 prints.

**`MAT-15` step 6a: `FEM_EM_MAT15_STEP5_DEGREE` on step 5's module, proven by a build-only window and a degree-1 re-gate.**
- **Executor:** implementer. **Tier:** heavy — W1 build-only `-n 2`, W2 item 361's closing window `-n 8`; W0 smoke `-n 1`. **Class:** lane 1 (median 6) / chain no-solve (median 7). **Cost:** ≈ 7 measured / ≈ 40 predicted min (15 fixed + ≈ 12 code + W0 ≈ 1 + W1 ≈ 2.5 + W2 ≈ 3 + spare). **If W1 has not started by minute 40, park** (rule (d)).
- **Depends on nothing open.** W1 and W2 read `logs/wf12/h12c-cu-s32-a2.npz` (the builder's (Z50) source, asserted present even on the build-only path, `test_wf10_f_human_tuned_b1_sar.py:300–301`) and W2 reads `logs/wf10/h9p-node-118.2437.json`; both are on the host at drafting. Absent ⇒ stop, BLOCKED; never re-run 358. **Family:** `MAT-15` step 6, landing.

**Code — `tests/validation/test_mat15_f_human_tuned_heating.py` only, additive. `src/`, H20a's / 359's modules, `test_mat20_tuned_heating_dose.py` and the other `MAT-15` modules are not edited.**
- **Knob.** `DEGREE_KNOB = "FEM_EM_MAT15_STEP5_DEGREE"`; `_step5_degree()` read at call time: unset or empty ⇒ `None` (the degree-1 path, untouched); `"2"` ⇒ 2; anything else, `"1"` included, ⇒ `ValueError`.
- **The guard.** `_step5_degree()` then reads H20a's `_degree_knob()` (private import from `test_wf10_f_human_tuned_b1_sar`, never a `test_*` name) and raises `ValueError` naming both variables unless the two agree (both unset, or both `2`). One degree, two switches that must be thrown together: entry 38's command (or entry 35's, if folded) sets both.
- **`_series_json_path(degree)`:** `None` ⇒ `SERIES_JSON` bit for bit; 2 ⇒ `logs/mat15_s5/theta-series-ccw-d2.json`; anything else ⇒ `ValueError`. The series write (`:590–594`) goes through it.
- **`STEP5_DEGREE2_GATES`:** `{"printed": ("(S_rep)",), "asserted": ("(A0)", "(P_RC)", "NEG node", "(P_w)", "(Lin)", "(Sigma_P)", "(V)", "(B)", "(E)", "(E')", "(S2)", "(M_theta-a)", "(M_theta-b)", "(iv)", "(d)", "(f')")}` — H20a's table in shape. (V) stays asserted: it is the tag-3 mesh volume against 359's logged value (ruling 361-3), degree-independent.
- **`_run()`:** call `_step5_degree()` first, before `:242`. At `None` the body is unchanged. At 2: read `_node_json_path(2)` (H20a, `:129–136`) before the builder if it exists and print the residual against it, else print "no degree-2 node record"; store (S_rep) as printed. `test_s_rep_source_wiring` at degree 2 calls `pytest.skip("(S_rep): degree-1 record comparand, printed only at degree 2 (STEP5_DEGREE2_GATES)")`; every other test is unchanged. The start print names the effective degree. `assert int(DEGREE) == 1` (`:242`) stays — it is the fixture's constant, which the knob does not change.
- **No `build_only` path is added to `_run()`.** The build-only proof goes through the builder's existing keyword in a new test (K1).

**Tests.**
- **(K0) `test_mat15_step6a_degree_knob_coding`** — ungated, numpy only, `monkeypatch`, both variables cleared first:
  - the parser: unset / `""` / `"  "` ⇒ `None`; `"2"` / `" 2 "` ⇒ 2 with `FEM_EM_WF10_DEGREE="2"`; `"1"`, `"3"`, `"x"`, `"2.0"` ⇒ `ValueError` — **the rejected values are the negative control**;
  - the guard: `MAT15=2` with `WF10` unset ⇒ `ValueError`; `MAT15` unset with `WF10=2` ⇒ `ValueError`;
  - `_series_json_path(None) == SERIES_JSON`; `_series_json_path(2)` ends `-d2.json` in the same directory and differs; `_series_json_path(1)` ⇒ `ValueError`;
  - the gate table: printed and asserted disjoint, (S_rep) printed, (V) and every thermal identity asserted;
  - the records' arithmetic `2_570_024 == 2 * 500_568 + 2 * 784_444` (H20a's `N_DOFS_D2_RECORD`, `E_RECORD`, `F_RECORD`, imported).
- **(K1) `test_mat15_step6a_build_only`** — selected by `FEM_EM_MAT15_STEP6A=1` (never `-k`); it does **not** request the `heated` fixture, so `FEM_EM_MAT15_STEP5` stays unset and none of the 25 step-5 tests runs. One build, both paths (H20a's K1 shape, `test_wf10_f_human_tuned_b1_sar.py:801–845`):
  - both knobs unset: `_step5_degree()` is `None`; `out = build_tuned_quadrature(comm, with_z50=False, build_only=True)`; asserted `out["degree"] == 1`, `out["ctx"]["solver"] is out["solver_built"]`, `n_dofs == E`;
  - `monkeypatch.setenv` both knobs to `"2"`: `_step5_degree() == 2`, then H20a's `_apply_degree(out["ctx"], comm)` on the same ctx; asserted degree 2 and `n_dofs == 2E + 2F`;
  - prints E, F, both counts, the build time and the summed RSS after the degree-2 space (`_rss_gib`).

**Anchors (asserted).**
- **(A0)** 369 387 cells within `CELL_COUNT_BAND`, (iv) none — the builder's own, in W1 and W2.
- **(Deg-1) W1:** unset ⇒ degree 1, the solver `_build_f_human` created (identity, no rebuild), unknowns = E. Record E = 500 568 (`20261011T140313Z_WF-10-step4a-W1.log:6217`, same mesh), printed beside.
- **(Deg-2) W1:** both knobs 2 ⇒ degree 2, unknowns = 2E + 2F exactly. Record 2 570 024 = 2·500 568 + 2·784 444 (`…140313Z…:6219`; `20261013T070008Z_WF-12-step4-H16.log:6532` at `-n 16`), printed beside; predicted equal.
- **(Bit) W2:** item 361's W1, re-run verbatim on the edited module with all knobs unset, is green — 25 passed, (S_rep) ≤ `REUSE_REPRODUCTION_RTOL` against 358's JSON read before the builder (the same comparison 361 read at 2.705e-13, `20261011T003520Z_MAT-15-step5-W1.log:6296`). **Printed beside** 361's readings (`:6296–6383`): (M_θ-a) 1.1599 %, (M_θ-b) 0.7623 %, mean θ(t\*) 3.553e-2 K and peak 10 g θ̄ 4.926e-2 K per accepted W, sub cells 3 306. "Bitwise" means reproduced at round-off (H20a's W2 reading); a printed reading off its record beyond the last printed digit is reported, not asserted.

**Negative controls (rule (e)).**
- **Asserted (backed):** a degree-1 space has E unknowns, not 2E + 2F; both numbers are measured on this mesh through this builder (`…140313Z…:6217–6220`). Assert `n(2) != n(1)` and `n(2) * E == n(1) * (2E + 2F)` in exact integers; the ratio 5.134216 is also the ceiling.
- **Asserted, the parser and guard:** K0's rejected values and both mismatches raise.
- **Predicted, printed:** the RSS after the degree-2 space ≤ 2 GiB summed (H20a 1.10 GiB, `:6221`) — a second solver held would show here.

**Windows.** Foreground, bare, repo-relative; Bash `timeout: 660000`; `-s`; complex.
- **W0** (smoke, `-n 1`; K0 and the existing coding test green, the gated tests skip — 361's W0 read 5 passed, 20 skipped):
  `scripts/testing/run_and_log.sh MAT-15-step6a-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -v -s --tb=short'"`
- **W1** (build-only, `-n 2`). Predicted 115–180 s harness: H20a's K1 build 118.3 s, harness 121 s on the same builder at `-n 2` (`20261011T140313Z_WF-10-step4a-W1.log:6217`). §5.1 durable capture verbatim:
  `scripts/testing/run_and_log.sh MAT-15-step6a-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_MAT15_STEP6A=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -s -v --tb=short > /workspace/logs/mat15-step6a-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/mat15-step6a-w1-raw.log; cat /workspace/logs/mat15-step6a-w1-raw.log; exit \$rc'"`
- **Before W2:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W2** (the degree-1 re-gate, `-n 8`, every knob **unset** but the step-5 selector; item 361's W1 command (above) with the label and raw-log name changed; predicted 155–220 s, as 361's 162 s):
  `scripts/testing/run_and_log.sh MAT-15-step6a-W2-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf10 && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_MAT15_STEP5=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -s -v --tb=short > /workspace/logs/mat15-step6a-w2-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/mat15-step6a-w2-raw.log; cat /workspace/logs/mat15-step6a-w2-raw.log; exit \$rc'"`
- **Why 361's window and not `MAT-20` step 3's** (`20261011T183420Z_MAT-20-step3-W1.log`, the module's latest window, series on): that window rewrites `logs/mat15_s5/theta-series-ccw.json`, which `MAT-20` step 3's W2 and `EX-110` Part B compare against a logged print. W2 must write nothing but the builder's node JSON (reproduced values, as 359 / 361 / H20a's W2 did). The series path's unset behaviour is K0's `_series_json_path(None) == SERIES_JSON` bit for bit.
- No module edit after W0 starts (rule (i)). Any edit ⇒ W0, W1 and W2 again.

**Traps already paid for.**
- **The degree-default trap** (`WF-15` step 1); never monkeypatch `DEGREE` in `test_wf7_f_human_b1_quadrature`; never hold two solvers (H20a's `_apply_degree` deletes first).
- **The builder rewrites `h9p-node-118.2437.json`** (W2 does, with reproduced values). Read it before the builder; at degree 2 the builder writes the `-d2` sibling (H20a premise (a)).
- **`build_only` returns before any `fem.form`** (no orphan `.c` risk); K1 must not request `heated`.
- **Every rank enters every collective** (item 227's deadlock; the submesh may be owned by some ranks only, `MAT-19` D4).
- **Orphan `.c` of any size** (item 160): after a kill, `pgrep -c python3` first; at 0, the rubric-4 sweep under label `MAT-15-step6a-cache-sweep`.
- Never `-k`; no pipe; never delete a log; never path-checkout a record file; two heavy windows never run at once.

**Scope.** A knob, its guard, two `-d2` paths and the dof count; the degree-1 path re-gated. **Not claimed:** any degree-2 θ, SAR or power reading (entry 38), entry 38's command or price, whether entry 38 folds into entry 35 (the 10-14 weekly's), the ASTM box, tissue, an absolute temperature, `ANS-9`.

**Negative result.** Report and stop. Park on `attempt/MAT-15-step6a-*` and mark item 413 BLOCKED in the same commit (rule (d)).
- **(Deg-1) red:** the knob leaks into the unset path. Park.
- **(Deg-2) red:** the space is not N1curl degree 2 first kind on this mesh. Known-issues entry; entry 38 stays PENDING.
- **W2 red, or a 361 reading not reproduced at round-off:** the edit is not bitwise. Park, revert nothing on `main`; the review rules.
- **W1 / W2 killed (status 124):** not a result; park with the measured time.

No `PREDICTED RED:` line. A red anywhere stops the slot.

**Status it can move.** No glyph: `MAT-15`'s Done-when is "step 1 asserted green ⇒ 🟡; step 2 ⇒ ✅" (§7, PROJECT_PLAN.md:2596), already ✅, and this item touches neither (rule (j) checked). It moves **a prerequisite of `xl` entry 38, not a glyph**:
- On green the §7 row gains one sentence: "**Step 6a** (§9 item 413, P22a): `FEM_EM_MAT15_STEP5_DEGREE=2` (guarded with `FEM_EM_WF10_DEGREE`) on step 5's module; unknowns 2E + 2F = <n> exact on the hole mesh (build-only, `-n 2`, <W1 log>); unset path re-gated (<W2 log>)."
- §10 P22a gains a DONE line.
- **`docs/testing/xl-pending.md` entry 38, same commit,** one edit only: in the Status paragraph, after "(i) §10 P22a (lane 1, no solve):", insert "**MET 2026-10-13 (§9 item 413, `<commit>`; `<W1 log>`, `<W2 log>`)** — `FEM_EM_MAT15_STEP5_DEGREE=2` must be set together with `FEM_EM_WF10_DEGREE=2` (the module raises otherwise); at degree 2 the series goes to `logs/mat15_s5/theta-series-ccw-d2.json` and (S_rep) is printed only —". **The implementer does not write the command and does not mark READY**; (ii) is entry 32's reading, which the 10-14 weekly reads.

**Execution order.** Code → W0 → W1 → `pgrep` → W2 → commit (module, harness logs, the §7 sentence, §10 P22a DONE, the entry 38 edit) → journal.

**Spawn line:** "Use the implementer agent to execute `MAT-15` §9 item 413 per `docs/planning/chunks/MAT-15.md` (read the latest Review ruling paragraph first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run every window in the foreground with `timeout: 660000`, pytest with `-s`, container `timeout -k 30`. Run W0, W1, then W2; before W2, `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0. Edit only `tests/validation/test_mat15_f_human_tuned_heating.py` (no `src/`, not H20a's or 359's module); never import a `test_*` name. If `logs/wf10/h9p-node-118.2437.json` or `logs/wf12/h12c-cu-s32-a2.npz` is absent, stop (BLOCKED). In `docs/testing/xl-pending.md` make only the one entry-38 insertion the brief names; never write its command or mark it READY."

**Outcome (§9 item 413, 2026-10-13 06:00 slot, run 11:00–11:10 UTC): DONE.** `DEGREE_KNOB` / `_step5_degree()` (guarded against `FEM_EM_WF10_DEGREE`), `_series_json_path`, `STEP5_DEGREE2_GATES`, (S_rep) printed-only and skipped at degree 2, K0 `test_mat15_step6a_degree_knob_coding`, K1 `test_mat15_step6a_build_only`, all in `tests/validation/test_mat15_f_human_tuned_heating.py`. No `src/` edit. W0 `20261013T110145Z_MAT-15-step6a-W0.log` (`-n 1`, 6 passed, 21 skipped). W1 `20261013T110156Z_MAT-15-step6a-W1.log` (build-only, `-n 2`, 116 s, 7 passed): (A0) 369 387 cells, (iv) none; (Deg-1) unset ⇒ degree 1, `_build_f_human`'s solver held, unknowns 500 568 = E; (Deg-2) both knobs 2 ⇒ 2 570 024 = 2·500 568 + 2·784 444 exact (record equal); NEG n(2)/n(1) = 5.134216 = (2E + 2F)/E asserted; summed RSS 1.13 GiB (predicted ≤ 2, printed). W2 `20261013T110401Z_MAT-15-step6a-W2-regate.log` (unset, `-n 8`, 166 s, 26 passed, 1 skipped = K1): (S_rep) 4.747e-13 (361: 2.705e-13; round-off), (V) 7.749e-14, (E) worst 2.608e-10, (M_θ-a) 1.1599 %, (M_θ-b) 0.7623 %, mean θ(t\*) 3.553024e-2 K, peak 10 g θ̄ 4.926111e-2 K per accepted W, sub cells 3 306 — every 361 reading reproduced to its printed digits. Entry 38 (i) MET; its command is the 10-14 weekly's.

## §7 row moved verbatim 2026-10-10 (rotate_plan_archive.py chunks)

 Coil-driven heating of the gel phantom on F-small | ✅ **Audit PASS 2026-10-02 10:30 review; tier relabelled heavy (chunk file).** 2026-10-02 (§9 item 168): step 2(b)'s P1 FEM on the phantom submesh matches the finite Robin cylinder series — LSQ rate 2.1511 (bar 1.7), finest interior 7.374e-4 and CN interior 3.4945e-4 (≤ 1 %), Gauss ≤ 1.3e-14, (∞) 2.204e-14 (≤ 1e-8), both controls miss (`20261002T093152Z_MAT-15-step2b-W4-steady.log`, 252 s; `20261002T093610Z_MAT-15-step2b-W5-transient.log`, 273 s; `-n 2`, complex). Uniform source, Robin walls, P1 only: no tissue, perfusion or external check (`ANS-9`). **Step 3 record** (§9 item 236, 2026-10-04, chain P5): on `ANS-7`'s copper F-small at T21's C_mode,Cu, driven across C₁ at 64.0 MHz (**not** a claim the coil is tuned to 64.000 MHz), degree-2 EM on the ×1 exact-strip hole mesh (80 798 cells), κ off, the identities hold per accepted watt — (W) 1.29e-12, (P) 7.39e-13, (B) 7.7e-16, (E′) 3.17e-14, (S2) 4.87e-6, σ_y mirror (M1) 1.36 % / (M2) 1.84 % (≤ 5 %), both controls miss; printed η = 0.5777, mean θ(t*) 0.688 K and peak 10 g θ̄ 0.943 K at φ₁ = 0° per accepted watt, on a 537-cell phantom submesh (`20261005T040135Z_MAT-15-step3-W1.log`, 156 s, `-n 8`). A record, no glyph. **Step 4 record** (§9 item 315, 2026-10-08, chain P15): the Pennes sink w(T − T∞) as a uniform numerical term (gel has no perfusion) on step 1's 10 MHz σ = 800 S/m fixture — the energy identity with the sink ≤ 1.63e-14 (four drives × w ∈ {24, 240, 2400} W/m³/K), every balance ≤ 8.7e-14, w → 0 (1e-8) reproducing P1's ⟨θ⟩(t\*) to 1.53e-12, ⟨θ⟩(t\*) and the 10 g θ̄ strictly decreasing in w, C4 pairs ≤ 0.049 %, both controls miss; printed ⟨θ⟩(t\*)/P1 = 0.9966 / 0.9664 / 0.7241 (`20261008T140322Z_MAT-15-step4-W1.log`, 303 s, `-n 2`). A record, no glyph. **Step 5 record** (§9 item 361, 2026-10-10, chain P17, provisional): H9′'s tuned copper F-human (degree-1 EM, 369 387-cell hole mesh, κ off, r = 0.2 Ω on all 32 gaps, two-port 90° drive at C_in,F-human = 118.2437 pF, 64 MHz) heats its own 0.015 m cylindrical gel phantom (3 306-cell P1 submesh), one G per accepted watt (P_acc,ccw = P_src − ΣP_sheet over the two feeds; the feed gaps' ESR term 1.27e-3 W printed beside it). The identities hold: (S_rep) 2.71e-13, (P_w) 4.97e-13 / 5.84e-13, (Σ_P) 0.0, (V) 7.7e-14 against 359's logged volume, (B) ≤ 7.9e-16 on ccw / cw / P17 / P21, (E) ≤ 3.3e-10 over seven runs, (E′) 2.2e-13 / 3.5e-13, (S2) 2.47e-6 both senses; the feed-bisector thermal mirror θ̄_cw(Mx) vs θ̄_ccw(x) reads 1.16 % at c\*_θ and 0.76 % rel L2 over the 18-point lattice (≤ 5 %); both controls miss ((d) 1.000000, (f′) 697 %). Printed: η = 0.2821, mean θ(t\*) 3.553e-2 K and peak 10 g θ̄ 4.926e-2 K at c\*_θ = r 45 mm, z_c, az₁₇ per accepted watt (peak-to-mean 1.39), κ_s 1.49 (`20261011T003520Z_MAT-15-step5-W1.log`, 162 s, `-n 8`). A record, no glyph; not the ASTM box, no absolute temperature or compliance claim. **Done-when:** step 1 asserted green ⇒ 🟡; step 2 ⇒ ✅. *History: `docs/planning/chunks/MAT-15.md`.* 

### Review ruling 2026-10-13 03:00 (daily) on the item 413 draft

Queued as drafted. `FEM_EM_WF10_DEGREE` already switches step 5's EM half to degree 2 (`test_wf10_f_human_tuned_b1_sar.py:307–309`), so the item adds a guard that the two knobs agree, a `-d2` series path (the degree-1 series that `mat:18` / `EX-110` read is never overwritten) and (S_rep) printed-only at degree 2; the re-gate is item 361's W1. The brief sits after item 361's outcome, ahead of the moved-§7-row section.

### §9 item 515 — 2026-10-18 (§10 P22b, `MAT-15` step 6b; lane 1, chain no-solve)

*Written for the 2026-10-18 03:00 daily review by a read-only drafter. Every file:line was read at drafting (`main` at `d94320db`); the review checks each command against its source before it queues the item.*

**Licence.** §10 chain P, "**P22b (lane 1, no solve).** Entry 42's prerequisite (rulings). Standard. Serial on nothing. **Then entry 42** (`xl` 10-19 or 10-20)" (`PROJECT_PLAN.md:4806–4807`). The 10-17 weekly's ruling defines the step: "the step-5 module takes its node through `FEM_EM_WF10_NODE` (H28's reader), build-only (Deg-2) at `-n 2`, the unset path bitwise" (`PROJECT_PLAN.md:4659–4669`). `xl` entry 42's Status gives the same text, plus "(the step-5 and 6a logs reproduced)" and "**The implementer landing P22b marks (i) MET in the same commit and does not write the command**" (`docs/testing/xl-pending.md:1999`). The pattern is H28, §9 item 511 (`docs/planning/chunks/WF-10.md` § "§9 item 511", landed `12df58ce`), applied to this module. The 2026-10-13 ruling on item 413 (above) binds. **No band, tolerance, Done-when or target is set or moved.**

**Deadline.** The item must be on `main` before the **2026-10-18 18:00 review**. If it is not, that review swaps `xl-queue.d/20261019-MAT-15-step6-P22.env` with the 10-20 file, so entry 37 runs Monday and entry 42 Tuesday (`PROJECT_PLAN.md:4667–4669`). A slot that cannot finish by then parks it, and does not hurry it through.

**`MAT-15` step 6b / P22b: the step-5 module reads its node from H28's `_node_pf()` in its power bookkeeping, its (S_rep) path and its prints; a build-only window at `-n 2` with the knob proves this; both degree-1 windows of step 6a are re-gated bitwise.**
- **Executor:** implementer.
- **Tier:**
  - W0 smoke (`-n 1`; 5 s measured, `20261013T110145Z_MAT-15-step6a-W0.log`);
  - W1 heavy by ceiling (`-n 2`, build-only, knobs on; 115 s measured for the same builder call at `d2root`, `20261017T170238Z_WF-10-H28-W1.log`);
  - W2 heavy by ceiling (`-n 2`, step 6a's W1 re-gate; 116 s measured, `20261013T110156Z_MAT-15-step6a-W1.log`);
  - W3 heavy (`-n 8`, step 6a's W2 re-gate; 166 s measured, `20261013T110401Z_MAT-15-step6a-W2-regate.log`).
  - §10 says "Standard". The windows still carry the 590 s container ceiling that items 413 and 511 used, because 166 s is within 8 % of the 180 s standard ceiling.
- **Class:** chain no-solve (median 18.6). **Cost:** ≈ 18.6 measured (class median; the weekly's own figure, `PROJECT_PLAN.md:4919`) / ≈ 28 predicted slot-minutes. The windows total ≈ 6.7 min (5 + 115 + 116 + 166 s), code ≈ 14 min, commit and journal ≈ 7. Item 511 ran in ≈ 13 min and item 413 in ≈ 10. **If W3 has not started by minute 40, park** (rule (d)).
- **Depends on:**
  - nothing in §9 (item 511 is DONE, `12df58ce`);
  - on disk, checked at drafting: `logs/wf12/h12c-cu-s32-a2.npz` (the builder asserts it even on the build-only path, `test_wf10_f_human_tuned_b1_sar.py:335–336`) and `logs/wf10/h9p-node-118.2437.json` (W3 reads it before the builder rewrites it, `test_mat15_f_human_tuned_heating.py:308–317`);
  - also on disk and **never written here**: `logs/wf10/h9p-node-118.2437-d2.json` (entry 35's, mtime 2026-10-16 02:19);
  - **absent, and correctly so:** `logs/wf10/h9p-node-115.5209-d2.json`. Entry 41 writes it on 10-22, and entry 42 writes it first, on 10-19 or 10-20 (see Traps). No window here needs it.
  - If either required file is absent, run W0 and W1, skip W2 and W3, and mark the item BLOCKED with "W2/W3 re-gate pending: <file> absent". Never re-solve to make one.
- **Family:** `MAT-15` step 6 (P22 → P22a → P22b); attempt 1 of 4.

**What the module is today** (`tests/validation/test_mat15_f_human_tuned_heating.py`, 873 lines). It imports `C_IN_F_HUMAN_PF`, `DEGREE_ENV`, `_degree_knob`, `_node_json_path` and `build_tuned_quadrature` from `test_wf10_f_human_tuned_b1_sar` (`:127–139`). It does **not** import `NODE_ENV`, `C_IN_F_HUMAN_D2_PF` or `_node_pf`.
- **The builder already follows the knob.** `build_tuned_quadrature` reads `c_pf = _node_pf()` at call time (`b1_sar.py:341–346`). It solves at that node, builds `out["parts"]` with it, and returns it as `out["c_in_pf"]` / `out["node_json"]` on both the full and the build-only path (`b1_sar.py:381–384`, `:639`). So **G = 1 / P_acc,ccw and every thermal number already follow the builder's node.** P_acc = P_src − ΣP_sheet over the two feeds, taken from `out["parts"]` (`:363–371`).
- **What does not follow the node:**
  - **(a) `:360` `z_rc = _z_rc(C_IN_F_HUMAN_PF * 1e-12, FREQUENCY_HZ)`, and `z_th` from it at `:361`.** These feed V_port = w − Z_th·I (`:368`), the printed cross-check ½ Re Σ V_port I\* (`:369`) and the feed-gap ESR term (`:370`). Both are printed (`:377–386`) and written to the series payload as `"p_vi_w"` / `"esr_w"` (`:663`, `:665`). At `d2root` both would use the wrong capacitor.
  - **(b) `:319` `_node_json_path(step5_degree)`, the one-argument form.** At degree 2 it reads entry 35's **118.2437** `-d2` record as the (S_rep) comparand (printed only), so a `d2root` run would print a meaningless residual against the wrong node.
  - **(c) the start print (`:301–303`) and the series payload (`:659–675`)** name no node.
  - **(d) no Y_in print and no SAR10g print.** Entry 42's readout asks for "Y_in at the node printed (an Im above 1e-3 S says the window is not at the root)" and "SAR10g(c\*) per accepted watt" (`xl-pending.md:2005`). The module prints neither today. The builder holds `out["y_in"][sense][i]` (`b1_sar.py:589`), and the module computes the 10 g ball mean of q at c\*_θ (`:542–547`, `:557`).
- **The existing guard.** `_step5_degree()` (`:166–182`) makes `FEM_EM_MAT15_STEP5_DEGREE` and `FEM_EM_WF10_DEGREE` agree. `_node_pf()` raises on `d2root` unless `FEM_EM_WF10_DEGREE=2` (`b1_sar.py:147–151`). The degree-1 path therefore cannot carry the knob.
- **The existing build-only K1** (`:823–873`) runs `monkeypatch.delenv(DEGREE_ENV)` and then calls the builder (`:834–838`). If `FEM_EM_WF10_NODE=d2root` is in the window environment, its `_node_pf()` raises, and K1 goes red. Item 511 fixed the same leak in WF-10's K1 with one `delenv(NODE_ENV)` line (`b1_sar.py:856`).

**Code. Edit only `tests/validation/test_mat15_f_human_tuned_heating.py`, additively.** Do not edit `src/`, `test_wf10_f_human_tuned_b1_sar.py`, `test_wf10_f_human_tuned_sar10g.py`, `test_mat20_tuned_heating_dose.py`, the root module, `tests/package/*`, any example, or `xl-queue.d/`. Never import a `test_*` name.
- **Imports (`:127–139`):** add `NODE_ENV`, `C_IN_F_HUMAN_D2_PF` and `_node_pf`.
- **New constant:** `STEP6B_KNOB = "FEM_EM_MAT15_STEP6B"`, the build-only selector. It is used in W1 only and never in entry 42.
- **New `_step6b_node()`:** returns `(degree, c_pf)` with `degree = _step5_degree()` called first (the guard keeps its order) and then `c_pf = _node_pf()`. When `degree is None` it asserts `c_pf is C_IN_F_HUMAN_PF`.
- **New `_s_rep_path(degree, c_pf)`:**
  - `None` ⇒ `NODE_JSON`;
  - `2` ⇒ `_node_json_path(2, c_pf=c_pf)`, which gives the 118.2437 `-d2` path character for character when `c_pf is C_IN_F_HUMAN_PF` (entry 38's behaviour, unchanged) and `h9p-node-115.5209-d2.json` at the root;
  - anything else ⇒ `ValueError`.
- **New `_node_impedances(c_pf, z0)`:** returns `(z_rc, z_th)` = `(_z_rc(c_pf * 1e-12, FREQUENCY_HZ), _thevenin_feed_impedance(z_rc, z0))`. It is the one place `:360–361` now live.
- **`_run()`:**
  - Replace `:300` with `step5_degree, c_pf = _step6b_node()`.
  - The start print adds the effective node, `NODE_ENV`'s raw value and the secant root.
  - On the degree-1 branch (`:308–317`) nothing changes; `:315` still compares against `C_IN_F_HUMAN_PF`.
  - On the degree-2 branch, `:319` becomes `node_json_d2 = _s_rep_path(step5_degree, c_pf)`. If the file is absent the branch prints "n/a" as today. **It never falls back to another node's record.**
  - Directly after the builder (`:333`), assert `out["c_in_pf"] == c_pf` exactly, and that `comm.allgather(c_pf)` holds a single value.
  - `:359–361` become `z0 = …; z_rc, z_th = _node_impedances(c_pf, z0)`. On the unset path `c_pf` is the same float object as `C_IN_F_HUMAN_PF`, so `z_rc` is bit-identical.
- **New prints (additive; every existing print stays as it is):**
  - **(Y) Y_in at the node** after the `G =` line (`:387`): `out["y_in"][s][i]` for both senses and both feeds, with the root reading named. That reading is P17 under cw = (1, +1j), the P17 (1, +1j) root (`b1_sar.py:74–78`). Print "|Im| > 1e-3 S ⇒ not at the root" as text only. **Printed, never asserted.**
  - **(SAR) SAR10g at c\*_θ per accepted watt** in the READING block (`:617–618`): the q-ball mean at c\*_θ divided by `PHANTOM_RHO_KG_PER_M3`, with per-W-of-P_src and per-deposited-W beside it as the θ readings have.
- **Series payload (`:659–675`):** add the additive keys `"c_in_pf"`, `"node_knob"` and `"degree"`. The path stays `_series_json_path(step5_degree)`: no `-d2` series exists on disk, and entry 38 (the detuned node) is superseded.
- **Docstring (`:43–50`):** one step-6b paragraph.
- **The existing K0 (`:778`) and K1 (`:823`)** each gain one `monkeypatch.delenv(NODE_ENV, raising=False)` line next to their existing delenvs, as item 511 did. Change nothing else in them.

**Tests added.**
- **(K0b) `test_mat15_step6b_node_coding(monkeypatch)`**, ungated and numpy only. It first runs `delenv` on `NODE_ENV`, `DEGREE_ENV` and `DEGREE_KNOB`, then asserts:
  - unset ⇒ `_step6b_node() == (None, C_IN_F_HUMAN_PF)`, with `is` on the float;
  - `d2root` with both degree knobs `"2"` ⇒ `(2, 115.5209)`, `is C_IN_F_HUMAN_D2_PF`;
  - `d2root` with both degree knobs unset ⇒ `ValueError` naming `DEGREE_ENV`;
  - `d2root` with only `FEM_EM_WF10_DEGREE=2` ⇒ `ValueError`, the step-6a guard firing first;
  - `"115.5209"`, `"x"` and `"d2"` ⇒ `ValueError`;
  - `_s_rep_path(None, C_IN_F_HUMAN_PF) == NODE_JSON`;
  - `_s_rep_path(2, C_IN_F_HUMAN_PF) == _node_json_path(2)`;
  - `_s_rep_path(2, C_IN_F_HUMAN_D2_PF)` ends `h9p-node-115.5209-d2.json` and differs from both;
  - `_s_rep_path(1, …)` ⇒ `ValueError`;
  - `_series_json_path(None) == SERIES_JSON`, so the `tests/package` importers' constant is unchanged;
  - **(N2) the closed form** at both nodes: `_node_impedances(c, z0)[0]` = `ESR_OHM − 1j/(2πf·c·1e-12)` to 1e-12 rel, and `[1]` = z0·z/(z0 + z) to 1e-12 rel.
- **(K1b) `test_mat15_step6b_build_only()`**, `@complex_only`, selected by `FEM_EM_MAT15_STEP6B=1` (never `-k`). It reads W1's command-line environment and does not monkeypatch. It does not request `heated`. It:
  - asserts `_step6b_node() == (2, C_IN_F_HUMAN_D2_PF)`, a wiring red if the window's environment is wrong (not a skip);
  - records the existence and mtime of `h9p-node-115.5209-d2.json` and `h9p-node-118.2437-d2.json` before the call;
  - calls `out = build_tuned_quadrature(comm, with_z50=False, build_only=True)` and asserts:
    - `out["degree"] == 2`;
    - **(Deg-2)** `out["n_dofs"] == 2E + 2F == N_DOFS_D2_RECORD` (2 570 024) exactly;
    - `out["c_in_pf"] == C_IN_F_HUMAN_D2_PF` exactly and `!= C_IN_F_HUMAN_PF`;
    - `os.path.basename(out["node_json"]) == "h9p-node-115.5209-d2.json"`;
    - **`_s_rep_path(2, _node_pf()) == out["node_json"]`**: the record `_run` would read is the one the builder writes;
    - `comm.allgather(_node_pf())` is identical on every rank;
    - `_node_impedances(out["c_in_pf"], z0)[0]` equals the closed form at 115.5209 pF;
    - both files' existence and mtime are **unchanged** after the call. Do not assert "absent": see Traps.
  - It prints the build time and summed RSS (H28 W1: 113.1 s, 1.10 GiB).

**Anchors (asserted; every band imported).**
- **(A0)** 369 387 cells within `CELL_COUNT_BAND` (`tests/mesh/test_birdcage_port_sheet_prerequisite.py:57`), in the builder, in W1–W3.
- **(Deg-2) W1:** 2 570 024 = 2·500 568 + 2·784 444 (`N_DOFS_D2_RECORD`, `E_RECORD`, `F_RECORD`, `b1_sar.py:109–111`). The same count was measured through this builder at `d2root` (`20261017T170238Z_WF-10-H28-W1.log:6230`) and on the unset-node path through this module (`20261013T110156Z_MAT-15-step6a-W1.log:6252`).
- **(N2)** the series-RC closed form through `_z_rc` (`test_wf10_f_human_copper_esr_root.py:188–191`) and `_thevenin_feed_impedance` (`:144–146`).
- **(Bit-W2)** step 6a's W1 re-run on the edited module with every knob unset is green, predicted 8 passed / 21 skipped (6a's 7 / 20 plus K0b and K1b). It reproduces (Deg-1) unknowns 500 568 = E with `_build_f_human`'s solver held (`…110156Z…:6250`), (Deg-2) 2 570 024 (`:6252`) and NEG n(2)/n(1) = 5.134216 (`:6253`).
- **(Bit-W3)** step 6a's W2 re-run with every knob unset but `FEM_EM_MAT15_STEP5=1` is green, predicted 27 passed / 2 skipped (6a's 26 / 1, `…110401Z_MAT-15-step6a-W2-regate.log:6581`, plus K0b passing and K1b skipping). (S_rep) ≤ `REUSE_REPRODUCTION_RTOL` (`test_port_birdcage_ring_matrix.py:171`) against 358's JSON read before the builder. The following are **printed beside, to the printed digit**:
  - **the lines that run through the edited `:360`:** ½ Re Σ V_port I\* 2.931095758904e-02 / 2.932256559685e-02 W and feed-gap ESR term 1.265338e-03 / 1.265081e-03 W, ccw / cw (`:6299–6300`). These are the bitwise test of the one line this item changes;
  - G = 3.441551785626e+01 1/W (`:6301`);
  - (Σ_P) 0.000e+00 (`:6302`);
  - (V) 7.749e-14 (`:6307`);
  - (E) worst 2.608e-10 (`:6320`);
  - (M_θ-a) 1.1599 % and (M_θ-b) 0.7623 % at `C4_COVARIANCE_BAND` (`test_birdcage_b1_plus_map.py:120`; `:6381–6382`);
  - mean θ(t\*) 3.553024e-02 K and peak 10 g θ̄ 4.926111e-02 K per accepted W (`:6386–6387`). Step 5's own window reads the same (`20261011T003520Z_MAT-15-step5-W1.log.gz:6379–6385`).
  - The new (Y) print at degree 1 is printed beside H28 W2's four Y_in readings: P17 cw +1.172678791e-02 +2.783852943e-06j S, and the others (`20261017T170443Z_WF-10-H28-W2-regate.log:6328–6331`).
  - "Bitwise" means reproduced at round-off, as item 413 read it. A printed reading that misses its record beyond the last printed digit is reported, not asserted.

**Negative controls (rule (e)).**
- **Asserted, by construction (closed form):** a knob-blind `:360` computes Z_RC at 118.2437 pF. At the root the separation Im Z_RC(115.5209) / Im Z_RC(118.2437) = 118.2437 / 115.5209 = **1.023570**, asserted in K0b and K1b against 1 for a blind module. The same comparison on the same closed form is measured at 1.023570 (`20261017T170226Z_WF-10-H28-W0.log`, item 511's NEG). The ceiling is the ratio itself.
- **Asserted, wiring:** K1b's `out["c_in_pf"] != C_IN_F_HUMAN_PF` and `_s_rep_path(2, _node_pf()) != _node_json_path(2)`. A module that kept the one-argument `:319` call reads the 118.2437 `-d2` path, which differs.
- **Asserted, coding:** K0b's rejected values and the two guard mismatches.
- **Predicted, printed only:** W1's summed RSS ≤ 2 GiB (`PREDICTED_RSS_D2_GIB`), against H28 W1's 1.10 GiB.

**Windows.** Each runs in the foreground with Bash `timeout: 660000`, with `-s`, in the complex build, with `tests/environment/test_complex_mode.py` first, and never with `-k`. **The harness call is the whole Bash command.**
- **W0** (smoke, `-n 1`; K0, K0b and the existing coding test green, the gated tests skip; predicted 7 passed / 22 skipped, against step 6a's 6 / 21, since K0b passes and K1b skips):
  `scripts/testing/run_and_log.sh MAT-15-step6b-W0 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 60 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -v -s --tb=short'"`
- **W1** (build-only, `-n 2`, all four knobs on the command line; predicted 110–200 s):
  `scripts/testing/run_and_log.sh MAT-15-step6b-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_WF10_DEGREE=2 FEM_EM_MAT15_STEP5_DEGREE=2 FEM_EM_WF10_NODE=d2root FEM_EM_MAT15_STEP6B=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -s -v --tb=short > /workspace/logs/mat15-step6b-w1-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/mat15-step6b-w1-raw.log; cat /workspace/logs/mat15-step6b-w1-raw.log; exit \$rc'"`
- **W2** (step 6a's W1 re-gate, `-n 2`, knobs **unset**; that command with only the label and raw-log name changed; predicted 110–200 s):
  `scripts/testing/run_and_log.sh MAT-15-step6b-W2-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_MAT15_STEP6A=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -s -v --tb=short > /workspace/logs/mat15-step6b-w2-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/mat15-step6b-w2-raw.log; cat /workspace/logs/mat15-step6b-w2-raw.log; exit \$rc'"`
- **Before W3:** `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` must read 0.
- **W3** (step 6a's W2 re-gate, item 361's window, `-n 8`, knobs **unset**; that command with only the label and raw-log name changed; predicted 155–230 s):
  `scripts/testing/run_and_log.sh MAT-15-step6b-W3-regate "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && mkdir -p /workspace/logs/wf10 && source /usr/local/bin/dolfinx-complex-mode && FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_MAT15_STEP5=1 PYTHONPATH=/workspace/src:/workspace timeout -k 30 590 mpiexec -n 8 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat15_f_human_tuned_heating.py -s -v --tb=short > /workspace/logs/mat15-step6b-w3-raw.log 2>&1; rc=\$?; echo \"[capture] rc=\$rc\" >> /workspace/logs/mat15-step6b-w3-raw.log; cat /workspace/logs/mat15-step6b-w3-raw.log; exit \$rc'"`
- No window sets `FEM_EM_MAT15_STEP5_SERIES`. The degree-1 series `logs/mat15_s5/theta-series-ccw.json`, which `MAT-20` step 3, `mat:18` and `tests/package/test_wf1_yaml_f_human_thermal*.py` read, is never rewritten.
- No module edit after W0 starts (rule (i)). Any edit means W0 to W3 again.

**Traps already paid for.**
- **The builder writes `out["node_json"]` on every full call** (`b1_sar.py:626–630`). W3 rewrites only the degree-1 `h9p-node-118.2437.json`, with reproduced values, the overwrite items 361, 413 and 511 accepted. W1 and W2 are build-only and return before the writer (`:378–384`). **No window here writes any `-d2` file.** Entry 35's `h9p-node-118.2437-d2.json` must keep its 10-16 02:19 mtime, and K1b asserts that.
- **Entry 42 writes `h9p-node-115.5209-d2.json` before entry 41 does.** Entry 42 runs on 10-19 or 10-20, and entry 41 is queued for 10-22 (`xl-queue.d/20261022-WF-10-step4-H28.env`). Entry 41's degree-2 b1_sar path reads nothing before its builder (`b1_sar.py:672–696`) and overwrites the file, so its readings are unaffected. But **H28's (K1n) asserts the file is absent after a build-only call** (`b1_sar.py:983`, `assert not exists_after`), so it goes red on any re-run after 10-19. K1b therefore asserts "unchanged", not "absent". The WF-10 module is out of scope here (flagged for the review).
- **The node-knob leak into K1.** W1's environment carries `d2root`, and K1 deletes `DEGREE_ENV` before its build. Without the new `delenv(NODE_ENV)` line, `_node_pf()` raises inside K1 if `FEM_EM_MAT15_STEP6A` were ever set beside it. W1 does not set it, but the line goes in (item 511's precedent).
- **The guard order.** `_step5_degree()` runs before `_node_pf()`, so a degree mismatch is reported as one before any node error.
- **The degree-default trap** (`WF-15` step 1). Never monkeypatch `DEGREE`, and keep `assert int(DEGREE) == 1` (`:304`). **Never hold two solvers**: `_apply_degree` deletes first.
- **The circular import.** Add no module-level import from `b1_sar` anywhere except this module's existing block.
- **Every rank enters every collective** (item 227); the `allgather` sits outside any `rank == 0` branch.
- **A killed window:** run `pgrep -c python3` first; at 0, run the rubric-4 cache sweep as `MAT-15-step6b-cache-sweep`.
- Never `-k`; no pipe; never delete a log; never path-checkout a record file; two heavy windows never run at once; durable capture ends in `; exit \$rc`.

**Scope.** The node knob reaches this module's Z_RC / Z_th, its (S_rep) path, its prints and its series payload. The item adds the (Y) and (SAR) prints entry 42's readout needs, a build-only proof at degree 2, and re-gates the degree-1 path. **Not claimed:** any degree-2 θ, SAR or Y_in reading (that is entry 42's); entry 42's command, tier or date; any `cem43` / six-minute reading on the degree-2 series; the WF-10 modules; f₀ / Q / match; the ASTM box; tissue; any absolute or compliance figure; any external comparison.

**Negative result.** Report and stop. Park on `attempt/MAT-15-step6b-*`, mark item 515 BLOCKED in the same commit, and leave entry 42 PENDING. The 18:00 review then applies the 10-19 / 10-20 swap.
- A K0b red means the knob, guard or path logic is wrong; W0 catches it before W1.
- A K1b red on `c_in_pf`, the path or the closed form means the module does not see the builder's node through `mpiexec`.
- A (Deg-2) red means the import or the call site is wrong; it also gets a known-issues entry.
- A W2 / W3 red, or a 6a or 361 reading not reproduced at round-off (above all the ½ Re Σ V_port I\* and ESR lines), means the edit is not bitwise. Revert nothing on `main`; the review rules.
- Status 124 is not a result; park with the measured time.

There is no `PREDICTED RED:` line, so a red anywhere stops the slot.

**Status it can move.** No glyph. `MAT-15`'s Done-when is "step 1 asserted green ⇒ 🟡; step 2 ⇒ ✅" (§7 row, `PROJECT_PLAN.md:2614`); it is already ✅, and step 6b is not in it (rule (j) checked). On green:
- The `MAT-15` row gains: "**Step 6b** (§9 item 515, P22b): `FEM_EM_WF10_NODE=d2root` (with both degree knobs at 2) reaches step 5's Z_RC / Z_th, its (S_rep) path `h9p-node-115.5209-d2.json` and its prints, with Y_in and SAR10g(c\*_θ) printed; (Deg-2) 2 570 024 at `-n 2` with the knob (<W1 log>); unset path re-gated (<W2 log>, <W3 log>)."
- §10 P22b gains a DONE line.
- **`docs/testing/xl-pending.md` entry 42, in the same commit:** in the Status paragraph, after "PENDING PREREQUISITE — §10 **P22b** (lane 1, no solve):", insert "**(i) MET <date> (§9 item 515, `<commit>`; <W1 log>, <W2 log>, <W3 log>)** —". **Do not write or edit the command, do not mark READY, and do not touch `xl-queue.d/`**: the replacement is the review's (`cmp` against entry 42). Report the exact knob and constant names as landed (`FEM_EM_WF10_NODE`, `d2root`, `C_IN_F_HUMAN_D2_PF`).

**Execution order.** Code → W0 → W1 → W2 → `pgrep` → W3 → commit (the module, four harness logs, the §7 sentence, §10 DONE, the entry-42 edit) → journal.

**Spawn line (if delegated):** "Use the implementer agent to execute `MAT-15` §9 item 515 (P22b) per `docs/planning/chunks/MAT-15.md` (read the brief and any review ruling beneath it first). The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped. Run W0, W1, W2, W3 in order, in the foreground, with `timeout: 660000`, pytest with `-s`, and container `timeout -k 30`. Before W3, run `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` and require 0. Edit only `tests/validation/test_mat15_f_human_tuned_heating.py`. Make no `src/` edit, never edit either WF-10 module, the `MAT-20` module or `tests/package/*`, and never import a `test_*` name. In `docs/testing/xl-pending.md`, make only the entry-42 (i) MET insertion the brief names; never write its command, never mark READY, never touch `xl-queue.d/`. You are the executor: do not spawn agents, and never background a harness window. If the docker socket is denied, stop and report the exact command string. On a red, park on `attempt/MAT-15-step6b-*`, mark item 515 BLOCKED in the same commit, and report. The item must land before the 2026-10-18 18:00 review."

**Outcome (2026-10-18, implementer, item 515): DONE, green.** W0 7 passed / 22 skipped (`20261018T094401Z_MAT-15-step6b-W0.log`, K0b's (N2) closed form at both nodes ≤ 2.3e-16 rel, NEG Im ratio 1.023570). W1 (knob on, `-n 2`, 114 s): 8 / 21; K1b (Deg-2) 2 570 024 exact, c_in 115.5209 pF on both ranks, node_json `h9p-node-115.5209-d2.json` = `_s_rep_path`, `-d2` records unchanged (115.5209 absent before and after, 118.2437 mtime held), build 110.4 s, RSS 1.10 GiB (`20261018T094411Z_MAT-15-step6b-W1.log`). W2 (unset, `-n 2`, 113 s): 8 / 21, (Deg-1) 500 568, (Deg-2) 2 570 024, NEG 5.134216 (`20261018T094613Z_MAT-15-step6b-W2-regate.log`). W3 (unset, `-n 8`, 162 s): 27 / 2; (S_rep) 4.955e-13; ½ Re Σ V_port I\* 2.931095758904e-02 / 2.932256559685e-02 W, ESR 1.265338e-03 / 1.265081e-03 W, G 3.441551785626e+01, (Σ_P) 0, (V) 7.749e-14, (M_θ) 1.1599 % / 0.7623 %, mean θ 3.553024e-02 K, peak 4.926111e-02 K: all to the printed digit. **Reported, not asserted:** (E) worst 3.469e-10 vs the 2.608e-10 record (band 1e-9). It does not pass through the edited line, and it already reads 3.310e-10 in `20261011T183420Z_MAT-20-step3-W1.log`, so it is thermal-solve round-off. The new (Y) P17 cw reading +1.172678791e-02 +2.783852940e-06j S sits beside H28 W2's ...943e-06j (last digits). SAR10g(c\*_θ) at degree 1 is 1.698396e-01 W/kg per accepted W (printed) (`20261018T094815Z_MAT-15-step6b-W3-regate.log`).

### Review ruling 2026-10-18 03:00 (item 515, P22b)

Accepted as drafted. Spot-check: `test_mat15_f_human_tuned_heating.py:319` (`node_json_d2 = _node_json_path(step5_degree)`, one argument) and `:360` (`_z_rc(C_IN_F_HUMAN_PF * 1e-12, …)`) were re-read.
- **Tier: heavy** (W3 measured 166 s against standard's 180 s). §10's "Standard" is a price label, and the 590 s ceiling is the safer wrapper.
- **The print-only Y_in and SAR10g lines are accepted.** They are additive, print only, and are what entry 42's Readout names.
- **Not this item's:**
  - Entry 42's `cem43` / six-minute row. The command cannot produce it, because the `MAT-20` dose module reads only the degree-1 series (`test_mat20_tuned_heating_dose.py:47`). It is noted in entry 42 for the weekly.
  - H28's build-only test, which asserts `h9p-node-115.5209-d2.json` absent after its call (`test_wf10_f_human_tuned_b1_sar.py:983`). It goes red only if `FEM_EM_WF10_STEP4C_BUILD` is set again after entry 42 has written the file. It is named for the weekly, and no slot sets that selector.
- **Deadline.** If this item is not on `main` by the 18:00 review, that review swaps the 10-19 / 10-20 queue files (entry 42's rule).
