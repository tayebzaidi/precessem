=== OPS-47
Closed 2026-09-13: the four open-chunk narratives (`POST-6`, `PORT-14`, `TH-15`, `WF-6`) were moved byte for byte to docs/planning/chunks/ with all anchors green, shrinking the plan to 5 546 lines. The 4 000-line guide is not met, a miss recorded with the guide not widened.

=== OPS-48
Closed 2026-09-16: the `EX-14` VTX reader was ported to the adios2 2.12 bindings and both examples now execute the read-back gate at relative difference 0 with the wrong-comparison control at 1. Audited PASS 2026-09-18; the writer-side hole that survey found went to `OPS-50`.
Done-when: the reader ported to the 2.12 bindings (`adios2.bindings.ADIOS` or `FileReader`, the VTX local-block walk kept), both examples run through the harness at `-n 2` with `relative difference ≤ VTX_ROUNDTRIP_RTOL` (1e-10, unchanged) **asserted and executed** (the log shows the `✓ written .bp reproduces` line, not the `⚠ unavailable` line), a wrong-comparison control (read-back vs 0.5 × in-memory) printed at rel ≈ 1, the guides' read-back sections re-cited to the new logs, censuses `broken=0` / docrefs `exit != 1`, elapsed recorded, the known-issues entry retired in the same commit.

=== OPS-49
Closed 2026-09-16: `write_xdmf_time_series` landed beside the byte-identical `write_xdmf_with_tags`, with the one-collection count identity, an exact h5py round trip and the collapse control asserted at -n 2. Audited PASS 2026-09-18; `EX-56` cannot adopt it and no follow-up is queued.
Done-when: additive `write_xdmf_time_series(filename, mesh, cell_tags, steps, comm)` with `steps = [(t, {name: Function}), …]`, one temporal collection whose `n` children each carry every field and CellTags; unit test at `-n 2`: (a) the XDMF parses to exactly one collection with `n` `<Time>` values equal to the `t`s written and every attribute on every child (in = out), (b) each step's arrays read back through `h5py` equal the gathered function arrays at rel ≤ 1e-12, (c) negative control — the same steps through the existing `write_xdmf_with_tags` at the pre-change commit (pinned `<sha>^`) collapse to ≤ 1 time value; `write_xdmf_with_tags` / `consolidate_xdmf_grids` byte-identical; elapsed recorded;

=== OPS-50
Closed 2026-09-21: (a) a B writer failure now raises in both examples with the exit-status control asserted, and (b) the post4 step-5 probe runs on adios2.bindings and prints PROBE_RESULT PASS after the licensed 0.11-image re-record and the two-sided `STEP4_SEP` ruling. The B-median nondeterminism known-issues entry stays open.
Done-when: (a) a `B` writer failure **raises** in both examples (the `A` writer may stay tolerant — no gate reads it); asserted control: a harness window that runs `mag:1` with `dolfinx.io.VTXWriter` monkeypatched to raise exits **non-zero**, and the same wrapper on the example file taken from the pinned pre-change commit exits **0** (the defect); the unpatched `mag:1` and `mag:2` reproduce `OPS-48`'s printed digits (`relative difference = 0.000e+00`, control `1.000e+00`) at `-n 2`; `PARAVIEW_GUIDE.md:166–170` reworded (both failure paths now exit non-zero).

=== OPS-51
Closed 2026-09-18: the six hardcoded rank-width and frequency log labels now format the real communicator size and frequency, with a unit guard counting 6 to 0 against the pre-change commit and the step-3a flag-off ladder control green at -n 8. No band, record or assertion was touched.
Done-when: the six literals replaced by the communicator's size and the rung's actual frequency; a `tests/unit` guard (real build, `-n 2`) asserting the count of the two literal patterns in `tests/validation/*.py` is **0**, with the asserted control that the same count on the files at the pinned pre-change commit is **6** (count identity, in = out); the ladder's step-3a flag-off control re-run at `-n 8` (its record width — `20260914T093325Z_ANS-4-step3a-w1.log`'s command verbatim) green with the 2a″ record assert **executed** at `STEP3A_RECORD_RTOL` (1e-6, unchanged) and the rung header reading `-n 8`;

=== OPS-52
Closed 2026-09-19: the WF-7 cost probe's phase-4 label is now a pure function that says flag-off control only at degree 1, the bracket dicts are re-dated from measured points, and the -n 8 flag-off control window passed with the S_driven record asserted.
Done-when: the phase-4 label is a pure function that says `flag-off control` only at degree 1 and otherwise prints the relative move beside the degree-1 record, asserted by a `tests/unit` test (including 5.50 % on the 09-19 values) with a source-text count control 1 → 0 against the pinned pre-change sha; the two non-gating bracket dicts re-dated from the measured points with their logs; the `-n 8` flag-off control window green with the probe's own executed `S_driven` assert against an untouched `S_DRIVEN_STEP0_RECORD`; the 2026-09-19 known-issues entry retired in the same commit.

=== OPS-53
Closed 2026-09-19: `write_setup_figure` now raises a ValueError on every rank for titles over `MAX_TITLE_CHARS` = 130, with the unit test, call-site scan and pre-change control green and the census unchanged. The clipped-legends half of the known-issues entry stays open for other items.
Done-when: a `MAX_TITLE_CHARS` constant set by measurement (bracket [126 whole, 150 clipped] ⇒ 130) and a `ValueError` raised before plotting on every rank; a `tests/unit` test asserting 131 raises / 130 does not, an `ast` scan asserting every constant `title=` at the 18 `write_setup_figure(` call sites under `examples/` is inside the limit, and the absent-guard control against the pinned pre-change sha; setup-figure census unchanged (`ok=18 missing=36 broken=0`) and one flagged `mesh:5` re-run digit-identical; the 2026-09-19 known-issues entry retired in the same commit.

=== OPS-55
Closed 2026-09-21: `evaluate_vector_field_parallel` now folds a point digest into its allgather and raises on every rank when equal-count point lists differ, with the pre-change control showing the silent defect and the caller subset unchanged. No current caller was passing rank-local points.
Done-when: in `tests/post/test_evaluation_collective_guard.py`, at `-n 1` and `-n 2`: equal counts with coordinates shifted by `0.01 * rank` inside the unit cube raise on every rank (allgathered flags `== [True] * size`), the message naming each rank's first differing row and deviation, and at `-n 1` the same construction does **not** raise (negative control); a pre-change control — the function at the pinned pre-change sha, `git show` into a scratch module under `/logs/` as `OPS-53` did — returns an all-valid mask on that input **and** values deviating from the closed form `(x, 2y, 3z)` at rank 0's points by the amount the shift predicts (quantitative: the defect was silent, and by how much);

=== OPS-56
Closed 2026-09-21: the sweep planner now decides endpoint inclusion and de-duplication relative to the step rather than the carrier frequency and rejects non-finite inputs, with the pre-change control reproducing all three failures. An auditor DEMOTE at the 2026-09-23 review was overruled and the closure stands.
Done-when: unit tests asserting exact grids for the three reproduced cases (endpoint present, length as computed by hand), a nondivisible span, the float-noise case `(0.1, 0.3, 0.1)` → 3 points, an on-grid endpoint bitwise equal to `stop_hz`, a coarse/refined near-duplicate merging to one point, and the intended `ValueError` for each non-finite input; a pre-change control showing the three reproduced cases fail on the pinned pre-change sha; the two existing planner tests green unchanged.

=== OPS-57
Closed 2026-09-21: the three determinant singularity tests became LinAlgError-to-ValueError wraps, `s_to_z` and `z_to_s` refuse non-finite input, and `sparameters_from_impedance` delegates to `z_to_s`. The unit round trip and the PORT gate re-runs at -n 2 are green with no tolerance touched.
Done-when: unit tests — round trip `z_to_s(s_to_z(S))` to ≤ 1e-12 on a 16-port matrix scaled so its determinant underflows (the pre-change control: the pinned pre-change sha raises on it); an exactly singular `I − S` still raises `ValueError` with the same message at all three sites; `sparameters_from_impedance` and `z_to_s` agree to ≤ 1e-14 on a random passive Z.

=== OPS-58
Closed 2026-09-22: the private-leak checker now scans commit messages, fails closed on inspection errors, reports SKIPPED distinct from clean, and has an exit-code test (24/24) plus a narrow shared-input allowlist for the one false positive. Audited PASS 2026-09-23; the operator still needs to re-run `install_git_hooks.sh` once.
Done-when: that test asserts exit codes exactly — planted figure staged → 1; planted figure in the message → 1; a truncated spelling of it → 1; clean → 0 with `clean:`; no reference data → 0 with `SKIPPED`; broken git → 2 — plus the pre-change control (the pinned pre-change script passes the message plant and the broken-git case); the CI `private-leak` job runs the same test, which becomes the one part of that job that tests something on a runner; `--audit` on this clone still reports clean (the count of figures printed, never a figure).

=== OPS-59
Closed 2026-09-21: (b) replaced the literal call-site count with the call-sites equal setup-PNGs equal census-ok identity, and (a) now admits setup figures by rule plus a declared list after the ans6 `metrics.json` was declared on the review's ruling. Both modules are green at -n 2 with the artifact identity reading 28 and the negative controls passing.
Done-when: (a) the artifact pin admits exactly the files matching `examples/*/figures/*_setup.png` that have a same-stem example script, by rule rather than by list, and nothing else — asserted negative control: a non-setup PNG and an orphan `*_setup.png` injected into the listed set both still fail the pin; count identity extras-before = committed setup PNGs, extras-after = 0.

=== OPS-60
Closed 2026-09-23: `scikit-rf` 2.1.0 is pinned in the image and its version test is green on the ordinary service at -n 2, unblocking `PORT-23`. An auditor DEMOTE on the tier word was overruled at the 2026-09-25 review and the tier relabelled standard.

=== OPS-61
Steps 0, 1a and 1b are done, and the first half of step 1c (item 159, 2026-10-01) shows ours at degree 1 agreeing with Palace on the exact f = 0.5 strip mesh to 0.0000 % in every class, so the step-1a gaps were port width. Open: item 166 (ours at degree 2), while the (a)/(b)/(c) reading against the private column is the weekly's.
Done-when: Palace `-np 8` solves `f05x` at 64 MHz: 117 206 elements, Rs 43.75 Ω/sq, port i = N_i, `+Z`, 50 Ω. Its 4×4 meets `RECIPROCITY_BAND`, σ_max and every column power ≤ 1 + `PASSIVITY_SIGMA_TOLERANCE`, and every S-class spread ≤ `ADJACENT_SPREAD_BAND` with pooled ≥ `POOLED_SEPARATION_FLOOR` × the worst — each asserted with the band imported.

=== OPS-62
Closed 2026-09-22 out of band: the `ANS-4` ladder's refinement guard now compares estimated unknowns instead of cells, strictly stronger, with pure-function controls green at -n 2. Audited PASS 2026-09-23; step 3e was ruled not re-run and 3f ran green under the new guard.

=== OPS-64
Closed 2026-09-28: `xl-run.sh` now records container rank counts before and after each window with an orphan verdict, and stages declared `XL_OUTPUTS` or refuses bad ones before the window, with the extended isolation test at 40/40 PASS. The tier was relabelled standard by the 2026-09-28 review.
Done-when: `xl-run.sh` prints the container-namespace rank count before and after every window into the log it commits, independent of the queued command; and either stages the output paths a queue file declares or refuses *before* the window starts when a declared output is neither ignored nor declared — asserted by extending `test_xl_window_isolation.sh` with one positive and one negative control each, run through the harness, elapsed recorded.

=== OPS-65
Closed 2026-09-27: the housekeeping citing set now includes chunk files, review journals, XL ledgers and xl-pending, 210 deleted logs were restored, and cited-but-missing fell from 284 to 74, all unrecoverable and listed. Audited PASS at the 2026-09-28 review.
Done-when: a unit test asserts that `cited_basenames()` contains a basename cited only from a chunk file and one cited only from the `xl` ledger; its control, pinned to the pre-change commit (`<sha>^`), is red; the dry-run census reads cited-but-missing = 0 (bar the unrecoverable, listed in `docs/testing/housekeeping.md`); elapsed recorded.
=== OPS-63
Done 2026-09-25 with audit PASS 2026-09-26: the `ans:2` example now asserts every quantity it labels ASSERTED, importing each band from the gate modules, and the AED figure is replaced by a pointer to the private adjudication file. The knob-route run reproduced the tracked metrics within the identity tolerance, the phrase census reads 0 after the change, and the `OPS-58` checker's failure to flag the figure is recorded as a finding.
Done-when: (i) the knob route at the tracked rung (`FEM_EM_ANS2_PHANTOM_RESOLUTION=0.0025`, `-n 4`, the route step 4a proved digit-identical to the unset run) runs green with the new asserts, and its `metrics_h0.0025.json` reproduces the tracked `metrics.json` on every primary leaf within the imported `EXACT_IDENTITY_RTOL` (step 4a measured worst 8.19e-11) — no number moved;

=== ARCH-1
Not started: embargoed until the weekly of 2026-10-17 by operator directive, so no daily review may queue it and earlier weeklies may only annotate the row. It waits on `OPS-61` step 1, `PORT-21` and one F-human physics step so the solver-agnostic interface is designed against two real backends.
Done-when: digit-identical reproduction of every imported band and record on the DolfinX backend before any second backend is wired — a refactor that moves a record is a defect.

=== GEO-5
Done 2026-09-30 (§9 item 103) with audit PASS: the phantom volume deficit falls at LSQ rate 1.9707 inside [1.7, 2.3], the equivalent-edge and partition gates hold on all four rungs, and the coil-slope control misses as predicted. Only the phantom leg is gated (coil and air stay on `GEO-17`'s sign gate), and a finer rung would need the a_eq/h_p band re-examined first.
Done-when: every mesh's tags partition it to 1e-9; tag 3's volume never exceeds πr²H and its relative deficit falls strictly, with LSQ rate against the requested size in [1.7, 2.3]; tag 3's equivalent edge (6√2·V₃/N₃)^{1/3} lies in [0.8, 1.5]·h_p at every rung; the coil's deficit regressed on the same ladder misses the rate band; `test_mesh_tag_integrity.py` is green in the same window. Harness, real, `-n 2`, elapsed recorded ⇒ ✅.

=== GEO-17
Done: step 1 (2026-08-20, audited compliant) made the region sizes actually apply, raising coil meshed/CAD to 0.8356. The `mesh:5` inverted control that lost separation on 0.11 was re-chosen on 2026-08-25 at `CONTROL_RESOLUTION` 0.018 with margins above the 0.05 separation, so the ✅ stands.

=== GEO-33
Step 0a (cracked box lips decoupled and tagged) is done as of 2026-10-02 (§9 item 149), and step 0b, the two-sided Leontovich lip term against the Pozar half-cavity Q ratio, is next. After that the 32-leg cost probe (step 1) and a single priced degree-1 strip-sheet drive remain, sequenced after the F-human degree-order question.
Done-when: 0a + 0b green ⇒ 🟡 (the lip term exists and is gated); ✅ additionally needs step 1 (the 32-leg cost probe) and the strip-sheet rung's single priced drive.

=== GEO-34
Steps 1a (STEP import with tag census), 1b (through the congruent-sheet config, 2026-09-26) and 2 (the preflight writer, 2026-09-30) are done, and §9 item 119 then wrote the import's own wall group so the wall re-enters step 2's identity. Steps 1c (independent kernel) and 3 remain open.

=== TH-5
Steps 1a (first-order ABC reflection on the TEM line falling at LSQ rate 2.0021 by volume projection), 1b-0 and 1b (the dipole on a curved ball boundary) all landed 2026-09-30. ✅ still needs step 2 and the §2.1 line, and the weekly ratifies the sphere re-pose.

=== TH-15
Steps 1, 2a, 2d, 3a, 3b, 3 and 3c are done: the birdcage 4×4 as a PEC hole passes the port gates at 10 / 64 / 128 MHz (2026-09-13) and its terminal-power excess is attributed to the sheets' Cauchy–Schwarz deficit. The row stays open on step 2's unitarity gate, whose `src/` replacement for the point-sampled path voltage is specified but unwritten and awaits a weekly re-scope.

=== TH-14
Done: the Leontovich surface-impedance boundary passed the lossy-cavity Q_c anchor, the copper F-small birdcage port gates at 10 / 64 / 128 MHz, and on 2026-09-14 the Dodd–Deeds copper-slab floor step that re-closed the row after the audit demotion. The §2.1 conductor-model line is written.

=== TH-16
Done 2026-09-30 (§9 item 93) with audit PASS: the quarter birdcage's three solves reproduce the 4×4 of its own mirror-built full mesh to 1.520e-12, with the PMC-for-PEC negative control asserted. Scope is one fixture at 10 MHz and degree 1, and the re-anchored gate is for the weekly to ratify.

=== TH-17
Done 2026-10-01 (§9 item 131): gates (c) and (d′) passed, completing the amended Done-when together with (a), (b) and (e), and the 10:30 review audited it PASS. Gate (a) is closed as re-posed at `C_mode` rather than `C_tuned`, a reading the weekly confirms.
Done-when: (a)–(e) through the harness at `-n 4`, the fixed-point shift recorded (< 0.1 % or the band re-read), elapsed recorded ⇒ **`TH-17` → ✅ and Phase 6's first physics target ticked.**

=== MAT-1
Done 2026-09-29 (§9 item 83, attempt 2): all three presets match the lossy-sphere closed form at 64 MHz, and the phantom path assembles exactly the same material arrays and operator as the `material_map` path. The original solve-level 1e-12 premise was re-anchored upstream of the solve, so no solver-determinism claim is made.
Done-when: On `MAT-4` step 1's lossy sphere (h = R/10, 64 MHz), each preset passed as `phantom_material` gives mean SAR and Im/Re E_z within 10 % of `uniform_sphere_sar_closed_form`

=== MAT-5
Done: step 1 (the Picard fixed point on the Robin ball) and step 2 (per-cell σ(T) inside the EM solve on `MAT-4`'s sphere) both passed 2026-09-30, audit PASS. The (D′) transfer gate's planted failures were measured the same day (§9 item 109), closing the audit caveat.
Done-when: ✅ = steps 1 and 2.

=== MAT-9
Done 2026-09-29 (§9 item 76): the steady and transient conduction steps are green and the held root anchor passes in its scale-relative form. The weekly's lane-3 ratification is separate.
Done-when: ✅ = steps 1 and 2 green.

=== MAT-11
Done 2026-09-30 (§9 item 97) with audit PASS: the perfused Robin ball meets its closed-form, rate and balance gates, both controls miss, and the rule (c) re-runs reproduce. A complex-build leg (§9 item 120) reproduced every digit; the weekly's ratification of the row stays open.
Done-when: (i) P1 interior relative L2 on θ ≤ 1 % at the 0.01 m rung (`MAT-9` (i), imported); (ii) P1 whole-ball LSQ rate ≥ 1.7 on the 0.02 / 0.01 / 0.005 m ladder; (iii) ∫q dV = ∮hθ dS + ∫wθ dV to ≤ 1e-9 on every rung; the perfusion-dropped and sign-flipped controls each miss (i)'s band; `MAT-9`'s, `MAT-5`'s and `MAT-10`'s modules re-run green at the default w = 0 (rule (c)) ⇒ ✅.

=== MAT-12
Done 2026-09-30 (§9 item 100, first try) with audit PASS: the complex-build transient driven by FEM SAR meets its closed-form, calorimetric-slope, balance and steady-limit gates, and both controls miss. The weekly's ratification of the row stands open.
Done-when: (E) per-step balance ≤ 1e-9 on every run; (∞) backward Euler 20 × 5 R²ρc/k equals `solve_steady_conduction` with q_h on the same submesh to ≤ 1e-8; the Re-E-only and quasi-static controls each miss (A_t)'s band; `MAT-9`'s modules re-run green in the real build (rule (c)) ⇒ ✅.

=== MAT-13
Done 2026-09-30 (§9 item 108), audited PASS at the 18:00 review: the transient Pennes θ-scheme meets its closed-form, CN rate, balance, steady-limit and w → 0 gates with all controls missing. The complex-build leg (§9 item 120) reproduced every digit, and the weekly ratifies the lane-3 row.
Done-when: the perfusion-dropped and sign-flipped controls each miss (P_w)'s band and the perfused-term-dropped balance misses 1e-3; `MAT-9`'s modules (real build) and `MAT-12`'s (complex build) re-run green on their recorded digits (rule (c)) ⇒ ✅.

=== MAT-14
Done: the complex-build SAR-driven perfused transient was green on every gate with each control missing its band, audited PASS at the 2026-10-01 03:00 review. The row was ratified by the 2026-09-30 weekly.

=== MAT-15
Step 1 is green and step 2(a)'s numpy series is recorded; step 2(b)'s FEM re-run as §9 item 160 is blocked on a bug in the test's (∞) check, with the steady and transient physics otherwise green.
Done-when: step 1 asserted green ⇒ 🟡; step 2 ⇒ ✅.
=== IMP-2
Both steps went green on 2026-10-01: step 0's numpy closed form, then step 1's FEM at aspect 2 and 5 inside the imported bands with the transverse control missing as predicted, audited PASS at the 18:00 review. Quasi-static, PEC, imposed trace only; no rod, no coil and no 1 g SAR claim.
Done-when: step 0 green ⇒ 🟡; step 1 at both aspects ⇒ ✅.

=== IMP-1
Closed 2026-10-01 on the bands-only closure ratified by the 2026-09-30 weekly and audited PASS; ρ is kept as a record with no convergence claim. Step 3 records follow: 3(a) done, 3(b) parked and read as "the floor is neither comparand model error nor box trace alone; discretisation not separated", with an a/6 degree 3 discriminator queued as §9 item 165.
Done-when: ✅ on (i)–(iii) at three rungs, (M), (iv′) and the control when §9 item 123 lands them green on `main`

=== POST-2
Closed 2026-09-30 and audited PASS: the ratio and span diagnostics read the lossy plane-wave closed form on `TH-6` within 5 % with the σ-blind control missing, and the cross-unit balance metric was replaced the same day by a reference-speed form. No energy quantity is gated and detecting a wrong solve is not claimed.
Done-when: On `TH-6`'s plane wave (0.1 m, 127.74 MHz, σ 0.7, εᵣ 78, 24³, degree 1), with B from `magnetic_flux_density_from_e` and both fields sampled over every owned cell by `compute_tagged_vector_magnitude_stats`: `e_to_b_mean_ratio` is within 5 % of ω/\|k\|, and `e_span_ratio` / `b_span_ratio` are within 5 % of 1 − e^{−αΔx_c}.

=== POST-8
Step 1 done 2026-09-27: the raw N1curl field in each named point's single incident cell is width-invariant, so the N1curl to CG1 nodal interpolation carries the partition dependence. Open: the weekly must choose a width-invariant sampling rule, adopt it at every census site, and re-issue `ans:2`'s Row 1.
Done-when: a sampling rule whose output is width-invariant to 1e-9 relative, asserted at `-n 2` and `-n 4` on `ans:2`'s points, and adopted by every point-sampling site the census found.

=== POST-9
Future work for Phase 9, not dated: the weekly places it, and it is not on §9 rule (5)'s lanes unless the operator releases it. Steps 0 to 4 are written but none has started.

=== PORT-5
Closed 2026-09-30 and audited PASS at the 10:30 review: through the real assembler, synthesised port states reproduce known S matrices within 1e-13 and the sanity metrics return their closed forms, with the transposed-comparand control missing by exactly 1. The frequency-sweep ambition has no landed code and stays outside the closure.
Done-when: (a) Through `run_n_port_sparameter_sweep`, port states synthesised from a known S₀ (a lossless reciprocal line, an ideal 3-port circulator, a seeded non-symmetric 4-port) assemble to S₀ within 1e-13, and the circulator misses its own transpose by exactly 1.

=== PORT-6
Closed 2026-09-30: the planner equals the integer-lattice closed form bitwise on the seeded 500-case family and the frozen `c1d167d` planner reproduces the `OPS-56` failures; the validity floor was fixed the same day, raising for steps under 16 ulp of the carrier. Planning only: orchestration of solves over a plan is outside the closure.
Done-when: `plan_frequency_sweep` equals the integer-lattice closed form bitwise on a seeded 500-case family (integer-Hz carriers 1e6–3e8 Hz, steps 100 Hz–1 MHz, on- and off-grid endpoints, 0–3 refinement windows); `z0_ohm` / `port_order` pass through exactly; the frozen `c1d167d` planner reproduces the three `OPS-56` failures (`20260921T112323Z_OPS-56.log:34–36`); the existing 19 tests are green unchanged.

=== PORT-14
Closed 2026-09-13 at step 3: the κ-derived width route is registered at 64 MHz with corrected residuals well under 1e-3 and the `PORT-9` gate module green. κ is carried as the sheet's named systematic, and the 128 MHz pair stays a printed out-of-sample reading.

=== PORT-20
All three steps executed 2026-09-20: the current-drive route now reports S as `z_to_s(Z)`, the closed-form identity and both negative controls are asserted, and the S-derived records are re-recorded at unmoved bands. Left over: the `ans:3` generated outputs still carry the pre-correction S table, filed in known-issues with its unblock condition.
Done-when: steps 1–3 executed by the implementer, the closed-form identity and both negative controls asserted, elapsed recorded, the known-issues entry retired, §2.1's `PORT-1` line gains "S on the current route is `z_to_s(Z)` since `PORT-20`".

=== PORT-21
Steps 1, 3, 4 and 5 have landed (step 1 closed 2026-09-21; steps 4 and 5 on 2026-10-01 as records with no glyph move), and step 2 is the weekly's private ruling. The row closes only when T6 is decided either way and §2.2 is written; no absolute S₁₁ or Z_in at 10 or 64 MHz is licensed.
Done-when: `PORT-21` → ✅ when T6 is decided either way (F12 for clause (i), F9 for clause (ii)) and §2.2 is written.

=== PORT-22
Step 1 closed 2026-09-20: the driven sweep across 44–84 MHz matches the circuit-layer prediction at every grid point, finds one sign change of Im Z_in bracketing 64 MHz, and the half-C control misses by far more than the band. Re-run green 2026-10-01 through the in-run κ option under `PORT-21` step 5; still one fixture at degree 1 with no absolute input-match claim.
Done-when: (i)–(iii) and the control asserted by the implementer, elapsed recorded; §2.1 gains "the F-small tuned resonance is located by a driven sweep consistent with the circuit layer across 44–84 MHz" — still one fixture, degree 1, series resonance, no absolute `S₁₁` (known-issues 2026-09-19), no mode-spectrum claim (`TH-17`'s).

=== PORT-23
Step 1 done 2026-09-23: the `scikit-rf` termination reproduces the raw (C2) reduction to machine precision on the stored 64 and 10 MHz matrices, with the swapped-port control asserted. Open for step 2, re-posed 2026-10-01 as reproducing Route B on `ans:7`'s stored grid matrices (§9 item 155).
Done-when: on `ans:7`'s 51 stored grid 4×4 S (`logs/ans7/B1.json` / `B2.json`, Cu and PEC), `terminate` (P2–P4 in `C_F`) reproduces Route B's `Z_bare` (`route_b_z` on `z_from_s`) to `S_VS_Z_ROUTE_RTOL` at every point; `parallel_rows` on the skrf-side `Z_in` gives f₀, R_in, Q_loaded within 1e-9 of the raw route and within `test_ans7_sparse_route_b.BANDS` of `metrics.json`'s rows; harness, complex, `-n 2` ⇒ `PORT-23` ✅.

=== WF-3
Closed 2026-09-30 and audited PASS at the 10:30 review: the eight quick-look numbers parsed from the written markdown match the lossy plane-wave closed form within 5 % and the JSON summary, with the σ-blind and swap controls missing. The balance line, the status line on physical fields and SAR are not claimed.
Done-when: the eight numbers parsed from the written markdown (\|E\| and \|B\| min / max / mean, the two shape ratios) are each within 5 % of the closed form on the owned centroids (e^{−αx_c}; × \|k\|/ω; ω/\|k\|); each equals its JSON `summary` value to 5e-7 relative; the JSON `summary` equals the metrics dict exactly; coverage reads 82 944/82 944/82 944 with no drops and no report warning flag; the σ-blind comparand misses by > 5 %; `POST-2`'s test is green in the same window.

=== WF-5
All three steps were asserted or recorded on 2026-10-01 and the Done-when is met, audited PASS at the 18:00 review with caveats (a thin A3 margin, the PEC control missing its tighter prediction). No h-convergence, absolute-Q, AED or tuning claim; the weekly still reconciles §10 T15 with the printed-only Δf₀.
Done-when: step 1 asserted ⇒ 🟡; steps 1–2 asserted and 3 recorded, elapsed recorded ⇒ ✅.

=== WF-7
Closed at step 3 on 2026-09-25 (audited PASS 2026-09-26): the 10 g SAR hotspot is located at human scale with the coverage, C16 and mirror identities green on one phantom at degree 1, 64 MHz. The degree-2 full 32×32 ran on 2026-09-26 as a record, and no compliance or absolute SAR claim is made.
Done-when: (ii) the **10 g hotspot** — the argmax of `mass_averaged_sar` over a fixed candidate lattice — and its fifteen **C16 images** within the imported `C4_COVARIANCE_BAND` (5 %); (iii) the **mirror identity** — the cw drive's 10 g SAR at the mirrored hotspot vs the ccw drive's at the hotspot — within the same band; (iv) every ball's containment in the phantom.

=== WF-10
Step 1 closed 2026-09-27 and step 1c was recorded as mixed, but step 1d W1 went red on 2026-10-01 (identity above `REDUCTION_BAND`, feed-ring k = ±1 share below 0.9) and is blocked, its code parked on an attempt branch. The reading is that the feed geometry, not the feed loading, sets the share; step 2 (tuned B₁⁺ and the 10 g SAR hotspot, on `WF-12` step 3's copper matrix) waits on step 1d.
Done-when: steps 1, 1d and 2 asserted, elapsed recorded ⇒ ✅.

=== WF-12
Step 1 (H12a) went green 2026-10-01: CAD identities and the one Cu drive's exact power partition asserted, priced at 369 387 cells and 9.6 GiB at `-n 8`, with the wall-area fixture band awaiting the 10-03 weekly's ratification. Steps 2 (the copper 32×32 gates) and 3 (C_mode re-derived on the copper matrix) are open.
Done-when: steps 1–3 ⇒ ✅; *(amended 2026-10-01 03:00 review, rule (j))* step 1's CAD identities (sheet areas, the box-minus-coil partition, cavity-wall area) and the one Cu drive's exact power partition asserted in one logged window, with cells, unknowns, time and memory recorded ⇒ 🟡 (§9 item 129).
=== EX-55
Closed 2026-09-14 and audited PASS: the 32-port ccw quadrature drive on the 16-leg birdcage asserts the C16 spread, the mirror identity and the exact power identity at their imported bands, with the census and the rule (a) gate re-run green. Nothing is open; no homogeneity, absolute or Larmor claim is made.
Done-when: the three imported identities asserted, artifact named, census `exit != 1`, elapsed recorded.

=== EX-57
Closed 2026-09-28: the setup-figure census exits 0 on `main` with every example guide carrying a setup figure and docrefs clean. The re-render of the earlier figures with clipped legends stays open for a review to queue.
Done-when: closes when the census exits 0 on `main`

=== EX-61
Done 2026-09-21 and audited PASS: `ports:19` shows the tuned birdcage's resonance curve with all three imported anchors green, f0, R_in and loaded Q printed but never gated, and the rule (a) lift re-gated green. Nothing is open; no absolute `S₁₁` or `Z_in` claim is made.
Done-when: the three imported assertions green through `./run_examples.sh -e ports:19`, artifacts named in the guide by full filename, both censuses before / after, docrefs `exit != 1`, elapsed recorded.

=== EX-62
Done 2026-09-26 and audited PASS: `mri:4` reproduces the `WF-7` human-scale 10 g SAR hotspot gate (coverage identity, C16 and mirror spreads, contained averaging balls, negative control) at `-n 8`, with the rule (a) lift re-run green and censuses matched. Nothing is open; no absolute SAR or compliance claim is made.

=== EX-64
Done 2026-09-30 and audited PASS: `th:11` solves the PEC and `NATURAL` (PMC) walls on one box, every imported anchor and control green on the gate's digits, rule (a) re-gate green and censuses clean. Nothing is open.
Done-when: asserts green through `./run_examples.sh -e th:11 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== EX-65
Done 2026-09-30 and audited PASS: `mat:3` runs the three gelled-saline presets on the lossy sphere, SAR and Im/Re inside `MAT4_BAND` against the closed form and every cross-preset control missing, with the re-gate and censuses green. Nothing is open.
Done-when: asserts green through `./run_examples.sh -e mat:3 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== EX-66
Done 2026-09-30 and audited PASS: `th:12` cuts the driven box on PMC and PEC symmetry planes and reproduces the full model's port records within the imported band, with the pinned-cut control missing, the gate re-run green and censuses clean. Nothing is open.
Done-when: asserts green through `./run_examples.sh -e th:12 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== EX-67
Done 2026-09-30 and audited PASS: `mat:4` heats the gel ball from SAR, with the steady, transient and calorimetric-slope readings equal to the `MAT-10` and `MAT-12` records under their imported bands and both controls missing. Nothing is open beyond the weekly's lane-3 ratification of the demonstrated rows.
Done-when: asserts green through `./run_examples.sh -e mat:4 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== EX-68
Done 2026-09-30 and audited PASS: `mat:5` runs its own σ(T) fixed-point loop over the `MAT-5` gate helpers, all eight iterates equal to the gate record and the frozen-σ control missing, with no gate edit. Nothing is open.
Done-when: asserts green through `./run_examples.sh -e mat:5 -n 2 -t 600`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== EX-69
Done 2026-09-30: `mat:6` runs green with the steady and transient Pennes perfusion readings on the Robin ball inside their imported bands, the w = 0 controls missing, and census and docrefs clean. Nothing is open.
Done-when: steady (i) `<= INTERIOR_BAND`, (iii) `<= BALANCE_BAND`, transient (P_w) `<= INTERIOR_BAND`, (E) `<= BALANCE_BAND`, all imported; controls w = 0 steady / transient each miss `INTERIOR_BAND`; asserts green through `./run_examples.sh -e mat:6 -n 2`, artifacts by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== EX-70
Done 2026-10-01: `th:13` reproduces the `TH-18` layered-foil reflection record for every termination to the record's digit, with the bulk-term control missing as predicted, the rule (a) lift re-gated green and censuses clean. Nothing is open; no copper gate or oblique-incidence claim is made.
Done-when: asserts green through `./run_examples.sh -e th:13 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== EX-71
Done 2026-09-30: `th:14` shows the total-current divergence on the two-slab box at CG2 against its vacuous CG1 twin, matching the `POST-3` record to the digit with the control and the gate re-run green and censuses clean. Nothing is open.
Done-when: asserts green through `./run_examples.sh -e th:14 -n 2`, artifacts named by full filename, both censuses before / after (`examples` +1), docrefs `exit != 1`, elapsed recorded.

=== ANS-2
Adjudicated AGREE 2026-09-18 on coil-driven whole-phantom power and pointwise SAR, and the step 4 halved-phantom rung run 2026-09-25 was ruled branch (a): the residual is the feed, and the pointwise band tightens to the measured spread for the weekly to carry into the private file. Step 2 remains an optional record step; no compliance, homogeneity, B₁⁺, Larmor or quadrature claim is licensed.

=== ANS-3
Adjudicated 2026-09-19 by the weekly: AGREE on the primary row, the mutual coupling, while the S rows DISAGREE through a defect on our side, the current-drive route's S table, handed to `PORT-20`. The self-impedance rows are informative only, as pre-declared; nothing else is open.

=== ANS-4
The chunk's self-consistency gates are green on every rung and the 2026-09-19 weekly ruled the order-matched comparison AGREE on the couplings and at 128 MHz but DISAGREE on the self class at 10 and 64 MHz, tracked by `PORT-21`. Step 4 landed 2026-09-28, honouring the conductor-factor knob on the rung path and readying `xl` entry 20, while the degree-2 ladder follow-ups await the weekly's ruling.

=== ANS-5
Done 2026-08-31: every ANS spec and comparison table now pins our degree 1 to HFSS Zero Order with a First Order sensitivity column, per the 2026-08-30 weekly ruling, and both generators were re-run green. Nothing is open; step 3 was ruled not taken.

=== ANS-6
The runnable half, the degree and resolution knobs and the step 3 private-mode comparison writer (2026-09-28) have all landed, and the operator's AED numbers are in. The glyph stays 🟡 because the weekly's adjudication is still open.

=== ANS-7
Adjudicated 2026-09-30: AGREE on the mode, on Route A against Route B and on the copper contribution, but the absolute rows are INCONCLUSIVE at degree 1 on one rung. Step 3 completed 2026-10-01 as records (the 7-node sparse Route B licensed at degree 1 and degree 2, plus the degree-2 rows), and the row stays open until T11/T12 supply the order- and h-matched rung for the weekly's T13 re-read.
Done-when: The row stays 🟡 until §10 T11 / T12 supply the order- and h-matched rung and the weekly re-reads (T13).

=== ANS-8
The runnable half landed 2026-10-01 against the v1 rows and the operator's v1 AED column is in but not adjudicated. SPEC v2 makes the own-`C_mode` row 2b primary, and the pure-algebra v2 rows follow-up (§9 item 161) plus the copper and First Order AED columns are open.
Done-when: ✅ only after the AED half lands and the weekly adjudicates share and leak at a stated relative level.
