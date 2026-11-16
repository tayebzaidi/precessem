# OPS-78 — which identity gates compare two separate solves

*History file for the §7 `OPS-78` row (`OPS-46`). Opened by the 2026-11-02 10:30 daily review, lane 1, from `OPS-77`'s tail branch and known-issues' `OPS-43` observation.*

**The question.** `OPS-77` measured one gate on one fixture and found a tail: (S) reached 0.115 × `IDENTITY_BAND`, and `WF-16`'s four-port (I) landed at 1.004e-13 on another fixture, so the reading is not the rod's alone. The weekly decides tonight whether the gate's design changes, and it cannot decide that from one gate. The name `IDENTITY_BAND` appears on 147 lines in 27 files. Most uses compare a quantity with itself on a shared route, where the floor is 0.000e+00 (`OPS-77`'s control), and are not exposed. **Nobody has listed the ones that compare two separate factorisations.**

**Done-when.** A table of every use site of a band ≤ 1e-11 on a solve output, classified shared route / separate solve / not a solve output, with the count checked against `grep`; for each separate-solve gate, its largest recorded read and its headroom; the observation gains the row; no band moved.

### §9 item 761 — 2026-11-02 (lane 1, known-issues; `OPS-77`'s tail branch)

**Why.** A red that is noise costs a slot and teaches people to distrust reds. This gives the weekly the whole exposed set. **A band is never loosened to swallow this,** and this item proposes no band.

**Anchor.** (Census) `grep -rn` over `tests/`, `src/` and `examples/` for the names defined with a value ≤ 1e-11 that gate a solve output: `IDENTITY_BAND`, `ROUTE_IDENTITY_BAND`, `DECOMPOSITION_IDENTITY_BAND`, plus any further `*_BAND` constant of that size the first grep shows imported by a test that solves. **Every use site is classified, and the table's row count equals the grep's count for each name,** asserted by the reader script below. Mesh and geometry identities are classified "not a solve output" in one line and not itemised.

**Negative control.** `mri:38`'s sweep-vs-point (S) must land in "shared route", with its recorded 0.000e+00 (`20261102T060640Z_EX-170.log:5643`), and `mri:37`'s (S) in "separate solve", with `OPS-77`'s 1.147e-13. If the classification splits those two wrongly, the method is wrong: stop.

**How to classify.** Read the call: two factorisations in one window on the same operator are a separate solve; one solve read twice is a shared route. For each separate-solve gate, read the largest recorded value through `tests.log_records.read_log_text` or the `.gz`-safe route in `scripts/testing/read_s_identity_record.py`, extended with a pattern list if that is the shortest way. **Name the log and line of every value.** A gate with no printed read says so; it is not given a value.

**Execution.** Extend `scripts/testing/read_s_identity_record.py` or add a sibling, stdlib only. Table in this file, largest read first. Known-issues' `OPS-43` observation gains one row with the count and the top three gates. No test, band or example is edited. No solve.

**Tier, ranks, wall clock.** Smoke, one harness window ≈ 10 s, through `scripts/testing/run_and_log.sh OPS-78-census "timeout -k 30 30 python3 scripts/testing/<reader>.py"`. Most of the item is reading. **Class:** ops / docs. **Cost:** ≈ 5.5 measured / 25 predicted min. **Timebox:** minute 30 of the item; a partial table lands with the unread sites listed. **Depends on:** nothing.

**Traps already paid for.** **Never `gunzip` a record to get green.** A band name can be re-exported, so count definitions and uses apart. `tests/mesh/test_birdcage_leg_offset.py:101` defines `IDENTITY_BAND` at 1e-12, `tests/validation/test_port_circuit_skrf.py:83` defines another, and `tests/validation/test_geo38_f_clinical_build.py:39` defines a **different** one at 1e-9, so resolve each use to its import. `.eval` bytes depend on batch composition (item 537), a different cause: a gate of that kind is listed as such.

**Scope.** A census. **Not claimed:** a failure rate, a cause, a new band, or a design.

**Negative result.** A separate-solve gate with a recorded read above its band is a finding: name it, stop the census there, BLOCKED, and report. No `PREDICTED RED:` line.

**Status it can move.** `OPS-78` ✅; known-issues' `OPS-43` observation gains the exposed set.

**Outcome (2026-11-02T18:08Z): STOPPED at the negative control. The census is not complete, and `OPS-78` does not close.** Reader `20261102T180848Z_OPS-78-census.log`, `Status: 3` (`:95`), 8 s (`:96`), smoke, host `python3` through the harness as the brief writes it, no solve.

**The control split wrongly (`:85–86`).** By the brief's rule, "two factorisations in one window on the same operator are a separate solve", **`mri:38` lands in "separate solve"; the brief requires "shared route".** `mri:37` lands where required. The call: `examples/mri/38_implant_rod_on_a_sweep.py:100` runs `solve_sweep`, and `:107` runs a new `Study(...).solve(comm)` per frequency. `_study_solve_sweep` (`src/fem_em_solver/study.py:1092`, `:1114`) and `_study_solve_ports` (`:777`, `:799`) each call `geom.build(comm)` and then `run_n_port_sparameter_sweep` with no solver handed across, and the example's own docstring says each "builds its own mesh" (`:8`). `reuse_factorization` shares a factor among the drives of one sweep, not between the two calls. So the 0.000e+00 at `20261102T060640Z_EX-170.log:5643`, `:5646` is **not** one solve read twice: it is two factorisations, in one process, reading bit-identical, on all 6 recorded reads (`:91`). `OPS-77`'s line that this is "the floor when both sides share a route" rests on the same premise and was not checked against the call.

**Why it stops the census.** The method separates exposed gates by counting factorisations. On this reading that count does not separate the 0.000e+00 gate from the 1.147e-13 gate, since both have two. Whatever does separate them is not named here, and every classification would inherit that gap. **Not claimed:** a cause. Two facts for whoever re-poses the rule, neither tested here: `mri:37`'s two sides come from different builders (`Study` geometry against `_build_rod` plus `build_four_port_sweep`), and its sheet `h`, `w` and `area` are gated only to `IDENTITY_BAND` (`tests/package/test_wf15_study_implant_rod.py:298`), not to equality. Against that, `WF-17 s25`'s (P=Y) compares two `Study` routes and still reads non-zero.

**What landed before the stop.** Name counts only, no site classified. 817 modules parsed, 38 `*_BAND` names defined ≤ 1e-11, **336 use sites** resolved to a definition ≤ 1e-11, 21 unresolved, and every name's total equal to `grep -rnw`'s (`:36–74`). The largest are `S0_BAND` 50, `B1_BAND` 44, `IDENTITY_BAND` 39, `D_BAND` 38, `Y0_CODING_BAND` 31, `L_BAND` 26. Reads of the comparisons opened before the stop, **classification withheld**, none above its band (`:87–92`):

| comparison | largest read | log:line | × band | headroom | reads |
|---|---|---|---|---|---|
| `WF-17 s25` (P=Y), YAML vs block | 1.180e-13 | `20261101T134207Z_WF-17-s25-W1.log:10975` | 0.1180 × 1e-12 | 8.5 × | 2 |
| `WF-17 s25` (S), block vs direct | 1.172e-13 | `20261101T134207Z_WF-17-s25-W1.log:10975` | 0.1172 × 1e-12 | 8.5 × | 2 |
| `mri:37` / `WF-15 s10` (S), Study vs direct | 1.147e-13 | `20261102T121721Z_OPS-77-mri37-r2.log:5568` | 0.1147 × 1e-12 | 8.7 × | 5 |
| `WF-16` (I), Study vs direct | 2.038e-13 | `20261005T203015Z_WF-16-W2-deg2.log:3916` | 0.0002 × 1e-9 | 4906.8 × | 8 |
| `mri:38` / `WF-15 s11` (S), sweep vs point | 0.000e+00 | `20261102T060640Z_EX-170.log:5646` | 0 × 1e-12 | floor | 6, all 0 |
| `PORT-21` step 5, option vs hand-rolled | no printed read found | none | band 1e-12 | none | 0 |

**Unread sites.** All 336 are unclassified. Call sites read so far: `mri:37`, `mri:38`, `WF-15 s10`, `WF-15 s11`, `WF-16` tuned (`tests/package/test_wf16_study_tuned_f_human.py`), `WF-17 s25`, `tests/package/test_wf17_tuned_copper_tutorial.py`, `PORT-22`'s route identity. Everything else is unread. The 21 unresolved reads sit in `examples/materials/18_six_minute_window_coil_heating.py` (12) and `examples/mri/32_femem_run_ball_probes_and_duty_quadratic_files.py` (9), and the reader did not follow how those files bind the names.

**Disclosed.** (1) `WF-16` (I), the neighbour `OPS-77` set beside a 1e-12 gate, is asserted at `LEG_D0_REPRODUCTION_BAND` = 1e-9 (`tests/package/test_wf16_study_four_port.py:402`), so it is outside this census's ≤ 1e-11 limit. Its largest read is 2.038e-13, above the 1.004e-13 `OPS-77` quoted. (2) The brief's "147 lines in 27 files" counts the substring, which also catches `AREA_IDENTITY_BAND` and the like. Today that grep gives 156 lines in 27 files; the word grep gives 103. The 1e-9 definitions in `test_geo38_f_clinical_build.py` and `06_birdcage_leg_gaps_port_sheets.py` account for the reader's 35 out-of-scope lines of the name. (3) Bare literals such as `<= 1e-12`, and `*_RTOL` or `*_BAR` constants, are outside the brief's `*_BAND` scope and were not counted. (4) The reader's control test is a source check on the two `study.py` functions plus the example's lines, not a measurement of MUMPS factor calls. (5) The reader had unlogged host dry runs before its one harness window, and was last edited before that window. It exits 3 on a split control by design. **No band moved; no test, example or `src` edited; known-issues untouched,** because the row belongs to a finished census.

**For the review.** Re-pose the classification rule and the control together, then re-run. Known-issues' `OPS-77` row says the 0.000e+00 is "the floor when both sides share a route", and that wording needs the same correction.

**Status moved (12:00 CST slot).** None; §9 item 761 BLOCKED, unblock on a review's re-posed rule and control. Attempt 1. The reader is parked on `attempt/OPS-78-20261102T181007Z`. The slot checked the split against `38_implant_rod_on_a_sweep.py:100` and `:107` and the log's `:85–86`, and agrees with it. The `implementer` was told to add the known-issues row and declined, for the reason in (5) above; the slot leaves it that way, so the known-issues wording also stands until the review corrects it.

### Review ruling 2026-11-02 18:00 — the control was wrong, not the rule

**The brief's rule stands; the required label for `mri:38` flips.** The slot stopped correctly, and its reading of the call is confirmed here (`examples/mri/38_implant_rod_on_a_sweep.py:100`, `:107`). The error was the 10:30 review's: it took 0.000e+00 for proof of one solve read twice, and `OPS-77` had written the same premise.

**Why the rule stands.** The census exists to list what the drift can reach. Known-issues' `OPS-43` observation localises that drift after assembly, in the factor or solve, with the A / b fingerprints identical. So a gate is exposed by repeating a factorisation, whether or not a given pair of solves happened to drift. **A zero read is a read inside the exposed class, not a different class.** `OPS-43`'s own tally at `-n 2` was 2 drifted of 34 on another fixture, so most pairs there read bit-identical.

**What stays open, and is not claimed.** `mri:37` read non-zero 5 of 5 and `mri:38` zero 6 of 6, which is more of a split than an intermittent drift alone suggests. `mri:37`'s two sides come from different builders and `mri:38`'s from one, but `WF-17 s25` (P=Y), two `Study` routes, reads non-zero. So **the census records the builder as a column and draws no conclusion from it.** No cause is claimed, and this review does not assert that `mri:38` can drift.

**Done-when unchanged.** Only "How to classify" and the negative control change. `OPS-77`'s line and known-issues' row are corrected by this review.

### §9 item 767 — 2026-11-02 (lane 1, known-issues; `OPS-78`, attempt 2 of 4; re-poses item 761)

**Why first.** The 11-02 weekly decides the gate's design at 21:00 and asked for the exposed set. A partial table on `main` by then is worth more than a complete one after.

**Start from the parked reader:** `git checkout attempt/OPS-78-20261102T181007Z -- scripts/testing/census_identity_bands.py`. Code by path is allowed. The branch is deleted when this item lands.

**Anchor.** (Census) as item 761: every use site of a `*_BAND` name ≤ 1e-11 resolved, the count equal to `grep -rnw` per name, asserted by the reader. Item 761's window counted 336 sites; this review's unlogged re-run of the same parse reads 348 on today's tree, since new modules landed after it. **The reader's own count is the record.**

**How to classify (replaces item 761's paragraph).**
1. **Not a solve output, in bulk.** A site in a module with no solver call is classified mechanically, by a token list in the reader, and reported as one line per name. This review's unlogged pass found 45 such modules holding ≈ 184 sites, with the tokens `.solve(`, `solve_sweep`, `run_n_port_sparameter_sweep`, `build_four_port_sweep(`, `_sweep(`, `solve_time_harmonic`, `solve_transient`, `_solve_`, `solver.solve`, `.solve_ports`, `run_two_port`, `_run_drive`. Confirm the list against the code; it is a filter, not a proof.
2. **The rest by reading the call.** ≈ 164 sites in 45 modules. **Separate solve:** the two sides of the comparison come from two factorisations in one window, same operator or not. **Shared route:** one factorisation, whether one solved field read twice or several drives through one held factor. **Not a solve output:** a mesh, area or geometry identity, or a closed form with no solve, even in a solve-bearing module.
3. **For every separate-solve gate, two extra columns:** *builder* (`same` when both sides reach the mesh through one builder function, `different` otherwise) and *reads* (count, and how many exactly 0). A gate with no printed read says so and gets no value.

**Order, so a partial table is the useful one:** the port and S-matrix modules first (`tests/package/`, `test_port*`, `examples/mri/37` and `38`, `test_th16*`, `test_th20*`, `test_wf1*`, `test_geo38*`), then the thermal modules (`test_mat*`).

**Negative control (replaces item 761's).** All three required, or stop:
- `mri:38` / `WF-15 s11` (S) → **separate solve**, builder `same`, 6 reads, all 0 (`20261102T060640Z_EX-170.log:5643`, `:5646`).
- `mri:37` / `WF-15 s10` (S) → **separate solve**, builder `different`, largest 1.147e-13 (`20261102T121721Z_OPS-77-mri37-r2.log:5568`).
- `tests/validation/test_port_lumped_two_torus.py:709`, `v_lumped` against the independently assembled `v_sheet_average` on one solved field → **shared route**. This is the control the rule was missing: without a member, the class could be empty. **Read `lumped_run` first and confirm it is one solve;** if it is not, name another `*_BAND` site that is, and say so.

The reader's `control()` carries these three and exits 3 on a split, as before. Its source check is not a count of MUMPS factor calls, which stays disclosed.

**Execution.** Table in this file, separate-solve gates largest read first, with log and line per value. Known-issues' `OPS-43` observation gains one row: the counts per class and the top three. No test, band or example is edited. No solve.

**Tier, ranks, wall clock.** Smoke, one harness window ≈ 10 s, `scripts/testing/run_and_log.sh OPS-78-census "timeout -k 30 30 python3 scripts/testing/census_identity_bands.py"`. Most of the item is reading. **Class:** ops / docs. **Cost:** ≈ 5.5 measured / 30 predicted min. **Timebox:** minute 30 of the item. **A partial table lands at the timebox with the unread sites listed;** it counts as landed only if the three controls pass and each module it lists is complete. **Depends on:** nothing. **Family:** `OPS-78`, attempt 2 of 4.

**Traps already paid for.** All of item 761's. **Never `gunzip` a record to get green:** the `.gz`-safe route stays. A bare literal such as `<= 1e-12`, and `*_RTOL` / `*_BAR` constants, stay out of scope and stay disclosed as such. `WF-16` (I) gates at `LEG_D0_REPRODUCTION_BAND` 1e-9, outside the limit. Item 761's `PORT-21` step 5 read "no printed read found" with `20261001T171226Z_PORT-21-step5-W1.log.gz` and siblings present, so look at one of those before banking that; it may be the reader's pattern. The 21 unresolved reads sit in two example files (`mat:18`, `mri:32`), and the reader did not follow how they bind the names.

**Scope.** A census. **Not claimed:** a failure rate, a cause, a new band, a design, or a statement that a zero-read gate cannot drift.

**Negative result.** A separate-solve gate with a recorded read above its band is a finding: name it, stop, BLOCKED, report. A split control: stop. No `PREDICTED RED:` line.

**Status it can move.** `OPS-78` ✅ on a complete table; on a partial one it stays open, with the unread list; known-issues' `OPS-43` observation gains the exposed set.

**Outcome (2026-11-03T01:40Z, 19:30 CST slot): complete. The table is whole, 0 sites unread, `OPS-78` ✅.** Reader `20261103T013749Z_OPS-78-census.log`, `Status: 0` (`:141`), 8 s (`:142`), smoke, through the harness as the brief writes it, no solve, and the reader is as committed. Every one of the 38 names' totals equals `grep -rnw`'s (`:74`). **369 lines classified, 0 unread (`:113`):** 92 separate solve, 51 shared route, 217 not a solve output, 9 out of scope (`:75–78`). The 369 are the 348 resolved sites plus the 21 lines `ast` could not follow. The 9 out of scope are `mri:32`'s `P_BAND`, which `b.P_BAND` binds to `MAT-10`'s 1e-2. The other 12 are `mat:18`'s `D_BAND` and `L_BAND`, bound from `m17s3` at 1e-13 and 1e-12, and they are separate solve.

**All three controls land as required (`:114–117`).** `mri:38` → separate solve, builder same, 6 reads, all 0. `mri:37` → separate solve, builder different, 1.147e-13. `test_port_lumped_two_torus.py:709` → shared route: `lumped_run` has one `.solve(` call, and its docstring reads "One mesh, one solve, both routes read off the same field". No substitute control was needed.

**Separate-solve gates, largest read first (`:118–138`).** None is above its band. Reads are every printed value the reader's patterns found across the named logs, so "reads" is a count of prints, not of processes.

| gate | fixture | band | largest read | log:line | × band | headroom | builder | reads (exactly 0) |
|---|---|---|---:|---|---:|---:|---|---:|
| `MAT-21` steps 2, 4 (H0), duty vs Duhamel sum, per-step mean | gel box | 1e-12 | 3.801e-13 | `20261014T003128Z_MAT-21-step2-a2-W2-identities.log:123` | 0.3801 | 2.6 × | same | 3 (0) |
| `MAT-20` / `mat:18` (FEM-2d) | F-small coil-driven phantom | 1e-13 | 3.476e-14 | `20261012T170911Z_EX-110.log:1920` | 0.3476 | 2.9 × | same | 5 (0) |
| `MAT-20` / `mat:18` (FEM) | F-small coil-driven phantom | 1e-13 | 3.434e-14 | `20261011T172810Z_MAT-20-step2-fem-W1.log:1969` | 0.3434 | 2.9 × | same | 5 (0) |
| `MAT-21` steps 2, 4 (FEM-c), left / right of a one-sided bound | gel box | 1 | 3.011e-01 | `20261014T003128Z_MAT-21-step2-a2-W2-identities.log:130` | 0.3011 | 3.3 × | same | 2 (0) |
| `mat:18` (Cm) | F-small coil-driven phantom | 1e-13 | 1.306e-14 | `20261012T170911Z_EX-110.log:1935` | 0.1306 | 7.7 × | same | 2 (0) |
| `WF-17 s25` (P=Y), YAML vs block | rod | 1e-12 | 1.180e-13 | `20261101T134207Z_WF-17-s25-W1.log:10975` | 0.1180 | 8.5 × | same | 2 (0) |
| `WF-17 s25` (S), block vs direct | rod | 1e-12 | 1.172e-13 | `20261101T134207Z_WF-17-s25-W1.log:10975` | 0.1172 | 8.5 × | different | 2 (0) |
| `mri:37` / `WF-15 s10` (S), Study vs direct | rod | 1e-12 | 1.147e-13 | `20261102T121721Z_OPS-77-mri37-r2.log:5568` | 0.1147 | 8.7 × | different | 5 (0) |
| `MAT-21` steps 2, 4 (Δ) | gel box | 1e-12 | 7.868e-14 | `20261013T183523Z_MAT-21-step2-W2-identities.log:111` | 0.0787 | 12.7 × | same | 3 (0) |
| `MAT-20` / `mat:18` (H0) | F-small coil-driven phantom | 1e-12 | 5.908e-14 | `20261012T170416Z_EX-110-r2-mat20-regate.log:1968` | 0.0591 | 16.9 × | same | 5 (0) |
| `mat:15` (Hist) | F-small coil-driven phantom | 1e-13 | 5.035e-15 | `20261008T004928Z_EX-99-item306.log:1936` | 0.0503 | 19.9 × | same | 21 (0) |
| `MAT-21` step 13 (C100-h) | gel box | 1e-13 | 4.885e-15 | `20261020T141240Z_MAT-21-step13-W2-duty-dose-fem.log:89` | 0.0489 | 20.5 × | same | 4 (0) |
| `MAT-21` steps 2, 4, 11 (D′) | gel box | 1e-13 | 4.481e-15 | `20261014T184215Z_MAT-21-step4-W2-p2-identities.log:118` | 0.0448 | 22.3 × | same | 4 (0) |
| `MAT-17` step 3 / `mat:15` (Δ) | F-small coil-driven phantom | 1e-12 | 3.719e-14 | `20261007T183803Z_EX-99.log:1939` | 0.0372 | 26.9 × | same | 9 (0) |
| `MAT-21` steps 2, 4 (FEM′-2d), residual subtracted | gel box | 1e-13 | 2.759e-15 | `20261014T184215Z_MAT-21-step4-W2-p2-identities.log:128` | 0.0276 | 36.2 × | same | 2 (0) |
| `MAT-21` steps 2, 4 (FEM′), residual subtracted | gel box | 1e-13 | 2.579e-15 | `20261014T184215Z_MAT-21-step4-W2-p2-identities.log:128` | 0.0258 | 38.8 × | same | 2 (0) |
| `MAT-17` step 3 / `mat:15` (D′) | F-small coil-driven phantom | 1e-13 | 2.152e-15 | `20261007T184229Z_EX-99-control.log:1933` | 0.0215 | 46.5 × | same | 9 (0) |
| `mri:38` / `WF-15 s11` (S), sweep vs point | rod | 1e-12 | 0.000e+00 | `20261102T060640Z_EX-170.log:5646` | 0 | floor | same | 6 (6) |
| `PORT-8` anchor S′, flipped sheets vs DSD | F-small | 1e-12 | 0.000e+00 | `20260929T110818Z_PORT-8-W1.log:1903` | 0 | floor | same | 1 (1) |
| `PORT-8` anchor Z′, flipped sheets vs DZD | F-small | 1e-12 | 0.000e+00 | `20260929T110818Z_PORT-8-W1.log:1903` | 0 | floor | same | 1 (1) |
| `PORT-21` step 5, option vs hand-rolled | F-small | 1e-12 | no printed read | none | none | none | same | 0 |

**The other classes, one line each (`:79–103`).**
- **Separate solve, 92 sites:** the 21 rows' gates, plus 23 sites of misses asserted *above* a band (`control, a miss`). A miss has no headroom to lose, so it has no row.
- **Shared route, 51 sites:** (Lin) `Y_in` 18, all on one node pair that the fixture asserts as `["solve", "held"]`; `v_lumped` against the sheet average and siblings 10; reciprocity through one held factor 7; (X1) 4; (Im W) 3; S vs `z_to_s(Z)` 3; Y symmetry 2; (FEM-c-win) 2; (Y) 1; one printed route 1.
- **Not a solve output, 217 sites:** no solver token 149 in 34 modules; (S0) schedule census 18; recorded S through `skrf` algebra 8; synthetic networks 12; closed form and quadrature 8; sheet geometry and spec fields 9; `numpy` round-off 5; the rest 8.

**What the table says, with no conclusion drawn.** The two tightest gates by headroom are thermal, not port gates. Against a builder split, `mri:37` and `WF-17 s25` (S) are different and non-zero, but `WF-17 s25` (P=Y) is same and non-zero, and most "same" rows are non-zero, so builder does not sort the reads. `PORT-8` is a second all-zero gate beside `mri:38`, on one read. **Not claimed:** a failure rate, a cause, a new band, a design, or that a zero-read gate cannot drift.

**Disclosed.**
1. **A read above band that is not a stop.** `MAT-21` step 2's as-posed (FEM) read 1.464e-13 against `D_BAND` 1e-13 (`20261013T183523Z_MAT-21-step2-W2-identities.log:121`, (FEM-2d) `:123`). The brief's line "a recorded read above its band is a finding" does not reach it, because the 2026-10-13 18:00 ruling had already re-posed that assertion as (FEM′), (FEM-c) and printed-only (FEM-ctl), and no use site compares it to a band today. The slot decided that; the reader's `READS` carries only the live gates. The weekly should know it is the same form as `MAT-20`'s (FEM), which still gates at 1e-13 and reads 0.34 ×, on another fixture.
2. **Residual-subtracted forms.** (FEM′), (FEM′-2d) and (FEM-c) subtract or bound by the run's own residual r, so a drift between runs should largely cancel in them. They count as separate solve by the brief's literal rule, and 10 of the 50 `Duhamel` sites are these. Their low reads may not measure exposure.
3. **`PORT-21` step 5.** `20261001T171043Z_PORT-21-step5-W0.log.gz`, `20261001T171226Z_PORT-21-step5-W1.log.gz`, `20261001T171613Z_PORT-21-step5-W2.log.gz` and both `*_PORT-22.log.gz` were searched for `dS|`, `route identity`, `route_identity` and `hand-rolled`, and none carries the print, so "no printed read" is the logs' and not the reader's pattern. The print and its node sit behind `FEM_EM_PORT22_WIDTH_ROUTE`; whether any of those windows set it, or ran with `-s`, was not established, so the gate gets no value and it is not known here that it has ever run.
4. **The token filter is a filter.** It left 149 sites as N, and it missed real solves that run through helpers: `mat:15`, `mat:18`, `test_mat20_fem_history.py` and five `Y0_CODING_BAND` modules carry `RULES` entries. Eighteen no-token modules naming a band likely to gate a solve were read, not all 34 line by line.
5. **Thermal class by the rule.** The Duhamel gates repeat runs, so they are separate solve here. Whether their transient solver is the factorisation `OPS-43` localised was not read, and their reads are round-off-sized against bands set near round-off, so a tight headroom there may have another source.
6. **Held factors.** "One held factor" for reciprocity and Y symmetry rests on `reuse_factorization` defaulting on for the lumped-sheet route (`src/fem_em_solver/ports/sparameters.py:421`). (Lin) rests on the fixture's `["solve", "held"]` assert. Neither is a count of MUMPS factor calls, and the controls' source check is not one either.
7. **Scope.** Bare literals such as `<= 1e-12`, and `*_RTOL` / `*_BAR` constants, stay out, so `UNIT_WEIGHT_CONTROL_RTOL`, `REUSE_REPRODUCTION_RTOL` and `DISCRETE_IDENTITY_RTOL` are not here. `WF-16` (I) at 1e-9 is outside the limit. `RULES`' line numbers pin today's tree: if a listed module moves, a re-run can misfile a site, and the reader only flags one that no rule catches.
8. **The reader had unlogged host dry runs before its one window,** and was last edited before that window. `read_s_identity_record.py`'s docstring keeps the old "floor when both sides share a route" wording; the brief edits no other script, so it stands.

**Status moved (19:30 CST slot).** `OPS-78` ✅; the known-issues `OPS-43` observation gains the row; `attempt/OPS-78-20261102T181007Z` is deleted, as the brief says. **No band, test or example moved.**

### §9 item 815 — 2026-11-04 (ops / docs; `OPS-78` census refresh at HEAD; opened 10:30 daily review from the `EX-180` audit)

**Why now.** The census's last window is `20261103T013749Z_OPS-78-census.log` on `f038ae35` (2026-11-02 18:12). Since then 39 files under `tests/`, `examples/` and `src/` changed (+6 826 / −64). Most are new modules (I55 … I63, P63, P65, `mat:37`, U121 … U123), but **two files with line-pinned `RULES` entries were edited**: `examples/mri/38_implant_rod_on_a_sweep.py` (rule `(161, 163)`) and `tests/package/test_wf17_implant_rod_tutorial.py` (rules `(158,)`, `(264, 285)`, +146 lines). This is point 7's trap: "if a listed module moves, a re-run can misfile a site". The 10:30 review's `EX-180` auditor also read, by code only, that `mat:37` (`PRE_BAR` 1e-13, a `solve_transient` token, no `RULES` entry) would print `UNREAD`. The 11-02 weekly made this table "the reading of any future separate-solve red", so it is kept current at HEAD.

**Anchor.** The reader's own gates, unchanged: every `*_BAND` name ≤ 1e-11 resolved, each name's total equal to `grep -rnw` (exit 2 otherwise), **0 unread** (exit 4 otherwise), the three controls (exit 3 otherwise), and no read above its band. **Exit 0 on the closing window is the gate.**

**Execution.**
1. **W1, read-only.** `scripts/testing/run_and_log.sh OPS-78-census-r1 "timeout -k 30 30 python3 scripts/testing/census_identity_bands.py"`. Record its exit and every `UNREAD` line. A non-zero exit here is the expected starting point, not a red.
2. **Classify each `UNREAD` site** by item 767's three-class rule ("How to classify", above): not a solve output / separate solve / shared route, reading the call. **Re-check the four pinned lines in the two edited files** against today's text: a rule that now points at a different statement is moved to the statement it was written for, and that move is disclosed.
3. **Edit `RULES` only** (and `SEPARATE`'s tuples if a new separate-solve gate appears, with its log pattern). **No test, band, example or `src/` edit.** `mat:37` is expected to be `N` with the reason "(S0) census / precondition", as `mat:36`. Confirm it by reading `examples/materials/37_gel_box_wall_pre_equilibration_bias.py:94–110` (the imported bands) and `:268–282` (the assertion block).
4. **W2, closing.** The same command, labelled `OPS-78-census-r2`, must exit 0.

**Negative control.** The three in `control()` are unchanged and must still land as before (`mri:38` separate / same / 6 reads all 0; `mri:37` separate / different / 1.147e-13; `test_port_lumped_two_torus.py:709` shared). If `mri:38`'s rule moved, its control reading is the check that the move is right.

**Tier, ranks, wall clock.** Smoke, no solve, no container: two harness windows of ≈ 10 s each. **Class:** ops / docs. **Cost:** ≈ 3.1 measured (class median; item 767's own was 5.5) / 25 predicted min. **Timebox:** minute 25 of the item. **Depends on:** nothing. **Family:** `OPS-78` refresh, attempt 1 of 4.

**Traps already paid for.** All of item 767's. `ast` cannot follow `b.P_BAND`-style bindings (the 21 lines in `mat:18` / `mri:32`). **Never `gunzip` a record** to get green. A bare literal and `*_RTOL` / `*_BAR` constants stay out of scope. The reader is not on the headless allowlist bare: it runs only through the harness. **The harness call is the whole Bash command.**

**Scope.** A census refresh. **Not claimed:** a failure rate, a cause, a band, a design.

**Negative result.** A separate-solve gate with a recorded read above its band is a finding: name it, stop, BLOCKED, report. A split control: stop. Sites still unread at the timebox: land nothing, park the reader on `attempt/*`, and list them. No `PREDICTED RED:` line.

**Status it can move.** `OPS-78`'s row gains a record line (the table current at HEAD, with the new separate-solve gates, if any, tabled with log:line). Known-issues' `OPS-43` observation gains the refreshed counts.

**Outcome (2026-11-04T18:02Z, 12:00 CST slot): complete. The table is current at HEAD (`59784023`), 0 unread, exit 0.** Smoke, no solve, through the harness as the brief writes it, 10 s per window. **W1** `20261104T180021Z_OPS-78-census-r1.log`: `Status: 4` (`:151`), the expected start; 392 classified, **10 unread** (`:113–123`), controls as required (`:124–127`), no read above band. **W2** `20261104T180148Z_OPS-78-census-r2.log`: `Status: 0` (`:148`), 39 names, every total equal to `grep -rnw`'s (`:75`); **402 lines classified, 0 unread (`:117`): 100 separate solve, 51 shared route, 242 not a solve output, 9 out of scope** (`:76–79`); the three controls as required (`:118–121`). The reader is as committed: it was last edited before W2.

**The 10 unread sites, classified by reading the call.**
- `tests/validation/test_mat21_wall_pre_equilibration_fem.py:176, 178, 231, 278` (`S0_BAND`): R2's Σχ·Δt against t\*, a schedule census on one run → **N, "(S0) census"**, as its siblings `test_mat21_wall_duty*_fem.py`.
- `tests/package/test_wf15_study_rod_caps.py:261, 269` (`IDENTITY_BAND`): (Id) the capacitor-sheet `Study` 3×3 against `_direct_three_port(_build_rod(...))` → **S, builder different** (a `Study` build against `_build_rod`).
- `tests/package/test_wf16_sweep_fields_at.py:230, 232, 237, 239` (`IDENTITY_BAND`): (I) `|B1+|` and 10 g SAR from the `solve_sweep(fields_at=)` result against the private `_solve_port_fields` route, each its own build and factorisation (the module's docstring, `:14–19`) → **S, builder same** (both reach the mesh through the `Study` geometry's `build`).

**The four pinned rules, re-checked against today's text.** Both files shifted, so every pin had moved off its statement:
- `examples/mri/38_implant_rod_on_a_sweep.py`: `(161, 163)` → **`(162, 164)`**, the (c-f) miss (`gates["c_f"] = miss > IDENTITY_BAND` and its print). Before the move, `:162` and `:164` fell into `ALL`, "(S) sweep vs point" — the same class, the wrong gate. `control()` reads that file's class at the (S) gate line, which moved from `:143` to **`:144`**; that one-integer edit in `control()` is the same pin moved and is disclosed here. The control still lands as required (`:118`).
- `tests/package/test_wf17_implant_rod_tutorial.py` (+17 lines above the old pins): `(158,)` "page text" → **`(175, 362)`** — `:175` is the old statement, `:362` is step 28's new `bands ==` page-text assert, which the old rules had filed as S; `(264, 285)` CONTROL → **`(281, 302)`**, the (c-off) miss print and `off_misses`. Step 28's (Id) `:422, 423, 430` is new: the caps block's S against `_direct_three_port(_build_rod(...))` → **S, "(Id) caps, block vs direct"**, builder different. `:279`, `:298`, `:299` stay in `ALL`, "(S), (P=Y) rod".

**New separate-solve gates, with their reads (`READS` gained three tuples).** None reads above its band.

| gate | fixture | band | largest read | log:line | × band | headroom | builder | reads (exactly 0) |
|---|---|---|---:|---|---:|---:|---|---:|
| `WF-17 s28` (Id), caps block vs direct | rod, caps | 1e-12 | 6.377e-14 | `20261104T031801Z_WF-17-s28-W1.log:3749` | 0.0638 | 15.7 × | different | 1 (0) |
| `WF-15 s13` (Id), caps Study vs direct | rod / off, caps | 1e-12 | 4.722e-16 | `20261103T210522Z_WF-15-s13-W1.log:7345` | 0.0005 | 2117.7 × | different | 2 (1) |
| `WF-16 s18` (I), `fields_at` sweep vs private route | rod | 1e-12 | 0.000e+00 | `20261103T133533Z_WF-16-s18-W1.log:5701` | 0 | floor | same | 5 (5) |

**A changed row: `mri:38` / `WF-15 s11` (S) is no longer all-zero.** `EX-170`'s 2026-11-03 re-run reads **9.342e-14 at 64 MHz and 1.076e-13 at 128 MHz** (`20261103T225531Z_EX-170.log:5643`, `:5646`). The row now reads 1.076e-13 = 0.1076 × band, 9.3 × headroom, **8 reads, 6 exactly 0** (W2 `:130`). Item 767's table read 6 of 6 zero, and the 11-02 weekly ruling's line "`mri:38` zero 6 of 6" was true of the logs then. It is not true of the logs now. That is a record for the weekly, and it claims nothing. It bears on the recorded `mri:37` / `mri:38` difference: `mri:38`'s builder-same gate has now drifted on one window, to the same size as `mri:37`'s largest. The control still passes, because `control()` reads the original `20261102T060640Z_EX-170.log` by name, as the brief fixed it. Every other row of item 767's table reads as before (W2 `:122–145`).

**Disclosed.** (1) `mat:37` (`examples/materials/37_gel_box_wall_pre_equilibration_bias.py`) printed no `UNREAD` line. Its tight bar is `PRE_BAR` (`:100`), and `*_BAR` constants are outside the census's `*_BAND` scope (item 767's disclosure 7). Its two `*_BAND` imports, `E_BAND` (`:72`) and `INTERIOR_BAND` (`:73`), are not among the names the census resolves at ≤ 1e-11. So the `EX-180` auditor's code-only prediction did not apply, and the module needed no rule. The brief's `:268–282` assertion block was not read line by line. (2) `20261102T015645Z_WF-17-s28-W1.log` carries no (Id) print, so step 28's row rests on the one later window. (3) Rule line numbers still pin today's tree (item 767's disclosure 7). (4) The reader had no unlogged run. W1 ran on the committed reader, and W2 ran after the edits. **No band, test, example or `src/` edited.**

**Status moved (12:00 CST slot).** `OPS-78`'s row gains the record line; known-issues' `OPS-43` observation gains the refreshed row. §9 item 815 DONE.
