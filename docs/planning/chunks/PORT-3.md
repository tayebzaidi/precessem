# PORT-3 — history

Item briefs and execution history moved out of PROJECT_PLAN §9 (operator directive 2026-10-01). The §7 row is the status of record.

### §9 item 154 — 2026-10-01 (lane 1 revalidation, the 2026-09-30 weekly's lane-1 target; queued by the 10:30 review)

`PORT-3` — the `PORT-9` / `PORT-11` four-port gates as one callable (implementer; standard, `-n 1`, real; ≈ 3 measured / ≈ 47 predicted min; depends on nothing).

- **Code (new, additive): `src/fem_em_solver/ports/network_gates.py`** (it does not exist yet; the item's letter, not a rule (c) case). Band constants as defaults: `RECIPROCITY_BAND = 1e-3`, `PASSIVITY_SIGMA_TOLERANCE = 1e-9`, `C4_SPREAD_BAND = 0.005`, `POOLED_SEPARATION_FLOOR = 10.0`. `terminated_z_from_s(s, z0)`: `Z[:, j] = z0 (e_j + S[:, j]) / (1 − S_jj)` (the matched-load route's `z_matrix`; checked by hand at 4 digits against `Z_11`, `Z_21` at 10 MHz, `LEG_D0_Z_COLUMN`, and `Z_11` at 64 MHz, `…PORT-11-step2.log.gz:9259`). `circulant_classes` and `class_spread` copied verbatim from `test_port_birdcage_four_port.py:176–215`. `four_port_network_gates(s, *, z=None, z0_ohm=50.0, <bands>)` returns a frozen dataclass with reciprocity, `sigma_max`, `column_power`, spreads, pooled spread, separation, the four verdicts and `passed`, with the expressions of `test_port_birdcage_leg_offset_sweep.py:387–399`; `z=None` derives Z from S. Exported from `ports/__init__.py`.
- **Test: new `tests/ports/test_network_gates.py`.**
  - (A0) each default band equals its test constant (`RECIPROCITY_BAND`, `PASSIVITY_SIGMA_TOLERANCE`, `ADJACENT_SPREAD_BAND`, `POOLED_SEPARATION_FLOOR`).
  - (A1) for both stored S (`LEG_D_S_MATRIX_10MHZ`; `PORT11_S_64MHZ_RECORD`, which lives in `test_ans4_resolution_ladder.py:209` as a printed control and is input data here), every reading equals (`==`) the in-test expressions recomputed with the imported `_circulant_classes` / `_class_spread`.
  - (A2) records — 64 MHz: spreads 0.0573 / 0.0599 / 0.0370 %, pooled 40.1838 %, separation 671.0527×, σ_max 0.999721388, max column power 0.804704664 (`…PORT-11-step2.log.gz:9269–9271,13920–13928`); 10 MHz: 0.0553 / 0.0353 / 0.0214 %, separation 166.6766×, σ_max 0.999992805 — each at its print precision; all verdicts pass. Reciprocity on stored S reads 0.0 exactly (the printed S is symmetric): expected, not a finding (the (d3c) rule).
  - (A3) derived Z vs `LEG_D0_Z_COLUMN` and vs the logged 64 MHz Z (`:9259–9262`, restated as a version-tagged record), relative ≤ `FREQUENCY_CONTROL_BAND` (1e-6); predicted ≲ 1e-8.
- **Negative controls (planted; predicted, binary).** Each planted case passes the unplanted derived Z explicitly, so only the planted gate can move: (N1) `S[0, 1] *= 1.01` → only `reciprocal` False; (N2) `S *= (1 + 1e-6)/σ_max` → only `passive` False; (N3) `z[0, 0] *= 1.01` → self spread ≈ 1 %, only `c4` False; (N4) adjacent entries set equal to the opposite class mean → only `resolved` False.
- **Window.** W1 (≈ 15 s): `scripts/testing/run_and_log.sh PORT-3-W1 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment tests/ports/test_network_gates.py tests/ports/test_port_definition.py -s -v --tb=short'"` (`test_port_definition.py` holds D1's checklist checks, `4de76c5c`, the code this 🧪 row was opened on).
- **Traps.** Import only constants and helpers from test modules — never `test_*` names or the `larmor_rungs` fixture. Never `HEAD:` pins. Harness bare, foreground, `timeout: 660000`; no edit after the last window.
- **Scope.** Stored matrices only; no solve; the `PORT-9` / `PORT-11` modules are not rewired onto the callable.
- **Status it can move.** `PORT-3` 🧪 → ✅ on the Done-when written by the 10:30 review.
- **Negative result.** (A1) not bitwise: the expression order differs — fix once, then park. (A3) > 1e-6: the derived-Z relation is wrong — print every entry and stop, BLOCKED. A planted case flipping two gates: report and stop.

**Outcome (2026-10-02 07:30 slot, third item, `implementer`): BLOCKED, attempt 1 of 4.** Branch `attempt/PORT-3-20261002T125518Z`. W1 is `20261002T125356Z_PORT-3-W1.log`, `-n 1`: 25 passed, 1 failed, 4 skipped, 49 s.
- (A0) green.
- (A1) every reading `==` at 10 and 64 MHz. No expression-order fix was needed.
- (A2) spreads, pooled spread, σ_max and column power match at print precision. Reciprocity reads 0.0 on the stored S, as expected.
- (A3) derived Z: worst entry 2.5e-10 against `LEG_D0_Z_COLUMN` and 3.4e-10 against the 64 MHz Z, restated as `PORT11_Z_64MHZ_RECORD`.
- (N2) flips only `passive`, (N3) only `c4` (self spread 1.006 %), (N4) only `resolved` (separation 0.646×).
- **(N1) red, the brief's own stop condition.** `S[0,1] *= 1.01` gives reciprocity 3.3e-3, which flips `reciprocal`. It also lifts σ_max to 1.000695, so `passive` flips too: the 64 MHz record has only 2.8e-4 of passivity headroom, and a 1 % change to an entry of magnitude ≈ 0.42 adds ≈ 1e-3. No real factor isolates the gate: reciprocity crosses 1e-3 at ≈ 1.003, and σ_max crosses 1 at ≈ 1.0029. Re-brief with a plant that keeps σ_max, such as a phase rotation of `S[0,1]` alone or an antisymmetric perturbation, pre-computed on the record.
- **For the review (a band question):** the executor compares the (A2) separation to its record at 1e-6 relative instead of print precision. The readings are 671.05257 against 671.0527 and 166.67665 against 166.6766. Both records were printed from unrounded Z, while our input is S rounded to 10 digits, and separation is a ratio of spreads near 5e-4. Rule whether that comparison stands.
- The planted cases use the 64 MHz record, because at 10 MHz (N3) would also flip `resolved` (pooled 9.2 %).


### §9 item 178 — 2026-10-02 (lane 1; PORT-3 attempt 2 of 4)

*Drafted by a read-only `Plan` agent and reviewed by the 2026-10-02 10:30 daily review.*

**Outcome (2026-10-02 15:00 slot, `implementer`): DONE, attempt 2 of 4; `PORT-3` 🧪 → ✅.** W2 is `20261002T201726Z_PORT-3-W2.log`, real, `-n 1`: 26 passed, 0 failed, 4 skipped, Status 0, 51 s — exactly as predicted.
- (A0) green; (A1) every reading `==` at 10 and 64 MHz (`:69–93`).
- (A2) every named reading within ½ ulp (`:95–119`). Tightest: 64 MHz self spread 0.0573415 % vs 0.0573, gap 4.147e-5 ≤ 5e-5 (margin 8.5e-6). σ_max gaps 2.60e-10 / 3.73e-10 ≤ 5e-10 at 64 / 10 MHz; column power 1.78e-10. Separation printed, not asserted: 671.0525735× vs 671.0527 (2.53 half-ulps, 1.885e-7) and 166.6766550× vs 166.6766 (1.10 half-ulps, 3.298e-7), as predicted. All verdicts pass.
- (A3) derived Z within 1e-6 at both frequencies (`:121–143`).
- (N1) φ = 1°: reciprocity 1.878257e-2 (= the predicted 2|sin φ|K, K = 0.538108; 18.78× band); σ_max(S′) − σ_max(S) = 0.0; max |Δ column power| 1.1e-16; spreads / pooled / separation unchanged; only `reciprocal` False (`:145–150`).
- (N2) only `passive` (σ_max 1.000001); (N3) only `c4` (self 1.0058 %, separation 39.95×); (N4) only `resolved` (separation 0.6460×) (`:152–162`).
- Code by path from `attempt/PORT-3-20261002T125518Z`; `network_gates.py` / `__init__.py` unchanged; only the test edited. The branch was deleted locally after landing.

**Read by the review (item 154).** The diagnosis holds on the parked module. `PORT11_S_64MHZ_RECORD` has only 2.786e-4 of passivity headroom (σ_max 0.9997213883), so a real scale of one off-diagonal entry cannot isolate reciprocity: reciprocity crosses 1e-3 at a factor of 1.003, and σ_max crosses 1 at ≈ 1.0029. That is a brief defect, not a callable defect. The callable read every gate correctly: (N2)–(N4) were green, and (N1) printed the 3.3e-3 / 1.000695 that its plant implies. (N1) is now re-planted as a unitary similarity, which preserves σ_max and the column powers exactly and breaks only symmetry. **(A2) is ruled:** see "Ruling on (A2)" below. No gate band changes.

`PORT-3`: the `PORT-9` / `PORT-11` four-port gates as one callable. This is item 154's parked code with (N1) re-planted and (A2) re-ruled.
- **Executor:** implementer. **Tier:** standard, **`-n 1`**, real, as the Done-when says.
- **Cost:** ≈ 4 measured / ≈ 20 predicted min.
- **Depends on:** nothing. **Attempt 2 of 4.**

**Code (by path, no `src/` edit):**
- Run `git checkout attempt/PORT-3-20261002T125518Z -- src/fem_em_solver/ports/network_gates.py src/fem_em_solver/ports/__init__.py tests/ports/test_network_gates.py` (branch commit `057d0b67`). `git diff --stat` should list exactly those three paths.
- The branch carries no log or record row: its W1 log and `test-results.md` row are already on `main`. Take nothing else, and never path-checkout `docs/testing/test-results.md`.
- `main` has not touched `src/fem_em_solver/ports/` or `tests/ports/` since the branch base `54f9ce89`.
- `network_gates.py` and `__init__.py` land unchanged. Only `tests/ports/test_network_gates.py` is edited, in two places:
  1. **(N1) re-planted.**
     - Replace `test_n1_planted_asymmetry_fails_only_reciprocity` with `test_n1_planted_nonreciprocal_phase_fails_only_reciprocity`.
     - Set `N1_PHASE_RAD = np.deg2rad(1.0)`. On `_base()`, for k in 1..3: `s[0, k] *= np.exp(1j*N1_PHASE_RAD)` and `s[k, 0] *= np.exp(-1j*N1_PHASE_RAD)`. Write it as an elementwise loop and leave the diagonal untouched bitwise. Do not build it as a `D @ S @ D.conj().T` matmul.
     - The plant is S′ = D S Dᴴ with D = diag(e^{jφ}, 1, 1, 1): a lossless non-reciprocal phase shifter in port 1's feed.
     - Call `_planted("N1", s, z, "reciprocal")` with the **unplanted** derived `z`.
     - Print σ_max(S′) − σ_max(S) and the max |Δ column power| against the unplanted S.
     - The docstring cites this item and the identity.
  2. **(A2), per the ruling below.**
     - Store each record as the string the log printed (e.g. `"0.0573"`, `"0.999721388"`).
     - Set the tolerance to `_half_ulp(text) = 0.5·10^(−decimals after the point)`.
     - Assert `|reading − float(text)| ≤ _half_ulp(text)` for: the three class spreads at both frequencies, the pooled spread at 64 MHz, σ_max at both frequencies, and max column power at 64 MHz.
     - **Delete** the `sep_rel <= FREQUENCY_CONTROL_BAND` assertion and its comment. Print the separation beside its record, as the deviation in half-ulps and as a relative figure, labelled "printed, not asserted (item 178 ruling)".
     - `FREQUENCY_CONTROL_BAND` stays imported for (A3).

**Anchors (asserted):**
- **(A0)** Default bands `==` the test constants.
- **(A1)** Every reading `==` the in-test arithmetic at 10 and 64 MHz. Measured bitwise in W1, `20261002T125356Z_PORT-3-W1.log:69–92`.
- **(A2)** Each reading within ½ ulp of its printed record (W1 readings `:98–118`; every one is inside):
  - **64 MHz:**
    - spread self 0.057341 % vs 0.0573 %: gap 4.15e-5 ≤ 5e-5, the tightest margin;
    - spread adjacent 0.059882 % vs 0.0599 %;
    - spread opposite 0.037040 % vs 0.0370 %;
    - pooled 40.183832 % vs 40.1838 %;
    - σ_max 0.999721388260 vs 0.999721388;
    - max column power 0.804704664178 vs 0.804704664.
  - **10 MHz:**
    - spreads 0.055266 / 0.035303 / 0.021392 % vs 0.0553 / 0.0353 / 0.0214 %;
    - σ_max 0.999992804627 vs 0.999992805: gap 3.7e-10 ≤ 5e-10.
  - Every verdict passes. Reciprocity reads 0.0 exactly, because the stored S is symmetric.
  - Separation, printed only: predicted 2.53 half-ulps (1.885e-7 relative) at 64 MHz and 1.10 half-ulps (3.298e-7) at 10 MHz. Its `resolved` verdict (≥ 10) is asserted through all-pass, at margins of 67× and 16.7×.
- **(A3)** Derived Z vs `LEG_D0_Z_COLUMN` and `PORT11_Z_64MHZ_RECORD` ≤ 1e-6 relative. Predicted worst 2.5e-10 / 3.4e-10 (`:144`).
- **The identity behind (N1).** For a unitary diagonal D, σ(D S Dᴴ) = σ(S) and |S′_ik| = |S_ik|, so `passive` cannot move. The plant's reciprocity residual is ‖S′ − S′ᵀ‖_F / ‖S′‖_F = 2|sin φ|·K, with K = √(2 Σ_{k=1..3} |S_0k|²) / ‖S‖_F = 0.538108 on the 64 MHz record.

**Negative controls (the 64 MHz record, each passed the unplanted derived `z`).** Each verdict pattern is **asserted**: it is the Done-when's "each fail exactly their own gate", and N2–N4 are backed by W1 `:151–160`. Every numeric factor is **predicted** and printed beside its reading (rule (e)).
- **(N1)** φ = 1°.
  - Reciprocity predicted 1.8783e-2, i.e. 18.8× `RECIPROCITY_BAND`.
  - The crossing is at φ ≈ 0.0532°; every φ in (0.0532°, 179.95°) isolates the gate.
  - |Δσ_max| ≲ 1e-15 (round-off), and the max |Δ column power| ≤ 2e-16.
  - Spreads, pooled spread and separation are identical to the unplanted ones.
  - Only `reciprocal` reads False.
- **(N2)** `S *= (1 + 1e-6)/σ_max`. This puts the plant 1000× over the `passive` tolerance.
  - Reciprocity stays 0.0 exactly; max column power is 0.805155 < 1.
  - No other gate has a crossing on this plant (`:151`).
- **(N3)** `z[0,0] *= 1.01`, with S untouched.
  - `c4` crosses at a factor of ≈ 1.00492 and `resolved` at ≈ 1.0405.
  - Predicted self spread 1.0058 % (2.0× the 0.5 % band) and separation 39.95× (`:155`).
- **(N4)** The adjacent entries set to the opposite-class mean, with S untouched.
  - Predicted separation 0.6460×, 15.5× under the floor (`:159`).
  - The spreads only shrink, so `c4` cannot flip.
- The planted cases stay at 64 MHz. At 10 MHz the pooled spread is 9.2 %, so (N3) would also flip `resolved`.

**Window W2 (≈ 50 s; W1 took 49 s).** One bare, foreground harness call, with Bash `timeout: 660000`:

`scripts/testing/run_and_log.sh PORT-3-W2 "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src:/workspace timeout -k 30 120 mpiexec -n 1 python3 -m pytest tests/environment tests/ports/test_network_gates.py tests/ports/test_port_definition.py -s -v --tb=short'"`

Predicted: 30 collected, **26 passed, 0 failed, 4 skipped** (the complex-only environment checks), Status 0.

There is no PREDICTED RED line: any red stops the slot.

**Traps:**
- Pass the **unplanted** derived `z` to every plant. Re-deriving Z from the (N1) S′ rotates the off-diagonal Z by e^{±jφ}. That opens an adjacent-class gap of ≈ 3.5 % and flips `c4` too.
- Build (N1) elementwise, never by matmul.
- Do not assert the separation at any band. In particular, do not reuse `FREQUENCY_CONTROL_BAND`.
- Run no probe outside the harness. Item 154 ran a bare `docker compose exec … python3 -c`, and host `python3` is denied.
- Import only constants and helpers from the test modules, never `test_*` names or the `larmor_rungs` fixture. No `HEAD:` pins.
- The harness call is the **whole** Bash command: bare, repo-relative, nothing chained or piped.
- Run with `-s`. Make no edit after W2 (rule (i)).
- `git diff --stat` shows only insertions on record files.

**Scope.** Stored matrices only; no solve. The `PORT-9` / `PORT-11` modules are not rewired onto the callable. No band moves: 1e-3, 1 + 1e-9, 0.5 %, 10× and (A3)'s 1e-6 all stay as they are.

**Status it can move.** `PORT-3` 🧪 → ✅ on the §7 Done-when, as amended by this review with the (A2) parenthetical. Every clause has its check:
- the callable with test-constant defaults: (A0);
- bitwise readings: (A1);
- print precision and all-pass: (A2);
- Z to 1e-6: (A3);
- each plant failing alone: (N1)–(N4);
- `test_port_definition.py` in the same window, real, `-n 1`: W2.

**On landing green:**
- Commit the three code paths with the W2 log, and append the `test-results.md` row by hand if one is missing.
- Write the attempts.md entry and mark item 178 DONE in the same commit.
- Then delete the branch with `git branch -D attempt/PORT-3-20261002T125518Z` (local only) and journal the deletion.

**Negative result.**
- **(N1) flips `passive`, or |Δσ_max| > 1e-12:** the callable's σ does not come from the S it was passed. That is a code defect: print every σ and stop.
- **(N1) reciprocity ≤ 1e-3:** the gate is blind to phase non-reciprocity. Stop.
- **(N1) flips `c4` / `resolved`:** the plant re-derived Z. Fix it once in the test, re-run once, then park.
- **(A1) is not bitwise, or (A2) misses by more than ½ ulp:** the module drifted. Stop.
- **On any stop:** park on a new `attempt/PORT-3-<UTC>` and mark item 178 BLOCKED with its unblock condition in the same commit (rule (d)). That spends attempt 2 of 4.

**Ruling on (A2) (2026-10-02 10:30 review).** The Done-when's words are "reproduce the logged **spreads / σ_max** at print precision". So (A2) asserts |reading − record| ≤ ½ unit in the last printed digit, on exactly those quantities, plus max column power, which item 154 listed. All of them pass at that bound in W1 (`:98–118`).

The separation is not named in the Done-when. Its printed digit lies below the input's information floor. S is rounded to 10 digits, which gives ≈ 3e-10 relative in Z (`:144`); dividing by a worst spread of ≈ 6e-4 leaves the ratio uncertain to ≈ 1e-7 to 1e-6 relative, while the record prints it to 7.5e-8 / 3.0e-7.

The executor's `FREQUENCY_CONTROL_BAND` assertion on the separation is therefore **struck, not ratified**: it was a band chosen in-slot and borrowed from an unrelated gate. The separation is printed against its record, and its `resolved` gate stays asserted. This is the same class of ruling as `IMP-1` item 165's: a reproduction band against a printed value is never tighter than ½ ulp of the printed mantissa. No gate band changes.

**Audit, 2026-10-02 18:00 daily review — PASS** (`auditor`; the review re-read W2:111, gap 4.147e-05 ≤ 5e-05). Evidence: `20261002T201726Z_PORT-3-W2.log:69–93,95–119,121–143,145–162,163–170,234,236–238`; record strings `20260826T110434Z_PORT-11-step2.log.gz` (zcat) :4675, :9269–9271, :13920–13928; asserts `tests/ports/test_network_gates.py:180–184,211–223,244–245,253–264`. `git diff 057d0b67 ab060446` on the two `src/` files is empty. Caveats: the log header is the parent `0862918b` (dirty tree; timestamps consistent); the 64 MHz self-spread margin is thin (8.5e-6 of 5e-5); the (N1) Δσ_max 0.0 is printed, not asserted.
