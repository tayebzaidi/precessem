# MAT-3 — history file

*Created by the 2026-10-25 18:00 daily review for a §9 brief (`OPS-46`). The §7 `MAT-3` row (🟡, smoke) carries the chunk's status and its Done-when. Step 1 landed on 2026-09-29 as §9 item 80 (`a097bed2`; record in `docs/planning/plan-archive.md`). That step covered a single Debye pole on `TH-6`'s lossy plane wave, with α/β within 0.026 % / 0.055 % at 127.74 MHz and both controls missing α by 12.0 % / 22.0 %. Its example is `EX-145` (`th:20`). Step 2 waits on the operator's Gabriel table under `docs/references/` (dashboard item 13).*

### §9 item 646 — 2026-10-25 (§9 rule (5) lane 3, `MAT-3` step 1b: a two-pole Cole–Cole medium, α > 0, on `TH-6`'s plane wave)

*Written by the 2026-10-25 18:00 daily review (lane-3 fact sheet). Lane 3 lists `MAT-3` by name. The anchor is the closed-form plane wave k = k₀√ε̂(ω) (the same as `TH-6`'s and step 1's), with ε̂ the Cole–Cole sum (Cole & Cole 1941, as cited in `src/fem_em_solver/materials/dispersion.py`). The band is `TH-6`'s 5 %, imported. Nothing here sets or moves a band, tolerance, record, Done-when or target. **Why:** step 1's FEM gate runs one Debye pole (α = 0), so **no FEM run has exercised a broadened pole (α > 0) or more than one pole**. Those are the cases Gabriel tissues use, and the evaluator's (jωτ)^(1−α) polar branch (`dispersion.py:110–113`) has only been tested by numpy identities. The parameters are synthetic, as step 1's were, so the step needs no table.*

**The medium (fixed by this brief; the executor re-derives the predictions in W0).** `ColeColeMaterial(epsilon_inf=4.0, poles=((50.0, 1/(2π·100 MHz), 0.30), (4000.0, 1/(2π·1 MHz), 0.35)), sigma_ionic=0.2)`. This is a muscle-like shape: a γ-like pole near 100 MHz and a large β-like pole near 1 MHz, both broadened. Predictions from the formula, computed by this review with the standard library:

| f | ε̂ | α_k (Np/m) | β_k (rad/m) | Debye-reduced (α_n → 0) α_k miss | wrong-f α_k miss |
|---|---|---|---|---|---|
| 63.87 MHz | 181.505 − 283.983j | 11.8044 | 21.5541 | 17.15 % | — |
| 127.74 MHz | 118.420 − 182.618j | 18.8581 | 34.7047 | 19.91 % | 29.77 % (ε̂(63.87) frozen, solved at 127.74) |

λ at 127.74 MHz is 2π/34.70 ≈ 0.181 m. At 24³ on the 0.1 m box that is ≈ 43 cells per wavelength; step 1 had ≈ 61 at β = 24.4.

**Anchors (asserted, `TH-6`'s imported 5 % band `BAND` from step 1's module):**
- **(Ev)** `complex_permittivity(f)` equals an **independent standard-library evaluation** of ε∞ + Σ Δεₙ/(1 + (jωτₙ)^(1−αₙ)) − jσ/(ωε₀), written in the test with `cmath` and `(1j*w*tau)**(1-alpha)` (principal branch, which is the evaluator's polar form for ωτ > 0), to ≤ 1e-12 relative at both frequencies. Then step 1's `_truth` cross-checks: `_wavenumber` / `_analytic_alpha_beta` from `TH-6`'s module at the frozen (σ_eff, ε′) equal k₀√ε̂ to 1e-12 / 1e-10.
- **(i)** 127.74 MHz, 12³ → 24³: rel L2 decreasing, rate > 0.9, 24³ rel L2 < `BAND`, α_k and β_k within `BAND` at 24³.
- **(ii)** 63.87 MHz, 24³: α_k and β_k within `BAND`.

**Negative controls (asserted, 16³, each a self-consistent solve; ceilings computed above):**
- **(a′) Debye-reduced:** the same Δε, τ and σ with every αₙ = 0 (`ColeColeMaterial.debye(4.0, poles=((50.0, τ₁), (4000.0, τ₂)), sigma_ionic=0.2)`), solved at 127.74 MHz, misses the Cole–Cole α_k by > `BAND`. The closed-form ceiling is 19.91 %, ≈ 4× the band. This is the control that a broadening-blind evaluator (α ignored) would fail.
- **(b′) wrong frequency:** ε̂(63.87 MHz) frozen and solved at 127.74 MHz misses α_k by > `BAND` (ceiling 29.77 %).
- The β_k misses of both controls are printed only, as step 1 does.

**Route — a new module, nothing edited.** `tests/validation/test_mat3_cole_cole_plane_wave.py`, modelled on `tests/validation/test_mat3_debye_plane_wave.py`. Import by name from it: `F_HIGH`, `F_LOW`, `BAND`. From `tests/validation/test_lossy_plane_wave.py` import `MU_R`, `_analytic_alpha_beta`, `_solve_plane_wave`, `_wavenumber`. From `fem_em_solver.materials` import `ColeColeMaterial`. Write the module's own `_truth(f, material)` and `_solve(n, f, material, material_f=None)` taking the material as an argument; do not reuse step 1's, which close over `DEBYE`. Never edit step 1's module, `TH-6`'s, or `dispersion.py`.

**Windows.**
- **W0** (smoke, `-n 1`, complex): the (Ev) node alone, plus printing the table above (ε̂, α_k, β_k and both ceilings). No mesh. Assert the two ceilings > 2 × `BAND` before any FEM window, so a control that cannot miss is caught for free.
- **W1** (standard, `-n 2`, complex, `timeout -k 30 180`), after `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0: `scripts/testing/run_and_log.sh MAT-3-step1b-W1 "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment tests/validation/test_mat3_cole_cole_plane_wave.py tests/validation/test_mat3_debye_plane_wave.py -v -s --tb=short'"`. Step 1's module re-runs unedited in the same window as the regression.
- **Prediction:** step 1's window ≈ 33 s for five solves. This one has five solves plus step 1's five ⇒ ≈ 60–90 s. Over 180 s is a kill: park, never raise the timeout.

**The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped.**

**Executor:** implementer. **Tier:** standard, `-n 2`, complex. **Class:** lane-3 closed form. **Cost:** ≈ 7.5 measured / 30 predicted slot-minutes. **Depends on:** nothing. **Family:** step 1b, attempt 1 of 4. **Daytime:** no.

**Traps already paid for.**
- `(jωτ)^(1−α)` must be the principal branch with arg = π(1−α)/2. A `(w*tau)**(1-alpha) * 1j**(1-alpha)` spelling is the same. A `(1j*w*tau)**(-(alpha-1))` sign slip conjugates the pole, and (Ev) is there to catch it.
- e^{+jωt}: the ionic term is −jσ/(ωε₀), and Im k < 0 on the decaying wave (step 1's `_truth` asserts `k.imag < 0`).
- `ufl.max_value` and ordering comparisons on complex operands do not compile (`OPS-22`). The material enters only through `HomogeneousMaterial`, so no new UFL is needed.
- `assemble_scalar` is rank-local, so reduce before asserting. pytest captures prints without `-s`.
- Never pipe the harness command, and never delete a log.

**Scope.** A two-pole broadened Cole–Cole medium through the frozen-frequency route, against the closed form. **Not claimed:** any tissue's published σ / ε_r (step 2, the operator's table), a time-domain or broadband solve, temperature dependence. `MAT-3` stays 🟡: its ✅ needs step 2.

**Negative result.** If (Ev) is red, it is an evaluator defect: report it, park on `attempt/MAT-3-step1b-*`, item 646 BLOCKED, and open a known-issues entry. Never edit `dispersion.py` in this item. If (i) or (ii) is red with (Ev) green, the FEM misses the closed form on a broadened medium: record and park. If a control is inside the band, the fixture cannot see broadening: record it, stop. No `PREDICTED RED:` line.

**Status it can move.** The §7 `MAT-3` state line gains "step 1b: a two-pole Cole–Cole medium (α = 0.30 / 0.35) on `TH-6`'s plane wave, green" as a registered record (glyph stays 🟡), or a known-issues entry (red). An example (a `th:` or `mat:` broadened-tissue plane wave) is the next review's to open under the ramp.

**Outcome (2026-10-26 02:17 UTC, 21:00 CDT slot, attempt 1 of 4): BLOCKED, parked `attempt/MAT-3-step1b-20261026T0217Z`.** The module is written to the brief and has been unedited since W0.
- **W0** (`-n 1`, 5 s) is green. (Ev) evaluator vs `cmath`: 1.687e-16 / 1.306e-16. ε̂, α_k and β_k match the table to every printed digit. Ceilings: 19.9130 % (a′, 127.74 MHz), 17.1466 % (a′, 63.87 MHz) and 29.7715 % (b′).
- **W1** (as written) was killed at 181 s, exit 124.
- **Bisect:** the full `tests/environment` costs ≈ 145 s. `test_dolfinx_api_migration.py` alone takes 144.42 s at `-n 1` (known-issues 2026-10-26).
- **Probe** (`test_complex_mode.py` + this module, `-n 2`, 45 s, all green), not the closing window:
  - (i) rel L2 9.496e-02 → 4.750e-02, rate 0.9994; α 18.86589 (0.0415 %), β 34.73752 (0.0947 %).
  - (ii) α 11.80833 (0.0334 %), β 21.56536 (0.0520 %).
  - (a′) Debye-reduced misses α by 20.0464 % (β 42.84 %, printed only).
  - (b′) wrong-f misses α by 30.1096 % (β 10.67 %, printed only).
- **Step 1's module** reproduces its 09-29 record to the digit in 35 s.
- **Logs:** `20261026T020127Z_MAT-3-step1b-W0.log`, `…020141Z_…-W1.log`, `…020520Z_…-W1probe.log`, `…020620Z_…-W1probe2.log`, `…020945Z_…-W1probe3.log`, `…021046Z_…-bisectA.log`, `…021129Z_…-bisectB.log`, `…021215Z_…-bisectC.log`, `…021500Z_…-envprice.log` (all on the branch).

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 646 — `MAT-3` step 1b — per `docs/planning/chunks/MAT-3.md` "§9 item 646". Read the brief and any review ruling beneath it first; the rulings hold. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Write a new module; never edit `test_mat3_debye_plane_wave.py`, `test_lossy_plane_wave.py` or `dispersion.py`. Spawned executors run every window foreground with Bash `timeout: 660000`."

### Review ruling 2026-10-26 03:00 (item 646)

Binds over the brief above.
- **Not landed by this review.** An `auditor` read the parked windows (INCONCLUSIVE on one item). Every anchor is an executed assert in the parked module (`test_mat3_cole_cole_plane_wave.py:109`, `:113`, `:117`, `:132–133`, `:178–186` at `1ffc2790`). Each probe window is footered, Status 0 (`W1probe.log:99`, `:121`; `W1probe3.log:92`, `:111`; `envprice.log:119`). Step 1's record reproduces to the digit (`W1probe3.log:67–73` against `20260929T183252Z_MAT-3-step1.log.gz:159–165`). **But nothing pins the module's bytes at probe time.** The module was untracked when the probe ran (02:05 UTC) and was committed at 02:17. So the probe windows cannot close the step. Item 646 is re-posed as **§9 item 648**, attempt 2 of 4, and the parked module runs unedited.
- **The environment-first ruling** (journal `2026-10-26-daily.md` §4; known-issues 2026-10-26) applies. A solve window puts `tests/environment/test_complex_mode.py` first. The full `tests/environment` runs once per item, in its own window `W-env` (`-n 1`, `timeout -k 30 300`), and never shares a window with a solve.
- **Correction to the composition:** `envprice` ran `test_dolfinx_api_migration.py` only (`envprice.log:12`), not the full `tests/environment`. The other two environment modules were green at `-n 2` in bisectA / bisectB (`…021046Z_…-bisectA.log:101`, `…021129Z_…-bisectB.log:116`).

### §9 item 648 — 2026-10-26 (§9 rule (5) lane 3, `MAT-3` step 1b re-posed: the parked Cole–Cole module run unedited under the environment-first ruling)

*Written by the 2026-10-26 03:00 daily review. This re-poses item 646 (attempt 2 of 4). The anchor, the medium, the controls, the scope, the negative-result disposition and the status it can move are item 646's above, unchanged. Nothing here sets or moves a band, tolerance, record, Done-when or target.*

**Land the parked work first (code and logs by path, record rows by hand).**
1. `git checkout attempt/MAT-3-step1b-20261026T0217Z -- tests/validation/test_mat3_cole_cole_plane_wave.py`, plus the nine `docs/testing/logs/20261026T02*_MAT-3-step1b-*.log` files the branch adds (new files, so a path checkout is safe for them).
2. Append the branch's nine `test-results.md` rows **by hand**. `git diff main...attempt/MAT-3-step1b-20261026T0217Z -- docs/testing/test-results.md` shows them. **Never path-checkout `test-results.md`** (the `8d4cf58` trap). `git diff --stat` on it must show insertions only.
3. Do not edit the module. Its `sha256sum` is printed in W0 below, so the closing windows pin its bytes.

**Windows** (foreground, Bash `timeout: 660000`, `-s`; the harness call is the whole Bash command: bare, repo-relative, nothing chained or piped).
- **W0** (smoke, `-n 1`): `scripts/testing/run_and_log.sh MAT-3-step1b-W0b "docker compose exec -T fem-em-solver bash -lc 'cd /workspace && sha256sum tests/validation/test_mat3_cole_cole_plane_wave.py'"`. Then `git show attempt/MAT-3-step1b-20261026T0217Z:tests/validation/test_mat3_cole_cole_plane_wave.py` through `sha256sum` on the host must give the same digest; state both in the outcome.
- **W1** (standard, `-n 2`, complex, `timeout -k 30 180`), after `docker compose -f docker/docker-compose.yml exec -T fem-em-solver pgrep -c python3` reads 0: `scripts/testing/run_and_log.sh MAT-3-step1b-W1b "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 180 mpiexec -n 2 python3 -m pytest tests/environment/test_complex_mode.py tests/validation/test_mat3_cole_cole_plane_wave.py tests/validation/test_mat3_debye_plane_wave.py -v -s --tb=short'"`. **Prediction ≈ 70–90 s** (probe 45 s + probe3 35 s, measured). The pass count must be 4 + 2 + 3 nodes per rank (complex_mode's four, the two new nodes, step 1's three), and every node must be PASSED on both ranks. Count the tokens per rank, not the lines (the interleave trap below).
- **W-env** (heavy by measurement, `-n 1`, `timeout -k 30 300`): `scripts/testing/run_and_log.sh MAT-3-step1b-Wenv "docker compose exec -T fem-em-solver bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && FEM_EM_REQUIRE_COMPLEX=1 PYTHONPATH=/workspace/src timeout -k 30 300 python3 -m pytest tests/environment -v --tb=short --durations=5'"`. **Prediction ≈ 150–165 s** (`envprice` 146 s for the migration module alone, plus < 10 s for the rest). This window sweeps the new module, which W1 does not.

**Executor:** implementer. **Tier:** standard (W1) / heavy (W-env), `-n 2` / `-n 1`, complex. **Class:** lane-3 closed form. **Cost:** ≈ 7.5 measured / 20 predicted slot-minutes. **Depends on:** nothing. **Family:** step 1b, attempt 2 of 4. **Daytime:** no.

**Traps already paid for** (item 646's, plus):
- **The rank-interleave trap** (`log-pathologist`, 2026-10-26 03:00, on item 644's W1): in a 2-rank stream a `PASSED` token printed at the end of a line can belong to the other rank's node. Count `PASSED` tokens per test ID (one per rank) and read the two pytest summary lines. Never infer one rank's verdict from line adjacency.
- A W1 kill at 180 s with W1's composition is a real overrun (≈ 2× the measured sum). Park, never raise the timeout, and report the last node each rank printed.

**Status it can move:** as item 646. The §7 `MAT-3` state line gains the step-1b record (glyph stays 🟡). §10 lane 3's step-1b note → DONE. Delete `attempt/MAT-3-step1b-20261026T0217Z` in the closing commit, once its module, logs and rows are on `main`. The example (`EX-157`, `th:23`) is §9 item 652 and waits on this item.

**Executor spawn prompt (verbatim, if delegated).** "Execute §9 item 648 — `MAT-3` step 1b re-posed — per `docs/planning/chunks/MAT-3.md` "§9 item 648". Read the brief, item 646's brief above it and the 2026-10-26 03:00 review ruling first; the ruling holds. Land the parked module and logs by path and the test-results rows by hand; never edit the module. The harness call is the whole Bash command: bare, repo-relative, nothing chained or piped. Spawned executors run every window foreground with Bash `timeout: 660000`."

**Outcome (2026-10-26, 09:41 UTC / 04:30 CDT slot, attempt 2 of 4) — DONE, green; `MAT-3` stays 🟡 (step 2 waits on the Gabriel table).** Run inline by the slot. Landed by path from `attempt/MAT-3-step1b-20261026T0217Z`: the module and its nine logs; the branch's nine `test-results.md` rows appended by hand (`git diff --stat`: 9 insertions, 0 deletions). The module is unedited. **W0b** `20261026T094144Z_MAT-3-step1b-W0b.log`: container `sha256sum` 0fa956168ac16dedf4b0f7c38564fd76d20e50b4c98c82ab606dd38a38b816b1, and the host's `git show attempt/…:tests/validation/test_mat3_cole_cole_plane_wave.py | sha256sum` gives the same digest. **W1b** `20261026T094151Z_MAT-3-step1b-W1b.log` (`-n 2`, complex, 60 s, Status 0). Both rank summaries read "7 passed", and there are 14 `PASSED` tokens, two per test ID. **Composition note:** the brief predicted 4 + 2 + 3 nodes per rank, but step 1's module has **one** test function carrying three readings, so the count is 4 + 2 + 1 = 7. The brief's error, not a missing node. Readings, to the digit of the 02:05 probe: (Ev) 1.687e-16 / 1.306e-16; (i) rel L2 9.49618e-02 → 4.75018e-02, rate 0.9994, α 18.86589 (0.0415 %), β 34.73752 (0.0947 %); (ii) α 11.80833 (0.0334 %), β 21.56536 (0.0520 %); (a′) Debye-reduced α miss 20.0464 % (ceiling 19.91 %; β 42.84 % printed); (b′) wrong-f α miss 30.1096 % (ceiling 29.77 %; β 10.67 % printed). Step 1's regression reproduces its record (α 0.0262 % / β 0.0554 %, controls 11.9990 % / 21.9949 %). **W-env** `20261026T094304Z_MAT-3-step1b-Wenv.log` (`-n 1`, complex, 103 s against a predicted 150–165, Status 0): 12 passed. The attempt branch is deleted in this slot, now that its content is on `main`.

### Review audit 2026-10-26 18:00 (item 648, step 1b) — PASS

`auditor` PASS on `075879df`: the module's sha256 `0fa95616…b816b1` printed in `20261026T094144Z_MAT-3-step1b-W0b.log` equals the committed blob and the parked blob `1ffc2790`; (Ev) 1.687e-16 / 1.306e-16 (`…W1b.log:65–66`), (i) 0.0415 % / 0.0947 % rate 0.9994 (`:74`), (ii) 0.0334 % / 0.0520 % (`:75`), (a′) 20.0464 % and (b′) 30.1096 % missed (`:76–77`), both ranks "7 passed" (`:93`, `:96`), Status 0 / 60 s (`:98–101`); `test-results.md` 12 insertions, 0 deletions. Caveats recorded: the row's tier read **smoke** against measured 60–103 s windows (re-declared **standard** in the row by this review); `BAND` is a local `0.05` restating `TH-6`'s 5 %, not a symbol import (value equal); W-env sweeps the module statically, so W1b is its only execution evidence. Consumed next by `Study` (`WF-15` step 8, §9 item 653).

## §7 row moved verbatim 2026-10-26 (rotate_plan_archive.py chunks)

 Debye/Cole-Cole dispersion models — **feature ladder C3** (operator directive 2026-09-04; matters once tissue-property maps replace the gel, Phase 7). **Done-when (written 2026-09-29 10:30 review, rule (j); step 1 queued, §9 lane 3 — the lane-3 opening ratified by the 2026-09-30 weekly).** *Step 1:* `materials/dispersion.py` evaluates Cole–Cole / Debye ε̂(ω) = ε∞ + Σ Δεₙ/(1 + (jωτₙ)^(1−αₙ)) + σ/(jωε₀), e^{+jωt} (Gabriel, Lau & Gabriel 1996, *Phys. Med. Biol.* 41:2271). Cole–Cole at α = 0 equals Debye to ≤ 1e-14. The Debye semicircle and Cole–Cole depressed-arc identities hold to ≤ 1e-12·Δε. At Δε = 0 it equals `complex_permittivity`, and −Im ε̂ > 0 everywhere. `TH-6`'s lossy plane wave solved with a Debye material at 63.87 and 127.74 MHz matches k₀√ε̂(ω) in α and β within `TH-6`'s imported 5 % band. Wrong-frequency and dispersion-dropped controls each miss α by more than the band. *Step 2:* one Gabriel 1996 tissue reproduces the paper's σ and εr at 64 and 128 MHz to a stated band, gated only after the parameter table and published values are committed to `docs/references/`. ✅ = steps 1 and 2 green. **Out of scope:** time-domain (ADE / convolution) dispersion, spatially varying tissue maps, T-dependent σ (`MAT-5`). **Step 1 green 2026-09-29 (item 80):** N1–N6 hold to ≤ 2.2e-16 relative (1.2e-14 absolute on the arcs), and the Debye plane wave's α / β are within 0.026 % / 0.055 % at 127.74 MHz (24³, rate 0.9998) and 0.042 % / 0.048 % at 63.87 MHz, with the two controls missing α by 12.0 % / 22.0 % (`20260929T183243Z_MAT-3-step1-numpy.log`, `20260929T183252Z_MAT-3-step1.log`). Step 2 waits on the Gabriel table landing in `docs/references/`. Step 1b (a two-pole Cole–Cole medium, α > 0): item 646 BLOCKED 2026-10-26 on `tests/environment`'s cost, re-posed as §9 item 648 (attempt 2 of 4). **Step 1b green 2026-10-26 (item 648; a two-pole Cole–Cole medium, α = 0.30 / 0.35, on `TH-6`'s plane wave; registered record):** (Ev) evaluator vs `cmath` ≤ 1.7e-16; α / β within 0.042 % / 0.095 % at 127.74 MHz (24³, rate 0.9994) and 0.033 % / 0.052 % at 63.87 MHz; Debye-reduced and wrong-f controls miss α by 20.0 % / 30.1 % (`20261026T094151Z_MAT-3-step1b-W1b.log`, 60 s; module bytes pinned by `20261026T094144Z_MAT-3-step1b-W0b.log`; `tests/environment` `20261026T094304Z_MAT-3-step1b-Wenv.log`, 103 s). Consumed by `Study` (`WF-15` step 8, §9 item 653, 2026-10-27). *History: `docs/planning/chunks/MAT-3.md`.* | 🟡 
