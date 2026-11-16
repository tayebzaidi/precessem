# OPS-77 — how close a 1e-12 cross-solve gate sits to the known MUMPS drift

*History file for the §7 `OPS-77` row (`OPS-46`). Opened by the 2026-11-02 03:00 daily review, lane 1, from `OPS-74`'s disclosure and known-issues' `OPS-43` observation.*

**The question.** Known-issues' 2026-09-10 observation records that solves are not bit-reproducible, with a drift of up to ≈ 4.5e-11 on a lossy fixture at `-n 2` (`MAT-1` probe, its second-fixture row), and states that any solve-output record tighter than ≈ 5e-11 is exposed. Since then `IDENTITY_BAND` = 1e-12 has been used across two separate solves: the `WF-15` / `EX-167` (S) gate, a `Study` against the direct build. `OPS-74` read that (S) at 1.295e-15 and then 6.598e-14 with no code change, ≈ 51 ×, and `EX-172`'s 10 MHz `C_mode` moved rel ≈ 1.4e-13 between windows. **Those are the same phenomenon the observation names, so they need no new diagnosis.** What nobody has measured is the distribution on this gate: whether it is 1e-14 noise that clears the band forever, or a tail that will fail one window in N for no reason.

**Done-when.** A table of every recorded (S) reading on the fixture, plus two new identical windows, with the largest read against `IDENTITY_BAND`; the known-issues observation gains the row; the weekly is told which branch below fired.

### §9 item 755 — 2026-11-02 (lane 1, known-issues; `OPS-74`'s disclosure)

**Why.** Red gates that are noise cost slots and teach people to distrust reds. The rule stands regardless: **a band is never loosened to swallow this.** The measurement only tells the weekly whether the gate's design needs its attention, as `OPS-43` (d)'s did.

**Anchor.** `IDENTITY_BAND` 1e-12, imported, on (S) in each window, asserted by the example's own gate. New: `scripts/testing/read_s_identity_record.py`, stdlib only, reading each listed log through the `.gz`-safe route and printing log, value, multiple of the band.

**Records to read, no solve.** Every `EX-167*` log, the `EX-167-close` log and every `WF-15-s10*`, `WF-15-s11*` log carrying an (S) line. Grep first and list what you found. Each reading is one sample of the same comparison; say which share a process.

**New windows.** `mri:37` twice, identical, in the rule's form:
`scripts/testing/run_and_log.sh OPS-77-mri37-r1 "./scripts/run_examples.sh -e mri:37 -n 2 -t 300"`, and `-r2`. Each ≈ 203 s.

**Negative control.** In `mri:38`'s closing log (`20261102T060640Z_EX-170.log`) (S) reads exactly 0.000e+00 between the sweep and point solves at both frequencies: read and cite it. That is the floor when the comparison shares a route, so a non-zero read in `mri:37` is the separate build's and not the reader's.

**Tier, ranks, wall clock.** Heavy, `-n 2`, two windows ≈ 203 s each; the reader is seconds. **Class:** ops / docs. **Cost:** ≈ 5.5 measured / 20 predicted min; expect ≈ 10 from the two windows. **Depends on:** nothing. Pre-check `pgrep -c python3` reads 0.

**Traps already paid for.** **Never `gunzip` a record to get green;** read through `tests/log_records.read_log_text` or its stdlib twin. `run_examples.sh` needs no `-T` (`OPS-74`). Do not edit `37_implant_rod_through_a_study.py` or its guide: a window here closes nothing about `EX-167`. Count verdicts per rank by token. `FEM_EM_SETUP_FIGURES` stays unset in both windows, since the earlier pair differed in it.

**Scope.** One gate on one fixture. **Not claimed:** a cause (the observation's stays undiagnosed), any other gate's exposure, a new band, a re-run of any closure.

**Negative result.** Pre-registered branches, the review's and then the weekly's to read:
- **All reads ≤ 1e-13 (≥ 10 × clear):** the gate is safe in practice; the row says so and `OPS-77` ✅.
- **Any read in (1e-13, 1e-12]:** ✅ as a record, and the row names the tail; the weekly decides whether the gate's design changes.
- **Any read > 1e-12:** the window is red by the imported band: stop, keep the log, BLOCKED; **it is a finding about the gate, and nothing moves the band.**

No `PREDICTED RED:` line.

**Status it can move.** `OPS-77` ✅; known-issues' `OPS-43` observation gains a third-fixture row with a measured distribution.

**Outcome (2026-11-02T12:23Z, 06:00 CST slot, third item): DONE — `OPS-77` ✅ as a record; the second branch fired, any read in (1e-13, 1e-12].** Two identical windows, `FEM_EM_SETUP_FIGURES` unset, no example edited. Both Status 0: `20261102T121354Z_OPS-77-mri37-r1.log` 204 s (`:5586–5587`), `20261102T121721Z_OPS-77-mri37-r2.log` 201 s (`:5586–5587`). Reader `20261102T122110Z_OPS-77-reader.log`, Status 0, 1 s, through the harness in the container.

The table (`20261102T122110Z_OPS-77-reader.log:36–41`). 28 logs matched; six carry the comparison. Each sample is one window in its own process, so no two share a process.

| log | (S) | × `IDENTITY_BAND` | headroom |
|---|---|---|---|
| `20261101T104952Z_WF-15-s10-W1.log:5588` | 9.801e-14 | 0.0980 | 10.2 × |
| `20261101T212457Z_EX-167.log:5571` | 1.295e-15 | 0.0013 | 772.2 × |
| `20261102T013649Z_EX-167-close.log:5568` | 6.598e-14 | 0.0660 | 15.2 × |
| `20261102T121354Z_OPS-77-mri37-r1.log:5568` | 1.357e-15 | 0.0014 | 736.9 × |
| `20261102T121721Z_OPS-77-mri37-r2.log:5568` | **1.147e-13** | 0.1147 | 8.7 × |

**Negative control.** In `20261102T060640Z_EX-170.log` (S) reads 0.000e+00 at 64 MHz (`:5643`) and 128 MHz (`:5646`), the floor when both sides share a route. So the non-zero reads in `mri:37` are the separate build's, not the reader's.

**Reading.** r1 and r2 are the same command on the same tree and read ≈ 85 × apart. Two reads sit near 1.3e-15 and three near 1e-13, with nothing between; five samples do not show whether that is two modes. The largest stands 8.7 × below the band, so no window has been red. Five samples set no failure rate: this says the tail reaches ≈ 0.11 × the band, not how often it reaches 1 ×.

**Beside the table, another fixture:** `WF-16`'s four-port Study against direct, (I) 1.004e-13 (`20261101T105341Z_WF-15-s10-Wreg.log:3716`). It is outside this gate and did not set the branch, which r2 sets alone; it also lands above 1e-13, so the reading is not the rod's alone. Not claimed as a second sample.

**Disclosed.** (1) The brief's globs also caught `WF-15-s10-Wreg`'s `[WF-16]` line, reported apart for the reason above. (2) `WF-15-s10-W1`'s line reads the same comparison, printed by the gate module, not the example. (3) The reader takes the band from the line the gate printed, so it restates nothing. The `1e-13` in its `over_1e-13` count is the brief's branch line, not a band. (4) The reader has no test of its own; its output is the table above. (5) ≈ 85 × is hand arithmetic on the two prints. **No band moved; no closure re-run.**

**Status moved.** `OPS-77` ✅; known-issues' `OPS-43` observation gains the row; §9 item 755 DONE; the weekly decides whether the gate's design changes.

**Correction (2026-11-02 18:00 review).** The negative control above reads 0.000e+00 as "the floor when both sides share a route". **That premise was wrong:** `examples/mri/38_implant_rod_on_a_sweep.py:100`, `:107` are two factorisations, as item 761 found (`20261102T180848Z_OPS-78-census.log:85–86`). So "the non-zero reads in `mri:37` are the separate build's, not the reader's" lost its support. The table and the tail branch do not rest on that line and stand. `OPS-78.md`'s 18:00 ruling carries the rest.
